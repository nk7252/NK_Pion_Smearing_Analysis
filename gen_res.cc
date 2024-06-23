// Headers and Namespaces.
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <chrono>
// #include <numbers> //std::numbers
#include <cmath>
#include <random>
#include <memory>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include <TF1.h>
#include <TMath.h>
#include "TStopwatch.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "Pythia8/Pythia.h" // Include Pythia headers.
using namespace Pythia8;	// Let Pythia8:: be implicit;

// Forward declarations
TF1* ChooseSpectrumFunction(int weightmethod, int PT_Min, int PT_Max, const std::string& particleType);
Pythia8::Vec4 clusterPhoton(Pythia8::Vec4& originalPhoton, int method, double randomE);
Pythia8::Vec4 PositionResSmear(Pythia8::Vec4 photon, double smearingFactorx,double smearingFactory,double smearingFactorz);
bool DeltaRcut(Pythia8::Vec4& Photon1, Pythia8::Vec4& Photon2, float DeltaRcutMax);
bool pTCut(const Pythia8::Vec4& particle, float ptCut);
bool AsymmCutcheck(Pythia8::Vec4& Photon1, Pythia8::Vec4& Photon2, float AsymmCutoff, bool asymcutbool);
void parseArguments(int argc, char* argv[], std::map<std::string, std::string>& params);

int main(int argc, char* argv[]){
    // Define default parameters
    //input params
    std::string particleType = "pion";
    int nParticles = 8 * 1000000;
    int PT_Max = 50;
    float PT_Min = 0;
    double PT_ratio = PT_Min / PT_Max;
    int MassNBins = 50;
    int binres = 1;
    int n_bins = binres * PT_Max;
    int weightMethod = 2;
    bool applyAsymmCut = true;
    float asymmCutValue = 0.6;
    int clusterOverlap = 1;
    float clusterOverlapProb = 0.99;
    float DeltaRcut_MAX = 1.1;
    float pt1cut = 1.5;
    float pt2cut = 1.5;
    float comb_ptcut = 0;
    float ptMaxCut = 50;
    float nclus_ptCut = 0.0;
    //weighting params
    double t = 4.5;
    double w = 0.114;
    double A = 229.6;
    double B = 14.43;
    double n = 8.1028;
    double m_param = 10.654;
    double p0 = 1.466;
    //smearing params
    float smeared_lower_bin_limit = 0.0;
    float smeared_upper_bin_limit = 1.0;
    float smear_factor_a = 0;
    float smear_factor_d = 0.0;
    float posit_smearingFactor = 2.8;
    float smear_factor_basevalue = 0.168;
    float smear_factor_step = 0.001;
    int smear_factor_steps = 1;
    //output params
    bool saveToTree = false;

    // Parse command-line arguments
    std::map<std::string, std::string> params;
    parseArguments(argc, argv, params);

    if (params.find("particleType") != params.end()) particleType = params["particleType"];
    if (params.find("nParticles") != params.end()) nParticles = std::stoi(params["nParticles"]);
    if (params.find("PT_Max") != params.end()) PT_Max = std::stoi(params["PT_Max"]);
    if (params.find("PT_Min") != params.end()) PT_Min = std::stof(params["PT_Min"]);
    if (params.find("weightMethod") != params.end()) weightMethod = std::stoi(params["weightMethod"]);
    if (params.find("applyAsymmCut") != params.end()) applyAsymmCut = std::stoi(params["applyAsymmCut"]);
    if (params.find("asymmCutValue") != params.end()) asymmCutValue = std::stof(params["asymmCutValue"]);
    if (params.find("clusterOverlap") != params.end()) clusterOverlap = std::stoi(params["clusterOverlap"]);
    if (params.find("clusterOverlapProb") != params.end()) clusterOverlapProb = std::stof(params["clusterOverlapProb"]);
    if (params.find("DeltaRcut_MAX") != params.end()) DeltaRcut_MAX = std::stof(params["DeltaRcut_MAX"]);
    if (params.find("pt1cut") != params.end()) pt1cut = std::stof(params["pt1cut"]);
    if (params.find("pt2cut") != params.end()) pt2cut = std::stof(params["pt2cut"]);
    if (params.find("comb_ptcut") != params.end()) comb_ptcut = std::stof(params["comb_ptcut"]);
    if (params.find("ptMaxCut") != params.end()) ptMaxCut = std::stof(params["ptMaxCut"]);
    if (params.find("nclus_ptCut") != params.end()) nclus_ptCut = std::stof(params["nclus_ptCut"]);
    if (params.find("saveToTree") != params.end()) saveToTree = std::stoi(params["saveToTree"]);
    if (params.find("smear_factor_basevalue") != params.end()) smear_factor_basevalue = std::stof(params["smear_factor_basevalue"]);
    if (params.find("smear_factor_step") != params.end()) smear_factor_step = std::stof(params["smear_factor_step"]);
    if (params.find("smear_factor_steps") != params.end()) smear_factor_steps = std::stoi(params["smear_factor_steps"]);

    TStopwatch timer;
    timer.Start();

    std::map<double, std::vector<double>> mass_pt_map;
    TF1* myFunc = ChooseSpectrumFunction(weightMethod, PT_Min, PT_Max, particleType);

    std::vector<std::string> WeightNames = {"EXP", "POWER", "WSHP", "HAGEDORN"};

    for (int smear_factor_itt = 0; smear_factor_itt < smear_factor_steps; smear_factor_itt++) {
        float smear_factor_c = smear_factor_basevalue + smear_factor_step * smear_factor_itt;
        float smear_factor_b = 0.154;

        TFile* output = new TFile(Form("%sFastMC_%f_sqrte_%f_const.root", particleType.c_str(), smear_factor_b, smear_factor_c), "recreate");

        TTree* tree = nullptr;
        if (saveToTree) {
            tree = new TTree("tree", "tree");
            tree->SetMaxTreeSize(500 * 1024 * 1024);
        }

        TH1* h4 = new TH1F("h4", "PT, unweighted", n_bins, PT_Min, PT_Max);
        TH1* h5 = new TH1F("h5", "Photon Pt, unweighted", n_bins, PT_Min, PT_Max);
        TH1* h6 = new TH1F("h6", "inv mass of gamma pair", MassNBins, 0, 1);
        TH1* h8 = new TH1F("h8", "inv mass of Photon pair, smeared", MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
        TH2F* h9 = new TH2F("h9", "Smeared Pt vs Smeared Inv Mass", n_bins, 0, PT_Max, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
        TH1* h10 = new TH1F("h10", "Smeared PT", n_bins, PT_Min, PT_Max);
        TH1* h16 = new TH1F("h16", "Smeared Photon pT", n_bins, PT_Min, PT_Max);
        TH1* h17 = new TH1F("h17", "Photon pT", n_bins, PT_Min, PT_Max);
        TH1* hInvMass_Cutson = new TH1F("hInvMass_Cutson", "PT,nSmeared+no_weight+cuts+pr", MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);

        std::vector<TH1*> hpionpt(WeightNames.size());
        std::vector<TH1*> h2(WeightNames.size());
        std::vector<TH1*> h3(WeightNames.size());
        std::vector<TH1*> h12(WeightNames.size());
        std::vector<TH1*> h20(WeightNames.size());
        std::vector<TH1*> h21(WeightNames.size());
        std::vector<TH2F*> h18(WeightNames.size());
        std::vector<TH1F*> h18_1d(WeightNames.size());
        std::vector<TH2F*> h27(WeightNames.size());
        std::vector<TH2F*> h28(WeightNames.size());
        std::vector<TH2F*> h29(WeightNames.size());
        std::vector<TH2F*> h28_v2(WeightNames.size());
        std::vector<TH2F*> h29_v2(WeightNames.size());
        std::vector<TH2F*> h30(WeightNames.size());
        std::vector<TH1F*> h30_1d(WeightNames.size());
        std::vector<TH2F*> h31(WeightNames.size());
        std::vector<TH1F*> h31_1d(WeightNames.size());
        std::vector<TH1F*> h31_pionspectrum(WeightNames.size());
        std::vector<TH2F*> h34(WeightNames.size());
        std::vector<TH1F*> h34_1d(WeightNames.size());
        std::vector<TH2F*> h35(WeightNames.size());
        std::vector<TH1F*> h35_1d(WeightNames.size());
        std::vector<TH2F*> h100(WeightNames.size());
        std::vector<TH1F*> h100_1d(WeightNames.size());

        for (int p = 0; p < WeightNames.size(); p++) {
            hpionpt[p] = new TH1D(Form("hpionpt_%i", p), Form("Pt no smear + no weight:%s", WeightNames[p].c_str()), n_bins, PT_Min, PT_Max);
            h2[p] = new TH1D(Form("h2_%i", p), Form("Photon Pt:%s", WeightNames[p].c_str()), n_bins, PT_Min, PT_Max);
            h3[p] = new TH1D(Form("h3_%i", p), Form("PT, weighted:%s", WeightNames[p].c_str()), n_bins, PT_Min, PT_Max);
            h12[p] = new TH1F(Form("h12_%i", p), Form("Smeared PT, weighted:%s", WeightNames[p].c_str()), n_bins, PT_Min, PT_Max);
            h20[p] = new TH1F(Form("h20_%i", p), Form("Smeared Photon pT, weighted:%s", WeightNames[p].c_str()), n_bins, PT_Min, PT_Max);
            h21[p] = new TH1F(Form("h21_%i", p), Form("Photon pT, weighted:%s", WeightNames[p].c_str()), n_bins, PT_Min, PT_Max);
            h18[p] = new TH2F(Form("h18_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
            h18_1d[p] = new TH1F(Form("h18_1d_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted:%s", WeightNames[p].c_str()), MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
            h27[p] = new TH2F(Form("h27_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. cluster:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
            h28[p] = new TH2F(Form("h28_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. cluster and Asymm Cut:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
            h29[p] = new TH2F(Form("h29_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. Asymm Cut:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
            h28_v2[p] = new TH2F(Form("h28_v2_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. Cluster+Asymm+pos res:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
            h29_v2[p] = new TH2F(Form("h29_v2_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. Asymm+Pos res:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
            h30[p] = new TH2F(Form("h30_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. Blair's cuts, no pos.res no occupancy:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, 0, smeared_upper_bin_limit);
            h30_1d[p] = new TH1F(Form("h30_1d_%i", p), Form("Smeared Inv Mass, weighted. Blair's cuts, no pos.res no occupancy:%s", WeightNames[p].c_str()), MassNBins, 0, smeared_upper_bin_limit);
            h31[p] = new TH2F(Form("h31_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. Blair's cuts+pos res:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, 0, smeared_upper_bin_limit);
            h31_1d[p] = new TH1F(Form("h31_1d_%i", p), Form("Smeared Inv Mass, weighted. Blair's cuts+pos res:%s", WeightNames[p].c_str()), MassNBins, 0, smeared_upper_bin_limit);
            h31_pionspectrum[p] = new TH1F(Form("h31_ps_%i", p), Form("Cuts+pos res, Smeared PT, weighted:%s", WeightNames[p].c_str()), MassNBins, PT_Min, PT_Max);
            h34[p] = new TH2F(Form("h34_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. Position Smearing:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
            h34_1d[p] = new TH1F(Form("h34_1d_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. Position Smearing:%s", WeightNames[p].c_str()), MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
            h35[p] = new TH2F(Form("h35_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. Blair's cuts+cluster:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, 0, smeared_upper_bin_limit);
            h35_1d[p] = new TH1F(Form("h35_1d_%i", p), Form("Smeared Inv Mass, weighted. Blair's cuts+cluster:%s", WeightNames[p].c_str()), MassNBins, 0, smeared_upper_bin_limit);
            h100[p] = new TH2F(Form("h100_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. All Cuts+effects:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
            h100_1d[p] = new TH1F(Form("h100_1d_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. All Cuts+effects:%s", WeightNames[p].c_str()), MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
        }

        std::random_device rd;
        std::random_device rdgamma;
        std::random_device rdgammacluster;
        std::random_device rdgammapositsmr;
        std::mt19937_64 gen(rd());
        std::mt19937_64 gen_gamma(rdgamma());
        std::mt19937_64 gen_gammacluster(rdgammacluster());
        std::mt19937_64 gen_gammapositsmear(rdgammapositsmr());
        std::normal_distribution<double> gammadis(0.0, 1.0);
        std::uniform_real_distribution<> gammacluster(0, 1.0);
        std::normal_distribution<double> gamma_positsmear(0.0, 1.0);
        std::uniform_real_distribution<> pdis(PT_ratio, 1.0);
        std::uniform_real_distribution<> adis(0.0, 2 * M_PI);

        Pythia pythia;
        pythia.readString("PromptPhoton:all = on");
        if (particleType == "pion") {
            pythia.readString("111:oneChannel = 1 1.0 0 22 22");
        } else if (particleType == "eta") {
            pythia.readString("221:oneChannel = 1 1.0 0 22 22");
        }
        pythia.init();
        pythia.event.clear();

        int particleCount = 0;
        double particleMass = (particleType == "pion") ? 0.1349768 : 0.54786;
        double mT_Scaling = 1.0;

        for (int i = 0; i < nParticles; i++) {
            double azimuthal_ang = adis(gen);
            double Pt = PT_Max * pdis(gen);

            double weight_function[4];
            double inv_yield[4];
            int WeightScale[4] = {1e+14, 1e+5, 1e+14, 1e+14};

            if (particleType == "eta") {
                mT_Scaling = 0.5 * pow((1.2 + sqrt(pow(0.54786, 2) + pow(Pt, 2))) / (1.2 + sqrt(pow(0.1349768, 2) + pow(Pt, 2))), -10);
            }

            weight_function[0] = exp(-Pt / 0.3);
            inv_yield[0] = mT_Scaling * WeightScale[0] * Pt * weight_function[0];

            weight_function[1] = pow(Pt, -8.14);
            inv_yield[1] = mT_Scaling * WeightScale[1] * Pt * weight_function[1];

            weight_function[2] = ((1 / (1 + exp((Pt - t) / w))) * A / pow(1 + Pt / p0, m_param) + (1 - (1 / (1 + exp((Pt - t) / w)))) * B / (pow(Pt, n)));
            inv_yield[2] = mT_Scaling * WeightScale[2] * Pt * weight_function[2];

            weight_function[3] = A / pow(1 + Pt / p0, m_param);
            inv_yield[3] = mT_Scaling * WeightScale[3] * Pt * weight_function[3];

            double particle_px = Pt * cos(azimuthal_ang);
            double particle_py = Pt * sin(azimuthal_ang);
            double particle_E = sqrt(particleMass * particleMass + Pt * Pt + 0.0);

            if (pythia.event.append((particleType == "pion") ? 111 : 221, 23, 0, 0, particle_px, particle_py, 0.0, particle_E, particleMass)) {
                particleCount++;
            } else {
                std::cout << "Failed to append particle at iteration " << i << std::endl;
            }
        }

        pythia.moreDecays();

        for (int i = 0; i < pythia.event.size(); i++) {
            if (pythia.event[i].id() == ((particleType == "pion") ? 111 : 221)) {
                int Gamma_daughters[2] = {pythia.event[i].daughter1(), pythia.event[i].daughter2()};
                double Pt = pythia.event[i].pT();

                if (particleType == "eta") {
                    mT_Scaling = 0.5 * pow((1.2 + sqrt(pow(0.54786, 2) + pow(Pt, 2))) / (1.2 + sqrt(pow(0.1349768, 2) + pow(Pt, 2))), -10);
                }

                double weight_function[4];
                double inv_yield[4];
                int WeightScale[4] = {1e+14, 1e+5, 1e+14, 1e+14};

                weight_function[0] = exp(-Pt / 0.3);
                inv_yield[0] = mT_Scaling * WeightScale[0] * Pt * weight_function[0];

                weight_function[1] = pow(Pt, -8.14);
                inv_yield[1] = mT_Scaling * WeightScale[1] * Pt * weight_function[1];

                weight_function[2] = ((1 / (1 + exp((Pt - t) / w))) * A / pow(1 + Pt / p0, m_param) + (1 - (1 / (1 + exp((Pt - t) / w)))) * B / (pow(Pt, n)));
                inv_yield[2] = mT_Scaling * WeightScale[2] * Pt * weight_function[2];

                weight_function[3] = A / pow(1 + Pt / p0, m_param);
                inv_yield[3] = mT_Scaling * WeightScale[3] * Pt * weight_function[3];

                if (pythia.event[Gamma_daughters[0]].id() == 22 && pythia.event[Gamma_daughters[1]].id() == 22) {
                    Pythia8::Vec4 gamma_lorentz[3];
                    Pythia8::Vec4 gamma_smeared[3];
                    Pythia8::Vec4 gamma_cluster[3];
                    Pythia8::Vec4 gamma_cluster_asymm[3];
                    Pythia8::Vec4 gamma_Blair_Cuts[3];
                    Pythia8::Vec4 gamma_All_Cuts[3];
                    Pythia8::Vec4 gamma_position_smear[3];
                    Pythia8::Vec4 gamma_Blair_position[3];

                    gamma_lorentz[0] = pythia.event[Gamma_daughters[0]].p();
                    gamma_lorentz[1] = pythia.event[Gamma_daughters[1]].p();
                    gamma_lorentz[2] = gamma_lorentz[0] + gamma_lorentz[1];
                    double inv_mass = gamma_lorentz[2].mCalc();

                    double scale_factor1 = sqrt(pow(smear_factor_b, 2) / gamma_lorentz[0].e() + pow(smear_factor_c, 2) + pow(smear_factor_d, 2));
                    double scale_factor2 = sqrt(pow(smear_factor_b, 2) / gamma_lorentz[1].e() + pow(smear_factor_c, 2) + pow(smear_factor_d, 2));

                    double smear_factor1 = scale_factor1 * gammadis(gen_gamma) + 1;
                    double smear_factor2 = scale_factor2 * gammadis(gen_gamma) + 1;

                    gamma_smeared[0] = smear_factor1 * pythia.event[Gamma_daughters[0]].p();
                    gamma_smeared[1] = smear_factor2 * pythia.event[Gamma_daughters[1]].p();

                    for (int photclust = 0; photclust < 2; photclust++) {
                        if (gammacluster(gen_gammacluster) > clusterOverlapProb && clusterOverlap == 1) {
                            gamma_cluster[photclust] = clusterPhoton(gamma_smeared[photclust], 2, myFunc->GetRandom());
                            gamma_cluster_asymm[photclust] = gamma_cluster[photclust];
                        } else {
                            gamma_cluster[photclust] = gamma_smeared[photclust];
                            gamma_cluster_asymm[photclust] = gamma_cluster[photclust];
                        }
                    }

                    gamma_All_Cuts[0] = PositionResSmear(gamma_cluster_asymm[0], posit_smearingFactor * gamma_positsmear(gen_gammapositsmear), posit_smearingFactor * gamma_positsmear(gen_gammapositsmear), posit_smearingFactor * gamma_positsmear(gen_gammapositsmear));
                    gamma_All_Cuts[1] = PositionResSmear(gamma_cluster_asymm[1], posit_smearingFactor * gamma_positsmear(gen_gammapositsmear), posit_smearingFactor * gamma_positsmear(gen_gammapositsmear), posit_smearingFactor * gamma_positsmear(gen_gammapositsmear));
                    gamma_All_Cuts[2] = gamma_All_Cuts[0] + gamma_All_Cuts[1];

                    gamma_Blair_position[0] = PositionResSmear(gamma_smeared[0], posit_smearingFactor * gamma_positsmear(gen_gammapositsmear), posit_smearingFactor * gamma_positsmear(gen_gammapositsmear), posit_smearingFactor * gamma_positsmear(gen_gammapositsmear));
                    gamma_Blair_position[1] = PositionResSmear(gamma_smeared[1], posit_smearingFactor * gamma_positsmear(gen_gammapositsmear), posit_smearingFactor * gamma_positsmear(gen_gammapositsmear), posit_smearingFactor * gamma_positsmear(gen_gammapositsmear));
                    gamma_Blair_position[2] = gamma_Blair_position[0] + gamma_Blair_position[1];

                    gamma_smeared[2] = gamma_smeared[0] + gamma_smeared[1];
                    gamma_cluster[2] = gamma_cluster[0] + gamma_cluster[1];
                    gamma_cluster_asymm[2] = gamma_cluster_asymm[0] + gamma_cluster_asymm[1];

                    gamma_Blair_Cuts[0] = gamma_cluster_asymm[0];
                    gamma_Blair_Cuts[1] = gamma_cluster_asymm[1];
                    gamma_Blair_Cuts[2] = gamma_cluster_asymm[2];

                    gamma_position_smear[0] = PositionResSmear(gamma_smeared[0], posit_smearingFactor * gamma_positsmear(gen_gammapositsmear), posit_smearingFactor * gamma_positsmear(gen_gammapositsmear), posit_smearingFactor * gamma_positsmear(gen_gammapositsmear));
                    gamma_position_smear[1] = PositionResSmear(gamma_smeared[1], posit_smearingFactor * gamma_positsmear(gen_gammapositsmear), posit_smearingFactor * gamma_positsmear(gen_gammapositsmear), posit_smearingFactor * gamma_positsmear(gen_gammapositsmear));

                    gamma_position_smear[2] = gamma_position_smear[0] + gamma_position_smear[1];

                    double inv_mass_smeared = gamma_smeared[2].mCalc();
                    h4->Fill(Pt);
                    h6->Fill(inv_mass);
                    h8->Fill(inv_mass_smeared);
                    h9->Fill(gamma_smeared[2].pT(), inv_mass_smeared);
                    h10->Fill(gamma_smeared[2].pT());
                    h17->Fill(gamma_lorentz[0].pT());
                    h17->Fill(gamma_lorentz[1].pT());
                    h16->Fill(gamma_smeared[0].pT());
                    h16->Fill(gamma_smeared[1].pT());

                    for (int p = 0; p < WeightNames.size(); p++) {
                        h20[p]->Fill(gamma_smeared[0].pT(), inv_yield[p]);
                        h20[p]->Fill(gamma_smeared[1].pT(), inv_yield[p]);
                        h21[p]->Fill(gamma_lorentz[0].pT(), inv_yield[p]);
                        h21[p]->Fill(gamma_lorentz[1].pT(), inv_yield[p]);
                        h3[p]->Fill(Pt, inv_yield[p]);
                        h12[p]->Fill(gamma_smeared[2].pT(), inv_yield[p]);
                        h18[p]->Fill(gamma_smeared[2].pT(), inv_mass_smeared, inv_yield[p]);
                        h18_1d[p]->Fill(inv_mass_smeared, inv_yield[p]);
                        h27[p]->Fill(gamma_cluster[2].pT(), gamma_cluster[2].mCalc(), inv_yield[p]);
                        h34[p]->Fill(gamma_position_smear[2].pT(), gamma_position_smear[2].mCalc(), inv_yield[p]);
                        h34_1d[p]->Fill(gamma_position_smear[2].mCalc(), inv_yield[p]);

                        if (AsymmCutcheck(gamma_smeared[0], gamma_smeared[1], asymmCutValue, applyAsymmCut) == true) {
                            h29[p]->Fill(gamma_smeared[2].pT(), gamma_smeared[2].mCalc(), inv_yield[p]);
                        }
                        if (AsymmCutcheck(gamma_cluster_asymm[0], gamma_cluster_asymm[1], asymmCutValue, applyAsymmCut) == true) {
                            h28[p]->Fill(gamma_cluster_asymm[2].pT(), gamma_cluster_asymm[2].mCalc(), inv_yield[p]);
                        }
                        if (AsymmCutcheck(gamma_position_smear[0], gamma_position_smear[1], asymmCutValue, applyAsymmCut) == true) {
                            h29_v2[p]->Fill(gamma_position_smear[2].pT(), gamma_position_smear[2].mCalc(), inv_yield[p]);
                        }
                        if (AsymmCutcheck(gamma_All_Cuts[0], gamma_All_Cuts[1], asymmCutValue, applyAsymmCut) == true) {
                            h28_v2[p]->Fill(gamma_All_Cuts[2].pT(), gamma_All_Cuts[2].mCalc(), inv_yield[p]);
                        }

                        if (DeltaRcut(gamma_smeared[0], gamma_smeared[1], DeltaRcut_MAX) == false && AsymmCutcheck(gamma_smeared[0], gamma_smeared[1], asymmCutValue, applyAsymmCut) == true && pTCut(gamma_smeared[0], pt1cut) == true && pTCut(gamma_smeared[1], pt2cut) == true && nclus_ptCut < gamma_smeared[0].pT() && gamma_smeared[0].pT() < ptMaxCut && nclus_ptCut < gamma_smeared[1].pT() && gamma_smeared[1].pT() < ptMaxCut && gamma_smeared[2].pT() > comb_ptcut * (pt1cut + pt2cut)) {
                            h30[p]->Fill(gamma_smeared[2].pT(), gamma_smeared[2].mCalc(), inv_yield[p]);
                            h30_1d[p]->Fill(gamma_smeared[2].mCalc(), inv_yield[p]);
                        }

                        if (DeltaRcut(gamma_Blair_position[0], gamma_Blair_position[1], DeltaRcut_MAX) == false && AsymmCutcheck(gamma_Blair_position[0], gamma_Blair_position[1], asymmCutValue, applyAsymmCut) == true && pTCut(gamma_Blair_position[0], pt1cut) == true && pTCut(gamma_Blair_position[1], pt2cut) == true && nclus_ptCut < gamma_Blair_position[0].pT() && gamma_Blair_position[0].pT() < ptMaxCut && nclus_ptCut < gamma_Blair_position[1].pT() && gamma_Blair_position[1].pT() < ptMaxCut && gamma_Blair_position[2].pT() > comb_ptcut * (pt1cut + pt2cut)) {
                            h31[p]->Fill(gamma_Blair_position[2].pT(), gamma_Blair_position[2].mCalc(), inv_yield[p]);
                            h31_1d[p]->Fill(gamma_Blair_position[2].mCalc(), inv_yield[p]);
                            h31_pionspectrum[p]->Fill(gamma_Blair_position[2].pT(), inv_yield[p]);
                            hpionpt[p]->Fill(gamma_lorentz[2].pT());
                            hInvMass_Cutson->Fill(gamma_Blair_position[2].mCalc());
                        }

                        if (DeltaRcut(gamma_Blair_Cuts[0], gamma_Blair_Cuts[1], DeltaRcut_MAX) == false && AsymmCutcheck(gamma_Blair_Cuts[0], gamma_Blair_Cuts[1], asymmCutValue, applyAsymmCut) == true && pTCut(gamma_Blair_Cuts[0], pt1cut) == true && pTCut(gamma_Blair_Cuts[1], pt2cut) == true && nclus_ptCut < gamma_Blair_Cuts[0].pT() && gamma_Blair_Cuts[0].pT() < ptMaxCut && nclus_ptCut < gamma_Blair_Cuts[1].pT() && gamma_Blair_Cuts[1].pT() < ptMaxCut && gamma_Blair_Cuts[2].pT() > comb_ptcut * (pt1cut + pt2cut)) {
                            h35[p]->Fill(gamma_Blair_Cuts[2].pT(), gamma_Blair_Cuts[2].mCalc(), inv_yield[p]);
                            h35_1d[p]->Fill(gamma_Blair_Cuts[2].mCalc(), inv_yield[p]);
                        }

                        if (DeltaRcut(gamma_All_Cuts[0], gamma_All_Cuts[1], DeltaRcut_MAX) == false && AsymmCutcheck(gamma_All_Cuts[0], gamma_All_Cuts[1], asymmCutValue, applyAsymmCut) == true && pTCut(gamma_All_Cuts[0], pt1cut) == true && pTCut(gamma_All_Cuts[1], pt2cut) == true && nclus_ptCut < gamma_All_Cuts[0].pT() && gamma_All_Cuts[0].pT() < ptMaxCut && nclus_ptCut < gamma_All_Cuts[1].pT() && gamma_All_Cuts[1].pT() < ptMaxCut && gamma_All_Cuts[2].pT() > comb_ptcut * (pt1cut + pt2cut)) {
                            h100[p]->Fill(gamma_All_Cuts[2].pT(), gamma_All_Cuts[2].mCalc(), inv_yield[p]);
                            h100_1d[p]->Fill(gamma_All_Cuts[2].mCalc(), inv_yield[p]);
                        }
                    }

                    try {
                        mass_pt_map[floor(Pt) + 1].push_back(inv_mass_smeared);
                    } catch (...) {
                        mass_pt_map.insert({floor(Pt) + 1, std::vector<double>()});
                    }

                    for (int j : Gamma_daughters) {
                        h5->Fill(pythia.event[j].pT());
                    }
                }
            }
        }

        double_t realtime = timer.RealTime();
        double_t cputime = timer.CpuTime();

        printf("real time =%f\n", realtime);
        printf("CPU time=%f\n", cputime);

        output->Write();
        output->Close();

        if (saveToTree) {
            delete tree;
        }
        delete output;
    }

    return 0;
}

TF1* ChooseSpectrumFunction(int weightmethod, int PT_Min, int PT_Max, const std::string& particleType) {
    TF1 *myFunc = nullptr;
    if (particleType == "pion") {
        if (weightmethod == 0) {
            myFunc = new TF1("myFunc", [](double *x, double *par) {
                return x[0] * par[0] * TMath::Exp(-x[0] / 0.3);
            }, PT_Min, PT_Max, 1);
            Double_t initialParameters[1] = {1.0};
            myFunc->SetParameters(initialParameters);
        } else if (weightmethod == 1) {
            myFunc = new TF1("myFunc", [](double *x, double *par) {
                return x[0] * par[0] * pow(x[0], -8.14);
            }, PT_Min, PT_Max, 1);
            Double_t initialParameters[1] = {1.0};
            myFunc->SetParameters(initialParameters);
        } else if (weightmethod == 2) {
            myFunc = new TF1("myFunc", [](double *x, double *par) {
                return x[0] * ((1 / (1 + exp((x[0] - par[0]) / par[1]))) * par[2] / pow(1 + x[0] / par[3], par[4]) + (1 - (1 / (1 + exp((x[0] - par[0]) / par[1])))) * par[5] / (pow(x[0], par[6])));
            }, PT_Min, PT_Max, 7);
            myFunc->SetParameters(4.5, 0.114, 229.6, 1.466, 10.654, 14.43, 8.1028);
        } else if (weightmethod == 3) {
            myFunc = new TF1("myFunc", [](double *x, double *par) {
                return x[0] * (par[0] / pow(1 + x[0] / par[1], par[2]));
            }, PT_Min, PT_Max, 3);
            myFunc->SetParameters(229.6, 1.466, 10.654);
        } else {
            std::cout << "Error: No Weight function found" << std::endl;
        }
    } else if (particleType == "eta") {
        if (weightmethod == 0) {
            myFunc = new TF1("myFunc", [](double *x, double *par) {
                return x[0] * par[1] * pow((par[2] + sqrt(pow(par[3], 2) + pow(x[0], 2))) / (par[2] + sqrt(pow(par[4], 2) + pow(x[0], 2))), -par[5]) * par[0] * TMath::Exp(-x[0] / 0.3);
            }, PT_Min, PT_Max, 6);
            myFunc->SetParameters(1.0, 0.5, 1.2, 0.54786, 0.1349768, 10);
        } else if (weightmethod == 1) {
            myFunc = new TF1("myFunc", [](double *x, double *par) {
                return x[0] * par[1] * pow((par[2] + sqrt(pow(par[3], 2) + pow(x[0], 2))) / (par[2] + sqrt(pow(par[4], 2) + pow(x[0], 2))), -par[5]) * par[0] * pow(x[0], -8.14);
            }, PT_Min, PT_Max, 6);
            myFunc->SetParameters(1.0, 0.5, 1.2, 0.54786, 0.1349768, 10);
        } else if (weightmethod == 2) {
            myFunc = new TF1("myFunc", [](double *x, double *par) {
                return x[0] * par[7] * pow((par[8] + sqrt(pow(par[9], 2) + pow(x[0], 2))) / (par[8] + sqrt(pow(par[10], 2) + pow(x[0], 2))), -par[11]) * ((1 / (1 + exp((x[0] - par[0]) / par[1]))) * par[2] / pow(1 + x[0] / par[3], par[4]) + (1 - (1 / (1 + exp((x[0] - par[0]) / par[1])))) * par[5] / (pow(x[0], par[6])));
            }, PT_Min, PT_Max, 12);
            Double_t params[] = {4.5, 0.114, 229.6, 1.466, 10.654, 14.43, 8.1028, 0.5, 1.2, 0.54786, 0.1349768, 10.0};
            myFunc->SetParameters(params);
        } else if (weightmethod == 3) {
            myFunc = new TF1("myFunc", [](double *x, double *par) {
                return x[0] * par[3] * pow((par[4] + sqrt(pow(par[5], 2) + pow(x[0], 2))) / (par[4] + sqrt(pow(par[6], 2) + pow(x[0], 2))), -par[7]) * (par[0] / pow(1 + x[0] / par[1], par[2]));
            }, PT_Min, PT_Max, 8);
            myFunc->SetParameters(229.6, 1.466, 10.654, 0.5, 1.2, 0.54786, 0.1349768, 10);
        } else {
            std::cout << "Error: No Weight function found" << std::endl;
        }
    }
    return myFunc;
}

Pythia8::Vec4 clusterPhoton(Pythia8::Vec4& originalPhoton, int method, double randomE) {
    Pythia8::Vec4 newPhoton;
    Pythia8::Rndm rndm;
    if (method == 2) {
        newPhoton.e(randomE);
        newPhoton.px(originalPhoton.px() * (randomE / originalPhoton.e()));
        newPhoton.py(originalPhoton.py() * (randomE / originalPhoton.e()));
        newPhoton.pz(originalPhoton.pz() * (randomE / originalPhoton.e()));
    }
    return newPhoton + originalPhoton;
}

Pythia8::Vec4 PositionResSmear(Pythia8::Vec4 photon, double smearingFactorx, double smearingFactory, double smearingFactorz) {
    double energy = photon.e();
    double posx = 900 * photon.px() / photon.e();
    double posy = 900 * photon.py() / photon.e();
    double posz = 900 * photon.pz() / photon.e();
    posx += smearingFactorx / sqrt(energy);
    posy += smearingFactory / sqrt(energy);
    posz += smearingFactorz / sqrt(energy);
    double px_smear = posx * photon.e() / 900;
    double py_smear = posy * photon.e() / 900;
    double pz_smear = posz * photon.e() / 900;
    double E_New = sqrt(px_smear * px_smear + py_smear * py_smear + pz_smear * pz_smear);
    double energyscale = energy / E_New;
    double xnew = energyscale * px_smear;
    double ynew = energyscale * py_smear;
    double znew = energyscale * pz_smear;
    Pythia8::Vec4 smearedPhoton(xnew, ynew, znew, energy);
    return smearedPhoton;
}

bool DeltaRcut(Pythia8::Vec4& Photon1, Pythia8::Vec4& Photon2, float DeltaRcutMax) {
    double dEta = Photon1.eta() - Photon2.eta();
    double dPhi = acos(cos(Photon1.phi() - Photon2.phi()));
    double deltaR = sqrt(dEta * dEta + dPhi * dPhi);
    return deltaR > DeltaRcutMax;
}

bool pTCut(const Pythia8::Vec4& particle, float ptCut) {
    double pT = particle.pT();
    return pT > ptCut;
}

bool AsymmCutcheck(Pythia8::Vec4& Photon1, Pythia8::Vec4& Photon2, float AsymmCutoff, bool asymcutbool) {
    if (asymcutbool == false) {
        return true;
    }
    return abs(Photon1.e() - Photon2.e()) / (Photon1.e() + Photon2.e()) < AsymmCutoff;
}

void parseArguments(int argc, char* argv[], std::map<std::string, std::string>& params) {
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        size_t equalPos = arg.find('=');
        if (equalPos != std::string::npos) {
            std::string key = arg.substr(0, equalPos);
            std::string value = arg.substr(equalPos + 1);
            params[key] = value;
        }
    }
}
