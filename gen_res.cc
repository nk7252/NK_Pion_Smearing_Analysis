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
using namespace Pythia8;    // Let Pythia8:: be implicit;

// methods for cluster overlap scaling.
enum ScalingMethod
{
    EXPONENTIAL,
    RATIONAL,
    HYPERBOLIC_TANGENT
};
// Forward declarations
TF1 *ChooseSpectrumFunction(int weightmethod, int PT_Min, int PT_Max, const std::string &particleType);
Pythia8::Vec4 clusterPhoton(Pythia8::Vec4 &originalPhoton, int method, double randomE);
Pythia8::Vec4 PositionResSmear(Pythia8::Vec4 photon, double smearingFactorx, double smearingFactory, double smearingFactorz);
bool DeltaRcut(Pythia8::Vec4 &Photon1, Pythia8::Vec4 &Photon2, float DeltaRcutMax);
float DeltaR(Pythia8::Vec4 &Photon1, Pythia8::Vec4 &Photon2);
bool pTCut(const Pythia8::Vec4 &particle, float ptCut);
bool eTCut(const Pythia8::Vec4 &particle, float etCut);
bool EtaCut(const Pythia8::Vec4 &particle, float EtaCutValue, bool ApplyEtaCut, bool debug);
bool AsymmCutcheck(Pythia8::Vec4 &Photon1, Pythia8::Vec4 &Photon2, float AsymmCutoff, bool asymcutbool);
void parseArguments(int argc, char *argv[], std::map<std::string, std::string> &params, bool debug);
double DetectorPhotonDistance(Pythia8::Vec4 &photon1, Pythia8::Vec4 &photon2, bool debug);
std::pair<Pythia8::Vec4, Pythia8::Vec4> adjustPhotonEnergiesSymmetric(Pythia8::Vec4 photon1, Pythia8::Vec4 photon2, int method_int, bool debug);
std::pair<Pythia8::Vec4, Pythia8::Vec4> adjustPhotonEnergiesAsymmetric(Pythia8::Vec4 photon1, Pythia8::Vec4 photon2, int method_int, bool debug);
double EnergySharingScale(double distance, ScalingMethod method);
Pythia8::Vec4 ClusterLossSimple(Pythia8::Vec4 &originalPhoton, float randomRoll, float threshold);

int main(int argc, char *argv[])
{
    // compiled with:
    // g++ gen_res.cc -o gen_res -w -I /home/nik/pythia8307/include -O2 -std=c++11 -pedantic -W -Wall -Wshadow -fPIC -pthread  -L/home/nik/pythia8307/lib/ -Wl,-rpath,/home/nik/pythia8307/lib -lpythia8 -ldl -L /home/nik/root/lib -Wl,-rpath,/home/nik/root/lib -lCore `root-config --cflags --glibs`
    // code is run with:
    // ./gen_res -n 8000000 -p 50.0 -m 0.0 -w WSHP -a True -v 0.6 -c True -o 0.99 -d 1.1 -q 1.3 -r 0.7 -b 0.12 -x 0.0 -y 50.0 -z 0.0 -f 2.8 -t False -s 1 -i 0.001
    // Define default parameters
    // input params
    std::string particleType = "Pion";
    // std::string particleType = "Eta";
    int nParticles = 8 * 1000000;
    int PT_Max = 40;
    float PT_Min = 0;
    double PT_ratio = PT_Min / PT_Max;
    std::string weightMethodStr = "WSHP";
    int weightMethod = 2;
    bool applyAsymmCut = true; // need to remove all do this bool
    float asymmCutValue = 0.6;
    int clusterOverlap = 1; // need to remove all do this bool
    float clusterOverlapProb = 0.99;
    float DeltaRcut_MAX = 1.1;
    float pt1cut = 1.3;
    float pt2cut = 0.7;
    float comb_ptcut = 0;
    float ptMaxCut = 100;
    float nclus_ptCut = 0.0;

    int PT_Max_bin = 20; // normally this, but now we want to match fun4all PT_Max;
    int MassNBins = 600;
    int binres = 1;
    int n_bins = 8 * 10; // binres * PT_Max;
    bool Debug = false;  // default should be false ofc.
    float etCut = 1.0;
    bool Apply_Eta_Cut = false; // need to remove all do this bool
    float eta_cut_val = 0.6;
    // weighting params
    double t = 4.5;
    double w = 0.114;
    double A = 229.6;
    double B = 14.43;
    double n = 8.1028;
    double m_param = 10.654;
    double p0 = 1.466;
    // smearing params
    // add first three to gui
    float smeared_lower_bin_limit = 0.0; // mass limits
    float smeared_upper_bin_limit = 1.2;
    float smear_factor_sqrtE = 0.154; // 0.154;
    float smear_factor_a = 0;
    float smear_factor_d = 0.0;
    float posit_smearingFactor = 2.8;
    float smear_factor_const = 0.12;
    float smear_factor_const_step_size = 0.001;
    int smear_factor_const_num_steps = 1;
    float tower_cluster_threshold = 0.07;
    // output params
    bool saveToTree = false;
    // untracked general parameters
    bool Debug_Hists = false;
    // int scaling_method_int=1;

    // Parse command-line arguments
    std::map<std::string, std::string> params;
    parseArguments(argc, argv, params, Debug);

    // Override default parameters with command-line arguments
    if (params.find("particleType") != params.end())
        particleType = params["particleType"];
    if (params.find("nParticles") != params.end())
        nParticles = std::stoi(params["nParticles"]);
    if (params.find("PT_Max") != params.end())
        PT_Max = std::stoi(params["PT_Max"]);
    if (params.find("PT_Min") != params.end())
        PT_Min = std::stof(params["PT_Min"]);
    if (params.find("weightMethod") != params.end())
        weightMethodStr = params["weightMethod"];
    if (params.find("applyAsymmCut") != params.end())
        applyAsymmCut = (params["applyAsymmCut"] == "true");
    if (params.find("asymmCutValue") != params.end())
        asymmCutValue = std::stof(params["asymmCutValue"]);
    if (params.find("clusterOverlap") != params.end())
        clusterOverlap = (params["clusterOverlap"] == "true");
    if (params.find("clusterOverlapProb") != params.end())
        clusterOverlapProb = std::stof(params["clusterOverlapProb"]);
    if (params.find("DeltaRcut_MAX") != params.end())
        DeltaRcut_MAX = std::stof(params["DeltaRcut_MAX"]);
    if (params.find("pt1cut") != params.end())
        pt1cut = std::stof(params["pt1cut"]);
    if (params.find("pt2cut") != params.end())
        pt2cut = std::stof(params["pt2cut"]);
    if (params.find("comb_ptcut") != params.end())
        comb_ptcut = std::stof(params["comb_ptcut"]);
    if (params.find("ptMaxCut") != params.end())
        ptMaxCut = std::stof(params["ptMaxCut"]);
    if (params.find("nclus_ptCut") != params.end())
        nclus_ptCut = std::stof(params["nclus_ptCut"]);
    if (params.find("posit_smearingFactor") != params.end())
        posit_smearingFactor = std::stof(params["posit_smearingFactor"]);
    if (params.find("saveToTree") != params.end())
        saveToTree = (params["saveToTree"] == "true");
    if (params.find("smear_factor_const") != params.end())
        smear_factor_const = std::stof(params["smear_factor_const"]);
    if (params.find("smear_factor_const_step_size") != params.end())
        smear_factor_const_step_size = std::stof(params["smear_factor_const_step_size"]);
    if (params.find("smear_factor_const_num_steps") != params.end())
        smear_factor_const_num_steps = std::stoi(params["smear_factor_const_num_steps"]);
    if (params.find("PT_Max_bin") != params.end())
        PT_Max_bin = std::stoi(params["PT_Max_bin"]);
    if (params.find("MassNBins") != params.end())
        MassNBins = std::stoi(params["MassNBins"]);
    if (params.find("binres") != params.end())
        binres = std::stoi(params["binres"]);
    if (params.find("n_bins") != params.end())
        n_bins = std::stoi(params["n_bins"]);
    if (params.find("Debug") != params.end())
        Debug = (params["Debug"] == "true");
    if (params.find("etCut") != params.end())
        etCut = std::stof(params["etCut"]);
    if (params.find("Apply_Eta_Cut") != params.end())
        Apply_Eta_Cut = (params["Apply_Eta_Cut"] == "true");
    if (params.find("eta_cut_val") != params.end())
        eta_cut_val = std::stof(params["eta_cut_val"]);
    if (params.find("smeared_lower_bin_limit") != params.end())
        smeared_lower_bin_limit = std::stof(params["smeared_lower_bin_limit"]);
    if (params.find("smeared_upper_bin_limit") != params.end())
        smeared_upper_bin_limit = std::stof(params["smeared_upper_bin_limit"]);
    if (params.find("smear_factor_sqrtE") != params.end())
        smear_factor_sqrtE = std::stof(params["smear_factor_sqrtE"]);

    // Map weightMethodStr to weightMethod integer
    if (weightMethodStr == "EXP")
    {
        weightMethod = 0;
    }
    else if (weightMethodStr == "POWER")
    {
        weightMethod = 1;
    }
    else if (weightMethodStr == "WSHP")
    {
        weightMethod = 2;
    }
    else if (weightMethodStr == "HAGEDORN")
    {
        weightMethod = 3;
    }
    else
    {
        std::cerr << "Unknown weight method: " << weightMethodStr << std::endl;
        return 1;
    }

    // Debugging: Print all parameters after parsing
    std::cout << "\n Parameters after parsing: \n"
              << std::endl;
    std::cout << "particleType: " << particleType << std::endl;
    std::cout << "nParticles: " << nParticles << std::endl;
    std::cout << "PT_Max: " << PT_Max << std::endl;
    std::cout << "PT_Min: " << PT_Min << std::endl;
    std::cout << "weightMethod: " << weightMethod << std::endl;
    std::cout << "applyAsymmCut: " << applyAsymmCut << std::endl;
    std::cout << "asymmCutValue: " << asymmCutValue << std::endl;
    std::cout << "clusterOverlap: " << clusterOverlap << std::endl;
    std::cout << "clusterOverlapProb: " << clusterOverlapProb << std::endl;
    std::cout << "DeltaRcut_MAX: " << DeltaRcut_MAX << std::endl;
    std::cout << "pt1cut: " << pt1cut << std::endl;
    std::cout << "pt2cut: " << pt2cut << std::endl;
    std::cout << "comb_ptcut: " << comb_ptcut << std::endl;
    std::cout << "ptMaxCut: " << ptMaxCut << std::endl;
    std::cout << "nclus_ptCut: " << nclus_ptCut << std::endl;
    std::cout << "posit_smearingFactor: " << posit_smearingFactor << std::endl;
    std::cout << "saveToTree: " << saveToTree << std::endl;
    std::cout << "smear_factor_const: " << smear_factor_const << std::endl;
    std::cout << "smear_factor_const_step_size: " << smear_factor_const_step_size << std::endl;
    std::cout << "smear_factor_const_num_steps: " << smear_factor_const_num_steps << std::endl;
    std::cout << "PT_Max_bin: " << PT_Max_bin << std::endl;
    std::cout << "MassNBins: " << MassNBins << std::endl;
    std::cout << "binres: " << binres << std::endl;
    std::cout << "n_bins: " << n_bins << std::endl;
    std::cout << "Debug: " << std::boolalpha << Debug << std::endl;
    std::cout << "etCut: " << etCut << std::endl;
    std::cout << "Apply_Eta_Cut: " << std::boolalpha << Apply_Eta_Cut << std::endl;
    std::cout << "eta_cut_val: " << eta_cut_val << std::endl;
    std::cout << "smeared_lower_bin_limit: " << smeared_lower_bin_limit << std::endl;
    std::cout << "smeared_upper_bin_limit: " << smeared_upper_bin_limit << std::endl;
    std::cout << "smear_factor_sqrtE: " << smear_factor_sqrtE << std::endl;

    TH2F *h_original = nullptr;
    TH2F *h_ClusterERatio = nullptr;
    /*
    TFile* frw = TFile::Open("/pioncode/rootfiles/", "READ");
    h_original = (TH2F *)frw->Get("h_clus_ERatio_2d");
    std::cout << "Using eta spectrum" << std::endl;
    if (h_original == nullptr)
    {
        std::cerr << "Error: Original histogram not found!" << std::endl;
    }
    else
    {
        h_ClusterERatio = (TH2F *)h_original->Clone("h_ClusterERatio_clone");
        if (h_ClusterERatio == nullptr)
        {
            std::cerr << "Error: Cloning of histogram failed!" << std::endl;
        }
        else
        {
            h_ClusterERatio->SetDirectory(0);
        }
    }
    // Close the file and clean up
    frw->Close();
    delete frw;
    //*/
    TStopwatch timer;
    timer.Start();

    std::map<double, std::vector<double>> mass_pt_map;
    TF1 *myFunc = ChooseSpectrumFunction(weightMethod, PT_Min, PT_Max, particleType);

    std::vector<std::string> WeightNames = {"EXP", "POWER", "WSHP", "HAGEDORN"};  // weightMethodStr
    std::vector<std::string> ClusterScalingNames = {"EXP", "RATIONAL", "H_TANH"}; // different

    for (int smear_factor_itt = 0; smear_factor_itt < smear_factor_const_num_steps; smear_factor_itt++)
    {
        float smear_factor_c = smear_factor_const + smear_factor_const_step_size * smear_factor_itt;
        float smear_factor_b = smear_factor_sqrtE;

        TFile *output = new TFile(Form("pioncode/rootfiles/%sFastMC_%f_sqrte_%f_const.root", particleType.c_str(), smear_factor_b, smear_factor_c), "recreate");

        TTree *tree = nullptr;
        if (saveToTree)
        {
            tree = new TTree("tree", "tree");
            tree->SetMaxTreeSize(500 * 1024 * 1024);
        }

        TH1 *h4 = new TH1F("h4", "PT, unweighted", n_bins, PT_Min, PT_Max_bin);
        TH1 *h5 = new TH1F("h5", "Photon Pt, unweighted", n_bins, PT_Min, PT_Max_bin);
        TH1 *h6 = new TH1F("h6", "inv mass of gamma pair", MassNBins, 0, 1);
        TH1 *h8 = new TH1F("h8", "inv mass of Photon pair, smeared", MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
        TH2F *h9 = new TH2F("h9", "Smeared Pt vs Smeared Inv Mass", n_bins, 0, PT_Max_bin, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
        TH1 *h10 = new TH1F("h10", "Smeared PT", n_bins, PT_Min, PT_Max_bin);
        TH1 *h16 = new TH1F("h16", "Smeared Photon pT", n_bins, PT_Min, PT_Max_bin);
        TH1 *h17 = new TH1F("h17", "Photon pT", n_bins, PT_Min, PT_Max_bin);
        TH1 *hInvMass_Cutson = new TH1F("hInvMass_Cutson", "PT,nSmeared+no_weight+cuts+pr", MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);

        std::vector<TH1 *> hpionpt(WeightNames.size());
        std::vector<TH1 *> htruthphotondistance_1d(WeightNames.size());
        std::vector<TH2F *> htruthphotondistance(WeightNames.size());
        std::vector<TH1 *> h2(WeightNames.size());
        std::vector<TH1 *> h3(WeightNames.size());
        std::vector<TH1 *> h12(WeightNames.size());
        std::vector<TH1 *> h20(WeightNames.size());
        std::vector<TH1 *> h21(WeightNames.size());
        std::vector<TH2F *> h18(WeightNames.size());
        std::vector<TH1F *> h18_1d(WeightNames.size());
        std::vector<TH2F *> h27(WeightNames.size());
        std::vector<TH2F *> h28(WeightNames.size());
        std::vector<TH2F *> h29(WeightNames.size());
        std::vector<TH2F *> h28_v2(WeightNames.size());
        std::vector<TH2F *> h29_v2(WeightNames.size());
        std::vector<TH2F *> h30(WeightNames.size());
        std::vector<TH1F *> h30_1d(WeightNames.size());
        std::vector<TH2F *> h31(WeightNames.size());
        std::vector<TH1F *> h31_1d(WeightNames.size());
        std::vector<TH1F *> h31_pionspectrum(WeightNames.size());
        std::vector<TH2F *> h34(WeightNames.size());
        std::vector<TH1F *> h34_1d(WeightNames.size());
        std::vector<TH2F *> h35(WeightNames.size());
        std::vector<TH1F *> h35_1d(WeightNames.size());
        std::vector<TH2F *> h100(WeightNames.size());
        std::vector<TH1F *> h100_1d(WeightNames.size());
        std::vector<TH2F *> h101(WeightNames.size());
        std::vector<TH2F *> h102(WeightNames.size());
        std::vector<TH1F *> h101_1d(WeightNames.size());
        std::vector<TH1F *> h101_dr(WeightNames.size());
        std::vector<TH1F *> h101_photon_dist_1d(WeightNames.size());
        std::vector<TH2F *> h101_photon_dist(WeightNames.size());
        std::vector<std::vector<TH2F *>> h101_asymm(WeightNames.size(), std::vector<TH2F *>(ClusterScalingNames.size()));
        std::vector<std::vector<TH2F *>> h101_symm(WeightNames.size(), std::vector<TH2F *>(ClusterScalingNames.size()));
        // std::vector<TH2F *> h101_asymm(WeightNames.size());
        // std::vector<TH2F *> h101_symm(WeightNames.size());

        for (int p = 0; p < WeightNames.size(); p++)
        {
            htruthphotondistance_1d[p] = new TH1F(Form("htruthphotondistance_1d_%i", p), Form("Photon distance distribution, weighted:%s", WeightNames[p].c_str()), 10000, 0, 2000);
            htruthphotondistance[p] = new TH2F(Form("htruthphotondistance_%i", p), Form("pT vs Photon distance distribution, weighted:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max_bin, 10000, 0, 2000);
            if (Debug_Hists)
            {
                hpionpt[p] = new TH1D(Form("hpionpt_%i", p), Form("Pt no smear + no weight:%s", WeightNames[p].c_str()), n_bins, PT_Min, PT_Max_bin);
                h2[p] = new TH1D(Form("h2_%i", p), Form("Photon Pt:%s", WeightNames[p].c_str()), n_bins, PT_Min, PT_Max_bin);
                h3[p] = new TH1D(Form("h3_%i", p), Form("PT, weighted:%s", WeightNames[p].c_str()), n_bins, PT_Min, PT_Max_bin);
                h12[p] = new TH1F(Form("h12_%i", p), Form("Smeared PT, weighted:%s", WeightNames[p].c_str()), n_bins, PT_Min, PT_Max_bin);
                h20[p] = new TH1F(Form("h20_%i", p), Form("Smeared Photon pT, weighted:%s", WeightNames[p].c_str()), n_bins, PT_Min, PT_Max_bin);
                h21[p] = new TH1F(Form("h21_%i", p), Form("Photon pT, weighted:%s", WeightNames[p].c_str()), n_bins, PT_Min, PT_Max_bin);
                h18[p] = new TH2F(Form("h18_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max_bin, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
                h18_1d[p] = new TH1F(Form("h18_1d_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted:%s", WeightNames[p].c_str()), MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
                h27[p] = new TH2F(Form("h27_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. cluster:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max_bin, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
                h28[p] = new TH2F(Form("h28_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. cluster and Asymm Cut:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max_bin, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
                h29[p] = new TH2F(Form("h29_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. Asymm Cut:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max_bin, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
                h28_v2[p] = new TH2F(Form("h28_v2_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. Cluster+Asymm+pos res:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max_bin, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
                h29_v2[p] = new TH2F(Form("h29_v2_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. Asymm+Pos res:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max_bin, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
                h30[p] = new TH2F(Form("h30_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. Blair's cuts, no pos.res no occupancy:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max_bin, MassNBins, 0, smeared_upper_bin_limit);
                h30_1d[p] = new TH1F(Form("h30_1d_%i", p), Form("Smeared Inv Mass, weighted. Blair's cuts, no pos.res no occupancy:%s", WeightNames[p].c_str()), MassNBins, 0, smeared_upper_bin_limit);
                h31[p] = new TH2F(Form("h31_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. Blair's cuts+pos res:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max_bin, MassNBins, 0, smeared_upper_bin_limit);
                h31_1d[p] = new TH1F(Form("h31_1d_%i", p), Form("Smeared Inv Mass, weighted. Blair's cuts+pos res:%s", WeightNames[p].c_str()), MassNBins, 0, smeared_upper_bin_limit);
                h31_pionspectrum[p] = new TH1F(Form("h31_ps_%i", p), Form("Cuts+pos res, Smeared PT, weighted:%s", WeightNames[p].c_str()), MassNBins, PT_Min, PT_Max_bin);
                h34[p] = new TH2F(Form("h34_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. Position Smearing:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max_bin, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
                h34_1d[p] = new TH1F(Form("h34_1d_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. Position Smearing:%s", WeightNames[p].c_str()), MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
                h35[p] = new TH2F(Form("h35_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. Blair's cuts+cluster:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max_bin, MassNBins, 0, smeared_upper_bin_limit);
                h35_1d[p] = new TH1F(Form("h35_1d_%i", p), Form("Smeared Inv Mass, weighted. Blair's cuts+cluster:%s", WeightNames[p].c_str()), MassNBins, 0, smeared_upper_bin_limit);
            }

            h100[p] = new TH2F(Form("h100_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. All Cuts+effects:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max_bin, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
            h100_1d[p] = new TH1F(Form("h100_1d_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. All Cuts+effects:%s", WeightNames[p].c_str()), MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
            h101[p] = new TH2F(Form("h101_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. Everything+eT cuts:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max_bin, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
            h101_1d[p] = new TH1F(Form("h101_1d_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. Everything+eT cuts:%s", WeightNames[p].c_str()), MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
            h101_dr[p] = new TH1F(Form("h101_dr_%i", p), Form("dR distribution, weighted. Everything+eT cuts:%s", WeightNames[p].c_str()), 10000, 0, 2);
            h101_photon_dist_1d[p] = new TH1F(Form("h101_photon_dist_1d_%i", p), Form("Photon distance distribution, weighted. Everything+eT cuts:%s", WeightNames[p].c_str()), 10000, 0, 2000);
            h101_photon_dist[p] = new TH2F(Form("h101_photon_dist_%i", p), Form("pT vs Photon distance distribution, weighted. Everything+eT cuts:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max_bin, 10000, 0, 2000);
            for (int c = 0; c < ClusterScalingNames.size(); c++)
            {
                h101_asymm[p][c] = new TH2F(Form("h101_%i_asymm_%i", p, c), Form("More asymm(%s):Pt vs Inv Mass, weighted, eT cuts:%s", ClusterScalingNames[c].c_str(), WeightNames[p].c_str()), n_bins, 0, PT_Max_bin, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
                h101_symm[p][c] = new TH2F(Form("h101_%i_symm_%i", p, c), Form("More symm(%s):Pt vs Inv Mass, weighted, eT cuts:%s", ClusterScalingNames[c].c_str(), WeightNames[p].c_str()), n_bins, 0, PT_Max_bin, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
            }
            h102[p] = new TH2F(Form("h102_%i", p), Form("Smeared Pt vs Smeared Inv Mass, weighted. same as 100+towerloss:%s", WeightNames[p].c_str()), n_bins, 0, PT_Max_bin, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
        }

        std::random_device rd;
        std::random_device rdgamma;
        std::random_device rdgammacluster;
        std::random_device rdgammapositsmr;
        std::random_device rdclusterloss;
        std::mt19937_64 gen(rd());
        std::mt19937_64 gen_gamma(rdgamma());
        std::mt19937_64 gen_gammacluster(rdgammacluster());
        std::mt19937_64 gen_gammapositsmear(rdgammapositsmr());
        std::mt19937_64 gen_clusterloss(rdclusterloss());
        std::normal_distribution<double> gammadis(0.0, 1.0); // mean 0 and std dev 1
        std::uniform_real_distribution<> gammacluster(0, 1.0);
        std::normal_distribution<double> gamma_positsmear(0.0, 1.0);
        std::uniform_real_distribution<> pdis(PT_ratio, 1.0);
        std::uniform_real_distribution<> adis(0.0, 2 * M_PI);
        std::uniform_real_distribution<> clusterloss(0, 0.1);

        Pythia pythia;
        pythia.readString("PromptPhoton:all = on");
        if (particleType == "Pion")
        {
            pythia.readString("111:oneChannel = 1 1.0 0 22 22");
        }
        else if (particleType == "Eta")
        {
            pythia.readString("221:oneChannel = 1 1.0 0 22 22");
        }
        pythia.init();
        pythia.event.clear();

        int particleCount = 0;
        double particleMass = (particleType == "Pion") ? 0.1349768 : 0.54786;
        double mT_Scaling = 1.0;

        for (int i = 0; i < nParticles; i++)
        {
            double azimuthal_ang = adis(gen);
            double Pt = PT_Max * pdis(gen);

            double weight_function[4];
            double inv_yield[4];
            int WeightScale[4] = {1e+14, 1e+5, 1e+14, 1e+14};

            if (particleType == "Eta")
            {
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

            if (pythia.event.append((particleType == "Pion") ? 111 : 221, 23, 0, 0, particle_px, particle_py, 0.0, particle_E, particleMass))
            {
                particleCount++;
            }
            else
            {
                std::cout << "Failed to append particle at iteration " << i << std::endl;
            }
        }

        pythia.moreDecays();

        for (int i = 0; i < pythia.event.size(); i++)
        {
            if (pythia.event[i].id() == ((particleType == "Pion") ? 111 : 221))
            {
                int Gamma_daughters[2] = {pythia.event[i].daughter1(), pythia.event[i].daughter2()};
                double Pt = pythia.event[i].pT();

                if (particleType == "Eta")
                {
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

                if (pythia.event[Gamma_daughters[0]].id() == 22 && pythia.event[Gamma_daughters[1]].id() == 22)
                {
                    Pythia8::Vec4 gamma_lorentz[3];
                    Pythia8::Vec4 gamma_smeared[3];
                    Pythia8::Vec4 gamma_cluster[3];
                    Pythia8::Vec4 gamma_cluster_asymm[3];
                    Pythia8::Vec4 gamma_Blair_Cuts[3];
                    Pythia8::Vec4 gamma_All_Cuts[3];
                    Pythia8::Vec4 gamma_position_smear[3];
                    Pythia8::Vec4 gamma_Blair_position[3];
                    Pythia8::Vec4 gamma_All_Towerthreshold[3];

                    gamma_lorentz[0] = pythia.event[Gamma_daughters[0]].p();
                    gamma_lorentz[1] = pythia.event[Gamma_daughters[1]].p();
                    gamma_lorentz[2] = gamma_lorentz[0] + gamma_lorentz[1];
                    double inv_mass = gamma_lorentz[2].mCalc();
                    double truthphotondistance = DetectorPhotonDistance(gamma_lorentz[0], gamma_lorentz[1], Debug);

                    double scale_factor1 = sqrt(pow(smear_factor_b, 2) / gamma_lorentz[0].e() + pow(smear_factor_c, 2) + pow(smear_factor_d, 2));
                    double scale_factor2 = sqrt(pow(smear_factor_b, 2) / gamma_lorentz[1].e() + pow(smear_factor_c, 2) + pow(smear_factor_d, 2));

                    double smear_factor1 = scale_factor1 * gammadis(gen_gamma) + 1;
                    double smear_factor2 = scale_factor2 * gammadis(gen_gamma) + 1;

                    gamma_smeared[0] = smear_factor1 * pythia.event[Gamma_daughters[0]].p();
                    gamma_smeared[1] = smear_factor2 * pythia.event[Gamma_daughters[1]].p();

                    for (int photclust = 0; photclust < 2; photclust++)
                    {
                        if (gammacluster(gen_gammacluster) > clusterOverlapProb && clusterOverlap == 1)
                        {
                            gamma_cluster[photclust] = clusterPhoton(gamma_smeared[photclust], 2, myFunc->GetRandom());
                            gamma_cluster_asymm[photclust] = gamma_cluster[photclust];
                        }
                        else
                        {
                            gamma_cluster[photclust] = gamma_smeared[photclust];
                            gamma_cluster_asymm[photclust] = gamma_cluster[photclust];
                        }
                    }

                    gamma_All_Cuts[0] = PositionResSmear(gamma_cluster_asymm[0], posit_smearingFactor * gamma_positsmear(gen_gammapositsmear), posit_smearingFactor * gamma_positsmear(gen_gammapositsmear), posit_smearingFactor * gamma_positsmear(gen_gammapositsmear));
                    gamma_All_Cuts[1] = PositionResSmear(gamma_cluster_asymm[1], posit_smearingFactor * gamma_positsmear(gen_gammapositsmear), posit_smearingFactor * gamma_positsmear(gen_gammapositsmear), posit_smearingFactor * gamma_positsmear(gen_gammapositsmear));
                    gamma_All_Cuts[2] = gamma_All_Cuts[0] + gamma_All_Cuts[1];
                    gamma_All_Towerthreshold[0] = ClusterLossSimple(gamma_All_Cuts[0], clusterloss(gen_clusterloss),tower_cluster_threshold);
                    gamma_All_Towerthreshold[1] = ClusterLossSimple(gamma_All_Cuts[1], clusterloss(gen_clusterloss),tower_cluster_threshold);
                    gamma_All_Towerthreshold[2] = gamma_All_Towerthreshold[0] + gamma_All_Towerthreshold[1];

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
                    if (Debug_Hists)
                    {
                        h4->Fill(Pt);
                        h6->Fill(inv_mass);
                        h8->Fill(inv_mass_smeared);
                        h9->Fill(gamma_smeared[2].pT(), inv_mass_smeared);
                        h10->Fill(gamma_smeared[2].pT());
                        h17->Fill(gamma_lorentz[0].pT());
                        h17->Fill(gamma_lorentz[1].pT());
                        h16->Fill(gamma_smeared[0].pT());
                        h16->Fill(gamma_smeared[1].pT());
                    }

                    for (int p = 0; p < WeightNames.size(); p++)
                    {

                        htruthphotondistance_1d[p]->Fill(truthphotondistance, inv_yield[p]);
                        htruthphotondistance[p]->Fill(gamma_smeared[2].pT(), truthphotondistance, inv_yield[p]);
                        if (Debug_Hists)
                        {
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
                            if (AsymmCutcheck(gamma_smeared[0], gamma_smeared[1], asymmCutValue, applyAsymmCut) == true)
                            {
                                h29[p]->Fill(gamma_smeared[2].pT(), gamma_smeared[2].mCalc(), inv_yield[p]);
                            }
                            if (AsymmCutcheck(gamma_cluster_asymm[0], gamma_cluster_asymm[1], asymmCutValue, applyAsymmCut) == true)
                            {
                                h28[p]->Fill(gamma_cluster_asymm[2].pT(), gamma_cluster_asymm[2].mCalc(), inv_yield[p]);
                            }
                            if (AsymmCutcheck(gamma_position_smear[0], gamma_position_smear[1], asymmCutValue, applyAsymmCut) == true)
                            {
                                h29_v2[p]->Fill(gamma_position_smear[2].pT(), gamma_position_smear[2].mCalc(), inv_yield[p]);
                            }
                            if (AsymmCutcheck(gamma_All_Cuts[0], gamma_All_Cuts[1], asymmCutValue, applyAsymmCut) == true)
                            {
                                h28_v2[p]->Fill(gamma_All_Cuts[2].pT(), gamma_All_Cuts[2].mCalc(), inv_yield[p]);
                            }
                            if (DeltaRcut(gamma_smeared[0], gamma_smeared[1], DeltaRcut_MAX) == false && AsymmCutcheck(gamma_smeared[0], gamma_smeared[1], asymmCutValue, applyAsymmCut) == true && pTCut(gamma_smeared[0], pt1cut) == true && pTCut(gamma_smeared[1], pt2cut) == true && nclus_ptCut < gamma_smeared[0].pT() && gamma_smeared[0].pT() < ptMaxCut && nclus_ptCut < gamma_smeared[1].pT() && gamma_smeared[1].pT() < ptMaxCut && gamma_smeared[2].pT() > comb_ptcut * (pt1cut + pt2cut))
                            {
                                h30[p]->Fill(gamma_smeared[2].pT(), gamma_smeared[2].mCalc(), inv_yield[p]);
                                h30_1d[p]->Fill(gamma_smeared[2].mCalc(), inv_yield[p]);
                            }

                            if (DeltaRcut(gamma_Blair_position[0], gamma_Blair_position[1], DeltaRcut_MAX) == false && AsymmCutcheck(gamma_Blair_position[0], gamma_Blair_position[1], asymmCutValue, applyAsymmCut) == true && pTCut(gamma_Blair_position[0], pt1cut) == true && pTCut(gamma_Blair_position[1], pt2cut) == true && nclus_ptCut < gamma_Blair_position[0].pT() && gamma_Blair_position[0].pT() < ptMaxCut && nclus_ptCut < gamma_Blair_position[1].pT() && gamma_Blair_position[1].pT() < ptMaxCut && gamma_Blair_position[2].pT() > comb_ptcut * (pt1cut + pt2cut))
                            {
                                h31[p]->Fill(gamma_Blair_position[2].pT(), gamma_Blair_position[2].mCalc(), inv_yield[p]);
                                h31_1d[p]->Fill(gamma_Blair_position[2].mCalc(), inv_yield[p]);
                                h31_pionspectrum[p]->Fill(gamma_Blair_position[2].pT(), inv_yield[p]);
                                hpionpt[p]->Fill(gamma_lorentz[2].pT());
                                hInvMass_Cutson->Fill(gamma_Blair_position[2].mCalc());
                            }

                            if (DeltaRcut(gamma_Blair_Cuts[0], gamma_Blair_Cuts[1], DeltaRcut_MAX) == false && AsymmCutcheck(gamma_Blair_Cuts[0], gamma_Blair_Cuts[1], asymmCutValue, applyAsymmCut) == true && pTCut(gamma_Blair_Cuts[0], pt1cut) == true && pTCut(gamma_Blair_Cuts[1], pt2cut) == true && nclus_ptCut < gamma_Blair_Cuts[0].pT() && gamma_Blair_Cuts[0].pT() < ptMaxCut && nclus_ptCut < gamma_Blair_Cuts[1].pT() && gamma_Blair_Cuts[1].pT() < ptMaxCut && gamma_Blair_Cuts[2].pT() > comb_ptcut * (pt1cut + pt2cut))
                            {
                                h35[p]->Fill(gamma_Blair_Cuts[2].pT(), gamma_Blair_Cuts[2].mCalc(), inv_yield[p]);
                                h35_1d[p]->Fill(gamma_Blair_Cuts[2].mCalc(), inv_yield[p]);
                            }
                        }

                        if (DeltaRcut(gamma_All_Cuts[0], gamma_All_Cuts[1], DeltaRcut_MAX) == false && AsymmCutcheck(gamma_All_Cuts[0], gamma_All_Cuts[1], asymmCutValue, applyAsymmCut) == true && pTCut(gamma_All_Cuts[0], pt1cut) == true && pTCut(gamma_All_Cuts[1], pt2cut) == true && nclus_ptCut < gamma_All_Cuts[0].pT() && gamma_All_Cuts[0].pT() < ptMaxCut && nclus_ptCut < gamma_All_Cuts[1].pT() && gamma_All_Cuts[1].pT() < ptMaxCut && gamma_All_Cuts[2].pT() > comb_ptcut * (pt1cut + pt2cut))
                        {
                            h100[p]->Fill(gamma_All_Cuts[2].pT(), gamma_All_Cuts[2].mCalc(), inv_yield[p]);
                            h100_1d[p]->Fill(gamma_All_Cuts[2].mCalc(), inv_yield[p]);
                        }

                        if(DeltaRcut(gamma_All_Towerthreshold[0], gamma_All_Towerthreshold[1], DeltaRcut_MAX) == false && AsymmCutcheck(gamma_All_Towerthreshold[0], gamma_All_Towerthreshold[1], asymmCutValue, applyAsymmCut) == true && pTCut(gamma_All_Towerthreshold[0], pt1cut) == true && pTCut(gamma_All_Towerthreshold[1], pt2cut) == true && nclus_ptCut < gamma_All_Towerthreshold[0].pT() && gamma_All_Towerthreshold[0].pT() < ptMaxCut && nclus_ptCut < gamma_All_Towerthreshold[1].pT() && gamma_All_Towerthreshold[1].pT() < ptMaxCut && gamma_All_Towerthreshold[2].pT() > comb_ptcut * (pt1cut + pt2cut))
                        {
                            h102[p]->Fill(gamma_All_Towerthreshold[2].pT(), gamma_All_Towerthreshold[2].mCalc(), inv_yield[p]);
                        }

                        if (DeltaRcut(gamma_All_Cuts[0], gamma_All_Cuts[1], DeltaRcut_MAX) == false && AsymmCutcheck(gamma_All_Cuts[0], gamma_All_Cuts[1], asymmCutValue, applyAsymmCut) == true && eTCut(gamma_All_Cuts[0], etCut) == true && eTCut(gamma_All_Cuts[1], etCut) == true && nclus_ptCut < gamma_All_Cuts[0].pT() && gamma_All_Cuts[0].pT() < ptMaxCut && nclus_ptCut < gamma_All_Cuts[1].pT() && gamma_All_Cuts[1].pT() < ptMaxCut && gamma_All_Cuts[2].pT() > comb_ptcut * (pt1cut + pt2cut))
                        {
                            h101[p]->Fill(gamma_All_Cuts[2].pT(), gamma_All_Cuts[2].mCalc(), inv_yield[p]);
                            h101_1d[p]->Fill(gamma_All_Cuts[2].mCalc(), inv_yield[p]);
                            h101_dr[p]->Fill(DeltaR(gamma_All_Cuts[0], gamma_All_Cuts[1]), inv_yield[p]);
                        }
                        // fill photon distance hist
                        double photon_dist_All_Cuts = DetectorPhotonDistance(gamma_All_Cuts[0], gamma_All_Cuts[1], Debug);
                        h101_photon_dist_1d[p]->Fill(photon_dist_All_Cuts, inv_yield[p]);
                        h101_photon_dist[p]->Fill(gamma_All_Cuts[2].pT(), photon_dist_All_Cuts, inv_yield[p]);
                        // clustering algorithm check
                        for (int c = 0; c < ClusterScalingNames.size(); c++)
                        {
                            auto [symmetricPhoton1, symmetricPhoton2] = adjustPhotonEnergiesSymmetric(gamma_All_Cuts[0], gamma_All_Cuts[1], c, Debug);
                            auto [asymmetricPhoton1, asymmetricPhoton2] = adjustPhotonEnergiesAsymmetric(gamma_All_Cuts[0], gamma_All_Cuts[1], c, Debug);
                            auto asymmetricPion = asymmetricPhoton1 + asymmetricPhoton2;
                            auto symmetricPion = symmetricPhoton1 + symmetricPhoton2;

                            if (DeltaRcut(asymmetricPhoton1, asymmetricPhoton2, DeltaRcut_MAX) == false && AsymmCutcheck(asymmetricPhoton1, asymmetricPhoton2, asymmCutValue, applyAsymmCut) == true && eTCut(asymmetricPhoton1, etCut) == true && eTCut(asymmetricPhoton2, etCut) == true && nclus_ptCut < asymmetricPhoton1.pT() && asymmetricPhoton1.pT() < ptMaxCut && nclus_ptCut < asymmetricPhoton2.pT() && asymmetricPhoton2.pT() < ptMaxCut && asymmetricPhoton1.pT() > pt1cut && asymmetricPhoton2.pT() > pt2cut && asymmetricPhoton1.pT() + asymmetricPhoton2.pT() > comb_ptcut * (pt1cut + pt2cut))
                            {
                                h101_asymm[p][c]->Fill(asymmetricPion.pT(), asymmetricPion.mCalc(), inv_yield[p]);
                            }
                            if (DeltaRcut(symmetricPhoton1, symmetricPhoton2, DeltaRcut_MAX) == false && AsymmCutcheck(symmetricPhoton1, symmetricPhoton2, asymmCutValue, applyAsymmCut) == true && eTCut(symmetricPhoton1, etCut) == true && eTCut(symmetricPhoton2, etCut) == true && nclus_ptCut < symmetricPhoton1.pT() && symmetricPhoton1.pT() < ptMaxCut && nclus_ptCut < symmetricPhoton2.pT() && symmetricPhoton2.pT() < ptMaxCut && symmetricPhoton1.pT() > pt1cut && symmetricPhoton2.pT() > pt2cut && symmetricPhoton1.pT() + symmetricPhoton2.pT() > comb_ptcut * (pt1cut + pt2cut))
                            {
                                h101_symm[p][c]->Fill(symmetricPion.pT(), symmetricPion.mCalc(), inv_yield[p]);
                            }
                        }
                    }

                    if (Debug_Hists)
                    {
                        try
                        {
                            mass_pt_map[floor(Pt) + 1].push_back(inv_mass_smeared);
                        }
                        catch (...)
                        {
                            mass_pt_map.insert({floor(Pt) + 1, std::vector<double>()});
                        }

                        for (int j : Gamma_daughters)
                        {
                            h5->Fill(pythia.event[j].pT());
                        }
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

        if (saveToTree)
        {
            delete tree;
        }
        delete output;
    }
    std::cout << "All done!" << std::endl; // debug line
    return 0;
}

TF1 *ChooseSpectrumFunction(int weightmethod, int PT_Min, int PT_Max, const std::string &particleType)
{
    TF1 *myFunc = nullptr;
    if (particleType == "Pion")
    {
        if (weightmethod == 0)
        {
            myFunc = new TF1("myFunc", [](double *x, double *par)
                             { return x[0] * par[0] * TMath::Exp(-x[0] / 0.3); }, PT_Min, PT_Max, 1);
            Double_t initialParameters[1] = {1.0};
            myFunc->SetParameters(initialParameters);
        }
        else if (weightmethod == 1)
        {
            myFunc = new TF1("myFunc", [](double *x, double *par)
                             { return x[0] * par[0] * pow(x[0], -8.14); }, PT_Min, PT_Max, 1);
            Double_t initialParameters[1] = {1.0};
            myFunc->SetParameters(initialParameters);
        }
        else if (weightmethod == 2)
        {
            myFunc = new TF1("myFunc", [](double *x, double *par)
                             { return x[0] * ((1 / (1 + exp((x[0] - par[0]) / par[1]))) * par[2] / pow(1 + x[0] / par[3], par[4]) + (1 - (1 / (1 + exp((x[0] - par[0]) / par[1])))) * par[5] / (pow(x[0], par[6]))); }, PT_Min, PT_Max, 7);
            myFunc->SetParameters(4.5, 0.114, 229.6, 1.466, 10.654, 14.43, 8.1028);
        }
        else if (weightmethod == 3)
        {
            myFunc = new TF1("myFunc", [](double *x, double *par)
                             { return x[0] * (par[0] / pow(1 + x[0] / par[1], par[2])); }, PT_Min, PT_Max, 3);
            myFunc->SetParameters(229.6, 1.466, 10.654);
        }
        else
        {
            std::cout << "Error: No Weight function found" << std::endl;
        }
    }
    else if (particleType == "Eta")
    {
        if (weightmethod == 0)
        {
            myFunc = new TF1("myFunc", [](double *x, double *par)
                             { return x[0] * par[1] * pow((par[2] + sqrt(pow(par[3], 2) + pow(x[0], 2))) / (par[2] + sqrt(pow(par[4], 2) + pow(x[0], 2))), -par[5]) * par[0] * TMath::Exp(-x[0] / 0.3); }, PT_Min, PT_Max, 6);
            myFunc->SetParameters(1.0, 0.5, 1.2, 0.54786, 0.1349768, 10);
        }
        else if (weightmethod == 1)
        {
            myFunc = new TF1("myFunc", [](double *x, double *par)
                             { return x[0] * par[1] * pow((par[2] + sqrt(pow(par[3], 2) + pow(x[0], 2))) / (par[2] + sqrt(pow(par[4], 2) + pow(x[0], 2))), -par[5]) * par[0] * pow(x[0], -8.14); }, PT_Min, PT_Max, 6);
            myFunc->SetParameters(1.0, 0.5, 1.2, 0.54786, 0.1349768, 10);
        }
        else if (weightmethod == 2)
        {
            myFunc = new TF1("myFunc", [](double *x, double *par)
                             { return x[0] * par[7] * pow((par[8] + sqrt(pow(par[9], 2) + pow(x[0], 2))) / (par[8] + sqrt(pow(par[10], 2) + pow(x[0], 2))), -par[11]) * ((1 / (1 + exp((x[0] - par[0]) / par[1]))) * par[2] / pow(1 + x[0] / par[3], par[4]) + (1 - (1 / (1 + exp((x[0] - par[0]) / par[1])))) * par[5] / (pow(x[0], par[6]))); }, PT_Min, PT_Max, 12);
            Double_t params[] = {4.5, 0.114, 229.6, 1.466, 10.654, 14.43, 8.1028, 0.5, 1.2, 0.54786, 0.1349768, 10.0};
            myFunc->SetParameters(params);
        }
        else if (weightmethod == 3)
        {
            myFunc = new TF1("myFunc", [](double *x, double *par)
                             { return x[0] * par[3] * pow((par[4] + sqrt(pow(par[5], 2) + pow(x[0], 2))) / (par[4] + sqrt(pow(par[6], 2) + pow(x[0], 2))), -par[7]) * (par[0] / pow(1 + x[0] / par[1], par[2])); }, PT_Min, PT_Max, 8);
            myFunc->SetParameters(229.6, 1.466, 10.654, 0.5, 1.2, 0.54786, 0.1349768, 10);
        }
        else
        {
            std::cout << "Error: No Weight function found" << std::endl;
        }
    }
    return myFunc;
}
// clustering and psotion smearing functions
Pythia8::Vec4 clusterPhoton(Pythia8::Vec4 &originalPhoton, int method, double randomE)
{
    Pythia8::Vec4 newPhoton;
    Pythia8::Rndm rndm;
    if (method == 1)
    {
        newPhoton.e(randomE);
        newPhoton.px(originalPhoton.px() * (randomE / originalPhoton.e()));
        newPhoton.py(originalPhoton.py() * (randomE / originalPhoton.e()));
        newPhoton.pz(originalPhoton.pz() * (randomE / originalPhoton.e()));
    }
    if (method == 2)
    {
        // Set the energy
        newPhoton.e(randomE);

        // Scale the momentum components according to the new energy
        double scale = randomE / originalPhoton.e();
        newPhoton.px(originalPhoton.px() * scale);
        newPhoton.py(originalPhoton.py() * scale);
        newPhoton.pz(originalPhoton.pz() * scale);

        // Explicitly enforce the massless condition: E^2 = px^2 + py^2 + pz^2
        double newE = sqrt(newPhoton.px() * newPhoton.px() +
                           newPhoton.py() * newPhoton.py() +
                           newPhoton.pz() * newPhoton.pz());

        // Adjust the energy to match the massless condition
        newPhoton.e(newE);
    }
    // Check for any residual mass due to numerical issues
    if (fabs(newPhoton.mCalc()) > 1e-6)
    {
        // Allowing a tiny tolerance for floating-point precision
        std::cerr << "Error: Clustered photon has mass" << std::endl;
    }
    return newPhoton + originalPhoton;
}

Pythia8::Vec4 PositionResSmear(Pythia8::Vec4 photon, double smearingFactorx, double smearingFactory, double smearingFactorz)
{
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

// cut functions
bool DeltaRcut(Pythia8::Vec4 &Photon1, Pythia8::Vec4 &Photon2, float DeltaRcutMax)
{
    double dEta = Photon1.eta() - Photon2.eta();
    double dPhi = acos(cos(Photon1.phi() - Photon2.phi()));
    double deltaR = sqrt(dEta * dEta + dPhi * dPhi);
    return deltaR > DeltaRcutMax;
}

float DeltaR(Pythia8::Vec4 &Photon1, Pythia8::Vec4 &Photon2)
{
    double dEta = Photon1.eta() - Photon2.eta();
    double dPhi = acos(cos(Photon1.phi() - Photon2.phi()));
    return sqrt(dEta * dEta + dPhi * dPhi);
}

bool pTCut(const Pythia8::Vec4 &particle, float ptCut)
{
    double pT = particle.pT();
    return pT > ptCut; // true if greater than the cut. I use this for the min and max pT cuts
}

bool eTCut(const Pythia8::Vec4 &particle, float eTCut)
{
    double eT = sqrt(pow(particle.mCalc(), 2) + particle.pT2());
    return eT > eTCut; // true if greater than the cut.
}

bool EtaCut(const Pythia8::Vec4 &particle, float EtaCutValue, bool ApplyEtaCut, bool debug)
{
    if (ApplyEtaCut == false) // will fill the histogram if the cut is not applied
    {
        return false;
    }
    float Eta = abs(particle.eta());
    if (debug)
        std::cout << "Debug:  " << " Eta, abs(Eta) = " << particle.eta() << " , " << Eta << std::endl;
    return Eta > EtaCutValue; // true if greater than the cut, I.e. out of the allowed range
}

bool AsymmCutcheck(Pythia8::Vec4 &Photon1, Pythia8::Vec4 &Photon2, float AsymmCutoff, bool asymcutbool)
{
    if (asymcutbool == false)
    {
        return true;
    }
    return abs(Photon1.e() - Photon2.e()) / (Photon1.e() + Photon2.e()) < AsymmCutoff;
}

// Function to parse command-line arguments
void parseArguments(int argc, char *argv[], std::map<std::string, std::string> &params, bool debug)
{
    for (int i = 1; i < argc; i++)
    {
        std::string arg = argv[i];
        size_t equalPos = arg.find('=');
        if (equalPos != std::string::npos)
        {
            std::string key = arg.substr(1, equalPos - 1); // Remove the leading '-'
            std::string value = arg.substr(equalPos + 1);
            params[key] = value;
        }
    }

    // Debugging: Print all parsed parameters
    if (debug)
    {
        std::cout << "\n Parsed command-line parameters: \n"
                  << std::endl;
        for (const auto &param : params)
        {
            std::cout << param.first << " = " << param.second << std::endl;
        }
    }
}

// Functions to calculate the scaling factor for the energy sharing between two photons
// Function to calculate the Distance between two particles on the projected surface of the emcal
double DetectorPhotonDistance(Pythia8::Vec4 &photon1, Pythia8::Vec4 &photon2, bool debug)
{
    // Calculate the projected distance between the two photons
    // photon 1
    double energy1 = photon1.e();
    double posx1 = 900 * photon1.px() / photon1.e();
    double posy1 = 900 * photon1.py() / photon1.e();
    double posz1 = photon1.pz() / photon1.e();
    // photon 2
    double energy2 = photon2.e();
    double posx2 = 900 * photon2.px() / photon2.e();
    double posy2 = 900 * photon2.py() / photon2.e();
    double posz2 = photon2.pz() / photon2.e();
    double dz = posz1 - posz2;

    if (debug && dz < 0.1)
    {
        std::cerr << "Warning: Photons are too close in z direction: " << dz << std::endl;
    }

    double distance = sqrt(pow(posx1 - posx2, 2) + pow(posy1 - posy2, 2) + pow(dz, 2));

    return distance;
}

std::pair<Pythia8::Vec4, Pythia8::Vec4> adjustPhotonEnergiesSymmetric(Pythia8::Vec4 photon1, Pythia8::Vec4 photon2, int method_int, bool debug)
{
    double tolerance = 1e-6;
    double distance = DetectorPhotonDistance(photon1, photon2, debug);

    // three choices, inverse exponential, hyberbolic tangent, or rational
    ScalingMethod method; // = EXPONENTIAL;
    switch (method_int)
    {
    case 0:
        method = EXPONENTIAL;
        break;
    case 1:
        method = RATIONAL;
        break;
    case 2:
        method = HYPERBOLIC_TANGENT;
        break;
    default:
        std::cerr << "Error: Unknown scaling method." << std::endl;
        break;
    }

    double shiftFactor = EnergySharingScale(distance, method);

    double totalEnergy = photon1.e() + photon2.e();
    double avgEnergy = totalEnergy / 2.0;

    double originalPhoton1E = photon1.e();
    double originalPhoton2E = photon2.e();

    photon1.e((photon1.e() * (1.0 - shiftFactor)) + (avgEnergy * shiftFactor));
    photon2.e((photon2.e() * (1.0 - shiftFactor)) + (avgEnergy * shiftFactor));

    // Scale the momentum components according to the new energy
    double scale1 = photon1.e() / originalPhoton1E;
    double scale2 = photon2.e() / originalPhoton2E;

    photon1.px(photon1.px() * scale1);
    photon1.py(photon1.py() * scale1);
    photon1.pz(photon1.pz() * scale1);

    photon2.px(photon2.px() * scale2);
    photon2.py(photon2.py() * scale2);
    photon2.pz(photon2.pz() * scale2);

    if (debug)
    {
        if (std::abs(photon1.mCalc()) > tolerance)
        {
            std::cerr << "Warning: Photon1 has non-zero mass after scaling: " << photon1.mCalc() << std::endl;
            std::cerr << "Photon1: " << photon1.e() << " " << photon1.px() << " " << photon1.py() << " " << photon1.pz() << std::endl;
        }
        if (std::abs(photon2.mCalc()) > tolerance)
        {
            std::cerr << "Warning: Photon2 has non-zero mass after scaling: " << photon2.mCalc() << std::endl;
            std::cerr << "Photon2: " << photon2.e() << " " << photon2.px() << " " << photon2.py() << " " << photon2.pz() << std::endl;
        }
    }
    return std::make_pair(photon1, photon2);
}

std::pair<Pythia8::Vec4, Pythia8::Vec4> adjustPhotonEnergiesAsymmetric(Pythia8::Vec4 photon1, Pythia8::Vec4 photon2, int method_int, bool debug)
{
    double tolerance = 1e-6;
    double distance = DetectorPhotonDistance(photon1, photon2, debug);

    double totalEnergy = photon1.e() + photon2.e();
    // three choices, inverse exponential, hyberbolic tangent, or rational
    ScalingMethod method; // = EXPONENTIAL;
    switch (method_int)
    {
    case 0:
        method = EXPONENTIAL;
        break;
    case 1:
        method = RATIONAL;
        break;
    case 2:
        method = HYPERBOLIC_TANGENT;
        break;
    default:
        std::cerr << "Error: Unknown scaling method." << std::endl;
        break;
    }

    double shiftFactor = EnergySharingScale(distance, method);
    double energyShift = (photon1.e() - photon2.e()) * shiftFactor;

    // Store original energies for scaling momentum
    double originalPhoton1E = photon1.e();
    double originalPhoton2E = photon2.e();

    photon1.e(photon1.e() + energyShift);
    photon2.e(photon2.e() - energyShift);

    // Scale the momentum components according to the new energy
    double scale1 = photon1.e() / originalPhoton1E;
    double scale2 = photon2.e() / originalPhoton2E;

    photon1.px(photon1.px() * scale1);
    photon1.py(photon1.py() * scale1);
    photon1.pz(photon1.pz() * scale1);

    photon2.px(photon2.px() * scale2);
    photon2.py(photon2.py() * scale2);
    photon2.pz(photon2.pz() * scale2);

    if (debug)
    {
        if (std::abs(photon1.mCalc()) > tolerance)
        {
            std::cerr << "Warning: Photon1 has non-zero mass after scaling: " << photon1.mCalc() << std::endl;
            std::cerr << "Photon1: " << photon1.e() << " " << photon1.px() << " " << photon1.py() << " " << photon1.pz() << std::endl;
        }
        if (std::abs(photon2.mCalc()) > tolerance)
        {
            std::cerr << "Warning: Photon2 has non-zero mass after scaling: " << photon2.mCalc() << std::endl;
            std::cerr << "Photon2: " << photon2.e() << " " << photon2.px() << " " << photon2.py() << " " << photon2.pz() << std::endl;
        }
    }

    return std::make_pair(photon1, photon2);
}

double EnergySharingScale(double distance, ScalingMethod method)
{
    double scale = 0.0;

    switch (method)
    {
    case EXPONENTIAL:
    {
        // 0.3*exp(-distance / 0.2); // exponential
        // parameters for the Exponential method
        double expParam1 = 0.7;  // finite energy sharing at low distances
        double expParam2 = 50.0; // distance at which scale = expParam*1/e~expParam*0.37
        // higher=higher mass?
        //  Exponential scaling: scale = expParam1 * exp(-distance / expParam2)
        scale = expParam1 * exp(-distance / expParam2);
        break;
    }

    case RATIONAL:
    {
        // double shiftFactor = A / (1.0 + pow(distance / d0, n)); // rational
        // parameters for the Rational method
        double rationalParam1 = 0.7;  // finite energy sharing at low distances
        double rationalParam2 = 50.0; // distance scaling factor. no idea what this should be
        double rationalParam3 = 1.3;  // A parameter that controls the steepness of the decline
        // Rational scaling: scale = rationalParam1 / (distance + rationalParam2)
        scale = rationalParam1 / (1 + pow(distance / rationalParam2, rationalParam3));
        break;
    }

    case HYPERBOLIC_TANGENT:
    {
        // double shiftFactor = A * tanh((distance - d0) / d1); // hyperbolic tangent
        // parameters for the Hyperbolic Tangent method
        double tanhParam1 = 0.7;  // finite energy sharing at low distances
        double tanhParam2 = 10.0; // sets the scale of the distance at which the energy sharing is = param1. lower = less sharing at high distances
        scale = tanhParam1 * tanh(tanhParam2 / distance);
        break;
    }

    default:
        std::cerr << "Error: Unknown scaling method." << std::endl;
        break;
    }

    return scale;
}

//function to "lose" energy due to the clusterizer tower threshold
//two versions, one simplified case where a random number is rolled and if that number times the photon energy is below the threshold, the energy is lost from the photon.
Pythia8::Vec4 ClusterLossSimple(Pythia8::Vec4 &originalPhoton, float randomRoll, float threshold)
{
    Pythia8::Vec4 newPhoton = originalPhoton;

    float energyLoss = randomRoll * 0.014;
    newPhoton = originalPhoton - energyLoss;

    return newPhoton;
}
//second version is create a 3x3 grid. the photon lands randomly somewhere in the central cell, and the energy is shared between the central cell and the 8 surrounding cells. Any cell with energy below the threshold is lost. This is tricky because you have to model how much energy each cell gets based on the central position. 
