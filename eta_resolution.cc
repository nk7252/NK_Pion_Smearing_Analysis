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
using namespace Pythia8;	// Let Pythia8:: be implicit.

//forward declarators
TF1* ChooseSpectrumFunction(int weightmethod, int PT_Min, int PT_Max);
Pythia8::Vec4 clusterPhoton(Pythia8::Vec4& originalPhoton, int method, double randomE);
Pythia8::Vec4 PositionResSmear(Pythia8::Vec4 photon, double smearingFactorx,double smearingFactory,double smearingFactorz);
bool DeltaRcut(Pythia8::Vec4& Photon1, Pythia8::Vec4& Photon2, float DeltaRcutMax);
bool pTCut(const Pythia8::Vec4& particle, float ptCut);
bool AsymmCutcheck(Pythia8::Vec4& Photon1, Pythia8::Vec4& Photon2, float AsymmCutoff,bool asymcutbool);

int main(){ 
	TStopwatch timer;
	timer.Start();
	int NPions = 8 * 1000000;//
	// int n_bins=1+ceil(log2(Nevents));
	int PT_Max = 50; // 64 normally, [0.2,10] to compare to geant sim
	float PT_Min = 0;	 // cross check to elimate power law problems
	double PT_ratio = PT_Min / PT_Max;
	int MassNBins =50;//40
	// int n_bins=round((1/4)*PT_Max);
	// int n_bins=round(1+3.222*log(NPions));
	int binres=1;//number of divisions per GeV
	int n_bins = binres*PT_Max;//multiple by bin res.
	std::map<double, std::vector<double>> mass_pt_map; // we want to have keys of a pT range?


	//-----------------------------------set weighting method
	int weightmethod = 2;//0=exp,1=power,2=wshp, 3=hagedorn(not implemented)
	std::vector<std::string> WeightNames = {"EXP", "POWER", "WSHP","HAGEDORN"};
	//-----------------------------------
	//cuts
	bool asymcut=true;//apply asymm cut.
	float asymval= 0.6;
	int clusteroverlay = 1;//overlayed cluster check
	float coprob=0.99;//random numbers(0-1) greater than this value will have some smearing added.

	//Blair specific cuts
	float DeltaRcut_MAX = 1.1; 
	float pt1cut =1.5;   
	float pt2cut =1.5;
    float comb_ptcut = 0;
	float ptMaxCut = 50;
	float nclus_ptCut =0.0;
	//nclus is not something I can do here is it          

	//--------------------Alternative paramaterization, woods saxon+hagedorn+power law
    float R_EP_Inf=0.5;
    float mT_a=1.2;
    int mT_n=-10;
	double t = 4.5;
	double w = 0.114;
	double A = 229.6;
	double B = 14.43;
	double n = 8.1028;
	double m_param = 10.654;
	double p0 = 1.466;
	// static smear 0.04
	// float smeared_lower_bin_limit=0.11;
	// float smeared_upper_bin_limit=0.16;
	float smeared_lower_bin_limit = 0.0;
	float smeared_upper_bin_limit = 1.0;
	float smear_factor_a = 0;
	float smear_factor_d = 0.0;  // 0.02;// test trying to include the beam momentum resolution.. will set to zero for now
	//float smear_factor_c = 0.02; // first parameter in test beam parameterization? 2.8 in test beam.
	float posit_smearingFactor = 2.8; // Example smearing factor for position. use half of phnx pos res(simplified) try 2.8 mm. can set the scale to whatever I want, so I will use mm.

	//std::cout << "Processing: " << WeightNames[weightmethod] << std::endl;
	//----------------------pion spectrum function for clusteroverlay
	// reserve a TF1 for the chosen function just in case
	TF1 *myFunc;
	myFunc=ChooseSpectrumFunction(weightmethod, PT_Min, PT_Max);

	for (int smear_factor_itt = 0; smear_factor_itt < 0 + 1; smear_factor_itt++){// originally int smear_factor_itt = 0; smear_factor_itt < 24 + 1; smear_factor_itt++
		// only want .155
		float smear_factor_basevalue = 0.168; //0.498, I used 1.6% + 12.7%/sqrt(E) fig 22, but that is from a special beam cross section config. trying with fig 24 data i.e 2.8% + 15.5%
		//--------------------preliminaries to read from root
		float smear_factor_c = smear_factor_basevalue + 0.001 * smear_factor_itt;//constant term
		float smear_factor_b = 0.154;//0.154,sqrt(E) term
		//////////////////////New//0.155 loop from twice test beam data paramaterization to half? this is 15.5% from https://arxiv.org/pdf/1704.01461.pdf fig 24b so going from 6.5% to 30.5%// need 24 steps for 1% diff each
		//////////////////////OLD//0.127 loop from twice test beam data paramaterization to half? this is 12.7% from https://arxiv.org/pdf/1704.01461.pdf fig 22b so going from 6.35% to 25.4%
		TFile *output = new TFile(Form("pioncode/rootfiles/EtaFastMC_%f_sqrte_%f_const.root", smear_factor_b , smear_factor_c ), "recreate");//

		//TTree *tree = new TTree("tree", "tree");
		//tree->SetMaxTreeSize(500 * 1024 * 1024); // set max tree size to 500 mb
		//------------------------------ book histograms. Need to cull this list. some are useless/redundant
		// TH1* h1 = new TH1F("h1", "pi0 E",128, 0, PT_Max);

		
		TH1 *h4 = new TH1F("h4", "Pion PT, unweighted", n_bins, PT_Min, PT_Max);
		TH1 *h5 = new TH1F("h5", "Photon Pt, unweighted", n_bins, PT_Min, PT_Max);
		TH1 *h6 = new TH1F("h6", "inv mass of gamma pair", MassNBins, 0, 1);
		//TH1 *h7 = new TH1F("h7", "ratio of Photon/Pion pt", n_bins, PT_Min, PT_Max);
		TH1 *h8 = new TH1F("h8", "inv mass of Photon pair, smeared", MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
		TH2F *h9 = new TH2F("h9", "Smeared Pion Pt vs Smeared Inv Mass", n_bins, 0, PT_Max, MassNBins, smeared_lower_bin_limit,  smeared_upper_bin_limit);
		TH1 *h10 = new TH1F("h10", "Smeared Pion PT", n_bins, PT_Min, PT_Max);
		//TH1 *h11 = new TH1F("h11", "Smeared Pion PT/Pion PT ratio", n_bins, PT_Min, PT_Max);
		
		//TH1 *h13 = new TH1F("h13", "ratio of weighted  Smeared-Pion-PT/ weighted Pion PT ", n_bins, PT_Min, PT_Max);
		//------------------------photon pT
		TH1 *h16 = new TH1F("h16", "Smeared Photon pT", n_bins, PT_Min, PT_Max);
		TH1 *h17 = new TH1F("h17", "Photon pT", n_bins, PT_Min, PT_Max);
		
		//TH1 *h24 = new TH1F("h24", "Photon pT ratio, smeared/unsmeared", n_bins, PT_Min, PT_Max);
		//TH1 *h26 = new TH1F("h26", "weighted, Photon pT ratio, smeared/unsmeared", n_bins, PT_Min, PT_Max);
		TH1 *hInvMass_Cutson = new TH1F("hInvMass_Cutson", "Pion pT,nSmeared+no_weight+cuts+pr", MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);

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

		//std::vector<TH2F*> h32(WeightNames.size());//empty
		//std::vector<TH2F*> h33(WeightNames.size());//empty

		std::vector<TH2F*> h34(WeightNames.size());
		std::vector<TH1F*> h34_1d(WeightNames.size());

		std::vector<TH2F*> h35(WeightNames.size());
		std::vector<TH1F*> h35_1d(WeightNames.size());

		std::vector<TH2F*> h100(WeightNames.size());
		std::vector<TH1F*> h100_1d(WeightNames.size());


		for (int p=0; p < WeightNames.size(); p++){
			hpionpt[p] = new TH1D(Form("hpionpt_%i",p),Form("Pion Pt no smear + no weight:%s",WeightNames[p].c_str()) , n_bins, PT_Min, PT_Max);
			h2[p] = new TH1D(Form("h2_%i",p),Form("Photon Pt:%s",WeightNames[p].c_str()) , n_bins, PT_Min, PT_Max); // will be weighted, is this redundant?
			h3[p] = new TH1D(Form("h3_%i",p),Form("Pion PT, weighted:%s",WeightNames[p].c_str()) , n_bins, PT_Min, PT_Max);
			h12[p] = new TH1F(Form("h12_%i",p),Form("Smeared Pion PT, weighted:%s",WeightNames[p].c_str()) , n_bins, PT_Min, PT_Max);
			h20[p] = new TH1F(Form("h20_%i",p),Form("Smeared Photon pT, weighted:%s",WeightNames[p].c_str()) , n_bins, PT_Min, PT_Max);
			h21[p] = new TH1F(Form("h21_%i",p),Form("Photon pT, weighted:%s",WeightNames[p].c_str()) , n_bins, PT_Min, PT_Max);

			h18[p] = new TH2F(Form("h18_%i",p),Form("Smeared Pion Pt vs Smeared Inv Mass, weighted:%s",WeightNames[p].c_str()) , n_bins, 0, PT_Max, MassNBins, smeared_lower_bin_limit,  smeared_upper_bin_limit);
			h18_1d[p] = new TH1F(Form("h18_1d_%i",p),Form("Smeared Pion Pt vs Smeared Inv Mass, weighted:%s",WeightNames[p].c_str()),MassNBins, smeared_lower_bin_limit,  smeared_upper_bin_limit);

			h27[p] = new TH2F(Form("h27_%i",p), Form("Smeared Pion Pt vs Smeared Inv Mass, weighted. cluster:%s",WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, smeared_lower_bin_limit,  smeared_upper_bin_limit);
			h28[p] = new TH2F(Form("h28_%i",p), Form("Smeared Pion Pt vs Smeared Inv Mass, weighted. cluster and Asymm Cut:%s",WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, smeared_lower_bin_limit,  smeared_upper_bin_limit);
			h29[p] = new TH2F(Form("h29_%i",p), Form("Smeared Pion Pt vs Smeared Inv Mass, weighted. Asymm Cut:%s",WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, smeared_lower_bin_limit,  smeared_upper_bin_limit);

            h28_v2[p] = new TH2F(Form("h28_v2_%i",p), Form("Smeared Pion Pt vs Smeared Inv Mass, weighted. Cluster+Asymm+pos res:%s",WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, smeared_lower_bin_limit,  smeared_upper_bin_limit);
			h29_v2[p] = new TH2F(Form("h29_v2_%i",p), Form("Smeared Pion Pt vs Smeared Inv Mass, weighted. Asymm+Pos res:%s",WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, smeared_lower_bin_limit,  smeared_upper_bin_limit);

			h30[p] = new TH2F(Form("h30_%i",p), Form("Smeared Pion Pt vs Smeared Inv Mass, weighted. Blair's cuts, no pos.res no occupancy:%s",WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, 0, smeared_upper_bin_limit);
			h30_1d[p] = new TH1F(Form("h30_1d_%i",p), Form("Smeared Inv Mass, weighted. Blair's cuts, no pos.res no occupancy:%s",WeightNames[p].c_str()), MassNBins, 0, smeared_upper_bin_limit);

			h31[p] = new TH2F(Form("h31_%i",p), Form("Smeared Pion Pt vs Smeared Inv Mass, weighted. Blair's cuts+pos res:%s",WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, 0, smeared_upper_bin_limit);
			h31_1d[p] = new TH1F(Form("h31_1d_%i",p), Form("Smeared Inv Mass, weighted. Blair's cuts+pos res:%s",WeightNames[p].c_str()), MassNBins, 0, smeared_upper_bin_limit);
			h31_pionspectrum[p] = new TH1F(Form("h31_ps_%i",p),Form("Cuts+pos res, Smeared Pion PT, weighted:%s",WeightNames[p].c_str()) , MassNBins, PT_Min, PT_Max);

			h34[p] = new TH2F(Form("h34_%i",p), Form("Smeared Pion Pt vs Smeared Inv Mass, weighted. Position Smearing:%s",WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, smeared_lower_bin_limit,  smeared_upper_bin_limit);
			h34_1d[p] = new TH1F(Form("h34_1d_%i",p), Form("Smeared Pion Pt vs Smeared Inv Mass, weighted. Position Smearing:%s",WeightNames[p].c_str()), MassNBins, smeared_lower_bin_limit,  smeared_upper_bin_limit);

			h35[p] = new TH2F(Form("h35_%i",p), Form("Smeared Pion Pt vs Smeared Inv Mass, weighted. Blair's cuts+cluster:%s",WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, 0, smeared_upper_bin_limit);
			h35_1d[p] = new TH1F(Form("h35_1d_%i",p), Form("Smeared Inv Mass, weighted. Blair's cuts+cluster:%s",WeightNames[p].c_str()), MassNBins, 0, smeared_upper_bin_limit);



			h100[p] = new TH2F(Form("h100_%i",p), Form("Smeared Pion Pt vs Smeared Inv Mass, weighted. All Cuts+effects:%s",WeightNames[p].c_str()), n_bins, 0, PT_Max, MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
			h100_1d[p] = new TH1F(Form("h100_1d_%i",p), Form("Smeared Pion Pt vs Smeared Inv Mass, weighted. All Cuts+effects:%s",WeightNames[p].c_str()), MassNBins, smeared_lower_bin_limit, smeared_upper_bin_limit);
		}

		// things to add probability to add an exponentially scaled (small) energy
		//--------------------set up random number generation
		std::random_device rd;// generate a random number to seed random generation
		std::random_device rdgamma;// generate a random number to seed random generation of daughter gamma for smearing
		std::random_device rdgammacluster;// random number to test clustering
		std::random_device rdgammapositsmr;// random number to test position smearing
		// std::mt19937 gen(rd());          // mersenne_twister_engine seeded with rd
		std::mt19937_64 gen(rd());            // mersenne_twister_engine 64 bit seeded with rd
		std::mt19937_64 gen_gamma(rdgamma()); // mersenne_twister_engine 64 bit seeded with rdgamma for gamma smearing
		std::mt19937_64 gen_gammacluster(rdgammacluster());//generate random number to test clustering 
		std::mt19937_64 gen_gammapositsmear(rdgammapositsmr());//generate random number to test clustering 
		// std::ranlux48 gen(rd()); // ranlux48 seeded with rd
		//--------------------generate random: momentum
		// std::uniform_real_distribution<> pdis(0.0,PT_Max);
		std::normal_distribution<double> gammadis(0.0, 1.0);  // generate normal distribution for gamma smearing, mean zero, variance 1
		std::uniform_real_distribution<> gammacluster(0, 1.0); // random probability to have photons clustered together. if above certain value, add a random photon from the smeared+weighted photon spectrum. cutoff will be tuned.
		std::normal_distribution<double> gamma_positsmear(0.0, 1.0);//normally distributed position smearing
		std::uniform_real_distribution<> pdis(PT_ratio, 1.0); // alternative scheme with min PT to avoid power law complications.
		//--------------------generate random: angle
		// double tpi=2*std::numbers::pi;
		std::uniform_real_distribution<> adis(0.0, 2 * M_PI);
		int id, size, no;
		Vec4 gamma_lorentz[3];
		Vec4 gamma_smeared[3];
		Vec4 gamma_cluster[3];
		Vec4 gamma_cluster_asymm[3];
		Vec4 gamma_Blair_Cuts[3];
		Vec4 gamma_All_Cuts[3];
		Vec4 gamma_position_smear[3];
		Vec4 gamma_Blair_position[3];

		double m, E, px, py, pz, eta_px, eta_py, eta_E, scale_factor1, scale_factor2, smear_factor1, smear_factor2, mT_Scaling;
        double eta_rest = 0.0;
        double eta_pz = 0.0;
        double Eta_M = 0.54786;
        double Pion_M = 0.1349768;
		double inv_mass, inv_mass_smeared;

		std::vector<int> WeightScale(WeightNames.size());
		std::vector<double> weight_function(WeightNames.size()), inv_yield(WeightNames.size());

		// std::cout << Pi0_M <<std::endl;
		//tree->Branch("id", &id, "id/I");
		// tree->Branch("size",&size, "size/I");
		// tree->Branch("no",&no, "no/I");
		//tree->Branch("m", &m, "m/D");
		//tree->Branch("E", &E, "E/D"); // can reconstruct in root easilly
		// tree->Branch("Pt",&Pt, "Pt/D"); // can reconstruct in root easilly
		//tree->Branch("px", &px, "px/D");
		//tree->Branch("py", &py, "py/D");
		//tree->Branch("pz", &pz, "pz/D");

		// Set up generation.  event
		Pythia pythia;								// Declare Pythia object
		pythia.readString("PromptPhoton:all = on"); // Switch on process.
		//pythia.readString("ParticleDecays:allowPhotonRadiation = on");
		// pythia.readString("SoftQCD:all = on"); // Switch on process.
		// pythia.readString("Beams:eCM = 14.e3"); // 14 TeV CM energy.
		//pythia.readString("221:all = eta -> gamma gamma");
		pythia.readString("221:oneChannel = 1 1.0 0 22 22");// force eta->gammagamma with 100% branching ratio?
		pythia.init(); 
		pythia.event.clear();

		int etaCount = 0;
		for (int i = 0; i < NPions; i++){//create primaries
            //if ((i+1 % 100000) == 0) std::cout << i+1 << std::endl;
			// std::cout << "I reached here" <<" "<< i <<std::endl;// debug line
			double azimuthal_ang = adis(gen); // generate a random angle from 0 to 2pi
			double Pt = PT_Max * pdis(gen);
			//----------------------different possible weights
			//std::cout << "EXP Weight" <<std::endl;
			weight_function[0]=exp(-Pt/0.3);//originally dividing by 0.2
			WeightScale[0]=1e+14;
            mT_Scaling=R_EP_Inf*pow((mT_a+sqrt(pow(Eta_M,2)+pow(Pt,2)))/(mT_a+sqrt(pow(Pion_M,2)+pow(Pt,2))),mT_n);

			inv_yield[0] = mT_Scaling*WeightScale[0]* Pt * weight_function[0];

			//std::cout << "Power Weight" <<std::endl;
			weight_function[1]=pow(Pt,-8.14);
			WeightScale[1]=1e+5;
			inv_yield[1] = mT_Scaling*WeightScale[0]* Pt * weight_function[1];

			//std::cout << "WSHP Weight" <<std::endl;
			weight_function[2]=((1/(1+exp((Pt-t)/w)))*A/pow(1+Pt/p0,m_param)+(1-(1/(1+exp((Pt-t)/w))))*B/(pow(Pt,n)));
			WeightScale[2]=1e+14;
			inv_yield[2] = mT_Scaling*WeightScale[2]* Pt * weight_function[2];

			//std::cout << "Hagedorn Weight" <<std::endl;
			weight_function[3]=A/pow(1+Pt/p0,m_param);
			WeightScale[3]=1e+14;
			inv_yield[3] = mT_Scaling*WeightScale[3]* Pt * weight_function[3];

			//double inv_yield = WeightScale* Pt * weight_function;
			

			//std::cout << inv_yield <<std::endl;
			//h3->Fill(Pt, inv_yield); // fill pi0 pt, weighted
			//h4->Fill(Pt);						// fill pi0 pt, unweighted
			eta_px = Pt * cos(azimuthal_ang);
			eta_py = Pt * sin(azimuthal_ang);
			eta_E = sqrt(Eta_M * Eta_M + Pt * Pt + eta_pz * eta_pz);
			//pythia.event.append(221, 23, 0, 0, eta_px, eta_py, eta_pz, eta_E, Eta_M); // 111 is pi0 add a partice, 221 is eta
			if (pythia.event.append(221, 23, 0, 0, eta_px, eta_py, eta_pz, eta_E, Eta_M)) {
            etaCount++;
        	} else {
            	std::cout << "Failed to append eta meson at iteration " << i << std::endl;
        	}


			// append(pid, use of particle flag, mother?, daughter?, px, py, pz, E, m)
			//id = pythia.event[i].id();
			//m = pythia.event[i].m();
			//px = pythia.event[i].px();
			//py = pythia.event[i].py();
			//pz = pythia.event[i].pz();
			//E = pythia.event[i].e();
			// std::cout << i <<" "<<id<< " " << m << " " << " E " << " " << E <<" " << " P " << " " << px<<" "<<py<<" "<<pz<<std::endl;
			//tree->Fill();
		}
		std::cout << "Generated eta mesons: " << etaCount << std::endl;
		pythia.moreDecays();
		//std::cout << "I reached here" << std::endl; // debug line

		for (int i = 0; i < pythia.event.size(); i++){ // loop over all events(primary)
			if ((i+1 % 100000) == 0) std::cout << i+1 << std::endl;
            if (pythia.event[i].id() == 221){ // if the ith event is an eta

				int Gamma_daughters[2] = {pythia.event[i].daughter1(), pythia.event[i].daughter2()}; // make array of daughter particles(di gamma) event ids
				double Pt = pythia.event[i].pT();
				mT_Scaling=R_EP_Inf*pow((mT_a+sqrt(pow(Eta_M,2)+pow(Pt,2)))/(mT_a+sqrt(pow(Pion_M,2)+pow(Pt,2))),mT_n);
				//std::cout << "EXP Weight" <<std::endl;
				weight_function[0]=exp(-Pt/0.3);//originally dividing by 0.2
				WeightScale[0]=1e+14;
				inv_yield[0] =mT_Scaling*WeightScale[0]* Pt * weight_function[0];

				//std::cout << "Power Weight" <<std::endl;
				weight_function[1]=pow(Pt,-8.14);
				WeightScale[1]=1e+5;
				inv_yield[1] = mT_Scaling*WeightScale[0]* Pt * weight_function[1];

				//std::cout << "WSHP Weight" <<std::endl;
				weight_function[2]=((1/(1+exp((Pt-t)/w)))*A/pow(1+Pt/p0,m_param)+(1-(1/(1+exp((Pt-t)/w))))*B/(pow(Pt,n)));
				WeightScale[2]=1e+14;
				inv_yield[2] = mT_Scaling*WeightScale[2]* Pt * weight_function[2];

				//std::cout << "Hagedorn Weight" <<std::endl;
				weight_function[3]=A/pow(1+Pt/p0,m_param);
				WeightScale[3]=1e+14;
				inv_yield[3] = mT_Scaling*WeightScale[3]* Pt * weight_function[3];

				if (pythia.event[Gamma_daughters[0]].id() == 22 && pythia.event[Gamma_daughters[1]].id() == 22){//filling histograms
					gamma_lorentz[0] = pythia.event[Gamma_daughters[0]].p();
					gamma_lorentz[1] = pythia.event[Gamma_daughters[1]].p();
					gamma_lorentz[2] = gamma_lorentz[0] + gamma_lorentz[1];
					inv_mass = gamma_lorentz[2].mCalc();

					scale_factor1 = sqrt(pow(smear_factor_b, 2) / gamma_lorentz[0].e() + pow(smear_factor_c, 2) + pow(smear_factor_d, 2));
					scale_factor2 = sqrt(pow(smear_factor_b, 2) / gamma_lorentz[1].e() + pow(smear_factor_c, 2) + pow(smear_factor_d, 2));

					smear_factor1 = scale_factor1 * gammadis(gen_gamma) + 1;
					smear_factor2 = scale_factor2 * gammadis(gen_gamma) + 1;

					gamma_smeared[0] = smear_factor1 * pythia.event[Gamma_daughters[0]].p();
					gamma_smeared[1] = smear_factor2 * pythia.event[Gamma_daughters[1]].p();

                    //clustering loop
                    for(int photclust=0;photclust<2;photclust++ ){
                        if (gammacluster(gen_gammacluster)>coprob && clusteroverlay==1){
						gamma_cluster[photclust] = clusterPhoton(gamma_smeared[photclust], 2, myFunc->GetRandom());
						gamma_cluster_asymm[photclust]=gamma_cluster[photclust];
                        }
                        else{
                            gamma_cluster[photclust]=gamma_smeared[photclust];
                            gamma_cluster_asymm[photclust]=gamma_cluster[photclust];
                        }
                    }

					gamma_All_Cuts[0]=PositionResSmear(gamma_cluster_asymm[0], posit_smearingFactor*gamma_positsmear(gen_gammapositsmear), posit_smearingFactor*gamma_positsmear(gen_gammapositsmear), posit_smearingFactor*gamma_positsmear(gen_gammapositsmear));
					gamma_All_Cuts[1]=PositionResSmear(gamma_cluster_asymm[1], posit_smearingFactor*gamma_positsmear(gen_gammapositsmear), posit_smearingFactor*gamma_positsmear(gen_gammapositsmear), posit_smearingFactor*gamma_positsmear(gen_gammapositsmear));
					gamma_All_Cuts[2] = gamma_All_Cuts[0] + gamma_All_Cuts[1];

					gamma_Blair_position[0]=PositionResSmear(gamma_smeared[0], posit_smearingFactor*gamma_positsmear(gen_gammapositsmear), posit_smearingFactor*gamma_positsmear(gen_gammapositsmear), posit_smearingFactor*gamma_positsmear(gen_gammapositsmear));
					gamma_Blair_position[1]=PositionResSmear(gamma_smeared[1], posit_smearingFactor*gamma_positsmear(gen_gammapositsmear), posit_smearingFactor*gamma_positsmear(gen_gammapositsmear), posit_smearingFactor*gamma_positsmear(gen_gammapositsmear));
					gamma_Blair_position[2]= gamma_Blair_position[0] + gamma_Blair_position[1];

					gamma_smeared[2] = gamma_smeared[0] + gamma_smeared[1];
					gamma_cluster[2] = gamma_cluster[0] + gamma_cluster[1];
					gamma_cluster_asymm[2] = gamma_cluster_asymm[0] + gamma_cluster_asymm[1];

					gamma_Blair_Cuts[0]=gamma_cluster_asymm[0];
					gamma_Blair_Cuts[1]=gamma_cluster_asymm[1];
					gamma_Blair_Cuts[2]=gamma_cluster_asymm[2];

					gamma_position_smear[0]=PositionResSmear(gamma_smeared[0], posit_smearingFactor*gamma_positsmear(gen_gammapositsmear), posit_smearingFactor*gamma_positsmear(gen_gammapositsmear), posit_smearingFactor*gamma_positsmear(gen_gammapositsmear));
					gamma_position_smear[1]=PositionResSmear(gamma_smeared[1], posit_smearingFactor*gamma_positsmear(gen_gammapositsmear), posit_smearingFactor*gamma_positsmear(gen_gammapositsmear), posit_smearingFactor*gamma_positsmear(gen_gammapositsmear));

					gamma_position_smear[2]=gamma_position_smear[0]+gamma_position_smear[1];

					inv_mass_smeared = gamma_smeared[2].mCalc();
					h4->Fill(Pt);			// fill pion pt, unweighted
					h6->Fill(inv_mass);
					h8->Fill(inv_mass_smeared);
					h9->Fill(gamma_smeared[2].pT(), inv_mass_smeared);//change  to smeared pion pT. was previously unsmeared: pythia.event[i].pT()
					h10->Fill(gamma_smeared[2].pT());	
					h17->Fill(gamma_lorentz[0].pT());//unsmeared energy spectrum
					h17->Fill(gamma_lorentz[1].pT());
					h16->Fill(gamma_smeared[0].pT());//smeared photon energy spectrum
					h16->Fill(gamma_smeared[1].pT());

					for (int p=0; p < WeightNames.size(); p++){

						h20[p]->Fill(gamma_smeared[0].pT(), inv_yield[p]);
						h20[p]->Fill(gamma_smeared[1].pT(), inv_yield[p]);
						h21[p]->Fill(gamma_lorentz[0].pT(), inv_yield[p]);
						h21[p]->Fill(gamma_lorentz[1].pT(), inv_yield[p]);						
						h3[p]->Fill(Pt, inv_yield[p]); // fill pion pt, weighted
						h12[p]->Fill(gamma_smeared[2].pT(), inv_yield[p]);
						h18[p]->Fill(gamma_smeared[2].pT(), inv_mass_smeared, inv_yield[p]);//
						h18_1d[p]->Fill(inv_mass_smeared, inv_yield[p]);//
						h27[p]->Fill(gamma_cluster[2].pT(), gamma_cluster[2].mCalc(), inv_yield[p]);
						h34[p]->Fill(gamma_position_smear[2].pT(), gamma_position_smear[2].mCalc(), inv_yield[p]);//position smearing
						h34_1d[p]->Fill(gamma_position_smear[2].mCalc(), inv_yield[p]);//position smearing

                        if(AsymmCutcheck(gamma_smeared[0], gamma_smeared[1], asymval,asymcut)==true){// asymm
                            h29[p]->Fill(gamma_smeared[2].pT(), gamma_smeared[2].mCalc(), inv_yield[p]);
                        }
                        if(AsymmCutcheck(gamma_cluster_asymm[0], gamma_cluster_asymm[1], asymval,asymcut)==true){//cluster+asymm
                            h28[p]->Fill(gamma_cluster_asymm[2].pT(), gamma_cluster_asymm[2].mCalc(), inv_yield[p]);
                        }
                        if(AsymmCutcheck(gamma_position_smear[0], gamma_position_smear[1], asymval,asymcut)==true){//asymm+pos res
                            h29_v2[p]->Fill(gamma_position_smear[2].pT(), gamma_position_smear[2].mCalc(), inv_yield[p]);
                        }
                        if(AsymmCutcheck(gamma_All_Cuts[0], gamma_All_Cuts[1], asymval,asymcut)==true){//cluster+asymm+pos res
                            h28_v2[p]->Fill(gamma_All_Cuts[2].pT(), gamma_All_Cuts[2].mCalc(), inv_yield[p]);
                        }

                        if(DeltaRcut(gamma_smeared[0], gamma_smeared[1], DeltaRcut_MAX)==false && AsymmCutcheck(gamma_smeared[0], gamma_smeared[1], asymval,asymcut)==true && pTCut(gamma_smeared[0], pt1cut)==true && pTCut(gamma_smeared[1], pt2cut)==true && nclus_ptCut<gamma_smeared[0].pT() && gamma_smeared[0].pT()<ptMaxCut && nclus_ptCut<gamma_smeared[1].pT() && gamma_smeared[1].pT()<ptMaxCut && gamma_smeared[2].pT()>comb_ptcut*(pt1cut+pt2cut)){
                            h30[p]->Fill(gamma_smeared[2].pT(), gamma_smeared[2].mCalc(), inv_yield[p]);// data+mc cuts, no pos res or occupancy
                            h30_1d[p]->Fill(gamma_smeared[2].mCalc(), inv_yield[p]);
                        }
                        
                        if(DeltaRcut(gamma_Blair_position[0], gamma_Blair_position[1], DeltaRcut_MAX)==false && AsymmCutcheck(gamma_Blair_position[0], gamma_Blair_position[1], asymval,asymcut)==true  && pTCut(gamma_Blair_position[0], pt1cut)==true && pTCut(gamma_Blair_position[1], pt2cut)==true  && nclus_ptCut<gamma_Blair_position[0].pT() && gamma_Blair_position[0].pT()<ptMaxCut && nclus_ptCut<gamma_Blair_position[1].pT() && gamma_Blair_position[1].pT()<ptMaxCut && gamma_Blair_position[2].pT()>comb_ptcut*(pt1cut+pt2cut)){
                            h31[p]->Fill(gamma_Blair_position[2].pT(), gamma_Blair_position[2].mCalc(), inv_yield[p]);// data+mc cuts,  pos res on ,no occupancy
                            h31_1d[p]->Fill(gamma_Blair_position[2].mCalc(), inv_yield[p]);
                            h31_pionspectrum[p]->Fill(gamma_Blair_position[2].pT(), inv_yield[p]);
                            hpionpt[p]->Fill(gamma_lorentz[2].pT());
                            hInvMass_Cutson->Fill(gamma_Blair_position[2].mCalc());
                        }

                        if(DeltaRcut(gamma_Blair_Cuts[0], gamma_Blair_Cuts[1], DeltaRcut_MAX)==false && AsymmCutcheck(gamma_Blair_Cuts[0], gamma_Blair_Cuts[1], asymval,asymcut)==true && pTCut(gamma_Blair_Cuts[0], pt1cut)==true && pTCut(gamma_Blair_Cuts[1], pt2cut)==true  && nclus_ptCut<gamma_Blair_Cuts[0].pT() && gamma_Blair_Cuts[0].pT()<ptMaxCut && nclus_ptCut<gamma_Blair_Cuts[1].pT()&& gamma_Blair_Cuts[1].pT()<ptMaxCut && gamma_Blair_Cuts[2].pT()>comb_ptcut*(pt1cut+pt2cut)){
                            h35[p]->Fill(gamma_Blair_Cuts[2].pT(), gamma_Blair_Cuts[2].mCalc(), inv_yield[p]);// data+mc cuts,  pos res off , occupancy on
                            h35_1d[p]->Fill(gamma_Blair_Cuts[2].mCalc(), inv_yield[p]);
                        }

                        if(DeltaRcut(gamma_All_Cuts[0], gamma_All_Cuts[1], DeltaRcut_MAX)==false && AsymmCutcheck(gamma_All_Cuts[0], gamma_All_Cuts[1], asymval,asymcut)==true && pTCut(gamma_All_Cuts[0], pt1cut)==true && pTCut(gamma_All_Cuts[1], pt2cut)==true  && nclus_ptCut<gamma_All_Cuts[0].pT() && gamma_All_Cuts[0].pT()<ptMaxCut && nclus_ptCut<gamma_All_Cuts[1].pT() && gamma_All_Cuts[1].pT()<ptMaxCut && gamma_All_Cuts[2].pT()>comb_ptcut*(pt1cut+pt2cut)){
                            h100[p]->Fill(gamma_All_Cuts[2].pT(), gamma_All_Cuts[2].mCalc(), inv_yield[p]);// blair+pos+occupancy
                            h100_1d[p]->Fill(gamma_All_Cuts[2].mCalc(), inv_yield[p]);
                        }

					}

					try
					{
						mass_pt_map[floor(Pt) + 1].push_back(inv_mass_smeared);
					}
					catch (...)
					{ // ellipses in catch argument means handles ALL exceptions.
						mass_pt_map.insert({floor(Pt) + 1, std::vector<double>()});
					}
					//*/
					for (int j : Gamma_daughters)
					{ // loop over daughter particles
						//id = pythia.event[j].id();
						//m = pythia.event[j].m();
						//px = pythia.event[j].px();
						//py = pythia.event[j].py();
						//pz = pythia.event[j].pz();
						//E = pythia.event[j].e();
						// std::cout << i <<" "<<id<< " " << m << " " << " E " << " " << E <<" " << " P " << " " << px<<" "<<py<<" "<<pz<<std::endl;
						//tree->Fill();
						h5->Fill(pythia.event[j].pT()); // fill gamma pt, unweighted
					}
				}
				else
				{
					// printf("non gamma decay\n");
					// printf("daughter1 ID=%g\n",pythia.event[Gamma_daughters[0]].id());
					// printf("daughter2 ID=%g\n",pythia.event[Gamma_daughters[1]].id());
				}
			}			
		}
			
		double_t realtime = timer.RealTime();
		double_t cputime = timer.CpuTime();

		printf("real time =%f\n", realtime);
		printf("CPU time=%f\n", cputime);

		output->Write();
		output->Close();

		//clean up
		//delete tree; 

        delete output;
	}
    
    return 0;
}

TF1* ChooseSpectrumFunction(int weightmethod, int PT_Min, int PT_Max){
	TF1 *myFunc = nullptr;
	if(weightmethod==0){
		myFunc = new TF1("myFunc", [](double *x, double *par){
			return x[0]*par[1]*pow((par[2]+sqrt(pow(par[3],2)+pow(x[0],2)))/(par[2]+sqrt(pow(par[4],2)+pow(x[0],2))),-par[5])*par[0] * TMath::Exp(-x[0]/0.3);
		}, PT_Min, PT_Max, 6);
		//Double_t initialParameters[1] = {1.0};
    	myFunc->SetParameters(1.0, 0.5, 1.2, 0.54786, 0.1349768,10);

	}
	else if(weightmethod==1){
		myFunc = new TF1("myFunc", [](double *x, double *par){
			return x[0]*par[1]*pow((par[2]+sqrt(pow(par[3],2)+pow(x[0],2)))/(par[2]+sqrt(pow(par[4],2)+pow(x[0],2))),-par[5])*par[0] * pow(x[0],-8.14);
			}, PT_Min, PT_Max, 6);
		//Double_t initialParameters[1] = {1.0};
    	myFunc->SetParameters(1.0, 0.5, 1.2, 0.54786, 0.1349768,10);
	}
	else if(weightmethod==2){
		myFunc = new TF1("myFunc", [](double *x, double *par){
			return x[0]*
			par[7]*pow((par[8]+sqrt(pow(par[9],2)+pow(x[0],2)))/(par[8]+sqrt(pow(par[10],2)+pow(x[0],2))),-par[11])
			*((1/(1+exp((x[0]-par[0])/par[1])))*par[2]/pow(1+x[0]/par[3],par[4])+(1-(1/(1+exp((x[0]-par[0])/par[1]))))*par[5]/(pow(x[0],par[6])));
		}, PT_Min, PT_Max, 12);
		Double_t params[] = {4.5, 0.114, 229.6, 1.466, 10.654, 14.43, 8.1028, 0.5, 1.2, 0.54786, 0.1349768, 10.0};
    	myFunc->SetParameters(params);
		//myFunc->SetParameters(4.5,  0.114, 229.6, 1.466, 10.654, 14.43, 8.1028, 0.5, 1.2, 0.54786, 0.1349768, 10.0);
		//t,w,A,p0,m_param,B,n
	}
	else if(weightmethod==3){
		//std::cout << "Hagedorn Weight" <<std::endl;
		myFunc = new TF1("myFunc", [](double *x, double *par){
			return x[0]*par[3]*pow((par[4]+sqrt(pow(par[5],2)+pow(x[0],2)))/(par[4]+sqrt(pow(par[6],2)+pow(x[0],2))),-par[7])*(par[0]/pow(1+x[0]/par[1],par[2]) );
			}, PT_Min, PT_Max, 8);
		myFunc->SetParameters( 229.6, 1.466, 10.654, 0.5, 1.2, 0.54786, 0.1349768,10);
		//weight_function=A/pow(1+Pt/p0,m_param);
		//WeightScale=1e+14;
	}
	else{
		std::cout << "Error:No Weight function found" <<std::endl;
	}
	return myFunc;
}

Pythia8::Vec4 clusterPhoton(Pythia8::Vec4& originalPhoton, int method, double randomE) {
	// I see two methods to do this. 
	// 1) draw a random energy, pz and phi then construct px and py
	// this is not really realistic because photons with certain momenta could never reach certain towers
	// 2) draw a random energy, copy the momentum from the original 4 vector and scale them down appropriately
	// 2.1) add some randomization
	// also non realistic. I am not sure what level of realism is needed for this though
	// 3) do a proper analysis for each energy. I would have to think it through.
    // Create a Pythia8 random number generator
	//Pythia8::Vec4 photongen;
	Pythia8::Vec4 newPhoton;
	Pythia8::Rndm rndm;
	if (method==1){
		
		//double theta = rndm.flat(0.0, M_PI);
		//double phi = rndm.flat(0.0, 2.0 * M_PI);
		//double theta = acos(photonPz/photonEnergy);// theta should be related to the energy
		

		// Calculate the momentum components based on the angles and energy
		//double photonPx = randomE * sin(theta) * cos(phi);
		//double photonPy = randomE * sin(theta) * sin(phi);
		//double photonPz = randomE * cos(theta);

		// Create a 4-vector for the new photon
		//Pythia8::Vec4 photon(photonPx, photonPy, photonPz, randomE);
		//newPhoton = photon;
	}
	else if(method==2){
		newPhoton.e(randomE);
		
		newPhoton.px(originalPhoton.px() * (randomE / originalPhoton.e()));
		newPhoton.py(originalPhoton.py() * (randomE / originalPhoton.e()));	
		newPhoton.pz(originalPhoton.pz() * (randomE / originalPhoton.e()));		
	}
	else if(method==3){
		
	}
    



    // Return the sum of the original 4-vector and the new photon 4-vector
    //return originalVector + photon;
	return newPhoton+originalPhoton;
}

Pythia8::Vec4 PositionResSmear(Pythia8::Vec4 photon, double smearingFactorx,double smearingFactory,double smearingFactorz) {//, Pythia8::Rndm* rndm

	//calculate the pT to keep the length of the four momenta the same. then smear the px and py to match.
	//*
	//double xsmear = smearingFactorx + 1;
	//there is a position vector on the detector
	//new vector which is x/r, y/r, and z/r
	// multiply by e to get px py pz
	// I can do the reverse now px=e*x/r.
	// reverse is e*x/px=r. px/e=x/|r|. I know the inner radius of EMCALL. I can use that for the magnitude of r? 
	// increase length of vector until it intersects with the EMCAL barrel
	// for p.rapidity zero r=EMCAL radius. off
	// now I have a position vector. I can smear the position vector. and then go back to the momentum vector. 
	//EMCAL inner radius =900mm

	//phnx pos res was like 5.7 mm/sqrt(E). try multiplying by 2.5 as a first guess as to the pos res scale.
	// scale not by towersize in p rapid? but by abs tower size?
	//tower size is 2.5 cm and it is half the size of phenix
	// module size is approx moiliere rad for sphenix and phenix. sphenix emcal is at a smaller rad than phenix emcal. so sphenix needs to be made of a material that has a smaller m radius? it is. some tungsten powder. so we need to look at the tower size to find the pos res. The block size is 54 mm. so 5.4 cm. so the pos res should be roughly a factor 2 less than phenix.
	// spehnix should be ~2.8 mm

	double energy = photon.e();
	//construct the position vector

	double posx =900*photon.px()/photon.e();// 900 mm for EMCAL inner radius
	double posy =900*photon.py()/photon.e();
	double posz =900*photon.pz()/photon.e();

	// smear the position vectors. abs smearing coming from the pos res of phnx(for now)

	posx += smearingFactorx/sqrt(energy);
	posy += smearingFactory/sqrt(energy);
	posz += smearingFactorz/sqrt(energy);

	// return to momentum vectors.

	double px_smear=posx*photon.e()/900;
	double py_smear=posy*photon.e()/900;
	double pz_smear=posz*photon.e()/900;

    // Keep the energy constant
    
	double E_New = sqrt(px_smear*px_smear +py_smear*py_smear +pz_smear*pz_smear);
	double energyscale=energy/E_New;

	double xnew=energyscale*px_smear;
	double ynew=energyscale*py_smear;
	double znew=energyscale*pz_smear;


    Pythia8::Vec4 smearedPhoton(xnew, ynew, znew, energy);// sqrt(xnew*xnew +ynew*ynew +znew*znew)
    return smearedPhoton;
}
   	



// Function to check if Delta R between two four-vectors is greater than a specified cut
bool DeltaRcut(Pythia8::Vec4& Photon1, Pythia8::Vec4& Photon2, float DeltaRcutMax) {
    double dEta = Photon1.eta() - Photon2.eta();
    double dPhi = acos(cos(Photon1.phi() - Photon2.phi()));  // Correct way to handle the periodicity
    double deltaR = sqrt(dEta * dEta + dPhi * dPhi);
    return deltaR > DeltaRcutMax;
}

// Function to check if the transverse momentum (pT) of a particle is greater than a specified cut
bool pTCut(const Pythia8::Vec4& particle, float ptCut) {
    // Calculate the transverse momentum (pT)
    double pT = particle.pT();

    // Check if pT is greater than the cut value
    return pT > ptCut;
}

// Function to check if the photons pass the asymmetry cut.
bool AsymmCutcheck(Pythia8::Vec4& Photon1, Pythia8::Vec4& Photon2, float AsymmCutoff,bool asymcutbool) {
    if(asymcutbool==false){
        return true;//asymmcut is off, so proceeed with fill logic
    }
    return abs(Photon1.e()-Photon2.e())/(Photon1.e()+Photon2.e())<AsymmCutoff;
}
