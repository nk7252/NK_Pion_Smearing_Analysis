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
Pythia8::Vec4 PositionResSmear(Pythia8::Vec4 photon, double smearingFactor);
bool DeltaRcut(Pythia8::Vec4& Photon1, Pythia8::Vec4& Photon2, float DeltaRcutMax);
bool pTCut(const Pythia8::Vec4& particle, float ptCut);

int main(){ 
	TStopwatch timer;
	timer.Start();
	int NPions = 1 * 1000000;//
	// int n_bins=1+ceil(log2(Nevents));
	int PT_Max = 64; // 65
	int PT_Min = 0;	 // cross check to elimate power law problems
	double PT_ratio = PT_Min / PT_Max;
	// int n_bins=round((1/4)*PT_Max);
	// int n_bins=round(1+3.222*log(NPions));
	int binres=2;//number of divisions per GeV
	int n_bins = binres*PT_Max;//multiple by bin res.
	std::map<double, std::vector<double>> mass_pt_map; // we want to have keys of a pT range?


	//-----------------------------------set weighting method
	//int weightmethod = 3;//0=exp,1=power,2=wshp, 3=hagedorn(not implemented)
	std::vector<std::string> WeightNames = {"EXP", "POWER", "WSHP","HAGEDORN"};
	//-----------------------------------
	//cuts
	int asymcut=1;//apply asymm cut.
	float asymval= 0.6;
	int clusteroverlay = 1;//overlayed cluster check
	float coprob=0.8;//random numbers(0-1) greater than this value will have some smearing added.

	//Blair specific cuts
	float DeltaRcut_MAX = 1.1; 
	float pt1cut =1.3;   
	float pt2cut =0.7; 
	//nclus is not something I can do here is it          

	//--------------------Alternative paramaterization, woods saxon+hagedorn+power law
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
	float smeared_upper_bin_limit = 0.25;
	float smear_factor_a = 0;
	float smear_factor_d = 0.02;  // 0.02;// test trying to include the beam momentum resolution.. will set to zero for now
	float smear_factor_c = 0.028; // first parameter in test beam parameterization?
	float posit_smearingFactor = 0.1; // Example smearing factor for position

	//std::cout << "Processing: " << WeightNames[weightmethod] << std::endl;
	//----------------------pion spectrum function for clusteroverlay
	// reserve a TF1 for the chosen function just in case
	TF1 *myFunc;
	myFunc=ChooseSpectrumFunction(2, PT_Min, PT_Max);

	for (int smear_factor_itt = 9; smear_factor_itt < 9 + 1; smear_factor_itt++)
	{// originally int smear_factor_itt = 0; smear_factor_itt < 24 + 1; smear_factor_itt++
	// only want .155
		float smear_factor_basevalue = 0.065; // I used 1.6% + 12.7%/sqrt(E) fig 22, but that is from a special beam cross section config. trying with fig 24 data i.e 2.8% + 15.5%
		//--------------------preliminaries to read from root
		float smear_factor_b = smear_factor_basevalue + 0.01 * smear_factor_itt;
		//float smear_factor_b = 0.03;
		//////////////////////New//0.155 loop from twice test beam data paramaterization to half? this is 15.5% from https://arxiv.org/pdf/1704.01461.pdf fig 24b so going from 6.5% to 30.5%// need 24 steps for 1% diff each
		//////////////////////OLD//0.127 loop from twice test beam data paramaterization to half? this is 12.7% from https://arxiv.org/pdf/1704.01461.pdf fig 22b so going from 6.35% to 25.4%

		TFile *output = new TFile(Form("pioncode/rootfiles/Pi0FastMC_%f.root", smear_factor_b), "recreate");//

		//TFile *output = new TFile(Form("pioncode/rootfiles/Pi0FastMC_%f_%s_ac%i_co%i.root", smear_factor_b, WeightNames[weightmethod].c_str(), asymcut, clusteroverlay), "recreate");//
		TTree *tree = new TTree("tree", "tree");
		tree->SetMaxTreeSize(500 * 1024 * 1024); // set max tree size to 500 mb
		//------------------------------ book histograms. Need to cull this list. some are useless/redundant
		// TH1* h1 = new TH1F("h1", "pi0 E",128, 0, PT_Max);

		
		TH1 *h4 = new TH1F("h4", "Pion PT, unweighted", n_bins, PT_Min, PT_Max);
		TH1 *h5 = new TH1F("h5", "Photon Pt, unweighted", n_bins, PT_Min, PT_Max);
		TH1 *h6 = new TH1F("h6", "inv mass of gamma pair", 100, 0, 1);
		//TH1 *h7 = new TH1F("h7", "ratio of Photon/Pion pt", n_bins, PT_Min, PT_Max);
		TH1 *h8 = new TH1F("h8", "inv mass of Photon pair, smeared", 100, smeared_lower_bin_limit, smeared_upper_bin_limit);
		TH2F *h9 = new TH2F("h9", "Smeared Pion Pt vs Smeared Inv Mass", n_bins, 0, PT_Max, 100, smeared_lower_bin_limit, 2 * smeared_upper_bin_limit);
		TH1 *h10 = new TH1F("h10", "Smeared Pion PT", n_bins, PT_Min, PT_Max);
		//TH1 *h11 = new TH1F("h11", "Smeared Pion PT/Pion PT ratio", n_bins, PT_Min, PT_Max);
		
		//TH1 *h13 = new TH1F("h13", "ratio of weighted  Smeared-Pion-PT/ weighted Pion PT ", n_bins, PT_Min, PT_Max);
		//------------------------photon pT
		TH1 *h16 = new TH1F("h16", "Smeared Photon pT", n_bins, PT_Min, PT_Max);
		TH1 *h17 = new TH1F("h17", "Photon pT", n_bins, PT_Min, PT_Max);
		
		//TH1 *h24 = new TH1F("h24", "Photon pT ratio, smeared/unsmeared", n_bins, PT_Min, PT_Max);
		//TH1 *h26 = new TH1F("h26", "weighted, Photon pT ratio, smeared/unsmeared", n_bins, PT_Min, PT_Max);

		std::vector<TH1*> h2(WeightNames.size());
		std::vector<TH1*> h3(WeightNames.size());
		std::vector<TH1*> h12(WeightNames.size());
		std::vector<TH1*> h20(WeightNames.size());
		std::vector<TH1*> h21(WeightNames.size());

		std::vector<TH2F*> h18(WeightNames.size());
		std::vector<TH2F*> h27(WeightNames.size());
		std::vector<TH2F*> h28(WeightNames.size());
		std::vector<TH2F*> h29(WeightNames.size());

		std::vector<TH2F*> h30(WeightNames.size());
		std::vector<TH2F*> h31(WeightNames.size());
		std::vector<TH2F*> h32(WeightNames.size());
		std::vector<TH2F*> h33(WeightNames.size());

		std::vector<TH2F*> h34(WeightNames.size());
		std::vector<TH2F*> h35(WeightNames.size());
		std::vector<TH1F*> h35_1d(WeightNames.size());

		std::vector<TH2F*> h100(WeightNames.size());


		for (int p=0; p < WeightNames.size(); p++){
			h2[p] = new TH1D(Form("h2_%i",p),Form("Photon Pt:%s",WeightNames[p].c_str()) , n_bins, PT_Min, PT_Max); // will be weighted, is this redundant?
			h3[p] = new TH1D(Form("h3_%i",p),Form("Pion PT, weighted:%s",WeightNames[p].c_str()) , n_bins, PT_Min, PT_Max);
			h12[p] = new TH1F(Form("h12_%i",p),Form("Smeared Pion PT, weighted:%s",WeightNames[p].c_str()) , n_bins, PT_Min, PT_Max);
			h20[p] = new TH1F(Form("h20_%i",p),Form("Smeared Photon pT, weighted:%s",WeightNames[p].c_str()) , n_bins, PT_Min, PT_Max);
			h21[p] = new TH1F(Form("h21_%i",p),Form("Photon pT, weighted:%s",WeightNames[p].c_str()) , n_bins, PT_Min, PT_Max);
			h18[p] = new TH2F(Form("h18_%i",p),Form("Smeared Pion Pt vs Smeared Inv Mass, weighted:%s",WeightNames[p].c_str()) , n_bins, 0, PT_Max, 100, smeared_lower_bin_limit, 2 * smeared_upper_bin_limit);
			h27[p] = new TH2F(Form("h27_%i",p), Form("Smeared Pion Pt vs Smeared Inv Mass, weighted. cluster:%s",WeightNames[p].c_str()), n_bins, 0, PT_Max, 100, smeared_lower_bin_limit, 2 * smeared_upper_bin_limit);
			h28[p] = new TH2F(Form("h28_%i",p), Form("Smeared Pion Pt vs Smeared Inv Mass, weighted. cluster and asym cut:%s",WeightNames[p].c_str()), n_bins, 0, PT_Max, 100, smeared_lower_bin_limit, 2 * smeared_upper_bin_limit);
			h29[p] = new TH2F(Form("h29_%i",p), Form("Smeared Pion Pt vs Smeared Inv Mass, weighted. asym cut:%s",WeightNames[p].c_str()), n_bins, 0, PT_Max, 100, smeared_lower_bin_limit, 2 * smeared_upper_bin_limit);

			h30[p] = new TH2F(Form("h30_%i",p),Form("Smeared Pion Pt vs Smeared Energy, weighted:%s",WeightNames[p].c_str()) , n_bins, 0, PT_Max, 100, smeared_lower_bin_limit, 2 * smeared_upper_bin_limit);
			h31[p] = new TH2F(Form("31_%i",p), Form("Smeared Pion Pt vs Smeared Energy, weighted. cluster:%s",WeightNames[p].c_str()), n_bins, 0, PT_Max, 100, smeared_lower_bin_limit, 2 * smeared_upper_bin_limit);
			h32[p] = new TH2F(Form("h32_%i",p), Form("Smeared Pion Pt vs Smeared Energy, weighted. asym cut:%s",WeightNames[p].c_str()), n_bins, 0, PT_Max, 100, smeared_lower_bin_limit, 2 * smeared_upper_bin_limit);
			h33[p] = new TH2F(Form("h33_%i",p), Form("Smeared Pion Pt vs Smeared Energy, weighted. cluster and asym cut:%s",WeightNames[p].c_str()), n_bins, 0, PT_Max, 100, smeared_lower_bin_limit, 2 * smeared_upper_bin_limit);



			h34[p] = new TH2F(Form("h34_%i",p), Form("Smeared Pion Pt vs Smeared Inv Mass, weighted. Position Smearing:%s",WeightNames[p].c_str()), n_bins, 0, PT_Max, 100, smeared_lower_bin_limit, 2 * smeared_upper_bin_limit);

			h35[p] = new TH2F(Form("h35_%i",p), Form("Smeared Pion Pt vs Smeared Inv Mass, weighted. Blair's cuts:%s",WeightNames[p].c_str()), n_bins, 0, PT_Max, 120, 0, 1.2);
			h35_1d[p] = new TH1F(Form("h35_1d_%i",p), Form("Smeared Inv Mass, weighted. Blair's cuts:%s",WeightNames[p].c_str()), 120, 0, 1.2);



			h100[p] = new TH2F(Form("h100_%i",p), Form("Smeared Pion Pt vs Smeared Inv Mass, weighted. All Cuts+effects:%s",WeightNames[p].c_str()), n_bins, 0, PT_Max, 100, smeared_lower_bin_limit, 2 * smeared_upper_bin_limit);
		}
		/*
		TH2F *h18 = new TH2F("h18", "Smeared Pion Pt vs Smeared Inv Mass, weighted", n_bins, 0, PT_Max, 100, smeared_lower_bin_limit, 2 * smeared_upper_bin_limit);
		//h18->Sumw2();
		TH2F *h27 = new TH2F("h27", "Smeared Pion Pt vs Smeared Inv Mass, weighted. cluster", n_bins, 0, PT_Max, 100, smeared_lower_bin_limit, 2 * smeared_upper_bin_limit);

		TH2F *h28 = new TH2F("h28", "Smeared Pion Pt vs Smeared Inv Mass, weighted. cluster and asym cut", n_bins, 0, PT_Max, 100, smeared_lower_bin_limit, 2 * smeared_upper_bin_limit);

		TH2F *h29 = new TH2F("h29", "Smeared Pion Pt vs Smeared Inv Mass, weighted. asym cut", n_bins, 0, PT_Max, 100, smeared_lower_bin_limit, 2 * smeared_upper_bin_limit);
		*/
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

		double m, E, px, py, pz, pi0_px, pi0_py, pi0_E, scale_factor1, scale_factor2, smear_factor1, smear_factor2;
		double P0rest = 0.0;
		double pi0_pz = 0.0;
		double Pi0_M = 0.1349768; // 135 MeV
		double inv_mass, inv_mass_smeared;

		std::vector<int> WeightScale(WeightNames.size());
		std::vector<double> weight_function(WeightNames.size()), inv_yield(WeightNames.size());

		// std::cout << Pi0_M <<std::endl;
		tree->Branch("id", &id, "id/I");
		// tree->Branch("size",&size, "size/I");
		// tree->Branch("no",&no, "no/I");
		tree->Branch("m", &m, "m/D");
		tree->Branch("E", &E, "E/D"); // can reconstruct in root easilly
		// tree->Branch("Pt",&Pt, "Pt/D"); // can reconstruct in root easilly
		tree->Branch("px", &px, "px/D");
		tree->Branch("py", &py, "py/D");
		tree->Branch("pz", &pz, "pz/D");

		// Set up generation.  event
		Pythia pythia;								// Declare Pythia object
		pythia.readString("PromptPhoton:all = on"); // Switch on process.
		pythia.readString("ParticleDecays:allowPhotonRadiation = on");
		// pythia.readString("SoftQCD:all = on"); // Switch on process.
		// pythia.readString("Beams:eCM = 14.e3"); // 14 TeV CM energy.
		pythia.readString("111:all = pi0 -> gamma gamma");
		pythia.init(); 
		pythia.event.clear();

		for (int i = 0; i < NPions; i++){
			// std::cout << "I reached here" <<" "<< i <<std::endl;// debug line
			double azimuthal_ang = adis(gen); // generate a random angle from 0 to 2pi
			double Pt = PT_Max * pdis(gen);
			//----------------------different possible weights
			//std::cout << "EXP Weight" <<std::endl;
			weight_function[0]=exp(-Pt/0.3);//originally dividing by 0.2
			WeightScale[0]=1e+14;
			inv_yield[0] = WeightScale[0]* Pt * weight_function[0];

			//std::cout << "Power Weight" <<std::endl;
			weight_function[1]=pow(Pt,-8.14);
			WeightScale[1]=1e+5;
			inv_yield[1] = WeightScale[0]* Pt * weight_function[1];

			//std::cout << "WSHP Weight" <<std::endl;
			weight_function[2]=((1/(1+exp((Pt-t)/w)))*A/pow(1+Pt/p0,m_param)+(1-(1/(1+exp((Pt-t)/w))))*B/(pow(Pt,n)));
			WeightScale[2]=1e+14;
			inv_yield[2] = WeightScale[2]* Pt * weight_function[2];

			//std::cout << "Hagedorn Weight" <<std::endl;
			weight_function[3]=A/pow(1+Pt/p0,m_param);
			WeightScale[3]=1e+14;
			inv_yield[3] = WeightScale[3]* Pt * weight_function[3];

		//double inv_yield = WeightScale* Pt * weight_function;
			

			//std::cout << inv_yield <<std::endl;
			//h3->Fill(Pt, inv_yield); // fill pi0 pt, weighted
			//h4->Fill(Pt);						// fill pi0 pt, unweighted
			pi0_px = Pt * cos(azimuthal_ang);
			pi0_py = Pt * sin(azimuthal_ang);
			pi0_E = sqrt(Pi0_M * Pi0_M + Pt * Pt + pi0_pz * pi0_pz);
			pythia.event.append(111, 23, 0, 0, pi0_px, pi0_py, pi0_pz, pi0_E, Pi0_M); // 111 is pi0 add a partice
			// append(pid, use of particle flag, mother?, daughter?, px, py, pz, E, m)
			id = pythia.event[i].id();
			m = pythia.event[i].m();
			px = pythia.event[i].px();
			py = pythia.event[i].py();
			pz = pythia.event[i].pz();
			E = pythia.event[i].e();
			// std::cout << i <<" "<<id<< " " << m << " " << " E " << " " << E <<" " << " P " << " " << px<<" "<<py<<" "<<pz<<std::endl;
			tree->Fill();
		}

		pythia.moreDecays();
		std::cout << "I reached here" << std::endl; // debug line
		//std::ofstream mycsv2;
		//mycsv2.open(Form("pioncode/csvfiles/pT_IMass_IYield_diag_%f.csv", smear_factor_b));
		for (int i = 0; i < pythia.event.size(); i++){ // loop over all events(pions)
			if (pythia.event[i].id() == 111){ // if the ith event is a pion

				int Gamma_daughters[2] = {pythia.event[i].daughter1(), pythia.event[i].daughter2()}; // make array of daughter particles(di gamma) event ids
				double Pt = pythia.event[i].pT();
				
				//std::cout << "EXP Weight" <<std::endl;
				weight_function[0]=exp(-Pt/0.3);//originally dividing by 0.2
				WeightScale[0]=1e+14;
				inv_yield[0] = WeightScale[0]* Pt * weight_function[0];

				//std::cout << "Power Weight" <<std::endl;
				weight_function[1]=pow(Pt,-8.14);
				WeightScale[1]=1e+5;
				inv_yield[1] = WeightScale[0]* Pt * weight_function[1];

				//std::cout << "WSHP Weight" <<std::endl;
				weight_function[2]=((1/(1+exp((Pt-t)/w)))*A/pow(1+Pt/p0,m_param)+(1-(1/(1+exp((Pt-t)/w))))*B/(pow(Pt,n)));
				WeightScale[2]=1e+14;
				inv_yield[2] = WeightScale[2]* Pt * weight_function[2];

				//std::cout << "Hagedorn Weight" <<std::endl;
				weight_function[3]=A/pow(1+Pt/p0,m_param);
				WeightScale[3]=1e+14;
				inv_yield[3] = WeightScale[3]* Pt * weight_function[3];

				//----------------------different possible weights
				// I may not need this line. The variables at play are defined above the first loop where the pions are created.
				//I do this weight check in that loop. the values I assign should remain up to this point from that intial loop.
				// this version is a waste of time. It would also be good to offload this part of the code to a function. It would be a lot cleaner.
				/*
				if(weightmethod==0){
					//std::cout << "EXP Weight" <<std::endl;
					weight_function=exp(-Pt/0.3);//originally dividing by 0.2
					WeightScale=1e+14;
				}
				else if(weightmethod==1){
					//std::cout << "Power Weight" <<std::endl;
					weight_function=pow(Pt,-8.14);
					WeightScale=1e+5;
				}
				else if(weightmethod==2){
					//std::cout << "WSHP Weight" <<std::endl;
					weight_function=((1/(1+exp((Pt-t)/w)))*A/pow(1+Pt/p0,m_param)+(1-(1/(1+exp((Pt-t)/w))))*B/(pow(Pt,n)));
					WeightScale=1e+14;
				}
				else if(weightmethod==3){
				//std::cout << "Hagedorn Weight" <<std::endl;
					weight_function=A/pow(1+Pt/p0,m_param);
					WeightScale=1e+14;
				}
				else{
					std::cout << "Error:No Weight method found" <<std::endl;
				}

				double inv_yield = WeightScale* Pt * weight_function;
				*/
				if (pythia.event[Gamma_daughters[0]].id() == 22 && pythia.event[Gamma_daughters[1]].id() == 22){// check that the decays are photons
					// gammadis(gen_gamma(rdgamma()));
					
					gamma_lorentz[0] = pythia.event[Gamma_daughters[0]].p();
					gamma_lorentz[1] = pythia.event[Gamma_daughters[1]].p();
					gamma_lorentz[2] = gamma_lorentz[0] + gamma_lorentz[1];
					inv_mass = gamma_lorentz[2].mCalc();

					scale_factor1 = sqrt(pow(smear_factor_b, 2) / gamma_lorentz[0].e() + pow(smear_factor_c, 2) + pow(smear_factor_d, 2));
					scale_factor2 = sqrt(pow(smear_factor_b, 2) / gamma_lorentz[1].e() + pow(smear_factor_c, 2) + pow(smear_factor_d, 2));

					/* 
					they said  "A beam momentum spread (δp/p ≈ 2%) is quadratically subtracted from σ/μ of the fit, in order to unfolded beam momentum spread from the relative energy resolution. The Gauss function parameter of μ and energy resolution from each fit are plotted against the nominal beam energy as linearity and resolution." 
					*/

					smear_factor1 = scale_factor1 * gammadis(gen_gamma) + 1;
					smear_factor2 = scale_factor2 * gammadis(gen_gamma) + 1;

					// std::cout << "gamma gen" << " " <<gammadis(gen_gamma)<<std::endl;
					// std::cout << "pion E" << " " <<gamma_lorentz[2].e()<< " " << "smear_factor1" << " " <<smear_factor1<< " " << "smear_factor2" << " " <<smear_factor2<< " " <<std::endl;

					// position smearing. smear z phi
					// generate random numbers and move endpoint of the vector by smearing in z and p hi so ~1 tower size in z and rphi.
					// if spread over 2 towers you can fit a better position. to ~0.5 tower size.


					gamma_smeared[0] = smear_factor1 * pythia.event[Gamma_daughters[0]].p();//is px py pz recalculated? I assume so
					// std::cout << "E" << " " <<gamma_lorentz[0].e()<< " " << "smeared E" << " " <<gamma_smeared[0].e()<< " " <<std::endl; // debug, is the factor being applied?
					gamma_smeared[1] = smear_factor2 * pythia.event[Gamma_daughters[1]].p();
					
					gamma_position_smear[0]=PositionResSmear(gamma_smeared[0], posit_smearingFactor*gamma_positsmear(gen_gammapositsmear));
					gamma_position_smear[1]=PositionResSmear(gamma_smeared[1], posit_smearingFactor*gamma_positsmear(gen_gammapositsmear));

					///*

					if (gammacluster(gen_gammacluster)>coprob && clusteroverlay==1){//overlay with photon cluster 1
						//std::cout << "before cluster" << " " << gamma_smeared[0].e() <<std::endl;

						// Randomly choose an entry from the branch
						//TBranch* branch = tree->GetBranch("pz");
						//Long64_t nEntries = branch->GetEntries();
						//Long64_t entry = static_cast<Long64_t>(gRandom->Uniform(nEntries));
						// Get the value from the chosen entry
						//branch->GetEntry(entry);
						// Print or use the z_momentum value
						//std::cout << "Random z_momentum: " << z_momentum << std::endl;
						//gamma_smeared[0] = clusterPhoton(gamma_smeared[0], 2, myFunc->GetRandom())

						gamma_cluster[0] = clusterPhoton(gamma_smeared[0], 2, myFunc->GetRandom());

						gamma_cluster_asymm[0]=gamma_cluster[0];
						//gamma_smeared[0].e(gamma_smeared[0].e() + myFunc->GetRandom());
						//std::cout << "after cluster" << " " << gamma_cluster[0].e() <<std::endl;
					}
					else{
						gamma_cluster[0]=gamma_smeared[0];
						gamma_cluster_asymm[0]=gamma_cluster[0];
					}
					if (gammacluster(gen_gammacluster)>coprob && clusteroverlay==1){//overlay with photon cluster 2
						//std::cout << "before cluster" << " " << gamma_smeared[1].e() <<std::endl;

						gamma_cluster[1] = clusterPhoton(gamma_smeared[1], 2, myFunc->GetRandom());

						gamma_cluster_asymm[1]=gamma_cluster[1];
						//gamma_smeared[1].e(gamma_smeared[1].e() +myFunc->GetRandom());
						//std::cout << "after cluster" << " " << gamma_cluster[1].e() <<std::endl;
					}
					else{
						gamma_cluster[1]=gamma_smeared[1];
						gamma_cluster_asymm[1]=gamma_cluster[1];
					}//*/

					gamma_All_Cuts[2] = PositionResSmear(gamma_cluster_asymm[0], posit_smearingFactor*gamma_positsmear(gen_gammapositsmear)) + PositionResSmear(gamma_cluster_asymm[1], posit_smearingFactor*gamma_positsmear(gen_gammapositsmear));

					gamma_smeared[2] = gamma_smeared[0] + gamma_smeared[1];
					gamma_cluster[2] = gamma_cluster[0] + gamma_cluster[1];
					gamma_cluster_asymm[2] = gamma_cluster_asymm[0] + gamma_cluster_asymm[1];

					gamma_Blair_Cuts[0]=gamma_cluster_asymm[0];
					gamma_Blair_Cuts[1]=gamma_cluster_asymm[1];
					gamma_Blair_Cuts[2]=gamma_cluster_asymm[2];

					gamma_position_smear[2]=gamma_position_smear[0]+gamma_position_smear[1];
					inv_mass_smeared = gamma_smeared[2].mCalc();
					// diagnostic
					//std::cout << "photon 1, unsmr: px= " << gamma_lorentz[0].px() << " , py= " << gamma_lorentz[0].py() << " , pz= " <<gamma_lorentz[0].pz() << " , pT= " <<gamma_lorentz[0].pT() << " , E= " <<gamma_lorentz[0].e() <<std::endl;

					//std::cout << "photon 2, unsmr: px= " << gamma_lorentz[1].px() << " , py= " << gamma_lorentz[1].py() << " , pz= " <<gamma_lorentz[1].pz() << " , pT= " <<gamma_lorentz[1].pT() << " , E= " <<gamma_lorentz[1].e() <<std::endl;

					//std::cout << "pion, unsmr: px= " << gamma_lorentz[2].px() << " , py= " << gamma_lorentz[2].py() << " , pz= " <<gamma_lorentz[2].pz() << " , pT= " <<gamma_lorentz[2].pT() << " , E= " <<gamma_lorentz[2].e() << " , M= " <<gamma_lorentz[2].mCalc() <<std::endl;

					//std::cout << "smeared Pion: px= " << gamma_smeared[2].px() << " , py= " << gamma_smeared[2].py() << " , pz= " <<gamma_smeared[2].pz() << " , pT= " <<gamma_smeared[2].pT() << " , E= " <<gamma_smeared[2].e()  << " , M= " << gamma_smeared[2].mCalc() <<std::endl;

					//std::cout << "Posit smeared Pion: px= " << gamma_position_smear[2].px() << " , py= " << gamma_position_smear[2].py() << " , pz= " <<gamma_position_smear[2].pz() << " , pT= " <<gamma_position_smear[2].pT() << " , E= " <<gamma_position_smear[2].e()  << " , M= " << gamma_position_smear[2].mCalc() <<std::endl;


					/*
					if(abs(gamma_smeared[0].e()-gamma_smeared[1].e())/(gamma_smeared[0].e()+gamma_smeared[1].e())>0.8 &&asymcut==1){//asymmetry cut
					//std::cout << "Asymmetry Cut" << " " << abs(gamma_smeared[0].e()-gamma_smeared[1].e())/(gamma_smeared[0].e()+gamma_smeared[1].e())<<std::endl;
					// if I am to save both, maybe filling here would be appropriate.
					//gamma_asymm=;
					//gamma_cluster_asymm;
					//h28->Fill(gamma_smeared[2].pT(), inv_mass_smeared, inv_yield);
					//h29->Fill(gamma_smeared[2].pT(), inv_mass_smeared, inv_yield);
					//continue;
					}*/
					for (int p=0; p < WeightNames.size(); p++){//weighted histograms
						//h2[p]->Fill(pythia.event[Gamma_daughters[0]].pT(), inv_yield[p]);
						//h2[p]->Fill(pythia.event[Gamma_daughters[1]].pT(), inv_yield[p]);

						h20[p]->Fill(gamma_smeared[0].pT(), inv_yield[p]);
						h20[p]->Fill(gamma_smeared[1].pT(), inv_yield[p]);

						h21[p]->Fill(gamma_lorentz[0].pT(), inv_yield[p]);
						h21[p]->Fill(gamma_lorentz[1].pT(), inv_yield[p]);
						
						h3[p]->Fill(Pt, inv_yield[p]); // fill pion pt, weighted
						h12[p]->Fill(gamma_smeared[2].pT(), inv_yield[p]);

						h18[p]->Fill(gamma_smeared[2].pT(), inv_mass_smeared, inv_yield[p]);//
						h27[p]->Fill(gamma_cluster[2].pT(), gamma_cluster[2].mCalc(), inv_yield[p]);

						h30[p]->Fill(gamma_smeared[2].pT(), gamma_smeared[2].e(), inv_yield[p]);//
						h31[p]->Fill(gamma_cluster[2].pT(), gamma_cluster[2].e(), inv_yield[p]);//

						
						h34[p]->Fill(gamma_position_smear[2].pT(), gamma_position_smear[2].mCalc(), inv_yield[p]);//position smearing

						if(abs(gamma_smeared[0].e()-gamma_smeared[1].e())/(gamma_smeared[0].e()+gamma_smeared[1].e())<asymval &&asymcut==1){//asymmetry cut
							//std::cout << "Asymmetry Cut" << " " << abs(gamma_smeared[0].e()-gamma_smeared[1].e())/(gamma_smeared[0].e()+gamma_smeared[1].e())<<std::endl;
							// if I am to save both, maybe filling here would be appropriate.
							h29[p]->Fill(gamma_smeared[2].pT(), gamma_smeared[2].mCalc(), inv_yield[p]);// asymm
							h28[p]->Fill(gamma_cluster_asymm[2].pT(), gamma_cluster_asymm[2].mCalc(), inv_yield[p]);//
							//continue;
							h32[p]->Fill(gamma_smeared[2].pT(), gamma_smeared[2].e(), inv_yield[p]);//
							h33[p]->Fill(gamma_cluster_asymm[2].pT(), gamma_cluster_asymm[2].e(), inv_yield[p]);//


							if(DeltaRcut(gamma_Blair_Cuts[0], gamma_Blair_Cuts[1], DeltaRcut_MAX)==false && pTCut(gamma_Blair_Cuts[0], pt1cut)==true && pTCut(gamma_Blair_Cuts[2], pt2cut)==true){
								h35[p]->Fill(gamma_Blair_Cuts[2].pT(), gamma_Blair_Cuts[2].mCalc(), inv_yield[p]);
								h35_1d[p]->Fill(gamma_Blair_Cuts[2].mCalc(), inv_yield[p]);
							}

							h100[p]->Fill(gamma_All_Cuts[2].pT(), gamma_All_Cuts[2].mCalc(), inv_yield[p]);
						}
						//std::cout << "smeared energy " << gamma_smeared[2].e() <<". clustered energy " << gamma_cluster[2].e() << " . ratio c/s "<< gamma_cluster[2].e()/gamma_smeared[2].e()<< std::endl;
					}

					// std::cout << "inv mass" << " " <<inv_mass<<std::endl;

					h4->Fill(Pt);			// fill pion pt, unweighted
					h6->Fill(inv_mass);
					h8->Fill(inv_mass_smeared);
					h9->Fill(gamma_smeared[2].pT(), inv_mass_smeared);//change  to smeared pion pT. was previously unsmeared: pythia.event[i].pT()
					h10->Fill(gamma_smeared[2].pT());	
					h17->Fill(gamma_lorentz[0].pT());//unsmeared energy spectrum
					h17->Fill(gamma_lorentz[1].pT());
					h16->Fill(gamma_smeared[0].pT());//smeared photon energy spectrum
					h16->Fill(gamma_smeared[1].pT());


					// mass_pt_map.insert[]
					//std::cout << "event number" << " " << i << " " << "smeared pT" << " " << std::scientific << gamma_smeared[2].pT() << " " << "Smeared Mass" << " " << inv_mass_smeared << " " << "inv yield" << " " << std::scientific << inv_yield << std::endl;
					
					/*
					if(i==0){
						mycsv2 << "event number" << "," <<  "smeared pT" <<  "," << "Smeared Mass" << "," << "inv yield" << "\n";
						mycsv2 << i << "," << std::scientific << gamma_smeared[2].pT() << "," << inv_mass_smeared << "," << std::scientific << inv_yield << "\n";
						//mycsv2 << "\n";
					}
					else{
						//std::cout << "I reached here" <<" "<< i <<std::endl;// debug line
						//mycsv2 << "event number" << "," <<  "smeared pT" <<  "," << "Smeared Mass" << "," << "inv yield" << "\n";
						mycsv2 << i << "," << std::scientific << gamma_smeared[2].pT() << "," << inv_mass_smeared << "," << std::scientific << inv_yield << "\n";
					}*/

					
					///*
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

						// if(pythia.event[j].id() == 22){// check that the jth event is a photon

						id = pythia.event[j].id();
						m = pythia.event[j].m();
						px = pythia.event[j].px();
						py = pythia.event[j].py();
						pz = pythia.event[j].pz();
						E = pythia.event[j].e();
						// std::cout << i <<" "<<id<< " " << m << " " << " E " << " " << E <<" " << " P " << " " << px<<" "<<py<<" "<<pz<<std::endl;
						tree->Fill();
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
		
		float Smeared_Mean_array[n_bins], Smeared_Variance_array[n_bins], nbins_array[n_bins]; 
		//mycsv2.close();
		/*
		std::ofstream mycsv;
		mycsv.open(Form("pioncode/csvfiles/Inv_Mass_mean_variance%f.csv", smear_factor_b));
		for (int i = 1; i < n_bins + 1; i++){
			TH1 *htemp1 = new TH1D("htemp1", "temp1", n_bins, PT_Min, PT_Max); // unweighted

			for (int j = 0; j < mass_pt_map[i].size(); j++)
			{
				htemp1->Fill(mass_pt_map[i][j]);
			}


			Smeared_Mean_array[i] = htemp1->GetMean();
			nbins_array[i] = i;

			Smeared_Variance_array[i] = pow(htemp1->GetStdDev(), 2);

			mycsv << i << "," << Smeared_Mean_array[i] << "," << Smeared_Variance_array[i] << "," << Smeared_Variance_array[i] / Smeared_Mean_array[i] << "\n";
			//std::cout << "Smeared Mean = " << Smeared_Mean_array[i] << " , " << "Smeared Variance = " << Smeared_Variance_array[i] << std::endl;
			// std::cout << "Smeared Mean = " << Smeared_Mean_array[i] << " , " << "Smeared Variance = " << "blank" <<std::endl;
			delete htemp1;
		}
		mycsv.close();*/


		//h7->Divide(h2, h3);
		//h11->Divide(h10, h4);
		//h13->Divide(h12, h3);
		//h24->Divide(h16, h17);
		//h26->Divide(h20, h21);
		double_t realtime = timer.RealTime();
		double_t cputime = timer.CpuTime();

		printf("real time =%f\n", realtime);
		printf("CPU time=%f\n", cputime);

		output->Write();
		output->Close();
		//clean up
		//delete tree; 
		//delete h2;//,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12,h13,h16,h17,h18,h20,h21,h24,h26,h27,h28,h29;
		delete output;
		
		
	}
	return 0;
}

TF1* ChooseSpectrumFunction(int weightmethod, int PT_Min, int PT_Max){
	TF1 *myFunc = nullptr;
	if(weightmethod==0){
		myFunc = new TF1("myFunc", [](double *x, double *par){
			return par[0] * TMath::Exp(-x[0]/0.3);
		}, PT_Min, PT_Max, 1);
		Double_t initialParameters[1] = {1.0};
    	myFunc->SetParameters(initialParameters);

	}
	else if(weightmethod==1){
		myFunc = new TF1("myFunc", [](double *x, double *par){
			return par[0] * pow(x[0],-8.14);
			}, PT_Min, PT_Max, 1);
		Double_t initialParameters[1] = {1.0};
    	myFunc->SetParameters(initialParameters);
	}
	else if(weightmethod==2){
		myFunc = new TF1("myFunc", [](double *x, double *par){
			return ((1/(1+exp((x[0]-par[0])/par[1])))*par[2]/pow(1+x[0]/par[3],par[4])+(1-(1/(1+exp((x[0]-par[0])/par[1]))))*par[5]/(pow(x[0],par[6])));
			}, PT_Min, PT_Max, 7);
		myFunc->SetParameters(4.5,  0.114, 229.6, 1.466, 10.654, 14.43, 8.1028);
		//t,w,A,p0,m_param,B,n
	}
	else if(weightmethod==3){
		//std::cout << "Hagedorn Weight" <<std::endl;
		myFunc = new TF1("myFunc", [](double *x, double *par){
			return (par[0]/pow(1+x[0]/par[1],par[2]) );
			}, PT_Min, PT_Max, 3);
		myFunc->SetParameters( 229.6, 1.466, 10.654);
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

Pythia8::Vec4 PositionResSmear(Pythia8::Vec4 photon, double smearingFactor) {//, Pythia8::Rndm* rndm

	//calculate the pT to keep the length of the four momenta the same. then smear the px and py to match.
	//*
	double xsmear = smearingFactor + 1;
	double px_smear = photon.px()*xsmear;

	double ysmear = smearingFactor + 1;
	double py_smear = photon.py()*ysmear;

	double zsmear = smearingFactor + 1;
	double pz_smear = photon.pz()*zsmear;


    // Keep the energy constant
    double energy = photon.e();
	double E_New = sqrt(px_smear*px_smear +py_smear*py_smear +pz_smear*pz_smear);
	double energyscale=energy/E_New;

	double xnew=energyscale*px_smear;
	double ynew=energyscale*py_smear;
	double znew=energyscale*pz_smear;

   	//std::cout << "E_New/E_old" << sqrt(xnew*xnew +ynew*ynew +znew*znew)/energy << std::endl;
    Pythia8::Vec4 smearedPhoton(xnew, ynew, znew, sqrt(xnew*xnew +ynew*ynew +znew*znew));
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

