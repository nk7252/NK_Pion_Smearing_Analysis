// Headers and Namespaces.
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <chrono>
// #include <numbers> //std::numbers
#include <cmath>
#include <random>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TStopwatch.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "Pythia8/Pythia.h" // Include Pythia headers.
using namespace Pythia8;	// Let Pythia8:: be implicit.

int main()
{ // Begin main program.
	// auto start = std::chrono::steady_clock::now();
	// std::cout << "Start Time (sec) = " << start.count() << std::endl;
	// int prefactor =3 ;
	TStopwatch timer;
	timer.Start();
	int NPions = 1 * 1000000;
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
	int weightmethod = 2;//0=exp,1=power,2=wshp
	std::vector<std::string> WeightNames = {"EXP", "POWER", "WSHP"};
	//-----------------------------------
	int asymcut=0;//apply asymm cut.
	int clusteroverlay = 0;//overlayed cluster check

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

	for (int smear_factor_itt = 9; smear_factor_itt < 9 + 1; smear_factor_itt++)
	{// originally int smear_factor_itt = 0; smear_factor_itt < 24 + 1; smear_factor_itt++
	// only want .155
		float smear_factor_basevalue = 0.065; // I used 1.6% + 12.7%/sqrt(E) fig 22, but that is from a special beam cross section config. trying with fig 24 data i.e 2.8% + 15.5%
		//--------------------preliminaries to read from root
		float smear_factor_b = smear_factor_basevalue + 0.01 * smear_factor_itt;
		//////////////////////New//0.155 loop from twice test beam data paramaterization to half? this is 15.5% from https://arxiv.org/pdf/1704.01461.pdf fig 24b so going from 6.5% to 30.5%// need 24 steps for 1% diff each
		//////////////////////OLD//0.127 loop from twice test beam data paramaterization to half? this is 12.7% from https://arxiv.org/pdf/1704.01461.pdf fig 22b so going from 6.35% to 25.4%

		TFile *output = new TFile(Form("pioncode/rootfiles/Pi0FastMC_%f_%s.root", smear_factor_b, WeightNames[weightmethod].c_str()), "recreate");
		TTree *tree = new TTree("tree", "tree");
		tree->SetMaxTreeSize(500 * 1024 * 1024); // set max tree size to 500 mb

		// TH1* h1 = new TH1F("h1", "pi0 E",128, 0, PT_Max);
		TH1 *h2 = new TH1D("h2", "gamma Pt", n_bins, PT_Min, PT_Max); // will be weighted
		TH1 *h3 = new TH1D("h3", "pi0 Pt", n_bins, PT_Min, PT_Max);
		TH1 *h4 = new TH1F("h4", "pi0 Pt, unweighted", n_bins, PT_Min, PT_Max);
		TH1 *h5 = new TH1F("h5", "gamma Pt, unweighted", n_bins, PT_Min, PT_Max);
		TH1 *h6 = new TH1F("h6", "inv mass of gamma pair", 100, 0, 1);
		TH1 *h7 = new TH1F("h7", "ratio of gamma/pi0 pt", n_bins, PT_Min, PT_Max);
		TH1 *h8 = new TH1F("h8", "inv mass of gamma pair, smeared", 100, smeared_lower_bin_limit, smeared_upper_bin_limit);
		TH2F *h9 = new TH2F("h9", "Smeared Pion Pt vs Smeared Inv Mass", n_bins, 0, PT_Max, 100, smeared_lower_bin_limit, 2 * smeared_upper_bin_limit);
		TH1 *h10 = new TH1F("h10", "Smeared Pion PT", n_bins, PT_Min, PT_Max);
		TH1 *h11 = new TH1F("h11", "Smeared Pion PT/Pion PT ratio", n_bins, PT_Min, PT_Max);
		TH1 *h12 = new TH1F("h12", "Smeared Pion PT, weighted", n_bins, PT_Min, PT_Max);
		TH1 *h13 = new TH1F("h13", "ratio of weighted  Smeared-Pion-PT/ weighted Pion PT ", n_bins, PT_Min, PT_Max);
		TH1 *h14 = new TH1F("h14", "ratio of weighted and unweighted ratio ", n_bins, PT_Min, PT_Max);
		TH1 *h15 = new TH1F("h15", "ratio of unweighted and weighted ratio ", n_bins, PT_Min, PT_Max);
		//
		TH1 *h16 = new TH1F("h16", "Smeared Photon pT", n_bins, PT_Min, PT_Max);
		TH1 *h17 = new TH1F("h17", "Photon pT", n_bins, PT_Min, PT_Max);
		TH1 *h20 = new TH1F("h20", "Smeared Photon pT, weighted", n_bins, PT_Min, PT_Max);
		TH1 *h21 = new TH1F("h21", "Photon pT, weighted", n_bins, PT_Min, PT_Max);
		TH1 *h24 = new TH1F("h24", "Photon pT ratio, smeared/unsmeared", n_bins, PT_Min, PT_Max);
		TH1 *h26 = new TH1F("h26", "weighted, Photon pT ratio, smeared/unsmeared", n_bins, PT_Min, PT_Max);
		TH1 *h28 = new TH1F("h28", "weighted/unweighted ratio of ratios, photons", n_bins, PT_Min, PT_Max);
		TH2F *h18 = new TH2F("h18", "Smeared Pion Pt vs Smeared Inv Mass, weighted", n_bins, 0, PT_Max, 100, smeared_lower_bin_limit, 2 * smeared_upper_bin_limit);
		h18->Sumw2();
		// things to add probability to add an exponentially scaled (small) energy
		//--------------------set up random number generation
		std::random_device rd;		// generate a random number to seed random generation
		std::random_device rdgamma; // generate a random number to seed random generation of daughter gamma for smearing
		std::random_device rdgammacluster;// random number to test clustering
		// std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd
		std::mt19937_64 gen(rd());			  // mersenne_twister_engine 64 bit seeded with rd
		std::mt19937_64 gen_gamma(rdgamma()); // mersenne_twister_engine 64 bit seeded with rdgamma for gamma smearing
		std::mt19937_64 gen_gammacluster(rdgammacluster());//generate random number to test clustering 
		// std::ranlux48 gen(rd()); // ranlux48 seeded with rd
		//--------------------generate random: momentum
		// std::uniform_real_distribution<> pdis(0.0,PT_Max);
		std::normal_distribution<double> gammadis(0.0, 1.0);  // generate normal distribution for gamma smearing, mean zero, variance 1
		std::uniform_real_distribution<> gammacluster(0, 1.0); // random probability to have photons clustered together. if above certain value, add a random photon from the smeared+weighted photon spectrum. cutoff will be tuned.
		std::uniform_real_distribution<> pdis(PT_ratio, 1.0); // alternative scheme with min PT to avoid power law complications.
		//--------------------generate random: angle
		// double tpi=2*std::numbers::pi;
		std::uniform_real_distribution<> adis(0.0, 2 * M_PI);
		int id, size, no, WeightScale;
		// std::vector<int> vec_id;
		// std::vector<double> vec_E;
		// std::vector<double> vec_Pi0_Pt;
		// std::vector<double> vec_P;
		// std::vector<double> vec_inv_mass;
		Vec4 gamma_lorentz[3];
		Vec4 gamma_smeared[3];
		double m, E, px, py, pz, pi0_px, pi0_py, pi0_E, scale_factor1, scale_factor2, smear_factor1, smear_factor2;
		double P0rest = 0.0;
		double pi0_pz = 0.0;
		double Pi0_M = 0.1349768; // 135 MeV
		double inv_mass, inv_mass_smeared,weight_function;

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
		pythia.init(); // Initialize; incoming pp beams is default.
		// pythia.particleData.mayDecay(111,true);
		// pythia.particleData.mayDecay(22,false);
		// pythia.particleData.mayDecay(11,false);
		// pythia.particleData.mayDecay(-11,false);
		//  Generate event(s).
		// pythia.next(); // Generate an(other) event. Fill event record.
		pythia.event.clear();
		// bool readytoappend = true;

		for (int i = 0; i < NPions; i++){
			// std::cout << "I reached here" <<" "<< i <<std::endl;// debug line
			double azimuthal_ang = adis(gen); // generate a random angle from 0 to 2pi
			double Pt = PT_Max * pdis(gen);
			//----------------------different possible weights
			if(weightmethod==0){
				//std::cout << "EXP Weight" <<std::endl;
				double weight_function=exp(-Pt/0.3);//originally dividing by 0.2
				int WeightScale=1e+20;
			}
			else if(weightmethod==1){
				//std::cout << "Power Weight" <<std::endl;
				double weight_function=pow(Pt,-8.14);
				int WeightScale=1e+5;
			}
			else if(weightmethod==2){
				//std::cout << "WSHP Weight" <<std::endl;
				double weight_function=((1/(1+exp((Pt-t)/w)))*A/pow(1+Pt/p0,m_param)+(1-(1/(1+exp((Pt-t)/w))))*B/(pow(Pt,n)));
				int WeightScale=1e+14;
			}
			else{
				std::cout << "Error:No Weight method found" <<std::endl;
			}
			double inv_yield = WeightScale* Pt * weight_function;

			h3->Fill(Pt, inv_yield); // fill pi0 pt, weighted
			h4->Fill(Pt);						// fill pi0 pt, unweighted
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
		// pythia.forceHadronLevel();
		// Event& event = pythia.event;
		//  force decays
		// pythia.resetDecays();
		// pythia.decay();
		pythia.moreDecays();
		std::cout << "I reached here" << std::endl; // debug line
		std::ofstream mycsv2;
		mycsv2.open(Form("pioncode/csvfiles/pT_IMass_IYield_diag_%f.csv", smear_factor_b));
		for (int i = 0; i < pythia.event.size(); i++)
		{ // loop over all events(pions)
			if (pythia.event[i].id() == 111)
			{ // if the ith event is a pion

				int Gamma_daughters[2] = {pythia.event[i].daughter1(), pythia.event[i].daughter2()}; // make array of daughter particles(di gamma) event ids
				double Pt = pythia.event[i].pT();
				//----------------------different possible weights
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
				else{
					std::cout << "Error:No Weight method found" <<std::endl;
				}

				double inv_yield = WeightScale* Pt * weight_function;

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

					gamma_smeared[0] = smear_factor1 * pythia.event[Gamma_daughters[0]].p();
					// std::cout << "E" << " " <<gamma_lorentz[0].e()<< " " << "smeared E" << " " <<gamma_smeared[0].e()<< " " <<std::endl; // debug, is the factor being applied?
					gamma_smeared[1] = smear_factor2 * pythia.event[Gamma_daughters[1]].p();


					if(abs(gamma_smeared[0].e()-gamma_smeared[1].e())/(gamma_smeared[0].e()+gamma_smeared[1].e())>0.8 &&asymcut==1){//asymmetry cut
					//std::cout << "Asymmetry Cut" << " " << abs(gamma_smeared[0].e()-gamma_smeared[1].e())/(gamma_smeared[0].e()+gamma_smeared[1].e())<<std::endl;
					continue;
					}

					gamma_smeared[2] = gamma_smeared[0] + gamma_smeared[1];
					inv_mass_smeared = gamma_smeared[2].mCalc();
					// std::cout << "inv mass" << " " <<inv_mass<<std::endl;
					h6->Fill(inv_mass);
					h8->Fill(inv_mass_smeared);
					h9->Fill(gamma_smeared[2].pT(), inv_mass_smeared);//change  to smeared pion pT. was previously unsmeared: pythia.event[i].pT()
					h10->Fill(gamma_smeared[2].pT());
					h12->Fill(gamma_smeared[2].pT(), Pt * weight_function);
					h17->Fill(gamma_lorentz[0].pT());//unsmeared energy spectrum
					h17->Fill(gamma_lorentz[1].pT());
					h16->Fill(gamma_smeared[0].pT());//smeared photn energy spectrum
					h16->Fill(gamma_smeared[1].pT());
					///*
					if (gammacluster(gen_gammacluster)>0.8 && clusteroverlay==1){//overlay with photon cluster 2
					std::cout << "before cluster" << " " << gamma_smeared[0].e() <<std::endl;
					gamma_smeared[0].e(gamma_smeared[0].e() + h16->GetRandom());
					std::cout << "after cluster" << " " << gamma_smeared[0].e() <<std::endl;
					}
					if (gammacluster(gen_gammacluster)>0.8 && clusteroverlay==1){//overlay with photon cluster 2
					std::cout << "before cluster" << " " << gamma_smeared[1].e() <<std::endl;
					gamma_smeared[1].e(gamma_smeared[1].e() +h16->GetRandom());
					std::cout << "after cluster" << " " << gamma_smeared[1].e() <<std::endl;
					}//*/
					h20->Fill(gamma_smeared[0].pT(), Pt * weight_function);
					h21->Fill(gamma_lorentz[0].pT(), Pt * weight_function);
					h20->Fill(gamma_smeared[1].pT(), Pt * weight_function);
					h21->Fill(gamma_lorentz[1].pT(), Pt * weight_function);
					h18->Fill(gamma_smeared[2].pT(), inv_mass_smeared, inv_yield);//change x to smeared pion pT. was previously unsmeared: pythia.event[i].pT()
					// mass_pt_map.insert[]
					//std::cout << "event number" << " " << i << " " << "smeared pT" << " " << std::scientific << gamma_smeared[2].pT() << " " << "Smeared Mass" << " " << inv_mass_smeared << " " << "inv yield" << " " << std::scientific << inv_yield << std::endl;
					
					
					if(i==0){
						mycsv2 << "event number" << "," <<  "smeared pT" <<  "," << "Smeared Mass" << "," << "inv yield" << "\n";
						mycsv2 << i << "," << std::scientific << gamma_smeared[2].pT() << "," << inv_mass_smeared << "," << std::scientific << inv_yield << "\n";
						//mycsv2 << "\n";
					}
					else{
						//std::cout << "I reached here" <<" "<< i <<std::endl;// debug line
						//mycsv2 << "event number" << "," <<  "smeared pT" <<  "," << "Smeared Mass" << "," << "inv yield" << "\n";
						mycsv2 << i << "," << std::scientific << gamma_smeared[2].pT() << "," << inv_mass_smeared << "," << std::scientific << inv_yield << "\n";
					}

					
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
						h2->Fill(pythia.event[j].pT(), Pt * weight_function);
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
		mycsv2.close();
		std::ofstream mycsv;

		mycsv.open(Form("pioncode/csvfiles/Inv_Mass_mean_variance%f.csv", smear_factor_b));
		for (int i = 1; i < n_bins + 1; i++)
		{
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
		mycsv.close();
		// TCanvas* c1 = new TCanvas("c1","A Simple Graph Example",200,10,700,500);
		// TGraph* g1 = new TGraph(n_bins, nbins_array, Smeared_Mean_array);
		// auto g2 = new TGraph(64, nbins_array,Smeared_Variance_array);
		// g1->SetTitle("Mean;X title;Y title");
		// g1->Draw("AC*");
		// g2->SetTitle("Var;X title;Y title");
		// g2->Draw("AC*");
		// g1->Write("MyGraph1");
		// g2->Write("MyGraph2");
		// delete g1;
		// delete g2;
		// h10->Fit("gaus","Q"); // Q quiet mode, https://root.cern.ch/root/htmldoc/guides/users-guide/FittingHistograms.html

		// h11->Fit("gaus","Q");

		// for(int i=0; i<n_bins;i++){
		// auto mean_graph (n_bins, Smeared_Mean);
		// auto variance_graph (n_bins, Smeared_Variance);
		// };

		/*
		std::cout << i <<" "<<id<< " " << m << " " << " E " << " " << E <<" " << " P " << " " << px<<" "<<py<<" "<<pz<<std::endl;
		*/
		h7->Divide(h2, h3);
		h11->Divide(h10, h4);
		h13->Divide(h12, h3);
		h14->Divide(h13, h11);
		h15->Divide(h11, h13);
		//
		h24->Divide(h16, h17);
		// h25->Divide(h18,h19);
		h26->Divide(h20, h21);
		// h27->Divide(h22,h23);
		h28->Divide(h26, h24);
		// h29->Divide(h27,h25);
		// h8->Divide(h3,h2);
		double_t realtime = timer.RealTime();
		double_t cputime = timer.CpuTime();

		printf("real time =%f\n", realtime);
		printf("CPU time=%f\n", cputime);

		output->Write();
		output->Close();
	}
	return 0;
}

