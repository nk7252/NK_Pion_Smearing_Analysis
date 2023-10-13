// Headers and Namespaces.
#define _USE_MATH_DEFINES
#include <iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <random>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TStopwatch.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "Math/Vector4D.h"
#include "TLorentzVector.h"

void simplephotongen(){ // Begin main program.
	// auto start = std::chrono::steady_clock::now();
	// std::cout << "Start Time (sec) = " << start.count() << std::endl;
	// int prefactor =3 ;
	TStopwatch timer;
	timer.Start();
	int NPhotons = 1 * 1000000;
	int PT_Max = 64; // 65
	int PT_Min = 0;	 // cross check to elimate power law problems
	double PT_ratio = PT_Min / PT_Max;
	int n_bins = PT_Max;

	//--------------------3 types of weighting: power law, exponential, woods saxon+hagedorn+power law
	//-------------------- power law param
	float plaw_param = 7.1;
	//-------------------- exp param
	float exp_param = 200;
	//--------------------Alternative paramaterization, woods saxon+hagedorn+power law
	double t = 4.5;
	double w = 0.114;
	double A = 229.6;
	double B = 14.43;
	double n = 8.1028;
	double m_param = 10.654;
	double p0 = 1.466;

	// smearing params
	float smear_factor_d = 0.02;  // 0.02;// test trying to include the beam momentum resolution.. will set to zero for now
	float smear_factor_c = 0.028; // first parameter in test beam parameterization?
	float smear_factor_b = 0.155;

	//TFile *output = new TFile(pioncode / rootfiles / simplephotoncheck.root, "recreate");
	//TTree *tree = new TTree("tree", "tree");

	TH1 *h1 = new TH1D("h1", "gamma Pt, flat", n_bins, PT_Min, PT_Max);
	//weighted
	TH1 *h2 = new TH1D("h2", "gamma Pt, exp weight", n_bins, PT_Min, PT_Max);
	TH1 *h3 = new TH1D("h3", "gamma Pt, power weight", n_bins, PT_Min, PT_Max);
	TH1 *h4 = new TH1D("h4", "gamma Pt, w.s.h.p. weight", n_bins, PT_Min, PT_Max);
	//weighted+smeared
	TH1 *h5 = new TH1D("h5", "gamma Pt smeared, exp weight", n_bins, PT_Min, PT_Max);
	TH1 *h6 = new TH1D("h6", "gamma Pt smeared, power weight", n_bins, PT_Min, PT_Max);
	TH1 *h7 = new TH1D("h7", "gamma Pt smeared, w.s.h.p. weight", n_bins, PT_Min, PT_Max);
	//weighted+smeared/weighted ratio
	TH1 *h8 = new TH1D("h8", "gamma Pt smeared+weighted/weighted, exp weight", n_bins, PT_Min, PT_Max);
	TH1 *h9 = new TH1D("h9", "gamma Pt smeared+weighted/weighted, power weight", n_bins, PT_Min, PT_Max);
	TH1 *h10 = new TH1D("h10", "gamma Pt smeared+weighted/weighted, w.s.h.p. weight", n_bins, PT_Min, PT_Max);

	std::random_device rd;		// generate a random number to seed random generation
	std::random_device rdgamma; // generate a random number to seed random generation of daughter gamma for smearing
	// std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd
	std::mt19937_64 gen(rd());							 // mersenne_twister_engine 64 bit seeded with rd
	std::mt19937_64 gen_gamma(rdgamma());				 // mersenne_twister_engine 64 bit seeded
	std::normal_distribution<double> gammadis(0.0, 1.0); // generate normal distribution for gamma smearing, mean zero, variance 1
	std::uniform_real_distribution<> pdis(PT_ratio, 1.0);
	std::uniform_real_distribution<> adis(0.0, 2 * M_PI);
	float pz =0.0;
	for (int i = 0; i < NPhotons; i++){
		float azimuthal_ang = adis(gen); // generate a random angle from 0 to 2pi
		float Pt = PT_Max * pdis(gen);
		h1->Fill(Pt);
		h1->Fill(Pt);
		//------------------------------smearing photons
		float scale_factor1 = sqrt(pow(smear_factor_b, 2) / Pt + pow(smear_factor_c, 2) + pow(smear_factor_d, 2));
		//float scale_factor2 = sqrt(pow(smear_factor_b, 2) / Pt + pow(smear_factor_c, 2) + pow(smear_factor_d, 2));// these are naturally the same if I say pz=0
		float smear_factor1 = scale_factor1 * gammadis(gen_gamma) + 1;
		float smear_factor2 = scale_factor1 * gammadis(gen_gamma) + 1;
		float smeared_Pt1=Pt*smear_factor1;
		float smeared_Pt2=Pt*smear_factor2;

		//------------------------------find energies/ Pt
		double weight_function1 = exp(-Pt/0.200 );// exponential. scaled it to 200
		//note hagedorn should be dominant up to ~5 GeV from wf3. we are using an exponential instead
		double weight_function2 = 1 / (pow(Pt, 7.1)); //power law
		// power law is becomes important around 5 GeV
		double weight_function3 = ((1 / (1 + exp((Pt - t) / w))) * A / pow(1 + Pt / p0, m_param) + (1 - (1 / (1 + exp((Pt - t) / w)))) * B / (pow(Pt, n)));
		h2->Fill(Pt, weight_function1);
		h2->Fill(Pt, weight_function1);
		h5->Fill(smeared_Pt1, weight_function1);
		h5->Fill(smeared_Pt2, weight_function1);

		h3->Fill(Pt, weight_function2);
		h3->Fill(Pt, weight_function2);
		h6->Fill(smeared_Pt1, weight_function2);
		h6->Fill(smeared_Pt2, weight_function2);

		h4->Fill(Pt, Pt*weight_function3);
		h4->Fill(Pt, Pt*weight_function3);
		h7->Fill(smeared_Pt1, Pt * weight_function3);
		h7->Fill(smeared_Pt2, Pt * weight_function3);
	}


	TCanvas* c1 = new TCanvas("c1", "c1", 800, 800);//unsmeared pTs
	// c1->SetFillColor(42);
	c1->Divide(2, 2); // nx, ny
	c1->cd(1);
	h1->Draw("colz");
	c1->cd(2);
	gPad->SetLogy();
	h2->Draw("colz");
	c1->cd(3);
	gPad->SetLogy();
	h3->Draw("colz");
	c1->cd(4);
	gPad->SetLogy();
	h4->Draw("colz");

	TCanvas* c3 = new TCanvas("c3", "c3", 900, 300);//smeared Pts
	c3->Divide(3, 1); // nx, ny
	c3->cd(1);
	gPad->SetLogy();
	h5->Draw("colz");
	c3->cd(2);
	gPad->SetLogy();
	h6->Draw("colz");
	c3->cd(3);
	gPad->SetLogy();
	h7->Draw("colz");

	TCanvas* c2 = new TCanvas("c2", "c2", 900, 300);//ratios
	// c1->SetFillColor(42);
	c2->Divide(3, 1); // nx, ny
	c2->cd(1);
	//h8->Clone("h5");
	//TH1D* h8 = (TH1D*)h5->Clone("h8");
	//h5->Scale(pow(10, 14));
	//h2->Scale(pow(10, 14));
	//h2->SetBinContent(0, 344);
	//h2->SetBinContent(65,2.51738E-122);
	//h8->Divide(h2);
	double h8binconttemp;
	for(int j=1;j<n_bins+1;j++){
		h8binconttemp=h5->GetBinContent(j)/h2->GetBinContent(j);
		h8->SetBinContent(j,h8binconttemp);
		std::cout << h5->GetBinContent(j) << " Divided by " << h2->GetBinContent(j) << " gives " << h8->GetBinContent(j) << std::endl;
	}
	//h8->Scale(pow(10, 14));
	gPad->SetLogy();
	h8->Draw("colz");
	//h8->DrawNormalized("colz");
	c2->cd(2);
	//h9->Clone("h6");
	h9->Divide(h6,h3);
	h9->Draw("colz");
	c2->cd(3);
	//h10->Clone("h7");
	h10->Divide(h7,h4);
	h10->Draw("colz");

	

	double_t realtime = timer.RealTime();
	double_t cputime = timer.CpuTime();

	printf("real time =%f\n", realtime);
	printf("CPU time=%f\n", cputime);
}
