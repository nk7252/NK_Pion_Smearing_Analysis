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

int main() {   // Begin main program.  	
	//auto start = std::chrono::steady_clock::now();
	//std::cout << "Start Time (sec) = " << start.count() << std::endl;
	//int prefactor =3 ;
	TStopwatch timer;
	timer.Start();
	int NPhotons=1*1000000;
	int PT_Max=64;//65
	int PT_Min=0;// cross check to elimate power law problems
	double PT_ratio=PT_Min/PT_Max;
    int n_bins=PT_Max;

    //--------------------3 types of weighting: power law, exponential, woods saxon+hagedorn+power law
    //-------------------- power law param
    float plaw_param=7.1;
    //-------------------- exp param
    float exp_param=7.1;
    //--------------------Alternative paramaterization, woods saxon+hagedorn+power law
	double t=4.5;
	double w=0.114;
	double A=229.6;
	double B=14.43;
	double n=8.1028;
	double m_param=10.654;
	double p0=1.466;

	//smearing params
	float smear_factor_d = 0.02;//0.02;// test trying to include the beam momentum resolution.. will set to zero for now
	float smear_factor_c = 0.028;// first parameter in test beam parameterization?
	float smear_factor_b = 0.155;

    TFile *output = new TFile(pioncode/rootfiles/simplephotoncheck.root, "recreate");
	TTree *tree = new TTree("tree", "tree");


    TH1* h1 = new TH1D("h1", "gamma Pt, power weight",n_bins, PT_Min, PT_Max);
    TH1* h2 = new TH1D("h2", "gamma Pt, exp weight",n_bins, PT_Min, PT_Max);
    TH1* h3 = new TH1D("h3", "gamma Pt, w.s.h.p. weight",n_bins, PT_Min, PT_Max);


    std::random_device rd; // generate a random number to seed random generation
	std::random_device rdgamma; // generate a random number to seed random generation of daughter gamma for smearing
	//std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd
	std::mt19937_64 gen(rd()); // mersenne_twister_engine 64 bit seeded with rd
	std::mt19937_64 gen_gamma(rdgamma()); // mersenne_twister_engine 64 bit seeded  
    std::normal_distribution<double> gammadis(0.0,1.0);// generate normal distribution for gamma smearing, mean zero, variance 1
	std::uniform_real_distribution<> pdis(PT_ratio,1.0);

	for(int i=0; i<NPhotons; i++){

	float gamma1_4mom[4];	
	float gamma2_4mom[4];

    scale_factor1=sqrt(pow(smear_factor_b,2)/gamma1_4mom[1]+pow(smear_factor_c,2)+pow(smear_factor_d,2));
	scale_factor2=sqrt(pow(smear_factor_b,2)/gamma2_4mom[1]+pow(smear_factor_c,2)+pow(smear_factor_d,2));
    smear_factor1=scale_factor1*gammadis(gen_gamma)+1;
	smear_factor2=scale_factor2*gammadis(gen_gamma)+1;

	
	float gamma1_4mom_smeared[4]=smear_factor1*gamma1_4mom;
	float gamma2_4mom_smeared[4]=smear_factor1*gamma1_4mom;



	h1->Fill(pT1,weight_function1);
	h1->Fill(pT2,weight_function1);

	h2->Fill(pT1,weight_function2);
	h2->Fill(pT2,weight_function2);

	h3->Fill(pT1,weight_function3);
	h3->Fill(pT2,weight_function3);
	}
double_t realtime =timer.RealTime();
double_t cputime = timer.CpuTime();

printf("real time =%f\n",realtime);
printf("CPU time=%f\n", cputime);
}

