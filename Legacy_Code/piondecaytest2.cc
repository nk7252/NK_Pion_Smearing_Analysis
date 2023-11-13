// Headers and Namespaces.  
#define _USE_MATH_DEFINES
#include <iostream>
#include <chrono>
//#include <numbers> //std::numbers
#include <cmath>
#include <random>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
     
#include "Pythia8/Pythia.h" // Include Pythia headers.  
using namespace Pythia8;    // Let Pythia8:: be implicit.  
 
int main() {   // Begin main program.  	
	//auto start = std::chrono::steady_clock::now();
	//std::cout << "Start Time (sec) = " << start.count() << std::endl;
	int Nevents=3*30000;
	//int n_bins=1+ceil(log2(Nevents));
	int PT_Max=65;//65
	int PT_Min=0;// cross check to elimate power law problems
	double PT_ratio=PT_Min/PT_Max;
	int n_bins=10*PT_Max;
	//--------------------Alternative paramaterization, woods saxon
	double t=4.5;
	double w=0.114;
	double A=229.6;
	double B=14.43;
	double n=8.1028;
	double m_param=10.654;
	double p0=1.466;
        //--------------------preliminaries to read from root
	TFile *output = new TFile("Pi0FastMC.root", "recreate");
	TTree *tree = new TTree("tree", "tree");
	//TH1* h1 = new TH1F("h1", "pi0 E",128, 0, PT_Max);
	TH1* h2 = new TH1F("h2", "gamma Pt",n_bins, PT_Min, PT_Max);// will be weighted
	TH1* h3 = new TH1F("h3", "pi0 Pt",n_bins, PT_Min, PT_Max);
	TH1* h4 = new TH1F("h4", "pi0 Pt, unweighted",n_bins, PT_Min, PT_Max);
	TH1* h5 = new TH1F("h5", "gamma Pt, unweighted",n_bins, PT_Min, PT_Max);
	TH1* h6 = new TH1F("h6", "gamma inv mass of pair, unweighted",100, 0, 1);
	TH1* h7 = new TH1F("h7", "ratio of gamma/pi0 pt",n_bins, PT_Min, PT_Max);
	TH1* h8 = new TH1F("h8", "ratio of pi0 pt/gamma",n_bins, PT_Min, PT_Max);
	//--------------------set up random number generation
	std::random_device rd; // generate a random number to seed random generation
	//std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd
	std::mt19937_64 gen(rd()); // mersenne_twister_engine 64 bit seeded with rd
	//std::ranlux48 gen(rd()); // ranlux48 seeded with rd
	//--------------------generate random: momentum
	//std::uniform_real_distribution<> pdis(0.0,PT_Max);
	//std::uniform_real_distribution<> pdis(0.0,1.0);// alternative scheme
	std::uniform_real_distribution<> pdis(PT_ratio,1.0);// alternative scheme with min PT to avoid power law complications. 
	//--------------------generate random: angle
	//double tpi=2*std::numbers::pi;
	std::uniform_real_distribution<> adis(0.0,2*M_PI);
	int id, size, no;
	std::vector<int> vec_id;
	std::vector<double> vec_E;
	std::vector<double> vec_Pi0_Pt;
	std::vector<double> vec_P;
	std::vector<double> vec_inv_mass;
	Vec4 gamma_lorentz[3];
	double m, E, px, py ,pz, pi0_px, pi0_py, pi0_E;
        double P0rest= 0.0;
        double pi0_pz=0.0;
        double Pi0_M = 0.1349768;//135 MeV
        double inv_mass, E_sum, P_sum, E_sum_sq, P_sum_sq, weight_function;
	std::cout << Pi0_M <<std::endl;
 	//int idnone=0;
 	//int idg2=0;
 	tree->Branch("id",&id, "id/I");
 	//tree->Branch("size",&size, "size/I");
 	//tree->Branch("no",&no, "no/I");
 	tree->Branch("m",&m, "m/D");
 	tree->Branch("E",&E, "E/D"); // can reconstruct in root easilly
 	//tree->Branch("Pt",&Pt, "Pt/D"); // can reconstruct in root easilly
 	tree->Branch("px",&px, "px/D");
 	tree->Branch("py",&py, "py/D");
 	tree->Branch("pz",&pz, "pz/D");
 	
       // Set up generation.  event
       Pythia pythia;            // Declare Pythia object  
       //pythia.readString("PromptPhoton:all = on"); // Switch on process. 
       //pythia.readString("ParticleDecays:allowPhotonRadiation = on");
       pythia.readString("111:all = pi0 -> gamma gamma"); // Switch on pi0 decay
       pythia.readstring("111:onMode = off"); 
       //pythia.readString("Beams:eCM = 14.e3"); // 14 TeV CM energy.  
       pythia.init(); // Initialize; incoming pp beams is default.  
       //pythia.particleData.mayDecay(111,true);
       //pythia.particleData.mayDecay(22,false);
       //pythia.particleData.mayDecay(11,false);
       //pythia.particleData.mayDecay(-11,false);
       // Generate event(s).  
       //pythia.next(); // Generate an(other) event. Fill event record.  
       pythia.event.clear();
       bool readytoappend = true;
       
       for(int i=0; i<Nevents; i++){
       //std::cout << "I reached here" <<std::endl;// debug line
       			double azimuthal_ang=adis(gen);// generate a random angle from 0 to 2pi
	       		//double Pt=pdis(gen); // generate a random momentum for pi0 from 0 to 64 GeV/c
	       		double Pt=PT_Max*sqrt(pdis(gen));
	       		//double Pt=PT_Max*sqrt(pdis(gen)); // generate a random momentum for pi0 from 0 to 64 GeV/c. use other method. pdis 0 to 1 
	       		vec_Pi0_Pt.push_back(Pt);
	       		weight_function=(1/(1+exp((vec_Pi0_Pt.back()-t)/w)))*A/pow(1+vec_Pi0_Pt.back()/p0,m_param)+(1-(1/(1+exp((vec_Pi0_Pt.back()-t)/w))))*B/(pow(vec_Pi0_Pt.back(),n));//new weight method, need to multiply by Pt too
	       		h3->Fill(vec_Pi0_Pt.back(),vec_Pi0_Pt.back()*weight_function);//fill pi0 pt, weighted
        		h4->Fill(vec_Pi0_Pt.back());//fill pi0 pt, unweighted	
	       		//double Pt=10; // FIxed value Pt. in GeV
	      		pi0_px=Pt*cos(azimuthal_ang);
	       		pi0_py=Pt*sin(azimuthal_ang);
	       		pi0_E=sqrt(Pi0_M*Pi0_M + Pt*Pt + pi0_pz*pi0_pz);
	       		pythia.event.append(111,23, 0, 0 ,pi0_px,pi0_py,pi0_pz, pi0_E, Pi0_M);
	       		
  
        	id = pythia.event[i].id();
        	m  = pythia.event[i].m();
        	px = pythia.event[i].px();
        	py = pythia.event[i].py();
        	pz = pythia.event[i].pz();
        	E = pythia.event[i].e();
        	//if(id==111){
       		//pythia.moreDecays();
        	//}
       		vec_id.push_back(id);
       		vec_E.push_back(E);
       		vec_P.push_back(sqrt(px*px+py*py+pz*pz));
			}
			pythia.forceHadronLevel();
			pythia.event.resetDecays();
			pythia.event.decay();
			for(int i=0;i<pythia.event.size();i++){
				if(pythia.event[i].id() == 111){
				int daughter1 = pythia.event[i].daughter1();
				int daughter2 = pythia.event[i].daughter2();
					for(int j = daughter1; j<= daughter2; j++){
						if(pythia.event[j].id() == 22){	
							gamma_lorentz=
						}
					}
				}
       		        if(vec_id[vec_id.size()-1]==22 && vec_id[vec_id.size()-2]==22){//find inv mass
       		        	gamma_lorentz[1]=pythia.event[vec_id.size()-1].p();
       		        	gamma_lorentz[2]=pythia.event[vec_id.size()-2].p();
       		        	gamma_lorentz[3]=gamma_lorentz[1]+gamma_lorentz[2];
       		        	inv_mass=gamma_lorentz[3].mCalc();
       		        	//std::cout << "inv mass" << " " << gamma_lorentz[3].mCalc() <<std::endl;
			       
       		        	//inv_mass=gamma_lorentz[3].mCalc();
       		        	vec_inv_mass.push_back(inv_mass);
       		        	std::cout << "inv mass" << " " <<vec_inv_mass.back()<<std::endl;
       		        	//std::cout << "last two ids" <<" "<< vec_id[vec_id.size()-1] << " " << vec_id[vec_id.size()-2] << " " << "last two E" <<" "<< vec_E[vec_E.size()-1] << " " << vec_E[vec_E.size()-2] << " " << "last two P" <<" "<< vec_P[vec_P.size()-1] << " " << vec_P[vec_P.size()-2]<<std::endl;
       		        	h6->Fill(vec_inv_mass.back());
       		        	readytoappend = true;
       		        	
       		        	
       		        }
       		        
       		        ///*
       		if(abs(vec_id[vec_id.size()-1])==11 && abs(vec_id[vec_id.size()-2])==11){
       			readytoappend = true;
       		}//*/
       		        
       		         ///*
       		      
			if(vec_id[vec_id.size()-1]==22 && vec_id[vec_id.size()-2]==22){
        		//if(id==22){// need to do it like the above loop. add this and previous
        			gamma_lorentz[1]=pythia.event[vec_id.size()-1].p();
        			gamma_lorentz[2]=pythia.event[vec_id.size()-2].p();
        	        	double gPt1=gamma_lorentz[1].pT();
        	        	double gPt2=gamma_lorentz[2].pT();
        	        	//weight_function=(1/(1+exp((vec_Pi0_Pt.back()-t)/w)))*A/pow(1+vec_Pi0_Pt.back()/p0,m_param)+(1-(1/(1+exp((vec_Pi0_Pt.back()-t)/w))))*B/(pow(vec_Pi0_Pt.back(),n));
        	        	
        			//h2->Fill(gPt1,pow(vec_Pi0_Pt.back(),-7.1));//fill gamma pt, weighted
        			h2->Fill(gPt1,vec_Pi0_Pt.back()*weight_function);
        			h5->Fill(gPt1);//fill gamma pt, unweighted
        			//h2->Fill(gPt2,pow(vec_Pi0_Pt.back(),-7.1));//fill gamma pt, weighted
        			h2->Fill(gPt2,vec_Pi0_Pt.back()*weight_function);
        			h5->Fill(gPt2);//fill gamma pt, unweighted		
        			//std::cout << "gamma" <<std::endl;// test if the loop is working
       			}
       			//*/
       			
       			
       			
       			//if(idnone==22){
       			//idg2=1;
       			//if(idg2==1){
       			//idg2=0;
       			//}
        	           //double gPt=sqrt(px*px+py*py);
        		   //h2->Fill(gPt,pow(Pt,-8.1));//fill gamma pt, weighted
        		   //h5->Fill(gPt);//fill gamma pt, unweighted
        		   //std::cout << "double gamma" <<std::endl;// test if the loop is working
       			//}


       		tree->Fill();
       		std::cout << i <<" "<<id<< " " << m << " " << " E " << " " << E <<" " << " P " << " " << px<<" "<<py<<" "<<pz<<std::endl;


       		//idnone=id;
       		//std::cout << "random" <<" "<<pdis(gen)<< " " <<"random"<<" "<<adis(gen)<<std::endl;// check random number generation
       		}
       		h7->Divide(h2,h3);
       		h8->Divide(h3,h2);
       		output->Write();
       		output->Close();
	//auto end =std::chrono::steady_clock::now();
	//auto elapsed = end - start;
	//std::cout << "elapsed (sec) = " << elapsed.count() << std::endl;    
return 0;  
}                // End main program with error-free return. 

