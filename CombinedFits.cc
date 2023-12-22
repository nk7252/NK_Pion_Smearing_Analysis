#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <TFile.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TGraphMultiErrors.h>

class filename_object {//object to hold my file names and related strings. things like settings for the files are contained here.
    public:
    std::vector<std::string> fileNames;
    std::vector<std::string> legendnames;
    std::vector<std::string> weightnames;
    std::vector<std::string> histlist;
    string filenamemod;
    string canvasnamemod;
    std::vector<float> plotxlims;
    std::vector<float> plotylims;
    int pTcutoff;
    std::vector<float> sqrtEsmearing;
    int binres;


};

float extractNumber(const std::string& filePath);
filename_object choosecomparisontype(int choosetype);
void OverlayMeans(filename_object filenameobj);
void OverlaySigmaOverMean(filename_object filenameobj);
void ClusterOverlayTestFunc(filename_object filenameobj, const char* histname);
TH1D* getYProjectionof2DHist(const char* fileName, const char* histName, int firstxbin, int lastxbin);
void plotOverlayedHistograms(filename_object filenameobj, const char* histName);
void FitYProjectionsAndGraph(TCanvas *canvas, TH2F *hist2D, const char *pdfname, double binres, int weightmethod);
void SliceAndFit(filename_object filenameobj);
void OverlayMeansAIO(filename_object filenameobj, const char* histName);
void OverlaySigmaOverMeanAIO(filename_object filenameobj, const char* histName);

void CombinedFits() {
    int fileset = 6 ;
    filename_object choosenfilenameobj = choosecomparisontype(fileset);
    // 0=weight type, 1=ac on/off, 2=co on/off, 4 ac&co on/off
    //5=weight type, new files+ hagedorn
    // want to do sets 1 and 5 today.

    if (fileset<5){// different files for different cuts/changes
    OverlayMeans(choosenfilenameobj);
    OverlaySigmaOverMean(choosenfilenameobj);
    plotOverlayedHistograms(choosenfilenameobj, "h12");//h12 is smeared pion pT, Weighted. h3 is unsmeared pion pT, weighted
    //SliceAndFit(choosenfilenameobj);

    /*
    if (choosenfilenameobj.fileNames.size()==2){// for subtraction of inv mass profile
        ClusterOverlayTestFunc(choosenfilenameobj, "h18");//const char* histname
    }
    //*/
    }
    // hist list: h18- smeared pTs vs smeared inv mass, h27- cluster, h28 cluster+asymm, h29 asymm
    if (fileset>5){// all in one file
        //tbd
        OverlayMeansAIO(choosenfilenameobj, "h27");
        OverlaySigmaOverMeanAIO(choosenfilenameobj, "h27");
        //plotOverlayedHistogramsAIO(choosenfilenameobj, "h12");//h12 is smeared pion pT, Weighted. h3 is unsmeared pion pT, weighted
        //SliceAndFit(choosenfilenameobj);
    }
}    

float extractNumber(const std::string& filePath) {
    // Find the first underscore in the filepath
    size_t firstUnderscorePos = filePath.find('_');

    // Find the second underscore in the filepath
    size_t secondUnderscorePos = filePath.find('_', firstUnderscorePos + 1);
    float testnumber;
    // Check if both underscores are found and there is a number between them
    if (firstUnderscorePos != std::string::npos && secondUnderscorePos != std::string::npos && secondUnderscorePos > firstUnderscorePos + 1) {
        //std::cout<< "I reached here" << std::endl;
        // Extract the substring between the two underscores
        std::string numberStr = filePath.substr(firstUnderscorePos + 1, secondUnderscorePos - firstUnderscorePos - 1);
        //std::cout<< "substring" <<" "<< numberStr << std::endl;
        testnumber =std::stof(numberStr);
        //std::cout<< "float" <<" "<< testnumber << std::endl;
    } 
    return testnumber;
} 

filename_object choosecomparisontype(int choosetype){
    filename_object filename_object1;// 0=weight type, 1=ac on/off, 2=co on/off, 3=ac&co on/off
    if(choosetype==0){
        //filename_object weightfilenameobj;
        filename_object1.fileNames={"pioncode/rootfiles/Pi0FastMC_0.155000_EXP.root", "pioncode/rootfiles/Pi0FastMC_0.155000_POWER.root", "pioncode/rootfiles/Pi0FastMC_0.155000_WSHP.root"};
        filename_object1.legendnames={"EXP","POWER","WSHP"};
        filename_object1.filenamemod="weightmethod";
        filename_object1.canvasnamemod=" for various weighting methods";  
        filename_object1.plotxlims={0.9,6.4};//min, max
        filename_object1.plotylims={0.13,0.17,0.08,0.25,0.0, 2.0}; //mean_min, mean_max,sm_min, sm_max, min h12, max h12
        filename_object1.pTcutoff=6;
        filename_object1.binres=2;
    }
    else if(choosetype==1){
        filename_object1.fileNames={"pioncode/rootfiles/Pi0FastMC_0.155000_WSHP_ac0_co0.root", "pioncode/rootfiles/Pi0FastMC_0.155000_WSHP_ac1_co0.root","pioncode/rootfiles/Pi0FastMC_0.155000_EXP_ac0_co0.root", "pioncode/rootfiles/Pi0FastMC_0.155000_EXP_ac1_co0.root","pioncode/rootfiles/Pi0FastMC_0.155000_POWER_ac0_co0.root", "pioncode/rootfiles/Pi0FastMC_0.155000_POWER_ac1_co0.root","pioncode/rootfiles/Pi0FastMC_0.155000_HAGEDORN_ac0_co0.root", "pioncode/rootfiles/Pi0FastMC_0.155000_HAGEDORN_ac1_co0.root"};
        filename_object1.legendnames={"Asym. cut off, WSHP","Asym. cut on, WSHP","Asym. cut off, EXP","Asym. cut on, EXP","Asym. cut off, POWER","Asym. cut on, POWER","Asym. cut off, HAGEDORN","Asym. cut on, HAGEDORN"};
        filename_object1.weightnames={"WSHP","WSHP","EXP","EXP","POWER","POWER","HAGEDORN","HAGEDORN"};
        filename_object1.filenamemod="AsymCutTest";
        filename_object1.canvasnamemod=", Asymmetry Cut: on vs off";
        filename_object1.plotxlims={0.1,16.4};//min, max
        filename_object1.plotylims={0.12,0.17,0.0,0.4,0.0, 2.0}; //mean_min, mean_max,sm_min, sm_max, min h12, max h12
        filename_object1.pTcutoff=16;
        filename_object1.binres=2;
    }
    else if(choosetype==2){
        filename_object1.fileNames={"pioncode/rootfiles/Pi0FastMC_0.155000_WSHP_ac0_co0.root", "pioncode/rootfiles/Pi0FastMC_0.155000_WSHP_ac0_co1.root"};
        filename_object1.legendnames={"Cluster Overlay off","Cluster Overlay on"};
        filename_object1.filenamemod="ClusterOverlayTest";
        filename_object1.canvasnamemod=", Cluster Overlap: on vs off";
        filename_object1.plotxlims={0.1,16.4};//min, max
        filename_object1.plotylims={0.135,0.145,0.05,0.25}; //mean_min, mean_max,sm_min,sm_max 
        filename_object1.pTcutoff=16; 
        filename_object1.binres=2;
    }
    else if(choosetype==3){
        filename_object1.fileNames={"pioncode/rootfiles/Pi0FastMC_0.155000_WSHP_ac0_co0.root", "pioncode/rootfiles/Pi0FastMC_0.155000_WSHP_ac1_co0.root", "pioncode/rootfiles/Pi0FastMC_0.155000_WSHP_ac0_co1.root", "pioncode/rootfiles/Pi0FastMC_0.155000_WSHP_ac1_co1.root"};
        filename_object1.legendnames={"AC,CO:00","AC,CO:10","AC,CO:01","AC,CO:11"};
        filename_object1.filenamemod="AsymClusterTest";
        filename_object1.canvasnamemod=", Cluster Overlap & Asymmmetry Cut: on(1) vs off(0)";
        filename_object1.plotxlims={0.1,16.4};//min, max
        filename_object1.plotylims={0.125,0.15,0.01,0.4}; //mean_min, mean_max,sm_min,sm_max  
        filename_object1.pTcutoff=16; 
        filename_object1.binres=2;
    }
    else if(choosetype==4){//co for 3% smearing
        filename_object1.fileNames={"pioncode/rootfiles/Pi0FastMC_0.030000_WSHP_ac0_co0.root", "pioncode/rootfiles/Pi0FastMC_0.030000_WSHP_ac0_co1.root"};
        filename_object1.legendnames={"Cluster Overlay off","Cluster Overlay on"};
        filename_object1.filenamemod="ClusterOverlayTest";
        filename_object1.canvasnamemod=", Cluster Overlap: on vs off";
        filename_object1.plotxlims={0.1,16.4};//min, max
        filename_object1.plotylims={0.1,0.137,0.0,0.085}; //mean_min, mean_max,sm_min,sm_max 
        filename_object1.pTcutoff=16; 
        filename_object1.binres=2;
    }
    else if(choosetype==5){
        //filename_object weightfilenameobj;
        filename_object1.fileNames={"pioncode/rootfiles/Pi0FastMC_0.155000_EXP_ac0_co0.root", "pioncode/rootfiles/Pi0FastMC_0.155000_POWER_ac0_co0.root", "pioncode/rootfiles/Pi0FastMC_0.155000_WSHP_ac0_co0.root","pioncode/rootfiles/Pi0FastMC_0.155000_HAGEDORN_ac0_co0.root"};
        filename_object1.legendnames={"EXP","POWER","WSHP","HAGEDORN"};
        filename_object1.weightnames={"EXP","POWER","WSHP","HAGEDORN"};
        filename_object1.filenamemod="weightmethod";
        filename_object1.canvasnamemod=" for various weighting methods";  
        filename_object1.plotxlims={0.9,6.4};//min, max
        filename_object1.plotylims={0.13,0.17,0.08,0.25,0.0, 2.0}; //mean_min, mean_max,sm_min, sm_max, min h12, max h12
        filename_object1.pTcutoff=6;
        filename_object1.binres=2;
    }// AOI sets following
    else if(choosetype==6){
        //filename_object weightfilenameobj;
        filename_object1.fileNames={"pioncode/rootfiles/Pi0FastMC_0.155000_EXP.root", "pioncode/rootfiles/Pi0FastMC_0.155000_POWER.root", "pioncode/rootfiles/Pi0FastMC_0.155000_WSHP.root","pioncode/rootfiles/Pi0FastMC_0.155000_HAGEDORN.root"};
        filename_object1.legendnames={"EXP","POWER","WSHP","HAGEDORN"};
        filename_object1.weightnames={"EXP","POWER","WSHP","HAGEDORN"};
        filename_object1.filenamemod="weightmethod_co1_ac1";
        filename_object1.canvasnamemod=" for various weighting methods. asymm+clustering";  
        filename_object1.plotxlims={0.9,6.4};//min, max
        filename_object1.plotylims={0.13,0.17,0.08,0.25,0.0, 2.0}; //mean_min, mean_max,sm_min, sm_max, min h12, max h12
        filename_object1.pTcutoff=6;
        filename_object1.binres=2;
    }

    for(size_t i=0; i < filename_object1.fileNames.size(); i++){
        filename_object1.sqrtEsmearing.push_back(extractNumber(filename_object1.fileNames[i]));
    }

return filename_object1;
}

//-----------------------------------------split file functions
void OverlayMeans(filename_object filenameobj) {
//void OverlayMeans(const std::vector<std::string> & fileNames, string filenamemod) {
    // Create a TCanvas
    TCanvas* canvas1 = new TCanvas("canvas1", "Overlay Means", 800, 600);
    
    //canvas1->SetGrid();
    //gStyle->SetOptTitle ( 0 );
    //gStyle->SetOptStat ( 0 );
    //gStyle->SetOptFit ( 1 1 1 1 );
    //gStyle->SetStatX ( . 8 9 ); 
    //gStyle->SetStatY ( . 8 9 );
    //gStyle->SetStatBorderSize ( 0 );
    float errparam=filenameobj.sqrtEsmearing[0];
    double binres=2;//number of divisions per GeV
    // Create a legend
    TLegend* legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);//0.7, 0.4, 0.9, 0.6
    TMultiGraph *MultiGraphs = new TMultiGraph();//h18->GetNbinsX()
    // Loop over each file
    for (size_t i = 0; i < filenameobj.fileNames.size(); ++i) {
        // Open the root file
        TFile* pionfile = new TFile(filenameobj.fileNames[i].c_str(), "READ");

        // Check if the file is open
        if (!pionfile || pionfile->IsZombie()) {
            std::cerr << "Error: Unable to open file " << filenameobj.fileNames[i] << std::endl;
            continue;
        }

        // Retrieve the 2D histogram
        //TH2F *h18 =(TH2F *)pionfile->Get("h18");
        TH2F* h18 = nullptr;
        pionfile->GetObject("h18", h18);

        // Check if the histogram exists
        if (!h18) {
            std::cerr << "Error: Histogram 'myHist2D' not found in file " << filenameobj.fileNames[i] << std::endl;
            pionfile->Close();
            continue;
        }

        // Create a histogram for means
        TGraphErrors *meanGraph = new TGraphErrors(h18->GetNbinsX());//

        // Loop over each bin in the X direction
        for (int binX = 1; binX <= h18->GetNbinsX(); ++binX) {
            //std::cout << "Nbin" << " " << binX << std::endl;
            // Project along Y for each sbinX
            TH1D* yProjection = h18->ProjectionY(Form("YProjection_%zu_%d", i, binX), binX, binX, "");

            // Fit the Y projection with a Gaussian
            yProjection->Fit("gaus", "Q");

            // Access the fit parameters
            TF1* fitFunc = yProjection->GetFunction("gaus");

            // Check if the fit function is valid
            if (fitFunc) {
                // Fill the mean histogram with the mean value
                //meanHistogram->SetBinContent(binX, fitFunc->GetParameter(1));
                if (filenameobj.weightnames[i]=="EXP"){//exp &&  binX <= filenameobj.pTcutoff*binres
                    meanGraph->SetPoint(binX, binX/binres,fitFunc->GetParameter(1));
                    meanGraph->SetPointError(binX, 0,fitFunc->GetParError(1));
                    if(binX/binres==3){
                        std::cout << binX <<" "<<fitFunc->GetParameter(1) << std::endl; // debug line
                    }
                }
                else if (filenameobj.weightnames[i]=="POWER"){//power &&  3*binres < binX
                    meanGraph->SetPoint(binX, binX/binres,fitFunc->GetParameter(1));
                    meanGraph->SetPointError(binX, 0,fitFunc->GetParError(1));
                    if(binX/binres==3){
                        std::cout << binX <<" "<<fitFunc->GetParameter(1) << std::endl; // debug line
                    }
                    //std::cout << binX <<" "<<fitFunc->GetParameter(1) << std::endl; // debug line
                }
                else if (filenameobj.weightnames[i]=="WSHP"){//woods saxon
                    meanGraph->SetPoint(binX, binX/binres,fitFunc->GetParameter(1));
                    meanGraph->SetPointError(binX, 0,fitFunc->GetParError(1));
                }
                else if (filenameobj.weightnames[i]=="HAGEDORN"){//HAGEDORN
                    meanGraph->SetPoint(binX, binX/binres,fitFunc->GetParameter(1));
                    meanGraph->SetPointError(binX, 0,fitFunc->GetParError(1));
                }
                
                //meanHistogram->SetBinError(binX, fitFunc->GetParError(1));
                //std::cout << "bin number" << " " << binX << " " << "Mean" << " " <<  fitFunc->GetParameter(1) << " " << "Mean error" << " " << fitFunc->GetParError(1) << std::endl;
            }
                // Add an entry to the legend
                //legend1->AddEntry(yProjection, Form("Version %zu, BinX %d", i, binX), "L");
            
            //std::cout << "I reached here, pre delete proj" << std::endl; // debug line
            // Clean up Y projection
            delete yProjection;
        }
        // Set different line colors for each version
        int MarkerStyle = i + 24; // 
        int MarkerColor = i + 1;
        //meanGraph->SetLineColor(lineColor);
        meanGraph->SetMarkerStyle(MarkerStyle);
        meanGraph->SetMarkerColor(MarkerColor);
        

       
        // Overlay the mean histogram on the same canvas
        if (i == 0) {
            //meanGraph->Draw("AP"); // Draw histogram for the first version
            //canvas1->Print("OverlayMeanHistograms.pdf");
            MultiGraphs->Add(meanGraph,"PE");
            std::cout << "draw for the first file" << std::endl; // debug line
        } else {
            MultiGraphs->Add(meanGraph,"PE");
            //meanGraph->Draw("P SAME"); // Draw subsequent histograms on the same canvas
            std::cout << "draw for subsequent" << std::endl; // debug line
        }   
        
        MultiGraphs->SetTitle(Form("Smeared Inv. Mass%s;pT (GeV);Inv. Mass (GeV)",filenameobj.canvasnamemod.c_str()));
        MultiGraphs->Draw("APE");

        // Add an entry to the legend
        //std::vector<std::string> legendstring = {"EXP","POWER","WSHP"};
        legend1->AddEntry(meanGraph, filenameobj.legendnames[i].c_str(), "P");

        // Close the file
        pionfile->Close();
        delete pionfile;

        //std::cout << "I reached here, close+delete file" << std::endl; // debug line
    }

    //std::cout << "I reached here, done with all files" << std::endl; // debug line
    // Draw the legend
    legend1->Draw();

    // Show the canvas
    MultiGraphs->GetXaxis()->SetLimits(filenameobj.plotxlims[0],filenameobj.plotxlims[1]);
    MultiGraphs->SetMinimum(filenameobj.plotylims[0]);
    MultiGraphs->SetMaximum(filenameobj.plotylims[1]);
    canvas1->SetMargin(0.2,0.1,0.1,0.1);
    gPad->Modified();
    gPad->Update();
    //canvas1->Modified();
    canvas1->Update();
    canvas1->Print(Form("pioncode/canvas_pdf/%s_%f_OverlayMeanHistograms.pdf",filenameobj.filenamemod.c_str(),errparam));

    // Clean up
    delete canvas1;
    delete legend1;
}

void OverlaySigmaOverMean(filename_object filenameobj) {
//void OverlaySigmaOverMean(const std::vector<std::string>& fileNames, string filenamemod) {
    // Create a TCanvas
    TCanvas* canvas1 = new TCanvas("canvas1", "Overlay Means", 800, 600);
    //canvas1->SetGrid();
    //gStyle->SetOptTitle ( 0 );
    //gStyle->SetOptStat ( 0 );
    //gStyle->SetOptFit ( 1 1 1 1 );
    //gStyle->SetStatX ( . 8 9 ); 
    //gStyle->SetStatY ( . 8 9 );
    //gStyle->SetStatBorderSize ( 0 );

    float errparam=filenameobj.sqrtEsmearing[0];
    double binres=2;//number of divisions per GeV
    // Create a legend
    TLegend* legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    TMultiGraph *MultiGraphs = new TMultiGraph();//h18->GetNbinsX()
    // Loop over each file
    for (size_t i = 0; i < filenameobj.fileNames.size(); ++i) {
        // Open the root file
        TFile* pionfile = new TFile(filenameobj.fileNames[i].c_str(), "READ");

        // Check if the file is open
        if (!pionfile || pionfile->IsZombie()) {
            std::cerr << "Error: Unable to open file " << filenameobj.fileNames[i] << std::endl;
            continue;
        }

        // Retrieve the 2D histogram
        //TH2F *h18 =(TH2F *)pionfile->Get("h18");
        TH2F* h18 = nullptr;
        pionfile->GetObject("h18", h18);

        // Check if the histogram exists
        if (!h18) {
            std::cerr << "Error: Histogram 'myHist2D' not found in file " << filenameobj.fileNames[i] << std::endl;
            pionfile->Close();
            continue;
        }

        // Create a histogram for means
        TGraphErrors *meanGraph = new TGraphErrors(h18->GetNbinsX());
        //TH1F* meanHistogram = new TH1F(Form("MeanHistogram_%zu", i), Form("Version %zu", i), h18->GetNbinsX(), 0.5, h18->GetNbinsX() + 0.5);



        // Loop over each bin in the X direction
        for (int binX = 1; binX <= h18->GetNbinsX(); ++binX) {
            // Project along Y for each binX
            TH1D* yProjection = h18->ProjectionY(Form("YProjection_%zu_%d", i, binX), binX, binX, "");

            // Fit the Y projection with a Gaussian
            yProjection->Fit("gaus", "Q");

            // Access the fit parameters
            TF1* fitFunc = yProjection->GetFunction("gaus");

            // Check if the fit function is valid
            if (fitFunc) {
                // Fill the mean histogram with the mean value
                //meanHistogram->SetBinContent(binX, fitFunc->GetParameter(1));
                if (filenameobj.weightnames[i]=="EXP"){//exp && binX <= filenameobj.pTcutoff*binres
                    double meanoversigma =fitFunc->GetParameter(2)/fitFunc->GetParameter(1);
                    double meanoversigmaerr = meanoversigma*(fitFunc->GetParError(2)/fitFunc->GetParameter(2)+fitFunc->GetParError(1)/fitFunc->GetParameter(1));//m/s*(serr/s+merr/m)

                    meanGraph->SetPoint(binX, binX/binres,meanoversigma);
                    meanGraph->SetPointError(binX, 0,meanoversigmaerr);
                }
                else if (filenameobj.weightnames[i]=="POWER"){//power 
                    double meanoversigma =fitFunc->GetParameter(2)/fitFunc->GetParameter(1);
                    double meanoversigmaerr=meanoversigma*(fitFunc->GetParError(2)/fitFunc->GetParameter(2)+fitFunc->GetParError(1)/fitFunc->GetParameter(1));//m/s*(serr/s+merr/m)

                    meanGraph->SetPoint(binX, binX/binres,meanoversigma);
                    meanGraph->SetPointError(binX, 0,meanoversigmaerr);
                    std::cout << binX <<" "<<meanoversigma << std::endl; // debug line

                }
                else if (filenameobj.weightnames[i]=="WSHP") {//woods saxon
                    double meanoversigma =fitFunc->GetParameter(2)/fitFunc->GetParameter(1);
                    double meanoversigmaerr=meanoversigma*(fitFunc->GetParError(2)/fitFunc->GetParameter(2)+fitFunc->GetParError(1)/fitFunc->GetParameter(1));//(m/s)_err=m/s*(serr/s+merr/m)

                    meanGraph->SetPoint(binX, binX/binres,meanoversigma);
                    meanGraph->SetPointError(binX, 0,meanoversigmaerr);
                }
                else if (filenameobj.weightnames[i]=="HAGEDORN") {//HAGEDORN
                    double meanoversigma =fitFunc->GetParameter(2)/fitFunc->GetParameter(1);
                    double meanoversigmaerr=meanoversigma*(fitFunc->GetParError(2)/fitFunc->GetParameter(2)+fitFunc->GetParError(1)/fitFunc->GetParameter(1));//(m/s)_err=m/s*(serr/s+merr/m)

                    meanGraph->SetPoint(binX, binX/binres,meanoversigma);
                    meanGraph->SetPointError(binX, 0,meanoversigmaerr);
                }

                
                //meanHistogram->SetBinError(binX, fitFunc->GetParError(1));
                //std::cout << "bin number" << " " << binX << " " << "Mean" << " " <<  fitFunc->GetParameter(1) << " " << "Mean error" << " " << fitFunc->GetParError(1) << std::endl;
            }
                // Add an entry to the legend
                //legend1->AddEntry(yProjection, Form("Version %zu, BinX %d", i, binX), "L");
            
            //std::cout << "I reached here, pre delete proj" << std::endl; // debug line
            // Clean up Y projection
            delete yProjection;
        }
        // Set different line colors for each version
        int MarkerStyle = i + 24; // 
        int MarkerColor = i + 1;
        //meanGraph->SetLineColor(lineColor);
        meanGraph->SetMarkerStyle(MarkerStyle);
        meanGraph->SetMarkerColor(MarkerColor);

        //std::cout << "I reached here, done with loop over bins" << std::endl; // debug line

        // Overlay the mean histogram on the same canvas
        if (i == 0) {
            //meanGraph->Draw("AP"); // Draw histogram for the first version
            //canvas1->Print("OverlayMeanHistograms.pdf");
            MultiGraphs->Add(meanGraph,"PE");
            std::cout << "draw for the first file" << std::endl; // debug line
        } else {
            MultiGraphs->Add(meanGraph,"PE");
            //meanGraph->Draw("P SAME"); // Draw subsequent histograms on the same canvas
            std::cout << "draw for subsequent" << std::endl; // debug line
        }   

        MultiGraphs->SetTitle(Form(" Sigma/Mean (Smeared Inv. Mass)%s;pT (GeV);sigma/mean)",filenameobj.canvasnamemod.c_str()));
        MultiGraphs->Draw("APE");
        // Add an entry to the legend
        //std::vector<std::string> legendstring = {"EXP","POWER","WSHP"};
        legend1->AddEntry(meanGraph, filenameobj.legendnames[i].c_str(), "P");

        // Close the file
        pionfile->Close();
        delete pionfile;

        std::cout << "I reached here, close+delete file" << std::endl; // debug line
    }

    std::cout << "I reached here, done with all files" << std::endl; // debug line
    // Draw the legend
    legend1->Draw();
    
    MultiGraphs->GetXaxis()->SetLimits(filenameobj.plotxlims[0],filenameobj.plotxlims[1]);
    MultiGraphs->SetMinimum(filenameobj.plotylims[2]);
    MultiGraphs->SetMaximum(filenameobj.plotylims[3]);
    canvas1->SetMargin(0.2,0.1,0.1,0.1);
    gPad->Modified();
    gPad->Update();

   //canvas1->Update();
    canvas1->Modified();
    canvas1->SaveAs(Form("pioncode/canvas_pdf/%s_%f_OverlaySigma_Over_Mean.pdf",filenameobj.filenamemod.c_str(),errparam));//Print-> works too
    
    // Clean up
    delete canvas1;
    delete legend1;
}

void ClusterOverlayTestFunc(filename_object filenameobj, const char* histname){//only works if filenames.size()=2 !!
    TCanvas* canvas1 = new TCanvas("canvas1", "Overlay Means", 800, 600);
    TLegend* legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    
    float errparam=filenameobj.sqrtEsmearing[0];

    ///*
    //need number of lines in hist. temp until I setup the object to hold that info?
    TFile* file = new TFile(filenameobj.fileNames[0].c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file " << filenameobj.fileNames[0].c_str() << std::endl;
        return nullptr;
    }

    // Get the 2D histogram from the file
    TH2F* hist2D = dynamic_cast<TH2F*>(file->Get(histname));
    if (!hist2D) {
        std::cerr << "Error: Could not retrieve 2D histogram " << file->Get(histname) << " from file" << std::endl;
        file->Close();
        return nullptr;
    }
    int NX= hist2D->GetNbinsX();
    file->Close();//*/
    //open the pdf?
    //int NX=128;
    canvas1->Print(Form("pioncode/canvas_pdf/%s_%f_InvMassprojections.pdf[",filenameobj.filenamemod.c_str(), errparam));

    for (int i=1;i<NX+1;i++){
        TH1D* yProjection1 = getYProjectionof2DHist(filenameobj.fileNames[0].c_str(), histname,i,i);
        TH1D* yProjection2 = getYProjectionof2DHist(filenameobj.fileNames[1].c_str(), histname,i,i);
        TH1D *histClone = (TH1D *)yProjection2->Clone("histClone");
        TH1D *ratioClone = (TH1D *)yProjection2->Clone("ratioClone");
        histClone->Add(yProjection1, -1);
        canvas1->Divide(1,3);
        canvas1->cd(1);   
        yProjection1->Draw();
        yProjection1->SetLineColor(kRed);
        yProjection2->Draw("SAME");
        yProjection2->SetLineColor(kBlue);
        yProjection1->SetTitle(Form("Cluster Overlay on vs off. Bin %i;Invariant Mass (GeV); Counts",i));
        gPad->Modified();
        gPad->Update();
        canvas1->cd(2);
        histClone->Draw();
        histClone->SetTitle(Form("Cluster Overlay on - off. Bin %i;Invariant Mass (GeV); Counts",i));
        canvas1->cd(3);
        ratioClone->Divide(yProjection1);
        ratioClone->Draw();
        ratioClone->SetMaximum(5);
        ratioClone->SetTitle(Form("Cluster Overlay on/off. Bin %i;Invariant Mass (GeV); Counts",i));
        canvas1->Modified();
        legend1->AddEntry(yProjection1, filenameobj.legendnames[0].c_str(), "P");
        legend1->AddEntry(yProjection2, filenameobj.legendnames[1].c_str(), "P");
        canvas1->Print(Form("pioncode/canvas_pdf/%s_%f_InvMassprojections.pdf",filenameobj.filenamemod.c_str(), errparam));//Print-> works too
        canvas1->Clear();
        legend1->Clear();
    }
    TH1D* yProjection1 = getYProjectionof2DHist(filenameobj.fileNames[0].c_str(), histname,1,NX);
    TH1D* yProjection2 = getYProjectionof2DHist(filenameobj.fileNames[1].c_str(), histname,1,NX);
    TH1D *histClone = (TH1D *)yProjection2->Clone("histClone");
    TH1D *ratioClone = (TH1D *)yProjection2->Clone("ratioClone");
    histClone->Add(yProjection1, -1);
    canvas1->Divide(1,3);
    canvas1->cd(1);   
    yProjection1->Draw();
    yProjection1->SetLineColor(kRed);
    yProjection2->Draw("SAME");
    yProjection2->SetLineColor(kBlue);
    yProjection1->SetTitle("Cluster Overlay on vs off. All Bins;Invariant Mass (GeV); Counts");
    legend1->AddEntry(yProjection1, filenameobj.legendnames[0].c_str(), "P");
    legend1->AddEntry(yProjection2, filenameobj.legendnames[1].c_str(), "P");
    gPad->Modified();
    gPad->Update();
    canvas1->cd(2);
    histClone->Draw();
    histClone->SetTitle("Cluster Overlay on - off. All Bins;Invariant Mass (GeV); Counts");
    canvas1->cd(3);
    ratioClone->Divide(yProjection1);
    ratioClone->Draw();
    ratioClone->SetMaximum(5);
    ratioClone->SetTitle("Cluster Overlay on/off. All Bins;Invariant Mass (GeV); Counts");
    canvas1->Modified();
    canvas1->Print(Form("pioncode/canvas_pdf/%s_%f_InvMassprojections.pdf",filenameobj.filenamemod.c_str(), errparam));//Print-> works too
    //close the pdf?
    canvas1->Print(Form("pioncode/canvas_pdf/%s_%f_InvMassprojections.pdf]",filenameobj.filenamemod.c_str(), errparam));
    // Clean up
    delete canvas1;
    //delete legend1;
}

TH1D* getYProjectionof2DHist(const char* fileName, const char* histName, int firstxbin, int lastxbin) {
    // Open the root file
    TFile* file = new TFile(fileName, "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file " << fileName << std::endl;
        return nullptr;
    }

    // Get the 2D histogram from the file
    TH2F* hist2D = dynamic_cast<TH2F*>(file->Get(histName));
    if (!hist2D) {
        std::cerr << "Error: Could not retrieve 2D histogram " << histName << " from file" << std::endl;
        file->Close();
        return nullptr;
    }
    //int NX = hist2D->GetNbinsX();
    // Create a 1D histogram for the y projection
    TH1D* hist1D = hist2D->ProjectionY("_py", firstxbin, lastxbin,"");
    //TH1D* hist1D = hist2D->ProjectionY("_py",1, 20,"");
    hist1D->SetDirectory(0);
    // Close the file
    file->Close();

    return hist1D;
}

void plotOverlayedHistograms(filename_object filenameobj, const char* histName) {

    // Create a canvas
    TCanvas* canvas = new TCanvas("canvas", "Overlayed Histograms", 800, 600);
    gStyle->SetOptStat ( 0 );
    // create a legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    // Loop over each file
    for (size_t i = 0; i < filenameobj.fileNames.size(); ++i) {
        // Open the root file
        TFile* pionfile = new TFile(filenameobj.fileNames[i].c_str(), "READ");

        // Check if the file is open
        if (!pionfile || pionfile->IsZombie()) {
            std::cerr << "Error: Unable to open file " << filenameobj.fileNames[i] << std::endl;
            continue;
        }

        // Retrieve the histogram
        TH1F* hist1 = dynamic_cast<TH1F*>(pionfile->Get(histName));

        // Check if histogram was retrieved successfully
        if (!hist1) {
            std::cerr << "Error: Could not retrieve histogram" << " .file" << filenameobj.fileNames[i].c_str() << std::endl;
            pionfile->Close();
            return;
        }
        // Draw The current Histogram
        if(i==0){
            hist1->Draw("HIST");
        }
        else{
            hist1->Draw("SAME HIST");
        }
        //set marker options
        //int MarkerStyle = i + 24; // 
        //int MarkerColor = i + 1;
        //int LineStyle = i ; // 
        int LineColor = i + 1;
        //hist1->SetLineStyle(LineStyle);
        hist1->SetLineColor(LineColor);
        //hist1->GetXaxis()->SetLimits(0,0.1);
        hist1->SetMinimum(1e-6);
        hist1->SetMaximum(1e10); 
        //hist1->SetAxisRange(0., 10000,"Y");
        hist1->SetAxisRange(0., 6.4,"X");
        hist1->GetXaxis()->SetTitle("Pion pT (GeV)");// to properly do this I would need to save one for each histogram somewhere
        hist1->GetYaxis()->SetTitle("Counts");
        // add an entry to the legend  
        legend->AddEntry(hist1, filenameobj.legendnames[i].c_str());//, "P"

        //detach histogram from file
        hist1->SetDirectory(0);

        // Close the file and delete from memory
        pionfile->Close();
        delete pionfile;
    }
    // Save or display the canvas
    gPad->SetLogy(); 
    legend->Draw();
    canvas->Update();
    canvas->SaveAs(Form("pioncode/canvas_pdf/overlayed_histogram_%s.pdf",histName));
    // canvas->Print("overlayed_histograms.pdf"); // Uncomment to save as PDF

    // clean up
    delete legend; 
    delete canvas;
}

void FitYProjectionsAndGraph(TCanvas *canvas, TH2F *hist2D, const char *pdfname, double binres, int weightmethod){
    std::ofstream mycsv;
	//const int numXBins = 64; // Number of X bins
	int numXBins = hist2D->GetNbinsX();
	mycsv.open(Form("pioncode/csvfiles/Fit_Param_%s.csv", pdfname));
	// Create TGraphs to store the fit parameters
	//TGraph *ampGraph = new TGraph();
	//TGraph *meanGraph = new TGraph();
	//TGraph *sigmaGraph = new TGraph();
	TGraphErrors *ampGraph = new TGraphErrors(numXBins);
    TGraphErrors *meanGraph = new TGraphErrors(numXBins);
    TGraphErrors *sigmaGraph = new TGraphErrors(numXBins);
	// Loop over X bins
	int upper_limit, lower_limit;
	if (weightmethod==0){//exp
		lower_limit=1;
		upper_limit=binres*6;
	}
	else if (weightmethod==1){//power 
		lower_limit=binres*2;
		upper_limit=numXBins;
	}
	else if (weightmethod==2 ){//woods saxon
		lower_limit=1;
		upper_limit=numXBins;
	}
	for (int i = lower_limit; i <= upper_limit; ++i)// originally until numXBins// for exp use i < 4*6+1// for power go from 4*3+1 to numXBins+1
	{
		// Create a Y projection for the current X bin
		TH1D *yProjection = hist2D->ProjectionY(Form("YProjection_%d", i), i , i ,"");

		// Fit the Y projection with a Gaussian
		yProjection->Fit("gaus", "Q");
		
		// Access the fit parameters
		TF1 *fitFunc = yProjection->GetFunction("gaus");
		if(fitFunc){
			double amplitude = fitFunc->GetParameter(0);
			double mean = fitFunc->GetParameter(1);
			double sigma = fitFunc->GetParameter(2);
			double amplitudeError = fitFunc->GetParError(0);
			double meanError = fitFunc->GetParError(1);
			double sigmaError = fitFunc->GetParError(2);

			// Store the parameters in the TGraphs
			ampGraph->SetPoint(i, i/binres, amplitude);
			meanGraph->SetPoint(i, i/binres, mean);
			sigmaGraph->SetPoint(i, i/binres, sigma);
			ampGraph->SetPointError(i, 0, amplitudeError);// No error in X
			meanGraph->SetPointError(i, 0, meanError); // No error in X
			sigmaGraph->SetPointError(i, 0, sigmaError); // No error in 
			//*
			if(i==0){
			mycsv << "Bin Number" << "," << "," << "# entries" << "," << ","  << "Amplitude" << "," << "Error" << "," << "," << "Sigma" << "," << "Error" << "," << "," << "Mean" << "," << "Error" << "," << "," << "Sigma/Mean" << "\n";
			}
			else{
			mycsv << i << "," << "," << yProjection->GetEntries() << "," << "," << amplitude << ","<< amplitudeError << "," << "," << sigma << ","<< sigmaError<< "," << "," << mean << ","<< meanError << "," << "," << sigma/mean << "\n";
			}
			//*/
		}
        
	}
	mycsv.close();
	// Create a multi-panel canvas and plot the TGraphs
	canvas->Divide(2, 2); // Divide the canvas into 3 parts
	canvas->cd(1);
    ampGraph->SetTitle("Amplitude vs. Bin Number");
    ampGraph->GetXaxis()->SetTitle("Bin Number");
    ampGraph->GetYaxis()->SetTitle("Amplitude");
	//gPad->SetTopMargin(0.12);
    ampGraph->Draw("AP");
	canvas->cd(2);
    meanGraph->SetTitle("Mean vs. Bin Number");
    meanGraph->GetXaxis()->SetTitle("Bin Number");
    meanGraph->GetYaxis()->SetTitle("Mean");
	//gPad->SetTopMargin(0.12);
    meanGraph->Draw("AP");
	canvas->cd(3);
	sigmaGraph->SetTitle("Sigma vs. Bin Number");
    sigmaGraph->GetXaxis()->SetTitle("Bin Number");
    sigmaGraph->GetYaxis()->SetTitle("Sigma");
	//gPad->SetTopMargin(0.12);
    sigmaGraph->Draw("AP");
	canvas->cd(4);
	gPad->SetLogz();
	hist2D->Draw("colz");
	//return 0;

}

void SliceAndFit(filename_object filenameobj){

	time_t rawtime;
	struct tm *timeinfo;
	char Time[80]; // time string with size given
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(Time, 80, "%Y-%m-%d-%H:%M:%S", timeinfo);
	puts(Time);				   // print Time to root

	int rel_error_param = 155; // 65 for full set. I will do 145 now// switched from 67(old) to 65(new)


	std::vector<std::string> WeightNames = {"EXP", "POWER", "WSHP"};
	int weightmethod=2;//0=exp,1=power,2=wshp
	double binres=2;

	for (int i = 0; i < 1; i++)
	{ // 3 if doing 14.5-16.5//26 for full set(old)
		////////////////////////////////////
		if (rel_error_param > 320)
		{
			break;
		}

		// const char* canvasname;
		// cout << rel_error_param;//Printf(rel_error_param);
		const char *pionfilename;
		if (rel_error_param < 100)
		{
			pionfilename = Form("Pi0FastMC_0.0%d000_%s.root", rel_error_param,WeightNames[weightmethod].c_str());
		}
		else
		{
			pionfilename = Form("Pi0FastMC_0.%d000_%s.root", rel_error_param,WeightNames[weightmethod].c_str());
		}
		// const char* pionfilename = Form("Pi0FastMC_%d000.root", rel_error_param);
		cout << pionfilename << "\n";
		TFile *pionfile = new TFile(Form("pioncode/rootfiles/%s", pionfilename)); // TFile *_file0 = TFile::Open("A:/root_v6.26.06/bin/pioncode/rootfiles/Pi0FastMC_0.067000.root")  then .ls to see files

		// TFile* pionfile = new TFile(Form("pioncode/rootfiles/Pi0FastMC_0.067000.root", rel_error_param));//TFile *_file0 = TFile::Open("A:/root_v6.26.06/bin/pioncode/rootfiles/Pi0FastMC_0.067000.root")  then .ls to see files
		// TString dir = gROOT->GetTutorialDir();
		// dir.Append("/hsimple.C");
		// dir.ReplaceAll("/./", "/");
		// if (!gInterpreter->IsLoaded(pionfile.Data())) gInterpreter->LoadMacro(pionfile.Data());
		// TFile* hsimpleFile = (TFile*)gROOT->ProcessLineFast("hsimple(1)");
		// if (!hsimpleFile) return;
		TH2F *h9 = (TH2F *)pionfile->Get("h9"); // h18 weighted, h9 unweighted
		TH2F *h18 = (TH2F *)pionfile->Get("h18");
		// Create a canvas and divide it
		// int i = 0;
		// int canvas_param = rel_error_param + i;

		// delete pionfilename;

		// const char* canvasname = Form("weighted_unweighted_SliceFit_ErrParam_67_thousandths", Time);//_date_%s,rel_error_param
		TString canvasname = Form("%s_Sliced_%d_thousandths_%s",WeightNames[weightmethod].c_str(), rel_error_param, Time); //_date_%s
		const char *pdfname = canvasname;
		TCanvas *c1 = new TCanvas(canvasname, canvasname, 3000, 1200);
		// TCanvas* c1 = new TCanvas(Form("c_%d", rel_error_param), "c1", 3000, 1000);
		// TCanvas* c1 = new TCanvas("c1", "c1", 3000, 1000);
		// c1->SetFillColor(42);
		c1->Divide(2, 1); // nx, ny
		TPad *leftPad = (TPad *)c1->cd(1);
		;
		leftPad->Divide(2, 3); // nx, ny

		////////////////////////////////// unweighted
		leftPad->cd(1);
		// c1->cd(1);
		gPad->SetTopMargin(0.12);
		gPad->SetFillColor(33);
		h9->GetXaxis()->SetTitle("Pion Pt [GeV/c]");
		h9->GetYaxis()->SetTitle("Invariant Mass [GeV/c^2]");
		h9->Draw("colz");
		printf("error test code\n");
		h9->GetXaxis()->SetLabelSize(0.06);
		h9->GetYaxis()->SetLabelSize(0.06);
		h9->SetMarkerColor(kYellow);
		// Fit slices projected along Y fron bins in X [1,64] with more than 2 bins  in Y filled
		h9->FitSlicesY(0, 0, -1, 0);

		//

		leftPad->cd(2);
		gPad->SetFillColor(33);
		TH2F *h9_0 = (TH2F *)pionfile->Get("h9_0");
		h9_0->GetXaxis()->SetTitle("Pion Pt [GeV/c]");
		h9_0->Draw();
		// TPad* rightPad = (TPad*)c1->cd(2);
		// rightPad->Divide(1, 2);
		// rightPad->cd(1);
		// Show fitted "mean" for each slice //Double_t mm = h->GetMean();
		leftPad->cd(3);
		gPad->SetTopMargin(0.12);
		gPad->SetLeftMargin(0.15);
		gPad->SetFillColor(33);
		TH2F *h9_1 = (TH2F *)pionfile->Get("h9_1");
		h9_1->GetXaxis()->SetTitle("Pion Pt [GeV/c]");
		h9_1->GetYaxis()->SetTitle("Mean");
		h9_1->Draw();

		// Show fitted "sigma" for each slice
		// rightPad->cd(2);
		leftPad->cd(4);
		gPad->SetTopMargin(0.12);
		gPad->SetLeftMargin(0.15);
		gPad->SetFillColor(33);
		TH2F *h9_2 = (TH2F *)pionfile->Get("h9_2");
		h9_2->SetMinimum(0.8);
		h9_2->GetXaxis()->SetTitle("Pion Pt [GeV/c]");
		h9_2->GetYaxis()->SetTitle("Sigma");
		h9_2->Draw();
		///*
		// Show fitted variance(sigma^2) for each slice
		leftPad->cd(5);
		TH2F *h9_3 = (TH2F *)h9_2->Clone("h9_3"); // clone sigma
		gPad->SetTopMargin(0.12);
		gPad->SetLeftMargin(0.15);
		gPad->SetFillColor(33);
		h9_3->Multiply(h9_2); // multiply cloned sigma by sigma
		h9_3->SetTitle("Value of par[2]^2=Variance;Pion Pt [GeV/c];Variance");
		// h9_3->SetMinimum(0.8);
		h9_3->Draw();

		// Show fitted mean/variance for each slice
		leftPad->cd(6);
		TH2F *h9_4 = (TH2F *)h9_1->Clone("h9_4"); // clone mean
		gPad->SetTopMargin(0.12);
		gPad->SetLeftMargin(0.15);
		gPad->SetFillColor(33);
		h9_4->Divide(h9_2);
		h9_4->SetTitle("Value of par[1]/par[2]^2=Mean/Variance;Pion Pt [GeV/c];Mean/Variance");
		// h9_4->SetMinimum(0.8);
		h9_4->Draw();
		//*/
		/////////////////////////////////////////////weighted
		TPad *rightPad = (TPad *)c1->cd(2);
		rightPad->Divide(2, 3);
		rightPad->cd(1);
		// c1->cd(1);
		gPad->SetTopMargin(0.12);
		gPad->SetFillColor(33);
		gPad->SetLogz();
		h18->Draw("colz");
		printf("error test code\n");
		h18->GetXaxis()->SetLabelSize(0.06);
		h18->GetYaxis()->SetLabelSize(0.06);
		h18->GetXaxis()->SetTitle("Pion Pt [GeV/c]");
		h18->GetYaxis()->SetTitle("Invariant Mass [GeV/c^2]");
		h18->SetMarkerColor(kYellow);
		// Fit slices projected along Y fron bins in X [1,64] with more than 2 bins  in Y filled
		h18->FitSlicesY(0, 0, -1, 0);//, "EMW"

		//

		// Show fitted "mean" for each slice
		rightPad->cd(2);
		gPad->SetFillColor(33);
		TH2F *h18_0 = (TH2F *)pionfile->Get("h18_0");
		h18_0->GetXaxis()->SetTitle("Pion Pt [GeV/c]");
		h18_0->Draw();
		// TPad* rightPad = (TPad*)c1->cd(2);
		// rightPad->Divide(1, 2);
		// rightPad->cd(1);
		rightPad->cd(3);
		gPad->SetTopMargin(0.12);
		gPad->SetLeftMargin(0.15);
		gPad->SetFillColor(33);
		TH2F *h18_1 = (TH2F *)pionfile->Get("h18_1");
		h18_1->GetYaxis()->SetTitle("Mean");
		h18_1->GetXaxis()->SetTitle("Pion Pt [GeV/c]");
		// h18_1->SetAxisRange(0.1, 0.16,"Y");
		h18_1->Draw();

		// Show fitted "sigma" for each slice
		// rightPad->cd(2);
		rightPad->cd(4);
		gPad->SetTopMargin(0.12);
		gPad->SetLeftMargin(0.15);
		gPad->SetFillColor(33);
		TH2F *h18_2 = (TH2F *)pionfile->Get("h18_2");
		h18_2->SetMinimum(0.8);
		h18_2->GetYaxis()->SetTitle("Sigma");
		h18_2->GetXaxis()->SetTitle("Pion Pt [GeV/c]");
		h18_2->Draw();

		// Show fitted variance(sigma^2) for each slice
		rightPad->cd(5);
		TH2F *h18_3 = (TH2F *)h18_2->Clone("h18_3"); // clone sigma
		gPad->SetTopMargin(0.12);
		gPad->SetLeftMargin(0.15);
		gPad->SetFillColor(33);
		h18_3->Multiply(h18_2); // multiply cloned sigma by sigma
		// h9_3->SetMinimum(0.8);
		h18_3->SetTitle("Value of par[2]^2=Variance;Pion Pt [GeV/c];Variance");
		h18_3->Draw();

		// Show fitted mean/variance for each slice
		rightPad->cd(6);
		// TH2F* h18_4 = (TH2F*)h18_1->Clone("h18_4");//clone mean
		TH2F *h18_4 = (TH2F *)h18_2->Clone("h18_4"); // clone sigma instead
		gPad->SetTopMargin(0.12);
		gPad->SetLeftMargin(0.15);
		gPad->SetFillColor(33);
		h18_4->Divide(h18_3); // divide by variance. 2 is sigma
		h18_4->SetTitle("Value of par[1]/par[2]^2=Mean/Variance;Pion Pt [GeV/c];Mean/Variance");
		// h18_4->Divide(h18_1);//divide sigma by mean.
		// h18_4->SetTitle("Value of par[2]/par[1]=Sigma/Mean;Pion Pt [GeV/c];Sigma/Mean");
		// h18_4->SetMinimum(0.8);
		h18_4->Draw();

		// const char* canvasname = Form("weighted_unweighted_SliceFit_ErrParam_%d/1000_date_%s", rel_error_param, Time);
		// const char* canvasname = Form("pioncode/canvas_pdf/collected_canvas_ErrParam_%d/1000_Date_%s", rel_error_param, Time);
		cout << canvasname << "\n";
		// c1->SaveAs(Form("pioncode/canvas_pdf/Weighted_unweighted_SliceFit_ErrParam%s.pdf", canvasname));
		// c1->SaveAs(Form("pioncode/canvas_pdf/%s.pdf", canvasname));
		c1->SaveAs(Form("pioncode/canvas_pdf/%s.pdf", pdfname));
		// c1->SaveAs(Form("pioncode/canvas_pdf/Weighted_unweighted_SliceFit_ErrParam_%d/1000_Date_%s.pdf", rel_error_param, Time));
		delete c1;

		if (rel_error_param == 155)
		{
			TCanvas *c2 = new TCanvas("c2", "c2", 400, 900);
			c2->Divide(1, 2);
			c2->cd(1);
			TH2F *h18_5 = (TH2F *)h18_2->Clone("h18_5");
			h18_5->SetAxisRange(0., 16., "x");
			// h18_5->Scale(1000/135);
			// h18_5->SetAxisRange(0., 0.2,"y");
			h18_5->Draw();
			h18_5->GetYaxis()->SetTitle("Sigma");
			h18_5->GetXaxis()->SetTitle("Pion Pt [GeV/c]");

			c2->cd(2);
			/* they said  "A beam momentum spread (δp/p ≈ 2%) is quadratically subtracted from σ/μ of the fit, in order to unfolded beam momentum spread from the relative energy resolution. The Gauss function parameter of μ and energy resolution from each fit are plotted against the  nominal beam energy as linearity and resolution." //*/
			TH2F *h18_6 = (TH2F *)h18_2->Clone("h18_6"); // sigma/mean

			h18_6->Divide(h18_1); // sigma/mean

			// TH2F* h18_7 = (TH2F*)h18_6->Clone("h18_7");// clone (sigma/mean)to get (sigma/mean)^2
			// h18_7->Multiply(h18_6);//(sigma/mean)^2
			// h18_7->Add(-0.0004);//why?
			// TH2F* h18_8 = (TH2F*)h18_7->Clone("h18_8");
			// h18_8->Divide(h18_7);// call this the new sigma/mean
			// h18_8->Multiply(h18_1);//multiply by the mean to find a new Sigma
			// h18_8->Divide(135);//scale it to find sigma_M/M directly

			h18_6->SetAxisRange(0., 16., "x");
			h18_6->SetAxisRange(0., 0.2, "y");
			h18_6->Draw();
			h18_6->GetYaxis()->SetTitle("Sigma/Mean");
			h18_6->GetXaxis()->SetTitle("Pion Pt [GeV/c]");
			//*/
			c2->SaveAs(Form("pioncode/canvas_pdf/%s_truncatedsigma.pdf", pdfname));
			delete c2;
		}

		
		TCanvas *c3 = new TCanvas("c3", "c3", 3000, 3000);
		
		FitYProjectionsAndGraph(c3, h18, pdfname, binres, weightmethod);
		c3->SaveAs(Form("pioncode/canvas_pdf/Alt_Projection_%s.pdf", pdfname));
		
		pionfile->Close();

		rel_error_param = rel_error_param + 10;
		cout << rel_error_param << "\n";
		delete c3;
	}
	//return 0;
}

//----------------------------------all in one file functions
void OverlayMeansAIO(filename_object filenameobj, const char* histName) {
//void OverlayMeans(const std::vector<std::string> & fileNames, string filenamemod) {
    // Create a TCanvas
    TCanvas* canvas1 = new TCanvas("canvas1", "Overlay Means", 800, 600);
    
    //canvas1->SetGrid();
    //gStyle->SetOptTitle ( 0 );
    //gStyle->SetOptStat ( 0 );
    //gStyle->SetOptFit ( 1 1 1 1 );
    //gStyle->SetStatX ( . 8 9 ); 
    //gStyle->SetStatY ( . 8 9 );
    //gStyle->SetStatBorderSize ( 0 );
    float errparam=filenameobj.sqrtEsmearing[0];
    double binres=2;//number of divisions per GeV
    // Create a legend
    TLegend* legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);//0.7, 0.4, 0.9, 0.6
    TMultiGraph *MultiGraphs = new TMultiGraph();//h18->GetNbinsX()
    // Loop over each file
    for (size_t i = 0; i < filenameobj.fileNames.size(); ++i) {
        // Open the root file
        TFile* pionfile = new TFile(filenameobj.fileNames[i].c_str(), "READ");

        // Check if the file is open
        if (!pionfile || pionfile->IsZombie()) {
            std::cerr << "Error: Unable to open file " << filenameobj.fileNames[i] << std::endl;
            continue;
        }

        // Retrieve the 2D histogram
        //TH2F *h18 =(TH2F *)pionfile->Get("h18");
        TH2F* htemp = nullptr;
        pionfile->GetObject(histName, htemp);

        // Check if the histogram exists
        if (!htemp) {
            std::cerr << "Error: Histogram 'myHist2D' not found in file " << filenameobj.fileNames[i] << std::endl;
            pionfile->Close();
            continue;
        }

        // Create a histogram for means
        TGraphErrors *meanGraph = new TGraphErrors(htemp->GetNbinsX());//

        // Loop over each bin in the X direction
        for (int binX = 1; binX <= htemp->GetNbinsX(); ++binX) {
            //std::cout << "Nbin" << " " << binX << std::endl;
            // Project along Y for each sbinX
            TH1D* yProjection = htemp->ProjectionY(Form("YProjection_%zu_%d", i, binX), binX, binX, "");

            // Fit the Y projection with a Gaussian
            yProjection->Fit("gaus", "Q");

            // Access the fit parameters
            TF1* fitFunc = yProjection->GetFunction("gaus");

            // Check if the fit function is valid
            if (fitFunc) {
                // Fill the mean histogram with the mean value
                //meanHistogram->SetBinContent(binX, fitFunc->GetParameter(1));
                if (filenameobj.weightnames[i]=="EXP"){//exp &&  binX <= filenameobj.pTcutoff*binres
                    meanGraph->SetPoint(binX, binX/binres,fitFunc->GetParameter(1));
                    meanGraph->SetPointError(binX, 0,fitFunc->GetParError(1));
                    if(binX/binres==3){
                        std::cout << binX <<" "<<fitFunc->GetParameter(1) << std::endl; // debug line
                    }
                }
                else if (filenameobj.weightnames[i]=="POWER"){//power &&  3*binres < binX
                    meanGraph->SetPoint(binX, binX/binres,fitFunc->GetParameter(1));
                    meanGraph->SetPointError(binX, 0,fitFunc->GetParError(1));
                    if(binX/binres==3){
                        std::cout << binX <<" "<<fitFunc->GetParameter(1) << std::endl; // debug line
                    }
                    //std::cout << binX <<" "<<fitFunc->GetParameter(1) << std::endl; // debug line
                }
                else if (filenameobj.weightnames[i]=="WSHP"){//woods saxon
                    meanGraph->SetPoint(binX, binX/binres,fitFunc->GetParameter(1));
                    meanGraph->SetPointError(binX, 0,fitFunc->GetParError(1));
                }
                else if (filenameobj.weightnames[i]=="HAGEDORN"){//HAGEDORN
                    meanGraph->SetPoint(binX, binX/binres,fitFunc->GetParameter(1));
                    meanGraph->SetPointError(binX, 0,fitFunc->GetParError(1));
                }
                
                //meanHistogram->SetBinError(binX, fitFunc->GetParError(1));
                //std::cout << "bin number" << " " << binX << " " << "Mean" << " " <<  fitFunc->GetParameter(1) << " " << "Mean error" << " " << fitFunc->GetParError(1) << std::endl;
            }
                // Add an entry to the legend
                //legend1->AddEntry(yProjection, Form("Version %zu, BinX %d", i, binX), "L");
            
            //std::cout << "I reached here, pre delete proj" << std::endl; // debug line
            // Clean up Y projection
            delete yProjection;
        }
        // Set different line colors for each version
        int MarkerStyle = i + 24; // 
        int MarkerColor = i + 1;
        //meanGraph->SetLineColor(lineColor);
        meanGraph->SetMarkerStyle(MarkerStyle);
        meanGraph->SetMarkerColor(MarkerColor);
        

       
        // Overlay the mean histogram on the same canvas
        if (i == 0) {
            //meanGraph->Draw("AP"); // Draw histogram for the first version
            //canvas1->Print("OverlayMeanHistograms.pdf");
            MultiGraphs->Add(meanGraph,"PE");
            std::cout << "draw for the first file" << std::endl; // debug line
        } else {
            MultiGraphs->Add(meanGraph,"PE");
            //meanGraph->Draw("P SAME"); // Draw subsequent histograms on the same canvas
            std::cout << "draw for subsequent" << std::endl; // debug line
        }   
        
        MultiGraphs->SetTitle(Form("Smeared Inv. Mass%s;pT (GeV);Inv. Mass (GeV)",filenameobj.canvasnamemod.c_str()));
        MultiGraphs->Draw("APE");

        // Add an entry to the legend
        //std::vector<std::string> legendstring = {"EXP","POWER","WSHP"};
        legend1->AddEntry(meanGraph, filenameobj.legendnames[i].c_str(), "P");

        // Close the file
        pionfile->Close();
        delete pionfile;

        //std::cout << "I reached here, close+delete file" << std::endl; // debug line
    }

    //std::cout << "I reached here, done with all files" << std::endl; // debug line
    // Draw the legend
    legend1->Draw();

    // Show the canvas
    MultiGraphs->GetXaxis()->SetLimits(filenameobj.plotxlims[0],filenameobj.plotxlims[1]);
    MultiGraphs->SetMinimum(filenameobj.plotylims[0]);
    MultiGraphs->SetMaximum(filenameobj.plotylims[1]);
    canvas1->SetMargin(0.2,0.1,0.1,0.1);
    gPad->Modified();
    gPad->Update();
    //canvas1->Modified();
    canvas1->Update();
    canvas1->Print(Form("pioncode/canvas_pdf/%s_%s_%f_OverlayMeanHistograms.pdf",histName,filenameobj.filenamemod.c_str(),errparam));

    // Clean up
    delete canvas1;
    delete legend1;
}


void OverlaySigmaOverMeanAIO(filename_object filenameobj, const char* histName){
//void OverlaySigmaOverMean(const std::vector<std::string>& fileNames, string filenamemod) {
    // Create a TCanvas
    TCanvas* canvas1 = new TCanvas("canvas1", "Overlay Means", 800, 600);
    //canvas1->SetGrid();
    //gStyle->SetOptTitle ( 0 );
    //gStyle->SetOptStat ( 0 );
    //gStyle->SetOptFit ( 1 1 1 1 );
    //gStyle->SetStatX ( . 8 9 ); 
    //gStyle->SetStatY ( . 8 9 );
    //gStyle->SetStatBorderSize ( 0 );

    float errparam=filenameobj.sqrtEsmearing[0];
    double binres=2;//number of divisions per GeV
    // Create a legend
    TLegend* legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    TMultiGraph *MultiGraphs = new TMultiGraph();//h18->GetNbinsX()
    // Loop over each file
    for (size_t i = 0; i < filenameobj.fileNames.size(); ++i) {
        // Open the root file
        TFile* pionfile = new TFile(filenameobj.fileNames[i].c_str(), "READ");

        // Check if the file is open
        if (!pionfile || pionfile->IsZombie()) {
            std::cerr << "Error: Unable to open file " << filenameobj.fileNames[i] << std::endl;
            continue;
        }

        // Retrieve the 2D histogram
        //TH2F *h18 =(TH2F *)pionfile->Get("h18");
        TH2F* htemp = nullptr;
        pionfile->GetObject(histName, htemp);

        // Check if the histogram exists
        if (!htemp) {
            std::cerr << "Error: Histogram 'myHist2D' not found in file " << filenameobj.fileNames[i] << std::endl;
            pionfile->Close();
            continue;
        }

        // Create a histogram for means
        TGraphErrors *meanGraph = new TGraphErrors(htemp->GetNbinsX());
        //TH1F* meanHistogram = new TH1F(Form("MeanHistogram_%zu", i), Form("Version %zu", i), h18->GetNbinsX(), 0.5, h18->GetNbinsX() + 0.5);



        // Loop over each bin in the X direction
        for (int binX = 1; binX <= htemp->GetNbinsX(); ++binX) {
            // Project along Y for each binX
            TH1D* yProjection = htemp->ProjectionY(Form("YProjection_%zu_%d", i, binX), binX, binX, "");

            // Fit the Y projection with a Gaussian
            yProjection->Fit("gaus", "Q");

            // Access the fit parameters
            TF1* fitFunc = yProjection->GetFunction("gaus");

            // Check if the fit function is valid
            if (fitFunc) {
                // Fill the mean histogram with the mean value
                //meanHistogram->SetBinContent(binX, fitFunc->GetParameter(1));
                if (filenameobj.weightnames[i]=="EXP"){//exp && binX <= filenameobj.pTcutoff*binres
                    double meanoversigma =fitFunc->GetParameter(2)/fitFunc->GetParameter(1);
                    double meanoversigmaerr = meanoversigma*(fitFunc->GetParError(2)/fitFunc->GetParameter(2)+fitFunc->GetParError(1)/fitFunc->GetParameter(1));//m/s*(serr/s+merr/m)

                    meanGraph->SetPoint(binX, binX/binres,meanoversigma);
                    meanGraph->SetPointError(binX, 0,meanoversigmaerr);
                }
                else if (filenameobj.weightnames[i]=="POWER"){//power 
                    double meanoversigma =fitFunc->GetParameter(2)/fitFunc->GetParameter(1);
                    double meanoversigmaerr=meanoversigma*(fitFunc->GetParError(2)/fitFunc->GetParameter(2)+fitFunc->GetParError(1)/fitFunc->GetParameter(1));//m/s*(serr/s+merr/m)

                    meanGraph->SetPoint(binX, binX/binres,meanoversigma);
                    meanGraph->SetPointError(binX, 0,meanoversigmaerr);
                    std::cout << binX <<" "<<meanoversigma << std::endl; // debug line

                }
                else if (filenameobj.weightnames[i]=="WSHP") {//woods saxon
                    double meanoversigma =fitFunc->GetParameter(2)/fitFunc->GetParameter(1);
                    double meanoversigmaerr=meanoversigma*(fitFunc->GetParError(2)/fitFunc->GetParameter(2)+fitFunc->GetParError(1)/fitFunc->GetParameter(1));//(m/s)_err=m/s*(serr/s+merr/m)

                    meanGraph->SetPoint(binX, binX/binres,meanoversigma);
                    meanGraph->SetPointError(binX, 0,meanoversigmaerr);
                }
                else if (filenameobj.weightnames[i]=="HAGEDORN") {//HAGEDORN
                    double meanoversigma =fitFunc->GetParameter(2)/fitFunc->GetParameter(1);
                    double meanoversigmaerr=meanoversigma*(fitFunc->GetParError(2)/fitFunc->GetParameter(2)+fitFunc->GetParError(1)/fitFunc->GetParameter(1));//(m/s)_err=m/s*(serr/s+merr/m)

                    meanGraph->SetPoint(binX, binX/binres,meanoversigma);
                    meanGraph->SetPointError(binX, 0,meanoversigmaerr);
                }

                
                //meanHistogram->SetBinError(binX, fitFunc->GetParError(1));
                //std::cout << "bin number" << " " << binX << " " << "Mean" << " " <<  fitFunc->GetParameter(1) << " " << "Mean error" << " " << fitFunc->GetParError(1) << std::endl;
            }
                // Add an entry to the legend
                //legend1->AddEntry(yProjection, Form("Version %zu, BinX %d", i, binX), "L");
            
            //std::cout << "I reached here, pre delete proj" << std::endl; // debug line
            // Clean up Y projection
            delete yProjection;
        }
        // Set different line colors for each version
        int MarkerStyle = i + 24; // 
        int MarkerColor = i + 1;
        //meanGraph->SetLineColor(lineColor);
        meanGraph->SetMarkerStyle(MarkerStyle);
        meanGraph->SetMarkerColor(MarkerColor);

        //std::cout << "I reached here, done with loop over bins" << std::endl; // debug line

        // Overlay the mean histogram on the same canvas
        if (i == 0) {
            //meanGraph->Draw("AP"); // Draw histogram for the first version
            //canvas1->Print("OverlayMeanHistograms.pdf");
            MultiGraphs->Add(meanGraph,"PE");
            std::cout << "draw for the first file" << std::endl; // debug line
        } else {
            MultiGraphs->Add(meanGraph,"PE");
            //meanGraph->Draw("P SAME"); // Draw subsequent histograms on the same canvas
            std::cout << "draw for subsequent" << std::endl; // debug line
        }   

        MultiGraphs->SetTitle(Form(" Sigma/Mean (Smeared Inv. Mass)%s;pT (GeV);sigma/mean)",filenameobj.canvasnamemod.c_str()));
        MultiGraphs->Draw("APE");
        // Add an entry to the legend
        //std::vector<std::string> legendstring = {"EXP","POWER","WSHP"};
        legend1->AddEntry(meanGraph, filenameobj.legendnames[i].c_str(), "P");

        // Close the file
        pionfile->Close();
        delete pionfile;

        std::cout << "I reached here, close+delete file" << std::endl; // debug line
    }

    std::cout << "I reached here, done with all files" << std::endl; // debug line
    // Draw the legend
    legend1->Draw();
    
    MultiGraphs->GetXaxis()->SetLimits(filenameobj.plotxlims[0],filenameobj.plotxlims[1]);
    MultiGraphs->SetMinimum(filenameobj.plotylims[2]);
    MultiGraphs->SetMaximum(filenameobj.plotylims[3]);
    canvas1->SetMargin(0.2,0.1,0.1,0.1);
    gPad->Modified();
    gPad->Update();

   //canvas1->Update();
    canvas1->Modified();
    canvas1->SaveAs(Form("pioncode/canvas_pdf/%s_%s_%f_OverlaySigma_Over_Mean.pdf",histName,filenameobj.filenamemod.c_str(),errparam));//Print-> works too
    
    // Clean up
    delete canvas1;
    delete legend1;
}

