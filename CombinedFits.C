#include <iostream>
#include <vector>
#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TF1.h>
#include <TGraphMultiErrors.h>

class filename_object {//object to hold my file names and related strings. things like settings for the files are contained here.
    public:
    std::vector<std::string> fileNames;
    std::vector<std::string> legendnames;
    string filenamemod;
    string canvasnamemod;
    std::vector<float> plotxlims;
    std::vector<float> plotylims;
    int pTcutoff;
    //float minplotpT;

};

filename_object choosecomparisontype(int choosetype){
    filename_object filename_object1;// 0=weight type, 1=ac on/off, 2=co on/off, 3=ac&co on/off
    if(choosetype==0){
        //filename_object weightfilenameobj;
        filename_object1.fileNames={"pioncode/rootfiles/Pi0FastMC_0.155000_EXP.root", "pioncode/rootfiles/Pi0FastMC_0.155000_POWER.root", "pioncode/rootfiles/Pi0FastMC_0.155000_WSHP.root"};
        filename_object1.legendnames={"EXP","POWER","WSHP"};
        filename_object1.filenamemod="weightmethod";
        filename_object1.canvasnamemod=" for various weighting methods";  
        filename_object1.plotxlims={0.9,6.4};//min, max
        filename_object1.plotylims={0.13,0.17,0.08,0.25}; //mean_min, mean_max,sm_min, sm_max
        filename_object1.pTcutoff=6;
    }
    else if(choosetype==1){
        filename_object1.fileNames={"pioncode/rootfiles/Pi0FastMC_0.155000_WSHP_ac0_co0.root", "pioncode/rootfiles/Pi0FastMC_0.155000_WSHP_ac1_co0.root"};
        filename_object1.legendnames={"Asym. cut off","Asym. cut on"};
        filename_object1.filenamemod="AsymCutTest";
        filename_object1.canvasnamemod=", Asymmetry Cut: on vs off";
        filename_object1.plotxlims={0.1,16.4};//min, max
        filename_object1.plotylims={0.125,0.15,0.01,0.4}; //mean_min, mean_max,sm_min,sm_max  
        filename_object1.pTcutoff=16;
    }
    else if(choosetype==2){
        filename_object1.fileNames={"pioncode/rootfiles/Pi0FastMC_0.155000_WSHP_ac0_co0.root", "pioncode/rootfiles/Pi0FastMC_0.155000_WSHP_ac0_co1.root"};
        filename_object1.legendnames={"Cluster Overlay off","Cluster Overlay on"};
        filename_object1.filenamemod="ClusterOverlayTest";
        filename_object1.canvasnamemod=", Cluster Overlap: on vs off";
        filename_object1.plotxlims={0.1,16.4};//min, max
        filename_object1.plotylims={0.135,0.145,0.05,0.25}; //mean_min, mean_max,sm_min,sm_max 
        filename_object1.pTcutoff=16; 
    }
    else if(choosetype==3){
        filename_object1.fileNames={"pioncode/rootfiles/Pi0FastMC_0.155000_WSHP_ac0_co0.root", "pioncode/rootfiles/Pi0FastMC_0.155000_WSHP_ac1_co0.root", "pioncode/rootfiles/Pi0FastMC_0.155000_WSHP_ac0_co1.root", "pioncode/rootfiles/Pi0FastMC_0.155000_WSHP_ac1_co1.root"};
        filename_object1.legendnames={"AC,CO:00","AC,CO:10","AC,CO:01","AC,CO:11"};
        filename_object1.filenamemod="AsymClusterTest";
        filename_object1.canvasnamemod=", Cluster Overlap & Asymmmetry Cut: on(1) vs off(0)";
        filename_object1.plotxlims={0.1,16.4};//min, max
        filename_object1.plotylims={0.125,0.15,0.01,0.4}; //mean_min, mean_max,sm_min,sm_max  
        filename_object1.pTcutoff=16; 
    }

return filename_object1;
}

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
                if (i==0 &&  binX <= filenameobj.pTcutoff*binres){//exp
                    meanGraph->SetPoint(binX, binX/binres,fitFunc->GetParameter(1));
                    meanGraph->SetPointError(binX, 0,fitFunc->GetParError(1));
                    if(binX/binres==3){
                        std::cout << binX <<" "<<fitFunc->GetParameter(1) << std::endl; // debug line
                    }
                }
                else if (i==1){//power &&  3*binres < binX
                    meanGraph->SetPoint(binX, binX/binres,fitFunc->GetParameter(1));
                    meanGraph->SetPointError(binX, 0,fitFunc->GetParError(1));
                    if(binX/binres==3){
                        std::cout << binX <<" "<<fitFunc->GetParameter(1) << std::endl; // debug line
                    }
                    //std::cout << binX <<" "<<fitFunc->GetParameter(1) << std::endl; // debug line
                }
                else if (i==2){//woods saxon
                    meanGraph->SetPoint(binX, binX/binres,fitFunc->GetParameter(1));
                    meanGraph->SetPointError(binX, 0,fitFunc->GetParError(1));
                }
                else if (i==3){//woods saxon
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
    canvas1->Print(Form("pioncode/canvas_pdf/%s_OverlayMeanHistograms.pdf",filenameobj.filenamemod.c_str()));

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
                if (i==0 && binX <= filenameobj.pTcutoff*binres){//exp
                    double meanoversigma =fitFunc->GetParameter(2)/fitFunc->GetParameter(1);
                    double meanoversigmaerr = meanoversigma*(fitFunc->GetParError(2)/fitFunc->GetParameter(2)+fitFunc->GetParError(1)/fitFunc->GetParameter(1));//m/s*(serr/s+merr/m)

                    meanGraph->SetPoint(binX, binX/binres,meanoversigma);
                    meanGraph->SetPointError(binX, 0,meanoversigmaerr);
                }
                else if (i==1){//power 
                    double meanoversigma =fitFunc->GetParameter(2)/fitFunc->GetParameter(1);
                    double meanoversigmaerr=meanoversigma*(fitFunc->GetParError(2)/fitFunc->GetParameter(2)+fitFunc->GetParError(1)/fitFunc->GetParameter(1));//m/s*(serr/s+merr/m)

                    meanGraph->SetPoint(binX, binX/binres,meanoversigma);
                    meanGraph->SetPointError(binX, 0,meanoversigmaerr);
                    std::cout << binX <<" "<<meanoversigma << std::endl; // debug line

                }
                else if (i==2) {//woods saxon
                    double meanoversigma =fitFunc->GetParameter(2)/fitFunc->GetParameter(1);
                    double meanoversigmaerr=meanoversigma*(fitFunc->GetParError(2)/fitFunc->GetParameter(2)+fitFunc->GetParError(1)/fitFunc->GetParameter(1));//(m/s)_err=m/s*(serr/s+merr/m)

                    meanGraph->SetPoint(binX, binX/binres,meanoversigma);
                    meanGraph->SetPointError(binX, 0,meanoversigmaerr);
                }
                else if (i==3) {//woods saxon
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
        std::vector<std::string> legendstring = {"EXP","POWER","WSHP"};
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
    canvas1->SaveAs(Form("pioncode/canvas_pdf/%s_OverlaySigma_Over_Mean.pdf",filenameobj.filenamemod.c_str()));//Print-> works too
    
    // Clean up
    delete canvas1;
    delete legend1;
}

void CombinedFits() {
    filename_object choosenfilenameobj = choosecomparisontype(3);// 0=weight type, 1=ac on/off, 2=co on/off, 4 ac&co on/off
    OverlayMeans(choosenfilenameobj);
    OverlaySigmaOverMean(choosenfilenameobj);
}    


//return 0;
// List of root file names
    //std::vector<std::string> fileNames = {"pioncode/rootfiles/Pi0FastMC_0.155000_EXP.root", "pioncode/rootfiles/Pi0FastMC_0.155000_POWER.root", "pioncode/rootfiles/Pi0FastMC_0.155000_WSHP.root"};//{"pioncode/Pi0FastMC_0.155000WSHP.root"};
    // {"pioncode/Pi0FastMC_0.155000EXP.root", "pioncode/Pi0FastMC_0.155000POWER.root", "pioncode/Pi0FastMC_0.15500WSHP.root"};
    //std::string filenamemod = "weightmethod";
    //std::string filenamemod = "weightmethod";
    // Overlay the means for each bin
    //OverlayMeans(fileNames,filenamemod);
    //OverlaySigmaOverMean(fileNames,filenamemod);