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
    string filenamemod;
    string canvasnamemod;
    std::vector<float> plotxlims;
    std::vector<float> plotylims;
    int pTcutoff;
    std::vector<float> sqrtEsmearing;

};

float extractNumber(const std::string& filePath);
filename_object choosecomparisontype(int choosetype);
void OverlayMeans(filename_object filenameobj);
void OverlaySigmaOverMean(filename_object filenameobj);
void ClusterOverlayTestFunc(filename_object filenameobj, const char* histname);
TH1D* getYProjectionof2DHist(const char* fileName, const char* histName, int firstxbin, int lastxbin);
void plotOverlayedHistograms(filename_object filenameobj, const char* histName);


void CombinedFits() {
    filename_object choosenfilenameobj = choosecomparisontype(0);// 0=weight type, 1=ac on/off, 2=co on/off, 4 ac&co on/off
    OverlayMeans(choosenfilenameobj);
    OverlaySigmaOverMean(choosenfilenameobj);
    plotOverlayedHistograms(choosenfilenameobj, "h3");//h12 is smeared pion pT, Weighted. h3 is unsmeared pion pT, weighted
    ///*
    if (choosenfilenameobj.fileNames.size()==2){// for subtraction of inv mass profile
        ClusterOverlayTestFunc(choosenfilenameobj, "h18");//const char* histname
    }
    //*/
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
    else if(choosetype==4){//co for 3% smearing
        filename_object1.fileNames={"pioncode/rootfiles/Pi0FastMC_0.030000_WSHP_ac0_co0.root", "pioncode/rootfiles/Pi0FastMC_0.030000_WSHP_ac0_co1.root"};
        filename_object1.legendnames={"Cluster Overlay off","Cluster Overlay on"};
        filename_object1.filenamemod="ClusterOverlayTest";
        filename_object1.canvasnamemod=", Cluster Overlap: on vs off";
        filename_object1.plotxlims={0.1,16.4};//min, max
        filename_object1.plotylims={0.1,0.137,0.0,0.085}; //mean_min, mean_max,sm_min,sm_max 
        filename_object1.pTcutoff=16; 
    }

    for(size_t i=0; i < filename_object1.fileNames.size(); i++){
        filename_object1.sqrtEsmearing.push_back(extractNumber(filename_object1.fileNames[i]));
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
            std::cerr << "Error: Could not retrieve histogram" << std::endl;
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


