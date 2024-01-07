#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <functional>
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
#include <TPad.h>
#include <TImage.h>

class filename_object {//object to hold my file names and related strings. things like settings for the files are contained here.
    public:
    std::vector<std::string> fileNames;
    std::vector<std::string> legendnames;
    std::vector<std::string> weightnames;
    string filenamemod;
    string canvasnamemod;
    std::vector<float> plotxlims;
    std::vector<float> plotylims;
    int pTcutoff;
    std::vector<float> sqrtEsmearing;
    int binres;
};
// declarations
float extractNumber(const std::string& filePath);
filename_object choosecomparisontype(int choosetype);
TCanvas* FitMeanAndPlot(filename_object filenameobj, int legendInt, const std::string& fileName, std::vector<std::string> HistList,std::vector<std::string> HistLegend);
void GraphAndSaveToPDF(filename_object filenameobj, std::vector<std::string> HistList, std::vector<std::string> HistLegend);
void ClusterOverlayTestFunc(filename_object filenameobj, const std::string& fileName, const char* histName, const std::string& LegendName);
TH1D* getYProjectionof2DHist(const char* fileName, const char* histName, int firstxbin, int lastxbin);



void AIOFit() {
    // 0=weight type
    int fileset = 0 ;
    filename_object choosenfilenameobj = choosecomparisontype(fileset);
    std::vector<std::string> HistList={"h18_","h27_","h29_","h28_"};
    std::vector<std::string> HistLegend={"Smeared Pion pT vs Inv Mass","Smeared Pion pT vs Inv Mass. cluster","Smeared Pion pT vs Inv Mass. asymm cut","Smeared Pion pT vs Inv Mass. clust+asymm"};
    GraphAndSaveToPDF(choosenfilenameobj,  HistList, HistLegend);
    
    //OverlayMeans(choosenfilenameobj);
    //OverlaySigmaOverMean(choosenfilenameobj);
    //plotOverlayedHistograms(choosenfilenameobj, "h12");//h12 is smeared pion pT, Weighted. h3 is unsmeared pion pT, weighted
    //SliceAndFit(choosenfilenameobj);

    //ClusterOverlayTestFunc(choosenfilenameobj,"pioncode/rootfiles/Pi0FastMC_0.155000.root", "h27_2", "test");
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
        filename_object1.fileNames={"pioncode/rootfiles/Pi0FastMC_0.155000.root"};
        filename_object1.legendnames={"EXP","POWER","WSHP","HAGEDORN"};
        filename_object1.weightnames={"EXP","POWER","WSHP","HAGEDORN"};
        filename_object1.filenamemod="weightmethod_co1_ac1";
        //filename_object1.canvasnamemod=" for various weighting methods. asymm+clustering";  
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

void GraphAndSaveToPDF(filename_object filenameobj, std::vector<std::string> HistList,std::vector<std::string> HistLegend) {
    // create canvas and open the output pdf
    TCanvas *canvas = new TCanvas("canvas1", "Canvas", 800, 600);
    //TCanvas* canvas;    
    canvas->Print("pioncode/canvas_pdf/output.pdf[");
    //loop over files. save 
    //int legendInt=0;
    
    for (size_t l = 0; l < filenameobj.weightnames.size(); ++l) {
        std::vector<std::string> histogramName;
        int legendInt=0;
        // Construct the histogram name with the index appended
        for (size_t v = 0; v < HistList.size(); ++v) {
            histogramName.push_back(HistList[v] + std::to_string(l));
            std::cout << histogramName[v] << std::endl; // debug line
        }

        for (const auto& fileName : filenameobj.fileNames) {
                canvas = FitMeanAndPlot(filenameobj, legendInt, fileName, histogramName, HistLegend);
                // Save canvas to PDF as a page
                canvas->Print("pioncode/canvas_pdf/output.pdf");
                legendInt++;
                std::cout << "filename loop done" << std::endl; // debug line
        }
       // Clear the histogramName vector for the next iteration
        histogramName.clear();
        std::cout << "l loop done" << std::endl; // debug line
    }
    
    // Close the PDF file
    canvas->Print("pioncode/canvas_pdf/output.pdf]");
    //clean up
    delete canvas;
}

TCanvas* FitMeanAndPlot(filename_object filenameobj, int legendInt, const std::string& fileName, std::vector<std::string> HistList,std::vector<std::string> HistLegend) {

    // create canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    // Create a legend
    TLegend* legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);

    TFile *file = new TFile(fileName.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Unable to open file " << fileName << std::endl;
        delete c1;
        delete legend1;
        return nullptr;
    }

    //plotting operations.
    float errparam=filenameobj.sqrtEsmearing[0];
    double binres=2;//number of divisions per GeV

    //int numhists=HistList.size();
    // Declare Histograms
    std::vector<TH2F*> temphist(HistList.size());
    // TMultiGraph to hold graphs.
    TMultiGraph *MultiGraphs = new TMultiGraph();
    // Create Tgraphs to hold means for each histogram
    std::vector<TGraphErrors*> meanGraph(HistList.size());

    // loop over the two histograms being compared. probably a better way to do this since it is a 1d array.
    for (size_t j=0; j < HistList.size(); j++){
        std::cout << "file" << j << " of" << HistList.size() << std::endl; // debug line
        // Load histograms from the file
        temphist[j] = dynamic_cast<TH2F*>(file->Get(HistList[j].c_str()));
        meanGraph[j] = new TGraphErrors(temphist[j]->GetNbinsX());
        std::cout << "loaded hists? made graph" << std::endl; // debug line
        if (!temphist[j] ) {
            std::cerr << "Error: Unable to retrieve histogram "<< j << " from file " << fileName << std::endl;
            file->Close();
            delete c1;
            return nullptr;
        }
        std::cout << "pre bin loop" << std::endl; // debug line
        // Loop over each bin in the X direction
        for (int binX = 1; binX <= temphist[j]->GetNbinsX(); binX++) {
            
            // Project along Y for each binX
            TH1D* yProjection = temphist[j]->ProjectionY(Form("YProjection_%zu_%d", j, binX), binX, binX, "");

            // Fit the Y projection with a Gaussian
            yProjection->Fit("gaus", "Q");

            // Access the fit parameters
            TF1* fitFunc = yProjection->GetFunction("gaus");

            // Check if the fit function is valid
            if (fitFunc) {
                meanGraph[j]->SetPoint(binX, binX/binres,fitFunc->GetParameter(1));
                meanGraph[j]->SetPointError(binX, 0,fitFunc->GetParError(1));
                // Set the legend entry to the title of the original 2D histogram
                //meanGraph[j]->SetTitle(temphist[j]->GetTitle());
            }
            // Clean up Y projection
            delete fitFunc;
            delete yProjection;
        }
        std::cout << "bin loop done" << std::endl; // debug line
        // Set different line colors for each version
        int MarkerStyle = j + 24; // 
        int MarkerColor = j + 1;
        //meanGraph->SetLineColor(lineColor);
        meanGraph[j]->SetMarkerStyle(MarkerStyle);
        meanGraph[j]->SetMarkerColor(MarkerColor);
        // Overlay the mean histogram on the same canvas
        if (j == 0) {
            //meanGraph->Draw("AP"); // Draw histogram for the first version
            //canvas1->Print("OverlayMeanHistograms.pdf");
            MultiGraphs->Add(meanGraph[j],"PE");
            std::cout << "draw for the first file" << std::endl; // debug line
        } else {
            MultiGraphs->Add(meanGraph[j],"PE");
            //meanGraph->Draw("P SAME"); // Draw subsequent histograms on the same canvas
            std::cout << "draw for subsequent" << std::endl; // debug line
        }       
    // Add an entry to the legend
    legend1->AddEntry(meanGraph[j], HistLegend[j].c_str(), "P");//meanGraph[j]->GetTitle(),"P"
    }

    std::cout << "position 1" << std::endl; // debug line
    MultiGraphs->SetTitle(Form("Smeared Pion pT vs Inv Mass: %s weight;pT (GeV);Inv. Mass (GeV)",filenameobj.weightnames[legendInt].c_str()));
    MultiGraphs->Draw("APE");



    
    // Draw the legend
    legend1->Draw();

    // Show the canvas
    MultiGraphs->GetXaxis()->SetLimits(filenameobj.plotxlims[0],filenameobj.plotxlims[1]);
    MultiGraphs->SetMinimum(filenameobj.plotylims[0]);
    MultiGraphs->SetMaximum(filenameobj.plotylims[1]);
    c1->SetMargin(0.2,0.1,0.1,0.1);
    gPad->Modified();
    gPad->Update();
    //canvas1->Modified();
    std::cout << "position 2" << std::endl; // debug line
    c1->Update();
    //c1->Print(Form("pioncode/canvas_pdf/%s_%f_OverlayMeanHistograms.pdf",filenameobj.filenamemod.c_str(),errparam));
    
    // Close the file and clean up
    file->Close();
    delete file;
    //delete legend1;
    std::cout << "position 3" << std::endl; // debug line
    return c1;
}


void ClusterOverlayTestFunc(filename_object filenameobj, const std::string& fileName, const char* histName, const std::string& LegendName){//only works if filenames.size()=2 !!
    TCanvas* canvas1 = new TCanvas("canvas1", "Overlay Means", 800, 600);
    TLegend* legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
    
    float errparam=filenameobj.sqrtEsmearing[0];

    ///*
    //need number of lines in hist. temp until I setup the object to hold that info?
    TFile* file = new TFile(filenameobj.fileNames[0].c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file " << fileName.c_str() << std::endl;
        //return nullptr;
    }

    // Get the 2D histogram from the file
    TH2F* hist2D = dynamic_cast<TH2F*>(file->Get("h27_2"));
    if (!hist2D) {
        std::cerr << "Error: Could not retrieve 2D histogram " << file->Get("h27_2") << " from file" << std::endl;
        file->Close();
        //return nullptr;
    }
    int NX= hist2D->GetNbinsX();
    file->Close();//*/
    //open the pdf?
    //int NX=128;
    canvas1->Print(Form("pioncode/canvas_pdf/%s_%f_InvMassprojections.pdf[",filenameobj.filenamemod.c_str(), errparam));

    for (int i=1;i<NX+1;i++){
        TH1D* yProjection1 = getYProjectionof2DHist(fileName.c_str(), "h27_2",i,i);//on
        TH1D* yProjection2 = getYProjectionof2DHist(fileName.c_str(), "h18_2",i,i);//off
        TH1D *histClone = (TH1D *)yProjection2->Clone("histClone");
        TH1D *ratioClone = (TH1D *)yProjection1->Clone("ratioClone");
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
        ratioClone->Divide(yProjection2);
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
    TH1D* yProjection1 = getYProjectionof2DHist(fileName.c_str(), "h27_2",1,NX);
    TH1D* yProjection2 = getYProjectionof2DHist(fileName.c_str(), "h18_2",1,NX);
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
    // Draw the legend
    //legend1->Draw();
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

