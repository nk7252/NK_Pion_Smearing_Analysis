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
Canvas* PlotAndFit(filename_object filenameobj,const std::string& fileName, const std::string& histName1, const std::string& histName2);
void GraphAndSaveToPDF(filename_object filenameobj, const std::string& histName1, const std::string& histName2);

int main () {
    // 0=weight type
    int fileset = 0 ;
    filename_object choosenfilenameobj = choosecomparisontype(fileset);
    
    OverlayMeans(choosenfilenameobj);
    OverlaySigmaOverMean(choosenfilenameobj);
    plotOverlayedHistograms(choosenfilenameobj, "h12");//h12 is smeared pion pT, Weighted. h3 is unsmeared pion pT, weighted
    //SliceAndFit(choosenfilenameobj);

    /*
    if (choosenfilenameobj.fileNames.size()==2){// for subtraction of inv mass profile
        ClusterOverlayTestFunc(choosenfilenameobj, "h18");//const char* histname
    }//*/
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
        filename_object1.fileNames={"pioncode/rootfiles/Pi0FastMC_0.155000_EXP.root", "pioncode/rootfiles/Pi0FastMC_0.155000_POWER.root", "pioncode/rootfiles/Pi0FastMC_0.155000_WSHP.root","pioncode/rootfiles/Pi0FastMC_0.155000_HAGEDORN.root"};
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

void GraphAndSaveToPDF(filename_object filenameobj, const std::string& histName1, const std::string& histName2) {
    // create canvas and open the output pdf
    TCanvas* canvas;    
    canvas->Print("output.pdf[");
    //loop over files. save 
    for (const auto& fileName : fileNames) {
        canvas = PlotAndFit(fileName, histName1, histName2);
        // Save canvas to PDF as a page
        canvas->Print("output.pdf");
    }

    // Close the PDF file
    canvas->Print("output.pdf]");
    //clean up
    delete canvas;
}

Canvas* PlotAndFit(filename_object filenameobj,const std::string& fileName, const std::string& histName1, const std::string& histName2) {
    TCanvas *canvas = new TCanvas("canvas", "Canvas", 800, 600);

    TFile *file = new TFile(fileName.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Unable to open file " << fileName << std::endl;
        delete canvas;
        return nullptr;
    }

    // Load histograms from the file
    TH1F *hist1 = static_cast<TH1F*>(file->Get(histName1.c_str()));
    TH1F *hist2 = static_cast<TH1F*>(file->Get(histName2.c_str()));

    if (!hist1 || !hist2) {
        std::cerr << "Error: Unable to retrieve histograms from file " << fileName << std::endl;
        file->Close();
        delete canvas;
        return nullptr;
    }
    //plotting operations.
    float errparam=filenameobj.sqrtEsmearing[0];
    double binres=2;//number of divisions per GeV
    // Create a legend
    TLegend* legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);//0.7, 0.4, 0.9, 0.6
    TMultiGraph *MultiGraphs = new TMultiGraph();//h18->GetNbinsX()
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

    // Close the file
    file->Close();

    return canvas;
}

