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
void OverlayMeans(filename_object filenameobj);
void OverlaySigmaOverMean(filename_object filenameobj);


int main () {
    // 0=weight type
    int fileset = 6 ;
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

void GraphAndSaveToPDF(const std::vector<std::string>& fileNames, const std::string& histName1, const std::string& histName2) {
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

Canvas* PlotAndFit(const std::string& fileName, const std::string& histName1, const std::string& histName2) {
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

    // Plot histograms on the same canvas
    hist1->Draw();
    hist2->Draw("SAME");

    // Optionally, add fitting or other operations

    // Close the file
    file->Close();

    return canvas;
}

