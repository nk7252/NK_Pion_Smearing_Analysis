#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <stdexcept>
#include <functional>
#include <utility>
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
#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

class filename_object
{ // object to hold my file names and related strings. things like settings for the files are contained here.
public:
    std::vector<std::string> fileNames;
    std::vector<std::string> legendnames;
    std::vector<std::string> weightnames;
    string filenamemod;
    std::string particletype;
    string canvasnamemod;
    std::vector<float> FitRange;
    std::vector<float> plotxlims;
    std::vector<float> plotylims;
    int pTcutoff;
    std::vector<int> sqrtEsmearing;
    int binres;
};
// declarations
int extractNumber(const std::string &filepath, int option);
filename_object choosecomparisontype(int choosetype);
TCanvas *FitMeanAndPlot(filename_object filenameobj, int legendInt, const std::string &fileName, std::vector<std::string> HistList, std::vector<std::string> HistLegend);
TCanvas *FitSigmaMeanAndPlot(filename_object filenameobj, int legendInt, const std::string &fileName, std::vector<std::string> HistList, std::vector<std::string> HistLegend);
void GraphAndSaveToPDF(filename_object filenameobj, std::vector<std::string> HistList, std::vector<std::string> HistLegend);
void ClusterOverlayTestFunc(filename_object filenameobj, const std::string &fileName, const char *histName1, const char *histName2, const std::string &LegendName);
TH1D *getYProjectionof2DHist(const char *fileName, const char *histName, int firstxbin, int lastxbin);
void transferHistogram(const char *sourceFileName, const char *histogramName, const char *targetFileName, const char *NewhistogramName);
void SliceAndFit(filename_object filenameobj, const char *histName, const char *fileName);
void ScaleHistogramErrorsAndFit(int EsmearfactorB, int EsmearfactorA, const char *fileName, const char *histName, double errorScaleFactor, double fitRangeLow, double fitRangeHigh, int numBins, double maxXRange, int histtype);
void ProcessTH3IntoGraphs(const std::string &fileName, const std::string &histName, int nSlices, const std::string &pdfName, int sliceSize = -1);
std::vector<TH2 *> SliceTH3(const std::string &fileName, const std::string &histName, int nSlices, int sliceSize = -1);
TGraphErrors *FitAndGenerateGraph(TH2 *slice, int index);
TCanvas *CompareMeans(filename_object comparisonFilenameObj, int legendInt, const filename_object &filenameobj1, const filename_object &filenameobj2, std::vector<std::string> HistList, std::vector<std::string> HistLegend);
void PlotHistogramsToPDF(std::vector<filename_object> filenameobjs, filename_object comparisonFilenameObj, std::vector<std::string> HistList, std::vector<std::string> HistLegend);
TCanvas *FitAndProjectY(filename_object filenameobj, const std::string &fileName, const std::string &histName);

void AIOFit()
{

    ///*
    // const char* destinationfile="pioncode/rootfiles/Pi0FastMC_0.154000_sqrte_0.130000_const.root";
    const char *sourcehistfile = "pioncode/rootfiles/Pi0FastMC_0.154000_sqrte_0.130000_const.root"; // using the source for functions later as source here

    // transferHistogram("pioncode/rootfiles/data/Data_from_Blair_original/diClusMass_23726_23746_nomPi0CalibCuts.root", "h_InvMass", sourcehistfile, "h_InvMass_data");

    // transferHistogram("pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_200_10000MeV-0000000013-00000.root", "h_InvMass_badcalib_smear_weighted_125", sourcehistfile, "h_InvMass_Single_pi0_smear12_5");

    // transferHistogram("pioncode/rootfiles/geant/nclus75/smear_125_cutson_noposcorr/OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_200_10000MeV-0000000013-00000.root", "h_InvMass_badcalib_smear_weighted_125", sourcehistfile, "h_InvMass_Single_pi0_weighted_nclus75_noposcor_smear125");
    // transferHistogram("pioncode/rootfiles/geant/nclus75/smear_115_cutson_noposcorr/OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_200_10000MeV-0000000013-00000.root", "h_InvMass_badcalib_smear_weighted_115", sourcehistfile, "h_InvMass_Single_pi0_weighted_nclus75_noposcor_smear115");
    // transferHistogram("pioncode/rootfiles/geant/nclus75/smear_072_cutson_noposcorr/OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_200_10000MeV-0000000013-00000.root", "h_InvMass_badcalib_smear_weighted_072", sourcehistfile, "h_InvMass_Single_pi0_weighted_nclus75_noposcor_smear072");
    // transferHistogram("pioncode/rootfiles/geant/nclus75/nosmear_noposcorr/OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_200_10000MeV-0000000013-00000_v2.root", "h_InvMass_weighted", sourcehistfile, "h_InvMass_Single_pi0_weighted_nclus75_noposcor_nosmear");

    // transferHistogram("pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_200_10000MeV-0000000013-00000.root", "h_pTdiff_InvMass", destinationhist, "h_pTdiff_InvMass_Single_pi0");

    // transferHistogram("pioncode/rootfiles/geant/smear_132_cutson_poscorr/OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_200_10000MeV-0000000013-00000.root", "h_InvMass_badcalib_smear_weighted_132", sourcehistfile, "h_InvMass_Single_pi0_weighted_poscor_smear132");//13% with poscor

    // transferHistogram("pioncode/rootfiles/geant/nosmear_poscorr/OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_200_10000MeV-0000000013-00000.root", "h_InvMass_weighted", sourcehistfile, "h_InvMass_Single_pi0_weighted_poscor_nosmear");//no smear with poscor

    // transferHistogram("pioncode/rootfiles/geant/nclus37/smear_125_cutson_noposcorr/OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_200_10000MeV-0000000013-00000.root", "h_InvMass_badcalib_smear_weighted_125", sourcehistfile, "h_InvMass_Single_pi0_weighted_poscor_smear125_nclus37");//12.5 smear with noposcor, nclus<37

    //*/

    // 0=weight type
    int fileset = 1;
    filename_object choosenfilenameobj = choosecomparisontype(fileset);
    std::vector<std::string> HistList = {"h29_", "h28_", "h29_v2_", "h28_v2_"};
    std::vector<std::string> HistLegend = {"Smear+Asymm", "SA+Cluster", "SA+Pos_Res", "SA+CL+PS"};
    // GraphAndSaveToPDF(choosenfilenameobj,  HistList, HistLegend);

    HistList = {"h30_", "h35_", "h31_", "h100_"};
    HistLegend = {"Cuts", "Cuts+Cluster", "Cuts+Pos_Res", "Cuts+CL+PR"};
    // std::vector<std::string> HistList={"h30_","h35_","h31_","h100_"};
    // std::vector<std::string> HistLegend={"Cuts","Cuts+Cluster","Cuts+Pos_Res","Cuts+CL+PR"};
    int fileset2 = 2; // eta
    filename_object choosenfilenameobj2 = choosecomparisontype(fileset2);
    GraphAndSaveToPDF(choosenfilenameobj2, HistList, HistLegend);

    int fileset0 = 0; // pion
    filename_object choosenfilenameobj0 = choosecomparisontype(fileset0);
    GraphAndSaveToPDF(choosenfilenameobj0, HistList, HistLegend);

    std::vector<filename_object> filenameobjs = {choosenfilenameobj0, choosenfilenameobj2};
    filename_object combfilenameobj = choosecomparisontype(3);
    // Assuming filenameobj1 and filenameobj2 are defined
    PlotHistogramsToPDF(filenameobjs, combfilenameobj, HistList, HistLegend);

    //*/
    // OverlayMeans(choosenfilenameobj);
    // OverlaySigmaOverMean(choosenfilenameobj);
    // plotOverlayedHistograms(choosenfilenameobj, "h12");//accepts 1d hist?//h12 is smeared pion pT, Weighted. h3 is unsmeared pion pT, weighted

    //*
    // SliceAndFit(choosenfilenameobj, "h18_2", sourcehistfile);// smeared
    // SliceAndFit(choosenfilenameobj, "h34_2", sourcehistfile);// position res

    // SliceAndFit(choosenfilenameobj, "h30_2", sourcehistfile);// just cuts
    // SliceAndFit(choosenfilenameobj, "h31_2", sourcehistfile);// cuts + pos res
    // SliceAndFit(choosenfilenameobj, "h35_2", sourcehistfile);// cuts + occupancy

    // SliceAndFit(choosenfilenameobj, "h100_2", sourcehistfile);// cuts + pos + occupancy
    //*/
    // void ScaleHistogramErrorsAndFit(int Esmearfactor ,const char* fileName, const char* histName,  double errorScaleFactor, double fitRangeLow, double fitRangeHigh, int numBins, double maxXRange, int histtype)
    // int histtype= 0=fastmc, 1=geant, 2=data?

    // ScaleHistogramErrorsAndFit(extractNumber(sourcehistfile, 1),extractNumber(sourcehistfile, 2), sourcehistfile, "h31_1d_2",  1.0, 0.13, 0.17 , 40, 0.4, 0);
    // ScaleHistogramErrorsAndFit(extractNumber(sourcehistfile, 1),extractNumber(sourcehistfile, 2), sourcehistfile, "h31_1d_2",  1.0, 0.13, 0.19 , 40, 0.4, 0);
    // ScaleHistogramErrorsAndFit(extractNumber(sourcehistfile, 1),extractNumber(sourcehistfile, 2), sourcehistfile, "h31_1d_2",  1.0, 0.10, 0.21 , 40, 0.4, 0);

    // ScaleHistogramErrorsAndFit(extractNumber(sourcehistfile, 1),extractNumber(sourcehistfile, 2), sourcehistfile, "h_InvMass_Single_pi0_weighted_noposcor_smear125",  1.0, 0.12, 0.18 , 40, 0.4, 0);

    //.09-.17 ~5
    // float start_point=0.09;
    // ScaleHistogramErrorsAndFit(extractNumber(sourcehistfile, 1),extractNumber(sourcehistfile, 2), sourcehistfile, "h_InvMass_Single_pi0_weighted_nclus75_noposcor_smear125",  1.0, 0.1, 0.18 , 40, 0.4, 1);

    // ScaleHistogramErrorsAndFit(extractNumber(sourcehistfile, 1),extractNumber(sourcehistfile, 2), sourcehistfile, "h_InvMass_Single_pi0_weighted_nclus75_noposcor_smear115",  1.0, 0.09, 0.19 , 40, 0.4, 1);
    // ScaleHistogramErrorsAndFit(extractNumber(sourcehistfile, 1),extractNumber(sourcehistfile, 2), sourcehistfile, "h_InvMass_Single_pi0_weighted_nclus75_noposcor_smear135",  1.0, 0.1, 0.20 , 40, 0.4, 1);
    // ScaleHistogramErrorsAndFit(extractNumber(sourcehistfile, 1),extractNumber(sourcehistfile, 2), sourcehistfile, "h_InvMass_Single_pi0_weighted_nclus75_noposcor_nosmear",  1.0, 0.13, 0.20 , 40, 0.4, 1);

    // ScaleHistogramErrorsAndFit(extractNumber(sourcehistfile, 1),extractNumber(sourcehistfile, 2), sourcehistfile, "h_InvMass_Single_pi0_weighted_poscor_nosmear",  1.0, 0.10, 0.17 , 40, 0.4, 1);

    // ScaleHistogramErrorsAndFit(extractNumber(sourcehistfile, 1),extractNumber(sourcehistfile, 2), sourcehistfile, "h_InvMass_Single_pi0_weighted_poscor_smear125_nclus37",  1.0, 0.10, 0.18 , 40, 0.4, 1);

    // ScaleHistogramErrorsAndFit(extractNumber(sourcehistfile, 1), extractNumber(sourcehistfile, 2), sourcehistfile, "h_InvMass_Single_pi0_smear12_5",  1.0, 0.10, 0.18 , 40, 0.4, 1);//if error for this and not the next one, you probably have a file where the smearing isn't applied.

    // ScaleHistogramErrorsAndFit(extractNumber(sourcehistfile, 1),extractNumber(sourcehistfile, 2), sourcehistfile, "h_InvMass_Single_pi0_weighted",  1.0, 0.10, 0.2 , 40, 0.4, 1);

    // ScaleHistogramErrorsAndFit(extractNumber(sourcehistfile, 1),extractNumber(sourcehistfile, 2), sourcehistfile, "h_InvMass_data",  1.0, 0.05, 0.3 , 40, 0.4, 2);

    // ProcessTH3IntoGraphs("pioncode/rootfiles/data/pt_nclus_differential_data/pt05pt05.root", "h_pipT_Nclus_mass", 20, "results");

    // ClusterOverlayTestFunc(choosenfilenameobj,"pioncode/rootfiles/Pi0FastMC_0.155000.root", "h27_2", "test");
    //  Code to exit ROOT after running the macro
    gApplication->Terminate(0);
}

int extractNumber(const std::string &filepath, int option)
{
    // file format
    //  pioncode/rootfiles/Pi0FastMC_0.154000_sqrte_0.130000_const.root
    // Validate the option (must be positive)
    if (option < 1)
    {
        std::cerr << "Error: Option must be a positive integer." << std::endl;
        return -1;
    }

    // Extract the filename from the filepath
    std::string filename = filepath.substr(filepath.find_last_of("/\\") + 1);

    // Split the filename at every period
    std::vector<std::string> segments;
    std::stringstream ss(filename);
    std::string item;
    while (std::getline(ss, item, '.'))
    {
        segments.push_back(item);
    }

    // Check if there are enough segments to meet the requested option
    // Since we're interested in segments after periods, the option should be less than the size of segments
    if (option >= segments.size())
    {
        std::cerr << "Error: The specified option exceeds the number of available numbers in the filepath." << std::endl;
        return -1;
    }

    // Attempt to convert the specified segment into a number
    try
    {
        int number = std::stoi(segments[option]); // Directly access the segment based on option
        return number;
    }
    catch (const std::invalid_argument &e)
    {
        std::cerr << "Error: Invalid number format in the specified segment." << std::endl;
        return -1; // or handle error differently
    }
    catch (const std::out_of_range &e)
    {
        std::cerr << "Error: Number out of range." << std::endl;
        return -1; // or handle error differently
    }
}

filename_object choosecomparisontype(int choosetype)
{
    filename_object filename_object1; // 0=weight type, 1=ac on/off, 2=co on/off, 3=ac&co on/off
    if (choosetype == 0)
    {
        // filename_object weightfilenameobj;
        filename_object1.fileNames = {"pioncode/rootfiles/PionFastMC_0.154000_sqrte_0.093000_const.root"};
        filename_object1.particletype = "Pion";
        filename_object1.legendnames = {"EXP", "POWER", "WSHP", "HAGEDORN"};
        filename_object1.weightnames = {"EXP", "POWER", "WSHP", "HAGEDORN"};
        filename_object1.filenamemod = "weightmethod_co1_ac1";
        // filename_object1.canvasnamemod=" for various weighting methods. asymm+clustering";
        filename_object1.FitRange = {0.1, 0.19};
        filename_object1.plotxlims = {2.01, 10};                          // min, max6.4
        filename_object1.plotylims = {0.14, 0.146, 0.08, 0.125, 0.0, 2.0}; // mean_min, mean_max,sm_min, sm_max, min h12, max h12
        filename_object1.pTcutoff = 30;
        filename_object1.binres = 2;
    }
    if (choosetype == 1)
    {
        // filename_object weightfilenameobj;
        filename_object1.fileNames = {"pioncode/rootfiles/EtaFastMC_0.154000_sqrte_0.093000_const.root"};
        filename_object1.particletype = "Eta";
        filename_object1.legendnames = {"EXP", "POWER", "WSHP", "HAGEDORN"};
        filename_object1.weightnames = {"EXP", "POWER", "WSHP", "HAGEDORN"};
        // filename_object1.legendnames={"WSHP"};
        // filename_object1.weightnames={"WSHP"};
        filename_object1.filenamemod = "weightmethod_co1_ac1";
        // filename_object1.canvasnamemod=" for various weighting methods. asymm+clustering";
        filename_object1.FitRange = {0.4, 0.8};
        filename_object1.plotxlims = {2.01, 10};                             // min, max6.4
        filename_object1.plotylims = {0.56, 0.59, 0.08, 0.16, 0.0, 2.0}; // mean_min, mean_max,sm_min, sm_max, min h12, max h12
        filename_object1.pTcutoff = 30;
        filename_object1.binres = 2;
    }
    if (choosetype == 2)
    {
        // filename_object weightfilenameobj;
        filename_object1.fileNames = {"pioncode/rootfiles/EtaFastMC_0.154000_sqrte_0.093000_const.root"};
        filename_object1.particletype = "Eta";
        filename_object1.legendnames = {"EXP", "POWER", "WSHP", "HAGEDORN"};
        filename_object1.weightnames = {"EXP", "POWER", "WSHP", "HAGEDORN"};
        // filename_object1.legendnames={"WSHP"};
        // filename_object1.weightnames={"WSHP"};
        filename_object1.filenamemod = "weightmethod_co1_ac1";
        // filename_object1.canvasnamemod=" for various weighting methods. asymm+clustering";
        filename_object1.FitRange = {0.4, 0.80};
        filename_object1.plotxlims = {2.01, 10};                             // min, max6.4
        filename_object1.plotylims = {0.568, 0.588, 0.075, 0.115, 0.0, 2.0}; // mean_min, mean_max,sm_min, sm_max, min h12, max h12
        filename_object1.pTcutoff = 30;
        filename_object1.binres = 2;
    }
    if (choosetype == 3)
    {
        // filename_object weightfilenameobj;
        filename_object1.fileNames = {"pioncode/rootfiles/EtaFastMC_0.154000_sqrte_0.093000_const.root"};
        filename_object1.particletype = "Pion&Eta";
        filename_object1.legendnames = {"EXP", "POWER", "WSHP", "HAGEDORN"};
        filename_object1.weightnames = {"EXP", "POWER", "WSHP", "HAGEDORN"};
        // filename_object1.legendnames={"WSHP"};
        // filename_object1.weightnames={"WSHP"};
        filename_object1.filenamemod = "weightmethod_co1_ac1";
        // filename_object1.canvasnamemod=" for various weighting methods. asymm+clustering";
        // filename_object1.FitRange={0.4,0.8};
        filename_object1.plotxlims = {2.01, 10};                            // min, max6.4
        filename_object1.plotylims = {0.23, 0.28, 0.08, 0.11, 0.0, 2.0}; // mean_min, mean_max,sm_min, sm_max, min h12, max h12
        filename_object1.pTcutoff = 30;
        filename_object1.binres = 2;
    }
    for (size_t i = 0; i < filename_object1.fileNames.size(); i++)
    {
        filename_object1.sqrtEsmearing.push_back(extractNumber(filename_object1.fileNames[i], 1));
    }
    return filename_object1;
}

void GraphAndSaveToPDF(filename_object filenameobj, std::vector<std::string> HistList, std::vector<std::string> HistLegend)
{
    // create canvas and open the output pdf
    TCanvas *canvas = new TCanvas("canvas1", "Canvas", 800, 600);
    // TCanvas* canvas;
    std::string combinedhiststring;
    for (const auto &str : HistList)
    {
        combinedhiststring += str;
    }
    canvas->Print(Form("pioncode/canvas_pdf/%s_%d_%s_OverlayPlot.pdf[", combinedhiststring.c_str(), filenameobj.sqrtEsmearing[0], filenameobj.particletype.c_str()));
    // loop over files. save
    int legendInt = 0;
    // loop over weight methods
    for (size_t l = 0; l < filenameobj.weightnames.size(); ++l)
    {
        std::vector<std::string> histogramName;
        // int legendInt=0;

        // Construct the histogram name with the index appended
        for (size_t v = 0; v < HistList.size(); ++v)
        {
            histogramName.push_back(HistList[v] + std::to_string(l));
            std::cout << histogramName[v] << std::endl; // debug line
        }

        for (const auto &fileName : filenameobj.fileNames)
        {
            canvas = FitMeanAndPlot(filenameobj, legendInt, fileName, histogramName, HistLegend);
            // Save mean canvas to PDF as a page
            canvas->Print(Form("pioncode/canvas_pdf/%s_%d_%s_OverlayPlot.pdf", combinedhiststring.c_str(), filenameobj.sqrtEsmearing[0], filenameobj.particletype.c_str()));

            canvas = FitSigmaMeanAndPlot(filenameobj, legendInt, fileName, histogramName, HistLegend);
            // Save sigma/mean canvas to PDF as a page
            canvas->Print(Form("pioncode/canvas_pdf/%s_%d_%s_OverlayPlot.pdf", combinedhiststring.c_str(), filenameobj.sqrtEsmearing[0], filenameobj.particletype.c_str()));

            // std::cout << "filename loop done" << std::endl; // debug line
        }

        // Clear the histogramName vector for the next iteration
        histogramName.clear();
        // std::cout << "l loop done" << std::endl; // debug line
        legendInt++;
    }

    // Close the PDF file
    canvas->Print(Form("pioncode/canvas_pdf/%s_%d_%s_OverlayPlot.pdf]", combinedhiststring.c_str(), filenameobj.sqrtEsmearing[0], filenameobj.particletype.c_str()));
    // clean up
    delete canvas;
}

void PlotHistogramsToPDF(std::vector<filename_object> filenameobjs, filename_object comparisonFilenameObj, std::vector<std::string> HistList, std::vector<std::string> HistLegend)
{
    // Create canvas and open the output PDF
    TCanvas *canvas = new TCanvas("canvas1", "Canvas", 800, 600);
    std::string combinedhiststring;
    for (const auto &str : HistList)
    {
        combinedhiststring += str;
    }
    std::string pdfFileName = Form("pioncode/canvas_pdf/%s_%d_%s_OverlayPlot.pdf", combinedhiststring.c_str(), comparisonFilenameObj.sqrtEsmearing[0], comparisonFilenameObj.particletype.c_str());
    canvas->Print((pdfFileName + "[").c_str());

    // Loop over weight methods
    for (size_t l = 0; l < comparisonFilenameObj.weightnames.size(); ++l)
    {
        std::vector<std::string> histogramName;
        for (size_t v = 0; v < HistList.size(); ++v)
        {
            histogramName.push_back(HistList[v] + std::to_string(l));
            // std::cout << histogramName[v] << std::endl; // Debug line
        }

        const auto &filenameobj1 = filenameobjs[0];
        const auto &filenameobj2 = filenameobjs[1];

        TCanvas *meanCanvas = CompareMeans(comparisonFilenameObj, l, filenameobj1, filenameobj2, histogramName, HistLegend);
        if (meanCanvas)
        {
            meanCanvas->Print(pdfFileName.c_str());
            delete meanCanvas;
        }
        /*
    //for (const auto& histName : histogramName) {
        TCanvas* projectionCanvas1 = FitAndProjectY(filenameobj1, filenameobj1.fileNames[0], histogramName[l]);
        if (projectionCanvas1) {
            projectionCanvas1->Print(pdfFileName.c_str());
            delete projectionCanvas1;
        }

        TCanvas* projectionCanvas2 = FitAndProjectY(filenameobj2, filenameobj2.fileNames[0], histogramName[l]);
        if (projectionCanvas2) {
            projectionCanvas2->Print(pdfFileName.c_str());
            delete projectionCanvas2;
        }

        TCanvas* sigmaCanvas1 = FitSigmaMeanAndPlot(filenameobj1, l, filenameobj1.fileNames[0], histogramName[l], HistLegend[l]);
        if (sigmaCanvas1) {
            sigmaCanvas1->Print(pdfFileName.c_str());
            delete sigmaCanvas1;
        }


        TCanvas* sigmaCanvas2 = FitSigmaMeanAndPlot(filenameobj2, l, filenameobj2.fileNames[0], histogramName[l], HistLegend[l]);
        if (sigmaCanvas2) {
            sigmaCanvas2->Print(pdfFileName.c_str());
            delete sigmaCanvas2;
        }

*/

        // std::cout << "Finished processing files: " << filenameobj1.fileNames[0] << " and " << filenameobj2.fileNames[0] << std::endl; // Debug line
    }

    // histogramName.clear();
    // std::cout << "Finished processing weight method: " << l << std::endl; // Debug line
    //  Close the PDF file
    canvas->Print((pdfFileName + "]").c_str());
    delete canvas;
}

TCanvas *FitAndProjectY(filename_object filenameobj, const std::string &fileName, const std::string &histName)
{
    // Open the file
    TFile *file = TFile::Open(fileName.c_str());
    if (!file || file->IsZombie())
    {
        std::cerr << "Error opening file: " << fileName << std::endl;
        return nullptr;
    }

    // Retrieve the 2D histogram
    TH2F *hist2D = dynamic_cast<TH2F *>(file->Get(histName.c_str()));
    if (!hist2D)
    {
        std::cerr << "Error retrieving histogram: " << histName << " from file: " << fileName << std::endl;
        file->Close();
        return nullptr;
    }

    // Project the histogram onto the Y axis
    TH1D *projY = hist2D->ProjectionY("projY", 0, hist2D->GetXaxis()->GetNbins());
    projY->SetTitle(Form("Projection of %s on Y-axis", histName.c_str()));

    // Create a canvas to draw the projection
    TCanvas *c2 = new TCanvas(Form("canvas_%s", histName.c_str()), "Projection Canvas", 800, 600);
    projY->Draw();
    c2->Update(); // Ensure the canvas is updated

    // Close the file
    file->Close();

    return c2;
}

TCanvas *FitMeanAndPlot(filename_object filenameobj, int legendInt, const std::string &fileName, std::vector<std::string> HistList, std::vector<std::string> HistLegend)
{
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
    // create canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    // Create a legend
    TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);

    TFile *file = new TFile(fileName.c_str(), "READ");
    if (!file || file->IsZombie())
    {
        std::cerr << "Error: Unable to open file " << fileName << std::endl;
        delete c1;
        delete legend1;
        return nullptr;
    }

    // plotting operations.
    float errparam = filenameobj.sqrtEsmearing[0];
    int binres = filenameobj.binres; // number of divisions per GeV

    // int numhists=HistList.size();
    //  Declare Histograms
    std::vector<TH2F *> temphist(HistList.size());
    // TMultiGraph to hold graphs.
    TMultiGraph *MultiGraphs = new TMultiGraph();
    // Create Tgraphs to hold means for each histogram
    std::vector<TGraphErrors *> meanGraph(HistList.size());

    // loop over the two histograms being compared. probably a better way to do this since it is a 1d array.
    for (size_t j = 0; j < HistList.size(); j++)
    {
        std::cout << "file" << j << " of" << HistList.size() << std::endl; // debug line
        // Load histograms from the file
        temphist[j] = dynamic_cast<TH2F *>(file->Get(HistList[j].c_str()));
        meanGraph[j] = new TGraphErrors(temphist[j]->GetNbinsX());
        // std::cout << "loaded hists? made graph" << std::endl; // debug line
        if (!temphist[j])
        {
            std::cerr << "Error: Unable to retrieve histogram " << j << " from file " << fileName << std::endl;
            file->Close();
            delete c1;
            return nullptr;
        }
        // std::cout << "pre bin loop" << std::endl; // debug line
        //  Loop over each bin in the X direction
        for (int binX = 1; binX <= temphist[j]->GetNbinsX(); binX++)
        {

            // Project along Y for each binX
            TH1D *yProjection = temphist[j]->ProjectionY(Form("YProjection_%zu_%d", j, binX), binX, binX, "");

            TF1 *fitFunc = new TF1("gaus", "gaus", filenameobj.FitRange[0], filenameobj.FitRange[1]);

            // Fit the Y projection with a Gaussian
            yProjection->Fit("gaus", "QER");

            // Access the fit parameters
            // TF1* fitFunc = yProjection->GetFunction("gaus");

            // Check if the fit function is valid
            if (fitFunc)
            {
                //meanGraph[j]->SetPoint(binX, binX / binres, fitFunc->GetParameter(1));
                //meanGraph[j]->SetPointError(binX, 0, fitFunc->GetParError(1));
                meanGraph[j]->AddPointError (binX / binres, fitFunc->GetParameter(1), 0, fitFunc->GetParError(1));
                // Set the legend entry to the title of the original 2D histogram
                // meanGraph[j]->SetTitle(temphist[j]->GetTitle());
            }
            // Clean up Y projection
            delete fitFunc;
            delete yProjection;
        }
        // std::cout << "bin loop done" << std::endl; // debug line
        //  Set different line colors for each version
        int MarkerStyle = j + 24; //
        int MarkerColor = j + 1;
        // meanGraph->SetLineColor(lineColor);
        meanGraph[j]->SetMarkerStyle(MarkerStyle);
        meanGraph[j]->SetMarkerColor(MarkerColor);
        // Overlay the mean histogram on the same canvas
        if (j == 0)
        {
            // meanGraph->Draw("AP"); // Draw histogram for the first version
            // canvas1->Print("OverlayMeanHistograms.pdf");
            MultiGraphs->Add(meanGraph[j], "PE");
            // std::cout << "draw for the first file" << std::endl; // debug line
        }
        else
        {
            MultiGraphs->Add(meanGraph[j], "PE");
            // meanGraph->Draw("P SAME"); // Draw subsequent histograms on the same canvas
            // std::cout << "draw for subsequent" << std::endl; // debug line
        }
        // Add an entry to the legend
        legend1->AddEntry(meanGraph[j], HistLegend[j].c_str(), "P"); // meanGraph[j]->GetTitle(),"P"
    }

    std::cout << "position 1" << std::endl; // debug line
    MultiGraphs->SetTitle(Form("Smeared %s pT vs Inv Mass: %s weight;pT (GeV);Inv. Mass (GeV)", filenameobj.particletype.c_str(), filenameobj.weightnames[legendInt].c_str()));
    MultiGraphs->Draw("APE");

    // Draw the legend
    legend1->Draw();

    // Show the canvas
    MultiGraphs->GetXaxis()->SetLimits(filenameobj.plotxlims[0], filenameobj.plotxlims[1]);
    MultiGraphs->SetMinimum(filenameobj.plotylims[0]);
    MultiGraphs->SetMaximum(filenameobj.plotylims[1]);
    c1->SetMargin(0.2, 0.1, 0.1, 0.1);
    gPad->Modified();
    gPad->Update();
    // canvas1->Modified();
    // std::cout << "position 2" << std::endl; // debug line
    c1->Update();
    // c1->Print(Form("pioncode/canvas_pdf/%s_%f_OverlayMeanHistograms.pdf",filenameobj.filenamemod.c_str(),errparam));

    // Close the file and clean up
    file->Close();
    delete file;
    // delete legend1;
    // std::cout << "position 3" << std::endl; // debug line
    return c1;
}

TCanvas *FitSigmaMeanAndPlot(filename_object filenameobj, int legendInt, const std::string &fileName, std::vector<std::string> HistList, std::vector<std::string> HistLegend)
{
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
    // create canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    // Create a legend
    TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);

    TFile *file = new TFile(fileName.c_str(), "READ");
    if (!file || file->IsZombie())
    {
        std::cerr << "Error: Unable to open file " << fileName << std::endl;
        delete c1;
        delete legend1;
        return nullptr;
    }

    // plotting operations.
    float errparam = filenameobj.sqrtEsmearing[0];
    int binres = filenameobj.binres; // number of divisions per GeV

    // int numhists=HistList.size();
    //  Declare Histograms
    std::vector<TH2F *> temphist(HistList.size());
    // TMultiGraph to hold graphs.
    TMultiGraph *MultiGraphs = new TMultiGraph();
    // Create Tgraphs to hold means for each histogram
    std::vector<TGraphErrors *> meanGraph(HistList.size());

    // loop over the two histograms being compared. probably a better way to do this since it is a 1d array.
    for (size_t j = 0; j < HistList.size(); j++)
    {
        std::cout << "file" << j << " of" << HistList.size() << std::endl; // debug line
        // Load histograms from the file
        temphist[j] = dynamic_cast<TH2F *>(file->Get(HistList[j].c_str()));
        meanGraph[j] = new TGraphErrors(temphist[j]->GetNbinsX());
        // std::cout << "loaded hists? made graph" << std::endl; // debug line
        if (!temphist[j])
        {
            std::cerr << "Error: Unable to retrieve histogram " << j << " from file " << fileName << std::endl;
            file->Close();
            delete c1;
            return nullptr;
        }
        // std::cout << "pre bin loop" << std::endl; // debug line
        //  Loop over each bin in the X direction
        for (int binX = 1; binX <= temphist[j]->GetNbinsX(); binX++)
        {

            // Project along Y for each binX
            TH1D *yProjection = temphist[j]->ProjectionY(Form("YProjection_%zu_%d", j, binX), binX, binX, "");
            TF1 *fitFunc = new TF1("gaus", "gaus", filenameobj.FitRange[0], filenameobj.FitRange[1]);
            // Fit the Y projection with a Gaussian
            yProjection->Fit("gaus", "QER");

            // Access the fit parameters
            // TF1* fitFunc = yProjection->GetFunction("gaus");

            // Check if the fit function is valid
            if (fitFunc)
            {
                // double meanoversigma =fitFunc->GetParameter(2)/fitFunc->GetParameter(1);
                // double meanoversigmaerr = meanoversigma*(fitFunc->GetParError(2)/fitFunc->GetParameter(2)+fitFunc->GetParError(1)/fitFunc->GetParameter(1));//m/s*(serr/s+merr/m)

                //meanGraph[j]->SetPoint(binX, binX / binres, fitFunc->GetParameter(2) / fitFunc->GetParameter(1));
                //meanGraph[j]->SetPointError(binX, 0, (fitFunc->GetParameter(2) / fitFunc->GetParameter(1)) * (fitFunc->GetParError(2) / fitFunc->GetParameter(2) + fitFunc->GetParError(1) / fitFunc->GetParameter(1))); //(sigma/mean)*[sigmaerr/sigma+meanerr/mean]
                meanGraph[j]->AddPointError(binX / binres, fitFunc->GetParameter(1), 0, fitFunc->GetParError(1));
                // Set the legend entry to the title of the original 2D histogram
                // meanGraph[j]->SetTitle(temphist[j]->GetTitle());
            }
            // Clean up Y projection
            delete fitFunc;
            delete yProjection;
        }
        // std::cout << "bin loop done" << std::endl; // debug line
        //  Set different line colors for each version
        int MarkerStyle = j + 24; //
        int MarkerColor = j + 1;
        // meanGraph->SetLineColor(lineColor);
        meanGraph[j]->SetMarkerStyle(MarkerStyle);
        meanGraph[j]->SetMarkerColor(MarkerColor);
        // Overlay the mean histogram on the same canvas
        if (j == 0)
        {
            // meanGraph->Draw("AP"); // Draw histogram for the first version
            // canvas1->Print("OverlayMeanHistograms.pdf");
            MultiGraphs->Add(meanGraph[j], "PE");
            // std::cout << "draw for the first file" << std::endl; // debug line
        }
        else
        {
            MultiGraphs->Add(meanGraph[j], "PE");
            // meanGraph->Draw("P SAME"); // Draw subsequent histograms on the same canvas
            // std::cout << "draw for subsequent" << std::endl; // debug line
        }
        // Add an entry to the legend
        legend1->AddEntry(meanGraph[j], HistLegend[j].c_str(), "P"); // meanGraph[j]->GetTitle(),"P"
    }

    // std::cout << "position 1" << std::endl; // debug line
    MultiGraphs->SetTitle(Form("Smeared %s pT vs Sigma/Inv Mass Mean: %s weight;pT (GeV);Sigma/Inv. Mass ", filenameobj.particletype.c_str(), filenameobj.weightnames[legendInt].c_str()));
    MultiGraphs->Draw("APE");

    // Draw the legend
    legend1->Draw();

    // Show the canvas
    MultiGraphs->GetXaxis()->SetLimits(filenameobj.plotxlims[0], filenameobj.plotxlims[1]);
    MultiGraphs->SetMinimum(filenameobj.plotylims[2]);
    MultiGraphs->SetMaximum(filenameobj.plotylims[3]);
    c1->SetMargin(0.2, 0.1, 0.1, 0.1);
    gPad->Modified();
    gPad->Update();
    // canvas1->Modified();
    // std::cout << "position 2" << std::endl; // debug line
    c1->Update();
    // c1->Print(Form("pioncode/canvas_pdf/%s_%f_OverlayMeanHistograms.pdf",filenameobj.filenamemod.c_str(),errparam));

    // Close the file and clean up
    file->Close();
    delete file;
    // delete legend1;
    // std::cout << "position 3" << std::endl; // debug line
    return c1;
}

TCanvas *CompareMeans(filename_object comparisonFilenameObj, int legendInt, const filename_object &filenameobj1, const filename_object &filenameobj2, std::vector<std::string> HistList, std::vector<std::string> HistLegend)
{
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
    // Create canvas
    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    // Create a legend
    TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);

    TFile *file1 = new TFile(filenameobj1.fileNames[0].c_str(), "READ");
    TFile *file2 = new TFile(filenameobj2.fileNames[0].c_str(), "READ");
    if (!file1 || file1->IsZombie() || !file2 || file2->IsZombie())
    {
        std::cerr << "Error: Unable to open files " << filenameobj1.fileNames[0] << " or " << filenameobj2.fileNames[0] << std::endl;
        delete c1;
        delete legend1;
        return nullptr;
    }

    int binres = comparisonFilenameObj.binres; // Number of divisions per GeV
    std::vector<TH2F *> temphist1(HistList.size()), temphist2(HistList.size());
    TMultiGraph *MultiGraphs = new TMultiGraph();
    std::vector<TGraphErrors *> meanGraph(HistList.size());

    for (size_t j = 0; j < HistList.size(); j++)
    {
        temphist1[j] = dynamic_cast<TH2F *>(file1->Get(HistList[j].c_str()));
        temphist2[j] = dynamic_cast<TH2F *>(file2->Get(HistList[j].c_str()));
        meanGraph[j] = new TGraphErrors(temphist1[j]->GetNbinsX());

        if (!temphist1[j] || !temphist2[j])
        {
            std::cerr << "Error: Unable to retrieve histogram " << j << " from files " << filenameobj1.fileNames[0] << " or " << filenameobj2.fileNames[0] << std::endl;
            file1->Close();
            file2->Close();
            delete c1;
            return nullptr;
        }

        for (int binX = 1; binX <= temphist1[j]->GetNbinsX(); binX++)
        {
            TH1D *yProjection1 = temphist1[j]->ProjectionY(Form("YProjection1_%zu_%d", j, binX), binX, binX, "");
            TH1D *yProjection2 = temphist2[j]->ProjectionY(Form("YProjection2_%zu_%d", j, binX), binX, binX, "");

            yProjection1->Fit("gaus", "QE");
            yProjection2->Fit("gaus", "QE");

            TF1 *fitFunc1 = yProjection1->GetFunction("gaus");
            TF1 *fitFunc2 = yProjection2->GetFunction("gaus");

            if (fitFunc1 && fitFunc2)
            {
                double mu1 = fitFunc1->GetParameter(1);
                double mu2 = fitFunc2->GetParameter(1);
                double mu1Err = fitFunc1->GetParError(1);
                double mu2Err = fitFunc2->GetParError(1);

                double ratio = mu1 / mu2;
                double ratioErr = ratio * sqrt(pow(mu1Err / mu1, 2) + pow(mu2Err / mu2, 2));

                meanGraph[j]->SetPoint(binX, binX / binres, ratio);
                meanGraph[j]->SetPointError(binX, 0, ratioErr);
            }

            delete fitFunc1;
            delete fitFunc2;
            delete yProjection1;
            delete yProjection2;
        }

        int MarkerStyle = j + 24;
        int MarkerColor = j + 1;
        meanGraph[j]->SetMarkerStyle(MarkerStyle);
        meanGraph[j]->SetMarkerColor(MarkerColor);

        MultiGraphs->Add(meanGraph[j], "PE");
        legend1->AddEntry(meanGraph[j], HistLegend[j].c_str(), "P");
    }

    MultiGraphs->SetTitle(Form("Smeared %s pT vs Mean Ratio: %s weight;pT (GeV);Mean Ratio", comparisonFilenameObj.particletype.c_str(), comparisonFilenameObj.weightnames[legendInt].c_str()));
    MultiGraphs->Draw("APE");

    legend1->Draw();

    MultiGraphs->GetXaxis()->SetLimits(comparisonFilenameObj.plotxlims[0], comparisonFilenameObj.plotxlims[1]);
    MultiGraphs->SetMinimum(comparisonFilenameObj.plotylims[0]);
    MultiGraphs->SetMaximum(comparisonFilenameObj.plotylims[1]);
    c1->SetMargin(0.2, 0.1, 0.1, 0.1);
    gPad->Modified();
    gPad->Update();

    file1->Close();
    file2->Close();
    delete file1;
    delete file2;

    return c1;
}

void ClusterOverlayTestFunc(filename_object filenameobj, const std::string &fileName, const char *histName1, const char *histName2, const std::string &LegendName)
{ // only works if filenames.size()=2 !!
    TCanvas *canvas1 = new TCanvas("canvas1", "Overlay Means", 800, 600);
    TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);

    float errparam = filenameobj.sqrtEsmearing[0];

    ///*
    // need number of lines in hist. temp until I setup the object to hold that info?
    TFile *file = new TFile(filenameobj.fileNames[0].c_str(), "READ");
    if (!file || file->IsZombie())
    {
        std::cerr << "Error: Could not open file " << fileName.c_str() << std::endl;
        // return nullptr;
    }

    // Get the 2D histogram from the file
    TH2F *hist2D = dynamic_cast<TH2F *>(file->Get(histName2));
    if (!hist2D)
    {
        std::cerr << "Error: Could not retrieve 2D histogram " << file->Get(histName2) << " from file" << std::endl;
        file->Close();
        // return nullptr;
    }
    int NX = hist2D->GetNbinsX();
    file->Close(); //*/
    // open the pdf?
    // int NX=128;
    canvas1->Print(Form("pioncode/canvas_pdf/%s_%f_InvMassprojections.pdf[", filenameobj.filenamemod.c_str(), errparam));

    for (int i = 1; i < NX + 1; i++)
    {
        TH1D *yProjection1 = getYProjectionof2DHist(fileName.c_str(), "histName2_2", i, i); // on
        TH1D *yProjection2 = getYProjectionof2DHist(fileName.c_str(), "histName1_2", i, i); // off
        TH1D *histClone = (TH1D *)yProjection2->Clone("histClone");
        TH1D *ratioClone = (TH1D *)yProjection1->Clone("ratioClone");
        histClone->Add(yProjection1, -1);
        canvas1->Divide(1, 3);
        canvas1->cd(1);
        yProjection1->Draw();
        yProjection1->SetLineColor(kRed);
        yProjection2->Draw("SAME");
        yProjection2->SetLineColor(kBlue);
        yProjection1->SetTitle(Form("Cluster Overlay on vs off. Bin %i;Invariant Mass (GeV); Counts", i));
        gPad->Modified();
        gPad->Update();
        canvas1->cd(2);
        histClone->Draw();
        histClone->SetTitle(Form("Cluster Overlay on - off. Bin %i;Invariant Mass (GeV); Counts", i));
        canvas1->cd(3);
        ratioClone->Divide(yProjection2);
        ratioClone->Draw();
        ratioClone->SetMaximum(5);
        ratioClone->SetTitle(Form("Cluster Overlay on/off. Bin %i;Invariant Mass (GeV); Counts", i));
        canvas1->Modified();
        legend1->AddEntry(yProjection1, filenameobj.legendnames[0].c_str(), "P");
        legend1->AddEntry(yProjection2, filenameobj.legendnames[1].c_str(), "P");
        canvas1->Print(Form("pioncode/canvas_pdf/%s_%f_InvMassprojections.pdf", filenameobj.filenamemod.c_str(), errparam)); // Print-> works too
        canvas1->Clear();
        legend1->Clear();
    }
    TH1D *yProjection1 = getYProjectionof2DHist(fileName.c_str(), "histName2_2", 1, NX);
    TH1D *yProjection2 = getYProjectionof2DHist(fileName.c_str(), "histName1_2", 1, NX);
    TH1D *histClone = (TH1D *)yProjection2->Clone("histClone");
    TH1D *ratioClone = (TH1D *)yProjection2->Clone("ratioClone");
    histClone->Add(yProjection1, -1);
    canvas1->Divide(1, 3);
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
    // legend1->Draw();
    canvas1->Print(Form("pioncode/canvas_pdf/%s_%f_InvMassprojections.pdf", filenameobj.filenamemod.c_str(), errparam)); // Print-> works too
    // close the pdf?
    canvas1->Print(Form("pioncode/canvas_pdf/%s_%f_InvMassprojections.pdf]", filenameobj.filenamemod.c_str(), errparam));
    // Clean up
    delete canvas1;
    // delete legend1;
}

TH1D *getYProjectionof2DHist(const char *fileName, const char *histName, int firstxbin, int lastxbin)
{
    // Open the root file
    TFile *file = new TFile(fileName, "READ");
    if (!file || file->IsZombie())
    {
        std::cerr << "Error: Could not open file " << fileName << std::endl;
        return nullptr;
    }

    // Get the 2D histogram from the file
    TH2F *hist2D = dynamic_cast<TH2F *>(file->Get(histName));
    if (!hist2D)
    {
        std::cerr << "Error: Could not retrieve 2D histogram " << histName << " from file" << std::endl;
        file->Close();
        return nullptr;
    }
    // int NX = hist2D->GetNbinsX();
    //  Create a 1D histogram for the y projection
    TH1D *hist1D = hist2D->ProjectionY("_py", firstxbin, lastxbin, "");
    // TH1D* hist1D = hist2D->ProjectionY("_py",1, 20,"");
    hist1D->SetDirectory(0);
    // Close the file
    file->Close();

    return hist1D;
}

void ScaleHistogramErrorsAndFit(int EsmearfactorB, int EsmearfactorA, const char *fileName, const char *histName, double errorScaleFactor, double fitRangeLow, double fitRangeHigh, int numBins, double maxXRange, int histtype)
{
    // filename_object filenameobj,
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
    SetsPhenixStyle();

    cout << "\n"
         << "processing: " << fileName << " Histogram: " << histName << "\n"
         << " With Error scaled by: " << errorScaleFactor << "\n"
         << " Fit from " << fitRangeLow << " to " << fitRangeHigh << "\n";
    TFile *pionfile = new TFile(fileName, "READ");
    if (!pionfile || pionfile->IsZombie())
    {
        std::cerr << "Error opening file." << std::endl;
        // Handle error or return
    }
    // TH2F* hist2D = dynamic_cast<TH2F*>(pionfile->Get(histName));
    TH1F *hist1D = (TH1F *)pionfile->Get(histName);
    if (!hist1D)
    {
        std::cerr << "Histogram not found." << std::endl;
        // Handle error or return
    }
    // Rebin the histogram to have 'numBins' bins
    // First, calculate the rebin factor assuming the histogram's range is 0 to maxXRange
    int currentNumBins = hist1D->GetNbinsX();
    double currentXMax = hist1D->GetXaxis()->GetXmax();
    int rebinFactor = currentNumBins / numBins;
    if (rebinFactor > 1)
    { // Only rebin if the factor is greater than 1
        std::cout << "current nbins: " << currentNumBins << " requested nbins: " << numBins << " rebin by: " << rebinFactor << std::endl;
        hist1D->Rebin(rebinFactor);
        std::cout << "new nbin check: " << hist1D->GetNbinsX() << std::endl;
    }

    // Set the maximum range on the x-axis to maxXRange
    // Note: This should be done after rebinning to maintain consistent bin widths
    hist1D->GetXaxis()->SetRangeUser(hist1D->GetXaxis()->GetXmin(), maxXRange);

    // Scale the error bars
    if (errorScaleFactor > 1)
    {
        std::cout << "rescaling error bars by: " << errorScaleFactor << std::endl;
        for (int i = 1; i <= hist1D->GetNbinsX(); ++i)
        {
            double originalError = hist1D->GetBinError(i);
            hist1D->SetBinError(i, originalError * errorScaleFactor);
        }
    }

    // Fit a Gaussian to the histogram within the specified range
    TF1 *gaussFit = new TF1("gaussFit", "gaus", fitRangeLow, fitRangeHigh);
    // hist1D->Fit(gaussFit, "WRQM"); // initial fit paramaters
    hist1D->Fit(gaussFit, "R"); // ML what does M do? Improve algorithm of tminuit. E is error est by Minos technique, L log likelihood

    // Create a canvas to draw the histogram and fit
    TString canvasname = Form("%s_PeakFit_SmearB_%d_A_%d_Thousandths_escale_%f", histName, EsmearfactorB / 1000, EsmearfactorA / 1000, errorScaleFactor);
    TCanvas *c1 = new TCanvas(canvasname, canvasname, 800, 600);
    hist1D->SetMinimum(0.0);
    hist1D->SetTitle(Form("%s;Inv. Mass (GeV); Counts", histName));
    hist1D->GetYaxis()->SetTitleOffset(1.5);
    float xbleft = 0.42;
    float ybleft = 0.7;
    float xtright = 0.9;
    float ytright = 0.9;
    TLegend *leg = new TLegend(xbleft, ybleft, xtright, ytright); // SetX1, SetY1, SetX2, SetY2; 1:bottom-left corner of legend, 2: top-right corner
    // The coordinates are in the range from 0 to 1, where (0, 0) is the bottom-left corner and (1, 1) is the top-right corner of the pad.
    leg->SetFillStyle(0);

    leg->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
    hist1D->SetStats(0); // Turn off stat box
    if (histtype == 1)
    { // SPMC=
        // leg->AddEntry("", "#it{#bf{sPHENIX}} Simulation, Internal", "");
        leg->AddEntry("", "Geant4", "");
        leg->AddEntry("", "Single Particle Simulation", "");
        // gStyle->SetOptFit(112); // 1 (fit parameters) + 10 (errors) + 2 (Chi-square/ndf)
        // gStyle->SetOptStat(0); // Turn off default statistics
    }

    // leg->AddEntry("", "Work In Progress", "");
    hist1D->Draw("E");      // Draw histogram with error bars
    gaussFit->Draw("SAME"); // Draw the fit on the same canvas
    leg->Draw("SAME");
    gPad->Modified();
    gPad->Update();

    TPaveText *pt2 = new TPaveText(xbleft + .2, 0.5, xtright, 0.7, "NDC"); // Adjust coordinates as needed
    pt2->SetFillColor(0);                                                  // Set the fill color to 0 for transparency
    pt2->SetFillStyle(0);                                                  // Set fill style to 0 (solid) with color 0 for transparency
    pt2->AddText(Form("#chi^{2}/NDF = %.2f", gaussFit->GetChisquare() / gaussFit->GetNDF()));
    pt2->AddText(Form("Mean = %.4f", gaussFit->GetParameter(1)));
    pt2->AddText(Form("Sigma = %.4f", gaussFit->GetParameter(2)));
    pt2->AddText(Form("Relative Width: %.2f%%", gaussFit->GetParameter(2) * 100.0f / gaussFit->GetParameter(1)));
    pt2->Draw("SAME");
    gPad->Modified(); // Apply the changes to the pad

    std::cout << "Chi-squared: " << gaussFit->GetChisquare() << std::endl;
    std::cout << "Number of Degrees of Freedom: " << gaussFit->GetNDF() << std::endl;
    std::cout << "Chi-squared/NDF: " << gaussFit->GetChisquare() / gaussFit->GetNDF() << std::endl;
    std::cout << "Relative Width: " << gaussFit->GetParameter(2) * 100.0f / gaussFit->GetParameter(1) << std::endl;

    // const char *pdfname = canvasname;
    //  Print the canvas to a file
    c1->SaveAs(Form("pioncode/canvas_pdf/%s.pdf", canvasname.Data()));

    // Second canvas: Custom list of fit results
    TCanvas *c2 = new TCanvas("canvas2", "Fit Parameters", 800, 600);
    c2->cd();

    TPaveText *pt = new TPaveText(0.1, 0.1, 0.9, 0.9, "blNDC"); // blNDC: borderless, normalized coordinates
    pt->SetTextAlign(12);                                       // Align text to the left
    pt->SetFillColor(0);                                        // Transparent background

    // Adding custom text entries
    pt->AddText(Form("Hist = %s", histName));
    pt->AddText("Fit Parameters:");
    pt->AddText(Form("Fit Range = %f to %f", fitRangeLow, fitRangeHigh));
    pt->AddText(Form("Mean = %f +/- %f", gaussFit->GetParameter(1), gaussFit->GetParError(1)));
    pt->AddText(Form("Sigma = %f +/- %f", gaussFit->GetParameter(2), gaussFit->GetParError(2)));
    pt->AddText(Form("Relative Width: %f", gaussFit->GetParameter(2) * 100.0f / gaussFit->GetParameter(1)));
    pt->AddText(Form("Chi2/NDF = %f / %d= %f", gaussFit->GetChisquare(), gaussFit->GetNDF(), gaussFit->GetChisquare() / gaussFit->GetNDF()));
    if (histtype == 0)
    { // fastmc
        pt->AddText("smearing details:");
        pt->AddText(Form("E smear B/sqrt(E) = %.1f percent", EsmearfactorB / 10000.f));
        pt->AddText(Form("E smear A = %.1f percent", EsmearfactorA / 10000.f));
        pt->AddText("Pos res = 2.8 mm ");
        pt->AddText("Cluster prob=off for h31,\n 1 percent for h100");
    }
    if (histtype == 0 || histtype == 1)
        pt->AddText("Pions Generated from 0.2-10 GeV"); // fastmc and spmc

    // You can add more custom lines as needed
    // pt->AddText("Additional Info: ...");

    pt->Draw();
    c2->SaveAs(Form("pioncode/canvas_pdf/FitInfo_%s.pdf", canvasname.Data()));

    // Clean up
    delete c1; // Automatically deletes gaussFit as well since it's drawn on the canvas
    delete c2; // Automatically deletes gaussFit as well since it's drawn on the canvas
}

// 2d hist fitting and metrics
void SliceAndFit(filename_object filenameobj, const char *histName, const char *fileName)
{
    // ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
    // int E_error_param = 155;static_cast<int>(
    // int E_error_param = filenameobj.sqrtEsmearing[0];
    // std::vector<std::string> WeightNames = filenameobj.weightnames//{"EXP", "POWER", "WSHP"};
    // int weightmethod=2;//0=exp,1=power,2=wshp
    // double binres=2;

    cout << "\n"
         << "processing:" << fileName << " Histogram: " << histName << "\n";
    TFile *pionfile = new TFile(fileName, "READ");
    if (!pionfile || pionfile->IsZombie())
    {
        std::cerr << "Error opening file." << std::endl;
        // Handle error or return
    }
    // TH2F* hist2D = dynamic_cast<TH2F*>(pionfile->Get(histName));
    TH2F *hist2D = (TH2F *)pionfile->Get(histName);
    if (!hist2D)
    {
        std::cerr << "Histogram not found." << std::endl;
        // Handle error or return
    }

    TString canvasname = Form("%s_Esmear_%d_Thousandths", histName, filenameobj.sqrtEsmearing[0]); // remvod %s for Time
    const char *pdfname = canvasname;
    TCanvas *c1 = new TCanvas(canvasname, canvasname, 3000, 1200);

    // c1->SetFillColor(42);
    // c1->Divide(2, 1); // nx, ny
    /////////////////////////////////////////////weighted

    c1->Divide(2, 3);
    c1->cd(1);
    // c1->cd(1);
    gPad->SetTopMargin(0.12);
    gPad->SetFillColor(33);
    gPad->SetLogz();
    hist2D->Draw("colz");
    // printf("error test code\n");
    hist2D->GetXaxis()->SetLabelSize(0.06);
    hist2D->GetYaxis()->SetLabelSize(0.06);
    hist2D->GetXaxis()->SetTitle(Form("%s Pt [GeV/c]", filenameobj.particletype.c_str()));
    hist2D->GetYaxis()->SetTitle("Invariant Mass [GeV/c^2]");
    hist2D->SetMarkerColor(kYellow);
    // Fit slices projected along Y fron bins in X [1,64] with more than 2 bins in Y filled
    TObjArray aSlices; // owner of slices?
    hist2D->GetXaxis()->SetRangeUser(0, 30);
    hist2D->FitSlicesY(0, 0, -1, 0, "QNR", &aSlices); //, "EMW" "QNR",L-log likelihood

    aSlices.SetOwner(kTRUE);
    //

    // Show fitted "mean" for each slice
    c1->cd(2);
    gPad->SetFillColor(33);

    // Create a 1D histogram for the integrals
    int nBinsX = hist2D->GetXaxis()->GetNbins();
    TH1D *hIntegrals = new TH1D("hIntegrals", "Integrals over Y;X-axis;Integral", nBinsX, hist2D->GetXaxis()->GetXmin(), hist2D->GetXaxis()->GetXmax());

    // note this gives the weighted entries!!!!!
    // Loop over each bin in X
    for (int i = 1; i <= nBinsX; ++i)
    {
        // Integrate over all bins in Y for the current bin in X
        double integral = hist2D->Integral(i, i, 1, hist2D->GetYaxis()->GetNbins());
        hIntegrals->SetBinContent(i, integral);
    }
    hIntegrals->GetXaxis()->SetTitle(Form("%s Pt [GeV/c]", filenameobj.particletype.c_str()));
    gPad->SetLogy();
    hIntegrals->Draw();

    c1->cd(3);
    gPad->SetTopMargin(0.12);
    gPad->SetLeftMargin(0.15);
    gPad->SetFillColor(33);
    TH2F *hist2D_1 = (TH2F *)aSlices[1];
    // TH2F *hist2D_1 = (TH2F *)pionfile->Get(Form("%s_1",histName));
    hist2D_1->GetYaxis()->SetTitle("Mean");
    hist2D_1->GetXaxis()->SetTitle(Form("%s Pt [GeV/c]", filenameobj.particletype.c_str()));
    // hist2D_1->SetAxisRange(0.1, 0.16,"Y");
    hist2D_1->Draw();

    // Show fitted "sigma" for each slice
    // c1->cd(2);
    c1->cd(4);
    gPad->SetTopMargin(0.12);
    gPad->SetLeftMargin(0.15);
    gPad->SetFillColor(33);
    TH2F *hist2D_2 = (TH2F *)aSlices[2];
    // TH2F *hist2D_2 = (TH2F *)pionfile->Get(Form("%s_2",histName));
    hist2D_2->SetMinimum(0.8);
    hist2D_2->GetYaxis()->SetTitle("Sigma");
    hist2D_2->GetXaxis()->SetTitle(Form("%s Pt [GeV/c]", filenameobj.particletype.c_str()));
    hist2D_2->Draw();

    // Show chi^2/ndf
    c1->cd(5);
    gPad->SetTopMargin(0.12);
    gPad->SetLeftMargin(0.15);
    gPad->SetFillColor(33);
    // TH2F *hist2D_0 = (TH2F *)pionfile->Get(Form("%s_0",histName));
    TH2F *hist2D_chi2 = (TH2F *)pionfile->Get(Form("%s_chi2", histName));
    // TH2F *hist2D_chi2 = (TH2F *)aSlices[3];
    if (!hist2D_chi2)
    {
        std::cerr << "Histogram " << Form("%s_chi2", histName) << " not found. Exiting SliceAndFit" << std::endl;
        return;
        // Handle error or return
    }
    hist2D_chi2->GetXaxis()->SetTitle(Form("%s Pt [GeV/c]", filenameobj.particletype.c_str()));
    // gPad->SetLogy();
    hist2D_chi2->Draw();
    // hist2D_chi2->SetMinimum(0.8);
    hist2D_chi2->SetTitle(Form("chi^2/NDF;%s Pt [GeV/c];chi^2/NDF", filenameobj.particletype.c_str()));

    // Show fitted sigma/mean for each slice
    c1->cd(6);
    TH2F *hist2D_4 = (TH2F *)hist2D_2->Clone(Form("%s_4", histName)); // clone sigma
    gPad->SetTopMargin(0.12);
    gPad->SetLeftMargin(0.15);
    gPad->SetFillColor(33);
    hist2D_4->Divide(hist2D_1); // divide by mean
    hist2D_4->SetTitle(Form("Value of par[2]/par[1]=Sigma/Mean;%s Pt [GeV/c];Sigma/Mean", filenameobj.particletype.c_str()));
    // hist2D_4->SetMinimum(0.8);
    hist2D_4->Draw();

    c1->SaveAs(Form("pioncode/canvas_pdf/%s.pdf", canvasname.Data()));

    delete c1;
    delete hIntegrals;
    aSlices.Delete();
    pionfile->Close();
    delete pionfile;
}

// misc operations
void transferHistogram(const char *sourceFileName, const char *histogramName, const char *targetFileName, const char *NewhistogramName)
{
    // Open the source file
    TFile sourceFile(sourceFileName, "READ");
    if (!sourceFile.IsOpen())
    {
        std::cerr << "Error: Unable to open source file!" << std::endl;
        return;
    }

    // Retrieve the histogram
    TObject *obj = sourceFile.Get(histogramName);
    if (!obj)
    {
        std::cerr << "Error: Histogram not found in source file!" << std::endl;
        sourceFile.Close();
        return;
    }

    TH1 *histogram = dynamic_cast<TH1 *>(obj);
    if (!histogram)
    {
        std::cerr << "Error: Object is not a histogram!" << std::endl;
        sourceFile.Close();
        return;
    }

    // Detach the histogram from the source file
    histogram->SetDirectory(0);

    // Close the source file
    sourceFile.Close();

    // Change the name of the histogram
    histogram->SetName(NewhistogramName);

    // Open the target file
    TFile targetFile(targetFileName, "UPDATE");
    if (!targetFile.IsOpen())
    {
        std::cerr << "Error: Unable to open target file!" << std::endl;
        return;
    }

    // Write the histogram to the target file
    histogram->Write();

    // Close the target file
    targetFile.Close();
}
