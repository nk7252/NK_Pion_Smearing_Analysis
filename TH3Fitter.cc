#include <TH1.h>
#include <TH1F.h>
#include <TH3.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <Math/MinimizerOptions.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

struct FitConfig {
    int xBinStart, xBinEnd, yBinStart, yBinEnd;
    double minInvMass, maxInvMass;

    // Constructor
    FitConfig(int xs, int xe, int ys, int ye, double minM, double maxM)
        : xBinStart(xs), xBinEnd(xe), yBinStart(ys), yBinEnd(ye), minInvMass(minM), maxInvMass(maxM) {}

    // The serialize and deserialize methods remain unchanged
    std::string serialize() const {
        std::stringstream ss;
        ss << xBinStart << "," << xBinEnd << ","
           << yBinStart << "," << yBinEnd << ","
           << minInvMass << "," << maxInvMass;
        return ss.str();
    }

    static FitConfig deserialize(const std::string& line) {
        std::stringstream ss(line);
        std::vector<std::string> result;
        while (ss.good()) {
            std::string substr;
            getline(ss, substr, ',');
            result.push_back(substr);
        }
        return {std::stoi(result[0]), std::stoi(result[1]),
                std::stoi(result[2]), std::stoi(result[3]),
                std::stod(result[4]), std::stod(result[5])};
    }
};

class LeftRightPolynomial {
public:
    double minMass;
    double maxMass;

    LeftRightPolynomial(double minM, double maxM) : minMass(minM), maxMass(maxM) {}

    double operator()(double *x, double *par) {
        if (x[0] >= minMass && x[0] <= maxMass) {
            TF1::RejectPoint();
            return 0;
        }
        // Polynomial (4th degree) calculation
        return par[0] + par[1]*x[0] + par[2]*x[0]*x[0] + par[3]*x[0]*x[0]*x[0] + par[4]*x[0]*x[0]*x[0]*x[0];
    }
};

// Global pointer to an instance of LeftRightPolynomial
// allows for each call(sequential) to pass the corresponding min and max mass
static LeftRightPolynomial* gLeftRightPoly = nullptr;

//declarations
std::vector<FitConfig> LoadOrCreateConfig(const std::string& configFile);
void FitAndSaveProjection(TH3* h3, const FitConfig& config, const std::string& outputPDF, std::ofstream& csvFile);
void AnalyzeAndFit(const std::string& rootFileName, const std::string& histName);
// automated fitting
std::vector<double> FitAndGetParams(TH1D* hProjZ, double minMass, double maxMass);
std::vector<TCanvas*> OptimizeFitRange(TH3* h3, int xBinStart, int xBinEnd, int yBinStart, int yBinEnd);
void OptimizeHistogramFit(const std::string& rootFileName, const std::string& histogramName);
//for fitting gaus+poly4
Double_t combinedFunction(Double_t *x, Double_t *par);
double LeftRightPolynomialBridge(double *x, double *par);
void appendtextfile(TF1* fitFunc, const std::string& fitName);
std::vector<TCanvas*> DrawBestHistogram(TH1D* hProjZ, double minMass, double maxMass);



void TH3Fitter(){
std::string rootFileName="pioncode/rootfiles/data/pt_nclus_differential_data/pt05pt05.root";
std::string histName= "h_pipT_Nclus_mass";

//AnalyzeAndFit(rootFileName, histName);
OptimizeHistogramFit(rootFileName, histName);

}


std::vector<FitConfig> LoadOrCreateConfig(const std::string& configFile) {
    std::vector<FitConfig> configs;
    std::ifstream file(configFile);
    std::string line;

    if (file) {
        // Skip the first line which contains the header
        std::getline(file, line);

        // Load configurations
        while (std::getline(file, line)) {
            configs.push_back(FitConfig::deserialize(line));
        }
    } else {
        // Create default configuration with header
        std::ofstream outFile(configFile);
        // Write the header line first
        outFile << "pTBinStart,pTBinEnd,nclusBinStart,nclusBinEnd,minInvMass,maxInvMass\n";

        // Example: Assuming 10 bins in X and Y, and invariant mass range [0.1, 1.0]
        for (int xBinStart = 1; xBinStart <= 10; xBinStart += 2) {
            for (int yBinStart = 1; yBinStart <= 10; yBinStart += 2) {
                FitConfig config(xBinStart, xBinStart + 1, yBinStart, yBinStart + 1, 0.1, 0.3);
                outFile << config.serialize() << std::endl;
            }
        }
    }
    return configs;
}

// Function to fit a projection and save results
void FitAndSaveProjection(TH3* h3, const FitConfig& config, const std::string& outputPDF, std::ofstream& csvFile) {
    // Assume h3 is your TH3 histogram
    TH1D* hProjZ = h3->ProjectionZ("projZ", config.xBinStart, config.xBinEnd, config.yBinStart, config.yBinEnd);
    TF1* fitFunc = new TF1("fitFunc", "gaus", config.minInvMass, config.maxInvMass);
    hProjZ->Fit(fitFunc, "RQ"); // "RQ" option for Range and Quiet


    hProjZ->SetTitle("Pion Inv. Mass for: nclus ,pT ;Pion Inv. Mass GeV;Counts");
    // Extract and print fit parameters
    double mean = fitFunc->GetParameter(1);
    double sigma = fitFunc->GetParameter(2);
    double chi2 = fitFunc->GetChisquare();
    double ndf = fitFunc->GetNDF();
    double resolution = 100 * sigma / mean;

    // Append fit results to CSV
    csvFile << config.xBinStart << "," << config.xBinEnd << "," << config.yBinStart << "," << config.yBinEnd << "," << mean << "," << sigma << "," << resolution << "," << chi2/ndf << std::endl;

    // Draw and save to PDF
    TCanvas c;
    hProjZ->Draw();
    c.Print(outputPDF.c_str(), "pdf");
}

// Main analysis function
void AnalyzeAndFit(const std::string& rootFileName, const std::string& histName) {
    // Prepare file paths
    std::string configFilePath = "pioncode/csvfiles/fitConfig.csv";
    std::string outputPDFPath = "pioncode/canvas_pdf/fitResults.pdf";
    std::string outputCSVPath = "pioncode/csvfiles/fitResults.csv";

    // Open the ROOT file and retrieve the histogram
    TFile* file = TFile::Open(rootFileName.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << rootFileName << std::endl;
        return;
    }
    TH3* h3 = dynamic_cast<TH3*>(file->Get(histName.c_str()));
    if (!h3) {
        std::cerr << "Histogram " << histName << " not found in file." << std::endl;
        file->Close();
        return;
    }
    // Load or create the configuration
    auto configs = LoadOrCreateConfig(configFilePath);

    // Setup output files
    std::ofstream csvFile(outputCSVPath);
    csvFile << "pTBinStart,pTBinEnd,NclusBinStart,NclusBinEnd,Mean,Sigma,Resolution,Chi2/NDF\n";
    TCanvas* c3 = new TCanvas("c3", "Fits", 800, 600);
    std::string pdfName = outputPDFPath;
    c3->Print((pdfName + "[").c_str()); // Open the PDF document

    // Perform fits and save results
    for (const auto& config : configs) {
          std::cout << "Fitting with config: "
              << "xBinStart=" << config.xBinStart << ", "
              << "xBinEnd=" << config.xBinEnd << ", "
              << "yBinStart=" << config.yBinStart << ", "
              << "yBinEnd=" << config.yBinEnd << ", "
              << "minInvMass=" << config.minInvMass << ", "
              << "maxInvMass=" << config.maxInvMass << std::endl;
        FitAndSaveProjection(h3, config, pdfName, csvFile);
    }

    // Finalize PDF and clean up
    c3->Print((pdfName + "]").c_str()); // Close the PDF document
    delete c3; // Clean up the canvas
    file->Close(); // Close the ROOT file
}

// Function to perform the Gaussian fit and return chi^2/NDF
std::vector<double> FitAndGetParams(TH1D* hProjZ, double minMass, double maxMass) {
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);


    // Rebin the histogram to have 'numBins' bins
    // First, calculate the rebin factor assuming the histogram's range is 0 to maxXRange
    /*
    int numBins= 120;
    int currentNumBins = hist->GetNbinsX();
    double currentXMax = hist->GetXaxis()->GetXmax();
    int rebinFactor = currentNumBins / numBins;
    if (rebinFactor > 1) { // Only rebin if the factor is greater than 1
        std::cout << "current nbins: " << currentNumBins <<" requested nbins: " << numBins << " rebin by: " << rebinFactor << std::endl;
        hist->Rebin(rebinFactor);
        std::cout << "new nbin check: " << hist->GetNbinsX() << std::endl;
    }
    */

    float leftlimit =0.0;
    float rightlimit =0.5;

    //hProjZ->GetXaxis()->SetRangeUser(0, 0.4);

    // Fit left and right regions with a polynomial, excluding Gaussian region

    LeftRightPolynomial polyFunc(minMass, maxMass); //exclusion range
    gLeftRightPoly = &polyFunc; // Point the global pointer to your instance
    TF1 *leftRightFit = new TF1("leftRightFit", LeftRightPolynomialBridge, leftlimit, rightlimit, 5);

    hProjZ->Fit(leftRightFit, "RQ0");// "RQ" option for Range and Quiet, 0 for do not display fit on canvas.
    // Fit Gaussian in the specified range
    TF1 *gausFit = new TF1("gausFit", "gaus", minMass, maxMass);//leftpolylim, rightpolylim
    hProjZ->Fit(gausFit, "RQ0");
    // Combined Gaussian + Polynomial fit
    TF1 *combinedFit = new TF1("combinedFit", combinedFunction, leftlimit, rightlimit, 8);
    // Set initial parameters from previous fits
    for (int i = 0; i < 3; ++i) combinedFit->SetParameter(i, gausFit->GetParameter(i));
    for (int i = 3; i < 8; ++i) combinedFit->SetParameter(i, leftRightFit->GetParameter(i-3));
    //try to improve the fit.
    hProjZ->Fit(combinedFit, "RQ0");//L//M
    //-------------------------------------------show the poly4 part seperately
    // Create a new function for just the polynomial part
    TF1 *polyPart = new TF1("polyPart", "pol4", leftlimit, rightlimit);
    // Set the parameters of polyPart to those from the combined fit
    // Assuming the first 5 parameters of combinedFit are for the polynomial
    for (int i = 0; i < 5; ++i) polyPart->SetParameter(i, combinedFit->GetParameter(i+3));
    // Create a new histogram to store the subtracted data
    TH1F *histSubtracted = (TH1F*)hProjZ->Clone("histSubtracted");
    // Subtract the polynomial part
    for (int i = 1; i <= hProjZ->GetNbinsX(); ++i) {
        double x = hProjZ->GetBinCenter(i);
        double y = hProjZ->GetBinContent(i) - polyPart->Eval(x);
        histSubtracted->SetBinContent(i, y);
    }
    TF1 *gausFit2 = new TF1("gausFit2", "gaus", minMass, maxMass);//leftmost_limit, 0.25
    for (int i = 0; i < 3; ++i) gausFit2->SetParameter(i, combinedFit->GetParameter(i));
    histSubtracted->Fit(gausFit2, "RQ0");//L
    //double chi2_s = gausFit2->GetChisquare();
    //double ndf_s = gausFit2->GetNDF();
    //double chi2ndf_s = chi2_s / ndf_s;


    //append // Append fit results to text file
    appendtextfile(gausFit, "gausFit");
    appendtextfile(combinedFit, "Combined Fit");
    appendtextfile(gausFit2, "subpgaus fit");




    std::vector<double> fitparams;
    // add fit results to vector output
    fitparams.push_back(gausFit2->GetChisquare());
    fitparams.push_back(gausFit2->GetNDF());
    fitparams.push_back(gausFit2->GetChisquare()/gausFit2->GetNDF());
    fitparams.push_back(gausFit2->GetParameter(2));
    fitparams.push_back(gausFit2->GetParameter(1));
    fitparams.push_back(gausFit2->GetParameter(2)* 100.0f / gausFit2->GetParameter(1));


    // Clean up
    delete leftRightFit;
    delete gausFit;
    delete combinedFit;
    delete polyPart;
    delete gausFit2;
    
    return fitparams;
}

// Main function to explore fit ranges and find optimal parameters
std::vector<TCanvas*> OptimizeFitRange(TH3* h3, int xBinStart, int xBinEnd, int yBinStart, int yBinEnd) {
    TH1D* hProjZ = h3->ProjectionZ("projZ", xBinStart, xBinEnd, yBinStart, yBinEnd);

    double bestChi2NDF = TMath::Infinity();
    double bestMinMass = 0, bestMaxMass = 0;

    double startMass = 0.1; // start of mass range
    double endMass = 0.25; // end of mass range
    double step = 0.005; // step size

    for (double minMass = startMass; minMass < endMass; minMass += step) {
        for (double maxMass = minMass + step; maxMass <= endMass; maxMass += step) {
            std::vector<double> fitresults;
            fitresults= FitAndGetParams(hProjZ, minMass, maxMass);
            double chi2NDF = fitresults[2];// chi2/ndf
            

            if (chi2NDF >= 0.9 && chi2NDF < bestChi2NDF) {
                bestChi2NDF = chi2NDF;
                bestMinMass = minMass;
                bestMaxMass = maxMass;
                //if (chi2NDF >= 0 && chi2NDF <= 10){

                //}
            }
        }
    }


    // Output the best fit range and chi^2/NDF
    std::cout << "Best Fit Range: [" << bestMinMass << ", " << bestMaxMass << "] with chi^2/NDF = " << bestChi2NDF << std::endl;

    // Optionally, perform and visualize the final fit with the best parameters
    std::vector<TCanvas*> canvases2;
    canvases2=DrawBestHistogram(hProjZ, bestMinMass, bestMaxMass);
    return canvases2;
}

void OptimizeHistogramFit(const std::string& rootFileName, const std::string& histogramName) {
    // Open the ROOT file
    TFile* file = TFile::Open(rootFileName.c_str(), "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file: " << rootFileName << std::endl;
        return;
    }

    // Load the histogram
    TH3* h3 = dynamic_cast<TH3*>(file->Get(histogramName.c_str()));
    if (!h3) {
        std::cerr << "Histogram " << histogramName << " not found in file " << rootFileName << std::endl;
        file->Close();
        delete file;
        return;
    }

    // Call the function to optimize the fit range 
    float xBinStart = 1, xBinEnd = h3->GetXaxis()->GetNbins();
    //int yBinStart = 1, yBinEnd =1;
    int nclusbinwidth= h3->GetYaxis()->GetBinWidth(1);///h3->GetYaxis()->GetNbins();
    float pTbinwidth= h3->GetXaxis()->GetBinWidth(1);

    TCanvas *canvas = new TCanvas("canvas", "Canvas", 800, 600);

    //open pdf
    canvas->Print(Form("pioncode/canvas_pdf/%s_Fit.pdf[",histogramName.c_str()));

    // look at ncluster ranges(bins of width 20 nclus)
    for(int yBinEnd=1; yBinEnd<=h3->GetYaxis()->GetNbins(); yBinEnd++){
        for(int yBinStart=1; yBinStart<=yBinEnd; yBinStart++){
            TCanvas *textCanvas = new TCanvas("textCanvas", "Canvas Info", 800, 600);
            textCanvas->cd();
            TPaveText* pt0 = new TPaveText(0.1, 0.1, 0.9, 0.9, "blNDC"); // blNDC: borderless, normalized coordinates
            pt0->SetTextAlign(12); // Align text to the left
            pt0->SetFillColor(0); // Transparent background

            // Adding custom text entries
            //pt0->AddText(Form("Fits to File = %s", histogramName.c_str()));
            pt0->AddText(Form("Fits to Histogram = %s", histogramName.c_str()));
            pt0->AddText(Form("For pT (GeV) = %.1f - %.1f ",(xBinStart-1)*pTbinwidth ,xBinEnd*pTbinwidth));
            pt0->AddText(Form("For nClusters = %.1d - %.1d",(yBinStart-1)*nclusbinwidth ,yBinEnd*nclusbinwidth));
            pt0->Draw();
            textCanvas->Print(Form("pioncode/canvas_pdf/%s_Fit.pdf",histogramName.c_str()));
            delete textCanvas;
            delete pt0;

            std::vector<TCanvas*> canvases3;
            canvases3=OptimizeFitRange(h3, xBinStart, xBinEnd, yBinStart, yBinEnd);

            for (auto* canvascol : canvases3) {//(int i=0;i<canvases3.size();i++){
                canvascol->Print(Form("pioncode/canvas_pdf/%s_Fit.pdf",histogramName.c_str()));
                delete canvascol;
            }
        }
    }
    


    //close pdf
    canvas->Print(Form("pioncode/canvas_pdf/%s_Fit.pdf]",histogramName.c_str()));

    delete canvas;
    canvases3.clear();
    // Close the ROOT file
    file->Close();
    delete file;
}


//fitting functions 
// Combined function for Gaussian + Polynomial
Double_t combinedFunction(Double_t *x, Double_t *par) {
    // Gaussian part
    Double_t gauss = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));
    
    // Polynomial part
    Double_t poly = par[3] + par[4]*x[0] + par[5]*x[0]*x[0] + par[6]*x[0]*x[0]*x[0] + par[7]*x[0]*x[0]*x[0]*x[0];
    
    return gauss + poly;
}

// Bridge function for simultaneous left and right polynomial fit, excluding the Gaussian region
double LeftRightPolynomialBridge(double *x, double *par) {
    if (gLeftRightPoly) return (*gLeftRightPoly)(x, par);
    return 0; // Just in case gLeftRightPoly is not set
}

void appendtextfile(TF1* fitFunc, const std::string& fitName){//, Double_t scale_factor
    // Open a text file in append mode
    std::ofstream outfile;
    outfile.open("fit_parameters.txt", std::ios_base::app);

    // Check if the file is open (and thus valid)
    if (outfile.is_open()) {
        // Check if the file is empty
        outfile.seekp(0, std::ios::end);
        size_t size = outfile.tellp();

        if (size == 0) {
            // File is new, write header
            outfile << "Fit Name          Mean          Sigma          Sigma/Mean          chi2          NDF           chi2/NDF" << std::endl;
        }

        // Write the parameters to the file
        outfile << fitName << "     "  << fitFunc->GetParameter(1) << "     " << fitFunc->GetParameter(2) << "     " << fitFunc->GetParameter(2)/fitFunc->GetParameter(1) << "     " << fitFunc->GetChisquare() << "     " << fitFunc->GetNDF() <<"     " << fitFunc->GetChisquare()/fitFunc->GetNDF() << std::endl;//<< scale_factor << "     "

        outfile.close();
    } else {
        std::cerr << "Error: Unable to open file." << std::endl;
    }
}


std::vector<TCanvas*> DrawBestHistogram(TH1D* hProjZ, double minMass, double maxMass) {
    // more thorough minimizer for fit
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    // Set the global fit strategy
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);


    float leftlimit =0.0;
    float rightlimit =0.5;
    //hist->GetXaxis()->SetRangeUser(0, 0.4);

    LeftRightPolynomial polyFunc(minMass, maxMass); //exclusion range
    gLeftRightPoly = &polyFunc; // Point the global pointer to your instance
    TF1 *leftRightFit = new TF1("leftRightFit", LeftRightPolynomialBridge, leftlimit, rightlimit, 5);

    hProjZ->Fit(leftRightFit, "R");

    // Fit Gaussian in the specified range
    TF1 *gausFit = new TF1("gausFit", "gaus", minMass, maxMass);//leftpolylim, rightpolylim
    hProjZ->Fit(gausFit, "R");


    // Combined Gaussian + Polynomial fit
    TF1 *combinedFit = new TF1("combinedFit", combinedFunction, leftlimit, rightlimit, 8);
    // Set initial parameters from previous fits
    for (int i = 0; i < 3; ++i) combinedFit->SetParameter(i, gausFit->GetParameter(i));
    for (int i = 3; i < 8; ++i) combinedFit->SetParameter(i, leftRightFit->GetParameter(i-3));
    //try to improve the fit.
    hProjZ->Fit(combinedFit, "R");//M
    //double chi2 = combinedFit->GetChisquare();
    //double ndf = combinedFit->GetNDF();
    //double chi2ndf = chi2 / ndf;

    //std::cout << "Chi-squared: " << chi2 << std::endl;
    //std::cout << "Number of Degrees of Freedom: " << ndf << std::endl;
    //std::cout << "Chi-squared/NDF: " << chi2ndf << std::endl;

    //-------------------------------------------show the poly4 part seperately


    // Create a new function for just the polynomial part
    TF1 *polyPart = new TF1("polyPart", "pol4", 0, 0.5);
    //same for gaussian part
    TF1 *GaussPart = new TF1("GaussPart", "gaus", minMass, maxMass);

    // Set the parameters of polyPart to those from the combined fit
    // Assuming the last 5 parameters of combinedFit are for the polynomial
    for (int i = 0; i < 5; ++i) polyPart->SetParameter(i, combinedFit->GetParameter(i+3));

    for (int i = 0; i < 3; ++i) GaussPart->SetParameter(i, combinedFit->GetParameter(i));
    
    // Create a new histogram to store the subtracted data
    TH1F *histSubtracted = (TH1F*)hProjZ->Clone("histSubtracted");

    // Subtract the polynomial part
    for (int i = 1; i <= hProjZ->GetNbinsX(); ++i) {
        double x = hProjZ->GetBinCenter(i);
        double y = hProjZ->GetBinContent(i) - polyPart->Eval(x);
        histSubtracted->SetBinContent(i, y);
    }
    TF1 *gausFit2 = new TF1("gausFit2", "gaus", minMass, maxMass);//leftmost_limit, 0.25
    for (int i = 0; i < 3; ++i) gausFit2->SetParameter(i, combinedFit->GetParameter(i));
    histSubtracted->Fit(gausFit2, "R");
    //double chi2_s = gausFit2->GetChisquare();
    //double ndf_s = gausFit2->GetNDF();
    //double chi2ndf_s = chi2_s / ndf_s;

    //std::cout << "Chi-squared: " << chi2_s << std::endl;
    //std::cout << "Number of Degrees of Freedom: " << ndf_s  << std::endl;
    //std::cout << "Chi-squared/NDF: " << chi2ndf_s  << std::endl;
    //std::cout << "Relative Width: " << 100*gausFit2->GetParameter(2)/gausFit2->GetParameter(1)  << " %" << std::endl;

    //store 2 separate functions for visualization
    TF1 *fleft = new TF1("fleft",LeftRightPolynomialBridge, leftlimit, minMass, 5);
    fleft->SetParameters(leftRightFit->GetParameters());
    //hist->GetListOfFunctions()->Add(fleft);
    //gROOT->GetListOfFunctions()->Remove(fleft);
    TF1 *fright = new TF1("fright",LeftRightPolynomialBridge, maxMass, rightlimit, 5);
    fright->SetParameters(leftRightFit->GetParameters());
    //hist->GetListOfFunctions()->Add(fright);
    //gROOT->GetListOfFunctions()->Remove(fright);

    
    // Create graphs for residuals
    TGraph *residuals1 = new TGraph(hProjZ->GetNbinsX());
    TGraph *residuals2 = new TGraph(hProjZ->GetNbinsX());

    for (int i = 1; i < hProjZ->GetNbinsX(); ++i) {
        // Get the bin center
        double bincenter1 = hProjZ->GetBinCenter(i);
        // Calculate the residual (Data - Fit)
        double residual1 = hProjZ->GetBinContent(i) - combinedFit->Eval(bincenter1);
        residuals1->SetPoint(i, bincenter1, residual1); // Set the point in the residuals graph
        //
        double bincenter2 = histSubtracted->GetBinCenter(i);
        double residual2 = histSubtracted->GetBinContent(i) - gausFit2->Eval(bincenter2);
        residuals2->SetPoint(i, bincenter2, residual2); // Set the point in the residuals graph

    }


    // Draw everything and add canvases to vector of canvases
    std::vector<TCanvas*> canvases;

    //-------------------------------------------------------------------------------------------canvas 1
    TCanvas *c1 = new TCanvas("c1", "Fits", 800, 600);
    hProjZ->SetTitle("Combined Fit (Gaus+Poly4); Inv. Mass (GeV); Counts");
    //hProjZ->SetStats(0); // Turn off stat box
    hProjZ->Draw("E");

    gausFit->SetLineColor(kGreen);
    gausFit->Draw("SAME");// draw the gaussian fit
    polyPart->SetLineColor(kRed);
    polyPart->Draw("SAME");
    //GaussPart->SetLineColor(kGreen);
    //GaussPart->Draw("SAME");

    //leftRightFit->SetLineColor(kBlue);
    fleft->SetLineColor(kBlue);
    fright->SetLineColor(kBlue);

    fleft->Draw("SAME");
    fright->Draw("SAME");// turn off to see the inflection better.
    //leftRightFit->Draw("SAME"); // Draw the left and right polynomial fits

    combinedFit->SetLineColor(kBlack);
    combinedFit->Draw("SAME"); // draw the combined fit

    // Add a legend
    TLegend *leg = new TLegend(0.1, 0.7, 0.3, 0.9);//bot left x, bot left y, top right x, top right y 
    leg->AddEntry(gausFit, "Gaussian Fit");
    leg->AddEntry(polyPart, "Polynomial Fit");
    leg->AddEntry(fleft, "Left & Right Polynomial Fit");
    leg->AddEntry(combinedFit, "Combined Fit");
    leg->Draw();

    canvases.push_back(c1);
    c1->SaveAs("combined_fits.pdf");

    //-------------------------------------------------------------------------------------------canvas 2
    TCanvas *c2 = new TCanvas("c2", "Subtracted Peak", 800, 600);
    histSubtracted->SetTitle("Peak after Background Subtraction; Inv. Mass (GeV); Counts");
    histSubtracted->SetStats(0); // Turn off stat box
    histSubtracted->Draw("E");
    histSubtracted->SetMinimum(0.0);


    canvases.push_back(c2);
    //c2->SaveAs("Subtracted_Peak.pdf");

    // Second canvas: Custom list of fit results
    TCanvas* c3 = new TCanvas("c3", "Fit Parameters", 800, 600);
    c3->cd();

    TPaveText* pt = new TPaveText(0.1, 0.1, 0.9, 0.9, "blNDC"); // blNDC: borderless, normalized coordinates
    pt->SetTextAlign(12); // Align text to the left
    pt->SetFillColor(0); // Transparent background

    // Adding custom text entries
    pt->AddText("Data Fit");
    pt->AddText("Fit Parameters:");
    pt->AddText(Form("Combined Fit Range = %f to %f", leftlimit, rightlimit));
    pt->AddText(Form("Peak Mean = %f +/- %f", combinedFit->GetParameter(1), combinedFit->GetParError(1)));
    pt->AddText(Form("Peak Sigma = %f +/- %f", combinedFit->GetParameter(2), combinedFit->GetParError(2)));
    pt->AddText(Form("Peak Relative Width: %f",combinedFit->GetParameter(2)* 100.0f / combinedFit->GetParameter(1)));   
    pt->AddText(Form("Peak Chi2/NDF = %f / %d= %f", combinedFit->GetChisquare(), combinedFit->GetNDF(),combinedFit->GetChisquare()/combinedFit->GetNDF()));
    pt->AddText(Form("Background Subtracted Peak Fit = %f to %f", minMass, maxMass));
    pt->AddText(Form("Mean = %f +/- %f", gausFit2->GetParameter(1), gausFit2->GetParError(1)));
    pt->AddText(Form("Sigma = %f +/- %f", gausFit2->GetParameter(2), gausFit2->GetParError(2)));
    pt->AddText(Form("Relative Width: %f",gausFit2->GetParameter(2)* 100.0f / gausFit2->GetParameter(1)));   
    pt->AddText(Form("Chi2/NDF = %f / %d= %f", gausFit2->GetChisquare(), gausFit2->GetNDF(),gausFit2->GetChisquare()/gausFit2->GetNDF()));
    pt->Draw();

    canvases.push_back(c3);
    //c3->SaveAs("FitInfo.pdf");
    TCanvas* c4 = new TCanvas("c4", "CombinedFit Residuals", 800, 600);
    c4->cd();
    residuals1->SetTitle("Residuals for Combined Fit (Gaus+Poly4); Inv. Mass (GeV); Residual Counts");
    residuals1->Draw("AP");
    residuals1->SetMarkerSize(1); // Increase the marker size
    residuals1->SetMarkerStyle(27); // Increase the marker size
    TLine *line = new TLine(leftlimit, 0, rightlimit, 0);
    line->SetLineColor(kRed);
    line->Draw("same");
    c4->Update();


    TCanvas* c5 = new TCanvas("c5", "Background Subtracted Peak Residuals", 800, 600);
    c5->cd();
    residuals2->SetTitle("Residuals for Peak after Background Subtraction; Inv. Mass (GeV); Residual Counts");
    residuals2->Draw("AP");
    residuals2->SetMarkerSize(1); // Increase the marker size
    residuals2->SetMarkerStyle(28); // Increase the marker size
    line->Draw("same");
    c5->Update();

    canvases.push_back(c4);
    canvases.push_back(c5);




    //clean up
    delete gausFit;
    delete gausFit2;
    delete polyPart;
    delete leftRightFit;
    delete fleft;
    delete fright;
    delete combinedFit;
    delete leg;

    return canvases;
}


