#include <TH1.h>
#include <TH3.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
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
std::vector<double> FitAndGetParams(TH1D* hProjZ, double minMass, double maxMass);
void OptimizeFitRange(TH3* h3, int xBinStart, int xBinEnd, int yBinStart, int yBinEnd);
void OptimizeHistogramFit(const std::string& rootFileName, const std::string& histogramName);
//for fitting gaus+poly4
Double_t combinedFunction(Double_t *x, Double_t *par);
double LeftRightPolynomialBridge(double *x, double *par);
void appendtextfile(TF1* fitFunc, const std::string& fitName);
//void DrawBestHistogram( float leftmost_gauslimit, float rightmost_gauslimit);



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

    // limits on gauss and poly
    //float leftpolylim = 0.11;
    //float rightpolylim = 0.19;
    //hProjZ->GetXaxis()->SetRangeUser(0, 0.4);

    // Fit left and right regions with a polynomial, excluding Gaussian region

    LeftRightPolynomial polyFunc(minMass, maxMass); //exclusion range
    gLeftRightPoly = &polyFunc; // Point the global pointer to your instance
    TF1 *leftRightFit = new TF1("leftRightFit", LeftRightPolynomialBridge, 0, 0.5, 5);

    hProjZ->Fit(leftRightFit, "RQ0");// "RQ" option for Range and Quiet, 0 for do not display fit on canvas.
    // Fit Gaussian in the specified range
    TF1 *gausFit = new TF1("gausFit", "gaus", minMass, maxMass);//leftpolylim, rightpolylim
    hProjZ->Fit(gausFit, "RQ0");
    // Combined Gaussian + Polynomial fit
    TF1 *combinedFit = new TF1("combinedFit", combinedFunction, minMass, maxMass, 8);
    // Set initial parameters from previous fits
    for (int i = 0; i < 3; ++i) combinedFit->SetParameter(i, gausFit->GetParameter(i));
    for (int i = 3; i < 8; ++i) combinedFit->SetParameter(i, leftRightFit->GetParameter(i-3));
    //try to improve the fit.
    hProjZ->Fit(combinedFit, "RQL0");//L//M
    //-------------------------------------------show the poly4 part seperately
    // Create a new function for just the polynomial part
    TF1 *polyPart = new TF1("polyPart", "pol4", minMass, maxMass);
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
    histSubtracted->Fit(gausFit2, "RQ0");
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
void OptimizeFitRange(TH3* h3, int xBinStart, int xBinEnd, int yBinStart, int yBinEnd) {
    TH1D* hProjZ = h3->ProjectionZ("projZ", xBinStart, xBinEnd, yBinStart, yBinEnd);

    double bestChi2NDF = TMath::Infinity();
    double bestMinMass = 0, bestMaxMass = 0;

    double startMass = 0.05; // start of mass range
    double endMass = 0.4; // end of mass range
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
    //std::cout << "Chi-squared: " << gaussFit->GetChisquare() << std::endl;
    //std::cout << "Number of Degrees of Freedom: " << gaussFit->GetNDF() << std::endl;
    //std::cout << "Chi-squared/NDF: " << gaussFit->GetChisquare()/gaussFit->GetNDF() << std::endl;
    //std::cout << "Relative Width: " << gaussFit->GetParameter(2)* 100.0f / gaussFit->GetParameter(1) << std::endl;
    //FitAndGetParams(hProjZ, bestMinMass, bestMaxMass); // This will also draw the fit
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

    // Call the function to optimize the fit range (assuming it's defined as before)
    // You might need to adjust the bin ranges according to your histogram
    int xBinStart = 1, xBinEnd = h3->GetXaxis()->GetNbins();
    int yBinStart = 1, yBinEnd =1;//yBinEnd = h3->GetYaxis()->GetNbins();
    OptimizeFitRange(h3, xBinStart, xBinEnd, yBinStart, yBinEnd);

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

/*
void DrawBestHistogram( float leftmost_gauslimit, float rightmost_gauslimit) {
    // more thorough minimizer for fit
    //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    // Set the global fit strategy
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);

    
    // Open the ROOT file and get the histogram
    //old file->cluster dependent cuts. sig fraction implement that cut. cut on centrality? 
    //TFile *file = new TFile("diClusMass_23726_23746_nomPi0CalibCuts.root");
    //new file. cluster dependent cut removed.
    TFile *file = new TFile("diClusMass_allruns_peripheral.root");
    TH1F *hist = (TH1F*)file->Get("h_InvMass");
    
    // Rebin the histogram to have 'numBins' bins
    // First, calculate the rebin factor assuming the histogram's range is 0 to maxXRange
    int numBins= 120;
    int currentNumBins = hist->GetNbinsX();
    double currentXMax = hist->GetXaxis()->GetXmax();
    int rebinFactor = currentNumBins / numBins;
    if (rebinFactor > 1) { // Only rebin if the factor is greater than 1
        std::cout << "current nbins: " << currentNumBins <<" requested nbins: " << numBins << " rebin by: " << rebinFactor << std::endl;
        hist->Rebin(rebinFactor);
        std::cout << "new nbin check: " << hist->GetNbinsX() << std::endl;
    }

    // overall limits
    float rightmost_limit= 0.3;// fit range limit
    float leftmost_limit= 0.05; // fit range limit. normally 0.05
    //float rightmost_limit= 0.257;// fit range limit
    //float leftmost_limit= 0.07; // fit range limit
    // limits on gauss and poly
    float leftpolylim = 0.11;
    float rightpolylim = 0.19;
    hist->GetXaxis()->SetRangeUser(0, 0.4);
    //Double_t scale_factor = 2.5; // Replace with the factor by which you want to scale the errors
    //Double_t error_replace= 0.1;
    scale_histogram_errors(hist, scale_factor);
    //replace_histogram_errors(hist, error_replace);

    // Fit left and right regions with a polynomial, excluding Gaussian region
    TF1 *leftRightFit = new TF1("leftRightFit", leftRightPolynomial, leftmost_limit, rightmost_limit, 5);
    hist->Fit(leftRightFit, "R");

    // Fit Gaussian in the specified range
    TF1 *gausFit = new TF1("gausFit", "gaus", leftpolylim, rightpolylim);//leftpolylim, rightpolylim
    hist->Fit(gausFit, "R");


    // Combined Gaussian + Polynomial fit
    TF1 *combinedFit = new TF1("combinedFit", combinedFunction, leftmost_limit, rightmost_limit, 8);
    // Set initial parameters from previous fits
    for (int i = 0; i < 3; ++i) combinedFit->SetParameter(i, gausFit->GetParameter(i));
    for (int i = 3; i < 8; ++i) combinedFit->SetParameter(i, leftRightFit->GetParameter(i-3));
    //try to improve the fit.
    hist->Fit(combinedFit, "RL");//M
    double chi2 = combinedFit->GetChisquare();
    double ndf = combinedFit->GetNDF();
    double chi2ndf = chi2 / ndf;

    std::cout << "Chi-squared: " << chi2 << std::endl;
    std::cout << "Number of Degrees of Freedom: " << ndf << std::endl;
    std::cout << "Chi-squared/NDF: " << chi2ndf << std::endl;

    //-------------------------------------------show the poly4 part seperately


    // Create a new function for just the polynomial part
    TF1 *polyPart = new TF1("polyPart", "pol4", leftmost_limit, rightmost_limit);

    // Set the parameters of polyPart to those from the combined fit
    // Assuming the first 5 parameters of combinedFit are for the polynomial
    for (int i = 0; i < 5; ++i) polyPart->SetParameter(i, combinedFit->GetParameter(i+3));
    
    
    // Create a new histogram to store the subtracted data
    TH1F *histSubtracted = (TH1F*)hist->Clone("histSubtracted");

    // Subtract the polynomial part
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
        double x = hist->GetBinCenter(i);
        double y = hist->GetBinContent(i) - polyPart->Eval(x);
        histSubtracted->SetBinContent(i, y);
    }
    TF1 *gausFit2 = new TF1("gausFit2", "gaus", leftmost_gauslimit, rightmost_gauslimit);//leftmost_limit, 0.25
    for (int i = 0; i < 3; ++i) gausFit2->SetParameter(i, combinedFit->GetParameter(i));
    histSubtracted->Fit(gausFit2, "R");
    double chi2_s = gausFit2->GetChisquare();
    double ndf_s = gausFit2->GetNDF();
    double chi2ndf_s = chi2_s / ndf_s;

    std::cout << "Chi-squared: " << chi2_s << std::endl;
    std::cout << "Number of Degrees of Freedom: " << ndf_s  << std::endl;
    std::cout << "Chi-squared/NDF: " << chi2ndf_s  << std::endl;
    std::cout << "Relative Width: " << 100*gausFit2->GetParameter(2)/gausFit2->GetParameter(1)  << " %" << std::endl;

    //store 2 separate functions for visualization
    TF1 *fleft = new TF1("fleft",leftRightPolynomial, leftmost_limit, leftpolylim, 5);
    fleft->SetParameters(leftRightFit->GetParameters());
    //hist->GetListOfFunctions()->Add(fleft);
    //gROOT->GetListOfFunctions()->Remove(fleft);
    TF1 *fright = new TF1("fright",leftRightPolynomial, rightpolylim, rightmost_limit, 5);
    fright->SetParameters(leftRightFit->GetParameters());
    //hist->GetListOfFunctions()->Add(fright);
    //gROOT->GetListOfFunctions()->Remove(fright);



    // Draw everything
    //-------------------------------------------------------------------------------------------canvas 1
    TCanvas *c1 = new TCanvas("c1", "Fits", 800, 600);
    hist->Draw("E");

    //gausFit->SetLineColor(kRed);
    //gausFit->Draw("SAME");// draw the gaussian fit
    polyPart->SetLineColor(kRed);
    polyPart->Draw("SAME");


    //leftRightFit->SetLineColor(kBlue);
    fleft->SetLineColor(kBlue);
    fright->SetLineColor(kBlue);

    fleft->Draw("SAME");
    //fright->Draw("SAME");// turn off to see the inflection better.
    //leftRightFit->Draw("SAME"); // Draw the left and right polynomial fits

    combinedFit->SetLineColor(kBlack);
    combinedFit->Draw("SAME"); // draw the combined fit

    // Add a legend
    TLegend *leg = new TLegend(0.1, 0.7, 0.3, 0.9);//bot left x, bot left y, top right x, top right y 
    leg->AddEntry(gausFit, "Gaussian Fit");
    //leg->AddEntry(leftRightFit, "Left & Right Polynomial Fit");
    leg->AddEntry(fleft, "Left & Right Polynomial Fit");
    leg->AddEntry(combinedFit, "Combined Fit");
    leg->Draw();

    // Save the canvas
    c1->SaveAs("combined_fits.pdf");
    //-------------------------------------------------------------------------------------------canvas 2
    TCanvas *c2 = new TCanvas("c2", "Subtracted Peak", 800, 600);
    histSubtracted->Draw();
    histSubtracted->SetMinimum(0.0);
    c2->SaveAs("Subtracted_Peak.pdf");

    //append // Append fit parameters to text file
    appendtextfile(combinedFit, "Combined Fit",scale_factor);
    appendtextfile(gausFit2, "subpgaus fit",scale_factor);

// Second canvas: Custom list of fit results
    TCanvas* c3 = new TCanvas("canvas2", "Fit Parameters", 800, 600);
    c3->cd();

    TPaveText* pt = new TPaveText(0.1, 0.1, 0.9, 0.9, "blNDC"); // blNDC: borderless, normalized coordinates
    pt->SetTextAlign(12); // Align text to the left
    pt->SetFillColor(0); // Transparent background

    // Adding custom text entries
    pt->AddText("Data Fit");
    pt->AddText("Fit Parameters:");
    pt->AddText(Form("Combined Fit Range = %f to %f", leftmost_limit, rightmost_limit));
    pt->AddText(Form("Peak Mean = %f +/- %f", combinedFit->GetParameter(1), combinedFit->GetParError(1)));
    pt->AddText(Form("Peak Sigma = %f +/- %f", combinedFit->GetParameter(2), combinedFit->GetParError(2)));
    pt->AddText(Form("Background Subtracted Peak Fit = %f to %f", leftmost_gauslimit, rightmost_gauslimit));
    pt->AddText(Form("Mean = %f +/- %f", gausFit2->GetParameter(1), gausFit2->GetParError(1)));
    pt->AddText(Form("Sigma = %f +/- %f", gausFit2->GetParameter(2), gausFit2->GetParError(2)));
    pt->AddText(Form("Relative Width: %f",gausFit2->GetParameter(2)* 100.0f / gausFit2->GetParameter(1)));   
    pt->AddText(Form("Chi2/NDF = %f / %d= %f", gausFit2->GetChisquare(), gausFit2->GetNDF(),gausFit2->GetChisquare()/gausFit2->GetNDF()));

    pt->Draw();
    c3->SaveAs("FitInfo.pdf");



    delete file;
    delete c1;
    delete c2;
    delete c3;
    delete gausFit;
    delete gausFit2;
    delete polyPart;
    delete leftRightFit;
    delete fleft;
    delete fright;
    delete combinedFit;
    delete leg;
    
}*/
