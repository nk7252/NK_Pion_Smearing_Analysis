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



//declarations
std::vector<FitConfig> LoadOrCreateConfig(const std::string& configFile);
void FitAndSaveProjection(TH3* h3, const FitConfig& config, const std::string& outputPDF, std::ofstream& csvFile);
void AnalyzeAndFit(const std::string& rootFileName, const std::string& histName);
double FitAndGetChi2NDF(TH1D* hProjZ, double minMass, double maxMass);
void OptimizeFitRange(TH3* h3, int xBinStart, int xBinEnd, int yBinStart, int yBinEnd);
void OptimizeHistogramFit(const std::string& rootFileName, const std::string& histogramName);


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
double FitAndGetChi2NDF(TH1D* hProjZ, double minMass, double maxMass) {
    TF1* fitFunc = new TF1("fitFunc", "gaus", minMass, maxMass);
    hProjZ->Fit(fitFunc, "RQ"); // "RQ" option for Range and Quiet

    double chi2 = fitFunc->GetChisquare();
    double ndf = fitFunc->GetNDF();
    double chi2ndf = (ndf > 0) ? chi2 / ndf : -1; // Avoid division by zero

    delete fitFunc; // Clean up
    return chi2ndf;
}

// Main function to explore fit ranges and find optimal parameters
void OptimizeFitRange(TH3* h3, int xBinStart, int xBinEnd, int yBinStart, int yBinEnd) {
    TH1D* hProjZ = h3->ProjectionZ("projZ", xBinStart, xBinEnd, yBinStart, yBinEnd);

    double bestChi2NDF = TMath::Infinity();
    double bestMinMass = 0, bestMaxMass = 0;

    double startMass = 0.09; // start of mass range
    double endMass = 0.3; // end of mass range
    double step = 0.01; // step size

    for (double minMass = startMass; minMass < endMass; minMass += step) {
        for (double maxMass = minMass + step; maxMass <= endMass; maxMass += step) {
            double chi2NDF = FitAndGetChi2NDF(hProjZ, minMass, maxMass);

            if (chi2NDF >= 0.9 && chi2NDF < bestChi2NDF) {
                bestChi2NDF = chi2NDF;
                bestMinMass = minMass;
                bestMaxMass = maxMass;
            }
        }
    }


    // Output the best fit range and chi^2/NDF
    std::cout << "Best Fit Range: [" << bestMinMass << ", " << bestMaxMass << "] with chi^2/NDF = " << bestChi2NDF << std::endl;

    // Optionally, perform and visualize the final fit with the best parameters
    FitAndGetChi2NDF(hProjZ, bestMinMass, bestMaxMass); // This will also draw the fit
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




