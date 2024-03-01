#include <fstream>
#include <iostream>
#include <vector>
#include <sstream>
#include <TH3.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TFile.h>

struct FitConfig {
    int xBinStart, xBinEnd, yBinStart, yBinEnd;
    double minInvMass, maxInvMass;

    // Serialize the configuration for CSV output
    std::string serialize() const {
        std::stringstream ss;
        ss << xBinStart << "," << xBinEnd << ","
           << yBinStart << "," << yBinEnd << ","
           << minInvMass << "," << maxInvMass;
        return ss.str();
    }

    // Deserialize the configuration from CSV input
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
void AnalyzeAndFit(TH3* h3, const std::string& configFilePath, const std::string& outputPDF, const std::string& outputCSV);

void TH3Fitter(){
const char* sourcehistfile="pioncode/rootfiles/data/pt_nclus_differential_data/pt05pt05.root";
const char* histname= "h_pipT_Nclus_mass";
}

std::vector<FitConfig> LoadOrCreateConfig(const std::string& configFile) {
    std::vector<FitConfig> configs;
    std::ifstream file(configFile);
    std::string line;

    if (file) {
        // Load configurations
        while (getline(file, line)) {
            configs.push_back(FitConfig::deserialize(line));
        }
    } else {
        // Create default configuration
        // Example: Assuming 10 bins in X and Y, and invariant mass range [0.1, 1.0]
        for (int xBinStart = 1; xBinStart <= 10; xBinStart += 2) {
            for (int yBinStart = 1; yBinStart <= 10; yBinStart += 2) {
                configs.emplace_back(xBinStart, xBinStart + 1, yBinStart, yBinStart + 1, 0.1, 1.0);
            }
        }

        // Save to file
        std::ofstream outFile(configFile);
        for (const auto& config : configs) {
            outFile << config.serialize() << std::endl;
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

    // Extract and print fit parameters
    double mean = fitFunc->GetParameter(1);
    double sigma = fitFunc->GetParameter(2);
    double chi2 = fitFunc->GetChisquare();
    double ndf = fitFunc->GetNDF();
    double resolution = 100 * sigma / mean;

    // Append fit results to CSV
    csvFile << config.xBinStart << "," << config.xBinEnd << ","
            << config.yBinStart << "," << config.yBinEnd << ","
            << mean << "," << sigma << "," << resolution << "," << chi2/ndf << std::endl;

    // Draw and save to PDF
    TCanvas c;
    hProjZ->Draw();
    c.Print(outputPDF.c_str(), "pdf");
}

// Main analysis function
void AnalyzeAndFit(TH3* h3, const std::string& configFilePath, const std::string& outputPDF, const std::string& outputCSV) {
    auto configs = LoadOrCreateConfig(configFilePath);

    // Open CSV file for output
    std::ofstream csvFile(outputCSV);
    csvFile << "XBinStart,XBinEnd,YBinStart,YBinEnd,Mean,Sigma,Resolution,Chi2/NDF\n"; // Header

    // Initialize PDF
    TCanvas* c = new TCanvas();
    c->Print((outputPDF + "[").c_str()); // Open PDF

    for (const auto& config : configs) {
        FitAndSaveProjection(h3, config, outputPDF, csvFile);
    }

    c->Print((outputPDF + "]").c_str()); // Close PDF
    delete c; // Clean up
}


