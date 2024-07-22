#include <TFile.h>
#include <TH1F.h>

void CopyHistogram(const char* inputFile, const char* outputFileName, const char* histogramName) {
    // Open the input ROOT file
    TFile* inFile = new TFile(inputFile, "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Error: Could not open input file " << inputFile << std::endl;
        return;
    }

    // Get the histogram from the input file
    TH1F* inputHist = dynamic_cast<TH1F*>(inFile->Get(histogramName));
    if (!inputHist) {
        std::cerr << "Error: Could not retrieve histogram " << histogramName << " from " << inputFile << std::endl;
        inFile->Close();
        return;
    }

    // Create a new output ROOT file
    TFile* outFile = new TFile(Form("../pioncode/rootfiles/%s",outputFileName), "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cerr << "Error: Could not create output file " << outputFileName << std::endl;
        inFile->Close();
        return;
    }

    // Clone the input histogram and write it to the output file
    TH1F* outputHist = dynamic_cast<TH1F*>(inputHist->Clone());
    outputHist->SetDirectory(0);  // Detach histogram from the file
    outputHist->Write();

    // Close the input and output files
    inFile->Close();
    outFile->Close();

    std::cout << "Histogram '" << histogramName << "' copied from '" << inputFile << "' to '" << outputFileName << "'" << std::endl;
}

int photonspectrumhist() {
    // Example usage
    const char* inputFile = "../pioncode/rootfiles/Pi0FastMC_0.155000_WSHP.root";           // Replace with your input ROOT file
    const char* outputFileName = "Photon_spectrum_hist.root";
    const char* histogramName = "h16";    // Replace with the name of the histogram you want to copy

    CopyHistogram(inputFile, outputFileName, histogramName);

    return 0;
}