#include <TFile.h>
#include <TH1.h>
#include <TObject.h>
#include <iostream>

void detach_hist() {
    // Open the input ROOT file
    TFile *inputFile = TFile::Open("pioncode/rootfiles/OUTHIST_iter_G4Hits_single_eta_p_600_20000MeV_0000000017_00merged_V14.root", "READ");
    
    //pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_single_eta_p_600_20000MeV_0000000017_00merged_V47.root
    //pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_single_pi0_p_200_20000MeV_0000000017_00merged_V53.root
    //pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_single_pi0_p_200_20000MeV_0000000017_00merged_V53.root
    
    if (!inputFile || inputFile->IsZombie()) {
        std::cerr << "Error: Could not open input file." << std::endl;
    }

    // Get the specific histogram by name
    TH1 *hist = dynamic_cast<TH1*>(inputFile->Get("h_truth_spectrum2"));
    if (!hist) {
        std::cerr << "Error: Histogram not found in the input file." << std::endl;
        inputFile->Close();
    }

    // Clone the histogram to detach it from the input file
    TH1 *histCopy = dynamic_cast<TH1*>(hist->Clone());
    if (!histCopy) {
        std::cerr << "Error: Could not clone the histogram." << std::endl;
        inputFile->Close();
    }

    // Give the cloned histogram a new name
    histCopy->SetName("seta_pt_spectrum");

    // Open the output ROOT file to save the cloned histogram
    TFile *outputFile = TFile::Open("pioncode/rootfiles/seta_spectrum.root", "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error: Could not create output file." << std::endl;
        inputFile->Close();
    }

    // Write the cloned histogram to the new file
    outputFile->cd();
    histCopy->Write();

    // Clean up
    inputFile->Close();
    outputFile->Close();

    std::cout << "Histogram successfully copied to output.root" << std::endl;
}