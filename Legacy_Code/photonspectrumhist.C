#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TRandom3.h>

int main() {
    // Create a histogram
    TH1F *histogram = new TH1F("myhist", "My Histogram", 100, 0, 10);

    // Create a random number generator
    TRandom3 randomGenerator;

    // Fill the histogram with random data
    for (int i = 0; i < 10000; ++i) {
        double value = randomGenerator.Gaus(5, 1);  // Example: Gaussian distribution with mean=5, sigma=1
        histogram->Fill(value);
    }

    // Save the histogram in a ROOT file
    TFile *file = new TFile("output.root", "RECREATE");
    histogram->Write();

    // Write and close the file
    file->Write();
    file->Close();

    // Cleanup
    delete histogram;
    delete file;

    return 0;
}
