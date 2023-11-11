#include <iostream>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>

void PlotMultipleVersions(const std::vector<std::string>& fileNames) {
    // Create a TCanvas
    TCanvas* canvas = new TCanvas("canvas", "Multiple Versions Plot", 800, 600);
    canvas->SetGrid();
    gStyle->SetOptStat(0);

    // Create a legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    // Loop over each file
    for (size_t i = 0; i < fileNames.size(); ++i) {
        // Open the root file
        TFile* file = new TFile(fileNames[i].c_str(), "READ");

        // Check if the file is open
        if (!file || file->IsZombie()) {
            std::cerr << "Error: Unable to open file " << fileNames[i] << std::endl;
            continue;
        }

        // Retrieve the histogram
        TH1F* hist = nullptr;
        file->GetObject("myHist", hist);

        // Check if the histogram exists
        if (!hist) {
            std::cerr << "Error: Histogram 'myHist' not found in file " << fileNames[i] << std::endl;
            file->Close();
            continue;
        }

        // Set different marker styles for each version
        int markerStyle = 20 + i; // Marker style: 20, 21, 22, ...
        hist->SetMarkerStyle(markerStyle);

        // Plot the histogram
        if (i == 0) {
            hist->Draw("P"); // Draw with markers for the first version
        } else {
            hist->Draw("P SAME"); // Draw with markers and on the same canvas for subsequent versions
        }

        // Add an entry to the legend
        legend->AddEntry(hist, Form("Version %zu", i), "P");

        // Close the file
        file->Close();
    }

    // Draw the legend
    legend->Draw();

    // Show the canvas
    canvas->Update();
    canvas->Modified();
    canvas->Print("MultipleVersionsPlot.pdf");
    //canvas->SaveAs(Form("pioncode/canvas_pdf/Alt_Projection_%s.pdf", pdfname));
		
    // Clean up
    delete canvas;
    delete legend;
}

int main() {
    // List of root file names
    std::vector<std::string> fileNames = {"pioncode/Pi0FastMC_0.155000EXP.root", "pioncode/Pi0FastMC_0.155000POWER.root", "pioncode/Pi0FastMC_0.155000WSHP.root"};

    // Plot the histograms from different versions
    PlotMultipleVersions(fileNames);

    return 0;
}