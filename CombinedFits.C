#include <iostream>
#include <TFile.h>
#include <TH2F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TF1.h>

void OverlayMeans(const std::vector<std::string>& fileNames) {
    // Create a TCanvas
    TCanvas* canvas = new TCanvas("canvas", "Overlay Means", 800, 600);
    canvas->SetGrid();
    gStyle->SetOptStat(0);

    // Create a legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);

    // Loop over each file
    for (size_t i = 0; i < fileNames.size(); ++i) {
        // Open the root file
        TFile* pionfile = new TFile(fileNames[i].c_str(), "READ");

        // Check if the file is open
        if (!pionfile || pionfile->IsZombie()) {
            std::cerr << "Error: Unable to open file " << fileNames[i] << std::endl;
            continue;
        }

        // Retrieve the 2D histogram
        TH2F *h18 =(TH2F *)pionfile->Get("h18");
        //TH2F* hist2D = nullptr;
        //file->GetObject("h18", hist2D);

        // Check if the histogram exists
        if (!hist2D) {
            std::cerr << "Error: Histogram 'myHist2D' not found in file " << fileNames[i] << std::endl;
            file->Close();
            continue;
        }

        // Set different line colors for each version
        int lineColor = i + 1; // Line color: 1, 2, 3, ...
        hist2D->SetLineColor(lineColor);

        // Loop over each bin in the X direction
        for (int binX = 1; binX <= hist2D->GetNbinsX(); ++binX) {
            // Project along Y for each binX
            TH1D* yProjection = hist2D->ProjectionY(Form("YProjection_%d_%d", i, binX), binX, binX, "");

            // Fit the Y projection with a Gaussian
            yProjection->Fit("gaus", "Q");

            // Access the fit parameters
            TF1* fitFunc = yProjection->GetFunction("gaus");

            // Check if the fit function is valid
            if (fitFunc) {
                // Overlay the mean on the same canvas
                if (i == 0 && binX == 1) {
                    yProjection->Draw("HIST"); // Draw histogram for the first version and binX
                } else {
                    yProjection->Draw("HIST SAME"); // Draw subsequent histograms on the same canvas
                }

                // Add an entry to the legend
                legend->AddEntry(yProjection, Form("Version %zu, BinX %d", i, binX), "L");
            }

            // Clean up Y projection
            delete yProjection;
        }

        // Close the file
        pionfile->Close();
    }

    // Draw the legend
    legend->Draw();

    // Show the canvas
    canvas->Update();
    canvas->Modified();
    canvas->Print("OverlayMeansPlot.pdf");

    // Clean up
    delete canvas;
    delete legend;
}

int CombinedFits() {
    // List of root file names
    std::vector<std::string> fileNames = {"pioncode/Pi0FastMC_0.155000EXP.root", "pioncode/Pi0FastMC_0.155000POWER.root", "pioncode/Pi0FastMC_0.155000WSHP.root"};

    // Overlay the means for each bin
    OverlayMeans(fileNames);

    return 0;
}
