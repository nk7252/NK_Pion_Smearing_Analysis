#include <iostream>
#include <vector>
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
    TCanvas* canvas1 = new TCanvas("canvas1", "Overlay Means", 800, 600);
    //canvas1->SetGrid();
    gStyle->SetOptStat(0);

    // Create a legend
    //TLegend* legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);

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
        if (!h18) {
            std::cerr << "Error: Histogram 'myHist2D' not found in file " << fileNames[i] << std::endl;
            pionfile->Close();
            continue;
        }

        // Create a histogram for means
        TH1F* meanHistogram = new TH1F(Form("MeanHistogram_%zu", i), Form("Version %zu", i), h18D->GetNbinsX(), 0.5, h18->GetNbinsX() + 0.5);



        // Loop over each bin in the X direction
        for (int binX = 1; binX <= h18->GetNbinsX(); ++binX) {
            // Project along Y for each binX
            TH1D* yProjection = h18->ProjectionY(Form("YProjection_%zu_%d", i, binX), binX, binX, "");

            // Fit the Y projection with a Gaussian
            yProjection->Fit("gaus", "Q");

            // Access the fit parameters
            TF1* fitFunc = yProjection->GetFunction("gaus");

            // Check if the fit function is valid
            if (fitFunc) {
                // Fill the mean histogram with the mean value
                meanHistogram->SetBinContent(binX, fitFunc->GetParameter(1));
                meanHistogram->SetBinError(binX, fitFunc->GetParError(1));
            }
                // Add an entry to the legend
                //legend1->AddEntry(yProjection, Form("Version %zu, BinX %d", i, binX), "L");
            
            std::cout << "I reached here, pre delete proj" << std::endl; // debug line
            // Clean up Y projection
            delete yProjection;
        }
        // Set different line colors for each version
        int lineColor = i + 1; // Line color: 1, 2, 3, ...
        meanHistogram->SetLineColor(lineColor);
        std::cout << "I reached here, done with loop over bins" << std::endl; // debug line

        // Overlay the mean histogram on the same canvas
        if (i == 0) {
            meanHistogram->Draw("E"); // Draw histogram for the first version
        } else {
            meanHistogram->Draw("E SAME"); // Draw subsequent histograms on the same canvas
        }

        // Add an entry to the legend
        //legend->AddEntry(meanHistogram, Form("Version %zu", i), "L");

        // Close the file
        pionfile->Close();
        delete pionfile;

        std::cout << "I reached here, close+delete file" << std::endl; // debug line
    }

    std::cout << "I reached here, done with all files" << std::endl; // debug line
    // Draw the legend
    //legend1->Draw();

    // Show the canvas
    canvas1->Update();
    canvas1->Modified();
    //canvas1->Print("OverlayMeanHistograms.pdf");

    // Clean up
    //delete canvas1;
    //delete legend1;
}

void CombinedFits() {
    // List of root file names
    std::vector<std::string> fileNames = {"pioncode/Pi0FastMC_0.155000WSHP.root"};
    // {"pioncode/Pi0FastMC_0.155000EXP.root", "pioncode/Pi0FastMC_0.155000POWER.root", "pioncode/Pi0FastMC_0.15500WSHP.root"};
    // Overlay the means for each bin
    OverlayMeans(fileNames);

    //return 0;
}
