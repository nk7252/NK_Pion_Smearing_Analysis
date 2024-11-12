#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TGraphErrors.h>
#include <Math/MinimizerOptions.h>
#include <iostream>
#include <vector> // Include the vector library
#include <algorithm> // For std::min

// local includes
#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

// Main function to perform the fits and generate PDFs
void eresFit()
{
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
    ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(10000);
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.001);
    ROOT::Math::MinimizerOptions::SetDefaultPrecision(1e-12);
    // ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(3);
    SetsPhenixStyle();

    // Open the ROOT file containing your histogram
    TFile *file = TFile::Open("pioncode/rootfiles/OUTHIST_iter_G4Hits_single_eta_p_600_20000MeV_0000000017_00merged_V21.root");
    if (!file || file->IsZombie())
    {
        std::cout << "Error opening file!" << std::endl;
        return;
    }

    // Retrieve the 2D histogram (e.g., resolution vs. pT)
    TH2D *h_res_ptTr = (TH2D *)file->Get("h_res_ptTr");
    if (!h_res_ptTr)
    {
        std::cout << "Histogram h_res_ptTr not found!" << std::endl;
        file->Close();
        return;
    }

    // Create histograms to store mean and sigma values
    int nBinsX = h_res_ptTr->GetNbinsX();
    double xMin = h_res_ptTr->GetXaxis()->GetXmin();
    double xMax = h_res_ptTr->GetXaxis()->GetXmax();

    int startBin = 1;
    int endBin_global = -1; // -1=actual last bin
    int endBin = endBin_global;
    int projectionBins = 4;
    double scale_factor = 1.0;
    if (endBin_global == -1)
        endBin = nBinsX; // Default to the last bin if not specified

    TH1D *h_mean = new TH1D("h_mean", "Mean of Energy Resolution", nBinsX, xMin, xMax);
    TH1D *h_sigma = new TH1D("h_sigma", "Sigma of Energy Resolution", nBinsX, xMin, xMax);

    // Create histograms to store chi² and chi²/ndf values
    TH1D *h_chi2 = new TH1D("h_chi2", "Chi^{2} of Individual Fits", nBinsX, xMin, xMax);
    TH1D *h_chi2_ndf = new TH1D("h_chi2_ndf", "Chi^{2}/NDF of Individual Fits", nBinsX, xMin, xMax);

    // Create TGraphErrors to store the relative width (sigma/mean) with proper error propagation
    TGraphErrors *Eres_Graph = new TGraphErrors();
    int graphPoint = 0; // Index for TGraphErrors points

    // Vectors to store pT ranges, chi2/ndf values, and fit ranges for the table
    std::vector<double> pT_low_edges;
    std::vector<double> pT_high_edges;
    std::vector<double> chi2ndf_values;
    std::vector<double> fitMin_values;
    std::vector<double> fitMax_values;

    // Open multi-page PDF for viewing all pT bin fits
    TCanvas *c_summary = new TCanvas("c_summary", "Summary of Fits", 800, 600);
    c_summary->Print("pioncode/canvas_pdf/e_res_fit_monitoring.pdf[");

    // Loop over each pT bin and perform projection and fitting
    for (int i = startBin; i <= endBin; i += projectionBins)
    {
        // Project the histogram along the Y-axis
        int lastBin = std::min(i + projectionBins - 1, nBinsX);
        // Project the 2D histogram onto the Y-axis (resolution) for the current pT bin
        TH1D *h_proj = h_res_ptTr->ProjectionY(Form("h_proj_%d", i), i, lastBin);

        // Skip empty projections
        if (h_proj->GetEntries() == 0)
        {
            delete h_proj;
            continue;
        }
        h_proj->Scale(1, "width");

        // Get the pT bin limits and center for dynamic fit range
        Double_t x1 = h_res_ptTr->GetXaxis()->GetBinLowEdge(i);
        Double_t x2 = h_res_ptTr->GetXaxis()->GetBinUpEdge(lastBin);
        Double_t pTValue = (x1 + x2) / 2.0;
        Double_t pTError = (x2 - x1) / 2.0;

        // Set the fit range for gausFit dynamically based on the pT bin value
        Double_t fitMin, fitMax;

        if (projectionBins == 1)
        {
            //*
            if (pTValue < 0.25)
            {
                fitMin = 0.7;
                fitMax = 1.07;
            }
            else if (pTValue < 0.5)
            {
                fitMin = 0.60;
                fitMax = 1.07;
            }
            else if (pTValue < 1)
            {
                fitMin = 0.75;
                fitMax = 1.07;
            }
            else if (pTValue < 5)
            {
                fitMin = 0.89;
                fitMax = 1.07;
            }
            else if (pTValue < 10)
            {
                fitMin = 0.89;
                fitMax = 1.07;
            }
            else if (pTValue < 14)
            {
                fitMin = 0.89;
                fitMax = 1.07;
            }
            else if (pTValue < 16.2)
            {
                fitMin = 0.89;
                fitMax = 1.05;
            }
            else if (pTValue < 18)
            {
                fitMin = 0.86;
                fitMax = 1.02;
            }
            else
            {
                continue;
                // fitMin = 0.90; fitMax = 1.05;
            }
            //*/
        }
        else if (projectionBins == 2)
        {
            if (pTValue < 0.25)
            {
                fitMin = 0.7;
                fitMax = 1.07;
            }
            else if (pTValue < 0.5)
            {
                fitMin = 0.60;
                fitMax = 1.07;
            }
            else if (pTValue < 1)
            {
                fitMin = 0.75;
                fitMax = 1.07;
            }
            else if (pTValue < 5)
            {
                fitMin = 0.89;
                fitMax = 1.07;
            }
            else if (pTValue < 10)
            {
                fitMin = 0.89;
                fitMax = 1.07;
            }
            else if (pTValue < 14)
            {
                fitMin = 0.89;
                fitMax = 1.07;
            }
            else if (pTValue < 16.2)
            {
                fitMin = 0.89;
                fitMax = 1.05;
            }
            else if (pTValue < 18)
            {
                fitMin = 0.86;
                fitMax = 1.02;
            }
            else
            {
                continue;
                // fitMin = 0.90; fitMax = 1.05;
            }
        }
        else if (projectionBins == 4)
        {
            if (pTValue < 0.25)
            {
                fitMin = 0.7;
                fitMax = 1.07;
            }
            else if (pTValue < 0.5)
            {
                fitMin = 0.60;
                fitMax = 1.07;
            }
            else if (pTValue < 1)
            {
                fitMin = 0.75;
                fitMax = 1.07;
            }
            else if (pTValue < 5)
            {
                fitMin = 0.89;
                fitMax = 1.07;
            }
            else if (pTValue < 10)
            {
                fitMin = 0.89;
                fitMax = 1.07;
            }
            else if (pTValue < 14)
            {
                fitMin = 0.89;
                fitMax = 1.07;
            }
            else if (pTValue < 16.2)
            {
                fitMin = 0.89;
                fitMax = 1.05;
            }
            else if (pTValue < 18)
            {
                fitMin = 0.86;
                fitMax = 1.02;
            }
            else
            {
                continue;
                // fitMin = 0.90; fitMax = 1.05;
            }
        }

        // Define the Gaussian fit function with the dynamic range
        TF1 *gausFit = new TF1("gausFit", "gaus", fitMin, fitMax);

        // Perform the fit on the projection
        h_proj->Fit(gausFit, "REQ"); // "R" for range, "Q" for quiet mode

        // Extract fit parameters
        Double_t mean = gausFit->GetParameter(1);
        Double_t sigma = gausFit->GetParameter(2);
        Double_t meanError = gausFit->GetParError(1);
        Double_t sigmaError = gausFit->GetParError(2);
        Double_t chi2 = gausFit->GetChisquare();
        Double_t ndf = gausFit->GetNDF();
        Double_t chi2ndf = chi2 / ndf;

        // Calculate the relative width and its error
        Double_t width = sigma / mean;
        Double_t widthError = width * sqrt(pow(meanError / mean, 2) + pow(sigmaError / sigma, 2));

        // Fill the mean and sigma histograms
        h_mean->SetBinContent(i, mean);
        h_mean->SetBinError(i, meanError);
        h_sigma->SetBinContent(i, sigma);
        h_sigma->SetBinError(i, sigmaError);
        std::cout << "Mean: " << mean << " +/- " << meanError << std::endl;
        std::cout << "Sigma: " << sigma << " +/- " << sigmaError << std::endl;
        std::cout << "Bin: " << i << " ,pT range:  ( " << x1 << " , " << x2 << ") ,  pT Bin: " << pTValue << " +/- " << pTError << std::endl;

        // Store chi² and chi²/ndf values
        h_chi2->SetBinContent(i, chi2);
        h_chi2_ndf->SetBinContent(i, chi2ndf);

        // **Store pT ranges, chi²/ndf values, and fit ranges for the table**
        pT_low_edges.push_back(x1);
        pT_high_edges.push_back(x2);
        chi2ndf_values.push_back(chi2ndf);
        fitMin_values.push_back(fitMin);
        fitMax_values.push_back(fitMax);

        // Add the relative width data point to the TGraphErrors
        Eres_Graph->SetPoint(graphPoint, pTValue, width);
        Eres_Graph->SetPointError(graphPoint, pTError, widthError);
        graphPoint++;

        // Draw and save each bin’s fit in the PDF
        c_summary->cd();
        h_proj->SetTitle(Form("p_{T} Bin: [%.2f, %.2f] GeV/c - Fit", x1, x2));
        h_proj->GetXaxis()->SetTitle("Resolution");
        h_proj->GetYaxis()->SetTitle("Entries");
        h_proj->Draw("PE");
        gausFit->SetLineColor(kRed);
        gausFit->Draw("same");

        // Add text box with bin range
        TPaveText *pave = new TPaveText(0.15, 0.75, 0.35, 0.85, "NDC");
        pave->AddText(Form("Bin Range: %.2f - %.2f GeV/c", x1, x2));
        pave->SetFillColor(0);
        pave->SetTextSize(0.03);
        pave->Draw();

        // Save this fit result to the multi-page PDF
        c_summary->Print("pioncode/canvas_pdf/e_res_fit_monitoring.pdf");

        // Create additional canvas for fit result details
        TCanvas *c_fitDetails = new TCanvas("c_fitDetails", "Fit Details", 800, 600);
        TPaveText *fitResults = new TPaveText(0.1, 0.6, 0.4, 0.9, "NDC");
        fitResults->AddText(Form("Fit Results for p_{T} Bin [%.2f, %.2f] GeV/c", x1, x2));
        fitResults->AddText(Form("Mean = %.4f #pm %.4f", mean, meanError));
        fitResults->AddText(Form("Sigma = %.4f #pm %.4f", sigma, sigmaError));
        fitResults->AddText(Form("Chi^{2} = %.2f", chi2));
        fitResults->AddText(Form("Chi^{2}/NDF = %.2f", chi2ndf));
        fitResults->SetFillColor(0);
        fitResults->SetTextSize(0.04);
        fitResults->Draw();

        // Save the fit results page to the PDF
        c_fitDetails->Print("pioncode/canvas_pdf/e_res_fit_monitoring.pdf");

        // Clean up
        delete pave;
        delete fitResults;
        delete c_fitDetails;
        delete h_proj;
        delete gausFit;
    }

    // **Create a table of Chi^2/NDF values and fit ranges, limited to 10 entries per page**

    int entries_per_page = 10;
    int total_entries = chi2ndf_values.size();
    int num_pages = (total_entries + entries_per_page - 1) / entries_per_page;

    for (int page = 0; page < num_pages; ++page)
    {
        // Clear the c_summary canvas
        c_summary->cd();
        c_summary->Clear();

        // Create the table using TPaveText
        TPaveText *table = new TPaveText(0.05, 0.05, 0.95, 0.95, "NDC");
        table->SetBorderSize(1);
        table->SetFillColor(0);
        table->SetTextAlign(12);
        table->SetTextFont(42);
        table->SetTextSize(0.03);

        // Add a header
        table->AddText("Chi^{2}/NDF and Fit Range Summary Table");
        table->AddText(Form("Page %d", page + 1));
        table->AddText(" ");

        // Add column headers with fixed-width formatting
        TString header = Form("%-20s %-15s %-20s", "p_{T} Range (GeV/c)", "Chi^{2}/NDF", "Fit Range");
        table->AddText(header);

        // Loop over entries in this page
        for (int i = page * entries_per_page; i < std::min((page + 1) * entries_per_page, total_entries); ++i)
        {
            TString line = Form("%-20s %-15.2f %-20s", 
                                Form("%.2f - %.2f", pT_low_edges[i], pT_high_edges[i]), 
                                chi2ndf_values[i], 
                                Form("%.2f - %.2f", fitMin_values[i], fitMax_values[i]));
            table->AddText(line);
        }

        // Draw the table
        table->Draw();

        // Save the table to the PDF
        c_summary->Print("pioncode/canvas_pdf/e_res_fit_monitoring.pdf");

        // Clean up
        delete table;
    }

    // Close the multi-page PDF for pT bin fits
    c_summary->Print("pioncode/canvas_pdf/e_res_fit_monitoring.pdf]");

    // Final fit over the entire sigma/mean range using TGraphErrors
    TCanvas *c_result = new TCanvas("c_result", "Energy Resolution Fit", 800, 600);
    c_result->Print("pioncode/canvas_pdf/energy_resolution_fit_results.pdf["); // Open PDF

    // Define the fit function
    TF1 *fitFunc = new TF1("fitFunc", "sqrt(([0]*[0])/sqrt(x*x + 0.1349768*0.1349768) + [1]*[1])", 0.6, 9);
    fitFunc->SetParNames("sqrt(E) term", "Constant term");
    // fitFunc->SetParameters(0.15, 0.06);
    fitFunc->SetParLimits(0, 0.1, 0.2);
    fitFunc->SetParLimits(1, 0.02, 0.08);

    // Perform the fit on the TGraphErrors
    Eres_Graph->Fit(fitFunc, "RE");

    // Draw and save the final fit
    Eres_Graph->SetTitle("Relative Width vs. p_{T}");
    Eres_Graph->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    Eres_Graph->GetYaxis()->SetTitle("Relative Width (#sigma/#mu)");
    Eres_Graph->SetMarkerStyle(20);
    Eres_Graph->SetMarkerSize(1);
    Eres_Graph->Draw("AP"); // "A" for axes, "P" for points with error bars
    fitFunc->SetLineColor(kRed);
    fitFunc->Draw("same");

    // Save the final fit canvas to a PDF
    c_result->Print("pioncode/canvas_pdf/energy_resolution_fit_results.pdf");

    // Additional page with original histograms and chi² plots
    TCanvas *c_summary_page = new TCanvas("c_summary_page", "Fit Summary", 1200, 800);
    c_summary_page->Divide(2, 2);

    // Plot mean histogram
    c_summary_page->cd(1);
    h_mean->SetTitle("Mean of Energy Resolution");
    h_mean->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_mean->GetYaxis()->SetTitle("Mean");
    h_mean->Draw("E");

    // Plot sigma histogram
    c_summary_page->cd(2);
    h_sigma->SetTitle("Sigma of Energy Resolution");
    h_sigma->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_sigma->GetYaxis()->SetTitle("Sigma");
    h_sigma->Draw("E");

    // Plot chi² histogram with correct formatting for superscript
    c_summary_page->cd(3);
    h_chi2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_chi2->GetYaxis()->SetTitle("Chi^{2}");
    h_chi2->SetTitle("Chi^{2} of Individual Fits");
    h_chi2->Draw("E");

    // Plot chi²/ndf histogram with correct formatting
    c_summary_page->cd(4);
    h_chi2_ndf->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_chi2_ndf->GetYaxis()->SetTitle("Chi^{2}/NDF");
    h_chi2_ndf->SetTitle("Chi^{2}/NDF of Individual Fits");
    h_chi2_ndf->Draw("E");
    c_summary_page->Print("pioncode/canvas_pdf/energy_resolution_fit_results.pdf");

    // Save all plots to the main PDF
    c_result->Print("pioncode/canvas_pdf/energy_resolution_fit_results.pdf]");

    // Clean up
    delete fitFunc;
    delete c_summary_page;
    delete c_summary;
    delete c_result;
    file->Close();
    gApplication->Terminate(0);
}
