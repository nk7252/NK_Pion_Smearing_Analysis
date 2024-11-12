#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <iostream>

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
    // int nBinsX = h_res_ptTr->GetNbinsX();

    int startBin = 1;
    int endBin_global = -1; //-1=actual last bin
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

    // Open multi-page PDF for viewing all pT bin fits
    TCanvas *c_summary = new TCanvas("c_summary", "Summary of Fits", 800, 600);
    c_summary->Print("pioncode/canvas_pdf/e_res_fit_monitoring.pdf[");

    // Loop over each pT bin and perform projection and fitting
    // for (int i = 1; i <= nBinsX; ++i) {
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
        Double_t x2 = h_res_ptTr->GetXaxis()->GetBinUpEdge(i);
        Double_t pTValue = h_res_ptTr->GetXaxis()->GetBinCenter(i);

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

        // Fill the mean and sigma histograms
        h_mean->SetBinContent(i, mean);
        h_mean->SetBinError(i, meanError);
        h_sigma->SetBinContent(i, sigma);
        h_sigma->SetBinError(i, sigmaError);

        // Store chi² and chi²/ndf values
        h_chi2->SetBinContent(i, chi2);
        h_chi2_ndf->SetBinContent(i, chi2ndf);

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

    // Close the multi-page PDF for pT bin fits
    c_summary->Print("pioncode/canvas_pdf/e_res_fit_monitoring.pdf]");

    // Calculate the relative width (sigma divided by mean)
    TH1D *h_relwidth = (TH1D *)h_sigma->Clone("h_relwidth");
    h_relwidth->Divide(h_mean);

    // Final fit over the entire sigma/mean range
    TCanvas *c_result = new TCanvas("c_result", "Energy Resolution Fit", 800, 600);
    c_result->Print("pioncode/canvas_pdf/energy_resolution_fit_results.pdf["); // Open PDF
    TF1 *fitFunc = new TF1("fitFunc", "sqrt(([0]*[0])/sqrt(x*x + 0.1349768*0.1349768) + [1]*[1])", 0.6, 9);
    fitFunc->SetParNames("sqrt(E) term", "Constant term");
    // fitFunc->SetParameters(0.15, 0.06);
    fitFunc->SetParLimits(0, 0.1, 0.2);
    fitFunc->SetParLimits(1, 0.02, 0.08);

    h_relwidth->Fit(fitFunc, "RE");

    // Draw and save the final fit
    h_relwidth->SetTitle("Relative Width vs. p_{T}");
    h_relwidth->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    h_relwidth->GetYaxis()->SetTitle("Relative Width (#sigma/#mu)");
    h_relwidth->Draw();
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
    h_mean->Draw("E");

    // Plot sigma histogram
    c_summary_page->cd(2);
    h_sigma->SetTitle("Sigma of Energy Resolution");
    h_sigma->GetXaxis()->SetTitle("p_{T} (GeV/c)");
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
    file->Close();
    gApplication->Terminate(0);
}
