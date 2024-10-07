#pragma once
// c++ includes
#include <iostream>
#include <vector>
#include <string>
// root includes
#include <TFile.h>
#include <TH2.h>
#include <TF1.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TROOT.h>
#include <Math/MinimizerOptions.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

// local includes
#include "sPhenixStyle.h"
#include "sPhenixStyle.C"


void Spectrum_Fit()
{
    // Set the default minimizer options
    ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
    ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(10000);
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.001);
    ROOT::Math::MinimizerOptions::SetDefaultPrecision(1e-12);

    SetsPhenixStyle();  // Apply sPHENIX style settings

    // Open the ROOT file
    TFile *file = new TFile("pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_pythia8_pp_mb_0000000015_merged_V55.root", "READ");
    if (!file->IsOpen())
    {
        std::cerr << "Error opening file: " << "pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_pythia8_pp_mb_0000000015_merged_V55.root" << std::endl;
        return;
    }

    // Retrieve the histogram
    TH1F *hist = dynamic_cast<TH1F *>(file->Get("h_truth_spectrum2"));
    if (!hist)
    {
        std::cerr << "Error retrieving histogram: " << "h_truth_spectrum2" << " from file" << std::endl;
        file->Close();
        return;
    }

    // Scale and rebin the histogram
    hist->Scale(1/hist->Integral(), "width");
    hist->Rebin(2);

    // Scale bin contents by bin center (pT), avoiding division by zero
    for (int i = 1; i <= hist->GetNbinsX(); i++)
    {
        double binpT = hist->GetBinCenter(i); // Get the bin center (pT)
        if (binpT == 0) continue;             // Avoid division by zero
        double binContent = hist->GetBinContent(i); // Get the bin content
        double newContent = binContent / binpT;     // Scale the bin content
        hist->SetBinContent(i, newContent);         // Set the new bin content
        hist->SetBinError(i, hist->GetBinError(i) / binpT); // Scale the bin error
    }

    // Define fitting ranges and transition point
    int lowedge = 1;
    int highedge = 16;
    float transition = 4.8;
    bool showauxilliaryfits = false;

    // Define the low pT fit function (Hagedorn)
    TF1 *lowPtFunc = new TF1("lowPtFunc", "[0] / pow(1 + x / [1], [2])", lowedge, transition);
    lowPtFunc->SetParameters(1525, 0.64, 7.9);
    lowPtFunc->SetNpx(1000);

    // Define the high pT fit function (Power Law)
    TF1 *highPtFunc = new TF1("highPtFunc", "[0] / (pow(x, [1]))", transition, highedge);
    //highPtFunc->SetParameters(404.4, 10);
    highPtFunc->SetParameters(141.5, 9.9);
    highPtFunc->SetNpx(1000);

    // Fit the histogram with low and high pT functions and get fit results
    TFitResultPtr fitResultPtrL = hist->Fit(lowPtFunc, "SR");
    const ROOT::Fit::FitResult &fitResultL = *fitResultPtrL;

    TFitResultPtr fitResultPtrH = hist->Fit(highPtFunc, "SR");
    const ROOT::Fit::FitResult &fitResultH = *fitResultPtrH;

    // Define the combined fit function
    TF1  *myFunc = new TF1("myFunc", "((1 / (1 + exp((x - [0]) / [1]))) * [2] / pow(1 + x / [3], [4]) + (1 - (1 / (1 + exp((x - [0]) / [1])))) * [5] / (pow(x, [6])))", lowedge, highedge);

    // Initialize combined function parameters
    myFunc->SetParameter(0, transition);
    myFunc->SetParameter(1, 0.114);
    myFunc->SetParLimits(1, 0.114*0.8,  0.114*1.6);
    for (int j=0; j<3; j++) myFunc->SetParameter(j+2, lowPtFunc->GetParameter(j));
    for (int j=0; j<2; j++) myFunc->SetParameter(j+5, highPtFunc->GetParameter(j));
    //myFunc->SetParLimits(5, highPtFunc->GetParameter(0)*0.99, highPtFunc->GetParameter(0)*1.05);
    //myFunc->SetParLimits(6, highPtFunc->GetParameter(1)*0.999, highPtFunc->GetParameter(1)*1.001);
    myFunc->FixParameter(5, highPtFunc->GetParameter(0));
    myFunc->FixParameter(6, highPtFunc->GetParameter(1));
    myFunc->SetNpx(1000);
    // Perform the fit and retrieve the fit result pointer
    TFitResultPtr fitResultPtr = hist->Fit(myFunc, "SRE");
    const ROOT::Fit::FitResult &fitResult = *fitResultPtr;


    // Number of points (bins) to calculate the confidence intervals for
    unsigned int nPoints = hist->GetNbinsX();

    // Arrays to hold x values and the corresponding confidence intervals
    double *x = new double[nPoints];
    double *ciCombined = new double[nPoints];
    double *ciLow = new double[nPoints];
    double *ciHigh = new double[nPoints];

    // Fill the x array with bin centers
    for (unsigned int i = 1; i <= nPoints; ++i)
    {
        x[i - 1] = hist->GetBinCenter(i);
    }

    // Calculate the confidence intervals using the FitResult method
    fitResult.GetConfidenceIntervals(nPoints, 1, 1, x, ciCombined, 0.68, false);  // For combined fit
    fitResultL.GetConfidenceIntervals(nPoints, 1, 1, x, ciLow, 0.68, false);      // For lowPtFunc
    fitResultH.GetConfidenceIntervals(nPoints, 1, 1, x, ciHigh, 0.68, false);     // For highPtFunc

    // Continue fitting until chi^2/ndf is acceptable or maximum iterations reached
    double chi2=myFunc->GetChisquare();
    double ndf=myFunc->GetNDF();
    int time = 0;
    int maxtime = 5000;
    while (chi2/ndf > 10)
    {
        if (time < maxtime)
        {
            time += 1;
            hist->Fit(myFunc, "R");
            chi2=myFunc->GetChisquare();
            ndf=myFunc->GetNDF();
        }
        else
        {
            std::cout << "Chi^2/ndf is not good enough. Exiting." << std::endl;
            break;
        }
    }

    // Prepare TGraphErrors for relative deviation plots
    TGraphErrors *gRelDevCombined = new TGraphErrors();
    TGraphErrors *gRelDevLow = new TGraphErrors();
    TGraphErrors *gRelDevHigh = new TGraphErrors();
    // characterize statistical uncertainty of the fit
    std::vector<double> bestParams(nParams), bestParamsErr(nParams),yDefault(nPoints), yUpper(nPoints), yLower(nPoints), yRatioUpper(nPoints), yRatioLower(nPoints), ymaxratioU(nPoints), ymaxratioL(nPoints);
    double bestChi2 = myFunc->GetChisquare();
    double bestChi2NDF = myFunc->GetChisquare() / myFunc->GetNDF();
    // Store the current best parameters
    for (int i = 0; i < nParams; ++i) {
        bestParams[i] = myFunc->GetParameter(i);
    }




    int pointIndexCombined = 0;
    int pointIndexLow = 0;
    int pointIndexHigh = 0;

    for (int i = 1; i <= hist->GetNbinsX(); i++)
    {
        double binpT = hist->GetBinCenter(i);
        double data = hist->GetBinContent(i);
        double error = hist->GetBinError(i);

        if (data == 0 || binpT < lowedge || binpT > highedge)
            continue;

        // Combined fit
        double fitCombined = myFunc->Eval(binpT);
        double fitErrorCombined = ciCombined[i - 1];

        if (fitCombined == 0)
            continue;

        double deviationCombined = (data - fitCombined) / fitCombined;
        double deviationErrorCombined = sqrt(pow(data * fitErrorCombined / pow(fitCombined, 2), 2) + pow(error / data, 2));

        gRelDevCombined->SetPoint(pointIndexCombined, binpT, deviationCombined);
        gRelDevCombined->SetPointError(pointIndexCombined, 0, deviationErrorCombined);
        pointIndexCombined++;

        // Low pT fit
        double fitLow = lowPtFunc->Eval(binpT);
        double fitErrorLow = ciLow[i - 1];

        if (fitLow == 0)
            continue;

        double deviationLow = (data - fitLow) / fitLow;
        double deviationErrorLow = sqrt(pow(data * fitErrorLow / pow(fitLow, 2), 2) + pow(error / data, 2));

        gRelDevLow->SetPoint(pointIndexLow, binpT, deviationLow);
        gRelDevLow->SetPointError(pointIndexLow, 0, deviationErrorLow);
        pointIndexLow++;

        // High pT fit
        double fitHigh = highPtFunc->Eval(binpT);
        double fitErrorHigh = ciHigh[i - 1];

        if (fitHigh == 0)
            continue;

        double deviationHigh = (data - fitHigh) / fitHigh;
        double deviationErrorHigh = sqrt(pow(data * fitErrorHigh / pow(fitHigh, 2), 2) + pow(error / data, 2));

        gRelDevHigh->SetPoint(pointIndexHigh, binpT, deviationHigh);
        gRelDevHigh->SetPointError(pointIndexHigh, 0, deviationErrorHigh);
        pointIndexHigh++;
    }

    // Set the number of valid points explicitly
    gRelDevCombined->Set(pointIndexCombined);
    gRelDevLow->Set(pointIndexLow);
    gRelDevHigh->Set(pointIndexHigh);

    // Clean up dynamically allocated memory
    delete[] x;
    delete[] ciCombined;
    delete[] ciLow;
    delete[] ciHigh;
    //-------------------------------------------------------------------------------------------------------------- canvas
    // Woods-Saxon function for transition visualization
    float transitionwidth = 0.3;
    TF1 *transitionFunc = new TF1("transitionFunc", "[2] / (1 + exp((x - [0]) / [1])) - [2] / 2", 0, 20);
    transitionFunc->SetParameters(myFunc->GetParameter(0), myFunc->GetParameter(1), transitionwidth); // 0.15 sets arbitrary height
    transitionFunc->SetNpx(1000);

    // Plot the spectrum with fits
    TCanvas *c1 = new TCanvas("c1", "Canvas1", 800, 600);
    hist->SetTitle("Pion Spectrum; #pi_{0} p_{T} (GeV/c); 1/p_{T} #times dN/dp_{T}");
    gStyle->SetOptFit(112);
    gStyle->SetOptStat(0);

    // Set x-axis range
    hist->GetXaxis()->SetRangeUser(1, 20);
    hist->Draw();

    // Draw the fit functions
    myFunc->SetLineColor(kRed);
    myFunc->Draw("same");
    if (showauxilliaryfits)
    {
        lowPtFunc->SetLineColor(kGreen);
        lowPtFunc->Draw("same");
        highPtFunc->SetLineColor(kBlue);
        highPtFunc->Draw("same");
    }

    // Create a legend and add entries
    TLegend *legend1 = new TLegend(0.55, 0.75, 0.9, 0.9);
    legend1->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
    legend1->AddEntry(hist, "Pythia: p+p mb #sqrt{s_{NN}} = 200 GeV", "L");
    legend1->AddEntry(myFunc, "Fit (Hagedorn + Power Law)", "L");
    if (showauxilliaryfits)
    {
        legend1->AddEntry(lowPtFunc, "Low p_{T} Fit (Hagedorn)", "L");
        legend1->AddEntry(highPtFunc, "High p_{T} Fit (Power Law)", "L");
    }
    legend1->SetTextSize(0.03);
    legend1->Draw("same");

    c1->SetLogy();
    c1->SetLogx();
    c1->Print("pioncode/canvas_pdf/spectrum_fit.pdf");

    // Plot the relative deviations
    TCanvas *c2 = new TCanvas("c2", "Canvas2", 800, 600);
    float reldevrange = 1.0;

    // Create a TMultiGraph to hold all deviation graphs
    TMultiGraph *mgDev = new TMultiGraph();

    // Customize markers and colors
    //gRelDevCombined->SetMarkerStyle(20); // Circle
    //gRelDevCombined->SetMarkerColor(kRed);
    gRelDevCombined->SetMarkerStyle(24); // Medium open circle
    gRelDevCombined->SetMarkerSize(1.5); // Increase marker size
    mgDev->Add(gRelDevCombined, "P");

    gRelDevLow->SetMarkerStyle(21); // Square
    gRelDevLow->SetMarkerColor(kRed);
    mgDev->Add(gRelDevLow, "P");

    gRelDevHigh->SetMarkerStyle(22); // Triangle
    gRelDevHigh->SetMarkerColor(kGreen+1);
    mgDev->Add(gRelDevHigh, "P");

    mgDev->SetTitle("Relative Deviation of Fits from Data; #pi_{0} p_{T} (GeV/c); (Data - Fit) / Fit");

    mgDev->Draw("AP");
    mgDev->SetMinimum(-reldevrange);
    mgDev->SetMaximum(reldevrange);
    mgDev->GetXaxis()->SetLimits(0, 20);
    c2->Update();

    // Draw the transition function
    transitionFunc->SetLineColor(kGray);
    transitionFunc->Draw("same");

    // Add a TLine at deviation = 0
    double xmin = mgDev->GetXaxis()->GetXmin();
    double xmax = mgDev->GetXaxis()->GetXmax();
    TLine *line = new TLine(xmin, 0, xmax, 0);
    line->SetLineColor(kBlack);
    line->SetLineWidth(3);
    line->Draw("same");

    // Add a filled region to the right of pT = 10
    double ymin = -reldevrange;
    double ymax = reldevrange;
    TBox *box = new TBox(10, ymin, xmax, ymax);
    box->SetFillColorAlpha(kRed, 0.3);
    box->SetFillStyle(3004);   // Optional: Set the fill style (diagonal hatching, etc.)
    box->Draw("same");

    // Create and draw the legend
    TLegend *legend = new TLegend(0.2, 0.75, 0.4, 0.9);
    legend->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
    legend->AddEntry(gRelDevCombined, "Combined Fit ", "P");
    legend->AddEntry(gRelDevLow, "Hagedorn", "P");
    legend->AddEntry(gRelDevHigh, "Power Law", "P");
    legend->AddEntry(box, "Low Stats Region", "F");
    legend->SetFillStyle(0);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.03);
    legend->Draw("same");

    c2->Print("pioncode/canvas_pdf/relative_deviation.pdf");
}

