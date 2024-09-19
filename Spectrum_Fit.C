#pragma once
// c++ includes
#include <iostream>
#include <vector>
#include <string>
// root includes
#include <TFile.h>
#include <TH2.h>
#include <TF1.h>
#include <TH1F.h>
#include <TRandom.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TBox.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TApplication.h>
#include <TROOT.h>
#include <TFitResult.h>
#include <Math/MinimizerOptions.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

// local includes
#include "sPhenixStyle.h"
#include "sPhenixStyle.C"


void Spectrum_Fit()
{
    //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
    //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("GSLMultiMin");
    //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("GSLMultiFit");
    //ROOT::Math::MinimizerOptions::SetDefaultMinimizer("GSLSimAn");


    ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
    ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
    ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(10000);
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.001);
    ROOT::Math::MinimizerOptions::SetDefaultPrecision(1e-12);

    SetsPhenixStyle();

    TFile *file = new TFile("pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_pythia8_pp_mb_0000000015_merged_V55.root", "READ");
    if (!file->IsOpen())
    {
        std::cerr << "Error opening file: " << "pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_pythia8_pp_mb_0000000015_merged_V55.root" << std::endl;
        return;
    }

    TH1F *hist = dynamic_cast<TH1F *>(file->Get("h_truth_spectrum2"));
    //h_truth_spectrum2, 200 bins 0-20
    //h_FullTruth_pt, 100 bins 0-20

    if (!hist)
    {
        std::cerr << "Error retrieving histogram: " << "h_truth_spectrum2" << " from file" << std::endl;
        file->Close();
        return;
    }

    hist->Scale(1/hist->Integral(), "width");
    //fit the histogram with a power law function
    //rebin the histogram
    hist->Rebin(2);

    //make legends
    //top left (0.1, 0.75, 0.4, 0.9)
    //top right (0.6, 0.75, 0.9, 0.9)
    TLegend *legend = new TLegend(0.2, 0.75, 0.4, 0.9); // Position (x1, y1, x2, y2)
    TLegend *legend1 = new TLegend(0.55, 0.75, 0.9, 0.9); // Position (x1, y1, x2, y2)

    for (int i = 1; i <= hist->GetNbinsX(); i++) {
        double binpT = hist->GetBinCenter(i);       // Get the bin center (pT)
        double binContent = hist->GetBinContent(i); // Get the bin content
        double newContent = binContent / binpT;     // Scale the bin content by the bin center (pT)
        hist->SetBinContent(i, newContent);         // Set the new bin content
        hist->SetBinError(i, hist->GetBinError(i) / binpT); // Scale the bin error by the bin center (pT)
    }

    int lowedge = 1;
    int highedge = 17;
    float transition = 3.5;
    //fit low pt and high pt separately, then funnel parameters into a single fit
    //low pt is a hagedorn function, high pt is a power law
    bool showauxilliaryfits = true;
    TF1 *lowPtFunc = new TF1("lowPtFunc", "[0] / pow(1 + x / [1], [2])", lowedge, transition);
    lowPtFunc->SetParameters(1525, 0.64, 7.9);
    //lowPtFunc->SetParameters(53, 1.04, 7.5);
    TF1 *highPtFunc = new TF1("highPtFunc", "[0] / (pow(x, [1]))", transition, highedge);
    highPtFunc->SetParameters(404.4, 10);
    lowPtFunc->SetNpx(1000);
    highPtFunc->SetNpx(1000);
    TFitResultPtr fitResultPtrL = hist->Fit(lowPtFunc, "SRE");  // 'S' option returns
    TFitResultPtr fitResultPtrH = hist->Fit(highPtFunc, "SRE");  // 'S' option returns
    //hist->Fit(lowPtFunc, "R");
    //hist->Fit(highPtFunc, "R");
    // Dereference pointers
    const ROOT::Fit::FitResult &fitResultL = *fitResultPtrL;  
    const ROOT::Fit::FitResult &fitResultH = *fitResultPtrH;


    //combined function with hagedorn for low pt and power law for high pt. Transition using a woods saxon function

    TF1  *myFunc = new TF1("myFunc", "((1 / (1 + exp((x - [0]) / [1]))) * [2] / pow(1 + x / [3], [4]) + (1 - (1 / (1 + exp((x - [0]) / [1])))) * [5] / (pow(x, [6])))", lowedge, highedge);
    
    //myFunc->SetParameters(4.5, 0.114, 229.6, 1.466, 10.654, 14.43, 8.1028);
    //myFunc->SetParameters(265, 2.1, 207.2, 1.968, 12.91, 442, 2);
    //myFunc->SetParameters(265, 2.1, 53.14, 1.89, 12.11, 442, 2);
    //myFunc->SetParameters(229, 14.28, 264, 1.968, 13.08, 648.9, 4.919);
    // use parameters from the low and high pt fits as starting values
    myFunc->SetParameter(0, transition);
    myFunc->SetParameter(1, 0.114);
    myFunc->SetParLimits(1, 0.114*0.6,  0.114*1.4);
    for (int j=0; j<3; j++) myFunc->SetParameter(j+2, lowPtFunc->GetParameter(j));
    for (int j=0; j<2; j++) myFunc->SetParameter(j+5, highPtFunc->GetParameter(j));
    //set parameter limits for high pt power law
    //myFunc->SetParLimits(5, 0.5*highPtFunc->GetParameter(0), 2.1*highPtFunc->GetParameter(0));
    //myFunc->SetParLimits(6, 0.5*highPtFunc->GetParameter(1), 2.1*highPtFunc->GetParameter(1));
    //myFunc->SetNpx(1000);
    //hist->Fit(myFunc, "RE");

    // Perform the fit and retrieve the fit result pointer
    TFitResultPtr fitResultPtr = hist->Fit(myFunc, "SRE");  // 'S' option returns TFitResultPtr
    // Get the FitResult from the TFitResultPtr
    const ROOT::Fit::FitResult &fitResult = *fitResultPtr;  // Dereference TFitResultPtr

    // Number of points (bins) to calculate the confidence intervals for
    unsigned int nPoints = hist->GetNbinsX();

    // Arrays to hold x values and the corresponding confidence intervals
    std::vector<double> x(nPoints);
    std::vector<std::vector<double>> ci(nPoints, std::vector<double>(3));
    std::vector<std::string> ciNames = {"Low p_{T} Fit (Hagedorn)", "High p_{T} Fit (Power Law)", "Combined Fit (Hagedorn+Power Law)"};
    // Fill the x array with bin centers
    for (unsigned int i = 1; i <= nPoints; ++i)
    {
        x[i - 1] = hist->GetBinCenter(i);
    }

    // Calculate the confidence intervals using the FitResult method
    fitResult.GetConfidenceIntervals(nPoints, 1, 1, x.data(), &ci[0][0], 0.68, false);  // 68% confidence interval
    fitResultL.GetConfidenceIntervals(nPoints, 1, 1, x.data(), &ci[0][1], 0.68, false);  // 68% confidence interval
    fitResultH.GetConfidenceIntervals(nPoints, 1, 1, x.data(), &ci[0][2], 0.68, false);  // 68% confidence interval

    // Number of parameters to fit
    const int nParams = 7;
    double params[nParams];    // To hold initial parameters
    double errors[nParams];    // To hold errors of parameters

    // Fill arrays with the initial parameters and their errors
    for (int i = 0; i < nParams; ++i) {
        params[i] = myFunc->GetParameter(i);
        errors[i] = myFunc->GetParError(i);
    }

    // Fix one parameter (e.g., transition point)
    myFunc->FixParameter(0, params[0]);

    // Variables to store the best fit results
    double bestParams[nParams];
    double bestChi2 = myFunc->GetChisquare();
    // Store the current best parameters
    for (int i = 0; i < nParams; ++i) {
        bestParams[i] = params[i];
    }

    // Start iterating over parameter space
    int iterations = 100;  // Number of iterations to explore the parameter space
    for (int i = 0; i < iterations; ++i) {
        // Randomize parameters within their error bounds
        for (int j = 1; j < nParams; ++j) {  // Start from 1 to avoid changing the fixed parameter
            params[j] = bestParams[j] + gRandom->Gaus(0, errors[j]);

            // Optionally check that the parameters remain within limits
            /*
            if (j == 1) { // For parameter 1
                params[j] = TMath::Min(TMath::Max(params[j], 0.114 * 0.6), 0.114 * 1.4);  // Keep within set limits
            }
            */
        }

        // Set new parameters to the function
        for (int j = 1; j < nParams; ++j) {
            myFunc->SetParameter(j, params[j]);
        }

        // Perform fit and check chi2
        hist->Fit(myFunc, "RQN");  // Fit without printing and reinitialization
        double newChi2 = myFunc->GetChisquare();

        if (newChi2 < bestChi2) {
            bestChi2 = newChi2;
            for (int j = 1; j < nParams; ++j) {
                bestParams[j] = params[j];
            }
        }
    }

    // Set best parameters after exploring the parameter space
    for (int i = 1; i < nParams; ++i) {
        myFunc->SetParameter(i, bestParams[i]);
    }

    // Perform final fit with best parameters
    hist->Fit(myFunc, "RE");

    // Output the best parameters and chi-square
    std::cout << "Best Fit Parameters: " << std::endl;
    for (int i = 0; i < nParams; ++i) {
        std::cout << "Parameter " << i << ": " << bestParams[i] << std::endl;
    }
    std::cout << "Best Chi2: " << bestChi2 << std::endl;

    // characterize statistical uncertainty of the fit
    double yDefault[nPoints], yUpper[nPoints], yLower[nPoints],yRatioUpper[nPoints][nParams], yRatioLower[nPoints][nParams], ymaxratioU[nPoints], ymaxratioL[nPoints];
    // Store the default fit values for later ratio computation
    for (int i = 0; i < nPoints; ++i) {
        x[i] = hist->GetBinCenter(i);
        yDefault[i] = myFunc->Eval(x[i]);
    }


    // Loop over each parameter and vary it ±1 sigma to create upper/lower envelopes
    for (int paramIndex = 1; paramIndex < nParams; ++paramIndex) {
        // Vary the parameter by +1 sigma and -1 sigma
        for (int direction = -1; direction <= 1; direction += 2) {
            myFunc->SetParameter(paramIndex, params[paramIndex] + direction * errors[paramIndex]);
            hist->Fit(myFunc, "RQ");  // Fit without printing

            for (int i = 0; i < nPoints; ++i) {
                if(i==0) continue; // Skip the first point
                double fitValue = myFunc->Eval(x[i]);
                if (direction == 1) {
                    // Upper envelope
                    yUpper[i] = fitValue;
                } else {
                    // Lower envelope
                    yLower[i] = fitValue;
                }
                if (yDefault[i] == 0) {
                    yRatioUpper[i][paramIndex] = 0;
                    yRatioLower[i][paramIndex] = 0;
                } else {
                    yRatioUpper[i][paramIndex] = yUpper[i] / yDefault[i];
                    yRatioLower[i][paramIndex] = yLower[i] / yDefault[i];
                }
                // list for debugging
                std::cout << "Parameter " << paramIndex << ", Point " << i << ", Direction " << direction << ": " << yRatioUpper[i][paramIndex] << ", " << yRatioLower[i][paramIndex] << std::endl;
            }
        }
        // Reset the parameter to its original value
        myFunc->SetParameter(paramIndex, params[paramIndex]);
    }
    //loop over the parameters and calculate the maximum deviation from the default fit for each point

    for (int i = 0; i < nPoints; ++i) {
        for (int j = 0; j < nParams; ++j) {
            double upper = yRatioUpper[i][j];
            double lower = yRatioLower[i][j];

            if(upper > ymaxratioU[i]) ymaxratioU[i] = upper;
            if(lower < ymaxratioL[i]) ymaxratioL[i] = lower;
            //check if the current deviation is larger than the previous maximum
            //if so, update the maximum
            //double ymaxlocal = TMath::Max(abs(upper - 1), abs(lower - 1));
            //if (ymaxlocal > ymaxratio[i]) ymaxratio[i] = ymaxlocal;
        }
    }


    //characterize the relative deviation of the fit from the data
    // find fit-data/data for each bin and plot it on a graph as a function of pT
    //make a vector of TGraphErrors to store the relative deviation of each fit from the data
    std::vector<TGraphErrors *> gRelDev (3);
    //TGraphErrors *gRelDev = new TGraphErrors(hist->GetNbinsX());
    int pointIndex = 0;
    double fitvalue[3][hist->GetNbinsX()];// one for each fit
    for (int j=0; j<3; j++)
    {   
        pointIndex = 0;
        gRelDev[j] = new TGraphErrors(hist->GetNbinsX());
        for (int i = 1; i <= hist->GetNbinsX(); i++)
        {
            double binpT = hist->GetBinCenter(i);
            double data = hist->GetBinContent(i);
            double error = hist->GetBinError(i);
            double fit = myFunc->Eval(binpT);
            double sigma_param =0; //placeholder for parameter uncertainty
            double fitError = ci[i-1][j]; // Use GetErrorY for confidence interval error
            double fit_sigma_total = sqrt(pow(fitError, 2) + pow(sigma_param, 2));

            if(data == 0 || fit == 0 || binpT < lowedge || binpT > highedge) 
                continue; // Skip invalid points

            double deviation = (data - fit) / fit;
            gRelDev[j]->SetPoint(pointIndex, binpT, deviation);
            
            double sigma_deviation = sqrt(pow(data*fit_sigma_total/pow(fit,2), 2) + pow(error/data, 2));

            //need to add correlation term
            gRelDev[j]->SetPointError(pointIndex, 0, sigma_deviation);

            //total fit sigma debug line, only fitError sigma_param and sigma_total
            std::cout << "Fit error: " << fitError << ", parameter error: " << sigma_param << ", total fit error: " << fit_sigma_total << std::endl;
            //error debug line 
            std::cout << "binpT: " << binpT << ", fit: " << fit << ", fit error: " << fit_sigma_total << ", data: " << data << ", error: " << error << ", deviation: " << deviation << ", deviation error: " << sigma_deviation << std::endl;
            pointIndex++; // Only increment for valid points
            //relative errors debug line
            std::cout <<  " data relative error: " << error/data << ", fit relative error: " << fit_sigma_total/fit << ", deviation relative error: " << sigma_deviation/deviation << std::endl;
            //add a new line to break up the output
            std::cout << "------------------------------------------------" << std::endl;
        }
        gRelDev[j]->Set(pointIndex); // Sets the number of valid points explicitly
    }

    // woods saxon function for transition between low and high pt
    float transitionwidth = 0.3;
    TF1 *transitionFunc = new TF1("transitionFunc", "[2] / (1 + exp((x - [0]) / [1]))-[2]/2", 0, 20);
    transitionFunc->SetParameters(myFunc->GetParameter(0), myFunc->GetParameter(1), transitionwidth);//0.15 sets arbitrary height
    transitionFunc->SetNpx(1000);

    TCanvas *c1 = new TCanvas("c1", "Canvas1", 800, 600);
    hist->SetTitle("Pion Spectrum; #pi_{0} p_{T} (GeV/c); #frac{dN}{p_{T}dp_{T}}");
    //add statbox with fit information
    //hist->SetStats(0);
    gStyle->SetOptFit(112);
    gStyle->SetOptStat(0); 
    //set x range
    hist->GetXaxis()->SetRangeUser(1, 20);
    hist->Draw();
    // Draw the fit function on the same canvas
    myFunc->SetLineColor(kRed); // Set the color of the fit function
    myFunc->Draw("same");       // Draw the fit function on the canvas
    lowPtFunc->SetLineColor(kGreen);
    lowPtFunc->Draw("same");
    highPtFunc->SetLineColor(kBlue);
    highPtFunc->Draw("same");

    legend1->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
    legend1->AddEntry(hist, "Pythia: p+p mb #sqrt{s_{NN}} = 200 GeV", "L");        // Add histogram to the legend
    legend1->AddEntry(myFunc, "Fit (Hagedorn + Power Law)", "L"); // Add the fit to the legend
    if(showauxilliaryfits)
    {
        legend1->AddEntry(lowPtFunc, "Low p_{T} Fit (Hagedorn)", "L");
        legend1->AddEntry(highPtFunc, "High p_{T} Fit (Power Law)", "L");
    }
    // Set text size for the legend (smaller value means smaller text)
    legend1->SetTextSize(0.03);  // Default is 0.04, so decrease to fit better
    legend1->Draw("same");       // Draw the legend on the canvas
    c1->SetLogy();
    c1->SetLogx();
    c1->Print("pioncode/canvas_pdf/spectrum_fit.pdf");

    TCanvas *c2 = new TCanvas("c2", "Canvas2", 800, 600);
    float reldevrange = 0.6;
    TMultiGraph *mg0 = new TMultiGraph();
    for (int j=0; j<3; j++){
        mg0->Add(gRelDev[j]);
        legend->AddEntry(gRelDev[j], ciNames[j].c_str(), "P"); // Add graph to the legend
    }
    mg0->SetTitle("Relative Deviation of Fit from Data; #pi_{0} p_{T} (GeV/c); (Data - Fit) / Fit");
    mg0->SetMinimum(-reldevrange);
    mg0->SetMaximum(reldevrange);
    mg0->GetXaxis()->SetLimits(0, 20);
    mg0->Draw("APE");
    //draw the transition function
    transitionFunc->SetLineColor(kGray);
    transitionFunc->Draw("same");

    // Add a TLine at deviation = 0
    double xmin = gRelDev[0]->GetXaxis()->GetXmin();
    double xmax = gRelDev[0]->GetXaxis()->GetXmax();
    TLine *line = new TLine(xmin, 0, xmax, 0);
    line->SetLineColor(kBlack); // Optional: Set the color of the line
    //line->SetLineStyle(2); // Optional: Set the line style (e.g., dashed)
    line->SetLineWidth(3);        // Set the line width (increase value for thicker line)
    line->Draw("same");

    // Add a filled region to the right of pT = 10
    double ymin = gRelDev[0]->GetYaxis()->GetXmin(); // Get the y-axis min value
    double ymax = gRelDev[0]->GetYaxis()->GetXmax(); // Get the y-axis max value
    if(reldevrange){
        ymin = -reldevrange;
        ymax = reldevrange;
    }
    TBox *box = new TBox(10, ymin, xmax, ymax);   // Create a box from pT = 10 to xmax
    box->SetFillColor(kRed);   // Set the fill color (e.g., yellow)
    box->SetFillStyle(3004);   // Optional: Set the fill style (diagonal hatching, etc.)
    box->Draw("same");
    legend->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
    legend->AddEntry(box, "Low Stats Region", "F"); // Add the filled region to the legend
    legend->SetFillStyle(0);  // Set fill style to 0 (transparent)
    legend->SetBorderSize(0); // Set border size to 0 (remove the border)
    legend->SetTextSize(0.03);  // Default is 0.04, so decrease to fit better
    legend->Draw("same");

    c2->Print("pioncode/canvas_pdf/relative_deviation.pdf");

    // Draw the ratio graph
    TCanvas *c3 = new TCanvas("c3", "Fit Variation with Shaded Envelope", 800, 600);
     // Draw yRatioUpper and yRatioLower on the same canvas
    TMultiGraph *mg = new TMultiGraph();
    // add tgraphs for each parameter
    std::vector<TGraph*> graphUpper(nParams);
    std::vector<TGraph*> graphLower(nParams);
    for (int j = 0; j < nParams; ++j) {
        graphUpper[j] = new TGraph(nPoints, x.data(), yRatioUpper[j]);
        graphLower[j] = new TGraph(nPoints, x.data(), yRatioLower[j]);
        mg->Add(graphUpper[j]);
        mg->Add(graphLower[j]);
    }
    mg->SetMinimum(0.5);
    mg->SetMaximum(1.5);
    //ymaxratio
    mg->Draw("APE");
    // Draw the envelope using TBox to shade the area between yLower and yUpper
    /*
    for (int i = 0; i < nPoints - 1; ++i) {
        TBox *box = new TBox(x[i], ymaxratioL[i], x[i + 1], ymaxratioU[i]);//
        box->SetFillColorAlpha(kGray, 0.5);  // Semi-transparent gray
        box->Draw("same");
    }
    //*/

    c1->Print("FitWithShadedEnvelope.png");



    //delete[] x;
    //delete[] ci;
    //gApplication->Terminate(0);
}

