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

    for (int i = 1; i <= hist->GetNbinsX(); i++) {
        double binpT = hist->GetBinCenter(i);       // Get the bin center (pT)
        double binContent = hist->GetBinContent(i); // Get the bin content
        double newContent = binContent / binpT;     // Scale the bin content by the bin center (pT)
        hist->SetBinContent(i, newContent);         // Set the new bin content
        hist->SetBinError(i, hist->GetBinError(i) / binpT); // Scale the bin error by the bin center (pT)
    }

    int lowedge = 1;
    int highedge = 11;
    float transition = 3.9;
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
    hist->Fit(lowPtFunc, "R");
    hist->Fit(highPtFunc, "R");


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
    double *x = new double[nPoints];
    double *ci = new double[nPoints];

    // Fill the x array with bin centers
    for (unsigned int i = 1; i <= nPoints; ++i)
    {
        x[i - 1] = hist->GetBinCenter(i);
    }

    // Calculate the confidence intervals using the FitResult method
    fitResult.GetConfidenceIntervals(nPoints, 1, 1, x, ci, 0.68, false);  // 68% confidence interval


    // if chi^2/ndf is not good enough, continue to fit. break at time intervals so this is not indefinite
    double chi2=myFunc->GetChisquare();
    double ndf=myFunc->GetNDF();
    int time = 0;
    int maxtime = 5000;//5000
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
    //characterize the relative deviation of the fit from the data
    // find fit-data/data for each bin and plot it on a graph as a function of pT
    TGraphErrors *gRelDev = new TGraphErrors(hist->GetNbinsX());
    int pointIndex = 0;
    for (int i = 1; i <= hist->GetNbinsX(); i++)
    {
        double binpT = hist->GetBinCenter(i);
        double data = hist->GetBinContent(i);
        double error = hist->GetBinError(i);
        double fit = myFunc->Eval(binpT);

        double fitError = ci[i - 1]; // Use GetErrorY for confidence interval error

        //double fiterror = myFunc->GetParError(i);
        if(data == 0 || fit == 0 || binpT < lowedge || binpT > highedge) 
            continue; // Skip invalid points

        double deviation = (data - fit) / fit;
        gRelDev->SetPoint(pointIndex, binpT, deviation);
        
        double deviationerror = sqrt(pow(data*fitError/pow(fit,2), 2) + pow(error/data, 2));
        //need to add correlation term
        gRelDev->SetPointError(pointIndex, 0, deviationerror);

        //error debug line 
        std::cout << "binpT: " << binpT << ", fit: " << fit << ", fit error: " << fitError << ", data: " << data << ", error: " << error << ", deviation: " << deviation << ", deviation error: " << deviationerror << std::endl;
        pointIndex++; // Only increment for valid points
        //relative errors debug line
        std::cout <<  " data relative error: " << error/data << ", fit relative error: " << fitError/fit << ", deviation relative error: " << deviationerror/deviation << std::endl;
    }
    gRelDev->Set(pointIndex); // Sets the number of valid points explicitly
    // Clean up dynamically allocated memory
    delete[] x;
    delete[] ci;

    // woods saxon function for transition between low and high pt
    float transitionwidth = 0.3;
    TF1 *transitionFunc = new TF1("transitionFunc", "[2] / (1 + exp((x - [0]) / [1]))-[2]/2", 0, 20);
    transitionFunc->SetParameters(myFunc->GetParameter(0), myFunc->GetParameter(1), transitionwidth);//0.15 sets arbitrary height
    transitionFunc->SetNpx(1000);

    TCanvas *c1 = new TCanvas("c1", "Canvas1", 800, 600);
    hist->SetTitle("Pion Spectrum; #pi_{0} p_{T} (GeV/c); 1/p_{T}#times dN/dp_{T}");
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
    // Create a legend and add entries for the histogram and the fit
    TLegend *legend1 = new TLegend(0.55, 0.75, 0.9, 0.9); // Position (x1, y1, x2, y2)
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
    gRelDev->SetTitle("Relative Deviation of Fit from Data; #pi_{0} p_{T} (GeV/c); (Data - Fit) / Fit");
    //add a Tline at deviation=0, across the entire x-axis
    //set y range
    float reldevrange = 0.6;
    gRelDev->SetMinimum(-reldevrange);
    gRelDev->SetMaximum(reldevrange);
    gRelDev->GetXaxis()->SetLimits(0, 20);
    gRelDev->Draw("AP");
    //draw the transition function
    transitionFunc->SetLineColor(kGray);
    transitionFunc->Draw("same");

    // Add a TLine at deviation = 0
    double xmin = gRelDev->GetXaxis()->GetXmin();
    double xmax = gRelDev->GetXaxis()->GetXmax();
    TLine *line = new TLine(xmin, 0, xmax, 0);
    line->SetLineColor(kBlack); // Optional: Set the color of the line
    //line->SetLineStyle(2); // Optional: Set the line style (e.g., dashed)
    line->SetLineWidth(3);        // Set the line width (increase value for thicker line)
    line->Draw("same");

    // Add a filled region to the right of pT = 10
    double ymin = gRelDev->GetYaxis()->GetXmin(); // Get the y-axis min value
    double ymax = gRelDev->GetYaxis()->GetXmax(); // Get the y-axis max value
    if(reldevrange){
        ymin = -reldevrange;
        ymax = reldevrange;
    }
    TBox *box = new TBox(10, ymin, xmax, ymax);   // Create a box from pT = 10 to xmax
    box->SetFillColor(kRed);   // Set the fill color (e.g., yellow)
    box->SetFillStyle(3004);   // Optional: Set the fill style (diagonal hatching, etc.)
    box->Draw("same");

    //top left (0.1, 0.75, 0.4, 0.9)
    //top right (0.6, 0.75, 0.9, 0.9)
    TLegend *legend = new TLegend(0.2, 0.75, 0.4, 0.9); // Position (x1, y1, x2, y2)
    legend->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
    legend->AddEntry(gRelDev, "Relative Deviation", "P"); // Add graph to the legend
    legend->AddEntry(box, "Low Stats Region", "F"); // Add the filled region to the legend
    legend->SetFillStyle(0);  // Set fill style to 0 (transparent)
    legend->SetBorderSize(0); // Set border size to 0 (remove the border)
    legend->SetTextSize(0.03);  // Default is 0.04, so decrease to fit better
    legend->Draw("same");

    c2->Print("pioncode/canvas_pdf/relative_deviation.pdf");



    //gApplication->Terminate(0);
}

