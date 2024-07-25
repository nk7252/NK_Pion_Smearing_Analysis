#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraphErrors.h>

//background subtraction for geant based MC
double combinedFunctionDoubleGaussPoly5(double *x, double *par)
{
  // First Gaussian part (e.g., pion peak)
  //double gauss1 =0;
  //if (x[0] >= 0.7 && x[0] <= 0.22){
    double gauss1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));
  //}


  // Second Gaussian part (e.g., eta peak)
  //double gauss2 =0;
  //if (x[0] >= 0.3 && x[0] <= 0.8){
    double gauss2 = par[3] * exp(-0.5 * pow((x[0] - par[4]) / par[5], 2));
  //}
  // Polynomial part (5th degree)
  double poly = par[6] + par[7] * x[0] + par[8] * x[0] * x[0] + par[9] * x[0] * x[0] * x[0] + par[10] * x[0] * x[0] * x[0] * x[0] + par[11] * x[0] * x[0] * x[0] * x[0] * x[0];

  return gauss1 + gauss2 + poly;
}

double poly5BG(double *x, double *par)
{
  // 5th degree polynomial background
  // Check if x is in the range of any Gaussian fit
  if ((x[0] >= 0.1 && x[0] <= 0.2) || (x[0] >= 0.52 && x[0] <= 0.68))
  {
    TF1::RejectPoint();
    return 0;
  }
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0] + par[4] * x[0] * x[0] * x[0] * x[0] + par[5] * x[0] * x[0] * x[0] * x[0] * x[0];
}

// scale the histogram's error bars
void scale_histogram_errors(TH1D *hist_error_scale, double scale_factor)
{
  for (int i = 1; i <= hist_error_scale->GetNbinsX(); ++i)
  {
    // Get the current error
    double current_error = hist_error_scale->GetBinError(i);

    // Scale the error by the scale factor
    hist_error_scale->SetBinError(i, current_error * scale_factor);
    // std::cout << "orig bin cont: " << hist_error_scale->GetBinContent(i) << " . bin error: " << current_error << " . New bin error: " << hist_error_scale->GetBinError(i) <<std::endl;
  }
}

void AnalyzeHistograms(const std::vector<std::string>& GeantFileNames,const std::vector<std::string>& FastMCFileNames, const std::vector<std::string>& GeanthistNames, const std::vector<std::string>& FastMChistNames) {
    std::vector<double> Pion_Means, Pion_Sigmas, Eta_Means, Eta_Sigmas, Mass_Ratios;
    std::vector<double> Pion_errorsMean, Pion_errorsSigma, Eta_errorsSecondMean, Eta_errorsSecondSigma, Mass_errorsRatio;
    std::vector<double> pT_Bins, pT_Bins_Errors;
    //geant based loop
    for (size_t i = 0; i < GeantFileNames.size(); ++i) {
        TFile file(fileNames[i].c_str(), "READ");
        if (!file.IsOpen()) {
            std::cerr << "Error opening file: " << GeantFileNames[i] << std::endl;
            continue;
        }

        TH2* hist = dynamic_cast<TH2*>(file.Get(histNames[i].c_str()));
        if (!hist) {
            std::cerr << "Error getting histogram: " << histNames[i] << " from file: " << GeantFileNames[i] << std::endl;
            file.Close();
            continue;
        }

        // Fit the histogram
        TF1* fitFunc = new TF1("fitFunc", "gaus(0) + gaus(3) + pol2(6)", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
        hist->Fit(fitFunc);

        // Extract parameters
        double mean1 = fitFunc->GetParameter(1);
        double sigma1 = fitFunc->GetParameter(2);
        double mean2 = fitFunc->GetParameter(4);
        double sigma2 = fitFunc->GetParameter(5);

        means.push_back(mean1);
        sigmas.push_back(sigma1 / mean1);
        secondMeans.push_back(mean2);
        secondSigmas.push_back(sigma2 / mean2);
        ratios.push_back(mean1 / mean2);

        // Assuming constant errors for simplicity
        errorsMean.push_back(fitFunc->GetParError(1));
        errorsSigma.push_back(fitFunc->GetParError(2) / mean1);
        errorsSecondMean.push_back(fitFunc->GetParError(4));
        errorsSecondSigma.push_back(fitFunc->GetParError(5) / mean2);
        errorsRatio.push_back(mean1 / mean2 * sqrt(pow(fitFunc->GetParError(1)/mean1, 2) + pow(fitFunc->GetParError(4)/mean2, 2)));

        pT_Bins.push_back(i);
        pT_Bins_Errors.push_back(0);

        file.Close();
    }

    //FastMC based loop

    // Create TGraphErrors
    TGraphErrors* graphMean = new TGraphErrors(means.size(), pT_Bins.data(), means.data(), pT_Bins_Errors.data(), errorsMean.data());
    TGraphErrors* graphSigma = new TGraphErrors(sigmas.size(), pT_Bins.data(), sigmas.data(), pT_Bins_Errors.data(), errorsSigma.data());
    TGraphErrors* graphSecondMean = new TGraphErrors(secondMeans.size(), pT_Bins.data(), secondMeans.data(), pT_Bins_Errors.data(), errorsSecondMean.data());
    TGraphErrors* graphSecondSigma = new TGraphErrors(secondSigmas.size(), pT_Bins.data(), secondSigmas.data(), pT_Bins_Errors.data(), errorsSecondSigma.data());
    TGraphErrors* graphRatio = new TGraphErrors(ratios.size(), pT_Bins.data(), ratios.data(), pT_Bins_Errors.data(), errorsRatio.data());

    // Optionally save graphs to a file
    TFile outputFile("fitResults.root", "RECREATE");
    graphMean->Write("MeanGraph");
    graphSigma->Write("SigmaGraph");
    graphSecondMean->Write("SecondMeanGraph");
    graphSecondSigma->Write("SecondSigmaGraph");
    graphRatio->Write("RatioGraph");
    outputFile.Close();
}

void fit_comparison() {
    std::vector<std::string> Geant_fileNames = {"pioncode/OUTHIST_iter_DST_CALO_CLUSTER_pythia8_pp_mb_3MHz_0000000011_merged_V1.root", "pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_200_10000MeV_0000000013_00merged_V3.root"};
    std::vector<std::string> fastmc_fileNames = {"pioncode/rootfiles/PionFastMC_0.154000_sqrte_0.077000_const.root", "pioncode/rootfiles/EtaFastMC_0.154000_sqrte_0.077000_const.root"};
    std::vector<std::string> Geant_histNames = {"h_InvMass_2d", "h_InvMass_2d"};
    std::vector<std::string> fastmc_histNames = {"h100_2", "h100_2"};

    AnalyzeHistograms(Geant_fileNames,fastmc_fileNames, Geant_histNames,fastmc_histNames);

    return 0;
}
