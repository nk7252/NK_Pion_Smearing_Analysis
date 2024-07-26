#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TH2.h>
#include <TF1.h>
#include <TMultiGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <Math/MinimizerOptions.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

// global bin var
std::vector<double> nuBins = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.36, 0.40, 0.44, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58, 0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1.0}; //, 1.04, 1.08, 1.12, 1.16, 1.2

// background subtraction for geant based MC
double combinedFunctionDoubleGaussPoly5(double *x, double *par)
{
  // First Gaussian part (e.g., pion peak)
  // double gauss1 =0;
  // if (x[0] >= 0.7 && x[0] <= 0.22){
  double gauss1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));
  //}

  // Second Gaussian part (e.g., eta peak)
  // double gauss2 =0;
  // if (x[0] >= 0.3 && x[0] <= 0.8){
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

void AnalyzeHistograms(const std::vector<std::string> &GeantFileNames, const std::vector<std::string> &FastMCFileNames, const std::vector<std::string> &GeanthistNames, const std::vector<std::string> &FastMChistNames)
{
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2"); //,"Simplex", "Migrad", "Fumili", "Combined"
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
  ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(10000);
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.001);
  ROOT::Math::MinimizerOptions::SetDefaultPrecision(1e-12);
  // ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(3);
  SetsPhenixStyle();
  // create canvas and legend
  TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
  TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);
  // Initialize the vectors to store the fit results
  std::vector<double> Pion_Mean, Pion_Width, Eta_Mean, Eta_Width, Mass_Ratio;
  std::vector<double> Pion_Mean_errors, Pion_Width_errors, Eta_Mean_errors, Eta_Width_errors, Mass_Ratio_errors;
  std::vector<double> pT_Bins, pT_Bins_Errors;
  // TMultiGraph to hold graphs.
  TMultiGraph *MultiGraphs = new TMultiGraph();
  // Create Tgraphs to hold means for each histogram
  std::vector<TGraphErrors *> meanGraph(GeantFileNames.size());
  // Geant4 histogram loop
  for (size_t i = 0; i < GeantFileNames.size(); ++i)
  {
    TFile file(fileNames[i].c_str(), "READ");
    if (!file.IsOpen())
    {
      std::cerr << "Error opening file: " << GeantFileNames[i] << std::endl;
      continue;
    }

    TH2 *hist = dynamic_cast<TH2 *>(file.Get(histNames[i].c_str()));
    if (!hist)
    {
      std::cerr << "Error getting histogram: " << histNames[i] << " from file: " << GeantFileNames[i] << std::endl;
      file.Close();
      continue;
    }

    // Fit the histogram
    TF1 *fitFunc = new TF1("fitFunc", "gaus(0) + gaus(3) + pol2(6)", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    hist->Fit(fitFunc);

    // Extract parameters
    double mean1 = fitFunc->GetParameter(1);
    double sigma1 = fitFunc->GetParameter(2);
    double mean2 = fitFunc->GetParameter(4);
    double sigma2 = fitFunc->GetParameter(5);

    Pion_Mean.push_back(mean1);
    Pion_Width.push_back(sigma1 / mean1);
    Eta_Mean.push_back(mean2);
    Eta_Width.push_back(sigma2 / mean2);
    Mass_Ratio.push_back(mean1 / mean2);
    // Assuming constant errors for simplicity
    Pion_Mean_errors.push_back(fitFunc->GetParError(1));
    Pion_Width_errors.push_back(fitFunc->GetParError(2) / mean1);
    Eta_Mean_errors.push_back(fitFunc->GetParError(4));
    Eta_Width_errors.push_back(fitFunc->GetParError(5) / mean2);
    Mass_Ratio_errors.push_back(mean1 / mean2 * sqrt(pow(fitFunc->GetParError(1) / mean1, 2) + pow(fitFunc->GetParError(4) / mean2, 2)));

    pT_Bins.push_back(i);
    pT_Bins_Errors.push_back(0);

    file.Close();
  }

  
  std::cout << "position 1" << std::endl; // debug line
  MultiGraphs->SetTitle(Form("Smeared %s pT vs Inv Mass: %s weight;pT (GeV);Inv. Mass (GeV)", filenameobj.particletype.c_str(), filenameobj.weightnames[legendInt].c_str()));
  MultiGraphs->Draw("APE");

  // Create TGraphErrors
  TGraphErrors *gPionMean = new TGraphErrors(Pion_Mean.size(), pT_Bins.data(), Pion_Mean.data(), pT_Bins_Errors.data(), Pion_Mean_errors.data());
  TGraphErrors *gPionWidth = new TGraphErrors(Pion_Width.size(), pT_Bins.data(), Pion_Width.data(), pT_Bins_Errors.data(), Pion_Width_errors.data());
  TGraphErrors *gEtaMean = new TGraphErrors(Eta_Mean.size(), pT_Bins.data(), Eta_Mean.data(), pT_Bins_Errors.data(), Eta_Mean_errors.data());
  TGraphErrors *gEtaWidth = new TGraphErrors(Eta_Width.size(), pT_Bins.data(), Eta_Width.data(), pT_Bins_Errors.data(), Eta_Width_errors.data());
  TGraphErrors *gMassRatio = new TGraphErrors(Mass_Ratio.size(), pT_Bins.data(), Mass_Ratio.data(), pT_Bins_Errors.data(), Mass_Ratio_errors.data());
  std::vector<double> Pion_Mean, Pion_Width, Eta_Mean, Eta_Width, Mass_Ratio;
  std::vector<double> Pion_Mean_errors, Pion_Width_errors, Eta_Mean_errors, Eta_Width_errors, Mass_Ratio_errors;
  // Optionally save graphs to a file
  TFile outputFile("ptdifferential_overlay.root", "RECREATE");
  gPionMean->Write("gPionMean");
  gPionWidth->Write("gPionWidth");
  gEtaMean->Write("gEtaMean");
  gEtaWidth->Write("gEtaWidth");
  gMassRatio->Write("gMassRatio");
  outputFile.Close();
}

void fit_comparison()
{
  std::vector<std::string> Geant_fileNames = {"pioncode/OUTHIST_iter_DST_CALO_CLUSTER_pythia8_pp_mb_3MHz_0000000011_merged_V1.root", "pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_200_10000MeV_0000000013_00merged_V3.root"};
  std::vector<std::string> fastmc_fileNames = {"pioncode/rootfiles/PionFastMC_0.154000_sqrte_0.077000_const.root", "pioncode/rootfiles/EtaFastMC_0.154000_sqrte_0.077000_const.root"};
  std::vector<std::string> Geant_histNames = {"h_InvMass_2d", "h_InvMass_2d"};
  std::vector<std::string> fastmc_histNames = {"h100_2", "h100_2"};

  AnalyzeHistograms(Geant_fileNames, fastmc_fileNames, Geant_histNames, fastmc_histNames);

  return 0;
}
