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
#include "sPhenixStyle.h"
#include "sPhenixStyle.C"

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

TH1D *rebinHistogram(TH1D *h, const std::vector<double> &binEdges)
{
  const int Nbins = binEdges.size() - 1;
  double bins[Nbins + 1];
  for (int i = 0; i < Nbins + 1; ++i)
  {
    bins[i] = binEdges[i];
  }
  return (TH1D *)h->Rebin(Nbins, "hrb", bins);
}

void AnalyzeHistograms(const std::vector<std::string> &GeantFileNames, const std::vector<std::string> &FastMCFileNames, const std::vector<std::string> &GeanthistNames, const std::vector<std::string> &FastMChistNames)
{
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2");
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
  ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(10000);
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.001);
  ROOT::Math::MinimizerOptions::SetDefaultPrecision(1e-12);
  // ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(3);
  SetsPhenixStyle();
  //
  bool var_bins = true;
  int rebinFactor = 1;
  bool dynamic_left = true;
  int startBin = 1;
  int endBin = -1;
  int projectionBins = 1;
  double scale_factor = 1.0;
  double limits[8] = {0.1, 0.2, 0.11, 0.19, 0.52, 0.68, 0.50, 0.64};

  TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
  TLegend *legend1 = new TLegend(0.7, 0.7, 0.9, 0.9);

  std::vector<double> Pion_Mean, Pion_Width, Eta_Mean, Eta_Width, Mass_Ratio;
  std::vector<double> Pion_Mean_errors, Pion_Width_errors, Eta_Mean_errors, Eta_Width_errors, Mass_Ratio_errors;
  std::vector<double> pT_Bins, pT_Bins_Errors;

  TMultiGraph *gPionMeans = new TMultiGraph();
  TMultiGraph *gPionWidths = new TMultiGraph();
  TMultiGraph *gEtaMeans = new TMultiGraph();
  TMultiGraph *gEtaWidths = new TMultiGraph();
  TMultiGraph *gMassRatios = new TMultiGraph();
  int totalfiles = GeantFileNames.size(); //+ 1 + FastMCFileNames.size()
  std::vector<TGraphErrors *> pionmeanGraph(totalfiles);
  std::vector<TGraphErrors *> pionwidthGraph(totalfiles);
  std::vector<TGraphErrors *> etameanGraph(totalfiles);
  std::vector<TGraphErrors *> etawidthGraph(totalfiles);
  std::vector<TGraphErrors *> massRatioGraph(totalfiles);

  // Loop through Geant4 files
  for (size_t j = 0; j < GeantFileNames.size(); ++j)
  {
    TFile file(GeantFileNames[i].c_str(), "READ");
    if (!file.IsOpen())
    {
      std::cerr << "Error opening file: " << GeantFileNames[j] << std::endl;
      continue;
    }

    TH2 *hist2D = dynamic_cast<TH2 *>(file.Get(GeanthistNames[j].c_str()));
    if (!hist2D)
    {
      std::cerr << "Error getting histogram: " << GeanthistNames[j] << " from file: " << GeantFileNames[j] << std::endl;
      file.Close();
      continue;
    }
    int nBinsX = hist2D->GetNbinsX();
    if (endBin == -1)
      endBin = nBinsX; // Default to the last bin if not specified

    // Loop over the x-axis bins
    for (int i = startBin; i <= endBin; i += projectionBins)
    {
      // Project the histogram along the Y-axis
      int lastBin = std::min(i + projectionBins - 1, nBinsX);
      TH1D *yProjection = hist2D->ProjectionY(Form("proj_%d", i), i, lastBin);
      // Check if the projection has enough entries to perform a fit
      if (hist->GetEntries() < 1000)
      { // Adjust the threshold as needed
        delete hist;
        continue;
      }
      TH1D *histF = (TH1D *)yProjection;
      // re binning
      if (var_bins && !nuBins.empty())
      {
        std::cout << "Rebinning histogram with non-uniform edges" << std::endl;
        histF = rebinHistogram(histF, nuBins); // nuBins
      }
      else if (rebinFactor > 1)
      {
        histF->Rebin(rebinFactor);
      }

      histF->Scale(1. / 2, "width");

      // Determine the leftmost point with a value in the projection histograms
      float leftmost_limit = 0;
      if (dynamic_left)
      {
        for (int bin = 1; bin <= histF->GetNbinsX(); ++bin)
        {
          if (histF->GetBinContent(bin) > 0)
          {
            leftmost_limit = histF->GetBinLowEdge(bin);
            // limits[0] = leftmost_limit;
            break;
          }
        }
      }
      double pt_min = hist2D->GetXaxis()->GetBinLowEdge(i);
      double pt_max = hist2D->GetXaxis()->GetBinUpEdge(lastBin);
      TString ptRange = Form("pt_%.2f-%.2f_GeV", pt_min, pt_max);
      scale_histogram_errors(histF, scale_factor);
      // fitting background only
      TF1 *leftRightFit;
      leftRightFit = new TF1("leftRightFit", poly5BG, leftmost_limit, 1.2, 6);
      histF->Fit(leftRightFit, "RE");

      // Fit first Gaussian in the specified range
      TF1 *gausFit = new TF1("gausFit", "gaus", limits[2], limits[3]);
      gausFit->SetParLimits(1, 0.11, 0.19);
      gausFit->SetParLimits(2, 0.01, 0.25);
      histF->Fit(gausFit, "RE");
      // Fit second Gaussian in the specified range
      TF1 *gausFit2 = new TF1("gausFit2", "gaus", limits[6], limits[7]);
      gausFit2->SetParLimits(1, 0.50, 0.64);
      gausFit2->SetParLimits(2, 0.03, 0.25);
      histF->Fit(gausFit2, "RE");

      // combined fit setup
      TF1 *combinedFit;
      combinedFit = new TF1("combinedFit", combinedFunctionDoubleGaussPoly5, limits[0], limits[1], 12);
      for (int j = 0; j < 3; ++j)
      {
        combinedFit->SetParameter(j, gausFit->GetParameter(j));
      }
      combinedFit->SetParLimits(0, 0, gausFit->GetParameter(0) * 1.25); // gausFit->GetParameter(0) *0.5
      combinedFit->SetParLimits(1, 0.11, 0.19);
      combinedFit->SetParLimits(2, 0.01, 0.25);
      for (int j = 0; j < 6; ++j)
        combinedFit->SetParameter(j + 6, leftRightFit->GetParameter(j));
      for (int j = 0; j < 3; ++j)
        combinedFit->SetParameter(j + 3, gausFit2->GetParameter(j));
      combinedFit->SetParLimits(3, 0, gausFit2->GetParameter(0) * 1.1); // gausFit2->GetParameter(0) *0.5
      combinedFit->SetParLimits(4, 0.5, 0.64);
      combinedFit->SetParLimits(5, 0.03, 0.25);

      // Fit the combined function
      combinedFit->SetNpx(1000);
      histF->Fit(combinedFit, "RE");

      // Get the fit parameters
      double Pmean = fitFunc->GetParameter(1);
      double Psigma = fitFunc->GetParameter(2);
      double PmeanErr = fitFunc->GetParError(1);
      double PsigmaErr = fitFunc->GetParError(2);
      double PWidth = Psigma / Pmean;
      double PWidthErr = PWidth * sqrt(pow(PmeanErr / Pmean, 2) + pow(PsigmaErr / Psigma, 2));

      double Emean = fitFunc->GetParameter(4);
      double Esigma = fitFunc->GetParameter(5);
      double EmeanErr = fitFunc->GetParError(4);
      double EsigmaErr = fitFunc->GetParError(5);
      double EWidth = Esigma / Emean;
      double EWidthErr = EWidth * sqrt(pow(EmeanErr / Emean, 2) + pow(EsigmaErr / Esigma, 2));

      double MassRatio = Pmean / Emean;
      double MassRatioErr = MassRatio * sqrt(pow(PmeanErr / Pmean, 2) + pow(EmeanErr / Emean, 2));

      Pion_Mean.push_back(Pmean);
      Pion_Width.push_back(PWidth);
      Pion_Mean_errors.push_back(PmeanErr);
      Pion_Width_errors.push_back(PWidthErr);

      Eta_Mean.push_back(Emean);
      Eta_Width.push_back(EWidth);
      Eta_Mean_errors.push_back(EmeanErr);
      Eta_Width_errors.push_back(EWidthErr);

      pT_Bins.push_back(hist->GetXaxis()->GetBinCenter(binX));
      pT_Bins_Errors.push_back(0);

      pionmeanGraph[j]->SetPoint(binX, hist->GetXaxis()->GetBinCenter(binX), Pmean);
      pionmeanGraph[j]->SetPointError(binX, 0, PmeanErr);
      pionwidthGraph[j]->SetPoint(binX, hist->GetXaxis()->GetBinCenter(binX), PWidth);
      pionwidthGraph[j]->SetPointError(binX, 0, PWidthErr);
      etameanGraph[j]->SetPoint(binX, hist->GetXaxis()->GetBinCenter(binX), Emean);
      etameanGraph[j]->SetPointError(binX, 0, EmeanErr);
      etawidthGraph[j]->SetPoint(binX, hist->GetXaxis()->GetBinCenter(binX), EWidth);
      etawidthGraph[j]->SetPointError(binX, 0, EWidthErr);
      massRatioGraph[j]->SetPoint(binX, hist->GetXaxis()->GetBinCenter(binX), MassRatio);
      massRatioGraph[j]->SetPointError(binX, 0, MassRatioErr);
    }

    int MarkerStyle = j + 24;
    int MarkerColor = j + 1;
    pionmeanGraph[j]->SetMarkerStyle(MarkerStyle);
    pionmeanGraph[j]->SetMarkerColor(MarkerColor);
    // pionmeanGraph[j]->SetMarkerSize(1.5);
    // pionmeanGraph[j]->SetLineColor(MarkerColor);
    // pionmeanGraph[j]->SetLineWidth(2);
    // pionmeanGraph[j]->SetLineStyle(1);
    // pionmeanGraph[j]->SetFillColor(0);
    // pionmeanGraph[j]->SetFillStyle(0);

    pionwidthGraph[j]->SetMarkerStyle(MarkerStyle);
    pionwidthGraph[j]->SetMarkerColor(MarkerColor);

    etameanGraph[j]->SetMarkerStyle(MarkerStyle);
    etameanGraph[j]->SetMarkerColor(MarkerColor);

    etawidthGraph[j]->SetMarkerStyle(MarkerStyle);
    etawidthGraph[j]->SetMarkerColor(MarkerColor);

    massRatioGraph[j]->SetMarkerStyle(MarkerStyle);
    massRatioGraph[j]->SetMarkerColor(MarkerColor);

    gPionMeans->Add(pionmeanGraph[j], "PE");
    legend1->AddEntry(pionmeanGraph[j], HistLegend[j].c_str(), "P");

    gPionWidths->Add(pionwidthGraph[j], "PE");
    legend1->AddEntry(pionwidthGraph[j], HistLegend[j].c_str(), "P");

    gEtaMeans->Add(etameanGraph[j], "PE");
    legend1->AddEntry(etameanGraph[j], HistLegend[j].c_str(), "P");

    gEtaWidths->Add(etawidthGraph[j], "PE");
    legend1->AddEntry(etawidthGraph[j], HistLegend[j].c_str(), "P");

    gMassRatios->Add(massRatioGraph[j], "PE");
    legend1->AddEntry(massRatioGraph[j], HistLegend[j].c_str(), "P");

    file.Close();
  }

  // Repeat the same for FastMC files
  /*
  for (size_t i = 0; i < FastMCFileNames.size(); ++i)
  {
    TFile file(FastMCFileNames[i].c_str(), "READ");
    if (!file.IsOpen())
    {
      std::cerr << "Error opening file: " << FastMCFileNames[i] << std::endl;
      continue;
    }

    TH2 *hist = dynamic_cast<TH2 *>(file.Get(FastMChistNames[i].c_str()));
    if (!hist)
    {
      std::cerr << "Error getting histogram: " << FastMChistNames[i] << " from file: " << FastMCFileNames[i] << std::endl;
      file.Close();
      continue;
    }

    for (int binX = 1; binX <= hist->GetNbinsX(); ++binX)
    {
      TH1D *yProjection = hist->ProjectionY(Form("YProjection_%zu_%d", i, binX), binX, binX, "");
      yProjection->Fit("gaus", "QE");

      TF1 *fitFunc = yProjection->GetFunction("gaus");
      if (fitFunc)
      {
        double mean = fitFunc->GetParameter(1);
        double sigma = fitFunc->GetParameter(2);
        double meanErr = fitFunc->GetParError(1);
        double sigmaErr = fitFunc->GetParError(2);

        Eta_Mean.push_back(mean);
        Eta_Width.push_back(sigma / mean);
        Eta_Mean_errors.push_back(meanErr);
        Eta_Width_errors.push_back(sigmaErr / mean);

        pT_Bins.push_back(hist->GetXaxis()->GetBinCenter(binX));
        pT_Bins_Errors.push_back(0);
      }

      delete yProjection;
    }

    file.Close();
  }
*/
  // draw multigraphs
  gPionMeans->SetTitle("Pion: Smeared pT vs Inv. Mass;#it{pT}_{#gamma#gamma} (GeV); Pion Peak Position (GeV)");
  gPionMeans->Draw("APE");
  legend1->Draw();
  // gPionMeans->GetXaxis()->SetLimits(comparisonFilenameObj.plotxlims[0], comparisonFilenameObj.plotxlims[1]);
  // gPionMeans->SetMinimum(comparisonFilenameObj.plotylims[0]);
  // gPionMeans->SetMaximum(comparisonFilenameObj.plotylims[1]);
  // c1->SetMargin(0.2, 0.1, 0.1, 0.1);
  // gPad->Modified();
  // gPad->Update();
  // c1->SaveAs(Form("%s/%
  gPionWidths->SetTitle("Pion: Smeared pT vs Inv. Mass;#it{pT}_{#gamma#gamma} (GeV); Pion Peak Position (GeV)");
  gPionWidths->Draw("APE");
  legend1->Draw();

  gEtaMeans->SetTitle("Eta: Smeared pT vs Inv. Mass;#it{pT}_{#gamma#gamma} (GeV); Eta Peak Position (GeV)");
  gEtaMeans->Draw("APE");
  legend1->Draw();

  gEtaWidths->SetTitle("Eta: Smeared pT vs Inv. Mass;#it{pT}_{#gamma#gamma} (GeV); Eta Peak Position (GeV)");
  gEtaWidths->Draw("APE");
  legend1->Draw();

  gMassRatios->SetTitle("Mass Ratios: Smeared pT vs Inv. Mass;#it{pT}_{#gamma#gamma} (GeV); Pion/Eta Mass Ratio");
  gMassRatios->Draw("APE");
  legend1->Draw();

  TFile outputFile("ptdifferential_overlay.root", "RECREATE");
  gPionMean->Write("gPionMean");
  gPionWidth->Write("gPionWidth");
  gEtaMean->Write("gEtaMean");
  gEtaWidth->Write("gEtaWidth");
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
