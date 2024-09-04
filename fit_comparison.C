#pragma once
//c++ includes
#include <iostream>
#include <vector>
#include <string>
//root includes
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
//local includes
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

void AnalyzeHistograms(const std::vector<std::string> &unweightedFileNames, const std::vector<std::string> &unweightedhistNames, const std::vector<std::string> &unweighted_legendNames, const std::vector<std::string> &SPMC_FileNames, const std::vector<std::string> &SPMC_histNames, const std::vector<std::string> &SPMC_legendNames,std::vector<int> SPMC_FileTypes, const std::vector<std::string> &FastMC_FileNames, const std::vector<std::string> &FastMC_histNames, const std::vector<std::string> &FastMC_legendNames,std::vector<int> FastMC_FileTypes, const std::vector<std::string> &Run2024_FileNames, const std::vector<std::string> &Run2024_legendNames)
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
  bool var_bins = false;
  int rebinFactor = 1;
  bool dynamic_left = true;
  int startBin = 9;
  int endBin_global = -1;
  int projectionBins = 4;
  double scale_factor = 1.0;
  double limits[10] = {0.05, 1.0, 0.09, 0.25, 0.05,0.35, 0.52, 0.68, 0.35, 1.0};
  /*
  std::vector<float> limits = {
      polyL, polyR,              // 0,1 Polynomial fit range: left and right limits
      gauss1L, gauss1R,          // 2-3 First Gaussian fit range: left and right limits
      polygauss1L, polygauss1R,  // 4,5 Exclusion zone for left and right polynomials: first gaussian
      gauss2L, gauss2R,          // 6,7 Second Gaussian fit range (if fitting eta peak): left and right limits
      polygauss2L, polygauss2R   // 8,9 Exclusion zone for left and right polynomials: second gaussian
  };
  */
  // Create a PDF to save the canvases
  TCanvas *dummyCanvas = new TCanvas(); // to create pdf
  dummyCanvas->Print("pioncode/canvas_pdf/ptdifferentialcomparison.pdf[");
  dummyCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf[");
  //top right (0.66, 0.7, 0.90, 0.9)
  //top left (0.2, 0.7, 0.44, 0.9)
  //bottom left (0.2, 0.2, 0.44, 0.4)
  //bottom right (0.66, 0.2, 0.90, 0.4)
  TLegend *legend1 = new TLegend(0.2, 0.7, 0.44, 0.9);//pion mean
  TLegend *legend2 = new TLegend(0.2, 0.2, 0.44, 0.4);//pion width
  TLegend *legend3 = new TLegend(0.66, 0.7, 0.90, 0.9);//eta mean
  TLegend *legend4 = new TLegend(0.2, 0.7, 0.44, 0.9);//eta width
  TLegend *legend5 = new TLegend(0.2, 0.7, 0.44, 0.9);//mass ratio
  TLegend *legend6 = new TLegend(0.66, 0.7, 0.90, 0.9);//pion resolution
  TLegend *legend7 = new TLegend(0.66, 0.7, 0.90, 0.9);//eta resolution
  legend1->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
  legend2->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
  legend3->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
  legend4->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
  legend5->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
  legend6->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
  legend7->AddEntry("", "#it{#bf{sPHENIX}} Internal", "");
  

  std::vector<double> Pion_Mean, Pion_Width, Eta_Mean, Eta_Width, Mass_Ratio;
  std::vector<double> Pion_Mean_errors, Pion_Width_errors, Eta_Mean_errors, Eta_Width_errors, Mass_Ratio_errors;
  std::vector<double> pT_Bins, pT_Bins_Errors;
  int endBin = endBin_global;

  TMultiGraph *gPionMeans = new TMultiGraph();
  TMultiGraph *gPionWidths = new TMultiGraph();
  TMultiGraph *gEtaMeans = new TMultiGraph();
  TMultiGraph *gEtaWidths = new TMultiGraph();
  TMultiGraph *gMassRatios = new TMultiGraph();
  TMultiGraph *gPResolutions = new TMultiGraph();
  TMultiGraph *gEResolutions = new TMultiGraph();
  int totalfiles = unweightedFileNames.size()+ FastMC_FileNames.size() + SPMC_FileNames.size() + Run2024_FileNames.size();
  std::vector<TGraphErrors *> pionmeanGraph(totalfiles);
  std::vector<TGraphErrors *> pionwidthGraph(totalfiles);
  std::vector<TGraphErrors *> etameanGraph(totalfiles);
  std::vector<TGraphErrors *> etawidthGraph(totalfiles);
  std::vector<TGraphErrors *> massRatioGraph(totalfiles);
  std::vector<TGraphErrors *> PresolutionGraph(totalfiles);
  std::vector<TGraphErrors *> EresolutionGraph(totalfiles);

  // Initialize the graphs
  for (size_t j = 0; j < totalfiles; ++j)
  {
    pionmeanGraph[j] = new TGraphErrors();
    pionwidthGraph[j] = new TGraphErrors();
    etameanGraph[j] = new TGraphErrors();
    etawidthGraph[j] = new TGraphErrors();
    massRatioGraph[j] = new TGraphErrors();
    PresolutionGraph[j] = new TGraphErrors();
    EresolutionGraph[j] = new TGraphErrors();
  }

  int MarkerStyle = 24;
  int MarkerColor = 1; 
  int filecounter = 0;

  std::cout << "Unweighted files: " << unweightedFileNames.size() << std::endl;
  // Loop through unweighted files(pythia,data)
  for (size_t j = 0; j < unweightedFileNames.size(); ++j)
  {
    TFile file(unweightedFileNames[j].c_str(), "READ");
    if (!file.IsOpen())
    {
      std::cerr << "Error opening file: " << unweightedFileNames[j] << std::endl;
      continue;
    }

    TH2 *hist2D = dynamic_cast<TH2 *>(file.Get(unweightedhistNames[j].c_str()));
    if (!hist2D)
    {
      std::cerr << "Error getting histogram: " << unweightedhistNames[j] << " from file: " << unweightedFileNames[j] << std::endl;
      file.Close();
      continue;
    }
    int nBinsX = hist2D->GetNbinsX();
    if (endBin_global == -1)
      endBin = nBinsX; // Default to the last bin if not specified

    // Loop over the x-axis bins
    int bincounter = 1;
    for (int i = startBin; i <= endBin; i += projectionBins)
    {
      // Project the histogram along the Y-axis
      int lastBin = std::min(i + projectionBins - 1, nBinsX);
      TH1D *yProjection = hist2D->ProjectionY(Form("proj_%d", i), i, lastBin);
      // Check if the projection has enough entries to perform a fit
      if (yProjection->GetEntries() < 1000)
      { // Adjust the threshold as needed
        delete yProjection;
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
      //float leftmost_limit = 0;
      if (dynamic_left)
      {
        for (int bin = 1; bin <= histF->GetNbinsX(); ++bin)
        {
          if (histF->GetBinContent(bin) > 0)
          {
            float leftmost_limit = histF->GetBinLowEdge(bin);
            limits[0] = leftmost_limit;
            break;
          }
        }
      }

      double pt_min = hist2D->GetXaxis()->GetBinLowEdge(i);
      double pt_max = hist2D->GetXaxis()->GetBinUpEdge(lastBin);
      TString ptRange = Form("pt_%.2f-%.2f_GeV", pt_min, pt_max);
      double pion_pt = (pt_min + pt_max) / 2.0;
      scale_histogram_errors(histF, scale_factor);
      // fitting background only
      TF1 *leftRightFit;
      leftRightFit = new TF1("leftRightFit", poly5BG, limits[0], limits[1], 6);
      histF->Fit(leftRightFit, "REQ");

      // Fit first Gaussian in the specified range
      TF1 *gausFit = new TF1("gausFit", "gaus", limits[2], limits[3]);
      gausFit->SetParLimits(1, 0.11, 0.19);
      gausFit->SetParLimits(2, 0.01, 0.25);
      histF->Fit(gausFit, "REQ");
      // Fit second Gaussian in the specified range
      TF1 *gausFit2 = new TF1("gausFit2", "gaus", limits[6], limits[7]);
      gausFit2->SetParLimits(1, 0.50, 0.64);
      gausFit2->SetParLimits(2, 0.03, 0.25);
      histF->Fit(gausFit2, "REQ");

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
      histF->Fit(combinedFit, "REQ");

      // Check if the fit returns NaN or Inf
      bool fitFailed = false;
      for (int i = 0; i < combinedFit->GetNpar(); i++) {
          double param = combinedFit->GetParameter(i);
          if (std::isnan(param) || std::isinf(param)) {
              fitFailed = true;
              break;
          }
      }
      if (fitFailed) {
          std::cout << "Fit returned NaN or Inf for slice: " << i << std::endl;
          continue;
      }

      // Get the fit parameters
      double Pmean = combinedFit->GetParameter(1);
      double Psigma = combinedFit->GetParameter(2);
      double PmeanErr = combinedFit->GetParError(1);
      double PsigmaErr = combinedFit->GetParError(2);
      double PWidth = Psigma / Pmean;
      double PWidthErr = PWidth * sqrt(pow(PmeanErr / Pmean, 2) + pow(PsigmaErr / Psigma, 2));

      double Emean = combinedFit->GetParameter(4);
      double Esigma = combinedFit->GetParameter(5);
      double EmeanErr = combinedFit->GetParError(4);
      double EsigmaErr = combinedFit->GetParError(5);
      double EWidth = Esigma / Emean;
      double EWidthErr = EWidth * sqrt(pow(EmeanErr / Emean, 2) + pow(EsigmaErr / Esigma, 2));

      double MassRatio = Pmean / Emean;
      double MassRatioErr = MassRatio * sqrt(pow(PmeanErr / Pmean, 2) + pow(EmeanErr / Emean, 2));

      //Pion_Mean.push_back(Pmean);
      //Pion_Width.push_back(PWidth);
      //Pion_Mean_errors.push_back(PmeanErr);
      //Pion_Width_errors.push_back(PWidthErr);

      //Eta_Mean.push_back(Emean);
      //Eta_Width.push_back(EWidth);
      //Eta_Mean_errors.push_back(EmeanErr);
      //Eta_Width_errors.push_back(EWidthErr);

      //pT_Bins.push_back(pion_pt);
      //pT_Bins_Errors.push_back(0);

      pionmeanGraph[filecounter]->SetPoint(bincounter, pion_pt, Pmean);
      pionmeanGraph[filecounter]->SetPointError(bincounter, 0, PmeanErr);
      pionwidthGraph[filecounter]->SetPoint(bincounter, pion_pt, PWidth);
      pionwidthGraph[filecounter]->SetPointError(bincounter, 0, PWidthErr);
      etameanGraph[filecounter]->SetPoint(bincounter, pion_pt, Emean);
      etameanGraph[filecounter]->SetPointError(bincounter, 0, EmeanErr);
      etawidthGraph[filecounter]->SetPoint(bincounter, pion_pt, EWidth);
      etawidthGraph[filecounter]->SetPointError(bincounter, 0, EWidthErr);
      massRatioGraph[filecounter]->SetPoint(bincounter, pion_pt, MassRatio);
      massRatioGraph[filecounter]->SetPointError(bincounter, 0, MassRatioErr);
      PresolutionGraph[filecounter]->SetPoint(bincounter, pion_pt, PWidth);  
      PresolutionGraph[filecounter]->SetPointError(bincounter, 0, PWidthErr); 
      EresolutionGraph[filecounter]->SetPoint(bincounter, pion_pt, EWidth);  
      EresolutionGraph[filecounter]->SetPointError(bincounter, 0, EWidthErr);
      bincounter++;
    }
    MarkerStyle+=1;
    MarkerColor+=1;
    if(MarkerColor==5 || MarkerColor==10) MarkerColor+=1;//avoid yellow
    pionmeanGraph[filecounter]->SetMarkerStyle(MarkerStyle);
    pionmeanGraph[filecounter]->SetMarkerColor(MarkerColor);
    // pionmeanGraph[filecounter]->SetMarkerSize(1.5);
    // pionmeanGraph[filecounter]->SetLineColor(MarkerColor);
    // pionmeanGraph[filecounter]->SetLineWidth(2);
    // pionmeanGraph[filecounter]->SetLineStyle(1);
    // pionmeanGraph[filecounter]->SetFillColor(0);
    // pionmeanGraph[filecounter]->SetFillStyle(0);

    pionwidthGraph[filecounter]->SetMarkerStyle(MarkerStyle);
    pionwidthGraph[filecounter]->SetMarkerColor(MarkerColor);

    etameanGraph[filecounter]->SetMarkerStyle(MarkerStyle);
    etameanGraph[filecounter]->SetMarkerColor(MarkerColor);

    etawidthGraph[filecounter]->SetMarkerStyle(MarkerStyle);
    etawidthGraph[filecounter]->SetMarkerColor(MarkerColor);

    massRatioGraph[filecounter]->SetMarkerStyle(MarkerStyle);
    massRatioGraph[filecounter]->SetMarkerColor(MarkerColor);

    PresolutionGraph[filecounter]->SetMarkerStyle(MarkerStyle);
    PresolutionGraph[filecounter]->SetMarkerColor(MarkerColor);

    EresolutionGraph[filecounter]->SetMarkerStyle(MarkerStyle);
    EresolutionGraph[filecounter]->SetMarkerColor(MarkerColor);

    gPionMeans->Add(pionmeanGraph[filecounter], "PE");
    legend1->AddEntry(pionmeanGraph[filecounter], unweighted_legendNames[j].c_str(), "P");

    gPionWidths->Add(pionwidthGraph[filecounter], "PE");
    legend2->AddEntry(pionwidthGraph[filecounter], unweighted_legendNames[j].c_str(), "P");

    gEtaMeans->Add(etameanGraph[filecounter], "PE");
    legend3->AddEntry(etameanGraph[filecounter], unweighted_legendNames[j].c_str(), "P");

    gEtaWidths->Add(etawidthGraph[filecounter], "PE");
    legend4->AddEntry(etawidthGraph[filecounter], unweighted_legendNames[j].c_str(), "P");

    gMassRatios->Add(massRatioGraph[filecounter], "PE");
    legend5->AddEntry(massRatioGraph[filecounter], unweighted_legendNames[j].c_str(), "P");

    //------------------------------------------------------------------------------------------------

    // Define a function for the pion energy resolution fit
    TF1 *PresolutionFit = new TF1("PresolutionFit", "sqrt([0]*[0]/x + [1]*[1])", 0.1, 20);
    PresolutionFit->SetParameters(0.1, 0.02);  // Initial guesses for a, b

    // Fit the resolution graph
    PresolutionGraph[filecounter]->Fit(PresolutionFit, "R");  // Fit and constrain to the range of pT

    // Create a canvas to plot the resolution graph and fit
    TCanvas *PresCanvas = new TCanvas("resCanvas", "Resolution Fit", 800, 600);
    PresolutionGraph[filecounter]->SetTitle("Energy Resolution; p_{T} (GeV/c); Pion #sigma / #mu");
    PresolutionGraph[filecounter]->Draw("APE");
    PresolutionFit->Draw("same");

    // Print the fit parameters on a new canvas
    TCanvas *PfitParamsCanvas = new TCanvas("fitParamsCanvas", "Fit Parameters", 800, 600);
    TPaveText *PparamsText = new TPaveText(0.1, 0.7, 0.9, 0.9, "NDC");
    PparamsText->AddText("Pion Fitted Resolution Parameters:");
    PparamsText->AddText(unweighted_legendNames[j].c_str());
    PparamsText->AddText(Form("Stochastic term (a): %.4f", PresolutionFit->GetParameter(0)));
    //paramsText->AddText(Form("Noise term (b): %.4f", PresolutionFit->GetParameter(2)));
    PparamsText->AddText(Form("Constant term (c): %.4f", PresolutionFit->GetParameter(1)));
    //add goodness of fit
    PparamsText->AddText(Form("Chi2/ndf: %.4f", PresolutionFit->GetChisquare() / PresolutionFit->GetNDF()));
    PparamsText->Draw();

    
    gPResolutions->Add(PresolutionGraph[filecounter], "PE");
    legend6->AddEntry(PresolutionGraph[filecounter], unweighted_legendNames[j].c_str(), "P");

    // Save the plot to the PDF
    PresCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
    PresCanvas->Close();
    PfitParamsCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
    PfitParamsCanvas->Close();

    //------------------------------------------------------------------------------------------------
    
    // Define a function for the eta energy resolution fit
    TF1 *EresolutionFit = new TF1("EresolutionFit", "sqrt([0]*[0]/x + [1]*[1])", 0.1, 20);
    EresolutionFit->SetParameters(0.1, 0.02);  // Initial guesses for a, b
    EresolutionGraph[filecounter]->Fit(EresolutionFit, "R");

    TCanvas *EresCanvas = new TCanvas("EresCanvas", "Resolution Fit", 800, 600);
    EresolutionGraph[filecounter]->SetTitle("Energy Resolution; p_{T} (GeV/c); #sigma / #mu");
    EresolutionGraph[filecounter]->Draw("APE");
    EresolutionFit->Draw("same");

    TCanvas *EfitParamsCanvas = new TCanvas("EfitParamsCanvas", "Fit Parameters", 800, 600);
    TPaveText *EparamsText = new TPaveText(0.1, 0.7, 0.9, 0.9, "NDC");
    EparamsText->AddText("Eta Fitted Resolution Parameters:");
    EparamsText->AddText(unweighted_legendNames[j].c_str());
    EparamsText->AddText(Form("Stochastic term (a): %.4f", EresolutionFit->GetParameter(0)));
    //paramsText->AddText(Form("Noise term (b): %.4f", EresolutionFit->GetParameter(2)));
    EparamsText->AddText(Form("Constant term (c): %.4f", EresolutionFit->GetParameter(1)));
    //add goodness of fit
    EparamsText->AddText(Form("Chi2/ndf: %.4f", EresolutionFit->GetChisquare() / EresolutionFit->GetNDF()));
    EparamsText->Draw();

    gEResolutions->Add(EresolutionGraph[filecounter], "PE");
    legend7->AddEntry(EresolutionGraph[filecounter], unweighted_legendNames[j].c_str(), "P");

    EresCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
    EresCanvas->Close();
    EfitParamsCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
    EfitParamsCanvas->Close();

    file.Close();
    filecounter++;
  }
  
  std::cout << "SPMC files: " << SPMC_FileNames.size() << std::endl;
  // Loop through SPMC files
  for (size_t j = 0; j < SPMC_FileNames.size(); ++j)
  {
    TFile file(SPMC_FileNames[j].c_str(), "READ");
    if (!file.IsOpen())
    {
      std::cerr << "Error opening file: " << SPMC_FileNames[j] << std::endl;
      continue;
    }

    TH2 *hist2D = dynamic_cast<TH2 *>(file.Get(SPMC_histNames[j].c_str()));
    if (!hist2D)
    {
      std::cerr << "Error getting histogram: " << SPMC_histNames[j] << " from file: " << SPMC_FileNames[j] << std::endl;
      file.Close();
      continue;
    }
    int nBinsX = hist2D->GetNbinsX();
    if (endBin_global == -1)
      endBin = nBinsX; // Default to the last bin if not specified

    // Loop over the x-axis bins
    int bincounter = 1;
    if(SPMC_FileTypes[j]==0)//pion
    {
      for (int i = startBin; i <= endBin; i += projectionBins)
      {
        // Project the histogram along the Y-axis
        int lastBin = std::min(i + projectionBins - 1, nBinsX);
        TH1D *yProjection = hist2D->ProjectionY(Form("proj_%d", i), i, lastBin);
        // Check if the projection has enough entries to perform a fit
        if (yProjection->GetEntries() < 1000)
        { // Adjust the threshold as needed
          delete yProjection;
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
        //float leftmost_limit = 0;
        if (dynamic_left)
        {
          for (int bin = 1; bin <= histF->GetNbinsX(); ++bin)
          {
            if (histF->GetBinContent(bin) > 0)
            {
              float leftmost_limit = histF->GetBinLowEdge(bin);
              limits[0] = leftmost_limit;
              break;
            }
          }
        }

        double pt_min = hist2D->GetXaxis()->GetBinLowEdge(i);
        double pt_max = hist2D->GetXaxis()->GetBinUpEdge(lastBin);
        TString ptRange = Form("pt_%.2f-%.2f_GeV", pt_min, pt_max);
        double pion_pt = (pt_min + pt_max) / 2.0;
        scale_histogram_errors(histF, scale_factor);
        
        // Fit Pion Gaussian in the specified range
        TF1 *gausFit = new TF1("gausFit", "gaus", limits[2], limits[3]);
        gausFit->SetParLimits(1, 0.11, 0.19);
        gausFit->SetParLimits(2, 0.01, 0.25);
        gausFit->SetNpx(1000);
        histF->Fit(gausFit, "RE");

        // Get the fit parameters
        double Pmean = gausFit->GetParameter(1);
        double Psigma = gausFit->GetParameter(2);
        double PmeanErr = gausFit->GetParError(1);
        double PsigmaErr = gausFit->GetParError(2);
        double PWidth = Psigma / Pmean;
        double PWidthErr = PWidth * sqrt(pow(PmeanErr / Pmean, 2) + pow(PsigmaErr / Psigma, 2));

        //Pion_Mean.push_back(Pmean);
        //Pion_Width.push_back(PWidth);
        //Pion_Mean_errors.push_back(PmeanErr);
        //Pion_Width_errors.push_back(PWidthErr);
        //pT_Bins.push_back(pion_pt);
        //pT_Bins_Errors.push_back(0);

        pionmeanGraph[filecounter]->SetPoint(bincounter, pion_pt, Pmean);
        pionmeanGraph[filecounter]->SetPointError(bincounter, 0, PmeanErr);
        pionwidthGraph[filecounter]->SetPoint(bincounter, pion_pt, PWidth);
        pionwidthGraph[filecounter]->SetPointError(bincounter, 0, PWidthErr);
        PresolutionGraph[filecounter]->SetPoint(bincounter, pion_pt, PWidth);  
        PresolutionGraph[filecounter]->SetPointError(bincounter, 0, PWidthErr); 

        bincounter++;
    }

      MarkerStyle+=1;
      MarkerColor+=1;
      if(MarkerColor==5 || MarkerColor==10) MarkerColor+=1;//avoid yellow
      pionmeanGraph[filecounter]->SetMarkerStyle(MarkerStyle);
      pionmeanGraph[filecounter]->SetMarkerColor(MarkerColor);
      // pionmeanGraph[filecounter]->SetMarkerSize(1.5);
      // pionmeanGraph[filecounter]->SetLineColor(MarkerColor);
      // pionmeanGraph[filecounter]->SetLineWidth(2);
      // pionmeanGraph[filecounter]->SetLineStyle(1);
      // pionmeanGraph[filecounter]->SetFillColor(0);
      // pionmeanGraph[filecounter]->SetFillStyle(0);

      pionwidthGraph[filecounter]->SetMarkerStyle(MarkerStyle);
      pionwidthGraph[filecounter]->SetMarkerColor(MarkerColor);
      PresolutionGraph[filecounter]->SetMarkerStyle(MarkerStyle);
      PresolutionGraph[filecounter]->SetMarkerColor(MarkerColor);


      gPionMeans->Add(pionmeanGraph[filecounter], "PE");
      legend1->AddEntry(pionmeanGraph[filecounter], SPMC_legendNames[j].c_str(), "P");

      gPionWidths->Add(pionwidthGraph[filecounter], "PE");
      legend2->AddEntry(pionwidthGraph[filecounter], SPMC_legendNames[j].c_str(), "P");

      //------------------------------------------------------------------------------------------------

      // Define a function for the pion energy resolution fit
      TF1 *PresolutionFit = new TF1("PresolutionFit", "sqrt([0]*[0]/x + [1]*[1])", 0.1, 20);
      PresolutionFit->SetParameters(0.1, 0.02);  // Initial guesses for a, b

      // Fit the resolution graph
      PresolutionGraph[filecounter]->Fit(PresolutionFit, "R");  // Fit and constrain to the range of pT

      // Create a canvas to plot the resolution graph and fit
      TCanvas *PresCanvas = new TCanvas("resCanvas", "Resolution Fit", 800, 600);
      PresolutionGraph[filecounter]->SetTitle("Energy Resolution; p_{T} (GeV/c); Pion #sigma / #mu");
      PresolutionGraph[filecounter]->Draw("APE");
      PresolutionFit->Draw("same");

      // Print the fit parameters on a new canvas
      TCanvas *PfitParamsCanvas = new TCanvas("fitParamsCanvas", "Fit Parameters", 800, 600);
      TPaveText *PparamsText = new TPaveText(0.1, 0.7, 0.9, 0.9, "NDC");
      PparamsText->AddText("Pion Fitted Resolution Parameters:");
      PparamsText->AddText(unweighted_legendNames[j].c_str());
      PparamsText->AddText(Form("Stochastic term (a): %.4f", PresolutionFit->GetParameter(0)));
      //paramsText->AddText(Form("Noise term (b): %.4f", PresolutionFit->GetParameter(2)));
      PparamsText->AddText(Form("Constant term (c): %.4f", PresolutionFit->GetParameter(1)));
      //add goodness of fit
      PparamsText->AddText(Form("Chi2/ndf: %.4f", PresolutionFit->GetChisquare() / PresolutionFit->GetNDF()));
      PparamsText->Draw();

      
      gPResolutions->Add(PresolutionGraph[filecounter], "PE");
      legend6->AddEntry(PresolutionGraph[filecounter], unweighted_legendNames[j].c_str(), "P");

      // Save the plot to the PDF
      PresCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
      PresCanvas->Close();
      PfitParamsCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
      PfitParamsCanvas->Close();

    }

    else if(FastMC_FileTypes[j]==1)//eta
    {
      for (int i = startBin; i <= endBin; i += projectionBins)
      {
        // Project the histogram along the Y-axis
        int lastBin = std::min(i + projectionBins - 1, nBinsX);
        TH1D *yProjection = hist2D->ProjectionY(Form("proj_%d", i), i, lastBin);
        // Check if the projection has enough entries to perform a fit
        if (yProjection->GetEntries() < 1000)
        { // Adjust the threshold as needed
          delete yProjection;
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
        //float leftmost_limit = 0;
        if (dynamic_left)
        {
          for (int bin = 1; bin <= histF->GetNbinsX(); ++bin)
          {
            if (histF->GetBinContent(bin) > 0)
            {
              float leftmost_limit = histF->GetBinLowEdge(bin);
              limits[0] = leftmost_limit;
              break;
            }
          }
        }

        double pt_min = hist2D->GetXaxis()->GetBinLowEdge(i);
        double pt_max = hist2D->GetXaxis()->GetBinUpEdge(lastBin);
        TString ptRange = Form("pt_%.2f-%.2f_GeV", pt_min, pt_max);
        double Eta_pt = (pt_min + pt_max) / 2.0;
        scale_histogram_errors(histF, scale_factor);
        
        // Fit Eta Gaussian in the specified range
        TF1 *gausFit = new TF1("gausFit", "gaus", limits[6], limits[7]);
        gausFit->SetParLimits(1, 0.50, 0.64);
        gausFit->SetParLimits(2, 0.03, 0.25);
        gausFit->SetNpx(1000);
        histF->Fit(gausFit, "REQ");

        // Check if the fit returns NaN or Inf
        bool fitFailed = false;
        for (int i = 0; i < gausFit->GetNpar(); i++) {
            double param = gausFit->GetParameter(i);
            if (std::isnan(param) || std::isinf(param)) {
                fitFailed = true;
                break;
            }
        }
        if (fitFailed) {
            std::cout << "Fit returned NaN or Inf for slice: " << i << std::endl;
            continue;
        }
        double Emean = gausFit->GetParameter(1);
        double Esigma = gausFit->GetParameter(2);
        double EmeanErr = gausFit->GetParError(1);
        double EsigmaErr = gausFit->GetParError(2);
        double EWidth = Esigma / Emean;
        double EWidthErr = EWidth * sqrt(pow(EmeanErr / Emean, 2) + pow(EsigmaErr / Esigma, 2));

        //double MassRatio = Pion_Mean[j] / Emean;
        //double MassRatioErr = MassRatio * sqrt(pow(Pion_Mean_errors[j] / Pion_Mean[j], 2) + pow(EmeanErr / Emean, 2));

        //Eta_Mean.push_back(Emean);
        //Eta_Width.push_back(EWidth);
        //Eta_Mean_errors.push_back(EmeanErr);
        //Eta_Width_errors.push_back(EWidthErr);
        //pT_Bins.push_back(Eta_pt);
        //pT_Bins_Errors.push_back(0);

        etameanGraph[filecounter]->SetPoint(bincounter, Eta_pt, Emean);
        etameanGraph[filecounter]->SetPointError(bincounter, 0, EmeanErr);
        etawidthGraph[filecounter]->SetPoint(bincounter, Eta_pt, EWidth);
        etawidthGraph[filecounter]->SetPointError(bincounter, 0, EWidthErr);
        //massRatioGraph[filecounter]->SetPoint(bincounter, Eta_pt, MassRatio);
        //massRatioGraph[filecounter]->SetPointError(bincounter, 0, MassRatioErr);
        EresolutionGraph[filecounter]->SetPoint(bincounter, Eta_pt, EWidth);  
        EresolutionGraph[filecounter]->SetPointError(bincounter, 0, EWidthErr);
        bincounter++;
        //std::cout << "Eta_pt: " << Eta_pt << " Emean: " << Emean << " EWidth: " << EWidth << std::endl;
        std::cout << "Bincounter: " << bincounter << std::endl;
      }

      MarkerStyle+=1;
      MarkerColor+=1;
      if(MarkerColor==5 || MarkerColor==10) MarkerColor+=1;//avoid yellow
      etameanGraph[filecounter]->SetMarkerStyle(MarkerStyle);
      etameanGraph[filecounter]->SetMarkerColor(MarkerColor);  
      // etameanGraph[j]->SetMarkerSize(1.5);
      // etameanGraph[j]->SetLineColor(MarkerColor);
      // etameanGraph[j]->SetLineWidth(2);
      // etameanGraph[j]->SetLineStyle(1);
      // etameanGraph[j]->SetFillColor(0);
      // etameanGraph[j]->SetFillStyle(0);

      etawidthGraph[filecounter]->SetMarkerStyle(MarkerStyle);
      etawidthGraph[filecounter]->SetMarkerColor(MarkerColor);

      //massRatioGraph[filecounter]->SetMarkerStyle(MarkerStyle);
      //massRatioGraph[filecounter]->SetMarkerColor(MarkerColor);

      EresolutionGraph[filecounter]->SetMarkerStyle(MarkerStyle);
      EresolutionGraph[filecounter]->SetMarkerColor(MarkerColor);

      gEtaMeans->Add(etameanGraph[filecounter], "PE");
      legend3->AddEntry(etameanGraph[filecounter], FastMC_legendNames[j].c_str(), "P");

      gEtaWidths->Add(etawidthGraph[filecounter], "PE");
      legend4->AddEntry(etawidthGraph[filecounter], FastMC_legendNames[j].c_str(), "P");

      //gMassRatios->Add(massRatioGraph[filecounter], "PE");
      //legend5->AddEntry(massRatioGraph[filecounter], FastMC_legendNames[j].c_str(), "P");

      //------------------------------------------------------------------------------------------------
      
      // Define a function for the eta energy resolution fit
      TF1 *EresolutionFit = new TF1("EresolutionFit", "sqrt([0]*[0]/x + [1]*[1])", 0.1, 20);
      EresolutionFit->SetParameters(0.1, 0.02);  // Initial guesses for a, b
      EresolutionGraph[filecounter]->Fit(EresolutionFit, "R");

      TCanvas *EresCanvas = new TCanvas("EresCanvas", "Resolution Fit", 800, 600);
      EresolutionGraph[filecounter]->SetTitle("Energy Resolution; p_{T} (GeV/c); #sigma / #mu");
      EresolutionGraph[filecounter]->Draw("APE");
      EresolutionFit->Draw("same");

      TCanvas *EfitParamsCanvas = new TCanvas("EfitParamsCanvas", "Fit Parameters", 800, 600);
      TPaveText *EparamsText = new TPaveText(0.1, 0.7, 0.9, 0.9, "NDC");
      EparamsText->AddText("Eta Fitted Resolution Parameters:");
      EparamsText->AddText(unweighted_legendNames[j].c_str());
      EparamsText->AddText(Form("Stochastic term (a): %.4f", EresolutionFit->GetParameter(0)));
      //paramsText->AddText(Form("Noise term (b): %.4f", EresolutionFit->GetParameter(2)));
      EparamsText->AddText(Form("Constant term (c): %.4f", EresolutionFit->GetParameter(1)));
      //add goodness of fit
      EparamsText->AddText(Form("Chi2/ndf: %.4f", EresolutionFit->GetChisquare() / EresolutionFit->GetNDF()));
      EparamsText->Draw();

      gEResolutions->Add(EresolutionGraph[filecounter], "PE");
      legend7->AddEntry(EresolutionGraph[filecounter], unweighted_legendNames[j].c_str(), "P");

      EresCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
      EresCanvas->Close();
      EfitParamsCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
      EfitParamsCanvas->Close();
    }

    file.Close();
    filecounter++;
  }

  std::cout << "FastMC files: " << FastMC_FileNames.size() << std::endl;
  // Repeat the same for FastMC files

  for (size_t j = 0; j < FastMC_FileNames.size(); ++j)
  {
    TFile file(FastMC_FileNames[j].c_str(), "READ");
    if (!file.IsOpen())
    {
      std::cerr << "Error opening file: " << FastMC_FileNames[j] << std::endl;
      continue;
    }

    TH2 *hist2D = dynamic_cast<TH2 *>(file.Get(FastMC_histNames[j].c_str()));
    if (!hist2D)
    {
      std::cerr << "Error getting histogram: " << FastMC_histNames[j] << " from file: " << FastMC_FileNames[j] << std::endl;
      file.Close();
      continue;
    }
    int nBinsX = hist2D->GetNbinsX();
    //std::cout << "nBinsX: " << nBinsX << std::endl;
    if (endBin_global == -1)
      endBin = nBinsX; // Default to the last bin if not specified
    std::cout << "endBin: " << endBin << std::endl;
    // Loop over the x-axis bins
    int bincounter = 1;
    if(FastMC_FileTypes[j]==0){//pion
      for (int i = startBin; i <= endBin; i += projectionBins)
      {
        // Project the histogram along the Y-axis
        int lastBin = std::min(i + projectionBins - 1, nBinsX);
        TH1D *yProjection = hist2D->ProjectionY(Form("proj_%d", i), i, lastBin);
        // Check if the projection has enough entries to perform a fit
        if (yProjection->GetEntries() < 1000)
        { // Adjust the threshold as needed
          delete yProjection;
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
        //float leftmost_limit = 0;
        if (dynamic_left)
        {
          for (int bin = 1; bin <= histF->GetNbinsX(); ++bin)
          {
            if (histF->GetBinContent(bin) > 0)
            {
              float leftmost_limit = histF->GetBinLowEdge(bin);
              limits[0] = leftmost_limit;
              break;
            }
          }
        }

        double pt_min = hist2D->GetXaxis()->GetBinLowEdge(i);
        double pt_max = hist2D->GetXaxis()->GetBinUpEdge(lastBin);
        TString ptRange = Form("pt_%.2f-%.2f_GeV", pt_min, pt_max);
        double pion_pt = (pt_min + pt_max) / 2.0;
        scale_histogram_errors(histF, scale_factor);
        
        // Fit Pion Gaussian in the specified range
        TF1 *gausFit = new TF1("gausFit", "gaus", limits[2], limits[3]);
        gausFit->SetParLimits(1, 0.11, 0.19);
        gausFit->SetParLimits(2, 0.01, 0.25);
        gausFit->SetNpx(1000);
        histF->Fit(gausFit, "REQ");

        // Check if the fit returns NaN or Inf
        bool fitFailed = false;
        for (int i = 0; i < gausFit->GetNpar(); i++) {
            double param = gausFit->GetParameter(i);
            if (std::isnan(param) || std::isinf(param)) {
                fitFailed = true;
                break;
            }
        }
        if (fitFailed) {
            std::cout << "Fit returned NaN or Inf for slice: " << i << std::endl;
            continue;
        }

        // Get the fit parameters
        double Pmean = gausFit->GetParameter(1);
        double Psigma = gausFit->GetParameter(2);
        double PmeanErr = gausFit->GetParError(1);
        double PsigmaErr = gausFit->GetParError(2);
        double PWidth = Psigma / Pmean;
        double PWidthErr = PWidth * sqrt(pow(PmeanErr / Pmean, 2) + pow(PsigmaErr / Psigma, 2));

        Pion_Mean.push_back(Pmean);
        Pion_Width.push_back(PWidth);
        Pion_Mean_errors.push_back(PmeanErr);
        Pion_Width_errors.push_back(PWidthErr);

        //pT_Bins.push_back(pion_pt);
        //pT_Bins_Errors.push_back(0);

        pionmeanGraph[filecounter]->SetPoint(bincounter, pion_pt, Pmean);
        pionmeanGraph[filecounter]->SetPointError(bincounter, 0, PmeanErr);
        pionwidthGraph[filecounter]->SetPoint(bincounter, pion_pt, PWidth);
        pionwidthGraph[filecounter]->SetPointError(bincounter, 0, PWidthErr);
        PresolutionGraph[filecounter]->SetPoint(bincounter, pion_pt, PWidth);  
        PresolutionGraph[filecounter]->SetPointError(bincounter, 0, PWidthErr); 

        bincounter++;
      }

      MarkerStyle+=1;
      MarkerColor+=1;
      if(MarkerColor==5 || MarkerColor==10) MarkerColor+=1;//avoid yellow
      pionmeanGraph[filecounter]->SetMarkerStyle(MarkerStyle);
      pionmeanGraph[filecounter]->SetMarkerColor(MarkerColor);
      // pionmeanGraph[j]->SetMarkerSize(1.5);
      // pionmeanGraph[j]->SetLineColor(MarkerColor);
      // pionmeanGraph[j]->SetLineWidth(2);
      // pionmeanGraph[j]->SetLineStyle(1);
      // pionmeanGraph[j]->SetFillColor(0);
      // pionmeanGraph[j]->SetFillStyle(0);

      pionwidthGraph[filecounter]->SetMarkerStyle(MarkerStyle);
      pionwidthGraph[filecounter]->SetMarkerColor(MarkerColor);

      PresolutionGraph[filecounter]->SetMarkerStyle(MarkerStyle);
      PresolutionGraph[filecounter]->SetMarkerColor(MarkerColor);


      gPionMeans->Add(pionmeanGraph[filecounter], "PE");
      legend1->AddEntry(pionmeanGraph[filecounter], FastMC_legendNames[j].c_str(), "P");

      gPionWidths->Add(pionwidthGraph[filecounter], "PE");
      legend2->AddEntry(pionwidthGraph[filecounter], FastMC_legendNames[j].c_str(), "P");

      //------------------------------------------------------------------------------------------------

      // Define a function for the pion energy resolution fit
      TF1 *PresolutionFit = new TF1("PresolutionFit", "sqrt([0]*[0]/x + [1]*[1])", 0.1, 20);
      PresolutionFit->SetParameters(0.1, 0.02);  // Initial guesses for a, b

      // Fit the resolution graph
      PresolutionGraph[filecounter]->Fit(PresolutionFit, "R");  // Fit and constrain to the range of pT

      // Create a canvas to plot the resolution graph and fit
      TCanvas *PresCanvas = new TCanvas("resCanvas", "Resolution Fit", 800, 600);
      PresolutionGraph[filecounter]->SetTitle("Energy Resolution; p_{T} (GeV/c); Pion #sigma / #mu");
      PresolutionGraph[filecounter]->Draw("APE");
      PresolutionFit->Draw("same");

      // Print the fit parameters on a new canvas
      TCanvas *PfitParamsCanvas = new TCanvas("fitParamsCanvas", "Fit Parameters", 800, 600);
      TPaveText *PparamsText = new TPaveText(0.1, 0.7, 0.9, 0.9, "NDC");
      PparamsText->AddText("Pion Fitted Resolution Parameters:");
      PparamsText->AddText(unweighted_legendNames[j].c_str());
      PparamsText->AddText(Form("Stochastic term (a): %.4f", PresolutionFit->GetParameter(0)));
      //paramsText->AddText(Form("Noise term (b): %.4f", PresolutionFit->GetParameter(2)));
      PparamsText->AddText(Form("Constant term (c): %.4f", PresolutionFit->GetParameter(1)));
      //add goodness of fit
      PparamsText->AddText(Form("Chi2/ndf: %.4f", PresolutionFit->GetChisquare() / PresolutionFit->GetNDF()));
      PparamsText->Draw();

      
      gPResolutions->Add(PresolutionGraph[filecounter], "PE");
      legend6->AddEntry(PresolutionGraph[filecounter], unweighted_legendNames[j].c_str(), "P");

      // Save the plot to the PDF
      PresCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
      PresCanvas->Close();
      PfitParamsCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
      PfitParamsCanvas->Close();
    }
    else if(FastMC_FileTypes[j]==1){//eta
      for (int i = startBin; i <= endBin; i += projectionBins)
      {
        // Project the histogram along the Y-axis
        int lastBin = std::min(i + projectionBins - 1, nBinsX);
        TH1D *yProjection = hist2D->ProjectionY(Form("proj_%d", i), i, lastBin);
        /*
        // Check if the projection has enough entries to perform a fit
        if (yProjection->GetEntries() < 1000)
        { // Adjust the threshold as needed
          delete yProjection;
          continue;
        }
        */
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
        //float leftmost_limit = 0;
        if (dynamic_left)
        {
          for (int bin = 1; bin <= histF->GetNbinsX(); ++bin)
          {
            if (histF->GetBinContent(bin) > 0)
            {
              float leftmost_limit = histF->GetBinLowEdge(bin);
              limits[0] = leftmost_limit;
              break;
            }
          }
        }

        double pt_min = hist2D->GetXaxis()->GetBinLowEdge(i);
        double pt_max = hist2D->GetXaxis()->GetBinUpEdge(lastBin);
        TString ptRange = Form("pt_%.2f-%.2f_GeV", pt_min, pt_max);
        double Eta_pt = (pt_min + pt_max) / 2.0;
        scale_histogram_errors(histF, scale_factor);
        
        // Fit Eta Gaussian in the specified range
        TF1 *gausFit = new TF1("gausFit", "gaus", limits[6], limits[7]);
        gausFit->SetParLimits(1, 0.50, 0.64);
        gausFit->SetParLimits(2, 0.03, 0.25);
        gausFit->SetNpx(1000);
        histF->Fit(gausFit, "REQ");

        // Check if the fit returns NaN or Inf
        bool fitFailed = false;
        for (int i = 0; i < gausFit->GetNpar(); i++) {
            double param = gausFit->GetParameter(i);
            if (std::isnan(param) || std::isinf(param)) {
                fitFailed = true;
                break;
            }
        }
        if (fitFailed) {
            std::cout << "Fit returned NaN or Inf for slice: " << i << std::endl;
            continue;
        }
        double Emean = gausFit->GetParameter(1);
        double Esigma = gausFit->GetParameter(2);
        double EmeanErr = gausFit->GetParError(1);
        double EsigmaErr = gausFit->GetParError(2);
        double EWidth = Esigma / Emean;
        double EWidthErr = EWidth * sqrt(pow(EmeanErr / Emean, 2) + pow(EsigmaErr / Esigma, 2));

        double MassRatio = Pion_Mean[j] / Emean;
        double MassRatioErr = MassRatio * sqrt(pow(Pion_Mean_errors[j] / Pion_Mean[j], 2) + pow(EmeanErr / Emean, 2));

        //Eta_Mean.push_back(Emean);
        //Eta_Width.push_back(EWidth);
        //Eta_Mean_errors.push_back(EmeanErr);
        //Eta_Width_errors.push_back(EWidthErr);
        //pT_Bins.push_back(Eta_pt);
        //pT_Bins_Errors.push_back(0);

        etameanGraph[filecounter]->SetPoint(bincounter, Eta_pt, Emean);
        etameanGraph[filecounter]->SetPointError(bincounter, 0, EmeanErr);
        etawidthGraph[filecounter]->SetPoint(bincounter, Eta_pt, EWidth);
        etawidthGraph[filecounter]->SetPointError(bincounter, 0, EWidthErr);
        massRatioGraph[filecounter]->SetPoint(bincounter, Eta_pt, MassRatio);
        massRatioGraph[filecounter]->SetPointError(bincounter, 0, MassRatioErr);
        EresolutionGraph[filecounter]->SetPoint(bincounter, Eta_pt, EWidth);  
        EresolutionGraph[filecounter]->SetPointError(bincounter, 0, EWidthErr);
        bincounter++;
        //std::cout << "Eta_pt: " << Eta_pt << " Emean: " << Emean << " EWidth: " << EWidth << std::endl;
        std::cout << "Bincounter: " << bincounter << std::endl;
      }

      MarkerStyle+=1;
      MarkerColor+=1;
      if(MarkerColor==5 || MarkerColor==10) MarkerColor+=1;//avoid yellow
      etameanGraph[filecounter]->SetMarkerStyle(MarkerStyle);
      etameanGraph[filecounter]->SetMarkerColor(MarkerColor);  
      // etameanGraph[j]->SetMarkerSize(1.5);
      // etameanGraph[j]->SetLineColor(MarkerColor);
      // etameanGraph[j]->SetLineWidth(2);
      // etameanGraph[j]->SetLineStyle(1);
      // etameanGraph[j]->SetFillColor(0);
      // etameanGraph[j]->SetFillStyle(0);

      etawidthGraph[filecounter]->SetMarkerStyle(MarkerStyle);
      etawidthGraph[filecounter]->SetMarkerColor(MarkerColor);

      massRatioGraph[filecounter]->SetMarkerStyle(MarkerStyle);
      massRatioGraph[filecounter]->SetMarkerColor(MarkerColor);

      EresolutionGraph[filecounter]->SetMarkerStyle(MarkerStyle);
      EresolutionGraph[filecounter]->SetMarkerColor(MarkerColor);

      gEtaMeans->Add(etameanGraph[filecounter], "PE");
      legend3->AddEntry(etameanGraph[filecounter], SPMC_legendNames[j].c_str(), "P");

      gEtaWidths->Add(etawidthGraph[filecounter], "PE");
      legend4->AddEntry(etawidthGraph[filecounter], SPMC_legendNames[j].c_str(), "P");

      gMassRatios->Add(massRatioGraph[filecounter], "PE");
      legend5->AddEntry(massRatioGraph[filecounter], SPMC_legendNames[j].c_str(), "P");

      //------------------------------------------------------------------------------------------------
    
      // Define a function for the eta energy resolution fit
      TF1 *EresolutionFit = new TF1("EresolutionFit", "sqrt([0]*[0]/x + [1]*[1])", 0.1, 20);
      EresolutionFit->SetParameters(0.1, 0.02);  // Initial guesses for a, b
      EresolutionGraph[filecounter]->Fit(EresolutionFit, "R");

      TCanvas *EresCanvas = new TCanvas("EresCanvas", "Resolution Fit", 800, 600);
      EresolutionGraph[filecounter]->SetTitle("Energy Resolution; p_{T} (GeV/c); #sigma / #mu");
      EresolutionGraph[filecounter]->Draw("APE");
      EresolutionFit->Draw("same");

      TCanvas *EfitParamsCanvas = new TCanvas("EfitParamsCanvas", "Fit Parameters", 800, 600);
      TPaveText *EparamsText = new TPaveText(0.1, 0.7, 0.9, 0.9, "NDC");
      EparamsText->AddText("Eta Fitted Resolution Parameters:");
      EparamsText->AddText(unweighted_legendNames[j].c_str());
      EparamsText->AddText(Form("Stochastic term (a): %.4f", EresolutionFit->GetParameter(0)));
      //paramsText->AddText(Form("Noise term (b): %.4f", EresolutionFit->GetParameter(2)));
      EparamsText->AddText(Form("Constant term (c): %.4f", EresolutionFit->GetParameter(1)));
      //add goodness of fit
      EparamsText->AddText(Form("Chi2/ndf: %.4f", EresolutionFit->GetChisquare() / EresolutionFit->GetNDF()));
      EparamsText->Draw();

      gEResolutions->Add(EresolutionGraph[filecounter], "PE");
      legend7->AddEntry(EresolutionGraph[filecounter], unweighted_legendNames[j].c_str(), "P");

      EresCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
      EresCanvas->Close();
      EfitParamsCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
      EfitParamsCanvas->Close();

    }

    file.Close();
    filecounter++;
  }

  
  std::cout << "Run2024 files: " << Run2024_FileNames.size() << std::endl;
 //data, pre made graphs
  if(Run2024_FileNames.size()>0){
    for (size_t j = 0; j < Run2024_FileNames.size(); ++j){
      // Load the new file and retrieve the TGraphErrors
      TFile newFile(Run2024_FileNames[0].c_str(), "READ");
      if (!newFile.IsOpen())
      {
        std::cerr << "Error opening file: " << Run2024_FileNames[0] << std::endl;
        return;
      }

      TGraphErrors* newPionMean = dynamic_cast<TGraphErrors*>(newFile.Get("gr_mass_pi0"));
      TGraphErrors* newPionSigma = dynamic_cast<TGraphErrors*>(newFile.Get("gr_width_pi0"));
      TGraphErrors* newEtaMean = dynamic_cast<TGraphErrors*>(newFile.Get("gr_mass_eta"));
      TGraphErrors* newEtaSigma = dynamic_cast<TGraphErrors*>(newFile.Get("gr_width_eta"));
      TGraphErrors* pionRelativeWidthGraph = new TGraphErrors();
      TGraphErrors* etaRelativeWidthGraph = new TGraphErrors();
      if (!newPionMean || !newPionSigma || !newEtaMean || !newEtaSigma)
      {
        std::cerr << "Error getting TGraphErrors from new file." << std::endl;
        newFile.Close();
        return;
      }
      
      int nPionPoints = newPionMean->GetN();
      for (int i = 0; i < nPionPoints; ++i) {
          double pionPt, pionMean, pionSigma;
          double pionMeanErr, pionSigmaErr;
          newPionMean->GetPoint(i, pionPt, pionMean);
          pionMeanErr = newPionMean->GetErrorY(i);
          newPionSigma->GetPoint(i, pionPt, pionSigma);
          pionSigmaErr = newPionSigma->GetErrorY(i);

          double pionRelativeWidth = pionMean > 0 ? pionSigma / pionMean : 0;
          double pionRelativeWidthErr = pionRelativeWidth * sqrt(pow(pionMeanErr / pionMean, 2) + pow(pionSigmaErr / pionSigma, 2));

          pionRelativeWidthGraph->SetPoint(i, pionPt, pionRelativeWidth);
          pionRelativeWidthGraph->SetPointError(i, 0, pionRelativeWidthErr);
      }


      int nEtaPoints = newEtaMean->GetN();
      for (int i = 0; i < nEtaPoints; ++i) {
          double etaPt, etaMean, etaSigma;
          double etaMeanErr, etaSigmaErr;
          newEtaMean->GetPoint(i, etaPt, etaMean);
          etaMeanErr = newEtaMean->GetErrorY(i);
          newEtaSigma->GetPoint(i, etaPt, etaSigma);
          etaSigmaErr = newEtaSigma->GetErrorY(i);

          double etaRelativeWidth = etaMean > 0 ? etaSigma / etaMean : 0;
          double etaRelativeWidthErr = etaRelativeWidth * sqrt(pow(etaMeanErr / etaMean, 2) + pow(etaSigmaErr / etaSigma, 2));

          etaRelativeWidthGraph->SetPoint(i, etaPt, etaRelativeWidth);
          etaRelativeWidthGraph->SetPointError(i, 0, etaRelativeWidthErr);
      }
      
      MarkerStyle+=1;
      MarkerColor+=1;
      if(MarkerColor==5 || MarkerColor==10) MarkerColor+=1;//avoid yellow
      //with all off will default to sphenix style
      //newPionMean->SetMarkerStyle(MarkerStyle);
      //newPionMean->SetMarkerColor(MarkerColor);
      //newEtaMean->SetMarkerStyle(MarkerStyle);
      //newEtaMean->SetMarkerColor(MarkerColor);
      //pionRelativeWidthGraph->SetMarkerStyle(MarkerStyle);
      //pionRelativeWidthGraph->SetMarkerColor(MarkerColor);
      //etaRelativeWidthGraph->SetMarkerStyle(MarkerStyle);
      //etaRelativeWidthGraph->SetMarkerColor(MarkerColor);
      // Add the new graphs to the multigraphs
      gPionMeans->Add(newPionMean, "PE");
      legend1->AddEntry(newPionMean, Run2024_legendNames[j].c_str(), "PE");

      gPionWidths->Add(pionRelativeWidthGraph, "PE");
      legend2->AddEntry(pionRelativeWidthGraph, Run2024_legendNames[j].c_str(), "PE");

      gEtaMeans->Add(newEtaMean, "PE");
      legend3->AddEntry(newEtaMean, Run2024_legendNames[j].c_str(), "PE");

      gEtaWidths->Add(etaRelativeWidthGraph, "PE");
      legend4->AddEntry(etaRelativeWidthGraph, Run2024_legendNames[j].c_str(), "PE");

      newFile.Close();
    }
  }

  // draw multigraphs
  TCanvas *c1 = new TCanvas("c1", "Canvas1", 800, 600);
  //gPad->SetFillColor(33);
  gPionMeans->SetTitle("Pion: Smeared pT vs Inv. Mass;#it{pT}_{#gamma#gamma} (GeV); Pion Peak Position (GeV)");
  gPionMeans->GetXaxis()->SetLimits(0.01, 10);
  gPionMeans->SetMinimum(0.135);
  gPionMeans->SetMaximum(0.17);
  gPionMeans->Draw("APE");
  legend1->SetFillStyle(0);
  legend1->SetTextAlign(32);
  legend1->SetTextSize(0.02);
  legend1->Draw();
  c1->Print("pioncode/canvas_pdf/ptdifferentialcomparison.pdf");
  // gPionMeans->GetXaxis()->SetLimits(comparisonFilenameObj.plotxlims[0], comparisonFilenameObj.plotxlims[1]);
  // gPionMeans->SetMinimum(comparisonFilenameObj.plotylims[0]);
  // gPionMeans->SetMaximum(comparisonFilenameObj.plotylims[1]);
  // c1->SetMargin(0.2, 0.1, 0.1, 0.1);
  // gPad->Modified();
  // gPad->Update();
  // c1->SaveAs(Form("%s/%
  TCanvas *c2 = new TCanvas("c2", "Canvas2", 800, 600);
  //gPad->SetFillColor(33);
  gPionWidths->SetTitle("Pion: Smeared pT vs Relative Width;#it{pT}_{#gamma#gamma} (GeV); Pion Relative Width");
  gPionWidths->GetXaxis()->SetLimits(0.01, 10);
  gPionWidths->SetMinimum(0.05);
  gPionWidths->SetMaximum(0.2);
  gPionWidths->Draw("APE");
  legend2->SetFillStyle(0);
  legend2->SetTextAlign(32);
  legend2->SetTextSize(0.02);
  legend2->Draw();
  c2->Print("pioncode/canvas_pdf/ptdifferentialcomparison.pdf");

  TCanvas *c3 = new TCanvas("c3", "Canvas3", 800, 600);
  //gPad->SetFillColor(33);
  gEtaMeans->SetTitle("Eta: Smeared pT vs Inv. Mass;#it{pT}_{#gamma#gamma} (GeV); Eta Peak Position (GeV)");
  gEtaMeans->GetXaxis()->SetLimits(0.01, 17);
  gEtaMeans->SetMinimum(0.45);
  gEtaMeans->SetMaximum(0.7);
  gEtaMeans->Draw("APE");
  legend3->SetFillStyle(0);
  legend3->SetTextAlign(32);
  legend3->SetTextSize(0.02);
  legend3->Draw();
  c3->Print("pioncode/canvas_pdf/ptdifferentialcomparison.pdf");

  TCanvas *c4 = new TCanvas("c4", "Canvas4", 800, 600);
  //gPad->SetFillColor(33);
  gEtaWidths->SetTitle("Eta: Smeared pT vs Relative Width;#it{pT}_{#gamma#gamma} (GeV); Eta Relative Width");
  gEtaWidths->GetXaxis()->SetLimits(0.01, 17);
  gEtaWidths->SetMinimum(0.01);
  gEtaWidths->SetMaximum(0.3);
  gEtaWidths->Draw("APE");
  legend4->SetFillStyle(0);
  legend4->SetTextAlign(32);
  legend4->SetTextSize(0.02);
  legend4->Draw();
  c4->Print("pioncode/canvas_pdf/ptdifferentialcomparison.pdf");

  TCanvas *c5 = new TCanvas("c5", "Canvas5", 800, 600);
  //gPad->SetFillColor(33);
  gMassRatios->SetTitle("Mass Ratios: Smeared pT vs Inv. Mass;#it{pT}_{#gamma#gamma} (GeV); Pion/Eta Mass Ratio");
  gMassRatios->GetXaxis()->SetLimits(0.01, 10);
  gMassRatios->SetMinimum(0.24);
  gMassRatios->Draw("APE");
  legend5->SetFillStyle(0);
  legend5->SetTextAlign(32);
  legend5->SetTextSize(0.02);
  legend5->Draw();
  //c5->Print("pioncode/canvas_pdf/ptdifferentialcomparison.pdf");

  TCanvas *c6 = new TCanvas("c6", "Canvas6", 800, 600);
  //gPad->SetFillColor(33);
  gPResolutions->SetTitle("Pion: Smeared pT vs Resolution;#it{pT}_{#gamma#gamma} (GeV); Pion Resolution");
  gPResolutions->GetXaxis()->SetLimits(0.01, 10);
  gPResolutions->SetMinimum(0.1);
  //gPResolutions->SetMaximum(0.3);
  gPResolutions->Draw("APE");
  legend6->SetFillStyle(0);
  legend6->SetTextAlign(32);
  legend6->SetTextSize(0.02);
  legend6->Draw();
  c6->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");

  TCanvas *c7 = new TCanvas("c7", "Canvas7", 800, 600);
  //gPad->SetFillColor(33);
  gEResolutions->SetTitle("Eta: Smeared pT vs Resolution;#it{pT}_{#gamma#gamma} (GeV); Eta Resolution");
  gEResolutions->GetXaxis()->SetLimits(0.01, 17);
  gEResolutions->SetMinimum(0.04);
  //gEResolutions->SetMaximum(0.3);
  gEResolutions->Draw("APE");
  legend7->SetFillStyle(0);
  legend7->SetTextAlign(32);
  legend7->SetTextSize(0.02);
  legend7->Draw();
  c7->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");


  // Close the PDF file
  dummyCanvas->Print("pioncode/canvas_pdf/ptdifferentialcomparison.pdf]");
  dummyCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf]");

  TFile outputFile("pioncode/rootfiles/ptdifferential_overlay.root", "RECREATE");
  gPionMeans->Write("gPionMean");
  gPionWidths->Write("gPionWidth");
  gEtaMeans->Write("gEtaMean");
  gEtaWidths->Write("gEtaWidth");
  gMassRatios->Write("gMassRatio");
  outputFile.Close();
}

void fit_comparison()
{
  //-----------------------------------------
  std::vector<std::string> unweighted_fileNames = {

    "pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_pythia8_pp_mb_0000000015_merged_V36.root",
    "pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_pythia8_pp_mb_0000000015_merged_V37.root"};//    "pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_pythia8_pp_mb_3MHz_0000000011__merged_V1.root",
  std::vector<std::string> unweighted_histNames = {"h_InvMass_2d", "h_InvMass_smear_2d_100"};//"h_InvMass_2d",
  std::vector<std::string> unweighted_legendNames = {"Pythia_wvfm_E","Pythia_wvfm_E+10%smr"};//"Pythia",

  //-----------------------------------------
  std::vector<std::string> SPMC_FileNames = {
    "pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_single_pi0_p_200_20000MeV_0000000017_00merged_V38.root",
    "pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_single_eta_p_600_20000MeV_0000000017_00merged_V39.root",
    "pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_single_eta_p_600_20000MeV_0000000017_00merged_V40.root",
    "pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_single_pi0_p_200_20000MeV_0000000017_00merged_V41.root"
    };
  //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_single_pi0_p_200_20000MeV_0000000017_00merged_V38OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_200_10000MeV_0000000013_00merged_V13.root"
  //,"pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_200_10000MeV_0000000013_00merged_V14.root"
  std::vector<std::string> SPMC_histNames = {
    "h_InvMass_smear_weighted_2d_0", 
    "h_InvMass_smear_weighted_2d_0",
    "h_InvMass_smear_weighted_2d_125",
    "h_InvMass_smear_weighted_2d_125"};
  std::vector<std::string> SPMC_legend = {"SPi0+0sm","SEta+0sm","SEta+12.5sm","SPi0+12.5sm"};
  std::vector<int> SPMC_FileTypes ={0,1,1,0};//0 for pion, 1 for eta

  //-----------------------------------------
  std::vector<std::string> FastMC_fileNames = {};
  //    "pioncode/rootfiles/PionFastMC_0.154000_sqrte_0.180000_const.root",
  //  "pioncode/rootfiles/EtaFastMC_0.154000_sqrte_0.120000_const.root",
  //  "pioncode/rootfiles/EtaFastMC_0.154000_sqrte_0.150000_const.root",
  //  "pioncode/rootfiles/EtaFastMC_0.154000_sqrte_0.180000_const.root",
  //  "pioncode/rootfiles/PionFastMC_0.154000_sqrte_0.180000_const.root",
  //  "pioncode/rootfiles/PionFastMC_0.154000_sqrte_0.180000_const.root",
  //  "pioncode/rootfiles/PionFastMC_0.154000_sqrte_0.180000_const.root",
  //  "pioncode/rootfiles/PionFastMC_0.154000_sqrte_0.180000_const.root"    
  std::vector<std::string> FastMC_histNames = {"h101_2", "h101_2", "h101_2","h101_2","h101_2_symm_2","h101_2_asymm_2","h101_2_symm_0","h101_2_symm_1"};
  std::vector<std::string> FastMC_legendNames = {
    "FastMC: 15.4%/#sqrt{E} #oplus 18%",
    "FastMC: 15.4%/#sqrt{E} #oplus 12%",
    "FastMC: 15.4%/#sqrt{E} #oplus 15%",
    "FastMC: 15.4%/#sqrt{E} #oplus 18%",
    "FastMC_symm_E: 15.4%/#sqrt{E} #oplus 18%",
    "FastMC_asymm_E: 15.4%/#sqrt{E} #oplus 18%",
    "FastMC_symm_R: 15.4%/#sqrt{E} #oplus 18%",
    "FastMC_symm_T: 15.4%/#sqrt{E} #oplus 18%"
    };//"PionFastMC", "EtaFastMC"
  std::vector<int> FastMC_FileTypes ={0,1,1,1,0,0,0,0};//0 for pion, 1 for eta
  //
  //-----------------------------------------
  std::vector<std::string> Run2024_fileNames = {"pioncode/rootfiles/meson_graphs.root"};
  std::vector<std::string> Run2024_legendNames = {"Run2024"};

  //-----------------------------------------
  AnalyzeHistograms(unweighted_fileNames, unweighted_histNames, unweighted_legendNames, SPMC_FileNames, SPMC_histNames, SPMC_legend, SPMC_FileTypes, FastMC_fileNames, FastMC_histNames, FastMC_legendNames,FastMC_FileTypes,Run2024_fileNames, Run2024_legendNames);


  gApplication->Terminate(0);
  // return 0;
}
