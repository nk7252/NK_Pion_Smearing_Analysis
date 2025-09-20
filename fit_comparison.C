#pragma once
// c++ includes
#include <iostream>
#include <string>
#include <vector>
// root includes
#include <Math/Factory.h>
#include <Math/Functor.h>
#include <Math/MinimizerOptions.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <TLegend.h>
#include <TMultiGraph.h>
#include <TPaveText.h>
#include <TROOT.h>
#include <TStyle.h>
// local includes
#include "sPhenixStyle.C"
#include "sPhenixStyle.h"
// roofit includes
#include "RooFit.h"
#include "RooRealVar.h"
#include "RooArgList.h"
#include "RooArgSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooAddPdf.h"
#include "RooChebychev.h"
#include "RooCrystalBall.h"
#include "RooPlot.h"
#include "RooAbsReal.h"

// (optional) quiet RooFit a bit — top of function/main, once
RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
RooMsgService::instance().setSilentMode(true);

// global bin var
std::vector<double> nuBins = {0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2, 0.21, 0.22, 0.24, 0.26, 0.28, 0.3, 0.32, 0.36, 0.40, 0.44, 0.48, 0.50, 0.52, 0.54, 0.56, 0.58, 0.60, 0.62, 0.64, 0.66, 0.68, 0.70, 0.72, 0.76, 0.8, 0.84, 0.88, 0.92, 0.96, 1.0}; //, 1.04, 1.08, 1.12, 1.16, 1.2
// helepers for roofit fits, unweighted histograms
// ===== helper to fit one projection with RooFit =====
struct RooSliceResult
{
  double pMean = std::numeric_limits<double>::quiet_NaN();
  double pSigma = std::numeric_limits<double>::quiet_NaN();
  double pMeanErr = std::numeric_limits<double>::quiet_NaN();
  double pSigmaErr = std::numeric_limits<double>::quiet_NaN();
  double eMean = std::numeric_limits<double>::quiet_NaN();
  double eSigma = std::numeric_limits<double>::quiet_NaN();
  double eMeanErr = std::numeric_limits<double>::quiet_NaN();
  double eSigmaErr = std::numeric_limits<double>::quiet_NaN();
  double chi2ndf = std::numeric_limits<double>::quiet_NaN();
};

static RooSliceResult FitSliceWithRooFitAndPages(
    TH1D *h,
    double fitMin, double fitMax,
    double muPionSeed,
    double pionLow, double pionHigh,
    double etaLow, double etaHigh,
    int sliceId,
    const char *pdfPath,
    const char *runBanner,
    const char *legendLine,
    double pt_min, double pt_max, double pt_center,
    bool makeSpectrumPage = true, bool makeParamsPage = true, bool makePullPage = false)
{
  using namespace RooFit;
  RooSliceResult out;

  // -------- observable & data
  RooRealVar x(Form("x_%d", sliceId), "m_{#gamma#gamma}", fitMin, fitMax, "GeV");
  RooDataHist data(Form("data_%d", sliceId), "data", RooArgList(x), Import(*h));

  // -------- π0 CrystalBall (BLUE)
  RooRealVar meanP(Form("meanP_%d", sliceId), "#pi^{0} mean", muPionSeed, pionLow, pionHigh);
  RooRealVar sigmaP(Form("sigmaP_%d", sliceId), "#pi^{0} sigma", 0.02, 0.002, 0.09);
  RooRealVar alphaP(Form("alphaP_%d", sliceId), "#pi^{0} alpha", 1.5, 0.2, 8.0);
  RooRealVar nP(Form("nP_%d", sliceId), "#pi^{0} n", 3.0, 1.0, 25.0);
  RooCrystalBall cbP(Form("cbP_%d", sliceId), "CB(#pi^{0})", x, meanP, sigmaP, alphaP, nP);

  // -------- η Gaussian (GREEN)
  RooRealVar meanE(Form("meanE_%d", sliceId), "#eta mean", 0.55, etaLow, etaHigh);
  RooRealVar sigmaE(Form("sigmaE_%d", sliceId), "#eta sigma", 0.05, 0.01, 0.15);
  RooGaussian gE(Form("gE_%d", sliceId), "G(#eta)", x, meanE, sigmaE);

  // -------- background Cheby2 (RED dashed)
  RooRealVar c0(Form("c0_%d", sliceId), "c0", 0.0, -2.0, 2.0);
  RooRealVar c1(Form("c1_%d", sliceId), "c1", 0.0, -2.0, 2.0);
  RooChebychev bkg(Form("bkg_%d", sliceId), "bkg", x, RooArgList(c0, c1));

  // -------- extended yields
  const double nTot = h->GetSumOfWeights();
  RooRealVar Np(Form("Np_%d", sliceId), "N_{#pi^{0}}", 0.35 * nTot, 0.0, 2.0 * nTot);
  RooRealVar Ne(Form("Ne_%d", sliceId), "N_{#eta}", 0.05 * nTot, 0.0, 2.0 * nTot);
  RooRealVar Nb(Form("Nb_%d", sliceId), "N_{bkg}", 0.60 * nTot, 0.0, 2.0 * nTot);

  RooAddPdf model(Form("model_%d", sliceId), "cbP + gE + bkg",
                  RooArgList(cbP, gE, bkg), RooArgList(Np, Ne, Nb));

  // -------- fit
  model.fitTo(data, Extended(kTRUE), SumW2Error(kTRUE), PrintLevel(-1));

  // -------- χ2/ndf (rough)
  RooPlot *fchi = x.frame();
  data.plotOn(fchi);
  model.plotOn(fchi);
  out.chi2ndf = fchi->chiSquare(model.getParameters(data)->getSize());
  delete fchi;

  // -------- export params
  out.pMean = meanP.getVal();
  out.pMeanErr = meanP.getError();
  out.pSigma = sigmaP.getVal();
  out.pSigmaErr = sigmaP.getError();
  out.eMean = meanE.getVal();
  out.eMeanErr = meanE.getError();
  out.eSigma = sigmaE.getVal();
  out.eSigmaErr = sigmaE.getError();

  // -------- PAGE A: spectrum (colors preserved)
  if (makeSpectrumPage)
  {
    RooPlot *fr = x.frame(Title("CB(#pi^{0}) + G(#eta) + Cheby2(bkg)"));
    data.plotOn(fr);
    // order: draw components first (so total draws last on top)
    model.plotOn(fr, Components(bkg), LineColor(kRed), LineStyle(kDashed));
    model.plotOn(fr, Components(cbP), LineColor(kBlue), LineStyle(kSolid));
    model.plotOn(fr, Components(gE), LineColor(kGreen), LineStyle(kSolid));
    model.plotOn(fr, LineColor(kBlack), LineStyle(kSolid));

    TCanvas *c = new TCanvas(Form("cSpec_%d", sliceId), "", 800, 600);
    fr->GetXaxis()->SetTitle("#it{m}_{#gamma#gamma} (GeV)");
    fr->GetYaxis()->SetTitle("#frac{1}{2 #it{p}_{T}} #frac{dN}{d#it{m}_{#gamma#gamma}}");
    fr->Draw();

    TLatex lt;
    lt.SetNDC(true);
    lt.SetTextSize(0.035);
    lt.DrawLatex(0.55, 0.92, "#bf{sPHENIX} Internal");
    if (runBanner && *runBanner)
      lt.DrawLatex(0.55, 0.88, runBanner);
    lt.DrawLatex(0.55, 0.84, Form("p_{T}: %.2f - %.2f GeV  (center: %.2f)", pt_min, pt_max, pt_center));

    // Legend with dummy lines (stable across RooFit versions)
    TLegend *leg = new TLegend(0.55, 0.60, 0.90, 0.83);
    leg->SetFillStyle(0);
    if (legendLine && *legendLine)
      leg->AddEntry((TObject *)0, legendLine, "");
    auto lModel = new TLine();
    lModel->SetLineColor(kBlack);
    auto lPi = new TLine();
    lPi->SetLineColor(kBlue);
    auto lEta = new TLine();
    lEta->SetLineColor(kGreen);
    auto lBkg = new TLine();
    lBkg->SetLineColor(kRed);
    lBkg->SetLineStyle(kDashed);
    leg->AddEntry(lModel, "Model", "l");
    leg->AddEntry(lPi, "#pi^{0} CrystalBall", "l");
    leg->AddEntry(lEta, "#eta Gaussian", "l");
    leg->AddEntry(lBkg, "Background (Cheby2)", "l");
    leg->Draw();

    c->Print(pdfPath);
    delete leg;
    delete lModel;
    delete lPi;
    delete lEta;
    delete lBkg;
    delete c;
    delete fr;
  }

  // -------- PAGE B: parameter summary
  if (makeParamsPage)
  {
    const double PWidth = out.pSigma / out.pMean;
    const double EWidth = out.eSigma / out.eMean;

    TCanvas *c = new TCanvas(Form("cPar_%d", sliceId), "", 800, 600);
    TPaveText *txt = new TPaveText(0.10, 0.70, 0.90, 0.90, "NDC");
    txt->SetFillStyle(0);
    if (legendLine && *legendLine)
      txt->AddText(legendLine);
    txt->AddText(Form("p_{T}: %.2f - %.2f GeV  (center: %.2f)", pt_min, pt_max, pt_center));
    txt->AddText(Form("#chi^{2}/ndf: %.3f", out.chi2ndf));
    txt->AddText(Form("#pi^{0} mean = %.5f #pm %.5f", out.pMean, out.pMeanErr));
    txt->AddText(Form("#pi^{0} #sigma = %.5f #pm %.5f", out.pSigma, out.pSigmaErr));
    txt->AddText(Form("#pi^{0} rel. width = %.3f%%", 100.0 * PWidth));
    txt->AddText(Form("#eta mean = %.5f #pm %.5f", out.eMean, out.eMeanErr));
    txt->AddText(Form("#eta #sigma = %.5f #pm %.5f", out.eSigma, out.eSigmaErr));
    txt->AddText(Form("#eta rel. width = %.3f%%", 100.0 * EWidth));
    txt->Draw();
    c->Print(pdfPath);
    delete txt;
    delete c;
  }

  // -------- PAGE C: pull plot (optional)
  if (makePullPage)
  {
    // Build pulls using same model/data
    RooPlot *base = x.frame();
    data.plotOn(base);
    model.plotOn(base);
    RooHist *hpull = base->pullHist();

    RooPlot *fp = x.frame(Title("Pull distribution"));
    fp->addPlotable(hpull, "P");

    TCanvas *c = new TCanvas(Form("cPull_%d", sliceId), "", 800, 600);
    fp->GetYaxis()->SetTitle("pull");
    fp->Draw();
    c->Print(pdfPath);

    delete c;
    delete fp;
    delete hpull;
    delete base;
  }

  return out;
}

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

Double_t myGausPol2(Double_t *x, Double_t *par)
{
  // x[0] is the current x-value
  // par[0] = Gaussian amplitude
  // par[1] = Gaussian mean
  // par[2] = Gaussian sigma
  // par[3], par[4], par[5] = pol2 coefficients
  double gausAmp = par[0];
  double gausMean = par[1];
  double gausSigma = par[2];

  // Gaussian piece
  double arg = (x[0] - gausMean) / gausSigma;
  double gausVal = gausAmp * TMath::Exp(-0.5 * arg * arg);

  // pol2 piece:  par[3] + par[4]*x + par[5]*x^2
  double polyVal = par[3] + par[4] * x[0] + par[5] * x[0] * x[0];

  return gausVal + polyVal;
}

Double_t myGausPol3(Double_t *x, Double_t *par)
{
  // par[0] = amplitude
  // par[1] = mean
  // par[2] = sigma
  // par[3..6] = pol3 coefficients

  double gausAmp = par[0];
  double gausMean = par[1];
  double gausSigma = par[2];

  double arg = (x[0] - gausMean) / gausSigma;
  double gausVal = gausAmp * TMath::Exp(-0.5 * arg * arg);

  // pol3: par[3] + par[4]*x + par[5]*x^2 + par[6]*x^3
  double polyVal = par[3] + par[4] * x[0] + par[5] * x[0] * x[0] + par[6] * x[0] * x[0] * x[0];

  return gausVal + polyVal;
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

double poly2BG(double *x, double *par)
{
  // 2nd degree polynomial background
  // Check if x is in the range of any Gaussian fit
  /*
  if (x[0] >= par[3] && x[0] <= par[4])
  {
    TF1::RejectPoint();
    return 0;
  }
  //*/
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0];
}
double poly2BGCut(double *x, double *par)
{
  // 2nd degree polynomial background
  // Check if x is in the range of any Gaussian fit
  ///*
  if (x[0] >= par[3] && x[0] <= par[4])
  {
    TF1::RejectPoint();
    return 0;
  }
  //*/
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0];
}

double poly3BG(double *x, double *par)
{
  // 3rd degree polynomial background
  // Check if x is in the range of any Gaussian fit
  /*
  if (x[0] >= 0.52 && x[0] <= 0.68)
  {
    TF1::RejectPoint();
    return 0;
  }
  */
  return par[0] + par[1] * x[0] + par[2] * x[0] * x[0] + par[3] * x[0] * x[0] * x[0];
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

void AnalyzeHistograms(const std::vector<std::string> &unweightedFileNames, const std::vector<std::string> &unweightedhistNames, const std::vector<std::string> &unweighted_legendNames, const std::vector<std::string> &SPMC_FileNames, const std::vector<std::string> &SPMC_histNames, const std::vector<std::string> &SPMC_legendNames, std::vector<int> SPMC_FileTypes, const std::vector<std::string> &FastMC_FileNames, const std::vector<std::string> &FastMC_histNames, const std::vector<std::string> &FastMC_legendNames, std::vector<int> FastMC_FileTypes, const std::vector<std::string> &Run2024_FileNames, const std::vector<std::string> &Run2024_legendNames)
{
  // be ready to try different minimizers...
  //  possible choices are:
  //      minName                  algoName
  //  Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
  //   Minuit2                     Fumili2
  //   Fumili
  //   GSLMultiMin                ConjugateFR, ConjugatePR, BFGS,
  //                               BFGS2, SteepestDescent
  //   GSLMultiFit
  //    GSLSimAn
  //    Genetic
  ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit2", "Combined"); //,"Simplex", "Migrad", "Fumili"
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("GSLMultiMin");//, "ConjugateFR"
  // ROOT::Math::MinimizerOptions::SetDefaultMinimizer("GSLMultiFit");//, "LevenbergMarquardt"
  // ROOT::Math::MinimizerOptions::SetDefaultAlgorithm("Fumili2");

  ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);
  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(1000000);
  ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(10000);
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(0.001);
  ROOT::Math::MinimizerOptions::SetDefaultPrecision(1e-12);
  // ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(3);
  SetsPhenixStyle();
  //
  bool var_bins = true;
  bool var_bins_unw = false;
  int rebinFactor = 2;
  bool dynamic_left = true;
  int startBin = 5;
  int endBin_global = -1;
  int projectionBins = 4;
  double scale_factor = 1.0;
  double scale_factor_unw = 1.0; // error scale up factor
  double limits[10] = {0.05, 1.0, 0.09, 0.25, 0.05, 0.35, 0.52, 0.68, 0.35, 1.0};
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
  dummyCanvas->Print("pioncode/canvas_pdf/ptdifferential_PionFit_results.pdf[");
  dummyCanvas->Print("pioncode/canvas_pdf/ptdifferential_EtaFit_results.pdf[");
  dummyCanvas->Print("pioncode/canvas_pdf/ptdifferential_unw_Fit_results.pdf[");
  // top right (0.66, 0.7, 0.90, 0.9)
  // top left (0.2, 0.7, 0.44, 0.9)
  // bottom left (0.2, 0.2, 0.44, 0.4)
  // bottom right (0.66, 0.2, 0.90, 0.4)
  TLegend *legend1 = new TLegend(0.2, 0.7, 0.44, 0.9);  // pion mean
  TLegend *legend2 = new TLegend(0.66, 0.7, 0.90, 0.9); // pion width
  TLegend *legend3 = new TLegend(0.66, 0.7, 0.90, 0.9); // eta mean
  TLegend *legend4 = new TLegend(0.66, 0.7, 0.90, 0.9); // eta width
  TLegend *legend5 = new TLegend(0.2, 0.7, 0.44, 0.9);  // mass ratio
  TLegend *legend6 = new TLegend(0.66, 0.7, 0.90, 0.9); // pion resolution
  TLegend *legend7 = new TLegend(0.66, 0.7, 0.90, 0.9); // eta resolution
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
  int totalfiles = unweightedFileNames.size() + FastMC_FileNames.size() + SPMC_FileNames.size() + Run2024_FileNames.size();
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
  // Each entry corresponds to:
  //{pt_min, pt_max, pionFitLow, pionFitHigh, bgLeft, bgRight}
  std::vector<std::tuple<double, double, double, double, double, double>> pionFitRanges = {
      {1.0, 2.0, 0.11, 0.16, 0.07, 0.25}, // pT: 2.0 - 2.5 GeV
      {2.0, 3.0, 0.12, 0.16, 0.05, 0.22},
      {3.0, 4.0, 0.12, 0.16, 0.05, 0.22},
      {4.0, 5.0, 0.11, 0.16, 0.05, 0.25},
      {5.0, 6.0, 0.11, 0.16, 0.07, 0.25},
      {6.0, 7.0, 0.11, 0.16, 0.09, 0.21},
      {7.0, 8.0, 0.12, 0.16, 0.09, 0.20},
      {8.0, 9.0, 0.12, 0.16, 0.1, 0.19},
      {9.0, 10.0, 0.12, 0.17, 0.1, 0.19},
      {10.0, 11.0, 0.21, 0.3, 0.18, 0.33},
      {11.0, 12.0, 0.21, 0.33, 0.2, 0.35},
      {12.0, 13.0, 0.21, 0.35, 0.2, 0.35},
      {13.0, 14.0, 0.21, 0.36, 0.2, 0.35},
      {14.0, 15.0, 0.3, 0.4, 0.25, 0.45},
      {15.0, 16.0, 0.3, 0.4, 0.25, 0.45},
      {16.0, 17.0, 0.3, 0.4, 0.25, 0.45},
      {17.0, 18.0, 0.3, 0.4, 0.25, 0.45},
      {18.0, 19.0, 0.11, 0.19, 0.05, 0.17},
      {19.0, 20.0, 0.11, 0.19, 0.05, 0.17},
  };

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
      if (var_bins_unw && !nuBins.empty())
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
      // float leftmost_limit = 0;
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
      // Scale bin contents by pT, avoiding division by zero
      for (int i = 1; i <= histF->GetNbinsX(); i++)
      {
        if (pion_pt == 0)
          continue;                                             // Avoid division by zero
        double binContent = histF->GetBinContent(i);            // Get the bin content
        double newContent = binContent / pion_pt;               // Scale the bin content
        histF->SetBinContent(i, newContent);                    // Set the new bin content
        histF->SetBinError(i, histF->GetBinError(i) / pion_pt); // Scale the bin error
      }

      // estimate the peak
      double wideMin = 0.05;
      double wideMax = 0.45;

      // Quick estimate: find maximum bin in [wideMin, wideMax]
      int iMin = histF->FindBin(wideMin);
      int iMax = histF->FindBin(wideMax);

      int maxBin = iMin;
      double maxVal = 0.0;
      for (int b = iMin; b <= iMax; b++)
      {
        double c = histF->GetBinContent(b);
        if (c > maxVal)
        {
          maxVal = c;
          maxBin = b;
        }
      }
      // Approx mean
      double muPeak = histF->GetBinCenter(maxBin);

      double pionFitLow = 0.10, pionFitHigh = 0.16; // Default pion peak fit range
      double bgLeft = 0.05, bgRight = 0.25;         // Default background fit range

      for (const auto &range : pionFitRanges)
      {
        if (pion_pt >= std::get<0>(range) && pion_pt < std::get<1>(range))
        {
          pionFitLow = std::get<2>(range);
          pionFitHigh = std::get<3>(range);
          bgLeft = std::get<4>(range);
          bgRight = std::get<5>(range);
          break;
        }
      }
      // ---- choose a global fit window; you were using pieces, this is simpler:
      double fitMin = std::max(0.05, bgLeft); // left boundary you already compute
      double fitMax = 0.78;                   // as before (or limits[1])

      // Right-peak seed if you want to be fancy: quick scan in [0.45, 0.65]
      int iEtaMin = histF->FindBin(0.45);
      int iEtaMax = histF->FindBin(0.65);
      int maxBinEta = iEtaMin;
      double maxEta = 0;
      for (int b = iEtaMin; b <= iEtaMax; ++b)
      {
        double v = histF->GetBinContent(b);
        if (v > maxEta)
        {
          maxEta = v;
          maxBinEta = b;
        }
      }
      double muEtaSeed = histF->GetBinCenter(maxBinEta);

      // Call RooFit (sliceId = unique integer; use e.g. (filecounter*10000 + bincounter))

      int sliceId = filecounter * 10000 + bincounter;
      RooSliceResult R = FitSliceWithRooFitAndPages(
          histF, fitMin, fitMax,
          muPeak, pionFitLow, pionFitHigh,
          0.5, 0.65,
          sliceId,
          "pioncode/canvas_pdf/ptdifferential_unw_Fit_results.pdf",
          "run2024:p+p  #sqrt{s_{NN}} = 200 GeV",
          unweighted_legendNames[j].c_str(),
          pt_min, pt_max, pion_pt,
          /*makeSpectrumPage=*/true,
          /*makeParamsPage=*/true,
          /*makePullPage=*/false // toggle on if desired
      );

      // Pull results (just like you did with TF1s)
      double Pmean = R.pMean;
      double Psigma = R.pSigma;
      double PmeanErr = R.pMeanErr;
      double PsigmaErr = R.pSigmaErr;

      double Emean = R.eMean;
      double Esigma = R.eSigma;
      double EmeanErr = R.eMeanErr;
      double EsigmaErr = R.eSigmaErr;

      // Derived quantities (unchanged)
      double PWidth = Psigma / Pmean;
      double PWidthErr = PsigmaErr / Pmean + Psigma * PmeanErr / (Pmean * Pmean);
      double EWidth = Esigma / Emean;
      double EWidthErr = EsigmaErr / Emean + Esigma * EmeanErr / (Emean * Emean);
      double MassRatio = Pmean / Emean;
      double MassRatioErr = MassRatio * std::sqrt(std::pow(PmeanErr / Pmean, 2) + std::pow(EmeanErr / Emean, 2));

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
    MarkerStyle += 1;
    MarkerColor += 1;
    if (MarkerColor == 5 || MarkerColor == 10)
      MarkerColor += 1; // avoid yellow
    if (MarkerStyle == 26 || MarkerStyle == 32)
      MarkerStyle += 1; // avoid triangles
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
    TF1 *PresolutionFit = new TF1("PresolutionFit", "sqrt(2)*sqrt([0]*[0]/x + [1]*[1])", 2, 5.5);
    PresolutionFit->SetParameters(0.154, 0.02); // Initial guesses for a, b
    // PresolutionFit->SetParLimits(0, 0.1, 0.18);
    // PresolutionFit->SetParLimits(1, 0.01, 0.2);
    PresolutionFit->SetLineColor(MarkerColor);
    PresolutionGraph[filecounter]->Fit(PresolutionFit, "RE"); // Fit and constrain to the range of pT

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
    PparamsText->AddText(Form("A/#sqrt{E} term (a): %.4f", PresolutionFit->GetParameter(0)));
    // paramsText->AddText(Form("Noise term (b): %.4f", PresolutionFit->GetParameter(2)));
    PparamsText->AddText(Form("Constant term (c): %.4f", PresolutionFit->GetParameter(1)));
    // add goodness of fit
    PparamsText->AddText(Form("Chi2/ndf: %.4f", PresolutionFit->GetChisquare() / PresolutionFit->GetNDF()));
    PparamsText->Draw();

    gPResolutions->Add(PresolutionGraph[filecounter], "PE");
    legend6->AddEntry(PresolutionGraph[filecounter], unweighted_legendNames[j].c_str(), "P");

    // Save the plot to the PDF
    PresCanvas->SaveAs("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
    // PresCanvas->Close();
    PfitParamsCanvas->SaveAs("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
    // PfitParamsCanvas->Close();

    //------------------------------------------------------------------------------------------------

    // Define a function for the eta energy resolution fit
    TF1 *EresolutionFit = new TF1("EresolutionFit", "sqrt(2)*sqrt([0]*[0]/x + [1]*[1])", 2, 8);
    EresolutionFit->SetParameters(0.154, 0.02); // Initial guesses for a, b
    EresolutionFit->SetParLimits(0, 0.1, 0.18);
    EresolutionFit->SetParLimits(1, 0.01, 0.2);
    EresolutionFit->SetLineColor(MarkerColor);
    EresolutionGraph[filecounter]->Fit(EresolutionFit, "RE");

    TCanvas *EresCanvas = new TCanvas("EresCanvas", "Resolution Fit", 800, 600);
    EresolutionGraph[filecounter]->SetTitle("Energy Resolution; p_{T} (GeV/c); #sigma / #mu");
    EresolutionGraph[filecounter]->Draw("APE");
    EresolutionFit->Draw("same");

    TCanvas *EfitParamsCanvas = new TCanvas("EfitParamsCanvas", "Fit Parameters", 800, 600);
    TPaveText *EparamsText = new TPaveText(0.1, 0.7, 0.9, 0.9, "NDC");
    EparamsText->AddText("Eta Fitted Resolution Parameters:");
    EparamsText->AddText(unweighted_legendNames[j].c_str());
    EparamsText->AddText(Form("A/#sqrt{E} term (a): %.4f", EresolutionFit->GetParameter(0)));
    // paramsText->AddText(Form("Noise term (b): %.4f", EresolutionFit->GetParameter(2)));
    EparamsText->AddText(Form("Constant term (c): %.4f", EresolutionFit->GetParameter(1)));
    // add goodness of fit
    EparamsText->AddText(Form("Chi2/ndf: %.4f", EresolutionFit->GetChisquare() / EresolutionFit->GetNDF()));
    EparamsText->Draw();

    gEResolutions->Add(EresolutionGraph[filecounter], "PE");
    legend7->AddEntry(EresolutionGraph[filecounter], unweighted_legendNames[j].c_str(), "P");

    EresCanvas->SaveAs("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
    // EresCanvas->Close();
    EfitParamsCanvas->SaveAs("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
    // EfitParamsCanvas->Close();
    // finished file
    std::cout << "Finished file: " << unweightedFileNames[j] << std::endl;
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
    if (SPMC_FileTypes[j] == 0) // pion
    {
      for (int i = startBin; i <= endBin; i += projectionBins)
      {
        // Project the histogram along the Y-axis
        int lastBin = std::min(i + projectionBins - 1, nBinsX);
        TH1D *yProjection = hist2D->ProjectionY(Form("proj_%d", i), i, lastBin);
        // Check if the projection has enough entries to perform a fit
        if (yProjection->GetEntries() < 300)
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
        // float leftmost_limit = 0;
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

        double pionFitLow = 0.10, pionFitHigh = 0.16; // Default pion peak fit range
        double bgLeft = 0.05, bgRight = 0.25;         // Default background fit range

        for (const auto &range : pionFitRanges)
        {
          if (pion_pt >= std::get<0>(range) && pion_pt < std::get<1>(range))
          {
            pionFitLow = std::get<2>(range);
            pionFitHigh = std::get<3>(range);
            bgLeft = std::get<4>(range);
            bgRight = std::get<5>(range);
            break;
          }
        }

        // Fit Pion Gaussian in the specified range
        TF1 *gausFit = new TF1("gausFit", "gaus", pionFitLow, pionFitHigh); // limits[2], limits[3]
        // gausFit->SetParLimits(1, 0.11, 0.19);
        // gausFit->SetParLimits(2, 0.01, 0.25);
        gausFit->SetNpx(1000);
        histF->Fit(gausFit, "REQ");

        // Get the fit parameters
        double Pmean = gausFit->GetParameter(1);
        double Psigma = gausFit->GetParameter(2);
        double PmeanErr = gausFit->GetParError(1);
        double PsigmaErr = gausFit->GetParError(2);
        double PWidth = Psigma / Pmean;
        double PWidthErr = PsigmaErr / Pmean + Psigma * PmeanErr / pow(Pmean, 2);
        // double PWidthErr = PWidth * sqrt(pow(PmeanErr / Pmean, 2) + pow(PsigmaErr / Psigma, 2));

        // Pion_Mean.push_back(Pmean);
        // Pion_Width.push_back(PWidth);
        // Pion_Mean_errors.push_back(PmeanErr);
        // Pion_Width_errors.push_back(PWidthErr);
        // pT_Bins.push_back(pion_pt);
        // pT_Bins_Errors.push_back(0);

        TCanvas *tempcanvas = new TCanvas("tempcanvas", "tempcanvas", 800, 600);
        histF->Draw();
        // tempcanvas->Print("pioncode/canvas_pdf/ptdifferential_Fit_results.pdf");
        float xbleft = 0.42;
        float ybleft = 0.7;
        float xtright = 0.9;
        float ytright = 0.9;
        TPaveText *pt2 = new TPaveText(xbleft + .2, 0.5, xtright, 0.7, "NDC"); // Adjust coordinates as needed
        pt2->SetFillColor(0);                                                  // Set the fill color to 0 for transparency
        pt2->SetFillStyle(0);                                                  // Set fill style to 0 (solid) with color 0 for transparency
        pt2->AddText("Pion Fitted Resolution Parameters:");
        pt2->AddText(SPMC_legendNames[j].c_str());
        pt2->AddText(Form("pt region (bin center): %.2f-%.2f GeV (%.2f)", pt_min, pt_max, pion_pt));
        pt2->AddText(Form("#chi^{2}/NDF = %.2f", gausFit->GetChisquare() / gausFit->GetNDF()));
        pt2->AddText(Form("Mean = %.4f", gausFit->GetParameter(1)));
        pt2->AddText(Form("Sigma = %.4f", gausFit->GetParameter(2)));
        pt2->AddText(Form("Relative Width: %.2f%%", gausFit->GetParameter(2) * 100.0f / gausFit->GetParameter(1)));
        pt2->Draw("SAME");
        gPad->Modified(); // Apply the changes to the pad
        tempcanvas->Print("pioncode/canvas_pdf/ptdifferential_PionFit_results.pdf");

        pionmeanGraph[filecounter]->SetPoint(bincounter, pion_pt, Pmean);
        pionmeanGraph[filecounter]->SetPointError(bincounter, 0, PmeanErr);
        pionwidthGraph[filecounter]->SetPoint(bincounter, pion_pt, PWidth);
        pionwidthGraph[filecounter]->SetPointError(bincounter, 0, PWidthErr);
        PresolutionGraph[filecounter]->SetPoint(bincounter, pion_pt, PWidth);
        PresolutionGraph[filecounter]->SetPointError(bincounter, 0, PWidthErr);
        bincounter++;
      }

      MarkerStyle += 1;
      MarkerColor += 1;
      if (MarkerColor == 5 || MarkerColor == 10)
        MarkerColor += 1; // avoid yellow
      if (MarkerStyle == 26 || MarkerStyle == 32)
        MarkerStyle += 1; // avoid triangles
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
      TF1 *PresolutionFit = new TF1("PresolutionFit", "sqrt(2)*sqrt([0]*[0]/x + [1]*[1])", 2, 5.5);
      PresolutionFit->SetParameters(0.154, 0.02); // Initial guesses for a, b
      // PresolutionFit->SetParLimits(0, 0.1, 0.18);
      // PresolutionFit->SetParLimits(1, 0.01, 0.2);
      PresolutionFit->SetLineColor(MarkerColor);
      PresolutionGraph[filecounter]->Fit(PresolutionFit, "RE"); // Fit and constrain to the range of pT
      /*
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
      PparamsText->AddText(Form("A/#sqrt{E} term (a): %.4f", PresolutionFit->GetParameter(0)));
      // paramsText->AddText(Form("Noise term (b): %.4f", PresolutionFit->GetParameter(2)));
      PparamsText->AddText(Form("Constant term (c): %.4f", PresolutionFit->GetParameter(1)));
      // add goodness of fit
      PparamsText->AddText(Form("Chi2/ndf: %.4f", PresolutionFit->GetChisquare() / PresolutionFit->GetNDF()));
      PparamsText->Draw();
      //*/
      gPResolutions->Add(PresolutionGraph[filecounter], "PE");
      legend6->AddEntry(PresolutionGraph[filecounter], SPMC_legendNames[j].c_str(), "P");

      // Save the plot to the PDF
      // PresCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
      // PresCanvas->Close();
      // PfitParamsCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
      // PfitParamsCanvas->Close();
    }
    else if (SPMC_FileTypes[j] == 1) // eta
    {
      for (int i = startBin; i <= endBin; i += projectionBins)
      {
        // Project the histogram along the Y-axis
        int lastBin = std::min(i + projectionBins - 1, nBinsX);
        TH1D *yProjection = hist2D->ProjectionY(Form("proj_%d", i), i, lastBin);
        // Check if the projection has enough entries to perform a fit
        if (yProjection->GetEntries() < 300)
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
        // float leftmost_limit = 0;
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
        float lower_limit = 0.46; // 0.45,0.51
        float upper_limit = 0.61; // 0.70,0.66
        /*
        if (Eta_pt >= 5.5 && Eta_pt <= 6.5){
          lower_limit = 0.45;
          upper_limit = 0.85;
        }
        else if (Eta_pt > 16.5){
          lower_limit = 0.51;
          upper_limit = 0.65;
        }
        */

        // Fit Eta Gaussian in the specified range
        TF1 *EtagausFit = new TF1("EtagausFit", "gaus", lower_limit, upper_limit); // limits[6], limits[7]
        EtagausFit->SetParLimits(1, 0.50, 0.65);
        EtagausFit->SetParLimits(2, 0.03, 0.20);
        EtagausFit->SetNpx(1000);
        histF->Fit(EtagausFit, "REQ");

        // Check if the fit returns NaN or Inf
        bool fitFailed = false;
        for (int i = 0; i < EtagausFit->GetNpar(); i++)
        {
          double param = EtagausFit->GetParameter(i);
          if (std::isnan(param) || std::isinf(param))
          {
            fitFailed = true;
            break;
          }
        }
        if (fitFailed)
        {
          std::cout << "Fit returned NaN or Inf for slice: " << i << std::endl;
          continue;
        }
        TCanvas *tempcanvas = new TCanvas("tempcanvas", "tempcanvas", 800, 600);
        histF->Draw();
        // tempcanvas->Print("pioncode/canvas_pdf/ptdifferential_Fit_results.pdf");
        float xbleft = 0.42;
        float ybleft = 0.7;
        float xtright = 0.9;
        float ytright = 0.9;
        TPaveText *pt2 = new TPaveText(xbleft + .2, 0.5, xtright, 0.7, "NDC"); // Adjust coordinates as needed
        pt2->SetFillColor(0);                                                  // Set the fill color to 0 for transparency
        pt2->SetFillStyle(0);                                                  // Set fill style to 0 (solid) with color 0 for transparency
        pt2->AddText("Eta Fitted Resolution Parameters:");
        pt2->AddText(SPMC_legendNames[j].c_str());
        pt2->AddText(Form("pt region (bin center): %.2f-%.2f GeV (%.2f)", pt_min, pt_max, Eta_pt));
        pt2->AddText(Form("#chi^{2}/NDF = %.2f", EtagausFit->GetChisquare() / EtagausFit->GetNDF()));
        pt2->AddText(Form("Mean = %.4f", EtagausFit->GetParameter(1)));
        pt2->AddText(Form("Sigma = %.4f", EtagausFit->GetParameter(2)));
        pt2->AddText(Form("Relative Width: %.2f%%", EtagausFit->GetParameter(2) * 100.0f / EtagausFit->GetParameter(1)));
        pt2->Draw("SAME");
        gPad->Modified(); // Apply the changes to the pad
        tempcanvas->Print("pioncode/canvas_pdf/ptdifferential_EtaFit_results.pdf");
        double Emean = EtagausFit->GetParameter(1);
        double Esigma = EtagausFit->GetParameter(2);
        double EmeanErr = EtagausFit->GetParError(1);
        double EsigmaErr = EtagausFit->GetParError(2);
        double EWidth = Esigma / Emean;
        double EWidthErr = EsigmaErr / Emean + Esigma * EmeanErr / pow(Emean, 2);
        // EWidth * sqrt(pow(EmeanErr / Emean, 2) + pow(EsigmaErr / Esigma, 2));

        // double MassRatio = Pion_Mean[j] / Emean;
        // double MassRatioErr = MassRatio * sqrt(pow(Pion_Mean_errors[j] / Pion_Mean[j], 2) + pow(EmeanErr / Emean, 2));

        // Eta_Mean.push_back(Emean);
        // Eta_Width.push_back(EWidth);
        // Eta_Mean_errors.push_back(EmeanErr);
        // Eta_Width_errors.push_back(EWidthErr);
        // pT_Bins.push_back(Eta_pt);
        // pT_Bins_Errors.push_back(0);

        etameanGraph[filecounter]->SetPoint(bincounter, Eta_pt, Emean);
        etameanGraph[filecounter]->SetPointError(bincounter, 0, EmeanErr);
        etawidthGraph[filecounter]->SetPoint(bincounter, Eta_pt, EWidth);
        etawidthGraph[filecounter]->SetPointError(bincounter, 0, EWidthErr);
        // massRatioGraph[filecounter]->SetPoint(bincounter, Eta_pt, MassRatio);
        // massRatioGraph[filecounter]->SetPointError(bincounter, 0, MassRatioErr);
        EresolutionGraph[filecounter]->SetPoint(bincounter, Eta_pt, EWidth);
        EresolutionGraph[filecounter]->SetPointError(bincounter, 0, EWidthErr);

        bincounter++;
        // std::cout << "Eta_pt: " << Eta_pt << " Emean: " << Emean << " EWidth: " << EWidth << std::endl;
        // std::cout << "Bincounter: " << bincounter << std::endl;
      }

      MarkerStyle += 1;
      MarkerColor += 1;
      if (MarkerColor == 5 || MarkerColor == 10)
        MarkerColor += 1; // avoid yellow
      if (MarkerStyle == 26 || MarkerStyle == 32)
        MarkerStyle += 1; // avoid triangles
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

      // massRatioGraph[filecounter]->SetMarkerStyle(MarkerStyle);
      // massRatioGraph[filecounter]->SetMarkerColor(MarkerColor);

      EresolutionGraph[filecounter]->SetMarkerStyle(MarkerStyle);
      EresolutionGraph[filecounter]->SetMarkerColor(MarkerColor);

      gEtaMeans->Add(etameanGraph[filecounter], "PE");
      legend3->AddEntry(etameanGraph[filecounter], SPMC_legendNames[j].c_str(), "P");

      gEtaWidths->Add(etawidthGraph[filecounter], "PE");
      legend4->AddEntry(etawidthGraph[filecounter], SPMC_legendNames[j].c_str(), "P");

      // gMassRatios->Add(massRatioGraph[filecounter], "PE");
      // legend5->AddEntry(massRatioGraph[filecounter], FastMC_legendNames[j].c_str(), "P");

      //------------------------------------------------------------------------------------------------

      // Define a function for the eta energy resolution fit
      TF1 *EresolutionFit = new TF1("EresolutionFit", "sqrt(2)*sqrt([0]*[0]/x + [1]*[1])", 2, 15);
      EresolutionFit->SetParameters(0.154, 0.02); // Initial guesses for a, b
      EresolutionFit->SetParLimits(0, 0.01, 0.3);
      EresolutionFit->SetParLimits(1, 0.01, 0.3);
      EresolutionGraph[filecounter]->Fit(EresolutionFit, "RE");
      EresolutionFit->SetLineColor(MarkerColor);
      ///*
      TCanvas *EresCanvas = new TCanvas("EresCanvas", "Resolution Fit", 800, 600);
      EresolutionGraph[filecounter]->SetTitle("Energy Resolution; p_{T} (GeV/c); #sigma / #mu");
      EresolutionGraph[filecounter]->Draw("APE");
      EresolutionFit->Draw("same");

      TCanvas *EfitParamsCanvas = new TCanvas("EfitParamsCanvas", "Fit Parameters", 800, 600);
      TPaveText *EparamsText = new TPaveText(0.1, 0.7, 0.9, 0.9, "NDC");
      EparamsText->AddText("Eta Fitted Resolution Parameters:");
      EparamsText->AddText(SPMC_legendNames[j].c_str());
      EparamsText->AddText(Form("A/#sqrt{E} term (a): %.4f", EresolutionFit->GetParameter(0)));
      // paramsText->AddText(Form("Noise term (b): %.4f", EresolutionFit->GetParameter(2)));
      EparamsText->AddText(Form("Constant term (c): %.4f", EresolutionFit->GetParameter(1)));
      // add goodness of fit
      EparamsText->AddText(Form("Chi2/ndf: %.4f", EresolutionFit->GetChisquare() / EresolutionFit->GetNDF()));
      EparamsText->Draw();
      //*/
      gEResolutions->Add(EresolutionGraph[filecounter], "PE");
      legend7->AddEntry(EresolutionGraph[filecounter], SPMC_legendNames[j].c_str(), "P");

      EresCanvas->SaveAs("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
      // EresCanvas->Close();
      EfitParamsCanvas->SaveAs("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
      // EfitParamsCanvas->Close();
    }
    std::cout << "Finished file: " << SPMC_FileNames[j] << std::endl;
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
    // std::cout << "nBinsX: " << nBinsX << std::endl;
    if (endBin_global == -1)
      endBin = nBinsX; // Default to the last bin if not specified
    std::cout << "endBin: " << endBin << std::endl;
    // Loop over the x-axis bins
    int bincounter = 1;
    if (FastMC_FileTypes[j] == 0)
    { // pion
      for (int i = startBin; i <= endBin; i += projectionBins)
      {
        // Project the histogram along the Y-axis
        int lastBin = std::min(i + projectionBins - 1, nBinsX);
        TH1D *yProjection = hist2D->ProjectionY(Form("proj_%d", i), i, lastBin);
        // Check if the projection has enough entries to perform a fit
        if (yProjection->GetEntries() < 300)
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
        // float leftmost_limit = 0;
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
        for (int i = 0; i < gausFit->GetNpar(); i++)
        {
          double param = gausFit->GetParameter(i);
          if (std::isnan(param) || std::isinf(param))
          {
            fitFailed = true;
            break;
          }
        }
        if (fitFailed)
        {
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

        // pT_Bins.push_back(pion_pt);
        // pT_Bins_Errors.push_back(0);

        pionmeanGraph[filecounter]->SetPoint(bincounter, pion_pt, Pmean);
        pionmeanGraph[filecounter]->SetPointError(bincounter, 0, PmeanErr);
        pionwidthGraph[filecounter]->SetPoint(bincounter, pion_pt, PWidth);
        pionwidthGraph[filecounter]->SetPointError(bincounter, 0, PWidthErr);
        PresolutionGraph[filecounter]->SetPoint(bincounter, pion_pt, PWidth);
        PresolutionGraph[filecounter]->SetPointError(bincounter, 0, PWidthErr);
        bincounter++;
      }

      MarkerStyle += 1;
      MarkerColor += 1;
      if (MarkerColor == 5 || MarkerColor == 10)
        MarkerColor += 1; // avoid yellow
      if (MarkerStyle == 26 || MarkerStyle == 32)
        MarkerStyle += 1; // avoid triangles
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
      TF1 *PresolutionFit = new TF1("PresolutionFit", "sqrt(2)*sqrt([0]*[0]/x + [1]*[1])", 2, 10);
      PresolutionFit->SetParameters(0.154, 0.02); // Initial guesses for a, b
      // set parameter limits
      // PresolutionFit->SetParLimits(0, 0.1, 0.18);
      PresolutionFit->SetParLimits(1, 0.07, 0.1);
      PresolutionFit->SetLineColor(MarkerColor);
      // Fit the resolution graph
      PresolutionGraph[filecounter]->Fit(PresolutionFit, "RE"); // Fit and constrain to the range of pT
                                                                /*
                                                                      // Create a canvas to plot the resolution graph and fit
                                                                      TCanvas *PresCanvas = new TCanvas("resCanvas", "Resolution Fit", 800, 600);
                                                                      PresolutionGraph[filecounter]->SetTitle("Energy Resolution; p_{T} (GeV/c); Pion #sigma / #mu");
                                                                      PresolutionGraph[filecounter]->Draw("APE");
                                                                      PresolutionFit->Draw("same");
                                                                      // Print the fit parameters on a new canvas
                                                                      TCanvas *PfitParamsCanvas = new TCanvas("fitParamsCanvas", "Fit Parameters", 800, 600);
                                                                      TPaveText *PparamsText = new TPaveText(0.1, 0.7, 0.9, 0.9, "NDC");
                                                                      PparamsText->AddText("Pion Fitted Resolution Parameters:");
                                                                      PparamsText->AddText(FastMC_legendNames[j].c_str());
                                                                      PparamsText->AddText(Form("A/#sqrt{E} term (a): %.4f", PresolutionFit->GetParameter(0)));
                                                                      // paramsText->AddText(Form("Noise term (b): %.4f", PresolutionFit->GetParameter(2)));
                                                                      PparamsText->AddText(Form("Constant term (c): %.4f", PresolutionFit->GetParameter(1)));
                                                                      // add goodness of fit
                                                                      PparamsText->AddText(Form("Chi2/ndf: %.4f", PresolutionFit->GetChisquare() / PresolutionFit->GetNDF()));
                                                                      PparamsText->Draw();
                                                                */
      gPResolutions->Add(PresolutionGraph[filecounter], "PE");
      legend6->AddEntry(PresolutionGraph[filecounter], FastMC_legendNames[j].c_str(), "P");

      // Save the plot to the PDF
      // PresCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
      // PresCanvas->Close();
      // PfitParamsCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
      // PfitParamsCanvas->Close();
    }
    else if (FastMC_FileTypes[j] == 1)
    { // eta
      for (int i = startBin; i <= endBin; i += projectionBins)
      {
        // Project the histogram along the Y-axis
        int lastBin = std::min(i + projectionBins - 1, nBinsX);
        TH1D *yProjection = hist2D->ProjectionY(Form("proj_%d", i), i, lastBin);
        ///*
        // Check if the projection has enough entries to perform a fit
        if (yProjection->GetEntries() < 300)
        { // Adjust the threshold as needed
          delete yProjection;
          continue;
        }
        //*/
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
        std::cout << "Rebinning done" << std::endl;

        histF->Scale(1. / 2, "width");
        std::cout << "Scaling done" << std::endl;
        // Determine the leftmost point with a value in the projection histograms
        // float leftmost_limit = 0;
        if (dynamic_left)
        {
          std::cout << "Dynamic left" << std::endl;
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
        // std::cout << "Leftmost limit: " << limits[0] << std::endl;
        double pt_min = hist2D->GetXaxis()->GetBinLowEdge(i);
        double pt_max = hist2D->GetXaxis()->GetBinUpEdge(lastBin);
        TString ptRange = Form("pt_%.2f-%.2f_GeV", pt_min, pt_max);
        double Eta_pt = (pt_min + pt_max) / 2.0;
        scale_histogram_errors(histF, scale_factor);
        std::cout << "Pre fit Setup done" << std::endl;
        // Fit Eta Gaussian in the specified range
        TF1 *EtagausFit = new TF1("EtagausFit", "gaus", 0.45, 0.75); // limits[6], limits[7]
        EtagausFit->SetParLimits(1, 0.50, 0.64);
        EtagausFit->SetParLimits(2, 0.03, 0.25);
        EtagausFit->SetNpx(1000);
        histF->Fit(EtagausFit, "REQ");

        std::cout << "Fit done" << std::endl;

        // Check if the fit returns NaN or Inf
        bool fitFailed = false;
        for (int i = 0; i < EtagausFit->GetNpar(); i++)
        {
          double param = EtagausFit->GetParameter(i);
          if (std::isnan(param) || std::isinf(param))
          {
            fitFailed = true;
            break;
          }
        }
        if (fitFailed)
        {
          std::cout << "Fit returned NaN or Inf for slice: " << i << std::endl;
          continue;
        }
        std::cout << "Fit parameters are good" << std::endl;

        double Emean = EtagausFit->GetParameter(1);
        double Esigma = EtagausFit->GetParameter(2);
        double EmeanErr = EtagausFit->GetParError(1);
        double EsigmaErr = EtagausFit->GetParError(2);
        double EWidth = Esigma / Emean;
        double EWidthErr = EWidth * sqrt(pow(EmeanErr / Emean, 2) + pow(EsigmaErr / Esigma, 2));

        // double MassRatio = Pion_Mean[j] / Emean;
        // double MassRatioErr = MassRatio * sqrt(pow(Pion_Mean_errors[j] / Pion_Mean[j], 2) + pow(EmeanErr / Emean, 2));
        std::cout << "Fit parameters stored" << std::endl;
        // Eta_Mean.push_back(Emean);
        // Eta_Width.push_back(EWidth);
        // Eta_Mean_errors.push_back(EmeanErr);
        // Eta_Width_errors.push_back(EWidthErr);
        // pT_Bins.push_back(Eta_pt);
        // pT_Bins_Errors.push_back(0);

        etameanGraph[filecounter]->SetPoint(bincounter, Eta_pt, Emean);
        etameanGraph[filecounter]->SetPointError(bincounter, 0, EmeanErr);
        etawidthGraph[filecounter]->SetPoint(bincounter, Eta_pt, EWidth);
        etawidthGraph[filecounter]->SetPointError(bincounter, 0, EWidthErr);
        // massRatioGraph[filecounter]->SetPoint(bincounter, Eta_pt, MassRatio);
        // massRatioGraph[filecounter]->SetPointError(bincounter, 0, MassRatioErr);
        //  EresolutionGraph[filecounter]->SetPoint(bincounter, Eta_pt, EWidth);
        //  EresolutionGraph[filecounter]->SetPointError(bincounter, 0, EWidthErr);

        bincounter++;
        // std::cout << "Eta_pt: " << Eta_pt << " Emean: " << Emean << " EWidth: " << EWidth << std::endl;
        // std::cout << "Bincounter: " << bincounter << std::endl;
        // which file is done
      }

      MarkerStyle += 1;
      MarkerColor += 1;
      if (MarkerColor == 5 || MarkerColor == 10)
        MarkerColor += 1; // avoid yellow
      if (MarkerStyle == 26 || MarkerStyle == 32)
        MarkerStyle += 1; // avoid triangles
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

      // EresolutionGraph[filecounter]->SetMarkerStyle(MarkerStyle);
      // EresolutionGraph[filecounter]->SetMarkerColor(MarkerColor);

      gEtaMeans->Add(etameanGraph[filecounter], "PE");
      legend3->AddEntry(etameanGraph[filecounter], FastMC_legendNames[j].c_str(), "P");

      gEtaWidths->Add(etawidthGraph[filecounter], "PE");
      legend4->AddEntry(etawidthGraph[filecounter], FastMC_legendNames[j].c_str(), "P");

      // gMassRatios->Add(massRatioGraph[filecounter], "PE");
      // legend5->AddEntry(massRatioGraph[filecounter], FastMC_legendNames[j].c_str(), "P");

      //------------------------------------------------------------------------------------------------
      /*
            // Define a function for the eta energy resolution fit
            TF1 *EresolutionFit = new TF1("EresolutionFit", "sqrt(2)*sqrt([0]*[0]/x + [1]*[1])", 2, 20);
            EresolutionFit->SetParameters(0.154, 0.02); // Initial guesses for a, b
            EresolutionFit->SetParLimits(0, 0.1, 0.18);
            EresolutionFit->SetParLimits(1, 0.01, 0.2);
            EresolutionFit->SetLineColor(MarkerColor);
            EresolutionGraph[filecounter]->Fit(EresolutionFit, "RE");


                  TCanvas *EresCanvas = new TCanvas("EresCanvas", "Resolution Fit", 800, 600);
                  EresolutionGraph[filecounter]->SetTitle("Energy Resolution; p_{T} (GeV/c); #sigma / #mu");
                  EresolutionGraph[filecounter]->Draw("APE");
                  EresolutionFit->Draw("same");

                  TCanvas *EfitParamsCanvas = new TCanvas("EfitParamsCanvas", "Fit Parameters", 800, 600);
                  TPaveText *EparamsText = new TPaveText(0.1, 0.7, 0.9, 0.9, "NDC");
                  EparamsText->AddText("Eta Fitted Resolution Parameters:");
                  EparamsText->AddText(FastMC_legendNames[j].c_str());
                  EparamsText->AddText(Form("A/#sqrt{E} term (a): %.4f", EresolutionFit->GetParameter(0)));
                  //paramsText->AddText(Form("Noise term (b): %.4f", EresolutionFit->GetParameter(2)));
                  EparamsText->AddText(Form("Constant term (c): %.4f", EresolutionFit->GetParameter(1)));
                  //add goodness of fit
                  EparamsText->AddText(Form("Chi2/ndf: %.4f", EresolutionFit->GetChisquare() / EresolutionFit->GetNDF()));
                  EparamsText->Draw();

            gEResolutions->Add(EresolutionGraph[filecounter], "PE");
            legend7->AddEntry(EresolutionGraph[filecounter], FastMC_legendNames[j].c_str(), "P");
      */
      // EresCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
      // EresCanvas->Close();
      // EfitParamsCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");
      // EfitParamsCanvas->Close();
    }

    std::cout << "Finished processing file: " << FastMC_FileNames[j] << " done" << std::endl;
    file.Close();
    filecounter++;
  }

  std::cout << "Run2024 files: " << Run2024_FileNames.size() << std::endl;
  // data, pre made graphs
  if (Run2024_FileNames.size() > 0)
  {
    for (size_t j = 0; j < Run2024_FileNames.size(); ++j)
    {
      // Load the new file and retrieve the TGraphErrors
      TFile newFile(Run2024_FileNames[0].c_str(), "READ");
      if (!newFile.IsOpen())
      {
        std::cerr << "Error opening file: " << Run2024_FileNames[0] << std::endl;
        return;
      }

      TGraphErrors *newPionMean = dynamic_cast<TGraphErrors *>(newFile.Get("gr_mass_pi0"));
      TGraphErrors *newPionSigma = dynamic_cast<TGraphErrors *>(newFile.Get("gr_width_pi0"));
      TGraphErrors *newEtaMean = dynamic_cast<TGraphErrors *>(newFile.Get("gr_mass_eta"));
      TGraphErrors *newEtaSigma = dynamic_cast<TGraphErrors *>(newFile.Get("gr_width_eta"));
      TGraphErrors *pionRelativeWidthGraph = new TGraphErrors();
      TGraphErrors *etaRelativeWidthGraph = new TGraphErrors();
      if (!newPionMean || !newPionSigma || !newEtaMean || !newEtaSigma)
      {
        std::cerr << "Error getting TGraphErrors from new file." << std::endl;
        newFile.Close();
        return;
      }

      int nPionPoints = newPionMean->GetN();
      for (int i = 0; i < nPionPoints; ++i)
      {
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
      for (int i = 0; i < nEtaPoints; ++i)
      {
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

      MarkerStyle += 1;
      MarkerColor += 1;
      if (MarkerColor == 5 || MarkerColor == 10)
        MarkerColor += 1; // avoid yellow
      if (MarkerStyle == 26 || MarkerStyle == 32)
        MarkerStyle += 1; // avoid triangles
      // with all off will default to sphenix style
      // newPionMean->SetMarkerStyle(MarkerStyle);
      // newPionMean->SetMarkerColor(MarkerColor);
      // newEtaMean->SetMarkerStyle(MarkerStyle);
      // newEtaMean->SetMarkerColor(MarkerColor);
      // pionRelativeWidthGraph->SetMarkerStyle(MarkerStyle);
      // pionRelativeWidthGraph->SetMarkerColor(MarkerColor);
      // etaRelativeWidthGraph->SetMarkerStyle(MarkerStyle);
      // etaRelativeWidthGraph->SetMarkerColor(MarkerColor);
      //  Add the new graphs to the multigraphs
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
  // gPad->SetFillColor(33);
  gPionMeans->SetTitle("Pion: Smeared pT vs Inv. Mass;#it{pT}_{#gamma#gamma} (GeV); Pion Peak Position (GeV)");
  gPionMeans->GetXaxis()->SetLimits(0.01, 10);
  gPionMeans->SetMinimum(0.130);
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
  // gPad->SetFillColor(33);
  gPionWidths->SetTitle("Pion: Smeared pT vs Relative Width;#it{pT}_{#gamma#gamma} (GeV); Pion Relative Width");
  gPionWidths->GetXaxis()->SetLimits(0.01, 10);
  gPionWidths->SetMinimum(0.05);
  gPionWidths->SetMaximum(0.25);
  gPionWidths->Draw("APE");
  legend2->SetFillStyle(0);
  legend2->SetTextAlign(32);
  legend2->SetTextSize(0.02);
  legend2->Draw();
  c2->Print("pioncode/canvas_pdf/ptdifferentialcomparison.pdf");

  TCanvas *c3 = new TCanvas("c3", "Canvas3", 800, 600);
  // gPad->SetFillColor(33);
  gEtaMeans->SetTitle("Eta: Smeared pT vs Inv. Mass;#it{pT}_{#gamma#gamma} (GeV); Eta Peak Position (GeV)");
  gEtaMeans->GetXaxis()->SetLimits(0.01, 25); // 20
  gEtaMeans->SetMinimum(0.45);
  gEtaMeans->SetMaximum(0.7);
  gEtaMeans->Draw("APE");
  legend3->SetFillStyle(0);
  legend3->SetTextAlign(32);
  legend3->SetTextSize(0.02);
  legend3->Draw();
  c3->Print("pioncode/canvas_pdf/ptdifferentialcomparison.pdf");

  TCanvas *c4 = new TCanvas("c4", "Canvas4", 800, 600);
  // gPad->SetFillColor(33);
  gEtaWidths->SetTitle("Eta: Smeared pT vs Relative Width;#it{pT}_{#gamma#gamma} (GeV); Eta Relative Width");
  gEtaWidths->GetXaxis()->SetLimits(0.01, 25); // 20
  // gEtaWidths->SetMinimum(0.01);
  // gEtaWidths->SetMaximum(0.3);
  gEtaWidths->SetMinimum(0.01);
  gEtaWidths->SetMaximum(0.25); // 11, 25
  gEtaWidths->Draw("APE");
  legend4->SetFillStyle(0);
  legend4->SetTextAlign(32);
  legend4->SetTextSize(0.02);
  legend4->Draw();
  c4->Print("pioncode/canvas_pdf/ptdifferentialcomparison.pdf");

  TCanvas *c5 = new TCanvas("c5", "Canvas5", 800, 600);
  // gPad->SetFillColor(33);
  gMassRatios->SetTitle("Mass Ratios: Smeared pT vs Inv. Mass;#it{pT}_{#gamma#gamma} (GeV); Pion/Eta Mass Ratio");
  gMassRatios->GetXaxis()->SetLimits(0.01, 20);
  gMassRatios->SetMinimum(0.24);
  gMassRatios->Draw("APE");
  legend5->SetFillStyle(0);
  legend5->SetTextAlign(32);
  legend5->SetTextSize(0.02);
  legend5->Draw();
  // c5->Print("pioncode/canvas_pdf/ptdifferentialcomparison.pdf");

  TCanvas *c6 = new TCanvas("c6", "Canvas6", 800, 600);
  // gPad->SetFillColor(33);
  gPResolutions->SetTitle("Pion: Smeared pT vs Resolution;#it{pT}_{#gamma#gamma} (GeV); Pion Resolution");
  gPResolutions->GetXaxis()->SetLimits(0.01, 10);
  gPResolutions->SetMinimum(0.06);
  gPResolutions->SetMaximum(0.16);
  gPResolutions->Draw("APE");
  legend6->SetFillStyle(0);
  legend6->SetTextAlign(32);
  legend6->SetTextSize(0.02);
  legend6->Draw();
  c6->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");

  TCanvas *c7 = new TCanvas("c7", "Canvas7", 800, 600);
  // gPad->SetFillColor(33);
  gEResolutions->SetTitle("Eta: Smeared pT vs Resolution;#it{pT}_{#gamma#gamma} (GeV); Eta Resolution");
  gEResolutions->GetXaxis()->SetLimits(0.01, 25);
  gEResolutions->SetMinimum(0.04);
  // gEResolutions->SetMaximum(0.3);
  gEResolutions->Draw("APE");
  legend7->SetFillStyle(0);
  legend7->SetTextAlign(32);
  legend7->SetTextSize(0.02);
  legend7->Draw();
  c7->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf");

  // Close the PDF file
  dummyCanvas->Print("pioncode/canvas_pdf/ptdifferentialcomparison.pdf]");
  dummyCanvas->Print("pioncode/canvas_pdf/ptdifferential_Energyres_results.pdf]");
  dummyCanvas->Print("pioncode/canvas_pdf/ptdifferential_PionFit_results.pdf]");
  dummyCanvas->Print("pioncode/canvas_pdf/ptdifferential_EtaFit_results.pdf]");
  dummyCanvas->Print("pioncode/canvas_pdf/ptdifferential_unw_Fit_results.pdf]");

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
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_run2pp_ana450_2024p009_000merged_V2.root",
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_run2pp_ana450_2024p009_000merged_V5.root",//without extra cuts
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_run2pp_ana450_2024p009_000merged_V5.root",
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_run2pp_ana450_2024p009_000merged_V10.root",//with extra cuts
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_run2pp_ana450_2024p009_000merged_V10.root",
      "pioncode/rootfiles/OUTHIST_iter_DST_CALOFITTING_run2pp_ana509_2024p022_v001_0004merged_V14.root ",
      "pioncode/rootfiles/OUTHIST_iter_DST_CALOFITTING_run2pp_ana509_2024p022_v001_0004merged_V14.root",
      "pioncode/rootfiles/OUTHIST_iter_DST_CALOFITTING_run2pp_ana509_2024p022_v001_0004merged_V14.root",
      "pioncode/rootfiles/OUTHIST_iter_DST_CALOFITTING_run2pp_ana509_2024p022_v001_0004merged_V14.root",
      "pioncode/rootfiles/OUTHIST_iter_DST_CALOFITTING_run2pp_ana509_2024p022_v001_0004merged_V14.root",
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_pythia8_Detroit_0000000021_0merged_V85.root",// calib
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALOFITTING_run2pp_ana509_2024p022_v001_0004merged_V14.root",// wfm
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_pythia8_Jet30_0000000021_00merged_V78.root",
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_pythia8_70mev_0000000021_00merged_V79.root",
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_run2pp_ana450_2024p009_000merged_V6.root",
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_run2pp_ana450_2024p009_000merged_V6.root",
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_pythia8_pp_mb_0000000015_merged_V65.root",
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_pythia8_pp_mb_0000000015_merged_V66.root",
      ////"pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_pythia8_pp_mb_0000000015_merged_V67.root",
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_pythia8_pp_mb_0000000015_merged_V69.root",
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_pythia8_pp_mb_0000000015_merged_V70.root",
      ////"pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_pythia8_pp_mb_0000000015_merged_V71.root",
      ////"pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_pythia8_pp_mb_0000000015_merged_V72.root",
      ////"pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_pythia8_pp_mb_0000000015_merged_V73.root",
  }; //    "pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_pythia8_pp_mb_3MHz_0000000011__merged_V1.root",
  //
  std::vector<std::string> unweighted_histNames = {
      //"h_diPhotonMasspT",
      //"h_diPhotonMasspT_photon5",
      //"h_diPhotonMasspT_photon4",
      //"h_diPhotonMasspT_photon5_mbd",
      //"h_diPhotonMasspT_photon5",
      "h_diPhotonMasspT_photon4",
      "h_diPhotonMasspT_photon5",
      "h_diPhotonMasspT_jet12",
      "h_diPhotonMasspT_jet10",
      "h_diPhotonMasspT_jet8",
  }; //"h_InvMass_2d",
  std::vector<std::string> unweighted_legendNames = {
      //"run2024_12/21/24",
      //"run2024_1/9/25_p5",
      //"run2024_1/9/25_p4",
      //"ana450_1/9/25_p5+mbd",
      //"ana450_1/9/25_p5",
      "ana509_9/19/25_p4",
      "ana509_9/19/25_p5",
      "ana509_9/19/25_j12",
      "ana509_9/19/25_j10",
      "ana509_9/19/25_j8",
      //"run2024_1/9/25_p5+mbd_eta/zvtx_cuts",
      //"run2024_1/9/25_p4+mbd_eta/zvtx_cuts",
      //"Detroit_MB_calib+0%",
      //"Detroit_MB_wfm+0%",
      //"Pythia_10GeVJets",
      //"Pythia_30GeVJets",
      //"eT_prob","prob_eT",
      //"70mev_Pythia_wvfm_vtx+11.5%smr",
      //"Pythia_wvfm_vtx+0smr",
      //"70mev_Pythia_wvfm_vtx+0smr"
  }; //"Pythia",

  //-----------------------------------------
  std::vector<std::string> SPMC_FileNames = {
      ////"pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_single_pi0_p_200_20000MeV_0000000017_00merged_V43.root",
      ////"pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_single_eta_p_600_20000MeV_0000000017_00merged_V44.root",
      ////"pioncode/rootfiles/OUTHIST_iter_G4Hits_single_eta_p_600_20000MeV_0000000017_00merged_V14.root",
      //"pioncode/rootfiles/OUTHIST_iter_G4Hits_single_eta_p_600_20000MeV_0000000017_00merged_V15.root",
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_pt_200_40000MeV_0000000024_0merged_V84.root",// new calib sample 2/26/25
      "pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_pt_200_40000MeV_0000000024_0merged_V7.root", // new wfm sample
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_pt_200_40000MeV_0000000024_0merged_V3.root",
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_pt_200_40000MeV_0000000024_0merged_V4.root",
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_pt_200_40000MeV_0000000024_0merged_V4.root",
      //"pioncode/rootfiles/OUTHIST_iter_G4Hits_single_eta_p_600_20000MeV_0000000017_00merged_V21.root",
      //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_single_pi0_p_200_20000MeV_0000000017_00merged_V42.root",
      "pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_single_eta_pt_200_50000MeV_0000000024_0merged_V8.root", // new wfm sample
                                                                                                                //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_single_eta_pt_200_50000MeV_0000000024_00merged_V2.root",
                                                                                                                //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_single_eta_pt_200_50000MeV_0000000024_00merged_V6.root",
                                                                                                                //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_single_eta_pt_200_50000MeV_0000000024_00merged_V6.root",
                                                                                                                //"pioncode/rootfiles/OUTHIST_iter_G4Hits_single_pi0_p_200_20000MeV_0000000017_00merged_V6.root",
                                                                                                                //"pioncode/rootfiles/OUTHIST_iter_G4Hits_single_pi0_p_200_20000MeV_0000000017_00merged_V7.root",
  };
  //
  //"pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_single_pi0_p_200_20000MeV_0000000017_00merged_V38
  //,"pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_200_10000MeV_0000000013_00merged_V14.root"
  std::vector<std::string> SPMC_histNames = {
      //"h_InvMass_smear_weighted_2d_0",
      //"h_InvMass_smear_weighted_2d_0",
      //"h_InvMass_smear_weighted_2d_125",
      "h_InvMass_smear_weighted_2d_0",
      //"h_truthmatched_mass_weighted_2d",
      //"h_InvMass_smear_weighted_2d_125",
      //"h_truthmatched_mass_weighted_2d",
      //"h_truthmatched_mass_etameson_weighted_2d",
      //"h_InvMass_smear_weighted_2d_125",
      "h_InvMass_smear_weighted_2d_0",
      "h_truthmatched_mass_etameson_weighted_2d",
      "h_InvMass_smear_weighted_2d_125",
      "h_truthmatched_mass_etameson_weighted_2d",
  };
  std::vector<std::string> SPMC_legend = {
      //"SPi0+0sm", "SEta+0sm",
      //"SEta+12.5sm,70mev,lowcut",
      //"SEta+12.5sm,30mev,lowcut",
      //"SEta+12.5sm,30mev,lowcut",
      //"SEta+0sm,70mev,lowcut",
      //"calib, SPi0+0%",
      "wfm, SPi0+0%",
      //"New SPi0+12.5%",
      //"New SPi0+12.5%+match",
      //"SEta+0sm+match,70mev",
      //"SPi0+12.5sm, more_cuts",
      //"SPi0_Weight_pythia+12.5",
      //"SPi0_tbtzs_Weight_pythia+12.5"
      "wfm, SEta+0%%",
      "New SEta+0%+match",
      "New SEta+12.5%",
      "New SEta+12.5%+match"};

  std::vector<int> SPMC_FileTypes = {
      // 0,
      0,
      // 0,
      // 0,
      // 0,
      1,
      // 0,
      1,
      1,
      1,
      1}; // 0 for pion, 1 for eta

  //-----------------------------------------
  std::vector<std::string> FastMC_fileNames =
      //*
      {
          /*
            "pioncode/rootfiles/PionFastMC_0.150000_sqrte_0.120000_const.root",
            "pioncode/rootfiles/PionFastMC_0.150000_sqrte_0.120000_const.root",
            ////"pioncode/rootfiles/EtaFastMC_0.100000_sqrte_0.170000_const.root",
            //"pioncode/rootfiles/EtaFastMC_0.158000_sqrte_0.027000_const.root",
            ////"pioncode/rootfiles/EtaFastMC_0.110000_sqrte_0.060000_const.root",
            ////"pioncode/rootfiles/EtaFastMC_0.149000_sqrte_0.032000_const.root",
            "pioncode/rootfiles/EtaFastMC_0.149000_sqrte_0.020000_const.root",
            "pioncode/rootfiles/EtaFastMC_0.149000_sqrte_0.030000_const.root",
            "pioncode/rootfiles/EtaFastMC_0.149000_sqrte_0.050000_const.root",
            "pioncode/rootfiles/EtaFastMC_0.130000_sqrte_0.080000_const.root",
            ////"pioncode/rootfiles/PionFastMC_0.140000_sqrte_0.004000_const.root",
            ////"pioncode/rootfiles/PionFastMC_0.140000_sqrte_0.004000_const.root",
            ////"pioncode/rootfiles/PionFastMC_0.140000_sqrte_0.004000_const.root",
            ////"pioncode/rootfiles/PionFastMC_0.140000_sqrte_0.004000_const.root",
            */
      }; //*/
  std::vector<std::string> FastMC_legendNames =
      {
          "FastMC: 15%/#sqrt{E} #oplus 12% no tct",
          "FastMC: 15%/#sqrt{E} #oplus 12% + tct",
          //"FastMC: 15.8%/#sqrt{E} #oplus 2.7%",
          //"FastMC: 11.0%/#sqrt{E} #oplus 6.0%",
          //"FastMC: 14.9%/#sqrt{E} #oplus 3.2%",
          "FastMC: 14.9%/#sqrt{E} #oplus 2.0%",
          "FastMC: 14.9%/#sqrt{E} #oplus 3.0%",
          "FastMC: 14.9%/#sqrt{E} #oplus 5.0%",
          "FastMC: 13.0%/#sqrt{E} #oplus 8.0%",
          ////"FastMC_symm_E: 14%/#sqrt{E} #oplus 4%",
          ////"FastMC_asymm_E: 14%/#sqrt{E} #oplus 4%",
          ////"FastMC_symm_R: 14%/#sqrt{E} #oplus 4%",
          ////"FastMC_symm_T: 14%/#sqrt{E} #oplus 4%"
      }; //"PionFastMC", "EtaFastMC"
  std::vector<std::string> FastMC_histNames =
      {
          "h100_2",
          "h102_2",
          "h100_2",
          "h100_2",
          "h100_2",
          "h100_2",
          "h100_2",
          "h100_2",
          "h100_2",
          //"h101_2_symm_2",
          //"h101_2_asymm_2",
          //"h101_2_symm_0", "h101_2_symm_1"
      };
  std::vector<int> FastMC_FileTypes =
      {
          0,
          0,
          // 1,
          1,
          // 1,
          1,
          1,
          1,
          1,
          // 0,
          // 0,
          // 0, 0
      }; // 0 for pion, 1 for eta
  //
  //-----------------------------------------
  std::vector<std::string> Run2024_fileNames =
      {
          //"pioncode/rootfiles/meson_graphs.root"
      };
  std::vector<std::string> Run2024_legendNames = {"Run2024"};

  //-----------------------------------------
  AnalyzeHistograms(unweighted_fileNames, unweighted_histNames, unweighted_legendNames, SPMC_FileNames, SPMC_histNames, SPMC_legend, SPMC_FileTypes, FastMC_fileNames, FastMC_histNames, FastMC_legendNames, FastMC_FileTypes, Run2024_fileNames, Run2024_legendNames);
  gApplication->Terminate(0);
  // return 0;
}
