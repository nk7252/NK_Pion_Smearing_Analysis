#include <iostream>
#include <cstring>
#include <string>
#include <fstream>

#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TFitResult.h>

void Photon_Fit_comparison()
{
    //open file
    TFile *file = new TFile("pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_pythia8_pp_mb_0000000015_merged_V54.root", "READ");
    if (!file->IsOpen())
    {
        std::cerr << "Error opening file: " << "pioncode/rootfiles/OUTHIST_iter_DST_CALO_WAVEFORM_pythia8_pp_mb_0000000015_merged_V54.root" << std::endl;
        return;
    }

    TH2F *hist2D = (TH2F *)gDirectory->Get("h_truth_spectrum2");

    if (!hist2D)
    {
        std::cerr << "Error retrieving histogram: " << "h_reco_ALLphotonE_2d_2" << " from file" << std::endl;
        file->Close();
        return;
    }
    //fit slices in y

	hist2D->FitSlicesY(0, 0, -1, 0);
    TCanvas *canvas = new TCanvas("canvas", "canvas", 800, 800);
    canvas->Divide(2, 2);
    //make mu plots
    canvas->cd(1);
    TH2F *hMu = (TH2F *)gDirectory->Get("h_reco_ALLphotonE_2d_1");
    hMu->SetTitle("#mu");
    hMu->GetXaxis()->SetTitle("E_{#gamma} [GeV]");
    hMu->GetYaxis()->SetTitle("Bin Number");
    hMu->Draw("colz");
    //make sigma plots
    canvas->cd(2);
    TH2F *hSigma = (TH2F *)gDirectory->Get("h_reco_ALLphotonE_2d_2");
    hSigma->SetTitle("#sigma");
    hSigma->GetXaxis()->SetTitle("E_{#gamma} [GeV]");
    hSigma->GetYaxis()->SetTitle("Bin Number");
    hSigma->Draw("colz");
    //make sigma/mu plots
    canvas->cd(3);
    TH2F *hSigmaOverMu = (TH2F *)hMu->Clone("hSigmaOverMu"); // clone mean
    hSigmaOverMu->Divide(hSigma);
    hSigmaOverMu->SetTitle("#sigma/#mu");
    hSigmaOverMu->GetXaxis()->SetTitle("");
    hSigmaOverMu->GetYaxis()->SetTitle("Bin Number");
    hSigmaOverMu->Draw("colz");
}
