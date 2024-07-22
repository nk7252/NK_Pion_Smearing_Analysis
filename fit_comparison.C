#include <iostream>
#include <vector>
#include <string>
#include <TFile.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraphErrors.h>

void AnalyzeHistograms(const std::vector<std::string>& GeantFileNames,const std::vector<std::string>& FastMCFileNames, const std::vector<std::string>& GeanthistNames, const std::vector<std::string>& FastMChistNames) {
    std::vector<double> Pion_Means, Pion_Sigmas, Eta_Means, Eta_Sigmas, Mass_Ratios;
    std::vector<double> errorsMean, errorsSigma, errorsSecondMean, errorsSecondSigma, errorsRatio;
    std::vector<double> xValues, xErrors;

    for (size_t i = 0; i < fileNames.size(); ++i) {
        TFile file(fileNames[i].c_str(), "READ");
        if (!file.IsOpen()) {
            std::cerr << "Error opening file: " << fileNames[i] << std::endl;
            continue;
        }

        TH2* hist = dynamic_cast<TH2*>(file.Get(histNames[i].c_str()));
        if (!hist) {
            std::cerr << "Error getting histogram: " << histNames[i] << " from file: " << fileNames[i] << std::endl;
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

        xValues.push_back(i);
        xErrors.push_back(0);

        file.Close();
    }

    // Create TGraphErrors
    TGraphErrors* graphMean = new TGraphErrors(means.size(), xValues.data(), means.data(), xErrors.data(), errorsMean.data());
    TGraphErrors* graphSigma = new TGraphErrors(sigmas.size(), xValues.data(), sigmas.data(), xErrors.data(), errorsSigma.data());
    TGraphErrors* graphSecondMean = new TGraphErrors(secondMeans.size(), xValues.data(), secondMeans.data(), xErrors.data(), errorsSecondMean.data());
    TGraphErrors* graphSecondSigma = new TGraphErrors(secondSigmas.size(), xValues.data(), secondSigmas.data(), xErrors.data(), errorsSecondSigma.data());
    TGraphErrors* graphRatio = new TGraphErrors(ratios.size(), xValues.data(), ratios.data(), xErrors.data(), errorsRatio.data());

    // Optionally save graphs to a file
    TFile outputFile("fitResults.root", "RECREATE");
    graphMean->Write("MeanGraph");
    graphSigma->Write("SigmaGraph");
    graphSecondMean->Write("SecondMeanGraph");
    graphSecondSigma->Write("SecondSigmaGraph");
    graphRatio->Write("RatioGraph");
    outputFile.Close();
}

int main() {
    std::vector<std::string> Geant_fileNames = {"pioncode/OUTHIST_iter_DST_CALO_CLUSTER_pythia8_pp_mb_3MHz_0000000011_merged.root", "pioncode/rootfiles/OUTHIST_iter_DST_CALO_CLUSTER_single_pi0_200_10000MeV_0000000013_sm125_nopc.root"};
    std::vector<std::string> fastmc_fileNames = {"pioncode/rootfiles/PionFastMC_0.154000_sqrte_0.077000_const.root", "pioncode/rootfiles/EtaFastMC_0.154000_sqrte_0.077000_const.root"};
    std::vector<std::string> Geant_histNames = {"h_InvMass_2d", "h_InvMass_2d"};
    std::vector<std::string> fastmc_histNames = {"h100_2", "h100_2"};

    AnalyzeHistograms(Geant_fileNames,fastmc_fileNames, Geant_histNames,fastmc_histNames);

    return 0;
}
