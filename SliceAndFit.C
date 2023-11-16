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
// #include <boost/filesystem.hpp>
// #include <boost/multiprecision/cpp_dec_float.hpp>
//-------------------------------
//function declarations;
//-------------------------------
void FitYProjectionsAndGraph(TCanvas *canvas, TH2F *hist2D, const char *pdfname, double binres, int weightmethod);

void SliceAndFit()
{
	// Change some default parameters in the current style
	// gStyle->SetLabelSize(0.06, "x");
	// gStyle->SetLabelSize(0.06, "y");
	// gStyle->SetFrameFillColor(38);
	// gStyle->SetTitleW(0.6);
	// gStyle->SetTitleH(0.1);

	// Connect the input file and get the 2-d histogram in memory
	/*
	std::unique_ptr<TFile> pionfile(TFile::Open("pioncode/rootfiles/Pi0FastMC_0.067000.root"));
	if (!pionfile || pionfile->IsZombie()) {
		std::cerr << "Error opening file" << endl;
		exit(-1);
	}

	*/
	time_t rawtime;
	struct tm *timeinfo;
	char Time[80]; // time string with size given
	time(&rawtime);
	timeinfo = localtime(&rawtime);
	strftime(Time, 80, "%Y-%m-%d-%H:%M:%S", timeinfo);
	puts(Time);				   // print Time to root
	int rel_error_param = 155; // 65 for full set. I will do 145 now// switched from 67(old) to 65(new)


	std::vector<std::string> WeightNames = {"EXP", "POWER", "WSHP"};
	int weightmethod=2;//0=exp,1=power,2=wshp
	double binres=2;

	for (int i = 0; i < 1; i++)
	{ // 3 if doing 14.5-16.5//26 for full set(old)
		////////////////////////////////////
		if (rel_error_param > 320)
		{
			break;
		}

		// const char* canvasname;
		// cout << rel_error_param;//Printf(rel_error_param);
		const char *pionfilename;
		if (rel_error_param < 100)
		{
			pionfilename = Form("Pi0FastMC_0.0%d000_%s.root", rel_error_param,WeightNames[weightmethod].c_str());
		}
		else
		{
			pionfilename = Form("Pi0FastMC_0.%d000_%s.root", rel_error_param,WeightNames[weightmethod].c_str());
		}
		// const char* pionfilename = Form("Pi0FastMC_%d000.root", rel_error_param);
		cout << pionfilename << "\n";
		TFile *pionfile = new TFile(Form("pioncode/rootfiles/%s", pionfilename)); // TFile *_file0 = TFile::Open("A:/root_v6.26.06/bin/pioncode/rootfiles/Pi0FastMC_0.067000.root")  then .ls to see files

		// TFile* pionfile = new TFile(Form("pioncode/rootfiles/Pi0FastMC_0.067000.root", rel_error_param));//TFile *_file0 = TFile::Open("A:/root_v6.26.06/bin/pioncode/rootfiles/Pi0FastMC_0.067000.root")  then .ls to see files
		// TString dir = gROOT->GetTutorialDir();
		// dir.Append("/hsimple.C");
		// dir.ReplaceAll("/./", "/");
		// if (!gInterpreter->IsLoaded(pionfile.Data())) gInterpreter->LoadMacro(pionfile.Data());
		// TFile* hsimpleFile = (TFile*)gROOT->ProcessLineFast("hsimple(1)");
		// if (!hsimpleFile) return;
		TH2F *h9 = (TH2F *)pionfile->Get("h9"); // h18 weighted, h9 unweighted
		TH2F *h18 = (TH2F *)pionfile->Get("h18");
		// Create a canvas and divide it
		// int i = 0;
		// int canvas_param = rel_error_param + i;

		// delete pionfilename;

		// const char* canvasname = Form("weighted_unweighted_SliceFit_ErrParam_67_thousandths", Time);//_date_%s,rel_error_param
		TString canvasname = Form("%s_Sliced_%d_thousandths_%s",WeightNames[weightmethod].c_str(), rel_error_param, Time); //_date_%s
		const char *pdfname = canvasname;
		TCanvas *c1 = new TCanvas(canvasname, canvasname, 3000, 1200);
		// TCanvas* c1 = new TCanvas(Form("c_%d", rel_error_param), "c1", 3000, 1000);
		// TCanvas* c1 = new TCanvas("c1", "c1", 3000, 1000);
		// c1->SetFillColor(42);
		c1->Divide(2, 1); // nx, ny
		TPad *leftPad = (TPad *)c1->cd(1);
		;
		leftPad->Divide(2, 3); // nx, ny

		////////////////////////////////// unweighted
		leftPad->cd(1);
		// c1->cd(1);
		gPad->SetTopMargin(0.12);
		gPad->SetFillColor(33);
		h9->GetXaxis()->SetTitle("Pion Pt [GeV/c]");
		h9->GetYaxis()->SetTitle("Invariant Mass [GeV/c^2]");
		h9->Draw("colz");
		printf("error test code\n");
		h9->GetXaxis()->SetLabelSize(0.06);
		h9->GetYaxis()->SetLabelSize(0.06);
		h9->SetMarkerColor(kYellow);
		// Fit slices projected along Y fron bins in X [1,64] with more than 2 bins  in Y filled
		h9->FitSlicesY(0, 0, -1, 0);

		//

		leftPad->cd(2);
		gPad->SetFillColor(33);
		TH2F *h9_0 = (TH2F *)pionfile->Get("h9_0");
		h9_0->GetXaxis()->SetTitle("Pion Pt [GeV/c]");
		h9_0->Draw();
		// TPad* rightPad = (TPad*)c1->cd(2);
		// rightPad->Divide(1, 2);
		// rightPad->cd(1);
		// Show fitted "mean" for each slice //Double_t mm = h->GetMean();
		leftPad->cd(3);
		gPad->SetTopMargin(0.12);
		gPad->SetLeftMargin(0.15);
		gPad->SetFillColor(33);
		TH2F *h9_1 = (TH2F *)pionfile->Get("h9_1");
		h9_1->GetXaxis()->SetTitle("Pion Pt [GeV/c]");
		h9_1->GetYaxis()->SetTitle("Mean");
		h9_1->Draw();

		// Show fitted "sigma" for each slice
		// rightPad->cd(2);
		leftPad->cd(4);
		gPad->SetTopMargin(0.12);
		gPad->SetLeftMargin(0.15);
		gPad->SetFillColor(33);
		TH2F *h9_2 = (TH2F *)pionfile->Get("h9_2");
		h9_2->SetMinimum(0.8);
		h9_2->GetXaxis()->SetTitle("Pion Pt [GeV/c]");
		h9_2->GetYaxis()->SetTitle("Sigma");
		h9_2->Draw();
		///*
		// Show fitted variance(sigma^2) for each slice
		leftPad->cd(5);
		TH2F *h9_3 = (TH2F *)h9_2->Clone("h9_3"); // clone sigma
		gPad->SetTopMargin(0.12);
		gPad->SetLeftMargin(0.15);
		gPad->SetFillColor(33);
		h9_3->Multiply(h9_2); // multiply cloned sigma by sigma
		h9_3->SetTitle("Value of par[2]^2=Variance;Pion Pt [GeV/c];Variance");
		// h9_3->SetMinimum(0.8);
		h9_3->Draw();

		// Show fitted mean/variance for each slice
		leftPad->cd(6);
		TH2F *h9_4 = (TH2F *)h9_1->Clone("h9_4"); // clone mean
		gPad->SetTopMargin(0.12);
		gPad->SetLeftMargin(0.15);
		gPad->SetFillColor(33);
		h9_4->Divide(h9_2);
		h9_4->SetTitle("Value of par[1]/par[2]^2=Mean/Variance;Pion Pt [GeV/c];Mean/Variance");
		// h9_4->SetMinimum(0.8);
		h9_4->Draw();
		//*/
		/////////////////////////////////////////////weighted
		TPad *rightPad = (TPad *)c1->cd(2);
		rightPad->Divide(2, 3);
		rightPad->cd(1);
		// c1->cd(1);
		gPad->SetTopMargin(0.12);
		gPad->SetFillColor(33);
		gPad->SetLogz();
		h18->Draw("colz");
		printf("error test code\n");
		h18->GetXaxis()->SetLabelSize(0.06);
		h18->GetYaxis()->SetLabelSize(0.06);
		h18->GetXaxis()->SetTitle("Pion Pt [GeV/c]");
		h18->GetYaxis()->SetTitle("Invariant Mass [GeV/c^2]");
		h18->SetMarkerColor(kYellow);
		// Fit slices projected along Y fron bins in X [1,64] with more than 2 bins  in Y filled
		h18->FitSlicesY(0, 0, -1, 0);//, "EMW"

		//

		// Show fitted "mean" for each slice
		rightPad->cd(2);
		gPad->SetFillColor(33);
		TH2F *h18_0 = (TH2F *)pionfile->Get("h18_0");
		h18_0->GetXaxis()->SetTitle("Pion Pt [GeV/c]");
		h18_0->Draw();
		// TPad* rightPad = (TPad*)c1->cd(2);
		// rightPad->Divide(1, 2);
		// rightPad->cd(1);
		rightPad->cd(3);
		gPad->SetTopMargin(0.12);
		gPad->SetLeftMargin(0.15);
		gPad->SetFillColor(33);
		TH2F *h18_1 = (TH2F *)pionfile->Get("h18_1");
		h18_1->GetYaxis()->SetTitle("Mean");
		h18_1->GetXaxis()->SetTitle("Pion Pt [GeV/c]");
		// h18_1->SetAxisRange(0.1, 0.16,"Y");
		h18_1->Draw();

		// Show fitted "sigma" for each slice
		// rightPad->cd(2);
		rightPad->cd(4);
		gPad->SetTopMargin(0.12);
		gPad->SetLeftMargin(0.15);
		gPad->SetFillColor(33);
		TH2F *h18_2 = (TH2F *)pionfile->Get("h18_2");
		h18_2->SetMinimum(0.8);
		h18_2->GetYaxis()->SetTitle("Sigma");
		h18_2->GetXaxis()->SetTitle("Pion Pt [GeV/c]");
		h18_2->Draw();

		// Show fitted variance(sigma^2) for each slice
		rightPad->cd(5);
		TH2F *h18_3 = (TH2F *)h18_2->Clone("h18_3"); // clone sigma
		gPad->SetTopMargin(0.12);
		gPad->SetLeftMargin(0.15);
		gPad->SetFillColor(33);
		h18_3->Multiply(h18_2); // multiply cloned sigma by sigma
		// h9_3->SetMinimum(0.8);
		h18_3->SetTitle("Value of par[2]^2=Variance;Pion Pt [GeV/c];Variance");
		h18_3->Draw();

		// Show fitted mean/variance for each slice
		rightPad->cd(6);
		// TH2F* h18_4 = (TH2F*)h18_1->Clone("h18_4");//clone mean
		TH2F *h18_4 = (TH2F *)h18_2->Clone("h18_4"); // clone sigma instead
		gPad->SetTopMargin(0.12);
		gPad->SetLeftMargin(0.15);
		gPad->SetFillColor(33);
		h18_4->Divide(h18_3); // divide by variance. 2 is sigma
		h18_4->SetTitle("Value of par[1]/par[2]^2=Mean/Variance;Pion Pt [GeV/c];Mean/Variance");
		// h18_4->Divide(h18_1);//divide sigma by mean.
		// h18_4->SetTitle("Value of par[2]/par[1]=Sigma/Mean;Pion Pt [GeV/c];Sigma/Mean");
		// h18_4->SetMinimum(0.8);
		h18_4->Draw();

		// const char* canvasname = Form("weighted_unweighted_SliceFit_ErrParam_%d/1000_date_%s", rel_error_param, Time);
		// const char* canvasname = Form("pioncode/canvas_pdf/collected_canvas_ErrParam_%d/1000_Date_%s", rel_error_param, Time);
		cout << canvasname << "\n";
		// c1->SaveAs(Form("pioncode/canvas_pdf/Weighted_unweighted_SliceFit_ErrParam%s.pdf", canvasname));
		// c1->SaveAs(Form("pioncode/canvas_pdf/%s.pdf", canvasname));
		c1->SaveAs(Form("pioncode/canvas_pdf/%s.pdf", pdfname));
		// c1->SaveAs(Form("pioncode/canvas_pdf/Weighted_unweighted_SliceFit_ErrParam_%d/1000_Date_%s.pdf", rel_error_param, Time));
		delete c1;

		if (rel_error_param == 155)
		{
			TCanvas *c2 = new TCanvas("c2", "c2", 400, 900);
			c2->Divide(1, 2);
			c2->cd(1);
			TH2F *h18_5 = (TH2F *)h18_2->Clone("h18_5");
			h18_5->SetAxisRange(0., 16., "x");
			// h18_5->Scale(1000/135);
			// h18_5->SetAxisRange(0., 0.2,"y");
			h18_5->Draw();
			h18_5->GetYaxis()->SetTitle("Sigma");
			h18_5->GetXaxis()->SetTitle("Pion Pt [GeV/c]");

			c2->cd(2);
			/* they said  "A beam momentum spread (δp/p ≈ 2%) is quadratically subtracted from σ/μ of the fit, in order to unfolded beam momentum spread from the relative energy resolution. The Gauss function parameter of μ and energy resolution from each fit are plotted against the  nominal beam energy as linearity and resolution." //*/
			TH2F *h18_6 = (TH2F *)h18_2->Clone("h18_6"); // sigma/mean

			h18_6->Divide(h18_1); // sigma/mean

			// TH2F* h18_7 = (TH2F*)h18_6->Clone("h18_7");// clone (sigma/mean)to get (sigma/mean)^2
			// h18_7->Multiply(h18_6);//(sigma/mean)^2
			// h18_7->Add(-0.0004);//why?
			// TH2F* h18_8 = (TH2F*)h18_7->Clone("h18_8");
			// h18_8->Divide(h18_7);// call this the new sigma/mean
			// h18_8->Multiply(h18_1);//multiply by the mean to find a new Sigma
			// h18_8->Divide(135);//scale it to find sigma_M/M directly

			h18_6->SetAxisRange(0., 16., "x");
			h18_6->SetAxisRange(0., 0.2, "y");
			h18_6->Draw();
			h18_6->GetYaxis()->SetTitle("Sigma/Mean");
			h18_6->GetXaxis()->SetTitle("Pion Pt [GeV/c]");
			//*/
			c2->SaveAs(Form("pioncode/canvas_pdf/%s_truncatedsigma.pdf", pdfname));
			delete c2;
		}

		
		TCanvas *c3 = new TCanvas("c3", "c3", 3000, 3000);
		
		FitYProjectionsAndGraph(c3, h18, pdfname, binres, weightmethod);
		c3->SaveAs(Form("pioncode/canvas_pdf/Alt_Projection_%s.pdf", pdfname));
		
		pionfile->Close();

		rel_error_param = rel_error_param + 10;
		cout << rel_error_param << "\n";
		delete c3;
	}
	//return 0;
	
}

void FitYProjectionsAndGraph(TCanvas *canvas, TH2F *hist2D, const char *pdfname, double binres, int weightmethod)
{	
	std::ofstream mycsv;
	//const int numXBins = 64; // Number of X bins
	int numXBins = hist2D->GetNbinsX();
	mycsv.open(Form("pioncode/csvfiles/Fit_Param_%s.csv", pdfname));
	// Create TGraphs to store the fit parameters
	//TGraph *ampGraph = new TGraph();
	//TGraph *meanGraph = new TGraph();
	//TGraph *sigmaGraph = new TGraph();
	TGraphErrors *ampGraph = new TGraphErrors(numXBins);
    TGraphErrors *meanGraph = new TGraphErrors(numXBins);
    TGraphErrors *sigmaGraph = new TGraphErrors(numXBins);
	// Loop over X bins
	int upper_limit, lower_limit;
	if (weightmethod==0){//exp
		lower_limit=1;
		upper_limit=binres*6;
	}
	else if (weightmethod==1){//power 
		lower_limit=binres*2;
		upper_limit=numXBins;
	}
	else if (weightmethod==2 ){//woods saxon
		lower_limit=1;
		upper_limit=numXBins;
	}
	for (int i = lower_limit; i <= upper_limit; ++i)// originally until numXBins// for exp use i < 4*6+1// for power go from 4*3+1 to numXBins+1
	{
		// Create a Y projection for the current X bin
		TH1D *yProjection = hist2D->ProjectionY(Form("YProjection_%d", i), i , i ,"");

		// Fit the Y projection with a Gaussian
		yProjection->Fit("gaus", "Q");
		
		// Access the fit parameters
		TF1 *fitFunc = yProjection->GetFunction("gaus");
		if(fitFunc){
			double amplitude = fitFunc->GetParameter(0);
			double mean = fitFunc->GetParameter(1);
			double sigma = fitFunc->GetParameter(2);
			double amplitudeError = fitFunc->GetParError(0);
			double meanError = fitFunc->GetParError(1);
			double sigmaError = fitFunc->GetParError(2);

			// Store the parameters in the TGraphs
			ampGraph->SetPoint(i, i/binres, amplitude);
			meanGraph->SetPoint(i, i/binres, mean);
			sigmaGraph->SetPoint(i, i/binres, sigma);
			ampGraph->SetPointError(i, 0, amplitudeError);// No error in X
			meanGraph->SetPointError(i, 0, meanError); // No error in X
			sigmaGraph->SetPointError(i, 0, sigmaError); // No error in 
			//*
			if(i==0){
			mycsv << "Bin Number" << "," << "," << "# entries" << "," << ","  << "Amplitude" << "," << "Error" << "," << "," << "Sigma" << "," << "Error" << "," << "," << "Mean" << "," << "Error" << "," << "," << "Sigma/Mean" << "\n";
			}
			else{
			mycsv << i << "," << "," << yProjection->GetEntries() << "," << "," << amplitude << ","<< amplitudeError << "," << "," << sigma << ","<< sigmaError<< "," << "," << mean << ","<< meanError << "," << "," << sigma/mean << "\n";
			}
			//*/
		}
        
	}
	mycsv.close();
	// Create a multi-panel canvas and plot the TGraphs
	canvas->Divide(2, 2); // Divide the canvas into 3 parts
	canvas->cd(1);
    ampGraph->SetTitle("Amplitude vs. Bin Number");
    ampGraph->GetXaxis()->SetTitle("Bin Number");
    ampGraph->GetYaxis()->SetTitle("Amplitude");
	//gPad->SetTopMargin(0.12);
    ampGraph->Draw("AP");
	canvas->cd(2);
    meanGraph->SetTitle("Mean vs. Bin Number");
    meanGraph->GetXaxis()->SetTitle("Bin Number");
    meanGraph->GetYaxis()->SetTitle("Mean");
	//gPad->SetTopMargin(0.12);
    meanGraph->Draw("AP");
	canvas->cd(3);
	sigmaGraph->SetTitle("Sigma vs. Bin Number");
    sigmaGraph->GetXaxis()->SetTitle("Bin Number");
    sigmaGraph->GetYaxis()->SetTitle("Sigma");
	//gPad->SetTopMargin(0.12);
    sigmaGraph->Draw("AP");
	canvas->cd(4);
	gPad->SetLogz();
	hist2D->Draw("colz");
	//return 0;
}
