//sseeds 4.22.24

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TKey.h>
#include <TClass.h>
#include <iostream>
#include <vector>
#include <string>
#include <TLatex.h>
#include <algorithm>
#include <TRandom.h>
#include "../../include/gmn.h"
#include "../../src/jsonmgr.C"

//fitranges
//sbs9 70p: -1.8 to 0.7
//sbs8 50p: -1.4 to 0.7
//sbs8 70p: -1.8 to 0.7
//sbs8 100p: -2.1 to 0.7
//sbs4 30p: -1.6 to 0.8
//sbs4 50p: -2.1 to 0.7
//sbs7 85p: -1.4 to 0.5

//Fit range override options
double hcalfit_l = -2.1; //lower fit/bin limit for hcal dx plots (m)
double hcalfit_h = 0.7; //upper fit/bin limit for hcal dx plots (m)

double xrange_min_W2 = 0.0;
double xrange_max_W2 = 2.5;

//Fit options
std::string fitopt = "RMQ0";

//fit min entries
int minEvents = 500;

//Total fits using Interpolate with elastic signal histo and 4th order poly fit to bg
TH1D *hdx_p;
TH1D *hdx_n;
TH1D *hdx_inel;
TH1D *hdx_dyanti;
TH1D *hdx_coinanti;
TH1D *hW2;
TH1D *hW2_inel;
TH1D *hW2dy;
TH1D *hW2dy_inel;
TH1D *hW2spot;
TH1D *hW2spot_inel;

Double_t fitW2Shift(double *x, double *par){
  // MC float params
  double W2_scale = par[0];
  double bg_scale = par[1];
  double W2_shift = par[2];
  double bg_shift = par[3];  

  // Apply shifts before interpolation
  double W2 = W2_scale * hW2->Interpolate(x[0] - W2_shift);
  double bg = bg_scale * hW2_inel->Interpolate(x[0] - bg_shift);

  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return W2 + bg;
}

Double_t fitW2dyShift(double *x, double *par){
  // MC float params
  double W2_scale = par[0];
  double bg_scale = par[1];
  double W2_shift = par[2];
  double bg_shift = par[3];  

  // Apply shifts before interpolation
  double W2 = W2_scale * hW2dy->Interpolate(x[0] - W2_shift);
  double bg = bg_scale * hW2dy_inel->Interpolate(x[0] - bg_shift);

  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return W2 + bg;
}

Double_t fitW2spotShift(double *x, double *par){
  // MC float params
  double W2_scale = par[0];
  double bg_scale = par[1];
  double W2_shift = par[2];
  double bg_shift = par[3];  

  // Apply shifts before interpolation
  double W2 = W2_scale * hW2spot->Interpolate(x[0] - W2_shift);
  double bg = bg_scale * hW2spot_inel->Interpolate(x[0] - bg_shift);

  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return W2 + bg;
}

Double_t fitdxShift(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron + fits::g_p2fit_cd(x, &par[4]);
}

Double_t fitdxShiftInel(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  
  double bg_scale = par[4];

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);
  double bg = bg_scale * hdx_inel->Interpolate(x[0] - dx_shift_n);
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron + bg;
}

// Forward declarations
void handleError(TFile *file1, TFile *file2, std::string marker);
void handleError(TFile *file1, std::string marker);
std::vector<std::pair<double, double>> fitAndFineFit(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, double pshift, double nshift, const std::string& fitOptions = "RBMQ0");
TH1D* calculateResiduals(TH1D* hData, TF1* fit, const char* residualName);
void plotFitWithResiduals(TCanvas* canvas, TH1D* data, Double_t (*fitFunc)(double*, double*), TH1D* mcSignal, TH1D* mcBg, const char* fitOpt, const char* title);
void plotFitWithoutResiduals(TCanvas* canvas, TH1D* data, Double_t (*fitFunc)(double*, double*), TH1D* mcSignal, TH1D* mcBg, const char* fitOpt, const char* title);

//main. kine=kinematic, mag=fieldsetting, pass=pass#, sb_min/max=sidebandlimits, shiftX=shifttodxdata, N=cutvarsliceN, sliceCutMax=NCutsFromZeroTosliceCutMax
void W2datamc(int kine=4, 
	      int mag=30, 
	      int pass=2, 
	      bool bestclus=true, 
	      bool thin=true,
	      bool wide=false,
	      bool effz=true,
	      bool alt = true) {
  //set draw params
  gStyle->SetPalette(55);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);

  std::string bestclus_word = "";
  if(bestclus)
    bestclus_word = "_bc";

  std::string thin_word = "";
  if(thin)
    thin_word = "_thin";

  std::string alt_word = "";
  if(alt)
    alt_word = "_alt";

  std::string wide_word = "";
  if(wide)
    wide_word = "_widecut";

  std::string effz_word = "";
  if(effz)
    effz_word = "_effz";

  std::string basePath = "/lustre19/expphy/volatile/halla/sbs/seeds";
  std::string finPath = Form("%s/gmn_analysis/dx_correlations%s_sbs%d_mag%d_pass%d%s%s%s.root", basePath.c_str(), bestclus_word.c_str(), kine, mag, pass, thin_word.c_str(), wide_word.c_str(), effz_word.c_str());
  std::string fmcinPath = Form("%s/gmn_analysis/dx_mc_sbs%d_mag%d_pass%d%s%s%s%s.root", basePath.c_str(), kine, mag, pass, thin_word.c_str(), alt_word.c_str(), wide_word.c_str(), effz_word.c_str());
  std::string foutPath = Form("%s/gmn_analysis/W2datamc_gmn_sbs%d_mag%d_pass%d%s%s%s.root", basePath.c_str(), kine, mag, pass, bestclus_word.c_str(), wide_word.c_str(), effz_word.c_str());

  TFile* inputFile = new TFile(finPath.c_str());
  if (!inputFile || inputFile->IsZombie()) 
    handleError(inputFile,"inputFile");

  cout << endl << endl << "In path: " << finPath << endl;
  cout << "In MC path: " << fmcinPath << endl << endl;

  //Get all histograms
  TH1D *hdx_data = dynamic_cast<TH1D*>(inputFile->Get("hdx_allcut"));
  hdx_data->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  TH1D *hW2_data = dynamic_cast<TH1D*>(inputFile->Get("hW2_nodycut"));
  hW2_data->GetXaxis()->SetRangeUser(xrange_min_W2,xrange_max_W2);
  TH1D *hW2_datady = dynamic_cast<TH1D*>(inputFile->Get("hW2"));
  hW2_datady->GetXaxis()->SetRangeUser(xrange_min_W2,xrange_max_W2);
  TH1D *hW2_dataspot = dynamic_cast<TH1D*>(inputFile->Get("hW2_spotcut"));
  hW2_dataspot->GetXaxis()->SetRangeUser(xrange_min_W2,xrange_max_W2);

  TH2D *hdx_vs_W2_data = dynamic_cast<TH2D*>(inputFile->Get("hist_W2"));
  hdx_vs_W2_data->GetYaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  int hdx_vs_W2_data_N = hdx_vs_W2_data->GetEntries();
  std::string W2Cuts = hdx_vs_W2_data->GetTitle();
  cout << endl << "Opened dx vs W2 with cuts: " << W2Cuts << endl << endl; 

  hdx_dyanti = dynamic_cast<TH1D*>(inputFile->Get("hdx_dyanti"));
  hdx_dyanti->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  hdx_coinanti = dynamic_cast<TH1D*>(inputFile->Get("hdx_coinanti"));
  hdx_coinanti->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  if (!hdx_data || !hdx_dyanti || !hdx_coinanti) 
    handleError(inputFile,"hdx_data");

  TFile* inputFileMC = new TFile(fmcinPath.c_str(), "READ");
  if (!inputFileMC || inputFileMC->IsZombie())
    handleError(inputFile,inputFileMC,"inputFileMC");

  //fix the interpolate functions for W2
  hW2 = dynamic_cast<TH1D*>(inputFileMC->Get("hW2_nodycut"));
  hW2->GetXaxis()->SetRangeUser(xrange_min_W2,xrange_max_W2);
  hW2_inel = dynamic_cast<TH1D*>(inputFileMC->Get("hW2_inel_dy"));
  hW2_inel->Scale(10e33);
  hW2_inel->GetXaxis()->SetRangeUser(xrange_min_W2,xrange_max_W2);
  hW2dy = dynamic_cast<TH1D*>(inputFileMC->Get("hW2"));
  hW2dy->GetXaxis()->SetRangeUser(xrange_min_W2,xrange_max_W2);
  hW2dy_inel = dynamic_cast<TH1D*>(inputFileMC->Get("hW2_inel"));
  hW2dy_inel->Scale(10e33);
  hW2dy_inel->GetXaxis()->SetRangeUser(xrange_min_W2,xrange_max_W2);
  hW2spot = dynamic_cast<TH1D*>(inputFileMC->Get("hW2_spotcut"));
  hW2spot->GetXaxis()->SetRangeUser(xrange_min_W2,xrange_max_W2);
  hW2spot_inel = dynamic_cast<TH1D*>(inputFileMC->Get("hW2_inel_spotcut"));
  hW2spot_inel->Scale(10e33);
  hW2spot_inel->GetXaxis()->SetRangeUser(xrange_min_W2,xrange_max_W2);

  //make clones for quality checks
  TH1D *hW2_mc_inel = (TH1D*)hW2_inel->Clone("hW2_mc_inel");
  TH1D *hW2_mcdy_inel = (TH1D*)hW2dy_inel->Clone("hW2_mcdy_inel");
  TH1D *hW2_mcspot_inel = (TH1D*)hW2spot_inel->Clone("hW2_mcspot_inel");

  TH1D *hW2_mc = (TH1D*)hW2->Clone("hW2_mc");
  TH1D *hW2_mcdy = (TH1D*)hW2dy->Clone("hW2_mcdy");
  TH1D *hW2_mcspot = (TH1D*)hW2spot->Clone("hW2_mcspot");

  //fix the interpolate functions for dx
  hdx_p = dynamic_cast<TH1D*>(inputFileMC->Get("hdx_p"));
  hdx_p->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  hdx_n = dynamic_cast<TH1D*>(inputFileMC->Get("hdx_n"));
  hdx_n->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  TH2D *hdx_vs_W2_n = dynamic_cast<TH2D*>(inputFileMC->Get("hist_W2_n"));
  hdx_vs_W2_n->GetYaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  int hdx_vs_W2_n_N = hdx_vs_W2_n->GetEntries();

  TH2D *hdx_vs_W2_p = dynamic_cast<TH2D*>(inputFileMC->Get("hist_W2_p"));
  hdx_vs_W2_p->GetYaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  int hdx_vs_W2_p_N = hdx_vs_W2_p->GetEntries();

  hdx_inel = dynamic_cast<TH1D*>(inputFileMC->Get("hdx_inel"));
  hdx_inel->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  hdx_inel->Scale(10e33); //account for lack of overall normalization from g4sbs generator
  if (!hdx_p || !hdx_n || !hdx_inel) 
    handleError(inputFile,inputFileMC,"hdx");

  TFile* outputFile = new TFile(foutPath.c_str(), "RECREATE");

  TCanvas *canvas = new TCanvas("canvas", "Data and MC Histograms", 1800, 1200);
  canvas->Divide(3, 3);  // Adjust the grid size according to the number of histograms

  //Get the W2 histogram from data
  //TH1D* hW2_data = hdx_vs_W2_data->ProjectionX("hW2_data");

  // Pad 1: W2
  canvas->cd(1);
  hW2_data->SetLineColor(kBlack);
  hW2_data->SetLineWidth(2);
  hW2_data->Draw();
  hW2_mc->SetLineColor(kGreen);
  hW2_mc->Scale(0.47);
  hW2_mc->Draw("same");
  hW2_mc_inel->SetLineColor(kRed-5);
  hW2_mc_inel->Scale(23);
  hW2_mc_inel->Draw("same");

  // Pad 2: W2 dy cut
  canvas->cd(2);
  hW2_datady->SetLineColor(kBlack);
  hW2_datady->SetLineWidth(2);
  hW2_datady->Draw();
  hW2_mcdy->SetLineColor(kGreen);
  hW2_mcdy->Scale(0.47);
  hW2_mcdy->Draw("same");
  hW2_mcdy_inel->SetLineColor(kRed-5);
  hW2_mcdy_inel->Scale(20);
  hW2_mcdy_inel->Draw("same");

  // Pad 3: W2 spot cut
  canvas->cd(3);
  hW2_dataspot->SetLineColor(kBlack);
  hW2_dataspot->SetLineWidth(2);
  hW2_dataspot->Draw();
  hW2_mcspot->SetLineColor(kGreen);
  hW2_mcspot->Scale(1.7);
  hW2_mcspot->Draw("same");
  hW2_mcspot_inel->SetLineColor(kRed-5);
  hW2_mcspot_inel->Scale(20);
  hW2_mcspot_inel->Draw("same");

  // Pad 4: hdx_data
  canvas->cd(4);
  hW2_data->SetLineColor(kBlack);
  hW2_data->SetLineWidth(2);
  hW2_data->Draw();


  // Pad 5: hdx_vs_W2_data
  canvas->cd(5);
  hW2_mc->SetLineColor(kGreen);
  hW2_mc->Draw("");
  hW2_mc_inel->SetLineColor(kRed-5);
  hW2_mc_inel->Draw("same");

  // Pad 6: hdx_n
  canvas->cd(6);
  hdx_vs_W2_data->Draw("COLZ");

  // Pad 7: hdx_mc
  canvas->cd(7);
  hdx_p->SetLineColor(kRed);
  hdx_p->Draw("p");
  hdx_n->SetLineColor(kBlue);
  hdx_n->Draw("p same");

  // Pad 8: hdx_vs_W2
  canvas->cd(8);
  hdx_data->Draw();

  // Pad 9: hdx_inel
  canvas->cd(9);
  hdx_inel->Draw();

  // Update the canvas to reflect the drawings
  canvas->Update();

  // // First canvas
  // TCanvas *cFITT = new TCanvas("cFITT", "Data and MC Fitted (W2)", 1200, 800);
  // //TF1* W2mcfit = new TF1("W2mcfit", fitW2Shift, hW2_data->GetXaxis()->GetXmin(), hW2_data->GetXaxis()->GetXmax(), 4);
  // plotFitWithResiduals(cFITT, hW2_data, fitW2Shift, hW2, hW2_inel, fitopt.c_str(), "W^{2} Data Fit with MC");

  // // Second canvas
  // TCanvas *cFITTdy = new TCanvas("cFITTdy", "Data and MC Fitted (W2dy)", 1200, 800);
  // plotFitWithResiduals(cFITTdy, hW2_datady, fitW2dyShift, hW2dy, hW2dy_inel, fitopt.c_str(), "W^{2} Data Fit with MC (dy)");

  // // Third canvas
  // TCanvas *cFITTspot = new TCanvas("cFITTspot", "Data and MC Fitted (W2spot)", 1200, 800);
  // plotFitWithResiduals(cFITTspot, hW2_dataspot, fitW2spotShift, hW2spot, hW2spot_inel, fitopt.c_str(), "W^{2} Data Fit with MC (spot)");

  // Fourth canvas
  TCanvas *csmfit = new TCanvas("csmfit", "Data and MC Fitted", 1200, 800);
  plotFitWithoutResiduals(csmfit, hW2_data, fitW2Shift, hW2, hW2_inel, fitopt.c_str(), "W^{2} Data Fit with MC (elastic cuts)");

  // Fourth canvas
  TCanvas *csmfitdy = new TCanvas("csmfitdy", "Data and MC Fitted (dy cuts)", 1200, 800);
  plotFitWithoutResiduals(csmfitdy, hW2_datady, fitW2dyShift, hW2dy, hW2dy_inel, fitopt.c_str(), "W^{2} Data Fit with MC (elastic and dy cuts)");

  

}


/////////////
/////////////
/////////////
/////////////
/////////////


void handleError(TFile *file1, TFile *file2, std::string marker) {
    if (file1) file1->Close();
    if (file2) file2->Close();
    std::cerr << "Error: File opening or histogram retrieval failed at " << marker << "." << std::endl;
}

void handleError(TFile *file1, std::string marker) {
    if (file1) file1->Close();
    std::cerr << "Error: File opening or histogram retrieval failed at " << marker << "." << std::endl;
}

//Function to fit then fine fit and return all fit parameters
std::vector<std::pair<double, double>> fitAndFineFit(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, double pshift, double nshift, const std::string& fitOptions = "RBMQ0") {
  TF1* fit = new TF1(fitName.c_str(), fitFormula.c_str(), hcalfit_l, hcalfit_h, paramCount);
  //fit->SetNpx(5000);
  for (int i=0; i<paramCount; ++i){ //reset parameters/errors for this set
    fit->SetParameter(i,0);
    fit->SetParError(i,0);
  }
  if(!(pshift==0&&nshift==0)){
    fit->FixParameter(2,pshift);
    fit->FixParameter(3,nshift);
  }

  histogram->Fit(fit, fitOptions.c_str());

  std::vector<std::pair<double, double>> parametersAndErrors(paramCount);
  for (int i = 0; i < paramCount; ++i) {
    parametersAndErrors[i].first = fit->GetParameter(i); // Parameter value
    parametersAndErrors[i].second = fit->GetParError(i); // Parameter error
  }

  // Fine Fit
  TF1* fineFit = new TF1((fitName + "_fine").c_str(), fitFormula.c_str(), hcalfit_l, hcalfit_h, paramCount);
  std::vector<double> fineFitInitialParams(paramCount);
  for (int i = 0; i < paramCount; ++i) {
    fineFitInitialParams[i] = parametersAndErrors[i].first;
  }

  fineFit->SetParameters(fineFitInitialParams.data());
  if(!(pshift==0&&nshift==0)){
    fit->FixParameter(2,pshift);
    fit->FixParameter(3,nshift);
  }

  histogram->Fit(fineFit, fitOptions.c_str());

  cout << fineFit->GetParameter(6) << endl;

  // Update parameters and errors with fine fit results
  for (int i = 0; i < paramCount; ++i) {
    parametersAndErrors[i].first = fineFit->GetParameter(i); // Fine fit parameter value
    parametersAndErrors[i].second = fineFit->GetParError(i); // Fine fit parameter error
  }

  fitqual.first = fineFit->GetChisquare();
  fitqual.second = fineFit->GetNDF();

  delete fit; // Delete fit to avoid carry-over
  delete fineFit; // Clean up
  return parametersAndErrors;
}



TH1D* calculateResiduals(TH1D* hData, TF1* fit, const char* residualName) {
  // Clone the data histogram to create the residuals histogram
  TH1D* hResiduals = (TH1D*)(hData->Clone(residualName));
  hResiduals->Reset();  // Clear the histogram

  // Loop over each bin and calculate the residuals
  for (int bin = 1; bin <= hData->GetNbinsX(); ++bin) {
    double x = hData->GetXaxis()->GetBinCenter(bin);
    double dataValue = hData->GetBinContent(bin);
    double fitValue = fit->Eval(x);

    // Calculate the residual as (data - fit)
    double residual = dataValue - fitValue;
    double error = hData->GetBinError(bin); // Assuming Poisson statistics, you can also add fit and bg errors if needed

    hResiduals->SetBinContent(bin, residual);
    hResiduals->SetBinError(bin, error);
  }

  return hResiduals;
}


void plotFitWithResiduals(TCanvas* canvas, TH1D* data, Double_t (*fitFunc)(double*, double*), TH1D* mcSignal, TH1D* mcBg, const char* fitOpt, const char* title) {

  // Generate a unique identifier for the fit function name
  TString uniqueID = TString::Format("_%u", gRandom->Integer(1000000));

  canvas->cd();

  // Create a pad for the main plot
  TPad *pad1 = new TPad("pad1", "Main Plot", 0.0, 0.3, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "Residuals", 0.0, 0.0, 1.0, 0.3);

  pad1->SetBottomMargin(0);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad1->Draw();
  pad2->Draw();

  pad1->cd();

  gPad->SetGridx();
  gPad->SetGridy();

  data->SetTitle(title);
  data->SetTitleFont(132);
  //data->SetMarkerStyle(7);
  data->SetLineColor(kBlack);
  data->SetLineWidth(2);
  data->Draw("hist E");

  TF1* fit = new TF1("fit", fitFunc, data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax(), 4);

  // Set parameter limits
  //fit->SetParLimits(2, -0.1, 0.1);
  fit->SetParLimits(3, -0.1, 0.1);

  // Fit the data
  data->Fit(fit, fitOpt);

  std::vector<double> fit_pars;
  for (int i = 0; i < fit->GetNpar(); ++i)
    fit_pars.push_back(fit->GetParameter(i));

  fit->SetLineColor(kGreen);
  fit->SetNpx(1000);
  fit->Draw("same");

  TH1D *fit_signal = util::shiftHistogramX(mcSignal, fit_pars[2]);
  fit_signal->Scale(fit_pars[0]);
  fit_signal->SetLineColor(kBlue);
  fit_signal->Draw("same");

  TH1D *fit_bg = util::shiftHistogramX(mcBg, fit_pars[3]);
  fit_bg->Scale(fit_pars[1]);
  fit_bg->SetLineColor(kRed - 5);
  fit_bg->Draw("same");

  // Calculate Chi-Square/NDF
  double chi2 = fit->GetChisquare();
  int ndf = fit->GetNDF();

  // Create a legend
  TLegend *legend = new TLegend(0.11, 0.7, 0.5, 0.9);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(data, "Data", "lep");
  legend->AddEntry(fit, "Fit (MC Signal + BG)", "l");
  legend->AddEntry(fit_signal, TString::Format("MC Signal (Scale: %.2f, Shift: %.2f)", fit_pars[0], fit_pars[2]), "l");
  legend->AddEntry(fit_bg, TString::Format("MC Background (Scale: %.2f, Shift: %.2f)", fit_pars[1], fit_pars[3]), "l");
  legend->AddEntry((TObject*)0, TString::Format("#chi^{2}/ndf: %.2f/%d", chi2, ndf), "");
  legend->Draw();

  // Create a pad for the residuals
  pad2->cd();
  gPad->SetGridx();

  // Calculate residuals (Data - Fit)
  TH1D *hRes = calculateResiduals(data, fit, Form("hResiduals_%s", uniqueID.Data()));
  double W2MaxValue = data->GetMaximum();
  hRes->GetYaxis()->SetRangeUser(-W2MaxValue / 4, W2MaxValue / 4);

  // Draw the residuals
  hRes->SetTitle("");
  hRes->GetYaxis()->SetTitle("Residuals");
  hRes->GetYaxis()->SetTitleSize(0.1);
  hRes->GetYaxis()->SetTitleOffset(0.5);
  hRes->GetYaxis()->SetLabelSize(0.1);
  hRes->GetXaxis()->SetTitleSize(0.1);
  hRes->GetXaxis()->SetLabelSize(0.1);
  hRes->GetXaxis()->SetTitle("x_{hcal}-x_{exp}");
  hRes->GetXaxis()->SetTitleFont(132);
  hRes->GetXaxis()->SetLabelFont(132);
  hRes->GetYaxis()->SetTitleFont(132);
  hRes->GetYaxis()->SetLabelFont(132);
  hRes->SetLineColor(kBlack);
  hRes->SetLineWidth(1);
  hRes->Draw("E");

  // Fill background
  double y_min = hRes->GetMinimum();
  double y_max = hRes->GetMaximum();

  // Create a TPave to fill the region between the axes
  TPave *pave = new TPave(xrange_min_W2, y_min, xrange_max_W2, y_max, 0, "NB");
  pave->SetFillColorAlpha(kGray, 0.35);
  pave->SetFillStyle(3001);
  pave->Draw("same");

  // Draw the residuals
  hRes->Draw("E2 same");
  hRes->Draw("E same");

  // Draw a line at y=0 for residuals
  TF1 *zeroLine = new TF1("zeroLine", "0", xrange_min_W2, xrange_max_W2);
  zeroLine->SetLineColor(kRed);
  zeroLine->SetLineStyle(2);
  zeroLine->Draw("same");

  // Update the canvas
  canvas->cd();
  canvas->Update();
}



void plotFitWithoutResiduals(TCanvas* canvas, TH1D* data, Double_t (*fitFunc)(double*, double*), TH1D* mcSignal, TH1D* mcBg, const char* fitOpt, const char* title) {

  // Generate a unique identifier for the fit function name
  TString uniqueID = TString::Format("_%u", gRandom->Integer(1000000));

  canvas->cd();

  // Create a pad for the main plot
  TPad *pad1 = new TPad("pad1", "Main Plot", 0.0, 0.0, 1.0, 1.0);
  pad1->SetBottomMargin(0.15);
  pad1->SetLeftMargin(0.15);
  pad1->SetRightMargin(0.05);
  pad1->SetTopMargin(0.1);
  pad1->Draw();
  pad1->cd();

  gPad->SetGridx();
  gPad->SetGridy();

  data->SetTitle(title);
  data->SetTitleFont(132);
  data->SetMarkerStyle(20);
  data->SetMarkerSize(0.5);
  data->SetLineColor(kBlack);
  data->SetLineWidth(1);
  data->Draw("hist E");

  TF1* fit = new TF1("fit", fitFunc, data->GetXaxis()->GetXmin(), data->GetXaxis()->GetXmax(), 4);

  // Set parameter limits
  fit->SetParLimits(2, -0.1, 0.1);
  fit->SetParLimits(3, -0.1, 0.1);

  // Fit the data
  data->Fit(fit, fitOpt);

  std::vector<double> fit_pars;
  for (int i = 0; i < fit->GetNpar(); ++i)
    fit_pars.push_back(fit->GetParameter(i));

  // Shade the fit area in green
  fit->SetLineColor(kGreen);
  fit->SetFillColorAlpha(kGreen, 0.35);
  fit->SetFillStyle(1001);
  fit->SetNpx(1000);
  fit->Draw("same FC");

  TH1D *fit_signal = util::shiftHistogramX(mcSignal, fit_pars[2]);
  fit_signal->Scale(fit_pars[0]);
  fit_signal->SetMarkerStyle(21); // Different filled circle marker
  fit_signal->SetMarkerColor(kBlue);
  fit_signal->SetMarkerSize(0.5);
  fit_signal->SetLineColor(kBlue);
  fit_signal->Draw("P same");

  TH1D *fit_bg = util::shiftHistogramX(mcBg, fit_pars[3]);
  fit_bg->Scale(fit_pars[1]);
  fit_bg->SetMarkerStyle(22); // Different filled circle marker
  fit_bg->SetMarkerColor(kRed - 5);
  fit_bg->SetMarkerSize(0.5);
  fit_bg->SetLineColor(kRed - 5);
  fit_bg->Draw("P same");

  // Create a legend
  TLegend *legend = new TLegend(0.15, 0.7, 0.5, 0.9);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.04);
  legend->SetTextColor(kViolet);
  legend->AddEntry(data, "Data", "lep");
  legend->AddEntry(fit, "Fit (MC Signal + BG)", "f"); // "f" for filled area
  legend->AddEntry(fit_signal, "MC Signal", "p");
  legend->AddEntry(fit_bg, "MC Background", "p");
  legend->Draw();

  // Update the canvas
  canvas->cd();
  canvas->Update();
}

