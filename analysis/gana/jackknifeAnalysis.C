//sseeds 4.22.24: Script to run over jacknife samples and extract Rsf

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
#include "../../include/gmn.h"
#include "../../src/jsonmgr.C"

//fitranges
//sbs9 70p: -1.8 to 0.7
//sbs8 70p: -1.8 to 0.7
//sbs4 30p: -1.7 to 0.7, shiftX=0.0, neutronshift=-0.05
//sbs4 50p: -2.1 to 0.7, shiftX=0.05, neutronshift=0.0
//sbs7 85p: -1.4 to 0.5

//Fit range override options
double hcalfit_l = -1.8; //lower fit/bin limit for hcal dx plots (m) sbs4, 50p
double hcalfit_h = 0.7; //upper fit/bin limit for hcal dx plots (m) sbs4, 50p

//Plot range option shared by all fits
double hcalr_l = -1.8; //lower fit/bin limit for hcal dx plots (m) sbs4, 50p
double hcalr_h = 0.7; //upper fit/bin limit for hcal dx plots (m) sbs4, 50p

//Fit options
std::string fitopt = "RMQ0";

//Total fits using Interpolate with elastic signal histo and 4th order poly fit to bg
TH1D *hdx_p;
TH1D *hdx_n;

//Total datasets
double N_datasets = 11;

Double_t fitFullShift_p2(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron + fits::g_p2fit(x, &par[4]);
}

// Forward declarations
void handleError(TFile *file1, TFile *file2, std::string marker);
void handleError(TFile *file1, std::string marker);
std::vector<std::pair<double, double>> fitAndFineFit(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, const std::string& fitOptions = "RBMQ0");
double weightedMean(const std::vector<double>& values, const std::vector<double>& errors);
double weightedStdDev(const std::vector<double>& values, const std::vector<double>& errors, double weightedMean);
double calculateMean(const std::vector<double>& values);
double calculateStdDev(const std::vector<double>& values, double mean);
double calculateJKErr(const std::vector<double>& values, double mean, double num_datasets);

//main. kine=kinematic, mag=fieldsetting, pass=pass#
void jackknifeAnalysis(int kine=8, 
		       int mag=70, 
		       int pass=2) {

  //set draw params
  gStyle->SetPalette(55);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);

  // Set up configuration and tune objects to load analysis parameters
  SBSconfig config(kine,0);

  // Obtain configuration pars from config class
  double E_beam = config.GetEbeam();
  double BB_angle = config.GetBBtheta_rad(); // In radians

  std::string thin_word = "_thin";
  std::string alt_word = "_alt";
  //std::string alt_word = "";

  //Paths to files
  std::string basePath = "/lustre19/expphy/volatile/halla/sbs/seeds";
  std::string finPath = Form("%s/gmn_analysis/jacknife_samples_sbs%d_mag%d_pass%d.root", basePath.c_str(), kine, mag, pass);
  std::string fmcinPath = Form("%s/gmn_analysis/dx_mc_sbs%d_mag%d_pass%d%s%s.root", basePath.c_str(), kine, mag, pass, thin_word.c_str(), alt_word.c_str());
  std::string foutPath = Form("%s/gmn_analysis/jackknife_out_sbs%d_mag%d_pass%d.root", basePath.c_str(), kine, mag, pass);

  //get data file
  TFile* inputFile = new TFile(finPath.c_str());
  if (!inputFile || inputFile->IsZombie()) 
    handleError(inputFile,"inputFile");

  //get MC file
  TFile* inputFileMC = new TFile(fmcinPath.c_str(), "READ");
  if (!inputFileMC || inputFileMC->IsZombie())
    handleError(inputFile,inputFileMC,"inputFileMC");

  //fix the interpolate functions
  hdx_p = dynamic_cast<TH1D*>(inputFileMC->Get("hdx_p"));
  hdx_p->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  hdx_n = dynamic_cast<TH1D*>(inputFileMC->Get("hdx_n"));
  hdx_n->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  std::map<int, TH1D*> histograms_stat;
  std::map<int, TH1D*> histograms_syst;
  TIter next(inputFile->GetListOfKeys());
  TKey* key;

  while ((key = dynamic_cast<TKey*>(next()))) {
    TH1D* hist = dynamic_cast<TH1D*>(key->ReadObj());
    if (hist) {
      std::string histName = hist->GetName();
      hist->GetXaxis()->SetRangeUser(hcalr_l,hcalr_h);

      if (histName.find("hdx_stat") != std::string::npos) {
	histograms_stat[std::atoi(histName.substr(9).c_str())] = hist;
      } else if (histName.find("hdx_syst") != std::string::npos) {
	histograms_syst[std::atoi(histName.substr(9).c_str())] = hist;
      }
    }
  }

  //make output file
  TFile* outputFile = new TFile(foutPath.c_str(), "RECREATE");

  struct HistogramData {
    TH1D* histogram;
    double histogramNumber;
    double protonScale;
    double neutronScale;
    double R_sf;
    double R_sf_err;
    std::vector<double> params;
    std::vector<double> paramErrors;
    TH1D* pMCslice = nullptr;  // Proton MC histogram slice
    TH1D* nMCslice = nullptr;  // Neutron MC histogram slice
    TF1* bgfit = nullptr;      // Background fit function

    HistogramData() : histogram(nullptr), histogramNumber(0), protonScale(0), neutronScale(0), R_sf(0), R_sf_err(0) {}
    ~HistogramData() {
      delete pMCslice;
      delete nMCslice;
      delete bgfit;
    }
  };

  std::vector<HistogramData> statData;
  std::vector<HistogramData> systData;

  std::vector<HistogramData> JKData;

  //Build equal stats set
  for (auto& pair : histograms_stat) {
    TH1D* hist = pair.second;

    HistogramData data;
    data.histogram = hist;
    data.histogramNumber = pair.first;

    std::pair<double,double> fitParStatP2;
    auto fitParStatP2_vector = fitAndFineFit(hist, "statFit", "fitFullShift_p2", 7, hcalfit_l, hcalfit_h, fitParStatP2, fitopt.c_str());

    data.protonScale = fitParStatP2_vector[0].first;
    data.neutronScale = fitParStatP2_vector[1].first;
    data.R_sf = data.neutronScale / data.protonScale;
    data.R_sf_err = sqrt(pow(fitParStatP2_vector[1].second/data.neutronScale, 2) + pow(fitParStatP2_vector[0].second/data.protonScale, 2)) * data.R_sf;

    // Store parameters and errors
    for (const auto& p : fitParStatP2_vector) {
      data.params.push_back(p.first);
      data.paramErrors.push_back(p.second);
    }

    statData.push_back(data);
  }

  //Build set by run number
  for (auto& pair : histograms_syst) {
    TH1D* hist = pair.second;

    HistogramData data;
    data.histogram = hist;
    data.histogramNumber = pair.first;

    std::pair<double,double> fitParSystP2;
    auto fitParSystP2_vector = fitAndFineFit(hist, "statFit", "fitFullShift_p2", 7, hcalfit_l, hcalfit_h, fitParSystP2, fitopt.c_str());

    data.protonScale = fitParSystP2_vector[0].first;
    data.neutronScale = fitParSystP2_vector[1].first;
    data.R_sf = data.neutronScale / data.protonScale;
    data.R_sf_err = sqrt(pow(fitParSystP2_vector[1].second/data.neutronScale, 2) + pow(fitParSystP2_vector[0].second/data.protonScale, 2)) * data.R_sf;

    // Store parameters and errors
    for (const auto& p : fitParSystP2_vector) {
      data.params.push_back(p.first);
      data.paramErrors.push_back(p.second);
    }

    systData.push_back(data);
  }

  //Build jackknife set
  for (auto it = histograms_stat.begin(); it != histograms_stat.end(); ++it) {
    TH1D* combinedHistogram = nullptr; // This histogram will be the combination of all except one

    // Combine all histograms except the one pointed to by 'it'
    for (auto jt = histograms_stat.begin(); jt != histograms_stat.end(); ++jt) {
      if (jt == it) continue; // Skip the current histogram in the outer loop
      if (!combinedHistogram) {
	combinedHistogram = (TH1D*)jt->second->Clone(Form("combined_%d", jt->first));
      } else {
	combinedHistogram->Add(jt->second); // Ensure proper scaling if necessary
      }
    }

    // Now fit the combined histogram
    if (combinedHistogram) {
      HistogramData data;
      data.histogram = combinedHistogram;
      data.histogramNumber = it->first; // Use the number of the excluded histogram to identify this sample

      std::pair<double,double> fitPar;
      auto fitResults = fitAndFineFit(combinedHistogram, "jackknifeFit", "fitFullShift_p2", 7, hcalfit_l, hcalfit_h, fitPar, fitopt.c_str());

      data.protonScale = fitResults[0].first;
      data.neutronScale = fitResults[1].first;
      data.R_sf = data.neutronScale / data.protonScale;
      data.R_sf_err = sqrt(pow(fitResults[1].second/data.neutronScale, 2) + pow(fitResults[0].second/data.protonScale, 2)) * data.R_sf;

      // Store parameters and errors
      for (const auto& p : fitResults) {
	data.params.push_back(p.first);
	data.paramErrors.push_back(p.second);
      }

      JKData.push_back(data);
    }
  }

  // Plot fit results on stat set
  TCanvas* canvasStat = new TCanvas("canvasStat", "Fits to Equal Stats Histograms", 1800, 1200);
  int Nhist_stat = statData.size();
  canvasStat->Divide(sqrt(Nhist_stat), sqrt(Nhist_stat) + 1);  // Dynamically determine the number of pads

  for (int i = 0; i < Nhist_stat-1; ++i) { //Do not draw the last histogram with insufficient stats
    canvasStat->cd(i + 1);

    HistogramData& data = statData[i];
    data.histogram->Draw();  // Draw experimental histogram

    // Draw the corresponding MC histograms
    TH1D *pMCslice = util::shiftHistogramX(hdx_p, data.params[2]);
    pMCslice->Scale(data.protonScale);

    TH1D *nMCslice = util::shiftHistogramX(hdx_n, data.params[3]);
    nMCslice->Scale(data.neutronScale);

    pMCslice->Draw("same E");
    nMCslice->Draw("same E");

    // Add legend
    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(data.histogram, "Experimental Data", "l");
    legend->AddEntry(data.pMCslice, "Proton MC Fit", "l");
    legend->AddEntry(data.nMCslice, "Neutron MC Fit", "l");
    legend->AddEntry((TObject*)0, Form("Nev: %0.0f", data.histogram->GetEntries()), "");
    legend->AddEntry((TObject*)0, Form("Ratio: %0.4f", data.R_sf), "");
    legend->Draw();

    TF1 *mcfit = new TF1(Form("mcfit_%d", i), fitFullShift_p2, hcalfit_l, hcalfit_h, 7);
    for (int j = 0; j < 7; ++j) {
      mcfit->SetParameter(j, data.params[j]);
    }

    // Create a new TH1D or fill an existing one with the values from TF1
    int nbins = data.histogram->GetNbinsX();
    double x_low = data.histogram->GetXaxis()->GetXmin();
    double x_high = data.histogram->GetXaxis()->GetXmax();
    TH1D* hFromStatTF1 = new TH1D(Form("hFromStatTF1_%d", i), "Histogram from Stat TF1", nbins, x_low, x_high);

    for (int bin = 1; bin <= nbins; ++bin) {
      double binCenter = data.histogram->GetBinCenter(bin);
      double funcValue = mcfit->Eval(binCenter);
      hFromStatTF1->SetBinContent(bin, funcValue);
    }

    int transparentGreen = TColor::GetColorTransparent(kGreen, 0.3);
    hFromStatTF1->SetLineColor(kGreen);
    hFromStatTF1->SetLineWidth(0);
    hFromStatTF1->SetFillColor(transparentGreen);
    hFromStatTF1->SetFillStyle(1001);
    hFromStatTF1->Draw("SAME LF2");

    canvasStat->Update();
  }

  // Plot fit results on syst set (by run number)
  TCanvas* canvasSyst = new TCanvas("canvasSyst", "Fits to Run Number Histograms", 1800, 1200);
  int Nhist_syst = systData.size();
  canvasSyst->Divide(sqrt(Nhist_syst), sqrt(Nhist_syst) + 1);  // Dynamically determine the number of pads

  for (int i = 0; i < Nhist_syst-1; ++i) { //Do not draw the last histogram with insufficient stats
    canvasSyst->cd(i + 1);

    HistogramData& data = systData[i];
    data.histogram->SetTitle(Form("dx, run number %0.0f", data.histogramNumber));
    data.histogram->Draw();  // Draw experimental histogram

    // Draw the corresponding MC histograms
    TH1D *pMCslice = util::shiftHistogramX(hdx_p, data.params[2]);
    pMCslice->Scale(data.protonScale);

    TH1D *nMCslice = util::shiftHistogramX(hdx_n, data.params[3]);
    nMCslice->Scale(data.neutronScale);

    pMCslice->Draw("same E");
    nMCslice->Draw("same E");

    // Add legend
    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(data.histogram, "Experimental Data", "l");
    legend->AddEntry(data.pMCslice, "Proton MC Fit", "l");
    legend->AddEntry(data.nMCslice, "Neutron MC Fit", "l");
    legend->AddEntry((TObject*)0, Form("Nev: %0.0f", data.histogram->GetEntries()), "");
    legend->AddEntry((TObject*)0, Form("Ratio: %0.4f", data.R_sf), "");

    legend->Draw();

    TF1 *mcfit = new TF1(Form("mcfit_%d", i), fitFullShift_p2, hcalfit_l, hcalfit_h, 7);
    for (int j = 0; j < 7; ++j) {
      mcfit->SetParameter(j, data.params[j]);
    }

    // Create a new TH1D or fill an existing one with the values from TF1
    int nbins = data.histogram->GetNbinsX();
    double x_low = data.histogram->GetXaxis()->GetXmin();
    double x_high = data.histogram->GetXaxis()->GetXmax();
    TH1D* hFromSystTF1 = new TH1D(Form("hFromSystTF1_%d", i), "Histogram from Syst TF1", nbins, x_low, x_high);

    for (int bin = 1; bin <= nbins; ++bin) {
      double binCenter = data.histogram->GetBinCenter(bin);
      double funcValue = mcfit->Eval(binCenter);
      hFromSystTF1->SetBinContent(bin, funcValue);
    }

    int transparentMagenta = TColor::GetColorTransparent(kMagenta, 0.3);
    hFromSystTF1->SetLineColor(kMagenta);
    hFromSystTF1->SetLineWidth(0);
    hFromSystTF1->SetFillColor(transparentMagenta);
    hFromSystTF1->SetFillStyle(1001);
    hFromSystTF1->Draw("SAME LF2");

    canvasSyst->Update();
  }

  // Plot fit results on syst set
  TCanvas* canvasJK = new TCanvas("canvasJK", "Fits to Jackknife Histograms", 1800, 1200);
  int Nhist_JK = JKData.size();
  canvasJK->Divide(sqrt(Nhist_JK), sqrt(Nhist_JK) + 1);  // Dynamically determine the number of pads

  for (int i = 0; i < Nhist_JK; ++i) {
    canvasJK->cd(i + 1);

    HistogramData& data = JKData[i];
    data.histogram->SetTitle(Form("dx, remove range %0.0f", data.histogramNumber));
    data.histogram->Draw();  // Draw experimental histogram

    // Draw the corresponding MC histograms
    TH1D *pMCslice = util::shiftHistogramX(hdx_p, data.params[2]);
    pMCslice->Scale(data.protonScale);

    TH1D *nMCslice = util::shiftHistogramX(hdx_n, data.params[3]);
    nMCslice->Scale(data.neutronScale);

    pMCslice->Draw("same E");
    nMCslice->Draw("same E");

    // Add legend
    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->AddEntry(data.histogram, "Experimental Data", "l");
    legend->AddEntry(data.pMCslice, "Proton MC Fit", "l");
    legend->AddEntry(data.nMCslice, "Neutron MC Fit", "l");
    legend->AddEntry((TObject*)0, Form("Nev: %0.0f", data.histogram->GetEntries()), "");
    legend->AddEntry((TObject*)0, Form("Ratio: %0.4f", data.R_sf), "");
    legend->Draw();

    TF1 *mcfit = new TF1(Form("mcfit_%d", i), fitFullShift_p2, hcalfit_l, hcalfit_h, 7);
    for (int j = 0; j < 7; ++j) {
      mcfit->SetParameter(j, data.params[j]);
    }

    // Create a new TH1D or fill an existing one with the values from TF1
    int nbins = data.histogram->GetNbinsX();
    double x_low = data.histogram->GetXaxis()->GetXmin();
    double x_high = data.histogram->GetXaxis()->GetXmax();
    TH1D* hFromJKTF1 = new TH1D(Form("hFromJKTF1_%d", i), "Histogram from JK TF1", nbins, x_low, x_high);

    for (int bin = 1; bin <= nbins; ++bin) {
      double binCenter = data.histogram->GetBinCenter(bin);
      double funcValue = mcfit->Eval(binCenter);
      hFromJKTF1->SetBinContent(bin, funcValue);
    }

    int transparentMagenta = TColor::GetColorTransparent(kMagenta, 0.3);
    hFromJKTF1->SetLineColor(kMagenta);
    hFromJKTF1->SetLineWidth(0);
    hFromJKTF1->SetFillColor(transparentMagenta);
    hFromJKTF1->SetFillStyle(1001);
    hFromJKTF1->Draw("SAME LF2");

    canvasJK->Update();
  }

  // Graphing results for Stat data with errors
  TCanvas* c_stat = new TCanvas("c_stat", "R_{sf} Stat", 1200, 400);
  std::vector<double> numbers_stat;
  std::vector<double> values_stat;
  std::vector<double> errors_stat;
  for (const auto& data : statData) {
    numbers_stat.push_back(data.histogramNumber);
    values_stat.push_back(data.R_sf);
    errors_stat.push_back(data.R_sf_err);
  }

  //Remove small stats last bin
  // Check if the vectors are not empty to avoid underflow
  if (!numbers_stat.empty() && !values_stat.empty() && !errors_stat.empty()) {
    numbers_stat.pop_back(); // Remove the last element from numbers_stat
    values_stat.pop_back();  // Remove the last element from values_stat
    errors_stat.pop_back();  // Remove the last element from errors_stat
  }

  //Get RMSE method stat error from root fits
  double RMSEValue = 0;
  double RMSEError = 0;

  if (!numbers_stat.empty()) {
    double xMin = *std::min_element(numbers_stat.begin(), numbers_stat.end());
    double xMax = *std::max_element(numbers_stat.begin(), numbers_stat.end());
    //double xMax = 11;
    
    TF1* RMSEfit = new TF1("RMSEfit", "pol0", xMin, xMax); // xMin and xMax should cover the range of your data.
    TGraphErrors* graph_stat = new TGraphErrors(numbers_stat.size(), numbers_stat.data(), values_stat.data(), nullptr, errors_stat.data());
    graph_stat->SetTitle("R_{sf} Stat;Run Range Index;R_{sf}");
    graph_stat->SetMarkerStyle(20);
    graph_stat->SetMarkerColor(kBlue);
    graph_stat->Draw("AP");
    graph_stat->GetHistogram()->SetMinimum(0.85);
    graph_stat->GetHistogram()->SetMaximum(1.15);
    RMSEfit->SetParameter(0,0.95); //initial guess
    graph_stat->Fit(RMSEfit, "RWQ");  // "W" uses the errors as weights, "Q" quiet mode.
    double RMSEfit_g = RMSEfit->GetParameter(0);
    RMSEfit->SetParameter(0,RMSEfit_g);
    graph_stat->Fit(RMSEfit, "RWQ");  // "W" uses the errors as weights, "Q" quiet mode.

    RMSEValue = RMSEfit->GetParameter(0);   // Gets the mean
    RMSEError = RMSEfit->GetParError(0);    // Gets the error
  }

  // Graphing results for Syst data with errors
  TCanvas* c_syst = new TCanvas("c_syst", "R_{sf} Syst", 1200, 400);
  std::vector<double> numbers_syst;
  std::vector<double> values_syst;
  std::vector<double> errors_syst;
  for (const auto& data : systData) {
    numbers_syst.push_back(data.histogramNumber);
    values_syst.push_back(data.R_sf);
    errors_syst.push_back(data.R_sf_err);
  }

  if (!numbers_syst.empty()) {
    TGraphErrors* graph_syst = new TGraphErrors(numbers_syst.size(), numbers_syst.data(), values_syst.data(), nullptr, errors_syst.data());
    graph_syst->SetTitle("R_{sf} Syst;Run Number;R_{sf}");
    graph_syst->SetMarkerStyle(21);
    graph_syst->SetMarkerColor(kMagenta);
    graph_syst->Draw("AP");
    graph_syst->GetHistogram()->SetMinimum(0.85);
    graph_syst->GetHistogram()->SetMaximum(1.15);
  }

  // Graphing results for Jackknife data with errors
  TCanvas* c_JK = new TCanvas("c_JK", "R_{sf} Jackknife", 1200, 400);
  std::vector<double> numbers_JK;
  std::vector<double> values_JK;
  std::vector<double> errors_JK;
  for (const auto& data : JKData) {
    numbers_JK.push_back(data.histogramNumber);
    values_JK.push_back(data.R_sf);
    errors_JK.push_back(data.R_sf_err);
  }

  if (!numbers_JK.empty()) {
    TGraphErrors* graph_JK = new TGraphErrors(numbers_JK.size(), numbers_JK.data(), values_JK.data(), nullptr, errors_JK.data());
    graph_JK->SetTitle("R_{sf} Jackknife;Subtracted Run Range Index;R_{sf}");
    graph_JK->SetMarkerStyle(22);
    graph_JK->SetMarkerColor(kMagenta);
    graph_JK->Draw("AP");
    graph_JK->GetHistogram()->SetMinimum(0.85);
    graph_JK->GetHistogram()->SetMaximum(1.15);
  }
  
  //Until we have infinite MC stats, leave this alone
  double mean = weightedMean(values_JK, errors_JK);
  double stdDev = weightedStdDev(values_JK, errors_JK, mean);

  double mean_syst = weightedMean(values_syst, errors_syst);
  double stdDev_syst = weightedStdDev(values_syst, errors_syst, mean_syst);

  //Jackknife uncertainty
  double meanJK = calculateMean(values_JK);
  double errJK = calculateJKErr(values_JK, meanJK, N_datasets);

  std::cout << "Mean of R_sf: " << mean << std::endl;
  std::cout << "StdDev of R_sf: " << stdDev << std::endl;

  //Use extract to get GMn and add error from model dependent extractions
  auto GMn_and_error = extract::extract_GMn_from_simc( RMSEValue,
						       RMSEError,
						       0.,              //Q2 will be calculated when 0.
						       E_beam,
						       BB_angle,
						       true);

  
  cout << "GMn/GD/mun = " << GMn_and_error[1].first << " +/- " << GMn_and_error[1].second << " (stat)" << endl;

  // Create a new canvas to compare the contributions
  TCanvas* JKCanvas = new TCanvas("JKCanvas", "Statistical Error With Jackknife", 800, 600);
  JKCanvas->Divide(1, 1); // Adjust based on how you want to layout the histograms/text
  JKCanvas->cd();

  // Prepare a TPaveText to display the jack knife results
  TPaveText* JKText = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
  JKText->SetTextAlign(12); // Left align
  JKText->SetBorderSize(1);
  JKText->SetFillStyle(0); // Transparent

  char rsfText[255];
  sprintf(rsfText, "Jackknife R_{sf} mean, std dev: %0.3f, %0.3f", meanJK, errJK);
  JKText->AddText(rsfText);

  char rmseText[255];
  sprintf(rmseText, "RMSE R_{sf} mean, error: %0.3f, %0.3f", RMSEValue, RMSEError);
  JKText->AddText(rmseText);

  char rtorText[255];
  sprintf(rtorText, "Run to Run R_{sf} mean, std dev: %0.3f, %0.3f", mean_syst, stdDev_syst);
  JKText->AddText(rtorText);

  char gmnText[255];
  sprintf(gmnText, "G_{M}^{n}/G_{D}/#mu_{n} = %.4f #pm %.4f (RMSE + WD_{fit})", GMn_and_error[1].first, GMn_and_error[1].second);
  JKText->AddText(gmnText);

  // Draw the JKText on the canvas
  JKText->Draw();
  JKCanvas->Update();
  JKCanvas->Write();

  //inputFile->Close();
  //inputFileMC->Close();
  outputFile->Write();
  //outputFile->Close();

  std::cout << "Analysis complete. Plots and data saved to " << foutPath << std::endl;
  
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
std::vector<std::pair<double, double>> fitAndFineFit(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, const std::string& fitOptions = "RBMQ0") {
  TF1* fit = new TF1(fitName.c_str(), fitFormula.c_str(), hcalfit_l, hcalfit_h, paramCount);
  //fit->SetNpx(5000);
  for (int i=0; i<paramCount; ++i){ //reset parameters/errors for this set
    fit->SetParameter(i,0);
    fit->SetParError(i,0);
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
  histogram->Fit(fineFit, fitOptions.c_str());

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

double weightedMean(const std::vector<double>& values, const std::vector<double>& errors) {
  double sumWeightedValues = 0;
  double sumWeights = 0;
    
  for (size_t i = 0; i < values.size(); ++i) {
    double weight = 1.0 / (errors[i] * errors[i]);
    sumWeights += weight;
    sumWeightedValues += weight * values[i];
  }
    
  return sumWeightedValues / sumWeights;
}

double weightedStdDev(const std::vector<double>& values, const std::vector<double>& errors, double weightedMean) {
  double sumWeightedVariances = 0;
  double sumWeights = 0;
    
  for (size_t i = 0; i < values.size(); ++i) {
    double weight = 1.0 / (errors[i] * errors[i]);
    sumWeights += weight;
    sumWeightedVariances += weight * pow(values[i] - weightedMean, 2);
  }
    
  return sqrt(sumWeightedVariances / sumWeights);
}

double calculateMean(const std::vector<double>& values) {
    double sum = 0.0;
    for (double value : values) {
        sum += value;
    }
    return sum / values.size();
}

double calculateStdDev(const std::vector<double>& values, double mean) {
    double sum = 0.0;
    for (double value : values) {
        sum += (value - mean) * (value - mean);
    }
    return sqrt(sum / values.size());
}

double calculateJKErr(const std::vector<double>& values, double mean, double num_datasets) {
    double sum = 0.0;
    for (double value : values) {
        sum += (value - mean) * (value - mean);
    }

    sum *= ((num_datasets-1)/num_datasets);

    return sqrt(sum / values.size());
}
