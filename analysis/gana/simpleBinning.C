//sseeds - 5.15.24: Plots branch variable comparisons with several algorithmically generated bin ranges
#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TStyle.h>
#include <TCanvas.h>

// Function to calculate the standard deviation
double calculateStdDev(const std::vector<double>& data) {
  double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
  double sq_sum = std::inner_product(data.begin(), data.end(), data.begin(), 0.0);
  double stddev = std::sqrt(sq_sum / data.size() - mean * mean);
  return stddev;
}

// Function to calculate the interquartile range (IQR)
double calculateIQR(std::vector<double> data) {
  std::sort(data.begin(), data.end());
  size_t q1_index = data.size() / 4;
  size_t q3_index = 3 * data.size() / 4;
  double Q1 = data[q1_index];
  double Q3 = data[q3_index];
  return Q3 - Q1;
}

// Function to calculate the number of bins using Scott's rule
int calculateBinsScottsRule(const std::vector<double>& data) {
  double stddev = calculateStdDev(data);
  double binWidth = 3.5 * stddev / std::cbrt(data.size());
  double range = *max_element(data.begin(), data.end()) - *min_element(data.begin(), data.end());
  return std::ceil(range / binWidth);
}

// Function to calculate the number of bins using the Freedman-Diaconis rule
int calculateBinsFreedmanDiaconis(const std::vector<double>& data) {
  double IQR = calculateIQR(data);
  double binWidth = 2.0 * IQR / std::cbrt(data.size());
  double range = *max_element(data.begin(), data.end()) - *min_element(data.begin(), data.end());
  return std::ceil(range / binWidth);
}

// Function to calculate the number of bins using Sturges' rule
int calculateBinsSturges(int n) {
  return std::ceil(std::log2(n) + 1);
}

//MAIN. kine=kinematic, mag=magnetic field setting (percent), pass=reconstruction pass, effz=use plot file with effective z offset, isld2=use plot file with LD2 data, branchname=plot branch name, addcuts=all cut string, cut_l=plot lower limit, cut_h=plot upper limit
void simpleBinning(int kine=8, 
		   int mag=70, 
		   int pass=2, 
		   bool effz=true,
		   bool isld2=true,
		   const std::string& branchname = "dx", 
		   const std::string& addcuts = "&&W2<1.2&&abs(coin)<10&&abs(dy)<0.3", 
		   double cut_l = -2.0, 
		   double cut_h = 1.0) {
  // Set draw params
  gStyle->SetPalette(55);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetLabelSize(0.03, "XYZ");
  gStyle->SetTitleSize(0.04, "XYZ");
  gStyle->SetEndErrorSize(0);

  //set up files and paths
  std::string effz_word = "";
  if(effz)
    effz_word = "_effz";

  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string filename = outdir_path + Form("/parse/parse_sbs%d_pass%d_barebones%s_ld2.root",kine,pass,effz_word.c_str());

  // Open the ROOT file
  TFile* inputFile = new TFile(filename.c_str());
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }

  // Get the tree
  TTree* tree = dynamic_cast<TTree*>(inputFile->Get("P"));
  if (!tree) {
    std::cerr << "Tree not found in file: " << filename << std::endl;
    inputFile->Close();
    return;
  }

  std::string cuts = Form("mag==%d&&tar==%d",mag,isld2) + addcuts;

  // Create a vector to store the data
  std::vector<double> data;

  // Create a histogram with 500 bins to get the data
  TH1D* tempHist = new TH1D("tempHist", "Temporary Histogram;X;Y", 500, cut_l, cut_h);
  tree->Draw((branchname + ">>tempHist").c_str(), cuts.c_str(), "goff");

  // Extract data from the histogram to the vector
  for (int i = 1; i <= tempHist->GetNbinsX(); ++i) {
    double binContent = tempHist->GetBinContent(i);
    double binCenter = tempHist->GetBinCenter(i);
    for (int j = 0; j < binContent; ++j) {
      data.push_back(binCenter);
    }
  }

  // Calculate the number of bins using different rules
  int binsSturges = calculateBinsSturges(data.size());
  int binsScott = calculateBinsScottsRule(data);
  int binsFreedmanDiaconis = calculateBinsFreedmanDiaconis(data);

  // Print the binning recommendations
  std::cout << "Number of bins according to Sturges' rule: " << binsSturges << std::endl;
  std::cout << "Number of bins according to Scott's rule: " << binsScott << std::endl;
  std::cout << "Number of bins according to the Freedman-Diaconis rule: " << binsFreedmanDiaconis << std::endl;

  // Create a canvas with four partitions
  TCanvas* canvas = new TCanvas("canvas", "Histogram Binning", 1200, 800);
  canvas->Divide(2, 2);

  // Create the common title for the histograms
  std::string histTitle = branchname + ", cuts: " + cuts;

  // Draw the histogram with 500 bins
  canvas->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  TH1D* hist500 = new TH1D("hist500", (histTitle + ";" + branchname + ";").c_str(), 500, cut_l, cut_h);
  tree->Draw((branchname + ">>hist500").c_str(), cuts.c_str(), "goff");
  hist500->SetLineWidth(2);
  hist500->SetLineColor(kBlack);
  hist500->Draw();

  // Create and style the legend
  TLegend* legend500 = new TLegend(0.5, 0.7, 0.9, 0.8);
  legend500->SetTextColor(kRed);
  legend500->SetBorderSize(0);  // No border
  legend500->SetFillStyle(0);   // Transparent fill
  legend500->SetTextSize(0.04);
  legend500->SetMargin(0.06);   // Remove left margin
  legend500->AddEntry((TObject*)0, "Binning Method: Fixed", "");
  legend500->AddEntry((TObject*)0, "Total Bins: 500", "");
  legend500->Draw();

  // Draw the histogram with bins calculated using Sturges' rule
  canvas->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  TH1D* histSturges = new TH1D("histSturges", (histTitle + ";" + branchname + ";").c_str(), binsSturges, cut_l, cut_h);
  tree->Draw((branchname + ">>histSturges").c_str(), cuts.c_str(), "goff");
  histSturges->SetLineWidth(2);
  histSturges->SetLineColor(kBlack);
  histSturges->Draw();

  // Create and style the legend
  TLegend* legendSturges = new TLegend(0.5, 0.7, 0.9, 0.8);
  legendSturges->SetTextColor(kRed);
  legendSturges->SetBorderSize(0);  // No border
  legendSturges->SetFillStyle(0);   // Transparent fill
  legendSturges->SetTextSize(0.04);
  legendSturges->SetMargin(0.06);   // Remove left margin
  legendSturges->AddEntry((TObject*)0, "Binning Method: Sturges' rule", "");
  legendSturges->AddEntry((TObject*)0, Form("Total Bins: %d", binsSturges), "");
  legendSturges->Draw();

  // Draw the histogram with bins calculated using Scott's rule
  canvas->cd(3);
  gPad->SetGridx();
  gPad->SetGridy();
  TH1D* histScott = new TH1D("histScott", (histTitle + ";" + branchname + ";").c_str(), binsScott, cut_l, cut_h);
  tree->Draw((branchname + ">>histScott").c_str(), cuts.c_str(), "goff");
  histScott->SetLineWidth(2);
  histScott->SetLineColor(kBlack);
  histScott->Draw();

  // Create and style the legend
  TLegend* legendScott = new TLegend(0.5, 0.7, 0.9, 0.8);
  legendScott->SetTextColor(kRed);
  legendScott->SetBorderSize(0);  // No border
  legendScott->SetFillStyle(0);   // Transparent fill
  legendScott->SetTextSize(0.04);
  legendScott->SetMargin(0.06);   // Remove left margin
  legendScott->AddEntry((TObject*)0, "Binning Method: Scott's rule", "");
  legendScott->AddEntry((TObject*)0, Form("Total Bins: %d", binsScott), "");
  legendScott->Draw();

  // Draw the histogram with bins calculated using the Freedman-Diaconis rule
  canvas->cd(4);
  gPad->SetGridx();
  gPad->SetGridy();
  TH1D* histFreedmanDiaconis = new TH1D("histFreedmanDiaconis", (histTitle + ";" + branchname + ";").c_str(), binsFreedmanDiaconis, cut_l, cut_h);
  tree->Draw((branchname + ">>histFreedmanDiaconis").c_str(), cuts.c_str(), "goff");
  histFreedmanDiaconis->SetLineWidth(2);
  histFreedmanDiaconis->SetLineColor(kBlack);
  histFreedmanDiaconis->Draw();

  // Create and style the legend
  TLegend* legendFreedmanDiaconis = new TLegend(0.5, 0.7, 0.9, 0.8);
  legendFreedmanDiaconis->SetTextColor(kRed);
  legendFreedmanDiaconis->SetBorderSize(0);  // No border
  legendFreedmanDiaconis->SetFillStyle(0);   // Transparent fill
  legendFreedmanDiaconis->SetTextSize(0.04);
  legendFreedmanDiaconis->SetMargin(0.06);   // Remove left margin
  legendFreedmanDiaconis->AddEntry((TObject*)0, "Binning Method: Freedman-Diaconis rule", "");
  legendFreedmanDiaconis->AddEntry((TObject*)0, Form("Total Bins: %d", binsFreedmanDiaconis), "");
  legendFreedmanDiaconis->Draw();

}
