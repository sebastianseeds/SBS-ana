//sseeds 
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

const int tdcbins = 400;
const double tdcllim = -20.0;
const double tdculim = 20.0;

// Main
void tdcTimeRes() {
  gStyle->SetOptStat(0110);
  gStyle->SetStatW(0.2); // Set width of the stats box
  gStyle->SetStatH(0.05); // Set height of the stats box
  gStyle->SetStatX(0.8); // Set X position (0 to 1, where 1 is the right edge of the canvas)
  gStyle->SetStatY(0.8); // Set Y position (0 to 1, where 1 is the top edge of the canvas)

  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fin_path = outdir_path + "/parse/parse_sbs8_pass2_barebones_effz_lh2.root";
  std::string fout_path = outdir_path + "/gmn_analysis/tdcres/tdcres_out.root";
  cout << "Setting up input path: " << fin_path << endl;
  cout << "Setting up output path: " << fout_path << endl;

  // Open the ROOT file
  TFile* inputFile = new TFile(fin_path.c_str());
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Error opening file: " << fin_path << std::endl;
    return;
  }

  // Get the tree
  TTree* tree = dynamic_cast<TTree*>(inputFile->Get("P"));
  if (!tree) {
    std::cerr << "Tree not found in file: " << fin_path << std::endl;
    inputFile->Close();
    return;
  }

  // Open an output ROOT file
  TFile* outputFile = new TFile(fout_path.c_str(), "RECREATE");

  // Canvas for TDC plots
  TCanvas *c1 = new TCanvas("c1", "TDC res", 1600, 1600);
  c1->Divide(2, 2);

  std::vector<std::string> variables = {"hcaltdc", "hcaltdc-bb_sh_atime", "hcaltdctw", "hcaltdctof"};
  std::vector<std::string> titles = {"hcaltdc;ns", "hcaltdc - bb_sh_atime;ns", "hcaltdctw;ns", "hcaltdctof;ns"};

  for (int i = 0; i < variables.size(); ++i) {
    c1->cd(i+1);
    std::string histName = "h" + variables[i];
    TH1D* hist = new TH1D(histName.c_str(), titles[i].c_str(), tdcbins, tdcllim, tdculim);
    tree->Draw((variables[i] + ">>" + histName).c_str(), "variable_name != 0", "COLZ");

    // Find the bin with the maximum value
    int binMax = hist->GetMaximumBin();
    double maxVal = hist->GetBinContent(binMax);

    // Find the first bins at half maximum to the left and right of the maximum bin
    double halfMax = maxVal / 2.0;

    int binLeft = binMax;
    while (binLeft > 1 && hist->GetBinContent(binLeft) > halfMax) {
      binLeft--;
    }

    int binRight = binMax;
    while (binRight < hist->GetNbinsX() && hist->GetBinContent(binRight) > halfMax) {
      binRight++;
    }

    // Fit the distribution with a Gaussian
    TF1* gaussFit = new TF1("gaussFit", "gaus", hist->GetBinCenter(binLeft), hist->GetBinCenter(binRight));
    hist->Fit(gaussFit, "R");

    // Get the mean and standard deviation from the fit
    double mean = gaussFit->GetParameter(1);
    double sigma = gaussFit->GetParameter(2);

    // Draw the histogram and fit
    hist->SetLineColor(kBlack);
    hist->SetLineWidth(2);
    hist->Draw();
    gaussFit->Draw("same");

    // Create a legend and add the mean and standard deviation
    TLegend* legend = new TLegend(0.6, 0.5, 0.85, 0.65);
    legend->SetBorderSize(0); // Remove border
    legend->SetFillStyle(0);  // Make background transparent
    legend->SetMargin(0);     // Remove left margin
    legend->SetTextColor(kRed);
    legend->AddEntry((TObject*)0, Form("Mean = %.4f ns", mean), "");
    legend->AddEntry((TObject*)0, Form("Sigma = %.4f ns", sigma), "");
    legend->Draw("same");

    c1->Update();
  }

  c1->Update();
  c1->Write();

  cout << "All plots created. Output file located here: " << fout_path << endl;
}
