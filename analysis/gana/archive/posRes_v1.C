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

vector<int> kine = {4, 11, 14}; //only kinematics with sufficient LH2 elastics at zero field
vector<double> Q2 = {3.0, 13.6, 7.5}; // Q^2 values corresponding to kinematics

std::string globalcuts = "mag==0&&W2<1.2&&abs(coin)<6&&bb_ps_e>0.2&&bb_gem_track_nhits>3&&abs(bb_tr_vz)<0.072";

int dxbins = 200;
double dxllim = -1.0;
double dxulim = 1.0;

void fitGaussianWithLegend(TH1D* hist, TCanvas *can, double Q2Value) {
  // Find the bin with the maximum value
  int binMax = hist->GetMaximumBin();
  double maxVal = hist->GetBinContent(binMax);
    
  // Find the first bins at half maximum to the left and right of the maximum bin
  double halfMax = maxVal / 3.0;
    
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
  can->cd();
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
  legend->AddEntry((TObject*)0, Form("Mean = %.4f", mean), "");
  legend->AddEntry((TObject*)0, Form("Sigma = %.4f", sigma), "");
  legend->AddEntry((TObject*)0, Form("Q^{2} = %.1f", Q2Value), "");
  legend->Draw("same");

  can->Update();
}

void posRes(bool skipcorrelations=true, bool addresscorr=false, bool effz=true, bool wide=false) {
  gStyle->SetOptStat(0110);
  gStyle->SetStatW(0.2); // Set width of the stats box
  gStyle->SetStatH(0.05); // Set height of the stats box
  gStyle->SetStatX(0.8); // Set X position (0 to 1, where 1 is the right edge of the canvas)
  gStyle->SetStatY(0.8); // Set Y position (0 to 1, where 1 is the top edge of the canvas)

  //set up files and paths
  std::string effz_word = "";
  if(effz)
    effz_word = "_effz";

  std::string outdir_path = gSystem->Getenv("OUT_DIR");

  for(int i = 0; i < kine.size(); ++i) {
    std::string fin_path = outdir_path + Form("/parse/parse_sbs%d_pass2_barebones%s_lh2.root",kine[i],effz_word.c_str());
    cout << "Setting up input path: " << fin_path << endl;
    std::string fout_path = outdir_path + Form("/gmn_analysis/posres/posres_out_sbs%d.root", kine[i]);
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

    TCanvas *c1 = new TCanvas(Form("c1_sbs%d",kine[i]),Form("dx res sbs%d",kine[i]),800,800);
    c1->Divide(2, 2);

    for(int j = 1; j <= 4; ++j) {
      c1->cd(j);
      std::string hdxName = "hdx_cuts_" + std::to_string(kine[i]);
      std::string hdxCuts = globalcuts + "&&abs(dy)<0.3";
      TH1D* hdx = new TH1D(hdxName.c_str(), "dx;m", dxbins, dxllim, dxulim);
      tree->Draw(("dx>>" + hdxName).c_str(), hdxCuts.c_str(), "COLZ");
      fitGaussianWithLegend(hdx, c1, Q2[i]);
    }

    TCanvas *c2 = new TCanvas(Form("c2_sbs%d",kine[i]),Form("dy res sbs%d",kine[i]),800,800);
    c2->Divide(2, 2);

    for(int j = 1; j <= 4; ++j) {
      c2->cd(j);
      std::string hdyName = "hdy_cuts_" + std::to_string(kine[i]);
      std::string hdyCuts = globalcuts + "&&abs(dx)<0.3";
      TH1D* hdy = new TH1D(hdyName.c_str(), "dy;m", dxbins, dxllim, dxulim);
      tree->Draw(("dy>>" + hdyName).c_str(), hdyCuts.c_str(), "COLZ");
      fitGaussianWithLegend(hdy, c2, Q2[i]);
    }

    c1->Update();
    c1->Write();

    c2->Update();
    c2->Write();
  }
  cout << "All plots created. Output files are located in: " << outdir_path << endl;
}
