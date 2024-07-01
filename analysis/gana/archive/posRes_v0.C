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

vector<int> kine = {4,11,14,8}; //only kinematics with sufficient LH2 elastics at zero field

std::string globalcuts = "mag==0&&W2<1.2&&abs(coin)<6&&bb_ps_e>0.2&&bb_gem_track_nhits>3&&abs(bb_tr_vz)<0.072";

int dxbins = 200;
double dxllim = -1.0;
double dxulim = 1.0;

void fitGaussianWithLegend(TH1D* hist, TCanvas *can) {
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
  can->cd();
  //TCanvas* canvas = new TCanvas("canvas", "Gaussian Fit", 800, 600);
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
  legend->Draw("same");

  can->Update();
}

//main
void posRes(bool skipcorrelations=true, 
	    bool addresscorr=false,
	    bool effz=true,
	    bool wide=false) {

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
  std::string fin_path = outdir_path + Form("/parse/parse_sbs%d_pass2_barebones%s_lh2.root",kine[0],effz_word.c_str());

  cout << "Setting up input path: " << fin_path << endl;

  std::string fout_path = outdir_path + Form("/gmn_analysis/posres/posres_out.root");

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

  // Clean elastic cuts
  // Using best cluster variables now
  // Open an output ROOT file
  TFile* outputFile = new TFile(fout_path.c_str(), "RECREATE");

  TCanvas *c1 = new TCanvas("c1",Form("dx res sbs%d",kine[0]),800,800);
  c1->cd();

  std::string hdxName = "hdx_cuts";
  std::string hdxCuts = globalcuts + "&&abs(dy)<0.3";
  TH1D* hdx = new TH1D( hdxName.c_str(), "dx;m", dxbins, dxllim, dxulim );

  // Draw the plot using the created histogram
  tree->Draw(("dx>>" + hdxName).c_str(), hdxCuts.c_str(), "COLZ");

  // Fit
  fitGaussianWithLegend(hdx,c1);

  c1->Update();
  c1->Write();

  cout << "All plots created. Output file located here: " << fout_path << endl;

}
