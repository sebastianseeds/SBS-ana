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

std::string globalcuts = "mag==0&&W2<1.2&&bb_ps_e>0.2&&bb_gem_track_nhits>3&&abs(bb_tr_vz)<0.072";
std::string tglobalcuts = "bb_tr_chi2<10&&bb_gem_track_nhits>3&&abs(bb_tr_vz)<0.07&&mag==0&&W2<1.2&&bb_ps_e>0.2&&bb_gem_track_nhits>3&&abs(bb_tr_vz)<0.072&&(bb_tr_r_x-0.9*bb_tr_r_th)<0.25&&(bb_tr_r_x-0.9*bb_tr_r_th)>-0.3&&(bb_tr_r_y-0.9*bb_tr_r_ph)<0.1&&(bb_tr_r_y-0.9*bb_tr_r_ph)>-0.1";
const int idselect = 125; //select one channel with onechan

std::string spotcut = "&&(pow(dx + 0.04, 2) / pow(0.26, 2)) + (pow(dy + 0.05, 2) / pow(0.33, 2)) <= 1"; //assumes zero field!!

const vector<int> midselect = {123,124,125,126,127}; //select several close channels

const int tdcbins = 45;
const double tdcllim = -15.0;
const double tdculim = 10.0;

const int tofoffset = 18; //offset to return tof corrected tdc peak to zero

// Main function
void tdcTimeRes(bool tglob=false, bool onechan=false, bool mchan=false, bool spots=true) {
  gStyle->SetOptStat(0110);
  gStyle->SetStatW(0.2); // Set width of the stats box
  gStyle->SetStatH(0.05); // Set height of the stats box
  gStyle->SetStatX(0.35); // Set X position (0 to 1, where 1 is the right edge of the canvas)
  gStyle->SetStatY(0.8); // Set Y position (0 to 1, where 1 is the top edge of the canvas)

  if(onechan&&mchan){
    cout << "ERROR: Select one of onechan and mchan." << endl;
    return;
  }

  std::string cuts;
  if(tglob)
    cuts = tglobalcuts;
  else
    cuts = globalcuts;

  if(onechan)
    cuts += Form("&&hcalpid==%d",idselect);

  if(mchan){
    cuts += "&&(";
    for(size_t i=0; i<midselect.size(); ++i){
      if(i==0)
	cuts += Form("hcalpid==%d",midselect[i]);
      else
	cuts += Form("||hcalpid==%d",midselect[i]);
    }
    cuts += ")";
  }
    
  if(spots)
    cuts += spotcut;

  cout << endl << "All plots with cuts " << cuts << endl << endl;

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
  TCanvas *c1 = new TCanvas("c1", "TDC res", 2400, 1600);
  c1->Divide(3, 2);

  std::vector<std::string> variables = {"hcaltdc", 
					"hcaltdc-bb_sh_atime", 
					"hcaltdctw", 
					Form("hcaltdctof-%d",tofoffset), 
					Form("hcaltdccorr-%d",tofoffset), 
					Form("hcaltdccorr-bb_sh_atime-%d",tofoffset)};

  std::vector<std::string> titles = {"TDC (no correction);ns", 
				     "TDC (trigger corrected);ns", 
				     "TDC (timewalk corrected);ns", 
				     "TDC (p_{N} corrected);ns", 
				     "TDC (timewalk and p_{N} corrected);ns", 
				     "TDC (all corrections);ns"};

  double baseline_res = 1.0;

  for (int i = 0; i < variables.size(); ++i) {
    c1->cd(i+1);
    std::string histName = "h" + variables[i];
    TH1D* hist = new TH1D(histName.c_str(), titles[i].c_str(), tdcbins, tdcllim, tdculim);
    tree->Draw((variables[i] + ">>" + histName).c_str(), cuts.c_str(), "COLZ");

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

    //Get percent improvement
    double improve = (baseline_res-sigma)/baseline_res*100.;

    // Create a legend and add the mean and standard deviation
    TLegend* legend = new TLegend(0.15, 0.5, 0.4, 0.65);
    legend->SetBorderSize(0); // Remove border
    legend->SetFillStyle(0);  // Make background transparent
    legend->SetMargin(0);     // Remove left margin
    legend->SetTextColor(kRed);
    legend->AddEntry((TObject*)0, Form("Mean = %.4f ns", mean), "");
    legend->AddEntry((TObject*)0, Form("Sigma = %.4f ns", sigma), "");
    if(i>0)
      legend->AddEntry((TObject*)0, Form("Improvement = %0.2f%%", improve), "");
    if(onechan)
      legend->AddEntry((TObject*)0, Form("Channel = %d", idselect), "");
    legend->Draw("same");

    if(i==0)
      baseline_res = sigma;

    c1->Update();
  }

  c1->Update();
  c1->Write();

  cout << "All plots created. Output file located here: " << fout_path << endl;
}
