//seeds
//Gets central Q2 value from MC with wide elastic and fiducial cuts. The Median is used as the central value from the non-symmetric Q2 distribution.
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>

std::string gcut = "mc_weight_norm*(W2>0.5&&W2<1.2&&abs(coin)<6&&bb_ps_e>0.2&&hcale>0.05&&fiducial_sig_x>2.0&&fiducial_sig_y>1.0)";

vector<std::string> gcut_vect = {"W^{2}>0.5","W^{2}<1.2","abs(BBCal - HCal coincidence time)<6","BBCal PS>0.2","HCal E>0.05","fiducial cut x safety margin = 2.0#sigma","fiducial cut y safety margin = 1.0#sigma"};

int Q2bins = 500;
vector<double> Q2llim = {0.0,0.0,7.0,10.5,4.5,1.5,1.5,1.5,1.5};
vector<double> Q2ulim = {6.0,6.0,13.0,16.5,10.5,7.5,7.5,7.5,7.5};

// Function to compute the median directly from the histogram
double computeMedianFromHistogram(TH1D* hist) {
  int nBins = hist->GetNbinsX();
  double totalEntries = hist->Integral(hist->FindBin(0.5), hist->FindBin(20.0)); //avoid over/underflow bins
  double cumulativeEntries = 0;
  for (int i = 1; i <= nBins; ++i) {
    cumulativeEntries += hist->GetBinContent(i);
    if (cumulativeEntries >= totalEntries / 2.0) {
      return hist->GetBinCenter(i);
    }
  }
  return 0; // Fallback, should not happen
}

// Function to compute the mode directly from the histogram
double computeModeFromHistogram(TH1D* hist) {
  int maxBin = hist->GetMaximumBin();
  return hist->GetBinCenter(maxBin);
}

// Main function
void getMCQ2() {

  //set draw params
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetPalette(55);
  gStyle->SetOptStat(0110);
  gStyle->SetStatTextColor(kBlack);
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.08);
  gStyle->SetStatX(0.35);
  gStyle->SetStatY(0.9);

  // Directory with parse files
  std::string outdir_path = gSystem->Getenv("OUT_DIR");

  // File paths and titles
  std::vector<std::string> filePaths = {
    outdir_path +"/parse/parse_mc_sbs4_30p_barebones_alt_effz.root",
    outdir_path +"/parse/parse_mc_sbs4_50p_barebones_alt_effz.root",
    outdir_path +"/parse/parse_mc_sbs7_85p_barebones_effz.root",
    outdir_path +"/parse/parse_mc_sbs11_100p_barebones_effz.root",
    outdir_path +"/parse/parse_mc_sbs14_70p_barebones_effz.root",
    outdir_path +"/parse/parse_mc_sbs8_50p_barebones_alt_effz.root",
    outdir_path +"/parse/parse_mc_sbs8_70p_barebones_alt_effz.root",
    outdir_path +"/parse/parse_mc_sbs8_100p_barebones_alt_effz.root",
    outdir_path +"/parse/parse_mc_sbs9_70p_barebones_alt_effz.root"
  };

  std::vector<int> kineidx = {
    4,
    4,
    7,
    11,
    14,
    8,
    8,
    8,
    9
  };

  std::vector<int> fieldidx = {
    30,
    50,
    85,
    100,
    70,
    50,
    70,
    100,
    70
  };
  // Create a canvas to draw the histograms
  TCanvas *c = new TCanvas("c", "Q^2 Distribution for Various Kinematics", 1800, 1500);
  c->Divide(3, 3);
  c->cd(1);

  std::vector<double> means, medians, modes;

  for (size_t i = 0; i < filePaths.size(); ++i) {

    cout << "Working on kinematic " << kineidx[i] << endl;

    // Open the ROOT file
    TFile* inputFile = new TFile(filePaths[i].c_str());
    if (!inputFile || inputFile->IsZombie()) {
      std::cerr << "Error opening file: " << filePaths[i] << std::endl;
      return;
    }

    // Get the tree
    TTree* tree = dynamic_cast<TTree*>(inputFile->Get("P"));
    if (!tree) {
      std::cerr << "Tree not found in file: " << filePaths[i] << std::endl;
      inputFile->Close();
      return;
    }

    cout << "Got parse tree P from file " << filePaths[i] << endl;

    // Create a canvas to draw the histograms
    c->cd(i+1);

    gPad->SetGridx();
    gPad->SetGridy();

    // Create a histogram to hold the Q^2 values
    std::string hQ2Name = Form("Q2_hist_%d",kineidx[i]);
    std::string hQ2Cut = gcut;
    std::string hQ2Title = Form("Q^{2} (SBS-%d); GeV^{2}",kineidx[i]);
    TH1D* histQ2 = new TH1D(hQ2Name.c_str(), 
			    hQ2Title.c_str(), 
			    Q2bins, 
			    Q2llim[i], 
			    Q2ulim[i]); // Adjust binning and range as needed

    // Fill the histogram with Q^2 values from the tree
    tree->Draw(("Q2>>" + hQ2Name).c_str(), hQ2Cut.c_str(), "HIST");

    // Compute the mean, median, and mode Q^2
    double mean_Q2 = histQ2->GetMean();
    double median_Q2 = computeMedianFromHistogram(histQ2);
    double mode_Q2 = computeModeFromHistogram(histQ2);

    means.push_back(mean_Q2);
    medians.push_back(median_Q2);
    modes.push_back(mode_Q2);

    //Style the histogram
    histQ2->SetLineColor(kBlack);
    histQ2->SetLineWidth(2);
    histQ2->SetFillColor(kBlue - 10);  // Light blue fill color
    histQ2->SetFillStyle(3004);        // Hatch fill style

    // Draw the histogram on the canvas
    histQ2->Draw("HIST");

    // Create and style the legend
    TLegend *legend = new TLegend(0.1, 0.6, 0.4, 0.7);
    legend->SetTextColor(kRed);
    legend->SetBorderSize(0);  // No border
    legend->SetFillStyle(0);   // Transparent fill
    legend->SetTextSize(0.04);
    legend->SetMargin(0.06);   // Remove left margin

    legend->AddEntry((TObject*)0, Form("Mean Q^{2}: %.2f", mean_Q2), "");
    legend->AddEntry((TObject*)0, Form("Median Q^{2}: %.2f", median_Q2), "");
    legend->AddEntry((TObject*)0, Form("Mode Q^{2}: %.2f", mode_Q2), "");
    legend->Draw();

    // Update the canvas
    c->Update();

  }

  c->Update();

  c->Write();

  // Create a canvas for the itemized cuts from gcut_vect
  TCanvas *cCuts = new TCanvas("cCuts", "Itemized Cuts", 1200, 800);
  cCuts->cd();
  TLatex latex;
  latex.SetTextSize(0.04);
  latex.DrawLatexNDC(0.1, 0.9, "Itemized Cuts:");

  for (size_t j = 0; j < gcut_vect.size(); ++j) {
    latex.DrawLatexNDC(0.1, 0.9 - 0.07 * (j + 1), gcut_vect[j].c_str());
  }

  latex.DrawLatexNDC(0.1, 0.9 - 0.07 * (gcut_vect.size() + 1), " ");
  latex.DrawLatexNDC(0.1, 0.9 - 0.07 * (gcut_vect.size() + 2), "MC weights are applied to each histogram.");

  cCuts->Update();
  cCuts->Write();

  // Write LaTeX table to console
  std::cout << "\\begin{table}[ht]" << std::endl;
  std::cout << "\\centering" << std::endl;
  std::cout << "\\begin{tabular}{|c|c|c|c|c|c|}" << std::endl;
  std::cout << "\\hline" << std::endl;
  std::cout << "\\textbf{Kinematic} & \\textbf{SBS Field} & \\textbf{Mean} & \\textbf{Median} & \\textbf{Mode} & \\textbf{Central (Median)} \\\\" << std::endl;
  std::cout << "\\hline" << std::endl;

  for (size_t k = 0; k < kineidx.size(); ++k) {
    std::cout << kineidx[k] << " & " << fieldidx[k] << " & " << means[k] << " & " << medians[k] << " & " << modes[k] << " & \\textbf{" << medians[k] << "} \\\\" << std::endl;
    std::cout << "\\hline" << std::endl;
  }

  std::cout << "\\end{tabular}" << std::endl;
  std::cout << "\\caption{$Q^2$ distribution statistics for different kinematic configurations}" << std::endl;
  std::cout << "\\label{tab:q2_statistics}" << std::endl;
  std::cout << "\\end{table}" << std::endl;


}
