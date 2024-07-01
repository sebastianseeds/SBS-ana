//sseeds - 1.24.24: Script to plot data MC comparison on HCal energy
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>

//MAIN (no args)
void plotHCalEnergy() {

  gStyle->SetOptStat(0);

  // Open the ROOT files
  TFile *file1 = TFile::Open("/lustre19/expphy/volatile/halla/sbs/seeds/parse/parse_mc_sbs9_70p_barebones.root");
  TFile *file2 = TFile::Open("/lustre19/expphy/volatile/halla/sbs/seeds/parse/parse_sbs9_pass2_barebones.root");

  if (!file1 || !file2) {
    std::cerr << "Error opening files. Exiting." << std::endl;
    return;
  }

  // Make sure the Tree exists in both files
  TTree *tree1 = nullptr;
  TTree *tree2 = nullptr;
  file1->GetObject("P", tree1);
  file2->GetObject("P", tree2);

  if (!tree1 || !tree2) {
    std::cerr << "Tree not found. Exiting." << std::endl;
    return;
  }

  // Create histograms
  TH1D *h1 = new TH1D("h1", "Elastic MC proton", 200, 0, 1);
  TH1D *h2 = new TH1D("h2", "Elastic Data LH2", 200, 0, 1);

  // Draw with conditions
  tree1->Draw("hcale>>h1", "hcale>0&&nucleon==0&&W2<1.2&&bb_ps_e>0.2&&abs(dy)<0.8", "goff");
  tree2->Draw("hcale>>h2", "W2<1.2&&bb_ps_e>0.2&&abs(dy)<0.8&&tar==0", "goff");

  // Calculate the scale factor and scale the data histogram
  double scaleFactor = h1->Integral() / h2->Integral();
  h2->Scale(scaleFactor);

  // Create a canvas
  TCanvas *c1 = new TCanvas("c1", "HCal Energy SBS-9", 800, 600);
  c1->cd();

  // Draw histograms
  h1->SetLineColor(kRed);
  h1->SetTitle("HCal Energy SBS-9;GeV;");
  h1->Draw("HIST");

  h2->SetLineColor(kBlue);
  h2->Draw("HIST SAME");

  // Add a legend
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(h1, "Elastic MC proton", "l");
  legend->AddEntry(h2, "Elastic Data LH2", "l");
  legend->Draw();

  // Update the canvas
  c1->Update();

}
