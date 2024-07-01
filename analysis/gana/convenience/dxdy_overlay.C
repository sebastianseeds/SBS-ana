//sseeds 4.1.24: Plot overlay script for visualization
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
#include "../../../src/jsonmgr.C"
#include "../../../include/gmn.h"

//MAIN. kine=kinematic, mag=SBS magnetic field, pass=pass
void dxdy_overlay(int kine=8, int mag=70, int pass=2 ) {

  // reading json configuration file
  JSONManager *jmgr = new JSONManager("../../../config/syst.json");

  //Get tight elastic cuts
  std::string globalcuts = jmgr->GetValueFromSubKey_str( Form("post_cuts_p%d",pass), Form("sbs%d_%d",kine,mag) );

  std::cout << "Loaded tight cuts: " << globalcuts << std::endl;

  std::vector<std::string> cuts = util::parseCuts(globalcuts); //this makes a vector of all individual cuts in the single globalcut string.

  std::cout << "Parsed cuts: " << std::endl;

  // Remove any cuts that contain "dy"
  std::string cut_elas;
  for (const auto& cut : cuts) {
    if (cut.find("dy") == std::string::npos) {
      cut_elas += cut + ";";
    }
  }

  std::cout << "Parsed cuts without 'dy': " << cut_elas << std::endl;

  //set up files and paths
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fin_path = outdir_path + Form("/parse/parse_sbs%d_pass%d_barebones.root",kine,pass);

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

  std::string histBaseName = "hdxdy";
  TH2D* hist_base = new TH2D( histBaseName.c_str(), "dxdy;dy;dx", 200, -1, 1,300,-2,1);

  // Draw the plot using the created histogram
  tree->Draw(("dx:dy>>" + histBaseName).c_str(), cut_elas.c_str(), "colz");

  // Create x-axis projection
  TH1D* hist1d_x = hist_base->ProjectionX();

  // Create y-axis projection
  TH1D* hist1d_y = hist_base->ProjectionY();

  // Create a canvas
  TCanvas* canvas = new TCanvas("canvas", "Overlay TH1Ds on TH2D", 800, 600);

  // Draw the 2D histogram
  hist_base->Draw("colz");

  // Draw the x-axis projection on the right side of the 2D histogram
  hist1d_x->SetLineColor(kRed);
  hist1d_x->Draw("same");

  // Draw the y-axis projection on top of the 2D histogram
  hist1d_y->SetLineColor(kBlue);
  hist1d_y->Draw("same");

  // Draw legend
  TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(hist1d_x, "X Projection", "l");
  legend->AddEntry(hist1d_y, "Y Projection", "l");
  legend->Draw();

  // Save the canvas as an image
  std::string outputFileName = "output.png";
  canvas->SaveAs(outputFileName.c_str());

  // Keep the canvas open
  canvas->Draw();
}
