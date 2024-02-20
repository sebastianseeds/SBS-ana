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

//plot with cuts script, first step on systematics analysis. Configured for pass 2.
void plotWithCuts_p2(int kine, int mag, int pass, bool bestclus ) {

  // reading json configuration file
  JSONManager *jmgr = new JSONManager("../../config/syst.json");

  std::string globalcuts = jmgr->GetValueFromSubKey_str( Form("post_cuts_p%d",pass), Form("sbs%d_%d",kine,mag) );

  cout << "Loaded tight cuts: " << globalcuts << endl;

  std::vector<std::string> cuts = util::parseCuts(globalcuts);

  cout << "Parsed cuts: " << endl;
  for( size_t i=0; i<cuts.size(); ++i ) 
    cout << cuts[i] << endl;

  std::string postbranches = jmgr->GetValueFromKey_str( Form("post_branches_p%d",pass) );

  cout << "Loaded branches: " << postbranches << endl;

  std::vector<std::string> branches = util::parseCuts(postbranches);

  cout << "Parsed branches: " << endl;
  for( size_t i=0; i<branches.size(); ++i ) 
    cout << branches[i] << endl;

  //set up files and paths
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  //std::string fin_path = outdir_path + Form("/parse/parse_sbs%d_pass%d.root",kine,pass);
  std::string fin_path = outdir_path + Form("/parse/parse_sbs%d_pass%d_barebones_rd2.root",kine,pass);

  std::string bestclus_word = "";

  if(bestclus)
    bestclus_word = "_bestclus";

  std::string fout_path = outdir_path + Form("/gmn_analysis/dx_correlations_sbs%d_mag%d_pass%d%s_rd2.root",kine,mag,pass,bestclus_word.c_str());

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

  // Branches to examine
  // std::vector<std::string> branches = {"hcale",
  // 				       "bb_tr_vz", 
  // 				       "dy", 
  // 				       "bb_ps_e",
  // 				       "W2",
  // 				       "bb_gem_track_nhits",
  // 				       "bb_etot_over_p",
  // 				       "coin",
  // 				       "hcalnblk",
  // 				       "bb_sh_rowblk",
  // 				       "bb_sh_colblk",
  // 				       "hcalx",
  // 				       "hcaly",
  // 				       "fiducial_sig_x",
  // 				       "fiducial_sig_y"};


  // Clean elastic cuts
  // Using best cluster variables now
  // Open an output ROOT file
  TFile* outputFile = new TFile(fout_path.c_str(), "RECREATE");

  // Define fixed binning for x and y axes
  int nBins_dx = 300; // Number of bins in x
  int nBins_cut = 300; // Number of bins in y

  // Loop over branches
  for (const auto& branch : branches) {

    std::string cutString;
    for (const auto& cut : cuts) {
      if (cut.find(branch) == std::string::npos) {
	if (!cutString.empty()) cutString += " && ";
	cutString += cut;
      }
    }

    std::string histTitle = branch + " vs dx;" + branch + ";dx";
    cout << "Working on histogram " << branch << endl;
    std::string histName = "hist_" + branch;

    cout << histName << endl;

    TH2D* hist = new TH2D(histName.c_str(), histTitle.c_str(), nBins_dx, 0, 0, nBins_cut, 0, 0);

    // Draw the plot using the created histogram
    tree->Draw(("dx:" + branch + ">>" + histName).c_str(), cutString.c_str(), "COLZ");

    // Write the histogram to the output file
    hist->Write();

    // Clean up the histogram
    delete hist;

  }

  // Loop over all pairs of branches to check correlations
  for (size_t i = 0; i < branches.size(); ++i) {
    for (size_t j = i + 1; j < branches.size(); ++j) { //avoid reverse pairs
      std::string branchX = branches[i];
      std::string branchY = branches[j];

      // Construct cut string excluding cuts involving either of the pair branches
      std::string cutString;
      for (const auto& cut : cuts) {
	if (cut.find(branchX) == std::string::npos && cut.find(branchY) == std::string::npos) {
	  if (!cutString.empty()) cutString += " && ";
	  cutString += cut;
	}
      }

      std::string histName = "coorhist_" + branchX + "_vs_" + branchY;
      std::string histTitle = branchX + " vs " + branchY + ";" + branchY + ";" + branchX;
      cout << "Working on histogram " << histName << endl;
      TH2D* hist = new TH2D(histName.c_str(), histTitle.c_str(), nBins_dx, 0, 0, nBins_cut, 0, 0);

      // Draw the plot using the created histogram
      tree->Draw((branchX + ":" + branchY + ">>" + histName).c_str(), cutString.c_str(), "COLZ");

      // Write the histogram to the output file
      hist->Write();

      // Clean up the histogram
      delete hist;
    }
  }

  // Create a canvas for displaying cuts
  TCanvas* cutCanvas = new TCanvas("cutCanvas", "Cuts", 800, 600);
  TPaveText* cutsText = new TPaveText(0.1, 0.1, 0.9, 0.9, "NB NDC"); // coordinates are normalized

  cutsText->AddText("Cuts Used:");
  cutsText->AddLine();
  for (const auto& cut : cuts) {
    cutsText->AddText(cut.c_str());
  }

  cutsText->SetTextAlign(12); // Align text to left
  cutsText->SetFillColor(0);  // Transparent background
  cutsText->Draw();

  // Write the canvas to the output file
  cutCanvas->Write();

  // Close the files
  inputFile->Close();
  delete inputFile;

  outputFile->Close();
  delete outputFile;

  cout << "All plots created. Output file located here: " << fout_path << endl;

}
