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

//this configured for pass gmn 0/1
void plotWithCuts(int kine, int mag, int pass=1 ) {

  //set up files and paths
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  //std::string fin_path = outdir_path + Form("/parse/parse_sbs%d_pass%d.root",kine,pass);
  std::string fin_path = outdir_path + Form("/parse/parse_sbs%d_pass%d_update.root",kine,pass);

  std::string fout_path = outdir_path + Form("/gmn_analysis/dx_correlations_sbs%d_mag%d_pass%d.root",kine,mag,pass);

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

  // Branches and corresponding cuts
  std::vector<std::string> branches = {"sbs.hcal.e", 
				       "bb.tr.vz", 
				       "dy", 
				       "bb.ps.e",
				       "W2",
				       "bb.gem.track.nhits",
				       "bb.sh.nclus",
				       "sbs.hcal.nclus"};

  // Elastic cuts
  std::vector<std::string> cuts = {"sbs.hcal.e>0.03", 
				   "abs(bb.tr.vz)<0.075",  
				   "abs(dy)<0.25", 
				   "bb.ps.e>0.2", 
				   "abs(W2-0.88)<0.3",
				   "bb.gem.track.nhits>3",
				   "sbs.hcal.nclus>0",
				   "bb.sh.nclus>0",
				   "tar==1", 
				   Form("mag==%d", mag)};

  // Clean elastic cuts

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
