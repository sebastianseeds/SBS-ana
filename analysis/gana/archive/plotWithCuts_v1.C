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

//Known correlations. Elastic selection on dx:<branch> will ignore cuts on correlated variables.
std::map<std::string, std::string> correlatedCuts = {
  {"bb.sh.nclus", "bb.sh.nclus&&sbs.hcal.nclus"},
  {"sbs.hcal.nclus", "sbs.hcal.nclus&&bb.sh.nclus"},
  {"dy", "dy&&W2"},
  {"W2", "W2&&dy"}
};

// trim function for accurate parsing
std::string trim(const std::string& str) {
  size_t first = str.find_first_not_of(" \t\n\r\f\v");
  if (std::string::npos == first) {
    return str;
  }
  size_t last = str.find_last_not_of(" \t\n\r\f\v");
  return str.substr(first, (last - first + 1));
}

//plot with cuts script, first step on systematics analysis. Configured for pass 2.
//bestclus uses bestclusterinfo from parse file. skipcorrelations doesn't plot cut vs cut plots. addresscorr removes cuts from dx vs cut TH2Ds where the cut is correlated (like dy and W2)
void plotWithCuts(int kine=8, 
		  int mag=70, 
		  int pass=2, 
		  bool skipcorrelations=true, 
		  bool addresscorr=false,
		  bool effz=true,
		  bool wide=false) {

  // reading json configuration file
  JSONManager *jmgr = new JSONManager("../../config/syst.json");

  //Get plot details
  int hbins = jmgr->GetValueFromSubKey<int>( "hbins", Form("sbs%d",kine) ); //gets to 7cm HCal pos res
  double hcalfit_l = jmgr->GetValueFromSubKey<int>( "hcalfit_l", Form("sbs%d",kine) ); //gets to 7cm HCal pos res
  double hcalfit_h = jmgr->GetValueFromSubKey<int>( "hcalfit_h", Form("sbs%d",kine) ); //gets to 7cm HCal pos res
  int W2bins = jmgr->GetValueFromSubKey<int>( "W2bins", Form("sbs%d",kine) );
  double W2fit_l = jmgr->GetValueFromSubKey<int>( "W2fit_l", Form("sbs%d_%d",kine,mag) );
  double W2fit_h = jmgr->GetValueFromSubKey<int>( "W2fit_h", Form("sbs%d_%d",kine,mag) );

  //int Nbins_dx = hbins;

  cout << "Nbins_dx:" << hbins << " hcal llim: " << hcalfit_l << " hcal ulim: " << hcalfit_h << endl;

  vector<int> plot_bins;
  vector<double> plot_lims;
  vector<double> llims;
  vector<double> ulims;

  jmgr->GetVectorFromSubKey<double>(Form("plot_limits_p%d",pass),Form("sbs%d_%d",kine,mag),plot_lims);  
  jmgr->GetVectorFromSubKey<int>(Form("plot_bins_p%d",pass),Form("sbs%d_%d",kine,mag),plot_bins);  
  
  for( size_t i=0; i<plot_lims.size(); ++i ){
    if( i%2==0 )
      llims.push_back(plot_lims[i]);
    else
      ulims.push_back(plot_lims[i]);
  }

  //Get wide/tight elastic cuts
  std::string globalcuts;
  if(wide){
    globalcuts = jmgr->GetValueFromSubKey_str( Form("post_cuts_p%d",pass), Form("sbs%d_%d",kine,mag) );
    cout << "Loaded wide cuts: " << globalcuts << endl;
  }else{
    globalcuts = jmgr->GetValueFromSubKey_str( Form("post_tcuts_p%d",pass), Form("sbs%d_%d",kine,mag) );
    cout << "Loaded tight cuts: " << globalcuts << endl;
  }

  cout << "Loaded tight cuts: " << globalcuts << endl;

  //Get nucleon spot cuts for W2 histograms
  std::string spotcut_n = jmgr->GetValueFromSubKey_str( "h_dn_spot_cuts_tight", Form("sbs%d_%d",kine,mag) );
  std::string spotcut_p = jmgr->GetValueFromSubKey_str( "h_dp_spot_cuts_tight", Form("sbs%d_%d",kine,mag) );

  std::string spotcut = "&&(" + spotcut_n + ")||(" + spotcut_p + ")";

  std::vector<std::string> cuts = util::parseCuts(globalcuts); //this makes a vector of all individual cuts in the single globalcut string.
  std::vector<std::string> ccuts;

  std::cout << "Parsed cuts: " << std::endl;

  std::string parsedCutString;
  std::string parsedCutString_W2;
  std::string parsedCutString_W2dy;

  for (const auto& cut : cuts) {
    std::string trimmedCut = trim(cut);

    if (!parsedCutString.empty()) parsedCutString += "&&";
    parsedCutString += trimmedCut;

    // Include cut in parsedCutString_W2 if it doesn't contain "W2"
    if (trimmedCut.find("W2") == std::string::npos) {
      if (!parsedCutString_W2.empty()) parsedCutString_W2 += "&&";
      parsedCutString_W2 += trimmedCut;

      // Include cut in parsedCutString_W2dy if it doesn't contain "W2" or "dy"
      if (trimmedCut.find("dy") == std::string::npos) {
	if (!parsedCutString_W2dy.empty()) parsedCutString_W2dy += "&&";
	parsedCutString_W2dy += trimmedCut;
      }
    }

    ccuts.push_back(cut);
  }

  cout << endl << endl << parsedCutString_W2 << endl << endl;
  cout << endl << endl << parsedCutString_W2dy << endl << endl;

  //Set up dy anticut bg histogram
  std::string dyanticut_noelas;
  std::string dyanticut_elas;

  for (size_t i = 0; i < cuts.size(); ++i) {
    std::cout << cuts[i] << std::endl;
    
    // Check if the cut contains "W2"
    size_t found_W2 = cuts[i].find("W2");
    if (found_W2 != std::string::npos) {
      // If "W2" is found, skip adding this cut to dyanticut_elas
      continue; // Skip to the next iteration of the loop
    }

    size_t found_dy = cuts[i].find("dy");
    if (found_dy != std::string::npos) {
      // Found "dy", replace the first '<' with '>'
      size_t found_lt = cuts[i].find('<', found_dy);
      if (found_lt != std::string::npos) {
	std::string modified_cut = cuts[i];
	modified_cut[found_lt] = '>';
	dyanticut_noelas = modified_cut; // Overwrite with this modified cut
	// Append modified cut to dyanticut_elas, but exclude "W2" cuts
	dyanticut_elas += (dyanticut_elas.empty() ? "" : "&&") + modified_cut;
      }
    } else {
      // "dy" not found, just append the cut to dyanticut_elas if it doesn't contain "W2"
      dyanticut_elas += (dyanticut_elas.empty() ? "" : "&&") + cuts[i];
    }
  }

  std::cout << std::endl << "dyanticut_noelas: " << dyanticut_noelas << std::endl;
  std::cout << "dyanticut_elas: " << dyanticut_elas << std::endl;

  //Set up coin anticut bg histogram
  std::string coinanticut_noelas;
  std::string coinanticut_elas;
  for (size_t i = 0; i < cuts.size(); ++i) {

    size_t found_coin = cuts[i].find("coin");
    if (found_coin != std::string::npos) {
      // Found "coin", replace the first '<' with '>'
      size_t found_lt = cuts[i].find('<', found_coin);
      if (found_lt != std::string::npos) {
	std::string modified_cut = cuts[i];
	modified_cut[found_lt] = '>';
	coinanticut_noelas = modified_cut; // Overwrite with this modified cut
	// Append modified cut to coinanticut_elas, but exclude "W2" cuts
	coinanticut_elas += (coinanticut_elas.empty() ? "" : "&&") + modified_cut;
      }
    } else {
      // "coin" not found, just append the cut to coinanticut_elas
      coinanticut_elas += (coinanticut_elas.empty() ? "" : "&&") + cuts[i];
    }
  }

  std::cout << "coinanticut_noelas: " << coinanticut_noelas << std::endl;
  std::cout << "coinanticut_elas: " << coinanticut_elas << std::endl;

  //return;

  std::string postbranches = jmgr->GetValueFromKey_str( Form("post_branches_p%d",pass) );

  cout << "Loaded branches: " << postbranches << endl;

  std::vector<std::string> branches = util::parseCuts(postbranches);

  cout << "Total branches: " << branches.size() << endl;

  cout << "Parsed branches with plot limits: " << endl;
  for( size_t i=0; i<branches.size(); ++i ){ 
    cout << branches[i] << " llim=" << llims[i] << " ulim=" << ulims[i] << " bins=" << plot_bins[i] << endl;
  }

  //Make a map between branches and plot details
  // Verify that all vectors are the same size
  if (!(branches.size() == llims.size() && branches.size() == ulims.size() && branches.size() == plot_bins.size())) {
    std::cerr << "ERROR: Branch vectors must be of the same size. Check json." << std::endl;
    return;
  }

  std::map<std::string, std::tuple<double, double, int>> branchmap;

  for(size_t i = 0; i < branches.size(); ++i) {
    branchmap[branches[i]] = std::make_tuple(llims[i], ulims[i], plot_bins[i]);
  }

  //set up files and paths
  std::string effz_word = "";
  if(effz)
    effz_word = "_effz";

  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fin_path = outdir_path + Form("/parse/parse_sbs%d_pass%d_barebones%s_ld2.root",kine,pass,effz_word.c_str());

  std::string skipcorrelations_word = "";
  if(skipcorrelations)
    skipcorrelations_word = "_thin";

  std::string wide_word = "";
  if(wide)
    wide_word = "_widecut";

  std::string fout_path = outdir_path + Form("/gmn_analysis/dx_correlations_sbs%d_mag%d_pass%d%s%s%s.root",kine,mag,pass,skipcorrelations_word.c_str(),wide_word.c_str(),effz_word.c_str());

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

  std::string histBaseName = "hdx_allcut";
  TH1D* hist_base = new TH1D( histBaseName.c_str(), "dx;m", hbins, hcalfit_l, hcalfit_h );

  // Draw the plot using the created histogram
  tree->Draw(("dx>>" + histBaseName).c_str(), globalcuts.c_str(), "COLZ");

  // Write the histogram to the output file
  hist_base->Write();

  delete hist_base;

  std::string histAntidyName = "hdx_dyanti";
  TH1D* hist_antidy = new TH1D( histAntidyName.c_str(), "dx, dy anticut;m", hbins, hcalfit_l, hcalfit_h );

  // Draw the plot using the created histogram
  tree->Draw(("dx>>" + histAntidyName).c_str(), dyanticut_elas.c_str(), "COLZ");

  // Write the histogram to the output file
  hist_antidy->Write();

  delete hist_antidy;

  std::string histAnticoinName = "hdx_coinanti";
  TH1D* hist_anticoin = new TH1D( histAnticoinName.c_str(), "dx, coin anticut;m", hbins, hcalfit_l, hcalfit_h );

  // Draw the plot using the created histogram
  tree->Draw(("dx>>" + histAnticoinName).c_str(), coinanticut_elas.c_str(), "COLZ");

  // Write the histogram to the output file
  hist_anticoin->Write();

  delete hist_anticoin;

  /////////////////////////////////////////
  // Draw the W2 histogram with all cuts
  std::string hW2hist = "hW2";
  std::string hW2_cut = parsedCutString_W2 + spotcut;
  TH1D* hW2 = new TH1D(hW2hist.c_str(), "W2;GeV^{2}", W2bins, W2fit_l, W2fit_h);

  // Draw the plot using the created histogram with the W2-excluded cut
  tree->Draw(("W2>>" + hW2hist).c_str(), Form("mc_weight_norm * (%s)", hW2_cut.c_str()), "COLZ");

  // Write the histogram to the output file
  hW2->Write();

  /////////////////////////////////////////
  // Draw the W2 histogram excluding dy cuts
  std::string hW2dyhist = "hW2_nodycut";
  std::string hW2dy_cut = parsedCutString_W2dy + spotcut;
  TH1D* hW2dy = new TH1D(hW2dyhist.c_str(), "W2 without dy cuts;GeV^{2}", W2bins, W2fit_l, W2fit_h);

  // Draw the plot using the created histogram with the W2 and dy-excluded cut
  tree->Draw(("W2>>" + hW2dyhist).c_str(), Form("mc_weight_norm * (%s)", hW2dy_cut.c_str()), "COLZ");

  // Write the histogram to the output file
  hW2dy->Write();

  // Loop over branches
  for (const auto& branch : branches) {

    //get plot details from branch 
    double llim;
    double ulim;
    int bins;
    if(branchmap.find(branch) != branchmap.end()) {
      auto& tuple = branchmap[branch];
      llim = std::get<0>(tuple);
      ulim = std::get<1>(tuple);
      bins = std::get<2>(tuple);

    } else {
      std::cout << "ERROR: Branch key not found: " << branch << std::endl;
    }

    cout << "Working on histogram " << branch << ", bins " << bins << ", llim " << llim << ", ulim " << ulim << endl;

    std::string cutString;
    bool branchIsCorrelated = correlatedCuts.find(branch) != correlatedCuts.end();

    // Process correlated branches
    if (branchIsCorrelated && addresscorr) {
      std::cout << "Branch is correlated on " << branch << std::endl;
      std::string correlatedCutString = correlatedCuts[branch];
      std::vector<std::string> correlatedCutsVector = util::parseCuts(correlatedCutString);

      // Loop through all cuts and exclude the ones that are correlated
      for (const auto& cut : cuts) {
	std::string trimmedCut = trim(cut);
	bool cutFound = false;
	for (const auto& correlatedCut : correlatedCutsVector) {
	  std::string trimmedCorrelatedCut = trim(correlatedCut);
	  if (trimmedCut.find(trimmedCorrelatedCut) != std::string::npos) {
	    cutFound = true;
	    break;
	  }
	}
	if (!cutFound) {
	  if (!cutString.empty()) cutString += "&&";
	  cutString += trimmedCut;
	}
      }
    }else {
      for (const auto& cut : cuts) {
	std::string trimmedCut = trim(cut);
	if (trimmedCut.find(branch) == std::string::npos) {
	  if (!cutString.empty()) cutString += "&&";
	  cutString += trimmedCut;
	}
      }
    }
    
    cout << "Plotting dx versus branch " << branch << " with surviving uncorrelated cuts: " << cutString << endl << endl;

    //Create the histogram
    //std::string histTitle = branch + " vs dx;" + branch + ";dx";
    std::string histTitle = cutString + ";" + branch + ";dx";
    std::string histName = "hist_" + branch;

    TH2D* hist = new TH2D(histName.c_str(), histTitle.c_str(), bins, llim, ulim, hbins, hcalfit_l, hcalfit_h);

    // Draw the plot using the created histogram
    tree->Draw(("dx:" + branch + ">>" + histName).c_str(), cutString.c_str(), "COLZ");

    // Write the histogram to the output file
    hist->Write();

    // Clean up the histogram
    delete hist;

  }
  
  if(skipcorrelations){
    cout << "All plots created. Correlation plots skipped. Output file located here: " << fout_path << endl;
    return;
  }

  // Loop over all pairs of branches to check correlations
  for (size_t i = 0; i < branches.size(); ++i) {

    //get outer branch plot details
    std::string branchX = branches[i];
    double Xllim;
    double Xulim;
    int Xbins;
    if(branchmap.find(branchX) != branchmap.end()) {
      auto& tuple = branchmap[branchX];
      Xllim = std::get<0>(tuple);
      Xulim = std::get<1>(tuple);
      Xbins = std::get<2>(tuple);

    } else {
      std::cout << "ERROR: Branch key not found: " << branchX << std::endl;
    }
    
    for (size_t j = i + 1; j < branches.size(); ++j) { //avoid reverse pairs

      //get inner branch plot details
      std::string branchY = branches[j];
      double Yllim;
      double Yulim;
      int Ybins;
      if(branchmap.find(branchY) != branchmap.end()) {
	auto& tuple = branchmap[branchY];
	Yllim = std::get<0>(tuple);
	Yulim = std::get<1>(tuple);
	Ybins = std::get<2>(tuple);

      } else {
	std::cout << "ERROR: Branch key not found: " << branchY << std::endl;
      }

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
      TH2D* hist = new TH2D(histName.c_str(), histTitle.c_str(), Xbins, Xllim, Xulim, Ybins, Yllim, Yulim);

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
