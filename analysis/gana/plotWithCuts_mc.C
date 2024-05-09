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
void plotWithCuts_mc(int kine=8, 
		     int mag=100, 
		     int pass=2, 
		     bool skipcorrelations=true,
		     bool addresscorr=false,
		     bool wide=false,
		     bool effz=true,
		     const char *replay_type = "alt" ) {

  // reading json configuration file
  JSONManager *jmgr = new JSONManager("../../config/syst.json");

  // Set up exception to MC file structure
  std::string rtype = replay_type;
  std::string rootfile_type = "";
  bool jboyd = false;
  bool alt = false;
  bool alt2 = false;
  if( rtype.compare("jboyd")==0 ){
    rootfile_type = "_jboyd";
    jboyd = true;
  }else if( rtype.compare("alt")==0 ){
    rootfile_type = "_alt";
    alt = true;
  }else if( rtype.compare("alt2")==0 ){
      rootfile_type = "_alt2";
      alt2 = true;
  }else if( rtype.compare("")!=0 )
    cout << "WARNING: invalid argument at replay_type. Valid entries include jboyd, alt, or nothing. Defaulting to nothing." << endl;

  //Get plot details
  int hbins = jmgr->GetValueFromSubKey<int>( "hbins", Form("sbs%d",kine) ); //gets to 7cm HCal pos res
  double hcalfit_l = jmgr->GetValueFromSubKey<int>( "hcalfit_l", Form("sbs%d",kine) ); //gets to 7cm HCal pos res
  double hcalfit_h = jmgr->GetValueFromSubKey<int>( "hcalfit_h", Form("sbs%d",kine) ); //gets to 7cm HCal pos res
  //double hcalfit_l = econst::hcalposXi_mc; //lower fit/bin limit for hcal dx plots (m)
  //double hcalfit_h = econst::hcalposXf_mc; //upper fit/bin limit for hcal dx plots (m)
  //int Nbins_dx = hbins;

  cout << "Nbins_dx:" << hbins << " hcal llim: " << hcalfit_l << " hcal ulim: " << hcalfit_h << endl;

  vector<int> plot_bins;
  vector<double> plot_lims;
  vector<double> llims;
  vector<double> ulims;

  jmgr->GetVectorFromSubKey<double>("plot_limits_mc",Form("sbs%d_%d",kine,mag),plot_lims);  
  jmgr->GetVectorFromSubKey<int>("plot_bins_mc",Form("sbs%d_%d",kine,mag),plot_bins);  
  
  for( size_t i=0; i<plot_lims.size(); ++i ){
    if( i%2==0 )
      llims.push_back(plot_lims[i]);
    else
      ulims.push_back(plot_lims[i]);
  }

  //Get wide/tight elastic cuts
  std::string globalcuts;
  if(wide){
    globalcuts = jmgr->GetValueFromSubKey_str( "post_cuts_mc", Form("sbs%d_%d",kine,mag) );
    cout << "Loaded wide cuts: " << globalcuts << endl;
  }else{
    globalcuts = jmgr->GetValueFromSubKey_str( "post_tcuts_mc", Form("sbs%d_%d",kine,mag) );
    cout << "Loaded tight cuts: " << globalcuts << endl;
  }

  std::vector<std::string> wcuts = util::parseCuts(globalcuts);
  std::vector<std::string> cuts;

  //Remove timing and tar cuts from globalcuts as MC is unreliable in this department and tar==LD2
  std::string parsedCutString;
  for (const auto& cut : wcuts) {
    std::string trimmedCut = trim(cut);
    if (trimmedCut.find("coin") == std::string::npos && trimmedCut.find("tar") == std::string::npos && trimmedCut.find("mag") == std::string::npos && trimmedCut.find("bb_grinch_tdc_clus_size") == std::string::npos) {
      if (!parsedCutString.empty()) parsedCutString += "&&";
      parsedCutString += trimmedCut;
      cuts.push_back(cut);
    }
  }

  cout << endl << endl << parsedCutString << endl << endl;

  cout << "Parsed cuts: " << endl;
  for( size_t i=0; i<cuts.size(); ++i ) 
    cout << cuts[i] << endl;

  //Make general cuts to extract proton and neutron from parsed data
  std::string protoncut = parsedCutString + "&&nucleon==0";
  std::string neutroncut = parsedCutString + "&&nucleon==1";

  //Get branches to examine systematics
  std::string postbranches = jmgr->GetValueFromKey_str( "post_branches_mc" );

  cout << "Loaded branches: " << postbranches << endl;

  std::vector<std::string> branches = util::parseCuts(postbranches);

  cout << "Parsed branches with plot limits: " << endl;
  for( size_t i=0; i<branches.size(); ++i ){ 
    cout << branches[i] << " llim=" << llims[i] << " ulim=" << ulims[i] << " bins=" << plot_bins[i] << endl;
  }

  cout << endl << endl;

  //Make a map between branches and plot details
  // Verify that all vectors are the same size
  if (!(branches.size() == llims.size() && branches.size() == ulims.size() && branches.size() == plot_bins.size())) {
    std::cerr << "ERROR: Branch vectors must be of the same size. Check json." << std::endl;
    cout << branches.size() << endl;
    cout << llims.size() << endl;
    cout << ulims.size() << endl;
    cout << plot_bins.size() << endl;


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
  std::string fin_path = outdir_path + Form("/parse/parse_mc_sbs%d_%dp_barebones%s%s.root",kine,mag,rootfile_type.c_str(),effz_word.c_str());
  std::string fin_path_bg = outdir_path + Form("/parse/parse_mcbg_sbs%d_%dp_barebones.root",kine,mag);

  std::string skipcorrelations_word = "";
  if(skipcorrelations)
    skipcorrelations_word = "_thin";

  std::string wide_word = "";
  if(wide)
    wide_word = "_widecut";

  std::string fout_path = outdir_path + Form("/gmn_analysis/dx_mc_sbs%d_mag%d_pass%d%s%s%s%s.root",kine,mag,pass,skipcorrelations_word.c_str(),rootfile_type.c_str(),wide_word.c_str(),effz_word.c_str());

  //Get MC inelastic background from file
  TFile* inputFileBG = new TFile(fin_path_bg.c_str());
  if (!inputFileBG || inputFileBG->IsZombie()) {
    std::cerr << "Error opening file: " << fin_path_bg << std::endl;
    return;
  }
  
  // Get the tree from the new input file
  TTree* treeBG = dynamic_cast<TTree*>(inputFileBG->Get("P"));
  if (!treeBG) {
    std::cerr << "Tree not found in file: " << fin_path_bg << std::endl;
    inputFileBG->Close();
    return;
  }
  
  // Clean elastic cuts
  // Using best cluster variables now
  // Open an output ROOT file
  TFile* outputFile = new TFile(fout_path.c_str(), "RECREATE");

  /////////////////////////////////////////
  //Draw the background histogram with all cuts
  std::string bghist = "hdx_inel";
  TH1D* hdx_inel = new TH1D(bghist.c_str(), "dx inelastic bg with cuts;m", hbins, hcalfit_l, hcalfit_h);

  // Draw the plot using the created histogram with the new cut
  treeBG->Draw(("dx>>" + bghist).c_str(), Form("mcsigma*mcomega*(%s)", parsedCutString.c_str()), "COLZ");
    
  // Write the histogram to the output file
  hdx_inel->Write();
    
  // Look at dx inelastic correlations with cut variables
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

    cout << "Working on inelastic histogram " << branch << ", bins " << bins << ", llim " << llim << ", ulim " << ulim << endl;

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

    std::string weightedCutString_bg = Form("mcsigma * mcomega * (%s)", cutString.c_str());

    cout << "Plotting inelastic dx versus branch " << branch << " with surviving uncorrelated weighted cuts: " << weightedCutString_bg << endl << endl;

    //Create the histogram
    std::string histTitle_bg = branch + " vs inelastic dx;" + branch + ";dx";
    std::string histName_bg = "hist_inel_" + branch;

    TH2D* hist_bg = new TH2D(histName_bg.c_str(), histTitle_bg.c_str(), bins, llim, ulim, hbins, hcalfit_l, hcalfit_h);

    // Draw the plot using the created histogram
    treeBG->Draw(("dx:" + branch + ">>" + histName_bg).c_str(), weightedCutString_bg.c_str(), "COLZ");

    // Write the histogram to the output file
    hist_bg->Write();

    // Clean up the histogram
    delete hist_bg;

  }

  // Cleanup
  delete hdx_inel;
  inputFileBG->Close();
  outputFile->Close();

  //Now get information from SIMC MC elastic data
  // Open the mc data root file
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

  std::cout << std::endl << "Successfully opened tree from file " << fin_path << std::endl;

  //reopen output file and update
  outputFile = new TFile(fout_path.c_str(), "UPDATE");

  /////////////////////////////////////////
  //Draw the proton histogram with all cuts
  std::string protonhist = "hdx_p";
  TH1D* hdx_p = new TH1D( protonhist.c_str(), "dx proton;m", hbins, hcalfit_l, hcalfit_h );

  // Draw the plot using the created histogram
  tree->Draw(("dx>>" + protonhist).c_str(), Form("mc_weight_norm * (%s)", protoncut.c_str()), "COLZ");

  // Write the histogram to the output file
  hdx_p->Write();

  delete hdx_p;

  //////////////////////////////////////////
  //Draw the neutron histogram with all cuts
  std::string neutronhist = "hdx_n";
  TH1D* hdx_n = new TH1D( neutronhist.c_str(), "dx neutron;m", hbins, hcalfit_l, hcalfit_h );

  // Draw the plot using the created histogram
  tree->Draw(("dx>>" + neutronhist).c_str(), Form("mc_weight_norm * (%s)", neutroncut.c_str()), "COLZ");

  // Write the histogram to the output file
  hdx_n->Write();

  delete hdx_n;

  //////////////////////////////////////
  //Draw overall histogram with all cuts
  std::string hist = "hdx";
  TH1D* hdx = new TH1D( hist.c_str(), "dx;m", hbins, hcalfit_l, hcalfit_h );

  // Draw the plot using the created histogram
  tree->Draw(("dx>>" + hist).c_str(), Form("mc_weight_norm * (%s)", parsedCutString.c_str()), "COLZ");

  // Write the histogram to the output file
  hdx->Write();

  delete hdx;

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
    if (branchIsCorrelated) {
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
    
    std::string protonCutString = cutString + "&&nucleon==0";
    std::string neutronCutString = cutString + "&&nucleon==1";

    std::string weightedCutString = Form("mc_weight_norm * (%s)", cutString.c_str());
    std::string weightedCutString_p = Form("mc_weight_norm * (%s)", protonCutString.c_str());
    std::string weightedCutString_n = Form("mc_weight_norm * (%s)", neutronCutString.c_str());

    cout << "Plotting dx versus branch " << branch << " with surviving uncorrelated weighted cuts: " << weightedCutString << endl << endl;

    //Create the histogram
    std::string histTitle = branch + " vs dx;" + branch + ";dx";
    std::string histTitle_p = branch + " vs dx (proton);" + branch + ";dx";
    std::string histTitle_n = branch + " vs dx (neutron);" + branch + ";dx";
    std::string histName = "hist_" + branch;
    std::string histName_p = "hist_" + branch + "_p";
    std::string histName_n = "hist_" + branch + "_n";

    TH2D* hist = new TH2D(histName.c_str(), histTitle.c_str(), bins, llim, ulim, hbins, hcalfit_l, hcalfit_h);
    TH2D* hist_p = new TH2D(histName_p.c_str(), histTitle_p.c_str(), bins, llim, ulim, hbins, hcalfit_l, hcalfit_h);
    TH2D* hist_n = new TH2D(histName_n.c_str(), histTitle_n.c_str(), bins, llim, ulim, hbins, hcalfit_l, hcalfit_h);

    // Draw the plot using the created histogram
    //tree->Draw(("dx:" + branch + ">>" + histName).c_str(), cutString.c_str(), "COLZ");
    tree->Draw(("dx:" + branch + ">>" + histName).c_str(), weightedCutString.c_str(), "COLZ");
    tree->Draw(("dx:" + branch + ">>" + histName_p).c_str(), weightedCutString_p.c_str(), "COLZ");
    tree->Draw(("dx:" + branch + ">>" + histName_n).c_str(), weightedCutString_n.c_str(), "COLZ");

    // Write the histogram to the output file
    hist->Write();
    hist_p->Write();
    hist_n->Write();

    // Clean up the histogram
    delete hist;
    delete hist_p;
    delete hist_n;

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

      std::string weightedCutString = Form("mc_weight_norm * (%s)", cutString.c_str());

      // Draw the plot using the created histogram
      tree->Draw((branchX + ":" + branchY + ">>" + histName).c_str(), weightedCutString.c_str(), "COLZ");

      // Write the histogram to the output file
      hist->Write();

      // Clean up the histogram
      delete hist;
    }
  }

  // // Convert map to a vector of strings for easy display
  // std::vector<std::string> cuts;
  // for (const auto& pair : correlatedCuts) {
  //   std::string cutInfo = pair.first + " correlated with " + pair.second;
  //   cuts.push_back(cutInfo);
  // }

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
