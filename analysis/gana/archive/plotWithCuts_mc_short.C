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
void plotWithCuts_mc_short(int kine=4, int pass=2, int mag=50, bool bestclus=false ) {

  // reading json configuration file
  JSONManager *jmgr = new JSONManager("../../config/syst.json");

  //Get plot details and binning
  int hbins = jmgr->GetValueFromSubKey<int>( "hbins", Form("sbs%d",kine) ); //gets to 7cm HCal pos res
  double hcalfit_l = econst::hcalposXi_mc; //lower fit/bin limit for hcal dx plots (m)
  double hcalfit_h = econst::hcalposXf_mc; //upper fit/bin limit for hcal dx plots (m)

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

  //Get tight elastic cuts
  std::string globalcuts = jmgr->GetValueFromSubKey_str( Form("post_cuts_p%d",pass), Form("sbs%d_%d",kine,mag) );

  cout << "Loaded tight cuts: " << globalcuts << endl;

  std::vector<std::string> cuts = util::parseCuts(globalcuts);

  //Remove timing and tar cuts from globalcuts as MC is unreliable in this department and tar==LD2
  std::string cutString;
  for (const auto& cut : cuts) {
    std::string trimmedCut = trim(cut);
    if (trimmedCut.find("coin") == std::string::npos && trimmedCut.find("tar") == std::string::npos && trimmedCut.find("mag") == std::string::npos) {
      if (!cutString.empty()) cutString += "&&";
      cutString += trimmedCut;
    }
  }

  cout << "Parsed cuts: " << endl;
  for( size_t i=0; i<cuts.size(); ++i ) 
    cout << cuts[i] << endl;

  //Make general cuts to extract proton and neutron from parsed data
  std::string protoncut = cutString + "&&nucleon==0";
  std::string neutroncut = cutString + "&&nucleon==1";

  //set up files and paths
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  //std::string fin_path = outdir_path + Form("/parse/parse_sbs%d_pass%d.root",kine,pass);
  std::string fin_path = outdir_path + Form("/parse/parse_mc_sbs%d_%dp_barebones.root",kine,mag);

  std::string bestclus_word = "";

  //Would need to change dx to dx_bc to implement
  if(bestclus)
    bestclus_word = "_bestclus";

  std::string fout_path = outdir_path + Form("/gmn_analysis/dx_mc_sbs%d_mag%d_pass%d%s_rd2_short.root",kine,mag,pass,bestclus_word.c_str());

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
  tree->Draw(("dx>>" + hist).c_str(), Form("mc_weight_norm * (%s)", cutString.c_str()), "COLZ");

  // Write the histogram to the output file
  hdx->Write();

  delete hdx;

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
