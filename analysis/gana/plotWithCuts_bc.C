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
#include <regex>
#include "TStopwatch.h"
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

//Known correlations. Elastic selection on dx:<branch> will ignore cuts on correlated variables.
std::map<std::string, std::string> correlatedCuts = {
  {"dy_bc", "dy_bc&&W2"},
  {"W2", "W2&&dy_bc"}
};

std::string wide_globalcut = "hcale>0.015&&abs(bb_tr_vz)<0.08&&abs(dy)<0.60&&bb_ps_e>0.10&&abs(W2-0.8)<0.8&&bb_gem_track_nhits>2&&abs(bb_etot_over_p-1.0)<0.4&&abs(coin)<17.0&&hcalon==1&&tar==1&&fiducial_sig_x>1&&fiducial_sig_y>1";

std::vector<std::string> excludeCuts = {"bb_grinch_tdc_clus_size","bb_grinch_tdc_clus_trackindex","pspot","nspot"};

//forward declarations
std::string trim(const std::string& str);
std::string addbc(std::string input);
std::vector<std::string> split(const std::string &s, char delimiter);
std::map<std::string, std::string> getRowContents(const std::string &filePath, int kine, int mag, const std::string &target, const std::vector<std::string> &excludeKeys);

//plot with cuts script to work with best clusters, first step on systematics analysis. Configured for pass 2.
//bestclus uses bestclusterinfo from parse file. skipcorrelationplots doesn't plot cut vs cut plots. addresscorr removes cuts from dx vs cut TH2Ds where the cut is correlated (like dy and W2)
void plotWithCuts_bc(int kine=9, 
		     int mag=70, 
		     int pass=2, 
		     bool skipcorrelationplots=true, 
		     bool addresscorr=false,
		     bool wide=false,
		     bool effz=true) {

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // reading json configuration file
  JSONManager *jmgr = new JSONManager("../../config/syst.json");

  //Get plot details
  int hbins = jmgr->GetValueFromSubKey<int>( "hbins", Form("sbs%d",kine) ); //gets to 7cm HCal pos res
  int hbins_dy = jmgr->GetValueFromSubKey<int>( "hbins_dy", Form("sbs%d",kine) );
  double hcalfit_l = jmgr->GetValueFromSubKey<double>( "hcalfit_l_wide", Form("sbs%d_%d",kine,mag) ); //gets to 7cm HCal pos res
  double hcalfit_h = jmgr->GetValueFromSubKey<double>( "hcalfit_h_wide", Form("sbs%d_%d",kine,mag) ); //gets to 7cm HCal pos res
  double hcalfit_l_dy = jmgr->GetValueFromSubKey<double>( "hcalfit_l_dy_wide", Form("sbs%d_%d",kine,mag) ); //gets to 7cm HCal pos res
  double hcalfit_h_dy = jmgr->GetValueFromSubKey<double>( "hcalfit_h_dy_wide", Form("sbs%d_%d",kine,mag) ); //gets to 7cm HCal pos res

  int W2bins = jmgr->GetValueFromSubKey<int>( "W2bins", Form("sbs%d",kine) );
  double W2fit_l = jmgr->GetValueFromSubKey<double>( "W2fit_l_wide", Form("sbs%d_%d",kine,mag) );
  double W2fit_h = jmgr->GetValueFromSubKey<double>( "W2fit_h_wide", Form("sbs%d_%d",kine,mag) );

  cout << "Running best cluster dx plots.." << endl << endl;

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

  // std::string globalcuts_raw;
  // if(wide){
  //   globalcuts_raw = jmgr->GetValueFromSubKey_str( Form("post_cuts_p%d",pass), Form("sbs%d_%d",kine,mag) );
  //   cout << "Loaded wide cuts: ";
  // }else{
  //   globalcuts_raw = jmgr->GetValueFromSubKey_str( Form("post_tcuts_p%d",pass), Form("sbs%d_%d",kine,mag) );
  //   cout << "Loaded tight cuts: ";
  // }

  //std::string globalcuts = addbc(globalcuts_raw);

  //cout << globalcuts << endl;

  //get new cuts from .csv
  std::string cutsheet_path = "/w/halla-scshelf2102/sbs/seeds/ana/data/p2_cutset.csv";
  std::string target = "ld2";
  std::vector<std::string> excludeKeys = {};

  // Get the row contents
  std::map<std::string, std::string> rowContents = getRowContents(cutsheet_path, kine, mag, target, excludeCuts);
  // Print the row contents
  cout << endl << "Loading NEW cuts..." << endl;
  for (const auto &content : rowContents) {
    std::cout << content.first << ": " << content.second << std::endl;
  }
  cout << endl;

  std::string globalcuts_raw;
  for (const auto &content : rowContents) {
    if (!content.second.empty()) {
      if (!globalcuts_raw.empty()) {
	globalcuts_raw += "&&";
      }
      globalcuts_raw += content.second;
    }
  }

  cout << endl <<"Concatenated globalcuts_raw: " << globalcuts_raw << endl << endl;

  std::string globalcuts = addbc(globalcuts_raw);

  cout << endl <<"Concatenated globalcuts: " << globalcuts << endl << endl;

  if(wide){
    globalcuts = wide_globalcut;
    cout << "Wide cut override. New cut string: " << wide_globalcut << endl;
  }

  //Get nucleon spot cuts for W2 histograms
  std::string spotcut_n = jmgr->GetValueFromSubKey_str( "h_dn_spot_cuts_tight", Form("sbs%d_%d",kine,mag) );
  std::string spotcut_p = jmgr->GetValueFromSubKey_str( "h_dp_spot_cuts_tight", Form("sbs%d_%d",kine,mag) );

  std::string spotcutraw = "&&((" + spotcut_n + ")||(" + spotcut_p + "))";

  std::string spotcut = addbc(spotcutraw);

  cout << "Loaded spotcut: " << spotcut << endl;

  std::string spotanticut_n = jmgr->GetValueFromSubKey_str( "h_dn_spot_anticuts_tight", Form("sbs%d_%d",kine,mag) );
  std::string spotanticut_p = jmgr->GetValueFromSubKey_str( "h_dp_spot_anticuts_tight", Form("sbs%d_%d",kine,mag) );

  std::string spotanticutraw = "&&" + spotanticut_n + "&&" + spotanticut_p;

  std::string spotanticut = addbc(spotanticutraw);

  std::cout << std::endl << "Loaded best cluster cuts: " << globalcuts << std::endl;

  std::vector<std::string> cuts = util::parseCuts(globalcuts); //this makes a vector of all individual cuts in the single globalcut string.

  std::vector<std::string> cuts_raw = util::parseCuts(globalcuts_raw);
  
  std::vector<std::string> ccuts;

  std::cout << std::endl << "Parsed cuts: " << std::endl;

  std::string parsedCutString;
  std::string parsedCutString_W2;
  std::string parsedCutString_W2dy;
  std::string parsedCutString_dy;

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

    // Include cut in parsedCutString_dy if it doesn't contain "dy"
    if (trimmedCut.find("dy") == std::string::npos) {
      if (!parsedCutString_dy.empty()) parsedCutString_dy += "&&";
      parsedCutString_dy += trimmedCut;
    }

    ccuts.push_back(cut);
  }

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
  
  // for (size_t i = 0; i < cuts.size(); ++i) {

  //   size_t found_coin = cuts[i].find("coin");
  //   if (found_coin != std::string::npos) {
  //     // Found "coin", replace the first '<' with '>'
  //     size_t found_lt = cuts[i].find('<', found_coin);
  //     if (found_lt != std::string::npos) {
  // 	std::string modified_cut = cuts[i];
  // 	modified_cut[found_lt] = '>';
  // 	coinanticut_noelas = modified_cut; // Overwrite with this modified cut
  // 	// Append modified cut to coinanticut_elas, but exclude "W2" cuts
  // 	coinanticut_elas += (coinanticut_elas.empty() ? "" : "&&") + modified_cut;
  //     }
  //   } else {
  //     // "coin" not found, just append the cut to coinanticut_elas
  //     coinanticut_elas += (coinanticut_elas.empty() ? "" : "&&") + cuts[i];
  //   }
  // }


  //Fast and dirty replacement - will only work if coin time cut isn't improved.
  std::string to_remove = "coin>-10.0";
  std::string to_replace = "coin<10.0";
  std::string replacement = "abs(coin)>10.0";

  for (size_t i = 0; i < cuts_raw.size(); ++i) {
    std::string modified_cut = cuts_raw[i];

    // Remove the substring "coin_bc>-10.0"
    size_t found_remove = modified_cut.find(to_remove);
    if (found_remove != std::string::npos) {
      // Erase the substring and any preceding "&&"
      if (found_remove != 0 && modified_cut[found_remove - 1] == '&') {
	modified_cut.erase(found_remove - 1, to_remove.length() + 2);
      } else {
	modified_cut.erase(found_remove, to_remove.length());
      }
    }

    // Replace the substring "coin_bc<10.0" with "abs(coin_bc)>10.0"
    size_t found_replace = modified_cut.find(to_replace);
    if (found_replace != std::string::npos) {
      modified_cut.replace(found_replace, to_replace.length(), replacement);
    }

    coinanticut_noelas = modified_cut; // Overwrite with this modified cut
    // Append modified cut to coinanticut_elas
    coinanticut_elas += (coinanticut_elas.empty() ? "" : "&&") + modified_cut;
  }

  size_t pos = coinanticut_elas.find("&&&&");
  while (pos != std::string::npos) {
    coinanticut_elas.replace(pos, 4, "&&");
    pos = coinanticut_elas.find("&&&&", pos);
  }
  
  std::cout << "coinanticut_noelas: " << coinanticut_noelas << std::endl;
  std::cout << "coinanticut_elas: " << coinanticut_elas << std::endl;



  //coinanticut_elas = "bb_etot_over_p>0.81&&bb_etot_over_p<1.17&&bb_gem_track_nhits>2&&bb_ps_e>0.234&&abs(bb_tr_vz)<0.073&&abs(coin)>10.0&&dy>-0.36&&dy<0.264&&hcale>0.0224&&hcalon==1&&mag==30&&tar==1&&xexp<0.8&&xexp>-1.3&&yexp<0.45&&yexp>-0.25";



  
  
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

  cout << endl << "    Drawing data histograms from parsed file: " << fin_path << endl << endl << endl;

  std::string skipcorrelationplots_word = "";
  if(skipcorrelationplots)
    skipcorrelationplots_word = "_thin";

  std::string wide_word = "";
  if(wide)
    wide_word = "_widecut";

  std::string fout_path = outdir_path + Form("/gmn_analysis/dx_correlations_bc_sbs%d_mag%d_pass%d%s%s%s.root",kine,mag,pass,skipcorrelationplots_word.c_str(),wide_word.c_str(),effz_word.c_str());

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
  TH1D* hist_base = new TH1D( histBaseName.c_str(), "dx (best cluster);m", hbins, hcalfit_l, hcalfit_h );

  // Draw the plot using the created histogram
  tree->Draw(("dx_bc>>" + histBaseName).c_str(), globalcuts.c_str(), "COLZ");

  // Write the histogram to the output file
  hist_base->Write();

  delete hist_base;

  std::string histBasenodyName = "hdx_allcut_nody";
  TH1D* hist_basenody = new TH1D( histBasenodyName.c_str(), "dx (best cluster);m", hbins, hcalfit_l, hcalfit_h );

  // Draw the plot using the created histogram
  tree->Draw(("dx_bc>>" + histBasenodyName).c_str(), parsedCutString_dy.c_str(), "COLZ");

  // Write the histogram to the output file
  hist_basenody->Write();

  delete hist_basenody;

  std::string histAntidyName = "hdx_dyanti";
  TH1D* hist_antidy = new TH1D( histAntidyName.c_str(), "dx, dy anticut (best cluster);m", hbins, hcalfit_l, hcalfit_h );

  // Draw the plot using the created histogram
  tree->Draw(("dx_bc>>" + histAntidyName).c_str(), dyanticut_elas.c_str(), "COLZ");

  // Write the histogram to the output file
  hist_antidy->Write();

  delete hist_antidy;

  std::string histAnticoinName = "hdx_coinanti";
  TH1D* hist_anticoin = new TH1D( histAnticoinName.c_str(), "dx, coin anticut (best cluster);m", hbins, hcalfit_l, hcalfit_h );

  cout << "FROM" << fin_path << " COINANTICUT APPLIED: " << coinanticut_elas << endl << endl << endl;
  
  // Draw the plot using the created histogram
  //tree->Draw(("dx_bc>>" + histAnticoinName).c_str(), coinanticut_elas.c_str(), "COLZ");

  //tree->Draw("dx>>hdx_coinanti","bb_etot_over_p>0.81&&bb_etot_over_p<1.17&&bb_gem_track_nhits>2&&bb_ps_e>0.234&&abs(bb_tr_vz)<0.073&&abs(coin)>10.0&&dy>-0.36&&dy<0.264&&hcale>0.0224&&hcalon==1&&mag==30&&tar==1&&xexp<0.8&&xexp>-1.3&&yexp<0.45&&yexp>-0.25");

  tree->Draw("dx>>hdx_coinanti",coinanticut_elas.c_str());

    
  // Write the histogram to the output file
  hist_anticoin->Write();

  // TCanvas *c1 = new TCanvas("c1","test",800,600);
  // c1->cd();
  // hdx_coinanti->Draw();
  
  delete hist_anticoin;

  /////////////////////////////////////////
  // Draw the W2 histogram with all cuts
  std::string hW2hist = "hW2";
  std::string hW2_cut = parsedCutString_W2;
  TH1D* hW2 = new TH1D(hW2hist.c_str(), "W2;GeV^{2}", W2bins, W2fit_l, W2fit_h);

  cout << endl << "Plotting best cluster W2 data histogram with cuts: " << hW2_cut << endl << endl;

  // Draw the plot using the created histogram with the W2-excluded cut
  tree->Draw(("W2>>" + hW2hist).c_str(), hW2_cut.c_str(), "COLZ");

  // Write the histogram to the output file
  hW2->Write();

  /////////////////////////////////////////
  // Draw the W2 histogram excluding dy cuts
  std::string hW2dyhist = "hW2_nodycut";
  std::string hW2dy_cut = parsedCutString_W2dy;
  TH1D* hW2dy = new TH1D(hW2dyhist.c_str(), "W2 without dy cuts with spotcut;GeV^{2}", W2bins, W2fit_l, W2fit_h);

  cout << endl << "Plotting best cluster W2 (no dy cut) data histogram with cuts: " << hW2dy_cut << endl << endl;

  // Draw the plot using the created histogram with the W2 and dy-excluded cut
  tree->Draw(("W2>>" + hW2dyhist).c_str(), hW2dy_cut.c_str(), "COLZ");

  // Write the histogram to the output file
  hW2dy->Write();


  /////////////////////////////////////////
  // Draw the W2 histogram excluding dy cuts
  std::string hW2spothist = "hW2_spotcut";
  std::string hW2spot_cut = parsedCutString_W2dy + spotcut;
  TH1D* hW2spot = new TH1D(hW2spothist.c_str(), "W2 with spotcut;GeV^{2}", W2bins, W2fit_l, W2fit_h);

  cout << endl << "Plotting best cluster W2 (spotcut) data histogram with cuts: " << hW2spot_cut << endl << endl;

  // Draw the plot using the created histogram with the W2 and dy-excluded cut
  tree->Draw(("W2>>" + hW2spothist).c_str(), hW2spot_cut.c_str(), "COLZ");

  // Write the histogram to the output file
  hW2spot->Write();


  /////////////////////////////////////////
  // Draw the W2 histogram excluding all cuts
  std::string hW2nocuthist = "hW2_nocutcut";
  std::string hW2nocut_cut = "";
  TH1D* hW2nocut = new TH1D(hW2nocuthist.c_str(), "W2 with no cuts;GeV^{2}", W2bins, W2fit_l, W2fit_h);

  cout << endl << "Plotting best cluster W2 (no cuts) data histogram: " << hW2nocut_cut << endl << endl;

  // Draw the plot using the created histogram with the W2 and dy-excluded cut
  tree->Draw(("W2>>" + hW2nocuthist).c_str(), hW2nocut_cut.c_str(), "COLZ");

  /////////////////////////////////////////
  // Draw the W2 histogram with spot anticuts
  std::string hW2anticuthist = "hW2_anticutcut";
  std::string hW2anticut_cut = parsedCutString_W2dy + spotanticut;
  TH1D* hW2anticut = new TH1D(hW2anticuthist.c_str(), "W2 with elastic and spot anticuts;GeV^{2}", W2bins, W2fit_l, W2fit_h);

  cout << endl << "Plotting best cluster W2 (spot anticut) data histogram: " << hW2anticut_cut << endl << endl;

  // Draw the plot using the created histogram with the W2 and dy-excluded cut
  tree->Draw(("W2>>" + hW2anticuthist).c_str(), hW2anticut_cut.c_str(), "COLZ");

  // Write the histogram to the output file
  hW2anticut->Write();


  ////DO SOME DXDYs for QA
  std::string hdxdyhist = "hdxdy";
  std::string hdxdy_cut = parsedCutString_dy;
  TH2D* hdxdy = new TH2D(hdxdyhist.c_str(), "dx vs dy with elastic cuts;m;m", hbins_dy, hcalfit_l_dy, hcalfit_h_dy, hbins, hcalfit_l, hcalfit_h);

  cout << endl << "Plotting best cluster dx vs dy (elastic cuts) data histogram: " << hdxdy_cut << endl << endl;

  // Draw the plot using the created histogram with the W2 and dy-excluded cut
  tree->Draw(("dx_bc:dy_bc>>" + hdxdyhist).c_str(), hdxdy_cut.c_str(), "COLZ");

  // Write the histogram to the output file
  hdxdy->Write();

  //dxdy spot cuts
  std::string hdxdyspothist = "hdxdy_spotcut";
  std::string hdxdyspot_cut = parsedCutString_dy + spotcut;
  TH2D* hdxdyspot = new TH2D(hdxdyspothist.c_str(), "dx vs dy with elastic and spot cuts;m;m", hbins_dy, hcalfit_l_dy, hcalfit_h_dy, hbins, hcalfit_l, hcalfit_h);

  cout << endl << "Plotting best cluster dx vs dy (elastic + spot cuts) data histogram: " << hdxdyspot_cut << endl << endl;

  // Draw the plot using the created histogram with the W2 and dy-excluded cut
  tree->Draw(("dx_bc:dy_bc>>" + hdxdyspothist).c_str(), hdxdyspot_cut.c_str(), "COLZ");

  // Write the histogram to the output file
  hdxdyspot->Write();

  //dxdy spot anticuts
  std::string hdxdyspotantihist = "hdxdy_spotanticut";
  std::string hdxdyspotanti_cut = parsedCutString_dy + spotanticut;
  TH2D* hdxdyspotanti = new TH2D(hdxdyspotantihist.c_str(), "dx vs dy with elastic and spot anticuts;m;m", hbins_dy, hcalfit_l_dy, hcalfit_h_dy, hbins, hcalfit_l, hcalfit_h);

  cout << endl << "Plotting best cluster dx vs dy (elastic + spot anticut) data histogram: " << hdxdyspotanti_cut << endl << endl;

  // Draw the plot using the created histogram with the W2 and dy-excluded cut
  tree->Draw(("dx_bc:dy_bc>>" + hdxdyspotantihist).c_str(), hdxdyspotanti_cut.c_str(), "COLZ");

  // Write the histogram to the output file
  hdxdyspotanti->Write();

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
    
    //Add best cluster affix for branches that depend on hcal clusters
    std::string bcbranch = branch;
    if(branch.compare("hcale")==0||branch.compare("dy")==0||branch.compare("coin")==0)
      bcbranch += "_bc";

    cout << "Plotting dx versus branch " << bcbranch << " with surviving uncorrelated cuts: " << cutString << endl << endl;

    //Create the histogram
    //std::string histTitle = bcbranch + " vs dx;" + bcbranch + ";dx";
    std::string histTitle = cutString + ";" + bcbranch + ";dx (best cluster)";
    std::string histName = "hist_" + bcbranch;

    TH2D* hist = new TH2D(histName.c_str(), histTitle.c_str(), bins, llim, ulim, hbins, hcalfit_l, hcalfit_h);

    // Draw the plot using the created histogram
    tree->Draw(("dx_bc:" + bcbranch + ">>" + histName).c_str(), cutString.c_str(), "COLZ");

    // Write the histogram to the output file
    hist->Write();

    // Clean up the histogram
    delete hist;

  }
  
  // Create a canvas for displaying cuts
  TCanvas* cutCanvas = new TCanvas("cutCanvas", "Cuts", 800, 600);
  TPaveText* cutsText = new TPaveText(0.1, 0.1, 0.9, 0.9, "NB NDC"); // coordinates are normalized

  cutsText->AddText("Cuts Used:");
  //cutsText->AddLine();
  for (const auto& cut : cuts) {
    cutsText->AddText(cut.c_str());
  }

  cutsText->SetTextAlign(12); // Align text to left
  cutsText->SetFillColor(0);  // Transparent background
  cutsText->Write();

  // Write the canvas to the output file
  cutCanvas->Write();


  if(skipcorrelationplots){
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

  // Close the files
  inputFile->Close();
  delete inputFile;

  outputFile->Close();
  delete outputFile;

  cout << endl << "All plots created. Output file located here: " << fout_path << endl << endl;

  st->Stop();

  // Send time efficiency report to console
  std::cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << std::endl;

}

// trim function for accurate parsing
std::string trim(const std::string& str) {
  size_t first = str.find_first_not_of(" \t\n\r\f\v");
  if (std::string::npos == first) {
    return str;
  }
  size_t last = str.find_last_not_of(" \t\n\r\f\v");
  return str.substr(first, (last - first + 1));
}

// add _bc to get best cluster branches from parse file
std::string addbc(std::string input) {
  const std::string suffix = "_bc";
  const std::string targets[5] = {"coin", "dy", "dx", "hcale", "hcalon"}; //branches that depend on hcal clusters

  for (const auto& target : targets) {
    std::string::size_type pos = 0;
    std::string token = target + suffix;  // Create the new token with the suffix

    // Continue searching the string for the target and replacing it
    while ((pos = input.find(target, pos)) != std::string::npos) {
      // Ensure we match whole words by checking character before and after the match
      bool match = true;
      if (pos != 0 && isalnum(input[pos - 1])) {
	match = false; // Check for character before
      }
      size_t endPos = pos + target.length();
      if (endPos < input.length() && isalnum(input[endPos])) {
	match = false; // Check for character after
      }

      if (match) {
	input.replace(pos, target.length(), token);
	pos += token.length(); // Move past the newly added part to avoid infinite loops
      } else {
	pos += target.length(); // Move past the current word if it's part of a larger word
      }
    }
  }

  return input;
}

// Function to split a string by a delimiter and return a vector of tokens
std::vector<std::string> split(const std::string &s, char delimiter) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter)) {
    tokens.push_back(token);
  }
  return tokens;
}

// Function to get the cell content based on kine, mag, target, and column name
std::map<std::string, std::string> getRowContents(const std::string &filePath, int kine, int mag, const std::string &target, const std::vector<std::string> &excludeKeys) {
  std::ifstream file(filePath);
  std::string line;
  std::map<std::string, std::string> rowContents;
  std::vector<std::string> headers;
  std::set<std::string> excludeSet(excludeKeys.begin(), excludeKeys.end());

  // Read the header line
  if (std::getline(file, line)) {
    headers = split(line, ',');
  }

  // Read the rest of the lines
  while (std::getline(file, line)) {
    std::vector<std::string> values = split(line, ',');
    if (values.size() < 3) {
      continue;
    }

    // Check if the row matches the specified kine, mag, and target
    if (std::stoi(values[0]) == kine && std::stoi(values[1]) == mag && values[2] == target) {
      // Store all column values except the first three in a map
      for (size_t i = 3; i < values.size(); ++i) {
	if (i < headers.size() && excludeSet.find(headers[i]) == excludeSet.end()) {
	  rowContents[headers[i]] = values[i];
	}
      }
      break;
    }
  }

  return rowContents;
}
