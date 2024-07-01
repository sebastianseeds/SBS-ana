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
#include <vector>
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

//gets values for optics v x and y cuts for plotting
std::vector<double> extractValues(const std::string& expression) {
  std::vector<double> values;
  std::regex pattern("[<>]=?\\s*(-?\\d*\\.?\\d+)");
  std::smatch matches;

  std::string::const_iterator searchStart(expression.cbegin());
  while (std::regex_search(searchStart, expression.cend(), matches, pattern)) {
    values.push_back(std::stod(matches[1].str()));
    searchStart = matches.suffix().first;
  }

  return values;
}

//MAIN. Spots are 5x dy sig, 3x dx sig
void pnrHDEmap(int kine=4, 
	     int mag=30, 
	     int pass=2, 
	     bool ellipse = true,
	     bool pidcuts = false,
	     bool effz=true) {

  //set draw params
  gStyle->SetPalette(55);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  
  cout << "Processing analysis of HCal detection efficiency for kinematic " << kine << " at field " << mag << "." << endl << endl;

  // reading json configuration file
  JSONManager *jmgr = new JSONManager("../../config/hdemap.json");

  //Get cuts
  std::string magcut = jmgr->GetValueFromSubKey_str( "mag_cuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded magnetic field cuts: " << magcut << endl;

  std::string trackcut = jmgr->GetValueFromSubKey_str( "track_v_cuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded track validity cuts: " << trackcut << endl;

  std::string optxcut = jmgr->GetValueFromSubKey_str( "optx_v_cuts", Form("sbs%d_%d",kine,mag) );
  vector<double> optxcut_vec = extractValues(optxcut);
  cout << "Loaded optics x validity cuts: " << optxcut << " (values: " << optxcut_vec[0] << ", " << optxcut_vec[1] << ")" << endl;

  std::string optycut = jmgr->GetValueFromSubKey_str( "opty_v_cuts", Form("sbs%d_%d",kine,mag) );
  vector<double> optycut_vec = extractValues(optycut);
  cout << "Loaded optics y validity cuts: " << optycut << " (values: " << optycut_vec[0] << ", " << optycut_vec[1] << ")" << endl;

  std::string fidxcut = jmgr->GetValueFromSubKey_str( "fidx_cut", Form("sbs%d_%d",kine,mag) );
  vector<double> fidxcut_vec = extractValues(fidxcut);
  cout << "Loaded fiducial x cut: " << fidxcut << " (values: " << fidxcut_vec[0] << ", " << fidxcut_vec[1] << ")" << endl;

  std::string fidycut = jmgr->GetValueFromSubKey_str( "fidy_cut", Form("sbs%d_%d",kine,mag) );
  vector<double> fidycut_vec = extractValues(fidycut);
  cout << "Loaded fiducial y cut: " << fidycut << " (values: " << fidycut_vec[0] << ", " << fidycut_vec[1] << ")" << endl;

  std::string earmcut = jmgr->GetValueFromSubKey_str( "earm_cuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded electron arm cuts: " << earmcut << endl;

  std::string earmpidcut = jmgr->GetValueFromSubKey_str( "earmpid_cuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded electron arm with PID cuts: " << earmcut << endl;

  //Spots are for LD2 now
  std::string h_prot_spotcut;
  if(ellipse){
    h_prot_spotcut = jmgr->GetValueFromSubKey_str( "h_dp_spot_cuts_tight", Form("sbs%d_%d",kine,mag) );
    cout << "Loaded ellipse proton delta cuts: " << h_prot_spotcut << endl;
  }else{
    h_prot_spotcut = jmgr->GetValueFromSubKey_str( "h_dp_spot_cuts", Form("sbs%d_%d",kine,mag) );
    cout << "Loaded box proton delta cuts: " << h_prot_spotcut << endl;
  }

  std::string h_prot_spotanticut;
  if(ellipse){
    h_prot_spotanticut = jmgr->GetValueFromSubKey_str( "h_dp_spot_anticuts_tight", Form("sbs%d_%d",kine,mag) );
    cout << "Loaded ellipse proton delta anticuts: " << h_prot_spotcut << endl;
  }else{
    h_prot_spotanticut = jmgr->GetValueFromSubKey_str( "h_dp_spot_anticuts", Form("sbs%d_%d",kine,mag) );
    cout << "Loaded box proton delta anticuts: " << h_prot_spotcut << endl;
  }

  std::string h_neut_spotcut;
  if(ellipse){
    h_neut_spotcut = jmgr->GetValueFromSubKey_str( "h_dn_spot_cuts_tight", Form("sbs%d_%d",kine,mag) );
    cout << "Loaded ellipse neutron delta cuts: " << h_neut_spotcut << endl;
  }else{
    h_neut_spotcut = jmgr->GetValueFromSubKey_str( "h_dn_spot_cuts", Form("sbs%d_%d",kine,mag) );
    cout << "Loaded box neutron delta cuts: " << h_neut_spotcut << endl;
  }

  std::string h_neut_spotanticut;
  if(ellipse){
    h_neut_spotanticut = jmgr->GetValueFromSubKey_str( "h_dn_spot_anticuts_tight", Form("sbs%d_%d",kine,mag) );
    cout << "Loaded ellipse neutron delta anticuts: " << h_neut_spotcut << endl;
  }else{
    h_neut_spotanticut = jmgr->GetValueFromSubKey_str( "h_dn_spot_anticuts", Form("sbs%d_%d",kine,mag) );
    cout << "Loaded box neutron delta anticuts: " << h_neut_spotcut << endl;
  }

  //stitch together relevant cuts for later
  std::string h_spotcut = "((" + h_prot_spotcut + ")||(" + h_neut_spotcut + "))";
  std::string h_spotanticut = "(" + h_prot_spotanticut + ")&&(" + h_neut_spotanticut + ")";

  std::string tropt_cuts = magcut + "&&" + trackcut + "&&" + optxcut + "&&" + optycut;
  std::string e_cuts;
  if(pidcuts)
    e_cuts = tropt_cuts + "&&" + earmpidcut;
  else
    e_cuts = tropt_cuts + "&&" + earmcut;

  std::string exp_cuts = magcut + "&&" + trackcut + "&&" + optxcut + "&&" + optycut + "&&" + h_spotcut;
  
  //cuts for exracting efficiency
  std::string sig_cuts = e_cuts + "&&" + fidxcut + "&&" + fidycut + "&&" + h_spotcut;
  std::string norm_cuts = e_cuts + "&&" + fidxcut + "&&" + fidycut;
  std::string anti_cuts = e_cuts + "&&" + fidxcut + "&&" + fidycut + "&&" + h_spotanticut;

  std::string sig_xcuts = e_cuts + "&&" + fidycut + "&&" + h_spotcut;
  std::string norm_xcuts = e_cuts + "&&" + fidycut;
  std::string anti_xcuts = e_cuts + "&&" + fidycut + "&&" + h_spotanticut;

  std::string sig_ycuts = e_cuts + "&&" + fidxcut + "&&" + h_spotcut;
  std::string norm_ycuts = e_cuts + "&&" + fidxcut;
  std::string anti_ycuts = e_cuts + "&&" + fidxcut + "&&" + h_spotanticut;

  //select eff cuts by nucleon
  std::string sig_xpcuts = e_cuts + "&&" + fidycut + "&&" + h_prot_spotcut;
  std::string norm_xpcuts = e_cuts + "&&" + fidycut;
  std::string anti_xpcuts = e_cuts + "&&" + fidycut + "&&" + h_prot_spotanticut;

  std::string sig_ypcuts = e_cuts + "&&" + fidxcut + "&&" + h_prot_spotcut;
  std::string norm_ypcuts = e_cuts + "&&" + fidxcut;
  std::string anti_ypcuts = e_cuts + "&&" + fidxcut + "&&" + h_prot_spotanticut;

  std::string sig_xncuts = e_cuts + "&&" + fidycut + "&&" + h_neut_spotcut;
  std::string norm_xncuts = e_cuts + "&&" + fidycut;
  std::string anti_xncuts = e_cuts + "&&" + fidycut + "&&" + h_neut_spotanticut;

  std::string sig_yncuts = e_cuts + "&&" + fidxcut + "&&" + h_neut_spotcut;
  std::string norm_yncuts = e_cuts + "&&" + fidxcut;
  std::string anti_yncuts = e_cuts + "&&" + fidxcut + "&&" + h_neut_spotanticut;  

  //Get plot limits/bins
  int bbmidxbins = jmgr->GetValueFromSubKey<int>( "bbbinsx", Form("sbs%d_%d",kine,mag) );
  int bbmidybins = jmgr->GetValueFromSubKey<int>( "bbbinsy", Form("sbs%d_%d",kine,mag) );
  int hcalxbins = jmgr->GetValueFromSubKey<int>( "hbinsx", Form("sbs%d_%d",kine,mag) );
  int hcalybins = jmgr->GetValueFromSubKey<int>( "hbinsy", Form("sbs%d_%d",kine,mag) );
  int W2bins = jmgr->GetValueFromSubKey<int>( "W2bins", Form("sbs%d_%d",kine,mag) );

  std::pair<double,double> bbmidxlim = {
    jmgr->GetValueFromSubKey<double>( "bbmidx_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "bbmidx_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded bbmidx plot limits min: " << bbmidxlim.first << ", max: " << bbmidxlim.second << ", with bins: " << bbmidxbins << endl;

  std::pair<double,double> bbmidylim = {
    jmgr->GetValueFromSubKey<double>( "bbmidy_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "bbmidy_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded bbmidy plot limits min: " << bbmidylim.first << ", max: " << bbmidylim.second << ", with bins: " << bbmidybins <<  endl;

  std::pair<double,double> hcalxlim = {
    jmgr->GetValueFromSubKey<double>( "hcalx_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "hcalx_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded hcalx plot limits min: " << hcalxlim.first << ", max: " << hcalxlim.second << ", with bins: " << hcalxbins <<  endl;

  std::pair<double,double> hcalylim = {
    jmgr->GetValueFromSubKey<double>( "hcaly_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "hcaly_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded hcaly plot limits min: " << hcalylim.first << ", max: " << hcalylim.second << ", with bins: " << hcalybins <<  endl;

  std::pair<double,double> W2lim = {
    jmgr->GetValueFromSubKey<double>( "W2_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "W2_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded W2 plot limits min: " << W2lim.first << ", max: " << W2lim.second << ", with bins: " << W2bins <<  endl;
 
  //set up files and paths
  std::string effz_word = "";
  if(effz)
    effz_word = "_effz";
  std::string tar_word = "_ld2";

  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fin_path = outdir_path + Form("/parse/parse_sbs%d_pass%d_barebones%s%s.root",kine,pass,effz_word.c_str(),tar_word.c_str());
  std::string fout_path = outdir_path + Form("/gmn_analysis/pnrHDEmap_sbs%d_mag%d_pass%d%s.root",kine,mag,pass,effz_word.c_str());

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

  ///////////////////////
  ///////////////////////
  ///////////////////////
  ///////////////////////
  ///////////////////////
  ///////////////////////
  //Begin Plots
  ///////////////////////
  ///////////////////////
  ///////////////////////

  // Draw vertical projection to bb midplane with track v cut
  std::string W2optxName = "hW2vOptx";
  std::string W2optxcut = trackcut + "&&" + magcut;
  std::string W2optxTitle = "(hW2vOptx) " + W2optxcut;
  TH2D* W2optx = new TH2D( W2optxName.c_str(), 
			   (W2optxTitle + ";bb_tr_r_x-0.9*bb_tr_r_th; GeV^{2}").c_str(), 
			   bbmidxbins, 
			   bbmidxlim.first, 
			   bbmidxlim.second, 
			   W2bins, 
			   W2lim.first, 
			   W2lim.second);

  // Draw the plot using the created histogram
  tree->Draw(("W2:(bb_tr_r_x-0.9*bb_tr_r_th)>>" + W2optxName).c_str(), W2optxcut.c_str(), "COLZ");

  // Create a new canvas
  TCanvas *cW2optx = new TCanvas("cW2optx", "W2 vs Optics x", 1000, 600);
  cW2optx->cd();

  // Draw the histogram on the canvas
  W2optx->Draw("COLZ");

  // Get cuts from json and make TLines
  TLine *xline1 = new TLine(optxcut_vec[0], W2lim.first, optxcut_vec[0], W2lim.second);
  TLine *xline2 = new TLine(optxcut_vec[1], W2lim.first, optxcut_vec[1], W2lim.second);

  xline1->SetLineColor(kRed);
  xline1->SetLineStyle(2);   
  xline1->SetLineWidth(2);
  xline2->SetLineColor(kRed);
  xline2->SetLineStyle(2); 
  xline2->SetLineWidth(2); 

  xline1->Draw("same");
  xline2->Draw("same");

  cW2optx->Update();
  cW2optx->Write();

  // Create a new canvas
  TCanvas *cW2opty = new TCanvas("cW2opty", "W2 vs Optics y", 1000, 600);
  cW2opty->cd();

  // Draw vertical projection to bb midplane with track v cut vertical projection cut
  std::string W2optyName = "hW2vOpty";
  std::string W2optycut = W2optxcut + "&&" + optxcut;
  std::string W2optyTitle = "(hW2vOpty) " + W2optycut;
  TH2D* W2opty = new TH2D( W2optyName.c_str(), 
			   (W2optyTitle + ";bb_tr_r_y-0.9*bb_tr_r_ph; GeV^{2}").c_str(), 
			   bbmidybins, 
			   bbmidylim.first, 
			   bbmidylim.second, 
			   W2bins, 
			   W2lim.first, 
			   W2lim.second);

  // Draw the plot using the created histogram
  tree->Draw(("W2:(bb_tr_r_y-0.9*bb_tr_r_ph)>>" + W2optyName).c_str(), W2optycut.c_str(), "COLZ");

  // Draw the histogram on the canvas
  W2opty->Draw("COLZ");

  // Get cuts from json and make TLines
  TLine *yline1 = new TLine(optycut_vec[0], W2lim.first, optycut_vec[0], W2lim.second);
  TLine *yline2 = new TLine(optycut_vec[1], W2lim.first, optycut_vec[1], W2lim.second);

  yline1->SetLineColor(kRed);
  yline1->SetLineStyle(2);   
  yline1->SetLineWidth(2);
  yline2->SetLineColor(kRed);
  yline2->SetLineStyle(2); 
  yline2->SetLineWidth(2); 

  yline1->Draw("same");
  yline2->Draw("same");

  // Write the histogram to the output file
  cW2opty->Update();
  cW2opty->Write();

  // Create a new canvas
  TCanvas *cxyexp = new TCanvas("cxyexp", "xexp vs yexp", 1200, 600);
  cxyexp->Divide(2,1);
  cxyexp->cd(1);

  // Draw the xexp heat maps
  std::string xyexpName = "hxyexp";
  std::string xyexpcuts;
  if(pidcuts)
    xyexpcuts = tropt_cuts + "&&" + earmpidcut;
  else
    xyexpcuts = tropt_cuts + "&&" + earmcut;
  std::string xyexpTitle = "(hxyexp) " + xyexpcuts;
  TH2D* hxyexp = new TH2D( xyexpName.c_str(), 
			   (xyexpTitle + ";y_{exp} (m); x_{exp} (m)").c_str(), 
			   hcalybins, 
			   hcalylim.first, 
			   hcalylim.second, 
			   hcalxbins, 
			   hcalxlim.first, 
			   hcalxlim.second);

  // Draw the plot using the created histogram
  tree->Draw(("xexp:yexp>>" + xyexpName).c_str(), xyexpcuts.c_str(), "COLZ");

  std::string xyexpspotName = "hxyexpspot";
  std::string xyexpspotcuts;
  if(pidcuts)
    xyexpspotcuts = tropt_cuts + "&&" + earmpidcut + "&&" + h_spotcut;
  else
    xyexpspotcuts = tropt_cuts + "&&" + earmcut + "&&" + h_spotcut;

  std::string xyexpspotTitle = "(hxyexpspot) " + xyexpspotcuts;
  TH2D* hxyexpspot = new TH2D( xyexpspotName.c_str(), 
			       (xyexpspotTitle + ";y_{exp} (m); x_{exp} (m)").c_str(), 
			       hcalybins, 
			       hcalylim.first, 
			       hcalylim.second, 
			       hcalxbins, 
			       hcalxlim.first, 
			       hcalxlim.second);

  // Draw the plot using the created histogram
  tree->Draw(("xexp:yexp>>" + xyexpspotName).c_str(), xyexpspotcuts.c_str(), "COLZ");

  // Bottom line (horizontal) at y = fidxcut_vec[0]
  TLine *xexpline1 = new TLine(fidycut_vec[0], fidxcut_vec[0], fidycut_vec[1], fidxcut_vec[0]);

  // Bottom line (horizontal) at y = fidxcut_vec[0]
  TLine *xexpline2 = new TLine(fidycut_vec[0], fidxcut_vec[1], fidycut_vec[1], fidxcut_vec[1]);

  // Left side line (vertical) at x = fidycut_vec[0]
  TLine *yexpline1 = new TLine(fidycut_vec[0], fidxcut_vec[0], fidycut_vec[0], fidxcut_vec[1]);

  // Right side line (vertical) at x = fidycut_vec[1]
  TLine *yexpline2 = new TLine(fidycut_vec[1], fidxcut_vec[0], fidycut_vec[1], fidxcut_vec[1]);

  // Draw the histogram on the canvas
  int xyexp_nev = hxyexp->GetEntries();
  hxyexp->Draw("COLZ");

  xexpline1->SetLineColor(kRed);
  xexpline1->SetLineStyle(2);   
  xexpline1->SetLineWidth(2);
  xexpline2->SetLineColor(kRed);
  xexpline2->SetLineStyle(2);   
  xexpline2->SetLineWidth(2);

  xexpline1->Draw("same");
  xexpline2->Draw("same");

  yexpline1->SetLineColor(kRed);
  yexpline1->SetLineStyle(2);   
  yexpline1->SetLineWidth(2);
  yexpline2->SetLineColor(kRed);
  yexpline2->SetLineStyle(2); 
  yexpline2->SetLineWidth(2); 

  yexpline1->Draw("same");
  yexpline2->Draw("same");

  TLegend *leg_exp = new TLegend(0.11, 0.7, 0.5, 0.89);
  leg_exp->AddEntry((TObject*)0, Form("N Ev: %d", xyexp_nev), "");
  leg_exp->AddEntry(xexpline1, "Fiducial cuts", "l");
  leg_exp->Draw();

  cxyexp->cd(2);

  // Draw the histogram on the canvas
  int xyexpspot_nev = hxyexpspot->GetEntries();
  hxyexpspot->Draw("COLZ");

  xexpline1->Draw("same");
  xexpline2->Draw("same");

  yexpline1->Draw("same");
  yexpline2->Draw("same");

  TLegend *leg_expspot = new TLegend(0.11, 0.7, 0.5, 0.89);
  leg_expspot->AddEntry((TObject*)0, Form("N Ev: %d", xyexpspot_nev), "");
  leg_expspot->AddEntry(xexpline1, "Fiducial cuts", "l");
  leg_expspot->Draw(); 

  // Write the canvas to the output file
  cxyexp->Update();
  cxyexp->Write();

  // Create a new canvas
  TCanvas *cW2fidx = new TCanvas("cW2fidx", "W2 vs xexp", 1000, 600);
  cW2fidx->cd();

  // Draw the W2 fiducial x cut histogram
  std::string W2fidxName = "hW2vFidx";
  std::string W2fidxTitle = "(hW2vFidx) " + exp_cuts;
  TH2D* W2fidx = new TH2D( W2fidxName.c_str(), 
			   (W2fidxTitle + ";x_{exp} (m); GeV^{2}").c_str(), 
			   hcalxbins, 
			   hcalxlim.first, 
			   hcalxlim.second, 
			   W2bins, 
			   W2lim.first, 
			   W2lim.second);

  // Draw the plot using the created histogram
  tree->Draw(("W2:xexp>>" + W2fidxName).c_str(), exp_cuts.c_str(), "COLZ");

  // Draw the histogram on the canvas
  W2fidx->Draw("COLZ");

  // Get cuts from json and make TLines
  TLine *xfidline1 = new TLine(fidxcut_vec[0], W2lim.first, fidxcut_vec[0], W2lim.second);
  TLine *xfidline2 = new TLine(fidxcut_vec[1], W2lim.first, fidxcut_vec[1], W2lim.second);

  xfidline1->SetLineColor(kRed);
  xfidline1->SetLineStyle(2);   
  xfidline1->SetLineWidth(2);
  xfidline2->SetLineColor(kRed);
  xfidline2->SetLineStyle(2);   
  xfidline2->SetLineWidth(2);

  xfidline1->Draw("same");
  xfidline2->Draw("same");

  // Write the histogram to the output file
  cW2fidx->Update();
  cW2fidx->Write();

  // Create a new canvas
  TCanvas *cW2fidy = new TCanvas("cW2fidy", "W2 vs yexp", 1000, 600);
  cW2fidy->cd();

  // Draw the fiducial y cut histogram
  std::string W2fidyName = "hW2vFidy";
  std::string W2fidyTitle = "(hW2vFidy) " + exp_cuts;
  TH2D* W2fidy = new TH2D( W2fidyName.c_str(), 
			   (W2fidyTitle + ";y_{exp} (m); GeV^{2}").c_str(), 
			   hcalybins, 
			   hcalylim.first, 
			   hcalylim.second, 
			   W2bins, 
			   W2lim.first, 
			   W2lim.second);

  // Draw the plot using the created histogram
  tree->Draw(("W2:yexp>>" + W2fidyName).c_str(), exp_cuts.c_str(), "COLZ");

  // Draw the histogram on the canvas
  W2fidy->Draw("COLZ");

  // Get cuts from json and make TLines
  TLine *yfidline1 = new TLine(fidycut_vec[0], W2lim.first, fidycut_vec[0], W2lim.second);
  TLine *yfidline2 = new TLine(fidycut_vec[1], W2lim.first, fidycut_vec[1], W2lim.second);

  yfidline1->SetLineColor(kRed);
  yfidline1->SetLineStyle(2);   
  yfidline1->SetLineWidth(2);
  yfidline2->SetLineColor(kRed);
  yfidline2->SetLineStyle(2); 
  yfidline2->SetLineWidth(2); 

  yfidline1->Draw("same");
  yfidline2->Draw("same");

  // Write the histogram to the output file
  cW2fidy->Update();
  cW2fidy->Write();

  // Create a new canvas
  TCanvas *cdxdy = new TCanvas("cdxdy", "dx vs dy", 1200, 600);
  cdxdy->Divide(2,1);
  cdxdy->cd(1);

  // Draw dx vs dy plot for comparison with spot cut region
  std::string dxdyName = "hdxdy";
  std::string dxdycut = norm_cuts;
  std::string dxdyTitle = "(hdxdy) " + dxdycut;
  TH2D* hdxdy = new TH2D( dxdyName.c_str(), 
			  (dxdyTitle + ";#Delta y (m); #Delta x (m)").c_str(), 
			  hcalybins, 
			  hcalylim.first, 
			  hcalylim.second, 
			  hcalxbins, 
			  hcalxlim.first, 
			  hcalxlim.second);

  // Draw the histogram
  tree->Draw(("dx_bc:dy_bc>>" + dxdyName).c_str(), dxdycut.c_str(), "COLZ");

  // Draw dx vs dy plot for comparison with spot cut region
  std::string dxdyspotName = "hdxdyspot";
  std::string dxdyspotcut = sig_cuts;
  std::string dxdyspotTitle = "(hdxdyspot) " + dxdyspotcut;
  TH2D* hdxdyspot = new TH2D( dxdyspotName.c_str(),
			      (dxdyspotTitle + ";#Delta y (m); #Delta x (m)").c_str(), 
			      hcalybins, 
			      hcalylim.first, 
			      hcalylim.second, 
			      hcalxbins, 
			      hcalxlim.first, 
			      hcalxlim.second);

  // Draw the plot using the created histogram
  tree->Draw(("dx_bc:dy_bc>>" + dxdyspotName).c_str(), dxdyspotcut.c_str(), "COLZ");

  // Draw the dxdy plot without spot cut
  hdxdy->Draw("COLZ");

  // Draw the dxdy plot with the spot cut
  cdxdy->cd(2);
  hdxdyspot->Draw("colz");

  cdxdy->Update();
  cdxdy->Write();

  // Create a new canvas
  TCanvas *cW2 = new TCanvas("cW2", "W2", 1200, 600);
  cW2->cd();

  // Draw W2 plot to look at effects of all cuts
  std::string W2Name = "hW2";
  std::string W2cut = tropt_cuts;
  std::string W2Title = "(hW2) " + W2cut;
  TH1D* hW2 = new TH1D( W2Name.c_str(), 
			(W2Title + "; GeV^{2}").c_str(), 
			W2bins, 
			W2lim.first, 
			W2lim.second);

  // Draw the histogram
  tree->Draw(("W2>>" + W2Name).c_str(), W2cut.c_str(), "COLZ");

  // Draw dx vs dy plot for comparison with spot cut region
  std::string W2spotName = "hW2spot";
  std::string W2spotcut = tropt_cuts + "&&" + h_spotcut;
  std::string W2spotTitle = "(hW2spot) " + W2spotcut;
  TH1D* hW2spot = new TH1D( W2spotName.c_str(), 
			    (W2spotTitle + "; GeV^{2}").c_str(), 
			    W2bins, 
			    W2lim.first, 
			    W2lim.second);

  // Draw the plot using the created histogram
  tree->Draw(("W2>>" + W2spotName).c_str(), W2spotcut.c_str(), "COLZ");

  // Draw the W2 plot without spot cut
  hW2->SetLineColor(kBlack);
  hW2->Draw();

  // Draw the W2 plot with the spot cut
  int transparentBlue = TColor::GetColorTransparent(kBlue, 0.3);
  hW2spot->SetLineColor(kBlue);
  hW2spot->SetLineWidth(0);
  hW2spot->SetFillColor(transparentBlue);
  hW2spot->SetFillStyle(1001);
  hW2spot->Draw("SAME");

  TLegend *leg_w2 = new TLegend(0.5, 0.7, 0.89, 0.89);
  leg_w2->AddEntry(hW2, "Electron Arm Cuts", "L");
  leg_w2->AddEntry(hW2spot,"Electron Arm and Nucleon Spot Cuts" , "L");
  leg_w2->Draw();

  cW2->Update();
  cW2->Write();


  // Create a new canvas
  TCanvas *craw = new TCanvas("craw", "raw", 1200, 600);
  craw->Divide(2,1);
  craw->cd(1);

  // Draw expected x plot with all cuts, all cuts minus spot, and all cuts antispot (save fidx on all)
  std::string xexp_normName = "hxexp_norm";
  std::string xexp_normcut = norm_xcuts;
  std::string xexp_normTitle = "(hxexp_norm) " + xexp_normcut;
  TH1D* hxexp_norm = new TH1D( xexp_normName.c_str(), 
			(xexp_normTitle + "; x_{exp} (m)").c_str(), 
			hcalxbins, 
			hcalxlim.first, 
			hcalxlim.second);

  // Draw the histogram
  tree->Draw(("xexp>>" + xexp_normName).c_str(), xexp_normcut.c_str(), "COLZ");

  std::string xexp_sigName = "hxexp_sig";
  std::string xexp_sigcut = sig_xcuts;
  std::string xexp_sigTitle = "(hxexp_sig) " + xexp_sigcut;
  TH1D* hxexp_sig = new TH1D( xexp_sigName.c_str(), 
			(xexp_sigTitle + "; x_{exp} (m)").c_str(), 
			hcalxbins, 
			hcalxlim.first, 
			hcalxlim.second);

  // Draw the histogram
  tree->Draw(("xexp>>" + xexp_sigName).c_str(), xexp_sigcut.c_str(), "COLZ");

  std::string xexp_antiName = "hxexp_anti";
  std::string xexp_anticut = anti_xcuts;
  std::string xexp_antiTitle = "(hxexp_anti) " + xexp_anticut;
  TH1D* hxexp_anti = new TH1D( xexp_antiName.c_str(), 
			(xexp_antiTitle + "; x_{exp} (m)").c_str(), 
			hcalxbins, 
			hcalxlim.first, 
			hcalxlim.second);

  // Draw the histogram
  tree->Draw(("xexp>>" + xexp_antiName).c_str(), xexp_anticut.c_str(), "COLZ");


  // Draw expected y plot with all cuts, all cuts minus spot, and all cuts antispot (save fidy on all)
  std::string yexp_normName = "hyexp_norm";
  std::string yexp_normcut = norm_ycuts;
  std::string yexp_normTitle = "(hyexp_norm) " + yexp_normcut;
  TH1D* hyexp_norm = new TH1D( yexp_normName.c_str(), 
			(yexp_normTitle + "; y_{exp} (m)").c_str(), 
			hcalybins, 
			hcalylim.first, 
			hcalylim.second);

  // Draw the histogram
  tree->Draw(("yexp>>" + yexp_normName).c_str(), yexp_normcut.c_str(), "COLZ");

  std::string yexp_sigName = "hyexp_sig";
  std::string yexp_sigcut = sig_ycuts;
  std::string yexp_sigTitle = "(hyexp_sig) " + yexp_sigcut;
  TH1D* hyexp_sig = new TH1D( yexp_sigName.c_str(), 
			(yexp_sigTitle + "; y_{exp} (m)").c_str(), 
			hcalybins, 
			hcalylim.first, 
			hcalylim.second);

  // Draw the histogram
  tree->Draw(("yexp>>" + yexp_sigName).c_str(), yexp_sigcut.c_str(), "COLZ");

  std::string yexp_antiName = "hyexp_anti";
  std::string yexp_anticut = anti_ycuts;
  std::string yexp_antiTitle = "(hyexp_anti) " + yexp_anticut;
  TH1D* hyexp_anti = new TH1D( yexp_antiName.c_str(), 
			(yexp_antiTitle + "; y_{exp} (m)").c_str(), 
			hcalybins, 
			hcalylim.first, 
			hcalylim.second);

  // Draw the histogram
  tree->Draw(("yexp>>" + yexp_antiName).c_str(), yexp_anticut.c_str(), "COLZ");

  // Draw the xexp comparisons
  hxexp_norm->Draw();
  hxexp_sig->SetFillStyle(3002);
  hxexp_sig->SetFillColor(kGreen+2);
  hxexp_sig->SetLineColor(kGreen+2);
  hxexp_sig->Draw("same");
  hxexp_anti->SetLineColor(kRed);
  hxexp_anti->SetLineWidth(2);
  hxexp_anti->Draw("same");

  // Draw the yexp comparisons
  craw->cd(2);
  hyexp_norm->Draw();
  hyexp_sig->SetFillStyle(3002);
  hyexp_sig->SetFillColor(kGreen+2);
  hyexp_sig->SetLineColor(kGreen+2);
  hyexp_sig->Draw("same");
  hyexp_anti->SetLineColor(kRed);
  hyexp_anti->SetLineWidth(2);
  hyexp_anti->Draw("same");
  
  craw->Update();
  craw->Write();

  // Create canvas and divide it
  TCanvas *cEff = new TCanvas("cEff", "HCal Detection Efficiency", 1200, 600);
  cEff->Divide(2,1);
  cEff->cd(1);

  //Now divide the sig:norm for efficiency per dimension
  std::string yeffName = "hyexp_eff";
  std::string yeffTitle = "HDE nucleon yexp; y_{exp} (m); Nucleon Detection Efficiency";
  TH1D* hyexp_eff = new TH1D(yeffName.c_str(), 
			     yeffTitle.c_str(), 
			     hcalybins, 
			     hcalylim.first, 
			     hcalylim.second);

  // Use Divide() to compute the efficiency
  hyexp_eff->Divide(hyexp_sig, hyexp_norm, 1.0, 1.0, "B"); // "B" option to use binomial errors
  hyexp_eff->SetMarkerStyle(20);
  hyexp_eff->SetMarkerColor(kBlack);  

  std::string xeffName = "hxexp_eff";
  std::string xeffTitle = "HDE nucleon xexp; x_{exp} (m); Nucleon Detection Efficiency";
  TH1D* hxexp_eff = new TH1D(xeffName.c_str(), 
			     xeffTitle.c_str(), 
			     hcalxbins, 
			     hcalxlim.first, 
			     hcalxlim.second);

  // Use Divide() to compute the efficiency
  hxexp_eff->Divide(hxexp_sig, hxexp_norm, 1.0, 1.0, "B");
  hxexp_eff->SetMarkerStyle(20);
  hxexp_eff->SetMarkerColor(kBlack);

  // Fit hyexp_eff within bounds
  TF1 *fit_hy = new TF1("fit_hy", "pol0", fidycut_vec[0], fidycut_vec[1]);
  hyexp_eff->Fit(fit_hy, "R");

  // Fit hxexp_eff from the first bin with data to fidxcut_vec[0]
  int first_nonempty_bin = hxexp_eff->FindFirstBinAbove(0);
  double left_binlimit_x = hxexp_eff->GetBinLowEdge(first_nonempty_bin);
  double fidxcutfirst = (left_binlimit_x > fidxcut_vec[1]) ? left_binlimit_x : fidxcut_vec[1];
  TF1 *fit_hx = new TF1("fit_hx", "pol0", fidxcutfirst, fidxcut_vec[0]);
  hxexp_eff->Fit(fit_hx, "R");

  //cout << left_binlimit_x << " " << fidxcut_vec[0] << " " << fidxcutfirst << endl;

  // Drawing hyexp_eff with fit and legend
  hxexp_eff->Draw("P");
  fit_hx->SetLineColor(kRed);

  TLegend *leg_x = new TLegend(0.3, 0.3, 0.7, 0.6);
  leg_x->AddEntry(hxexp_eff, "x_{exp} Detection Efficiency", "P");
  leg_x->AddEntry(fit_hx, Form("Avg Eff = %.3f #pm %.3f", fit_hx->GetParameter(0), fit_hx->GetParError(0)), "L");
  leg_x->Draw();

  // Drawing hxexp_eff with fit and legend
  cEff->cd(2);
  hyexp_eff->Draw("P");
  fit_hy->SetLineColor(kRed);

  TLegend *leg_y = new TLegend(0.3, 0.3, 0.7, 0.6);
  leg_y->AddEntry(hyexp_eff, "y_{exp} Detection Efficiency", "P");
  leg_y->AddEntry(fit_hy, Form("Avg Eff = %.3f #pm %.3f", fit_hy->GetParameter(0), fit_hy->GetParError(0)), "L");
  leg_y->Draw();

  cEff->Update();
  cEff->Write();

  ////////////////////////
  ////////////////////////
  // Now look at n_det/N_exp / p_exp/N_exp and compare to MC sig_rn / sig_rp
  ////////////////////////
  ////////////////////////
  ////////////////////////

  //PROTON
   // Draw expected x plot with all cuts, all cuts minus spot, and all cuts antispot (save fidx on all)
  std::string xpexp_normName = "hxpexp_norm";
  std::string xpexp_normcut = norm_xpcuts;
  std::string xpexp_normTitle = "(hxpexp_norm) " + xpexp_normcut;
  TH1D* hxpexp_norm = new TH1D( xpexp_normName.c_str(), 
			(xpexp_normTitle + "; x_{exp} (m)").c_str(), 
			hcalxbins, 
			hcalxlim.first, 
			hcalxlim.second);

  // Draw the histogram
  tree->Draw(("xexp>>" + xpexp_normName).c_str(), xpexp_normcut.c_str(), "COLZ");

  std::string xpexp_sigName = "hxpexp_sig";
  std::string xpexp_sigcut = sig_xpcuts;
  std::string xpexp_sigTitle = "(hxpexp_sig) " + xpexp_sigcut;
  TH1D* hxpexp_sig = new TH1D( xpexp_sigName.c_str(), 
			(xpexp_sigTitle + "; x_{exp} (m)").c_str(), 
			hcalxbins, 
			hcalxlim.first, 
			hcalxlim.second);

  // Draw the histogram
  tree->Draw(("xexp>>" + xpexp_sigName).c_str(), xpexp_sigcut.c_str(), "COLZ");

  std::string xpexp_antiName = "hxpexp_anti";
  std::string xpexp_anticut = anti_xpcuts;
  std::string xpexp_antiTitle = "(hxpexp_anti) " + xpexp_anticut;
  TH1D* hxpexp_anti = new TH1D( xpexp_antiName.c_str(), 
			(xpexp_antiTitle + "; x_{exp} (m)").c_str(), 
			hcalxbins, 
			hcalxlim.first, 
			hcalxlim.second);

  // Draw the histogram
  tree->Draw(("xexp>>" + xpexp_antiName).c_str(), xpexp_anticut.c_str(), "COLZ");


  // Draw expected y plot with all cuts, all cuts minus spot, and all cuts antispot (save fidy on all)
  std::string ypexp_normName = "hypexp_norm";
  std::string ypexp_normcut = norm_ypcuts;
  std::string ypexp_normTitle = "(hypexp_norm) " + ypexp_normcut;
  TH1D* hypexp_norm = new TH1D( ypexp_normName.c_str(), 
			(ypexp_normTitle + "; y_{exp} (m)").c_str(), 
			hcalybins, 
			hcalylim.first, 
			hcalylim.second);

  // Draw the histogram
  tree->Draw(("yexp>>" + ypexp_normName).c_str(), ypexp_normcut.c_str(), "COLZ");

  std::string ypexp_sigName = "hypexp_sig";
  std::string ypexp_sigcut = sig_ypcuts;
  std::string ypexp_sigTitle = "(hypexp_sig) " + ypexp_sigcut;
  TH1D* hypexp_sig = new TH1D( ypexp_sigName.c_str(), 
			(ypexp_sigTitle + "; y_{exp} (m)").c_str(), 
			hcalybins, 
			hcalylim.first, 
			hcalylim.second);

  // Draw the histogram
  tree->Draw(("yexp>>" + ypexp_sigName).c_str(), ypexp_sigcut.c_str(), "COLZ");

  std::string ypexp_antiName = "hypexp_anti";
  std::string ypexp_anticut = anti_ypcuts;
  std::string ypexp_antiTitle = "(hypexp_anti) " + ypexp_anticut;
  TH1D* hypexp_anti = new TH1D( ypexp_antiName.c_str(), 
			(ypexp_antiTitle + "; y_{exp} (m)").c_str(), 
			hcalybins, 
			hcalylim.first, 
			hcalylim.second);

  // Draw the histogram
  tree->Draw(("yexp>>" + ypexp_antiName).c_str(), ypexp_anticut.c_str(), "COLZ");

  //NEUTRON
  // Draw expected x plot with all cuts, all cuts minus spot, and all cuts antispot (save fidx on all)
  std::string xnexp_normName = "hnpexp_norm";
  std::string xnexp_normcut = norm_xncuts;
  std::string xnexp_normTitle = "(hxnexp_norm) " + xnexp_normcut;
  TH1D* hxnexp_norm = new TH1D( xnexp_normName.c_str(), 
			(xnexp_normTitle + "; x_{exp} (m)").c_str(), 
			hcalxbins, 
			hcalxlim.first, 
			hcalxlim.second);

  // Draw the histogram
  tree->Draw(("xexp>>" + xnexp_normName).c_str(), xnexp_normcut.c_str(), "COLZ");

  std::string xnexp_sigName = "hxnexp_sig";
  std::string xnexp_sigcut = sig_xncuts;
  std::string xnexp_sigTitle = "(hxnexp_sig) " + xnexp_sigcut;
  TH1D* hxnexp_sig = new TH1D( xnexp_sigName.c_str(), 
			(xnexp_sigTitle + "; x_{exp} (m)").c_str(), 
			hcalxbins, 
			hcalxlim.first, 
			hcalxlim.second);

  // Draw the histogram
  tree->Draw(("xexp>>" + xnexp_sigName).c_str(), xnexp_sigcut.c_str(), "COLZ");

  std::string xnexp_antiName = "hxnexp_anti";
  std::string xnexp_anticut = anti_xncuts;
  std::string xnexp_antiTitle = "(hxnexp_anti) " + xnexp_anticut;
  TH1D* hxnexp_anti = new TH1D( xnexp_antiName.c_str(), 
			(xnexp_antiTitle + "; x_{exp} (m)").c_str(), 
			hcalxbins, 
			hcalxlim.first, 
			hcalxlim.second);

  // Draw the histogram
  tree->Draw(("xexp>>" + xnexp_antiName).c_str(), xnexp_anticut.c_str(), "COLZ");


  // Draw expected y plot with all cuts, all cuts minus spot, and all cuts antispot (save fidy on all)
  std::string ynexp_normName = "hynexp_norm";
  std::string ynexp_normcut = norm_yncuts;
  std::string ynexp_normTitle = "(hynexp_norm) " + ynexp_normcut;
  TH1D* hynexp_norm = new TH1D( ynexp_normName.c_str(), 
			(ynexp_normTitle + "; y_{exp} (m)").c_str(), 
			hcalybins, 
			hcalylim.first, 
			hcalylim.second);

  // Draw the histogram
  tree->Draw(("yexp>>" + ynexp_normName).c_str(), ynexp_normcut.c_str(), "COLZ");

  std::string ynexp_sigName = "hynexp_sig";
  std::string ynexp_sigcut = sig_yncuts;
  std::string ynexp_sigTitle = "(hynexp_sig) " + ynexp_sigcut;
  TH1D* hynexp_sig = new TH1D( ynexp_sigName.c_str(), 
			(ynexp_sigTitle + "; y_{exp} (m)").c_str(), 
			hcalybins, 
			hcalylim.first, 
			hcalylim.second);

  // Draw the histogram
  tree->Draw(("yexp>>" + ynexp_sigName).c_str(), ynexp_sigcut.c_str(), "COLZ");

  std::string ynexp_antiName = "hynexp_anti";
  std::string ynexp_anticut = anti_yncuts;
  std::string ynexp_antiTitle = "(hynexp_anti) " + ynexp_anticut;
  TH1D* hynexp_anti = new TH1D( ynexp_antiName.c_str(), 
			(ynexp_antiTitle + "; y_{exp} (m)").c_str(), 
			hcalybins, 
			hcalylim.first, 
			hcalylim.second);

  // Draw the yexp plot
  tree->Draw(("yexp>>" + ynexp_antiName).c_str(), ynexp_anticut.c_str(), "COLZ");

  // Create a new canvas
  TCanvas *crawall = new TCanvas("crawall", "raw", 1800, 1200);
  crawall->Divide(2,2);

  // Draw the proton xexp comparisons
  int transparentGreen = TColor::GetColorTransparent(kGreen, 0.3);
  crawall->cd(1);
  hxpexp_norm->Draw();
  hxpexp_sig->SetFillStyle(1001);
  hxpexp_sig->SetFillColor(transparentGreen);
  hxpexp_sig->SetLineColor(kGreen);
  hxpexp_sig->SetLineWidth(0);
  hxpexp_sig->Draw("same");
  hxpexp_anti->SetLineColor(kRed);
  hxpexp_anti->SetLineWidth(2);
  hxpexp_anti->Draw("same");

  // Draw the proton yexp comparisons
  crawall->cd(2);
  hypexp_norm->Draw();
  hypexp_sig->SetFillStyle(1001);
  hypexp_sig->SetFillColor(transparentGreen);
  hypexp_sig->SetLineColor(kGreen);
  hypexp_sig->SetLineWidth(0);
  hypexp_sig->Draw("same");
  hypexp_anti->SetLineColor(kRed);
  hypexp_anti->SetLineWidth(2);
  hypexp_anti->Draw("same");

  // Draw the neutron xexp comparisons
  crawall->cd(3);
  hxnexp_norm->Draw();
  hxnexp_sig->SetFillStyle(1001);
  hxnexp_sig->SetFillColor(transparentGreen);
  hxnexp_sig->SetLineColor(kGreen);
  hxnexp_sig->SetLineWidth(0);
  hxnexp_sig->Draw("same");
  hxnexp_anti->SetLineColor(kRed);
  hxnexp_anti->SetLineWidth(2);
  hxnexp_anti->Draw("same");

  // Draw the neutron yexp comparisons
  crawall->cd(4);
  hynexp_norm->Draw();
  hynexp_sig->SetFillStyle(1001);
  hynexp_sig->SetFillColor(transparentGreen);
  hynexp_sig->SetLineColor(kGreen);
  hynexp_sig->SetLineWidth(0);
  hynexp_sig->Draw("same");
  hynexp_anti->SetLineColor(kRed);
  hynexp_anti->SetLineWidth(2);
  hynexp_anti->Draw("same");
  
  crawall->Update();
  crawall->Write();

  //Get the ratios
  // Create canvas and divide it
  TCanvas *cnpr = new TCanvas("cnpr", "n:p detected ratios", 1200, 600);
  cnpr->Divide(2,1);
  cnpr->cd(1);

  //n:p detected, X
  std::string nprxName = "hnprx";
  std::string nprxTitle = "n_{det}:p_{det}, x_{exp}; x_{exp} (m); n:p ratio";
  TH1D* hnprx = new TH1D(nprxName.c_str(), 
			 nprxTitle.c_str(), 
			 hcalxbins, 
			 hcalxlim.first, 
			 hcalxlim.second);

  hnprx->Divide(hxnexp_sig, hxpexp_sig, 1.0, 1.0, "B");
  hnprx->GetYaxis()->SetRangeUser(0.0,1.5*hnprx->GetMaximum());
  hnprx->SetMarkerStyle(20);
  hnprx->SetMarkerColor(kBlack);

  //n:p detected, Y
  std::string npryName = "hnpry";
  std::string npryTitle = "n_{det}:p_{det}, y_{exp}; y_{exp} (m); n:p ratio";
  TH1D* hnpry = new TH1D(npryName.c_str(), 
			 npryTitle.c_str(), 
			 hcalybins, 
			 hcalylim.first, 
			 hcalylim.second);

  hnpry->Divide(hynexp_sig, hypexp_sig, 1.0, 1.0, "B");
  hnprx->GetYaxis()->SetRangeUser(0.0,1.5*hnprx->GetMaximum());
  hnpry->SetMarkerStyle(20);
  hnpry->SetMarkerColor(kBlack);

  // pol0 fit, x
  TF1 *fit_nprx = new TF1("fit_nprx", "pol0", fidxcutfirst, fidxcut_vec[0]);
  hnprx->Fit(fit_nprx, "R");

  // pol0 fit, y
  TF1 *fit_npry = new TF1("fit_npry", "pol0", fidycut_vec[0], fidycut_vec[1]);
  hnpry->Fit(fit_npry, "R");

  // Drawing hnprx with fit and legend
  hnprx->Draw("P");
  fit_nprx->SetLineColor(kRed);

  TLegend *leg_nprx = new TLegend(0.11, 0.75, 0.89, 0.89);
  leg_nprx->AddEntry(hnprx, "x_{exp} n:p detected ratio", "P");
  leg_nprx->AddEntry(fit_nprx, Form("n:p detected ratio mean = %.3f #pm %.3f", fit_nprx->GetParameter(0), fit_nprx->GetParError(0)), "L");
  leg_nprx->Draw();

  cnpr->Update();

  // Drawing hnprx with fit and legend
  cnpr->cd(2);
  hnpry->Draw("P");
  fit_npry->SetLineColor(kRed);

  TLegend *leg_npry = new TLegend(0.11, 0.75, 0.89, 0.89);
  leg_npry->AddEntry(hnpry, "y_{exp} n:p detected ratio", "P");
  leg_npry->AddEntry(fit_npry, Form("n:p detected ratio mean = %.3f #pm %.3f", fit_npry->GetParameter(0), fit_npry->GetParError(0)), "L");
  leg_npry->Draw();

  cnpr->Update();
  cnpr->Write();

  cout << "All plots created. Output file located here: " << fout_path << endl;

}
