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

double binErrLimit = 0.1; //Set limit on HDE error to use HDE prot vs pos exp point for corrections

// Function to concatenate cuts with their names
std::string concatenateCutsWithNames(const std::vector<std::pair<std::string, std::string>>& cutsWithNames) {
  std::string concatenatedCuts;
  for (const auto& cut : cutsWithNames) {
    concatenatedCuts += cut.first + "&&" + cut.second + "&&_____&&";
  }
  return concatenatedCuts;
}

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

Double_t SBpol0rej_b = -0.7; // Central-band fit reject begin, 8.0(-0.9), 8.5(-0.2), 8.7(0.0)
Double_t SBpol0rej_e = -0.4; // Central-band fit reject end, 8.0(-0.5), 8.5(0.1), 8.7(0.5)

Double_t SBpol0rej_b2 = 0.3;
Double_t SBpol0rej_e2 = 0.8;

// Zeroth order poly (constant) with sideband limits
Double_t sbBGfit_pol0(double *x, double *par) {
  Double_t yint = par[0];

  if (x[0] > SBpol0rej_b && x[0] < SBpol0rej_e) {
    TF1::RejectPoint();
    return 0;
  }

  return yint;
}

// Zeroth order poly (constant) with sideband limits
Double_t sbBGfit2_pol0(double *x, double *par) {
  Double_t yint = par[0];

  if (x[0] > SBpol0rej_b && x[0] < SBpol0rej_e) {
    TF1::RejectPoint();
    return 0;
  }

  if (x[0] > SBpol0rej_b2 && x[0] < SBpol0rej_e2) {
    TF1::RejectPoint();
    return 0;
  }

  return yint;
}

//MAIN
void pHDEmap(int kine=11, 
	     int mag=100, 
	     int pass=2, 
	     bool ellipse = true,
	     bool pidcuts = false,
	     bool effz=true,
	     bool mcopt=false) {

  //set draw params
  gStyle->SetPalette(55);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(1111);
  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0110);
  gStyle->SetStatTextColor(kBlack);
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.08);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  
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

  std::string h_spotcut;
  if(ellipse){
    h_spotcut = jmgr->GetValueFromSubKey_str( "h_spot_cuts_tight", Form("sbs%d_%d",kine,mag) );
    cout << "Loaded ellipse proton delta cuts: " << h_spotcut << endl;
  }else{
    h_spotcut = jmgr->GetValueFromSubKey_str( "h_spot_cuts", Form("sbs%d_%d",kine,mag) );
    cout << "Loaded box proton delta cuts: " << h_spotcut << endl;
  }

  std::string h_spotanticut;
  if(ellipse){
    h_spotanticut = jmgr->GetValueFromSubKey_str( "h_spot_anticuts_tight", Form("sbs%d_%d",kine,mag) );
    cout << "Loaded ellipse proton delta anticuts: " << h_spotcut << endl;
  }else{
    h_spotanticut = jmgr->GetValueFromSubKey_str( "h_spot_anticuts", Form("sbs%d_%d",kine,mag) );
    cout << "Loaded box proton delta anticuts: " << h_spotcut << endl;
  }

  //stitch together relevant cuts for later
  std::string e_cuts;
  std::string elas_cuts;
  if(pidcuts){
    e_cuts = magcut + "&&" + trackcut + "&&" + optxcut + "&&" + optycut + "&&" + earmpidcut;
    elas_cuts = earmpidcut;
  }else{
    e_cuts = magcut + "&&" + trackcut + "&&" + optxcut + "&&" + optycut + "&&" + earmcut;
    elas_cuts = earmcut;
  }
  
  std::string tra_cuts = magcut + "&&" + trackcut;
  std::string opt_cuts = magcut + "&&" + optxcut + "&&" + optycut;
  std::string tropt_cuts = magcut + "&&" + trackcut + "&&" + optxcut + "&&" + optycut;
  std::string exp_cuts = tropt_cuts + "&&" + h_spotcut;

  std::string fid_cuts = fidxcut + "&&" + fidycut;

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
  
  std::string mc_norm_cuts = "nucleon==0&&" + optxcut + "&&" + optycut + "&&" + earmcut + "&&" + fidxcut + "&&" + fidycut;
  std::string mc_norm_xcuts = "nucleon==0&&" + optxcut + "&&" + optycut + "&&" + earmcut + "&&" + fidycut;
  std::string mc_norm_ycuts = "nucleon==0&&" + optxcut + "&&" + optycut + "&&" + earmcut + "&&" + fidxcut;

  std::string mc_sig_cuts = "nucleon==0&&" + optxcut + "&&" + optycut + "&&" + earmcut + "&&" + fidxcut + "&&" + fidycut + "&&" + h_spotcut;
  std::string mc_sig_xcuts = "nucleon==0&&" + optxcut + "&&" + optycut + "&&" + earmcut + "&&" + fidycut + "&&" + h_spotcut;
  std::string mc_sig_ycuts = "nucleon==0&&" + optxcut + "&&" + optycut + "&&" + earmcut + "&&" + fidxcut + "&&" + h_spotcut;

  std::string mc_anti_cuts = "nucleon==0&&" + optxcut + "&&" + optycut + "&&" + earmcut + "&&" + fidxcut + "&&" + fidycut + "&&" + h_spotanticut;
  std::string mc_anti_xcuts = "nucleon==0&&" + optxcut + "&&" + optycut + "&&" + earmcut + "&&" + fidycut + "&&" + h_spotanticut;
  std::string mc_anti_ycuts = "nucleon==0&&" + optxcut + "&&" + optycut + "&&" + earmcut + "&&" + fidxcut + "&&" + h_spotanticut;

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
  std::string tar_word = "_lh2";
  std::string mc_word = "";
  if(mcopt)
    mc_word = "_mcana";

  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fin_path;
  if(mcopt)
    fin_path = outdir_path + Form("/parse/parse_mc_sbs%d_%dp_barebones_alt%s.root",kine,mag,effz_word.c_str());
  else
    fin_path = outdir_path + Form("/parse/parse_sbs%d_pass%d_barebones%s%s.root",kine,pass,effz_word.c_str(),tar_word.c_str());

  std::string fout_path = outdir_path + Form("/gmn_analysis/HDE/pHDEmap_sbs%d_mag%d_pass%d%s%s.root",kine,mag,pass,effz_word.c_str(),mc_word.c_str());

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

  // Create a new canvas
  TCanvas *cW2optx = new TCanvas("cW2optx", "W2 vs Optics x", 1000, 600);
  cW2optx->cd();

  // Draw vertical projection to bb midplane with track v cut
  std::string W2optxName = "hW2vOptx";
  std::string W2optxcut = trackcut + "&&" + magcut;
  std::string W2optxTitle = "W^{2} vs BB midplane x (Track Validity Cuts) ";
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

  // Write the histogram to the output file
  cW2optx->Update();
  cW2optx->Write();

  // Create a new canvas
  TCanvas *cW2opty = new TCanvas("cW2opty", "W2 vs Optics y", 1000, 600);
  cW2opty->cd();

  // Draw vertical projection to bb midplane with track v cut vertical projection cut
  std::string W2optyName = "hW2vOpty";
  std::string W2optycut = W2optxcut + "&&" + optxcut;
  std::string W2optyTitle = "W^{2} vs BB midplane y (Track and Optical x Validity Cuts)";
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
  std::string xyexpcuts = tropt_cuts + "&&" + earmcut;
  if(mcopt)
    xyexpcuts += "&&nucleon==1";
  std::string xyexpTitle = "x_{exp} vs y_{exp} (Validity and Elastic Cuts)";
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
  std::string xyexpspotcuts = tropt_cuts + "&&" + earmcut + "&&" + h_spotcut;
  if(mcopt)
    xyexpspotcuts += "&&nucleon==1";
  std::string xyexpspotTitle = "x_{exp} vs y_{exp} (Validity, Elastic, Spot Cuts)";
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
  std::string W2fidxTitle = "W^{2} vs x_{exp} (Validity and Spot Cuts)";
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
  std::string W2fidyTitle = "W^{2} vs y_{exp} (Validity and Spot Cuts)";
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
  if(mcopt)
    dxdycut = mc_norm_cuts;
  std::string dxdyTitle = "#Delta x vs #Delta y (Validity, Fiducial, and Elastic Cuts)";
  TH2D* hdxdy = new TH2D( dxdyName.c_str(), 
			  (dxdyTitle + ";#Delta y (m); #Delta x (m)").c_str(), 
			  hcalybins, 
			  hcalylim.first, 
			  hcalylim.second, 
			  hcalxbins, 
			  hcalxlim.first, 
			  hcalxlim.second);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("dx_bc:dy_bc>>" + dxdyName).c_str(), dxdycut.c_str(), "COLZ");

  // Draw dx vs dy plot for comparison with spot cut region
  std::string dxdyspotName = "hdxdyspot";
  std::string dxdyspotcut = sig_cuts;
  if(mcopt)
    dxdyspotcut = mc_sig_cuts;
  std::string dxdyspotTitle = "#Delta x vs #Delta y (Validity, Fiducial, Elastic, and Spot Cuts)";
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
  TCanvas *cmclus = new TCanvas("cmclus", "multicluster cut comparison", 1200, 600);
  cmclus->Divide(2,1);
  cmclus->cd(1);

  // Draw associated cluster plot (centroids within 60 cm radius) with e-arm and spot cuts
  std::string mclusName = "hmclus";
  std::string mcluscut = sig_cuts;
  if(mcopt)
    mcluscut = mc_sig_cuts;
  std::string mclusTitle = "#Delta x vs #Delta y (Validity, Fiducial, Elastic, and Spot Cuts)";
  TH1D* hmclus = new TH1D( mclusName.c_str(), 
			  (mclusTitle + ";N Associated Clusters").c_str(), 
			  5, 
			  0, 
			  5);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("hcalnclus_mclus>>" + mclusName).c_str(), mcluscut.c_str(), "COLZ");

  // Draw associated cluster plot (centroids within 60 cm radius) with e-arm and antispot cuts
  std::string mclusantiName = "hmclusanti";
  std::string mclusanticut = anti_cuts;
  if(mcopt)
    mclusanticut = mc_anti_cuts;
  std::string mclusantiTitle = "#Delta x vs #Delta y (Validity, Fiducial, Elastic, and Antispot Cuts)";
  TH1D* hmclusanti = new TH1D( mclusantiName.c_str(), 
			  (mclusantiTitle + ";N Associated Clusters").c_str(), 
			  5, 
			  0, 
			  5);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("hcalnclus_mclus>>" + mclusantiName).c_str(), mclusanticut.c_str(), "COLZ");

  // Draw dx vs dy plot for comparison with spot cut region
  std::string dxdymclusName = "hdxdymclus";
  std::string dxdymcluscut = norm_cuts + "&&hcalnclus_mclus>1";
  if(mcopt)
    dxdymcluscut = mc_norm_cuts + "&&hcalnclus_mclus>1";
  std::string dxdymclusTitle = "#Delta x vs #Delta y (Validity, Fiducial, Elastic, and Mclus Cuts)";
  TH2D* hdxdymclus = new TH2D( dxdymclusName.c_str(),
			      (dxdymclusTitle + ";#Delta y (m); #Delta x (m)").c_str(), 
			      hcalybins, 
			      hcalylim.first, 
			      hcalylim.second, 
			      hcalxbins, 
			      hcalxlim.first, 
			      hcalxlim.second);

  // Draw the plot using the created histogram
  tree->Draw(("dx_bc:dy_bc>>" + dxdymclusName).c_str(), dxdymcluscut.c_str(), "COLZ");


  // Draw the mclus plot without spot cut
  hmclus->SetLineColor(kBlack);
  hmclus->GetYaxis()->SetRangeUser(0,1.4*hmclus->GetMaximum());
  hmclus->Draw();

  hmclusanti->SetLineColor(kRed);
  hmclusanti->Draw("same");

  TLegend *leg_mclus = new TLegend(0.11, 0.7, 0.89, 0.89);
  leg_mclus->AddEntry(hmclus, "Validity, Fiducial, and Elastic Cuts", "L");
  leg_mclus->AddEntry(hmclusanti,"Add Antispot Cuts" , "L");
  leg_mclus->Draw();

  // Draw the mclus plot with the spot cut
  cmclus->cd(2);
  hdxdymclus->Draw("colz");

  cmclus->Update();
  cmclus->Write();

  // Create a new canvas
  TCanvas *cW2 = new TCanvas("cW2", "W2", 1200, 600);
  cW2->cd();

  // Draw W2 plot to look at effects of all cuts
  std::string W2Name = "hW2";
  std::string W2cut = tropt_cuts;
  std::string W2Title = "W^{2}";
  TH1D* hW2 = new TH1D( W2Name.c_str(), 
			(W2Title + "; GeV^{2}").c_str(), 
			W2bins, 
			W2lim.first, 
			W2lim.second);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("W2>>" + W2Name).c_str(), W2cut.c_str(), "COLZ");

  // Draw dx vs dy plot for comparison with spot cut region
  std::string W2spotName = "hW2spot";
  std::string W2spotcut = exp_cuts;
  std::string W2spotTitle = "W^{2} (Validity and Spot Cuts)";
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

  TLegend *leg_w2 = new TLegend(0.11, 0.7, 0.49, 0.89);
  leg_w2->AddEntry(hW2, "Validity Cuts", "L");
  leg_w2->AddEntry(hW2spot,"Validity and Spot Cuts" , "L");
  leg_w2->Draw();

  cW2->Update();
  cW2->Write();

  // Draw expected x plot with all cuts, all cuts minus spot, and all cuts antispot (save fidx on all)
  // Create a new canvas
  TCanvas *craw = new TCanvas("craw", "raw", 1200, 600);
  craw->Divide(2,1);
  craw->cd(1);

  std::string xexp_normName = "hxexp_norm";
  std::string xexp_normcut = norm_xcuts;
  if(mcopt)
    xexp_normcut = mc_norm_xcuts;
  std::string xexp_normTitle = "x_{exp}";
  TH1D* hxexp_norm = new TH1D( xexp_normName.c_str(), 
			(xexp_normTitle + "; x_{exp} (m)").c_str(), 
			hcalxbins, 
			hcalxlim.first, 
			hcalxlim.second);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("xexp>>" + xexp_normName).c_str(), xexp_normcut.c_str(), "COLZ");
  hxexp_norm->GetYaxis()->SetRangeUser(0.0,1.5*hxexp_norm->GetMaximum());

  std::string xexp_sigName = "hxexp_sig";
  std::string xexp_sigcut = sig_xcuts;
  if(mcopt)
    xexp_sigcut = mc_sig_xcuts;
  std::string xexp_sigTitle = "";
  TH1D* hxexp_sig = new TH1D( xexp_sigName.c_str(), 
			(xexp_sigTitle + "; x_{exp} (m)").c_str(), 
			hcalxbins, 
			hcalxlim.first, 
			hcalxlim.second);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("xexp>>" + xexp_sigName).c_str(), xexp_sigcut.c_str(), "COLZ");

  std::string xexp_antiName = "hxexp_anti";
  std::string xexp_anticut = anti_xcuts;
  if(mcopt)
    xexp_anticut = mc_anti_xcuts;
  std::string xexp_antiTitle = "";
  TH1D* hxexp_anti = new TH1D( xexp_antiName.c_str(), 
			(xexp_antiTitle + "; x_{exp} (m)").c_str(), 
			hcalxbins, 
			hcalxlim.first, 
			hcalxlim.second);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("xexp>>" + xexp_antiName).c_str(), xexp_anticut.c_str(), "COLZ");


  // Draw expected y plot with all cuts, all cuts minus spot, and all cuts antispot (save fidy on all)
  std::string yexp_normName = "hyexp_norm";
  std::string yexp_normcut = norm_ycuts;
  if(mcopt)
    yexp_normcut = mc_norm_ycuts;
  std::string yexp_normTitle = "y_{exp}";
  TH1D* hyexp_norm = new TH1D( yexp_normName.c_str(), 
			(yexp_normTitle + "; y_{exp} (m)").c_str(), 
			hcalybins, 
			hcalylim.first, 
			hcalylim.second);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("yexp>>" + yexp_normName).c_str(), yexp_normcut.c_str(), "COLZ");
  hyexp_norm->GetYaxis()->SetRangeUser(0.0,1.5*hyexp_norm->GetMaximum());

  std::string yexp_sigName = "hyexp_sig";
  std::string yexp_sigcut = sig_ycuts;
  if(mcopt)
    yexp_sigcut = mc_sig_ycuts;
  std::string yexp_sigTitle = "";
  TH1D* hyexp_sig = new TH1D( yexp_sigName.c_str(), 
			(yexp_sigTitle + "; y_{exp} (m)").c_str(), 
			hcalybins, 
			hcalylim.first, 
			hcalylim.second);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("yexp>>" + yexp_sigName).c_str(), yexp_sigcut.c_str(), "COLZ");

  std::string yexp_antiName = "hyexp_anti";
  std::string yexp_anticut = anti_ycuts;
  if(mcopt)
    yexp_anticut = mc_anti_ycuts;
  std::string yexp_antiTitle = "";
  TH1D* hyexp_anti = new TH1D( yexp_antiName.c_str(), 
			(yexp_antiTitle + "; y_{exp} (m)").c_str(), 
			hcalybins, 
			hcalylim.first, 
			hcalylim.second);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("yexp>>" + yexp_antiName).c_str(), yexp_anticut.c_str(), "COLZ");

  int transparentGreen = TColor::GetColorTransparent(kGreen, 0.3);
  double hxexpmax = hxexp_norm->GetMaximum();
  double hyexpmax = hyexp_norm->GetMaximum();

  TLine *xexpline3 = new TLine(fidxcut_vec[0], 0.0, fidxcut_vec[0], hxexpmax);
  TLine *xexpline4 = new TLine(fidxcut_vec[1], 0.0, fidxcut_vec[1], hxexpmax);
  xexpline3->SetLineColor(kRed);
  xexpline3->SetLineStyle(2);   
  xexpline3->SetLineWidth(2);
  xexpline4->SetLineColor(kRed);
  xexpline4->SetLineStyle(2);   
  xexpline4->SetLineWidth(2);

  TLine *yexpline3 = new TLine(fidycut_vec[0], 0.0, fidycut_vec[0], hyexpmax);
  TLine *yexpline4 = new TLine(fidycut_vec[1], 0.0, fidycut_vec[1], hyexpmax);
  yexpline3->SetLineColor(kRed);
  yexpline3->SetLineStyle(2);   
  yexpline3->SetLineWidth(2);
  yexpline4->SetLineColor(kRed);
  yexpline4->SetLineStyle(2); 
  yexpline4->SetLineWidth(2); 

  // Draw the xexp comparisons
  hxexp_norm->SetLineColor(kBlack);
  hxexp_norm->Draw();
  hxexp_sig->SetFillStyle(1001);
  hxexp_sig->SetFillColor(transparentGreen);
  hxexp_sig->SetLineColor(kGreen);
  hxexp_sig->SetLineWidth(0);
  hxexp_sig->Draw("same");
  hxexp_anti->SetLineColor(kRed);
  hxexp_anti->SetLineWidth(2);
  hxexp_anti->Draw("same");

  //xexpline3->Draw("same");
  //xexpline4->Draw("same");

  TLegend *lxexp = new TLegend(0.11, 0.75, 0.89, 0.89);
  lxexp->AddEntry(hxexp_norm, "Validity, Fiducial y, Elastic Cuts", "l");
  lxexp->AddEntry(hxexp_sig, "Validity, Fiducial y, Elastic, Spot Cuts", "l");
  lxexp->AddEntry(hxexp_anti, "Difference", "l");
  lxexp->Draw();

  craw->Update();

  // Draw the yexp comparisons
  craw->cd(2);
  hyexp_norm->SetLineColor(kBlack);
  hyexp_norm->Draw();
  hyexp_sig->SetFillStyle(1001);
  hyexp_sig->SetFillColor(transparentGreen);
  hyexp_sig->SetLineColor(kGreen);
  hyexp_sig->SetLineWidth(0);
  hyexp_sig->Draw("same");
  hyexp_anti->SetLineColor(kRed);
  hyexp_anti->SetLineWidth(2);
  hyexp_anti->Draw("same");
  
  //yexpline3->Draw("same");
  //yexpline4->Draw("same");

  TLegend *lyexp = new TLegend(0.11, 0.75, 0.89, 0.89);
  lyexp->AddEntry(hyexp_norm, "Validity, Fiducial x, Elastic Cuts", "l");
  lyexp->AddEntry(hyexp_sig, "Validity, Fiducial x, Elastic, Spot Cuts", "l");
  lyexp->AddEntry(hyexp_anti, "Difference", "l");
  lyexp->Draw();

  craw->Update();
  craw->Write();

  // Create canvas and divide it
  TCanvas *cEff = new TCanvas("cEff", "HCal Detection Efficiency", 1200, 600);
  cEff->Divide(2,1);
  cEff->cd(1);

  //Now divide the sig:norm for efficiency per dimension
  std::string yeffName = "hyexp_eff";
  std::string yeffTitle = "HDE proton yexp; y_{exp} (m); Proton Detection Efficiency";
  TH1D* hyexp_eff = new TH1D(yeffName.c_str(), 
			     yeffTitle.c_str(), 
			     hcalybins, 
			     hcalylim.first, 
			     hcalylim.second);

  // Use Divide() to compute the efficiency
  hyexp_eff->Divide(hyexp_sig, hyexp_norm, 1.0, 1.0, "B"); // "B" option to use binomial errors
  hyexp_eff->SetMarkerStyle(20);
  hyexp_eff->SetMarkerColor(kBlack);  

  TH1D* hyexp_corr = (TH1D*)hyexp_eff->Clone("hyexp_corr"); //Clone for later

  std::string xeffName = "hxexp_eff";
  std::string xeffTitle = "HDE proton xexp; x_{exp} (m); Proton Detection Efficiency";
  TH1D* hxexp_eff = new TH1D(xeffName.c_str(), 
			     xeffTitle.c_str(), 
			     hcalxbins, 
			     hcalxlim.first, 
			     hcalxlim.second);

  // Use Divide() to compute the efficiency
  hxexp_eff->Divide(hxexp_sig, hxexp_norm, 1.0, 1.0, "B");
  hxexp_eff->SetMarkerStyle(20);
  hxexp_eff->SetMarkerColor(kBlack);

  TH1D* hxexp_corr = (TH1D*)hxexp_eff->Clone("hxexp_corr"); //Clone for later

  // Fit hyexp_eff within specified bounds
  TF1 *fit_hy = new TF1("fit_hy", "pol0", fidycut_vec[0], fidycut_vec[1]);
  //TF1 *fit_hy = new TF1("fit_hy", "sbBGfit_pol0", fidycut_vec[0], fidycut_vec[1], 1);

  //ensure that no unity no error bins affect the fit
  for (int i = 1; i <= hyexp_eff->GetNbinsX(); ++i) {
    if (hyexp_eff->GetBinError(i) == 0 || hyexp_eff->GetBinContent(i) == 1) {
      hyexp_eff->SetBinContent(i, 0);
      hyexp_eff->SetBinError(i, 0);
    }
  }

  hyexp_eff->Fit(fit_hy, "R");  // "R" option for fit range

  // Fit hxexp_eff from the first bin with data or fidxcut_vec[1] to fidxcut_vec[0]
  int first_nonempty_bin = hxexp_eff->FindFirstBinAbove(0);
  double left_binlimit_x = hxexp_eff->GetBinLowEdge(first_nonempty_bin);
  double fidxcutfirst = (left_binlimit_x > fidxcut_vec[1]) ? left_binlimit_x : fidxcut_vec[1];
  //TF1 *fit_hx = new TF1("fit_hx", "pol0", fidxcutfirst, fidxcut_vec[0]);
  TF1 *fit_hx = new TF1("fit_hx", "sbBGfit_pol0", fidxcutfirst, fidxcut_vec[0], 1);
  //TF1 *fit_hx = new TF1("fit_hx", "pol0", 0.5, fidxcut_vec[0]);

  //ensure that no unity no error bins affect the fit
  for (int i = 1; i <= hxexp_eff->GetNbinsX(); ++i) {
    if (hxexp_eff->GetBinError(i) == 0 || hxexp_eff->GetBinContent(i) == 1) {
      hxexp_eff->SetBinContent(i, 0);
      hxexp_eff->SetBinError(i, 0);
    }
  }

  hxexp_eff->Fit(fit_hx, "R");

  cout << left_binlimit_x << " " << fidxcut_vec[0] << " " << fidxcutfirst << endl;

  // Drawing hyexp_eff with fit and legend
  hxexp_eff->Draw("P");
  fit_hx->SetLineColor(kRed);

  TLegend *leg_x = new TLegend(0.22, 0.12, 0.7, 0.3);
  double xMeanEff = fit_hx->GetParameter(0);
  double xErrEff = fit_hx->GetParError(0);
  leg_x->AddEntry(hxexp_eff, "x_{exp} Detection Efficiency", "P");
  leg_x->AddEntry(fit_hx, Form("Mean Eff = %.3f #pm %.3f", xMeanEff, xErrEff), "L");
  leg_x->Draw();

  // Drawing hxexp_eff with fit and legend
  cEff->cd(2);
  hyexp_eff->Draw("P");
  fit_hy->SetLineColor(kRed);

  TLegend *leg_y = new TLegend(0.22, 0.12, 0.7, 0.3);
  double yMeanEff = fit_hy->GetParameter(0);
  double yErrEff = fit_hy->GetParError(0);
  leg_y->AddEntry(hyexp_eff, "y_{exp} Detection Efficiency", "P");
  leg_y->AddEntry(fit_hy, Form("Mean Eff = %.3f #pm %.3f", yMeanEff, yErrEff), "L");
  leg_y->Draw();

  cEff->Update();
  cEff->Write();

  // Create canvas for map display. Assume mean eff is norm, get correction factor vs x/y exp
  TCanvas *cCorr = new TCanvas("cCorr", "Correction Factors", 1200, 600);
  cCorr->Divide(2,1);
  cCorr->cd(1);

  TLine *xexpline5 = new TLine(fidxcut_vec[0], 0.0, fidxcut_vec[0], 1.1);
  TLine *xexpline6 = new TLine(fidxcut_vec[1], 0.0, fidxcut_vec[1], 1.1);
  xexpline5->SetLineColor(kRed);
  xexpline5->SetLineStyle(2);   
  xexpline5->SetLineWidth(2);
  xexpline6->SetLineColor(kRed);
  xexpline6->SetLineStyle(2);   
  xexpline6->SetLineWidth(2);

  TLine *yexpline5 = new TLine(fidycut_vec[0], 0.0, fidycut_vec[0], 1.1);
  TLine *yexpline6 = new TLine(fidycut_vec[1], 0.0, fidycut_vec[1], 1.1);
  yexpline5->SetLineColor(kRed);
  yexpline5->SetLineStyle(2);   
  yexpline5->SetLineWidth(2);
  yexpline6->SetLineColor(kRed);
  yexpline6->SetLineStyle(2); 
  yexpline6->SetLineWidth(2); 

  //Now divide the eff by mean eff
  std::string xcorrTitle = "Correction factor xexp; x_{exp} (m); Correction factor";

  hxexp_corr->Scale( 1.0 / xMeanEff );
  hxexp_corr->SetTitle(xcorrTitle.c_str());
  hxexp_corr->SetMarkerStyle(20);
  hxexp_corr->SetMarkerColor(kBlue);

  //Fit with pol6 to extract analytic fit to data
  // TF1* xp6map = new TF1("xp6map", "pol6", fidxcutfirst, fidxcut_vec[0]);
  // hxexp_corr->Fit(xp6map, "R");
  // xp6map->SetLineColor(kMagenta);

  // vector<double> xpars;
  // cout << endl << "Correction vs xexp pol6 parameters (p0->p6): " << endl;
  // for( int i; i<7; ++i ){
  //   xpars.push_back(xp6map->GetParameter(i));
  //   cout << xpars[i] << endl;
  // }

  //cout << endl << "Full range of corr map on x: " << fidxcutfirst << " " << fidxcut_vec[0] << endl << endl;

  std::map<double, double> correctionMapx;
  int nBinsx = hxexp_corr->GetNbinsX();
  for (int i = 1; i <= nBinsx; ++i) {
    double binLeftEdge = hxexp_corr->GetXaxis()->GetBinLowEdge(i);
    double binRightEdge = hxexp_corr->GetXaxis()->GetBinUpEdge(i);
    double correctionValue = hxexp_corr->GetBinContent(i);
    // Only add bins within the specified range
    if (binRightEdge >= fidxcutfirst && binLeftEdge <= fidxcut_vec[0]) {
      if(binRightEdge == fidxcutfirst||binLeftEdge == fidxcut_vec[0]) //set to 1 if out of range
	correctionMapx[binRightEdge] = 1.0;
      else if(hxexp_corr->GetBinError(i) > binErrLimit) //set to 1 if the stats low
	correctionMapx[binRightEdge] = 1.0;
      else if(correctionValue>1.05) //set to 1 if HDE artefact
	correctionMapx[binRightEdge] = 1.0;
      else
	correctionMapx[binRightEdge] = correctionValue; //set to the det eff otherwise
      //cout << binRightEdge << " " << correctionValue << endl;
    }
    hxexp_corr->SetBinError(i,0);
  }

  hxexp_corr->Draw("P");

  xexpline5->Draw("same");
  xexpline6->Draw("same");

  cCorr->Update();
  cCorr->cd(2);

  std::string ycorrTitle = "Correction factor yexp; x_{exp} (m); Correction factor";

  hyexp_corr->Scale( 1.0 / yMeanEff );
  hyexp_corr->SetTitle(ycorrTitle.c_str());
  hyexp_corr->SetMarkerStyle(20);
  hyexp_corr->SetMarkerColor(kBlue);

  //Fit with pol6 to extract analytic fit to data
  // TF1* yp6map = new TF1("yp6map", "pol6", fidycut_vec[0], fidycut_vec[1]);
  // hyexp_corr->Fit(yp6map, "R");
  // yp6map->SetLineColor(kMagenta);

  // vector<double> ypars;
  // cout << endl << "Correction vs yexp pol6 parameters (p0->p6): " << endl;
  // for( int i; i<7; ++i ){
  //   ypars.push_back(yp6map->GetParameter(i));
  //   cout << ypars[i] << endl;
  // }

  //cout << endl << "Full range of corr map on y: " << fidycut_vec[0] << " " << fidycut_vec[1] << endl << endl;

  std::map<double, double> correctionMapy;
  int nBinsy = hyexp_corr->GetNbinsX();
  for (int i = 1; i <= nBinsy; ++i) {
    double binLeftEdge = hyexp_corr->GetXaxis()->GetBinLowEdge(i);
    double binRightEdge = hyexp_corr->GetXaxis()->GetBinUpEdge(i);
    double correctionValue = hyexp_corr->GetBinContent(i);
    // Only add bins within the specified range
    if (binRightEdge >= fidycut_vec[1] && binLeftEdge <= fidycut_vec[0]) {
      if(binRightEdge == fidycut_vec[1]||binLeftEdge == fidycut_vec[0]) //set to 1 if out of range
	correctionMapy[binRightEdge] = 1.0;
      else if(hyexp_corr->GetBinError(i) > binErrLimit) //set to 1 if stats low
	correctionMapy[binRightEdge] = 1.0;
      else if(correctionValue>1.05) //set to 1 if HDE artefact
	correctionMapy[binRightEdge] = 1.0;
      else
	correctionMapy[binRightEdge] = correctionValue; //set to the det eff otherwise
      //cout << binRightEdge << " " << correctionValue << endl;
    }
    hyexp_corr->SetBinError(i,0);
  }

  hyexp_corr->Draw("P");

  yexpline5->Draw("same");
  yexpline6->Draw("same");

  //Write correction map out for xexp to file
  std::string xcorrtxt = "effmaps/xcorr_sbs" + std::to_string(kine) + "_mag" + std::to_string(mag) + ".txt";
  std::ofstream xcorrFile(xcorrtxt);
  if (!xcorrFile) {
    std::cerr << "Error opening file for writing: " << xcorrtxt << std::endl;
    return;
  }
  xcorrFile << std::fixed << std::setprecision(6); // Set precision
  xcorrFile << "#HCal Efficiency correction map. List of pairs follow: <x_exp, next upper limit> <correction factor>" << std::endl;
  xcorrFile << "#Created from SBS" << kine << " mag" << mag << std::endl;
  for (const auto& entry : correctionMapx) {
    xcorrFile << entry.first << " " << entry.second << std::endl;
  }
  xcorrFile.close();

  //Write correction map out for yexp to file
  std::string ycorrtxt = "effmaps/ycorr_sbs" + std::to_string(kine) + "_mag" + std::to_string(mag) + ".txt";
  std::ofstream ycorrFile(ycorrtxt);
  if (!ycorrFile) {
    std::cerr << "Error opening file for writing: " << ycorrtxt << std::endl;
    return;
  }
  ycorrFile << std::fixed << std::setprecision(6); // Set precision
  ycorrFile << "#HCal Efficiency correction map. List of pairs follow: <y_exp, next upper limit> <correction factor>" << std::endl;
  ycorrFile << "#Created from SBS" << kine << " mag" << mag << std::endl;
  for (const auto& entry : correctionMapy) {
    ycorrFile << entry.first << " " << entry.second << std::endl;
  }
  ycorrFile.close();

  // //Display all cuts
  // //Track validity only
  // std::string disname = "Track Validity Cuts";
  // TCanvas *cdis1 = new TCanvas("cdis1",disname.c_str(),800,800);
  // util::parseAndDisplayCuts(disname.c_str(), tra_cuts.c_str(), cdis1);

  // cdis1->Write();

  // // //Track + optical x validity
  // // std::string dis2name = "Track and Optical x Validity Cuts";
  // // TCanvas *cdis2 = new TCanvas("cdis2",dis2name.c_str(),800,800);
  // // util::parseAndDisplayCuts(dis2name.c_str(), W2optycut.c_str(), cdis2);

  // // cdis2->Write();

  // //Validity Cuts
  // std::string dis3name = "Optical Validity Cuts";
  // TCanvas *cdis3 = new TCanvas("cdis3",dis3name.c_str(),800,800);
  // util::parseAndDisplayCuts(dis3name.c_str(), opt_cuts.c_str(), cdis3);

  // cdis3->Write();

  // //Elastic Cuts
  // std::string dis4name = "Elastic Cuts";
  // TCanvas *cdis4 = new TCanvas("cdis4",dis4name.c_str(),800,800);
  // util::parseAndDisplayCuts(dis4name.c_str(), elas_cuts.c_str(), cdis4);

  // cdis4->Write();

  // //Proton Spot Cuts
  // std::string dis5name = "Spot Cuts";
  // TCanvas *cdis5 = new TCanvas("cdis5",dis5name.c_str(),800,800);
  // util::parseAndDisplayCuts(dis5name.c_str(), h_spotcut.c_str(), cdis5);

  // cdis5->Write();

  // //Fiducial Cuts
  // std::string dis6name = "Fiducial Cuts";
  // TCanvas *cdis6 = new TCanvas("cdis6",dis6name.c_str(),800,800);
  // util::parseAndDisplayCuts(dis6name.c_str(), fid_cuts.c_str(), cdis6);

  // cdis6->Write();

  // Calculate efficiency over both directions and estimate error
  double teff = (xMeanEff/pow(xErrEff,2)+yMeanEff/pow(yErrEff,2)) / (1/pow(xErrEff,2)+1/pow(yErrEff,2));
  double terr = sqrt(1 / (1/pow(xErrEff,2)+1/pow(yErrEff,2)));


  // Define all cuts and their corresponding names
  std::vector<std::pair<std::string, std::string>> cutsWithNames = {
    {"Track Validity Cuts", tra_cuts},
    // {"Track and Optical x Validity Cuts", W2optycut}, // Uncomment if needed
    {"Optical Validity Cuts", opt_cuts},
    {"Elastic Cuts", elas_cuts},
    {"Spot Cuts", h_spotcut},
    {"Fiducial Cuts", fid_cuts}
  };
  
  // Concatenate all cuts with their names
  std::string allCuts = concatenateCutsWithNames(cutsWithNames);

  std::string efficiencypart = "&&_____&&Efficiency&&xexp: " + 
    std::to_string(xMeanEff) + " +/- " + 
    std::to_string(xErrEff) + "&&yexp: " + 
    std::to_string(yMeanEff) + " +/- " + 
    std::to_string(yErrEff) + "&&combined: " +
    std::to_string(teff) + " +/- " + 
    std::to_string(terr);

  allCuts += efficiencypart;

  // Display all concatenated cuts on a single canvas
  std::string disnameall = "All Cuts";
  TCanvas *cdisall = new TCanvas("cdis", disnameall.c_str(), 1000, 800);
  util::parseAndDisplayCuts(disnameall.c_str(), allCuts.c_str(), cdisall);

  cout << "All plots created. Output file located here: " << fout_path << endl;

}
