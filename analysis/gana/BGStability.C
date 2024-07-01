//sseeds 5.22.23: Script to vary backgrounds and extract GMn with systematic error estimation from the spread of these extractions.

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TKey.h>
#include <TClass.h>
#include <iostream>
#include <vector>
#include <string>
#include <TLatex.h>
#include <algorithm>
#include <TRandom.h>
#include "../../include/gmn.h"
#include "../../src/jsonmgr.C"

//Fit range
double hcalfit_l; //default defined by json
double hcalfit_h; //default defined by json

double wideband_offset = 0.3; //adds and subtracts this value from the fit limits for pol4 and sideband fits

double HDE_syst_N = 94.316; //percent error on nucleon yield (from proton HDE analysis, assuming similar)
double HDE_systerr_N = 0.583; //percent error on neutron yield (from proton HDE analysis, assuming similar)

//Fit options
std::string fitopt = "RMQ0";

//get random seed
string passkey = "hytp";

//Total fits using Interpolate with elastic signal histo and 4th order poly fit to bg
TH1D *hdx_p;
TH1D *hdx_n;
TH1D *hdx_p_wide;
TH1D *hdx_n_wide;
TH1D *hdx_inel;
TH1D *hdx_antidy;
TH1D *hdx_anticoin;

Double_t fit_gausBG(double *x, double *par) {
  // MC float params for scaling and shifting
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p_wide->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n_wide->Interpolate(x[0] - dx_shift_n);
  
  // Sum of components
  return proton + neutron + fits::g_gfit(x, &par[4]);
}

Double_t fit_pol4BG(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p_wide->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n_wide->Interpolate(x[0] - dx_shift_n);
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron + fits::g_p4fit(x, &par[4]);
}

// Double_t fit_pol4BG(double *x, double *par){
//   // MC float params
//   double proton_scale = par[0];
//   double neutron_scale = par[1];
//   double dx_shift_p = par[2]; // Shift for proton histogram
//   double dx_shift_n = par[3]; // Shift for neutron histogram  

//   // Apply shifts before interpolation
//   double proton = proton_scale * hdx_p_wide->Interpolate(x[0] - dx_shift_p);
//   double neutron = neutron_scale * hdx_n_wide->Interpolate(x[0] - dx_shift_n);
  
//   // Use the remaining parameters for fits::g_p4fit, starting from par[4]
//   return proton + neutron + 
//     abs(par[4]) + abs(par[5]) * x[0] + abs(par[6]) * x[0] * x[0] + abs(par[7]) * x[0] * x[0] * x[0] + abs(par[8]) * x[0] * x[0] * x[0] * x[0];
// }

double fit_xmin;
double fit_xmax;

// Double_t fit_pol4BG(double *x, double *par){
//   // MC float params
//   double proton_scale = par[0];
//   double neutron_scale = par[1];
//   double dx_shift_p = par[2]; // Shift for proton histogram
//   double dx_shift_n = par[3]; // Shift for neutron histogram  

//   // Apply shifts before interpolation
//   double proton = proton_scale * hdx_p_wide->Interpolate(x[0] - dx_shift_p);
//   double neutron = neutron_scale * hdx_n_wide->Interpolate(x[0] - dx_shift_n);
  
//   // Calculate the polynomial background
//   double bg = fits::g_p4fit(x, &par[4]);

//   // Compute the baseline (straight line connecting the endpoints)
//   double x_start = fit_xmin;  // Define the start of the fit range
//   double x_end = fit_xmax;    // Define the end of the fit range
//   double y_start = fits::g_p4fit(&x_start, &par[4]); // Value at start
//   double y_end = fits::g_p4fit(&x_end, &par[4]);     // Value at end

//   double slope = (y_end - y_start) / (x_end - x_start);
//   double intercept = y_start - slope * x_start;

//   double baseline = slope * x[0] + intercept;

//   // Ensure the fit result is never below the baseline
//   double fit_result = proton + neutron + bg;
//   if (fit_result < baseline) {
//     fit_result = baseline;
//   }

//   return fit_result;
// }

Double_t fit_pol2BG(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron + fits::g_p2fit_cd(x, &par[4]);
}

Double_t fit_pol1BG(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron + fits::g_p1fit(x, &par[4]);
}

Double_t fit_pol3BG(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron + fits::g_p3fit(x, &par[4]);
}

Double_t fit_inelBG(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  
  double bg_scale = par[4];

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);
  double bg = bg_scale * hdx_inel->Interpolate(x[0] - dx_shift_n);
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron + bg;
}

Double_t fit_dyantiBG(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  
  double bg_scale = par[4];

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);
  double bg = bg_scale * hdx_antidy->Interpolate(x[0]); //no shift for data BG
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron + bg;
}

Double_t fit_coinantiBG(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  
  double bg_scale = par[4];

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);
  double bg = bg_scale * hdx_anticoin->Interpolate(x[0]); //no shift for data BG
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron + bg;
}

Double_t fit_noBG(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);
  
  return proton + neutron;
}

Double_t SBpol4rej_b; //Central-band fit begin
Double_t SBpol4rej_e; //Central-band fit end

//Fourth order poly with sideband limits
Double_t sbBGfit(double *x, double *par){

  Double_t yint = par[0];
  Double_t p1 = par[1];
  Double_t p2 = par[2];
  Double_t p3 = par[3];
  Double_t p4 = par[4];

  if(x[0]>SBpol4rej_b && x[0]<SBpol4rej_e) { 
    TF1::RejectPoint();
    return 0;
  }

  return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3)+p4*pow(x[0],4);
}

// Forward declarations
void handleError(TFile *file1, TFile *file2, std::string marker);
void handleError(TFile *file1, std::string marker);
void createAndSaveCanvasWithFitsAndResiduals(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, TF1* bg, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, int pol, double &elastics, std::pair<double,double> &rsf, int kine, int mag);
void createAndSaveCanvasWithFitsAndResiduals_nobg(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, double &elastics, std::pair<double,double> &rsf, int kine, int mag);
void createAndSaveCanvasWithFitsAndResiduals_altbg(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, TH1D* hdxinel, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, std::string bgtype, double &elastics, std::pair<double,double> &rsf, int kine, int mag);
std::vector<std::pair<double, double>> fitAndFineFit(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, const std::string& fitOptions = "RBMQ0");
std::vector<std::pair<double, double>> fitAndFineFit_fixshift(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, double pshift, double nshift, const std::string& fitOptions = "RBMQ0");
void subtractBackground(TH1D* hist, TF1* bgFit);
int FindCutIndex(const std::string& branches, const std::string& searchStr);
std::string RemovePrefix(const std::string& originalString, const std::string& prefix);
TH1D* calculateResiduals(TH1D* hData, TF1* fit, const char* residualName);
void chi2_constrained(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
std::pair<double, double> calculateChiSquareWithZeroMCError(TH1D* histogram, TF1* fitFunc);
void plotBackgrounds(const std::vector<TF1*>& bgFunctions);

/////////////////
/////////////////
/////////////////
/////////////////
/////////////////
/////////////////
/////////////////

bool s4f50 = false;

//MAIN. kine=kinematic, mag=fieldsetting, pass=pass#, *_override=p/n shift par override, blind=add blinding, bestclus=use best cluster plots, thin=use plots without correlations, widecut=use plots with wide cuts, effz=use plots with effective z implemented, alt=use plots from alternate MC files
void BGStability(int kine=8, 
		 int mag=100, 
		 int pass=2, 
		 double nshift_override=0.0, 
		 double pshift_override=0.0, 
		 bool blind=false, 
		 bool bestclus=true, 
		 bool thin=true,
		 bool widecut=false,
		 bool effz=true,
		 bool alt = true) {
  //set draw params
  gStyle->SetPalette(55);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetLabelSize(0.03, "XYZ");
  gStyle->SetTitleSize(0.04, "XYZ");
  gStyle->SetEndErrorSize(0);

  // Set up configuration and tune objects to load analysis parameters
  SBSconfig config(kine,0);

  // Obtain configuration pars from config class
  double E_beam = config.GetEbeam();
  double BB_angle = config.GetBBtheta_rad(); // In radians

  //load json file
  JSONManager *jmgr = new JSONManager("../../config/syst.json");

  //Get dx plot details
  hcalfit_l = jmgr->GetValueFromSubKey<double>( "hcalfit_l", Form("sbs%d_%d",kine,mag) ); //gets to 7cm HCal pos res
  hcalfit_h = jmgr->GetValueFromSubKey<double>( "hcalfit_h", Form("sbs%d_%d",kine,mag) ); //gets to 7cm HCal pos res

  std::cout << "Loaded hcal fit limits lower(upper): " << hcalfit_l << "(" << hcalfit_h << ")" << std::endl;

  fit_xmin = hcalfit_l;
  fit_xmax = hcalfit_h;
  
  //set sideband
  SBpol4rej_b = -1.8; 
  SBpol4rej_e = 0.7;

  //Get tight elastic cuts
  std::string globalcuts;
  if(!widecut){
    globalcuts = jmgr->GetValueFromSubKey_str( Form("post_tcuts_p%d",pass), Form("sbs%d_%d",kine,mag) );
    cout << "Loaded wide cuts: " << globalcuts << endl;

  }else{
    globalcuts = jmgr->GetValueFromSubKey_str( Form("post_cuts_p%d",pass), Form("sbs%d_%d",kine,mag) );
    cout << "Loaded tight cuts: " << globalcuts << endl;

  }

  std::vector<std::string> cuts = util::parseCuts(globalcuts); //this makes a vector of all individual cuts in the single globalcut string.

  std::cout << "Parsed cuts: " << std::endl;
  for ( size_t i=0; i<cuts.size(); ++i )
    cout << cuts[i] << endl;

  //get paths and set random seed
  Double_t blind_factor = 1.;
  if(blind)
    blind_factor = util::generateRandomNumber(passkey);

  std::string bestclus_word = "";
  if(bestclus)
    bestclus_word = "_bc";

  std::string thin_word = "";
  if(thin)
    thin_word = "_thin";

  std::string alt_word = "";
  if(alt)
    alt_word = "_alt";

  std::string wide_word = "";
  if(widecut)
    wide_word = "_widecut";

  std::string effz_word = "";
  if(effz)
    effz_word = "_effz";

  std::string basePath = "/lustre19/expphy/volatile/halla/sbs/seeds";

  std::string finPath = Form("%s/gmn_analysis/dx_correlations%s_sbs%d_mag%d_pass%d%s%s%s.root", basePath.c_str(), bestclus_word.c_str(), kine, mag, pass, thin_word.c_str(),wide_word.c_str(), effz_word.c_str());
  std::string fmcinPath = Form("%s/gmn_analysis/dx_mc_sbs%d_mag%d_pass%d%s%s%s%s.root", basePath.c_str(), kine, mag, pass, thin_word.c_str(), alt_word.c_str(),wide_word.c_str(), effz_word.c_str());
  std::string foutPath = Form("%s/gmn_analysis/bgstab_sbs%d_mag%d_pass%d%s%s.root", basePath.c_str(), kine, mag, pass, bestclus_word.c_str(), effz_word.c_str());

  std::string fmcinelinPath = Form("%s/gmn_analysis/dx_mc_sbs%d_mag%d_pass%d_rd2.root", basePath.c_str(), kine, mag, pass);

  if(kine>4)
    fmcinelinPath = fmcinPath;
  
  TFile* inputFile = new TFile(finPath.c_str());
  if (!inputFile || inputFile->IsZombie()) 
    handleError(inputFile,"inputFile");

  //Get histograms from data plot file
  TH1D *hdx_data_raw = dynamic_cast<TH1D*>(inputFile->Get("hdx_allcut"));
  TH1D *hdx_data = (TH1D*)(hdx_data_raw->Clone("hdx_data"));
  TH1D *hdx_data_pol3test = (TH1D*)(hdx_data_raw->Clone("hdx_data_pol3test"));
  //TH1D *hdx_data = (TH1D*)(hdx_projection->Clone("hdx_data"));
  hdx_data->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  TH1D *hdx_data_wide = (TH1D*)(hdx_data_raw->Clone("hdx_data_wide"));
  hdx_data_wide->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  TH1D *hdx_data_set = (TH1D*)(hdx_data_raw->Clone("hdx_data_set"));
  hdx_data_set->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  TH1D *hdx_data_sb = (TH1D*)(hdx_data_raw->Clone("hdx_data_sb"));
  hdx_data_sb->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  TH1D *hdx_sb_nobg = (TH1D*)(hdx_data_raw->Clone("hdx_sb_nobg"));
  hdx_sb_nobg->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  //Get anticut histograms
  TH1D *hdx_antidy_raw = dynamic_cast<TH1D*>(inputFile->Get("hdx_dyanti"));
  hdx_antidy = (TH1D*)(hdx_antidy_raw->Clone("hdx_antidy"));
  hdx_antidy->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  TH1D *hdx_antidy_clone = (TH1D*)(hdx_antidy_raw->Clone("hdx_antidy_clone"));
  hdx_antidy_clone->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  TH1D *hdx_anticoin_raw = dynamic_cast<TH1D*>(inputFile->Get("hdx_coinanti"));
  hdx_anticoin = (TH1D*)(hdx_anticoin_raw->Clone("hdx_anticoin"));
  hdx_anticoin->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  TH1D *hdx_anticoin_clone = (TH1D*)(hdx_anticoin_raw->Clone("hdx_anticoin_clone"));
  hdx_anticoin_clone->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  if (!hdx_data_raw || !hdx_antidy_raw || !hdx_anticoin_raw) 
    handleError(inputFile,"hdx_data");

  //get histograms from MC plot file
  TFile* inputFileMC = new TFile(fmcinPath.c_str(), "READ");
  if (!inputFileMC || inputFileMC->IsZombie())
    handleError(inputFile,inputFileMC,"inputFileMC");

  //Get dxdy for spot check on MC
  TH2D *hdxdy_raw = dynamic_cast<TH2D*>(inputFileMC->Get("hdxdy"));
  TH2D *hdxdy = (TH2D*)(hdxdy_raw->Clone("hdxdy"));

  //fix the interpolate functions
  TH1D *hdx_p_raw = dynamic_cast<TH1D*>(inputFileMC->Get("hdx_p"));
  hdx_p = (TH1D*)(hdx_p_raw->Clone("hdx_p"));
  TH1D *hdx_p_clone = (TH1D*)(hdx_p_raw->Clone("hdx_p_clone"));
  hdx_p->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  hdx_p_wide = (TH1D*)(hdx_p_raw->Clone("hdx_p"));
  hdx_p_wide->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  TH1D *hdx_n_raw = dynamic_cast<TH1D*>(inputFileMC->Get("hdx_n"));
  hdx_n = (TH1D*)(hdx_n_raw->Clone("hdx_n"));
  TH1D *hdx_n_clone = (TH1D*)(hdx_n_raw->Clone("hdx_n_clone"));
  hdx_n->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  hdx_n_wide = (TH1D*)(hdx_n_raw->Clone("hdx_n"));
  hdx_n_wide->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  //get histograms from MC plot file
  TFile* inputInelFileMC = new TFile(fmcinelinPath.c_str(), "READ");
  
  TH1D *hdx_inel_raw = dynamic_cast<TH1D*>(inputInelFileMC->Get("hdx_inel"));
  hdx_inel = (TH1D*)(hdx_inel_raw->Clone("hdx_inel"));
  hdx_inel->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  TH1D *hdx_inel_clone = (TH1D*)(hdx_inel_raw->Clone("hdx_inel_clone"));
  hdx_inel_clone->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  //Handle missing histograms
  if (!hdx_p_raw || !hdx_n_raw || !hdx_inel_raw || !hdx_data_raw) 
    handleError(inputFile,inputFileMC,"hdx");

  hdx_inel->Scale(10e33); //account for lack of overall normalization from g4sbs generator

  TFile* outputFile = new TFile(foutPath.c_str(), "RECREATE");

  //test canvas
  TH1D *hdx_data_test = (TH1D*)(hdx_data_raw->Clone("hdx_data_test"));
  TH1D *hdx_p_test = (TH1D*)(hdx_p_raw->Clone("hdx_p_test"));
  TH1D *hdx_n_test = (TH1D*)(hdx_n_raw->Clone("hdx_n_test"));
  TH1D *hdx_inel_test = (TH1D*)(hdx_inel_raw->Clone("hdx_inel_test"));
  TH1D *hdx_antidy_test = (TH1D*)(hdx_antidy_raw->Clone("hdx_antidy_test"));
  TH1D *hdx_anticoin_test = (TH1D*)(hdx_anticoin_raw->Clone("hdx_anticoin_test"));

  // Normalize histograms to unity
  hdx_p_test->Scale(1.5 / hdx_p_test->Integral());
  hdx_n_test->Scale(0.5 / hdx_n_test->Integral());
  hdx_inel_test->Scale(1.0 / hdx_inel_test->Integral());
  hdx_antidy_test->Scale(1.0 / hdx_antidy_test->Integral());
  hdx_anticoin_test->Scale(1.0 / hdx_anticoin_test->Integral());
  hdx_data_test->Scale(3.0 / hdx_data_test->Integral());
    
  // Create a canvas to plot histograms
  TCanvas *canvas = new TCanvas("canvas", "Histograms Normalized to Unity", 800, 600);
    
  // Set histogram drawing options
  hdx_p_test->SetLineColor(kRed);
  hdx_n_test->SetLineColor(kBlue);
  hdx_inel_test->SetLineColor(kGreen);
  hdx_antidy_test->SetLineColor(kMagenta);
  hdx_anticoin_test->SetLineColor(kOrange);
  hdx_data_test->SetLineWidth(2);
  hdx_data_test->SetLineColor(kBlack);
    
  // Draw histograms on the same canvas
  hdx_p_test->GetYaxis()->SetRangeUser(0,0.04);
  hdx_p_test = util::shiftHistogramX(hdx_p_test, 0.2);
  hdx_p_test->SetTitle("dx;m");
  hdx_p_test->Draw("HIST");
  hdx_n_test = util::shiftHistogramX(hdx_n_test, 0.2);
  hdx_n_test->Draw("HIST SAME");
  hdx_inel_test->Draw("HIST SAME");
  hdx_antidy_test->DrawNormalized("HIST SAME");
  hdx_anticoin_test->DrawNormalized("HIST SAME");
  hdx_data_test->Draw("HIST SAME");
    
  // Add a legend
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(hdx_p_test, "simc proton", "l");
  legend->AddEntry(hdx_n_test, "simc neutron", "l");
  legend->AddEntry(hdx_inel_test, "inel bg distribution", "l");
  legend->AddEntry(hdx_antidy_test, "antidy bg", "l");
  legend->AddEntry(hdx_anticoin_test, "anticoin bg", "l");
  legend->AddEntry(hdx_data_test, "data", "l");
  legend->Draw();
    
  // Update the canvas
  canvas->Update();
    
  TCanvas *cdxdy = new TCanvas("cdxdy", "dx vs dy", 800, 600);
  cdxdy->cd();
  hdxdy->Draw();

  //Run through various fits with different background methods

  //For now, best fit to data. Use to get shifts for all other fits.
  std::pair<double,double> qualset;
  auto setpar_vector = fitAndFineFit(hdx_data_set, "pol2BG", "fit_pol2BG", 7, hcalfit_l, hcalfit_h, qualset, fitopt.c_str());

  cout << "Set vector: ";
  for (auto par : setpar_vector )
    cout << "(" << par.first << ", " << par.second << ") ";
  cout << endl << endl;

  //Get the same shift for all BG comparisons
  double pshift = setpar_vector[2].first;
  double nshift = setpar_vector[3].first;

  std::pair<double,double> qualp2;
  auto p2par_vector = fitAndFineFit_fixshift(hdx_data, "pol2BG", "fit_pol2BG", 7, hcalfit_l, hcalfit_h, qualp2, pshift, nshift, fitopt.c_str());
  
  cout << "Pol 2 vector: ";
  for (auto par : p2par_vector )
    cout << "(" << par.first << ", " << par.second << ") ";
  cout << endl << endl;

  //If a set is given as override arguments, take those
  if(pshift_override!=0)
    pshift = pshift_override;
  if(nshift_override!=0)
    nshift = nshift_override;

  //Redo the pol2 fit if overrides exist
  if(pshift_override!=0 || nshift_override!=0){
    p2par_vector = fitAndFineFit_fixshift(hdx_data, "pol2BG", "fit_pol2BG", 7, hcalfit_l, hcalfit_h, qualp2, pshift, nshift, fitopt.c_str());
    cout << "Overriding p and n shift parameters!" << endl;
  }

  //Pol1 fit for concavity comparisons only
  std::pair<double,double> qualp1;
  auto p1par_vector = fitAndFineFit_fixshift(hdx_data, "pol1BG", "fit_pol1BG", 6, hcalfit_l, hcalfit_h, qualp1, pshift, nshift, fitopt.c_str());

  cout << "Pol 1 vector: ";
  for (auto par : p1par_vector )
    cout << "(" << par.first << ", " << par.second << ") ";
  cout << endl << endl;

  //pol3 fit
  std::pair<double,double> qualp3;
  auto p3par_vector = fitAndFineFit_fixshift(hdx_data, "pol3BG", "fit_pol3BG", 8, hcalfit_l, hcalfit_h, qualp3, pshift, nshift, fitopt.c_str());

  cout << "Pol 3 vector: ";
  for (auto par : p3par_vector )
    cout << "(" << par.first << ", " << par.second << ") ";
  cout << endl << endl;

  //pol3 fit an pin (test)
  // Create a canvas
  TCanvas *c = new TCanvas("c", "Minuet test to pol3 endpoint constraint", 800, 600);
    
  // Draw the histogram
  hdx_data_pol3test->Draw("E");
    
  // Set up Minuit for constrained fit
  TMinuit *minuit = new TMinuit(8); // 4 parameters for proton+neutron + 4 for pol3
  minuit->SetFCN(chi2_constrained);
  minuit->SetObjectFit(hdx_data_pol3test);

  minuit->DefineParameter(0, "proton_scale", 1, 0.1, 0, 10);
  minuit->DefineParameter(1, "neutron_scale", 1, 0.1, 0, 10);
  minuit->DefineParameter(2, "dx_shift_p", 0, 0.01, -1, 1);
  minuit->DefineParameter(3, "dx_shift_n", 0, 0.01, -1, 1);
  minuit->DefineParameter(4, "a", 1, 0.1, -10, 10);
  minuit->DefineParameter(5, "b", 1, 0.1, -10, 10);
  minuit->DefineParameter(6, "c", 1, 0.1, -10, 10);
  minuit->DefineParameter(7, "d", 1, 0.1, -10, 10);
    
  // Perform the minimization
  minuit->Migrad();
    
  // Get the parameters
  double par[8];
  for (int i = 0; i < 8; ++i) {
    minuit->GetParameter(i, par[i], par[i]);
  }
    
  // Create and draw the fitted function
  TF1 *fitFunc = new TF1("fitFunc", fit_pol3BG, -3, 6, 8);
  fitFunc->SetParameters(par);
  fitFunc->SetLineColor(kRed);
  fitFunc->Draw("same");
    
  // Update the canvas
  c->Update();

  std::pair<double,double> qualp4;
  auto p4par_vector = fitAndFineFit_fixshift(hdx_data_wide, "pol4BG", "fit_pol4BG", 9, hcalfit_l, hcalfit_h, qualp4, pshift, nshift, fitopt.c_str());

  cout << "Pol 4 vector: ";
  for (auto par : p4par_vector )
    cout << "(" << par.first << ", " << par.second << ") ";
  cout << endl << endl;

  std::pair<double,double> qualgaus;
  auto gauspar_vector = fitAndFineFit_fixshift(hdx_data_wide, "gausBG", "fit_gausBG", 7, hcalfit_l-wideband_offset, hcalfit_h+wideband_offset, qualgaus, pshift, nshift, fitopt.c_str());

  cout << "Gaus vector: ";
  for (auto par : gauspar_vector )
    cout << "(" << par.first << ", " << par.second << ") ";
  cout << endl << endl;

  //This returns parameters for the sideband pol4 background only
  std::pair<double,double> qualsb;
  auto sbpar_vector = fitAndFineFit_fixshift(hdx_data_wide, "sbBG", "sbBGfit", 5, hcalfit_l-wideband_offset, hcalfit_h+wideband_offset, qualsb, pshift, nshift, fitopt.c_str());

  std::pair<double,double> qualinel;
  auto inelpar_vector = fitAndFineFit_fixshift(hdx_data, "inelBG", "fit_inelBG", 5, hcalfit_l, hcalfit_h, qualinel, pshift, nshift, fitopt.c_str());

  std::pair<double,double> qualantidy;
  auto antidypar_vector = fitAndFineFit_fixshift(hdx_data, "antidyBG", "fit_dyantiBG", 5, hcalfit_l, hcalfit_h, qualantidy, pshift, nshift, fitopt.c_str());

  std::pair<double,double> qualanticoin;
  auto anticoinpar_vector = fitAndFineFit_fixshift(hdx_data, "anticoinBG", "fit_coinantiBG", 5, hcalfit_l, hcalfit_h, qualanticoin, pshift, nshift, fitopt.c_str());

  // Make background functions
  TF1 *bg_gausFit = new TF1("bg_gausFit",fits::g_gfit,hcalfit_l,hcalfit_h,3);
  for (int i=0; i<3; ++i){
    bg_gausFit->SetParameter(i,gauspar_vector[i+4].first);
    bg_gausFit->SetParError(i,gauspar_vector[i+4].second);
  }

  TF1 *bg_p4Fit = new TF1("bg_pol4Fit",fits::g_p4fit,hcalfit_l,hcalfit_h,5);
  for (int i=0; i<5; ++i){
    bg_p4Fit->SetParameter(i,p4par_vector[i+4].first);
    bg_p4Fit->SetParError(i,p4par_vector[i+4].second);
  }

  TF1 *bg_p2Fit = new TF1("bg_pol2Fit",fits::g_p2fit_cd,hcalfit_l,hcalfit_h,3);
  for (int i=0; i<3; ++i){
    bg_p2Fit->SetParameter(i,p2par_vector[i+4].first);
    bg_p2Fit->SetParError(i,p2par_vector[i+4].second);
  }

  TF1 *bg_p1Fit = new TF1("bg_pol1Fit",fits::g_p1fit,hcalfit_l,hcalfit_h,2);
  for (int i=0; i<2; ++i){
    bg_p1Fit->SetParameter(i,p1par_vector[i+4].first);
    bg_p1Fit->SetParError(i,p1par_vector[i+4].second);
  }

  TF1 *bg_p3Fit = new TF1("bg_pol3Fit",fits::g_p3fit,hcalfit_l,hcalfit_h,4);
  for (int i=0; i<2; ++i){
    bg_p3Fit->SetParameter(i,p3par_vector[i+4].first);
    bg_p3Fit->SetParError(i,p3par_vector[i+4].second);
  }

  TF1 *bg_sbFit = new TF1("bg_sbFit",fits::g_p4fit,hcalfit_l,hcalfit_h,5);
  for (int i=0; i<5; ++i){
    bg_sbFit->SetParameter(i,sbpar_vector[i].first);
    bg_sbFit->SetParError(i,sbpar_vector[i].second);
  }

  std::vector<TF1*> bgFunctions = {bg_gausFit, bg_p4Fit, bg_p2Fit, bg_p1Fit, bg_p3Fit, bg_sbFit};
  plotBackgrounds(bgFunctions);

  //Print the canvas
  //pol4 fit to bg
  double elastics_gaus;
  std::pair<double,double> rsf_gaus;
  createAndSaveCanvasWithFitsAndResiduals(hdx_data, hdx_p_clone, hdx_n_clone, bg_gausFit, gauspar_vector, qualgaus, "gaus BG", "~/gmn_plots/new_gausfit.pdf", -1, elastics_gaus, rsf_gaus, kine, mag);

  //pol4 fit to bg
  double elastics_p4;
  std::pair<double,double> rsf_p4;
  createAndSaveCanvasWithFitsAndResiduals(hdx_data, hdx_p_clone, hdx_n_clone, bg_p4Fit, p4par_vector, qualp4, "pol4 BG", "~/gmn_plots/new_p4fit.pdf", 4, elastics_p4, rsf_p4, kine, mag);

  //pol2 fit to bg
  double elastics_p2;
  std::pair<double,double> rsf_p2;
  createAndSaveCanvasWithFitsAndResiduals(hdx_data, hdx_p_clone, hdx_n_clone, bg_p2Fit, p2par_vector, qualp2, "pol2 BG", "~/gmn_plots/new_p2fit.pdf", 2, elastics_p2, rsf_p2, kine, mag);

  //pol1 fit to bg
  double elastics_p1;
  std::pair<double,double> rsf_p1;
  //createAndSaveCanvasWithFitsAndResiduals(hdx_data, hdx_p_clone, hdx_n_clone, bg_p1Fit, p1par_vector, qualp1, "pol1 BG", "~/gmn_plots/new_p1fit.pdf", 1, elastics_p1, rsf_p1, kine, mag);

  //pol3 fit to bg
  double elastics_p3;
  std::pair<double,double> rsf_p3;
  createAndSaveCanvasWithFitsAndResiduals(hdx_data, hdx_p_clone, hdx_n_clone, bg_p3Fit, p3par_vector, qualp3, "pol3 BG", "~/gmn_plots/new_p3fit.pdf", 3, elastics_p3, rsf_p3, kine, mag);

  // //inel bg
  double elastics_inel;
  std::pair<double,double> rsf_inel;
  createAndSaveCanvasWithFitsAndResiduals_altbg(hdx_data, hdx_p_clone, hdx_n_clone, hdx_inel_clone, inelpar_vector, qualinel, "inel BG", "~/gmn_plots/new_inelfit.pdf", "inel", elastics_inel, rsf_inel, kine, mag);
  
  //dy anticut bg
  double elastics_antidy;
  std::pair<double,double> rsf_antidy;
  createAndSaveCanvasWithFitsAndResiduals_altbg(hdx_data, hdx_p_clone, hdx_n_clone, hdx_antidy_clone, antidypar_vector, qualantidy, "antidy BG", "~/gmn_plots/new_antidy.pdf", "dy", elastics_antidy, rsf_antidy, kine, mag);
  
  // //coin anticut bg
  double elastics_anticoin;
  std::pair<double,double> rsf_anticoin;
  createAndSaveCanvasWithFitsAndResiduals_altbg(hdx_data, hdx_p_clone, hdx_n_clone, hdx_anticoin_clone, anticoinpar_vector, qualanticoin, "anticoin BG", "~/gmn_plots/new_anticoin.pdf", "coin", elastics_anticoin, rsf_anticoin, kine, mag);

  //subtract bg fits
  subtractBackground(hdx_sb_nobg,bg_sbFit);

  std::pair<double,double> qualnobg;
  auto sbpar_nobg_vector = fitAndFineFit(hdx_sb_nobg, "noBG", "fit_noBG", 4, hcalfit_l, hcalfit_h, qualnobg, fitopt.c_str());
  
  TCanvas *ctest = new TCanvas("ctest","test canvas",1200,800);
  ctest->cd();
  hdx_data_sb->Draw();
  bg_sbFit->Draw("same");
  hdx_sb_nobg->Draw("same");

  double elastics_sb_nobg;
  std::pair<double,double> rsf_nobg;
  createAndSaveCanvasWithFitsAndResiduals_nobg(hdx_sb_nobg, hdx_p_clone, hdx_n_clone, sbpar_nobg_vector, qualnobg, "Sideband Order-4 Poly BG subtracted", "~/gmn_plots/new_sbnobg.pdf", elastics_sb_nobg, rsf_nobg, kine, mag);


  ///////////////////
  ///////////////////
  ///////////////////
  ///////////////////
  ///////////////////
  ///////////////////
  ///////////////////
  ///////////////////
  //Get Error estimates

  cout << endl << "ALL RSF by BG:" << endl;
  cout << "Antidy: " << rsf_antidy.first << " pm " << rsf_antidy.second << endl;
  cout << "Anticoin: " << rsf_anticoin.first << " pm " << rsf_anticoin.second << endl;
  cout << "pol2: " << rsf_p2.first << " pm " << rsf_p2.second << endl;
  cout << "pol3: " << rsf_p3.first << " pm " << rsf_p3.second << endl;
  cout << "pol4: " << rsf_p4.first << " pm " << rsf_p4.second << endl;
  cout << "Gaussian: " << rsf_gaus.first << " pm " << rsf_gaus.second << endl;
  cout << "Inelastic: " << rsf_inel.first << " pm " << rsf_inel.second << endl;

  cout << endl << endl;
  
  
  //get arb error from analytical bg fits
  //double Rsf_reported = (rsf_gaus.first + rsf_p4.first + rsf_p2.first) / 3;
  //double Rsf_reported = (rsf_gaus.first + rsf_p4.first + rsf_p2.first) / 3;

  //double Rsf_reported = (rsf_gaus.first + rsf_p3.first + rsf_anticoin.first + rsf_p2.first + rsf_antidy.first) / 5;

  double Rsf_reported = (rsf_gaus.first + rsf_anticoin.first + rsf_p2.first + rsf_antidy.first) / 4;

  //double Rsf_reported = (rsf_anticoin.first + rsf_p2.first + rsf_antidy.first) / 3;

  // Calculate the standard deviation of the first elements
  // double variance = (std::pow(rsf_gaus.first - Rsf_reported, 2) +
  // 		     std::pow(rsf_p4.first - Rsf_reported, 2) +
  // 		     std::pow(rsf_p2.first - Rsf_reported, 2)) / 3;

  // double variance = (std::pow(rsf_gaus.first - Rsf_reported, 2) +
  // 		     std::pow(rsf_anticoin.first - Rsf_reported, 2) +
  // 		     std::pow(rsf_antidy.first - Rsf_reported, 2) +
  // 		     std::pow(rsf_p3.first - Rsf_reported, 2) +
  // 		     std::pow(rsf_p2.first - Rsf_reported, 2)) / 5;

  double variance = (std::pow(rsf_gaus.first - Rsf_reported, 2) +
		     std::pow(rsf_anticoin.first - Rsf_reported, 2) +
		     std::pow(rsf_antidy.first - Rsf_reported, 2) +
		     std::pow(rsf_p2.first - Rsf_reported, 2)) / 4;
  
  // double variance = (std::pow(rsf_anticoin.first - Rsf_reported, 2) +
  // 		     std::pow(rsf_antidy.first - Rsf_reported, 2) +
  // 		     std::pow(rsf_p2.first - Rsf_reported, 2)) / 3;

  
  double systerr_inel = std::sqrt(variance);

  double staterr = rsf_p2.second;

  double systerr_hde = Rsf_reported*sqrt(2)*(HDE_systerr_N/HDE_syst_N); //assuming p and n have similar HDE error. Follows from (dRsf/Rsf)^2=(dNn/Nn)^2+(dNp/Np)^2 and dN/dHDE = N/HDE 

  double systerr = sqrt(pow(systerr, 2) + pow(systerr_hde, 2)); //std quadrature

  cout << endl << endl << "Rsf error budget (syst_inel, syst_hde, syst_total, stat): " << systerr_inel << ", " << systerr_hde << ", " << systerr << ", " << staterr << endl << endl;

  cout << Rsf_reported << " " << systerr << " " << staterr << " " << systerr_hde << endl;

  auto GMn_and_error_budget = extract::extract_GMn_from_simc( Rsf_reported,
  							      staterr,
							      systerr,
  							      0.,          
  							      E_beam,
  							      BB_angle,
  							      true);

  auto antidy_ex = extract::extract_GMn_from_simc( rsf_antidy.first,
						   rsf_antidy.second,
						   0.,
						   0.,          
						   E_beam,
						   BB_angle,
						   true);


  auto nobg_ex = extract::extract_GMn_from_simc( rsf_nobg.first,
						   rsf_nobg.second,
						   0.,
						   0.,          
						   E_beam,
						   BB_angle,
						   true);

  cout << endl << endl << "//////////////////////////" << endl;


  auto& nobg_ex_tuple = nobg_ex[1]; // Access the tuple for corrected GMn

  // Extracting values from the tuple for clarity
  double nobg_ex_gmn = std::get<0>(nobg_ex_tuple);
  double nobg_ex_stat = std::get<1>(nobg_ex_tuple);
  double nobg_ex_syst = std::get<2>(nobg_ex_tuple);
  double nobg_ex_model = std::get<3>(nobg_ex_tuple);
  double nobg_ex_total = std::get<4>(nobg_ex_tuple);

  // Output the values
  std::cout << "Corrected GMn (sideband only): " << nobg_ex_gmn << " \n"
	    << "Statistical Error (sideband only): " << nobg_ex_stat << " \n"
	    << "Systematic Error (sideband only): " << nobg_ex_syst << " \n"
	    << "Model Error (sideband only): " << nobg_ex_model << " \n"
	    << "Total Error (sideband only): " << nobg_ex_total << " \n";

  cout << endl << endl << "//////////////////////////" << endl;


  auto& antidy_ex_tuple = antidy_ex[1]; // Access the tuple for corrected GMn

  // Extracting values from the tuple for clarity
  double antidy_ex_gmn = std::get<0>(antidy_ex_tuple);
  double antidy_ex_stat = std::get<1>(antidy_ex_tuple);
  double antidy_ex_syst = std::get<2>(antidy_ex_tuple);
  double antidy_ex_model = std::get<3>(antidy_ex_tuple);
  double antidy_ex_total = std::get<4>(antidy_ex_tuple);

  // Output the values
  std::cout << "Corrected GMn (antidy only): " << antidy_ex_gmn << " \n"
	    << "Statistical Error (antidy only): " << antidy_ex_stat << " \n"
	    << "Systematic Error (antidy only): " << antidy_ex_syst << " \n"
	    << "Model Error (antidy only): " << antidy_ex_model << " \n"
	    << "Total Error (antidy only): " << antidy_ex_total << " \n";

  cout << endl << endl << "//////////////////////////" << endl;

  // Check if the vector has the expected number of elements to avoid out-of-range errors
  auto& corrected_GMn_tuple = GMn_and_error_budget[1]; // Access the tuple for corrected GMn

  // Extracting values from the tuple for clarity
  double corrected_GMn = std::get<0>(corrected_GMn_tuple);
  double error_stat = std::get<1>(corrected_GMn_tuple);
  double error_syst = std::get<2>(corrected_GMn_tuple);
  double error_model = std::get<3>(corrected_GMn_tuple);
  double error_total = std::get<4>(corrected_GMn_tuple);

  // Output the values
  std::cout << "Corrected GMn: " << corrected_GMn << " \n"
	    << "Statistical Error: " << error_stat << " \n"
	    << "Systematic Error: " << error_syst << " \n"
	    << "Model Error: " << error_model << " \n"
	    << "Total Error: " << error_total << " \n";


  // Create a canvas for displaying parameters
  TCanvas* paramCanvas = new TCanvas("paramCanvas", "GMn and Error", 800, 600);
  paramCanvas->cd();
  TPaveText* paramText = new TPaveText(0.1, 0.1, 0.9, 0.9); // Normalized coordinates

  // Add a header
  paramText->AddText("GMn and Error:");
  paramText->AddText(" ");
  //paramText->AddLine();
    
  // Adding each parameter to the TPaveText
  paramText->AddText(Form("Corrected GMn: %0.3f ", corrected_GMn));
  paramText->AddText(Form("Statistical Error: %0.3f ", error_stat));
  paramText->AddText(Form("Systematic Error: %0.3f ", error_syst));
  paramText->AddText(Form("Model Error: %0.3f ", error_model));
  paramText->AddText(Form("Total Error: %0.3f ", error_total));

  // Set text alignment and fill color
  //paramText->SetTextAlign(12); // Align text to left
  //paramText->SetFillColor(0);  // Transparent background

  // Draw the TPaveText on the canvas
  paramText->Draw();

  cout << endl << endl;

  // // Create a canvas for displaying cuts
  // TCanvas* cutCanvas = new TCanvas("cutCanvas", "Cuts", 800, 600);
  // cutCanvas->cd();
  // TPaveText* cutsText = new TPaveText(0.1, 0.1, 0.9, 0.9); // coordinates are normalized

  // cutsText->AddText("Cuts Used:");
  // cutsText->AddLine();
  // for (const auto& cut : cuts) {
  //   cutsText->AddText(cut.c_str());
  // }

  // cutsText->SetTextAlign(12); // Align text to left
  // cutsText->SetFillColor(0);  // Transparent background
  // cutsText->Draw();

  // // Write the canvas to the output file
  // cutCanvas->Write();
  
  // Create a canvas for the itemized cuts from gcut_vect
  TCanvas *cCuts = new TCanvas("cCuts", "Itemized Cuts", 1200, 800);
  cCuts->cd();
  TLatex latex;
  latex.SetTextSize(0.04);
  latex.DrawLatexNDC(0.1, 0.9, "Itemized Cuts:");

  for (size_t j = 0; j < cuts.size(); ++j) {
    latex.DrawLatexNDC(0.1, 0.9 - 0.07 * (j + 1), cuts[j].c_str());
  }

  //latex.DrawLatexNDC(0.1, 0.9 - 0.07 * (gcut_vect.size() + 1), " ");
  //latex.DrawLatexNDC(0.1, 0.9 - 0.07 * (gcut_vect.size() + 2), "MC weights are applied to each histogram.");

  cCuts->Update();
  cCuts->Write();


  outputFile->Write();
  //outputFile->Close();

  cout << "All plots created. Output file located here: " << foutPath << endl;

}

/////////////
/////////////
//functions
/////////////
/////////////
/////////////


void handleError(TFile *file1, TFile *file2, std::string marker) {
    if (file1) file1->Close();
    if (file2) file2->Close();
    std::cerr << "Error: File opening or histogram retrieval failed at " << marker << "." << std::endl;
}

void handleError(TFile *file1, std::string marker) {
    if (file1) file1->Close();
    std::cerr << "Error: File opening or histogram retrieval failed at " << marker << "." << std::endl;
}

void createAndSaveCanvasWithFitsAndResiduals(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, TF1* bg, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, int pol, double &elastics, std::pair<double,double> &rsf, int kine, int mag) {

  // Generate a unique identifier for the fit function name
  TString uniqueID = TString::Format("_%u", gRandom->Integer(1000000));

  //Clone the dx histo for bg sub
  TH1D *hdx_histplot = (TH1D*)(hdx->Clone("hdx_histplot"));
  TH1D *hdx_clone = (TH1D*)(hdx->Clone("hdx_clone"));
  TH1D *hdx_p_clone = (TH1D*)(hdxp->Clone("hdx_p_clone"));
  TH1D *hdx_n_clone = (TH1D*)(hdxn->Clone("hdx_n_clone"));

  //Recreate the fit
  TF1* fit;
  std::string bgword = "";
  if (pol == 4){
    fit = new TF1(TString::Format("fitp4%s", uniqueID.Data()), fit_pol4BG, hcalfit_l - wideband_offset, hcalfit_h + wideband_offset, 9);
    bgword = "Order-4 Polynomial ";
  }else if (pol == 2){
    fit = new TF1(TString::Format("fitp2%s", uniqueID.Data()), fit_pol2BG, hcalfit_l, hcalfit_h, 7);
    bgword = "Order-2 Polynomial ";
  }else if (pol == 1){
    fit = new TF1(TString::Format("fitp1%s", uniqueID.Data()), fit_pol1BG, hcalfit_l, hcalfit_h, 6);
    bgword = "Order-1 Polynomial ";
  }else if (pol == 3){
    fit = new TF1(TString::Format("fitp3%s", uniqueID.Data()), fit_pol3BG, hcalfit_l, hcalfit_h, 8);
    bgword = "Order-3 Polynomial ";
  }else if (pol == -1){
    fit = new TF1(TString::Format("fitgaus%s", uniqueID.Data()), fit_gausBG, hcalfit_l, hcalfit_h, 7);
    bgword = "Gaussian ";
  }

  for (size_t i = 0; i < pars.size(); ++i) {
    fit->SetParameter(i, pars[i].first);
    fit->SetParError(i, pars[i].second);
  }

  std::pair<double,double> altchisqrndf = calculateChiSquareWithZeroMCError(hdx, fit);
  
  //fit->SetChisquare(altchisqrndf.first);
  //fit->SetNDF(altchisqrndf.second);

  fit->SetChisquare(qual.first);
  fit->SetNDF(qual.second);

  // Create a canvas with two pads
  TCanvas *cTotal = new TCanvas(savePath, Form("Data with Fits and Residuals %s",type), 1200, 800);
  TPad *pad1 = new TPad("pad1", "Pad with the fit", 0.0, 0.3, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "Pad with the residuals", 0.0, 0.0, 1.0, 0.3);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad1->Draw();
  pad2->Draw();

  // First pad: Draw the histogram with the fits
  pad1->cd();

  gPad->SetGridx();
  gPad->SetGridy();

  //hdx_histplot->SetLineColor(kBlack);
  //hdx_histplot->SetLineWidth(2);
  hdx_histplot->SetTitle(Form("dx, %s (SBS-%d, %d%% field);m", type, kine, mag));
  hdx_histplot->SetTitleFont(132);
  hdx_histplot->GetXaxis()->SetTitleFont(132);
  hdx_histplot->GetXaxis()->SetLabelFont(132);
  hdx_histplot->GetYaxis()->SetTitleFont(132);
  hdx_histplot->GetYaxis()->SetLabelFont(132);
  hdx_histplot->SetMarkerStyle(21);
  hdx_histplot->SetMarkerSize(0.5);
  hdx_histplot->Draw("P");
  hdx->SetLineColor(kBlack);
  hdx->SetLineWidth(1);
  hdx->Draw("E same");
  fit->SetLineColor(kGreen);
  fit->SetLineWidth(1);
  fit->SetFillColorAlpha(kGreen,0.35);
  fit->SetFillStyle(1001);
  fit->SetNpx(1000);
  fit->Draw("sameFC");

  hdx_p_clone = util::shiftHistogramX(hdx_p, pars[2].first);
  hdx_p_clone->Scale(pars[0].first);
  hdx_p_clone->SetLineColor(kRed);
  hdx_p_clone->SetLineWidth(2);
  hdx_p_clone->SetMarkerStyle(20);  // Filled circle marker
  hdx_p_clone->SetMarkerColor(kRed);
  hdx_p_clone->SetMarkerSize(0.5);
  hdx_p_clone->Draw("P same");

  hdx_n_clone = util::shiftHistogramX(hdx_n, pars[3].first);
  hdx_n_clone->Scale(pars[1].first);
  hdx_n_clone->SetLineColor(kBlue);
  hdx_n_clone->SetLineWidth(2);
  hdx_n_clone->SetMarkerStyle(20);  // Filled circle marker
  hdx_n_clone->SetMarkerColor(kBlue);
  hdx_n_clone->SetMarkerSize(0.5);
  hdx_n_clone->Draw("P same");

  // Set background function
  bg->SetLineColor(kBlack);
  bg->SetFillColorAlpha(kBlack,0.35);
  bg->SetFillStyle(3013);
  bg->Draw("same");

  double psumerror;
  double psum = hdx_p_clone->IntegralAndError(0, hdx_p_clone->GetNbinsX() + 1, psumerror, "");

  double pscaleerror = pars[0].second;
  double pscale = pars[0].first;

  double nsumerror;
  double nsum = hdx_n_clone->IntegralAndError(0, hdx_n_clone->GetNbinsX() + 1, nsumerror, "");

  double nscaleerror = pars[1].second;
  double nscale = pars[1].first;

  double np_sum_ratio = nsum/psum;
  double np_par_ratio = nscale/pscale;
  double np_sum_ratio_error = np_sum_ratio * sqrt(pow(psumerror/psum, 2) + pow(nsumerror/nsum, 2));
  double np_par_ratio_error = np_par_ratio * sqrt(pow(pscaleerror/pscale, 2) + pow(nscaleerror/nscale, 2));

  rsf.first = np_par_ratio;
  rsf.second = np_par_ratio_error;

  // Construct legend within function
  TLegend* leg = new TLegend(0.5, 0.5, 0.89, 0.89);
  leg->SetBorderSize(0); // Remove border
  //leg->SetFillStyle(0);  // Make background transparent
  leg->SetTextColor(kBlack);
  leg->SetTextFont(132);
  leg->SetTextSize(0.05);
  leg->AddEntry(hdx_histplot, "Data", "p");
  leg->AddEntry(fit, "Fit (Proton + Neutron + BG)", "f");
  leg->AddEntry(bg, Form("%s Background",bgword.c_str()), "f");
  leg->AddEntry(hdx_p_clone, "Proton SIMC MC", "p");
  leg->AddEntry(hdx_n_clone, "Neutron SIMC MC", "p");
  leg->AddEntry( (TObject*)0, "", "");
  leg->AddEntry( (TObject*)0, Form("n/p scale ratio R_{sf} : %0.3f #pm %0.3f",np_par_ratio, np_par_ratio_error), "");
  leg->AddEntry( (TObject*)0, "", "");
  //leg->AddEntry( (TObject*)0, Form("data N events : %0.0f",hdx->GetEntries()), "");
  //leg->AddEntry( (TObject*)0, Form("n/p yield ratio R'' : %0.3f #pm %0.3f",np_sum_ratio, np_sum_ratio_error), "");
  //leg->AddEntry( (TObject*)0, Form("dx shift pars, n/p : %0.3f / %0.3f ", pars[3].first, pars[2].first), "");
  leg->AddEntry( (TObject*)0, Form("#chi^{2}/ndf: %0.3f/%d",fit->GetChisquare(),fit->GetNDF()), "");
  leg->Draw("same");

  //Background subraction on dx_clone and total elastic statistics
  elastics = 0;
  for (int bin = 1; bin <= hdx->GetNbinsX(); ++bin) {
    double bgValue = bg->Eval(hdx_clone->GetXaxis()->GetBinCenter(bin));
    if( hdx_clone->GetBinContent(bin) > bgValue )
      elastics += hdx_clone->GetBinContent(bin) - bgValue;
    hdx_clone->SetBinContent(bin, hdx_clone->GetBinContent(bin) - bgValue);
    hdx_clone->SetBinError(bin, hdx_clone->GetBinError(bin));
  }

  TH1D *hRes = calculateResiduals(hdx, fit, Form("hResiduals_%s",uniqueID.Data()));
  TH1D *hRes_histplot = (TH1D*)(hRes->Clone(Form("hRes_histplot_%s",type)));

  hRes_histplot->SetTitle("");
  hRes_histplot->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  double dxMaxValue = hdx_clone->GetMaximum();
  hRes_histplot->GetYaxis()->SetRangeUser(-dxMaxValue/4,dxMaxValue/4);

  hRes_histplot->SetLineColor(kBlack);
  hRes_histplot->SetLineWidth(1);
  //hRes_histplot->SetFillStyle(1001);
  //hRes_histplot->SetFillColorAlpha(kRed,0.75);

  // Second pad: Draw the residuals
  pad2->cd();
  gPad->SetGridx();

  hRes_histplot->GetYaxis()->SetTitle("Residuals");
  hRes_histplot->GetYaxis()->SetTitleSize(0.1);
  hRes_histplot->GetYaxis()->SetTitleOffset(0.5);
  hRes_histplot->GetYaxis()->SetLabelSize(0.1);
  hRes_histplot->GetXaxis()->SetTitleSize(0.1);
  hRes_histplot->GetXaxis()->SetLabelSize(0.1);
  hRes_histplot->GetXaxis()->SetTitle("x_{hcal}-x_{exp}");
  hRes_histplot->GetXaxis()->SetTitleFont(132);
  hRes_histplot->GetXaxis()->SetLabelFont(132);
  hRes_histplot->GetYaxis()->SetTitleFont(132);
  hRes_histplot->GetYaxis()->SetLabelFont(132);
  hRes_histplot->Draw("E");

  // Get the axis range
  double y_min = hRes_histplot->GetMinimum();
  double y_max = hRes_histplot->GetMaximum();

  // Create a TPave to fill the region between the axes
  TPave *pave = new TPave(hcalfit_l, y_min, hcalfit_h, y_max, 0, "NB");
  pave->SetFillColorAlpha(kGray,0.4); // Light grey RGB
  pave->SetFillStyle(3001); // 3001: transparent fill style
  pave->Draw("same");

  // Draw a line at y=0 for residuals
  TF1 *zeroLine = new TF1("zeroLine", "0", hRes->GetXaxis()->GetXmin(), hRes->GetXaxis()->GetXmax());
  zeroLine->SetLineColor(kRed);
  zeroLine->SetLineStyle(2);
  zeroLine->Draw("same");

  // Redraw the histogram axes and title over the TPave
  hRes_histplot->Draw("AXIS SAME");

  // Update the canvas
  cTotal->cd();
  cTotal->Update();

  // Save the canvas to a file
  cTotal->SaveAs(savePath);
}

//Create canvas for reporting
void createAndSaveCanvasWithFitsAndResiduals_nobg(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, double &elastics, std::pair<double,double> &rsf, int kine, int mag) {

  // Generate a unique identifier for the fit function name
  TString uniqueID = TString::Format("_%u", gRandom->Integer(1000000));

  //Clone the dx histos
  TH1D *hdx_histplot = (TH1D*)(hdx->Clone("hdx_histplot"));
  TH1D *hdx_clone = (TH1D*)(hdx->Clone("hdx_clone"));
  TH1D *hdx_p_clone = (TH1D*)(hdxp->Clone("hdx_p_clone"));
  TH1D *hdx_n_clone = (TH1D*)(hdxn->Clone("hdx_n_clone"));

  //Recreate the fit
  TF1* fit = new TF1("fit_noBG", fit_noBG, hcalfit_l-wideband_offset, hcalfit_h+wideband_offset, 4);
  
  for (size_t i = 0; i < pars.size(); ++i) {
    fit->SetParameter(i, pars[i].first);
    fit->SetParError(i, pars[i].second);
  }

  std::pair<double,double> altchisqrndf = calculateChiSquareWithZeroMCError(hdx, fit);
  
  //fit->SetChisquare(altchisqrndf.first);
  //fit->SetNDF(altchisqrndf.second);

  fit->SetChisquare(qual.first);
  fit->SetNDF(qual.second);

  // Create a canvas with two pads
  TCanvas *cTotal = new TCanvas(savePath, Form("No BG Data with Fits and Residuals %s",type), 1200, 800);
  TPad *pad1 = new TPad("pad1", "Pad with the fit", 0.0, 0.3, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "Pad with the residuals", 0.0, 0.0, 1.0, 0.3);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad1->Draw();
  pad2->Draw();

  // First pad: Draw the histogram with the fits
  pad1->cd();

  gPad->SetGridx();
  gPad->SetGridy();

  //hdx_histplot->SetLineColor(kBlack);
  //hdx_histplot->SetLineWidth(2);
  hdx_histplot->SetTitle(Form("dx, %s (SBS-%d, %d%% field);m", type, kine, mag));
  hdx_histplot->SetTitleFont(132);
  hdx_histplot->GetXaxis()->SetTitleFont(132);
  hdx_histplot->GetXaxis()->SetLabelFont(132);
  hdx_histplot->GetYaxis()->SetTitleFont(132);
  hdx_histplot->GetYaxis()->SetLabelFont(132);
  hdx_histplot->SetMarkerStyle(21);
  hdx_histplot->SetMarkerSize(0.5);
  hdx_histplot->Draw("P");
  hdx->SetLineColor(kBlack);
  hdx->SetLineWidth(1);
  hdx->Draw("E same");
  fit->SetNpx(1000);
  fit->SetLineColor(kGreen);
  fit->SetLineWidth(1);
  fit->SetFillColorAlpha(kGreen,0.35);
  fit->SetFillStyle(1001);
  fit->Draw("sameFC");

  hdx_p_clone = util::shiftHistogramX(hdx_p, pars[2].first);
  hdx_p_clone->Scale(pars[0].first);
  hdx_p_clone->SetLineColor(kRed);
  hdx_p_clone->SetLineWidth(2);
  hdx_p_clone->SetMarkerStyle(20);  // Filled circle marker
  hdx_p_clone->SetMarkerColor(kRed);
  hdx_p_clone->SetMarkerSize(0.5);
  hdx_p_clone->Draw("P same");

  hdx_n_clone = util::shiftHistogramX(hdx_n, pars[3].first);
  hdx_n_clone->Scale(pars[1].first);
  hdx_n_clone->SetLineColor(kBlue);
  hdx_n_clone->SetLineWidth(2);
  hdx_n_clone->SetMarkerStyle(20);  // Filled circle marker
  hdx_n_clone->SetMarkerColor(kBlue);
  hdx_n_clone->SetMarkerSize(0.5);
  hdx_n_clone->Draw("P same");

  double psumerror;
  double psum = hdx_p_clone->IntegralAndError(0, hdx_p_clone->GetNbinsX() + 1, psumerror, "");

  double pscaleerror = pars[0].second;
  double pscale = pars[0].first;

  double nsumerror;
  double nsum = hdx_n_clone->IntegralAndError(0, hdx_n_clone->GetNbinsX() + 1, nsumerror, "");

  double nscaleerror = pars[1].second;
  double nscale = pars[1].first;

  double np_sum_ratio = nsum/psum;
  double np_par_ratio = nscale/pscale;
  double np_sum_ratio_error = np_sum_ratio * sqrt(pow(psumerror/psum, 2) + pow(nsumerror/nsum, 2));
  double np_par_ratio_error = np_par_ratio * sqrt(pow(pscaleerror/pscale, 2) + pow(nscaleerror/nscale, 2));

  rsf.first = np_par_ratio;
  rsf.second = np_par_ratio_error;

  // Construct legend within function
  TLegend* leg = new TLegend(0.45, 0.5, 0.89, 0.89);
  leg->SetBorderSize(0); // Remove border
  //leg->SetFillStyle(0);  // Make background transparent
  leg->SetTextColor(kBlack);
  leg->SetTextFont(132);
  leg->SetTextSize(0.04);
  leg->AddEntry(hdx, "Data", "p");
  leg->AddEntry(fit, "Fit (Proton + Neutron)", "f");
  leg->AddEntry(hdx_p_clone, "Proton SIMC MC", "p");
  leg->AddEntry(hdx_n_clone, "Neutron SIMC MC", "p");
  leg->AddEntry( (TObject*)0, "", "");
  leg->AddEntry( (TObject*)0, Form("n/p scale ratio R_{sf} : %0.3f #pm %0.3f",np_par_ratio, np_par_ratio_error), "");
  leg->AddEntry( (TObject*)0, "", "");
  //leg->AddEntry( (TObject*)0, Form("data N events (bg sub) : %0.0f",hdx->GetEntries()), "");
  //leg->AddEntry( (TObject*)0, Form("n/p yield ratio R'' : %0.3f #pm %0.3f",np_sum_ratio, np_sum_ratio_error), "");
  //leg->AddEntry( (TObject*)0, Form("dx shift pars, n/p : %0.3f / %0.3f ", pars[3].first, pars[2].first), "");
  leg->AddEntry( (TObject*)0, Form("#chi^{2}/ndf: %0.3f/%d",fit->GetChisquare(),fit->GetNDF()), "");
  leg->Draw("same");

  //Total elastic statistics
  elastics = 0;
  for (int bin = 1; bin <= hdx->GetNbinsX(); ++bin)
    elastics += hdx_clone->GetBinContent(bin);

  TH1D *hRes = calculateResiduals(hdx_clone, fit, Form("hResiduals_%s",uniqueID.Data()));
  TH1D *hRes_histplot = (TH1D*)(hRes->Clone(Form("hRes_histplot_%s",type)));

  hRes_histplot->SetTitle("");
  hRes_histplot->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  double dxMaxValue = hdx_clone->GetMaximum();
  hRes_histplot->GetYaxis()->SetRangeUser(-dxMaxValue/4,dxMaxValue/4);

  hRes_histplot->SetLineColor(kBlack);
  hRes_histplot->SetLineWidth(1);
  //hRes_histplot->SetFillStyle(1001);
  //hRes_histplot->SetFillColorAlpha(kRed,0.75);

  // Second pad: Draw the residuals
  pad2->cd();
  gPad->SetGridx();

  hRes_histplot->GetYaxis()->SetTitle("Residuals");
  hRes_histplot->GetYaxis()->SetTitleSize(0.1);
  hRes_histplot->GetYaxis()->SetTitleOffset(0.5);
  hRes_histplot->GetYaxis()->SetLabelSize(0.1);
  hRes_histplot->GetXaxis()->SetTitleSize(0.1);
  hRes_histplot->GetXaxis()->SetLabelSize(0.1);
  hRes_histplot->GetXaxis()->SetTitle("x_{hcal}-x_{exp}");
  hRes_histplot->GetXaxis()->SetTitleFont(132);
  hRes_histplot->GetXaxis()->SetLabelFont(132);
  hRes_histplot->GetYaxis()->SetTitleFont(132);
  hRes_histplot->GetYaxis()->SetLabelFont(132);
  hRes_histplot->Draw("E");

  // Get the axis range
  double y_min = hRes_histplot->GetMinimum();
  double y_max = hRes_histplot->GetMaximum();

  // Create a TPave to fill the region between the axes
  TPave *pave = new TPave(hcalfit_l, y_min, hcalfit_h, y_max, 0, "NB");
  pave->SetFillColorAlpha(kGray,0.4); // Light grey RGB
  pave->SetFillStyle(3001); // 3001: transparent fill style
  pave->Draw("same");

  // Draw a line at y=0 for residuals
  TF1 *zeroLine = new TF1("zeroLine", "0", hRes->GetXaxis()->GetXmin(), hRes->GetXaxis()->GetXmax());
  zeroLine->SetLineColor(kRed);
  zeroLine->SetLineStyle(2);
  zeroLine->Draw("same");

  // Redraw the histogram axes and title over the TPave
  hRes_histplot->Draw("AXIS SAME");

  // Update the canvas
  cTotal->cd();
  cTotal->Update();

  // Save the canvas to a file
  cTotal->SaveAs(savePath);
}

//Create canvas for reporting, expects a that the peaks are allowed to slide
void createAndSaveCanvasWithFitsAndResiduals_altbg(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, TH1D* hdxaltbg, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, std::string bgtype, double &elastics, std::pair<double,double> &rsf, int kine, int mag) {

  // Generate a unique identifier for the fit function name
  TString uniqueID = TString::Format("_%u", gRandom->Integer(1000000));

  //Clone the dx histos
  TH1D *hdx_histplot = (TH1D*)(hdx->Clone("hdx_histplot"));
  TH1D *hdx_clone = (TH1D*)(hdx->Clone("hdx_clone"));
  TH1D *hdx_p_clone = (TH1D*)(hdxp->Clone("hdx_p_clone"));
  TH1D *hdx_n_clone = (TH1D*)(hdxn->Clone("hdx_n_clone"));

  //Recreate the fit
  TF1* fit;
  std::string bgword = "";
  if(bgtype.compare("inel")==0){
    fit = new TF1(Form("fitinel_%s",uniqueID.Data()), fit_inelBG, hcalfit_l, hcalfit_h, 5);
    bgword = "MC Inelastic ";
  }else if(bgtype.compare("dy")==0){
    fit = new TF1(Form("fitdy_%s",uniqueID.Data()), fit_dyantiBG, hcalfit_l, hcalfit_h, 5);
    bgword = "dy Anticut ";
  }else if(bgtype.compare("coin")==0){
    fit = new TF1(Form("fitcoin_%s",uniqueID.Data()), fit_coinantiBG, hcalfit_l, hcalfit_h, 5);
    bgword = "Coin Time Anticut ";
  }else{
    cout << "ERROR: Invalid MC background type. Must be inel, antidy, or anticoin." << endl;
    return;
  }

  for (size_t i = 0; i < pars.size(); ++i) {
    fit->SetParameter(i, pars[i].first);
    fit->SetParError(i, pars[i].second);
  }

  std::pair<double,double> altchisqrndf = calculateChiSquareWithZeroMCError(hdx, fit);
  
  //fit->SetChisquare(altchisqrndf.first);
  //fit->SetNDF(altchisqrndf.second);

  fit->SetChisquare(qual.first);
  fit->SetNDF(qual.second);

  // Create a canvas with two pads
  TCanvas *cTotal = new TCanvas(savePath, Form("MC BG Data with Fits and Residuals %s",type), 1200, 800);
  TPad *pad1 = new TPad("pad1", "Pad with the fit", 0.0, 0.3, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "Pad with the residuals", 0.0, 0.0, 1.0, 0.3);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad1->Draw();
  pad2->Draw();

  // First pad: Draw the histogram with the fits
  pad1->cd();

  gPad->SetGridx();
  gPad->SetGridy();

  //hdx_histplot->SetLineColor(kBlack);
  //hdx_histplot->SetLineWidth(1);
  hdx_histplot->SetTitle(Form("dx, %s (SBS-%d, %d%% field);m", type, kine, mag));
  hdx_histplot->SetTitleFont(132);
  hdx_histplot->GetXaxis()->SetTitleFont(132);
  hdx_histplot->GetXaxis()->SetLabelFont(132);
  hdx_histplot->GetYaxis()->SetTitleFont(132);
  hdx_histplot->GetYaxis()->SetLabelFont(132);
  hdx_histplot->SetMarkerStyle(21);
  hdx_histplot->SetMarkerSize(0.5);
  hdx_histplot->Draw("P");
  hdx->SetLineColor(kBlack);
  hdx->SetLineWidth(1);
  hdx->Draw("E same");
  fit->SetNpx(1000);
  fit->SetLineColor(kGreen);
  fit->SetLineWidth(1);
  fit->SetFillColorAlpha(kGreen,0.35);
  fit->SetFillStyle(1001);
  fit->Draw("sameFC");

  hdx_p_clone = util::shiftHistogramX(hdx_p, pars[2].first);
  hdx_p_clone->Scale(pars[0].first);
  hdx_p_clone->SetLineColor(kRed-5);
  hdx_p_clone->SetLineWidth(2);
  hdx_p_clone->SetMarkerStyle(20);  // Filled circle marker
  hdx_p_clone->SetMarkerColor(kRed);
  hdx_p_clone->SetMarkerSize(0.5);
  hdx_p_clone->Draw("P same");

  hdx_n_clone = util::shiftHistogramX(hdx_n, pars[3].first);
  hdx_n_clone->Scale(pars[1].first);
  hdx_n_clone->SetLineColor(kBlue-5);
  hdx_n_clone->SetLineWidth(2);
  hdx_n_clone->SetMarkerStyle(20);  // Filled circle marker
  hdx_n_clone->SetMarkerColor(kBlue);
  hdx_n_clone->SetMarkerSize(0.5);
  hdx_n_clone->Draw("P same");

  // Set background function
  TH1D *hdx_altbg_clone = util::shiftHistogramX(hdxaltbg,pars[3].first);
  hdx_altbg_clone->Scale(pars[4].first);

  hdx_altbg_clone->SetLineColor(kBlack);
  hdx_altbg_clone->SetFillColorAlpha(kBlack,0.35);
  hdx_altbg_clone->SetFillStyle(3013);
  hdx_altbg_clone->Draw("hist same");

  TH1D *hdx_altbg_clone_clone = util::shiftHistogramX(hdxaltbg,pars[3].first);
  hdx_altbg_clone_clone->Scale(pars[4].first);
  hdx_altbg_clone_clone->SetLineColor(kBlack);
  hdx_altbg_clone_clone->SetLineWidth(1);
  hdx_altbg_clone_clone->Draw("E same");

  double psumerror;
  double psum = hdx_p_clone->IntegralAndError(0, hdx_p_clone->GetNbinsX() + 1, psumerror, "");

  double pscaleerror = pars[0].second;
  double pscale = pars[0].first;

  double nsumerror;
  double nsum = hdx_n_clone->IntegralAndError(0, hdx_n_clone->GetNbinsX() + 1, nsumerror, "");

  double nscaleerror = pars[1].second;
  double nscale = pars[1].first;

  double np_sum_ratio = nsum/psum;
  double np_par_ratio = nscale/pscale;
  double np_sum_ratio_error = np_sum_ratio * sqrt(pow(psumerror/psum, 2) + pow(nsumerror/nsum, 2));
  double np_par_ratio_error = np_par_ratio * sqrt(pow(pscaleerror/pscale, 2) + pow(nscaleerror/nscale, 2));

  rsf.first = np_par_ratio;
  rsf.second = np_par_ratio_error;

  // Construct legend within function
  TLegend* leg = new TLegend(0.5, 0.5, 0.89, 0.89);
  leg->SetBorderSize(0); // Remove border
  //leg->SetFillStyle(0);  // Make background transparent
  leg->SetTextColor(kBlack);
  leg->SetTextFont(132);
  leg->SetTextSize(0.04);
  leg->AddEntry( hdx_histplot, "Data", "p" );
  leg->AddEntry( fit, "Fit (Proton + Neutron + BG)", "f" );
  leg->AddEntry( hdx_altbg_clone, Form("%s Background",bgword.c_str()), "f" );
  leg->AddEntry( hdx_p_clone, "Proton SIMC MC", "p" );
  leg->AddEntry( hdx_n_clone, "Neutron SIMC MC", "p" );
  leg->AddEntry( (TObject*)0, "", "");
  leg->AddEntry( (TObject*)0, Form("n/p scale ratio R_{sf} : %0.3f #pm %0.3f",np_par_ratio, np_par_ratio_error), "");
  leg->AddEntry( (TObject*)0, "", "");
  //leg->AddEntry( (TObject*)0, Form("data N events : %0.0f",hdx->GetEntries()), "");
  //leg->AddEntry( (TObject*)0, Form("n/p yield ratio R'' : %0.3f #pm %0.3f",np_sum_ratio, np_sum_ratio_error), "");
  //leg->AddEntry( (TObject*)0, Form("dx shift pars, n/p : %0.3f / %0.3f ", pars[3].first, pars[2].first), "");
  //leg->AddEntry( (TObject*)0, Form("bg scale: %0.3f", pars[4].first), "");
  leg->AddEntry( (TObject*)0, Form("#chi^{2}/ndf: %0.3f/%d",fit->GetChisquare(),fit->GetNDF()), "");
  leg->Draw("same");

  //subtract the background from dx_clone and get elastic stats
  elastics = 0;
  for (int bin = 1; bin <= hdx_clone->GetNbinsX(); ++bin) {
    // Get the value at the current bin in bg_hist
    double bgValue = hdx_altbg_clone->GetBinContent(bin);
    if( hdx_clone->GetBinContent(bin) > bgValue )
      elastics += hdx_clone->GetBinContent(bin) - bgValue;
    // Subtract bgValue from the current bin content of hdx_clone and update it
    hdx_clone->SetBinContent(bin, hdx_clone->GetBinContent(bin) - bgValue);
    // Assuming you want to combine errors in quadrature
    double hdxError = hdx_clone->GetBinError(bin);
    double bgError = hdx_altbg_clone->GetBinError(bin);
    hdx_clone->SetBinError(bin, sqrt(hdxError*hdxError + bgError*bgError));
  }

  TH1D *hRes = calculateResiduals(hdx, fit, Form("hResiduals_%s",uniqueID.Data()));
  TH1D *hRes_histplot = (TH1D*)(hRes->Clone(Form("hRes_histplot_%s_%s",type,uniqueID.Data())));

  hRes_histplot->SetTitle("");
  hRes_histplot->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  double dxMaxValue = hdx_clone->GetMaximum();
  hRes_histplot->GetYaxis()->SetRangeUser(-dxMaxValue/4,dxMaxValue/4);

  hRes_histplot->SetLineColor(kBlack);
  hRes_histplot->SetLineWidth(1);
  //hRes_histplot->SetFillStyle(1001);
  //hRes_histplot->SetFillColorAlpha(kRed,0.75);

  // Second pad: Draw the residuals
  pad2->cd();
  gPad->SetGridx();

  hRes_histplot->GetYaxis()->SetTitle("Residuals");
  hRes_histplot->GetYaxis()->SetTitleSize(0.1);
  hRes_histplot->GetYaxis()->SetTitleOffset(0.5);
  hRes_histplot->GetYaxis()->SetLabelSize(0.1);
  hRes_histplot->GetXaxis()->SetTitleSize(0.1);
  hRes_histplot->GetXaxis()->SetLabelSize(0.1);
  hRes_histplot->GetXaxis()->SetTitle("x_{hcal}-x_{exp}");
  hRes_histplot->GetXaxis()->SetTitleFont(132);
  hRes_histplot->GetXaxis()->SetLabelFont(132);
  hRes_histplot->GetYaxis()->SetTitleFont(132);
  hRes_histplot->GetYaxis()->SetLabelFont(132);
  hRes_histplot->Draw("E");

  // Get the axis range
  double y_min = hRes_histplot->GetMinimum();
  double y_max = hRes_histplot->GetMaximum();

  // Create a TPave to fill the region between the axes
  TPave *pave = new TPave(hcalfit_l, y_min, hcalfit_h, y_max, 0, "NB");
  pave->SetFillColorAlpha(kGray,0.4); // Light grey RGB
  pave->SetFillStyle(3001); // 3001: transparent fill style
  pave->Draw("same");

  // Draw a line at y=0 for residuals
  TF1 *zeroLine = new TF1("zeroLine", "0", hRes->GetXaxis()->GetXmin(), hRes->GetXaxis()->GetXmax());
  zeroLine->SetLineColor(kRed);
  zeroLine->SetLineStyle(2);
  zeroLine->Draw("same");

  // Redraw the histogram axes and title over the TPave
  hRes_histplot->Draw("AXIS SAME");

  // Update the canvas
  cTotal->cd();
  cTotal->Update();

  // Save the canvas to a file
  cTotal->SaveAs(savePath);
}

//Function to fit then fine fit and return all fit parameters
std::vector<std::pair<double, double>> fitAndFineFit(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, const std::string& fitOptions = "RBMQ0") {
  TF1* fit = new TF1(fitName.c_str(), fitFormula.c_str(), hcalfit_l, hcalfit_h, paramCount);
  //fit->SetNpx(5000);
  for (int i=0; i<paramCount; ++i){ //reset parameters/errors for this set
    fit->SetParameter(i,0);
    fit->SetParError(i,0);
  }

  //fix bad peak fit for sbs4-50p with tight cuts. Informed from loose cut scales and residuals
  //fit->SetParLimits(0,0.0755,1.0);
  
  histogram->Fit(fit, fitOptions.c_str());

  std::vector<std::pair<double, double>> parametersAndErrors(paramCount);
  for (int i = 0; i < paramCount; ++i) {
    parametersAndErrors[i].first = fit->GetParameter(i); // Parameter value
    parametersAndErrors[i].second = fit->GetParError(i); // Parameter error
  }

  // Fine Fit
  TF1* fineFit = new TF1((fitName + "_fine").c_str(), fitFormula.c_str(), hcalfit_l, hcalfit_h, paramCount);
  std::vector<double> fineFitInitialParams(paramCount);
  for (int i = 0; i < paramCount; ++i) {
    fineFitInitialParams[i] = parametersAndErrors[i].first;
  }
  fineFit->SetParameters(fineFitInitialParams.data());

  //fix bad peak fit for sbs4-50p with tight cuts. Informed from loose cut scales and residuals
  //fineFit->SetParLimits(0,0.0755,1.0);
  
  histogram->Fit(fineFit, fitOptions.c_str());

  // Update parameters and errors with fine fit results
  for (int i = 0; i < paramCount; ++i) {
    parametersAndErrors[i].first = fineFit->GetParameter(i); // Fine fit parameter value
    parametersAndErrors[i].second = fineFit->GetParError(i); // Fine fit parameter error
  }

  cout << "P0!!! " << fineFit->GetParameter(0) << endl;
  
  fitqual.first = fineFit->GetChisquare();
  fitqual.second = fineFit->GetNDF();

  delete fit; // Delete fit to avoid carry-over
  delete fineFit; // Clean up
  return parametersAndErrors;
}

//Function to fit then fine fit and return all fit parameters
std::vector<std::pair<double, double>> fitAndFineFit_fixshift(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, double pshift, double nshift, const std::string& fitOptions = "RBMQ0") {
  TF1* fit = new TF1(fitName.c_str(), fitFormula.c_str(), hcalfit_l, hcalfit_h, paramCount);
  //fit->SetNpx(5000);
  for (int i=0; i<paramCount; ++i){ //reset parameters/errors for this set
    fit->SetParameter(i,0);
    fit->SetParError(i,0);
  }

  // double histfitmin = histogram->GetMaximum()/100;
  // double histfitmax = histogram->GetMaximum();
  // double histfitrmin = histogram->GetXaxis()->GetXmin();
  // double histfitrmax = histogram->GetXaxis()->GetXmax();
  // double histfitrange = histogram->GetXaxis()->GetXmax()-histogram->GetXaxis()->GetXmin();

  //fix bad peak fit for sbs4-50p with tight cuts. Informed from loose cut scales and residuals.
  //fit->SetParLimits(0,0.0755,1.0);
  
  double histfitmin = histogram->GetMaximum()/100;
  double histfitmax = histogram->GetMaximum();
  double histfitrmin = hcalfit_l;
  double histfitrmax = hcalfit_h;
  double histfitrange = histfitrmax-histfitrmin;

  //wooohoooo magic numbers!!!
  if(fitName.compare("gausBG")==0){
    fit->SetParLimits(4,histfitmin,histfitmax);
    fit->SetParLimits(5,histfitrmin+histfitrange/3,histfitrmax-histfitrange/3);
    fit->SetParLimits(6,histfitrange/10,histfitrange);

    cout << "GAUS MEAN HIST LIMS: " << histfitrmin << "  " << histfitrmax << " " << histfitrange/4 << endl;

  }

  fit->FixParameter(2,pshift); //fix x shift for proton
  fit->FixParameter(3,nshift); //fix x shift for neutron

  histogram->Fit(fit, fitOptions.c_str());

  std::vector<std::pair<double, double>> parametersAndErrors(paramCount);
  for (int i = 0; i < paramCount; ++i) {
    parametersAndErrors[i].first = fit->GetParameter(i); // Parameter value
    parametersAndErrors[i].second = fit->GetParError(i); // Parameter error
  }

  // Fine Fit
  TF1* fineFit = new TF1((fitName + "_fine").c_str(), fitFormula.c_str(), hcalfit_l, hcalfit_h, paramCount);
  std::vector<double> fineFitInitialParams(paramCount);
  for (int i = 0; i < paramCount; ++i) {
    fineFitInitialParams[i] = parametersAndErrors[i].first;
  }

  //wooohoooo magic numbers!!!
  if(fitName.compare("gausBG")==0){
    fineFit->SetParLimits(4,histfitmin,histfitmax);
    fineFit->SetParLimits(5,histfitrmin+histfitrange*0.4,histfitrmax-histfitrange*0.4);
    fineFit->SetParLimits(6,histfitrange/10,histfitrange);
  }

  fineFit->FixParameter(2,pshift); //keep fixed x shift
  fineFit->FixParameter(3,nshift); //keep fixed x shift

  fineFit->SetParameters(fineFitInitialParams.data());

  //fix bad peak fit for sbs4-50p with tight cuts. Informed from loose cut scales and residuals.
  //fineFit->SetParLimits(0,0.0755,1.0);
  
  histogram->Fit(fineFit, fitOptions.c_str());

  // Update parameters and errors with fine fit results
  for (int i = 0; i < paramCount; ++i) {
    parametersAndErrors[i].first = fineFit->GetParameter(i); // Fine fit parameter value
    parametersAndErrors[i].second = fineFit->GetParError(i); // Fine fit parameter error
  }

  fitqual.first = fineFit->GetChisquare();
  fitqual.second = fineFit->GetNDF();

  if(fitName.compare("gausBG")==0)
    for( auto par : parametersAndErrors )
      cout << par.first << endl;

  delete fit; // Delete fit to avoid carry-over
  delete fineFit; // Clean up
  return parametersAndErrors;
}

//simple background subtraction function
void subtractBackground(TH1D* hist, TF1* bgFit) {
  if (!hist || !bgFit) {
    std::cerr << "Error: Null pointer passed to subtractBackground." << std::endl;
    return;
  }

  for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
    double bgValue = bgFit->Eval(hist->GetXaxis()->GetBinCenter(bin));
    double content = hist->GetBinContent(bin) - bgValue;
    double error = hist->GetBinError(bin); // Assuming you want to keep the original error

    hist->SetBinContent(bin, content);
    hist->SetBinError(bin, error);

  }
}

int FindCutIndex(const std::string& branches, const std::string& searchStr) {
  std::istringstream iss(branches);
  std::string token;
  int count = 0;

  // Split the branches string by "&&" and search for the smaller string
  while (std::getline(iss, token, '&')) {
    // Check if the next character is also '&' to correctly skip "&&" delimiters
    if(iss.peek() == '&') iss.ignore();

    // Check if the token contains the searchStr
    if (token.find(searchStr) != std::string::npos) {
      return count; // Return the current count (index) if found
    }
    ++count; // Increment count for each passed "&&"
  }
  return -1; // Return -1 if the searchStr is not found
}

std::string RemovePrefix(const std::string& originalString, const std::string& prefix) {
    if (originalString.substr(0, prefix.length()) == prefix) {
        // If the string starts with the prefix, remove it
        return originalString.substr(prefix.length());
    }
    // If the string does not start with the prefix, return it unchanged
    return originalString;
}

TH1D* calculateResiduals(TH1D* hData, TF1* fit, const char* residualName) {
    // Clone the data histogram to create the residuals histogram
    TH1D* hResiduals = (TH1D*)(hData->Clone(residualName));
    hResiduals->Reset();  // Clear the histogram

    // Loop over each bin and calculate the residuals
    for (int bin = 1; bin <= hData->GetNbinsX(); ++bin) {
        double x = hData->GetXaxis()->GetBinCenter(bin);
        double dataValue = hData->GetBinContent(bin);
        double fitValue = fit->Eval(x);

        // Calculate the residual as (data - fit)
        double residual = dataValue - fitValue;
        double error = hData->GetBinError(bin); // Assuming Poisson statistics, you can also add fit and bg errors if needed

        hResiduals->SetBinContent(bin, residual);
        hResiduals->SetBinError(bin, error);
    }

    return hResiduals;
}

void chi2_constrained(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
    // Retrieve the histogram to fit
    TH1D *hist = (TH1D*)gMinuit->GetObjectFit();
    
    Double_t chi2 = 0;
    Int_t nbins = hist->GetNbinsX();
    
    // Loop over bins to calculate chi2
    for (int i = 1; i <= nbins; ++i) {
        Double_t x = hist->GetBinCenter(i);
        Double_t y = hist->GetBinContent(i);
        Double_t yerr = hist->GetBinError(i);
        Double_t yfit = fit_pol3BG(&x, par);
        
        if (yerr > 0) {
            chi2 += TMath::Power((y - yfit) / yerr, 2);
        }
    }
    
    // Constraints: Fit should pass through (x1, y1) and (x2, y2)
    Double_t x1 = hist->GetXaxis()->GetXmin(); // Example endpoint x1
    Double_t y1 = hist->GetBinContent(hist->FindBin(x1));  // Example endpoint y1
    Double_t x2 = hist->GetXaxis()->GetXmax(); // Example endpoint x2
    Double_t y2 = hist->GetBinContent(hist->FindBin(x2));  // Example endpoint y2
    
    Double_t yfit1 = fit_pol3BG(&x1, par);
    Double_t yfit2 = fit_pol3BG(&x2, par);
    
    chi2 += 1000 * TMath::Power(y1 - yfit1, 2); // Large weight to enforce the constraint
    chi2 += 1000 * TMath::Power(y2 - yfit2, 2); // Large weight to enforce the constraint
    
    f = chi2;
}

std::pair<double, double> calculateChiSquareWithZeroMCError(TH1D* histogram, TF1* fitFunc) {
  int nbins = histogram->GetNbinsX();
  double chi2 = 0;
  int ndf = 0;

  for (int i = 1; i <= nbins; ++i) {
    double x = histogram->GetBinCenter(i);
    double y_data = histogram->GetBinContent(i);
    double yerr_data = histogram->GetBinError(i);
    //double yerr_data = 1; //check to see if the by-eye agreement is better when not accounting for data err
    double y_fit = fitFunc->Eval(x);

    if (yerr_data > 0) {
      chi2 += TMath::Power((y_data - y_fit) / yerr_data, 2);
      ndf++;
    }
  }

  ndf -= fitFunc->GetNpar(); // Subtract the number of fit parameters from the number of degrees of freedom

  std::cout << "Chi-square (zero MC error): " << chi2 << std::endl;
  std::cout << "Degrees of freedom: " << ndf << std::endl;
  std::cout << "Chi-square/ndf: " << chi2 / ndf << std::endl;

  return std::make_pair(chi2, ndf);
}

void plotBackgrounds(const std::vector<TF1*>& bgFunctions) {
  // Create a canvas
  TCanvas *c1 = new TCanvas("c1", "Background Fits", 800, 600);

  // Define colors and line styles
  std::vector<int> colors = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange, kBlack, kViolet, kSpring};
  std::vector<int> styles = {1, 2, 3, 4, 5, 6, 7, 8, 9};

  // Ensure we have enough colors and styles
  int numFunctions = bgFunctions.size();
  if (numFunctions > colors.size()) {
    for (int i = colors.size(); i < numFunctions; ++i) {
      colors.push_back(kBlack); // Add default color if not enough
    }
  }
  if (numFunctions > styles.size()) {
    for (int i = styles.size(); i < numFunctions; ++i) {
      styles.push_back(1); // Add default style if not enough
    }
  }

  // Draw each background function with different colors and styles
  for (int i = 0; i < numFunctions; ++i) {
    bgFunctions[i]->SetLineColor(colors[i]);
    bgFunctions[i]->SetLineStyle(styles[i]);
    if (i == 0) {
      bgFunctions[i]->Draw();
    } else {
      bgFunctions[i]->Draw("same");
    }
  }

  // Create a legend
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  for (int i = 0; i < numFunctions; ++i) {
    legend->AddEntry(bgFunctions[i], bgFunctions[i]->GetName(), "l");
  }
  legend->Draw();

  // Update the canvas
  c1->Update();
}
