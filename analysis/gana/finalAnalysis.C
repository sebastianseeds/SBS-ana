//sseeds 1.22.23

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
#include "../../include/gmn.h"
#include "../../src/jsonmgr.C"

//Fit range
double hcalfit_l; //default defined by json
double hcalfit_h; //default defined by json

//Fit range override options
double hcalfit_l_or = -1.8; //lower fit/bin limit for hcal dx plots (m) sbs4, 50p
double hcalfit_h_or = 0.7; //upper fit/bin limit for hcal dx plots (m) sbs4, 50p

//Plot range option shared by all fits
double hcalr_l = -1.8; //lower fit/bin limit for hcal dx plots (m) sbs4, 50p
double hcalr_h = 0.7; //upper fit/bin limit for hcal dx plots (m) sbs4, 50p

//Fit options
std::string fitopt = "RMQ0";
//std::string fitopt = "RBMQ0"; //R:setrange,B:predefinedTF1,M:improvefit,Q:quiet,0:nofitline
//std::string fitopt = "RLQ0"; //L:loglikelihood for counts
//std::string fitopt = "RLEMQ0"; //E:bettererrorest
//std::string fitopt = "RWLEMQ0"; //WL:weightedloglikelihood for weighted counts (N/A here)
//std::string fitopt = "RPEMQ0"; //P:pearsonloglikelihood for expected error (N/A here)

//fit min entries
int minEvents = 500;

//get random seed
string passkey = "fltx";

//set up exclusions
const vector<std::string> exclusions = {"hist_bb_tr_n",
					"hist_bb_sh_nclus",
					"hist_bb_sh_colblk",
					"hist_bb_sh_rowblk",
					"hist_sbs_hcal_nclus",
					"hist_hcalx",
					"hist_hcaly",
					"hist_thetapq_n",
					"hist_thetapq_p",
					"hist_hcalnblk",
					"hist_fiducial_sig_y",
					"hist_fiducial_sig_x",
					"hist_bb_gem_track_chi2ndf"};

//Total fits using Interpolate with elastic signal histo and 4th order poly fit to bg
TH1D *hdx_p;
TH1D *hdx_n;
TH1D *hdx_inel;
TH1D *hdx_dyanti;
TH1D *hdx_coinanti;

Double_t fitFull(double *x, double *par){
  double dx = x[0];
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double proton = proton_scale * hdx_p->Interpolate(dx);
  double neutron = neutron_scale * hdx_n->Interpolate(dx);
  return proton + neutron + fits::g_p4fit(x,&par[2]);
}

Double_t fitFullInel(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double bg_scale = par[2];

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0]);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0]);
  double bg = bg_scale * hdx_inel->Interpolate(x[0]);
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron + bg;
}

Double_t fitFull_nobg(double *x, double *par){
  double dx = x[0];
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double proton = proton_scale * hdx_p->Interpolate(dx);
  double neutron = neutron_scale * hdx_n->Interpolate(dx);
  return proton + neutron;
}

Double_t fitFullShift(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron + fits::g_p4fit(x, &par[4]);
}

Double_t fitFullShift_p2(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron + fits::g_p2fit(x, &par[4]);
}

Double_t fitFullShiftInel(double *x, double *par){
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

Double_t fitFullShiftDyAnti(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  
  double bg_scale = par[4];

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);
  double bg = bg_scale * hdx_dyanti->Interpolate(x[0]); //no shift for data BG
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron + bg;
}

Double_t fitFullShiftCoinAnti(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  
  double bg_scale = par[4];

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);
  double bg = bg_scale * hdx_coinanti->Interpolate(x[0]); //no shift for data BG
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron + bg;
}

Double_t fitFullShift_nobg(double *x, double *par){
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

TH1D *hdx_slice_p;
TH1D *hdx_slice_n;

Double_t fitSlices(double *x, double *par){
  double dx = x[0];
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron + fits::g_p4fit(x, &par[4]);
}

Double_t fitSlices_p2(double *x, double *par){
  double dx = x[0];
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_slice_p->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_slice_n->Interpolate(x[0] - dx_shift_n);
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron + fits::g_p2fit(x, &par[4]);
}

Double_t fitSlices_nobg(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron;
}

TH1D *hdx_p_full;
TH1D *hdx_n_full;

Double_t fitSliceShift(double *x, double *par){

  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p_full->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n_full->Interpolate(x[0] - dx_shift_n);
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron + fits::g_p4fit(x, &par[4]);
}

Double_t fitSliceShift_nobg(double *x, double *par){

  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p_full->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n_full->Interpolate(x[0] - dx_shift_n);
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron;
}

Double_t SBpol4rej_b; //Central-band fit begin
Double_t SBpol4rej_e; //Central-band fit end

Double_t BGfit(double *x, double *par){

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
void createAndSaveCanvasWithFitsAndResiduals(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, TF1* bg, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, bool shiftfit, double &elastics);
void createAndSaveCanvasWithFitsAndResiduals_nobg(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, bool shiftfit, double &elastics);
void createAndSaveCanvasWithFitsAndResiduals_altbg(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, TH1D* hdxinel, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, std::string bgtype, double &elastics);
//std::vector<double> fitAndFineFit(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, const std::string& fitOptions);
std::vector<std::pair<double, double>> fitAndFineFit(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, const std::string& fitOptions = "RBMQ0");
std::vector<std::pair<double, double>> fitAndFineFit_slice(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, double pshift, double nshift, const std::string& fitOptions = "RBMQ0");
void subtractBackground(TH1D* hist, TF1* bgFit);
int FindCutIndex(const std::string& branches, const std::string& searchStr);
std::string RemovePrefix(const std::string& originalString, const std::string& prefix);

//fitranges
//sbs9 70p: -1.8 to 0.7
//sbs8 70p: -1.8 to 0.7
//sbs4 30p: -1.7 to 0.7, shiftX=0.0, neutronshift=-0.05
//sbs4 50p: -2.1 to 0.7, shiftX=0.05, neutronshift=0.0
//sbs7 85p: -1.4 to 0.5

//main. kine=kinematic, mag=fieldsetting, pass=pass#, sb_min/max=sidebandlimits, shiftX=shifttodxdata, N=cutvarsliceN
void finalAnalysis(int kine=9, 
		   int mag=70, 
		   int pass=2, 
		   double sb_min=-1.6, 
		   double sb_max=0.5, 
		   double shiftX=0.0, 
		   double neutronshift=0.0, 
		   int N=12, 
		   bool detail=true, 
		   bool mc_slices=true, 
		   bool blind=false, 
		   bool bestclus=true, 
		   bool thin=true, 
		   bool fitrangeoverride=true, 
		   bool alt = false) {
  //set draw params
  gStyle->SetPalette(55);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);

  //set sideband
  SBpol4rej_b = sb_min; 
  SBpol4rej_e = sb_max;

  // Set up configuration and tune objects to load analysis parameters
  SBSconfig config(kine,0);

  // Obtain configuration pars from config class
  double E_beam = config.GetEbeam();
  double BB_angle = config.GetBBtheta_rad(); // In radians

  //load json file
  JSONManager *jmgr = new JSONManager("../../config/syst.json");

  //Get dx plot details
  hcalfit_l = jmgr->GetValueFromSubKey<int>( "hcalfit_l", Form("sbs%d",kine) ); //gets to 7cm HCal pos res
  hcalfit_h = jmgr->GetValueFromSubKey<int>( "hcalfit_h", Form("sbs%d",kine) ); //gets to 7cm HCal pos res

  if(fitrangeoverride){
    hcalfit_l = hcalfit_l_or;
    hcalfit_h = hcalfit_h_or;
  }

  //Get tight elastic cuts
  std::string globalcuts = jmgr->GetValueFromSubKey_str( Form("post_cuts_p%d",pass), Form("sbs%d_%d",kine,mag) );

  cout << "Loaded tight cuts: " << globalcuts << endl;

  std::vector<std::string> cuts = util::parseCuts(globalcuts); //this makes a vector of all individual cuts in the single globalcut string.

  std::cout << "Parsed cuts: " << std::endl;

  //Get tight elastic cut strings and limits
  std::string branches = jmgr->GetValueFromKey_str( Form("post_branches_p%d",pass) );
  int hbins = jmgr->GetValueFromSubKey<int>( "hbins", Form("sbs%d",kine) ); //gets to 7cm HCal pos res

  vector<double> cut_lims;
  vector<double> llims;
  vector<double> ulims;

  jmgr->GetVectorFromSubKey<double>(Form("cut_limits_p%d",pass),Form("sbs%d_%d",kine,mag),cut_lims);  
  
  for( size_t i=0; i<cut_lims.size(); ++i ){
    if( i%2==0 )
      llims.push_back(cut_lims[i]);
    else
      ulims.push_back(cut_lims[i]);
  }

  //get paths and set random seed
  Double_t blind_factor = 1.;
  if(blind)
    blind_factor = util::generateRandomNumber(passkey);

  std::string detailword = "";
  if(detail)
    detailword = "_detail";

  std::string bestclus_word = "";
  if(bestclus)
    bestclus_word = "_bestclus";

  std::string thin_word = "";
  if(thin)
    thin_word = "_thin";

  std::string alt_word = "";
  if(alt)
    alt_word = "_alt";

  std::string basePath = "/lustre19/expphy/volatile/halla/sbs/seeds";
  //std::string finPath = Form("%s/gmn_analysis/dx_correlations_sbs%d_mag%d_pass%d_rd2.root", basePath.c_str(), kine, mag, pass);
  std::string finPath = Form("%s/gmn_analysis/dx_correlations_sbs%d_mag%d_pass%d%s%s.root", basePath.c_str(), kine, mag, pass, bestclus_word.c_str(), thin_word.c_str());
  //std::string fmcinPath = Form("%s/gmn_analysis/dx_mc_sbs%d_mag%d_pass%d_rd2.root", basePath.c_str(), kine, mag, pass);
  std::string fmcinPath = Form("%s/gmn_analysis/dx_mc_sbs%d_mag%d_pass%d%s%s.root", basePath.c_str(), kine, mag, pass, thin_word.c_str(), alt_word.c_str());
  std::string foutPath = Form("%s/gmn_analysis/analyze_gmn_sbs%d_mag%d_pass%d%s%s.root", basePath.c_str(), kine, mag, pass, detailword.c_str(), bestclus_word.c_str());

  TFile* inputFile = new TFile(finPath.c_str());
  if (!inputFile || inputFile->IsZombie()) 
    handleError(inputFile,"inputFile");

  //Get histograms from 
  TH1D *hdx_data = dynamic_cast<TH1D*>(inputFile->Get("hdx_allcut"));
  hdx_data->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  hdx_dyanti = dynamic_cast<TH1D*>(inputFile->Get("hdx_dyanti"));
  hdx_dyanti->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  hdx_coinanti = dynamic_cast<TH1D*>(inputFile->Get("hdx_coinanti"));
  hdx_coinanti->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  if (!hdx_data || !hdx_dyanti || !hdx_coinanti) 
    handleError(inputFile,"hdx_data");

  TFile* inputFileMC = new TFile(fmcinPath.c_str(), "READ");
  if (!inputFileMC || inputFileMC->IsZombie())
    handleError(inputFile,inputFileMC,"inputFileMC");

  //fix the interpolate functions
  hdx_p = dynamic_cast<TH1D*>(inputFileMC->Get("hdx_p"));
  hdx_p->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  hdx_n = dynamic_cast<TH1D*>(inputFileMC->Get("hdx_n"));
  hdx_n->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  hdx_inel = dynamic_cast<TH1D*>(inputFileMC->Get("hdx_inel"));
  hdx_inel->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  hdx_inel->Scale(10e33); //account for lack of overall normalization from g4sbs generator
  if (!hdx_p || !hdx_n || !hdx_inel) 
    handleError(inputFile,inputFileMC,"hdx");

  TFile* outputFile = new TFile(foutPath.c_str(), "RECREATE");

  //Set up some histogram clones for analysis
  TH1D *hdx_shifted = util::shiftHistogramX(hdx_data, shiftX);
  TH1D *hdx_shifted_clone = (TH1D*)(hdx_shifted->Clone("hdx_shifted_clone"));
  TH1D *hdx_shifted_nobg = (TH1D*)(hdx_shifted->Clone("hdx_shifted_nobg"));
  TH1D *hdx_raw = (TH1D*)(hdx_data->Clone("hdx_raw"));
  TH1D *hdx_raw_clone = (TH1D*)(hdx_data->Clone("hdx_raw_clone"));
  TH1D *hdx_raw_nobg = (TH1D*)(hdx_data->Clone("hdx_raw_nobg"));
  TH1D *hdx_raw_inel = (TH1D*)(hdx_data->Clone("hdx_raw_inel"));
  TH1D *hdx_raw_inel_clone = (TH1D*)(hdx_data->Clone("hdx_raw_inel_clone"));
  TH1D *hdx_raw_antidy = (TH1D*)(hdx_data->Clone("hdx_raw_antidy"));
  TH1D *hdx_raw_antidy_clone = (TH1D*)(hdx_data->Clone("hdx_raw_antidy_clone"));
  TH1D *hdx_raw_anticoin = (TH1D*)(hdx_data->Clone("hdx_raw_anticoin"));
  TH1D *hdx_raw_anticoin_clone = (TH1D*)(hdx_data->Clone("hdx_raw_anticoin_clone"));
  TH1D *hdx_raw_p2 = (TH1D*)(hdx_data->Clone("hdx_raw_p2"));
  TH1D *hdx_raw_p2_nobg = (TH1D*)(hdx_data->Clone("hdx_raw_p2_nobg"));
  TH1D *hdx_raw_p2_clone = (TH1D*)(hdx_data->Clone("hdx_raw_p2_clone"));

  TH1D *hdx_sb_nobg = (TH1D*)(hdx_data->Clone("hdx_sb_nobg"));
  TH1D *hdx_shifted_sb = (TH1D*)(hdx_shifted->Clone("hdx_shifted_sb"));
  TH1D *hdx_shifted_inel = (TH1D*)(hdx_shifted->Clone("hdx_shifted_inel"));
  TH1D *hdx_shifted_inel_clone = (TH1D*)(hdx_shifted->Clone("hdx_shifted_inel_clone"));

  //clone the mc signal 
  hdx_p_full = (TH1D*)(hdx_p->Clone("hdx_p_full"));
  hdx_n_full = (TH1D*)(hdx_n->Clone("hdx_n_full"));
  TH1D *hdx_p_clone = (TH1D*)(hdx_p->Clone("hdx_p_clone"));
  hdx_n = util::shiftHistogramX(hdx_n, neutronshift); 
  TH1D *hdx_n_clone = (TH1D*)(hdx_n->Clone("hdx_n_clone"));

  //clone the alt background
  TH1D *hdx_inel_clone = (TH1D*)(hdx_inel->Clone("hdx_inel_clone"));
  TH1D *hdx_antidy_clone = (TH1D*)(hdx_dyanti->Clone("hdx_antidy_clone"));
  TH1D *hdx_anticoin_clone = (TH1D*)(hdx_coinanti->Clone("hdx_anticoin_clone"));

  //test canvas
  TH1D *hdx_data_test = (TH1D*)(hdx_data->Clone("hdx_data_test"));
  TH1D *hdx_p_test = (TH1D*)(hdx_p->Clone("hdx_p_test"));
  TH1D *hdx_n_test = (TH1D*)(hdx_n->Clone("hdx_n_test"));
  TH1D *hdx_inel_test = (TH1D*)(hdx_inel->Clone("hdx_inel_test"));
  TH1D *hdx_antidy_test = (TH1D*)(hdx_dyanti->Clone("hdx_antidy_test"));
  TH1D *hdx_anticoin_test = (TH1D*)(hdx_coinanti->Clone("hdx_anticoin_test"));

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
  //hdx_antidy_test->DrawNormalized("HIST SAME");
  //hdx_anticoin_test->DrawNormalized("HIST SAME");
  hdx_data_test->Draw("HIST SAME");
    
  // Add a legend
  TLegend *legend = new TLegend(0.7, 0.7, 0.9, 0.9);
  legend->AddEntry(hdx_p_test, "simc proton", "l");
  legend->AddEntry(hdx_n_test, "simc neutron", "l");
  legend->AddEntry(hdx_inel_test, "bg distribution", "l");
  //legend->AddEntry(hdx_antidy_test, "hdx_antidy_test", "l");
  //legend->AddEntry(hdx_anticoin_test, "hdx_anticoin_test", "l");
  legend->AddEntry(hdx_data_test, "data", "l");
  legend->Draw();
    
  // Update the canvas
  canvas->Update();
    
  // Optionally, save the canvas to a file
  canvas->SaveAs("normalized_histograms_test.png");

  //Run through various fits and background methods
  std::pair<double,double> fullQual;
  auto fullpar_vector = fitAndFineFit(hdx_shifted, "fullFit", "fitFull", 7, hcalfit_l, hcalfit_h, fullQual, fitopt.c_str());

  std::pair<double,double> shiftQual;
  auto shiftpar_vector = fitAndFineFit(hdx_raw, "shiftFit", "fitFullShift", 9, hcalfit_l, hcalfit_h, shiftQual, fitopt.c_str());

  ////for now, best fit to data
  std::pair<double,double> shiftQualp2;
  auto shiftp2par_vector = fitAndFineFit(hdx_raw_p2, "shiftFit", "fitFullShift_p2", 7, hcalfit_l, hcalfit_h, shiftQualp2, fitopt.c_str());

  std::pair<double,double> sbQual;
  auto sbpar_vector = fitAndFineFit(hdx_shifted_sb, "sbFit", "BGfit", 5, hcalfit_l, hcalfit_h, sbQual, fitopt.c_str());

  std::pair<double,double> shiftInelQual;
  auto shiftinel_vector = fitAndFineFit(hdx_raw_inel, "shiftFitInel", "fitFullShiftInel", 5, hcalfit_l, hcalfit_h, shiftInelQual, fitopt.c_str());

  std::pair<double,double> shiftAntidyQual;
  auto shiftantidy_vector = fitAndFineFit(hdx_raw_antidy, "shiftFitAntidy", "fitFullShiftDyAnti", 5, hcalfit_l, hcalfit_h, shiftAntidyQual, fitopt.c_str());

  std::pair<double,double> shiftAnticoinQual;
  auto shiftanticoin_vector = fitAndFineFit(hdx_raw_anticoin, "shiftFitAnticoin", "fitFullShiftCoinAnti", 5, hcalfit_l, hcalfit_h, shiftAnticoinQual, fitopt.c_str());


  // Make background functions
  TF1 *bg_fullFit = new TF1("bg_fullFit",fits::g_p4fit,hcalfit_l,hcalfit_h,5);
  for (int i=0; i<5; ++i){
    bg_fullFit->SetParameter(i,fullpar_vector[i+2].first);
    bg_fullFit->SetParError(i,fullpar_vector[i+2].second);
  }

  TF1 *bg_shiftFit = new TF1("bg_shiftFit",fits::g_p4fit,hcalfit_l,hcalfit_h,5);
  for (int i=0; i<5; ++i){
    bg_shiftFit->SetParameter(i,shiftpar_vector[i+4].first);
    bg_shiftFit->SetParError(i,shiftpar_vector[i+4].second);
  }

  TF1 *bg_shiftp2Fit = new TF1("bg_shiftp2Fit",fits::g_p2fit,hcalfit_l,hcalfit_h,3);
  for (int i=0; i<3; ++i){
    bg_shiftp2Fit->SetParameter(i,shiftp2par_vector[i+4].first);
    bg_shiftp2Fit->SetParError(i,shiftp2par_vector[i+4].second);
  }

  TF1 *bg_sbFit = new TF1("bg_sbFit",fits::g_p4fit,hcalfit_l,hcalfit_h,5);
  for (int i=0; i<5; ++i){
    bg_sbFit->SetParameter(i,sbpar_vector[i].first);
    bg_sbFit->SetParError(i,sbpar_vector[i].second);
  }

  //Print the canvas
  //pol4 fit to bg fixed peak locs
  //createAndSaveCanvasWithFitsAndResiduals(hdx_shifted_clone, hdx_p_clone, hdx_n_clone, bg_fullFit, fullpar_vector, fullQual, "pol4 BG", "~/gmn_plots/combined_dx_fits_residuals_shift.pdf", false);
  
  //pol4 fit to bg
  double elastics_p4;
  createAndSaveCanvasWithFitsAndResiduals(hdx_raw_clone, hdx_p_clone, hdx_n_clone, bg_shiftFit, shiftpar_vector, shiftQual, "shiftfit pol4 BG", "~/gmn_plots/combined_dx_fits_residuals_fitshift.pdf", true, elastics_p4);

  //pol2 fit to bg
  double elastics_p2;
  createAndSaveCanvasWithFitsAndResiduals(hdx_raw_p2_clone, hdx_p_clone, hdx_n_clone, bg_shiftp2Fit, shiftp2par_vector, shiftQual, "shiftfit pol2 BG", "~/gmn_plots/combined_dx_fits_residuals_fitshiftp2.pdf", true, elastics_p2);

  //inel bg
  double elastics_inel;
  createAndSaveCanvasWithFitsAndResiduals_altbg(hdx_raw_inel_clone, hdx_p_clone, hdx_n_clone, hdx_inel_clone, shiftinel_vector, shiftInelQual, "shiftfit MC BG", "~/gmn_plots/combined_dx_fits_residuals_inel_fitshift.pdf", "inel", elastics_inel);
  
  //dy anticut bg
  double elastics_dy;
  createAndSaveCanvasWithFitsAndResiduals_altbg(hdx_raw_antidy_clone, hdx_p_clone, hdx_n_clone, hdx_antidy_clone, shiftantidy_vector, shiftAntidyQual, "shiftfit antidy BG", "~/gmn_plots/combined_dx_fits_residuals_antidy_fitshift.pdf", "dy", elastics_dy);
  
  //coin anticut bg
  double elastics_coin;
  createAndSaveCanvasWithFitsAndResiduals_altbg(hdx_raw_anticoin_clone, hdx_p_clone, hdx_n_clone, hdx_anticoin_clone, shiftanticoin_vector, shiftAnticoinQual, "shiftfit anticoin BG", "~/gmn_plots/combined_dx_fits_residuals_anticoin_fitshift.pdf", "coin", elastics_coin);

  //subtract bg fits
  subtractBackground(hdx_shifted_nobg,bg_fullFit);
  std::pair<double,double> fullnobgQual;
  auto fullpar_nobg_vector = fitAndFineFit(hdx_shifted_nobg, "fullFit_nobg", "fitFull_nobg", 2, hcalfit_l, hcalfit_h, fullnobgQual, fitopt.c_str());
  
  subtractBackground(hdx_raw_nobg,bg_shiftFit);
  std::pair<double,double> shiftnobgQual;
  auto shiftpar_nobg_vector = fitAndFineFit(hdx_raw_nobg, "shiftFit_nobg", "fitFullShift_nobg", 4, hcalfit_l, hcalfit_h, shiftnobgQual, fitopt.c_str());

  subtractBackground(hdx_raw_p2_nobg,bg_shiftp2Fit);
  std::pair<double,double> shiftp2nobgQual;
  auto shiftp2par_nobg_vector = fitAndFineFit(hdx_raw_p2_nobg, "shiftp2Fit_nobg", "fitFullShift_nobg", 4, hcalfit_l, hcalfit_h, shiftp2nobgQual, fitopt.c_str());

  subtractBackground(hdx_sb_nobg,bg_sbFit);
  std::pair<double,double> sbnobgQual;
  auto sbpar_nobg_vector = fitAndFineFit(hdx_sb_nobg, "sbFit_nobg", "fitFullShift_nobg", 4, hcalfit_l, hcalfit_h, sbnobgQual, fitopt.c_str());

  //Print the canvas
  //createAndSaveCanvasWithFitsAndResiduals_nobg(hdx_shifted_nobg, hdx_p_clone, hdx_n_clone, fullpar_nobg_vector, fullnobgQual, "pol4 BG subtracted", "~/gmn_plots/combined_dx_fits_residuals_full_nobg.pdf", false);
  
  double elastics_p4_nobg;
  createAndSaveCanvasWithFitsAndResiduals_nobg(hdx_raw_nobg, hdx_p_clone, hdx_n_clone, shiftpar_nobg_vector, shiftnobgQual, "shiftfit pol4 BG subtracted", "~/gmn_plots/combined_dx_fits_residuals_fitshift_nobg.pdf", true, elastics_p4_nobg);
  
  double elastics_p2_nobg;  
  createAndSaveCanvasWithFitsAndResiduals_nobg(hdx_raw_p2_nobg, hdx_p_clone, hdx_n_clone, shiftp2par_nobg_vector, shiftp2nobgQual, "shiftfit pol2 BG subtracted", "~/gmn_plots/combined_dx_fits_residuals_fitshiftp2_nobg.pdf", true, elastics_p2_nobg);
  
  double elastics_sb_nobg;
  createAndSaveCanvasWithFitsAndResiduals_nobg(hdx_sb_nobg, hdx_p_clone, hdx_n_clone, sbpar_nobg_vector, sbnobgQual, "shiftfit sideband BG subtracted", "~/gmn_plots/combined_dx_fits_residuals_fitsb_nobg.pdf", true, elastics_sb_nobg);


  /////////////
  /////////////
  /////////////
  /////////////
  /////////////
  /////////////
  /////////////
  //....proceed to detailed analysis

  if(!detail)
    return;

  // Setup report struct
  struct ReportData {
    TH1D* sliceHistogram = nullptr;
    TH1D* slicePHistogram = nullptr;
    TH1D* sliceNHistogram = nullptr;
    double fitParams[7] = {}; // Array for fit parameters
    double fitErrors[7] = {}; // Array for fit errors
    double winLow = 0;
    double winHigh = 0;
    double nev = 0;
    double mcpnev = 0;
    double mcnnev = 0;
    double scaleratio = 0;
    double scaleratioerr = 0;
    double chisqrndf = 0;
    std::string type;

    // Constructor for easier initialization (optional)
    ReportData(TH1D* hist = nullptr, TH1D* histp = nullptr, TH1D* histn = nullptr, double wl = 0, double wh = 0, double nv = 0, double mvp = 0, double mvn = 0, double sr = 0, double se = 0, double cs = 0, std::string t = "")
      : sliceHistogram(hist), slicePHistogram(histp), sliceNHistogram(histn), winLow(wl), winHigh(wh), nev(nv), mcpnev(mvp), mcnnev(mvn), scaleratio(sr), scaleratioerr(se), chisqrndf(cs), type(t) {
      // Initialize fitParams and fitErrors with default values if needed
      for (int i = 0; i < 7; ++i) {
    	fitParams[i] = 0.0;
    	fitErrors[i] = 0.0;
      }
    }
  };

  std::vector<std::vector<ReportData>> reports;

  std::vector<double> relativeErrors;
  std::vector<string> relativeErrorSources;

  // Loop over all TH2D histograms
  int iteridx = 0;
  TIter next(inputFile->GetListOfKeys());
  TKey* key;

  while ((key = dynamic_cast<TKey*>(next()))) {

    if (strcmp(key->GetClassName(), "TH2D") != 0) continue; // Process only TH2D objects

    //Get 2D histogram from data file
    TH2D* hist = dynamic_cast<TH2D*>(key->ReadObj());
    if (!hist) continue;

    //Get corresponding histogram from MC file
    std::string histName = hist->GetName();
    std::string histName_mc_p = histName + "_p";
    std::string histName_mc_n = histName + "_n";

    cout << "Loaded histogram " << histName << " for analysis." << endl;

    //Get proton and neutron histograms for slice analysis
    TH2D* hist_mc_p = dynamic_cast<TH2D*>(inputFileMC->Get(histName_mc_p.c_str()));
    TH2D* hist_mc_n = dynamic_cast<TH2D*>(inputFileMC->Get(histName_mc_n.c_str()));
    if (!hist_mc_p || !hist_mc_n) continue;

    // Check if the histogram name starts with "coorhist" and skip if so
    if (histName.find("coorhist") == 0) continue; // Skip if name starts with coorhist

    if (histName.find("hist_inel") == 0) continue; // Skip if name contains inel

    bool skip = false;

    for( size_t name=0; name<exclusions.size(); ++name){

      if (histName.find(exclusions[name]) == 0){
	cout << "Found exclusion: " << exclusions[name] << endl;
	skip=true;
	break;
      }
    }

    //skip on exclusions
    if(skip)
      continue;

    ///////////////////////////////
    // Start variable width binning
    ///////////////////////////////

    // Get llim and ulim and set the range for analysis
    std::string current_branch = RemovePrefix(histName, "hist_");
    int branch_index = FindCutIndex(branches,current_branch);
    double current_llim = llims[branch_index];
    cout << "For branch " << current_branch << " index " << branch_index << ", loaded the lower limit from json: " << current_llim << endl;

    double current_ulim = ulims[branch_index];
    cout << "For branch " << current_branch << " index " << branch_index << ", loaded the upper limit from json: " << current_ulim << endl;

    double histStart = hist->GetXaxis()->GetBinLowEdge(1);
    double histEnd = hist->GetXaxis()->GetBinUpEdge(hist->GetNbinsX());

    if( current_llim < histStart )
      current_llim = histStart;
    if( current_ulim > histEnd )
      current_ulim = histEnd;

    cout << "Working on dx vs " << current_branch << " with lower branch limit " << current_llim << " and upper branch limit " << current_ulim << endl;

    hist->GetXaxis()->SetRangeUser(current_llim, current_ulim);
    hist->GetYaxis()->SetRangeUser(hcalfit_l, hcalfit_h);

    // Assuming current_llim and current_ulim are set to the desired x-axis range
    int xbinLow = hist->GetXaxis()->FindBin(current_llim);
    int xbinHigh = hist->GetXaxis()->FindBin(current_ulim);

    // Project the TH2D onto a TH1D for the specified x-axis range
    TH1D* wideProjY = hist->ProjectionY("_py", xbinLow, xbinHigh);

    // Use Integral() to get the total number of events within the specified range
    int totalEvents = wideProjY->Integral();

    //set M = N at this point to introduce errors where loops shouldn't proceed.
    int M = N;

    // Make sure we have more events than the number of slices
    if (totalEvents < M) {
      cout << "ERROR: Total events are less than the number of requested slices M. Cannot proceed." << endl;
      return;
    }

    // Calculate the approximate number of events per slice
    int eventsPerSlice = totalEvents / M;

    cout << "Events per slice: " << eventsPerSlice << endl;

    int eventsAccumulated = 0;
    int currentBin = 1;
    std::vector<int> sliceBoundaries; // Stores the upper boundary of each slice
    std::vector<double> exactSliceEvents; //Stores the total number of events in each slice

    // Before starting the loop, ensure you add the first boundary
    sliceBoundaries.push_back(xbinLow);

    cout << "Getting slice boundaries for " << histName << " with " << eventsPerSlice << " events per slice over " << M << " slices." << endl;

    cout << "   lower bound: " << current_llim << endl;

    double sum_events_on_slice = 0;

    // Start the loop from the first bin within your range of interest
    for (int i = xbinLow; i <= xbinHigh; ++i) {
      // Accumulate events in this bin
      int eventsInBin = hist->ProjectionY("_py", i, i)->GetEntries();
      eventsAccumulated += eventsInBin;

      // Define a new slice if conditions are met
      if (eventsAccumulated >= eventsPerSlice || i == xbinHigh) { // Adjusted condition to include xbinHigh
        double sliceXUBound = hist->GetXaxis()->GetBinLowEdge(i + 1); // For logging purposes
	sliceBoundaries.push_back(i);
	exactSliceEvents.push_back(eventsAccumulated);
	sum_events_on_slice += eventsAccumulated;

        cout << "   acc events: " << eventsAccumulated << endl;
        cout << "   next bound: " << sliceXUBound << endl;

        eventsAccumulated = 0; // Reset for the next slice
      }
    }

    if( sum_events_on_slice != totalEvents ){
      cout << endl << endl << "WARNING: Total number of events summed up over all sliced != to total events used to calculate intervals." << endl;
      cout << "  sum ev: " << sum_events_on_slice << endl;
      cout << "  total events on slice: " << totalEvents << endl << endl;
    }

    double xLow = current_llim;
    double xHigh = current_ulim;

    // Create a canvas for this TH2D histogram
    TCanvas* th2dCanvas = new TCanvas(Form("canvas_%s", hist->GetName()), hist->GetTitle(), 800, 600);
    th2dCanvas->cd();
    TH2D *histout = (TH2D*)(hist->Clone("hdx_shifted_clone"));
    std::string histTitle = hist->GetTitle();
    std::string histoutTitle = histTitle + Form(". Appx ev/slice: %d", eventsPerSlice);
    histout->SetTitle(histoutTitle.c_str());
    histout->GetXaxis()->SetRangeUser(xLow,xHigh);
    histout->GetYaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
    histout->Draw("COLZ");  // Draw the histogram

    double maxEventsPerSlice = *std::max_element(exactSliceEvents.begin(), exactSliceEvents.end());

    for (size_t i = 0; i < sliceBoundaries.size() - 1; ++i) {
      int binLow = sliceBoundaries[i];
      int binHigh = sliceBoundaries[i + 1];
      double x1 = hist->GetXaxis()->GetBinLowEdge(binLow);
      double x2 = hist->GetXaxis()->GetBinUpEdge(binHigh);
    
      // Assuming the y-axis of the histogram only spans the non-empty bins
      double y1 = hist->GetYaxis()->GetBinLowEdge(hist->GetYaxis()->GetFirst());
      double y2 = hist->GetYaxis()->GetBinUpEdge(hist->GetYaxis()->GetLast());

      double sliceEvents = exactSliceEvents[i];
      double alpha = ( sliceEvents / maxEventsPerSlice ) / 2;

      // Draw a black line at the end of each slice
      TLine *line = new TLine(x2, y1, x2, y2);
      line->SetLineColor(kBlack);
      line->Draw("same");

    }

    th2dCanvas->Update();
    th2dCanvas->Write(); 

    // Set up TH2D report hists
    std::string ratioHistName = std::string("ratio: ") + hist->GetName() + ";" + hist->GetName() + ";scale factor n:p ratio";
    std::string ratioHistName_w = std::string("R, weighted proj: ") + hist->GetName() + "; scale ratio";;

    TH1D* ratioHist = new TH1D(ratioHistName.c_str(), ratioHistName.c_str(), N, xLow, xHigh);
    TH1D* ratio_w = new TH1D(ratioHistName_w.c_str(), ratioHistName_w.c_str(), 10*N, 0.5, 1.5);


    // Resize the vectors to fit the number of slices
    std::vector<std::vector<double>> scalefactor(M, std::vector<double>(7, 0));
    std::vector<std::vector<double>> scaleerr(M, std::vector<double>(7, 0));
    std::vector<double> scaleratio(M, 0);
    std::vector<double> scaleratio_err(M, 0);
    std::vector<double> totalevents(M, 0);
    std::vector<double> bin_low(M, 0);
    std::vector<double> bin_high(M, 0);

    // Ensure there's a slot for the current histogram's reports
    if (iteridx >= reports.size()) {
      reports.resize(iteridx + 1);
    }

    // Now, for the current TH2D histogram indexed by iteridx,
    // resize the inner vector to hold all M slices.
    reports[iteridx].resize(M);

    // Begin loop over slices
    for (int i = 0; i < M; ++i) {
      int binLow = sliceBoundaries[i]; // The lower boundary of the current slice
      int binHigh = (i < M - 1) ? sliceBoundaries[i + 1] - 1 : xbinHigh; // The upper boundary

      bin_low[i] = hist->GetXaxis()->GetBinLowEdge(binLow);
      bin_high[i] = hist->GetXaxis()->GetBinUpEdge(binHigh);
      cout << "Slicing on TH2D X axis range: " << bin_low[i] << " to " << bin_high[i] << endl;

      //Do slices for fits to MC
      TH1D *cellslice = hist->ProjectionY(Form("cellslice_%d", i+1), binLow, binHigh);

      //Default set the slice MC comparison equal to the total hist MC projY
      hdx_slice_p = hdx_p;
      hdx_slice_n = hdx_n;

      //Get slices from corresponding MC hist
      if(mc_slices){
	hdx_slice_p = hist_mc_p->ProjectionY(Form("%s_slicemc_p_%d",histName.c_str(), i+1), binLow, binHigh);
	hdx_slice_n = hist_mc_n->ProjectionY(Form("%s_slicemc_n_%d",histName.c_str(), i+1), binLow, binHigh);
      }

      //Reporting histograms
      TH1D *reportslice = hist->ProjectionY(Form("reportslice_%d_%s", i+1, histName.c_str()), binLow, binHigh);
      TH1D *reportslice_p = (TH1D*)(hdx_slice_p->Clone(Form("reportslice_p_%d_%s", i+1, histName.c_str())));
      TH1D *reportslice_n = (TH1D*)(hdx_slice_n->Clone(Form("reportslice_n_%d_%s", i+1, histName.c_str())));

      //int slice_nEntries = cellslice->Integral(1,sliceNBins); // Total number of entries in the slice
      int slice_nEntries = reportslice->Integral(1,reportslice->GetNbinsX()); // Total number of entries in the slice
      int slicemcp_nEntries = reportslice_p->Integral(1,reportslice_p->GetNbinsX());
      int slicemcn_nEntries = reportslice_n->Integral(1,reportslice_n->GetNbinsX());
      double slice_error = 0;
      if (slice_nEntries > 0) {
	slice_error = 1/sqrt(slice_nEntries);
      }

      //Rsf extraction using second order poly fit
      std::pair<double,double> sliceQualp2;

      cout << "Max value for mc scaled proton, neutron: ";
      auto slicep2par_vector = fitAndFineFit_slice(cellslice, "sliceFit", "fitSlices_p2", 7, hcalfit_l, hcalfit_h, sliceQualp2, shiftp2par_vector[2].first, shiftp2par_vector[3].first, fitopt.c_str());

      // //Rsf extraction using coin anticut fit
      // std::pair<double,double> sliceQualanticoin;

      // cout << "Max value for mc scaled proton, neutron: ";
      // auto slicep2par_vector = fitAndFineFit_slice(cellslice, "sliceFit", "fitSlices_p2", 7, hcalfit_l, hcalfit_h, sliceQualp2, shiftp2par_vector[2].first, shiftp2par_vector[3].first, fitopt.c_str());

      //fix x shift errors from overall fit to data. fitAndFineFit_slice() fixes this already.
      slicep2par_vector[2].second = shiftp2par_vector[2].second;
      slicep2par_vector[3].second = shiftp2par_vector[3].second;

      reportslice_p = (TH1D*)(hdx_slice_p->Clone(Form("reportslice_p_%d_%s", i+1, histName.c_str())));
      reportslice_n = (TH1D*)(hdx_slice_n->Clone(Form("reportslice_n_%d_%s", i+1, histName.c_str())));

      double csndf = sliceQualp2.first/sliceQualp2.second;

      for(int k = 0; k < 7; ++k) {
	double param = slicep2par_vector[k].first;
	double error = slicep2par_vector[k].second;

	cout << "On slice M=" << i << " par " << k << "=" << param << " with error " << error << endl;  

	// Directly populate the scalefactor and scaleerr arrays
	scalefactor[i][k] = param;
	scaleerr[i][k] = error;
      }

      // Calculate the scale factor ratio and its error
      if (scalefactor[i][0] != 0) { // Avoid division by zero
	scaleratio[i] = scalefactor[i][1] / scalefactor[i][0];
	scaleratio_err[i] = sqrt(pow(scaleerr[i][1]/scalefactor[i][1], 2) + pow(scaleerr[i][0]/scalefactor[i][0], 2)) * scaleratio[i];
      } else {
	scaleratio[i] = 0;
	scaleratio_err[i] = 0;
      }

      reports[iteridx][i] = ReportData( reportslice, 
					reportslice_p,
					reportslice_n,
					bin_low[i],
					bin_high[i],
					slice_nEntries,
					slicemcp_nEntries,
					slicemcn_nEntries,
					scaleratio[i],
					scaleratio_err[i],
					csndf,
					histName );

      if( bin_low[i]==bin_high[i] || bin_high[i] > xbinHigh || bin_high[i] < bin_low[i] ){

	TH1D* nullHist = new TH1D("nullHist", "Null Histogram", 500, 0, 1);

	reports[iteridx][i] = ReportData( nullHist, 
					  nullHist,
					  nullHist,
					  0,
					  1,
					  0,
					  0,
					  0,
					  0,
					  0,
					  0,
					  histName );

      }

      for (int par=0; par<7; ++par){
	reports[iteridx][i].fitParams[par] = scalefactor[i][par];
	reports[iteridx][i].fitErrors[par] = scaleerr[i][par];
      }

      // Fill the TH1Ds with ratios and yields
      ratioHist->SetBinContent(i, scaleratio[i]);
      ratioHist->SetBinError(i, scaleratio_err[i]);

    } //endloop over N slices

    //Now get weighted histograms
    // Loop over all bins and print the value in each bin
    int ratioHist_nBins = ratioHist->GetXaxis()->GetNbins();
    for (int k = 1; k <= ratioHist_nBins; k++) { // bin indices start from 1
      double ratioHist_binValue = ratioHist->GetBinContent(k);
      double ratioHist_binErr = ratioHist->GetBinError(k);
    
      // Use the square of the error as weight
      double weight = 1.0 / (ratioHist_binErr * ratioHist_binErr);
  
      // Check for division by zero
      if(ratioHist_binErr > 0) {
	ratio_w->Fill(ratioHist_binValue, weight);
      }

    }

    // Before writing the histogram, check for invalid values
    for (int i = 1; i <= ratioHist->GetNbinsX(); ++i) {
      if (std::isinf(ratioHist->GetBinContent(i)) || std::isnan(ratioHist->GetBinContent(i))) {
	std::cerr << "Invalid value found in histogram " << ratioHist->GetName() << " bin " << i << std::endl;
	ratioHist->SetBinContent(i, 0); // Setting to 0 or some default value
      }
    }

    //Write the raw histograms to file
    ratioHist->Write();
    ratio_w->Write();

    //Write out canvas for ratio histograms
    TCanvas *Rcanvas = new TCanvas(Form("R_{sf} %s",histName.c_str()), Form("R_{sf} %s",histName.c_str()), 800, 600);
    Rcanvas->cd();

    // Set the marker style to stars
    ratio_w->SetMarkerStyle(29); // 29 is the marker code for stars in ROOT
    ratio_w->SetMarkerSize(1.4); // Adjust the marker size as needed
    
    // Draw the histogram with the option "P" (draws the histogram using the current marker)
    ratio_w->Draw("hist P");
    // Retrieve the Standard Deviation (SD) of the histogram and mean. Get contribution to syst error.
    double sd = ratio_w->GetStdDev();
    double mean = ratio_w->GetMean();
    
    //relativeErrors.push_back(pow(sd/mean,2));
    relativeErrors.push_back(pow(sd/mean,2));

    std::string prefixToRemove = "hist_";

    if (histName.substr(0, prefixToRemove.length()) == prefixToRemove) {
        histName = histName.substr(prefixToRemove.length());
    }

    relativeErrorSources.push_back(histName.c_str());

    // If you also need the Standard Deviation error, you can retrieve it as follows
    double sdError = ratio_w->GetStdDevError();

    // Create a legend to display the SD
    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9); // Adjust the coordinates to fit your canvas
    legend->SetHeader(Form("R_{sf} %s slices",histName.c_str()), "C"); // Optional header

    // Format the SD text and add it to the legend
    char sdText[255];
    sprintf(sdText, "SD = %.2f #pm %.2f", sd, sdError);
    char evText[255];
    sprintf(evText, "entries = %.2f", ratio_w->GetEntries());
    legend->AddEntry((TObject*)0, sdText, "");
    legend->AddEntry((TObject*)0, evText, "");

    // Draw the legend on the canvas
    legend->Draw();

    // Write the histogram to the current directory or file
    Rcanvas->Write();

    cout << endl << iteridx << endl;

    iteridx++;

  }//endloop over TH2D histograms

  //set up vectors for later tgrapherrors
  std::vector<std::vector<double>> Rsf_vec;
  std::vector<std::vector<double>> Rsferr_vec;
  std::vector<std::vector<double>> xval_vec;
  
  Rsf_vec.resize(iteridx);
  Rsferr_vec.resize(iteridx);
  xval_vec.resize(iteridx);

  //write out reports per cut
  for (int r = 0; r < iteridx; ++r) {

    TCanvas* canvasSlices = new TCanvas(Form("dx slices over %s", reports[r][0].type.c_str()), Form("dx slices over %s", reports[r][0].type.c_str()), 1800, 1200);

    // Assuming N is the total number of histograms to be plotted on the canvas
    int optimalRows = 1;
    int optimalCols = 1;

    // Find the nearest square root for N to determine the grid size
    int sqrtN = (int)std::sqrt(N);
    for (int cols = sqrtN; cols <= N; ++cols) {
      if (N % cols == 0) { // If cols is a factor of N
	optimalCols = cols;
	optimalRows = N / cols;
	break; // Found the optimal layout
      }
    }

    canvasSlices->Divide(optimalCols, optimalRows); // Adjust the division based on reportSamples or desired layout

    TF1* bgfit[N];

    // loop over slices
    for (int i = 0; i < N; ++i) {
      if (reports[r][i].sliceHistogram == nullptr)
        continue;

      canvasSlices->cd(i + 1);

      //catch bad fit ranges and do not write
      std::string tempTitle = reports[r][i].sliceHistogram->GetTitle();      
      if( tempTitle.compare("Null Histogram")==0 )
	continue;

      //get mc nucleon parameters
      double p_scale = reports[r][i].fitParams[0];
      double n_scale = reports[r][i].fitParams[1];
      double p_shift = reports[r][i].fitParams[2];
      double n_shift = reports[r][i].fitParams[3];

      double MCev = reports[r][i].mcpnev + reports[r][i].mcnnev;

      // Simplify the histogram title
      reports[r][i].sliceHistogram->SetTitle(Form("%s %0.2f to %0.2f", reports[r][0].type.c_str(), reports[r][i].winLow, reports[r][i].winHigh));
      reports[r][i].sliceHistogram->SetLineWidth(1);
      reports[r][i].sliceHistogram->Draw();

      // Retrieve the number of bins, and the x-axis limits of the slice histogram
      int nbins = reports[r][i].sliceHistogram->GetNbinsX();
      double x_low = reports[r][i].sliceHistogram->GetXaxis()->GetXmin();
      double x_high = reports[r][i].sliceHistogram->GetXaxis()->GetXmax();

      //get scale ratio, error, and midpoint for slice
      double ratio = reports[r][i].scaleratio;
      double error = reports[r][i].scaleratioerr;
      double xval_midpoint = (reports[r][i].winLow + reports[r][i].winHigh)/2;

      //fill vectors for later tgrapherrors
      Rsf_vec[r].push_back(ratio);
      Rsferr_vec[r].push_back(error);
      xval_vec[r].push_back(xval_midpoint);

      // Draw the corresponding MC histograms
      hdx_slice_p = reports[r][i].slicePHistogram; //update the MC slice for the later overall fit
      TH1D *pMCslice = util::shiftHistogramX(reports[r][i].slicePHistogram, p_shift);
      pMCslice->Scale(p_scale);

      hdx_slice_n = reports[r][i].sliceNHistogram; //update the MC slice for the later overall fit
      TH1D *nMCslice = util::shiftHistogramX(reports[r][i].sliceNHistogram, n_shift);
      nMCslice->Scale(n_scale);
      
      pMCslice->Draw("same");
      nMCslice->Draw("same");

      // Create a legend and add the remaining parameters
      TLegend* legend = new TLegend(0.59, 0.69, 0.89, 0.89); // Adjust the position as needed
      legend->AddEntry((TObject*)0, Form("Nev: %0.0f", reports[r][i].nev), "");
      legend->AddEntry((TObject*)0, Form("MCev: %0.0f", MCev), "");
      legend->AddEntry((TObject*)0, Form("Ratio: %0.4f #pm %0.4f", ratio, error), "");
      legend->AddEntry((TObject*)0, Form("#chi^{2}/ndf: %0.4f", reports[r][i].chisqrndf), "");
      legend->Draw();

      bgfit[i] = new TF1(Form("bgfit_%d",i),fits::g_p2fit,hcalfit_l,hcalfit_h,3);
      //cout << "Background fit parameter" << endl;
      for (int j=0; j<3; ++j){
	bgfit[i]->SetParameter(j,reports[r][i].fitParams[j+4]);
	//cout << "   " << j << " = " << reports[r][i].fitParams[j+4] << endl;
      }

      bgfit[i]->SetLineWidth(2);
      bgfit[i]->Draw("same");

      TF1 *mcfit = new TF1(Form("mcfit_%d_%d", r, i), fitSlices_p2, hcalfit_l, hcalfit_h, 7);
      mcfit->SetParameters(reports[r][i].fitParams);
      // Can set the error here, might be relevant later
      // mcfit->SetParErrors(reports[r][i].fitErrors);

      // Create a new TH1D or fill an existing one with the values from TF1
      TH1D* hFromTF1 = new TH1D(Form("hFromTF1_%d_%d", r, i), "Histogram from TF1", nbins, x_low, x_high);

      for (int bin = 1; bin <= nbins; ++bin) {
	double binCenter = reports[r][i].sliceHistogram->GetBinCenter(bin);
	double funcValue = mcfit->Eval(binCenter);
	hFromTF1->SetBinContent(bin, funcValue);
      }

      int transparentGreen = TColor::GetColorTransparent(kGreen, 0.3);
      hFromTF1->SetLineColor(kGreen);
      hFromTF1->SetLineWidth(0);
      hFromTF1->SetFillColor(transparentGreen);
      hFromTF1->SetFillStyle(1001);
      hFromTF1->Draw("SAME LF2");

      canvasSlices->Update();
    }//endloop over N (slices)

    canvasSlices->Write();
  }//endloop over r (cuts)

  // Now, let's create and draw TGraphErrors for each cut
  for (int r = 0; r < iteridx; ++r) {
    // The number of valid points for this cut might be different from N if some were skipped
    int numValidPoints = Rsf_vec[r].size();
    std::string cutname = RemovePrefix(reports[r][0].type,"hist_");

    // Creating a TGraphErrors for the current cut
    TGraphErrors* graphErrors = new TGraphErrors(numValidPoints);
    graphErrors->SetTitle(Form("Ratio vs slice midpoint for cut %s;%s;R_{sf}", cutname.c_str(), cutname.c_str()));

    // Filling the TGraphErrors with data
    for (int i = 0; i < numValidPoints; ++i) {
      graphErrors->SetPoint(i, xval_vec[r][i], Rsf_vec[r][i]);
      graphErrors->SetPointError(i, 0, Rsferr_vec[r][i]); // Assuming no error in x
    }

    // Set some graphical attributes
    graphErrors->SetMarkerStyle(21);
    graphErrors->SetMarkerColor(kBlue);
    graphErrors->SetLineColor(kBlue);

    // Drawing the graph
    TCanvas* graphCanvas = new TCanvas(Form("GraphCanvas_%d", r), Form("Ratio vs slice midpoint for cut %s", cutname.c_str()), 1600, 600);
    graphCanvas->cd();

    double minY = 0.7; // Minimum y-axis value
    double maxY = 1.2; // Maximum y-axis value
    graphErrors->SetMinimum(minY);
    graphErrors->SetMaximum(maxY);

    graphErrors->Draw("AP"); // Draw with markers and a line connecting points

    // Save or Write the canvas as needed
    // graphCanvas->SaveAs(Form("TGraphErrors_Cut_%d.png", r)); // Save to a file
    graphCanvas->Write(); // Or write to an open ROOT file
  }

  // Create a new canvas to compare the contributions
  TCanvas* comparisonCanvas = new TCanvas("comparisonCanvas", "Systematic Error Contributions", 800, 600);
  comparisonCanvas->Divide(1, 1); // Adjust based on how you want to layout the histograms/text
  comparisonCanvas->cd();

  // Prepare a TPaveText to display the contributions
  TPaveText* contributionsText = new TPaveText(0.1, 0.1, 0.9, 0.9, "NDC");
  contributionsText->SetTextAlign(12); // Left align
  contributionsText->SetBorderSize(1);
  contributionsText->SetFillStyle(0); // Transparent

  for (size_t i = 0; i < relativeErrors.size(); ++i) {
    char contributionText[255];
    sprintf(contributionText, "%s: e^{2} = %.4f", relativeErrorSources[i].c_str(), relativeErrors[i]);
    contributionsText->AddText(contributionText);
  }

  // Add a blank line for spacing before the total
  contributionsText->AddText(" ");

  //////////////////////////////////////
  //Rsf and error contributions with p2 bg model
  double Telastics = elastics_p2;
  double pscale_bestfit = shiftp2par_vector[0].first;
  double nscale_bestfit = shiftp2par_vector[1].first;
  double sum_of_par_err_squares = 0.0;
  // Assuming bg and shift pars equally and independently effect Rsf
  // for ( size_t p=0; p<shiftp2par_vector.size(); ++p) {
  //   double term = shiftp2par_vector[p].second/shiftp2par_vector[p].first;

  //   sum_of_par_err_squares += term * term;
  //   cout << shiftp2par_vector[p].first << " " << term << endl;
  // }

  // Assuming only the scale parameters meaningfully and independently effect Rsf
  for ( size_t p=0; p<2; ++p) {
    double term = shiftp2par_vector[p].second/shiftp2par_vector[p].first;

    sum_of_par_err_squares += term * term;
    //cout << shiftp2par_vector[p].first << " " << term << endl;
  }
  //////////////////////////////////////

  // // Overall error calculation where the errors are independent and uncorrelated
  // double sum_of_squares = 0.0;
  // for (double error : relativeErrors) {
  //   sum_of_squares += error * error;
  // }
  // Overall error calculation where the errors are independent and uncorrelated
  double sum_of_squares = 0.0;
  for (double error : relativeErrors) {
    sum_of_squares += error;
  }

  double totalSystematicError = std::sqrt(sum_of_squares);

  double totalStatisticalError = sqrt(1/Telastics);

  double Rsf_reported = nscale_bestfit / pscale_bestfit;

  double totalFitError = Rsf_reported * std::sqrt(sum_of_par_err_squares);

  double totalError = std::sqrt(totalSystematicError*totalSystematicError + totalStatisticalError*totalStatisticalError + totalFitError*totalFitError);

  // Add a line for the total systematic error to the TPaveText
  char totalSystematicErrorText[255];
  sprintf(totalSystematicErrorText, "Estimated Systematic Error: #sigma_{syst} = %.4f", totalSystematicError);
  contributionsText->AddText(totalSystematicErrorText);

  // Add a line for the total statistical error to the TPaveText
  char totalStatisticalErrorText[255];
  sprintf(totalStatisticalErrorText, "Estimated Statistical Error: #sigma_{stat} = %.4f", totalStatisticalError);
  contributionsText->AddText(totalStatisticalErrorText);

  // Add a line for the total statistical error to the TPaveText
  char totalFitErrorText[255];
  sprintf(totalFitErrorText, "Estimated Fit Error: #sigma_{fit} = %.4f", totalFitError);
  contributionsText->AddText(totalFitErrorText);

  // Add a line for the total error to the TPaveText
  char totalErrorText[255];
  sprintf(totalErrorText, "Estimated Total Quadrature Error: #sigma_{total} = %.4f", totalError);
  contributionsText->AddText(totalErrorText);

  // Add a blank line for spacing before the overall
  contributionsText->AddText(" ");

  // Add a line for the total error to the TPaveText
  char totalText[255];
  sprintf(totalText, "R_{sf} = %.4f #pm %.4f", Rsf_reported, totalError);
  contributionsText->AddText(totalText);

  //Use extract to get GMn and add error from model dependent extractions
  auto GMn_and_error = extract::extract_GMn_from_simc( Rsf_reported,
						       totalError,
						       0.,              //Q2 will be calculated when 0.
						       E_beam,
						       BB_angle,
						       true);

  
  // Add a blank line for spacing before the overall
  contributionsText->AddText(" ");

  // Add a line for proton model RCS to the TPaveText
  char prcsText[255];
  sprintf(prcsText, "Arrington07 proton RCS = %.4f #pm %.4f", GMn_and_error[4].first, GMn_and_error[4].second);
  contributionsText->AddText(prcsText);

  // Add a line for neutron model-dependent RCS to the TPaveText
  char nrcsText[255];
  sprintf(nrcsText, "Model-dependent neutron RCS = %.4f #pm %.4f", GMn_and_error[3].first, GMn_and_error[3].second);
  contributionsText->AddText(nrcsText);
  
  // Add a blank line for spacing before the overall
  contributionsText->AddText(" ");

  // Add a line for gmn extracted to the TPaveText
  char gmnText[255];
  sprintf(gmnText, "G_{M}^{n} = %.4f #pm %.4f", GMn_and_error[0].first, GMn_and_error[0].second);
  contributionsText->AddText(gmnText);

  // Add a line for gmn extracted and normalized to the TPaveText
  char gmncText[255];
  sprintf(gmncText, "G_{M}^{n}/G_{D}/#mu_{n} = %.4f #pm %.4f", GMn_and_error[1].first, GMn_and_error[1].second);
  contributionsText->AddText(gmncText);

  // Draw the contributionsText on the canvas
  contributionsText->Draw();
  comparisonCanvas->Update();
  comparisonCanvas->Write();
  
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


  outputFile->Write();
  //outputFile->Close();

  cout << "All plots created. Output file located here: " << foutPath << endl;

}

/////////////
/////////////
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

void createAndSaveCanvasWithFitsAndResiduals(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, TF1* bg, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, bool shiftfit, double &elastics) {

  //Clone the dx histo for bg sub
  TH1D *hdx_histplot = (TH1D*)(hdx->Clone("hdx_histplot"));
  TH1D *hdx_clone = (TH1D*)(hdx->Clone("hdx_clone"));
  TH1D *hdx_p_clone = (TH1D*)(hdxp->Clone("hdx_p_clone"));
  TH1D *hdx_n_clone = (TH1D*)(hdxn->Clone("hdx_n_clone"));

  //Recreate the fit
  TF1* fit;
  if(shiftfit)
    fit = new TF1("fit", fitFullShift, hcalfit_l, hcalfit_h, 9);
  else
    fit = new TF1("fit", fitFull, hcalfit_l, hcalfit_h, 7);
  
  for (size_t i = 0; i < pars.size(); ++i) {
    fit->SetParameter(i, pars[i].first);
    fit->SetParError(i, pars[i].second);
  }

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
  hdx_histplot->SetLineColor(kBlack);
  hdx_histplot->SetLineWidth(1);
  hdx_histplot->SetTitle(Form("dx, %s;m",type));
  hdx_histplot->Draw("hist");
  hdx->SetLineColor(kBlack);
  hdx->SetLineWidth(1);
  hdx->Draw("E same");
  fit->SetLineColor(kGreen);
  fit->SetLineWidth(1);
  fit->Draw("same");
  if(shiftfit)
    hdx_p_clone = util::shiftHistogramX(hdx_p,pars[2].first);
  hdx_p_clone->Scale(pars[0].first);
  hdx_p_clone->SetLineColor(kRed-5);
  hdx_p_clone->SetLineWidth(2);
  hdx_p_clone->Draw("E same");
  if(shiftfit)
    hdx_n_clone = util::shiftHistogramX(hdx_n,pars[3].first);
  hdx_n_clone->Scale(pars[1].first);
  hdx_n_clone->SetLineColor(kBlue-5);
  hdx_n_clone->SetLineWidth(2);
  hdx_n_clone->Draw("E same");

  // Set background function
  bg->SetLineColor(kBlack);
  bg->SetFillColorAlpha(kGray,0.35);
  bg->SetFillStyle(1001);
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

  // Construct legend within function
  TLegend* leg = new TLegend(0.5, 0.5, 0.89, 0.89);
  leg->AddEntry(hdx, "Data", "l");
  leg->AddEntry(fit, "Fit", "l");
  leg->AddEntry(bg, "Background", "l");
  leg->AddEntry(hdx_p_clone, "Proton SIMC MC", "l");
  leg->AddEntry(hdx_n_clone, "Neutron SIMC MC", "l");
  leg->AddEntry( (TObject*)0, "", "");
  leg->AddEntry( (TObject*)0, Form("data N events : %0.0f",hdx->GetEntries()), "");
  //leg->AddEntry( (TObject*)0, Form("MC N events : %0.0f",(hdx_p_clone->GetEntries()+hdx_n_clone->GetEntries())), "");
  // leg->AddEntry( (TObject*)0, Form("n/p scale ratio R_{sf} : %0.3f #pm %0.3f (fit error)",np_par_ratio, np_par_ratio_error), "");
  // leg->AddEntry( (TObject*)0, Form("n/p yield ratio R'' : %0.3f #pm %0.3f (fit and sum statistical error)",np_sum_ratio, np_sum_ratio_error), "");
  leg->AddEntry( (TObject*)0, Form("n/p scale ratio R_{sf} : %0.3f #pm %0.3f",np_par_ratio, np_par_ratio_error), "");
  leg->AddEntry( (TObject*)0, Form("n/p yield ratio R'' : %0.3f #pm %0.3f",np_sum_ratio, np_sum_ratio_error), "");
  if(shiftfit)
    //leg->AddEntry( (TObject*)0, Form("dx shift pars, n/p : %0.3f / %0.3f ", pars[2].first, pars[3].first), "");
  // leg->AddEntry( (TObject*)0, Form("#chi^{2}/ndf: %0.3f",fit->GetChisquare()/fit->GetNDF()), "");
  leg->AddEntry( (TObject*)0, Form("#chi^{2}/ndf: %0.3f/%d",fit->GetChisquare(),fit->GetNDF()), "");
  leg->Draw("same");

  // Construct residuals
  //add the p and n histograms together for comparison
  TH1D* sumHistogram = new TH1D("sumHistogram", "Sum of Histograms", hdx_p_clone->GetNbinsX(), hdx_p_clone->GetXaxis()->GetXmin(), hdx_p_clone->GetXaxis()->GetXmax());
  sumHistogram->Add(hdx_p_clone, hdx_n_clone);
  //sumHistogram->SetName(Form("sh_%s",savePath)); //stop the warning messages

  //Background subraction and total elastic statistics
  elastics = 0;
  for (int bin = 1; bin <= hdx->GetNbinsX(); ++bin) {
    double bgValue = bg->Eval(hdx_clone->GetXaxis()->GetBinCenter(bin));
    if( hdx_clone->GetBinContent(bin) > bgValue )
      elastics += hdx_clone->GetBinContent(bin) - bgValue;
    hdx_clone->SetBinContent(bin, hdx_clone->GetBinContent(bin) - bgValue);
    hdx_clone->SetBinError(bin, hdx_clone->GetBinError(bin));
  }

  TH1D *hRes = util::makeResidualHistoWithError("dx",hdx_clone,sumHistogram,true,false);
  TH1D *hRes_histplot = (TH1D*)(hRes->Clone("hRes_histplot"));

  hRes_histplot->SetTitle("");
  //hRes->SetName(Form("hr_%s",savePath));
  hRes_histplot->GetXaxis()->SetRangeUser(hcalr_l,hcalr_h);

  double dxMaxValue = hdx_clone->GetMaximum();
  hRes->GetYaxis()->SetRangeUser(-dxMaxValue/4,dxMaxValue/4);
  hRes_histplot->GetYaxis()->SetRangeUser(-dxMaxValue/4,dxMaxValue/4);

  hRes_histplot->SetLineColor(kRed);
  hRes_histplot->SetLineWidth(2);
  hRes->SetLineColor(kBlack);
  hRes->SetMarkerStyle(1);
  hRes->SetLineWidth(1);

  // Second pad: Draw the residuals
  pad2->cd();
  hRes_histplot->Draw("hist");
  hRes_histplot->GetYaxis()->SetTitle("Residuals");
  hRes_histplot->GetYaxis()->SetTitleSize(0.07);
  hRes_histplot->GetYaxis()->SetTitleOffset(0.5);
  hRes_histplot->GetYaxis()->SetLabelSize(0.07);
  hRes_histplot->GetXaxis()->SetTitleSize(0.07);
  hRes_histplot->GetXaxis()->SetLabelSize(0.07);
  hRes_histplot->GetXaxis()->SetTitle("x_{hcal}-x_{exp}");
  hRes->Draw("E same");

  // Draw a line at y=0 for residuals
  TF1 *zeroLine = new TF1("zeroLine", "0", hRes->GetXaxis()->GetXmin(), hRes->GetXaxis()->GetXmax());
  zeroLine->SetLineColor(kBlack);
  zeroLine->SetLineStyle(2);
  zeroLine->Draw("same");

  // Update the canvas
  cTotal->cd();
  cTotal->Update();

  // Save the canvas to a file
  cTotal->SaveAs(savePath);
}

//Create canvas for reporting
void createAndSaveCanvasWithFitsAndResiduals_nobg(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, bool shiftfit, double &elastics) {

  //Clone the dx histos
  TH1D *hdx_histplot = (TH1D*)(hdx->Clone("hdx_histplot"));
  TH1D *hdx_clone = (TH1D*)(hdx->Clone("hdx_clone"));
  TH1D *hdx_p_clone = (TH1D*)(hdxp->Clone("hdx_p_clone"));
  TH1D *hdx_n_clone = (TH1D*)(hdxn->Clone("hdx_n_clone"));

  //Recreate the fit
  TF1* fit;
  if(shiftfit)
    fit = new TF1("fit", fitFullShift_nobg, hcalfit_l, hcalfit_h, 4);
  else
    fit = new TF1("fit", fitFull_nobg, hcalfit_l, hcalfit_h, 2);
  
  for (size_t i = 0; i < pars.size(); ++i) {
    fit->SetParameter(i, pars[i].first);
    fit->SetParError(i, pars[i].second);
  }

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
  hdx_histplot->SetLineColor(kBlack);
  hdx_histplot->SetLineWidth(1);
  hdx_histplot->SetTitle(Form("dx, %s;m",type));
  hdx_histplot->Draw("hist");
  hdx->SetLineColor(kBlack);
  hdx->SetLineWidth(1);
  hdx->Draw("E same");
  fit->SetLineColor(kGreen);
  fit->SetLineWidth(1);
  fit->Draw("same");
  if(shiftfit)
    hdx_p_clone = util::shiftHistogramX(hdx_p,pars[2].first);
  hdx_p_clone->Scale(pars[0].first);
  hdx_p_clone->SetLineColor(kRed-5);
  hdx_p_clone->SetLineWidth(2);
  hdx_p_clone->Draw("E same");
  if(shiftfit)
    hdx_n_clone = util::shiftHistogramX(hdx_n,pars[3].first);
  hdx_n_clone->Scale(pars[1].first);
  hdx_n_clone->SetLineColor(kBlue-5);
  hdx_n_clone->SetLineWidth(2);
  hdx_n_clone->Draw("E same");

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

  // Construct legend within function
  TLegend* leg = new TLegend(0.45, 0.5, 0.89, 0.89);
  leg->AddEntry(hdx, "Data", "l");
  leg->AddEntry(fit, "Fit", "l");
  leg->AddEntry(hdx_p_clone, "Proton SIMC MC", "l");
  leg->AddEntry(hdx_n_clone, "Neutron SIMC MC", "l");
  leg->AddEntry( (TObject*)0, "", "");
  leg->AddEntry( (TObject*)0, Form("data N events : %0.0f",hdx->GetEntries()), "");
  //leg->AddEntry( (TObject*)0, Form("MC N events : %0.0f",(hdx_p_clone->GetEntries()+hdx_n_clone->GetEntries())), "");
  // leg->AddEntry( (TObject*)0, Form("n/p scale ratio R_{sf} : %0.3f #pm %0.3f (fit error)",np_par_ratio, np_par_ratio_error), "");
  // leg->AddEntry( (TObject*)0, Form("n/p yield ratio R'' : %0.3f #pm %0.3f (fit and sum statistical error)",np_sum_ratio, np_sum_ratio_error), "");
  // if(shiftfit)
  leg->AddEntry( (TObject*)0, Form("n/p scale ratio R_{sf} : %0.3f #pm %0.3f",np_par_ratio, np_par_ratio_error), "");
  leg->AddEntry( (TObject*)0, Form("n/p yield ratio R'' : %0.3f #pm %0.3f",np_sum_ratio, np_sum_ratio_error), "");
  //leg->AddEntry( (TObject*)0, Form("dx shift pars, n/p : %0.3f / %0.3f ", pars[2].first, pars[3].first), "");
  // leg->AddEntry( (TObject*)0, Form("#chi^{2}/ndf: %0.3f",fit->GetChisquare()/fit->GetNDF()), "");
  leg->AddEntry( (TObject*)0, Form("#chi^{2}/ndf: %0.3f/%d",fit->GetChisquare(),fit->GetNDF()), "");
  leg->Draw("same");

  // Construct residuals
  //add the p and n histograms together for comparison
  TH1D* sumHistogram = new TH1D("sumHistogram", "Sum of Histograms", hdx_p_clone->GetNbinsX(), hdx_p_clone->GetXaxis()->GetXmin(), hdx_p_clone->GetXaxis()->GetXmax());
  sumHistogram->Add(hdx_p_clone, hdx_n_clone);
  //sumHistogram->SetName(Form("nobg_sh_%s",savePath)); //stop the warning messages

  //Total elastic statistics
  elastics = 0;
  for (int bin = 1; bin <= hdx->GetNbinsX(); ++bin)
    elastics += hdx_clone->GetBinContent(bin);
  
  TH1D *hRes = util::makeResidualHistoWithError("dx",hdx_clone,sumHistogram,true,false);
  TH1D *hRes_histplot = (TH1D*)(hRes->Clone("hRes_histplot"));

  hRes_histplot->SetTitle("");
  //hRes->SetName(Form("hr_%s",savePath));
  hRes_histplot->GetXaxis()->SetRangeUser(hcalr_l,hcalr_h);

  double dxMaxValue = hdx_clone->GetMaximum();
  hRes->GetYaxis()->SetRangeUser(-dxMaxValue/4,dxMaxValue/4);
  hRes_histplot->GetYaxis()->SetRangeUser(-dxMaxValue/4,dxMaxValue/4);

  hRes_histplot->SetLineColor(kRed);
  hRes_histplot->SetLineWidth(2);
  hRes->SetMarkerStyle(1);
  hRes->SetLineColor(kBlack);
  hRes->SetLineWidth(1);

  // Second pad: Draw the residuals
  pad2->cd();
  hRes_histplot->Draw("hist");
  hRes_histplot->GetYaxis()->SetTitle("Residuals");
  hRes_histplot->GetYaxis()->SetTitleSize(0.07);
  hRes_histplot->GetYaxis()->SetTitleOffset(0.5);
  hRes_histplot->GetYaxis()->SetLabelSize(0.07);
  hRes_histplot->GetXaxis()->SetTitleSize(0.07);
  hRes_histplot->GetXaxis()->SetLabelSize(0.07);
  hRes_histplot->GetXaxis()->SetTitle("x_{hcal}-x_{exp}");
  hRes->Draw("E same");

  // Draw a line at y=0 for residuals
  TF1 *zeroLine = new TF1("zeroLine", "0", hRes->GetXaxis()->GetXmin(), hRes->GetXaxis()->GetXmax());
  zeroLine->SetLineColor(kBlack);
  zeroLine->SetLineStyle(2);
  zeroLine->Draw("same");

  // Update the canvas
  cTotal->cd();
  cTotal->Update();

  // Save the canvas to a file
  cTotal->SaveAs(savePath);
}

//Create canvas for reporting, expects a that the peaks are allowed to slide
void createAndSaveCanvasWithFitsAndResiduals_altbg(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, TH1D* hdxaltbg, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, std::string bgtype, double &elastics) {

  //Clone the dx histos
  TH1D *hdx_histplot = (TH1D*)(hdx->Clone("hdx_histplot"));
  TH1D *hdx_clone = (TH1D*)(hdx->Clone("hdx_clone"));
  TH1D *hdx_p_clone = (TH1D*)(hdxp->Clone("hdx_p_clone"));
  TH1D *hdx_n_clone = (TH1D*)(hdxn->Clone("hdx_n_clone"));
  //TH1D *hdx_altbg_clone = (TH1D*)(hdxaltbg->Clone("hdx_altbg_clone"));

  //Recreate the fit
  TF1* fit;
  if(bgtype.compare("inel")==0){
    fit = new TF1("fit", fitFullShiftInel, hcalfit_l, hcalfit_h, 5);
  }else if(bgtype.compare("dy")==0){
    fit = new TF1("fit", fitFullShiftDyAnti, hcalfit_l, hcalfit_h, 5);
  }else if(bgtype.compare("coin")==0){
    fit = new TF1("fit", fitFullShiftCoinAnti, hcalfit_l, hcalfit_h, 5);
  }else{
    cout << "ERROR: Invalid MC background type. Must be inel, antidy, or anticoin." << endl;
    return;
  }

  for (size_t i = 0; i < pars.size(); ++i) {
    fit->SetParameter(i, pars[i].first);
    fit->SetParError(i, pars[i].second);
  }

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
  hdx_histplot->SetLineColor(kBlack);
  hdx_histplot->SetLineWidth(1);
  hdx_histplot->SetTitle(Form("dx, %s;m",type));
  hdx_histplot->Draw("hist");
  hdx->SetLineColor(kBlack);
  hdx->SetLineWidth(1);
  hdx->Draw("E same");
  fit->SetLineColor(kGreen);
  fit->SetLineWidth(1);
  fit->Draw("same");
  hdx_p_clone = util::shiftHistogramX(hdx_p,pars[2].first);
  hdx_p_clone->Scale(pars[0].first);
  hdx_p_clone->SetLineColor(kRed-5);
  hdx_p_clone->SetLineWidth(2);
  hdx_p_clone->Draw("E same");
  hdx_n_clone = util::shiftHistogramX(hdx_n,pars[3].first);
  hdx_n_clone->Scale(pars[1].first);
  hdx_n_clone->SetLineColor(kBlue-5);
  hdx_n_clone->SetLineWidth(2);
  hdx_n_clone->Draw("E same");

  // Set background function
  TH1D *hdx_altbg_clone = util::shiftHistogramX(hdxaltbg,pars[3].first);
  hdx_altbg_clone->Scale(pars[4].first);
  hdx_altbg_clone->SetLineColor(kBlack);
  hdx_altbg_clone->SetFillColorAlpha(kGray,0.35);
  hdx_altbg_clone->SetFillStyle(1001);
  hdx_altbg_clone->Draw("same");

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

  // Construct legend within function
  TLegend* leg = new TLegend(0.5, 0.5, 0.89, 0.89);
  leg->AddEntry( hdx, "Data", "l" );
  leg->AddEntry( fit, "Fit", "l" );
  leg->AddEntry( hdx_p_clone, "Proton SIMC MC", "l" );
  leg->AddEntry( hdx_n_clone, "Neutron SIMC MC", "l" );
  leg->AddEntry( hdx_altbg_clone, Form("BG, %s",bgtype.c_str()), "l" );
  leg->AddEntry( (TObject*)0, "", "");
  leg->AddEntry( (TObject*)0, Form("data N events : %0.0f",hdx->GetEntries()), "");
  //leg->AddEntry( (TObject*)0, Form("MC N events : %0.0f",(hdx_p_clone->GetEntries()+hdx_n_clone->GetEntries())), "");
  // leg->AddEntry( (TObject*)0, Form("n/p scale ratio R_{sf} : %0.3f #pm %0.3f (fit error)",np_par_ratio, np_par_ratio_error), "");
  // leg->AddEntry( (TObject*)0, Form("n/p yield ratio R'' : %0.3f #pm %0.3f (fit and sum statistical error)",np_sum_ratio, np_sum_ratio_error), "");
  leg->AddEntry( (TObject*)0, Form("n/p scale ratio R_{sf} : %0.3f #pm %0.3f",np_par_ratio, np_par_ratio_error), "");
  leg->AddEntry( (TObject*)0, Form("n/p yield ratio R'' : %0.3f #pm %0.3f",np_sum_ratio, np_sum_ratio_error), "");
  //leg->AddEntry( (TObject*)0, Form("dx shift pars, n/p : %0.3f / %0.3f ", pars[2].first, pars[3].first), "");
  leg->AddEntry( (TObject*)0, Form("bg scale: %0.3f", pars[4].first), "");
  //leg->AddEntry( (TObject*)0, Form("#chi^{2}/ndf: %0.3f",fit->GetChisquare()/fit->GetNDF()), "");
  leg->AddEntry( (TObject*)0, Form("#chi^{2}/ndf: %0.3f/%d",fit->GetChisquare(),fit->GetNDF()), "");
  leg->Draw("same");

  // Construct residuals
  //add the p and n histograms together for comparison
  TH1D* sumHistogram = new TH1D("sumHistogram", "Sum of Histograms", hdx_p_clone->GetNbinsX(), hdx_p_clone->GetXaxis()->GetXmin(), hdx_p_clone->GetXaxis()->GetXmax());
  sumHistogram->Add(hdx_p_clone, hdx_n_clone);
  //sumHistogram->SetName(Form("nobg_sh_%s",savePath)); //stop the warning messages

  //subtract the background and get elastic stats
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

  TH1D *hRes = util::makeResidualHistoWithError("dx",hdx_clone,sumHistogram,true,false);
  TH1D *hRes_histplot = (TH1D*)(hRes->Clone("hRes_histplot"));

  hRes_histplot->SetTitle("");
  //hRes->SetName(Form("hr_%s",savePath));
  hRes_histplot->GetXaxis()->SetRangeUser(hcalr_l,hcalr_h);

  double dxMaxValue = hdx_clone->GetMaximum();
  hRes->GetYaxis()->SetRangeUser(-dxMaxValue/4,dxMaxValue/4);
  hRes_histplot->GetYaxis()->SetRangeUser(-dxMaxValue/4,dxMaxValue/4);

  hRes_histplot->SetLineColor(kRed);
  hRes_histplot->SetLineWidth(2);
  hRes->SetMarkerStyle(1);
  hRes->SetLineColor(kBlack);
  hRes->SetLineWidth(1);

  // Second pad: Draw the residuals
  pad2->cd();
  hRes_histplot->Draw("hist");
  hRes_histplot->GetYaxis()->SetTitle("Residuals");
  hRes_histplot->GetYaxis()->SetTitleSize(0.07);
  hRes_histplot->GetYaxis()->SetTitleOffset(0.5);
  hRes_histplot->GetYaxis()->SetLabelSize(0.07);
  hRes_histplot->GetXaxis()->SetTitleSize(0.07);
  hRes_histplot->GetXaxis()->SetLabelSize(0.07);
  hRes_histplot->GetXaxis()->SetTitle("x_{hcal}-x_{exp}");
  hRes->Draw("E same");

  // Draw a line at y=0 for residuals
  TF1 *zeroLine = new TF1("zeroLine", "0", hRes->GetXaxis()->GetXmin(), hRes->GetXaxis()->GetXmax());
  zeroLine->SetLineColor(kBlack);
  zeroLine->SetLineStyle(2);
  zeroLine->Draw("same");

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
  histogram->Fit(fineFit, fitOptions.c_str());

  // Update parameters and errors with fine fit results
  for (int i = 0; i < paramCount; ++i) {
    parametersAndErrors[i].first = fineFit->GetParameter(i); // Fine fit parameter value
    parametersAndErrors[i].second = fineFit->GetParError(i); // Fine fit parameter error
  }

  fitqual.first = fineFit->GetChisquare();
  fitqual.second = fineFit->GetNDF();

  delete fit; // Delete fit to avoid carry-over
  delete fineFit; // Clean up
  return parametersAndErrors;
}

//Function to fit then fine fit and return all fit parameters
std::vector<std::pair<double, double>> fitAndFineFit_slice(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, double pshift, double nshift, const std::string& fitOptions = "RBMQ0") {
  TF1* fit = new TF1(fitName.c_str(), fitFormula.c_str(), hcalfit_l, hcalfit_h, paramCount);
  //fit->SetNpx(5000);
  for (int i=0; i<paramCount; ++i){ //reset parameters/errors for this set
    fit->SetParameter(i,0);
    fit->SetParError(i,0);
  }
  fit->FixParameter(2,pshift); //fix x shift for proton
  fit->FixParameter(3,nshift); //fix x shift for neutron

  //Make certain that the bg fit is concave down when pol2 bg.
  if(fitFormula.compare("fitSlices_p2")==0){
    fit->SetParLimits(4,0,10e8);
    fit->SetParLimits(5,-10e8,0);
    fit->SetParLimits(6,-10e8,0);
  }

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
  fineFit->FixParameter(2,pshift); //keep fixed x shift
  fineFit->FixParameter(3,nshift); //keep fixed x shift

  //Make certain that the bg fit is concave down when pol2 bg.
  if(fitFormula.compare("fitSlices_p2")==0){
    fineFit->SetParLimits(4,0,10e8);
    fineFit->SetParLimits(5,-10e8,0);
    fineFit->SetParLimits(6,-10e8,0);
  }

  fineFit->SetParameters(fineFitInitialParams.data());
  histogram->Fit(fineFit, fitOptions.c_str());

  //Write out max mc values for verification
  cout << hdx_slice_p->GetMaximum() << ", " << hdx_slice_n->GetMaximum();

  // Update parameters and errors with fine fit results
  for (int i = 0; i < paramCount; ++i) {
    parametersAndErrors[i].first = fineFit->GetParameter(i); // Fine fit parameter value
    parametersAndErrors[i].second = fineFit->GetParError(i); // Fine fit parameter error
  }

  fitqual.first = fineFit->GetChisquare();
  fitqual.second = fineFit->GetNDF();

  delete fit; // Delete fit to avoid carry-over
  delete fineFit; // Clean up
  return parametersAndErrors;
}

//simply background subtraction function
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
