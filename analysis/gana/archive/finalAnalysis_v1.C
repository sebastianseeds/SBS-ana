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
// double hcalfit_l_or = -2.1; //lower fit/bin limit for hcal dx plots (m) sbs4, 50p
// double hcalfit_h_or = 0.7; //upper fit/bin limit for hcal dx plots (m) sbs4, 50p
double hcalfit_l_or = -1.7; //lower fit/bin limit for hcal dx plots (m) sbs4, 30p
double hcalfit_h_or = 0.7; //upper fit/bin limit for hcal dx plots (m) sbs4, 30p

// double hcalfit_l_or = econst::hcalposXi_mc; //lower fit/bin limit for hcal dx plots (m)
// double hcalfit_h_or = econst::hcalposXf_mc; //upper fit/bin limit for hcal dx plots (m)
//double hcalfit_l_or = -3.0; //lower fit/bin limit for hcal dx plots (m)
//double hcalfit_h_or = 2.0; //upper fit/bin limit for hcal dx plots (m)

//Plot range option
// double hcalr_l = -2.1; //lower fit/bin limit for hcal dx plots (m) sbs4, 50p
// double hcalr_h = 0.7; //upper fit/bin limit for hcal dx plots (m) sbs4, 50p
double hcalr_l = -1.7; //lower fit/bin limit for hcal dx plots (m) sbs4, 30p
double hcalr_h = 0.7; //upper fit/bin limit for hcal dx plots (m) sbs4, 30p

// double hcalr_l = econst::hcalposXi_mc; //lower fit/bin limit for hcal dx plots (m)
// double hcalr_h = econst::hcalposXf_mc; //upper fit/bin limit for hcal dx plots (m)

//Fit options
std::string fitopt = "RMQ0";
//std::string fitopt = "RBMQ0"; //R:setrange,B:predefinedTF1,M:improvefit,Q:quiet,0:nofitline
//std::string fitopt = "RLQ0"; //L:loglikelihood for counts
//std::string fitopt = "RLEMQ0"; //E:bettererrorest
//std::string fitopt = "RWLEMQ0"; //WL:weightedloglikelihood for weighted counts (N/A here)
//std::string fitopt = "RPEMQ0"; //P:pearsonloglikelihood for expected error (N/A here)

//fit min entries
int minEntries = 500;

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
					"hist_fiducial_sig_x"};

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
void createAndSaveCanvasWithFitsAndResiduals(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, TF1* bg, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, bool shiftfit);
void createAndSaveCanvasWithFitsAndResiduals_nobg(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, bool shiftfit);
void createAndSaveCanvasWithFitsAndResiduals_altbg(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, TH1D* hdxinel, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, std::string bgtype);
//std::vector<double> fitAndFineFit(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, const std::string& fitOptions);
std::vector<std::pair<double, double>> fitAndFineFit(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, const std::string& fitOptions = "RBMQ0");
void subtractBackground(TH1D* hist, TF1* bgFit);
int FindCutIndex(const std::string& branches, const std::string& searchStr);
std::string RemovePrefix(const std::string& originalString, const std::string& prefix);

//shiftX, SBS-4 = 0.05
//Side-band pol4 fit
//sbs9: -1.8 to 0.8
//sbs4 30p: -1.6 to 0.7, shiftX=0.0, neutronshift=-0.05
//sbs4 50p: -2.0 to 0.7, shiftX=0.05, neutronshift=0.0

//main. kine=kinematic, mag=fieldsetting, pass=pass#, sb_min/max=sidebandlimits, shiftX=shifttodxdata, N=cutvarsliceN
void finalAnalysis(int kine=4, int mag=50, int pass=2, double sb_min=-2.0, double sb_max=0.7, double shiftX=0.0, double neutronshift=0.0, int N=6, bool detail=true, bool mc_slices=true, bool blind=false, bool bestclus=true, bool thin=true, bool fitrangeoverride=true, bool alt = true) {
  //set draw params
  gStyle->SetPalette(55);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);

  //set sideband
  SBpol4rej_b = sb_min; 
  SBpol4rej_e = sb_max;

  //load json file
  JSONManager *jmgr = new JSONManager("../../config/syst.json");

  //Get dx plot details
  hcalfit_l = jmgr->GetValueFromSubKey<int>( "hcalfit_l", Form("sbs%d",kine) ); //gets to 7cm HCal pos res
  hcalfit_h = jmgr->GetValueFromSubKey<int>( "hcalfit_h", Form("sbs%d",kine) ); //gets to 7cm HCal pos res

  if(fitrangeoverride){
    hcalfit_l = hcalfit_l_or;
    hcalfit_h = hcalfit_h_or;
  }

  //Get tight elastic cut strings and limits
  std::string branches = jmgr->GetValueFromKey_str( Form("post_branches_p%d",pass) );

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

  // for (size_t i=0; i<llims.size(); ++i)
  //   cout << llims[i] << " " << ulims[i] << endl;

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
  createAndSaveCanvasWithFitsAndResiduals(hdx_shifted_clone, hdx_p_clone, hdx_n_clone, bg_fullFit, fullpar_vector, fullQual, "pol4 BG", "~/gmn_plots/combined_dx_fits_residuals_shift.pdf", false);
  
  //pol4 fit to bg
  createAndSaveCanvasWithFitsAndResiduals(hdx_raw_clone, hdx_p_clone, hdx_n_clone, bg_shiftFit, shiftpar_vector, shiftQual, "shiftfit pol4 BG", "~/gmn_plots/combined_dx_fits_residuals_fitshift.pdf", true);

  //pol2 fit to bg
  createAndSaveCanvasWithFitsAndResiduals(hdx_raw_p2_clone, hdx_p_clone, hdx_n_clone, bg_shiftp2Fit, shiftp2par_vector, shiftQual, "shiftfit pol2 BG", "~/gmn_plots/combined_dx_fits_residuals_fitshiftp2.pdf", true);
  
  //inel bg
  createAndSaveCanvasWithFitsAndResiduals_altbg(hdx_raw_inel_clone, hdx_p_clone, hdx_n_clone, hdx_inel_clone, shiftinel_vector, shiftInelQual, "shiftfit MC BG", "~/gmn_plots/combined_dx_fits_residuals_inel_fitshift.pdf", "inel");
  
  //dy anticut bg
  createAndSaveCanvasWithFitsAndResiduals_altbg(hdx_raw_antidy_clone, hdx_p_clone, hdx_n_clone, hdx_antidy_clone, shiftantidy_vector, shiftAntidyQual, "shiftfit antidy BG", "~/gmn_plots/combined_dx_fits_residuals_antidy_fitshift.pdf", "dy");
  
  //coin anticut bg
  createAndSaveCanvasWithFitsAndResiduals_altbg(hdx_raw_anticoin_clone, hdx_p_clone, hdx_n_clone, hdx_anticoin_clone, shiftanticoin_vector, shiftAnticoinQual, "shiftfit anticoin BG", "~/gmn_plots/combined_dx_fits_residuals_anticoin_fitshift.pdf", "coin");

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
  createAndSaveCanvasWithFitsAndResiduals_nobg(hdx_shifted_nobg, hdx_p_clone, hdx_n_clone, fullpar_nobg_vector, fullnobgQual, "pol4 BG subtracted", "~/gmn_plots/combined_dx_fits_residuals_full_nobg.pdf", false);
  
  createAndSaveCanvasWithFitsAndResiduals_nobg(hdx_raw_nobg, hdx_p_clone, hdx_n_clone, shiftpar_nobg_vector, shiftnobgQual, "shiftfit pol4 BG subtracted", "~/gmn_plots/combined_dx_fits_residuals_fitshift_nobg.pdf", true);
  
  createAndSaveCanvasWithFitsAndResiduals_nobg(hdx_raw_p2_nobg, hdx_p_clone, hdx_n_clone, shiftp2par_nobg_vector, shiftp2nobgQual, "shiftfit pol2 BG subtracted", "~/gmn_plots/combined_dx_fits_residuals_fitshiftp2_nobg.pdf", true);
  
  createAndSaveCanvasWithFitsAndResiduals_nobg(hdx_sb_nobg, hdx_p_clone, hdx_n_clone, sbpar_nobg_vector, sbnobgQual, "shiftfit sideband BG subtracted", "~/gmn_plots/combined_dx_fits_residuals_fitsb_nobg.pdf", true);


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

  //Setup report struct
  struct ReportData {
    TH1D* sliceHistogram = nullptr;
    double pscale = 0;
    double nscale = 0;
    double perr = 0;
    double nerr = 0;
    double winLow = 0;
    double winHigh = 0;
    double nev = 0;
    double scaleratio = 0;
    std::string type;

    // Constructor for easier initialization (optional)
    ReportData(TH1D* hist = nullptr, double ps = 0, double ns = 0, double pe = 0, double ne = 0, 
               double wl = 0, double wh = 0, double nv = 0, double sr = 0, std::string t = "") 
      : sliceHistogram(hist), pscale(ps), nscale(ns), perr(pe), nerr(ne), 
	winLow(wl), winHigh(wh), nev(nv), scaleratio(sr), type(t) {}
  };

  ReportData reports[nregions][reportSamples];

  // Loop over all TH2D histograms
  int iteridx;
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

    // Create a canvas for this TH2D histogram
    TCanvas* th2dCanvas = new TCanvas(Form("canvas_%s", hist->GetName()), hist->GetTitle(), 800, 600);
    th2dCanvas->cd();
    hist->Draw("COLZ");  // Draw the histogram
    th2dCanvas->Write();

    //Determine the first and last bin with data along the x-axis for all data TH2Ds
    int histBE[4] = {0};
    util::checkTH2DBinsAndEntries(hist,minEntries,histBE[0],histBE[1],histBE[2],histBE[3]);

    //Determine the first and last bin with data along the x-axis for all data TH2Ds
    int histBE_mcp[4] = {0};
    util::checkTH2DBinsAndEntries(hist_mc_p,minEntries,histBE_mcp[0],histBE_mcp[1],histBE_mcp[2],histBE_mcp[3]);

    //Determine the first and last bin with data along the x-axis for all data TH2Ds
    int histBE_mcn[4] = {0};
    util::checkTH2DBinsAndEntries(hist_mc_n,minEntries,histBE_mcn[0],histBE_mcn[1],histBE_mcn[2],histBE_mcn[3]);

    //get first bin where sufficient data exists on all histograms if mc_slices enabled
    int first_bin;
    int last_bin;
    if(mc_slices){
      first_bin = std::max({histBE[0],histBE_mcp[0],histBE_mcn[0]});
      last_bin = std::min({histBE[1],histBE_mcp[1],histBE_mcn[1]});
    }else{
      first_bin = histBE[0];
      last_bin = histBE[1];
    }

    //finally check llim and ulim on this histogram and apply appropriate elastic cuts for this SF R value
    std::string prefix = "hist_";
    std::string modHistName = RemovePrefix(histName, prefix);

    int limitIndex = FindCutIndex(branches,modHistName);
    
    if(limitIndex<0){
      cout << "ERROR: limit index < 0" << endl;
      return;
    }

    double cutllim = llims[limitIndex];
    double cutulim = ulims[limitIndex];

    double xLow = hist->GetXaxis()->GetBinLowEdge(first_bin);
    double xHigh = hist->GetXaxis()->GetBinUpEdge(last_bin);

    xLow = std::max(xLow,cutllim);
    xHigh = std::min(xHigh,cutulim);

    // Create TH1D to store ratios and yields with the same x-axis range as the TH2D histogram
    std::string ratioHistName = std::string("ratio: ") + hist->GetName() + ";" + hist->GetName() + ";scale factor n:p ratio";
    std::string ratioHistNameBG = std::string("ratio_bgsub: ") + hist->GetName() + ";" + hist->GetName() + ";scale factor n:p ratio";
    std::string yieldHistName = std::string("yield: ") + hist->GetName() + ";" + hist->GetName() + ";p+n yield";
    std::string ratioHistName_w = std::string("R, weighted proj: ") + hist->GetName() + "; scale ratio";;
    std::string ratioHistNameBG_w = std::string("R_bg, weighted proj: ") + hist->GetName() + "; scale ratio background subtracted";
 
    TH1D* ratioHist = new TH1D(ratioHistName.c_str(), ratioHistName.c_str(), N, xLow, xHigh);
    TH1D* ratioHist_bg = new TH1D(ratioHistNameBG.c_str(), ratioHistNameBG.c_str(), N, xLow, xHigh);
    TH1D* yieldHist = new TH1D(yieldHistName.c_str(), yieldHistName.c_str(), N, xLow, xHigh);
    TH1D* ratio_w = new TH1D(ratioHistName_w.c_str(), ratioHistName_w.c_str(), 3*N, -1, 3);
    TH1D* ratio_w_bg = new TH1D(ratioHistNameBG_w.c_str(), ratioHistNameBG_w.c_str(), 3*N, -1, 3);


    cout << endl << endl << "Working on correlation histogram " << ratioHistName << ". first_bin last_bin first_N last_N : " << histBE[0] << " " << histBE[1] << " " << histBE[2] << " " << histBE[3] << endl;

    if(mc_slices){
      cout << "MC p. First_bin last_bin first_N last_N : " << histBE_mcp[0] << " " << histBE_mcp[1] << " " << histBE_mcp[2] << " " << histBE_mcp[3] << endl;
      
      cout << "MC n. First_bin last_bin first_N last_N : " << histBE_mcn[0] << " " << histBE_mcn[1] << " " << histBE_mcn[2] << " " << histBE_mcn[3] << endl;
    }

    //return;

    //data slices
    TH1D *cellslice[N];
    TH1D *cellslice_rd2[N];
    TH1D *cellslice_bgsub[N];
    TH1D *sliceclone[N];

    //mc slices
    TH1D *cellslice_mc[N];
    TH1D *cellslice_mc_rd2[N];
    TH1D *sliceclone_mc[N];

    double scalefactor[N][7];
    double scalefactor_bg[N][2];
    double sbpar[N][5];

    double sMin[N];
    double sMax[N];

    double scaleratio[N];
    double scaleratio_err[N];
    double scaleratio_bgsub[N];
    double scaleratio_bgsub_err[N];
    double totalevents[N];

    // Arrays for scale factors for hdx_p and hdx_n for each slice after pol4 bg subtract
    double pScale[N];
    double nScale[N];
    double scaleratio_bg[N];
    double scaleratio_bg_err[N];

    // store slice window for reporting
    double bin_low[N];
    double bin_high[N];

    // Make N slices and fit each
    for (int i = 0; i < N; ++i) {
      
      sMin[i] = 0;
      sMax[i] = 0;
      scaleratio[i] = 0.;
      scaleratio_err[i] = 0.;
      scaleratio_bgsub[i] = 0.;
      scaleratio_bgsub_err[i] = 0.;
      totalevents[i] = 0.;
      pScale[i] = 0.;
      nScale[i] = 0.;
      scaleratio_bg[i] = 0.;
      scaleratio_bg_err[i] = 0.;
      bin_low[i] = 0.;
      bin_high[i] = 0.;

      //initialize this slice array of scale factor ratios
      for (int sf=0; sf<7; ++sf)
	scalefactor[i][sf] = 0.;

      int binLow = first_bin + i * (last_bin - first_bin + 1) / N;
      int binHigh = first_bin + (i + 1) * (last_bin - first_bin + 1) / N - 1;

      if(fitrangeoverride){
	binLow = hist->GetXaxis()->FindBin(hcalfit_l);
	binHigh = hist->GetXaxis()->FindBin(hcalfit_h);
      }

      bin_low[i] = hist->GetXaxis()->GetBinCenter(binLow);
      bin_high[i] = hist->GetXaxis()->GetBinCenter(binHigh);

      //Do slices for fits to MC
      cellslice[i] = hist->ProjectionY(Form("%s_cellslice_%d", histName.c_str(), i+1), binLow, binHigh);
      cellslice_rd2[i] = hist->ProjectionY(Form("cellslice_rd2_%d", i+1), binLow, binHigh);
      
      if( cellslice[i] == nullptr || cellslice_mc[i] == nullptr )
	continue;

      //Get slices from corresponding MC hist
      if(mc_slices){
	hdx_slice_p = hist_mc_p->ProjectionY(Form("%s_slicemc_p_%d",histName.c_str(), i+1), binLow, binHigh);
	hdx_slice_n = hist_mc_n->ProjectionY(Form("%s_slicemc_n_%d",histName.c_str(), i+1), binLow, binHigh);
	if( cellslice[i]->GetEntries()<minEntries || 
	    hdx_slice_p->GetEntries()<minEntries || 
	    hdx_slice_n->GetEntries()<minEntries )
	  continue;
      
      }else{
	if( cellslice[i]->GetEntries()<minEntries )	
	  continue;
      }

      //Diagnostic for now - write out all histograms. N MUST BE SMALL
      cellslice[i]->Draw();
      cellslice[i]->Write();

      if(mc_slices){
      	hdx_slice_p->Draw();
      	hdx_slice_p->Write();
	
      	hdx_slice_n->Draw();
      	hdx_slice_n->Write();
      }

      cout << "....chugging through at: " << i << endl;

      int sliceNBins = cellslice[i]->GetNbinsX();
      double sliceMin = 0;
      double sliceMax = 0;

      int slice_nEntries = cellslice[i]->Integral(); // Total number of entries in the slice
      double slice_error = 1/sqrt(slice_nEntries); // Poisson error on histogram

      // Find first bin with data on TH1D slice
      for (int bin = 1; bin <= sliceNBins; ++bin) {
        if (cellslice[i]->GetBinContent(bin) != 0) {
	  sliceMin = cellslice[i]->GetXaxis()->GetBinCenter(bin);
	  break;
        }
      }
      sMin[i] = sliceMin;

      // Find last bin with data on TH1D slice
      for (int bin = sliceNBins; bin >= 1; --bin) {
        if (cellslice[i]->GetBinContent(bin) != 0) {
	  sliceMax = cellslice[i]->GetXaxis()->GetBinCenter(bin);
	  break;
        }
      }
      sMax[i] = sliceMax;

      //Get Scale Factor ratio from this slice
      TF1* fitFunc = new TF1("fitFunc", fitSliceShift, sliceMin, sliceMax, 9);
      fitFunc->SetNpx(5000);

      //reset fit parameters to eliminate hold-over effect
      for( int i=0; i<9; ++i )
	fitFunc->SetParameter(i,0.);
      
      cellslice[i]->Fit(fitFunc, fitopt.c_str());  // Fitting the slice

      double mcParams_array[9]={0.};
      for( int i=0; i<9; ++i )
	mcParams_array[i] = fitFunc->GetParameter(i);
      
      //Fit again for fine tuning
      TF1* fitFunc_rd2 = new TF1("fitFunc_rd2", fitSliceShift, sliceMin, sliceMax, 9);
      fitFunc_rd2->SetNpx(5000);
      fitFunc_rd2->SetParameters(mcParams_array);

      cellslice_rd2[i]->Fit(fitFunc_rd2, fitopt.c_str());  // Fitting the slice second with fine tuned params

      double mcParams_rd2_array[9]={0.};
      for( int i=0; i<9; ++i ){
	mcParams_rd2_array[i] = fitFunc_rd2->GetParameter(i);
      }

      // Get the error on the fit parameters
      int nParams = fitFunc->GetNpar(); // Gets the number of parameters
      double paramErrors[9]; // Array to store the errors

      for (int j = 0; j < nParams; ++j) 
	paramErrors[j] = fitFunc->GetParError(j); // Get the error for each parameter
      
      // Populate the scalefactor array
      for (int p = 0; p < 7; ++p) scalefactor[i][p] = mcParams_rd2_array[p];

      // Calculate the scale factor ratio and its error
      if (scalefactor[i][0] != 0) { // Check for division by zero for the denominator
	scaleratio[i] = scalefactor[i][1] / scalefactor[i][0];
	scaleratio_err[i] = sqrt(pow(paramErrors[1]/scalefactor[i][1], 2) + pow(paramErrors[0]/scalefactor[i][0], 2)) * scaleratio[i]; // quadrature error
      } else {
	scaleratio[i] = 0;
	scaleratio_err[i] = 0;
      }

      //So slices for sideband fits
      sliceclone[i] = hist->ProjectionY(Form("slice_clone_%d", i+1), binLow, binHigh);

      //BG subtract with sideband fit and get yield for this slice
      TF1 *sbfit = new TF1("sbfit",BGfit,sliceMin,sliceMax,5);
      sliceclone[i]->Fit("sbfit","RBMQ");

      const double* sbParams = sbfit->GetParameters();
      for (int p = 0; p < 5; ++p) sbpar[i][p] = sbParams[p];

      double *sbpars = sbfit->GetParameters();

      TF1 *sb = new TF1("sb",fits::g_p4fit,sliceMin,sliceMax,5);
      sb->SetParameters(&sbpars[0]);

      double sberr; util::subtractFunctionAndGetTotalAndError(sliceclone[i],
							      sb,
							      sb_min,
							      sb_max,
							      totalevents[i],
							      sliceMin,
							      sliceMax,
							      sberr);

      //Fill the TH1Ds with ratios and yields
      ratioHist->SetBinContent(i, scaleratio[i]);
      ratioHist->SetBinError(i, scaleratio_err[i]);
      yieldHist->SetBinContent(i, totalevents[i]);
      yieldHist->SetBinError(i, sberr);

      //cout << "For bin " << i << ", ratio " << scaleratio[i] << " and yield " << totalevents[i] << endl;

      //Add scale factor after BG subtract for comparison.

      // Create polynomial background function from scalefactor parameters
      TF1 *pol4 = new TF1(Form("pol4_%d", i), fits::g_p4fit, sMin[i], sMax[i], 5);
      pol4->SetParameters(&scalefactor[i][4]); // Assuming indices 2-6 are pol4 parameters

      // Background subtract cellslice
      cellslice_bgsub[i] = (TH1D*)cellslice[i]->Clone(Form("bgSubtractedSlice_%d", i));
      for (int bin = 1; bin <= cellslice_bgsub[i]->GetNbinsX(); ++bin) {
	double bgValue = pol4->Eval(cellslice_bgsub[i]->GetXaxis()->GetBinCenter(bin));
	cellslice_bgsub[i]->SetBinContent(bin, cellslice_bgsub[i]->GetBinContent(bin) - bgValue);
	cellslice_bgsub[i]->SetBinError(bin, cellslice_bgsub[i]->GetBinError(bin));
      }

      // Define fit function for BG subtracted histograms
      TF1 *mcFit = new TF1(Form("mcFit_%d", i), fitSliceShift_nobg, sliceMin, sliceMax,4);

      //reset fit parameters to eliminate hold-over effect
      for( int i=0; i<4; ++i )
	mcFit->SetParameter(i,0.);

      // Fit the background subtracted cellslice with the BG subtracted histograms
      cellslice_bgsub[i]->Fit(mcFit, fitopt.c_str());

      // Get the error on the fit parameters
      int nParams_bg = mcFit->GetNpar(); // Gets the number of parameters
      double paramErrors_bg[4]; // Array to store the errors

      for (int j = 0; j < nParams_bg; ++j) {
	paramErrors_bg[j] = mcFit->GetParError(j); // Get the error for each parameter
      }

      // Store the scale factors
      pScale[i] = mcFit->GetParameter(0);
      nScale[i] = mcFit->GetParameter(1);
      //scaleratio_bg[i] = nScale[i]/pScale[i];
      
      // Calculate the background subtracted scale factor ratio and its error
      if (pScale[i] != 0) { // Check for division by zero for the denominator
	scaleratio_bg[i] = nScale[i]/pScale[i];
	scaleratio_bg_err[i] = sqrt(pow(paramErrors_bg[1]/pScale[i], 2) + pow(paramErrors_bg[0]/nScale[i], 2)) * scaleratio_bg[i]; // Error propagation formula
      } else {
	scaleratio_bg[i] = 0;
	scaleratio_bg_err[i] = 0;
      }

      scalefactor_bg[i][0] = pScale[i];
      scalefactor_bg[i][1] = nScale[i];

      ratioHist_bg->SetBinContent(i, scaleratio_bg[i]);
      ratioHist_bg->SetBinError(i, scaleratio_bg_err[i]);

    } //endloop over N slices


    ////////////////////////////////////////////////////////
    cout << "Filling ratio histograms..." << endl;
    
    //Now get weighted histograms
    // Loop over all bins and print the value in each bin
    int ratioHist_nBins = ratioHist->GetXaxis()->GetNbins();
    for (int k = 1; k <= ratioHist_nBins; k++) { // bin indices start from 1
      double ratioHist_binValue = ratioHist->GetBinContent(k);
      double ratioHist_bg_binValue = ratioHist_bg->GetBinContent(k);
      double yieldHist_binValue = yieldHist->GetBinContent(k);
      
      //Should use error for weights. 1/err^2
      double weight = ratioHist->GetBinError(k);
      double weight_bg = ratioHist_bg->GetBinError(k);

      // ratio_w->Fill(ratioHist_binValue,yieldHist_binValue);
      // ratio_w_bg->Fill(ratioHist_bg_binValue,yieldHist_binValue);

      ratio_w->Fill(ratioHist_binValue,weight);
      ratio_w_bg->Fill(ratioHist_bg_binValue,weight_bg);
    }

    ////////////////////////////////////////////////////////
    cout << "Checking ratioHist for division by zero..." << endl;

    // Before writing the histogram, check for invalid values
    for (int i = 1; i <= ratioHist->GetNbinsX(); ++i) {
      if (std::isinf(ratioHist->GetBinContent(i)) || std::isnan(ratioHist->GetBinContent(i))) {
	std::cerr << "Invalid value found in histogram " << ratioHist->GetName() << " bin " << i << std::endl;
	ratioHist->SetBinContent(i, 0); // Setting to 0 or some default value
      }
    }

    ////////////////////////////////////////////////////////
    cout << "Writing the MC canvas with example fits..." << endl;

    //Now write out the mc canvas
    TCanvas* canvasSlices = new TCanvas(Form("MC Fits %s", hist->GetName()), hist->GetTitle(), 1200, 600);
    canvasSlices->Divide(3, 2); // Divide canvas into sub-pads for slices
    for( int i=0; i<6; ++i ){
      canvasSlices->cd(i + 1); // Go to the ith sub-pad

      //Get even sample
      int j = i*(int)(N/6);
      //int j = i+20;

      if( cellslice[j] == nullptr || j>N )
	continue;

      cellslice[j]->SetTitle(Form("slice window: %0.2f - %0.2f, ratio: %0.2f",bin_low[j],bin_high[j],scaleratio[j]));
      cellslice[j]->Draw();
      TF1 *mcfit;
      if(mc_slices)
	mcfit = new TF1("mcfit", fitSlices, sMin[i], sMax[i], 7);
      else
	mcfit = new TF1("mcfit", fitSliceShift, sMin[i], sMax[i], 9);
      mcfit->SetParameters(&scalefactor[j][0]);
      mcfit->SetLineColor(kGreen);
      mcfit->Draw("SAME");
      canvasSlices->Update();
    }	
    canvasSlices->Write();

    //continue;
    
    /*
    ////////////////////////////////////////////////////////
    cout << "Writing the BG sub MC canvas with example fits..." << endl;

    //Now write out the mc bg sub canvas
    TCanvas* canvasSlices_bg = new TCanvas(Form("MC Fits BG Sub %s", hist->GetName()), hist->GetTitle(), 800, 600);
    canvasSlices_bg->Divide(3, 2); // Divide canvas into sub-pads for slices
    for( int i=0; i<6; ++i ){
      canvasSlices_bg->cd(i + 1); // Go to the ith sub-pad

      //Get even sample
      int j = i*(int)(N/6);
      //int j = i+20;


      cout << i << " " << j << " " << sliceclone[j]->GetEntries() << endl;

      if(cellslice_bgsub[j] == nullptr || j>N)
	continue;

      cellslice_bgsub[j]->SetTitle(Form("slice window: %0.2f - %0.2f, ratio bg: %0.2f",bin_low[j],bin_high[j],scaleratio_bg[j]));
      cellslice_bgsub[j]->Draw();
      TF1 *mcfit = new TF1("mcfit", fitSlice_nobg, sMin[i], sMax[i], 7);
      mcfit->SetParameters(&scalefactor_bg[j][0]);
      mcfit->SetLineColor(kRed);
      mcfit->Draw("SAME");
      canvasSlices_bg->Update();
    }	
    canvasSlices_bg->Write();

    ////////////////////////////////////////////////////////
    cout << "Writing the sideband BG fit canvas with example fits..." << endl;

    //Now write out the sideband canvas
    TCanvas* canvasClones = new TCanvas(Form("Sideband Fits %s", hist->GetName()), hist->GetTitle(), 800, 600);
    canvasClones->Divide(3, 2); // Divide canvas into sub-pads for slices
    for( int i=0; i<6; ++i ){
      canvasClones->cd(i + 1); // Go to the ith sub-pad

      //Get even sample
      int j = i*(int)(N/6);
      //int j = i+20;

      if(sliceclone[j] == nullptr || j>N)
	continue;

      sliceclone[j]->SetTitle(Form("slice window: %0.2f - %0.2f, p+n yield: %0.2f",bin_low[j],bin_high[j],totalevents[j]));
      sliceclone[j]->Draw();
      TF1 *sbfit = new TF1("sbfit",fits::g_p4fit,sMin[i],sMax[i],5);
      sbfit->SetParameters(&sbpar[j][0]);
      sbfit->SetLineColor(kGreen);
      sbfit->Draw("SAME");
      canvasClones->Update();
    }	
    */

    //canvasClones->Write();

    ratioHist->Draw();
    ratioHist->Write();

    ratioHist_bg->Draw();
    ratioHist_bg->Write();

    yieldHist->Draw();
    yieldHist->Write();

    // Ratio histogram with RMS
    TCanvas *c1 = new TCanvas(Form("SF R %s", hist->GetName()), Form("SF R %s", hist->GetTitle()), 1200, 600);
    ratio_w->SetLineColor(kBlack);
    ratio_w->SetLineWidth(1);
    ratio_w->Draw("hist"); // Draw histogram as a line
    ratio_w->Draw("E same"); // Draw histogram errors on the same plot

    // Creating a legend for the first histogram
    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend->SetHeader("Full pol4 BG Fit", "C");
    char rmsText[255];
    double rms = ratio_w->GetRMS();
    double rmsError = ratio_w->GetRMSError();
    sprintf(rmsText, "RMS = %.2f #pm %.2f", rms, rmsError);
    legend->AddEntry((TObject*)0, rmsText, "");
    legend->Draw();
    c1->Update(); // Update the canvas to reflect changes

    // Ratio histogram with RMS (Background Subtracted Version)
    TCanvas *c2 = new TCanvas(Form("SF R BG sub %s", hist->GetName()), Form("SF R BG sub %s", hist->GetTitle()), 1200, 600);
    ratio_w_bg->SetLineColor(kBlack);
    ratio_w_bg->SetLineWidth(1);
    ratio_w_bg->Draw("hist"); 
    ratio_w_bg->Draw("E same");

    // Creating a legend for the second histogram
    TLegend *legend_bg = new TLegend(0.6, 0.7, 0.9, 0.9);
    legend_bg->SetHeader("pol4 BG Subtracted", "C");
    char rmsText_bg[255];
    double rms_bg = ratio_w_bg->GetRMS();
    double rmsError_bg = ratio_w_bg->GetRMSError();
    sprintf(rmsText_bg, "RMS = %.2f #pm %.2f", rms_bg, rmsError_bg);
    legend_bg->AddEntry((TObject*)0, rmsText_bg, "");
    legend_bg->Draw();
    c2->Update(); // Update the canvas to reflect changes

    // Optionally, write the canvases to a file
    c1->Write();
    c2->Write();

  }//endloop over TH2D histograms

  outputFile->Write();
  outputFile->Close();

  cout << "All plots created. Output file located here: " << foutPath << endl;

}

void handleError(TFile *file1, TFile *file2, std::string marker) {
    if (file1) file1->Close();
    if (file2) file2->Close();
    std::cerr << "Error: File opening or histogram retrieval failed at " << marker << "." << std::endl;
}

void handleError(TFile *file1, std::string marker) {
    if (file1) file1->Close();
    std::cerr << "Error: File opening or histogram retrieval failed at " << marker << "." << std::endl;
}

void createAndSaveCanvasWithFitsAndResiduals(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, TF1* bg, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, bool shiftfit) {

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

  // Construct legend within function if preferred
  TLegend* leg = new TLegend(0.6, 0.5, 0.89, 0.89);
  leg->AddEntry(hdx, "Data", "l");
  leg->AddEntry(fit, "Fit", "l");
  leg->AddEntry(bg, "Background", "l");
  leg->AddEntry(hdx_p_clone, "Proton SIMC MC", "l");
  leg->AddEntry(hdx_n_clone, "Neutron SIMC MC", "l");
  leg->AddEntry( (TObject*)0, "", "");
  leg->AddEntry( (TObject*)0, Form("data N events : %0.0f",hdx->GetEntries()), "");
  leg->AddEntry( (TObject*)0, Form("MC N events : %0.0f",(hdx_p_clone->GetEntries()+hdx_n_clone->GetEntries())), "");
  leg->AddEntry( (TObject*)0, Form("n/p scale ratio R_{sf} : %0.3f #pm %0.3f (fit error)",np_par_ratio, np_par_ratio_error), "");
  leg->AddEntry( (TObject*)0, Form("n/p yield ratio R'' : %0.3f #pm %0.3f (fit and sum statistical error)",np_sum_ratio, np_sum_ratio_error), "");
  if(shiftfit)
    leg->AddEntry( (TObject*)0, Form("dx shift pars, n/p : %0.3f / %0.3f ", pars[2].first, pars[3].first), "");
  leg->AddEntry( (TObject*)0, Form("#chi^{2}/ndf: %0.3f",fit->GetChisquare()/fit->GetNDF()), "");
  leg->Draw("same");

  // Construct residuals
  //add the p and n histograms together for comparison
  TH1D* sumHistogram = new TH1D("sumHistogram", "Sum of Histograms", hdx_p_clone->GetNbinsX(), hdx_p_clone->GetXaxis()->GetXmin(), hdx_p_clone->GetXaxis()->GetXmax());
  sumHistogram->Add(hdx_p_clone, hdx_n_clone);
  //sumHistogram->SetName(Form("sh_%s",savePath)); //stop the warning messages

  for (int bin = 1; bin <= hdx->GetNbinsX(); ++bin) {
    double bgValue = bg->Eval(hdx_clone->GetXaxis()->GetBinCenter(bin));
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
void createAndSaveCanvasWithFitsAndResiduals_nobg(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, bool shiftfit) {

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

  // Construct legend within function if preferred
  TLegend* leg = new TLegend(0.6, 0.5, 0.89, 0.89);
  leg->AddEntry(hdx, "Data", "l");
  leg->AddEntry(fit, "Fit", "l");
  leg->AddEntry(hdx_p_clone, "Proton SIMC MC", "l");
  leg->AddEntry(hdx_n_clone, "Neutron SIMC MC", "l");
  leg->AddEntry( (TObject*)0, "", "");
  leg->AddEntry( (TObject*)0, Form("data N events : %0.0f",hdx->GetEntries()), "");
  leg->AddEntry( (TObject*)0, Form("MC N events : %0.0f",(hdx_p_clone->GetEntries()+hdx_n_clone->GetEntries())), "");
  leg->AddEntry( (TObject*)0, Form("n/p scale ratio R_{sf} : %0.3f #pm %0.3f (fit error)",np_par_ratio, np_par_ratio_error), "");
  leg->AddEntry( (TObject*)0, Form("n/p yield ratio R'' : %0.3f #pm %0.3f (fit and sum statistical error)",np_sum_ratio, np_sum_ratio_error), "");
  if(shiftfit)
    leg->AddEntry( (TObject*)0, Form("dx shift pars, n/p : %0.3f / %0.3f ", pars[2].first, pars[3].first), "");
  leg->AddEntry( (TObject*)0, Form("#chi^{2}/ndf: %0.3f",fit->GetChisquare()/fit->GetNDF()), "");
  leg->Draw("same");

  // Construct residuals
  //add the p and n histograms together for comparison
  TH1D* sumHistogram = new TH1D("sumHistogram", "Sum of Histograms", hdx_p_clone->GetNbinsX(), hdx_p_clone->GetXaxis()->GetXmin(), hdx_p_clone->GetXaxis()->GetXmax());
  sumHistogram->Add(hdx_p_clone, hdx_n_clone);
  //sumHistogram->SetName(Form("nobg_sh_%s",savePath)); //stop the warning messages

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
void createAndSaveCanvasWithFitsAndResiduals_altbg(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, TH1D* hdxaltbg, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, std::string bgtype) {

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

  // Construct legend within function if preferred
  TLegend* leg = new TLegend(0.6, 0.5, 0.89, 0.89);
  leg->AddEntry( hdx, "Data", "l" );
  leg->AddEntry( fit, "Fit", "l" );
  leg->AddEntry( hdx_p_clone, "Proton SIMC MC", "l" );
  leg->AddEntry( hdx_n_clone, "Neutron SIMC MC", "l" );
  leg->AddEntry( hdx_altbg_clone, Form("BG, %s",bgtype.c_str()), "l" );
  leg->AddEntry( (TObject*)0, "", "");
  leg->AddEntry( (TObject*)0, Form("data N events : %0.0f",hdx->GetEntries()), "");
  leg->AddEntry( (TObject*)0, Form("MC N events : %0.0f",(hdx_p_clone->GetEntries()+hdx_n_clone->GetEntries())), "");
  leg->AddEntry( (TObject*)0, Form("n/p scale ratio R_{sf} : %0.3f #pm %0.3f (fit error)",np_par_ratio, np_par_ratio_error), "");
  leg->AddEntry( (TObject*)0, Form("n/p yield ratio R'' : %0.3f #pm %0.3f (fit and sum statistical error)",np_sum_ratio, np_sum_ratio_error), "");
  leg->AddEntry( (TObject*)0, Form("dx shift pars, n/p : %0.3f / %0.3f ", pars[2].first, pars[3].first), "");
  leg->AddEntry( (TObject*)0, Form("bg scale: %0.3f", pars[4].first), "");
  leg->AddEntry( (TObject*)0, Form("#chi^{2}/ndf: %0.3f",fit->GetChisquare()/fit->GetNDF()), "");
  leg->Draw("same");

  // Construct residuals
  //add the p and n histograms together for comparison
  TH1D* sumHistogram = new TH1D("sumHistogram", "Sum of Histograms", hdx_p_clone->GetNbinsX(), hdx_p_clone->GetXaxis()->GetXmin(), hdx_p_clone->GetXaxis()->GetXmax());
  sumHistogram->Add(hdx_p_clone, hdx_n_clone);
  //sumHistogram->SetName(Form("nobg_sh_%s",savePath)); //stop the warning messages

  //subtract the background
  for (int bin = 1; bin <= hdx_clone->GetNbinsX(); ++bin) {
    // Get the value at the current bin in bg_hist
    double bgValue = hdx_altbg_clone->GetBinContent(bin);
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
