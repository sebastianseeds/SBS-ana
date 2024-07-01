//sseeds 4.22.24: script to vary course fiducial cut and extract Rsf. Comparison on dependent ranges determines fine cut.

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

//fitranges
//sbs9 70p: -1.8 to 0.7
//sbs8 70p: -1.8 to 0.7
//sbs4 30p: -1.7 to 0.7, shiftX=0.0, neutronshift=-0.05
//sbs4 50p: -2.1 to 0.7, shiftX=0.05, neutronshift=0.0
//sbs7 85p: -1.4 to 0.5

//Fit range override options
double hcalfit_l = -2.1; //lower fit/bin limit for hcal dx plots (m) sbs4, 50p
double hcalfit_h = 0.7; //upper fit/bin limit for hcal dx plots (m) sbs4, 50p

//Plot range option shared by all fits
double hcalr_l = -2.1; //lower fit/bin limit for hcal dx plots (m) sbs4, 50p
double hcalr_h = 0.7; //upper fit/bin limit for hcal dx plots (m) sbs4, 50p

double xrange_min_fidx = 0.;
double xrange_max_fidx = 23.;
double xrange_min_fidy = 0.;
double xrange_max_fidy = 13.;

//Fit options
std::string fitopt = "RMQ0";
//std::string fitopt = "RBMQ0"; //R:setrange,B:predefinedTF1,M:improvefit,Q:quiet,0:nofitline
//std::string fitopt = "RLQ0"; //L:loglikelihood for counts
//std::string fitopt = "RLEMQ0"; //E:bettererrorest
//std::string fitopt = "RWLEMQ0"; //WL:weightedloglikelihood for weighted counts (N/A here)
//std::string fitopt = "RPEMQ0"; //P:pearsonloglikelihood for expected error (N/A here)

//fit min entries
int minEvents = 500;

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
  return proton + neutron + fits::g_p2fit_cd(x, &par[4]);
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

// Forward declarations
void handleError(TFile *file1, TFile *file2, std::string marker);
void handleError(TFile *file1, std::string marker);
std::vector<std::pair<double, double>> fitAndFineFit(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, double pshift, double nshift, const std::string& fitOptions = "RBMQ0");

//main. kine=kinematic, mag=fieldsetting, pass=pass#, N=total range plots, sliceCutMin=minimum cut on slices, sliceCutMax=max cut on slices, BG=background function, bestclus=use plot file with best cluster selection, thin=use plot file without correlations, wide=use plots with wide cuts, effz=use plots with effective z applied, alt=use plot file using MC alt files
void fiducialStability(int kine=8, 
		       int mag=100, 
		       int pass=2, 
		       int N=12,
		       double sliceCutMin=0.05,
		       double sliceCutMax=5.0,
		       std::string BG="pol2",
		       bool dxopt=true,
		       bool bestclus=true, 
		       bool thin=true,
		       bool wide=false,
		       bool effz=true,
		       bool alt = true) {
  //set draw params
  gStyle->SetPalette(55);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);

  //load json file
  JSONManager *jmgr = new JSONManager("../../config/syst.json");

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
  if(wide)
    wide_word = "_widecut";

  std::string effz_word = "";
  if(effz)
    effz_word = "_effz";

  std::string fidy_word = "";
  if(!dxopt)
    fidy_word = "_fidy";

  std::string basePath = "/lustre19/expphy/volatile/halla/sbs/seeds";
  std::string finPath = Form("%s/gmn_analysis/dx_correlations%s_sbs%d_mag%d_pass%d%s%s%s.root", basePath.c_str(), bestclus_word.c_str(), kine, mag, pass, thin_word.c_str(), wide_word.c_str(), effz_word.c_str());
  std::string fmcinPath = Form("%s/gmn_analysis/dx_mc_sbs%d_mag%d_pass%d%s%s%s%s.root", basePath.c_str(), kine, mag, pass, thin_word.c_str(), alt_word.c_str(), wide_word.c_str(), effz_word.c_str());
  std::string foutPath = Form("%s/gmn_analysis/fidstab_gmn_sbs%d_mag%d_pass%d%s%s%s%s.root", basePath.c_str(), kine, mag, pass, bestclus_word.c_str(), wide_word.c_str(), effz_word.c_str(), fidy_word.c_str());

  TFile* inputFile = new TFile(finPath.c_str());
  if (!inputFile || inputFile->IsZombie()) 
    handleError(inputFile,"inputFile");

  //Get all histograms
  TH1D *hdx_data = dynamic_cast<TH1D*>(inputFile->Get("hdx_allcut"));
  hdx_data->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  TH2D *hdx_vs_fidx_data;
  if(dxopt)
    hdx_vs_fidx_data = dynamic_cast<TH2D*>(inputFile->Get("hist_fiducial_sig_x"));
  else
    hdx_vs_fidx_data = dynamic_cast<TH2D*>(inputFile->Get("hist_fiducial_sig_y"));

  //TH2D *hdx_vs_fidx_data = dynamic_cast<TH2D*>(inputFile->Get("hist_fiducial_sig_x"));
  hdx_vs_fidx_data->GetYaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  int hdx_vs_fidx_data_N = hdx_vs_fidx_data->GetEntries();
  std::string fidxCuts = hdx_vs_fidx_data->GetTitle();
  cout << endl << "Opened dx vs fiducial with cuts: " << fidxCuts << endl << endl;

  TH2D *hdx_vs_fidy_data = dynamic_cast<TH2D*>(inputFile->Get("hist_fiducial_sig_y"));
  hdx_vs_fidy_data->GetYaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  int hdx_vs_fidy_data_N = hdx_vs_fidy_data->GetEntries();
  cout << endl << "Opened dx vs fid y with cuts: " << hdx_vs_fidx_data->GetTitle() << endl << endl;  

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

  TH2D *hdx_vs_fidx_n;
  if(dxopt)
    hdx_vs_fidx_n = dynamic_cast<TH2D*>(inputFileMC->Get("hist_fiducial_sig_x_n"));
  else
    hdx_vs_fidx_n = dynamic_cast<TH2D*>(inputFileMC->Get("hist_fiducial_sig_y_n"));
  
  //TH2D *hdx_vs_fidx_n = dynamic_cast<TH2D*>(inputFileMC->Get("hist_fiducial_sig_x_n"));
  hdx_vs_fidx_n->GetYaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  int hdx_vs_fidx_n_N = hdx_vs_fidx_n->GetEntries();

  TH2D *hdx_vs_fidy_n = dynamic_cast<TH2D*>(inputFileMC->Get("hist_fiducial_sig_y_n"));
  hdx_vs_fidy_n->GetYaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  int hdx_vs_fidy_n_N = hdx_vs_fidy_n->GetEntries();

  TH2D *hdx_vs_fidx_p;
  if(dxopt)
    hdx_vs_fidx_p = dynamic_cast<TH2D*>(inputFileMC->Get("hist_fiducial_sig_x_p"));
  else
    hdx_vs_fidx_p = dynamic_cast<TH2D*>(inputFileMC->Get("hist_fiducial_sig_y_p"));

  //TH2D *hdx_vs_fidx_p = dynamic_cast<TH2D*>(inputFileMC->Get("hist_fiducial_sig_x_p"));
  hdx_vs_fidx_p->GetYaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  int hdx_vs_fidx_p_N = hdx_vs_fidx_p->GetEntries();

  TH2D *hdx_vs_fidy_p = dynamic_cast<TH2D*>(inputFileMC->Get("hist_fiducial_sig_y_p"));
  hdx_vs_fidy_p->GetYaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  int hdx_vs_fidy_p_N = hdx_vs_fidy_p->GetEntries();

  hdx_inel = dynamic_cast<TH1D*>(inputFileMC->Get("hdx_inel"));
  hdx_inel->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  hdx_inel->Scale(10e33); //account for lack of overall normalization from g4sbs generator
  if (!hdx_p || !hdx_n || !hdx_inel) 
    handleError(inputFile,inputFileMC,"hdx");

  TFile* outputFile = new TFile(foutPath.c_str(), "RECREATE");

  TCanvas *canvas = new TCanvas("canvas", "Data and MC Histograms", 1800, 1200);
  canvas->Divide(4, 3);  // Adjust the grid size according to the number of histograms

  // Pad 1: hdx_data
  canvas->cd(1);
  hdx_data->Draw();

  // Pad 2: hdx_vs_fidx_data
  canvas->cd(2);
  hdx_vs_fidx_data->Draw("COLZ");

  // Pad 3: hdx_vs_fidy_data
  canvas->cd(3);
  hdx_vs_fidy_data->Draw("COLZ");

  // Pad 4: hdx_dyanti
  canvas->cd(4);
  hdx_dyanti->Draw();

  // Pad 5: hdx_coinanti
  canvas->cd(5);
  hdx_coinanti->Draw();

  // Pad 6: hdx_p
  canvas->cd(6);
  hdx_p->Draw();

  // Pad 7: hdx_n
  canvas->cd(7);
  hdx_n->Draw();

  // Pad 8: hdx_vs_fidx_n
  canvas->cd(8);
  hdx_vs_fidx_n->Draw("COLZ");

  // Pad 9: hdx_vs_fidy_n
  canvas->cd(9);
  hdx_vs_fidy_n->Draw("COLZ");

  // Pad 10: hdx_vs_fidx_p
  canvas->cd(10);
  hdx_vs_fidx_p->Draw("COLZ");

  // Pad 11: hdx_vs_fidy_p
  canvas->cd(11);
  hdx_vs_fidy_p->Draw("COLZ");

  // Pad 12: hdx_inel
  canvas->cd(12);
  hdx_inel->Draw();

  // Update the canvas to reflect the drawings
  canvas->Update();

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

  std::vector<ReportData> fidxreports;
  fidxreports.resize(N);
  // std::vector<ReportData> fidyreports;
  // fidyreports.resize(N);

  // Loop over dx vs fiducial x
  hdx_vs_fidx_data->GetXaxis()->SetRangeUser(xrange_min_fidx, xrange_max_fidx);
  hdx_vs_fidx_data->GetYaxis()->SetRangeUser(hcalfit_l, hcalfit_h);

  // Assuming current_llim and current_ulim are set to the desired x-axis range
  int xbinLow_fidx = hdx_vs_fidx_data->GetXaxis()->FindBin(xrange_min_fidx);
  int xbinWideLow_fidx = hdx_vs_fidx_data->GetXaxis()->FindBin(0.5); //Set beneath the cut to map the region
  int xbinHigh_fidx = hdx_vs_fidx_data->GetXaxis()->FindBin(xrange_max_fidx);
  int xbinCutMin_fidx = hdx_vs_fidx_data->GetXaxis()->FindBin(sliceCutMin);
  int xbinCutMax_fidx = hdx_vs_fidx_data->GetXaxis()->FindBin(sliceCutMax);
  int xbinBegin;

  //cout << sliceCutMax << " " << xbinCutMax_fidx << " " << xrange_max_fidx << " " << xbinHigh_fidx << endl;

  // Bins per range
  int xbinRanges_fidx;
  if(sliceCutMax>0 && sliceCutMin==0){
    int modrange = xbinCutMax_fidx - xbinLow_fidx;
    xbinBegin = xbinLow_fidx;
    if(modrange<N)
      N=modrange;
    xbinRanges_fidx = modrange / N;
  }else if(sliceCutMax==0 && sliceCutMin>0){
    int modrange = xbinHigh_fidx - xbinCutMin_fidx;
    xbinBegin = xbinCutMin_fidx;
    if(modrange<N)
      N=modrange;
    xbinRanges_fidx = modrange / N;
  }else if(sliceCutMax>0 && sliceCutMin>0){
    int modrange = xbinCutMax_fidx - xbinCutMin_fidx;
    xbinBegin = xbinCutMin_fidx;
    if(modrange<N)
      N=modrange;
    xbinRanges_fidx = modrange / N;
  }else{
    int modrange = xbinHigh_fidx - xbinLow_fidx;
    xbinBegin = xbinLow_fidx;
    if(modrange<N)
      N=modrange;
    xbinRanges_fidx = modrange / N;
  }

  cout << "Bins in x per cut: " << xbinRanges_fidx << " on N = " << N << " cuts." << endl;

  // Project the TH2D onto a TH1D for the specified x-axis range
  TH1D* fullProjY_fidx = hdx_vs_fidx_data->ProjectionY("_py", xbinLow_fidx, xbinHigh_fidx);
  TH1D* wideProjY_fidx = hdx_vs_fidx_data->ProjectionY("_wpy", xbinWideLow_fidx, xbinHigh_fidx);

  // Use Integral() to get the total number of events within the specified range
  int totalEvents_fidx = fullProjY_fidx->Integral();

  std::pair<double,double> wideQualp2;

  auto widep2par_vector = fitAndFineFit(wideProjY_fidx, "sliceFitWide", "fitFullShift_p2", 7, hcalfit_l, hcalfit_h, wideQualp2, 0, 0, fitopt.c_str());

  double fixed_pshift = widep2par_vector[2].first;
  double fixed_nshift = widep2par_vector[3].first;

  cout << "Fit over wide range of fiducial x (0.5-x_max) yield shifts p:n -> " << fixed_pshift << ":" << fixed_nshift << endl;

  vector<double> cutx_fidx;

  cout << "Looping over fiducial y slices.." << endl;
  for (int i = 0; i < N; ++i) {

    int binStart = xbinBegin + i*xbinRanges_fidx;
    int binEnd = xbinHigh_fidx;

    cutx_fidx.push_back(hdx_vs_fidx_data->GetXaxis()->GetBinLowEdge(binStart));
    //cout << "Fiducial cut at x=" << cutx_fidx[i] << endl;

    //Get the data slices
    TH1D *hdx_slice = hdx_vs_fidx_data->ProjectionY(Form("hdx_fidx_slice_%d", i+1), binStart, binEnd);
    TH1D *dx_fidx_slice = hdx_vs_fidx_data->ProjectionY(Form("dxfidcut_%d", i+1), binStart, binEnd);

    //Get the MC slices
    TH1D *dx_fidx_p_slice = hdx_vs_fidx_p->ProjectionY(Form("dxfidcut_p_%d", i+1), binStart, binEnd);
    TH1D *dx_fidx_n_slice = hdx_vs_fidx_n->ProjectionY(Form("dxfidcut_n_%d", i+1), binStart, binEnd);
    hdx_p = hdx_vs_fidx_p->ProjectionY(Form("hdx_fidx_slice_p_%d", i+1), binStart, binEnd);
    hdx_n = hdx_vs_fidx_n->ProjectionY(Form("hdx_fidx_slice_n_%d", i+1), binStart, binEnd);

    //Get total entries on slice
    int slice_nEntries = dx_fidx_slice->Integral();
    int slicemcp_nEntries = dx_fidx_p_slice->Integral();
    int slicemcn_nEntries = dx_fidx_n_slice->Integral();

    //Rsf extraction using second order poly fit
    std::pair<double,double> sliceQualp2;

    auto slicep2par_vector = fitAndFineFit(hdx_slice, "sliceFit", "fitFullShift_p2", 7, hcalfit_l, hcalfit_h, sliceQualp2, fixed_pshift, fixed_nshift, fitopt.c_str());

    double csndf = sliceQualp2.first/sliceQualp2.second;

    double scaleratio = slicep2par_vector[1].first / slicep2par_vector[0].first;
    double scaleratio_err = sqrt(pow(slicep2par_vector[1].second/slicep2par_vector[1].first, 2) + pow(slicep2par_vector[0].second/slicep2par_vector[0].first, 2)) * scaleratio;

    double xshift_p = slicep2par_vector[2].first;
    double xshift_n = slicep2par_vector[3].first;

    //Write fit results to struct
    fidxreports[i] = ReportData( dx_fidx_slice, 
				 dx_fidx_p_slice,
				 dx_fidx_n_slice,
				 cutx_fidx[i],
				 xrange_max_fidx,
				 slice_nEntries,
				 slicemcp_nEntries,
				 slicemcn_nEntries,
				 scaleratio,
				 scaleratio_err,
				 csndf,
				 "dx vs fiducial x" );

    //cout << endl << endl << "Writing out all dx vs fid fit pars: " << endl;
    for (int par=0; par<7; ++par){
      fidxreports[i].fitParams[par] = slicep2par_vector[par].first;
      fidxreports[i].fitErrors[par] = slicep2par_vector[par].second;

      //cout << "p" << par << " " << dxreports[i].fitParams[par] << " ";
      
    }
    
    //cout << endl << endl;

  }
  /*
  
  // Loop over dx vs fiducial y
  hdx_vs_fidy_data->GetXaxis()->SetRangeUser(xrange_min_fidy, xrange_max_fidy);
  hdx_vs_fidy_data->GetYaxis()->SetRangeUser(hcalfit_l, hcalfit_h);
  
  // Assuming current_llim and current_ulim are set to the desired x-axis range
  int xbinLow_fidy = hdx_vs_fidy_data->GetXaxis()->FindBin(xrange_min_fidy);
  int xbinHigh_fidy = hdx_vs_fidy_data->GetXaxis()->FindBin(xrange_max_fidy);
  
  // Bins per range
  int xbinRanges_fidy = (xbinHigh_fidy - xbinLow_fidy) / N;
  
  // Project the TH2D onto a TH1D for the specified x-axis range
  TH1D* wideProjY_fidy = hdx_vs_fidy_data->ProjectionY("_py", xbinLow_fidy, xbinHigh_fidy);
  
  // Use Integral() to get the total number of events within the specified range
  int totalEvents_fidy = wideProjY_fidy->Integral();
  vector<double> cutx_fidy;

  cout << "Looping over fiducial y slices.." << endl;
  for (int i = 0; i < N; ++i) {
    
    int binStart = xbinLow_fidy + i*xbinRanges_fidy;
    int binEnd = xbinHigh_fidy;
    
    cutx_fidy.push_back(hdx_vs_fidy_data->GetXaxis()->GetBinLowEdge(binStart));
    cout << "Fiducial cut at x=" << cutx_fidy[i] << endl;
    
    //Get the data slices
    TH1D *hdx_slice = hdx_vs_fidy_data->ProjectionY(Form("hdx_fidy_slice_%d", i+1), binStart, binEnd);
    TH1D *dx_fidy_slice = hdx_vs_fidy_data->ProjectionY(Form("dxfidcut_%d", i+1), binStart, binEnd);
    
    //Get the MC slices
    TH1D *dx_fidy_p_slice = hdx_vs_fidy_p->ProjectionY(Form("dxfidcut_p_%d", i+1), binStart, binEnd);
    TH1D *dx_fidy_n_slice = hdx_vs_fidy_n->ProjectionY(Form("dxfidcut_n_%d", i+1), binStart, binEnd);
    hdx_p = hdx_vs_fidy_p->ProjectionY(Form("hdx_fidy_slice_p_%d", i+1), binStart, binEnd);
    hdx_n = hdx_vs_fidy_n->ProjectionY(Form("hdx_fidy_slice_n_%d", i+1), binStart, binEnd);
    
    //Get total entries on slice
    int slice_nEntries = dx_fidy_slice->Integral(binStart,binEnd);
    int slicemcp_nEntries = dx_fidy_p_slice->Integral(binStart,binEnd);
    int slicemcn_nEntries = dx_fidy_n_slice->Integral(binStart,binEnd);
    
    //Rsf extraction using second order poly fit
    std::pair<double,double> sliceQualp2;
    
    auto slicep2par_vector = fitAndFineFit(hdx_slice, "sliceFit", "fitFullShift_p2", 7, hcalfit_l, hcalfit_h, sliceQualp2, fitopt.c_str());
    
    double csndf = sliceQualp2.first/sliceQualp2.second;
    
    double scaleratio = slicep2par_vector[1].first / slicep2par_vector[0].first;
    double scaleratio_err = sqrt(pow(slicep2par_vector[1].second/slicep2par_vector[1].first, 2) + pow(slicep2par_vector[0].second/slicep2par_vector[0].first, 2)) * scaleratio;
    
    double xshift_p = slicep2par_vector[2].first;
    double xshift_n = slicep2par_vector[3].first;
    
    //Write fit results to struct
    fidyreports[i] = ReportData( dx_fidy_slice, 
				 dx_fidy_p_slice,
				 dx_fidy_n_slice,
				 cutx_fidy[i],
				 xrange_max_fidy,
				 slice_nEntries,
				 slicemcp_nEntries,
				 slicemcn_nEntries,
				 scaleratio,
				 scaleratio_err,
				 csndf,
				 "dx vs fiducial y" );
    
    for (int par=0; par<7; ++par){
      fidyreports[i].fitParams[par] = slicep2par_vector[par].first;
      fidyreports[i].fitErrors[par] = slicep2par_vector[par].second;
    }
  
 
  }
  
  */

  //Write out the dx histograms
  
  //get a clone for displays
  TH2D* clonedFidx = (TH2D*)hdx_vs_fidx_data->Clone("clonedFidx");
  // double fidx_ymin = clonedFidx->GetYaxis()->GetXmin();
  // double fidx_ymax = clonedFidx->GetYaxis()->GetXmax();
  double fidx_ymin = hcalfit_l;
  double fidx_ymax = hcalfit_h;

  TCanvas* c0 = new TCanvas("c0", "Slice Locations", 800,600);
  c0->cd();
  clonedFidx->Draw("colz");

  for (auto& cut : cutx_fidx){
    TLine* line = new TLine(cut, fidx_ymin, cut, fidx_ymax);
    line->SetLineWidth(2);
    line->SetLineColor(kRed); // Set line color to red
    line->Draw();
  }
  c0->Update();

  //set up vectors for later tgrapherrors
  std::vector<double> Rsf_vec_fidx;
  std::vector<double> Rsferr_vec_fidx;
  std::vector<double> xval_vec_fidx;
  std::vector<double> nev_vec_fidx;

  TCanvas* canvasSlices_fidx = new TCanvas("dx slices over fiducial x", "dx slices over fiducial x", 1800, 1200);
  
  // Assuming N is the total number of histograms to be plotted on the canvas
  int optimalRows_fidx = 1;
  int optimalCols_fidx = 1;
  
  // Find the nearest square root for N to determine the grid size
  int sqrtN_fidx = (int)std::sqrt(N);
  for (int cols = sqrtN_fidx; cols <= N; ++cols) {
    if (N % cols == 0) { // If cols is a factor of N
      optimalCols_fidx = cols;
      optimalRows_fidx = N / cols;
      break; // Found the optimal layout
    }
  }
  
  canvasSlices_fidx->Divide(optimalCols_fidx, optimalRows_fidx); // Adjust the division based on reportSamples or desired layout
  
  TF1* bgfit_fidx[N];
  
  // loop over slices
  for (int i = 0; i < N; ++i) {
    if (fidxreports[i].sliceHistogram == nullptr)
      continue;
    
    canvasSlices_fidx->cd(i + 1);
    
    //catch bad fit ranges and do not write
    std::string tempTitle = fidxreports[i].sliceHistogram->GetTitle();      
    if( tempTitle.compare("Null Histogram")==0 )
      continue;
    
    //get mc nucleon parameters
    double p_scale = fidxreports[i].fitParams[0];
    double n_scale = fidxreports[i].fitParams[1];
    double p_shift = fidxreports[i].fitParams[2];
    double n_shift = fidxreports[i].fitParams[3];
    
    double MCev = fidxreports[i].mcpnev + fidxreports[i].mcnnev;
    
    // Simplify the histogram title
    fidxreports[i].sliceHistogram->SetTitle(Form("%s %0.2f to %0.2f", fidxreports[0].type.c_str(), fidxreports[i].winLow, fidxreports[i].winHigh));
    fidxreports[i].sliceHistogram->SetLineWidth(1);
    fidxreports[i].sliceHistogram->Draw();
    
    // Retrieve the number of bins, and the x-axis limits of the slice histogram
    int nbins = fidxreports[i].sliceHistogram->GetNbinsX();
    double x_low = fidxreports[i].sliceHistogram->GetXaxis()->GetXmin();
    double x_high = fidxreports[i].sliceHistogram->GetXaxis()->GetXmax();
    
    //get scale ratio, error, and midpoint for slice
    double ratio = fidxreports[i].scaleratio;
    double error = fidxreports[i].scaleratioerr;
    double xval = fidxreports[i].winLow;
    double nev = (totalEvents_fidx - fidxreports[i].nev)/totalEvents_fidx*100;
    
    //fill vectors for later tgrapherrors
    Rsf_vec_fidx.push_back(ratio);
    Rsferr_vec_fidx.push_back(error);
    xval_vec_fidx.push_back(xval);
    nev_vec_fidx.push_back(nev); //Total remaining stats

    // Draw the corresponding MC histograms
    hdx_p = fidxreports[i].slicePHistogram; //update the MC slice for the later overall fit
    TH1D *pMCslice = util::shiftHistogramX(fidxreports[i].slicePHistogram, p_shift);
    pMCslice->Scale(p_scale);
    
    hdx_n = fidxreports[i].sliceNHistogram; //update the MC slice for the later overall fit
    TH1D *nMCslice = util::shiftHistogramX(fidxreports[i].sliceNHistogram, n_shift);
    nMCslice->Scale(n_scale);
    
    pMCslice->Draw("same");
    nMCslice->Draw("same");

    //cout << "bin " << i << " pshift " << p_shift << " nshift " << n_shift << endl;
    
    // Create a legend and add the remaining parameters
    TLegend* legend = new TLegend(0.49, 0.59, 0.89, 0.89); // Adjust the position as needed
    legend->SetMargin(0.1);
    legend->AddEntry((TObject*)0, Form("Nev: %0.0f", fidxreports[i].nev), "");
    legend->AddEntry((TObject*)0, Form("MCev: %0.0f", MCev), "");
    legend->AddEntry((TObject*)0, Form("Ratio: %0.2f", ratio), "");
    legend->AddEntry((TObject*)0, Form("Shift p/n: %0.2f/%0.2f", p_shift, n_shift), "");
    legend->AddEntry((TObject*)0, Form("#chi^{2}/ndf: %0.2f", fidxreports[i].chisqrndf), "");
    legend->Draw();
    
    bgfit_fidx[i] = new TF1(Form("bgfit_fidx_%d",i),fits::g_p2fit_cd,hcalfit_l,hcalfit_h,3);
    //cout << "Background fit parameter" << endl;
    for (int j=0; j<3; ++j){
      bgfit_fidx[i]->SetParameter(j,fidxreports[i].fitParams[j+4]);
      //cout << "   " << j << " = " << dxreports[i].fitParams[j+4] << endl;
    }
    
    bgfit_fidx[i]->SetLineWidth(2);
    bgfit_fidx[i]->Draw("same");
    
    TF1 *mcfit = new TF1(Form("mcfit_%d", i), fitFullShift_p2, hcalfit_l, hcalfit_h, 7);
    mcfit->SetParameters(fidxreports[i].fitParams);
    // Can set the error here, might be relevant later
    // mcfit->SetParErrors(dxreports[i].fitErrors);
    
    // Create a new TH1D or fill an existing one with the values from TF1
    TH1D* hFromTF1_fidx = new TH1D(Form("hFromTF1_fidx_%d", i), "Histogram from fid x TF1", nbins, x_low, x_high);
    
    for (int bin = 1; bin <= nbins; ++bin) {
      double binCenter = fidxreports[i].sliceHistogram->GetBinCenter(bin);
      double funcValue = mcfit->Eval(binCenter);
      hFromTF1_fidx->SetBinContent(bin, funcValue);
    }
    
    int transparentSpring = TColor::GetColorTransparent(kSpring+10, 0.3);
    hFromTF1_fidx->SetLineColor(kSpring+10);
    hFromTF1_fidx->SetLineWidth(0);
    hFromTF1_fidx->SetFillColor(transparentSpring);
    hFromTF1_fidx->SetFillStyle(1001);
    hFromTF1_fidx->Draw("SAME LF2");
    
    canvasSlices_fidx->Update();
  }//endloop over N (slices)
  
  canvasSlices_fidx->Write();
  
  //Create Graphs
  int numValidPoints_fidx = Rsf_vec_fidx.size();

  TGraphErrors* graphErrors_fidx = new TGraphErrors(numValidPoints_fidx);
  for (int i = 0; i < numValidPoints_fidx; ++i) {
    graphErrors_fidx->SetPoint(i, xval_vec_fidx[i], Rsf_vec_fidx[i]);
    graphErrors_fidx->SetPointError(i, 0, Rsferr_vec_fidx[i]);
  }
  graphErrors_fidx->SetMarkerStyle(21);
  graphErrors_fidx->SetMarkerColor(kSpring+10);
  graphErrors_fidx->SetLineColor(kSpring+10);
  graphErrors_fidx->SetTitle("R_{sf} vs. N sig; N sig; R_{sf}");

  TGraph* graph_nev = new TGraph(numValidPoints_fidx);
  for (int i = 0; i < numValidPoints_fidx; ++i) {
    graph_nev->SetPoint(i, xval_vec_fidx[i], nev_vec_fidx[i]);
  }
  graph_nev->SetMarkerStyle(20);
  graph_nev->SetMarkerColor(kGreen-5);
  graph_nev->SetLineColor(kGreen-5);
  graph_nev->SetTitle("Nev Left in Wide Cut vs. N sig; N sig; Nev");

  // // Create a canvas and divide it into 2 sub-pads
  // TCanvas* c1 = new TCanvas("c1", "Stacked Graphs", 800, 600);
  // c1->Divide(1, 2); // Divide canvas into 1 column and 2 rows

  // // Draw the first graph in the first pad
  // c1->cd(1);
  // graphErrors_fidx->Draw("AP");

  // // Draw the second graph in the second pad
  // c1->cd(2);
  // graph_nev->Draw("AP");

  // // Set the Y-axis of the second graph to start at zero
  // Double_t ymax = 1.1 * *std::max_element(nev_vec_fidx.begin(), nev_vec_fidx.end());
  // graph_nev->GetHistogram()->SetMinimum(0); // Set the minimum y-value to zero
  // graph_nev->GetHistogram()->SetMaximum(ymax); // Adjust the maximum y-value

  // // Update the canvas
  // c1->Update();

  //Create overlayed tgraph object
  gROOT->Reset();
  TCanvas *c2 = new TCanvas("c2","",800,500);
  TPad *pad = new TPad("pad","",0,0,1,1);
  //pad->SetFillColor(42);
  pad->SetGrid();
  pad->Draw();
  pad->cd();

  pad->GetFrame()->SetBorderSize(12);

  graphErrors_fidx->GetYaxis()->SetLabelSize(0.04);
  graphErrors_fidx->GetYaxis()->SetLabelFont(42);
  graphErrors_fidx->GetYaxis()->SetTitleOffset(1.0);
  graphErrors_fidx->GetYaxis()->SetTitleSize(0.04);
  graphErrors_fidx->GetYaxis()->SetTitleFont(42); 

  graphErrors_fidx->Draw("AP");

  //create a transparent pad drawn on top of the main pad
  c2->cd();
  TPad *overlay = new TPad("overlay","",0,0,1,1);
  overlay->SetFillStyle(4000);
  overlay->SetFillColor(0);
  overlay->SetFrameFillStyle(4000);
  overlay->Draw();
  overlay->cd();

  graph_nev->GetHistogram()->GetYaxis()->SetLabelOffset(999);
  graph_nev->GetHistogram()->GetYaxis()->SetTickLength(0);
  graph_nev->GetHistogram()->GetYaxis()->SetTitleOffset(999);
  graph_nev->SetTitle("");

  graph_nev->Draw("AP");

  Double_t xmin = graphErrors_fidx->GetXaxis()->GetXmin();
  Double_t xmax = graphErrors_fidx->GetXaxis()->GetXmax();

  Double_t ymin = graph_nev->GetHistogram()->GetMinimum();
  Double_t ymax_nev = graph_nev->GetHistogram()->GetMaximum();

  //Draw an axis on the right side
  TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax_nev, ymin, ymax_nev, 510,"+L");
  axis->SetTitle("ev cut (%)");
  axis->SetTitleColor(kGreen-5);  
  axis->SetTitleOffset(1.0);  
  axis->SetLineColor(kGreen-5);
  axis->SetLabelColor(kGreen-5);
  axis->Draw();

  // Add a legend
  TLegend *leg = new TLegend(0.71, 0.45, 0.86, 0.55);
  leg->AddEntry(graphErrors_fidx, "R_{sf}", "p");
  leg->AddEntry(graph_nev, "N events", "p");
  leg->AddEntry((TObject*)0, Form("ev tot: %d", totalEvents_fidx), "");
  leg->Draw();

  std::string displayname = "Cuts on dx vs fidx (all cuts)";
  util::parseAndDisplayCuts(displayname.c_str(), fidxCuts.c_str());


  /*
  // Now, let's create and draw TGraphErrors for each cut
  // The number of valid points for this cut might be different from N if some were skipped
  int numValidPoints_fidx = Rsf_vec_fidx.size();
  
  // Creating a TGraphErrors for the current cut
  TGraphErrors* graphErrors_fidx = new TGraphErrors(numValidPoints_fidx);
  graphErrors_fidx->SetTitle("Ratio vs fiducial x cut; N prot sig;R_{sf}");
			
  // Filling the TGraphErrors with data
  for (int i = 0; i < numValidPoints_fidx; ++i) {
    graphErrors_fidx->SetPoint(i, xval_vec_fidx[i], Rsf_vec_fidx[i]);
    graphErrors_fidx->SetPointError(i, 0, Rsferr_vec_fidx[i]); // Assuming no error in x
  }

  // Set some graphical attributes
  graphErrors_fidx->SetMarkerStyle(21);
  graphErrors_fidx->SetMarkerColor(kOrange);
  graphErrors_fidx->SetLineColor(kOrange);

  // Drawing the graph
  TCanvas* graphCanvas_fidx = new TCanvas("GraphCanvas_fidx", "Ratio vs fiducial x cut", 1600, 600);
  graphCanvas_fidx->cd();

  // double minY = 0.7; // Minimum y-axis value
  // double maxY = 1.2; // Maximum y-axis value
  // graphErrors->SetMinimum(minY);
  // graphErrors->SetMaximum(maxY);

  graphErrors_fidx->Draw("AP"); // Draw with markers and a line connecting points

  // Save or Write the canvas as needed
  // graphCanvas->SaveAs(Form("TGraphErrors_Cut_%d.png", r)); // Save to a file
  graphCanvas_fidx->Write(); // Or write to an open ROOT file

  //TCanvas *c1; util::parseAndDisplayCuts(c1,fidxCuts.c_str());

  std::string displayname = "Cuts on dx vs fidx (all slices)";
  util::parseAndDisplayCuts(displayname.c_str(), fidxCuts.c_str());

*/

/*
  //Write out the fidy histograms

  //set up vectors for later tgrapherrors
  std::vector<double> Rsf_vec_fidy;
  std::vector<double> Rsferr_vec_fidy;
  std::vector<double> xval_vec_fidy;

  TCanvas* canvasSlices_fidy = new TCanvas("dx slices over fiducial y", "dx slices over fiducial y", 1800, 1200);

  // Assuming N is the total number of histograms to be plotted on the canvas
  int optimalRows_fidy = 1;
  int optimalCols_fidy = 1;

  // Find the nearest square root for N to determine the grid size
  int sqrtN_fidy = (int)std::sqrt(N);
  for (int cols = sqrtN_fidy; cols <= N; ++cols) {
    if (N % cols == 0) { // If cols is a factor of N
      optimalCols_fidy = cols;
      optimalRows_fidy = N / cols;
      break; // Found the optimal layout
    }
  }

  canvasSlices_fidy->Divide(optimalCols_fidy, optimalRows_fidy); // Adjust the division based on reportSamples or desired layout

  TF1* bgfit_fidy[N];

  // loop over slices
  for (int i = 0; i < N; ++i) {
    if (fidyreports[i].sliceHistogram == nullptr)
      continue;

    canvasSlices_fidy->cd(i + 1);

    //catch bad fit ranges and do not write
    std::string tempTitle = fidyreports[i].sliceHistogram->GetTitle();      
    if( tempTitle.compare("Null Histogram")==0 )
      continue;

    //get mc nucleon parameters
    double p_scale = fidyreports[i].fitParams[0];
    double n_scale = fidyreports[i].fitParams[1];
    double p_shift = fidyreports[i].fitParams[2];
    double n_shift = fidyreports[i].fitParams[3];

    double MCev = fidyreports[i].mcpnev + fidyreports[i].mcnnev;

    // Simplify the histogram title
    fidyreports[i].sliceHistogram->SetTitle(Form("%s %0.2f to %0.2f", fidyreports[0].type.c_str(), fidyreports[i].winLow, fidyreports[i].winHigh));
    fidyreports[i].sliceHistogram->SetLineWidth(1);
    fidyreports[i].sliceHistogram->Draw();

    // Retrieve the number of bins, and the x-axis limits of the slice histogram
    int nbins = fidyreports[i].sliceHistogram->GetNbinsX();
    double x_low = fidyreports[i].sliceHistogram->GetXaxis()->GetXmin();
    double x_high = fidyreports[i].sliceHistogram->GetXaxis()->GetXmax();

    //get scale ratio, error, and midpoint for slice
    double ratio = fidyreports[i].scaleratio;
    double error = fidyreports[i].scaleratioerr;
    double xval = fidyreports[i].winLow;

    //fill vectors for later tgrapherrors
    Rsf_vec_fidy.push_back(ratio);
    Rsferr_vec_fidy.push_back(error);
    xval_vec_fidy.push_back(xval);

    // Draw the corresponding MC histograms
    hdx_p = fidyreports[i].slicePHistogram; //update the MC slice for the later overall fit
    TH1D *pMCslice = util::shiftHistogramX(fidyreports[i].slicePHistogram, p_shift);
    pMCslice->Scale(p_scale);

    hdx_n = fidyreports[i].sliceNHistogram; //update the MC slice for the later overall fit
    TH1D *nMCslice = util::shiftHistogramX(fidyreports[i].sliceNHistogram, n_shift);
    nMCslice->Scale(n_scale);
      
    pMCslice->Draw("same");
    nMCslice->Draw("same");

    // Create a legend and add the remaining parameters
    TLegend* legend = new TLegend(0.59, 0.69, 0.89, 0.89); // Adjust the position as needed
    legend->AddEntry((TObject*)0, Form("Nev: %0.0f", fidyreports[i].nev), "");
    legend->AddEntry((TObject*)0, Form("MCev: %0.0f", MCev), "");
    legend->AddEntry((TObject*)0, Form("Ratio: %0.4f #pm %0.4f", ratio, error), "");
    legend->AddEntry((TObject*)0, Form("#chi^{2}/ndf: %0.4f", fidyreports[i].chisqrndf), "");
    legend->Draw();

    bgfit_fidy[i] = new TF1(Form("bgfit_fidy_%d",i),fits::g_p2fit,hcalfit_l,hcalfit_h,3);
    //cout << "Background fit parameter" << endl;
    for (int j=0; j<3; ++j){
      bgfit_fidy[i]->SetParameter(j,fidyreports[i].fitParams[j+4]);
      //cout << "   " << j << " = " << fidyreports[i].fitParams[j+4] << endl;
    }

    bgfit_fidy[i]->SetLineWidth(2);
    bgfit_fidy[i]->Draw("same");

    TF1 *mcfit = new TF1(Form("mcfit_%d", i), fitFullShift_p2, hcalfit_l, hcalfit_h, 7);
    mcfit->SetParameters(fidyreports[i].fitParams);
    // Can set the error here, might be relevant later
    // mcfit->SetParErrors(fidyreports[i].fitErrors);

    // Create a new TH1D or fill an existing one with the values from TF1
    TH1D* hFromTF1_fidy = new TH1D(Form("hFromTF1_fidy_%d", i), "Histogram from fid y TF1", nbins, x_low, x_high);

    for (int bin = 1; bin <= nbins; ++bin) {
      double binCenter = fidyreports[i].sliceHistogram->GetBinCenter(bin);
      double funcValue = mcfit->Eval(binCenter);
      hFromTF1_fidy->SetBinContent(bin, funcValue);
    }

    int transparentBlue = TColor::GetColorTransparent(kBlue, 0.3);
    hFromTF1_fidy->SetLineColor(kBlue);
    hFromTF1_fidy->SetLineWidth(0);
    hFromTF1_fidy->SetFillColor(transparentBlue);
    hFromTF1_fidy->SetFillStyle(1001);
    hFromTF1_fidy->Draw("SAME LF2");

    canvasSlices_fidy->Update();
  }//endloop over N (slices)

  canvasSlices_fidy->Write();

  // Now, let's create and draw TGraphErrors for each cut
  // The number of valid points for this cut might be different from N if some were skipped
  int numValidPoints_fidy = Rsf_vec_fidy.size();

  // Creating a TGraphErrors for the current cut
  TGraphErrors* graphErrors_fidy = new TGraphErrors(numValidPoints_fidy);
  graphErrors_fidy->SetTitle("Ratio vs fiducial y cut; N prot sig;R_{sf}");

  // Filling the TGraphErrors with data
  for (int i = 0; i < numValidPoints_fidy; ++i) {
    graphErrors_fidy->SetPoint(i, xval_vec_fidy[i], Rsf_vec_fidy[i]);
    graphErrors_fidy->SetPointError(i, 0, Rsferr_vec_fidy[i]); // Assuming no error in x
  }

  // Set some graphical attributes
  graphErrors_fidy->SetMarkerStyle(21);
  graphErrors_fidy->SetMarkerColor(kBlue);
  graphErrors_fidy->SetLineColor(kBlue);

  // Drawing the graph
  TCanvas* graphCanvas_fidy = new TCanvas("GraphCanvas_fidy", "Ratio vs fiducial y cut", 1600, 600);
  graphCanvas_fidy->cd();

  // double minY = 0.7; // Minimum y-axis value
  // double maxY = 1.2; // Maximum y-axis value
  // graphErrors->SetMinimum(minY);
  // graphErrors->SetMaximum(maxY);

  graphErrors_fidy->Draw("AP"); // Draw with markers and a line connecting points

  // Save or Write the canvas as needed
  // graphCanvas->SaveAs(Form("TGraphErrors_Cut_%d.png", r)); // Save to a file
  graphCanvas_fidy->Write(); // Or write to an open ROOT file
*/

  cout << "Analysis complete. Output file written to " << foutPath  << endl;

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

//Function to fit then fine fit and return all fit parameters
std::vector<std::pair<double, double>> fitAndFineFit(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, double pshift, double nshift, const std::string& fitOptions = "RBMQ0") {
  TF1* fit = new TF1(fitName.c_str(), fitFormula.c_str(), hcalfit_l, hcalfit_h, paramCount);
  //fit->SetNpx(5000);
  for (int i=0; i<paramCount; ++i){ //reset parameters/errors for this set
    fit->SetParameter(i,0);
    fit->SetParError(i,0);
  }
  if(!(pshift==0&&nshift==0)){
    fit->FixParameter(2,pshift);
    fit->FixParameter(3,nshift);
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
  if(!(pshift==0&&nshift==0)){
    fit->FixParameter(2,pshift);
    fit->FixParameter(3,nshift);
  }

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
