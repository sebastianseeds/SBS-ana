//sseeds 4.22.24

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

//Fit range override options
double hcalfit_l; //lower fit/bin limit for hcal dx plots (m)
double hcalfit_h; //upper fit/bin limit for hcal dx plots (m)

double xrange_min_dy = -1.5;
double xrange_max_dy = 1.5;

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
std::string addbc(std::string input);
std::vector<std::string> split(const std::string &s, char delimiter);
std::map<std::string, std::string> getRowContents(const std::string &filePath, int kine, int mag, const std::string &target, const std::vector<std::string> &excludeKeys);

//main. kine=kinematic, mag=fieldsetting, pass=pass#, sb_min/max=sidebandlimits, shiftX=shifttodxdata, N=cutvarsliceN, sliceCutMax=NCutsFromZeroTosliceCutMax
void dyStability(int kine=9, 
		 int mag=70, 
		 int pass=2, 
		 int N=12,
		 double mean=-0.05,
		 double sigma=0.1,
		 double start_sigma=1,
		 double max_sigma=10,
		 std::string BG="pol2",
		 bool bestclus=true, 
		 bool thin=true,
		 bool wide=false,
		 bool effz=true,
		 bool alt = true) {
  //set draw params
  gStyle->SetPalette(kViridis);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);

  if(start_sigma>max_sigma){
    cout << "ERROR: starting sigma value for analysis greater than max sigma value." << endl;
    return;
  }

  //load json file
  JSONManager *jmgr = new JSONManager("../../config/syst.json");

  // Get dx plot details
  hcalfit_l = jmgr->GetValueFromSubKey<double>("hcalfit_l", Form("sbs%d_%d", kine, mag));
  hcalfit_h = jmgr->GetValueFromSubKey<double>("hcalfit_h", Form("sbs%d_%d", kine, mag));
  
  //get new cuts from .csv
  std::string cutsheet_path = "/w/halla-scshelf2102/sbs/seeds/ana/data/p2_cutset.csv";
  std::string target = "ld2";
  std::vector<std::string> excludeKeys = {};

  // Get the row contents
  std::map<std::string, std::string> rowContents = getRowContents(cutsheet_path, kine, mag, target, excludeKeys);
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
  
  // //Get tight elastic cuts
  // std::string globalcuts = jmgr->GetValueFromSubKey_str( Form("post_cuts_p%d",pass), Form("sbs%d_%d",kine,mag) );

  // cout << "Loaded tight cuts: " << globalcuts << endl;

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

  std::string basePath = "/lustre19/expphy/volatile/halla/sbs/seeds";
  std::string finPath = Form("%s/gmn_analysis/dx_correlations%s_sbs%d_mag%d_pass%d%s%s%s.root", basePath.c_str(), bestclus_word.c_str(), kine, mag, pass, thin_word.c_str(), wide_word.c_str(), effz_word.c_str());
  std::string fmcinPath = Form("%s/gmn_analysis/dx_mc_sbs%d_mag%d_pass%d%s%s%s%s.root", basePath.c_str(), kine, mag, pass, thin_word.c_str(), alt_word.c_str(), wide_word.c_str(), effz_word.c_str());
  std::string foutPath = Form("%s/gmn_analysis/dystab_gmn_sbs%d_mag%d_pass%d%s%s%s.root", basePath.c_str(), kine, mag, pass, bestclus_word.c_str(), wide_word.c_str(), effz_word.c_str());

  TFile* inputFile = new TFile(finPath.c_str());
  if (!inputFile || inputFile->IsZombie()) 
    handleError(inputFile,"inputFile");

  //Get all histograms
  TH1D *hdx_data = dynamic_cast<TH1D*>(inputFile->Get("hdx_allcut"));
  hdx_data->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  std::string dyhist_s = "hist_dy";
  if(bestclus)
    dyhist_s += "_bc";

  TH2D *hdx_vs_dy_data = dynamic_cast<TH2D*>(inputFile->Get(dyhist_s.c_str()));
  hdx_vs_dy_data->GetYaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  int hdx_vs_dy_data_N = hdx_vs_dy_data->GetEntries();
  std::string dyCuts = hdx_vs_dy_data->GetTitle();
  cout << endl << "Opened dx vs dy with cuts: " << dyCuts << endl << endl; 

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

  TH2D *hdx_vs_dy_n = dynamic_cast<TH2D*>(inputFileMC->Get("hist_dy_n"));
  hdx_vs_dy_n->GetYaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  int hdx_vs_dy_n_N = hdx_vs_dy_n->GetEntries();

  TH2D *hdx_vs_dy_p = dynamic_cast<TH2D*>(inputFileMC->Get("hist_dy_p"));
  hdx_vs_dy_p->GetYaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  int hdx_vs_dy_p_N = hdx_vs_dy_p->GetEntries();

  hdx_inel = dynamic_cast<TH1D*>(inputFileMC->Get("hdx_inel"));
  hdx_inel->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  hdx_inel->Scale(10e33); //account for lack of overall normalization from g4sbs generator
  if (!hdx_p || !hdx_n || !hdx_inel) 
    handleError(inputFile,inputFileMC,"hdx");

  TFile* outputFile = new TFile(foutPath.c_str(), "RECREATE");

  TCanvas *canvas = new TCanvas("canvas", "Data and MC Histograms", 1800, 1200);
  canvas->Divide(3, 3);  // Adjust the grid size according to the number of histograms

  // Pad 1: hdx_data
  canvas->cd(1);
  hdx_data->Draw();

  // Pad 2: hdx_vs_dy_data
  canvas->cd(2);
  hdx_vs_dy_data->Draw("COLZ");

  // Pad 4: hdx_dyanti
  canvas->cd(3);
  hdx_dyanti->Draw();

  // Pad 5: hdx_coinanti
  canvas->cd(4);
  hdx_coinanti->Draw();

  // Pad 6: hdx_p
  canvas->cd(5);
  hdx_p->Draw();

  // Pad 7: hdx_n
  canvas->cd(6);
  hdx_n->Draw();

  // Pad 8: hdx_vs_dy_n
  canvas->cd(7);
  hdx_vs_dy_n->Draw("COLZ");

  // Pad 10: hdx_vs_dy_p
  canvas->cd(8);
  hdx_vs_dy_p->Draw("COLZ");

  // Pad 12: hdx_inel
  canvas->cd(9);
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

  std::vector<ReportData> dyreports;
  dyreports.resize(N);

  // Loop over dx vs dy
  hdx_vs_dy_data->GetXaxis()->SetRangeUser(xrange_min_dy, xrange_max_dy);
  hdx_vs_dy_data->GetYaxis()->SetRangeUser(hcalfit_l, hcalfit_h);

  // Get plot limits
  int xbinLow_dy = hdx_vs_dy_data->GetXaxis()->FindBin(xrange_min_dy);
  int xbinHigh_dy = hdx_vs_dy_data->GetXaxis()->FindBin(xrange_max_dy); 

  // Get cut region
  int xbinmean_dy = hdx_vs_dy_data->GetXaxis()->FindBin(mean);
  double xCutMax = mean + (max_sigma*sigma);
  double xCutMin = mean - (max_sigma*sigma);
  int xbinCutMax_dy = hdx_vs_dy_data->GetXaxis()->FindBin(xCutMax);
  int xbinCutMin_dy = hdx_vs_dy_data->GetXaxis()->FindBin(xCutMin);
  double startcut = start_sigma * sigma;
  int xbinCutStart_low_dy = hdx_vs_dy_data->GetXaxis()->FindBin(mean-startcut);
  int xbinCutStart_high_dy = hdx_vs_dy_data->GetXaxis()->FindBin(mean+startcut);

  cout << "Bin at starting sigma RHS = " << xbinCutStart_high_dy << endl;

  int NbinsRight = xbinCutMax_dy - xbinCutStart_high_dy;
  
  cout << "Mean bin " << xbinmean_dy << " with rightmost bin " << xbinCutMax_dy << ".." << endl;

  cout << "N bins right of mean " << NbinsRight << ".." << endl;

  //Adjust total number of bins if not enough available on the TH2D
  if(NbinsRight<N)
    N=NbinsRight;

  int xbinRanges_dy = NbinsRight / N;

  cout << "Each added range per cut: " << xbinRanges_dy << ".." << endl;

  cout << "Added bins in x per cut: " << 2*xbinRanges_dy << " on N = " << N << " cuts." << endl;

  // Project the TH2D onto a TH1D for the specified x-axis range
  TH1D* fullProjY_dy = hdx_vs_dy_data->ProjectionY("_py", xbinLow_dy, xbinHigh_dy);
  TH1D* wideProjY_dy = hdx_vs_dy_data->ProjectionY("_wpy");

  // Use Integral() to get the total number of events within the specified range
  int totalEvents_dy = fullProjY_dy->Integral();

  std::pair<double,double> wideQualp2;

  auto widep2par_vector = fitAndFineFit(wideProjY_dy, "sliceFitWide", "fitFullShift_p2", 7, hcalfit_l, hcalfit_h, wideQualp2, 0, 0, fitopt.c_str());

  double fixed_pshift = widep2par_vector[2].first;
  double fixed_nshift = widep2par_vector[3].first;

  cout << "Fit over full dy range yield shifts p:n -> " << fixed_pshift << ":" << fixed_nshift << endl;

  vector<double> cutx_low_dy;
  vector<double> cutx_high_dy;

  cout << "Looping over dy ranges.." << endl;
  for (int i = 0; i < N; ++i) {

    int binStart = xbinCutStart_low_dy - (i+1)*xbinRanges_dy;
    int binEnd = xbinCutStart_high_dy + (i+1)*xbinRanges_dy;

    if(binStart<xbinLow_dy || binEnd>xbinHigh_dy){
      cout << "ERROR: cut ranges exceed length of TH2D in x." << endl;
      return;

    }

    cutx_low_dy.push_back(hdx_vs_dy_data->GetXaxis()->GetBinLowEdge(binStart));
    cutx_high_dy.push_back(hdx_vs_dy_data->GetXaxis()->GetBinLowEdge(binEnd+1));

    cout << "Cut range from x=" << hdx_vs_dy_data->GetXaxis()->GetBinLowEdge(binStart) << " to " << hdx_vs_dy_data->GetXaxis()->GetBinLowEdge(binEnd+1) << " (bin " << binStart << " to bin " << binEnd << ")" << endl;

    //Get the data slices
    TH1D *hdx_slice = hdx_vs_dy_data->ProjectionY(Form("hdx_dy_slice_%d", i+1), binStart, binEnd);
    TH1D *dx_dy_slice = hdx_vs_dy_data->ProjectionY(Form("dxfidcut_%d", i+1), binStart, binEnd);

    //Get the MC slices
    TH1D *dx_dy_p_slice = hdx_vs_dy_p->ProjectionY(Form("dxfidcut_p_%d", i+1), binStart, binEnd);
    TH1D *dx_dy_n_slice = hdx_vs_dy_n->ProjectionY(Form("dxfidcut_n_%d", i+1), binStart, binEnd);
    hdx_p = hdx_vs_dy_p->ProjectionY(Form("hdx_dy_slice_p_%d", i+1), binStart, binEnd);
    hdx_n = hdx_vs_dy_n->ProjectionY(Form("hdx_dy_slice_n_%d", i+1), binStart, binEnd);

    //Get total entries on slice
    int slice_nEntries = dx_dy_slice->Integral();
    int slicemcp_nEntries = dx_dy_p_slice->Integral();
    int slicemcn_nEntries = dx_dy_n_slice->Integral();

    //Rsf extraction using second order poly fit
    std::pair<double,double> sliceQualp2;

    auto slicep2par_vector = fitAndFineFit(hdx_slice, "sliceFit", "fitFullShift_p2", 7, hcalfit_l, hcalfit_h, sliceQualp2, fixed_pshift, fixed_nshift, fitopt.c_str());

    double csndf = sliceQualp2.first/sliceQualp2.second;

    double scaleratio = slicep2par_vector[1].first / slicep2par_vector[0].first;
    double scaleratio_err = sqrt(pow(slicep2par_vector[1].second/slicep2par_vector[1].first, 2) + pow(slicep2par_vector[0].second/slicep2par_vector[0].first, 2)) * scaleratio;

    double xshift_p = slicep2par_vector[2].first;
    double xshift_n = slicep2par_vector[3].first;

    //Write fit results to struct
    dyreports[i] = ReportData( dx_dy_slice, 
				 dx_dy_p_slice,
				 dx_dy_n_slice,
				 cutx_low_dy[i],
				 cutx_high_dy[i],
				 slice_nEntries,
				 slicemcp_nEntries,
				 slicemcn_nEntries,
				 scaleratio,
				 scaleratio_err,
				 csndf,
				 "dx vs dy" );

    for (int par=0; par<7; ++par){
      dyreports[i].fitParams[par] = slicep2par_vector[par].first;
      dyreports[i].fitErrors[par] = slicep2par_vector[par].second;

      
    }
   
  }
 
  //Write out the dx histograms
  
  //get a clone for displays
  TH2D* cloneddy = (TH2D*)hdx_vs_dy_data->Clone("cloneddy");
  double dy_ymin = hcalfit_l;
  double dy_ymax = hcalfit_h;

  TCanvas* c0 = new TCanvas("c0", "Slice Locations", 800,600);
  c0->cd();
  cloneddy->SetTitle("dy vs dx");
  cloneddy->Draw("colz");

  int numberOfLines = cutx_low_dy.size(); // Assuming cutx_low_dy and cutx_high_dy are same size
  int gradientSteps = numberOfLines;
  int* colorGradient = new int[gradientSteps];
  double r1 = 1.0, g1 = 0.0, b1 = 0.0; // Starting color (red)
  double r2 = 0.0, g2 = 0.0, b2 = 1.0; // Ending color (blue)
  for (int i = 0; i < gradientSteps; ++i) {
    colorGradient[i] = TColor::GetColor(
  					static_cast<Float_t>(r1 + (r2 - r1) * i / (gradientSteps - 1)),
  					static_cast<Float_t>(g1 + (g2 - g1) * i / (gradientSteps - 1)),
  					static_cast<Float_t>(b1 + (b2 - b1) * i / (gradientSteps - 1))
  					);
  }

  int lineThickness = 1; // Start line thickness

  // Drawing lines with increasing thickness for lower limits
  for (int i = 0; i < numberOfLines; ++i) {
    TLine* line = new TLine(cutx_low_dy[i], dy_ymin, cutx_low_dy[i], dy_ymax);
    line->SetLineColor(colorGradient[i]); // Set line color dynamically
    line->SetLineWidth(2); // Set line thickness dynamically
    line->Draw();
    ++lineThickness; // Increase line thickness for next line
  }

  lineThickness = 1; // Reset line thickness for upper limits

  // Drawing lines with increasing thickness for upper limits
  for (int i = 0; i < numberOfLines; ++i) {
    TLine* line = new TLine(cutx_high_dy[i], dy_ymin, cutx_high_dy[i], dy_ymax);
    line->SetLineColor(colorGradient[i]); // Set line color dynamically
    line->SetLineWidth(2); // Set line thickness dynamically
    line->Draw();
    ++lineThickness; // Increase line thickness for next line
  }

  c0->Update();

  // for (auto& cut : cutx_low_dy){
  //   TLine* line = new TLine(cut, dy_ymin, cut, dy_ymax);
  //   line->SetLineColor(kRed); // Set line color to red
  //   line->Draw();
  // }

  // for (auto& cut : cutx_high_dy){
  //   TLine* line = new TLine(cut, dy_ymin, cut, dy_ymax);
  //   line->SetLineColor(kRed); // Set line color to red
  //   line->Draw();
  // }

  c0->Update();

  //set up vectors for later tgrapherrors
  std::vector<double> Rsf_vec_dy;
  std::vector<double> Rsferr_vec_dy;
  std::vector<double> xval_vec_dy;
  std::vector<double> sig_vec_dy;
  std::vector<double> nev_vec_dy;

  TCanvas* canvasSlices_dy = new TCanvas("dx slices over dy", "dx slices over dy", 1800, 1200);
  
  // Assuming N is the total number of histograms to be plotted on the canvas
  int optimalRows_dy = 1;
  int optimalCols_dy = 1;
  
  // Find the nearest square root for N to determine the grid size
  int sqrtN_dy = (int)std::sqrt(N);
  for (int cols = sqrtN_dy; cols <= N; ++cols) {
    if (N % cols == 0) { // If cols is a factor of N
      optimalCols_dy = cols;
      optimalRows_dy = N / cols;
      break; // Found the optimal layout
    }
  }
  
  canvasSlices_dy->Divide(optimalCols_dy, optimalRows_dy); // Adjust the division based on reportSamples or desired layout
  
  TF1* bgfit_dy[N];
  
  // loop over slices
  for (int i = 0; i < N; ++i) {
    if (dyreports[i].sliceHistogram == nullptr)
      continue;
    
    canvasSlices_dy->cd(i + 1);
    
    //catch bad fit ranges and do not write
    std::string tempTitle = dyreports[i].sliceHistogram->GetTitle();      
    if( tempTitle.compare("Null Histogram")==0 )
      continue;
    
    //get mc nucleon parameters
    double p_scale = dyreports[i].fitParams[0];
    double n_scale = dyreports[i].fitParams[1];
    double p_shift = dyreports[i].fitParams[2];
    double n_shift = dyreports[i].fitParams[3];
    
    double MCev = dyreports[i].mcpnev + dyreports[i].mcnnev;
    
    // Simplify the histogram title
    dyreports[i].sliceHistogram->SetTitle(Form("%s %0.4f to %0.4f", dyreports[0].type.c_str(), dyreports[i].winLow, dyreports[i].winHigh));
    dyreports[i].sliceHistogram->SetLineWidth(1);
    dyreports[i].sliceHistogram->Draw();
    
    // Retrieve the number of bins, and the x-axis limits of the slice histogram
    int nbins = dyreports[i].sliceHistogram->GetNbinsX();
    double x_low = dyreports[i].sliceHistogram->GetXaxis()->GetXmin();
    double x_high = dyreports[i].sliceHistogram->GetXaxis()->GetXmax();
    
    //get scale ratio, error, and midpoint for slice
    double ratio = dyreports[i].scaleratio;
    double error = dyreports[i].scaleratioerr;
    double xval = dyreports[i].winLow;
    double sig = (i+1)*((max_sigma-start_sigma)/N)+start_sigma;
    double nev = (totalEvents_dy - dyreports[i].nev)/totalEvents_dy*100;

    //fill vectors for later tgrapherrors
    Rsf_vec_dy.push_back(ratio-0.025); //offset correction to account for included W2 cut (+0.0155)
    Rsferr_vec_dy.push_back(error);
    xval_vec_dy.push_back(xval);
    sig_vec_dy.push_back(sig);
    nev_vec_dy.push_back(nev); //Total remaining stats

    // Draw the corresponding MC histograms
    hdx_p = dyreports[i].slicePHistogram; //update the MC slice for the later overall fit
    TH1D *pMCslice = util::shiftHistogramX(dyreports[i].slicePHistogram, p_shift);
    pMCslice->Scale(p_scale);
    
    hdx_n = dyreports[i].sliceNHistogram; //update the MC slice for the later overall fit
    TH1D *nMCslice = util::shiftHistogramX(dyreports[i].sliceNHistogram, n_shift);
    nMCslice->Scale(n_scale);
    
    pMCslice->Draw("same");
    nMCslice->Draw("same");

    //cout << "bin " << i << " pshift " << p_shift << " nshift " << n_shift << endl;
    
    // Create a legend and add the remaining parameters
    TLegend* legend = new TLegend(0.49, 0.59, 0.89, 0.89); // Adjust the position as needed
    legend->SetMargin(0.1);
    legend->AddEntry((TObject*)0, Form("Nev: %0.0f", dyreports[i].nev), "");
    legend->AddEntry((TObject*)0, Form("MCev: %0.0f", MCev), "");
    legend->AddEntry((TObject*)0, Form("Ratio: %0.2f", ratio), "");
    legend->AddEntry((TObject*)0, Form("Shift p/n: %0.2f/%0.2f", p_shift, n_shift), "");
    legend->AddEntry((TObject*)0, Form("#chi^{2}/ndf: %0.2f", dyreports[i].chisqrndf), "");
    legend->Draw();
    
    bgfit_dy[i] = new TF1(Form("bgfit_dy_%d",i),fits::g_p2fit_cd,hcalfit_l,hcalfit_h,3);
    //cout << "Background fit parameter" << endl;
    for (int j=0; j<3; ++j){
      bgfit_dy[i]->SetParameter(j,dyreports[i].fitParams[j+4]);
      //cout << "   " << j << " = " << dxreports[i].fitParams[j+4] << endl;
    }
    
    bgfit_dy[i]->SetLineWidth(2);
    bgfit_dy[i]->Draw("same");
    
    TF1 *mcfit = new TF1(Form("mcfit_%d", i), fitFullShift_p2, hcalfit_l, hcalfit_h, 7);
    mcfit->SetParameters(dyreports[i].fitParams);
    // Can set the error here, might be relevant later
    // mcfit->SetParErrors(dxreports[i].fitErrors);
    
    // Create a new TH1D or fill an existing one with the values from TF1
    TH1D* hFromTF1_dy = new TH1D(Form("hFromTF1_dy_%d", i), "Histogram from fid x TF1", nbins, x_low, x_high);
    
    for (int bin = 1; bin <= nbins; ++bin) {
      double binCenter = dyreports[i].sliceHistogram->GetBinCenter(bin);
      double funcValue = mcfit->Eval(binCenter);
      hFromTF1_dy->SetBinContent(bin, funcValue);
    }
    
    int transparentBlue = TColor::GetColorTransparent(kBlue-9, 0.3);
    hFromTF1_dy->SetLineColor(kBlue-9);
    hFromTF1_dy->SetLineWidth(0);
    hFromTF1_dy->SetFillColor(transparentBlue);
    hFromTF1_dy->SetFillStyle(1001);
    hFromTF1_dy->Draw("SAME LF2");
    
    canvasSlices_dy->Update();
  }//endloop over N (slices)
  
  canvasSlices_dy->Write();
  
  // Now, let's create and draw TGraphErrors for each cut
  // The number of valid points for this cut might be different from N if some were skipped
  int numValidPoints_dy = Rsf_vec_dy.size();
  
  TGraphErrors* graphErrors_dy = new TGraphErrors(numValidPoints_dy);
  for (int i = 0; i < numValidPoints_dy; ++i) {
    graphErrors_dy->SetPoint(i, sig_vec_dy[i], Rsf_vec_dy[i]);
    graphErrors_dy->SetPointError(i, 0, Rsferr_vec_dy[i]);
  }
  graphErrors_dy->SetMarkerStyle(21);
  graphErrors_dy->SetMarkerColor(kBlue-9);
  graphErrors_dy->SetLineColor(kBlue-9);
  graphErrors_dy->SetTitle(Form("R_{sf} vs. N sig (#sigma = %.3f m, mean = %.3f m); N sig; R_{sf}",sigma,mean));

  TGraph* graph_nev = new TGraph(numValidPoints_dy);
  for (int i = 0; i < numValidPoints_dy; ++i) {
    graph_nev->SetPoint(i, sig_vec_dy[i], nev_vec_dy[i]);
  }
  graph_nev->SetMarkerStyle(20);
  graph_nev->SetMarkerColor(kGreen-5);
  graph_nev->SetLineColor(kGreen-5);
  graph_nev->SetTitle(Form("Nev Left in Wide Cut vs. N sig (#sigma = %.3f m, mean = %.3f m); N sig; Nev",sigma,mean));

  // // Create a canvas and divide it into 2 sub-pads
  // TCanvas* c1 = new TCanvas("c1", "Stacked Graphs", 800, 600);
  // c1->Divide(1, 2); // Divide canvas into 1 column and 2 rows

  // // Draw the first graph in the first pad
  // c1->cd(1);
  // graphErrors_dy->Draw("AP");

  // // Draw the second graph in the second pad
  // c1->cd(2);
  // graph_nev->Draw("AP");

  // // Set the Y-axis of the second graph to start at zero
  // Double_t ymax = 1.1 * *std::max_element(nev_vec_dy.begin(), nev_vec_dy.end());
  // graph_nev->GetHistogram()->SetMinimum(0); // Set the minimum y-value to zero
  // graph_nev->GetHistogram()->SetMaximum(ymax); // Adjust the maximum y-value

  // // Update the canvas
  // c1->Update();

  //Create overlayed tgraph object
  gROOT->Reset();
  TCanvas *c2 = new TCanvas("c2","",800,600);
  TPad *pad = new TPad("pad","",0,0,1,1);
  //pad->SetFillColor(42);
  pad->SetGrid();
  pad->Draw();
  pad->cd();

  pad->GetFrame()->SetBorderSize(12);

  graphErrors_dy->GetYaxis()->SetLabelSize(0.04);
  graphErrors_dy->GetYaxis()->SetLabelFont(42);
  graphErrors_dy->GetYaxis()->SetTitleOffset(1.2);
  graphErrors_dy->GetYaxis()->SetTitleSize(0.04);
  graphErrors_dy->GetYaxis()->SetTitleFont(42); 

  graphErrors_dy->Draw("AP");

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

  Double_t xmin = graphErrors_dy->GetXaxis()->GetXmin();
  Double_t xmax = graphErrors_dy->GetXaxis()->GetXmax();

  Double_t ymin = graph_nev->GetHistogram()->GetMinimum();
  Double_t ymax_nev = graph_nev->GetHistogram()->GetMaximum();

  //Draw an axis on the right side
  TGaxis *axis = new TGaxis(xmax, ymin, xmax, ymax_nev, ymin, ymax_nev, 510,"+L");
  axis->SetTitle("ev cut (%)");
  axis->SetTitleColor(kGreen-5);  
  axis->SetTitleOffset(0.8);  
  axis->SetLineColor(kGreen-5);
  axis->SetLabelColor(kGreen-5);
  axis->Draw();

  // Add a legend
  TLegend *leg = new TLegend(0.71, 0.45, 0.86, 0.55);
  leg->AddEntry(graphErrors_dy, "R_{sf}", "p");
  leg->AddEntry(graph_nev, "N events", "p");
  leg->AddEntry((TObject*)0, Form("ev tot: %d", totalEvents_dy), "");
  leg->Draw();

  std::string displayname = "Cuts on dx vs dy (all sigma ranges)";
  util::parseAndDisplayCuts(displayname.c_str(), dyCuts.c_str());

  cout << "Analysis complete. Output file written to " << foutPath  << endl;



  /*

  // Creating a TGraphErrors for the current cut
  TGraphErrors* graphErrors_dy = new TGraphErrors(numValidPoints_dy);
  graphErrors_dy->SetTitle("Ratio vs dy cut range; N sig;R_{sf}");
			
  // Filling the TGraphErrors with data
  for (int i = 0; i < numValidPoints_dy; ++i) {
    graphErrors_dy->SetPoint(i, sig_vec_dy[i], Rsf_vec_dy[i]);
    graphErrors_dy->SetPointError(i, 0, Rsferr_vec_dy[i]); // Assuming no error in x
  }

  // Set some graphical attributes
  graphErrors_dy->SetMarkerStyle(21);
  graphErrors_dy->SetMarkerColor(kBlue);
  graphErrors_dy->SetLineColor(kBlue);

  // Drawing the graph
  TCanvas* graphCanvas_dy = new TCanvas("GraphCanvas_dy", "Ratio vs dy cut range", 1600, 600);
  graphCanvas_dy->cd();

  // double minY = 0.7; // Minimum y-axis value
  // double maxY = 1.2; // Maximum y-axis value
  // graphErrors->SetMinimum(minY);
  // graphErrors->SetMaximum(maxY);

  graphErrors_dy->Draw("AP"); // Draw with markers and a line connecting points

  // Save or Write the canvas as needed
  // graphCanvas->SaveAs(Form("TGraphErrors_Cut_%d.png", r)); // Save to a file
  graphCanvas_dy->Write(); // Or write to an open ROOT file

  //TCanvas *c1; util::parseAndDisplayCuts(c1,dyCuts.c_str());

  std::string displayname = "Cuts on dx vs dy (all ranges)";
  util::parseAndDisplayCuts(displayname.c_str(), dyCuts.c_str());


  cout << "Analysis complete. Output file written to " << foutPath  << endl;

  */

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

  //fit->SetParLimits(6, -FLT_MAX, 0);

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

  //fit->SetParLimits(6, -FLT_MAX, 0);

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
