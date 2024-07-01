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
double hcalfit_l; //lower fit/bin limit for hcal dx plots (m) sbs4, 50p
double hcalfit_h; //upper fit/bin limit for hcal dx plots (m) sbs4, 50p

double xrange_min_bb_ps_e = 0.;
double xrange_max_bb_ps_e = 10;

//Fit options
std::string fitopt = "RMQ0";

//fit min entries
int minEvents = 500;

//Total fits using Interpolate with elastic signal histo and 4th order poly fit to bg
TH1D *hdx_p;
TH1D *hdx_n;

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


// Forward declarations
void handleError(TFile *file1, TFile *file2, std::string marker);
void handleError(TFile *file1, std::string marker);
std::vector<std::pair<double, double>> fitAndFineFit(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, double pshift, double nshift, const std::string& fitOptions = "RBMQ0");
std::string addbc(std::string input);
std::vector<std::string> split(const std::string &s, char delimiter);
std::map<std::string, std::string> getRowContents(const std::string &filePath, int kine, int mag, const std::string &target, const std::vector<std::string> &excludeKeys);

//main. kine=kinematic, mag=fieldsetting, pass=pass#, sb_min/max=sidebandlimits, shiftX=shifttodxdata, N=cutvarsliceN, sliceCutMax=NCutsFromZeroTosliceCutMax
void pseStability(int kine=9, 
		  int mag=70, 
		  int pass=2, 
		  int N=9,
		  double sliceCutMin=0.20,
		  double sliceCutMax=0.34,
		  std::string BG="pol2",
		  bool backwards=false,
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
  std::string foutPath = Form("%s/gmn_analysis/psestab_gmn_sbs%d_mag%d_pass%d%s%s%s.root", basePath.c_str(), kine, mag, pass, bestclus_word.c_str(), wide_word.c_str(), effz_word.c_str());

  TFile* inputFile = new TFile(finPath.c_str());
  if (!inputFile || inputFile->IsZombie()) 
    handleError(inputFile,"inputFile");

  //Get all histograms
  TH1D *hdx_data = dynamic_cast<TH1D*>(inputFile->Get("hdx_allcut"));
  hdx_data->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  TH2D *hdx_vs_bb_ps_e_data = dynamic_cast<TH2D*>(inputFile->Get("hist_bb_ps_e"));
  hdx_vs_bb_ps_e_data->GetYaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  int hdx_vs_bb_ps_e_data_N = hdx_vs_bb_ps_e_data->GetEntries();
  std::string bb_ps_eCuts = hdx_vs_bb_ps_e_data->GetTitle();
  cout << endl << "Opened dx vs hcal E with cuts: " << bb_ps_eCuts << endl << endl;

  if (!hdx_data) 
    handleError(inputFile,"hdx_data");

  TFile* inputFileMC = new TFile(fmcinPath.c_str(), "READ");
  if (!inputFileMC || inputFileMC->IsZombie())
    handleError(inputFile,inputFileMC,"inputFileMC");

  //fix the interpolate functions
  hdx_p = dynamic_cast<TH1D*>(inputFileMC->Get("hdx_p"));
  hdx_p->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  hdx_n = dynamic_cast<TH1D*>(inputFileMC->Get("hdx_n"));
  hdx_n->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  TH2D *hdx_vs_bb_ps_e_n = dynamic_cast<TH2D*>(inputFileMC->Get("hist_bb_ps_e_n"));
  hdx_vs_bb_ps_e_n->GetYaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  int hdx_vs_bb_ps_e_n_N = hdx_vs_bb_ps_e_n->GetEntries();

  TH2D *hdx_vs_bb_ps_e_p = dynamic_cast<TH2D*>(inputFileMC->Get("hist_bb_ps_e_p"));
  hdx_vs_bb_ps_e_p->GetYaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  int hdx_vs_bb_ps_e_p_N = hdx_vs_bb_ps_e_p->GetEntries();

  if (!hdx_p || !hdx_n) 
    handleError(inputFile,inputFileMC,"hdx");

  TFile* outputFile = new TFile(foutPath.c_str(), "RECREATE");

  TCanvas *canvas = new TCanvas("canvas", "Data and MC Histograms", 1800, 1200);
  canvas->Divide(3, 2);  // Adjust the grid size according to the number of histograms

  // Pad 1: hdx_data
  canvas->cd(1);
  hdx_data->Draw();

  // Pad 2: hdx_vs_bb_ps_e_data
  canvas->cd(2);
  hdx_vs_bb_ps_e_data->Draw("COLZ");

  // Pad 6: hdx_p
  canvas->cd(3);
  hdx_p->Draw();

  // Pad 7: hdx_n
  canvas->cd(4);
  hdx_n->Draw();

  // Pad 8: hdx_vs_bb_ps_e_n
  canvas->cd(5);
  hdx_vs_bb_ps_e_n->Draw("COLZ");

  // Pad 10: hdx_vs_bb_ps_e_p
  canvas->cd(6);
  hdx_vs_bb_ps_e_p->Draw("COLZ");

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

  /////////////////
  /////////////////
  /////////////////
  //dependent ranges
  /////////////////
  /////////////////

  std::vector<ReportData> bb_ps_ereports;
  bb_ps_ereports.resize(N);

  // Loop over dx vs fiducial x
  hdx_vs_bb_ps_e_data->GetXaxis()->SetRangeUser(xrange_min_bb_ps_e, xrange_max_bb_ps_e);
  hdx_vs_bb_ps_e_data->GetYaxis()->SetRangeUser(hcalfit_l, hcalfit_h);

  // Assuming current_llim and current_ulim are set to the desired x-axis range
  int xbinLow_bb_ps_e = hdx_vs_bb_ps_e_data->GetXaxis()->FindBin(xrange_min_bb_ps_e);
  int xbinWideLow_bb_ps_e = hdx_vs_bb_ps_e_data->GetXaxis()->FindBin(0.5); //Set beneath the cut to map the region
  int xbinHigh_bb_ps_e = hdx_vs_bb_ps_e_data->GetXaxis()->FindBin(xrange_max_bb_ps_e);
  int xbinCutMin_bb_ps_e = hdx_vs_bb_ps_e_data->GetXaxis()->FindBin(sliceCutMin);
  int xbinCutMax_bb_ps_e = hdx_vs_bb_ps_e_data->GetXaxis()->FindBin(sliceCutMax);
  if(sliceCutMax==0.)
    xbinCutMax_bb_ps_e = xbinHigh_bb_ps_e;
  int xbinBegin;

  // Bins per range
  int xbinRanges_bb_ps_e;
  if(sliceCutMax>0 && sliceCutMin==0){
    int modrange = xbinCutMax_bb_ps_e - xbinLow_bb_ps_e;
    xbinBegin = xbinLow_bb_ps_e;
    if(modrange<N)
      N=modrange;
    xbinRanges_bb_ps_e = modrange / N;
  }else if(sliceCutMax==0 && sliceCutMin>0){
    int modrange = xbinHigh_bb_ps_e - xbinCutMin_bb_ps_e;
    xbinBegin = xbinCutMin_bb_ps_e;
    if(modrange<N)
      N=modrange;
    xbinRanges_bb_ps_e = modrange / N;
  }else if(sliceCutMax>0 && sliceCutMin>0){
    int modrange = xbinCutMax_bb_ps_e - xbinCutMin_bb_ps_e;
    xbinBegin = xbinCutMin_bb_ps_e;
    if(modrange<N)
      N=modrange;
    xbinRanges_bb_ps_e = modrange / N;
  }else{
    int modrange = xbinHigh_bb_ps_e - xbinLow_bb_ps_e;
    xbinBegin = xbinLow_bb_ps_e;
    if(modrange<N)
      N=modrange;
    xbinRanges_bb_ps_e = modrange / N;
  }

  cout << "Bins in x per cut: " << xbinRanges_bb_ps_e << " on N = " << N << " cuts." << endl;

  // Project the TH2D onto a TH1D for the specified x-axis range
  TH1D* fullProjY_bb_ps_e = hdx_vs_bb_ps_e_data->ProjectionY("_py", xbinLow_bb_ps_e, xbinHigh_bb_ps_e);
  TH1D* wideProjY_bb_ps_e = hdx_vs_bb_ps_e_data->ProjectionY("_wpy", xbinWideLow_bb_ps_e, xbinHigh_bb_ps_e);

  // Use Integral() to get the total number of events within the specified range
  int totalEvents_bb_ps_e = fullProjY_bb_ps_e->Integral();

  std::pair<double,double> wideQualp2;

  auto widep2par_vector = fitAndFineFit(wideProjY_bb_ps_e, "sliceFitWide", "fitFullShift_p2", 7, hcalfit_l, hcalfit_h, wideQualp2, 0, 0, fitopt.c_str());

  double fixed_pshift = widep2par_vector[2].first;
  double fixed_nshift = widep2par_vector[3].first;

  cout << "Fit over full range of PS E yield shifts p:n -> " << fixed_pshift << ":" << fixed_nshift << endl;

  vector<double> cutx_bb_ps_e;

  cout << "Looping over PS E slices, full range.." << endl;
  for (int i = 0; i < N; ++i) {

    int binStart;
    int binEnd;

    if(backwards){
      binStart = xbinBegin;
      binEnd = xbinHigh_bb_ps_e - i*xbinRanges_bb_ps_e;
      cutx_bb_ps_e.push_back(hdx_vs_bb_ps_e_data->GetXaxis()->GetBinLowEdge(binEnd));

    }else{
      binStart = xbinBegin + i*xbinRanges_bb_ps_e;
      binEnd = xbinHigh_bb_ps_e;
      cutx_bb_ps_e.push_back(hdx_vs_bb_ps_e_data->GetXaxis()->GetBinLowEdge(binStart));

    }

    //cout << "Fiducial cut at x=" << cutx_bb_ps_e[i] << endl;

    cout << "Cut from " << hdx_vs_bb_ps_e_data->GetXaxis()->GetBinLowEdge(binStart) << " to " << hdx_vs_bb_ps_e_data->GetXaxis()->GetBinLowEdge(binEnd) << endl;

    //Get the data slices
    TH1D *hdx_slice = hdx_vs_bb_ps_e_data->ProjectionY(Form("hdx_bb_ps_e_slice_%d", i+1), binStart, binEnd);
    TH1D *dx_bb_ps_e_slice = hdx_vs_bb_ps_e_data->ProjectionY(Form("dxfidcut_%d", i+1), binStart, binEnd);

    //Get the MC slices
    TH1D *dx_bb_ps_e_p_slice = hdx_vs_bb_ps_e_p->ProjectionY(Form("dxfidcut_p_%d", i+1), binStart, binEnd);
    TH1D *dx_bb_ps_e_n_slice = hdx_vs_bb_ps_e_n->ProjectionY(Form("dxfidcut_n_%d", i+1), binStart, binEnd);
    hdx_p = hdx_vs_bb_ps_e_p->ProjectionY(Form("hdx_bb_ps_e_slice_p_%d", i+1), binStart, binEnd);
    hdx_n = hdx_vs_bb_ps_e_n->ProjectionY(Form("hdx_bb_ps_e_slice_n_%d", i+1), binStart, binEnd);

    //Get total entries on slice
    int slice_nEntries = dx_bb_ps_e_slice->Integral();
    int slicemcp_nEntries = dx_bb_ps_e_p_slice->Integral();
    int slicemcn_nEntries = dx_bb_ps_e_n_slice->Integral();

    //Rsf extraction using second order poly fit
    std::pair<double,double> sliceQualp2;

    auto slicep2par_vector = fitAndFineFit(hdx_slice, "sliceFit", "fitFullShift_p2", 7, hcalfit_l, hcalfit_h, sliceQualp2, fixed_pshift, fixed_nshift, fitopt.c_str());

    double csndf = sliceQualp2.first/sliceQualp2.second;

    double scaleratio = slicep2par_vector[1].first / slicep2par_vector[0].first;
    double scaleratio_err = sqrt(pow(slicep2par_vector[1].second/slicep2par_vector[1].first, 2) + pow(slicep2par_vector[0].second/slicep2par_vector[0].first, 2)) * scaleratio;

    double xshift_p = slicep2par_vector[2].first;
    double xshift_n = slicep2par_vector[3].first;

    //Write fit results to struct
    bb_ps_ereports[i] = ReportData( dx_bb_ps_e_slice, 
				 dx_bb_ps_e_p_slice,
				 dx_bb_ps_e_n_slice,
				 cutx_bb_ps_e[i],
				 xrange_max_bb_ps_e,
				 slice_nEntries,
				 slicemcp_nEntries,
				 slicemcn_nEntries,
				 scaleratio,
				 scaleratio_err,
				 csndf,
				 "dx vs bb_ps_e" );

    //cout << endl << endl << "Writing out all dx vs fid fit pars: " << endl;
    for (int par=0; par<7; ++par){
      bb_ps_ereports[i].fitParams[par] = slicep2par_vector[par].first;
      bb_ps_ereports[i].fitErrors[par] = slicep2par_vector[par].second;

      //cout << "p" << par << " " << dxreports[i].fitParams[par] << " ";
      
    }
    
    //cout << endl << endl;

  }
 
  //Write out the dx histograms
  
  //get a clone for displays
  TH2D* clonedbb_ps_e = (TH2D*)hdx_vs_bb_ps_e_data->Clone("clonedBb_Ps_E");
  double bb_ps_e_ymin = hcalfit_l;
  double bb_ps_e_ymax = hcalfit_h;

  TCanvas* c0 = new TCanvas("c0", "Slice Locations", 800,600);
  c0->cd();
  clonedbb_ps_e->SetTitle("BBCal Preshower Energy vs dx");
  clonedbb_ps_e->Draw("colz");

  for (auto& cut : cutx_bb_ps_e){
    TLine* line = new TLine(cut, bb_ps_e_ymin, cut, bb_ps_e_ymax);
    line->SetLineWidth(2);
    line->SetLineColor(kRed); // Set line color to red
    line->Draw();
  }
  c0->Update();

  //set up vectors for later tgrapherrors
  std::vector<double> Rsf_vec_bb_ps_e;
  std::vector<double> Rsferr_vec_bb_ps_e;
  std::vector<double> xval_vec_bb_ps_e;
  std::vector<double> nev_vec_bb_ps_e;

  TCanvas* canvasSlices_bb_ps_e = new TCanvas("dx slices over hcal E", "dx slices over hcal E", 1800, 1200);
  
  // Assuming N is the total number of histograms to be plotted on the canvas
  int optimalRows_bb_ps_e = 1;
  int optimalCols_bb_ps_e = 1;
  
  // Find the nearest square root for N to determine the grid size
  int sqrtN_bb_ps_e = (int)std::sqrt(N);
  for (int cols = sqrtN_bb_ps_e; cols <= N; ++cols) {
    if (N % cols == 0) { // If cols is a factor of N
      optimalCols_bb_ps_e = cols;
      optimalRows_bb_ps_e = N / cols;
      break; // Found the optimal layout
    }
  }
  
  canvasSlices_bb_ps_e->Divide(optimalCols_bb_ps_e, optimalRows_bb_ps_e); // Adjust the division based on reportSamples or desired layout
  
  TF1* bgfit_bb_ps_e[N];
  
  // loop over slices
  for (int i = 0; i < N; ++i) {
    if (bb_ps_ereports[i].sliceHistogram == nullptr)
      continue;
    
    canvasSlices_bb_ps_e->cd(i + 1);
    
    //catch bad fit ranges and do not write
    std::string tempTitle = bb_ps_ereports[i].sliceHistogram->GetTitle();      
    if( tempTitle.compare("Null Histogram")==0 )
      continue;
    
    //get mc nucleon parameters
    double p_scale = bb_ps_ereports[i].fitParams[0];
    double n_scale = bb_ps_ereports[i].fitParams[1];
    double p_shift = bb_ps_ereports[i].fitParams[2];
    double n_shift = bb_ps_ereports[i].fitParams[3];
    
    double MCev = bb_ps_ereports[i].mcpnev + bb_ps_ereports[i].mcnnev;
    
    // Simplify the histogram title
    bb_ps_ereports[i].sliceHistogram->SetTitle(Form("%s %0.4f to %0.4f", bb_ps_ereports[0].type.c_str(), bb_ps_ereports[i].winLow, bb_ps_ereports[i].winHigh));
    bb_ps_ereports[i].sliceHistogram->SetLineWidth(1);
    bb_ps_ereports[i].sliceHistogram->Draw();
    
    // Retrieve the number of bins, and the x-axis limits of the slice histogram
    int nbins = bb_ps_ereports[i].sliceHistogram->GetNbinsX();
    double x_low = bb_ps_ereports[i].sliceHistogram->GetXaxis()->GetXmin();
    double x_high = bb_ps_ereports[i].sliceHistogram->GetXaxis()->GetXmax();
    
    //get scale ratio, error, and midpoint for slice
    double ratio = bb_ps_ereports[i].scaleratio;
    double error = bb_ps_ereports[i].scaleratioerr;
    double xval = bb_ps_ereports[i].winLow;
    double nev = (totalEvents_bb_ps_e - bb_ps_ereports[i].nev) / totalEvents_bb_ps_e * 100; //each cut removes data
    
    //fill vectors for later tgrapherrors
    Rsf_vec_bb_ps_e.push_back(ratio);
    Rsferr_vec_bb_ps_e.push_back(error);
    xval_vec_bb_ps_e.push_back(xval);
    nev_vec_bb_ps_e.push_back(nev); //Total remaining stats

    // Draw the corresponding MC histograms
    hdx_p = bb_ps_ereports[i].slicePHistogram; //update the MC slice for the later overall fit
    TH1D *pMCslice = util::shiftHistogramX(bb_ps_ereports[i].slicePHistogram, p_shift);
    pMCslice->Scale(p_scale);
    
    hdx_n = bb_ps_ereports[i].sliceNHistogram; //update the MC slice for the later overall fit
    TH1D *nMCslice = util::shiftHistogramX(bb_ps_ereports[i].sliceNHistogram, n_shift);
    nMCslice->Scale(n_scale);
    
    pMCslice->Draw("same");
    nMCslice->Draw("same");

    //cout << "bin " << i << " pshift " << p_shift << " nshift " << n_shift << endl;
    
    // Create a legend and add the remaining parameters
    TLegend* legend = new TLegend(0.49, 0.59, 0.89, 0.89); // Adjust the position as needed
    legend->SetMargin(0.1);
    legend->AddEntry((TObject*)0, Form("Nev: %0.0f", bb_ps_ereports[i].nev), "");
    legend->AddEntry((TObject*)0, Form("MCev: %0.0f", MCev), "");
    legend->AddEntry((TObject*)0, Form("Ratio: %0.2f", ratio), "");
    legend->AddEntry((TObject*)0, Form("Shift p/n: %0.2f/%0.2f", p_shift, n_shift), "");
    legend->AddEntry((TObject*)0, Form("#chi^{2}/ndf: %0.2f", bb_ps_ereports[i].chisqrndf), "");
    legend->Draw();
    
    bgfit_bb_ps_e[i] = new TF1(Form("bgfit_bb_ps_e_%d",i),fits::g_p2fit_cd,hcalfit_l,hcalfit_h,3);
    //cout << "Background fit parameter" << endl;
    for (int j=0; j<3; ++j){
      bgfit_bb_ps_e[i]->SetParameter(j,bb_ps_ereports[i].fitParams[j+4]);
      //cout << "   " << j << " = " << dxreports[i].fitParams[j+4] << endl;
    }
    
    bgfit_bb_ps_e[i]->SetLineWidth(2);
    bgfit_bb_ps_e[i]->Draw("same");
    
    TF1 *mcfit = new TF1(Form("mcfit_%d", i), fitFullShift_p2, hcalfit_l, hcalfit_h, 7);
    mcfit->SetParameters(bb_ps_ereports[i].fitParams);
    // Can set the error here, might be relevant later
    // mcfit->SetParErrors(dxreports[i].fitErrors);
    
    // Create a new TH1D or fill an existing one with the values from TF1
    TH1D* hFromTF1_bb_ps_e = new TH1D(Form("hFromTF1_bb_ps_e_%d", i), "Histogram from fid x TF1", nbins, x_low, x_high);
    
    for (int bin = 1; bin <= nbins; ++bin) {
      double binCenter = bb_ps_ereports[i].sliceHistogram->GetBinCenter(bin);
      double funcValue = mcfit->Eval(binCenter);
      hFromTF1_bb_ps_e->SetBinContent(bin, funcValue);
    }
    
    int transparentPink = TColor::GetColorTransparent(kPink-4, 0.3);
    hFromTF1_bb_ps_e->SetLineColor(kPink-4);
    hFromTF1_bb_ps_e->SetLineWidth(0);
    hFromTF1_bb_ps_e->SetFillColor(transparentPink);
    hFromTF1_bb_ps_e->SetFillStyle(1001);
    hFromTF1_bb_ps_e->Draw("SAME LF2");
    
    canvasSlices_bb_ps_e->Update();
  }//endloop over N (slices)
  
  canvasSlices_bb_ps_e->Write();
  
  //Create Graphs
  int numValidPoints_bb_ps_e = Rsf_vec_bb_ps_e.size();

  TGraphErrors* graphErrors_bb_ps_e = new TGraphErrors(numValidPoints_bb_ps_e);
  for (int i = 0; i < numValidPoints_bb_ps_e; ++i) {
    graphErrors_bb_ps_e->SetPoint(i, xval_vec_bb_ps_e[i], Rsf_vec_bb_ps_e[i]);
    graphErrors_bb_ps_e->SetPointError(i, 0, Rsferr_vec_bb_ps_e[i]);
  }
  graphErrors_bb_ps_e->SetMarkerStyle(21);
  graphErrors_bb_ps_e->SetMarkerColor(kPink-4);
  graphErrors_bb_ps_e->SetLineColor(kPink-4);
  graphErrors_bb_ps_e->SetTitle("R_{sf} vs. PS E_{cut}; GeV; R_{sf}");

  TGraph* graph_nev = new TGraph(numValidPoints_bb_ps_e);
  for (int i = 0; i < numValidPoints_bb_ps_e; ++i) {
    graph_nev->SetPoint(i, xval_vec_bb_ps_e[i], nev_vec_bb_ps_e[i]);
  }
  graph_nev->SetMarkerStyle(20);
  graph_nev->SetMarkerColor(kGreen-5);
  graph_nev->SetLineColor(kGreen-5);
  graph_nev->SetTitle("Nev Left in Wide Cut vs. PS E_{cut}; GeV; Nev");

  //Create overlayed tgraph object
  gROOT->Reset();
  TCanvas *c2 = new TCanvas("c2","",800,600);
  TPad *pad = new TPad("pad","",0,0,1,1);
  //pad->SetFillColor(42);
  pad->SetGrid();
  pad->Draw();
  pad->cd();

  pad->GetFrame()->SetBorderSize(12);

  graphErrors_bb_ps_e->GetYaxis()->SetLabelSize(0.04);
  graphErrors_bb_ps_e->GetYaxis()->SetLabelFont(42);
  graphErrors_bb_ps_e->GetYaxis()->SetTitleOffset(1.2);
  graphErrors_bb_ps_e->GetYaxis()->SetTitleSize(0.04);
  graphErrors_bb_ps_e->GetYaxis()->SetTitleFont(42); 

  graphErrors_bb_ps_e->Draw("AP");

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

  Double_t xmin = graphErrors_bb_ps_e->GetXaxis()->GetXmin();
  Double_t xmax = graphErrors_bb_ps_e->GetXaxis()->GetXmax();

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
  leg->AddEntry(graphErrors_bb_ps_e, "R_{sf}", "p");
  leg->AddEntry(graph_nev, "N events", "p");
  leg->AddEntry((TObject*)0, Form("ev tot: %d", totalEvents_bb_ps_e), "");
  leg->Draw();

  std::string displayname = "Cuts on dx vs bb_ps_e (all cuts)";
  util::parseAndDisplayCuts(displayname.c_str(), bb_ps_eCuts.c_str());


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
  if(pshift==0&&nshift==0){
    fit->SetParLimits(2,-7,7);
    fit->SetParLimits(3,-7,7);
  }else{
    fit->FixParameter(2,pshift);
    fit->FixParameter(3,nshift);
  }

  //fix bad peak fit for sbs4-50p with tight cuts. Informed from loose cut scales.
  //fit->SetParLimits(0,0.048,1.0);
  
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
  if(pshift==0&&nshift==0){
    fineFit->SetParLimits(2,-7,7);
    fineFit->SetParLimits(3,-7,7);
  }else{
    fineFit->FixParameter(2,pshift);
    fineFit->FixParameter(3,nshift);
    cout << "Fixing fineFit shift parameters at " << pshift << ", " << nshift << endl;
  }

  //fix bad peak fit for sbs4-50p with tight cuts. Informed from loose cut scales.
  //fineFit->SetParLimits(0,0.048,1.0);
  
  histogram->Fit(fineFit, fitOptions.c_str());

  // Update parameters and errors with fine fit results
  for (int i = 0; i < paramCount; ++i) {
    parametersAndErrors[i].first = fineFit->GetParameter(i); // Fine fit parameter value
    parametersAndErrors[i].second = fineFit->GetParError(i); // Fine fit parameter error
  }

  cout << endl << "P0: " << fineFit->GetParameter(0) << endl << endl;
  
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
