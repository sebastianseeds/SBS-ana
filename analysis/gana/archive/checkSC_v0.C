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

//main. kine=kinematic, mag=fieldsetting, pass=pass#, fixshift=fix shift pars for comparisons, bestclus=use best cluster plots, thin=use plot file without correlations, wide=use plots with wide cuts, alt=use plot file using MC alt files
void checkSupercluster(int kine=4, 
		       int mag=30, 
		       int pass=2, 
		       bool fixshift=true,
		       bool bestclus=true, 
		       bool thin=true,
		       bool wide=true,
		       bool alt = true) {
  //set draw params
  gStyle->SetPalette(55);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);


  // reading json configuration file
  JSONManager *jmgr = new JSONManager("../../config/syst.json");

  //Get plot details
  int hbins = jmgr->GetValueFromSubKey<int>( "hbins", Form("sbs%d",kine) ); //gets to 7cm HCal pos res
  double hcalfit_l = jmgr->GetValueFromSubKey<int>( "hcalfit_l", Form("sbs%d",kine) ); //gets to 7cm HCal pos res
  double hcalfit_h = jmgr->GetValueFromSubKey<int>( "hcalfit_h", Form("sbs%d",kine) ); //gets to 7cm HCal pos res
  //double hcalfit_l = econst::hcalposXi_mc; //lower fit/bin limit for hcal dx plots (m)
  //double hcalfit_h = econst::hcalposXf_mc; //upper fit/bin limit for hcal dx plots (m)

  //int Nbins_dx = hbins;

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

  //Get wide/tight elastic cuts
  std::string globalcuts;
  if(wide){
    globalcuts = jmgr->GetValueFromSubKey_str( Form("post_cuts_p%d",pass), Form("sbs%d_%d",kine,mag) );
    cout << "Loaded wide cuts: " << globalcuts << endl;
  }else{
    globalcuts = jmgr->GetValueFromSubKey_str( Form("post_tcuts_p%d",pass), Form("sbs%d_%d",kine,mag) );
    cout << "Loaded tight cuts: " << globalcuts << endl;
  }

  cout << "Loaded tight cuts: " << globalcuts << endl;

  std::vector<std::string> cuts = util::parseCuts(globalcuts); //this makes a vector of all individual cuts in the single globalcut string.

  std::cout << "Parsed cuts: " << std::endl;

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
  for (size_t i = 0; i < cuts.size(); ++i) {

    size_t found_coin = cuts[i].find("coin");
    if (found_coin != std::string::npos) {
      // Found "coin", replace the first '<' with '>'
      size_t found_lt = cuts[i].find('<', found_coin);
      if (found_lt != std::string::npos) {
	std::string modified_cut = cuts[i];
	modified_cut[found_lt] = '>';
	coinanticut_noelas = modified_cut; // Overwrite with this modified cut
	// Append modified cut to coinanticut_elas, but exclude "W2" cuts
	coinanticut_elas += (coinanticut_elas.empty() ? "" : "&&") + modified_cut;
      }
    } else {
      // "coin" not found, just append the cut to coinanticut_elas
      coinanticut_elas += (coinanticut_elas.empty() ? "" : "&&") + cuts[i];
    }
  }

  std::cout << "coinanticut_noelas: " << coinanticut_noelas << std::endl;
  std::cout << "coinanticut_elas: " << coinanticut_elas << std::endl;

  //return;

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
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fin_path = outdir_path + Form("/parse/parse_sbs%d_pass%d_barebones.root",kine,pass);

  std::string bestclus_word = "";
  if(bestclus)
    bestclus_word = "_bestclus";

  std::string skipcorrelations_word = "";
  if(skipcorrelations)
    skipcorrelations_word = "_thin";

  std::string wide_word = "";
  if(wide)
    wide_word = "_widecut";

  std::string fout_path = outdir_path + Form("/gmn_analysis/dx_correlations_sbs%d_mag%d_pass%d%s%s%s.root",kine,mag,pass,bestclus_word.c_str(),skipcorrelations_word.c_str(),wide_word.c_str());

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
  TH1D* hist_base = new TH1D( histBaseName.c_str(), "dx;m", hbins, hcalfit_l, hcalfit_h );

  // Draw the plot using the created histogram
  tree->Draw(("dx>>" + histBaseName).c_str(), globalcuts.c_str(), "COLZ");

  // Write the histogram to the output file
  hist_base->Write();

  delete hist_base;

  std::string histAntidyName = "hdx_dyanti";
  TH1D* hist_antidy = new TH1D( histAntidyName.c_str(), "dx, dy anticut;m", hbins, hcalfit_l, hcalfit_h );

  // Draw the plot using the created histogram
  tree->Draw(("dx>>" + histAntidyName).c_str(), dyanticut_elas.c_str(), "COLZ");

  // Write the histogram to the output file
  hist_antidy->Write();


















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
    bestclus_word = "_bestclus";

  std::string thin_word = "";
  if(thin)
    thin_word = "_thin";

  std::string alt_word = "";
  if(alt)
    alt_word = "_alt";

  std::string wide_word = "";
  if(wide)
    wide_word = "_widecut";

  std::string basePath = "/lustre19/expphy/volatile/halla/sbs/seeds";
  std::string finPath = Form("%s/gmn_analysis/dx_correlations_sbs%d_mag%d_pass%d%s%s%s.root", basePath.c_str(), kine, mag, pass, bestclus_word.c_str(), thin_word.c_str(), wide_word.c_str());
  std::string fmcinPath = Form("%s/gmn_analysis/dx_mc_sbs%d_mag%d_pass%d%s%s%s.root", basePath.c_str(), kine, mag, pass, thin_word.c_str(), alt_word.c_str(), wide_word.c_str());
  std::string foutPath = Form("%s/gmn_analysis/fidstab_gmn_sbs%d_mag%d_pass%d%s%s.root", basePath.c_str(), kine, mag, pass, bestclus_word.c_str(), wide_word.c_str());

  TFile* inputFile = new TFile(finPath.c_str());
  if (!inputFile || inputFile->IsZombie()) 
    handleError(inputFile,"inputFile");

  //Get all histograms
  TH1D *hdx_data = dynamic_cast<TH1D*>(inputFile->Get("hdx_allcut"));
  hdx_data->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  TH2D *hdx_vs_fidx_data = dynamic_cast<TH2D*>(inputFile->Get("hist_fiducial_sig_x"));
  hdx_vs_fidx_data->GetYaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  int hdx_vs_fidx_data_N = hdx_vs_fidx_data->GetEntries();
  std::string fidxCuts = hdx_vs_fidx_data->GetTitle();
  cout << endl << "Opened dx vs fid x with cuts: " << fidxCuts << endl << endl;

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

  TH2D *hdx_vs_fidx_n = dynamic_cast<TH2D*>(inputFileMC->Get("hist_fiducial_sig_x_n"));
  hdx_vs_fidx_n->GetYaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  int hdx_vs_fidx_n_N = hdx_vs_fidx_n->GetEntries();

  TH2D *hdx_vs_fidy_n = dynamic_cast<TH2D*>(inputFileMC->Get("hist_fiducial_sig_y_n"));
  hdx_vs_fidy_n->GetYaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  int hdx_vs_fidy_n_N = hdx_vs_fidy_n->GetEntries();

  TH2D *hdx_vs_fidx_p = dynamic_cast<TH2D*>(inputFileMC->Get("hist_fiducial_sig_x_p"));
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
  int xbinCutMax_fidx = hdx_vs_fidx_data->GetXaxis()->FindBin(sliceCutMax);

  cout << sliceCutMax << " " << xbinCutMax_fidx << " " << xrange_max_fidx << " " << xbinHigh_fidx << endl;

  // Bins per range
  int xbinRanges_fidx = (xbinHigh_fidx - xbinLow_fidx) / N;
  if(sliceCutMax>0)
    xbinRanges_fidx = (xbinCutMax_fidx - xbinLow_fidx) / N;

  cout << "Bins in x per cut: " << xbinRanges_fidx << endl;

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

    int binStart = xbinLow_fidx + i*xbinRanges_fidx;
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
    int slice_nEntries = dx_fidx_slice->Integral(binStart,binEnd);
    int slicemcp_nEntries = dx_fidx_p_slice->Integral(binStart,binEnd);
    int slicemcn_nEntries = dx_fidx_n_slice->Integral(binStart,binEnd);

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
    line->SetLineColor(kRed); // Set line color to red
    line->Draw();
  }
  c0->Update();

  //set up vectors for later tgrapherrors
  std::vector<double> Rsf_vec_fidx;
  std::vector<double> Rsferr_vec_fidx;
  std::vector<double> xval_vec_fidx;
  
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
    
    //fill vectors for later tgrapherrors
    Rsf_vec_fidx.push_back(ratio);
    Rsferr_vec_fidx.push_back(error);
    xval_vec_fidx.push_back(xval);
    
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
    
    bgfit_fidx[i] = new TF1(Form("bgfit_fidx_%d",i),fits::g_p2fit,hcalfit_l,hcalfit_h,3);
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
    
    int transparentOrange = TColor::GetColorTransparent(kOrange, 0.3);
    hFromTF1_fidx->SetLineColor(kOrange);
    hFromTF1_fidx->SetLineWidth(0);
    hFromTF1_fidx->SetFillColor(transparentOrange);
    hFromTF1_fidx->SetFillStyle(1001);
    hFromTF1_fidx->Draw("SAME LF2");
    
    canvasSlices_fidx->Update();
  }//endloop over N (slices)
  
  canvasSlices_fidx->Write();
  
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
