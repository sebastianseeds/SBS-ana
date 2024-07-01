//sseeds 10.23.23 - Script to run over pgun and ngun digitized and replayed data and return histograms for further analysis

#include <vector>
#include <iostream>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TCanvas.h"
#include "TF2.h"
#include "TH1.h"
#include "TVirtualFitter.h"
#include "TRandom.h"
#include "TFitResultPtr.h"
#include "TMinuit.h"
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

const int nbin = 150;
const double xmin1 = 1.4;
const double xmin2 = 9.0;

void FitAndDrawGraphWithErrorBand(TGraph* graph1, 
				  TGraph* graph2, 
				  TCanvas* canvas, 
				  double fitMin, 
				  double fitMax, 
				  const vector<double>& N1_at_bin, 
				  const vector<double>& N2_at_bin, 
				  double tfac,
				  vector<double> pN,
				  vector<double> hde,
				  vector<double> hdeerr,
				  TFile *f1) {
  
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetEndErrorSize(0);

  // Create and set up the fit function
  TF1* fitFunc = new TF1("fitFunc", "pol4", fitMin, fitMax);
    
  // Fit the first graph within the specified range
  graph1->Fit(fitFunc, "R");

  cout << endl << endl;

  cout << "Value of proton MC fit at pN = 2.36 GeV (SBS4): " << fitFunc->Eval(2.36) << endl;
  cout << "Value of proton MC fit at pN = 3.17 GeV (SBS8): " << fitFunc->Eval(3.17) << endl;
  cout << "Value of proton MC fit at pN = 3.20 GeV (SBS9): " << fitFunc->Eval(3.20) << endl;

  cout << endl;

  // Create a TGraphErrors to represent the error band for the first graph
  TGraphErrors* errorBand1 = new TGraphErrors();
  errorBand1->SetName("errorBand1");
    
  int N1_min = 100000;
  int N2_min = 100000;

  int nPoints1 = graph1->GetN();
  for (int i = 0; i < nPoints1; ++i) {
    double x = graph1->GetX()[i];
    double fitValue = fitFunc->Eval(x);
    double error = TMath::Sqrt(fitValue/100 * (1 - fitValue/100) / N1_at_bin[i])*100.;
    errorBand1->SetPoint(i, x, fitValue);
    errorBand1->SetPointError(i, 0, error);

    if( N1_at_bin[i]<N1_min )
      N1_min=N1_at_bin[i];
    
    cout << "proton (x,fit,err,N): " << x << " " << fitValue << " " << error << " " << N1_at_bin[i] << endl;
  }

  // Fit the second graph within the specified range
  graph2->Fit(fitFunc, "R+");

  // Create a TGraphErrors to represent the error band for the second graph
  TGraphErrors* errorBand2 = new TGraphErrors();
  errorBand2->SetName("errorBand2");
    
  int nPoints2 = graph2->GetN();
  for (int i = 0; i < nPoints2; ++i) {
    double x = graph2->GetX()[i];
    double fitValue = fitFunc->Eval(x);
    double error = TMath::Sqrt(fitValue/100 * (1 - fitValue/100) / N2_at_bin[i])*100.;
    errorBand2->SetPoint(i, x, fitValue);
    errorBand2->SetPointError(i, 0, error);

    if( N2_at_bin[i]<N2_min )
      N2_min=N2_at_bin[i];

    //cout << "neutron: " << fitValue << " " << error << " " << N2_at_bin[i] << endl;
  }

  // Set fill colors for the error bands
  errorBand1->SetFillStyle(1001);
  errorBand1->SetFillColorAlpha(kRed, 0.35);
  errorBand1->SetLineWidth(1);
  errorBand2->SetFillStyle(1001);
  errorBand2->SetFillColorAlpha(kBlue, 0.35);
  errorBand2->SetLineWidth(1);

  //Add data points
  //SBS4 dx sideband
  auto graph3 = new TGraphErrors();
  graph3->SetPoint(0,pN[0],hde[0]);
  graph3->SetPointError(0,0,hdeerr[0]);
  graph3->SetMarkerColor(kGreen-3);
  graph3->SetMarkerStyle(20);
  graph3->SetMarkerSize(1);
  graph3->SetLineColor(kGreen-3);
  graph3->SetLineWidth(2);
  //SBS4 W2 anticut
  auto graph4 = new TGraphErrors();
  graph4->SetPoint(0,pN[1],hde[1]);
  graph4->SetPointError(0,0,hdeerr[1]);
  graph4->SetMarkerColor(kGreen-3);
  graph4->SetMarkerStyle(21);
  graph4->SetMarkerSize(1);
  graph4->SetLineColor(kGreen-3);
  graph4->SetLineWidth(2);

  //SBS8 dx sideband
  auto graph5 = new TGraphErrors();
  graph5->SetPoint(0,pN[2],hde[2]);
  graph5->SetPointError(0,0,hdeerr[2]);
  graph5->SetMarkerColor(kGreen-1);
  graph5->SetMarkerStyle(20);
  graph5->SetMarkerSize(1);
  graph5->SetLineColor(kGreen-1);
  graph5->SetLineWidth(2);

  //SBS8 W2 anticut
  auto graph6 = new TGraphErrors();
  graph6->SetPoint(0,pN[3],hde[3]);
  graph6->SetPointError(0,0,hdeerr[3]);
  graph6->SetMarkerColor(kGreen-1);
  graph6->SetMarkerStyle(21);
  graph6->SetMarkerSize(1);
  graph6->SetLineColor(kGreen-1);
  graph6->SetLineWidth(2);

  //SBS4 30% positional spot
  auto graph7 = new TGraphErrors();
  graph7->SetPoint(0,pN[4],hde[4]);
  graph7->SetPointError(0,0,hdeerr[4]);
  graph7->SetMarkerColor(kBlack);
  graph7->SetMarkerStyle(71);
  graph7->SetMarkerSize(2);
  graph7->SetLineColor(kBlack);
  //graph7->SetLineWidth(2);

  //SBS8 70% positional spot
  auto graph8 = new TGraphErrors();
  graph8->SetPoint(0,pN[5],hde[5]);
  graph8->SetPointError(0,0,hdeerr[5]);
  graph8->SetMarkerColor(kRed);
  graph8->SetMarkerStyle(21);
  graph8->SetMarkerSize(2);
  graph8->SetLineColor(kRed);
  //graph8->SetLineWidth(2);

  //SBS9 70% positional spot
  auto graph9 = new TGraphErrors();
  graph9->SetPoint(0,pN[6],hde[6]);
  graph9->SetPointError(0,0,hdeerr[6]);
  graph9->SetMarkerColor(kRed);
  graph9->SetMarkerStyle(22);
  graph9->SetMarkerSize(2);
  graph9->SetLineColor(kRed);
  //graph9->SetLineWidth(2);

  //SBS8 0% positional spot, dip correction
  auto graph10 = new TGraphErrors();
  graph10->SetPoint(0,pN[7],hde[7]);
  graph10->SetPointError(0,0,hdeerr[7]);
  graph10->SetMarkerColor(kBlack);
  graph10->SetMarkerStyle(33);
  graph10->SetMarkerSize(2);
  graph10->SetLineColor(kBlack);
  //graph10->SetLineWidth(2);

  //SBS8 50% positional spot, dip correction
  auto graph11 = new TGraphErrors();
  graph11->SetPoint(0,pN[8],hde[8]);
  graph11->SetPointError(0,0,hdeerr[8]);
  graph11->SetMarkerColor(kBlack);
  graph11->SetMarkerStyle(34);
  graph11->SetMarkerSize(2);
  graph11->SetLineColor(kBlack);
  //graph10->SetLineWidth(2);

  //SBS8 70% positional spot, dip correction
  auto graph12 = new TGraphErrors();
  graph12->SetPoint(0,pN[9],hde[9]);
  graph12->SetPointError(0,0,hdeerr[9]);
  graph12->SetMarkerColor(kBlack);
  graph12->SetMarkerStyle(29);
  graph12->SetMarkerSize(2);
  graph12->SetLineColor(kBlack);
  //graph10->SetLineWidth(2);

  //SBS11 70% positional spot, dip correction
  auto graph13 = new TGraphErrors();
  graph13->SetPoint(0,pN[10],hde[10]);
  graph13->SetPointError(0,0,hdeerr[10]);
  graph13->SetMarkerColor(kBlack);
  graph13->SetMarkerStyle(22);
  graph13->SetMarkerSize(2);
  graph13->SetLineColor(kBlack);
  //graph11->SetLineWidth(2);

  TMultiGraph* mg = new TMultiGraph();
  mg->SetTitle("HCAL efficiency, MC/data Comparison");
  mg->Add(graph1, "AP");
  mg->Add(graph2, "AP");
  //mg->Add(graph3, "AP");
  //mg->Add(graph4, "AP");
  //mg->Add(graph5, "AP");
  //mg->Add(graph6, "AP");
   mg->Add(graph7, "AP");
  // mg->Add(graph8, "AP");
  // mg->Add(graph9, "AP");
  mg->Add(graph10, "AP");
  mg->Add(graph11, "AP");
  mg->Add(graph12, "AP");
  // mg->Add(graph13, "AP");

  mg->Add(errorBand1, "E3");
  mg->Add(errorBand2, "E3");

  //Add a legend to the canvas
  //auto l1 = new TLegend(0.11,0.65,0.45,0.89);
  auto l1 = new TLegend(0.51,0.11,0.89,0.49);
  //l1->SetTextSize( 0.03 );
  l1->AddEntry( graph1, Form("Proton, min ev/cell: %d",N1_min), "p");
  l1->AddEntry( graph2, Form("Neutron, min ev/cell: %d",N2_min), "p");
  //l1->AddEntry( graph3, "Data: LH2, SBS4, dx sideband", "p");
  //l1->AddEntry( graph4, "Data: LH2, SBS4, W^{2} anticut", "p");
  //l1->AddEntry( graph5, "Data: LH2, SBS8, dx sideband", "p");
  //l1->AddEntry( graph6, "Data: LH2, SBS8, W^{2} anticut", "p");
  //l1->AddEntry( graph4, "Data: LH2, SBS8, Anticut Method", "p");
   l1->AddEntry( graph7, "Data: LH2, SBS4 30%, proton", "p");
  // l1->AddEntry( graph8, "Data: LH2, SBS8 70%, proton", "p");
  // l1->AddEntry( graph9, "Data: LH2, SBS9 70%, proton", "p");
  l1->AddEntry( graph10, "Data: LH2, SBS8 0%, proton, dip correction", "p");
  l1->AddEntry( graph11, "Data: LH2, SBS8 50%, proton, dip correction", "p");
  l1->AddEntry( graph12, "Data: LH2, SBS8 70%, proton, dip correction", "p");
  // l1->AddEntry( graph13, "Data: LH2, SBS9 70%, proton, dip correction", "p");
  l1->AddEntry( (TObject*)0, "", "");
  //l1->AddEntry( (TObject*)0, "Binomial Error on MC Fits", "");
  l1->AddEntry( (TObject*)0, Form("Threshold Energy E_{T} = %0.2f*E_{Peak}",1/(double)tfac), "");
  l1->AddEntry( (TObject*)0, "4x4 clusters", "");
  //l1->AddEntry( errorBand1, "Total fit", "l");
  //l1->AddEntry( errorBand2, "Total fit", "l");

  // Draw the TMultiGraph
  mg->Draw("AP3");

  l1->Draw("same");

  mg->GetYaxis()->SetTitle("Efficiency (%)");
  mg->GetYaxis()->SetRangeUser(80,100);
  mg->GetXaxis()->SetLimits(1.,8.);

  //mg->GetXaxis()->SetRangeUser(0.,10.);
  //mg->GetXaxis()->SetRange(0.,10.);
  // mg->SetMaximum(10);
  // mg->SetMinimum(0);
  mg->GetXaxis()->SetTitle("Nucleon Momentum (GeV/c)");

  // Update the canvas
  canvas->Update();
  canvas->Modified();
  gPad->Modified();
  //mg->GetXaxis()->SetRange(0.,10.);

  canvas->SaveAs("/work/halla/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/quality_plots/hcal_de_new.png");
  canvas->SaveAs("/work/halla/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/quality_plots/hcal_de_new.pdf");
  canvas->SaveAs("/work/halla/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/quality_plots/hcal_de_new.jpeg");

  canvas->Write();
}

void FitAndDrawGraph(TGraph* graph, TCanvas* canvas, double fitMin, double fitMax) {
    // Create and set up the fit function
    TF1* fitFunc = new TF1("fitFunc", "pol4", fitMin, fitMax);
    
    // Fit the graph within the specified range
    graph->Fit(fitFunc, "QR");

    // Draw the graph and the fit curve
    graph->Draw("AP");
    fitFunc->Draw("SAME");

    // Update the canvas
    canvas->Update();
}


//Uses g4sbs replays of simulated data set containing pgun/ngun, zero field, SBS4 geometry
void hde_mc_data( int iter = 1 ) //iteration 0 gets mean values of hcalE vs nucleonp; 1 gets eff
{ //main  
  
  if( iter!=0 && iter!=1 ){
    cerr << "ERROR: Enter only iteration 0 (get E mean) or iteration 1 (get detection efficiency)" << endl;
    return;
  }
    
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  std::string date = util::getDate();

  //gStyle->SetEndErrorSize(0);

  //set draw params
  gStyle->SetPalette(55);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  //gStyle->SetEndErrorSize(0);

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/smchde.json");
  
  string proton_rfile = jmgr->GetValueFromSubKey_str( "rootfiles", "proton" );
  string neutron_rfile = jmgr->GetValueFromSubKey_str( "rootfiles", "neutron" );
  string prot_param_path = jmgr->GetValueFromSubKey_str( "params", "proton" );
  string neut_param_path = jmgr->GetValueFromSubKey_str( "params", "neutron" );
  string outfilename;
  if(iter==0)
    outfilename = jmgr->GetValueFromKey_str( "outfile_0" );
  else
    outfilename = jmgr->GetValueFromKey_str( "outfile_1" );

  double pmin = jmgr->GetValueFromKey<double>( "pmin" );
  double pmax = jmgr->GetValueFromKey<double>( "pmax" );
  double Emin = jmgr->GetValueFromKey<double>( "Emin" );
  double Emax = jmgr->GetValueFromKey<double>( "Emax" );
  double fit_fac = jmgr->GetValueFromKey<double>( "fit_fac" );
  double tfac = jmgr->GetValueFromKey<double>( "tfac" );
  int jmgr_nbin = jmgr->GetValueFromKey<int>( "nbin" );
  // double pN_cent = jmgr->GetValueFromSubKey<Double_t>( "pN_cent", Form("sbs%d",kine) );
  // double HDE_sb = jmgr->GetValueFromSubKey<Double_t>( "HDE_sb", Form("sbs%d",kine) );
  // double HDE_sb_err = jmgr->GetValueFromSubKey<Double_t>( "HDE_sb_err", Form("sbs%d",kine) );
  // double HDE_anticut = jmgr->GetValueFromSubKey<Double_t>( "HDE_anticut", Form("sbs%d",kine) );
  // double HDE_anticut_err = jmgr->GetValueFromSubKey<Double_t>( "HDE_anticut_err", Form("sbs%d",kine) );
  vector<double> pN; jmgr->GetVectorFromKey<double>( "pN", pN );
  vector<double> hde; jmgr->GetVectorFromKey<double>( "hde", hde );
  vector<double> hdeerr; jmgr->GetVectorFromKey<double>( "hdeerr", hdeerr );
  vector<double> hdeerr_b; jmgr->GetVectorFromKey<double>( "hdeerr_b", hdeerr_b );

  if( jmgr_nbin!=nbin ){
    cout << "ERROR: const int nbin not equal to common configuration json file value." << endl;
    return;
  }

  //general static parameters for this analysis
  double p_step = (pmax-pmin)/nbin; //Amount of momentum traversed per bin
  double oEmean_p[nbin] = {0.}; //hcal proton E mean values obtained from iter 0
  double oEmean_n[nbin] = {0.}; //hcal neutron E mean values obtained from iter 0
  
  if( iter==1 ){
    //Get proton E means
    ifstream prot_param_file( prot_param_path );
    if( !prot_param_file ){
      cerr << endl << "ERROR: No input constant file present -> path to meanE_proton.txt expected." << endl;
      return;
    }
  
    Int_t n1=0;
    Double_t d1;
    string line1;
    
    while( getline( prot_param_file, line1 ) ){
      if( line1.at(0) == '#' ) {
  	continue;
      }
      stringstream ss( line1 );
      ss >> d1;     
      oEmean_p[n1] = d1;      
      n1++;
    }

    //Get neutron E means
    ifstream neut_param_file( neut_param_path );
    if( !neut_param_file ){
      cerr << endl << "ERROR: No input constant file present -> path to meanE_proton.txt expected." << endl;
      return;
    }
  
    n1=0;
    d1=0;
    string line2;
    
    while( getline( neut_param_file, line2 ) ){
      if( line2.at(0) == '#' ) {
  	continue;
      }
      stringstream ss( line2 );
      ss >> d1;     
      oEmean_n[n1] = d1;      
      n1++;
    }

    //console out to check
    cout << endl << "proton E vs p means: " << endl;
    for( int i=0; i<nbin; ++i )
      cout << oEmean_p[i] << endl;
    
    cout << endl << endl << "neutron E vs p means: " << endl;
    for( int i=0; i<nbin; ++i )
      cout << oEmean_n[i] << endl;
  }

  //set up output files
  //TFile *fout = new TFile( outfilename.c_str(), "RECREATE" );
  TFile *fout = new TFile( "test.root", "RECREATE" );
  
  //set up diagnostic histograms
  TH2D *hEdepvP_p = new TH2D("hEdepvP_p","HCal E dep vs proton momentum; p_{proton} (GeV); E_{hcal} (GeV)", nbin, 1, 9, 200, Emin, Emax);
  TH2D *hEdepvP_n = new TH2D("hEdepvP_n","HCal E dep vs neutron momentum; p_{neutron} (GeV); E_{hcal} (GeV)", nbin, 1, 9, 200, Emin, Emax);
  TH2D *hEdepvP_p_Ecut = new TH2D("hEdepvP_p_Ecut","HCal E dep vs proton momentum, Mean E / 4 cut; p_{proton} (GeV); E_{hcal} (GeV)", nbin, 1, 9, 200, Emin, Emax);
  TH2D *hEdepvP_n_Ecut = new TH2D("hEdepvP_n_Ecut","HCal E dep vs neutron momentum, Mean E / 4 cut; p_{neutron} (GeV); E_{hcal} (GeV)", nbin, 1, 9, 200, Emin, Emax);
  TH1D *hdx_p = new TH1D("hdx_p","dx proton (sd track);x_{HCAL}-x_{expect} (m)", 400, -2, 2);
  TH1D *hdy_p = new TH1D("hdy_p","dy proton (sd track);y_{HCAL}-y_{expect} (m)", 400, 3.8, 7.8);
  TH1D *hdx_n = new TH1D("hdx_n","dx neutron (sd track);x_{HCAL}-x_{expect} (m)", 400, -2, 2);
  TH1D *hdy_n = new TH1D("hdy_n","dy neutron (sd track);y_{HCAL}-y_{expect} (m)", 400, 3.8, 7.8);
  TH1D *hdx_p_v2 = new TH1D("hdx_p_v2","dx proton (angles);x_{HCAL}-x_{expect} (m)", 400, -2, 2);
  TH1D *hdy_p_v2 = new TH1D("hdy_p_v2","dy proton (angles);y_{HCAL}-y_{expect} (m)", 400, -2, 2);
  TH1D *hdx_n_v2 = new TH1D("hdx_n_v2","dx neutron (angles);x_{HCAL}-x_{expect} (m)", 400, -2, 2);
  TH1D *hdy_n_v2 = new TH1D("hdy_n_v2","dy neutron (angles);y_{HCAL}-y_{expect} (m)", 400, -2, 2);
  TH1D *hxexp = new TH1D("hxexp","x exp (angles);x_{expect} (m)", 400, -2, 2);
  TH1D *hyexp = new TH1D("hyexp","y exp (angles);y_{expect} (m)", 400, -2, 2);
  
  TH2D *hdxvp_p = new TH2D("hdxvp_p","dx vs proton p; p_{p} (GeV); x_{HCAL}-x_{expect} (m)", nbin, 1, 9, 400, -2, 2);
  TH2D *hdyvp_p = new TH2D("hdyvp_p","dy vs proton p; p_{p} (GeV); y_{HCAL}-y_{expect} (m)", nbin, 1, 9, 400, -2, 2);
  TH2D *hdxvp_n = new TH2D("hdxvp_n","dx vs neutron p; p_{n} (GeV); x_{HCAL}-x_{expect} (m)", nbin, 1, 9, 400, -2, 2);
  TH2D *hdyvp_n = new TH2D("hdyvp_n","dy vs neutron p; p_{n} (GeV); y_{HCAL}-y_{expect} (m)", nbin, 1, 9, 400, -2, 2);

  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;
  std::string nuc;
  
  // declare ints to hold number of nucleons that pass/fail energy check on iter 1
  static int P_pass[nbin] = {0};
  static int P_tot[nbin] = {0};
  static int N_pass[nbin] = {0};
  static int N_tot[nbin] = {0};

  //loop over nucleons
  for (Int_t n=0; n<2; n++) {
    
    C = new TChain("T");
    
    if( n==0 ){ 
      nuc = "proton";
      C->Add(proton_rfile.c_str());
    }else if( n==1 ){ 
      nuc = "neutron";
      C->Add(neutron_rfile.c_str());
    }else
      break;
    
    if( C->GetEntries()==0 ){
      cout << "ERROR: Missing MC file or empty MC file." << endl;
      return;
    } 

    // setting up ROOT tree branch addresses
    C->SetBranchStatus("*",0);    
    
    // HCal general
    Double_t hcalid, hcale, hcalx, hcaly, hcalr, hcalc, hcaltdc, hcalatime;
    std::vector<std::string> hcalvar = {"idblk","e","x","y","rowblk","colblk","tdctimeblk","atimeblk"};
    std::vector<void*> hcalvarlink = {&hcalid,&hcale,&hcalx,&hcaly,&hcalr,&hcalc,&hcaltdc,&hcalatime};
    rvars::setbranch(C, "sbs.hcal", hcalvar, hcalvarlink);
    
    // MC nucleon
    Double_t mc_p, mc_px, mc_py, mc_pz, mc_vx, mc_vy, mc_vz, mc_nucl, mc_posx, mc_posy;
    std::vector<std::string> mcvar = {"mc_ep","mc_epx","mc_epy","mc_epz","mc_vx","mc_vy","mc_vz","mc_nucl","sdtrack_posx","sdtrack_posy"};
    std::vector<void*> mcvarlink = {&mc_p,&mc_px,&mc_py,&mc_pz,&mc_vx,&mc_vy,&mc_vz,&mc_nucl,&mc_posx,&mc_posy};
    rvars::setbranch(C, "MC", mcvar, mcvarlink);
    
    // event indices
    long nevent = 0, nevents = C->GetEntries(); 

    while (C->GetEntry(nevent++)) {
      
      std::cout << "Processing " << nuc << " MC data, event " << nevent << "/" << nevents << "\r";
      std::cout.flush();
      
      Int_t bin = (mc_p-pmin)/p_step; //will fail if MC isn't configured to match pmin and pmax
      Double_t E_thresh;

      if( bin>=nbin || bin < 0 ){
	std::cout << "Warning: nucleon momentum bin out of bounds at " << bin << endl;
	continue;
      }

      if( n==0 ){
	P_tot[bin]++;
	E_thresh = oEmean_p[bin]/tfac;
	if( hcale > E_thresh ){
	  P_pass[bin]++;
	  hEdepvP_p_Ecut->Fill(mc_p,hcale);
	}
	hEdepvP_p->Fill(mc_p,hcale);
      }
      if( n==1 ){
	N_tot[bin]++;
	E_thresh = oEmean_n[bin]/tfac;
	if( hcale > E_thresh ){
	  N_pass[bin]++;
	  hEdepvP_n_Ecut->Fill(mc_p,hcale);
	}
	hEdepvP_n->Fill(mc_p,hcale);
      }

    } //end event loop

    // reset chain for the next run config
    C->Reset();

  } //end run loop
  
  if(iter==0){
    //finish up
    fout->Write();

    cout << "Analysis complete. Output file written to " << outfilename << endl;

    // Send time efficiency report to console
    st->Stop();
    cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;
    
    return;
  }

  //get a couple of vectors for later
  vector<double> Nval_p;
  vector<double> Nval_n;

  for( int i=0; i<nbin; ++i ){
    Nval_p.push_back((double)P_tot[i]);
    Nval_n.push_back((double)N_tot[i]);
  }

  //Reporting on iteration 1
  gStyle->SetOptFit();
  //gStyle->SetEndErrorSize(0);
  TCanvas *c1 = new TCanvas("c1","HCal E vs Nucleon P",1600,1200);
  c1->SetGrid();
  c1->cd();

  //Iter 1 arrays
  double binp[nbin] = {0.};
  double hde_proton[nbin] = {0.};
  double hde_neutron[nbin] = {0.};
  double hde_npratio[nbin] = {0.};

  auto mg = new TMultiGraph();

  //loop over bins and extract ratio of detected / expected per nucleon and super ratio
  for(int b=0; b<nbin; ++b){
    double p = b*p_step+pmin;
    binp[b] = p;
      
    hde_proton[b] = (double)P_pass[b] / (double)P_tot[b] *100.;
    hde_neutron[b] = (double)N_pass[b] / (double)N_tot[b] *100.;
      
    hde_npratio[b] = hde_neutron[b]/hde_proton[b];

  }
    
  //Draw simple graphs
  auto grp = new TGraph(nbin,binp,hde_proton);
  grp->SetTitle("Proton");
  grp->SetMarkerColor(kRed);
  grp->SetMarkerStyle(20);
  grp->SetMarkerSize(0.5);
  grp->SetLineColor(kRed);
  grp->SetLineWidth(0);
  mg->Add(grp);

  //auto grn = new TGraphErrors(nbin,binp,hde_neutron,binerr_n,binerr_n);
  auto grn = new TGraph(nbin,binp,hde_neutron);
  grn->SetTitle("Neutron");
  grn->SetMarkerColor(kBlue);
  grn->SetMarkerStyle(21);
  grn->SetMarkerSize(0.5);
  grn->SetLineColor(kBlue);
  grn->SetLineWidth(0);
  mg->Add(grn);

  mg->SetTitle(Form("HCAL Efficiency (E_{T}=1/%0.0f E_{Peak}) (4x4 cluster)",tfac));
  mg->GetXaxis()->SetTitle("Nucleon Momentum (GeV/c)");
  mg->GetYaxis()->SetTitle("Efficiency (%)");
  mg->Draw("AP");

  mg->GetYaxis()->SetRangeUser(80.,105.);

  c1->Modified();

  c1->BuildLegend();

  c1->Write();

  TCanvas *c2 = new TCanvas("c2","HCal Efficiency Ratio (N/P)(E_{T}=1/4 E_{peak})",1600,1200);
  c2->SetGrid();
  c2->cd();

  auto grr = new TGraph(nbin,binp,hde_npratio);
  grr->SetTitle("HCal Efficiency Ratio (N/P)(E_{T}=1/4 E_{peak})");
  grr->SetMarkerColor(kBlack);
  grr->SetMarkerStyle(20);
  grr->SetMarkerSize(2);
  grr->SetLineColor(kMagenta);
  grr->SetLineWidth(0);
  grr->GetXaxis()->SetTitle("Nucleon Momentum (GeV/c)");
  grr->GetYaxis()->SetTitle("Efficiency Ratio (N/P)");
  grr->Draw("AP");

  grr->GetYaxis()->SetRangeUser(0.9,1.05);

  c2->Modified();

  c2->Write();

  TCanvas *c3 = new TCanvas("c3","HCal Efficiency Ratio (N/P)(E_{T}=1/4 E_{peak}) w/Error",1600,1200);
  c3->SetGrid();
  c3->cd();

  //FitTwoGraphsWithPolynomial(grp,grn,c3);
  //FitAndDrawGraph(grp,c3,xmin1,xmin2);
  //FitAndDrawGraphWithErrorBand(grp,grn,c3,xmin1,xmin2,Nval_p,Nval_n,tfac,fout);
  FitAndDrawGraphWithErrorBand(grp,grn,c3,xmin1,xmin2,Nval_p,Nval_n,tfac,pN,hde,hdeerr_b,fout);

  c3->Modified();
  c3->Write();

  //finish up
  fout->Write();

  cout << "Analysis complete. Output file written to " << outfilename << endl;

  // Send time efficiency report to console
  st->Stop();
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
