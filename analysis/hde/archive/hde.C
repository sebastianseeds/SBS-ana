//SSeeds 11.03.21 - Production - Code designed to project scattered proton during zero-field running to face of HCal and determine detection efficiency. Currently measures relative to BB elastic track detection and sets a limit on good hits in HCal at 2*sigma the proton peak at zero field.
//Updates 1.11.23 - Analysis - Updates to use gmna framework and to assess various methods for obtaining detection efficiency from data for HCal
//Update 2.23.23 - Corrections and improvements to reporting and fits.
//Update 3.9.23 - Added small changes to allow for efficiency over all field settings. Added various methods for extracting efficiencies.
//Update 3.23.23 - Generalized fits with header gmna.h and added efficiency updates
//Update 3.28.23 - Updates for analysis class upgrade and jupyter notebook applications

#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TString.h"
#include <iostream>
#include <fstream>
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include "TMath.h"
#include "../include/gmn.h"

//Configure total fits for this application. Component fits in gmna.h
Double_t Tfit_hcal(Double_t *x, Double_t *par){ //4 order poly fit for HCal bg
  return g_p4fit(x,&par[0])+g_gfit(x,&par[5]);
}
Double_t Tfit_expo(Double_t *x, Double_t *par){
  return g_expofit(x,&par[0])+g_gfit(x,&par[3]);  //expo fit
}
Double_t Tfit_scexpo(Double_t *x, Double_t *par){
  return g_scexpofit(x,&par[0])+g_gfit(x,&par[4]);  //offset expo fit
}
Double_t Tfit_p3(Double_t *x, Double_t *par){
  return g_p3fit(x,&par[0])+g_gfit(x,&par[4]);  //3rd order poly fit
}
Double_t Tfit_p4(Double_t *x, Double_t *par){
  return g_p4fit(x,&par[0])+g_gfit(x,&par[5]);  //4th order poly fit
}
Double_t Tfit_p5(Double_t *x, Double_t *par){
  return g_p5fit(x,&par[0])+g_gfit(x,&par[6]);  //5th order poly fit
}
Double_t Tfit_p6(Double_t *x, Double_t *par){
  return g_p6fit(x,&par[0])+g_gfit(x,&par[7]);  //6th order poly fit
}
Double_t Tfit_p8(Double_t *x, Double_t *par){
  return g_p8fit(x,&par[0])+g_gfit(x,&par[9]);  //8th order poly fit
}

//Configured only for LH2 data. Use -2(expo alt), -1(expo), 0(gaus), 1-8(poly)
void hde( Int_t kine=-1, Int_t mag=-1, Int_t W2fit=-1  ){ //main

  if( W2fit == 0 || 
      W2fit == 1 || 
      W2fit == 2 || 
      W2fit == 7 || 
      W2fit > 8 || 
      W2fit < -2 ){
    cout << "ERROR: W2fit not configured for setting. Try again." << endl;
    return;
  }

  const char *tar = "LH2"; //Should only consider LH2 for HCal detection efficiency for now

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  TChain *C = new TChain("T");

  // First obtain index for nfset constants
  Int_t kIdx=-1;
  for( Int_t i=0; i<nkine; i++ ){
    if( gIdx[i]==kine ) kIdx=i;
  }
  if( kIdx==-1 ){
    cout << "Error: Input parameters out of bounds. Please execute with the format:" << endl;
    cout << "  root -l \"hcal_detEff.C( <kine>,  <mag> )\" " << endl;
    cout << "  ..where kine = { 4, 7, 8, 9, 11, 14 }" << endl;
    cout << "  ..where mag = { 0, 30, 50, 70, 85, 100 }" << endl;
    return;
  }
  
  // Declare simple variables and strings to keep track of file extensions
  const char *config_prefix = gSystem->Getenv("CONFIG_DIR");

  TString configfilename = Form("%s/hcal_conf/SBS%d/secal_lh2_sbs%d_f%d.cfg",config_prefix,kine,kine,mag);
  TString outputfilename = Form("outfiles/HDE_lh2data_SBS%d_mag%d.root", kine, mag );
  TString reportPath = Form("efficiencyReports/effReport_sbs%d_mag%d.txt",kine,mag);

  // Declare general physics parameters to be modified by input config file
  Double_t E_e = -1.; // Energy of beam (incoming electrons from accelerator)
  Double_t HCal_d = -1.; // Distance to HCal from scattering chamber for comm1
  Double_t HCal_th = -1.; // Angle HCal center makes with exit beamline  
  Double_t BB_th = -1.; // Angle BB spectrometer makes with exit beamline
  Double_t W2_mean = -1.; // Mean of W at current kinematic
  Double_t W2_sig = -1.; // Width of W at current kinematic
  Double_t dx0_n = -1; // Peak location of neutron delta x. NOT USED
  Double_t dx0_p = -1; // Peak location of proton delta x
  Double_t dy0 = -1; // Peak location of delta y
  Double_t dx_sig_n = -1; // Sigma of neutron delta x peak. NOT USED
  Double_t dx_sig_p = -1; // Sigma of proton delta x peak
  Double_t dy_sig = -1; // Sigma of delta y peak
  Double_t atime0 = -1; // Peak position of ADC time signal
  Double_t atime_sig = -1; // Sigma of ADC time signal peak
  Int_t useAlshield = -1.; //Use 1/8" al shield on scattering chamber exit? 1:yes 0:no

  //EARM FIT TUNING
  ///////////////////////////////////////////////////////////
  Double_t W2fitmax = 1.4; //Max value in W2 plot range
  Double_t binfac = 400.; //Number of bins/unit W2
  Double_t confac = 2.; //Number of sigma in delta x/y considered a "detection"
  Double_t W2confac = 2.; //Number of sigma in W2 considered elastic
  Double_t dconfac = 10.; //Number of sigma outside of which inelastics only
  ///////////////////////////////////////////////////////////

  // Reading config file
  ifstream configfile(configfilename);
  TString currentline;
  cout << endl << "Chaining the following runs: " << endl;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endlist") ){
    if( !currentline.BeginsWith("#") ){
      if(!currentline) cout << "WARNING: No file exists at " << currentline << "." << endl;
      C->Add(currentline);
      if(currentline){
	cout << currentline << " ..check" << endl;
      }
    }    
  }

  TCut globalcut = "";
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("endcut") ){
    if( !currentline.BeginsWith("#") ){
      globalcut += currentline;
    }    
    cout << "Loading globalcut: " << globalcut << endl;
  }

  cout << endl << "Loading the following parameters from " << configfilename << ":" << endl;
  while( currentline.ReadLine( configfile ) && !currentline.BeginsWith("#") ){
    TObjArray *tokens = currentline.Tokenize(" ");
    Int_t ntokens = tokens->GetEntries();
    if( ntokens>1 ){
      TString skey = ( (TObjString*)(*tokens)[0] )->GetString();
      if( skey == "E_e" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	E_e = sval.Atof();
	cout << "Loading beam energy: " << E_e << endl;
      }
      if( skey == "HCal_d" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_d = sval.Atof();
	cout << "Loading HCal distance: " << HCal_d << endl;
      }
      if( skey == "HCal_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	HCal_th = sval.Atof() * TMath::DegToRad();	
	cout << "Loading HCal angle: " << HCal_th << endl;
      }
      if( skey == "BB_th" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	BB_th = sval.Atof() * TMath::DegToRad();	
	cout << "Loading BBCal angle: " << BB_th << endl;
      }
      if( skey == "W2_mean" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	W2_mean = sval.Atof();
	cout << "Loading W2 mean cut: " << W2_mean << endl;
      }
      if( skey == "W2_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	W2_sig = sval.Atof();
	cout << "Loading W2 sigma cut: " << W2_sig << endl;
      }
      if( skey == "dx0_n" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx0_n = sval.Atof();
	cout << "Loading x position of neutron spot: " << dx0_n << endl;
      }
      if( skey == "dx0_p" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx0_p = sval.Atof();
	cout << "Loading y position of proton spot: " << dx0_p << endl;
      }
      if( skey == "dy0" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dy0 = sval.Atof();
	cout << "Loading y position of both hadron spots: " << dy0 << endl;
      }
      if( skey == "dx_sig_n" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx_sig_n = sval.Atof();
	cout << "Loading x sigma of neutron spot: " << dx_sig_n << endl;
      }
      if( skey == "dx_sig_p" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dx_sig_p = sval.Atof();
	cout << "Loading x sigma of proton spot: " << dx_sig_p << endl;
      }
      if( skey == "dy_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	dy_sig = sval.Atof();
	cout << "Loading y sigma of both hadron spots: " << dy_sig << endl;
      }
      if( skey == "atime0" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	atime0 = sval.Atof();
	cout << "Loading ADC time mean: " << atime0 << endl;
      }
      if( skey == "atime_sig" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	atime_sig = sval.Atof();
	cout << "Loading ADC time sigma: " << atime_sig << endl;
      }
      if( skey == "useAlshield" ){
	TString sval = ( (TObjString*)(*tokens)[1] )->GetString();
	useAlshield = sval.Atoi();
	cout << "Loading Aluminum absorber option: " << useAlshield << endl;
      }
    }
    delete tokens;
  }
  
  if( E_e==-1 || 
      BB_th==-1 || 
      HCal_d==-1 || 
      HCal_th==-1 || 
      W2_mean==-1 || 
      W2_sig==-1 || 
      dx0_p==-1 || 
      dy0==-1 || 
      dx_sig_p==-1 || 
      dy_sig==-1 || 
      atime0==-1 || 
      atime_sig==-1 || 
      useAlshield==-1 ){
    cout << "Error: Unable to load setup parameters correctly. Please source setup_gmna.sh and try again. If problem persists, check configuration file." << endl;
    return;
  }

  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut, C );

  cout << "Setup parameters loaded." << endl;

  cout << endl;

  // Declare general detector and physics parameters. NOT USED (until tdc signals better understood)
  Double_t TDCT_id[maxTDCTrigChan], TDCT_tdc[maxTDCTrigChan]; 
  Int_t TDCTndata;

  // BB params
  Double_t BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  Double_t BBtr_vz[maxTracks], BBtr_chi2[maxTracks], BBtr_tgth[maxTracks], BBtr_tgph[maxTracks];
  Double_t BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e, GEMtr_hits, BBeoverp;

  // HCal params
  Double_t cid[maxHCalChan], crow[maxHCalRows], ccol[maxHCalCols];
  Double_t ce[maxHCalChan], catime[maxHCalChan], ctdc[maxHCalChan];
  Double_t nclus, nblk;
  Double_t hcalx, hcaly, hcale, kineW2;

  //Switch off all unnecessary branches
  C->SetBranchStatus( "*", 0 );
  C->SetBranchStatus( "sbs.hcal.x", 1 );
  C->SetBranchStatus( "sbs.hcal.y", 1 );
  C->SetBranchStatus( "sbs.hcal.e", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.row", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.col", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.tdctime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.atime", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.id", 1 );
  C->SetBranchStatus( "sbs.hcal.clus_blk.e", 1 );
  C->SetBranchStatus( "sbs.hcal.nblk", 1 );
  C->SetBranchStatus( "sbs.hcal.nclus", 1 );
  C->SetBranchStatus( "bb.tr.n", 1 );
  C->SetBranchStatus( "bb.tr.p", 1 );
  C->SetBranchStatus( "bb.tr.px", 1 );
  C->SetBranchStatus( "bb.tr.py", 1 );
  C->SetBranchStatus( "bb.tr.pz", 1 );
  C->SetBranchStatus( "bb.tr.vx", 1 );
  C->SetBranchStatus( "bb.tr.vy", 1 );
  C->SetBranchStatus( "bb.tr.vz", 1 );
  C->SetBranchStatus( "bb.tr.tg_th", 1 );
  C->SetBranchStatus( "bb.tr.tg_ph", 1 );
  C->SetBranchStatus( "bb.ps.e", 1 );
  C->SetBranchStatus( "bb.ps.x", 1 );
  C->SetBranchStatus( "bb.ps.y", 1 );
  C->SetBranchStatus( "bb.sh.e", 1 );
  C->SetBranchStatus( "bb.sh.x", 1 );
  C->SetBranchStatus( "bb.sh.y", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdc", 1 );
  C->SetBranchStatus( "bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus( "Ndata.bb.tdctrig.tdcelemID", 1 );
  C->SetBranchStatus( "e.kine.W2", 1 );
  C->SetBranchStatus( "bb.gem.track.nhits", 1 );
  C->SetBranchStatus( "bb.etot_over_p", 1 );

  //Link the branches to vars
  C->SetBranchAddress( "sbs.hcal.x", &hcalx );
  C->SetBranchAddress( "sbs.hcal.y", &hcaly );
  C->SetBranchAddress( "sbs.hcal.e", &hcale );
  C->SetBranchAddress( "sbs.hcal.clus_blk.row", crow );
  C->SetBranchAddress( "sbs.hcal.clus_blk.col", ccol );
  C->SetBranchAddress( "sbs.hcal.clus_blk.tdctime", ctdc );
  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", catime );
  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", cid );
  C->SetBranchAddress( "sbs.hcal.clus_blk.atime", ce );
  C->SetBranchAddress( "sbs.hcal.nblk", &nblk );
  C->SetBranchAddress( "sbs.hcal.nclus", &nclus );
  C->SetBranchAddress( "bb.tr.n", &BBtr_n );
  C->SetBranchAddress( "bb.tr.p", BBtr_p );
  C->SetBranchAddress( "bb.tr.px", BBtr_px );
  C->SetBranchAddress( "bb.tr.py", BBtr_py );
  C->SetBranchAddress( "bb.tr.pz", BBtr_pz );
  C->SetBranchAddress( "bb.tr.vz", BBtr_vz );
  C->SetBranchAddress( "bb.tr.tg_th", BBtr_tgth );
  C->SetBranchAddress( "bb.tr.tg_ph", BBtr_tgph );
  C->SetBranchAddress( "bb.ps.x", &BBps_x );
  C->SetBranchAddress( "bb.ps.y", &BBps_y );
  C->SetBranchAddress( "bb.ps.e", &BBps_e );
  C->SetBranchAddress( "bb.sh.x", &BBsh_x );
  C->SetBranchAddress( "bb.sh.y", &BBsh_y );
  C->SetBranchAddress( "bb.sh.e", &BBsh_e );
  C->SetBranchAddress( "bb.tdctrig.tdc", TDCT_tdc );
  C->SetBranchAddress( "bb.tdctrig.tdcelemID", TDCT_id );
  C->SetBranchAddress( "Ndata.bb.tdctrig.tdcelemID", &TDCTndata );
  C->SetBranchAddress( "e.kine.W2", &kineW2 );
  C->SetBranchAddress( "bb.gem.track.nhits", &GEMtr_hits );
  C->SetBranchAddress( "bb.etot_over_p", &BBeoverp );

  Int_t nentries = C->GetEntries();
  cout << endl;

  cout << "Opened tree with " << nentries << " entries." << endl;
  
  // Define some global variables
  Double_t emin = 0.0; // Minimum HCal energy to be considered as a hit, added for hard check on 1D distributions
  Double_t efficiency_rel; // Relative detection efficiency of HCAL (elastic events detected by HCal) / (elastic events as defined by BB tracking)
  Int_t hits_elBB = 0; // Count of total elastic hits detected in BB
  Int_t hits_gHCAL = 0; // Count of all elastic events that HCal detected
  Double_t pBeam = E_e/(1.+E_e/M_p*(1.-cos(BB_th))); // Momentum of beam calculated from theta
  Double_t Eloss_outgoing = celldiameter/2.0/sin(BB_th) * rho_tgt * dEdx_tgt; //Approximately 1 MeV, could correct further with raster position
  Double_t W2min_elastic = W2_mean - 2.*W2_sig;
  Double_t W2max_elastic = W2_mean + 2.*W2_sig;
  Int_t hcalbins = 500;
  //Fit limits for HCal/W2 elastics
  Double_t hcalfit_l = posHCalXi_MC-2*HCalblk_l_v_MC;
  Double_t hcalfit_h = posHCalXf_MC+2*HCalblk_l_v_MC;
  Double_t harmrange = (hcalfit_h) - (hcalfit_l);
  Double_t fit_l = 0.0;
  Double_t fit_h = W2fitmax;

  //Outfiles
  TFile *fout = new TFile(outputfilename,"RECREATE");

  //Histograms
  TH2D *hxy_nocut = new TH2D("hxy_nocut",";dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_hactivecut = new TH2D("hxy_hactivecut",";dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH1D *hvz = new TH1D("hvz","",250,-0.15,0.15);
  TH1D *hdpel_nocut = new TH1D("hdpel_nocut",";p/p_{elastic}(#theta)-1;", 250, -1.0, 0.5);
  TH1D *hdpel_BBcut = new TH1D("hdpel_BBcut",";p/p_{elastic}(#theta)-1;", 250, -1.0, 0.5);
  TH1D *hW2_nocut = new TH1D( "hW2_nocut","W2, no cut;W^2 (GeV);", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_nocut_wide = new TH1D( "hW2_nocut_wide","W2, no cut;W^2 (GeV);", binfac*6.0, 0.0, 6.0 ); //Make a large upper limit
  TH1D *hW2_BBcut = new TH1D( "hW2_BBcut","W2, elastic cut;W^2 (GeV);", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_BBcut_norange = new TH1D( "hW2_BBcut_norange",";W^2 (GeV);", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_BBcut_HCalcut = new TH1D( "hW2_BBcut_HCalcut",";W^2 (GeV);", binfac*W2fitmax, 0.0, W2fitmax );
  TH2D *hdxdy_HCAL_BBcut = new TH2D("hdxdy_HCAL_BBcut",";dy_{HCAL} (m); dx_{HCAL} (m)", 250, posHCalYi_MC-2*HCalblk_l_h_MC, posHCalYf_MC+2*HCalblk_l_h_MC, 250, hcalfit_l, hcalfit_h );
  TH2D *hdxdy_HCAL_nocut = new TH2D("hdxdy_HCAL_nocut",";dy_{HCAL} (m); dx_{HCAL} (m)", 250, posHCalYi_MC-2*HCalblk_l_h_MC, posHCalYf_MC+2*HCalblk_l_h_MC, 250, hcalfit_l, hcalfit_h );
  TH1D *hdx_HCAL_nocut = new TH1D("hdx_HCAL_nocut",";x_{HCAL}-x_{expect} (m);", hcalbins, hcalfit_l, hcalfit_h);
  TH1D *hdy_HCAL_nocut = new TH1D("hdy_HCAL_nocut",";y_{HCAL}-y_{expect} (m);", hcalbins, posHCalYi_MC-2*HCalblk_l_h_MC, posHCalYf_MC+2*HCalblk_l_h_MC);
  TH1D *hdx_HCAL_BBcut = new TH1D("hdx_HCAL_BBcut",";x_{HCAL}-x_{expect} (m);", hcalbins, hcalfit_l, hcalfit_h);
  TH1D *hdy_HCAL_BBcut = new TH1D("hdy_HCAL_BBcut",";y_{HCAL}-y_{expect} (m);", hcalbins, posHCalYi_MC-2*HCalblk_l_h_MC, posHCalYf_MC+2*HCalblk_l_h_MC);
  TH1D *hdx_HCAL_BBcut_dycut = new TH1D("hdx_HCAL_BBcut_dycut",";x_{HCAL}-x_{expect} (m);", hcalbins, hcalfit_l, hcalfit_h);

  TH1D *hW2_BBYcut = new TH1D( "hW2_BBYcut","W2, hcal dy cut;W^2 (GeV);", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_badYcut = new TH1D( "hW2_badYcut","W2, hcal dy anticut;W^2 (GeV);", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hdx_badW2cut = new TH1D( "hdx_badW2cut","hcal dx, W2 anticut;x_{HCAL}-x_{expect} (m);", hcalbins, hcalfit_l, hcalfit_h);
  TH1D *hdx_badYcut = new TH1D( "hdx_badYcut","hcal dx, hcal dy anticut;x_{HCAL}-x_{expect} (m);", hcalbins, hcalfit_l, hcalfit_h);

  //Accounting and diagnostic variables
  Long64_t nevent = 0;
  Int_t treenum = 0, currenttreenum = 0, currentrun = 0;

  cout << "Processing events.." << endl;

  Double_t progress = 0.;
  
  while( progress<1.0 ){
    
    Int_t barwidth = 70;
    Int_t step = 1;
    if( nevent >= nentries ) break;
    
    while( C->GetEntry( nevent++ ) ){

      //Single-loop globalcut
      currenttreenum = C->GetTreeNumber();
      if( nevent == 1 || currenttreenum != treenum ){
	treenum = currenttreenum; 
	GlobalCut->UpdateFormulaLeaves();
      }
      bool failedglobal = GlobalCut->EvalInstance(0) == 0;

      //Correct the beam energy with energy loss in target using vertex position
      Double_t Eloss = (BBtr_vz[0]+l_tgt/2.0) * rho_tgt * dEdx_tgt + uwallthick_LH2 * rho_Al * dEdx_Al; //approximately 3 MeV
      Double_t E_corr = E_e - Eloss;

      //Physics Calculations
      Double_t p_corr = BBtr_p[0] - Eloss_outgoing; //Neglecting the mass of e'
      Double_t etheta = acos( BBtr_pz[0]/BBtr_p[0] ); //Will use the uncorrected track momentum to reconstruct e' theta
      Double_t ephi = atan2( BBtr_py[0], BBtr_px[0] );   
      TVector3 vertex( 0, 0, BBtr_vz[0] ); // z location of vertex in hall coordinates
      TLorentzVector Pbeam( 0, 0, E_corr, E_corr ); //Mass of e negligable
      TLorentzVector kprime( BBtr_px[0], BBtr_py[0], BBtr_pz[0], BBtr_p[0] );
      TLorentzVector Ptarg( 0, 0, 0, M_p );   
      TLorentzVector q = Pbeam - kprime;
      TLorentzVector PgammaN = Ptarg + q; //(-px, -py, ebeam - pz, Mp + ebeam - p)     
      Double_t pel = E_corr/( 1.+E_corr/M_p*( 1.-cos(etheta) ) );
      Double_t nu = E_corr - BBtr_p[0];
      Double_t pp = sqrt( pow(nu,2)+2.*M_p*nu );
      Double_t phinucleon = ephi + PI; //assume coplanarity
      Double_t thetanucleon = acos( (E_corr - BBtr_pz[0])/pp ); //use elastic constraint on nucleon kinematics     
      TVector3 pNhat( sin(thetanucleon)*cos(phinucleon), sin(thetanucleon)*sin(phinucleon), cos(thetanucleon) );
      
      //Define HCal coordinate system
      TVector3 HCAL_zaxis( sin(-HCal_th), 0, cos(-HCal_th) );
      TVector3 HCAL_xaxis( 0, -1, 0 );
      TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();      
      TVector3 HCAL_origin = HCal_d * HCAL_zaxis + HCALHeight * HCAL_xaxis;     
      TVector3 HCALpos = HCAL_origin + hcalx * HCAL_xaxis + hcaly * HCAL_yaxis;
      Double_t sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis ) / (pNhat.Dot( HCAL_zaxis ) ); 
      TVector3 HCAL_intersect = vertex + sintersect * pNhat;   

      //Calculate quantities of interest
      Double_t yexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_yaxis );
      Double_t xexpect_HCAL = (HCAL_intersect - HCAL_origin).Dot( HCAL_xaxis ); 
      Double_t dx = hcalx - xexpect_HCAL;
      Double_t dy = hcaly - yexpect_HCAL;
      //Double_t dr = sqrt( pow( dx, 2 )+pow( dy, 2 ));
      Double_t E_ep = sqrt( pow(M_e,2) + pow(BBtr_p[0],2) ); // Obtain the scattered electron energy 
      Double_t p_ep = BBtr_p[0];     
      Double_t Q2 = 2*E_corr*E_ep*( 1-(BBtr_pz[0]/p_ep) ); // Obtain Q2 from beam energy, outgoing electron energy, and momenta  
      Double_t W2 = kineW2; //Get invarient mass transfer W from the four-momentum of the scattered nucleon
      //Double_t W = PgammaN.M();
      //Double_t W = sqrt(kineW2);
      Double_t E_pp = nu+M_p; // Get energy of the proton
      Double_t Enucleon = sqrt(pow(pp,2)+pow(M_p,2)); // Check on E_pp, same     
      Double_t KE_p = nu; // For elastics
      //Double_t dx = hcalx - xexpect_HCAL;
      //Double_t dy = hcaly - yexpect_HCAL;  
      Double_t pelastic = E_corr/(1.+(E_corr/M_p)*(1.0-cos(etheta))); 	
      Double_t precon = BBtr_p[0] + Eloss_outgoing; //reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
      Double_t nu_recon = E_corr - precon;
      Double_t Q2recon = 2.0*E_corr*precon*(1.0-cos(etheta));
      Double_t W2recon = pow(M_p,2) + 2.0*M_p*nu_recon - Q2recon;      

      //Calculate q vector
      TLorentzVector Kprime_recon(precon*BBtr_px[0]/BBtr_p[0],precon*BBtr_py[0]/BBtr_p[0],precon*BBtr_pz[0]/BBtr_p[0],precon);
      TLorentzVector q_recon = Pbeam - Kprime_recon;
      TVector3 qvect_recon = q_recon.Vect();

      //Get theta pq for neutron and proton
      //Calculate expected neutron direction
      TVector3 NeutronDirection = (HCALpos - vertex).Unit();

      //Calculate expected proton direction with SBS deflection
      double BdL = SBSfield * maxSBSfield * Dgap;
      double proton_deflection = tan( 0.3 * BdL / qvect_recon.Mag() ) * (HCal_d - (SBSdist + Dgap/2.0) );
      TVector3 ProtonDirection = (HCALpos + proton_deflection * HCAL_xaxis - vertex).Unit();
	
      bool inboundsW2 = W2 < W2fitmax && W2 > 0.;

      hxy_nocut->Fill(yexpect_HCAL,xexpect_HCAL);

      //No events which project hadrons off the face of HCal should be considered
      bool offhcal = 
	yexpect_HCAL > (posHCalYf_MC-HCalblk_l_h_MC) ||
	yexpect_HCAL < (posHCalYi_MC+HCalblk_l_h_MC) ||
	xexpect_HCAL > (posHCalXf_MC-HCalblk_l_v_MC) ||
	xexpect_HCAL < (posHCalXi_MC+HCalblk_l_v_MC);
      
      // if(  yexpect_HCAL > (posHCalYf_MC-HCalblk_l_h_MC) ||
      // 	   yexpect_HCAL < (posHCalYi_MC+HCalblk_l_h_MC) ||
      // 	   xexpect_HCAL > (posHCalXf_MC-HCalblk_l_v_MC) ||
      // 	   xexpect_HCAL < (posHCalXi_MC+HCalblk_l_v_MC)){
      // 	continue;
      // }
      if( offhcal ) continue;

      hxy_hactivecut->Fill(yexpect_HCAL,xexpect_HCAL);

      //Fill histograms without any non-hcal-active-area cuts
      hdpel_nocut->Fill( BBtr_p[0]/pel - 1.0 );
      if( inboundsW2 ) hW2_nocut->Fill( W2 ); //Subtract dy anticut version from this
      hW2_nocut_wide->Fill( W2 );
      hvz->Fill( vertex.Z() );
      hdxdy_HCAL_nocut->Fill( dy, dx );
      hdx_HCAL_nocut->Fill( dx );  //Subtract W2/dy anticut version from this
      hdy_HCAL_nocut->Fill( dy );
	  	
      //Set up bools for e-arm cuts by kinematic
      // bool earmcut[nkine] = {
      // 	BBtr_n!=1||BBps_e<=0.141||abs(BBtr_vz[0])>=0.08||GEMtr_hits<3||abs(BBeoverp-0.99)>0.27,
      // 	BBtr_n<=0||BBps_e<=0.2||abs(BBtr_vz[0])>=0.08||GEMtr_hits<3,
      // 	BBtr_n<=0||BBps_e<=0.2||abs(BBtr_vz[0])>=0.08||GEMtr_hits<3||BBtr_p[0]<=2.0,
      // 	BBtr_n<=0||BBps_e<=0.2||abs(BBtr_vz[0])>=0.08||GEMtr_hits<3||BBtr_p[0]<=2.0,
      // 	BBtr_n!=1||BBps_e<=0.2||abs(BBtr_vz[0])>=0.08||GEMtr_hits<3||abs(BBtr_tgth[0])>=0.15||abs(BBtr_tgph[0])>=0.3,	
      // 	BBtr_n!=1||BBps_e<=0.2||abs(BBtr_vz[0])>=0.08||GEMtr_hits<3||abs(BBtr_tgth[0])>=0.15||abs(BBtr_tgph[0])>=0.3,
      //};
      
      //E-arm only cuts to gain tight selection on elastics for comparison with HCal (not generally used)
      ////////////////////////////////////////////////////////////////////////////////////////////////////
      //if( BBtr_n!=1||BBps_e<0.141||abs(BBtr_vz[0])>0.08||GEMtr_hits<3||abs(BBeoverp-0.99)>0.27 ) continue;
      //if( earmcut[kIdx] ) continue;
      //if( failedglobal ) continue;
      ////////////////////////////////////////////////////////////////////////////////////////////////////

      bool elastic = abs( W2 - W2_mean ) < W2confac*W2_sig && !failedglobal;

      //Cuts useful for anticut analysis
      bool gW2 = elastic; //Pull all the stops to find W2 background
      bool gY = abs(dy-dy0) < confac*dy_sig && !failedglobal; //Pull all the stops to find dx background
      //bool gY = abs(dy-dy0) < HCalblk_l_h_MC && !failedglobal; //consider pos res of one block width

      if( W2<=W2fitmax && W2>=0.0 ) hW2_BBcut->Fill( W2 );

      //Consider efficiency only on HCal Y slice corresponding to position resolution HCal

      // if( dy<-0.6 & dy>0.6 ){
      // 	hW2_badYcut->Fill( W2 ); //Anticut W2 histo
      // 	hdx_badYcut->Fill( dx ); //Alt anticut dx histo
      // }
      if( dy>0.7 ){
      	hW2_badYcut->Fill( W2 ); //Anticut W2 histo
      	hdx_badYcut->Fill( dx ); //Alt anticut dx histo
      }

      // if( !gY ){
      // 	hW2_badYcut->Fill( W2 ); //Anticut W2 histo
      // 	hdx_badYcut->Fill( dx ); //Alt anticut dx histo
      // }

      // if( !gW2 ){
      // 	hdx_badW2cut->Fill( dx ); //Anticut dx histo
      // }

      if( W2<0.5 && W2>0. ){
	hdx_badW2cut->Fill( dx ); //Anticut dx histo
      }

      //Can continue here if we want elastic cuts
      if( !gW2 ) continue;

      //Fill histograms with BigBite and HCal active area cuts
      hdpel_BBcut->Fill( BBtr_p[0]/pel - 1.0 );
      hW2_BBcut_norange->Fill( W2 );
      if( inboundsW2 ) hW2_BBYcut->Fill( W2 );
      hdxdy_HCAL_BBcut->Fill( dy, dx );
      hdx_HCAL_BBcut->Fill( dx );
      hdy_HCAL_BBcut->Fill( dy );

      if( !gY ) continue;

      hdx_HCAL_BBcut_dycut->Fill( dx ); //Harm elastics (efficiency numerator)
      if( abs(dx-dx0_p)<confac*dx_sig_p ){
	hW2_BBcut_HCalcut->Fill( W2 );
      }

      cout << "[";
      Int_t pos = barwidth*progress;
      for( Int_t i=0; i<barwidth; ++i ){
	if( i<pos ) cout << "_";
	else if( i==pos ){ 
	  
	  if( step%4==0 ){
	    cout << "(>0_0)>";
	  }
	  if( step%4==1 ){
	    cout << "<(0_0)>";
	  }
	  if( step%4==2 ){
	    cout << "<(0_0<)";
	  }
	  if( step%4==3 ){
	    cout << "<( ; )>";
	  }
	  
	}
	else cout << " ";
      }
      progress = (Double_t)( ( nevent+1. )/nentries );
      
      cout << "]" << Int_t( progress*100 ) << "%\r";
      cout.flush();
      if( nevent%10000==0 ) step++;
      
    }
    
  }

  //dy anti-cut BG subtraction methods
  TH1D *hdx_HCAL_nocut_clone = (TH1D*)(hdx_HCAL_nocut->Clone("hdx_HCAL_nocut_clone")); //for later fit option

  //Make canvas for W2 dy anticut analysis
  TCanvas *c1 = new TCanvas("c1",Form("hcal_sbs%d%smag%d_W2_dyanticut",kine, tar, mag),1200,500);
  gStyle->SetPalette(53);
  c1->Divide(2,1);
  c1->cd(1);

  //Get difference between gY and !gY W2 after scaling to normalize
  hW2_nocut->SetLineColor(kBlack);
  hW2_nocut->SetTitle("W2");
  hW2_nocut->Draw();
  Double_t gY_lastbin = hW2_nocut->GetBinContent(binfac*W2fitmax);
  Double_t bY_lastbin = hW2_badYcut->GetBinContent(binfac*W2fitmax);

  Double_t BGfac = gY_lastbin/bY_lastbin;

  cout << endl << "e-arm dy anticut background factor: " << BGfac << endl;

  //Scale W2 "background" dy anticut histo up/down to match last bin of W2 nocut histo
  TH1D *hW2_badYcut_scaled = (TH1D*)(hW2_badYcut->Clone("hW2_badYcut_scaled"));
  hW2_badYcut_scaled->Scale(BGfac);
  hW2_badYcut_scaled->SetLineColor(kRed);
  hW2_badYcut_scaled->Draw("same hist");
  TH1D *hW2_BGsub = (TH1D*)(hW2_nocut->Clone("hW2_BGsub"));
  hW2_BGsub->Add(hW2_badYcut_scaled,-1);
  for( Int_t i=1; i<=binfac*W2fitmax; i++ ){
    Double_t binc = hW2_BGsub->GetBinContent(i);
    if( binc<0. ) hW2_BGsub->SetBinContent(i,0.);
  }

  //Add a legend to the canvas
  auto W2bgsublegend = new TLegend(0.1,0.8,0.5,0.9);
  W2bgsublegend->SetTextSize( 0.03 );
  W2bgsublegend->AddEntry( hW2_nocut, "No Cut", "l");
  W2bgsublegend->AddEntry( hW2_badYcut_scaled, Form("dy<%d*sigma anticut (scaled)",(int)confac), "l");
  W2bgsublegend->Draw();

  c1->cd(2);
  hW2_BGsub->SetLineColor(kBlue);
  hW2_BGsub->SetTitle("W2, BG Subtracted");
  hW2_BGsub->Draw("hist");

  Double_t elastics_alt = hW2_BGsub->Integral(); //integrate all events in histogram after subtraction
  
  Double_t eel_bincontents = 0.;
  for( Int_t i=1; i<=binfac*W2fitmax; i++ ){ //bins start at 1
    eel_bincontents += hW2_BGsub->GetBinContent(i);
  }

  c1->Write();

  //Make canvas for hcal W2 anticut analysis
  TCanvas *c2 = new TCanvas("c2",Form("hcal_sbs%d%smag%d_dx_W2anticut",kine, tar, mag),1200,500);
  gStyle->SetPalette(53);
  c2->Divide(2,1);
  c2->cd(1);

  //Get diff gW2 and !gW2 dx after scaling to normalize
  hdx_HCAL_nocut->SetLineColor(kBlack);
  hdx_HCAL_nocut->SetTitle("HCal dx");
  hdx_HCAL_nocut->Draw();
  // hdx_HCAL_BBcut_dycut->SetLineColor(kBlack);
  // hdx_HCAL_BBcut_dycut->SetTitle("HCal dx");
  // hdx_HCAL_BBcut_dycut->Draw();

  Double_t gW2_hlastbin = hdx_HCAL_nocut->GetBinContent(hcalbins);
  //Double_t gW2_hlastbin = hdx_HCAL_BBcut_dycut->GetBinContent(hcalbins);
  Double_t bW2_hlastbin = hdx_badW2cut->GetBinContent(hcalbins);

  Double_t hcalBGfac = gW2_hlastbin/bW2_hlastbin;

  cout << endl << "h-arm W2 anticut background factor: " << hcalBGfac << endl;

  //Scale dx "background" W2 anticut histo up/down to match last bin of dx nocut histo
  TH1D *hdx_badW2cut_scaled = (TH1D*)(hdx_badW2cut->Clone("hdx_badW2cut_scaled"));
  hdx_badW2cut_scaled->Scale(hcalBGfac);
  hdx_badW2cut_scaled->SetLineColor(kRed);
  hdx_badW2cut_scaled->Draw("same hist");
  TH1D *hdx_BGsub = (TH1D*)(hdx_HCAL_nocut->Clone("hdx_BGsub"));
  //TH1D *hdx_BGsub = (TH1D*)(hdx_HCAL_BBcut_dycut->Clone("hdx_BGsub"));
  hdx_BGsub->Add(hdx_badW2cut_scaled,-1);
  for( Int_t i=1; i<=hcalbins; i++ ){
    Double_t binc = hdx_BGsub->GetBinContent(i);
    if( binc<0. ) hdx_BGsub->SetBinContent(i,0.);
  }

  //Add a legend to the canvas
  auto dxbgsublegend = new TLegend(0.1,0.8,0.5,0.9);
  dxbgsublegend->SetTextSize( 0.03 );
  dxbgsublegend->AddEntry( hdx_HCAL_nocut, "No Cut", "l");
  //dxbgsublegend->AddEntry( hdx_HCAL_BBcut_dycut, "BB elastic cut, dy cut", "l");
  dxbgsublegend->AddEntry( hdx_badW2cut_scaled, Form("W2<%d*sig anticut (scaled)",(int)W2confac), "l");
  dxbgsublegend->Draw();

  c2->cd(2);
  hdx_BGsub->SetLineColor(kBlue);
  hdx_BGsub->SetTitle("HCal dx, BG subtracted (W2 anticut)");
  hdx_BGsub->Draw("hist");

  Double_t hcalelastics_antiW2 = hdx_BGsub->Integral();

  Double_t hel_bincontents = 0.;
  for( Int_t i=1; i<=hcalbins; i++ ){ //bins start at 1
    hel_bincontents += hdx_BGsub->GetBinContent(i);
  }

  c2->Write();

  //Make canvas for hcal dy anticut analysis
  TCanvas *c0 = new TCanvas("c0",Form("hcal_sbs%d%smag%d_dx_dyanticut",kine, tar, mag),1200,500);
  gStyle->SetPalette(53);
  c0->Divide(2,1);
  c0->cd(1);

  //Get diff gY and !gY dx after scaling to normalize
  hdx_HCAL_nocut->SetLineColor(kBlack);
  hdx_HCAL_nocut->SetTitle("HCal dx");
  hdx_HCAL_nocut->Draw();
  // hdx_HCAL_BBcut_dycut->SetLineColor(kBlack);
  // hdx_HCAL_BBcut_dycut->SetTitle("HCal dx");
  // hdx_HCAL_BBcut_dycut->Draw();


  Double_t gY_hlastbin = hdx_HCAL_nocut->GetBinContent(hcalbins);
  //Double_t gY_hlastbin = hdx_HCAL_BBcut_dycut->GetBinContent(hcalbins);
  Double_t bY_hlastbin = hdx_badYcut->GetBinContent(hcalbins);

  Double_t hcalYBGfac = gY_hlastbin/bY_hlastbin;

  cout << endl << "h-arm dy anticut background factor: " << hcalYBGfac << endl;

  //Scale dx "background" W2 anticut histo up/down to match last bin of dx nocut histo
  TH1D *hdx_badYcut_scaled = (TH1D*)(hdx_badYcut->Clone("hdx_badYcut_scaled"));
  hdx_badYcut_scaled->Scale(hcalYBGfac);
  hdx_badYcut_scaled->SetLineColor(kRed);
  hdx_badYcut_scaled->Draw("same hist");
  TH1D *hdx_YBGsub = (TH1D*)(hdx_HCAL_nocut->Clone("hdx_YBGsub"));
  //TH1D *hdx_YBGsub = (TH1D*)(hdx_HCAL_BBcut_dycut->Clone("hdx_YBGsub"));
  hdx_YBGsub->Add(hdx_badYcut_scaled,-1);
  for( Int_t i=1; i<=hcalbins; i++ ){
    Double_t binc = hdx_YBGsub->GetBinContent(i);
    if( binc<0. ) hdx_YBGsub->SetBinContent(i,0.);
  }


  //Add a legend to the canvas
  auto dxabgsublegend = new TLegend(0.1,0.8,0.5,0.9);
  dxabgsublegend->SetTextSize( 0.03 );
  dxabgsublegend->AddEntry( hdx_HCAL_nocut, "No Cut", "l");
  //dxabgsublegend->AddEntry( hdx_HCAL_BBcut_dycut, "BB elastic cut, dy cut", "l");
  dxabgsublegend->AddEntry( hdx_badYcut_scaled, Form("dy<%d*sig anticut (scaled)",(int)confac), "l");
  dxabgsublegend->Draw();

  c0->cd(2);
  hdx_YBGsub->SetLineColor(kBlue);
  hdx_YBGsub->SetTitle("Hcal dx, BG subtracted (dy anticut)");
  hdx_YBGsub->Draw("hist");

  Double_t hcalelastics_antidy = hdx_YBGsub->Integral();

  Double_t hely_bincontents = 0.;
  for( Int_t i=1; i<=hcalbins; i++ ){ //bins start at 1
    hely_bincontents += hdx_YBGsub->GetBinContent(i);
  }

  c0->Write();

  //Make canvas for dxdy to show dy cut
  TCanvas *c3 = new TCanvas("c3",Form("hcal_sbs%d%smag%d",kine, tar, mag),450,750);
  gStyle->SetPalette(53);
  c3->cd();
  hdxdy_HCAL_nocut->GetYaxis()->SetTitleOffset(1.3); //can use BBcut
  hdxdy_HCAL_nocut->Draw("colz");
  TLine l_pl; //left cut
  l_pl.SetLineColor(kBlue);
  l_pl.SetLineWidth(3);
  l_pl.DrawLine(dy0-confac*dy_sig,hcalfit_l,dy0-confac*dy_sig,hcalfit_h);
  TLine l_pr; //right cut
  l_pr.SetLineColor(kBlue);
  l_pr.SetLineWidth(3);
  l_pr.DrawLine(dy0+confac*dy_sig,hcalfit_l,dy0+confac*dy_sig,hcalfit_h);
  c3->Write();

  //Make canvas to hold fitted W2 and HCal dx histograms
  TCanvas *c4 = new TCanvas("c4",Form("W2_sbs%d%smag%d",kine, tar, mag),1200,500);
  c4->Divide(2,1);
  gStyle->SetPalette(53);

  ///////////////
  ////HCal dx fit
  c4->cd(1);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  //Declare initial fit params for hcal dx
  vector<Double_t> hcalsetpar; //keep consistent with W2 for later alt fits
  hcalsetpar.push_back(0.0); //pol4 p0
  hcalsetpar.push_back(1.0); //pol4 p1
  hcalsetpar.push_back(1.0); //pol4 p2
  hcalsetpar.push_back(1.0); //pol4 p3
  hcalsetpar.push_back(1.0); //pol4 p4
  //hcalsetpar.push_back(2924.74); //gaussian amplitude
  hcalsetpar.push_back(15000); //gaussian amplitude
  hcalsetpar.push_back(dx0_p); //gaussian mean
  //hcalsetpar.push_back(dx_sig_p); //gaussian sigma
  hcalsetpar.push_back(0.074); //gaussian sigma

  //Declare total fit function
  TF1 *hcaltotalfit = new TF1("hcaltotalfit",Tfit_hcal,hcalfit_l,hcalfit_h,8); //total fit

  //Set the initial fit params
  hcaltotalfit->SetParameter(0,hcalsetpar[0]);
  hcaltotalfit->SetParameter(1,hcalsetpar[1]);
  hcaltotalfit->SetParameter(2,hcalsetpar[2]);
  hcaltotalfit->SetParameter(3,hcalsetpar[3]);
  hcaltotalfit->SetParameter(4,hcalsetpar[4]);
  hcaltotalfit->SetParameter(5,hcalsetpar[5]);
  hcaltotalfit->SetParameter(6,hcalsetpar[6]);
  hcaltotalfit->SetParameter(7,hcalsetpar[7]);
  
  hcaltotalfit->SetParLimits(5,10000,100000);
  hcaltotalfit->SetParLimits(6,-1.,0.0);
  hcaltotalfit->SetParLimits(7,0.05,0.5);

  //Fit the dy cut HCal dx histogram with total fit
  hcaltotalfit->SetLineColor(kGreen);
  hcaltotalfit->SetLineWidth(4);
  hdx_HCAL_BBcut_dycut->Fit("hcaltotalfit","RBM"); //Consider dx distribution in HCal with full elastic selection including dy
  //hdx_HCAL_BBcut->Fit("hcaltotalfit","RBM"); //Consider dx distribution in HCal with elastic selection no dy
  //hdx_HCAL_nocut_clone->Fit("hcaltotalfit","RBM"); //Consider dx distribution in HCal no cuts

  //Get parameters from total fit and generate identical bg and elastic functions to write on canvas
  Double_t *hcaltpar = hcaltotalfit->GetParameters();

  TF1 *hcalbg = new TF1("hcalbg",g_p4fit,hcalfit_l,hcalfit_h,5); //inelastic pol4 background
  hcalbg->SetParameters(&hcaltpar[0]);
  hcalbg->SetLineColor(kRed);
  hcalbg->Draw("same");

  TF1 *p = new TF1("p",g_gfit,hcalfit_l,hcalfit_h,3); //proton elastic gaussian
  p->SetParameters(&hcaltpar[5]);
  p->SetLineColor(kBlack);
  p->Draw("same");

  //Get elastics/inelastics from normalized integrals of fit functions
  Double_t hcalinelastics = (hcalbins/harmrange)*hcalbg->Integral(hcalfit_l,hcalfit_h);
  Double_t hcalelastics = (hcalbins/harmrange)*p->Integral(hcalfit_l,hcalfit_h);

  //Add a legend to the canvas
  auto hcallegend = new TLegend(0.1,0.7,0.48,0.9);
  hcallegend->SetTextSize( 0.03 );
  hcallegend->SetHeader( "HCal dx, E Arm Cuts and dy Cut", "C" );
  hcallegend->AddEntry( hcalbg, Form("Inelastic Fit: %d ev", (int)hcalinelastics), "l");
  hcallegend->AddEntry( p, Form("Elastic Fit: %d ev", (int)hcalelastics), "l");
  hcallegend->AddEntry( hcaltotalfit, "Total Fit", "l");
  hcallegend->Draw();

  //////////////////
  ////Bigbite W2 fit
  c4->cd(2);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);

  //Declare initial fit params for W2
  vector<Double_t> setpar;
  Double_t glimh[3];
  Double_t gliml[3];
  Double_t earmEvents = hW2_BBcut->GetEntries();

  //Set some reasonable limits for gaussian fit to elastic peak
  gliml[0] = earmEvents/binfac/2; //Amplitude constrained by number of events and binfac
  glimh[0] = earmEvents;
  gliml[1] = 0.75; //Mean constrained by W2 peak for elastics
  glimh[1] = 1.25;
  gliml[2] = 0.05; //Empirical limits
  glimh[2] = 0.50;

  //W2 overall fits
  TF1 *totalfit;
  TF1 *bg;
  TF1 *e = new TF1("e",g_gfit,fit_l,fit_h,3); //The elastic peak always gaussian

  //expo alt fit
  if( W2fit==-2 ){
    setpar.push_back(35);
    setpar.push_back(465.);
    setpar.push_back(-3.5);
    setpar.push_back(3.5);
    setpar.push_back(3057.);
    setpar.push_back(0.8077);
    setpar.push_back(0.0961);

    totalfit = new TF1("totalfit",Tfit_scexpo,fit_l,fit_h,7);
    bg = new TF1("bg",g_scexpofit,fit_l,fit_h,4);

    totalfit->FixParameter(0,setpar[0]);
    totalfit->FixParameter(1,setpar[1]);
    totalfit->FixParameter(2,setpar[2]);
    totalfit->FixParameter(3,setpar[3]);
    totalfit->SetParameter(4,setpar[4]);
    totalfit->SetParameter(5,setpar[5]);
    totalfit->SetParameter(6,setpar[6]);
    totalfit->SetParLimits(4,gliml[0],glimh[0]);
    totalfit->SetParLimits(5,gliml[1],glimh[1]);
    totalfit->SetParLimits(6,gliml[2],glimh[2]);
  }

  //expo fit
  if( W2fit==-1 ){
    setpar.push_back(465.0);
    setpar.push_back(-3.5);
    setpar.push_back(3.5);
    setpar.push_back(3057.);
    setpar.push_back(0.8077);
    setpar.push_back(0.0961);

    totalfit = new TF1("totalfit",Tfit_expo,fit_l,fit_h,6);
    bg = new TF1("bg",g_expofit,fit_l,fit_h,3);

    totalfit->SetParameter(0,setpar[0]);
    totalfit->SetParLimits(0,0.,10000.);
    totalfit->SetParameter(1,setpar[1]);
    totalfit->SetParLimits(1,-10000.,0.);
    totalfit->SetParameter(2,setpar[2]);
    totalfit->SetParLimits(2,0.,10000.);
    totalfit->SetParameter(3,setpar[3]);
    totalfit->SetParameter(4,setpar[4]);
    totalfit->SetParameter(5,setpar[5]);
    totalfit->SetParLimits(3,gliml[0],glimh[0]);
    totalfit->SetParLimits(4,gliml[1],glimh[1]);
    totalfit->SetParLimits(5,gliml[2],glimh[2]);
  }

  //3rd order poly fit
  if( W2fit==3 ){
    setpar.push_back(0.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(3057.);
    setpar.push_back(0.8077);
    setpar.push_back(0.0961);

    totalfit = new TF1("totalfit",Tfit_p3,fit_l,fit_h,7);
    bg = new TF1("bg",g_p3fit,fit_l,fit_h,4);

    totalfit->SetParameter(0,setpar[0]);
    totalfit->SetParameter(1,setpar[1]);
    totalfit->SetParameter(2,setpar[2]);
    totalfit->SetParameter(3,setpar[3]);
    totalfit->SetParameter(4,setpar[4]);
    totalfit->SetParameter(5,setpar[5]);
    totalfit->SetParameter(6,setpar[6]);
    totalfit->SetParLimits(4,gliml[0],glimh[0]);
    totalfit->SetParLimits(5,gliml[1],glimh[1]);
    totalfit->SetParLimits(6,gliml[2],glimh[2]);
  }

  //4th order poly fit
  if( W2fit==4 ){
    setpar.push_back(0.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(3057.);
    setpar.push_back(0.8077);
    setpar.push_back(0.0961);

    totalfit = new TF1("totalfit",Tfit_p4,fit_l,fit_h,8);
    bg = new TF1("bg",g_p4fit,fit_l,fit_h,5);

    totalfit->SetParameter(0,setpar[0]);
    totalfit->SetParameter(1,setpar[1]);
    totalfit->SetParameter(2,setpar[2]);
    totalfit->SetParameter(3,setpar[3]);
    totalfit->SetParameter(4,setpar[4]);
    totalfit->SetParameter(5,setpar[5]);
    totalfit->SetParameter(6,setpar[6]);
    totalfit->SetParameter(7,setpar[7]);
    totalfit->SetParLimits(5,gliml[0],glimh[0]);
    totalfit->SetParLimits(6,gliml[1],glimh[1]);
    totalfit->SetParLimits(7,gliml[2],glimh[2]);
  }

  //5th order poly fit
  if( W2fit==5 ){
    setpar.push_back(0.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(3057.);
    setpar.push_back(0.8077);
    setpar.push_back(0.0961);

    totalfit = new TF1("totalfit",Tfit_p5,fit_l,fit_h,9);
    bg = new TF1("bg",g_p5fit,fit_l,fit_h,6);

    totalfit->SetParameter(0,setpar[0]);
    totalfit->SetParameter(1,setpar[1]);
    totalfit->SetParameter(2,setpar[2]);
    totalfit->SetParameter(3,setpar[3]);
    totalfit->SetParameter(4,setpar[4]);
    totalfit->SetParameter(5,setpar[5]);
    totalfit->SetParameter(6,setpar[6]);
    totalfit->SetParameter(7,setpar[7]);
    totalfit->SetParameter(8,setpar[8]);
    totalfit->SetParLimits(6,gliml[0],glimh[0]);
    totalfit->SetParLimits(7,gliml[1],glimh[1]);
    totalfit->SetParLimits(8,gliml[2],glimh[2]);
  }

  //6th order poly fit
  if( W2fit==6 ){
    setpar.push_back(0.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(3057.);
    setpar.push_back(0.8077);
    setpar.push_back(0.0961);

    totalfit = new TF1("totalfit",Tfit_p6,fit_l,fit_h,10);
    bg = new TF1("bg",g_p6fit,fit_l,fit_h,7);

    totalfit->SetParameter(0,setpar[0]);
    totalfit->SetParameter(1,setpar[1]);
    totalfit->SetParameter(2,setpar[2]);
    totalfit->SetParameter(3,setpar[3]);
    totalfit->SetParameter(4,setpar[4]);
    totalfit->SetParameter(5,setpar[5]);
    totalfit->SetParameter(6,setpar[6]);
    totalfit->SetParameter(7,setpar[7]);
    totalfit->SetParameter(8,setpar[8]);
    totalfit->SetParameter(9,setpar[9]);
    totalfit->SetParLimits(7,gliml[0],glimh[0]);
    totalfit->SetParLimits(8,gliml[1],glimh[1]);
    totalfit->SetParLimits(9,gliml[2],glimh[2]);
  }

  //8th order poly fit
  if( W2fit==8 ){
    setpar.push_back(0.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(1.0);
    setpar.push_back(3057.);
    setpar.push_back(0.8077);
    setpar.push_back(0.0961);

    totalfit = new TF1("totalfit",Tfit_p8,fit_l,fit_h,12);
    bg = new TF1("bg",g_p8fit,fit_l,fit_h,9);

    totalfit->SetParameter(0,setpar[0]);
    totalfit->SetParameter(1,setpar[1]);
    totalfit->SetParameter(2,setpar[2]);
    totalfit->SetParameter(3,setpar[3]);
    totalfit->SetParameter(4,setpar[4]);
    totalfit->SetParameter(5,setpar[5]);
    totalfit->SetParameter(6,setpar[6]);
    totalfit->SetParameter(7,setpar[7]);
    totalfit->SetParameter(8,setpar[8]);
    totalfit->SetParameter(9,setpar[9]);
    totalfit->SetParameter(10,setpar[10]);
    totalfit->SetParameter(11,setpar[11]);
    totalfit->SetParLimits(9,gliml[0],glimh[0]);
    totalfit->SetParLimits(10,gliml[1],glimh[1]);
    totalfit->SetParLimits(11,gliml[2],glimh[2]);
  }

  //Fit the W2 histogram with total fit
  totalfit->SetLineColor(kGreen);
  totalfit->SetLineWidth(4);
  hW2_BBcut->Fit("totalfit","RBM"); //Fit the W2 distribution with total fit

  //Get all parameters
  Double_t *tpar = totalfit->GetParameters();

  //Get fit parameters for bg function and draw identical function on canvas
  bg->SetParameters(&tpar[0]);
  bg->SetLineColor(kRed);
  bg->Draw("same");

  //Get fit parameters for e gaussian and draw identical gaussian on canvas
  e->SetParameters(&tpar[3]);
  e->SetLineColor(kBlack);
  e->Draw("same");

  //Integrate the elastic peak to get total elastics detected in BB
  Double_t inelastics = binfac*bg->Integral(0.0,W2fitmax); //scale with binfac
  Double_t elastics = binfac*e->Integral(0.0,W2fitmax);

  //Add a legend
  auto legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->SetTextSize( 0.03 );
  legend->SetHeader( "W2, Electron Arm Cuts", "C" );
  legend->AddEntry( bg, Form("Inelastic Fit: %d ev", (int)inelastics), "l");
  legend->AddEntry( e, Form("Elastic Fit: %d ev", (int)elastics), "l");
  legend->AddEntry( totalfit, "Total Fit", "l");
  legend->Draw();

  c4->Write();
  fout->Write();

  //Calculate error assuming fit error is negligable
  //Binomial
  Double_t eff = hcalelastics/elastics;
  Double_t eff_alt = hcalelastics_antiW2/elastics_alt;
  Double_t eff_alt2 = hcalelastics_antidy/elastics_alt;
  //Double_t eff_althybrid = hcalelastics/elastics_alt;
  Double_t efferr = sqrt(eff*(1-eff)/elastics);
  Double_t efferr_alt = sqrt(eff_alt*(1-eff_alt)/elastics_alt);
  Double_t efferr_alt2 = sqrt(eff_alt2*(1-eff_alt2)/elastics_alt);
  //Double_t efferr_althybrid = sqrt(eff_althybrid*(1-eff_althybrid)/elastics_alt);

  //Declare outfile
  ofstream report;

  cout << "Chi^2 for total fit: " << totalfit->GetChisquare() << endl;
  report << "Chi^2 for total fit: " << totalfit->GetChisquare() << endl;

  report.open( reportPath );
  report << "HCal detection efficiency report for SBS" << kine << " LH2 data at " << mag << " percent field" << endl << endl;

  cout << "Total elastics detected in HCal: " << hcalelastics << endl;
  report << "Total elastics detected in HCal: " << hcalelastics << endl << endl;

  cout << "Total elastics detected in Bigbite with e-arm cuts only: " << elastics << endl;
  report << "Total elastics detected in Bigbite with e-arm cuts only: " << elastics << endl << endl;

  cout << "Total elastics detected in HCal (antidy): " << hcalelastics_antidy << endl;
  report << "Total elastics detected in HCal (antidy): " << hcalelastics_antidy << endl << endl;

  cout << "Total elastics detected in HCal (alt): " << hcalelastics_antiW2 << endl;
  report << "Total elastics detected in HCal (alt): " << hcalelastics_antiW2 << endl << endl;

  cout << "Total elastics detected in Bigbite with e-arm cuts only (alt): " << elastics_alt << endl;
  report << "Total elastics detected in Bigbite with e-arm cuts only (alt): " << elastics_alt << endl << endl;

  cout << "Total entries in W2 plot: " << earmEvents << endl << endl;
  report << "Total entries in W2 plot: " << earmEvents << endl << endl;

  cout << "==============================================================" << endl;
  cout << "Total HCal elastic detection efficiency (fits): " << eff << " +/- " << efferr << endl;
  cout << "==============================================================" << endl << endl;

  cout << "===================================================================" << endl;
  cout << "Total HCal elastic detection efficiency (anticut W2): " << eff_alt << " +/- " << efferr_alt << endl;
  cout << "===================================================================" << endl << endl;

  cout << "===================================================================" << endl;
  cout << "Total HCal elastic detection efficiency (anticut dy): " << eff_alt2 << " +/- " << efferr_alt2 << endl;
  cout << "===================================================================" << endl << endl;

  // cout << "===================================================================" << endl;
  // cout << "Total HCal elastic detection efficiency (althybrid): " << eff_althybrid << " +/- " << efferr_althybrid << endl;
  // cout << "===================================================================" << endl << endl;

  report << "==============================================================" << endl;
  report << "Total HCal elastic detection efficiency (fits): " << eff << " +/- " << efferr <<  endl;
  report << "==============================================================" << endl;

  report << "===================================================================" << endl;
  report << "Total HCal elastic detection efficiency (anticut W2): " << eff_alt << " +/- " << efferr_alt << endl;
  report << "===================================================================" << endl << endl;

  report << "===================================================================" << endl;
  report << "Total HCal elastic detection efficiency (anticut dy): " << eff_alt2 << " +/- " << efferr_alt2 << endl;
  report << "===================================================================" << endl << endl;


  cout << endl << endl << "raw bincontent efficiency anticut dy: " << hely_bincontents / eel_bincontents << endl;
  cout << endl << endl << "raw bincontent efficiency anticut W2: " << hel_bincontents / eel_bincontents << endl;

  report.close();

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
