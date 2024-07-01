//SSeeds 9.1.23 Script to aggregate simc generator mc output by nucleon and build histograms for further analysis

#include <vector>
#include <iostream>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "TMatrixD.h"
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

//////////////////////////////
//Manual override data options
bool norm_override = false;
Double_t Ntried_override = 100000;
Double_t luminosity_override = 3.8475e+36;
Double_t genvol_override = 12.566;

//main (type should either be "hpc" for jlab-HPC or "jboyd" for John Boyd's personal replays or "pdatta" for Provakar Datta's personal replays
void plotWithCuts_mc( Int_t kine=9, Int_t mag=70, const char *replay_type = "jboyd", const char *nucleon = "both", bool alt_infile = false, bool sync_jobs=true ){

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/gmn.json");

  std::string rootfile_dir = jmgr->GetValueFromSubKey_str( Form("rootfile_dir_mc_rc_%dp",mag), Form("sbs%d",kine) );

  std::string rootfile_dir_alt = "null";
  if(kine==4&&mag==30&&alt_infile)
    rootfile_dir_alt = jmgr->GetValueFromSubKey_str( Form("rootfile_dir_mc_rc_%dp_alt",mag), Form("sbs%d",kine) );

  std::string histfile_dir = jmgr->GetValueFromSubKey_str( Form("rootfile_dir_mc_rc_%dp_hist",mag), Form("sbs%d",kine) );
  std::string histfile_dir_alt = "null";
  if(kine==4&&mag==30&&alt_infile)
    histfile_dir_alt = jmgr->GetValueFromSubKey_str( "rootfile_dir_mc_rc_hist_alt", Form("sbs%d",kine) );

  Double_t W2fitmax = jmgr->GetValueFromSubKey<Double_t>( "W2fitmax", Form("sbs%d",kine) );
  Double_t binfac = jmgr->GetValueFromSubKey<Double_t>( "binfac", Form("sbs%d",kine) );
  Double_t hbinfac = jmgr->GetValueFromSubKey<Double_t>( "hbinfac", Form("sbs%d",kine) ); //gets to 7cm HCal pos res
  std::string target = jmgr->GetValueFromKey_str( "target" ); //Must always be LD2 for GMn analysis
  Int_t pass = jmgr->GetValueFromSubKey<Int_t>( "pass", Form("sbs%d",kine) );
  Double_t minE = jmgr->GetValueFromSubKey<Double_t>( "minE", Form("sbs%d",kine) );
  Double_t ebound_l = jmgr->GetValueFromSubKey<Double_t>( "ebound_l", Form("sbs%d",kine) );
  Double_t ebound_h = jmgr->GetValueFromSubKey<Double_t>( "ebound_h", Form("sbs%d",kine) );
  Double_t coin_sigma_factor = jmgr->GetValueFromSubKey<Double_t>( "coin_sigma_factor", Form("sbs%d",kine) );
  std::string partialName_p = jmgr->GetValueFromSubKey_str( "partial_name_p", Form("sbs%d",kine) );
  std::string partialName_n = jmgr->GetValueFromSubKey_str( "partial_name_n", Form("sbs%d",kine) );

  if( pass>1 ){
    std::cout << "As of 3.31.23, the highest GMn replay pass is 1. Enter a valid pass." << endl;
    return;
  }

  //set up default parameters for all analysis
  std::string runsheet_dir = "/w/halla-scshelf2102/sbs/seeds/ana/data"; //unique to my environment for now
  Int_t nruns = -1; //Always analyze all available runs
  Int_t verb = 0; //Don't print diagnostic info by default
  Double_t hcalfit_l = econst::hcalposXi_mc; //lower fit/bin limit for hcal dx plots (m)
  Double_t hcalfit_h = econst::hcalposXf_mc; //upper fit/bin limit for hcal dx plots (m)
  Double_t harmrange = econst::hcal_vrange; //Full range of hcal dx plots (m)

  //set up output files
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fout_path = outdir_path + Form("/gmn_analysis/mc_rc_%s_out_sbs%d_mag%d.root",nucleon,kine,mag);
  TFile *fout = new TFile( fout_path.c_str(), "RECREATE" );
  //TFile *fout = new TFile( "test.root", "RECREATE" );


  //file search on .root and .csv/.hist, get extensions. Proton first
  std::vector<std::string> rootFileNames_p, histFileNames_p; 
  std::vector<std::string> rootFileNames_p_alt, histFileNames_p_alt; 
  std::string histfile_dir_p;
  std::string rootfile_dir_p;
 
  std::string rtype = replay_type;
  if( rtype.compare("jboyd")==0 )
    util::FindMatchingFiles(histfile_dir,rootfile_dir,partialName_p,histFileNames_p,rootFileNames_p,true);
  else if( rtype.compare("pdatta")==0 ){
    util::FindMatchingFiles(histfile_dir,rootfile_dir,partialName_p,histFileNames_p,rootFileNames_p,false);
  }else{

    histfile_dir_p = histfile_dir + "deep";
    rootfile_dir_p = rootfile_dir + "deep";

    util::FindMatchingFiles(histfile_dir_p,rootfile_dir_p,partialName_p,histFileNames_p,rootFileNames_p,false);

  }
  
  if(kine==4&&mag==30&&alt_infile){
    histfile_dir_p = histfile_dir_alt + "deep";
    rootfile_dir_p = rootfile_dir_alt + "deep";
    util::FindMatchingFiles(histfile_dir_p,rootfile_dir_p,partialName_p,histFileNames_p_alt,rootFileNames_p_alt,false);

  }
    
  //Check neutron root/hist files

  std::vector<std::string> rootFileNames_n, histFileNames_n;
  std::vector<std::string> rootFileNames_n_alt, histFileNames_n_alt;
  std::string histfile_dir_n;
  std::string rootfile_dir_n;

  if( rtype.compare("jboyd")==0 )
    util::FindMatchingFiles(histfile_dir,rootfile_dir,partialName_n,histFileNames_n,rootFileNames_n,true);
  else if( rtype.compare("pdatta")==0 ){
    util::FindMatchingFiles(histfile_dir,rootfile_dir,partialName_n,histFileNames_n,rootFileNames_n,false);
  }else{
    histfile_dir_n = histfile_dir + "deen";
    rootfile_dir_n = rootfile_dir + "deen";
    util::FindMatchingFiles(histfile_dir_n,rootfile_dir_n,partialName_n,histFileNames_n,rootFileNames_n,false);

  }

  if(kine==4&&mag==30&&alt_infile){
    histfile_dir_n = histfile_dir_alt + "deen";
    rootfile_dir_n = rootfile_dir_alt + "deen";
    util::FindMatchingFiles(histfile_dir_n,rootfile_dir_n,partialName_n,histFileNames_n_alt,rootFileNames_n_alt,false);

  }

  //Check to make sure the number of rootfiles and number of histfiles match on protons and neutrons
  if( rootFileNames_p.size() != histFileNames_p.size() || 
      rootFileNames_n.size() != histFileNames_n.size() ||
      rootFileNames_p_alt.size() != histFileNames_p_alt.size() || 
      rootFileNames_n_alt.size() != histFileNames_n_alt.size())
    cerr << "ERROR: FindMatchingFiles() failure, vector sizes mismatch" << endl;

  //if two sets, add together
  if(kine==4&&mag==30&&alt_infile){

    rootFileNames_p.insert(rootFileNames_p.end(), rootFileNames_p_alt.begin(), rootFileNames_p_alt.end());
    rootFileNames_n.insert(rootFileNames_n.end(), rootFileNames_n_alt.begin(), rootFileNames_n_alt.end());
    histFileNames_p.insert(histFileNames_p.end(), histFileNames_p_alt.begin(), histFileNames_p_alt.end());
    histFileNames_n.insert(histFileNames_n.end(), histFileNames_n_alt.begin(), histFileNames_n_alt.end());

  }

  //Throw away jobs for which there does not exist both a proton and neutron rootfile
  if(sync_jobs)
    util::synchronizeJobNumbers(rootFileNames_p,rootFileNames_n);

  //Check to be sure same number of stats between p and n
  int pNfiles = rootFileNames_p.size();
  int nNfiles = rootFileNames_n.size();
  cout << endl << "Number of proton files: " << pNfiles << ", Number of neutron files: " << nNfiles << endl;
  if(rootFileNames_p.size() != rootFileNames_n.size()){
    cout << endl << "WARNING: Number of proton and neutron files loaded do not match. Check MC file directory" << endl << endl;
    if(sync_jobs)
      return;
  }

  cout << endl << endl << "CHECKING ROOTFILE PATHS..." << endl;
  for( size_t f=0; f<rootFileNames_p.size(); ++f )
    cout << rootFileNames_p[f] << endl;

  cout << endl << endl;


  ////////////
  //HISTOGRAMS

  //basic
  TH2D *hxy_nocut = new TH2D("hxy_nocut","HCal x vs HCal y, no cut; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );

  //E-arm
  
  TH1D *hW2_nocut = new TH1D( "hW2_nocut", "W^{2}, no cut; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_cut = new TH1D( "hW2_cut", "W^{2}, accmatch/global/atime cuts; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  TH2D *hxy_exp_n = new TH2D("hxy_exp_n","HCal x vs y expected from e', elastic cuts neutron; y_{HCAL} (m); x_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_exp_p = new TH2D("hxy_exp_p","HCal x vs y expected from e', elastic cuts proton; y_{HCAL} (m); x_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_exp_n_fid = new TH2D("hxy_exp_n_fid","HCal x vs y expected from e', all cuts neutron; y_{HCAL} (m); x_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_exp_p_fid = new TH2D("hxy_exp_p_fid","HCal x vs y expected from e', all cuts proton; y_{HCAL} (m); x_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );


  //H-arm
  TH2D *hdxdy_nocut = new TH2D("hdxdy_nocut","HCal dxdy, no cut; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hdxdy_cut = new TH2D("hdxdy_cut","HCal dxdy, all cuts; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );

  TH2D *hdxdy_nocut_p = new TH2D("hdxdy_nocut_p","proton dxdy, no cut; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hdxdy_cut_p = new TH2D("hdxdy_cut_p","proton dxdy, all cuts; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );

  TH2D *hdxdy_nocut_n = new TH2D("hdxdy_nocut_n","neutron dxdy, no cut; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hdxdy_cut_n = new TH2D("hdxdy_cut_n","neutron dxdy, all cuts; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );

  TH1D *hdx_nocut = new TH1D( "hdx_nocut", "dx, no cut;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);
  TH1D *hdx_cut = new TH1D( "hdx_cut", "dx, e-arm cut, no dy;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);
  TH1D *hdx_cut_nofid = new TH1D( "hdx_cut_nofid", "dx, all cuts sans fiducial;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);
  TH1D *hdx_cut_failfid = new TH1D( "hdx_cut_failfid", "dx, all cuts fail fiducial;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);

  TH1D *hdx_nocut_p = new TH1D( "hdx_nocut_p", "proton dx, no cut;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);
  TH1D *hdx_cut_p = new TH1D( "hdx_cut_p", "proton dx, e-arm cut, no dy;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);

  TH1D *hdx_nocut_n = new TH1D( "hdx_nocut_n", "neutron dx, no cut;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);
  TH1D *hdx_cut_n = new TH1D( "hdx_cut_n", "neutron dx, e-arm cut, no dy;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);

  std::string tar = "LD2";

  //set up configuration and tune objects to load analysis parameters
  SBSconfig config(kine,mag);

  //Obtain configuration pars from config file
  Double_t hcaltheta = config.GetHCALtheta_rad();
  Double_t hcaldist = config.GetHCALdist();
  Double_t sbsdist = config.GetSBSdist();
  Double_t bbthr = config.GetBBtheta_rad(); //in radians
  Double_t ebeam = config.GetEbeam();

  //Set up hcal active area with bounds that match database on pass
  vector<Double_t> hcalaa = cut::hcalaa_data(1,1);

  //SBStune *tune = new SBStune(kine,mag);
  SBStune tune(kine,mag);
    
  //Reporting. tar should always equal curtar as categorized by good run list
  std::cout << "Settings are.." << endl;
  cout << endl << config << endl;
  cout << endl << tune << endl;

  //Obtain cuts from tune class
  std:string gcut   = tune.Getglobcut_earm();
  Double_t W2mean   = tune.GetW2mean();
  Double_t W2sig    = tune.GetW2sig();
  Double_t dx0_n    = tune.Getdx0_n();
  Double_t dx0_p    = tune.Getdx0_p();
  Double_t dy0      = tune.Getdy0();
  Double_t dxsig_n  = tune.Getdxsig_n();
  Double_t dxsig_p  = tune.Getdxsig_p();
  Double_t dysig    = tune.Getdysig();
  Double_t atime0   = tune.Getatime0();
  Double_t atimesig = tune.Getatimesig();
  Double_t atimediff0   = tune.Getatimediff0();
  Double_t atimediffsig = tune.Getatimediffsig();

  Double_t W2min = W2mean - 2*W2sig;
  Double_t W2max = W2mean + 2*W2sig;

  //Set minimum energy for expected recoil nucleon to be considered
  Double_t hminE = 0.0;  //0.05 GeV or nothing
    
  //minimize reporting
  std::string nuc = "none";
  std::string snuc = nucleon;
  //std::string load_file = "";
  int nnuc = 2; //two nucleons, proton and neutron
  int nfiles = 0;

  TChain *C = nullptr;

  //Main loop over nucleons (r==0, proton; r==1, neutron)
  for (int r=0; r<nnuc; r++) {
    
    //validate nucleon argument
    if( snuc.compare("proton")==0 && r==0 )
      continue;
    else if( snuc.compare("neutron")==0 && r==1 )
      continue;
    else if( snuc.compare("both")!=0 ){
      cerr << "ERROR: Enter a valid nucleon. Valid arguments are {proton,neutron,both}" << endl;
      return;
    }
    
    //update current nucleon
    if( r==0 ){
      nuc = "p";
      nfiles = rootFileNames_p.size();
    }else if( r==1 ){
      nuc = "n";
      nfiles = rootFileNames_n.size();
    }

    //loop over simulation files
    for( int f=0; f<nfiles; ++f ){

      std::string load_file;
      std::string hist_file;
      
      if(r==0){
	load_file = rootFileNames_p[f];
	hist_file = histFileNames_p[f];
	
      }else if(r==1){
	load_file = rootFileNames_n[f];
	hist_file = histFileNames_n[f];
      }
      
      //Get metadata for MC weights
      Int_t Ntried;
      Double_t luminosity;
      Double_t genvol;
      if( norm_override ){
	cout << "Normalization override enabled. ";
	Ntried = (Int_t)Ntried_override;
	luminosity = luminosity_override;
	genvol = genvol_override;
      }else if( rtype.compare("jboyd")==0 ){
	Ntried = (Int_t)util::searchSimcHistFile("Ntried", hist_file);
	luminosity = util::searchSimcHistFile("luminosity", hist_file);
	genvol = util::searchSimcHistFile("genvol", hist_file);
      }else{
	Ntried = (Int_t)util::searchSimpleCSVForValue(hist_file, "N_tries");
	luminosity = util::searchSimpleCSVForValue(hist_file, "Luminosity_s-1_cm-2");
	genvol = util::searchSimpleCSVForValue(hist_file, "Generation_Volume");
      }

      std::cout << "For this run: Ntried=" << Ntried << " luminosity=" << luminosity << " genvol=" << genvol << endl; 
      
      if( !norm_override && (Ntried==-1 || luminosity==-1 || genvol==-1) ){
	cout << "Failed to read normalization data, skipping file.." << endl;
	continue;
      }

      //first attempt to resolve segmentation fault on large data sets
      if (C != nullptr) {
	delete C;
      }

      C = new TChain("T"); 
      
      //add the right replayed mc output file
      C->Add(load_file.c_str());
    
      long nevents = C->GetEntries();

      cout << "Opening " << nuc << " MC file at " << load_file << " with " << nevents << " entries for analysis." << endl;

      // setting up ROOT tree branch addresses
      C->SetBranchStatus("*",0);    

      C->SetBranchStatus("bb.sh.nclus",1);

      // HCal general
      Double_t hcalid, hcale, hcalx, hcaly, hcalr, hcalc, hcaltdc, hcalatime;
      std::vector<std::string> hcalvar = {"idblk","e","x","y","rowblk","colblk","tdctimeblk","atimeblk"};
      std::vector<void*> hcalvarlink = {&hcalid,&hcale,&hcalx,&hcaly,&hcalr,&hcalc,&hcaltdc,&hcalatime};
      rvars::setbranch(C, "sbs.hcal", hcalvar, hcalvarlink);

      // HCal cluster branches (primary block information for id, tdc, and atime)
      Double_t hcalcid[econst::maxclus], hcalce[econst::maxclus], hcalcx[econst::maxclus], hcalcy[econst::maxclus], hcalctdctime[econst::maxclus], hcalcatime[econst::maxclus];
      Int_t Nhcalcid;
      std::vector<std::string> hcalcvar = {"id","e","x","y","tdctime","atime","id"};
      std::vector<void*> hcalcvarlink = {&hcalcid,&hcalce,&hcalcx,&hcalcy,&hcalctdctime,&hcalcatime,&Nhcalcid};
      rvars::setbranch(C, "sbs.hcal.clus", hcalcvar, hcalcvarlink,6);

      // Monte Carlo variables
      Double_t mcweight;
      std::vector<std::string> mcvar = {"simc_Weight"};
      std::vector<void*> mcvarlink = {&mcweight};
      rvars::setbranch(C, "MC", mcvar, mcvarlink);
   
      // BBCal shower timing
      Int_t Natime_sh; 
      Double_t atime_sh[econst::maxclus];
      std::vector<std::string> shvar = {"atime","atime"};
      std::vector<void*> shvarlink = {&Natime_sh,&atime_sh};
      rvars::setbranch(C, "bb.sh.clus", shvar, shvarlink, 0);  
    
      // track branches
      Double_t ntrack, p[econst::maxtrack],px[econst::maxtrack],py[econst::maxtrack],pz[econst::maxtrack],xtr[econst::maxtrack],ytr[econst::maxtrack],thtr[econst::maxtrack],phtr[econst::maxtrack];
      Double_t vx[econst::maxtrack],vy[econst::maxtrack],vz[econst::maxtrack];
      Double_t xtgt[econst::maxtrack],ytgt[econst::maxtrack],thtgt[econst::maxtrack],phtgt[econst::maxtrack];
      std::vector<std::string> trvar = {"n","p","px","py","pz","x","y","th","ph","vx","vy","vz","tg_x","tg_y","tg_th","tg_ph"};
      std::vector<void*> trvarlink = {&ntrack,&p,&px,&py,&pz,&xtr,&ytr,&thtr,&phtr,&vx,&vy,&vz,&xtgt,&ytgt,&thtgt,&phtgt};
      rvars::setbranch(C,"bb.tr",trvar,trvarlink);
   
      // ekine branches
      Double_t ekineQ2, ekineW2, ekineeps, ekinenu, ekineqx, ekineqy, ekineqz;
      std::vector<std::string> ekinevar = {"Q2","W2","epsilon","nu","q_x","q_y","q_z"};
      std::vector<void*> ekinevarlink = {&ekineQ2,&ekineW2,&ekineeps,&ekinenu,&ekineqx,&ekineqy,&ekineqz};
      rvars::setbranch(C, "e.kine", ekinevar, ekinevarlink);
    
      // fEvtHdr branches
      UInt_t rnum, gevnum, trigbits;
      std::vector<std::string> evhdrvar = {"fRun","fEvtNum","fTrigBits"};
      std::vector<void*> evhdrlink = {&rnum,&gevnum,&trigbits};
      rvars::setbranch(C,"fEvtHdr",evhdrvar,evhdrlink);
    
      // globalcut branches
      C->SetBranchStatus("bb.gem.track.nhits", 1);
      C->SetBranchStatus("bb.etot_over_p", 1);
      C->SetBranchStatus("bb.tr.n", 1);
      C->SetBranchStatus("bb.ps.e", 1);
      C->SetBranchStatus("bb.sh.e", 1);

      TCut GCut = gcut.c_str();
      TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );

      // get experimental quantities by run
      std::cout << "Uncorrected average beam energy on " << tar << " for run: " << ebeam << std::endl;
      //set up hcal coordinate system with hcal angle wrt exit beamline
      vector<TVector3> hcalaxes; vars::sethcalaxes( hcaltheta, hcalaxes );
      //TVector3 hcalorigin = hcaldist*hcalaxes[2] + econst::hcalvoff*hcalaxes[0];
      TVector3 hcalorigin = hcaldist*hcalaxes[2];
      Double_t BdL = econst::sbsmaxfield * econst::sbsdipolegap * (mag/100); //scale crudely by magnetic field
      Double_t Eloss_outgoing = econst::celldiameter/2.0/sin(bbthr) * econst::lh2tarrho * econst::lh2dEdx;
      //vector<Double_t> hcalaa = cut::hcalaa_data(1,1);
      //vector<Double_t> hcalaa = cut::hcalaa_data_alt(1,1); //ostensibly for pass 0

      // event indices
      long nevent = 0; 
      Int_t treenum = 0, currenttreenum = 0;
    
      while (C->GetEntry(nevent++)) {
      
	std::cout << "Processing " << nuc << " file " << f << "/" << nfiles << ", event " << nevent << "/" << nevents << "\r";
	std::cout.flush();
      
	///////
	//Single-loop globalcut method. Save pass/fail for output tree.
	currenttreenum = C->GetTreeNumber();
	if( nevent == 1 || currenttreenum != treenum ){
	  treenum = currenttreenum; 
	  GlobalCut->UpdateFormulaLeaves();
	}
	bool failedglobal = GlobalCut->EvalInstance(0) == 0;
      
	///////
	//E-arm physics calculations
	//correct beam energy with vertex information
	Double_t ebeam_c = vars::ebeam_c( ebeam, vz[0], target.c_str() );
	TVector3 vertex( 0., 0., vz[0] );
      
	//reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
	Double_t precon = p[0] + Eloss_outgoing;
      
	//set up four-momenta with some empty for various calculation methods
	TLorentzVector pbeam( 0., 0., ebeam_c, ebeam_c ); //beam momentum
	//TLorentzVector pe( px[0], py[0], pz[0], p[0] ); //e' plvect
	TLorentzVector pe( precon*px[0]/p[0], precon*py[0]/p[0], precon*pz[0]/p[0], precon ); //e' recon plvect
	TLorentzVector ptarg; vars::setPN(nuc,ptarg); //target momentum
	TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
	TVector3 qv = q.Vect();
	TLorentzVector pN; //N' momentum
	TVector3 pNhat; //Unit N' 3-vector
      
	//simple calculations for e' and N'
	Double_t etheta = vars::etheta(pe); 
	Double_t ephi = vars::ephi(pe);
	Double_t pcent = vars::pcentral(ebeam,etheta,nuc); //e' p reconstructed by angles
	Double_t phNexp = ephi + physconst::pi;
	Double_t Q2, W2, nu, thNexp, pNexp;
	Double_t ebeam_o = vars::ebeam_o( ebeam_c, etheta, target.c_str() ); //Second energy correction accounting for energy loss leaving target

	//Using default e' momentum calculation method which reconstructs via track angles
	nu = pbeam.E() - pcent;
	pNexp = vars::pN_expect( nu, nuc );
	thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	pNhat = vars::pNhat_track( thNexp, phNexp );
	pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	Q2 = vars::Q2( pbeam.E(), pe.E(), etheta );
	W2 = vars::W2( pbeam.E(), pe.E(), Q2, nuc );

	/////////////////////
	//W2 elastic cut bool
	bool failedW2 = W2<W2min || W2>W2max;

	///////
	//HCal active area cut (acceptance matching). Save pass/fail for output tree.
	vector<Double_t> xyhcalexp; vars::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
      
	//Get expected proton deflection for thetapq (this needs work to be useful)
	Double_t protdeflect = tan( 0.3 * BdL / qv.Mag() ) * (hcaldist - (sbsdist + econst::sbsdipolegap/2.0) );
      
	//Get some primary cluster primary block variables
	Double_t dx = hcalx - xyhcalexp[0];
	Double_t dy = hcaly - xyhcalexp[1];

	//check the acceptance limits
	hxy_nocut->Fill(hcaly,hcalx);

	//Primary cluster atime cut
	bool failedpclusatime =  abs( hcalcatime[0] - atime0 ) > coin_sigma_factor*atimesig;

	//caculate final weight for this event
	Double_t final_mc_weight = mcweight*luminosity*genvol/Ntried;

	//Fill nocut histograms for both
	hW2_nocut->Fill(W2,final_mc_weight);
	hdxdy_nocut->Fill(dy,dx,final_mc_weight);
	hdx_nocut->Fill(dx,final_mc_weight);

	if( nuc.compare("p")==0 ){
	  hdxdy_nocut_p->Fill(dy,dx,final_mc_weight);
	  hdx_nocut_p->Fill(dx,final_mc_weight);
	}

	if( nuc.compare("n")==0 ){
	  hdxdy_nocut_n->Fill(dy,dx,final_mc_weight);
	  hdx_nocut_n->Fill(dx,final_mc_weight);
	}

	//H-arm fiducial cuts
	bool hcalON = cut::hcalaaON(hcalx,hcaly,hcalaa);
	vector<Double_t> fid = cut::hcalfid(dxsig_p,dysig,hcalaa);
	bool passed_fid = cut::hcalfidIN(xyhcalexp[0],xyhcalexp[1],dx0_p,fid);

	//E-arm cuts
	if( failedglobal || failedW2 )
	  continue;
    
	hxy_exp_n->Fill(xyhcalexp[1],xyhcalexp[0]);
	hxy_exp_p->Fill(xyhcalexp[1],xyhcalexp[0]+dx0_p); //expected proton position from average difference between dxdy neutron and proton spots

	if(passed_fid){
	  hxy_exp_n_fid->Fill(xyhcalexp[1],xyhcalexp[0]);
	  hxy_exp_p_fid->Fill(xyhcalexp[1],xyhcalexp[0]+dx0_p); //expected proton position from average difference between dxdy neutron and proton spots
	}

	//dy cut
	bool faileddy = abs( dy - dy0 ) > 3*dysig;
	
	if( !faileddy && hcalON ){
	  hdx_cut_nofid->Fill(dx,final_mc_weight); //primary dx histo, no fiducial cut
	  if(!passed_fid)
	    hdx_cut_failfid->Fill(dx,final_mc_weight); //primary dx histo, 
	}

	if( !passed_fid || !hcalON || faileddy)
	  continue;

	//Fill cut histograms for both
	hW2_cut->Fill(W2,final_mc_weight);
	hdxdy_cut->Fill(dy,dx,final_mc_weight);
	hdx_cut->Fill(dx,final_mc_weight);

	if( nuc.compare("p")==0 ){
	  hdxdy_cut_p->Fill(dy,dx,final_mc_weight);
	  hdx_cut_p->Fill(dx,final_mc_weight);
	}

	if( nuc.compare("n")==0 ){
	  hdxdy_cut_n->Fill(dy,dx,final_mc_weight);
	  hdx_cut_n->Fill(dx,final_mc_weight);
	}

      } //end event loop

    } //end file loop

  } //end nucleon loop
   
   //Make fiducial check on neutron hypothesis
  vector<Double_t> fid = cut::hcalfid(dxsig_p,dysig,hcalaa);
  Double_t hcalx_ft = fid[2];
  Double_t hcalx_fb = fid[3];
  Double_t hcaly_fr = fid[0];
  Double_t hcaly_fl = fid[1];

  //Make Fiducial Safety Margin TLines
  TLine* l1_f = new TLine(hcalx_ft, hcaly_fr, hcalx_ft, hcaly_fl);
  l1_f->SetLineColor(kGreen);
  TLine* l2_f = new TLine(hcalx_fb, hcaly_fr, hcalx_fb, hcaly_fl);
  l2_f->SetLineColor(kGreen);
  TLine* l3_f = new TLine(hcalx_ft, hcaly_fr, hcalx_fb, hcaly_fr);
  l3_f->SetLineColor(kGreen);
  TLine* l4_f = new TLine(hcalx_ft, hcaly_fl, hcalx_fb, hcaly_fl);
  l4_f->SetLineColor(kGreen);

  //Make Fiducial Line
  // TLine* l1_fl = new TLine(hcalx_t-2*dysig, hcaly_fl, hcalx_t-2*dysig, hcaly_fl+dx0_p);
  // l1_fl->SetLineColor(kMagenta);

  TCanvas *c0 = new TCanvas("c0","HCal Fiducial cuts, neutron hypothesis",1200,500);
  c0->Divide(2,1);

  c0->cd(1);
  c0->SetLogz();
  hxy_exp_n->Draw("colz");
  l1_f->Draw();
  l2_f->Draw();
  l3_f->Draw();
  l4_f->Draw();

  auto leg0a = new TLegend(0.1,0.8,0.5,0.9);
  leg0a->AddEntry( l1_f, "HCal Safety Margin", "l");
  leg0a->Draw();

  c0->Update();

  c0->cd(2);
  c0->SetLogz();
  hxy_exp_n_fid->Draw("colz");
  l1_f->Draw();
  l2_f->Draw();
  l3_f->Draw();
  l4_f->Draw();
  //l1_fl->Draw();

  auto leg0b = new TLegend(0.1,0.8,0.5,0.9);
  leg0b->AddEntry( l1_f, "HCal Safety Margin", "l");
  leg0b->Draw();

  c0->Update();
  c0->Write();

  //Make fiducial check on proton hypothesis

  TCanvas *c1 = new TCanvas("c1","HCal Fiducial cuts, proton hypothesis",1200,500);
  gStyle->SetPalette(53);
  c1->Divide(2,1);

  c1->cd(1);
  c1->SetLogz();
  hxy_exp_p->Draw("colz");
  l1_f->Draw();
  l2_f->Draw();
  l3_f->Draw();
  l4_f->Draw();

  auto leg1a = new TLegend(0.1,0.8,0.5,0.9);
  leg1a->AddEntry( l1_f, "HCal Safety Margin", "l");
  leg1a->Draw();

  c1->Update();

  c1->cd(2);
  c1->SetLogz();
  hxy_exp_p_fid->Draw("colz");
  l1_f->Draw();
  l2_f->Draw();
  l3_f->Draw();
  l4_f->Draw();
  //l1_fl->Draw();

  auto leg1b = new TLegend(0.1,0.8,0.5,0.9);
  leg1b->AddEntry( l1_f, "HCal Safety Margin", "l");
  leg1b->Draw();

  c1->Update();
  c1->Write();

  fout->Write();
  
  cout << endl << "Analysis complete. Outfile located at " << fout_path << endl;

  // Send time efficiency report to console
  st->Stop();
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
