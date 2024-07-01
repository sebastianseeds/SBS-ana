//SSeeds 9.1.23 Script to extract the p/n yield ratio via data/MC comparison. Intent is to loop over all data from a given set/field setting, apply elastic cuts, and build a dx histogram. From this point, read in MC (with RC) distributions for quasi-elastic protons and neutrons (independently), fit these distributions with several distributions to get a best fit, then compose a sum of these fits allowing for a scaling parameter which maps to the total MC yields. Next, do the same with g4sbs inelastic generator and add this to the total fit function (sum) and apply the now three floating parameter fit to the data and check residuals and chi-square. Varying the fit function may yield better results, but the ratio will be extracted. See gmn_calc.C for application of this ratio to gmn extraction.
//sseeds update 9.23.23 Parsed cluster sorting and e' momentum strategies. Renamed and refocused script to produce output histograms for later analysis.
//sseeds update 9.24.23 Debug and prevent memory leaks
//sseeds update 11.14.23 remove assignScore from here and add to util function list
//sseeds update 12.19.23 added option to parse globalcuts and output analysis tree to examine effects of cuts on dx

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

//Manual override info
bool debug = false;

//MAIN
void data_elastic( Int_t kine=9, Int_t magset=70, Int_t pass=2, bool systematicInfo=true )
{

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/gmn.json");

  std::string rootfile_word = "";
  if(pass>1)
    rootfile_word = Form("_pass%d",pass);

  std::string rootfile_dir = jmgr->GetValueFromSubKey_str( Form("rootfile_dir%s",rootfile_word.c_str()), Form("sbs%d",kine) );
  Double_t W2fitmax = jmgr->GetValueFromSubKey<Double_t>( "W2fitmax", Form("sbs%d",kine) );
  Double_t binfac = jmgr->GetValueFromSubKey<Double_t>( "binfac", Form("sbs%d",kine) );
  Double_t hbinfac = jmgr->GetValueFromSubKey<Double_t>( "hbinfac", Form("sbs%d",kine) ); //gets to 7cm HCal pos res
  std::string target = jmgr->GetValueFromKey_str( "target" ); //Must always be LD2 for GMn analysis
  Int_t epm = jmgr->GetValueFromSubKey<Int_t>( "epm", Form("sbs%d",kine) );
  Double_t minE = jmgr->GetValueFromSubKey<Double_t>( "minE", Form("sbs%d",kine) );
  Int_t cluster_idx = jmgr->GetValueFromSubKey<Int_t>( "cluster_idx", Form("sbs%d",kine) );
  //Int_t pass = jmgr->GetValueFromSubKey<Int_t>( "pass", Form("sbs%d",kine) );
  Double_t ebound_l = jmgr->GetValueFromSubKey<Double_t>( "ebound_l", Form("sbs%d",kine) );
  Double_t ebound_h = jmgr->GetValueFromSubKey<Double_t>( "ebound_h", Form("sbs%d",kine) );
  Double_t coin_sigma_factor = jmgr->GetValueFromSubKey<Double_t>( "coin_sigma_factor", Form("sbs%d",kine) );
  vector<Double_t> coin_profile;
  jmgr->GetVectorFromSubKey<Double_t>("coin_profile",Form("sbs%d",kine),coin_profile);
  Double_t hcal_v_offset = jmgr->GetValueFromSubKey<Double_t>( Form("hcal_offset%s",rootfile_word.c_str()), Form("sbs%d",kine) );

  //Necessary for loop over directories in sbs8
  std::vector<TString> directories = {
    rootfile_dir + "/SBS0percent",
    rootfile_dir + "/SBS100percent",
    rootfile_dir + "/SBS50percent",
    rootfile_dir + "/SBS70percent_part1",
    rootfile_dir + "/SBS70percent_part2",
    rootfile_dir + "/SBS70percent_part3",
    rootfile_dir + "/SBS70percent_part4"
  };

  std::cout << "Loaded HCal vertical offset: " << hcal_v_offset << std::endl;

  //set up configuration and tune objects to load analysis parameters
  SBSconfig config(kine,magset);

  //Obtain configuration pars from config file
  Double_t hcaltheta = config.GetHCALtheta_rad();
  Double_t hcaldist = config.GetHCALdist();
  Double_t sbsdist = config.GetSBSdist();
  Double_t bbthr = config.GetBBtheta_rad(); //in radians

  //Set up hcal active area with bounds that match database on pass
  vector<Double_t> hcalaa;
  if(pass<1)
    hcalaa = cut::hcalaa_data_alt(1,1);
  else
    hcalaa = cut::hcalaa_mc(1,1); //verified 2.10.24

  //SBStune *tune = new SBStune(kine,mag);
  SBStune tune(kine,magset);
    
  //Reporting. tar should always equal curtar as categorized by good run list
  std::cout << "Settings are.." << endl;
  cout << endl << config << endl;
  cout << endl << tune << endl;

  //Obtain cuts from tune class
  std::string gcut   = tune.Getglobcut();
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

  Double_t W2min = W2mean - 3*W2sig;
  Double_t W2max = W2mean + 3*W2sig;

  cout << "Loaded globalcut: " << gcut << endl;

  //set up default parameters for all analysis
  std::string runsheet_dir = "/w/halla-scshelf2102/sbs/seeds/ana/data"; //unique to my environment for now
  Int_t nruns = -1; //Always analyze all available runs
  Int_t verb = 0; //Don't print diagnostic info by default
  Double_t hcalfit_l = econst::hcalposXi_mc; //lower fit/bin limit for hcal dx plots (m)
  Double_t hcalfit_h = econst::hcalposXf_mc; //upper fit/bin limit for hcal dx plots (m) 
  Double_t harmrange = econst::hcal_vrange; //Full range of hcal dx plots (m)

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<crun> corun; 
  util::ReadRunList(runsheet_dir,nruns,kine,target.c_str(),pass,verb,corun); //modifies nruns to be very large when -1

  //set up output files
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fout_path;
  if(systematicInfo)
    fout_path = outdir_path + Form("/gmn_analysis/gmn_elastic_fid_sys_out_sbs%d_mag%d_clusteridx%d_epm%d.root",kine,magset,cluster_idx,epm);
  else
    fout_path = outdir_path + Form("/gmn_analysis/gmn_elastic_fid_out_sbs%d_mag%d_clusteridx%d_epm%d.root",kine,magset,cluster_idx,epm);

  TFile *fout = new TFile( fout_path.c_str(), "RECREATE" );

  ////////////
  //HISTOGRAMS

  //basic H-arm
  TH2D *hxy_nocut = new TH2D("hxy_nocut","HCal x vs HCal y, no cut; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_aacut = new TH2D("hxy_aacut","HCal x vs HCal y, acceptance matching cut no sigma; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_acccut = new TH2D("hxy_acccut","HCal x vs HCal y, acceptance matching cut; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );

  //E-arm
  TH1D *hW2_nocut = new TH1D( "hW2_nocut", "W^{2}, no cut; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_cut = new TH1D( "hW2_cut", "W^{2}, accmatch/global/atime cuts; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  TH2D *hxy_exp_n = new TH2D("hxy_exp_n","HCal x vs y expected from e', elastic cuts neutron; y_{HCAL} (m); x_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_exp_p = new TH2D("hxy_exp_p","HCal x vs y expected from e', elastic cuts proton; y_{HCAL} (m); x_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_exp_n_fid = new TH2D("hxy_exp_n_fid","HCal x vs y expected from e', all cuts neutron; y_{HCAL} (m); x_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_exp_p_fid = new TH2D("hxy_exp_p_fid","HCal x vs y expected from e', all cuts proton; y_{HCAL} (m); x_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );

  //Both arms
  TH2D *hdxdy_nocut = new TH2D("hdxdy_nocut","dxdy, no cut; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hdxdy_cut = new TH2D("hdxdy_cut","dxdy, all cuts; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hdxdy_cut_nofid = new TH2D("hdxdy_cut_nofid","dxdy, all cuts sans fiducial; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hdxvE = new TH2D("hdxvE","HCal dx vs HCal E, all cuts; E_{HCAL} (GeV); dx_{HCAL} (m)", 400, 0.0, 4.0, 600, -3.0, 3.0 );
  TH1D *hdx_nocut = new TH1D( "hdx_nocut", "dx, earm cut only (W2 and global);x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);
  TH1D *hdx_cut = new TH1D( "hdx_cut", "dx, all cuts;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);
  TH1D *hdx_cut_nofid = new TH1D( "hdx_cut_nofid", "dx, all cuts sans fiducial;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);
  TH1D *hdx_cut_failfid = new TH1D( "hdx_cut_failfid", "dx, all cuts fail fiducial;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);
  TH1D *hdx_cut_globalonly = new TH1D( "hdx_cut_globalonly", "dx, global cuts only;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);
  TH1D *hdx_cut_W2only = new TH1D( "hdx_cut_W2only", "dx, W2 cuts only;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);
  TH1D *hdx_cut_earm = new TH1D( "hdx_cut_earm", "dx, earm cuts only;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);

  TH1D *hcoin_pclus = new TH1D( "hcoin_pclus", "HCal ADCt - BBCal ADCt, pclus, no cut; ns", 200, 0, 100 );
  TH1D *hcoin_pclus_cut = new TH1D( "hcoin_pclus_cut", "HCal ADCt - BBCal ADCt, pclus; ns", 200, 0, 100 );
  TH1D *hcoin = new TH1D( "hcoin", "HCal ADCt - BBCal ADCt, no cut; ns", 200, 0, 100 );
  TH1D *hcoin_cut = new TH1D( "hcoin_cut", "HCal ADCt - BBCal ADCt; ns", 200, 0, 100 );

  TH2D *hHcalXY = new TH2D( "hHcalXY", "HCal X vs HCal Y; Y (m); X (m)",600,-3,3,600,-3,3);
  TH2D *hHcalXY_allcuts = new TH2D( "hHcalXY_allcuts", "HCal X vs HCal Y, all cuts, best cluster; Y (m); X (m)",600,-3,3,600,-3,3);

  //general
  TH1D *hMott_cs = new TH1D( "hMott_cs", "Mott Cross Section, no cut; (GeV/c)^{-2}", 200, 0, 0.0001 );

  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;

  // create output tree
  TTree *P = new TTree("P","Analysis Tree"); 

  // uncut output tree vars
  Double_t dx_out;
  Double_t dy_out;
  Double_t xexp_out;
  Double_t yexp_out;
  Double_t hcalx_out;
  Double_t hcaly_out;
  Double_t W2_out;
  Double_t Q2_out;
  Double_t mott_out;
  Double_t hcalE_out;
  Double_t bbEtot_out;
  Double_t bbshE_out;
  Double_t bbpsE_out;
  Double_t hcalatime_out;
  Double_t bbshatime_out;
  Double_t bbpsatime_out;
  Double_t bbgemNhits_out;
  Double_t bbtrx_out;
  Double_t bbtry_out;
  Double_t bbtrp_out;
  Double_t coinP1_out;
  Double_t coinP2_out;
  Double_t dy0_out;
  Double_t dx0_p_out;
  Double_t dx0_n_out;
  Double_t dysig_out;
  Double_t dxsig_p_out;
  Double_t dxsig_n_out;
  Double_t W2mean_out;
  Double_t W2sig_out;

  Int_t hcalnclus_out;
  Int_t hcalnblk_out;
  Int_t bbshnclus_out;
  Int_t bbshnblk_out;
  Int_t bbpsnclus_out;
  Int_t bbpsnblk_out;
  Int_t bbtrN_out;
  Int_t passedfid_out;  
  Int_t run_out;
  Int_t mag_out;

  // set new output tree branches
  P->Branch( "dx", &dx_out, "dx/D" );
  P->Branch( "dy", &dy_out, "dy/D" );
  P->Branch( "xexp", &xexp_out, "xexp/D" );
  P->Branch( "yexp", &yexp_out, "yexp/D" );
  P->Branch( "hcalx", &hcalx_out, "hcalx/D" );
  P->Branch( "hcaly", &hcaly_out, "hcaly/D" );
  P->Branch( "W2", &W2_out, "W2/D" );
  P->Branch( "Q2", &Q2_out, "Q2/D" );
  P->Branch( "mott", &mott_out, "mott/D" );
  P->Branch( "hcalE", &hcalE_out, "hcalE/D" );
  P->Branch( "bbEtot", &bbEtot_out, "bbEtot/D" );
  P->Branch( "bbshE", &bbshE_out, "bbshE/D" );
  P->Branch( "bbpsE", &bbpsE_out, "bbpsE/D" );
  P->Branch( "hcalatime", &hcalatime_out, "hcalatime/D" );
  P->Branch( "bbshatime", &bbshatime_out, "bbshatime/D" );
  P->Branch( "bbpsatime", &bbpsatime_out, "bbpsatime/D" );
  P->Branch( "bbgemNhits", &bbgemNhits_out, "bbgemNhits/D" );
  P->Branch( "bbtrx", &bbtrx_out, "bbtrx/D" );
  P->Branch( "bbtry", &bbtry_out, "bbtry/D" );
  P->Branch( "bbtrp", &bbtrp_out, "bbtrp/D" );
  P->Branch( "coinP1", &coinP1_out, "coinP1/D" );
  P->Branch( "coinP2", &coinP2_out, "coinP2/D" );
  P->Branch( "dy0", &dy0_out, "dy0/D" );
  P->Branch( "dx0_p", &dx0_p_out, "dx0_p/D" );
  P->Branch( "dx0_n", &dx0_n_out, "dx0_n/D" );
  P->Branch( "dysig", &dysig_out, "dysig/D" );
  P->Branch( "dxsig_p", &dxsig_p_out, "dxsig_p/D" );
  P->Branch( "dxsig_n", &dxsig_n_out, "dxsig_n/D" );
  P->Branch( "W2mean", &W2mean_out, "W2mean/D" );
  P->Branch( "W2sig", &W2sig_out, "W2sig/D" );

  P->Branch( "hcalnclus", &hcalnclus_out, "hcalnclus/I" );
  P->Branch( "hcalnblk", &hcalnblk_out, "hcalnblk/I" );
  P->Branch( "bbshnclus", &bbshnclus_out, "bbshnclus/I" );
  P->Branch( "bbpsnclus", &bbpsnclus_out, "bbpsnclus/I" );
  P->Branch( "bbshnblk", &bbshnblk_out, "bbshnblk/I" );
  P->Branch( "bbpsnblk", &bbpsnblk_out, "bbpsnblk/I" );
  P->Branch( "bbtrN", &bbtrN_out, "bbtrN/I" );
  P->Branch( "passedfid", &passedfid_out, "passedfid/I" );
  P->Branch( "run", &run_out, "run/I" );
  P->Branch( "mag", &mag_out, "mag/I" );

  //HCal edge boundaries
  Double_t left;
  Double_t right;
  Double_t top;
  Double_t bottom;

  //SBS-4 SBS-7 (pass0)
  if( pass==0 && (kine==4 || kine==7) ){
    left = econst::hcalposYi_p0;
    right = econst::hcalposYf_p0;
    top = econst::hcalposXf_p0;
    bottom = econst::hcalposXi_p0;
  }else{
    //pass1 and >pass1
    left = econst::hcalposYi_mc;
    right = econst::hcalposYf_mc;
    top = econst::hcalposXf_mc;
    bottom = econst::hcalposXi_mc;
  }

  //Active area boundaries
  Double_t leftAA;
  Double_t rightAA;
  Double_t topAA;
  Double_t bottomAA;

  //SBS-4 SBS-7 (pass0)
  if( pass==0 && (kine==4 || kine==7) ){
    leftAA = (econst::hcalposYi_p0+econst::hcalblk_w_p0);
    rightAA = (econst::hcalposYf_p0-econst::hcalblk_w_p0);
    topAA = (econst::hcalposXf_p0-econst::hcalblk_h_p0);
    bottomAA = (econst::hcalposXi_p0+econst::hcalblk_h_p0);
  }else{
    //pass1 and >pass1
    leftAA = (econst::hcalposYi_mc+econst::hcalblk_div_h);
    rightAA = (econst::hcalposYf_mc-econst::hcalblk_div_h);
    topAA = (econst::hcalposXf_mc-econst::hcalblk_div_v);
    bottomAA = (econst::hcalposXi_mc+econst::hcalblk_div_v);
  }
  
  //Safety Margin boundaries
  Double_t leftSM;
  Double_t rightSM;
  Double_t topSM;
  Double_t bottomSM;

  //Cut all events that are projected to outermost edge of HCal
  //Add 3sigma proton peak safety margin (x and y) to ensure no expected detections lie one boundary of HCal

  //SBS-4 SBS-7 (pass0)
  if( pass==0 && (kine==4 || kine==7) ){
    leftSM = (econst::hcalposYi_p0+econst::hcalblk_w_p0+3*dysig);
    rightSM = (econst::hcalposYf_p0-econst::hcalblk_w_p0-3*dysig);
    topSM = (econst::hcalposXf_p0-econst::hcalblk_h_p0-3*dxsig_p-dx0_p);
    bottomSM = (econst::hcalposXi_p0+econst::hcalblk_h_p0+3*dxsig_p-dx0_p);
  }else{
    //pass1 and >pass1
    leftSM = (econst::hcalposYi_mc+econst::hcalblk_div_h+3*dysig);
    rightSM = (econst::hcalposYf_mc-econst::hcalblk_div_h-3*dysig);
    topSM = (econst::hcalposXf_mc-econst::hcalblk_div_v-3*dxsig_p-dx0_p);
    bottomSM = (econst::hcalposXi_mc+econst::hcalblk_div_v+3*dxsig_p-dx0_p);
  }

  // setup reporting indices
  std::cout << std::endl << std::endl;

  for (Int_t irun=0; irun<nruns; irun++) {

    //DEBUG
    if(debug&&irun>5) continue;

    // accessing run info
    Int_t runnum = corun[irun].runnum;
    Int_t mag = corun[irun].sbsmag / 21; //convert to percent
    Double_t ebeam = corun[irun].ebeam; //get beam energy per run
    std::string tar = corun[irun].target;

    //Only proceed for deuterium at desired field setting
    if( mag!=magset || tar.compare("LD2")!=0 ) 
      continue;

    std::cout << "Analyzing run " << runnum << ".." << std::endl;

    std::string rfname;
    std::string rfname_sbs8;
    if( kine==8 ){
      TString pattern = Form("*%d*",corun[irun].runnum);
      TString foundDir = util::FindFileInDirectories(pattern, directories);
      if (!foundDir.IsNull()) {
	std::cout << "SBS 8 ld2 file found in: " << foundDir << std::endl;
	rfname_sbs8 = foundDir;
      } else {
	std::cout << "SBS 8 ld2 file not found in " << foundDir << std::endl;
      }
      rfname = rfname_sbs8 + Form("/*%d*",corun[irun].runnum);
    }else
      rfname = rootfile_dir + Form("/*%d*",corun[irun].runnum);

    //first attempt to resolve segmentation fault on large data sets
    if (C != nullptr) {
      delete C;
    }

    C = new TChain("T");
    C->Add(rfname.c_str());

    // setting up ROOT tree branch addresses
    C->SetBranchStatus("*",0);    

    // HCal general
    Double_t hcalid, hcale, hcalx, hcaly, hcalr, hcalc, hcaltdc, hcalatime, hcalidx, hcalnclus;
    std::vector<std::string> hcalvar = {"idblk","e","x","y","rowblk","colblk","tdctimeblk","atimeblk","index","nclus"};
    std::vector<void*> hcalvarlink = {&hcalid,&hcale,&hcalx,&hcaly,&hcalr,&hcalc,&hcaltdc,&hcalatime,&hcalidx,&hcalnclus};
    rvars::setbranch(C, "sbs.hcal", hcalvar, hcalvarlink);

    // HCal cluster branches (primary block information for id, tdc, and atime)
    Double_t hcalcid[econst::maxclus], hcalce[econst::maxclus], hcalcx[econst::maxclus], hcalcy[econst::maxclus], hcalctdctime[econst::maxclus], hcalcatime[econst::maxclus], hcalcnblk[econst::maxclus];
    //Double_t hcalcnblk;
    Int_t Nhcalcid;
    std::vector<std::string> hcalcvar = {"id","e","x","y","tdctime","atime","nblk","id"};
    std::vector<void*> hcalcvarlink = {&hcalcid,&hcalce,&hcalcx,&hcalcy,&hcalctdctime,&hcalcatime,&hcalcnblk,&Nhcalcid};
    rvars::setbranch(C, "sbs.hcal.clus", hcalcvar, hcalcvarlink, 7);

    // hodoscope cluster mean time
    Int_t Nhodotmean; 
    Double_t hodotmean[econst::maxclus];
    std::vector<std::string> hodovar = {"clus.tmean","clus.tmean"};
    std::vector<void*> hodovarlink = {&Nhodotmean,&hodotmean};
    rvars::setbranch(C, "bb.hodotdc", hodovar, hodovarlink, 0); 

    // BBCal shower timing
    Double_t atime_sh, e_sh, nclus_sh, nblk_sh;
    std::vector<std::string> shvar = {"atimeblk","e","nclus","nblk"};
    std::vector<void*> shvarlink = {&atime_sh,&e_sh,&nclus_sh,&nblk_sh};
    rvars::setbranch(C, "bb.sh", shvar, shvarlink);  

    // BBCal preshower timing
    Double_t atime_ps, e_ps, nclus_ps, nblk_ps;
    std::vector<std::string> psvar = {"atimeblk","e","nclus","nblk"};
    std::vector<void*> psvarlink = {&atime_ps,&e_ps,&nclus_ps,&nblk_ps};
    rvars::setbranch(C, "bb.ps", psvar, psvarlink); 

    // BB GEM hits
    Double_t gem_hits[econst::maxtrack];
    std::vector<std::string> gemvar = {"track.nhits"};
    std::vector<void*> gemvarlink = {&gem_hits};
    rvars::setbranch(C, "bb.gem", gemvar, gemvarlink);      

    // track branches
    Double_t ntrack, p[econst::maxtrack],px[econst::maxtrack],py[econst::maxtrack],pz[econst::maxtrack],xtr[econst::maxtrack],ytr[econst::maxtrack],thtr[econst::maxtrack],phtr[econst::maxtrack];
    Double_t vx[econst::maxtrack],vy[econst::maxtrack],vz[econst::maxtrack];
    Double_t xtgt[econst::maxtrack],ytgt[econst::maxtrack],thtgt[econst::maxtrack],phtgt[econst::maxtrack];
    std::vector<std::string> trvar = {"n","p","px","py","pz","x","y","th","ph","vx","vy","vz","tg_x","tg_y","tg_th","tg_ph"};
    std::vector<void*> trvarlink = {&ntrack,&p,&px,&py,&pz,&xtr,&ytr,&thtr,&phtr,&vx,&vy,&vz,&xtgt,&ytgt,&thtgt,&phtgt};
    rvars::setbranch(C,"bb.tr",trvar,trvarlink);
    
    // tdctrig branches
    Int_t Ntdctrigid;
    Double_t tdctrig[econst::maxtrack], tdctrigid[econst::maxtrack];
    std::vector<std::string> tdcvar = {"tdcelemID","tdcelemID","tdc"};
    std::vector<void*> tdcvarlink = {&tdctrigid,&Ntdctrigid,&tdctrig};
    rvars::setbranch(C,"bb.tdctrig",tdcvar,tdcvarlink,1);
  
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
    //C->SetBranchStatus("bb.gem.track.nhits", 1);
    C->SetBranchStatus("bb.etot_over_p", 1);
    C->SetBranchStatus("bb.tr.n", 1);
    //C->SetBranchStatus("bb.ps.e", 1);

    //set up the global cut formula
    TCut GCut = gcut.c_str();
    TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );

    //get each globalcut
    vector<string> allGlobalCuts = util::parseGlobalCut(gcut);

    // get experimental quantities by run
    std::cout << "Uncorrected average beam energy on " << tar << " for run: " << ebeam << std::endl;
    //set up hcal coordinate system with hcal angle wrt exit beamline
    vector<TVector3> hcalaxes; vars::sethcalaxes( hcaltheta, hcalaxes );
    //TVector3 hcalorigin = hcaldist*hcalaxes[2] + econst::hcalvoff*hcalaxes[0];
    TVector3 hcalorigin = hcaldist*hcalaxes[2] + hcal_v_offset*hcalaxes[0];
    Double_t BdL = econst::sbsmaxfield * econst::sbsdipolegap * (mag/100); //scale crudely by magnetic field
    Double_t Eloss_outgoing = econst::celldiameter/2.0/sin(bbthr) * econst::lh2tarrho * econst::lh2dEdx;

    // set nucleon for LD2
    std::string nucleon = "np";

    // event indices
    long nevent = 0, nevents = C->GetEntries(); 

    //ttree formula markers
    Int_t treenum = 0, currenttreenum = 0;

    while (C->GetEntry(nevent++)) {
      
      std::cout << "Processing run (" << irun << "/" << nruns << ") " << corun[irun].runnum << " event " << nevent << "/" << nevents << "\r";
      std::cout.flush();
      
      ///////
      //Single-loop globalcut method. Save pass/fail for output tree.
      currenttreenum = C->GetTreeNumber();
      if( nevent == 1 || currenttreenum != treenum ){
	treenum = currenttreenum; 
	GlobalCut->UpdateFormulaLeaves();
      }
      bool failedglobal = GlobalCut->EvalInstance(0) == 0;

      if(currenttreenum != treenum){
	cout << "Multiple tree numbers for a single run" << endl;
	return;
      }

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
      TLorentzVector ptarg; vars::setPN(nucleon,ptarg); //target momentum
      TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
      TVector3 qv = q.Vect();
      TLorentzVector pN; //N' momentum
      TVector3 pNhat; //Unit N' 3-vector
      
      //simple calculations for e' and N'
      Double_t etheta = vars::etheta(pe); 
      Double_t ephi = vars::ephi(pe);
      Double_t pcent = vars::pcentral(ebeam,etheta,nucleon); //e' p reconstructed by angles
      Double_t phNexp = ephi + physconst::pi;
      Double_t Q2, W2, nu, thNexp, pNexp;
      Double_t ebeam_o = vars::ebeam_o( ebeam_c, etheta, target.c_str() ); //Second energy correction accounting for energy loss leaving target

      //Calculate Mott cross section for this event
      Double_t MCS = (pow(physconst::alpha,2)*pow(cos(etheta/2),2)*precon)/(4*pow(ebeam_c,3)*pow(sin(etheta/2),4));

      /* Can reconstruct e' momentum for downstream calculations differently:
      v1 - Use four-momentum member functions
      v2 - Use all available ekine (tree) vars and calculate vectors (should be the same as v1)
      v3 - Use reconstructed angles as independent qty (usually preferable given GEM precision at most kinematics)
      v4 - Use reconstructed momentum as independent qty */
      
      if( epm==1 ){
	//v1
	pN = q + ptarg;
	pNhat = pN.Vect().Unit();
	Q2 = -q.M2();
	W2 = pN.M2();
	nu = q.E();
      }else if( epm==2 ){
	//v2
	Q2 = ekineQ2;
	W2 = ekineW2;
	nu = ekinenu;
	pNexp = vars::pN_expect( nu, nucleon );
	thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	pNhat = vars::pNhat_track( thNexp, phNexp );
	pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
      }else if( epm==3 ){
	//v3
	nu = pbeam.E() - pcent;
	//nu = pbeam.E() - precon;
	pNexp = vars::pN_expect( nu, nucleon );
	thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	pNhat = vars::pNhat_track( thNexp, phNexp );
	pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	Q2 = vars::Q2( pbeam.E(), pe.E(), etheta );
	W2 = vars::W2( pbeam.E(), pe.E(), Q2, nucleon );
      }else if( epm==4 ){
	//v4
	nu = pbeam.E() - pe.E();
	pNexp = vars::pN_expect( nu, nucleon );
	thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	pNhat = vars::pNhat_track( thNexp, phNexp );
	pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	Q2 = vars::Q2( pbeam.E(), pe.E(), etheta );
	W2 = vars::W2( pbeam.E(), pe.E(), Q2, nucleon );
      }else{
	Q2 = ekineQ2;
	W2 = ekineW2;
	nu = ekinenu;
	pNexp = vars::pN_expect( nu, nucleon );
	thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	pNhat = vars::pNhat_track( thNexp, phNexp );
	pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	std::cout << "Warning: epm version incorrect. Defaulting to version 2." << endl;
      }

      /////////////////////
      //W2 elastic cut bool
      bool failedW2 = W2<W2min || W2>W2max;

      ///////
      //HCal active area cut (acceptance matching). Save pass/fail for output tree.
      vector<Double_t> xyhcalexp; vars::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );

      //Fill for active area check
      hHcalXY->Fill(hcaly,hcalx);								    

      //////////////////////
      //ALL CLUSTER ANALYSIS
      
      //Set up clone clusters for selection analysis.
      vector<double> clone_cluster_intime;
      vector<double> clone_cluster_score;

      //loop through all clusters and select without HCal position information
      for( int c=0; c<Nhcalcid; c++ ){
	
	//calculate h-arm physics quantities per cluster
	double atime = hcalcatime[c];
	double atime_diff = atime - atime_sh; //Assuming best shower time on primary cluster
	double ce = hcalce[c];
	
	//using hcal atime until after pass2, wide cut around 5 sigma
	bool passedCoin = abs(atime_diff-coin_profile[1])<5*coin_profile[2];
	bool passedE = ce>minE;

	//Replicate the in-time algorithm with new cluster to be sorted later
	clone_cluster_intime.push_back(ce);
	if( !passedCoin )
	  clone_cluster_intime[c] = 0;
	
	//Get score (no position info). Will be sorted later
	double cascore = util::assignScore( ce, atime_diff, hcalce[(Int_t)hcalidx], coin_profile);
	clone_cluster_score.push_back(cascore);
	
      }//endloop over cluster elements

      Int_t cidx_e = (Int_t)hcalidx;
      
      if( hcale != hcalce[(Int_t)hcalidx] ){
	cerr << "ERROR: Sorting failure. Debug and rerun." << endl;
	cout << "Index: " << (Int_t)hcalidx << ", hcale:" << hcale << endl;
      }

      //Get best score/intime indices from clone clusters
      Int_t score_idx = -1;
      Int_t intime_idx = -1;
      Double_t score = 0.;
      Double_t intime = 0.;
      for( int c=0; c<Nhcalcid; c++ ){
	if( clone_cluster_score[c]>score ){
	  score = clone_cluster_score[c];
	  score_idx = c;
	}
	if( clone_cluster_intime[c]>intime ){
	  intime = clone_cluster_intime[c];
	  intime_idx = c;
	}
      }

      Int_t cidx_intime = intime_idx;
      Int_t cidx_score = score_idx;

      //Switch between best clusters for systematic analysis
      Int_t cidx_best;
      
      switch (cluster_idx) {
      case 1:
	cidx_best = 0;
	break;
      case 2:
	cidx_best = cidx_e;
	break;
      case 3:
	cidx_best = cidx_intime;
	break;
      case 4:
	cidx_best = cidx_score;
	break;
      default:
	cidx_best = 3;
      }
      
      //Calculations from the best cluster
      Double_t dx_bestcluster = hcalcx[cidx_best] - xyhcalexp[0];
      Double_t dy_bestcluster = hcalcy[cidx_best] - xyhcalexp[1];
      Double_t hatime_bestcluster = hcalcatime[cidx_best];
      Double_t hcoin_bestcluster = hcalcatime[cidx_best] - atime_sh;
      Double_t hcoin_pcluster = hcalcatime[0] - atime_sh;
      Double_t ce_bestcluster = hcalce[cidx_best];
      Int_t cnblk = (Int_t)hcalcnblk[cidx_best];
      //bool passedcoin = abs( hatime_bestcluster - atimediff0 ) < coin_sigma_factor*atimediffsig;

      //H-arm fiducial cuts
      bool hcalON = cut::hcalaaON(hcalcx[cidx_best],hcalcy[cidx_best],hcalaa);
      vector<Double_t> fid = cut::hcalfid(dxsig_p,dysig,hcalaa);
      bool passed_fid = cut::hcalfidIN(xyhcalexp[0],xyhcalexp[1],dx0_p,fid);

      //Fill analysis tree variables before making cuts on systematicInfo
      if(systematicInfo){
	dx_out = dx_bestcluster;
	dy_out = dy_bestcluster;
	xexp_out = xyhcalexp[0];
	yexp_out = xyhcalexp[1];
	hcalx_out = hcalcx[cidx_best];
	hcaly_out = hcalcy[cidx_best];
	W2_out = W2;
	Q2_out = Q2;
	mott_out = MCS;
	hcalE_out = ce_bestcluster;
	bbEtot_out = e_sh + e_ps;
	bbshE_out = e_sh;
	bbpsE_out = e_ps;
	hcalatime_out = hatime_bestcluster;
	bbshatime_out = atime_sh;
	bbpsatime_out = atime_ps;
	bbgemNhits_out = gem_hits[0];
	bbtrx_out = xtr[0];
	bbtry_out = ytr[0];
	bbtrp_out = p[0];
	coinP1_out = coin_profile[1];
	coinP2_out = coin_profile[2];
	dy0_out = dy0;
	dx0_p_out = dx0_p;
	dx0_n_out = dx0_n;
	dysig_out = dysig;
	dxsig_p_out = dxsig_p;
	dxsig_n_out = dxsig_n;
	W2mean_out = W2mean;
	W2sig_out = W2sig;

	hcalnclus_out = Nhcalcid;
	hcalnblk_out = cnblk;
	bbshnclus_out = (Int_t)nclus_sh;
	bbshnblk_out = (Int_t)nblk_sh;
	bbpsnclus_out = (Int_t)nclus_ps;
	bbpsnblk_out = (Int_t)nblk_ps;
	bbtrN_out = (Int_t)ntrack;
	passedfid_out = (Int_t)passed_fid;
	run_out = runnum;
	mag_out = mag;
      }      

      if(!failedglobal)
	hdx_cut_globalonly->Fill(dx_bestcluster);

      if(!failedW2)
	hdx_cut_W2only->Fill(dx_bestcluster);

      //E-arm only cuts first
      if( failedglobal || failedW2 )
	continue;

      //Fill mott cross section after elastic cut
      hMott_cs->Fill(MCS);

      //quick check on hcal geometry
      hxy_nocut->Fill(hcalcy[cidx_best],hcalcx[cidx_best]);
      hxy_nocut->Fill(0.,0.);
      if(hcalON)
	hxy_aacut->Fill(hcalcy[cidx_best],hcalcx[cidx_best]);
      if(passed_fid)
	hxy_acccut->Fill(hcalcy[cidx_best],hcalcx[cidx_best]); //this is largely useless

      //Fill e-arm only cut histograms
      hcoin->Fill(hcoin_bestcluster);
      hcoin_pclus->Fill(hcoin_pcluster);
      hW2_nocut->Fill(W2);
      hdxdy_nocut->Fill(dy_bestcluster,dx_bestcluster);
      hdx_nocut->Fill(dx_bestcluster);
      hxy_exp_n->Fill(xyhcalexp[1],xyhcalexp[0]);
      hxy_exp_p->Fill(xyhcalexp[1],xyhcalexp[0]+dx0_p); //expected proton position from average difference between dxdy neutron and proton spots

      if(passed_fid){
	hxy_exp_n_fid->Fill(xyhcalexp[1],xyhcalexp[0]);
	hxy_exp_p_fid->Fill(xyhcalexp[1],xyhcalexp[0]+dx0_p); //expected proton position from average difference between dxdy neutron and proton spots
      }

      //Both arm cuts
      bool failedcoin = abs( hcoin_bestcluster - coin_profile[1] ) > coin_sigma_factor*coin_profile[2];
      bool faileddy = abs( dy_bestcluster - dy0 ) > 3*dysig;

      if( !faileddy && !failedcoin && hcalON ){
	hdx_cut_nofid->Fill(dx_bestcluster); //primary dx histo, no fiducial cut
	hdxdy_cut_nofid->Fill(dy_bestcluster,dx_bestcluster);
	if(!passed_fid)
	  hdx_cut_failfid->Fill(dx_bestcluster); //primary dx histo, 
      }
      
      //Make both arm cuts
      if( !passed_fid || !hcalON || failedcoin || faileddy )
	continue;

      hdxvE->Fill(ce_bestcluster,dx_bestcluster);
      hcoin_cut->Fill(hcoin_bestcluster);
      hcoin_pclus_cut->Fill(hcoin_pcluster);
      hW2_cut->Fill(W2);
      hdxdy_cut->Fill(dy_bestcluster,dx_bestcluster);
      hdx_cut->Fill(dx_bestcluster); //primary dx histo, fiducial cut
      hHcalXY_allcuts->Fill(hcalcy[cidx_best],hcalcx[cidx_best]);

      P->Fill();

    } //end event loop
    
    // reset chain for the next run config
    C->Reset();
    
  } //end run loop
  
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
