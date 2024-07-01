//sseeds 03.31.23 - test script to use analysis framework to produce hcal detection efficiency output more efficiently
//04.10.23 Update - Added W2 interpolate method for obtaining background fits to both total W2 distribution and W2 with HCal anticut 
//04.11.23 Update - Added back direct dx "detected" yield method for comparison. Fixed thetapq calculation and included in elastic cuts.
//5.20.23 Update - broke from general script and focused this script on extraction of hcal detection efficiency directly from dx after dy cuts normalized by strong earm elastic cuts. 
//5.30.23 Update - NOTE that fits (and resulting yields) are very sensitive to fit ranges. Background shape varies considerably and fit range should reflect 3sigma about the dx elastic peak to avoid overcounting

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

const Int_t atimeSigFac = 5;
const Double_t dx_fitlow = -0.3;
const Double_t dx_fithigh = 0.3;

//Total fits using Interpolate with elastic signal histo and 6th order poly fit to bg
TH1D *hW2elas;
Double_t W2total(Double_t *x, Double_t *par){
  Double_t W2 = x[0];
  Double_t sig_scale = par[0];
  Double_t signal = sig_scale * hW2elas->Interpolate(W2);
  return signal + fits::g_p6fit(x,&par[1]);
}

TH1D *hW2elasres;
Double_t W2totalres(Double_t *x, Double_t *par){
  Double_t W2 = x[0];
  Double_t sig_scale = par[0];
  Double_t signal = sig_scale * hW2elasres->Interpolate(W2);
  return signal + fits::g_p6fit(x,&par[1]);
}

TH1D *hdxelastic;
Double_t dxtotal(Double_t *x, Double_t *par){ //get poly fit to bg with scaled fit to "pure elastics"
  Double_t dx= x[0];
  Double_t sig_scale = par[0];
  Double_t signal = sig_scale * hdxelastic->Interpolate(dx);
  return signal + fits::g_p6fit(x,&par[1]);
}

//Pass only kinematic and sbs magnetic field setting. Config file for all kinematics located ../../config/shde.json
//cluster_idx refers to cluster selection method:
//// 1) Highest Energy cluster
//// 2) Cluster which minimizes the time between expected adc time and primary block adc time
//// 3) Cluster which minimizes proton thetapq
//// 4) Cluster which minimizes distance between cluster center and expected location of scattered nucleon
//// 5) Cluster which passes 5 sigma cut on adc time, then minimizes proton theta pq.
//// 6) Cluster which passes 5 sigma cut on adc time, then maximizes energy.
//// 7) Cluster which maximizes score (see algorithm in util.C)
void hde_dxdirect( Int_t kine=4, Int_t magset=30, Int_t cluster_idx=7)
{ //main  

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  //Some comparison specific variables
  const Double_t pcsigfac = 0.58; //sbs4
  const Double_t sampfrac = 0.0641; //sbs4
  const Double_t dxhmin = -10; //hcal, GeV
  const Double_t dxhmax = 10; //hcal, GeV
  const Double_t dyhmin = -10; //hcal, GeV
  const Double_t dyhmax = 10; //hcal, GeV
  const Int_t nfac = 100; //Number of bins for hcal E vs nucleon p fits obtained to ensure 1000 events per bin

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/shde.json");

  std::string rootfile_dir = jmgr->GetValueFromSubKey_str( "rootfile_dir", Form("sbs%d",kine) );
  Double_t W2fitmax = jmgr->GetValueFromSubKey<Double_t>( "W2fitmax", Form("sbs%d",kine) );
  Double_t binfac = jmgr->GetValueFromSubKey<Double_t>( "binfac", Form("sbs%d",kine) );
  Double_t hbinfac = jmgr->GetValueFromSubKey<Double_t>( "hbinfac", Form("sbs%d",kine) );
  Double_t thetapqcut = jmgr->GetValueFromSubKey<Double_t>( "thetapqcut", Form("sbs%d",kine) );
  std::string target = jmgr->GetValueFromKey_str( "target" ); //Must always be LH2 data for proton hde
  Int_t epm = jmgr->GetValueFromSubKey<Int_t>( "epm", Form("sbs%d",kine) );
  Int_t pass = jmgr->GetValueFromSubKey<Int_t>( "pass", Form("sbs%d",kine) );
  Double_t ebound_l = jmgr->GetValueFromSubKey<Double_t>( "ebound_l", Form("sbs%d",kine) );
  Double_t ebound_h = jmgr->GetValueFromSubKey<Double_t>( "ebound_h", Form("sbs%d",kine) );
  vector<Double_t> coin_profile;
  jmgr->GetVectorFromSubKey<Double_t>("coin_profile",Form("sbs%d",kine),coin_profile);
  Double_t hcal_v_offset = jmgr->GetValueFromSubKey<Double_t>( "hcal_offset", Form("sbs%d",kine) );

  if( pass>1 ){
    std::cout << "As of 3.31.23, the highest GMn replay pass is 1. Enter a valid pass." << endl;
    return;
  }

  //set up default parameters for all analysis
  std::string runsheet_dir = "/w/halla-scshelf2102/sbs/seeds/ana/data"; //unique to my environment for now
  Int_t nruns = -1; //Always analyze all available runs
  Int_t verb = 0; //Don't print diagnostic info by default
  Double_t hcalfit_l = econst::hcalposXi_p0-2*econst::hcalblk_w_p0; //lower fit/bin limit for hcal dx plots (m)
  Double_t hcalfit_h = econst::hcalposXf_p0+2*econst::hcalblk_h_p0; //upper fit/bin limit for hcal dx plots (m)
  Double_t harmrange = (hcalfit_h) - (hcalfit_l); //Full range of hcal dx plots (m)

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<crun> corun; 
  util::ReadRunList(runsheet_dir,nruns,kine,target.c_str(),pass,verb,corun); //modifies nruns to be very large when -1. Select run list for LH2 only for hde

  //set up output files
  string outfilename = Form("outfiles/hde_p1_sbs%d_%s_epm%d_magset%d.root",kine,target.c_str(),epm,magset);
  TFile *fout = new TFile( outfilename.c_str(), "RECREATE" );

  ////////////
  //HISTOGRAMS

  //Basic
  TH2D *hxy_nocut = new TH2D("hxy_nocut","HCal x vs HCal y, no cut; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_aacut = new TH2D("hxy_aacut","HCal x vs HCal y, acceptance matching cut no sigma; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxy_acccut = new TH2D("hxy_acccut","HCal x vs HCal y, acceptance matching cut; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );

  TH2D *hdxdy_nocut = new TH2D("hdxdy_nocut","HCal dxdy, no cut; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hdxdy_spotcut = new TH2D("hdxdy_spotcut","HCal dxdy, cut on proton spot; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hdxdy_earm_cut = new TH2D("hdxdy_earm_cut","HCal dxdy, e-arm elastic cut; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH1D *hW2_earm_cut = new TH1D( "hW2_earm_cut", "W^{2}, global/dy/dx/thetapq cut; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  

  //set up diagnostic histograms
  TH2D *hxyexp_nocut = new TH2D("hxyexp_nocut","HCal X vs Y Expected from BB, No Cut;dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxyexp_aacut = new TH2D("hxyexp_aacut","HCal X vs Y Expected from BB, Acceptance Match Cut;dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );

  //set up detection efficiency histograms
  TH1D *hW2_allcut = new TH1D( "hW2_allcut", "W^{2}, global/dy/dx/thetapq cut;W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_anticut = new TH1D( "hW2_anticut", "W^{2}, elastic anticut (global e-arm and hcal dy/dx);W^{2} (GeV^{2});", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_nocut = new TH1D( "hW2_nocut", "W^{2}, no cut;W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hdx_allcut = new TH1D( "hdx_allcut", "dx, W^{2} and dy cut;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);
  TH1D *hdx_nocut = new TH1D( "hdx_nocut", "dx, no cut;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);
  TH1D *hdx_anticut = new TH1D( "hdx_anticut", "dx, e-arm elastic anticut;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);

  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;

  // create output tree
  TTree *P = new TTree("P","Analysis Tree"); 

  // output tree vars
  Double_t W2_out;
  Double_t Q2_out;
  Double_t nu_out;
  Double_t ep_out; //track reconstructed e' momentum
  Double_t hcalpexp_out;
  Int_t mag_out;
  Int_t run_out; //run number
  Int_t failedglobal_out; //earm global variable cuts only
  Int_t failedaccmatch_out;
  Int_t cidx_best_atime_out;
  Int_t cidx_best_tpq_out;
  Int_t cidx_best_dxdy_out;
  Int_t cidx_best_tpqatime_out;
  Double_t bc_atime_out;
  Double_t bc_hcale_out;
  Double_t bc_thpq_out;
  Double_t bc_dx_out;
  Double_t bc_dy_out;

  //cluster tree vars
  Double_t cpblkid_out[econst::maxclus];
  Double_t chatime_out[econst::maxclus];
  Double_t chcale_out[econst::maxclus];
  Double_t cthetapq_p_out[econst::maxclus];
  Double_t cthetapq_n_out[econst::maxclus];
  Double_t cdx_out[econst::maxclus];
  Double_t cdy_out[econst::maxclus];
  Int_t cpid_out[econst::maxclus]; //-1:neither,0:ambiguous,1:proton,2:neutron
  Int_t cfailedhcaltime_out[econst::maxclus];
  Int_t cinspot_p_out[econst::maxclus];

  // set output tree branches
  P->Branch( "W2", &W2_out, "W2/D" );
  P->Branch( "Q2", &Q2_out, "Q2/D" );
  P->Branch( "nu", &nu_out, "nu/D" );
  P->Branch( "ep", &ep_out, "ep/D" );
  P->Branch( "hcalpexp", &hcalpexp_out, "hcalpexp/D" );
  P->Branch( "mag", &mag_out, "mag/I" );
  P->Branch( "run", &run_out, "run_out/I" );
  P->Branch( "failedglobal", &failedglobal_out, "failedglobal/I" );
  P->Branch( "failedaccmatch", &failedaccmatch_out, "failedaccmatch/I" );
  P->Branch( "cidx_best_atime", &cidx_best_atime_out, "cidx_best_atime/I" );
  P->Branch( "cidx_best_tpq", &cidx_best_tpq_out, "cidx_best_tpq/I" );
  P->Branch( "cidx_best_dxdy", &cidx_best_dxdy_out, "cidx_best_dxdy/I" );
  P->Branch( "cidx_best_tpqatime", &cidx_best_tpqatime_out, "cidx_best_tpqatime/I" );
  P->Branch( "bc_atime", &bc_atime_out, "bc_atime/D" );
  P->Branch( "bc_hcale", &bc_hcale_out, "bc_hcale/D" );
  P->Branch( "bc_thpq", &bc_thpq_out, "bc_thpq/D" );
  P->Branch( "bc_dx", &bc_dx_out, "bc_dx/D" );
  P->Branch( "bc_dy", &bc_dy_out, "bc_dy/D" );

  P->Branch( "cpblkid", &cpblkid_out, Form("cpblkid[%d]/D",econst::maxclus) );
  P->Branch( "chatime", &chatime_out, Form("chatime[%d]/D",econst::maxclus) );
  P->Branch( "chcale", &chcale_out, Form("chcale[%d]/D",econst::maxclus) );
  P->Branch( "cthetapq_p", &cthetapq_p_out, Form("cthetapq_p[%d]/D",econst::maxclus) );
  P->Branch( "cthetapq_n", &cthetapq_n_out, Form("cthetapq_n[%d]/D",econst::maxclus) );
  P->Branch( "cdx", &cdx_out, Form("cdx[%d]/D",econst::maxclus) );
  P->Branch( "cdy", &cdy_out, Form("cdy[%d]/D",econst::maxclus) );
  P->Branch( "cpid", &cpid_out, Form("cpid[%d]/I",econst::maxclus) );
  P->Branch( "cfailedhcaltime", &cfailedhcaltime_out, Form("cfailedhcaltime[%d]/I",econst::maxclus) );
  P->Branch( "cinspot_p", &cinspot_p_out, Form("cinspot_p[%d]/I",econst::maxclus) );

  // setup reporting indices
  Int_t curmag = -1;
  std::string curtar = "";
  std::cout << std::endl << std::endl;
  Double_t dxmean;
  Double_t dxsigma;

  Double_t dy0;
  Double_t dysig;
  Double_t dx0_p;
  Double_t dxsig_p;

  for (Int_t irun=0; irun<nruns; irun++) {
    // accessing run info
    Int_t runnum = corun[irun].runnum;
    Int_t mag = corun[irun].sbsmag / 21; //convert to percent
    Double_t ebeam = corun[irun].ebeam; //get beam energy per run
    std::string tar = corun[irun].target;

    //Do not continue if current run has undesired sbs magnet config
    if( magset!=-1 && mag!=magset ) continue;
    std::cout << "Analyzing run " << runnum << ".." << std::endl;

    //set up configuration and tune objects to load analysis parameters
    //SBSconfig *config = new SBSconfig(kine,mag);
    SBSconfig config(kine,mag);

    //Obtain configuration pars from config file
    Double_t hcaltheta = config.GetHCALtheta_rad();
    Double_t hcaldist = config.GetHCALdist();
    Double_t sbsdist = config.GetSBSdist();
    Double_t bbthr = config.GetBBtheta_rad(); //in radians
    
    //SBStune *tune = new SBStune(kine,mag);
    SBStune tune(kine,mag);
    
    //Reporting. tar should always equal curtar as categorized by good run list
    if( tar.compare(curtar)!=0 || mag!=curmag ){
      std::cout << "Settings change.." << endl;
      cout << config;
      cout << tune;
      curtar = tar;
      curmag = mag;
    }

    //Obtain cuts from tune class
    std:string gcut   = tune.Getglobcut_earm();
    Double_t W2mean   = tune.GetW2mean();
    Double_t W2sig    = tune.GetW2sig();
    Double_t dx0_n    = tune.Getdx0_n();
    dx0_p    = tune.Getdx0_p();
    dy0      = tune.Getdy0();
    Double_t dxsig_n  = tune.Getdxsig_n();
    dxsig_p  = tune.Getdxsig_p();
    dysig    = tune.Getdysig();
    Double_t atime0   = tune.Getatime0();
    Double_t atimesig = tune.Getatimesig();

    //For this analysis, only one magnetic field
    dxmean = dx0_p;
    dxsigma = dxsig_p;

    //Set cut/anticut limits for proton
    Double_t dxmin = dx0_p - 3*dxsig_p; //Obtain 99.73% of elastic sample
    Double_t dxmax = dx0_p + 3*dxsig_p; //Obtain 99.73% of elastic sample
    Double_t dymin = dy0 - 2*dysig;
    Double_t dymax = dy0 + 2*dysig;
    Double_t W2min = W2mean - 2*W2sig;
    Double_t W2max = W2mean + 2*W2sig;
    Double_t atmin = atime0 - 6*atimesig;
    Double_t atmax = atime0 + 6*atimesig;

    //Set minimum energy for expected recoil nucleon to be considered
    Double_t hminE = 0.0;  //0.05 GeV or nothing

    //Acceptance matching limits for analysis
    // Double_t topAcc = (econst::hcalposYi_p0+econst::hcalblk_w_p0-dysig);
    // Double_t bottomAcc = (econst::hcalposYf_p0-econst::hcalblk_w_p0+dysig);
    // Double_t leftAcc = (econst::hcalposXf_p0-econst::hcalblk_h_p0+dxsig_p-dx0_p);
    // Double_t rightAcc = (econst::hcalposXi_p0+econst::hcalblk_h_p0-dxsig_p-dx0_p);

    Double_t leftAcc = (econst::hcalposYi_p0+econst::hcalblk_w_p0-dysig);
    Double_t rightAcc = (econst::hcalposYf_p0-econst::hcalblk_w_p0+dysig);
    Double_t topAcc = (econst::hcalposXf_p0-econst::hcalblk_h_p0+dxsig_p-dx0_p);
    Double_t bottomAcc = (econst::hcalposXi_p0+econst::hcalblk_h_p0-dxsig_p-dx0_p);

    std::string rfname = rootfile_dir + Form("/*%d*",corun[irun].runnum);
    //std::cout << "Switching to run " << rfname << ".." << std::endl;
    C = new TChain("T");
    C->Add(rfname.c_str());

    // setting up ROOT tree branch addresses
    C->SetBranchStatus("*",0);    

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
    
    // hodoscope cluster mean time
    Int_t Nhodotmean; 
    Double_t hodotmean[econst::maxclus];
    std::vector<std::string> hodovar = {"clus.tmean","clus.tmean"};
    std::vector<void*> hodovarlink = {&Nhodotmean,&hodotmean};
    rvars::setbranch(C, "bb.hodotdc", hodovar, hodovarlink, 0);  
    
    // BBCal shower timing
    Double_t atime_sh, nclus_sh;
    std::vector<std::string> shvar = {"atimeblk","nclus"};
    std::vector<void*> shvarlink = {&atime_sh,&nclus_sh};
    rvars::setbranch(C, "bb.sh", shvar, shvarlink);  

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
    TVector3 hcalorigin = hcaldist*hcalaxes[2] + hcal_v_offset*hcalaxes[0];
    Double_t BdL = econst::sbsmaxfield * econst::sbsdipolegap * (mag/100); //scale crudely by magnetic field
    Double_t Eloss_outgoing = econst::celldiameter/2.0/sin(bbthr) * econst::lh2tarrho * econst::lh2dEdx;
    //vector<Double_t> hcalaa_p0 = cut::hcalaa_data(1,1); //exclude one block on periphery

    // set nucleon defaults by target
    std::string nucleon;
    if( tar.compare("LH2")==0 )
      nucleon = "p"; 
    else if( tar.compare("LD2")==0 )      
      nucleon = "np";
    else{
      nucleon = "np";
      std::cout << "Error: incorrect target information loaded from good run list. Check and verify." << std::endl;
    }

    // event indices
    long nevent = 0, nevents = C->GetEntries(); 
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
	//cout << "Updating formula leaves and switching segment at event: " << nevent << endl;
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

      //Fill earm physics analysis tree variables
      W2_out = W2;
      Q2_out = Q2;
      nu_out = nu;
      ep_out = precon;
      hcalpexp_out = pNexp;

      /////////////////////
      //W2 elastic cut bool
      bool failedW2 = W2<W2min || W2>W2max;

      ///////
      //HCal active area cut (acceptance matching). Save pass/fail for output tree. Good only for zero field
      vector<Double_t> xyhcalexp; vars::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
      bool failedaccmatch = 
	xyhcalexp[1] > bottomAcc ||
	xyhcalexp[1] < topAcc ||
	xyhcalexp[0] > leftAcc ||
	xyhcalexp[0] < rightAcc;

      hxyexp_nocut->Fill( xyhcalexp[1], xyhcalexp[0] );
      if( !failedaccmatch )
	hxyexp_aacut->Fill( xyhcalexp[1], xyhcalexp[0] );
     
      //Get expected proton deflection for thetapq (zero field, so this shouldn't matter)
      Double_t protdeflect = tan( 0.3 * BdL / qv.Mag() ) * (hcaldist - (sbsdist + econst::sbsdipolegap/2.0) );

      //////////////////////
      //ALL CLUSTER ANALYSIS

      //Set up indexes and variables for best cluster analysis
      Int_t cidx_atime = 0;
      Double_t c_atimediff = 1000.;
      Int_t cidx_thetapq = 0;
      Double_t c_thetapqdiff = 1000.;
      Int_t cidx_dxdy = 0;
      Double_t c_dxdydiff = 1000.;
      Int_t cidx_tpq_atime = 0;
      Double_t c_tpq_atimediff = 1000.;
      Int_t cidx_e_coin = 0;
      Double_t c_e_coindiff = 0.;
      Int_t cidx_score = 0;
      Double_t c_score = 0.;

      Double_t tpq_p[econst::maxclus];

      Double_t max_e_among_clus = hcale;

      //////////////////////////
      //CLUSTERING AND SELECTION
      //loop through all clusters Ndata.sbs.hcal.clus.id
      for( Int_t c=0; c<Nhcalcid; c++ ){

	//calculate h-arm physics quantities per cluster
	Double_t cid = hcalcid[c]; //sbs.hcal.clus.id
	Double_t cdx = hcalcx[c] - xyhcalexp[0];  //sbs.hcal.clus.x
	Double_t cdy = hcalcy[c] - xyhcalexp[1];  //sbs.hcal.clus.y
	Double_t ce = hcalce[c];  //sbs.hcal.clus.e
	TVector3 chcalpos = hcalorigin + hcalcx[c]*hcalaxes[0] + hcalcy[c]*hcalaxes[1];  
	TVector3 cneutdir = ( chcalpos - vertex ).Unit();
	TVector3 cprotdir = ( chcalpos + protdeflect * hcalaxes[0] - vertex ).Unit();
	Double_t cthetapq_p = acos( cprotdir.Dot( qv.Unit() ) );
	tpq_p[c]=cthetapq_p;
	Double_t cthetapq_n = acos( cneutdir.Dot( qv.Unit() ) );
	Int_t cpid;
	bool cis_p = abs(cdx - dx0_p)<3*dxsig_p && abs(cdy - dy0)<3*dysig;
	bool cis_n = abs(cdx - dx0_n)<3*dxsig_n && abs(cdy - dy0)<3*dysig;
	Double_t c_cointime = hcalcatime[c] - atime_sh;  //sbs.hcal.clus.atime , bb.sh.atimeblk
	bool c_failedhcaltime = abs( hcalcatime[c] - atime0 ) > atimeSigFac*atimesig;
	bool c_failedcoin = abs( c_cointime - coin_profile[1] ) > atimeSigFac*coin_profile[2];
	bool c_spotcheck_p = util::Nspotcheck(cdy,cdx,dy0,dx0_p,dysig,dxsig_p,0);

	//determine pid from dx and dy information
	if( cis_p && !cis_n )
	  cpid=1;
	else if( cis_n && !cis_p )
	  cpid=2;
	else if( cis_p && cis_n )
	  cpid=0;
	else
	  cpid=-1;

	//////////////////////////////////////////////
	//get best cluster indexes for various methods

	//get cluster which minimizes coin diff (METHOD 2)
	if( abs(hcalcatime[c]-atime0)<c_atimediff ){
	  cidx_atime = c;
	  c_atimediff = hcalcatime[c];
	}
	//get cluster which minimizes thetapq (proton for hydrogen analysis) (METHOD 3)
	if( cthetapq_p<c_thetapqdiff ){
	  cidx_thetapq = c;
	  c_thetapqdiff = cthetapq_p;
	}
	//get cluster which minimizes distance between cluster center and expected loc dxdy (METHOD 4)
	Double_t dxdydist = sqrt( pow(cdx-dx0_p,2)+pow(cdy-dy0,2));
	if( dxdydist<c_dxdydiff ){
	  cidx_dxdy = c;
	  c_dxdydiff = dxdydist;
	}
	//get cluster which passes wide cut on atime, then minimizes thetapq (METHOD 5)
	if( !c_failedhcaltime && cthetapq_p<c_tpq_atimediff ){
	  cidx_tpq_atime = c;
	  c_tpq_atimediff = cthetapq_p; 
	}

	//get cluster which passes wide cut on coin, then maximizes energy (METHOD 6)
	if( !c_failedcoin && ce>c_e_coindiff ){
	  cidx_e_coin = c;
	  c_e_coindiff = ce;
	}
	//get cluster with best time/energy score (METHOD 7)
	Double_t score = util::assignScore(ce,c_cointime,max_e_among_clus,coin_profile);
	if(score>c_score){
	  cidx_score = c;
	  c_score = score;
	}

	//fill analysis tree variables
	chcale_out[c] = hcalce[c];
	cpblkid_out[c] = hcalcid[c];
	chatime_out[c] = hcalcatime[c];
	cthetapq_p_out[c] = cthetapq_p;
	cthetapq_n_out[c] = cthetapq_n;
	cdx_out[c] = cdx;
	cdy_out[c] = cdy;
	cpid_out[c] = cpid;
	cfailedhcaltime_out[c] = (Int_t)c_failedhcaltime;
	cinspot_p_out[c] = (Int_t)c_spotcheck_p;
       
      }
	
      //Switch between best clusters for systematic analysis
      Int_t cidx_best;

      switch (cluster_idx) {
      case 1:
	cidx_best = 0;
	break;
      case 2:
	cidx_best = cidx_atime;
	break;
      case 3:
	cidx_best = cidx_thetapq;
	break;
      case 4:
	cidx_best = cidx_dxdy;
	break;
      case 5:
	cidx_best = cidx_tpq_atime;
	break;
      case 6:
	cidx_best = cidx_e_coin;
	break;
      case 7:
	cidx_best = cidx_score;
	break;
      default:
	cidx_best = 0;
      }
      
      //Calculations from the best cluster for later use
      Double_t dx_bestcluster = hcalcx[cidx_best] - xyhcalexp[0];
      Double_t dy_bestcluster = hcalcy[cidx_best] - xyhcalexp[1];
      Double_t hatime_bestcluster = hcalcatime[cidx_best];
      Double_t thpq_bestcluster = tpq_p[cidx_best];
      Double_t ce_bestcluster = hcalce[cidx_best];

      //check if best cluster position is within the proton spot
      bool spotcheck_p = util::Nspotcheck(dy_bestcluster,dx_bestcluster,dy0,3*dx0_p,dysig,3*dxsig_p,0);
      //check if best cluster passes wide coin cut
      bool passedcoin = abs( hatime_bestcluster - atime0 ) < atimeSigFac*atimesig;

      //quick check on acceptance
      hxy_nocut->Fill(hcalcy[cidx_best],hcalcx[cidx_best]);
      hxy_nocut->Fill(0.,0.);

      bool aaregion = 
	hcalcy[cidx_best] > (econst::hcalposYf_p0-econst::hcalblk_w_p0) ||
	hcalcy[cidx_best] < (econst::hcalposYi_p0+econst::hcalblk_w_p0) ||
	hcalcx[cidx_best] > (econst::hcalposXf_p0-econst::hcalblk_h_p0) ||
	hcalcx[cidx_best] < (econst::hcalposXi_p0+econst::hcalblk_h_p0);

      if(aaregion)
	hxy_aacut->Fill(hcalcy[cidx_best],hcalcx[cidx_best]);

      bool accregion = 
      	hcalcy[cidx_best] > bottomAcc ||
	hcalcy[cidx_best] < topAcc ||
	hcalcx[cidx_best] > leftAcc ||
	hcalcx[cidx_best] < rightAcc;

      if(accregion)
	hxy_acccut->Fill(hcalcy[cidx_best],hcalcx[cidx_best]);

      hdxdy_nocut->Fill(dy_bestcluster,dx_bestcluster);
      if( spotcheck_p )
	hdxdy_spotcut->Fill(dy_bestcluster,dx_bestcluster);

      //Fill final output tree variables    
      mag_out = mag;
      run_out = runnum;
      failedglobal_out = (Int_t)failedglobal;
      failedaccmatch_out = (Int_t)failedaccmatch;
      cidx_best_atime_out = cidx_atime;
      cidx_best_tpq_out = cidx_thetapq;
      cidx_best_dxdy_out = cidx_dxdy;
      cidx_best_tpqatime_out = cidx_tpq_atime;
      bc_atime_out = hatime_bestcluster;
      bc_hcale_out = ce_bestcluster;
      bc_thpq_out = thpq_bestcluster;
      bc_dx_out = dx_bestcluster;
      bc_dy_out = dy_bestcluster;

      P->Fill();

      ////////////////////////////////////////////////
      //HCal detection efficiency - dx direct approach
      ////////////////////////////////////////////////
      if( !failedaccmatch ){ //cut on acceptance matching for all analysis
	
	//Obtain distributions with only acceptance matching
	hW2_nocut->Fill( W2 );
	hdx_nocut->Fill( dx_bestcluster );
	    
	//Make tightest elastic cut possible for elastic SHAPE only (W2 and dx)
	if( spotcheck_p &&
	    passedcoin &&
	    !failedglobal ){

	  if( thpq_bestcluster<0.04 ){
	    hW2_allcut->Fill( W2 ); //Don't cut on W2 for W2 distribution
	    //if( !failedW2 && ce_bestcluster>0.01 && ce_bestcluster<0.14 )
	    if( !failedW2 && ce_bestcluster>0.035 && ce_bestcluster<0.105)
	      hdx_allcut->Fill(dx_bestcluster); //Cut on dx in spotcheck may be too tight
	  }

	}

	//Make cut on e-arm only for dx distribution
	if( !failedglobal &&
	    !failedW2 )
	  hdxdy_earm_cut->Fill(dy_bestcluster,dx_bestcluster);
	  
	//Make hcal only elastic anticut
	if( thpq_bestcluster>0.04 &&
	    !spotcheck_p &&
	    !passedcoin )
	  hW2_anticut->Fill( W2 );
	  
	//Make e-arm only elastic anticut
	if( failedglobal &&
	    failedW2 )
	  hdx_anticut->Fill( dx_bestcluster );

      }
	 
    } //end event loop

    // reset chain for the next run config
    C->Reset();

  } //end run loop

  ///////////////////
  //W2 ANTICUT METHOD

  TH1D *hW2 = (TH1D*)(hW2_nocut->Clone("hW2"));
  TH1D *hW2res = (TH1D*)(hW2_anticut->Clone("hW2res"));
  TH1D *hdx = (TH1D*)(hdx_nocut->Clone("hdx")); 
  TH1D *hdx_2 = (TH1D*)(hdx_nocut->Clone("hdx")); 
  TH1D *hdx_3 = (TH1D*)(hdx_nocut->Clone("hdx")); 

  Double_t W2fullint = hW2->Integral(0,W2fitmax);
  Double_t W2resfullint = hW2res->Integral(0,W2fitmax);

  //Make canvas for (expected-residuals)/expected
  TCanvas *c1 = new TCanvas("c1",Form("HDE sbs%d mag%d",kine,magset),1200,500);
  gStyle->SetPalette(53);
  c1->Divide(2,1);
  c1->cd(1);

  //Acceptance matching cut only first with "pure elastic" hW2_allcut
  hW2elas = (TH1D*)(hW2_allcut->Clone("hW2elas"));
  TF1 *tfit = new TF1("tfit",W2total,0,W2fitmax,8); //tfit npar = 1+pNfit_npar+1
  TF1 *bg = new TF1("bg",fits::g_p6fit,0.,W2fitmax,7);
  tfit->SetLineColor(kGreen);
  hW2->SetTitle("W^{2}, Acceptance Cut Only");
  hW2->Fit("tfit","RBM");

  Double_t *tpar = tfit->GetParameters();

  //Get fit parameters for bg function and draw identical function on canvas
  bg->SetParameters(&tpar[1]);
  bg->SetLineColor(kRed);
  bg->SetFillColor(kRed);
  bg->SetFillStyle(3005);
  bg->Draw("same");
  hW2elas->SetLineColor(kBlue);
  hW2elas->SetLineWidth(2);
  hW2elas->SetFillColor(kBlue);
  hW2elas->SetFillStyle(3003);
  hW2elas->Draw("same");

  //get error params
  Int_t W2elasb = ebound_l*binfac;
  Int_t W2elase = ebound_h*binfac;
  TFitResultPtr s = hW2->Fit("tfit","S");
  //TMatrixD *bgcov = tfit->GetCovarianceMatrix(); 
  Double_t bgerr = tfit->IntegralError(ebound_l,ebound_h,s->GetParams(),s->GetCovarianceMatrix().GetMatrixArray())*binfac;

  //get expected elastics (elastic divergence from bg begin: ebound_l, end: ebound_h)
  Double_t bgint = bg->Integral(ebound_l,ebound_h)*binfac;
  //Double_t bgerr = bg->IntegralError(ebound_l,ebound_h,bg->GetParameters(),bgcov->GetMatrixArray())*binfac;
  Double_t W2nocutint = hW2->Integral(W2elasb,W2elase);
  Double_t W2nocutinterr = sqrt(W2nocutint);
  Double_t W2elas = W2nocutint - bgint;
  Double_t serr = sqrt( pow(bgerr,2) + pow(W2nocutinterr,2) );

  //Add a legend to the canvas
  auto nocutlegend = new TLegend(0.1,0.6,0.5,0.9);
  //nocutlegend->SetTextSize( 0.03 );
  nocutlegend->AddEntry( hW2elas, "Tight Elastic Cut", "l");
  nocutlegend->AddEntry( bg, "Interpolated Background (scaled)", "l");
  nocutlegend->AddEntry( tfit, "Total fit", "l");
  nocutlegend->AddEntry( (TObject*)0, "", "");
  nocutlegend->AddEntry( (TObject*)0, Form("Number of Elastics Expected: %d",(Int_t)W2elas), "");
  nocutlegend->Draw();

  c1->cd(2);
  
  //HCal anticut to obtain missed elastics
  hW2elasres = (TH1D*)(hW2_allcut->Clone("hW2elasres"));
  TF1 *tfitres = new TF1("tfitres",W2totalres,0,W2fitmax,8);
  tfitres->SetNpx(10000);  // Default is usually 100
  TF1 *bgres = new TF1("bgres",fits::g_p6fit,0.,W2fitmax,7);
  tfitres->SetLineColor(kGreen);
  hW2res->SetTitle("W^{2}, Acceptance Cut and HCal best cluster elastic anticut");
  hW2res->Fit("tfitres","RBM");

  Double_t *tparres = tfitres->GetParameters();

  //Get fit parameters for bg function and draw identical function on canvas
  bgres->SetParameters(&tparres[1]);
  bgres->SetLineColor(kRed);
  bgres->SetFillColor(kRed);
  bgres->SetFillStyle(3005);
  bgres->Draw("same");
  hW2elasres->SetLineColor(kBlue);
  hW2elasres->SetLineWidth(2);
  hW2elasres->SetFillColorAlpha(kBlue,0.35);
  hW2elasres->SetFillStyle(3003);
  hW2elasres->Draw("same");

  //get error
  TFitResultPtr r = hW2res->Fit("tfitres","S");
  //TMatrixD *bgrescov = r->GetCovarianceMatrix();
  Double_t bgreserr = tfitres->IntegralError(ebound_l,ebound_h,r->GetParams(),r->GetCovarianceMatrix().GetMatrixArray())*binfac;

  //get elastics missing from hcal
  Double_t bgresint = bgres->Integral(ebound_l,ebound_h)*binfac; //from fit
  //Double_t bgreserr = bgres->IntegralError(ebound_l,ebound_h,bgres->GetParameters(),bgrescov->GetMatrixArray())*binfac;
  Double_t W2resint = hW2res->Integral(W2elasb,W2elase); //from histo
  Double_t W2reserr = sqrt(W2resint);
  Double_t W2elasres = W2resint - bgresint;
  Double_t reserr = sqrt( pow(bgreserr,2) + pow(W2reserr,2) );
  
  Double_t effres = ( (W2elas-W2elasres) / W2elas )*100.;
  Double_t effreserr = effres*sqrt(pow(serr/W2elas,2)+pow(reserr/W2elasres,2));
  Double_t effhist =  ( (W2fullint-W2resfullint) / W2elas )*100.; //Nonsense - BG will not remain the same between W2nocut and W2residuals
  //Double_t efferr = sqrt( pow(sqrt(effres/100.*(1-effres/100.)/W2elas)*100.,2) + effreserr);
  Double_t efferr = sqrt(effres/100.*(1-effres/100.)/W2elas)*100.;

  //Add a legend to the canvas
  auto reslegend = new TLegend(0.1,0.6,0.5,0.9);
  //reslegend->SetTextSize( 0.03 );
  reslegend->AddEntry( hW2elasres, "Tight Elastic Cut", "l");
  reslegend->AddEntry( bgres, "Interpolated Background (scaled)", "l");
  reslegend->AddEntry( tfitres, "Total fit", "l");
  reslegend->AddEntry( (TObject*)0, "", "");
  reslegend->AddEntry( (TObject*)0, Form("Number of Elastics Detected: %d",(Int_t)(W2elas-W2elasres)), "");
  reslegend->AddEntry( (TObject*)0, Form("HDE via residuals: %0.3f%% +/- %0.3f%%",effres,efferr), "");
  // reslegend->AddEntry( (TObject*)0, "", "");
  // reslegend->AddEntry( (TObject*)0, Form("HDE via histo subtraction: %0.3f%%",effhist), "");
  reslegend->Draw();

  c1->Write();

  ///////////////////////////////
  //DX ELASTIC INTERPOLATE METHOD

  //Make canvas for hcal dx detected / expected 
  TCanvas *c2 = new TCanvas("c2",Form("HDE v2 sbs%d mag%d",kine,magset),1200,500);
  gStyle->SetPalette(53);
  c2->Divide(2,1);
  c2->cd(1);

  //Draw the expected elastic extraction first
  hW2->Draw();
  bg->Draw("same");
  hW2elas->Draw("same");

  auto wlegend = new TLegend(0.1,0.7,0.5,0.9);
  wlegend->AddEntry( hW2elas, "Tight Elastic Cut", "l");
  wlegend->AddEntry( bg, "Interpolated Background (scaled)", "l");
  wlegend->AddEntry( tfit, "Total fit", "l");
  wlegend->AddEntry( (TObject*)0, "", "");
  wlegend->AddEntry( (TObject*)0, Form("Number of Elastics Expected: %d",(Int_t)W2elas), "");
  wlegend->Draw();

  c2->cd(2);

  Double_t lowdxcut = dxmean-3*dxsigma;
  Double_t highdxcut = dxmean+3*dxsigma;
  // Double_t lowdxcut = dx_fitlow;
  // Double_t highdxcut = dx_fithigh;
  Double_t lowdxrange = dxmean-3*dxsigma;
  Double_t highdxrange = dxmean+3*dxsigma;
  Int_t dxelasb = hcalfit_l*hbinfac;
  Int_t dxelase = hcalfit_h*hbinfac;

  hdx_2->SetMinimum(0.0);
  hdx_2->Draw();

  //hdx_2->GetXaxis()->SetRangeUser(lowdxcut,highdxcut);
  //hdx_2->GetXaxis()->SetRangeUser(-0.3,0.3);
  //hdx_2->GetXaxis()->SetRangeUser(-0.5,0.5);
  hdx_2->GetXaxis()->SetRangeUser(dx_fitlow,dx_fithigh);

  //Then draw the dx interpolated signal fit
  hdxelastic = (TH1D*)(hdx_allcut->Clone("hdxelastic"));
  // TF1 *tfitdx = new TF1("tfitdx",dxtotal,hcalfit_l,hcalfit_h,4); //tfit npar = 1+pNfit_npar+1
  // TF1 *bgdx = new TF1("bgdx",fits::g_p2fit,hcalfit_l,hcalfit_h,3); //poly pN fit Nparam = N+1
  TF1 *tfitdx = new TF1("tfitdx",dxtotal,dx_fitlow,dx_fithigh,8); //tfit npar = 1+pNfit_npar+1
  TF1 *bgdx = new TF1("bgdx",fits::g_p6fit,dx_fitlow,dx_fithigh,7); //poly pN fit Nparam = N+1
  tfitdx->SetLineColor(kGreen);
  hdx_2->SetTitle("dx best cluster, Acceptance Cut Only");
  hdx_2->Fit("tfitdx","RBM");

  Double_t *tpardx = tfitdx->GetParameters();

  //Get fit parameters for bg function and draw identical function on canvas
  bgdx->SetParameters(&tpardx[1]);
  bgdx->SetLineColor(kRed);
  bgdx->SetFillColor(kRed);
  bgdx->SetFillStyle(3005);
  bgdx->Draw("same");
  hdxelastic->SetLineColor(kBlue);
  hdxelastic->SetLineWidth(2);
  hdxelastic->SetFillColorAlpha(kBlue,0.35);
  hdxelastic->SetFillStyle(3003);
  hdxelastic->Draw("same");

  // TLine *line1 = new TLine(lowdxcut, hdx_2->GetMinimum(), lowdxcut, hdx_2->GetMaximum());
  // line1->SetLineColor(kRed);
  // line1->SetLineWidth(2);
  // line1->Draw("same");

  // TLine *line2 = new TLine(highdxcut, hdx_2->GetMinimum(), highdxcut, hdx_2->GetMaximum());
  // line2->SetLineColor(kRed);
  // line2->SetLineWidth(2);
  // line2->Draw("same");

  //Get elastics detected in hcal
  Int_t dxbinmin = hdx_2->FindBin(lowdxcut);
  Int_t dxbinmax = hdx_2->FindBin(highdxcut);

  //get error params

  TFitResultPtr t = hdx_2->Fit("tfitdx","S");
  //TMatrixD *bgcov = tfit->GetCovarianceMatrix(); 
  //Double_t dxbgerr = tfitdx->IntegralError(hcalfit_l,hcalfit_h,t->GetParams(),t->GetCovarianceMatrix().GetMatrixArray())*hbinfac;
  Double_t dxbgerr = tfitdx->IntegralError(lowdxcut,highdxcut,t->GetParams(),t->GetCovarianceMatrix().GetMatrixArray())*hbinfac;

  cout << "dxbgerr: " << dxbgerr << endl;

  Double_t bgdxint = bgdx->Integral(lowdxcut,highdxcut)*hbinfac;

  // Double_t int_C = bgdx->Integral(hcalfit_l,hcalfit_h);
  // Double_t int_D = bgdx->Integral(hcalfit_l,hcalfit_h)*hbinfac;
  // Double_t int_E = bgdx->Integral(lowdxcut,highdxcut);
  // Double_t int_F = bgdx->Integral(lowdxcut,highdxcut)*hbinfac;
  // Double_t int_G = hdx_2->Integral();
  // Double_t int_H = hdx_2->Integral(lowdxcut,highdxcut);
  // Double_t int_I = hdx_2->Integral(dxbinmin,dxbinmax);

  Double_t int_direct = 0;
  for( int i=dxbinmin; i<=dxbinmax; i++ ){
    int_direct+=hdx_2->GetBinContent(i);
  }

  //Double_t bgdxint = bgdx->Integral(lowdxcut,highdxcut)*hbinfac;
  Double_t dxnocutint = int_direct;
  //Double_t dxnocutint = hdx_2->Integral(dxbinmin,dxbinmax);
  Double_t dxerr = sqrt(dxnocutint);

  cout << "dxerr: " << dxerr << endl;

  Double_t dxTerr = sqrt( pow(dxbgerr,2) + pow(dxerr,2) );


  Double_t dxelas = dxnocutint - bgdxint;  //hcal elastics dx direct

  Double_t effdx =  ( dxelas / W2elas )*100.; 

  Double_t effdxerr = effdx*sqrt(pow(dxbgerr/dxelas,2)+pow(dxerr/dxTerr,2));

  //Double_t efferr2 = sqrt( pow(sqrt(effdx/100.*(1-effdx/100.)/W2elas)*100.,2) + effdxerr);
  Double_t efferr2 = sqrt(effdx/100.*(1-effdx/100.)/W2elas)*100.;


  // cout << "///////////////////////////////" << endl;
  // cout << "Let's figure this out:" << endl;
  // cout << "bgdx fit values!" << endl;
  // cout << "hbinfac: " << hbinfac << endl;
  // cout << "hcalfit_l: " << hcalfit_l << endl;
  // cout << "hcalfit_h: " << hcalfit_h << endl;
  // cout << "lowdxcut: " << lowdxcut << endl;
  // cout << "highdxcut: " << highdxcut << endl;
  // cout << "dxbinmin: " << dxbinmin << endl;
  // cout << "dxbinmax: " << dxbinmax << endl;
  // cout << "bgdx full integral on full range: " << int_C << endl;
  // cout << "bgdx full integral on full range with binfac: " << int_D << endl;
  // cout << "bgdx full integral on restricted range: " << int_E << endl;
  // cout << "bgdx full integral on restricted range with binfac: " << int_F << endl;
  // cout << endl;
  // cout << "dxnocutint TH1D values!" << endl;
  // cout << "dxelasb: " << dxelasb << endl;
  // cout << "dxelase: " << dxelase << endl;
  // cout << "hdx_2 full integral: " << int_G << endl;
  // cout << "hdx_2 full integral on full range: " << dxnocutint << endl;
  // cout << "hdx_2 full integral on restricted range: " << int_H << endl;
  // cout << "hdx_2 full integral on FindBin restricted range: " << int_I << endl;
  // cout << "hdx_2 full integral on restricted range manual: " << int_direct << endl;
  
  //Add a legend to the canvas
  auto dxlegend = new TLegend(0.1,0.7,0.5,0.9);
  dxlegend->AddEntry( hdxelastic, "Tight Elastic Cut", "l");
  dxlegend->AddEntry( bgdx, Form("Interpolated Background (scaled): %d",(Int_t)bgdxint), "l");
  dxlegend->AddEntry( tfitdx, "Total fit", "l");
  dxlegend->AddEntry( (TObject*)0, "", "");
  dxlegend->AddEntry( (TObject*)0, Form("Total Events on Range: %d",(Int_t)dxnocutint), "");
  dxlegend->AddEntry( (TObject*)0, Form("Number of Elastics Detected: %d",(Int_t)dxelas), "");
  dxlegend->AddEntry( (TObject*)0, Form("HDE via dx: %0.3f%% +/- %0.3f%%",effdx,efferr2), "");
  dxlegend->Draw();

  c2->Write();

  /////////////////////
  //DX SIDEBAND METHOD

  //Make canvas for hcal dx detected / expected 
  TCanvas *c3 = new TCanvas("c3",Form("HDE sideband sbs%d mag%d",kine,magset),1200,500);
  gStyle->SetPalette(53);
  c3->Divide(2,1);
  c3->cd(1);

  //Draw the expected elastic extraction first
  hW2->Draw();
  bg->Draw("same");
  hW2elas->Draw("same");

  auto sbwlegend = new TLegend(0.1,0.7,0.5,0.9);
  sbwlegend->AddEntry( hW2elas, "Tight Elastic Cut", "l");
  sbwlegend->AddEntry( bg, "Interpolated Background (scaled)", "l");
  sbwlegend->AddEntry( tfit, "Total fit", "l");
  sbwlegend->AddEntry( (TObject*)0, "", "");
  sbwlegend->AddEntry( (TObject*)0, Form("Number of Elastics Expected: %d",(Int_t)W2elas), "");
  sbwlegend->Draw();

  c3->cd(2);

  // Perform side-band analysis
  auto [totalEvents, polyParams] = util::performSideBandAnalysis(hdx_3, lowdxcut, highdxcut);

  // Create a TF1 with the pol4 parameters
  TF1* sbbgdx = new TF1("sbbgdx", "pol4", hdx_3->GetXaxis()->GetXmin(), hdx_3->GetXaxis()->GetXmax());
  sbbgdx->SetParameters(polyParams.data()); // Set all parameters in one call

  hdx_3->SetMinimum(0.0);
  hdx_3->Draw("E");

  sbbgdx->SetLineColor(kRed);
  sbbgdx->Draw("SAME");

  // Retrieve the minimum and maximum y-values of the histogram's y-axis for line placement
  double yMin = hdx_3->GetMinimum();
  double yMax = hdx_3->GetMaximum();

  // Create lines at xLow and xHigh
  TLine* lineLow = new TLine(lowdxcut, yMin, lowdxcut, yMax);
  TLine* lineHigh = new TLine(highdxcut, yMin, highdxcut, yMax);

  // Set line color to kMagenta
  lineLow->SetLineColor(kMagenta);
  lineHigh->SetLineColor(kMagenta);

  // Draw the lines on the same canvas
  lineLow->Draw("SAME");
  lineHigh->Draw("SAME");

  TLegend* sblegend = new TLegend(0.7, 0.7, 0.89, 0.89); // Adjust these coordinates as needed
  sblegend->SetBorderSize(1);
  sblegend->SetFillColor(0);
  
  // Add the histogram, fit function, and lines to the legend
  sblegend->AddEntry( hdx_3, "Data", "lpe" );
  sblegend->AddEntry( sbbgdx, "Background Fit", "l" );
  sblegend->AddEntry( lineLow, "xLow", "l" );
  sblegend->AddEntry( lineHigh, "xHigh", "l" );
  sblegend->AddEntry( (TObject*)0, "", "" );
  sblegend->AddEntry( (TObject*)0, Form("Number of Elastics Detected: %d",(Int_t)totalEvents), "" );

  // Draw the legend
  sblegend->Draw();
  c3->Write();

  ///////////////////////////
  //ACCEPTANCE MATCHING CHECK

  //Make canvas for hcal dx detected / expected 
  TCanvas *c4 = new TCanvas("c4",Form("HDE acc match sbs%d mag%d",kine,magset),1200,500);
  gStyle->SetPalette(53);
  c4->Divide(2,1);
  c4->cd(1);

  // Create lines for hcal boundary
  // Double_t top = (econst::hcalposYi_p0+econst::hcalblk_w_p0);
  // Double_t bottom = (econst::hcalposYf_p0-econst::hcalblk_w_p0);
  // Double_t left = (econst::hcalposXf_p0-econst::hcalblk_h_p0);
  // Double_t right = (econst::hcalposXi_p0+econst::hcalblk_h_p0); 

  Double_t left = (econst::hcalposYi_p0+econst::hcalblk_w_p0);
  Double_t right = (econst::hcalposYf_p0-econst::hcalblk_w_p0);
  Double_t top = (econst::hcalposXf_p0-econst::hcalblk_h_p0);
  Double_t bottom = (econst::hcalposXi_p0+econst::hcalblk_h_p0);  

  TLine* lineTop = new TLine(left, top, right, top);
  lineTop->SetLineColor(kGray);
  TLine* lineBottom = new TLine(left, bottom, right, bottom);
  lineBottom->SetLineColor(kGray);
  TLine* lineLeft = new TLine(left, bottom, left, top);
  lineLeft->SetLineColor(kGray);
  TLine* lineRight = new TLine(right, bottom, right, top);
  lineRight->SetLineColor(kGray);

  // Create lines for acceptance to match
  // Double_t topAcc = (econst::hcalposYi_p0+econst::hcalblk_w_p0-dysig);
  // Double_t bottomAcc = (econst::hcalposYf_p0-econst::hcalblk_w_p0+dysig);
  // Double_t leftAcc = (econst::hcalposXf_p0-econst::hcalblk_h_p0+dxsig_p-dx0_p);
  // Double_t rightAcc = (econst::hcalposXi_p0+econst::hcalblk_h_p0-dxsig_p-dx0_p);

  Double_t leftAcc = (econst::hcalposYi_p0+econst::hcalblk_w_p0-dysig);
  Double_t rightAcc = (econst::hcalposYf_p0-econst::hcalblk_w_p0+dysig);
  Double_t topAcc = (econst::hcalposXf_p0-econst::hcalblk_h_p0+dxsig_p-dx0_p);
  Double_t bottomAcc = (econst::hcalposXi_p0+econst::hcalblk_h_p0-dxsig_p-dx0_p);

  TLine* lineTopAcc = new TLine(leftAcc, topAcc, rightAcc, topAcc);
  lineTopAcc->SetLineColor(kGreen);
  TLine* lineBottomAcc = new TLine(leftAcc, bottomAcc, rightAcc, bottomAcc);
  lineBottomAcc->SetLineColor(kGreen);
  TLine* lineLeftAcc = new TLine(leftAcc, bottomAcc, leftAcc, topAcc);
  lineLeftAcc->SetLineColor(kGreen);
  TLine* lineRightAcc = new TLine(rightAcc, bottomAcc, rightAcc, topAcc);
  lineRightAcc->SetLineColor(kGreen);

  // Draw the histogram
  hxyexp_nocut->Draw();

  // Draw the lines on top of the histogram
  lineTop->Draw("SAME");
  lineBottom->Draw("SAME");
  lineLeft->Draw("SAME");
  lineRight->Draw("SAME");
  lineTopAcc->Draw("SAME");
  lineBottomAcc->Draw("SAME");
  lineLeftAcc->Draw("SAME");
  lineRightAcc->Draw("SAME");

  auto accleg = new TLegend(0.1,0.7,0.5,0.9);
  accleg->AddEntry( lineTop, "HCal Boundary", "l");
  accleg->AddEntry( lineTopAcc, "Projected Acceptance", "l");
  accleg->Draw();

  c4->cd(2);

  hxyexp_aacut->Draw();

  // Draw the lines on top of the histogram
  lineTop->Draw("SAME");
  lineBottom->Draw("SAME");
  lineLeft->Draw("SAME");
  lineRight->Draw("SAME");
  lineTopAcc->Draw("SAME");
  lineBottomAcc->Draw("SAME");
  lineLeftAcc->Draw("SAME");
  lineRightAcc->Draw("SAME");

  auto acclegcut = new TLegend(0.1,0.7,0.5,0.9);
  acclegcut->AddEntry( lineTop, "HCal Boundary", "l");
  acclegcut->AddEntry( lineTopAcc, "Projected Acceptance", "l");
  acclegcut->Draw();

  fout->Write();

  cout << "Analysis complete. Outfile written to " << outfilename << endl;

  // Send time efficiency report to console
  st->Stop();
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
