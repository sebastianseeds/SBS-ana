//SSeeds 9.1.23 Script to extract the p/n yield ratio via data/MC comparison. Intent is to loop over all data from a given set/field setting, apply elastic cuts, and build a dx histogram. From this point, read in MC (with RC) distributions for quasi-elastic protons and neutrons (independently), fit these distributions with several distributions to get a best fit, then compose a sum of these fits allowing for a scaling parameter which maps to the total MC yields. Next, do the same with g4sbs inelastic generator and add this to the total fit function (sum) and apply the now three floating parameter fit to the data and check residuals and chi-square. Varying the fit function may yield better results, but the ratio will be extracted. See gmn_calc.C for application of this ratio to gmn extraction.

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

const int atimeSigFac = 5;
const int atimediffSigFac = 5;
const double dx_fitlow = -0.3;
const double dx_fithigh = 0.3;

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
//// 5) Cluster which passes atimeSigFac sigma cut on adc time, then minimizes proton theta pq
//// 6) Cluster which passes atimeSigFac sigma cut on adc time, then minimizes dxdy distance from either p dx mean or n dx mean
void gmn_ratio( Int_t kine=4, Int_t magset=50, Int_t cluster_idx=6)
{ //main  

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // //Some comparison specific variables
  // const Double_t pcsigfac = 0.58; //sbs4
  // const Double_t sampfrac = 0.0641; //sbs4
  // const Double_t dxhmin = -10; //hcal, GeV
  // const Double_t dxhmax = 10; //hcal, GeV
  // const Double_t dyhmin = -10; //hcal, GeV
  // const Double_t dyhmax = 10; //hcal, GeV
  // const Int_t nfac = 100; //Number of bins for hcal E vs nucleon p fits obtained to ensure 1000 events per bin

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/ratio.json");

  std::string rootfile_dir = jmgr->GetValueFromSubKey_str( "rootfile_dir", Form("sbs%d",kine) );
  Double_t W2fitmax = jmgr->GetValueFromSubKey<Double_t>( "W2fitmax", Form("sbs%d",kine) );
  Double_t binfac = jmgr->GetValueFromSubKey<Double_t>( "binfac", Form("sbs%d",kine) );
  Double_t hbinfac = jmgr->GetValueFromSubKey<Double_t>( "hbinfac", Form("sbs%d",kine) ); //gets to 7cm HCal pos res
  Double_t thetapqcut = jmgr->GetValueFromSubKey<Double_t>( "thetapqcut", Form("sbs%d",kine) );
  std::string target = jmgr->GetValueFromKey_str( "target" ); //Must always be LD2 for GMn analysis
  Int_t epm = jmgr->GetValueFromSubKey<Int_t>( "epm", Form("sbs%d",kine) );
  Int_t pass = jmgr->GetValueFromSubKey<Int_t>( "pass", Form("sbs%d",kine) );
  Double_t ebound_l = jmgr->GetValueFromSubKey<Double_t>( "ebound_l", Form("sbs%d",kine) );
  Double_t ebound_h = jmgr->GetValueFromSubKey<Double_t>( "ebound_h", Form("sbs%d",kine) );

  if( pass>1 ){
    std::cout << "As of 3.31.23, the highest GMn replay pass is 1. Enter a valid pass." << endl;
    return;
  }

  //set up default parameters for all analysis
  std::string runsheet_dir = "/w/halla-scshelf2102/sbs/seeds/ana/data"; //unique to my environment for now
  Int_t nruns = -1; //Always analyze all available runs
  Int_t verb = 0; //Don't print diagnostic info by default
  // Double_t hcalfit_l = econst::hcalposXi_mc-2*econst::hcalblk_div_v; //lower fit/bin limit for hcal dx plots (m)
  // Double_t hcalfit_h = econst::hcalposXf_mc+2*econst::hcalblk_div_v; //upper fit/bin limit for hcal dx plots (m)
  Double_t hcalfit_l = econst::hcalposXi_mc; //lower fit/bin limit for hcal dx plots (m)
  Double_t hcalfit_h = econst::hcalposXf_mc; //upper fit/bin limit for hcal dx plots (m)
  Double_t harmrange = econst::hcal_vrange; //Full range of hcal dx plots (m)

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<crun> corun; 
  util::ReadRunList(runsheet_dir,nruns,kine,target.c_str(),pass,verb,corun); //modifies nruns to be very large when -1

  //set up output files
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  //std::string fout_path = outdir_path + Form("/gmn_analysis/gmn_ratio_out_sbs%d_mag%d_clusteridx%d.root",kine,magset,cluster_idx);
  std::string fout_path = "test.root";
  TFile *fout = new TFile( fout_path.c_str(), "RECREATE" );

  ////////////
  //HISTOGRAMS

  //Basic
  TH2D *hdxdy_nocut = new TH2D("hdxdy_nocut","HCal dxdy, no cut; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hdxdy_cut = new TH2D("hdxdy_cut","HCal dxdy, all cuts; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hdxdy_cut_newclus = new TH2D("hdxdy_cut_newclus","HCal dxdy, all cuts, new cluster selection; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );

  TH1D *hatime_nocut = new TH1D("hatime_nocut","HCal ADC time, no cut; ns",400,-200,200);
  TH1D *hatime_cut = new TH1D("hatime_cut","HCal ADC time, all cut, no coin; ns",400,-200,200);
  TH1D *hatime_cut_newclus = new TH1D("hatime_cut_newclus","HCal ADC time, all cut, best cluster; ns",400,-200,200);

  TH1D *hatimediff_nocut = new TH1D("hatimediff_nocut","HCal ADC time - BBCal time, no cut; ns",400,-200,200);
  TH1D *hatimediff_cut = new TH1D("hatimediff_cut","HCal ADC time - BBCal time, all cut, no coin; ns",400,-200,200);
  TH1D *hatimediff_cut_newclus = new TH1D("hatimediff_cut_newclus","HCal ADC time - BBCal time, all cut, best cluster; ns",400,-200,200);
  
  TH1D *hW2_nocut = new TH1D( "hW2_nocut", "W^{2}, no cut; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_cut = new TH1D( "hW2_cut", "W^{2}, accmatch/global/atime cuts; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_cut_newclus = new TH1D( "hW2_cut_newclus", "W^{2}, accmatch/global/atime cuts, best cluster; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );

  TH1D *hdx_nocut = new TH1D( "hdx_nocut", "dx, no cut;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);
  TH1D *hdx_cut = new TH1D( "hdx_cut", "dx, e-arm cut, no dy;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);
  TH1D *hdx_allcut = new TH1D( "hdx_allcut", "dx, e-arm/dy cut;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);
  TH1D *hdx_allcut_newclus = new TH1D( "hdx_allcut_newclus", "dx, e-arm/dy cut, new cluster selection;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);

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
  Int_t passall_out;
  Int_t cidx_best_atime_out;
  Int_t cidx_best_tpq_out;
  Int_t cidx_best_dxdy_out;
  Int_t cidx_best_tpqatime_out;
  Double_t bc_atime_out;
  Double_t bc_atimediff_out;
  Double_t bc_hcale_out;
  Double_t bc_thpq_out;
  Double_t bc_dx_out;
  Double_t bc_dy_out;

  //cluster tree vars
  Double_t cpblkid_out[econst::maxclus];
  Double_t chatime_out[econst::maxclus];
  Double_t chatimediff_out[econst::maxclus];
  Double_t chcale_out[econst::maxclus];
  Double_t cthetapq_p_out[econst::maxclus];
  Double_t cthetapq_n_out[econst::maxclus];
  Double_t cdx_out[econst::maxclus];
  Double_t cdy_out[econst::maxclus];
  Int_t cpid_out[econst::maxclus]; //-1:neither,0:ambiguous,1:proton,2:neutron
  Int_t cfailedcoin_out[econst::maxclus];
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
  P->Branch( "passall", &passall_out, "passall/I" );
  P->Branch( "cidx_best_atime", &cidx_best_atime_out, "cidx_best_atime/I" );
  P->Branch( "cidx_best_tpq", &cidx_best_tpq_out, "cidx_best_tpq/I" );
  P->Branch( "cidx_best_dxdy", &cidx_best_dxdy_out, "cidx_best_dxdy/I" );
  P->Branch( "cidx_best_tpqatime", &cidx_best_tpqatime_out, "cidx_best_tpqatime/I" );
  P->Branch( "bc_atime", &bc_atime_out, "bc_atime/D" );
  P->Branch( "bc_atimediff", &bc_atimediff_out, "bc_atimediff/D" );
  P->Branch( "bc_hcale", &bc_hcale_out, "bc_hcale/D" );
  P->Branch( "bc_thpq", &bc_thpq_out, "bc_thpq/D" );
  P->Branch( "bc_dx", &bc_dx_out, "bc_dx/D" );
  P->Branch( "bc_dy", &bc_dy_out, "bc_dy/D" );

  P->Branch( "cpblkid", &cpblkid_out, Form("cpblkid[%d]/D",econst::maxclus) );
  P->Branch( "chatime", &chatime_out, Form("chatime[%d]/D",econst::maxclus) );
  P->Branch( "chatimediff", &chatimediff_out, Form("chatimediff[%d]/D",econst::maxclus) );
  P->Branch( "chcale", &chcale_out, Form("chcale[%d]/D",econst::maxclus) );
  P->Branch( "cthetapq_p", &cthetapq_p_out, Form("cthetapq_p[%d]/D",econst::maxclus) );
  P->Branch( "cthetapq_n", &cthetapq_n_out, Form("cthetapq_n[%d]/D",econst::maxclus) );
  P->Branch( "cdx", &cdx_out, Form("cdx[%d]/D",econst::maxclus) );
  P->Branch( "cdy", &cdy_out, Form("cdy[%d]/D",econst::maxclus) );
  P->Branch( "cpid", &cpid_out, Form("cpid[%d]/I",econst::maxclus) );
  P->Branch( "cfailedcoin", &cfailedcoin_out, Form("cfailedcoin[%d]/I",econst::maxclus) );
  P->Branch( "cinspot_p", &cinspot_p_out, Form("cinspot_p[%d]/I",econst::maxclus) );

  // setup reporting indices
  std::cout << std::endl << std::endl;
  Double_t dxmean;
  Double_t dxsigma;

  for (Int_t irun=0; irun<nruns; irun++) {
    // accessing run info
    Int_t runnum = corun[irun].runnum;
    Int_t mag = corun[irun].sbsmag / 21; //convert to percent
    Double_t ebeam = corun[irun].ebeam; //get beam energy per run
    std::string tar = corun[irun].target;

    //Only proceed for deuterium at desired field setting
    if( mag!=magset || tar.compare("LD2")!=0 ) 
      continue;

    std::cout << "Analyzing run " << runnum << ".." << std::endl;

    //set up configuration and tune objects to load analysis parameters
    SBSconfig config(kine,mag);

    //Obtain configuration pars from config file
    Double_t hcaltheta = config.GetHCALtheta_rad();
    Double_t hcaldist = config.GetHCALdist();
    Double_t sbsdist = config.GetSBSdist();
    Double_t bbthr = config.GetBBtheta_rad(); //in radians
    
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

    //For this analysis, only one magnetic field
    dxmean = dx0_p;
    dxsigma = dxsig_p;

    Double_t W2min = W2mean - 2*W2sig;
    Double_t W2max = W2mean + 2*W2sig;

    //Set minimum energy for expected recoil nucleon to be considered
    Double_t hminE = 0.0;  //0.05 GeV or nothing
    
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
    rvars::setbranch(C, "sbs.hcal.clus", hcalcvar, hcalcvarlink, 6);
    
    // hodoscope cluster mean time
    Int_t Nhodotmean; 
    Double_t hodotmean[econst::maxclus];
    std::vector<std::string> hodovar = {"clus.tmean","clus.tmean"};
    std::vector<void*> hodovarlink = {&Nhodotmean,&hodotmean};
    rvars::setbranch(C, "bb.hodotdc", hodovar, hodovarlink, 0); 

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
    TVector3 hcalorigin = hcaldist*hcalaxes[2] + econst::hcalvoff*hcalaxes[0];
    Double_t BdL = econst::sbsmaxfield * econst::sbsdipolegap * (mag/100); //scale crudely by magnetic field
    Double_t Eloss_outgoing = econst::celldiameter/2.0/sin(bbthr) * econst::lh2tarrho * econst::lh2dEdx;

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
    Int_t clusterBadSort = 0;

    while (C->GetEntry(nevent++)) {
      
      std::cout << "Bad Cluster Sorts: " << clusterBadSort << ". Processing run " << corun[irun].runnum << " event " << nevent << "/" << nevents << "\r";
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
      //HCal active area cut (acceptance matching). Save pass/fail for output tree.
      vector<Double_t> xyhcalexp; vars::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
      bool failedaccmatch = 
	xyhcalexp[1] > (econst::hcalposYf_mc-econst::hcalblk_div_h) ||
	xyhcalexp[1] < (econst::hcalposYi_mc+econst::hcalblk_div_h) ||
	xyhcalexp[0] > (econst::hcalposXf_mc-econst::hcalblk_div_v) ||
	xyhcalexp[0] < (econst::hcalposXi_mc+econst::hcalblk_div_v);
      
      //Get expected proton deflection for thetapq (this needs work to be useful)
      Double_t protdeflect = tan( 0.3 * BdL / qv.Mag() ) * (hcaldist - (sbsdist + econst::sbsdipolegap/2.0) );
      
      //Get some primary cluster primary block variables
      Double_t dx = hcalcx[0] - xyhcalexp[0];
      Double_t dy = hcalcy[0] - xyhcalexp[1];

      //Primary cluster atime cut
      bool failedpclusatime =  abs( hcalcatime[0] - atime0 ) > atimeSigFac*atimesig;

      //dy cut
      bool faileddy = abs(dy-dy0)>dysig;

      //////////////////////
      //ALL CLUSTER ANALYSIS
      
      //Set up indexes and variables for best cluster analysis. All index defaults are zero (highest E)
      Int_t cidx_atime = 0;
      Double_t c_atimediff = 1000.;
      Int_t cidx_e = 0;
      Double_t c_e = 0.;
      Int_t cidx_thetapq = 0;
      Double_t c_thetapqdiff = 1000.;
      Int_t cidx_dxdy = 0;
      Double_t c_dxdydiff = 1000.;
      Int_t cidx_tpq_atime = 0;
      Double_t c_tpq_atimediff = 1000.;
      Int_t cidx_dxdy_p = 0;
      Double_t c_dxdydiff_p = 1000.;
      Int_t cidx_dxdy_n = 0;
      Double_t c_dxdydiff_n = 1000.;

      Double_t tpq_p[econst::maxclus];

      //loop through all clusters
      for( Int_t c=0; c<Nhcalcid; c++ ){

	//calculate h-arm physics quantities per cluster
	Double_t atime_diff = hcalcatime[c] - atime_sh[0]; //Assuming best shower time on primary cluster
	Double_t cdx = hcalcx[c] - xyhcalexp[0];
	Double_t cdy = hcalcy[c] - xyhcalexp[1];	  
	TVector3 chcalpos = hcalorigin + hcalcx[c]*hcalaxes[0] + hcalcy[c]*hcalaxes[1];
	TVector3 cneutdir = ( chcalpos - vertex ).Unit();
	TVector3 cprotdir = ( chcalpos + protdeflect * hcalaxes[0] - vertex ).Unit();
	Double_t cthetapq_p = acos( cprotdir.Dot( qv.Unit() ) );
	tpq_p[c]=cthetapq_p;
	Double_t cthetapq_n = acos( cneutdir.Dot( qv.Unit() ) );
	Int_t cpid;
	bool cis_p = abs(cdx - dx0_p)<3*dxsig_p && abs(cdy - dy0)<3*dysig;
	bool cis_n = abs(cdx - dx0_n)<3*dxsig_n && abs(cdy - dy0)<3*dysig;
	bool c_failedcoin = abs( hcalcatime[c] - atime0 ) > atimeSigFac*atimesig;
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
	//get cluster which maximizes energy (METHOD 1)
	if( hcalce[c]>c_e ){
	  cidx_e = c;
	  c_e = hcalce[c];
	}
	//get cluster which minimizes adct diff from signal (METHOD 2)
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
	if( !c_failedcoin && cthetapq_p<c_tpq_atimediff ){
	  cidx_tpq_atime = c;
	  c_tpq_atimediff = cthetapq_p; 
	}

	//get cluster which passes wide cut on atime, then minimizes dxdy for either p or n hypothesis (METHOD 6)
	Double_t dxdydist_p = sqrt( pow(cdx-dx0_p,2)+pow(cdy-dy0,2));
	Double_t dxdydist_n = sqrt( pow(cdx-dx0_n,2)+pow(cdy-dy0,2));
	if( !c_failedcoin && dxdydist_p<c_dxdydiff_p ){
	  cidx_dxdy_p = c;
	  c_dxdydiff_p = dxdydist_p;
	}
	if( !c_failedcoin && dxdydist_n<c_dxdydiff_n ){
	  cidx_dxdy_n = c;
	  c_dxdydiff_n = dxdydist_n;
	}

	//fill analysis tree variables
	chcale_out[c] = hcalce[c];
	cpblkid_out[c] = hcalcid[c];
	chatime_out[c] = hcalcatime[c];
	chatimediff_out[c] = hcalcatime[c] - atime_sh[0];
	cthetapq_p_out[c] = cthetapq_p;
	cthetapq_n_out[c] = cthetapq_n;
	cdx_out[c] = cdx;
	cdy_out[c] = cdy;
	cpid_out[c] = cpid;
	cfailedcoin_out[c] = (Int_t)c_failedcoin;
	cinspot_p_out[c] = (Int_t)c_spotcheck_p;
       
      }

      if( cidx_e!=0 )
	clusterBadSort++;

      //Switch between best clusters for systematic analysis
      Int_t cidx_best;
      
      switch (cluster_idx) {
      case 1:
	cidx_best = cidx_e;
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
	if( c_dxdydiff_p < c_dxdydiff_n )
	  cidx_best = cidx_dxdy_p;
	else
	  cidx_best = cidx_dxdy_n;
	break;
      default:
	cidx_best = 0;
      }
      
      //Calculations from the best cluster for later use
      Double_t dx_bestcluster = hcalcx[cidx_best] - xyhcalexp[0];
      Double_t dy_bestcluster = hcalcy[cidx_best] - xyhcalexp[1];
      Double_t hatime_bestcluster = hcalcatime[cidx_best];
      Double_t hatimediff_bestcluster = hcalcatime[cidx_best] - atime_sh[0];
      Double_t thpq_bestcluster = tpq_p[cidx_best];
      Double_t ce_bestcluster = hcalce[cidx_best];
      bool passedatime = abs( hatime_bestcluster - atime0 ) < atimeSigFac*atimesig;
      bool passedcoin = abs( hatime_bestcluster - atimediff0 ) < atimediffSigFac*atimediffsig;
      
      //Fill some output analysis tree variables
      mag_out = mag;
      run_out = runnum;
      failedglobal_out = (Int_t)failedglobal;
      failedaccmatch_out = (Int_t)failedaccmatch;

      //Fill some nocut histograms
      hW2_nocut->Fill(W2);
      hdxdy_nocut->Fill(dy,dx);
      hdx_nocut->Fill(dx);
      hatime_nocut->Fill(hcalcatime[0]);
      hatimediff_nocut->Fill(hcalcatime[0]-atime_sh[0]);

      if( !failedpclusatime )
	hW2_cut->Fill(W2);

      if( passedatime )
	hW2_cut_newclus->Fill(W2);

      //E-arm cuts first
      if( failedglobal || failedaccmatch || failedW2 )
	continue;
    
      //Fill final output tree variables    
      cidx_best_atime_out = cidx_atime;
      cidx_best_tpq_out = cidx_thetapq;
      cidx_best_dxdy_out = cidx_dxdy;
      cidx_best_tpqatime_out = cidx_tpq_atime;
      bc_atime_out = hatime_bestcluster;
      bc_atimediff_out = hatimediff_bestcluster;
      bc_hcale_out = ce_bestcluster;
      bc_thpq_out = thpq_bestcluster;
      bc_dx_out = dx_bestcluster;
      bc_dy_out = dy_bestcluster;
      
      P->Fill();
      
      //Fill all cut timing histograms
      hatime_cut->Fill(hcalcatime[0]);
      hatimediff_cut->Fill(hcalcatime[0]-atime_sh[0]);
	  
      //Fill all cut best cluster timing histograms
      hatime_cut_newclus->Fill(hatime_bestcluster);
      hatimediff_cut->Fill(hatimediff_bestcluster);
	

      //Fill all cut histograms
      hdxdy_cut->Fill(dy,dx);
      hdx_cut->Fill(dx);
      
      //Fill all best cluster histograms
      hdxdy_cut_newclus->Fill(dy_bestcluster,dx_bestcluster);

      if( !faileddy ){
	hdx_allcut->Fill(dx);
	hdx_allcut_newclus->Fill(dx_bestcluster);
      }

    } //end event loop
    
    // reset chain for the next run config
    C->Reset();
    
  } //end run loop
  
  fout->Write();
  
  cout << endl << "Analysis complete. Outfile located at " << fout_path << endl;

  // Send time efficiency report to console
  st->Stop();
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
