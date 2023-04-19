//sseeds 03.31.23 - test script to use analysis framework to produce hcal detection efficiency output more efficiently
//04.10.23 Update - Added W2 interpolate method for obtaining background fits to both total W2 distribution and W2 with HCal anticut 
//04.11.23 Update - Added back direct dx "detected" yield method for comparison. Fixed thetapq calculation and included in elastic cuts.

#include <vector>
#include <iostream>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

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
  return signal + fits::g_p12fit(x,&par[1]);
}

//Pass only kinematic and sbs magnetic field setting. Config file for all kinematics located ../../config/shde.json
void hde( Int_t kine=4, Int_t magset=0, bool pclusonly=false )
{ //main  

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  Double_t pcsigfac = 0.58;
  Double_t sampfrac = 0.0641;

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

  if( pass>1 ){
    std::cout << "As of 3.31.23, the highest GMn replay pass is 1. Enter a valid pass." << endl;
    return;
  }

  //set up default parameters for all analysis
  std::string runsheet_dir = "/w/halla-scshelf2102/sbs/seeds/ana/data"; //unique to my environment for now
  Int_t nruns = -1; //Always analyze all available runs
  Int_t verb = 0; //Don't print diagnostic info by default
  Double_t hcalfit_l = econst::hcalposXi_mc-2*econst::hcalblk_div_v; //lower fit/bin limit for hcal dx plots (m)
  Double_t hcalfit_h = econst::hcalposXf_mc+2*econst::hcalblk_div_v; //upper fit/bin limit for hcal dx plots (m)
  Double_t harmrange = (hcalfit_h) - (hcalfit_l); //Full range of hcal dx plots (m)

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<crun> corun; 
  util::ReadRunList(runsheet_dir,nruns,kine,target.c_str(),pass,verb,corun); //modifies nruns to be very large when -1. Select run list for LH2 only for hde

  //set up output files
  TFile *fout = new TFile( Form("outfiles/hde_sbs%d_%s_epm%d_magset%d.root",kine,target.c_str(),epm,magset), "RECREATE" );

  //set up diagnostic histograms
  TH2D *hxy_nocut = new TH2D("hxy_nocut","HCal xy, no cut;dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxyexp_nocut = new TH2D("hxyexp_nocut","HCal xy expected from BB, no cut;dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hxyexp_aacut = new TH2D("hxyexp_aacut","HCal xy expected from BB, hcal active area cut;dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH1D *hbclus_atime = new TH1D("hbclus_atime","Best coin cluster idx", 10, 0, 10 );
  TH1D *hbclus_thetapq = new TH1D("hbclus_thetapq","Best thetapq cluster idx", 10, 0, 10 );
  TH1D *hbclus_dxdy = new TH1D("hbclus_dxdy","Best dxdy cluster idx", 10, 0, 10 );

  //set up detection efficiency histograms
  TH1D *hW2_allcut = new TH1D( "hW2_allcut", "W^{2}, global/dy/dx/thetapq cut;W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_anticut = new TH1D( "hW2_anticut", "W^{2}, elastic anticut (global e-arm and hcal dy/dx);W^{2} (GeV^{2});", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_nocut = new TH1D( "hW2_nocut", "W^{2}, no cut;W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hdx_allcut = new TH1D( "hdx_allcut", "dx, W^{2} and dy cut;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);
  TH1D *hdx_nocut = new TH1D( "hdx_nocut", "dx, no cut;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);

  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;

  // create output tree
  TTree *P = new TTree("P","Analysis Tree"); 

  // output tree vars
  Double_t dx_out;
  Double_t dy_out;
  Double_t W2_out;
  Double_t Q2_out;
  Double_t nu_out;
  Double_t ep_out; //track reconstructed e' momentum
  Double_t hcale_out;
  Double_t hcaleexp_out;
  Double_t hcalatime_out;
  Double_t thetapq_pout;
  Double_t thetapq_nout;
  Int_t pblkid_out;
  Int_t pid_out; //-1:neither,0:ambiguous,1:proton,2:neutron
  Int_t mag_out;
  Int_t run_out; //run number
  Int_t failedglobal_out;
  Int_t failedaccmatch_out;
  Int_t failedcoin_out;

  //cluster tree vars
  Double_t cpblkid_out[econst::maxclus];
  Double_t chatime_out[econst::maxclus];
  Double_t cthetapq_p_out[econst::maxclus];
  Double_t cthetapq_n_out[econst::maxclus];
  Double_t cdx_out[econst::maxclus];
  Double_t cdy_out[econst::maxclus];
  Int_t cpid_out[econst::maxclus];

  // set output tree branches
  P->Branch( "dx", &dx_out, "dx/D" );
  P->Branch( "dy", &dy_out, "dy/D" );
  P->Branch( "W2", &W2_out, "W2/D" );
  P->Branch( "Q2", &Q2_out, "Q2/D" );
  P->Branch( "nu", &nu_out, "nu/D" );
  P->Branch( "hcale", &hcale_out, "hcale/D" );
  P->Branch( "hcaleexp", &hcaleexp_out, "hcaleexp/D" );
  P->Branch( "hcalatime", &hcalatime_out, "hcalatime/D" );
  P->Branch( "thetapq_p", &thetapq_pout, "thetapq_p/D" );
  P->Branch( "thetapq_n", &thetapq_nout, "thetapq_n/D" );
  P->Branch( "pblkid", &pblkid_out, "pblkid/I" );
  P->Branch( "pid", &pid_out, "pid_out/I" );
  P->Branch( "mag", &mag_out, "mag/I" );
  P->Branch( "run", &run_out, "run_out/I" );
  P->Branch( "failedglobal", &failedglobal_out, "failedglobal/I" );
  P->Branch( "failedaccmatch", &failedaccmatch_out, "failedaccmatch/I" );
  P->Branch( "failedcoin", &failedcoin_out, "failedcoin/I" );

  P->Branch( "cpblkid", &cpblkid_out, Form("cpblkid[%d]/D",econst::maxchan) );
  P->Branch( "chatime", &chatime_out, Form("chatime[%d]/D",econst::maxchan) );
  P->Branch( "cthetapq_p", &cthetapq_p_out, Form("cthetapq_p[%d]/D",econst::maxchan) );
  P->Branch( "cthetapq_n", &cthetapq_n_out, Form("cthetapq_n[%d]/D",econst::maxchan) );
  P->Branch( "cdx", &cdx_out, Form("cdx[%d]/D",econst::maxchan) );
  P->Branch( "cdy", &cdy_out, Form("cdy[%d]/D",econst::maxchan) );
  P->Branch( "cpid", &cpid_out, Form("cpid[%d]/I",econst::maxchan) );

  // setup reporting indices
  Int_t curmag = -1;
  std::string curtar = "";
  std::cout << std::endl << std::endl;
  Double_t dxmean;
  Double_t dxsigma;

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
    Double_t dx0_p    = tune.Getdx0_p();
    Double_t dy0      = tune.Getdy0();
    Double_t dxsig_n  = tune.Getdxsig_n();
    Double_t dxsig_p  = tune.Getdxsig_p();
    Double_t dysig    = tune.Getdysig();
    Double_t atime0   = tune.Getatime0();
    Double_t atimesig = tune.Getatimesig();

    //For this analysis, only one magnetic field
    dxmean = dx0_p;
    dxsigma = dxsig_p;

    //Set cut/anticut limits for proton
    Double_t dxmin = dx0_p - 3*dxsig_p; //Obtain 99.5% of elastic sample
    Double_t dxmax = dx0_p + 3*dxsig_p; //Obtain 99.5% of elastic sample
    Double_t dymin = dy0 - 2*dysig;
    Double_t dymax = dy0 + 2*dysig;
    Double_t W2min = W2mean - 3*W2sig;
    Double_t W2max = W2mean + 3*W2sig;
    Double_t atmin = atime0 - 3*atimesig;
    Double_t atmax = atime0 + 3*atimesig;

    //Set minimum energy for expected recoil nucleon to be considered
    Double_t hminE = 0.05;  //GeV
    
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

    // HCal cluster branches
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

    while (C->GetEntry(nevent++)) {

      std::cout << "Processing run " << corun[irun].runnum << " event " << nevent << "/" << nevents << "\r";
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

      //if( failedglobal ) continue;
      Int_t fglobalint = 0;
      if( failedglobal ) fglobalint = 1;
      ///////
      
      ///////
      //HCal coin cut. Save pass/fail for output tree.
      bool failedcoin = abs( hcalatime - atime0 ) > 3*atimesig;

      //if( failedcoin ) continue;
      Int_t fcoinint = 0;
      if( failedcoin ) fcoinint = 1;
      ///////

      ///////
      //Physics calculations
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

      //Calculate h-arm quantities
      vector<Double_t> xyhcalexp; vars::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
      TVector3 hcalpos = hcalorigin + hcalx*hcalaxes[0] + hcaly*hcalaxes[1]; //-,+ before
      Double_t dx = hcalx - xyhcalexp[0];
      Double_t dy = hcaly - xyhcalexp[1];
      TVector3 neutdir = (hcalpos - vertex).Unit();
      Double_t protdeflect = tan( 0.3 * BdL / qv.Mag() ) * (hcaldist - (sbsdist + econst::sbsdipolegap/2.0) );
      TVector3 protdir = ( hcalpos + protdeflect * hcalaxes[0] - vertex ).Unit();
      Double_t thetapq_p = acos( protdir.Dot( qv.Unit() ) );
      Double_t thetapq_n = acos( neutdir.Dot( qv.Unit() ) );

      Double_t KE_exp = nu*sampfrac/pcsigfac; //Expected measured E in HCal

      //Get PID
      Int_t pid;
      bool is_p = abs(dx - dx0_p)<3*dxsig_p && abs(dy - dy0)<3*dysig;
      bool is_n = abs(dx - dx0_n)<3*dxsig_n && abs(dy - dy0)<3*dysig;
      if( is_p && !is_n )
	pid=1;
      else if( is_n && !is_p )
	pid=2;
      else if( is_p && is_n )
	pid=0;
      else
	pid=-1;

      //Fill diagnostic histos
      hxy_nocut->Fill( hcaly, hcalx );
      hxyexp_nocut->Fill( xyhcalexp[1], xyhcalexp[0] );

      ///////
      //HCal active area cut (acceptance matching). Save pass/fail for output tree.
      bool failedaccmatch = 
	xyhcalexp[1] > (econst::hcalposYf_mc-econst::hcalblk_div_h) ||
	xyhcalexp[1] < (econst::hcalposYi_mc+econst::hcalblk_div_h) ||
	xyhcalexp[0] > (econst::hcalposXf_mc-econst::hcalblk_div_v) ||
	xyhcalexp[0] < (econst::hcalposXi_mc+econst::hcalblk_div_v);

      //if( failedaccmatch ) continue;
      Int_t faccmatchint = 0;
      if( failedaccmatch ) faccmatchint = 1;
      ///////

      if( !failedaccmatch ){
	hxyexp_aacut->Fill( xyhcalexp[1], xyhcalexp[0] );
      }

      //Fill all block physics output - this is inefficient, should get primary cluster info from here first
      Int_t cidx_atime = 0;
      Double_t c_atimediff = 1000.;
      Int_t cidx_thetapq = 0;
      Double_t c_thetapqdiff = 1000.;
      Int_t cidx_dxdy = 0;
      Double_t c_dxdydiff = 1000.;
      Double_t tpq_best;

      for( Int_t c=0; c<Nhcalcid; c++ ){
	Double_t cdx = hcalcx[c] - xyhcalexp[0];
	Double_t cdy = hcalcy[c] - xyhcalexp[1];	  
	TVector3 chcalpos = hcalorigin + hcalcx[c]*hcalaxes[0] + hcalcy[c]*hcalaxes[1];
	TVector3 cneutdir = ( chcalpos - vertex ).Unit();
	TVector3 cprotdir = ( chcalpos + protdeflect * hcalaxes[0] - vertex ).Unit();
	Double_t cthetapq_p = acos( cprotdir.Dot( qv.Unit() ) );
	Double_t cthetapq_n = acos( cneutdir.Dot( qv.Unit() ) );
	Int_t cpid;
	bool cis_p = abs(cdx - dx0_p)<3*dxsig_p && abs(cdy - dy0)<3*dysig;
	bool cis_n = abs(cdx - dx0_n)<3*dxsig_n && abs(cdy - dy0)<3*dysig;
	if( cis_p && !cis_n )
	  cpid=1;
	else if( cis_n && !cis_p )
	  cpid=2;
	else if( cis_p && cis_n )
	  cpid=0;
	else
	  cpid=-1;

	//get best cluster indexes
	//get cluster which minimizes coin diff
	if( abs(hcalcatime[c]-atime0)<c_atimediff ){
	  cidx_atime = c;
	  c_atimediff = hcalcatime[c];
	}
	//get cluster which minimizes thetapq (proton for this analysis)
	if( cthetapq_p<c_thetapqdiff ){
	  cidx_thetapq = c;
	  c_thetapqdiff = cthetapq_p;
	}
	//get cluster which minimizes distance between cluster center and expected loc dxdy
	Double_t dxdydist = sqrt( pow(cdx-dx0_p,2)+pow(cdy-dy0,2));
	if( dxdydist<c_dxdydiff ){
	  cidx_dxdy = c;
	  c_dxdydiff = dxdydist;
	  tpq_best = cthetapq_p; //record the "best" cluster thetapq for later
	}

	cpblkid_out[c] = hcalcid[c];
	chatime_out[c] = hcalcatime[c];
	cthetapq_p_out[c] = cthetapq_p;
	cthetapq_n_out[c] = cthetapq_n;
	cdx_out[c] = cdx;
	cdy_out[c] = cdy;
	  
	cpid_out[c] = cpid;
      }
	
      //Calculate dx, dy, thetapq from the best cluster for later use (USING dxdy for now)
      Double_t dx_bestcluster = hcalcx[cidx_dxdy] - xyhcalexp[0];
      Double_t dy_bestcluster = hcalcy[cidx_dxdy] - xyhcalexp[1];
      Double_t hatime_bestcluster = hcalcatime[cidx_dxdy];
      
      hbclus_atime->Fill(cidx_atime);
      hbclus_thetapq->Fill(cidx_thetapq);
      hbclus_dxdy->Fill(cidx_dxdy);

      //Fill physics output tree     
      dx_out = dx;
      dy_out = dy;
      W2_out = W2;
      Q2_out = Q2;
      nu_out = nu;
      ep_out = precon;
      hcale_out = hcale;
      hcaleexp_out = KE_exp;
      hcalatime_out = hcalatime;
      thetapq_pout = thetapq_p;
      thetapq_nout = thetapq_n;
      pblkid_out = hcalid;
      mag_out = mag;
      run_out = runnum;
      pid_out = pid;
      failedglobal_out = fglobalint;
      failedaccmatch_out = faccmatchint;
      failedcoin_out = fcoinint;

      P->Fill();

      /////////////////////////////////////////////
      //HCal Detection Efficiency - primary cluster
      /////////////////////////////////////////////
      if( pclusonly ){
	//All analyses must first require elastic acceptance matching between e and h arms per event
	if( !failedaccmatch&&KE_exp>hminE ){
	  //Obtain W2 distribution with only acceptance matching
	  hW2_nocut->Fill( W2 );
	  hdx_nocut->Fill( dx );
	  //Make strongest possible hcal elastic anticut (failedglobal left out)
	  if( thetapq_p>0.04&&(dx<dxmin||dx>dxmax)&&(dy<dymin||dy>dymax) ){
	    hW2_anticut->Fill( W2 );
	  }
	  //And tightest elastic cut possible for elastic shape (both W2 and dx)
	  if( dy>dymin&&dy<dymax&&failedglobal==0&&hcalatime>atmin&&hcalatime<atmax ){
	    if( thetapq_p<0.04&&dx>dxmin&&dx<dxmax ){
	      hW2_allcut->Fill( W2 );
	    }
	    if( W2>W2min&&W2<W2max ){
	      hdx_allcut->Fill( dx );
	    }
	  }
	}
      }
      /////////////////////////////////////////////

      /////////////////////////////////////////////
      //HCal Detection Efficiency - "best" cluster
      /////////////////////////////////////////////
      if( !pclusonly ){
	//All analyses must first require elastic acceptance matching between e and h arms per event
	if( !failedaccmatch&&KE_exp>hminE ){
	  //Obtain W2 distribution with only acceptance matching
	  hW2_nocut->Fill( W2 );
	  hdx_nocut->Fill( dx_bestcluster );
	  //Make strongest possible hcal elastic anticut (failedglobal left out)
	  if( tpq_best>0.04&&
	      (dx_bestcluster<dxmin||dx_bestcluster>dxmax)&&
	      (dy_bestcluster<dymin||dy_bestcluster>dymax) ){
	    hW2_anticut->Fill( W2 );
	  }
	  //And tightest elastic cut possible for elastic shape (both W2 and dx)
	  if( dy_bestcluster>dymin&&dy_bestcluster<dymax&&
	      failedglobal==0&&
	      hatime_bestcluster>atmin&&
	      hatime_bestcluster<atmax ){
	    if( tpq_best<0.04&&
		dx_bestcluster>dxmin&&
		dx_bestcluster<dxmax ){
	      hW2_allcut->Fill( W2 );
	    }
	    if( tpq_best<0.04&&
		W2>W2min&&
		W2<W2max ){
	      hdx_allcut->Fill( dx_bestcluster );
	    }
	  }
	}
      }
      /////////////////////////////////////////////


    } //end event loop

    // reset chain for the next run config
    C->Reset();

  } //end run loop

  ////////////////////////////////////////////////////
  //HCal detection efficiency post-event-loop analysis
  ////////////////////////////////////////////////////

  TH1D *hW2 = (TH1D*)(hW2_nocut->Clone("hW2"));
  TH1D *hW2res = (TH1D*)(hW2_anticut->Clone("hW2res"));
  TH1D *hdx = (TH1D*)(hdx_nocut->Clone("hdx")); 

  Double_t W2fullint = hW2->Integral(0,W2fitmax);
  Double_t W2resfullint = hW2res->Integral(0,W2fitmax);

  //Make canvas for (expected-residuals)/expected
  TCanvas *c1 = new TCanvas("c1",Form("HDE sbs%d mag%d",kine,magset),1200,500);
  gStyle->SetPalette(53);
  c1->Divide(1,2);
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

  //get expected elastics (elastic divergence from bg begin: ebound_l, end: ebound_h)
  Double_t bgint = bg->Integral(ebound_l,ebound_h)*binfac;
  Int_t W2elasb = ebound_l*binfac;
  Int_t W2elase = ebound_h*binfac;
  Double_t W2nocutint = hW2->Integral(W2elasb,W2elase);
  Double_t W2elas = W2nocutint - bgint;

  //Add a legend to the canvas
  auto nocutlegend = new TLegend(0.1,0.6,0.5,0.9);
  nocutlegend->SetTextSize( 0.03 );
  nocutlegend->AddEntry( hW2elas, "Tight Elastic Cut", "l");
  nocutlegend->AddEntry( bg, "Interpolated Background (scaled)", "l");
  nocutlegend->AddEntry( tfit, "Total fit (scaled background + interpolated signal)", "l");
  nocutlegend->Draw();

  c1->cd(2);
  
  //HCal anticut to obtain missed elastics
  hW2elasres = (TH1D*)(hW2_allcut->Clone("hW2elasres"));
  TF1 *tfitres = new TF1("tfitres",W2totalres,0,W2fitmax,8);
  TF1 *bgres = new TF1("bgres",fits::g_p6fit,0.,W2fitmax,7);
  tfitres->SetLineColor(kGreen);
  if( pclusonly )
    hW2res->SetTitle("W^{2}, Acceptance Cut and HCal highest E cluster elastic anticut");
  else
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
  hW2elasres->SetFillColor(kBlue);
  hW2elasres->SetFillStyle(3003);
  hW2elasres->Draw("same");

  //get elastics missing from hcal
  Double_t bgresint = bgres->Integral(ebound_l,ebound_h)*binfac;
  Double_t W2resint = hW2res->Integral(W2elasb,W2elase);
  Double_t W2elasres = W2resint - bgresint;

  Double_t effres = ( (W2elas-W2elasres) / W2elas )*100.;
  Double_t effhist =  ( (W2fullint-W2resfullint) / W2elas )*100.; //Nonsense - BG will not remain the same between W2nocut and W2residuals
  Double_t efferr = sqrt(effres/100.*(1-effres/100.)/W2elas)*100.;

  //Add a legend to the canvas
  auto reslegend = new TLegend(0.1,0.6,0.5,0.9);
  reslegend->SetTextSize( 0.03 );
  reslegend->AddEntry( hW2elasres, "HCal dx/dy Elastic Anticut", "l");
  reslegend->AddEntry( bgres, "Interpolated Background (scaled)", "l");
  reslegend->AddEntry( tfitres, "Total fit (scaled background + interpolated signal)", "l");
  reslegend->AddEntry( (TObject*)0, "", "");
  reslegend->AddEntry( (TObject*)0, Form("HDE via residuals: %0.3f%% +/- %0.3f%%",effres,efferr), "");
  // reslegend->AddEntry( (TObject*)0, "", "");
  // reslegend->AddEntry( (TObject*)0, Form("HDE via histo subtraction: %0.3f%%",effhist), "");
  reslegend->Draw();

  c1->Write();

  //Make canvas for hcal dx detected / expected 
  TCanvas *c2 = new TCanvas("c2",Form("HDE v2 sbs%d mag%d",kine,magset),1200,500);
  gStyle->SetPalette(53);
  c2->Divide(1,2);
  c2->cd(1);

  //Draw the expected elastic extraction first
  hW2->Draw();
  bg->Draw("same");
  hW2elasres->Draw("same");

  c2->cd(2);

  //Then draw the dx interpolated signal fit
  hdxelastic = (TH1D*)(hdx_allcut->Clone("hdxelastic"));
  TF1 *tfitdx = new TF1("tfitdx",dxtotal,hcalfit_l,hcalfit_h,14); //tfit npar = 1+pNfit_npar+1
  TF1 *bgdx = new TF1("bgdx",fits::g_p12fit,hcalfit_l,hcalfit_h,13);
  tfitdx->SetLineColor(kGreen);
  if( pclusonly )
    hdx->SetTitle("dx highest E cluster, Acceptance Cut Only");
  else
    hdx->SetTitle("dx best cluster, Acceptance Cut Only");
  hdx->Fit("tfitdx","RBM");

  Double_t *tpardx = tfitdx->GetParameters();

  //Get fit parameters for bg function and draw identical function on canvas
  bgdx->SetParameters(&tpardx[1]);
  bgdx->SetLineColor(kRed);
  bgdx->SetFillColor(kRed);
  bgdx->SetFillStyle(3005);
  bgdx->Draw("same");
  hdxelastic->SetLineColor(kBlue);
  hdxelastic->SetLineWidth(2);
  hdxelastic->SetFillColor(kBlue);
  hdxelastic->SetFillStyle(3003);
  hdxelastic->Draw("same");

  //Get elastics detected in hcal
  //Double_t bgdxint = bgdx->Integral(hcalfit_l,hcalfit_h)*hbinfac;
  Double_t bgdxint = bgdx->Integral(dxmean-2*dxsigma,dxmean+2*dxsigma)*hbinfac;
  //Double_t dxnocutint = hdx->Integral(1,harmrange*hbinfac);
  Double_t dxnocutint = hdx->Integral((-hcalfit_l+dxmean-2*dxsigma)*hbinfac,(-hcalfit_l+dxmean+2*dxsigma)*hbinfac);
  Double_t dxelas = dxnocutint - bgdxint;  //hcal elastics alt

  Double_t effdx =  ( dxelas / W2elas )*100.; 

  //Add a legend to the canvas
  auto dxlegend = new TLegend(0.1,0.7,0.5,0.9);
  dxlegend->SetTextSize( 0.03 );
  dxlegend->AddEntry( hdxelastic, "Tight Elastic Cut", "l");
  dxlegend->AddEntry( bgdx, "Interpolated Background (scaled)", "l");
  dxlegend->AddEntry( tfitdx, "Total fit (background + signal)", "l");
  dxlegend->AddEntry( (TObject*)0, "", "");
  dxlegend->AddEntry( (TObject*)0, Form("HDE via dx: %0.3f%%",effdx), "");
  dxlegend->Draw();

  c2->Write();
  fout->Write();

  // Send time efficiency report to console
  st->Stop();
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
