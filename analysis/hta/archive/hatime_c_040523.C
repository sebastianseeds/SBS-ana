//sseeds 04.03.23 - Script to extract improved ADC time via application of ToF corrections and waveform corrections

#include <vector>
#include <iostream>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include <TF1Convolution.h>
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

//Passing kine==-1 will run all kinematics
void hatime_c( Int_t kine=4 )
{ //main
  
  Int_t epm = 2; //Set e' momentum reconstruction to simplest method
  Int_t pass = 1; //fix for pass 0/1 for now
  Double_t adcbinw = 4.; //adc sample bin width (4*40=160ns over whole waveform)

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/shatime_c.json");
  std::string rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");

  vector<Int_t> lh2r;
  jmgr->GetVectorFromSubKey<Int_t>( "lh2", "sbs4_runs", lh2r );
  jmgr->GetVectorFromSubKey<Int_t>( "lh2", "sbs8_runs", lh2r );
  jmgr->GetVectorFromSubKey<Int_t>( "lh2", "sbs9_runs", lh2r );

  vector<Int_t> ld2r;
  jmgr->GetVectorFromSubKey<Int_t>( "ld2", "sbs11_runs", ld2r );

  //set up default parameters for all analysis
  std::string runsheet_dir = "/w/halla-scshelf2102/sbs/seeds/ana/data"; //unique to my environment for now
  Int_t nruns = -1; //Always analyze all available runs
  Int_t verb = 0; //Don't print diagnostic info by default

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<crun> crunh; 
  vector<crun> crund; 
  util::ReadRunList(runsheet_dir,nruns,kine,"LH2",pass,verb,crunh); //modifies nruns to be very large when -1
  util::ReadRunList(runsheet_dir,nruns,kine,"LD2",pass,verb,crund); //modifies nruns to be very large when -1

  //std::string rootfile_dir = jmgr->GetValueFromKey_str(Form("rootfile_dir_%s_sbs%d",target.c_str(),kine));
  
  //set up output files
  TFile *fout = new TFile( Form("outfiles/ohatime_sbs%d.root",kine), "RECREATE" );

  //set up waveform histograms
  TH1D *landhist[econst::hcalrow][econst::hcalcol];
  TH1D *sghist[econst::hcalrow][econst::hcalcol];
  TH1D *lghist[econst::hcalrow][econst::hcalcol];
  for(Int_t r = 0; r < econst::hcalrow; r++) {
    for(Int_t c = 0; c < econst::hcalcol; c++) {
      landhist[r][c] = util::hhsamps(r,c,econst::maxsamp);
      sghist[r][c] = (TH1D*)landhist[r][c]->Clone(Form("sghist %d-%d",r,c));
      lghist[r][c] = (TH1D*)landhist[r][c]->Clone(Form("lghist %d-%d",r,c));
    }
  }

  //set up general histograms
  TH2D *hatimeID = new TH2D("hatimeID","HCal tree ADC time vs ID; channel;ns",0,econst::hcalchan,econst::hcalchan,econst::maxsamp*4.,econst::minsamp,econst::maxsamp*4.);
  TH2D *hatimeclID = new TH2D("hatimeclID","HCal waveform ADC time (landau) vs ID; channel;ns",0,econst::hcalchan,econst::hcalchan,econst::maxsamp*4.,econst::minsamp,econst::maxsamp*4.);
  TH2D *hatimecsgID = new TH2D("hatimecsgID","HCal waveform ADC time (skew gaus) vs ID; channel;ns",0,econst::hcalchan,econst::hcalchan,econst::maxsamp*4.,econst::minsamp,econst::maxsamp*4.);
  TH2D *hatimeclgID = new TH2D("hatimeclgID","HCal waveform ADC time (gaus landau conv) vs ID; channel;ns",0,econst::hcalchan,econst::hcalchan,econst::maxsamp*4.,econst::minsamp,econst::maxsamp*4.);
  TH2D *hatimecbID = new TH2D("hatimecbID","HCal waveform ADC time (min ch2) vs ID; channel;ns",0,econst::hcalchan,econst::hcalchan,econst::maxsamp*4.,econst::minsamp,econst::maxsamp*4.);

  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;
  
  // create output tree
  TTree *P = new TTree("P","Analysis Tree"); 
  
  // output tree vars
  Double_t pblkid_out;
  Double_t latime_out;
  Double_t sgatime_out;
  Double_t batime_out;
  Double_t atime_out;
  
  Double_t dx_out;
  Double_t dy_out;
  Double_t W2_out;
  Double_t Q2_out;
  Double_t hcale_out;
  Double_t hcalatime_out;
  Double_t hodotmean_out;
  Double_t thetapq_pout;
  Double_t thetapq_nout;
  Int_t mag_out;
  Int_t run_out;
  Int_t tar_out;
  Int_t failedglobal_out;
  Int_t failedaccmatch_out;
  Int_t failedcoin_out;

  // set output tree branches
  P->Branch( "pblkid", &pblkid_out, "dx/D" );
  P->Branch( "latime", &latime_out, "dx/D" );
  P->Branch( "sgatime", &sgatime_out, "dx/D" );
  P->Branch( "batime", &batime_out, "dx/D" );
  P->Branch( "atime", &atime_out, "dx/D" );

  P->Branch( "dx", &dx_out, "dx/D" );
  P->Branch( "dy", &dy_out, "dy/D" );
  P->Branch( "W2", &W2_out, "W2/D" );
  P->Branch( "Q2", &Q2_out, "Q2/D" );
  P->Branch( "hcale", &hcale_out, "hcale/D" );
  P->Branch( "hcalatime", &hcalatime_out, "hcalatime/D" );
  P->Branch( "hodotmean", &hodotmean_out, "hodotmean/D" );
  P->Branch( "thetapq_p", &thetapq_pout, "thetapq_pout/D" );
  P->Branch( "thetapq_n", &thetapq_nout, "thetapq_nout/D" );
  P->Branch( "mag", &mag_out, "mag_out/I" );
  P->Branch( "run", &run_out, "run_out/I" );
  P->Branch( "tar", &tar_out, "tar_out/I" );
  P->Branch( "failedglobal", &failedglobal_out, "failedglobal/I" );
  P->Branch( "failedaccmatch", &failedaccmatch_out, "failedaccmatch/I" );
  P->Branch( "failedcoin", &failedcoin_out, "failedcoin/I" );

  // setup reporting indices
  Int_t curmag = -1;
  std::string curtar = "";
  //bool testover = false;

  for ( Int_t t=0; t<2; t++ ){
    //for ( Int_t t=0; t<2; t++ ){  //replace for testing
    // t==0, lh2; t==1, ld2

    for (Int_t irun=0; irun<nruns; irun++) {
      //if( testover==true ) continue;
      // accessing run info
      Int_t runnum;
      std::string targ;
      if( t==0 ){
	targ = "LH2";
	runnum = crunh[irun].runnum;
	bool skip = true;
	for( Int_t el=0; el<lh2r.size(); el++ ){     
	  if( runnum==lh2r[el] ) skip=false;
	}
	if( skip==true ) continue;
      }
      if( t==1 ){
	targ = "LD2";
	runnum = crund[irun].runnum;
	bool skip = true;
	for( Int_t el=0; el<ld2r.size(); el++ ){     
	  if( runnum==ld2r[el] ) skip=false;
	}
	if( skip==true ) continue;
      }

      //testover=true;
      std::cout << "Analyzing run " << runnum << ".." << std::endl;
      

      std::string rfname = rootfile_dir + Form("/*%d*",crunh[irun].runnum);

      Int_t mag = crunh[irun].sbsmag / 21; //convert to percent
      Double_t ebeam = crunh[irun].ebeam; //get beam energy per run
      std::string tar = crunh[irun].target; //get target for current run
      Int_t conf = crunh[irun].sbsconf; //get config of current run

      if( conf!=kine && kine!=-1) continue;

      if( targ.compare( tar )!=0 ){
	std::cout << "Error: target from grl and json mismatch. Check json file and resubmit." << std::endl;
	return;
      }

      //set up configuration and tune objects to load analysis parameters
      //SBSconfig *config = new SBSconfig(kine,mag);
      SBSconfig config(kine,mag);

      //Obtain configuration pars from config file
      Double_t hcaltheta = config.GetHCALtheta_rad();
      Double_t hcaldist = config.GetHCALdist();
      Double_t sbsdist = config.GetSBSdist();
    
      //SBStune *tune = new SBStune(kine,mag);
      SBStune tune(kine,mag);
    
      if( tar.compare(curtar)!=0 || mag!=curmag ){
	std::cout << "Settings change.." << endl;
	cout << config;
	cout << tune;
	curtar = tar;
	curmag = mag;
      }

      //Obtain cuts from tune class
      std:string gcut   = tune.Getglobcut_hexp();
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

      //std::string rfname = rootfile_dir + Form("/*%d*",crunh[irun].runnum);
      C = new TChain("T");
      C->Add(rfname.c_str());

      // setting up ROOT tree branch addresses
      C->SetBranchStatus("*",0);    

      // HCal general
      Double_t hcalid, hcale, hcalx, hcaly, hcalr, hcalc, hcala, hcalamp, hcaltdc, hcalatime;
      std::vector<std::string> hcalvar = {"idblk","e","x","y","rowblk","colblk","a_p","a_amp_p","tdctimeblk","atimeblk"};
      std::vector<void*> hcalvarlink = {&hcalid,&hcale,&hcalx,&hcaly,&hcalr,&hcalc,&hcala,&hcalamp,&hcaltdc,&hcalatime};
      rvars::setbranch(C, "sbs.hcal", hcalvar, hcalvarlink);

      // Double_t nsamps[econst::maxchan], samps[econst::maxsamps], sampsrow[econst::maxchan], sampscol[econst::maxchan], a_p[econst::maxchan], a_amp_p[econst::maxchan], sampsidx[econst::maxchan];
      // Int_t Nssave;
      // C->SetBranchStatus("sbs.hcal.a_p",1);
      // C->SetBranchAddress("sbs.hcal.a_p",a_p);
      // C->SetBranchStatus("sbs.hcal.a_amp_p",1);
      // C->SetBranchAddress("sbs.hcal.a_amp_p",a_amp_p);
      // C->SetBranchStatus("sbs.hcal.nsamps",1);
      // C->SetBranchAddress("sbs.hcal.nsamps",nsamps);
      // C->SetBranchStatus("sbs.hcal.samps",1);
      // C->SetBranchAddress("sbs.hcal.samps",samps);
      // C->SetBranchStatus("sbs.hcal.samps_idx",1);
      // C->SetBranchAddress("sbs.hcal.samps_idx",sampsidx);
      // C->SetBranchStatus("sbs.hcal.adcrow",1);
      // C->SetBranchAddress("sbs.hcal.adcrow",sampsrow);
      // C->SetBranchStatus("sbs.hcal.adccol",1);
      // C->SetBranchAddress("sbs.hcal.adccol",sampscol);
      // C->SetBranchStatus("Ndata.sbs.hcal.adcrow",1);
      // C->SetBranchAddress("Ndata.sbs.hcal.adcrow",&Nssave);

      //HCal samples
      Double_t nsamps[econst::maxchan], samps[econst::maxsamps], sampsrow[econst::maxchan], sampscol[econst::maxchan], a_p[econst::maxchan], a_amp_p[econst::maxchan], sampsidx[econst::maxchan];
      Int_t Nssave;
      std::vector<std::string> hcalsvar = {"nsamps","samps","adcrow","adccol","samps_idx","a_p","a_amp_p","adcrow"};
      std::vector<void*> hcalsvarlink = {&nsamps,&samps,&sampsrow,&sampscol,&sampsidx,&a_p,&a_amp_p,&Nssave};
      rvars::setbranch(C, "sbs.hcal", hcalsvar, hcalsvarlink, 7);

      // HCal cluster branches
      Double_t hcalcid[econst::maxclus], hcalce[econst::maxclus], hcalcx[econst::maxclus], hcalcy[econst::maxclus], hcalctdctime[econst::maxclus], hcalcatime[econst::maxclus];
      Int_t Nhcalcid;
      std::vector<std::string> hcalcvar = {"id","e","x","y","tdctime","atime","id"};
      std::vector<void*> hcalcvarlink = {&hcalcid,&hcalce,&hcalcx,&hcalcy,&hcalctdctime,&hcalcatime,&Nhcalcid};
      rvars::setbranch(C, "sbs.hcal.clus", hcalcvar, hcalcvarlink, 6);

      // HCal primary clust block branches
      Double_t hcalcbid[econst::maxclus], hcalcbe[econst::maxclus], hcalcbx[econst::maxclus], hcalcby[econst::maxclus], hcalcbrow[econst::maxclus], hcalcbcol[econst::maxclus], hcalcbtdctime[econst::maxclus], hcalcbatime[econst::maxclus];
      Int_t Nhcalcbid;
      std::vector<std::string> hcalcbvar = {"id","e","row","col","x","y","tdctime","atime","id"};
      std::vector<void*> hcalcbvarlink = {&hcalcbid,&hcalcbe,&hcalcbrow,&hcalcbcol,&hcalcbx,&hcalcby,&hcalcbtdctime,&hcalcbatime,&Nhcalcbid};
      rvars::setbranch(C, "sbs.hcal.clus_blk", hcalcbvar, hcalcbvarlink, 8);

      // hodoscope cluster mean time
      Int_t Nhodotmean; 
      Double_t hodotmean[econst::maxclus];
      std::vector<std::string> hodovar = {"clus.tmean","clus.tmean"};
      std::vector<void*> hodovarlink = {&hodotmean,&Nhodotmean};
      rvars::setbranch(C, "bb.hodotdc", hodovar, hodovarlink, 1);  
    
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
      Int_t pltcntr = 0;

      while (C->GetEntry(nevent++)) {

	cout << nevent << "/" << nevents << " \r";
	cout.flush();

	//if( nevent>1000 ) break; //testing

	///////
	//Single-loop globalcut method. Save pass/fail for output tree.
	currenttreenum = C->GetTreeNumber();
	if( nevent == 1 || currenttreenum != treenum ){
	  treenum = currenttreenum; 
	  GlobalCut->UpdateFormulaLeaves();
	  cout << "Updating formula leaves and switching segment at event: " << nevent << endl;
	}
	bool failedglobal = GlobalCut->EvalInstance(0) == 0;

	//if( failedglobal ) continue;
	Int_t fglobalint = 0;
	if( failedglobal ) fglobalint = 1;
	///////
      
	///////
	//HCal coin cut. Save pass/fail for output tree.
	bool failedcoin = abs( hcalcbatime[0] - atime0 ) > 3*atimesig;

	//if( failedcoin ) continue;
	Int_t fcoinint = 0;
	if( failedcoin ) fcoinint = 1;
	///////

	///////
	//Physics calculations
	//correct beam energy with vertex information
	Double_t ebeam_c = vars::ebeam_c( ebeam, vz[0], targ );
	TVector3 vertex( 0., 0., vz[0] );

	//set up four-momenta with some empty for various calculation methods
	TLorentzVector pbeam( 0., 0., ebeam_c, ebeam_c ); //beam momentum
	TLorentzVector pe( px[0], py[0], pz[0], p[0] ); //e' momentum
	TLorentzVector ptarg; vars::setPN(nucleon,ptarg); //target momentum
	TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
	TLorentzVector pN; //N' momentum
	TVector3 pNhat; //Unit N' 3-vector
      
	//simple calculations for e' and N'
	Double_t etheta = vars::etheta(pe); 
	Double_t ephi = vars::ephi(pe);
	Double_t pcent = vars::pcentral(ebeam,etheta,nucleon); //e' p reconstructed by angles
	Double_t phNexp = ephi + physconst::pi;
	Double_t Q2, W2, nu, thNexp, pNexp;
	Double_t ebeam_o = vars::ebeam_o( ebeam_c, etheta, targ ); //Second energy correction accounting for energy loss leaving target

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
	  pNhat = vars::qVect_unit( thNexp, phNexp );
	  pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	}else if( epm==3 ){
	  //v3
	  nu = pbeam.E() - pcent;
	  pNexp = vars::pN_expect( nu, nucleon );
	  thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	  pNhat = vars::qVect_unit( thNexp, phNexp );
	  pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	  Q2 = vars::Q2( pbeam.E(), pe.E(), etheta );
	  W2 = vars::W2( pbeam.E(), pe.E(), W2, nucleon );
	}else if( epm==4 ){
	  //v4
	  nu = pbeam.E() - pe.E();
	  pNexp = vars::pN_expect( nu, nucleon );
	  thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	  pNhat = vars::qVect_unit( thNexp, phNexp );
	  pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	  Q2 = vars::Q2( pbeam.E(), pe.E(), etheta );
	  W2 = vars::W2( pbeam.E(), pe.E(), W2, nucleon );
	}else{
	  Q2 = ekineQ2;
	  W2 = ekineW2;
	  nu = ekinenu;
	  pNexp = vars::pN_expect( nu, nucleon );
	  thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	  pNhat = vars::qVect_unit( thNexp, phNexp );
	  pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	  std::cout << "Warning: epm version incorrect. Defaulting to version 2." << endl;
	}

	//Calculate h-arm quantities
	vector<Double_t> xyhcalexp; vars::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
	TVector3 hcalpos = hcalorigin - hcalcbx[0]*hcalaxes[0] + hcalcby[0]*hcalaxes[1]; //from primary blk
	Double_t dx = hcalcbx[0] - xyhcalexp[0];
	Double_t dy = hcalcby[0] - xyhcalexp[1];
	TVector3 neutdir = (hcalpos - vertex).Unit();
	Double_t protdeflect = tan( 0.3 * BdL / q.Mag() ) * (hcaldist - (sbsdist + econst::sbsdipolegap/2.0) );
	TVector3 protdir = ( hcalpos + protdeflect * hcalaxes[0] - vertex ).Unit();

	//Fill diagnostic histos

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

	///////
	//extract atime with waveforms from primary cluster primary block for now
	Int_t pblkid = hcalcbid[0];
	Int_t r,c,idx,n,sub;
	Double_t adc[econst::hcalrow][econst::hcalcol];
	Double_t amp[econst::hcalrow][econst::hcalcol];
	Float_t peak[econst::hcalrow][econst::hcalcol];
	//Double_t adc_p[econst::hcalrow][econst::hcalcol];

	for( r = 0; r < econst::hcalrow; r++ ) {
	  for( c = 0; c < econst::hcalcol; c++ ) {
	    landhist[r][c]->Reset( "ICES M" );
	    sghist[r][c]->Reset( "ICES M" );
	    lghist[r][c]->Reset( "ICES M" );
	    peak[r][c] = 0.0;
	    adc[r][c] = 0.0;
	    amp[r][c] = 0.0;
	  }
	}

	//all identical -> cout << "hcalid: " << hcalid << " hcalcid[0]: " << hcalcid[0] << " hcalcbid[0]:" << hcalcbid[0] << endl;
	
	//Primary cluster element
	Int_t pblkrow = (int)hcalcbid[0]/econst::hcalcol;
	Int_t pblkcol = (int)hcalcbid[0]%econst::hcalcol;
	//cout << "EVENT: " << nevent << " pblkrow: " << pblkrow << " pblkcol: " << pblkcol << endl;

	//atime fit parameters
	Double_t sgmpv = -1000.;
	Double_t lmpv = -1000.;
	Double_t brising_edge = -1000.;
	Double_t lrising_edge = -1000.;
	Double_t lrising_edge_v2 = -1000.;
	Double_t sgrising_edge = -1000.;
	Double_t sgrising_edge_v2 = -1000.;
	Double_t lgrising_edge = -1000.;

	bool pblkfill = false;

	for(Int_t m = 0; m < Nssave; m++) {

	  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
	  //gSystem->RedirectOutput("/dev/null");

	  //cout << "entering samps.." << endl;
	  r = sampsrow[m];
	  c = sampscol[m];
	  Int_t chan = econst::hcalcol*r+c;

	  //Check indices and continue if problem or end of samples exists
	  if(r < 0 || c < 0) {
	    std::cerr << "Warning: row/col negative. Check indices." << std::endl;
	    continue;
	  }
	  if(r >= econst::hcalrow || c >= econst::hcalcol){
	    cout << "Error: r/c out of bounds" << endl;
	    break; 
	  }
	  idx = sampsidx[m];

	  //cout << "chan: " << chan << " idx: " << idx << endl;

	  //if( chan!=hcalcbid[0] ) continue; //continue if samples do not correspond to primary block

	  n = nsamps[m];
	  
	  adc[r][c] = a_p[m];
	  amp[r][c] = a_amp_p[m];

	  bool saturated = amp[r][c]>3500;
	  bool mansat = false;
	  bool negped = adc[r][c]<-5;
	  bool skip = true;
	  for( int b = 0; b < Nhcalcbid; b++ ){
	    if( chan==(hcalcbid[0]-1) ){ 
	      skip=false;
	      //cout << chan << " " << hcalcbid[0]-1 << " " << hcalcbid[b]-1 << endl;
	    }
	  }

	  if( !saturated && !negped && skip==true ) continue;
      
	  //cout << "IN mLOOP event: " << nevent << " pblkrow: " << pblkrow << " pblkcol: " << pblkcol << endl;
	  //cout << "IN mLOOP event: " << nevent << " post cut r: " << r << " c: " << c << endl;

	  bool displayed = false;
	  for(Int_t s = econst::minsamp; s < econst::maxsamp && s < n; s++) {
	    displayed = true;
	    landhist[r][c]->SetBinContent(s+1-econst::minsamp,samps[idx+s]);
	    sghist[r][c]->SetBinContent(s+1-econst::minsamp,samps[idx+s]);
	    lghist[r][c]->SetBinContent(s+1-econst::minsamp,samps[idx+s]);
	    if(peak[r][c]<samps[idx+s])
	      peak[r][c]=samps[idx+s];
	    if(peak[r][c]>2) {
	      mansat = true;
	    }
	  }
	  //landhist[r][c]->SetLineColor(kBlue+1);
	  // if( saturated==true || negped==true ){
	  //   cout << "Bad block found at r:" << r+1 << " c:" << c+1;
	  //   if( saturated==true ) cout << " -saturated";
	  //   if( negped==true ) cout << " -negative pedestal";
	  //   cout << endl;
	  //   //subCanv[2]->cd(errIdx);
	  //   //subCanv[2]->SetGrid();
	  //   landhist[r][c]->SetTitle(TString::Format("%d-%d (TDC=%g,TMu=%g)",r+1,c+1,tdc[r][c],tdcmult[r][c]));
	  //   landhist[r][c]->SetLineColor(kRed+1);
	  //   landhist[r][c]->Draw();
	  //   gPad->SetFillColor(18);
	  //   gPad->Update();
	  //   errIdx++;
	  // }
	  if(!displayed) {
	    std::cerr << "Skipping empty module: " << m << std::endl;
	    for(Int_t s = 0;  s < econst::maxsamp; s++) {
	      landhist[r][c]->SetBinContent(s+1,-404);
	      sghist[r][c]->SetBinContent(s+1,-404);
	      lghist[r][c]->SetBinContent(s+1,-404);
	    }
	  }
	  
	  /////////////////////////////
	  //if( pblkrow!=r || pblkcol!=c ) continue;
	  /////////////////////////////

	  TF1 *lfit;
	  TF1 *sgfit = new TF1(Form("sg r:%d c:%d",r,c),fits::g_sgfit, econst::minsamp, econst::maxsamp-1, 4);
	  //TF1 *f = new TF1("f","2.*gaus(x,[0],[1],[2])*ROOT::Math::normal_cdf([3]*x,1,0)",econst::minsamp, econst::maxsamp);
	  //TF1 *pfit = new TF1(Form("p r:%d c:%d",r,c),fits::g_p3fit, econst::minsamp, econst::maxsamp, 4);
	  //TH1D *hsclone;
	  //TH1D *fclone;
	  //TH1D *pclone;

	  //hsclone = (TH1D*)landhist[r][c]->Clone(Form("hsclone %d-%d",r,c));
	  //fclone = (TH1D*)landhist[r][c]->Clone(Form("fclone %d-%d",r,c));
	  //pclone = (TH1D*)landhist[r][c]->Clone(Form("pclone %d-%d",r,c));

	  
	  //landau fit
	  landhist[r][c]->Fit("landau", "Q", "", econst::minsamp, econst::maxsamp);
	  lfit = landhist[r][c]->GetFunction("landau");
	  lmpv = lfit->GetMaximumX();
	  Double_t lch2 = lfit->GetChisquare();
	  Double_t lmean = lfit->GetParameter(1);
	  Double_t lsig = lfit->GetParameter(2);
	  lrising_edge = lmean-3*lsig;
	  lrising_edge_v2 = lfit->GetX(0.00001,econst::minsamp+1,econst::maxsamp-1);
	  //landhist[r][c]->SetTitle(Form("landau ev:%ld r:%d c:%d adc:%0.2f amp:%0.2f mpv:%0.2f mean:%0.2f sig:%0.2f chi2:%0.2f",nevent,r,c,adc[r][c],amp[r][c],lmpv,lmean,lsig,lch2));
	  
	  

	  //landhist[r][c]->Write();

	  //sgfit
	  Double_t sgamean = sghist[r][c]->GetMean();
	  Double_t sgasig = (econst::maxsamp-1)*0.15;
	  //sgfit->SetParameters(hsclone->GetMaximum(),hsclone->GetMean(),2.,1.);
	  sgfit->SetParameters(sghist[r][c]->GetMaximum(),sgamean,2.,1.);
	  sgfit->SetParLimits(0,0.005,2.0); //bad limits
	  sgfit->SetParLimits(1,econst::minsamp,econst::maxsamp-1); //bad limits
	  // //sgfit->SetParLimits(1,sgamean-sgasig,sgamean+sgasig);
	  // //sgfit->SetParLimits(2,0.,(econst::maxsamp-1)*0.2);
	  sgfit->SetParLimits(2,4.,sgasig); //good limits
	  sgfit->SetParLimits(3,5.,13.); //bad limits
	  sghist[r][c]->Fit(sgfit,"RQ","",econst::minsamp, econst::maxsamp);
	  sgmpv = sgfit->GetMaximumX();
	  Double_t sgsig = sgfit->GetParameter(2);
	  Double_t sskew = sgfit->GetParameter(3);
	  sgrising_edge = sgmpv+sgsig*sqrt(2)*TMath::ErfInverse(-sskew); //not very effective
	  //TF1 *threshold = new TF1("threshold",fits::g_thresh,econst::minsamp,econst::maxsamp,1);
	  //threshold->SetParameter(0,0.0005);
	  sgrising_edge_v2 = sgfit->GetX(0.00001,econst::minsamp+1,econst::maxsamp-1);
	  //cout << endl << sgrising_edge << endl;
	  Double_t sgch2 = sgfit->GetChisquare();
	  //Double_t sgmean = sgfit->GetParameter(1);
	  //Double_t sgsig = sgfit->GetParameter(2);
	  //Double_t sskew = sgfit->GetParameter(3);
	  //sghist[r][c]->SetTitle(Form("sgaus ev:%ld r:%d c:%d adc:%0.2f amp:%0.2f mpv:%0.2f mean:%0.2f sig:%0.2f skew:%0.2f chi2:%0.2f",nevent,r,c,adc[r][c],amp[r][c],sgmpv,sgmean,sgsig,sskew,sgch2));

	  //cout << "skewg_mpv: " << sgmpv << " skewg_ch2: " << sgch2 << " land_mpv: " << lmpv << " land_ch2: " << lch2 << endl;

	  //cout << "event: " << nevent << " chan: " << chan << " lmpv*adcbinw: " << lmpv*adcbinw << " sgmpv*adcbinw: " << sgmpv*adcbinw << " hcalcbatime[0]: " << hcalcbatime[0] << endl;
	  
	  //cout << endl << lrising_edge_v2 << endl;

	  hatimeclID->Fill(chan,lrising_edge_v2*adcbinw);
	  hatimecsgID->Fill(chan,sgrising_edge_v2*adcbinw);
	  //hatimeclgID->Fill(chan,sgmpv*adcbinw);
	  hatimeID->Fill(chan,hcalcbatime[0]);

	  //Get best rising edge
	  bool glre = lrising_edge_v2>econst::minsamp && lrising_edge_v2<econst::maxsamp; //good landau fit re
	  bool gsgre = sgrising_edge_v2>econst::minsamp && sgrising_edge_v2<econst::maxsamp; //good skew gaus fit re
	  if( lch2 < sgch2 && glre ){ //if landau chi2 better than skew gaus chi2 and good landau fit re
	    brising_edge = lrising_edge_v2*adcbinw;
	  }else if( gsgre ){ //if skew gaus chi2 better than landau chi2 and good skew gaus fit re
	    brising_edge = sgrising_edge_v2*adcbinw;
	  }else{ //if no good fit re, use tree atime
	    brising_edge = hcalcbatime[0];
	  }
	  hatimecbID->Fill(chan,brising_edge);

	  //gaussian landau convolution fit. This takes very long..
	  // TF1Convolution *f_conv = new TF1Convolution("landau","gaus",econst::minsamp-1,econst::maxsamp,true);
	  // f_conv->SetRange(econst::minsamp-1,econst::maxsamp);
	  // f_conv->SetNofPointsFFT(1000);
	  // TF1 *f = new TF1("f",*f_conv,econst::minsamp,econst::maxsamp-1,f_conv->GetNpar());
	  // //5 fit parameters, 
	  // //  p[0]:landau mean, 
	  // //  p[1]:landau sigma, 
	  // //  p[2]:gaus amp, 
	  // //  p[3]:gaus mean,
	  // //  p[4]:gaus sigma
	  // f->SetParameters( sgamean, sgasig, lghist[r][c]->GetMaximum(), sgamean, sgasig );
	  
	  // lghist[r][c]->Fit("f","RQ","",econst::minsamp, econst::maxsamp-1);
	  // Double_t lgp1 = f->GetParameter(0);
	  // Double_t lgp2 = f->GetParameter(1);
	  // Double_t lgp3 = f->GetParameter(2);
	  // Double_t lgp4 = f->GetParameter(3);
	  // Double_t lgp5 = f->GetParameter(4);
	  // lgrising_edge = f->GetX(0.002,econst::minsamp+1,econst::maxsamp-1);
	  // Double_t lgch2 = f->GetChisquare();

	  //diagnostic plots
	  // if( pltcntr<20 && nevent%5==0 && pblkfill==false){

	  //   lghist[r][c]->SetTitle(Form("lgconv ev:%ld r:%d c:%d at0:%f re:%0.2f p1:%0.2f p2:%0.2f p3:%0.2f p4:%0.2f p5:%0.2f chi2:%0.2f",nevent,r,c,hcalatime/4.,lgrising_edge,lgp1,lgp2,lgp3,lgp4,lgp5,lgch2));
	  //   lghist[r][c]->SetName(Form("lgev%ld",nevent));
	  //   lghist[r][c]->Write();

	  //   //if( lch2 < sgch2 ){
	  //     Double_t lmean = lfit->GetParameter(1);
	  //     landhist[r][c]->SetTitle(Form("landau ev:%ld r:%d c:%d at0:%f re:%0.2f rev2:%0.2f adc:%0.2f amp:%0.2f mpv:%0.2f mean:%0.2f sig:%0.2f chi2:%0.2f",nevent,r,c,hcalatime/4.,lrising_edge,lrising_edge_v2,adc[r][c],amp[r][c],lmpv,lmean,lsig,lch2));
	  //     landhist[r][c]->SetName(Form("lev%ld",nevent));
	  //     landhist[r][c]->Write();
	  //     //}else{
	  //     Double_t sgmean = sgfit->GetParameter(1);
	  //     sghist[r][c]->SetTitle(Form("sgaus ev:%ld r:%d c:%d at0:%f re:%0.2f rev2:%0.2f adc:%0.2f amp:%0.2f mpv:%0.2f mean:%0.2f sig:%0.2f skew:%0.2f chi2:%0.2f",nevent,r,c,hcalatime/4.,sgrising_edge,sgrising_edge_v2,adc[r][c],amp[r][c],sgmpv,sgmean,sgsig,sskew,sgch2));
	  //     sghist[r][c]->SetName(Form("sgev%ld",nevent));
	  //     sghist[r][c]->Write();
	  //     //}
	  //   pltcntr++;
	  //   pblkfill=true;
	  // }

	  // if( lgrising_edge==-1000. || lrising_edge==-1000. || sgrising_edge==-1000. ) cout << endl << "lgrising_edge:" << lgrising_edge << " lrising_edge:" << lrising_edge << " sgrising_edge:" << sgrising_edge << endl;

	  //f - Results in singular fits and emphasis on gaussians
	  // f->SetParameters(hsclone->GetMaximum(),hsclone->GetMean(),2.,1.);
	  // f->SetParLimits(0,0.005,2.0);
	  // f->SetParLimits(1,econst::minsamp,econst::maxsamp);
	  // f->SetParLimits(2,0.,(econst::maxsamp-1)*0.25);
	  // f->SetParLimits(3,0.05,15.);
	  // fclone->Fit(f,"RQ","",econst::minsamp, econst::maxsamp-1);
	  // Double_t fmean = f->GetParameter(1);
	  // Double_t fskew = f->GetParameter(3);
	  // Double_t fmpv = f->GetMaximumX();
	  // fclone->SetTitle(Form("sgaus-f r:%d c:%d adc:%f amp:%f mpv:%0.2f mean:%0.2f skew:%f",r,c,adc[r][c],amp[r][c],fmpv,fmean,fskew));

	  //hsclone->Write();

	  // pclone->Fit(pfit,"RQ","",econst::minsamp, econst::maxsamp);
	  // //Double_t pmean = pfit->GetParameter(1);
	  // Double_t pmpv = pfit->GetMaximumX();
	  // Double_t pch2 = pfit->GetChisquare();
	  // pclone->SetTitle(Form("pfit r:%d c:%d adc:%f amp:%f mpv:%0.2f chi2:%0.2f",r,c,adc[r][c],amp[r][c],pmpv,pch2));

	  // pclone->Write();
	  // if( nevent>6000 && nevent<6100 && amp[r][c]>80. ){
	    
	  //   landhist[r][c]->Write();
	  //   sghist[r][c]->Write();
	  //   //fclone->Write();

	  //   // if( lch2 < sgch2 ){
	  //   //   landhist[r][c]->Write();
	  //   // }else{
	  //   //   hsclone->Write();
	  //   // }

	  // }
	}



	//Fill physics output tree  
	if( pblkfill==false ){

	  pblkid_out = (double)hcalcbid[0];
	  latime_out = lrising_edge_v2*adcbinw;
	  sgatime_out = sgrising_edge_v2*adcbinw;
	  batime_out = brising_edge*adcbinw;
	  atime_out = hcalcbatime[0];
   
	  dx_out = dx;
	  dy_out = dy;
	  W2_out = W2;
	  Q2_out = Q2;
	  hcale_out = hcale;
	  hcalatime_out = hcalcbatime[0];
	  hodotmean_out = hodotmean[0];
	  thetapq_pout = acos( protdir.Dot( pNhat ) );
	  thetapq_nout = acos( neutdir.Dot( pNhat ) );
	  mag_out = mag;
	  run_out = runnum;
	  tar_out = t;
	  failedglobal_out = fglobalint;
	  failedaccmatch_out = faccmatchint;
	  failedcoin_out = fcoinint;

	  P->Fill();
	  pblkfill=true;
	}

      }

      // getting ready for the next run
      C->Reset();

    }//endloop over runs

  }//endloop over targets

  fout->Write();

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
