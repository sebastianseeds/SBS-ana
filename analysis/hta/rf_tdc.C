//sseeds 5.22.23 - script to compare rftime and tdc time from hcal reference tdc and tdc respectively

//NOTE: currently configured to function only on HCal expert replays

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

//Passing kine==-1 will run all kinematics, pass is replay pass
void rf_tdc( Int_t kine=4, Int_t pass=1 )
{ //main

  Double_t adcbinw = econst::hcaladc_binw; //adc sample bin width (4*40=160ns over whole waveform)

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/shatime_c.json");
  std::string rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");

  vector<Int_t> lh2r;
  jmgr->GetVectorFromSubKey<Int_t>( "lh2", "sbs4_runs", lh2r );
  jmgr->GetVectorFromSubKey<Int_t>( "lh2", "sbs9_runs", lh2r );

  vector<Int_t> ld2r;
  jmgr->GetVectorFromSubKey<Int_t>( "ld2", "sbs8_runs", ld2r );
  jmgr->GetVectorFromSubKey<Int_t>( "ld2", "sbs11_runs", ld2r );

  //set up default parameters for all analysis
  std::string runsheet_dir = "/w/halla-scshelf2102/sbs/seeds/ana/data"; //unique to my environment for now
  Int_t nruns = -1; //Always analyze all available runs for this application. If not, need to pass nruns for lh2 and ld2 separately.
  Int_t verb = 0; //Don't print diagnostic info by default

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<crun> crunh; 
  vector<crun> crund; 
  util::ReadRunList(runsheet_dir,nruns,kine,"LH2",pass,verb,crunh); //modifies nruns to be very large when -1
  util::ReadRunList(runsheet_dir,nruns,kine,"LD2",pass,verb,crund); //modifies nruns to be very large when -1

  //set up output files
  TFile *fout = new TFile( Form("outfiles/rftdc_sbs%d.root",kine), "RECREATE" );

  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;
  
  // create output tree
  TTree *P = new TTree("P","Analysis Tree"); 
  
  // output tree vars
  Double_t pblkid_out; //hcal primary cluster, primary block id
  Double_t atime_out; //hcal primary cluster adc time, tree

  Double_t bcatime_dxdy_out; //hcal best cluster adc time, tree
  Double_t bcatime_atime_out; //hcal best cluster adc time, tree
  Double_t bcatime_thpq_out; //hcal best cluster adc time, tree
  Double_t bcpblkid_dxdy_out; //hcal best cluster, primary block id
  Double_t bcpblkid_atime_out; //hcal best cluster, primary block id
  Double_t bcpblkid_thpq_out; //hcal best cluster, primary block id

  Double_t dx_out; //hcal primary cluster dx
  Double_t dy_out; //hcal primary cluster dy 
  Double_t W2_out;
  Double_t Q2_out;
  Double_t hcale_out; //hcal primary cluster energy
  Double_t pse_out; //bbcal preshower primary cluster energy
  Double_t she_out; //bbcal shower primary cluster energy
  Double_t ep_out; //track reconstructed e' momentum
  Double_t hcalatime_out; //hcal primary cluster adc time, tree
  Double_t hodotmean_out; //hodoscope primary cluster mean tdc time
  Double_t thetapq_pout; //proton
  Double_t thetapq_nout; //neutron
  Int_t pid_out; //-1:neither,0:ambiguous,1:proton,2:neutron
  Int_t mag_out; //sbs magnetic field strength (percent)
  Int_t run_out; //run number
  Int_t tar_out; //target, LH2 or LD2
  Int_t failedglobal_out; //failed global cut (=1)
  Int_t failedaccmatch_out; //failed acceptance matching, hcal active area (=1)
  Int_t failedcoin_out; //failed coincidence time, hcal adc time (=1)
  Int_t failedw2_out; //failed w2 cut on elastic peak

  //cluster tree vars
  Double_t cpblkid_out[econst::maxclus];
  Double_t chatime_out[econst::maxclus];
  Double_t cthetapq_p_out[econst::maxclus];
  Double_t cthetapq_n_out[econst::maxclus];
  Double_t cdx_out[econst::maxclus];
  Double_t cdy_out[econst::maxclus];
  Int_t cpid_out[econst::maxclus];

  //reference tdc vars
  Double_t reftdc_out[econst::maxclus];
  Double_t refid_out[econst::maxclus];
  Double_t refmult_out[econst::maxclus];

  // set output tree branches
  P->Branch( "pblkid", &pblkid_out, "pblkid/D" );
  P->Branch( "atime", &atime_out, "atime/D" );

  P->Branch( "bcatime_dxdy", &bcatime_dxdy_out, "bcatime_dxdy/D" );
  P->Branch( "bcatime_atime", &bcatime_atime_out, "bcatime_atime/D" );
  P->Branch( "bcatime_thpq", &bcatime_thpq_out, "bcatime_thpq/D" );
  P->Branch( "bcpblkid_dxdy", &bcpblkid_dxdy_out, "bcpblkid_dxdy/D" );
  P->Branch( "bcpblkid_atime", &bcpblkid_atime_out, "bcpblkid_atime/D" );
  P->Branch( "bcpblkid_thpq", &bcpblkid_thpq_out, "bcpblkid_thpq/D" );

  P->Branch( "dx", &dx_out, "dx/D" );
  P->Branch( "dy", &dy_out, "dy/D" );
  P->Branch( "W2", &W2_out, "W2/D" );
  P->Branch( "Q2", &Q2_out, "Q2/D" );
  P->Branch( "hcale", &hcale_out, "hcale/D" );
  P->Branch( "pse", &pse_out, "pse/D" );
  P->Branch( "she", &pse_out, "she/D" );
  P->Branch( "ep", &ep_out, "ep/D" );
  P->Branch( "hcalatime", &hcalatime_out, "hcalatime/D" );
  P->Branch( "hodotmean", &hodotmean_out, "hodotmean/D" );
  P->Branch( "thetapq_p", &thetapq_pout, "thetapq_pout/D" );
  P->Branch( "thetapq_n", &thetapq_nout, "thetapq_nout/D" );
  P->Branch( "pid", &pid_out, "pid_out/I" );
  P->Branch( "mag", &mag_out, "mag_out/I" );
  P->Branch( "run", &run_out, "run_out/I" );
  P->Branch( "tar", &tar_out, "tar_out/I" );
  P->Branch( "failedglobal", &failedglobal_out, "failedglobal/I" );
  P->Branch( "failedaccmatch", &failedaccmatch_out, "failedaccmatch/I" );
  P->Branch( "failedcoin", &failedcoin_out, "failedcoin/I" );
  P->Branch( "failedw2", &failedw2_out, "failedw2/I" );

  P->Branch( "cpblkid", &cpblkid_out, Form("cpblkid[%d]/D",econst::maxclus) );
  P->Branch( "chatime", &chatime_out, Form("chatime[%d]/D",econst::maxclus) );
  P->Branch( "cthetapq_p", &cthetapq_p_out, Form("cthetapq_p[%d]/D",econst::maxclus) );
  P->Branch( "cthetapq_n", &cthetapq_n_out, Form("cthetapq_n[%d]/D",econst::maxclus) );
  P->Branch( "cdx", &cdx_out, Form("cdx[%d]/D",econst::maxclus) );
  P->Branch( "cdy", &cdy_out, Form("cdy[%d]/D",econst::maxclus) );
  P->Branch( "cpid", &cpid_out, Form("cpid[%d]/I",econst::maxclus) );

  P->Branch( "reftdc", &reftdc_out, Form("reftdc[%d]/D",econst::maxclus) );
  P->Branch( "refid", &refid_out, Form("refid[%d]/D",econst::maxclus) );
  P->Branch( "refmult", &refmult_out, Form("refmult[%d]/D",econst::maxclus) );

  //Set up histogram
  TH2D *htdcvrf = new TH2D("htdcvrf",
			   "hcal tdc vs hcal rftime; tdc (ns); rftime (ns)",
			   400,
			   -300,
			   100,
			   320,
			   -20,
			   300);

  TH2D *htdcvrf_pclus = new TH2D("htdcvrf_pclus",
			   "hcal tdc vs hcal rftime (primary cluster, primary block); tdc (ns); rftime (ns)",
			   400,
			   -300,
			   100,
			   320,
			   -20,
			   300);

  TH2D *htdcvrf_bclus = new TH2D("htdcvrf_bclus",
			   "hcal tdc vs hcal rftime (best cluster, primary block); tdc (ns); rftime (ns)",
			   400,
			   -300,
			   100,
			   320,
			   -20,
			   300);



  // setup reporting indices
  Int_t curmag = -1;
  std::string curtar = "";

  for ( Int_t t=0; t<2; t++ ){ //loop over targets
    // t==0, lh2; t==1, ld2

    for (Int_t irun=0; irun<nruns; irun++) { //loop over runs in each target category
      // accessing run info
      Int_t runnum;
      Int_t mag;
      Double_t ebeam;
      Int_t conf;
      std::string targ;
      std::string rfname;
      if( t==0 ){
	runnum = crunh[irun].runnum;
	mag = crunh[irun].sbsmag / 21; //convert to percent
	ebeam = crunh[irun].ebeam; //get beam energy per run
	conf = crunh[irun].sbsconf; //get config of current run
	targ = crunh[irun].target;
	rfname = rootfile_dir + Form("/*%d*",crunh[irun].runnum);
	bool skip = true;
	for( Int_t el=0; el<lh2r.size(); el++ ){     
	  if( runnum==lh2r[el] ) skip=false;
	}
	if( skip==true ) continue;
      }
      if( t==1 ){
	runnum = crund[irun].runnum;
	mag = crund[irun].sbsmag / 21;
	ebeam = crund[irun].ebeam;
	conf = crund[irun].sbsconf; //get config of current run
	targ = crund[irun].target;
	rfname = rootfile_dir + Form("/*%d*",crund[irun].runnum);
	bool skip = true;
	for( Int_t el=0; el<ld2r.size(); el++ ){     
	  if( runnum==ld2r[el] ) skip=false;
	}
	if( skip==true ) continue;
      }

      std::cout << "Analyzing run " << runnum << ".." << std::endl;

      if( conf!=kine && kine!=-1) continue; //do not proceed if kinematic is constrained

      //set up configuration and tune objects to load analysis parameters
      SBSconfig config(kine,mag);

      //Obtain configuration pars from config file
      Double_t hcaltheta = config.GetHCALtheta_rad();
      Double_t hcaldist = config.GetHCALdist();
      Double_t sbsdist = config.GetSBSdist();
      Double_t bbthr = config.GetBBtheta_rad(); //in radians

      SBStune tune(kine,mag);
    
      if( targ.compare(curtar)!=0 || mag!=curmag ){
	std::cout << "Settings change.." << endl;
	cout << config;
	cout << tune;
	curtar = targ;
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

      // HCal all
      Double_t hcalalltdc[econst::hcalchan], hcalallid[econst::hcalchan];
      Int_t Nhcalallid;
      // std::vector<std::string> hcalallvar = {"tdcelemID","tdc","tdcelemID"};
      // std::vector<void*> hcalallvarlink = {&hcalallid,&hcalalltdc,&hcalallid};
      // rvars::setbranch(C, "sbs.hcal", hcalallvar, hcalallvarlink, 2);

      // HCal general
      Double_t hcalid, hcale, hcalx, hcaly, hcalr, hcalc, hcala, hcalamp, hcaltdc, hcalatime;
      std::vector<std::string> hcalvar = {"idblk","e","x","y","rowblk","colblk","a_p","a_amp_p","tdctimeblk","atimeblk"};
      std::vector<void*> hcalvarlink = {&hcalid,&hcale,&hcalx,&hcaly,&hcalr,&hcalc,&hcala,&hcalamp,&hcaltdc,&hcalatime};
      rvars::setbranch(C, "sbs.hcal", hcalvar, hcalvarlink);

      // HCal cluster branches
      Double_t hcalcid[econst::maxclus], hcalce[econst::maxclus], hcalcx[econst::maxclus], hcalcy[econst::maxclus], hcalctdctime[econst::maxclus], hcalcatime[econst::maxclus];
      Int_t Nhcalcid;
      std::vector<std::string> hcalcvar = {"id","e","x","y","tdctime","atime","id"};
      std::vector<void*> hcalcvarlink = {&hcalcid,&hcalce,&hcalcx,&hcalcy,&hcalctdctime,&hcalcatime,&Nhcalcid};
      rvars::setbranch(C, "sbs.hcal.clus", hcalcvar, hcalcvarlink, 6);

      // HCal reference tdc branches
      Double_t reftdc[econst::maxclus], refid[econst::maxclus], refmult[econst::maxclus];
      Int_t Nrefid;
      std::vector<std::string> hcalrvar = {"tdcelemID","tdc","tdc_mult","tdcelemID"};
      std::vector<void*> hcalrvarlink = {&refid,&reftdc,&refmult,&Nrefid};
      rvars::setbranch(C, "sbs.hcal.Ref", hcalrvar, hcalrvarlink, 3);

      // HCal primary clust block branches
      Double_t hcalcbid[econst::maxclus], hcalcbe[econst::maxclus], hcalcbx[econst::maxclus], hcalcby[econst::maxclus], hcalcbrow[econst::maxclus], hcalcbcol[econst::maxclus], hcalcbtdctime[econst::maxclus], hcalcbatime[econst::maxclus];
      Int_t Nhcalcbid;
      std::vector<std::string> hcalcbvar = {"id","e","row","col","x","y","tdctime","atime","id"};
      std::vector<void*> hcalcbvarlink = {&hcalcbid,&hcalcbe,&hcalcbrow,&hcalcbcol,&hcalcbx,&hcalcby,&hcalcbtdctime,&hcalcbatime,&Nhcalcbid};
      rvars::setbranch(C, "sbs.hcal.clus_blk", hcalcbvar, hcalcbvarlink, 8);

      // bbcal clus var
      Double_t eSH, xSH, ySH, rblkSH, cblkSH, idblkSH, atimeSH, ePS, rblkPS, cblkPS, idblkPS, atimePS;
      std::vector<std::string> bbcalclvar = {"sh.e","sh.x","sh.y","sh.rowblk","sh.colblk","sh.idblk","sh.atimeblk","ps.e","ps.rowblk","ps.colblk","ps.idblk","ps.atimeblk"};
      std::vector<void*> bbcalclvarlink = {&eSH,&xSH,&ySH,&rblkSH,&cblkSH,&idblkSH,&atimeSH,&ePS,&rblkPS,&cblkPS,&idblkPS,&atimePS};
      rvars::setbranch(C, "bb", bbcalclvar, bbcalclvarlink);

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
      std::cout << "Uncorrected average beam energy on " << targ << " for run: " << ebeam << std::endl;
      //set up hcal coordinate system with hcal angle wrt exit beamline
      vector<TVector3> hcalaxes; vars::sethcalaxes( hcaltheta, hcalaxes );
      TVector3 hcalorigin = hcaldist*hcalaxes[2] + econst::hcalvoff*hcalaxes[0];
      Double_t BdL = econst::sbsmaxfield * econst::sbsdipolegap * (mag/100); //scale crudely by magnetic field
      Double_t Eloss_outgoing = econst::celldiameter/2.0/sin(bbthr) * econst::lh2tarrho * econst::lh2dEdx;

      // set nucleon defaults by target
      std::string nucleon;
      if( targ.compare("LH2")==0 )
	nucleon = "p"; 
      else if( targ.compare("LD2")==0 )      
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
	Double_t ebeam_o = vars::ebeam_o( ebeam_c, etheta, targ ); //Second energy correction accounting for energy loss leaving target

	//v3
	nu = pbeam.E() - pcent;
	pNexp = vars::pN_expect( nu, nucleon );
	thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	pNhat = vars::pNhat_track( thNexp, phNexp );
	pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	Q2 = vars::Q2( pbeam.E(), pe.E(), etheta );
	W2 = vars::W2( pbeam.E(), pe.E(), W2, nucleon );


	//Calculate h-arm quantities
	vector<Double_t> xyhcalexp; vars::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
	TVector3 hcalpos = hcalorigin + hcalcbx[0]*hcalaxes[0] + hcalcby[0]*hcalaxes[1]; //from primary blk
	Double_t dx = hcalcbx[0] - xyhcalexp[0];
	Double_t dy = hcalcby[0] - xyhcalexp[1];
	TVector3 neutdir = ( hcalpos - vertex ).Unit();
	Double_t protdeflect = tan( 0.3 * BdL / q.Mag() ) * (hcaldist - (sbsdist + econst::sbsdipolegap/2.0) );
	TVector3 protdir = ( hcalpos + protdeflect * hcalaxes[0] - vertex ).Unit();
	Double_t thetapq_p = acos( protdir.Dot( qv.Unit() ) );
	Double_t thetapq_n = acos( neutdir.Dot( qv.Unit() ) );

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

	/////////////////////////////////////////////////
	//Primary W2 cut on elastics and HCal active area
	bool failedw2 = fabs(W2-W2mean)>W2sig; //Observed mean W2 cut on elastic peak
	Int_t fw2int = 0;
	if( failedw2 ) fw2int = 1;
      
	//Fill all block physics output - this is inefficient, should get primary cluster info from here first
	Int_t cidx_atime = 0;
	Double_t c_atimediff = 1000.;
	Int_t cidx_thetapq = 0;
	Double_t c_thetapqdiff = 1000.;
	Int_t cidx_dxdy = 0;
	Double_t c_dxdydiff = 1000.;
	Double_t tpq_best;

	///Check clusters for best
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

	  //get cluster which minimizes thetapq (proton for this analysis), passes coin, and pid well defined
	  if( cthetapq_p<c_thetapqdiff && abs(hcalcatime[c]-atime0)<3*atimesig && (cpid==1||cpid==2) ){
	    cidx_thetapq = c;
	    if( cpid==1 )
	      c_thetapqdiff = cthetapq_p;
	    if( cpid==2 )
	      c_thetapqdiff = cthetapq_n;

	  }

	  cpblkid_out[c] = hcalcid[c];
	  chatime_out[c] = hcalcatime[c];
	  cthetapq_p_out[c] = cthetapq_p;
	  cthetapq_n_out[c] = cthetapq_n;
	  cdx_out[c] = cdx;
	  cdy_out[c] = cdy;
	  
	  cpid_out[c] = cpid;
	}

	//HCal RF time in tdc reference channel 3
	Double_t hcal_rf_time = reftdc[3];
	//loop over tdc reference channels
	for( Int_t c=0; c<Nrefid; c++ ){
	  reftdc_out[c] = reftdc[c];
	  refid_out[c] = refid[c];
	  refmult_out[c] = refmult[c];

	}

	//fill simple analysis histogram for block 161
	if( !failedglobal && !failedaccmatch && !failedcoin && !failedw2 ){

	  //Fill from all blocks (tdc does not require inclusion in cluster)
	  htdcvrf->Fill(hcalalltdc[161],hcal_rf_time);
	  //Fill from primary cluster on tree
	  if( hcalid==161 )
	    htdcvrf_pclus->Fill(hcaltdc,hcal_rf_time);
	  //Fill from best cluster (atime cut, then minimize thetapq)
	  if( hcalcid[cidx_thetapq]==161 )
	    htdcvrf_bclus->Fill(hcalctdctime[cidx_thetapq],hcal_rf_time);

	}

	//Calculate dx, dy, thetapq from the best cluster for later use
	Double_t dx_bestcluster = hcalcx[cidx_thetapq] - xyhcalexp[0];
	Double_t dy_bestcluster = hcalcy[cidx_thetapq] - xyhcalexp[1];
	Double_t hatime_bestcluster = hcalcatime[cidx_thetapq];

	pblkid_out = (double)hcalcbid[0];
	atime_out = hcalcbatime[0];
   
	bcatime_dxdy_out = hcalcatime[cidx_dxdy];
	bcatime_atime_out = hcalcatime[cidx_atime];
	bcatime_thpq_out = hcalcatime[cidx_thetapq];
	bcpblkid_dxdy_out = (double)hcalcbid[cidx_dxdy];
	bcpblkid_atime_out = (double)hcalcbid[cidx_atime];
	bcpblkid_thpq_out = (double)hcalcbid[cidx_thetapq];

	dx_out = dx;
	dy_out = dy;
	W2_out = W2;
	Q2_out = Q2;
	hcale_out = hcale;
	pse_out = ePS;
	she_out = eSH;
	ep_out = p[0];
	hcalatime_out = hcalcbatime[0];
	hodotmean_out = hodotmean[0];
	thetapq_pout = thetapq_p ;
	thetapq_nout = thetapq_n;
	mag_out = mag;
	run_out = runnum;
	tar_out = t;
	pid_out = pid;
	failedglobal_out = fglobalint;
	failedaccmatch_out = faccmatchint;
	failedcoin_out = fcoinint;
	failedw2_out = fw2int;


	P->Fill();

      }//end loop over event

      // getting ready for the next run
      C->Reset();

    }//endloop over runs

  }//endloop over targets

  fout->Write();

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
