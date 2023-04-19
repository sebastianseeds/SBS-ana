//sseeds 4.10.23 - Short script to parse a kinematic and produce a physics output tree

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

//
void parse_sh( Int_t kine=14, Int_t epm=2 )
{ //main  

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/parse_sh.json");

  std::string rootfile_dir = jmgr->GetValueFromSubKey_str( "rootfile_dir", Form("sbs%d",kine));
  Int_t pass = 1;

  if( pass>1 ){
    std::cout << "As of 3.31.23, the highest GMn replay pass is 1. Enter a valid pass." << endl;
    return;
  }

  //set up default parameters for all analysis
  std::string runsheet_dir = "/w/halla-scshelf2102/sbs/seeds/ana/data"; //unique to my environment for now
  Int_t nhruns = -1; //Always analyze all available runs
  Int_t ndruns = -1; //Always analyze all available runs
  Int_t verb = 0; //Don't print diagnostic info by default
  Double_t binfac = 400.;

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<crun> crunh; 
  vector<crun> crund;
  util::ReadRunList(runsheet_dir,nhruns,kine,"LH2",pass,verb,crunh); //modifies nruns to be very large when -1
  util::ReadRunList(runsheet_dir,ndruns,kine,"LD2",pass,verb,crund);

  //set up output files
  TFile *fout = new TFile( Form("outfiles/shortparse_sbs%d.root",kine), "RECREATE" );

  //set up diagnostic histograms
  TH2D *hW2mag = new TH2D( "hW2mag", "W^{2} vs sbsmag; \%; GeV^{2}", 20, 0, 100, 200, 0, 2 );
  TH2D *hQ2mag = new TH2D( "hQ2mag", "Q^{2} vs sbsmag; \%; GeV^{2}", 20, 0, 100, 200, 0, 4 ); //can generalize
  TH2D *hdxmag = new TH2D( "hdxmag","dx vs sbsmag; \%; x_{HCAL}-x_{expect} (m)", 20, 0, 100, 800, -4, 4 );
  TH2D *hdymag = new TH2D( "hdymag","dy vs sbsmag; \%; y_{HCAL}-y_{expect} (m)", 20, 0, 100, 800, -4, 4 );

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
  Double_t hcale_out;
  Double_t pse_out;
  Double_t she_out;
  Double_t ep_out;
  Double_t eoverp_out;
  Double_t hcalatime_out;
  Double_t hodotmean_out;
  Double_t thetapq_pout;
  Double_t thetapq_nout;
  Int_t run_out;
  Int_t tar_out; //0: LH2, 1:LD2
  Int_t mag_out;
  Int_t failedglobal_out;
  Int_t failedaccmatch_out;
  Int_t failedcoin_out;

  // set output tree branches
  P->Branch( "dx", &dx_out, "dx/D" );
  P->Branch( "dy", &dy_out, "dy/D" );
  P->Branch( "W2", &W2_out, "W2/D" );
  P->Branch( "Q2", &Q2_out, "Q2/D" );
  P->Branch( "nu", &nu_out, "nu/D" );
  P->Branch( "hcale", &hcale_out, "hcale/D" );
  P->Branch( "pse", &pse_out, "pse/D" );
  P->Branch( "she", &she_out, "she/D" );
  P->Branch( "ep", &ep_out, "ep/D" );
  P->Branch( "eoverp", &eoverp_out, "eoverp/D" );
  P->Branch( "hcalatime", &hcalatime_out, "hcalatime/D" );
  P->Branch( "hodotmean", &hodotmean_out, "hodotmean/D" );
  P->Branch( "thetapq_p", &thetapq_pout, "thetapq_p/D" );
  P->Branch( "thetapq_n", &thetapq_nout, "thetapq_n/D" );
  P->Branch( "mag", &mag_out, "mag/I" );
  P->Branch( "run", &run_out, "run/I" );
  P->Branch( "tar", &tar_out, "tar/I" );
  P->Branch( "failedglobal", &failedglobal_out, "failedglobal/I" );
  P->Branch( "failedaccmatch", &failedaccmatch_out, "failedaccmatch/I" );
  P->Branch( "failedcoin", &failedcoin_out, "failedcoin/I" );

  // setup reporting indices
  Int_t curmag = -1;
  std::string curtar = "";
  
  for ( Int_t t=0; t<2; t++ ){ //loop over targets
    // t==0, lh2; t==1, ld2
    Int_t nruns;
    if( t==0 ) nruns = nhruns;
    if( t==1 ) nruns = ndruns;

    for (Int_t irun=0; irun<nruns; irun++) {
      // accessing run info
      Int_t runnum;
      Int_t mag;
      Int_t ebeam;
      std::string targ;
      std::string rfname;
      if( t==0 ){
	runnum = crunh[irun].runnum;
	mag = crunh[irun].sbsmag / 21; //convert to percent
	ebeam = crunh[irun].ebeam; //get beam energy per run
	targ = crunh[irun].target;
	rfname = rootfile_dir + Form("/*%d*",crunh[irun].runnum);
      }
      if( t==1 ){
	runnum = crund[irun].runnum;
	mag = crund[irun].sbsmag / 21;
	ebeam = crund[irun].ebeam;
	targ = crund[irun].target;
	rfname = rootfile_dir + Form("/*%d*",crund[irun].runnum);
      }
      std::cout << "Analyzing " << targ << " run " << runnum << ".." << std::endl;

      //set up configuration and tune objects to load analysis parameters
      SBSconfig config(kine,mag);

      //Obtain configuration pars from config file
      Double_t hcaltheta = config.GetHCALtheta_rad();
      Double_t hcaldist = config.GetHCALdist();
      Double_t sbsdist = config.GetSBSdist();
    
      //SBStune *tune = new SBStune(kine,mag);
      SBStune tune(kine,mag);
    
      //Reporting. tar should always equal curtar as categorized by good run list
      if( targ.compare(curtar)!=0 || mag!=curmag ){
	std::cout << "Settings change.." << endl;
	cout << config;
	cout << tune;
	curtar = targ;
	curmag = mag;
      }
      
      //Obtain cuts from tune class
    std:string gcut   = tune.Getglobcut();
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
    
      // bbcal clus var
      Double_t eSH, xSH, ySH, rblkSH, cblkSH, idblkSH, atimeSH, ePS, rblkPS, cblkPS, idblkPS, atimePS;
      std::vector<std::string> bbcalclvar = {"sh.e","sh.x","sh.y","sh.rowblk","sh.colblk","sh.idblk","sh.atimeblk","ps.e","ps.rowblk","ps.colblk","ps.idblk","ps.atimeblk"};
      std::vector<void*> bbcalclvarlink = {&eSH,&xSH,&ySH,&rblkSH,&cblkSH,&idblkSH,&atimeSH,&ePS,&rblkPS,&cblkPS,&idblkPS,&atimePS};
      rvars::setbranch(C, "bb", bbcalclvar, bbcalclvarlink);

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

      TCut GCut = gcut.c_str();
      TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );

      // get experimental quantities by run
      std::cout << "Uncorrected average beam energy on " << targ << " for run: " << ebeam << std::endl;
      //set up hcal coordinate system with hcal angle wrt exit beamline
      vector<TVector3> hcalaxes; vars::sethcalaxes( hcaltheta, hcalaxes );
      TVector3 hcalorigin = hcaldist*hcalaxes[2] + econst::hcalvoff*hcalaxes[0];
      Double_t BdL = econst::sbsmaxfield * econst::sbsdipolegap * (mag/100); //scale crudely by magnetic field

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

      while (C->GetEntry(nevent++)) {

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
	bool failedcoin = abs( hcalatime - atime0 ) > 3*atimesig;

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
	  W2 = vars::W2( pbeam.E(), pe.E(), Q2, nucleon );
	}else if( epm==4 ){
	  //v4
	  nu = pbeam.E() - pe.E();
	  pNexp = vars::pN_expect( nu, nucleon );
	  thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	  pNhat = vars::qVect_unit( thNexp, phNexp );
	  pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	  Q2 = vars::Q2( pbeam.E(), pe.E(), etheta );
	  W2 = vars::W2( pbeam.E(), pe.E(), Q2, nucleon );
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
	TVector3 hcalpos = hcalorigin - hcalx*hcalaxes[0] + hcaly*hcalaxes[1];
	Double_t dx = hcalx - xyhcalexp[0];
	Double_t dy = hcaly - xyhcalexp[1];
	TVector3 neutdir = (hcalpos - vertex).Unit();
	Double_t protdeflect = tan( 0.3 * BdL / q.Mag() ) * (hcaldist - (sbsdist + econst::sbsdipolegap/2.0) );
	TVector3 protdir = ( hcalpos + protdeflect * hcalaxes[0] - vertex ).Unit();
	Double_t thetapq_p = acos( protdir.Dot( pNhat ) );
	Double_t thetapq_n = acos( neutdir.Dot( pNhat ) );
	Double_t eoverp = (ePS + eSH) / p[0];

	//Fill diagnostic histos
	hQ2mag->Fill( mag, Q2 );
	hW2mag->Fill( mag, W2 );
	hdxmag->Fill( mag, dx );
	hdymag->Fill( mag, dy );

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

	//Fill physics output tree     
	dx_out = dx;
	dy_out = dy;
	W2_out = W2;
	Q2_out = Q2;
	nu_out = nu;
	hcale_out = hcale;
	pse_out = ePS;
	she_out = eSH;
	ep_out = p[0];
	eoverp_out = eoverp;
	hcalatime_out = hcalatime;
	hodotmean_out = hodotmean[0];
	thetapq_pout = thetapq_p;
	thetapq_nout = thetapq_n;
	mag_out = mag;
	run_out = runnum;
	tar_out = t;
	failedglobal_out = fglobalint;
	failedaccmatch_out = faccmatchint;
	failedcoin_out = fcoinint;

	P->Fill();

      }//end event loop

      // getting ready for the next run
      C->Reset();

    }//end run loop

  }//end target loop


  fout->Write();

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
