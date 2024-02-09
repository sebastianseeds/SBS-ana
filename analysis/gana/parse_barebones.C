//sseeds 10.23.23 - Updated parsing script to cut both inelastic events and unused branches. Configured only to parse data files, not MC. Event parsing using wide globalcuts, wide W2 cuts, and wide coin (HCal/BBCal) cuts. Branch parsing includes only branches that sseeds is using for his gmn analysis
//Update 2.2.24 - Same method, made simpler without class references for troubleshooting

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

const Int_t maxClus = 35; //larger than max per tree
const Int_t maxBlk = 25;
const Double_t W2max = 3.0; //very large compared to nucleon mass
const Double_t coin_sigma_factor = 5.; //Wide coincidence timing cut

//Specific wide cut for all parsing
const std::string gcut = "bb.ps.e>0.1&&abs(bb.tr.vz[0])<0.12&&bb.sh.nclus>0&&sbs.hcal.nclus>0";

//MAIN
void parse_barebones( Int_t kine=4, Int_t pass=2, Int_t cluster_method=4, bool verbose=false )
{   

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  //One set of data files in json config shared between pass0/1 per kinematic
  if( pass==0 )
    pass=1;

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/parse.json");

  std::string rootfile_dir_lh2 = jmgr->GetValueFromSubKey_str( Form("rootfile_dir_lh2_p%d",pass), Form("sbs%d",kine) );;
  std::string rootfile_dir_ld2 = jmgr->GetValueFromSubKey_str( Form("rootfile_dir_ld2_p%d",pass), Form("sbs%d",kine) );
  Double_t minE = jmgr->GetValueFromSubKey<Double_t>( Form("minE_p%d",pass), Form("sbs%d",kine) );
  vector<Double_t> coin_profile;
  jmgr->GetVectorFromSubKey<Double_t>(Form("coin_profile_p%d",pass),Form("sbs%d",kine),coin_profile);
  Double_t hcal_v_offset = jmgr->GetValueFromSubKey<Double_t>( Form("hcal_offset_p%d",pass), Form("sbs%d",kine) );

  //set up default parameters for all analysis
  std::string runsheet_dir = "/w/halla-scshelf2102/sbs/seeds/ana/data"; //unique to my environment for now
  Int_t nhruns = -1; //Always analyze all available runs
  Int_t ndruns = -1; //Always analyze all available runs
  Int_t verb = 0; //Don't print diagnostic info by default
  Double_t binfac = 400.;

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<crun> crunh; 
  vector<crun> crund;
  util::ReadRunList(runsheet_dir,nhruns,kine,"LH2",pass,verb,crunh); //modifies nruns
  util::ReadRunList(runsheet_dir,ndruns,kine,"LD2",pass,verb,crund); //modifies nruns

  // outfile path
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string parse_path = outdir_path + Form("/parse/parse_sbs%d_pass%d_barebones.root",kine,pass);

  //set up output files
  TFile *fout = new TFile( parse_path.c_str(), "RECREATE" );

  //set up diagnostic histograms
  TH2D *hW2mag = new TH2D( "hW2mag", "W^{2} vs sbsmag; \%; GeV^{2}", 20, 0, 100, 200, 0, 2 );
  TH2D *hQ2mag = new TH2D( "hQ2mag", "Q^{2} vs sbsmag; \%; GeV^{2}", 20, 0, 100, 200, 0, 4 ); //can generalize
  TH2D *hdxmag = new TH2D( "hdxmag","dx vs sbsmag; \%; x_{HCAL}-x_{expect} (m)", 20, 0, 100, 800, -4, 4 );
  TH2D *hdymag = new TH2D( "hdymag","dy vs sbsmag; \%; y_{HCAL}-y_{expect} (m)", 20, 0, 100, 800, -4, 4 );
  TH1D *hcidxfail_d = new TH1D( "hcidxfail_d","Cluster Index Failed vs Run; runnum", ndruns, 0, ndruns );
  TH1D *hcidxfail_h = new TH1D( "hcidxfail_h","Cluster Index Failed vs Run; runnum", nhruns, 0, nhruns );

  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;

  // create output tree
  TTree *P = new TTree("P","Analysis Tree"); 

  // new output tree vars
  //Universals
  Int_t run_out;
  Int_t tar_out; //0:LH2, 1:LD2
  Int_t mag_out;
  Int_t event_out;
  Int_t trig_out;
  Double_t xexp_out;
  Double_t yexp_out;
  Double_t W2_out;
  Double_t Q2_out;
  Double_t nu_out;
  Double_t precon_out;

  //Fiducial slices
  Int_t failedfid_0_0_out;
  Int_t failedfid_0_5_out;
  Int_t failedfid_1_0_out;
  Int_t failedfid_1_5_out;
  Int_t failedfid_2_0_out;
  Int_t failedfid_2_5_out;
  Int_t failedfid_3_0_out;
  Int_t failedfid_3_5_out;

  //Primary cluster
  Double_t dx_out;
  Double_t dy_out;
  Double_t coin_out;
  Double_t thetapq_pout;
  Double_t thetapq_nout;
  Int_t nucleon_out; //via data ellipse inclusion. -1=outside both, 0=inside both, 1=proton, 2=neutron
  Int_t hcalon_out;
  Double_t hcalnblk_out;
  Double_t hcalpid_out;
  Double_t hcalx_out;
  Double_t hcaly_out;
  Double_t hcale_out;  

  //Best cluster
  Double_t dx_bc_out;
  Double_t dy_bc_out;
  Double_t coin_bc_out;
  Double_t thetapq_bc_pout;
  Double_t thetapq_bc_nout;
  Int_t nucleon_bc_out; //via data ellipse inclusion. -1=outside both, 0=inside both, 1=proton, 2=neutron
  Int_t hcalon_bc_out;
  Double_t hcalnblk_bc_out;
  Double_t hcalpid_bc_out;
  Double_t hcalx_bc_out;
  Double_t hcaly_bc_out;
  Double_t hcale_bc_out;  

  // relevant old output tree vars
  Double_t bb_tr_vz_out;
  Double_t bb_tr_p_out;
  Double_t bb_ps_e_out;
  Double_t bb_ps_rowblk_out;
  Double_t bb_ps_colblk_out;
  Double_t bb_sh_e_out;
  Double_t bb_sh_rowblk_out;
  Double_t bb_sh_colblk_out;
  Double_t bb_hodotdc_clus_tmean_out;
  Double_t bb_gem_track_nhits_out;
  Double_t bb_etot_over_p_out;

  P->Branch("run_out", &run_out, "run_out/I");
  P->Branch("tar_out", &tar_out, "tar_out/I");
  P->Branch("mag_out", &mag_out, "mag_out/I");
  P->Branch("event_out", &event_out, "event_out/I");
  P->Branch("trig_out", &trig_out, "trig_out/I");
  P->Branch("xexp_out", &xexp_out, "xexp_out/D");
  P->Branch("yexp_out", &yexp_out, "yexp_out/D");
  P->Branch("W2_out", &W2_out, "W2_out/D");
  P->Branch("Q2_out", &Q2_out, "Q2_out/D");
  P->Branch("nu_out", &nu_out, "nu_out/D");
  P->Branch("precon_out", &precon_out, "precon_out/D");

  P->Branch("failedfid_0_0_out", &failedfid_0_0_out, "failedfid_0_0_out/I");
  P->Branch("failedfid_0_5_out", &failedfid_0_5_out, "failedfid_0_5_out/I");
  P->Branch("failedfid_1_0_out", &failedfid_1_0_out, "failedfid_1_0_out/I");
  P->Branch("failedfid_1_5_out", &failedfid_1_5_out, "failedfid_1_5_out/I");
  P->Branch("failedfid_2_0_out", &failedfid_2_0_out, "failedfid_2_0_out/I");
  P->Branch("failedfid_2_5_out", &failedfid_2_5_out, "failedfid_2_5_out/I");
  P->Branch("failedfid_3_0_out", &failedfid_3_0_out, "failedfid_3_0_out/I");
  P->Branch("failedfid_3_5_out", &failedfid_3_5_out, "failedfid_3_5_out/I");

  P->Branch("dx_out", &dx_out, "dx_out/D");
  P->Branch("dy_out", &dy_out, "dy_out/D");
  P->Branch("coin_out", &coin_out, "coin_out/D");
  P->Branch("thetapq_pout", &thetapq_pout, "thetapq_pout/D");
  P->Branch("thetapq_nout", &thetapq_nout, "thetapq_nout/D");
  P->Branch("nucleon_out", &nucleon_out, "nucleon_out/I");
  P->Branch("hcalon_out", &hcalon_out, "hcalon_out/I");
  P->Branch("hcalnblk_out", &hcalnblk_out, "hcalnblk_out/D");
  P->Branch("hcalpid_out", &hcalpid_out, "hcalpid_out/D");
  P->Branch("hcalx_out", &hcalx_out, "hcalx_out/D");
  P->Branch("hcaly_out", &hcaly_out, "hcaly_out/D");
  P->Branch("hcale_out", &hcale_out, "hcale_out/D");

  P->Branch("dx_bc_out", &dx_bc_out, "dx_bc_out/D");
  P->Branch("dy_bc_out", &dy_bc_out, "dy_bc_out/D");
  P->Branch("coin_bc_out", &coin_bc_out, "coin_bc_out/D");
  P->Branch("thetapq_bc_pout", &thetapq_bc_pout, "thetapq_bc_pout/D");
  P->Branch("thetapq_bc_nout", &thetapq_bc_nout, "thetapq_bc_nout/D");
  P->Branch("nucleon_bc_out", &nucleon_bc_out, "nucleon_bc_out/I");
  P->Branch("hcalon_bc_out", &hcalon_bc_out, "hcalon_bc_out/I");
  P->Branch("hcalnblk_bc_out", &hcalnblk_bc_out, "hcalnblk_bc_out/D");
  P->Branch("hcalpid_bc_out", &hcalpid_bc_out, "hcalpid_bc_out/D");
  P->Branch("hcalx_bc_out", &hcalx_bc_out, "hcalx_bc_out/D");
  P->Branch("hcaly_bc_out", &hcaly_bc_out, "hcaly_bc_out/D");
  P->Branch("hcale_bc_out", &hcale_bc_out, "hcale_bc_out/D");

  P->Branch("bb_tr_vz_out", &bb_tr_vz_out, "bb_tr_vz_out/D");
  P->Branch("bb_tr_p_out", &bb_tr_p_out, "bb_tr_p_out/D");
  P->Branch("bb_ps_e_out", &bb_ps_e_out, "bb_ps_e_out/D");
  P->Branch("bb_ps_rowblk_out", &bb_ps_rowblk_out, "bb_ps_rowblk_out/D");
  P->Branch("bb_ps_colblk_out", &bb_ps_colblk_out, "bb_ps_colblk_out/D");
  P->Branch("bb_sh_e_out", &bb_sh_e_out, "bb_sh_e_out/D");
  P->Branch("bb_sh_rowblk_out", &bb_sh_rowblk_out, "bb_sh_rowblk_out/D");
  P->Branch("bb_sh_colblk_out", &bb_sh_colblk_out, "bb_sh_colblk_out/D");
  P->Branch("bb_hodotdc_clus_tmean_out", &bb_hodotdc_clus_tmean_out, "bb_hodotdc_clus_tmean_out/D");
  P->Branch("bb_gem_track_nhits_out", &bb_gem_track_nhits_out, "bb_gem_track_nhits_out/D");
  P->Branch("bb_etot_over_p_out", &bb_etot_over_p_out, "bb_etot_over_p_out/D");

  // setup reporting indices
  Int_t curmag = -1;
  std::string curtar = "";

  for ( Int_t t=0; t<2; t++ ){ //loop over targets
    // t==0, lh2; t==1, ld2
    Int_t nruns = 1;
    if( t==0 ){
      nruns = nhruns;
      std::cout << std::endl << "Proceeding to LH2 data.." << std::endl << std::endl;
    }
    if( t==1 ){
      nruns = ndruns;
      std::cout << std::endl << "Proceeding to LD2 data.." << std::endl << std::endl;
    }

    for (Int_t irun=0; irun<nruns; irun++) {
      
      cout << "irun/nruns " << irun << "/" << nruns << endl;

      // accessing run info
      Int_t runnum;
      Int_t mag;
      Double_t ebeam;
      Double_t charge;
      std::string targ;
      std::string rfname;
      if( t==0 ){
	//std::cout << crunh;
	runnum = crunh[irun].runnum;
	mag = crunh[irun].sbsmag / 21; //convert to percent
	ebeam = crunh[irun].ebeam; //get beam energy per run
	targ = crunh[irun].target;
	rfname = rootfile_dir_lh2 + Form("/*%d*",crunh[irun].runnum);
	charge = crunh[irun].charge;
      }
      if( t==1 ){
	//std::cout << crund;
	runnum = crund[irun].runnum;
	mag = crund[irun].sbsmag / 21;
	ebeam = crund[irun].ebeam;
	targ = crund[irun].target;
	rfname = rootfile_dir_ld2 + Form("/*%d*",crund[irun].runnum);
	charge = crund[irun].charge;
      }

      //set up configuration and tune objects to load analysis parameters
      SBSconfig config(kine,mag);

      //Obtain configuration pars from config file
      Double_t hcaltheta = config.GetHCALtheta_rad();
      Double_t hcaldist = config.GetHCALdist();
      Double_t sbsdist = config.GetSBSdist();
      Double_t bbthr = config.GetBBtheta_rad(); //in radians

      SBStune tune(kine,mag);

      //Reporting. tar should always equal curtar as categorized by good run list
      if( targ.compare(curtar)!=0 || mag!=curmag ){
	if(verbose){
	  std::cout << "Settings change.." << std::endl;
	  std::cout << config;
	  //std::cout << tune;
	}
	curtar = targ;
	curmag = mag;
      }

      //Obtain cuts from tune class
      //std::string gcut   = tune.Getglobcut_wide();
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

      //first attempt to resolve segmentation fault on large data sets
      if (C != nullptr) {
	delete C;
      }

      C = new TChain("T");
      C->Add(rfname.c_str());

      // setting up ROOT tree branch addresses
      C->SetBranchStatus("*",0);    

      // HCal general
      Double_t hcalid, hcale, hcalx, hcaly, hcalr, hcalc, hcaltdc, hcalatime, hcalidx, nclus, nblk;
      std::vector<std::string> hcalvar = {"idblk","e","x","y","rowblk","colblk","tdctimeblk","atimeblk","index","nclus","nblk"};
      std::vector<void*> hcalvarlink = {&hcalid,&hcale,&hcalx,&hcaly,&hcalr,&hcalc,&hcaltdc,&hcalatime,&hcalidx,&nclus,&nblk};
      rvars::setbranch(C, "sbs.hcal", hcalvar, hcalvarlink);
      
      // HCal cluster branches
      Double_t hcalcid[econst::maxclus], hcalce[econst::maxclus], hcalcx[econst::maxclus], hcalcy[econst::maxclus], hcalctdctime[econst::maxclus], hcalcatime[econst::maxclus], hcalcrow[econst::maxclus], hcalcnblk[econst::maxclus], hcalccol[econst::maxclus];
      Int_t Nhcalcid;
      std::vector<std::string> hcalcvar = {"id","e","x","y","tdctime","atime","row","col","nblk","id"};
      std::vector<void*> hcalcvarlink = {&hcalcid,&hcalce,&hcalcx,&hcalcy,&hcalctdctime,&hcalcatime,&hcalcrow,&hcalccol,&hcalcnblk,&Nhcalcid};
      rvars::setbranch(C, "sbs.hcal.clus", hcalcvar, hcalcvarlink, 9);

      // HCal cluster blk branches
      Double_t hcalcbid[econst::maxclus], hcalcbe[econst::maxclus], hcalcbx[econst::maxclus], hcalcby[econst::maxclus], hcalcbtdctime[econst::maxclus], hcalcbatime[econst::maxclus];
      Int_t Nhcalcbid;
      std::vector<std::string> hcalcbvar = {"id","e","x","y","tdctime","atime","id"};
      std::vector<void*> hcalcbvarlink = {&hcalcbid,&hcalcbe,&hcalcbx,&hcalcby,&hcalcbtdctime,&hcalcbatime,&Nhcalcbid};
      rvars::setbranch(C, "sbs.hcal.clus_blk", hcalcbvar, hcalcbvarlink, 6);

      // bbcal clus var
      Double_t eSH, xSH, ySH, rblkSH, cblkSH, idblkSH, atimeSH, nclusSH, ePS, rblkPS, cblkPS, idblkPS, atimePS;
      std::vector<std::string> bbcalclvar = {"sh.e","sh.x","sh.y","sh.rowblk","sh.colblk","sh.idblk","sh.atimeblk","sh.nclus","ps.e","ps.rowblk","ps.colblk","ps.idblk","ps.atimeblk"};
      std::vector<void*> bbcalclvarlink = {&eSH,&xSH,&ySH,&rblkSH,&cblkSH,&idblkSH,&atimeSH,&nclusSH,&ePS,&rblkPS,&cblkPS,&idblkPS,&atimePS};
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

      // other bb branches
      Double_t gemNhits, eop;
      std::vector<std::string> miscbbvar = {"gem.track.nhits","etot_over_p"};
      std::vector<void*> miscbbvarlink = {&gemNhits,&eop};
      rvars::setbranch(C, "bb", miscbbvar, miscbbvarlink);

      TCut GCut = gcut.c_str();

      TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );

      // get experimental quantities by run
      //set up hcal coordinate system with hcal angle wrt exit beamline
      vector<TVector3> hcalaxes; vars::sethcalaxes( hcaltheta, hcalaxes );
      TVector3 hcalorigin = hcaldist*hcalaxes[2] + hcal_v_offset*hcalaxes[0];

      cout << "ATTENTION: Loaded vertical offset of " << hcal_v_offset << endl;

      Double_t BdL = econst::sbsmaxfield * econst::sbsdipolegap * (mag/100); //scale crudely by magnetic field
      Double_t Eloss_outgoing;
      if(t==0)
	Eloss_outgoing = econst::celldiameter/2.0/sin(bbthr) * econst::lh2tarrho * econst::lh2dEdx;
      if(t==1)
	Eloss_outgoing = econst::celldiameter/2.0/sin(bbthr) * econst::ld2tarrho * econst::ld2dEdx;

      //Set up hcal active area with bounds that match database on pass
      vector<Double_t> hcalaa;
      if(pass<1)
	hcalaa = cut::hcalaa_data_alt(1,1);
      else
	hcalaa = cut::hcalaa_data(1,1);

      // set nucleon defaults by target
      std::string nucleon;
      if( targ.compare("LH2")==0 )
	nucleon = "p"; 
      else if( targ.compare("LD2")==0 )      
	nucleon = "np";
      else{
	nucleon = "np";
	std::cout << "ERROR: incorrect target information loaded from good run list. Check and verify." << std::endl;
      }

      // event indices
      long nevent = 0, npassed = 0, nevents = C->GetEntries(); 
      Int_t treenum = 0, currenttreenum = 0;

      if(verbose)
	std::cout << "Beginning analysis of run " << runnum << ", target " << targ << ", magnetic field " << mag << "%, total accumulated charge " << charge << " C." << std::endl;
      
      while (C->GetEntry(nevent++)) {
	
	std::cout << "Processing run " << runnum << " event " << nevent << " / " << nevents << ", total passed cuts " << npassed << "\r";
	std::cout.flush();

	if((Int_t)hcalidx>9){
	  if(t==0)
	    hcidxfail_h->Fill(irun);
	  if(t==1)
	    hcidxfail_d->Fill(irun);
	  continue;
	}

	///////
	//Single-loop globalcut method. Save pass/fail for output tree.
	bool failedglobal = false;

	currenttreenum = C->GetTreeNumber();
	if( nevent == 1 || currenttreenum != treenum ){
	  treenum = currenttreenum; 
	  GlobalCut->UpdateFormulaLeaves();
	  //cout << "Updating formula leaves and switching segment at event: " << nevent << endl;
	}
	failedglobal = GlobalCut->EvalInstance(0) == 0;
	if( failedglobal ){
	  continue;
	}

	///////
	//Physics calculations
	//correct beam energy with vertex information
	Double_t ebeam_c;
	ebeam_c = vars::ebeam_c( ebeam, vz[0], targ );

	TVector3 vertex( 0., 0., vz[0] );

	//reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
	Double_t precon = p[0] + Eloss_outgoing;

	//set up four-momenta with some empty for various calculation methods
	TLorentzVector pbeam( 0., 0., ebeam_c, ebeam_c ); //beam momentum
	//TLorentzVector pe( px[0], py[0], pz[0], p[0] ); //e' momentum
	TLorentzVector pe( precon*px[0]/p[0], precon*py[0]/p[0], precon*pz[0]/p[0], precon ); //e' recon plvect
	TLorentzVector ptarg; vars::setPN(nucleon,ptarg); //target momentum
	TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
	TLorentzVector pN; //N' momentum
	TVector3 pNhat; //Unit N' 3-vector
      
	//simple calculations for e' and N'
	Double_t etheta = vars::etheta(pe); 
	Double_t ephi = vars::ephi(pe);
	Double_t pcent = vars::pcentral(ebeam,etheta,nucleon); //e' p reconstructed by angles
	Double_t phNexp = ephi + physconst::pi;
	Double_t Q2, W2, nu, thNexp, pNexp, ebeam_o;
	ebeam_o = vars::ebeam_o( ebeam_c, etheta, targ ); //Second energy correction accounting for energy loss leaving target
       
	//reconstruct track with angles (better resolution with GEMs)
	nu = pbeam.E() - pcent;
	pNexp = vars::pN_expect( nu, nucleon );
	thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	pNhat = vars::pNhat_track( thNexp, phNexp );
	pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	Q2 = vars::Q2( pbeam.E(), pe.E(), etheta );
	W2 = vars::W2( pbeam.E(), pe.E(), Q2, nucleon );
	
	/////////////////////
	//W2 elastic cut
	bool failedW2 = W2>W2max;
	if(failedW2)
	  continue;

	npassed++;

	Double_t comp_ev_fraction = (Double_t)npassed/(Double_t)nevent;
	Double_t ev_fraction = (Double_t)npassed/(Double_t)nevents;

	//Calculate h-arm quantities
	vector<Double_t> xyhcalexp; vars::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
	TVector3 hcalpos = hcalorigin + hcalx*hcalaxes[0] + hcaly*hcalaxes[1];
	Double_t dx = hcalx - xyhcalexp[0];
	Double_t dy = hcaly - xyhcalexp[1];
	TVector3 neutdir = ( hcalpos - vertex ).Unit();
	Double_t protdeflect = tan( 0.3 * BdL / q.Mag() ) * (hcaldist - (sbsdist + econst::sbsdipolegap/2.0) );
	TVector3 protdir = ( hcalpos + protdeflect * hcalaxes[0] - vertex ).Unit();
	Double_t thetapq_p = acos( protdir.Dot( pNhat ) );
	Double_t thetapq_n = acos( neutdir.Dot( pNhat ) );
	Double_t eoverp = ( ePS + eSH ) / p[0];

	double pclus_diff = hcalatime - atimeSH;
	double pclus_nblk = hcalcnblk[0];
	double pblkid = hcalcid[0];

	//Determine rough PID
	bool in_p = util::Nspotcheck(dy,dx,dy0,dx0_p,dysig,dxsig_p,0);
	bool in_n = util::Nspotcheck(dy,dx,dy0,dx0_n,dysig,dxsig_n,0);

	Int_t PID = -1;
	if(in_p&&in_n)
	  PID=0;
	else if(in_p)
	  PID=1;
	else if(in_n)
	  PID=2;

	//////////////////////
	//ALL CLUSTER ANALYSIS
      
	//Set up clone clusters for selection analysis.
	vector<double> clone_cluster_intime;
	vector<double> clone_cluster_score;

	//loop through all clusters and select without HCal position information
	for( int c=0; c<Nhcalcid; c++ ){
	
	  //calculate h-arm physics quantities per cluster
	  double atime = hcalcatime[c];
	  double atime_diff = atime - atimeSH; //Assuming best shower time on primary cluster
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
      
	if( hcale != hcalce[(Int_t)hcalidx] && verbose )
	  cerr << "WARNING: Sorting failure on cluster index " << (Int_t)hcalidx << ". Check max clusters allowed on tree." << endl;
	

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
      
	switch (cluster_method) {
	case 1:
	  cidx_best = 0;
	  break;
	case 2:
	  cidx_best = cidx_e;
	  break;
	case 3:
	  cidx_best = cidx_intime;
	  break;
	case 4:  //default, best elastic yield without bias
	  cidx_best = cidx_score;
	  break;
	default:
	  cidx_best = 3;
	}

	//Calculations from the best cluster
	TVector3 hcalpos_bestcluster = hcalorigin + 
	  hcalcx[cidx_best]*hcalaxes[0] + 
	  hcalcy[cidx_best]*hcalaxes[1];

	TVector3 protdir_bc = ( hcalpos_bestcluster + protdeflect * hcalaxes[0] - vertex ).Unit();
	TVector3 neutdir_bc = ( hcalpos_bestcluster - vertex ).Unit();
	Double_t thetapq_bc_p = acos( protdir_bc.Dot( pNhat ) );
	Double_t thetapq_bc_n = acos( neutdir_bc.Dot( pNhat ) );

	Double_t dx_bestcluster = hcalcx[cidx_best] - xyhcalexp[0];
	Double_t dy_bestcluster = hcalcy[cidx_best] - xyhcalexp[1];
	Double_t hatime_bestcluster = hcalcatime[cidx_best];
	Double_t hcoin_bestcluster = hcalcatime[cidx_best] - atimeSH;
	Double_t ce_bestcluster = hcalce[cidx_best];
	Double_t x_bestcluster = hcalcx[cidx_best];
	Double_t y_bestcluster = hcalcy[cidx_best];
	Double_t row_bestcluster = hcalcrow[cidx_best];
	Double_t col_bestcluster = hcalccol[cidx_best];
	Double_t nblk_bestcluster = hcalcnblk[cidx_best];
	Double_t pblkid_bestcluster = hcalcid[cidx_best];

	//Determine rough PID
	bool in_p_bc = util::Nspotcheck(dy_bestcluster,dx_bestcluster,dy0,dx0_p,dysig,dxsig_p,0);
	bool in_n_bc = util::Nspotcheck(dy_bestcluster,dx_bestcluster,dy0,dx0_n,dysig,dxsig_n,0);

	Int_t PID_bc = -1;
	if(in_p_bc&&in_n_bc)
	  PID_bc=0;
	else if(in_p_bc)
	  PID_bc=1;
	else if(in_n_bc)
	  PID_bc=2;


	//H-arm fiducial cuts (best_cluster)
	bool hcalON_bc = cut::hcalaaON(hcalcx[cidx_best],hcalcy[cidx_best],hcalaa);

	//Set up multiple fiducial cut points
	//Nsig=0
	vector<Double_t> fid_0_0 = cut::hcalfid(dxsig_p,dysig,hcalaa,0);
	bool passed_fid_0_0 = cut::hcalfidIN(xyhcalexp[0],xyhcalexp[1],dx0_p,fid_0_0);
	//Nsig=0.5
	vector<Double_t> fid_0_5 = cut::hcalfid(dxsig_p,dysig,hcalaa,0.5);
	bool passed_fid_0_5 = cut::hcalfidIN(xyhcalexp[0],xyhcalexp[1],dx0_p,fid_0_5);
	//Nsig=1.
	vector<Double_t> fid_1_0 = cut::hcalfid(dxsig_p,dysig,hcalaa,1.);
	bool passed_fid_1_0 = cut::hcalfidIN(xyhcalexp[0],xyhcalexp[1],dx0_p,fid_1_0);
	//Nsig=1.5
	vector<Double_t> fid_1_5 = cut::hcalfid(dxsig_p,dysig,hcalaa,1.5);
	bool passed_fid_1_5 = cut::hcalfidIN(xyhcalexp[0],xyhcalexp[1],dx0_p,fid_1_5);
	//Nsig=2.
	vector<Double_t> fid_2_0 = cut::hcalfid(dxsig_p,dysig,hcalaa,2.);
	bool passed_fid_2_0 = cut::hcalfidIN(xyhcalexp[0],xyhcalexp[1],dx0_p,fid_2_0);
	//Nsig=2.5
	vector<Double_t> fid_2_5 = cut::hcalfid(dxsig_p,dysig,hcalaa,2.5);
	bool passed_fid_2_5 = cut::hcalfidIN(xyhcalexp[0],xyhcalexp[1],dx0_p,fid_2_5);
	//Nsig=3
	vector<Double_t> fid_3_0 = cut::hcalfid(dxsig_p,dysig,hcalaa,3.);
	bool passed_fid_3_0 = cut::hcalfidIN(xyhcalexp[0],xyhcalexp[1],dx0_p,fid_3_0);
	//Nsig=3.5
	vector<Double_t> fid_3_5 = cut::hcalfid(dxsig_p,dysig,hcalaa,3.5);
	bool passed_fid_3_5 = cut::hcalfidIN(xyhcalexp[0],xyhcalexp[1],dx0_p,fid_3_5);

	//H-arm acceptance cut (primary cluster)
	bool hcalON = cut::hcalaaON(hcalx,hcaly,hcalaa);

	//H-arm acceptance cut (primary cluster)
	//bool hcalON_bc = cut::hcalaaON(x_bestcluster,y_bestcluster,hcalaa);

	//Fill diagnostic histos
	hQ2mag->Fill( mag, Q2 );
	hW2mag->Fill( mag, W2 );
	hdxmag->Fill( mag, dx );
	hdymag->Fill( mag, dy );
	
	////////////////////
	//coincidence time BBCal/HCal cut
	bool failedwidecoin_bc = abs( hatime_bestcluster - atimediff0 ) > coin_sigma_factor*atimediffsig;

	//Fill new output tree  
	//Universals
	mag_out = mag;
	run_out = runnum;
	tar_out = t;
	event_out = (Int_t)gevnum;
	trig_out = (Int_t)trigbits;
	xexp_out = xyhcalexp[0];
	yexp_out = xyhcalexp[1];
	W2_out = W2;
	Q2_out = Q2;
	nu_out = nu;
	precon_out = precon;

	//Fiducial slices
	failedfid_0_0_out = (Int_t)!passed_fid_0_0;
	failedfid_0_5_out = (Int_t)!passed_fid_0_5;
	failedfid_1_0_out = (Int_t)!passed_fid_1_0;
	failedfid_1_5_out = (Int_t)!passed_fid_1_5;
	failedfid_2_0_out = (Int_t)!passed_fid_2_0;
	failedfid_2_5_out = (Int_t)!passed_fid_2_5;
	failedfid_3_0_out = (Int_t)!passed_fid_3_0;
	failedfid_3_5_out = (Int_t)!passed_fid_3_5;

	//Primary Cluster
	dx_out = dx;
	dy_out = dy;
	coin_out = pclus_diff;
	thetapq_pout = thetapq_p;
	thetapq_nout = thetapq_n;
	nucleon_out = PID;
	hcalnblk_out = pclus_nblk;
	hcalon_out = hcalON;
	hcalpid_out = pblkid;
	hcalx_out = hcalx;
	hcaly_out = hcaly;
	hcale_out = hcale;

	//Best Cluster
	dx_bc_out = dx_bestcluster;
	dy_bc_out = dy_bestcluster;
	coin_bc_out = hcoin_bestcluster;
	thetapq_bc_pout = thetapq_bc_p;
	thetapq_bc_nout = thetapq_bc_n;
	nucleon_bc_out = PID_bc;
	hcalnblk_bc_out = nblk_bestcluster;
	hcalon_bc_out = hcalON_bc;
	hcalpid_bc_out = pblkid_bestcluster;
	hcalx_bc_out = x_bestcluster;
	hcaly_bc_out = y_bestcluster;
	hcale_bc_out = ce_bestcluster;

	//Fill old output tree
	bb_tr_p_out = p[0];
	bb_tr_vz_out = vz[0];
	bb_ps_e_out = ePS;
	bb_ps_rowblk_out = rblkPS;
	bb_ps_colblk_out = cblkPS;
	bb_sh_e_out = eSH;
	bb_sh_rowblk_out = rblkSH;
	bb_sh_colblk_out = cblkSH;
	bb_hodotdc_clus_tmean_out = hodotmean[0];
	bb_gem_track_nhits_out = gemNhits;
	bb_etot_over_p_out = eop;

	P->Fill();	
	
      }//end event loop

      // getting ready for the next run
      C->Reset();

    }//end run loop

  }//end target loop

  fout->Write();

  std::cout << std::endl << "Barebones parsing complete. Output written to " << parse_path << std::endl << std::endl;

  st->Stop();

  // Send time efficiency report to console
  std::cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << std::endl;

}
