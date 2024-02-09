//sseeds 10.23.23 - Updated parsing script to cut both inelastic events and unused branches. Configured only to parse data files, not MC. Event parsing using wide globalcuts, wide W2 cuts, and wide coin (HCal/BBCal) cuts. Branch parsing includes only branches that sseeds is using for his gmn analysis

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

//MAIN
void parse( Int_t kine=7, Int_t epm=3, Int_t cluster_method = 4, Int_t pass=0, bool verbose=false )
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

  // if( pass>1 ){
  //   std::cout << "As of 10.31.23, the highest GMn replay pass is 1. Enter a valid pass." << endl;
  //   return;
  // }

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
  std::string parse_path = outdir_path + Form("/parse/parse_sbs%d_pass%d_update.root",kine,pass);

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
  TChain *TS = nullptr;

  // create output tree
  TTree *P = new TTree("P","Analysis Tree"); 

  // new output tree vars
  Double_t dx_out;
  Double_t dy_out;
  Double_t xexp_out;
  Double_t yexp_out;
  Double_t W2_out;
  Double_t Q2_out;
  Double_t nu_out;
  Double_t eoverp_out;
  Double_t thetapq_pout;
  Double_t thetapq_nout;
  Double_t charge_ts_out;
  Int_t nucleon_bc_out; //via data ellipse inclusion. -1=outside both, 0=inside both, 1=proton, 2=neutron
  Int_t run_out;
  Int_t tar_out; //0:LH2, 1:LD2
  Int_t mag_out;
  Int_t event_out;
  Int_t trig_out;
  Int_t idx_bc_out;
  Int_t dx_bc_out;
  Int_t dy_bc_out;
  Int_t failedwidecoin_bc_out;
  Int_t failedfid_0_0_out;
  Int_t failedfid_0_5_out;
  Int_t failedfid_1_0_out;
  Int_t failedfid_1_5_out;
  Int_t failedfid_2_0_out;
  Int_t failedfid_2_5_out;
  Int_t failedfid_3_0_out;
  Int_t failedfid_3_5_out;
  Int_t hcalon_out;
  Int_t hcalon_bc_out;
  //Int_t cluster_idx_fail_out;
  Double_t charge_avg_out;
  Double_t comp_event_frac_out;
  Double_t event_frac_out;
  Double_t acc_charge_out;

  // relevant old output tree vars
  Double_t bb_tr_vx_out;
  Double_t bb_tr_vy_out;
  Double_t bb_tr_vz_out;
  Double_t bb_tr_p_out;
  Double_t bb_tr_px_out;
  Double_t bb_tr_py_out;
  Double_t bb_tr_pz_out;
  Double_t bb_tr_n_out;
  Double_t e_kine_Q2_out;
  Double_t e_kine_W2_out;
  Double_t e_kine_nu_out;
  Double_t bb_ps_e_out;
  Double_t bb_ps_rowblk_out;
  Double_t bb_ps_colblk_out;
  Double_t bb_sh_e_out;
  Double_t bb_sh_rowblk_out;
  Double_t bb_sh_colblk_out;
  Double_t bb_sh_atimeblk_out;
  Double_t bb_sh_nclus_out;
  Double_t bb_hodotdc_clus_tmean_out;
  Double_t bb_gem_track_nhits_out;
  Double_t bb_etot_over_p_out;
  Double_t sbs_hcal_e_out;
  Double_t sbs_hcal_x_out;
  Double_t sbs_hcal_y_out;
  Double_t sbs_hcal_rowblk_out;
  Double_t sbs_hcal_colblk_out;
  Double_t sbs_hcal_atimeblk_out;
  Double_t sbs_hcal_tdctimeblk_out;
  Double_t sbs_hcal_clus_id_out[maxClus];
  Double_t sbs_hcal_clus_e_out[maxClus];
  Double_t sbs_hcal_clus_x_out[maxClus];
  Double_t sbs_hcal_clus_y_out[maxClus];
  Double_t sbs_hcal_clus_tdctime_out[maxClus];
  Double_t sbs_hcal_clus_atime_out[maxClus];
  Double_t sbs_hcal_nclus_out;
  Double_t sbs_hcal_nblk_out;
  Double_t sbs_hcal_clus_blk_id_out[maxBlk];
  Double_t sbs_hcal_clus_blk_e_out[maxBlk];
  Double_t sbs_hcal_clus_blk_x_out[maxBlk];
  Double_t sbs_hcal_clus_blk_y_out[maxBlk];
  Double_t sbs_hcal_clus_blk_atime_out[maxBlk];
  Double_t sbs_hcal_clus_blk_tdctime_out[maxBlk];
  Int_t Ndata_sbs_hcal_clus_blk_id_out;
  Int_t Ndata_sbs_hcal_clus_id_out;
  Double_t sbs_hcal_bclus_x_out;
  Double_t sbs_hcal_bclus_y_out;
  Double_t sbs_hcal_bclus_e_out;
  Double_t sbs_hcal_bclus_row_out;
  Double_t sbs_hcal_bclus_col_out;
  Double_t sbs_hcal_bclus_atime_out;

  // set new output tree branches
  P->Branch( "dx", &dx_out, "dx/D" );
  P->Branch( "dy", &dy_out, "dy/D" );
  P->Branch( "xexp", &xexp_out, "xexp/D" );
  P->Branch( "yexp", &yexp_out, "yexp/D" );
  P->Branch( "W2", &W2_out, "W2/D" );
  P->Branch( "Q2", &Q2_out, "Q2/D" );
  P->Branch( "nu", &nu_out, "nu/D" );
  P->Branch( "eoverp", &eoverp_out, "eoverp/D" );
  P->Branch( "thetapq_p", &thetapq_pout, "thetapq_p/D" );
  P->Branch( "thetapq_n", &thetapq_nout, "thetapq_n/D" );
  P->Branch( "charge_ts", &charge_ts_out, "charge_ts/D" );
  P->Branch( "nucleon_bc", &nucleon_bc_out, "nucleon_bc/I" );
  P->Branch( "mag", &mag_out, "mag/I" );
  P->Branch( "run", &run_out, "run/I" );
  P->Branch( "tar", &tar_out, "tar/I" );
  P->Branch( "event", &event_out, "event/I" );
  P->Branch( "trig", &trig_out, "trig/I" );
  P->Branch( "idx_bc", &idx_bc_out, "idx_bc/I" );
  P->Branch( "dx_bc", &dx_bc_out, "dx_bc/D" );
  P->Branch( "dy_bc", &dy_bc_out, "dy_bc/D" );
  P->Branch( "failedwidecoin_bc", &failedwidecoin_bc_out, "failedwidecoin_bc/I" );
  P->Branch( "failedfid_0_0", &failedfid_0_0_out, "failedfid_0_0/I" );
  P->Branch( "failedfid_0_5", &failedfid_0_5_out, "failedfid_0_5/I" );
  P->Branch( "failedfid_1_0", &failedfid_1_0_out, "failedfid_1_0/I" );
  P->Branch( "failedfid_1_5", &failedfid_1_5_out, "failedfid_1_5/I" );
  P->Branch( "failedfid_2_0", &failedfid_2_0_out, "failedfid_2_0/I" );
  P->Branch( "failedfid_2_5", &failedfid_2_5_out, "failedfid_2_5/I" );
  P->Branch( "failedfid_3_0", &failedfid_3_0_out, "failedfid_3_0/I" );
  P->Branch( "failedfid_4_5", &failedfid_3_5_out, "failedfid_3_5/I" );
  P->Branch( "hcalon", &hcalon_out, "hcalon/I" );
  P->Branch( "hcalon_bc", &hcalon_bc_out, "hcalon_bc/I" );
  //P->Branch( "cluster_idx_fail", &cluster_idx_fail_out, "cluster_idx_fail/I" );
  P->Branch( "charge_avg", &charge_avg_out, "charge_avg/D" );
  P->Branch( "comp_event_frac", &comp_event_frac_out, "comp_event_frac/D" );
  P->Branch( "event_frac", &event_frac_out, "event_frac/D" );
  P->Branch( "sbs.hcal.bclus.x", &sbs_hcal_bclus_x_out, "sbs.hcal.bclus.x/D" );
  P->Branch( "sbs.hcal.bclus.y", &sbs_hcal_bclus_y_out, "sbs.hcal.bclus.y/D" );
  P->Branch( "sbs.hcal.bclus.e", &sbs_hcal_bclus_e_out, "sbs.hcal.bclus.e/D" );
  P->Branch( "sbs.hcal.bclus.row", &sbs_hcal_bclus_row_out, "sbs.hcal.bclus.row/D" );
  P->Branch( "sbs.hcal.bclus.col", &sbs_hcal_bclus_col_out, "sbs.hcal.bclus.col/D" );
  P->Branch( "sbs.hcal.bclus.atime", &sbs_hcal_bclus_atime_out, "sbs.hcal.bclus.atime/D" );

  // set relevant old output tree branches
  P->Branch( "bb.tr.vx", &bb_tr_vx_out, "bb.tr.vx/D" );
  P->Branch( "bb.tr.vy", &bb_tr_vy_out, "bb.tr.vy/D" );
  P->Branch( "bb.tr.vz", &bb_tr_vz_out, "bb.tr.vz/D" );
  P->Branch( "bb.tr.p", &bb_tr_p_out, "bb.tr.p/D" );
  P->Branch( "bb.tr.px", &bb_tr_px_out, "bb.tr.px/D" );
  P->Branch( "bb.tr.py", &bb_tr_py_out, "bb.tr.py/D" );
  P->Branch( "bb.tr.pz", &bb_tr_pz_out, "bb.tr.pz/D" );
  P->Branch( "bb.tr.n", &bb_tr_n_out, "bb.tr.n/D" );
  P->Branch( "e.kine.Q2", &e_kine_Q2_out, "e.kine.Q2/D" );
  P->Branch( "e.kine.W2", &e_kine_W2_out, "e.kine.W2/D" );
  P->Branch( "e.kine.nu", &e_kine_nu_out, "e.kine.nu/D" );
  P->Branch( "bb.ps.e", &bb_ps_e_out, "bb.ps.e/D" );
  P->Branch( "bb.ps.rowblk", &bb_ps_rowblk_out, "bb.ps.rowblk/D" );
  P->Branch( "bb.ps.colblk", &bb_ps_colblk_out, "bb.ps.colblk/D" );
  P->Branch( "bb.sh.e", &bb_sh_e_out, "bb.sh.e/D" );
  P->Branch( "bb.sh.rowblk", &bb_sh_rowblk_out, "bb.sh.rowblk/D" );
  P->Branch( "bb.sh.colblk", &bb_sh_colblk_out, "bb.sh.colblk/D" );
  P->Branch( "bb.sh.atimeblk", &bb_sh_atimeblk_out, "bb.sh.atimeblk/D" );
  P->Branch( "bb.sh.nclus", &bb_sh_nclus_out, "bb.sh.nclus/D" );
  P->Branch( "bb.hodotdc.clus.tmean", &bb_hodotdc_clus_tmean_out, "bb.hodotdc.clus.tmean/D" );
  P->Branch( "bb.gem.track.nhits", &bb_gem_track_nhits_out, "bb.gem.track.nhits/D" );
  P->Branch( "bb.etot_over_p", &bb_etot_over_p_out, "bb.etot_over_p/D" );
  P->Branch( "sbs.hcal.e", &sbs_hcal_e_out, "sbs.hcal.e/D" );
  P->Branch( "sbs.hcal.x", &sbs_hcal_x_out, "sbs.hcal.x/D" );
  P->Branch( "sbs.hcal.y", &sbs_hcal_y_out, "sbs.hcal.y/D" );
  P->Branch( "sbs.hcal.rowblk", &sbs_hcal_rowblk_out, "sbs.hcal.rowblk/D" );
  P->Branch( "sbs.hcal.colblk", &sbs_hcal_colblk_out, "sbs.hcal.colblk/D" );
  P->Branch( "sbs.hcal.atimeblk", &sbs_hcal_atimeblk_out, "sbs.hcal.atimeblk/D" );
  P->Branch( "sbs.hcal.clus.id", &sbs_hcal_clus_id_out, "sbs.hcal.clus.id/D" );
  P->Branch( "sbs.hcal.clus.e", &sbs_hcal_clus_e_out, "sbs.hcal.clus.e/D" );
  P->Branch( "sbs.hcal.clus.x", &sbs_hcal_clus_x_out, "sbs.hcal.clus.x/D" );
  P->Branch( "sbs.hcal.clus.y", &sbs_hcal_clus_y_out, "sbs.hcal.clus.y/D" );
  P->Branch( "sbs.hcal.clus.tdctime", &sbs_hcal_clus_tdctime_out, "sbs.hcal.clus.tdctime/D" );
  P->Branch( "sbs.hcal.clus.atime", &sbs_hcal_clus_atime_out, "sbs.hcal.clus.atime/D" );
  P->Branch( "sbs.hcal.clus_blk.id", &sbs_hcal_clus_blk_id_out, "sbs.hcal.clus_blk.id/D" );
  P->Branch( "sbs.hcal.clus_blk.e", &sbs_hcal_clus_blk_e_out, "sbs.hcal.clus_blk.e/D" );
  P->Branch( "sbs.hcal.clus_blk.x", &sbs_hcal_clus_blk_x_out, "sbs.hcal.clus_blk.x/D" );
  P->Branch( "sbs.hcal.clus_blk.y", &sbs_hcal_clus_blk_y_out, "sbs.hcal.clus_blk.y/D" );
  P->Branch( "sbs.hcal.clus_blk.tdctime", &sbs_hcal_clus_blk_tdctime_out, "sbs.hcal.clus_blk.tdctime/D" );
  P->Branch( "sbs.hcal.clus_blk.atime", &sbs_hcal_clus_blk_atime_out, "sbs.hcal.clus_blk.atime/D" );
  P->Branch( "sbs.hcal.nclus", &sbs_hcal_nclus_out, "sbs.hcal.nclus/D" );
  P->Branch( "sbs.hcal.nblk", &sbs_hcal_nblk_out, "sbs.hcal.nblk/D" );
  P->Branch( "Ndata.sbs.hcal.clus.id", &Ndata_sbs_hcal_clus_id_out, "Ndata.sbs.hcal.clus.id/I" );
  P->Branch( "Ndata.sbs.hcal.clus_blk.id", &Ndata_sbs_hcal_clus_blk_id_out, "Ndata.sbs.hcal.clus_blk.id/I" );

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
      
      

      cout << "here 1 with irun/nruns " << irun << "/" << nruns << endl;

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
     
      std::string label = Form("Run: %d",runnum);

      cout << "here 2" << endl;


      // if(!verbose)
      // 	util::updateProgress(irun + 1, nruns, label);

      //analyze bcm from TSsbs for charge accumulated per recorded event number
      typedef struct{
	Double_t enumber;
	Double_t clock_diff;
	Double_t avg_current;
	Double_t acc_charge;
      } SCALARAGG;

      vector<SCALARAGG> Sagg;

      cout << "here 3" << endl;


      //first attempt to resolve segmentation fault on large data sets
      if (TS != nullptr) {
	delete TS;
      }

      cout << "here 4" << endl;


      //Create chain for epics/scalar tree
      TS = new TChain("TSsbs");
      TS->Add( rfname.c_str() );

      cout << "here 5" << endl;


      //Declare dummies and counters for needed vars
      Double_t Sevnum = 0., Sdnew = 0., Sclk = 0., clkLast = 0., curLast = 0., curAvg = 0., clkDiff = 0.;
      Int_t DumTN = 0, DumCurTN = 0, DumEvLast = 0, DumEvOff = 0, DumRun = 0, evDiff = 0, evTCoff = 0;
      Long64_t CurEv = 0;
      Long64_t NumEv = TS->GetEntries();
      bool treechange = false;
      bool begin = true;
      Double_t chargeLast = 0;
      Double_t chargeAvgLast = 0;
      //Switch on and set addresses for branches
      TS->SetBranchStatus( "*", 0 );
      TS->SetBranchStatus( "evNumber", 1 );
      TS->SetBranchStatus( "sbs.bcm.dnew.cnt", 1 );
      TS->SetBranchStatus( "sbs.104kHz_CLK.cnt", 1 );
      TS->SetMakeClass(1);
      TS->SetBranchAddress( "evNumber", &Sevnum );
      TS->SetBranchAddress( "sbs.bcm.dnew.cnt", &Sdnew );
      TS->SetBranchAddress( "sbs.104kHz_CLK.cnt", &Sclk );
      //Loop over events (dt=2s for epics tree)
      //Get rate (Sdnew/3318/clkLast) then multiply by the time (clkLast), without the offset, 3318 is cnts/C
      while( TS->GetEntry( CurEv++ ) ){

	Double_t clk = Sclk/103700;
	if(begin)
	  clkLast=clk;
	clkDiff = clk-clkLast;
	clkLast = clk;

	Double_t cur = Sdnew/3318;
	if(begin){
	  curLast=cur;
	  begin=false;
	}
	curAvg = (cur+curLast)/2;
	curLast = cur;

	chargeLast = curAvg*clkDiff*10E-9; //charge = current*time, convert from ns->s

	//Fill the scalar data structure with charge and clock corrections applied
	SCALARAGG thisSCALAR = { Sevnum, clkDiff, curAvg, chargeLast };
	Sagg.push_back( thisSCALAR );

	//cout << "SCALARAGG vals on CurEv " << CurEv << ": " << Sagg.back().enumber << " " << Sagg.back().aggcharge << " " << Sagg.back().clock << endl;

      }

      cout << "here 6" << endl;


      //set up configuration and tune objects to load analysis parameters
      SBSconfig config(kine,mag);

      //Obtain configuration pars from config file
      Double_t hcaltheta = config.GetHCALtheta_rad();
      Double_t hcaldist = config.GetHCALdist();
      Double_t sbsdist = config.GetSBSdist();
      Double_t bbthr = config.GetBBtheta_rad(); //in radians

      //SBStune *tune = new SBStune(kine,mag);
      SBStune tune(kine,mag);
    
      cout << "here 7" << endl;


      //Reporting. tar should always equal curtar as categorized by good run list
      if( targ.compare(curtar)!=0 || mag!=curmag ){
	if(verbose){
	  std::cout << "Settings change.." << std::endl;
	  std::cout << config;
	  std::cout << tune;
	}
	curtar = targ;
	curmag = mag;
      }
      
      cout << "here 7" << endl;


      //Obtain cuts from tune class
      std:string gcut   = tune.Getglobcut_wide();
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

      cout << "here 8" << endl;


      //first attempt to resolve segmentation fault on large data sets
      if (C != nullptr) {
	delete C;
      }

      cout << "here 9" << endl;


      C = new TChain("T");
      C->Add(rfname.c_str());

      // setting up ROOT tree branch addresses
      C->SetBranchStatus("*",0);    

      cout << "here 10" << endl;


      // HCal general
      Double_t hcalid, hcale, hcalx, hcaly, hcalr, hcalc, hcaltdc, hcalatime, hcalidx, nclus, nblk;
      std::vector<std::string> hcalvar = {"idblk","e","x","y","rowblk","colblk","tdctimeblk","atimeblk","index","nclus","nblk"};
      std::vector<void*> hcalvarlink = {&hcalid,&hcale,&hcalx,&hcaly,&hcalr,&hcalc,&hcaltdc,&hcalatime,&hcalidx,&nclus,&nblk};
      rvars::setbranch(C, "sbs.hcal", hcalvar, hcalvarlink);

      // HCal cluster branches
      Double_t hcalcid[econst::maxclus], hcalce[econst::maxclus], hcalcx[econst::maxclus], hcalcy[econst::maxclus], hcalctdctime[econst::maxclus], hcalcatime[econst::maxclus], hcalcrow[econst::maxclus], hcalccol[econst::maxclus];
      Int_t Nhcalcid;
      std::vector<std::string> hcalcvar = {"id","e","x","y","tdctime","atime","row","col","id"};
      std::vector<void*> hcalcvarlink = {&hcalcid,&hcalce,&hcalcx,&hcalcy,&hcalctdctime,&hcalcatime,&hcalcrow,&hcalccol,&Nhcalcid};
      rvars::setbranch(C, "sbs.hcal.clus", hcalcvar, hcalcvarlink, 8);

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


      cout << "here 11" << endl;


      // globalcut branches
      //C->SetBranchStatus("bb.gem.track.nhits", 1);
      //C->SetBranchStatus("bb.etot_over_p", 1);
      //C->SetBranchStatus("bb.tr.n", 1);

      TCut GCut = gcut.c_str();

      TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );

      // get experimental quantities by run
      //std::cout << "Uncorrected average beam energy on " << targ << " for run: " << ebeam << std::endl;
      //set up hcal coordinate system with hcal angle wrt exit beamline
      vector<TVector3> hcalaxes; vars::sethcalaxes( hcaltheta, hcalaxes );
      //TVector3 hcalorigin = hcaldist*hcalaxes[2] + econst::hcalvoff*hcalaxes[0];
      TVector3 hcalorigin = hcaldist*hcalaxes[2] + hcal_v_offset*hcalaxes[0];
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

      Double_t TSlastEv = -1.;
      Double_t TSnextEv = Sagg[0].enumber;
      Double_t TSlastt = Sagg[0].clock_diff;
      Double_t TSlastCur = Sagg[0].avg_current;
      Double_t TSlastC = chargeAvgLast; //Just use last average for first few events wherever necessary.
      Double_t TSstructidx=-1;

      if(verbose)
	std::cout << "Beginning analysis of run " << runnum << ", target " << targ << ", magnetic field " << mag << "%, total accumulated charge " << charge << " C." << std::endl;

      cout << "here 12" << endl;

      //cout << "N events " << nevents << endl;

      while (C->GetEntry(nevent++)) {
	
	//cout << "processing event " << nevent << endl;

	//std::cout << "Processing run " << runnum << " event " << nevent << " / " << nevents << "\r";
	//std::cout.flush();
	
	//std::cout << "Processing run " << runnum << " event " << nevent << " / " << nevents << endl;

	// if((Int_t)hcalidx>9){
	//   if(t==0)
	//     hcidxfail_h->Fill(irun);
	//   if(t==1)
	//     hcidxfail_d->Fill(irun);
	//   continue;
	// }

	// if(!verbose)
	//   util::updateProgress(irun + 1, nevent +1, nruns, nevents, label);

	//Access TS charge struct members
	if( gevnum>=TSnextEv ){
	  TSstructidx++;
	  TSlastEv = Sagg[TSstructidx].enumber;
	  TSnextEv = Sagg[TSstructidx+1].enumber;
	  TSlastt = Sagg[TSstructidx].clock_diff;
	  TSlastCur = Sagg[TSstructidx].avg_current;
	  TSlastC = Sagg[TSstructidx].acc_charge;
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
	if( failedglobal )
	  continue;


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
	  //v3 (uses BB track angles, default)
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
	  std::cout << "WARNING: epm version incorrect. Defaulting to version 2." << endl;
	}


	/////////////////////
	//W2 elastic cut
	bool failedW2 = W2>W2max;
	if(failedW2)
	  continue;

	npassed++;

	Double_t comp_ev_fraction = (Double_t)npassed/(Double_t)nevent;
	Double_t ev_fraction = (Double_t)npassed/(Double_t)nevents;
	Double_t accumulated_charge = charge*ev_fraction;

	//cout << npassed << " " << nevent << " " << nevents << " " << ev_fraction << " " << comp_ev_fraction << " " << accumulated_charge << endl;


	//Calculate h-arm quantities
	vector<Double_t> xyhcalexp; vars::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
	TVector3 hcalpos = hcalorigin + hcalx*hcalaxes[0] + hcaly*hcalaxes[1];
	Double_t dx = hcalx - xyhcalexp[0];
	Double_t dy = hcaly - xyhcalexp[1];
	TVector3 neutdir = (hcalpos - vertex).Unit();
	Double_t protdeflect = tan( 0.3 * BdL / q.Mag() ) * (hcaldist - (sbsdist + econst::sbsdipolegap/2.0) );
	TVector3 protdir = ( hcalpos + protdeflect * hcalaxes[0] - vertex ).Unit();
	Double_t thetapq_p = acos( protdir.Dot( pNhat ) );
	Double_t thetapq_n = acos( neutdir.Dot( pNhat ) );
	Double_t eoverp = (ePS + eSH) / p[0];


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
	Double_t dx_bestcluster = hcalcx[cidx_best] - xyhcalexp[0];
	Double_t dy_bestcluster = hcalcy[cidx_best] - xyhcalexp[1];
	Double_t hatime_bestcluster = hcalcatime[cidx_best];
	Double_t hcoin_bestcluster = hcalcatime[cidx_best] - atimeSH;
	Double_t ce_bestcluster = hcalce[cidx_best];
	Double_t x_bestcluster = hcalcx[cidx_best];
	Double_t y_bestcluster = hcalcy[cidx_best];
	Double_t row_bestcluster = hcalcrow[cidx_best];
	Double_t col_bestcluster = hcalccol[cidx_best];


	//Determine rough PID
	bool in_p = util::Nspotcheck(dy_bestcluster,dx_bestcluster,dy0,dx0_p,dysig,dxsig_p,0);
	bool in_n = util::Nspotcheck(dy_bestcluster,dx_bestcluster,dy0,dx0_n,dysig,dxsig_n,0);

	Int_t PID = -1;
	if(in_p&&in_n)
	  PID=0;
	else if(in_p)
	  PID=1;
	else if(in_n)
	  PID=2;


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


	//Fill diagnostic histos
	hQ2mag->Fill( mag, Q2 );
	hW2mag->Fill( mag, W2 );
	hdxmag->Fill( mag, dx );
	hdymag->Fill( mag, dy );
	

	////////////////////
	//coincidence time BBCal/HCal cut
	bool failedwidecoin_bc = abs( hatime_bestcluster - atimediff0 ) > coin_sigma_factor*atimediffsig;


	//Fill new output tree     
	dx_out = dx;
	dy_out = dy;
	xexp_out = xyhcalexp[0];
	yexp_out = xyhcalexp[1];
	W2_out = W2;
	Q2_out = Q2;
	nu_out = nu;
	eoverp_out = eoverp;
	thetapq_pout = thetapq_p;
	thetapq_nout = thetapq_n;
	charge_ts_out = TSlastC;
	nucleon_bc_out = PID;
	mag_out = mag;
	run_out = runnum;
	tar_out = t;
	event_out = (Int_t)gevnum;
	trig_out = (Int_t)trigbits;
	idx_bc_out = cidx_best;
	dx_bc_out = dx_bestcluster;
	dy_bc_out = dy_bestcluster;
	failedwidecoin_bc_out = (Int_t)failedwidecoin_bc;
	failedfid_0_0_out = (Int_t)!passed_fid_0_0;
	failedfid_0_5_out = (Int_t)!passed_fid_0_5;
	failedfid_1_0_out = (Int_t)!passed_fid_1_0;
	failedfid_1_5_out = (Int_t)!passed_fid_1_5;
	failedfid_2_0_out = (Int_t)!passed_fid_2_0;
	failedfid_2_5_out = (Int_t)!passed_fid_2_5;
	failedfid_3_0_out = (Int_t)!passed_fid_3_0;
	failedfid_3_5_out = (Int_t)!passed_fid_3_5;
	hcalon_out = (Int_t)!hcalON;
	hcalon_bc_out = (Int_t)!hcalON_bc;
	charge_avg_out = charge;
	comp_event_frac_out = comp_ev_fraction;
	event_frac_out = ev_fraction;
	acc_charge_out = accumulated_charge;
	sbs_hcal_bclus_x_out = x_bestcluster;
	sbs_hcal_bclus_y_out = y_bestcluster;
	sbs_hcal_bclus_e_out = ce_bestcluster;
	sbs_hcal_bclus_row_out = row_bestcluster;
	sbs_hcal_bclus_col_out = col_bestcluster;

	//cout <<acc_charge_out<<"/"<<charge<< endl;


	//Fill old output tree
	bb_tr_vx_out = vx[0];
	bb_tr_vy_out = vy[0];
	bb_tr_vz_out = vz[0];
	bb_tr_p_out = p[0];
	bb_tr_px_out = px[0];
	bb_tr_py_out = py[0];
	bb_tr_pz_out = pz[0];
	bb_tr_n_out = ntrack;
	e_kine_Q2_out = ekineQ2;
	e_kine_W2_out = ekineW2;
	e_kine_nu_out = ekinenu;
	bb_ps_e_out = ePS;
	bb_ps_rowblk_out = rblkPS;
	bb_ps_colblk_out = cblkPS;
	bb_sh_e_out = eSH;
	bb_sh_rowblk_out = rblkSH;
	bb_sh_colblk_out = cblkSH;
	bb_sh_atimeblk_out = atimeSH;
	bb_sh_nclus_out = nclusSH;
	bb_hodotdc_clus_tmean_out = hodotmean[0];
	bb_gem_track_nhits_out = gemNhits;
	bb_etot_over_p_out = eop;
	sbs_hcal_e_out = hcale;
	sbs_hcal_x_out = hcalx;
	sbs_hcal_y_out = hcaly;
	sbs_hcal_rowblk_out = hcalr;
	sbs_hcal_colblk_out = hcalc;
	sbs_hcal_atimeblk_out = hcalatime;
	sbs_hcal_tdctimeblk_out = hcaltdc;
	sbs_hcal_nclus_out = nclus;
	Ndata_sbs_hcal_clus_id_out = Nhcalcid;
	for( Int_t c=0; c<Nhcalcid; ++c ){
	  Double_t cid = hcalcid[c];
	  Double_t ce = hcalce[c];
	  Double_t cx = hcalcx[c];
	  Double_t cy = hcalcy[c];
	  Double_t ctdc = hcalctdctime[c];
	  Double_t catime = hcalcatime[c];

	  sbs_hcal_clus_id_out[c] = cid;
	  sbs_hcal_clus_e_out[c] = ce;
	  sbs_hcal_clus_x_out[c] = cx;	  
	  sbs_hcal_clus_y_out[c] = cy;
	  sbs_hcal_clus_tdctime_out[c] = ctdc;
	  sbs_hcal_clus_atime_out[c] = catime;
	}
	sbs_hcal_nblk_out = nblk;
	Ndata_sbs_hcal_clus_blk_id_out = Nhcalcbid;
	for( Int_t b=0; b<Nhcalcbid; ++b ){
	  Double_t cbid = hcalcbid[b];
	  Double_t cbe = hcalcbe[b];
	  Double_t cbx = hcalcbx[b];
	  Double_t cby = hcalcby[b];
	  Double_t cbtdc = hcalcbtdctime[b];
	  Double_t cbatime = hcalcbatime[b];

	  sbs_hcal_clus_blk_id_out[b] = cbid;
	  sbs_hcal_clus_blk_e_out[b] = cbe;
	  sbs_hcal_clus_blk_x_out[b] = cbx;	  
	  sbs_hcal_clus_blk_y_out[b] = cby;
	  sbs_hcal_clus_blk_tdctime_out[b] = cbtdc;
	  sbs_hcal_clus_blk_atime_out[b] = cbatime;

	}

	P->Fill();

      }//end event loop

      cout << "Here 01" << endl;

      // getting ready for the next run
      C->Reset();

      cout << "Here 02" << endl;

    }//end run loop

  }//end target loop

  fout->Write();

  std::cout << std::endl << "Parsing complete. Output written to " << parse_path << std::endl << std::endl;

  st->Stop();

  // Send time efficiency report to console
  std::cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << std::endl;

}
