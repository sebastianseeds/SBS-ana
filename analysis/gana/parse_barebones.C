//sseeds 10.23.23 - Updated parsing script to cut both inelastic events and unused branches. Configured only to parse data files, not MC. Event parsing using wide globalcuts, wide W2 cuts, and wide coin (HCal/BBCal) cuts. Branch parsing includes only branches that sseeds is using for his gmn analysis
//Update 2.2.24 - Same method, made simpler without class references for troubleshooting
//Added multicluster analysis
//Update 5.13.24 - Added track variables to assess validity and set LD2 OR LH2

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
const Double_t intime_cut_override = 10.; //override the intime cut for kinematics with many fields settings
const Double_t nsig_step = 0.1; //number of sigma to step through in dx until fiducial cut failure written to output tree
const Double_t R_mclus = 0.6; //search region for clusters to add to highest energy cluster on multicluster analysis
const Double_t Nsig_fid_qual = 1.; //number of proton sigma to add to safety margin for fid cut quality plots

//Pion production cluster check
const double spatial_threshold = 60.0; // cm, informed by shower size max in MC
const double time_threshold = 10.0; // ns, wide, informed by shower formation limits (light propagation and decay time in WLS/scintillator HCal)

//Specific wide cut for all parsing. Account for possible survival of additional triggers in reconstructed data set.
const std::string gcut = "bb.ps.e>0.1&&abs(bb.tr.vz[0])<0.12&&fEvtHdr.fTrigBits==1";

//Forward declarations
void checkClusters(
		   const std::vector<double>& cluster_centroids_x, 
		   const std::vector<double>& cluster_centroids_y, 
		   const std::vector<double>& cluster_times
		   );

//MAIN. kine=kinematic, pass=reconstruction pass, cluster_method=best cluster selection method, coin_override=use intime_cut_override global instead of json, verbose=more error messages, lh2opt=use LH2 data, ld2opt=use LD2 data, effz=apply effective z offsets, ep_fourvec=reconstruct e' track using root fourvectors, debug=include debugging messages
void parse_barebones( Int_t kine = 7, 
		      Int_t pass = 2, 
		      Int_t cluster_method = 3,
		      bool coin_override = true,
		      bool verbose = false, 
		      bool lh2opt = true, 
		      bool ld2opt = false, 
		      bool effz = true, 
		      bool ep_fourvec = false, 
		      bool debug = false )
{   

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  //One set of data files in json config shared between pass0/1 per kinematic
  if( pass==0 )
    pass=1;

  //Use the intime algorithm if the intime cut is to be used
  if( coin_override )
    cluster_method=3;

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/parse.json");

  std::string rootfile_dir_lh2 = jmgr->GetValueFromSubKey_str( Form("rootfile_dir_lh2_p%d",pass), Form("sbs%d",kine) );
  std::string rootfile_dir_ld2 = jmgr->GetValueFromSubKey_str( Form("rootfile_dir_ld2_p%d",pass), Form("sbs%d",kine) );
  Double_t minE = jmgr->GetValueFromSubKey<Double_t>( Form("minE_p%d",pass), Form("sbs%d",kine) );
  Double_t hcal_v_offset = jmgr->GetValueFromSubKey<Double_t>( Form("hcal_offset_p%d",pass), Form("sbs%d",kine) );

  //get fit profiles
  vector<Double_t> coin_profile;
  jmgr->GetVectorFromSubKey<Double_t>(Form("coin_profile_p%d",pass),Form("sbs%d",kine),coin_profile);
  vector<Double_t> tof_profile;
  jmgr->GetVectorFromSubKey<Double_t>(Form("tof_profile_p%d",pass),Form("sbs%d",kine),tof_profile);
  vector<Double_t> tw_profile;
  jmgr->GetVectorFromSubKey<Double_t>(Form("tw_profile_p%d",pass),Form("sbs%d",kine),tw_profile);

  if( cluster_method==4 )
    cout << "Loaded coincidence time profile. Mean " << coin_profile[1] << ", sigma " << coin_profile[2] << ". Intime around " << coin_sigma_factor*coin_profile[2] << "." << endl;
  else if( cluster_method==3 )
    cout << "Using intime method for cluster selection." << endl;

  //Necessary for loop over directories in sbs8
  std::vector<TString> directories_h = {
    rootfile_dir_lh2 + "/SBS0percent",
    rootfile_dir_lh2 + "/SBS100percent",
    rootfile_dir_lh2 + "/SBS50percent",
    rootfile_dir_lh2 + "/SBS70percent_part1",
    rootfile_dir_lh2 + "/SBS70percent_part2",
    rootfile_dir_lh2 + "/SBS70percent_part3"
  };

  std::vector<TString> directories_d = {
    rootfile_dir_ld2 + "/SBS0percent",
    rootfile_dir_ld2 + "/SBS100percent",
    rootfile_dir_ld2 + "/SBS50percent",
    rootfile_dir_ld2 + "/SBS70percent_part1",
    rootfile_dir_ld2 + "/SBS70percent_part2",
    rootfile_dir_ld2 + "/SBS70percent_part3",
    rootfile_dir_ld2 + "/SBS70percent_part4"
  };

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
  std::string debug_word = "";
  if(debug)
    debug_word = "_debug";
  std::string effz_word = "";
  if(effz)
    effz_word = "_effz";
  std::string four_word = "";
  if(ep_fourvec)
    four_word = "_fourvec";
  std::string tar_word = "";
  if(lh2opt && !ld2opt)
    tar_word = "_lh2";
  else if(!lh2opt && ld2opt)
    tar_word = "_ld2";

  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string parse_path = outdir_path + Form("/parse/bunchedClusterTest_parse_sbs%d_pass%d_barebones%s%s%s%s.root",kine,pass,debug_word.c_str(),effz_word.c_str(),four_word.c_str(),tar_word.c_str());

  //set up output files
  TFile *fout = new TFile( parse_path.c_str(), "RECREATE" );

  //set up diagnostic histograms
  TH2D *hW2mag = new TH2D( "hW2mag", "W^{2} vs sbsmag; \%; GeV^{2}", 20, 0, 100, 200, 0, 2 );
  TH2D *hQ2mag = new TH2D( "hQ2mag", "Q^{2} vs sbsmag; \%; GeV^{2}", 20, 0, 100, 200, 0, 4 ); //can generalize
  TH2D *hdxmag = new TH2D( "hdxmag","dx vs sbsmag; \%; x_{HCAL}-x_{expect} (m)", 20, 0, 100, 800, -4, 4 );
  TH2D *hdymag = new TH2D( "hdymag","dy vs sbsmag; \%; y_{HCAL}-y_{expect} (m)", 20, 0, 100, 800, -4, 4 );
  TH1D *hcidxfail_d = new TH1D( "hcidxfail_d","Cluster Index Failed vs Run; runnum", ndruns, 0, ndruns );
  TH1D *hcidxfail_h = new TH1D( "hcidxfail_h","Cluster Index Failed vs Run; runnum", nhruns, 0, nhruns );

  TH2D *hexpxy_p = new TH2D( "hexpxy_p","proton hcal exp x vs hcal exp y; #sigma (m); x_{expect} (m)", 400, -2, 2, 600, -4, 2 );
  TH2D *hexpxy_p_fid = new TH2D( "hexpxy_p_fid",Form("proton hcal exp x vs hcal exp y %0.1f#sigma fid cut; #sigma (m); x_{expect} (m)",Nsig_fid_qual), 400, -2, 2, 600, -4, 2 );

  TH2D *hexpxy_n = new TH2D( "hexpxy_n","neutron hcal exp x vs hcal exp y; #sigma (m); x_{expect} (m)", 400, -2, 2, 600, -4, 2 );
  TH2D *hexpxy_n_fid = new TH2D( "hexpxy_n_fid",Form("neutron hcal exp x vs hcal exp y %0.1f#sigma fid cut; #sigma (m); x_{expect} (m)",Nsig_fid_qual), 400, -2, 2, 600, -4, 2 );

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
  Double_t W2_uc_out;
  Double_t Q2_out;
  Double_t Q2_uc_out;
  Double_t nu_out;
  Double_t precon_out;
  Double_t tau_out;
  Double_t epsilon_out;
  Double_t dttof_out;
  Double_t dttw_out;
  Double_t dttw_bc_out;
  Double_t dttw_b1_out;
  Double_t dttw_b2_out;
  Double_t dttw_b3_out;

  //Fiducial slices
  Double_t fiducial_sig_x_out;
  Double_t fiducial_sig_y_out;

  //Primary cluster
  Double_t dx_out;
  Double_t dy_out;
  Double_t coin_out;
  Double_t thetapq_pout;
  Double_t thetapq_nout;
  double sigsep_pout; //Total scale N, where N*sig_dx_p and N*sig_dy_p (and mean dx_p, dy_p) defines the minimum area necessary to include point on event (dx,dy)
  double sigsep_nout; //Total scale N, where N*sig_dx_n and N*sig_dy_n (and mean dx_n, dy_n) defines the minimum area necessary to include point on event (dx,dy)
  int hcalon_out;
  double accsep_xout; //Total scale N, where N*block_width_x necessary to reject event on active area cut
  double accsep_yout; //Total scale N, where N*block_width_y necessary to reject event on active area cut
  Double_t hcalnblk_out;
  Double_t hcalpid_out;
  Double_t hcalx_out;
  Double_t hcaly_out;
  Double_t hcale_out;  
  Double_t hcaltdc_out;  
  Double_t hcaltdc_b1_out;  
  Double_t hcaltdc_b2_out;  
  Double_t hcaltdc_b3_out;
  Double_t hcaltdctw_out;  
  Double_t hcaltdctof_out;  
  Double_t hcaltdccorr_out;  
  Double_t hcalatime_out;  
  Double_t hcalatime_b1_out;  
  Double_t hcalatime_b2_out;  
  Double_t hcalatime_b3_out; 
  Double_t hcal_index_out;
  Double_t hcalcbid_b1_out;  
  Double_t hcalcbid_b2_out;  
  Double_t hcalcbid_b3_out;  

  //Best cluster
  Double_t dx_bc_out;
  Double_t dy_bc_out;
  Double_t coin_bc_out;
  Double_t thetapq_bc_pout;
  Double_t thetapq_bc_nout;
  double sigsep_bc_pout; //Total scale N, where N*sig_dx_p and N*sig_dy_p (and mean dx_p, dy_p) defines the minimum area necessary to include point on event (dx,dy)
  double sigsep_bc_nout; //Total scale N, where N*sig_dx_n and N*sig_dy_n (and mean dx_n, dy_n) defines the minimum area necessary to include point on event (dx,dy)  int hcalon_bc_out;
  int hcalon_bc_out;
  double accsep_bc_xout; //Total scale N, where N*block_width_x necessary to reject event on active area cut
  double accsep_bc_yout; //Total scale N, where N*block_width_y necessary to reject event on active area cut
  Double_t hcalnblk_bc_out;
  Double_t hcalpid_bc_out;
  Double_t hcalx_bc_out;
  Double_t hcaly_bc_out;
  Double_t hcale_bc_out;  
  Double_t hcaltdc_bc_out;  
  Double_t hcaltdctw_bc_out;  
  Double_t hcaltdctof_bc_out;  
  Double_t hcaltdccorr_bc_out;  
  Double_t hcalatime_bc_out; 

  //Multi cluster
  double dx_mclus_out;
  double dy_mclus_out;
  double coin_mclus_out;
  double thetapq_mclus_pout;
  double thetapq_mclus_nout;
  double sigsep_mclus_pout; //Total scale N, where N*sig_dx_p and N*sig_dy_p (and mean dx_p, dy_p) defines the minimum area necessary to include point on event (dx,dy)
  double sigsep_mclus_nout; //Total scale N, where N*sig_dx_n and N*sig_dy_n (and mean dx_n, dy_n) defines the minimum area necessary to include point on event (dx,dy)  int hcalon_mclus_out;
  int hcalon_mclus_out;
  double accsep_mclus_xout; //Total scale N, where N*block_width_x necessary to reject event on active area cut
  double accsep_mclus_yout; //Total scale N, where N*block_width_y necessary to reject event on active area cut
  double hcalnclus_mclus_out;
  double hcalpid_mclus_out;
  double hcalx_mclus_out;
  double hcaly_mclus_out;
  double hcale_mclus_out;  

  //Additional cluster energy ranked (first)
  double dx_fc_out;
  double dy_fc_out;
  double coin_fc_out;
  double thetapq_fc_pout;
  double thetapq_fc_nout;
  double sigsep_fc_pout; //Total scale N, where N*sig_dx_p and N*sig_dy_p (and mean dx_p, dy_p) defines the minimum area necessary to include point on event (dx,dy)
  double sigsep_fc_nout; //Total scale N, where N*sig_dx_n and N*sig_dy_n (and mean dx_n, dy_n) defines the minimum area necessary to include point on event (dx,dy)  int hcalon_fc_out;
  int hcalon_fc_out;
  double accsep_fc_xout; //Total scale N, where N*block_width_x necessary to reject event on active area cut
  double accsep_fc_yout; //Total scale N, where N*block_width_y necessary to reject event on active area cut
  double hcalnblk_fc_out;
  double hcalpid_fc_out;
  double hcalx_fc_out;
  double hcaly_fc_out;
  double hcale_fc_out;

  //Additional cluster energy ranked (second)
  double dx_sc_out;
  double dy_sc_out;
  double coin_sc_out;
  double thetapq_sc_pout;
  double thetapq_sc_nout;
  double sigsep_sc_pout; //Total scale N, where N*sig_dx_p and N*sig_dy_p (and mean dx_p, dy_p) defines the minimum area necessary to include point on event (dx,dy)
  double sigsep_sc_nout; //Total scale N, where N*sig_dx_n and N*sig_dy_n (and mean dx_n, dy_n) defines the minimum area necessary to include point on event (dx,dy)  int hcalon_sc_out;
  int hcalon_sc_out;
  double accsep_sc_xout; //Total scale N, where N*block_width_x necessary to reject event on active area cut
  double accsep_sc_yout; //Total scale N, where N*block_width_y necessary to reject event on active area cut
  double hcalnblk_sc_out;
  double hcalpid_sc_out;
  double hcalx_sc_out;
  double hcaly_sc_out;
  double hcale_sc_out;

  //Additional cluster energy ranked (third)
  double dx_tc_out;
  double dy_tc_out;
  double coin_tc_out;
  double thetapq_tc_pout;
  double thetapq_tc_nout;
  double sigsep_tc_pout; //Total scale N, where N*sig_dx_p and N*sig_dy_p (and mean dx_p, dy_p) defines the minimum area necessary to include point on event (dx,dy)
  double sigsep_tc_nout; //Total scale N, where N*sig_dx_n and N*sig_dy_n (and mean dx_n, dy_n) defines the minimum area necessary to include point on event (dx,dy)  int hcalon_tc_out;
  int hcalon_tc_out;
  double accsep_tc_xout; //Total scale N, where N*block_width_x necessary to reject event on active area cut
  double accsep_tc_yout; //Total scale N, where N*block_width_y necessary to reject event on active area cut
  double hcalnblk_tc_out;
  double hcalpid_tc_out;
  double hcalx_tc_out;
  double hcaly_tc_out;
  double hcale_tc_out;

  // relevant old output tree vars
  Double_t hcal_nclus_out;
  Double_t bb_tr_vz_out;
  Double_t bb_tr_n_out;
  Double_t bb_tr_p_out;
  Double_t bb_tr_th_out;
  Double_t bb_tr_ph_out;
  Double_t bb_tr_r_th_out;
  Double_t bb_tr_r_x_out;
  Double_t bb_tr_r_ph_out;
  Double_t bb_tr_r_y_out;
  Double_t bb_tr_chi2_out;
  Double_t bb_ps_e_out;
  Double_t bb_ps_rowblk_out;
  Double_t bb_ps_colblk_out;
  Double_t bb_sh_e_out;
  Double_t bb_sh_rowblk_out;
  Double_t bb_sh_colblk_out;
  Double_t bb_sh_atime_out;
  Double_t bb_hodotdc_clus_tmean_out;
  Double_t bb_grinch_tdc_clus_size_out;
  Double_t bb_grinch_tdc_clus_trackindex_out;
  Double_t bb_grinch_tdc_clus_adc_out;
  Double_t bb_gem_track_nhits_out;
  Double_t bb_gem_track_ngoodhits_out;
  Double_t bb_gem_track_chi2ndf_out;
  Double_t bb_etot_over_p_out;

  P->Branch("run", &run_out, "run/I");
  P->Branch("tar", &tar_out, "tar/I");
  P->Branch("mag", &mag_out, "mag/I");
  P->Branch("event", &event_out, "event/I");
  P->Branch("trig", &trig_out, "trig/I");
  P->Branch("xexp", &xexp_out, "xexp/D");
  P->Branch("yexp", &yexp_out, "yexp/D");
  P->Branch("W2", &W2_out, "W2/D");
  P->Branch("W2_uc", &W2_uc_out, "W2_uc/D");
  P->Branch("Q2", &Q2_out, "Q2/D");
  P->Branch("Q2_uc", &Q2_uc_out, "Q2_uc/D");
  P->Branch("nu", &nu_out, "nu/D");
  P->Branch("tau", &tau_out, "tau/D");
  P->Branch("epsilon", &epsilon_out, "epsilon/D");
  P->Branch("precon", &precon_out, "precon/D");
  P->Branch("dttof", &dttof_out, "dttof/D");
  P->Branch("dttw", &dttw_out, "dttw/D");
  P->Branch("dttw_bc", &dttw_bc_out, "dttw_bc/D");
  P->Branch("dttw_b1", &dttw_b1_out, "dttw_b1/D");
  P->Branch("dttw_b2", &dttw_b2_out, "dttw_b2/D");
  P->Branch("dttw_b3", &dttw_b3_out, "dttw_b3/D");

  P->Branch("fiducial_sig_x", &fiducial_sig_x_out, "fiducial_sig_x/D");
  P->Branch("fiducial_sig_y", &fiducial_sig_y_out, "fiducial_sig_y/D");

  P->Branch("dx", &dx_out, "dx/D");
  P->Branch("dy", &dy_out, "dy/D");
  P->Branch("coin", &coin_out, "coin/D");
  P->Branch("thetapq_p", &thetapq_pout, "thetapq_p/D");
  P->Branch("thetapq_n", &thetapq_nout, "thetapq_n/D");
  P->Branch("sigsep_p", &sigsep_pout, "sigsep_p/D");
  P->Branch("sigsep_n", &sigsep_nout, "sigsep_n/D");
  P->Branch("hcalon", &hcalon_out, "hcalon/I");
  P->Branch("accsep_x", &accsep_xout, "accsep_x/D");
  P->Branch("accsep_y", &accsep_yout, "accsep_y/D");
  P->Branch("hcalnblk", &hcalnblk_out, "hcalnblk/D");
  P->Branch("hcalpid", &hcalpid_out, "hcalpid/D");
  P->Branch("hcalx", &hcalx_out, "hcalx/D");
  P->Branch("hcaly", &hcaly_out, "hcaly/D");
  P->Branch("hcale", &hcale_out, "hcale/D");
  P->Branch("hcaltdc", &hcaltdc_out, "hcaltdc/D");
  P->Branch("hcaltdc_b1", &hcaltdc_b1_out, "hcaltdc_b1/D");
  P->Branch("hcaltdc_b2", &hcaltdc_b2_out, "hcaltdc_b2/D");
  P->Branch("hcaltdc_b3", &hcaltdc_b3_out, "hcaltdc_b3/D");
  P->Branch("hcaltdctw", &hcaltdctw_out, "hcaltdctw/D");
  P->Branch("hcaltdctof", &hcaltdctof_out, "hcaltdctof/D");
  P->Branch("hcaltdccorr", &hcaltdccorr_out, "hcaltdccorr/D");
  P->Branch("hcalatime", &hcalatime_out, "hcalatime/D");
  P->Branch("hcalatime_b1", &hcalatime_b1_out, "hcalatime_b1/D");
  P->Branch("hcalatime_b2", &hcalatime_b2_out, "hcalatime_b2/D");
  P->Branch("hcalatime_b3", &hcalatime_b3_out, "hcalatime_b3/D");
  P->Branch("hcal_index", &hcal_index_out, "hcal_index/D");
  P->Branch("hcalcbid_b1", &hcalcbid_b1_out, "hcalcbid_b1/D");
  P->Branch("hcalcbid_b2", &hcalcbid_b2_out, "hcalcbid_b2/D");
  P->Branch("hcalcbid_b3", &hcalcbid_b3_out, "hcalcbid_b3/D");

  P->Branch("dx_bc", &dx_bc_out, "dx_bc/D");
  P->Branch("dy_bc", &dy_bc_out, "dy_bc/D");
  P->Branch("coin_bc", &coin_bc_out, "coin_bc/D");
  P->Branch("thetapq_bc_p", &thetapq_bc_pout, "thetapq_bc_p/D");
  P->Branch("thetapq_bc_n", &thetapq_bc_nout, "thetapq_bc_n/D");
  P->Branch("sigsep_bc_p", &sigsep_bc_pout, "sigsep_bc_p/D");
  P->Branch("sigsep_bc_n", &sigsep_bc_nout, "sigsep_bc_n/D");
  P->Branch("accsep_bc_x", &accsep_bc_xout, "accsep_bc_x/D");
  P->Branch("accsep_bc_y", &accsep_bc_yout, "accsep_bc_y/D");
  P->Branch("hcalon_bc", &hcalon_bc_out, "hcalon_bc/I");
  P->Branch("hcalnblk_bc", &hcalnblk_bc_out, "hcalnblk_bc/D");
  P->Branch("hcalpid_bc", &hcalpid_bc_out, "hcalpid_bc/D");
  P->Branch("hcalx_bc", &hcalx_bc_out, "hcalx_bc/D");
  P->Branch("hcaly_bc", &hcaly_bc_out, "hcaly_bc/D");
  P->Branch("hcale_bc", &hcale_bc_out, "hcale_bc/D");
  P->Branch("hcaltdc_bc", &hcaltdc_bc_out, "hcaltdc_bc/D");
  P->Branch("hcaltdctw_bc", &hcaltdctw_bc_out, "hcaltdctw_bc/D");
  P->Branch("hcaltdctof_bc", &hcaltdctof_bc_out, "hcaltdctof_bc/D");
  P->Branch("hcaltdccorr_bc", &hcaltdccorr_bc_out, "hcaltdccorr_bc/D");
  P->Branch("hcalatime_bc", &hcalatime_bc_out, "hcalatime_bc/D");

  P->Branch("dx_mclus", &dx_mclus_out, "dx_mclus/D");
  P->Branch("dy_mclus", &dy_mclus_out, "dy_mclus/D");
  P->Branch("coin_mclus", &coin_mclus_out, "coin_mclus/D");
  P->Branch("thetapq_mclus_p", &thetapq_mclus_pout, "thetapq_mclus_p/D");
  P->Branch("thetapq_mclus_n", &thetapq_mclus_nout, "thetapq_mclus_n/D");
  P->Branch("sigsep_mclus_p", &sigsep_mclus_pout, "sigsep_mclus_p/D");
  P->Branch("sigsep_mclus_n", &sigsep_mclus_nout, "sigsep_mclus_n/D");
  P->Branch("accsep_mclus_x", &accsep_mclus_xout, "accsep_mclus_x/D");
  P->Branch("accsep_mclus_y", &accsep_mclus_yout, "accsep_mclus_y/D");
  P->Branch("hcalon_mclus", &hcalon_mclus_out, "hcalon_mclus/I");
  P->Branch("hcalnclus_mclus", &hcalnclus_mclus_out, "hcalnclus_mclus/D");
  P->Branch("hcalpid_mclus", &hcalpid_mclus_out, "hcalpid_mclus/D");
  P->Branch("hcalx_mclus", &hcalx_mclus_out, "hcalx_mclus/D");
  P->Branch("hcaly_mclus", &hcaly_mclus_out, "hcaly_mclus/D");
  P->Branch("hcale_mclus", &hcale_mclus_out, "hcale_mclus/D");

  P->Branch("dx_fc", &dx_fc_out, "dx_fc/D");
  P->Branch("dy_fc", &dy_fc_out, "dy_fc/D");
  P->Branch("coin_fc", &coin_fc_out, "coin_fc/D");
  P->Branch("thetapq_fc_p", &thetapq_fc_pout, "thetapq_fc_p/D");
  P->Branch("thetapq_fc_n", &thetapq_fc_nout, "thetapq_fc_n/D");
  P->Branch("sigsep_fc_p", &sigsep_fc_pout, "sigsep_fc_p/D");
  P->Branch("sigsep_fc_n", &sigsep_fc_nout, "sigsep_fc_n/D");
  P->Branch("accsep_fc_x", &accsep_fc_xout, "accsep_fc_x/D");
  P->Branch("accsep_fc_y", &accsep_fc_yout, "accsep_fc_y/D");
  P->Branch("hcalon_fc", &hcalon_fc_out, "hcalon_fc/I");
  P->Branch("hcalnblk_fc", &hcalnblk_fc_out, "hcalnblk_fc/D");
  P->Branch("hcalpid_fc", &hcalpid_fc_out, "hcalpid_fc/D");
  P->Branch("hcalx_fc", &hcalx_fc_out, "hcalx_fc/D");
  P->Branch("hcaly_fc", &hcaly_fc_out, "hcaly_fc/D");
  P->Branch("hcale_fc", &hcale_fc_out, "hcale_fc/D");

  P->Branch("dx_sc", &dx_sc_out, "dx_sc/D");
  P->Branch("dy_sc", &dy_sc_out, "dy_sc/D");
  P->Branch("coin_sc", &coin_sc_out, "coin_sc/D");
  P->Branch("thetapq_sc_p", &thetapq_sc_pout, "thetapq_sc_p/D");
  P->Branch("thetapq_sc_n", &thetapq_sc_nout, "thetapq_sc_n/D");
  P->Branch("sigsep_sc_p", &sigsep_sc_pout, "sigsep_sc_p/D");
  P->Branch("sigsep_sc_n", &sigsep_sc_nout, "sigsep_sc_n/D");
  P->Branch("accsep_sc_x", &accsep_sc_xout, "accsep_sc_x/D");
  P->Branch("accsep_sc_y", &accsep_sc_yout, "accsep_sc_y/D");
  P->Branch("hcalon_sc", &hcalon_sc_out, "hcalon_sc/I");
  P->Branch("hcalnblk_sc", &hcalnblk_sc_out, "hcalnblk_sc/D");
  P->Branch("hcalpid_sc", &hcalpid_sc_out, "hcalpid_sc/D");
  P->Branch("hcalx_sc", &hcalx_sc_out, "hcalx_sc/D");
  P->Branch("hcaly_sc", &hcaly_sc_out, "hcaly_sc/D");
  P->Branch("hcale_sc", &hcale_sc_out, "hcale_sc/D");

  P->Branch("dx_tc", &dx_tc_out, "dx_tc/D");
  P->Branch("dy_tc", &dy_tc_out, "dy_tc/D");
  P->Branch("coin_tc", &coin_tc_out, "coin_tc/D");
  P->Branch("thetapq_tc_p", &thetapq_tc_pout, "thetapq_tc_p/D");
  P->Branch("thetapq_tc_n", &thetapq_tc_nout, "thetapq_tc_n/D");
  P->Branch("sigsep_tc_p", &sigsep_tc_pout, "sigsep_tc_p/D");
  P->Branch("sigsep_tc_n", &sigsep_tc_nout, "sigsep_tc_n/D");
  P->Branch("accsep_tc_x", &accsep_tc_xout, "accsep_tc_x/D");
  P->Branch("accsep_tc_y", &accsep_tc_yout, "accsep_tc_y/D");
  P->Branch("hcalon_tc", &hcalon_tc_out, "hcalon_tc/I");
  P->Branch("hcalnblk_tc", &hcalnblk_tc_out, "hcalnblk_tc/D");
  P->Branch("hcalpid_tc", &hcalpid_tc_out, "hcalpid_tc/D");
  P->Branch("hcalx_tc", &hcalx_tc_out, "hcalx_tc/D");
  P->Branch("hcaly_tc", &hcaly_tc_out, "hcaly_tc/D");
  P->Branch("hcale_tc", &hcale_tc_out, "hcale_tc/D");

  P->Branch("hcal_nclus", &hcal_nclus_out, "hcal_nclus/D");
  P->Branch("bb_tr_n", &bb_tr_n_out, "bb_tr_n/D");
  P->Branch("bb_tr_vz", &bb_tr_vz_out, "bb_tr_vz/D");
  P->Branch("bb_tr_p", &bb_tr_p_out, "bb_tr_p/D");
  P->Branch("bb_tr_th", &bb_tr_th_out, "bb_tr_th/D");
  P->Branch("bb_tr_ph", &bb_tr_ph_out, "bb_tr_ph/D");
  P->Branch("bb_tr_r_x", &bb_tr_r_x_out, "bb_tr_r_x/D");
  P->Branch("bb_tr_r_th", &bb_tr_r_th_out, "bb_tr_r_th/D");
  P->Branch("bb_tr_r_y", &bb_tr_r_y_out, "bb_tr_r_y/D");
  P->Branch("bb_tr_r_ph", &bb_tr_r_ph_out, "bb_tr_r_ph/D");
  P->Branch("bb_tr_chi2", &bb_tr_chi2_out, "bb_tr_chi2/D");
  P->Branch("bb_ps_e", &bb_ps_e_out, "bb_ps_e/D");
  P->Branch("bb_ps_rowblk", &bb_ps_rowblk_out, "bb_ps_rowblk/D");
  P->Branch("bb_ps_colblk", &bb_ps_colblk_out, "bb_ps_colblk/D");
  P->Branch("bb_sh_e", &bb_sh_e_out, "bb_sh_e/D");
  P->Branch("bb_sh_rowblk", &bb_sh_rowblk_out, "bb_sh_rowblk/D");
  P->Branch("bb_sh_colblk", &bb_sh_colblk_out, "bb_sh_colblk/D");
  P->Branch("bb_sh_atime", &bb_sh_atime_out, "bb_sh_atime/D");
  P->Branch("bb_hodotdc_clus_tmean", &bb_hodotdc_clus_tmean_out, "bb_hodotdc_clus_tmean/D");
  P->Branch("bb_grinch_tdc_clus_size", &bb_grinch_tdc_clus_size_out, "bb_grinch_tdc_clus_size/D");
  P->Branch("bb_grinch_tdc_clus_trackindex", &bb_grinch_tdc_clus_trackindex_out, "bb_grinch_tdc_clus_trackindex/D");
  P->Branch("bb_grinch_tdc_clus_adc", &bb_grinch_tdc_clus_adc_out, "bb_grinch_tdc_clus_adc/D");
  P->Branch("bb_gem_track_nhits", &bb_gem_track_nhits_out, "bb_gem_track_nhits/D");
  P->Branch("bb_gem_track_ngoodhits", &bb_gem_track_ngoodhits_out, "bb_gem_track_ngoodhits/D");
  P->Branch("bb_gem_track_chi2ndf", &bb_gem_track_chi2ndf_out, "bb_gem_track_chi2ndf/D");
  P->Branch("bb_etot_over_p", &bb_etot_over_p_out, "bb_etot_over_p/D");

  // setup reporting indices
  Int_t curmag = -1;
  std::string curtar = "";

  //Set up hcal active area with bounds that match database on pass
  vector<Double_t> hcalaa;
  if(pass<2)
    hcalaa = cut::hcalaa_data_alt(1,1);
  else
    hcalaa = cut::hcalaa_mc(1,1); //verified 2.10.24

  for ( Int_t t=0; t<2; t++ ){ //loop over targets
    // t==0, lh2; t==1, ld2
    Int_t nruns = 1;

    // only process target of interest
    if( !lh2opt && t==0 )
      continue;
    if( !ld2opt && t==1 )
      continue;

    if( t==0 ){
      nruns = nhruns;
      std::cout << std::endl << "Proceeding to LH2 data.." << std::endl << std::endl;
    }
    if( t==1 ){
      nruns = ndruns;
      std::cout << std::endl << "Proceeding to LD2 data.." << std::endl << std::endl;
    }

    int debug_total = 0;

    for (Int_t irun=0; irun<nruns; irun++) {
      
      cout << endl << "irun/nruns " << irun << "/" << nruns << endl << endl;

      if(debug){
	cout << "Debug mode enabled. Running roughly 10k events..." << endl;
	if( debug_total>10000 )
	  continue;
      }

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
	
	//search for sbs8 file in subdirectories
	std::string rfname_sbs8_h;
	if( kine==8 ){
	  TString pattern_h = Form("*%d*",crunh[irun].runnum);
	  TString foundDir_h = util::FindFileInDirectories(pattern_h, directories_h);
	  if (!foundDir_h.IsNull()) {
	    std::cout << "SBS 8 lh2 file found in: " << foundDir_h << std::endl;
	    rfname_sbs8_h = foundDir_h;
	  } else {
	    std::cout << "SBS 8 lh2 file not found in " << foundDir_h << std::endl;
	  }

	  rfname = rfname_sbs8_h + Form("/*%d*",crunh[irun].runnum);
	  
	}else
	  rfname = rootfile_dir_lh2 + Form("/*%d*",crunh[irun].runnum);
	charge = crunh[irun].charge;
      }
      if( t==1 ){
	//std::cout << crund;
	runnum = crund[irun].runnum;
	mag = crund[irun].sbsmag / 21;
	ebeam = crund[irun].ebeam;
	targ = crund[irun].target;

	//search for sbs8 file in subdirectories
	std::string rfname_sbs8_d;
	if( kine==8 ){
	  TString pattern_d = Form("*%d*",crund[irun].runnum);
	  TString foundDir_d = util::FindFileInDirectories(pattern_d, directories_d);
	  if (!foundDir_d.IsNull()) {
	    std::cout << "SBS 8 ld2 file found in: " << foundDir_d << std::endl;
	    rfname_sbs8_d = foundDir_d;
	  } else {
	    std::cout << "SBS 8 ld2 file not found in " << foundDir_d << std::endl;
	  }
	  rfname = rfname_sbs8_d + Form("/*%d*",crund[irun].runnum);
	}else
	  rfname = rootfile_dir_ld2 + Form("/*%d*",crund[irun].runnum);

	charge = crund[irun].charge;
      }

      //set up configuration and tune objects to load analysis parameters
      SBSconfig config(kine,mag);

      //Obtain configuration pars from config file
      Double_t hcaltheta = config.GetHCALtheta_rad();
      Double_t hcaldist;
      if(effz){
	hcaldist = config.GetHCALeffdist();
	cout << "Loading effective z offset " << hcaldist << "..." << endl;
      }else
	hcaldist = config.GetHCALdist();
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
      Double_t dx_del   = tune.Getdx_del();
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

      cout << "Loaded file for analysis: " << rfname << endl;

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
      Double_t r_x[econst::maxtrack],r_th[econst::maxtrack],r_y[econst::maxtrack],r_ph[econst::maxtrack],tr_chi2[econst::maxtrack];
      std::vector<std::string> trvar = {"n","p","px","py","pz","x","y","th","ph","vx","vy","vz","tg_x","tg_y","tg_th","tg_ph","r_x","r_th","r_y","r_ph","chi2"};
      std::vector<void*> trvarlink = {&ntrack,&p,&px,&py,&pz,&xtr,&ytr,&thtr,&phtr,&vx,&vy,&vz,&xtgt,&ytgt,&thtgt,&phtgt,&r_x,&r_th,&r_y,&r_ph,&tr_chi2};
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
      Double_t gemNhits, gemNgoodhits, gemChiSqr, grinchClusSize, grinchClusTrIndex, grinchClusADC, eop;
      std::vector<std::string> miscbbvar = {"gem.track.nhits","gem.track.ngoodhits","gem.track.chi2ndf","grinch_tdc.clus.size","grinch_tdc.clus.trackindex","grinch_tdc.clus.adc","etot_over_p"};
      std::vector<void*> miscbbvarlink = {&gemNhits,&gemNgoodhits,&gemChiSqr,&grinchClusSize,&grinchClusTrIndex,&grinchClusADC,&eop};
      rvars::setbranch(C, "bb", miscbbvar, miscbbvarlink);

      TCut GCut = gcut.c_str();

      TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );

      // get experimental quantities by run
      //set up hcal coordinate system with hcal angle wrt exit beamline
      vector<TVector3> hcalaxes; vars::sethcalaxes( hcaltheta, hcalaxes );
      TVector3 hcalorigin = hcaldist*hcalaxes[2] + hcal_v_offset*hcalaxes[0];

      Double_t BdL = econst::sbsmaxfield * econst::sbsdipolegap * (mag/100); //scale crudely by magnetic field
      Double_t Eloss_outgoing;
      if(t==0)
	Eloss_outgoing = econst::celldiameter/2.0/sin(bbthr) * econst::lh2tarrho * econst::lh2dEdx;
      if(t==1)
	Eloss_outgoing = econst::celldiameter/2.0/sin(bbthr) * econst::ld2tarrho * econst::ld2dEdx;

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
      
      //Event loop
      while (C->GetEntry(nevent++)) {
	
	std::cout << "Processing run " << runnum << " event " << nevent << " / " << nevents << ", total passed cuts " << npassed << "\r";
	std::cout.flush();

	if((Int_t)hcalidx>9){
	  if(t==0)
	    hcidxfail_h->Fill(irun);
	  if(t==1)
	    hcidxfail_d->Fill(irun);
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
	ebeam_c = vars::ebeam_c( ebeam, vz[0], targ ); //correct for energy loss of beam electron as it passes into the target material towards the vertex position.

	TVector3 vertex( 0., 0., vz[0] );

	//reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
	Double_t precon = p[0] + Eloss_outgoing; //correct the momentum of e' after the vertex as it passes out of the target material.

	//set up four-momenta with some empty for various calculation methods
	TLorentzVector pbeam( 0., 0., ebeam_c, ebeam_c ); //beam momentum
	//TLorentzVector pe( px[0], py[0], pz[0], p[0] ); //e' momentum
	TLorentzVector pe( precon*px[0]/p[0], precon*py[0]/p[0], precon*pz[0]/p[0], precon ); //e' recon plvect

	//TLorentzVector pe( px[0], py[0], pz[0], p[0] ); //e' recon plvect
	TLorentzVector ptarg; vars::setPN(nucleon,ptarg); //target momentum
	TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
	TLorentzVector pN; //N' momentum
	TVector3 pNhat; //Unit N' 3-vector
      
	//simple calculations for e' and N'
	Double_t etheta = vars::etheta(pe); 
	Double_t ephi = vars::ephi(pe);
	Double_t pcent = vars::pcentral(ebeam,etheta,nucleon); //e' p reconstructed by angles
	Double_t phNexp = ephi + physconst::pi;
	Double_t Q2, Q2_uc, W2, W2_uc, nu, thNexp, pNexp, ebeam_o, tau, epsilon;
	ebeam_o = vars::ebeam_o( ebeam_c, etheta, targ ); //Second energy correction accounting for energy loss leaving target

	//Use reconstructed angles as independent qty (usually preferable given GEM precision at most kinematics)
	nu = pbeam.E() - pcent;
	pNexp = vars::pN_expect( nu, nucleon );
	//thNexp = acos((ebeam-pz[0])/pNexp);
	thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	pNhat = vars::pNhat_track( thNexp, phNexp );
	pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	Q2 = vars::Q2( pbeam.E(), pe.E(), etheta );
	Q2_uc = ekineQ2;
	W2 = vars::W2( pbeam.E(), pe.E(), Q2, nucleon );
	W2_uc = ekineW2;
	tau = vars::tau( Q2, nucleon );
	epsilon = vars::epsilon( tau, etheta );

	//Use four-momentum member functions for higher Q2 points (7,11)
	if(ep_fourvec){
	  pN = q + ptarg;
	  pNhat = pN.Vect().Unit();
	  Q2 = -q.M2();
	  W2 = pN.M2();
	  nu = q.E();
	}

	/////////////////////
	//W2 elastic cut
	bool failedW2 = W2>W2max;
	if(failedW2)
	  continue;

	npassed++;
	debug_total++;

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

	//Determine separation from data proton and neutron in dx,dy
	double separation_p = util::NspotScaleFactor(dy, dx, dy0, dx0_p, dysig, dxsig_p, 0);
	double separation_n = util::NspotScaleFactor(dy, dx, dy0, dx0_n, dysig, dxsig_n, 0);

	//////////////////////
	//ALL CLUSTER ANALYSIS
      
	//Set up clone clusters for selection analysis.
	vector<double> clone_cluster_intime;
	vector<double> clone_cluster_score;

	vector<double> cluster_centroids_x;
	vector<double> cluster_centroids_y;
	vector<double> cluster_times;
	vector<double> cluster_energies;

	// Set up the structure to hold energy and index, then sort based on energy
	std::vector<std::pair<double, int>> energy_index_pairs;

	//Add multicluster analysis for higher Q2 where hcal clusters may be separated
	vector<int> intime_cluster_indices;

	//loop through all clusters and select without HCal position information
	for( int c=0; c<Nhcalcid; c++ ){
	
	  //add to all cluster vectors
	  cluster_centroids_x.push_back(hcalcx[c]);
	  cluster_centroids_y.push_back(hcalcy[c]);
	  cluster_times.push_back(hcalcatime[c]);
	  cluster_energies.push_back(hcalce[c]);

	  //calculate h-arm physics quantities per cluster
	  double atime = hcalcatime[c];
	  double atime_diff = atime - atimeSH; //Assuming best shower time on primary cluster
	  double ce = hcalce[c];

	  //using hcal atime until after pass2, wide cut around sigma factor
	  bool passedCoin;
	    
	  if( coin_override )
	    passedCoin = abs(atime_diff)<intime_cut_override;
	  else
	    passedCoin = abs(atime_diff-coin_profile[1])<coin_sigma_factor*coin_profile[2];
	  
	  //Replicate the in-time algorithm with new cluster to be sorted later. All cluster elements included
	  clone_cluster_intime.push_back(ce);
	  if( !passedCoin ){
	    clone_cluster_intime[c] = 0;
	  }else{
	    intime_cluster_indices.push_back(c);
	  }

	  //Get score (no position info). Will be sorted later
	  double cascore = util::assignScore( ce, atime_diff, hcalce[(Int_t)hcalidx], coin_profile );
	  clone_cluster_score.push_back(cascore);
	  
	  //Add energy cluster index and energy to pair
	  energy_index_pairs.push_back(std::make_pair(ce, c));

	}//endloop over cluster elements

	// Check for the existence of associated clusters
	checkClusters(cluster_centroids_x, cluster_centroids_y, cluster_times);

	// Sort in descending order based on the cluster energy
	std::sort(energy_index_pairs.begin(), energy_index_pairs.end(), [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
	    return a.first > b.first; // Using the first element of pair (energy) for comparison
	  });

	///////////////////////
	//Multicluster analysis

	// Step 1: Find the highest energy cluster from intime_cluster_indices
	int highest_energy_idx = -1;
	double highest_energy = -1;
	for (int idx : intime_cluster_indices) {
	  if (hcalce[idx] > highest_energy) {
	    highest_energy = hcalce[idx];
	    highest_energy_idx = idx;
	  }
	}

	// Variables for combined cluster properties
	double mclus_e = highest_energy; // Start with highest energy cluster
	double mclus_x = hcalcx[highest_energy_idx] * highest_energy; // Weighted by energy
	double mclus_y = hcalcy[highest_energy_idx] * highest_energy; // Weighted by energy
	double mclus_total_energy = highest_energy; // Total energy for weighted position calculation
	double nclus_multicluster = 0;

	if(Nhcalcid>0)
	  nclus_multicluster++;

	// Step 2: Identify nearby clusters and combine
	for (int idx : intime_cluster_indices) {
	  if (idx == highest_energy_idx) continue; // Skip the highest energy cluster itself

	  double difx = hcalcx[idx] - hcalcx[highest_energy_idx];
	  double dify = hcalcy[idx] - hcalcy[highest_energy_idx];
	  double distance = sqrt(difx*difx + dify*dify);

	  if (distance <= R_mclus) {
	    // This cluster is within the radius, combine it
	    nclus_multicluster++;
	    mclus_e += hcalce[idx]; // Add the energy of the nearby cluster
	    mclus_x += hcalcx[idx] * hcalce[idx]; // Add weighted position
	    mclus_y += hcalcy[idx] * hcalce[idx]; // Add weighted position
	    mclus_total_energy += hcalce[idx]; // Increase total energy for weighting
	  }
	}

	// Step 3: Recompute the combined position based on the total combined energy
	if (mclus_total_energy > 0) { // To avoid division by zero
	  mclus_x /= mclus_total_energy;
	  mclus_y /= mclus_total_energy;
	}

	// if( nclus_multicluster>1 )
	//   cout << "  Multiple correlated clusters = " << nclus_multicluster << " with variables: energy = " << mclus_e << ", x " << mclus_x << ", y " << mclus_y << endl;

	/////////////////////////

	// Get cluster indices for top three energy clusters. Default to 0 if they don't exist.
	int cidxe_fc = 0;
	int cidxe_sc = 0;
	int cidxe_tc = 0;
	if(Nhcalcid>0)
	  cidxe_fc = energy_index_pairs[0].second;
	if(Nhcalcid>1)
	  cidxe_sc = energy_index_pairs[1].second;
	if(Nhcalcid>2)
	  cidxe_tc = energy_index_pairs[2].second;
      
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
	    intime_idx = c; //The dimension of clone_cluster_intime is always the whole cluster
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
	  cidx_best = cidxe_fc;
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

	//Calculations from the multi cluster
	TVector3 hcalpos_multicluster = hcalorigin + 
	  mclus_x*hcalaxes[0] + 
	  mclus_y*hcalaxes[1];

	TVector3 protdir_mclus = ( hcalpos_multicluster + protdeflect * hcalaxes[0] - vertex ).Unit();
	TVector3 neutdir_mclus = ( hcalpos_multicluster - vertex ).Unit();
	double thetapq_mclus_p = acos( protdir_mclus.Dot( pNhat ) );
	double thetapq_mclus_n = acos( neutdir_mclus.Dot( pNhat ) );

	//Calculations from the first cluster
	TVector3 hcalpos_firstcluster = hcalorigin + 
	  hcalcx[cidxe_fc]*hcalaxes[0] + 
	  hcalcy[cidxe_fc]*hcalaxes[1];

	TVector3 protdir_fc = ( hcalpos_firstcluster + protdeflect * hcalaxes[0] - vertex ).Unit();
	TVector3 neutdir_fc = ( hcalpos_firstcluster - vertex ).Unit();
	double thetapq_fc_p = acos( protdir_fc.Dot( pNhat ) );
	double thetapq_fc_n = acos( neutdir_fc.Dot( pNhat ) );

	//Calculations from the second cluster
	TVector3 hcalpos_secondcluster = hcalorigin + 
	  hcalcx[cidxe_sc]*hcalaxes[0] + 
	  hcalcy[cidxe_sc]*hcalaxes[1];

	TVector3 protdir_sc = ( hcalpos_secondcluster + protdeflect * hcalaxes[0] - vertex ).Unit();
	TVector3 neutdir_sc = ( hcalpos_secondcluster - vertex ).Unit();
	double thetapq_sc_p = acos( protdir_sc.Dot( pNhat ) );
	double thetapq_sc_n = acos( neutdir_sc.Dot( pNhat ) );

	//Calculations from the third cluster
	TVector3 hcalpos_thirdcluster = hcalorigin + 
	  hcalcx[cidxe_tc]*hcalaxes[0] + 
	  hcalcy[cidxe_tc]*hcalaxes[1];

	TVector3 protdir_tc = ( hcalpos_thirdcluster + protdeflect * hcalaxes[0] - vertex ).Unit();
	TVector3 neutdir_tc = ( hcalpos_thirdcluster - vertex ).Unit();
	double thetapq_tc_p = acos( protdir_tc.Dot( pNhat ) );
	double thetapq_tc_n = acos( neutdir_tc.Dot( pNhat ) );

	//best cluster variables
	double dx_bestcluster = hcalcx[cidx_best] - xyhcalexp[0];
	double dy_bestcluster = hcalcy[cidx_best] - xyhcalexp[1];
	double hatime_bestcluster = hcalcatime[cidx_best];
	double hcoin_bestcluster = hcalcatime[cidx_best] - atimeSH;
	double ce_bestcluster = hcalce[cidx_best];
	double x_bestcluster = hcalcx[cidx_best];
	double y_bestcluster = hcalcy[cidx_best];
	double row_bestcluster = hcalcrow[cidx_best];
	double col_bestcluster = hcalccol[cidx_best];
	double nblk_bestcluster = hcalcnblk[cidx_best];
	double pblkid_bestcluster = hcalcid[cidx_best];

	//multi cluster variables
	double dx_multicluster = mclus_x - xyhcalexp[0];
	double dy_multicluster = mclus_y - xyhcalexp[1];
	double hatime_multicluster = hcalcatime[cidxe_fc]; //Use adc time from highest e clus
	double hcoin_multicluster = hcalcatime[cidxe_fc] - atimeSH;
	double ce_multicluster = mclus_e;
	double x_multicluster = mclus_x;
	double y_multicluster = mclus_y;

	//first cluster variables
	double dx_firstcluster = hcalcx[cidxe_fc] - xyhcalexp[0];
	double dy_firstcluster = hcalcy[cidxe_fc] - xyhcalexp[1];
	double hatime_firstcluster = hcalcatime[cidxe_fc];
	double hcoin_firstcluster = hcalcatime[cidxe_fc] - atimeSH;
	double ce_firstcluster = hcalce[cidxe_fc];
	double x_firstcluster = hcalcx[cidxe_fc];
	double y_firstcluster = hcalcy[cidxe_fc];
	double row_firstcluster = hcalcrow[cidxe_fc];
	double col_firstcluster = hcalccol[cidxe_fc];
	double nblk_firstcluster = hcalcnblk[cidxe_fc];
	double pblkid_firstcluster = hcalcid[cidxe_fc];

	//second cluster variables
	double dx_secondcluster = hcalcx[cidxe_sc] - xyhcalexp[0];
	double dy_secondcluster = hcalcy[cidxe_sc] - xyhcalexp[1];
	double hatime_secondcluster = hcalcatime[cidxe_sc];
	double hcoin_secondcluster = hcalcatime[cidxe_sc] - atimeSH;
	double ce_secondcluster = hcalce[cidxe_sc];
	double x_secondcluster = hcalcx[cidxe_sc];
	double y_secondcluster = hcalcy[cidxe_sc];
	double row_secondcluster = hcalcrow[cidxe_sc];
	double col_secondcluster = hcalccol[cidxe_sc];
	double nblk_secondcluster = hcalcnblk[cidxe_sc];
	double pblkid_secondcluster = hcalcid[cidxe_sc];

	//third cluster variables
	double dx_thirdcluster = hcalcx[cidxe_tc] - xyhcalexp[0];
	double dy_thirdcluster = hcalcy[cidxe_tc] - xyhcalexp[1];
	double hatime_thirdcluster = hcalcatime[cidxe_tc];
	double hcoin_thirdcluster = hcalcatime[cidxe_tc] - atimeSH;
	double ce_thirdcluster = hcalce[cidxe_tc];
	double x_thirdcluster = hcalcx[cidxe_tc];
	double y_thirdcluster = hcalcy[cidxe_tc];
	double row_thirdcluster = hcalcrow[cidxe_tc];
	double col_thirdcluster = hcalccol[cidxe_tc];
	double nblk_thirdcluster = hcalcnblk[cidxe_tc];
	double pblkid_thirdcluster = hcalcid[cidxe_tc];

	//Determine separation from data proton and neutron in dx,dy (best cluster)
	double separation_bc_p = util::NspotScaleFactor(dy_bestcluster, dx_bestcluster, dy0, dx0_p, dysig, dxsig_p, 0);
	double separation_bc_n = util::NspotScaleFactor(dy_bestcluster, dx_bestcluster, dy0, dx0_n, dysig, dxsig_n, 0);

	//Determine separation from data proton and neutron in dx,dy (multi cluster)
	double separation_mclus_p = util::NspotScaleFactor(dy_multicluster, dx_multicluster, dy0, dx0_p, dysig, dxsig_p, 0);
	double separation_mclus_n = util::NspotScaleFactor(dy_multicluster, dx_multicluster, dy0, dx0_n, dysig, dxsig_n, 0);

	//Determine separation from data proton and neutron in dx,dy (first cluster)
	double separation_fc_p = util::NspotScaleFactor(dy_firstcluster, dx_firstcluster, dy0, dx0_p, dysig, dxsig_p, 0);
	double separation_fc_n = util::NspotScaleFactor(dy_firstcluster, dx_firstcluster, dy0, dx0_n, dysig, dxsig_n, 0);

	//Determine separation from data proton and neutron in dx,dy (second cluster)
	double separation_sc_p = util::NspotScaleFactor(dy_secondcluster, dx_secondcluster, dy0, dx0_p, dysig, dxsig_p, 0);
	double separation_sc_n = util::NspotScaleFactor(dy_secondcluster, dx_secondcluster, dy0, dx0_n, dysig, dxsig_n, 0);

	//Determine separation from data proton and neutron in dx,dy (third cluster)
	double separation_tc_p = util::NspotScaleFactor(dy_thirdcluster, dx_thirdcluster, dy0, dx0_p, dysig, dxsig_p, 0);
	double separation_tc_n = util::NspotScaleFactor(dy_thirdcluster, dx_thirdcluster, dy0, dx0_n, dysig, dxsig_n, 0);

	//Find fiducial cut sigma factor
	std::pair<Double_t, Double_t> fiducial_factors = cut::findFidFailure(dxsig_p, 
									     dysig, 
									     xyhcalexp[0],
									     xyhcalexp[1], 
									     dx_del,
									     hcalaa);

	//Fill fiducial quality histograms
	hexpxy_p->Fill(xyhcalexp[1],xyhcalexp[0]-dx_del);
	hexpxy_n->Fill(xyhcalexp[1],xyhcalexp[0]);
	
	if(fiducial_factors.first>Nsig_fid_qual && fiducial_factors.second>Nsig_fid_qual){
	  hexpxy_p_fid->Fill(xyhcalexp[1],xyhcalexp[0]-dx_del);
	  hexpxy_n_fid->Fill(xyhcalexp[1],xyhcalexp[0]);
	}

	//H-arm active area cut (primary cluster)
	bool hcalON = cut::hcalaaON(hcalx,hcaly,hcalaa);

	//H-arm active area cuts (best_cluster)
	bool hcalON_bc = cut::hcalaaON(hcalcx[cidx_best],hcalcy[cidx_best],hcalaa);
	bool hcalON_mclus = cut::hcalaaON(mclus_x,mclus_y,hcalaa);
	bool hcalON_fc = cut::hcalaaON(hcalcx[cidxe_fc],hcalcy[cidxe_fc],hcalaa);
	bool hcalON_sc = cut::hcalaaON(hcalcx[cidxe_sc],hcalcy[cidxe_sc],hcalaa);
	bool hcalON_tc = cut::hcalaaON(hcalcx[cidxe_tc],hcalcy[cidxe_tc],hcalaa);

	//Get min active area cut in x and y to reject event (in terms of block width x and y)
	std::pair<double,double> minaacut = util::minaa(hcalx,hcaly);
	std::pair<double,double> minaacut_bc = util::minaa(hcalcx[cidx_best],hcalcy[cidx_best]);
	std::pair<double,double> minaacut_mclus = util::minaa(mclus_x,mclus_y);
	std::pair<double,double> minaacut_fc = util::minaa(hcalcx[cidxe_fc],hcalcy[cidxe_fc]);
	std::pair<double,double> minaacut_sc = util::minaa(hcalcx[cidxe_sc],hcalcy[cidxe_sc]);
	std::pair<double,double> minaacut_tc = util::minaa(hcalcx[cidxe_tc],hcalcy[cidxe_tc]);

	//Get timing variables
	double hcalTDC = hcaltdc;
	double hcalTDC_bc = hcalctdctime[cidx_best];
	double hcalADCt = hcalatime;
	double hcalADCt_bc = hcalcatime[cidx_best];
	
	//Get timing variables with nucleon momentum ToF and timewalk corrections
	double dt_tw = tw_profile[0]/pow(hcale,tw_profile[1]);
	double dt_tw_bc = tw_profile[0]/pow(ce_bestcluster,tw_profile[1]);
	double dt_tof = tof_profile[0]*pNexp + tof_profile[1]*pow(pNexp,2) + tof_profile[2]*pow(pNexp,3);
	double dt_tw_b1 = tw_profile[0]/pow(hcalcbe[0],tw_profile[1]);
	double dt_tw_b2 = tw_profile[0]/pow(hcalcbe[1],tw_profile[1]);
	double dt_tw_b3 = tw_profile[0]/pow(hcalcbe[2],tw_profile[1]);

	double hcalTDCtof = hcaltdc - dt_tof;
	double hcalTDCtof_bc = hcalctdctime[cidx_best] - dt_tof;
	double hcalTDCtw = hcaltdc - dt_tw;
	double hcalTDCtw_bc = hcalctdctime[cidx_best] - dt_tw;

	// Apply tof and tw corrections
	double hcalTDC_corrected = hcaltdc - dt_tof - dt_tw;
	double hcalTDC_corrected_bc = hcalctdctime[cidx_best] - dt_tof - dt_tw_bc;

	// Get internal resolution
	double hcalTDC_block1 = hcalcbtdctime[0];
	double hcalTDC_block2 = hcalcbtdctime[1];
	double hcalTDC_block3 = hcalcbtdctime[2];

	double hcalADCt_block1 = hcalcbatime[0];
	double hcalADCt_block2 = hcalcbatime[1];
	double hcalADCt_block3 = hcalcbatime[2];

	double hcal_blockid1 = hcalcbid[0];
	double hcal_blockid2 = hcalcbid[1];
	double hcal_blockid3 = hcalcbid[2];

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
	W2_uc_out = W2_uc;
	Q2_out = Q2;
	Q2_uc_out = Q2_uc;
	nu_out = nu;
	tau_out = tau;
	epsilon_out = epsilon;
	precon_out = precon;
	dttof_out = dt_tof;
	dttw_out = dt_tw;
	dttw_bc_out = dt_tw_bc;
	dttw_b1_out = dt_tw_b1;
	dttw_b2_out = dt_tw_b2;
	dttw_b3_out = dt_tw_b3;

	//Fiducial slices
	fiducial_sig_x_out = fiducial_factors.first;
	fiducial_sig_y_out = fiducial_factors.second;

	//Primary Cluster
	dx_out = dx;
	dy_out = dy;
	coin_out = pclus_diff;
	thetapq_pout = thetapq_p;
	thetapq_nout = thetapq_n;
	sigsep_pout = separation_p;
	sigsep_nout = separation_n;
	hcalnblk_out = pclus_nblk;
	hcalon_out = hcalON;
	accsep_xout = minaacut.first;
	accsep_yout = minaacut.second;
	hcalpid_out = pblkid;
	hcalx_out = hcalx;
	hcaly_out = hcaly;
	hcale_out = hcale;
	hcaltdc_out = hcalTDC;
	hcaltdc_b1_out = hcalTDC_block1;
	hcaltdc_b2_out = hcalTDC_block2;
	hcaltdc_b3_out = hcalTDC_block3;
	hcaltdctw_out = hcalTDCtw;
	hcaltdctof_out = hcalTDCtof;
	hcaltdccorr_out = hcalTDC_corrected;
	hcalatime_out = hcalADCt;
	hcalatime_b1_out = hcalADCt_block1;
	hcalatime_b2_out = hcalADCt_block2;
	hcalatime_b3_out = hcalADCt_block3;
	hcal_index_out = hcalidx;
	hcalcbid_b1_out = hcal_blockid1;
	hcalcbid_b2_out = hcal_blockid2;
	hcalcbid_b3_out = hcal_blockid3;

	//Best Cluster
	dx_bc_out = dx_bestcluster;
	dy_bc_out = dy_bestcluster;
	coin_bc_out = hcoin_bestcluster;
	thetapq_bc_pout = thetapq_bc_p;
	thetapq_bc_nout = thetapq_bc_n;
	sigsep_bc_pout = separation_bc_p;
	sigsep_bc_nout = separation_bc_n;
	hcalnblk_bc_out = nblk_bestcluster;
	hcalon_bc_out = hcalON_bc;
	accsep_bc_xout = minaacut_bc.first;
	accsep_bc_yout = minaacut_bc.second;
	hcalpid_bc_out = pblkid_bestcluster;
	hcalx_bc_out = x_bestcluster;
	hcaly_bc_out = y_bestcluster;
	hcale_bc_out = ce_bestcluster;
	hcaltdc_bc_out = hcalTDC_bc;
	hcaltdctw_bc_out = hcalTDCtw_bc;
	hcaltdctof_bc_out = hcalTDCtof_bc;
	hcaltdccorr_bc_out = hcalTDC_corrected_bc;
	hcalatime_bc_out = hcalADCt_bc;

	//Multi Cluster
	dx_mclus_out = dx_multicluster;
	dy_mclus_out = dy_multicluster;
	coin_mclus_out = hcoin_multicluster;
	thetapq_mclus_pout = thetapq_mclus_p;
	thetapq_mclus_nout = thetapq_mclus_n;
	sigsep_mclus_pout = separation_mclus_p;
	sigsep_mclus_nout = separation_mclus_n;
	hcalon_mclus_out = hcalON_bc;
	accsep_mclus_xout = minaacut_bc.first;
	accsep_mclus_yout = minaacut_bc.second;
	hcalnclus_mclus_out = nclus_multicluster;
	hcalpid_mclus_out = pblkid_firstcluster; //use highest e cluster highest e blk id
	hcalx_mclus_out = x_multicluster;
	hcaly_mclus_out = y_multicluster;
	hcale_mclus_out = ce_multicluster;

	//First Cluster
	dx_fc_out = dx_firstcluster;
	dy_fc_out = dy_firstcluster;
	coin_fc_out = hcoin_firstcluster;
	thetapq_fc_pout = thetapq_fc_p;
	thetapq_fc_nout = thetapq_fc_n;
	sigsep_fc_pout = separation_fc_p;
	sigsep_fc_nout = separation_fc_n;
	hcalnblk_fc_out = nblk_firstcluster;
	hcalon_fc_out = hcalON_fc;
	accsep_fc_xout = minaacut_fc.first;
	accsep_fc_yout = minaacut_fc.second;
	hcalpid_fc_out = pblkid_firstcluster;
	hcalx_fc_out = x_firstcluster;
	hcaly_fc_out = y_firstcluster;
	hcale_fc_out = ce_firstcluster;

	//Second Cluster
	dx_sc_out = dx_secondcluster;
	dy_sc_out = dy_secondcluster;
	coin_sc_out = hcoin_secondcluster;
	thetapq_sc_pout = thetapq_sc_p;
	thetapq_sc_nout = thetapq_sc_n;
	sigsep_sc_pout = separation_sc_p;
	sigsep_sc_nout = separation_sc_n;
	hcalnblk_sc_out = nblk_secondcluster;
	hcalon_sc_out = hcalON_bc;
	accsep_sc_xout = minaacut_bc.first;
	accsep_sc_yout = minaacut_bc.second;
	hcalpid_sc_out = pblkid_secondcluster;
	hcalx_sc_out = x_secondcluster;
	hcaly_sc_out = y_secondcluster;
	hcale_sc_out = ce_secondcluster;

	//Third Cluster
	dx_tc_out = dx_thirdcluster;
	dy_tc_out = dy_thirdcluster;
	coin_tc_out = hcoin_thirdcluster;
	thetapq_tc_pout = thetapq_tc_p;
	thetapq_tc_nout = thetapq_tc_n;
	sigsep_tc_pout = separation_tc_p;
	sigsep_tc_nout = separation_tc_n;
	hcalnblk_tc_out = nblk_thirdcluster;
	hcalon_tc_out = hcalON_bc;
	accsep_tc_xout = minaacut_bc.first;
	accsep_tc_yout = minaacut_bc.second;
	hcalpid_tc_out = pblkid_thirdcluster;
	hcalx_tc_out = x_thirdcluster;
	hcaly_tc_out = y_thirdcluster;
	hcale_tc_out = ce_thirdcluster;

	//Fill old output tree
	hcal_nclus_out = nclus;
	bb_tr_p_out = p[0];
	bb_tr_n_out = ntrack;
	bb_tr_vz_out = vz[0];
	bb_tr_th_out = thtr[0];
	bb_tr_ph_out = phtr[0];
	bb_tr_r_x_out = r_x[0];
	bb_tr_r_th_out = r_th[0];
	bb_tr_r_y_out = r_y[0];
	bb_tr_r_ph_out = r_ph[0];
	bb_tr_chi2_out = tr_chi2[0];
	bb_ps_e_out = ePS;
	bb_ps_rowblk_out = rblkPS;
	bb_ps_colblk_out = cblkPS;
	bb_sh_e_out = eSH;
	bb_sh_rowblk_out = rblkSH;
	bb_sh_colblk_out = cblkSH;
	bb_sh_atime_out = atimeSH;
	bb_hodotdc_clus_tmean_out = hodotmean[0];
	bb_grinch_tdc_clus_size_out = grinchClusSize;
	bb_grinch_tdc_clus_trackindex_out = grinchClusTrIndex;
	bb_grinch_tdc_clus_adc_out = grinchClusADC;
	bb_gem_track_nhits_out = gemNhits;
	bb_gem_track_ngoodhits_out = gemNgoodhits;
	bb_gem_track_chi2ndf_out = gemChiSqr;
	bb_etot_over_p_out = eop;

	P->Fill();	
	
      }//end event loop

      // getting ready for the next run
      C->Reset();

    }//end run loop

  }//end target loop

  // Draw fiducial quality histograms

  // Add TLines to form the rectangle representing the active area
  TLine *line1a = new TLine(econst::hcalposYi_mc, econst::hcalposXi_mc, econst::hcalposYf_mc, econst::hcalposXi_mc); // bottom line
  TLine *line2a = new TLine(econst::hcalposYi_mc, econst::hcalposXf_mc, econst::hcalposYf_mc, econst::hcalposXf_mc); // top line
  TLine *line3a = new TLine(econst::hcalposYi_mc, econst::hcalposXi_mc, econst::hcalposYi_mc, econst::hcalposXf_mc); // left line
  TLine *line4a = new TLine(econst::hcalposYf_mc, econst::hcalposXi_mc, econst::hcalposYf_mc, econst::hcalposXf_mc); // right line
  TLine *line1 = new TLine(hcalaa[2], hcalaa[0], hcalaa[3], hcalaa[0]); // bottom line
  TLine *line2 = new TLine(hcalaa[2], hcalaa[1], hcalaa[3], hcalaa[1]); // top line
  TLine *line3 = new TLine(hcalaa[2], hcalaa[0], hcalaa[2], hcalaa[1]); // left line
  TLine *line4 = new TLine(hcalaa[3], hcalaa[0], hcalaa[3], hcalaa[1]); // right line

  // Set the line color to black for visibility
  line1a->SetLineColor(kBlack);
  line2a->SetLineColor(kBlack);
  line3a->SetLineColor(kBlack);
  line4a->SetLineColor(kBlack);
  line1->SetLineColor(kGreen);
  line2->SetLineColor(kGreen);
  line3->SetLineColor(kGreen);
  line4->SetLineColor(kGreen);

  // Create a canvas to draw the raw expected pos histograms
  TCanvas *c1 = new TCanvas("c1", "Raw expected elastic positions", 1200, 600);
  c1->Divide(2,1);
  c1->cd(1);

  hexpxy_p->Draw("colz");
  line1a->Draw("same");
  line2a->Draw("same");
  line3a->Draw("same");
  line4a->Draw("same");
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  line4->Draw("same");

  c1->cd(2);
  hexpxy_n->Draw("colz");
  line1a->Draw("same");
  line2a->Draw("same");
  line3a->Draw("same");
  line4a->Draw("same");
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  line4->Draw("same");

  c1->Update();

  c1->Write();

  // Create a canvas to draw the raw expected pos histograms
  TCanvas *c2 = new TCanvas("c2", Form("%0.1f sigma fid cut expected elastic positions",Nsig_fid_qual), 1200, 600);
  c2->Divide(2,1);
  c2->cd(1);

  hexpxy_p_fid->Draw("colz");
  line1a->Draw("same");
  line2a->Draw("same");
  line3a->Draw("same");
  line4a->Draw("same");
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  line4->Draw("same");

  c2->cd(2);
  hexpxy_n_fid->Draw("colz");
  line1a->Draw("same");
  line2a->Draw("same");
  line3a->Draw("same");
  line4a->Draw("same");
  line1->Draw("same");
  line2->Draw("same");
  line3->Draw("same");
  line4->Draw("same");

  c2->Update();

  c2->Write();

  fout->Write();

  std::cout << std::endl << "Barebones parsing complete. Output written to " << parse_path << std::endl << std::endl;

  st->Stop();

  // Send time efficiency report to console
  std::cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << std::endl;

}

void checkClusters(
		   const std::vector<double>& cluster_centroids_x, 
		   const std::vector<double>& cluster_centroids_y, 
		   const std::vector<double>& cluster_times
		   ) {

  size_t num_clusters = cluster_centroids_x.size();
  for (size_t i = 0; i < num_clusters; ++i) {
    for (size_t j = i + 1; j < num_clusters; ++j) {
      double dx = cluster_centroids_x[i] - cluster_centroids_x[j];
      double dy = cluster_centroids_y[i] - cluster_centroids_y[j];
      double distance = std::sqrt(dx * dx + dy * dy);

      if (distance < spatial_threshold) {
	double time_difference = std::abs(cluster_times[i] - cluster_times[j]);
	if (time_difference < time_threshold) {
	  // Output centroid difference, time difference, and confirmation
	  std::cout << "Centroid difference: (" << dx << ", " << dy << "), "
		    << "Time difference: " << time_difference << " ns, "
		    << "Pair exists.\n";
	}
      }
    }
  }
}
