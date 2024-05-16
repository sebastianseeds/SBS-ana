//sseeds 10.23.23 - Updated parsing script to cut both inelastic events and unused branches. Configured only to parse data files, not MC. Event parsing using wide globalcuts, wide W2 cuts, and wide coin (HCal/BBCal) cuts. Branch parsing includes only branches that sseeds is using for his gmn analysis
//Update 2.2.24 - Same method, made simpler without class references for troubleshooting
//Version of parse_barebones.C updated to handle monte carlo data.
//csv structure: jobid,Nthrown,Ntried,genvol(MeV*sr^2),luminosity(ub^-1),ebeam(GeV),charge(mC),RndmSeed. jobid is dropped during metadata formation such that each element of metadata objects are as listed - 1. 
//Update 4.2.24 - Added multicluster analysis for higher Q2

#include <vector>
#include <iostream>
#include <regex>
#include <algorithm>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

const int maxClus = 35; //larger than max per tree
const int maxBlk = 25;
const double W2max = 3.0; //very large compared to nucleon mass
const double coin_sigma_factor = 5.; //Wide coincidence timing cut
const double nsig_step = 0.1; //number of sigma to step through in dx until fiducial cut failure written to output tree
const int cluster_method = 4; //Generally useless since we just take the primary cluster in MC
const double hcal_v_offset = 0.; //Should be no offset after pass1 
const double R_mclus = 0.3; //Total radius outside of primary tree cluster where additional clusters are absorbed
const Double_t Nsig_fid_qual = 1.; //number of proton sigma to add to safety margin for fid cut quality plots

//norm_override data
const double Ntried_override = 100000;
const double luminosity_override = 3.8475e+36;
const double genvol_override = 12.566;

//limit_size data
const int maxEvents = 5e6; //5M
Long64_t maxFileSize = 8e9; //8 GB

//Specific wide cut for all parsing
const std::string gcut = "bb.ps.e>0.1&&abs(bb.tr.vz[0])<0.12";

//For now, use alt for all kine except 7, 14, 11

//MAIN
void parse_mc_barebones( int kine = 4, 
			 int mag = 30, 
			 const char *replay_type = "alt", 
			 bool verbose=false, 
			 bool norm_override=false, 
			 bool effz=true, 
			 bool ep_fourvec=false, 
			 bool limit_size=false )
{   

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // Set up exception to MC file structure
  std::string rtype = replay_type;
  std::string rootfile_type = "";
  bool jboyd = false;
  bool alt = false;
  bool alt2 = false;
  if( rtype.compare("jboyd")==0 ){
      rootfile_type = "_jboyd";
      jboyd = true;
  }else if( rtype.compare("alt")==0 ){
      rootfile_type = "_alt";
      alt = true;
  }else if( rtype.compare("alt2")==0 ){
      rootfile_type = "_alt2";
      alt2 = true;
  }else if( rtype.compare("")!=0 )
    cout << "WARNING: invalid argument at replay_type. Valid entries include jboyd, alt, or nothing. Defaulting to nothing." << endl;
    
  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/gmn_mc.json");

  //Get rootfiles and metadata accounting for exceptions
  //All MC for analysis is LD2
  std::string rootfile_dir = jmgr->GetValueFromSubKey_str( Form("filedir_sbs%d%s",kine,rootfile_type.c_str()), Form("%dp",mag) );

  //Check if trailing forward slash is in json, if not add it
  if (!rootfile_dir.empty() && rootfile_dir.back() != '/') {
    rootfile_dir += '/';
  }

  std::string histfile_dir = rootfile_dir;

  if( jboyd ){
    rootfile_dir += Form("simc/SBS%d/",kine);
    histfile_dir += "hist/";
  }else
    histfile_dir += "simcout/";

  std::cout << "Rootfile directory: " << rootfile_dir <<  std::endl;
  std::cout << "Histfile directory: " << histfile_dir <<  std::endl;

  //Get search parameters for MC files
  std::string effz_word = "";
  if(effz)
    effz_word = "_effz";
  std::string four_word = "";
  if(ep_fourvec)
    four_word = "_fourvec";

  std::string partialName_p = jmgr->GetValueFromSubKey_str( Form("p_string_sbs%d%s",kine,rootfile_type.c_str()), Form("%dp",mag) );
  std::string partialName_n = jmgr->GetValueFromSubKey_str( Form("n_string_sbs%d%s",kine,rootfile_type.c_str()), Form("%dp",mag) );

  //set up default parameters for all analysis
  double binfac = 400.;

  // outfile path
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string parse_path = outdir_path + Form("/parse/parse_mc_sbs%d_%dp_barebones%s%s%s.root",kine,mag,rootfile_type.c_str(),effz_word.c_str(),four_word.c_str());

  cout << "Creating parse file and parse log file at " << parse_path << endl;

  //set up output files
  TFile *fout = new TFile( parse_path.c_str(), "RECREATE" );
  std::ofstream logfile(outdir_path + Form("/parse/parselog_mc_sbs%d_%dp_barebones%s.root",kine,mag,rootfile_type.c_str()));

  // Obtain current time
  std::time_t now = std::time(nullptr);
    
  // Convert time_t to tm struct for local time
  std::tm *ltm = std::localtime(&now);

  logfile << "Logfile created for " << parse_path << " at " 
	  << std::setfill('0') << std::setw(2) << 1 + ltm->tm_mon << "."  // Month (01-12)
	  << std::setfill('0') << std::setw(2) << ltm->tm_mday << "."    // Day (01-31)
	  << 1900 + ltm->tm_year << " "                                  // Year
	  << std::setfill('0') << std::setw(2) << ltm->tm_hour << ":"    // Hour (00-23)
	  << std::setfill('0') << std::setw(2) << ltm->tm_min << endl << endl;

  logfile << "Rootfile type: " << rootfile_type << endl;

  //Check neutron root/hist files for generic replay output type
  std::vector<std::string> rootFileNames_n;

  //std::map<int, std::pair<std::string, std::vector<float>>> csvData_n;
  std::vector<std::pair<std::string, std::vector<float>>> metadata_n;

  if(!jboyd)
    util::SyncFilesWithCsv(histfile_dir, rootfile_dir, partialName_n, rootFileNames_n, metadata_n);

  //Check proton root/hist files for generic replay output type
  std::vector<std::string> rootFileNames_p;

  //std::map<int, std::pair<std::string, std::vector<float>>> csvData_p;
  std::vector<std::pair<std::string, std::vector<float>>> metadata_p;

  if(!jboyd)
    util::SyncFilesWithCsv(histfile_dir, rootfile_dir, partialName_p, rootFileNames_p, metadata_p);

  //Check neutron root/hist files for jboyd replay output type
  std::vector<std::string> histFileNames_n;

  if(jboyd)
    util::FindMatchingFiles(histfile_dir,rootfile_dir,partialName_n,histFileNames_n,rootFileNames_n,true);

  //Check proton root/hist files for jboyd replay output type
  std::vector<std::string> histFileNames_p;

  if(jboyd)
    util::FindMatchingFiles(histfile_dir,rootfile_dir,partialName_p,histFileNames_p,rootFileNames_p,true);

  int size_rfiles_raw_p = rootFileNames_p.size();
  int size_rfiles_raw_n = rootFileNames_n.size();  
  int size_hfiles_raw_p = histFileNames_p.size();
  int size_hfiles_raw_n = histFileNames_n.size();

  //Double check that FindMatchingFiles is give 1:1 results
  if(jboyd)
    if( (size_rfiles_raw_p != size_hfiles_raw_p) || 
	(size_rfiles_raw_n != size_hfiles_raw_n) )
      cerr << "ERROR: FindMatchingFiles() failure, vector sizes mismatch" << endl;

  //Throw away jobs for which there does not exist both a proton and neutron rootfile and check if synchronizeJobNumbers() is working.
  int protonfiles = 0;
  int neutronfiles = 0;

  //Make certain that there is a 1:1 correspondance between proton and neutron files. Throw away any unpaired files from analysis
  if(jboyd){
    util::synchronizeJobNumbers(rootFileNames_p,rootFileNames_n);
    protonfiles = rootFileNames_p.size();
    neutronfiles = rootFileNames_n.size();
    
    if( protonfiles!=neutronfiles )
      cout << "ERROR: synchronizeJobNumbers() for jboyd file not working. pfiles/nfiles: " << protonfiles << "/" << neutronfiles << endl;

    for( int i=0; i<protonfiles; ++i ){
      cout << "idx: " << i << " protonfile: " << rootFileNames_p[i] << endl;
      if(jboyd)
	cout << "      protonhistfile: " << histFileNames_p[i] << endl;
      cout << "      neutronfile: " << rootFileNames_n[i] << endl;
      if(jboyd)
	cout << "      neutronhistfile: " << histFileNames_n[i] << endl;

    }

  }else{

    // for( size_t p=0; p<metadata_p.size(); ++p )
    //   cout << metadata_p[p].first << endl;
    // for( size_t n=0; n<metadata_n.size(); ++n )
    //   cout << metadata_n[n].first << endl;


    util::synchronizeJobNumbers(metadata_p,metadata_n);
    protonfiles = metadata_p.size();
    neutronfiles = metadata_n.size();
    
    if( protonfiles!=neutronfiles ){
      logfile << "ERROR: synchronizeJobNumbers() for generic file not working. pfiles/nfiles: " << protonfiles << "/" << neutronfiles << endl;
      cout << "ERROR: synchronizeJobNumbers() for generic file not working. pfiles/nfiles: " << protonfiles << "/" << neutronfiles << endl;
    }
    for( int i=0; i<protonfiles; ++i ){
      logfile << "idx: " << i << " protonfile: " << metadata_p[i].first << endl;
      logfile << "      neutronfile: " << metadata_n[i].first << endl << endl;

      cout << "idx: " << i << " protonfile: " << metadata_p[i].first << endl;
      cout << "      neutronfile: " << metadata_n[i].first << endl << endl;
    }
  }

  //set up diagnostic objects
  int totalNtried_p = 0;
  int totalNtried_n = 0;

  TH1D *hNtried_p = new TH1D("hNtried_p","Number of proton MC events tried by job index;job index",size_rfiles_raw_p,0,size_rfiles_raw_p);
  TH1D *hNtried_n = new TH1D("hNtried_n","Number of neutron MC events tried by job index;job index",size_rfiles_raw_n,0,size_rfiles_raw_n);

  //set up experimental constants
  vector<double> coin_profile = {1000., 0., 5.};

  //Set up hcal active area with bounds that match database on pass
  vector<double> hcalaa = cut::hcalaa_mc(1,1); //verified 2.10.24

  //set up configuration and tune objects to load analysis parameters
  SBSconfig config(kine,mag);

  //Obtain configuration pars from config file
  double hcaltheta = config.GetHCALtheta_rad();
  double hcaldist;
  if(effz){
    hcaldist = config.GetHCALeffdist();
    cout << "Loading effective z offset " << hcaldist << "..." << endl;
  }else
    hcaldist = config.GetHCALdist();  
  double sbsdist = config.GetSBSdist();
  double bbthr = config.GetBBtheta_rad(); //in radians
  double ebeam_tune = config.GetEbeam();

  SBStune tune(kine,mag);

  //Obtain cuts from tune class
  double W2mean   = tune.GetW2mean();
  double W2sig    = tune.GetW2sig();
  double dx0_n    = tune.Getdx0_n();
  double dx0_p    = tune.Getdx0_p();
  double dx_del   = tune.Getdx_del();
  double dy0      = tune.Getdy0();
  double dxsig_n  = tune.Getdxsig_n();
  double dxsig_p  = tune.Getdxsig_p();
  double dysig    = tune.Getdysig();
  double atime0   = tune.Getatime0();
  double atimesig = tune.Getatimesig();
  double atimediff0   = tune.Getatimediff0();
  double atimediffsig = tune.Getatimediffsig();

  //Make fiducial cut quality histograms
  TH2D *hexpxy_p = new TH2D( "hexpxy_p","proton hcal exp x vs hcal exp y; #sigma (m); x_{expect} (m)", 400, -2, 2, 600, -4, 2 );
  TH2D *hexpxy_p_fid = new TH2D( "hexpxy_p_fid",Form("proton hcal exp x vs hcal exp y %0.1f#sigma fid cut; #sigma (m); x_{expect} (m)",Nsig_fid_qual), 400, -2, 2, 600, -4, 2 );

  TH2D *hexpxy_n = new TH2D( "hexpxy_n","neutron hcal exp x vs hcal exp y; #sigma (m); x_{expect} (m)", 400, -2, 2, 600, -4, 2 );
  TH2D *hexpxy_n_fid = new TH2D( "hexpxy_n_fid",Form("neutron hcal exp x vs hcal exp y %0.1f#sigma fid cut; #sigma (m); x_{expect} (m)",Nsig_fid_qual), 400, -2, 2, 600, -4, 2 );

  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;

  // create output tree
  TTree *P = new TTree("P","Analysis Tree"); 

  // limit the output so that root files don't grow to ridiculous sizes.
  // Set maximum tree size to 1 GB
  TTree::SetMaxTreeSize(maxFileSize);

  // new output tree vars
  //Universals
  int job_out;
  int event_out;
  int trig_out;
  double xexp_out;
  double yexp_out;
  double W2_out;
  double Q2_out;
  double nu_out;
  double tau_out;
  double epsilon_out;
  double precon_out;

  //MC data
  int nucleon_out;
  int ntried_out;
  double genvol_out;
  double lumi_out;
  double ebeam_out;
  double charge_out;
  double seed_out;
  double mc_weight_out;
  double mc_weight_norm_out;

  //Fiducial slices
  double fiducial_sig_x_out;
  double fiducial_sig_y_out;

  //Primary cluster (default, containing highest energy block)
  double dx_out;
  double dy_out;
  double coin_out;
  double thetapq_pout;
  double thetapq_nout;
  double sigsep_pout; //Total scale N, where N*sig_dx_p and N*sig_dy_p (and mean dx_p, dy_p) defines the minimum area necessary to include point on event (dx,dy)
  double sigsep_nout; //Total scale N, where N*sig_dx_n and N*sig_dy_n (and mean dx_n, dy_n) defines the minimum area necessary to include point on event (dx,dy)
  int hcalon_out;
  double accsep_xout; //Total scale N, where N*block_width_x necessary to reject event on active area cut
  double accsep_yout; //Total scale N, where N*block_width_y necessary to reject event on active area cut
  double hcalnblk_out;
  double hcalnclus_out;
  double hcalpid_out;
  double hcalx_out;
  double hcaly_out;
  double hcale_out;  
  double hcal_index_out;

  //Best cluster
  double dx_bc_out;
  double dy_bc_out;
  double coin_bc_out;
  double thetapq_bc_pout;
  double thetapq_bc_nout;
  double sigsep_bc_pout; //Total scale N, where N*sig_dx_p and N*sig_dy_p (and mean dx_p, dy_p) defines the minimum area necessary to include point on event (dx,dy)
  double sigsep_bc_nout; //Total scale N, where N*sig_dx_n and N*sig_dy_n (and mean dx_n, dy_n) defines the minimum area necessary to include point on event (dx,dy)  int hcalon_bc_out;
  int hcalon_bc_out;
  double accsep_bc_xout; //Total scale N, where N*block_width_x necessary to reject event on active area cut
  double accsep_bc_yout; //Total scale N, where N*block_width_y necessary to reject event on active area cut
  double hcalnblk_bc_out;
  double hcalpid_bc_out;
  double hcalx_bc_out;
  double hcaly_bc_out;
  double hcale_bc_out; 

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
  double bb_tr_vz_out;
  double bb_tr_n_out;
  double bb_tr_p_out;
  double bb_tr_th_out;
  double bb_tr_ph_out;
  double bb_tr_r_th_out;
  double bb_tr_r_x_out;
  double bb_tr_r_ph_out;
  double bb_tr_r_y_out;
  double bb_tr_chi2_out;
  double bb_ps_e_out;
  double bb_ps_rowblk_out;
  double bb_ps_colblk_out;
  double bb_sh_e_out;
  double bb_sh_rowblk_out;
  double bb_sh_colblk_out;
  double bb_hodotdc_clus_tmean_out;
  double bb_grinch_tdc_clus_size_out;
  double bb_grinch_tdc_clus_trackindex_out;
  double bb_gem_track_nhits_out;
  double bb_gem_track_ngoodhits_out;
  double bb_gem_track_chi2ndf_out;
  double bb_etot_over_p_out;

  P->Branch("job", &job_out, "job/I");
  P->Branch("event", &event_out, "event/I");
  P->Branch("trig", &trig_out, "trig/I");
  P->Branch("xexp", &xexp_out, "xexp/D");
  P->Branch("yexp", &yexp_out, "yexp/D");
  P->Branch("W2", &W2_out, "W2/D");
  P->Branch("Q2", &Q2_out, "Q2/D");
  P->Branch("nu", &nu_out, "nu/D");
  P->Branch("tau", &tau_out, "tau/D");
  P->Branch("epsilon", &epsilon_out, "epsilon/D");
  P->Branch("precon", &precon_out, "precon/D");

  P->Branch("nucleon", &nucleon_out, "nucleon/I");
  P->Branch("ntried", &ntried_out, "ntried/I");
  P->Branch("genvol", &genvol_out, "genvol/D");
  P->Branch("lumi", &lumi_out, "lumi/D");
  P->Branch("ebeam", &ebeam_out, "ebeam/D");
  P->Branch("charge", &charge_out, "charge/D");
  P->Branch("seed", &seed_out, "seed/D");
  P->Branch("mc_weight", &mc_weight_out, "mc_weight/D");
  P->Branch("mc_weight_norm", &mc_weight_norm_out, "mc_weight_norm/D");

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
  P->Branch("hcalnclus", &hcalnclus_out, "hcalnclus/D");
  P->Branch("hcalpid", &hcalpid_out, "hcalpid/D");
  P->Branch("hcalx", &hcalx_out, "hcalx/D");
  P->Branch("hcaly", &hcaly_out, "hcaly/D");
  P->Branch("hcale", &hcale_out, "hcale/D");
  P->Branch("hcal_index", &hcal_index_out, "hcal_index/D");

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
  P->Branch("bb_hodotdc_clus_tmean", &bb_hodotdc_clus_tmean_out, "bb_hodotdc_clus_tmean/D");
  P->Branch("bb_grinch_tdc_clus_size", &bb_grinch_tdc_clus_size_out, "bb_grinch_tdc_clus_size/D");
  P->Branch("bb_grinch_tdc_clus_trackindex", &bb_grinch_tdc_clus_trackindex_out, "bb_grinch_tdc_clus_trackindex/D");
  P->Branch("bb_gem_track_nhits", &bb_gem_track_nhits_out, "bb_gem_track_nhits/D");
  P->Branch("bb_gem_track_ngoodhits", &bb_gem_track_ngoodhits_out, "bb_gem_track_ngoodhits/D");
  P->Branch("bb_gem_track_chi2ndf", &bb_gem_track_chi2ndf_out, "bb_gem_track_chi2ndf/D");
  P->Branch("bb_etot_over_p", &bb_etot_over_p_out, "bb_etot_over_p/D");

  std::string nuc = "none";
  int nnuc = 2; //two nucleons, proton and neutron
  int nfiles = 0;

  int protonevents = 0;
  int neutronevents = 0;

  //Main loop over nucleons (r==0, proton; r==1, neutron)
  for (int r=0; r<nnuc; r++) {

    //update current nucleon
    if( r==0 ){
      nuc = "p";
      nfiles = protonfiles;
      if( limit_size && protonevents>maxEvents ) //if over max limit on events, continue
	continue;
    }else if( r==1 ){
      nuc = "n";
      nfiles = neutronfiles;
      if( limit_size && neutronevents>maxEvents )
	continue;
    }

    //loop over simulation files
    for( int f=0; f<nfiles; ++f ){

      std::string load_file;
      int job_number = -1;
      int Ntried = -1;
      int Nthrown = -1;
      double luminosity = -1.;
      double genvol = -1.;
      double ebeam = -1.;
      double charge = -1.;
      double seed = -1.;
      bool using_rej_samp = 0;
      double max_weight = -1;
      double weight_gt_max = 0.;
      double obs_max_weight = 0.;

      //Get rootfile paths and metadata
      if(r==0){  //proton
	//for jboyd structure
	if(jboyd){
	  load_file = rootFileNames_p[f];
	  std::string hist_file = histFileNames_p[f];

	  Ntried = (int)util::searchSimcHistFile("Ntried", hist_file);
	  luminosity = util::searchSimcHistFile("luminosity", hist_file);
	  genvol = util::searchSimcHistFile("genvol", hist_file);

	}else{
	  load_file = metadata_p[f].first;
	  
	  job_number = metadata_p[f].second[0];
	  Nthrown = metadata_p[f].second[1];
	  Ntried = metadata_p[f].second[2]; //see csv structure in preamble
	  genvol = metadata_p[f].second[3];
	  luminosity = metadata_p[f].second[4];
	  ebeam =  metadata_p[f].second[5];
	  charge =  metadata_p[f].second[6];
	  seed =  metadata_p[f].second[7];
	  if(metadata_p[f].second.size() > 8) {  //added just in case old mc file is loaded
	    using_rej_samp = metadata_p[f].second[8];
	    max_weight = metadata_p[f].second[9];
	    weight_gt_max = metadata_p[f].second[10];
	    obs_max_weight = metadata_p[f].second[11];
	  }
	}
	
	cout << "Loaded proton file at " << load_file << endl;
	logfile << "Loaded PROTON file at " << load_file << endl;
	

      }else if(r==1){  //neutron
	//jboyd structure
	if(jboyd){
	  load_file = rootFileNames_n[f];
	  std::string hist_file = histFileNames_n[f];

	  Ntried = (int)util::searchSimcHistFile("Ntried", hist_file);
	  luminosity = util::searchSimcHistFile("luminosity", hist_file);
	  genvol = util::searchSimcHistFile("genvol", hist_file);

	}else{
	  load_file = metadata_n[f].first;
	  
	  job_number = metadata_n[f].second[0];
	  Nthrown = metadata_n[f].second[1];
	  Ntried = metadata_n[f].second[2]; //see csv structure in preamble
	  genvol = metadata_n[f].second[3];
	  luminosity = metadata_n[f].second[4];
	  ebeam =  metadata_n[f].second[5];
	  charge =  metadata_n[f].second[6];
	  seed =  metadata_n[f].second[7];
	  if(metadata_n[f].second.size() > 8) {  //added just in case old mc file is loaded
	    using_rej_samp = metadata_n[f].second[8];
	    max_weight = metadata_n[f].second[9];
	    weight_gt_max = metadata_n[f].second[10];
	    obs_max_weight = metadata_n[f].second[11];
	  }
	}

	cout << "Loaded neutron file at " << load_file << endl;
	logfile << "Loaded NEUTRON file at " << load_file << endl;

      }

      //if override is enabled, change metadata
      if(norm_override){
	Ntried = Ntried_override; //see csv structure in preamble
	genvol = genvol_override;
	luminosity = luminosity_override;
      }

      cout << "For job " << job_number <<  ", Ntried=" << Ntried << " luminosity=" << luminosity << " genvol=" << genvol <<  " max weight=" << max_weight << endl; 
      logfile << "For job " << job_number <<  ", Ntried=" << Ntried << " luminosity=" << luminosity << " genvol=" << genvol << " max weight=" << max_weight <<  endl; 

      //fill diagnostics
      if(r==0){
	totalNtried_p += Ntried;
	hNtried_p->SetBinContent(job_number,Ntried);
      }else if(r==1){
	totalNtried_n += Ntried;
	hNtried_n->SetBinContent(job_number,Ntried);
      }

      if( !norm_override && (Ntried==-1 || luminosity==-1 || genvol==-1) ){
	cout << "WARNING: Failed to read normalization data, skipping file.." << endl;
	continue;
      }

      //send values of ebeam to console for comparison
      ebeam /= 1000.; //convert to GeV

      if( ebeam == -1 )
	ebeam = ebeam_tune;

      //first attempt to resolve segmentation fault on large data sets
      if (C != nullptr) {
	delete C;
      }

      C = new TChain("T");
      C->Add(load_file.c_str());

      //Reporting on farm running
      // Get current system time
      auto now = std::chrono::system_clock::now();
    
      // Convert to time_t which is needed for converting to tm (time structure)
      std::time_t now_time_t = std::chrono::system_clock::to_time_t(now);
    
      // Convert to local time
      std::tm* now_tm = std::localtime(&now_time_t);
    
      // Print the time in "HH:MM AM/PM" format
      std::cout << "Switched to file " << load_file << " at " << std::put_time(now_tm, "%I:%M %p") << std::endl;
      logfile << "Switched to file " << load_file << " at " << std::put_time(now_tm, "%I:%M %p") << std::endl;

      // setting up ROOT tree branch addresses
      C->SetBranchStatus("*",0);    

      // HCal general
      double hcalid, hcale, hcalx, hcaly, hcalr, hcalc, hcaltdc, hcalatime, hcalidx, nclus, nblk;
      std::vector<std::string> hcalvar = {"idblk","e","x","y","rowblk","colblk","tdctimeblk","atimeblk","index","nclus","nblk"};
      std::vector<void*> hcalvarlink = {&hcalid,&hcale,&hcalx,&hcaly,&hcalr,&hcalc,&hcaltdc,&hcalatime,&hcalidx,&nclus,&nblk};
      rvars::setbranch(C, "sbs.hcal", hcalvar, hcalvarlink);
      
      // HCal cluster branches
      double hcalcid[econst::maxclus], hcalce[econst::maxclus], hcalcx[econst::maxclus], hcalcy[econst::maxclus], hcalctdctime[econst::maxclus], hcalcatime[econst::maxclus], hcalcrow[econst::maxclus], hcalcnblk[econst::maxclus], hcalccol[econst::maxclus];
      int Nhcalcid;
      std::vector<std::string> hcalcvar = {"id","e","x","y","tdctime","atime","row","col","nblk","id"};
      std::vector<void*> hcalcvarlink = {&hcalcid,&hcalce,&hcalcx,&hcalcy,&hcalctdctime,&hcalcatime,&hcalcrow,&hcalccol,&hcalcnblk,&Nhcalcid};
      rvars::setbranch(C, "sbs.hcal.clus", hcalcvar, hcalcvarlink, 9);

      // HCal cluster blk branches
      double hcalcbid[econst::maxclus], hcalcbe[econst::maxclus], hcalcbx[econst::maxclus], hcalcby[econst::maxclus], hcalcbtdctime[econst::maxclus], hcalcbatime[econst::maxclus];
      int Nhcalcbid;
      std::vector<std::string> hcalcbvar = {"id","e","x","y","tdctime","atime","id"};
      std::vector<void*> hcalcbvarlink = {&hcalcbid,&hcalcbe,&hcalcbx,&hcalcby,&hcalcbtdctime,&hcalcbatime,&Nhcalcbid};
      rvars::setbranch(C, "sbs.hcal.clus_blk", hcalcbvar, hcalcbvarlink, 6);

      // bbcal clus var
      double eSH, xSH, ySH, rblkSH, cblkSH, idblkSH, atimeSH, nclusSH, ePS, rblkPS, cblkPS, idblkPS, atimePS;
      std::vector<std::string> bbcalclvar = {"sh.e","sh.x","sh.y","sh.rowblk","sh.colblk","sh.idblk","sh.atimeblk","sh.nclus","ps.e","ps.rowblk","ps.colblk","ps.idblk","ps.atimeblk"};
      std::vector<void*> bbcalclvarlink = {&eSH,&xSH,&ySH,&rblkSH,&cblkSH,&idblkSH,&atimeSH,&nclusSH,&ePS,&rblkPS,&cblkPS,&idblkPS,&atimePS};
      rvars::setbranch(C, "bb", bbcalclvar, bbcalclvarlink);

      // hodoscope cluster mean time
      int Nhodotmean; 
      double hodotmean[econst::maxclus];
      std::vector<std::string> hodovar = {"clus.tmean","clus.tmean"};
      std::vector<void*> hodovarlink = {&Nhodotmean,&hodotmean};
      rvars::setbranch(C, "bb.hodotdc", hodovar, hodovarlink, 0); 

      // track branches
      double ntrack, p[econst::maxtrack],px[econst::maxtrack],py[econst::maxtrack],pz[econst::maxtrack],xtr[econst::maxtrack],ytr[econst::maxtrack],thtr[econst::maxtrack],phtr[econst::maxtrack];
      double vx[econst::maxtrack],vy[econst::maxtrack],vz[econst::maxtrack];
      double xtgt[econst::maxtrack],ytgt[econst::maxtrack],thtgt[econst::maxtrack],phtgt[econst::maxtrack];
      double r_x[econst::maxtrack],r_th[econst::maxtrack],r_y[econst::maxtrack],r_ph[econst::maxtrack],tr_chi2[econst::maxtrack];
      std::vector<std::string> trvar = {"n","p","px","py","pz","x","y","th","ph","vx","vy","vz","tg_x","tg_y","tg_th","tg_ph","r_x","r_th","r_y","r_ph","chi2"};
      std::vector<void*> trvarlink = {&ntrack,&p,&px,&py,&pz,&xtr,&ytr,&thtr,&phtr,&vx,&vy,&vz,&xtgt,&ytgt,&thtgt,&phtgt,&r_x,&r_th,&r_y,&r_ph,&tr_chi2};
      rvars::setbranch(C,"bb.tr",trvar,trvarlink);

      // tdctrig branches
      // int Ntdctrigid;
      // double tdctrig[econst::maxtrack], tdctrigid[econst::maxtrack];
      // std::vector<std::string> tdcvar = {"tdcelemID","tdcelemID","tdc"};
      // std::vector<void*> tdcvarlink = {&tdctrigid,&Ntdctrigid,&tdctrig};
      // rvars::setbranch(C,"bb.tdctrig",tdcvar,tdcvarlink,1);

      // ekine branches
      double ekineQ2, ekineW2, ekineeps, ekinenu, ekineqx, ekineqy, ekineqz;
      std::vector<std::string> ekinevar = {"Q2","W2","epsilon","nu","q_x","q_y","q_z"};
      std::vector<void*> ekinevarlink = {&ekineQ2,&ekineW2,&ekineeps,&ekinenu,&ekineqx,&ekineqy,&ekineqz};
      rvars::setbranch(C, "e.kine", ekinevar, ekinevarlink);
    
      // fEvtHdr branches
      UInt_t rnum, gevnum, trigbits;
      std::vector<std::string> evhdrvar = {"fRun","fEvtNum","fTrigBits"};
      std::vector<void*> evhdrlink = {&rnum,&gevnum,&trigbits};
      rvars::setbranch(C,"fEvtHdr",evhdrvar,evhdrlink);

      // other bb branches
      Double_t gemNhits, gemNgoodhits, gemChiSqr, grinchClusSize, grinchClusTrIndex, eop;
      std::vector<std::string> miscbbvar = {"gem.track.nhits","gem.track.ngoodhits","gem.track.chi2ndf","grinch_tdc.clus.size","grinch_tdc.clus.trackindex","etot_over_p"};
      std::vector<void*> miscbbvarlink = {&gemNhits,&gemNgoodhits,&gemChiSqr,&grinchClusSize,&grinchClusTrIndex,&eop};
      rvars::setbranch(C, "bb", miscbbvar, miscbbvarlink);

      // Monte Carlo variables
      double mcweight;
      std::vector<std::string> mcvar = {"simc_Weight"};
      std::vector<void*> mcvarlink = {&mcweight};
      rvars::setbranch(C, "MC", mcvar, mcvarlink);

      TCut GCut = gcut.c_str();

      TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );

      // get experimental quantities by run
      //set up hcal coordinate system with hcal angle wrt exit beamline
      vector<TVector3> hcalaxes; vars::sethcalaxes( hcaltheta, hcalaxes );
      TVector3 hcalorigin = hcaldist*hcalaxes[2] + hcal_v_offset*hcalaxes[0];

      double BdL = econst::sbsmaxfield * econst::sbsdipolegap * (mag/100); //scale crudely by magnetic field
      double Eloss_outgoing = econst::celldiameter/2.0/sin(bbthr) * econst::ld2tarrho * econst::ld2dEdx;

      // set nucleon defaulted to LD2
      std::string nucleon = "np";

      // event indices
      long nevent = 0, npassed = 0, nevents = C->GetEntries(); 
      int treenum = 0, currenttreenum = 0;

      cout << "Beginning analysis of " << nuc << " rootfile " << load_file << endl;
      logfile << "Beginning analysis of " << nuc << " rootfile " << load_file << endl;
      
      //Event loop
      while (C->GetEntry(nevent++)) {
	
	// std::cout << "Processing job " << f << " event " << nevent << " / " << nevents << ", total passed cuts " << npassed << ". Total proton events: " << protonevents << ". Total neutron events: " << neutronevents << "\r";
	// std::cout.flush();

	//
	if( r==0 && limit_size && protonevents>maxEvents )
	  break;
	if( r==1 && limit_size && neutronevents>maxEvents )
	  break;

	///////
	//Single-loop globalcut method. Save pass/fail for output tree.
	bool failedglobal = false;

	currenttreenum = C->GetTreeNumber();
	if( nevent == 1 || currenttreenum != treenum ){
	  treenum = currenttreenum; 
	  GlobalCut->UpdateFormulaLeaves();
	}

	failedglobal = GlobalCut->EvalInstance(0) == 0;
	if( failedglobal ){
	  continue;
	}

	///////
	//Physics calculations
	//correct beam energy with vertex information
	double ebeam_c;
	ebeam_c = vars::ebeam_c( ebeam, vz[0], "LD2" );

	TVector3 vertex( 0., 0., vz[0] );

	//reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
	double precon = p[0] + Eloss_outgoing;

	//set up four-momenta with some empty for various calculation methods
	TLorentzVector pbeam( 0., 0., ebeam_c, ebeam_c ); //beam momentum
	//TLorentzVector pe( px[0], py[0], pz[0], p[0] ); //e' momentum
	TLorentzVector pe( precon*px[0]/p[0], precon*py[0]/p[0], precon*pz[0]/p[0], precon ); //e' recon plvect
	TLorentzVector ptarg; vars::setPN(nucleon,ptarg); //target momentum
	TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
	TLorentzVector pN; //N' momentum
	TVector3 pNhat; //Unit N' 3-vector
      
	//simple calculations for e' and N'
	double etheta = vars::etheta(pe); 
	double ephi = vars::ephi(pe);
	double pcent = vars::pcentral(ebeam,etheta,nucleon); //e' p reconstructed by angles
	double phNexp = ephi + physconst::pi;
	double Q2, W2, nu, thNexp, pNexp, ebeam_o, tau, epsilon;
	ebeam_o = vars::ebeam_o( ebeam_c, etheta, "LD2" ); //Second energy correction accounting for energy loss leaving target
       
	//reconstruct track with angles (better resolution with GEMs)
	nu = pbeam.E() - pcent;
	pNexp = vars::pN_expect( nu, nucleon );
	//thNexp = acos((ebeam-pz[0])/pNexp);
	thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	pNhat = vars::pNhat_track( thNexp, phNexp );
	pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	Q2 = vars::Q2( pbeam.E(), pe.E(), etheta );
	W2 = vars::W2( pbeam.E(), pe.E(), Q2, nucleon );
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

	// if(ep_fourvec)
	//   thNexp = acos((ebeam-pz[0])/pNexp);


	/////////////////////
	//W2 elastic cut
	bool failedW2 = W2>W2max;
	if(failedW2)
	  continue;

	npassed++;

	double comp_ev_fraction = (double)npassed/(double)nevent;
	double ev_fraction = (double)npassed/(double)nevents;

	//Calculate h-arm quantities
	vector<double> xyhcalexp; vars::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
	TVector3 hcalpos = hcalorigin + hcalx*hcalaxes[0] + hcaly*hcalaxes[1];
	double dx = hcalx - xyhcalexp[0];
	double dy = hcaly - xyhcalexp[1];
	TVector3 neutdir = ( hcalpos - vertex ).Unit();
	double protdeflect = tan( 0.3 * BdL / q.Mag() ) * (hcaldist - (sbsdist + econst::sbsdipolegap/2.0) );
	TVector3 protdir = ( hcalpos + protdeflect * hcalaxes[0] - vertex ).Unit();
	double thetapq_p = acos( protdir.Dot( pNhat ) );
	double thetapq_n = acos( neutdir.Dot( pNhat ) );
	double eoverp = ( ePS + eSH ) / p[0];

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

	// Set up the structure to hold energy and index, then sort based on energy
	std::vector<std::pair<double, int>> energy_index_pairs;
	
	//Add multicluster analysis for higher Q2 where hcal clusters may be separated
	vector<int> intime_cluster_indices;

	//loop through all clusters and select without HCal position information
	for( int c=0; c<Nhcalcid; c++ ){
	
	  //calculate h-arm physics quantities per cluster
	  double atime = hcalcatime[c];
	  double atime_diff = atime - atimeSH; //Assuming best shower time on primary cluster
	  double ce = hcalce[c];

	  //using hcal atime until after pass2
	  bool passedCoin = abs(atime_diff-coin_profile[1])<coin_sigma_factor*coin_profile[2];

	  //Replicate the in-time algorithm with new cluster to be sorted later
	  clone_cluster_intime.push_back(ce);
	  if( !passedCoin ){
	    clone_cluster_intime[c] = 0;
	  }else{
	    intime_cluster_indices.push_back(c);
	  }

	  //Get score (no position info). Will be sorted later
	  double cascore = util::assignScore( ce, atime_diff, hcalce[(int)hcalidx], coin_profile);
	  clone_cluster_score.push_back(cascore);
	  
	  //Add energy cluster index and energy to pair
	  energy_index_pairs.push_back(std::make_pair(ce, c));
	
	}//endloop over cluster elements

	// Sort in descending order based on the cluster energy
	std::sort(energy_index_pairs.begin(), energy_index_pairs.end(), [](const std::pair<double, int>& a, const std::pair<double, int>& b) {
	    return a.first > b.first; // Using the first element of pair (energy) for comparison
	  });
	
	// cout << "total clusters on this event " << Nhcalcid << endl;

	// for (size_t i=0; i<energy_index_pairs.size(); ++i)
	//   cout << "Sorted cluster index " << energy_index_pairs[i].second << " with energy " << energy_index_pairs[i].first << endl;

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

	if( nclus_multicluster>1 )
	  cout << "  Multiple correlated clusters = " << nclus_multicluster << " with variables: energy = " << mclus_e << ", x " << mclus_x << ", y " << mclus_y << endl;

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

	if( hcale != hcalce[(int)hcalidx] && verbose )
	  cerr << "WARNING: Sorting failure on cluster index " << (int)hcalidx << ". Check max clusters allowed on tree." << endl;
	
	//Get best score/intime indices from clone clusters
	int score_idx = -1;
	int intime_idx = -1;
	double score = 0.;
	double intime = 0.;
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

	int cidx_intime = intime_idx;
	int cidx_score = score_idx;

	//Switch between best clusters for systematic analysis
	int cidx_best;
      
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
	double thetapq_bc_p = acos( protdir_bc.Dot( pNhat ) );
	double thetapq_bc_n = acos( neutdir_bc.Dot( pNhat ) );

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
	std::pair<double, double> fiducial_factors = cut::findFidFailure(dxsig_p, 
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

	////////////////////
	//coincidence time BBCal/HCal cut
	bool failedwidecoin_bc = abs( hatime_bestcluster - atimediff0 ) > coin_sigma_factor*atimediffsig;

	//caculate final weight for this event
	Double_t final_mc_weight;

	if(using_rej_samp)
	  final_mc_weight = max_weight*luminosity*genvol/Ntried;
	else
	  final_mc_weight = mcweight*luminosity*genvol/Ntried;
	  
	//Fill new output tree  
	//Universals
	job_out = job_number;
	event_out = (int)gevnum;
	trig_out = (int)trigbits;
	xexp_out = xyhcalexp[0];
	yexp_out = xyhcalexp[1];
	W2_out = W2;
	Q2_out = Q2;
	nu_out = nu;
	tau_out = tau;
	epsilon_out = epsilon;
	precon_out = precon;

	//MC data
	nucleon_out = r;
	ntried_out = Ntried;
	genvol_out = genvol;
	lumi_out = luminosity;
	ebeam_out = ebeam_c;
	charge_out = charge;
	seed_out = seed;
	mc_weight_out = mcweight;
	mc_weight_norm_out = final_mc_weight;

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
	hcal_index_out = hcalidx;

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
	bb_hodotdc_clus_tmean_out = hodotmean[0];
	bb_grinch_tdc_clus_size_out = grinchClusSize;
	bb_grinch_tdc_clus_trackindex_out = grinchClusTrIndex;
	bb_gem_track_nhits_out = gemNhits;
	bb_gem_track_ngoodhits_out = gemNgoodhits;
	bb_gem_track_chi2ndf_out = gemChiSqr;
	bb_etot_over_p_out = eop;
	hcalnclus_out = Nhcalcid;

	P->Fill();	
	
	if(r==0)
	  protonevents++;
	else
	  neutronevents++;

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

  logfile.close();

  std::cout << std::endl << "Barebones parsing complete. Output written to " << parse_path << std::endl << std::endl;
  logfile << std::endl << "Barebones parsing complete. Output written to " << parse_path << std::endl << std::endl;

  st->Stop();

  // Send time efficiency report to console
  std::cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << std::endl;

}
