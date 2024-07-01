//sseeds 10.23.23 - Updated parsing script to cut both inelastic events and unused branches. Configured only to parse data files, not MC. Event parsing using wide globalcuts, wide W2 cuts, and wide coin (HCal/BBCal) cuts. Branch parsing includes only branches that sseeds is using for his gmn analysis
//Update 2.2.24 - Same method, made simpler without class references for troubleshooting
//Version of parse_barebones.C updated to handle monte carlo data.
//csv structure: jobid,Nthrown,Ntried,genvol(MeV*sr^2),luminosity(ub^-1),ebeam(GeV),charge(mC),RndmSeed. jobid is dropped during metadata formation such that each element of metadata objects are as listed - 1. 
//This version focused on g4sbs inelastic events only

#include <vector>
#include <iostream>
#include <regex>
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
const int cluster_method = 4; //Hard coded for MC where cluster selection isn't essential
const double hcal_v_offset = 0.; //Should be no offset in MC data

//norm_override data
const double Ntried_override = 100000;
const double luminosity_override = 3.8475e+36;
const double genvol_override = 12.566;

//limit_size data
const int maxEvents = 5000000; //5M
Long64_t maxFileSize = 8000000000; //8 GB

//Specific wide cut for all parsing
const std::string gcut = "bb.ps.e>0.1&&abs(bb.tr.vz[0])<0.12";

//MAIN. kine=kinematic, mag=magnetic field setting (percent), verbose=added console output, norm_override=use single normalization value from global var for all events, effz=use effective z offset, limit_size=segment output to limit size of output files with global var
void parse_mcbg_barebones( int kine=8, int mag = 100, bool verbose=false, bool norm_override=false, bool effz=true, bool limit_size=false )
{   

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/gmn_mc.json");

  //Get rootfiles and metadata accounting for exceptions
  //All MC for analysis is LD2
  std::string rootfile_dir = jmgr->GetValueFromSubKey_str( Form("filedir_inel_sbs%d",kine), Form("%dp",mag) );

  //Get search parameters for MC files
  std::string partialName = jmgr->GetValueFromSubKey_str( Form("inel_string_sbs%d",kine), Form("%dp",mag) );

  std::string fileset = rootfile_dir + "*" + partialName + "*.root";

  //set up default parameters for all analysis
  double binfac = 400.;

  // outfile path
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string parse_path = outdir_path + Form("/parse/parse_mcbg_sbs%d_%dp_barebones.root",kine,mag);

  //set up output files
  TFile *fout = new TFile( parse_path.c_str(), "RECREATE" );

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
  double ebeam = config.GetEbeam();

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

  // create output tree
  TTree *P = new TTree("P","Analysis Tree"); 

  // limit the output so that root files don't grow to ridiculous sizes.
  // Set maximum tree size to 1 GB
  TTree::SetMaxTreeSize(maxFileSize);

  // new output tree vars
  //Universals
  int event_out;
  int trig_out;
  int failedglobal_out;
  int failedW2_out;
  double xexp_out;
  double yexp_out;
  double W2_out;
  double Q2_out;
  double nu_out;
  double precon_out;

  //MC data
  double mcsigma_out;
  double mcomega_out;
  double mcweight_out;
  double mcnucl_out;

  //Fiducial slices
  double fiducial_sig_x_out;
  double fiducial_sig_y_out;

  //Primary cluster
  double dx_out;
  double dy_out;
  double coin_out;
  double thetapq_pout;
  double thetapq_nout;
  int nucleon_spot_out; //via data ellipse inclusion. -1=outside both, 0=inside both, 1=proton, 2=neutron
  int hcalon_out;
  double hcalnblk_out;
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
  int nucleon_bc_out; //via data ellipse inclusion. -1=outside both, 0=inside both, 1=proton, 2=neutron
  int hcalon_bc_out;
  double hcalnblk_bc_out;
  double hcalpid_bc_out;
  double hcalx_bc_out;
  double hcaly_bc_out;
  double hcale_bc_out;  

  // relevant old output tree vars
  double bb_tr_vz_out;
  double bb_tr_n_out;
  double bb_tr_p_out;
  double bb_tr_th_out;
  double bb_tr_ph_out;
  double bb_tr_r_th_out;
  double bb_tr_r_x_out;
  double bb_ps_e_out;
  double bb_ps_rowblk_out;
  double bb_ps_colblk_out;
  double bb_sh_e_out;
  double bb_sh_rowblk_out;
  double bb_sh_colblk_out;
  //double bb_hodotdc_clus_tmean_out;
  //double bb_grinch_tdc_clus_size_out;
  //double bb_grinch_tdc_clus_trackindex_out;
  double bb_gem_track_nhits_out;
  double bb_gem_track_ngoodhits_out;
  double bb_gem_track_chi2ndf_out;
  double bb_etot_over_p_out;

  P->Branch("event", &event_out, "event/I");
  P->Branch("trig", &trig_out, "trig/I");
  P->Branch("failedglobal", &failedglobal_out, "failedglobal/I");
  P->Branch("failedW2", &failedW2_out, "failedW2/I");
  P->Branch("xexp", &xexp_out, "xexp/D");
  P->Branch("yexp", &yexp_out, "yexp/D");
  P->Branch("W2", &W2_out, "W2/D");
  P->Branch("Q2", &Q2_out, "Q2/D");
  P->Branch("nu", &nu_out, "nu/D");
  P->Branch("precon", &precon_out, "precon/D");

  P->Branch("mcsigma", &mcsigma_out, "mcsigma/D");
  P->Branch("mcomega", &mcomega_out, "mcomega/D");
  P->Branch("mcweight", &mcweight_out, "mcweight/D");
  P->Branch("mcnucl", &mcnucl_out, "mcnucl/D");

  P->Branch("fiducial_sig_x", &fiducial_sig_x_out, "fiducial_sig_x/D");
  P->Branch("fiducial_sig_y", &fiducial_sig_y_out, "fiducial_sig_y/D");

  P->Branch("dx", &dx_out, "dx/D");
  P->Branch("dy", &dy_out, "dy/D");
  P->Branch("coin", &coin_out, "coin/D");
  P->Branch("thetapq_p", &thetapq_pout, "thetapq_p/D");
  P->Branch("thetapq_n", &thetapq_nout, "thetapq_n/D");
  P->Branch("nucleon_spot", &nucleon_spot_out, "nucleon_spot/I");
  P->Branch("hcalon", &hcalon_out, "hcalon/I");
  P->Branch("hcalnblk", &hcalnblk_out, "hcalnblk/D");
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
  P->Branch("nucleon_bc", &nucleon_bc_out, "nucleon_bc/I");
  P->Branch("hcalon_bc", &hcalon_bc_out, "hcalon_bc/I");
  P->Branch("hcalnblk_bc", &hcalnblk_bc_out, "hcalnblk_bc/D");
  P->Branch("hcalpid_bc", &hcalpid_bc_out, "hcalpid_bc/D");
  P->Branch("hcalx_bc", &hcalx_bc_out, "hcalx_bc/D");
  P->Branch("hcaly_bc", &hcaly_bc_out, "hcaly_bc/D");
  P->Branch("hcale_bc", &hcale_bc_out, "hcale_bc/D");

  P->Branch("bb_tr_n", &bb_tr_n_out, "bb_tr_n/D");
  P->Branch("bb_tr_vz", &bb_tr_vz_out, "bb_tr_vz/D");
  P->Branch("bb_tr_p", &bb_tr_p_out, "bb_tr_p/D");
  P->Branch("bb_tr_th", &bb_tr_th_out, "bb_tr_th/D");
  P->Branch("bb_tr_ph", &bb_tr_ph_out, "bb_tr_ph/D");
  P->Branch("bb_tr_r_x", &bb_tr_r_x_out, "bb_tr_r_x/D");
  P->Branch("bb_tr_r_th", &bb_tr_r_th_out, "bb_tr_r_th/D");
  P->Branch("bb_ps_e", &bb_ps_e_out, "bb_ps_e/D");
  P->Branch("bb_ps_rowblk", &bb_ps_rowblk_out, "bb_ps_rowblk/D");
  P->Branch("bb_ps_colblk", &bb_ps_colblk_out, "bb_ps_colblk/D");
  P->Branch("bb_sh_e", &bb_sh_e_out, "bb_sh_e/D");
  P->Branch("bb_sh_rowblk", &bb_sh_rowblk_out, "bb_sh_rowblk/D");
  P->Branch("bb_sh_colblk", &bb_sh_colblk_out, "bb_sh_colblk/D");
  //P->Branch("bb_hodotdc_clus_tmean", &bb_hodotdc_clus_tmean_out, "bb_hodotdc_clus_tmean/D");
  //P->Branch("bb_grinch_tdc_clus_size", &bb_grinch_tdc_clus_size_out, "bb_grinch_tdc_clus_size/D");
  //P->Branch("bb_grinch_tdc_clus_trackindex", &bb_grinch_tdc_clus_trackindex_out, "bb_grinch_tdc_clus_trackindex/D");
  P->Branch("bb_gem_track_nhits", &bb_gem_track_nhits_out, "bb_gem_track_nhits/D");
  P->Branch("bb_gem_track_ngoodhits", &bb_gem_track_ngoodhits_out, "bb_gem_track_ngoodhits/D");
  P->Branch("bb_gem_track_chi2ndf", &bb_gem_track_chi2ndf_out, "bb_gem_track_chi2ndf/D");
  P->Branch("bb_etot_over_p", &bb_etot_over_p_out, "bb_etot_over_p/D");

  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = new TChain("T");

  cout << "Adding files " << fileset << " to chain..." << endl; 

  C->Add(fileset.c_str());

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
  // int Nhodotmean; 
  // double hodotmean[econst::maxclus];
  // std::vector<std::string> hodovar = {"clus.tmean","clus.tmean"};
  // std::vector<void*> hodovarlink = {&Nhodotmean,&hodotmean};
  // rvars::setbranch(C, "bb.hodotdc", hodovar, hodovarlink, 0); 

  // track branches
  double ntrack, p[econst::maxtrack],px[econst::maxtrack],py[econst::maxtrack],pz[econst::maxtrack],xtr[econst::maxtrack],ytr[econst::maxtrack],thtr[econst::maxtrack],phtr[econst::maxtrack];
  double vx[econst::maxtrack],vy[econst::maxtrack],vz[econst::maxtrack];
  double xtgt[econst::maxtrack],ytgt[econst::maxtrack],thtgt[econst::maxtrack],phtgt[econst::maxtrack];
  double r_x[econst::maxtrack],r_th[econst::maxtrack];
  std::vector<std::string> trvar = {"n","p","px","py","pz","x","y","th","ph","vx","vy","vz","tg_x","tg_y","tg_th","tg_ph","r_x","r_th"};
  std::vector<void*> trvarlink = {&ntrack,&p,&px,&py,&pz,&xtr,&ytr,&thtr,&phtr,&vx,&vy,&vz,&xtgt,&ytgt,&thtgt,&phtgt,&r_x,&r_th};
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
  double gemNhits, gemNgoodhits, gemChiSqr, eop;
  std::vector<std::string> miscbbvar = {"gem.track.nhits","gem.track.ngoodhits","gem.track.chi2ndf","etot_over_p"};
  std::vector<void*> miscbbvarlink = {&gemNhits,&gemNgoodhits,&gemChiSqr,&eop};
  rvars::setbranch(C, "bb", miscbbvar, miscbbvarlink);

  // Monte Carlo variables
  double mcsigma, mcomega, mcnucl;
  std::vector<std::string> mcvar = {"mc_sigma","mc_omega","mc_nucl"};
  std::vector<void*> mcvarlink = {&mcsigma,&mcomega,&mcnucl};
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
      
  //Event loop
  while (C->GetEntry(nevent++)) {
	
    std::cout << "Processing event " << nevent << " / " << nevents << ", total passed cuts " << npassed << "\r";
    std::cout.flush();

    ///////
    //Single-loop globalcut method. Save pass/fail for output tree.
    bool failedglobal = false;

    currenttreenum = C->GetTreeNumber();
    if( nevent == 1 || currenttreenum != treenum ){
      treenum = currenttreenum; 
      GlobalCut->UpdateFormulaLeaves();
    }

    failedglobal = GlobalCut->EvalInstance(0) == 0;

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
    double Q2, W2, nu, thNexp, pNexp, ebeam_o;
    ebeam_o = vars::ebeam_o( ebeam_c, etheta, "LD2" ); //Second energy correction accounting for energy loss leaving target
       
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

    if(!failedglobal && !failedW2){
      npassed++;
      if(npassed>maxEvents && limit_size){
	cout << "N Events exceeded maximum. Check limit_size if unexpected." << endl;
	return;
      }
    }

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

    //Determine rough PID
    bool in_p = util::Nspotcheck(dy,dx,dy0,dx0_p,dysig,dxsig_p,0);
    bool in_n = util::Nspotcheck(dy,dx,dy0,dx0_n,dysig,dxsig_n,0);

    int PID = -1;
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
      bool passedCoin = abs(atime_diff-coin_profile[1])<coin_sigma_factor*coin_profile[2];

      //Replicate the in-time algorithm with new cluster to be sorted later
      clone_cluster_intime.push_back(ce);
      if( !passedCoin )
	clone_cluster_intime[c] = 0;
	
      //Get score (no position info). Will be sorted later
      double cascore = util::assignScore( ce, atime_diff, hcalce[(int)hcalidx], coin_profile);
      clone_cluster_score.push_back(cascore);
	
    }//endloop over cluster elements

    int cidx_e = (int)hcalidx;
      
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
    double thetapq_bc_p = acos( protdir_bc.Dot( pNhat ) );
    double thetapq_bc_n = acos( neutdir_bc.Dot( pNhat ) );

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

    //Determine rough PID
    bool in_p_bc = util::Nspotcheck(dy_bestcluster,dx_bestcluster,dy0,dx0_p,dysig,dxsig_p,0);
    bool in_n_bc = util::Nspotcheck(dy_bestcluster,dx_bestcluster,dy0,dx0_n,dysig,dxsig_n,0);

    int PID_bc = -1;
    if(in_p_bc&&in_n_bc)
      PID_bc=0;
    else if(in_p_bc)
      PID_bc=1;
    else if(in_n_bc)
      PID_bc=2;

    //Find fiducial cut sigma factor
    std::pair<double, double> fiducial_factors = cut::findFidFailure(dxsig_p, 
								     dysig, 
								     xyhcalexp[0],
								     xyhcalexp[1], 
								     dx_del,
								     hcalaa);
	

    //cout << "Fiducial cut failed at x:y " << fiducial_factors.first << ":" << fiducial_factors.second << endl;

    //H-arm active area cuts (best_cluster)
    bool hcalON_bc = cut::hcalaaON(hcalcx[cidx_best],hcalcy[cidx_best],hcalaa);

    //H-arm active area cut (primary cluster)
    bool hcalON = cut::hcalaaON(hcalx,hcaly,hcalaa);

    ////////////////////
    //coincidence time BBCal/HCal cut
    bool failedwidecoin_bc = abs( hatime_bestcluster - atimediff0 ) > coin_sigma_factor*atimediffsig;

    //caculate final weight for this event
    Double_t mcweight = mcsigma/mcomega;

    //Fill new output tree  
    //Universals
    event_out = (int)gevnum;
    trig_out = (int)trigbits;
    xexp_out = xyhcalexp[0];
    yexp_out = xyhcalexp[1];
    W2_out = W2;
    Q2_out = Q2;
    nu_out = nu;
    precon_out = precon;
    failedglobal_out = failedglobal;
    failedW2_out = failedW2;

    //MC data
    mcsigma_out = mcsigma;
    mcomega_out = mcomega;
    mcweight_out = mcweight;
    mcnucl_out = mcnucl;

    //Fiducial slices
    fiducial_sig_x_out = fiducial_factors.first;
    fiducial_sig_y_out = fiducial_factors.second;

    //Primary Cluster
    dx_out = dx;
    dy_out = dy;
    coin_out = pclus_diff;
    thetapq_pout = thetapq_p;
    thetapq_nout = thetapq_n;
    nucleon_spot_out = PID;
    hcalnblk_out = pclus_nblk;
    hcalon_out = hcalON;
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
    nucleon_bc_out = PID_bc;
    hcalnblk_bc_out = nblk_bestcluster;
    hcalon_bc_out = hcalON_bc;
    hcalpid_bc_out = pblkid_bestcluster;
    hcalx_bc_out = x_bestcluster;
    hcaly_bc_out = y_bestcluster;
    hcale_bc_out = ce_bestcluster;

    //Fill old output tree
    bb_tr_p_out = p[0];
    bb_tr_n_out = ntrack;
    bb_tr_vz_out = vz[0];
    bb_tr_th_out = thtr[0];
    bb_tr_ph_out = phtr[0];
    bb_tr_r_x_out = r_x[0];
    bb_tr_r_th_out = r_th[0];
    bb_ps_e_out = ePS;
    bb_ps_rowblk_out = rblkPS;
    bb_ps_colblk_out = cblkPS;
    bb_sh_e_out = eSH;
    bb_sh_rowblk_out = rblkSH;
    bb_sh_colblk_out = cblkSH;
    //bb_hodotdc_clus_tmean_out = hodotmean[0];
    //bb_grinch_tdc_clus_size_out = grinchClusSize;
    //bb_grinch_tdc_clus_trackindex_out = grinchClusTrIndex;
    bb_gem_track_nhits_out = gemNhits;
    bb_gem_track_ngoodhits_out = gemNgoodhits;
    bb_gem_track_chi2ndf_out = gemChiSqr;
    bb_etot_over_p_out = eop;

    P->Fill();	

  }//end event loop

  fout->Write();

  std::cout << std::endl << "Barebones mcbg parsing complete. Output written to " << parse_path << std::endl << std::endl;

  st->Stop();

  // Send time efficiency report to console
  std::cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << std::endl;

}
