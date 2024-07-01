//sseeds 04.20.23 - Script to extract single-channel best possible timing resolutions by kinematic/field

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

//Passing kine==-1 will run all kinematics, pass is replay pass, epm is e' momentum calculation method
void htres( Int_t kine=8, Int_t epm=3, bool waveform=false, Int_t pass=1 )
{ //main

  Double_t adcbinw = econst::hcaladc_binw; //adc sample bin width (4*40=160ns over whole waveform)

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // reading config file
  JSONManager *jmgr = new JSONManager("../../config/shtres.json");
  std::string rootfile_dir;
  if(waveform)
    rootfile_dir = jmgr->GetValueFromKey_str("expert_rootfile_dir");
  else
    rootfile_dir = jmgr->GetValueFromSubKey_str("general_rootfile_dir",Form("sbs%d",kine));

  //define structure to encode parameters for continuous ep fit tdc/adct corrections
  typedef struct{
    vector<Int_t> kIdx;
    Int_t nfset_lh2; //total number of different field settings for lh2
    Int_t nfset_ld2; //total number of different field settings for ld2
    vector<Int_t> fset_lh2; //all field settings for lh2 (percent, -1 indicates no setting)
    vector<Int_t> fset_ld2; //all field settings for ld2 (percent, -1 indicates no setting)
    Int_t TOFfitp_mset; //field settings with TOF vs nucleon p available
    vector<Double_t> TOFfitp_p; //TOF vs proton p fit parameters
    vector<Double_t> TOFfitp_n; //TOF vs neutron p fit parameters
    Double_t TOFfitp_ulim; //upper nucleon momentum limit on fits, beyond which fits diverge quickly
  } SBSSET;

  SBSSET fitparams[econst::nkine]; //One structure per kinematic

  Int_t gIdx; //get the global index for this kinematic
  for( Int_t k=0; k<econst::nkine; k++ ){
    jmgr->GetVectorFromKey<Int_t>("kIdx",fitparams[k].kIdx);
    Int_t ki = fitparams[k].kIdx[k];

    fitparams[k].nfset_lh2 = jmgr->GetValueFromSubKey<Int_t>("nfset_lh2",Form("sbs%d",ki));    
    fitparams[k].nfset_ld2 = jmgr->GetValueFromSubKey<Int_t>("nfset_ld2",Form("sbs%d",ki));
    jmgr->GetVectorFromSubKey<Int_t>("fset_lh2",Form("sbs%d",ki),fitparams[k].fset_lh2);
    jmgr->GetVectorFromSubKey<Int_t>("fset_ld2",Form("sbs%d",ki),fitparams[k].fset_ld2);
    fitparams[k].TOFfitp_mset = jmgr->GetValueFromSubKey<Int_t>("TOFfitp_mset",Form("sbs%d",ki));
    jmgr->GetVectorFromSubKey<Double_t>("TOFfitp_p",Form("sbs%d",ki),fitparams[k].TOFfitp_p);
    jmgr->GetVectorFromSubKey<Double_t>("TOFfitp_n",Form("sbs%d",ki),fitparams[k].TOFfitp_n);
    fitparams[k].TOFfitp_ulim = jmgr->GetValueFromSubKey<Double_t>("TOFfitp_ulim",Form("sbs%d",ki));

    if( fitparams[0].kIdx[k]==kine )
      gIdx=k;
  }

  //Read in timewalk parameters
  typedef struct{
    vector<Double_t> tdcp;
    vector<Double_t> adcp;
  } TW;

  TW timewalk[3]; //three timewalk parameters of the form p0*exp(-p1*x[0])+p2

  for( Int_t p=0; p<3; p++ ){
    std::string tdcpath = Form("params/TW/TDC/tdctwP%d_sbs%d.txt",p,kine);
    util::readParam(tdcpath,timewalk[p].tdcp);
    std::string adctpath = Form("params/TW/ADC/adcttwP%d_sbs%d.txt",p,kine);
    util::readParam(adctpath,timewalk[p].adcp);
  }

  //Read in gain parameters for more accurate timewalk corrections
  vector<Double_t> oldadcgain;
  vector<Double_t> adcgain;
  std::string oldadcgainpath = "params/gain/oldcoeff.txt";
  std::string adcgainpath = Form("params/gain/coeff_sbs%d.txt",kine);
  std::string gmoniker = "sbs.hcal.adc.gain";
  util::readParam(oldadcgainpath,gmoniker,econst::hcalchan,oldadcgain);
  util::readParam(adcgainpath,gmoniker,econst::hcalchan,adcgain);

  vector<Int_t> lh2r; //jmgr adds (not replace) when called
  jmgr->GetVectorFromSubKey<Int_t>( "lh2runs", "sbs4_runs", lh2r );
  jmgr->GetVectorFromSubKey<Int_t>( "lh2runs", "sbs9_runs", lh2r );

  vector<Int_t> ld2r;
  jmgr->GetVectorFromSubKey<Int_t>( "ld2runs", "sbs8_runs", ld2r );
  jmgr->GetVectorFromSubKey<Int_t>( "ld2runs", "sbs11_runs", ld2r );

  //set up default parameters for all analysis
  std::string runsheet_dir = "/w/halla-scshelf2102/sbs/seeds/ana/data"; //unique to my environment for now
  Int_t nruns = -1; //Always analyze all available runs for this application. If not, need to pass nruns for lh2 and ld2 separately.
  Int_t nruns_h = -1;
  Int_t nruns_d = -1;
  Int_t verb = 0; //Don't print diagnostic info by default

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<crun> crunh; 
  vector<crun> crund; 
  util::ReadRunList(runsheet_dir,nruns_h,kine,"LH2",pass,verb,crunh); //modifies nruns to be very large when -1
  util::ReadRunList(runsheet_dir,nruns_d,kine,"LD2",pass,verb,crund); //modifies nruns to be very large when -1
  
  //set up output files
  TFile *fout = new TFile( Form("outfiles/htresout_sbs%d_epm%d.root",kine,epm), "RECREATE" );

  //set up waveform histograms
  TH1D *landhist[econst::hcalrow][econst::hcalcol]; //For landau fits
  TH1D *sghist[econst::hcalrow][econst::hcalcol]; //For skewed gaussian fits
  TH1D *lghist[econst::hcalrow][econst::hcalcol]; //For landau gaussian convolution fits
  for(Int_t r = 0; r < econst::hcalrow; r++) {
    for(Int_t c = 0; c < econst::hcalcol; c++) {
      landhist[r][c] = util::hhsamps(r,c,econst::maxsamp);
      sghist[r][c] = (TH1D*)landhist[r][c]->Clone(Form("sghist %d-%d",r,c));
      lghist[r][c] = (TH1D*)landhist[r][c]->Clone(Form("lghist %d-%d",r,c));
    }
  }

  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;
  
  // create output tree
  TTree *P = new TTree("P","Analysis Tree"); 

  // timing
  Double_t pblkid_out; //hcal primary cluster, primary block id
  Double_t latime_out; //hcal primary cluster landau fit adc time
  Double_t sgatime_out; //hcal primary cluster skewed gaussian fit adc time
  Double_t batime_out; //hcal primary cluster "best" fit adc time
  Double_t atime_out; //hcal primary cluster adc time, tree
  Double_t tdc_out; //hcal primary cluster tdc time, tree

  Double_t bcatime_dxdy_out; //hcal dxdy best cluster adc time, tree
  Double_t bcatime_atime_out; //hcal atime best cluster adc time, tree
  Double_t bcatime_thpq_out; //hcal thetapq best cluster adc time, tree
  Double_t bctdc_dxdy_out; //hcal dxdy best cluster adc time, tree
  Double_t bctdc_atime_out; //hcal atime best cluster adc time, tree
  Double_t bctdc_thpq_out; //hcal thetapq best cluster adc time, tree
  Double_t bcpblkid_dxdy_out; //hcal dxdy best cluster, primary block id
  Double_t bcpblkid_atime_out; //hcal atime best cluster, primary block id
  Double_t bcpblkid_thpq_out; //hcal thetapq best cluster, primary block id

  //"best" timing (atime cut around 3sig, then min thetapq)
  Double_t btdc_out;
  Double_t badct_out;
  Double_t tofcorr_out;
  Double_t bclustwcorr_out;
  Double_t bclusatwcorr_out;
  Double_t bcpblkid_cuts_out; //channel of best cluster primary block

  //physics
  Double_t dx_out; //hcal primary cluster dx
  Double_t dy_out; //hcal primary cluster dy 
  Double_t W2_out;
  Double_t Q2_out;
  Double_t pN_out;
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

  //cluster tree vars
  Double_t cpblkid_out[econst::maxclus];
  Double_t chatime_out[econst::maxclus];
  Double_t chtdc_out[econst::maxclus];
  Double_t cthetapq_p_out[econst::maxclus];
  Double_t cthetapq_n_out[econst::maxclus];
  Double_t cdx_out[econst::maxclus];
  Double_t cdy_out[econst::maxclus];
  Int_t cpid_out[econst::maxclus];

  // set output tree branches
  P->Branch( "pblkid", &pblkid_out, "pblkid/D" );
  P->Branch( "latime", &latime_out, "latime/D" );
  P->Branch( "sgatime", &sgatime_out, "sgatime/D" );
  P->Branch( "batime", &batime_out, "batime/D" );
  P->Branch( "atime", &atime_out, "atime/D" );

  P->Branch( "bcatime_dxdy", &bcatime_dxdy_out, "bcatime_dxdy/D" );
  P->Branch( "bcatime_atime", &bcatime_atime_out, "bcatime_atime/D" );
  P->Branch( "bcatime_thpq", &bcatime_thpq_out, "bcatime_thpq/D" );
  P->Branch( "bctdc_dxdy", &bctdc_dxdy_out, "bctdc_dxdy/D" );
  P->Branch( "bctdc_atime", &bctdc_atime_out, "bctdc_atime/D" );
  P->Branch( "bctdc_thpq", &bctdc_thpq_out, "bctdc_thpq/D" );
  P->Branch( "bcpblkid_dxdy", &bcpblkid_dxdy_out, "bcpblkid_dxdy/D" );
  P->Branch( "bcpblkid_atime", &bcpblkid_atime_out, "bcpblkid_atime/D" );
  P->Branch( "bcpblkid_thpq", &bcpblkid_thpq_out, "bcpblkid_thpq/D" );
  
  P->Branch( "btdc", &btdc_out, "btdc/D" );
  P->Branch( "badct", &badct_out, "badct/D" );
  P->Branch( "tofcorr", &tofcorr_out, "tofcorr/D" );
  P->Branch( "bclustwcorr", &bclustwcorr_out, "bclustwcorr/D" );
  P->Branch( "bclusatwcorr", &bclusatwcorr_out, "bclusatwcorr/D" );
  P->Branch( "bcpblkid_cuts", &bcpblkid_cuts_out, "bcpblkid_cuts/D" );

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

  P->Branch( "cpblkid", &cpblkid_out, Form("cpblkid[%d]/D",econst::maxchan) );
  P->Branch( "chatime", &chatime_out, Form("chatime[%d]/D",econst::maxchan) );
  P->Branch( "chtdc", &chtdc_out, Form("chtdc[%d]/D",econst::maxchan) );
  P->Branch( "cthetapq_p", &cthetapq_p_out, Form("cthetapq_p[%d]/D",econst::maxchan) );
  P->Branch( "cthetapq_n", &cthetapq_n_out, Form("cthetapq_n[%d]/D",econst::maxchan) );
  P->Branch( "cdx", &cdx_out, Form("cdx[%d]/D",econst::maxchan) );
  P->Branch( "cdy", &cdy_out, Form("cdy[%d]/D",econst::maxchan) );
  P->Branch( "cpid", &cpid_out, Form("cpid[%d]/I",econst::maxchan) );

  // setup reporting indices
  Int_t curmag = -1;
  std::string roottar = "";

  for ( Int_t t=0; t<2; t++ ){ //loop over targets
    // t==0, lh2; t==1, ld2
    std::string rootfile_fdir;

    if( t==0 ){
      roottar = "LH2";
      nruns = nruns_h;
    }else{
      roottar = "LD2";
      nruns = nruns_d;
    }

    if( !waveform )
      rootfile_fdir = rootfile_dir + Form("/%s/rootfiles",roottar.c_str());
    else
      rootfile_fdir = rootfile_dir;

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
	rfname = rootfile_fdir + Form("/*%d*",crunh[irun].runnum);

	//selects from the runsheet only those runs for which expert replayes exist
	bool skip = true;
	for( Int_t el=0; el<lh2r.size(); el++ ){     
	  if( runnum==lh2r[el] ) skip=false;
	}
	if( skip==true && waveform ) continue;
      }
      if( t==1 ){
	runnum = crund[irun].runnum;
	mag = crund[irun].sbsmag / 21;
	ebeam = crund[irun].ebeam;
	conf = crund[irun].sbsconf; //get config of current run
	targ = crund[irun].target;
	rfname = rootfile_fdir + Form("/*%d*",crund[irun].runnum);

	//selects from the runsheet only those runs for which expert replayes exist
	bool skip = true;
	for( Int_t el=0; el<ld2r.size(); el++ ){     
	  if( runnum==ld2r[el] ) skip=false;
	}
	if( skip==true && waveform ) continue;
      }      

      if( conf!=kine && kine!=-1) continue; //do not proceed if kinematic is constrained
      if( !waveform && mag!=fitparams[gIdx].TOFfitp_mset ) continue; //continue if TOF corrections not available

      std::cout << "Analyzing run " << runnum << ".." << std::endl;

      //set up configuration and tune objects to load analysis parameters
      SBSconfig config(kine,mag);

      //Obtain configuration pars from config file
      Double_t hcaltheta = config.GetHCALtheta_rad();
      Double_t hcaldist = config.GetHCALdist();
      Double_t sbsdist = config.GetSBSdist();
      Double_t bbthr = config.GetBBtheta_rad(); //in radians

      SBStune tune(kine,mag);
    
      //if( targ.compare(curtar)!=0 || mag!=curmag ){
      if( mag!=curmag ){
	std::cout << "Settings change.." << endl;
	cout << config;
	cout << tune;
	//curtar = targ;
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

      std::cout << "Added " << rfname << " to the chain with " << C->GetEntries() << " entries.." << std::endl;

      // setting up ROOT tree branch addresses
      C->SetBranchStatus("*",0);    

      // HCal general
      Double_t hcalid, hcale, hcalx, hcaly, hcalr, hcalc, hcala, hcalamp, hcaltdc, hcalatime;
      std::vector<std::string> hcalvar = {"idblk","e","x","y","rowblk","colblk","a_p","a_amp_p","tdctimeblk","atimeblk"};
      std::vector<void*> hcalvarlink = {&hcalid,&hcale,&hcalx,&hcaly,&hcalr,&hcalc,&hcala,&hcalamp,&hcaltdc,&hcalatime};
      rvars::setbranch(C, "sbs.hcal", hcalvar, hcalvarlink);

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
	Double_t Q2, W2, nu, thNexp, pNexp, pNcalc, pp, pn;
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
	  pNexp = vars::pN_expect( nu, nucleon );
	  thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
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
	  W2 = vars::W2( pbeam.E(), pe.E(), W2, nucleon );
	}else if( epm==4 ){
	  //v4
	  nu = pbeam.E() - pe.E();
	  pNexp = vars::pN_expect( nu, nucleon );
	  thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	  pNhat = vars::pNhat_track( thNexp, phNexp );
	  pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	  Q2 = vars::Q2( pbeam.E(), pe.E(), etheta );
	  W2 = vars::W2( pbeam.E(), pe.E(), W2, nucleon );
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
	pp = sqrt(pow(nu,2)+2.*physconst::Mp*nu); //momentum of proton
	pn = sqrt(pow(nu,2)+2.*physconst::Mn*nu); //momentum of proton

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
	if( (is_p && !is_n) || targ.compare("LH2")==0 ){
	  pid=1;
	  pNcalc = vars::pN_expect( nu, "p" );
	}else if( is_n && !is_p ){
	  pid=2;
	  pNcalc = vars::pN_expect( nu, "n" );
	}else if( is_p && is_n ){
	  pid=0;
	  pNcalc = vars::pN_expect( nu, "n" );
	}else{
	  pid=-1;
	  pNcalc = pNexp;
	}

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
	//extract atime with waveforms from primary cluster primary block
	Int_t pblkid = hcalcbid[0];
	Int_t r,c,idx,n,sub;
	Double_t adc[econst::hcalrow][econst::hcalcol];
	Double_t amp[econst::hcalrow][econst::hcalcol];
	Float_t peak[econst::hcalrow][econst::hcalcol];

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
	
	//Primary cluster element
	Int_t pblkrow = (int)hcalcbid[0]/econst::hcalcol;
	Int_t pblkcol = (int)hcalcbid[0]%econst::hcalcol;

	//atime fit parameters
	Double_t sgmpv = -1000.;
	Double_t lmpv = -1000.;
	Double_t brising_edge = -1000.;
	Double_t lrising_edge = -1000.;
	Double_t lrising_edge_v2 = -1000.;
	Double_t sgrising_edge = -1000.;
	Double_t sgrising_edge_v2 = -1000.;
	Double_t lgrising_edge = -1000.;

	//bool pblkfill = false;

	//Loop over hcal samples for primary cluster primary block waveform extraction and analysis
	for(Int_t m = 0; m < Nssave; m++) {

	  if( !waveform ) continue;

	  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
	  //gSystem->RedirectOutput("/dev/null");

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

	  n = nsamps[m];
	  
	  adc[r][c] = a_p[m];
	  amp[r][c] = a_amp_p[m];

	  bool saturated = amp[r][c]>3500; //in RAU
	  bool mansat = false;
	  bool negped = adc[r][c]<-5; //Should be fixed and never occur
	  bool skip = true;
	  for( int b = 0; b < Nhcalcbid; b++ ){
	    if( chan==(hcalcbid[0]-1) ){ 
	      skip=false;
	    }
	  }

	  if( skip==true ) continue; //Only continue if the current channel is the primary blk channel

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

	  if(!displayed) {
	    std::cerr << "Skipping empty module: " << m << std::endl;
	    for(Int_t s = 0;  s < econst::maxsamp; s++) {
	      landhist[r][c]->SetBinContent(s+1,-404);
	      sghist[r][c]->SetBinContent(s+1,-404);
	      lghist[r][c]->SetBinContent(s+1,-404);
	    }
	  }
	  
	  //landau fit
	  TF1 *lfit;
	  landhist[r][c]->Fit("landau", "Q", "", econst::minsamp, econst::maxsamp);
	  lfit = landhist[r][c]->GetFunction("landau");
	  lmpv = lfit->GetMaximumX();
	  Double_t lch2 = lfit->GetChisquare();
	  Double_t lmean = lfit->GetParameter(1);
	  Double_t lsig = lfit->GetParameter(2);
	  lrising_edge = lmean-3*lsig; //Not as effective as brute force
	  lrising_edge_v2 = lfit->GetX(0.00001,econst::minsamp+1,econst::maxsamp-1);

	  //sgfit
	  TF1 *sgfit = new TF1(Form("sg r:%d c:%d",r,c),fits::g_sgfit, econst::minsamp, econst::maxsamp-1, 4);
	  Double_t sgamean = sghist[r][c]->GetMean();
	  Double_t sgasig = (econst::maxsamp-1)*0.15;
	  sgfit->SetParameters(sghist[r][c]->GetMaximum(),sgamean,2.,1.);
	  sgfit->SetParLimits(0,0.005,2.0);
	  sgfit->SetParLimits(1,econst::minsamp,econst::maxsamp-1);
	  sgfit->SetParLimits(2,4.,sgasig);
	  sgfit->SetParLimits(3,5.,13.);
	  sghist[r][c]->Fit(sgfit,"RQ","",econst::minsamp, econst::maxsamp);
	  sgmpv = sgfit->GetMaximumX();
	  Double_t sgsig = sgfit->GetParameter(2);
	  Double_t sskew = sgfit->GetParameter(3);
	  sgrising_edge = sgmpv+sgsig*sqrt(2)*TMath::ErfInverse(-sskew); //Not as effective as brute force
	  sgrising_edge_v2 = sgfit->GetX(0.00001,econst::minsamp+1,econst::maxsamp-1);
	  Double_t sgch2 = sgfit->GetChisquare();

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

	  //gaussian landau convolution fit. This takes very long. Will leave commented out.
	  /*
	  TF1 *glf;
	  Double_t land_mean=sgamean, land_sig=sgasig; 
	  Double_t gaus_amp=lghist[r][c]->GetMaximum(), gaus_mean=sgamean, gaus_sig=sgasig;
	  Int_t no_FFT_pts=1000;

	  fits::g_conv_gausland( lghist[r][c], glf, land_mean, land_sig, gaus_amp, gaus_mean, gaus_sig, no_FFT_pts );
	  */

	  //diagnostic plots
	  /*
	  if( pltcntr<20 && nevent%5==0 && pblkfill==false){

	    lghist[r][c]->SetTitle(Form("lgconv ev:%ld r:%d c:%d at0:%f re:%0.2f p1:%0.2f p2:%0.2f p3:%0.2f p4:%0.2f p5:%0.2f chi2:%0.2f",nevent,r,c,hcalatime/4.,lgrising_edge,lgp1,lgp2,lgp3,lgp4,lgp5,lgch2));
	    lghist[r][c]->SetName(Form("lgev%ld",nevent));
	    lghist[r][c]->Write();

	    //if( lch2 < sgch2 ){
	      Double_t lmean = lfit->GetParameter(1);
	      landhist[r][c]->SetTitle(Form("landau ev:%ld r:%d c:%d at0:%f re:%0.2f rev2:%0.2f adc:%0.2f amp:%0.2f mpv:%0.2f mean:%0.2f sig:%0.2f chi2:%0.2f",nevent,r,c,hcalatime/4.,lrising_edge,lrising_edge_v2,adc[r][c],amp[r][c],lmpv,lmean,lsig,lch2));
	      landhist[r][c]->SetName(Form("lev%ld",nevent));
	      landhist[r][c]->Write();
	      //}else{
	      Double_t sgmean = sgfit->GetParameter(1);
	      sghist[r][c]->SetTitle(Form("sgaus ev:%ld r:%d c:%d at0:%f re:%0.2f rev2:%0.2f adc:%0.2f amp:%0.2f mpv:%0.2f mean:%0.2f sig:%0.2f skew:%0.2f chi2:%0.2f",nevent,r,c,hcalatime/4.,sgrising_edge,sgrising_edge_v2,adc[r][c],amp[r][c],sgmpv,sgmean,sgsig,sskew,sgch2));
	      sghist[r][c]->SetName(Form("sgev%ld",nevent));
	      sghist[r][c]->Write();
	      //}
	    pltcntr++;
	    pblkfill=true;
	  }
	  */

	} //end loop over hcal samples

	//Fill all cluster primary block physics output - this is inefficient, should get primary cluster info from here first
	Int_t cidx_atime = 0;
	Double_t c_atimediff = 1000.;
	Int_t cidx_thetapq = 0;
	Double_t c_thetapqdiff = 1000.;
	Int_t cidx_dxdy = 0;
	Double_t c_dxdydiff = 1000.;
	Double_t tpq_best;

	Int_t cidx_cuts = 0;
	Double_t c_cutsdiff = 1000.;
	bool cuts_skip = true;

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
	    c_atimediff = abs(hcalcatime[c]-atime0);
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

	  //get cluster which passes atime cut first, then minimizes thetapq ("best")
	  if( abs(hcalcatime[c]-atime0)<3*atimesig ){
	    if( cthetapq_p<c_cutsdiff ){
	      cuts_skip = false;
	      cidx_cuts = c;
	      c_cutsdiff = cthetapq_p;
	    }
	  }

	  cpblkid_out[c] = hcalcid[c];
	  chatime_out[c] = hcalcatime[c];
	  chtdc_out[c] = hcalctdctime[c];
	  cthetapq_p_out[c] = cthetapq_p;
	  cthetapq_n_out[c] = cthetapq_n;
	  cdx_out[c] = cdx;
	  cdy_out[c] = cdy;
	  
	  cpid_out[c] = cpid;
	}

	//Calculate dx, dy, thetapq from the best cluster for later use (USING dxdy for now)
	Double_t dx_bestcluster = hcalcx[cidx_cuts] - xyhcalexp[0];
	Double_t dy_bestcluster = hcalcy[cidx_cuts] - xyhcalexp[1];
	Double_t hatime_bestcluster = hcalcatime[cidx_cuts];
	Int_t bcluschan = hcalcbid[cidx_cuts];

	//For best clusters, apply time of flight to all
	Double_t protmom = pp;
	Double_t neutmom = pn;
	Double_t p0_p = fitparams[gIdx].TOFfitp_p[0];
	Double_t p1_p = fitparams[gIdx].TOFfitp_p[1];
	Double_t p2_p = fitparams[gIdx].TOFfitp_p[2];
	Double_t p3_p = fitparams[gIdx].TOFfitp_p[3];
	Double_t p0_n = fitparams[gIdx].TOFfitp_n[0];
	Double_t p1_n = fitparams[gIdx].TOFfitp_n[1];
	Double_t p2_n = fitparams[gIdx].TOFfitp_n[2];
	Double_t p3_n = fitparams[gIdx].TOFfitp_n[3];

	//time of flight proton correction from momentum fit
	Double_t tofcorr = 0.;
	if( protmom>fitparams[gIdx].TOFfitp_ulim ) protmom = fitparams[gIdx].TOFfitp_ulim; //do not pass corrections for momenta > fit limit in MC
	Double_t tofpcorr_p = p0_p+p1_p*protmom+p2_p*pow(protmom,2)+p3_p*pow(protmom,3);

	//time of flight neutron correction from momentum fit
	if( neutmom>fitparams[gIdx].TOFfitp_ulim ) neutmom = fitparams[gIdx].TOFfitp_ulim;
	Double_t tofpcorr_n = p0_n+p1_n*neutmom+p2_n*pow(neutmom,2)+p3_n*pow(neutmom,3);

	if( pid==1 )
	  tofcorr = tofpcorr_p;
	else if( pid==2 || pid==0)
	  tofcorr = tofpcorr_n;
	else if( pid==-1 )
	  tofcorr = 0.;

	//timewalk corrections - leave off overall offset corrections (p2)
	Double_t bcE = hcalcbe[cidx_cuts]/oldadcgain[bcluschan]*adcgain[bcluschan]; //best here is closest to expected adctime. Correct with pass2 gain coefficients not already included in pass0/1 replays
	Double_t twp0tdc = timewalk[0].tdcp[bcluschan];
	Double_t twp1tdc = timewalk[1].tdcp[bcluschan];
	Double_t twp0adc = timewalk[0].adcp[bcluschan];
	Double_t twp1adc = timewalk[1].adcp[bcluschan];

	Double_t twcorr = twp0tdc*exp(-twp1tdc*bcE); //tdc
	Double_t atwcorr = twp0adc*exp(-twp1adc*bcE); //adct

	//trigger corrections (hodoscope primary (idx 0) cluster meantime)
	Double_t trcorr = hodotmean[0];

	pblkid_out = (double)hcalcbid[0];
	latime_out = lrising_edge_v2*adcbinw;
	sgatime_out = sgrising_edge_v2*adcbinw;
	batime_out = brising_edge;
	atime_out = hcalcbatime[0];
	tdc_out = hcalcbtdctime[0];

	//mark events where no "best" cluster was available (failed atime cut or no passable thetapq )
	Double_t skipmarker = 0.;
	if( cuts_skip )
	  skipmarker = -1000.;
	btdc_out = hcalctdctime[cidx_cuts]-trcorr-twcorr-tofcorr+skipmarker;
	badct_out = hcalcatime[cidx_cuts]-trcorr-tofcorr+skipmarker; //not clear that tw is effective with adct

	tofcorr_out = tofcorr;
	bclustwcorr_out = twcorr;
	bclusatwcorr_out = atwcorr;
	bcpblkid_cuts_out = (double)hcalcbid[cidx_cuts];

	bcatime_dxdy_out = hcalcatime[cidx_dxdy];
	bcatime_atime_out = hcalcatime[cidx_atime];
	bcatime_thpq_out = hcalcatime[cidx_thetapq];
	bctdc_dxdy_out = hcalctdctime[cidx_dxdy];
	bctdc_atime_out = hcalctdctime[cidx_atime];
	bctdc_thpq_out = hcalctdctime[cidx_thetapq];
	bcpblkid_dxdy_out = (double)hcalcbid[cidx_dxdy];
	bcpblkid_atime_out = (double)hcalcbid[cidx_atime];
	bcpblkid_thpq_out = (double)hcalcbid[cidx_thetapq];

	dx_out = dx;
	dy_out = dy;
	W2_out = W2;
	Q2_out = Q2;
	pN_out = pNcalc;
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
