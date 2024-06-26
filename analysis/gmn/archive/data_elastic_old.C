//SSeeds 9.1.23 Script to extract the p/n yield ratio via data/MC comparison. Intent is to loop over all data from a given set/field setting, apply elastic cuts, and build a dx histogram. From this point, read in MC (with RC) distributions for quasi-elastic protons and neutrons (independently), fit these distributions with several distributions to get a best fit, then compose a sum of these fits allowing for a scaling parameter which maps to the total MC yields. Next, do the same with g4sbs inelastic generator and add this to the total fit function (sum) and apply the now three floating parameter fit to the data and check residuals and chi-square. Varying the fit function may yield better results, but the ratio will be extracted. See gmn_calc.C for application of this ratio to gmn extraction.
//sseeds update 9.23.23 Parsed cluster sorting and e' momentum strategies. Renamed and refocused script to produce output histograms for later analysis.
//sseeds update 9.24.23 Debug and remove memory leaks

#include <vector>
#include <iostream>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "TMatrixD.h"
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

const int coinSigFac = 5;
const double dx_fitlow = -0.3;
const double dx_fithigh = 0.3;

//assign score to potential cluster based on probability density of gaussian fit to atime and maximum energy, exclude position info. Here, gfit is a vector with coin atime fit parameters, val_2 is the cluster coin atime, val_1 is the cluster energy, and max1 is the max cluster energy for the event
double assignScore( double val_1, double val_2, double max1, const std::vector<double> &gfit ) {

  if( gfit.size()!=3 ){
    cout << "ERROR: size of dxfit vector not equal to expected number of gaussian parameters (3)." << endl;
    return 0.;
  }

  // Build a TF1 with fit provided fit parameters
  TF1 *gauss = new TF1("gauss", "gaus", 0, 100);
  gauss->SetParameter(0,gfit[0]);
  gauss->SetParameter(1,gfit[1]);
  gauss->SetParameter(2,gfit[2]);
    
  double x_value = val_2;
  double density = gauss->Eval(x_value);

  // Compute the score based on val 2
  double score = density / gauss->GetParameter(0);  // Normalize by the peak of the Gaussian

  // Include a product component based on val 1 (linear weight)
  score *= val_1 / max1;

  if(score==0){
    //cout << "WARNING: score is zero. val_1=" << val_1 << ", val_2= " << val_2 << " density= " << density << " 1 comp= " << val_1 / max1 << endl;
    score=1e-38; //write very small number to score to exclude it without breaking the sorting later
  }

  return score;
}


//cluster non-zero sorting. if two single valued vectors equal to zero are returned, no non-zero values were found
std::pair<std::vector<double>, std::vector<size_t>> sortNonZero(const std::vector<double> &arr) {
  // Pair of value and its index in the original array
  std::vector<std::pair<double, size_t>> nonZeroPairs;

  // Populate the pairs
  bool empty = true;
  for (size_t i = 0; i < arr.size(); ++i) {
    if (arr[i] != 0.0) {
      nonZeroPairs.push_back({arr[i], i});
      empty = false;
    }
  }

  //prevent the seg fault!
  if( empty )
    nonZeroPairs.push_back({0.,0.});

  // Sort the pairs based on value in descending order
  std::sort(nonZeroPairs.begin(), nonZeroPairs.end(), [](const auto& a, const auto& b) {
      return a.first > b.first;
    });

  std::vector<double> sortedValues;
  std::vector<size_t> sortedIndices;

  for (const auto& pair : nonZeroPairs) {
    sortedValues.push_back(pair.first);
    sortedIndices.push_back(pair.second);
  }

  return {sortedValues, sortedIndices};
}

//MAIN
void data_elastic( Int_t kine=8, Int_t magset=70 )
{

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/gmn.json");

  std::string rootfile_dir = jmgr->GetValueFromSubKey_str( "rootfile_dir", Form("sbs%d",kine) );
  Double_t W2fitmax = jmgr->GetValueFromSubKey<Double_t>( "W2fitmax", Form("sbs%d",kine) );
  Double_t binfac = jmgr->GetValueFromSubKey<Double_t>( "binfac", Form("sbs%d",kine) );
  Double_t hbinfac = jmgr->GetValueFromSubKey<Double_t>( "hbinfac", Form("sbs%d",kine) ); //gets to 7cm HCal pos res
  std::string target = jmgr->GetValueFromKey_str( "target" ); //Must always be LD2 for GMn analysis
  Int_t epm = jmgr->GetValueFromSubKey<Int_t>( "epm", Form("sbs%d",kine) );
  Double_t minE = jmgr->GetValueFromSubKey<Double_t>( "minE", Form("sbs%d",kine) );
  Int_t cluster_idx = jmgr->GetValueFromSubKey<Int_t>( "cluster_idx", Form("sbs%d",kine) );
  Int_t pass = jmgr->GetValueFromSubKey<Int_t>( "pass", Form("sbs%d",kine) );
  Double_t ebound_l = jmgr->GetValueFromSubKey<Double_t>( "ebound_l", Form("sbs%d",kine) );
  Double_t ebound_h = jmgr->GetValueFromSubKey<Double_t>( "ebound_h", Form("sbs%d",kine) );
  vector<Double_t> coin_profile;
  jmgr->GetVectorFromSubKey<Double_t>("coin_profile",Form("sbs%d",kine),coin_profile);

  if( pass>1 ){
    std::cout << "As of 3.31.23, the highest GMn replay pass is 1. Enter a valid pass." << endl;
    return;
  }

  //set up default parameters for all analysis
  std::string runsheet_dir = "/w/halla-scshelf2102/sbs/seeds/ana/data"; //unique to my environment for now
  Int_t nruns = -1; //Always analyze all available runs
  Int_t verb = 0; //Don't print diagnostic info by default
  Double_t hcalfit_l = econst::hcalposXi_mc; //lower fit/bin limit for hcal dx plots (m)
  Double_t hcalfit_h = econst::hcalposXf_mc; //upper fit/bin limit for hcal dx plots (m)
  Double_t harmrange = econst::hcal_vrange; //Full range of hcal dx plots (m)

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<crun> corun; 
  util::ReadRunList(runsheet_dir,nruns,kine,target.c_str(),pass,verb,corun); //modifies nruns to be very large when -1

  //set up output files
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fout_path = outdir_path + Form("/gmn_analysis/gmn_elastic_out_sbs%d_mag%d_clusteridx%d_epm%d.root",kine,magset,cluster_idx,epm);
  TFile *fout = new TFile( fout_path.c_str(), "RECREATE" );

  ////////////
  //HISTOGRAMS

  //E-arm
  TH1D *hW2_nocut = new TH1D( "hW2_nocut", "W^{2}, no cut; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_cut = new TH1D( "hW2_cut", "W^{2}, accmatch/global/atime cuts; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );

  //Basic
  TH2D *hdxdy_nocut = new TH2D("hdxdy_nocut","dxdy, no cut; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );
  TH2D *hdxdy_cut = new TH2D("hdxdy_cut","dxdy, all cuts; dy_{HCAL} (m); dx_{HCAL} (m)", 400, -2.0, 2.0, 600, -3.0, 3.0 );

  TH1D *hdx_nocut = new TH1D( "hdx_nocut", "dx, no cut;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);
  TH1D *hdx_cut = new TH1D( "hdx_cut", "dx, e-arm cut, no dy;x_{HCAL}-x_{expect} (m)", hbinfac*harmrange, hcalfit_l, hcalfit_h);

  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;

  // setup reporting indices
  std::cout << std::endl << std::endl;

  for (Int_t irun=0; irun<nruns; irun++) {

    //DEBUG
    //if( irun > 1 ) break;

    // accessing run info
    Int_t runnum = corun[irun].runnum;
    Int_t mag = corun[irun].sbsmag / 21; //convert to percent
    Double_t ebeam = corun[irun].ebeam; //get beam energy per run
    std::string tar = corun[irun].target;

    //Only proceed for deuterium at desired field setting
    if( mag!=magset || tar.compare("LD2")!=0 ) 
      continue;

    std::cout << "Analyzing run " << runnum << ".." << std::endl;

    //set up configuration and tune objects to load analysis parameters
    SBSconfig config(kine,mag);

    //Obtain configuration pars from config file
    Double_t hcaltheta = config.GetHCALtheta_rad();
    Double_t hcaldist = config.GetHCALdist();
    Double_t sbsdist = config.GetSBSdist();
    Double_t bbthr = config.GetBBtheta_rad(); //in radians
    
    //SBStune *tune = new SBStune(kine,mag);
    SBStune tune(kine,mag);
    
    //Reporting. tar should always equal curtar as categorized by good run list
    std::cout << "Settings are.." << endl;
    cout << endl << config << endl;
    cout << endl << tune << endl;

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
    Double_t atimediff0   = tune.Getatimediff0();
    Double_t atimediffsig = tune.Getatimediffsig();

    Double_t W2min = W2mean - 2*W2sig;
    Double_t W2max = W2mean + 2*W2sig;

    std::string rfname = rootfile_dir + Form("/*%d*",corun[irun].runnum);
    //std::cout << "Switching to run " << rfname << ".." << std::endl;
    C = new TChain("T");
    C->Add(rfname.c_str());

    // setting up ROOT tree branch addresses
    C->SetBranchStatus("*",0);    

    // HCal general
    Double_t hcalid, hcale, hcalx, hcaly, hcalr, hcalc, hcaltdc, hcalatime, hcalidx, hcalnclus;
    std::vector<std::string> hcalvar = {"idblk","e","x","y","rowblk","colblk","tdctimeblk","atimeblk","index","nclus"};
    std::vector<void*> hcalvarlink = {&hcalid,&hcale,&hcalx,&hcaly,&hcalr,&hcalc,&hcaltdc,&hcalatime,&hcalidx,&hcalnclus};
    rvars::setbranch(C, "sbs.hcal", hcalvar, hcalvarlink);

    // HCal cluster branches (primary block information for id, tdc, and atime)
    Double_t hcalcid[econst::maxclus], hcalce[econst::maxclus], hcalcx[econst::maxclus], hcalcy[econst::maxclus], hcalctdctime[econst::maxclus], hcalcatime[econst::maxclus];
    Int_t Nhcalcid;
    std::vector<std::string> hcalcvar = {"id","e","x","y","tdctime","atime","id"};
    std::vector<void*> hcalcvarlink = {&hcalcid,&hcalce,&hcalcx,&hcalcy,&hcalctdctime,&hcalcatime,&Nhcalcid};
    rvars::setbranch(C, "sbs.hcal.clus", hcalcvar, hcalcvarlink, 6);
    
    // hodoscope cluster mean time
    Int_t Nhodotmean; 
    Double_t hodotmean[econst::maxclus];
    std::vector<std::string> hodovar = {"clus.tmean","clus.tmean"};
    std::vector<void*> hodovarlink = {&Nhodotmean,&hodotmean};
    rvars::setbranch(C, "bb.hodotdc", hodovar, hodovarlink, 0); 

    // BBCal shower timing
    Double_t atime_sh;
    std::vector<std::string> shvar = {"atimeblk"};
    std::vector<void*> shvarlink = {&atime_sh};
    rvars::setbranch(C, "bb.sh", shvar, shvarlink);  
    
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
    C->SetBranchStatus("bb.sh.nclus", 1);
    C->SetBranchStatus("sbs.hcal.nclus", 1);

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
      
      //DEBUG
      //if( nevent > 100000 ) break;

      ///////
      //Single-loop globalcut method. Save pass/fail for output tree.
      currenttreenum = C->GetTreeNumber();
      if( nevent == 1 || currenttreenum != treenum ){
	treenum = currenttreenum; 
	GlobalCut->UpdateFormulaLeaves();
	//cout << "Updating formula leaves and switching segment at event: " << nevent << endl;
      }
      bool failedglobal = GlobalCut->EvalInstance(0) == 0;
      
      ///////
      //E-arm physics calculations
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

      /////////////////////
      //W2 elastic cut bool
      bool failedW2 = W2<W2min || W2>W2max;

      ///////
      //HCal active area cut (acceptance matching). Save pass/fail for output tree.
      vector<Double_t> xyhcalexp; vars::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
      bool failedaccmatch = 
	xyhcalexp[1] > (econst::hcalposYf_mc-econst::hcalblk_div_h) ||
	xyhcalexp[1] < (econst::hcalposYi_mc+econst::hcalblk_div_h) ||
	xyhcalexp[0] > (econst::hcalposXf_mc-econst::hcalblk_div_v) ||
	xyhcalexp[0] < (econst::hcalposXi_mc+econst::hcalblk_div_v);
      
      //Get expected proton deflection for thetapq (this needs work to be useful)
      Double_t protdeflect = tan( 0.3 * BdL / qv.Mag() ) * (hcaldist - (sbsdist + econst::sbsdipolegap/2.0) );
      
      //Get some primary cluster primary block variables
      Double_t dx = hcalx - xyhcalexp[0];
      Double_t dy = hcaly - xyhcalexp[1];

      //Primary cluster atime cut
      bool failedpclusatime =  abs( hcalcatime[0] - atime0 ) > coinSigFac*atimesig;

      //dy cut
      bool failedpclusdy = abs(dy-dy0)>dysig;

      //E-arm cuts first
      if( failedglobal || failedaccmatch || failedW2 )
	continue;

      //////////////////////
      //ALL CLUSTER ANALYSIS
      
      //Set up clone clusters for selection analysis.
      vector<double> clone_cluster_e;
      vector<double> clone_cluster_intime;
      vector<double> clone_cluster_score;

      //loop through all clusters and select without HCal position information
      for( int c=0; c<Nhcalcid; c++ ){
	
	//calculate h-arm physics quantities per cluster
	double atime = hcalcatime[c];
	double atime_diff = atime - atime_sh; //Assuming best shower time on primary cluster
	double ce = hcalce[c];
	
	//using hcal atime until after pass2, wide cut around 5 sigma
	bool passedCoin = abs(atime_diff-coin_profile[1])<5*coin_profile[2];
	bool passedE = ce>minE;
	
	//Build clone for later sorting
	clone_cluster_e.push_back(ce);

	//Replicate the in-time algorithm with new cluster to be sorted later
	clone_cluster_intime.push_back(ce);
	if( !passedCoin )
	  clone_cluster_intime[c] = 0;
	
	//Get score (no position info). Will be sorted later
	double cascore = assignScore( ce, atime_diff, hcalce[(Int_t)hcalidx], coin_profile);
	clone_cluster_score.push_back(cascore);
	
      }//endloop over cluster elements
      
      // auto [sortedEValues, sortedEIndices] = sortNonZero(clone_cluster_e);
      // auto [sortedIntimeValues, sortedIntimeIndices] = sortNonZero(clone_cluster_intime);
      // auto [sortedAltValues, sortedAltIndices] = sortNonZero(clone_cluster_score);

      // Int_t cidx_e = sortedEIndices[0];
      // Int_t cidx_intime = sortedIntimeIndices[0];
      // Int_t cidx_score = sortedAltIndices[0];

      Int_t cidx_e = (Int_t)hcalidx;
      
      if( hcale != hcalce[(Int_t)hcalidx] )
	cerr << "ERROR: Sorting failure. Debug and rerun." << endl;

      Int_t score_idx = -1;
      Int_t intime_idx = -1;
      //Int_t e_idx = -1;
      Double_t score = 0.;
      Double_t intime = 0.;
      //Double_t energy = 0.;
      for( int c=0; c<Nhcalcid; c++ ){
	if( clone_cluster_score[c]>score ){
	  score = clone_cluster_score[c];
	  score_idx = c;
	}
	if( clone_cluster_intime[c]>intime ){
	  intime = clone_cluster_intime[c];
	  intime_idx = c;
	}
	// if( hcalce[c]>energy ){
	//   energy = hcalce[c];
	//   e_idx = c;
	// }
      }

      Int_t cidx_intime = intime_idx;
      Int_t cidx_score = score_idx;
      //Int_t cidx_e = e_idx;

      //Switch between best clusters for systematic analysis
      Int_t cidx_best;
      
      switch (cluster_idx) {
      case 1:
	cidx_best = 0;
	break;
      case 2:
	cidx_best = cidx_e;
	break;
      case 3:
	cidx_best = cidx_intime;
	break;
      case 4:
	cidx_best = cidx_score;
	break;
      default:
	cidx_best = 3;
      }
      
      //Calculations from the best cluster for later use
      Double_t dx_bestcluster = hcalcx[cidx_best] - xyhcalexp[0];
      Double_t dy_bestcluster = hcalcy[cidx_best] - xyhcalexp[1];
      Double_t hatime_bestcluster = hcalcatime[cidx_best];
      Double_t hatimediff_bestcluster = hcalcatime[cidx_best] - atime_sh;
      Double_t ce_bestcluster = hcalce[cidx_best];
      bool passedcoin = abs( hatime_bestcluster - atimediff0 ) < coinSigFac*atimediffsig;
      bool faileddy = abs(dy_bestcluster-dy0)>dysig;

      //Fill nocut histograms for both
      hW2_nocut->Fill(W2);
      hdxdy_nocut->Fill(dy_bestcluster,dx_bestcluster);
      hdx_nocut->Fill(dx_bestcluster);

      //H-arm cuts last
      if( faileddy )
	continue;

      //Fill cut histograms for both
      hW2_cut->Fill(W2);
      hdxdy_cut->Fill(dy_bestcluster,dx_bestcluster);
      hdx_cut->Fill(dx_bestcluster);

    } //end event loop
    
    // reset chain for the next run config
    C->Reset();
    
  } //end run loop
  
  fout->Write();
  
  cout << endl << "Analysis complete. Outfile located at " << fout_path << endl;

  // Send time efficiency report to console
  st->Stop();
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
