//SSeeds 1.22.23

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

const Int_t maxClus = 35; //larger than max per tree
const Int_t maxBlk = 25;

//Side-band pol4 fit requires beginning and end of rejection region declared here updated later
double SBpol4rej_b = -1.7; //Central-band fit begin
double SBpol4rej_e = 0.7; //Central-band fit end

//Side-band background fit, pol4 with RejectPoint()
double BGfit(double *x, double *par){
  double yint = par[0];
  double p1 = par[1];
  double p2 = par[2];
  double p3 = par[3];
  double p4 = par[4];

  //TF1::RejectPoint() for sideband analysis
  if(x[0]>SBpol4rej_b && x[0]<SBpol4rej_e) { 
    TF1::RejectPoint();
    return 0;
  }

  return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3)+p4*pow(x[0],4);
}

//MAIN
void systematics( Int_t kine=4, Int_t magset=30, Int_t pass=1 )
{

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // Catch first set together for consistency
  if(pass==0)
    pass=1;

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
  Double_t coin_sigma_factor = jmgr->GetValueFromSubKey<Double_t>( "coin_sigma_factor", Form("sbs%d",kine) );
  vector<Double_t> coin_profile;
  jmgr->GetVectorFromSubKey<Double_t>("coin_profile",Form("sbs%d",kine),coin_profile);
  Double_t hcal_v_offset = jmgr->GetValueFromSubKey<Double_t>( "hcal_offset", Form("sbs%d",kine) );

  std::cout << "Loaded HCal vertical offset prior to pass2: " << hcal_v_offset << std::endl;

  //set up configuration and tune objects to load analysis parameters
  SBSconfig config(kine,magset);

  //Obtain configuration pars from config file
  Double_t hcaltheta = config.GetHCALtheta_rad();
  Double_t hcaldist = config.GetHCALdist();
  Double_t sbsdist = config.GetSBSdist();
  Double_t bbthr = config.GetBBtheta_rad(); //in radians

  //Set up hcal active area with bounds that match database on pass
  vector<Double_t> hcalaa;
  if(pass<2)
    hcalaa = cut::hcalaa_data_alt(1,1);
  else
    hcalaa = cut::hcalaa_data(1,1);

  //SBStune *tune = new SBStune(kine,mag);
  SBStune tune(kine,magset);
    
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

  Double_t W2min = W2mean - 3*W2sig;
  Double_t W2max = W2mean + 3*W2sig;

  Double_t hcalfit_l = econst::hcalposXi_mc; //lower fit/bin limit for hcal dx plots (m)
  Double_t hcalfit_h = econst::hcalposXf_mc; //upper fit/bin limit for hcal dx plots (m) 
  Double_t harmrange = econst::hcal_vrange; //Full range of hcal dx plots (m)

  //set up files and paths
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fout_path = outdir_path + Form("/gmn_analysis/systematics_sbs%d_mag_%d.root",kine,magset);
  std::string fin_path = outdir_path + Form("/parse/parse_sbs%d_pass%d_update.root",kine,pass);

  if(!util::checkFile){
    cerr << "ERROR: Input file does not exist." << endl;
    return;
  }

  TFile *fout = new TFile( fout_path.c_str(), "RECREATE" );

  ////////////
  //HISTOGRAMS

  



  //HCal edge boundaries
  Double_t left;
  Double_t right;
  Double_t top;
  Double_t bottom;

  //SBS-4 SBS-7 (pass0)
  if( kine==4 || kine==7 ){
    left = econst::hcalposYi_p0;
    right = econst::hcalposYf_p0;
    top = econst::hcalposXf_p0;
    bottom = econst::hcalposXi_p0;
  }else{
    //pass1 and >pass1
    left = econst::hcalposYi;
    right = econst::hcalposYf;
    top = econst::hcalposXf;
    bottom = econst::hcalposXi;
  }

  //Active area boundaries
  Double_t leftAA;
  Double_t rightAA;
  Double_t topAA;
  Double_t bottomAA;

  //SBS-4 SBS-7 (pass0)
  if( kine==4 || kine==7 ){
    leftAA = (econst::hcalposYi_p0+econst::hcalblk_w_p0);
    rightAA = (econst::hcalposYf_p0-econst::hcalblk_w_p0);
    topAA = (econst::hcalposXf_p0-econst::hcalblk_h_p0);
    bottomAA = (econst::hcalposXi_p0+econst::hcalblk_h_p0);
  }else{
    //pass1 and >pass1
    leftAA = (econst::hcalposYi+econst::hcalblk_w);
    rightAA = (econst::hcalposYf-econst::hcalblk_w);
    topAA = (econst::hcalposXf-econst::hcalblk_h);
    bottomAA = (econst::hcalposXi+econst::hcalblk_h);
  }
  
  //Safety Margin boundaries
  Double_t leftSM;
  Double_t rightSM;
  Double_t topSM;
  Double_t bottomSM;

  //Cut all events that are projected to outermost edge of HCal
  //Add 3sigma proton peak safety margin (x and y) to ensure no expected detections lie one boundary of HCal

  //SBS-4 SBS-7 (pass0)
  if( kine==4 || kine==7 ){
    leftSM = (econst::hcalposYi_p0+econst::hcalblk_w_p0+3*dysig);
    rightSM = (econst::hcalposYf_p0-econst::hcalblk_w_p0-3*dysig);
    topSM = (econst::hcalposXf_p0-econst::hcalblk_h_p0-3*dxsig_p-dx0_p);
    bottomSM = (econst::hcalposXi_p0+econst::hcalblk_h_p0+3*dxsig_p-dx0_p);
  }else{
    //pass1 and >pass1
    leftSM = (econst::hcalposYi+econst::hcalblk_w+3*dysig);
    rightSM = (econst::hcalposYf-econst::hcalblk_w-3*dysig);
    topSM = (econst::hcalposXf-econst::hcalblk_h-3*dxsig_p-dx0_p);
    bottomSM = (econst::hcalposXi+econst::hcalblk_h+3*dxsig_p-dx0_p);
  }


  TChain *C = new TChain("P");
  C->Add(fin_path.c_str());
  
  // set up branch variables
  double hcalx, hcaly, hcale;
  double hcalcx[maxClus], hcalcy[maxClus], hcalce[maxClus], hcalcid[maxClus], hcalctdc[maxClus], hcalcatime[maxClus];
  double hcalnclus, hcalnblk;
  double hcalcbid[maxBlk], hcalcbe[maxBlk], hcalcbx[maxBlk], hcalcby[maxBlk], 

  // setting up ROOT tree branch addresses
  C->SetBranchStatus("*",1);    

  C->SetBranchAddress("sbs.hcal.x", &hcalx);
  C->SetBranchAddress("sbs.hcal.y", &hcaly);
  C->SetBranchAddress("sbs.hcal.e", &hcale);
    

    //set up the global cut formula
    TCut GCut = gcut.c_str();
    TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );

    //get each globalcut
    vector<string> allGlobalCuts = util::parseGlobalCut(gcut);

    // get experimental quantities by run
    std::cout << "Uncorrected average beam energy on " << tar << " for run: " << ebeam << std::endl;
    //set up hcal coordinate system with hcal angle wrt exit beamline
    vector<TVector3> hcalaxes; vars::sethcalaxes( hcaltheta, hcalaxes );
    //TVector3 hcalorigin = hcaldist*hcalaxes[2] + econst::hcalvoff*hcalaxes[0];
    TVector3 hcalorigin = hcaldist*hcalaxes[2] + hcal_v_offset*hcalaxes[0];
    Double_t BdL = econst::sbsmaxfield * econst::sbsdipolegap * (mag/100); //scale crudely by magnetic field
    Double_t Eloss_outgoing = econst::celldiameter/2.0/sin(bbthr) * econst::lh2tarrho * econst::lh2dEdx;

    // set nucleon for LD2
    std::string nucleon = "np";

    // event indices
    long nevent = 0, nevents = C->GetEntries(); 

    //ttree formula markers
    Int_t treenum = 0, currenttreenum = 0;

    while (C->GetEntry(nevent++)) {
      
      std::cout << "Processing run (" << irun << "/" << nruns << ") " << corun[irun].runnum << " event " << nevent << "/" << nevents << "\r";
      std::cout.flush();
      
      ///////
      //Single-loop globalcut method. Save pass/fail for output tree.
      currenttreenum = C->GetTreeNumber();
      if( nevent == 1 || currenttreenum != treenum ){
	treenum = currenttreenum; 
	GlobalCut->UpdateFormulaLeaves();
      }
      bool failedglobal = GlobalCut->EvalInstance(0) == 0;

      if(currenttreenum != treenum){
	cout << "Multiple tree numbers for a single run" << endl;
	return;
      }

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

      //Calculate Mott cross section for this event
      Double_t MCS = (pow(physconst::alpha,2)*pow(cos(etheta/2),2)*precon)/(4*pow(ebeam_c,3)*pow(sin(etheta/2),4));

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
	//nu = pbeam.E() - precon;
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

      //Fill for active area check
      hHcalXY->Fill(hcaly,hcalx);								    

      //////////////////////
      //ALL CLUSTER ANALYSIS
      
      //Set up clone clusters for selection analysis.
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

	//Replicate the in-time algorithm with new cluster to be sorted later
	clone_cluster_intime.push_back(ce);
	if( !passedCoin )
	  clone_cluster_intime[c] = 0;
	
	//Get score (no position info). Will be sorted later
	double cascore = util::assignScore( ce, atime_diff, hcalce[(Int_t)hcalidx], coin_profile);
	clone_cluster_score.push_back(cascore);
	
      }//endloop over cluster elements

      Int_t cidx_e = (Int_t)hcalidx;
      
      if( hcale != hcalce[(Int_t)hcalidx] ){
	cerr << "ERROR: Sorting failure. Debug and rerun." << endl;
	cout << "Index: " << (Int_t)hcalidx << ", hcale:" << hcale << endl;
      }

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
      
      //Calculations from the best cluster
      Double_t dx_bestcluster = hcalcx[cidx_best] - xyhcalexp[0];
      Double_t dy_bestcluster = hcalcy[cidx_best] - xyhcalexp[1];
      Double_t hatime_bestcluster = hcalcatime[cidx_best];
      Double_t hcoin_bestcluster = hcalcatime[cidx_best] - atime_sh;
      Double_t hcoin_pcluster = hcalcatime[0] - atime_sh;
      Double_t ce_bestcluster = hcalce[cidx_best];
      Int_t cnblk = (Int_t)hcalcnblk[cidx_best];
      //bool passedcoin = abs( hatime_bestcluster - atimediff0 ) < coin_sigma_factor*atimediffsig;

      //H-arm fiducial cuts
      bool hcalON = cut::hcalaaON(hcalcx[cidx_best],hcalcy[cidx_best],hcalaa);
      vector<Double_t> fid = cut::hcalfid(dxsig_p,dysig,hcalaa);
      bool passed_fid = cut::hcalfidIN(xyhcalexp[0],xyhcalexp[1],dx0_p,fid);

      //Fill analysis tree variables before making cuts on systematicInfo
      if(systematicInfo){
	dx_out = dx_bestcluster;
	dy_out = dy_bestcluster;
	xexp_out = xyhcalexp[0];
	yexp_out = xyhcalexp[1];
	hcalx_out = hcalcx[cidx_best];
	hcaly_out = hcalcy[cidx_best];
	W2_out = W2;
	Q2_out = Q2;
	mott_out = MCS;
	hcalE_out = ce_bestcluster;
	bbEtot_out = e_sh + e_ps;
	bbshE_out = e_sh;
	bbpsE_out = e_ps;
	hcalatime_out = hatime_bestcluster;
	bbshatime_out = atime_sh;
	bbpsatime_out = atime_ps;
	bbgemNhits_out = gem_hits[0];
	bbtrx_out = xtr[0];
	bbtry_out = ytr[0];
	bbtrp_out = p[0];
	coinP1_out = coin_profile[1];
	coinP2_out = coin_profile[2];
	dy0_out = dy0;
	dx0_p_out = dx0_p;
	dx0_n_out = dx0_n;
	dysig_out = dysig;
	dxsig_p_out = dxsig_p;
	dxsig_n_out = dxsig_n;
	W2mean_out = W2mean;
	W2sig_out = W2sig;

	hcalnclus_out = Nhcalcid;
	hcalnblk_out = cnblk;
	bbshnclus_out = (Int_t)nclus_sh;
	bbshnblk_out = (Int_t)nblk_sh;
	bbpsnclus_out = (Int_t)nclus_ps;
	bbpsnblk_out = (Int_t)nblk_ps;
	bbtrN_out = (Int_t)ntrack;
	passedfid_out = (Int_t)passed_fid;
	run_out = runnum;
	mag_out = mag;
      }      

      if(!failedglobal)
	hdx_cut_globalonly->Fill(dx_bestcluster);

      if(!failedW2)
	hdx_cut_W2only->Fill(dx_bestcluster);

      //E-arm only cuts first
      if( failedglobal || failedW2 )
	continue;

      //Fill mott cross section after elastic cut
      hMott_cs->Fill(MCS);

      //quick check on hcal geometry
      hxy_nocut->Fill(hcalcy[cidx_best],hcalcx[cidx_best]);
      hxy_nocut->Fill(0.,0.);
      if(hcalON)
	hxy_aacut->Fill(hcalcy[cidx_best],hcalcx[cidx_best]);
      if(passed_fid)
	hxy_acccut->Fill(hcalcy[cidx_best],hcalcx[cidx_best]); //this is largely useless

      //Fill e-arm only cut histograms
      hcoin->Fill(hcoin_bestcluster);
      hcoin_pclus->Fill(hcoin_pcluster);
      hW2_nocut->Fill(W2);
      hdxdy_nocut->Fill(dy_bestcluster,dx_bestcluster);
      hdx_nocut->Fill(dx_bestcluster);
      hxy_exp_n->Fill(xyhcalexp[1],xyhcalexp[0]);
      hxy_exp_p->Fill(xyhcalexp[1],xyhcalexp[0]+dx0_p); //expected proton position from average difference between dxdy neutron and proton spots

      if(passed_fid){
	hxy_exp_n_fid->Fill(xyhcalexp[1],xyhcalexp[0]);
	hxy_exp_p_fid->Fill(xyhcalexp[1],xyhcalexp[0]+dx0_p); //expected proton position from average difference between dxdy neutron and proton spots
      }

      //Both arm cuts
      bool failedcoin = abs( hcoin_bestcluster - coin_profile[1] ) > coin_sigma_factor*coin_profile[2];
      bool faileddy = abs( dy_bestcluster - dy0 ) > 3*dysig;

      if( !faileddy && !failedcoin && hcalON ){
	hdx_cut_nofid->Fill(dx_bestcluster); //primary dx histo, no fiducial cut
	hdxdy_cut_nofid->Fill(dy_bestcluster,dx_bestcluster);
	if(!passed_fid)
	  hdx_cut_failfid->Fill(dx_bestcluster); //primary dx histo, 
      }
      
      //Make both arm cuts
      if( !passed_fid || !hcalON || failedcoin || faileddy )
	continue;

      hdxvE->Fill(ce_bestcluster,dx_bestcluster);
      hcoin_cut->Fill(hcoin_bestcluster);
      hcoin_pclus_cut->Fill(hcoin_pcluster);
      hW2_cut->Fill(W2);
      hdxdy_cut->Fill(dy_bestcluster,dx_bestcluster);
      hdx_cut->Fill(dx_bestcluster); //primary dx histo, fiducial cut
      hHcalXY_allcuts->Fill(hcalcy[cidx_best],hcalcx[cidx_best]);

      P->Fill();

    } //end event loop
    
    // reset chain for the next run config
    C->Reset();
    
  } //end run loop
  
  //Make fiducial check on neutron hypothesis
  vector<Double_t> fid = cut::hcalfid(dxsig_p,dysig,hcalaa);
  Double_t hcalx_ft = fid[2];
  Double_t hcalx_fb = fid[3];
  Double_t hcaly_fr = fid[0];
  Double_t hcaly_fl = fid[1];

  //Make Fiducial Safety Margin TLines
  TLine* l1_f = new TLine(hcalx_ft, hcaly_fr, hcalx_ft, hcaly_fl);
  l1_f->SetLineColor(kGreen);
  TLine* l2_f = new TLine(hcalx_fb, hcaly_fr, hcalx_fb, hcaly_fl);
  l2_f->SetLineColor(kGreen);
  TLine* l3_f = new TLine(hcalx_ft, hcaly_fr, hcalx_fb, hcaly_fr);
  l3_f->SetLineColor(kGreen);
  TLine* l4_f = new TLine(hcalx_ft, hcaly_fl, hcalx_fb, hcaly_fl);
  l4_f->SetLineColor(kGreen);

  //Make Fiducial Line
  // TLine* l1_fl = new TLine(hcalx_t-2*dysig, hcaly_fl, hcalx_t-2*dysig, hcaly_fl+dx0_p);
  // l1_fl->SetLineColor(kMagenta);

  TCanvas *c0 = new TCanvas("c0","HCal Fiducial cuts, neutron hypothesis",1200,500);
  c0->Divide(2,1);

  c0->cd(1);
  c0->SetLogz();
  hxy_exp_n->Draw("colz");
  l1_f->Draw();
  l2_f->Draw();
  l3_f->Draw();
  l4_f->Draw();

  auto leg0a = new TLegend(0.1,0.8,0.5,0.9);
  leg0a->AddEntry( l1_f, "HCal Safety Margin", "l");
  leg0a->Draw();

  c0->Update();

  c0->cd(2);
  c0->SetLogz();
  hxy_exp_n_fid->Draw("colz");
  l1_f->Draw();
  l2_f->Draw();
  l3_f->Draw();
  l4_f->Draw();
  //l1_fl->Draw();

  auto leg0b = new TLegend(0.1,0.8,0.5,0.9);
  leg0b->AddEntry( l1_f, "HCal Safety Margin", "l");
  leg0b->Draw();

  c0->Update();
  c0->Write();

  //Make fiducial check on proton hypothesis

  TCanvas *c1 = new TCanvas("c1","HCal Fiducial cuts, proton hypothesis",1200,500);
  gStyle->SetPalette(53);
  c1->Divide(2,1);

  c1->cd(1);
  c1->SetLogz();
  hxy_exp_p->Draw("colz");
  l1_f->Draw();
  l2_f->Draw();
  l3_f->Draw();
  l4_f->Draw();

  auto leg1a = new TLegend(0.1,0.8,0.5,0.9);
  leg1a->AddEntry( l1_f, "HCal Safety Margin", "l");
  leg1a->Draw();

  c1->Update();

  c1->cd(2);
  c1->SetLogz();
  hxy_exp_p_fid->Draw("colz");
  l1_f->Draw();
  l2_f->Draw();
  l3_f->Draw();
  l4_f->Draw();
  //l1_fl->Draw();

  auto leg1b = new TLegend(0.1,0.8,0.5,0.9);
  leg1b->AddEntry( l1_f, "HCal Safety Margin", "l");
  leg1b->Draw();

  c1->Update();
  c1->Write();

  fout->Write();
  
  cout << endl << "Analysis complete. Outfile located at " << fout_path << endl;

  // Send time efficiency report to console
  st->Stop();
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
