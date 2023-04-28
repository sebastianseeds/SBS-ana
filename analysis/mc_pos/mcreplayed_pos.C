//sseeds 04.19.23 - Script to extract hcal position resolution from MC

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

//Uses g4sbs replays of LH2 data
void mcreplayed_pos( Int_t kine=4, Int_t epm=1, bool wide=false )
{ //main  
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  std::string date = util::getDate();

  //general static parameters for this analysis
  const Double_t detxmin = -3.0; //nucleon, GeV
  const Double_t detxmax = 3.0; //nucleon, GeV
  const Double_t detymin = -3.0; //nucleon, GeV
  const Double_t detymax = 3.0; //nucleon, GeV
  Double_t dxmin, dxmax, dymin, dymax;
  if(wide){
    dxmin = -15;
    dxmax = 15;
    dymin = -15;
    dymax = 15;
  }else{
    dxmin = -4; //hcal, GeV
    dxmax = 4; //hcal, GeV
    dymin = -1.5; //hcal, GeV
    dymax = 1.5; //hcal, GeV
  }

  const Int_t nfac = 100; //Number of bins for hcal E vs nucleon p fits obtained to ensure 1000 events per bin

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/smcpos.json");
  
  std::string rdir = jmgr->GetValueFromKey_str( "rootdir" );

  //return all sbs magnetic field settings for the kinematic given
  vector<Int_t> magset;
  jmgr->GetVectorFromSubKey<Int_t>( "magset", Form("sbs%d",kine), magset );

  Int_t Nmagset = magset.size();

  //set up output files
  TFile *fout = new TFile( Form("outfiles/mcrep_pos_sbs%d_epm%d.root",kine,epm), "RECREATE" );
  
  //set up diagnostic histograms

  TH1D *hnu[Nmagset];

  TH1D *hHCALx[Nmagset];
  TH1D *hHCALy[Nmagset];

  TH1D *hMCx[Nmagset];
  TH1D *hMCy[Nmagset];

  TH1D *hdx_p[Nmagset];
  TH1D *hdy_p[Nmagset];

  TH1D *hdx_p_v2[Nmagset];
  TH1D *hdy_p_v2[Nmagset];

  TH2D *hdxdy_p[Nmagset];
  TH2D *hdxdy_p_v2[Nmagset];

  TH1D *hxexp[Nmagset];
  TH1D *hyexp[Nmagset];

  for( Int_t m=0; m<Nmagset; m++ ){
    hnu[m] = new TH1D(Form("hnu_mag%d",magset[m]),Form("HCal elastic KE (nu) mag%d;nu (GeV)",magset[m]), 300,0,10);

    hHCALx[m] = new TH1D(Form("hHCALx_mag%d",magset[m]),Form("HCALx mag%d;x_{HCAL} (m)",magset[m]), nfac*(detxmax-detxmin), detxmin, detxmax);
    hHCALy[m] = new TH1D(Form("hHCALy_mag%d",magset[m]),Form("HCALy mag%d;y_{HCAL} (m)",magset[m]), nfac*(detymax-detymin), detymin, detymax);

    hMCx[m] = new TH1D(Form("hMCx_mag%d",magset[m]),Form("MCx mag%d;x_{MC} (m)",magset[m]), nfac*(detymax-detymin), detymin, detymax);
    hMCy[m] = new TH1D(Form("hMCy_mag%d",magset[m]),Form("MCy mag%d;y_{MC} (m)",magset[m]), nfac*(detxmax-detxmin), detxmin, detxmax);

    hdx_p[m] = new TH1D(Form("hdx_p_mag%d",magset[m]),Form("dx proton (sd track) mag%d;x_{HCAL}-x_{expect} (m)",magset[m]), nfac*(dxmax-dxmin), dxmin, dxmax);
    hdy_p[m] = new TH1D(Form("hdy_p_mag%d",magset[m]),Form("dy proton (sd track) mag%d;y_{HCAL}-y_{expect} (m)",magset[m]), nfac*(dymax-dymin), dymin, dymax);

    hdx_p_v2[m] = new TH1D(Form("hdx_p_v2_mag%d",magset[m]),Form("dx proton (angles) mag%d;x_{HCAL}-x_{expect} (m)",magset[m]), nfac*(dxmax-dxmin), dxmin, dxmax);
    hdy_p_v2[m] = new TH1D(Form("hdy_p_v2_mag%d",magset[m]),Form("dy proton (angles) mag%d;y_{HCAL}-y_{expect} (m)",magset[m]), nfac*(dymax-dymin), dymin, dymax);

    hdxdy_p[m] = new TH2D(Form("hdxdy_p_mag%d",magset[m]),Form("dxdy proton (sd track) mag%d;y_{HCAL}-y_{expect} (m);x_{HCAL}-x_{expect} (m)",magset[m]), nfac*(dymax-dymin), dymin+6, dymax+6, nfac*(dxmax-dxmin), dxmin, dxmax);
    hdxdy_p_v2[m] = new TH2D(Form("hdxdy_p_v2_mag%d",magset[m]),Form("dxdy proton (angles) mag%d;y_{HCAL}-y_{expect} (m);x_{HCAL}-x_{expect} (m)",magset[m]), nfac*(dymax-dymin), dymin, dymax, nfac*(dxmax-dxmin), dxmin, dxmax);

    hxexp[m] = new TH1D(Form("hxexp_mag%d",magset[m]),Form("x exp (angles) mag%d;x_{expect} (m)",magset[m]), nfac*(detxmax-detxmin), detxmin, detxmax);
    hyexp[m] = new TH1D(Form("hyexp_mag%d",magset[m]),Form("y exp (angles) mag%d;y_{expect} (m)",magset[m]), nfac*(detymax-detymin), detymin, detymax);

  }
  
  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;

  // create output tree
  TTree *P = new TTree("P","Analysis Tree"); 

  // output tree vars
  Double_t dx_out;
  Double_t dy_out;
  Double_t dxsd_out;
  Double_t dysd_out;
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
  P->Branch( "dxsd", &dxsd_out, "dxsd/D" );
  P->Branch( "dysd", &dysd_out, "dysd/D" );
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

  //loop over magnetic field settings
  for (Int_t n=0; n<Nmagset; n++) {
        
    string rfile = rdir + Form("replayed_gmn_sbs%d_lh2_%dp_job*",kine,magset[n]);
  
    C = new TChain("T");
    
    C->Add(rfile.c_str());
    
    //setting up ROOT tree branch addresses
    C->SetBranchStatus("*",0);    
    
    // HCal general
    Double_t hcalid, hcale, hcalx, hcaly, hcalr, hcalc, hcaltdc, hcalatime;
    std::vector<std::string> hcalvar = {"idblk","e","x","y","rowblk","colblk","tdctimeblk","atimeblk"};
    std::vector<void*> hcalvarlink = {&hcalid,&hcale,&hcalx,&hcaly,&hcalr,&hcalc,&hcaltdc,&hcalatime};
    rvars::setbranch(C, "sbs.hcal", hcalvar, hcalvarlink);
    
    // MC truth e'
    Double_t mc_p, mc_px, mc_py, mc_pz, mc_vx, mc_vy, mc_vz, mc_nucl;
    std::vector<std::string> mcvar = {"mc_ep","mc_epx","mc_epy","mc_epz","mc_vx","mc_vy","mc_vz","mc_nucl"};
    std::vector<void*> mcvarlink = {&mc_p,&mc_px,&mc_py,&mc_pz,&mc_vx,&mc_vy,&mc_vz,&mc_nucl};
    rvars::setbranch(C, "MC", mcvar, mcvarlink);

    // MC sdtrack
    Double_t mc_posx[econst::maxtrack], mc_posy[econst::maxtrack], mc_tid[econst::maxtrack], mc_pid[econst::maxtrack];
    Int_t nsdtrack;
    std::vector<std::string> mcsdvar = {"sdtrack_posx","sdtrack_posx","sdtrack_posy","sdtrack_tid","sdtrack_pid"};
    std::vector<void*> mcsdvarlink = {&nsdtrack,&mc_posx,&mc_posy,&mc_tid,&mc_pid};
    rvars::setbranch(C, "MC", mcsdvar, mcsdvarlink,0);
    
    // bbcal clus var
    Double_t eSH, xSH, ySH, rblkSH, cblkSH, idblkSH, atimeSH, ePS, rblkPS, cblkPS, idblkPS, atimePS;
    std::vector<std::string> bbcalclvar = {"sh.e","sh.x","sh.y","sh.rowblk","sh.colblk","sh.idblk","sh.atimeblk","ps.e","ps.rowblk","ps.colblk","ps.idblk","ps.atimeblk"};
    std::vector<void*> bbcalclvarlink = {&eSH,&xSH,&ySH,&rblkSH,&cblkSH,&idblkSH,&atimeSH,&ePS,&rblkPS,&cblkPS,&idblkPS,&atimePS};
    rvars::setbranch(C, "bb", bbcalclvar, bbcalclvarlink);

    // track branches
    Double_t ntrack, p[econst::maxtrack],px[econst::maxtrack],py[econst::maxtrack],pz[econst::maxtrack],xtr[econst::maxtrack],ytr[econst::maxtrack],thtr[econst::maxtrack],phtr[econst::maxtrack];
    Double_t vx[econst::maxtrack],vy[econst::maxtrack],vz[econst::maxtrack];
    Double_t xtgt[econst::maxtrack],ytgt[econst::maxtrack],thtgt[econst::maxtrack],phtgt[econst::maxtrack];
    std::vector<std::string> trvar = {"n","p","px","py","pz","x","y","th","ph","vx","vy","vz","tg_x","tg_y","tg_th","tg_ph"};
    std::vector<void*> trvarlink = {&ntrack,&p,&px,&py,&pz,&xtr,&ytr,&thtr,&phtr,&vx,&vy,&vz,&xtgt,&ytgt,&thtgt,&phtgt};
    rvars::setbranch(C,"bb.tr",trvar,trvarlink);

    // set up static params
    SBSconfig config(4,magset[n]); //Simulation files set up for sbs4 at zero field
    Double_t hcaltheta = config.GetHCALtheta_rad();
    Double_t hcaldist = config.GetHCALdist();
    Double_t ebeam = config.GetEbeam();
    Double_t sbsdist = config.GetSBSdist();

    std::string nucleon = "p";

    vector<TVector3> hcalaxes; vars::sethcalaxes( hcaltheta, hcalaxes );
    //TVector3 hcalorigin = hcaldist*hcalaxes[2] + econst::hcalvoff*hcalaxes[0];
    TVector3 hcalorigin = hcaldist*hcalaxes[2];
    Double_t mag = magset[n];
    Double_t BdL = econst::sbsmaxfield * econst::sbsdipolegap * (mag/100); //scale crudely by magnetic field

    // event indices
    long nevent = 0, nevents = C->GetEntries(); 

    while (C->GetEntry(nevent++)) {
      
      if( nevent>nevents-10 ) break;
      std::cout << "Processing sbs" << kine << " field " << magset[n] << "% MC data, event " << nevent << "/" << nevents << "\r";
      std::cout.flush();

      Double_t mcposy=-1000;
      Double_t mcposx=-1000;
      
      for( Int_t i=0; i<nsdtrack; i++ ){
      	//if( mc_tid[i] == 1 && mc_pid[i]==physconst::IDXp){
      	if( mc_pid[i]==physconst::IDXp){
      	  mcposy = mc_posy[i];
      	  mcposx = mc_posx[i];
      	  break;
      	}
      }

      // mcposy = mc_posy[0];
      // mcposx = mc_posx[0];

      //calculated from sd tracks
      Double_t dx = -hcalx - mcposy; //MC convention in hall coordinates, hcal in transport coordinates
      Double_t dy = -hcaly - mcposx;

      hHCALx[n]->Fill(hcalx);
      hHCALy[n]->Fill(hcaly);
      hMCx[n]->Fill(mcposx);
      hMCy[n]->Fill(mcposy);

      //calculated using MC e' momenta using epm1 (MC ep vector is flawless)
      //TVector3 vertex( mc_vx, mc_vy, mc_vz ); //With MC truth
      TVector3 vertex( 0., 0., vz[0] ); //With e-arm tracks
      TLorentzVector pbeam( 0., 0., ebeam, ebeam ); //beam momentum
      //TLorentzVector pe( mc_px, mc_py, mc_pz, mc_p ); //e' plvect, MC truth
      TLorentzVector pe( px[0], py[0], pz[0], p[0] ); //e' plvect e-arm tracks
      TLorentzVector ptarg; vars::setPN(nucleon,ptarg); //target momentum
      TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
      TVector3 qv = q.Vect();
      // TLorentzVector pN = q + ptarg;
      // TVector3 pNhat = pN.Vect().Unit();
      // Double_t Q2 = -q.M2();
      // Double_t W2 = pN.M2();
      // Double_t nu = q.E();
      TLorentzVector pN; //N' momentum
      TVector3 pNhat; //Unit N' 3-vector
      Double_t Q2, W2, nu, thNexp, pNexp;;

      //simple calculations for e' and N'
      Double_t etheta = vars::etheta(pe); 
      Double_t ephi = vars::ephi(pe);
      Double_t pcent = vars::pcentral(ebeam,etheta,nucleon); //e' p reconstructed by angles
      Double_t phNexp = ephi + physconst::pi;
      std::string nucleon = "p"; //For LH2 MC

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
	Q2 = -q.M2();
	W2 = pN.M2();
	nu = q.E();
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
	Q2 = -q.M2();
	W2 = pN.M2();
	nu = q.E();
	pNexp = vars::pN_expect( nu, nucleon );
	thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	pNhat = vars::pNhat_track( thNexp, phNexp );
	pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	std::cout << "Warning: epm version incorrect. Defaulting to version 2." << endl;
      }

      vector<Double_t> xyhcalexp; vars::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
      TVector3 hcalpos = hcalorigin + hcalx*hcalaxes[0] + hcaly*hcalaxes[1];
      // Double_t dx = hcalx - xyhcalexp[0];
      // Double_t dy = hcaly - xyhcalexp[1];
      TVector3 neutdir = (hcalpos - vertex).Unit();
      Double_t protdeflect = tan( 0.3 * BdL / q.Mag() ) * (hcaldist - (sbsdist + econst::sbsdipolegap/2.0) );
      TVector3 protdir = ( hcalpos + protdeflect * hcalaxes[0] - vertex ).Unit();
      Double_t thetapq_p = acos( protdir.Dot( pNhat ) );
      Double_t thetapq_n = acos( neutdir.Dot( pNhat ) );
      // Double_t thNexp = acos( mc_pz / mc_p );
      // Double_t phNexp = atan2( mc_py, mc_px );
      // TVector3 pNhat = vars::pNhat_track( thNexp, phNexp );
      Double_t eoverp = (ePS + eSH) / p[0];
      hnu[n]->Fill(nu);

      //vector<Double_t> xyhcalexp; vars::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
      hxexp[n]->Fill(xyhcalexp[0]);
      hyexp[n]->Fill(xyhcalexp[1]);

      Double_t dx_v2 = hcalx - xyhcalexp[0];
      Double_t dy_v2 = hcaly - xyhcalexp[1];

      hdx_p[n]->Fill(dx);
      hdy_p[n]->Fill(dy);
      hdx_p_v2[n]->Fill(dx_v2);
      hdy_p_v2[n]->Fill(dy_v2);

      hdxdy_p[n]->Fill(dy,dx);
      hdxdy_p_v2[n]->Fill(dy_v2,dx_v2);

      //Fill physics output tree     
      dx_out = dx;
      dy_out = dy;
      dxsd_out = dx_v2;
      dysd_out = dy_v2;
      W2_out = W2;
      Q2_out = Q2;
      nu_out = nu;
      hcale_out = hcale;
      pse_out = ePS;
      she_out = eSH;
      ep_out = p[0];
      eoverp_out = eoverp;
      hcalatime_out = 0;
      hodotmean_out = 0;
      thetapq_pout = thetapq_p;
      thetapq_nout = thetapq_n;
      mag_out = mag;
      run_out = 0;
      tar_out = 0;
      failedglobal_out = 0;
      failedaccmatch_out = 0;
      failedcoin_out = 0;

      P->Fill();

    } //end event loop

    // reset chain for the next run config
    C->Reset();    

  } //end magset loop
  
  //Make canvas for (expected-residuals)/expected
  TCanvas *c[Nmagset];

  for( Int_t m=0; m<Nmagset; m++ ){

    c[m] = new TCanvas(Form("c_m%d",magset[m]),Form("dxdy sbs%d mag%d",kine,magset[m]),1200,500);
    gStyle->SetPalette(53);
    c[m]->Divide(2,1);

    c[m]->cd(1);

    hdxdy_p[m]->Draw("colz");
    
    c[m]->cd(2);
    
    hdxdy_p_v2[m]->Draw("colz");

    c[m]->Write();

  }

  fout->Write();

  // Send time efficiency report to console
  st->Stop();
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
