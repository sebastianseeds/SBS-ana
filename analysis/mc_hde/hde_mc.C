//sseeds 04.18.23 - Script to extract hcal detection efficiency from Juan Carlos Cornejo method detailed by Provakar datta well here: https://sbs.jlab.org/DocDB/0003/000339/001/wSoft_gmnAna_11_18_22.pdf.
//04.25.23 Update - Adds position efficiency information and efficiency ratio
//10.23.23 Update - Fixing memory issues and adding error/fits

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

const double pmin = 1.0; //nucleon, GeV
const double pmax = 9.0; //nucleon, GeV
const double Emin = 0.; //hcal, GeV
const double Emax = 2.; //hcal, GeV
const int nbin = 150; //Number of bins for hcal E vs nucleon p fits obtained to ensure 1000 events per bin
const double fit_fac = 0.5;

// Function to fit slices of a TH2D histogram with Gaussian
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> FitSlices(TH2D* hist2D, int nbins) {
    std::vector<double> x_values;
    std::vector<double> fit_means;
    std::vector<double> fit_sigmas;
    
    //double fit_range = fit_fac*;  // Adjust this value for the desired fit range
    double fit_range=1.;

    // Loop over the specified number of bins
    for (int bin = 1; bin <= nbins; ++bin) {
        // Calculate the center of the current bin
        double x_center = hist2D->GetXaxis()->GetBinCenter(bin);
        x_values.push_back(x_center);

	fit_range = fit_fac*x_center;

        // Get the projection of the current bin along the Y-axis
        TH1D* slice = hist2D->ProjectionY("", bin, bin);

        // Find the maximum bin content and its bin number
        int max_bin = slice->GetMaximumBin();
        double max_content = slice->GetBinContent(max_bin);

        // Define the fit range symmetrically around the maximum bin
        int min_fit_bin = std::max(1, max_bin - int(fit_range * slice->GetNbinsX()));
        int max_fit_bin = std::min(slice->GetNbinsX(), max_bin + int(fit_range * slice->GetNbinsX()));

        double fit_min = slice->GetXaxis()->GetBinCenter(min_fit_bin);
        double fit_max = slice->GetXaxis()->GetBinCenter(max_fit_bin);

        // Create and fit a Gaussian to the slice within the defined range
        TF1* fit_function = new TF1("gaussian", "gaus", fit_min, fit_max);
        fit_function->SetParameters(max_content, x_center, 0.1);  // Initial parameters

        TCanvas* c = new TCanvas();
        slice->Fit(fit_function, "RQ");

        // Retrieve the fit results (mean and sigma)
        double mean = fit_function->GetParameter(1);
        double sigma = fit_function->GetParameter(2);

        fit_means.push_back(mean);
        fit_sigmas.push_back(sigma);

        // Clean up
        delete fit_function;
        delete slice;
        delete c;
    }

    return std::make_tuple(x_values, fit_means, fit_sigmas);
}

//Uses g4sbs replays of simulated data set containing pgun/ngun, zero field, SBS4 geometry
void hde_mc( Int_t iter = 1, Double_t tfac = 4. ) //iteration 0 gets mean values of hcalE vs nucleonp; 1 gets eff
{ //main  
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  std::string date = util::getDate();

  //general static parameters for this analysis
  Double_t p_step = (pmax-pmin)/nbin; //Amount of momentum traversed per bin
  Double_t oEmean_p[nbin] = {0.}; //hcal proton E mean values obtained from iter 0
  Double_t oEmean_n[nbin] = {0.}; //hcal neutron E mean values obtained from iter 0

  // Create files to store hcal E fit parameters for iter 1 cuts
  ofstream E_prot;
  ofstream E_neut;

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/smchde.json");
  
  std::string proton_rfile = jmgr->GetValueFromSubKey_str( "rootfiles", "proton" );
  std::string neutron_rfile = jmgr->GetValueFromSubKey_str( "rootfiles", "neutron" );
  std::string prot_param_path = jmgr->GetValueFromSubKey_str( "params", "proton" );
  std::string neut_param_path = jmgr->GetValueFromSubKey_str( "params", "neutron" );
  Double_t pp0 = jmgr->GetValueFromSubKey<Double_t>( "Efit_p0", "proton" );
  Double_t pp1 = jmgr->GetValueFromSubKey<Double_t>( "Efit_p1", "proton" );
  Double_t np0 = jmgr->GetValueFromSubKey<Double_t>( "Efit_p0", "neutron" );
  Double_t np1 = jmgr->GetValueFromSubKey<Double_t>( "Efit_p1", "neutron" );

  Int_t nruns = -1; //Always analyze all available runs
  
  if( iter==1 ){
    //Get proton E means
    ifstream prot_param_file( prot_param_path );
    if( !prot_param_file ){
      cerr << endl << "ERROR: No input constant file present -> path to meanE_proton.txt expected." << endl;
      return 0;
    }
  
    Int_t n1=0;
    Double_t d1;
    string line1;
    
    while( getline( prot_param_file, line1 ) ){
      if( line1.at(0) == '#' ) {
	continue;
      }
      stringstream ss( line1 );
      ss >> d1;     
      oEmean_p[n1] = d1;      
      n1++;
    }

    //Get neutron E means
    ifstream neut_param_file( neut_param_path );
    if( !neut_param_file ){
      cerr << endl << "ERROR: No input constant file present -> path to meanE_proton.txt expected." << endl;
      return 0;
    }
  
    n1=0;
    d1=0;
    string line2;
    
    while( getline( neut_param_file, line2 ) ){
      if( line2.at(0) == '#' ) {
	continue;
      }
      stringstream ss( line2 );
      ss >> d1;     
      oEmean_n[n1] = d1;      
      n1++;
    }

    //console out to check
    cout << endl << "proton means: " << endl;
    for( int i=0; i<nbin; ++i )
      cout << oEmean_p[i] << endl;
    
    cout << endl << endl << "neutron means: " << endl;
    for( int i=0; i<nbin; ++i )
      cout << oEmean_n[i] << endl;
  }

  //set up output files
  TFile *fout = new TFile( Form("outfiles/hdemc_i%d_tfac%0.0f.root",iter,tfac), "RECREATE" );
  
  //set up diagnostic histograms
  TH2D *hEdepvP_p = new TH2D("hEdepvP_p","HCal E dep vs proton momentum; p_{proton} (GeV); E_{hcal} (GeV)", nbin, 1, 9, 200, Emin, Emax);
  TH2D *hEdepvP_n = new TH2D("hEdepvP_n","HCal E dep vs neutron momentum; p_{neutron} (GeV); E_{hcal} (GeV)", nbin, 1, 9, 200, Emin, Emax);
  TH2D *hEdepvP_p_Ecut = new TH2D("hEdepvP_p_Ecut","HCal E dep vs proton momentum, Mean E / 4 cut; p_{proton} (GeV); E_{hcal} (GeV)", nbin, 1, 9, 200, Emin, Emax);
  TH2D *hEdepvP_n_Ecut = new TH2D("hEdepvP_n_Ecut","HCal E dep vs neutron momentum, Mean E / 4 cut; p_{neutron} (GeV); E_{hcal} (GeV)", nbin, 1, 9, 200, Emin, Emax);
  TH1D *hdx_p = new TH1D("hdx_p","dx proton (sd track);x_{HCAL}-x_{expect} (m)", 400, -2, 2);
  TH1D *hdy_p = new TH1D("hdy_p","dy proton (sd track);y_{HCAL}-y_{expect} (m)", 400, 3.8, 7.8);
  TH1D *hdx_n = new TH1D("hdx_n","dx neutron (sd track);x_{HCAL}-x_{expect} (m)", 400, -2, 2);
  TH1D *hdy_n = new TH1D("hdy_n","dy neutron (sd track);y_{HCAL}-y_{expect} (m)", 400, 3.8, 7.8);
  TH1D *hdx_p_v2 = new TH1D("hdx_p_v2","dx proton (angles);x_{HCAL}-x_{expect} (m)", 400, -2, 2);
  TH1D *hdy_p_v2 = new TH1D("hdy_p_v2","dy proton (angles);y_{HCAL}-y_{expect} (m)", 400, -2, 2);
  TH1D *hdx_n_v2 = new TH1D("hdx_n_v2","dx neutron (angles);x_{HCAL}-x_{expect} (m)", 400, -2, 2);
  TH1D *hdy_n_v2 = new TH1D("hdy_n_v2","dy neutron (angles);y_{HCAL}-y_{expect} (m)", 400, -2, 2);
  TH1D *hxexp = new TH1D("hxexp","x exp (angles);x_{expect} (m)", 400, -2, 2);
  TH1D *hyexp = new TH1D("hyexp","y exp (angles);y_{expect} (m)", 400, -2, 2);
  
  TH2D *hdxvp_p = new TH2D("hdxvp_p","dx vs proton p; p_{p} (GeV); x_{HCAL}-x_{expect} (m)", nbin, 1, 9, 400, -2, 2);
  TH2D *hdyvp_p = new TH2D("hdyvp_p","dy vs proton p; p_{p} (GeV); y_{HCAL}-y_{expect} (m)", nbin, 1, 9, 400, -2, 2);
  TH2D *hdxvp_n = new TH2D("hdxvp_n","dx vs neutron p; p_{n} (GeV); x_{HCAL}-x_{expect} (m)", nbin, 1, 9, 400, -2, 2);
  TH2D *hdyvp_n = new TH2D("hdyvp_n","dy vs neutron p; p_{n} (GeV); y_{HCAL}-y_{expect} (m)", nbin, 1, 9, 400, -2, 2);

  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;
  std::string nuc;
  
  // declare ints to hold number of nucleons that pass/fail energy check on iter 1
  Int_t P_pass[nbin] = {0};
  Int_t P_tot[nbin] = {0};
  Int_t N_pass[nbin] = {0};
  Int_t N_tot[nbin] = {0};

  //loop over nucleons
  for (Int_t n=0; n<2; n++) {
    
    C = new TChain("T");
    
    if( n==0 ){ 
      nuc = "proton";
      C->Add(proton_rfile.c_str());
    }else if( n==1 ){ 
      nuc = "neutron";
      C->Add(neutron_rfile.c_str());
    }else
      break;
    
    if( C->GetEntries()==0 ){
      cout << "ERROR: Missing MC file or empty MC file." << endl;
      return;
    } 

    // setting up ROOT tree branch addresses
    C->SetBranchStatus("*",0);    
    
    // HCal general
    Double_t hcalid, hcale, hcalx, hcaly, hcalr, hcalc, hcaltdc, hcalatime;
    std::vector<std::string> hcalvar = {"idblk","e","x","y","rowblk","colblk","tdctimeblk","atimeblk"};
    std::vector<void*> hcalvarlink = {&hcalid,&hcale,&hcalx,&hcaly,&hcalr,&hcalc,&hcaltdc,&hcalatime};
    rvars::setbranch(C, "sbs.hcal", hcalvar, hcalvarlink);
    
    // MC nucleon
    Double_t mc_p, mc_px, mc_py, mc_pz, mc_vx, mc_vy, mc_vz, mc_nucl, mc_posx, mc_posy;
    std::vector<std::string> mcvar = {"mc_ep","mc_epx","mc_epy","mc_epz","mc_vx","mc_vy","mc_vz","mc_nucl","sdtrack_posx","sdtrack_posy"};
    std::vector<void*> mcvarlink = {&mc_p,&mc_px,&mc_py,&mc_pz,&mc_vx,&mc_vy,&mc_vz,&mc_nucl,&mc_posx,&mc_posy};
    rvars::setbranch(C, "MC", mcvar, mcvarlink);
    
    // set up hcal
    SBSconfig config(4,0); //Simulation files set up for sbs4 at zero field
    Double_t hcaltheta = config.GetHCALtheta_rad();
    Double_t hcaldist = config.GetHCALdist();

    vector<TVector3> hcalaxes; vars::sethcalaxes( hcaltheta, hcalaxes );
    //TVector3 hcalorigin = hcaldist*hcalaxes[2] + econst::hcalvoff*hcalaxes[0];
    TVector3 hcalorigin = hcaldist*hcalaxes[2];

    // event indices
    long nevent = 0, nevents = C->GetEntries(); 

    while (C->GetEntry(nevent++)) {
      
      std::cout << "Processing " << nuc << " MC data, event " << nevent << "/" << nevents << "\r";
      std::cout.flush();
      
      Int_t bin = (mc_p-pmin)/p_step; //will fail if MC isn't configured to match pmin and pmax
      Double_t E_thresh;

      if( bin>nbin || bin < 0 ){
	std::cout << "Warning: nucleon momentum bin out of bounds at " << bin << endl;
	continue;
      }

      //calculated from sd tracks
      Double_t dx = -hcalx - mc_posy; //MC convention in hall coordinates, hcal in transport coordinates
      Double_t dy = hcaly - mc_posx;

      //calculated using MC nucleon momenta
      TVector3 vertex( mc_vx, mc_vy, mc_vz );
      Double_t thNexp = acos( mc_pz / mc_p );
      Double_t phNexp = atan2( mc_py, mc_px );
      TVector3 pNhat = vars::pNhat_track( thNexp, phNexp );

      vector<Double_t> xyhcalexp; vars::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
      hxexp->Fill(xyhcalexp[0]);
      hyexp->Fill(xyhcalexp[1]);

      Double_t dx_v2 = hcalx - xyhcalexp[0];
      Double_t dy_v2 = hcaly - xyhcalexp[1];

      if( n==0 ){
	P_tot[bin]++;
	E_thresh = oEmean_p[bin]/tfac;
	if( hcale > E_thresh ){
	  P_pass[bin]++;
	  hEdepvP_p_Ecut->Fill(mc_p,hcale);
	}
	hEdepvP_p->Fill(mc_p,hcale);
	hdx_p->Fill(dx);
	hdy_p->Fill(dy);
	hdxvp_p->Fill(mc_p,dx_v2);
	hdyvp_p->Fill(mc_p,dy_v2);
	hdx_p_v2->Fill(dx_v2);
	hdy_p_v2->Fill(dy_v2);
      }
      if( n==1 ){
	N_tot[bin]++;
	E_thresh = oEmean_n[bin]/tfac;
	if( hcale > E_thresh ){
	  N_pass[bin]++;
	  hEdepvP_n_Ecut->Fill(mc_p,hcale);
	}
	hEdepvP_n->Fill(mc_p,hcale);
	hdx_n->Fill(dx);
	hdy_n->Fill(dy);
	hdxvp_n->Fill(mc_p,dx_v2);
	hdyvp_n->Fill(mc_p,dy_v2);
	hdx_n_v2->Fill(dx_v2);
	hdy_n_v2->Fill(dy_v2);
      }

    } //end event loop

    // reset chain for the next run config
    C->Reset();

  } //end run loop

  //Reporting
  gStyle->SetOptFit();
  gStyle->SetEndErrorSize(0);
  TCanvas *c1 = new TCanvas("c1","HCal E vs Nucleon P",1600,1200);
  c1->SetGrid();
  c1->cd();

  //Iter 0 arrays
  Double_t binp[nbin] = {0.};

  Double_t Emean_p[nbin] = {0.};
  Double_t Esig_p[nbin] = {0.};
  Double_t Emean_n[nbin] = {0.};
  Double_t Esig_n[nbin] = {0.};
  Double_t binerr_p[nbin] = {0.};
  Double_t binerr_n[nbin] = {0.};

  //Iter 1 arrays
  Double_t hde_proton[nbin] = {0.};
  Double_t hde_neutron[nbin] = {0.};
  Double_t hde_npratio[nbin] = {0.};

  auto mg = new TMultiGraph();

  for(Int_t b=0; b<nbin; b++){
    Double_t p = b*p_step+pmin;
    binp[b] = p;
      
    hde_proton[b] = (double)P_pass[b] / (double)P_tot[b] *100.;
    hde_neutron[b] = (double)N_pass[b] / (double)N_tot[b] *100.;
      
    hde_npratio[b] = hde_neutron[b]/hde_proton[b];

  }
    
  //Draw graphs
  //auto grp = new TGraphErrors(nbin,binp,hde_proton,binerr_p,binerr_p);
  auto grp = new TGraph(nbin,binp,hde_proton);
  grp->SetTitle("Proton");
  grp->SetMarkerColor(kRed);
  grp->SetMarkerStyle(20);
  grp->SetMarkerSize(1);
  grp->SetLineColor(kRed);
  grp->SetLineWidth(0);
  mg->Add(grp);

  //auto grn = new TGraphErrors(nbin,binp,hde_neutron,binerr_n,binerr_n);
  auto grn = new TGraph(nbin,binp,hde_neutron);
  grn->SetTitle("Neutron");
  grn->SetMarkerColor(kBlue);
  grn->SetMarkerStyle(21);
  grn->SetMarkerSize(1);
  grn->SetLineColor(kBlue);
  grn->SetLineWidth(0);
  mg->Add(grn);

  mg->SetTitle(Form("HCAL Efficiency (E_{T}=1/%0.0f E_{Peak}) (4x4 cluster)",tfac));
  mg->GetXaxis()->SetTitle("Nucleon Momentum (GeV/c)");
  mg->GetYaxis()->SetTitle("Efficiency (%)");
  mg->Draw("AP");

  mg->GetYaxis()->SetRangeUser(80.,105.);

  c1->Modified();

  c1->BuildLegend();

  c1->Write();

  
  cout << "1" << endl;

  TCanvas *c0 = new TCanvas("c0","HCal Efficiency Ratio (N/P)(E_{T}=1/4 E_{peak})",1600,1200);
  //c2->Divide(2,1);
  c0->SetGrid();
  
  c0->cd();

  auto grr = new TGraph(nbin,binp,hde_npratio);
  grr->SetTitle("HCal Efficiency Ratio (N/P)(E_{T}=1/4 E_{peak})");
  grr->SetMarkerColor(kMagenta);
  grr->SetMarkerStyle(20);
  grr->SetMarkerSize(1);
  grr->SetLineColor(kMagenta);
  grr->SetLineWidth(0);
  grr->GetXaxis()->SetTitle("Nucleon Momentum (GeV/c)");
  grr->GetYaxis()->SetTitle("Efficiency Ratio (N/P)");
  //mg->Add(grr);
  grr->Draw("AP");

  grr->GetYaxis()->SetRangeUser(0.9,1.05);

  c0->Modified();

  c0->Write();

  cout << "2" << endl;


  TCanvas *c2 = new TCanvas("c2","HCal dx Sigma vs Nucleon p (MC)",1600,1200);
  //c2->Divide(2,1);
  c2->SetGrid();
  
  c2->cd();

  Double_t dbinp[nbin] = {0.};

  Double_t dxsig_p[nbin] = {0.};
  Double_t dxsig_n[nbin] = {0.};
  Double_t dysig_p[nbin] = {0.};
  Double_t dysig_n[nbin] = {0.};

  auto dxmg = new TMultiGraph();
  auto dymg = new TMultiGraph();
  auto allmg = new TMultiGraph();

  //loop over nucleons for slices, n==0 proton, n==1 neutron
  for(Int_t n=0; n<2; n++){
      
    Double_t p0;
    Double_t p1;
    TH2D *hdxvp;
    TH2D *hdyvp;
      
    if( n==0 ){ //if proton
      hdxvp = (TH2D*)(hdxvp_p->Clone("hdxvp"));
      hdyvp = (TH2D*)(hdyvp_p->Clone("hdyvp"));
    }else if( n==1 ){ //if neutron
      hdxvp = (TH2D*)(hdxvp_n->Clone("hdxvp"));
      hdyvp = (TH2D*)(hdyvp_n->Clone("hdyvp"));
    }else
      break;
    
    Double_t fitl;
    Double_t fith;

    TH1D *pbindxslice[nbin]; 
    TH1D *pbindyslice[nbin];

    for(Int_t b=0; b<nbin; b++){

      //Get expected mean from fit to hcale vs nucleon p
      Double_t p = b*p_step+pmin;
      Double_t fitp1exp = 0;
      Double_t fitp2exp = 0.10;

      pbindxslice[b] = hdxvp->ProjectionY(Form("pbindxslice_%d",b+1),b+1,b+1);
      pbindyslice[b] = hdyvp->ProjectionY(Form("pbindyslice_%d",b+1),b+1,b+1);

      fitl = fitp1exp - 3*fitp2exp;
      fith = fitp1exp + 3*fitp2exp;
	  
      TF1 *gausfitdx = new TF1("gausfitdx",fits::g_gfit,fitl,fith,3);
      gausfitdx->SetParameter(0,800);
      gausfitdx->SetParameter(1,fitp1exp);
      gausfitdx->SetParLimits(1,fitl,fith);
      gausfitdx->SetParameter(2,fitp2exp);
      gausfitdx->SetParLimits(2,0,0.25);

      TF1 *gausfitdy = new TF1("gausfitdy",fits::g_gfit,fitl,fith,3);
      gausfitdy->SetParameter(0,800);
      gausfitdy->SetParameter(1,fitp1exp);
      gausfitdy->SetParLimits(1,fitl,fith);
      gausfitdy->SetParameter(2,fitp2exp);
      gausfitdy->SetParLimits(2,0,0.25);

      pbindxslice[b]->Fit("gausfitdx","RBMQ");
      pbindxslice[b]->Draw();
      pbindyslice[b]->Fit("gausfitdy","RBMQ");
      pbindyslice[b]->Draw();
	
      dbinp[b] = p;
      if( n==0 ){
  	dxsig_p[b] = gausfitdx->GetParameter(2)*100; //convert to cm
  	dysig_p[b] = gausfitdy->GetParameter(2)*100;
  	pbindxslice[b]->SetTitle(Form("dxslice Loop:%d Np:%f Nuc:%d sig:%f",b,p,n,dxsig_p[b]));
  	pbindyslice[b]->SetTitle(Form("dyslice Loop:%d Np:%f Nuc:%d sig:%f",b,p,n,dysig_p[b]));   

      }else if( n==1 ){
  	dxsig_n[b] = gausfitdx->GetParameter(2)*100;
  	dysig_n[b] = gausfitdy->GetParameter(2)*100;
  	pbindxslice[b]->SetTitle(Form("dxslice Loop:%d Np:%f Nuc:%d sig:%f",b,p,n,dxsig_p[b]));
  	pbindyslice[b]->SetTitle(Form("dyslice Loop:%d Np:%f Nuc:%d sig:%f",b,p,n,dysig_p[b]));
   
      }
	
    } //end loop over bins

  } //end loop over nucleons

  cout << "3" << endl;

  //Draw graphs
  auto dxgrp = new TGraph(nbin,binp,dxsig_p);
  dxgrp->SetTitle("Proton X (RMS)");
  dxgrp->SetMarkerColor(kRed);
  dxgrp->SetMarkerStyle(21);
  dxgrp->SetMarkerSize(1);
  dxgrp->SetLineColor(kRed);
  dxgrp->SetLineWidth(0);
  dxmg->Add(dxgrp);
  allmg->Add(dxgrp);

  auto dxgrn = new TGraph(nbin,binp,dxsig_n);
  dxgrn->SetTitle("Neutron X (RMS)");
  dxgrn->SetMarkerColor(kBlue);
  dxgrn->SetMarkerStyle(20);
  dxgrn->SetMarkerSize(1);
  dxgrn->SetLineColor(kBlue);
  dxgrn->SetLineWidth(0);
  dxmg->Add(dxgrn);
  allmg->Add(dxgrn);

  dxmg->SetTitle("HCal X Res vs Nucleon p (MC)");
  dxmg->GetXaxis()->SetTitle("Nucleon p (GeV)");
  dxmg->GetYaxis()->SetTitle("dx sigma (m)");
  dxmg->Draw("AP");

  c2->BuildLegend();

  c2->Write();

  cout << "4" << endl;


  TCanvas *c3 = new TCanvas("c3","HCal dy Sigma vs Nucleon P (MC)",1600,1200);
  //c2->Divide(2,1);
  c3->SetGrid();

  c3->cd();

  //Draw graphs
  auto dygrp = new TGraph(nbin,dbinp,dysig_p);
  dygrp->SetTitle("Proton Y (RMS)");
  dygrp->SetMarkerColor(kRed);
  dygrp->SetMarkerStyle(25);
  dygrp->SetMarkerSize(1);
  dygrp->SetLineColor(kRed);
  dygrp->SetLineWidth(0);
  dymg->Add(dygrp);
  allmg->Add(dygrp);

  auto dygrn = new TGraphErrors(nbin,dbinp,dysig_n);
  dygrn->SetTitle("Neutron Y (RMS)");
  dygrn->SetMarkerColor(kBlue);
  dygrn->SetMarkerStyle(24);
  dygrn->SetMarkerSize(1);
  dygrn->SetLineColor(kBlue);
  dygrn->SetLineWidth(0);
  dymg->Add(dygrn);
  allmg->Add(dygrn);

  dymg->SetTitle("HCal Y Res vs Nucleon p (MC)");
  dymg->GetXaxis()->SetTitle("Nucleon p (GeV)");
  dymg->GetYaxis()->SetTitle("dy sigma (m)");
  dymg->Draw("AP");

  c3->BuildLegend();

  c3->Write();

  cout << "5" << endl;


  TCanvas *c4 = new TCanvas("c4","HCal dx and dy Sigma vs Nucleon P (MC)",1600,1200);
  //c2->Divide(2,1);
  c4->SetGrid();

  c4->cd();

  allmg->SetTitle("HCal Spatial Resolution (4x4 cluster)");
  allmg->GetXaxis()->SetTitle("Nucleon momentum (GeV/c)");
  allmg->GetYaxis()->SetTitle("X and Y Resolution (RMS) (cm)");
  allmg->Draw("AP");
  

  allmg->GetYaxis()->SetRangeUser(0,10);

  c4->Modified();

  c4->BuildLegend();

  c4->Write();


  cout << "6" << endl;


  fout->Write();

  // Send time efficiency report to console
  st->Stop();
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
