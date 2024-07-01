//sseeds 5.10.23 - Script designed to run over all parsed data and return fitted location of neutron peak as a function of magnetic field setting per kinematic.

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

const Int_t nkine = 6; //total kinematics
const Int_t nmags = 6;
const Int_t kidx[6] = {4,7,11,14,8,9};
const Int_t midx[6] = {0,30,50,70,85,100}; //magnetic field index for this analysis
const Double_t dxmax = 2.;
const Double_t dxmin = -2.;
const Double_t dxbinfac = 100.;
const Double_t dxfitmax = 0.2;
const Double_t dxfitmin = -0.1;
const Int_t fset_ld2[6][4] = {{0,30,50,-1},
			      {85,-1,-1,-1},
			      {0,100,-1,-1},
			      {70,-1,-1,-1},
			      {0,50,70,100},
			      {70,-1,-1,-1}}; //All field settings for all kinematics LH2. -1 indicates no setting

void neutron_peak_diagnostic(){ //main
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // outfile path
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string npeak_path = outdir_path + "/shortparse/neutron_peak_diagnostic_out.root";

  TFile *fout = new TFile( npeak_path.c_str(), "RECREATE" );

  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;

  //set up histograms
  TH1D *hdx[nkine][nmags];
  TH2D *hdxvindex = new TH2D( "hdxvindex", "HCal dx vs kinematic.magset; index", 150, 0, 15, dxbinfac*(dxmax-dxmin), dxmin, dxmax);

  //loop over kinematics
  for( Int_t k=0; k<nkine; k++ ){

    Int_t kine = kidx[k];

    for( Int_t nm=0; nm<nmags; nm++ )
      hdx[k][nm] = new TH1D( Form("hdx_%d_%d",kine,midx[nm]), Form("HCal dx SBS%d mag%d; m",kine,midx[nm]), dxbinfac*(dxmax-dxmin), dxmin, dxmax);

    C = new TChain("P"); //P designated tree name for output analysis tree
    std::string Pfile = Form("outfiles/shortparse_sbs%d.root",kine);
    C->Add(Pfile.c_str());

    //Elastics and LD2 selected
    string cut = "failedglobal==0&&failedaccmatch==0&&failedcoin==0&&abs(dy)<0.2&&abs(W2-1)<0.4&&tar==1";
    TCut GCut = cut.c_str();
    TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );

    // setting up ROOT tree branches
    Int_t mag, run, tar, failedglobal, failedaccmatch, failedcoin;
    Double_t dx, dy, W2, Q2, nu, hcale, pse, she, ep, eoverp, hcalatime, hodotmean, thetapq_p, thetapq_n;

    //C->SetMakeClass(1);

    // adding all branches from P tree to be ready for additional plots
    C->SetBranchStatus("*",1);

    C->SetBranchStatus( "dx", 1 );
    C->SetBranchStatus( "dy", 1 );
    C->SetBranchStatus( "W2", 1 );
    C->SetBranchStatus( "Q2", 1 );
    C->SetBranchStatus( "nu", 1 );
    C->SetBranchStatus( "hcale", 1 );
    C->SetBranchStatus( "pse", 1 );
    C->SetBranchStatus( "she", 1 );
    C->SetBranchStatus( "ep", 1 );
    C->SetBranchStatus( "eoverp", 1 );
    C->SetBranchStatus( "hcalatime", 1 );
    C->SetBranchStatus( "hodotmean", 1 );
    C->SetBranchStatus( "thetapq_p", 1 );
    C->SetBranchStatus( "thetapq_n", 1 );
    C->SetBranchStatus( "mag", 1 );
    C->SetBranchStatus( "run", 1 );
    C->SetBranchStatus( "tar", 1 );
    C->SetBranchStatus( "failedglobal", 1 );
    C->SetBranchStatus( "failedaccmatch", 1 );
    C->SetBranchStatus( "failedcoin", 1 );

    C->SetBranchAddress( "dx", &dx );
    C->SetBranchAddress( "dy", &dy );
    C->SetBranchAddress( "W2", &W2 );
    C->SetBranchAddress( "Q2", &Q2 );
    C->SetBranchAddress( "nu", &nu );
    C->SetBranchAddress( "hcale", &hcale );
    C->SetBranchAddress( "pse", &pse );
    C->SetBranchAddress( "she", &she );
    C->SetBranchAddress( "ep", &ep );
    C->SetBranchAddress( "eoverp", &eoverp );
    C->SetBranchAddress( "hcalatime", &hcalatime );
    C->SetBranchAddress( "hodotmean", &hodotmean );
    C->SetBranchAddress( "thetapq_p", &thetapq_p );
    C->SetBranchAddress( "thetapq_n", &thetapq_n );
    C->SetBranchAddress( "mag", &mag );
    C->SetBranchAddress( "run", &run );
    C->SetBranchAddress( "tar", &tar );
    C->SetBranchAddress( "failedglobal", &failedglobal );
    C->SetBranchAddress( "failedaccmatch", &failedaccmatch );
    C->SetBranchAddress( "failedcoin", &failedcoin );

    // event indices
    Long64_t nevent = 0, nevents = C->GetEntries(); 
    Int_t treenum = 0, currenttreenum = 0;

    // loop over entries
    while (C->GetEntry(nevent++)) {
   
      // if( nevent>10000 )
      // 	break;

      Double_t pdone = (Double_t)nevent / (Double_t)nevents * 100.;

      cout << "Processing kinematic " << kine << ", " << Form("%0.2f",pdone) << "% complete.. \r";
      
      cout.flush();

      //Single-loop globalcut on elastics
      currenttreenum = C->GetTreeNumber();
      if( nevent == 1 || currenttreenum != treenum ){
      	treenum = currenttreenum; 
      	GlobalCut->UpdateFormulaLeaves();
      }
      bool failgcut = GlobalCut->EvalInstance(0) == 0;
      
      if( failgcut ) continue;

      Int_t mag_index = -1;
      for( Int_t nm=0; nm<nmags; nm++ )
	if( mag==midx[nm] ) mag_index = nm;

      hdx[k][mag_index]->Fill(dx);

      //Set up kine.mag sub index for single plot comparison
      Double_t magsubidx = (double)mag*0.01;
      if( mag==100 )
	magsubidx = 0.9;

      Double_t tdidx = (double)kine+magsubidx;

      hdxvindex->Fill(tdidx,dx);

    }//end loop over entries

    C->Reset();

  }//end loop over kinematics
  hdxvindex->Write();

  cout << "Got here" << endl;

  // Loop over all histos and fit neutron peaks
  Double_t meandx[nkine][nmags] = {-10.};
  Double_t sigdx[nkine][nmags] = {0.};
  Double_t cell[nkine][nmags] = {0.};
  Double_t cellerr[nkine][nmags] = {0.};

  Double_t meandx_m[nmags][nkine] = {-10.};
  Double_t sigdx_m[nmags][nkine] = {0.};
  Double_t cell_m[nmags][nkine] = {0.};
  Double_t cellerr_m[nmags][nkine] = {0.};

  Double_t cell_a[nkine*nmags] = {-10.};
  Double_t meandx_a[nkine*nmags] = {-10.}; 
  Double_t sigdx_a[nmags*nkine] = {0.};
  Double_t cellerr_a[nkine*nmags] = {0.};


  for( Int_t k=0; k<nkine; k++ ){

    Int_t kine = kidx[k];
    
    for( Int_t nm=0; nm<nmags; nm++ ){
      
      Int_t mag = midx[nm];

      cell[k][nm] = midx[nm];
      cell_m[nm][k] = kine;

      Double_t magsubidx = (double)mag*0.01;
      if( mag==100 )
	magsubidx = 0.9;

      Double_t tdidx = (double)kine+magsubidx;

      Int_t caidx = k*nmags+nm;

      cell_a[caidx] = tdidx;

      if( hdx[k][nm]->GetEntries() < 100 ) //Only evaluate fit for field settings where data was taken
	continue;

      TF1 *gfit = new TF1("gfit",fits::g_gfit,dxfitmin,dxfitmax,3);
      gfit->SetLineWidth(4);
      gfit->SetParameter(0,1000);
      gfit->SetParameter(1,0);
      gfit->SetParLimits(1,dxfitmin,dxfitmax);
      gfit->SetParameter(2,0.2);
      gfit->SetParLimits(2,0.01,1.0);

      hdx[k][nm]->Fit("gfit","RBM");

      meandx[k][nm] = gfit->GetParameter(1);
      meandx_m[nm][k] = gfit->GetParameter(1);
      meandx_a[caidx] = gfit->GetParameter(1);
      sigdx[k][nm] = gfit->GetParameter(2);
      sigdx_m[nm][k] = gfit->GetParameter(2);
      sigdx_a[caidx] = gfit->GetParameter(2);

    }
  }
      
  TGraphErrors *gnp_kine[nkine];
  TGraphErrors *gnp_mag[nmags];
  TGraphErrors *gnp = new TGraphErrors( nmags*nkine, cell_a, meandx_a, cellerr_a, sigdx_a ); 

  gnp->GetXaxis()->SetLimits(0,15);  
  gnp->GetYaxis()->SetLimits(dxmin,dxmax);
  gnp->SetTitle( "Neutron Peak Mean; kine.mag; dx (m)" );
  gnp->GetXaxis()->SetTitle("SBS Field");
  gnp->GetYaxis()->SetTitle("dx (m)");
  gnp->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
  gnp->Draw();
  gnp->Write("gnp");

  for( Int_t k=0; k<nkine; k++ ){

    Int_t kine = kidx[k];
    
    gnp_kine[k] = new TGraphErrors( nmags, cell[k], meandx[k], cellerr[k], sigdx[k] ); 
    gnp_kine[k]->GetXaxis()->SetLimits(0,100);  
    gnp_kine[k]->GetYaxis()->SetLimits(dxmin,dxmax);
    gnp_kine[k]->SetTitle(Form("Neutron Peak Mean, kine: %d",kine ) );
    gnp_kine[k]->GetXaxis()->SetTitle("SBS Field");
    gnp_kine[k]->GetYaxis()->SetTitle("dx (m)");
    gnp_kine[k]->SetMarkerStyle(20); // idx 20 Circles, idx 21 Boxes
    gnp_kine[k]->Draw();
    gnp_kine[k]->Write(Form("gnp_kine%d",kine));
  }

  for( Int_t m=0; m<nmags; m++ ){

    Int_t mag = midx[m];
    
    gnp_mag[m] = new TGraphErrors( nkine, cell_m[m], meandx_m[m], cellerr_m[m], sigdx_m[m] ); 
    gnp_mag[m]->GetXaxis()->SetLimits(0,10);  
    gnp_mag[m]->GetYaxis()->SetLimits(dxmin,dxmax);
    gnp_mag[m]->SetTitle(Form("Neutron Peak Mean, mag: %d%%",mag ) );
    gnp_mag[m]->GetXaxis()->SetTitle("SBS Field");
    gnp_mag[m]->GetYaxis()->SetTitle("dx (m)");
    gnp_mag[m]->SetMarkerStyle(21); // idx 20 Circles, idx 21 Boxes
    gnp_mag[m]->Draw();
    gnp_mag[m]->Write(Form("gnp_mag%d",mag));
  }

     

  // Save and send time efficiency report to console
  fout->Write();
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}//end main
