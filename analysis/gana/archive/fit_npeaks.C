//sseeds 5.11.23: backup script to tune neutron peaks from output histograms neutron_peak_diagnostic.C

#include <vector>
#include <iostream>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <string>
#include <TChain.h>
#include <TSystem.h>
#include <TStopwatch.h>
#include <unistd.h>
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TCut.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TLegend.h"
#include "TVector3.h"
#include "TGraphErrors.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TF1.h"
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

const Int_t nkine = 6; //total kinematics
const Int_t mmags = 4; //maximum field settings per kinematic
const Int_t nmag = 12; //total magnetic field setting / kinematic combos
const Double_t fitl = -0.1;
const Double_t fitu = 0.2;
const Int_t kidx[6] = {4,7,11,14,8,9};

const Int_t fset[6][4] = {{0,30,50,-1},
			  {85,-1,-1,-1},
			  {0,100,-1,-1},
			  {70,-1,-1,-1},
			  {0,50,70,100},
			  {70,-1,-1,-1}}; //All field settings for all kinematics LH2. -1 indicates no setting

void fit_npeaks( Int_t kine = 4 ){

  // outfile path
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string file_path = outdir_path + "/shortparse/neutron_peak_diagnostic_out.root";
  std::string out_path = outdir_path + "/shortparse/npeakfits_allkine.root";

  TCanvas *c1 = new TCanvas("c1","LD2 dx Kinematics 4,7,11",2200,1200);
  c1->Divide(3,4);
  c1->SetGrid();

  TCanvas *c2 = new TCanvas("c2","LD2 dx Kinematics 14,8,9",2200,1200);
  c2->Divide(3,4);
  c2->SetGrid();

  gStyle->SetOptStat(0);
  gStyle->SetPalette(53);
  gStyle->SetOptFit();
  gStyle->SetEndErrorSize(0);

  
  TH1D *hnull = new TH1D("hnull","hnull",100,-1,1);

  TH1D *hdx[nkine*mmags];

  vector<Double_t> mean;
  vector<Double_t> sig;
  vector<Double_t> cell;
  vector<Double_t> cellerr;

  for( Int_t k=0; k<nkine; k++ ){
    Int_t kine = kidx[k];
    
    for( Int_t m=0; m<mmags; m++ ){
      Int_t kmidx = k*mmags+m;

      Int_t mfset = fset[k][m];

      Double_t magsubidx = (double)mfset*0.01;
      if( mfset==100 )
	magsubidx = 0.9;

      Double_t tdidx = (double)kine+magsubidx;
      
      //index through canvas
      c1->cd(kmidx+1);
      if( kmidx>=12 )
	c2->cd(kmidx-11);
      
      //get the histogram
      TFile *f1 = TFile::Open(file_path.c_str());
      hdx[kmidx] = (TH1D*)f1->Get(Form("hdx_%d_%d",kine,mfset));

      if( mfset==-1 || hdx[kmidx]->GetEntries()<100 ){ //avoid failed fits and empty points
	hnull->Draw();
	continue;
      }

      cell.push_back(tdidx);
      cellerr.push_back(0.);
      
      hdx[kmidx]->Fit("gaus","Q","",fitl,fitu);
      hdx[kmidx]->Draw();

      TF1 *gfit;
      gfit = hdx[kmidx]->GetFunction("gaus");
      mean.push_back( gfit->GetParameter(1) );
      sig.push_back( gfit->GetParameter(2) );
      hdx[kmidx]->SetTitle( Form("k:%d s:%d m:%0.2f s:%0.2f",kine,mfset,mean.back(),sig.back()) );
      hdx[kmidx]->Write();

    }
    
  }

  TCanvas *c3 = new TCanvas("c3","All Kine dx n fits",2200,1200);
  c3->SetGrid();
  c3->cd();

  //Make graphs with errors for reporting
  TGraphErrors *gnp = new TGraphErrors( cell.size(), &(cell[0]), &(mean[0]), &(cellerr[0]), &(sig[0]) );

  gnp->SetTitle("dx vs kine.mag");
  gnp->GetXaxis()->SetTitle("kine.mag");
  gnp->GetYaxis()->SetTitle("dx (m)");
  gnp->SetMarkerStyle(33); // idx 20 Circles, idx 21 Boxes
  gnp->SetMarkerColor(kSpring-6);
  gnp->SetMarkerSize(3);
  gnp->SetLineWidth(3);
  gnp->SetLineColor(kSpring-6);

  gnp->Draw();

  TFile *fout = new TFile( out_path.c_str(), "RECREATE" );
  gnp->Write("gnp");

  // for( Int_t i=0; i<nkine*mmags; i++ )
  //   hdx[i]->Write();

  fout->Write();
}
