//sseeds 4.25.23 - plot position comparison with MC, sbs 4

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

const Int_t nkine = 6; //total kinematics
const Int_t kIdx[6] = {4,7,11,14,8,9}; //indexing kinematics for processing
//const Int_t pIdx[6] = {4,5,0,3,2,1};
const Int_t pIdx[6] = {4,0,5,3,2,1};
const Int_t pIdxmc[6] = {2,1,5,4,3,0};
const Double_t fmean[6] = {0.15,0.1,0.15,0.28,0.38,0.54};

//MAIN. (no args)
void plot_sbs4_comparison( ){
  
  TCanvas *c1 = new TCanvas("c1","SBS4 MC vs Data dx/dy",1600,1200);
  c1->Divide(2,1);

  TFile *fmc[nkine];
  TH1D *hmc[nkine];

  TFile *fmc = TFile::Open("outfiles/mcrep_pos_sbs4_epm3.root");
  TFile *fdata = TFile::Open("../gana/outfiles/shortparse_sbs4.root");

  auto hdx = new THStack("hdx","HCal dx");
  auto hdy = new THStack("hdy","HCal dy");

  Double_t fitl;
  Double_t fith;
  Double_t fits = 0.5;
  for( Int_t i=0; i<nkine; i++ ){
    fitl = fmean[i]-fits;
    fith = fmean[i]+fits;

    fmc[i] = TFile::Open(Form("../MC/output/MChcalE_sbs%d.root",kIdx[i]));
    hmc[i]= (TH1D*)fmc[i]->Get("hHCALe_clus");
    //hmc[i]->GetXaxis()->SetRange(0,0.17);
    hmc[i]->SetLineWidth(0);
    
    hmc[i]->Fit("gaus","0","",fitl,fith);
    fitmc[i] = hmc[i]->GetFunction("gaus");
    fitmcp1[i] = fitmc[i]->GetParameter(1);

    hmc[i]->SetTitle(Form("SBS%d, peak %0.2f GeV",kIdx[i],fitmcp1[i]));

  }

  c1->cd(1);

  for( Int_t i=0; i<nkine; i++ ){
    //Int_t j = pIdxmc[i];
    Int_t j = pIdx[i];

    hsmc->Add(hmc[j]);
  }

  hsmc->Draw("pfc nostack");

  hsmc->GetXaxis()->SetTitle("E_{dep} (GeV)");
  hsmc->GetXaxis()->SetTitleOffset(1.4);

  c1->Modified();
  
  gPad->BuildLegend(0.45,0.65,0.89,0.89,"");


  gStyle->SetOptStat(0);
  //gStyle->SetPalette(kOcean);
  gStyle->SetPalette(kCividis);
  //Get MC samp frac all kine
  TFile *f[nkine];
  TH1D *h[nkine];
  TF1 *fit[nkine];
  Double_t fitp1[nkine];

  auto hs = new THStack("hs","HCal Cluster Energy (Calibrated Data)");

  for( Int_t i=0; i<nkine; i++ ){
    fitl = fmean[i]-fits;
    fith = fmean[i]+fits;

    f[i] = TFile::Open(Form("qreplay_sbs%d.root",kIdx[i]));
    h[i]= (TH1D*)f[i]->Get("hHCALe");
    //h[i]->GetXaxis()->SetRange(0,0.17);
    h[i]->SetLineWidth(0);
    
    h[i]->Fit("gaus","0","",fitl,fith);
    fit[i] = h[i]->GetFunction("gaus");
    fitp1[i] = fit[i]->GetParameter(1);

    h[i]->SetTitle(Form("SBS%d, peak %0.2f GeV",kIdx[i],fitp1[i]));


  }

  c1->cd(2);
  

  for( Int_t i=0; i<nkine; i++ ){
    Int_t j = pIdx[i];

    hs->Add(h[j]);
  }
  
  hs->Draw("pfc nostack");

  hs->GetXaxis()->SetTitle("E_{dep} (GeV)");
  hs->GetXaxis()->SetTitleOffset(1.4);

  c1->Modified();
  
  gPad->BuildLegend(0.45,0.65,0.89,0.89,"");

  //c1->SaveAs("/u/home/sseeds/Plots/MC_SF_allkine.png");

}
