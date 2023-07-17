//sseeds 4.10.23 Short script to create and aggregate general physics plots over all kinematics after short parsing. Intended for use on parsed root files with output "P" tree from parse_sh.C

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
//const Int_t kIdx[6] = {4,7,11,14,8,9}; //indexing kinematics for processing, chronological order
const Int_t kIdx[6] = {8,11,4,9,7,14}; //indexing kinematics for processing, highest to lowest nevents
//const Int_t kIdx[6] = {7,14,9,4,11,8}; //indexing kinematics for processing, lowest to highest nevents

void aggplots(){ //main
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // outfile path
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string out_path = outdir_path + "/shortparse/psE_shE_plots_allkine.root";

  TFile *fout = new TFile( out_path.c_str(), "RECREATE" );

  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;

  //set up histograms
  TH1D *hpse[nkine];
  TH1D *hshe[nkine];

  //loop over kinematics
  for( Int_t k=0; k<nkine; k++ ){
    
    Int_t kine = kIdx[k];
    hpse[k] = new TH1D( Form("hpse_%d",kine), Form("BBCal Preshower E, SBS%d; GeV",kine), 240, 0., 1.2);
    hshe[k] = new TH1D( Form("hshe_%d",kine), Form("BBCal Shower E, SBS%d; GeV",kine), 350, 0., 3.5);

    C = new TChain("P"); //P designated tree name for output analysis tree
    std::string Pfile = outdir_path + Form("/shortparse/shortparse_sbs%d.root",kine);

    C->Add(Pfile.c_str());

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

    // loop over entries
    while (C->GetEntry(nevent++)) {
      
      if( nevent > 10000 ) break;

      Double_t pdone = (Double_t)nevent / (Double_t)nevents * 100.;

      cout << "Processing kinematic " << kine << ", " << Form("%0.2f",pdone) << "% complete.. \r";
      cout.flush();

      hpse[k]->Fill(pse);
      hshe[k]->Fill(she);

    }//end loop over entries

    C->Reset();

  }//end loop over kinematics

  //Build canvas to write histograms
  TCanvas *c1 = new TCanvas("c1","Preshower Energy, All Kinematics",1600,1200);
  c1->Divide(1,2);
  gStyle->SetPalette(53);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  c1->cd(1);

  for( Int_t k; k<nkine; k++ ){
    Int_t kine = kIdx[k];
    hpse[k]->SetLineColor(kPink+4-2*k);
    hpse[k]->SetLineWidth(2);
    hpse[k]->SetFillColor(kPink+4-2*k);
    hpse[k]->SetFillStyle(3001+k);
    if( k==0 ){
      hpse[k]->Draw();
    }else{
      hpse[k]->Draw("same");
    }
    //gPad->Modified(); gPad->Update();
  }

  // auto pselegend = new TLegend(0.1,0.6,0.3,0.9);
  // pselegend->SetTextSize( 0.03 );
  // for( Int_t k; k<nkine; k++ ){
  //   Int_t kine = kIdx[k];
  //   pselegend->AddEntry( hpse[k], Form("SBS%d",kine), "l");
  // }

  //pselegend->Draw();
  //gPad->Modified(); gPad->Update();

  // c1->Write();

  // TCanvas *c2 = new TCanvas("c2","Shower Energy, All Kinematics",1600,1200);
  // gStyle->SetPalette(53);
  // gStyle->SetOptFit(0);
  // gStyle->SetOptStat(0);
  // c2->cd();

  c1->cd(2);

  for( Int_t k; k<nkine; k++ ){
    Int_t kine = kIdx[k];
    hshe[k]->SetLineColor(kBlue+4-2*k);
    hshe[k]->SetLineWidth(2);
    hshe[k]->SetFillColor(kBlue+4-2*k);
    hshe[k]->SetFillStyle(3001+k);
    if( k==0 ){
      hshe[k]->Draw();
    }else{
      hshe[k]->Draw("same");
    }
    //gPad->Modified(); gPad->Update();
  }

  // auto shelegend = new TLegend(0.1,0.6,0.3,0.9);
  // shelegend->SetTextSize( 0.03 );
  // for( Int_t k; k<nkine; k++ ){
  //   Int_t kine = kIdx[k];
  //   shelegend->AddEntry( hshe[k], Form("SBS%d",kine), "l");
  // }

  // shelegend->Draw();
  // gPad->Modified(); gPad->Update();

  //c2->Write();

  // TCanvas *c1 = new TCanvas("c1","Preshower Energy, All Kinematics",1600,1200);
  // TCanvas *c2 = new TCanvas("c2","Shower Energy, All Kinematics",1600,1200);
  // auto c1leg = new TLegend(0.1,0.6,0.3,0.9);
  // auto c2leg = new TLegend(0.1,0.6,0.3,0.9);

  // for( Int_t k; k<nkine; k++ ){
  //   Int_t kine = kIdx[k];
    
  //   c1->cd();
  //   hpse[k]->SetLineColor(kPink+4-2*k);
  //   hpse[k]->SetLineWidth(2);
  //   hpse[k]->SetFillColor(kPink+4-2*k);
  //   hpse[k]->SetFillStyle(3001+k);
  //   hpse[k]->Draw("same");
  //   c1leg->AddEntry( hpse[k], Form("SBS%d",kine), "l");
  //   c1leg->Draw();
  //   gPad->Modified(); gPad->Update();

  //   c2->cd();
  //   hshe[k]->SetLineColor(kBlue+4-2*k);
  //   hshe[k]->SetLineWidth(2);
  //   hshe[k]->SetFillColor(kBlue+4-2*k);
  //   hshe[k]->SetFillStyle(3001+k);
  //   hshe[k]->Draw("same");
  //   c2leg->AddEntry( hshe[k], Form("SBS%d",kine), "l");
  //   c2leg->Draw();
  //   gPad->Modified(); gPad->Update();

  // }

  // c1->Write();
  // c2->Write();

  // Save and send time efficiency report to console
  fout->Write();
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}//end main
