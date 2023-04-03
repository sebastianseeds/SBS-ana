//sseeds 03.31.23 - test script to explore the utility functions

#include <vector>
#include <iostream>

#include "../../include/gmn.h"
#include "../../src/jsonmgr.C"

void testana( ) 
{
  
  std::string runsheet_dir = "/w/halla-scshelf2102/sbs/seeds/ana/data";
  Int_t nruns = -1;
  Int_t kine = 7;
  std::string target = "LH2";
  Int_t pass = 1;
  Int_t verb = 0;
  Int_t mag = 85;

  vector<crun> corun; 
  util::ReadRunList(runsheet_dir,nruns,kine,target,pass,verb,corun);
 
  std::string rootfile_dir;
  TFile *fout = new TFile( Form("test_sbs%d_%s.root",kine,target.c_str()), "RECREATE" );

  if( pass<2 ){
    if( kine==4 ) rootfile_dir = Form("/work/halla/sbs/sbs-gmn/pass0/SBS4/%s/rootfiles",target.c_str());
    if( kine==7 ) rootfile_dir = Form("/work/halla/sbs/sbs-gmn/pass0/SBS7/%s/rootfiles",target.c_str());
    if( kine==11 ) rootfile_dir = Form("/work/halla/sbs/sbs-gmn/pass1/SBS11/%s/rootfiles",target.c_str());
    if( kine==14 ) rootfile_dir = Form("/work/halla/sbs/sbs-gmn/pass1/SBS14/%s/rootfiles",target.c_str());
    if( kine==8 ) rootfile_dir = Form("/work/halla/sbs/sbs-gmn/pass1/SBS8/%s/rootfiles",target.c_str());
    if( kine==9 ) rootfile_dir = Form("/work/halla/sbs/sbs-gmn/pass1/SBS9/%s/rootfiles",target.c_str());
  }else{
    std::cout << "As of 3.31.23, the highest GMn replay pass is 1. Enter a valid pass." << endl;
    return;
  }
  
  SBSconfig config( kine, mag );
  cout << config;

  SBStune tune( kine, mag );
  cout << tune;

  std:string gcut = tune.Getglobcut();

  // parsing trees depending on target type
  TChain *C = nullptr;

  TH1D *hcale = new TH1D( "hcale", "", 400, 0, 5 );

  for (Int_t irun=0; irun<nruns; irun++) {
    std::cout << "Analyzing run " << corun[irun].runnum << std::endl;

    std::string rfname = rootfile_dir + Form("/*%d*",corun[irun].runnum);
    C = new TChain("T");
    C->Add(rfname.c_str());

    // setting up ROOT tree branch addresses
    C->SetBranchStatus("*",0);
    // HCal general
    Double_t hcalid, hcale, hcalx, hcaly, hcalr, hcalc, hcaltdc, hcalatime;
    std::vector<std::string> hcalvar = {"idblk","e","x","y","rowblk","colblk","tdctimeblk","atimeblk"};
    std::vector<void*> hcalvarlink = {&hcalid,&hcale,&hcalx,&hcaly,&hcalr,&hcalc,&hcaltdc,&hcalatime};
    rvars::setbranch(C, "sbs.hcal", hcalvar, hcalvarlink);


    // HCal
    Double_t sbs_hcal_e; 
    rvars::setbranch(C, "sbs.hcal.e", "", &sbs_hcal_e);
    // BB tracks
    Double_t bb_tr_p; 
    rvars::setbranch(C, "bb.tr.p", "", &bb_tr_p);
    
    // looping through EPICS events
    long nevent = 0, nevents = C->GetEntries();  
    while (C->GetEntry(nevent++)) {
      hcale->Fill( sbs_hcal_e );
      
    }

    // getting ready for the next run
    C->Reset();

  }

  fout->Write();
}
