//sseeds 03.28.23 - script designed to loop over all events in E tree and generate avg E per run/evt

#include <vector>
#include <iostream>

#include "../include/gmn.h"
#include "../src/jsonmgr.C"

Int_t tstE( Int_t kine = 9 ) 
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  
  // reading input config file
  JSONManager *jmgr = new JSONManager(Form("../config/SBS%d/sgetmeanE.json",kine));

  // set up the desired SBS configuration
  Int_t conf = jmgr->GetValueFromKey<Int_t>("SBS_config");
  std::string target = jmgr->GetValueFromKey_str("target");
  Int_t pass = jmgr->GetValueFromKey<Int_t>("replay_pass"); 
  Int_t verb = 0;

  // read relevant run info from relevant good runlist spreadsheet
  std::string runsheet_dir = jmgr->GetValueFromKey_str("runsheet_dir");
  Int_t nruns = jmgr->GetValueFromKey<Int_t>("Nruns_to_ana"); // # runs to analyze

  //cout << "NOPROBLEMO" << endl;

  cout << runsheet_dir << " " << nruns << " " << conf << " " << target << " " << pass << " " << verb << endl;

  //-----------test

  std::string fst = "/grl_"; 
  std::string mid = "_pass";
  std::string lst = ".csv";
  if (pass < 2) pass = 1; // single spreadsheet exists for pass 0 & 1
  std::string run_spreadsheet = runsheet_dir + fst + target + mid + std::to_string(pass) + lst;

  cout << run_spreadsheet << endl;

  // Reading the spreadsheet
  if (nruns < 0) nruns = 1e6;         // replay all runs if nruns < 0
  ifstream run_data; run_data.open(run_spreadsheet);
  string readline;
  if(run_data.is_open()){
    cout << "YOYUOUY" << endl;
  }
  //=====================testend

  //vector<crun> corun; util::ReadRunList(runsheet_dir,nruns,conf,target,pass,verb,corun); //set verbosity to 1
 
  //cout << corun[0] << endl;

  // // parsing trees depending on target type
  // TChain *C = nullptr;
  // std::string rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");

  // TString outFile; outFile = Form("epout/get_HALLA_p_prun_SBS%d_%s.csv",conf,target.c_str());
  // ofstream outFile_data; outFile_data.open(outFile);
  // outFile_data << "runnum," << "evnum," << "HALLA_p" << std::endl;
  
  // for (Int_t irun=0; irun<nruns; irun++) {
  //   std::cout << "Analyzing run " << corun[irun].runnum << std::endl;

  //   std::string rfname = rootfile_dir + Form("/*%d*",corun[irun].runnum);
  //   C = new TChain("E");
  //   C->Add(rfname.c_str());

  //   // setting up ROOT tree branch addresses
  //   C->SetBranchStatus("*",0);
  //   // beam energy
  //   Double_t HALLA_p; rvars::setbranch(C, "HALLA_p", "", &HALLA_p);
  //   // global enent number
  //   long evnum;   setrootvar::setbranch(C, "evnum", "", &evnum);
    
  //   // looping through EPICS events
  //   long nevent = 0, nevents = C->GetEntries();  
  //   while (C->GetEntry(nevent++)) {
  //     outFile_data << corun[irun].runnum << "," << evnum << "," << HALLA_p << std::endl;
  //   }

  //   // getting ready for the next run
  //   C->Reset();
  // }

  return 0;
}
