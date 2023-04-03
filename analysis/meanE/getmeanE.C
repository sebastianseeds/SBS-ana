//sseeds 03.28.23 - script adapted from P.Datta designed to loop over all events in E tree and generate avg beam E per run/evt

#include <vector>
#include <iostream>

#include "../include/gmn.h"
#include "../src/jsonmgr.C"

Int_t getmeanE( Int_t kine = 9, Int_t nruns = -1, const char* targ = "lh2" ) 
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  
  // reading input config file
  JSONManager *jmgr = new JSONManager("../config/sgetmeanE.json");

  // set up the desired SBS configuration;
  std::string target = targ;
  Int_t pass = -1;
  if( kine==4 || kine==7 ){
    pass=0;
  }else if( kine==11 || kine==14 || kine==8 || kine==9 ){
    pass=1;
  }else{
    std::cout << "Error: enter a valid kinematic." << std::endl;
    return 0;
  }
  Int_t verb = 1;

  // read relevant run info from relevant good runlist spreadsheet
  std::string runsheet_dir = jmgr->GetValueFromKey_str("runsheet_dir");

  vector<crun> corun; util::ReadRunList(runsheet_dir,nruns,kine,target,pass,verb,corun); //set verbosity to 1
 
  // parsing trees depending on target type
  TChain *C = nullptr;
 
  std::string rootfile_dir = jmgr->GetValueFromKey_str(Form("rootfile_dir_%s_sbs%d",target.c_str(),kine));

  TString outFile; outFile = Form("epout/getmeanE_sbs%d_%s.csv",kine,target.c_str());
  ofstream outFile_data; outFile_data.open(outFile);
  outFile_data << "runnum," << "evnum," << "HALLA_p" << std::endl;
  
  for (Int_t irun=0; irun<nruns; irun++) {
    std::cout << "Analyzing run " << corun[irun].runnum << std::endl;

    std::string rfname = rootfile_dir + Form("/*%d*",corun[irun].runnum);
    C = new TChain("E");
    C->Add(rfname.c_str());

    // setting up ROOT tree branch addresses
    C->SetBranchStatus("*",0);
    // beam energy
    Double_t HALLA_p; rvars::setbranch(C, "HALLA_p", "", &HALLA_p);
    // global enent number
    long evnum; rvars::setbranch(C, "evnum", "", &evnum);
    
    // looping through EPICS events
    long nevent = 0, nevents = C->GetEntries();  
    while (C->GetEntry(nevent++)) {
      outFile_data << corun[irun].runnum << "," << evnum << "," << HALLA_p << std::endl;
    }

    // getting ready for the next run
    C->Reset();
  }

  return 0;
}
