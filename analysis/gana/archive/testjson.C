//sseeds test script to verify function of json configuration load
#include <TSystem.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TCollection.h>
#include <TRegexp.h>
#include <TString.h>
#include <iostream>
#include <vector>
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

void testjson( int kine = 9, int mag = 70, int pass = 2 ) {

  //access json file
  JSONManager *jmgr = new JSONManager("../../config/syst.json");

  //access first key and read contents
  cout << endl << "Accessing first key " << Form("post_cuts_p%d",pass) << ", subkey " << Form("sbs%d_%d",kine,mag) << endl;

  std::string postcuts = jmgr->GetValueFromSubKey_str( Form("post_cuts_p%d",pass), Form("sbs%d_%d",kine,mag) );

  cout << "Loaded tight cuts: " << postcuts << endl;

  std::vector<std::string> cuts = util::parseCuts(postcuts);

  cout << "Parsed cuts: " << endl;
  for( size_t i=0; i<cuts.size(); ++i ) 
    cout << cuts[i] << endl;

  //access second key and read contents
  cout << endl << "Accessing second key " << Form("post_branches_p%d",pass)  << endl;

  std::string postbranches = jmgr->GetValueFromKey_str( Form("post_branches_p%d",pass) );

  cout << "Loaded branches: " << postbranches << endl;

  std::vector<std::string> branches = util::parseCuts(postbranches);

  cout << "Parsed branches: " << endl;
  for( size_t i=0; i<branches.size(); ++i ) 
    cout << branches[i] << endl;

  //access second key and read contents
  cout << endl << "Accessing third key " << Form("cut_limits_p%d",pass)  << endl;

  std::vector<double> cutlimits; jmgr->GetVectorFromSubKey<Double_t>( Form("cut_limits_p%d",pass), Form("sbs%d_%d",kine,mag), cutlimits );

  int limitidx = 0;
  cout << "Cut limits: " << endl;
  for( size_t i=0; i<branches.size(); ++i ){ 
    cout << "For branch " << branches[i] << ": ";
    cout << cutlimits[limitidx] << " to ";
    limitidx++;
    cout << cutlimits[limitidx] << endl;
    limitidx++;
  }

  return;
}
