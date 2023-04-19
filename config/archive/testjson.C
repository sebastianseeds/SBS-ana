//sseeds 04.03.23 - Script to test the added functionality of jsonmgr.C

#include <iostream>
#include "TStopwatch.h"
#include "../src/jsonmgr.C"
#include "../include/gmn.h"

void testjson( )
{ //main

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // reading input config file
  JSONManager *jmgr = new JSONManager("test.json");

  int val0 = jmgr->GetValueFromKey<int>("key1");

  std::cout << val0 << std::endl;

  double val1 = jmgr->GetValueFromKey<double>("key2");

  std::cout << val1 << std::endl;

  vector<int> vec0; jmgr->GetVectorFromKey<int>("keyvec1",vec0);

  for( Int_t el=0; el<vec0.size(); el++ ){
    std::cout << vec0[el] << std::endl;
  }
  
  int val2 = jmgr->GetValueFromSubKey<int>("keyset1","a");

  std::cout << val2 << endl;

  vector<int> vec1; jmgr->GetVectorFromSubKey<int>("keyset2","a",vec1);

  for( Int_t el=0; el<vec1.size(); el++ ){
    std::cout << vec1[el] << std::endl;
  }

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}//endmain
