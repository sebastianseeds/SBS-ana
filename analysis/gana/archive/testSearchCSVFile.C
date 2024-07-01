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


std::string fileloc = "/lustre19/expphy/volatile/halla/sbs/jboyd/simulation/out_dir/MC_REPLAY_OUT_DIR/hist/jb_FARM_SIMC_gmn_SBS9_LD2_proton_mag70port0656_1200uA_elas_250k_06_01_2024_job644.hist";

//MAIN
void testSearchCSVFile( )
{   

  int Ntried;
  Ntried = (int)util::searchSimcHistFile("Ntried", fileloc);
  double genvol;
  genvol = (double)util::searchSimcHistFile("genvol", fileloc);  
  double lumi;
  lumi = (double)util::searchSimcHistFile("luminosity", fileloc);

  cout << "Ntried: " << Ntried << ", genvol: " << genvol << ", luminosity: " << lumi << endl;

  return 0;

}
