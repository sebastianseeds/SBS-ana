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

const Int_t maxClus = 35; //larger than max per tree
const Int_t maxBlk = 25;
const Double_t W2max = 3.0; //very large compared to nucleon mass
const Double_t coin_sigma_factor = 5.; //Wide coincidence timing cut
const Double_t nsig_step = 0.1; //number of sigma to step through in dx until fiducial cut failure written to output tree

bool norm_override = false;
double Ntried_override = 100000;
double luminosity_override = 3.8475e+36;
double genvol_override = 12.566;

//Specific wide cut for all parsing
const std::string gcut = "bb.ps.e>0.1&&abs(bb.tr.vz[0])<0.12";

//MAIN
void testcheck_mcfiles( Int_t kine=4, Int_t mag = 30, const char *replay_type = "", bool verbose=false )
{   

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // Set up exception to MC file structure
  std::string rtype = replay_type;
  std::string rootfile_type = "";
  bool jboyd = false;
  if( rtype.compare("jboyd")==0 ){
      rootfile_type = "_jboyd";
      jboyd = true;
  }

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/gmn_mc.json");

  //Get rootfiles and metadata accounting for exceptions
  //All MC for analysis is LD2
  std::string rootfile_dir = jmgr->GetValueFromSubKey_str( Form("filedir_sbs%d%s",kine,rootfile_type.c_str()), Form("%dp",mag) );

  if( jboyd )
    rootfile_dir += Form("simc/SBS%d/",kine);

  std::string histfile_dir = rootfile_dir;

  if( jboyd )
    histfile_dir += "hist/";
  else
    histfile_dir += "simcout/";

  std::string partialName_p = jmgr->GetValueFromSubKey_str( Form("p_string_sbs%d",kine), Form("%dp",mag) );
  std::string partialName_n = jmgr->GetValueFromSubKey_str( Form("n_string_sbs%d",kine), Form("%dp",mag) );

  //set up default parameters for all analysis
  double binfac = 400.;

  // outfile path
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string parse_path = outdir_path + Form("/parse/parse_mc_sbs%d_barebones.root",kine);

  //set up output files
  TFile *fout = new TFile( parse_path.c_str(), "RECREATE" );

  //Check neutron root/hist files
  std::vector<std::string> rootFileNames_n;

  std::map<int, std::pair<std::string, std::vector<float>>> csvData_n;

  util::SyncFilesWithCsv(histfile_dir, rootfile_dir, partialName_n, rootFileNames_n, csvData_n);

  // Write verified data to a file
  //std::ofstream outFile("verified_data.txt");
  for (const auto& [jobId, dataPair] : csvData) {
    cout << "Neutron Job ID: " << jobId << endl;
    cout << "  rootfile name: " << dataPair.first << endl;
    cout << "  data:";
    for (auto& val : dataPair.second) {
      cout << val << " ";
    }
    cout << endl << endl;
    //marker++;
  }

  return 0;

}
