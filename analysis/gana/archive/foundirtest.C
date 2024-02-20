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

int foundirtest( int kine = 8, int pass = 2 ) {

  JSONManager *jmgr = new JSONManager("../../config/parse.json");

  std::string rootfile_dir_lh2 = jmgr->GetValueFromSubKey_str( Form("rootfile_dir_lh2_p%d",pass), Form("sbs%d",kine) );
  std::string rootfile_dir_ld2 = jmgr->GetValueFromSubKey_str( Form("rootfile_dir_ld2_p%d",pass), Form("sbs%d",kine) );

  //Necessary for loop over directories in sbs8
  std::vector<TString> directories_h = {
    rootfile_dir_lh2 + "/SBS0percent",
    rootfile_dir_lh2 + "/SBS100percent",
    rootfile_dir_lh2 + "/SBS50percent",
    rootfile_dir_lh2 + "/SBS70percent_part1",
    rootfile_dir_lh2 + "/SBS70percent_part2",
    rootfile_dir_lh2 + "/SBS70percent_part3"
  };

  std::vector<TString> directories_d = {
    rootfile_dir_ld2 + "/SBS0percent",
    rootfile_dir_ld2 + "/SBS100percent",
    rootfile_dir_ld2 + "/SBS50percent",
    rootfile_dir_ld2 + "/SBS70percent_part1",
    rootfile_dir_ld2 + "/SBS70percent_part2",
    rootfile_dir_ld2 + "/SBS70percent_part3",
    rootfile_dir_ld2 + "/SBS70percent_part4"
  };

  TString pattern_h = "*13574*";  // Replace with the wildcard pattern you're searching for

  TString foundDir_h = util::FindFileInDirectories(pattern_h, directories_h);
  if (!foundDir_h.IsNull()) {
    std::cout << "File found in: " << foundDir_h << std::endl;
  } else {
    std::cout << "File not found in: " << foundDir_h << std::endl;
  }

  TString pattern_d = "*13570*";  // Replace with the wildcard pattern you're searching for

  TString foundDir_d = util::FindFileInDirectories(pattern_d, directories_d);
  if (!foundDir_d.IsNull()) {
    std::cout << "File found in: " << foundDir_d << std::endl;
  } else {
    std::cout << "File not found in: " << foundDir_d << std::endl;
  }

  return 0;
}
