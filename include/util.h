#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <dirent.h>

#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"

#include "../include/crun.h"
#include "../include/econst.h"

namespace util {

  //Console
  string getDate();

  //General
  void readParam(std::string const_path,vector<Double_t> &param); //reads parameters from txt file
  void readParam(std::string const_path,std::string type,Int_t setsize,vector<Double_t> &param); //overload to read subset of data from sbs db-like file

  //HCal histograms
  TH2D *hhcalrowcol(std::string name); // returns hcal row/col 2d histo
  TH2D *hhcalxy(std::string name); // returns hcal xy 2d histo (data coordinates)
  TH2D *hhcalxy_mc(std::string name); // returns hcal xy 2d histo (mc coordinates)
  TH2D *hdxdy(std::string name); // returns hcal dxdy 2d histo (wide coordinates)
  TH1D *hhsamps(Int_t row,Int_t col,Int_t bins); // returns hcal waveform histogram

  // draws rectangular cut regions
  void drawarea(vector<Double_t> dimensions,      // a vector with extreme points
		Int_t lcolor,  // Default = 2 
		Int_t lwidth,  // Default = 4
		Int_t lstyle); // Default = 9


  //Kinematic histograms
  TH1D *hW2(std::string name);   // returns W histogram
  TH1D *hQ2(std::string name,   // returns Q2 histogram
		Int_t conf);   // SBS config

  //Functions to read csv files
  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   Int_t &nruns,                // No. of runs to analyze
		   Int_t sbsconf,               // SBS configuration
		   std::string target,        // target type
		   Int_t replay_pass,           // replay pass
		   Int_t verbose,               // verbosity
		   vector<crun> &corun);    // Output: Vector of crun structs

  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   Int_t &nruns,                // No. of runs to analyze
		   Int_t sbsconf,               // SBS configuration
		   std::string target,        // target type
		   Int_t replay_pass,           // replay pass
		   Int_t sbsmag,                // SBS magnet current (in %)
		   Int_t verbose,               // verbosity
		   vector<crun> &corun);    // Output: Vector of crun structs

  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   Int_t &nruns,                // No. of runs to analyze
		   Int_t sbsconf,               // SBS configuration
		   std::string target,        // target type
		   Int_t replay_pass,           // replay pass
		   Int_t sbsmag,                // SBS magnet current (in %)
		   Int_t bbmag,                 // BB magnet current (in %)
		   Int_t verbose,               // verbosity
		   vector<crun> &corun);    // Output: Vector of crun structs

  //Functions to read pd sort files
  Int_t LoadROOTTree(std::string path,          // ROOT file directory path
		   std::vector<crun> corun, // crun objects with run info
		   bool sort,                 // Sort by segments before parsing?
		   bool verbose,              // verbosity
		   TChain* &C);               // Output: TChain with data

  Int_t LoadROOTTree(std::string path,          // ROOT file directory path
		   crun corun,              // crun object with run info
		   bool sort,                 // Sort by segments before parsing?
		   bool verbose,              // verbosity
		   TChain* &C);               // Output: TChain with data
}

#endif
