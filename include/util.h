#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <dirent.h>
#include <filesystem>

#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"

#include "../include/crun.h"
#include "../include/econst.h"

namespace fs = std::filesystem;

namespace util {

  //Console
  string getDate();

  //send progress in for loop to console
  void updateProgress(Int_t currentIteration, 
		      Int_t totalIterations);

  //send progress in for loop to console and print a label (for nests)
  void updateProgress(Int_t currentIteration, 
		      Int_t totalIterations,
		      std::string label);
  
  //send progress in nested for loop to console
  void updateProgress(Int_t currentIteration_outer, 
		      Int_t currentIteration_inner, 
		      Int_t totalIterations_outer, 
		      Int_t totalIterations_inner,
		      std::string label);

  //General
  void readParam(std::string const_path,vector<Double_t> &param); //reads parameters from txt file
  void readParam(std::string const_path,std::string type,Int_t setsize,vector<Double_t> &param); //overload to read subset of data from sbs db-like file

  //HCal histograms
  TH2D *hhcalrowcol(std::string name); // returns hcal row/col 2d histo
  TH2D *hhcalxy(std::string name); // returns hcal xy 2d histo (data coordinates)
  TH2D *hhcalxy_mc(std::string name); // returns hcal xy 2d histo (mc coordinates)
  TH2D *hdxdy(std::string name); // returns hcal dxdy 2d histo (wide coordinates)
  TH1D *hhsamps(Int_t row,Int_t col,Int_t bins); // returns hcal waveform histogram

  // checks if point is within proton/neutron spot
  bool Nspotcheck(Double_t dy,               //hcaly
		  Double_t dx,               //hcalx
		  Double_t dy_mean,          //elastic peak location in dy
		  Double_t dx_mean,          //elastic peak location in dx
		  Double_t dy_sigma,         //elastic peak sigma in dy
		  Double_t dx_sigma,         //elastic peak sigma in dx
		  Double_t rotationAngle);   //rotation angle in radians, if ever applicable
  
  // draws rectangular cut regions
  void drawarea(vector<Double_t> dimensions,      // a vector with extreme points
		Int_t lcolor,  // Default = 2 
		Int_t lwidth,  // Default = 4
		Int_t lstyle); // Default = 9


  //Kinematic histograms
  TH1D *hW2(std::string name);  // returns W histogram
  TH1D *hQ2(std::string name,   // returns Q2 histogram
	    Int_t conf);        // SBS config

  //Cluster finding score-based algorithm
  double assignScore( Double_t val_1,   //cluster energy
		      Double_t val_2,   //cluster coin time (via primary block)
		      Double_t max1,    //max cluster energy among clusters
		      const std::vector<Double_t> &gfit );  //coin atime gaus fit parameters

  //Functions to read csv files
  void ReadRunList(std::string runsheet_dir,    // Dir. path containing CSV files with run info
		   Int_t &nruns,                // No. of runs to analyze
		   Int_t sbsconf,               // SBS configuration
		   std::string target,          // target type
		   Int_t replay_pass,           // replay pass
		   Int_t verbose,               // verbosity
		   vector<crun> &corun);        // Output: Vector of crun structs

  void ReadRunList(std::string runsheet_dir,    // Dir. path containing CSV files with run info
		   Int_t &nruns,                // No. of runs to analyze
		   Int_t sbsconf,               // SBS configuration
		   std::string target,          // target type
		   Int_t replay_pass,           // replay pass
		   Int_t sbsmag,                // SBS magnet current (in %)
		   Int_t verbose,               // verbosity
		   vector<crun> &corun);        // Output: Vector of crun structs

  void ReadRunList(std::string runsheet_dir,    // Dir. path containing CSV files with run info
		   Int_t &nruns,                // No. of runs to analyze
		   Int_t sbsconf,               // SBS configuration
		   std::string target,          // target type
		   Int_t replay_pass,           // replay pass
		   Int_t sbsmag,                // SBS magnet current (in %)
		   Int_t bbmag,                 // BB magnet current (in %)
		   Int_t verbose,               // verbosity
		   vector<crun> &corun);        // Output: Vector of crun structs

  //Functions to read pd sort files
  Int_t LoadROOTTree(std::string path,          // ROOT file directory path
		     std::vector<crun> corun,   // crun objects with run info
		     bool sort,                 // Sort by segments before parsing?
		     bool verbose,              // verbosity
		     TChain* &C);               // Output: TChain with data

  Int_t LoadROOTTree(std::string path,          // ROOT file directory path
		     crun corun,                // crun object with run info
		     bool sort,                 // Sort by segments before parsing?
		     bool verbose,              // verbosity
		     TChain* &C);               // Output: TChain with data

  //Functions for simc and gmn analysis
  Double_t searchSimpleCSVForValue(const std::string& filePath,   //Path to the .csv file
			     const std::string& key);       //key

  Double_t searchSimcHistFile(const TString &searchPattern, // moniker to search file for ("luminosity","genvol", "Ntried")
			      const TString &filename);     // full path to name

  std::vector<Double_t> fitGaussianAndGetFineParams(TH1D* hist,       // 1D histo passed for fitting
						    Double_t sig,     // Estimated sigma for course fit
						    Double_t low,     // Lower bound for first fit, default first bin
						    Double_t high);   // Upper bound for first fit, default last bin

  std::vector<Double_t> fitSkewedGaussianAndGetFineParams(TH1D* hist,       // 1D histo passed for fitting
							  Double_t sig,     // Estimated sigma for course fit
							  Double_t low,     // Lower bound for first fit, default first bin
							  Double_t high);   // Upper bound for first fit, default last bin
  
  TH1D *makeResidualHisto(TString identifier,
			  TH1D *histo_1,
			  TH1D *histo_2,
			  bool match_bins,
			  bool match_x_range);

  std::pair<double, std::vector<double>> performSideBandAnalysis(TH1D* hist,    //TH1D for sideband analysis
								 double xLow,   //Low bound for band boundary
								 double xHigh); //High bound for band boundary

  TH1D *MirrorHistogram(TH1D *originalHist); //TH1D distribution to mirror

  void FindMatchingFiles(const std::string& directory1,      //path to .csv/.hist
			 const std::string& directory2,      //path to .root
			 const std::string& partialName,     //search string
			 std::vector<std::string>& vector1,  //vector of .csv/.hist file extensions
			 std::vector<std::string>& vector2,  //vector of .root file extensions
			 bool jboyd);                        //bool to alter file name structure to jboyd convention
    
}

#endif
