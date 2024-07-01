#ifndef UTIL_H
#define UTIL_H

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <dirent.h>
#include <filesystem>

#include <string>
#include <regex>
#include <unordered_set>
#include <algorithm>

#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"
#include <TRandom.h>

#include <TSystemDirectory.h>
#include <TList.h>
#include <TCollection.h>
#include <TRegexp.h>

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

  //parse global cut string delimited by && into individual cut strings
  std::vector<std::string> parseCuts(const std::string& cutsString);         //global cut string delimited by &&

  //Find files for pass2
  TString FindFileInDirectories(const TString& pattern,                      //look for runnum
				const std::vector<TString>& directories);    //pass possible directories

  //General
  void readParam(std::string const_path,vector<Double_t> &param); //reads parameters from txt file
  void readParam(std::string const_path,std::string type,Int_t setsize,vector<Double_t> &param); //overload to read subset of data from sbs db-like file

  Double_t generateRandomNumber(const std::string& passkey);      //four-element string which generates a random number between 0.95 and 1.05

  std::vector<std::string> parseGlobalCut(const std::string& input);  //Take in global cut and parse all cut elements for systematic analysis

  bool checkFile(const std::string& filePath);              //Full path to file

  bool checkTH1D(TH1D* histogram,                           //Input histogram to check
		 const std::string& name);                  //Name of histogram to indicate which failed

  //HCal histograms
  TH2D *hhcalrowcol(std::string name); // returns hcal row/col 2d histo
  TH2D *hhcalxy(std::string name); // returns hcal xy 2d histo (data coordinates)
  TH2D *hhcalxy_mc(std::string name); // returns hcal xy 2d histo (mc coordinates)
  TH2D *hdxdy(std::string name); // returns hcal dxdy 2d histo (wide coordinates)
  TH1D *hhsamps(Int_t row,Int_t col,Int_t bins); // returns hcal waveform histogram

  // Check position per event against hcal active area and return minimum N_x(N_y)*blockwidth_x(y) to include cluster center. X= .first, y= .second
  std::pair<double, double> minaa(double hcalx, double hcaly);

  // checks if point is within proton/neutron spot
  bool Nspotcheck(Double_t dy,               //hcaly
		  Double_t dx,               //hcalx
		  Double_t dy_mean,          //elastic peak location in dy
		  Double_t dx_mean,          //elastic peak location in dx
		  Double_t dy_sigma,         //elastic peak sigma in dy
		  Double_t dx_sigma,         //elastic peak sigma in dx
		  Double_t rotationAngle);   //rotation angle in radians, if ever applicable

  // calculates and returns N where N*sig_dx(dy) is necessary to include a dx,dy point in ellipse
  double NspotScaleFactor(Double_t dy,               //hcaly
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

  double CalculateCorrelationFactor(TH2D* hist,              //TH2D to get correlation factor from
				    double xmin,             //xmin for constraint
				    double xmax,             //xmax for constraint
				    double ymin,             //ymin for constraint
				    double ymax);            //ymax for constraint

  void checkTH2DBinsAndEntries(TH2D *hist,                        //TH2D to check bins and entries
			       int minEntries,                    //min_entries
			       int &firstBinWithEntries,          //get first bin with min_entries in x
			       int &lastBinWithEntries,           //get last bin with min_entries in x
			       int &firstBinEntries,              //get first bin with min_entries in x entries
			       int &lastBinEntries);              //get last bin with min_entries in x entries

  void fitAndCalculateRatio(TH1D* hist,                           //dx histogram for LD2
			    TCanvas* canvas,                      //QA histogram
			    double &ratio,                        //yield ratio
			    std::vector<double>& params,          //pass 3(gaus)+3(gaus)+5(pol4)=11 set params
			    std::vector<double>& fitRange,        //pass 2 params (xfitmin, xfitmax)
			    bool fixpol=false);                   //option to fix the pol4 params

  void fitSkewedGaussianAndGetFineParams(TH1D* hist,                  // 1D histo passed for fitting
					 Double_t sig,                   // Estimated sigma for course fit
					 Double_t low,                   // Lower bound for first fit, default first bin
					 Double_t high,                  // Upper bound for first fit, default last bin
					 std::vector<Double_t> &params,  //vector to store parameters
					 Double_t alpha_ll);             // Lower Limit on alpha (optional)
  
  TH1D *makeResidualHisto(TString identifier,                     //string to identify the variable analyzed
			  TH1D *histo_1,                          //histogram 1 (base)
			  TH1D *histo_2,                          //histogram 2 (subtracted)
			  bool match_bins,                        //bool to error on N bins mismatch
			  bool match_x_range);                    //bool to error on x range mismatch
  
  TH1D *makeResidualHistoWithError(TString identifier,                     //string to identify the variable analyzed
				   TH1D *histo_1,                          //histogram 1 (base)
				   TH1D *histo_2,                          //histogram 2 (subtracted)
				   bool match_bins,                        //bool to error on N bins mismatch
				   bool match_x_range);                    //bool to error on x range mismatch

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
    
  /* void SyncFilesWithCsv(const std::string& directory1,                  //path to .csv/.hist */
  /* 			const std::string& directory2,                  //path to .root */
  /* 			const std::string& partialName,                 //search string */
  /* 			std::vector<std::string>& vector1,              //vector of .root file extensions */
  /* 			std::map<int, std::pair<std::string, std::vector<float>>>& csvData);   //map jobid->csv data on root file and rootfile extension */

    
  void SyncFilesWithCsv(const std::string& directory1,                  //path to .csv/.hist
			const std::string& directory2,                  //path to .root
			const std::string& partialName,                 //search string
			std::vector<std::string>& vector1,              //vector of .root file extensions
			std::pair<std::string, std::vector<float>>& csvData);   //rootfile path and metadata vector

  void synchronizeJobNumbers(std::vector<std::string>& vecA,   //Vector of strings A of form "bla_job<int>bla"
			     std::vector<std::string>& vecB);  //Vector of strings B of form "bla_job<int>bla"

  //overload for pairs
  void synchronizeJobNumbers(std::vector<std::pair<std::string, std::vector<float>>>& vecA,
			     std::vector<std::pair<std::string, std::vector<float>>>& vecB);

  //overload for maps
  void synchronizeJobNumbers(std::map<int, std::pair<std::string, std::vector<float>>>& csvData_n,
			     std::map<int, std::pair<std::string, std::vector<float>>>& csvData_p);

  std::vector<TH1D*> createBootstrapSamples(TH1D* originalHist,             //Original dataset (full)
					    int nBootstrapSamples);         //Number of sample datasets desired

  TH1D* shiftHistogramX(TH1D* originalHist,                  //original dataset (full)
			double shiftValue);                  //amount to shift each bin
  
  TH1D* cloneAndCutHistogram(TH1D* originalHist,             //original dataset (full) 
			     double xMin,                    //constraint min
			     double xMax);                   //constraint max

  double GetTotalQuadratureError(TH1D* hist);                //input histogram 

  TH1D* subtractFunctionFromHistogramWithinRange(const TH1D* hist, //input histogram
						 const TF1* func,  //input tuned functional fit to bg
						 double xmin,      //<xmin, zero
						 double xmax);     //>xmax, zero

  void subtractFunctionAndGetTotalAndError(TH1D* hist,         //data histogram
					   TF1* func,          //fit to background
					   double xMin,        //reject point region else full histo, x minimum
					   double xMax,        //reject point region else full histo, x maximum
					   double &total,      //total difference, modifying input
					   double xRangeLow,   //full histo region, x minimum
					   double xRangeHigh,  //full histo region, x maximum
					   double &error,      //total error
					   bool abs_sub);      //bool to choose if hist-fit diff >0 else 0

  void plotNormalizedTH1DsOnCanvas(TCanvas* canvas,                              //passed canvas for viewing
				   TH1D* hist1,                                  //first TH1D
				   TH1D* hist2,                                  //second TH1D
				   const char* title = "Normalized Histograms",  //Title for overlay
				   const char* x_label = "X-axis",               //xaxis title for overlay
				   const char* y_label = "Frequency");           //yaxis title for overlay

  
  void parseAndDisplayCuts(const char* name, const char* cuts);
  void parseAndDisplayCuts(const char* name, const char* cuts, TCanvas *c1);     //overload to modify c1

}

#endif
