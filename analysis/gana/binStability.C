//sseeds 6.14.24: Script to vary the number of bins in dx, extract Rsf, and evaluate bin stability for fits

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TKey.h>
#include <TClass.h>
#include <iostream>
#include <vector>
#include <string>
#include <TLatex.h>
#include <algorithm>
#include <TRandom.h>
#include "../../include/gmn.h"
#include "../../src/jsonmgr.C"

//Fit range
double hcalfit_l; //default defined by json
double hcalfit_h; //default defined by json

double HDE_syst_N = 94.316; //percent error on nucleon yield (from proton HDE analysis, assuming similar)
double HDE_systerr_N = 0.583; //percent error on neutron yield (from proton HDE analysis, assuming similar)

//Fit options
std::string fitopt = "RMQ0";

//Exclude these globalcuts until MC grinch ToT is ready
std::vector<std::string> excludeCuts = {"bb_grinch_tdc_clus_size","bb_grinch_tdc_clus_trackindex"};

//Total fits using Interpolate with elastic signal histo and 4th order poly fit to bg
TH1D *hdx_p;
TH1D *hdx_n;
TH1D *hdx_p_wide;
TH1D *hdx_n_wide;

Double_t fit_gausBG(double *x, double *par) {
  // MC float params for scaling and shifting
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p_wide->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n_wide->Interpolate(x[0] - dx_shift_n);
  
  // Sum of components
  return proton + neutron + fits::g_gfit(x, &par[4]);
}

Double_t fit_pol2BG(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron + fits::g_p2fit_cd(x, &par[4]);
}

// Forward declarations
void handleError(TFile *file1, TFile *file2, std::string marker);
void handleError(TFile *file1, std::string marker);
void createAndSaveCanvasWithFitsAndResiduals(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, TF1* bg, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, int pol, double &elastics, std::pair<double,double> &rsf, int kine, int mag);
std::vector<std::pair<double, double>> fitAndFineFit(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, const std::string& fitOptions = "RBMQ0");
std::vector<std::pair<double, double>> fitAndFineFit_fixshift(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, double pshift, double nshift, const std::string& fitOptions = "RBMQ0");
TH1D* calculateResiduals(TH1D* hData, TF1* fit, const char* residualName);
std::vector<std::string> split(const std::string &s, char delimiter);
std::map<std::string, std::string> getRowContents(const std::string &filePath, int kine, int mag, const std::string &target, const std::vector<std::string> &excludeKeys);
std::string addbc(std::string input);

/////////////////
/////////////////
/////////////////
/////////////////
/////////////////
/////////////////
/////////////////

//MAIN. kine=kinematic, mag=fieldsetting, pass=pass#, binsets=total plots, binranges=bin step size, bestclus=use best cluster plots, thin=use plot file without correlations, widecut=use plots with wide cuts, effz=use plots with effective z applied, alt=use plot file using MC alt files
void binStability(int kine=4, 
		  int mag=30, 
		  int pass=2, 
		  int binsets=8,
		  int binranges=50,
		  bool bestclus=true, 
		  bool thin=true,
		  bool widecut=false,
		  bool effz=true,
		  bool alt = true) {
  //set draw params
  gStyle->SetPalette(55);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetLabelSize(0.03, "XYZ");
  gStyle->SetTitleSize(0.04, "XYZ");
  //gStyle->SetEndErrorSize(0);

  // Set up configuration and tune objects to load analysis parameters
  SBSconfig config(kine,0);

  // Obtain configuration pars from config class
  double E_beam = config.GetEbeam();
  double BB_angle = config.GetBBtheta_rad(); // In radians

  //load json file
  JSONManager *jmgr = new JSONManager("../../config/syst.json");

  //Get dx plot details
  hcalfit_l = jmgr->GetValueFromSubKey<double>( "hcalfit_l", Form("sbs%d_%d",kine,mag) );
  hcalfit_h = jmgr->GetValueFromSubKey<double>( "hcalfit_h", Form("sbs%d_%d",kine,mag) );

  cout << "Loaded hcal plot x limits: " << hcalfit_l << " to " << hcalfit_h << endl;

  int def_bins = jmgr->GetValueFromSubKey<int>( "hbins", Form("sbs%d_%d",kine,mag) );
  
  // Get tight elastic cuts
  std::string globalcutsbla = jmgr->GetValueFromSubKey_str(Form("post_tcuts_p%d", pass), Form("sbs%d_%d", kine, mag));
  cout << "Loaded OLD cuts: " << globalcutsbla << endl;

  //get new cuts from .csv
  std::string cutsheet_path = "/w/halla-scshelf2102/sbs/seeds/ana/data/p2_cutset.csv";
  std::string target = "ld2";

  // Get the row contents
  std::map<std::string, std::string> rowContents = getRowContents(cutsheet_path, kine, mag, target, excludeCuts);
  // Print the row contents
  cout << endl << "Loading NEW cuts..." << endl;
  for (const auto &content : rowContents) {
    std::cout << content.first << ": " << content.second << std::endl;
  }
  cout << endl;

  std::string globalcuts_raw;
  for (const auto &content : rowContents) {
    if (!content.second.empty()) {
      if (!globalcuts_raw.empty()) {
	globalcuts_raw += "&&";
      }
      globalcuts_raw += content.second;
    }
  }

  cout << endl <<"Concatenated globalcuts_raw: " << globalcuts_raw << endl << endl;

  std::string globalcuts = addbc(globalcuts_raw);

  cout << endl <<"Concatenated globalcuts: " << globalcuts << endl << endl;

  std::vector<std::string> cuts = util::parseCuts(globalcuts); //this makes a vector of all individual cuts in the single globalcut string.

  std::cout << "Parsed cuts: " << std::endl;
  for ( size_t i=0; i<cuts.size(); ++i )
    cout << cuts[i] << endl;

  std::string bestclus_word = "";
  if(bestclus)
    bestclus_word = "_bc";

  std::string thin_word = "";
  if(thin)
    thin_word = "_thin";

  std::string alt_word = "";
  if(alt)
    alt_word = "_alt";

  std::string wide_word = "";
  if(widecut)
    wide_word = "_widecut";

  std::string effz_word = "";
  if(effz)
    effz_word = "_effz";


  //Set up paths
  std::string basePath = "/lustre19/expphy/volatile/halla/sbs/seeds";

  std::string fmcinPath = Form("%s/gmn_analysis/dx_mc_sbs%d_mag%d_pass%d%s%s%s%s.root", basePath.c_str(), kine, mag, pass, thin_word.c_str(), alt_word.c_str(),wide_word.c_str(), effz_word.c_str());

  std::string foutPath = Form("%s/gmn_analysis/binstab_sbs%d_mag%d_pass%d%s%s.root", basePath.c_str(), kine, mag, pass, bestclus_word.c_str(), effz_word.c_str());

  std::string fin_path = basePath + Form("/parse/parse_sbs%d_pass%d_barebones%s_ld2.root",kine,pass,effz_word.c_str());

  cout << endl << "    Drawing data histograms from parsed file: " << fin_path << endl << endl << endl;

  // Open the parse file
  TFile* inputFile = new TFile(fin_path.c_str());
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Error opening file: " << fin_path << std::endl;
    return;
  }

  // Get the tree
  TTree* tree = dynamic_cast<TTree*>(inputFile->Get("P"));
  if (!tree) {
    std::cerr << "Tree not found in file: " << fin_path << std::endl;
    inputFile->Close();
    return;
  }

  //Create binsets histograms
  std::vector<TH1D*> histograms;
  std::vector<int> histobins;

  int minbins = 100;

  for (int i = 1; i <= binsets; ++i) {
    int numBins = i * binranges + minbins;
    histobins.push_back(numBins);
    std::string histName = Form("hdx_allcut_%d_bins", numBins);
    TH1D* hist = new TH1D(histName.c_str(), Form("dx (best cluster) with %d bins;m", numBins), numBins, hcalfit_l, hcalfit_h);
    tree->Draw(("dx_bc>>" + histName).c_str(), globalcuts.c_str(), "COLZ");
    hist->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
    histograms.push_back(hist);
  }
  
  // Create one more histogram with default bins
  TH1D* controlHist = new TH1D("hdx_allcut_def_bins", "dx (best cluster) with default bins;m", 400, hcalfit_l, hcalfit_h);
  tree->Draw("dx_bc>>hdx_allcut_def_bins", globalcuts.c_str(), "COLZ");
  controlHist->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  histograms.push_back(controlHist);
  histobins.push_back(400);
  
  //get histograms from MC plot file
  TFile* inputFileMC = new TFile(fmcinPath.c_str(), "READ");
  if (!inputFileMC || inputFileMC->IsZombie())
    handleError(inputFile,inputFileMC,"inputFileMC");

  //fix the interpolate functions
  TH1D *hdx_p_raw = dynamic_cast<TH1D*>(inputFileMC->Get("hdx_p"));
  hdx_p = (TH1D*)(hdx_p_raw->Clone("hdx_p"));
  TH1D *hdx_p_clone = (TH1D*)(hdx_p_raw->Clone("hdx_p_clone"));
  hdx_p->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  hdx_p_wide = (TH1D*)(hdx_p_raw->Clone("hdx_p"));
  hdx_p_wide->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  TH1D *hdx_n_raw = dynamic_cast<TH1D*>(inputFileMC->Get("hdx_n"));
  hdx_n = (TH1D*)(hdx_n_raw->Clone("hdx_n"));
  TH1D *hdx_n_clone = (TH1D*)(hdx_n_raw->Clone("hdx_n_clone"));
  hdx_n->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  hdx_n_wide = (TH1D*)(hdx_n_raw->Clone("hdx_n"));
  hdx_n_wide->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  //Handle missing histograms
  if (!hdx_p_raw || !hdx_n_raw ) 
    handleError(inputFile,inputFileMC,"hdx");

  TFile* outputFile = new TFile(foutPath.c_str(), "RECREATE");

  //For now, best fit to data. Use to get shifts for all other fits.
  std::pair<double,double> qualset;
  auto setpar_vector = fitAndFineFit(controlHist, "pol2BG_set", "fit_pol2BG", 7, hcalfit_l, hcalfit_h, qualset, fitopt.c_str());

  cout << "Parameter set vector: ";
  for (auto par : setpar_vector )
    cout << "(" << par.first << ", " << par.second << ") ";
  cout << endl << endl;

  //Get the same shift for all BG comparisons
  double pshift = setpar_vector[2].first;
  double nshift = setpar_vector[3].first;

  // Make background function for control
  TF1 *bg_p2Fit = new TF1("bg_pol2Fit",fits::g_p2fit_cd,hcalfit_l,hcalfit_h,3);
  for (int i=0; i<3; ++i){
    bg_p2Fit->SetParameter(i,setpar_vector[i+4].first);
    bg_p2Fit->SetParError(i,setpar_vector[i+4].second);
  }

  //print the control to canvas
  double control_p2;
  std::pair<double,double> control_rsf_p2;
  createAndSaveCanvasWithFitsAndResiduals(controlHist, hdx_p_clone, hdx_n_clone, bg_p2Fit, setpar_vector, qualset, "pol2 BG", "~/gmn_plots/control_p2.pdf", 2, control_p2, control_rsf_p2, kine, mag);

  // Fit the histograms
  std::vector<std::pair<double, double>> qualp2_vec;
  std::vector<std::vector<std::pair<double,double>>> p2par_vector_vec;
  std::vector<double> Rsf_values, Rsf_errors;

  for (size_t i = 0; i < histograms.size(); ++i){
    std::pair<double, double> qualp2;
    auto p2par_vector = fitAndFineFit_fixshift(histograms[i], "pol2BG", "fit_pol2BG", 7, hcalfit_l, hcalfit_h, qualp2, pshift, nshift, fitopt.c_str());
    cout << "Pol 2 vector for histogram " << histograms[i]->GetName() << ": ";
    for (auto par : p2par_vector)
      cout << "(" << par.first << ", " << par.second << ") ";
    cout << endl << endl;

    p2par_vector_vec.push_back(p2par_vector);
    qualp2_vec.push_back(qualp2);

    // Calculate Rsf and its error for each histogram
    double psum = p2par_vector[0].first;
    double psum_error = p2par_vector[0].second;
    double nsum = p2par_vector[1].first;
    double nsum_error = p2par_vector[1].second;
    double Rsf = nsum / psum;
    double Rsf_error = Rsf * sqrt(pow(psum_error / psum, 2) + pow(nsum_error / nsum, 2));
    
    Rsf_values.push_back(Rsf);
    Rsf_errors.push_back(Rsf_error);
  }

  // Create TGraphs for chi-square/ndf vs. nbins and Rsf vs. nbins. Color thing does nothing.
  TGraphErrors* grChi2NDF = new TGraphErrors();
  TGraphErrors* grRsf = new TGraphErrors();
  for (size_t i = 0; i < histograms.size() - 1; ++i) {
    int nbins = histobins[i];
    double chi2ndf = qualp2_vec[i].first / qualp2_vec[i].second;
    
    // Set the point for grChi2NDF
    grChi2NDF->SetPoint(i, nbins, chi2ndf);
    grChi2NDF->SetPointError(i, 0, 0); // No error on nbins or chi2/ndf
    
    // Set the color of the point if nbins == 400
    if (nbins == 400) {
      grChi2NDF->SetMarkerColor(kGreen);
      grChi2NDF->SetMarkerStyle(21);
    } else {
      grChi2NDF->SetMarkerColor(kBlack); // Default color
      grChi2NDF->SetMarkerStyle(20);
    }
    
    double Rsf = Rsf_values[i];
    double Rsf_error = Rsf_errors[i];
    
    // Set the point for grRsf
    grRsf->SetPoint(i, nbins, Rsf);
    grRsf->SetPointError(i, 0, Rsf_error);
    
    // Set the color of the point if nbins == 400
    if (nbins == 400) {
      grRsf->SetMarkerColor(kGreen);
      grRsf->SetMarkerStyle(21);
    } else {
      grRsf->SetMarkerColor(kBlack); // Default color
      grRsf->SetMarkerStyle(20);
    }
  }

  // Draw the control histogram
  //TCanvas* cControl = new TCanvas("cControl", "Control Histogram", 800, 600);
  //controlHist->Draw();
  //controlHist->Fit("fit_pol2BG");

  // Determine number of rows and columns dynamically
  int nPads = binsets + 1; // Include control histogram
  int nCols = std::ceil(std::sqrt(nPads));
  int nRows = std::ceil(static_cast<double>(nPads) / nCols);

  // Draw the other histograms
  TCanvas* c1 = new TCanvas("c1", "Bin Stability Histograms", 1800, 1200);
  // int nPads = binsets + 1; // Include control histogram
  // int nCols = 3;
  // int nRows = (nPads + nCols - 1) / nCols;
  c1->Divide(nCols, nRows);

  for (size_t i = 0; i < histograms.size(); ++i) {
    c1->cd(i + 1);

    //get the fit pars and chisqr/ndf for this histogram
    std::vector<std::pair<double,double>> pars = p2par_vector_vec[i];
    std::pair<double,double> qual = qualp2_vec[i];

    TF1* fit = new TF1(TString::Format("fitp2_binset%zu", i), fit_pol2BG, hcalfit_l, hcalfit_h, 7);

    for (size_t i = 0; i < pars.size(); ++i) {
      fit->SetParameter(i, pars[i].first);
      fit->SetParError(i, pars[i].second);
    }

    fit->SetChisquare(qual.first);
    fit->SetNDF(qual.second);

    gPad->SetGridx();
    gPad->SetGridy();

    std::string controlword = "";
    int lasthist = histograms.size()-1;
    if(i==lasthist)
      controlword = " (Shift Parameter Control)";
    
    histograms[i]->SetTitle(Form("dx (SBS-%d, %d%% field)%s;m", kine, mag, controlword.c_str()));    histograms[i]->SetTitleFont(132);
    histograms[i]->GetXaxis()->SetTitleFont(132);
    histograms[i]->GetXaxis()->SetLabelFont(132);
    histograms[i]->GetYaxis()->SetTitleFont(132);
    histograms[i]->GetYaxis()->SetLabelFont(132);
    histograms[i]->SetMarkerStyle(24);
    histograms[i]->SetMarkerSize(0.5);
    histograms[i]->Draw("P");

    fit->SetLineColor(kGreen);
    fit->SetLineWidth(1);
    fit->SetFillColorAlpha(kGreen,0.35);
    fit->SetFillStyle(1001);
    fit->SetNpx(1000);
    fit->Draw("sameFC");

    TLegend* leg = new TLegend(0.5, 0.5, 0.89, 0.89);
    leg->SetBorderSize(0); // Remove border
    leg->SetFillStyle(0);
    leg->SetTextColor(kRed);
    leg->SetTextFont(132);
    leg->SetTextSize(0.05);
    leg->AddEntry(histograms[i], Form("Data, %d bins", histobins[i]), "p");
    leg->AddEntry(fit, "Fit (MC + BG)", "f");
    leg->AddEntry( (TObject*)0, Form("R_{sf} : %0.4f",Rsf_values[i]), "");
    leg->AddEntry( (TObject*)0, Form("#chi^{2}/ndf: %0.0f/%d",fit->GetChisquare(),fit->GetNDF()), "");
    leg->Draw("same");

  }

  // Create a canvas for the graphs
  TCanvas* cGraphs = new TCanvas("cGraphs", "Chi-Square/NDF and Rsf vs. Fit Range", 1800, 800);
  cGraphs->Divide(2, 1);

  // Draw chi-square/ndf vs. fit range
  cGraphs->cd(1);
  gPad->SetGridx();
  gPad->SetGridy();
  grChi2NDF->SetTitle("#chi^{2}/NDF vs. N bins;N bins;#chi^{2}/NDF");
  grChi2NDF->SetMarkerStyle(21);
  grChi2NDF->GetYaxis()->SetRangeUser(0.0,4.0);
  grChi2NDF->Draw("AP");

  // Draw Rsf vs. fit range
  cGraphs->cd(2);
  gPad->SetGridx();
  gPad->SetGridy();
  grRsf->SetTitle("R_{sf} vs. N bins;N bins;R_{sf}");
  grRsf->SetMarkerStyle(21);
  grRsf->GetYaxis()->SetRangeUser(0.8,1.1);
  grRsf->Draw("AP");
  
  // Define the maximum number of cuts per page and the margin
  const int maxCutsPerPage = 20; // Adjust this as needed
  const double margin = 0.1;
  const double textHeight = 0.04;

  // Create a canvas with height adjusted for the number of cuts
  int numCuts = cuts.size();
  int canvasHeight = 800 + (numCuts - maxCutsPerPage) * 40;
  if (canvasHeight < 800) canvasHeight = 800; // Ensure minimum height

  TCanvas *cCuts = new TCanvas("cCuts", "Itemized Cuts", 1200, canvasHeight);
  cCuts->cd();
  TLatex latex;
  latex.SetTextSize(textHeight);

  // Draw the title
  latex.DrawLatexNDC(0.1, 0.9, "Itemized Cuts:");

  // Adjust the text size if necessary
  double adjustedTextSize = textHeight;
  if (numCuts > maxCutsPerPage) {
    adjustedTextSize = (0.9 - margin) / numCuts;
    latex.SetTextSize(adjustedTextSize);
  }

  // Draw the cuts
  for (size_t j = 0; j < cuts.size(); ++j) {
    latex.DrawLatexNDC(0.1, 0.9 - adjustedTextSize * (j + 1), cuts[j].c_str());
  }

  cCuts->Update();
  cCuts->Write();

  outputFile->Write();

  cout << "All plots created. Output file located here: " << foutPath << endl;

}

/////////////
/////////////
//functions
/////////////
/////////////
/////////////


void handleError(TFile *file1, TFile *file2, std::string marker) {
    if (file1) file1->Close();
    if (file2) file2->Close();
    std::cerr << "Error: File opening or histogram retrieval failed at " << marker << "." << std::endl;
}

void handleError(TFile *file1, std::string marker) {
    if (file1) file1->Close();
    std::cerr << "Error: File opening or histogram retrieval failed at " << marker << "." << std::endl;
}

void createAndSaveCanvasWithFitsAndResiduals(TH1D* hdx, TH1D* hdxp, TH1D* hdxn, TF1* bg, const std::vector<std::pair<double,double>> pars, std::pair<double,double> qual, const char* type, const char* savePath, int pol, double &elastics, std::pair<double,double> &rsf, int kine, int mag) {

  // Generate a unique identifier for the fit function name
  TString uniqueID = TString::Format("_%u", gRandom->Integer(1000000));

  //Clone the dx histo for bg sub
  TH1D *hdx_histplot = (TH1D*)(hdx->Clone("hdx_histplot"));
  TH1D *hdx_clone = (TH1D*)(hdx->Clone("hdx_clone"));
  TH1D *hdx_p_clone = (TH1D*)(hdxp->Clone("hdx_p_clone"));
  TH1D *hdx_n_clone = (TH1D*)(hdxn->Clone("hdx_n_clone"));

  //Recreate the fit
  TF1* fit;
  std::string bgword = "";
  if (pol == 2){
    fit = new TF1(TString::Format("fitp2%s", uniqueID.Data()), fit_pol2BG, hcalfit_l, hcalfit_h, 7);
    bgword = "Order-2 Polynomial ";
  }else{
    cout << "ERROR: bin stability supports only pol2 bg fits." << endl;
    return;
  }

  for (size_t i = 0; i < pars.size(); ++i) {
    fit->SetParameter(i, pars[i].first);
    fit->SetParError(i, pars[i].second);
  }

  //std::pair<double,double> altchisqrndf = calculateChiSquareWithZeroMCError(hdx, fit);
  
  //fit->SetChisquare(altchisqrndf.first);
  //fit->SetNDF(altchisqrndf.second);

  fit->SetChisquare(qual.first);
  fit->SetNDF(qual.second);

  // Create a canvas with two pads
  TCanvas *cTotal = new TCanvas(savePath, Form("Data with Fits and Residuals %s",type), 1200, 800);
  TPad *pad1 = new TPad("pad1", "Pad with the fit", 0.0, 0.3, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "Pad with the residuals", 0.0, 0.0, 1.0, 0.3);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad1->Draw();
  pad2->Draw();

  // First pad: Draw the histogram with the fits
  pad1->cd();

  gPad->SetGridx();
  gPad->SetGridy();

  //hdx_histplot->SetLineColor(kBlack);
  //hdx_histplot->SetLineWidth(2);
  hdx_histplot->SetTitle(Form("dx, %s (SBS-%d, %d%% field);m", type, kine, mag));
  hdx_histplot->SetTitleFont(132);
  hdx_histplot->GetXaxis()->SetTitleFont(132);
  hdx_histplot->GetXaxis()->SetLabelFont(132);
  hdx_histplot->GetYaxis()->SetTitleFont(132);
  hdx_histplot->GetYaxis()->SetLabelFont(132);
  hdx_histplot->SetMarkerStyle(21);
  hdx_histplot->SetMarkerSize(0.5);
  hdx_histplot->Draw("P");
  hdx->SetLineColor(kBlack);
  hdx->SetLineWidth(1);
  hdx->Draw("E same");
  fit->SetLineColor(kGreen);
  fit->SetLineWidth(1);
  fit->SetFillColorAlpha(kGreen,0.35);
  fit->SetFillStyle(1001);
  fit->SetNpx(1000);
  fit->Draw("sameFC");

  hdx_p_clone = util::shiftHistogramX(hdx_p, pars[2].first);
  hdx_p_clone->Scale(pars[0].first);
  hdx_p_clone->SetLineColor(kRed);
  hdx_p_clone->SetLineWidth(2);
  hdx_p_clone->SetMarkerStyle(20);  // Filled circle marker
  hdx_p_clone->SetMarkerColor(kRed);
  hdx_p_clone->SetMarkerSize(0.5);
  hdx_p_clone->Draw("P same");

  hdx_n_clone = util::shiftHistogramX(hdx_n, pars[3].first);
  hdx_n_clone->Scale(pars[1].first);
  hdx_n_clone->SetLineColor(kBlue);
  hdx_n_clone->SetLineWidth(2);
  hdx_n_clone->SetMarkerStyle(20);  // Filled circle marker
  hdx_n_clone->SetMarkerColor(kBlue);
  hdx_n_clone->SetMarkerSize(0.5);
  hdx_n_clone->Draw("P same");

  // Set background function
  bg->SetLineColor(kBlack);
  bg->SetFillColorAlpha(kBlack,0.35);
  bg->SetFillStyle(3013);
  bg->Draw("same");

  double psumerror;
  double psum = hdx_p_clone->IntegralAndError(0, hdx_p_clone->GetNbinsX() + 1, psumerror, "");

  double pscaleerror = pars[0].second;
  double pscale = pars[0].first;

  double nsumerror;
  double nsum = hdx_n_clone->IntegralAndError(0, hdx_n_clone->GetNbinsX() + 1, nsumerror, "");

  double nscaleerror = pars[1].second;
  double nscale = pars[1].first;

  double np_sum_ratio = nsum/psum;
  double np_par_ratio = nscale/pscale;
  double np_sum_ratio_error = np_sum_ratio * sqrt(pow(psumerror/psum, 2) + pow(nsumerror/nsum, 2));
  double np_par_ratio_error = np_par_ratio * sqrt(pow(pscaleerror/pscale, 2) + pow(nscaleerror/nscale, 2));

  rsf.first = np_par_ratio;
  rsf.second = np_par_ratio_error;

  // Construct legend within function
  TLegend* leg = new TLegend(0.5, 0.5, 0.89, 0.89);
  leg->SetBorderSize(0); // Remove border
  //leg->SetFillStyle(0);  // Make background transparent
  leg->SetTextColor(kBlack);
  leg->SetTextFont(132);
  leg->SetTextSize(0.05);
  leg->AddEntry(hdx_histplot, "Data", "p");
  leg->AddEntry(fit, "Fit (Proton + Neutron + BG)", "f");
  leg->AddEntry(bg, Form("%s Background",bgword.c_str()), "f");
  leg->AddEntry(hdx_p_clone, "Proton SIMC MC", "p");
  leg->AddEntry(hdx_n_clone, "Neutron SIMC MC", "p");
  leg->AddEntry( (TObject*)0, "", "");
  leg->AddEntry( (TObject*)0, Form("n/p scale ratio R_{sf} : %0.3f #pm %0.3f",np_par_ratio, np_par_ratio_error), "");
  leg->AddEntry( (TObject*)0, "", "");
  //leg->AddEntry( (TObject*)0, Form("data N events : %0.0f",hdx->GetEntries()), "");
  //leg->AddEntry( (TObject*)0, Form("n/p yield ratio R'' : %0.3f #pm %0.3f",np_sum_ratio, np_sum_ratio_error), "");
  //leg->AddEntry( (TObject*)0, Form("dx shift pars, n/p : %0.3f / %0.3f ", pars[3].first, pars[2].first), "");
  leg->AddEntry( (TObject*)0, Form("#chi^{2}/ndf: %0.3f/%d",fit->GetChisquare(),fit->GetNDF()), "");
  leg->Draw("same");

  //Background subraction on dx_clone and total elastic statistics
  elastics = 0;
  for (int bin = 1; bin <= hdx->GetNbinsX(); ++bin) {
    double bgValue = bg->Eval(hdx_clone->GetXaxis()->GetBinCenter(bin));
    if( hdx_clone->GetBinContent(bin) > bgValue )
      elastics += hdx_clone->GetBinContent(bin) - bgValue;
    hdx_clone->SetBinContent(bin, hdx_clone->GetBinContent(bin) - bgValue);
    hdx_clone->SetBinError(bin, hdx_clone->GetBinError(bin));
  }

  TH1D *hRes = calculateResiduals(hdx, fit, Form("hResiduals_%s",uniqueID.Data()));
  TH1D *hRes_histplot = (TH1D*)(hRes->Clone(Form("hRes_histplot_%s",type)));

  hRes_histplot->SetTitle("");
  hRes_histplot->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  double dxMaxValue = hdx_clone->GetMaximum();
  hRes_histplot->GetYaxis()->SetRangeUser(-dxMaxValue/4,dxMaxValue/4);

  hRes_histplot->SetLineColor(kBlack);
  hRes_histplot->SetLineWidth(1);
  //hRes_histplot->SetFillStyle(1001);
  //hRes_histplot->SetFillColorAlpha(kRed,0.75);

  // Second pad: Draw the residuals
  pad2->cd();
  gPad->SetGridx();

  hRes_histplot->GetYaxis()->SetTitle("Residuals");
  hRes_histplot->GetYaxis()->SetTitleSize(0.1);
  hRes_histplot->GetYaxis()->SetTitleOffset(0.5);
  hRes_histplot->GetYaxis()->SetLabelSize(0.1);
  hRes_histplot->GetXaxis()->SetTitleSize(0.1);
  hRes_histplot->GetXaxis()->SetLabelSize(0.1);
  hRes_histplot->GetXaxis()->SetTitle("x_{hcal}-x_{exp}");
  hRes_histplot->GetXaxis()->SetTitleFont(132);
  hRes_histplot->GetXaxis()->SetLabelFont(132);
  hRes_histplot->GetYaxis()->SetTitleFont(132);
  hRes_histplot->GetYaxis()->SetLabelFont(132);
  hRes_histplot->Draw("E");

  // Get the axis range
  double y_min = hRes_histplot->GetMinimum();
  double y_max = hRes_histplot->GetMaximum();

  // Create a TPave to fill the region between the axes
  TPave *pave = new TPave(hcalfit_l, y_min, hcalfit_h, y_max, 0, "NB");
  pave->SetFillColorAlpha(kGray,0.4); // Light grey RGB
  pave->SetFillStyle(3001); // 3001: transparent fill style
  pave->Draw("same");

  // Draw a line at y=0 for residuals
  TF1 *zeroLine = new TF1("zeroLine", "0", hRes->GetXaxis()->GetXmin(), hRes->GetXaxis()->GetXmax());
  zeroLine->SetLineColor(kRed);
  zeroLine->SetLineStyle(2);
  zeroLine->Draw("same");

  // Redraw the histogram axes and title over the TPave
  hRes_histplot->Draw("AXIS SAME");

  // Update the canvas
  cTotal->cd();
  cTotal->Update();

  // Save the canvas to a file
  cTotal->SaveAs(savePath);
}

//Function to fit then fine fit and return all fit parameters
std::vector<std::pair<double, double>> fitAndFineFit(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, const std::string& fitOptions = "RBMQ0") {
  TF1* fit = new TF1(fitName.c_str(), fitFormula.c_str(), hcalfit_l, hcalfit_h, paramCount);
  //fit->SetNpx(5000);
  for (int i=0; i<paramCount; ++i){ //reset parameters/errors for this set
    fit->SetParameter(i,0);
    fit->SetParError(i,0);
  }

  histogram->Fit(fit, fitOptions.c_str());

  std::vector<std::pair<double, double>> parametersAndErrors(paramCount);
  for (int i = 0; i < paramCount; ++i) {
    parametersAndErrors[i].first = fit->GetParameter(i); // Parameter value
    parametersAndErrors[i].second = fit->GetParError(i); // Parameter error
  }

  // Fine Fit
  TF1* fineFit = new TF1((fitName + "_fine").c_str(), fitFormula.c_str(), hcalfit_l, hcalfit_h, paramCount);
  std::vector<double> fineFitInitialParams(paramCount);
  for (int i = 0; i < paramCount; ++i) {
    fineFitInitialParams[i] = parametersAndErrors[i].first;
  }
  fineFit->SetParameters(fineFitInitialParams.data());
  histogram->Fit(fineFit, fitOptions.c_str());

  // Update parameters and errors with fine fit results
  for (int i = 0; i < paramCount; ++i) {
    parametersAndErrors[i].first = fineFit->GetParameter(i); // Fine fit parameter value
    parametersAndErrors[i].second = fineFit->GetParError(i); // Fine fit parameter error
  }

  fitqual.first = fineFit->GetChisquare();
  fitqual.second = fineFit->GetNDF();

  delete fit; // Delete fit to avoid carry-over
  delete fineFit; // Clean up
  return parametersAndErrors;
}

//Function to fit then fine fit and return all fit parameters
std::vector<std::pair<double, double>> fitAndFineFit_fixshift(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, double pshift, double nshift, const std::string& fitOptions = "RBMQ0") {
  TF1* fit = new TF1(fitName.c_str(), fitFormula.c_str(), hcalfit_l, hcalfit_h, paramCount);
  for (int i=0; i<paramCount; ++i){ //reset parameters/errors for this set
    fit->SetParameter(i,0);
    fit->SetParError(i,0);
  }

  double histfitmin = histogram->GetMaximum()/100;
  double histfitmax = histogram->GetMaximum();
  double histfitrmin = hcalfit_l;
  double histfitrmax = hcalfit_h;
  double histfitrange = histfitrmax-histfitrmin;

  fit->FixParameter(2,pshift); //fix x shift for proton
  fit->FixParameter(3,nshift); //fix x shift for neutron

  histogram->Fit(fit, fitOptions.c_str());

  std::vector<std::pair<double, double>> parametersAndErrors(paramCount);
  for (int i = 0; i < paramCount; ++i) {
    parametersAndErrors[i].first = fit->GetParameter(i); // Parameter value
    parametersAndErrors[i].second = fit->GetParError(i); // Parameter error
  }

  // Fine Fit
  TF1* fineFit = new TF1((fitName + "_fine").c_str(), fitFormula.c_str(), hcalfit_l, hcalfit_h, paramCount);
  std::vector<double> fineFitInitialParams(paramCount);
  for (int i = 0; i < paramCount; ++i) {
    fineFitInitialParams[i] = parametersAndErrors[i].first;
  }

  fineFit->FixParameter(2,pshift); //keep fixed x shift
  fineFit->FixParameter(3,nshift); //keep fixed x shift

  fineFit->SetParameters(fineFitInitialParams.data());
  histogram->Fit(fineFit, fitOptions.c_str());

  // Update parameters and errors with fine fit results
  for (int i = 0; i < paramCount; ++i) {
    parametersAndErrors[i].first = fineFit->GetParameter(i); // Fine fit parameter value
    parametersAndErrors[i].second = fineFit->GetParError(i); // Fine fit parameter error
  }

  fitqual.first = fineFit->GetChisquare();
  fitqual.second = fineFit->GetNDF();

  delete fit; // Delete fit to avoid carry-over
  delete fineFit; // Clean up
  return parametersAndErrors;
}

TH1D* calculateResiduals(TH1D* hData, TF1* fit, const char* residualName) {
  // Clone the data histogram to create the residuals histogram
  TH1D* hResiduals = (TH1D*)(hData->Clone(residualName));
  hResiduals->Reset();  // Clear the histogram

  // Loop over each bin and calculate the residuals
  for (int bin = 1; bin <= hData->GetNbinsX(); ++bin) {
    double x = hData->GetXaxis()->GetBinCenter(bin);
    double dataValue = hData->GetBinContent(bin);
    double fitValue = fit->Eval(x);

    // Calculate the residual as (data - fit)
    double residual = dataValue - fitValue;
    double error = hData->GetBinError(bin); // Assuming Poisson statistics, you can also add fit and bg errors if needed

    hResiduals->SetBinContent(bin, residual);
    hResiduals->SetBinError(bin, error);
  }

  return hResiduals;
}

// add _bc to get best cluster branches from parse file
std::string addbc(std::string input) {
  const std::string suffix = "_bc";
  const std::string targets[5] = {"coin", "dy", "dx", "hcale", "hcalon"}; //branches that depend on hcal clusters

  for (const auto& target : targets) {
    std::string::size_type pos = 0;
    std::string token = target + suffix;  // Create the new token with the suffix

    // Continue searching the string for the target and replacing it
    while ((pos = input.find(target, pos)) != std::string::npos) {
      // Ensure we match whole words by checking character before and after the match
      bool match = true;
      if (pos != 0 && isalnum(input[pos - 1])) {
	match = false; // Check for character before
      }
      size_t endPos = pos + target.length();
      if (endPos < input.length() && isalnum(input[endPos])) {
	match = false; // Check for character after
      }

      if (match) {
	input.replace(pos, target.length(), token);
	pos += token.length(); // Move past the newly added part to avoid infinite loops
      } else {
	pos += target.length(); // Move past the current word if it's part of a larger word
      }
    }
  }

  return input;
}


// Function to split a string by a delimiter and return a vector of tokens
std::vector<std::string> split(const std::string &s, char delimiter) {
  std::vector<std::string> tokens;
  std::string token;
  std::istringstream tokenStream(s);
  while (std::getline(tokenStream, token, delimiter)) {
    tokens.push_back(token);
  }
  return tokens;
}

// Function to get the cell content based on kine, mag, target, and column name
std::map<std::string, std::string> getRowContents(const std::string &filePath, int kine, int mag, const std::string &target, const std::vector<std::string> &excludeKeys) {
  std::ifstream file(filePath);
  std::string line;
  std::map<std::string, std::string> rowContents;
  std::vector<std::string> headers;
  std::set<std::string> excludeSet(excludeKeys.begin(), excludeKeys.end());

  // Read the header line
  if (std::getline(file, line)) {
    headers = split(line, ',');
  }

  // Read the rest of the lines
  while (std::getline(file, line)) {
    std::vector<std::string> values = split(line, ',');
    if (values.size() < 3) {
      continue;
    }

    // Check if the row matches the specified kine, mag, and target
    if (std::stoi(values[0]) == kine && std::stoi(values[1]) == mag && values[2] == target) {
      // Store all column values except the first three in a map
      for (size_t i = 3; i < values.size(); ++i) {
	if (i < headers.size() && excludeSet.find(headers[i]) == excludeSet.end()) {
	  rowContents[headers[i]] = values[i];
	}
      }
      break;
    }
  }

  return rowContents;
}
