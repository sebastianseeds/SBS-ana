//sseeds 
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include <regex>
#include <vector>
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

double percent_of_data = 0.97;

//sbs4,30 sbs8,70, sbs9,70
vector<double> dx_mean = {0.7405,0.8776,0.9049};
vector<double> dy_mean = {0.0462,0.0503,0.0292};

double binErrLimit = 0.1; //Set limit on HDE error to use HDE prot vs pos exp point for corrections

// Function to calculate N sigma for percent_of_data% coverage
auto findNSigma = [](TH1D* hist, double mean, double sigma) {
  int nbins = hist->GetNbinsX();
  double totalIntegral = hist->Integral();

  // Iterate over values of N until percent_of_data% coverage is achieved
  for (double N = 1; N < 10; N += 0.01) {
    int leftBin = hist->FindBin(mean - N * sigma);
    int rightBin = hist->FindBin(mean + N * sigma);
    double integral = hist->Integral(leftBin, rightBin);

    cout << leftBin << " " << rightBin << " " << integral / totalIntegral << endl;

    if (integral / totalIntegral >= percent_of_data) {
      return N;
    }
  }

  return 0.0;  // Fallback in case percent_of_data% is not found within 10 sigma
};

// Function to fit the peak region with a Gaussian
auto fitPeakRegion = [](TH1D* hist) {
  int maxBin = hist->GetMaximumBin();
  double maxValue = hist->GetBinContent(maxBin);
  double mean = hist->GetBinCenter(maxBin);

  // Find half maximum points
  double halfMax = maxValue / 2.0;
  int leftBin = maxBin, rightBin = maxBin;

  // Scan left to find where the histogram falls to half max
  while (leftBin > 0 && hist->GetBinContent(leftBin) > halfMax) {
    --leftBin;
  }

  // Scan right to find where the histogram falls to half max
  while (rightBin < hist->GetNbinsX() && hist->GetBinContent(rightBin) > halfMax) {
    ++rightBin;
  }

  double leftEdge = hist->GetBinCenter(leftBin);
  double rightEdge = hist->GetBinCenter(rightBin);

  // Fit Gaussian in this range
  TF1* gaus = new TF1("gaus", "gaus", leftEdge, rightEdge);
  hist->Fit(gaus, "RQ");

  return gaus;
};


// Function to concatenate cuts with their names
std::string concatenateCutsWithNames(const std::vector<std::pair<std::string, std::string>>& cutsWithNames) {
    std::string concatenatedCuts;
    for (const auto& cut : cutsWithNames) {
        concatenatedCuts += cut.first + "&&" + cut.second + "&&_____&&";
    }
    return concatenatedCuts;
}

//gets values for optics v x and y cuts for plotting
std::vector<double> extractValues(const std::string& expression) {
  std::vector<double> values;
  std::regex pattern("[<>]=?\\s*(-?\\d*\\.?\\d+)");
  std::smatch matches;

  std::string::const_iterator searchStart(expression.cbegin());
  while (std::regex_search(searchStart, expression.cend(), matches, pattern)) {
    values.push_back(std::stod(matches[1].str()));
    searchStart = matches.suffix().first;
  }

  return values;
}

//MAIN
void spotSizeCal(int kine=9, 
		 int mag=70, 
		 int pass=2, 
		 bool ellipse = true,
		 bool pidcuts = false,
		 bool effz=true) {

  //set draw params
  gStyle->SetPalette(55);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);
  
  cout << "Processing analysis of HCal detection efficiency for kinematic " << kine << " at field " << mag << "." << endl << endl;

  // reading json configuration file
  JSONManager *jmgr = new JSONManager("../../config/hdemap.json");

  //Get cuts
  std::string magcut = jmgr->GetValueFromSubKey_str( "mag_cuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded magnetic field cuts: " << magcut << endl;

  std::string trackcut = jmgr->GetValueFromSubKey_str( "track_v_cuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded track validity cuts: " << trackcut << endl;

  std::string optxcut = jmgr->GetValueFromSubKey_str( "optx_v_cuts", Form("sbs%d_%d",kine,mag) );
  vector<double> optxcut_vec = extractValues(optxcut);
  cout << "Loaded optics x validity cuts: " << optxcut << " (values: " << optxcut_vec[0] << ", " << optxcut_vec[1] << ")" << endl;

  std::string optycut = jmgr->GetValueFromSubKey_str( "opty_v_cuts", Form("sbs%d_%d",kine,mag) );
  vector<double> optycut_vec = extractValues(optycut);
  cout << "Loaded optics y validity cuts: " << optycut << " (values: " << optycut_vec[0] << ", " << optycut_vec[1] << ")" << endl;

  std::string fidxcut = jmgr->GetValueFromSubKey_str( "fidx_cut", Form("sbs%d_%d",kine,mag) );
  vector<double> fidxcut_vec = extractValues(fidxcut);
  cout << "Loaded fiducial x cut: " << fidxcut << " (values: " << fidxcut_vec[0] << ", " << fidxcut_vec[1] << ")" << endl;

  std::string fidycut = jmgr->GetValueFromSubKey_str( "fidy_cut", Form("sbs%d_%d",kine,mag) );
  vector<double> fidycut_vec = extractValues(fidycut);
  cout << "Loaded fiducial y cut: " << fidycut << " (values: " << fidycut_vec[0] << ", " << fidycut_vec[1] << ")" << endl;

  std::string earmcut = jmgr->GetValueFromSubKey_str( "earm_cuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded electron arm cuts: " << earmcut << endl;

  std::string earmpidcut = jmgr->GetValueFromSubKey_str( "earmpid_cuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded electron arm with PID cuts: " << earmcut << endl;

  //stitch together relevant cuts for later
  std::string e_cuts;
  std::string elas_cuts;
  if(pidcuts){
    e_cuts = magcut + "&&" + trackcut + "&&" + optxcut + "&&" + optycut + "&&" + earmpidcut;
    elas_cuts = earmpidcut;
  }else{
    e_cuts = magcut + "&&" + trackcut + "&&" + optxcut + "&&" + optycut + "&&" + earmcut;
    elas_cuts = earmcut;
  }
  
  std::string tra_cuts = magcut + "&&" + trackcut;
  std::string opt_cuts = magcut + "&&" + optxcut + "&&" + optycut;

  std::string fid_cuts = fidxcut + "&&" + fidycut;

  //cuts for exracting efficiency
  std::string norm_cuts = e_cuts + "&&" + fidxcut + "&&" + fidycut;

  std::string norm_xcuts = e_cuts + "&&" + fidycut;

  std::string norm_ycuts = e_cuts + "&&" + fidxcut;
  
  //Get plot limits/bins
  int bbmidxbins = jmgr->GetValueFromSubKey<int>( "bbbinsx", Form("sbs%d_%d",kine,mag) );
  int bbmidybins = jmgr->GetValueFromSubKey<int>( "bbbinsy", Form("sbs%d_%d",kine,mag) );
  int hcalxbins = jmgr->GetValueFromSubKey<int>( "hbinsx", Form("sbs%d_%d",kine,mag) );
  int hcalybins = jmgr->GetValueFromSubKey<int>( "hbinsy", Form("sbs%d_%d",kine,mag) );
  int W2bins = jmgr->GetValueFromSubKey<int>( "W2bins", Form("sbs%d_%d",kine,mag) );

  std::pair<double,double> bbmidxlim = {
    jmgr->GetValueFromSubKey<double>( "bbmidx_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "bbmidx_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded bbmidx plot limits min: " << bbmidxlim.first << ", max: " << bbmidxlim.second << ", with bins: " << bbmidxbins << endl;

  std::pair<double,double> bbmidylim = {
    jmgr->GetValueFromSubKey<double>( "bbmidy_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "bbmidy_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded bbmidy plot limits min: " << bbmidylim.first << ", max: " << bbmidylim.second << ", with bins: " << bbmidybins <<  endl;

  std::pair<double,double> hcalxlim = {
    jmgr->GetValueFromSubKey<double>( "hcalx_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "hcalx_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded hcalx plot limits min: " << hcalxlim.first << ", max: " << hcalxlim.second << ", with bins: " << hcalxbins <<  endl;

  std::pair<double,double> hcalylim = {
    jmgr->GetValueFromSubKey<double>( "hcaly_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "hcaly_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded hcaly plot limits min: " << hcalylim.first << ", max: " << hcalylim.second << ", with bins: " << hcalybins <<  endl;

  std::pair<double,double> W2lim = {
    jmgr->GetValueFromSubKey<double>( "W2_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "W2_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded W2 plot limits min: " << W2lim.first << ", max: " << W2lim.second << ", with bins: " << W2bins <<  endl;
 
  //set up files and paths
  std::string effz_word = "";
  if(effz)
    effz_word = "_effz";
  std::string tar_word = "_lh2";

  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fin_path = outdir_path + Form("/parse/parse_sbs%d_pass%d_barebones%s%s.root",kine,pass,effz_word.c_str(),tar_word.c_str());
  std::string fout_path = outdir_path + Form("/gmn_analysis/misc/spotSizeCal_sbs%d_mag%d_pass%d%s.root",kine,mag,pass,effz_word.c_str());

  cout << "Setting up output path: " << fout_path << endl;


  // Open the ROOT file
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

  // Open an output ROOT file
  TFile* outputFile = new TFile(fout_path.c_str(), "RECREATE");

  ///////////////////////
  ///////////////////////
  ///////////////////////
  ///////////////////////
  ///////////////////////
  ///////////////////////
  ////Begin Plots////////
  ///////////////////////
  ///////////////////////
  ///////////////////////


  // Create a new canvas
  TCanvas *cdxdy = new TCanvas("cdxdy", "dx vs dy proj", 1200, 600);
  cdxdy->Divide(2,1);
  cdxdy->cd(1);

  // Draw dx vs dy plot for comparison with spot cut region
  std::string dxdyName = "hdxdy";
  std::string dxdycut = norm_cuts;
  std::string dxdyTitle = "#Delta x vs #Delta y (Validity, Fiducial, and Elastic Cuts)";
  TH2D* hdxdy = new TH2D( dxdyName.c_str(), 
			  (dxdyTitle + ";#Delta y (m);#Delta x (m)").c_str(), 
			  hcalybins, 
			  hcalylim.first, 
			  hcalylim.second, 
			  hcalxbins, 
			  hcalxlim.first, 
			  hcalxlim.second);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("dx_bc:dy_bc>>" + dxdyName).c_str(), dxdycut.c_str(), "COLZ");

  // Create projections onto x and y axes
  TH1D *hdx_proj = hdxdy->ProjectionX("hdx_proj");
  TH1D *hdy_proj = hdxdy->ProjectionY("hdy_proj");

  // Fit the peak regions
  TF1 *gaus_x = fitPeakRegion(hdx_proj);
  double mean_x = gaus_x->GetParameter(1);
  double sigma_x = gaus_x->GetParameter(2);

  TF1 *gaus_y = fitPeakRegion(hdy_proj);
  double mean_y = gaus_y->GetParameter(1);
  double sigma_y = gaus_y->GetParameter(2);

  hdx_proj->Draw("hist");
  gaus_x->Draw("same");

  // Add legend for x projection
  TLegend *leg_x = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg_x->AddEntry(gaus_x, Form("Mean: %.4f", mean_x), "l");
  leg_x->AddEntry(gaus_x, Form("Sigma: %.4f", sigma_x), "l");
  leg_x->Draw();

  cdxdy->cd(2);
  hdy_proj->Draw("hist");
  gaus_y->Draw("same");

  // Add legend for y projection
  TLegend *leg_y = new TLegend(0.6, 0.7, 0.9, 0.9);
  leg_y->AddEntry(gaus_y, Form("Mean: %.4f", mean_y), "l");
  leg_y->AddEntry(gaus_y, Form("Sigma: %.4f", sigma_y), "l");
  leg_y->Draw();

  // Create a new canvas
  TCanvas *cdxdy2 = new TCanvas("cdxdy2", "dx vs dy", 1200, 600);
  cdxdy2->Divide(2,1);
  cdxdy2->cd(1);

  // Find N sigma for 99% coverage
  double nsigma_x = findNSigma(hdx_proj, mean_x, sigma_x);
  double nsigma_y = findNSigma(hdy_proj, mean_y, sigma_y);

  // cout << endl << "mean_x = " << mean_x << endl;
  // cout << "mean_y = " << mean_x << endl << endl;

  // cout << "nsigma_x = " << nsigma_x << " (sigma_x = " << sigma_x << ")" << endl;
  // cout << "nsigma_y = " << nsigma_x << " (sigma_y = " << sigma_x << ")" << endl;

  // Draw the dxdy plot without spot cut
  hdxdy->Draw("COLZ");

  // Draw dx vs dy plot for comparison with spot cut region
  std::string dxdyspotName = "hdxdyspot";
  //std::string dxdyspotcut = norm_cuts + "&&((pow(dx_bc + " + std::to_string(mean_x) + ", 2) / pow(" + std::to_string(nsigma_x*sigma_x) + ", 2)) + (pow(dy_bc + " + std::to_string(mean_x) + ", 2) / pow(" + std::to_string(nsigma_y*sigma_y) + ", 2)) <= 1)";
  std::string dxdyspotcut = norm_cuts + "&&((pow(dx_bc - " + std::to_string(mean_y) + ", 2) / pow(" + std::to_string(nsigma_y*sigma_y) + ", 2)) + (pow(dy_bc - " + std::to_string(mean_x) + ", 2) / pow(" + std::to_string(nsigma_x*sigma_x) + ", 2)) <= 1)";
  std::string dxdyspotTitle = "#Delta x vs #Delta y (Validity, Fiducial, Elastic, and Spot Cuts)";
  TH2D* hdxdyspot = new TH2D( dxdyspotName.c_str(),
  			      (dxdyspotTitle + ";#Delta y (m); #Delta x (m)").c_str(), 
  			      hcalybins, 
  			      hcalylim.first, 
  			      hcalylim.second, 
  			      hcalxbins, 
  			      hcalxlim.first, 
  			      hcalxlim.second);

  // Draw the plot using the created histogram
  tree->Draw(("dx_bc:dy_bc>>" + dxdyspotName).c_str(), dxdyspotcut.c_str(), "COLZ");

  // Draw the dxdy plot without spot cut
  hdxdy->Draw("COLZ");

  // Draw the dxdy plot with the spot cut
  cdxdy2->cd(2);
  hdxdyspot->Draw("colz");

  cdxdy2->Update();
  cdxdy2->Write();

  // Define all cuts and their corresponding names
  std::vector<std::pair<std::string, std::string>> cutsWithNames = {
    {"Track Validity Cuts", tra_cuts},
    // {"Track and Optical x Validity Cuts", W2optycut}, // Uncomment if needed
    {"Optical Validity Cuts", opt_cuts},
    {"Elastic Cuts", elas_cuts},
    //{"Spot Cuts", h_spotcut},
    {"Fiducial Cuts", fid_cuts}
  };
  
  // Concatenate all cuts with their names
  std::string allCuts = concatenateCutsWithNames(cutsWithNames);

  std::string spotoptpart = "&&_____&&Optimized Spot Sigma&&_____&&percent of data in cut: " + std::to_string(100.*percent_of_data) +
    "&&dx mean/sigma/N_sigma: " + std::to_string(mean_y) + "/" + std::to_string(sigma_y) + "/" + std::to_string(nsigma_y) + 
    "&&dy mean/sigma/N_sigma: " + std::to_string(mean_x) + "/" + std::to_string(sigma_x) + "/" + std::to_string(nsigma_x);

  allCuts += spotoptpart;

  // Display all concatenated cuts on a single canvas
  std::string disnameall = "All Cuts";
  TCanvas *cdisall = new TCanvas("cdis", disnameall.c_str(), 1000, 800);
  util::parseAndDisplayCuts(disnameall.c_str(), allCuts.c_str(), cdisall);

  cout << "All plots created. Output file located here: " << fout_path << endl;

}
