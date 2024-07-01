//sseeds
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
#include <regex>
#include "../../include/gmn.h"
#include "../../src/jsonmgr.C"

//fit params
double xrange_l = -3.0;
double xrange_h = 4.0;
int xbins = 600;

double yrange_l = -2.0;
double yrange_h = 2.0;
int ybins = 400;

//dxdy fit params
double dxrange_l = -2.0;
double dxrange_h = 1.0;
int dxbins = 300;

double dyrange_l = -1;
double dyrange_h = 1;
int dybins = 200;

//HCal dims
double ymin = econst::hcalposXi_mc;
double ymax = econst::hcalposXf_mc;
double xmin = econst::hcalposYi_mc;
double xmax = econst::hcalposYf_mc;

// Utility function to extract boundary values from cut strings using regular expressions
std::vector<double> extractBoundaries(const std::string& cut) {
    std::vector<double> boundaries;
    std::regex re("[-+]?[0-9]*\\.?[0-9]+");
    std::sregex_iterator next(cut.begin(), cut.end(), re);
    std::sregex_iterator end;
    while (next != end) {
        std::smatch match = *next;
        try {
            boundaries.push_back(std::stod(match.str()));
        } catch (const std::invalid_argument& e) {
            std::cerr << "Error: Unable to convert boundary value from string: " << match.str() << std::endl;
        }
        next++;
    }
    return boundaries;
}


// Utility function to extract parameters from the spot cut string
void extractEllipseParams(const std::string& cut, double& dx_center, double& dy_center, double& dx_sigma, double& dy_sigma) {
  std::regex re_dx_center("dx_bc \\+ ([\\d\\.-]+)");
  std::regex re_dy_center("dy_bc \\+ ([\\d\\.-]+)");
  std::regex re_dx_sigma("/ pow\\((\\d+\\.\\d+), 2\\)");
  std::regex re_dy_sigma("/ pow\\((\\d+\\.\\d+), 2\\)");

  std::smatch match;

  if (std::regex_search(cut, match, re_dx_center) && match.size() > 1) {
    dx_center = -std::stod(match.str(1)); // negate the center value as per the provided cut
  }
  if (std::regex_search(cut, match, re_dy_center) && match.size() > 1) {
    dy_center = -std::stod(match.str(1)); // negate the center value as per the provided cut
  }
  if (std::regex_search(cut, match, re_dx_sigma) && match.size() > 1) {
    dx_sigma = std::stod(match.str(1));
  }
  if (std::regex_search(cut, match, re_dy_sigma) && match.size() > 1) {
    dy_sigma = std::stod(match.str(1));
  }
}

//main
void fidnew(int kine=8, 
	    int mag=100, 
	    int pass=2, 
	    bool effz=true,
	    bool alt = true) {
  //set draw params
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(0);
  gStyle->SetEndErrorSize(0);
  gStyle->SetPalette(55);
  gStyle->SetOptStat(0110);
  gStyle->SetStatTextColor(kBlack);
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.08);
  gStyle->SetStatX(0.35);
  gStyle->SetStatY(0.9);

  //get globalcut. Select magnetic field setting and wide elastics
  std::string gcut = Form("W2<1.2&&abs(coin)<6&&mag==%d",mag);

  //set up files and paths
  std::string effz_word = "";
  if(effz)
    effz_word = "_effz";
  std::string tar_word = "_ld2";

  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fin_path = outdir_path + Form("/parse/parse_sbs%d_pass%d_barebones%s%s.root",kine,pass,effz_word.c_str(),tar_word.c_str());

  std::string fout_path = outdir_path + Form("/gmn_analysis/fidnew_sbs%d_mag%d_pass%d%s.root",kine,mag,pass,effz_word.c_str());

  cout << "Setting up output path: " << fout_path << endl;

  // reading json configuration file
  //JSONManager *jmgr = new JSONManager("../../config/hdemap.json");
  JSONManager *jmgr = new JSONManager("../../config/nphde.json");

  //proton spotcut
  std::string p_spotcut = jmgr->GetValueFromSubKey_str( "h_dp_spot_cuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded ellipse proton delta cuts: " << p_spotcut << endl;

  //neutron spotcut
  std::string n_spotcut = jmgr->GetValueFromSubKey_str( "h_dn_spot_cuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded ellipse neutron delta cuts: " << n_spotcut << endl;

  //fiducial cut xexp
  std::string fidcut_x = jmgr->GetValueFromSubKey_str( "fidx_cut", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded fiducial xexp cut: " << fidcut_x << endl;

  //fiducial cut yexp
  std::string fidcut_y = jmgr->GetValueFromSubKey_str( "fidy_cut", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded fiducial yexp cut: " << fidcut_y << endl;

  //get fiducial boundaries
  std::vector<double> yBoundaries = extractBoundaries(fidcut_x);
  std::vector<double> xBoundaries = extractBoundaries(fidcut_y);

  cout << xBoundaries[0] << " " << xBoundaries[1] << endl;
  
  // Ensure that there are exactly two boundaries for each cut
  if (xBoundaries.size() != 2 || yBoundaries.size() != 2) {
    std::cerr << "Error: Invalid fiducial cut boundaries." << std::endl;
    return;
  }

  //both spotcut
  std::string N_spotcut = "((" + n_spotcut + ")||(" + p_spotcut + "))";

  cout << "Loaded ellipse proton and neutron delta cuts: " << N_spotcut << endl;

  // Extract parameters for the ellipses from the spot cuts
  double p_dx_center, p_dy_center, p_dx_sigma, p_dy_sigma;
  double n_dx_center, n_dy_center, n_dx_sigma, n_dy_sigma;

  extractEllipseParams(p_spotcut, p_dx_center, p_dy_center, p_dx_sigma, p_dy_sigma);
  extractEllipseParams(n_spotcut, n_dx_center, n_dy_center, n_dx_sigma, n_dy_sigma);

  // Create TEllipse objects
  TEllipse *p_ellipse = new TEllipse(p_dy_center, p_dx_center, p_dy_sigma, p_dx_sigma);
  p_ellipse->SetLineColor(kRed);
  p_ellipse->SetLineWidth(2);
  p_ellipse->SetFillStyle(0); // Hollow ellipse

  TEllipse *n_ellipse = new TEllipse(n_dy_center, n_dx_center, n_dy_sigma, n_dx_sigma);
  n_ellipse->SetLineColor(kBlue);
  n_ellipse->SetLineWidth(2);
  n_ellipse->SetFillStyle(0); // Hollow ellipse

  // Make TLines that form a box with HCal dims.
  TLine *line5 = new TLine(xmin, ymin, xmax, ymin);
  TLine *line6 = new TLine(xmax, ymin, xmax, ymax);
  TLine *line7 = new TLine(xmax, ymax, xmin, ymax);
  TLine *line8 = new TLine(xmin, ymax, xmin, ymin);

  line5->SetLineColor(kRed);
  line6->SetLineColor(kRed);
  line7->SetLineColor(kRed);
  line8->SetLineColor(kRed);
  line5->SetLineWidth(2);
  line6->SetLineWidth(2);
  line7->SetLineWidth(2);
  line8->SetLineWidth(2);
  line5->SetLineStyle(2);
  line6->SetLineStyle(2);
  line7->SetLineStyle(2);
  line8->SetLineStyle(2);

  // More Tlines for the fiducial cuts
  TLine* line10 = new TLine(xBoundaries[0], yBoundaries[0], xBoundaries[1], yBoundaries[0]);
  TLine* line20 = new TLine(xBoundaries[1], yBoundaries[0], xBoundaries[1], yBoundaries[1]);
  TLine* line30 = new TLine(xBoundaries[1], yBoundaries[1], xBoundaries[0], yBoundaries[1]);
  TLine* line40 = new TLine(xBoundaries[0], yBoundaries[1], xBoundaries[0], yBoundaries[0]);

  line10->SetLineColor(kBlue);
  line20->SetLineColor(kBlue);
  line30->SetLineColor(kBlue);
  line40->SetLineColor(kBlue);
  line10->SetLineWidth(2);
  line20->SetLineWidth(2);
  line30->SetLineWidth(2);
  line40->SetLineWidth(2);
  line10->SetLineStyle(2);
  line20->SetLineStyle(2);
  line30->SetLineStyle(2);
  line40->SetLineStyle(2);

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

  /////////////////////
  // Create Control
  ////////////////////

  // Create a new canvas
  TCanvas *cControl = new TCanvas("cControl", "Fiducial Control (no cuts)", 2100, 600);
  cControl->Divide(3,1);
  cControl->cd(1);

  // Draw dx vs dy with no cuts
  std::string controldxdyName = "control_dxdy";
  std::string controldxdyCut = gcut;
  std::string controldxdyTitle = "dx vs dy; dy (m); dx (m)";
  TH2D* controldxdy = new TH2D( controldxdyName.c_str(), 
			   controldxdyTitle.c_str(), 
			   dybins, 
			   dyrange_l, 
			   dyrange_h, 
			   dxbins, 
			   dxrange_l, 
			   dxrange_h);

  // Draw the plot using the created histogram
  tree->Draw(("dx_bc:dy_bc>>" + controldxdyName).c_str(), controldxdyCut.c_str(), "COLZ");

  // Draw the histogram on the canvas
  controldxdy->Draw("COLZ");

  // Draw the ellipses on the histogram
  p_ellipse->Draw("same");
  n_ellipse->Draw("same");

  // Write the histogram to the output file
  cControl->Update();
  cControl->Write();

  cControl->cd(2);

  // Draw expected x vs expected y (from qvector) with both spot cuts
  std::string controlexpName = "control_exp";
  std::string controlexpCut = gcut + "&&" + N_spotcut;
  std::string controlexpTitle = "x_{exp} vs y_{exp}; y_{exp} (m); x_{exp} (m)";
  TH2D* controlexp = new TH2D( controlexpName.c_str(), 
			   controlexpTitle.c_str(), 
			   ybins, 
			   yrange_l, 
			   yrange_h, 
			   xbins, 
			   xrange_l, 
			   xrange_h);

  // Draw the plot using the created histogram
  tree->Draw(("xexp:yexp>>" + controlexpName).c_str(), controlexpCut.c_str(), "COLZ");

  // Draw the histogram on the canvas
  controlexp->Draw("COLZ");

  line5->Draw("same");
  line6->Draw("same");
  line7->Draw("same");
  line8->Draw("same");

  cControl->Update();

  cControl->cd(3);

  gPad->SetGridx();
  gPad->SetGridy();

  // Draw expected x with both spot cuts
  std::string controlxexpName = "control_xexp";
  std::string controlxexpCut = gcut + "&&" + N_spotcut;
  std::string controlxexpTitle = "x_{exp}; x_{exp} (m)";
  TH1D* controlxexp = new TH1D( controlxexpName.c_str(), 
			       controlxexpTitle.c_str(), 
			       xbins, 
			       xrange_l, 
			       xrange_h);

  // Draw the plot using the created histogram
  tree->Draw(("xexp>>" + controlxexpName).c_str(), controlxexpCut.c_str(), "COLZ");

  // Draw the histogram on the canvas
  controlxexp->Draw("COLZ");

  double controlxexp_max = controlxexp->GetMaximum();
  TLine *vline5 = new TLine(ymin,0,ymin,controlxexp_max);
  vline5->SetLineColor(kRed);
  vline5->SetLineWidth(2);
  vline5->SetLineStyle(2);
  
  TLine *vline6 = new TLine(ymax,0,ymax,controlxexp_max);
  vline6->SetLineColor(kRed);
  vline6->SetLineWidth(2);
  vline6->SetLineStyle(2);

  TLine *vline50 = new TLine(yBoundaries[0],0,yBoundaries[0],controlxexp_max);
  vline50->SetLineColor(kBlue);
  vline50->SetLineWidth(2);
  vline50->SetLineStyle(2);
  
  TLine *vline60 = new TLine(yBoundaries[1],0,yBoundaries[1],controlxexp_max);
  vline60->SetLineColor(kBlue);
  vline60->SetLineWidth(2);
  vline60->SetLineStyle(2);

  vline5->Draw("same");
  vline6->Draw("same");

  // Write the histogram to the output file
  cControl->Update();
  cControl->Write();

  ///////////////////////
  // Create Proton fid
  ////////////////////////

  // Create a new canvas
  TCanvas *cProton = new TCanvas("cProton", "Fiducial Proton (no cuts)", 2100, 600);
  cProton->Divide(3,1);
  cProton->cd(1);

  // Draw dx vs dy with no cuts
  std::string protondxdyName = "proton_dxdy";
  std::string protondxdyCut = gcut + "&&" + p_spotcut;
  std::string protondxdyTitle = "dx vs dy (proton spotcut); dy (m); dx (m)";
  TH2D* protondxdy = new TH2D( protondxdyName.c_str(), 
			   protondxdyTitle.c_str(), 
			   dybins, 
			   dyrange_l, 
			   dyrange_h, 
			   dxbins, 
			   dxrange_l, 
			   dxrange_h);

  // Draw the plot using the created histogram
  tree->Draw(("dx_bc:dy_bc>>" + protondxdyName).c_str(), protondxdyCut.c_str(), "COLZ");

  // Draw the histogram on the canvas
  protondxdy->Draw("COLZ");

  // Write the histogram to the output file
  cProton->Update();
  cProton->Write();

  cProton->cd(2);

  // Draw expected x vs expected y (from qvector) with no cuts
  std::string protonexpName = "proton_exp";
  std::string protonexpCut = gcut + "&&" + p_spotcut;
  std::string protonexpTitle = "x_{exp} vs y_{exp}; y_{exp} (m); x_{exp} (m)";
  TH2D* protonexp = new TH2D( protonexpName.c_str(), 
			   protonexpTitle.c_str(), 
			   ybins, 
			   yrange_l, 
			   yrange_h, 
			   xbins, 
			   xrange_l, 
			   xrange_h);

  // Draw the plot using the created histogram
  tree->Draw(("xexp:yexp>>" + protonexpName).c_str(), protonexpCut.c_str(), "COLZ");

  // Draw the histogram on the canvas
  protonexp->Draw("COLZ");

  // Make TLines that form a box on the histogram with HCal dims.
  line5->Draw("same");
  line6->Draw("same");
  line7->Draw("same");
  line8->Draw("same");
  line10->Draw("same");
  line20->Draw("same");
  line30->Draw("same");
  line40->Draw("same");

  cProton->Update();

  cProton->cd(3);

  gPad->SetGridx();
  gPad->SetGridy();

  // Draw expected x with both spot cuts
  std::string protonxexpName = "proton_xexp";
  std::string protonxexpCut = gcut + "&&" + p_spotcut;
  std::string protonxexpTitle = "x_{exp}; x_{exp} (m)";
  TH1D* protonxexp = new TH1D( protonxexpName.c_str(), 
			       protonxexpTitle.c_str(), 
			       xbins, 
			       xrange_l, 
			       xrange_h);

  // Draw the plot using the created histogram
  tree->Draw(("xexp>>" + protonxexpName).c_str(), protonxexpCut.c_str(), "COLZ");

  // Draw the histogram on the canvas
  protonxexp->Draw("COLZ");

  double protonxexp_max = protonxexp->GetMaximum();
  TLine *vline1 = new TLine(ymin,0,ymin,protonxexp_max);
  vline1->SetLineColor(kRed);
  vline1->SetLineWidth(2);
  vline1->SetLineStyle(2);
  
  TLine *vline2 = new TLine(ymax,0,ymax,protonxexp_max);
  vline2->SetLineColor(kRed);
  vline2->SetLineWidth(2);
  vline2->SetLineStyle(2);

  TLine *vline10 = new TLine(yBoundaries[0],0,yBoundaries[0],protonxexp_max);
  vline10->SetLineColor(kBlue);
  vline10->SetLineWidth(2);
  vline10->SetLineStyle(2);
  
  TLine *vline20 = new TLine(yBoundaries[1],0,yBoundaries[1],protonxexp_max);
  vline20->SetLineColor(kBlue);
  vline20->SetLineWidth(2);
  vline20->SetLineStyle(2);

  vline1->Draw("same");
  vline2->Draw("same");

  // Write the histogram to the output file
  cProton->Update();
  cProton->Write();

  ///////////////////////
  // Create Neutron fid
  ////////////////////////

  // Create a new canvas
  TCanvas *cNeutron = new TCanvas("cNeutron", "Fiducial Neutron (no cuts)", 2100, 600);
  cNeutron->Divide(3,1);
  cNeutron->cd(1);

  // Draw dx vs dy with no cuts
  std::string neutrondxdyName = "neutron_dxdy";
  std::string neutrondxdyCut = gcut + "&&" + n_spotcut;
  std::string neutrondxdyTitle = "dx vs dy (neutron spotcut); dy (m); dx (m)";
  TH2D* neutrondxdy = new TH2D( neutrondxdyName.c_str(), 
			   neutrondxdyTitle.c_str(), 
			   dybins, 
			   dyrange_l, 
			   dyrange_h, 
			   dxbins, 
			   dxrange_l, 
			   dxrange_h);

  // Draw the plot using the created histogram
  tree->Draw(("dx_bc:dy_bc>>" + neutrondxdyName).c_str(), neutrondxdyCut.c_str(), "COLZ");

  // Draw the histogram on the canvas
  neutrondxdy->Draw("COLZ");

  // Write the histogram to the output file
  cNeutron->Update();
  cNeutron->Write();

  cNeutron->cd(2);

  // Draw expected x vs expected y (from qvector) with no cuts
  std::string neutronexpName = "neutron_exp";
  std::string neutronexpCut = gcut + "&&" + n_spotcut;
  std::string neutronexpTitle = "x_{exp} vs y_{exp}; y_{exp} (m); x_{exp} (m)";
  TH2D* neutronexp = new TH2D( neutronexpName.c_str(), 
			   neutronexpTitle.c_str(), 
			   ybins, 
			   yrange_l, 
			   yrange_h, 
			   xbins, 
			   xrange_l, 
			   xrange_h);

  // Draw the plot using the created histogram
  tree->Draw(("xexp:yexp>>" + neutronexpName).c_str(), neutronexpCut.c_str(), "COLZ");

  // Draw the histogram on the canvas
  neutronexp->Draw("COLZ");

  // Make TLines that form a box on the histogram with HCal dims.
  line5->Draw("same");
  line6->Draw("same");
  line7->Draw("same");
  line8->Draw("same");
  line10->Draw("same");
  line20->Draw("same");
  line30->Draw("same");
  line40->Draw("same");

  cNeutron->Update();

  cNeutron->cd(3);

  gPad->SetGridx();
  gPad->SetGridy();

  // Draw expected x with both spot cuts
  std::string neutronxexpName = "neutron_xexp";
  std::string neutronxexpCut = gcut + "&&" + n_spotcut;
  std::string neutronxexpTitle = "x_{exp}; x_{exp} (m)";
  TH1D* neutronxexp = new TH1D( neutronxexpName.c_str(), 
			       neutronxexpTitle.c_str(), 
			       xbins, 
			       xrange_l, 
			       xrange_h);

  // Draw the plot using the created histogram
  tree->Draw(("xexp>>" + neutronxexpName).c_str(), neutronxexpCut.c_str(), "COLZ");

  // Draw the histogram on the canvas
  neutronxexp->Draw("COLZ");

  double neutronxexp_max = neutronxexp->GetMaximum();
  TLine *vline3 = new TLine(ymin,0,ymin,neutronxexp_max);
  vline3->SetLineColor(kRed);
  vline3->SetLineWidth(2);
  vline3->SetLineStyle(2);
  
  TLine *vline4 = new TLine(ymax,0,ymax,neutronxexp_max);
  vline4->SetLineColor(kRed);
  vline4->SetLineWidth(2);
  vline4->SetLineStyle(2);

  TLine *vline30 = new TLine(yBoundaries[0],0,yBoundaries[0],neutronxexp_max);
  vline30->SetLineColor(kBlue);
  vline30->SetLineWidth(2);
  vline30->SetLineStyle(2);
  
  TLine *vline40 = new TLine(yBoundaries[1],0,yBoundaries[1],neutronxexp_max);
  vline40->SetLineColor(kBlue);
  vline40->SetLineWidth(2);
  vline40->SetLineStyle(2);

  vline3->Draw("same");
  vline4->Draw("same");

  // Write the histogram to the output file
  cNeutron->Update();
  cNeutron->Write();

  //////////////////////
  //Draw them all together
  ////////////////////////


  // Create a new canvas
  TCanvas *cTotal = new TCanvas("cTotal", "Fiducial Total (no cuts)", 2100, 1800);
  cTotal->Divide(3,3);
  cTotal->cd(1);

  // Draw control dxdy
  controldxdy->Draw("COLZ");

  cTotal->Update();

  cTotal->cd(2);

  // Draw control xexp vs yexp
  controlexp->Draw("COLZ");

  line5->Draw("same");
  line6->Draw("same");
  line7->Draw("same");
  line8->Draw("same");

  // Create and draw the legend
  TLegend *l1 = new TLegend(0.5, 0.7, 0.9, 0.8); // Adjust the position as needed
  l1->SetBorderSize(0);
  l1->SetFillStyle(0);
  l1->AddEntry(line5, "HCal Full Acceptance", "l");
  l1->Draw("same");

  cTotal->Update();

  cTotal->cd(3);

  gPad->SetGridx();
  gPad->SetGridy();

  // Draw control xexp
  controlxexp->SetLineColor(kBlack);
  controlxexp->SetLineWidth(2);
  controlxexp->SetFillColor(kBlue-10);
  controlxexp->SetFillStyle(3003);
  controlxexp->GetYaxis()->SetRangeUser(0,1.3*controlxexp->GetMaximum());
  controlxexp->Draw("COLZ");

  vline5->Draw("same");
  vline6->Draw("same");

  l1->Draw("same");

  cTotal->Update();
  cTotal->cd(4);

  // Draw proton dxdy
  protondxdy->Draw("COLZ");

  cTotal->Update();

  cTotal->cd(5);

  // Draw proton xexp vs yexp
  protonexp->Draw("COLZ");

  line5->Draw("same");
  line6->Draw("same");
  line7->Draw("same");
  line8->Draw("same");
  line10->Draw("same");
  line20->Draw("same");
  line30->Draw("same");
  line40->Draw("same");

  // Create and draw the legend
  TLegend *l2 = new TLegend(0.5, 0.8, 0.9, 0.9); // Adjust the position as needed
  l2->SetBorderSize(0);
  l2->SetFillStyle(0);
  l2->AddEntry(line5, "HCal Full Acceptance", "l");
  l2->AddEntry(line10, "Fiducial Cuts", "l");
  l2->Draw("same");

  cTotal->Update();

  cTotal->cd(6);

  gPad->SetGridx();
  gPad->SetGridy();

  // Draw proton xexp
  protonxexp->SetLineColor(kBlack);
  protonxexp->SetLineWidth(2);
  protonxexp->SetFillColor(kBlue-10);
  protonxexp->SetFillStyle(3003);
  protonxexp->GetYaxis()->SetRangeUser(0,1.3*protonxexp->GetMaximum());
  protonxexp->Draw("COLZ");

  vline1->Draw("same");
  vline2->Draw("same");
  vline10->Draw("same");
  vline20->Draw("same");

  l2->Draw("same");

  cTotal->Update();
  cTotal->cd(7);

  // Draw neutron dxdy
  neutrondxdy->Draw("COLZ");

  cTotal->Update();

  cTotal->cd(8);

  // Draw neutron xexp vs yexp
  neutronexp->Draw("COLZ");

  line5->Draw("same");
  line6->Draw("same");
  line7->Draw("same");
  line8->Draw("same");
  line10->Draw("same");
  line20->Draw("same");
  line30->Draw("same");
  line40->Draw("same");

  l2->Draw("same");

  cTotal->Update();

  cTotal->cd(9);

  gPad->SetGridx();
  gPad->SetGridy();

  // Draw neutron xexp
  neutronxexp->SetLineColor(kBlack);
  neutronxexp->SetLineWidth(2);
  neutronxexp->SetFillColor(kBlue-10);
  neutronxexp->SetFillStyle(3003);
  neutronxexp->GetYaxis()->SetRangeUser(0,1.3*neutronxexp->GetMaximum());
  neutronxexp->Draw("COLZ");

  vline3->Draw("same");
  vline4->Draw("same");
  vline30->Draw("same");
  vline40->Draw("same");

  l2->Draw("same");

  cTotal->Update();


  cTotal->Write();

  

}
