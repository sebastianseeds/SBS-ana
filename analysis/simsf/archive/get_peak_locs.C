//sseeds 11.9.23 Script to run over several outputs from sf_dxdy.C, extract proton peak locations, plot means on TGraph, fit TGraph, and get fit params with multiple fits. Also plots all proton peaks on same histogram along with data proton peak in different color.

#include <vector>
#include <iostream>
#include <algorithm>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TVector3.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TF1.h"
#include "TStopwatch.h"

const int NFiles = 8;
const double gen_sig = 0.18;

void createAndFitGraph(const std::vector<double>& gx, const std::vector<double>& gy, TCanvas* canvas) {
    if (gx.size() != gy.size() || gx.empty()) {
        std::cerr << "Vectors must be of the same size and not empty." << std::endl;
        return;
    }
    if (!canvas) {
        std::cerr << "Canvas is null!" << std::endl;
        return;
    }

    // Create TGraph
    TGraph *graph = new TGraph(gx.size(), gx.data(), gy.data());

    // Fit with a straight line
    TF1 *fitFunc = new TF1("fitFunc", "pol1", gx.front(), gx.back()); // pol1 for linear fit
    graph->Fit(fitFunc, "Q"); // "Q" for quiet mode (no print)

    // Print slope and y-intercept
    double slope = fitFunc->GetParameter(1);
    double yIntercept = fitFunc->GetParameter(0);
    std::cout << "Slope: " << slope << ", Y-intercept: " << yIntercept << std::endl;

    // Draw graph
    canvas->cd(); // Change to the provided canvas
    graph->Draw("AP"); // "A" to draw the axes, "P" to draw the graph as points

    canvas->Update();
    //canvas->Update();

    // Clean up
    //delete fitFunc;
    //delete graph;
}

double fitGaussianAndGetMean(TH1D* hist, double sig) {
    if (!hist) {
        std::cerr << "Histogram is null!" << std::endl;
        return -1;
    }

    // Find the bin with the maximum value
    int binMax = hist->GetMaximumBin();
    double xMax = hist->GetXaxis()->GetBinCenter(binMax);

    // Define the fit range
    double fitMin = xMax - sig;
    double fitMax = xMax + sig;

    // Fit the histogram within the specified range
    TF1 *gausFit = new TF1("gausFit", "gaus", fitMin, fitMax);
    hist->Fit(gausFit, "RQ"); // "R" for fit range, "Q" for quiet mode (no print)

    // Get the mean from the fit
    double mean = gausFit->GetParameter(1); // Parameter 1 is the mean of the Gaussian

    // Clean up
    delete gausFit;

    return mean;
}

void get_peak_locs(){  
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  //define input file directory
  string inputfiledir = "/volatile/halla/sbs/seeds/scale_field_sbs9/";

  //define output file name
  string outputfilename = inputfiledir + "analysis/sf_graph.root";

  //set up output file
  TFile *fout = new TFile( outputfilename.c_str(), "RECREATE" );

  //get dx plots for proton
  TFile *sf_file[NFiles];
  string inputfiles[NFiles];
  TH1D *hdx_p[NFiles];
  vector<double> means_p;
  vector<double> field_scales;

  //set up the holding canvas
  TCanvas *c0 = new TCanvas("c0","all proton peaks",1600,1200);
  c0->cd();
  gStyle->SetOptStat(0);

  //read in the files, retrieve some histograms, get means and field scales
  int sf = 0;
  for( int i=0; i<NFiles; ++i ){
    field_scales.push_back(i);
    inputfiles[i] = inputfiledir + Form("plots/sf%d_dxdy_out.root",sf*10);
    sf++;
    sf_file[i] = new TFile(inputfiles[i].c_str(), "READ"); // Open each ROOT file in read mode
    hdx_p[i] = (TH1D*)sf_file[i]->Get("hdx_HCAL_p");
    means_p.push_back(fitGaussianAndGetMean(hdx_p[i],gen_sig));
    hdx_p[i]->SetLineColor(kRed+i);
    hdx_p[i]->Draw("same");
  }  

  for( int i=0; i<NFiles; ++i )
    hdx_p[i]->Write();
  
  //set up the canvas
  TCanvas *c1 = new TCanvas("c1","scale field comparisons",1600,1200);
  c1->cd();
  gStyle->SetOptStat(0);

  //create tgraph
  createAndFitGraph(field_scales,means_p,c1);

  c1->Write();

  fout->Write();

  cout << "Analysis complete. Outfile located at " << outputfilename << endl;

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;    
 

}
