//sseeds 10.23.23 - Script to take output file from hde_mc_data.C, fit each bin in p, and produce energy means for efficiency analysis later.

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

const int nbin = 150; //Number of bins for hcal E vs nucleon p fits obtained to ensure 1000 events per bin

// Function to fit slices of a TH2D histogram with Gaussian
std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> FitSlices(TH2D* hist2D, 
										    int nbins, 
										    double fit_fac, 
										    string nucleon,
										    TFile *outputfile ) {
  std::vector<double> x_values;
  std::vector<double> fit_means;
  std::vector<double> fit_sigmas;
    
  //double fit_range = fit_fac*;  // Adjust this value for the desired fit range
  double fit_range=1.;

  std::vector<TH1D*> fittedSlices;

  // Loop over the specified number of bins
  for (int bin = 1; bin <= nbins; ++bin) {

    std::cout << "Processing p bin " << bin << "/" << nbins << "\r";
    std::cout.flush();

    // Calculate the center of the current bin
    double x_center = hist2D->GetXaxis()->GetBinCenter(bin);
    x_values.push_back(x_center);

    fit_range = fit_fac*x_center;

    // Get the projection of the current bin along the Y-axis
    TH1D* slice = hist2D->ProjectionY(Form("pSlice_%s_%d", nucleon.c_str(), bin), bin, bin);

    // Find the maximum bin content and its bin number
    int max_bin = slice->GetMaximumBin();
    double max_content = slice->GetBinContent(max_bin);
    double max_center = slice->GetXaxis()->GetBinCenter(max_bin);
    double first_center = slice->GetXaxis()->GetBinCenter(2);
    double last_center = slice->GetXaxis()->GetBinCenter(slice->GetNbinsX());

    // Define the fit range symmetrically around the maximum bin with E det eff in mind (sig/E=const)
    double fit_min = std::max(first_center, max_center - fit_range);
    double fit_max = std::min(last_center, max_center + fit_range);

    // Create and fit a Gaussian to the slice within the defined range
    TF1* fit_function = new TF1("gaussian", "gaus", fit_min, fit_max);
    fit_function->SetParameters(max_content, max_center, 0.1);  // Initial parameters

    TCanvas* c = new TCanvas();
    slice->Fit(fit_function, "RQ");

    fittedSlices.push_back(slice);

    // Retrieve the fit results (mean and sigma)
    double mean = fit_function->GetParameter(1);
    double sigma = fit_function->GetParameter(2);

    fit_means.push_back(mean);
    fit_sigmas.push_back(sigma);

    // Write to file
    slice->Write();

    // Clean up
    delete fit_function;
    delete slice;
    delete c;
  }

  return std::make_tuple(x_values, fit_means, fit_sigmas);
}

//Uses g4sbs replays of simulated data set containing pgun/ngun, zero field, SBS4 geometry
//MAIN. (no args)
void hde_mc_getEmean( )
{  
  
  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  std::string date = util::getDate();

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/smchde.json");
  
  string proton_rfile = jmgr->GetValueFromSubKey_str( "rootfiles", "proton" );
  string neutron_rfile = jmgr->GetValueFromSubKey_str( "rootfiles", "neutron" );
  string prot_param_path = jmgr->GetValueFromSubKey_str( "params", "proton" );
  string neut_param_path = jmgr->GetValueFromSubKey_str( "params", "neutron" );
  string infilename = jmgr->GetValueFromKey_str( "outfile_0" ); //This produces E means. Read in only first iteration

  double pmin = jmgr->GetValueFromKey<double>( "pmin" );
  double pmax = jmgr->GetValueFromKey<double>( "pmax" );
  double Emin = jmgr->GetValueFromKey<double>( "Emin" );
  double Emax = jmgr->GetValueFromKey<double>( "Emax" );
  double fit_fac = jmgr->GetValueFromKey<double>( "fit_fac" );
  double tfac = jmgr->GetValueFromKey<double>( "tfac" );
  int jmgr_nbin = jmgr->GetValueFromKey<int>( "nbin" );

  if( jmgr_nbin!=nbin ){
    cout << "ERROR: const int nbin not equal to common configuration json file value." << endl;
    return;
  }
 
  //set up output file
  //TFile *fout = new TFile( "/lustre19/expphy/volatile/halla/sbs/seeds/hde_analysis/hdemc_emean.root", "RECREATE" );

  gStyle->SetOptFit();
  gStyle->SetEndErrorSize(0);

  // open files and get histograms
  TFile *i0_file = TFile::Open(infilename.c_str(), "READ");
  TH2D *hEdepvP_p = (TH2D*)i0_file->Get("hEdepvP_p");
  TH2D *hEdepvP_n = (TH2D*)i0_file->Get("hEdepvP_n");

  //clone the Evp histos
  TH2D *cEvP_p = (TH2D*)(hEdepvP_p->Clone("cEvP_p"));
  TH2D *cEvP_n = (TH2D*)(hEdepvP_n->Clone("cEvP_n"));

  //i0_file->Close();

  TFile *fout = new TFile( "outfiles/hdemc_emean.root", "RECREATE");

  //make multigraph for both p and n Evp means/sigma graphs
  auto mgb = new TMultiGraph();

  int Nbins_p = cEvP_p->GetNbinsX();
  int Nbins_n = cEvP_n->GetNbinsX();

  if( Nbins_p!=nbin || Nbins_n!=nbin ){
    cerr << "ERROR: Nucleon vs P histogram bins should match general configuration parameter nbin, but do not" << endl;
    return;
  }

  cout << "Processing proton E vs p.." << endl;

  //Call slice and fit function for both p and n
  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> fitResults_p =
    FitSlices(cEvP_p, Nbins_p, fit_fac, "p", fout);

  cout << "Processing neutron E vs p.." << endl;

  std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> fitResults_n =
    FitSlices(cEvP_n, Nbins_n, fit_fac, "n", fout);

  //Get unique vectors for TGraphErrors
  std::vector<double> x_values_p = std::get<0>(fitResults_p);
  std::vector<double> fit_means_p = std::get<1>(fitResults_p);
  std::vector<double> fit_sigmas_p = std::get<2>(fitResults_p);

  std::vector<double> x_values_n = std::get<0>(fitResults_n);
  std::vector<double> fit_means_n = std::get<1>(fitResults_n);
  std::vector<double> fit_sigmas_n = std::get<2>(fitResults_n);

  //Build canvas for Evp p/n overlay
  TCanvas *c1 = new TCanvas("c1","HCal E vs Nucleon P",1600,1200);
  c1->SetGrid();
  c1->cd();

  //make first tgraph for proton and add to multigraph
  TGraphErrors* gp = new TGraphErrors(Nbins_p, &(x_values_p[0]), &(fit_means_p[0]), nullptr, &(fit_sigmas_p[0]));
  gp->SetTitle("Proton");
  gp->SetMarkerColor(kRed);
  gp->SetMarkerStyle(20);
  gp->SetMarkerSize(2);
  gp->SetLineColor(kRed);
  gp->SetLineWidth(2);
  //gp->Draw("AP");
  mgb->Add(gp);

  // //make second tgraph for neutron and add to multigraph
  TGraphErrors* gn = new TGraphErrors(Nbins_n, &(x_values_n[0]), &(fit_means_n[0]), nullptr, &(fit_sigmas_n[0]));
  gn->SetTitle("Neutron");
  gn->SetMarkerColor(kBlue);
  gn->SetMarkerStyle(22);
  gn->SetMarkerSize(2);
  gn->SetLineColor(kBlue);
  gn->SetLineWidth(2);
  //gn->Draw("AP");
  mgb->Add(gn);

  mgb->SetTitle("HCal E vs Nucleon p");
  mgb->GetXaxis()->SetTitle("Nucleon p (GeV)");
  mgb->GetYaxis()->SetTitle("E_{hcal}");
  mgb->Draw("AP");

  c1->Update();
  c1->Write();

  //write the E means to params/ for next iteration of hde_mc_data.C
  ofstream E_prot;
  E_prot.open( prot_param_path );
  E_prot << "#HCal E mean via gaussian fit to proton data obtained " << date.c_str() << endl;
  E_prot << "#" << endl;

  for( Int_t b=0; b<nbin; b++ ){   
    E_prot << fit_means_p[b] << endl;
  }
  E_prot.close();

  ofstream E_neut;
  E_neut.open( neut_param_path );
  E_neut << "#HCal E mean via gaussian fit to neutron data obtained " << date.c_str() << endl;
  E_neut << "#" << endl;

  for( Int_t b=0; b<nbin; b++ ){   
    E_neut << fit_means_n[b] << endl;
  }
  E_neut.close();

  //finish up
  fout->Write();

  cout << "Analysis complete. Output histograms sent to file outfiles/hdemc_emean.root" <<  endl;

  // Send time efficiency report to console
  st->Stop();
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
