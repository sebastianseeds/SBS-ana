//sseeds 11.9.23 Script to run over several outputs from sf_dxdy.C and extract the maximum variance between proton peak location and fiducial cut proton peak location among field settings

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
#include "../../include/gmn.h"

const int NFiles = 11;
const double gen_sig = 0.18;

//MAIN. (no args)
void sf_getFidVar(){  

  //define input file directory
  string inputfiledir = "/volatile/halla/sbs/seeds/scale_field_sbs9/";

  //define output file name
  string outputfilename = inputfiledir + "analysis/sf_getFitVar_out.root";

  //set up output file
  TFile *fout = new TFile( outputfilename.c_str(), "RECREATE" );

  double diff[NFiles];
  double sf[NFiles];

  double maxDiff = 0;

  //set up output for scale = 0 
  string inputfile_0 = inputfiledir + "plots/sf0_dxdy_out.root";
  TFile *sf0_file = new TFile(inputfile_0.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_0 = (TH1D*)sf0_file->Get("hdx_HCAL_p");
  TH1D *hdx_p_0_fid = (TH1D*)sf0_file->Get("hdx_HCAL_p_fid");

  diff[0] = util::fitGaussianAndGetFineMean(hdx_p_0_fid,gen_sig) - util::fitGaussianAndGetFineMean(hdx_p_0,gen_sig);
  sf[0] = 0;

  if( abs(diff[0])>maxDiff )
    maxDiff = abs(diff[0]);
  
  //set up output for scale = 10 
  string inputfile_10 = inputfiledir + "plots/sf10_dxdy_out.root";
  TFile *sf10_file = new TFile(inputfile_10.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_10 = (TH1D*)sf10_file->Get("hdx_HCAL_p");
  TH1D *hdx_p_10_fid = (TH1D*)sf10_file->Get("hdx_HCAL_p_fid");

  diff[1] = util::fitGaussianAndGetFineMean(hdx_p_10_fid,gen_sig) - util::fitGaussianAndGetFineMean(hdx_p_10,gen_sig);
  sf[1] = 10;

  if( abs(diff[1])>maxDiff )
    maxDiff = abs(diff[1]);

  //set up output for scale = 20 
  string inputfile_20 = inputfiledir + "plots/sf20_dxdy_out.root";
  TFile *sf20_file = new TFile(inputfile_20.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_20 = (TH1D*)sf20_file->Get("hdx_HCAL_p");
  TH1D *hdx_p_20_fid = (TH1D*)sf20_file->Get("hdx_HCAL_p_fid");

  diff[2] = util::fitGaussianAndGetFineMean(hdx_p_20_fid,gen_sig) - util::fitGaussianAndGetFineMean(hdx_p_20,gen_sig);
  sf[2] = 20;

  if( abs(diff[2])>maxDiff )
    maxDiff = abs(diff[2]);

  //set up output for scale = 30 
  string inputfile_30 = inputfiledir + "plots/sf30_dxdy_out.root";
  TFile *sf30_file = new TFile(inputfile_30.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_30 = (TH1D*)sf30_file->Get("hdx_HCAL_p");
  TH1D *hdx_p_30_fid = (TH1D*)sf30_file->Get("hdx_HCAL_p_fid");

  diff[3] = util::fitGaussianAndGetFineMean(hdx_p_30_fid,gen_sig) - util::fitGaussianAndGetFineMean(hdx_p_30,gen_sig);
  sf[3] = 30;

  if( abs(diff[3])>maxDiff )
    maxDiff = abs(diff[3]);

  //set up output for scale = 40 
  string inputfile_40 = inputfiledir + "plots/sf40_dxdy_out.root";
  TFile *sf40_file = new TFile(inputfile_40.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_40 = (TH1D*)sf40_file->Get("hdx_HCAL_p");
  TH1D *hdx_p_40_fid = (TH1D*)sf40_file->Get("hdx_HCAL_p_fid");

  diff[4] = util::fitGaussianAndGetFineMean(hdx_p_40_fid,gen_sig) - util::fitGaussianAndGetFineMean(hdx_p_40,gen_sig);
  sf[4] = 40;

  if( abs(diff[4])>maxDiff )
    maxDiff = abs(diff[4]);

  //set up output for scale = 50 
  string inputfile_50 = inputfiledir + "plots/sf50_dxdy_out.root";
  TFile *sf50_file = new TFile(inputfile_50.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_50 = (TH1D*)sf50_file->Get("hdx_HCAL_p");
  TH1D *hdx_p_50_fid = (TH1D*)sf50_file->Get("hdx_HCAL_p_fid");

  diff[5] = util::fitGaussianAndGetFineMean(hdx_p_50_fid,gen_sig) - util::fitGaussianAndGetFineMean(hdx_p_50,gen_sig);
  sf[5] = 50;

  if( abs(diff[5])>maxDiff )
    maxDiff = abs(diff[5]);

  //set up output for scale = 60 
  string inputfile_60 = inputfiledir + "plots/sf60_dxdy_out.root";
  TFile *sf60_file = new TFile(inputfile_60.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_60 = (TH1D*)sf60_file->Get("hdx_HCAL_p");
  TH1D *hdx_p_60_fid = (TH1D*)sf60_file->Get("hdx_HCAL_p_fid");

  diff[6] = util::fitGaussianAndGetFineMean(hdx_p_60_fid,gen_sig) - util::fitGaussianAndGetFineMean(hdx_p_60,gen_sig);
  sf[6] = 60;

  if( abs(diff[6])>maxDiff )
    maxDiff = abs(diff[6]);

  //set up output for scale = 70 
  string inputfile_70 = inputfiledir + "plots/sf70_dxdy_out.root";
  TFile *sf70_file = new TFile(inputfile_70.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_70 = (TH1D*)sf70_file->Get("hdx_HCAL_p");
  TH1D *hdx_p_70_fid = (TH1D*)sf70_file->Get("hdx_HCAL_p_fid");

  diff[7] = util::fitGaussianAndGetFineMean(hdx_p_70_fid,gen_sig) - util::fitGaussianAndGetFineMean(hdx_p_70,gen_sig);
  sf[7] = 70;

  if( abs(diff[7])>maxDiff )
    maxDiff = abs(diff[7]);

  //set up output for scale = 80 
  string inputfile_80 = inputfiledir + "plots/sf80_dxdy_out.root";
  TFile *sf80_file = new TFile(inputfile_80.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_80 = (TH1D*)sf80_file->Get("hdx_HCAL_p");
  TH1D *hdx_p_80_fid = (TH1D*)sf80_file->Get("hdx_HCAL_p_fid");

  diff[8] = util::fitGaussianAndGetFineMean(hdx_p_80_fid,gen_sig) - util::fitGaussianAndGetFineMean(hdx_p_80,gen_sig);
  sf[8] = 80;

  if( abs(diff[8])>maxDiff )
    maxDiff = abs(diff[8]);

  //set up output for scale = 90 
  string inputfile_90 = inputfiledir + "plots/sf90_dxdy_out.root";
  TFile *sf90_file = new TFile(inputfile_90.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_90 = (TH1D*)sf90_file->Get("hdx_HCAL_p");
  TH1D *hdx_p_90_fid = (TH1D*)sf90_file->Get("hdx_HCAL_p_fid");

  diff[9] = util::fitGaussianAndGetFineMean(hdx_p_90_fid,gen_sig) - util::fitGaussianAndGetFineMean(hdx_p_90,gen_sig);
  sf[9] = 90;

  if( abs(diff[9])>maxDiff )
    maxDiff = abs(diff[9]);

  //set up output for scale = 100 
  string inputfile_100 = inputfiledir + "plots/sf100_dxdy_out.root";
  TFile *sf100_file = new TFile(inputfile_100.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_100 = (TH1D*)sf100_file->Get("hdx_HCAL_p");
  TH1D *hdx_p_100_fid = (TH1D*)sf100_file->Get("hdx_HCAL_p_fid");

  diff[10] = util::fitGaussianAndGetFineMean(hdx_p_100_fid,gen_sig) - util::fitGaussianAndGetFineMean(hdx_p_100,gen_sig);
  sf[10] = 100;

  if( abs(diff[10])>maxDiff )
    maxDiff = abs(diff[10]);

  //set up output for data
  string inputfile_data = inputfiledir + "plots/data_elastic_out_sbs9.root";
  TFile *sfdata_file = new TFile(inputfile_data.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx = (TH1D*)sfdata_file->Get("hdx_cut");
  TH1D *hdx_nofid = (TH1D*)sfdata_file->Get("hdx_cut_nofid");

  double dataDiff = util::fitGaussianAndGetFineMean(hdx,gen_sig) - util::fitGaussianAndGetFineMean(hdx_nofid,gen_sig);
  double dataMean = util::fitGaussianAndGetFineMean(hdx,gen_sig);

  //set up the canvas
  TCanvas *c0 = new TCanvas("c0","data fiducial p peak diff",1600,1200);
  c0->cd();
  //gStyle->SetOptStat(0);

  hdx_nofid->SetTitle("dx sbs9");
  hdx_nofid->SetLineColor(kMagenta);
  hdx_nofid->Draw();
  
  //hdx->SetTitle("dx sbs9");
  hdx->SetLineColor(kBlack);
  hdx->Draw("SAME");
  
  c0->Update();

  // Create and populate the legend
  TLegend *leg = new TLegend(0.7, 0.7, 0.89, 0.89); // Adjust these coordinates as needed
  TString data_difference;
  data_difference.Form("Abs Diff: %0.5f m", dataDiff);  
  TString data_mean;
  data_mean.Form("Mean with Fid Cut: %0.5f m", dataMean);
  leg->AddEntry(hdx, "Elastic and Fiducial Cut", "l");
  leg->AddEntry(hdx_nofid, "Elastic Cut", "l");
  leg->AddEntry((TObject*)0, "", "");
  leg->AddEntry((TObject*)0, data_mean, "");
  leg->AddEntry((TObject*)0, data_difference, "");
  leg->Draw();
  c0->Write();

  // Create TGraph
  TGraph *graph = new TGraph(NFiles, sf, diff);
  graph->SetMarkerStyle(20); // e.g., 20 is a full circle
  graph->SetMarkerSize(1.0); // Size of the marker
  graph->GetXaxis()->SetTitle("Scale Field (%)");
  graph->GetYaxis()->SetTitle("Proton Peak Pos Diff (m)");
  graph->SetTitle("Proton Peak Pos Diff MC (fidcut-nocut) vs Scale Field");

  //set up the canvas
  TCanvas *c1 = new TCanvas("c1","scale field comparisons",1600,1200);
  c1->cd();

  graph->Draw("AP");

  // Create and populate the legend
  TLegend *legend = new TLegend(0.5, 0.6, 0.8, 0.8); // Adjust these coordinates as needed
  TString max_difference;
  max_difference.Form("Max Abs Diff: %0.5f m", maxDiff);
  legend->AddEntry(graph, max_difference, "");
  legend->Draw();
  c1->Write();

  c1->SaveAs("/work/halla/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/quality_plots/sf_getFidVar.png");

  fout->Write();

  cout << "Analysis complete. Outfile located at " << outputfilename << endl;

}
