//sseeds 11.9.23 Script to run over several outputs from sf_dxdy.C, extract proton peak locations, plot means on TGraph, fit TGraph, and get fit params with multiple fits. Also plots all proton peaks on same histogram along with data proton peak in different color. Setting added to include fine adjustment around pass0/1 data proton peak location using fiducial cut variation to set bounds.

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

const bool verbose = true;

const int NFiles = 11;
const double gen_sig = 0.18;

const string field_scale_course[11] = {"0","10","20","30","40","50","60","70","80","90","100"};
const double fs_course[11] = {0,10,20,30,40,50,60,70,80,90,100};
const string field_scale_fine[11] = {"62p2653","63p5395","64p8137","66p0879","67p3621","68p6363","69p9106","70","71p1848","72p4590","73p7332"};
const double fs_fine[11] = {62.2653,63.5395,64.8137,66.0879,67.3621,68.6363,69.9106,70.,71.1848,72.4590,73.7332};

void sf_graph(bool course = false){  

  string field_scale[NFiles];
  string course_string;
  double x_min = -10.;
  double x_max = 10.;
  for( int i=0; i<NFiles; ++i ){
    if(course){
      field_scale[i] = field_scale_course[i];
      course_string = "course";
    }else{
      field_scale[i] = field_scale_fine[i];
      course_string = "fine";
      x_max = -0.5; //Don't select on protons only for fine processing. Explicitely cut out neutron peak.
    }
  }

  //define input file directory
  string inputfiledir = "/volatile/halla/sbs/seeds/scale_field_sbs9/";

  //define output file name
  string outputfilename = inputfiledir + Form("analysis/sf_%sgraph_out.root",course_string.c_str());

  //set up output file
  TFile *fout = new TFile( outputfilename.c_str(), "RECREATE" );

  //set up output for scale = 0 
  string inputfile_0 = inputfiledir + Form("plots/sf%s_dxdy_out.root",field_scale[0].c_str());
  TFile *sf0_file = new TFile(inputfile_0.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_0;
  if(course)
    hdx_p_0 = (TH1D*)sf0_file->Get("hdx_HCAL_p");
  else
    hdx_p_0 = (TH1D*)sf0_file->Get("hdx_HCAL_fid");

  //set up output for scale = 10 
  string inputfile_10 = inputfiledir + Form("plots/sf%s_dxdy_out.root",field_scale[1].c_str());
  TFile *sf10_file = new TFile(inputfile_10.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_10;
  if(course)
    hdx_p_10 = (TH1D*)sf10_file->Get("hdx_HCAL_p");
  else
    hdx_p_10 = (TH1D*)sf10_file->Get("hdx_HCAL_fid");

  //set up output for scale = 20 
  string inputfile_20 = inputfiledir + Form("plots/sf%s_dxdy_out.root",field_scale[2].c_str());
  TFile *sf20_file = new TFile(inputfile_20.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_20;
  if(course)
    hdx_p_20 = (TH1D*)sf20_file->Get("hdx_HCAL_p");
  else
    hdx_p_20 = (TH1D*)sf20_file->Get("hdx_HCAL_fid");

  //set up output for scale = 30 
  string inputfile_30 = inputfiledir + Form("plots/sf%s_dxdy_out.root",field_scale[3].c_str());
  TFile *sf30_file = new TFile(inputfile_30.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_30;
  if(course)
    hdx_p_30 = (TH1D*)sf30_file->Get("hdx_HCAL_p");
  else
    hdx_p_30 = (TH1D*)sf30_file->Get("hdx_HCAL_fid");

  //set up output for scale = 40 
  string inputfile_40 = inputfiledir + Form("plots/sf%s_dxdy_out.root",field_scale[4].c_str());
  TFile *sf40_file = new TFile(inputfile_40.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_40;
  if(course)
    hdx_p_40 = (TH1D*)sf40_file->Get("hdx_HCAL_p");
  else
    hdx_p_40 = (TH1D*)sf40_file->Get("hdx_HCAL_fid");

  //set up output for scale = 50 
  string inputfile_50 = inputfiledir + Form("plots/sf%s_dxdy_out.root",field_scale[5].c_str());
  TFile *sf50_file = new TFile(inputfile_50.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_50;
  if(course)
    hdx_p_50 = (TH1D*)sf50_file->Get("hdx_HCAL_p");
  else
    hdx_p_50 = (TH1D*)sf50_file->Get("hdx_HCAL_fid");

  //set up output for scale = 60 
  string inputfile_60 = inputfiledir + Form("plots/sf%s_dxdy_out.root",field_scale[6].c_str());
  TFile *sf60_file = new TFile(inputfile_60.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_60;
  if(course)
    hdx_p_60 = (TH1D*)sf60_file->Get("hdx_HCAL_p");
  else
    hdx_p_60 = (TH1D*)sf60_file->Get("hdx_HCAL_fid");

  //set up output for scale = 70 
  string inputfile_70 = inputfiledir + Form("plots/sf%s_dxdy_out.root",field_scale[7].c_str());
  TFile *sf70_file = new TFile(inputfile_70.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_70;
  if(course)
    hdx_p_70 = (TH1D*)sf70_file->Get("hdx_HCAL_p");
  else
    hdx_p_70 = (TH1D*)sf70_file->Get("hdx_HCAL_fid");

  //set up output for scale = 80 
  string inputfile_80 = inputfiledir + Form("plots/sf%s_dxdy_out.root",field_scale[8].c_str());
  TFile *sf80_file = new TFile(inputfile_80.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_80;
  if(course)
    hdx_p_80 = (TH1D*)sf80_file->Get("hdx_HCAL_p");
  else
    hdx_p_80 = (TH1D*)sf80_file->Get("hdx_HCAL_fid");

  //set up output for scale = 90 
  string inputfile_90 = inputfiledir + Form("plots/sf%s_dxdy_out.root",field_scale[9].c_str());
  TFile *sf90_file = new TFile(inputfile_90.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_90;
  if(course)
    hdx_p_90 = (TH1D*)sf90_file->Get("hdx_HCAL_p");
  else
    hdx_p_90 = (TH1D*)sf90_file->Get("hdx_HCAL_fid");

  //set up output for scale = 100 
  string inputfile_100 = inputfiledir + Form("plots/sf%s_dxdy_out.root",field_scale[10].c_str());
  TFile *sf100_file = new TFile(inputfile_100.c_str(), "READ"); // Open each ROOT file in read mode
  TH1D *hdx_p_100;
  if(course)
    hdx_p_100 = (TH1D*)sf100_file->Get("hdx_HCAL_p");
  else
    hdx_p_100 = (TH1D*)sf100_file->Get("hdx_HCAL_fid");

  double gy[NFiles];
  double gx[NFiles];
  double gdy[NFiles];
  double gdx[NFiles] = {0}; //no error in x (field scale set by MC)

  //set up the holding canvas
  TCanvas *c0 = new TCanvas("c0","all peaks",1600,1200);
  c0->cd();
  gStyle->SetOptStat(0);

  //fit 0
  hdx_p_0->SetLineColor(kRed);
  hdx_p_0->GetYaxis()->SetRangeUser(0,1000);
  hdx_p_0->GetXaxis()->SetTitle("dx (m)");
  if(!course)
    hdx_p_0->GetXaxis()->SetRangeUser(fs_fine[0]-gen_sig,fs_fine[NFiles-1]+gen_sig);
  hdx_p_0->SetTitle("dx, scale field comparisons");
  hdx_p_0->Draw();

  vector<double> par0 = util::fitGaussianAndGetFineParams(hdx_p_0,gen_sig,x_min,x_max);
  gy[0]=par0[1];
  gdy[0]=par0[2];
  if(course)
    gx[0]=fs_course[0];
  else
    gx[0]=fs_fine[0];
  cout << "fit 0: p0_x=" << par0[0] << " p1_x=" << par0[1] << " p2_x=" << par0[2] << endl;

  //fit 1
  hdx_p_10->SetLineColor(kOrange);
  hdx_p_10->Draw("same");

  vector<double> par1 = util::fitGaussianAndGetFineParams(hdx_p_10,gen_sig,x_min,x_max);
  gy[1]=par1[1];
  gdy[1]=par1[2];
  if(course)
    gx[1]=fs_course[1];
  else
    gx[1]=fs_fine[1];
  cout << "fit 1: p0_x=" << par1[0] << " p1_x=" << par1[1] << " p2_x=" << par1[2] << endl;

  //fit 2
  hdx_p_20->SetLineColor(kGreen);
  hdx_p_20->Draw("same");

  vector<double> par2 = util::fitGaussianAndGetFineParams(hdx_p_20,gen_sig,x_min,x_max);
  gy[2]=par2[1];
  gdy[2]=par2[2];
  if(course)
    gx[2]=fs_course[2];
  else
    gx[2]=fs_fine[2];
  cout << "fit 2: p0_x=" << par2[0] << " p1_x=" << par2[1] << " p2_x=" << par2[2] << endl;

  //fit 3
  hdx_p_30->SetLineColor(kMagenta);
  hdx_p_30->Draw("same");

  vector<double> par3 = util::fitGaussianAndGetFineParams(hdx_p_30,gen_sig,x_min,x_max);
  gy[3]=par3[1];
  gdy[3]=par3[2];
  if(course)
    gx[3]=fs_course[3];
  else
    gx[3]=fs_fine[3];
  cout << "fit 3: p0_x=" << par3[0] << " p1_x=" << par3[1] << " p2_x=" << par3[2] << endl;

  //fit 4
  hdx_p_40->SetLineColor(kBlue);
  hdx_p_40->Draw("same");

  vector<double> par4 = util::fitGaussianAndGetFineParams(hdx_p_40,gen_sig,x_min,x_max);
  gy[4]=par4[1];
  gdy[4]=par4[2];
  if(course)
    gx[4]=fs_course[4];
  else
    gx[4]=fs_fine[4];
  cout << "fit 4: p0_x=" << par4[0] << " p1_x=" << par4[1] << " p2_x=" << par4[2] << endl;

  //fit 5
  hdx_p_50->SetLineColor(kCyan);
  hdx_p_50->Draw("same");

  vector<double> par5 = util::fitGaussianAndGetFineParams(hdx_p_50,gen_sig,x_min,x_max);
  gy[5]=par5[1];
  gdy[5]=par5[2];
  if(course)
    gx[5]=fs_course[5];
  else
    gx[5]=fs_fine[5];
  cout << "fit 5: p0_x=" << par5[0] << " p1_x=" << par5[1] << " p2_x=" << par5[2] << endl;

  //fit 6
  hdx_p_60->SetLineColor(kPink);
  hdx_p_60->Draw("same");

  vector<double> par6 = util::fitGaussianAndGetFineParams(hdx_p_60,gen_sig,x_min,x_max);
  gy[6]=par6[1];
  gdy[6]=par6[2];
  if(course)
    gx[6]=fs_course[6];
  else
    gx[6]=fs_fine[6];
  cout << "fit 6: p0_x=" << par6[0] << " p1_x=" << par6[1] << " p2_x=" << par6[2] << endl;

  //fit 7
  hdx_p_70->SetLineColor(kViolet);
  hdx_p_70->Draw("same");

  vector<double> par7 = util::fitGaussianAndGetFineParams(hdx_p_70,gen_sig,x_min,x_max);
  gy[7]=par7[1];
  gdy[7]=par7[2];
  if(course)
    gx[7]=fs_course[7];
  else
    gx[7]=fs_fine[7];
  cout << "fit 7: p0_x=" << par7[0] << " p1_x=" << par7[1] << " p2_x=" << par7[2] << endl;

  //fit 8
  hdx_p_80->SetLineColor(kBlack);
  hdx_p_80->Draw("same");

  vector<double> par8 = util::fitGaussianAndGetFineParams(hdx_p_80,gen_sig,x_min,x_max);
  gy[8]=par8[1];
  gdy[8]=par8[2];
  if(course)
    gx[8]=fs_course[8];
  else
    gx[8]=fs_fine[8];
  cout << "fit 8: p0_x=" << par8[0] << " p1_x=" << par8[1] << " p2_x=" << par8[2] << endl;

  //fit 9
  hdx_p_90->SetLineColor(kAzure);
  hdx_p_90->Draw("same");

  vector<double> par9 = util::fitGaussianAndGetFineParams(hdx_p_90,gen_sig,x_min,x_max);
  gy[9]=par9[1];
  gdy[9]=par9[2];
  if(course)
    gx[9]=fs_course[9];
  else
    gx[9]=fs_fine[9];
  cout << "fit 9: p0_x=" << par9[0] << " p1_x=" << par9[1] << " p2_x=" << par9[2] << endl;

  //fit 10
  hdx_p_100->SetLineColor(kGray);
  hdx_p_100->Draw("same");

  vector<double> par10 = util::fitGaussianAndGetFineParams(hdx_p_100,gen_sig,x_min,x_max);
  gy[10]=par10[1];
  gdy[10]=par10[2];
  if(course)
    gx[10]=fs_course[10];
  else
    gx[10]=fs_fine[10];
  cout << "fit 10: p0_x=" << par10[0] << " p1_x=" << par10[1] << " p2_x=" << par10[2] << endl;

  // Create and populate the legend
  TLegend *legend0 = new TLegend(0.59, 0.7, 0.89, 0.9); // Adjust these coordinates as needed
  legend0->AddEntry(hdx_p_0, Form("scale field = %s %%",field_scale[0].c_str()), "l");
  legend0->AddEntry(hdx_p_10, Form("scale field = %s %%",field_scale[1].c_str()), "l");
  legend0->AddEntry(hdx_p_20, Form("scale field = %s %%",field_scale[2].c_str()), "l");
  legend0->AddEntry(hdx_p_30, Form("scale field = %s %%",field_scale[3].c_str()), "l");
  legend0->AddEntry(hdx_p_40, Form("scale field = %s %%",field_scale[4].c_str()), "l");
  legend0->AddEntry(hdx_p_50, Form("scale field = %s %%",field_scale[5].c_str()), "l");
  legend0->AddEntry(hdx_p_60, Form("scale field = %s %%",field_scale[6].c_str()), "l");
  legend0->AddEntry(hdx_p_70, Form("scale field = %s %%",field_scale[7].c_str()), "l");
  legend0->AddEntry(hdx_p_80, Form("scale field = %s %%",field_scale[8].c_str()), "l");
  legend0->AddEntry(hdx_p_90, Form("scale field = %s %%",field_scale[9].c_str()), "l");
  legend0->AddEntry(hdx_p_100, Form("scale field = %s %%",field_scale[10].c_str()), "l");
  legend0->Draw();

  c0->Write();

  c0->SaveAs(Form("/work/halla/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/quality_plots/sf_%soverlay.png",course_string.c_str()));

  // Create TGraph
  gStyle->SetEndErrorSize(0);

  TGraph *graph = new TGraphErrors(NFiles, gx, gy, gdx, gdy);
  graph->SetMarkerStyle(20); // e.g., 20 is a full circle
  graph->SetMarkerSize(1.0); // Size of the marker
  graph->SetTitle("Proton Peak Means vs Scale Field");

  // Fit with a straight line
  TF1 *fitFunc = new TF1("fitFunc", "pol1"); // pol1 for linear fit
  graph->Fit(fitFunc, "Q"); // "Q" for quiet mode (no print)

  // Print slope and y-intercept
  double slope = fitFunc->GetParameter(1);
  double yIntercept = fitFunc->GetParameter(0);

  //set up the canvas
  TCanvas *c1 = new TCanvas("c1","scale field comparisons",1600,1200);
  c1->cd();

  graph->Draw("AP");

  // Create and populate the legend
  TLegend *legend = new TLegend(0.5, 0.7, 0.8, 0.9); // Adjust these coordinates as needed
  TString slopeStr;
  slopeStr.Form("Slope: %0.5f", slope);
  TString interceptStr;
  interceptStr.Form("Y-Intercept: %0.5f", yIntercept);
  legend->AddEntry(graph, slopeStr, "l");
  legend->AddEntry(graph, interceptStr, "l");
  legend->Draw();

  std::cout << "Slope: " << slope << ", Y-intercept: " << yIntercept << std::endl;

  c1->Write();

  c1->SaveAs(Form("/work/halla/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/quality_plots/sf_%sgraph.png",course_string.c_str()));

  fout->Write();

  cout << "Analysis complete. Outfile located at " << outputfilename << endl;

  if(verbose){
    cout << endl << endl << "Proton Peak Mean Locations" << endl;
    for( int i=0; i<NFiles; ++i )
      cout << gy[i] << endl;

    cout << endl << endl <<"Proton Peak Std Dev" << endl;
    for( int i=0; i<NFiles; ++i )
      cout << gdy[i] << endl;
  }

}
