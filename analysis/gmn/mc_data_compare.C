//SSeeds 9.1.23 Script to extract the p/n yield ratio via data/MC comparison. Intent is to loop over all data from a given set/field setting, apply elastic cuts, and build a dx histogram. From this point, read in MC (with RC) distributions for quasi-elastic protons and neutrons (independently), fit these distributions with several distributions to get a best fit, then compose a sum of these fits allowing for a scaling parameter which maps to the total MC yields. Next, do the same with g4sbs inelastic generator and add this to the total fit function (sum) and apply the now three floating parameter fit to the data and check residuals and chi-square. Varying the fit function may yield better results, but the ratio will be extracted. See gmn_calc.C for application of this ratio to gmn extraction.

#include <vector>
#include <iostream>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "TMatrixD.h"
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

const int atimeSigFac = 5;
const int atimediffSigFac = 5;
const double dx_fitlow = -0.3;
const double dx_fithigh = 0.3;

// Set up fit functions for comparisons
// const double hcalfit_l = econst::hcalposXi_mc; //lower fit/bin limit for hcal dx plots (m)
// const double hcalfit_h = econst::hcalposXf_mc; //upper fit/bin limit for hcal dx plots (m)
const double hcalfit_l = -2; //lower fit/bin limit for hcal dx plots (m)
const double hcalfit_h = 0.8; //upper fit/bin limit for hcal dx plots (m)

//Total fits using Interpolate with elastic signal histo and 6th order poly fit to bg
TH1D *hdx_p;
TH1D *hdx_p_clone;
TH1D *hdx_n;
TH1D *hdx_n_clone;
TH1D *hdx_bg;

Double_t dxtotal(Double_t *x, Double_t *par){
  Double_t dx = x[0];
  Double_t proton_scale = par[0];
  Double_t neutron_scale = par[1];
  Double_t bg_scale = par[2];
  Double_t proton = proton_scale * hdx_p->Interpolate(dx);
  Double_t neutron = neutron_scale * hdx_n->Interpolate(dx);
  Double_t bg = bg_scale * hdx_bg->Interpolate(dx);
  return proton + neutron + bg;
}

Double_t dxtotal_bgfit(Double_t *x, Double_t *par){
  Double_t dx = x[0];
  Double_t proton_scale = par[0];
  Double_t neutron_scale = par[1];
  Double_t proton = proton_scale * hdx_p_clone->Interpolate(dx);
  Double_t neutron = neutron_scale * hdx_n_clone->Interpolate(dx);
  return proton + neutron + fits::g_p4fit(x,&par[2]);
}

double IntegrateTH1D(TH1D* hist) {
    double integral = 0.0;
    
    for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
        integral += hist->GetBinContent(bin) * hist->GetBinWidth(bin);
    }
    
    return integral;
}

//Side-band pol4 fit requires beginning and end of rejection region declared here updated later
double SBpol4rej_b = -1.7; //Central-band fit begin
double SBpol4rej_e = 0.7; //Central-band fit end

//Side-band background fit, pol4 with RejectPoint()
double BGfit(double *x, double *par){
  double yint = par[0];
  double p1 = par[1];
  double p2 = par[2];
  double p3 = par[3];
  double p4 = par[4];

  //TF1::RejectPoint() for sideband analysis
  if(x[0]>SBpol4rej_b && x[0]<SBpol4rej_e) { 
    TF1::RejectPoint();
    return 0;
  }

  return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3)+p4*pow(x[0],4);
}

// Function to fit a bootstrap sample and return the scale factors
std::pair<double, double> fitBootstrapSampleBGPoly(TH1D* bootstrapHist, Double_t xMin, Double_t xMax) {
    TF1 *fitFunc = new TF1("totalFit", dxtotal_bgfit, xMin, xMax, 7);
    fitFunc->SetNpx(5000);

    // Perform the fit
    bootstrapHist->Fit(fitFunc, "RBMQ0"); // R for range, Q for quiet mode

    // Extract scale factors
    double scaleFactor1 = fitFunc->GetParameter(0);
    double scaleFactor2 = fitFunc->GetParameter(1);

    //delete fitFunc; // Clean up

    return std::make_pair(scaleFactor1, scaleFactor2);
}

// Function to fit a bootstrap sample and return the scale factors
std::pair<double, double> fitBootstrapSampleBGInter(TH1D* bootstrapHist, Double_t xMin, Double_t xMax) {
    TF1 *fitFunc = new TF1("totalFit", dxtotal, xMin, xMax, 3);
    fitFunc->SetNpx(5000);

    // Perform the fit
    bootstrapHist->Fit(fitFunc, "RBMQ0"); // R for range, Q for quiet mode

    // Extract scale factors
    double scaleFactor1 = fitFunc->GetParameter(0);
    double scaleFactor2 = fitFunc->GetParameter(1);

    //delete fitFunc; // Clean up

    return std::make_pair(scaleFactor1, scaleFactor2);
}

//main (shiftX resulting from fiducial cut: -0.02 SBS9,-0.04 SBS8, 0.0 SBS4)
void mc_data_compare( Int_t kine=4, 
		      Int_t magset=30, 
		      Int_t nBootstrapSamples=10,
		      bool fiducialCut = true,
		      string passkey="abcd",
		      double shiftX=0.1,
		      bool systematic_check=true)
{  
  // Define blinding factor
  Double_t blind_factor = util::generateRandomNumber(passkey);

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  //gStyle->SetOptFit();
  gStyle->SetEndErrorSize(0);

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/gmn.json");
  Int_t epm = jmgr->GetValueFromSubKey<Int_t>( "epm", Form("sbs%d",kine) );
  Int_t cluster_idx = jmgr->GetValueFromSubKey<Int_t>( "cluster_idx", Form("sbs%d",kine) );
  Double_t hcal_v_offset = jmgr->GetValueFromSubKey<Double_t>( "hcal_offset", Form("sbs%d",kine) );

  // set up paths
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fout_path = outdir_path + Form("/gmn_analysis/mc_data_compare_out_sbs%d_mag%d.root",kine,magset);

  TFile *fout = new TFile( fout_path.c_str(), "RECREATE" );

  std::string mc_path = outdir_path + Form("/gmn_analysis/mc_rc_both_out_sbs%d_mag%d.root",kine,magset);
  std::string sys_word = "";
  if(systematic_check)
    sys_word = "_sys";
  std::string data_path = outdir_path + Form("/gmn_analysis/gmn_elastic_fid%s_out_sbs%d_mag%d_clusteridx%d_epm%d.root",sys_word.c_str(),kine,magset,cluster_idx,epm);

  std::string bg_path = outdir_path + "/gmn_analysis/mc_inelastic_out_sbs8_mag0.root";

  // open files and get histograms
  //elastic mc file
  TFile *mc_file = TFile::Open(mc_path.c_str());
  hdx_p = (TH1D*)mc_file->Get("hdx_cut_p"); //key histo for interpolate fits
  hdx_p_clone = (TH1D*)(hdx_p->Clone("hdx_p_clone")); //key histo for interpolate fits
  TH1D *hdx_p_clone2 = (TH1D*)(hdx_p->Clone("hdx_p_clone2"));
  TH1D *hdx_p_clone3 = (TH1D*)(hdx_p->Clone("hdx_p_clone3"));
  hdx_n = (TH1D*)mc_file->Get("hdx_cut_n"); //key histo for interpolate fits
  hdx_n_clone = (TH1D*)(hdx_n->Clone("hdx_n_clone")); //key histo for interpolate fits
  TH1D *hdx_n_clone2 = (TH1D*)(hdx_n->Clone("hdx_n_clone2"));
  TH1D *hdx_n_clone3 = (TH1D*)(hdx_n->Clone("hdx_n_clone3"));
  

  TH1D *hdx_mc = (TH1D*)mc_file->Get("hdx_cut");
  TH1D *hdx_mc_nofid = (TH1D*)mc_file->Get("hdx_cut_nofid");
  TH1D *hdx_mc_failfid = (TH1D*)mc_file->Get("hdx_cut_failfid");
  TH1D *hdx_mc_scaled = (TH1D*)(hdx_mc->Clone("hdx_mc_scaled"));

  //background mc file
  TFile *bg_file = TFile::Open(bg_path.c_str());
  hdx_bg = (TH1D*)bg_file->Get("hdx_cut"); //key histo for bg interpolate fits

  //data file
  TFile *data_file = TFile::Open(data_path.c_str());

  //dx histograms
  TH1D *hdx_data_core;
  if(fiducialCut)
    hdx_data_core = (TH1D*)data_file->Get("hdx_cut");
  else
    hdx_data_core = (TH1D*)data_file->Get("hdx_cut_nofid");
  TH1D *hdx_data_nocut_core = (TH1D*)data_file->Get("hdx_nocut");

  //nocut dx histos
  TH1D *hdx_data_nocut = util::shiftHistogramX(hdx_data_nocut_core,shiftX);
  TH1D *hdx_data_nocut_B = (TH1D*)(hdx_data_nocut->Clone("hdx_data_nocut_B"));

  //cut dx histos
  TH1D *hdx_data = util::shiftHistogramX(hdx_data_core,shiftX);
  TH1D *hdx_data_clone = (TH1D*)(hdx_data->Clone("hdx_data_clone")); 
  TH1D *hdx_data_bsclone = (TH1D*)(hdx_data->Clone("hdx_data_bsclone"));
  TH1D *hdx_data_scaled = (TH1D*)(hdx_data->Clone("hdx_data_scaled"));
  TH1D *hdx_data_clone_B = (TH1D*)(hdx_data->Clone("hdx_data_clone_B"));

  //get mott cross section plot
  TH1D *hMott = (TH1D*)data_file->Get("hMott_cs");

  //get QA plots
  TH2F *hxy = (TH2F*)data_file->Get("hHcalXY");
  TH2F *hdxdy = (TH2F*)data_file->Get("hdxdy_cut");
  TH1D *hdx_data_nofid_core = (TH1D*)data_file->Get("hdx_cut_nofid");
  TH1D *hdx_data_nofid = util::shiftHistogramX(hdx_data_nofid_core,shiftX);
  TH1D *hdx_data_nofid_B = (TH1D*)(hdx_data_nofid->Clone("hdx_data_nofid_B"));
  TH1D *hdx_data_failfid_core = (TH1D*)data_file->Get("hdx_cut_failfid");
  TH1D *hdx_data_failfid = util::shiftHistogramX(hdx_data_failfid_core,shiftX);

  //create scale ratio histogram for error
  TH1D *hErr = new TH1D("hErr","n/p scale ratio",400,0.5,0.9);

  //Check to see if all plots are loaded correctly
  bool allLoaded = true;

  allLoaded &= util::checkTH1D(hdx_p,"hdx_p");
  allLoaded &= util::checkTH1D(hdx_n,"hdx_n");
  allLoaded &= util::checkTH1D(hdx_mc,"hdx_mc");
  allLoaded &= util::checkTH1D(hdx_mc_nofid,"hdx_mc_nofid");
  allLoaded &= util::checkTH1D(hdx_mc_failfid,"hdx_mc_failfid");
  allLoaded &= util::checkTH1D(hdx_data_core,"hdx_data_core");
  allLoaded &= util::checkTH1D(hdx_data_nofid_core,"hdx_data_nofid_core");
  allLoaded &= util::checkTH1D(hdx_data_failfid_core,"hdx_data_failfid_core");
  allLoaded &= util::checkTH1D(hMott,"hMott");
  
  if (!allLoaded) {
    cout << endl << endl << "Failure to load all histograms. Check data and mc files." << endl;
    cout << "Data: " << data_path << endl;
    cout << "MC: " << mc_path << endl;
    return 1;
  }

  //Canvas to view effects of fiducial cuts on data
  TCanvas *cFD = new TCanvas("cFD","Effects of Fiducial Cut, Data",800,600);
  cFD->cd();

  hdx_data_nofid->SetLineColor(kBlack);
  hdx_data_nofid->Draw("hist");

  hdx_data->SetLineColor(kBlue);
  hdx_data->Draw("hist same");

  hdx_data_failfid->SetLineColor(kRed);
  hdx_data_failfid->Draw("hist same");

  //Add a legend to the canvas
  auto legFD = new TLegend(0.11,0.6,0.41,0.89);
  legFD->AddEntry( hdx_data_nofid, "No Fiducial Cut", "l");
  legFD->AddEntry( hdx_data, "Fiducial Cut", "l");
  legFD->AddEntry( hdx_data_failfid, "Failed Fiducial Cut", "l");
  legFD->Draw();

  cFD->Update();
  cFD->Write();

  //Canvas to view effects of fiducial cuts on mc
  TCanvas *cFMC = new TCanvas("cFMC","Effects of Fiducial Cut, MC",800,600);
  cFMC->cd();

  hdx_mc_nofid->SetLineColor(kBlack);
  hdx_mc_nofid->Draw("hist");

  hdx_mc->SetLineColor(kBlue);
  hdx_mc->Draw("hist same");

  hdx_mc_failfid->SetLineColor(kRed);
  hdx_mc_failfid->Draw("hist same");

  //Add a legend to the canvas
  auto legFMC = new TLegend(0.11,0.6,0.41,0.89);
  legFMC->AddEntry( hdx_mc_nofid, "No Fiducial Cut", "l");
  legFMC->AddEntry( hdx_mc, "Fiducial Cut", "l");
  legFMC->AddEntry( hdx_mc_failfid, "Failed Fiducial Cut", "l");
  legFMC->Draw();

  cFMC->Update();
  cFMC->Write();

  //Canvas to view Mott cros section calculated by event
  TCanvas *cQ = new TCanvas("cQ","Mott Cross Section",800,600);
  cQ->cd();

  //Fit and extract weighted Mott cross section
  std::vector<Double_t> mottP = util::fitGaussianAndGetFineParams(hMott, 1e-5);
  hMott->Draw();

  cQ->Update();
  cQ->Write();

  //Canvas to view bootstrap error samples
  TCanvas *c0 = new TCanvas("c0","bootstrap Hist Example",800,600);

  // Set up bootstrap samples for error analysis
  std::vector<TH1D*> bootstrapHists = util::createBootstrapSamples(hdx_data_bsclone, nBootstrapSamples, c0);
  for( int s=0; s<nBootstrapSamples; ++s ){
    std::pair<double, double> scaleFactors = fitBootstrapSampleBGPoly(bootstrapHists[s],hcalfit_l,hcalfit_h);
    double scale_n = scaleFactors.first;
    double scale_p = scaleFactors.second;
    double scale_ratio = blind_factor * scale_n/scale_p;

    hErr->Fill(scale_ratio);
  }

  TCanvas *c0a = new TCanvas("c0a","bootstrap mean hist",800,600);
  c0a->cd();

  //No fits necessary for pure normal distribution (assuming no neglected correlations)
  double ratio_mean = hErr->GetMean();
  double ratio_stdDev = hErr->GetRMS();

  hErr->SetTitle(Form("Bootstrap Ratio Mean Hist (mean: %0.3f RMS: %0.3f)",ratio_mean,ratio_stdDev));
  hErr->Draw();
  c0a->Update();
  c0a->Write();

  // Set up main fit to data for QA
  TF1 *tfitdx = new TF1("tfitdx",dxtotal,hcalfit_l,hcalfit_h,3);
  tfitdx->SetNpx(5000);  // Default is usually 100
  tfitdx->SetLineColor(kGreen);
  TF1 *tfitdx_bgfit = new TF1("tfitdx_bgfit",dxtotal_bgfit,hcalfit_l,hcalfit_h,7); //tfit npar = 1+pNfit_npar+1
  tfitdx_bgfit->SetNpx(5000);  // Default is usually 100

  //Make canvas for direct Data/ MC scale Compare
  TCanvas *c1 = new TCanvas("c1","Data/ MC scale Compare",1200,500);
  c1->cd();

  //TODO: Set maximum dynamically
  double dxpMax = hdx_p->GetMaximum();
  double dxdataMax = hdx_data->GetMaximum();
  
  hdx_p->SetLineColor(kRed);
  hdx_p->SetFillStyle(3005);
  //hdx_p->Draw();

  hdx_data->SetTitle("dx");
  hdx_data->SetLineColor(kBlack);
  hdx_data->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);

  if( dxpMax > dxdataMax ){
    hdx_p->Draw("");
    hdx_data->Draw("same hist");
  }else{
    hdx_data->Draw("hist");
    hdx_p->Draw("same");
  }

  hdx_n->SetLineColor(kBlue);
  hdx_n->SetFillStyle(3003);
  hdx_n->Draw("same");
  
  //Add a legend to the canvas
  auto legdmc = new TLegend(0.11,0.6,0.41,0.89);
  legdmc->AddEntry( hdx_data, "Elastic, Fiducial Cut", "l");
  legdmc->AddEntry( hdx_p, "Elastic proton, MC", "l");
  legdmc->AddEntry( hdx_n, "Elastic neutron, MC", "l");
  legdmc->Draw();

  c1->Update();

  c1->Write();

  hdx_data->Fit("tfitdx_bgfit","RBM0");

  double *tpar = tfitdx_bgfit->GetParameters();

  //Make canvas for interpolate fit with fourth order poly background
  TCanvas *c2 = new TCanvas("c2","Interpolate/4th Order Poly BG",1200,500);
  c2->cd();

  tfitdx_bgfit->SetLineColor(kGreen);
  hdx_data_clone->GetXaxis()->SetRangeUser(hcalfit_l,hcalfit_h);
  hdx_data_clone->Draw("hist");

  hdx_data_clone->SetTitle("dx poly bg fit");
  hdx_data_clone->Fit("tfitdx_bgfit","RBM");
  tfitdx_bgfit->Draw("same");

  double *tpar_clone = tfitdx_bgfit->GetParameters();

  hdx_p_clone2->Scale(tpar_clone[0]);
  hdx_p_clone2->SetLineColor(kOrange);
  hdx_p_clone2->SetLineWidth(2);
  hdx_p_clone->SetFillColorAlpha(kOrange,0.3);
  hdx_p_clone->SetFillStyle(3005);
  hdx_p_clone2->Draw("same");

  hdx_n_clone2->Scale(tpar_clone[1]);
  hdx_n_clone2->SetLineColor(kBlue);
  hdx_n_clone2->SetLineWidth(2);
  hdx_n_clone->SetFillColorAlpha(kBlue,0.3);
  hdx_n_clone->SetFillStyle(3003);
  hdx_n_clone2->Draw("same");

  TF1 *bg = new TF1("bg",fits::g_p4fit,hcalfit_l,hcalfit_h,5);
  bg->SetParameters(&tpar_clone[2]);
  bg->SetLineColor(kRed);
  bg->SetFillColor(kRed);
  bg->SetFillStyle(3004);
  bg->Draw("same");
  
  //Get integrals of scaled proton and neutron histograms
  double Np = IntegrateTH1D(hdx_p_clone2);
  double Nn = IntegrateTH1D(hdx_n_clone2);

  double pn_ratio = blind_factor * Np/Nn;

  double pscale = tpar[0];
  double nscale = tpar[1];

  //double np_par_ratio = blind_factor * nscale/pscale;
  double np_par_ratio = nscale/pscale;

  //Add a legend to the canvas
  auto leg = new TLegend(0.6,0.6,0.89,0.89);
  leg->AddEntry( hdx_data_clone, "Tight Elastic Cut, Data", "l");
  leg->AddEntry( hdx_p_clone2, "MC, proton", "l");
  leg->AddEntry( hdx_n_clone2, "MC, neutron", "l");
  leg->AddEntry( bg, "poly fit to BG", "l");
  leg->AddEntry( tfitdx_bgfit, "Total fit", "l");
  leg->AddEntry( (TObject*)0, "", "");
  //leg->AddEntry( (TObject*)0, Form("proton scale: %0.3f",pscale), "");
  //leg->AddEntry( (TObject*)0, Form("neutron scale: %0.3f",nscale), "");
  leg->AddEntry( (TObject*)0, Form("n/p scale ratio: %0.3f +/- %0.3f",ratio_mean,ratio_stdDev), "");
  leg->Draw();

  c2->Update();

  c2->SaveAs("/work/halla/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/quality_plots/gmn_data_mc.pdf");


  c2->Write();

  //Make canvas for residual analysis
  TCanvas *c3 = new TCanvas("c3","Data/MC Residuals",1200,500);
  c3->cd();

  //add the p and n histograms together for comparison
  TH1D* sumHistogram = new TH1D("sumHistogram", "Sum of Histograms", hdx_p_clone2->GetNbinsX(), hdx_p_clone2->GetXaxis()->GetXmin(), hdx_p_clone2->GetXaxis()->GetXmax());
  sumHistogram->Add(hdx_p_clone2, hdx_n_clone2);

  TH1D *hRes = util::makeResidualHisto("dx",hdx_data_clone,sumHistogram,true,false);

  //hRes
  //hRes->Scale(0.1);
  double dxMaxValue = hdx_data->GetMaximum();
  hRes->GetYaxis()->SetRangeUser(-dxMaxValue/2,dxMaxValue/2);

  hRes->SetLineColor(kRed);
  hRes->SetLineWidth(2);
  hRes->Draw("same");

  c3->Update();
  c3->Write();

  //Make canvas for normalized data and mc comparisons before ratio
  TCanvas *c4 = new TCanvas("c4","Data/MC Direct Compare",1200,500);
  c4->cd();

  //Set the names
  hdx_data_scaled->SetName("Data");
  hdx_mc_scaled->SetName("MC");

  //Normalize both to unity and plot together for reference
  util::plotNormalizedTH1DsOnCanvas(c4,hdx_data_scaled,hdx_mc_scaled,"dx");

  c4->Update();
  c4->Write();

  //Make canvas for normalized data and mc comparisons before ratio
  TCanvas *c5 = new TCanvas("c5","Data-Cut/Data-SB Compare",1200,500);
  c5->Divide(2,1);
  c5->cd(1);

  double nocut_xmin = hdx_data_nocut->GetXaxis()->GetXmin();
  double nocut_xmax = hdx_data_nocut->GetXaxis()->GetXmax();

  TF1 *sbfit = new TF1("sbfit",BGfit,nocut_xmin,nocut_xmax,5);
  hdx_data_nocut->Draw("hist");
  hdx_data_nocut->SetTitle("dx, no cuts, sideband fit");
  hdx_data_nocut->GetXaxis()->SetTitle("m");
  hdx_data_nocut->Fit("sbfit","RBM");
  
  double *sbpar = sbfit->GetParameters();

  TF1 *sb = new TF1("sb",fits::g_p4fit,nocut_xmin,nocut_xmax,5);
  sb->SetParameters(&sbpar[0]);
  sb->SetLineColor(kRed);
  sb->SetFillColor(kRed);
  sb->SetFillStyle(3005);
  sb->Draw("same");

  c5->cd(2);

  hdx_data_nocut_B->SetTitle("dx, data");
  hdx_data_nocut_B->Draw("hist");

  TH1D *hdx_data_sb = util::subtractFunctionFromHistogramWithinRange(hdx_data_nocut_B,sb,hcalfit_l,hcalfit_h);

  hdx_data_sb->SetLineColor(kBlue);
  hdx_data_sb->Draw("hist same");

  hdx_data_clone_B->SetLineColor(kGreen);
  hdx_data_clone_B->Draw("hist same");

  hdx_data_nofid_B->SetLineColor(kGreen-5);
  hdx_data_nofid_B->Draw("hist same");

  auto leg5a = new TLegend(0.6,0.6,0.89,0.89);
  leg5a->AddEntry( hdx_data_nocut_B, "No Cut, Data", "l");
  leg5a->AddEntry( hdx_data_sb, "BG subtract by SB, Data", "l");
  leg5a->AddEntry( hdx_data_clone_B, "Elastic Cut, Data", "l");
  leg5a->AddEntry( hdx_data_nofid_B, "Elastic Cut No Fid, Data", "l");
  leg5a->AddEntry( (TObject*)0, "", "");
  leg5a->Draw();

  c5->Update();
  c5->Write();

  //Make canvas for normalized data and mc comparisons before ratio
  TCanvas *c6 = new TCanvas("c6","Raw fit to data/mc yields",800,500);
  c6->cd();

  double pN_bla = hdx_p_clone3->Integral();
  double nN_bla = hdx_n_clone3->Integral();

  cout << endl << endl << endl << "n:p ratio MC -> " << nN_bla/pN_bla << endl;

  hdx_p_clone3->Add(hdx_n_clone3);
  vector<double> dgsetpar = {1, -0.9, 0.12, 1, 0.0, 0.12, 0.0, 0.0, 0.0, 0.0, 0.0};
  vector<double> dgfitrange = {hcalfit_l,hcalfit_h};
  double yield_ratio_mc;
  hdx_p_clone3->SetTitle("dx, elastic MC");
  util::fitAndCalculateRatio(hdx_p_clone3,c6,yield_ratio_mc,dgsetpar,dgfitrange,true);

  c6->Update();
  c6->Write();

  //Make canvas for normalized data and mc comparisons before ratio
  TCanvas *c7 = new TCanvas("c7","Raw fit to data/mc yields",800,500);
  c7->cd();

  double yield_ratio_data;
  hdx_data_sb->SetTitle("dx earm cuts, sideband bg subtracted");
  util::fitAndCalculateRatio(hdx_data_sb,c7,yield_ratio_data,dgsetpar,dgfitrange,true);

  c7->Update();
  c7->Write();

  fout->Write();
  
  cout << endl << "Analysis complete. Outfile located at " << fout_path << endl;

  // Send time efficiency report to console
  st->Stop();
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
