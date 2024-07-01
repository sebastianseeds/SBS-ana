//sseeds 03.31.23 - test script to use analysis framework to produce hcal detection efficiency output more efficiently
//04.10.23 Update - Added W2 interpolate method for obtaining background fits to both total W2 distribution and W2 with HCal anticut 
//04.11.23 Update - Added back direct dx "detected" yield method for comparison. Fixed thetapq calculation and included in elastic cuts.
//5.20.23 Update - broke from general script and focused this script on extraction of hcal detection efficiency directly from dx after dy cuts normalized by strong earm elastic cuts. 
//5.30.23 Update - NOTE that fits (and resulting yields) are very sensitive to fit ranges. Background shape varies considerably and fit range should reflect 3sigma about the dx elastic peak to avoid overcounting
//10.10.23 Update - split to do data loop first then use hde_analysis.C for fitting and comparisons

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

const Int_t atimeSigFac = 5;
const Double_t dx_fitlow = -0.3;
const Double_t dx_fithigh = 0.3;

double gmean, gsigma, gamp, galpha;

bool verbose = true;

//Total fits using Interpolate with elastic signal histo and 6th order poly fit to bg
TH1D *hW2elas;
Double_t W2total(Double_t *x, Double_t *par){
  Double_t W2 = x[0];
  Double_t sig_scale = par[0];
  Double_t signal = sig_scale * hW2elas->Interpolate(W2);
  return signal + fits::g_p6fit(x,&par[1]);
}

std::vector<Double_t> W2elas_fitpars;
Double_t W2total_sgfit(Double_t *x, Double_t *par){
  Double_t W2 = x[0];
  Double_t amp = W2elas_fitpars[0];
  //Double_t amp = par[0];
  Double_t offset = W2elas_fitpars[1];
  Double_t sigma = W2elas_fitpars[2];
  Double_t alpha = W2elas_fitpars[3];

  return amp*exp( -pow( x[0]-offset,2. )/( 2.*pow(sigma,2.) ) )*( 1+erf( (x[0]-offset)*alpha/sigma*sqrt(2.) ) ) + fits::g_p6fit(x,&par[0]);
}

Double_t sgfit(Double_t *x, Double_t *par){
  Double_t amp = par[0];
  Double_t offset = par[1];
  Double_t sigma = par[2];
  Double_t alpha = par[3];
  //return amp*exp(-0.5*pow((x[0]-offset)/sigma,2.))*0.5*(1+erf(alpha*((x[0]-offset)/sigma)));
  return amp*exp( -pow( x[0]-offset,2. )/( 2.*pow(sigma,2.) ) )*( 1+erf( (x[0]-offset)*alpha/sigma*sqrt(2.) ) );
}

// Double_t W2total_gfit(Double_t *x, Double_t *par){
//   Double_t W2 = x[0];
//   Double_t amp = par[0];
//   Double_t offset = W2elas_fitpars[0];
//   Double_t sigma = W2elas_fitpars[1];

//   return amp*exp(-0.5*pow((x[0]-offset)/sigma,2.)) + fits::g_p6fit(x,&par[1]);
// }

Double_t W2total_gfit(Double_t *x, Double_t *par){
  double gauss = par[0] * TMath::Gaus(x[0], gmean, gsigma);
  double poly = par[1] + par[2] * x[0] + par[3] * x[0] * x[0] + 
    par[4] * x[0] * x[0] * x[0] + par[5] * x[0] * x[0] * x[0] * x[0] + par[6] * x[0] * x[0] * x[0] * x[0] * x[0] + par[7] * x[0] * x[0] * x[0] * x[0] * x[0] * x[0];

  return gauss + poly;
}


TH1D *hW2elasres;
Double_t W2totalres(Double_t *x, Double_t *par){
  Double_t W2 = x[0];
  Double_t sig_scale = par[0];
  Double_t signal = sig_scale * hW2elasres->Interpolate(W2);
  return signal + fits::g_p6fit(x,&par[1]);
}

Double_t W2totalres_sgfit(Double_t *x, Double_t *par){
  Double_t W2 = x[0];
  Double_t amp = par[0];
  Double_t offset = W2elas_fitpars[1];
  Double_t sigma = W2elas_fitpars[2];
  Double_t alpha = W2elas_fitpars[3];

  return amp*exp( -pow( x[0]-offset,2. )/( 2.*pow(sigma,2.) ) )*( 1+erf( (x[0]-offset)*alpha/sigma*sqrt(2.) ) ) + fits::g_p6fit(x,&par[1]);
}

Double_t W2totalres_gfit(Double_t *x, Double_t *par){
  Double_t W2 = x[0];
  Double_t amp = par[0];
  Double_t offset = W2elas_fitpars[0];
  Double_t sigma = W2elas_fitpars[1];

  return amp*exp(-0.5*pow((x[0]-offset)/sigma,2.)) + fits::g_p6fit(x,&par[1]);
}

double W2_signalfit(double *x, double *par) {
    double gauss = par[0] * TMath::Gaus(x[0], gmean, gsigma);
    double poly = par[1] + par[2] * x[0] + par[3] * x[0] * x[0] + 
                  par[4] * x[0] * x[0] * x[0] + par[5] * x[0] * x[0] * x[0] * x[0];

    return gauss + poly;
}

TH1D *hdxelastic;
Double_t dxtotal(Double_t *x, Double_t *par){ //get poly fit to bg with scaled fit to "pure elastics"
  Double_t dx= x[0];
  Double_t sig_scale = par[0];
  Double_t signal = sig_scale * hdxelastic->Interpolate(dx);
  return signal + fits::g_p6fit(x,&par[1]);
}

//Side-band pol4 fit
Double_t SBpol4rej_b; //Central-band fit begin
Double_t SBpol4rej_e; //Central-band fit end

Double_t BGfit(double *x, double *par){

  Double_t yint = par[0];
  Double_t p1 = par[1];
  Double_t p2 = par[2];
  Double_t p3 = par[3];
  Double_t p4 = par[4];

  if(x[0]>SBpol4rej_b && x[0]<SBpol4rej_e) { 
    TF1::RejectPoint();
    return 0;
  }

  return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3)+p4*pow(x[0],4);
}

//to rescale dx x for better sideband analysis
TH1D* cloneAndCutHistogram(TH1D* originalHist, double xMin, double xMax) {
  if (!originalHist) {
    std::cerr << "Original histogram is null!" << std::endl;
    return nullptr;
  }

  // Get bin info
  int binMin = originalHist->FindBin(xMin);
  int binMax = originalHist->FindBin(xMax);
  int dnBins = binMax - binMin + 1; // +1 to include both binMin and binMax

  TH1D *scaledHist = new TH1D("scaledHist", "scaledHist", dnBins, xMin, xMax);

  for (int i = binMin, j = 1; i <= binMax; ++i, ++j) {
    double binContent = originalHist->GetBinContent(i);
    double binError = originalHist->GetBinError(i);
    scaledHist->SetBinContent(j, binContent);
    scaledHist->SetBinError(j, binError);

  }

  return scaledHist;
}

double GetTotalQuadratureError(TH1D* hist) {
  if (!hist) {
    std::cerr << "GetTotalQuadratureError: Histogram pointer is null." << std::endl;
    return 0.0;
  }

  double totalErrorSquared = 0.0;

  // Loop over all bins (excluding underflow and overflow)
  for (int i = 1; i <= hist->GetNbinsX(); ++i) {
    double binError = hist->GetBinError(i);
    totalErrorSquared += binError * binError;
  }

  return sqrt(totalErrorSquared);
}


double GetGrossFitError(TH1D* hist, TF1* func) {
  if (!hist) {
    std::cerr << "GetTotalQuadratureError: Histogram pointer is null." << std::endl;
    return 0.0;
  }

  double totalErrorSquared = 0.0;

  // Loop over all bins (excluding underflow and overflow)
  for (int i = 1; i <= hist->GetNbinsX(); ++i) {
    double binContent = hist->GetBinContent(i);
    double binCenter = hist->GetBinCenter(i);
    double funcValue = func->Eval(binCenter);
    double binError = hist->GetBinError(i);
    double subtractedValue = binContent - funcValue;
    totalErrorSquared += subtractedValue * subtractedValue + binError * binError;
  }

  return sqrt(totalErrorSquared);
}

void subtractFunctionAndGetTotalAndError(TH1D* hist, TF1* func, double xMin, double xMax, double &total, double xRangeLow, double xRangeHigh, double &error) {
  if (!hist || !func) {
    std::cerr << "Histogram or function is null!" << std::endl;
    return -1.0;
  }

  // Initialize total and fit error
  total = 0.0;
  error = 0.0;

  // Loop over the specified range in the histogram
  int binMin = hist->FindBin(xMin);
  int binMax = hist->FindBin(xMax);

  int firstbin = hist->FindBin(xRangeLow);
  int lastbin = hist->FindBin(xRangeHigh);

  // set up placeholders for background and signal
  double error_sig = 0.;
  double error_bg = 0.;
  int sb_bins = 0;
  int sig_bins = 0;

  for(int i = firstbin; i<=lastbin; ++i) {
    double binCenter = hist->GetBinCenter(i);
    double funcValue = func->Eval(binCenter);
    double binContent = hist->GetBinContent(i);
    double binErr = hist->GetBinError(i);
        
    // Subtract function value from bin content
    double subtractedValue = binContent - funcValue;

    if( i>=binMin && i<=binMax ){
      error_sig += binErr;
      sig_bins++;
      // Add the subtracted value to the total, if it's positive
      if (subtractedValue > 0) {
	total += subtractedValue; 
      }
    }else{
      error_bg += subtractedValue;
      sb_bins++;
    }
  }

  double error_bin_bg = error_bg/(double)sb_bins;

  double totalErrorSquared = 0.0;
  for(int i = binMin; i<=binMax; ++i) {
    double binError = hist->GetBinError(i);
    totalErrorSquared += binError * binError + error_bin_bg * error_bin_bg;
  }
  
  error = sqrt(totalErrorSquared);

  return;
}

void hde_dx( Int_t kine=8, Int_t magset=70, Int_t spot_sigma=1, Double_t W2max = 1.5 )
{ //main  

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/shde.json");

  std::string rootfile_dir = jmgr->GetValueFromSubKey_str( "rootfile_dir", Form("sbs%d",kine) );
  Double_t W2fitmax = W2max;
  Double_t binfac = jmgr->GetValueFromSubKey<Double_t>( "binfac", Form("sbs%d",kine) );
  Double_t hbinfac = jmgr->GetValueFromSubKey<Double_t>( "hbinfac", Form("sbs%d",kine) );
  Double_t thetapqcut = jmgr->GetValueFromSubKey<Double_t>( "thetapqcut", Form("sbs%d",kine) );
  std::string target = jmgr->GetValueFromKey_str( "target" ); //Must always be LH2 data for proton hde
  Int_t epm = jmgr->GetValueFromSubKey<Int_t>( "epm", Form("sbs%d",kine) );
  Int_t cluster_idx = jmgr->GetValueFromSubKey<Int_t>( "cluster_idx", Form("sbs%d",kine) );
  Int_t pass = jmgr->GetValueFromSubKey<Int_t>( "pass", Form("sbs%d",kine) );
  Double_t ebound_l = jmgr->GetValueFromSubKey<Double_t>( "ebound_l", Form("sbs%d",kine) );
  Double_t ebound_h = W2max;
  vector<Double_t> coin_profile;
  jmgr->GetVectorFromSubKey<Double_t>("coin_profile",Form("sbs%d",kine),coin_profile);
  Double_t hcal_v_offset = jmgr->GetValueFromSubKey<Double_t>( "hcal_offset", Form("sbs%d",kine) );

  //get tune params
  SBStune tune(kine,magset);

  Double_t dxmean = tune.Getdx0_p();
  Double_t dxsigma = tune.Getdxsig_p();
  Double_t hcalfit_l = econst::hcalposXi_p0-2*econst::hcalblk_w_p0; //lower fit/bin limit for hcal dx plots (m)
  Double_t hcalfit_h = econst::hcalposXf_p0+2*econst::hcalblk_h_p0; //upper fit/bin limit for hcal dx plots (m)
  Double_t harmrange = (hcalfit_h) - (hcalfit_l); //Full range of hcal dx plots (m)


  string infilename = Form("outfiles/hde_dataloop_sbs%d_magset%d_spotsig%d.root",kine,magset,spot_sigma);
  string outfilename = Form("outfiles/hdedx_analysis_sbs%d_magset%d_spotsig%d.root",kine,magset,spot_sigma);

  // Open file
  TFile *f1 = TFile::Open(infilename.c_str(), "READ");

  if (!f1 || f1->IsZombie()) {
    std::cerr << "Error opening file: " << infilename << std::endl;
    return;
  }

  // Retrieve the histograms
  TH1D *h1 = (TH1D*)f1->Get("hW2_nocut");
  TH1D *h2 = (TH1D*)f1->Get("hW2_anticut");
  TH1D *h3 = (TH1D*)f1->Get("hW2_anticut_dym");
  TH1D *h4 = (TH1D*)f1->Get("hW2_allcut");
  TH1D *h5 = (TH1D*)f1->Get("hW2_allcut_dym");
  TH1D *h6 = (TH1D*)f1->Get("hdx_nocut");
  TH1D *h7 = (TH1D*)f1->Get("hdx_allcut");
  TH1D *h8 = (TH1D*)f1->Get("hdx_allcut_dym");

  if (!h1 || !h2 || !h3 || !h4 || !h5 || !h6 || !h7 || !h8) {
    std::cerr << "Error retrieving histograms." << std::endl;
    f1->Close();
    delete f1;
    return;
  }

  TH1D *hW2_nocut = (TH1D*)h1->Clone("hW2_nocut");
  TH1D *hW2_anticut = (TH1D*)h2->Clone("hW2_anticut");
  TH1D *hW2_anticut_dym = (TH1D*)h3->Clone("hW2_anticut_dym");
  TH1D *hW2_allcut = (TH1D*)h4->Clone("hW2_allcut");
  TH1D *hW2_allcut_dym = (TH1D*)h5->Clone("hW2_allcut_dym");
  TH1D *hdx_nocut = (TH1D*)h6->Clone("hdx_nocut");
  TH1D *hdx_allcut = (TH1D*)h7->Clone("hdx_allcut");
  TH1D *hdx_allcut_dym = (TH1D*)h8->Clone("hdx_allcut_dym");

  TFile *fout = new TFile( outfilename.c_str(), "RECREATE" );
  if (!fout || fout->IsZombie()) {
    std::cerr << "Error opening output file: " << outfilename << std::endl;
    return;
  }

  hW2_nocut->Draw();
  hW2_anticut->Draw();
  hW2_anticut_dym->Draw();
  hW2_allcut->Draw();
  hW2_allcut_dym->Draw();
  hdx_nocut->Draw();
  hdx_allcut->Draw();
  hdx_allcut_dym->Draw();

  //Make canvas for (expected-residuals)/expected
  TCanvas *c0 = new TCanvas("c0",Form("SG fit sbs%d mag%d",kine,magset),1200,500);
  gStyle->SetPalette(53);
  c0->cd();

  TH1D *hW2elas_gfit = (TH1D*)(hW2_allcut->Clone("hW2elas2_gfit"));
  W2elas_fitpars = util::fitSkewedGaussianAndGetFineParams(hW2elas_gfit,0.1,0.75,1.0);
  hW2elas_gfit->Draw();
  gamp = W2elas_fitpars[0];
  gmean = W2elas_fitpars[1];
  gsigma = W2elas_fitpars[2];
  galpha = W2elas_fitpars[3];

  c0->Write();

  /////////////////////
  //DX SIDEBAND METHOD

  TH1D *hW2 = (TH1D*)(hW2_nocut->Clone("hW2"));
  TH1D *hW2_2 = (TH1D*)(hW2_nocut->Clone("hW2"));
  TH1D *hW2res = (TH1D*)(hW2_anticut->Clone("hW2res"));
  TH1D *hW2res_2 = (TH1D*)(hW2_anticut->Clone("hW2res"));
  TH1D *hdx = (TH1D*)(hdx_nocut->Clone("hdx")); 
  TH1D *hdx_2 = (TH1D*)(hdx_nocut->Clone("hdx")); 
  TH1D *hdx_3 = (TH1D*)(hdx_nocut->Clone("hdx")); 
  TH1D *hdx_4 = (TH1D*)(hdx_nocut->Clone("hdx")); 
  TH1D *hW2elas_2 = (TH1D*)(hW2_allcut->Clone("hW2elas"));
  
  Double_t W2fullint = hW2->Integral(0,W2fitmax);
  Double_t W2resfullint = hW2res->Integral(0,W2fitmax);
  
  Double_t lowdxcut = dxmean-3*dxsigma;
  Double_t highdxcut = dxmean+3*dxsigma;
  Double_t lowdxrange = dxmean-6*dxsigma;
  Double_t highdxrange = dxmean+6*dxsigma;
  Int_t dxelasb = hcalfit_l*hbinfac;
  Int_t dxelase = hcalfit_h*hbinfac;

  //Make canvas for (expected-residuals)/expected
  TCanvas *c1 = new TCanvas("c1","HDE Total",1200,500);
  //gStyle->SetPalette(53);
  //c1->SetGridx();
  c1->Divide(2,1);
  c1->cd(1);
  
  //hW2elas_2->Draw();

  //Acceptance matching cut only first with "pure elastic" hW2_allcut
  //hW2elas = (TH1D*)(hW2_allcut->Clone("hW2elas"));
  
  hW2->Draw();

  //TH1D *hW2elas_temp = (TH1D*)(hW2_allcut->Clone("hW2elas_temp"));
  //hW2elas = util::MirrorHistogram(hW2elas_temp);

  //hW2elas2->GetXaxis()->SetRangeUser(ebound_l,ebound_h);

  //Do skewed gaussian fit to signal
  TF1 *tfit = new TF1("tfit",W2total_sgfit,0,W2fitmax,7); //tfit npar = 1+pNfit_npar+1
  //TF1 *bg = new TF1("bg",fits::g_p6fit,0.,W2fitmax,7);
  //TF1 *signalfit = new TF1("signalfit",fits::g_sgfit,0.,W2fitmax,4);
  tfit->SetLineColor(kGreen);
  //hW2->GetXaxis()->SetRangeUser(ebound_l,ebound_h);
  hW2->SetTitle("W^{2}, Acceptance Cut Only");
  hW2->Fit("tfit","RBM");

  // Double_t *tpar = tfit->GetParameters();
  
  // //Get fit parameters for bg function and draw identical function on canvas
  // bg->SetParameters(&tpar[0]);
  // bg->SetLineColor(kRed);
  // bg->SetFillColor(kRed);
  // bg->SetFillStyle(3005);
  // bg->Draw("same");
  // signalfit->SetParameters(&W2elas_fitpars[0]);
  // signalfit->SetLineColor(kBlue);
  // signalfit->SetFillColor(kBlue);
  // signalfit->SetFillStyle(3003);
  // signalfit->Draw("same");
  
  // // hW2elas->SetLineColor(kBlue);
  // // hW2elas->SetLineWidth(2);
  // // hW2elas->SetFillColor(kBlue);
  // // hW2elas->SetFillStyle(3003);
  // // hW2elas->Draw("same");
  
  // //get error params
  // Int_t W2elasb = ebound_l*binfac;
  // Int_t W2elase = ebound_h*binfac;
  
  // //get integral error from fits to total W2
  // TFitResultPtr s = hW2->Fit("tfit","S");
  // //TFitResultPtr t = hW2->Fit("bg","S");
  // Double_t tfiterr = tfit->IntegralError(ebound_l,ebound_h,s->GetParams(),s->GetCovarianceMatrix().GetMatrixArray())*binfac;
  
  // Double_t signal_error = GetTotalQuadratureError(hW2elas_2);
  
  // Double_t tot_error = sqrt(pow(tfiterr,2)+pow(signal_error,2));

  // //get expected elastics (elastic divergence from bg begin: ebound_l, end: ebound_h)
  // Double_t bgint = bg->Integral(ebound_l,ebound_h)*binfac;
  // Double_t W2nocutint = hW2->Integral(W2elasb,W2elase);
  // Double_t W2elas = W2nocutint - bgint;
  
  // //Add a legend to the canvas
  // auto nocutlegend = new TLegend(0.1,0.6,0.5,0.9);
  // //nocutlegend->SetTextSize( 0.03 );
  // nocutlegend->AddEntry( hW2elas, "Tight elastic cut (raw)", "l");
  // nocutlegend->AddEntry( bg, "Background (polynomial)", "l");
  // nocutlegend->AddEntry( tfit, "Total interpolated fit", "l");
  // nocutlegend->AddEntry( (TObject*)0, "", "");
  // nocutlegend->AddEntry( (TObject*)0, Form("Number of Elastics Expected: %d",(Int_t)W2elas), "");
  // nocutlegend->Draw();

  // //Make canvas for hcal dx detected / expected 
  // // TCanvas *c3 = new TCanvas("c3",Form("HDE sideband sbs%d mag%d",kine,magset),1200,500);
  // // gStyle->SetPalette(53);
  // // c3->Divide(2,1);
  // // c3->cd(1);

  // // //Draw the expected elastic extraction first
  // // hW2->Draw();
  // // bg->Draw("same");
  // // hW2elas->Draw("same");

  // // auto sbwlegend = new TLegend(0.1,0.7,0.5,0.9);
  // // sbwlegend->AddEntry( hW2elas, "Tight Elastic Cut", "l");
  // // sbwlegend->AddEntry( bg, "Interpolated Background (scaled)", "l");
  // // sbwlegend->AddEntry( tfit, "Total fit", "l");
  // // sbwlegend->AddEntry( (TObject*)0, "", "");
  // // sbwlegend->AddEntry( (TObject*)0, Form("Number of Elastics Expected: %d",(Int_t)W2elas), "");
  // // sbwlegend->Draw();

  // c1->cd(2);

  // TH1D *hdx_3_scaled = cloneAndCutHistogram(hdx_3,lowdxrange,highdxrange);

  // //set reject point range
  // SBpol4rej_b = lowdxcut;
  // SBpol4rej_e = highdxcut;

  // TF1 *bgrpfit = new TF1("bgrpfit",BGfit,lowdxrange,highdxrange,5);
  // bgrpfit->SetLineColor(kRed);
  // hdx_3_scaled->SetTitle("dx, Sideband BG Subtraction");
  // hdx_3_scaled->GetXaxis()->SetTitle("m");
  // hdx_3_scaled->Fit("bgrpfit","RBM");

  // Double_t *bgrppar = bgrpfit->GetParameters();

  // TF1 *bgrp = new TF1("bgrp",fits::g_p4fit,lowdxrange,highdxrange,5);
  // bgrp->SetParameters(&bgrppar[0]);
  // bgrp->SetLineColor(kRed);
  // bgrp->SetFillColor(kRed);
  // bgrp->SetFillStyle(3005);
  // bgrp->Draw("same");

  // hdx_4->SetTitle("x_{exp}-x_{hcal} (m)");
  // hdx_4->SetMinimum(0.0);
  // hdx_4->SetLineColor(kBlue);
  // hdx_4->Draw("same E");

  // //get total between bg and dx total
  // double totalEvents, bgrperr; subtractFunctionAndGetTotalAndError(hdx_4,bgrp,
  // 								   lowdxcut,
  // 								   highdxcut,
  // 								   totalEvents,
  // 								   lowdxrange,
  // 								   highdxrange,
  // 								   bgrperr);

  // //Total the error for quadrature sum
  // Double_t miss_error_dx = bgrperr;
  // //Double_t num_error_dx = sqrt(pow(tot_error,2)+pow(miss_error_dx,2));
  // Double_t num_error_dx = miss_error_dx;

  // // Retrieve the minimum and maximum y-values of the histogram's y-axis for line placement
  // double yMin = hdx_3_scaled->GetMinimum();
  // double yMax = hdx_3_scaled->GetMaximum();

  // // Create lines at xLow and xHigh
  // TLine* lineLow = new TLine(lowdxcut, yMin, lowdxcut, yMax*1.04);
  // TLine* lineHigh = new TLine(highdxcut, yMin, highdxcut, yMax*1.04);

  // // Set line color to kMagenta
  // lineLow->SetLineColor(kMagenta);
  // lineHigh->SetLineColor(kMagenta);

  // // Draw the lines on the same canvas
  // lineLow->Draw("SAME");
  // lineHigh->Draw("SAME");

  // TLegend* sblegend = new TLegend(0.1,0.7,0.5,0.9); // Adjust these coordinates as needed
  // sblegend->SetBorderSize(1);
  // sblegend->SetFillColor(0);
  
  // //calculate hde
  // Double_t hde_sb = totalEvents/W2elas*100.;
  // //Double_t hde_sb_err = hde_sb*sqrt(pow(tot_error/W2elas,2)+pow(num_error_dx/totalEvents,2));
  // Double_t Derr_hist = GetTotalQuadratureError(hW2);
  // Double_t Derr_fit = GetGrossFitError(hW2,tfit);
  // Double_t Derr = sqrt(pow(Derr_hist,2)+pow(Derr_fit,2));

  // Double_t hde_sb_err = hde_sb*sqrt(pow(Derr/W2elas,2)+pow(num_error_dx/totalEvents,2));

  // // Add the histogram, fit function, and lines to the legend
  // sblegend->AddEntry( hdx_3_scaled, "Data", "lpe" );
  // sblegend->AddEntry( bgrpfit, "Background Fit", "l" );
  // sblegend->AddEntry( lineLow, "Sideband Limits", "l" );
  // sblegend->AddEntry( (TObject*)0, "", "" );
  // sblegend->AddEntry( (TObject*)0, Form("Proton spot cut: %d sigma",spot_sigma), "" );
  // sblegend->AddEntry( (TObject*)0, Form("Number of elastics detected: %d",(Int_t)totalEvents), "" );
  // sblegend->AddEntry( (TObject*)0, Form("Detection efficiency: %0.3f%% +/- %0.3f%%",hde_sb,hde_sb_err), "" );

  // // Draw the legend
  // sblegend->Draw();
  // c1->Write();

  cout << "Analysis complete. Outfile written to " << outfilename << endl;

  // Send time efficiency report to console
  st->Stop();
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
