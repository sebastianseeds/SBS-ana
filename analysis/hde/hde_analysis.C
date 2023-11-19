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

//Total fits using Interpolate with elastic signal histo and 6th order poly fit to bg
TH1D *hW2elas;
Double_t W2total(Double_t *x, Double_t *par){
  Double_t W2 = x[0];
  Double_t sig_scale = par[0];
  Double_t signal = sig_scale * hW2elas->Interpolate(W2);
  return signal + fits::g_p6fit(x,&par[1]);
}

TH1D *hW2elasres;
Double_t W2totalres(Double_t *x, Double_t *par){
  Double_t W2 = x[0];
  Double_t sig_scale = par[0];
  Double_t signal = sig_scale * hW2elasres->Interpolate(W2);
  return signal + fits::g_p6fit(x,&par[1]);
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

    // if(j<40||j>(dnBins-40)){
    //   scaledHist->SetBinContent(j, binContent);
    //   scaledHist->SetBinError(j, binError);
    // }else{
    //   scaledHist->SetBinContent(j, 0.);
    //   scaledHist->SetBinError(j, 0.);
    // }
  }

  return scaledHist;
}

double subtractFunctionAndGetTotal(TH1D* hist, TF1* func, double xMin, double xMax) {
  if (!hist || !func) {
    std::cerr << "Histogram or function is null!" << std::endl;
    return -1.0;
  }

  // Initialize total content after subtraction
  double totalContent = 0.0;

  // Loop over the specified range in the histogram
  int binMin = hist->FindBin(xMin);
  int binMax = hist->FindBin(xMax);

  for (int i = binMin; i <= binMax; ++i) {
    double binCenter = hist->GetBinCenter(i);
    double funcValue = func->Eval(binCenter);
    double binContent = hist->GetBinContent(i);
        
    // Subtract function value from bin content
    double subtractedValue = binContent - funcValue;

    // Add the subtracted value to the total, if it's positive
    if (subtractedValue > 0) {
      totalContent += subtractedValue;
    }
  }

  return totalContent;
}

void hde_analysis( Int_t kine=4, Int_t magset=30, Int_t spot_sigma=2)
{ //main  

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/shde.json");

  std::string rootfile_dir = jmgr->GetValueFromSubKey_str( "rootfile_dir", Form("sbs%d",kine) );
  Double_t W2fitmax = jmgr->GetValueFromSubKey<Double_t>( "W2fitmax", Form("sbs%d",kine) );
  Double_t binfac = jmgr->GetValueFromSubKey<Double_t>( "binfac", Form("sbs%d",kine) );
  Double_t hbinfac = jmgr->GetValueFromSubKey<Double_t>( "hbinfac", Form("sbs%d",kine) );
  Double_t thetapqcut = jmgr->GetValueFromSubKey<Double_t>( "thetapqcut", Form("sbs%d",kine) );
  std::string target = jmgr->GetValueFromKey_str( "target" ); //Must always be LH2 data for proton hde
  Int_t epm = jmgr->GetValueFromSubKey<Int_t>( "epm", Form("sbs%d",kine) );
  Int_t cluster_idx = jmgr->GetValueFromSubKey<Int_t>( "cluster_idx", Form("sbs%d",kine) );
  Int_t pass = jmgr->GetValueFromSubKey<Int_t>( "pass", Form("sbs%d",kine) );
  Double_t ebound_l = jmgr->GetValueFromSubKey<Double_t>( "ebound_l", Form("sbs%d",kine) );
  Double_t ebound_h = jmgr->GetValueFromSubKey<Double_t>( "ebound_h", Form("sbs%d",kine) );
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

  //load histograms

  //string target = "LH2";

  //string infilename = Form("outfiles/hde_dataloop_sbs%d_%s_epm%d_magset%d.root",kine,target.c_str(),epm,magset);
  string infilename = Form("outfiles/hde_dataloop_sbs%d_magset%d_spotsig%d.root",kine,magset,spot_sigma);
  //string outfilename = Form("outfiles/hde_analysis_sbs%d_%s_epm%d_magset%d.root",kine,target.c_str(),epm,magset);
  string outfilename = Form("outfiles/hde_analysis_sbs%d_magset%d_spotsig%d.root",kine,magset,spot_sigma);

  // Open file
  TFile *f1 = TFile::Open(infilename.c_str(), "READ");

  if (!f1 || f1->IsZombie()) {
    std::cerr << "Error opening file: " << infilename << std::endl;
    return;
  }

  // Retrieve the histograms
  TH1D *h1 = (TH1D*)f1->Get("hW2_nocut");
  TH1D *h2 = (TH1D*)f1->Get("hW2_anticut");
  TH1D *h3 = (TH1D*)f1->Get("hW2_allcut");
  TH1D *h4 = (TH1D*)f1->Get("hdx_nocut");
  TH1D *h5 = (TH1D*)f1->Get("hdx_allcut");

  if (!h1 || !h2 || !h3 || !h4 || !h5) {
    std::cerr << "Error retrieving histograms." << std::endl;
    f1->Close();
    delete f1;
    return;
  }


  TH1D *hW2_nocut = (TH1D*)h1->Clone("hW2_nocut");
  TH1D *hW2_anticut = (TH1D*)h2->Clone("hW2_anticut");
  TH1D *hW2_allcut = (TH1D*)h3->Clone("hW2_allcut");
  TH1D *hdx_nocut = (TH1D*)h4->Clone("hdx_nocut");
  TH1D *hdx_allcut = (TH1D*)h5->Clone("hdx_allcut");
  
  // f1->Close();
  // delete f1;

  TFile *fout = new TFile( outfilename.c_str(), "RECREATE" );
  if (!fout || fout->IsZombie()) {
    std::cerr << "Error opening output file: " << outfilename << std::endl;
    return;
  }

  hW2_nocut->Draw();
  hW2_anticut->Draw();
  hW2_allcut->Draw();
  hdx_nocut->Draw();
  hdx_allcut->Draw();

  /////////////////////////////////////////
  //W2 INCLUSIVE ANTICUT INTERPOLATE METHOD

  TH1D *hW2 = (TH1D*)(hW2_nocut->Clone("hW2"));
  TH1D *hW2res = (TH1D*)(hW2_anticut->Clone("hW2res"));
  TH1D *hdx = (TH1D*)(hdx_nocut->Clone("hdx")); 
  TH1D *hdx_2 = (TH1D*)(hdx_nocut->Clone("hdx")); 
  TH1D *hdx_3 = (TH1D*)(hdx_nocut->Clone("hdx")); 
  TH1D *hdx_4 = (TH1D*)(hdx_nocut->Clone("hdx")); 

  Double_t W2fullint = hW2->Integral(0,W2fitmax);
  Double_t W2resfullint = hW2res->Integral(0,W2fitmax);

  //Make canvas for (expected-residuals)/expected
  TCanvas *c1 = new TCanvas("c1",Form("HDE sbs%d mag%d",kine,magset),1200,500);
  gStyle->SetPalette(53);
  c1->Divide(2,1);
  c1->cd(1);

  //Acceptance matching cut only first with "pure elastic" hW2_allcut
  hW2elas = (TH1D*)(hW2_allcut->Clone("hW2elas"));
  TF1 *tfit = new TF1("tfit",W2total,0,W2fitmax,8); //tfit npar = 1+pNfit_npar+1
  TF1 *bg = new TF1("bg",fits::g_p6fit,0.,W2fitmax,7);
  tfit->SetLineColor(kGreen);
  hW2->SetTitle("W^{2}, Acceptance Cut Only");
  hW2->Fit("tfit","RBM");

  Double_t *tpar = tfit->GetParameters();

  //Get fit parameters for bg function and draw identical function on canvas
  bg->SetParameters(&tpar[1]);
  bg->SetLineColor(kRed);
  bg->SetFillColor(kRed);
  bg->SetFillStyle(3005);
  bg->Draw("same");
  hW2elas->SetLineColor(kBlue);
  hW2elas->SetLineWidth(2);
  hW2elas->SetFillColor(kBlue);
  hW2elas->SetFillStyle(3003);
  hW2elas->Draw("same");

  //get error params
  Int_t W2elasb = ebound_l*binfac;
  Int_t W2elase = ebound_h*binfac;
  TFitResultPtr s = hW2->Fit("tfit","S");
  //TMatrixD *bgcov = tfit->GetCovarianceMatrix(); 
  Double_t bgerr = tfit->IntegralError(ebound_l,ebound_h,s->GetParams(),s->GetCovarianceMatrix().GetMatrixArray())*binfac;

  //get expected elastics (elastic divergence from bg begin: ebound_l, end: ebound_h)
  Double_t bgint = bg->Integral(ebound_l,ebound_h)*binfac;
  //Double_t bgerr = bg->IntegralError(ebound_l,ebound_h,bg->GetParameters(),bgcov->GetMatrixArray())*binfac;
  Double_t W2nocutint = hW2->Integral(W2elasb,W2elase);
  Double_t W2nocutinterr = sqrt(W2nocutint);
  Double_t W2elas = W2nocutint - bgint;
  Double_t serr = sqrt( pow(bgerr,2) + pow(W2nocutinterr,2) );

  //Add a legend to the canvas
  auto nocutlegend = new TLegend(0.1,0.6,0.5,0.9);
  //nocutlegend->SetTextSize( 0.03 );
  nocutlegend->AddEntry( hW2elas, "Tight Elastic Cut", "l");
  nocutlegend->AddEntry( bg, "Interpolated Background (scaled)", "l");
  nocutlegend->AddEntry( tfit, "Total fit", "l");
  nocutlegend->AddEntry( (TObject*)0, "", "");
  nocutlegend->AddEntry( (TObject*)0, Form("Number of Elastics Expected: %d",(Int_t)W2elas), "");
  nocutlegend->Draw();

  c1->cd(2);
  
  //HCal anticut to obtain missed elastics
  hW2elasres = (TH1D*)(hW2_allcut->Clone("hW2elasres"));
  TF1 *tfitres = new TF1("tfitres",W2totalres,0,W2fitmax,8);
  tfitres->SetNpx(10000);  // Default is usually 100
  TF1 *bgres = new TF1("bgres",fits::g_p6fit,0.,W2fitmax,7);
  tfitres->SetLineColor(kGreen);
  hW2res->SetTitle("W^{2}, Acceptance Cut and HCal best cluster elastic anticut");
  hW2res->Fit("tfitres","RBM");

  Double_t *tparres = tfitres->GetParameters();

  //Get fit parameters for bg function and draw identical function on canvas
  bgres->SetParameters(&tparres[1]);
  bgres->SetLineColor(kRed);
  bgres->SetFillColor(kRed);
  bgres->SetFillStyle(3005);
  bgres->Draw("same");
  hW2elasres->SetLineColor(kBlue);
  hW2elasres->SetLineWidth(2);
  hW2elasres->SetFillColorAlpha(kBlue,0.35);
  hW2elasres->SetFillStyle(3003);
  hW2elasres->Draw("same");

  //get error
  TFitResultPtr r = hW2res->Fit("tfitres","S");
  //TMatrixD *bgrescov = r->GetCovarianceMatrix();
  Double_t bgreserr = tfitres->IntegralError(ebound_l,ebound_h,r->GetParams(),r->GetCovarianceMatrix().GetMatrixArray())*binfac;

  //get elastics missing from hcal
  Double_t bgresint = bgres->Integral(ebound_l,ebound_h)*binfac; //from fit
  //Double_t bgreserr = bgres->IntegralError(ebound_l,ebound_h,bgres->GetParameters(),bgrescov->GetMatrixArray())*binfac;
  Double_t W2resint = hW2res->Integral(W2elasb,W2elase); //from histo
  Double_t W2reserr = sqrt(W2resint);
  Double_t W2elasres = W2resint - bgresint;
  Double_t reserr = sqrt( pow(bgreserr,2) + pow(W2reserr,2) );
  
  Double_t effres = ( (W2elas-W2elasres) / W2elas )*100.;
  Double_t effreserr = effres*sqrt(pow(serr/W2elas,2)+pow(reserr/W2elasres,2));
  Double_t effhist =  ( (W2fullint-W2resfullint) / W2elas )*100.; //Nonsense - BG will not remain the same between W2nocut and W2residuals
  //Double_t efferr = sqrt( pow(sqrt(effres/100.*(1-effres/100.)/W2elas)*100.,2) + effreserr);
  Double_t efferr = sqrt(effres/100.*(1-effres/100.)/W2elas)*100.;

  //Add a legend to the canvas
  auto reslegend = new TLegend(0.1,0.6,0.5,0.9);
  //reslegend->SetTextSize( 0.03 );
  reslegend->AddEntry( hW2elasres, "Tight Elastic Cut", "l");
  reslegend->AddEntry( bgres, "Interpolated Background (scaled)", "l");
  reslegend->AddEntry( tfitres, "Total fit", "l");
  reslegend->AddEntry( (TObject*)0, "", "");
  reslegend->AddEntry( (TObject*)0, Form("Number of Elastics Detected: %d",(Int_t)(W2elas-W2elasres)), "");
  reslegend->AddEntry( (TObject*)0, Form("HDE via residuals: %0.3f%% +/- %0.3f%%",effres,efferr), "");
  // reslegend->AddEntry( (TObject*)0, "", "");
  // reslegend->AddEntry( (TObject*)0, Form("HDE via histo subtraction: %0.3f%%",effhist), "");
  reslegend->Draw();

  c1->Write();

  ///////////////////////////////
  //DX ELASTIC INTERPOLATE METHOD

  //Make canvas for hcal dx detected / expected 
  TCanvas *c2 = new TCanvas("c2",Form("HDE v2 sbs%d mag%d",kine,magset),1200,500);
  gStyle->SetPalette(53);
  c2->Divide(2,1);
  c2->cd(1);

  //Draw the expected elastic extraction first
  hW2->Draw();
  bg->Draw("same");
  hW2elas->Draw("same");

  auto wlegend = new TLegend(0.1,0.7,0.5,0.9);
  wlegend->AddEntry( hW2elas, "Tight Elastic Cut", "l");
  wlegend->AddEntry( bg, "Interpolated Background (scaled)", "l");
  wlegend->AddEntry( tfit, "Total fit", "l");
  wlegend->AddEntry( (TObject*)0, "", "");
  wlegend->AddEntry( (TObject*)0, Form("Number of Elastics Expected: %d",(Int_t)W2elas), "");
  wlegend->Draw();

  c2->cd(2);

  Double_t lowdxcut = dxmean-3*dxsigma;
  Double_t highdxcut = dxmean+3*dxsigma;
  // Double_t lowdxcut = dx_fitlow;
  // Double_t highdxcut = dx_fithigh;
  Double_t lowdxrange = dxmean-6*dxsigma;
  Double_t highdxrange = dxmean+6*dxsigma;
  Int_t dxelasb = hcalfit_l*hbinfac;
  Int_t dxelase = hcalfit_h*hbinfac;

  hdx_2->SetMinimum(0.0);
  hdx_2->Draw();

  //hdx_2->GetXaxis()->SetRangeUser(lowdxcut,highdxcut);
  //hdx_2->GetXaxis()->SetRangeUser(-0.3,0.3);
  //hdx_2->GetXaxis()->SetRangeUser(-0.5,0.5);
  hdx_2->GetXaxis()->SetRangeUser(dx_fitlow,dx_fithigh);

  //Then draw the dx interpolated signal fit
  hdxelastic = (TH1D*)(hdx_allcut->Clone("hdxelastic"));
  // TF1 *tfitdx = new TF1("tfitdx",dxtotal,hcalfit_l,hcalfit_h,4); //tfit npar = 1+pNfit_npar+1
  // TF1 *bgdx = new TF1("bgdx",fits::g_p2fit,hcalfit_l,hcalfit_h,3); //poly pN fit Nparam = N+1
  TF1 *tfitdx = new TF1("tfitdx",dxtotal,dx_fitlow,dx_fithigh,8); //tfit npar = 1+pNfit_npar+1
  TF1 *bgdx = new TF1("bgdx",fits::g_p6fit,dx_fitlow,dx_fithigh,7); //poly pN fit Nparam = N+1
  tfitdx->SetLineColor(kGreen);
  hdx_2->SetTitle("dx best cluster, Acceptance Cut Only");
  hdx_2->Fit("tfitdx","RBM");

  Double_t *tpardx = tfitdx->GetParameters();

  //Get fit parameters for bg function and draw identical function on canvas
  bgdx->SetParameters(&tpardx[1]);
  bgdx->SetLineColor(kRed);
  bgdx->SetFillColor(kRed);
  bgdx->SetFillStyle(3005);
  bgdx->Draw("same");
  hdxelastic->SetLineColor(kBlue);
  hdxelastic->SetLineWidth(2);
  hdxelastic->SetFillColorAlpha(kBlue,0.35);
  hdxelastic->SetFillStyle(3003);
  hdxelastic->Draw("same");

  // TLine *line1 = new TLine(lowdxcut, hdx_2->GetMinimum(), lowdxcut, hdx_2->GetMaximum());
  // line1->SetLineColor(kRed);
  // line1->SetLineWidth(2);
  // line1->Draw("same");

  // TLine *line2 = new TLine(highdxcut, hdx_2->GetMinimum(), highdxcut, hdx_2->GetMaximum());
  // line2->SetLineColor(kRed);
  // line2->SetLineWidth(2);
  // line2->Draw("same");

  //Get elastics detected in hcal
  Int_t dxbinmin = hdx_2->FindBin(lowdxcut);
  Int_t dxbinmax = hdx_2->FindBin(highdxcut);

  //get error params

  TFitResultPtr t = hdx_2->Fit("tfitdx","S");
  //TMatrixD *bgcov = tfit->GetCovarianceMatrix(); 
  //Double_t dxbgerr = tfitdx->IntegralError(hcalfit_l,hcalfit_h,t->GetParams(),t->GetCovarianceMatrix().GetMatrixArray())*hbinfac;
  Double_t dxbgerr = tfitdx->IntegralError(lowdxcut,highdxcut,t->GetParams(),t->GetCovarianceMatrix().GetMatrixArray())*hbinfac;

  cout << "dxbgerr: " << dxbgerr << endl;

  Double_t bgdxint = bgdx->Integral(lowdxcut,highdxcut)*hbinfac;

  // Double_t int_C = bgdx->Integral(hcalfit_l,hcalfit_h);
  // Double_t int_D = bgdx->Integral(hcalfit_l,hcalfit_h)*hbinfac;
  // Double_t int_E = bgdx->Integral(lowdxcut,highdxcut);
  // Double_t int_F = bgdx->Integral(lowdxcut,highdxcut)*hbinfac;
  // Double_t int_G = hdx_2->Integral();
  // Double_t int_H = hdx_2->Integral(lowdxcut,highdxcut);
  // Double_t int_I = hdx_2->Integral(dxbinmin,dxbinmax);

  Double_t int_direct = 0;
  for( int i=dxbinmin; i<=dxbinmax; i++ ){
    int_direct+=hdx_2->GetBinContent(i);
  }

  //Double_t bgdxint = bgdx->Integral(lowdxcut,highdxcut)*hbinfac;
  Double_t dxnocutint = int_direct;
  //Double_t dxnocutint = hdx_2->Integral(dxbinmin,dxbinmax);
  Double_t dxerr = sqrt(dxnocutint);

  cout << "dxerr: " << dxerr << endl;

  Double_t dxTerr = sqrt( pow(dxbgerr,2) + pow(dxerr,2) );


  Double_t dxelas = dxnocutint - bgdxint;  //hcal elastics dx direct

  Double_t effdx =  ( dxelas / W2elas )*100.; 

  Double_t effdxerr = effdx*sqrt(pow(dxbgerr/dxelas,2)+pow(dxerr/dxTerr,2));

  //Double_t efferr2 = sqrt( pow(sqrt(effdx/100.*(1-effdx/100.)/W2elas)*100.,2) + effdxerr);
  Double_t efferr2 = sqrt(effdx/100.*(1-effdx/100.)/W2elas)*100.;


  // cout << "///////////////////////////////" << endl;
  // cout << "Let's figure this out:" << endl;
  // cout << "bgdx fit values!" << endl;
  // cout << "hbinfac: " << hbinfac << endl;
  // cout << "hcalfit_l: " << hcalfit_l << endl;
  // cout << "hcalfit_h: " << hcalfit_h << endl;
  // cout << "lowdxcut: " << lowdxcut << endl;
  // cout << "highdxcut: " << highdxcut << endl;
  // cout << "dxbinmin: " << dxbinmin << endl;
  // cout << "dxbinmax: " << dxbinmax << endl;
  // cout << "bgdx full integral on full range: " << int_C << endl;
  // cout << "bgdx full integral on full range with binfac: " << int_D << endl;
  // cout << "bgdx full integral on restricted range: " << int_E << endl;
  // cout << "bgdx full integral on restricted range with binfac: " << int_F << endl;
  // cout << endl;
  // cout << "dxnocutint TH1D values!" << endl;
  // cout << "dxelasb: " << dxelasb << endl;
  // cout << "dxelase: " << dxelase << endl;
  // cout << "hdx_2 full integral: " << int_G << endl;
  // cout << "hdx_2 full integral on full range: " << dxnocutint << endl;
  // cout << "hdx_2 full integral on restricted range: " << int_H << endl;
  // cout << "hdx_2 full integral on FindBin restricted range: " << int_I << endl;
  // cout << "hdx_2 full integral on restricted range manual: " << int_direct << endl;
  
  //Add a legend to the canvas
  auto dxlegend = new TLegend(0.1,0.7,0.5,0.9);
  dxlegend->AddEntry( hdxelastic, "Tight Elastic Cut", "l");
  dxlegend->AddEntry( bgdx, Form("Interpolated Background (scaled): %d",(Int_t)bgdxint), "l");
  dxlegend->AddEntry( tfitdx, "Total fit", "l");
  dxlegend->AddEntry( (TObject*)0, "", "");
  dxlegend->AddEntry( (TObject*)0, Form("Total Events on Range: %d",(Int_t)dxnocutint), "");
  dxlegend->AddEntry( (TObject*)0, Form("Number of Elastics Detected: %d",(Int_t)dxelas), "");
  dxlegend->AddEntry( (TObject*)0, Form("HDE via dx: %0.3f%% +/- %0.3f%%",effdx,efferr2), "");
  dxlegend->Draw();

  c2->Write();

  /////////////////////
  //DX SIDEBAND METHOD

  //Make canvas for hcal dx detected / expected 
  TCanvas *c3 = new TCanvas("c3",Form("HDE sideband sbs%d mag%d",kine,magset),1200,500);
  gStyle->SetPalette(53);
  c3->Divide(2,1);
  c3->cd(1);

  //Draw the expected elastic extraction first
  hW2->Draw();
  bg->Draw("same");
  hW2elas->Draw("same");

  auto sbwlegend = new TLegend(0.1,0.7,0.5,0.9);
  sbwlegend->AddEntry( hW2elas, "Tight Elastic Cut", "l");
  sbwlegend->AddEntry( bg, "Interpolated Background (scaled)", "l");
  sbwlegend->AddEntry( tfit, "Total fit", "l");
  sbwlegend->AddEntry( (TObject*)0, "", "");
  sbwlegend->AddEntry( (TObject*)0, Form("Number of Elastics Expected: %d",(Int_t)W2elas), "");
  sbwlegend->Draw();

  c3->cd(2);

  TH1D *hdx_3_scaled = cloneAndCutHistogram(hdx_3,lowdxrange,highdxrange);

  // Perform side-band analysis
  //auto [totalEvents, polyParams] = util::performSideBandAnalysis(hdx_3_scaled, lowdxcut, highdxcut);

  //set reject point range
  SBpol4rej_b = lowdxcut;
  SBpol4rej_e = highdxcut;

  TF1 *bgrpfit = new TF1("bgrpfit",BGfit,lowdxrange,highdxrange,5);
  bgrpfit->SetLineColor(kRed);
  hdx_3_scaled->Fit("bgrpfit","RBM");

  Double_t *bgrppar = bgrpfit->GetParameters();

  TF1 *bgrp = new TF1("bgrp",fits::g_p4fit,lowdxrange,highdxrange,5);
  bgrp->SetParameters(&bgrppar[0]);
  bgrp->SetLineColor(kRed);
  bgrp->SetFillColor(kRed);
  bgrp->SetFillStyle(3005);
  bgrp->Draw("same");

  // Create a TF1 with the pol4 parameters
  // TF1* sbbgdx = new TF1("sbbgdx", "pol4", hdx_3_scaled->GetXaxis()->GetXmin(), hdx_3_scaled->GetXaxis()->GetXmax());
  // sbbgdx->SetParameters(polyParams.data()); // Set all parameters in one call

  hdx_4->SetTitle("x_{exp}-x_{hcal} (m)");
  hdx_4->SetMinimum(0.0);
  hdx_4->SetLineColor(kBlue);
  hdx_4->Draw("same E");

  //get total between bg and dx total
  double totalEvents = subtractFunctionAndGetTotal(hdx_4,bgrp,lowdxcut,highdxcut);
  
  // sbbgdx->SetLineColor(kRed);
  // sbbgdx->Draw("SAME");

  // Retrieve the minimum and maximum y-values of the histogram's y-axis for line placement
  double yMin = hdx_3_scaled->GetMinimum();
  double yMax = hdx_3_scaled->GetMaximum();

  // Create lines at xLow and xHigh
  TLine* lineLow = new TLine(lowdxcut, yMin, lowdxcut, yMax);
  TLine* lineHigh = new TLine(highdxcut, yMin, highdxcut, yMax);

  // Set line color to kMagenta
  lineLow->SetLineColor(kMagenta);
  lineHigh->SetLineColor(kMagenta);

  // Draw the lines on the same canvas
  lineLow->Draw("SAME");
  lineHigh->Draw("SAME");

  TLegend* sblegend = new TLegend(0.1,0.7,0.5,0.9); // Adjust these coordinates as needed
  sblegend->SetBorderSize(1);
  sblegend->SetFillColor(0);
  
  //calculate hde
  double hde_sb = totalEvents/W2elas*100.;

  // Add the histogram, fit function, and lines to the legend
  sblegend->AddEntry( hdx_3_scaled, "Data", "lpe" );
  sblegend->AddEntry( bgrpfit, "Background Fit", "l" );
  sblegend->AddEntry( lineLow, "Sideband Limits", "l" );
  sblegend->AddEntry( (TObject*)0, "", "" );
  sblegend->AddEntry( (TObject*)0, Form("Number of Elastics Detected: %d",(Int_t)totalEvents), "" );
  sblegend->AddEntry( (TObject*)0, Form("Detection Efficiency: %0.2f%%",hde_sb), "" );

  // Draw the legend
  sblegend->Draw();
  c3->Write();

  fout->Write();

  cout << "Analysis complete. Outfile written to " << outfilename << endl;

  // Send time efficiency report to console
  st->Stop();
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
