//sseeds 03.31.23 - test script to use analysis framework to produce hcal detection efficiency output more efficiently
//04.10.23 Update - Added W2 interpolate method for obtaining background fits to both total W2 distribution and W2 with HCal anticut 
//04.11.23 Update - Added back direct dx "detected" yield method for comparison. Fixed thetapq calculation and included in elastic cuts.
//5.20.23 Update - broke from general script and focused this script on extraction of hcal detection efficiency directly from dx after dy cuts normalized by strong earm elastic cuts. 
//5.30.23 Update - NOTE that fits (and resulting yields) are very sensitive to fit ranges. Background shape varies considerably and fit range should reflect 3sigma about the dx elastic peak to avoid overcounting
//10.10.23 Update - split to do data loop first then use hde_analysis.C for fitting and comparisons
//1.11.24 Update - rewrite to simplify analysis. Considering only W2 skewed gaussian fits for denominator, thus totalfit anticut and totalfit dx sideband methods

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

//Set up necessary local fit functions
//Several fits used are global and located in fits.C

//W2 FITS
std::vector<double> W2elas_sgfitpars;

//Skewed gaussian + sixth order polynomial total fit for W2. Function only scales the skewed gaussian whose shape is determined with tight elastic cuts on e-arm and h-arm
double W2total_sgfit(double *x, double *par){
  double W2 = x[0];
  double amp = par[0];
  double offset = W2elas_sgfitpars[1];
  double sigma = W2elas_sgfitpars[2];
  double alpha = W2elas_sgfitpars[3];

  return amp*exp( -pow( x[0]-offset,2. )/( 2.*pow(sigma,2.) ) )*( 1+erf( (x[0]-offset)*alpha/sigma*sqrt(2.) ) ) + fits::g_p6fit(x,&par[1]);
}

//Duplicate fit here for anticut histogram
double W2totalanticut_sgfit(double *x, double *par){
  double W2 = x[0];
  double amp = par[0];
  double offset = W2elas_sgfitpars[1];
  double sigma = W2elas_sgfitpars[2];
  double alpha = W2elas_sgfitpars[3];

  return amp*exp( -pow( x[0]-offset,2. )/( 2.*pow(sigma,2.) ) )*( 1+erf( (x[0]-offset)*alpha/sigma*sqrt(2.) ) ) + fits::g_p6fit(x,&par[1]);
}

//DX FITS
//Side-band pol4 fit requires beginning and end of rejection region declared here updated later
double SBpol4rej_b; //Central-band fit begin
double SBpol4rej_e; //Central-band fit end

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

//set up semi-static spot cuts for elastic shape (exp) and allowed detection events (det)
double det_spot_sigma = 3.;
double exp_spot_sigma = 1.;
double W2_fit_width_sigma = 2.; //Set to get fine param fit on W2 elastic signal
//double dx_fit_width_sigma = 4.1; //Set to minimize variance arising from BB shadow
//double dx_sb_width_sigma = 2.3; //Tuned to best chi-squared sb fit to data with fourth order poly
//double dx_fit_width_sigma = 4.1; //SBS4, SBS8 (anticut)
//double dx_sb_width_sigma = 2.3; //SBS4, SBS8 (anticut)
double dx_fit_width_sigma = 3.2; //SBS8 (sideband)
double dx_sb_width_sigma = 2.3; //SBS8 (sideband)


void hde_ana( Int_t kine=4, Int_t magset=30, bool scoopclus = false )
{ //main  

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/shde.json");
  double binfac = jmgr->GetValueFromSubKey<double>( "binfac", Form("sbs%d",kine) );
  double hbinfac = jmgr->GetValueFromSubKey<double>( "hbinfac", Form("sbs%d",kine) );
  double ebound_l = jmgr->GetValueFromSubKey<double>( "ebound_l", Form("sbs%d",kine) );
  double ebound_h = jmgr->GetValueFromSubKey<double>( "ebound_h", Form("sbs%d",kine) );
 
  //get tune params
  SBStune tune(kine,magset);
  double W2mean = tune.GetW2mean();
  double W2sig = tune.GetW2sig();
  double dxmean = tune.Getdx0_p();
  double dxsigma = tune.Getdxsig_p();
  double hcalfit_l = econst::hcalposXi_p0-2*econst::hcalblk_w_p0; //lower fit/bin limit for hcal dx plots (m)
  double hcalfit_h = econst::hcalposXf_p0+2*econst::hcalblk_h_p0; //upper fit/bin limit for hcal dx plots (m)
  double harmrange = (hcalfit_h) - (hcalfit_l); //Full range of hcal dx plots (m)

  // Get file names
  string infilename = Form("outfiles/hde_dataloop_sbs%d_magset%d_dspotsig%0.1f_espotsig%0.1f.root",kine,magset,det_spot_sigma,exp_spot_sigma);
  string outfilename = Form("outfiles/hde_analysis_sbs%d_magset%d_dspotsig%0.1f_espotsig%0.1f.root",kine,magset,det_spot_sigma,exp_spot_sigma);

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
  TH1D *h9 = (TH1D*)f1->Get("hW2_gcut");
  TH1D *h10 = (TH1D*)f1->Get("hdx_gcut");
  TH1D *h11 = (TH1D*)f1->Get("hdx_nocut_scoop");
  TH1D *h12 = (TH1D*)f1->Get("hdx_gcut_scoop");
  TH1D *h13 = (TH1D*)f1->Get("hW2_allcut_scoop");
  TH1D *h14 = (TH1D*)f1->Get("hW2_anticut_scoop");

  if (!h1 || !h2 || !h3 || !h4 || !h5 || !h6 || !h7 || !h8 || !h9) {
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
  TH1D *hW2_gcut = (TH1D*)h9->Clone("hW2_gcut");
  TH1D *hdx_gcut = (TH1D*)h10->Clone("hdx_gcut");
  TH1D *hdx_nocut_scoop = (TH1D*)h11->Clone("hdx_nocut_scoop");
  TH1D *hdx_gcut_scoop = (TH1D*)h12->Clone("hdx_gcut_scoop");
  TH1D *hW2_allcut_scoop = (TH1D*)h13->Clone("hW2_allcut_scoop");
  TH1D *hW2_anticut_scoop = (TH1D*)h14->Clone("hW2_anticut_scoop");

  TFile *fout = new TFile( outfilename.c_str(), "RECREATE" );
  if (!fout || fout->IsZombie()) {
    std::cerr << "Error opening output file: " << outfilename << std::endl;
    return;
  }

  /////////////////////////////////
  //W2 INCLUSIVE ANTICUT FIT METHOD

  //Make canvas to check "pure" elastic fit
  TCanvas *c0 = new TCanvas("c0",Form("Elastic Sample sbs%d mag%d",kine,magset),600,500);
  gStyle->SetPalette(55);
  c0->cd();

  //Clone and scale the W2 histograms
  TH1D *hW2pure_scaled;
  if(scoopclus)
    hW2pure_scaled = util::cloneAndCutHistogram(hW2_allcut_scoop,ebound_l,ebound_h);
  else
    hW2pure_scaled = util::cloneAndCutHistogram(hW2_allcut,ebound_l,ebound_h);

  TF1 *pure = new TF1("pure",fits::g_sgfit,ebound_l,ebound_h,4);

  std::vector<double> W2pure_sgfitpars; util::fitSkewedGaussianAndGetFineParams(hW2pure_scaled,W2_fit_width_sigma,ebound_l,ebound_h,W2pure_sgfitpars,0.1);
  pure->SetParameters(&W2pure_sgfitpars[0]);

  hW2pure_scaled->Draw();
  pure->Draw("same");

  //Add a legend to the canvas
  auto purelegend = new TLegend(0.1,0.6,0.5,0.9);
  purelegend->AddEntry( (TObject*)0, Form("Alpha: %0.4f",W2pure_sgfitpars[3]), "");
  purelegend->Draw();

  //Make canvas for (expected-residuals)/expected
  TCanvas *c1 = new TCanvas("c1",Form("HDE sbs%d mag%d",kine,magset),1200,500);
  gStyle->SetPalette(55);
  c1->SetGridx();
  c1->Divide(2,1);
  c1->cd(1);

  //Clone and scale the W2 histograms
  TH1D *hW2_scaled = util::cloneAndCutHistogram(hW2_gcut,ebound_l,ebound_h); //should use globalcut here to compare only the effects of the hadron arm
  
  TH1D *hW2elas_scaled;
  if(scoopclus)
    hW2elas_scaled = util::cloneAndCutHistogram(hW2_allcut_scoop,ebound_l,ebound_h);
  else
    hW2elas_scaled = util::cloneAndCutHistogram(hW2_allcut,ebound_l,ebound_h);

  //Fit the "pure" elastic W2 distribution with a skewed gaussian
  util::fitSkewedGaussianAndGetFineParams(hW2elas_scaled,W2_fit_width_sigma,ebound_l,ebound_h,W2elas_sgfitpars,0.1);
  
  //Set up fits and fit
  TF1 *tfit = new TF1("tfit",W2total_sgfit,ebound_l,ebound_h,8); //tfit npar = 1+pNfit_npar+1
  tfit->SetNpx(10000);  // Default is usually 100
  TF1 *bg = new TF1("bg",fits::g_p6fit,ebound_l,ebound_h,7);
  bg->SetNpx(10000);  // Default is usually 100
  TF1 *elastics = new TF1("elastics",fits::g_sgfit,ebound_l,ebound_h,4);
  elastics->SetNpx(10000);  // Default is usually 100
  tfit->SetLineColor(kGreen);
  hW2_scaled->SetTitle("W^{2}, Acceptance cut only");
  hW2_scaled->Fit("tfit","RBM");

  double *tpar = tfit->GetParameters();
  
  //Get fit parameters for bg function and draw identical function on canvas
  bg->SetParameters(&tpar[1]);
  bg->SetLineColor(kRed);
  bg->SetFillColor(kRed);
  bg->SetFillStyle(3005);
  bg->Draw("same");
  elastics->SetParameters(&W2pure_sgfitpars[0]);
  elastics->SetLineColor(kBlue);
  elastics->SetLineWidth(2);
  elastics->SetFillColor(kBlue);
  elastics->SetFillStyle(3003);
  elastics->Draw("same");

  //get expected elastics
  double totaldiff = 0.;
  double totalerr = 0.;
  util::subtractFunctionAndGetTotalAndError(hW2_scaled,bg,ebound_l,ebound_h,totaldiff,ebound_l,ebound_h,totalerr,false);

  //Add a legend to the canvas
  auto nocutlegend = new TLegend(0.1,0.6,0.5,0.9);
  nocutlegend->AddEntry( elastics, "Tight elastic cut (sk gaus)", "l");
  nocutlegend->AddEntry( bg, "Background (poly)", "l");
  nocutlegend->AddEntry( tfit, "Total sum fit", "l");
  nocutlegend->AddEntry( (TObject*)0, "", "");
  nocutlegend->AddEntry( (TObject*)0, Form("Exp Elastics: %d",(Int_t)totaldiff), "");
  nocutlegend->Draw();
  
  c1->cd(2);
  
  //HCal anticut to obtain missed elastics
  
  //Clone some additional histograms to total anticut events
  TH1D *hW2anticut_scaled;
  if(scoopclus)
    hW2anticut_scaled = util::cloneAndCutHistogram(hW2_anticut_scoop,ebound_l,ebound_h);
  else
    hW2anticut_scaled = util::cloneAndCutHistogram(hW2_anticut,ebound_l,ebound_h);

  //Get total number of events in W2 anticut histogram
  double W2anticutfullint = hW2anticut_scaled->Integral(ebound_l,ebound_h);

  //Set up fits and fit
  TF1 *tfitanticut = new TF1("tfitanticut",W2totalanticut_sgfit,ebound_l,ebound_h,8);
  tfitanticut->SetNpx(10000);  // Default is usually 100
  TF1 *bganticut = new TF1("bganticut",fits::g_p6fit,ebound_l,ebound_h,7);
  bganticut->SetNpx(10000);  // Default is usually 100
  TF1 *missed = new TF1("missed",fits::g_sgfit,ebound_l,ebound_h,4);
  missed->SetNpx(10000);  // Default is usually 100
  tfitanticut->SetLineColor(kGreen);
  hW2anticut_scaled->SetTitle("W^{2}, Acceptance cut and HCal best cluster elastic anticut");
  hW2anticut_scaled->Fit("tfitanticut","RBM");

  double *tparanticut = tfitanticut->GetParameters();

  //Get fit parameters for bg function and draw identical function on canvas
  bganticut->SetParameters(&tparanticut[1]);
  bganticut->SetLineColor(kRed);
  bganticut->SetFillColor(kRed);
  bganticut->SetFillStyle(3005);
  bganticut->Draw("same");
  missed->SetParameter(0,tparanticut[0]);
  missed->SetParameter(1,W2elas_sgfitpars[1]);
  missed->SetParameter(2,W2elas_sgfitpars[2]);
  missed->SetParameter(3,W2elas_sgfitpars[3]);
  missed->SetLineColor(kBlue);
  missed->SetLineWidth(2);
  missed->SetFillColorAlpha(kBlue,0.35);
  missed->SetFillStyle(3003);
  missed->Draw("same");
  
  //get missed elastics from anticut
  double anticut_totaldiff = 0.;
  double anticut_totalerr = 0.;
  util::subtractFunctionAndGetTotalAndError(hW2anticut_scaled,bganticut,ebound_l,ebound_h,anticut_totaldiff,ebound_l,ebound_h,anticut_totalerr,false);

  //Calculate efficiency and final error
  double effanticut = ( (totaldiff-anticut_totaldiff) / totaldiff )*100.;
  double effanticuterr = effanticut*sqrt(pow(totalerr/totaldiff,2)+pow(anticut_totalerr/anticut_totaldiff,2));

  //Add a legend to the canvas
  auto anticutlegend = new TLegend(0.1,0.6,0.5,0.9);
  anticutlegend->AddEntry( missed, "Tight elastic cut shape (scaled)", "l");
  anticutlegend->AddEntry( bganticut, "Background (polynomial)", "l");
  anticutlegend->AddEntry( tfitanticut, "Total sum fit", "l");
  anticutlegend->AddEntry( (TObject*)0, "", "");
  anticutlegend->AddEntry( (TObject*)0, Form("Proton spot elastic shape: %0.1f sigma",exp_spot_sigma), "" );
  anticutlegend->AddEntry( (TObject*)0, Form("Number of Missed Elastics: %d",(Int_t)anticut_totaldiff), "");
  anticutlegend->AddEntry( (TObject*)0, Form("Detection Efficiency: %0.3f%% +/- %0.3f%%",effanticut,effanticuterr), "");
  anticutlegend->Draw();

  //simply binomial error for comparison
  double W2_binomial = sqrt((effanticut/100.*(1-effanticut/100.))/totaldiff)*100.;
  cout << "Binomial error on W2 inclusive anticut: " << W2_binomial << endl;

  c1->Write();

  /////////////////////
  //DX SIDEBAND METHOD

  //Make canvas for hcal dx detected / expected 
  TCanvas *c2 = new TCanvas("c2",Form("HDE dx sideband sbs%d mag%d",kine,magset),1200,500);
  gStyle->SetPalette(55);
  gStyle->SetOptStat(0);
  c2->Divide(2,1);
  c2->cd(1);

  //Draw the expected elastic extraction first
  hW2_scaled->Draw();
  bg->Draw("same");
  elastics->Draw("same");

  auto sbwlegend = new TLegend(0.1,0.7,0.5,0.9);
  sbwlegend->AddEntry( elastics, "Tight elastic cut (sk gaus)", "l");
  sbwlegend->AddEntry( bg, "Background (poly)", "l");
  sbwlegend->AddEntry( tfit, "Total sum fit", "l");
  sbwlegend->AddEntry( (TObject*)0, "", "");
  sbwlegend->AddEntry( (TObject*)0, Form("Exp Elastics: %d",(Int_t)totaldiff), "");
  sbwlegend->Draw();

  c2->cd(2);

  Double_t lowdxcut = dxmean-dx_sb_width_sigma*dxsigma;
  Double_t highdxcut = dxmean+dx_sb_width_sigma*dxsigma;
  //Double_t lowdxcut = dxmean-0.01-dx_sb_width_sigma*dxsigma;
  //Double_t highdxcut = dxmean-0.01+dx_sb_width_sigma*dxsigma;
  Double_t lowdxrange = dxmean-dx_fit_width_sigma*dxsigma;
  Double_t highdxrange = dxmean+dx_fit_width_sigma*dxsigma;

  TH1D *hdx_gcut_scaled;
  if(scoopclus)
    hdx_gcut_scaled = util::cloneAndCutHistogram(hdx_gcut_scoop,lowdxrange,highdxrange);
  else
    hdx_gcut_scaled = util::cloneAndCutHistogram(hdx_gcut,lowdxrange,highdxrange);

  //set reject point range
  SBpol4rej_b = lowdxcut;
  SBpol4rej_e = highdxcut;

  TF1 *bgrpfit = new TF1("bgrpfit",BGfit,lowdxrange,highdxrange,5);
  bgrpfit->SetLineColor(kRed);
  hdx_gcut_scaled->SetTitle("dx, no cuts");
  hdx_gcut_scaled->GetXaxis()->SetTitle("m");
  hdx_gcut_scaled->Fit("bgrpfit","RBM");

  double *bgrppar = bgrpfit->GetParameters();

  TF1 *bgrp = new TF1("bgrp",fits::g_p4fit,lowdxrange,highdxrange,5);
  bgrp->SetParameters(&bgrppar[0]);
  bgrp->SetLineColor(kRed);
  bgrp->SetFillColor(kRed);
  bgrp->SetFillStyle(3005);
  bgrp->Draw("same");

  //get total between bg and dx total
  double sb_totaldiff, sb_totalerr; util::subtractFunctionAndGetTotalAndError(hdx_gcut_scaled,
									      bgrp,
									      lowdxcut,
									      highdxcut,
									      sb_totaldiff,
									      lowdxrange,
									      highdxrange,
									      sb_totalerr,
									      true);

  // Retrieve the minimum and maximum y-values of the histogram's y-axis for line placement
  double yMin = hdx_gcut_scaled->GetMinimum();
  double yMax = hdx_gcut_scaled->GetMaximum();

  // Create lines at xLow and xHigh
  TLine* lineLow = new TLine(lowdxcut, yMin, lowdxcut, yMax*1.04);
  TLine* lineHigh = new TLine(highdxcut, yMin, highdxcut, yMax*1.04);

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
  double hde_sb = sb_totaldiff/totaldiff*100.;
  double hde_sb_err = hde_sb*sqrt(pow(totalerr/totaldiff,2)+pow(sb_totalerr/sb_totaldiff,2));

  // Add the histogram, fit function, and lines to the legend
  sblegend->AddEntry( hdx_gcut_scaled, "Data", "lpe" );
  sblegend->AddEntry( bgrpfit, "Background Fit", "l" );
  sblegend->AddEntry( lineLow, "Sideband Limits", "l" );
  sblegend->AddEntry( (TObject*)0, "", "" );
  //sblegend->AddEntry( (TObject*)0, Form("Proton spot elastic shape: %0.1f sigma",exp_spot_sigma), "" );
  //sblegend->AddEntry( (TObject*)0, Form("Sideband limits: low %0.2f, high %0.2f",lowdxcut,highdxcut), "" );
  sblegend->AddEntry( (TObject*)0, Form("Total Elas: %d",(Int_t)sb_totaldiff), "" );
  sblegend->AddEntry( (TObject*)0, "Proton Det Eff", "" );
  sblegend->AddEntry( (TObject*)0, Form("%0.3f%% +/- %0.3f%%",hde_sb,hde_sb_err), "" );
  // Draw the legend
  sblegend->Draw();
  c2->Write();

  //simply binomial error for comparison
  double W2_sideband = sqrt((hde_sb/100.*(1-hde_sb/100.))/totaldiff)*100.;
  cout << "Binomial error on W2 inclusive anticut: " << W2_sideband << endl;

  cout << sqrt((0.96923*(1-0.96923))/totaldiff)*100. << endl;

  cout << endl << "Analysis complete. Outfile written to " << outfilename << endl;

  // Send time efficiency report to console
  st->Stop();
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
