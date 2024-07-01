//sseeds 1.22.23. Script to compare expected cut correlations from MC with those from data

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
#include "../../include/gmn.h"

double hcalfit_l = econst::hcalposXi_mc; //lower fit/bin limit for hcal dx plots (m)
double hcalfit_h = econst::hcalposXf_mc; //upper fit/bin limit for hcal dx plots (m)
double hcalr_l = -2.2; //lower fit/bin limit for hcal dx plots (m)
double hcalr_h = 0.8; //upper fit/bin limit for hcal dx plots (m)

//Total fits using Interpolate with elastic signal histo and 4th order poly fit to bg
TH1D *hdx_p;
TH1D *hdx_n;
string passkey = "flpx";

const vector<std::string> exclusions = {"hist_bb.tr.n",
					"hist_bb.gem.track.nhits",
					"hist_bb.sh.nclus",
					"hist_sbs.hcal.nclus"};

Double_t fitFull(double *x, double *par){
  double dx = x[0];
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double proton = proton_scale * hdx_p->Interpolate(dx);
  double neutron = neutron_scale * hdx_n->Interpolate(dx);
  return proton + neutron + fits::g_p4fit(x,&par[2]);
}

Double_t fitFull_rd2(double *x, double *par){
  double dx = x[0];
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double proton = proton_scale * hdx_p->Interpolate(dx);
  double neutron = neutron_scale * hdx_n->Interpolate(dx);
  return proton + neutron + fits::g_p4fit(x,&par[2]);
}

Double_t fitSlice(double *x, double *par){
  double dx = x[0];
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double proton = proton_scale * hdx_p->Interpolate(dx);
  double neutron = neutron_scale * hdx_n->Interpolate(dx);
  return proton + neutron + fits::g_p4fit(x,&par[2]);
}

Double_t fitSlice_rd2(double *x, double *par){
  double dx = x[0];
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double proton = proton_scale * hdx_p->Interpolate(dx);
  double neutron = neutron_scale * hdx_n->Interpolate(dx);
  return proton + neutron + fits::g_p4fit(x,&par[2]);
}

Double_t fitSlice_nobg(double *x, double *par){
  double dx = x[0];
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double proton = proton_scale * hdx_p->Interpolate(dx);
  double neutron = neutron_scale * hdx_n->Interpolate(dx);
  return proton + neutron;
}

//Side-band pol4 fit
//sbs9: -1.8 to 0.8
//sbs4 30p: -1.0 to 0.5 (careful)
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

//main. kine=kinematic, mag=fieldsetting, pass=pass#, sb_min/max=sidebandlimits, shiftX=shifttodxdata, N=cutvarsliceN
void analyzeCorrelations(int kine, int mag, int pass, double sb_min, double sb_max, double shiftX, int N) {

  //set up files and paths
  //std::string outdir_path = gSystem->Getenv("OUT_DIR");

  gStyle->SetPalette(55);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);

  // Define blinding factor
  Double_t blind_factor = util::generateRandomNumber(passkey);

  std::string outdir_path = "/lustre19/expphy/volatile/halla/sbs/seeds";
  std::string fin_path = outdir_path + Form("/gmn_analysis/dx_correlations_sbs%d_mag%d_pass%d_rd2.root",kine,mag,pass);
  std::string fin_path_mc = outdir_path + Form("/gmn_analysis/dx_mc_sbs%d_mag%d_pass%d_rd2.root",kine,mag,pass);
  std::string fout_path = outdir_path + Form("/gmn_analysis/analyze_correlations_sbs%d_mag%d_pass%d_rd2_finefit.root",kine,mag,pass);

  //std::string fout_path = outdir_path + "/garbage.root";

  // Open the correlations root file
  TFile* inputFile = new TFile(fin_path.c_str());
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Error opening file: " << fin_path << std::endl;
    return;
  }

  ///////////////////////////////////
  // Get the overall dx histogram
  TH1D *hdx_data = dynamic_cast<TH1D*>(inputFile->Get("hdx_allcut"));
  if (!hdx_data) {
    std::cerr << "Error: hdx_allcut not found in " << fin_path << std::endl;
    inputFile->Close();
    return;
  }

  // Open the mc root file
  TFile* inputFile_mc = new TFile(fin_path_mc.c_str());
  if (!inputFile_mc || inputFile_mc->IsZombie()) {
    std::cerr << "Error opening file: " << fin_path_mc << std::endl;
    return;
  }

  // Get hdx_p and hdx_n histograms
  hdx_p = dynamic_cast<TH1D*>(inputFile_mc->Get("hdx_p"));
  hdx_n = dynamic_cast<TH1D*>(inputFile_mc->Get("hdx_n"));
  if (!hdx_p || !hdx_n) {
    std::cerr << "Error: hdx_p or hdx_n not found in " << fin_path_mc << std::endl;
    inputFile->Close();
    inputFile_mc->Close();
    return;
  }

  // Open an output ROOT file
  TFile* outputFile = new TFile(fout_path.c_str(), "RECREATE");

  // Set up the sideband fit
  SBpol4rej_b = sb_min;
  SBpol4rej_e = sb_max;

  
  // Obtain overall fit to dx 
  gStyle->SetOptStat(0);
  TH1D *hdx_shifted = util::shiftHistogramX(hdx_data,shiftX);
  TH1D *hdx_shifted_clone = (TH1D*)(hdx_shifted->Clone("hdx_shifted_clone"));
  TH1D *hdx_shifted_clone2 = (TH1D*)(hdx_shifted->Clone("hdx_shifted_clone2"));
  TH1D *hdx_shifted_clone3 = (TH1D*)(hdx_shifted->Clone("hdx_shifted_clone3"));
  TH1D *hdx_p_clone = (TH1D*)(hdx_p->Clone("hdx_p_clone"));
  TH1D *hdx_n_clone = (TH1D*)(hdx_n->Clone("hdx_n_clone"));

  //Get Scale Factor ratio from this slice
  TF1* fullFit = new TF1("fullFit_rd1", fitFull, hcalfit_l, hcalfit_h, 7);
  fullFit->SetNpx(5000);
  hdx_shifted->Fit(fullFit, "RBMQ0");  // Fitting the hist

  double fullpar_array[7]={0.};
  for( int i=0; i<7; ++i ){
    fullpar_array[i] = fullFit->GetParameter(i);
  }

  //Fit again for fine tuning
  TF1* fullFit_rd2 = new TF1("fullFit_rd2", fitFull_rd2, hcalfit_l, hcalfit_h, 7);
  //fullFit_rd2->SetNpx(5000);
  fullFit_rd2->SetParameters(fullpar_array);

  hdx_shifted_clone->Fit(fullFit_rd2, "RBMQ0");  // Fitting the slice second with fine tuned params

  double fullpar_rd2_array[7]={0.};
  for( int i=0; i<7; ++i ){
    fullpar_rd2_array[i] = fullFit_rd2->GetParameter(i);
  }
  
  
  TCanvas *c2 = new TCanvas("c2","Interpolate/4th Order Poly BG",800,600);
  c2->cd();

  hdx_shifted_clone2->GetXaxis()->SetRangeUser(hcalr_l,hcalr_h);
  hdx_shifted_clone2->SetLineColor(kBlack);
  hdx_shifted_clone2->SetTitle("dx, Q^{2}=3.0 GeV^{2}");

  hdx_shifted_clone2->Draw("hist");

  fullFit_rd2->SetLineColor(kGreen);
  fullFit_rd2->SetLineWidth(1);
  fullFit_rd2->Draw("same");

  hdx_p_clone->Scale(fullpar_rd2_array[0]);
  hdx_p_clone->SetLineColor(kRed-5);
  hdx_p_clone->SetLineWidth(2);
  //hdx_p_clone->SetFillColorAlpha(kRed,0.35);
  //hdx_p_clone->SetFillStyle(3005);
  hdx_p_clone->Draw("E same");

  hdx_n_clone->Scale(fullpar_rd2_array[1]);
  hdx_n_clone->SetLineColor(kBlue-5);
  hdx_n_clone->SetLineWidth(2);
  //hdx_n_clone->SetFillColorAlpha(kBlue,0.35);
  //hdx_n_clone->SetFillStyle(3003);
  hdx_n_clone->Draw("E same");

  TF1 *bg = new TF1("bg",fits::g_p4fit,hcalfit_l,hcalfit_h,5);
  bg->SetParameters(&fullpar_rd2_array[2]);
  bg->SetLineColor(kBlack);
  bg->SetFillColorAlpha(kGray,0.35);
  //bg->SetFillStyle(3004);
  bg->SetFillStyle(1001);
  bg->Draw("same");

  double pscale = fullpar_rd2_array[0];
  double nscale = fullpar_rd2_array[1];

  //double np_par_ratio = blind_factor * nscale/pscale;
  double np_par_ratio = nscale/pscale*blind_factor;

  //Add a legend to the canvas
  auto leg = new TLegend(0.6,0.6,0.89,0.89);
  leg->AddEntry( hdx_shifted_clone2, "tight elastic cut, LD2 data", "l");
  leg->AddEntry( hdx_p_clone, "SIMC MC, proton", "l");
  leg->AddEntry( hdx_n_clone, "SIMC MC, neutron", "l");
  leg->AddEntry( bg, "pol4 fit to BG", "l");
  leg->AddEntry( fullFit_rd2, "Total fit", "l");
  leg->AddEntry( (TObject*)0, "", "");
  leg->AddEntry( (TObject*)0, Form("data N events : %0.0f",hdx_shifted_clone2->GetEntries()), "");
  leg->AddEntry( (TObject*)0, Form("MC N events : %0.0f",(hdx_p_clone->GetEntries()+hdx_n_clone->GetEntries())), "");
  leg->AddEntry( (TObject*)0, Form("n/p scale ratio R_{sf} : %0.3f",np_par_ratio), "");
  leg->AddEntry( (TObject*)0, Form("#chi^{2}: %0.3f",fullFit_rd2->GetChisquare()), "");
  leg->Draw("same");

  c2->Update();
  c2->Write();

  c2->SaveAs("~/gmn_plots/analyzeCorrelations_dx_out.pdf");
  
  TCanvas *c3 = new TCanvas("c3","Interpolate/4th Order Poly BG Residuals",800,600);
  c3->cd();
  //add the p and n histograms together for comparison
  TH1D* sumHistogram = new TH1D("sumHistogram", "Sum of Histograms", hdx_p_clone->GetNbinsX(), hdx_p_clone->GetXaxis()->GetXmin(), hdx_p_clone->GetXaxis()->GetXmax());
  sumHistogram->Add(hdx_p_clone, hdx_n_clone);

  TF1 *pol4 = new TF1("pol4", fits::g_p4fit, hcalfit_l, hcalfit_h, 5);
  pol4->SetParameters(&fullpar_rd2_array[2]); // Assuming indices 2-6 are pol4 parameters
  
  for (int bin = 1; bin <= hdx_shifted_clone3->GetNbinsX(); ++bin) {
    double bgValue = pol4->Eval(hdx_shifted_clone3->GetXaxis()->GetBinCenter(bin));
    hdx_shifted_clone3->SetBinContent(bin, hdx_shifted_clone3->GetBinContent(bin) - bgValue);
    hdx_shifted_clone3->SetBinError(bin, hdx_shifted_clone3->GetBinError(bin));
  }

  TH1D *hRes = util::makeResidualHisto("dx",hdx_shifted_clone3,sumHistogram,true,false);
  
  hRes->SetTitle("");
  hRes->GetXaxis()->SetRangeUser(hcalr_l,hcalr_h);

  double dxMaxValue = hdx_shifted_clone3->GetMaximum();
  hRes->GetYaxis()->SetRangeUser(-dxMaxValue/2,dxMaxValue/2);

  hRes->SetLineColor(kRed);
  hRes->SetLineWidth(2);
  hRes->Draw();

  c3->Update();
  c3->Write();

  //end dx full plot
  ////////////////////

  // Create a canvas with two pads
  TCanvas *cTotal = new TCanvas("cTotal", "Data with Fits and Residuals", 800, 800);
  TPad *pad1 = new TPad("pad1", "Pad with the fit", 0.0, 0.3, 1.0, 1.0);
  TPad *pad2 = new TPad("pad2", "Pad with the residuals", 0.0, 0.0, 1.0, 0.3);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.2);
  pad1->Draw();
  pad2->Draw();

  // First pad: Draw the histogram with the fits
  pad1->cd();
  hdx_shifted_clone2->Draw("hist");
  fullFit_rd2->Draw("same");
  hdx_p_clone->Draw("E same");
  hdx_n_clone->Draw("E same");
  bg->Draw("same");
  leg->Draw("same"); // Assuming 'leg' is your TLegend object

  // Second pad: Draw the residuals
  pad2->cd();
  hRes->Draw();

  // Adjust the y-axis of the residuals to be more visually appealing
  hRes->GetYaxis()->SetTitle("Residuals");
  hRes->GetYaxis()->SetTitleSize(0.07);
  hRes->GetYaxis()->SetTitleOffset(0.5);
  hRes->GetYaxis()->SetLabelSize(0.07);
  hRes->GetXaxis()->SetTitleSize(0.07);
  hRes->GetXaxis()->SetLabelSize(0.07);
  hRes->GetXaxis()->SetTitle("x_{hcal}-x_{exp}"); // Set the appropriate title

  // Draw a line at y=0 for residuals
  TF1 *zeroLine = new TF1("zeroLine", "0", hcalr_l, hcalr_h);
  zeroLine->SetLineColor(kBlack);
  zeroLine->SetLineStyle(2);
  zeroLine->Draw("same");

  // Update the canvases
  cTotal->cd();
  cTotal->Update();

  // Save the canvas to a file
  cTotal->SaveAs("~/gmn_plots/combined_dx_fits_residuals.pdf");

  return;

  // Loop over all TH2D histograms
  TIter next(inputFile->GetListOfKeys());
  TKey* key;
  while ((key = dynamic_cast<TKey*>(next()))) {

    if (strcmp(key->GetClassName(), "TH2D") != 0) continue; // Process only TH2D objects

    TH2D* hist = dynamic_cast<TH2D*>(key->ReadObj());
    if (!hist) continue;

    // Check if the histogram name starts with "coorhist" and skip if so
    std::string histName = hist->GetName();
    if (histName.find("coorhist") == 0) continue; // Skip if name doesn't start with "hist"

    bool skip = false;

    for( size_t name=0; name<exclusions.size(); ++name){

      if (histName.find(exclusions[name]) == 0){
	cout << "Found exclusion: " << exclusions[name] << endl;
	skip=true;
	break;
      }
    }

    if(skip)
      continue;

    // Create a canvas for this TH2D histogram
    TCanvas* th2dCanvas = new TCanvas(Form("canvas_%s", hist->GetName()), hist->GetTitle(), 800, 600);
    th2dCanvas->cd();
    hist->Draw("COLZ");  // Draw the histogram
    th2dCanvas->Write();

    // Determine the first and last bin with data along the x-axis
    int firstBinWithData = 1;
    int lastBinWithData = hist->GetNbinsX();
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
      if (hist->ProjectionY("_py", i, i)->GetEntries() > 0) {
	firstBinWithData = i;
	break;
      }
    }
    for (int i = hist->GetNbinsX(); i >= 1; --i) {
      if (hist->ProjectionY("_py", i, i)->GetEntries() > 0) {
	lastBinWithData = i;
	break;
      }
    }

    // Get x-axis values corresponding to these bins
    double xLow = hist->GetXaxis()->GetBinLowEdge(firstBinWithData);
    double xHigh = hist->GetXaxis()->GetBinUpEdge(lastBinWithData);

    // Create TH1D to store ratios and yields with the same x-axis range as the TH2D histogram
    std::string ratioHistName = std::string("ratio: ") + hist->GetName() + ";" + hist->GetName() + ";scale factor n:p ratio";
    std::string ratioHistNameBG = std::string("ratio_bgsub: ") + hist->GetName() + ";" + hist->GetName() + ";scale factor n:p ratio";
    std::string yieldHistName = std::string("yield: ") + hist->GetName() + ";" + hist->GetName() + ";p+n yield";
    std::string ratioHistName_w = std::string("R, weighted proj: ") + hist->GetName() + "; scale ratio";;
    std::string ratioHistNameBG_w = std::string("R_bg, weighted proj: ") + hist->GetName() + "; scale ratio background subtracted";
 
    TH1D* ratioHist = new TH1D(ratioHistName.c_str(), ratioHistName.c_str(), N, xLow, xHigh);
    TH1D* ratioHist_bg = new TH1D(ratioHistNameBG.c_str(), ratioHistNameBG.c_str(), N, xLow, xHigh);
    TH1D* yieldHist = new TH1D(yieldHistName.c_str(), yieldHistName.c_str(), N, xLow, xHigh);
    TH1D* ratio_w = new TH1D(ratioHistName_w.c_str(), ratioHistName_w.c_str(), N, -1, 3);
    TH1D* ratio_w_bg = new TH1D(ratioHistNameBG_w.c_str(), ratioHistNameBG_w.c_str(), N, -1, 3);

    cout << endl << endl << "Working on correlation histogram " << ratioHistName << ". First bin with data: " << firstBinWithData << " (x = " << xLow << "), last bin with data: " << lastBinWithData << " (x = " << xHigh << ")" << endl << endl;

    TH1D *cellslice[N];
    TH1D *cellslice_rd2[N];
    TH1D *cellslice_bgsub[N];
    TH1D *sliceclone[N];

    double scalefactor[N][7];
    double scalefactor_bg[N][2];
    double sbpar[N][5];

    double sMin[N];
    double sMax[N];

    double scaleratio[N];
    double scaleratio_err[N];
    double scaleratio_bgsub[N];
    double scaleratio_bgsub_err[N];
    double totalevents[N];

    // Arrays for scale factors for hdx_p and hdx_n for each slice after pol4 bg subtract
    double pScale[N];
    double nScale[N];
    double scaleratio_bg[N];
    double scaleratio_bg_err[N];

    // store slice window for reporting
    double bin_low[N];
    double bin_high[N];

    // Make N slices and fit each
    for (int i = 0; i < N; ++i) {
      
      //cout << "Slice " << i << endl;

      int binLow = firstBinWithData + i * (lastBinWithData - firstBinWithData + 1) / N;
      int binHigh = firstBinWithData + (i + 1) * (lastBinWithData - firstBinWithData + 1) / N - 1;

      bin_low[i] = hist->GetXaxis()->GetBinCenter(binLow);
      bin_high[i] = hist->GetXaxis()->GetBinCenter(binHigh);

      //Do slices for fits to MC
      cellslice[i] = hist->ProjectionY(Form("cellslice_%d", i+1), binLow, binHigh);
      cellslice_rd2[i] = hist->ProjectionY(Form("cellslice_rd2_%d", i+1), binLow, binHigh);

      if( cellslice[i] == nullptr )
	continue;

      int sliceNBins = cellslice[i]->GetNbinsX();
      double sliceMin = 0;
      double sliceMax = 0;

      int slice_nEntries = cellslice[i]->Integral(); // Total number of entries in the slice
      double slice_error = 1/sqrt(slice_nEntries); // Poisson error on histogram

      // Find first bin with data
      for (int bin = 1; bin <= sliceNBins; ++bin) {
        if (cellslice[i]->GetBinContent(bin) != 0) {
	  sliceMin = cellslice[i]->GetXaxis()->GetBinCenter(bin);
	  break;
        }
      }
      sMin[i] = sliceMin;

      // Find last bin with data
      for (int bin = sliceNBins; bin >= 1; --bin) {
        if (cellslice[i]->GetBinContent(bin) != 0) {
	  sliceMax = cellslice[i]->GetXaxis()->GetBinCenter(bin);
	  break;
        }
      }
      sMax[i] = sliceMax;

      //Get Scale Factor ratio from this slice
      TF1* fitFunc = new TF1("fitFunc", fitSlice, sliceMin, sliceMax, 7);
      fitFunc->SetNpx(5000);
      cellslice[i]->Fit(fitFunc, "RBMQ0");  // Fitting the slice

      //double* mcParams = fitFunc->GetParameters(); //0, proton; 1, neutron, 2-6, pol4

      double mcParams_array[7]={0.};
      for( int i=0; i<7; ++i ){
	mcParams_array[i] = fitFunc->GetParameter(i);
      }

      //Fit again for fine tuning
      TF1* fitFunc_rd2 = new TF1("fitFunc_rd2", fitSlice_rd2, sliceMin, sliceMax, 7);
      fitFunc_rd2->SetNpx(5000);
      fitFunc_rd2->SetParameters(mcParams_array);

      cellslice_rd2[i]->Fit(fitFunc_rd2, "RBMQ0");  // Fitting the slice second with fine tuned params

      //double* mcParams_rd2 = fitFunc_rd2->GetParameters();

      double mcParams_rd2_array[7]={0.};
      for( int i=0; i<7; ++i ){
	mcParams_rd2_array[i] = fitFunc_rd2->GetParameter(i);
      }

      // Get the error on the fit parameters
      int nParams = fitFunc->GetNpar(); // Gets the number of parameters
      double paramErrors[7]; // Array to store the errors

      for (int j = 0; j < nParams; ++j) {
	paramErrors[j] = fitFunc->GetParError(j); // Get the error for each parameter
      }

      // Populate the scalefactor array
      for (int p = 0; p < 7; ++p) scalefactor[i][p] = mcParams_rd2_array[p];

      // Calculate the scale factor ratio and its error
      if (scalefactor[i][0] != 0) { // Check for division by zero for the denominator
	scaleratio[i] = scalefactor[i][1] / scalefactor[i][0];
	scaleratio_err[i] = sqrt(pow(paramErrors[1]/scalefactor[i][1], 2) + pow(paramErrors[0]/scalefactor[i][0], 2)) * scaleratio[i]; // Error propagation formula
      } else {
	scaleratio[i] = 0;
	scaleratio_err[i] = 0;
      }

      //So slices for sideband fits
      sliceclone[i] = hist->ProjectionY(Form("slice_clone_%d", i+1), binLow, binHigh);

      //BG subtract with sideband fit and get yield for this slice
      TF1 *sbfit = new TF1("sbfit",BGfit,sliceMin,sliceMax,5);
      sliceclone[i]->Fit("sbfit","RBMQ");

      const double* sbParams = sbfit->GetParameters();
      for (int p = 0; p < 5; ++p) sbpar[i][p] = sbParams[p];

      double *sbpars = sbfit->GetParameters();

      TF1 *sb = new TF1("sb",fits::g_p4fit,sliceMin,sliceMax,5);
      sb->SetParameters(&sbpars[0]);

      double sberr; util::subtractFunctionAndGetTotalAndError(sliceclone[i],
							      sb,
							      sb_min,
							      sb_max,
							      totalevents[i],
							      sliceMin,
							      sliceMax,
							      sberr);

      //Fill the TH1Ds with ratios and yields
      ratioHist->SetBinContent(i, scaleratio[i]);
      ratioHist->SetBinError(i, scaleratio_err[i]);
      yieldHist->SetBinContent(i, totalevents[i]);
      yieldHist->SetBinError(i, sberr);

      //cout << "For bin " << i << ", ratio " << scaleratio[i] << " and yield " << totalevents[i] << endl;

      //Add scale factor after BG subtract for comparison.

      // Create polynomial background function from scalefactor parameters
      TF1 *pol4 = new TF1(Form("pol4_%d", i), fits::g_p4fit, sMin[i], sMax[i], 5);
      pol4->SetParameters(&scalefactor[i][2]); // Assuming indices 2-6 are pol4 parameters

      // Background subtract cellslice
      cellslice_bgsub[i] = (TH1D*)cellslice[i]->Clone(Form("bgSubtractedSlice_%d", i));
      for (int bin = 1; bin <= cellslice_bgsub[i]->GetNbinsX(); ++bin) {
	double bgValue = pol4->Eval(cellslice_bgsub[i]->GetXaxis()->GetBinCenter(bin));
	cellslice_bgsub[i]->SetBinContent(bin, cellslice_bgsub[i]->GetBinContent(bin) - bgValue);
	cellslice_bgsub[i]->SetBinError(bin, cellslice_bgsub[i]->GetBinError(bin));
      }

      // Define fit function for MC histograms
      TF1 *mcFit = new TF1(Form("mcFit_%d", i), fitSlice_nobg, sliceMin, sliceMax,2);
      mcFit->SetParameters(1, 1); // Initial guess for scaling factors

      // Fit the background subtracted cellslice with the MC histograms
      cellslice_bgsub[i]->Fit(mcFit, "RBMQ0");

      // Get the error on the fit parameters
      int nParams_bg = mcFit->GetNpar(); // Gets the number of parameters
      double paramErrors_bg[7]; // Array to store the errors

      for (int j = 0; j < nParams_bg; ++j) {
	paramErrors_bg[j] = mcFit->GetParError(j); // Get the error for each parameter
      }

      // Store the scale factors
      pScale[i] = mcFit->GetParameter(0);
      nScale[i] = mcFit->GetParameter(1);
      //scaleratio_bg[i] = nScale[i]/pScale[i];
      
      // Calculate the background subtracted scale factor ratio and its error
      if (pScale[i] != 0) { // Check for division by zero for the denominator
	scaleratio_bg[i] = nScale[i]/pScale[i];
	scaleratio_bg_err[i] = sqrt(pow(paramErrors_bg[1]/pScale[i], 2) + pow(paramErrors_bg[0]/nScale[i], 2)) * scaleratio_bg[i]; // Error propagation formula
      } else {
	scaleratio_bg[i] = 0;
	scaleratio_bg_err[i] = 0;
      }

      scalefactor_bg[i][0] = pScale[i];
      scalefactor_bg[i][1] = nScale[i];

      ratioHist_bg->SetBinContent(i, scaleratio_bg[i]);
      ratioHist_bg->SetBinError(i, scaleratio_bg_err[i]);

    }

    //Now get weighted histograms
    // Loop over all bins and print the value in each bin
    int ratioHist_nBins = ratioHist->GetXaxis()->GetNbins();
    for (int k = 1; k <= ratioHist_nBins; k++) { // bin indices start from 1
      double ratioHist_binValue = ratioHist->GetBinContent(k);
      double ratioHist_bg_binValue = ratioHist_bg->GetBinContent(k);
      double yieldHist_binValue = yieldHist->GetBinContent(k);
      
      ratio_w->Fill(ratioHist_binValue,yieldHist_binValue);
      ratio_w_bg->Fill(ratioHist_bg_binValue,yieldHist_binValue);
    }

    // Before writing the histogram, check for invalid values
    for (int i = 1; i <= ratioHist->GetNbinsX(); ++i) {
      if (std::isinf(ratioHist->GetBinContent(i)) || std::isnan(ratioHist->GetBinContent(i))) {
	std::cerr << "Invalid value found in histogram " << ratioHist->GetName() << " bin " << i << std::endl;
	ratioHist->SetBinContent(i, 0); // Setting to 0 or some default value
      }
    }

    

    //Now write out the mc canvas
    TCanvas* canvasSlices = new TCanvas(Form("MC Fits %s", hist->GetName()), hist->GetTitle(), 800, 600);
    canvasSlices->Divide(3, 2); // Divide canvas into sub-pads for slices
    for( int i=0; i<6; ++i ){
      canvasSlices->cd(i + 1); // Go to the ith sub-pad

      //Get even sample
      int j = i*(int)(N/6);
      //int j = i+20;

      if( cellslice[j] == nullptr )
	continue;

      cellslice[j]->SetTitle(Form("slice window: %0.2f - %0.2f, ratio: %0.2f",bin_low[j],bin_high[j],scaleratio[j]));
      cellslice[j]->Draw();
      TF1 *mcfit = new TF1("mcfit", fitSlice, sMin[i], sMax[i], 7);
      mcfit->SetParameters(&scalefactor[j][0]);
      mcfit->SetLineColor(kRed);
      mcfit->Draw("SAME");
      canvasSlices->Update();
    }	
    canvasSlices->Write();

    // //Now write overall fit to dx
    // TCanvas* dxfull = new TCanvas("dxfull","Interpolate/4th Order Poly BG", 800, 600);
    // dxfull->cd();
    // hdx_shifted_clone2->SetTitle("dx;x_{hcal}-x_{exp}");
    // hdx_shifted_clone2->Draw();
    // TF1 *fullfit = new TF1("fullfit", fitFull, hcalfit_l, hcalfit_h, 7);
    // fullfit->SetParameters(&fullpar_rd2_array[0]);
    // fullfit->SetLineColor(kRed);
    // fullfit->Draw("SAME");
    // dxfull->Update();

    // dxfull->Write();

    //Now write out the mc bg sub canvas
    TCanvas* canvasSlices_bg = new TCanvas(Form("MC Fits BG Sub %s", hist->GetName()), hist->GetTitle(), 800, 600);
    canvasSlices_bg->Divide(3, 2); // Divide canvas into sub-pads for slices
    for( int i=0; i<6; ++i ){
      canvasSlices_bg->cd(i + 1); // Go to the ith sub-pad

      //Get even sample
      int j = i*(int)(N/6);
      //int j = i+20;
      cellslice_bgsub[j]->SetTitle(Form("slice window: %0.2f - %0.2f, ratio bg: %0.2f",bin_low[j],bin_high[j],scaleratio_bg[j]));
      cellslice_bgsub[j]->Draw();
      TF1 *mcfit = new TF1("mcfit", fitSlice_nobg, sMin[i], sMax[i], 7);
      mcfit->SetParameters(&scalefactor_bg[j][0]);
      mcfit->SetLineColor(kRed);
      mcfit->Draw("SAME");
      canvasSlices_bg->Update();
    }	
    canvasSlices_bg->Write();

    //Now write out the sideband canvas
    TCanvas* canvasClones = new TCanvas(Form("Sideband Fits %s", hist->GetName()), hist->GetTitle(), 800, 600);
    canvasClones->Divide(3, 2); // Divide canvas into sub-pads for slices
    for( int i=0; i<6; ++i ){
      canvasClones->cd(i + 1); // Go to the ith sub-pad

      //Get even sample
      int j = i*(int)(N/6);
      //int j = i+20;
      sliceclone[j]->SetTitle(Form("slice window: %0.2f - %0.2f, p+n yield: %0.2f",bin_low[j],bin_high[j],totalevents[j]));
      sliceclone[j]->Draw();
      TF1 *sbfit = new TF1("sbfit",fits::g_p4fit,sMin[i],sMax[i],5);
      sbfit->SetParameters(&sbpar[j][0]);
      sbfit->SetLineColor(kGreen);
      sbfit->Draw("SAME");
      canvasClones->Update();
    }	
    canvasClones->Write();

    ratioHist->Draw();
    ratioHist->Write();

    ratioHist_bg->Draw();
    ratioHist_bg->Write();

    yieldHist->Draw();
    yieldHist->Write();

    ratio_w->Draw();
    ratio_w->Write();

    ratio_w_bg->Draw();
    ratio_w_bg->Write();
    
  }

  cout << "All plots created. Output file located here: " << fout_path << endl;

}
