//sseeds 1.22.23

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TKey.h>
#include <TClass.h>
#include <iostream>
#include <vector>
#include <string>
#include "../../include/gmn.h"

double hcalfit_l = econst::hcalposXi_mc; //lower fit/bin limit for hcal dx plots (m)
double hcalfit_h = econst::hcalposXf_mc; //upper fit/bin limit for hcal dx plots (m)

const int nregions = 7; //total cut variables to look at
const int reportSamples = 6;

//Fit options
std::string fitopt = "RLQ0";

//set up exclusions
const vector<std::string> exclusions = {"hist_bb_tr_n",
					"hist_bb_sh_nclus",
					"hist_bb_sh_colblk",
					"hist_bb_sh_rowblk",
					"hist_sbs_hcal_nclus",
					"hist_hcalx",
					"hist_hcaly",
					"hist_thetapq_n",
					"hist_thetapq_p",
					"hist_hcalnblk",
					"hist_fiducial_sig_y",
					"hist_fiducial_sig_x"};

//Total fits using Interpolate with elastic signal histo and 4th order poly fit to bg
TH1D *hdx_p;
TH1D *hdx_n;

Double_t fitSlice_nobg(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0]);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0]);
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron;
}

Double_t fitSliceShift_nobg(double *x, double *par){
  // MC float params
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double dx_shift_p = par[2]; // Shift for proton histogram
  double dx_shift_n = par[3]; // Shift for neutron histogram  

  // Apply shifts before interpolation
  double proton = proton_scale * hdx_p->Interpolate(x[0] - dx_shift_p);
  double neutron = neutron_scale * hdx_n->Interpolate(x[0] - dx_shift_n);
  
  // Use the remaining parameters for fits::g_p4fit, starting from par[4]
  return proton + neutron;
}

std::vector<std::pair<double, double>> fitAndFineFit(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, const std::string& fitOptions = "RBMQ0");

//main. kine=kinematic, mag=fieldsetting, pass=pass#, N=cutvarsliceN
void mc_stability(int kine=4, int mag=50, int pass=2, int minEvents=500, int N=6, bool save=true) {

  //set draw params
  gStyle->SetPalette(55);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0);

  std::string outdir_path = "/lustre19/expphy/volatile/halla/sbs/seeds";
  std::string fin_path = outdir_path + Form("/gmn_analysis/dx_mc_sbs%d_mag%d_pass%d_rd2.root",kine,mag,pass);
  std::string fout_path = outdir_path + Form("/gmn_analysis/mc_stability_sbs%d_mag%d_pass%d_Nslices%d.root",kine,mag,pass,N);

  // Open the correlations root file
  TFile* inputFile = new TFile(fin_path.c_str());
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Error opening file: " << fin_path << std::endl;
    return;
  }

  // Get hdx_p and hdx_n histograms
  hdx_p = dynamic_cast<TH1D*>(inputFile->Get("hdx_p"));
  TH1D *hdxp = (TH1D*)(hdx_p->Clone("hdxp"));
  hdx_n = dynamic_cast<TH1D*>(inputFile->Get("hdx_n"));
  TH1D *hdxn = (TH1D*)(hdx_n->Clone("hdxn"));

  if (!hdx_p || !hdx_n) {
    std::cerr << "Error: hdx_p or hdx_n not found in " << fin_path << std::endl;
    inputFile->Close();
    return;
  }

  // Open an output ROOT file
  TFile* outputFile = new TFile(fout_path.c_str(), "RECREATE");

  struct ReportData {
    TH1D* sliceHistogram = nullptr;
    double pscale = 0;
    double nscale = 0;
    double perr = 0;
    double nerr = 0;
    double winLow = 0;
    double winHigh = 0;
    double nev = 0;
    double scaleratio = 0;
    std::string type;

    // Constructor for easier initialization (optional)
    ReportData(TH1D* hist = nullptr, double ps = 0, double ns = 0, double pe = 0, double ne = 0, 
               double wl = 0, double wh = 0, double nv = 0, double sr = 0, std::string t = "") 
      : sliceHistogram(hist), pscale(ps), nscale(ns), perr(pe), nerr(ne), 
	winLow(wl), winHigh(wh), nev(nv), scaleratio(sr), type(t) {}
  };

  ReportData reports[nregions][reportSamples];

  // Loop over all TH2D histograms
  int iteridx = 0;
  TIter next(inputFile->GetListOfKeys());
  TKey* key;
  while ((key = dynamic_cast<TKey*>(next()))) {

    if (strcmp(key->GetClassName(), "TH2D") != 0) continue; // Process only TH2D objects

    TH2D* hist = dynamic_cast<TH2D*>(key->ReadObj());
    if (!hist) continue;

    // Skip if the histogram name starts with "coorhist"
    std::string histName = hist->GetName();
    if (histName.find("coorhist") == 0) continue;

    if (histName.find("hist_inel") == 0) continue; // Skip if name contains inel

    bool skip = false;

    for( size_t name=0; name<exclusions.size(); ++name){

      if (histName.find(exclusions[name]) == 0){
	cout << "Found exclusion: " << exclusions[name] << ", skipping..." << endl;
	skip=true;
	break;
      }
    }

    //skip on exclusions
    if(skip)
      continue;

    // Skip all histograms that are exclusively proton (_p) or neutron (_n)
    // Calculate the positions for "_n" and "_p" if they are at the end of the histName
    size_t pos_n = histName.length() - 2; // Length of "_n" is 2
    size_t pos_p = histName.length() - 2; // Length of "_p" is also 2

    // Check if the histogram name ends with "_n" or "_p"
    bool endsWithN = histName.rfind("_n") == pos_n;
    bool endsWithP = histName.rfind("_p") == pos_p;

    if (endsWithN || endsWithP) continue;

    // Create a canvas for this TH2D histogram
    TCanvas* th2dCanvas = new TCanvas(Form("Full 2D %s", hist->GetName()), Form("Full 2D %s", hist->GetName()), 800, 600);
    th2dCanvas->cd();
    hist->Draw("COLZ");  // Draw the histogram
    th2dCanvas->Write();
    
    delete th2dCanvas;

    // Determine the first and last bin with data along the x-axis
    int firstBinWithData = 1;
    int lastBinWithData = hist->GetNbinsX();
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
      if (hist->ProjectionY("_py", i, i)->GetEntries() > minEvents) {
	firstBinWithData = i;
	break;
      }
    }
    for (int i = hist->GetNbinsX(); i >= 1; --i) {
      if (hist->ProjectionY("_py", i, i)->GetEntries() > minEvents) {
	lastBinWithData = i;
	break;
      }
    }

    if(firstBinWithData == lastBinWithData || firstBinWithData > lastBinWithData){
      cout << "WARNING: first and last bin with data don't leave data to analyze. Skipping..." << endl;
      continue;
    }
      
    if(lastBinWithData - firstBinWithData <= N){
      cout << "WARNING: N exceeds total bins with data. Proceeding anyway..." << endl;
    }

    // Get x-axis values corresponding to these bins
    double xLow = hist->GetXaxis()->GetBinLowEdge(firstBinWithData);
    double xHigh = hist->GetXaxis()->GetBinUpEdge(lastBinWithData);

    // Create TH1D to store ratios and yields with the same x-axis range as the TH2D histogram
    std::string ratioHistName = std::string("ratio: ") + hist->GetName() + ";" + hist->GetName() + ";scale factor n:p ratio";
    std::string yieldHistName = std::string("yield: ") + hist->GetName() + ";" + hist->GetName() + ";p+n yield";
    std::string ratioHistName_w = std::string("R, weighted proj: ") + hist->GetName() + "; scale ratio";;
 
    TH1D* ratioHist = new TH1D(ratioHistName.c_str(), ratioHistName.c_str(), N, xLow, xHigh);
    TH1D* yieldHist = new TH1D(yieldHistName.c_str(), yieldHistName.c_str(), N, xLow, xHigh);
    TH1D* ratio_w = new TH1D(ratioHistName_w.c_str(), ratioHistName_w.c_str(), 10*N, 0.5, 1.5);

    cout << endl << endl << "Working on correlation histogram " << ratioHistName << ". First bin with data: " << firstBinWithData << " (x = " << xLow << "), last bin with data: " << lastBinWithData << " (x = " << xHigh << ")" << endl << endl;

    TH1D *cellslice[N];

    double scalefactor[N][7];
    double scaleerr[N][7];

    double sMin[N];
    double sMax[N];

    double scaleratio[N];
    double scaleratio_err[N];
    double totalevents[N];

    // Arrays for scale factors for hdx_p and hdx_n for each slice after pol4 bg subtract
    double pScale[N];
    double nScale[N];

    // store slice window for reporting, x value for bin locations, TH2D
    double bin_low[N];
    double bin_high[N];

    int samplingInterval = max(1, N / reportSamples); // Ensure the interval is at least 1

    // Make N slices and fit each
    for (int i = 0; i < N; ++i) {
      
      // Calculate the sample index for the current slice
      int currentSampleIndex = i / samplingInterval;
    
      // Determine if the current slice is one of the slices to report
      bool isReportSlice = (i % samplingInterval == 0) && (currentSampleIndex < reportSamples);

      if(isReportSlice)
	cout << "Reporting out on hist " << histName << " on slice " << i << endl;

      //Get range for slice in TH2D (bins)
      int binLow = firstBinWithData + i * (lastBinWithData - firstBinWithData + 1) / N;
      int binHigh = firstBinWithData + (i + 1) * (lastBinWithData - firstBinWithData + 1) / N - 1;

      //(x values)
      bin_low[i] = hist->GetXaxis()->GetBinCenter(binLow);
      bin_high[i] = hist->GetXaxis()->GetBinCenter(binHigh);

      //Do slices for fits to MC
      cellslice[i] = hist->ProjectionY(Form("cellslice_%d", i+1), binLow, binHigh);
      TH1D *reportslice = hist->ProjectionY(Form("reportslice_%d_%s", i+1, histName.c_str()), binLow, binHigh);

      // if(isReportSlice)
      // 	reportslices[iteridx][currentSampleIndex] = (TH1D*)(cellslice[i]->Clone(Form("reportslice_%d_%d", iteridx, i)));

      //Get TH1D x axis window parameters
      int sliceNBins = cellslice[i]->GetNbinsX();
      double sliceMin = 0;
      double sliceMax = 0;

      int slice_nEntries = cellslice[i]->Integral(1,sliceNBins); // Total number of entries in the slice
      double slice_error = 0;
      if (slice_nEntries > 0) {
	slice_error = 1/sqrt(slice_nEntries);
      }
      
      // Find first bin with data on this slice
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

      if( sMin[i]>sMax[i] ){
	cout << "ERROR: address sMin/sMax logic." << endl;
	return;
      }

      //Get Scale Factor ratio from this slice
      std::pair<double,double> fullnobgQual;
      auto fullpar_nobg_vector = fitAndFineFit(cellslice[i], "sliceFit_nobg", "fitSlice_nobg", 2, hcalfit_l, hcalfit_h, fullnobgQual, fitopt.c_str());

      for(int k = 0; k < 2; ++k) {
	double param = fullpar_nobg_vector[k].first;
	double error = fullpar_nobg_vector[k].second;

	// Directly populate the scalefactor and scaleerr arrays
	scalefactor[i][k] = param;
	scaleerr[i][k] = error;
      }

      // Calculate the scale factor ratio and its error
      if (scalefactor[i][0] != 0) { // Avoid division by zero
	scaleratio[i] = scalefactor[i][1] / scalefactor[i][0];
	scaleratio_err[i] = sqrt(pow(scaleerr[i][1]/scalefactor[i][1], 2) + pow(scaleerr[i][0]/scalefactor[i][0], 2)) * scaleratio[i];
      } else {
	scaleratio[i] = 0;
	scaleratio_err[i] = 0;
      }

      if(isReportSlice){

	reports[iteridx][currentSampleIndex] = ReportData( reportslice, 
							   fullpar_nobg_vector[0].first,
							   fullpar_nobg_vector[1].first,
							   fullpar_nobg_vector[0].second,
							   fullpar_nobg_vector[1].second,
							   bin_low[i],
							   bin_high[i],
							   slice_nEntries,
							   scaleratio[i],
							   histName );

      }
  

      // Fill the TH1Ds with ratios and yields
      ratioHist->SetBinContent(i, scaleratio[i]);
      ratioHist->SetBinError(i, scaleratio_err[i]);

    } //endloop over N slices

    //Now get weighted histograms
    // Loop over all bins and print the value in each bin
    int ratioHist_nBins = ratioHist->GetXaxis()->GetNbins();
    for (int k = 1; k <= ratioHist_nBins; k++) { // bin indices start from 1
      double ratioHist_binValue = ratioHist->GetBinContent(k);
      double ratioHist_binErr = ratioHist->GetBinError(k);
    
      ratio_w->Fill(ratioHist_binValue,ratioHist_binErr);

    }

    // Before writing the histogram, check for invalid values
    for (int i = 1; i <= ratioHist->GetNbinsX(); ++i) {
      if (std::isinf(ratioHist->GetBinContent(i)) || std::isnan(ratioHist->GetBinContent(i))) {
	std::cerr << "Invalid value found in histogram " << ratioHist->GetName() << " bin " << i << std::endl;
	ratioHist->SetBinContent(i, 0); // Setting to 0 or some default value
      }
    }

    //Write the raw histograms to file
    ratioHist->Write();
    ratio_w->Write();

    //Write out canvas for ratio hisotgrams
    TCanvas *Rcanvas = new TCanvas(Form("R_{sf} %s",histName.c_str()), Form("R_{sf} %s",histName.c_str()), 800, 600);
    Rcanvas->cd();

    //ratio_w->Draw("hist");
    ratio_w->Draw("E");

    // Retrieve the RMS (Root Mean Square) of the histogram
    double rms = ratio_w->GetRMS();
    // If you also need the RMS error, you can retrieve it as follows
    double rmsError = ratio_w->GetRMSError();

    // Create a legend to display the RMS
    TLegend *legend = new TLegend(0.6, 0.7, 0.9, 0.9); // Adjust the coordinates to fit your canvas
    legend->SetHeader(Form("R_{sf} %s slices",histName.c_str()), "C"); // Optional header

    // Format the RMS text and add it to the legend
    char rmsText[255];
    sprintf(rmsText, "RMS = %.2f #pm %.2f", rms, rmsError);
    legend->AddEntry((TObject*)0, rmsText, "");

    // Draw the legend on the canvas
    legend->Draw();

    // Write the histogram to the current directory or file
    Rcanvas->Write();

    iteridx++;

  } //endloop over cut variables

  for (int r = 0; r < nregions; ++r) {
    TCanvas* canvasSlices = new TCanvas(Form("dx slices over %s", reports[r][0].type.c_str()), Form("dx slices over %s", reports[r][0].type.c_str()), 1600, 800);
    canvasSlices->Divide(3, 2); // Adjust the division based on reportSamples or desired layout

    for (int i = 0; i < reportSamples; ++i) {
      // Check if the histogram is valid; assuming 'i > N' check is no longer needed due to struct use
      if (reports[r][i].sliceHistogram == nullptr)
	continue;

      canvasSlices->cd(i + 1); // Go to the ith sub-pad

      TH1D* hdx_p_clone = (TH1D*)(hdxp->Clone(Form("hdx_p_clone_%d_%d", r, i)));
      TH1D* hdx_n_clone = (TH1D*)(hdxn->Clone(Form("hdx_n_clone_%d_%d", r, i)));

      reports[r][i].sliceHistogram->SetTitle(Form("slice window: %0.2f to %0.2f, Nev: %0.0f, ratio: %0.2f", reports[r][i].winLow, reports[r][i].winHigh, reports[r][i].nev, reports[r][i].scaleratio));
      reports[r][i].sliceHistogram->SetLineWidth(1);
      reports[r][i].sliceHistogram->Draw();

      hdx_p_clone->Scale(reports[r][i].pscale);
      hdx_p_clone->SetLineColor(kRed - 5);
      hdx_p_clone->SetLineWidth(1);
      hdx_p_clone->Draw("E same");

      hdx_n_clone->Scale(reports[r][i].nscale);
      hdx_n_clone->SetLineColor(kBlue - 5);
      hdx_n_clone->SetLineWidth(1);
      hdx_n_clone->Draw("E same");

      TF1* mcfit = new TF1(Form("mcfit_%d_%d", r, i), fitSlice_nobg, hcalfit_l, hcalfit_h, 2);
      mcfit->SetParameter(0, reports[r][i].pscale);
      mcfit->SetParameter(1, reports[r][i].nscale);
      mcfit->SetParError(0, reports[r][i].perr);
      mcfit->SetParError(1, reports[r][i].nerr); // Fixed to set error for parameter 1

      int transparentGreen = TColor::GetColorTransparent(kGreen, 0.3);

      if(!save){
	mcfit->SetLineColor(kGreen);
	mcfit->SetLineWidth(0);
	mcfit->SetFillColor(transparentGreen); // Use the transparent green color
	mcfit->SetFillStyle(1001); // Solid fill
	mcfit->Draw("SAME LF2"); // Draw the fit with line and filled area
      }else{
	mcfit->SetLineColor(kGreen);
	mcfit->Draw("same");
      }
	
      canvasSlices->Update();
    }

    // Optional: Save or further process the canvas
    if(save)
      canvasSlices->Write();
  }

  cout << "All plots created. Output file located here: " << fout_path << endl;

}

//Function to fit then fine fit and return all fit parameters
std::vector<std::pair<double, double>> fitAndFineFit(TH1D* histogram, const std::string& fitName, const std::string& fitFormula, int paramCount, double hcalfit_l, double hcalfit_h, std::pair<double,double>& fitqual, const std::string& fitOptions = "RBMQ0") {

  TF1* fit = new TF1(fitName.c_str(), fitFormula.c_str(), hcalfit_l, hcalfit_h, paramCount);
  //fit->SetNpx(5000);
  for (int i=0; i<paramCount; ++i) //reset parameters for this set
    fit->SetParameter(i,0);
  histogram->Fit(fit, fitOptions.c_str());

  std::vector<std::pair<double, double>> parametersAndErrors(paramCount);
  for (int i = 0; i < paramCount; ++i) {
    parametersAndErrors[i].first = fit->GetParameter(i); // Parameter value
    parametersAndErrors[i].second = fit->GetParError(i); // Parameter error
  }

  // Fine Fit
  TF1* fineFit = new TF1((fitName + "_fine").c_str(), fitFormula.c_str(), hcalfit_l, hcalfit_h, paramCount);
  std::vector<double> fineFitInitialParams(paramCount);
  for (int i = 0; i < paramCount; ++i) {
    fineFitInitialParams[i] = parametersAndErrors[i].first;
  }
  fineFit->SetParameters(fineFitInitialParams.data());
  histogram->Fit(fineFit, fitOptions.c_str());

  // Update parameters and errors with fine fit results
  for (int i = 0; i < paramCount; ++i) {
    parametersAndErrors[i].first = fineFit->GetParameter(i); // Fine fit parameter value
    parametersAndErrors[i].second = fineFit->GetParError(i); // Fine fit parameter error
  }

  fitqual.first = fineFit->GetChisquare();
  fitqual.second = fineFit->GetNDF();

  delete fit; // Delete fit to avoid carry-over
  delete fineFit; // Clean up
  return parametersAndErrors;
}
