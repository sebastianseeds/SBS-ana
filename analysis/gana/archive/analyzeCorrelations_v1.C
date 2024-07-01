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

//Total fits using Interpolate with elastic signal histo and 4th order poly fit to bg
TH1D *hdx_p;
TH1D *hdx_n;

Double_t fitSlice(double *x, double *par){
  double dx = x[0];
  double proton_scale = par[0];
  double neutron_scale = par[1];
  double proton = proton_scale * hdx_p->Interpolate(dx);
  double neutron = neutron_scale * hdx_n->Interpolate(dx);
  return proton + neutron + fits::g_p4fit(x,&par[2]);
}

//Side-band pol4 fit
//sbs9: -1.2 to 0.8
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

void analyzeCorrelations(int kine, int mag, int pass, double sb_min, double sb_max, int N) {

  //set up files and paths
  //std::string outdir_path = gSystem->Getenv("OUT_DIR");

  std::string outdir_path = "/lustre19/expphy/volatile/halla/sbs/seeds";
  std::string fin_path = outdir_path + Form("/gmn_analysis/dx_correlations_sbs%d_mag%d_pass%d.root",kine,mag,pass);
  std::string fin_path_mc = outdir_path + Form("/gmn_analysis/mc_rc_both_out_sbs%d_mag%d.root",kine,mag);
  std::string fout_path = outdir_path + Form("/gmn_analysis/analyze_correlations_sbs%d_pass%d.root",kine,pass);

  // Open the correlations root file
  TFile* inputFile = new TFile(fin_path.c_str());
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Error opening file: " << fin_path << std::endl;
    return;
  }

  // Open the mc root file
  TFile* inputFile_mc = new TFile(fin_path_mc.c_str());
  if (!inputFile_mc || inputFile_mc->IsZombie()) {
    std::cerr << "Error opening file: " << fin_path_mc << std::endl;
    return;
  }

  // Get hdx_p and hdx_n histograms from B.root
  hdx_p = dynamic_cast<TH1D*>(inputFile_mc->Get("hdx_cut_p"));
  hdx_n = dynamic_cast<TH1D*>(inputFile_mc->Get("hdx_cut_n"));
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

  // Loop over all TH2D histograms in A.root
  TIter next(inputFile->GetListOfKeys());
  TKey* key;
  while ((key = dynamic_cast<TKey*>(next()))) {

    if (strcmp(key->GetClassName(), "TH2D") != 0) continue; // Process only TH2D objects

    TH2D* hist = dynamic_cast<TH2D*>(key->ReadObj());
    if (!hist) continue;

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
    std::string ratioHistName = std::string("ratio: ") + hist->GetName() + ";" + hist->GetName() + ";n:p ratio";
    std::string yieldHistName = std::string("yield: ") + hist->GetName() + ";" + hist->GetName() + ";p+n yield";
    TH1D* ratioHist = new TH1D(ratioHistName.c_str(), ratioHistName.c_str(), N, xLow, xHigh);
    TH1D* yieldHist = new TH1D(yieldHistName.c_str(), yieldHistName.c_str(), N, xLow, xHigh);

    cout << endl << endl << "Working on correlation histogram " << ratioHistName << ". First bin with data: " << firstBinWithData << " (x = " << xLow << "), last bin with data: " << lastBinWithData << " (x = " << xHigh << ")" << endl << endl;

    // Create canvases for this histogram
    TCanvas* canvasSlices = new TCanvas(Form("canvasSlices_%s", hist->GetName()), hist->GetTitle(), 800, 600);
    TCanvas* canvasClones = new TCanvas(Form("canvasClones_%s", hist->GetName()), hist->GetTitle(), 800, 600);
    canvasSlices->Divide(3, 2); // Divide canvas into sub-pads for slices
    canvasClones->Divide(3, 2); // Divide canvas into sub-pads for slice clones

    // TH1D *cellslice[N];
    // TH1D *sliceclone[N];

    // Make N slices and fit each
    for (int i = 0; i < N; ++i) {
      
      cout << "Slice " << i << endl;

      int binLow = firstBinWithData + i * (lastBinWithData - firstBinWithData + 1) / N;
      int binHigh = firstBinWithData + (i + 1) * (lastBinWithData - firstBinWithData + 1) / N - 1;

      //Do slices for fits to MC
      TH1D *cellslice = hist->ProjectionY(Form("cellslice_%d", i+1), binLow, binHigh);
	
      //debugging
      std::cout << "Slice " << i << ": binLow = " << binLow << ", binHigh = " << binHigh
		<< ", Entries = " << cellslice->GetEntries() << std::endl;	       


      int sliceNBins = cellslice->GetNbinsX();
      double sliceMin = 0;
      double sliceMax = 0;

      // Find first bin with data
      for (int i = 1; i <= sliceNBins; ++i) {
        if (cellslice->GetBinContent(i) != 0) {
	  sliceMin = cellslice->GetXaxis()->GetBinCenter(i);
	  break;
        }
      }

      // Find last bin with data
      for (int i = sliceNBins; i >= 1; --i) {
        if (cellslice->GetBinContent(i) != 0) {
	  sliceMax = cellslice->GetXaxis()->GetBinCenter(i);
	  break;
        }
      }

      //Get Scale Factor ratio from this slice
      TF1 *fitFunc = new TF1("fitFunc", fitSlice, sliceMin, sliceMax, 7);
      fitFunc->SetNpx(5000);
      cellslice->Fit(fitFunc, "RBMQ0");  // Fitting the slice

      double scaleFactor1 = fitFunc->GetParameter(0);
      double scaleFactor2 = fitFunc->GetParameter(1);

      // Check for division by zero
      double scaleRatio = (scaleFactor1 != 0) ? scaleFactor2 / scaleFactor1 : 0;

      // Draw the slice on the slices canvas
      if( i<6 ){
	canvasSlices->cd(i + 1); // Go to the ith sub-pad
	cellslice->Draw();
	fitFunc->Draw("SAME");
	canvasSlices->Update();
      }	

      //So slices for sideband fits
      TH1D *slice_clone = hist->ProjectionY(Form("slice_clone_%d", i+1), binLow, binHigh);
      TH1D *slice_clone_2 = hist->ProjectionY(Form("slice_clone_2_%d", i+1), binLow, binHigh);

      //BG subtract with sideband fit and get yield for this slice
      TF1 *sbfit = new TF1("sbfit",BGfit,sliceMin,sliceMax,5);
      slice_clone->Fit("sbfit","RBMQ");

      Double_t *sbpar = sbfit->GetParameters();

      TF1 *sb = new TF1("sb",fits::g_p4fit,sliceMin,sliceMax,5);
      sb->SetLineColor(kGreen);
      sb->SetLineWidth(2);
      sb->SetParameters(&sbpar[0]);

      double totalEvents, sberr; util::subtractFunctionAndGetTotalAndError(slice_clone,
									   sb,
									   sb_min,
									   sb_max,
									   totalEvents,
									   sliceMin,
									   sliceMax,
									   sberr);

      // Draw the slice clone on the clones canvas
      if( i<6 ){
	canvasClones->cd(i + 1); // Go to the ith sub-pad of the clones canvas
	slice_clone_2->Draw();
	sb->Draw("SAME");
	canvasClones->Update();
      }

      //Fill the TH1Ds with ratios and yields
      ratioHist->SetBinContent(i, scaleRatio);
      yieldHist->SetBinContent(i, totalEvents);

      cout << "For bin " << i << ", ratio " << scaleRatio << " and yield " << totalEvents << endl;

    }

    // Before writing the histogram, check for invalid values
    for (int i = 1; i <= ratioHist->GetNbinsX(); ++i) {
      if (std::isinf(ratioHist->GetBinContent(i)) || std::isnan(ratioHist->GetBinContent(i))) {
	std::cerr << "Invalid value found in histogram " << ratioHist->GetName() << " bin " << i << std::endl;
	ratioHist->SetBinContent(i, 0); // Setting to 0 or some default value
      }
    }

    // Write the canvases to the output ROOT file
    canvasSlices->Write();
    canvasClones->Write();

    cout << "Got here 1" << endl;

    ratioHist->Draw();
    ratioHist->Write();

    yieldHist->Draw();
    yieldHist->Write();
  }

  cout << "Got here 1" << endl;

  cout << "All plots created. Output file located here: " << fout_path << endl;

}
