//seeds
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TLegend.h>
#include <TH2D.h>
#include <iostream>
#include <vector>
#include <string>

//wide elastic
std::string globalcuts = "W2<1.2&&bb_ps_e>0.2&&bb_gem_track_nhits>3&&abs(bb_tr_vz)<0.072";

//remove randoms as best as possible and pions
std::string lglobalcuts = "bb_ps_e>0.2&&bb_gem_track_nhits>3&&abs(bb_tr_vz)<0.072";

const int minev = 150;
const int tdcbins = 45;
const double tdcllim = -15.0;
const double tdculim = 5.0;

// Helper function to fit histogram with a Gaussian and extract mean and sigma
void fitHistogram(TH1D* hist, double &mean, double &sigma) {
  int binMax = hist->GetMaximumBin();
  double maxVal = hist->GetBinContent(binMax);
  double fifthmax = maxVal / 5.0;
  double halfmax = maxVal / 2.0;
  
  int binLeft = binMax;
  while (binLeft > 1 && hist->GetBinContent(binLeft) > halfmax) {
    binLeft--;
  }
  
  int binRight = binMax;
  while (binRight < hist->GetNbinsX() && hist->GetBinContent(binRight) > fifthmax) {
    binRight++;
  }
  
  TF1* gaussFit = new TF1("gaussFit", "gaus", hist->GetBinCenter(binLeft), hist->GetBinCenter(binRight));
  hist->Fit(gaussFit, "R");
  
  mean = gaussFit->GetParameter(1);
  sigma = gaussFit->GetParameter(2);
}

void tdcInternalRes(int kine=4, bool loose=false) {
  gStyle->SetPalette(55);
  gStyle->SetOptStat(0110);
  gStyle->SetStatTextColor(kBlack);
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.08);
  gStyle->SetStatX(0.4);
  gStyle->SetStatY(0.8);

  std::string cuts;
  if(loose)
    cuts = lglobalcuts;
  else
    cuts = globalcuts;

  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fout_path = outdir_path + "/gmn_analysis/tdcres/timeinternalres_out.root";
  std::cout << "Setting up output path: " << fout_path << std::endl;

  TFile* outputFile = new TFile(fout_path.c_str(), "RECREATE");

  std::string fin_path = outdir_path + Form("/parse/parse_sbs%d_pass2_barebones_effz_lh2.root", kine);
  std::cout << "Setting up input path: " << fin_path << std::endl;
  TFile* inputFile = new TFile(fin_path.c_str());
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Error opening file: " << fin_path << std::endl;
    return;
  }

  TTree* tree = dynamic_cast<TTree*>(inputFile->Get("P"));
  if (!tree) {
    std::cerr << "Tree not found in file: " << fin_path << std::endl;
    inputFile->Close();
    return;
  }

  // Create base histograms for internal timing resolution
  std::vector<TH1D*> base_hists(288, nullptr);
  for (int i = 0; i < 288; ++i) {
    base_hists[i] = new TH1D(Form("base_hist_%d", i), Form("Base Histogram %d;ns;Counts", i), tdcbins, tdcllim, tdculim);
  }

  // Create a 2D histogram for the heat map
  TH2D* heat_map = new TH2D("heat_map", "Event Counts per Channel;Column;Row;Counts", 12, 0.5, 12.5, 24, 0.5, 24.5);

  // Account for timewalk with ADC energy per block
  for (int i = 0; i < 288; ++i) {
    std::string cut = cuts + Form("&&hcalcbid_b1==%d&&hcalcbid_b2==%d", i + 1, i + 2);
    tree->Draw(Form("(hcaltdc_b1-dttw_b1-hcaltdc_b2)>>base_hist_%d", i), cut.c_str(), "COLZ");

    // Fill the heat map
    int col = (i % 12) + 1;
    int row = (i / 12) + 1;
    heat_map->Fill(col, row, base_hists[i]->GetEntries());
  }

  // Create overall distribution of std devs
  TH1D* hcal_internal_timing_res = new TH1D("hcal_internal_timing_res", "HCal Internal Timing Resolution;ns;Counts", 20, 0, 2);

  for (auto hist : base_hists) {
    if (hist->GetEntries() > minev) {
      double mean, sigma;
      fitHistogram(hist, mean, sigma);
      hcal_internal_timing_res->Fill(sigma);
    }
  }

  // Create overall distribution of RMS
  TH1D* hcal_internal_timing_rms = new TH1D("hcal_internal_timing_rms", "HCal Internal Timing Resolution;ns;Counts", 60, 0, 6);

  for (auto hist : base_hists) {
    if (hist->GetEntries() > 50) {
      double rms = hist->GetRMS();
      hcal_internal_timing_rms->Fill(rms);
    }
  }

  TCanvas *c1 = new TCanvas("c1", "Internal Timing Resolution", 800, 600);
  c1->cd();
  double jmean, jsigma, jrms;
  //fitHistogram(hcal_internal_timing_res, jmean, jsigma);

  TF1* gausFit = new TF1("gausFit", "gaus",0,2);
  hcal_internal_timing_res->Fit(gausFit, "R");
  jmean = gausFit->GetParameter(1);
  jsigma = gausFit->GetParameter(2);
  jrms = hcal_internal_timing_res->GetRMS();

  hcal_internal_timing_res->SetLineColor(kBlack);
  hcal_internal_timing_res->SetLineWidth(2);
  hcal_internal_timing_res->Draw();

  gPad->SetGridx();
  gPad->SetGridy();

  // Create a legend and add the mean and standard deviation
  TLegend* jlegend = new TLegend(0.16, 0.45, 0.5, 0.65);
  jlegend->SetBorderSize(0); // Remove border
  jlegend->SetFillStyle(0);  // Make background transparent
  jlegend->SetMargin(0);     // Remove left margin
  jlegend->SetTextColor(kRed);
  jlegend->AddEntry((TObject*)0, Form("Mean = %.2f ns", jmean), "");

  jlegend->Draw();
  
  c1->Update();
  c1->Write();

  TCanvas *c2 = new TCanvas("c2", "Heat Map", 800, 600);
  c2->cd();
  heat_map->Draw("COLZ");
  c2->Write();

  // Select four histograms to check the fit quality
  std::vector<int> selected_indices;
  int marker = 0;
  for (size_t i=0; i<base_hists.size(); ++i) {
    if(base_hists[i]->GetEntries()>minev){
      selected_indices.push_back(i);
      marker++;
    }
    if(marker>3)
      break;
  }

  // Create a canvas with 6 columns and 4 rows
  TCanvas *c3 = new TCanvas("c3", "Fit Quality Check", 2400, 1600);
  c3->Divide(6, 4);

  int maxN[24] = {0};

  // Ensure there are 288 histograms, with each element representing a single cell
  for (size_t i = 0; i < base_hists.size(); ++i) {
    int col = i % 12;  // 12 cells per row
    int row = i / 12;

    // Skip plotting if the histogram has too few entries
    if (base_hists[i]->GetEntries() < minev || base_hists[i]->GetEntries() < maxN[row])
      continue;

    maxN[row] = base_hists[i]->GetEntries();

    // Select the appropriate pad in the canvas and plot the histogram
    c3->cd(row + 1);  // +1 because pad numbers start at 1
    gPad->SetGridx();
    gPad->SetGridy();

    base_hists[i]->SetLineColor(kBlack);
    base_hists[i]->SetLineWidth(2);
    base_hists[i]->Draw();
    if (base_hists[i]->GetFunction("gaussFit")) {
      base_hists[i]->GetFunction("gaussFit")->Draw("same");
    }

    double sigma = base_hists[i]->GetFunction("gaussFit")->GetParameter(2);

    // Create a legend and add the mean and standard deviation
    TLegend* legend = new TLegend(0.16, 0.45, 0.5, 0.65);
    legend->SetBorderSize(0); // Remove border
    legend->SetFillStyle(0);  // Make background transparent
    legend->SetMargin(0);     // Remove left margin
    legend->SetTextColor(kRed);
    legend->AddEntry((TObject*)0, Form("Sigma = %.4f ns", sigma), "");

    legend->Draw();

    c3->Update();

  }

  // Write the canvas to the file
  c3->Write();

  
  TCanvas *c4 = new TCanvas("c4", "Internal Timing Resolution RMS", 800, 600);
  c4->cd();
  double imean, isigma, irms;

  TF1* gausFitt = new TF1("gausFitt", "gaus",0,6);
  hcal_internal_timing_rms->Fit(gausFitt, "R");
  imean = gausFitt->GetParameter(1);
  isigma = gausFitt->GetParameter(2);
  irms = hcal_internal_timing_rms->GetRMS();

  hcal_internal_timing_rms->Draw();

  gPad->SetGridx();
  gPad->SetGridy();


  // Create a legend and add the mean and standard deviation
  TLegend* legend = new TLegend(0.16, 0.45, 0.5, 0.65);
  legend->SetBorderSize(0); // Remove border
  legend->SetFillStyle(0);  // Make background transparent
  legend->SetMargin(0);     // Remove left margin
  //legend->SetTextSize(0.05);
  legend->SetTextColor(kRed);
  legend->AddEntry((TObject*)0, Form("Mean = %.2f ns", imean), "");
  //legend->AddEntry((TObject*)0, Form("Sigma = %.4f ns", isigma), "");
  //legend->AddEntry((TObject*)0, Form("RMS = %.4f ns", irms), "");

  legend->Draw();
  
  c4->Update();

  c4->Write();


  std::cout << "All plots created. Output files are located in: " << outdir_path << std::endl;
}
