//sseeds - makes several dx plots with all cuts equalizing statistics, but keeping run numbers separate.
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include <TLegend.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <string>
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

//main. kine=kinematic, mag=magnetic field (percent), pass=pass#, bestclus=use plot file with best cluster selection, Nsamps=number of total jacknife samples
void jacknifeSamples(int kine=4, int mag=50, int pass=2, bool bestclus=true, int Nsamps=10 ) {

  // Configuration and initial setup (assuming you have this from your original code)
  JSONManager *jmgr = new JSONManager("../../config/syst.json");
  int hbins = jmgr->GetValueFromSubKey<int>("hbins", Form("sbs%d",kine));
  double hcalfit_l = jmgr->GetValueFromSubKey<double>("hcalfit_l", Form("sbs%d",kine));
  double hcalfit_h = jmgr->GetValueFromSubKey<double>("hcalfit_h", Form("sbs%d",kine));

  //Get tight elastic cuts
  std::string globalcuts = jmgr->GetValueFromSubKey_str( Form("post_cuts_p%d",pass), Form("sbs%d_%d",kine,mag) );
  cout << "Loaded tight cuts: " << globalcuts << endl;
  TCut globalcut = globalcuts.c_str();

  // Output paths and files
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fin_path = outdir_path + Form("/parse/parse_sbs%d_pass%d_barebones.root", kine, pass);
  std::string fout_path = outdir_path + Form("/gmn_analysis/jacknife_samples_sbs%d_mag%d_pass%d.root", kine, mag, pass);

  TFile* inputFile = new TFile(fin_path.c_str(), "READ");
  TTree* tree = dynamic_cast<TTree*>(inputFile->Get("P"));
  if (!tree) {
    std::cerr << "Tree not found in file: " << fin_path << std::endl;
    inputFile->Close();
    return;
  }

  // Open an output ROOT file
  TFile* outputFile = new TFile(fout_path.c_str(), "RECREATE");

  // Determine the number of runs and events
  int run_number;
  double delta_x;
  tree->SetBranchAddress("run", &run_number);
  tree->SetBranchAddress("dx", &delta_x);

  // Put cuts on tree
  TEventList *elist = new TEventList("elist","Elastic Event List");
  tree->Draw(">>elist",globalcut);

  std::map<int, int> run_event_count;
  Long64_t nentries = elist->GetN();
  for (Long64_t i=0; i < nentries; i++) {
    tree->GetEntry( elist->GetEntry(i) );
    run_event_count[run_number]++;
  }

  // Output the run numbers and their corresponding event counts
  std::cout << std::endl;
  for (const auto& pair : run_event_count) {
    std::cout << "Run Number: " << pair.first << ", Total Events: " << pair.second << std::endl;
  }

  int num_events_per_sample = nentries / Nsamps;

  std::cout << std::endl << "Number of Events per Sample: " << num_events_per_sample << std::endl;

  // Creating histograms for each segment
  std::vector<TH1D*> histograms_stat;
  int current_events = 0;
  int current_run = -1;
  int hist_count = 0;
  bool start_new_hist = true;  // Flag to start a new histogram

  // Equalize statistics for statisical analysis
  for (Long64_t i = 0; i < nentries; i++) {
    tree->GetEntry( elist->GetEntry(i) ); 

    // Check if we need to start a new histogram
    if (start_new_hist) {
      // Create a new histogram for a new run or because the previous one is full
      hist_count++;
      std::cout << "Run-to-run statistical analysis. Starting new histogram " << hist_count << " for run number " << run_number << std::endl;
      histograms_stat.push_back(new TH1D(Form("hdx_stat_%d", hist_count), "dx; m", hbins, hcalfit_l, hcalfit_h));
      current_events = 0;
      start_new_hist = false;
    }

    // Add the event to the current histogram
    histograms_stat.back()->Fill(delta_x);
    current_events++;
    current_run = run_number;

    // Check if the current histogram has reached or exceeded the desired number of events
    if (current_events >= num_events_per_sample) {
      start_new_hist = true;  // Set the flag to start a new histogram at the next run change
    }
  }

  // Creating histograms for each systematic study segment
  std::vector<TH1D*> histograms_syst;
  current_events = 0;
  current_run = -1;

  // Do not combine experimental runs for systematic analysis
  for (Long64_t i=0; i < nentries; i++) {
    tree->GetEntry( elist->GetEntry(i) );    

    if (current_run != run_number) {
      // Save and reset histogram when changing run or exceeding event limit
      std::cout << "Run-to-run systematics. Run number: " << run_number << std::endl;

      histograms_syst.push_back(new TH1D(Form("hdx_syst_%d", run_number), "dx;m", hbins, hcalfit_l, hcalfit_h));
      current_events = 0;
      current_run = run_number;
    }
    histograms_syst.back()->Fill(delta_x);
    current_events++;
  }

  for (auto& hist : histograms_stat) {
    hist->Write(); // Make sure each histogram is written if not already written.
    delete hist;
  }
  histograms_stat.clear();

  for (auto& hist : histograms_syst) {
    hist->Write(); // Make sure each histogram is written if not already written.
    delete hist;
  }
  histograms_syst.clear();

  outputFile->Close(); // Make sure to close the file.
  inputFile->Close();

  std::cout << std::endl << "Analysis complete. Output file located " << fout_path << std::endl;

}
