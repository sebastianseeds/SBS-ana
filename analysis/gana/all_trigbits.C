//sseeds 4.15.24, script to check trigbits on all files to examine available trigger information on pass2 replays. Script also asses the TDC and ADCt windows and minima/maxima per run
#include <TFile.h>
#include <TSystem.h>
#include <TChain.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <iostream>
#include <TString.h>
#include <TChain.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>

const double tdc_ul = 300.;
const double tdc_ll = -300.;
const double adct_ul = 200.;
const double adct_ll = -200.;

// Forward Declarations
void writeDataToFile(const std::map<int, std::set<int>>& runToTrigBits, int kine);
void writeDataCtsToFile(const std::map<int, std::set<int>>& runToTrigBits,
			const std::map<int, std::map<int, int>>& runToTrigBitCounts,
			int kine);
// Main. kine=kinematic
void all_trigbits(int kine = 9) {
  std::vector<std::string> directories;

  switch(kine) {
  case 4:
    directories = {"/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take4/SBS4/LD2/rootfiles",
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take4/SBS4/LH2/rootfiles"};
    break;
  case 7:
    directories = {"/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take3/SBS7/LD2/rootfiles",
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take3/SBS7/LH2/rootfiles"};
    break;
  case 14:
    directories = {"/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take4/SBS14/LD2/rootfiles",
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take4/SBS14/LH2/rootfiles"};
    break;
  case 9:
    directories = {"/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take3/SBS9/LD2/rootfiles",
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take3/SBS9/LH2/rootfiles"};
    break;
  case 11:
    directories = {"/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take4/SBS11/LD2/part1/rootfiles",
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take4/SBS11/LD2/part2/rootfiles",
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take4/SBS11/LD2/part3/rootfiles",
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take4/SBS11/LH2/rootfiles"};
    break;
  case 8:
    directories = {"/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take3/SBS8/LD2/rootfiles/SBS0percent", 
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take3/SBS8/LD2/rootfiles/SBS100percent", 
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take3/SBS8/LD2/rootfiles/SBS50percent", 
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take3/SBS8/LD2/rootfiles/SBS70percent_part1", 
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take3/SBS8/LD2/rootfiles/SBS70percent_part2", 
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take3/SBS8/LD2/rootfiles/SBS70percent_part3", 
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take3/SBS8/LD2/rootfiles/SBS70percent_part4", 
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take3/SBS8/LH2/rootfiles/SBS0percent", 
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take3/SBS8/LH2/rootfiles/SBS100percent", 
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take3/SBS8/LH2/rootfiles/SBS50percent", 
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take3/SBS8/LH2/rootfiles/SBS70percent_part1", 
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take3/SBS8/LH2/rootfiles/SBS70percent_part2", 
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take3/SBS8/LH2/rootfiles/SBS70percent_part3"};
    break;
  default:
    std::cout << "ERROR: Enter a valid kinematic setting." << std::endl;
    return;
  }

  // Outfile path
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string file_path = outdir_path + Form("/gmn_analysis/trigbits_and_timingwindow_out_sbs%d.root",kine);
  TFile *fout = new TFile( file_path.c_str(), "RECREATE" ); 

  TChain *C = new TChain("T");

  // Loop over directories and add files to the chain
  std::cout << "Adding SBS-" << kine << " files to chain.." << std::endl;

  for (const auto& dir : directories) {
    TString pattern = TString::Format("%s/*.root", dir.c_str());
    C->Add(pattern);
  }

  //C->Add("/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take4/SBS4/LD2/rootfiles/e1209019_fullreplay_11595*.root");

  std::cout << "Setting up all needed branches.." << std::endl;

  // Set up variables to hold data from branches
  UInt_t trigBits, runNumber;
  Double_t tdc, adct;

  // Switch on only the branches we need
  C->SetMakeClass(1);
  C->SetBranchStatus("*",0);
  C->SetBranchStatus("Event_Branch.fEvtHdr.fTrigBits",1);
  C->SetBranchStatus("Event_Branch.fEvtHdr.fRun",1);
  C->SetBranchStatus("sbs.hcal.tdctimeblk",1);
  C->SetBranchStatus("sbs.hcal.atimeblk",1);

  // Branches are of type UInt_t
  C->SetBranchAddress("Event_Branch.fEvtHdr.fTrigBits", &trigBits);
  C->SetBranchAddress("Event_Branch.fEvtHdr.fRun", &runNumber);
  C->SetBranchAddress("sbs.hcal.tdctimeblk", &tdc);
  C->SetBranchAddress("sbs.hcal.atimeblk", &adct);

  std::map<int, std::set<int>> runToTrigBits;
  std::map<int, std::map<int, int>> runToTrigBitCounts;
  std::set<int> uniqueRuns;  // To track unique run numbers
  std::set<int> processedRuns;  // To track processed run numbers

  std::map<int, std::pair<double, double>> runToTDCRange; // map to store min and max TDC values per run
  std::map<int, std::pair<double, double>> runToADCRange; // map to store min and max ADC values per run

  // Loop over all entries in the chain
  std::cout << "Looping over entries in chain.." << std::endl;
  Long64_t nentries = C->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {

    std::cout << "Processing event " << i << "/" << nentries << "\r";
    std::cout.flush();

    C->GetEntry(i);

    runToTrigBits[runNumber].insert(trigBits);
    uniqueRuns.insert(runNumber);  // Track unique run numbers
    runToTrigBitCounts[runNumber][trigBits]++; // Increment the count for this trigbit

    // Check if the run number has been printed before
    if (processedRuns.find(runNumber) == processedRuns.end()) {
      std::cout << std::endl << "Processing New Run Number: " << runNumber << std::endl;
      std::cout.flush();
      processedRuns.insert(runNumber);
    }

    // Update the TDC and ADC ranges per run
    if (runToTDCRange.find(runNumber) == runToTDCRange.end()) {
      if(tdc<tdc_ll||tdc>tdc_ul)
	runToTDCRange[runNumber] = {0, 0};
      else
	runToTDCRange[runNumber] = {tdc, tdc}; // Initialize min and max to the current tdc value
    } else {
      if(tdc<tdc_ll||tdc>tdc_ul) // Eliminate values for the TDC written where no data exists
	continue;
      runToTDCRange[runNumber].first = std::min(runToTDCRange[runNumber].first, tdc);
      runToTDCRange[runNumber].second = std::max(runToTDCRange[runNumber].second, tdc);
    }

    if (runToADCRange.find(runNumber) == runToADCRange.end()) {
      if(adct<adct_ll||adct>adct_ul)
	runToADCRange[runNumber] = {0, 0}; // Initialize min and max to the current adct value
      else
	runToADCRange[runNumber] = {adct, adct}; // Initialize min and max to the current adct value
    } else {
      if(adct<adct_ll||adct>adct_ul)
	continue;
      runToADCRange[runNumber].first = std::min(runToADCRange[runNumber].first, adct);
      runToADCRange[runNumber].second = std::max(runToADCRange[runNumber].second, adct);
    }

  }

  // Output the total number of unique run numbers
  std::cout << "Total unique run numbers encountered: " << uniqueRuns.size() << std::endl;

  // Prepare variables to store the differences
  std::vector<double> tdcDiff, adctDiff;

  // Calculate differences and prepare for graphing
  for (const auto& pair : runToTDCRange) {
    double diff = fabs(pair.second.second - pair.second.first); // Absolute difference for TDC
    tdcDiff.push_back(diff);
  }

  for (const auto& pair : runToADCRange) {
    double diff = fabs(pair.second.second - pair.second.first); // Absolute difference for ADCT
    adctDiff.push_back(diff);
  }

  // Draw the graph
  TCanvas* c1 = new TCanvas("c1", "Data Analysis", 1600, 600);
  c1->SetGrid(1,1);
  c1->cd();

  // Now prepare to plot the data using TGraph
  std::vector<int> x, y;
  for (const auto& pair : runToTrigBits) {
    for (int bit : pair.second) {
      x.push_back(pair.first);
      y.push_back(bit);
    }
  }

  // Create the graph
  TGraph* graph = new TGraph(x.size(), x.data(), y.data());
  graph->SetTitle(Form("Trigger Bits by Run, SBS-%d;Run Number;Trig Bits",kine));
  graph->SetMarkerStyle(29);
  graph->GetYaxis()->SetRangeUser(0,40);
  graph->Draw("AP");
  c1->Update();
  c1->Write();

  // Write a text file to record these trigbits as well
  writeDataToFile(runToTrigBits,kine);
  writeDataCtsToFile(runToTrigBits,runToTrigBitCounts,kine);

  // Display the graphs on a canvas
  TCanvas* c2 = new TCanvas("c2", "TDC and ADCT Analysis", 1600, 800);
  c2->Divide(2, 1);
  c2->cd(1);

  // Plot TGraphErrors for TDC
  std::vector<int> runs;
  std::vector<double> tdcMeans, tdcErrors, adctMeans, adctErrors;

  for (const auto& pair : runToTDCRange) {
    runs.push_back(pair.first);
    double mean = (pair.second.first + pair.second.second) / 2;
    double error = (pair.second.second - pair.second.first) / 2;
    tdcMeans.push_back(mean);
    tdcErrors.push_back(error);
  }

  // Convert 'runs' from int to double
  std::vector<double> runs_d;
  for (int run : runs) {
    runs_d.push_back(static_cast<double>(run));
  }

  TGraphErrors* grTDC = new TGraphErrors(runs.size(), runs_d.data(), tdcMeans.data(), nullptr, tdcErrors.data());

  // Configure the graphs (titles, marker styles etc.) and draw them
  grTDC->SetTitle("TDC vs Run;Run Number;TDC");
  grTDC->SetMarkerStyle(0);
  grTDC->Draw("AP");
  c2->cd(2);

  // Plot TGraphErrors for ADCt
  for (const auto& pair : runToADCRange) {
    adctMeans.push_back((pair.second.first + pair.second.second) / 2);
    adctErrors.push_back((pair.second.second - pair.second.first) / 2);
  }

  TGraphErrors* grADCT = new TGraphErrors(runs.size(), runs_d.data(), adctMeans.data(), nullptr, adctErrors.data());

  grADCT->SetTitle("ADCT vs Run;Run Number;ADCT");
  grADCT->SetMarkerStyle(0);

  grADCT->Draw("AP");
  c2->Update();
  c2->Write();

  // Plot the TDC and ADCt differences
  TCanvas* c3 = new TCanvas("c3", "Differences in Measurements", 1600, 800);
  c3->Divide(2, 1);
  c3->cd(1);

  // Plot TDC differences
  TGraph* graphTDCDiff = new TGraph(runs.size(), runs_d.data(), tdcDiff.data());

  // Set titles and styles for the new graphs
  graphTDCDiff->SetTitle("Absolute Difference in TDC vs Run;Run Number;|TDC Max - TDC Min|");
  graphTDCDiff->SetMarkerStyle(29);
  graphTDCDiff->SetMarkerColor(kBlue);
  graphTDCDiff->Draw("AP");
  c3->cd(2);

  // Plot ADCt differences
  TGraph* graphADCTDiff = new TGraph(runs.size(), runs_d.data(), adctDiff.data());

  // Set titles and styles for the new graphs
  graphADCTDiff->SetTitle("Absolute Difference in ADCT vs Run;Run Number;|ADCT Max - ADCT Min|");
  graphADCTDiff->SetMarkerStyle(29);
  graphADCTDiff->SetMarkerColor(kRed);
  graphADCTDiff->Draw("AP");
  c3->Update();
  c3->Write();

  // Plot the TDC and ADCt differences

  // Find min and max for TDC
  double tdcMin = *std::min_element(tdcDiff.begin(), tdcDiff.end());
  double tdcMax = *std::max_element(tdcDiff.begin(), tdcDiff.end());
  
  // Find min and max for ADC
  double adcMin = *std::min_element(adctDiff.begin(), adctDiff.end());
  double adcMax = *std::max_element(adctDiff.begin(), adctDiff.end());

  // Calculate number of bins using Sturges' formula
  int tdcBins = std::ceil(std::log2(tdcDiff.size()) + 1);
  int adcBins = std::ceil(std::log2(adctDiff.size()) + 1);

  TCanvas* c4 = new TCanvas("c4", "Measured Window Distribution", 1600, 800);
  c4->Divide(2, 1);
  c4->cd(1);

  TH1D *TDCdiff = new TH1D("TDCdiff","TDCdiff",tdcBins,tdcMin,tdcMax);
  for( size_t i=0; i<tdcDiff.size(); ++i )
    TDCdiff->Fill(tdcDiff[i]);

  TDCdiff->SetTitle("Absolute Difference in TDC (Measured TDC Window);ns");
  TDCdiff->Draw();
  c4->cd(2);

  TH1D *ADCtdiff = new TH1D("ADCtdiff","ADCtdiff",adcBins,adcMin,adcMax);
  for( size_t i=0; i<adctDiff.size(); ++i )
    ADCtdiff->Fill(adctDiff[i]);

  ADCtdiff->SetTitle("Absolute Difference in ADC Time (Measured ADCt Window);ns");
  ADCtdiff->Draw();
  c4->Update();
  c4->Write();

  fout->Write();

}

void writeDataToFile(const std::map<int, std::set<int>>& runToTrigBits, int kine) {
  std::ofstream outFile(Form("trigbits_per_run/trigger_bits_by_run_sbs%d.txt",kine));  // Create and open the output file

  if (!outFile.is_open()) {
    std::cerr << "Failed to open the file for writing." << std::endl;
    return;
  }

  // Iterate over each pair in the map
  for (const auto& pair : runToTrigBits) {
    int runNumber = pair.first;
    const std::set<int>& trigBits = pair.second;

    outFile << runNumber;  // Write the run number to the file

    for (int bit : trigBits) {
      outFile << ", " << bit;  // Write each trigger bit, preceded by a comma
    }
        
    outFile << std::endl;  // End the line after all trigger bits for a run have been written
  }

  outFile.close();  // Close the file
}

void writeDataCtsToFile(const std::map<int, std::set<int>>& runToTrigBits,
			const std::map<int, std::map<int, int>>& runToTrigBitCounts,
			int kine) {
  std::ofstream outFile(Form("trigbits_per_run/trigger_bits_cts_by_run_sbs%d.txt", kine)); // Create and open the output file

  if (!outFile.is_open()) {
    std::cerr << "Failed to open the file for writing." << std::endl;
    return;
  }

  outFile << "## <run number>, <fEvtHdr.fTrigBits>, <N Events>" << std::endl;

  // Iterate over each pair in the runToTrigBits map
  for (const auto& pair : runToTrigBits) {
    int runNumber = pair.first;
    const std::set<int>& trigBits = pair.second;

    // Check if we have the count data for this run
    auto countIt = runToTrigBitCounts.find(runNumber);
    if (countIt != runToTrigBitCounts.end()) {
      // We have counts for this run, write them to the file
      const std::map<int, int>& trigBitCounts = countIt->second;

      for (int bit : trigBits) {
	outFile << runNumber << ", " << bit << ", " << trigBitCounts.at(bit) << std::endl; // Write run number, trig bit, and number of events
      }
    } else {
      // No counts found for this run, could log an error or write out just the run number and trig bits
      for (int bit : trigBits) {
	outFile << runNumber << ", " << bit << ", " << "N/A" << std::endl; // Write run number, trig bit, and N/A for the counts
      }
    }
  }

  outFile.close(); // Close the file
}
