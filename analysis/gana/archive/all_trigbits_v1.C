//sseeds 4.15.24, script to check trigbits on all files to examine available trigger information on pass2 replays
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

void all_trigbits(int kine = 7) {
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
    directories = {"/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take4/SBS11/LD2/rootfiles/part1",
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take4/SBS11/LH2/rootfiles/part1",
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take4/SBS11/LD2/rootfiles/part2",
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take4/SBS11/LH2/rootfiles/part2",
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take4/SBS11/LD2/rootfiles/part3",
		   "/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass2_take4/SBS11/LH2/rootfiles/part3"};
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
  UInt_t trigBits;
  UInt_t runNumber;
  double tdc;
  double adct;

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
  std::set<int> uniqueRuns;  // To track unique run numbers
  std::set<int> processedRuns;  // To track processed run numbers

  // Loop over all entries in the chain
  std::cout << "Looping over entries in chain.." << std::endl;
  Long64_t nentries = C->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    C->GetEntry(i);

    runToTrigBits[runNumber].insert(trigBits);
    uniqueRuns.insert(runNumber);  // Track unique run numbers

    // Check if the run number has been printed before
    if (processedRuns.find(runNumber) == processedRuns.end()) {
      std::cout << "Processing New Run Number: " << runNumber << std::endl;
      std::cout.flush();
      processedRuns.insert(runNumber);
    }

  }

  // Output the total number of unique run numbers
  std::cout << "Total unique run numbers encountered: " << uniqueRuns.size() << std::endl;

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

  // Draw the graph
  TCanvas* c1 = new TCanvas("c1", "Data Analysis", 1600, 600);
  c1->SetGrid(1,1);
  graph->Draw("AP");
}
