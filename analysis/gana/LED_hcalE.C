//sseeds 4.15.24, script to check use trigbits 32 (HCal DVCS pulser trigger only) to see hcal E distribution throughout runs and plot the means.
//Better analysis will entail full HCal block readout, this will only return primary cluster spectra.
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

// Main. kine=kinematic
void LED_hcalE(int kine = 7) {
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
  std::string file_path = outdir_path + Form("/gmn_analysis/LED_hcale_out_sbs%d.root",kine);
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
  Double_t tdc, adct, hcale, hcalid, hcalx, hcaly;

  // Switch on only the branches we need
  C->SetMakeClass(1);
  C->SetBranchStatus("*",0);
  C->SetBranchStatus("Event_Branch.fEvtHdr.fTrigBits",1);
  C->SetBranchStatus("Event_Branch.fEvtHdr.fRun",1);
  C->SetBranchStatus("sbs.hcal.tdctimeblk",1);
  C->SetBranchStatus("sbs.hcal.atimeblk",1);
  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);

  // Branches are of type UInt_t
  C->SetBranchAddress("Event_Branch.fEvtHdr.fTrigBits", &trigBits);
  C->SetBranchAddress("Event_Branch.fEvtHdr.fRun", &runNumber);
  C->SetBranchAddress("sbs.hcal.tdctimeblk", &tdc);
  C->SetBranchAddress("sbs.hcal.atimeblk", &adct);
  C->SetBranchAddress("sbs.hcal.e", &hcale);
  C->SetBranchAddress("sbs.hcal.x", &hcalx);
  C->SetBranchAddress("sbs.hcal.y", &hcaly);

  // Create histograms
  std::map<int, TH1D*> histograms_hcale;
  std::map<int, TH2D*> histograms_hcalpos;

  TH1D *hcale_ol = new TH1D("hhcale_ol", "HCal PClus E LED Trig Only Overall;Energy (GeV);Entries", 200, 0., 1.0);
  TH2D *hcalpos_ol = new TH2D("hhcalpos_ol", "HCal PClus Pos LED Trig Only Overall;Energy (GeV);Entries", 400, -2, 2, 600, -3, 2);

  // Loop over all entries in the chain, get the run number, and add a TH1D for hcale where trigbits==32 exists
  std::cout << "Looping over entries in chain.." << std::endl;
  Long64_t nentries = C->GetEntries();
  for (Long64_t i = 0; i < nentries; ++i) {
    C->GetEntry(i);

    cout << "working on run " << runNumber << " event " << i << "/" << nentries << "\r";
    cout.flush();

    if (trigBits != 32) continue; //skip if the trigger wasn't DVCS pulser alone

    cout << endl << "LED event!" << endl;

    if (histograms_hcale.find(runNumber) == histograms_hcale.end()) {
      histograms_hcale[runNumber] = new TH1D(Form("hhcale_%d", runNumber), Form("HCal PClus E LED Trig Only for Run %d;Energy (GeV);Entries", runNumber), 200, 0., 1.0);
    }

    if (histograms_hcalpos.find(runNumber) == histograms_hcalpos.end()) {
      histograms_hcalpos[runNumber] = new TH2D(Form("hhcalpos_%d", runNumber), Form("HCal PClus Pos LED Trig Only for Run %d;Energy (GeV);Entries", runNumber), 400, -2, 2, 600, -3, 2);
    }

    histograms_hcale[runNumber]->Fill(hcale);
    hcale_ol->Fill(hcale);
    histograms_hcalpos[runNumber]->Fill(hcaly,hcalx);
    hcalpos_ol->Fill(hcaly,hcalx);
  }
			      
  // Save histograms
  for (auto& hist : histograms_hcale) {
    hist.second->Write();
    delete hist.second;
  }
  for (auto& hist : histograms_hcalpos) {
    hist.second->Write();
    delete hist.second;
  }

  hcale_ol->Write();
  hcalpos_ol->Write();

  fout->Close();

  cout << "Analysis complete and written to " << file_path << endl;

}
