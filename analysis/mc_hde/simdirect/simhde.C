//SSeeds 4.14.23 - script adapted from monte-carlo TOF to extract HCal detection efficiencies from material properties and detector geometry from g4sbs

#include <vector>
#include <iostream>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "gmn_tree.C"
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

//Add limits for detection efficiency based on active area in hcal (1 block center to center division from edge, x and y)
const Double_t xmin = -1.30125;
const Double_t xmax = 2.19125;
const Double_t ymin = -6.72934;
const Double_t ymax = -5.17994;

//Main, kine->kinematic, field->sbs magnetic field setting in percent
void simhde( Int_t kine = 4, Int_t field = 30 ){
  
  // Define a clock to check overall time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );
  
  // Get the date
  string date = util::getDate();

  //Set processing params
  TChain *C = new TChain("T");
  TFile *fout = new TFile( Form("outfiles/simhde_sbs%d_mag%d.root",kine,field), "RECREATE" );
  Int_t nevent=1;

  C->Add(Form("/lustre19/expphy/volatile/halla/sbs/seeds/sim032123/gmn_sbs%d_ld2_%dp_job*",kine,field));

  //Use root makeclass() functionality, no rush
  gmn_tree *T = new gmn_tree(C);

  if( C->GetEntries()!=0 ){
    cout << "Opened file successfully." << endl;
  }else{
    cout << "Error: No file found." << endl;
    return;
  }

  //Diagnostic histograms
  TH1D *hSDx = new TH1D("hSDx","HCal Elastic Hadron Sensitive Detector Boundary Crossing, x; x_{HCAL} (m)",600,-3.,3.);
  TH1D *hSDy = new TH1D("hSDy","HCal Elastic Hadron Sensitive Detector Boundary Crossing, y; y_{HCAL} (m)",400,-7.5,-5.0);  
  TH2D *hSDxy = new TH2D("hSDxy","Elastic Hadron Sensitive Detector Boundary Crossing, x:y;y_{HCAL} (m); x_{HCAL} (m)", 400, -7.5, -5., 600, -3.0, 3.0 );
  TH1D *hevphit = new TH1D("hevphit","Proton hits per event; N",100,100,100);
  TH1D *hevnhit = new TH1D("hevnhit","Neutron hits per event; N",100,100,100);
  TH1D *hevhitposdiff = new TH1D("hevhitposdiff","Distance between SD BC location and hit position; (m)",200,0.,10.);
  TH1D *hMID = new TH1D("hMID","Mother Particle ID; idx",200,0.,10.);
  TH1D *hsumE = new TH1D("hsumE","Total E Dep HCal; GeV",200,0.,10.);

  //Define loop parameters and counters
  Long64_t Nevents = C->GetEntries();
  Int_t NphitEdep = 0;
  Int_t NprotEdep = 0;
  Int_t NprotSDbc = 0;
  Int_t NnhitEdep = 0;
  Int_t NneutEdep = 0;
  Int_t NneutSDbc = 0;

  Int_t sumEp = 0;
  Int_t sumEn = 0;
  Int_t pev = 0;
  Int_t nev = 0;

  cout << "Opened tree with " << Nevents << " simulated events." << endl;

  while( T->GetEntry( nevent++ ) ){ 
    
    cout << "Processing event: " << nevent << "/" << Nevents << "\r";
    cout.flush();
    
    Int_t prothits = 0;
    Int_t neuthits = 0;

    Double_t sumE = T->Harm_HCalScint_det_esum;
    hsumE->Fill(sumE);

    //extract hde from SD track information and corresponding energy deposited in scintillator
    for( Int_t ihit=0; ihit<T->Harm_HCalScint_hit_nhits; ihit++ ){ //loop over all hits in HCal
      Int_t sdMID = ( *(T->SDTrack_MID) )[( *(T->Harm_HCalScint_hit_sdtridx) )[ ihit ]]; //sd hit particle ID
      hMID->Fill(sdMID);

      //skip tracks which are not made by the primary particle
      if( sdMID!=0 ) continue;

      Double_t sdBCtime = ( *(T->SDTrack_T) )[( *(T->Harm_HCalScint_hit_sdtridx) )[ ihit ]]; //boundary crossing time
      Int_t sdPID = ( *(T->SDTrack_PID) )[( *(T->Harm_HCalScint_hit_sdtridx) )[ ihit ]]; //sd hit particle ID
      Double_t sdhitPosY = ( *(T->SDTrack_posx) )[( *(T->Harm_HCalScint_hit_sdtridx) )[ ihit ]]; //sd hit position Y
      Double_t sdhitPosX = ( *(T->SDTrack_posy) )[( *(T->Harm_HCalScint_hit_sdtridx) )[ ihit ]]; //sd hit position X
      Double_t sdhitpx = ( *(T->SDTrack_momx) )[( *(T->Harm_HCalScint_hit_sdtridx) )[ ihit ]]; //sd hit momentum x
      Double_t sdhitpy = ( *(T->SDTrack_momy) )[( *(T->Harm_HCalScint_hit_sdtridx) )[ ihit ]]; //sd hit momentum y
      Double_t sdhitpz = ( *(T->SDTrack_momz) )[( *(T->Harm_HCalScint_hit_sdtridx) )[ ihit ]]; //sd hit momentum z
      Double_t hitPosY = ( *(T->Harm_HCalScint_hit_xcellg) )[ ihit ]; //hit position Y
      Double_t hitPosX = ( *(T->Harm_HCalScint_hit_ycellg) )[ ihit ]; //hit position X
      Double_t hitEdep = ( *(T->Harm_HCalScint_hit_sumedep) )[ ihit ]; //hit energy deposited in scint

      Int_t cell = (*(T->Harm_HCalScint_hit_cell))[ihit];

      //calculate distance between sensitive detector boundary crossing and hit location
      Double_t hitdiff = sqrt(pow(hitPosX-sdhitPosX,2)+pow(hitPosY-sdhitPosY,2));
      hevhitposdiff->Fill(hitdiff);
      //cout << hitdiff << endl;

      hSDx->Fill(sdhitPosX);
      hSDy->Fill(sdhitPosY);
      hSDxy->Fill(sdhitPosX,sdhitPosY);

      if( sdhitPosX>xmin&&
      	  sdhitPosX<xmax&&
      	  sdhitPosY>ymin&&
      	  sdhitPosY<ymax )
      	{
	  if( sdPID == physconst::IDXp ){ //particle id is proton
	    if( sumE>0 && hitdiff<econst::hcalblk_div_hyp ) NprotEdep++; //"detected" if energy deposited by sd event
	    if( hitEdep>0 ) NphitEdep++;
	    NprotSDbc++; //"expected" if proton crossed within active area
	    prothits++;
	  }else if( sdPID == physconst::IDXn ){ //particle id is neutron
	    if( sumE>0 && hitdiff<econst::hcalblk_div_hyp ) NneutEdep++; //"detected" if energy deposited by sd event
	    if( hitEdep>0 ) NnhitEdep++;
	    NneutSDbc++; //"expected" if proton crossed within active area
	    neuthits++;
	  }
	}
    } //end loop over hits in hcal scintillator
    
    if( prothits>0 ){
      pev++;
      if( sumE>0 ) sumEp++;
    }
    if( neuthits>0 ){
      nev++;
      if( sumE>0 ) sumEn++;
    }
    hevphit->Fill(prothits);
    hevnhit->Fill(neuthits);

  } //end loop over events

  cout << endl << endl;

  cout << "NprotEdep : " << NprotEdep << endl;
  cout << "NneutEdep : " << NneutEdep << endl;
  cout << "NprotSDbc : " << NprotSDbc << endl;
  cout << "NneutSDbc : " << NneutSDbc << endl;
  cout << "NphitEdep : " << NphitEdep << endl;
  cout << "NnhitEdep : " << NnhitEdep << endl;
  cout << "sumEp     : " << sumEp << endl;
  cout << "sumEn     : " << sumEn << endl;
  cout << "pev       : " << pev << endl;
  cout << "nev       : " << nev << endl;


  //Efficiency is "detected" over "expected"
  Double_t hde_prot = ((double)NprotEdep / (double)NprotSDbc )*100.;
  Double_t hde_neut = ((double)NneutEdep / (double)NneutSDbc )*100.;

  cout << "HCal proton detection efficiency (all tracks, loc cut): " << hde_prot << "%" << endl;
  cout << "HCal neutron detection efficiency (all tracks, loc cut): " << hde_neut << "%" << endl; 

  Double_t hde_phit = ((double)NphitEdep / (double)NprotSDbc )*100.;
  Double_t hde_nhit = ((double)NnhitEdep / (double)NneutSDbc )*100.;

  cout << "HCal proton detection efficiency (all tracks): " << hde_phit << "%" << endl;
  cout << "HCal neutron detection efficiency (all tracks): " << hde_phit << "%" << endl;

  Double_t hde_pev = ((double)sumEp / (double)pev )*100.;
  Double_t hde_nev = ((double)sumEn / (double)nev )*100.;

  cout << "HCal proton detection efficiency: " << hde_pev << "%" << endl;
  cout << "HCal neutron detection efficiency: " << hde_nev << "%" << endl; 

  fout->Write();

  // Send time efficiency report to console
  st->Stop();
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
