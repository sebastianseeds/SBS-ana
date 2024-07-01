//sseeds 04.20.23 - Script to extract single-channel best possible timing resolutions by kinematic/field

#include <vector>
#include <iostream>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include <TF1Convolution.h>
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

const Double_t atime_approx_FWHM = 3; //HCal adctime approximate full-width-half-max of elastic peak
const Double_t tdc_approx_FWHM = 3; //HCal tdctime approximate full-width-half-max of elastic peak

void overlayWithGaussianFits(TH2D *h2d, TCanvas *canvas) {
  // switch to canvas location
  canvas->cd();

  // Create a TGraphErrors to hold the means and standard deviations
  TGraphErrors *graphErrors = new TGraphErrors(h2d->GetNbinsX());

  // Loop over the bins in the TH2D
  for (int i = 1; i <= h2d->GetNbinsX(); ++i) {
    // Get the projection of this bin along the y-axis
    TH1D *projY = h2d->ProjectionY("_py", i, i);

    //Skip the fit if less than 100 entries exist in the bin
    if( projY->GetEntries()<100 ){
      //graphErrors->SetPoint(i - 1, h2d->GetXaxis()->GetBinCenter(i), 0);
      //graphErrors->SetPointError(i - 1, 0, 0); // Assuming no error on x
      continue;
    }

    //declare some dynamic fit variables
    Int_t binMax = projY->GetMaximumBin();
    Double_t binCenter = projY->GetBinCenter( binMax );
    Double_t fitLowerLim = binCenter - atime_approx_FWHM;
    Double_t fitUpperLim = binCenter + atime_approx_FWHM;

    // Fit a Gaussian to this projection
    TF1 *gausFit = new TF1("gausFit", "gaus");
    projY->Fit(gausFit, "Q", "", fitLowerLim, fitUpperLim ); // "Q" for quiet mode

    // Get the mean and standard deviation from the fit
    double mean = gausFit->GetParameter(1);
    double stdDev = gausFit->GetParameter(2);

    // Get the error on the mean and standard deviation
    double meanError = gausFit->GetParError(1);
    double stdDevError = gausFit->GetParError(2);

    // Set this point in the TGraphErrors
    graphErrors->SetPoint(i - 1, h2d->GetXaxis()->GetBinCenter(i), mean);
    graphErrors->SetPointError(i - 1, 0, stdDev); // Assuming no error on x
  }

  // Draw the TH2D
  h2d->Draw("COLZ");

  // Customize the appearance of the TGraphErrors
  graphErrors->SetMarkerStyle(20);
  graphErrors->SetMarkerColor(kBlack);
  graphErrors->SetLineColor(kBlack);
  //graphErrors->SetLineStyle(0);
  //graphErrors->SetEndErrorSize(0);

  // Draw the TGraphErrors on the same canvas
  graphErrors->Draw("P SAME");

  // Update the canvas
  gPad->Update();
}


//MAIN kine is kinematic, pass is replay pass, epm is e' momentum calculation method
void hclusblk_res( Int_t kine=8, Int_t epm=3, Int_t pass=1 )
{ //main

  Double_t adcbinw = econst::hcaladc_binw; //adc sample bin width (4*40=160ns over whole waveform)

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // reading config file
  JSONManager *jmgr = new JSONManager("../../config/shtres.json");
  std::string rootfile_dir;
  rootfile_dir = jmgr->GetValueFromSubKey_str("general_rootfile_dir",Form("sbs%d",kine));

  cout << rootfile_dir << endl;

  //define structure to encode parameters for continuous ep fit tdc/adct corrections
  typedef struct{
    vector<Int_t> kIdx;
    Int_t nfset_lh2; //total number of different field settings for lh2
    Int_t nfset_ld2; //total number of different field settings for ld2
    vector<Int_t> fset_lh2; //all field settings for lh2 (percent, -1 indicates no setting)
    vector<Int_t> fset_ld2; //all field settings for ld2 (percent, -1 indicates no setting)
    Int_t TOFfitp_mset; //field settings with TOF vs nucleon p available
    vector<Double_t> TOFfitp_p; //TOF vs proton p fit parameters
    vector<Double_t> TOFfitp_n; //TOF vs neutron p fit parameters
    Double_t TOFfitp_ulim; //upper nucleon momentum limit on fits, beyond which fits diverge quickly
  } SBSSET;

  SBSSET fitparams[econst::nkine]; //One structure per kinematic

  Int_t gIdx; //get the global index for this kinematic
  for( Int_t k=0; k<econst::nkine; k++ ){
    jmgr->GetVectorFromKey<Int_t>("kIdx",fitparams[k].kIdx);
    Int_t ki = fitparams[k].kIdx[k];

    fitparams[k].nfset_lh2 = jmgr->GetValueFromSubKey<Int_t>("nfset_lh2",Form("sbs%d",ki));    
    fitparams[k].nfset_ld2 = jmgr->GetValueFromSubKey<Int_t>("nfset_ld2",Form("sbs%d",ki));
    jmgr->GetVectorFromSubKey<Int_t>("fset_lh2",Form("sbs%d",ki),fitparams[k].fset_lh2);
    jmgr->GetVectorFromSubKey<Int_t>("fset_ld2",Form("sbs%d",ki),fitparams[k].fset_ld2);
    fitparams[k].TOFfitp_mset = jmgr->GetValueFromSubKey<Int_t>("TOFfitp_mset",Form("sbs%d",ki));
    jmgr->GetVectorFromSubKey<Double_t>("TOFfitp_p",Form("sbs%d",ki),fitparams[k].TOFfitp_p);
    jmgr->GetVectorFromSubKey<Double_t>("TOFfitp_n",Form("sbs%d",ki),fitparams[k].TOFfitp_n);
    fitparams[k].TOFfitp_ulim = jmgr->GetValueFromSubKey<Double_t>("TOFfitp_ulim",Form("sbs%d",ki));

    if( fitparams[0].kIdx[k]==kine )
      gIdx=k;
  }

  cout << fitparams[4].nfset_lh2 << endl;

  //Read in timewalk parameters
  typedef struct{
    vector<Double_t> tdcp;
    vector<Double_t> adcp;
  } TW;

  TW timewalk[3]; //three timewalk parameters of the form p0*exp(-p1*x[0])+p2

  for( Int_t p=0; p<3; p++ ){
    std::string tdcpath = Form("params/TW/TDC/tdctwP%d_sbs%d.txt",p,kine);
    util::readParam(tdcpath,timewalk[p].tdcp);
    std::string adctpath = Form("params/TW/ADC/adcttwP%d_sbs%d.txt",p,kine);
    util::readParam(adctpath,timewalk[p].adcp);
  }

  //Read in gain parameters for more accurate timewalk corrections
  vector<Double_t> oldadcgain;
  vector<Double_t> adcgain;
  std::string oldadcgainpath = "params/gain/oldcoeff.txt";
  std::string adcgainpath = Form("params/gain/coeff_sbs%d.txt",kine);
  std::string gmoniker = "sbs.hcal.adc.gain";
  util::readParam(oldadcgainpath,gmoniker,econst::hcalchan,oldadcgain);
  util::readParam(adcgainpath,gmoniker,econst::hcalchan,adcgain);

  vector<Int_t> lh2r; //jmgr adds (not replace) when called
  jmgr->GetVectorFromSubKey<Int_t>( "lh2runs", "sbs4_runs", lh2r );
  jmgr->GetVectorFromSubKey<Int_t>( "lh2runs", "sbs9_runs", lh2r );

  vector<Int_t> ld2r;
  jmgr->GetVectorFromSubKey<Int_t>( "ld2runs", "sbs8_runs", ld2r );
  jmgr->GetVectorFromSubKey<Int_t>( "ld2runs", "sbs11_runs", ld2r );

  //set up default parameters for all analysis
  std::string runsheet_dir = "/w/halla-scshelf2102/sbs/seeds/ana/data"; //unique to my environment for now
  Int_t nruns = -1; //Always analyze all available runs for this application. If not, need to pass nruns for lh2 and ld2 separately.
  Int_t nruns_h = -1;
  Int_t nruns_d = -1;
  Int_t verb = 0; //Don't print diagnostic info by default

  //read the run list to parse run numbers associated to input parameters and set up for loop over runs
  vector<crun> crunh; 
  vector<crun> crund; 
  util::ReadRunList(runsheet_dir,nruns_h,kine,"LH2",pass,verb,crunh); //modifies nruns to be very large when -1
  util::ReadRunList(runsheet_dir,nruns_d,kine,"LD2",pass,verb,crund); //modifies nruns to be very large when -1
  
  //set up output files

  std::string test_file = "test.root";

  //TFile *fout = new TFile( Form("outfiles/test_htresout_sbs%d_epm%d.root",kine,epm), "RECREATE" );
  //TFile *fout = new TFile( Form("/lustre19/expphy/volatile/halla/sbs/seeds/hcal_internal_timing/htresout_sbs%d_epm%d.root",kine,epm), "RECREATE" );

  TFile *fout = new TFile( test_file.c_str(), "RECREATE" );

  // re-allocate memory at each run to load different cuts/parameters
  TChain *C = nullptr;

  //Set up analysis histograms and arrays
  Int_t Nevents[288][9]={0};
  Int_t Neventstdc[288][9]={0};

  TH1D *hdx_HCAL = new TH1D("hdx_HCAL ", " ; x_{HCAL} - x_{exp} (m)  ", 250,-4,3);
  TH1D *hdy_HCAL = new TH1D("hdy_HCAL","; y_{HCAL} - y_{exp} (m);", 250, -1.25,1.25);
  TH1D *hW2 = new TH1D("hW2", " ;GeV2  ", 500,0,5);
  TH1D *hpse = new TH1D("hpse"," ;GeV", 200, 0, 2);
  
  TH2D *hadiff_id[9];
  TH2D *htdiff_id[9];

  Int_t Nevents_gtdc[288][9]={0};

  TH1D *hpp = new TH1D("hpp",
		       "Expected Scattered Proton Momentum; GeV",
		       500,
		       0,
		       5);

  
  TH1D *hpblke = new TH1D("hpblke",
		       "Primary block E; GeV",
		       200,
		       0,
		       0.5);

  
  TH1D *hsblke = new TH1D("hsblke",
		       "Secondary block E; GeV",
		       200,
		       0,
		       0.5);

  TH1D *hadctsig = new TH1D("hadctsig",
			    "ADC time sigma, channel doubles > 500 events; ns",
			    100,
			    0,
			    2);

  TH1D *htdcsig = new TH1D("htdcsig",
			    "TDC sig, chan pairs > 250 events, min block E > 30 MeV, tw corr; ns",
			    100,
			    0,
			    2);

  TH1D *hgtdc = new TH1D("hgtdc",
			 "TDC time pblk-blk diff, chan 126, idx 1; ns",
			 100,
			 -5,
			 5);

  TH2D *hadct_E = new TH2D("hadct_E",
			      "ADC time (all channels, primary block) vs E",
			      100,
			      0,
			      0.5,
			      280,
			      -20,
			      120);

  TH2D *htdc_E = new TH2D("htdc_E",
			      "TDC time (primary cluster channels, primary block) vs E",
			      100,
			      0,
			      0.5,
			      280,
			      -20,
			      120);


  TH2D *htdc_E_corr = new TH2D("htdc_E_corr",
			      "TDC time corrected (primary cluster channels, primary block) vs E",
			      100,
			      0,
			      0.5,
			      280,
			      -20,
			      120);
  

  TH2D *hElastic_doubles[9];
  TH1D *hadcTdiff_doubles[288][9];
  TH1D *htdcdiff_doubles[288][9];

  for( Int_t n=0; n<9; n++ ){
    hElastic_doubles[n] = new TH2D(Form("hElastic_doubles_%d",n),
				   Form("Number Elastic Events, Nidx:%d; col; row",n),
				   12,
				   0,
				   12,
				   24,
				   0,
				   24);

    hadiff_id[n] = new TH2D(Form("hadiff_id_i%d",n),
			    Form("ADC time diff vs ID, Nidx:%d; channel; ns",n),
			    econst::hcalchan,
			    0,
			    econst::hcalchan,
			    160,
			    -20,
			    20);

    htdiff_id[n] = new TH2D(Form("htdiff_id_i%d",n),
			    Form("TDC diff vs ID, Nidx:%d; channel; ns",n),
			    econst::hcalchan,
			    0,
			    econst::hcalchan,
			    160,
			    -20,
			    20);

    for( Int_t b=0; b<econst::hcalchan; b++ ){
      hadcTdiff_doubles[b][n] = new TH1D(Form("hadcTdiff_doubles_b%d_i%d",b,n),
					 Form("ADCt, Primary Block:%d, Nidx:%d",b,n),
					 160,
					 -20,
					 20);

      htdcdiff_doubles[b][n] = new TH1D(Form("htdcdiff_doubles_b%d_i%d",b,n),
					 Form("TDC, Primary Block:%d, Nidx:%d",b,n),
					 160,
					 -20,
					 20);
    }
  }


  // setup reporting indices
  Int_t curmag = -1;
  std::string roottar = "";
  Int_t good_events = 0;

  Int_t test_idx = 0;

  for ( Int_t t=0; t<2; t++ ){ //loop over targets
    // t==0, lh2; t==1, ld2
    std::string rootfile_fdir;

    //for now, select only LH2 for analysis
    if( t==1 )
      continue;

    if( t==0 ){
      roottar = "LH2";
      nruns = nruns_h;
    }else{
      roottar = "LD2";
      nruns = nruns_d;
    }

    rootfile_fdir = rootfile_dir + Form("/%s/rootfiles",roottar.c_str());

    for (Int_t irun=0; irun<nruns; irun++) { //loop over runs in each target category
      // accessing run info
      test_idx++;
     
      if( test_idx>2 )
	continue;
      

      Int_t runnum;
      Int_t mag;
      Double_t ebeam;
      Int_t conf;
      std::string targ;
      std::string rfname;
      if( t==0 ){
	runnum = crunh[irun].runnum;
	mag = crunh[irun].sbsmag / 21; //convert to percent
	ebeam = crunh[irun].ebeam; //get beam energy per run
	conf = crunh[irun].sbsconf; //get config of current run
	targ = crunh[irun].target;
	rfname = rootfile_fdir + Form("/*%d*",crunh[irun].runnum);

	cout << "rootfile name: " << rfname << endl;
      }
      if( t==1 ){
	runnum = crund[irun].runnum;
	mag = crund[irun].sbsmag / 21;
	ebeam = crund[irun].ebeam;
	conf = crund[irun].sbsconf; //get config of current run
	targ = crund[irun].target;
	rfname = rootfile_fdir + Form("/*%d*",crund[irun].runnum);

      }      

      if( conf!=kine && kine!=-1) continue; //do not proceed if kinematic is constrained
      //if( mag!=fitparams[gIdx].TOFfitp_mset ) continue; //continue if TOF corrections not available

      std::cout << "Analyzing run " << runnum << ".." << std::endl;

      //set up configuration and tune objects to load analysis parameters
      SBSconfig config(kine,mag);

      //Obtain configuration pars from config file
      Double_t hcaltheta = config.GetHCALtheta_rad();
      Double_t hcaldist = config.GetHCALdist();
      Double_t sbsdist = config.GetSBSdist();
      Double_t bbthr = config.GetBBtheta_rad(); //in radians

      SBStune tune(kine,mag);
    
      //if( targ.compare(curtar)!=0 || mag!=curmag ){
      if( mag!=curmag ){
	std::cout << "Settings change.." << endl;
	cout << config;
	cout << tune;
	//curtar = targ;
	curmag = mag;
      }

      //Obtain cuts from tune class
      std:string gcut   = tune.Getglobcut_hexp();
      Double_t W2mean   = tune.GetW2mean();
      Double_t W2sig    = tune.GetW2sig();
      Double_t dx0_n    = tune.Getdx0_n();
      Double_t dx0_p    = tune.Getdx0_p();
      Double_t dy0      = tune.Getdy0();
      Double_t dxsig_n  = tune.Getdxsig_n();
      Double_t dxsig_p  = tune.Getdxsig_p();
      Double_t dysig    = tune.Getdysig();
      Double_t atime0   = tune.Getatime0();
      Double_t atimesig = tune.Getatimesig();

      Double_t dymin = dy0 - 2*dysig;
      Double_t dymax = dy0 + 2*dysig;
      Double_t W2min = W2mean - 2*W2sig;
      Double_t W2max = W2mean + 2*W2sig;

      //std::string rfname = rootfile_dir + Form("/*%d*",crunh[irun].runnum);
      C = new TChain("T");
      C->Add(rfname.c_str());

      std::cout << "Added " << rfname << " to the chain with " << C->GetEntries() << " entries.." << std::endl;

      // setting up ROOT tree branch addresses
      C->SetBranchStatus("*",0);    

      // HCal general
      Double_t hcalid, hcale, hcalx, hcaly, hcalr, hcalc, hcaltdc, hcalatime;
      std::vector<std::string> hcalvar = {"idblk","e","x","y","rowblk","colblk","tdctimeblk","atimeblk"};
      std::vector<void*> hcalvarlink = {&hcalid,&hcale,&hcalx,&hcaly,&hcalr,&hcalc,&hcaltdc,&hcalatime};
      rvars::setbranch(C, "sbs.hcal", hcalvar, hcalvarlink);

      // HCal cluster branches
      Double_t hcalcid[econst::maxclus], hcalce[econst::maxclus], hcalcx[econst::maxclus], hcalcy[econst::maxclus], hcalctdctime[econst::maxclus], hcalcatime[econst::maxclus];
      Int_t Nhcalcid;
      std::vector<std::string> hcalcvar = {"id","e","x","y","tdctime","atime","id"};
      std::vector<void*> hcalcvarlink = {&hcalcid,&hcalce,&hcalcx,&hcalcy,&hcalctdctime,&hcalcatime,&Nhcalcid};
      rvars::setbranch(C, "sbs.hcal.clus", hcalcvar, hcalcvarlink, 6);

      // HCal primary clust block branches
      Double_t hcalcbid[econst::maxclus], hcalcbe[econst::maxclus], hcalcbx[econst::maxclus], hcalcby[econst::maxclus], hcalcbrow[econst::maxclus], hcalcbcol[econst::maxclus], hcalcbtdctime[econst::maxclus], hcalcbatime[econst::maxclus];
      Int_t Nhcalcbid;
      std::vector<std::string> hcalcbvar = {"id","e","row","col","x","y","tdctime","atime","id"};
      std::vector<void*> hcalcbvarlink = {&hcalcbid,&hcalcbe,&hcalcbrow,&hcalcbcol,&hcalcbx,&hcalcby,&hcalcbtdctime,&hcalcbatime,&Nhcalcbid};
      rvars::setbranch(C, "sbs.hcal.clus_blk", hcalcbvar, hcalcbvarlink, 8);

      // bbcal clus var
      Double_t eSH, xSH, ySH, rblkSH, cblkSH, idblkSH, atimeSH, ePS, rblkPS, cblkPS, idblkPS, atimePS;
      std::vector<std::string> bbcalclvar = {"sh.e","sh.x","sh.y","sh.rowblk","sh.colblk","sh.idblk","sh.atimeblk","ps.e","ps.rowblk","ps.colblk","ps.idblk","ps.atimeblk"};
      std::vector<void*> bbcalclvarlink = {&eSH,&xSH,&ySH,&rblkSH,&cblkSH,&idblkSH,&atimeSH,&ePS,&rblkPS,&cblkPS,&idblkPS,&atimePS};
      rvars::setbranch(C, "bb", bbcalclvar, bbcalclvarlink);

      // hodoscope cluster mean time
      Int_t Nhodotmean; 
      Double_t hodotmean[econst::maxclus];
      std::vector<std::string> hodovar = {"clus.tmean","clus.tmean"};
      std::vector<void*> hodovarlink = {&hodotmean,&Nhodotmean};
      rvars::setbranch(C, "bb.hodotdc", hodovar, hodovarlink, 1);  
    
      // track branches
      Double_t ntrack, p[econst::maxtrack],px[econst::maxtrack],py[econst::maxtrack],pz[econst::maxtrack],xtr[econst::maxtrack],ytr[econst::maxtrack],thtr[econst::maxtrack],phtr[econst::maxtrack];
      Double_t vx[econst::maxtrack],vy[econst::maxtrack],vz[econst::maxtrack];
      Double_t xtgt[econst::maxtrack],ytgt[econst::maxtrack],thtgt[econst::maxtrack],phtgt[econst::maxtrack];
      std::vector<std::string> trvar = {"n","p","px","py","pz","x","y","th","ph","vx","vy","vz","tg_x","tg_y","tg_th","tg_ph"};
      std::vector<void*> trvarlink = {&ntrack,&p,&px,&py,&pz,&xtr,&ytr,&thtr,&phtr,&vx,&vy,&vz,&xtgt,&ytgt,&thtgt,&phtgt};
      rvars::setbranch(C,"bb.tr",trvar,trvarlink);
    
      // tdctrig branches
      Int_t Ntdctrigid;
      Double_t tdctrig[econst::maxtrack], tdctrigid[econst::maxtrack];
      std::vector<std::string> tdcvar = {"tdcelemID","tdcelemID","tdc"};
      std::vector<void*> tdcvarlink = {&tdctrigid,&Ntdctrigid,&tdctrig};
      rvars::setbranch(C,"bb.tdctrig",tdcvar,tdcvarlink,1);
  
      // ekine branches
      Double_t ekineQ2, ekineW2, ekineeps, ekinenu, ekineqx, ekineqy, ekineqz;
      std::vector<std::string> ekinevar = {"Q2","W2","epsilon","nu","q_x","q_y","q_z"};
      std::vector<void*> ekinevarlink = {&ekineQ2,&ekineW2,&ekineeps,&ekinenu,&ekineqx,&ekineqy,&ekineqz};
      rvars::setbranch(C, "e.kine", ekinevar, ekinevarlink);
    
      // fEvtHdr branches
      UInt_t rnum, gevnum, trigbits;
      std::vector<std::string> evhdrvar = {"fRun","fEvtNum","fTrigBits"};
      std::vector<void*> evhdrlink = {&rnum,&gevnum,&trigbits};
      rvars::setbranch(C,"fEvtHdr",evhdrvar,evhdrlink);
    
      // globalcut branches
      C->SetBranchStatus("bb.tr.n", 1);
      C->SetBranchStatus("bb.ps.e", 1);
      C->SetBranchStatus("bb.sh.e", 1);

      TCut GCut = gcut.c_str();
      TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", GCut, C );

      // get experimental quantities by run
      std::cout << "Uncorrected average beam energy on " << targ << " for run: " << ebeam << std::endl;
      //set up hcal coordinate system with hcal angle wrt exit beamline
      vector<TVector3> hcalaxes; vars::sethcalaxes( hcaltheta, hcalaxes );
      TVector3 hcalorigin = hcaldist*hcalaxes[2] + econst::hcalvoff*hcalaxes[0];
      Double_t BdL = econst::sbsmaxfield * econst::sbsdipolegap * (mag/100); //scale crudely by magnetic field
      Double_t Eloss_outgoing = econst::celldiameter/2.0/sin(bbthr) * econst::lh2tarrho * econst::lh2dEdx;

      // set nucleon defaults by target
      std::string nucleon;
      if( targ.compare("LH2")==0 )
	nucleon = "p"; 
      else if( targ.compare("LD2")==0 )      
	nucleon = "np";
      else{
	nucleon = "np";
	std::cout << "Error: incorrect target information loaded from good run list. Check and verify." << std::endl;
      }

      // event indices
      long nevent = 0, nevents = C->GetEntries(); 
      Int_t treenum = 0, currenttreenum = 0;
      Int_t pltcntr = 0;

      while (C->GetEntry(nevent++)) {

	cout << "Processing: " << nevent << "/" << nevents << " - good events: " << good_events << " \r";
	cout.flush();

	///////
	//Single-loop globalcut method. Save pass/fail for output tree.
	currenttreenum = C->GetTreeNumber();
	if( nevent == 1 || currenttreenum != treenum ){
	  treenum = currenttreenum; 
	  GlobalCut->UpdateFormulaLeaves();
	  cout << "Updating formula leaves and switching segment at event: " << nevent << endl;
	}
	bool failedglobal = GlobalCut->EvalInstance(0) == 0;
      
	///////
	//HCal coin cut. Save pass/fail for output tree.
	bool failedcoin = abs( hcalcbatime[0] - atime0 ) > 3*atimesig;


	///////
	//Physics calculations
	//correct beam energy with vertex information
	Double_t ebeam_c = vars::ebeam_c( ebeam, vz[0], targ );
	TVector3 vertex( 0., 0., vz[0] );

	//reconstructed momentum, corrected for mean energy loss exiting the target (later we'll also want to include Al shielding on scattering chamber)
	Double_t precon = p[0] + Eloss_outgoing;

	//set up four-momenta with some empty for various calculation methods
	TLorentzVector pbeam( 0., 0., ebeam_c, ebeam_c ); //beam momentum
	//TLorentzVector pe( px[0], py[0], pz[0], p[0] ); //e' plvect
	TLorentzVector pe( precon*px[0]/p[0], precon*py[0]/p[0], precon*pz[0]/p[0], precon ); //e' recon plvect
	TLorentzVector ptarg; vars::setPN(nucleon,ptarg); //target momentum
	TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
	TVector3 qv = q.Vect();
	TLorentzVector pN; //N' momentum
	TVector3 pNhat; //Unit N' 3-vector
      
	//simple calculations for e' and N'
	Double_t etheta = vars::etheta(pe); 
	Double_t ephi = vars::ephi(pe);
	Double_t pcent = vars::pcentral(ebeam,etheta,nucleon); //e' p reconstructed by angles
	Double_t phNexp = ephi + physconst::pi;
	Double_t Q2, W2, nu, thNexp, pNexp, pNcalc, pp, pn;
	Double_t ebeam_o = vars::ebeam_o( ebeam_c, etheta, targ ); //Second energy correction accounting for energy loss leaving target

	/* Can reconstruct e' momentum for downstream calculations differently:
	   v1 - Use four-momentum member functions
	   v2 - Use all available ekine (tree) vars and calculate vectors (should be the same as v1)
	   v3 - Use reconstructed angles as independent qty (usually preferable given GEM precision at most kinematics)
	   v4 - Use reconstructed momentum as independent qty */
      
	if( epm==1 ){
	  //v1
	  pN = q + ptarg;
	  pNhat = pN.Vect().Unit();
	  Q2 = -q.M2();
	  W2 = pN.M2();
	  nu = q.E();
	  pNexp = vars::pN_expect( nu, nucleon );
	  thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	}else if( epm==2 ){
	  //v2
	  Q2 = ekineQ2;
	  W2 = ekineW2;
	  nu = ekinenu;
	  pNexp = vars::pN_expect( nu, nucleon );
	  thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	  pNhat = vars::pNhat_track( thNexp, phNexp );
	  pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	}else if( epm==3 ){
	  //v3
	  nu = pbeam.E() - pcent;
	  pNexp = vars::pN_expect( nu, nucleon );
	  thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	  pNhat = vars::pNhat_track( thNexp, phNexp );
	  pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	  Q2 = vars::Q2( pbeam.E(), pe.E(), etheta );
	  W2 = vars::W2( pbeam.E(), pe.E(), Q2, nucleon );
	}else if( epm==4 ){
	  //v4
	  nu = pbeam.E() - pe.E();
	  pNexp = vars::pN_expect( nu, nucleon );
	  thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	  pNhat = vars::pNhat_track( thNexp, phNexp );
	  pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	  Q2 = vars::Q2( pbeam.E(), pe.E(), etheta );
	  W2 = vars::W2( pbeam.E(), pe.E(), Q2, nucleon );
	}else{
	  Q2 = ekineQ2;
	  W2 = ekineW2;
	  nu = ekinenu;
	  pNexp = vars::pN_expect( nu, nucleon );
	  thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
	  pNhat = vars::pNhat_track( thNexp, phNexp );
	  pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
	  std::cout << "Warning: epm version incorrect. Defaulting to version 2." << endl;
	}
	pp = sqrt(pow(nu,2)+2.*physconst::Mp*nu); //momentum of proton
	pn = sqrt(pow(nu,2)+2.*physconst::Mn*nu); //momentum of neutron

	//Calculate h-arm quantities
	vector<Double_t> xyhcalexp; vars::getxyhcalexpect( vertex, pNhat, hcalorigin, hcalaxes, xyhcalexp );
	TVector3 hcalpos = hcalorigin + hcalcbx[0]*hcalaxes[0] + hcalcby[0]*hcalaxes[1]; //from primary blk
	Double_t dx = hcalcbx[0] - xyhcalexp[0];
	Double_t dy = hcalcby[0] - xyhcalexp[1];
	TVector3 neutdir = ( hcalpos - vertex ).Unit();
	Double_t protdeflect = tan( 0.3 * BdL / q.Mag() ) * (hcaldist - (sbsdist + econst::sbsdipolegap/2.0) );
	TVector3 protdir = ( hcalpos + protdeflect * hcalaxes[0] - vertex ).Unit();
	Double_t thetapq_p = acos( protdir.Dot( qv.Unit() ) );
	Double_t thetapq_n = acos( neutdir.Dot( qv.Unit() ) );

	// //Get PID
	// Int_t pid;
	// bool is_p = abs(dx - dx0_p)<3*dxsig_p && abs(dy - dy0)<3*dysig;
	// bool is_n = abs(dx - dx0_n)<3*dxsig_n && abs(dy - dy0)<3*dysig;
	// if( (is_p && !is_n) || targ.compare("LH2")==0 ){
	//   pid=1;
	//   pNcalc = vars::pN_expect( nu, "p" );
	// }else if( is_n && !is_p ){
	//   pid=2;
	//   pNcalc = vars::pN_expect( nu, "n" );
	// }else if( is_p && is_n ){
	//   pid=0;
	//   pNcalc = vars::pN_expect( nu, "n" );
	// }else{
	//   pid=-1;
	//   pNcalc = pNexp;
	// }

	///////
	//HCal active area cut (acceptance matching). Save pass/fail for output tree.
	// bool failedaccmatch = 
	//   xyhcalexp[1] > (econst::hcalposYf_mc-econst::hcalblk_div_h) ||
	//   xyhcalexp[1] < (econst::hcalposYi_mc+econst::hcalblk_div_h) ||
	//   xyhcalexp[0] > (econst::hcalposXf_mc-econst::hcalblk_div_v) ||
	//   xyhcalexp[0] < (econst::hcalposXi_mc+econst::hcalblk_div_v);

	/////////////////////
	//W2 elastic cut bool
	bool failedW2 = W2<W2min || W2>W2max;

	/////////////////////
	//dy cut
	bool faileddy = dy<dymin || dy>dymax;

	/////////////////////
	//dx cut
	bool faileddx = dx<-1.2 || dx>1.0; //empirical

	/////////////////////
	//pp cut
	bool failedpp = pp<2.4 || pp>2.6;

	/////////////////////
	//hcal primary bloc energy cut
	Double_t pblke = hcalcbe[0];
	bool failede = pblke<0.2 || pblke>0.4;

	///////////////////////////////////////////
	///////////////////////////////////////////
	//Make elastic cut here to constrain timing
	//if( failedglobal || failedW2 || faileddy || faileddx || failedpp || failede)
	//continue;



	good_events++;
	
	//Primary cluster primary block
	Int_t pblkrow = (int)hcalcbid[0]/econst::hcalcol;
	Int_t pblkcol = (int)hcalcbid[0]%econst::hcalcol;
	Int_t pblkid = (int)hcalcbid[0];
	Int_t nblk = Nhcalcbid;
	Double_t padct = hcalcbatime[0];
	Double_t ptdc = hcalcbtdctime[0];
	//Double_t ptdc_tw = ptdc + 15.6373*pblke; //via pol1 fit to tdc vs E min 0.05
	Double_t ptdc_tw = ptdc + 32.4746*pblke; //via pol1 fit to tdc vs E min 0.03

	//Skip events where the primary block is on the edge of the acceptance
	if( pblkrow == 0 || pblkrow == 23 || pblkcol == 0 || pblkcol == 11 )
	  continue;

	hpblke->Fill(pblke);
	hpp->Fill(pp);
	hdx_HCAL->Fill(dx);
	hdy_HCAL->Fill(dy);
	hW2->Fill(W2);
	hadct_E->Fill(pblke,padct);
	htdc_E->Fill(pblke,ptdc+100);
	htdc_E_corr->Fill(pblke,ptdc_tw+100);
	hpse->Fill(ePS);


	if( failedglobal || failedW2 )
	  continue;


	//cout << "event:" << nevent << " pblkid:" << pblkid << endl;

	//Loop over remaining blocks skipping the first element to fill diagnostic/histograms
	for( int b=1; b<nblk; b++){
	  Int_t blkid = (int)hcalcbid[b];
	  //cout << " nblk:" << nblk << " element:" << b << " id:" << blkid;;

	  Double_t adct = hcalcbatime[b];
	  Double_t tdc = hcalcbtdctime[b];
	  Int_t blkrow = blkid/econst::hcalcol;
	  Int_t blkcol = blkid%econst::hcalcol;
	  Double_t blke = hcalcbe[b];
	  //Double_t tdc_tw = tdc + 15.6373*blke; //fit min at 0.05
	  Double_t tdc_tw = tdc + 32.4746*blke; //fit min at 0.03

	  bool blke_cut = hcalcbe[b]>0.03 && hcalcbe[b]<0.15;

	  Double_t diff_adct = padct - adct;
	  Double_t diff_tdc = ptdc - tdc;
	  Double_t diff_tdc_tw = ptdc_tw - tdc_tw;

	  bool row_neighbor = false;
	  if( blkrow == pblkrow-1 || 
	      blkrow == pblkrow+1 ||
	      blkrow == pblkrow )
	    row_neighbor = true;

	  bool col_neighbor = false;
	  if( blkcol == pblkcol-1 || 
	      blkcol == pblkcol+1 ||
	      blkcol == pblkcol )
	    col_neighbor = true;

	  //only check neighbors and throw out cluster members that aren't
	  if( !row_neighbor || !col_neighbor ){
	    //cout << " NOT a neighbor" << endl;
	    continue;
	  }
	  Int_t idx_neighbor;
	  if( blkrow == pblkrow-1 )
	    idx_neighbor = -(pblkid-blkid)+13;
	  if( blkrow == pblkrow )
	    idx_neighbor = -(pblkid-blkid)+4;
	  if( blkrow == pblkrow+1 )
	    idx_neighbor = -(pblkid-blkid)-5;

	  //cout << " idx:" << idx_neighbor << endl;

	  // if( row_neighbor && col_neighbor ){
	  //   tdiff_adct[pblkid][idx_neighbor]->Fill( diff_adct );
	  //   tdiff_tdc[pblkid][idx_neighbor]->Fill( diff_tdc );
	  // }

	  //fill array and histos with available events and data
	  Nevents[pblkid][idx_neighbor]++;
	  hElastic_doubles[idx_neighbor]->Fill(pblkcol,pblkrow);
	  hadcTdiff_doubles[pblkid][idx_neighbor]->Fill(diff_adct);
	  hadiff_id[idx_neighbor]->Fill(pblkid,diff_adct);
	  

	  if( blke_cut ){
	    Neventstdc[pblkid][idx_neighbor]++;
	    htdiff_id[idx_neighbor]->Fill(pblkid,diff_tdc);
	    htdcdiff_doubles[pblkid][idx_neighbor]->Fill(diff_tdc_tw);
	  }
	  //cout << hcale << endl;

	  

	  if( b==1 && hcale>0.1 && abs(tdc)<100 && abs(ptdc)<100 ){
	    //cout << diff_tdc << " = " << ptdc << " - " << tdc << endl;
	    hsblke->Fill(blke);

	    Nevents_gtdc[pblkid][idx_neighbor]++;
	    if( pblkid==126 && idx_neighbor==1 )
	      hgtdc->Fill(diff_tdc);
	  }
	}

      }//end loop over event

      // getting ready for the next run
      C->Reset();

    }//endloop over runs

  }//endloop over targets

  //diagnostic, write out 
  Int_t chan = 0;
  Int_t tdc_best_chan = -1;
  Int_t tdc_best_idx = -1;
  Int_t tdc_best_N = 0;

  TCanvas *c1[9];
  TCanvas *c2[9];

  gStyle->SetEndErrorSize(0);

  for( Int_t n=0; n<9; n++ ){
    //Skip the index which corresponds to the difference of pblkid with itself
    if(n==4)
      continue;

    c1[n] = new TCanvas(Form("c1_%d",n),Form("ADCt Index %d",n),1200,1200);
    c1[n]->cd();

    overlayWithGaussianFits(hadiff_id[n],c1[n]);

    c1[n]->Write();

    c2[n] = new TCanvas(Form("c2_%d",n),Form("TDC Index %d",n),1200,1200);
    c2[n]->cd();

    overlayWithGaussianFits(htdiff_id[n],c2[n]);

    c2[n]->Write();

    cout << "Idx set " << n << " fitted and complete" << endl;
  }


  for( Int_t r=0; r<econst::hcalrow; r++ ){
    for( Int_t c=0; c<econst::hcalcol; c++ ){
      for( Int_t n=0; n<9; n++ ){
	if( Nevents_gtdc[chan][n] > tdc_best_N ){
	  tdc_best_N = Nevents_gtdc[chan][n];
	  tdc_best_chan = chan;
	  tdc_best_idx = n;
	}
	//cout << Nevents[chan][n] << ",";

	Int_t adctN = hadcTdiff_doubles[chan][n]->GetEntries();
	Int_t tdcN = htdcdiff_doubles[chan][n]->GetEntries();

	hadcTdiff_doubles[chan][n]->SetName(Form("hadctDiff_b%d_i%d_e%d",chan,n,adctN));
	htdcdiff_doubles[chan][n]->SetName(Form("htdcDiff_b%d_i%d_e%d",chan,n,tdcN));

	  cout << "ADCt " << chan << ":" << n << " " << hadcTdiff_doubles[chan][n]->GetEntries() << ":" << hadcTdiff_doubles[chan][n]->GetMaximumBin() << endl;

	  
	  // if( chan==266 )
	  //   continue;

	  //declare a 1D fit function to fit no cut adc time histogram
	  TF1 *f1;
	  TF1 *f2;

	  //ADCt pair fits
	  if( Nevents[chan][n] > 500 && hadcTdiff_doubles[chan][n]->GetMaximumBin() > 20 && hadcTdiff_doubles[chan][n]->GetMaximumBin() < 140 ){
	    //declare some dynamic fit variables
	    Int_t atimeBinMax = hadcTdiff_doubles[chan][n]->GetMaximumBin();
	    Double_t atimeBinCenter = hadcTdiff_doubles[chan][n]->GetBinCenter( atimeBinMax );
	    Double_t atime_fitLowerLim = atimeBinCenter - atime_approx_FWHM;
	    Double_t atime_fitUpperLim = atimeBinCenter + atime_approx_FWHM;
	    
	    //Do the adct fit
	    if( hadcTdiff_doubles[chan][n]->GetEntries()>100 && 
		atimeBinCenter>atime_fitLowerLim &&
		atimeBinCenter<atime_fitUpperLim ){
	      hadcTdiff_doubles[chan][n]->Fit( "gaus", "Q", "", atime_fitLowerLim, atime_fitUpperLim ); //fit the no cut adc time with a root predifined gaussian with the dynamic limits defined above
	      f1 = hadcTdiff_doubles[chan][n]->GetFunction( "gaus" ); //set the fit function equal to this fit
	      
	      //extract fit parameters. gaussian has three parameters, f(x) = p0 * exp(-0.5 * ((x - p1) / p2)^2)
	      Double_t atime_mean = f1->GetParameter(1); //where p1 is the mean
	      Double_t atime_sig = f1->GetParameter(2); //where p2 is the std dev
	      hadctsig->Fill(atime_sig);
	      
	      hadcTdiff_doubles[chan][n]->SetTitle(Form("ADC time elastic cut channel %d index %d. Mean:%f Sigma:%f",chan,n,atime_mean,atime_sig)); //Write the relevant fit parameters to the title of the histogram.
	    }else{
	      cout << "ADCt ERROR on chan " << chan << " index " << n << " entries:bincenter " << hadcTdiff_doubles[chan][n]->GetEntries() << ":" << atimeBinCenter << endl;
	    }
	  }

	  cout << "TDC " << chan << ":" << n << " " << htdcdiff_doubles[chan][n]->GetEntries() << ":" << htdcdiff_doubles[chan][n]->GetMaximumBin() << endl;

	  if( Neventstdc[chan][n] > 250 && htdcdiff_doubles[chan][n]->GetMaximumBin() > 20 && htdcdiff_doubles[chan][n]->GetMaximumBin() < 140){
	    Int_t tdcBinMax = htdcdiff_doubles[chan][n]->GetMaximumBin();
	    Double_t tdcBinCenter = htdcdiff_doubles[chan][n]->GetBinCenter( tdcBinMax );
	    Double_t tdc_fitLowerLim = tdcBinCenter - tdc_approx_FWHM;
	    Double_t tdc_fitUpperLim = tdcBinCenter + tdc_approx_FWHM;
	    
	    if( htdcdiff_doubles[chan][n]->GetEntries()>100 && 
		tdcBinCenter>tdc_fitLowerLim &&
		tdcBinCenter<tdc_fitUpperLim ){
	      //Do the tdc fit
	      htdcdiff_doubles[chan][n]->Fit( "gaus", "Q", "", tdc_fitLowerLim, tdc_fitUpperLim );
	      f2 = htdcdiff_doubles[chan][n]->GetFunction( "gaus" );
	      
	      Double_t tdc_mean = f2->GetParameter(1); //where p1 is the mean
	      Double_t tdc_sig = f2->GetParameter(2); //where p2 is the std dev
	      htdcsig->Fill(tdc_sig);
	      string lowsig = "";
	      if( tdc_sig<1.1 )
		lowsig=std::to_string(tdc_sig);
	      
	      htdcdiff_doubles[chan][n]->SetTitle(Form("TDC time elastic cut channel %d index %d. Mean:%f Sigma:%f",chan,n,tdc_mean,tdc_sig)); //Write the relevant fit parameters to the title of the histogram.
	      htdcdiff_doubles[chan][n]->SetName(Form("htdcDiff_b%d_i%d_e%0.1f_%s",chan,n,htdcdiff_doubles[chan][n]->GetEntries(),lowsig.c_str()));
	      
	    }else{
	      cout << "TDC ERROR on chan " << chan << " index " << n << " entries:bincenter " << htdcdiff_doubles[chan][n]->GetEntries() << ":" << tdcBinCenter << endl;
	    }
	    
	    //hadctsig->Fill(atime_sig);
	    //htdcsig->Fill(tdc_sig);
	  
	  }
      }
      chan++;
    }
    cout << endl;
    
  }

  cout << endl << endl << "TDC best channel with high E cuts: " << tdc_best_chan << " idx: " << tdc_best_idx << " with " << tdc_best_N << " events." << endl << endl;

  fout->Write();

  st->Stop();

  // Send time efficiency report to console
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
