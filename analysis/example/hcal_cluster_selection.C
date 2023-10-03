//SSeeds 9.1.23 Standalone script to sort through hcal cluster and evaluate them based on separate metrics

#include <vector>
#include <iostream>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TVector.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"

//global variables
const int maxTracks = 100; 
const int maxClus = 100;
const double coinSigFac = 4;
const double dSigFac = 3;
const double PI = TMath::Pi();
const double M_e = 0.0051;
const double M_p = 0.938272;
const double M_n = 0.939565;
const double hcalheight = -0.2897; //relevant to correct for bad HCal height in DB on pass0/1
const double binfac = 400;
//const double hbinfac = 14.286; //results in 7cm hcal resolution
//const double hbinfac = 28.572; //results in 3.5cm hcal resolution, 7cm inherent averaged
const double hbinfac = 100.; 
const double W2fitmax = 1.4;
// LH2
const double lh2tarrho = 0.0723;     //g/cc, target density
const double lh2cthick = 0.02;       //cm, target cell thickness
const double lh2uwallthick = 0.0145; //cm, upstream wall thickness
const double lh2dwallthick = 0.015;  //cm, downstream wall thickness
const double lh2dEdx = 0.00574; 
// LD2
const double ld2tarrho = 0.169;      //g/cc, target density
const double ld2dEdx = 0.0114;      
const double ld2uwallthick = 0.0145; //cm, assume same as hydrogen for now
const double ld2dwallthick = 0.015;  //cm, assume same as hydrogen for now
// Aluminum
const double rho_al = 2.7; //g/cc
const double aldEdx = 0.0021;
//all targ
const double l_tgt = 0.15;
//hcal dims
const double hcalposXi = -2.655;   //m, distance from beam center to top of HCal w/75cm offset
const double hcalposXf = 1.155;   //m, distance from beam center to bottom of HCal w/75cm offset
const double hcalposYi = -0.92964; //m, distance from beam center to opposite-beam side of HCal
const double hcalposYf = 0.92964;  //m, distance from beam center to beam side of HCal
const double hcal_vrange = 3.81;
const double hcal_hrange = 1.85928;
//cluster pass weights - may not need these after implementation of PDF-based assignments
const double eweight = 1.00;
const double tweight = 0.57;
const double dweight = 0.81;
//score assigments based on dx probability density
const double dx_fit7[3] = {800,-0.693,0.062}; //sbs7, 70% field
//const double dx_fit4_0[3] = {40600,0.019,0.059}; //sbs4, zero field
const double dx_fit4[3] = {9300,-1.092,0.093}; //sbs4, 50% field
const double dx_fit14[3] = {800,-0.693,0.062}; //sbs14
const double dx_fit11[3] = {800,-0.693,0.062}; //sbs11
const double dx_fit8[3] = {800,-0.693,0.062}; //sbs8
const double dx_fit9[3] = {800,-0.693,0.062}; //sbs9

//using atime here where coin is non-gaussian. will move to coin everywhere after pass2
const double atime_fit4[3] = {2446,51.69,3.96}; //hcal
const double atime_fit7[3] = {124.8,65.52,1.70}; //coin
const double atime_fit11[3] = {47.7,60.28,1.89}; //coin
const double atime_fit14[3] = {694,59.30,1.52}; //coin
const double atime_fit8[3] = {21682,47.33,1.61}; //coin
const double atime_fit9[3] = {5620,52.00,1.57}; //coin

/////DIAGNOSTIC
vector<int> sbs4_0p_runs = {11573,11587,11588}; //zero field for diagnostics
bool diagnostic = false; //enable to run the above three runs at zero field for diagnostics
/////

vector<int> sbs4_50p_runs = {11589,11590,11592};
vector<int> sbs7_85p_runs = {11989,11991,11992,11993,11994,12000};
vector<int> sbs11_100p_runs = {12355,12363,12367,12368,12369,12400};
vector<int> sbs14_70p_runs = {13242,13243,13244,13312,13320,13321};
vector<int> sbs8_70p_runs = {13482,13483,13484,13485,13486,13487};
vector<int> sbs9_70p_runs = {13656,13657,13663,13676,13683,13696};

//assign score to potential cluster based on probability density of gaussian fit to same data
double assignScore( double val, string type, const std::vector<double> &dxfit ) {

  if( dxfit.size()!=3 ){
    cout << "ERROR: size of dxfit vector not equal to expected number of gaussian parameters (3)." << endl;
    return 0.;
  }

  double weight;
  if( type.compare("e")==0 )
    weight = eweight;
  if( type.compare("t")==0 )
    weight = tweight;
  if( type.compare("d")==0 )
    weight = dweight;

  // Fit the histogram with a Gaussian function
  TF1 *gauss = new TF1("gauss", "gaus", -5, 5);
  gauss->SetParameter(0,dxfit[0]);
  gauss->SetParameter(1,dxfit[1]);
  gauss->SetParameter(2,dxfit[2]);
    
  double x_value = val;
  double density = gauss->Eval(x_value);

  // Compute the score
  double score = weight * density / gauss->GetParameter(0);  // Normalize by the peak of the Gaussian

  if(score==0){
    //cout << "WARNING: score is zero. val=" << val << " type=" << type << " density= " << density << endl;
    score=1e-38;
  }

  return score;
}

//assign score to potential cluster based on probability density of gaussian fit to atime and maximum energy, exclude position info. Here, gfit is a vector with coin atime fit parameters, val_2 is the cluster coin atime, val_1 is the cluster energy, and max1 is the max cluster energy for the event
double assignScoreAlt( double val_1, double val_2, double max1, const std::vector<double> &gfit ) {

  if( gfit.size()!=3 ){
    cout << "ERROR: size of dxfit vector not equal to expected number of gaussian parameters (3)." << endl;
    return 0.;
  }

  // Build a TF1 with fit provided fit parameters
  TF1 *gauss = new TF1("gauss", "gaus", 0, 100);
  gauss->SetParameter(0,gfit[0]);
  gauss->SetParameter(1,gfit[1]);
  gauss->SetParameter(2,gfit[2]);
    
  double x_value = val_2;
  double density = gauss->Eval(x_value);

  // Compute the score based on val 2
  double score = density / gauss->GetParameter(0);  // Normalize by the peak of the Gaussian

  // Include a product component based on val 1 (linear weight)
  score *= val_1 / max1;

  if(score==0){
    //cout << "WARNING: score is zero. val_1=" << val_1 << ", val_2= " << val_2 << " density= " << density << " 1 comp= " << val_1 / max1 << endl;
    score=1e-38; //write very small number to score to exclude it without breaking the sorting later
  }

  return score;
}


//cluster non-zero sorting. if two single valued vectors equal to zero are returned, no non-zero values were found
std::pair<std::vector<double>, std::vector<size_t>> sortNonZero(const std::vector<double> &arr) {
    // Pair of value and its index in the original array
    std::vector<std::pair<double, size_t>> nonZeroPairs;

    // Populate the pairs
    bool empty = true;
    for (size_t i = 0; i < arr.size(); ++i) {
        if (arr[i] != 0.0) {
            nonZeroPairs.push_back({arr[i], i});
	    empty = false;
        }
    }

    //prevent the seg fault!
    if( empty )
      nonZeroPairs.push_back({0.,0.});

    // Sort the pairs based on value in descending order
    std::sort(nonZeroPairs.begin(), nonZeroPairs.end(), [](const auto& a, const auto& b) {
        return a.first > b.first;
    });

    std::vector<double> sortedValues;
    std::vector<size_t> sortedIndices;

    for (const auto& pair : nonZeroPairs) {
        sortedValues.push_back(pair.first);
        sortedIndices.push_back(pair.second);
    }

    return {sortedValues, sortedIndices};
}

//implement cut on dxdy
bool IsInsideEllipse(double x, double y, double x_mean, double x_sigma, double y_mean, double y_sigma) {
    double term1 = (x - x_mean) * (x - x_mean) / (x_sigma * x_sigma);
    double term2 = (y - y_mean) * (y - y_mean) / (y_sigma * y_sigma);

    return (term1 + term2) <= 1.0;
}

//NOTE: this analysis should be performed on LH2 to constrain dx position
//main
void hcal_cluster_selection( int kine = 4, int pass = 0, string target = "LH2", int run = -1, bool shorty = true, bool coin = false ){  
  
  //Set up the output file
  string short_word = "";
  if(shorty)
    short_word = "_short";
  string diag_word = "";
  if(diagnostic)
    diag_word = "_diagnostic";

  //string outfile_directory = Form("/volatile/halla/sbs/seeds/hcal_cluster_selection_out_sbs%d.root",kine);
  string outfile_directory = Form( "outfiles/hcs_out_sbs%d%s%s.root", kine, short_word.c_str(), diag_word.c_str() );
  TFile *fout = new TFile( outfile_directory.c_str(), "RECREATE" );

  //Switch between settings based on kinematic
  string default_extension;
  string short_extension;
  double E_e;
  double BB_d;
  double BB_th;
  double HCal_d;
  double HCal_th;
  double W2_mean;
  double W2_sig;
  double dy_mean;
  double dy_sig;
  double dx_mean;
  double dx_sig;
  double coin_mean;
  double coin_sig;
  double minE;
  vector<int> short_extension_runs;

  switch (kine) {
  case 4:
    short_extension = "/w/halla-scshelf2102/sbs/sbs-gmn/pass0/SBS4/LH2/rootfiles/e1209019_fullreplay_";
    default_extension = "/w/halla-scshelf2102/sbs/sbs-gmn/pass0/SBS4/LH2/rootfiles/*.root";
    E_e = 3.7393; //Beam energy in GeV
    BB_d = 1.7988; //Distance to bigbite from target in m
    BB_th = TMath::DegToRad() * 36.0; //Angle wrt beamline for bigbite arm in radians
    HCal_d =11.0; //Distance to HCal from target in m
    HCal_th = TMath::DegToRad() * 31.9; //Angle wrt beamline for HCal in radians
    W2_mean = 1.00; //W^2 location of the mean of the elastic peak
    W2_sig = 0.24; //W^2 width of the elastic peak
    dy_mean = -1.687; //Elastic peak in dy mean
    dy_sig = 0.69; //Elastic peak in dy sigma
    dx_mean = -1.09359;
    dx_sig = 0.0879841;
    coin_mean = 56.;
    coin_sig = 4.;
    minE = 0.020;
    for( int r=0; r<sbs4_50p_runs.size(); ++r )
      short_extension_runs.push_back(sbs4_50p_runs[r]);
    break;
  case 7:
    short_extension = "/w/halla-scshelf2102/sbs/sbs-gmn/pass0/SBS7/LH2/rootfiles/e1209019_fullreplay_";
    default_extension = "/w/halla-scshelf2102/sbs/sbs-gmn/pass0/SBS7/LH2/rootfiles/*.root";
    E_e = 7.9072;
    BB_d = 1.850;
    BB_th = 40.0*TMath::DegToRad();
    HCal_d =14.0;
    HCal_th = 16.1*TMath::DegToRad();
    W2_mean = 0.98;
    W2_sig = 0.35;
    dy_mean = -0.0276561;
    dy_sig = 0.283612;
    dx_mean = -0.691916;
    dx_sig = 0.0636087;
    coin_mean = 51.;
    coin_sig = 6.;
    minE = 0.050;
    for( int r=0; r<sbs7_85p_runs.size(); ++r )
      short_extension_runs.push_back(sbs7_85p_runs[r]);
    break;
  case 11:
    short_extension = "/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS11/LH2/rootfiles/e1209019_fullreplay_";
    default_extension = "/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS11/LH2/rootfiles/*.root";
    E_e = 9.8594;
    BB_d = 1.550;
    BB_th = 42.0*TMath::DegToRad();
    HCal_d =14.5;
    HCal_th = 13.3*TMath::DegToRad();
    W2_mean = 0.98;
    W2_sig = 0.35;
    dy_mean = 0.015;
    dy_sig = 0.083;
    dx_mean = -0.631161;
    dx_sig = 0.0772046;
    coin_mean = 51.;
    coin_sig = 6.;
    minE = 0.050;
    for( int r=0; r<sbs11_100p_runs.size(); ++r )
      short_extension_runs.push_back(sbs11_100p_runs[r]);
    break;
  case 14:
    short_extension = "/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS14/LH2/rootfiles/e1209019_fullreplay_";
    default_extension = "/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS14/LH2/rootfiles/*.root";
    E_e = 5.9649;
    BB_d = 1.850;
    BB_th = 46.5*TMath::DegToRad();
    HCal_d =14.0;
    HCal_th = 17.3*TMath::DegToRad();
    W2_mean = 0.98;
    W2_sig = 0.35;
    dy_mean = 0.001;
    dy_sig = 0.045;
    dx_mean = -0.722565;
    dx_sig = 0.0761898;
    coin_mean = 51.;
    coin_sig = 6.;
    minE = 0.050;
    for( int r=0; r<sbs14_70p_runs.size(); ++r )
      short_extension_runs.push_back(sbs14_70p_runs[r]);
    break;
  case 8:
    short_extension = "/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS8/LH2/rootfiles/e1209019_fullreplay_";
    default_extension = "/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS8/LH2/rootfiles/*.root";
    E_e = 5.965;
    BB_d = 1.97473;
    BB_th = TMath::DegToRad() * 26.5;
    HCal_d =11.0;
    HCal_th = TMath::DegToRad() * 29.4;
    W2_mean = 1.00;
    W2_sig = 0.24;
    dy_mean = -1.687;
    dy_sig = 0.69;
    dx_mean = -0.769036;
    dx_sig = 0.139637;
    coin_mean = 51.;
    coin_sig = 6.;
    minE = 0.050;
    for( int r=0; r<sbs8_70p_runs.size(); ++r )
      short_extension_runs.push_back(sbs8_70p_runs[r]);
    break;
  case 9:
    short_extension = "/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS9/LH2/rootfiles/e1209019_fullreplay_";
    default_extension = "/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS9/LH2/rootfiles/*.root";
    E_e = 4.013;
    BB_d = 1.550;
    BB_th = 49.0*TMath::DegToRad();
    HCal_d =11.0;
    HCal_th = 22.0*TMath::DegToRad();
    W2_mean = 0.98;
    W2_sig = 0.35;
    dy_mean = -0.0276561;
    dy_sig = 0.283612;
    dx_mean = -0.819153;
    dx_sig = 0.0824273;
    coin_mean = 51.;
    coin_sig = 6.;
    minE = 0.050;
    for( int r=0; r<sbs9_70p_runs.size(); ++r )
      short_extension_runs.push_back(sbs9_70p_runs[r]);
    break;
  default:
    kine = 4;
  }
  
  //Select the rootfile or just use the default
  string rootfile_parent_directory = Form( "/w/halla-scshelf2102/sbs/sbs-gmn/pass%d/SBS%d/", pass, kine );
  string rootfile_path;
  if( run == -1 && shorty )
    rootfile_path = short_extension;
  else if( run == -1 && !shorty )
    rootfile_path = default_extension;
  else{
    rootfile_path = rootfile_parent_directory + Form("%s/rootfiles/*%d*",target.c_str(),run);
    cout << "No default run selected. Please enter dx location of proton peak." << endl;
    cin >> dx_mean;
    cout << "Please enter dx sigma of proton peak." << endl;
    cin >> dx_sig;
  }

  //Very generic globalcut. Keep events where no clusters exist in HCal for analysis
  TCut globalcut = "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.075&&bb.sh.nclus>0&&sbs.hcal.nclus>0";

  ////////////
  //HISTOGRAMS

  //Basic
  TH1D *hW2_nocut = new TH1D( "hW2_nocut", "W^{2}, no cut; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hHCALadiff = new TH1D( "hHCALadiff", "HCal ADCt - BBCal SH ADCt, pblkclus; ns", 100, 0., 100. );
  TH1D *hHCALe_W2cut = new TH1D( "hHCALe_W2cut", "HCal Energy, W2 Cut; GeV", 150, 0., 1.5 );

  //Scoring (best cluster)
  // TH1D *hscore_crude = new TH1D( "hscore_crude", "Best Cluster Weighted Score; score", 100, 0.0, 5. );
  // TH1D *hscore_allWeights = new TH1D( "hscore_allWeights", "Best Cluster and dx Location Weighted Score; score", 100, 0.0, 5. );

  //Basic cluster dx
  TH1D *hdx_pclus = new TH1D( "hdx_pclus", "dx, primary cluster;x_{HCAL}-x_{expect} (m)", hbinfac*hcal_vrange, hcalposXi, hcalposXf);
  TH1D *hdx_eclus = new TH1D( "hdx_eclus", "dx, highest e cluster;x_{HCAL}-x_{expect} (m)", hbinfac*hcal_vrange, hcalposXi, hcalposXf);
  TH1D *hdx_dclus = new TH1D( "hdx_dclus", "dx, closest dxdy cluster;x_{HCAL}-x_{expect} (m)", hbinfac*hcal_vrange, hcalposXi, hcalposXf);
  TH1D *hdx_tclus = new TH1D( "hdx_tclus", "dx, lowest time difference cluster;x_{HCAL}-x_{expect} (m)", hbinfac*hcal_vrange, hcalposXi, hcalposXf);
  TH1D *hdy_pclus = new TH1D( "hdy_pclus", "dy, primary cluster;y_{HCAL}-y_{expect} (m)", hbinfac*hcal_hrange, hcalposYi, hcalposYf);
  TH1D *hdy_eclus = new TH1D( "hdy_eclus", "dy, highest e cluster;y_{HCAL}-y_{expect} (m)", hbinfac*hcal_hrange, hcalposYi, hcalposYf);
  TH1D *hdy_dclus = new TH1D( "hdy_dclus", "dy, closest dxdy cluster;y_{HCAL}-y_{expect} (m)", hbinfac*hcal_hrange, hcalposYi, hcalposYf);
  TH1D *hdy_tclus = new TH1D( "hdy_tclus", "dy, lowest time difference cluster;y_{HCAL}-y_{expect} (m)", hbinfac*hcal_hrange, hcalposYi, hcalposYf);
  
  TH1D *hdx_bclus = new TH1D( "hdx_bclus", "dx, cluster with best weighted score;x_{HCAL}-x_{expect} (m)", hbinfac*hcal_vrange, hcalposXi, hcalposXf);
  TH1D *hdy_bclus = new TH1D( "hdy_bclus", "dy, cluster with best weighted score;y_{HCAL}-y_{expect} (m)", hbinfac*hcal_hrange, hcalposYi, hcalposYf);

  TH1D *hdx_stepclus = new TH1D( "hdx_stepclus", "dx, cluster wide coin cut, min dxdy;x_{HCAL}-x_{expect} (m)", hbinfac*hcal_vrange, hcalposXi, hcalposXf);
  TH1D *hdy_stepclus = new TH1D( "hdy_stepclus", "dy, cluster wide coin cut, min dxdy;y_{HCAL}-y_{expect} (m)", hbinfac*hcal_hrange, hcalposYi, hcalposYf);

  TH1D *hdx_intime = new TH1D( "hdx_intime", "dx, cluster wide coin cut, max E;x_{HCAL}-x_{expect} (m)", hbinfac*hcal_vrange, hcalposXi, hcalposXf);
  TH1D *hdy_intime = new TH1D( "hdy_intime", "dy, cluster wide coin cut, max E;y_{HCAL}-y_{expect} (m)", hbinfac*hcal_hrange, hcalposYi, hcalposYf);

  TH1D *hdx_altscore = new TH1D( "hdx_altscore", "dx, cluster BEST (atime PDF val * E/Emax);x_{HCAL}-x_{expect} (m)", hbinfac*hcal_vrange, hcalposXi, hcalposXf);
  TH1D *hdy_altscore = new TH1D( "hdy_altscore", "dy, cluster BEST (atime PDF val * E/Emax);y_{HCAL}-y_{expect} (m)", hbinfac*hcal_hrange, hcalposYi, hcalposYf);

  // TH1D *hW2_nocut = new TH1D( "hW2_nocut", "W^{2}, no cut; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  // TH1D *hW2_pclus = new TH1D( "hW2_pclus", "W^{2}, primary cluster; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  // TH1D *hW2_eclus = new TH1D( "hW2_eclus", "W^{2}, highest e cluster; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  // TH1D *hW2_dclus = new TH1D( "hW2_dclus", "W^{2}, closest dxdy cluster; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  // TH1D *hW2_tclus = new TH1D( "hW2_tclus", "W^{2}, lowest time difference cluster; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );

  //Pass Cluster
  TH1D *hdx_passe = new TH1D( "hdx_passe", "dx, clusters that pass E cut;x_{HCAL}-x_{expect} (m)", hbinfac*hcal_vrange, hcalposXi, hcalposXf);
  TH1D *hdx_passd = new TH1D( "hdx_passd", "dx, clusters that pass dxdy cut;x_{HCAL}-x_{expect} (m)", hbinfac*hcal_vrange, hcalposXi, hcalposXf);
  TH1D *hdx_passt = new TH1D( "hdx_passt", "dx, clusters that pass time diff cut;x_{HCAL}-x_{expect} (m)", hbinfac*hcal_vrange, hcalposXi, hcalposXf);
  TH1D *hdx_passed = new TH1D( "hdx_passed", "dx, clusters that pass E and dxdy cut;x_{HCAL}-x_{expect} (m)", hbinfac*hcal_vrange, hcalposXi, hcalposXf);
  TH1D *hdx_passet = new TH1D( "hdx_passet", "dx, clusters that pass E and time diff cut;x_{HCAL}-x_{expect} (m)", hbinfac*hcal_vrange, hcalposXi, hcalposXf);
  TH1D *hdx_passdt = new TH1D( "hdx_passdt", "dx, clusters that pass dxdy and time diff cut;x_{HCAL}-x_{expect} (m)", hbinfac*hcal_vrange, hcalposXi, hcalposXf);
  TH1D *hdx_passall = new TH1D( "hdx_passall", "dx, clusters that pass all cuts;x_{HCAL}-x_{expect} (m)", hbinfac*hcal_vrange, hcalposXi, hcalposXf);

  TH1D *hW2_passe = new TH1D( "hW2_passe", "W^{2}, clusters that pass E cut; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_passd = new TH1D( "hW2_passd", "W^{2}, clusters that pass dxdy cut; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_passt = new TH1D( "hW2_passt", "W^{2}, clusters that pass time diff cut; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_passed = new TH1D( "hW2_passed", "W^{2}, clusters that pass E and dxdy cut; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_passet = new TH1D( "hW2_passet", "W^{2}, clusters that pass E and time diff cut; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_passdt = new TH1D( "hW2_passdt", "W^{2}, clusters that pass dxdy and time diff cut; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );
  TH1D *hW2_passall = new TH1D( "hW2_passall", "W^{2}, clusters that pass all cuts; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );

  TH1D *hW2_noposcut = new TH1D( "hW2_noposcut", "W^{2}, all cuts, no hcal position; W^{2} (GeV^{2})", binfac*W2fitmax, 0.0, W2fitmax );

  // create output tree for analysis
  TTree *P = new TTree("P","Analysis Tree"); 

  // output tree vars
  double W2_out; //W2
  double Q2_out; //Q2
  double nu_out; //N KE
  double ep_out; //track reconstructed e' momentum
  double hcalpexp_out; //HCal expected momentum
  int failedglobal_out; //earm global variable cuts only
  int failedaccmatch_out; //failed to project elastic nucleon onto HCal
  int cidx_best_atime_out; //cluster index best atime difference
  int cidx_best_dxdy_out; //cluster index best dxdy
  int cidx_best_e_out; //cluster index highest energy
  double bec_atime_out; //best energy cluster atime
  double bec_atimediff_out; //best energy cluster atime diff w/BBCal
  double bec_hcale_out; //best energy cluster energy
  double bec_dx_out; //best energy cluster dx
  double bec_dy_out; //best energy cluster dy
  double btc_atime_out; //best time diff w/bbcal clusters
  double btc_atimediff_out;
  double btc_hcale_out;
  double btc_dx_out;
  double btc_dy_out;
  double bdc_atime_out; //best dxdy clusters
  double bdc_atimediff_out;
  double bdc_hcale_out;
  double bdc_dx_out;
  double bdc_dy_out;
  double bpc_atime_out; //default primary clusters
  double bpc_atimediff_out;
  double bpc_hcale_out;
  double bpc_blk_hcale_out; //primary block energy
  double bpc_dx_out;
  double bpc_dy_out;

  //Overall best cluster set by index
  int cluster_idx = 4; //set to hybrid method by default (1=high E, 2=low dt, 3=low dxdy, 4=hybrid)

  double bc_atime_out; //overall best cluster atime
  double bc_atimediff_out; //overall best  cluster atime diff w/BBCal
  double bc_hcale_out; //overall best  cluster energy
  double bc_dx_out; //overall best  cluster dx
  double bc_dy_out; //overall best  cluster dy

  //cluster tree vars
  int passall_idx_out[maxClus]; //array of indices for clusters that pass all cuts
  double cpblkid_out[maxClus]; //primary block id for each cluster
  double chatime_out[maxClus]; //atime
  double chatimediff_out[maxClus]; //atime diff with bbcal
  double chcale_out[maxClus]; //cluster energy 
  double cdx_out[maxClus]; //
  double cdy_out[maxClus];

  // set output tree branches
  P->Branch( "W2", &W2_out, "W2/D" );
  P->Branch( "Q2", &Q2_out, "Q2/D" );
  P->Branch( "nu", &nu_out, "nu/D" );
  P->Branch( "ep", &ep_out, "ep/D" );
  P->Branch( "hcalpexp", &hcalpexp_out, "hcalpexp/D" );
  P->Branch( "failedglobal", &failedglobal_out, "failedglobal/I" );
  P->Branch( "failedaccmatch", &failedaccmatch_out, "failedaccmatch/I" );
  P->Branch( "cidx_best_atime", &cidx_best_atime_out, "cidx_best_atime/I" );
  P->Branch( "cidx_best_e", &cidx_best_e_out, "cidx_best_e/I" );
  P->Branch( "cidx_best_dxdy", &cidx_best_dxdy_out, "cidx_best_dxdy/I" );
  //best energy cluster vars
  P->Branch( "bec_atime", &bec_atime_out, "bec_atime/D" );
  P->Branch( "bec_atimediff", &bec_atimediff_out, "bec_atimediff/D" );
  P->Branch( "bec_hcale", &bec_hcale_out, "bec_hcale/D" );
  P->Branch( "bec_dx", &bec_dx_out, "bec_dx/D" );
  P->Branch( "bec_dy", &bec_dy_out, "bec_dy/D" );
  //best time cluster vars
  P->Branch( "btc_atime", &btc_atime_out, "btc_atime/D" );
  P->Branch( "btc_atimediff", &btc_atimediff_out, "btc_atimediff/D" );
  P->Branch( "btc_hcale", &btc_hcale_out, "btc_hcale/D" );
  P->Branch( "btc_dx", &btc_dx_out, "btc_dx/D" );
  P->Branch( "btc_dy", &btc_dy_out, "btc_dy/D" );
  //best distance cluster vars
  P->Branch( "bdc_atime", &bdc_atime_out, "bdc_atime/D" );
  P->Branch( "bdc_atimediff", &bdc_atimediff_out, "bdc_atimediff/D" );
  P->Branch( "bdc_hcale", &bdc_hcale_out, "bdc_hcale/D" );
  P->Branch( "bdc_dx", &bdc_dx_out, "bdc_dx/D" );
  P->Branch( "bdc_dy", &bdc_dy_out, "bdc_dy/D" );
  //best primary cluster vars
  P->Branch( "bpc_atime", &bpc_atime_out, "bpc_atime/D" );
  P->Branch( "bpc_atimediff", &bpc_atimediff_out, "bpc_atimediff/D" );
  P->Branch( "bpc_hcale", &bpc_hcale_out, "bpc_hcale/D" );
  P->Branch( "bpc_dx", &bpc_dx_out, "bpc_dx/D" );
  P->Branch( "bpc_dy", &bpc_dy_out, "bpc_dy/D" );
  //best overall cluster vars
  P->Branch( "bc_atime", &bc_atime_out, "bc_atime/D" );
  P->Branch( "bc_atimediff", &bc_atimediff_out, "bc_atimediff/D" );
  P->Branch( "bc_hcale", &bc_hcale_out, "bc_hcale/D" );
  P->Branch( "bc_dx", &bc_dx_out, "bc_dx/D" );
  P->Branch( "bc_dy", &bc_dy_out, "bc_dy/D" );
  //all cluster vars
  P->Branch( "passall_idx", &passall_idx_out, Form("passall_idx[%d]/I",maxClus) );
  P->Branch( "cpblkid", &cpblkid_out, Form("cpblkid[%d]/D",maxClus) );
  P->Branch( "chatime", &chatime_out, Form("chatime[%d]/D",maxClus) );
  P->Branch( "chatimediff", &chatimediff_out, Form("chatimediff[%d]/D",maxClus) );
  P->Branch( "chcale", &chcale_out, Form("chcale[%d]/D",maxClus) );
  P->Branch( "cdx", &cdx_out, Form("cdx[%d]/D",maxClus) );
  P->Branch( "cdy", &cdy_out, Form("cdy[%d]/D",maxClus) );

  //Set up chain
  TChain *C = new TChain("T");
  
  if( diagnostic ){
    for( int r=0; r<sbs4_0p_runs.size(); ++r ){
      string rootfile_full_path = Form("/work/halla/sbs/sbs-gmn/pass0/SBS4/LH2/rootfiles/e1209019_fullreplay_%d*.root",sbs4_0p_runs[r]);
      C->Add(rootfile_full_path.c_str());
    }
  }else if( run == -1 && shorty )
    for( int r=0; r<short_extension_runs.size(); ++r ){
      string rootfile_full_path = rootfile_path + Form("%d*.root",short_extension_runs[r]);
      C->Add(rootfile_full_path.c_str());
    }
  else
    C->Add(rootfile_path.c_str());
  
  TTreeFormula *GlobalCut = new TTreeFormula( "GlobalCut", globalcut, C );
    
  //track vars
  double BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks], BBtr_vz[maxTracks];
  double BBtr_n;

  //bbcal vars
  double BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e, BBsh_atime, BBsh_nclus;

  //hcal vars
  double hcalx, hcaly, hcale, Nclus, atimeblk;
  double hcalcx[maxClus], hcalcy[maxClus], hcalce[maxClus], hcalcbe[maxClus], hcalcatime[maxClus], hcalcid[maxClus];
  int Nhcalcid, Nhcalcbe;

  //switch on the branches
  C ->SetBranchStatus("*",0);
  //hcal
  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus("sbs.hcal.atimeblk");
  C->SetBranchStatus("sbs.hcal.nclus",1);
  C->SetBranchStatus("sbs.hcal.clus.e",1);
  C->SetBranchStatus("sbs.hcal.clus_blk.e",1);
  C->SetBranchStatus("Ndata.sbs.hcal.clus_blk.e",1);
  C->SetBranchStatus("sbs.hcal.clus.atime",1);
  C->SetBranchStatus("sbs.hcal.clus.x",1);
  C->SetBranchStatus("sbs.hcal.clus.y",1);
  C->SetBranchStatus("sbs.hcal.clus.id",1);
  C->SetBranchStatus("Ndata.sbs.hcal.clus.id",1);
  //track
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.p",1);
  C->SetBranchStatus("bb.tr.vz",1);
  //bbcal
  C->SetBranchStatus("bb.ps.e",1);
  C->SetBranchStatus("bb.ps.x",1);
  C->SetBranchStatus("bb.ps.y",1);
  C->SetBranchStatus("bb.sh.e",1);
  C->SetBranchStatus("bb.sh.x",1);
  C->SetBranchStatus("bb.sh.y",1);
  C->SetBranchStatus("bb.sh.atimeblk",1);
  C->SetBranchStatus("bb.sh.nclus",1);

  //link the branches
  //hcal
  C->SetBranchAddress("sbs.hcal.x", &hcalx);
  C->SetBranchAddress("sbs.hcal.y", &hcaly);
  C->SetBranchAddress("sbs.hcal.e", &hcale);
  C->SetBranchAddress("sbs.hcal.atimeblk", &atimeblk);
  C->SetBranchAddress("sbs.hcal.nclus", &Nclus);
  C->SetBranchAddress("sbs.hcal.clus.e", hcalce);
  C->SetBranchAddress("sbs.hcal.clus_blk.e", hcalcbe);
  C->SetBranchAddress("Ndata.sbs.hcal.clus_blk.e", &Nhcalcbe);
  C->SetBranchAddress("sbs.hcal.clus.atime", hcalcatime);
  C->SetBranchAddress("sbs.hcal.clus.x", hcalcx);
  C->SetBranchAddress("sbs.hcal.clus.y", hcalcy);
  C->SetBranchAddress("sbs.hcal.clus.id", hcalcid);
  C->SetBranchAddress("Ndata.sbs.hcal.clus.id", &Nhcalcid);

  //track
  C->SetBranchAddress("bb.tr.n", &BBtr_n);
  C->SetBranchAddress("bb.tr.px", BBtr_px);
  C->SetBranchAddress("bb.tr.py", BBtr_py);
  C->SetBranchAddress("bb.tr.pz", BBtr_pz);
  C->SetBranchAddress("bb.tr.p", BBtr_p);
  C->SetBranchAddress("bb.tr.vz", BBtr_vz);
  //bbcal
  C->SetBranchAddress("bb.ps.e", &BBps_e);
  C->SetBranchAddress("bb.ps.x", &BBps_x);
  C->SetBranchAddress("bb.ps.y", &BBps_y);
  C->SetBranchAddress("bb.sh.e", &BBsh_e);
  C->SetBranchAddress("bb.sh.x", &BBsh_x);
  C->SetBranchAddress("bb.sh.y", &BBsh_y);
  C->SetBranchAddress("bb.sh.atimeblk", &BBsh_atime);
  C->SetBranchAddress("bb.sh.nclus", &BBsh_nclus);

  //Set up HCal coordinate system
  vector<TVector3> hcalaxes;
  TVector3 hcal_zaxis( sin(-HCal_th), 0, cos(-HCal_th) ); // Clock-wise rotation about Y axis
  TVector3 hcal_xaxis( 0, -1, 0 ); // -Y axis of Hall coordinate system = X axis of hcal coordinate system
  TVector3 hcal_yaxis = hcal_zaxis.Cross(hcal_xaxis).Unit();
  hcalaxes.push_back(hcal_xaxis);
  hcalaxes.push_back(hcal_yaxis);
  hcalaxes.push_back(hcal_zaxis);

  TVector3 hcalorigin = HCal_d*hcalaxes[2] + hcalheight*hcalaxes[0];

  // event indices
  long nevent = 0, nevents = C->GetEntries(); 
  int treenum = 0, currenttreenum = 0;

  // elastic sorting markers
  int clusterBadSort = 0;
  int clusterBadAllSort = 0;
  int NElasClusters = 0;

  // all cluster sorting markers
  int allClusterBadAllSort = 0;
  int allClusterBadSort = 0;
  int NClusters = 0;

  while (C->GetEntry(nevent++)) {
      
    cout << "Number of elastic clusters: " << NElasClusters << ". Bad Cluster Sorts: " << clusterBadSort << ". Processing event " << nevent << "/" << nevents << "\r";
    cout.flush();
      
    //Single-loop globalcut method. Save pass/fail for output tree.
    currenttreenum = C->GetTreeNumber();
    if( nevent == 1 || currenttreenum != treenum ){
      treenum = currenttreenum; 
      GlobalCut->UpdateFormulaLeaves();
    }
    bool failedglobal = GlobalCut->EvalInstance(0) == 0;
      
    //E-arm physics calculations
    //correct beam energy with vertex information - leaving out run to run ebeam averaging for this standalone calculation
    double ebeam_c = 0.;
    if( target.compare("LH2") == 0 )
      ebeam_c = E_e - ( (BBtr_vz[0]+l_tgt/2.0) * lh2tarrho * lh2dEdx + lh2uwallthick * rho_al * aldEdx );
    else if( target.compare("LD2") == 0 )
      ebeam_c = E_e - ( (BBtr_vz[0]+l_tgt/2.0) * ld2tarrho * ld2dEdx + ld2uwallthick * rho_al * aldEdx );

    //vertex position
    TVector3 vertex( 0., 0., BBtr_vz[0] );

    //track momentum
    double trackp = BBtr_p[0];

    //set up four-momenta and physics variables per run
    TLorentzVector pbeam( 0., 0., ebeam_c, ebeam_c ); //beam momentum
    //TLorentzVector pe( px[0], py[0], pz[0], p[0] ); //e' plvect
    TLorentzVector pe( trackp*BBtr_px[0]/BBtr_p[0], 
		       trackp*BBtr_py[0]/BBtr_p[0], 
		       trackp*BBtr_pz[0]/BBtr_p[0], 
		       trackp ); //e' recon plvect
    TLorentzVector ptarg;
    if(target.compare("LH2") == 0) 
      ptarg.SetPxPyPzE( 0., 0., 0., M_p );
    else if(target.compare("LD2") == 0) 
      ptarg.SetPxPyPzE( 0., 0., 0., 0.5*(M_n+M_p) );
    TLorentzVector q = pbeam - pe; //virtual photon mom. or mom. transferred to scatter nucleon (N')
    TVector3 qv = q.Vect();
    double etheta = acos( pe.Pz() / pe.E() ); 
    double ephi = atan2( pe.Py(), pe.Px() );
    double pcent;
    if( target.compare("LH2") == 0 ) 
      pcent = E_e/( 1. + ( E_e/M_p )*( 1.0 - cos(etheta) ) );
    else if( target.compare("LD2") == 0 ) {
      double Nmass = 0.5*( M_n + M_p );
      pcent = E_e/( 1. + ( E_e/Nmass )*( 1.0 - cos(etheta) ) );
    }
    double phNexp = ephi + PI;
    TLorentzVector pN; //N' momentum
    //TVector3 pNhat; //Unit N' 3-vector
    double Q2, W2, nu, thNexp, pNexp;

    //Calculate using reconstructed track angles.
    nu = pbeam.E() - pcent;
    if( target.compare("LH2") == 0 ) 
      pNexp = sqrt( pow(nu, 2.) + 2. * M_p * nu );
    else if( target.compare("LD2") == 0 ) 
      pNexp = sqrt( pow(nu, 2.) + 2. * 0.5*(M_n+M_p) * nu );
    thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
    TVector3 pNhat( sin(thNexp) * cos(phNexp), 
		    sin(thNexp) * sin(phNexp), 
		    cos(thNexp) );
    pN.SetPxPyPzE( pNexp*pNhat.X(), pNexp*pNhat.Y(), pNexp*pNhat.Z(), nu+ptarg.E() );
    Q2 = 2.0 * pbeam.E() * pe.E() * ( 1.0 - cos(etheta) );
    if( target.compare("LH2") == 0 ) 
      W2 = pow( M_p, 2.0 ) + 2.0*M_p * (pbeam.E()-pe.E()) - Q2;
    else if( target.compare("LD2") == 0 ) 
      W2 = pow( 0.5*(M_n + M_p), 2.0 ) + 2.0 * 0.5*(M_n + M_p) * (pbeam.E()-pe.E()) - Q2;
    
    //Fill earm physics analysis tree variables and basic histograms
    hW2_nocut->Fill(W2);
    W2_out = W2;
    Q2_out = Q2;
    nu_out = nu;
    ep_out = trackp;
    hcalpexp_out = pNexp;

    /////////////////////
    //W2 elastic cut bool
    bool failedW2 = abs(W2-W2_mean)>W2_sig;

    if(!failedW2)
      hHCALe_W2cut->Fill(hcale);

    //HCal active area cut (acceptance matching). Save pass/fail for output tree.
    vector<double> xyhcalexp; 
    double sintersect = ( hcalorigin - vertex).Dot(hcalaxes[2] ) / ( pNhat.Dot(hcalaxes[2]) );
    TVector3 hcal_intersect = vertex + sintersect*pNhat; 
    double xexpect_hcal = ( hcal_intersect - hcalorigin ).Dot( hcalaxes[0] );
    double yexpect_hcal = ( hcal_intersect - hcalorigin ).Dot( hcalaxes[1] );
    xyhcalexp.push_back( xexpect_hcal );
    xyhcalexp.push_back( yexpect_hcal );

    bool failedaccmatch = 
	xyhcalexp[1] > hcalposYf ||
	xyhcalexp[1] < hcalposYi ||
	xyhcalexp[0] > hcalposXf ||
	xyhcalexp[0] < hcalposXi;

    //Get some primary cluster primary block variables
    double dx = hcalcx[0] - xyhcalexp[0];
    double dy = hcalcy[0] - xyhcalexp[1];
    double atime_pdiff = hcalcatime[0] - BBsh_atime; //Assuming best shower time on primary cluster
    hHCALadiff->Fill(atime_pdiff);

    //dy cut
    bool faileddy = abs(dy-dy_mean)>dy_sig;

    //Fill some output analysis tree variables
    failedglobal_out = (Int_t)failedglobal;
    failedaccmatch_out = (Int_t)failedaccmatch;
    NClusters++;

    // Sort the entries by energy before e-arm cuts for later
    vector<double> cearr;
    for( int c=0; c<Nhcalcid; c++ )
      cearr.push_back(hcalce[c]);

    // Create a vector of indices
    vector<size_t> indices(cearr.size());
    for (size_t i = 0; i < cearr.size(); ++i)
      indices[i] = i;
    
    // Sort the indices based on the values in arr
    std::sort(indices.begin(), indices.end(),
	      [&cearr](size_t i1, size_t i2) { return cearr[i1] < cearr[i2]; });

    // cout << "raw indices post sort: ";
    // for( int i=0; i<indices.size(); ++i )
    //   cout << indices[i] << " ";
    // cout << endl;

    // adjust the sorting markers
    int marker = 0;
    bool not_descending_energy = false;
    int high_e_index = 0;
    for (size_t index : indices) {
      
      int reverse_index = cearr.size()-index-1;
      if( marker == 0 )
	high_e_index = reverse_index;
      //cout << "marker:index:energy " << marker << ":" << reverse_index << ":" << cearr[marker] << endl;
      if( marker==0 && reverse_index!=marker )
	allClusterBadSort++;
      
      if( reverse_index!=marker ){
      	allClusterBadAllSort++;
	not_descending_energy = true;
      	break;
      }
      marker++;
    }
    //cout << endl;

    //primary cluster primary block
    double atimepdiff = atimeblk - BBsh_atime; //Assuming best shower time on primary cluster

    if( !failedglobal && !failedaccmatch && hcale>0.02 && abs(atimepdiff-52)<10 )
      hW2_noposcut->Fill(W2);

    //Make broad elastic cuts now
    if( failedW2 || failedglobal || failedaccmatch )
      continue;

    //cout << endl << "NEVENT (pass elastic): " << nevent << endl << endl;

    //Fill some histograms and increment some markers
    hdx_pclus->Fill(dx);
    hdy_pclus->Fill(dy);
    NElasClusters+=Nhcalcid;
    if( not_descending_energy )
      clusterBadAllSort++;

    //////////////////////
    //ALL CLUSTER ANALYSIS

    //Get appropriate fits
    vector<double> dxfit;
    vector<double> atimefit;

    if( kine==4 ){
      for( int i=0; i<3; ++i){
	dxfit.push_back(dx_fit4[i]);
	atimefit.push_back(atime_fit4[i]);
      }
    }else if( kine==7 ){
      for( int i=0; i<3; ++i){
	dxfit.push_back(dx_fit7[i]);
	atimefit.push_back(atime_fit7[i]);
      }
    }else if( kine==14 ){
      for( int i=0; i<3; ++i){
	dxfit.push_back(dx_fit14[i]);
	atimefit.push_back(atime_fit14[i]);
      }
    }else if( kine==11 ){
      for( int i=0; i<3; ++i){
	dxfit.push_back(dx_fit11[i]);
	atimefit.push_back(atime_fit11[i]);
      }
    }else if( kine==8 ){
      for( int i=0; i<3; ++i){
	dxfit.push_back(dx_fit8[i]);
	atimefit.push_back(atime_fit8[i]);
      }
    }else if( kine==9 ){
      for( int i=0; i<3; ++i){
	dxfit.push_back(dx_fit9[i]);
	atimefit.push_back(atime_fit9[i]);
      }
    }else{
      cout << "ERROR: Enter a valid GMn kinematic." << endl;
      return;
    }

    //Set up indexes and variables for best cluster analysis. All index defaults are zero (highest E)
    int cidx_e = 0;
    double c_e = 0.;    
    int cidx_atime = 0;
    double c_atimediff = 1000.;
    int cidx_dxdy = 0;
    double c_dxdydiff = 1000.;
    int cidx_dxdy_coin = 0;
    double c_dxdydiff_coin = 1000.;
    int passallindex = 0;
    vector<int> passed_coin;
    vector<int> passed_energy;
    vector<int> passed_dxdy;
    double cluster_energy = 0;
    vector<double> clone_cluster_e;
    vector<double> clone_cluster_altscore;

    //loop through the elements of the primary cluster and get total energy for comparison
    for( int b=0; b<Nhcalcbe; ++b )
      cluster_energy += hcalcbe[b];

    //cout << "max energy " << cearr[0] << endl;

    // if( cearr[0] != hcalce[high_e_index] )
    //   cout << "WARNING!: Bad sort on energy " << cearr[0] << ":" << hcalce[high_e_index] << endl;

    //loop through all clusters
    for( int c=0; c<Nhcalcid; c++ ){

      //calculate h-arm physics quantities per cluster
      double atime = hcalcatime[c];
      double atime_diff = atime - BBsh_atime; //Assuming best shower time on primary cluster
      double cdx = hcalcx[c] - xyhcalexp[0];
      double cdy = hcalcy[c] - xyhcalexp[1];
      double ce = hcalce[c];
      TVector3 chcalpos = hcalorigin + hcalcx[c]*hcalaxes[0] + hcalcy[c]*hcalaxes[1];
      double dxdydist = sqrt( pow(cdx-dx_mean,2)+pow(cdy-dy_mean,2));

      //using hcal atime until after pass2, wide cut around 5 sigma
      bool passedCoin = false;
      if( coin==false )
	passedCoin = abs(atime-atimefit[1])<5*atimefit[2];
      else
	passedCoin = abs(atime_diff-atimefit[1])<5*atimefit[2];

      //cluster cuts
      bool c_failedcoin = abs(atime_diff-coin_mean)>coinSigFac*coin_sig;
      bool c_failede = ce<minE;
      bool c_faileddxdy = !IsInsideEllipse(cdx,cdy,dx_mean,dSigFac*dx_sig,dy_mean,dSigFac*dy_sig);

      //Replicate the in-time algorithm
      clone_cluster_e.push_back(ce);
      if( !passedCoin )
	clone_cluster_e[c] = 0;

      //Do alt score (no position info)
      double cascore = 0;
      if( coin==false )
	cascore = assignScoreAlt( ce, atime, hcalce[high_e_index], atimefit);
      else
	cascore = assignScoreAlt( ce, atime_diff, hcalce[high_e_index], atimefit);

      clone_cluster_altscore.push_back(cascore);

      //////////////////////////////////////////////
      //get best cluster indexes for various methods
      //get cluster which maximizes energy (ENERGY METHOD, index 1)
      if( hcalce[c]>c_e ){
	cidx_e = c;
	c_e = ce;
      }
      //get cluster which minimizes adct diff from bbcal (TIME METHOD, index 2)
      if( abs(atime_diff)<c_atimediff ){
	cidx_atime = c;
	c_atimediff = atime_diff;
      }
      //get cluster which minimizes distance between cluster center and expected loc dxdy (POSITION METHOD, index 3)
      if( dxdydist<c_dxdydiff ){
	cidx_dxdy = c;
	c_dxdydiff = dxdydist;
      }
      //get cluster which passes wide cut on atime, then minimizes dxdy for either p or n hypothesis (HYBRID METHOD, index 4)
      if( !c_failedcoin && dxdydist<c_dxdydiff_coin ){
	cidx_dxdy_coin = c;
	c_dxdydiff_coin = dxdydist;
      }

      //Fill cut pass histograms
      if( !c_failede ){
	hdx_passe->Fill(cdx);
	hW2_passe->Fill(W2);
      }
      if( !c_failedcoin ){
	hdx_passt->Fill(cdx);
	hW2_passt->Fill(W2);
      }
      if( !c_faileddxdy ){
	hdx_passd->Fill(cdx);
	hW2_passd->Fill(W2);
      }
      if( !c_failede && !c_failedcoin ){
	hdx_passet->Fill(cdx);
	hW2_passet->Fill(W2);
      }
      if( !c_failede && !c_faileddxdy ){
	hdx_passed->Fill(cdx);
	hW2_passed->Fill(W2);
      }
      if( !c_faileddxdy && !c_failedcoin ){
	hdx_passdt->Fill(cdx);
	hW2_passdt->Fill(W2);
      }
      if( !c_failede && !c_faileddxdy && !c_failedcoin ){
	hdx_passall->Fill(cdx);
	hW2_passall->Fill(W2);
      }

      //fill analysis tree variables
      //cluster cut tree vars
      if( !c_failedcoin && !c_failede && !c_faileddxdy ){
	passall_idx_out[passallindex] = c;
	passallindex++;
      }
      //primary cluster
      if( c==0 ){
	bpc_atime_out = atime;
	bpc_atimediff_out = atime_diff;
	bpc_hcale_out = ce;
	bpc_blk_hcale_out = hcalcbe[0];
	bpc_dx_out = cdx;
	bpc_dy_out = cdy;
      }
      //general
      chcale_out[c] = hcalce[c];
      cpblkid_out[c] = hcalcid[c];
      chatime_out[c] = hcalcatime[c];
      chatimediff_out[c] = hcalcatime[c] - BBsh_atime;
      cdx_out[c] = cdx;
      cdy_out[c] = cdy;
       
    }//endloop over cluster elements

    //Compute cluster dx vars
    double cdx_e = hcalcx[cidx_e] - xyhcalexp[0];
    double cdx_d = hcalcx[cidx_dxdy] - xyhcalexp[0];
    double cdx_t = hcalcx[cidx_atime] - xyhcalexp[0];

    //Errors and warnings to check that sbs.hcal.e and sbs.hcal.clus_blk.e reference the same (highest e) cluster necessarily, but sbs.hcal.clus.e does not necessarily.
    if( hcale != cluster_energy )
      cout << "ERROR: sbs.hcal.e != SUM(sbs.hcal.clus_blk.e[*]) on event " << nevent << endl;
    // if( hcalce[0] != hcale )
    //   cout << "WARNING: sbs.hcal.clus.e != sbs.hcal.e on event " << nevent << endl;

    // for( int i=0; i<clone_cluster_altscore.size(); ++i)
    //   cout << "clone_cluster_altscore[" << i << "] = " << clone_cluster_altscore[i] << endl;
    // cout << endl;

    //Fill best cluster tree variables
    bec_atime_out = hcalcatime[cidx_e];
    bec_atimediff_out = hcalcatime[cidx_e]-BBsh_atime;
    bec_hcale_out = hcalce[cidx_e];
    bec_dx_out = hcalcx[cidx_e] - xyhcalexp[0];
    bec_dy_out = hcalcy[cidx_e] - xyhcalexp[1];

    btc_atime_out = hcalcatime[cidx_atime];
    btc_atimediff_out = hcalcatime[cidx_atime]-BBsh_atime;
    btc_hcale_out = hcalce[cidx_atime];
    btc_dx_out = hcalcx[cidx_atime] - xyhcalexp[0];
    btc_dy_out = hcalcy[cidx_atime] - xyhcalexp[1];

    bdc_atime_out = hcalcatime[cidx_dxdy];
    bdc_atimediff_out = hcalcatime[cidx_dxdy]-BBsh_atime;
    bdc_hcale_out = hcalce[cidx_dxdy];
    bdc_dx_out = hcalcx[cidx_dxdy] - xyhcalexp[0];
    bdc_dy_out = hcalcy[cidx_dxdy] - xyhcalexp[1];

    //Fill best cluster histograms
    hdx_eclus->Fill(cdx_e);
    hdy_eclus->Fill(hcalcy[cidx_e] - xyhcalexp[1]);
    hdx_dclus->Fill(cdx_d);
    hdy_dclus->Fill(hcalcy[cidx_dxdy] - xyhcalexp[1]);
    hdx_tclus->Fill(cdx_t);
    hdy_tclus->Fill(hcalcy[cidx_atime] - xyhcalexp[1]);

    if( cidx_e!=0 )
      clusterBadSort++;

    //Switch between best clusters for systematic analysis
    int cidx_best;
      
    //implement John's ranking scheme with added weights
    vector<double> crude_score;
    vector<double> allweights_score;

    for( int c=0; c<Nhcalcid; ++c ){
      crude_score.push_back(0.);
      allweights_score.push_back(0.);

      if( c == cidx_e ){
      	crude_score[c] += eweight;
      	allweights_score[c] += assignScore(cdx_e, "e", dxfit);

      }
      if( c == cidx_atime ){
      	crude_score[c] += tweight;
      	allweights_score[c] += assignScore(cdx_t, "t", dxfit);

      }
      if( c == cidx_dxdy ){
      	crude_score[c] += dweight;
      	allweights_score[c] += assignScore(cdx_d, "d", dxfit);

      }
    }

    // for( int i=0; i<clone_cluster_e.size(); ++i )
    //   cout << "clone_cluster_e[" << i << "] = " << clone_cluster_e[i] << " ";

    // cout << endl << endl;

    auto [sortedIntimeValues, sortedIntimeIndices] = sortNonZero(clone_cluster_e);
    auto [sortedCrudeValues, sortedCrudeIndices] = sortNonZero(crude_score);
    auto [sortedAllWValues, sortedAllWIndices] = sortNonZero(allweights_score);
    auto [sortedAltValues, sortedAltIndices] = sortNonZero(clone_cluster_altscore);
    
    // for( int i=0; i<sortedIntimeValues.size(); ++i )
    //   cout << sortedIntimeValues[i] << " ";

    // cout << endl << endl << endl;

    //proceed to default
    if( sortedAllWValues[0] == 0 )

    hdx_bclus->Fill( hcalcx[sortedAllWIndices[0]] - xyhcalexp[0]);
    hdy_bclus->Fill( hcalcy[sortedAllWIndices[0]] - xyhcalexp[1]);

    int cidx_intime = sortedIntimeIndices[0];
    int cidx_altscore = sortedAltIndices[0];

    //cout << "sorted intime pidx: " << cidx_intime << " with e: " << sortedIntimeValues[0] << endl;

    switch (cluster_idx) {
    case 1:
      cidx_best = cidx_e; //highest energy cluster
      break;
    case 2:
      cidx_best = cidx_atime; //best coin time cluster
      break;
    case 3:
      cidx_best = cidx_dxdy; //best position cluster
      break;
    case 4:
      cidx_best = cidx_dxdy_coin; //cut coin, then best position cluster
      break;
    case 5:
      cidx_best = cidx_intime; //cut coin, the hightest energy cluster (intime algorithm)
      break;
    default:
      cidx_best = 0; //where 0 is primary cluster on tree
    }
      
    //Best overall cluster variables/histos
    bc_atime_out = hcalcatime[cidx_best];
    bc_atimediff_out = hcalcatime[cidx_best]-BBsh_atime;
    bc_hcale_out = hcalce[cidx_best];
    bc_dx_out = hcalcx[cidx_best] - xyhcalexp[0];
    bc_dy_out = hcalcy[cidx_best] - xyhcalexp[1];
    hdx_stepclus->Fill(hcalcx[cidx_best] - xyhcalexp[0]);
    hdy_stepclus->Fill(hcalcx[cidx_best] - xyhcalexp[1]);
    hdx_intime->Fill(hcalcx[cidx_intime] - xyhcalexp[0]);
    hdy_intime->Fill(hcalcx[cidx_intime] - xyhcalexp[1]);
    hdx_altscore->Fill(hcalcx[cidx_altscore] - xyhcalexp[0]);
    hdy_altscore->Fill(hcalcx[cidx_altscore] - xyhcalexp[1]);

    //Fill final output tree variables    
    cidx_best_atime_out = cidx_atime;
    cidx_best_dxdy_out = cidx_dxdy;
    cidx_best_e_out = cidx_e;

    P->Fill();

  } //end event loop
  
  fout->Write();
  
  cout << endl << "Analysis complete. Outfile located at " << outfile_directory << endl;

  double percent_bad_sorts = (double)clusterBadSort/(double)NElasClusters*100.;
  double percent_bad_sorts_allclus = (double)allClusterBadSort/(double)NClusters*100.;
  double percent_bad_sorts_any = (double)clusterBadAllSort/(double)NElasClusters*100.;
  double percent_bad_sorts_any_allclus = (double)allClusterBadAllSort/(double)NClusters*100.;


  cout << endl << "Total elastic bad sorts / clusters ratio: " << clusterBadSort << "/" << NElasClusters << Form(" (%0.3f%%)",percent_bad_sorts) << endl << endl;
  cout << endl << "Total bad sorts / clusters ratio: " << allClusterBadSort << "/" << NClusters << Form(" (%0.3f%%)",percent_bad_sorts_allclus) << endl << endl;
  cout << endl << "Total elastic bad sorts / clusters ratio (any mismatch): " << clusterBadAllSort << "/" << NElasClusters << Form(" (%0.3f%%)",percent_bad_sorts_any) << endl << endl;
  cout << endl << "Total bad sorts / clusters ratio (any mismatch): " << allClusterBadAllSort << "/" << NClusters << Form(" (%0.3f%%)",percent_bad_sorts_any_allclus) << endl << endl;

}
