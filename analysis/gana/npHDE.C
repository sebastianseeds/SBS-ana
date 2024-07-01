//sseeds 
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
#include <regex>
#include <vector>
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

double binErrLimit = 0.1; //Set limit on HDE error to use HDE prot vs pos exp point for corrections

// Function to concatenate cuts with their names
std::string concatenateCutsWithNames(const std::vector<std::pair<std::string, std::string>>& cutsWithNames) {
    std::string concatenatedCuts;
    for (const auto& cut : cutsWithNames) {
        concatenatedCuts += cut.first + "&&" + cut.second + "&&_____&&";
    }
    return concatenatedCuts;
}

//gets values for optics v x and y cuts for plotting
std::vector<double> extractValues(const std::string& expression) {
  std::vector<double> values;
  std::regex pattern("[<>]=?\\s*(-?\\d*\\.?\\d+)");
  std::smatch matches;

  std::string::const_iterator searchStart(expression.cbegin());
  while (std::regex_search(searchStart, expression.cend(), matches, pattern)) {
    values.push_back(std::stod(matches[1].str()));
    searchStart = matches.suffix().first;
  }

  return values;
}

Double_t SBpol0rej_b = -0.2; // Central-band fit reject begin
Double_t SBpol0rej_e = 0.1; // Central-band fit reject end

// Zeroth order poly (constant) with sideband limits
Double_t sbBGfit_pol0(double *x, double *par) {
    Double_t yint = par[0];

    if (x[0] > SBpol0rej_b && x[0] < SBpol0rej_e) {
      TF1::RejectPoint();
      return 0;
    }

    return yint;
}

//MAIN
void npHDE(int kine=8, 
	   int mag=70, 
	   int pass=2, 
	   bool effz=true) {

  //set draw params
  gStyle->SetPalette(55);
  gStyle->SetCanvasPreferGL(kTRUE);
  gStyle->SetOptFit(111);
  gStyle->SetEndErrorSize(0);
  gStyle->SetOptStat(0110);
  gStyle->SetStatTextColor(kBlack);
  gStyle->SetStatW(0.15);
  gStyle->SetStatH(0.08);
  gStyle->SetStatX(0.9);
  gStyle->SetStatY(0.9);
  
  cout << "Processing analysis of HCal detection efficiency for kinematic " << kine << " at field " << mag << "." << endl << endl;

  // reading json configuration file
  JSONManager *jmgr = new JSONManager("../../config/nphde.json");

  //Get cuts
  std::string magcut = jmgr->GetValueFromSubKey_str( "mag_cuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded magnetic field cuts: " << magcut << endl;

  std::string trackcut = jmgr->GetValueFromSubKey_str( "track_v_cuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded track validity cuts: " << trackcut << endl;

  std::string trackcut_novz = jmgr->GetValueFromSubKey_str( "track_v_cuts_novz", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded track validity cuts, no vertex: " << trackcut_novz << endl;

  std::string trackcut_nochi2 = jmgr->GetValueFromSubKey_str( "track_v_cuts_nochi2", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded track validity cuts, no chi2: " << trackcut_nochi2 << endl;

  std::string optxcut = jmgr->GetValueFromSubKey_str( "optx_v_cuts", Form("sbs%d_%d",kine,mag) );
  vector<double> optxcut_vec = extractValues(optxcut);
  cout << "Loaded optics x validity cuts: " << optxcut << " (values: " << optxcut_vec[0] << ", " << optxcut_vec[1] << ")" << endl;

  std::string optycut = jmgr->GetValueFromSubKey_str( "opty_v_cuts", Form("sbs%d_%d",kine,mag) );
  vector<double> optycut_vec = extractValues(optycut);
  cout << "Loaded optics y validity cuts: " << optycut << " (values: " << optycut_vec[0] << ", " << optycut_vec[1] << ")" << endl;

  std::string fidxcut = jmgr->GetValueFromSubKey_str( "fidx_cut", Form("sbs%d_%d",kine,mag) );
  vector<double> fidxcut_vec = extractValues(fidxcut);
  cout << "Loaded fiducial x cut: " << fidxcut << " (values: " << fidxcut_vec[0] << ", " << fidxcut_vec[1] << ")" << endl;

  std::string fidycut = jmgr->GetValueFromSubKey_str( "fidy_cut", Form("sbs%d_%d",kine,mag) );
  vector<double> fidycut_vec = extractValues(fidycut);
  cout << "Loaded fiducial y cut: " << fidycut << " (values: " << fidycut_vec[0] << ", " << fidycut_vec[1] << ")" << endl;

  std::string earmcut = jmgr->GetValueFromSubKey_str( "earm_cuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded electron arm cuts: " << earmcut << endl;

  std::string earmcut_noW2 = jmgr->GetValueFromSubKey_str( "earm_cuts_noW2", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded electron arm cuts with no W2: " << earmcut_noW2 << endl;

  std::string earmcut_noEp = jmgr->GetValueFromSubKey_str( "earm_cuts_noEp", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded electron arm cuts with no E/p: " << earmcut_noEp << endl;

  std::string earmcut_nopsE = jmgr->GetValueFromSubKey_str( "earm_cuts_nopsE", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded electron arm cuts with no PS E: " << earmcut_nopsE << endl;

  std::string harmcut = jmgr->GetValueFromSubKey_str( "harm_cuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded hadron arm cuts: " << harmcut << endl;

  std::string harmcut_noE = jmgr->GetValueFromSubKey_str( "harm_cuts_noE", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded hadron arm cuts with no HCal E: " << harmcut_noE << endl;

  std::string harmcut_nocoin = jmgr->GetValueFromSubKey_str( "harm_cuts_nocoin", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded hadron arm cuts with no coin: " << harmcut_nocoin << endl;

  std::string harmanticut = jmgr->GetValueFromSubKey_str( "harm_anticuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded hadron arm anticuts: " << harmcut << endl;

  std::string p_dxcut = jmgr->GetValueFromSubKey_str( "h_dp_dx_cuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded dx proton delta cuts: " << p_dxcut << endl;

  std::string n_dxcut = jmgr->GetValueFromSubKey_str( "h_dn_dx_cuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded dx neutron delta cuts: " << n_dxcut << endl;

  std::string p_spotcut = jmgr->GetValueFromSubKey_str( "h_dp_spot_cuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded ellipse proton delta cuts: " << p_spotcut << endl;

  std::string n_spotcut = jmgr->GetValueFromSubKey_str( "h_dn_spot_cuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded ellipse neutron delta cuts: " << n_spotcut << endl;

  std::string p_spotanticut = jmgr->GetValueFromSubKey_str( "h_dp_spot_anticuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded ellipse proton delta cuts: " << p_spotcut << endl;

  std::string n_spotanticut = jmgr->GetValueFromSubKey_str( "h_dn_spot_anticuts", Form("sbs%d_%d",kine,mag) );
  cout << "Loaded ellipse neutron delta cuts: " << n_spotcut << endl;

  //stitch together relevant cuts for later
  std::string nucleon_spotcut = "((" + n_spotcut + ")||(" + p_spotcut + "))";
  std::string nucleon_spotanticut = n_spotanticut + "&&" + p_spotanticut;
  std::string harm_fullanticut = "((" + nucleon_spotanticut + ")||(" + harmanticut + "))";

  std::string tropt_cuts = magcut + "&&" + trackcut + "&&" + optxcut + "&&" + optycut;
  std::string p_exp_cuts = tropt_cuts + "&&" + p_spotcut;
  std::string n_exp_cuts = tropt_cuts + "&&" + n_spotcut;
  std::string nucleon_exp_cuts = tropt_cuts + "&&" + nucleon_spotcut;
  std::string e_cuts = magcut + "&&" + trackcut + "&&" + optxcut + "&&" + optycut + "&&" + earmcut;
  std::string e_cuts_novz = magcut + "&&" + trackcut_novz + "&&" + optxcut + "&&" + optycut + "&&" + earmcut;
  std::string e_cuts_nochi2 = magcut + "&&" + trackcut_nochi2 + "&&" + optxcut + "&&" + optycut + "&&" + earmcut;
  std::string e_cuts_noW2 = magcut + "&&" + trackcut + "&&" + optxcut + "&&" + optycut + "&&" + earmcut_noW2;
  std::string e_cuts_noEp = magcut + "&&" + trackcut + "&&" + optxcut + "&&" + optycut + "&&" + earmcut_noEp;
  std::string e_cuts_nopsE = magcut + "&&" + trackcut + "&&" + optxcut + "&&" + optycut + "&&" + earmcut_nopsE;

  std::string tra_cuts = magcut + "&&" + trackcut;
  std::string opt_cuts = magcut + "&&" + optxcut + "&&" + optycut;

  std::string fid_cuts = fidxcut + "&&" + fidycut;

  //cuts for exracting efficiency
  std::string sig_cuts = e_cuts + "&&" + harmcut + "&&" + fidxcut + "&&" + fidycut + "&&" + nucleon_spotcut;
  std::string norm_cuts = e_cuts + "&&" + harmcut + "&&" + fidxcut + "&&" + fidycut;
  std::string anti_cuts = e_cuts + "&&" + harmcut + "&&" + fidxcut + "&&" + fidycut + "&&" + nucleon_spotanticut;

  //for nucleon of either type
  std::string sig_xcuts = e_cuts + "&&" + harmcut + "&&" + fidycut + "&&" + nucleon_spotcut;
  std::string norm_xcuts = e_cuts + "&&" + fidycut;
  std::string anti_xcuts = e_cuts + "&&" + fidycut + "&&" + harm_fullanticut;

  std::string sig_ycuts = e_cuts + "&&" + harmcut + "&&" + fidxcut + "&&" + nucleon_spotcut;
  std::string norm_ycuts = e_cuts + "&&" + fidxcut;
  std::string anti_ycuts = e_cuts + "&&" + fidxcut + "&&" + harm_fullanticut;

  //for proton (ratio)
  std::string sig_pxcuts = e_cuts + "&&" + harmcut + "&&" + fidycut + "&&" + p_spotcut;
  std::string sig_pycuts = e_cuts + "&&" + harmcut + "&&" + fidxcut + "&&" + p_spotcut;
  std::string sig_pcuts = e_cuts + "&&" + harmcut + "&&" + fidxcut + "&&" + fidycut + "&&" + p_spotcut;
  std::string sig_pcuts_nody = e_cuts + "&&" + harmcut + "&&" + fidxcut + "&&" + fidycut + "&&" + p_dxcut;
  std::string sig_pcuts_nochi2 = e_cuts_nochi2 + "&&" + harmcut + "&&" + fidxcut + "&&" + fidycut + "&&" + p_spotcut;
  std::string sig_pcuts_novz = e_cuts_novz + "&&" + harmcut + "&&" + fidxcut + "&&" + fidycut + "&&" + p_spotcut;
  std::string sig_pcuts_noW2 = e_cuts_noW2 + "&&" + harmcut + "&&" + fidxcut + "&&" + fidycut + "&&" + p_spotcut;
  std::string sig_pcuts_noEp = e_cuts_noEp + "&&" + harmcut + "&&" + fidxcut + "&&" + fidycut + "&&" + p_spotcut;
  std::string sig_pcuts_nopsE = e_cuts_nopsE + "&&" + harmcut + "&&" + fidxcut + "&&" + fidycut + "&&" + p_spotcut;
  std::string sig_pcuts_noE = e_cuts + "&&" + harmcut_noE + "&&" + fidxcut + "&&" + fidycut + "&&" + p_spotcut;
  std::string sig_pcuts_nocoin = e_cuts + "&&" + harmcut_nocoin + "&&" + fidxcut + "&&" + fidycut + "&&" + p_spotcut;

  //for neutron (ratio)
  std::string sig_nxcuts = e_cuts + "&&" + harmcut + "&&" + fidycut + "&&" + n_spotcut;
  std::string sig_nycuts = e_cuts + "&&" + harmcut + "&&" + fidxcut + "&&" + n_spotcut;
  std::string sig_ncuts = e_cuts + "&&" + harmcut + "&&" + fidxcut + "&&" + fidycut + "&&" + n_spotcut;
  std::string sig_ncuts_nody = e_cuts + "&&" + harmcut + "&&" + fidxcut + "&&" + fidycut + "&&" + n_dxcut;
  std::string sig_ncuts_nochi2 = e_cuts_nochi2 + "&&" + harmcut + "&&" + fidxcut + "&&" + fidycut + "&&" + n_spotcut;
  std::string sig_ncuts_novz = e_cuts_novz + "&&" + harmcut + "&&" + fidxcut + "&&" + fidycut + "&&" + n_spotcut;
  std::string sig_ncuts_noW2 = e_cuts_noW2 + "&&" + harmcut + "&&" + fidxcut + "&&" + fidycut + "&&" + n_spotcut;
  std::string sig_ncuts_noEp = e_cuts_noEp + "&&" + harmcut + "&&" + fidxcut + "&&" + fidycut + "&&" + n_spotcut;
  std::string sig_ncuts_nopsE = e_cuts_nopsE + "&&" + harmcut + "&&" + fidxcut + "&&" + fidycut + "&&" + n_spotcut;
  std::string sig_ncuts_noE = e_cuts + "&&" + harmcut_noE + "&&" + fidxcut + "&&" + fidycut + "&&" + n_spotcut;
  std::string sig_ncuts_nocoin = e_cuts + "&&" + harmcut_nocoin + "&&" + fidxcut + "&&" + fidycut + "&&" + n_spotcut;

  //Get plot limits/bins
  int bbmidxbins = jmgr->GetValueFromSubKey<int>( "bbbinsx", Form("sbs%d_%d",kine,mag) );
  int bbmidybins = jmgr->GetValueFromSubKey<int>( "bbbinsy", Form("sbs%d_%d",kine,mag) );
  int hcalxbins = jmgr->GetValueFromSubKey<int>( "hbinsx", Form("sbs%d_%d",kine,mag) );
  int hcalybins = jmgr->GetValueFromSubKey<int>( "hbinsy", Form("sbs%d_%d",kine,mag) );
  int W2bins = jmgr->GetValueFromSubKey<int>( "W2bins", Form("sbs%d_%d",kine,mag) );
  int Epbins = jmgr->GetValueFromSubKey<int>( "Epbins", Form("sbs%d_%d",kine,mag) );
  int psEbins = jmgr->GetValueFromSubKey<int>( "psEbins", Form("sbs%d_%d",kine,mag) );
  int Ebins = jmgr->GetValueFromSubKey<int>( "Ebins", Form("sbs%d_%d",kine,mag) );
  int coinbins = jmgr->GetValueFromSubKey<int>( "coinbins", Form("sbs%d_%d",kine,mag) );
  int dybins = jmgr->GetValueFromSubKey<int>( "dybins", Form("sbs%d_%d",kine,mag) );
  int vzbins = jmgr->GetValueFromSubKey<int>( "vzbins", Form("sbs%d_%d",kine,mag) );
  int chi2bins = jmgr->GetValueFromSubKey<int>( "chi2bins", Form("sbs%d_%d",kine,mag) );

  std::pair<double,double> bbmidxlim = {
    jmgr->GetValueFromSubKey<double>( "bbmidx_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "bbmidx_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded bbmidx plot limits min: " << bbmidxlim.first << ", max: " << bbmidxlim.second << ", with bins: " << bbmidxbins << endl;

  std::pair<double,double> bbmidylim = {
    jmgr->GetValueFromSubKey<double>( "bbmidy_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "bbmidy_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded bbmidy plot limits min: " << bbmidylim.first << ", max: " << bbmidylim.second << ", with bins: " << bbmidybins <<  endl;

  std::pair<double,double> hcalxlim = {
    jmgr->GetValueFromSubKey<double>( "hcalx_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "hcalx_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded hcalx plot limits min: " << hcalxlim.first << ", max: " << hcalxlim.second << ", with bins: " << hcalxbins <<  endl;

  std::pair<double,double> hcalylim = {
    jmgr->GetValueFromSubKey<double>( "hcaly_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "hcaly_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded hcaly plot limits min: " << hcalylim.first << ", max: " << hcalylim.second << ", with bins: " << hcalybins <<  endl;

  std::pair<double,double> W2lim = {
    jmgr->GetValueFromSubKey<double>( "W2_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "W2_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded W2 plot limits min: " << W2lim.first << ", max: " << W2lim.second << ", with bins: " << W2bins <<  endl;

  std::pair<double,double> Eplim = {
    jmgr->GetValueFromSubKey<double>( "Ep_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "Ep_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded e/p plot limits min: " << Eplim.first << ", max: " << Eplim.second << ", with bins: " << Epbins <<  endl;

  std::pair<double,double> psElim = {
    jmgr->GetValueFromSubKey<double>( "psE_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "psE_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded PS E plot limits min: " << psElim.first << ", max: " << psElim.second << ", with bins: " << psEbins <<  endl;

  std::pair<double,double> Elim = {
    jmgr->GetValueFromSubKey<double>( "E_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "E_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded HCal E plot limits min: " << Elim.first << ", max: " << Elim.second << ", with bins: " << Ebins <<  endl;

  std::pair<double,double> coinlim = {
    jmgr->GetValueFromSubKey<double>( "coin_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "coin_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded coin plot limits min: " << coinlim.first << ", max: " << coinlim.second << ", with bins: " << coinbins <<  endl;

  std::pair<double,double> dylim = {
    jmgr->GetValueFromSubKey<double>( "dy_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "dy_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded dy plot limits min: " << dylim.first << ", max: " << dylim.second << ", with bins: " << dybins <<  endl;

  std::pair<double,double> vzlim = {
    jmgr->GetValueFromSubKey<double>( "vz_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "vz_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded vz plot limits min: " << vzlim.first << ", max: " << vzlim.second << ", with bins: " << vzbins <<  endl;

  std::pair<double,double> chi2lim = {
    jmgr->GetValueFromSubKey<double>( "chi2_llim", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "chi2_ulim", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded chi2 plot limits min: " << chi2lim.first << ", max: " << chi2lim.second << ", with bins: " << chi2bins <<  endl;
 
  //Get fit ranges
  std::pair<double,double> W2fit = {
    jmgr->GetValueFromSubKey<double>( "W2_lfit", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "W2_ufit", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded W2 plot fit limits min: " << W2fit.first << ", max: " << W2fit.second << endl;

  std::pair<double,double> Epfit = {
    jmgr->GetValueFromSubKey<double>( "Ep_lfit", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "Ep_ufit", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded E/p plot fit limits min: " << Epfit.first << ", max: " << Epfit.second << endl;

  std::pair<double,double> psEfit = {
    jmgr->GetValueFromSubKey<double>( "psE_lfit", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "psE_ufit", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded PS E plot fit limits min: " << psEfit.first << ", max: " << psEfit.second << endl;

  std::pair<double,double> Efit = {
    jmgr->GetValueFromSubKey<double>( "E_lfit", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "E_ufit", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded HCal E plot fit limits min: " << Efit.first << ", max: " << Efit.second << endl;

  std::pair<double,double> coinfit = {
    jmgr->GetValueFromSubKey<double>( "coin_lfit", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "coin_ufit", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded coin plot fit limits min: " << coinfit.first << ", max: " << coinfit.second << endl;

  std::pair<double,double> dyfit = {
    jmgr->GetValueFromSubKey<double>( "dy_lfit", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "dy_ufit", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded dy plot fit limits min: " << dyfit.first << ", max: " << dyfit.second << endl;

  std::pair<double,double> vzfit = {
    jmgr->GetValueFromSubKey<double>( "vz_lfit", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "vz_ufit", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded vz plot fit limits min: " << vzfit.first << ", max: " << vzfit.second << endl;

  std::pair<double,double> chi2fit = {
    jmgr->GetValueFromSubKey<double>( "chi2_lfit", Form("sbs%d_%d",kine,mag) ), 
    jmgr->GetValueFromSubKey<double>( "chi2_ufit", Form("sbs%d_%d",kine,mag) )
  };
  cout << "Loaded chi2 plot fit limits min: " << chi2fit.first << ", max: " << chi2fit.second << endl;


  //set up files and paths
  std::string effz_word = "";
  if(effz)
    effz_word = "_effz";
  std::string tar_word = "_ld2";
  std::string mc_word = "";

  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fin_path = outdir_path + Form("/parse/parse_sbs%d_pass%d_barebones%s%s.root",kine,pass,effz_word.c_str(),tar_word.c_str());

  std::string fout_path = outdir_path + Form("/gmn_analysis/HDE/npHDE_sbs%d_mag%d_pass%d%s%s_alt2.root",kine,mag,pass,effz_word.c_str(),mc_word.c_str());

  cout << "Setting up output path: " << fout_path << endl;


  // Open the ROOT file
  TFile* inputFile = new TFile(fin_path.c_str());
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "Error opening file: " << fin_path << std::endl;
    return;
  }

  // Get the tree
  TTree* tree = dynamic_cast<TTree*>(inputFile->Get("P"));
  if (!tree) {
    std::cerr << "Tree not found in file: " << fin_path << std::endl;
    inputFile->Close();
    return;
  }

  // Open an output ROOT file
  TFile* outputFile = new TFile(fout_path.c_str(), "RECREATE");

  ///////////////////////
  ///////////////////////
  ///////////////////////
  ///////////////////////
  ///////////////////////
  ///////////////////////
  //Begin Plots
  ///////////////////////
  ///////////////////////
  ///////////////////////

  // Create a new canvas
  TCanvas *cW2optx = new TCanvas("cW2optx", "W2 vs Optics x", 1000, 600);
  cW2optx->cd();

  // Draw vertical projection to bb midplane with track v cut
  std::string W2optxName = "hW2vOptx";
  std::string W2optxcut = trackcut + "&&" + magcut;
  std::string W2optxTitle = "W^{2} vs BB midplane x (Track Validity Cuts) ";
  TH2D* W2optx = new TH2D( W2optxName.c_str(), 
			   (W2optxTitle + ";bb_tr_r_x-0.9*bb_tr_r_th; GeV^{2}").c_str(), 
			   bbmidxbins, 
			   bbmidxlim.first, 
			   bbmidxlim.second, 
			   W2bins, 
			   W2lim.first, 
			   W2lim.second);

  // Draw the plot using the created histogram
  tree->Draw(("W2:(bb_tr_r_x-0.9*bb_tr_r_th)>>" + W2optxName).c_str(), W2optxcut.c_str(), "COLZ");

  // Draw the histogram on the canvas
  W2optx->Draw("COLZ");

  // Get cuts from json and make TLines
  TLine *xline1 = new TLine(optxcut_vec[0], W2lim.first, optxcut_vec[0], W2lim.second);
  TLine *xline2 = new TLine(optxcut_vec[1], W2lim.first, optxcut_vec[1], W2lim.second);

  xline1->SetLineColor(kRed);
  xline1->SetLineStyle(2);   
  xline1->SetLineWidth(2);
  xline2->SetLineColor(kRed);
  xline2->SetLineStyle(2); 
  xline2->SetLineWidth(2); 

  xline1->Draw("same");
  xline2->Draw("same");

  // Write the histogram to the output file
  cW2optx->Update();
  cW2optx->Write();

  // Create a new canvas
  TCanvas *cW2opty = new TCanvas("cW2opty", "W2 vs Optics y", 1000, 600);
  cW2opty->cd();

  // Draw vertical projection to bb midplane with track v cut vertical projection cut
  std::string W2optyName = "hW2vOpty";
  std::string W2optycut = W2optxcut + "&&" + optxcut;
  std::string W2optyTitle = "W^{2} vs BB midplane y (Track and Optical x Validity Cuts)";
  TH2D* W2opty = new TH2D( W2optyName.c_str(), 
			   (W2optyTitle + ";bb_tr_r_y-0.9*bb_tr_r_ph; GeV^{2}").c_str(), 
			   bbmidybins, 
			   bbmidylim.first, 
			   bbmidylim.second, 
			   W2bins, 
			   W2lim.first, 
			   W2lim.second);

  // Draw the plot using the created histogram
  tree->Draw(("W2:(bb_tr_r_y-0.9*bb_tr_r_ph)>>" + W2optyName).c_str(), W2optycut.c_str(), "COLZ");

  // Draw the histogram on the canvas
  W2opty->Draw("COLZ");

  // Get cuts from json and make TLines
  TLine *yline1 = new TLine(optycut_vec[0], W2lim.first, optycut_vec[0], W2lim.second);
  TLine *yline2 = new TLine(optycut_vec[1], W2lim.first, optycut_vec[1], W2lim.second);

  yline1->SetLineColor(kRed);
  yline1->SetLineStyle(2);   
  yline1->SetLineWidth(2);
  yline2->SetLineColor(kRed);
  yline2->SetLineStyle(2); 
  yline2->SetLineWidth(2); 

  yline1->Draw("same");
  yline2->Draw("same");

  // Write the histogram to the output file
  cW2opty->Update();
  cW2opty->Write();

  // Create a new canvas
  TCanvas *cp_xyexp = new TCanvas("cp_xyexp", "xexp vs yexp", 1200, 600);
  cp_xyexp->Divide(2,1);
  cp_xyexp->cd(1);

  // Draw the xexp heat maps
  std::string p_xyexpName = "hp_xyexp";
  std::string p_xyexpcuts = tropt_cuts + "&&" + earmcut;
  std::string p_xyexpTitle = "x_{exp} vs y_{exp} (Validity and Elastic Cuts)";
  TH2D* hp_xyexp = new TH2D( p_xyexpName.c_str(), 
			   (p_xyexpTitle + ";y_{exp} (m); x_{exp} (m)").c_str(), 
			   hcalybins, 
			   hcalylim.first, 
			   hcalylim.second, 
			   hcalxbins, 
			   hcalxlim.first, 
			   hcalxlim.second);

  // Draw the plot using the created histogram
  tree->Draw(("xexp:yexp>>" + p_xyexpName).c_str(), p_xyexpcuts.c_str(), "COLZ");

  std::string p_xyexpspotName = "hp_xyexpspot";
  std::string p_xyexpspotcuts = tropt_cuts + "&&" + earmcut + "&&" + p_spotcut;
  std::string p_xyexpspotTitle = "x_{exp} vs y_{exp} (Validity, Elastic, Proton Spot Cuts)";
  TH2D* hp_xyexpspot = new TH2D( p_xyexpspotName.c_str(), 
			       (p_xyexpspotTitle + ";y_{exp} (m); x_{exp} (m)").c_str(), 
			       hcalybins, 
			       hcalylim.first, 
			       hcalylim.second, 
			       hcalxbins, 
			       hcalxlim.first, 
			       hcalxlim.second);

  // Draw the plot using the created histogram
  tree->Draw(("xexp:yexp>>" + p_xyexpspotName).c_str(), p_xyexpspotcuts.c_str(), "COLZ");

  // Bottom line (horizontal) at y = fidxcut_vec[0]
  TLine *xexpline1 = new TLine(fidycut_vec[0], fidxcut_vec[0], fidycut_vec[1], fidxcut_vec[0]);

  // Bottom line (horizontal) at y = fidxcut_vec[0]
  TLine *xexpline2 = new TLine(fidycut_vec[0], fidxcut_vec[1], fidycut_vec[1], fidxcut_vec[1]);

  // Left side line (vertical) at x = fidycut_vec[0]
  TLine *yexpline1 = new TLine(fidycut_vec[0], fidxcut_vec[0], fidycut_vec[0], fidxcut_vec[1]);

  // Right side line (vertical) at x = fidycut_vec[1]
  TLine *yexpline2 = new TLine(fidycut_vec[1], fidxcut_vec[0], fidycut_vec[1], fidxcut_vec[1]);

  // Draw the histogram on the canvas
  int p_xyexp_nev = hp_xyexp->GetEntries();
  hp_xyexp->Draw("COLZ");

  xexpline1->SetLineColor(kRed);
  xexpline1->SetLineStyle(2);   
  xexpline1->SetLineWidth(2);
  xexpline2->SetLineColor(kRed);
  xexpline2->SetLineStyle(2);   
  xexpline2->SetLineWidth(2);

  xexpline1->Draw("same");
  xexpline2->Draw("same");

  yexpline1->SetLineColor(kRed);
  yexpline1->SetLineStyle(2);   
  yexpline1->SetLineWidth(2);
  yexpline2->SetLineColor(kRed);
  yexpline2->SetLineStyle(2); 
  yexpline2->SetLineWidth(2); 

  yexpline1->Draw("same");
  yexpline2->Draw("same");

  TLegend *leg_exp = new TLegend(0.11, 0.7, 0.5, 0.89);
  leg_exp->AddEntry((TObject*)0, Form("N Ev: %d", p_xyexp_nev), "");
  leg_exp->AddEntry(xexpline1, "Fiducial cuts", "l");
  leg_exp->Draw();

  cp_xyexp->cd(2);

  // Draw the histogram on the canvas
  int p_xyexpspot_nev = hp_xyexpspot->GetEntries();
  hp_xyexpspot->Draw("COLZ");

  xexpline1->Draw("same");
  xexpline2->Draw("same");

  yexpline1->Draw("same");
  yexpline2->Draw("same");

  TLegend *legp_expspot = new TLegend(0.11, 0.7, 0.5, 0.89);
  legp_expspot->AddEntry((TObject*)0, Form("N Ev: %d", p_xyexpspot_nev), "");
  legp_expspot->AddEntry(xexpline1, "Fiducial cuts", "l");
  legp_expspot->Draw(); 

  // Write the canvas to the output file
  cp_xyexp->Update();
  cp_xyexp->Write();

  // Create a new canvas for neutron expected x v y
  TCanvas *cn_xyexp = new TCanvas("cn_xyexp", "xexp vs yexp", 1200, 600);
  cn_xyexp->Divide(2,1);
  cn_xyexp->cd(1);

  // Draw the xexp heat maps
  std::string n_xyexpName = "hn_xyexp";
  std::string n_xyexpcuts = tropt_cuts + "&&" + earmcut;
  std::string n_xyexpTitle = "x_{exp} vs y_{exp} (Validity and Elastic Cuts)";
  TH2D* hn_xyexp = new TH2D( n_xyexpName.c_str(), 
			   (n_xyexpTitle + ";y_{exp} (m); x_{exp} (m)").c_str(), 
			   hcalybins, 
			   hcalylim.first, 
			   hcalylim.second, 
			   hcalxbins, 
			   hcalxlim.first, 
			   hcalxlim.second);

  // Draw the plot using the created histogram
  tree->Draw(("xexp:yexp>>" + n_xyexpName).c_str(), n_xyexpcuts.c_str(), "COLZ");

  std::string n_xyexpspotName = "hn_xyexpspot";
  std::string n_xyexpspotcuts = tropt_cuts + "&&" + earmcut + "&&" + n_spotcut;
  std::string n_xyexpspotTitle = "x_{exp} vs y_{exp} (Validity, Elastic, Neutron Spot Cuts)";
  TH2D* hn_xyexpspot = new TH2D( n_xyexpspotName.c_str(), 
			       (n_xyexpspotTitle + ";y_{exp} (m); x_{exp} (m)").c_str(), 
			       hcalybins, 
			       hcalylim.first, 
			       hcalylim.second, 
			       hcalxbins, 
			       hcalxlim.first, 
			       hcalxlim.second);

  // Draw the plot using the created histogram
  tree->Draw(("xexp:yexp>>" + n_xyexpspotName).c_str(), n_xyexpspotcuts.c_str(), "COLZ");

  // Draw the histogram on the canvas
  int n_xyexp_nev = hn_xyexp->GetEntries();
  hn_xyexp->Draw("COLZ");

  xexpline1->Draw("same");
  xexpline2->Draw("same");

  yexpline1->Draw("same");
  yexpline2->Draw("same");

  TLegend *legn_exp = new TLegend(0.11, 0.7, 0.5, 0.89);
  legn_exp->AddEntry((TObject*)0, Form("N Ev: %d", n_xyexp_nev), "");
  legn_exp->AddEntry(xexpline1, "Fiducial cuts", "l");
  legn_exp->Draw();

  cn_xyexp->cd(2);

  // Draw the histogram on the canvas
  int n_xyexpspot_nev = hn_xyexpspot->GetEntries();
  hn_xyexpspot->Draw("COLZ");

  xexpline1->Draw("same");
  xexpline2->Draw("same");

  yexpline1->Draw("same");
  yexpline2->Draw("same");

  TLegend *leg_expspot = new TLegend(0.11, 0.7, 0.5, 0.89);
  leg_expspot->AddEntry((TObject*)0, Form("N Ev: %d", n_xyexpspot_nev), "");
  leg_expspot->AddEntry(xexpline1, "Fiducial cuts", "l");
  leg_expspot->Draw(); 

  // Write the canvas to the output file
  cn_xyexp->Update();
  cn_xyexp->Write();

  // Create a new canvas
  TCanvas *cW2fidx_p = new TCanvas("cW2fidx_p", "W2 vs xexp proton", 1000, 600);
  cW2fidx_p->cd();

  // Draw the W2 fiducial x cut histogram
  std::string W2fidx_pName = "hW2vFidx_P";
  std::string W2fidx_pTitle = "W^{2} vs x_{exp} (Validity and Proton Spot Cuts)";
  TH2D* W2fidx_p = new TH2D( W2fidx_pName.c_str(), 
			   (W2fidx_pTitle + ";x_{exp} (m); GeV^{2}").c_str(), 
			   hcalxbins, 
			   hcalxlim.first, 
			   hcalxlim.second, 
			   W2bins, 
			   W2lim.first, 
			   W2lim.second);

  // Draw the plot using the created histogram
  tree->Draw(("W2:xexp>>" + W2fidx_pName).c_str(), p_exp_cuts.c_str(), "COLZ");

  // Draw the histogram on the canvas
  W2fidx_p->Draw("COLZ");

  // Get cuts from json and make TLines
  TLine *xfidline1 = new TLine(fidxcut_vec[0], W2lim.first, fidxcut_vec[0], W2lim.second);
  TLine *xfidline2 = new TLine(fidxcut_vec[1], W2lim.first, fidxcut_vec[1], W2lim.second);

  xfidline1->SetLineColor(kRed);
  xfidline1->SetLineStyle(2);   
  xfidline1->SetLineWidth(2);
  xfidline2->SetLineColor(kRed);
  xfidline2->SetLineStyle(2);   
  xfidline2->SetLineWidth(2);

  xfidline1->Draw("same");
  xfidline2->Draw("same");

  // Write the histogram to the output file
  cW2fidx_p->Update();
  cW2fidx_p->Write();

  // Create a new canvas
  TCanvas *cW2fidy_p = new TCanvas("cW2fidy_p", "W2 vs yexp", 1000, 600);
  cW2fidy_p->cd();

  // Draw the fiducial y cut histogram
  std::string W2fidy_pName = "hW2vFidy_P";
  std::string W2fidy_pTitle = "W^{2} vs y_{exp} (Validity and Proton Spot Cuts)";
  TH2D* W2fidy_p = new TH2D( W2fidy_pName.c_str(), 
			   (W2fidy_pTitle + ";y_{exp} (m); GeV^{2}").c_str(), 
			   hcalybins, 
			   hcalylim.first, 
			   hcalylim.second, 
			   W2bins, 
			   W2lim.first, 
			   W2lim.second);

  // Draw the plot using the created histogram
  tree->Draw(("W2:yexp>>" + W2fidy_pName).c_str(), p_exp_cuts.c_str(), "COLZ");

  // Draw the histogram on the canvas
  W2fidy_p->Draw("COLZ");

  // Get cuts from json and make TLines
  TLine *yfidline1 = new TLine(fidycut_vec[0], W2lim.first, fidycut_vec[0], W2lim.second);
  TLine *yfidline2 = new TLine(fidycut_vec[1], W2lim.first, fidycut_vec[1], W2lim.second);

  yfidline1->SetLineColor(kRed);
  yfidline1->SetLineStyle(2);   
  yfidline1->SetLineWidth(2);
  yfidline2->SetLineColor(kRed);
  yfidline2->SetLineStyle(2); 
  yfidline2->SetLineWidth(2); 

  yfidline1->Draw("same");
  yfidline2->Draw("same");

  // Write the histogram to the output file
  cW2fidy_p->Update();
  cW2fidy_p->Write();

  // Create a new canvas for neutron
  TCanvas *cW2fidx_n = new TCanvas("cW2fidx_n", "W2 vs xexp neutron", 1000, 600);
  cW2fidx_n->cd();

  // Draw the W2 fiducial x cut histogram
  std::string W2fidx_nName = "hW2vFidx_N";
  std::string W2fidx_nTitle = "W^{2} vs x_{exp} (Validity and Neutron Spot Cuts)";
  TH2D* W2fidx_n = new TH2D( W2fidx_nName.c_str(), 
			   (W2fidx_nTitle + ";x_{exp} (m); GeV^{2}").c_str(), 
			   hcalxbins, 
			   hcalxlim.first, 
			   hcalxlim.second, 
			   W2bins, 
			   W2lim.first, 
			   W2lim.second);

  // Draw the plot using the created histogram
  tree->Draw(("W2:xexp>>" + W2fidx_nName).c_str(), n_exp_cuts.c_str(), "COLZ");

  // Draw the histogram on the canvas
  W2fidx_n->Draw("COLZ");

  // Get cuts from json and make TLines
  xfidline1->Draw("same");
  xfidline2->Draw("same");

  // Write the histogram to the output file
  cW2fidx_n->Update();
  cW2fidx_n->Write();

  // Create a new canvas
  TCanvas *cW2fidy_n = new TCanvas("cW2fidy_n", "W2 vs yexp", 1000, 600);
  cW2fidy_n->cd();

  // Draw the fiducial y cut histogram
  std::string W2fidy_nName = "hW2vFidy_N";
  std::string W2fidy_nTitle = "W^{2} vs y_{exp} (Validity and Neutron Spot Cuts)";
  TH2D* W2fidy_n = new TH2D( W2fidy_nName.c_str(), 
			   (W2fidy_nTitle + ";y_{exp} (m); GeV^{2}").c_str(), 
			   hcalybins, 
			   hcalylim.first, 
			   hcalylim.second, 
			   W2bins, 
			   W2lim.first, 
			   W2lim.second);

  // Draw the plot using the created histogram
  tree->Draw(("W2:yexp>>" + W2fidy_nName).c_str(), n_exp_cuts.c_str(), "COLZ");

  // Draw the histogram on the canvas
  W2fidy_n->Draw("COLZ");

  // Get cuts from json and make TLines
  yfidline1->Draw("same");
  yfidline2->Draw("same");

  // Write the histogram to the output file
  cW2fidy_n->Update();
  cW2fidy_n->Write();

  // Create a new canvas
  TCanvas *cdxdy = new TCanvas("cdxdy", "dx vs dy", 1200, 600);
  cdxdy->Divide(2,1);
  cdxdy->cd(1);

  // Draw dx vs dy plot for comparison with spot cut region
  std::string dxdyName = "hdxdy";
  std::string dxdycut = norm_cuts;
  std::string dxdyTitle = "#Delta x vs #Delta y (Validity, Fiducial, and Elastic Cuts)";
  TH2D* hdxdy = new TH2D( dxdyName.c_str(), 
			  (dxdyTitle + ";#Delta y (m); #Delta x (m)").c_str(), 
			  hcalybins, 
			  hcalylim.first, 
			  hcalylim.second, 
			  hcalxbins, 
			  hcalxlim.first, 
			  hcalxlim.second);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("dx_bc:dy_bc>>" + dxdyName).c_str(), dxdycut.c_str(), "COLZ");

  // Draw dx vs dy plot for comparison with spot cut region
  std::string dxdyspotName = "hdxdyspot";
  std::string dxdyspotcut = sig_cuts;
  std::string dxdyspotTitle = "#Delta x vs #Delta y (Validity, Fiducial, Elastic, and Nucleon Spot Cuts)";
  TH2D* hdxdyspot = new TH2D( dxdyspotName.c_str(),
			      (dxdyspotTitle + ";#Delta y (m); #Delta x (m)").c_str(), 
			      hcalybins, 
			      hcalylim.first, 
			      hcalylim.second, 
			      hcalxbins, 
			      hcalxlim.first, 
			      hcalxlim.second);

  // Draw the plot using the created histogram
  tree->Draw(("dx_bc:dy_bc>>" + dxdyspotName).c_str(), dxdyspotcut.c_str(), "COLZ");

  // Draw the dxdy plot without spot cut
  hdxdy->Draw("COLZ");

  // Draw the dxdy plot with the spot cut
  cdxdy->cd(2);
  hdxdyspot->Draw("colz");

  cdxdy->Update();
  cdxdy->Write();

  // Create a new canvas
  TCanvas *cmclus = new TCanvas("cmclus", "multicluster cut comparison", 1200, 600);
  cmclus->Divide(2,1);
  cmclus->cd(1);

  // Draw associated cluster plot (centroids within 60 cm radius) with e-arm and spot cuts
  std::string mclusName = "hmclus";
  std::string mcluscut = sig_cuts;
  std::string mclusTitle = "#Delta x vs #Delta y (Validity, Fiducial, Elastic, and Spot Cuts)";
  TH1D* hmclus = new TH1D( mclusName.c_str(), 
			  (mclusTitle + ";N Associated Clusters").c_str(), 
			  5, 
			  0, 
			  5);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("hcalnclus_mclus>>" + mclusName).c_str(), mcluscut.c_str(), "COLZ");

  // Draw associated cluster plot (centroids within 60 cm radius) with e-arm and antispot cuts
  std::string mclusantiName = "hmclusanti";
  std::string mclusanticut = anti_cuts;
  std::string mclusantiTitle = "#Delta x vs #Delta y (Validity, Fiducial, Elastic, and Antispot Cuts)";
  TH1D* hmclusanti = new TH1D( mclusantiName.c_str(), 
			  (mclusantiTitle + ";N Associated Clusters").c_str(), 
			  5, 
			  0, 
			  5);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("hcalnclus_mclus>>" + mclusantiName).c_str(), mclusanticut.c_str(), "COLZ");

  // Draw dx vs dy plot for comparison with spot cut region
  std::string dxdymclusName = "hdxdymclus";
  std::string dxdymcluscut = norm_cuts + "&&hcalnclus_mclus>1";
  std::string dxdymclusTitle = "#Delta x vs #Delta y (Validity, Fiducial, Elastic, and Mclus Cuts)";
  TH2D* hdxdymclus = new TH2D( dxdymclusName.c_str(),
			      (dxdymclusTitle + ";#Delta y (m); #Delta x (m)").c_str(), 
			      hcalybins, 
			      hcalylim.first, 
			      hcalylim.second, 
			      hcalxbins, 
			      hcalxlim.first, 
			      hcalxlim.second);

  // Draw the plot using the created histogram
  tree->Draw(("dx_bc:dy_bc>>" + dxdymclusName).c_str(), dxdymcluscut.c_str(), "COLZ");


  // Draw the mclus plot without spot cut
  hmclus->SetLineColor(kBlack);
  hmclus->GetYaxis()->SetRangeUser(0,1.4*hmclus->GetMaximum());
  hmclus->Draw();

  hmclusanti->SetLineColor(kRed);
  hmclusanti->Draw("same");

  TLegend *leg_mclus = new TLegend(0.11, 0.7, 0.89, 0.89);
  leg_mclus->AddEntry(hmclus, "Validity, Fiducial, and Elastic Cuts", "L");
  leg_mclus->AddEntry(hmclusanti,"Add Antispot Cuts" , "L");
  leg_mclus->Draw();

  // Draw the mclus plot with the spot cut
  cmclus->cd(2);
  hdxdymclus->Draw("colz");

  cmclus->Update();
  cmclus->Write();

  // Create a new canvas
  TCanvas *cW2_p = new TCanvas("cW2_p", "W2_p", 1200, 600);
  cW2_p->cd();

  // Draw W2 plot to look at effects of all cuts
  std::string W2_pName = "hW2_p";
  std::string W2_pcut = tropt_cuts;
  std::string W2_pTitle = "W^{2}";
  TH1D* hW2_p = new TH1D( W2_pName.c_str(), 
			(W2_pTitle + "; GeV^{2}").c_str(), 
			W2bins, 
			W2lim.first, 
			W2lim.second);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("W2>>" + W2_pName).c_str(), W2_pcut.c_str(), "COLZ");

  // Draw dx vs dy plot for comparison with spot cut region
  std::string W2_pspotName = "hW2_pspot";
  std::string W2_pspotcut = p_exp_cuts;
  std::string W2_pspotTitle = "W^{2} (Validity and Proton Spot Cuts)";
  TH1D* hW2_pspot = new TH1D( W2_pspotName.c_str(), 
			    (W2_pspotTitle + "; GeV^{2}").c_str(), 
			    W2bins, 
			    W2lim.first, 
			    W2lim.second);

  // Draw the plot using the created histogram
  tree->Draw(("W2>>" + W2_pspotName).c_str(), W2_pspotcut.c_str(), "COLZ");

  // Draw the W2 plot without spot cut
  hW2_p->SetLineColor(kBlack);
  hW2_p->Draw();

  // Draw the W2 plot with the spot cut
  int transparentRed = TColor::GetColorTransparent(kRed, 0.3);
  hW2_pspot->SetLineColor(kRed);
  hW2_pspot->SetLineWidth(0);
  hW2_pspot->SetFillColor(transparentRed);
  hW2_pspot->SetFillStyle(1001);
  hW2_pspot->Draw("SAME");

  TLegend *leg_pw2 = new TLegend(0.11, 0.7, 0.49, 0.89);
  leg_pw2->AddEntry(hW2_p, "Validity Cuts", "L");
  leg_pw2->AddEntry(hW2_pspot,"Validity and Proton Spot Cuts" , "f");
  leg_pw2->Draw();

  cW2_p->Update();
  cW2_p->Write();

  // Create a new canvas for neutron
  TCanvas *cW2_n = new TCanvas("cW2_n", "W2_n", 1200, 600);
  cW2_n->cd();

  // Draw W2 plot to look at effects of all cuts
  std::string W2_nName = "hW2_n";
  std::string W2_ncut = tropt_cuts;
  std::string W2_nTitle = "W^{2}";
  TH1D* hW2_n = new TH1D( W2_nName.c_str(), 
			(W2_nTitle + "; GeV^{2}").c_str(), 
			W2bins, 
			W2lim.first, 
			W2lim.second);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("W2>>" + W2_nName).c_str(), W2_ncut.c_str(), "COLZ");

  // Draw dx vs dy plot for comparison with spot cut region
  std::string W2_nspotName = "hW2_nspot";
  std::string W2_nspotcut = n_exp_cuts;
  std::string W2_nspotTitle = "W^{2} (Validity and Neutron Spot Cuts)";
  TH1D* hW2_nspot = new TH1D( W2_nspotName.c_str(), 
			    (W2_nspotTitle + "; GeV^{2}").c_str(), 
			    W2bins, 
			    W2lim.first, 
			    W2lim.second);

  // Draw the plot using the created histogram
  tree->Draw(("W2>>" + W2_nspotName).c_str(), W2_nspotcut.c_str(), "COLZ");

  // Draw the W2 plot without spot cut
  hW2_n->SetLineColor(kBlack);
  hW2_n->Draw();

  // Draw the W2 plot with the spot cut
  int transparentBlue = TColor::GetColorTransparent(kBlue, 0.3);
  hW2_nspot->SetLineColor(kBlue);
  hW2_nspot->SetLineWidth(0);
  hW2_nspot->SetFillColor(transparentBlue);
  hW2_nspot->SetFillStyle(1001);
  hW2_nspot->Draw("SAME");

  TLegend *leg_nw2 = new TLegend(0.11, 0.7, 0.49, 0.89);
  leg_nw2->AddEntry(hW2_n, "Validity Cuts", "L");
  leg_nw2->AddEntry(hW2_nspot,"Validity and Neutron Spot Cuts" , "f");
  leg_nw2->Draw();

  cW2_n->Update();
  cW2_n->Write();

  // Draw expected x plot with all cuts, all cuts minus spot, and all cuts antispot (save fidx on all)
  // Create a new canvas for both nucleons
  TCanvas *craw = new TCanvas("craw", "raw", 1200, 600);
  craw->Divide(2,1);
  craw->cd(1);

  std::string xexp_normName = "hxexp_norm";
  std::string xexp_normcut = norm_xcuts;
  std::string xexp_normTitle = "x_{exp}";
  TH1D* hxexp_norm = new TH1D( xexp_normName.c_str(), 
			(xexp_normTitle + "; x_{exp} (m)").c_str(), 
			hcalxbins, 
			hcalxlim.first, 
			hcalxlim.second);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("xexp>>" + xexp_normName).c_str(), xexp_normcut.c_str(), "COLZ");
  hxexp_norm->GetYaxis()->SetRangeUser(0.0,1.5*hxexp_norm->GetMaximum());

  std::string xexp_sigName = "hxexp_sig";
  std::string xexp_sigcut = sig_xcuts;
  std::string xexp_sigTitle = "";
  TH1D* hxexp_sig = new TH1D( xexp_sigName.c_str(), 
			(xexp_sigTitle + "; x_{exp} (m)").c_str(), 
			hcalxbins, 
			hcalxlim.first, 
			hcalxlim.second);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("xexp>>" + xexp_sigName).c_str(), xexp_sigcut.c_str(), "COLZ");

  std::string xexp_antiName = "hxexp_anti";
  std::string xexp_anticut = anti_xcuts;
  std::string xexp_antiTitle = "";
  TH1D* hxexp_anti = new TH1D( xexp_antiName.c_str(), 
			(xexp_antiTitle + "; x_{exp} (m)").c_str(), 
			hcalxbins, 
			hcalxlim.first, 
			hcalxlim.second);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("xexp>>" + xexp_antiName).c_str(), xexp_anticut.c_str(), "COLZ");


  // Draw expected y plot with all cuts, all cuts minus spot, and all cuts antispot (save fidy on all)
  std::string yexp_normName = "hyexp_norm";
  std::string yexp_normcut = norm_ycuts;
  std::string yexp_normTitle = "y_{exp}";
  TH1D* hyexp_norm = new TH1D( yexp_normName.c_str(), 
			(yexp_normTitle + "; y_{exp} (m)").c_str(), 
			hcalybins, 
			hcalylim.first, 
			hcalylim.second);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("yexp>>" + yexp_normName).c_str(), yexp_normcut.c_str(), "COLZ");
  hyexp_norm->GetYaxis()->SetRangeUser(0.0,1.5*hyexp_norm->GetMaximum());

  std::string yexp_sigName = "hyexp_sig";
  std::string yexp_sigcut = sig_ycuts;
  std::string yexp_sigTitle = "";
  TH1D* hyexp_sig = new TH1D( yexp_sigName.c_str(), 
			(yexp_sigTitle + "; y_{exp} (m)").c_str(), 
			hcalybins, 
			hcalylim.first, 
			hcalylim.second);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("yexp>>" + yexp_sigName).c_str(), yexp_sigcut.c_str(), "COLZ");

  std::string yexp_antiName = "hyexp_anti";
  std::string yexp_anticut = anti_ycuts;
  std::string yexp_antiTitle = "";
  TH1D* hyexp_anti = new TH1D( yexp_antiName.c_str(), 
			(yexp_antiTitle + "; y_{exp} (m)").c_str(), 
			hcalybins, 
			hcalylim.first, 
			hcalylim.second);

  // Draw the best cluster plot using the created histogram
  tree->Draw(("yexp>>" + yexp_antiName).c_str(), yexp_anticut.c_str(), "COLZ");

  int transparentGreen = TColor::GetColorTransparent(kGreen, 0.3);
  double hxexpmax = hxexp_norm->GetMaximum();
  double hyexpmax = hyexp_norm->GetMaximum();

  TLine *xexpline3 = new TLine(fidxcut_vec[0], 0.0, fidxcut_vec[0], hxexpmax);
  TLine *xexpline4 = new TLine(fidxcut_vec[1], 0.0, fidxcut_vec[1], hxexpmax);
  xexpline3->SetLineColor(kRed);
  xexpline3->SetLineStyle(2);   
  xexpline3->SetLineWidth(2);
  xexpline4->SetLineColor(kRed);
  xexpline4->SetLineStyle(2);   
  xexpline4->SetLineWidth(2);

  TLine *yexpline3 = new TLine(fidycut_vec[0], 0.0, fidycut_vec[0], hyexpmax);
  TLine *yexpline4 = new TLine(fidycut_vec[1], 0.0, fidycut_vec[1], hyexpmax);
  yexpline3->SetLineColor(kRed);
  yexpline3->SetLineStyle(2);   
  yexpline3->SetLineWidth(2);
  yexpline4->SetLineColor(kRed);
  yexpline4->SetLineStyle(2); 
  yexpline4->SetLineWidth(2); 

  // Draw the xexp comparisons
  hxexp_norm->SetLineColor(kBlack);
  hxexp_norm->Draw();
  hxexp_sig->SetFillStyle(1001);
  hxexp_sig->SetFillColor(transparentGreen);
  hxexp_sig->SetLineColor(kGreen);
  hxexp_sig->SetLineWidth(0);
  hxexp_sig->Draw("same");
  hxexp_anti->SetLineColor(kRed);
  hxexp_anti->SetLineWidth(2);
  hxexp_anti->Draw("same");

  TLegend *lxexp = new TLegend(0.11, 0.75, 0.89, 0.89);
  lxexp->AddEntry(hxexp_norm, "Validity, Fiducial y, Elastic Cuts", "l");
  lxexp->AddEntry(hxexp_sig, "Validity, Fiducial y, Elastic, Nucleon Spot Cuts", "l");
  lxexp->AddEntry(hxexp_anti, "Difference", "l");
  lxexp->Draw();

  craw->Update();

  // Draw the yexp comparisons
  craw->cd(2);
  hyexp_norm->SetLineColor(kBlack);
  hyexp_norm->Draw();
  hyexp_sig->SetFillStyle(1001);
  hyexp_sig->SetFillColor(transparentGreen);
  hyexp_sig->SetLineColor(kGreen);
  hyexp_sig->SetLineWidth(0);
  hyexp_sig->Draw("same");
  hyexp_anti->SetLineColor(kRed);
  hyexp_anti->SetLineWidth(2);
  hyexp_anti->Draw("same");

  TLegend *lyexp = new TLegend(0.11, 0.75, 0.89, 0.89);
  lyexp->AddEntry(hyexp_norm, "Validity, Fiducial x, Elastic Cuts", "l");
  lyexp->AddEntry(hyexp_sig, "Validity, Fiducial x, Elastic, Nucleon Spot Cuts", "l");
  lyexp->AddEntry(hyexp_anti, "Difference", "l");
  lyexp->Draw();

  craw->Update();
  craw->Write();

  // Create canvas and divide it
  TCanvas *cEff = new TCanvas("cEff", "HCal Detection Efficiency", 1200, 600);
  cEff->Divide(2,1);
  cEff->cd(1);

  //Now divide the sig:norm for efficiency per dimension
  std::string yeffName = "hyexp_eff";
  std::string yeffTitle = "HDE nucleon yexp; y_{exp} (m); Nucleon Detection Efficiency";
  TH1D* hyexp_eff = new TH1D(yeffName.c_str(), 
			     yeffTitle.c_str(), 
			     hcalybins, 
			     hcalylim.first, 
			     hcalylim.second);

  // Use Divide() to compute the efficiency
  hyexp_eff->Divide(hyexp_sig, hyexp_norm, 1.0, 1.0, "B"); // "B" option to use binomial errors
  //hyexp_eff->SetMarkerStyle(6);
  //hyexp_eff->SetMarkerColor(kBlack);  
  hyexp_eff->SetLineColor(kBlack);

  TH1D* hyexp_corr = (TH1D*)hyexp_eff->Clone("hyexp_corr"); //Clone for later

  std::string xeffName = "hxexp_eff";
  std::string xeffTitle = "HDE nucleon xexp; x_{exp} (m); Nucleon Detection Efficiency";
  TH1D* hxexp_eff = new TH1D(xeffName.c_str(), 
			     xeffTitle.c_str(), 
			     hcalxbins, 
			     hcalxlim.first, 
			     hcalxlim.second);

  // Use Divide() to compute the efficiency
  hxexp_eff->Divide(hxexp_sig, hxexp_norm, 1.0, 1.0, "B");
  //hxexp_eff->SetMarkerStyle(6);
  //hxexp_eff->SetMarkerColor(kBlack);
  hxexp_eff->SetLineColor(kBlack);

  TH1D* hxexp_corr = (TH1D*)hxexp_eff->Clone("hxexp_corr"); //Clone for later

  // Fit hyexp_eff within specified bounds
  TF1 *fit_hy = new TF1("fit_hy", "pol0", fidycut_vec[0], fidycut_vec[1]);
  hyexp_eff->Fit(fit_hy, "R");  // "R" option for fit range

  // Fit hxexp_eff from the first bin with data or fidxcut_vec[1] to fidxcut_vec[0]
  int first_nonempty_bin = hxexp_eff->FindFirstBinAbove(0);
  double left_binlimit_x = hxexp_eff->GetBinLowEdge(first_nonempty_bin);
  double fidxcutfirst = (left_binlimit_x > fidxcut_vec[1]) ? left_binlimit_x : fidxcut_vec[1];
  TF1 *fit_hx = new TF1("fit_hx", "pol0", fidxcutfirst, fidxcut_vec[0]);

  hxexp_eff->Fit(fit_hx, "R");

  //cout << left_binlimit_x << " " << fidxcut_vec[0] << " " << fidxcutfirst << endl;

  // Drawing hyexp_eff with fit and legend
  hxexp_eff->GetYaxis()->SetRangeUser(0.0,1.0);
  hxexp_eff->Draw("L");
  fit_hx->SetLineColor(kRed);

  TLegend *leg_x = new TLegend(0.22, 0.12, 0.7, 0.3);
  double xMeanEff = fit_hx->GetParameter(0);
  double xErrEff = fit_hx->GetParError(0);
  leg_x->AddEntry(hxexp_eff, "x_{exp} Detection Efficiency", "P");
  leg_x->AddEntry(fit_hx, Form("Mean Eff = %.3f #pm %.3f", xMeanEff, xErrEff), "L");
  leg_x->Draw();

  // Drawing hxexp_eff with fit and legend
  cEff->cd(2);
  hyexp_eff->GetYaxis()->SetRangeUser(0.0,1.0);
  hyexp_eff->Draw("L");
  fit_hy->SetLineColor(kRed);

  TLegend *leg_y = new TLegend(0.22, 0.12, 0.7, 0.3);
  double yMeanEff = fit_hy->GetParameter(0);
  double yErrEff = fit_hy->GetParError(0);
  leg_y->AddEntry(hyexp_eff, "y_{exp} Detection Efficiency", "P");
  leg_y->AddEntry(fit_hy, Form("Mean Eff = %.3f #pm %.3f", yMeanEff, yErrEff), "L");
  leg_y->Draw();

  cEff->Update();
  cEff->Write();

  /////////////////////////
  /////////////////////////
  /////////////////////////
  /////////////////////////
  //Ratio
  /////////////////////////
  /////////////////////////
  /////////////////////////
  /////////////////////////

  // Create a new canvas for both nucleons
  TCanvas *cratio = new TCanvas("cratio", "ratio", 2100, 1800);
  cratio->Divide(2,2);
  cratio->cd(1);

  //Get xexp total after harm and earm cuts for proton
  std::string xexp_ratpName = "hxexp_ratp";
  std::string xexp_ratpcut = sig_pcuts;
  std::string xexp_ratpTitle = "x_{exp}";

  double ratxllim = fidxcut_vec[1]*1.1;
  double ratxulim = fidxcut_vec[0]*1.1;

  TH1D* hxexp_ratp = new TH1D( xexp_ratpName.c_str(), 
			      (xexp_ratpTitle + "; x_{exp} (m)").c_str(), 
			      hcalybins, 
			      ratxllim, 
			      ratxulim);

  tree->Draw(("xexp>>" + xexp_ratpName).c_str(), xexp_ratpcut.c_str(), "COLZ");
  hxexp_ratp->GetYaxis()->SetRangeUser(0.0,1.0);

  //Get xexp total after harm and earm cuts for neutron
  std::string xexp_ratnName = "hxexp_ratn";
  std::string xexp_ratncut = sig_ncuts;
  std::string xexp_ratnTitle = "x_{exp}";
  TH1D* hxexp_ratn = new TH1D( xexp_ratnName.c_str(), 
			      (xexp_ratnTitle + "; x_{exp} (m)").c_str(), 
			      hcalybins, 
			      ratxllim, 
			      ratxulim);

  tree->Draw(("xexp>>" + xexp_ratnName).c_str(), xexp_ratncut.c_str(), "COLZ");
  hxexp_ratn->GetYaxis()->SetRangeUser(0.0,1.0);

  //Get n:p ratio for xexp
  std::string xratName = "hxexp_rat";
  std::string xratTitle = "n:p ratio xexp; x_{exp} (m); n:p";
  TH1D* hratx = new TH1D(xratName.c_str(), 
			 xratTitle.c_str(), 
			 hcalybins, 
			 ratxllim, 
			 ratxulim);

  hratx->Divide(hxexp_ratn, hxexp_ratp, 1.0, 1.0, "B"); // "B" option to use binomial errors

  // Handle points where the ratio is unity with no error
  for (int i = 1; i <= hratx->GetNbinsX(); ++i) {
    if (hratx->GetBinContent(i) == 1.0) {
      hratx->SetBinContent(i, 0); // Set to 0
      hratx->SetBinError(i, 0);   // Set error to 0
    }
  }


  TF1 *fit_hratx = new TF1("fit_hratx", "pol0", fidxcut_vec[0], fidxcut_vec[1]);
  hratx->GetYaxis()->SetRangeUser(0.0,1.0);
  hratx->Fit(fit_hratx, "R");  // "R" option for fit range

  double npratio_xexp = fit_hratx->GetParameter(0);
  double npratioerr_xexp = fit_hratx->GetParError(0);

  hratx->SetLineColor(kBlack);
  hratx->Draw("L");

  cratio->Update();
  cratio->cd(2);

  //Get yexp total after harm and earm cuts for proton
  std::string yexp_ratpName = "hyexp_ratp";
  std::string yexp_ratpcut = sig_pcuts;
  std::string yexp_ratpTitle = "x_{exp}";

  double ratyllim = fidycut_vec[1]*1.1;
  double ratyulim = fidycut_vec[0]*1.1;

  TH1D* hyexp_ratp = new TH1D( yexp_ratpName.c_str(), 
			      (yexp_ratpTitle + "; x_{exp} (m)").c_str(), 
			      hcalybins, 
			      ratyllim, 
			      ratyulim);

  tree->Draw(("yexp>>" + yexp_ratpName).c_str(), yexp_ratpcut.c_str(), "COLZ");
  hyexp_ratp->GetYaxis()->SetRangeUser(0.0,1.0);

  //Get yexp total after harm and earm cuts for neutron
  std::string yexp_ratnName = "hyexp_ratn";
  std::string yexp_ratncut = sig_ncuts;
  std::string yexp_ratnTitle = "x_{exp}";
  TH1D* hyexp_ratn = new TH1D( yexp_ratnName.c_str(), 
			      (yexp_ratnTitle + "; x_{exp} (m)").c_str(), 
			      hcalybins, 
			      ratyllim, 
			      ratyulim);

  tree->Draw(("yexp>>" + yexp_ratnName).c_str(), yexp_ratncut.c_str(), "COLZ");
  hyexp_ratn->GetYaxis()->SetRangeUser(0.0,1.0);

  //Get n:p ratio for yexp
  std::string yratName = "hyexp_rat";
  std::string yratTitle = "n:p ratio yexp; y_{exp} (m); n:p";
  TH1D* hraty = new TH1D(yratName.c_str(), 
			 yratTitle.c_str(), 
			 hcalybins, 
			 ratyllim, 
			 ratyulim);

  hraty->Divide(hyexp_ratn, hyexp_ratp, 1.0, 1.0, "B"); // "B" option to use binomial errors

  // Handle points where the ratio is unity with no error
  for (int i = 1; i <= hraty->GetNbinsX(); ++i) {
    if (hraty->GetBinContent(i) == 1.0) {
      hraty->SetBinContent(i, 0); // Set to 0
      hraty->SetBinError(i, 0);   // Set error to 0
    }
  }

  TF1 *fit_hraty = new TF1("fit_hraty", "pol0", fidycut_vec[0], fidycut_vec[1]);
  hraty->GetYaxis()->SetRangeUser(0.0,1.0);
  hraty->Fit(fit_hraty, "R");  // "R" option for fit range

  double npratio_yexp = fit_hraty->GetParameter(0);
  double npratioerr_yexp = fit_hraty->GetParError(0);

  hraty->SetLineColor(kBlack);
  hraty->Draw("L");

  cratio->Update();
  cratio->cd(3);

  // Create a 2D histogram with appropriate binning
  std::string vsnratName = "hvsnexp_rat";
  std::string vsnratcut = sig_ncuts;
  std::string vsnratTitle = "n:p ratio xexp vs yexp; y_{exp} (m); x_{exp} (m)";
  TH2D *hratvsn = new TH2D(vsnratName.c_str(), 
			  vsnratTitle.c_str(), 
			  60, 
			  hcalylim.first, 
			  hcalylim.second, 
			  120, 
			  hcalxlim.first, 
			  hcalxlim.second);

  tree->Draw(("xexp:yexp>>" + vsnratName).c_str(), vsnratcut.c_str(), "COLZ");

  std::string vspratName = "hvspexp_rat";
  std::string vspratcut = sig_pcuts;
  std::string vspratTitle = "n:p ratio xexp vs yexp; y_{exp} (m); x_{exp} (m)";
  TH2D *hratvsp = new TH2D(vspratName.c_str(), 
			  vspratTitle.c_str(), 
			  60, 
			  hcalylim.first, 
			  hcalylim.second, 
			  120, 
			  hcalxlim.first, 
			  hcalxlim.second);

  tree->Draw(("xexp:yexp>>" + vspratName).c_str(), vspratcut.c_str(), "COLZ");

  //ratio
  std::string vsratName = "hvsexp_rat";
  std::string vsratTitle = "n:p ratio xexp vs yexp; y_{exp} (m); x_{exp} (m)";
  TH2D *hratvs = new TH2D(vsratName.c_str(), 
			  vsratTitle.c_str(), 
			  60, 
			  hcalylim.first, 
			  hcalylim.second, 
			  120, 
			  hcalxlim.first, 
			  hcalxlim.second);

  hratvs->Divide(hratvsn, hratvsp, 1.0, 1.0, "B"); // "B" option to use binomial errors

  hratvs->Draw("colz");

  cratio->Update();
  cratio->cd(4);

  // Create a ratio vs W2 histogram
  std::string vsnw2ratName = "hvsnexp_w2rat";
  std::string vsnw2ratcut = sig_ncuts_noW2;
  std::string vsnw2ratTitle = "neutron W^{2}; W^{2} (GeV^{2})";
  TH1D *hw2ratvsn = new TH1D(vsnw2ratName.c_str(), 
			  vsnw2ratTitle.c_str(), 
			  W2bins, 
			  W2lim.first, 
			  W2lim.second);

  tree->Draw(("W2>>" + vsnw2ratName).c_str(), vsnw2ratcut.c_str(), "COLZ");

  std::string vspw2ratName = "hvspexp_w2rat";
  std::string vspw2ratcut = sig_pcuts_noW2;
  std::string vspw2ratTitle = "proton W^{2}; W^{2} (GeV^{2})";
  TH1D *hw2ratvsp = new TH1D(vspw2ratName.c_str(), 
			  vspw2ratTitle.c_str(), 
			  W2bins, 
			  W2lim.first, 
			  W2lim.second);

  tree->Draw(("W2>>" + vspw2ratName).c_str(), vspw2ratcut.c_str(), "COLZ");

  //w2ratio
  std::string vsw2ratName = "hvsexp_w2rat";
  std::string vsw2ratTitle = "n:p ratio W^{2}; W^{2} (GeV^{2}; n:p ratio)";
  TH1D *hw2ratvs = new TH1D(vsw2ratName.c_str(), 
			  vsw2ratTitle.c_str(), 
			  W2bins, 
			  W2lim.first, 
			  W2lim.second);

  hw2ratvs->Divide(hw2ratvsn, hw2ratvsp, 1.0, 1.0, "B"); // "B" option to use binomial errors

  cout << "Checking for source of ratio vs W2 unity values without errors.." << endl;

  for (int i = 1; i <= hw2ratvs->GetNbinsX(); ++i) {
    if (hw2ratvs->GetBinError(i) == 0) {
      hw2ratvs->SetBinContent(i, 0);
      hw2ratvs->SetBinError(i, 0);
    }
  }

  TF1 *fit_ratvs = new TF1("fit_ratvs", "pol0", W2fit.first, W2fit.second);
  hw2ratvs->Fit(fit_ratvs, "R");

  hw2ratvs->SetLineColor(kBlack);
  hw2ratvs->GetYaxis()->SetRangeUser(0.0,1.0);
  hw2ratvs->Draw("L");
  fit_ratvs->SetLineColor(kRed);

  cratio->Update();
  cratio->Write();


  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  //dy Stability canvas
  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  //////////////////////////////


  // Create a new canvas for both nucleons
  TCanvas *cstabdy = new TCanvas("cstabdy", "stabdy", 1200, 800);
  //cstabdy->Divide(2,2);
  cstabdy->cd();

 // Create a ratio vs dy histogram
  std::string vsndyratName = "hvsnexp_dyrat";
  std::string vsndyratcut = sig_ncuts_nody;
  std::string vsndyratTitle = "neutron E/p; E/p";
  TH1D *hdyratvsn = new TH1D(vsndyratName.c_str(), 
			  vsndyratTitle.c_str(), 
			  dybins, 
			  dylim.first, 
			  dylim.second);

  tree->Draw(("dy_bc>>" + vsndyratName).c_str(), vsndyratcut.c_str(), "COLZ");

  std::string vspdyratName = "hvspexp_dyrat";
  std::string vspdyratcut = sig_pcuts_nody;
  std::string vspdyratTitle = "proton E/p; E/p";
  TH1D *hdyratvsp = new TH1D(vspdyratName.c_str(), 
			  vspdyratTitle.c_str(), 
			  dybins, 
			  dylim.first, 
			  dylim.second);

  tree->Draw(("dy_bc>>" + vspdyratName).c_str(), vspdyratcut.c_str(), "COLZ");

  //dyratio
  std::string vsdyratName = "hvsexp_dyrat";
  std::string vsdyratTitle = "n:p ratio dy; dy; n:p ratio)";
  TH1D *hdyratvs = new TH1D(vsdyratName.c_str(), 
			  vsdyratTitle.c_str(), 
			  dybins, 
			  dylim.first, 
			  dylim.second);

  hdyratvs->Divide(hdyratvsn, hdyratvsp, 1.0, 1.0, "B"); // "B" option to use binomial errors

  cout << "Checking for source of ratio vs dy unity values without errors.." << endl;

  for (int i = 1; i <= hdyratvs->GetNbinsX(); ++i) {
    if (hdyratvs->GetBinError(i) == 0) {
      hdyratvs->SetBinContent(i, 0);
      hdyratvs->SetBinError(i, 0);
    }
  }

  TF1 *fit_ratvs_dy = new TF1("fit_ratvs_dy", "pol0", dyfit.first, dyfit.second);
  hdyratvs->Fit(fit_ratvs_dy, "R");

  hdyratvs->SetLineColor(kBlack);
  hdyratvs->GetYaxis()->SetRangeUser(0.0,1.0);
  hdyratvs->Draw("L");
  fit_ratvs_dy->SetLineColor(kRed);

  cstabdy->Update();
  cstabdy->Write();

  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  //Stability canvas
  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  //////////////////////////////
  //////////////////////////////

  // Create a new canvas for both nucleons
  TCanvas *cstab = new TCanvas("cstab", "stab", 2100, 1800);
  cstab->Divide(2,4);
  cstab->cd(1);

  // Create a ratio vs Ep histogram
  std::string vsnEpratName = "hvsnexp_Eprat";
  std::string vsnEpratcut = sig_ncuts_noEp;
  std::string vsnEpratTitle = "neutron E/p; E/p";
  TH1D *hEpratvsn = new TH1D(vsnEpratName.c_str(), 
			  vsnEpratTitle.c_str(), 
			  Epbins, 
			  Eplim.first, 
			  Eplim.second);

  tree->Draw(("bb_etot_over_p>>" + vsnEpratName).c_str(), vsnEpratcut.c_str(), "COLZ");

  std::string vspEpratName = "hvspexp_Eprat";
  std::string vspEpratcut = sig_pcuts_noEp;
  std::string vspEpratTitle = "proton E/p; E/p";
  TH1D *hEpratvsp = new TH1D(vspEpratName.c_str(), 
			  vspEpratTitle.c_str(), 
			  Epbins, 
			  Eplim.first, 
			  Eplim.second);

  tree->Draw(("bb_etot_over_p>>" + vspEpratName).c_str(), vspEpratcut.c_str(), "COLZ");

  //Epratio
  std::string vsEpratName = "hvsexp_Eprat";
  std::string vsEpratTitle = "n:p ratio E/p; E/p; n:p ratio";
  TH1D *hEpratvs = new TH1D(vsEpratName.c_str(), 
			  vsEpratTitle.c_str(), 
			  Epbins, 
			  Eplim.first, 
			  Eplim.second);

  hEpratvs->Divide(hEpratvsn, hEpratvsp, 1.0, 1.0, "B"); // "B" option to use binomial errors

  cout << "Checking for source of ratio vs E/p unity values without errors.." << endl;

  for (int i = 1; i <= hEpratvs->GetNbinsX(); ++i) {
    if (hEpratvs->GetBinError(i) == 0) {
      hEpratvs->SetBinContent(i, 0);
      hEpratvs->SetBinError(i, 0);
    }
  }

  TF1 *fit_ratvs_Ep = new TF1("fit_ratvs_Ep", "pol0", Epfit.first, Epfit.second);
  hEpratvs->Fit(fit_ratvs_Ep, "R");

  hEpratvs->SetLineColor(kBlack);
  hEpratvs->GetYaxis()->SetRangeUser(0.0,1.0);
  hEpratvs->GetXaxis()->SetLabelSize(0.05);
  hEpratvs->GetYaxis()->SetLabelSize(0.05);
  hEpratvs->GetXaxis()->SetTitleSize(0.05);
  hEpratvs->GetYaxis()->SetTitleSize(0.05);
  hEpratvs->Draw("L");
  fit_ratvs_Ep->SetLineColor(kRed);
  fit_ratvs_Ep->SetLineWidth(2);
  
  cstab->Update();
  cstab->cd(2);

  // Create a ratio vs psE histogram
  std::string vsnpsEratName = "hvsnexp_psErat";
  std::string vsnpsEratcut = sig_ncuts_nopsE;
  std::string vsnpsEratTitle = "neutron PS E; GeV";
  TH1D *hpsEratvsn = new TH1D(vsnpsEratName.c_str(), 
			  vsnpsEratTitle.c_str(), 
			  psEbins, 
			  psElim.first, 
			  psElim.second);

  tree->Draw(("bb_ps_e>>" + vsnpsEratName).c_str(), vsnpsEratcut.c_str(), "COLZ");

  std::string vsppsEratName = "hvspexp_psErat";
  std::string vsppsEratcut = sig_pcuts_nopsE;
  std::string vsppsEratTitle = "proton PS E; GeV";
  TH1D *hpsEratvsp = new TH1D(vsppsEratName.c_str(), 
			  vsppsEratTitle.c_str(), 
			  psEbins, 
			  psElim.first, 
			  psElim.second);

  tree->Draw(("bb_ps_e>>" + vsppsEratName).c_str(), vsppsEratcut.c_str(), "COLZ");

  //psEratio
  std::string vspsEratName = "hvsexp_psErat";
  std::string vspsEratTitle = "n:p ratio PS E; GeV; n:p ratio";
  TH1D *hpsEratvs = new TH1D(vspsEratName.c_str(), 
			  vspsEratTitle.c_str(), 
			  psEbins, 
			  psElim.first, 
			  psElim.second);

  hpsEratvs->Divide(hpsEratvsn, hpsEratvsp, 1.0, 1.0, "B"); // "B" option to use binomial errors

  cout << "Checking for source of ratio vs PS E unity values without errors.." << endl;

  for (int i = 1; i <= hpsEratvs->GetNbinsX(); ++i) {
    if (hpsEratvs->GetBinError(i) == 0) {
      hpsEratvs->SetBinContent(i, 0);
      hpsEratvs->SetBinError(i, 0);
    }
  }

  TF1 *fit_ratvs_psE = new TF1("fit_ratvs_psE", "pol0", psEfit.first, psEfit.second);
  hpsEratvs->Fit(fit_ratvs_psE, "R");

  hpsEratvs->SetLineColor(kBlack);
  hpsEratvs->GetYaxis()->SetRangeUser(0.0,1.0);
  hpsEratvs->Draw("L");
  fit_ratvs_psE->SetLineColor(kRed);

  cstab->Update();
  cstab->cd(3);

  // Create a ratio vs HCal E histogram
  std::string vsnEratName = "hvsnexp_Erat";
  std::string vsnEratcut = sig_ncuts_noE;
  std::string vsnEratTitle = "neutron HCal E; GeV";
  TH1D *hEratvsn = new TH1D(vsnEratName.c_str(), 
			  vsnEratTitle.c_str(), 
			  Ebins, 
			  Elim.first, 
			  Elim.second);

  tree->Draw(("hcale_bc>>" + vsnEratName).c_str(), vsnEratcut.c_str(), "COLZ");

  std::string vspEratName = "hvspexp_Erat";
  std::string vspEratcut = sig_pcuts_noE;
  std::string vspEratTitle = "proton HCal E; GeV";
  TH1D *hEratvsp = new TH1D(vspEratName.c_str(), 
			  vspEratTitle.c_str(), 
			  Ebins, 
			  Elim.first, 
			  Elim.second);

  tree->Draw(("hcale_bc>>" + vspEratName).c_str(), vspEratcut.c_str(), "COLZ");

  //HCal E ratio
  std::string vsEratName = "hvsexp_Erat";
  std::string vsEratTitle = "n:p ratio HCal E; GeV; n:p ratio";
  TH1D *hEratvs = new TH1D(vsEratName.c_str(), 
			  vsEratTitle.c_str(), 
			  Ebins, 
			  Elim.first, 
			  Elim.second);

  hEratvs->Divide(hEratvsn, hEratvsp, 1.0, 1.0, "B"); // "B" option to use binomial errors

  cout << "Checking for source of ratio vs HCal E unity values without errors.." << endl;

  for (int i = 1; i <= hEratvs->GetNbinsX(); ++i) {
    if (hEratvs->GetBinError(i) == 0) {
      hEratvs->SetBinContent(i, 0);
      hEratvs->SetBinError(i, 0);
    }
  }

  TF1 *fit_ratvs_E = new TF1("fit_ratvs_E", "pol0", Efit.first, Efit.second);
  hEratvs->Fit(fit_ratvs_E, "R");

  hEratvs->SetLineColor(kBlack);
  hEratvs->GetYaxis()->SetRangeUser(0.0,1.0);
  hEratvs->Draw("L");
  fit_ratvs_E->SetLineColor(kRed);

  cstab->Update();
  cstab->cd(4);

  // Create a ratio vs coin histogram
  std::string vsncoinratName = "hvsnexp_coinrat";
  std::string vsncoinratcut = sig_ncuts_nocoin;
  std::string vsncoinratTitle = "neutron coin time; ns; n:p ratio";
  TH1D *hcoinratvsn = new TH1D(vsncoinratName.c_str(), 
			  vsncoinratTitle.c_str(), 
			  coinbins, 
			  coinlim.first, 
			  coinlim.second);

  tree->Draw(("coin>>" + vsncoinratName).c_str(), vsncoinratcut.c_str(), "COLZ");

  std::string vspcoinratName = "hvspexp_coinrat";
  std::string vspcoinratcut = sig_pcuts_nocoin;
  std::string vspcoinratTitle = "proton coin time; ns; n:p ratio";
  TH1D *hcoinratvsp = new TH1D(vspcoinratName.c_str(), 
			  vspcoinratTitle.c_str(), 
			  coinbins, 
			  coinlim.first, 
			  coinlim.second);

  tree->Draw(("coin>>" + vspcoinratName).c_str(), vspcoinratcut.c_str(), "COLZ");

  //coinratio
  std::string vscoinratName = "hvsexp_coinrat";
  std::string vscoinratTitle = "n:p ratio coin time; ns; n:p ratio";
  TH1D *hcoinratvs = new TH1D(vscoinratName.c_str(), 
			  vscoinratTitle.c_str(), 
			  coinbins, 
			  coinlim.first, 
			  coinlim.second);

  hcoinratvs->Divide(hcoinratvsn, hcoinratvsp, 1.0, 1.0, "B"); // "B" option to use binomial errors

  cout << "Checking for source of ratio vs coin unity values without errors.." << endl;

  for (int i = 1; i <= hcoinratvs->GetNbinsX(); ++i) {
    if (hcoinratvs->GetBinError(i) == 0) {
      hcoinratvs->SetBinContent(i, 0);
      hcoinratvs->SetBinError(i, 0);
    }
  }

  TF1 *fit_ratvs_coin = new TF1("fit_ratvs_coin", "pol0", coinfit.first, coinfit.second);
  hcoinratvs->Fit(fit_ratvs_coin, "R");

  hcoinratvs->SetLineColor(kBlack);
  hcoinratvs->GetYaxis()->SetRangeUser(0.0,3.0);
  hcoinratvs->Draw("L");
  fit_ratvs_coin->SetLineColor(kRed);

  cstab->Update();
  cstab->cd(5);

  hw2ratvs->Draw("L");

  cstab->Update();
  cstab->cd(6);

  hdyratvs->Draw("L");
  cstab->cd(7);

  // Create a ratio vs vz histogram
  std::string vsnvzratName = "hvsnexp_vzrat";
  std::string vsnvzratcut = sig_ncuts_novz;
  std::string vsnvzratTitle = "neutron vz; ns; n:p ratio";
  TH1D *hvzratvsn = new TH1D(vsnvzratName.c_str(), 
			  vsnvzratTitle.c_str(), 
			  vzbins, 
			  vzlim.first, 
			  vzlim.second);

  tree->Draw(("bb_tr_vz>>" + vsnvzratName).c_str(), vsnvzratcut.c_str(), "COLZ");

  std::string vspvzratName = "hvspexp_vzrat";
  std::string vspvzratcut = sig_pcuts_novz;
  std::string vspvzratTitle = "proton vz; ns; n:p ratio";
  TH1D *hvzratvsp = new TH1D(vspvzratName.c_str(), 
			  vspvzratTitle.c_str(), 
			  vzbins, 
			  vzlim.first, 
			  vzlim.second);

  tree->Draw(("bb_tr_vz>>" + vspvzratName).c_str(), vspvzratcut.c_str(), "COLZ");

  //vzratio
  std::string vsvzratName = "hvsexp_vzrat";
  std::string vsvzratTitle = "n:p ratio vz; ns; n:p ratio";
  TH1D *hvzratvs = new TH1D(vsvzratName.c_str(), 
			  vsvzratTitle.c_str(), 
			  vzbins, 
			  vzlim.first, 
			  vzlim.second);

  hvzratvs->Divide(hvzratvsn, hvzratvsp, 1.0, 1.0, "B"); // "B" option to use binomial errors

  cout << "Checking for source of ratio vs vz unity values without errors.." << endl;

  for (int i = 1; i <= hvzratvs->GetNbinsX(); ++i) {
    if (hvzratvs->GetBinError(i) == 0) {
      hvzratvs->SetBinContent(i, 0);
      hvzratvs->SetBinError(i, 0);
    }
  }

  TF1 *fit_ratvs_vz = new TF1("fit_ratvs_vz", "pol0", vzfit.first, vzfit.second);
  hvzratvs->Fit(fit_ratvs_vz, "R");

  hvzratvs->SetLineColor(kBlack);
  hvzratvs->GetYaxis()->SetRangeUser(0.0,1.0);
  hvzratvs->Draw("L");
  fit_ratvs_vz->SetLineColor(kRed);

  cstab->Update();
  cstab->cd(8);

  // Create a ratio vs chi2 histogram
  std::string vsnchi2ratName = "hvsnexp_chi2rat";
  std::string vsnchi2ratcut = sig_ncuts_nochi2;
  std::string vsnchi2ratTitle = "neutron chi2; ns; n:p ratio";
  TH1D *hchi2ratvsn = new TH1D(vsnchi2ratName.c_str(), 
			  vsnchi2ratTitle.c_str(), 
			  chi2bins, 
			  chi2lim.first, 
			  chi2lim.second);

  tree->Draw(("bb_tr_chi2>>" + vsnchi2ratName).c_str(), vsnchi2ratcut.c_str(), "COLZ");

  std::string vspchi2ratName = "hvspexp_chi2rat";
  std::string vspchi2ratcut = sig_pcuts_nochi2;
  std::string vspchi2ratTitle = "proton chi2; ns; n:p ratio";
  TH1D *hchi2ratvsp = new TH1D(vspchi2ratName.c_str(), 
			  vspchi2ratTitle.c_str(), 
			  chi2bins, 
			  chi2lim.first, 
			  chi2lim.second);

  tree->Draw(("bb_tr_chi2>>" + vspchi2ratName).c_str(), vspchi2ratcut.c_str(), "COLZ");

  //chi2ratio
  std::string vschi2ratName = "hvsexp_chi2rat";
  std::string vschi2ratTitle = "n:p ratio chi2; ns; n:p ratio";
  TH1D *hchi2ratvs = new TH1D(vschi2ratName.c_str(), 
			  vschi2ratTitle.c_str(), 
			  chi2bins, 
			  chi2lim.first, 
			  chi2lim.second);

  hchi2ratvs->Divide(hchi2ratvsn, hchi2ratvsp, 1.0, 1.0, "B"); // "B" option to use binomial errors

  cout << "Checking for source of ratio vs chi2 unity values without errors.." << endl;

  for (int i = 1; i <= hchi2ratvs->GetNbinsX(); ++i) {
    if (hchi2ratvs->GetBinError(i) == 0) {
      hchi2ratvs->SetBinContent(i, 0);
      hchi2ratvs->SetBinError(i, 0);
    }
  }

  TF1 *fit_ratvs_chi2 = new TF1("fit_ratvs_chi2", "pol0", chi2fit.first, chi2fit.second);
  hchi2ratvs->Fit(fit_ratvs_chi2, "R");

  hchi2ratvs->SetLineColor(kBlack);
  hchi2ratvs->GetYaxis()->SetRangeUser(0.0,1.0);
  hchi2ratvs->Draw("L");
  fit_ratvs_chi2->SetLineColor(kRed);

  cstab->Update();
  cstab->Write();

  /////////////////////
  /////////////////////
  /////////////////////
  /////////////////////
  //reporting
  /////////////////////
  /////////////////////
  /////////////////////

  // Calculate efficiency over both directions and estimate error
  double teff = (xMeanEff/pow(xErrEff,2)+yMeanEff/pow(yErrEff,2)) / (1/pow(xErrEff,2)+1/pow(yErrEff,2));
  double terr = sqrt(1 / (1/pow(xErrEff,2)+1/pow(yErrEff,2)));

  double trat = (npratio_xexp/pow(npratioerr_xexp,2)+npratio_yexp/pow(npratioerr_yexp,2)) / (1/pow(npratioerr_xexp,2)+1/pow(npratioerr_yexp,2));
  double traterr = sqrt(1 / (1/pow(npratioerr_xexp,2)+1/pow(npratioerr_yexp,2)));

  // Define all cuts and their corresponding names
  std::vector<std::pair<std::string, std::string>> cutsWithNames = {
    {"Track Validity Cuts", tra_cuts},
    // {"Track and Optical x Validity Cuts", W2optycut}, // Uncomment if needed
    {"Optical Validity Cuts", opt_cuts},
    {"Elastic Electron Arm Cuts", earmcut},
    {"Elastic Hadron Arm Cuts", harmcut},
    {"Proton Spot Cuts", p_spotcut},
    {"Neutron Spot Cuts", n_spotcut},
    {"Fiducial Cuts", fid_cuts}
  };
  
  // Concatenate all cuts with their names
  std::string allCuts = concatenateCutsWithNames(cutsWithNames);

  std::string efficiencypart = "&&_____&&Efficiency&&xexp: " + 
    std::to_string(xMeanEff) + " +/- " + 
    std::to_string(xErrEff) + "&&yexp: " + 
    std::to_string(yMeanEff) + " +/- " + 
    std::to_string(yErrEff) + "&&combined: " +
    std::to_string(teff) + " +/- " + 
    std::to_string(terr) + "&&";

  allCuts += efficiencypart;

  std::string npratiopart = "&&_____&&neutron:proton Ratio&&xexp: " + 
    std::to_string(npratio_xexp) + " +/- " + 
    std::to_string(npratioerr_xexp) + "&&yexp: " + 
    std::to_string(npratio_yexp) + " +/- " + 
    std::to_string(npratioerr_yexp) + "&&combined: " +
    std::to_string(trat) + " +/- " + 
    std::to_string(traterr) + "&&";

  allCuts += npratiopart;

  // Display all concatenated cuts on a single canvas
  std::string disnameall = "All Cuts";
  TCanvas *cdisall = new TCanvas("cdis", disnameall.c_str(), 1000, 800);
  util::parseAndDisplayCuts(disnameall.c_str(), allCuts.c_str(), cdisall);

  cdisall->Write();
  
  cout << "All plots created. Output file located here: " << fout_path << endl;

}
