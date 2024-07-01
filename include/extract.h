#ifndef EXTRACT_H
#define EXTRACT_H

#include <iostream>
#include <cmath>
#include <vector>
#include <TCanvas.h>
#include <TLatex.h>
#include <TString.h>
#include <sstream>
#include <chrono>
#include <iomanip>
#include <TPaveText.h>
#include <vector>
#include <string>

namespace extract {

  //Gets {{GMn, dGMn}, {GMn/GD/mun, dGMn/GD/mun}, {GD, abs(mun)}, {NRCS, dNRCS}} where NCS is the neutron reduced cross section
  std::vector<std::pair<double,double>> extract_GMn_from_simc(double Rsf,       //Ratio scale factor
							      double Rsf_err,   //Total error on scale factor
							      double Q2,        //Q2 for this kinematic
							      double E_beam,    //Beamn energy for this kine
							      double BB_angle,  //Big Bite angle in radians
							      bool tpecorr);    //bool to switch off TPE corrections from model-dependent extraction

  //Same parameters as above with error components separated
  std::vector<std::tuple<double, double, double, double, double>> extract_GMn_from_simc_full(double Rsf, 
											     double Rsf_err_stat, 
											     double Rsf_err_syst, 
											     double Q2, 
											     double E_beam, 
											     double BB_angle, 
											     bool tpecorr = true);


}

#endif
