// S.Seeds. 3.14.24 This script will extract GMn for SBS experiment 12-09-019 using the following:
//  Riordan parameterization for GEn in simc
//  Kelly parameterization for GEp, GMn, GMp in simc
//  Arrington07 for GMp, GEp for proton model dependent extraction (TPE corrected)
//  Kelly parameterization for GEn (for now, pending TPE corrected alternative)
//
// Parameterizations absorb mu_n and mu_p such that magnetic FF output -> GMp(n)/mu_p(n)
// All simc parameterization from https://github.com/MarkKJones/simc_gfortran, branch bigbite, commit de9c643 accessed 3.18.24
// Extraction parameters from Ye, Arrington, Hill, & Lee Parameterization, 2018 (https://doi.org/10.1016/j.physletb.2017.11.023) and
// Arrington07 - Global Analysis of proton elastic form factor data with tpe corrections (https://journals.aps.org/prc/pdf/10.1103/PhysRevC.76.035205)
// As of 3.21.24 this script assumes negligable error in epsilon, Q2, tau, and all physical constants

#include <iostream>
#include <cmath>
#include <vector>
#include "../../../include/gmn.h"

using namespace std;

// Constants
const double M_p = 0.938272; // Proton mass in GeV
const double M_n = 0.939565; // Neutron mass in GeV
const double Pi = acos(-1);
const double alpha = 1.0 / 137.0; // Fine-structure constant
const double GD_const = 0.71; // Dipole form factor range factor
const double mu_p = 2.79284735; //Proton magnetic moment
const double mu_n = -1.9130427; //Neutron anomalous magnetic moment
const double M_pi = 0.13957; // Charged pion mass in GeV
const double ye_tnot = -0.7; // Ye parameterization constant t0

// Interpolated error from Arrington 07 in the form {GEp_err, GMp_err}
const std::pair<double, double> arrerr_sbs4 = {0.04797657, 0.018205331};
const std::pair<double, double> arrerr_sbs8 = {0.066452257, 0.02521617};
const std::pair<double, double> arrerr_sbs9 = {0.066452257, 0.02521617};

const double GE_over_GD_errslope = 0.0184; //via interpolated fit to Arrington07 list
const double GM_over_GDmup_errslope = 0.0025; //via interpolated fit to Arrington07 list

// Forward declarations
double calculateQ2(double E_beam, double BB_angle, double M);
double calculateYeZ(double Q2, double tcut, double tnot);
double calculateYeX(double Q2);
double calcRiordan(double tau, double a1, double a2, double a3, double b1, double b2, double b3);
double calcYe(double yez);
double calcYe_err(double yex, double GD);
double getRCSsimc(double tau_value, double epsilon_value, double GD, bool isproton);
std::pair<double,double> getRCSmdp_with_error(double tau_p, double epsilon_p, double GD);
std::pair<double, double> getRCSmdp_borntpe_with_error(double tau_p, double epsilon_p, double GD);
std::pair<double, double> getGMn_and_error_from_nRCS(double nRCS, double tau_value, double epsilon_value, double GD, double err_nRCS, double GEn, double err_GEn);

// Main
void extract_GMn_from_simc(int kine=4, double Rsf=0.953, double Rsf_err= 0.0951, double Q2=3.0335, bool tpecorr = true) {

  cout << "Extracting GMn from SIMC reduced cross section and global fits to GEp, GEn, and GMn..." << endl << endl;

  // Set up configuration and tune objects to load analysis parameters
  SBSconfig config(kine,0);

  // Obtain configuration pars from config class
  double E_beam = config.GetEbeam();
  double BB_angle = config.GetBBtheta_rad(); // In radians

  // Calculate average Q2
  double avg_mass = (M_p+M_n)/2;

  double Q2_central = calculateQ2(E_beam, BB_angle, avg_mass);

  cout << "Loaded configuration parameters for kinematic " << kine << "." << endl;
  cout << "   Beam energy = " << E_beam << endl;
  cout << "   Bigbite angle = " << BB_angle << endl;
  cout << "   Measured Q2 = " << Q2 << endl;
  cout << "   Central Q2 = " << Q2_central << endl;
  cout << "   Rsf = " << Rsf << " +/- " << Rsf_err << endl << endl;

  cout << "Calculating shared analysis parameters..." << endl;

  // Calculate taus
  double tau_p = Q2 / (4.0 * M_p * M_p);
  double tau_n = Q2 / (4.0 * M_n * M_n);

  // Calculate epsilons
  double epsilon_n = 1.0 / (1.0 + 2.0 * (1.0 + tau_p) * pow(tan(BB_angle / 2.0), 2.0));
  double epsilon_p = 1.0 / (1.0 + 2.0 * (1.0 + tau_p) * pow(tan(BB_angle / 2.0), 2.0));

  // Get dipole form factor
  double GD = pow(1+Q2/GD_const,-2);

  cout << "   Tau proton = " << tau_p << endl;
  cout << "   Tau neutron = " << tau_n << endl;
  cout << "   Epsilon proton = " << epsilon_p << endl;
  cout << "   Epsilon neutron = " << epsilon_n << endl;
  cout << "   Dipole FF = " << GD << endl << endl;

  cout << "Calculating reduced cross section with simc parameterizations..." << endl;

  double simc_RCS_proton = getRCSsimc(tau_p,epsilon_p,GD,true);
  double simc_RCS_neutron = getRCSsimc(tau_n,epsilon_n,GD,false);

  double simc_RCSR = simc_RCS_neutron / simc_RCS_proton;

  double RCSR = Rsf * simc_RCSR;
  double RCSR_error = Rsf_err * simc_RCSR; //No added error to convert back to model independent RCSR

  cout << "   Reduced cross section, proton simc = " << simc_RCS_proton << endl;
  cout << "   Reduced cross section, neutron simc = " << simc_RCS_neutron << endl;
  cout << "   Reduced cross section ratio simc = " << simc_RCSR << endl;
  cout << "   Model independent Born + TPE cross section ratio using Rsf = " << Rsf << ":" << endl;
  cout << "      " << RCSR << endl << endl;

  double m_RCS_proton;
  double m_RCS_proton_error;

  std::pair m_pRCS_with_err = getRCSmdp_with_error(tau_p, epsilon_p, GD);
  std::pair m_pRCS_with_err_tpe = getRCSmdp_borntpe_with_error(tau_p, epsilon_p, GD);

  if(tpecorr){
    cout << "Calculating Arrington07 Born + TPE proton cross section..." << endl;
    m_RCS_proton = m_pRCS_with_err_tpe.first;
    m_RCS_proton_error = m_pRCS_with_err_tpe.second;
  }else{
    cout << "Calculating Arrington07 Born proton cross section..." << endl;
    m_RCS_proton = m_pRCS_with_err.first;
    m_RCS_proton_error = m_pRCS_with_err.second;
  }
  
  cout << "   Model reduced cross section, proton = " << m_RCS_proton << endl;
  cout << "   Model reduced cross section error, proton = " << m_RCS_proton_error << endl << endl;

  cout << "Calculating model-dependent measured Born + TPE neutron reduced cross section..." << endl;

  double md_RCS_neutron = RCSR * m_RCS_proton;
  double md_RCS_neutron_error = md_RCS_neutron * sqrt(pow(Rsf_err/Rsf,2)+pow(m_RCS_proton_error/m_RCS_proton,2));

  cout << "   Model dependent reduced cross section, neutron = " << md_RCS_neutron << endl << endl;

  cout << "Calculating Ye parameterization constants and Ye GEn..." << endl;

  double ye_tcut = 4 * M_pi * M_pi;
  double ye_z = calculateYeZ(Q2,ye_tcut,ye_tnot); //Base for GEn expansion
  double ye_x = calculateYeX(Q2); //Base for GEn error expansion
  double ye_GEn = calcYe(ye_z);
  double ye_GEn_error = calcYe_err(ye_x, GD);

  cout << "   Ye tcut = " << ye_tcut << endl;
  cout << "   Ye tnot = " << ye_tnot << endl;
  cout << "   Ye z = " << ye_z << endl;
  cout << "   Ye x = " << ye_x << endl;
  cout << "      Ye GEn = " << ye_GEn << endl;
  cout << "      Ye GEn error = " << ye_GEn_error << endl << endl;

  cout << "Calculating Ye GEn with 5%% TPE correction..." << endl; 
  
  //NEEDS WORK - this is just a placeholder from J. Arr private comm. Should verify necessary correction (if any) to Ye parameterization
  double tpe_ye_GEn = ye_GEn - 0.05*GD;

  cout << "   Ye TPE corrected GEn = " << tpe_ye_GEn << endl << endl;

  cout << "Extracting GMn..." << endl;

  // With model for neutron reduced cross section, solve for GMn
  //double GMn = getGMn_from_nRCS(md_RCS_neutron, tau_n, epsilon_n, GD, tpe_ye_GEn);
  std::pair<double,double> GMn_and_error = getGMn_and_error_from_nRCS(md_RCS_neutron, tau_n, epsilon_n, GD, md_RCS_neutron_error, tpe_ye_GEn, ye_GEn_error);
  double abs_mu_n = abs(mu_n);

  cout << "   GMn = " << GMn_and_error.first << endl;
  cout << "   dGMn = " << GMn_and_error.second << endl;
  cout << "   GMn/GD/mun = " << GMn_and_error.first/GD/abs_mu_n << endl;
  cout << "   dGMn/GD/mun = " << GMn_and_error.second/GD/abs_mu_n << endl;
  
  return;
}

// Helper function to calculate Q2
double calculateQ2(double E_beam, double BB_angle, double M) {
  double ePeak = E_beam / (1.0 + (2.0 * E_beam / M) * pow(sin(BB_angle / 2.0), 2.0));
  return 4 * E_beam * ePeak * pow(sin(BB_angle / 2.0), 2.0);
}

double calculateYeZ(double Q2, double tcut, double tnot){
  double numerator = sqrt(tcut+Q2) - sqrt(tcut-tnot);
  double denominator = sqrt(tcut+Q2) + sqrt(tcut-tnot);
  return numerator / denominator;
}

double calculateYeX(double Q2){
  double yex;
  if (Q2 == 1)
    yex = 0.00000001;
  else
    yex = log10(Q2);
  
  return yex;
}

// Function to calculate GEn using provided coefficients
double calcRiordan(double tau, double a1, double a2, double a3, double b1, double b2, double b3) {
  double numerator = (a1 * tau) + (a2 * tau * tau) + (a3 * tau * tau * tau);
  double denominator = 1 + (b1 * tau) + (b2 * tau * tau) + (b3 * tau * tau * tau);
  return numerator / denominator;
}

// double calcYe(double yez){
//   double GEn = 0.048919981 - 0.064525054*yez - 0.240825897*yez*yez + 0.392108745*yez*yez*yez + 0.300445259*yez*yez*yez*yez - 0.661888687*yez*yez*yez*yez*yez - 0.17563977*yez*yez*yez*yez*yez*yez + 0.624691724*yez*yez*yez*yez*yez*yez*yez - 0.077684299*yez*yez*yez*yez*yez*yez*yez*yez - 0.236003975*yez*yez*yez*yez*yez*yez*yez*yez*yez*yez + 0.090401973*yez*yez*yez*yez*yez*yez*yez*yez*yez*yez*yez;
//   return GEn;
// }

double calcYe(double yez){
  double GEn = 0.048919981 +
    -0.064525054 * yez +
    -0.240825897 * pow(yez, 2) +
    0.392108745 * pow(yez, 3) +
    0.300445259 * pow(yez, 4) +
    -0.661888687 * pow(yez, 5) +
    -0.17563977 * pow(yez, 6) +
    0.624691724 * pow(yez, 7) +
    -0.077684299 * pow(yez, 8) +
    -0.236003975 * pow(yez, 9) +
    0.090401973 * pow(yez, 10);
  return GEn;
}

double calcYe_err(double yex, double GD){
  double GEn_err_exp = -2.07194073 +
    1.13809127 * yex +
    1.01431277 * pow(yex, 2) +
    -0.31330138 * pow(yex, 3) +
    -0.273293676 * pow(yex, 4) +
    0.257350595 * pow(yex, 5) +
    -0.206042113 * pow(yex, 6) +
    -0.168497322 * pow(yex, 7) +
    0.137784515 * pow(yex, 8) +
    0.075759196 * pow(yex, 9) +
    -0.02675113 * pow(yex, 10) +
    -0.017525731 * pow(yex, 11) +
    0.000703582 * pow(yex, 12) +
    0.001479621 * pow(yex, 13) +
    0.000197375 * pow(yex, 14);

  double GEn_err = pow(10, GEn_err_exp) * GD;

  return GEn_err;
}

// Function to return reduced cross section for a given tau as calculated in simc
double getRCSsimc(double tau_value, double epsilon_value, double GD, bool isproton) {
  // GEp calculation (JJ Kelly)
  double GEp = (1.0 - 0.24 * tau_value) / (1.0 + 10.98 * tau_value + 12.82 * tau_value * tau_value + 21.97 * tau_value * tau_value * tau_value);

  // GMp/mup calculation (JJ Kelly)
  double GMp = (2.793 * (1 + 0.12 * tau_value)) / (1.0 + 10.97 * tau_value + 18.86 * tau_value * tau_value + 6.55 * tau_value * tau_value * tau_value);

  // GMn/mun calculation (JJ Kelly)
  double GMn = (-1.913 * (1.0 + 2.33 * tau_value)) / (1.0 + 14.72 * tau_value + 24.2 * tau_value * tau_value + 84.1 * tau_value * tau_value * tau_value);

  // GEn calculation using the provided Riordan function coefficients
  double GEn = GD * calcRiordan(tau_value, 1.52, 2.629, 3.055, 5.222, 0.04, 11.438);

  double RCS;

  if(isproton)
    RCS = (epsilon_value * GEp * GEp + tau_value * GMp * GMp) / (epsilon_value * (1.0 + tau_value));
  else
    RCS = (epsilon_value * GEn * GEn + tau_value * GMn * GMn) / (epsilon_value * (1.0 + tau_value));

  return RCS;
}

// Function to return model reduced cross section for proton using Arrington07 parameterization
std::pair<double,double> getRCSmdp_with_error(double tau_p, double epsilon_p, double GD) {
  // GEp calculation (Arrington07)
  double GEp = (1 + 3.439*tau_p - 1.602*tau_p*tau_p + 0.068*tau_p*tau_p*tau_p) / 
    (1 + 15.055*tau_p + 48.061*tau_p*tau_p + 99.304*tau_p*tau_p*tau_p + 0.012*tau_p*tau_p*tau_p*tau_p + 8.65*tau_p*tau_p*tau_p*tau_p*tau_p);

  // GMp/mup calculation (Arrington07)
  double GMp_over_mup = (1 - 1.465*tau_p + 1.26*tau_p*tau_p + 0.262*tau_p*tau_p*tau_p) / 
    (1 + 9.627*tau_p + 0.*tau_p*tau_p + 0.*tau_p*tau_p*tau_p + 11.179*tau_p*tau_p*tau_p*tau_p + 13.245*tau_p*tau_p*tau_p*tau_p*tau_p);

  double RCSmdp = (epsilon_p * GEp * GEp + tau_p * GMp_over_mup * GMp_over_mup * mu_p * mu_p) / (epsilon_p * (1.0 + tau_p));

  // Error calculation on model
  double Q2 = tau_p * (4.0 * M_p * M_p);
  double GEp_error = GE_over_GD_errslope * Q2 * GD;
  double GMp_error = GM_over_GDmup_errslope * Q2 * GD * mu_p;
  
  // Simplified error calculation (Assuming independent variables and linear propagation)
  double RCSmdp_error = std::sqrt(std::pow(2 * epsilon_p * GEp * GEp_error, 2) + 
                                  std::pow(2 * tau_p * GMp_over_mup * mu_p * GMp_error, 2)) / 
                                  (epsilon_p * (1.0 + tau_p));

  return {RCSmdp, RCSmdp_error};
}

// Function to return model reduced cross section for proton using Arrington07 parameterization and its error
std::pair<double, double> getRCSmdp_borntpe_with_error(double tau_p, double epsilon_p, double GD) {
  // Constants
  const double mu_p = 2.79284735; // Proton magnetic moment, assuming it's defined somewhere
  
  // GEp calculation (Arrington07)
  double GEp = (1 - 1.651*tau_p + 1.287*tau_p*tau_p - 0.185*tau_p*tau_p*tau_p) / 
    (1 + 9.531*tau_p + 0.591*tau_p*tau_p + 0.0*tau_p*tau_p*tau_p + 0.0*tau_p*tau_p*tau_p*tau_p + 4.994*tau_p*tau_p*tau_p*tau_p*tau_p);

  // GMp/mup calculation (Arrington07)
  double GMp_over_mup = (1 - 2.151*tau_p + 4.261*tau_p*tau_p + 0.159*tau_p*tau_p*tau_p) / 
    (1 + 8.647*tau_p + 0.001*tau_p*tau_p + 5.245*tau_p*tau_p*tau_p + 82.817*tau_p*tau_p*tau_p*tau_p + 14.191*tau_p*tau_p*tau_p*tau_p*tau_p);

  double RCSmdp = (epsilon_p * std::pow(GEp, 2) + tau_p * std::pow(GMp_over_mup * mu_p, 2)) / (epsilon_p * (1.0 + tau_p));

  // Error calculation on model
  double Q2 = tau_p * (4.0 * M_p * M_p);
  double GEp_error = GE_over_GD_errslope * Q2 * GD;
  double GMp_error = GM_over_GDmup_errslope * Q2 * GD * mu_p;
  
  // Simplified error calculation (Assuming independent variables and linear propagation)
  double RCSmdp_error = std::sqrt(std::pow(2 * epsilon_p * GEp * GEp_error, 2) + 
                                  std::pow(2 * tau_p * GMp_over_mup * mu_p * GMp_error, 2)) / 
                                  (epsilon_p * (1.0 + tau_p));

  return {RCSmdp, RCSmdp_error};
}

std::pair<double, double> getGMn_and_error_from_nRCS(double nRCS, double tau_value, double epsilon_value, double GD, double err_nRCS, double GEn, double err_GEn) {
    double GEn_riordan = GD * calcRiordan(tau_value, 1.52, 2.629, 3.055, 5.222, 0.04, 11.438) * 1.05;
    if(GEn == 0) {
        GEn = GEn_riordan;
    }

    double GMn = sqrt((nRCS * epsilon_value * (1.0 + tau_value) - epsilon_value * GEn * GEn) / tau_value);

    // Propagate errors
    double dGMn_dNRCs = (epsilon_value * (1.0 + tau_value)) / (2 * tau_value * GMn);
    double dGMn_dGEn = (-epsilon_value * 2 * GEn) / (tau_value * GMn);
    
    double err_GMn = sqrt(pow(dGMn_dNRCs * err_nRCS, 2) + pow(dGMn_dGEn * err_GEn, 2));

    return std::make_pair(GMn, err_GMn);
}
