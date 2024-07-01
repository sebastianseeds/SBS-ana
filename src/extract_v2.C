// S.Seeds. 3.24.24 This namespace will extract GMn for SBS experiment 12-09-019 using the following:
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

#include "../include/extract.h"

#include <cmath>
#include <vector>
#include <tuple>

namespace extract {

  const double M_p = 0.938272; // Proton mass in GeV
  const double M_n = 0.939565; // Neutron mass in GeV
  const double Pi = acos(-1);
  const double alpha = 1.0 / 137.0; // Fine-structure constant
  const double GD_const = 0.71; // Dipole form factor range factor
  const double mu_p = 2.79284735; // Proton magnetic moment
  const double mu_n = -1.9130427; // Neutron anomalous magnetic moment
  const double M_pi = 0.13957; // Charged pion mass in GeV
  const double ye_tnot = -0.7; // Ye parameterization constant t0

  double calculateQ2(double E_beam, double BB_angle, double M) {
    return 4 * E_beam * E_beam / (1.0 + (2.0 * E_beam / M) * pow(sin(BB_angle / 2.0), 2.0)) * pow(sin(BB_angle / 2.0), 2.0);
  }

  double calculateYeZ(double Q2, double tcut, double tnot) {
    return (sqrt(tcut + Q2) - sqrt(tcut - tnot)) / (sqrt(tcut + Q2) + sqrt(tcut - tnot));
  }

  double calculateYeX(double Q2) {
    return log10(Q2);
  }

  double calcRiordan(double tau, double a1, double a2, double a3, double b1, double b2, double b3) {
    return ((a1 * tau) + (a2 * pow(tau, 2)) + (a3 * pow(tau, 3))) / (1 + (b1 * tau) + (b2 * pow(tau, 2)) + (b3 * pow(tau, 3)));
  }

  double calcYe(double yez) {
    return 0.048919981 - 0.064525054 * yez - 0.240825897 * pow(yez, 2) + 0.392108745 * pow(yez, 3) +
           0.300445259 * pow(yez, 4) - 0.661888687 * pow(yez, 5) - 0.17563977 * pow(yez, 6) +
           0.624691724 * pow(yez, 7) - 0.077684299 * pow(yez, 8) - 0.236003975 * pow(yez, 9) +
           0.090401973 * pow(yez, 10);
  }

  double calcYe_err(double yex, double GD) {
    double GEn_err_exp = -2.07194073 + 1.13809127 * yex + 1.01431277 * pow(yex, 2) - 0.31330138 * pow(yex, 3) -
                         0.273293676 * pow(yex, 4) + 0.257350595 * pow(yex, 5) - 0.206042113 * pow(yex, 6) -
                         0.168497322 * pow(yex, 7) + 0.137784515 * pow(yex, 8) + 0.075759196 * pow(yex, 9) -
                         0.02675113 * pow(yex, 10) - 0.017525731 * pow(yex, 11) + 0.000703582 * pow(yex, 12) +
                         0.001479621 * pow(yex, 13) + 0.000197375 * pow(yex, 14);
    return pow(10, GEn_err_exp) * GD;
  }

  // Here you should ensure the function correctly calculates and returns the model cross section and its error
  std::pair<double, double> getRCSmdp_with_error(double tau, double epsilon, double GD, bool isProton) {
    // Simplified model: actual implementation should correctly calculate GEp, GMp or GEn, GMn based on isProton
    double formFactor = 1.0 / (1.0 + GD * tau); // Placeholder for actual calculation
    double error = formFactor * 0.05; // Placeholder for error calculation
    return {formFactor, error};
  }

  // This function should correctly decompose and handle the cross-section model calculations and error propagations
  std::vector<std::tuple<double, double, double, double, double>> extract_GMn_from_simc(double Rsf, double Rsf_err_stat, double Rsf_err_syst, double Q2, double E_beam, double BB_angle, bool tpecorr) {
    double Q2_calculated = calculateQ2(E_beam, BB_angle, (M_p + M_n) / 2);
    if (Q2 == 0) Q2 = Q2_calculated;

    double tau = Q2 / (4 * M_p * M_p); // Assuming proton, adjust for neutron as needed
    double epsilon = 1 / (1 + 2 * (1 + tau) * pow(tan(BB_angle / 2.0), 2));

    double GD = pow(1 + Q2 / GD_const, -2);

    // Retrieve cross-section and its error from model
    auto [RCS, RCS_error] = getRCSmdp_with_error(tau, epsilon, GD, true); // true for proton, adjust as needed

    // Calculate corrected values and errors
    double corrected_RCS = RCS * Rsf;
    double corrected_RCS_error_stat = RCS * Rsf_err_stat;
    double corrected_RCS_error_syst = RCS * Rsf_err_syst;
    double corrected_RCS_error_total = sqrt(pow(corrected_RCS_error_stat, 2) + pow(corrected_RCS_error_syst, 2) + pow(RCS_error, 2));

    std::vector<std::tuple<double, double, double, double, double>> results;
    results.push_back(std::make_tuple(corrected_RCS, corrected_RCS_error_stat, corrected_RCS_error_syst, RCS_error, corrected_RCS_error_total));

    return results;
  }
}
