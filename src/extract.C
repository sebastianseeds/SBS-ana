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

namespace extract{

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
    const double mu_p = 2.79284735; // Proton magnetic moment
  
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
  
    // // Simplified error calculation (Assuming independent variables and linear propagation)
    // double RCSmdp_error = std::sqrt(std::pow(2 * epsilon_p * GEp * GEp_error, 2) + 
    // 				    std::pow(2 * tau_p * GMp_over_mup * mu_p * GMp_error, 2)) / 
    //   (epsilon_p * (1.0 + tau_p));

    double RCSmdp_error = sqrt(
			       pow((2 * GEp * GEp_error) / (1.0 + tau_p), 2) +
			       pow((2 * tau_p * GMp_over_mup * mu_p * GMp_error) / (epsilon_p * (1.0 + tau_p)), 2)
			       );

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
  
  std::vector<std::pair<double,double>> extract_GMn_from_simc(double Rsf, double Rsf_err, double Q2, double E_beam, double BB_angle, bool tpecorr = true) {

    // Calculate average Q2
    double avg_mass = (M_p+M_n)/2;

    double Q2_central = calculateQ2(E_beam, BB_angle, avg_mass);

    if(Q2==0.)
      Q2 = Q2_central;

    // Calculate taus
    double tau_p = Q2 / (4.0 * M_p * M_p);
    double tau_n = Q2 / (4.0 * M_n * M_n);

    // Calculate epsilons
    double epsilon_n = 1.0 / (1.0 + 2.0 * (1.0 + tau_p) * pow(tan(BB_angle / 2.0), 2.0));
    double epsilon_p = 1.0 / (1.0 + 2.0 * (1.0 + tau_p) * pow(tan(BB_angle / 2.0), 2.0));

    // Get dipole form factor
    double GD = pow(1+Q2/GD_const,-2);

    double simc_RCS_proton = getRCSsimc(tau_p,epsilon_p,GD,true);
    double simc_RCS_neutron = getRCSsimc(tau_n,epsilon_n,GD,false);

    double simc_RCSR = simc_RCS_neutron / simc_RCS_proton;

    double RCSR = Rsf * simc_RCSR;
    double RCSR_error = Rsf_err * simc_RCSR; //No added error to convert back to model independent RCSR

    double m_RCS_proton;
    double m_RCS_proton_error;

    std::pair m_pRCS_with_err = getRCSmdp_with_error(tau_p, epsilon_p, GD);
    std::pair m_pRCS_with_err_tpe = getRCSmdp_borntpe_with_error(tau_p, epsilon_p, GD);

    if(tpecorr){
      m_RCS_proton = m_pRCS_with_err_tpe.first;
      m_RCS_proton_error = m_pRCS_with_err_tpe.second;
    }else{
      m_RCS_proton = m_pRCS_with_err.first;
      m_RCS_proton_error = m_pRCS_with_err.second;
    }
  
    double md_RCS_neutron = RCSR * m_RCS_proton;
    double md_RCS_neutron_error = md_RCS_neutron * sqrt(pow(Rsf_err/Rsf,2)+pow(m_RCS_proton_error/m_RCS_proton,2));

    double ye_tcut = 4 * M_pi * M_pi;
    double ye_z = calculateYeZ(Q2,ye_tcut,ye_tnot); //Base for GEn expansion
    double ye_x = calculateYeX(Q2); //Base for GEn error expansion
    double ye_GEn = calcYe(ye_z);
    double ye_GEn_error = calcYe_err(ye_x, GD);
  
    //NEEDS WORK - this is just a placeholder from J. Arr private comm. Should verify necessary correction (if any) to Ye parameterization
    double tpe_ye_GEn = ye_GEn - 0.05*GD;

    // With model for neutron reduced cross section, solve for GMn
    //double GMn = getGMn_from_nRCS(md_RCS_neutron, tau_n, epsilon_n, GD, tpe_ye_GEn);
    std::pair<double,double> GMn_and_error = getGMn_and_error_from_nRCS(md_RCS_neutron, tau_n, epsilon_n, GD, md_RCS_neutron_error, tpe_ye_GEn, ye_GEn_error);
    double abs_mu_n = abs(mu_n);

    //return the vector of pairs with GMn (first) GMn error (second), base then with dipole and mun corr
    std::pair<double,double> GMn_and_error_corr = std::make_pair(GMn_and_error.first/GD/abs_mu_n, GMn_and_error.second/GD/abs_mu_n);

    std::pair<double,double> GD_and_mun = std::make_pair(GD, abs_mu_n);

    std::pair<double,double> NRCS_and_error = std::make_pair(md_RCS_neutron, md_RCS_neutron_error);

    std::pair<double,double> model_PRCS_and_error = std::make_pair(m_RCS_proton, m_RCS_proton_error);

    std::vector<std::pair<double,double>> allPar;
    allPar.push_back(GMn_and_error);
    allPar.push_back(GMn_and_error_corr);
    allPar.push_back(GD_and_mun);
    allPar.push_back(NRCS_and_error);
    allPar.push_back(model_PRCS_and_error);
    
    return allPar;
  }

  // * Function to extract GMn from SIMC with detailed error components.
  // * This function calculates GMn, GMn corrected for dipole and mu_n,
  // * and their respective errors including statistical, systematic, and model errors.
  // *
  // * @param Rsf Reported scale factor from experimental results.
  // * @param Rsf_err_stat Statistical error in the reported scale factor.
  // * @param Rsf_err_syst Systematic error in the reported scale factor.
  // * @param Q2 Four-momentum transfer squared. If zero, it will be calculated.
  // * @param E_beam Beam energy in GeV.
  // * @param BB_angle BigBite angle in radians.
  // * @param tpecorr Flag to use TPE-corrected model or not.
  // * @return Vector of tuples. Each tuple contains a measurement and its error components.
   
  std::vector<std::tuple<double, double, double, double, double>> extract_GMn_from_simc(double Rsf, double Rsf_err_stat, double Rsf_err_syst, double Q2, double E_beam, double BB_angle, bool tpecorr = true) {
    // Calculate average Q2 if not provided
    double avg_mass = (M_p + M_n) / 2;
    double Q2_central = calculateQ2(E_beam, BB_angle, avg_mass);
    if (Q2 == 0.) Q2 = Q2_central;

    // Calculate tau for proton and neutron
    double tau_p = Q2 / (4.0 * M_p * M_p);
    double tau_n = Q2 / (4.0 * M_n * M_n);
    
    // Calculate epsilon (virtual photon polarization)
    double epsilon = 1.0 / (1.0 + 2.0 * (1.0 + tau_p) * pow(tan(BB_angle / 2.0), 2.0));
    
    // Calculate dipole form factor GD
    double GD = pow(1 + Q2 / GD_const, -2);

    // Calculate RCS for proton and neutron using SIMC model
    double simc_RCS_proton = getRCSsimc(tau_p, epsilon, GD, true);
    double simc_RCS_neutron = getRCSsimc(tau_n, epsilon, GD, false);
    
    // Scale MC reduced cross section ratio by reported scale factor Rsf to get model independent ratio
    double simc_RCSR = simc_RCS_neutron / simc_RCS_proton;
    double RCSR = Rsf * simc_RCSR; //Rsf * simc_RCSR * (epsilon_p*(1+tau_p)) / (epsilon_n*(1+tau_n)) 

    // Calculate errors for RCSR
    double RCSR_error_stat = Rsf_err_stat * simc_RCSR;
    double RCSR_error_syst = Rsf_err_syst * simc_RCSR;
    double RCSR_error_model = sqrt(RCSR_error_stat * RCSR_error_stat + RCSR_error_syst * RCSR_error_syst);

    // Get model-dependent proton RCS and its errors using Arrington07 parameterization
    auto [m_RCS_proton, m_RCS_proton_error] = tpecorr ? getRCSmdp_borntpe_with_error(tau_p, epsilon, GD) : getRCSmdp_with_error(tau_p, epsilon, GD);
    double md_RCS_neutron = RCSR * m_RCS_proton;
    double md_RCS_neutron_error_stat = md_RCS_neutron * (RCSR_error_stat / RCSR);
    double md_RCS_neutron_error_syst = md_RCS_neutron * (RCSR_error_syst / RCSR);
    double md_RCS_neutron_error_model = sqrt(pow(m_RCS_proton_error / m_RCS_proton, 2)) * md_RCS_neutron;
    double md_RCS_neutron_error_total = sqrt(md_RCS_neutron_error_stat * md_RCS_neutron_error_stat + md_RCS_neutron_error_syst * md_RCS_neutron_error_syst + md_RCS_neutron_error_model * md_RCS_neutron_error_model);

    // Get GEn and model error
    double md_ye_GEn = calcYe(calculateYeZ(Q2, 4 * M_pi * M_pi, ye_tnot));
    double md_ye_GEn_err = calcYe_err(calculateYeX(Q2), GD);

    // Calculate GMn from neutron RCS, proton RCS, GEn, tau, and epsilon
    double GMn = sqrt((md_RCS_neutron * epsilon * (1.0 + tau_n) - epsilon * md_ye_GEn * md_ye_GEn) / tau_n);

    // Derivative of GMn w.r.t. md_RCS_neutron
    double partial_GMn_md_RCS_neutron = (epsilon * (1.0 + tau_n)) / (2 * tau_n * GMn);

    // Derivative of GMn w.r.t. md_ye_GEn
    double partial_GMn_md_ye_GEn = (-2 * epsilon * md_ye_GEn) / (2 * tau_n * GMn);

    // Statistical, systematic, and model errors for GMn
    double GMn_error_stat = abs(partial_GMn_md_RCS_neutron * md_RCS_neutron_error_stat);
    double GMn_error_syst = abs(partial_GMn_md_RCS_neutron * md_RCS_neutron_error_syst);
    double GMn_error_model = sqrt(pow(partial_GMn_md_RCS_neutron * m_RCS_proton_error, 2) + pow(partial_GMn_md_ye_GEn * md_ye_GEn_err, 2));

    // Total error for GMn using quadrature sum of all components
    double GMn_error_total = sqrt(GMn_error_stat * GMn_error_stat + GMn_error_syst * GMn_error_syst + GMn_error_model * GMn_error_model);

    // Organize the results into a tuple
    std::tuple<double, double, double, double, double> GMn_results = std::make_tuple(
										     GMn, GMn_error_stat, GMn_error_syst, GMn_error_model, GMn_error_total);

    // Optionally, calculate and add corrected GMn
    double GMn_corrected = GMn / (GD * fabs(mu_n));
    double GMn_corrected_error = GMn_error_total / (GD * fabs(mu_n));  // Simplified error correction assuming mu_n and GD are constants

    std::tuple<double, double, double, double, double> GMn_corrected_results = std::make_tuple(
											       GMn_corrected, GMn_error_stat / (GD * fabs(mu_n)), GMn_error_syst / (GD * fabs(mu_n)), GMn_error_model / (GD * fabs(mu_n)), GMn_corrected_error);

    // Return a vector of tuples containing the results
    std::vector<std::tuple<double, double, double, double, double>> results;
    results.push_back(GMn_results);
    results.push_back(GMn_corrected_results);

    return results;

    // // Calculate GMn using neutron RCS model and propagate errors
    // auto [GMn, GMn_error_stat, GMn_error_syst, GMn_error_model] = getGMn_and_error_from_nRCS(md_RCS_neutron, tau_n, epsilon, GD, md_RCS_neutron_error_total, calcYe(calculateYeZ(Q2, 4 * M_pi * M_pi, ye_tnot)), calcYe_err(calculateYeX(Q2), GD));
    // double GMn_corrected = GMn / (GD * fabs(mu_n));
    // double GMn_corrected_error_stat = GMn_error_stat / (GD * fabs(mu_n));
    // double GMn_corrected_error_syst = GMn_error_syst / (GD * fabs(mu_n));
    // double GMn_corrected_error_model = GMn_error_model / (GD * fabs(mu_n));
    // double GMn_corrected_error_total = sqrt(GMn_corrected_error_stat * GMn_corrected_error_stat + GMn_corrected_error_syst * GMn_corrected_error_syst + GMn_corrected_error_model * GMn_corrected_error_model);

    // Organize results into a vector of tuples
    // std::vector<std::tuple<double, double, double, double, double>> results;
    // results.push_back(std::make_tuple(GMn, GMn_error_stat, GMn_error_syst, GMn_error_model, sqrt(GMn_error_stat * GMn_error_stat + GMn_error_syst * GMn_error_syst + GMn_error_model * GMn_error_model)));  // GMn
    // results.push_back(std::make_tuple(GMn_corrected, GMn_corrected_error_stat, GMn_corrected_error_syst, GMn_corrected_error_model, GMn_corrected_error_total)); // Corrected GMn
    // results.push_back(std::make_tuple(md_RCS_neutron, 0.0, md_RCS_neutron_error, m_RCS_proton_error, md_RCS_neutron_error)); // Neutron RCS
    // results.push_back(std::make_tuple(m_RCS_proton, 0.0, m_RCS_proton_error, m_RCS_proton_error, m_RCS_proton_error)); // Model Proton RCS
    // results.push_back(std::make_tuple(GD, 0.0, 0.0, 0.0, 0.0));  // GD
    // results.push_back(std::make_tuple(fabs(mu_n), 0.0, 0.0, 0.0, 0.0)); // mu_n

    // return results;
  }
}
