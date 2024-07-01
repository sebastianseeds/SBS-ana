// S.Seeds. 3.1.24 This script will estimate Q2 for p/n with bb central angle and calculate GEp, GMp, GEn, and GMn with parameterizations in simc.

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
const double GD_const = 0.71; // Dipole form factor constant
const double mu_p = 2.79284735; //Proton magnetic moment
const double mu_n = -1.9130427; //Neutron anomalous magnetic moment

// Forward declarations
double CalculateQ2(double E_beam, double BB_angle, double M);
double calcRiordan(double tau, double a1, double a2, double a3, double b1, double b2, double b3);
//void CalculateFormFactors(double E_beam, double BB_angle);
void CalculateFormFactors(double E_beam, double BB_angle, double Q2); // Overload

// Main. kine:kinematic, Q2p:Q2 override
void calc_simc_FF(int kine, double Q2p=0) {

  // Set up configuration and tune objects to load analysis parameters
  SBSconfig config(kine,0);

  // Obtain configuration pars from config class
  double E_beam = config.GetEbeam();
  double BB_angle = config.GetBBtheta_rad(); // In radians

  // Calculate average Q2 for proton
  double avg_mass = (M_p+M_n)/2;

  double Q2 = CalculateQ2(E_beam, BB_angle, avg_mass);

  if( Q2p!=0 )
    Q2 = Q2p;

  cout << "Loaded configuration parameters for kinematic " << kine << "." << endl;
  cout << "   Beam energy = " << E_beam << endl;
  cout << "   Bigbite angle = " << BB_angle << endl;
  cout << "   ";
  if( Q2p==0 )
    cout << "Central ";
  cout << "Q2 = " << Q2 << endl << endl;

  cout << "Calculating Sachs form factors and Born cross section with simc parameterizations..." << endl;

  cout << "   Central Q2 values:" << endl;
  CalculateFormFactors(E_beam, BB_angle, Q2);

  return;
}

// Helper function to calculate Q2
double CalculateQ2(double E_beam, double BB_angle, double M) {
  double ePeak = E_beam / (1.0 + (2.0 * E_beam / M) * pow(sin(BB_angle / 2.0), 2.0));
  return 4 * E_beam * ePeak * pow(sin(BB_angle / 2.0), 2.0);
}

// Function to calculate GEn using provided coefficients
double calcRiordan(double tau, double a1, double a2, double a3, double b1, double b2, double b3) {
  double numerator = (a1 * tau) + (a2 * tau * tau) + (a3 * tau * tau * tau);
  double denominator = 1 + (b1 * tau) + (b2 * tau * tau) + (b3 * tau * tau * tau);
  return numerator / denominator;
}

// Function to calculate form factors as they are calculated in simc with given Q2
void CalculateFormFactors(double E_beam, double BB_angle, double Q2) {

  // Calculate taus
  double tau_p = Q2 / (4.0 * M_p * M_p);
  double tau_n = Q2 / (4.0 * M_n * M_n);

  // Calculate epsilons
  double epsilon_n = 1.0 / (1.0 + 2.0 * (1.0 + tau_p) * pow(tan(BB_angle / 2.0), 2.0));
  double epsilon_p = 1.0 / (1.0 + 2.0 * (1.0 + tau_p) * pow(tan(BB_angle / 2.0), 2.0));

  // Get dipole form factor
  double GD = pow(1+Q2/GD_const,-2);

  // GEp calculation (JJ Kelly)
  double GEp = (1.0 - 0.24 * tau_p) / (1.0 + 10.98 * tau_p + 12.82 * tau_p * tau_p + 21.97 * tau_p * tau_p * tau_p);

  // GMp calculation (JJ Kelly)
  double GMp = (2.79 * (1 + 0.12 * tau_p)) / (1.0 + 10.97 * tau_p + 18.86 * tau_p * tau_p + 6.55 * tau_p * tau_p * tau_p);

  // GMn calculation (JJ Kelly)
  double GMn = (-1.913 * (1.0 + 2.33 * tau_n)) / (1.0 + 14.72 * tau_n + 24.2 * tau_n * tau_n + 84.1 * tau_n * tau_n * tau_n);

  // GEn calculation using the provided Riordan function coefficients
  double GEn = GD * calcRiordan(tau_n, 1.52, 2.629, 3.055, 5.222, 0.04, 11.438);

  // Calculate reduced cross sections defined with 1/(1+tau) factor included
  double reducedCrossSection_n = (epsilon_n * GEn * GEn + tau_n * GMn * GMn) / (epsilon_n * (1.0 + tau_n));
  double reducedCrossSection_p = (epsilon_p * GEp * GEp + tau_p * GMp * GMp) / (epsilon_p * (1.0 + tau_p));

  // Output the results
  cout << "   GEp / GD: " << GEp / GD << ", GMp / (mu_p * GD): " << GMp / (mu_p * GD) << ", GEn / GD: " << GEn / GD << ", GMn / (mu_n * GD): " << GMn / (mu_n * GD) << endl;
  cout << "   GD: " << GD << endl;  
  cout << "      Reduced cross section neutron: " << reducedCrossSection_n << endl;
  cout << "      Reduced cross section proton: " << reducedCrossSection_p << endl;
  cout << endl;
  cout << "      Reduced cross section n:p ratio, simc: " << reducedCrossSection_n/reducedCrossSection_p << endl << endl;

}

// // Overload function to calculate form factors as they are calculated in simc with central Q2 from kinematic settings
// void CalculateFormFactors(double E_beam, double BB_angle) {


//   // Calculate taus
//   double tau_p = Q2 / (4.0 * M_p * M_p);
//   double tau_n = Q2 / (4.0 * M_n * M_n);

//   // Calculate epsilons
//   double epsilon_n = 1.0 / (1.0 + 2.0 * (1.0 + tau_p) * pow(tan(BB_angle / 2.0), 2.0));
//   double epsilon_p = 1.0 / (1.0 + 2.0 * (1.0 + tau_p) * pow(tan(BB_angle / 2.0), 2.0));

//   // Get dipole form factor
//   double GD = pow(1+Q2/GD_const,-2);

//   // GEp calculation (JJ Kelly)
//   double GEp = (1.0 - 0.24 * tau_p) / (1.0 + 10.98 * tau_p + 12.82 * tau_p * tau_p + 21.97 * tau_p * tau_p * tau_p);

//   // GMp calculation (JJ Kelly)
//   double GMp = (2.79 * (1.0 + 0.12 * tau_p)) / (1.0 + 10.97 * tau_p + 18.86 * tau_p * tau_p + 6.55 * tau_p * tau_p * tau_p);

//   // GMn calculation (JJ Kelly)
//   double GMn = (-1.913 * (1.0 + 2.33 * tau_n)) / (1.0 + 14.72 * tau_n + 24.2 * tau_n * tau_n + 84.1 * tau_n * tau_n * tau_n);

//   // GEn calculation using the provided Riordan function coefficients
//   double GEn = calcRiordan(tau_n, 1.52, 2.629, 3.055, 5.222, 0.04, 11.438);

//   // Calculate reduced cross sections defined with 1/(1+tau) factor included
//   double reducedCrossSection_n = (epsilon_n * GEn * GEn + tau_n * GMn * GMn) / (epsilon_n * (1.0 + tau_n));
//   double reducedCrossSection_p = (epsilon_p * GEp * GEp + tau_p * GMp * GMp) / (epsilon_p * (1.0 + tau_p));

//   // Output the results
//   cout << "   GEp / GD: " << GEp / GD << ", GMp / GD: " << GMp / GD << ", GEn / GD: " << GEn / GD << ", GMn / GD: " << GMn / GD << endl;
//   cout << "   GD: " << GD << endl;  
//   cout << "      Reduced cross section neutron: " << reducedCrossSection_n << endl;
//   cout << "      Reduced cross section proton: " << reducedCrossSection_p << endl;
//   cout << endl;
//   cout << "      Reduced cross section n:p ratio, simc: " << reducedCrossSection_n/reducedCrossSection_p << endl << endl;

// }
