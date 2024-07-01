// S.Seeds. 11.15.22 This script will estimate the energy of the elastically scattered electron and proton (e-p) taking into consideration the angular acceptance of BBSH for a proton or neutron target. Made for GMn, but generalizable for all SBS experiments

#include <TMath.h>
#include <iostream>
#include <iomanip>
#include <vector>

using namespace std;

const Double_t M_p = 0.938272; // GeV
const Double_t M_n = 0.939565; // GeV
const Double_t M_e = 0.00051; // GeV
const Double_t ff = 0.05; // Added arbitrary smearing factor to account for beam energy fluctuations and fermi motion in downstream estimates
const Double_t hcal_width = 1.70434; // m
const Double_t hcal_sampfrac = 0.0795; // m
const Double_t hcal_threshconv = 6.914; // MeV/mV
const Double_t bbcal_threshconv = 7.2; // MeV/mV

//MAIN. E_beam=beam energy, BB_angle=angle of BigBite spectrometer wrt downstream beamline, BB_distance=distance to BigBite from target, Q2=Q2, hcal_angle=angle of HCal wrt downstream beamline, hcal_distance=distance to the HCal from target
void GMnElasPeak( Double_t E_beam, Double_t BB_angle, Double_t BB_distance, Double_t Q2, Double_t hcal_angle, Double_t hcal_distance ){

  hcal_angle = hcal_angle*(TMath::Pi()/180.);
  Double_t hcal_minang = hcal_angle - ( (hcal_width/2)/hcal_distance ); //approx with arclength
  Double_t hcal_maxang = hcal_angle + ( (hcal_width/2)/hcal_distance ); //approx with arclength

  Double_t sh_faceDist = 3.1 + BB_distance; //1.2m to GEMs, another 1.9m to BBCal
  
  Double_t sh_ypos[7] = {-0.2565, -0.171, -.0855, 0.0, 0.0855, 0.171, 0.2565}; //Relative positions of shower columns
  Double_t effective_BBang[7] = {0.};
  Double_t hcal_p_projang[7] = {0.};
  Double_t hcal_n_projang[7] = {0.};

  Double_t deltaBBang = 0.;
  for(int shcol=0; shcol<7; shcol++){
    effective_BBang[shcol] = (sh_ypos[shcol]/sh_faceDist) + BB_angle*(TMath::Pi()/180.); //with deg to rad
  }

  cout << endl;
  
  Double_t elasPeak_e[7] = {0.};
  Double_t elasPeak_p[7] = {0.};
  Double_t elasPeak_n[7] = {0.};
  Double_t epsilon_p[7] = {0.};
  Double_t epsilon_n[7] = {0.};
  Double_t nu_p[7] = {0.}; //keeping around in arrays for possible readout in future versions. Just KE where KE = E_beam - E_e'
  Double_t nu_n[7] = {0.};
  Double_t p_p[7] = {0.};
  Double_t p_n[7] = {0.};
  Double_t p_e[7] = {0.};
  Double_t hcal_minKE = 12345.;
  bool hcalON[7] = {false};

  for(int shcol=0; shcol<7; shcol++){
    //Proton
    Double_t ePeak_e = E_beam/( 1. + (2.*E_beam/M_p)*pow(sin(effective_BBang[shcol]/2.), 2.) );

    elasPeak_e[shcol] = ePeak_e;
    nu_p[shcol] = E_beam-ePeak_e;
    p_p[shcol] = sqrt(pow(nu_p[shcol],2.)+2.*M_p*nu_p[shcol]);

    Double_t phel = sqrt( pow(nu_p[shcol],2.)+2.*M_p*nu_p[shcol] );

    epsilon_p[shcol] = pow( 1+2*( 1+Q2/(4*M_p) )*pow( tan(effective_BBang[shcol]/2),2 ), -1 );
    elasPeak_p[shcol] = nu_p[shcol]+M_p;
    hcal_p_projang[shcol] = acos( (E_beam - ePeak_e*cos(effective_BBang[shcol]))/phel );
    
    //Neutron
    ePeak_e = E_beam/( 1. + ( 2.*E_beam/M_n )*pow( sin(effective_BBang[shcol]/2. ), 2.) );
    nu_n[shcol] = E_beam-ePeak_e;
    p_n[shcol] = sqrt(pow(nu_n[shcol],2.)+2.*M_n*nu_n[shcol]);

    phel = sqrt( pow(nu_n[shcol],2.)+2.*M_n*nu_n[shcol] );

    epsilon_n[shcol] = pow( 1+2*( 1+Q2/(4*M_n) )*pow( tan(effective_BBang[shcol]/2),2 ), -1) ;
    elasPeak_n[shcol] = nu_n[shcol]+M_n;
    hcal_n_projang[shcol] = acos( (E_beam - ePeak_e*cos(effective_BBang[shcol]))/phel );

    hcalON[shcol] = hcal_p_projang[shcol]>hcal_minang&&hcal_p_projang[shcol]<hcal_maxang;

    //Check for lowest KE in hcal
    if( nu_p[shcol]<hcal_minKE && hcalON[shcol] ) hcal_minKE = nu_p[shcol];
    if( nu_n[shcol]<hcal_minKE && hcalON[shcol] ) hcal_minKE = nu_n[shcol];

  }

  cout << endl;

  cout << "///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "//   BBSh col  |  e' angle   |  e' p      |  elas proton KE  |  elas neutron KE |  had proj ang  |  on hcal? | p epsilon | n epsilon //" << endl;

  for(int shcol=0; shcol<7; shcol++){
    
    //bool hcalON = hcal_p_projang[shcol]>hcal_minang&&hcal_p_projang[shcol]<hcal_maxang;

    cout << "//___________________________________________________________________________________________________________________________________//" << endl;

    std::cout << std::setprecision(5) << std::fixed;

    cout << "//   " << shcol <<  "         |  " << effective_BBang[shcol]*(180./TMath::Pi())  << "   |  " << elasPeak_e[shcol]  << "   |  " << nu_p[shcol]  << "         |  " << nu_n[shcol] << "         |  " << hcal_p_projang[shcol]*(180./TMath::Pi()) << "      |  " << hcalON[shcol] << "        |  " << epsilon_p[shcol] << "  |  " << epsilon_n[shcol] << "  //" << endl; 

  }

  cout << "///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////" << endl;

  Double_t hcal_ms = hcal_minKE*hcal_sampfrac*1000.; //convert from GeV to MeV

  cout << endl << endl << "Lowest energy sampled in HCal with estimated smearing (MeV): " << hcal_ms - hcal_ms*ff << "." << endl;
  cout << endl << endl << "Suggested hcal threshold with estimated smearing: " << (hcal_ms - hcal_ms*ff)/hcal_threshconv << "." << endl;

  cout << endl << endl << "Lowest energy sampled in HCal (MeV): " << hcal_ms << "." << endl;
  cout << endl << endl << "Suggested hcal threshold: " << hcal_ms/hcal_threshconv << "." << endl << endl;
  cout << endl << endl << "Central proton momentum: " << p_p[3] << endl;
  cout << endl << endl << "Central neutron momentum: " << p_n[3] << endl;
  cout << endl << endl << "Q^2: " << Q2 << endl << endl << endl;

}

