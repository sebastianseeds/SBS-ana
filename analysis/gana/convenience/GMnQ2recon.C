// S.Seeds. 11.15.22 This script will estimate Q2 for p/n with bb angle 

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
void GMnQ2recon( Double_t E_beam, Double_t BB_angle, Double_t BB_distance, Double_t Q2, Double_t hcal_angle, Double_t hcal_distance ){

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
  
  Double_t Q2_recon_p[7] = {0.};
  Double_t Q2_recon_n[7] = {0.};

  for(int shcol=0; shcol<7; shcol++){
    Double_t ePeak_ep = E_beam/( 1. + (2.*E_beam/M_p)*pow(sin(effective_BBang[shcol]/2.), 2.) );
    Double_t ePeak_en = E_beam/( 1. + (2.*E_beam/M_n)*pow(sin(effective_BBang[shcol]/2.), 2.) );

    //////////////////
    Q2_recon_p[shcol] = 4*E_beam*ePeak_ep*(pow(sin(effective_BBang[shcol]/2.), 2.));
    Q2_recon_n[shcol] = 4*E_beam*ePeak_en*(pow(sin(effective_BBang[shcol]/2.), 2.));
    ///////////////////

  }

  cout << endl << endl;

  cout << "Q2 by SH col (from proton) = ";
  for ( int i=0; i<7; ++i )
    cout << Q2_recon_p[i] << " ";

  cout << endl;

  cout << "Q2 by SH col (from neutron)= ";
  for ( int i=0; i<7; ++i )
    cout << Q2_recon_n[i] << " ";

  cout << endl;


}

