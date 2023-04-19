#ifndef VARS_H
#define VARS_H

#include <cmath>
#include "TVector3.h"
#include "TLorentzVector.h"

#include "physconst.h"
#include "econst.h"

namespace vars{

  // central e' momentum
  Double_t pcentral( Double_t ebeam, Double_t etheta, std::string Ntype );       // GeV 
  Double_t etheta( TLorentzVector Peprime );  // scattering angle (rad)
  Double_t ephi( TLorentzVector Peprime );    // angle of scattering plane (rad)

  // construct target nucleon 4-momentum (assuming at rest)
  void setPN( std::string Ntype, TLorentzVector &PN );

  // expected recoil nucleon momentum
  Double_t pN_expect( Double_t nu,             // energy of the virtual photon (GeV)
		      std::string Ntype );    // type of nucleon in the reaction
		   
  // projected nucleon 3 vector using track angles (pNhat)
  TVector3 pNhat_track( Double_t Ntheta,      // recoil nucleon theta (rad)
		       Double_t Nphi );       // recoil nucleon phi (rad) 

  // construct hcal coordinate system in hall coordinate system
  void sethcalaxes( Double_t sbstheta_rad,           // SBS angle (rad)
		    vector<TVector3> &hcal_axes );  // hcal axes x,y,z

  // get the expected x and y positions of the recoil nucleon on face of hcal
  void getxyhcalexpect( TVector3 vertex,                 // vertex vector (in hall system)
		        TVector3 pNhat,                  // projected q vector
		        TVector3 hcal_origin,            // hcal origin vector (in hall system)
		        vector<TVector3> hcal_axes,      // hcal axes (in hall system)
		        vector<Double_t> &xyhcalexpect );   // expected x and y positions       

  // physics
  Double_t Q2( Double_t ebeam, Double_t eeprime, Double_t etheta );                // GeV, GeV, rad
  Double_t W2( Double_t ebeam, Double_t eeprime, Double_t Q2, std::string Ntype ); // GeV, GeV, GeV2
  Double_t W( Double_t ebeam, Double_t eeprime, Double_t Q2, std::string Ntype );  // GeV, GeV, GeV2

  // beam luminosity
  Double_t luminosity( Double_t ibeam, std::string targetType );               // A

  // beam energy vertex correction (energy loss through target)
  Double_t ebeam_c( Double_t ebeam, Double_t vz, std::string target ); //GeV  
  Double_t ebeam_o( Double_t ebeam_c, Double_t etheta, std::string target ); //GeV  
  void setpbeam( Double_t ebeam, TLorentzVector &pbeam ); //Beam electron four momentum

}

#endif
