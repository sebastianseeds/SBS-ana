//Seeds - Step through calculation of xyexpect_hcal

#include <vector>
#include <iostream>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

const Double_t Mn = 0.939565413;
const Double_t Mp = 0.938272081;
const Double_t PI = 3.14159;

//main
void calculateXYexp(Double_t sbstheta_rad, //hcal arm theta
		    Double_t hcaldist,      //hcal distance
		    Double_t bbthr,         //electron arm theta
		    Double_t ebeam,         //beam energy
		    Double_t l_tgt,         //target length
		    Double_t ld2tarrho,     //rho for LD2, g/cc
		    Double_t ld2dEdx,       //LD2 energy loss per unit length GeV/cm
		    Double_t ld2uwallthick, //LD2 beer can upstream wall thickness in cm
		    Double_t rho_al,        //rho for aluminum, g/cc
		    Double_t aldEdx,        //Aluminum energy loss per unit length GeV/cm
		    Double_t celldiameter,  //Diameter of target cell, cm
		    Double_t hcal_v_offset, //HCal vertical offset
		    Double_t vz,            //Event vertex position
		    Double_t px,            //Event e' track momentum, x component
		    Double_t py,            //Event e' track momentum, y component
		    Double_t pz,            //Event e' track momentum, z component
		    Double_t p){            //Event e' track momentum magnitude

  //Establish HCal Coordinates
  vector<TVector3> hcal_axes;
  TVector3 hcal_zaxis( sin(-sbstheta_rad), 0, cos(-sbstheta_rad) ); // Clock-wise rotation about Y axis
  TVector3 hcal_xaxis( 0, -1, 0 );                                  // -Y axis of Hall CoS = X axis of hcal CoS
  TVector3 hcal_yaxis = hcal_zaxis.Cross(hcal_xaxis).Unit();
  hcal_axes.push_back(hcal_xaxis);
  hcal_axes.push_back(hcal_yaxis);
  hcal_axes.push_back(hcal_zaxis);

  TVector3 hcal_origin = hcaldist*hcal_axes[2] + hcal_v_offset*hcal_axes[0];

  //Correct the beam energy with target information
  Double_t ebeam_c = ebeam - ( (vz+l_tgt/2.0) * ld2tarrho * ld2dEdx + ld2uwallthick * rho_al * aldEdx );

  //Get the vertex position from tracking variables
  TVector3 vertex( 0., 0., vz );

  //Compute the momentum of the beam electrons
  TLorentzVector pbeam( 0., 0., ebeam_c, ebeam_c );

  //Reconstruct e' momentum with estimated energy loss from target and get e' theta/phi
  Double_t Eloss_outgoing = celldiameter/2.0/sin(bbthr) * ld2tarrho * ld2dEdx;
  Double_t precon = p + Eloss_outgoing;
  TLorentzVector pe( precon*px/p, precon*py/p, precon*pz/p, precon );
  Double_t etheta = acos( pe.Pz() / pe.E() ); 
  Double_t ephi = atan2( pe.Py(), pe.Px() );

  //could put etheta back into Eloss_outgoing to refine the correction

  //Reconstruct the central momentum of scattered nucleons from etheta
  Double_t Nmass = 0.5*( Mn + Mp );
  Double_t pcent = ebeam/( 1. + ( ebeam/Nmass )*( 1.0 - cos(etheta) ) );

  //Calculate nucleon physical quantities from elastic kinematics
  //Double_t nu = pbeam.E() - pcent; //CHECK THIS
  Double_t nu = pbeam.E() - precon;
  Double_t pNexp = sqrt( pow(nu, 2.) + 2. * Nmass * nu );
  Double_t thNexp = acos( (pbeam.E()-pcent*cos(etheta)) / pNexp );
  Double_t phNexp = ephi + PI;

  //Use them to get the unit vector for an elastic nucleon pointing from the target to HCal
  TVector3 pNhat( sin(thNexp) * cos(phNexp), sin(thNexp) * sin(phNexp), cos(thNexp) );

  //Build the two component vector containing expected x and y HCal positions
  vector<Double_t> xyhcalexpect;

  // intersection of a ray with a plane in hcal coordinates
  Double_t sintersect = ( hcal_origin - vertex).Dot(hcal_axes[2] ) / ( pNhat.Dot(hcal_axes[2]) );
  // ray from Hall origin onto the face of hcal where the nucleon hit
  TVector3 hcal_intersect = vertex + sintersect*pNhat; 

  //Calculate the quantities
  Double_t xexpect_hcal = ( hcal_intersect - hcal_origin ).Dot( hcal_axes[0] );
  Double_t yexpect_hcal = ( hcal_intersect - hcal_origin ).Dot( hcal_axes[1] );

  //Push them back into the vector
  xyhcalexpect.push_back( xexpect_hcal );
  xyhcalexpect.push_back( yexpect_hcal );

  //console out
  std::cout << "HCal X expected : HCal Y expected -> " << xyhcalexpect[0] << ":" << xyhcalexpect[1] << std::endl;

}
