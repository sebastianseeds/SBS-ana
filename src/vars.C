#include "../include/vars.h"

namespace vars {
  
  //e' momentum reconstructed w/track angles and beam E
  Double_t pcentral( Double_t ebeam, Double_t etheta, std::string Ntype ) {
    Double_t temp = 0.;
    if( Ntype.compare("p") == 0 ) 
      temp = ebeam/( 1. + ( ebeam/physconst::Mp )*( 1.0 - cos(etheta) ) );
    else if( Ntype.compare("n") == 0 ) 
      temp = ebeam/( 1. + ( ebeam/physconst::Mn )*( 1.0 - cos(etheta) ) );
    else if( Ntype.compare("np") == 0 ) {
      Double_t Nmass = 0.5*( physconst::Mn + physconst::Mp );
      temp = ebeam/( 1. + ( ebeam/Nmass )*( 1.0 - cos(etheta) ) );
    }
    else
      std::cerr << "Error: [vars::pcentral] cannot be calculated. Enter a valid nucleon type." << std::endl;
    return temp;
  }

  //e' theta wrt exit beamline
  Double_t etheta( TLorentzVector Peprime ) {
    return acos( Peprime.Pz() / Peprime.E() );
  }

  //e' phi wrt exit beamline
  Double_t ephi( TLorentzVector Peprime ) {
    return atan2( Peprime.Py(), Peprime.Px() );
  }

  //rest nucleon four-vector (hall coordinates)
  void setPN( std::string Ntype, TLorentzVector &PN ) {
    if(Ntype.compare( "p" ) == 0) 
      PN.SetPxPyPzE( 0., 0., 0., physconst::Mp );
    else if(Ntype.compare( "n" ) == 0) 
      PN.SetPxPyPzE( 0., 0., 0., physconst::Mn );
    else if(Ntype.compare( "np" ) == 0) 
      PN.SetPxPyPzE( 0., 0., 0., 0.5*(physconst::Mn+physconst::Mp) );
    else
      std::cerr << "Error: [vars::setPN] cannot be calculated. Enter a valid nucleon type." << std::endl;
  } 

  //scattered nucleon expected momentum
  Double_t pN_expect( Double_t nu, std::string Ntype ) {
    if(Ntype.compare( "p" ) == 0)                      
      return sqrt( pow(nu, 2.) + 2. * physconst::Mp * nu );
    else if( Ntype.compare( "n" ) == 0 )      
      return sqrt( pow(nu, 2.) + 2. * physconst::Mn * nu );
    else if( Ntype.compare( "np" ) == 0 )      
      return sqrt( pow(nu, 2.) + 2. * 0.5*(physconst::Mn+physconst::Mp) * nu );
    else {
      std::cerr << "Error: [vars::pN_expect] cannot be calculated. Enter a valid nucleon type." << std::endl;
      return -1;
    }
  }

  //virtual photon q unit vector
  TVector3 qVect_unit( Double_t Ntheta, Double_t Nphi ) {
    TVector3 pNhat( sin(Ntheta) * cos(Nphi), sin(Ntheta) * sin(Nphi), cos(Ntheta) );
    return pNhat;
  }

  //sets hcal axes by kinematic (from hall coordinates)
  void sethcalaxes( Double_t sbstheta_rad,                           // SBS angle in radian 
		    vector<TVector3> &hcal_axes ) {
    TVector3 hcal_zaxis( sin(-sbstheta_rad), 0, cos(-sbstheta_rad) ); // Clock-wise rotation about Y axis
    TVector3 hcal_xaxis( 0, -1, 0 );                                  // -Y axis of Hall CoS = X axis of hcal CoS
    TVector3 hcal_yaxis = hcal_zaxis.Cross(hcal_xaxis).Unit();
    hcal_axes.push_back(hcal_xaxis);
    hcal_axes.push_back(hcal_yaxis);
    hcal_axes.push_back(hcal_zaxis);
  }

  //gets expected location of scattered nucleon assuming straight line projections from BB track
  void getxyhcalexpect( TVector3 vertex, TVector3 pNhat, TVector3 hcal_origin, 
			vector<TVector3> hcal_axes, vector<Double_t> &xyhcalexpect ) {
    // intersection of a ray with a plane in hcal coordinates
    Double_t sintersect = ( hcal_origin - vertex).Dot(hcal_axes[2] ) / ( pNhat.Dot(hcal_axes[2]) );
    // ray from Hall origin onto the face of hcal where the nucleon hit
    TVector3 hcal_intersect = vertex + sintersect*pNhat; 

    Double_t xexpect_hcal = ( hcal_intersect - hcal_origin ).Dot( hcal_axes[0] );
    Double_t yexpect_hcal = ( hcal_intersect - hcal_origin ).Dot( hcal_axes[1] );

    xyhcalexpect.push_back( xexpect_hcal );
    xyhcalexpect.push_back( yexpect_hcal );
  }

  //calculates inverse momentum transfer Q^2
  Double_t Q2( Double_t ebeam, Double_t eeprime, Double_t etheta ) {
    return 2.0 * ebeam * eeprime*( 1.0 - cos(etheta) );
  }

  //calculates invariant mass squared W^2
  Double_t W2( Double_t ebeam, Double_t eeprime, Double_t Q2, std::string Ntype ) {
    Double_t temp = 0.;
    if( Ntype.compare("p") == 0 ) 
      temp = pow( physconst::Mp, 2.0 ) + 2.0*physconst::Mp * (ebeam-eeprime) - Q2;
    else if( Ntype.compare( "n" ) == 0 ) 
      temp = pow( physconst::Mn, 2.0 ) + 2.0*physconst::Mn * (ebeam-eeprime) - Q2;
    else if( Ntype.compare( "np" ) == 0 ) 
      temp = pow( 0.5*(physconst::Mn + physconst::Mp), 2.0 ) + 2.0 * 0.5*(physconst::Mn + physconst::Mp) * (ebeam-eeprime) - Q2;
    else
      std::cerr << "Error: [vars::W2] cannot be calculated. Enter a valid nucleon type." << std::endl;
    return temp;
  }

  //calculates invariant mass W
  Double_t W( Double_t ebeam, Double_t eeprime, Double_t Q2, std::string Ntype ) {
    return max( 0., sqrt( vars::W2( ebeam, eeprime, Q2, Ntype ) ) );
  }

  //calculates luminosity
  Double_t luminosity( Double_t ibeam, std::string targetType ) {
    Double_t lumi = 0.;
    if(targetType.compare("LH2") == 0)
      lumi = ( (ibeam/physconst::qe) * econst::l_tgt * econst::lh2tarrho * (physconst::N_A/physconst::H2_Amass) );
    else if(targetType.compare("LD2") == 0)
      lumi = ( (ibeam/physconst::qe) * econst::l_tgt * econst::ld2tarrho * (physconst::N_A/physconst::D2_Amass) );
    else
      std::cerr << "Error: [vars::Luminosity] cannot be calculated. Enter a valid nucleon type." << std::endl;
    return lumi;
  }

  //beam energy loss correction to vertex
  Double_t ebeam_c( Double_t ebeam, Double_t vz, std::string target ){
    Double_t ebeam_c = 0.;
    if( target.compare("LH2") == 0 )
      ebeam_c = ebeam - ( (vz+econst::l_tgt/2.0) * econst::lh2tarrho * econst::lh2dEdx + econst::lh2uwallthick * econst::rho_al * econst::aldEdx );
    else if( target.compare("LD2") == 0 )
      ebeam_c = ebeam - ( (vz+econst::l_tgt/2.0) * econst::ld2tarrho * econst::ld2dEdx + econst::ld2uwallthick * econst::rho_al * econst::aldEdx );
    else
      std::cerr << "Warning: [vars::ebeam_c] cannot be calculated. Enter a valid nucleon type." << std::endl;
    return ebeam;
  }

  //beam energy loss correction from vertex, estimate for now
  Double_t ebeam_o( Double_t ebeam_c, Double_t etheta, std::string target ){
    Double_t ebeam_o = 0.;
    if( target.compare("LH2") == 0 )
      ebeam_o = ebeam_c - econst::celldiameter/2.0/sin(etheta) * econst::lh2tarrho * econst::lh2dEdx;
    else if( target.compare("LD2") == 0 )
      ebeam_o = ebeam_c - econst::celldiameter/2.0/sin(etheta) * econst::ld2tarrho * econst::ld2dEdx;
    else
      std::cerr << "Warning: [vars::ebeam_o] cannot be calculated. Enter a valid nucleon type." << std::endl;
    return ebeam_c;
  }

  //beam electron four-vector (hall coordinates)
  void setpbeam( Double_t ebeam, TLorentzVector &pbeam ){
    pbeam.SetPxPyPzE( 0., 0., ebeam, ebeam );
  } 

} //::vars
