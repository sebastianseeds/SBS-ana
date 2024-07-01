#include "../include/cut.h"

namespace cut {

  // Establish hcal active area excluding N blks from edge, MC DB
  std::vector<Double_t> hcalaa_mc (int exblkN_x=1, int exblkN_y=1) {
    std::vector<Double_t> hcalaa;
    Double_t hcalaaXi = econst::hcalposXi_mc + exblkN_x*econst::hcalblk_div_v;
    Double_t hcalaaXf = econst::hcalposXf_mc - exblkN_x*econst::hcalblk_div_v;
    Double_t hcalaaYi = econst::hcalposYi_mc + exblkN_y*econst::hcalblk_div_h;
    Double_t hcalaaYf = econst::hcalposYf_mc - exblkN_y*econst::hcalblk_div_h;
    hcalaa.push_back( hcalaaXi ); 
    hcalaa.push_back( hcalaaXf );
    hcalaa.push_back( hcalaaYi );
    hcalaa.push_back( hcalaaYf );
    return hcalaa;
  }

  // Establish hcal active area excluding N blks from edge, >pass0 DB
  std::vector<Double_t> hcalaa_data (int exblkN_x=1, int exblkN_y=1) {
    std::vector<Double_t> hcalaa;
    Double_t hcalaaXi = econst::hcalposXi_mc + exblkN_x*econst::hcalblk_div_v;
    Double_t hcalaaXf = econst::hcalposXf_mc - exblkN_x*econst::hcalblk_div_v;
    Double_t hcalaaYi = econst::hcalposYi_mc + exblkN_y*econst::hcalblk_div_h;
    Double_t hcalaaYf = econst::hcalposYf_mc - exblkN_y*econst::hcalblk_div_h;
    hcalaa.push_back( hcalaaXi ); 
    hcalaa.push_back( hcalaaXf );
    hcalaa.push_back( hcalaaYi );
    hcalaa.push_back( hcalaaYf );
    return hcalaa;
  }

  // Establish hcal active area excluding N blks from edge, pass 0 DB
  std::vector<Double_t> hcalaa_data_alt (int exblkN_x=1, int exblkN_y=1) {
    std::vector<Double_t> hcalaa;
    Double_t hcalaaXi = econst::hcalposXi_p0 + exblkN_x*econst::hcalblk_h_p0;
    Double_t hcalaaXf = econst::hcalposXf_p0 - exblkN_x*econst::hcalblk_h_p0;
    Double_t hcalaaYi = econst::hcalposYi_p0 + exblkN_y*econst::hcalblk_w_p0;
    Double_t hcalaaYf = econst::hcalposYf_p0 - exblkN_y*econst::hcalblk_w_p0;
    hcalaa.push_back( hcalaaXi ); 
    hcalaa.push_back( hcalaaXf );
    hcalaa.push_back( hcalaaYi );
    hcalaa.push_back( hcalaaYf );
    return hcalaa;
  }

  // Check position per event against hcal active area (TRUE if detection on active area)
  bool hcalaaON (Double_t hcalx, Double_t hcaly, std::vector<Double_t> hcalaa) {
    bool on = false;
    // active area dimensions
    Double_t hcalx_t = hcalaa[0];
    Double_t hcalx_b = hcalaa[1];
    Double_t hcaly_l = hcalaa[2];
    Double_t hcaly_r = hcalaa[3];
    on = hcaly>hcaly_l && hcaly<hcaly_r && hcalx>hcalx_t && hcalx<hcalx_b;
    return on;
  } 

  // Establish coordinates of fiducial cut with x/y margin
  std::vector<Double_t> hcalfid (Double_t dxsig, Double_t dysig, std::vector<Double_t> hcalaa, Double_t Nsig=1.) {
    std::vector<Double_t> fid;
    Double_t hcalx_t = hcalaa[0] + Nsig*dxsig;  // -X margin
    Double_t hcalx_b = hcalaa[1] - Nsig*dxsig;  // +X margin
    Double_t hcaly_l = hcalaa[2] + Nsig*dysig;  // -Y margin
    Double_t hcaly_r = hcalaa[3] - Nsig*dysig;  // +Y margin
    fid.push_back( hcalx_t ); 
    fid.push_back( hcalx_b );
    fid.push_back( hcaly_l );
    fid.push_back( hcaly_r );
    return fid;
  }

  // Overload for most cases where dy margin negligable and dx top/bottom configurable by p/n
  std::vector<Double_t> hcalfid (Double_t dxsig_p, Double_t dxsig_n, Double_t dysig, std::vector<Double_t> hcalaa, Double_t Nsig=1.) {
    std::vector<Double_t> fid;
    Double_t hcalx_t = hcalaa[0] + Nsig*dxsig_p;  // -X margin (relevant for proton)
    Double_t hcalx_b = hcalaa[1] - Nsig*dxsig_n;  // +X margin (relevant for neutron)
    Double_t hcaly_l = hcalaa[2] + Nsig*dysig;    // -Y margin
    Double_t hcaly_r = hcalaa[3] - Nsig*dysig;    // +Y margin
    fid.push_back( hcalx_t ); 
    fid.push_back( hcalx_b );
    fid.push_back( hcaly_l );
    fid.push_back( hcaly_r );
    return fid;
  }

  // Overload for case where nsig varies between disp and trans directions
  std::vector<Double_t> hcalfid_both (Double_t dxsig, Double_t dysig, std::vector<Double_t> hcalaa, Double_t Nsigx=1., Double_t Nsigy=1.) {
    std::vector<Double_t> fid;
    Double_t hcalx_t = hcalaa[0] + Nsigx*dxsig;  // -X margin (relevant for proton)
    Double_t hcalx_b = hcalaa[1] - Nsigx*dxsig;  // +X margin (relevant for neutron)
    Double_t hcaly_l = hcalaa[2] + Nsigy*dysig;  // -Y
    Double_t hcaly_r = hcalaa[3] - Nsigy*dysig;  // +Y
    fid.push_back( hcalx_t ); 
    fid.push_back( hcalx_b );
    fid.push_back( hcaly_l );
    fid.push_back( hcaly_r );
    return fid;
  }

  // Check position by event and verify in hcal fiducial area
  bool hcalfidIN (Double_t hcalx_exp, Double_t hcaly_exp, Double_t dx_pn, vector<Double_t> fid) {
    Double_t hcalx_t = fid[0]; //hcal -X fid
    Double_t hcalx_b = fid[1]; //hcal +X fid
    Double_t hcaly_l = fid[2]; //hcal -Y fid
    Double_t hcaly_r = fid[3]; //hcal +Y fid

    //proton hypothesis
    Double_t hcalx_exp_p = hcalx_exp - dx_pn;      //define the exp pos of a proton from obs dx peak diff

    bool infid = (hcaly_exp>hcaly_l) && (hcaly_exp<hcaly_r) &&      //dy same for protons and neutrons
		 (hcalx_exp>hcalx_t) && (hcalx_exp<hcalx_b) &&      //dx acceptance check for neutrons
                 (hcalx_exp_p>hcalx_t) && (hcalx_exp_p<hcalx_b);	//dx acceptance check for for protons							     
    return infid;
  } 

  // Check position by event and verify in hcal fiducial area. 
  // If pass both dispersive and transverse fiducial cuts = 0
  // If pass only transverse cut = 1
  // If pass only dispersive cut = 2
  // If fail both = 3
  Int_t hcalFidIN_index(Double_t hcalx_exp, Double_t hcaly_exp, Double_t dx_pn, vector<Double_t> fid) {
    Double_t hcalx_t = fid[0]; // hcal -X fid
    Double_t hcalx_b = fid[1]; // hcal +X fid
    Double_t hcaly_l = fid[2]; // hcal -Y fid
    Double_t hcaly_r = fid[3]; // hcal +Y fid

    // proton hypothesis
    Double_t hcalx_exp_p = hcalx_exp - dx_pn; // define the exp pos of a proton from obs dx peak diff

    // Check if hcalx_exp and hcaly_exp are within the fiducial boundaries
    bool infid_disp = (hcalx_exp_p > hcalx_t) && (hcalx_exp_p < hcalx_b) &&    //proton cut
                      (hcalx_exp > hcalx_t) && (hcalx_exp < hcalx_b);          //neutron cut

    bool infid_trans = (hcaly_exp > hcaly_l) && (hcaly_exp < hcaly_r);

    // Determine the index based on which cuts are passed
    if (infid_disp && infid_trans) {
      return 0; // Passes both cuts
    } else if (!infid_disp && infid_trans) {
      return 1; // Only passes transverse cut
    } else if (infid_disp && !infid_trans) {
      return 2; // Only passes dispersive cut
    } else {
      return 3; // Fails both cuts
    }
  }

  //Calculates the number of sigma in x and y after which the event fails the fiducial cut
  std::pair<Double_t, Double_t> findFidFailure(Double_t dxsig, 
					       Double_t dysig, 
					       Double_t hcalx_exp, 
					       Double_t hcaly_exp, 
					       Double_t dx_pn,
					       std::vector<Double_t> hcalaa) {

    //Calulate both hypotheses
    double hcalx_exp_n = hcalx_exp; //neutron
    double hcalx_exp_p = hcalx_exp - dx_pn; //proton

    //cout << "expected neutron x: " << hcalx_exp_n << ", expected proton x: " << hcalx_exp_p << endl;

    //Get active area
    Double_t hcalx_t = hcalaa[0]; //hcal -X
    Double_t hcalx_b = hcalaa[1]; //hcal +X 
    Double_t hcaly_l = hcalaa[2]; //hcal -Y 
    Double_t hcaly_r = hcalaa[3]; //hcal +Y 

    //cout << "hcal top boundary: " << hcalx_t << ", hcal bottom boundary: " << hcalx_b << endl;

    //Calculate factor of sigma necessary to cause fiducial cut to fail. Leave 0 where hcalaa fail
    Double_t Nsigy_l = std::max(0.0, (hcaly_exp - hcaly_l) / dysig);
    Double_t Nsigy_r = std::max(0.0, (hcaly_r - hcaly_exp) / dysig);
    Double_t Nsigx_t_n = std::max(0.0, (hcalx_exp_n - hcalx_t) / dxsig);
    Double_t Nsigx_b_n = std::max(0.0, (hcalx_b - hcalx_exp_n) / dxsig);
    Double_t Nsigx_t_p = std::max(0.0, (hcalx_exp_p - hcalx_t) / dxsig);
    Double_t Nsigx_b_p = std::max(0.0, (hcalx_b - hcalx_exp_p) / dxsig);
    
    //Get min sigma to fail (both hypothesis in x)
    
    //Return zero sigma where projection is off of hcal aa in y
    Double_t Nsigy = 0.;
    if( Nsigy_l!=0. && Nsigy_r!=0. )
      Nsigy = std::min( Nsigy_l, Nsigy_r);

    //Return zero sigma where projection is off of hcal aa in x
    Double_t Nsigx_n = 0.;
    if( Nsigx_b_n!=0. && Nsigx_t_n!=0. )
      Nsigx_n = std::min( Nsigx_b_n, Nsigx_t_n);
    Double_t Nsigx_p = 0.;
    if( Nsigx_b_p!=0. && Nsigx_t_p!=0. )
      Nsigx_p = std::min( Nsigx_b_p, Nsigx_t_p);

    //Get hypothesis which fails first in x
    Double_t Nsigx = 0.;
    if( Nsigx_n!=0. && Nsigx_p!=0. )
      Nsigx = std::min( Nsigx_n, Nsigx_p);

    // Return the Nsigx and Nsigy values that caused the failure
    return std::make_pair(Nsigx, Nsigy);
  }


}
