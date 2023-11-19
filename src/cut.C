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

  // Establish hcal active area excluding N blks from edge, Pass0/1 DB
  std::vector<Double_t> hcalaa_data (int exblkN_x=1, int exblkN_y=1) {
    std::vector<Double_t> hcalaa;
    Double_t hcalaaXi = econst::hcalposXi_p0 + exblkN_x*econst::hcalblk_w_p0;
    Double_t hcalaaXf = econst::hcalposXf_p0 - exblkN_x*econst::hcalblk_w_p0;
    Double_t hcalaaYi = econst::hcalposYi_p0 + exblkN_y*econst::hcalblk_h_p0;
    Double_t hcalaaYf = econst::hcalposYf_p0 - exblkN_y*econst::hcalblk_h_p0;
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
    Double_t hcaly_r = hcalaa[2];
    Double_t hcaly_l = hcalaa[3];
    on = hcaly>hcaly_r && hcaly<hcaly_l && hcalx>hcalx_t && hcalx<hcalx_b;
    return on;
  } 

  // Establish coordinates of fiducial cut with x/y margin
  std::vector<Double_t> hcalfid (Double_t dxsig, Double_t dysig, std::vector<Double_t> hcalaa) {
    std::vector<Double_t> fid;
    Double_t hcalx_t = hcalaa[0] + dxsig;  // top margin
    Double_t hcalx_b = hcalaa[1] - dxsig;  // bottom margin
    Double_t hcaly_r = hcalaa[2] + dysig;  // right margin
    Double_t hcaly_l = hcalaa[3] - dysig;  // left margin
    fid.push_back( hcalx_t ); 
    fid.push_back( hcalx_b );
    fid.push_back( hcaly_r );
    fid.push_back( hcaly_l );
    return fid;
  }

  // Overload for most cases where dy margin negligable and dx top/bottom configurable by p/n
  std::vector<Double_t> hcalfid (Double_t dxsig_p, Double_t dxsig_n, Double_t dysig, std::vector<Double_t> hcalaa) {
    std::vector<Double_t> fid;
    Double_t hcalx_t = hcalaa[0] + dxsig_p;  // top margin (relevant for proton)
    Double_t hcalx_b = hcalaa[1] - dxsig_n;  // bottom margin (relevant for neutron)
    Double_t hcaly_r = hcalaa[2] + dysig;    // right margin
    Double_t hcaly_l = hcalaa[3] - dysig;    // left margin
    fid.push_back( hcalx_t ); 
    fid.push_back( hcalx_b );
    fid.push_back( hcaly_r );
    fid.push_back( hcaly_l );
    return fid;
  }

  // Check position by event and verify in hcal fiducial area
  bool hcalfidIN (Double_t hcalx_exp, Double_t hcaly_exp, Double_t dx_pn, vector<Double_t> fid) {
    Double_t hcalx_t = fid[0];
    Double_t hcalx_b = fid[1];
    Double_t hcaly_r = fid[2];
    Double_t hcaly_l = fid[3];

    Double_t hcalx_exp_p = hcalx_exp - dx_pn;      //define the exp pos of a proton from obs dx peak diff

    bool infid = hcaly_exp>hcaly_r && hcaly_exp<hcaly_l &&      //dy same for protons and neutrons
	         hcalx_exp>hcalx_t && hcalx_exp<hcalx_b &&      //dx for neutrons
	         hcalx_exp_p>hcalx_t && hcalx_exp_p<hcalx_b;	//dx for protons							     
    return infid;
  } 

}
