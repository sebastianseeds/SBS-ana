#ifndef CUT_H
#define CUT_H

#include <iostream>
#include "TString.h"
#include "econst.h"

namespace cut {

  //// HCal
  // Active area cut
  // MC DB: exblkN_x = exblkN_y = 1 (default)
  std::vector<Double_t> hcalaa_mc (int exblkN_x ,  // No. of blocks to exclude form top and bottom (Default=1) 
				 int exblkN_y);  // No. of blocks to exclude form left and right (Default=1)

  // Data DB pass 0: exblkN_x = exblkN_y = 1 (default); uses data fit to sbs 8 hcal x y data to establish dispersive bounds
  std::vector<Double_t> hcalaa_data (int exblkN_x ,  // No. of blocks to exclude form top and bottom (Default=1) 
				       int exblkN_y);  // No. of blocks to exclude form left and right (Default=1)

  // Data DB pass 1: exblkN_x = exblkN_y = 1 (default)
  std::vector<Double_t> hcalaa_data_alt (int exblkN_x,   // No. of blocks to exclude form top and bottom (Default=1)
				   int exblkN_y);  // No. of blocks to exclude form left and right (Default=1)

  // Check position by event and verify on hcal active area
  bool hcalaaON (Double_t hcalx,                 // vertical (x) pos of hcal detection
		 Double_t hcaly,                 // horizontal (y) pos of hcal detection
		 std::vector<Double_t> hcalaa);  // HCAL active area coordinates

  // Defines fiducial cut region
  std::vector<Double_t> hcalfid (Double_t dxsig,                      // vertical (x) pos of hcal detection
				 Double_t dysig,                      // horizontal (y) pos of hcal detection
				 std::vector<Double_t> hcalaa,        // HCAL active area coordinates
				 Double_t Nsig);                         // Number of sigma of each of dx/dy to cut on
  // Overloaded for most cases
  std::vector<Double_t> hcalfid (Double_t dxsig_p,                    // vertical (x) pos of exp recoil p at HCAL
				 Double_t dxsig_n,                    // vertical (x) pos of exp recoil n at HCAL
				 Double_t dysig,                      // horizontal (y) pos of hcal detection
				 std::vector<Double_t> hcalaa,        // HCAL active area coordinates
				 Double_t Nsig);                         // Number of sigma of each of dx/dy to cut on
  
  std::vector<Double_t> hcalfid_both (Double_t dxsig_p,                    // vertical (x) pos of exp recoil p at HCAL
				      Double_t dysig,                      // horizontal (y) pos of hcal detection
				      std::vector<Double_t> hcalaa,        // HCAL active area coordinates
				      Double_t Nsigx,                      // Number of sigma dx to cut on
				      Double_t Nsigy);                     // Number of sigma dy to cut on

  // Check position by event and verify in hcal fiducial area
  bool hcalfidIN (Double_t hcalx_exp,                                 // m, expected vert. (x) pos of hcal detection
		  Double_t hcaly_exp,                                 // m, expected horiz. (y) pos of hcal detection
		  Double_t dx_pn,                                     // m, obs dx peak diff p/n
		  std::vector<Double_t> fid);                         // coordinates of fiducial cut

  // Check position by event and verify in hcal fiducial area. Return bool indicating how the fid cut failed
  Int_t hcalfidIN_index (Double_t hcalx_exp,                                 // m, expected vert. (x) pos of hcal detection
			 Double_t hcaly_exp,                                 // m, expected horiz. (y) pos of hcal detection
			 Double_t dx_pn,                                     // m, obs dx peak diff p/n
			 std::vector<Double_t> fid);                         // coordinates of fiducial cut


  std::pair<Double_t, Double_t> findFidFailure(Double_t dxsig,                // sigma in x of proton peak
					       Double_t dysig,                // sigma in y of proton peak
					       Double_t hcalx_exp,            // vertical (x) exp pos of hcal detection
					       Double_t hcaly_exp,            // horizontal (y) exp pos of hcal detection
					       Double_t dx_pn,                // m, obs dx peak diff p/n
					       std::vector<Double_t> hcalaa);  // HCAL active area coordinates

  ////BBCal

  ////Hodo

  ////GEM

  ////Grinch

}

#endif
