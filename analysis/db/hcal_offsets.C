//sseeds 4.21.23 - Short script to generate db params with vertical and horizontal offsets adapted from old

#include <vector>
#include <iostream>
#include <iomanip>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

//v_offset should be positive for towards the sky, h_offset should be positive for away from beamline
//v_offset for vertical offset (x, dispersive)
//h_offset for horizontal offset (y, non-dispersive) 
void hcal_offsets( Double_t v_offset, Double_t h_offset){ //main
  
  //Require that user pass offsets
  if( v_offset==-1000. || h_offset==-1000. ){
    std::cout << "Error: User must pass two args <vertical offset (+,sky)> <horizontal offset (+,away-from-beam)>" << std::endl;
    return;
  }

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // Get the date
  string date = util::getDate();

  // Create files to store hcal positions
  ofstream hcal_db_offsets;
  std::string offPath = "outfiles/hcal_db_offset.txt";

  hcal_db_offsets.open( offPath );
  
  // Declare arrays for storing calculated positions
  Double_t vpos[econst::hcalrow];
  Double_t hpos[econst::hcalcol];

  //Calculate positions for blocks. Each position in db refers to the center of the block
  for( Int_t r=0; r<econst::hcalrow; r++ ){ //loop over rows
    vpos[r] = -( ( econst::hcal_vrange - econst::hcalblk_div_v ) / 2 ) - v_offset + ( r * econst::hcalblk_div_v );
  }

  for( Int_t c=0; c<econst::hcalcol; c++ ){ //loop over cols
    hpos[c] = -( ( econst::hcal_hrange - econst::hcalblk_div_h ) / 2 ) - h_offset + ( c * econst::hcalblk_div_h );
  }
  
  Double_t test = -(( 3.81 - 0.15875 ) / 2 ) - 0.365;

  //Console/txt outs  
  std::cout << std::setprecision(7); //std::cout default precision is 5 sig figs, set to 7
  hcal_db_offsets << std::setprecision(7); //std::cout default precision is 5 sig figs, set to 7
  std::cout << std::endl;
  
  hcal_db_offsets << "#HCal position offsets obtained " << date.c_str() << std::endl;
  hcal_db_offsets << "sbs.hcal.ypos =" << std::endl;
  std::cout << "sbs.hcal.ypos =" << std::endl;

  for( Int_t r=0; r<econst::hcalrow; r++ ){ //loop over rows
    
    for( Int_t c=0; c<econst::hcalcol; c++ ){ //loop over cols

      hcal_db_offsets << hpos[c] << "  ";
      std::cout << hpos[c] << "  ";

    }
    hcal_db_offsets << std::endl;
    std::cout << std::endl;
  }

  hcal_db_offsets << std::endl;
  std::cout << std::endl;
  hcal_db_offsets << "sbs.hcal.xpos =" << std::endl;
  std::cout << "sbs.hcal.xpos =" << std::endl;

  for( Int_t r=0; r<econst::hcalrow; r++ ){ //loop over rows
    
    for( Int_t c=0; c<econst::hcalcol; c++ ){ //loop over cols

      hcal_db_offsets << vpos[r] << "  ";
      std::cout << vpos[r] << "  ";

    }
    hcal_db_offsets << std::endl;
    std::cout << std::endl;
  }
  
  hcal_db_offsets.close();

  std::cout << std::endl << "Done. CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << std::endl;

}//end main
