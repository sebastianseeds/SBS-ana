//sseeds 4.21.23 - Short script to generate db params with vertical and horizontal offsets adapted to work without sbs-ana framework

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

const Int_t kNrows = 24;                 // Total number of HCal rows
const Int_t kNcols = 12;                 // Total number of HCal columns
const Double_t hcalblk_div_h = 0.15494;  //m, horizontal center-to-center dist.
const Double_t hcalblk_div_v = 0.15875;  //m, vertical center-to-center dist.
const Double_t hcal_hrange = 1.85928;    //m, total range in horizontal direction of HCal (end-to-end)
const Double_t hcal_vrange = 3.81;       //m, total range in vertical direction of HCal (end-to-end)

// Get today's date
string getDate(){
  time_t now = time(0);
  tm ltm = *localtime(&now);
  
  string yyyy = to_string(1900 + ltm.tm_year);
  string mm = to_string(1 + ltm.tm_mon);
  string dd = to_string(ltm.tm_mday);
  string date = mm + '_' + dd + '_' + yyyy;
  
  return date;
}

//v_offset should be positive for towards the sky, h_offset should be positive for away from beamline
//v_offset for vertical offset (x, dispersive)
//h_offset for horizontal offset (y, non-dispersive) 
void HCal_genDBoffsets( Double_t v_offset = -1000., Double_t h_offset = -1000.){ //main

  //Require that user pass offsets
  if( v_offset==-1000. || h_offset==-1000. ){
    std::cout << "Error: User must pass two args <vertical offset (+,sky)> <horizontal offset (+,away-from-beam)>" << std::endl;
    return;
  }

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  // Get the date
  string date = getDate();

  // Create files to store hcal positions
  ofstream hcal_db_offsets;
  std::string offPath = "hcal_db_offset.txt";

  hcal_db_offsets.open( offPath );
  
  // Declare arrays for storing calculated positions
  Double_t vpos[kNrows];
  Double_t hpos[kNcols];

  //Calculate positions for blocks. Each position in db refers to the center of the block
  for( Int_t r=0; r<kNrows; r++ ){ //loop over rows
    vpos[r] = -( ( hcal_vrange - hcalblk_div_v ) / 2 ) - v_offset + ( r * hcalblk_div_v );
  }

  for( Int_t c=0; c<kNcols; c++ ){ //loop over cols
    hpos[c] = -( ( hcal_hrange - hcalblk_div_h ) / 2 ) - h_offset + ( c * hcalblk_div_h );
  }
  
  Double_t test = -(( 3.81 - 0.15875 ) / 2 ) - 0.365;

  //Console/txt outs  
  std::cout << std::setprecision(7); //std::cout default precision is 5 sig figs, set to 7
  hcal_db_offsets << std::setprecision(7); //std::cout default precision is 5 sig figs, set to 7
  std::cout << std::endl;
  
  hcal_db_offsets << "#HCal position offsets obtained " << date.c_str() << std::endl;
  hcal_db_offsets << "sbs.hcal.ypos =" << std::endl;
  std::cout << "sbs.hcal.ypos =" << std::endl;

  for( Int_t r=0; r<kNrows; r++ ){ //loop over rows
    
    for( Int_t c=0; c<kNcols; c++ ){ //loop over cols

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

  for( Int_t r=0; r<kNrows; r++ ){ //loop over rows
    
    for( Int_t c=0; c<kNcols; c++ ){ //loop over cols

      hcal_db_offsets << vpos[r] << "  ";
      std::cout << vpos[r] << "  ";

    }
    hcal_db_offsets << std::endl;
    std::cout << std::endl;
  }
  
  hcal_db_offsets.close();

  std::cout << std::endl << "Done. CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << std::endl;

  std::cout << "Database parameters generated and stored in working directory: " << offPath << std::endl;

}//end main
