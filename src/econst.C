#include "../include/econst.h"

namespace econst {

  //Wide cuts using all subsystems
  std::string globcut(Int_t config) {
    if(config==1)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&abs(bb.etot_over_p-0.92)<0.2&&sbs.hcal.e>0.01&&bb.ps.e+bb.sh.e>1.7";
    else if(config==4)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&abs(bb.etot_over_p-0.92)<0.2&&sbs.hcal.e>0.01&&bb.ps.e+bb.sh.e>1.7";
    else if(config==7)
      return "bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>2&&bb.tr.p[0]>2.0&&sbs.hcal.e>0.01&&bb.ps.e>0.2";
    else if(config==11)
      return "bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>2&&bb.tr.p[0]>2.0&&sbs.hcal.e>0.01&&bb.ps.e>0.2";
    else if(config==14)
      return "bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>2&&bb.tr.p[0]>1.6&&sbs.hcal.e>0.01&&bb.ps.e>0.2";
    else if(config==8)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&sbs.hcal.e>0.01&&abs(bb.tr.tg_th[0])<0.15&&abs(bb.tr.tg_ph[0])<0.3"; 
    else if(config==9)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&sbs.hcal.e>0.01&&abs(bb.tr.tg_th[0])<0.15&&abs(bb.tr.tg_ph[0])<0.3";
    else if(config==4363)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&abs(bb.etot_over_p-0.92)<0.2&&sbs.hcal.e>0.01&&bb.ps.e+bb.sh.e>1.7";
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return "";
    }
  }

  //Wide cuts using Bigbite subsystems
  std::string globcut_earm(Int_t config) {
    if(config==1)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&abs(bb.etot_over_p-0.92)<0.2&&bb.ps.e+bb.sh.e>1.7";
    else if(config==4)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&abs(bb.etot_over_p-0.92)<0.2&&bb.ps.e+bb.sh.e>1.7";
    else if(config==7)
      return "bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>2&&bb.tr.p[0]>2.0&&bb.ps.e>0.2";
    else if(config==11)
      return "bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>2&&bb.tr.p[0]>2.0&&bb.ps.e>0.2";
    else if(config==14)
      return "bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>2&&bb.tr.p[0]>1.6&&bb.ps.e>0.2";
    else if(config==8)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&abs(bb.tr.tg_th[0])<0.15&&abs(bb.tr.tg_ph[0])<0.3"; 
    else if(config==9)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&abs(bb.tr.tg_th[0])<0.15&&abs(bb.tr.tg_ph[0])<0.3";
    else if(config==4363)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&abs(bb.etot_over_p-0.92)<0.2&&bb.ps.e+bb.sh.e>1.7";
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return "";
    }
  }

  //Beam energy in GeV
  Double_t ebeam(Int_t config) {
    // first: lh2, second: ld2
    if(config==1)
      return 1.916;
    else if(config==4)
      //return 3.7278;
      return 3.7393;
    else if(config==7)
      //return 7.906;
      return 7.9308;
    else if(config==11)
      //return 9.91;
      return 9.889;
    else if(config==14)
      //return 5.965;
      return 5.9827;
    else if(config==8)
      //return 5.965;
      return 5.9826;
    else if(config==9)
      //return 4.013;
      return 4.0268;
    else if(config==4363)
      return 6.373;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

  //Angle of the BigBite (electron) arm wrt exit beamline in degrees
  Double_t bbtheta(Int_t config){
    if(config==1)
      return 51.0;
    else if(config==4)
      return 36.0;
    else if(config==7)
      return 40.0;
    else if(config==11)
      return 42.0;
    else if(config==14)
      return 46.5;
    else if(config==8)
      return 26.5;
    else if(config==9)
      return 49.0;
    else if(config==4363)
      return 36.5;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

  //Distance from the target to the BigBite magnet in m
  Double_t bbdist(Int_t config){
    if(config==1)
      return 1.8518;
    else if(config==4)
      return 1.7988;
    else if(config==7)
      return 1.84896;
    else if(config==11)
      return 1.55146;
    else if(config==14)
      return 1.84787;
    else if(config==8)
      return 1.97473;
    else if(config==9)
      return 1.550;
    else if(config==4363)
      return 1.63;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

  //Angle of the SBS magnet (in hadron arm) wrt exit beamline in degrees
  Double_t sbstheta(Int_t config){
    if(config==1)
      return 33.5;
    else if(config==4)
      return 31.9;
    else if(config==7)
      return 16.1;
    else if(config==11)
      return 13.3;
    else if(config==14)
      return 17.3;
    else if(config==8)
      return 29.9;
    else if(config==9)
      return 22.5;
    else if(config==4363)
      return 22.1;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

  //Distance from target to the SBS magnet in m
  Double_t sbsdist(Int_t config){
    if(config==1||config==4||config==7||config==11
       ||config==14||config==8||config==9)
      return 2.25;
    else if(config==4363)
      return 2.8;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

  //Effective distance from target to the hadron calorimeter in m to eliminate dx:hcalx correlation
  Double_t hcaleffdist(Int_t config){
    if(config==1)
      return 13.5;
    else if(config==4)
      return 11.65;
    else if(config==7)
      return 14.925;
    else if(config==11)
      return 15.77;
    else if(config==14)
      return 15.13;
    else if(config==8)
      return 11.77;
    else if(config==9)
      return 11.7;
    else if(config==4363)
      return 17.0;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

  //Distance from target to the hadron calorimeter in m
  Double_t hcaldist(Int_t config){
    if(config==1)
      return 13.5;
    else if(config==4||config==8||config==9)
      return 11.0;
    else if(config==7||config==14)
      return 14.0;
    else if(config==11)
      return 14.5;
    else if(config==4363)
      return 17.0;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

  //Angle hadron calorimeter makes wrt exit beamline in degrees
  Double_t hcaltheta(Int_t config){
    if(config==1)
      return 33.5;
    else if(config==4)
      return 31.9;
    else if(config==7)
      return 16.1;
    else if(config==11)
      return 13.3;
    else if(config==14)
      return 17.3;
    else if(config==8)
      return 29.4;
    else if(config==9)
      return 22.0;
    else if(config==4363)
      return 21.6; 
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

} //::econst
