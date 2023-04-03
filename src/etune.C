#include "../include/etune.h"

namespace etune {

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

  //Elastic invariant mass squared peak GeV^2
  Double_t W2mean(Int_t config) {
    if(config==1)
      return 0.92;
    else if(config==4)
      return 0.917994;
    else if(config==7)
      return 0.88;
    else if(config==11)
      return 0.925;
    else if(config==14)
      return 0.870819;
    else if(config==8)
      return 0.91;
    else if(config==9)
      return 0.91;
    else if(config==4363)
      return 0.92;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

  //Elastic invariant mass squared sigma GeV^2
  Double_t W2sig(Int_t config) {
    if(config==1)
      return 0.325;
    else if(config==4)
      return 0.167922;
    else if(config==7)
      return 0.5;
    else if(config==11)
      return 0.325;
    else if(config==14)
      return 0.19443;
    else if(config==8)
      return 0.21;
    else if(config==9)
      return 0.17;
    else if(config==4363)
      return 0.325;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

  //Location of neutron elastic peak in HCal dx (m)
  Double_t dx0_n(Int_t config,Int_t mag) {
    if(config==1)
      return 0.0;
    else if(config==4){
      if(mag==0) return 0.01914;
      else if(mag==30) return -0.00632;
      else if(mag==50) return 0.0129126;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS4 (hint: 0, 30, or 50)." << std::endl;
      return -1;
      }
    }else if(config==7){
      if(mag==85) return 0.02437;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS7 (hint: it's 85)." << std::endl;
      return -1;
      }
    }else if(config==11){
      if(mag==0) return 0.0405545;
      else if(mag==100) return 0.0728046;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS11 (hint: 0 or 100)." << std::endl;
      return -1;
      }
    }else if(config==14){
      if(mag==0) return 0.078859;
      else if(mag==70) return 0.078859;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: 0 or 70)." << std::endl;
      return -1;
      }
    }else if(config==8){
      if(mag==0) return 0.100465;
      else if(mag==50) return 0.0876799;
      else if(mag==70) return 0.0920908;
      else if(mag==100) return 0.0904357;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS8 (hint: 0, 50, 70, or 100)." << std::endl;
      return -1;
      }
    }else if(config==9){
      if(mag==70) return 0.08268;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: it's 70)." << std::endl;
      return -1;
      }
    }else if(config==4363)
      return 0.0;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

  //Location of proton elastic peak in HCal dx (m)
  Double_t dx0_p(Int_t config,Int_t mag) {
  if(config==1)
      return 0.0;
    else if(config==4){
      if(mag==0) return 0.01914;
      else if(mag==30) return -0.6449;
      else if(mag==50) return -1.08927;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS4 (hint: 0, 30, or 50)." << std::endl;
      return -1;
      }
    }else if(config==7){
      if(mag==85) return -0.6914;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS7 (hint: it's 85)." << std::endl;
      return -1;
      }
    }else if(config==11){
      if(mag==0) return 0.0405545;
      else if(mag==100) return -0.633216;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS11 (hint: 0 or 100)." << std::endl;
      return -1;
      }
    }else if(config==14){
      if(mag==0) return 0.078859;
      else if(mag==70) return -0.719273;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: 0 or 70)." << std::endl;
      return -1;
      }
    }else if(config==8){
      if(mag==0) return 0.100465;
      else if(mag==50) return -0.516198;
      else if(mag==70) return -0.765124;
      else if(mag==100) return -1.13518;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS8 (hint: 0, 50, 70, or 100)." << std::endl;
      return -1;
      }
    }else if(config==9){
      if(mag==70) return -0.82026;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: it's 70)." << std::endl;
      return -1;
      }
    }else if(config==4363)
      return -0.5;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

  //Location of elastic peaks in HCal dy (m)
  Double_t dy0(Int_t config) {
      if(config==1)
      return 0.;
    else if(config==4)
      return -0.0270143;
    else if(config==7)
      return 0.0158817;
    else if(config==11)
      return -0.0119;
    else if(config==14)
      return 0.00220912;
    else if(config==8)
      return -0.0119;
    else if(config==9)
      return -0.0175041;
    else if(config==4363)
      return 0.;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

  //Sigma of neutron elastic peak in HCal dx (m)
  Double_t dxsig_n(Int_t config,Int_t mag) {
  if(config==1)
      return 0.0;
    else if(config==4){
      if(mag==0) return 0.05911;
      else if(mag==30) return 0.2128;
      else if(mag==50) return 0.193991;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS4 (hint: 0, 30, or 50)." << std::endl;
      return -1;
      }
    }else if(config==7){
      if(mag==85) return 0.107;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS7 (hint: it's 85)." << std::endl;
      return -1;
      }
    }else if(config==11){
      if(mag==0) return 0.0635993;
      else if(mag==100) return 0.187977;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS11 (hint: 0 or 100)." << std::endl;
      return -1;
      }
    }else if(config==14){
      if(mag==0) return 0.171009;
      else if(mag==70) return 0.171009;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: 0 or 70)." << std::endl;
      return -1;
      }
    }else if(config==8){
      if(mag==0) return 0.163139;
      else if(mag==50) return 0.195652;
      else if(mag==70) return 0.185411;
      else if(mag==100) return 0.179498;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS8 (hint: 0, 50, 70, or 100)." << std::endl;
      return -1;
      }
    }else if(config==9){
      if(mag==70) return 0.173228;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: it's 70)." << std::endl;
      return -1;
      }
    }else if(config==4363)
      return 0.1;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

  //Sigma of proton elastic peak in HCal dx (m)
  Double_t dxsig_p(Int_t config,Int_t mag) {
  if(config==1)
      return 0.0;
    else if(config==4){
      if(mag==0) return 0.05911;
      else if(mag==30) return 0.1806;
      else if(mag==50) return 0.193991;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS4 (hint: 0, 30, or 50)." << std::endl;
      return -1;
      }
    }else if(config==7){
      if(mag==85) return 0.1246;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS7 (hint: it's 85)." << std::endl;
      return -1;
      }
    }else if(config==11){
      if(mag==0) return 0.0635993;
      else if(mag==100) return 0.191482;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS11 (hint: 0 or 100)." << std::endl;
      return -1;
      }
    }else if(config==14){
      if(mag==0) return 0.171009;
      else if(mag==70) return 0.161225;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: 0 or 70)." << std::endl;
      return -1;
      }
    }else if(config==8){
      if(mag==0) return 0.163139;
      else if(mag==50) return 0.186157;
      else if(mag==70) return 0.2023366;
      else if(mag==100) return 0.212879;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS8 (hint: 0, 50, 70, or 100)." << std::endl;
      return -1;
      }
    }else if(config==9){
      if(mag==70) return 0.172697;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: it's 70)." << std::endl;
      return -1;
      }
    }else if(config==4363)
      return 0.1;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

  //Sigma of elastic peaks in HCal dy (m)
  Double_t dysig(Int_t config) {
       if(config==1)
      return 0.1;
    else if(config==4)
      return 0.0869443;
    else if(config==7)
      return 0.155968;
    else if(config==11)
      return 0.08675;
    else if(config==14)
      return 0.0824202;
    else if(config==8)
      return 0.08675;
    else if(config==9)
      return 0.134702;
    else if(config==4363)
      return 0.1;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

  //HCal ADCt elastic peak (ns)
  Double_t atime0(Int_t config) {
    if(config==1)
      return 51.5;
    else if(config==4)
      return 51.466;
    else if(config==7)
      return 63.1685;
    else if(config==11)
      return 50.36;
    else if(config==14)
      return 58.7825;
    else if(config==8)
      return 50.36;
    else if(config==9)
      return 50.6158;
    else if(config==4363)
      return 51.5;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

  //HCal ADCt elastic sigma (ns)
  Double_t atimesig(Int_t config) {
       if(config==1)
      return 3.0;
    else if(config==4)
      return 3.67744;
    else if(config==7)
      return 3.69617;
    else if(config==11)
      return 3.73;
    else if(config==14)
      return 3.68101;
    else if(config==8)
      return 3.73;
    else if(config==9)
      return 3.52261;
    else if(config==4363)
      return 3.0;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

} //::etune
