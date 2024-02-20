#include "../include/etune.h"

namespace etune {

  //Wide cuts using all subsystems
  std::string globcut(Int_t config) {
    if(config==1)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.075&&bb.gem.track.nhits>3&&abs(bb.etot_over_p-0.92)<0.2&&sbs.hcal.e>0.01&&bb.sh.nclus>0&&sbs.hcal.nclus>0";
    else if(config==4)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.075&&bb.gem.track.nhits>3&&abs(bb.etot_over_p-0.92)<0.2&&sbs.hcal.e>0.02&&bb.sh.nclus>0&&sbs.hcal.nclus>0";
    else if(config==7)
      return "bb.tr.n>0&&abs(bb.tr.vz[0])<0.075&&bb.gem.track.nhits>2&&bb.tr.p[0]>2.0&&sbs.hcal.e>0.06&&bb.ps.e>0.2&&bb.sh.nclus>0&&sbs.hcal.nclus>0";
    else if(config==11)
      return "bb.tr.n>0&&abs(bb.tr.vz[0])<0.075&&bb.gem.track.nhits>2&&bb.tr.p[0]>2.0&&sbs.hcal.e>0.11&&bb.ps.e>0.2&&bb.sh.nclus>0&&sbs.hcal.nclus>0";
    else if(config==14)
      return "bb.tr.n>0&&abs(bb.tr.vz[0])<0.075&&bb.gem.track.nhits>2&&bb.tr.p[0]>1.6&&sbs.hcal.e>0.05&&bb.ps.e>0.2&&bb.sh.nclus>0&&sbs.hcal.nclus>0";
    else if(config==8)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.075&&bb.gem.track.nhits>3&&sbs.hcal.e>0.02&&bb.sh.nclus>0&&sbs.hcal.nclus>0"; 
    else if(config==9)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.075&&bb.gem.track.nhits>3&&sbs.hcal.e>0.02&&bb.sh.nclus>0&&sbs.hcal.nclus>0";
    else if(config==4363)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.gem.track.nhits>3&&abs(bb.etot_over_p-0.92)<0.2&&sbs.hcal.e>0.01&&bb.sh.nclus>0&&sbs.hcal.nclus>0";
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return "";
    }
  }

  //Wide cuts using all subsystems
  std::string globcut_wide(Int_t config) {
    if(config==1)
      return "bb.tr.n>0&&bb.ps.e>0.1&&abs(bb.tr.vz[0])<0.12&&bb.gem.track.nhits>1&&bb.sh.nclus>0&&sbs.hcal.nclus>0";
    else if(config==4)
      return "bb.tr.n>0&&bb.ps.e>0.1&&abs(bb.tr.vz[0])<0.12&&bb.gem.track.nhits>1&&bb.sh.nclus>0&&sbs.hcal.nclus>0";
    else if(config==7)
      return "bb.tr.n>0&&bb.ps.e>0.1&&abs(bb.tr.vz[0])<0.12&&bb.gem.track.nhits>1&&bb.sh.nclus>0&&sbs.hcal.nclus>0";
    else if(config==11)
      return "bb.tr.n>0&&bb.ps.e>0.1&&abs(bb.tr.vz[0])<0.12&&bb.gem.track.nhits>1&&bb.sh.nclus>0&&sbs.hcal.nclus>0";
    else if(config==14)
      return "bb.tr.n>0&&bb.ps.e>0.1&&abs(bb.tr.vz[0])<0.12&&bb.gem.track.nhits>1&&bb.sh.nclus>0&&sbs.hcal.nclus>0";
    else if(config==8)
      return "bb.tr.n>0&&bb.ps.e>0.1&&abs(bb.tr.vz[0])<0.12&&bb.gem.track.nhits>1&&bb.sh.nclus>0&&sbs.hcal.nclus>0";
    else if(config==9)
      return "bb.tr.n>0&&bb.ps.e>0.1&&abs(bb.tr.vz[0])<0.12&&bb.gem.track.nhits>1&&bb.sh.nclus>0&&sbs.hcal.nclus>0";
    else if(config==4363)
      return "bb.tr.n>0&&bb.ps.e>0.1&&abs(bb.tr.vz[0])<0.12&&bb.gem.track.nhits>1&&bb.sh.nclus>0&&sbs.hcal.nclus>0";
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return "";
    }
  }

  //Wide cuts using Bigbite subsystems
  std::string globcut_earm(Int_t config) {
    if(config==1)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.075&&bb.gem.track.nhits>3&&abs(bb.etot_over_p-0.92)<0.2&&bb.sh.nclus>0";
    else if(config==4)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.075&&bb.gem.track.nhits>3&&abs(bb.etot_over_p-0.92)<0.2&&bb.sh.nclus>0";
    else if(config==7)
      return "bb.tr.n>0&&abs(bb.tr.vz[0])<0.075&&bb.gem.track.nhits>2&&bb.ps.e>0.2&&bb.sh.nclus>0";
    else if(config==11)
      return "bb.tr.n>0&&abs(bb.tr.vz[0])<0.075&&bb.gem.track.nhits>2&&bb.ps.e>0.2&&bb.sh.nclus>0";
    else if(config==14)
      return "bb.tr.n>0&&abs(bb.tr.vz[0])<0.075&&bb.gem.track.nhits>2&&bb.tr.p[0]>1.6&&bb.ps.e>0.2&&bb.sh.nclus>0";
    else if(config==8)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.075&&bb.gem.track.nhits>3&&abs(bb.tr.tg_th[0])<0.15&&abs(bb.tr.tg_ph[0])<0.3&&bb.sh.nclus>0"; 
    else if(config==9)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.075&&bb.gem.track.nhits>3&&abs(bb.tr.tg_th[0])<0.15&&abs(bb.tr.tg_ph[0])<0.3&&bb.sh.nclus>0";
    else if(config==4363)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.075&&bb.gem.track.nhits>3&&abs(bb.etot_over_p-0.92)<0.2&&bb.ps.e+bb.sh.e>1.7&&bb.sh.nclus>0";
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return "";
    }
  }

  //Wide cuts using branches available to hcal expert replays
  std::string globcut_hexp(Int_t config) {
    if(config==1)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.ps.e+bb.sh.e>1.7";
    else if(config==4)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.ps.e+bb.sh.e>1.7";
    else if(config==7)
      return "bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.tr.p[0]>2.0&&bb.ps.e>0.2";
    else if(config==11)
      return "bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.tr.p[0]>2.0&&bb.ps.e>0.2";
    else if(config==14)
      return "bb.tr.n>0&&abs(bb.tr.vz[0])<0.08&&bb.tr.p[0]>1.6&&bb.ps.e>0.2";
    else if(config==8)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&abs(bb.tr.tg_th[0])<0.15&&abs(bb.tr.tg_ph[0])<0.3"; 
    else if(config==9)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&abs(bb.tr.tg_th[0])<0.15&&abs(bb.tr.tg_ph[0])<0.3";
    else if(config==4363)
      return "bb.tr.n==1&&bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.08&&bb.ps.e+bb.sh.e>1.7";
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return "";
    }
  }

  //Elastic invariant mass squared peak GeV^2
  Double_t W2mean(Int_t config,Int_t mag) {
    if(config==1)
      return 0.92;
    else if(config==4){
      if(mag==0) return 0.845610;
      else if(mag==30) return 0.913681;
      else if(mag==50) return 0.959941;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS4 (hint: 0, 30, or 50)." << std::endl;
      return -1;
      }
    }else if(config==7){
      if(mag==85) return 0.88;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS7 (hint: it's 85)." << std::endl;
      return -1;
      }
    }else if(config==11){
      if(mag==0) return 0.925;
      else if(mag==100) return 0.925;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS11 (hint: 0 or 100)." << std::endl;
      return -1;
      }
    }else if(config==14){
      if(mag==0) return 0.870819;
      else if(mag==70) return 0.870819;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: 0 or 70)." << std::endl;
      return -1;
      }
    }else if(config==8){
      if(mag==0) return 0.91;
      else if(mag==50) return 0.91;
      else if(mag==70) return 0.91;
      else if(mag==100) return 0.91;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS8 (hint: 0, 50, 70, or 100)." << std::endl;
      return -1;
      }
    }else if(config==9){
      if(mag==70) return 0.91;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: it's 70)." << std::endl;
      return -1;
      }
    }else if(config==4363){
      return 0.92;
    }else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }
  
  //Elastic invariant mass squared sigma GeV^2
  Double_t W2sig(Int_t config,Int_t mag) {
    if(config==1)
      return 0.325;
    else if(config==4){
      if(mag==0) return 0.0909989;
      else if(mag==30) return 0.0911511;
      else if(mag==50) return 0.0873467;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS4 (hint: 0, 30, or 50)." << std::endl;
      return -1;
      }
    }else if(config==7){
      if(mag==85) return 0.5;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS7 (hint: it's 85)." << std::endl;
      return -1;
      }
    }else if(config==11){
      if(mag==0) return 0.325;
      else if(mag==100) return 0.325;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS11 (hint: 0 or 100)." << std::endl;
      return -1;
      }
    }else if(config==14){
      if(mag==0) return 0.19443;
      else if(mag==70) return 0.19443;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: 0 or 70)." << std::endl;
      return -1;
      }
    }else if(config==8){
      if(mag==0) return 0.21;
      else if(mag==50) return 0.21;
      else if(mag==70) return 0.21;
      else if(mag==100) return 0.21;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS8 (hint: 0, 50, 70, or 100)." << std::endl;
      return -1;
      }
    }else if(config==9){
      if(mag==70) return 0.17;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: it's 70)." << std::endl;
      return -1;
      }
    }else if(config==4363){
      return 0.325;
    }else{
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
      else if(mag==30) return -2.91822e-03;
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
      else if(mag==70) return -0.00276996;
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
      if(mag==0) return 0.02764;
      else if(mag==30) return -6.31051e-01;
      else if(mag==50) return -1.09994;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS4 (hint: 0, 30, or 50)." << std::endl;
      return -1;
      }
    }else if(config==7){
      if(mag==85) return -0.6782;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS7 (hint: it's 85)." << std::endl;
      return -1;
      }
    }else if(config==11){
      if(mag==0) return 0.0405545;
      else if(mag==100) return -0.6737;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS11 (hint: 0 or 100)." << std::endl;
      return -1;
      }
    }else if(config==14){
      if(mag==0) return 0.078859;
      else if(mag==70) return -0.7840;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: 0 or 70)." << std::endl;
      return -1;
      }
    }else if(config==8){
      if(mag==0) return 0.0;
      else if(mag==50) return -0.5881;
      else if(mag==70) return -8.15032e-01;
      else if(mag==100) return -1.1973;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS8 (hint: 0, 50, 70, or 100)." << std::endl;
      return -1;
      }
    }else if(config==9){
      if(mag==70) return -0.841747;
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
  Double_t dy0(Int_t config,Int_t mag) {
      if(config==1)
      return 0.;
    else if(config==4){
      if(mag==0) return -0.0424;
      else if(mag==30) return -4.51123e-02;
      else if(mag==50) return -0.0270143;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS4 (hint: 0, 30, or 50)." << std::endl;
      return -1;
      }
    }else if(config==7){
      if(mag==85) return 0.0107;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS7 (hint: it's 85)." << std::endl;
      return -1;
      }
    }else if(config==11){
      if(mag==0) return 0.0428;
      else if(mag==100) return 0.0428;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS11 (hint: 0 or 100)." << std::endl;
      return -1;
      }
    }else if(config==14){
      if(mag==0) return 0.00220912;
      else if(mag==70) return 0.00177;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: 0 or 70)." << std::endl;
      return -1;
      }
    }else if(config==8){
      if(mag==0) return -0.0190;
      else if(mag==50) return -0.0171;
      else if(mag==70) return -5.11974e-02;
      else if(mag==100) return -0.00988;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS8 (hint: 0, 50, 70, or 100)." << std::endl;
      return -1;
      }
    }else if(config==9){
      if(mag==70) return -0.0334186;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: it's 70)." << std::endl;
      return -1;
      }
    }else if(config==4363)
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
      else if(mag==30) return 1.63299e-01;
      else if(mag==50) return 1.54280e-01;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS4 (hint: 0, 30, or 50)." << std::endl;
      return -1;
      }
    }else if(config==7){
      if(mag==85) return 0.183;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS7 (hint: it's 85)." << std::endl;
      return -1;
      }
    }else if(config==11){
      if(mag==0) return 0.0635993;
      else if(mag==100) return 0.324;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS11 (hint: 0 or 100)." << std::endl;
      return -1;
      }
    }else if(config==14){
      if(mag==0) return 0.171009;
      else if(mag==70) return 0.222;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: 0 or 70)." << std::endl;
      return -1;
      }
    }else if(config==8){
      if(mag==0) return 0.196;
      else if(mag==50) return 0.238;
      else if(mag==70) return 0.236;
      else if(mag==100) return 0.217;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS8 (hint: 0, 50, 70, or 100)." << std::endl;
      return -1;
      }
    }else if(config==9){
      if(mag==70) return 0.218;
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
      if(mag==0) return 0.05652;
      else if(mag==30) return 6.85202e-02;
      else if(mag==50) return 1.54280e-01;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS4 (hint: 0, 30, or 50)." << std::endl;
      return -1;
      }
    }else if(config==7){
      if(mag==85) return 0.151;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS7 (hint: it's 85)." << std::endl;
      return -1;
      }
    }else if(config==11){
      if(mag==0) return 0.0635993;
      else if(mag==100) return 0.277;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS11 (hint: 0 or 100)." << std::endl;
      return -1;
      }
    }else if(config==14){
      if(mag==0) return 0.171009;
      else if(mag==70) return 0.197;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: 0 or 70)." << std::endl;
      return -1;
      }
    }else if(config==8){
      if(mag==0) return 0.196;
      else if(mag==50) return 0.220;
      else if(mag==70) return 1.33254e-01;
      else if(mag==100) return 0.246;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS8 (hint: 0, 50, 70, or 100)." << std::endl;
      return -1;
      }
    }else if(config==9){
      if(mag==70) return 0.0826679;
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
  Double_t dysig(Int_t config,Int_t mag) {
       if(config==1)
      return 0.1;
    else if(config==4){
      if(mag==0) return 0.06364;
      else if(mag==30) return 6.08594e-02;
      else if(mag==50) return 6.08594e-02;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS4 (hint: 0, 30, or 50)." << std::endl;
      return -1;
      }
    }else if(config==7){
      if(mag==85) return 0.247;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS7 (hint: it's 85)." << std::endl;
      return -1;
      }
    }else if(config==11){
      if(mag==0) return 0.189;
      else if(mag==100) return 0.189;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS11 (hint: 0 or 100)." << std::endl;
      return -1;
      }
    }else if(config==14){
      if(mag==0) return 0.255;
      else if(mag==70) return 0.255;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: 0 or 70)." << std::endl;
      return -1;
      }
    }else if(config==8){
      if(mag==0) return 0.402;
      else if(mag==50) return 7.15760e-02;
      else if(mag==70) return 7.15760e-02;
      else if(mag==100) return 7.15760e-02;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS8 (hint: 0, 50, 70, or 100)." << std::endl;
      return -1;
      }
    }else if(config==9){
      if(mag==70) return 0.0584762;
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

  //HCal ADCt elastic peak (ns)
  Double_t atime0(Int_t config,Int_t mag) {
    if(config==1)
      return 51.5;
    else if(config==4){
      if(mag==0) return 51.466;
      else if(mag==30) return 51.466;
      else if(mag==50) return 51.466;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS4 (hint: 0, 30, or 50)." << std::endl;
      return -1;
      }
    }else if(config==7){
      if(mag==85) return 63.1685;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS7 (hint: it's 85)." << std::endl;
      return -1;
      }
    }else if(config==11){
      if(mag==0) return 50.36;
      else if(mag==100) return 50.36;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS11 (hint: 0 or 100)." << std::endl;
      return -1;
      }
    }else if(config==14){
      if(mag==0) return 58.7825;
      else if(mag==70) return 58.7825;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: 0 or 70)." << std::endl;
      return -1;
      }
    }else if(config==8){
      if(mag==0) return 50.36;
      else if(mag==50) return 50.36;
      else if(mag==70) return 50.36;
      else if(mag==100) return 50.36;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS8 (hint: 0, 50, 70, or 100)." << std::endl;
      return -1;
      }
    }else if(config==9){
      if(mag==70) return 50.6158;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: it's 70)." << std::endl;
      return -1;
      }
    }else if(config==4363)
      return 51.5;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

  //HCal ADCt elastic sigma (ns)
  Double_t atimesig(Int_t config,Int_t mag) {
       if(config==1)
      return 3.0;
    else if(config==4){
      if(mag==0) return 3.67744;
      else if(mag==30) return 3.67744;
      else if(mag==50) return 3.67744;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS4 (hint: 0, 30, or 50)." << std::endl;
      return -1;
      }
    }else if(config==7){
      if(mag==85) return 3.69617;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS7 (hint: it's 85)." << std::endl;
      return -1;
      }
    }else if(config==11){
      if(mag==0) return 3.73;
      else if(mag==100) return 3.73;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS11 (hint: 0 or 100)." << std::endl;
      return -1;
      }
    }else if(config==14){
      if(mag==0) return 3.68101;
      else if(mag==70) return 3.68101;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: 0 or 70)." << std::endl;
      return -1;
      }
    }else if(config==8){
      if(mag==0) return 3.73;
      else if(mag==50) return 3.73;
      else if(mag==70) return 3.73;
      else if(mag==100) return 3.73;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS8 (hint: 0, 50, 70, or 100)." << std::endl;
      return -1;
      }
    }else if(config==9){
      if(mag==70) return 3.52261;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: it's 70)." << std::endl;
      return -1;
      }
    }else if(config==4363)
      return 3.0;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

  //HCal ADCt - BBCal Shower ADCt elastic peak (ns)
  Double_t atimediff0(Int_t config,Int_t mag) {
    if(config==1)
      return 51.5;
    else if(config==4){
      if(mag==0) return 51.466;
      else if(mag==30) return 51.466;
      else if(mag==50) return 51.466;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS4 (hint: 0, 30, or 50)." << std::endl;
      return -1;
      }
    }else if(config==7){
      if(mag==85) return 63.1685;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS7 (hint: it's 85)." << std::endl;
      return -1;
      }
    }else if(config==11){
      if(mag==0) return 50.36;
      else if(mag==100) return 50.36;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS11 (hint: 0 or 100)." << std::endl;
      return -1;
      }
    }else if(config==14){
      if(mag==0) return 58.7825;
      else if(mag==70) return 58.7825;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: 0 or 70)." << std::endl;
      return -1;
      }
    }else if(config==8){
      if(mag==0) return 50.36;
      else if(mag==50) return 50.36;
      else if(mag==70) return 50.36;
      else if(mag==100) return 50.36;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS8 (hint: 0, 50, 70, or 100)." << std::endl;
      return -1;
      }
    }else if(config==9){
      if(mag==70) return 50.6158;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: it's 70)." << std::endl;
      return -1;
      }
    }else if(config==4363)
      return 51.5;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

  //HCal ADCt - BBCal Shower ADCt elastic sigma (ns)
  Double_t atimediffsig(Int_t config,Int_t mag) {
       if(config==1)
      return 3.0;
    else if(config==4){
      if(mag==0) return 3.67744;
      else if(mag==30) return 3.67744;
      else if(mag==50) return 3.67744;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS4 (hint: 0, 30, or 50)." << std::endl;
      return -1;
      }
    }else if(config==7){
      if(mag==85) return 3.69617;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS7 (hint: it's 85)." << std::endl;
      return -1;
      }
    }else if(config==11){
      if(mag==0) return 3.73;
      else if(mag==100) return 3.73;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS11 (hint: 0 or 100)." << std::endl;
      return -1;
      }
    }else if(config==14){
      if(mag==0) return 3.68101;
      else if(mag==70) return 3.68101;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: 0 or 70)." << std::endl;
      return -1;
      }
    }else if(config==8){
      if(mag==0) return 3.73;
      else if(mag==50) return 3.73;
      else if(mag==70) return 3.73;
      else if(mag==100) return 3.73;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS8 (hint: 0, 50, 70, or 100)." << std::endl;
      return -1;
      }
    }else if(config==9){
      if(mag==70) return 3.52261;
      else{
      std::cerr << "Error: enter a valid magnetic field \% for SBS14 (hint: it's 70)." << std::endl;
      return -1;
      }
    }else if(config==4363)
      return 3.0;
    else{
      std::cerr << "Error: enter a valid SBS kinematic." << std::endl;
      return -1;
    }
  }

} //::etune
