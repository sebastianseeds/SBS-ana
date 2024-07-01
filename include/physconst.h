#ifndef PHYSCONST_H
#define PHYSCONST_H

#include "TMath.h"

//Via pdg and pdg.lbl.gov/2020/reviews/rpp2020-rev-monte-carlo-numbering.pdf
namespace physconst {

  // math
  static const Double_t pi = TMath::Pi();

  // light
  static const Double_t c = 299792458; // m/s
  static const Int_t IDXc = 22; //mc index
  
  // electron
  static const Double_t qe = 1.602176634E-19; // C
  static const Double_t Me = 0.5109989461E-03; // +/- 31E-13 GeV
  static const Int_t IDXe = 11; //mc index

  // muon
  static const Double_t Mmu = 0.1056583745; // +/- 24E-10 GeV
  static const Int_t IDXmu = 13; //mc index
  
  // proton
  static const Double_t Mp = 0.938272081; // +/- 6E-9 GeV
  static const Double_t mup = 2.7928473446; // +/- 7E-11 muN
  static const Int_t IDXp = 2212; //mc index

  // neutron
  static const Double_t Mn = 0.939565413; // +/- 6E-9 GeV
  static const Double_t mun = -1.9130427; // +/- 5E-7 muN
  static const Int_t IDXn = 2112; //mc index

  // atomic 
  static const Double_t N_A = 6.02214076E23; // 1/mol
  static const Double_t H2_Amass = 1.00784; // u
  static const Double_t D2_Amass = 2.013553212745; // u

  // other
  static const Double_t alpha = 0.00729927; //unitless, fine structure constant

}

#endif
