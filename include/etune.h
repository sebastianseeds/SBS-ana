#ifndef ETUNE_H
#define ETUNE_H

#include <iostream>
#include <string>
#include <fstream>
#include "TMath.h"
#include "TString.h"

namespace etune {

  // Any tuneable constants for all configurations go here

  // Following quantities vary with configuration
  std::string  globcut(Int_t config);      //List of wide cuts over all subsystems
  std::string  globcut_earm(Int_t config); //List of wide cuts over BigBite subsystems
  Double_t     W2mean(Int_t config);
  Double_t     W2sig(Int_t config);
  Double_t     dx0_n(Int_t config,Int_t mag);
  Double_t     dx0_p(Int_t config,Int_t mag);
  Double_t     dy0(Int_t config);
  Double_t     dxsig_n(Int_t config,Int_t mag);
  Double_t     dxsig_p(Int_t config,Int_t mag);
  Double_t     dysig(Int_t config);
  Double_t     atime0(Int_t config);
  Double_t     atimesig(Int_t config);
}

// a class for SBS config
class SBStune {
 public:

  Int_t        GetSBStconf()       const { return fSBStconf; }
  Int_t        GetSBStmag()        const { return fSBStmag; }
  std::string  Getglobcut()        const { return fglobcut; }
  std::string  Getglobcut_earm()   const { return fglobcut_earm; }
  Double_t     GetW2mean()         const { return fW2mean; }
  Double_t     GetW2sig()          const { return fW2sig; }
  Double_t     Getdx0_n()          const { return fdx0_n; }
  Double_t     Getdx0_p()          const { return fdx0_p; }
  Double_t     Getdy0()            const { return fdy0; }
  Double_t     Getdxsig_n()        const { return fdxsig_n; }
  Double_t     Getdxsig_p()        const { return fdxsig_p; }
  Double_t     Getdysig()          const { return fdysig; }
  Double_t     Getatime0()         const { return fatime0; }
  Double_t     Getatimesig()       const { return fatimesig; }

  // constructor
  SBStune(Int_t conf, Int_t sbsmag) {
    fSBStconf       = conf;
    fSBStmag        = sbsmag;
    fglobcut        = etune::globcut(conf);
    fglobcut_earm   = etune::globcut_earm(conf);
    fW2mean         = etune::W2mean(conf);
    fW2sig          = etune::W2sig(conf);
    fdx0_n          = etune::dx0_n(conf,sbsmag);
    fdx0_p          = etune::dx0_p(conf,sbsmag);
    fdy0            = etune::dy0(conf);
    fdxsig_n        = etune::dxsig_n(conf,sbsmag);
    fdxsig_p        = etune::dxsig_p(conf,sbsmag);
    fdysig          = etune::dysig(conf);
    fatime0         = etune::atime0(conf);
    fatimesig       = etune::atimesig(conf);
  }

  // define an ostream operator to print to screen conveniently
  friend ostream& operator <<(ostream &out, const SBStune& sbstune) {
    out  << " -------------------------- "                                                        << std::endl
	 << Form(" SBS Tune Config: %d, "                         , sbstune.fSBStconf)             << std::endl
	 << Form(" SBS Tune Magnet Settings: %d (p), "            , sbstune.fSBStmag)              << std::endl
    	 << Form(" Global Cut: %s,"                               , sbstune.fglobcut.c_str())     << std::endl
    	 << Form(" Global Cut Earm: %s,"                          , sbstune.fglobcut_earm.c_str())<< std::endl
    	 << Form(" W2 Elastic Peak Mean : %0.5f (GeV),"           , sbstune.fW2mean)              << std::endl
      	 << Form(" W2 Elastic Peak Sigma: %0.5f (GeV),"           , sbstune.fW2sig)               << std::endl
    	 << Form(" HCal dx Elastic Neutron Peak Mean: %0.5f (m)," , sbstune.fdx0_n)               << std::endl
    	 << Form(" HCal dx Elastic Proton Peak Mean: %0.5f (m),"  , sbstune.fdx0_p)               << std::endl
    	 << Form(" HCal dy Elastic Peak Mean: %0.5f (m),"         , sbstune.fdy0)                 << std::endl
    	 << Form(" HCal dx Elastic Neutron Peak Sigma: %0.5f (m),", sbstune.fdxsig_n)             << std::endl
    	 << Form(" HCal dx Elastic Proton Peak Sigma: %0.5f (m)," , sbstune.fdxsig_p)             << std::endl
    	 << Form(" HCal dy Elastic Peak Sigma: %0.5f (m),"        , sbstune.fdysig)               << std::endl
    	 << Form(" HCal Elastic ADCt Mean: %0.5f (m),"            , sbstune.fatime0)              << std::endl
    	 << Form(" HCal Elastic ADCt Sigma: %0.5f (m),"           , sbstune.fatimesig)            << std::endl
	 << " -------------------------- "                        << std::endl                    << std::endl;
    return out;
  }

 private:
  Int_t        fSBStconf;         // SBS configuration number
  Int_t        fSBStmag;          // SBS magnet settings (%)
  std::string  fglobcut;          // wide elastic global cut using all subsystems
  std::string  fglobcut_earm;     // wide elastic global cut using BigBite subsystems only
  Double_t     fW2mean;           // Location of elastic peak in invariant mass distribution (GeV)
  Double_t     fW2sig;            // Width of elastic peak in invariant mass distribution 1-sigma (GeV)
  Double_t     fdx0_n;            // Location of neutron elastic peak in HCal dx distribution (m)
  Double_t     fdx0_p;            // Location of proton elastic peak in HCal dx distribution (m)
  Double_t     fdy0;              // Location of elastic peak in HCal dy distribution (m)
  Double_t     fdxsig_n;          // Width of neutron elastic peak in HCal dx distribution (m)
  Double_t     fdxsig_p;          // Width of proton elastic peak in HCal dx distribution (m)
  Double_t     fdysig;            // Width of elastic peak in HCal dy distribution (m)
  Double_t     fatime0;           // ADC time wrt BB trigger for elastic events (mean of distribution) (ns)
  Double_t     fatimesig;         // ADC time wrt BB trigger for elastic events (sigma of distribution) (ns)
};

#endif
