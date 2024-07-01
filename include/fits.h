#ifndef FITS_H
#define FITS_H

#include <cmath>
//#include "TF1.h"

//Double_t *g_pNfit[7]() = {g_p2fit,g_p3fit,g_p4fit,g_p5fit,g_p6fit,g_p7fit,g_p8fit};


namespace fits{

  Double_t g_gfit( Double_t *x, Double_t *par ); //gaussian fit
  
  Double_t g_sgfit( Double_t *x, Double_t *par ); //skewed gaussian fit

  Double_t g_doubleGausPlusPol( Double_t *x, Double_t *par); //two gaussians plus 4th order bg fit

  Double_t g_expofit( Double_t *x, Double_t *par ); //exponential fit

  Double_t g_scexpofit( Double_t *x, Double_t *par ); //exponential fit with offset

  Double_t g_p0fit( Double_t *x, Double_t *par ); //offset fit

  Double_t g_p1fit( Double_t *x, Double_t *par ); //linear fit

  Double_t g_p2fit( Double_t *x, Double_t *par ); //2nd order polynomial fit
  Double_t g_p2fit_cd( Double_t *x, Double_t *par ); //2nd order polynomial fit concave down only

  Double_t g_p3fit( Double_t *x, Double_t *par ); //3rd order polynomial fit

  Double_t g_p4fit( Double_t *x, Double_t *par ); //4th order polynomial fit

  Double_t g_p5fit( Double_t *x, Double_t *par ); //5th order polynomial fit

  Double_t g_p6fit( Double_t *x, Double_t *par ); //6th order polynomial fit

  Double_t g_p7fit( Double_t *x, Double_t *par ); //7th order polynomial fit

  Double_t g_p8fit( Double_t *x, Double_t *par ); //8th order polynomial fit

  Double_t g_p9fit( Double_t *x, Double_t *par ); //9th order polynomial fit

  Double_t g_p10fit( Double_t *x, Double_t *par ); //10th order polynomial fit

  Double_t g_p11fit( Double_t *x, Double_t *par ); //11th order polynomial fit

  Double_t g_p12fit( Double_t *x, Double_t *par ); //12th order polynomial fit

  /* Double_t ftotal[7]( Double_t *x, Double_t *par ); //Total fit function that expects g_pNfit[] and called by g_interp() */

  /* void g_interp( TH1D *puresamp, TH1D *fulldist, TF1 *tfit, Int_t orderp ); //N-order poly fit to BG + interpolated signal */

  /* void g_interp( TH1D *puresamp, TH1D *fulldist, TF1 *tfit, Int_t orderp, std::string fitname ); //N-order poly fit to BG + interpolated signal (overloaded) */

  void g_conv_gausland( TH1D *h, TF1 *f, Double_t land_mean, Double_t land_sig, Double_t gaus_amp, Double_t gaus_mean, Double_t gaus_sig, Int_t no_FFT_pts ); //gaussian landau convolution

  void g_conv_gausland( TH1D *h, TF1 *f ); //gaussian landau convolution (overloaded)

}

#endif
