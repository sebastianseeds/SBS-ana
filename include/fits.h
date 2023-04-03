#ifndef FITS_H
#define FITS_H

#include <cmath>
//#include "TF1.h"

namespace fits{

  Double_t g_gfit( Double_t *x, Double_t *par );  //gaussian fit

  Double_t g_expofit(Double_t *x, Double_t *par); //exponential fit

  Double_t g_scexpofit(Double_t *x, Double_t *par); //exponetial fit with offset

  Double_t g_lfit(Double_t *x, Double_t *par); //linear fit

  Double_t g_p2fit(Double_t *x, Double_t *par); //2nd order polynomial fit

  Double_t g_p3fit(Double_t *x, Double_t *par); //3rd order polynomial fit

  Double_t g_p4fit(Double_t *x, Double_t *par); //4th order polynomial fit

  Double_t g_p5fit(Double_t *x, Double_t *par); //5th order polynomial fit

  Double_t g_p6fit(Double_t *x, Double_t *par); //6th order polynomial fit

  Double_t g_p7fit(Double_t *x, Double_t *par); //7th order polynomial fit

  Double_t g_p8fit(Double_t *x, Double_t *par); //8th order polynomial fit

}

#endif
