#include "../include/fits.h"

//Total fits should be made of these and configured by script.
namespace fits {

  //gaussian fit
  Double_t g_gfit(Double_t *x, Double_t *par){
    Double_t amp = par[0];
    Double_t offset = par[1];
    Double_t sigma = par[2];
    return amp*exp(-0.5*pow((x[0]-offset)/sigma,2.));
  }

  //expo fit
  Double_t g_expofit(Double_t *x, Double_t *par){
    Double_t amp = par[0];
    Double_t offset = par[1];
    Double_t str = par[2];
    return amp*exp(offset+str*x[0]);
  }

  //expo fit with offset
  Double_t g_scexpofit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t amp = par[1];
    Double_t offset = par[2];
    Double_t str = par[3];
    return yint+amp*exp(offset+str*x[0]);
  }

  //linear fit
  Double_t g_lfit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];;
    return yint+p1*x[0];
  }

  //2nd order poly fit
  Double_t g_p2fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    return yint+p1*x[0]+p2*pow(x[0],2);
  }

  //3rd order poly fit
  Double_t g_p3fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3);
  }

  //4th order poly fit
  Double_t g_p4fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    Double_t p4 = par[4];
    return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3)+p4*pow(x[0],4);
  }

  //5th order poly fit
  Double_t g_p5fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    Double_t p4 = par[4];
    Double_t p5 = par[5];
    return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3)+p4*pow(x[0],4)+p5*pow(x[0],5);
  }

  //6th order poly fit
  Double_t g_p6fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    Double_t p4 = par[4];
    Double_t p5 = par[5];
    Double_t p6 = par[6];
    return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3)+p4*pow(x[0],4)+p5*pow(x[0],5)+p6*pow(x[0],6);
  }

  //7th order poly fit
  Double_t g_p7fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    Double_t p4 = par[4];
    Double_t p5 = par[5];
    Double_t p6 = par[6];
    Double_t p7 = par[7];
    return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3)+p4*pow(x[0],4)+p5*pow(x[0],5)+p6*pow(x[0],6)+p7*pow(x[0],7);
  }

  //8th order poly fit
  Double_t g_p8fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    Double_t p4 = par[4];
    Double_t p5 = par[5];
    Double_t p6 = par[6];
    Double_t p7 = par[7];
    Double_t p8 = par[8];
    return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3)+p4*pow(x[0],4)+p5*pow(x[0],5)+p6*pow(x[0],6)+p7*pow(x[0],7)+p8*pow(x[0],8);
  }

  //Reject Point fits. Limits configured by script.
  /*
  //Side-band gaussian fit (fits only wings/ends of histogram)
  Double_t SBgrej_b; //Side-band reject begin
  Double_t SBgrej_e; //Side-band reject end

  Double_t g_SBgfit(double *x, double *par){
  Double_t amp = par[0];
  Double_t loc = par[1];
  Double_t sigma = par[2];
  if(x[0]>SBgrej_b && x[0]<SBgrej_e) { 
  TF1::RejectPoint();
  return 0;
  }
  return amp*exp(-0.5*pow((x[0]-loc)/sigma),2.);
  }

  //Side-band expo fit (fits only wings/ends of histogram)
  Double_t SBexporej_b; //Side-band reject begin
  Double_t SBexporej_e; //Side-band reject end

  Double_t g_SBexpofit(double *x, double *par){
  Double_t amp = par[0];
  Double_t offset = par[1];
  Double_t str = par[2];
  if(x[0]>SBexporej_b && x[0]<SBexporej_e) { 
  TF1::RejectPoint();
  return 0;
  }
  return amp*exp(offset+str*x[0]);
  }

  //Central-band gaussian fit (fits only center of histogram)
  Double_t CBgrej_b; //Central-band fit begin
  Double_t CBgrej_e; //Central-band fit end

  Double_t g_CBgfit(double *x, double *par){
  Double_t amp = par[0];
  Double_t loc = par[1];
  Double_t sigma = par[2];
  if(x[0]<CBgrej_b || x[0]>CBgrej_e) { 
  TF1::RejectPoint();
  return 0;
  }
  return amp*exp(-0.5*pow((x[0]-loc)/sigma),2.);
  }

  //Central-band expo fit (fits only center of histogram)
  Double_t CBexporej_b; //Central-band fit begin
  Double_t CBexporej_e; //Central-band fit end

  Double_t g_CBexpofit(double *x, double *par){
  Double_t amp = par[0];
  Double_t offset = par[1];
  Double_t str = par[2];
  if(x[0]>CBexporej_b && x[0]<CBexporej_e) { 
  TF1::RejectPoint();
  return 0;
  }
  return amp*exp(offset+str*x[0]);
  }
  */

} //::fits
