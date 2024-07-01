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

  //skewed gaussian fit
  Double_t g_sgfit(Double_t *x, Double_t *par){
    Double_t amp = par[0];
    Double_t offset = par[1];
    Double_t sigma = par[2];
    Double_t alpha = par[3];

    return amp*exp( -pow( x[0]-offset,2. )/( 2.*pow(sigma,2.) ) )*( 1+erf( (x[0]-offset)*alpha/sigma*sqrt(2.) ) );
  }

  // Combined Gaussian and fourth-order polynomial fit function
  Double_t g_doubleGaussPlusPol( Double_t *x, Double_t *par ) {
    // First Gaussian
    Double_t gauss1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));

    // Second Gaussian
    Double_t gauss2 = par[3] * exp(-0.5 * pow((x[0] - par[4]) / par[5], 2));

    // Fourth-order polynomial
    Double_t poly = par[6] + par[7] * x[0] + par[8] * pow(x[0], 2) + par[9] * pow(x[0], 3) + par[10] * pow(x[0], 4);

    return gauss1 + gauss2 + poly;
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

  //offset fit
  Double_t g_p0fit(Double_t *x, Double_t *par) {
    Double_t yint = par[0];
    return yint;
  }

  //linear fit
  Double_t g_p1fit(Double_t *x, Double_t *par){
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

  //2nd order poly fit concave down only for small background slices
  Double_t g_p2fit_cd(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = -fabs(par[2]);
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

  //9th order poly fit
  Double_t g_p9fit(Double_t *x, Double_t *par){
    Double_t fit = 0.;
    for( Int_t p=0; p<9; p++ ){
      fit = par[p]*pow(x[0],p);
    }
    return fit;
  }


  //10th order poly fit
  Double_t g_p10fit(Double_t *x, Double_t *par){
    Double_t fit = 0.;
    for( Int_t p=0; p<10; p++ ){
      fit = par[p]*pow(x[0],p);
    }
    return fit;
  }


  //11th order poly fit
  Double_t g_p11fit(Double_t *x, Double_t *par){
    Double_t fit = 0.;
    for( Int_t p=0; p<11; p++ ){
      fit = par[p]*pow(x[0],p);
    }
    return fit;
  }


  //12th order poly fit
  Double_t g_p12fit(Double_t *x, Double_t *par){
    Double_t yint = par[0];
    Double_t p1 = par[1];
    Double_t p2 = par[2];
    Double_t p3 = par[3];
    Double_t p4 = par[4];
    Double_t p5 = par[5];
    Double_t p6 = par[6];
    Double_t p7 = par[7];
    Double_t p8 = par[8];
    Double_t p9 = par[9];
    Double_t p10 = par[10];
    Double_t p11 = par[11];
    Double_t p12 = par[12];
    return yint+p1*x[0]+p2*pow(x[0],2)+p3*pow(x[0],3)+p4*pow(x[0],4)+p5*pow(x[0],5)+p6*pow(x[0],6)+p7*pow(x[0],7)+p8*pow(x[0],8)+p9*pow(x[0],9)+p10*pow(x[0],10)+p11*pow(x[0],11)+p12*pow(x[0],12);
  }
  
  

  // //N-order poly fit to signal+background on histogram with signal interpolated from "pure" signal histogram (same limits)
  // //Double_t *g_pNfit[7]() = {g_p2fit,g_p3fit,g_p4fit,g_p5fit,g_p6fit,g_p7fit,g_p8fit};
  // Double_t ftotal[7](Double_t *x, Double_t *par) =
  // {
  //   {//get poly fit to bg with scaled fit to "pure elastics"
  //     Double_t distX = x[0];
  //     Double_t sig_scale = par[0];
  //     Double_t signal = sig_scale * puresamp->Interpolate(distX);
  //     return signal + g_p2fit(x,&par[1]);
  //   },
  //   {//get poly fit to bg with scaled fit to "pure elastics"
  //     Double_t distX = x[0];
  //     Double_t sig_scale = par[0];
  //     Double_t signal = sig_scale * puresamp->Interpolate(distX);
  //     return signal + g_p3fit(x,&par[1]);
  //   },
  //   {//get poly fit to bg with scaled fit to "pure elastics"
  //     Double_t distX = x[0];
  //     Double_t sig_scale = par[0];
  //     Double_t signal = sig_scale * puresamp->Interpolate(distX);
  //     return signal + g_p4fit(x,&par[1]);
  //   },
  //   {//get poly fit to bg with scaled fit to "pure elastics"
  //     Double_t distX = x[0];
  //     Double_t sig_scale = par[0];
  //     Double_t signal = sig_scale * puresamp->Interpolate(distX);
  //     return signal + g_p5fit(x,&par[1]);
  //   },
  //   {//get poly fit to bg with scaled fit to "pure elastics"
  //     Double_t distX = x[0];
  //     Double_t sig_scale = par[0];
  //     Double_t signal = sig_scale * puresamp->Interpolate(distX);
  //     return signal + g_p6fit(x,&par[1]);
  //   },
  //   {//get poly fit to bg with scaled fit to "pure elastics"
  //     Double_t distX = x[0];
  //     Double_t sig_scale = par[0];
  //     Double_t signal = sig_scale * puresamp->Interpolate(distX);
  //     return signal + g_p7fit(x,&par[1]);
  //   },
  //   {//get poly fit to bg with scaled fit to "pure elastics"
  //     Double_t distX = x[0];
  //     Double_t sig_scale = par[0];
  //     Double_t signal = sig_scale * puresamp->Interpolate(distX);
  //     return signal + g_p8fit(x,&par[1]);
  //   }
  // }

  // void g_interp( TH1D *puresamp, TH1D *fulldist, TF1 *tfit, Int_t orderp ){
  //   //Double_t *g_pNfit[7]() = {g_p2fit,g_p3fit,g_p4fit,g_p5fit,g_p6fit,g_p7fit,g_p8fit};
  //   Double_t xmin = fulldist->GetXaxis()->GetXmin();
  //   Double_t xmax = fulldist->GetXaxis()->GetXmax();

  //   // Double_t ftotal(Double_t *x, Double_t *par){ //get poly fit to bg with scaled fit to "pure elastics"
  //   //   Double_t distX = x[0];
  //   //   Double_t sig_scale = par[0];
  //   //   Double_t signal = sig_scale * puresamp->Interpolate(distX);
  //   //   return signal + g_pNfit[orderp-2](x,&par[1]);
  //   // }
    
  //   tfit = new TF1("tfit",ftotal[orderp-2],xmin,xmax,orderp+2);
  //   tfit->SetLineColor(kGreen);

  //   fulldist->Fit("tfit","RBM");
  // }

  // //N-order poly fit to signal+background on histogram with signal interpolated from "pure" signal histogram (same limits, overloaded to allow user control over fit function name)
  // void g_interp( TH1D *puresamp, TH1D *fulldist, TF1 *tfit, Int_t orderp, std::string fitname ){
  //   //Double_t *g_pNfit[7]() = {g_p2fit,g_p3fit,g_p4fit,g_p5fit,g_p6fit,g_p7fit,g_p8fit};
  //   Double_t xmin = fulldist->GetXaxis()->GetXmin();
  //   Double_t xmax = fulldist->GetXaxis()->GetXmax();

  //   // Double_t ftotal(Double_t *x, Double_t *par){ //get poly fit to bg with scaled fit to "pure elastics"
  //   //   Double_t distX = x[0];
  //   //   Double_t sig_scale = par[0];
  //   //   Double_t signal = sig_scale * puresamp->Interpolate(distX);
  //   //   return signal + g_pNfit[orderp-2](x,&par[1]);
  //   // }
    
  //   tfit = new TF1(fitname,ftotal,xmin,xmax,orderp+2);
  //   tfit->SetLineColor(kGreen);

  //   fulldist->Fit(fitname,"RBM");
  // }

  //Gaussian Landau convolution
  void g_conv_gausland( TH1D *h, TF1 *f, Double_t land_mean, Double_t land_sig, Double_t gaus_amp, Double_t gaus_mean, Double_t gaus_sig, Int_t no_FFT_pts ){
    Double_t xmin = h->GetXaxis()->GetXmin();
    Double_t xmax = h->GetXaxis()->GetXmax();
    TF1Convolution *f_conv = new TF1Convolution("landau","gaus",xmin,xmax);
    f_conv->SetRange(xmin,xmax);
    f_conv->SetNofPointsFFT(no_FFT_pts);
    f = new TF1("glc_fit",*f_conv,xmin,xmax,f_conv->GetNpar());
    f->SetParameters( land_mean, land_sig, gaus_amp, gaus_mean, gaus_sig );

    h->Fit( "f", "RQ", "", xmin, xmax);
  }

  //Guassian Landau convolution (overloaded)
  void g_conv_gausland( TH1D *h, TF1 *f ){
    Double_t land_mean=1., land_sig=1., gaus_amp=1., gaus_mean=1., gaus_sig=1.;
    Int_t no_FFT_pts=1000;
    Double_t xmin = h->GetXaxis()->GetXmin();
    Double_t xmax = h->GetXaxis()->GetXmax();
    TF1Convolution *f_conv = new TF1Convolution("landau","gaus",xmin,xmax);
    f_conv->SetRange(xmin,xmax);
    f_conv->SetNofPointsFFT(no_FFT_pts);
    f = new TF1("glc_fit",*f_conv,xmin,xmax,f_conv->GetNpar());
    f->SetParameters( land_mean, land_sig, gaus_amp, gaus_mean, gaus_sig );

    h->Fit( "f", "RQ", "", xmin, xmax);
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
