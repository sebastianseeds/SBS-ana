//SSeeds 9.1.23 Script to extract the p/n yield ratio via data/MC comparison. Intent is to loop over all data from a given set/field setting, apply elastic cuts, and build a dx histogram. From this point, read in MC (with RC) distributions for quasi-elastic protons and neutrons (independently), fit these distributions with several distributions to get a best fit, then compose a sum of these fits allowing for a scaling parameter which maps to the total MC yields. Next, do the same with g4sbs inelastic generator and add this to the total fit function (sum) and apply the now three floating parameter fit to the data and check residuals and chi-square. Varying the fit function may yield better results, but the ratio will be extracted. See gmn_calc.C for application of this ratio to gmn extraction.

#include <vector>
#include <iostream>
#include "TCut.h"
#include "TLorentzVector.h"
#include "TTreeFormula.h"
#include "TChain.h"
#include "TLatex.h"
#include "TH1D.h"
#include "TStopwatch.h"
#include "TMatrixD.h"
#include "../../src/jsonmgr.C"
#include "../../include/gmn.h"

const int atimeSigFac = 5;
const int atimediffSigFac = 5;
const double dx_fitlow = -0.3;
const double dx_fithigh = 0.3;

const double bfac = 0.0842;

//Total fits using Interpolate with elastic signal histo and 6th order poly fit to bg
TH1D *hdx_p;
TH1D *hdx_p_clone;
TH1D *hdx_n;
TH1D *hdx_n_clone;
TH1D *hdx_bg;

Double_t dxtotal(Double_t *x, Double_t *par){
  Double_t dx = x[0];
  Double_t proton_scale = par[0];
  Double_t neutron_scale = par[1];
  Double_t bg_scale = par[2];
  Double_t proton = proton_scale * hdx_p->Interpolate(dx);
  Double_t neutron = neutron_scale * hdx_n->Interpolate(dx);
  Double_t bg = bg_scale * hdx_bg->Interpolate(dx);
  return proton + neutron + bg;
}

Double_t dxtotal_bgfit(Double_t *x, Double_t *par){
  Double_t dx = x[0];
  Double_t proton_scale = par[0];
  Double_t neutron_scale = par[1];
  Double_t proton = proton_scale * hdx_p_clone->Interpolate(dx);
  Double_t neutron = neutron_scale * hdx_n_clone->Interpolate(dx);
  return proton + neutron + fits::g_p4fit(x,&par[2]);
}

double IntegrateTH1D(TH1D* hist) {
    double integral = 0.0;
    
    for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
        integral += hist->GetBinContent(bin) * hist->GetBinWidth(bin);
    }
    
    return integral;
}

//main
void mc_data_compare( Int_t kine=8, Int_t magset=70 )
{  

  // Define a clock to check macro processing time
  TStopwatch *st = new TStopwatch();
  st->Start( kTRUE );

  gStyle->SetOptFit();
  gStyle->SetEndErrorSize(0);

  // reading input config file
  JSONManager *jmgr = new JSONManager("../../config/gmn.json");
  Int_t epm = jmgr->GetValueFromSubKey<Int_t>( "epm", Form("sbs%d",kine) );
  Int_t cluster_idx = jmgr->GetValueFromSubKey<Int_t>( "cluster_idx", Form("sbs%d",kine) );
  Double_t hcal_v_offset = jmgr->GetValueFromSubKey<Double_t>( "hcal_offset", Form("sbs%d",kine) );

  // set up paths
  std::string outdir_path = gSystem->Getenv("OUT_DIR");
  std::string fout_path = outdir_path + Form("/gmn_analysis/mc_data_compare_out_sbs%d_mag%d.root",kine,magset);

  //TFile *fout = new TFile( fout_path.c_str(), "RECREATE" );
  TFile *fout = new TFile( "testout.root", "RECREATE" );

  std::string mc_path = outdir_path + Form("/gmn_analysis/mc_rc_both_out_sbs%d_mag%d.root",kine,magset);
  //std::string data_path = outdir_path + Form("/gmn_analysis/gmn_elastic_out_sbs%d_mag%d_clusteridx%d_epm%d.root",kine,magset,cluster_idx,epm);
  std::string data_path = outdir_path + Form("/gmn_analysis/gmn_elastic_fid_out_sbs%d_mag%d_clusteridx%d_epm%d.root",kine,magset,cluster_idx,epm);
  //std::string bg_path = outdir_path + Form("/gmn_analysis/mc_inelastic_out_sbs%d_mag0.root",kine);

  std::string bg_path = outdir_path + "/gmn_analysis/mc_inelastic_out_sbs8_mag0.root";

  // open files and get histograms
  TFile *mc_file = TFile::Open(mc_path.c_str());
  hdx_p = (TH1D*)mc_file->Get("hdx_cut_p");
  hdx_p_clone = (TH1D*)(hdx_p->Clone("hdx_p_clone"));
  TH1D *hdx_p_clone2 = (TH1D*)(hdx_p->Clone("hdx_p_clone2"));
  hdx_n = (TH1D*)mc_file->Get("hdx_cut_n");
  hdx_n_clone = (TH1D*)(hdx_n->Clone("hdx_n_clone"));
  TH1D *hdx_n_clone2 = (TH1D*)(hdx_n->Clone("hdx_n_clone2"));

  TFile *bg_file = TFile::Open(bg_path.c_str());
  hdx_bg = (TH1D*)bg_file->Get("hdx_cut");

  TFile *data_file = TFile::Open(data_path.c_str());
  TH1D *hdx_data = (TH1D*)data_file->Get("hdx_cut");
  TH1D *hdx_data_clone = (TH1D*)(hdx_data->Clone("hdx_data_clone")); 
  TH2F *hxy = (TH2F*)data_file->Get("hHcalXY");
  //TH2D *hxy_clone = (TH2D*)(hdx_data->Clone("hxy_clone")); 
  TH2F *hxy_fcut = (TH2F*)data_file->Get("hHcalXY_allcuts");
  //TH2D *hxy_fcut_clone = (TH2D*)(hdx_data->Clone("hxy_fcut_clone")); 

  // Set up fit functions for comparisons
  Double_t hcalfit_l = econst::hcalposXi_mc; //lower fit/bin limit for hcal dx plots (m)
  Double_t hcalfit_h = econst::hcalposXf_mc; //upper fit/bin limit for hcal dx plots (m)

  TF1 *tfitdx = new TF1("tfitdx",dxtotal,hcalfit_l,hcalfit_h,3);
  tfitdx->SetNpx(5000);  // Default is usually 100
  TF1 *tfitdx_bgfit = new TF1("tfitdx_bgfit",dxtotal_bgfit,hcalfit_l,hcalfit_h,7); //tfit npar = 1+pNfit_npar+1
  tfitdx_bgfit->SetNpx(5000);  // Default is usually 100

  // Check if histograms are loaded
  if (!hdx_p || !hdx_n || !hdx_bg) {
    std::cerr << "Failed to load histograms." << std::endl;
    return;
  }

  //Make canvas for all interpolate fit
  TCanvas *c1 = new TCanvas("c1","All interpolate",1200,500);
  c1->cd();

  tfitdx->SetLineColor(kGreen);
  hdx_data->SetTitle("dx");
  hdx_data->Draw();


  //hdx_p->Scale(tpar[0]);
  hdx_p->SetLineColor(kRed);
  //hdx_p->SetFillColor(kRed);
  hdx_p->SetFillStyle(3005);
  hdx_p->Draw("same");

  // hdx_n->Scale(tpar[1]);
  hdx_n->SetLineColor(kBlue);
  //hdx_n->SetFillColor(kBlue);
  hdx_n->SetFillStyle(3003);
  hdx_n->Draw("same");

  // hdx_bg->Scale(tpar[2]);
  // hdx_bg->SetLineColor(kBlack);
  // hdx_bg->SetFillColor(kBlack);
  // hdx_bg->SetFillStyle(3004);
  // hdx_bg->Draw("same");
  
  c1->Update();

  c1->Write();

  hdx_data->Fit("tfitdx","RBM0");

  double *tpar = tfitdx->GetParameters();


  //Make canvas for interpolate fit with fourth order poly background
  TCanvas *c2 = new TCanvas("c2","Interpolate/4th Order Poly BG",1200,500);
  c2->cd();

  tfitdx_bgfit->SetLineColor(kGreen);
  hdx_data_clone->SetTitle("dx poly bg fit");
  hdx_data_clone->Fit("tfitdx_bgfit","RBM");
  hdx_data_clone->Draw();

  double *tpar_clone = tfitdx_bgfit->GetParameters();

  hdx_p_clone2->Scale(tpar_clone[0]);
  hdx_p_clone2->SetLineColor(kOrange);
  hdx_p_clone2->SetLineWidth(2);
  hdx_p_clone->SetFillColorAlpha(kOrange,0.3);
  hdx_p_clone->SetFillStyle(3005);
  hdx_p_clone2->Draw("same");

  hdx_n_clone2->Scale(tpar_clone[1]);
  hdx_n_clone2->SetLineColor(kBlue);
  hdx_n_clone2->SetLineWidth(2);
  hdx_n_clone->SetFillColorAlpha(kBlue,0.3);
  hdx_n_clone->SetFillStyle(3003);
  hdx_n_clone2->Draw("same");

  TF1 *bg = new TF1("bg",fits::g_p4fit,hcalfit_l,hcalfit_h,5);
  bg->SetParameters(&tpar_clone[2]);
  bg->SetLineColor(kRed);
  bg->SetFillColor(kRed);
  bg->SetFillStyle(3004);
  bg->Draw("same");
  
  //Get integrals of scaled proton and neutron histograms
  double Np = IntegrateTH1D(hdx_p_clone2);
  double Nn = IntegrateTH1D(hdx_n_clone2);

  //double pn_ratio = Np/Nn * bfac;

  //blind for now
  double pn_ratio = 0.;

  double pscale = tpar[0];
  double nscale = tpar[1];
  double np_par_ratio = nscale/pscale * bfac;

  // double pscale = 0;
  // double nscale = 0;
  // double np_par_ratio = 0;

  //Add a legend to the canvas
  auto leg = new TLegend(0.1,0.6,0.5,0.9);
  leg->AddEntry( hdx_data_clone, "Tight Elastic Cut, Data", "l");
  leg->AddEntry( hdx_p_clone2, "MC, proton", "l");
  leg->AddEntry( hdx_n_clone2, "MC, neutron", "l");
  leg->AddEntry( bg, "poly fit to BG", "l");
  leg->AddEntry( tfitdx_bgfit, "Total fit", "l");
  leg->AddEntry( (TObject*)0, "", "");
  //leg->AddEntry( (TObject*)0, Form("proton scale: %0.3f",pscale), "");
  //leg->AddEntry( (TObject*)0, Form("neutron scale: %0.3f",nscale), "");
  //leg->AddEntry( (TObject*)0, Form("n/p scale ratio: %0.3f",np_par_ratio), "");
  leg->AddEntry( (TObject*)0, Form("proton/neutron ratio: %0.3f",0.), "");
  leg->Draw();


  c2->Update();

  c2->SaveAs("/work/halla/sbs/seeds/HCal_replay/hcal/hcalCalibration/SBS/quality_plots/gmn_data_mc.pdf");


  c2->Write();

  //Make canvas for residual analysis
  TCanvas *c3 = new TCanvas("c3","Data/MC Residuals",1200,500);
  c3->cd();

  //add the p and n histograms together for comparison
  TH1D* sumHistogram = new TH1D("sumHistogram", "Sum of Histograms", hdx_p_clone2->GetNbinsX(), hdx_p_clone2->GetXaxis()->GetXmin(), hdx_p_clone2->GetXaxis()->GetXmax());
  sumHistogram->Add(hdx_p_clone2, hdx_n_clone2);

  TH1D *hRes = util::makeResidualHisto("dx",hdx_data_clone,sumHistogram,true,false);

  //hRes
  hRes->Scale(0.1);

  hRes->SetLineColor(kRed);
  hRes->SetLineWidth(2);
  // hdx_n_clone->SetFillColor(kBlue);
  // hdx_n_clone->SetFillStyle(3003);
  hRes->Draw("same");

  c3->Update();
  c3->Write();

  //Make canvas for data fiducial and active area plots
  //Get cuts from objects
  SBStune tune(kine,magset);
  Double_t dxsig_p  = tune.Getdxsig_p();
  Double_t dysig    = tune.Getdysig();
  Double_t dx0_p    = tune.Getdx0_p();

  vector<Double_t> hcaledge = cut::hcalaa_mc(0,0);
  Double_t hcalx_t = hcaledge[2];
  Double_t hcalx_b = hcaledge[3];
  Double_t hcaly_r = hcaledge[0]-hcal_v_offset;
  Double_t hcaly_l = hcaledge[1]-hcal_v_offset;
  vector<Double_t> hcalaa = cut::hcalaa_mc(1,1);
  Double_t hcalx_aat = hcalaa[2];
  Double_t hcalx_aab = hcalaa[3];
  Double_t hcaly_aar = hcalaa[0]-hcal_v_offset;
  Double_t hcaly_aal = hcalaa[1]-hcal_v_offset;
  vector<Double_t> fid = cut::hcalfid(dxsig_p,dysig,hcalaa);
  Double_t hcalx_ft = fid[2];
  Double_t hcalx_fb = fid[3];
  Double_t hcaly_fr = fid[0]-hcal_v_offset;
  Double_t hcaly_fl = fid[1]-hcal_v_offset;

  //Make Edge TLines
  TLine* l1_e = new TLine(hcalx_t, hcaly_r, hcalx_t, hcaly_l);
  l1_e->SetLineColor(kRed);
  TLine* l2_e = new TLine(hcalx_b, hcaly_r, hcalx_b, hcaly_l);
  l2_e->SetLineColor(kRed);
  TLine* l3_e = new TLine(hcalx_t, hcaly_r, hcalx_b, hcaly_r);
  l3_e->SetLineColor(kRed);
  TLine* l4_e = new TLine(hcalx_t, hcaly_l, hcalx_b, hcaly_l);
  l4_e->SetLineColor(kRed);

  //Make Active Area TLines
  TLine* l1_aa = new TLine(hcalx_aat, hcaly_aar, hcalx_aat, hcaly_aal);
  l1_aa->SetLineColor(kBlue);
  TLine* l2_aa = new TLine(hcalx_aab, hcaly_aar, hcalx_aab, hcaly_aal);
  l2_aa->SetLineColor(kBlue);
  TLine* l3_aa = new TLine(hcalx_aat, hcaly_aar, hcalx_aab, hcaly_aar);
  l3_aa->SetLineColor(kBlue);
  TLine* l4_aa = new TLine(hcalx_aat, hcaly_aal, hcalx_aab, hcaly_aal);
  l4_aa->SetLineColor(kBlue);

  //Make Fiducial Safety Margin TLines
  TLine* l1_f = new TLine(hcalx_ft, hcaly_fr, hcalx_ft, hcaly_fl);
  l1_f->SetLineColor(kGreen);
  TLine* l2_f = new TLine(hcalx_fb, hcaly_fr, hcalx_fb, hcaly_fl);
  l2_f->SetLineColor(kGreen);
  TLine* l3_f = new TLine(hcalx_ft, hcaly_fr, hcalx_fb, hcaly_fr);
  l3_f->SetLineColor(kGreen);
  TLine* l4_f = new TLine(hcalx_ft, hcaly_fl, hcalx_fb, hcaly_fl);
  l4_f->SetLineColor(kGreen);

  //Make Fiducial Line
  TLine* l1_fl = new TLine(hcalx_t-2*dysig, hcaly_fl, hcalx_t-2*dysig, hcaly_fl+dx0_p);
  l1_fl->SetLineColor(kMagenta);

  TCanvas *c4 = new TCanvas("c4","HCal Active Area and Fiducial Region",1200,500);
  c4->Divide(2,1);

  c4->cd(1);
  c4->SetLogz();
  hxy->Draw("colz");
  l1_e->Draw();
  l2_e->Draw();
  l3_e->Draw();
  l4_e->Draw();
  l1_aa->Draw();
  l2_aa->Draw();
  l3_aa->Draw();
  l4_aa->Draw();
  l1_f->Draw();
  l2_f->Draw();
  l3_f->Draw();
  l4_f->Draw();
  l1_fl->Draw();

  c4->Update();

  c4->cd(2);
  c4->SetLogz();
  hxy_fcut->Draw("colz");
  l1_e->Draw();
  l2_e->Draw();
  l3_e->Draw();
  l4_e->Draw();
  l1_aa->Draw();
  l2_aa->Draw();
  l3_aa->Draw();
  l4_aa->Draw();
  l1_f->Draw();
  l2_f->Draw();
  l3_f->Draw();
  l4_f->Draw();
  l1_fl->Draw();

  c4->Update();
  c4->Write();

  fout->Write();
  
  cout << endl << "Analysis complete. Outfile located at " << fout_path << endl;

  // Send time efficiency report to console
  st->Stop();
  cout << "CPU time elapsed = " << st->CpuTime() << " s = " << st->CpuTime()/60.0 << " min. Real time = " << st->RealTime() << " s = " << st->RealTime()/60.0 << " min." << endl;

}
