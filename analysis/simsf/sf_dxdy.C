//sseeds - short script to test various field scales for post-pass2 g4sbs data-simulation peak matching (sbs9)
//11.14.23 - added fiducial cut and fiducial region check

#include <TSystem.h>
#include <TChain.h>
#include <TString.h>
#include "TFile.h"
#include "TTree.h"
#include <TNtuple.h>
#include "TCanvas.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include "TMath.h"
#include "TH1F.h"
#include <TH2.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>
#include <TPaveLabel.h>
#include <TProfile.h>
#include <TPolyLine.h>
#include <TObjArray.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <math.h>
#include <stack>
#include "TLorentzVector.h"
#include "TCut.h"
#include "TLatex.h"
#include <iomanip>
#include "../../include/gmn.h"

//global variables
const int maxTracks = 16; 
const double PI = TMath::Pi();
const double M_e = 0.0051;
const double M_p = 0.938272;
const double M_n = 0.939565;
const double rough_sig_x = 0.182;
const double rough_sig_y = 0.328;

const double hcalheight = 0.;

//Hack to get input fine scale to read the right file without changing all the extensions
string intToFourDigitString(int number) {
    if (number >= 0 && number <= 9999) {
        std::ostringstream oss;
        oss << std::setw(4) << std::setfill('0') << number;
        return oss.str();
    } else {
        return "Invalid input";
    }
}

//MAIN. scale=nominal kinematic magnetic field (percent), course=load course MC separations, fine_scale=MC field scale under study
void sf_dxdy( int scale = 100, bool course = false, int fine_scale = 0 ){//main 

  TChain *C = new TChain("T");

  string type = "fine";
  if(course)
    type = "course";

  string inputfilename = Form("/volatile/halla/sbs/seeds/scale_field_sbs9/course_data/replayed_sf%d*.root",scale);
  double scale_field = (double)scale;
  string fine_string = intToFourDigitString(fine_scale);

  if( course ){
    inputfilename = Form("/volatile/halla/sbs/seeds/scale_field_sbs9/fine_data/replayed_sf%dp%s*.root",scale,fine_string.c_str());
    double fine_scale_field = (double)fine_scale/10000;
    scale_field = (double)scale + fine_scale_field;
  }

  cout << "Reading in " << type << " scale field " << scale_field << " data file at: " << inputfilename << endl;

  C->Add(inputfilename.c_str());

  long Nentries = C->GetEntries();

  if( Nentries<1 ){
    cout << "ERROR: File is empty or missing." << endl;
    return;
  }

  TCut globalcut = "bb.ps.e>0.15&&abs(bb.tr.vz[0])<0.075&&sbs.hcal.nclus>0&&bb.sh.nclus>0&&bb.tr.p[0]>1.5";

  string outputfilename = Form("/volatile/halla/sbs/seeds/scale_field_sbs9/plots/sf%d_dxdy_out.root",scale);

  if( course )
    outputfilename = Form("/volatile/halla/sbs/seeds/scale_field_sbs9/plots/sf%dp%s_dxdy_out.root",scale,fine_string.c_str());

  //string outputfilename = "test.root";

  //SBS 9
  
  double E_e = 4.0268;
  double BB_d = 1.550;
  double BB_th = 49.0*TMath::DegToRad();
  double HCal_d = 11.0;
  double HCal_th = 22.0*TMath::DegToRad();
  double W2_mean = 0.91;
  double W2_sig = 0.4;
  double dy_mean = 0.0;
  double dy_sig = 0.3;
  
  TEventList *elist = new TEventList("elist","Elastic Event List");
  C->Draw(">>elist",globalcut);

  double BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  double BBtr_vz[maxTracks];
  double BBtr_n, BBps_x, BBps_y, BBps_e, BBsh_x, BBsh_y, BBsh_e;
  double HCALx, HCALy, HCALe, ekineW2;

  double nucleon;

  C ->SetBranchStatus("*",0);

  C->SetBranchStatus("sbs.hcal.x",1);
  C->SetBranchStatus("sbs.hcal.y",1);
  C->SetBranchStatus("sbs.hcal.e",1);
  C->SetBranchStatus("bb.tr.n",1);
  C->SetBranchStatus("bb.tr.px",1);
  C->SetBranchStatus("bb.tr.py",1);
  C->SetBranchStatus("bb.tr.pz",1);
  C->SetBranchStatus("bb.tr.p",1);
  C->SetBranchStatus("bb.tr.vz",1);
  C->SetBranchStatus("bb.ps.e",1);
  C->SetBranchStatus("bb.ps.x",1);
  C->SetBranchStatus("bb.ps.y",1);
  C->SetBranchStatus("bb.sh.e",1);
  C->SetBranchStatus("bb.sh.x",1);
  C->SetBranchStatus("bb.sh.y",1);
  C->SetBranchStatus("bb.gem.track.nhits",1);
  C->SetBranchStatus("bb.etot_over_p",1);
  C->SetBranchStatus("MC.mc_nucl",1);

  C->SetBranchAddress("sbs.hcal.x", &HCALx);
  C->SetBranchAddress("sbs.hcal.y", &HCALy);
  C->SetBranchAddress("sbs.hcal.e", &HCALe);
  C->SetBranchAddress("bb.tr.n", &BBtr_n);
  C->SetBranchAddress("bb.tr.px", BBtr_px);
  C->SetBranchAddress("bb.tr.py", BBtr_py);
  C->SetBranchAddress("bb.tr.pz", BBtr_pz);
  C->SetBranchAddress("bb.tr.p", BBtr_p);
  C->SetBranchAddress("bb.tr.vz", BBtr_vz);
  C->SetBranchAddress("bb.ps.e", &BBps_e);
  C->SetBranchAddress("bb.ps.x", &BBps_x);
  C->SetBranchAddress("bb.ps.y", &BBps_y);
  C->SetBranchAddress("bb.sh.e", &BBsh_e);
  C->SetBranchAddress("bb.sh.x", &BBsh_x);
  C->SetBranchAddress("bb.sh.y", &BBsh_y);
  C->SetBranchAddress("MC.mc_nucl", &nucleon);

  TFile *fout = new TFile(outputfilename.c_str(),"RECREATE");

  TH2D *hxy_HCAL = new TH2D("hxy_HCAL", " ; y_{HCAL} (m); x_{HCAL} (m)  ", 250, -1.25, 1.25, 250,-4,3);
  TH2D *hxy_HCAL_fid = new TH2D("hxy_HCAL_fid", " ; y_{HCAL} (m); x_{HCAL} (m)  ", 250, -1.25, 1.25, 250,-4,3);

  TH1D *hdx_HCAL = new TH1D("hdx_HCAL", " ; x_{HCAL} - x_{exp} (m)  ", 250,-4,3);
  TH1D *hdx_HCAL_money = new TH1D("hdx_HCAL_money", " ; x_{HCAL} - x_{exp} (m)  ", 250,-4,3);
  TH1D *hdy_HCAL = new TH1D("hdy_HCAL","; y_{HCAL} - y_{exp} (m);", 250, -1.25,1.25);
  TH1D *hy = new TH1D("hy","; y_{HCAL} (m);", 250, -1.25,1.25);
  TH1D *hyexp = new TH1D("hyexp","; y_{exp} (m);", 250, -1.25,1.25);
  TH1D *hW2 = new TH1D("hW2", " ;GeV2  ", 100,0,5);
  TH1D *hQ2 = new TH1D("hQ2", " ;GeV2  ", 150,0,15);

  TH1D *hdx_HCAL_p = new TH1D("hdx_HCAL_p", " ; x_{HCAL} - x_{exp} (m)  ", 250,-4,3);
  TH1D *hdx_HCAL_n = new TH1D("hdx_HCAL_n", " ; x_{HCAL} - x_{exp} (m)  ", 250,-4,3);

  TH1D *hdx_HCAL_fid = new TH1D("hdx_HCAL_fid", " ; x_{HCAL} - x_{exp} (m)  ", 250,-4,3);
  TH1D *hdy_HCAL_fid = new TH1D("hdy_HCAL_fid","; y_{HCAL} - y_{exp} (m);", 250, -1.25,1.25);
  TH1D *hdx_HCAL_p_fid = new TH1D("hdx_HCAL_p_fid", " ; x_{HCAL} - x_{exp} (m)  ", 250,-4,3);
  TH1D *hdx_HCAL_n_fid = new TH1D("hdx_HCAL_n_fid", " ; x_{HCAL} - x_{exp} (m)  ", 250,-4,3);

  Long64_t elastic = 0;

  Long64_t Nevents = elist->GetN();

  vector<Double_t> hcalaa = cut::hcalaa_mc(1,1); //exclude one block on periphery
  
  cout << Nevents << " passed global cut." << endl;

  for( Long64_t nevent = 1; nevent <Nevents; nevent++){

    cout << " Entry = " << nevent << "/" << Nevents << " (elastics = " << elastic << ")" << "\r";
    cout.flush();
    
    C->GetEntry(elist ->GetEntry(nevent)); 
    
    double etheta = acos( BBtr_pz[0]/BBtr_p[0]);
    double ephi = atan2(BBtr_py[0],BBtr_px[0]);
    
    TVector3 vertex(0,0,BBtr_vz[0]);
    TLorentzVector Pbeam(0,0,E_e,E_e);
    TLorentzVector kprime(BBtr_px[0], BBtr_py[0], BBtr_pz[0], BBtr_p[0]);
    TLorentzVector Ptarg(0, 0, 0, M_p);

    TLorentzVector q = Pbeam - kprime;
    TLorentzVector PgammaN = Ptarg + q; //should go through and write this out. Momentum of virtual photon

    double pel = E_e/ (1. +E_e/M_p*(1.-cos(etheta)));//momentum of elastically scattered electron 
    double nu = E_e -BBtr_p[0]; //kinetic energy of the elastically scattered electron 
    double pp = sqrt(pow(nu,2)+2 *M_p*nu); 
    double phinucleon = ephi + PI; //coplanar 
    double thetanucleon = acos((E_e - BBtr_pz[0])/pp);

    TVector3 pNhat( sin(thetanucleon)*cos(phinucleon), sin(thetanucleon)*sin(phinucleon), cos(thetanucleon) );
    
    TVector3 HCAL_zaxis(sin(-HCal_th),0,cos(-HCal_th));
    TVector3 HCAL_xaxis(0,-1,0);
    TVector3 HCAL_yaxis = HCAL_zaxis.Cross(HCAL_xaxis).Unit();

    TVector3 HCAL_origin = HCal_d * HCAL_zaxis + hcalheight * HCAL_xaxis;

    double sintersect = ( HCAL_origin - vertex ).Dot( HCAL_zaxis) / (pNhat.Dot( HCAL_zaxis ) );
    
    TVector3 HCAL_intersect = vertex + sintersect * pNhat; 

    double yexpect_HCAL = ( HCAL_intersect - HCAL_origin ).Dot( HCAL_yaxis );
    double xexpect_HCAL = ( HCAL_intersect - HCAL_origin ).Dot( HCAL_xaxis );

    hy->Fill(HCALy);
    hyexp->Fill(yexpect_HCAL);
    
    double Q2 = 2*E_e *BBtr_p[0]*( 1-cos(etheta));

    double W2 = pow( M_p,2 )+2*M_p*nu-Q2;
   
    hW2->Fill(W2);
    hQ2->Fill(Q2);
    
    double dx = HCALx - xexpect_HCAL;
    double dy = HCALy - yexpect_HCAL;

    //H-arm fiducial cuts
    double est_peak_loc = -0.01308*scale-0.00157;
 
    //cout << est_peak_loc << endl;

    bool hcalON = cut::hcalaaON(HCALx,HCALy,hcalaa);
    vector<Double_t> fid = cut::hcalfid(rough_sig_x,rough_sig_y,hcalaa);
    bool passed_fid = cut::hcalfidIN(xexpect_HCAL,yexpect_HCAL,est_peak_loc,fid);
    
    bool on_and_fid = hcalON && passed_fid;

    hdy_HCAL->Fill(dy);
    hdx_HCAL->Fill(dx);
    hxy_HCAL->Fill(HCALy, HCALx);
    if(nucleon==0)
      hdx_HCAL_n->Fill(dx);
    if(nucleon==1)
      hdx_HCAL_p->Fill(dx);

    if( on_and_fid ){
      hdy_HCAL_fid->Fill(dy);
      hdx_HCAL_fid->Fill(dx);
      hxy_HCAL_fid->Fill(HCALy, HCALx);
      if(nucleon==0)
	hdx_HCAL_n_fid->Fill(dx);
      if(nucleon==1)
	hdx_HCAL_p_fid->Fill(dx);
    }

    if(abs(dy-dy_mean)>dy_sig) continue;

    hdx_HCAL_money->Fill(dx);
    
    elastic++;

  }

  cout << endl << endl << "Total elastics = " << elastic << endl;

  cout << "Analysis complete. Out file written to " << outputfilename << endl;

  fout->Write();

}// end main







