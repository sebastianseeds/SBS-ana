//seeds - script to analyze downbending data and compare against MC wapp generator. Making use of A.Rathnayake's analysis root files employing AJRPs downbending optics for BigBite.
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TEllipse.h>
#include <iostream>
#include <TChain.h>

static const double Mp = 0.938272081; // +/- 6E-9 GeV
static const double Mn = 0.939565413; // +/- 6E-9 GeV

static const double Mpip = 0.13957;
static const double Mpiz = 0.13498;

static const double EBeam = 4.0268;

static const int maxTracks = 100;

static const double fitmin = 0.05;
static const double fitmax = 0.4;

// Energy fit function
double skewedGaussian(Double_t *x, Double_t *par) {
    Double_t arg = (x[0] - par[1]) / par[2];
    Double_t skew = par[3] * TMath::Erf(par[4] * arg);
    return par[0] * TMath::Exp(-0.5 * arg * arg) * (1 + skew);
}

// Forward declarations
bool pass_EndPointCut( double m_bbtrpz, double m_bbtrp);
double getMPV( TF1 *fit, double xmin, double xmax );

// Main
void analyze_dwnbending() {
  // Set up display
  gStyle->SetOptStat(0);

  // Open the ROOT file
  TFile *file = TFile::Open("/volatile/halla/sbs/seeds/wapp_analysis/downbending_data/dwnbendana_pass2_sbs9_sbs70p_lh2_2.root");
  if (!file || file->IsZombie()) {
    std::cerr << "File could not be opened!" << std::endl;
    return;
  }

  // Get the tree from the file
  TTree *tree = dynamic_cast<TTree*>(file->Get("T"));
  if (!tree) {
    std::cerr << "Tree not found in the file!" << std::endl;
    file->Close();
    return;
  }

  // Prepare histograms for hcaldx and hcaldy
  TH1D *hdx = new TH1D("hdx", "HCAL dx;dx (m)", 100, -4, 4);
  TH1D *hdy = new TH1D("hdy", "HCAL dy;dy (m)", 100, -4, 4);

  // Draw hcaldx and hcaldy from the tree
  tree->Draw("sbs.neutronrecon.hcaldx>>hdx");
  tree->Draw("sbs.neutronrecon.hcaldy>>hdy");

  // Setting up the TChain
  std::cout << "Adding wapp mc files to chain from /lustre19/expphy/volatile/halla/sbs/seeds/p657sf_sbs9_sbs70p_wapp_heep" << std::endl;

  TChain *chain = new TChain("T");
  chain->Add("/lustre19/expphy/volatile/halla/sbs/seeds/wapp_analysis/p657sf_sbs9_sbs70p_wapp_heep/replayed_*.root");

  // Check if the chain has any entries
  if (chain->GetEntries() == 0) {
    std::cerr << "The chain contains no entries. Check if the files have the correct tree structure." << std::endl;
    return; // Exit if the chain is empty
  }

  TH1D *h_e_noendpoint = new TH1D("h_e_noendpoint", "Histogram of HCal Energy; Energy (GeV); Entries (norm unity)", 100, 0, 1);

  chain->Draw("sbs.hcal.e>>h_e_noendpoint","sbs.hcal.e>0.02");

  TH1D *h_e = new TH1D("h_e", "Histogram of HCal Energy, BB track endpoint mom cut; Energy (GeV); Entries (norm unity)", 100, 0, 1);

  // Loop over entries and return hcal neutron energy with endpoint cuts on BB tracks
  double BBtr_p[maxTracks], BBtr_px[maxTracks], BBtr_py[maxTracks], BBtr_pz[maxTracks];
  double hcale;

  chain->SetBranchStatus("*",0);
  chain->SetBranchStatus( "bb.tr.px", 1 );
  chain->SetBranchStatus( "bb.tr.py", 1 );
  chain->SetBranchStatus( "bb.tr.pz", 1 );
  chain->SetBranchStatus( "bb.tr.p", 1 );
  chain->SetBranchStatus( "sbs.hcal.e", 1 );

  chain->SetBranchAddress( "bb.tr.px", BBtr_px );
  chain->SetBranchAddress( "bb.tr.py", BBtr_py );
  chain->SetBranchAddress( "bb.tr.pz", BBtr_pz );
  chain->SetBranchAddress( "bb.tr.p", BBtr_p );
  chain->SetBranchAddress( "sbs.hcal.e", &hcale );

  Long64_t Nevents = chain->GetEntries();

  cout << "Looping over MC chain and applying endpoint cut..." << endl;
  int passcuts = 0;
  for(Long64_t nevent = 0; nevent<Nevents; nevent++){

    cout << "Running MC event " << nevent << "/" << Nevents << ", with " << passcuts << " passing cuts.\r";
    cout.flush();

    chain->GetEntry(  nevent );

    //Check endpoint momentum
    bool passEndPoint = pass_EndPointCut(BBtr_pz[0],BBtr_p[0]);

    if(!passEndPoint)
      continue;

    if(hcale<0.02)
      continue;

    h_e->Fill(hcale);
    passcuts++;

  }

  cout << endl << "Finished loop. Plotting results..." << endl;

  if(h_e->GetEntries() == 0){
    cout << "No entries in h_e, endpoint cut problem likely." << endl;
    return;
  }

  // Fit histograms with a gaussian to find mean and std dev
  TF1 *fitdx = new TF1("fitdx", "gaus", -4, 4);
  TF1 *fitdy = new TF1("fitdy", "gaus", -4, 4);

  hdx->Fit(fitdx, "R");
  hdy->Fit(fitdy, "R");

  // Extract mean and std dev
  double mean_dx = fitdx->GetParameter(1);
  double sigma_dx = fitdx->GetParameter(2);
  double mean_dy = fitdy->GetParameter(1);
  double sigma_dy = fitdy->GetParameter(2);

  // Plot dx vs dy with an ellipse cut at 3 sigma around the combined peak
  TH2D *hdx_dy = new TH2D("hdx_dy", "HCAL dx vs dy;dy;dx", 100, -2, 2, 100, -4, 2);
  tree->Draw("sbs.neutronrecon.hcaldx:sbs.neutronrecon.hcaldy>>hdx_dy", "", "COLZ");

  // Draw ellipse at 3 sigma
  TEllipse *ellipse = new TEllipse(mean_dy, mean_dx, 3*sigma_dy, 3*sigma_dx);
  ellipse->SetFillStyle(0);
  ellipse->SetLineWidth(2);
  ellipse->SetLineColor(kRed);
  ellipse->Draw("same");

  // Plot hcal energy with 3 sigma cut
  TH1D *he = new TH1D("he", "Histogram of Neutron Recon Energy with 3 Sigma Spot Cut;Energy (GeV);Entries (norm unity)", 100, 0, 1);

  TCut cut = Form("sqrt(pow((sbs.neutronrecon.hcaldx-%f)/%f, 2) + pow((sbs.neutronrecon.hcaldy-%f)/%f, 2)) < 3", mean_dx, sigma_dx, mean_dy, sigma_dy);
  tree->Draw("sbs.hcal.e>>he", cut);
  TH1D *he_clone = (TH1D*)he->Clone("he_clone");
  he_clone->SetTitle("HCal E with 3 Sigma Cut; Energy (GeV)");

  // Display everything on a canvas
  TCanvas *c1 = new TCanvas("c1", "Data Analysis", 800, 600);
  c1->Divide(2,2);
  c1->cd(1); hdx->Draw();
  c1->cd(2); hdy->Draw();
  c1->cd(3); hdx_dy->Draw("COLZ"); ellipse->Draw("same");
  c1->cd(4); he_clone->Draw("E");

  cout << he_clone->GetEntries() << endl;

  // Normalize the histograms
  he->Scale(1.0 / he->Integral());
  h_e->Scale(1.0 / h_e->Integral());
  h_e_noendpoint->Scale(1.0 / h_e_noendpoint->Integral());

  // Draw the histograms on the same canvas for visual comparison
  TCanvas *c2 = new TCanvas("c2", "Comparison of Distributions", 800, 600);
  c2->SetLeftMargin(0.15); // Increase the left margin
  he->SetLineColor(kGreen);
  h_e->SetLineColor(kBlue);
  h_e_noendpoint->SetLineColor(kRed);
  he->Draw("hist");
  h_e->Draw("hist same");
  h_e_noendpoint->Draw("hist same");

  TF1 *skewGausFit_mc = new TF1("skewGausFit_mc", skewedGaussian, fitmin, fitmax, 5);
  skewGausFit_mc->SetParameters(50, 0.2, fitmin, 0.1, 0.2); // Initial parameters: Amplitude, Mean, StdDev, SkewAmplitude, SkewWidth
  skewGausFit_mc->SetParNames("Amplitude", "Mean", "StdDev", "SkewAmplitude", "SkewWidth");
  skewGausFit_mc->SetLineColor(kBlue-2);

  h_e->Fit(skewGausFit_mc, "R");
  //double mpv_mc = skewGausFit_mc->GetParameter(1); // Get the mean value from the fit
  double mpv_mc = getMPV(skewGausFit_mc,fitmin,fitmax);

  // Display the fit
  skewGausFit_mc->Draw("Same");

  TF1 *skewGausFit_data = new TF1("skewGausFit_data", skewedGaussian, fitmin, fitmax, 5);
  skewGausFit_data->SetParameters(50, 0.2, fitmin, 0.1, 0.2); // Initial parameters: Amplitude, Mean, StdDev, SkewAmplitude, SkewWidth
  skewGausFit_data->SetParNames("Amplitude", "Mean", "StdDev", "SkewAmplitude", "SkewWidth");
  skewGausFit_data->SetLineColor(kGreen-2);

  he->Fit(skewGausFit_data, "R");
  //double mpv_data = skewGausFit_data->GetParameter(1); // Get the mean value from the fit
  double mpv_data = getMPV(skewGausFit_data,fitmin,fitmax);

  // Display the fit
  skewGausFit_data->Draw("Same");

  TF1 *skewGausFit_noep = new TF1("skewGausFit_noep", skewedGaussian, fitmin, fitmax, 5);
  skewGausFit_noep->SetParameters(50, 0.2, fitmin, 0.1, 0.2); // Initial parameters: Amplitude, Mean, StdDev, SkewAmplitude, SkewWidth
  skewGausFit_noep->SetParNames("Amplitude", "Mean", "StdDev", "SkewAmplitude", "SkewWidth");
  skewGausFit_noep->SetLineColor(kRed-2);

  h_e_noendpoint->Fit(skewGausFit_noep, "R");
  //double mpv_noep = skewGausFit_noep->GetParameter(1); // Get the mean value from the fit
  double mpv_noep = getMPV(skewGausFit_noep,fitmin,fitmax);

  // Display the fit
  skewGausFit_noep->Draw("Same");

  // Adding legend
  TLegend *leg = new TLegend(0.4, 0.7, 0.89, 0.89);
  leg->AddEntry(he, Form("Downbending Track Data; Nev = %.0f; MPV = %.2f",he->GetEntries(),mpv_data), "l");
  leg->AddEntry(h_e, Form("WAPP MC Neutrons, endpoint BB mom cut; Nev = %.0f, MPV = %.2f",h_e->GetEntries(),mpv_mc), "l");
  leg->AddEntry(h_e_noendpoint, Form("WAPP MC Neutrons; Nev = %.0f, MPV = %.2f",h_e_noendpoint->GetEntries(),mpv_noep), "l");
  leg->Draw();

  // Update the canvas to display changes
  c2->Update();

}

//Adapted from ANR https://github.com/anuruddhadilshan/GMn-adr-ana/blob/3845b2523a0e2f30fb84cdccfb77c6dc4105e8a7/HCal-NDE-ana/eventclass.h#L152C2-L200C3
//Endpoint momentum cut to remove multi-pion channel response
bool pass_EndPointCut( double m_bbtrpz, double m_bbtrp) {
  double E_gamma = EBeam;

  // Mandelstam s variable for the photon(with beam energy) + proton(at rest = target) pair.
  double s = std::pow(Mp,2) + 2*E_gamma*Mp;

  // Establish the CM frame.
  double v_CMframe = E_gamma/(E_gamma + Mp); // Velocity of the CM frame w.r.t the lab frame (in the direction of the beam i.e. +Z).
  double beta = v_CMframe;
  double gamm_fac = 1/std::sqrt(1 - std::pow(beta,2));
  double E_tot_CM = std::sqrt(s);

  // photon + p --> pi+ + n
  // Calculating pi+ production thresholds (maximum) in terms of momentum and energy for the given beam energy incident on a stationary hydrogen (proton) target.
  double p_piplus_CM_1pi = std::sqrt( std::pow((s + std::pow(Mn,2) - std::pow(Mp,2)),2)/(4*s) - std:: pow(Mn,2) );
  // Transform the pi+ momentum into lab frame.
  double E_piplus_CM_1pi = std::sqrt(std::pow(p_piplus_CM_1pi,2) + std::pow(Mpip,2));
  double p_piplus_lab_1pi = gamm_fac*(p_piplus_CM_1pi + v_CMframe*E_piplus_CM_1pi); // Maximum possible pi+ momentum in the lab frame.
  double E_piplus_lab_1pi = gamm_fac*(E_piplus_CM_1pi + v_CMframe*p_piplus_CM_1pi);

  // photon + p --> pi+ + pi0 + n 
  // Calculating pi+ production thresholds (maximum) in terms of momentum and energy for the given beam energy incident on a stationary hydrogen (proton) target.
  double m_T = Mn + Mpiz; // Treat the pi0 and n as a composite system.
  double p_piplus_CM_2pi = std::sqrt( std::pow((s + std::pow(m_T,2) - std::pow(Mpip,2)),2)/(4*s) - std::pow(m_T,2) ); // Maximum pi+ momentum possible for the given beam energy in the CM frame.
  // Transform the pi+ momentum into lab frame.	
  double E_piplus_CM_2pi = std::sqrt(std::pow(p_piplus_CM_2pi,2) + std::pow(Mpip,2));
  double p_piplus_lab_2pi = gamm_fac*(p_piplus_CM_2pi + v_CMframe*E_piplus_CM_2pi); // Maximum possible pi+ momentum in the lab frame.
  double E_piplus_lab_2pi = gamm_fac*(E_piplus_CM_2pi + v_CMframe*p_piplus_CM_2pi);

  // Now we want to extend the calculation to the case where we have non-zero polar angle for theta/theta'
  double m_bbpoltgt = acos(m_bbtrpz/m_bbtrp); //Polar scattering angle of pi+ in radians.
  double tan_theta = tan(m_bbpoltgt); 

  double cos_thetaprime_1pi_pos = ( -v_CMframe*E_piplus_CM_1pi*std::pow(gamm_fac*tan_theta,2) + std::sqrt(std::pow(gamm_fac*tan_theta,2)*(std::pow(p_piplus_CM_1pi,2) - std::pow(v_CMframe*E_piplus_CM_1pi,2)) + std::pow(p_piplus_CM_1pi,2))) / ( p_piplus_CM_1pi*(std::pow(gamm_fac*tan_theta,2) + 1) );
  double p_piplus_lab_1pi_theta_pos = std::sqrt( (std::pow(gamm_fac,2)-1)*std::pow(p_piplus_CM_1pi,2)*std::pow(cos_thetaprime_1pi_pos,2) + 2*std::pow(gamm_fac,2)*p_piplus_CM_1pi*v_CMframe*E_piplus_CM_1pi*cos_thetaprime_1pi_pos + std::pow(p_piplus_CM_1pi,2)+std::pow(gamm_fac*v_CMframe*E_piplus_CM_1pi,2));
  double p_piplus_lab_1pi_theta_max = p_piplus_lab_1pi_theta_pos;  

  double cos_thetaprime_2pi_pos = ( -v_CMframe*E_piplus_CM_2pi*std::pow(gamm_fac*tan_theta,2) + std::sqrt(std::pow(gamm_fac*tan_theta,2)*(std::pow(p_piplus_CM_2pi,2) - std::pow(v_CMframe*E_piplus_CM_2pi,2)) + std::pow(p_piplus_CM_2pi,2))) / ( p_piplus_CM_2pi*(std::pow(gamm_fac*tan_theta,2) + 1) );
  double p_piplus_lab_2pi_theta_pos = std::sqrt( (std::pow(gamm_fac,2)-1)*std::pow(p_piplus_CM_2pi,2)*std::pow(cos_thetaprime_2pi_pos,2) + 2*std::pow(gamm_fac,2)*p_piplus_CM_2pi*v_CMframe*E_piplus_CM_2pi*cos_thetaprime_2pi_pos + std::pow(p_piplus_CM_2pi,2)+std::pow(gamm_fac*v_CMframe*E_piplus_CM_2pi,2));
  double p_piplus_lab_2pi_theta_max = p_piplus_lab_2pi_theta_pos;
  double piplus_momentum_limit = p_piplus_lab_2pi_theta_max*(1 + 1.5/100.0) - 0.2; //Substract 200 MeV/c from the momentum threshold give some room. 

  if ( m_bbtrp < piplus_momentum_limit ) return false; // End-point cut.

  if ( m_bbtrp > p_piplus_lab_1pi_theta_max+0.2 ) return false; // Max possible momentum for single pion production.

  return true; 
}

double getMPV( TF1 *fit, double xmin, double xmax ){

  double mpv = fit->GetMaximumX(xmin, xmax);

  return mpv;

}
