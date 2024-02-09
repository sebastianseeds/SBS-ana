//sseeds 2.2.24 GMn/nTPE plot of Q^2 coverage.

#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TH2D.h>
#include <TMarker.h>
#include <TBox.h>
#include <TLatex.h>

void plot_epsilon_vs_Q2() {
  
  gStyle->SetOptStat(0);

  // Kinematics data
  const int nPoints = 6;
  double Q2[nPoints] = {3.0, 4.5, 4.5, 7.4, 9.9, 13.5};
  double E[nPoints] = {3.73, 5.98, 4.03, 5.98, 7.93, 9.89};
  double epsilon[nPoints] = {0.719, 0.798, 0.514, 0.466, 0.497, 0.412};
  int kinematics[nPoints] = {4, 8, 9, 14, 7, 11};
  int colors[nPoints] = {kRed, kBlue, kGreen, kMagenta, kCyan, kOrange};

  // Create canvas
  TCanvas *c1 = new TCanvas("c1", "Epsilon vs Q^2", 800, 600);

  // Create and draw the TH2D histogram
  TH2D *h2 = new TH2D("h2", "#it{#epsilon}\ vsQ^{2};Q^{2} GeV^{2};#it{#epsilon}", 100, 0, 15, 100, 0, 1);
  h2->Draw();

  // Add grid lines to the canvas
  c1->SetGridx(1);
  c1->SetGridy(1);

  // Create a TGraph for each kinematic point
  TGraph *graphs[nPoints];
  TLegend *leg = new TLegend(0.6, 0.6, 0.89, 0.89);
  leg->SetBorderSize(0); // Remove border of legend
  //leg->SetFillStyle(0); // Transparent legend background
  //leg->SetFillColor(0);

  for (int i = 0; i < nPoints; ++i) {
    graphs[i] = new TGraph(1, &Q2[i], &epsilon[i]);
    graphs[i]->SetMarkerColor(colors[i]);
    graphs[i]->SetMarkerStyle(21);
    graphs[i]->SetMarkerSize(0.8); // Reduced marker size by 1/4
    graphs[i]->Draw("P same");

    // Add points to legend
    leg->AddEntry(graphs[i], Form("Kine %d, E = %0.2f GeV", kinematics[i], E[i]), "p");
  }

  // Shaded region for nTPE
  TBox *box = new TBox(4, 0, 5, 1);
  box->SetFillColorAlpha(kGray, 0.5);
  box->Draw();
  leg->AddEntry(box, "nTPE", "f");

  // Draw legend
  leg->Draw();

  // Update canvas
  c1->Update();

  // Create a folder named 'plots' in the current directory if it doesn't exist
  gSystem->mkdir("plots", true);

  // Save the canvas as a PDF in the 'plots' folder with the script's name
  c1->SaveAs("plots/plot_epsilon_vs_Q2.pdf");

}
