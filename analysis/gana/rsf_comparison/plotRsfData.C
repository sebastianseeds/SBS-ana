#include <fstream>
#include <iostream>
#include <string>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TMarker.h>

void plotRsfData() {
  std::ifstream infile("rsfdata.txt"); // Adjust the path to your data file as needed
  if (!infile.is_open()) {
    std::cout << "Error: Failed to open the data file." << std::endl;
    return; // Exit the function if the file cannot be opened
  }

  std::string bg_fit_type;
  float Rsf, Rsf_err, chisqr;
  int field, wf, alt;

  std::vector<float> x, y, ex, ey;
  std::vector<int> markerStyle;
  std::vector<int> colors;
  std::vector<bool> filled;
  std::vector<std::string> types;

  while (infile >> bg_fit_type >> Rsf >> Rsf_err >> chisqr >> field >> wf >> alt) {
    x.push_back(chisqr);
    y.push_back(Rsf);
    ex.push_back(0); // No error on chisqr
    ey.push_back(Rsf_err);
    types.push_back(bg_fit_type);

    cout << chisqr << " " << Rsf << " " << Rsf_err << " " << field << " " << wf << " " << alt << endl;

    // Define marker style based on 'field'
    markerStyle.push_back((field == 50) ? 20 : 21); // Circle for 50, square for 30

    // Color code based on 'bg_fit_type'
    if (bg_fit_type == "pol2") colors.push_back(kRed);
    else if (bg_fit_type == "pol4") colors.push_back(kBlue);
    else if (bg_fit_type == "antidy") colors.push_back(kGreen);
    else if (bg_fit_type == "anticoin") colors.push_back(kMagenta);
    else if (bg_fit_type == "sideband") colors.push_back(kCyan);
    else colors.push_back(kBlack); // Default color

    // Fill style based on 'alt'
    filled.push_back((alt == 0));
  }

  cout << endl;

  // Before plotting, calculate the ranges
  float minX = *std::min_element(x.begin(), x.end()) - 0.1; // Subtract a small value for padding
  float maxX = *std::max_element(x.begin(), x.end()) + 0.1; // Add a small value for padding
  float minY = *std::min_element(y.begin(), y.end()) - 0.1; // Subtract a small value for padding
  float maxY = *std::max_element(y.begin(), y.end()) + 0.1; // Add a small value for padding

  TCanvas *c1 = new TCanvas("c1", "Rsf vs chisqr", 200, 10, 700, 500);
  c1->cd();

  // Initialize the legend
  TLegend *leg = new TLegend(0.1, 0.7, 0.3, 0.9);
    
  // Adding legend entries for each bg_fit_type with a colored line
  TLine *line = new TLine(); // Dummy line for legend
  line->SetLineWidth(2);
  line->SetLineColor(kRed);
  leg->AddEntry(line, "pol2", "l");
  TLine *line2 = new TLine(); // Dummy line for legend
  line2->SetLineWidth(2);
  line2->SetLineColor(kBlue);
  leg->AddEntry(line2, "pol4", "l");
  TLine *line3 = new TLine(); // Dummy line for legend
  line3->SetLineWidth(2);
  line3->SetLineColor(kGreen);
  leg->AddEntry(line3, "antidy", "l");
  TLine *line4 = new TLine(); // Dummy line for legend
  line4->SetLineWidth(2);
  line4->SetLineColor(kMagenta);
  leg->AddEntry(line4, "anticoin", "l");
  TLine *line5 = new TLine(); // Dummy line for legend
  line5->SetLineWidth(2);
  line5->SetLineColor(kCyan);
  leg->AddEntry(line5, "sideband", "l");

  // Adding legend entries for 50p, 30p, and 30p alt with symbols
  TMarker *marker50p = new TMarker(0, 0, 20); // Solid circle for 50p
  marker50p->SetMarkerColor(kBlack);
  leg->AddEntry(marker50p, "50p", "p");

  TMarker *marker30p = new TMarker(0, 0, 21); // Solid square for 30p
  marker30p->SetMarkerColor(kBlack);
  leg->AddEntry(marker30p, "30p", "p");

  TMarker *marker30pAlt = new TMarker(0, 0, 25); // Hollow square for 30p alt
  marker30pAlt->SetMarkerColor(kBlack);
  leg->AddEntry(marker30pAlt, "30p alt", "p");
  bool firstGraphDrawn = false;

  for (size_t i = 0; i < x.size(); ++i) {

    cout << x[i] << " " << y[i] << " " << ex[i] << " " << ey[i] << " " << markerStyle[i] << " " << colors[i] << " " << filled[i] << " " << types[i] << endl;

    TGraphErrors *gr = new TGraphErrors(1, &x[i], &y[i], &ex[i], &ey[i]);
    gr->SetMarkerStyle(markerStyle[i]);
    gr->SetMarkerColor(colors[i]);
    if (!filled[i]) gr->SetMarkerStyle(markerStyle[i]+4); // Open marker if alt != 0

    // Add to legend
    //leg->AddEntry(gr, types[i].c_str(), "p");
    
    // Draw the graph
    if (!firstGraphDrawn) {
      gr->Draw("AP");
      gr->GetHistogram()->SetMinimum(minY); // Set y-axis range
      gr->GetHistogram()->SetMaximum(maxY); // Set y-axis range
      gr->GetXaxis()->SetLimits(minX, maxX); // Set x-axis range
      firstGraphDrawn = true;
    } else {
      gr->Draw("P SAME");
    }
  }

  leg->Draw();
  c1->Update();
}
