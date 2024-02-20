#include <TH1D.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>

void test() {
  std::vector<double> A = {2.1, 4.5, 6.7, 4.6, 1.1}; // Bin contents
  std::vector<double> B = {1.2, 2.4, 3.6, 4.8, 6.0}; // Bin centers

  if (B.size() != A.size()) {
    std::cerr << "Error: The size of vectors A and B must be the same." << std::endl;
    return;
  }

  // Calculate bin edges
  std::vector<double> binEdges;
  binEdges.push_back(B.front() - (B[1] - B[0]) / 2); // First edge
  for (size_t i = 0; i < B.size() - 1; ++i) {
    binEdges.push_back((B[i] + B[i + 1]) / 2); // Middle edges
  }
  binEdges.push_back(B.back() + (B.back() - B[B.size() - 2]) / 2); // Last edge

  // Create a TH1D histogram
  TH1D* hist = new TH1D("hist", "Histogram with Custom Bin Centers", B.size(), &binEdges[0]);

  // Fill the histogram
  for (size_t i = 0; i < A.size(); ++i) {
    hist->SetBinContent(i + 1, A[i]);
  }

  // Draw the histogram
  TCanvas* c = new TCanvas("c", "Canvas", 800, 600);
  hist->Draw();

  // Optionally save the canvas or histogram here

}
