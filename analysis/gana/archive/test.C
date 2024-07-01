#include "TROOT.h"
#include "TFile.h"
#include "TH2D.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TMarker.h"

void createGaussianTH2D() {
    // Define the number of bins and range for the histogram
    const Int_t nBinsX = 1000;
    const Int_t nBinsY = 1000;
    const Double_t xMin = -5.0;
    const Double_t xMax = 5.0;
    const Double_t yMin = -5.0;
    const Double_t yMax = 5.0;

    // Create a 2D histogram
    TH2D *h2 = new TH2D("h2", "2D Gaussian Distribution;X;Y", nBinsX, xMin, xMax, nBinsY, yMin, yMax);

    // Define the number of entries
    const Int_t nEntries = 10000;

    // Create random number generators for Gaussian distribution
    TRandom *rand = new TRandom();

    // Fill the histogram with Gaussian distributed values
    for (Int_t i = 0; i < nEntries; ++i) {
        Double_t x = rand->Gaus(-2.9, 0.01); // mean = 0, sigma = 1
        Double_t y = rand->Gaus(-2.9, 0.01); // mean = 0, sigma = 1
        h2->Fill(x, y);
    }

    // Draw the histogram
    TCanvas *c1 = new TCanvas("c1", "2D Gaussian Distribution", 800, 600);
    h2->Draw("COLZ");

    // Define the square properties
    Double_t squareSide = 6.0;
    Double_t squareHalfSide = squareSide / 2.0;
    //Double_t squareHalfSide = squareSide;

    // Draw the square (removing the bottom right corner)
    TLine *line1 = new TLine(-squareHalfSide, -squareHalfSide, squareHalfSide, -squareHalfSide);
    TLine *line2 = new TLine(squareHalfSide, -squareHalfSide, squareHalfSide, 0);
    TLine *line3 = new TLine(squareHalfSide, squareHalfSide, -squareHalfSide, squareHalfSide);
    TLine *line4 = new TLine(-squareHalfSide, squareHalfSide, -squareHalfSide, -squareHalfSide);

    line1->SetLineColor(kRed);
    line2->SetLineColor(kRed);
    line3->SetLineColor(kRed);
    line4->SetLineColor(kRed);

    // line1->SetLineWidth(50);
    // line2->SetLineWidth(50);
    // line3->SetLineWidth(50);
    // line4->SetLineWidth(50);

    line1->SetLineWidth(2);
    line2->SetLineWidth(2);
    line3->SetLineWidth(2);
    line4->SetLineWidth(2);

    line1->Draw();
    line2->Draw();
    line3->Draw();
    line4->Draw();

    // Add a red dot at the center of the Gaussian
    TMarker *centerMarker = new TMarker(3, 3, 20); // 20 is the marker style (full circle)
    centerMarker->SetMarkerColor(kRed);
    centerMarker->SetMarkerSize(2); // Adjust the size as needed
    centerMarker->Draw();

    // Save the histogram to a ROOT file
    TFile *outputFile = new TFile("GaussianTH2D_with_Square_and_Center.root", "RECREATE");
    h2->Write();
    outputFile->Close();


}

int test() {
    createGaussianTH2D();
    return 0;
}
