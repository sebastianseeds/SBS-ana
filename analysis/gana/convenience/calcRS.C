//sseeds
#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>
#include <TCanvas.h>
#include <TLatex.h>

// Function to calculate the Rosenbluth Slope (RS) and FFR
std::pair<std::pair<double, double>, std::pair<double, double>> calculateRosenbluthSlopeAndFFR(
											       double R_epsilon1, double R_epsilon1_err,
											       double R_epsilon2, double R_epsilon2_err,
											       double epsilon1, double epsilon1_err,
											       double epsilon2, double epsilon2_err,
											       double RS_p, double RS_p_err)
{
  // Calculate delta_epsilon
  double delta_epsilon = epsilon1 - epsilon2;
  if (delta_epsilon == 0) {
    throw std::invalid_argument("epsilon1 and epsilon2 cannot be equal");
  }

  // Calculate R_Mott ratios (assuming they are given or calculated separately)
  double R_Mott_epsilon1 = 1.0;  // Placeholder value, replace with actual calculation
  double R_Mott_epsilon2 = 1.0;  // Placeholder value, replace with actual calculation

  // Calculate A and B
  double A = (R_epsilon1 / R_epsilon2);
  double B = (R_Mott_epsilon1 / R_Mott_epsilon2) * ((1.0 + epsilon2 * RS_p) / (1.0 + epsilon1 * RS_p));

  // Calculate RS_n
  double RS_n = (A - B) / (B * delta_epsilon);

  // Calculate the error on RS_n
  // Propagate errors for A and B
  double A_err = A * sqrt(pow(R_epsilon1_err / R_epsilon1, 2) + pow(R_epsilon2_err / R_epsilon2, 2));
  double B_err = B * sqrt(pow(RS_p_err / RS_p, 2) + pow(epsilon1_err * epsilon2 / pow(1.0 + epsilon1 * RS_p, 2), 2) + pow(epsilon2_err * epsilon1 / pow(1.0 + epsilon2 * RS_p, 2), 2));

  // Propagate errors for RS_n
  double RS_n_err = RS_n * sqrt(pow(A_err / (A - B), 2) + pow(B_err / B, 2) + pow(epsilon1_err / delta_epsilon, 2) + pow(epsilon2_err / delta_epsilon, 2));

  // Calculate FFR
  double FFR = sqrt(RS_n);
  double FFR_err = RS_n_err / (2 * FFR); // Error propagation for sqrt

  return std::make_pair(std::make_pair(RS_n, RS_n_err), std::make_pair(FFR, FFR_err));
}

int calcRS() {
  // SBS8 and SBS9 values
  double R_epsilon1 = 0.4016; // sbs8 0.3991 (1 sig), 0.4002 (2 sig), 0.4016 (0.5 sig)
  double R_epsilon1_err = 0.0033; // sbs8 0.0012 (1 sig), 0.0012 (2 sig), 0.0033 (0.5 sig)
  double R_epsilon2 = 0.3867; // sbs9 0.3907 (1 sig), 0.3988 (2 sig), 0.3867 (0.5 sig)
  double R_epsilon2_err = 0.0033; // sbs9 0.0012 (1 sig), 0.0012 (2 sig), 0.0033 (0.5 sig)
  double epsilon1 = 0.803; // sbs8, MC mean
  double epsilon1_err = 0.005; // estimated from data spread
  double epsilon2 = 0.516; // sbs9, MC mean
  double epsilon2_err = 0.005; // estimated from data spread
  double RS_p = 0.107; //from nTPE proposal
  double RS_p_err = 0.010; //from nTPE proposal

  try {
    std::pair<std::pair<double, double>, std::pair<double, double>> result = calculateRosenbluthSlopeAndFFR(
													    R_epsilon1, R_epsilon1_err,
													    R_epsilon2, R_epsilon2_err,
													    epsilon1, epsilon1_err,
													    epsilon2, epsilon2_err,
													    RS_p, RS_p_err);

    double RS_n = result.first.first;
    double RS_n_err = result.first.second;
    double FFR = result.second.first;
    double FFR_err = result.second.second;

    std::cout << "RS_n: " << RS_n << " ± " << RS_n_err << std::endl;
    std::cout << "FFR: " << FFR << " ± " << FFR_err << std::endl;

    // Create a canvas to display the results
    TCanvas *c = new TCanvas("c", "RS and FFR Results", 800, 600);

    // Create TLatex objects to display the results
    TLatex latex;
    latex.SetTextSize(0.03);
    latex.DrawLatex(0.1, 0.8, Form("RS_n: %.4f #pm %.4f", RS_n, RS_n_err));
    latex.DrawLatex(0.1, 0.7, Form("FFR: %.4f #pm %.4f", FFR, FFR_err));

    // Save the canvas as an image or PDF
    c->SaveAs("RS_and_FFR_results.png");
    c->SaveAs("RS_and_FFR_results.pdf");

  } catch (const std::invalid_argument& e) {
    std::cerr << "Error: " << e.what() << std::endl;
  }

  return 0;
}
