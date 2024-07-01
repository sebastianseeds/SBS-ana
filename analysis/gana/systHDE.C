//sseeds - 5.10.24: Script to extract systematic error from HDE variations across HCAL. Proton and neutron are assumed to be similar for this estimation.
#include <iostream>
#include <vector>
#include <cmath>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TText.h>

// Struct to hold measurement data
struct Measurement {
    double efficiency;
    double error;
};

// Function to calculate the weighted mean
double calculateWeightedMean(const std::vector<Measurement>& measurements) {
    double sumWeights = 0;
    double sumWeightedEff = 0;
    for (const auto& m : measurements) {
        double weight = 1 / (m.error * m.error);
        sumWeights += weight;
        sumWeightedEff += m.efficiency * weight;
    }
    return sumWeightedEff / sumWeights;
}

// Function to calculate the statistical error
double calculateStatisticalError(const std::vector<Measurement>& measurements) {
    double sumWeights = 0;
    for (const auto& m : measurements) {
        double weight = 1 / (m.error * m.error);
        sumWeights += weight;
    }
    return std::sqrt(1 / sumWeights);
}

// Function to calculate the systematic error
double calculateSystematicError(const std::vector<Measurement>& measurements, double weightedMean) {
    double sumDiffsSquared = 0;
    for (const auto& m : measurements) {
        sumDiffsSquared += std::pow(m.efficiency - weightedMean, 2);
    }
    return std::sqrt(sumDiffsSquared / (measurements.size() - 1));
}

// Function to calculate the Z-Score
double calculateZScore(double weightedMean, double mcValue, double combinedError) {
    return (weightedMean - mcValue) / combinedError;
}

// Function to calculate the chi-squared
double calculateChiSquared(const std::vector<Measurement>& measurements, double mcValue) {
    double chi2 = 0.0;
    for (const auto& m : measurements) {
        chi2 += std::pow((m.efficiency - mcValue) / m.error, 2);
    }
    return chi2;
}

// Function to create and fill a canvas with the results
void createCanvas(const std::string& title, const std::vector<std::string>& results) {
    TCanvas *c1 = new TCanvas("c1", "HCal SBS8 proton detection eff. mean and error", 800, 600);
    TPaveText *pt = new TPaveText(0.1, 0.1, 0.9, 0.9);
    
    pt->AddText(title.c_str());
    for (const auto& result : results) {
        pt->AddText(result.c_str());
    }
    
    pt->Draw();
    c1->SaveAs("~/Plots/systematic_error_analysis.png");
}

//MAIN. (no args)
int systHDE() {
    std::vector<Measurement> measurements = {
        {0.936, 0.000696}, // 0% field
        {0.945885, 0.00112}, // 50% field
        {0.944370, 0.000301} // 70% field
    };
    double mcValue = 0.943; // MC fit to scaled efficiency at p_N = 3.17 GeV
    double mcBinomialError = 0.00195; // Binomial error, roughly 14000 ev per MC pN bin 

    double weightedMean = calculateWeightedMean(measurements);
    double statError = calculateStatisticalError(measurements);
    double sysError = calculateSystematicError(measurements, weightedMean);
    double combinedError = std::sqrt(statError * statError + sysError * sysError + mcBinomialError * mcBinomialError);
    //double combinedError = std::sqrt(statError * statError + sysError * sysError);
    

    double zScore = calculateZScore(weightedMean, mcValue, combinedError);
    double chiSquared = calculateChiSquared(measurements, mcValue);

    // Create a vector of results for display
    std::vector<std::string> results;
    results.push_back("Weighted Mean: " + std::to_string(weightedMean));
    results.push_back("Statistical Error: " + std::to_string(statError));
    results.push_back("Systematic Error: " + std::to_string(sysError));
    //results.push_back("MC Binomial Error: " + std::to_string(mcBinomialError));
    results.push_back("Combined Error: " + std::to_string(combinedError));
    results.push_back("Z-Score: " + std::to_string(zScore));
    //results.push_back("Chi-Squared: " + std::to_string(chiSquared));

    // Print results to console
    for (const auto& result : results) {
        std::cout << result << std::endl;
    }

    // Create and save the canvas with results
    createCanvas("HCal SBS8 proton detection eff. mean and error", results);

    return 0;
}
