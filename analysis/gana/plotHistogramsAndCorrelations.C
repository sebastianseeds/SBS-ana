#include <TFile.h>
#include <TCanvas.h>
#include <TH2D.h>
#include <TKey.h>
#include <TIterator.h>
#include <TClass.h>
#include <TLegend.h>
#include <TPRegexp.h>
#include <TObjArray.h>
#include <TObjString.h>

void plotHistogramsAndCorrelations(int kine=4, int mag=30, int pass=1){

  //"/volatile/halla/sbs/seeds/gmn_analysis/dx_correlations_sbs4_mag30_pass1.root"
  
  //set up files and paths
  std::string outdir_path = gSystem->Getenv("OUT_DIR");

  std::string fin_path = outdir_path + Form("/gmn_analysis/dx_correlations_sbs%d_mag%d_pass%d.root",kine,mag,pass);

  TFile* file = new TFile(fin_path.c_str());
  // Check if the file is open and good to go
  if (!file || !file->IsOpen()) {
    std::cerr << "File is not open!" << std::endl;
    return;
  }

  // Prepare an iterator to go through all objects in the file
  TIter next(file->GetListOfKeys());
  TKey *key;
  TCanvas *canvas = nullptr;
  int histCount = 0;
  const int histsPerCanvas = 4;

  // Get vectors to store correlations
  std::vector<double> correlations;
  std::vector<double> newCorrelations;
  std::vector<std::string> names;

  //Make map to store restricted ranges
  std::map<std::string, std::pair<double, double>> axisLimits;

  //Just get limits manually
  axisLimits["sbs.hcal.e"]={0.03,2.0};
  axisLimits["bb.tr.vz"]={-0.075,0.075};
  axisLimits["dy"]={-0.25,0.25};
  axisLimits["bb.ps.e"]={0.2,5};
  axisLimits["W2"]={0.58,1.18};
  axisLimits["bb.gem.track.nhits"]={3,6};
  axisLimits["sbs.hcal.nclus"]={0,20};
  axisLimits["bb.sh.nclus"]={0,20};

  // Loop over all keys in the file
  while ((key = (TKey*)next())) {
    // Check for TH2D class
    if (strcmp(key->GetClassName(), "TH2D") != 0) continue;

    // Get the histogram
    TH2D *hist = (TH2D*)key->ReadObj();
    // Check for histograms starting with "coorhist"
    if (TString(hist->GetName()).BeginsWith("coorhist")) {
      // Create a new canvas every 4 histograms
      if (histCount % histsPerCanvas == 0) {
	// If there's an existing canvas, update it
	if (canvas) canvas->Update();
	// Create a new canvas
	canvas = new TCanvas(Form("c%d",histCount), Form("Histograms and Correlations Hist %d",histCount), 1200, 800);
	canvas->Divide(2, 2); // Divide the canvas into 4 pads
      }
            
      // Increment the histogram count
      histCount++;

      // Draw the histogram on the appropriate pad
      canvas->cd(histCount % histsPerCanvas == 0 ? histsPerCanvas : histCount % histsPerCanvas);

      // Determine x and y axis limits from the histogram name or title if necessary
      double xmin = axisLimits[hist->GetXaxis()->GetTitle()].first;
      double xmax = axisLimits[hist->GetXaxis()->GetTitle()].second;
      double ymin = axisLimits[hist->GetYaxis()->GetTitle()].first;
      double ymax = axisLimits[hist->GetYaxis()->GetTitle()].second;

      // Clone and restrict the range of the histogram
      TH2D* newHist = (TH2D*)hist->Clone(TString::Format("%s_restricted", hist->GetName()));
      newHist->GetXaxis()->SetRangeUser(xmin, xmax);
      newHist->GetYaxis()->SetRangeUser(ymin, ymax);

      hist->Draw("COLZ"); // Draw the histogram in color

      // Retrieve the actual axes limits of the histogram
      double histXmin = hist->GetXaxis()->GetXmin();
      double histXmax = hist->GetXaxis()->GetXmax();
      double histYmin = hist->GetYaxis()->GetXmin();
      double histYmax = hist->GetYaxis()->GetXmax();

      // Adjust the box coordinates to not exceed the histogram boundaries
      double boxXmin = std::max(xmin, histXmin);
      double boxXmax = std::min(xmax, histXmax);
      double boxYmin = std::max(ymin, histYmin);
      double boxYmax = std::min(ymax, histYmax);

      // Draw lines to form a box around the cut region within the histogram boundaries
      TLine *lineTop = new TLine(boxXmin, boxYmax, boxXmax, boxYmax);
      lineTop->SetLineColor(kRed);
      lineTop->Draw("same");

      TLine *lineBottom = new TLine(boxXmin, boxYmin, boxXmax, boxYmin);
      lineBottom->SetLineColor(kRed);
      lineBottom->Draw("same");

      TLine *lineLeft = new TLine(boxXmin, boxYmin, boxXmin, boxYmax);
      lineLeft->SetLineColor(kRed);
      lineLeft->Draw("same");

      TLine *lineRight = new TLine(boxXmax, boxYmin, boxXmax, boxYmax);
      lineRight->SetLineColor(kRed);
      lineRight->Draw("same");

      // Calculate the correlation coefficient
      Double_t correlation = hist->GetCorrelationFactor();
      Double_t newCorrelation = newHist->GetCorrelationFactor();

      correlations.push_back(correlation);
      newCorrelations.push_back(newCorrelation);
      names.push_back(hist->GetName());

      // Create a legend and add the correlation info
      TLegend *leg = new TLegend(0.1, 0.7, 0.5, 0.9);
      TString corrInfo;
      TString newCorrInfo;
      corrInfo.Form("Correlation: %.2f", correlation);
      leg->AddEntry((TObject*)0, corrInfo, "");
      newCorrInfo.Form("Restricted Correlation: %.2f", newCorrelation);
      leg->AddEntry((TObject*)0, newCorrInfo, "");
      leg->Draw();
    }
  }

  // Update the last canvas
  if (canvas) canvas->Update();

  // Assuming you've collected correlations, newCorrelations, and names already

  // Check if we have data to plot
  if (!correlations.empty() && correlations.size() == newCorrelations.size()) {
    // Create a new canvas for the summary plot
    TCanvas* summaryCanvas = new TCanvas("summaryCanvas", "Correlation Factors", 1200, 800);
    summaryCanvas->cd();

    const int n = correlations.size();
    TGraph* graphCorrelations = new TGraph(n); // Original correlations
    TGraph* graphNewCorrelations = new TGraph(n); // New (restricted) correlations

    for (int i = 0; i < n; ++i) {
      graphCorrelations->SetPoint(i, i + 1, correlations[i]);
      graphNewCorrelations->SetPoint(i, i + 1, newCorrelations[i]);
    }

    // Customize the original correlations graph
    graphCorrelations->SetMarkerStyle(21); // Square marker
    graphCorrelations->SetMarkerColor(kBlue); // Blue color for original correlations
    graphCorrelations->SetTitle("Correlation Factors;;Correlation Factor");

    // Customize the new correlations graph
    graphNewCorrelations->SetMarkerStyle(22); // Triangle marker
    graphNewCorrelations->SetMarkerColor(kRed); // Red color for new correlations


    graphCorrelations->GetHistogram()->SetMaximum(0.25); // Set the maximum y value to 0.25

    // Drawing both graphs
    graphCorrelations->Draw("AP");
    graphNewCorrelations->Draw("P same"); // Draw on the same canvas without axes

    // Setting up the legend
    TLegend* leg = new TLegend(0.1, 0.7, 0.3, 0.9);
    leg->AddEntry(graphCorrelations, "Original Correlation", "p");
    leg->AddEntry(graphNewCorrelations, "Restricted Correlation", "p");
    leg->Draw();

    // Adjusting x-axis labels
    for (int i = 0; i < n; ++i) {
      TString label = names[i].c_str();
      label.ReplaceAll("coorhist_", ""); // Clean up the label
      graphCorrelations->GetXaxis()->SetBinLabel(graphCorrelations->GetXaxis()->FindBin(i + 1), label);
    }
    graphCorrelations->GetXaxis()->LabelsOption("v"); // Optional: Rotate labels for readability

    summaryCanvas->Modified();
    summaryCanvas->Update();
  }



  // //Plot all restricted correlations on tgraph
  // if (!correlations.empty()) {
  //   // Create a new canvas for the summary plot
  //   TCanvas* summaryCanvas = new TCanvas("summaryCanvas", "Restricted Correlation Factors", 1200, 800);
  //   summaryCanvas->cd();

  //   const int n = correlations.size();
  //   TGraph* graph = new TGraph(n);
  //   graph->SetTitle("Restricted Correlation Factors;;Correlation Factor");

  //   for (int i = 0; i < n; ++i) {
  //     graph->SetPoint(i, i + 1, correlations[i]);
  //   }

  //   graph->GetXaxis()->Set(n, 0.5, n + 0.5);
  //   for (int i = 0; i < n; ++i) {
  //       TString label = names[i].c_str(); // Convert std::string to TString
  //       label.ReplaceAll("coorhist_", ""); // Remove the "coorhist_" part
  //       graph->GetXaxis()->SetBinLabel(graph->GetXaxis()->FindBin(i + 1), label);
  //   }

  //   graph->SetMarkerStyle(21);
  //   graph->SetMarkerColor(kRed);
  //   graph->Draw("AP");

  //   summaryCanvas->Modified();
  //   summaryCanvas->Update();

  // }


}
