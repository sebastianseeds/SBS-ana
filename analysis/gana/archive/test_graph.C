#include <TGraph.h>
#include <TCanvas.h>
#include <vector>

void test_graph() {
  // Initialize vectors with 10 elements each for x and y coordinates
  std::vector<double> xpts = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
  std::vector<double> ypts = {1.0, 2.1, 1.9, 3.0, 2.5, 3.5, 3.0, 4.0, 3.9, 4.1};

  // Ensure the vectors are not empty and have the same size
  if (!xpts.empty() && xpts.size() == ypts.size()) {
    TGraph* graph = new TGraph(xpts.size(), xpts.data(), ypts.data());
    graph->SetTitle("Simple TGraph;X values;Y values");

    // Create a canvas to draw the graph
    TCanvas* c1 = new TCanvas("c1", "A Simple Graph Example", 800, 600);
    graph->SetMarkerStyle(21);  // Set marker style to filled squares
    graph->Draw("AP");  // Draw the graph with axis and points

    // Optionally, save the canvas as a file
    c1->SaveAs("graph.png");  // Save the plot as a PNG file
  } else {
    std::cerr << "Error: Vectors are empty or do not have the same size." << std::endl;
  }
}
