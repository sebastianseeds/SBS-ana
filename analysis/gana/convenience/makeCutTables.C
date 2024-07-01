//sseeds
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>

std::string sanitizeLaTeX(const std::string& input) {
  std::string sanitized;
  for (char c : input) {
    if (c == '_') {
      sanitized += "\\_";
    } else {
      sanitized += c;
    }
  }
  return sanitized;
}

std::string verbizeLaTeX(const std::string& input) {
  
  std::string verbized = "\\verb|" + input + "|";
  
  return verbized;
}

void makeCutTables() {

  //get new cuts from .csv
  std::string cutsheet_path = "/w/halla-scshelf2102/sbs/seeds/ana/data/p2_cutset.csv";
  
  std::ifstream inFile(cutsheet_path);
  std::ofstream outFile("manycutsheets.txt");

  if (!inFile.is_open() || !outFile.is_open()) {
    std::cerr << "Error opening file!" << std::endl;
    return;
  }

  std::string line;
  std::vector<std::vector<std::string>> data;
  std::vector<std::vector<std::string>> data2;
  while (std::getline(inFile, line)) {
    std::stringstream ss(line);
    std::string cell;
    std::vector<std::string> row;
    std::vector<std::string> row2;
    while (std::getline(ss, cell, ',')) {
      std::string sancell = sanitizeLaTeX(cell);
      std::string verbcell = verbizeLaTeX(cell);
      
      row.push_back(verbcell);
      row2.push_back(sancell);
    }
    data.push_back(row);
    data2.push_back(row2);
  }

  std::vector<std::string> headers = data[0];
  std::vector<std::string> headers2 = data2[0];
  data.erase(data.begin()); // Remove headers from data

  for (size_t col = 3; col < headers.size(); ++col) {
    outFile << "\\begin{table}[ht]\n";
    outFile << "\\centering\n";
    outFile << "{\\footnotesize\n";
    outFile << "\\begin{tabular}{|c|c|c|c|}\n";
    outFile << "\\hline\n";
    outFile << headers[0] << " & " << headers[1] << " & " << headers[2] << " & " << headers[col] << " \\\\ \\hline\n";

    for (const auto& row : data) {
      outFile << row[0] << " & " << row[1] << " & " << row[2] << " & " << row[col] << " \\\\ \\hline\n";
    }

    outFile << "\\end{tabular}\n";
    outFile << "}\n";
    outFile << "\\captionsetup{width=0.7\\textwidth}\n";
    outFile << "\\caption[Elastic Selection Cuts: " << headers2[col] << "]{\\small Elastic Selection Cuts: " << headers2[col] << "}\n";
    outFile << "\\label{tab:" << headers2[col] << "}\n";
    outFile << "\\end{table}\n\n";

    outFile << "NEW TABLE \n\n";
    
  }

  inFile.close();
  outFile.close();
}
