//sseeds - 6.15.24: Generates one long cut table with all cut parameters.
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

//MAIN (no arguments)
void makeLongCutTable() {

  //get new cuts from .csv
  std::string cutsheet_path = "/w/halla-scshelf2102/sbs/seeds/ana/data/p2_cutset.csv";
  std::string outfile_path = "cuttable.txt";
  
  std::ifstream inFile(cutsheet_path);
  std::ofstream outFile(outfile_path);

  if (!inFile.is_open() || !outFile.is_open()) {
    std::cerr << "Error opening file!" << std::endl;
    return;
  }

  std::string line;
  std::vector<std::vector<std::string>> data;
  while (std::getline(inFile, line)) {
    std::stringstream ss(line);
    std::string cell;
    std::vector<std::string> row;
    while (std::getline(ss, cell, ',')) {
      std::string modcell = "$" + cell + "$";
      row.push_back(modcell);
    }
    data.push_back(row);
  }

  outFile << "\\begin{landscape}\n";
  outFile << "\\scriptsize\n";
  outFile << "\\begin{singlespace}\n";
  outFile << "\\begin{longtable}{|";
  for (size_t i = 0; i < data[0].size(); ++i) {
    outFile << "p{1.5cm}|";
  }
  outFile << "}\n\\caption[HCal PMT Plateaus and Alphas]{HCal PMT Plateaus and Alphas}\n";
  outFile << "\\label{tab:hcalalphas} \\\\\n\\hline\n";

  outFile << "\\endfirsthead\n\n";
  outFile << "\\multicolumn{" << data[0].size() << "}{c}%\n";
  outFile << "{{\\bfseries \\tablename\\ \\thetable{} -- continued from previous page}} \\\\\n";
  outFile << "\\hline\n";
  outFile << "\\endhead\n\n";
  outFile << "\\hline\n";
  outFile << "\\multicolumn{" << data[0].size() << "}{|r|}{{Continued on next page}} \\\\ \\hline\n";
  outFile << "\\endfoot\n\n";
  outFile << "\\hline \\hline\n";
  outFile << "\\endlastfoot\n";

  for (const auto& row : data) {
    for (const auto& cell : row) {
      outFile << cell << " &";
    }
    outFile.seekp(-1, std::ios_base::end);  // Remove the last '&'
    outFile << "\\\\ \\hline\n";
  }

  outFile << "\\end{longtable}\n";
  outFile << "\\end{singlespace}\n";
  outFile << "\\end{landscape}\n";

  inFile.close();
  outFile.close();
}
