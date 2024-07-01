//sseeds - 7.24.25: script to sort output SIMC MC .csv files in terms of job number
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <sstream>

struct CsvRow {
  int jobid;
  std::vector<std::string> data;
};

bool compareJobid(const CsvRow& a, const CsvRow& b) {
  return a.jobid < b.jobid;
}

//MAIN. nucleon=nucleon files to sort (0=neutron, 1=proton)
int sortcsv(int nucleon = 1) {

  std::string nucword = "n";
  if(nucleon==1)
    nucword = "p";

  std::string inputFile = Form("p390sf_sbs4_sbs30p_simc_dee%s_summary.csv",nucword.c_str());
  std::string outputFile = Form("sorted_p390sf_sbs4_sbs30p_simc_dee%s_summary.csv",nucword.c_str());

  std::ifstream infile(inputFile);
  std::ofstream outfile(outputFile);

  if (!infile.is_open() || !outfile.is_open()) {
    std::cerr << "Error opening file!" << std::endl;
    return 1;
  }

  std::string line;
  std::vector<CsvRow> rows;

  // Read the header
  std::getline(infile, line);
  outfile << line << std::endl;

  // Read the data rows
  while (std::getline(infile, line)) {
    std::istringstream ss(line);
    std::string item;
    CsvRow row;
        
    std::getline(ss, item, ',');
    row.jobid = std::stoi(item);
    row.data.push_back(item);
        
    while (std::getline(ss, item, ',')) {
      row.data.push_back(item);
    }
        
    rows.push_back(row);
  }

  // Sort the rows by jobid
  std::sort(rows.begin(), rows.end(), compareJobid);

  // Write the sorted rows to the output file
  for (const auto& row : rows) {
    for (size_t i = 0; i < row.data.size(); ++i) {
      outfile << row.data[i];
      if (i < row.data.size() - 1) {
	outfile << ",";
      }
    }
    outfile << std::endl;
  }

  infile.close();
  outfile.close();

  return 0;
}
