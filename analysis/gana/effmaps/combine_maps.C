//sseeds script to combine maps over same direction
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <cstdlib>
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"

using namespace std;

void combine_maps(const char* dir, int kine, int mag_1, int mag_2) {
    stringstream file_1_path, file_2_path, output_file_path, combined_file_path, latex_file_path;
    file_1_path << dir << "corr_sbs" << kine << "_mag" << mag_1 << ".txt";
    file_2_path << dir << "corr_sbs" << kine << "_mag" << mag_2 << ".txt";
    output_file_path << dir << "corr_sbs" << kine << "_mag" << mag_1 << "_xmag" << mag_2 << ".txt";
    combined_file_path << dir << "corr_sbs" << kine << "_mag" << mag_1 << "_xmag" << mag_2 << "_combined.txt";
    latex_file_path << dir << "corr_sbs" << kine << "_mag" << mag_1 << "_xmag" << mag_2 << ".tex";

    map<double, double> map_1, map_2;

    // Load first map
    ifstream file_1(file_1_path.str().c_str());
    if (!file_1.is_open()) {
        cerr << "Error opening file: " << file_1_path.str() << endl;
        return;
    }
    string line;
    while (getline(file_1, line)) {
        if (line[0] == '#') continue;
        stringstream ss(line);
        double x, correction;
        ss >> x >> correction;
        map_1[x] = correction;
    }
    file_1.close();

    // Load second map
    ifstream file_2(file_2_path.str().c_str());
    if (!file_2.is_open()) {
        cerr << "Error opening file: " << file_2_path.str() << endl;
        return;
    }
    while (getline(file_2, line)) {
        if (line[0] == '#') continue;
        stringstream ss(line);
        double x, correction;
        ss >> x >> correction;
        map_2[x] = correction;
    }
    file_2.close();

    // Multiply maps and write to output files
    ofstream output_file(output_file_path.str().c_str());
    ofstream combined_file(combined_file_path.str().c_str());
    ofstream latex_file(latex_file_path.str().c_str());

    if (!output_file.is_open()) {
        cerr << "Error opening file: " << output_file_path.str() << endl;
        return;
    }
    if (!combined_file.is_open()) {
        cerr << "Error opening file: " << combined_file_path.str() << endl;
        return;
    }
    if (!latex_file.is_open()) {
        cerr << "Error opening file: " << latex_file_path.str() << endl;
        return;
    }

    output_file << "#HCal Efficiency correction map. List of pairs follow: <x_exp, next upper limit> <correction factor>\n";
    output_file << "#Created from SBS" << kine << " mag" << mag_1 << " and mag" << mag_2 << "\n";

    combined_file << "#HCal Efficiency correction map. Combined with original values\n";
    combined_file << "#Created from SBS" << kine << " mag" << mag_1 << " and mag" << mag_2 << "\n";

    latex_file << "\\begin{table}[h!]\n\\centering\n";
    latex_file << "\\begin{tabular}{|c|c|c|c|}\n";
    latex_file << "\\hline\n";
    latex_file << "x_exp & Correction & Original 1 & Original 2 \\\\\n";
    latex_file << "\\hline\n";

    for (const auto& pair : map_1) {
        double x = pair.first;
        double correction_1 = pair.second;
        double correction_2 = (map_2.find(x) != map_2.end()) ? map_2[x] : 1.0;
        double correction_product = correction_1 * correction_2;

        output_file << x << " " << correction_product << "\n";
        combined_file << x << " " << correction_product << " " << correction_1 << " " << correction_2 << "\n";
        latex_file << x << " & " << correction_product << " & " << correction_1 << " & " << correction_2 << " \\\\\n";
    }

    output_file.close();
    combined_file.close();

    latex_file << "\\hline\n";
    latex_file << "\\end{tabular}\n";
    latex_file << "\\caption{HCal Efficiency correction map.}\n";
    latex_file << "\\label{tab:hcal_efficiency_correction}\n";
    latex_file << "\\end{table}\n";
    latex_file.close();
}
