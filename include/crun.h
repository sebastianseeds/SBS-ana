#ifndef CRUN_H
#define CRUN_H

typedef struct crun {

  Int_t runnum;
  Int_t sbsconf;
  std::string target;
  Int_t sbsmag;           // SBS magnet current (A)
  Int_t bbmag;            // BB magnet current (A)
  Double_t ebeam;         // GeV, avg. over entire run 
  Double_t ebeam_std;     // GeV, std. dev. over entire run 
  Double_t charge;        // C, total charge collected by the run
  Double_t DAQltime;      // %

  // constructor 
  crun(): 
  runnum(0),sbsconf(0),target("NONE"),sbsmag(0),bbmag(0),ebeam(0),ebeam_std(0),charge(0),DAQltime(0)
  {}

  // define an ostream operator to print to screen conveniently
  friend ostream& operator <<(ostream &out, const crun& corun) {
    out << " -----------------------------------------" << std::endl;
    out << " Run number             : " << corun.runnum << std::endl;
    out << " SBS config             : " << corun.sbsconf << std::endl;
    out << " Target                 : " << corun.target << std::endl;
    out << " SBS mag. cur. (A)      : " << corun.sbsmag << std::endl;
    out << " BB mag. cur. (A)       : " << corun.bbmag << std::endl;
    out << " Avg. ebeam (GeV)       : " << corun.ebeam << std::endl;
    out << " Ebeam std. dev. (GeV)  : " << corun.ebeam_std << std::endl;
    out << " Charge (C)             : " << corun.charge << std::endl;
    out << " -----------------------------------------" << std::endl << std::endl;
    return out;
  }

  // sets data by reading runsheet (exclusively for util::ReadRunList functions)
  void SetDataRunSheet(std::vector<std::string> data) {
    sbsconf = stoi(data[0]);
    runnum = stoi(data[1]);
    target = data[2];
    sbsmag = stoi(data[3]);
    bbmag = stoi(data[4]);
    ebeam = stod(data[5]);
    ebeam_std = stod(data[6]);
    charge = stod(data[7]);
  }

} crun_t;  

#endif
