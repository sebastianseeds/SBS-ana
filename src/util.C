#include "../include/util.h"

namespace util {

  ////Console
  // returns string with current month day and year (date)
  string getDate(){
    time_t now = time(0);
    tm ltm = *localtime(&now);
  
    string yyyy = to_string(1900 + ltm.tm_year);
    string mm = to_string(1 + ltm.tm_mon);
    string dd = to_string(ltm.tm_mday);
    string date = mm + '_' + dd + '_' + yyyy;
  
    return date;
  }

  ////HCal
  // returns TH2D for hcal face (row,col), note that row/col start at 1
  TH2D *hhcalrowcol(std::string name) {
    TH2D *h = new TH2D(name.c_str(), ";HCAL columns;HCAL rows",
		       econst::hcalcol, 0, econst::hcalcol,
		       econst::hcalrow, 0, econst::hcalrow);
    return h;
  }

  // returns TH2D for hcal face (x,y), data
  TH2D *hhcalxy(std::string name) {
    Double_t ymin = econst::hcalposYi;
    Double_t ymax = econst::hcalposYf;
    Double_t xmin = econst::hcalposXi;
    Double_t xmax = econst::hcalposXf;

    TH2D *h = new TH2D(name.c_str(), ";hcaly_{exp} (m);hcalx_{exp} (m)",
		       econst::hcalcol, ymin, ymax,
		       econst::hcalrow, xmin, xmax);
    return h;
  }

  // returns TH2D for hcal face (x,y), mc
  TH2D *hhcalxy_mc(std::string name) {
    Double_t ymin = econst::hcalposYi_mc;
    Double_t ymax = econst::hcalposYf_mc;
    Double_t xmin = econst::hcalposXi_mc;
    Double_t xmax = econst::hcalposXf_mc;

    TH2D *h = new TH2D(name.c_str(), ";hcaly_{exp} (m);hcalx_{exp} (m)",
		       econst::hcalcol, ymin, ymax,
		       econst::hcalrow, xmin, xmax);
    return h;
  }

  // returns TH2D for hcal dxdy, wide coordinates
  TH2D *hdxdy(std::string name) {
    TH2D *h = new TH2D(name.c_str(), "; hcaly_{obs} - hcaly_{exp} (m); hcalx_{obs} - hcalx_{exp} (m)",
		       250, -1.25, 1.25, 250, -3.5, 2);
    return h;
  }

  // returns TH1D for hcal waveforms using hcal ADC bin limits
  TH1D *hhsamps(Int_t row, Int_t col, Int_t bins)
  {
    TH1D *h = new TH1D(TString::Format("h%02d%02d",row,col),
		       TString::Format("%d-%d",row+1,col+1),bins,econst::minsamp,econst::maxsamp);
    h->SetStats(0);
    h->SetLineWidth(2);
    return h;
  }

  // returns lines which draw out the edges of hcal or cut area. Pass dimensions {xmin, xmax, ymin, ymax}
  void drawarea(vector<Double_t> dimensions, Int_t lcolor=2, Int_t lwidth=4, Int_t lstyle=9) {
    Double_t top = dimensions[0];                 // -X axis
    Double_t bottom = dimensions[1];              // +X axis
    Double_t right = dimensions[2];               // -Y axis
    Double_t left = dimensions[3];                // +Y axis
    TLine line;
    line.SetLineColor(lcolor); 
    line.SetLineWidth(lwidth); 
    line.SetLineStyle(lstyle);
    line.DrawLine(right, bottom, left, bottom); // bottom margin
    line.DrawLine(right, top, left, top);       // top margin
    line.DrawLine(right, top, right, bottom);   // right margin
    line.DrawLine(left, top, left, bottom);     // left margin
  }


  ////kinematic histograms
  // returns W2 histogram
  TH1D *hW2(std::string name) {
    TH1D *h = new TH1D(name.c_str(), "W^{2} Distribution (GeV^{2})", 250,0,2);
    return h;
  }
  // returns Q2 histogram
  TH1D *hQ2(std::string name,       // Name of histogram
		Int_t conf) {             // SBS config
    Int_t nbin=0; Double_t hmin=-100, hmax=-100;
    if (conf==4) { nbin=100; hmin=1.; hmax=4.; } 
    else if (conf==14) { nbin=100; hmin=5.; hmax=10.; }
    else if (conf==7) { nbin=120; hmin=6.; hmax=12.; }
    else if (conf==11) { nbin=200; hmin=8.; hmax=18.; }
    else if (conf==8 || conf==9) { nbin=120; hmin=3.; hmax=9.; }
    else cerr << "Error on [util::hQ2], enter valid kinematic." << endl;
    TH1D *h = new TH1D(name.c_str(), "Q^{2} Distribution (GeV^{2})", 
		       nbin, hmin, hmax);
    return h;
  }

  // functions to read csv files
  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   Int_t &nruns,                // No. of runs to analyze
		   Int_t sbsconf,               // SBS configuration
		   std::string target,        // target type
		   Int_t replay_pass,           // replay pass
		   Int_t verbose,               // verbosity
		   vector<crun> &corun)     // Output: Vector of crun structs
  {
    // Define the name of the relevant run spreadsheet
    std::string fst = "/grl_"; 
    std::string mid = "_pass";
    std::string lst = ".csv";
    if (replay_pass < 2) replay_pass = 1; // single spreadsheet exists for pass 0 & 1
    std::string run_spreadsheet = runsheet_dir + fst + target + mid + std::to_string(replay_pass) + lst;

    // Reading the spreadsheet
    if (nruns < 0) nruns = 1e6;         // replay all runs if nruns < 0
    ifstream run_data; run_data.open(run_spreadsheet);
    string readline;
    if(run_data.is_open()){
      std::cout << "Reading run info from: "<< run_spreadsheet 
		<< std::endl << std::endl;
      string skip_header; getline(run_data, skip_header);  // skipping column header
      corun.clear();
      while(getline(run_data,readline)){                   // reading each line
	istringstream tokenStream(readline);
	string token;
	char delimiter = ',';
	vector<string> temp;
	while(getline(tokenStream,token,delimiter)){       // reading each element of a line
	  string temptoken=token;
	  temp.push_back(temptoken);
	}
	// add relevant info to crun objects
	if (stoi(temp[0]) == sbsconf) {
	  if (corun.size() >= nruns) break;
	  crun temp_cr;
	  temp_cr.SetDataRunSheet(temp);
	  corun.push_back(temp_cr);

	}

	temp.clear();
      }
      // let's update nruns with total no. of runs to analyze
      nruns = corun.size();
      if (verbose == 1) {
	std::cout << "First run info:" << std::endl << corun[0];
	std::cout << "Last run info:" << std::endl << corun[nruns-1];
      }
    }else
      throw "Error on [util::ReadRunList] Run spreadsheet doesn't exist";
    run_data.close();
  }
  //overload 1
  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   Int_t &nruns,                // No. of runs to analyze
		   Int_t sbsconf,               // SBS configuration
		   std::string target,        // target type
		   Int_t replay_pass,           // replay pass
		   Int_t sbsmag,                // SBS magnet current (in %)
		   Int_t verbose,               // verbosity
		   vector<crun> &corun)     // Output: Vector of crun structs
  {
    // Define the name of the relevant run spreadsheet
    std::string fst = "/grl_"; 
    std::string mid = "_pass";
    std::string lst = ".csv";
    if (replay_pass < 2) replay_pass = 1; // single spreadsheet exists for pass 0 & 1
    std::string run_spreadsheet = runsheet_dir + fst + target + mid + std::to_string(replay_pass) + lst;

    // convert magnet field values from % to A
    sbsmag *= 21;    // 100% SBS magnet current = 2100 A

    // Reading the spreadsheet
    if (nruns < 0) nruns = 1e6;         // replay all runs if nruns < 0
    ifstream run_data; run_data.open(run_spreadsheet);
    string readline;
    if(run_data.is_open()){
      std::cout << "Reading run info from: "<< run_spreadsheet 
		<< std::endl << std::endl;
      string skip_header; getline(run_data, skip_header); // skipping column header
      while(getline(run_data,readline)){                  // reading each line
	istringstream tokenStream(readline);
	string token;
	char delimiter = ',';
	vector<string> temp;
	while(getline(tokenStream,token,delimiter)){      // reading each element of a line
	  string temptoken=token;
	  temp.push_back(temptoken);
	}
	// add relevant info to crun objects
	if (stoi(temp[0]) == sbsconf && 
	    stoi(temp[3]) == sbsmag) {
	  if (corun.size() >= nruns) break;
	  crun temp_cr;
	  temp_cr.SetDataRunSheet(temp);
	  corun.push_back(temp_cr);
	}

	temp.clear();
      }
      // let's update nruns with total no. of runs to analyze
      nruns = corun.size();
      if (verbose == 1) {
	std::cout << "First run info:" << std::endl << corun[0];
	std::cout << "Last run info:" << std::endl << corun[nruns-1];
      }
    }else
      throw "Error on [util::ReadRunList], run spreadsheet doesn't exist.";
    run_data.close();
  }
  //overload 2
  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   Int_t &nruns,                // No. of runs to analyze
		   Int_t sbsconf,               // SBS configuration
		   std::string target,        // target type
		   Int_t replay_pass,           // replay pass
		   Int_t sbsmag,                // SBS magnet current (in %)
		   Int_t bbmag,                 // BB magnet current (in %)
		   Int_t verbose,               // verbosity
		   vector<crun> &corun)     // Output: Vector of crun structs
  {
    // Define the name of the relevant run spreadsheet
    std::string fst = "/grl_"; 
    std::string mid = "_pass";
    std::string lst = ".csv";
    if (replay_pass < 2) replay_pass = 1; // single spreadsheet exists for pass 0 & 1
    std::string run_spreadsheet = runsheet_dir + fst + target + mid + std::to_string(replay_pass) + lst;

    // convert magnet field values from % to A
    sbsmag *= 21;    // 100% SBS magnet current = 2100 A
    bbmag *= 7.5;    // 100% BB magnet current = 750 A

    // Reading the spreadsheet
    if (nruns < 0) nruns = 1e6;         // replay all runs if nruns < 0
    ifstream run_data; run_data.open(run_spreadsheet);
    string readline;
    if(run_data.is_open()){
      std::cout << "Reading run info from: "<< run_spreadsheet 
		<< std::endl << std::endl;     
      string skip_header; getline(run_data, skip_header); // skipping column header
      while(getline(run_data,readline)){                  // reading each line
	istringstream tokenStream(readline);
	string token;
	char delimiter = ',';
	vector<string> temp;
	while(getline(tokenStream,token,delimiter)){      // reading each element of a line
	  string temptoken=token;
	  temp.push_back(temptoken);
	}
	// add relevant info to crun objects
	if (stoi(temp[0]) == sbsconf 
	    && stoi(temp[3]) == sbsmag 
	    && stoi(temp[4]) == bbmag) {
	  if (corun.size() >= nruns) break;
	  crun temp_cr;
	  temp_cr.SetDataRunSheet(temp);
	  corun.push_back(temp_cr);
	}

	temp.clear();
      }
      // let's update nruns with total no. of runs to analyze
      nruns = corun.size();
      if (verbose == 1) {
	std::cout << "First run info:" << std::endl << corun[0];
	std::cout << "Last run info:" << std::endl << corun[nruns-1];
      }
    }else
      throw "Error on [util::ReadRunList], run spreadsheet doesn't exist";
    run_data.close();
  }

  // functions to read pd sort files
  void ListDirectory(const char *path,std::vector<std::string> &list){
    struct dirent *entry;
    DIR *dir = opendir(path);
    if(dir==NULL){
      return;
    }
    std::string myStr;
    while ((entry = readdir(dir)) != NULL) {
      // std::cout << " mystr= " << myStr << std::endl;
      myStr = entry->d_name;
      list.push_back(myStr);
    }
    closedir(dir);
  }

  Int_t SplitString(const char delim, const std::string myStr, std::vector<std::string> &out){
    // split a string by a delimiter
    std::stringstream ss(myStr);
    std::vector<std::string> result;
    while( ss.good() ){
      std::string substr;
      std::getline(ss, substr, delim);
      out.push_back(substr);
    }
    return 0;
  }

  bool sortbyval(const pair<Int_t, Int_t> &a, const pair<Int_t, Int_t> &b) {
    return (a.first < b.first);
  } 

  Int_t GetROOTFileMetaData(const char *rfDirPath, Int_t run,
			  std::vector<Int_t> &data, 
			  std::vector<pair<Int_t, Int_t>> &segB_segE,
			  Int_t verbose){
    // determine the beginning and end segment number for a CODA run
    // input: 
    // - rfDirPath: path to where the ROOT file(s) are located 
    // - run: CODA run number 
    // output: vector containing:  
    // - stream: stream number of the EVIO file associated with the run  
    // - begSeg: start segment number of the ROOT file associated with run 
    // - endSeg: ending segment number of the ROOT file associated with run
    // - number of files associated with the run

    Int_t rc=-1; // assume fail 

    // first get list of files in the directory 
    std::vector<std::string> fileList;
    ListDirectory(rfDirPath,fileList); 

    if (verbose > 1) {
      std::cout << "----" << std::endl;
      std::cout << " fileList.size() " << fileList.size() << std::endl;
      for (Int_t i = 0; i < fileList.size(); i++) {
	std::cout << " fileList[" << i << "] = " << fileList[i] << std::endl;
      }
      std::cout << "----" << std::endl;
    }

    // skip first two entries since we get . and .. at top of list
    for(Int_t i=0;i<2;i++) fileList.erase(fileList.begin()); 

    // identify the right index, and find number of files 
    Int_t j=-1, fileCnt=0, theRun=0, stream=0;
    const Int_t NF = fileList.size();
    std::vector<std::string> o1;
    for(Int_t i=0; i<NF; i++){ 

      // lets implement a filter 
      if (fileList[i].find("e1209019_fullreplay") == 0) {

	// split each entry based on a sufficient delimiter 
	SplitString('_', fileList[i], o1);
	// 3rd entry (index 2) is the one that has the run number 
	theRun = std::atoi(o1[2].c_str());
	if(theRun==run){
	  fileCnt++; 

	  // split each entry based on a sufficient delimiter 
	  std::vector<std::string> o2, o3; 
	  // determine the stream (index 3) 
	  std::string theStr = o1[3]; 
	  SplitString('m', theStr, o2); 
	  stream = std::atoi(o2[1].c_str());
	  o2.clear();
	  // determine the segment numbers (index 4 and 5) 
	  theStr = o1[4]; 
	  SplitString('g', theStr, o2); 
	  Int_t bseg = std::atoi(o2[1].c_str());
	  Int_t eseg = std::atoi(o1[5].c_str());
	  segB_segE.push_back( make_pair(bseg, eseg) );
	  o2.clear();
	}
	o1.clear();
      } 
    }

    // lets sort the segments by beginning segment number 
    // this is necessary for proper beam charge calculation
    sort(segB_segE.begin(), segB_segE.end(), sortbyval);  

    if(fileCnt==0){
      throw "Error on [util::GetROOTFileMetaData], root file directory is empty.";
    }

    data.push_back(stream); 
    data.push_back(fileCnt);
 
    return 0;
  } 

  Int_t LoadROOTTree(std::string path,
		   std::vector<crun> corun,    // multiple CODA runs
		   bool sort,
		   Int_t verbose,
		   TChain* &C) 
  {
    
    if (!corun.empty()) {
      const Int_t nruns = corun.size();
      std::cout << "Parsing ROOT files from " << nruns << " runs.." << std::endl;
      if (!sort) {
	for (Int_t i=0; i<nruns; i++) {
	  std::string rfname = Form("%s/*%d*",path.c_str(),corun[i].runnum);
	  if (verbose > 1) std::cout << rfname << std::endl;
	  C->Add(rfname.c_str());
	}
      } else {
	std::cout << "Sorting by segments.." << std::endl;
	Int_t aRun, aNumFiles, aStream, rc;
	std::vector<Int_t> md;
	std::vector<pair<Int_t, Int_t>> segB_segE;
	// Looping through unique runs
	for (Int_t irun=0; irun<nruns; irun++) {
	  aRun = corun[irun].runnum;
	  rc = GetROOTFileMetaData(path.c_str(),aRun,md,segB_segE,0);
	  aStream   = md[0]; //stream no.
	  aNumFiles = md[1]; //total # segments
	  if (verbose > 0) {
	    std::cout << "----" << std::endl;
	    std::cout << Form(" Run %d, Total # of segments %d",aRun,aNumFiles) << std::endl;
	    std::cout << Form(" Beg seg %d-%d, End seg %d-%d",segB_segE[0].first,segB_segE[0].second,
			      segB_segE[aNumFiles-1].first,segB_segE[aNumFiles-1].second) << std::endl;
	    std::cout << "----" << std::endl;
	  }
	  // Looping through segments and adding to tree
	  for (Int_t iseg=0; iseg<aNumFiles; iseg++) {
	    std::string rfname = Form("%s/e1209019_fullreplay_%d_stream%d_seg%d_%d.root",path.c_str(),
				      aRun,aStream,segB_segE[iseg].first,segB_segE[iseg].second);
	    if (verbose > 1) std::cout << rfname << std::endl;
	    C->Add(rfname.c_str());
	  }
	  // getting ready for next run
	  md.clear();
	  segB_segE.clear();
	}
      }
      if (C->GetEntries()==0) 
	throw "Error on [util::LoadROOTTree], empty/nonexistant root tree.";
    }else 
      throw "Error on [util::LoadROOTTree], empty coda run list.";

    return 0;
  }

  Int_t LoadROOTTree(std::string path,
		   crun corun,           // single CODA run
		   bool sort,
		   Int_t verbose,
		   TChain* &C) 
  {
    
    if (corun.runnum != 0) {
       std::cout << "Parsing ROOT files from run " << corun.runnum << std::endl;
      if (!sort) {
	std::string rfname = Form("%s/*%d*",path.c_str(),corun.runnum);
	if (verbose > 1) std::cout << rfname << std::endl;
	C->Add(rfname.c_str());
      } else {
	std::cout << "Sorting by segments.." << std::endl;
	Int_t aRun, aNumFiles, aStream, rc;
	std::vector<Int_t> md;
	std::vector<pair<Int_t, Int_t>> segB_segE;
	// Looping through unique runs
	aRun = corun.runnum;
	rc = GetROOTFileMetaData(path.c_str(),aRun,md,segB_segE,0);
	aStream   = md[0]; //stream no.
	aNumFiles = md[1]; //total # segments
	if (verbose > 0) {
	  std::cout << "----" << std::endl;
	  std::cout << Form(" Run %d, Total # of segments %d",aRun,aNumFiles) << std::endl;
	  std::cout << Form(" Beg seg %d-%d, End seg %d-%d",segB_segE[0].first,segB_segE[0].second,
			    segB_segE[aNumFiles-1].first,segB_segE[aNumFiles-1].second) << std::endl;
	  std::cout << "----" << std::endl;
	}
	// Looping through segments and adding to tree
	for (Int_t iseg=0; iseg<aNumFiles; iseg++) {
	  std::string rfname = Form("%s/e1209019_fullreplay_%d_stream%d_seg%d_%d.root",path.c_str(),
				    aRun,aStream,segB_segE[iseg].first,segB_segE[iseg].second);
	  if (verbose > 1) std::cout << rfname << std::endl;
	  C->Add(rfname.c_str());
	}
	// getting ready for next run
	md.clear();
	segB_segE.clear();
      }
      if (C->GetEntries()==0) 
	throw "Error on [util::LoadROOTTree], empty/missing root file.";
    }else 
      throw "Error on [util::LoadROOTTree], crun object is empty.";
    
    return 0;
  }
}
