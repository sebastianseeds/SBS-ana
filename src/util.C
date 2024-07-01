#include "../include/util.h"
#include "../include/fits.h"

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

  //sends current progress to console
  void updateProgress(Int_t currentIteration, Int_t totalIterations) {
    const Int_t barWidth = 50;

    // Calculate the percentage completion
    double progress = static_cast<double>(currentIteration) / totalIterations;

    // Calculate the number of characters to fill in the progress bar
    Int_t barLength = static_cast<Int_t>(barWidth * progress);

    // Calculate estimated time to completion
    static auto startTime = std::chrono::high_resolution_clock::now();
    auto currentTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedTime = currentTime - startTime;
    double estimatedTimeRemaining = (elapsedTime.count() / progress) - elapsedTime.count();

    // Print the progress bar
    std::cout << "[" << std::setw(barWidth) << std::left << std::string(barLength, '=') << "] ";
    
    // Print the percentage complete
    std::cout << std::fixed << std::setprecision(1) << progress * 100.0 << "% ";

    // Print the estimated time to completion
    Int_t minutes = static_cast<Int_t>(estimatedTimeRemaining / 60);
    Int_t seconds = static_cast<Int_t>(estimatedTimeRemaining) % 60;
    std::cout << "ETC: " << minutes << "m " << seconds << "s\r";
    
    // Flush the output to update the terminal
    std::cout.flush();
  }

  //sends current progress to console with a label (overload)
  void updateProgress(Int_t currentIteration, Int_t totalIterations, std::string label) {
    const Int_t barWidth = 50;

    // Calculate the percentage completion
    double progress = static_cast<double>(currentIteration) / totalIterations;

    // Calculate the number of characters to fill in the progress bar
    Int_t barLength = static_cast<Int_t>(barWidth * progress);

    // Calculate estimated time to completion
    static auto startTime = std::chrono::high_resolution_clock::now();
    auto currentTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedTime = currentTime - startTime;
    double estimatedTimeRemaining = (elapsedTime.count() / progress) - elapsedTime.count();

    // Print the progress bar
    std::cout << "[" << std::setw(barWidth) << std::left << std::string(barLength, '=') << "] ";
    
    // Print the percentage complete
    std::cout << std::fixed << std::setprecision(1) << progress * 100.0 << "% ";

    // Print the estimated time to completion
    Int_t minutes = static_cast<Int_t>(estimatedTimeRemaining / 60);
    Int_t seconds = static_cast<Int_t>(estimatedTimeRemaining) % 60;
    std::cout << label << ". ETC: " << minutes << "m " << seconds << "s\r";
    
    // Flush the output to update the terminal
    std::cout.flush();
  }

  void updateProgress(int currentIteration_outer, int currentIteration_inner, int totalIterations_outer, int totalIterations_inner,std::string label) {
    const int barWidth = 50;


    // Define whirlwind characters
    const std::string whirlwindChars = "/-\\|";

    // Calculate the index for whirlwind animation
    int whirlwindIndex = (currentIteration_outer + currentIteration_inner) % whirlwindChars.length();

    // Calculate the percentage completion for the outer loop
    double progressOuter = static_cast<double>(currentIteration_outer) / totalIterations_outer;

    // Calculate the percentage completion for the inner loop
    double progressInner = static_cast<double>(currentIteration_inner) / totalIterations_inner;

    // Calculate the combined progress based on the nested loops
    double combinedProgress = progressOuter + (progressInner / totalIterations_outer);

    // Calculate the number of characters to fill in the progress bar
    int barLength = static_cast<int>(barWidth * combinedProgress);

    // Calculate estimated time to completion
    static auto startTime = std::chrono::high_resolution_clock::now();
    auto currentTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedTime = currentTime - startTime;
    double estimatedTimeRemaining = (elapsedTime.count() / combinedProgress) - elapsedTime.count();

    // Print the progress bar
    std::cout << "[" << std::setw(barWidth) << std::left << std::string(barLength, '=') << whirlwindChars[whirlwindIndex] << "] ";

    // Print the percentage complete
    std::cout << std::fixed << std::setprecision(1) << combinedProgress * 100.0 << "% ";

    // Print the estimated time to completion
    int minutes = static_cast<int>(estimatedTimeRemaining / 60);
    int seconds = static_cast<int>(estimatedTimeRemaining) % 60;
    std::cout << label << ". ETC: " << minutes << "m " << seconds << "s\r";

    // Flush the output to update the terminal
    std::cout.flush();
  }

  TString FindFileInDirectories(const TString& pattern, const std::vector<TString>& directories) {
    TRegexp re(pattern, kTRUE);

    // Loop over the list of directories
    for (const auto& dir : directories) {
      void* dirp = gSystem->OpenDirectory(dir);
      const char* entry;
        
      // Iterate over all files in the directory
      while ((entry = gSystem->GetDirEntry(dirp))) {
	TString file = entry;
	// Check if the current file name matches the pattern
	if (file.Contains(re)) {
	  gSystem->FreeDirectory(dirp);
	  return dir; // Return the directory where the file is found
	}
      }
      gSystem->FreeDirectory(dirp);
    }
    // File not found in any of the directories; return empty string
    return "";
  }

  //parse global cut string delimited by && into individual cut strings
  std::vector<std::string> parseCuts(const std::string& cutsString) {
    std::vector<std::string> cuts;
    std::string delimiter = "&&";
    size_t pos = 0;
    std::string token;
    std::string s = cutsString + delimiter; // Ensure the string ends with delimiter to catch the last condition

    while ((pos = s.find(delimiter)) != std::string::npos) {
      token = s.substr(0, pos);
      cuts.push_back(token);
      s.erase(0, pos + delimiter.length());
    }

    return cuts;
  }

  // reads a set of parameters into a vector from a txt file at <const_path>. Assumes endl delimiter
  void readParam( std::string const_path, vector<Double_t> &param ) {

    ifstream const_file; const_file.open( const_path );

    //Get all params
    if(const_file.is_open()){
  
      Int_t n1=0;
      Double_t d1;
      string readline;
    
      while( getline( const_file, readline ) ){
	if( readline.at(0) == '#' ) {
	  continue;
	}
	stringstream ss( readline );
	ss >> d1;     
	param.push_back(d1);
	n1++;
      }
      
    }else
      throw "Error on [util::readConst] Parameter file doesn't exist";

    const_file.close();
  }

  // overload 1. Allows for selection of subset of data with type and expected set size
  void readParam( std::string const_path, std::string type, Int_t setsize, vector<Double_t> &param ) {
    ifstream const_file; const_file.open( const_path );
    
    //Get all params at datatype
    if(const_file.is_open()){

      Int_t n1=0;
      Double_t d1=0;
      string readline;
      bool skip_line = true;

      while( getline( const_file, readline ) ){
      
	if( n1==( setsize ) ) break;
      
	TString Tline = (TString)readline;
      
	if( Tline.BeginsWith( type.c_str() ) && skip_line==true ){	
	  skip_line = false;
	  continue;
	}
      
	if( skip_line==false ){
	  istringstream ss( readline );
	  while( ss >> d1 ){
	    param.push_back(d1);
	    n1++;
	  }
	}
      }
    }else
      throw "Error on [util::readConst] Parameter file doesn't exist";

    const_file.close();
  }

  //Generates random number with passkey intended for blinding yields/scale-factors
  Double_t generateRandomNumber(const std::string& passkey) {
    if (passkey.length() != 4) {
      throw std::invalid_argument("Passkey must be a four-letter string");
    }

    // Seed the random number generator
    unsigned int seed = 0;
    for (char ch : passkey) {
      seed += static_cast<unsigned int>(ch);
    }
    srand(seed);

    // Generate a random number between 0.95 and 1.05
    Double_t randomFraction = static_cast<Double_t>(rand()) / RAND_MAX; // Range [0, 1)
    return 0.95 + randomFraction * 0.1; // Scale to range [0.95, 1.05)
  }

  std::vector<std::string> parseGlobalCut(const std::string& input) {
    std::vector<std::string> elements;
    std::string token;
    std::istringstream tokenStream(input);

    while (std::getline(tokenStream, token, '&')) {
      if (!token.empty() && token.back() != '>' && token.back() != '<' && token.back() != '=') {
	// Remove any leading '&' left from splitting
	if (token.front() == '&') {
	  token.erase(0, 1);
	}
	elements.push_back(token);
      } else {
	// Handle case where '&&' is part of a comparison
	std::string nextToken;
	std::getline(tokenStream, nextToken, '&');
	token += "&" + nextToken;
	elements.push_back(token);
      }
    }

    return elements;
  }

  //Function to see if file exists at path
  bool checkFile(const std::string& filePath) {
    std::ifstream file(filePath.c_str());
    return file.good();
  }

  //Function to check TH1D histogram loaded from dataloop script output and return those that don't exist
  bool checkTH1D(TH1D* histogram, const std::string& name) {
    if (!histogram) {
        std::cerr << "Failed to load histogram: " << name << std::endl;
        return false;
    }
    return true;
  }
  
  ////HCal

  // Check position per event against hcal active area and return minimum N_x(N_y)*blockwidth_x(y) to include cluster center. X= .first, y= .second
  std::pair<double, double> minaa(double hcalx, double hcaly) {
    // HCAL dimensions from some constants namespace or directly defined
    double hcalaaXi = econst::hcalposXi_mc; // Begin of acceptance in x
    double hcalaaXf = econst::hcalposXf_mc; // End of acceptance in x
    double hcalaaYi = econst::hcalposYi_mc; // Begin of acceptance in y
    double hcalaaYf = econst::hcalposYf_mc; // End of acceptance in y

    double block_wx = econst::hcalblk_div_v; // Block width in x
    double block_wy = econst::hcalblk_div_h; // Block width in y

    // Calculate the minimum difference between the cluster center and the edges of the acceptance area in units of block width. Ensure the cluster center is within the acceptance area.

    // For x:
    double minDistX = std::min(std::abs(hcalx - hcalaaXi), std::abs(hcalx - hcalaaXf)) / block_wx;
    // For y:
    double minDistY = std::min(std::abs(hcaly - hcalaaYi), std::abs(hcaly - hcalaaYf)) / block_wy;

    // Create a pair to return both minimum distances (scales)
    std::pair<double, double> block_scales = std::make_pair(minDistX, minDistY);

    return block_scales;
  }

  // checks if a point is within proton or neutron spot. rotation angle in rad, if ever applicable
  bool Nspotcheck(Double_t dy, Double_t dx, Double_t dy_mean, Double_t dx_mean, Double_t dy_sigma, Double_t dx_sigma, Double_t rotationAngle=0) {

    // Caculate semimajor and semiminor axes
    Double_t dyAxis = dy_sigma;
    Double_t dxAxis = dx_sigma;

    // Apply rotation angle if applicable
    Double_t cosAngle = std::cos(rotationAngle);
    Double_t sinAngle = std::sin(rotationAngle);
    Double_t dyRot = (dy - dy_mean) * cosAngle + (dx - dx_mean) * sinAngle;
    Double_t dxRot = (dx - dx_mean) * cosAngle - (dy - dy_mean) * sinAngle;

    // Check if point is within the ellipse equation
    Double_t result = ((dyRot * dyRot) / (dyAxis * dyAxis)) + ((dxRot * dxRot) / (dxAxis * dxAxis));
    return result <= 1.0;
  }

  // Function that calculates and returns the scale factor N
  double NspotScaleFactor(Double_t dy, Double_t dx, Double_t dy_mean, Double_t dx_mean, Double_t dy_sigma, Double_t dx_sigma, Double_t rotationAngle=0) {
    // Calculate semimajor and semiminor axes
    Double_t dyAxis = dy_sigma;
    Double_t dxAxis = dx_sigma;

    // Apply rotation angle if applicable
    Double_t cosAngle = std::cos(rotationAngle);
    Double_t sinAngle = std::sin(rotationAngle);
    Double_t dyRot = (dy - dy_mean) * cosAngle + (dx - dx_mean) * sinAngle;
    Double_t dxRot = (dx - dx_mean) * cosAngle - (dy - dy_mean) * sinAngle;

    // Calculate the scale factor N such that the point is on the boundary of the ellipse
    // For a point (x, y) to lie on the ellipse defined by axes a*N and b*N, where a and b are the original semiaxes, we must have (x^2)/(a*N)^2 + (y^2)/(b*N)^2 = 1. Rearranging for N gives us the following:
    Double_t scaleFactor_dy = std::sqrt((dyRot * dyRot) / (dyAxis * dyAxis));
    Double_t scaleFactor_dx = std::sqrt((dxRot * dxRot) / (dxAxis * dxAxis));

    // The actual scale factor N needed is the maximum of these two, to ensure the point lies within or on the boundary
    Double_t N = std::max(scaleFactor_dy, scaleFactor_dx);

    return N;
  }

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

  //assign score to potential cluster based on probability density of gaussian fit to atime and maximum energy, exclude position info. Here, gfit is a vector with coin atime fit parameters, val_2 is the cluster coin atime, val_1 is the cluster energy, and max1 is the max cluster energy for the event
  Double_t assignScore( double val_1, double val_2, double max1, const std::vector<double> &gfit ) {

    if( gfit.size()!=3 ){
      cout << "ERROR: size of dxfit vector not equal to expected number of gaussian parameters (3)." << endl;
      return 0.;
    }

    // Build a TF1 with fit provided fit parameters
    TF1 *gauss = new TF1("gauss", "gaus", 0, 100);
    gauss->SetParameter(0,gfit[0]);
    gauss->SetParameter(1,gfit[1]);
    gauss->SetParameter(2,gfit[2]);
    
    Double_t x_value = val_2;
    Double_t density = gauss->Eval(x_value);

    // Compute the score based on val 2
    Double_t score = density / gauss->GetParameter(0);  // Normalize by the peak of the Gaussian

    // Include a product component based on val 1 (linear weight)
    score *= val_1 / max1;

    if(score==0){
      //cout << "WARNING: score is zero. val_1=" << val_1 << ", val_2= " << val_2 << " density= " << density << " 1 comp= " << val_1 / max1 << endl;
      score=1e-38; //write very small number to score to exclude it without breaking the sorting later
    }

    delete gauss; //cleanup

    return score;
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

  //utility function to search two column csv (typical simc/g4sbs output summary file) for value
  Double_t searchSimpleCSVForValue(const std::string& filePath, const std::string& key) {
    std::ifstream file(filePath);
    std::string line;

    if (!file.is_open()) {
      std::cerr << "Error opening file: " << filePath << std::endl;
      return -1.0; // Error code, assuming all valid values are positive
    }

    while (std::getline(file, line)) {
      std::istringstream lineStream(line);
      std::string currentKey;
      std::string valueStr;

      // Get the key and the value as strings
      if (std::getline(lineStream, currentKey, ',') && std::getline(lineStream, valueStr)) {
	// Check if the current key matches the desired key
	if (currentKey == key) {
	  // Convert the value to double and return it
	  try {
	    double value = std::stod(valueStr);
	    return value;
	  } catch (const std::invalid_argument& ia) {
	    std::cerr << "Invalid argument: " << ia.what() << std::endl;
	    return -1.0;
	  }
	}
      }
    }

    // Key not found
    std::cerr << "Key not found: " << key << std::endl;
    return -1.0;
  }

  //utility function to parse simc hist file for JBoyd simc replayed .root files Fall 2023
  Double_t searchSimcHistFile( const TString &searchPattern, const TString &filename ){
    std::ifstream inputFile(filename.Data()); // Open the file for reading

    if (!inputFile.is_open()) {
      std::cerr << "Failed to open the file." << std::endl;
      std::cerr << filename.Data() << std::endl;
      return -99;
    }

    TString line;
    Double_t foundValue = -1.0; // Default value if not found, using a Double_t
    TString previousLine; // Store the previous line for exceptions
    while (line.ReadLine(inputFile)) {
      // Search for lines containing the specified search pattern
      if (line.Contains(searchPattern)) {
	// Exception for "luminosity" and other patterns with scientific notation
	if (line.EndsWith("ub^-1") || line.EndsWith("geV^2")) {
	  // Extract the numerical value using TString::Tokenize
	  TString delim("=");
	  TObjArray* tokens = line.Tokenize(delim);
	  if (tokens->GetEntries() >= 2) {
	    TObjString* valueStr = dynamic_cast<TObjString*>(tokens->At(1));
	    if (valueStr) {
	      TString valueString = valueStr->String();
	      // Convert the value to a Double_t
	      foundValue = valueString.Atof();
	      delete tokens; // Clean up the tokens
	      return foundValue;
	    }
	  }
	  delete tokens; // Clean up the tokens if no valid value was found
	} else {
	  // For other patterns, try to parse as a Double_t
	  TString valueString;
	  TString delim("=");
	  TObjArray* tokens = line.Tokenize(delim);
	  if (tokens->GetEntries() >= 2) {
	    TObjString* valueStr = dynamic_cast<TObjString*>(tokens->At(1));
	    if (valueStr) {
	      valueString = valueStr->String();
	    }
	  }
	  // Attempt to convert to a Double_t
	  if (valueString.IsFloat()) {
	    foundValue = valueString.Atof();
	    delete tokens; // Clean up the tokens
	    return foundValue;
	  }
	  delete tokens; // Clean up the tokens if no valid value was found
	}
      }
      // Store the previous line for exceptions
      previousLine = line;
    }

    // Handle the "luminosity" and other exceptions when no appropriate line is found
    if ((searchPattern == "luminosity" && previousLine.EndsWith("ub^-1")) ||
        (searchPattern == "genvol" )) {
      // Extract the numerical value using TString::Tokenize
      TString delim("=");
      TObjArray* tokens = previousLine.Tokenize(delim);
      if (tokens->GetEntries() >= 2) {
	TObjString* valueStr = dynamic_cast<TObjString*>(tokens->At(1));
	if (valueStr) {
	  TString valueString = valueStr->String();
	  // Convert the value to a Double_t
	  foundValue = valueString.Atof();
	}
      }
      delete tokens; // Clean up the tokens if no valid value was found
    }

    inputFile.close(); // Close the file

    return foundValue; // Return the default value if not found
  }

  //function to create a residual histogram courtesy jboyd
  TH1D *makeResidualHisto(TString identifier, TH1D *histo_1, TH1D *histo_2, bool match_bins = true, bool match_x_range = false ) {
    // Extract information from the first histogram
    TString histo_1_name = histo_1->GetName();
    Int_t Nbins_histo_1 = histo_1->GetNbinsX();
    Double_t MinX_histo_1 = histo_1->GetXaxis()->GetXmin();
    Double_t MaxX_histo_1 = histo_1->GetXaxis()->GetXmax();

    // Extract information from the second histogram
    TString histo_2_name = histo_2->GetName();
    Int_t Nbins_histo_2 = histo_2->GetNbinsX();
    Double_t MinX_histo_2 = histo_2->GetXaxis()->GetXmin();
    Double_t MaxX_histo_2 = histo_2->GetXaxis()->GetXmax();

    // Check for matching bin numbers if required
    if(match_bins && (Nbins_histo_1 != Nbins_histo_2)) {
      cout << "Histograms need to have matching bin numbers for accurate residual calculation." << endl;
      return nullptr; // Return null to indicate failure
    }

    // Check for matching x ranges if required
    if(match_x_range && ((MinX_histo_1 != MinX_histo_2) || (MaxX_histo_1 != MaxX_histo_2))) {
      cout << "Histograms need to have matching x ranges for accurate residual calculation." << endl;
      return nullptr; // Return null to indicate failure
    }

    // Create the residual histogram
    TH1D *resid_histo = new TH1D("resid_histo", Form("Residual Histogram - %s: %s and %s", identifier.Data(), histo_1_name.Data(), histo_2_name.Data()), Nbins_histo_1, MinX_histo_1, MaxX_histo_1);

    // Calculate residuals for each bin
    for(int bin = 1; bin <= Nbins_histo_1; bin++) { // Note: ROOT bin indices start from 1
      double content1 = histo_1->GetBinContent(bin);
      double content2 = histo_2->GetBinContent(bin);
      double residual = content1 - content2;

      resid_histo->SetBinContent(bin, residual);
    }

    return resid_histo;
  }

  TH1D *makeResidualHistoWithError(TString identifier, TH1D *histo_1, TH1D *histo_2, bool match_bins = true, bool match_x_range = false ){

    // Generate a unique identifier for the fit function name
    TString uniqueID = TString::Format("_%u", gRandom->Integer(1000000));

    TString histo_1_name = histo_1->GetName();
    Int_t Nbins_histo_1 = histo_1->GetNbinsX();
    Double_t MinX_histo_1 = histo_1->GetXaxis()->GetXmin();
    Double_t MaxX_histo_1 = histo_1->GetXaxis()->GetXmax();

    TString histo_2_name = histo_2->GetName();
    Int_t Nbins_histo_2 = histo_2->GetNbinsX();
    Double_t MinX_histo_2 = histo_2->GetXaxis()->GetXmin();
    Double_t MaxX_histo_2 = histo_2->GetXaxis()->GetXmax();

    // Ensure histograms have the same binning and x range if required
    if( match_bins && (Nbins_histo_1 != Nbins_histo_2) ){
      cout << "Histograms need to have matching bin numbers for accurate residual calculation." << endl;
      return nullptr; // Exiting if bin numbers do not match
    }

    if( match_x_range && ((MinX_histo_1 != MinX_histo_2) || (MaxX_histo_1 != MaxX_histo_2))){
      cout << "Histograms need to have matching x ranges for accurate residual calculation." << endl;
      return nullptr; // Exiting if x ranges do not match
    }

    // Creating the residual histogram
    TH1D *resid_histo = new TH1D(Form("resid_histo_%s",uniqueID.Data()), Form("Residual Histogram with Errors - %s: %s and %s", identifier.Data(), histo_1_name.Data(), histo_2_name.Data()), Nbins_histo_1, MinX_histo_1, MaxX_histo_1);

    for( int bin = 1; bin <= Nbins_histo_1; bin++ ){ // Note: bin indices in ROOT start from 1
      double value1 = histo_1->GetBinContent(bin);
      double error1 = histo_1->GetBinError(bin);
      double value2 = histo_2->GetBinContent(bin);
      double error2 = histo_2->GetBinError(bin);

      double residual = value1 - value2;
      double error = sqrt(error1*error1 + error2*error2); // Calculating the error for the residual

      resid_histo->SetBinContent(bin, residual);
      resid_histo->SetBinError(bin, error); // Setting the error for each bin in the residual histogram
    }

    return resid_histo;
  }

  //Fits a gaussian to a distribution twice, first course, then fine
  std::vector<Double_t> fitGaussianAndGetFineParams(TH1D* hist, Double_t sig, Double_t low = -1e38, Double_t high = 1e38) {
    
    std::vector<Double_t> params(3); // Vector to store amplitude, mean, and sigma

    if (!hist) {
      std::cerr << "Histogram is null!" << std::endl;
      return params;
    }

    // Find the bin numbers corresponding to the specified range. If default values are passed, real edge bins returned by FindBin().
    Int_t binLow = hist->FindBin(low);
    Int_t binHigh = hist->FindBin(high);

    // Initial values for maximum content and bin
    Double_t maxContent = 0;
    Int_t maxBin = -1;

    // Iterate over the bins in the range
    for (Int_t i = binLow; i <= binHigh; ++i) {
      Double_t content = hist->GetBinContent(i);
      if (content > maxContent) {
	maxContent = content;
	maxBin = i;
      }
    }

    Double_t xMax = hist->GetXaxis()->GetBinCenter(maxBin);

    // Define the fit range
    Double_t fitMin = xMax - sig;
    Double_t fitMax = xMax + sig;

    // Fit the histogram within the specified range
    TF1 *gausFit = new TF1("gausFit", "gaus", fitMin, fitMax);
    hist->Fit(gausFit, "RQ"); // "R" for fit range, "Q" for quiet mode (no print)

    // Get the mean from the fit
    Double_t amplitude = gausFit->GetParameter(0); // Parameter 0 is the amplitude of the Gaussian 
    Double_t mean = gausFit->GetParameter(1); // Parameter 1 is the mean of the Gaussian
    Double_t sigma = gausFit->GetParameter(2); // Parameter 2 is the std dev of the Gaussian

    // Clean up
    delete gausFit;

    // Fit the histogram within the fine range
    Double_t fit2Min = mean - 1*sigma;
    Double_t fit2Max = mean + 1*sigma;
    TF1 *gausFit_fine = new TF1("gausFit_fine", "gaus", fit2Min, fit2Max);
    gausFit_fine->SetParameters(amplitude,mean,sigma);
    hist->Fit(gausFit_fine, "RQ"); // "R" for fit range, "Q" for quiet mode (no print)

    // Store the parameters in the vector
    params[0] = gausFit_fine->GetParameter(0); // Amplitude
    //cout << gausFit_fine->GetParameter(0) << endl;
    params[1] = gausFit_fine->GetParameter(1); // Mean
    //cout << gausFit_fine->GetParameter(1) << endl;
    params[2] = gausFit_fine->GetParameter(2); // Sigma
    //cout << gausFit_fine->GetParameter(2) << endl;

    return params;
  }

  //Get Pearson correlation factor via covariance analysis on constrained range of TH2D
  double CalculateCorrelationFactor(TH2D* hist, double xmin, double xmax, double ymin, double ymax) {
    if (!hist) {
      std::cerr << "Histogram is null." << std::endl;
      return 0; // Error indication
    }

    double sumX = 0, sumY = 0, sumXY = 0, sumX2 = 0, sumY2 = 0, n = 0;
    double x, y, binContent;

    // First, compute the means of x and y within the specified bounds
    for (int ix = 1; ix <= hist->GetNbinsX(); ++ix) {
      x = hist->GetXaxis()->GetBinCenter(ix);
      if (x < xmin || x > xmax) continue;

      for (int iy = 1; iy <= hist->GetNbinsY(); ++iy) {
	y = hist->GetYaxis()->GetBinCenter(iy);
	if (y < ymin || y > ymax) continue;

	binContent = hist->GetBinContent(ix, iy);
	sumX += x * binContent;
	sumY += y * binContent;
	n += binContent;
      }
    }

    double meanX = sumX / n;
    double meanY = sumY / n;

    // Reset sums for the correlation calculation
    sumX = sumY = sumXY = sumX2 = sumY2 = 0;

    // Now, compute the components of the correlation formula
    for (int ix = 1; ix <= hist->GetNbinsX(); ++ix) {
      x = hist->GetXaxis()->GetBinCenter(ix);
      if (x < xmin || x > xmax) continue;

      for (int iy = 1; iy <= hist->GetNbinsY(); ++iy) {
	y = hist->GetYaxis()->GetBinCenter(iy);
	if (y < ymin || y > ymax) continue;

	binContent = hist->GetBinContent(ix, iy);
	double deltaX = x - meanX;
	double deltaY = y - meanY;
	sumXY += deltaX * deltaY * binContent;
	sumX2 += deltaX * deltaX * binContent;
	sumY2 += deltaY * deltaY * binContent;
      }
    }

    // Finally, calculate and return the correlation coefficient
    return sumXY / (sqrt(sumX2) * sqrt(sumY2));
  }

  void checkTH2DBinsAndEntries(TH2D *hist,
			       int minEntries,
			       int &firstBinWithEntries, 
			       int &lastBinWithEntries, 
			       int &firstBinEntries,
			       int &lastBinEntries){

    firstBinWithEntries = 1;
    lastBinWithEntries = hist->GetNbinsX();
    firstBinEntries = 0;
    lastBinEntries = 0;
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
      int entries = hist->ProjectionY("_py", i, i)->GetEntries();
      if (entries > minEntries) {
	firstBinWithEntries = i;
	firstBinEntries = entries;
	break;
      }
    }
    for (int i = hist->GetNbinsX(); i >= 1; --i) {
      int entries = hist->ProjectionY("_py", i, i)->GetEntries();
      if (entries > minEntries) {
	lastBinWithEntries = i;
	lastBinEntries = entries;
	break;
      }
    }

    return;
  }

  Double_t d_doubleGausPlusPol( Double_t *x, Double_t *par ) {
    // First Gaussian
    Double_t gauss1 = par[0] * exp(-0.5 * pow((x[0] - par[1]) / par[2], 2));

    // Second Gaussian
    Double_t gauss2 = par[3] * exp(-0.5 * pow((x[0] - par[4]) / par[5], 2));

    // Fourth-order polynomial
    Double_t poly = par[6] + par[7] * x[0] + par[8] * pow(x[0], 2) + par[9] * pow(x[0], 3) + par[10] * pow(x[0], 4);

    return gauss1 + gauss2 + poly;
  }

  // Function to perform 2x gaus + pol4 fit, calculate the gaus event ratio, and show the fit on a canvas
  void fitAndCalculateRatio(TH1D* hist, TCanvas* canvas, double& ratio, std::vector<double>& params, std::vector<double>& fitRange, bool fixpol=false) {
    if (!hist || !canvas) {
      std::cerr << "Error on fitAndCalculateRatio: Null pointer passed to the function." << std::endl;
      ratio = -1.0; // Indicate error
      return;
    }

    if (params.size() != 11) {
      std::cerr << "Error on fitAndCalculateRatio: Parameter vector must have 11 elements." << std::endl;
      ratio = -1.0; // Indicate error
      return;
    }

    if (fitRange.size() != 2) {
      std::cerr << "Error on fitAndCalculateRatio: Fit range vector must have 2 elements." << std::endl;
      ratio = -1.0; // Indicate error
      return;
    }

    // Create the combined Gaussian and polynomial fit function
    TF1* fitFunc = new TF1("fitFunc", d_doubleGausPlusPol, fitRange[0], fitRange[1], 11);

    // Set or fix parameters
    for (size_t i = 0; i < params.size(); ++i) {
      fitFunc->SetParameter(i, params[i]);
      // Prevent root from fitting a single peak twice
      if(i==1||i==4){
	fitFunc->SetParLimits(i,params[i]-params[i+1],params[i]+params[i+1]);
	cout << i << " " << params[i]-params[i+1] << " " << params[i]+params[i+1] << endl;
      }
      if(i==2||i==5)
	fitFunc->SetParLimits(i,0.10,0.16);
      if(i>5 && fixpol)
	fitFunc->FixParameter(i, params[i]);
    }

    // Perform the fit
    hist->Fit(fitFunc, "R");

    cout << "Actual fits: " << endl;

    for( int i=0; i<11; ++i )
      cout << fitFunc->GetParameter(i) << endl;

    // Calculate the total number of events under each Gaussian
    double totalLeft = 0.0;
    double totalRight = 0.0;
    int nBins = hist->GetNbinsX();
    for (int i = 1; i <= nBins; ++i) {
      double binCenter = hist->GetBinCenter(i);
      double binWidth = hist->GetBinWidth(i);

      // Evaluate each Gaussian component at the bin center
      double gauss1Value = fitFunc->EvalPar(&binCenter, &fitFunc->GetParameters()[0]);
      double gauss2Value = fitFunc->EvalPar(&binCenter, &fitFunc->GetParameters()[3]);

      // Subtract the polynomial contribution to isolate the Gaussian components
      double polyValue = fitFunc->EvalPar(&binCenter, &fitFunc->GetParameters()[6]);
      gauss1Value -= polyValue;
      gauss2Value -= polyValue;

      // Accumulate the total events considering the bin width
      totalLeft += gauss1Value * binWidth;
      totalRight += gauss2Value * binWidth;
    }

    // Calculate ratio
    ratio = (totalRight / totalLeft);

    // Draw the histogram and the fit
    canvas->cd();
    hist->Draw("hist");
    fitFunc->Draw("SAME");

    // Create separate functions for each Gaussian and the polynomial
    TF1* gauss1Func = new TF1("gauss1", "gaus", fitRange[0], fitRange[1]);
    gauss1Func->SetParameters(fitFunc->GetParameter(0), fitFunc->GetParameter(1), fitFunc->GetParameter(2));
    gauss1Func->SetLineColor(kRed);
    gauss1Func->Draw("SAME");

    TF1* gauss2Func = new TF1("gauss2", "gaus", fitRange[0], fitRange[1]);
    gauss2Func->SetParameters(fitFunc->GetParameter(3), fitFunc->GetParameter(4), fitFunc->GetParameter(5));
    gauss2Func->SetLineColor(kBlue);
    gauss2Func->Draw("SAME");

    TF1* polyFunc = new TF1("poly", "pol4", fitRange[0], fitRange[1]);
    polyFunc->SetParameters(&fitFunc->GetParameters()[6]);
    polyFunc->SetLineColor(kGreen);
    if(!fixpol)
      polyFunc->Draw("SAME");

    // Create a legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry( gauss1Func, "Proton Gaus Fit", "l");
    legend->AddEntry( gauss2Func, "Neutron Gaus Fit", "l");
    if(!fixpol)
      legend->AddEntry( polyFunc, "4th order poly BG", "l");
    std::ostringstream stream;
    stream << "Yield Left: " << totalLeft;
    legend->AddEntry((TObject*)0, stream.str().c_str(), "");
    stream.str(""); // Clear the stream
    stream << "Yield Right: " << totalRight;
    legend->AddEntry((TObject*)0, stream.str().c_str(), "");
    stream.str(""); // Clear the stream
    stream << "Ratio R/L: " << ratio;
    legend->AddEntry((TObject*)0, stream.str().c_str(), "");
    legend->Draw();

  }

  //skewed gaussian fit - duplicated from fits.C for use below
  Double_t d_sgfit(Double_t *x, Double_t *par){
    Double_t amp = par[0];
    Double_t offset = par[1];
    Double_t sigma = par[2];
    Double_t alpha = par[3];

    return amp*exp( -pow( x[0]-offset,2. )/( 2.*pow(sigma,2.) ) )*( 1+erf( (x[0]-offset)*alpha/sigma*sqrt(2.) ) );
  }

  void fitSkewedGaussianAndGetFineParams(TH1D* hist, Double_t sig, Double_t low, Double_t high, std::vector<Double_t> &params, Double_t alpha_ll = 0.0) {
    if (!hist) {
      throw std::runtime_error("Histogram is null!");
    }

    params.resize(4); // Resize to store amplitude, mean, sigma, and alpha

    // Coarse Fit
    hist->GetXaxis()->SetRangeUser(low, high);
    TF1 *gausFit = new TF1("gausFit", "gaus", low, high);
    hist->Fit(gausFit, "RQ");

    // Fine Fit
    Double_t amp = gausFit->GetParameter(0);
    Double_t mean = gausFit->GetParameter(1);
    Double_t sigma = gausFit->GetParameter(2);
    Double_t fit2Min = mean - sig * sigma;
    Double_t fit2Max = mean + sig * sigma;

    delete gausFit;

    TF1 *sg_fine = new TF1("sg_fine", d_sgfit, fit2Min, fit2Max, 4);
    sg_fine->SetParameters(amp, mean, sigma, 0); // Initial guess for alpha is 0
    sg_fine->SetParLimits(3,alpha_ll,1);
    hist->Fit(sg_fine, "RQ");

    // Store parameters
    for (int i = 0; i < 4; ++i) {
      params[i] = sg_fine->GetParameter(i);
    }

    delete sg_fine;
  }

  //takes bounds for side-band analysis and fits side-bands to pol4, subtracts this BG, returns events left and fit params
  std::pair<double, std::vector<double>> performSideBandAnalysis(TH1D* hist, double xLow, double xHigh) {
    std::vector<double> polyParams;
    double totalEvents = 0;

    if (!hist) {
      std::cerr << "Histogram is null!" << std::endl;
      return std::make_pair(-1, polyParams);
    }

    // Fit a 4th order polynomial excluding the range between xLow and xHigh
    TF1* background = new TF1("background", "pol4", hist->GetXaxis()->GetXmin(), hist->GetXaxis()->GetXmax());
    for (int i = hist->FindBin(xLow); i <= hist->FindBin(xHigh); ++i) {
      background->RejectPoint(i);
    }
    hist->Fit(background, "RQ");

    // Store the parameters of the polynomial
    for (int i = 0; i <= 4; ++i) {
      polyParams.push_back(background->GetParameter(i));
    }

    // Subtract the background and total the events
    for (int i = hist->FindBin(xLow); i <= hist->FindBin(xHigh); ++i) {
      double binCenter = hist->GetBinCenter(i);
      double binContent = hist->GetBinContent(i);
      double backgroundValue = background->Eval(binCenter);
      double subtractedValue = binContent - backgroundValue;

      if (subtractedValue > 0) {
	totalEvents += subtractedValue;
      }
    }

    delete background;

    return std::make_pair(totalEvents, polyParams);
  }

  //function that mirrors first half of a distribution onto the second half
  TH1D* MirrorHistogram(TH1D* originalHist) {
    if (!originalHist) {
      std::cerr << "Null histogram passed!" << std::endl;
      return nullptr;
    }

    // Get the number of bins and the maximum bin
    int nBins = originalHist->GetNbinsX();
    int maxBin = originalHist->GetMaximumBin();

    // Create a new histogram
    TH1D* mirroredHist = new TH1D(*originalHist); // Copy constructor
    mirroredHist->SetNameTitle("mirroredHist", "Mirrored Histogram");

    // Loop over the bins and mirror the first part onto the second part
    for (int i = maxBin + 1; i <= nBins; ++i) {
      int mirrorBinIndex = maxBin - (i - maxBin);
      double mirroredValue = originalHist->GetBinContent(mirrorBinIndex);
      mirroredHist->SetBinContent(i, mirroredValue);
    }

    return mirroredHist;
  }

  //script to look through directory and return one-to-one matches for .hist files and .root files. Intended for jboyd file structure, but lef generic for otherwise.
  void FindMatchingFiles(const std::string& directory1, 
			 const std::string& directory2, 
			 const std::string& partialName, 
			 std::vector<std::string>& vector1, 
			 std::vector<std::string>& vector2, 
			 bool jboyd) {
    try {
      if (!fs::is_directory(directory1) || !fs::is_directory(directory2)) {
	std::cerr << "Error: One or both directories are invalid." << std::endl;
	return;
      }

      // Step 1: Populate vector1 with all matching files from directory1
      for (const auto& entry : fs::directory_iterator(directory1)) {
	if (entry.is_regular_file() && entry.path().filename().string().find(partialName) != std::string::npos) {
	  vector1.push_back(entry.path().string());
	}
      }

      // Temporary vector to track indices of vector1 elements with no match
      std::vector<int> unmatchedIndices;
      int matches = 0;
      // Step 2 & 3: Try to find a match in directory2 for each file in vector1
      for (size_t i = 0; i < vector1.size(); ++i) {

	cout << "Looking for dir1 to dir2 matches, " << i << "/" << vector1.size() << " total matches " << matches << "\r";
	cout.flush();

	fs::path filePath1 = vector1[i];
	std::string filename1 = filePath1.filename().string();
            
	// Modify filename based on jboyd flag
	std::string expectedFilename2 = (jboyd ? "replayed_digitized_" : "replayed_") + filename1;
	expectedFilename2 = expectedFilename2.substr(0, expectedFilename2.find_last_of('.')) + ".root";

	bool matchFound = false;
	for (const auto& entry : fs::directory_iterator(directory2)) {
	  if (entry.is_regular_file() && entry.path().filename() == expectedFilename2) {
	    vector2.push_back(entry.path().string());
	    matchFound = true;
	    matches++;
	    break; // Match found, no need to continue searching
	  }
	}

	if (!matchFound) {
	  // Mark index for removal if no match is found
	  unmatchedIndices.push_back(i);
	}
      }

      // Step 4: Remove unpaired elements from vector1
      // Iterate in reverse to avoid messing up indices after removal
      for (auto it = unmatchedIndices.rbegin(); it != unmatchedIndices.rend(); ++it) {
	vector1.erase(vector1.begin() + *it);
      }
    } catch (const fs::filesystem_error& e) {
      std::cerr << "Error: " << e.what() << std::endl;
    }
  }

  void SyncFilesWithCsv(const std::string& directory1, // Path to .csv
			const std::string& directory2, // Path to .root
			const std::string& partialName, // Search word
			std::vector<std::string>& vector1, // Vector for .root file paths
			std::vector<std::pair<std::string, std::vector<float>>>& csvData) { // Vector to store CSV data with root file path
    try {
      // Populate vector1 with digitized and replayed .root files
      std::string replayed_partialName = "replayed_" + partialName;
      for (const auto& entry : fs::directory_iterator(directory2)) {
      	if (entry.is_regular_file() && entry.path().filename().string().find(replayed_partialName) != std::string::npos && entry.path().extension() == ".root") {
      	  vector1.push_back(entry.path().string());
      	}
      }

      // Populate vector1 with digitized and replayed .root files
      // for (const auto& entry : fs::directory_iterator(directory2)) {
      // 	if (entry.is_regular_file() && entry.path().filename().string().find(partialName) != std::string::npos && entry.path().extension() == ".root") {
      // 	  vector1.push_back(entry.path().string());
      // 	}
      // }

      // Read the CSV file
      fs::path csvPath;
      bool csvFound = false;
      for (const auto& entry : fs::directory_iterator(directory1)) {
	if (entry.is_regular_file() && entry.path().filename().string().find(partialName) != std::string::npos && entry.path().filename().string().find("_summary.csv") != std::string::npos) {
	  csvPath = entry.path();
	  csvFound = true;
	  break;
	}
      }

      if (!csvFound) {
	std::cerr << "CSV file matching " << partialName << " not found in directory1." << std::endl;
	return;
      }

      std::ifstream csvFile(csvPath);
      std::string line;
      while (getline(csvFile, line)) {
	if (line.empty() || line[0] == 'j') continue; // Skip header or empty lines
	std::istringstream ss(line);
	float jobid;
	ss >> jobid; // First column is jobid, now stored as a float
	char comma;
	std::vector<float> rowData = {jobid}; // Initialize vector with jobid as the first element
	float data;

	while (ss >> comma >> data) {
	  rowData.push_back(data);
	}
	// Add placeholder pair; actual path to be matched later
	csvData.push_back({"", rowData});
      }

      // Attempt to match .root files in vector1 with csvData based on job ID
      for (auto& dataPair : csvData) {
	std::regex jobIdPattern("job_([0-9]+)\\.root$");
	for (const auto& path : vector1) {
	  std::smatch matches;
	  if (std::regex_search(path, matches, jobIdPattern) && matches.size() > 1) {
	    int jobId = std::stoi(matches[1].str());
	    if (jobId == static_cast<int>(dataPair.second[0])) {
	      dataPair.first = path; // Match found, update path in csvData
	      break; // Stop searching after finding a match
	    }
	  }
	}
      }

      // Optionally, remove entries from csvData without a matched path
      csvData.erase(std::remove_if(csvData.begin(), csvData.end(), [](const auto& pair) { return pair.first.empty(); }), csvData.end());

    } catch (const fs::filesystem_error& e) {
      std::cerr << "Error: " << e.what() << std::endl;
    }
  }

  //Function to look through two vectors of strings and remove entries which don't have a corresponding job pair.
  void synchronizeJobNumbers(std::vector<std::string>& vecA, std::vector<std::string>& vecB) {
    std::regex jobNumberRegex("job([0-9]+)");
    std::smatch match;
    std::unordered_set<std::string> jobNumbersA, jobNumbersB;

    // Extract job numbers from vecA
    for (const auto& str : vecA) {
      if (std::regex_search(str, match, jobNumberRegex)) {
	jobNumbersA.insert(match[0]);
      }
    }

    // Extract job numbers from vecB
    for (const auto& str : vecB) {
      if (std::regex_search(str, match, jobNumberRegex)) {
	jobNumbersB.insert(match[0]);
      }
    }

    // Remove unpaired entries from vecA
    vecA.erase(std::remove_if(vecA.begin(), vecA.end(), 
			      [&](const std::string& str) {
				if (std::regex_search(str, match, jobNumberRegex) && jobNumbersB.find(match[0]) == jobNumbersB.end()) {
				  std::cout << "Removing from vecA: " << str << std::endl;
				  return true;
				}
				return false;
			      }), vecA.end());

    // Remove unpaired entries from vecB
    vecB.erase(std::remove_if(vecB.begin(), vecB.end(), 
			      [&](const std::string& str) {
				if (std::regex_search(str, match, jobNumberRegex) && jobNumbersA.find(match[0]) == jobNumbersA.end()) {
				  std::cout << "Removing from vecB: " << str << std::endl;
				  return true;
				}
				return false;
			      }), vecB.end());
  }

  //overload
  void synchronizeJobNumbers(std::vector<std::pair<std::string, std::vector<float>>>& metadata_p,
			     std::vector<std::pair<std::string, std::vector<float>>>& metadata_n) {
    std::regex jobNumberRegex("job_([0-9]+)\\.root$");
    std::smatch match;
    std::unordered_set<std::string> jobNumbersP, jobNumbersN;

    // Extract job numbers from proton metadata
    for (const auto& dataPair : metadata_p) {
      if (std::regex_search(dataPair.first, match, jobNumberRegex) && match.size() > 1) {
	jobNumbersP.insert(match[1].str());
      }
    }

    // Extract job numbers from neutron metadata
    for (const auto& dataPair : metadata_n) {
      if (std::regex_search(dataPair.first, match, jobNumberRegex) && match.size() > 1) {
	jobNumbersN.insert(match[1].str());
      }
    }

    // Remove unpaired entries from proton metadata
    metadata_p.erase(std::remove_if(metadata_p.begin(), metadata_p.end(),
				    [&](const std::pair<std::string, std::vector<float>>& dataPair) {
				      if (std::regex_search(dataPair.first, match, jobNumberRegex) && match.size() > 1) {
					if (jobNumbersN.find(match[1].str()) == jobNumbersN.end()) {
					  std::cout << "Removing from proton metadata: " << dataPair.first << std::endl;
					  return true;
					}
				      }
				      return false;
				    }), metadata_p.end());

    // Remove unpaired entries from neutron metadata
    metadata_n.erase(std::remove_if(metadata_n.begin(), metadata_n.end(),
				    [&](const std::pair<std::string, std::vector<float>>& dataPair) {
				      if (std::regex_search(dataPair.first, match, jobNumberRegex) && match.size() > 1) {
					if (jobNumbersP.find(match[1].str()) == jobNumbersP.end()) {
					  std::cout << "Removing from neutron metadata: " << dataPair.first << std::endl;
					  return true;
					}
				      }
				      return false;
				    }), metadata_n.end());
  }


  /*
  void synchronizeJobNumbers(std::vector<std::pair<std::string, std::vector<float>>>& vecA,
			     std::vector<std::pair<std::string, std::vector<float>>>& vecB) {
    // Extract job numbers from vecA
    std::unordered_set<float> jobNumbersA;
    for (const auto& item : vecA) {
      if (!item.second.empty()) {
	jobNumbersA.insert(item.second[0]);
      }
    }

    // Extract job numbers from vecB
    std::unordered_set<float> jobNumbersB;
    for (const auto& item : vecB) {
      if (!item.second.empty()) {
	jobNumbersB.insert(item.second[0]);
      }
    }

    // Remove unpaired entries from vecA
    vecA.erase(std::remove_if(vecA.begin(), vecA.end(),
			      [&](const std::pair<std::string, std::vector<float>>& item) {
				return item.second.empty() || jobNumbersB.find(item.second[0]) == jobNumbersB.end();
			      }), vecA.end());

    // Remove unpaired entries from vecB
    vecB.erase(std::remove_if(vecB.begin(), vecB.end(),
			      [&](const std::pair<std::string, std::vector<float>>& item) {
				return item.second.empty() || jobNumbersA.find(item.second[0]) == jobNumbersA.end();
			      }), vecB.end());
  }

  //overload
  void synchronizeJobNumbers(std::map<int, std::pair<std::string, std::vector<float>>>& csvData_n,
			     std::map<int, std::pair<std::string, std::vector<float>>>& csvData_p) {
    // Step 1: Extract job numbers (keys) from both maps
    std::unordered_set<int> jobNumbersA, jobNumbersB;
    for (const auto& pair : csvData_n) {
      jobNumbersA.insert(pair.first);
    }
    for (const auto& pair : csvData_p) {
      jobNumbersB.insert(pair.first);
    }

    // Step 2: Remove unpaired entries from csvData_n
    for (auto it = csvData_n.begin(); it != csvData_n.end();) {
      if (jobNumbersB.find(it->first) == jobNumbersB.end()) {
	std::cout << "Removing from csvData_n: Job ID " << it->first << std::endl;
	it = csvData_n.erase(it); // Erase and move to the next element
      } else {
	++it; // Move to the next element
      }
    }

    // Step 3: Remove unpaired entries from csvData_p
    for (auto it = csvData_p.begin(); it != csvData_p.end();) {
      if (jobNumbersA.find(it->first) == jobNumbersA.end()) {
	std::cout << "Removing from csvData_p: Job ID " << it->first << std::endl;
	it = csvData_p.erase(it); // Erase and move to the next element
      } else {
	++it; // Move to the next element
      }
    }
  }
  */
  //Function to take TH1D, randomly sample from it, and return statistically fluctuated sample histograms
  std::vector<TH1D*> createBootstrapSamples(TH1D* originalHist, int nBootstrapSamples, TCanvas *canvas) {
    std::vector<TH1D*> bootstrapHists;

    if (!originalHist) return bootstrapHists;

    // Get the number of bins and range of the histogram
    int nBins = originalHist->GetNbinsX();
    double xMin = originalHist->GetXaxis()->GetXmin();
    double xMax = originalHist->GetXaxis()->GetXmax();

    // Initialize random number generator
    TRandom3 rnd(0); // 0 to use the machine clock

    // Prepare the canvas
    canvas->Divide(2, 2); // Divide into 4 pads for the first 4 histograms

    for (int i = 0; i < nBootstrapSamples; ++i) {
      // Create a new histogram for each bootstrap sample
      TH1D* bootstrapHist = new TH1D(Form("bootstrapHist%d", i), "Bootstrap Sample;X;Frequency", nBins, xMin, xMax);

      // Fill the bootstrap histogram by randomly sampling the original histogram
      for (int j = 0; j < originalHist->GetEntries(); ++j) {
	double val = originalHist->GetRandom();
	bootstrapHist->Fill(val);
      }

      // Store the bootstrap histogram for further analysis
      bootstrapHists.push_back(bootstrapHist);

      // Plot the first four histograms
      if (i < 4) {
	canvas->cd(i + 1);
	bootstrapHist->Draw();
      }
    }

    return bootstrapHists;
  }

  // Utility function to shift every bin of a TH1D along the x-axis
  TH1D* shiftHistogramX(TH1D* originalHist, double shiftValue) {
    if (!originalHist) return nullptr;

    // Preserve the total number of entries in the histogram for further analysis
    double totalEntries = originalHist->GetEntries();

    // Create a new histogram with the same binning as the original
    TH1D *shiftedHist = (TH1D*)(originalHist->Clone("shiftedHist"));

    // Clear the contents of the cloned histogram
    shiftedHist->Reset();

    // Shift each bin
    for (int i = 1; i <= originalHist->GetNbinsX(); ++i) {
      // Calculate new bin center
      double oldBinCenter = originalHist->GetBinCenter(i);
      double oldBinContent = originalHist->GetBinContent(i);
      double oldBinError = originalHist->GetBinError(i);

      // Find the bin in the new histogram that corresponds to the new bin center
      int newBin = shiftedHist->FindBin(oldBinCenter + shiftValue);

      // Add the content and error to the new bin
      // Note: If multiple old bins shift into the same new bin, their contents and errors are added
      double newBinContent = shiftedHist->GetBinContent(newBin) + oldBinContent;
      double newBinError = sqrt(pow(shiftedHist->GetBinError(newBin), 2) + pow(oldBinError, 2));

      shiftedHist->SetBinContent(newBin, newBinContent);
      shiftedHist->SetBinError(newBin, newBinError);
    }

    // Restore the total number of entries
    shiftedHist->SetEntries(totalEntries);

    return shiftedHist;
  }

  //Function to constrain x range of histogram and remove underflow and overflow for cleaner analysis
  TH1D* cloneAndCutHistogram(TH1D* originalHist, double xMin, double xMax) {
    if (!originalHist) {
      std::cerr << "Original histogram is null!" << std::endl;
      return nullptr;
    }

    // Get bin info
    int binMin = originalHist->FindBin(xMin);
    int binMax = originalHist->FindBin(xMax);
    int dnBins = binMax - binMin + 1; // +1 to include both binMin and binMax

    TH1D *scaledHist = new TH1D("scaledHist", "scaledHist", dnBins, xMin, xMax);

    for (int i = binMin, j = 1; i <= binMax; ++i, ++j) {
      double binContent = originalHist->GetBinContent(i);
      double binError = originalHist->GetBinError(i);
      scaledHist->SetBinContent(j, binContent);
      scaledHist->SetBinError(j, binError);

    }

    return scaledHist;
  }

  //gets total error where measurements are independent
  double GetTotalQuadratureError(TH1D* hist) {
    if (!hist) {
      std::cerr << "GetTotalQuadratureError: Histogram pointer is null." << std::endl;
      return 0.0;
    }

    double totalErrorSquared = 0.0;

    // Loop over all bins (excluding underflow and overflow)
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
      double binError = hist->GetBinError(i);
      totalErrorSquared += binError * binError;
    }

    return sqrt(totalErrorSquared);
  }

  //defines a new TH1D which is the same as the input TH1D - input TF1 bin by bin within xrange, else zero
  TH1D* subtractFunctionFromHistogramWithinRange(const TH1D* hist, const TF1* func, double xmin, double xmax) {
    if (!hist || !func) {
      std::cerr << "Error: Null pointer passed to the function." << std::endl;
      return nullptr;
    }

    // Clone the input histogram
    TH1D* resultHist = (TH1D*)hist->Clone();

    // Loop over each bin of the cloned histogram
    int nBins = resultHist->GetNbinsX();
    for (int i = 1; i <= nBins; ++i) {
      double binCenter = resultHist->GetBinCenter(i);

      if (binCenter < xmin || binCenter > xmax) {
	// Set bin content to zero outside the specified range
	resultHist->SetBinContent(i, 0);
      } else {
	// Subtract the function value from the bin content within the specified range
	double funcValue = func->Eval(binCenter);
	double binContent = resultHist->GetBinContent(i);
	resultHist->SetBinContent(i, binContent - funcValue);
      }
    }

    return resultHist;
  }

  //function that subtracts a fit from a histogram and returns the total error 
  void subtractFunctionAndGetTotalAndError(TH1D* hist, TF1* func, double xMin, double xMax, double &total, double xRangeLow, double xRangeHigh, double &error, bool abs_sub = false) {
    if (!hist || !func) {
      throw std::runtime_error("Histogram or function is null!");
    }

    total = 0.0;
    error = 0.0;

    int binMin = hist->FindBin(xMin);
    int binMax = hist->FindBin(xMax);

    int firstbin = hist->FindBin(xRangeLow);
    int lastbin = hist->FindBin(xRangeHigh);

    double error_sig = 0.;
    double error_bg = 0.;
    int sig_bins = 0;
    int bg_bins = 0;

    for(int i = firstbin; i <= lastbin; ++i) {
      double binCenter = hist->GetBinCenter(i);
      double funcValue = func->Eval(binCenter);
      double binContent = hist->GetBinContent(i);
      double binErr = hist->GetBinError(i);

      double subtractedValue = binContent - funcValue;

      if( i >= binMin && i <= binMax ){
	error_sig += pow(binErr, 2);
	sig_bins++;
	if (subtractedValue > 0 && abs_sub==true) {
	  total += subtractedValue; 
	}else if( abs_sub==false){
	  total += subtractedValue;
	}
      } else {
	error_bg += pow(binErr, 2);
	bg_bins++;
      }
    }

    double error_bin_bg = bg_bins > 0 ? sqrt(error_bg) / bg_bins : 0;

    error = sqrt(error_sig) + error_bin_bg;

    return;
  }

  //Function that takes two TH1D histograms, normalizes them, and plots them together on a passed canvas
  void plotNormalizedTH1DsOnCanvas(TCanvas* canvas, TH1D* hist1, TH1D* hist2, const char* title = "Normalized Histograms", const char* x_label = "X-axis", const char* y_label = "Frequency") {
    // Check for null pointers
    if (!canvas || !hist1 || !hist2) {
      std::cerr << "Error: Null pointer passed." << std::endl;
      return;
    }

    // Clone the histograms to avoid modifying the original ones
    TH1D* h1 = (TH1D*)hist1->Clone();
    TH1D* h2 = (TH1D*)hist2->Clone();

    // Normalize both histograms to unity
    h1->Scale(1.0 / h1->Integral());
    h2->Scale(1.0 / h2->Integral());

    // Determine the maximum y-value between both histograms
    double max_y = std::max(h1->GetMaximum(), h2->GetMaximum());

    // Set different line colors for visibility
    h1->SetLineColor(kRed);
    h2->SetLineColor(kBlue);

    // Draw the histograms on the provided canvas
    canvas->cd();
    h1->Draw("HIST");
    h2->Draw("HIST SAME");

    // Adjust the y-axis range to include the maximum y-value
    h1->GetYaxis()->SetRangeUser(0, max_y * 1.1); // Adding 10% for padding

    // Set titles and labels
    h1->SetTitle(title);
    h1->GetXaxis()->SetTitle(x_label);
    h1->GetYaxis()->SetTitle(y_label);

    // Add a legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    legend->AddEntry(h1, hist1->GetName(), "l");
    legend->AddEntry(h2, hist2->GetName(), "l");
    legend->Draw();

    // Update the canvas to display the plots
    canvas->Update();
  }

  // Takes name and globalcut, parses cuts, displays them on canvas with name
  void parseAndDisplayCuts(const char* name, const char* cuts) {
    // Convert the input string to a ROOT TString for easy manipulation
    TString strCuts(cuts);
    
    // Parse the conditions separated by "&&" into a vector of TStrings
    TObjArray* tokens = strCuts.Tokenize("&&");
    if (!tokens) return;
    
    // Create a canvas to display the conditions
    TCanvas* c1 = new TCanvas("c1", name, 800, 600);
    
    // Create a TPaveText to hold and display the conditions
    TPaveText* pt = new TPaveText(0.1, 0.1, 0.9, 0.9);
    pt->AddText(name);
    pt->AddText(" "); // Adding an empty line for better readability
    
    // Loop over the parsed conditions and add each to the TPaveText
    for (int i = 0; i < tokens->GetEntries(); ++i) {
      TObjString* token = (TObjString*)tokens->At(i);
      TString condition = token->GetString();
      //condition.ReplaceAll("abs(", ""); // Optional: Simplify display by removing 'abs('
      //condition.ReplaceAll(")", "");    // Optional: Simplify display by removing ')'
      pt->AddText(condition.Data());
    }
    
    // Draw the TPaveText on the canvas
    pt->Draw();
    
    // Update the canvas to display the changes
    c1->Update();

    // Clean up
    //delete tokens; // Free memory used by the TObjArray
  }

  // Takes name and globalcut, parses cuts, displays them on canvas with name
  void parseAndDisplayCuts(const char* name, const char* cuts, TCanvas *c1) {
    // Convert the input string to a ROOT TString for easy manipulation
    TString strCuts(cuts);
    
    // Parse the conditions separated by "&&" into a vector of TStrings
    TObjArray* tokens = strCuts.Tokenize("&&");
    if (!tokens) return;
    
    // Create a TPaveText to hold and display the conditions
    TPaveText* pt = new TPaveText(0.1, 0.1, 0.9, 0.9);
    pt->AddText(name);
    pt->AddText(" "); // Adding an empty line for better readability
    
    // Loop over the parsed conditions and add each to the TPaveText
    for (int i = 0; i < tokens->GetEntries(); ++i) {
      TObjString* token = (TObjString*)tokens->At(i);
      TString condition = token->GetString();
      // Optional: Simplify display by removing 'abs(' and ')'
      //condition.ReplaceAll("abs(", "");
      //condition.ReplaceAll(")", "");
      pt->AddText(condition.Data());
    }
    
    // Draw the TPaveText on the canvas
    pt->Draw();
    
    // Update the canvas to display the changes
    c1->Update();

    // Clean up
    //delete pt;      // Free memory used by TPaveText
    //delete tokens;  // Free memory used by the TObjArray
  }

}
