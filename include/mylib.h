#ifndef mylib_h
#define mylib_h

ofstream outfile;

using namespace std;
const double alpha = 1 - 0.6827;
const double signal_scale = 0.01;
// == Debugging Mode
bool debug = false;
bool blind_SR = false;
int tag_year = 0;
TString Cycle_name;
TString Cycle_name_signal;
TString SingleMuon;
TString DoubleEG;
// == Call all needed maps
map<TString, TH1*> maphist;
map<TString, TH2D*> maphist2D;
map<TString, TGraph*> map_gr;
map<TString, TGraphAsymmErrors*> map_asym_gr;
map<TString, TFile*> mapfile;
map<TString, TCanvas*> mapcanvas;
map<TString, TPad*> mappad;
map<TString, THStack*> maphstack;
map<TString, TLegend*> maplegend;
map<TString, TH1F*> mapfunc;
map<TString, double> map_overflow;
map<TString, TLine*> mapline;
map<TString, TKey*> maphistcheck;
map<TString, TList*> maplist;
map<TString, std::vector<double> > map_bin_vector;
map<TString, std::vector<TString> > map_sample_names;
map<TString, std::vector<double> > map_syst_array;
map<TString, std::vector<double> > map_syst_table;


const int N_syst = 35;
TString systematics[N_syst] = {"central",
			       "ElectronScaleDown", "ElectronScaleUp",
			       "ElectronSmearDown", "ElectronSmearUp",
			       "JetsResDown", "JetsResUp",
			       "JetsScaleDown", "JetsScaleUp",
			       //"SD_JMR_Down", "SD_JMR_Up", // -- Because there is no W/Z tagging in the analysis
			       "SD_JMS_Down", "SD_JMS_Up",
			       "PUReweight_Down", "PUReweight_Up",
			       "MuonScaleDown", "MuonScaleUp",
			       "MuonSmearDown", "MuonSmearUp",
			       //"PDFNormDown", "PDFNormUp" 
			       "Prefire_Up", "Prefire_Down",
			       "MuonRecoSFUp", "MuonRecoSFDown",
			       "MuonIDSFUp", "MuonIDSFDown",
			       "MuonISOSFUp", "MuonISOSFDown",
			       "MuonTriggerSFUp", "MuonTriggerSFDown",
			       "ElectronRecoSFUp", "ElectronRecoSFDown",
			       "ElectronIDSFUp", "ElectronIDSFDown",
			       "ElectronTriggerSFUp", "ElectronTriggerSFDown",
			       "ZPtRwUp", "ZPtRwDown",
};

const int N_syst_comparison = 16;
TString systematics_comparison[N_syst] = {"ElectronScale",
					  "ElectronSmear",
					  "JetsRes",
					  "JetsScale",
					  //"SD_JMR_", // -- Because there is no W/Z tagging in the analysis
					  "SD_JMS_",
					  "PUReweight_",
					  "MuonScale",
					  //"MuonSmear",
					  "Prefire_",
					  "MuonRecoSF",
					  "MuonIDSF",
					  "MuonISOSF",
					  "MuonTriggerSF",
					  "ElectronRecoSF",
					  "ElectronIDSF",
					  "ElectronTriggerSF",
					  "ZPtRw",
};

vector<double> vx, vy, vexl, vexh, veyl, veyh;

TH1F * GetHist(TString hname){

  TH1F *h = NULL;
  std::map<TString, TH1F*>::iterator mapit = mapfunc.find(hname);
  if(mapit != mapfunc.end()) return mapit-> second;

  return h;

}

double get_DY_norm_SF(TString nameofhistogram){

  double DY_norm = 1.;

  // -- Flavor dependent
  /*
  if(tag_year == 2016){
    if(nameofhistogram.Contains("0AK8")){
      if(nameofhistogram.Contains("DiMu")){
	DY_norm = 0.812;
      }
      if(nameofhistogram.Contains("DiEle")){
	DY_norm = 0.962;
      }
    }
    if(nameofhistogram.Contains("1AK8")){
      if(nameofhistogram.Contains("DiMu")){
	DY_norm = 0.709;
      }
      if(nameofhistogram.Contains("DiEle")){
	DY_norm = 0.689;
      }
    }
    if(nameofhistogram.Contains("2AK8")){
      if(nameofhistogram.Contains("DiMu")){
	DY_norm = 0.604;
      }
      if(nameofhistogram.Contains("DiEle")){
	DY_norm = 0.511;
      }
    }
  }
  if(tag_year == 2017){
    if(nameofhistogram.Contains("0AK8")){
      if(nameofhistogram.Contains("DiMu")){
	DY_norm = 1.095;
      }
      if(nameofhistogram.Contains("DiEle")){
	DY_norm = 1.093;
      }
    }
    if(nameofhistogram.Contains("1AK8")){
      if(nameofhistogram.Contains("DiMu")){
	DY_norm = 1.056;
      }
      if(nameofhistogram.Contains("DiEle")){
	DY_norm = 1.062;
      }
    }
    if(nameofhistogram.Contains("2AK8")){
      if(nameofhistogram.Contains("DiMu")){
	DY_norm = 0.747;
      }
      if(nameofhistogram.Contains("DiEle")){
	DY_norm = 0.881;
      }
    }
  }
  if(tag_year == 2018){
    if(nameofhistogram.Contains("0AK8")){
      if(nameofhistogram.Contains("DiMu")){
	DY_norm = 1.127;
      }
      if(nameofhistogram.Contains("DiEle")){
	DY_norm = 1.087;
      }
    }
    if(nameofhistogram.Contains("1AK8")){
      if(nameofhistogram.Contains("DiMu")){
	DY_norm = 0.947;
      }
      if(nameofhistogram.Contains("DiEle")){
	DY_norm = 0.952;
      }
    }
    if(nameofhistogram.Contains("2AK8")){
      if(nameofhistogram.Contains("DiMu")){
	DY_norm = 0.997;
      }
      if(nameofhistogram.Contains("DiEle")){
	DY_norm = 0.957;
      }
    }
  }
  */
  // -- Flavor Independent
  if(tag_year == 2016){
    if(nameofhistogram.Contains("0AK8")){
      DY_norm = 0.875;
    }
    if(nameofhistogram.Contains("1AK8")){
      DY_norm = 0.701;
    }
    if(nameofhistogram.Contains("2AK8")){
      DY_norm = 0.583;
    }
  }
  if(tag_year == 2017){
    if(nameofhistogram.Contains("0AK8")){
      DY_norm = 1.090;
    }
    if(nameofhistogram.Contains("1AK8")){
      DY_norm = 1.042;
    }
    if(nameofhistogram.Contains("2AK8")){
      DY_norm = 0.837;
    }
  }
  if(tag_year == 2018){
    if(nameofhistogram.Contains("0AK8")){
      DY_norm = 1.093;
    }
    if(nameofhistogram.Contains("1AK8")){
      DY_norm = 0.946;
    }
    if(nameofhistogram.Contains("2AK8")){
      DY_norm = 0.959;
    }
  }

  return DY_norm;
}

vector<double> Get_Syst_Error(TString nameofhistogram, TString sample){
  
  vector<double> error_vector;
  error_vector.clear();
  if(debug) cout << nameofhistogram << ", " << sample << endl;  
  Int_t N_rebinned_x = mapfunc[nameofhistogram + Cycle_name + sample + "rebin"] -> GetNbinsX();
  
  for(Int_t j = 1; j < N_rebinned_x; j++){
    double current_bin_error =  mapfunc[nameofhistogram + Cycle_name + sample + "rebin"] -> GetBinError(j);
    for(int i_syst = 1; i_syst < N_syst; i_syst++){
      if(GetHist(nameofhistogram + systematics[i_syst] + sample + "rebin")){
	double diff = fabs(mapfunc[nameofhistogram + Cycle_name + sample + "rebin"] -> GetBinContent(j) - mapfunc[nameofhistogram + systematics[i_syst] + sample + "rebin"] -> GetBinContent(j));
	if(debug) cout << systematics[i_syst] << " : " << diff << endl;;
	current_bin_error = pow(current_bin_error, 2);
	current_bin_error += pow(diff, 2);
	current_bin_error = pow(current_bin_error, 0.5);
      }
    }
    error_vector.push_back(current_bin_error);
    if(debug) cout << mapfunc[nameofhistogram + Cycle_name + sample + "rebin"] -> GetBinContent(j) << ", syst : " << current_bin_error << endl;
  }

  return error_vector;
}

void Rebin_with_overflow(TString histname, int N_bin, double binx[]){

  if(debug) cout << "[Rebin_with_overflow]" << endl;

  mapfunc[histname + "rebin"] = (TH1F*)mapfunc[histname] -> Rebin(N_bin - 1, histname + "rebin", binx);

  double last_bin = 0.;
  last_bin =  mapfunc[histname + "rebin"] -> GetBinContent(N_bin - 1) + mapfunc[histname + "rebin"] -> GetBinContent(N_bin);
  mapfunc[histname + "rebin"] -> SetBinContent(N_bin - 1, last_bin);

}

void Rebin_with_overflow_limit(TString histname, int N_bin, double binx[]){
  //if(debug) cout << "[Rebin_with_overflow_limit]" << endl;
  mapfunc[histname] -> Integral();
  mapfunc[histname + "rebin"] = (TH1F*)mapfunc[histname] -> Clone();
  mapfunc[histname + "rebin"] = (TH1F*)mapfunc[histname + "rebin"] -> Rebin(N_bin - 1, histname + "rebin", binx);
  //if(debug) cout << "[Rebin_with_overflow_limit] " << histname + "rebin" << endl;
  
  //mapfunc[histname + "rebin"] -> Integral();
}

void Remove_Negative_Bins(TString current_histname){
  int N_bin = mapfunc[current_histname] -> GetNbinsX();
  for(int i = 1; i <= N_bin; i++){
    double current_bin_content = mapfunc[current_histname] -> GetBinContent(i);
    if (current_bin_content < 0.000001){
      //current_bin_content = 0.; // -- For table
      current_bin_content = 0.000001;
      double current_error = current_bin_content / 1.0;
      mapfunc[current_histname] -> SetBinContent(i, current_bin_content);
      //mapfunc[current_histname] -> SetBinError(i, 0.); // -- For table
    }
  }
}


void Save_syst_array(TString current_histname, TString systematics, TString cycle_and_sample, int N_bin, double binx[]){
  
  cout << "[Save_syst_array]" << endl;
  
  TString histname_Up   = current_histname + "_" + systematics + "Up"   + cycle_and_sample;
  TString histname_Down = current_histname + "_" + systematics + "Down" + cycle_and_sample;
  TString histname_central = current_histname + "_central" + cycle_and_sample + "rebin"; // -- central hists are already rebinned

  // -- booleans to check if there are up/down hists
  bool up_hist = false;
  bool down_hist = false;
  
  cout << "[Save_syst_array] : " << histname_Up << endl;
  if(!mapfunc[histname_Up]) return;
  Rebin_with_overflow(histname_Up, N_bin, binx);
  Rebin_with_overflow(histname_Down, N_bin, binx);
  Remove_Negative_Bins(histname_Up + "rebin");
  Remove_Negative_Bins(histname_Down + "rebin");

  //mapfunc[histname_Up + "rebin"] = (TH1F*)mapfunc[histname_Up] ->  Rebin(N_bin - 1, histname_Up + "tmp", binx);
  //mapfunc[histname_Down + "rebin"] = (TH1F*)mapfunc[histname_Down] ->  Rebin(N_bin - 1, histname_Down + "tmp", binx);

  double last_bin_up = 0.;
  double last_bin_down = 0.;

  last_bin_up = mapfunc[histname_Up + "rebin"] -> GetBinContent(N_bin - 1) + mapfunc[histname_Up + "rebin"] -> GetBinContent(N_bin);
  last_bin_down = mapfunc[histname_Down + "rebin"] -> GetBinContent(N_bin - 1) + mapfunc[histname_Down + "rebin"] -> GetBinContent(N_bin);

  
  map_syst_array[current_histname + systematics + "Up"].clear();
  map_syst_array[current_histname + systematics + "Down"].clear();
  //map_syst_array
  for(Int_t j = 1; j < N_bin; j++){
    double current_entry_Up = mapfunc[histname_Up + "rebin"] -> GetBinContent(j);
    double current_entry_Down =  mapfunc[histname_Down + "rebin"] ->GetBinContent(j);

    if(current_entry_Up > current_entry_Down){
      map_syst_array[current_histname + cycle_and_sample + systematics + "Up"].push_back(current_entry_Up);
      map_syst_array[current_histname + cycle_and_sample + systematics + "Down"].push_back(current_entry_Down);
    }
    else{
      map_syst_array[current_histname + cycle_and_sample + systematics + "Up"].push_back(current_entry_Down);
      map_syst_array[current_histname + cycle_and_sample + systematics + "Down"].push_back(current_entry_Up);
    }
  }
  cout << "[Save_syst_array] " << current_histname + cycle_and_sample + systematics + "Up" << endl;
  cout << "[Save_syst_array] " << current_histname + cycle_and_sample + systematics << "Up size : " << map_syst_array[current_histname + cycle_and_sample + systematics + "Up"].size() << endl;
  cout << "[Save_syst_array] " << current_histname + cycle_and_sample + systematics << "Down size : " << map_syst_array[current_histname + cycle_and_sample + systematics + "Down"].size() << endl;
}

void sum_syst_error(TString current_histname, TString cycle_and_sample, int N_bin){
  
  // -- get central entries
  std::vector<double> central_content;
  central_content.clear();
  
  for(int j = 1; j < N_bin; j++){
    double current_content = mapfunc[current_histname + "_central" + cycle_and_sample + "rebin"] -> GetBinContent(j);
    central_content.push_back(current_content);
  }
  
  cout << "[sum_syst_error] " << endl;
  map_syst_array[current_histname + cycle_and_sample + "Up"].clear();
  map_syst_array[current_histname + cycle_and_sample + "Down"].clear();

  for(int j = 1; j < N_bin; j++){
    double current_error_up = 0.;
    double current_error_down = 0.;
    
    double current_stat_error = mapfunc[current_histname + "_central" + cycle_and_sample + "rebin"] -> GetBinError(j);
    current_error_up = current_error_up + pow(current_stat_error, 2);
    current_error_down = current_error_down+ pow(current_stat_error, 2);
    
    for(int i_syst = 0; i_syst < N_syst_comparison; i_syst++){
      cout << "[sum_syst_error] for j = " << j << ", i_syst = " << i_syst << endl;

      TString current_syst = systematics_comparison[i_syst];
      cout << "[sum_syst_error] map_syst_array[" << current_histname + cycle_and_sample + current_syst << "Up].size() : " << map_syst_array[current_histname + cycle_and_sample + current_syst + "Up"].size() << endl;

      double current_up = fabs(central_content.at(j-1) - map_syst_array[current_histname + cycle_and_sample + current_syst + "Up"].at(j-1) );
      double current_down = fabs(central_content.at(j-1) - map_syst_array[current_histname + cycle_and_sample + current_syst + "Down"].at(j-1) );
      
      cout << "[sum_syst_error] current_content : " << central_content.at(j-1) << ", up : " << map_syst_array[current_histname + cycle_and_sample + current_syst + "Up"].at(j-1) << ", down : " << map_syst_array[current_histname + cycle_and_sample + current_syst + "Down"].at(j-1) << endl;
      
      current_error_up = current_error_up + pow(current_up, 2);
      current_error_down = current_error_down + pow(current_down, 2);
    }
    
    current_error_up = pow(current_error_up, 0.5);
    current_error_down = pow(current_error_down, 0.5);
    
    map_syst_array[current_histname + cycle_and_sample + "Up"].push_back(current_error_up);
    map_syst_array[current_histname + cycle_and_sample + "Down"].push_back(current_error_down);
  }

  cout << "[sum_syst_error] Up size : " << map_syst_array[current_histname + cycle_and_sample + "Up"].size() << endl;
  cout << "[sum_syst_error] Down size : " << map_syst_array[current_histname + cycle_and_sample + "Down"].size() << endl;

}

void Write_syst_error(TString current_histname,  TString systematics, TString Cycle_name, TString current_sample, int N_bin, double binx[]){
  
  cout << "[Write_syst_error]" << endl;

  TString cycle_and_sample = Cycle_name + current_sample;
  TString histname_Up   = current_histname + "_" + systematics + "Up"   + cycle_and_sample;
  TString histname_Down = current_histname + "_" + systematics + "Down" + cycle_and_sample;
  TString histname_central = current_histname + "_central" + cycle_and_sample + "rebin"; // -- central hists are already rebinned

  // -- booleans to check if there are up/down hists
  bool up_hist = false;
  bool down_hist = false;

  if(!mapfunc[histname_Up]) return;

  mapfunc[histname_Up + "rebin"] = (TH1F*)mapfunc[histname_Up] ->  Rebin(N_bin - 1, histname_Up + "rebin", binx);
  mapfunc[histname_Down + "rebin"] = (TH1F*)mapfunc[histname_Down] ->  Rebin(N_bin - 1, histname_Down + "rebin", binx);
  
  double last_bin_up = 0.;
  double last_bin_down = 0.;

  last_bin_up = mapfunc[histname_Up + "rebin"] -> GetBinContent(N_bin - 1) + mapfunc[histname_Up + "rebin"] -> GetBinContent(N_bin);
  last_bin_down = mapfunc[histname_Down + "rebin"] -> GetBinContent(N_bin - 1) + mapfunc[histname_Down + "rebin"] -> GetBinContent(N_bin);

  mapfunc[histname_Up + "rebin"] -> SetBinContent(N_bin - 1, last_bin_up);
  mapfunc[histname_Down + "rebin"] -> SetBinContent(N_bin - 1, last_bin_down);
  
  // -- ADD year on syst name for correlated systematics over years 
  TString name_up = current_histname + "_" + current_sample + "_" + systematics + "Up";
  TString name_down = current_histname + "_" + current_sample + "_" + systematics + "Down";
  // -- Uncorrelated Systematics
  if(systematics.Contains("TriggerSF") || systematics.Contains("JetsRes") || systematics.Contains("Prefire")){
    name_up = current_histname + "_" + current_sample + "_" + systematics + TString::Itoa(tag_year,10) + "Up";
    name_down = current_histname + "_" + current_sample + "_" + systematics +  TString::Itoa(tag_year,10) + "Down";
  }
  
  mapfunc[histname_Up + "rebin"] -> SetName(name_up);
  mapfunc[histname_Down + "rebin"] -> SetName(name_down);
  

  cout << "[Write_syst_error] histname_Up : " << histname_Up << endl;
  //mapfunc[histname_Up + "rebin"] -> Write();
  //mapfunc[histname_Down + "rebin"] -> Write();
  
}

void Write_syst_error_limit(TString current_histname,  TString systematics, TString Cycle_name, TString current_sample, int N_bin, double binx[]){
  
  //cout << "[Write_syst_error_limit]" << endl;

  TString cycle_and_sample = Cycle_name + current_sample;
  TString histname_Up   = current_histname + "_" + systematics + "Up"   + cycle_and_sample;
  TString histname_Down = current_histname + "_" + systematics + "Down" + cycle_and_sample;
  TString histname_central = current_histname + "_central" + cycle_and_sample + "rebin"; // -- central hists are already rebinned
  
  if(!mapfunc[histname_Up]) return;
  
  mapfunc[histname_Up + "rebin"] = (TH1F*)mapfunc[histname_central] -> Clone();
  mapfunc[histname_Down + "rebin"] = (TH1F*)mapfunc[histname_central] -> Clone();
  
  if(mapfunc[histname_Up]){
    mapfunc[histname_Up + "rebin"] = (TH1F*)mapfunc[histname_Up] -> Clone();
    mapfunc[histname_Up + "rebin"] = (TH1F*)mapfunc[histname_Up + "rebin"] -> Rebin(N_bin - 1, histname_Up + "rebin", binx); 
  }
  if(mapfunc[histname_Down]){
    mapfunc[histname_Down + "rebin"] = (TH1F*)mapfunc[histname_Down] -> Clone();
    mapfunc[histname_Down + "rebin"] = (TH1F*)mapfunc[histname_Down + "rebin"] -> Rebin(N_bin - 1, histname_Down + "rebin", binx);
  }

  // -- ADD year on syst name for correlated systematics over years
  TString name_up = current_histname + "_" + current_sample + "_" + systematics + "Up";
  TString name_down = current_histname + "_" + current_sample + "_" + systematics + "Down";
  // -- Uncorrelated Systematics
  if(systematics.Contains("TriggerSF")){
    name_up = current_histname + "_" + current_sample + "_" + systematics + TString::Itoa(tag_year,10) + "Up";
    name_down = current_histname + "_" + current_sample + "_" + systematics +  TString::Itoa(tag_year,10) + "Down";
  }
  
  mapfunc[histname_Up + "rebin"] -> SetName(name_up);
  mapfunc[histname_Down + "rebin"] -> SetName(name_down);
  
  //cout << "[Write_syst_error_limit] histname_Up : " << histname_Up << endl;
}


void Remove_Syst_bin_for_Empty_central(TString histname, TString systematics, TString cycle, TString sample){
  
  TString nameofhistogram = histname + "_central";
  TString cycle_and_sample = cycle + sample;
  TString histname_central = nameofhistogram + cycle_and_sample + "rebin";
  TString histname_Up = histname + "_" + systematics + "Up"   + cycle_and_sample + "rebin";
  TString histname_Down = histname + "_" + systematics + "Down"   + cycle_and_sample + "rebin";
  if(systematics == "PDF" || systematics == "Scale"){
    histname_Up = histname + "_" + sample + "_" + systematics + "Up";
    histname_Down = histname + "_" + sample + "_" + systematics + "Down";
  }

  int N_bin = mapfunc[histname_central] -> GetNbinsX();
  for(int i = 1; i <= N_bin; i++){
    double central_content = mapfunc[histname_central] -> GetBinContent(i);
    if(central_content < 0.0000011){

      mapfunc[histname_Up] -> SetBinContent(i, 0.000001);
      mapfunc[histname_Up] -> SetBinError(i, 0.);
      mapfunc[histname_Down] -> SetBinContent(i, 0.000001);
      mapfunc[histname_Down] -> SetBinError(i, 0.);

      // -- For table
      /*
      mapfunc[histname_Up] -> SetBinContent(i, 0.);
      mapfunc[histname_Up] -> SetBinError(i, 0.);
      mapfunc[histname_Down] -> SetBinContent(i, 0.);
      mapfunc[histname_Down] -> SetBinError(i, 0.);
      */    
    }
  }
}

void Remove_Syst_bin_for_Empty_syst(TString histname, TString systematics, TString cycle, TString sample){
  TString nameofhistogram = histname + "_central";
  TString cycle_and_sample = cycle + sample;
  TString histname_central = nameofhistogram + cycle_and_sample + "rebin";
  TString histname_Up = histname + "_" + systematics + "Up"   + cycle_and_sample + "rebin";
  TString histname_Down = histname + "_" + systematics + "Down"   + cycle_and_sample + "rebin";
  if(systematics == "PDF" || systematics == "Scale"){
    histname_Up = histname + "_" + sample + "_" + systematics + "Up";
    histname_Down = histname + "_" + sample + "_" + systematics + "Down";
  }
  
  int N_bin = mapfunc[histname_central] -> GetNbinsX();
  for(int i = 1; i <= N_bin; i++){
    double central_content = mapfunc[histname_central] -> GetBinContent(i);
    double central_error = mapfunc[histname_central] -> GetBinError(i); 
    if(central_content > 0.0000011){
      double up_content = mapfunc[histname_Up] -> GetBinContent(i);
      double down_content = mapfunc[histname_Down] -> GetBinContent(i);
      
      if(up_content < 0.0000011){
	mapfunc[histname_Up] -> SetBinContent(i, central_content);
	mapfunc[histname_Up] -> SetBinError(i, central_error);
      }
      if(down_content < 0.0000011){
	mapfunc[histname_Down] -> SetBinContent(i, central_content);
	mapfunc[histname_Down] -> SetBinError(i, central_error);
      }
    }
  }
}


void change_to_pseudo_data(TString current_histname){
  
  int N_bin = mapfunc[current_histname] -> GetNbinsX();
  for(int i = 1; i <= N_bin; i++){
    int current_bin_int =  mapfunc[current_histname] -> GetBinContent(i);
    double current_error = pow(current_bin_int, 0.5);
    mapfunc[current_histname] -> SetBinContent(i, current_bin_int);
    mapfunc[current_histname] -> SetBinError(i, current_error);
  }
}


void Estimate_PDF_Error(TString current_histname, TString Cycle_name, TString current_sample, int N_bin, double binx[]){
  
  cout << "[Estimate_PDF_Error]" << endl;
  TString cycle_and_sample = Cycle_name + current_sample;
  TString histname_central = current_histname + "_central" + cycle_and_sample + "rebin"; // -- central hists are already rebinned
  mapfunc[histname_central] -> SetBinContent(0, 0.);
  mapfunc[histname_central] -> SetBinContent(N_bin, 0.);
  double integ_central_hist = mapfunc[histname_central] -> Integral();
  
  // -- Rebin Hessian
  for(int i_err = 1; i_err < 101; i_err++){
    TString current_err_hist = current_histname + "_central_Hessian_" + TString::Itoa(i_err, 10) + cycle_and_sample;
    //cout << "current_err_hist : " << current_err_hist << endl;
    //mapfunc[current_err_hist + "rebin"] = (TH1F*)mapfunc[current_err_hist] -> Clone();
    //mapfunc[current_err_hist + "rebin"] = (TH1F*)mapfunc[current_err_hist + "rebin"] -> Rebin(N_bin, current_err_hist + "rebin", binx);
    
    mapfunc[current_err_hist + "rebin"] = (TH1F*)mapfunc[current_err_hist] -> Rebin(N_bin - 1, current_err_hist + "rebin", binx);
    double last_bin = 0.;
    last_bin = mapfunc[current_err_hist + "rebin"] -> GetBinContent(N_bin - 1) + mapfunc[current_err_hist + "rebin"] -> GetBinContent(N_bin);
    mapfunc[current_err_hist + "rebin"] -> SetBinContent(N_bin - 1, last_bin);
  }
  
  // -- Rebin Alpha_s
  for(int i_err = 0; i_err < 2; i_err++){
    TString current_err_hist = current_histname + "_central_AlphaS_" + TString::Itoa(i_err, 10) + cycle_and_sample;
    //mapfunc[current_err_hist + "rebin"]= (TH1F*)mapfunc[current_err_hist] -> Clone();
    //mapfunc[current_err_hist + "rebin"] = (TH1F*)mapfunc[current_err_hist + "rebin"] ->Rebin(N_bin, current_err_hist + "rebin", binx);
    
    mapfunc[current_err_hist + "rebin"] = (TH1F*)mapfunc[current_err_hist] -> Rebin(N_bin - 1, current_err_hist + "rebin", binx);
    double last_bin = 0.;
    last_bin = mapfunc[current_err_hist + "rebin"] -> GetBinContent(N_bin - 1) + mapfunc[current_err_hist + "rebin"] -> GetBinContent(N_bin);
    mapfunc[current_err_hist + "rebin"] -> SetBinContent(N_bin - 1, last_bin);
  }

  // -- Save PDF errors in a vector
  std::vector<double> PDF_errors;
  PDF_errors.clear();
  for(int i_bin = 1; i_bin < N_bin; i_bin++){
    double current_Hessian_err = 0.;
    double current_AlphaS_err = 0.;
    double current_PDF_error = 0.;
    
    // -- Sum Hessian Errors
    for(int i_err = 1; i_err < 101; i_err++){
      TString current_err_hist = current_histname + "_central_Hessian_" + TString::Itoa(i_err, 10) + cycle_and_sample + "rebin";
      double current_diff = mapfunc[histname_central] -> GetBinContent(i_bin) - mapfunc[current_err_hist] -> GetBinContent(i_bin);
      double current_relative_diff = fabs(current_diff / mapfunc[histname_central] -> GetBinContent(i_bin));
      current_diff = pow(current_diff, 2);
      if(current_relative_diff > 0.5) current_diff = 0;
      current_Hessian_err += current_diff;
    }
    
    // -- Sum AlphaS Errors
    TString hist_alpha_0 = current_histname + "_central_AlphaS_0" + cycle_and_sample + "rebin";
    TString hist_alpha_1 = current_histname + "_central_AlphaS_1" + cycle_and_sample + "rebin";
    current_AlphaS_err = 0.5 * ( mapfunc[hist_alpha_0] -> GetBinContent(i_bin) - mapfunc[hist_alpha_1] -> GetBinContent(i_bin) );
    current_AlphaS_err = pow(current_AlphaS_err, 2);

    current_PDF_error = current_Hessian_err + current_AlphaS_err;
    current_PDF_error = pow(current_PDF_error, 0.5);

    PDF_errors.push_back(current_PDF_error);
  }
  
  // -- Make PDF syst Up/Down Shape
  TString PDF_Up = current_histname + "_" + current_sample + "_PDFUp";
  TString PDF_Down = current_histname + "_" + current_sample + "_PDFDown";
  mapfunc[PDF_Up] = (TH1F*)mapfunc[histname_central] -> Clone();
  mapfunc[PDF_Down] = (TH1F*)mapfunc[histname_central] -> Clone();
  for(int i_bin = 1; i_bin < N_bin; i_bin++){
    double current_Up_content = mapfunc[histname_central] -> GetBinContent(i_bin) + PDF_errors.at(i_bin-1);
    double current_Down_content = mapfunc[histname_central] -> GetBinContent(i_bin) - PDF_errors.at(i_bin-1);
    current_Down_content = max(0.00000001, current_Down_content);
    mapfunc[PDF_Up] -> SetBinContent(i_bin, current_Up_content);
    mapfunc[PDF_Down] -> SetBinContent(i_bin, current_Down_content);
  }

  // -- Normalize syst hists
  /*
  double integ_up_hist = mapfunc[PDF_Up] -> Integral();
  double integ_down_hist = mapfunc[PDF_Down] -> Integral();
  mapfunc[PDF_Up] -> Scale(integ_central_hist / integ_up_hist);
  mapfunc[PDF_Down] -> Scale(integ_central_hist / integ_down_hist);

  cout << "cetral : " << integ_central_hist << ", PDF up : " << mapfunc[PDF_Up] -> Integral() << ", PDF down : " << mapfunc[PDF_Down] -> Integral() << endl; 
  */
  // -- Set Syst hist names
  mapfunc[PDF_Up] -> SetName(PDF_Up);
  mapfunc[PDF_Down] -> SetName(PDF_Down);
}

void Estimate_Scale_Error(TString current_histname, TString Cycle_name, TString current_sample, int N_bin, double binx[]){
  
  cout << "[Estimate_Scale_Error]" << endl;
  TString cycle_and_sample = Cycle_name + current_sample;
  TString histname_central = current_histname + "_central" + cycle_and_sample + "rebin"; // -- central hists are already rebinned
  double integ_central_hist = mapfunc[histname_central] -> Integral();
  
  // -- Rebin Scale Hists
  for(int i_err = 1; i_err < 9; i_err++){
    if(i_err == 5 || i_err == 7) continue;
    TString current_err_hist = current_histname + "_central_Scale_" + TString::Itoa(i_err, 10) + cycle_and_sample;
    //cout << "current_err_hist : " << current_err_hist << endl;
    //mapfunc[current_err_hist + "rebin"] = (TH1F*)mapfunc[current_err_hist] -> Clone();
    //mapfunc[current_err_hist + "rebin"] = (TH1F*)mapfunc[current_err_hist + "rebin"] -> Rebin(N_bin, current_err_hist + "rebin", binx);
    mapfunc[current_err_hist + "rebin"] = (TH1F*)mapfunc[current_err_hist] -> Rebin(N_bin - 1, current_err_hist + "rebin", binx);
    double last_bin = 0.;
    last_bin = mapfunc[current_err_hist + "rebin"] -> GetBinContent(N_bin - 1) + mapfunc[current_err_hist + "rebin"] -> GetBinContent(N_bin);
    mapfunc[current_err_hist + "rebin"] -> SetBinContent(N_bin - 1, last_bin);
  }

  // -- Save Scale errors in two vectors
  std::vector<double> scale_errors_up;
  std::vector<double> scale_errors_down;
  scale_errors_up.clear();
  scale_errors_down.clear();
  for(int i_bin = 1; i_bin < N_bin; i_bin++){
    std::vector<double> current_scale_errors;
    current_scale_errors.clear();
    
    // -- Save Bin Contents from Scale error hists
    for(int i_err = 1; i_err < 9; i_err++){
      if(i_err == 5 || i_err == 7) continue;
      TString current_err_hist = current_histname + "_central_Scale_" + TString::Itoa(i_err, 10) + cycle_and_sample + "rebin";
      double current_content = mapfunc[current_err_hist] -> GetBinContent(i_bin);
      current_scale_errors.push_back(current_content);
    }
    
    // -- Get Max/Min (Envelope) among Scale Errors
    double current_max_scale = 0.;
    double current_min_scale = (mapfunc[histname_central] -> GetMaximum()) * 100.;
    
    for(unsigned int i_scale = 0; i_scale < current_scale_errors.size(); i_scale++){
      double current_scale = current_scale_errors.at(i_scale);
      
      if(current_scale > current_max_scale) current_max_scale = current_scale;
      if(current_scale < current_min_scale) current_min_scale = current_scale;
    } 

    scale_errors_up.push_back(current_max_scale);
    scale_errors_down.push_back(current_min_scale);
  }// -- for N_bin

  // -- Make Scale syst Up/Down Shape
  TString Scale_Up = current_histname + "_" + current_sample + "_ScaleUp";
  TString Scale_Down = current_histname + "_" + current_sample + "_ScaleDown";
  mapfunc[Scale_Up] = (TH1F*)mapfunc[histname_central] -> Clone();
  mapfunc[Scale_Down] = (TH1F*)mapfunc[histname_central] -> Clone();
  for(int i_bin = 1; i_bin < N_bin; i_bin++){
    double current_Up_content = scale_errors_up.at(i_bin - 1);
    double current_Down_content = scale_errors_down.at(i_bin - 1);
    current_Down_content = max(0.00000001, current_Down_content);
    mapfunc[Scale_Up] -> SetBinContent(i_bin, current_Up_content);
    mapfunc[Scale_Down] -> SetBinContent(i_bin,current_Down_content);
  }

  // -- Normalize
  /*
  double integ_up_hist = mapfunc[Scale_Up] -> Integral();
  double integ_down_hist = mapfunc[Scale_Down] -> Integral();
  mapfunc[Scale_Up] -> Scale(integ_central_hist / integ_up_hist);
  mapfunc[Scale_Down] -> Scale(integ_central_hist / integ_down_hist);
  */
  mapfunc[Scale_Up] -> SetName(Scale_Up);
  mapfunc[Scale_Down] -> SetName(Scale_Down);

}

void Estimate_PDF_accep_Error(TString current_histname, TString Cycle_name, TString current_sample){
  cout << "[Estimate_PDF_Error]" << endl;
  TString cycle_and_sample = Cycle_name + current_sample;
  TString histname_central = current_histname + "_central" + cycle_and_sample;
  TString norm_central = "Norm_central" + cycle_and_sample;
  double central_SR_norm = mapfunc[histname_central] -> Integral();
  double central_norm = mapfunc[norm_central] -> GetBinContent(1);
  double central_accep = central_SR_norm/central_norm;
  cout << "[[Estimate_PDF_accep_Error]] central_norm : " << central_norm << endl;
  for(int i_err = 1; i_err < 101; i_err++){
    
    TString current_err_hist = current_histname + "_central_Hessian_" + TString::Itoa(i_err, 10) + cycle_and_sample;
    double current_norm = mapfunc[current_err_hist] -> Integral();
  }
  

}

void Proper_error_data(TString nameofhistogram, TString current_data, int N_bin, double binx[]){
  
  //vector<double> vx, vy, vexl, vexh, veyl, veyh;
  
  map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] = new TGraphAsymmErrors(GetHist(nameofhistogram + Cycle_name + current_data + "rebin"));
  
  for(int i = 0; i < N_bin; i++){
    int N = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetBinContent(i + 1);
    double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
    double U =  (N==0) ? ( ROOT::Math::gamma_quantile_c(alpha,N+1,1) ) : ( ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) );
    if( N!=0 ){
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEYlow(i, (N-L) ); // / (binx[i + 1] - binx[i]) );
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEXlow(i, 0);
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEYhigh(i, (U-N) );// / (binx[i + 1] - binx[i]) );
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEXhigh(i, 0);
      
      double current_x = -1., current_y = -1.;
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> GetPoint(i, current_x, current_y);
      double modified_y = -1.;
      modified_y = current_y; // ( binx[i + 1] - binx[i] );
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPoint(i, current_x, modified_y);

      
      //if(debug) cout << "i : " << i << ", current_x : " << current_x << ", current_y : " << current_y << ", modified_y : " << modified_y << endl;

      veyl.push_back( (N-L)); // / (binx[i + 1] - binx[i]) );
      veyh.push_back( (U-N)); // / (binx[i + 1] - binx[i]));
    }
    else{
      double zerodata_err_low = 0.1 ;// / (binx[i + 1] - binx[i]);
      double zerodata_err_high = 1.8 ;// / (binx[i + 1] - binx[i]);

      veyl.push_back(zerodata_err_low);
      veyh.push_back(zerodata_err_high);
      
      double current_x = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetBinCenter(i + 1);
      double current_ex = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetBinWidth(i + 1) /2.;

      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEYlow(i, zerodata_err_low);
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEXlow(i, 0.);
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEYhigh(i, zerodata_err_high);
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEXhigh(i, 0.);

      vx.push_back(current_x);
      vexl.push_back(current_ex);
      vexh.push_back(current_ex);
    }
  }//end for i on N_bin
    
}

TGraphAsymmErrors* hist_to_graph(TH1D* hist, bool YErrorZero=false){

  TH1::SetDefaultSumw2(true);

  int Nbins = hist->GetXaxis()->GetNbins();
  double x[Nbins], y[Nbins], xlow[Nbins], xup[Nbins], ylow[Nbins], yup[Nbins];
  TAxis *xaxis = hist->GetXaxis();
  for(Int_t i=0; i<Nbins; i++){
    x[i] = xaxis->GetBinCenter(i+1);
    y[i] = hist->GetBinContent(i+1);
    xlow[i] = xaxis->GetBinCenter(i+1)-xaxis->GetBinLowEdge(i+1);
    xup[i] = xaxis->GetBinUpEdge(i+1)-xaxis->GetBinCenter(i+1);
    ylow[i] = hist->GetBinError(i+1);
    yup[i] = hist->GetBinError(i+1);
    if(YErrorZero){
      ylow[i] = 0;
      yup[i] = 0;
    }
    //ylow[i] = 0;
    //yup[i] = 0;
    //cout << "x = " << x[i] << ", y = " << y[i] << ", x_low = " << xlow[i] << ", xup = " << xup[i] << ", ylow = " << ylow[i] << ", yup = " << yup[i] << endl;
  }
  TGraphAsymmErrors *out = new TGraphAsymmErrors(Nbins, x, y, xlow, xup, ylow, yup);
  out->SetTitle("");
  return out;

}

TGraphAsymmErrors* hist_to_graph(TH1D* hist, int n_skip_x_left){

  TH1::SetDefaultSumw2(true);

  int Nbins = hist->GetXaxis()->GetNbins()-n_skip_x_left;
  double x[Nbins], y[Nbins], xlow[Nbins], xup[Nbins], ylow[Nbins], yup[Nbins];
  TAxis *xaxis = hist->GetXaxis();
  for(Int_t i=1; i<=Nbins; i++){
    x[i-1] = xaxis->GetBinCenter(i+n_skip_x_left);
    y[i-1] = hist->GetBinContent(i+n_skip_x_left);
    xlow[i-1] = xaxis->GetBinCenter(i+n_skip_x_left)-xaxis->GetBinLowEdge(i+n_skip_x_left);
    xup[i-1] = xaxis->GetBinUpEdge(i+n_skip_x_left)-xaxis->GetBinCenter(i+n_skip_x_left);
    ylow[i-1] = hist->GetBinError(i+n_skip_x_left);
    yup[i-1] = hist->GetBinError(i+n_skip_x_left);
    //ylow[i-1] = 0;
    //yup[i-1] = 0;
    //cout << "x = " << x[i-1] << ", y = " << y[i-1] << ", x_low = " << xlow[i-1] << ", xup = " << xup[i-1] << ", ylow = " << ylow[i-1] << ", yup = " << yup[i-1] << endl;
  }
  TGraphAsymmErrors *out = new TGraphAsymmErrors(Nbins, x, y, xlow, xup, ylow, yup);
  out->SetTitle("");
  return out;

}

TGraphAsymmErrors* hist_to_graph(TH1D* hist, int n_skip_x_left, int n_x_shift, int i_x_shift){

  TH1::SetDefaultSumw2(true);

  int Nbins = hist->GetXaxis()->GetNbins()-n_skip_x_left;
  double x[Nbins], y[Nbins], xlow[Nbins], xup[Nbins], ylow[Nbins], yup[Nbins];
  TAxis *xaxis = hist->GetXaxis();
  for(Int_t i=1; i<=Nbins; i++){
    x[i-1] = xaxis->GetBinCenter(i+n_skip_x_left);
    y[i-1] = hist->GetBinContent(i+n_skip_x_left);
    xlow[i-1] = xaxis->GetBinCenter(i+n_skip_x_left)-xaxis->GetBinLowEdge(i+n_skip_x_left);
    xup[i-1] = xaxis->GetBinUpEdge(i+n_skip_x_left)-xaxis->GetBinCenter(i+n_skip_x_left);
    ylow[i-1] = hist->GetBinError(i+n_skip_x_left);
    yup[i-1] = hist->GetBinError(i+n_skip_x_left);

    double dx = (xaxis->GetBinUpEdge(i+n_skip_x_left)-xaxis->GetBinLowEdge(i+n_skip_x_left))/double(n_x_shift+1);
    x[i-1] = xaxis->GetBinLowEdge(i+n_skip_x_left) + double(i_x_shift+1) * dx;
    xlow[i-1] = double(i_x_shift+1) * dx;
    xup[i-1] = xaxis->GetBinUpEdge(i+n_skip_x_left)-x[i-1];

    //ylow[i-1] = 0;
    //yup[i-1] = 0;
    //cout << "x = " << x[i-1] << ", y = " << y[i-1] << ", x_low = " << xlow[i-1] << ", xup = " << xup[i-1] << ", ylow = " << ylow[i-1] << ", yup = " << yup[i-1] << endl;
  }
  TGraphAsymmErrors *out = new TGraphAsymmErrors(Nbins, x, y, xlow, xup, ylow, yup);
  out->SetTitle("");
  return out;

}


TGraphAsymmErrors* GraphSubtract(TGraphAsymmErrors *a, TGraphAsymmErrors *b, bool Rel){

  //==== do a-b

  int NX = a->GetN();
  TGraphAsymmErrors* gr_out = (TGraphAsymmErrors*)a->Clone();

  for(int i=0; i<NX; i++){

    double a_x, a_y, b_x, b_y;

    a->GetPoint(i, a_x, a_y);
    b->GetPoint(i, b_x, b_y);

    if(Rel==true){
      gr_out->SetPoint( i, a_x, (a_y-b_y)/b_y );
    }
    else{
      gr_out->SetPoint( i, a_x, a_y-b_y );
    }

  }

  return gr_out;

}

void RemoveLargeError(TGraphAsymmErrors *a){

  int NX = a->GetN();

  for(int i=0; i<NX; i++){

    double x, y, yerr_low, yerr_high;

    a->GetPoint(i, x, y);
    yerr_low  = a->GetErrorYlow(i);
    yerr_high = a->GetErrorYhigh(i);

    if(y+yerr_high==1.){
      a->SetPointEYhigh( i, yerr_low );
    }

  }

}

void ScaleGraph(TGraphAsymmErrors *a, double c){

  int NX = a->GetN();

  for(int i=0; i<NX; i++){

    double x, y, yerr_low, yerr_high;

    a->GetPoint(i, x, y);
    yerr_low  = a->GetErrorYlow(i);
    yerr_high = a->GetErrorYhigh(i);

    a->SetPoint(i, x, c*y);
    a->SetPointEYlow(i, c*yerr_low);
    a->SetPointEYhigh(i, c*yerr_high);

  }

}



double GetMaximum(TH1D* hist){

  TAxis *xaxis = hist->GetXaxis();

  double maxval(-1.);
  for(int i=1; i<=xaxis->GetNbins(); i++){
    if( hist->GetBinContent(i) + hist->GetBinError(i) > maxval ){
      maxval = hist->GetBinContent(i) + hist->GetBinError(i);
    }
  }

  return maxval;

}

double GetMaximum(TGraphAsymmErrors *a){

  int NX = a->GetN();

  double maxval(-9999.);
  for(int i=0; i<NX; i++){

    double x, y, yerr_low, yerr_high;

    a->GetPoint(i, x, y);
    yerr_low  = a->GetErrorYlow(i);
    yerr_high = a->GetErrorYhigh(i);

    if( y+yerr_high > maxval ){
      maxval = y+yerr_high;
    }

  }

  return maxval;

}

double GetYieldSystematics(TH1D *hist){

  int n_syst = hist->GetXaxis()->GetNbins();
  int n_source = (n_syst-1)/2;

  //==== Bin index
  //==== i=1 : central
  //==== i=2 : _MuonEn_up
  //==== i=3 : _MuonEn_down
  //==== -> n_syst = 3
  //==== -> n_source = (n_syst-1)/2 = (3-1)/2 = 1

  double y_central = hist->GetBinContent(1);

  double sum_syst = 0.;
  for(int i=1; i<=n_source; i++){
    double yield_up = hist->GetBinContent(i*2);
    double yield_down = hist->GetBinContent(i*2+1);

    double syst_up = fabs(yield_up-y_central);
    double syst_down = fabs(yield_down-y_central);

    sum_syst += 0.5*(syst_up*syst_up+syst_down*syst_down);

  }

  return sqrt(sum_syst);

}

TDirectory *MakeTemporaryDirectory(){

  gROOT->cd();
  TDirectory* tempDir = 0;
  int counter = 0;
  while (not tempDir) {
    // First, let's find a directory name that doesn't exist yet:
    std::stringstream dirname;
    dirname << "HNCommonLeptonFakes_%i" << counter;
    if (gROOT->GetDirectory((dirname.str()).c_str())) {
      ++counter;
      continue;
    }
    // Let's try to make this directory:
    tempDir = gROOT->mkdir((dirname.str()).c_str());

  }

  return tempDir;

}

void AddPhantomZero(double a, TString align, int digit_int, int digit_frac){

  if(align=="r"){
    int number_maxdigit = 0;
    for(int i=10; i>=0; i--){
      if(a/pow(10,i)>=1.){
        number_maxdigit = i;
        break;
      }
    }
    //cout << "number_maxdigit = " << number_maxdigit << endl;
    for(int i=0; i<digit_int-(number_maxdigit+1); i++) cout << "\\phantom{0}";
    cout << std::fixed<<std::setprecision(digit_frac) << a;
  }

  else if(align=="l"){
    int target_total_digit = digit_int+digit_frac;
    int number_maxdigit = 0;
    for(int i=10; i>=0; i--){
      if(a/pow(10,i)>=1.){
        number_maxdigit = i;
        break;
      }
    }
    //cout << "target_total_digit = " << target_total_digit << endl;
    //cout << "number_maxdigit = " << number_maxdigit << endl;
    //cout << "--> # of \\phantom{0} = " << target_total_digit-digit_int-(number_maxdigit+1) << endl;

    cout << std::fixed<<std::setprecision(digit_frac) << a;
    for(int i=0; i<target_total_digit-(number_maxdigit+1)-digit_frac; i++) cout << "\\phantom{0}";
  }

}

vector<int> GetZPMassRange(int mN){

  vector<int> this_masses;

  if(mN==100){
    this_masses = {
      400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==200){
    this_masses = {
      600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==300){
    this_masses = {
      800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==400){
    this_masses = {
      1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==500){
    this_masses = {
      1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==600){
    this_masses = {
      1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==700){
    this_masses = {
      1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==800){
    this_masses = {
      1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==900){
    this_masses = {
      2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1000){
    this_masses = {
      2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1100){
    this_masses = {
      2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1200){
    this_masses = {
      2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1300){
    this_masses = {
      2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1400){
    this_masses = {
      3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1500){
    this_masses = {
      3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1600){
    this_masses = {
      3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1700){
    this_masses = {
      3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1800){
    this_masses = {
      3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1900){
    this_masses = {
      4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==2000){
    this_masses = {
      4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==2100){
    this_masses = {
      4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==2200){
    this_masses = {
      4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==2300){
    this_masses = {
      4800, 5000,
    };
    return this_masses;
  }
  if(mN==2400){
    this_masses = {
      5000,
    };
    return this_masses;
  }
  
  return this_masses;

}

vector<int> GetHNMassRange(int mZP, bool useallmass=false){

  vector<int> this_masses;

  if(useallmass){
    if(mZP==400){
      this_masses = {
        100,
      };
      return this_masses;
    }
    if(mZP==600){
      this_masses = {
        100, 200,
      };
      return this_masses;
    }
    if(mZP==800){
      this_masses = {
        100, 200, 300,
      };
      return this_masses;
    }
    if(mZP==1000){
      this_masses = {
        100, 200, 300, 400,
      };
      return this_masses;
    }
    if(mZP==1200){
      this_masses = {
        100, 200, 300, 400, 500,
      };
      return this_masses;
    }
    if(mZP==1400){
      this_masses = {
        100, 200, 300, 400, 500, 600,
      };
      return this_masses;
    }
    if(mZP==1600){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700,
      };
      return this_masses;
    }
    if(mZP==1800){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800,
      };
      return this_masses;
    }
    if(mZP==2000){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900,
      };
      return this_masses;
    }
    if(mZP==2200){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
      };
      return this_masses;
    }if(mZP==1400){
      this_masses = {
        100, 200, 300, 400, 500, 600,
      };
      return this_masses;
    }
    if(mZP==1600){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700,
      };
      return this_masses;
    }
    if(mZP==1800){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800,
      };
      return this_masses;
    }
    if(mZP==2000){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900,
      };
      return this_masses;
    }
    if(mZP==2200){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
      };
      return this_masses;
    }
    if(mZP==2400){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100,
      };
      return this_masses;
    }
    if(mZP==2600){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200,
      };
      return this_masses;
    }
    if(mZP==2800){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300,
      };
      return this_masses;
    }
    if(mZP==3000){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400,
      };
      return this_masses;
    }
    if(mZP==3200){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500,
      };
      return this_masses;
    }
    if(mZP==3400){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600,
      };
      return this_masses;
    }
    if(mZP==3600){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700,
      };
      return this_masses;
    }
    if(mZP==3800){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800,
      };
      return this_masses;
    }
    if(mZP==4000){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900,
      };
      return this_masses;
    }
    if(mZP==4200){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
      };
      return this_masses;
    }
    if(mZP==4400){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100,
      };
      return this_masses;
    }
    if(mZP==4600){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200,
      };
      return this_masses;
    }
    if(mZP==4800){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300,
      };
      return this_masses;
    }
    if(mZP==5000){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400,
      };
      return this_masses;
    }
  }
  else{
    if(mZP==400){
      this_masses = {
        100,
      };
      return this_masses;
    }
    if(mZP==800){
      this_masses = {
	100, 300,
      };
      return this_masses;
    }
    if(mZP==1200){
      this_masses = {
	//100, 300, 500,
	100, 500,
      };
      return this_masses;
    }
    if(mZP==1600){
      this_masses = {
	//100, 300, 500, 700,
	100, 700,
      };
      return this_masses;
    }
    if(mZP==2000){
      this_masses = {
	//100, 300, 500, 700, 900,
	100, 900,
      };
      return this_masses;
    }
    if(mZP==2400){
      this_masses = {
	//100, 300, 500, 700, 900, 1100,
	100, 1100,
      };
      return this_masses;
    }
    if(mZP==2800){
      this_masses = {
	//100, 300, 500, 700, 900, 1100, 1300,
	100, 1300,
      };
      return this_masses;
    }
    if(mZP==3200){
      this_masses = {
	//100, 300, 500, 700, 900, 1100, 1300, 1500,
	100, 1500,
      };
      return this_masses;
    }
    if(mZP==3600){
      this_masses = {
	//100, 300, 500, 700, 900, 1100, 1300, 1500, 1700,
	100, 700, 1700,
      };
      return this_masses;
    }
    if(mZP==4000){
      this_masses = {
	//100, 300, 500, 700, 900, 1100, 1300, 1500, 1700, 1900,
	100, 900, 1900,
      };
      return this_masses;
    }
    if(mZP==4400){
      this_masses = {
	//100, 300, 500, 700, 900, 1100, 1300, 1500, 1700, 1900, 2100,
	100, 1100, 2100,
      };
      return this_masses;
    }
    if(mZP==4800){
      this_masses = {
        //100, 300, 500, 700, 900, 1100, 1300, 1500, 1700, 1900, 2100, 2300,
	100, 1300, 2300,
      };
      return this_masses;
    }


  }

  cout << "[GetHNMassRange(int mZP)] Wrong mZP!!" << endl;
  return this_masses;
}

bool IsSignalCATSamaple(TString samplename){

  if(samplename.Contains("HNpair")) return true;

  return false;

}

int GetLowMN(double current_mN){
  current_mN = current_mN / 100.;
  int low_mN = current_mN;
  low_mN = low_mN * 100;

  return low_mN;
}



#endif
