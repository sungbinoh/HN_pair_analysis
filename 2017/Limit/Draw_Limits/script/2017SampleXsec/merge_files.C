#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TChain.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TKey.h>
#include <TList.h>
#include <TObject.h>
#include <TCollection.h>
#include <TLatex.h>
#include <Math/DistFunc.h>

using namespace std;

map<TString, TH1*> maphist;
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

//cycle name
TString Cycle_name = "HN_pair_all";

//sample names 
TString SingleMuon = "SingleMuon";
TString DoubleEG = "DoubleEG";

//DY
TString DY_high = "DYJets";

//VV 
TString WW = "WWTo2L2Nu_powheg";
TString WZ_2L = "WZTo2L2Q";
TString WZ_3L = "WZTo3LNu";
TString ZZ_2L = "ZZTo2L2Q";
TString ZZ_4L = "ZZTo4L_powheg";

TString ttbar = "TTLL_powheg";


const double alpha = 1 - 0.6827;

bool debug = true;

///////////////////////
// -- Getting Histogram
/////////////////////// 
TH1F * GetHist(TString hname){

  TH1F *h = NULL;
  std::map<TString, TH1F*>::iterator mapit = mapfunc.find(hname);
  if(mapit != mapfunc.end()) return mapit-> second;

  return h;
  
}

////////////////////
// -- Open ROOT file
////////////////////
void openfile(TString cyclename, TString samplename, TString dir, TString histname){

  TString filename = cyclename + "_" + samplename + ".root";
  cout << "[[openfile ]]opening : " << filename << endl;

  TFile *current_file = new TFile ((filename)) ;

  cout << "[[openfile ]] Cd : " << dir << endl;
  gDirectory->cd(dir);
  gDirectory->ls();
  TH1F * current_hist = (TH1F*)gDirectory -> Get(histname);
  current_hist -> SetDirectory(0);
  TH1::AddDirectory(kFALSE);

  mapfunc[histname + cyclename + samplename] = current_hist;

  current_file -> Close();
  delete current_file;

}

void make_histogram(TString nameofhistogram, int N_bin, double binx[]){
  
  TString current_data;
  if(nameofhistogram.Contains("EMu") || nameofhistogram.Contains("DiMu")) current_data = SingleMuon;
  else if(nameofhistogram.Contains("DiEle")) current_data = DoubleEG;
  else return;

  TString hstack = nameofhistogram;
  TString func = nameofhistogram;
  TString clone = nameofhistogram;

  maphstack[hstack] = new THStack(hstack, "Stacked_" + nameofhistogram);
  

  int n_kind = 3;
  TString samples_array[] = {WW, ttbar, DY_high};
  Int_t colour_array[] = {419, 416, 400};

  TString name_cycle = nameofhistogram + Cycle_name;
  if(debug) cout << "check1" << endl;
  TString WW = "WWTo2L2Nu_powheg";
  TString WZ_2L = "WZTo2L2Q";
  TString WZ_3L = "WZTo3LNu";
  TString ZZ_2L = "ZZTo2L2Q";
  TString ZZ_4L = "ZZTo4L_powheg";


  GetHist(name_cycle +  WW) -> Add(GetHist(name_cycle + WZ_2L));
  GetHist(name_cycle +  WW) -> Add(GetHist(name_cycle + WZ_3L));
  GetHist(name_cycle +  WW) -> Add(GetHist(name_cycle + ZZ_2L));
  GetHist(name_cycle +  WW) -> Add(GetHist(name_cycle + ZZ_4L));
  
  TString overflow = "overflow";
  Int_t nx_func    = GetHist(nameofhistogram + Cycle_name + current_data) -> GetNbinsX()+1;
  Double_t x1_func = GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinLowEdge(1);
  Double_t bw_func = GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinWidth(nx_func);
  Double_t x2_func = GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinLowEdge(nx_func)+bw_func;
  binx[N_bin] = binx[N_bin - 1] + bw_func;
  if(debug) cout << "[[make_histogram]] binx[N_bin] ; " << binx[N_bin]  << endl;
  mapfunc[func] = new TH1F("", "", nx_func, x1_func, x2_func);
  if(debug) cout << "func rebin rebinning" << endl;
  mapfunc[func + "rebin"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "rebin", binx);
  
  for(int i = 0; i < n_kind; i++){
    if(debug) cout << samples_array[i] << endl;
    if(mapfunc[nameofhistogram + Cycle_name + samples_array[i]]){
      Int_t nx    = GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetNbinsX()+1;
      Double_t x1 = GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetBinLowEdge(1);
      Double_t bw = GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetBinWidth(nx);
      Double_t x2 = GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetBinLowEdge(nx)+bw;

      TH1F *htmp = new TH1F("", "", nx, x1, x2);

      for (Int_t j = 1; j <= nx; j++) {
        htmp -> SetBinContent(j, GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetBinContent(j) );
        htmp -> SetBinError(j, GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetBinError(j) );
      }
      
      if(debug) cout << "rebinning" << endl;
      
      mapfunc[nameofhistogram + Cycle_name + samples_array[i] + "rebin"] = dynamic_cast<TH1F*>(htmp -> Rebin(N_bin, nameofhistogram + Cycle_name + samples_array[i] + "rebin", binx));
      
      GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin") -> SetFillColor(colour_array[i]);
      GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin") -> SetLineColor(colour_array[i]);
      if(debug) cout << samples_array[i] << endl;
      
      // -- divide bin content and error by bin width
      for(int j = 1; j <= N_bin; j++){
	double current_content = GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin") -> GetBinContent(j);
        double current_error = GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin") -> GetBinError(j);
        GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin") -> SetBinContent(j, current_content / (binx[j] - binx[j - 1]) );
        GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin") -> SetBinError(j, current_error / (binx[j] - binx[j - 1]) );
      }
      if(debug) cout << "divide done" << endl;

      maphstack[hstack] -> Add(GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin"));
      mapfunc[func + "rebin"] -> Add(GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin"));
      cout << func + "rebin" << endl;
      cout << "hstack : " << hstack << endl;
    }
  }//for loop 
  
  

}


void signal_rebin(TString nameofhistogram, TString sample_name, int N_bin, double binx[]){
 
  TString hstack = nameofhistogram + sample_name;
  TString func = nameofhistogram + sample_name;
  TString clone = nameofhistogram + sample_name;

 
  TString overflow = "overflow";
  Int_t nx_func    = GetHist(nameofhistogram + Cycle_name + sample_name) -> GetNbinsX()+1;
  Double_t x1_func = GetHist(nameofhistogram + Cycle_name + sample_name) -> GetBinLowEdge(1);
  Double_t bw_func = GetHist(nameofhistogram + Cycle_name + sample_name) -> GetBinWidth(nx_func);
  Double_t x2_func = GetHist(nameofhistogram + Cycle_name + sample_name) -> GetBinLowEdge(nx_func)+bw_func;
  binx[N_bin] = binx[N_bin - 1] + bw_func;
  if(debug) cout << "[[make_histogram]] binx[N_bin] ; " << binx[N_bin]  << endl;
  mapfunc[func] = new TH1F("", "", nx_func, x1_func, x2_func);
  if(debug) cout << "func rebin rebinning" << endl;
  mapfunc[func + "rebin"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "rebin", binx);
 


  Int_t nx    = GetHist(nameofhistogram + Cycle_name + sample_name) -> GetNbinsX()+1;
  Double_t x1 = GetHist(nameofhistogram + Cycle_name + sample_name) -> GetBinLowEdge(1);
  Double_t bw = GetHist(nameofhistogram + Cycle_name + sample_name) -> GetBinWidth(nx);
  Double_t x2 = GetHist(nameofhistogram + Cycle_name + sample_name) -> GetBinLowEdge(nx)+bw;

  TH1F *htmp = new TH1F("", "", nx, x1, x2);

  for (Int_t j = 1; j <= nx; j++) {
    htmp -> SetBinContent(j, GetHist(nameofhistogram + Cycle_name + sample_name) -> GetBinContent(j) );
    htmp -> SetBinError(j, GetHist(nameofhistogram + Cycle_name + sample_name) -> GetBinError(j) );
  }

  if(debug) cout << "rebinning" << endl;

  mapfunc[nameofhistogram + Cycle_name + sample_name + "rebin"] = dynamic_cast<TH1F*>(htmp -> Rebin(N_bin, nameofhistogram + Cycle_name + sample_name+ "rebin", binx));

  if(debug) cout << sample_name << endl;

  // -- divide bin content and error by bin width
  for(int j = 1; j <= N_bin; j++){
    double current_content = GetHist(nameofhistogram + Cycle_name + sample_name + "rebin") -> GetBinContent(j);
    double current_error = GetHist(nameofhistogram + Cycle_name + sample_name + "rebin") -> GetBinError(j);
    GetHist(nameofhistogram + Cycle_name + sample_name + "rebin") -> SetBinContent(j, current_content / (binx[j] - binx[j - 1]) );
    GetHist(nameofhistogram + Cycle_name + sample_name + "rebin") -> SetBinError(j, current_error / (binx[j] - binx[j - 1]) );
  }

  
  mapfunc[func + "rebin"] -> Add(GetHist(nameofhistogram + Cycle_name + sample_name + "rebin"));
  cout << "[[signal_rebin]] func + rebin : " << func + "rebin" << endl;
  

}

void draw_histogram(TString histname){

  if(debug) cout << "[[draw_histogram]] histname : " << histname << endl;

  unsigned int N_bin = map_bin_vector["mZp"].size();
  double current_bins[100];
  for(unsigned int i_bin = 0; i_bin < N_bin; i_bin++){
    current_bins[i_bin] = map_bin_vector["mZp"].at(i_bin);
  }
  
  int N_regions = 1;
  TString regions[] = {"SR"};

  int N_channels = 2;
  TString channels[] = {"DiEle",
                        "DiMu",
  };

  TString current_dir = "empty";
  TString current_channel = "empty";

  for(int i_region = 0; i_region < N_regions; i_region++){
    if(histname.Contains(regions[i_region])) current_dir = regions[i_region];
  }
  for(int i_channel = 0; i_channel < N_channels; i_channel++){
    if(histname.Contains(channels[i_channel])) current_channel = channels[i_channel];
  }

  current_dir = current_dir + "_" + current_channel;

  if(debug) cout << "current_dir : " << current_dir << "/" << histname << endl;
  if(current_dir.Contains("empty") || current_channel.Contains("empty")) return;

  TString WW = "WWTo2L2Nu_powheg";
  TString WZ_2L = "WZTo2L2Q";
  TString WZ_3L = "WZTo3LNu";
  TString ZZ_2L = "ZZTo2L2Q";
  TString ZZ_4L = "ZZTo4L_powheg";
  
  openfile(Cycle_name, DoubleEG, current_dir, histname);
  openfile(Cycle_name, SingleMuon, current_dir, histname);
  openfile(Cycle_name, DY_high, current_dir, histname);
  openfile(Cycle_name, WW, current_dir,  histname);
  openfile(Cycle_name, WZ_2L, current_dir, histname);
  openfile(Cycle_name, WZ_3L, current_dir, histname);
  openfile(Cycle_name, ZZ_2L, current_dir, histname);
  openfile(Cycle_name, ZZ_4L, current_dir, histname);
  openfile(Cycle_name, ttbar, current_dir, histname);
  
  make_histogram(histname, N_bin, current_bins);
  
  TString signal_text_name = "";
  if(histname.Contains("DiMu")) signal_text_name = "HN_pair_MM.txt";
  if(histname.Contains("DiEle")) signal_text_name = "HN_pair_EE.txt";
  
  ifstream data_file;
  data_file.open(signal_text_name);

  char line[500];
  if(data_file.is_open()){
    while(!data_file.eof()){
      data_file.getline(line, 500);
      if(debug) cout << line << endl;
      TString this_line = line;
      if(this_line.Contains("HNPairToJJJJ")){
	openfile(Cycle_name, this_line, current_dir, histname);
	signal_rebin(histname, this_line, N_bin, current_bins);
      }
    }//while not end of file
  }//end if file open

  data_file.close();

}


void get_denom(TString channel){
  
  TString signal_text_name = "";
  if(channel.Contains("DiMu")) signal_text_name = "HN_pair_MM.txt";
  if(channel.Contains("DiEle")) signal_text_name = "HN_pair_EE.txt";

  ifstream data_file;
  data_file.open(signal_text_name);

  char line[500];
  if(data_file.is_open()){
    while(!data_file.eof()){
      data_file.getline(line, 500);
      if(debug) cout << line << endl;
      TString this_line = line;
      if(this_line.Contains("HNPairToJJJJ")){
        openfile(Cycle_name, this_line, "", "sumW");
	mapfunc["denom" + this_line] = new TH1F("", "", 1, 0., 1.);
	
	double current_denom = mapfunc["sumW" + Cycle_name + this_line] -> GetBinContent(2);
	mapfunc["denom" + this_line] -> Fill(0.5, current_denom * 0.005);
	
      }
    }
  }
  data_file.close();
  

}




///////////////////
// -- Main Function
///////////////////
void plot(){
  
  ofstream myfile;
  myfile.open ("masses.txt");
  
  

  for(int i = 0;  i < 26; i++){
    map_bin_vector["mZp"].push_back(0. + 200. * i);
  }
  
  cout << "map_bin_vector[mZp].size() : " << map_bin_vector["mZp"].size() << endl; 
  for(unsigned int i = 0; i < map_bin_vector["mZp"].size(); i++){
    cout << "map_bin_vector[mZp].at(" << i << ") : " << map_bin_vector["mZp"].at(i) << endl;
  }
  
  draw_histogram("mZp_SR_DiMu");
  draw_histogram("mZp_SR_DiEle");
  
  get_denom("DiMu");
  get_denom("DiEle");


  
  
  TFile *MyFile = new TFile("Bkg_VS_signal_cutbased_e.root","RECREATE"); 
  gDirectory -> mkdir("MuMu");
  gDirectory -> cd("MuMu");
  maphstack["mZp_SR_DiMu"] -> SetName("hs_h_OS_lljjjjmass_SR1_DiMu");
  maphstack["mZp_SR_DiMu"] -> Write();

  ifstream mumu_signal_list;
  mumu_signal_list.open("HN_pair_MM.txt");

  char line[500];
  if(mumu_signal_list.is_open()){
    while(!mumu_signal_list.eof()){
      mumu_signal_list.getline(line, 500);
      
      TString this_line = line;
      if(this_line.Contains("HNPairToJJJJ")){
	TObjArray *tx = this_line.Tokenize("_");
	TString mZp_str = ((TObjString *)(tx->At(2)))->String();
	TString mN_str = ((TObjString *)(tx->At(3)))->String();
	
	mZp_str.Remove(0,2);
	mZp_str.Prepend("Zp");
	mN_str.Prepend("H");
	
	myfile << mZp_str << "_" << mN_str << endl;
	
	mapfunc["mZp_SR_DiMu" + this_line + "rebin"] -> SetName("h_OS_lljjjjmass_SR1_DiMu_" + mZp_str + "_" + mN_str);
	mapfunc["mZp_SR_DiMu" + this_line + "rebin"] -> Write();
	mapfunc["denom" + this_line] -> SetName("Den_" + mZp_str + "_" + mN_str + "_h_OS_lljjjjmass_SR1_DiMu");
	mapfunc["denom" + this_line] -> Write();
      }
    }
  }
  mumu_signal_list.close();
  gDirectory ->cd("../");

  
  gDirectory -> mkdir("ElEl");
  gDirectory ->cd("ElEl");
  maphstack["mZp_SR_DiEle"] -> SetName("hs_h_OS_lljjjjmass_SR1_DiEle");
  maphstack["mZp_SR_DiEle"] -> Write();
  
  ifstream elel_signal_list;
  elel_signal_list.open("HN_pair_EE.txt");
  if(elel_signal_list.is_open()){
    while(!elel_signal_list.eof()){
      elel_signal_list.getline(line, 500);

      TString this_line = line;
      if(this_line.Contains("HNPairToJJJJ")){
        TObjArray *tx = this_line.Tokenize("_");
        TString mZp_str = ((TObjString *)(tx->At(2)))->String();
        TString mN_str = ((TObjString *)(tx->At(3)))->String();

	mZp_str.Remove(0,2);
        mZp_str.Prepend("Zp");
	mN_str.Prepend("H");

	mapfunc["mZp_SR_DiEle" + this_line + "rebin"] ->SetName("h_OS_lljjjjmass_SR1_DiEle_" + mZp_str + "_" + mN_str);
	mapfunc["mZp_SR_DiEle" + this_line + "rebin"] ->Write();
	mapfunc["denom" + this_line] -> SetName("Den_" + mZp_str + "_" + mN_str + "_h_OS_lljjjjmass_SR1_DiEle");
	mapfunc["denom" + this_line] ->Write();
	
      }
    }
  }
  
  elel_signal_list.close();
  gDirectory ->cd("../");

  
  
  myfile.close(); 
  
}
