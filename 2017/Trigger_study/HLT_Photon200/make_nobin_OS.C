/*
Code to call root files containing histograms and to stack them
-----------------------------------------------------------------------------------------------
openfile : open root file with certain name and call histograms in it
gethist : call histogram saved in maphist
makehistogram : make stacked histogrtam with applied xmin, xmax, rebin, ymax of nameofhistogram
plot() : main working function
*/

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
map<TString, ofstream> map_outfile;
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

//cycle name
TString Cycle_name = "HN_pair_all";
//sample names
//TString data = "data_DoubleMuon";
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

double rel_stat_error_cutoff = 0.05;
double basis_bin_fraction = 0.05;
int zero_n_cutoff = 10;


/// Getting Histogram Function ///////////////////////////////////////////////
TH1F * GetHist(TString hname){
  
  TH1F *h = NULL;
  std::map<TString, TH1F*>::iterator mapit = mapfunc.find(hname);
  if(mapit != mapfunc.end()) return mapit-> second;
  
  return h;
  
}
//////////////////////////////////////////////////////////////////////////////


/// Open ROOT file ///////////////////////////////////////////////////////////
void openfile(TString cyclename, TString samplename, TString dir,TString histname){

  TString filename = cyclename + "_" + samplename + ".root";
  cout << "[[openfile ]]opening : " << filename << endl;
  
  TFile *current_file = new TFile ((filename)) ;

  cout << "[[openfile ]] Cd : " << dir << endl;
  TDirectory *current_dir = (TDirectory *)current_file -> Get(dir); 
  if(!current_dir) return;
  gDirectory->cd(dir);
  gDirectory -> ls();

  cout << "[[openfile ]] Get : " << histname << endl;

  //if(debug && samplename.Contains("SingleMuon")) gDirectory -> ls();
  
  TH1F * current_hist = (TH1F*)gDirectory -> Get(histname);
  if(current_hist){
    current_hist -> SetDirectory(0);
  }
  TH1::AddDirectory(kFALSE);
  
  mapfunc[histname + cyclename + samplename] = current_hist;
  
  current_file -> Close();
  delete current_file;
  
  
}
//////////////////////////////////////////////////////////////////////////////


/// -- Get Info from Hist
double cal_rel_stat_error(TString nameofhistogram, int i_start, int i_end){
  
  TString hstack = "hstack" + nameofhistogram;
  double sum_error_sq = 0.;
  double sum_content = 0.;
  for(int i_stat = i_start; i_stat < i_end; i_stat++){
    double current_content = mapfunc[hstack] -> GetBinContent(i_stat + 1);
    double current_error = mapfunc[hstack] -> GetBinError(i_stat + 1);

    sum_content = sum_content + current_content * current_content;
    sum_error_sq = sum_error_sq + current_error * current_error;
  }

  double rel_stat_error = pow(sum_error_sq, 0.5) / sum_content;

  return rel_stat_error;
}

double cal_content(TString nameofhistogram, int i_start, int i_end){
  
  TString hstack = "hstack" + nameofhistogram;
  double sum_content = 0.;
  for(int i_stat = i_start; i_stat < i_end; i_stat++){
    double current_content = mapfunc[hstack] -> GetBinContent(i_stat + 1);
    sum_content = sum_content +	current_content * current_content;
  }

  return sum_content;
}
//////////////////////////////////////////////////////////////////////////////


void find_good_binning(TString nameofhistogram, int basis_bin, double bin_width, double x_min, double x_max){

  //if(nameofhistogram.Contains("h_OS_lljjjjmass_AK8_1_CR3_OS_DiEle")) debug = true;
  
  TString hstack = "hstack" + nameofhistogram;

  // -- Start from zero : find the bin when MC starts to have non-zero value
  double bin_content = 0.;
  int first_bin = 0;
  while(bin_content < 0.00001 && first_bin * basis_bin * bin_width + x_min < x_max ){
    bin_content = cal_content(nameofhistogram, first_bin * basis_bin, (first_bin + 1) * basis_bin);
    first_bin ++;
  }

  // -- Warning empty MC & return
  if(debug) cout << "first bin_content : " << bin_content << ", first_bin : " << first_bin * basis_bin * bin_width + x_min << endl;
  if(first_bin * basis_bin * bin_width + 1. > x_max ){
    map_outfile["outfile"] << "############################## EMPTY MC, resign! ##############################" << endl;
    return;
  }
  
  if(debug) cout << "x_max : " << x_max << endl;
  
  if(nameofhistogram.Contains("Eta")) map_outfile["outfile"] << "{-3.0";
  else map_outfile["outfile"] << "{0.";

  for(int i_low = 1; i_low < first_bin; i_low++){
    map_outfile["outfile"] << ", " << i_low * basis_bin * bin_width + x_min ;
  }
  
  if(debug) cout << "low end "<< endl;
  
  // -- Run upto xmax to get good binning
  int i_current_x_low = first_bin;
  int i_current_x_high = first_bin + 1;
  int zero_n = 0;

  while(i_current_x_high * basis_bin * bin_width + x_min < x_max){
    double current_rel_error = fabs(cal_rel_stat_error(nameofhistogram, i_current_x_low * basis_bin , i_current_x_high * basis_bin));
    double next_bin_value = 0.;
    zero_n = 0; // -- stop merging bins if zero bin happens more than 4 times
    while( ((current_rel_error > rel_stat_error_cutoff && zero_n < zero_n_cutoff) || next_bin_value < 0 ) && (i_current_x_high + 1) * basis_bin * bin_width + x_min < x_max){
      if(debug) cout << "current_rel_error : " << current_rel_error << ", zero_n : " << zero_n << endl;
      
      i_current_x_high ++;
      current_rel_error = fabs(cal_rel_stat_error(nameofhistogram, i_current_x_low * basis_bin, i_current_x_high * basis_bin));
      double content_next_bin = 0.;
      content_next_bin = cal_content(nameofhistogram, (i_current_x_high - 1) * basis_bin, i_current_x_high * basis_bin);
      if( fabs(content_next_bin) < 0.000001 ) zero_n++;
      next_bin_value = content_next_bin;
    }
    
    if(zero_n == zero_n_cutoff){
      if(debug) cout << "break" << endl;
      break;
    }
    
    map_outfile["outfile"] << ", " << i_current_x_high * basis_bin * bin_width + x_min;
    
    i_current_x_low = i_current_x_high;
    i_current_x_high = i_current_x_high + 1;
  }

  if(debug) cout << "while end" << endl;
  
  if(zero_n > zero_n_cutoff - 1){
    while(i_current_x_high * basis_bin * bin_width + x_min < x_max){
      map_outfile["outfile"] << ", " << i_current_x_high * basis_bin * bin_width + x_min;
      i_current_x_high++;
    }
  }
  
  map_outfile["outfile"] << ", " << i_current_x_high * basis_bin * bin_width + x_min << " };" << endl;
  
  debug = false;
  
}


//////////////////////////////////////////////////////////////////////////////
void makehistogram(TString nameofhistogram, TString channel){
  
  bool blind = false;
  blind = (channel.Contains("SR")) && (!channel.Contains("EMu"));
  
  TString current_data;
  if(channel.Contains("EMu") || channel.Contains("DiMu")) current_data = SingleMuon;
  if(channel.Contains("DiEle")) current_data = DoubleEG;

  TString hstack = "hstack" + nameofhistogram;
  
  gStyle->SetOptTitle(0);
  
  int n_kind = 3;
  TString samples_array[] = {WW, ttbar, DY_high};
  bool hist_true[] = {false, false, false};
  
  if(debug) cout << "2" << endl;
  
  TString name_cycle = nameofhistogram + Cycle_name;
  if(debug) cout << "check1" << endl;
    
  if(debug) cout << "2.1" << endl;
  
  if(!GetHist(name_cycle + current_data)){
    map_outfile["outfile"]  << "############################## EMPTY DATA, resign! ##############################" << endl;
    return;
  }
  
  Int_t nx_func    = GetHist(name_cycle + current_data) -> GetNbinsX();
  Double_t x1_func = GetHist(name_cycle + current_data) -> GetBinLowEdge(1);
  Double_t bw_func = GetHist(name_cycle + current_data) -> GetBinWidth(nx_func);
  Double_t x2_func = GetHist(name_cycle + current_data) -> GetBinLowEdge(nx_func) + bw_func;
   
  for(int i_sample = 0; i_sample < n_kind; i_sample++){
    mapfunc[samples_array[i_sample] + hstack] = new TH1F("", "", nx_func, x1_func, x2_func);
  }
  
  TString WWs[] ={WW, WZ_2L, WZ_3L, ZZ_2L, ZZ_4L};
  int N_WWs = 5;
  for(int i_WW = 0; i_WW < N_WWs; i_WW++){
    if(GetHist(name_cycle + WWs[i_WW])){
      mapfunc[WW + hstack] -> Add(GetHist(name_cycle + WWs[i_WW]));
      hist_true[0] = true;
    }
  }

  if(GetHist(name_cycle + ttbar)){
    mapfunc[ttbar + hstack] -> Add(GetHist(name_cycle + ttbar));
    hist_true[1] = true;
  }
  if(GetHist(name_cycle + DY_high)){
    mapfunc[DY_high + hstack] -> Add(GetHist(name_cycle + DY_high));
    hist_true[2] = true;
  }
  
  if(debug) cout << "nx_func : " << nx_func << ", x1_func : " << x1_func << ", x2_func : " << x2_func << endl;
  
  int bin_basis = nx_func * basis_bin_fraction;
  if(debug) cout << "bin_basis : " << bin_basis << endl;

  mapfunc[hstack] = new TH1F("", "", nx_func, x1_func, x2_func);
  
  for(int i_stack = 0; i_stack < n_kind; i_stack++){
    if(hist_true[i_stack]) mapfunc[hstack] -> Add( GetHist(name_cycle + samples_array[i_stack]) );
  }
  
  find_good_binning(nameofhistogram, bin_basis, bw_func, x1_func, x2_func);
  
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

void run_for_hist(TString dir, TString histname, TString dir_2nd_tier){

  openfile(Cycle_name, DoubleEG, dir, histname);
  openfile(Cycle_name, SingleMuon, dir, histname);
  openfile(Cycle_name, DY_high, dir, histname);
  openfile(Cycle_name, WW, dir, histname);
  openfile(Cycle_name, WZ_2L, dir, histname);
  openfile(Cycle_name, WZ_3L, dir, histname);
  openfile(Cycle_name, ZZ_2L, dir, histname);
  openfile(Cycle_name, ZZ_4L, dir, histname);
  openfile(Cycle_name, ttbar, dir, histname);
  
  cout << "[[run_for_hist]] File Open Finished" << endl;
  
  map_outfile["outfile"] << histname << "\t";
  makehistogram(histname, dir_2nd_tier);
  
}


/// Main nction ////////////////////////////////////////////////////////////
void plot(){

  map_outfile["outfile"].open("binning_OS.txt");
    
  TDirectory* origDir = gDirectory;
  
  int n_hist = 17;
  TString histograms[] = {"mZp",
			  "mZp_0AK8",
			  "mZp_1AK8",
			  "mZp_2AK8",
			  "mN",
			  "mll",
			  //"mll_tight",
			  "Nvtx_69p2",
			  "Nvtx_noCorr",
			  "Lepton_0_Eta",
                          "Lepton_0_Pt",
			  "Lepton_0_dXY",
			  "Lepton_0_TrkRelIso",
			  "Lepton_1_Eta",
                          "Lepton_1_Pt",
			  "Lepton_1_dXY",
			  "Lepton_1_TrkRelIso",
			  "MET"
  };

  int N_regions = 5;
  TString regions[] = {"CR_Zmass",
                       "CR_ttbar",
		       "tight_CR_Zmass",
		       "tight_CR_ttbar",//"CR_inv_mll",
		       "SR"
  };
  
  int N_channels = 3;
  TString channels[] = {"DiEle_central",
			"DiMu_central",
			"EMu_central"
  };
  
  
  for(int i_hist = 0; i_hist < n_hist; i_hist++){
    for(int i = 0; i < N_regions; i++){
      for(int j = 0; j < N_channels; j++){//channel
	TString dir = regions[i] + "_" + channels[j];
	
	if(histograms[i_hist].Contains("MET")){
	  if(regions[i].Contains("CR_Zmass") || regions[i].Contains("CR_ttbar")){
	    run_for_hist(dir, histograms[i_hist] + "_" + dir, channels[j]);
	  }
	  else continue;;
	}
	

	if(regions[i].Contains("SR")){
	  if(channels[j].Contains("DiEle") || channels[j].Contains("DiMu")){
	    if(histograms[i_hist].Contains("mZp") || histograms[i_hist].Contains("mN") || histograms[i_hist].Contains("Pt")) run_for_hist(dir, histograms[i_hist] + "_" + dir, channels[j]);
	  }
	}
	else if(histograms[i_hist].Contains("mN")){
	  if(!regions[i].Contains("SR")) continue;
	}
	else if(regions[i].Contains("tight")){
	  if(histograms[i_hist].Contains("Lepton")) run_for_hist(dir, histograms[i_hist] + "_" + dir, channels[j]);
	} 
	else run_for_hist(dir, histograms[i_hist] + "_" + dir, channels[j]);
	
	mapfunc.clear();
	maphist.clear();
      }
    }
  }

  map_outfile["outfile"].close();

}// End of Main Function ////////////////////////////////////////////////////// 
