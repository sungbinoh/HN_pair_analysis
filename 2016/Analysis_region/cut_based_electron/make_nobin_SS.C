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
TString SingleMuon = "data_SingleMuon";
TString DoubleEG = "data_DoubleEG";

//DY
TString DY_low = "SKDYJets_10to50";
TString DY_high = "SKDYJets";

TString Wjets = "SKWJets";

//single top
TString SingleTbar_tW = "SKSingleTbar_tW";
TString SingleTbar_t = "SKSingleTbar_t";
TString SingleTop_s = "SKSingleTop_s";
TString SingleTop_tW = "SKSingleTop_tW";
TString SingleTop_t = "SKSingleTop_t";

TString WGtoLNuG = "SKWGtoLNuG";

//VV
TString WW = "SKWW";
TString WZ = "SKWZ";
TString ZZ = "SKZZ";

//VVV
TString WWW = "SKWWW";
TString WWZ = "SKWWZ";
TString WZZ = "SKWZZ";
TString ZZZ = "SKZZZ";

//TString ttbar = "SKTTJets_aMC";
TString ttbar = "SKTT_powheg";
//TString ttbar = "SKTTLL_powheg";
TString ZGto2LG = "SKZGto2LG";

//data driven
TString emu_method = "SKEMu_method";
TString CF_bkg = "SKCF_bkg";

const double alpha = 1 - 0.6827;

bool debug = true;


/// Getting Histogram Function ///////////////////////////////////////////////
TH1F * GetHist(TString hname){
  
  TH1F *h = NULL;
  std::map<TString, TH1F*>::iterator mapit = mapfunc.find(hname);
  if(mapit != mapfunc.end()) return mapit-> second;
  
  return h;
  
}
//////////////////////////////////////////////////////////////////////////////


/// Open ROOT file ///////////////////////////////////////////////////////////
void openfile(TString cyclename, TString samplename, TString dir_1st_tier, TString dir_2nd_tier, TString histname){

  if(debug) cout << "[[openfile]] " << histname + cyclename + samplename << endl;    
  
  TString underbar = "_";
  TString version = "dilep_cat_v8-0-7.root";
  TString filename = cyclename + underbar + samplename + underbar + version;
  if(samplename.Contains("data") || samplename.Contains("NNPDF30_13TeV-powheg")) filename = cyclename + underbar + samplename + "_cat_v8-0-7.root";
  
  TFile *current_file = new TFile ((filename)) ;

  gDirectory->cd(dir_1st_tier);
  gDirectory->cd(dir_1st_tier + dir_2nd_tier);

  //if(debug && samplename.Contains("SingleMuon")) gDirectory -> ls();
  
  TH1F * current_hist = (TH1F*)gDirectory -> Get(histname);
  current_hist -> SetDirectory(0);
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

  TString hstack = "hstack" + nameofhistogram;

  // -- Start from zero : find the bin when MC starts to have non-zero value
  double bin_content = 0.;
  int first_bin = 0;
  while(bin_content < 0.00001 && first_bin * basis_bin * bin_width < x_max ){
    bin_content = cal_content(nameofhistogram, first_bin * basis_bin, (first_bin + 1) * basis_bin);
    first_bin ++;
  }

  // -- Warning empty MC & return
  if(debug) cout << "first bin_content : " << bin_content << ", first_bin : " << first_bin * basis_bin * bin_width << endl;
  if(first_bin * basis_bin * bin_width + 1. > x_max ){
    map_outfile["outfile"] << "############################## EMPTY MC, resign! ##############################" << endl;
    return;
  }
  
  if(debug) cout << "x_max : " << x_max << endl;
  
  map_outfile["outfile"] << "{0";
  for(int i_low = 1; i_low < first_bin; i_low++){
    map_outfile["outfile"] << ", " << i_low * basis_bin * bin_width ;
  }
  
  if(debug) cout << "low end "<< endl;
  
  // -- Run upto xmax to get good binning
  int i_current_x_low = first_bin;
  int i_current_x_high = first_bin + 1;
  int zero_n = 0;
  
  while(i_current_x_high * basis_bin * bin_width < x_max){
    double current_rel_error = cal_rel_stat_error(nameofhistogram, i_current_x_low * basis_bin, i_current_x_high * basis_bin);
    
    zero_n = 0; // -- stop merging bins if zero bin happens more than 4 times
    while(current_rel_error > 0.1 && zero_n < 5 && (i_current_x_high + 1) * basis_bin * bin_width < x_max){
      if(debug) cout << "current_rel_error : " << current_rel_error << ", zero_n : " << zero_n << endl;
      
      i_current_x_high ++;
      current_rel_error = cal_rel_stat_error(nameofhistogram, i_current_x_low * basis_bin, i_current_x_high * basis_bin);
      double content_next_bin = 0.;
      content_next_bin = cal_content(nameofhistogram, (i_current_x_high - 1) * basis_bin, i_current_x_high * basis_bin);
      if( content_next_bin < 0.000001 ) zero_n++;
    }
    
    if(zero_n == 5){
      if(debug) cout << "break" << endl;
      break;
    }
    
    map_outfile["outfile"] << ", " << i_current_x_high * basis_bin * bin_width;
    
    i_current_x_low = i_current_x_high;
    i_current_x_high = i_current_x_high + 1;
  }

  if(debug) cout << "while end" << endl;
  
  if(zero_n > 4){
    while(i_current_x_high * basis_bin * bin_width < x_max){
      map_outfile["outfile"] << ", " << i_current_x_high * basis_bin * bin_width;
      i_current_x_high++;
    }
  }
  
  map_outfile["outfile"] << ", " << i_current_x_high * basis_bin * bin_width << " };" << endl;
  
  
  
}


//////////////////////////////////////////////////////////////////////////////
void makehistogram(TString nameofhistogram, TString channel){
  
  bool blind = false;
  blind = (channel.Contains("SR")) && (!channel.Contains("EMu"));
  
  TString current_data;
  if(channel.Contains("EMu") || channel.Contains("DiMu")) current_data = SingleMuon;
  if(channel.Contains("DiEle")) current_data = DoubleEG;

  TString hstack = "hstack" + nameofhistogram;
  
  maphstack[hstack] = new THStack(hstack, "Stacked_" + nameofhistogram);
  gStyle->SetOptTitle(0);
  
  int n_kind = 3;
  TString samples_array[] = {WZ, emu_method, DY_high};
  Int_t colour_array[] = {419, 416, 400};
  TString samples_legend[] = {"Other backgrounds", "Data Driven Flavour symm. bkg", "Z/#gamma + jets"};
  
  if(debug) cout << "2" << endl;
  
  TString name_cycle = nameofhistogram + Cycle_name;
  if(debug) cout << "check1" << endl;
  GetHist(name_cycle +  DY_high) -> Add(GetHist(name_cycle + DY_low));
  
  if(debug) cout << "2.1" << endl;
  
  GetHist(name_cycle +  SingleTop_s) -> Add(GetHist(name_cycle + SingleTbar_t));
  GetHist(name_cycle +  SingleTop_s) -> Add(GetHist(name_cycle + SingleTop_t));
  GetHist(name_cycle +  WZ) -> Add(GetHist(name_cycle + ZZ));
  GetHist(name_cycle +  WGtoLNuG) -> Add(GetHist(name_cycle + ZGto2LG));
  GetHist(name_cycle +  WWW) -> Add(GetHist(name_cycle + WWZ));
  GetHist(name_cycle +  WWW) -> Add(GetHist(name_cycle + WZZ));
  GetHist(name_cycle +  WWW) -> Add(GetHist(name_cycle + ZZZ));
  GetHist(name_cycle + WZ) -> Add(GetHist(name_cycle +  SingleTop_s));
  GetHist(name_cycle + WZ) -> Add(GetHist(name_cycle +  WWW));
  GetHist(name_cycle + WZ) -> Add(GetHist(name_cycle +  WGtoLNuG));

  Int_t nx_func    = GetHist(name_cycle + current_data) -> GetNbinsX();
  Double_t x1_func = GetHist(name_cycle + current_data) -> GetBinLowEdge(1);
  Double_t bw_func = GetHist(name_cycle + current_data) -> GetBinWidth(nx_func);
  Double_t x2_func = GetHist(name_cycle + current_data) -> GetBinLowEdge(nx_func) + bw_func;
  
  if(debug) cout << "nx_func : " << nx_func << ", x1_func : " << x1_func << ", x2_func : " << x2_func << endl;
  
  int bin_basis = nx_func / 50.;
  if(debug) cout << "bin_basis : " << bin_basis << endl;

  mapfunc[hstack] = new TH1F("", "", nx_func, x1_func, x2_func);
  
  for(int i_stack = 0; i_stack < n_kind; i_stack++){
    mapfunc[hstack] -> Add( GetHist(name_cycle + samples_array[i_stack]) );
  }
  
  find_good_binning(nameofhistogram, bin_basis, bw_func, x1_func, x2_func);

  if(debug) cout << "[[makehistogram]] end" << endl;
  
}
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

void run_for_hist(TString dir_1st_tier, TString histname, TString dir_2nd_tier){

  if(debug) cout << "[[run_for_hist]] " << histname << endl;
  
  openfile(Cycle_name, DoubleEG, dir_1st_tier, dir_2nd_tier, histname);
  openfile(Cycle_name, SingleMuon, dir_1st_tier, dir_2nd_tier, histname);
  openfile(Cycle_name, DY_low, dir_1st_tier, dir_2nd_tier, histname);
  openfile(Cycle_name, DY_high, dir_1st_tier, dir_2nd_tier, histname);
  openfile(Cycle_name, SingleTbar_tW, dir_1st_tier, dir_2nd_tier, histname);
  openfile(Cycle_name, SingleTbar_t, dir_1st_tier, dir_2nd_tier, histname);
  openfile(Cycle_name, SingleTop_s, dir_1st_tier, dir_2nd_tier, histname);
  openfile(Cycle_name, SingleTop_tW, dir_1st_tier, dir_2nd_tier, histname);
  openfile(Cycle_name, SingleTop_t, dir_1st_tier, dir_2nd_tier, histname);
  openfile(Cycle_name, WW, dir_1st_tier, dir_2nd_tier, histname);
  openfile(Cycle_name, WZ, dir_1st_tier, dir_2nd_tier, histname);
  openfile(Cycle_name, ZZ, dir_1st_tier, dir_2nd_tier, histname);
  openfile(Cycle_name, WWW, dir_1st_tier, dir_2nd_tier, histname);
  openfile(Cycle_name, WWZ, dir_1st_tier, dir_2nd_tier, histname);
  openfile(Cycle_name, WZZ, dir_1st_tier, dir_2nd_tier, histname);
  openfile(Cycle_name, ZZZ, dir_1st_tier, dir_2nd_tier, histname);
  openfile(Cycle_name, WGtoLNuG, dir_1st_tier, dir_2nd_tier, histname);
  openfile(Cycle_name, ZGto2LG, dir_1st_tier, dir_2nd_tier, histname);
  openfile(Cycle_name, ttbar, dir_1st_tier, dir_2nd_tier, histname);
  openfile(Cycle_name, CF_bkg, dir_1st_tier, dir_2nd_tier, histname);

  map_outfile["outfile"] << histname << "\t";
  makehistogram(histname, dir_2nd_tier);
  
}


/// Main nction ////////////////////////////////////////////////////////////
void plot(){

  map_outfile["outfile"].open("binning_SS.txt");

  
  TDirectory* origDir = gDirectory;
  
  int n_hist = 11;
  TString histograms[] = {"h_SS_lljjjjmass",
                          "h_SS_lljjjjmass_AK8_0",
                          "h_SS_lljjjjmass_AK8_1",
                          "h_SS_lljjjjmass_AK8_2",
                          "h_SS_leadingljjmass",
                          "h_SS_secondljjmass",
                          "h_SS_llmass",
                          "h_SS_leadingLeptonPt",
                          "h_SS_secondLeptonPt",
                          "h_SS_leadingjet_pt",
                          "h_SS_secondjet_pt"
  };

  int n_dir_1st = 2;
  TString dir_1st_tier[] = {"CR3",
			    "SR1"
  };

  int n_dir_2nd = 1;
  TString dir_2nd_tier[] = {"SS_DiEle",
  };

  int n_channel = 1;
  TString channel[] = {"_SS_DiEle"};
  for(int i_hist = 0; i_hist < n_hist; i_hist++){
    for(int i = 0; i < n_dir_1st; i++){
      for(int j = 0; j < n_channel; j++){//channel
        run_for_hist(dir_1st_tier[i], histograms[i_hist] + "_" + dir_1st_tier[i] + channel[j], channel[j]);
        /*
        origDir -> cd();
        gDirectory -> Cd(dir_1st_tier[i]);
        gDirectory -> Cd(dir_1st_tier[i] + channel[j]);
        mapfunc[histograms[i_hist] + "_" + dir_1st_tier[i] + "_OS_EMu"+ "emu"] -> SetName(histograms[i_hist] + "_" + dir_1st_tier[i] + channel[j]);
        mapfunc[histograms[i_hist] + "_" + dir_1st_tier[i] + "_OS_EMu"+ "emu"] -> Write();
	*/
        mapfunc.clear();
	maphist.clear();
      }
    }
  }

  map_outfile["outfile"].close();

}// End of Main Function ////////////////////////////////////////////////////// 
