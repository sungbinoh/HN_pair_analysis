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


//emu ratio results
double emu_over_ee = 0.450;
double emu_over_mumu = 0.647;

const double alpha = 1 - 0.6827;
bool debug = false;

/// Getting Histogram Function ///////////////////////////////////////////////
TH1F * GetHist(TString hname){
  
  TH1F *h = NULL;
  std::map<TString, TH1F*>::iterator mapit = mapfunc.find(hname);
  if(mapit != mapfunc.end()) return mapit-> second;
  
  return h;
  
}
//////////////////////////////////////////////////////////////////////////////


/// Open ROOT file ///////////////////////////////////////////////////////////
void openfile(TString cyclename, TString samplename, TString dir_1st_tier, TString histname){
  
  TString underbar = "_";
  TString version = "dilep_cat_v8-0-7.root";
  TString filename = cyclename + underbar + samplename + underbar + version;
  if(samplename.Contains("data") || samplename.Contains("NNPDF30_13TeV-powheg")) filename = cyclename + underbar + samplename + "_cat_v8-0-7.root";
  
  //cout << "opening : " << filename << endl;
  
  TFile *current_file = new TFile ((filename));

  // -- Get EMu CF hist
  gDirectory->cd(dir_1st_tier);
  gDirectory->cd(dir_1st_tier + "_SS_EMu_CF");
  
  if(debug && samplename.Contains("SingleMuon")) gDirectory -> ls();
  
  TH1F * current_hist_EMu = (TH1F*)gDirectory -> Get(histname + "_" + dir_1st_tier + "_SS_EMu_CF");
  current_hist_EMu -> SetDirectory(0);
  TH1::AddDirectory(kFALSE);
  
  mapfunc[histname + "_" + dir_1st_tier + "_EMu_CF" + cyclename + samplename] = current_hist_EMu;


  // -- Get DiEle CF hist
  gDirectory -> cd("../");
  gDirectory->cd(dir_1st_tier + "_SS_DiEle_CF");
  
  if(debug && samplename.Contains("SingleMuon")) gDirectory -> ls();

  TH1F * current_hist_DiEle = (TH1F*)gDirectory -> Get(histname + "_" + dir_1st_tier + "_SS_DiEle_CF");
  current_hist_DiEle -> SetDirectory(0);
  TH1::AddDirectory(kFALSE);

  mapfunc[histname + "_" + dir_1st_tier + "_DiEle_CF" + cyclename + samplename] = current_hist_DiEle;
    
  //if(debug) current_hist -> Draw();
  //return;
  //cout << "[[openfile]] " << histname + cyclename + samplename << endl;
  
  //if(debug) mapfunc[histname + cyclename + samplename] -> Draw();
  //return;
  current_file -> Close();
  delete current_file;

  //if(debug) mapfunc[histname + cyclename + samplename] -> Draw();
  //return;
}
//////////////////////////////////////////////////////////////////////////////


void make_hist_after_emu_method_and_CF(TString nameofhistogram, TString dir_1st_tier, TString channel){
  if(debug) cout << "[[make_hist_after_emu_method_and_CF]] making " << nameofhistogram << "_" << channel << endl;
  
  double current_emu_ratio = 1.;
  if(channel.Contains("DiEle")) current_emu_ratio = emu_over_ee;
  else if(channel.Contains("DiMu")) current_emu_ratio = emu_over_mumu;
  
  TString current_data = SingleMuon;
  TString cyclename = Cycle_name;
  
  if(debug) cout << "[[make_hist_after_emu_method_and_CF]] making data driven plot of " + nameofhistogram + "_" + channel << endl;
  
  TString histname_EMu = nameofhistogram + "_" + dir_1st_tier + "_EMu_CF";
  TString histname_DiEle = nameofhistogram + "_" + dir_1st_tier + "_DiEle_CF";

  
  //clone emu data plot
  //mapfunc[histname_EMu + cyclename + current_data] -> Draw();
  //return;
  if(debug) cout << "returned?" << endl;
  mapfunc[nameofhistogram + "CF"] = (TH1F*)GetHist(histname_EMu + cyclename + current_data) -> Clone(nameofhistogram);
  //mapfunc[nameofhistogram + "CF"] -> Draw();
  //return;
  
  if(debug) cout << "cloned emu data hist" << endl;

  TString name_cycle = histname_EMu + cyclename;
  //subtract non flavour symmetric bkgs
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  WZ), -1);
  if(debug) cout << "minus WZ" << endl;
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  ZZ), -1);
  if(debug) cout << "minus ZZ" << endl;
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  SingleTop_s), -1);
  if(debug) cout << "minus SingleTop_s" << endl;
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  SingleTbar_t), -1);
  if(debug) cout << "minus SingleTbar_t" << endl;
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  SingleTop_t), -1);
  if(debug) cout << "minus SingleTop_t" << endl;
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  WGtoLNuG), -1);
  if(debug) cout << "minus WGtoLNuG" << endl;
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  ZGto2LG), -1);
  if(debug) cout << "minus ZGto2LG" << endl;
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  WWW), -1);
  if(debug) cout << "minus WWW" << endl;
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  WWZ), -1);
  if(debug) cout << "minus WWZ" << endl;
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  WZZ), -1);
  if(debug) cout << "minus WZZ" << endl;
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  ZZZ), -1);
  if(debug) cout << "minus ZZZ" << endl;
  mapfunc[nameofhistogram + "CF"] -> Scale(current_emu_ratio);
  if(debug) cout << "scaled" << endl;
 
  name_cycle = histname_DiEle + cyclename;;
  
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle + DY_low));
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle + DY_high));
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  WZ));
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  ZZ));
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  SingleTop_s));
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  SingleTbar_t));
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  SingleTop_t));
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  WGtoLNuG));
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  ZGto2LG));
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  WWW));
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  WWZ));
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  WZZ));
  mapfunc[nameofhistogram + "CF"] -> Add(GetHist(name_cycle +  ZZZ));

  if(debug) cout << "CF added" << endl;
  if(debug) cout << nameofhistogram + "emu" + " finished" << endl;

  
}

//////////////////////////////////////////////////////////////////////////////

void save_hist_after_emu_method(TString dir_1st_tier, TString histname, TString channel){
  
  if(debug) cout << "open files for " << dir_1st_tier + "/" + histname << endl;
  
  openfile(Cycle_name, SingleMuon, dir_1st_tier, histname);
  openfile(Cycle_name, DY_low, dir_1st_tier, histname);
  openfile(Cycle_name, DY_high, dir_1st_tier, histname);
  openfile(Cycle_name, SingleTbar_tW, dir_1st_tier, histname);
  openfile(Cycle_name, SingleTbar_t, dir_1st_tier, histname);
  openfile(Cycle_name, SingleTop_s, dir_1st_tier, histname);
  openfile(Cycle_name, SingleTop_tW, dir_1st_tier, histname);
  openfile(Cycle_name, SingleTop_t, dir_1st_tier, histname);
  openfile(Cycle_name, WW, dir_1st_tier, histname);
  openfile(Cycle_name, WZ, dir_1st_tier, histname);
  openfile(Cycle_name, ZZ, dir_1st_tier, histname);
  openfile(Cycle_name, WWW, dir_1st_tier, histname);
  openfile(Cycle_name, WWZ, dir_1st_tier, histname);
  openfile(Cycle_name, WZZ, dir_1st_tier, histname);
  openfile(Cycle_name, ZZZ, dir_1st_tier, histname);
  openfile(Cycle_name, WGtoLNuG, dir_1st_tier, histname);
  openfile(Cycle_name, ZGto2LG, dir_1st_tier, histname);
  openfile(Cycle_name, ttbar, dir_1st_tier, histname);
  
  if(debug) cout << "[[save_hist_after_emu_method]] making " << histname << " for " << channel << endl;
  make_hist_after_emu_method_and_CF(histname, dir_1st_tier, channel);
  
}



/// Main nction ////////////////////////////////////////////////////////////
void plot(){
  
  //TDirectory*  = gDirectory;
  cout << "defining arrays" << endl;
  
  int n_hist = 11;
  TString histograms[] = {"lljjjjmass",
			  "lljjjjmass_AK8_0",
			  "lljjjjmass_AK8_1",
			  "lljjjjmass_AK8_2",
			  "leadingljjmass",
			  "secondljjmass",
			  "llmass",
			  "leadingLeptonPt",
			  "secondLeptonPt",
			  "leadingjet_pt",
			  "secondjet_pt"
  };

  int n_dir_1st = 5;
  TString dir_1st_tier[] = {"CR1",
			    "CR2",
			    "CR3",
			    "CR4",
			    "SR1"
  };

  int n_dir_2nd = 3;
  TString dir_2nd_tier[] = {"SS_DiEle",
			    "SS_DiMu",
			    "SS_EMu"
  };
  
  cout << "making output root file and set directories" << endl;
  
  TFile *MyFile = new TFile("HN_pair_all_SKCF_bkg_dilep_cat_v8-0-7.root","RECREATE");
  TDirectory* origDir = gDirectory;
  for(int i = 0; i < n_dir_1st; i++){
    gDirectory -> mkdir( dir_1st_tier[i] );
    gDirectory -> Cd( dir_1st_tier[i] );
    for(int j = 0; j < n_dir_2nd; j++){
      gDirectory -> mkdir( dir_1st_tier[i] + "_" + dir_2nd_tier[j] );
    }
    gDirectory -> Cd("../");
  }

  cout << "making emu method hists" << endl;
  
  TString channel[] = {"_SS_DiEle"};
  for(int i_hist = 0; i_hist < n_hist; i_hist++){
    for(int i = 0; i < n_dir_1st; i++){
      for(int j = 0; j < 1; j++){  
	save_hist_after_emu_method(dir_1st_tier[i], "h_OS_" +  histograms[i_hist], channel[j]);
	//return;
	origDir -> cd();
	gDirectory -> Cd(dir_1st_tier[i]);
	gDirectory -> Cd(dir_1st_tier[i] + "_SS_DiEle");
	mapfunc["h_OS_" + histograms[i_hist] + "CF"] -> SetName("h_SS_" + histograms[i_hist] + "_" + dir_1st_tier[i] + channel[j]);
	mapfunc["h_OS_" + histograms[i_hist] + "CF"] -> Write();
	
	//mapfunc.clear();
	//maphist.clear();
      }
    }
  }

  
  

  

}// End of Main Function ////////////////////////////////////////////////////// 
