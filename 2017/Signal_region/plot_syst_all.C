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
#include "setTDRStyle.C"

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
map<TString, std::vector<double>> map_bin_vector;

//cycle name
TString Cycle_name = "HN_pair_all";
//sample names
//TString data = "data_DoubleMuon";
TString SingleMuon = "SingleMuon";
TString DoubleEG = "DoubleEG";

//DY
TString DY_high = "DYJets";


//Wjets
TString WJets = "WJets_MG";


//VV
TString VV = "VV";
TString WW = "WWTo2L2Nu_powheg";
TString WZ_2L = "WZTo2L2Q";
TString WZ_3L = "WZTo3LNu";
TString ZZ_2L = "ZZTo2L2Q";
TString ZZ_4L = "ZZTo4L_powheg";

TString ttbar = "TTLL_powheg";

const double alpha = 1 - 0.6827;

bool debug = false;

// -- Set TDR Style
void rootlogon()
{
  // Load CMS style
  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();
}

/// Getting Histogram Function ///////////////////////////////////////////////
TH1F * GetHist(TString hname){
  
  TH1F *h = NULL;
  std::map<TString, TH1F*>::iterator mapit = mapfunc.find(hname);
  if(mapit != mapfunc.end()) return mapit-> second;
  
  return h;
  
}
//////////////////////////////////////////////////////////////////////////////


/// Open ROOT file ///////////////////////////////////////////////////////////
void openfile(TString cyclename, TString samplename, TString dir, TString histname){
  
  TString filename = cyclename + "_" + samplename + ".root";
  //cout << "[[openfile ]]opening : " << filename << endl;

  TFile *current_file = new TFile ((filename)) ;
  
  //cout << "[[openfile ]] Cd : " << dir << endl;
  gDirectory->cd(dir);
      
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


// -- Maing drawing function
void make_histogram(TString nameofhistogram, TString current_histname, int N_bin, double binx[]){

  if(debug){
    cout << "[[make_histogram]] checking binning arrary" << endl;
    for(int i = 0; i < N_bin; i++){
      cout << binx[i] << ", ";
    }
    cout << "" << endl;
  }

  const int N_syst_x = N_bin - 1;
  double syst_x_bin[N_syst_x];
  double syst_x_bin_error[N_syst_x];
  double syst_Stat[N_syst_x];
  double syst_ElectronScale[N_syst_x];
  double syst_ElectronSmear[N_syst_x];
  double syst_JetsRes[N_syst_x];
  double syst_JetsScale[N_syst_x];
  double syst_SD_JMR[N_syst_x];
  double syst_SD_JMS[N_syst_x];
  double syst_PUReweight[N_syst_x];
    
  for(int i = 0; i < N_syst_x; i++){
    syst_x_bin[i] = (binx[i+1] + binx[i]) / 2.0;
    syst_x_bin_error[i] = (binx[i+1] - binx[i]) / 2.0;
  }
  for(int i = 0; i < N_syst_x; i++){
    cout << syst_x_bin[i] << ", ";
  }
  cout << "" << endl;
  for(int i = 0; i < N_syst_x; i++){
    cout << syst_x_bin_error[i] << ", ";
  }
  cout << "" << endl;
  
  
  const int N_syst = 15;
  TString systematics[N_syst] = {"central",
                                 "ElectronScaleDown", "ElectronScaleUp",
                                 "ElectronSmearDown", "ElectronSmearUp",
                                 "JetsResDown", "JetsResUp",
                                 "JetsScaleDown", "JetsScaleUp",
				 "SD_JMR_Down", "SD_JMR_Up",
                                 "SD_JMS_Down", "SD_JMS_Up",
				 //"PDFNormDown", "PDFNormUp"
				 "PUReweight_Down", "PUReweight_Up",
  };

  
  
  TString current_data;
  if(nameofhistogram.Contains("EMu") || nameofhistogram.Contains("DiMu")) current_data = SingleMuon;
  else if(nameofhistogram.Contains("DiEle")) current_data = DoubleEG;
  else return;
  
  TString title_y = "Rel. Syst";;

  
  bool blind = false;
  blind = (nameofhistogram.Contains("SR")) && (!nameofhistogram.Contains("EMu"));
  
  TString pad1 = nameofhistogram;
  TString pad2 = nameofhistogram;
  TString canvas = nameofhistogram;
  TString hstack = nameofhistogram;
  TString legend = nameofhistogram;
  TString func = nameofhistogram;
  TString clone = nameofhistogram;
  TString line = nameofhistogram;
  pad1.Insert(0, "pad1_");
  pad2.Insert(0, "pad2_");
  canvas.Insert(0, "c_");
  hstack.Insert(0, "hs_");
  legend.Insert(0, "legend_");
  func.Insert(0, "ratio_");
  clone.Insert(0, "h3_");
  line.Insert(0, "l_");
  
  mapcanvas[canvas] = new TCanvas(canvas,"",800,800);
  gStyle -> SetOptStat(1111);
  mapcanvas[canvas] -> SetTopMargin( 0.05 );
  mapcanvas[canvas] -> SetBottomMargin( 0.13 );
  mapcanvas[canvas] -> SetRightMargin( 0.05 );
  mapcanvas[canvas] -> SetLeftMargin( 0.16 );
  //mapcanvas[canvas] -> SetLogy();
  
  maplegend[legend] = new TLegend(0.20, 0.60, 0.50, 0.92);
  
  if(debug) cout << "1" << endl;
  
  int n_kind = 4;
  TString samples_array[] = {VV, WJets, ttbar, DY_high};
  Int_t colour_array[] = {419, 600, 416, 400};
  TString samples_legend[] = {"Other backgrounds", "WJets", "ttbar", "Z/#gamma + jets"};
  
  if(debug) cout << "2" << endl;
  
  TString name_cycle = nameofhistogram + Cycle_name;
  if(debug) cout << "check1" << endl;
    
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
    
  mapfunc[func + "syst_ElectronScale"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "rebin", binx);
  mapfunc[func + "syst_ElectronSmear"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "rebin", binx);
  mapfunc[func + "syst_JetsRes"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "rebin", binx);
  mapfunc[func + "syst_JetsScale"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "rebin", binx);
  mapfunc[func + "syst_SD_JMR"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "rebin", binx);
  mapfunc[func + "syst_SD_JMS"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "rebin", binx);
  mapfunc[func + "syst_PUReweight"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "rebin", binx);

  
  for(int i = 0; i < n_kind; i++){
    if(debug) cout << samples_array[i] << endl;
    if(mapfunc[nameofhistogram + Cycle_name + samples_array[i]]){
      Int_t nx    = GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetNbinsX()+1;
      Double_t x1 = GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetBinLowEdge(1);
      Double_t bw = GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetBinWidth(nx);
      Double_t x2 = GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetBinLowEdge(nx)+bw;
      
      TH1F *htmp = new TH1F("", "", nx, x1, x2);
      TH1F *htmp_ElectronScale = new TH1F("", "", nx, x1, x2);
      TH1F *htmp_ElectronSmear = new TH1F("", "", nx, x1, x2);
      TH1F *htmp_JetsRes = new TH1F("", "", nx, x1, x2);
      TH1F *htmp_JetsScale = new TH1F("", "", nx, x1, x2);
      TH1F *htmp_SD_JMR = new TH1F("", "", nx, x1, x2);
      TH1F *htmp_SD_JMS = new TH1F("", "", nx, x1, x2);
      TH1F *htmp_PUReweight = new TH1F("", "", nx, x1, x2);


      
      for (Int_t j = 1; j <= nx; j++) {
	double current_bin_content = GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetBinContent(j);
        htmp -> SetBinContent(j, current_bin_content);
	htmp_ElectronScale -> SetBinContent(j, current_bin_content);
        htmp_ElectronSmear -> SetBinContent(j, current_bin_content);
        htmp_JetsRes -> SetBinContent(j, current_bin_content);
        htmp_JetsScale -> SetBinContent(j, current_bin_content);
	htmp_SD_JMR -> SetBinContent(j, current_bin_content);
        htmp_SD_JMS -> SetBinContent(j, current_bin_content);
	htmp_PUReweight -> SetBinContent(j, current_bin_content);

	double stat_error, ElectronScale_error, ElectronSmear_error, JetsRes_error, JetsScale_error, SD_JMR_error, SD_JMS_error, PUReweight_error;
	
	stat_error = GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetBinError(j); //stat error
	if(GetHist(current_histname + "_ElectronScaleDown"   + Cycle_name + samples_array[i]) || GetHist(current_histname + "_ElectronScaleUp" + Cycle_name + samples_array[i]) ){
	  ElectronScale_error = 0.;
	  if(GetHist(current_histname + "_ElectronScaleDown" + Cycle_name + samples_array[i])) ElectronScale_error += pow(current_bin_content - GetHist(current_histname + "_ElectronScaleDown" + Cycle_name + samples_array[i]) -> GetBinContent(j), 2);
	  if(GetHist(current_histname + "_ElectronScaleUp"   + Cycle_name + samples_array[i])) ElectronScale_error += pow(current_bin_content - GetHist(current_histname + "_ElectronScaleUp" + Cycle_name + samples_array[i]) -> GetBinContent(j), 2);
	  ElectronScale_error = pow(ElectronScale_error, 0.5);
	}
	if(GetHist(current_histname + "_ElectronSmearDown"   + Cycle_name + samples_array[i]) || GetHist(current_histname + "_ElectronSmearUp" + Cycle_name + samples_array[i]) ){
	  ElectronSmear_error = 0.;
	  if(GetHist(current_histname + "_ElectronSmearDown" + Cycle_name + samples_array[i])) ElectronSmear_error += pow(current_bin_content - GetHist(current_histname + "_ElectronSmearDown" + Cycle_name + samples_array[i]) -> GetBinContent(j), 2);
	  if(GetHist(current_histname + "_ElectronSmearUp"   + Cycle_name + samples_array[i])) ElectronSmear_error += pow(current_bin_content - GetHist(current_histname + "_ElectronSmearUp" + Cycle_name + samples_array[i]) -> GetBinContent(j), 2);
	  ElectronSmear_error = pow(ElectronSmear_error, 0.5);
	}
	if(GetHist(current_histname + "_JetsResDown"         + Cycle_name + samples_array[i]) || GetHist(current_histname + "_JetsResUp"       + Cycle_name + samples_array[i]) ){
	  JetsRes_error = 0.;
	  if(GetHist(current_histname + "_JetsResDown"       + Cycle_name + samples_array[i])) JetsRes_error += pow(current_bin_content - GetHist(current_histname + "_JetsResDown"       + Cycle_name + samples_array[i]) -> GetBinContent(j), 2);
	  if(GetHist(current_histname + "_JetsResUp"         + Cycle_name + samples_array[i])) JetsRes_error += pow(current_bin_content - GetHist(current_histname + "_JetsResUp"       + Cycle_name + samples_array[i]) -> GetBinContent(j), 2);
	  JetsRes_error = pow(JetsRes_error, 0.5);
	}
	if(GetHist(current_histname + "_JetsScaleDown"       + Cycle_name + samples_array[i]) || GetHist(current_histname + "_JetsScaleUp"     + Cycle_name + samples_array[i]) ){
	  JetsScale_error = 0.;
	  if(GetHist(current_histname + "_JetsScaleDown"     + Cycle_name + samples_array[i])) JetsScale_error += pow(current_bin_content - GetHist(current_histname + "_JetsScaleDown"     + Cycle_name + samples_array[i]) -> GetBinContent(j), 2);
	  if(GetHist(current_histname + "_JetsScaleUp"       + Cycle_name + samples_array[i])) JetsScale_error += pow(current_bin_content - GetHist(current_histname + "_JetsScaleUp"     + Cycle_name + samples_array[i]) -> GetBinContent(j), 2);
	  JetsScale_error = pow(JetsScale_error, 0.5);
	}
	if(GetHist(current_histname + "_SD_JMR_Down"         + Cycle_name + samples_array[i]) || GetHist(current_histname + "_SD_JMR_Up"       + Cycle_name + samples_array[i]) ){
	  SD_JMR_error = 0.;
	  if(GetHist(current_histname + "_SD_JMR_Down"       + Cycle_name + samples_array[i])) SD_JMR_error += pow(current_bin_content - GetHist(current_histname + "_SD_JMR_Down"       + Cycle_name + samples_array[i]) -> GetBinContent(j), 2);
	  if(GetHist(current_histname + "_SD_JMR_Up"         + Cycle_name + samples_array[i])) SD_JMR_error += pow(current_bin_content - GetHist(current_histname + "_SD_JMR_Up"       + Cycle_name + samples_array[i]) -> GetBinContent(j), 2);
	  SD_JMR_error = pow(SD_JMR_error, 0.5);
	}
	if(GetHist(current_histname + "_SD_JMS_Down"         + Cycle_name + samples_array[i]) || GetHist(current_histname + "_SD_JMS_Up"       + Cycle_name + samples_array[i]) ){
	  SD_JMS_error = 0.;
	  if(GetHist(current_histname + "_SD_JMS_Down"       + Cycle_name + samples_array[i])) SD_JMS_error += pow(current_bin_content - GetHist(current_histname + "_SD_JMS_Down"       + Cycle_name + samples_array[i]) -> GetBinContent(j), 2);
	  if(GetHist(current_histname + "_SD_JMS_Up"         + Cycle_name + samples_array[i])) SD_JMS_error += pow(current_bin_content - GetHist(current_histname + "_SD_JMS_Up"       + Cycle_name + samples_array[i]) -> GetBinContent(j), 2);
	  SD_JMS_error = pow(SD_JMS_error, 0.5);
	}
	if(GetHist(current_histname + "_PUReweight_Down"         + Cycle_name + samples_array[i]) || GetHist(current_histname + "_PUReweight_Up"       + Cycle_name + samples_array[i]) ){
          PUReweight_error = 0.;
          if(GetHist(current_histname + "_PUReweight_Down"       + Cycle_name + samples_array[i])) PUReweight_error += pow(current_bin_content - GetHist(current_histname + "_PUReweight_Down"       + Cycle_name + samples_array[i]) -> GetBinContent(j), 2);
          if(GetHist(current_histname + "_PUReweight_Up"         + Cycle_name + samples_array[i])) PUReweight_error += pow(current_bin_content - GetHist(current_histname + "_PUReweight_Up"       + Cycle_name + samples_array[i]) -> GetBinContent(j), 2);
          PUReweight_error = pow(PUReweight_error, 0.5);
        }
	
	htmp -> SetBinError(j, stat_error);
	htmp_ElectronScale -> SetBinError(j, ElectronScale_error);
	htmp_ElectronSmear -> SetBinError(j, ElectronSmear_error);
	htmp_JetsRes -> SetBinError(j, JetsRes_error);
	htmp_JetsScale -> SetBinError(j, JetsScale_error);
	htmp_SD_JMR -> SetBinError(j, SD_JMR_error);
	htmp_SD_JMS -> SetBinError(j, SD_JMS_error);
	htmp_PUReweight -> SetBinError(j, PUReweight_error);

      }
      
      if(debug) cout << "rebinning" << endl;
      
      mapfunc[nameofhistogram + Cycle_name + samples_array[i] + "rebin"] = dynamic_cast<TH1F*>(htmp -> Rebin(N_bin, nameofhistogram + Cycle_name + samples_array[i] + "rebin", binx));
      GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin") -> SetFillColor(colour_array[i]);
      GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin") -> SetLineColor(colour_array[i]);
      
      mapfunc[nameofhistogram + Cycle_name + samples_array[i] + "rebin_ElectronScale"] = dynamic_cast<TH1F*>(htmp_ElectronScale -> Rebin(N_bin, nameofhistogram + Cycle_name + samples_array[i] + "rebin_ElectronScale", binx));
      mapfunc[nameofhistogram + Cycle_name + samples_array[i] + "rebin_ElectronSmear"] = dynamic_cast<TH1F*>(htmp_ElectronSmear -> Rebin(N_bin, nameofhistogram + Cycle_name + samples_array[i] + "rebin_ElectronSmear", binx));
      mapfunc[nameofhistogram + Cycle_name + samples_array[i] + "rebin_JetsRes"]       = dynamic_cast<TH1F*>(htmp_JetsRes       -> Rebin(N_bin, nameofhistogram + Cycle_name + samples_array[i] + "rebin_JetsRes"      , binx));
      mapfunc[nameofhistogram + Cycle_name + samples_array[i] + "rebin_JetsScale"]     = dynamic_cast<TH1F*>(htmp_JetsScale     -> Rebin(N_bin, nameofhistogram + Cycle_name + samples_array[i] + "rebin_JetsScale"    , binx));
      mapfunc[nameofhistogram + Cycle_name + samples_array[i] + "rebin_SD_JMR"]        = dynamic_cast<TH1F*>(htmp_SD_JMR        -> Rebin(N_bin, nameofhistogram + Cycle_name + samples_array[i] + "rebin_SD_JMR"       , binx));
      mapfunc[nameofhistogram + Cycle_name + samples_array[i] + "rebin_SD_JMS"]        = dynamic_cast<TH1F*>(htmp_SD_JMS        -> Rebin(N_bin, nameofhistogram + Cycle_name + samples_array[i] + "rebin_SD_JMS"       , binx));
      mapfunc[nameofhistogram + Cycle_name + samples_array[i] + "rebin_PUReweight"]    = dynamic_cast<TH1F*>(htmp_PUReweight    -> Rebin(N_bin, nameofhistogram + Cycle_name + samples_array[i] + "rebin_PUReweight"   , binx));

      
      if(debug) cout << samples_array[i] << endl;
      
      mapfunc[func + "rebin"] -> Add(GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin"));
      mapfunc[func + "syst_ElectronScale"] -> Add(GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin_ElectronScale"));
      mapfunc[func + "syst_ElectronSmear"] -> Add(GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin_ElectronSmear"));
      mapfunc[func + "syst_JetsRes"]       -> Add(GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin_JetsRes"));
      mapfunc[func + "syst_JetsScale"]     -> Add(GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin_JetsScale"));
      mapfunc[func + "syst_SD_JMR"]        -> Add(GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin_SD_JMR"));
      mapfunc[func + "syst_SD_JMS"]        -> Add(GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin_SD_JMS"));
      mapfunc[func + "syst_PUReweight"]    -> Add(GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin_PUReweight"));
    }
  }//for loop
  

  if(debug) cout << "4" << endl;
  
  for(int i = 0; i < N_syst_x; i++){
    double current_denom = mapfunc[func + "rebin"]              -> GetBinContent(i+1);
    cout << "current_denom : " << current_denom << endl;
    if(current_denom > 10e-6){ 
      syst_Stat[i]          = mapfunc[func + "rebin"]              -> GetBinError(i+1) / mapfunc[func + "rebin"]              -> GetBinContent(i+1);
      syst_ElectronScale[i] = mapfunc[func + "syst_ElectronScale"] -> GetBinError(i+1) / mapfunc[func + "syst_ElectronScale"] -> GetBinContent(i+1);
      syst_ElectronSmear[i] = mapfunc[func + "syst_ElectronSmear"] -> GetBinError(i+1) / mapfunc[func + "syst_ElectronSmear"] -> GetBinContent(i+1);
      syst_JetsRes[i]       = mapfunc[func + "syst_JetsRes"]       -> GetBinError(i+1) / mapfunc[func + "syst_JetsRes"]       -> GetBinContent(i+1);
      syst_JetsScale[i]     = mapfunc[func + "syst_JetsScale"]     -> GetBinError(i+1) / mapfunc[func + "syst_JetsScale"]     -> GetBinContent(i+1);
      syst_SD_JMR[i]        = mapfunc[func + "syst_SD_JMR"]        -> GetBinError(i+1) / mapfunc[func + "syst_SD_JMR"]        -> GetBinContent(i+1);
      syst_SD_JMS[i]        = mapfunc[func + "syst_SD_JMS"]        -> GetBinError(i+1) / mapfunc[func + "syst_SD_JMS"]        -> GetBinContent(i+1);
      syst_PUReweight[i]    = mapfunc[func + "syst_PUReweight"]    -> GetBinError(i+1) / mapfunc[func + "syst_PUReweight"]    -> GetBinContent(i+1);
    }
    else{
      syst_Stat[i]          = 0.;
      syst_ElectronScale[i] = 0.;
      syst_ElectronSmear[i] = 0.;
      syst_JetsRes[i]       = 0.;
      syst_JetsScale[i]     = 0.;
      syst_SD_JMR[i]        = 0.;
      syst_SD_JMS[i]        = 0.;
      syst_PUReweight[i]    = 0.;
    }
  }
  
  cout << "stat errpr" << endl;
  for(int i = 0; i < N_syst_x; i++){
    cout << syst_Stat[i] << ", ";
  }
  cout << "" << endl;

  cout << "syst_JetsScale" << endl;
  for(int i = 0; i < N_syst_x; i++){
    cout << syst_JetsScale[i] << ", ";
  }
  cout << "" << endl;
  
  TString name_x = nameofhistogram;
  TGraph *gr_stat_error         = new TGraph(N_syst_x, syst_x_bin, syst_Stat);
  gStyle->SetOptTitle(0);
  gr_stat_error -> SetMarkerStyle(20);
  gr_stat_error -> SetMarkerColor(kBlack);
  gr_stat_error -> SetLineColor(kBlack);
  gr_stat_error -> GetYaxis() -> SetLabelSize(0.03);
  gr_stat_error -> GetYaxis() -> SetTitle(title_y);
  gr_stat_error -> GetYaxis() -> SetTitleSize(0.07);
  gr_stat_error -> GetYaxis() -> SetTitleOffset(1.02);
  gr_stat_error -> GetXaxis() -> SetTitle(current_histname);
  gr_stat_error -> GetXaxis() -> SetTitleSize(0.15);
  gr_stat_error -> GetXaxis() -> SetLabelSize(0.03);
  

  gr_stat_error -> GetXaxis() -> SetRangeUser(x1_func, x2_func);
  gr_stat_error -> GetYaxis() -> SetRangeUser(0., 1.);

  gr_stat_error -> Draw("ALP");
  
  TGraph *gr_syst_ElectronScale = new TGraph(N_syst_x, syst_x_bin, syst_ElectronScale);
  gr_syst_ElectronScale -> SetMarkerStyle(20);
  gr_syst_ElectronScale -> SetMarkerColor(kBlue);
  gr_syst_ElectronScale -> SetLineColor(kBlue);
  gr_syst_ElectronScale -> Draw("same LP");
  
  TGraph *gr_syst_ElectronSmear = new TGraph(N_syst_x, syst_x_bin, syst_ElectronSmear);
  gr_syst_ElectronSmear -> SetMarkerStyle(20);
  gr_syst_ElectronSmear -> SetMarkerColor(kCyan);
  gr_syst_ElectronSmear -> SetLineColor(kCyan);
  gr_syst_ElectronSmear -> Draw("same LP");
  
  TGraph *gr_syst_JetsRes       = new TGraph(N_syst_x, syst_x_bin, syst_JetsRes);
  gr_syst_JetsRes -> SetMarkerStyle(20);
  gr_syst_JetsRes -> SetMarkerColor(kRed);
  gr_syst_JetsRes -> SetLineColor(kRed);
  gr_syst_JetsRes -> Draw("same LP");

  TGraph *gr_syst_JetsScale     = new TGraph(N_syst_x, syst_x_bin, syst_JetsScale);
  gr_syst_JetsScale -> SetMarkerStyle(20);
  gr_syst_JetsScale -> SetMarkerColor(kOrange);
  gr_syst_JetsScale -> SetLineColor(kOrange);
  gr_syst_JetsScale -> Draw("same LP");

  TGraph *gr_syst_SD_JMR     = new TGraph(N_syst_x, syst_x_bin, syst_SD_JMR);
  gr_syst_SD_JMR -> SetMarkerStyle(20);
  gr_syst_SD_JMR -> SetMarkerColor(kGreen);
  gr_syst_SD_JMR -> SetLineColor(kGreen);
  gr_syst_SD_JMR -> Draw("same LP");
  
  TGraph *gr_syst_SD_JMS     = new TGraph(N_syst_x, syst_x_bin, syst_SD_JMS);
  gr_syst_SD_JMS -> SetMarkerStyle(20);
  gr_syst_SD_JMS -> SetMarkerColor(kYellow);
  gr_syst_SD_JMS -> SetLineColor(kYellow);
  gr_syst_SD_JMS -> Draw("same LP");
  
  TGraph *gr_syst_PUReweight     = new TGraph(N_syst_x, syst_x_bin, syst_PUReweight);
  gr_syst_PUReweight -> SetMarkerStyle(28);
  gr_syst_PUReweight -> SetMarkerColor(kBlack);
  gr_syst_PUReweight -> SetLineColor(kBlack);
  gr_syst_PUReweight -> SetLineStyle(6);
  gr_syst_PUReweight -> Draw("same LP");
  
  maplegend[legend] -> AddEntry(gr_stat_error, "Stat error", "lp");
  maplegend[legend] -> AddEntry(gr_syst_PUReweight, "PUReweight error", "lp");
  maplegend[legend] -> AddEntry(gr_syst_ElectronScale, "ElectronScal error", "lp");
  maplegend[legend] -> AddEntry(gr_syst_ElectronSmear, "ElectronSmear error", "lp");
  maplegend[legend] -> AddEntry(gr_syst_JetsRes, "JetsRes error", "lp");
  maplegend[legend] -> AddEntry(gr_syst_JetsScale, "JetsScale error", "lp");
  maplegend[legend] -> AddEntry(gr_syst_SD_JMS, "SDMassScale error", "lp");
  maplegend[legend] -> AddEntry(gr_syst_SD_JMR, "SDMassResolution error", "lp");
  
  maplegend[legend] -> Draw("same");
  
  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_Lumi.SetTextSize(0.035);
  latex_Lumi.DrawLatex(0.7, 0.96, "41.3 fb^{-1} (13 TeV)");
  
  if(debug) cout << "9" << endl;


  TString pdfname;
  if(nameofhistogram.Contains("EMu")){
    pdfname = "./plots/EMu/";

  }
  else if(nameofhistogram.Contains("DiEle")){
    pdfname = "./plots/DiEle/";

  }
  else if(nameofhistogram.Contains("DiMu")){
    pdfname = "./plots/DiMu/";
  }
  else return;

  pdfname.Append(nameofhistogram);
  pdfname.Append("_syst.pdf");
  if(debug) cout << "9.1" << endl;

  mapcanvas[canvas] -> SaveAs(pdfname);
  
  if(debug) cout << "10" << endl;
 
}

// -- Setup hist name and open files
void draw_histogram(TString histname){
  
  TString current_histname = "";
  
  TObjArray *token = histname.Tokenize("_");
  int N_token = token -> GetEntries();
  current_histname = ((TObjString *)(token->At(0)))->String();
  for(int i_token = 1; i_token < N_token - 1; i_token++){
    current_histname = current_histname + "_" + ((TObjString *)(token->At(i_token)))->String();
  }
  
  if(debug) cout << "[[draw_histogram]] histname : " << current_histname << endl;
  
  unsigned int N_bin = map_bin_vector[histname].size();
  double current_bins[100];
  for(unsigned int i_bin = 0; i_bin < N_bin; i_bin++){
    current_bins[i_bin] = map_bin_vector[histname].at(i_bin);
  }
  
  
  // -- Get which region(CR1, CR2 , ...) and channel (DiEle, DiMu, ...)
  int N_regions = 6;

  TString regions[] = {"CR_Zmass",
		       "CR_ttbar",
		       "tight_CR_Zmass",
                       "tight_CR_ttbar",
		       "CR_inv_mll",
		       "SR"
  };
  /*
  TString regions[] = {                       "SR"
  };
  */
  int N_channels = 3;
  TString channels[] = {"DiEle",
                        "DiMu",
                        "EMu"
  };
  
  const int N_syst = 15;
  TString systematics[N_syst] = {"central",
				 "ElectronScaleDown", "ElectronScaleUp",
				 "ElectronSmearDown", "ElectronSmearUp",
				 "JetsResDown", "JetsResUp",
				 "JetsScaleDown", "JetsScaleUp",
				 "SD_JMR_Down", "SD_JMR_Up",
                                 "SD_JMS_Down", "SD_JMS_Up",
				 "PUReweight_Down", "PUReweight_Up",
				 //"PDFNormDown", "PDFNormUp"
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
  
  if(debug) cout << "current_dir : " << current_dir << "/" << current_histname << "+systflags" << endl;
  if(current_dir.Contains("empty") || current_channel.Contains("empty")) return;

  for(int i_syst = 0; i_syst < N_syst; i_syst++){
    TString current_dir_syst = current_dir + "_" + systematics[i_syst];
    TString current_hist_syst = current_histname + "_" + systematics[i_syst];

    openfile(Cycle_name, DoubleEG, current_dir_syst, current_hist_syst);
    openfile(Cycle_name, SingleMuon, current_dir_syst, current_hist_syst);
    openfile(Cycle_name, DY_high, current_dir_syst, current_hist_syst);
    openfile(Cycle_name, VV, current_dir_syst, current_hist_syst);
    openfile(Cycle_name, WJets, current_dir_syst, current_hist_syst);
    openfile(Cycle_name, ttbar, current_dir_syst, current_hist_syst);
  }
  make_histogram(histname, current_histname, N_bin, current_bins);
  
}

// -- Open binning txt file
void open_binning_file(TString filename){
  
  ifstream data_file;
  data_file.open(filename);
  
  char line[500];
  if(data_file.is_open()){
    while(!data_file.eof()){
      data_file.getline(line, 500);
      if(debug) cout << line << endl;
      TString this_line = line;
      if(this_line.Contains("resign")) continue;
      
      TObjArray *tx = this_line.Tokenize("\t");
      int N_line_part = tx -> GetEntries();
      if(N_line_part < 1) continue;
      if(debug){
	cout << "[[open_binning_file]] histname : " << ((TObjString *)(tx->At(0)))->String() << endl;
	cout << "[[open_binning_file]] binning : " << ((TObjString *)(tx->At(1)))->String() << endl;
      }
      
      // -- saving current histname for mapping : TString && clear the map vector
      TString current_histname = ((TObjString *)(tx->At(0)))->String();
      if(debug) cout << "[[open_binning_file]] current_histname : " << current_histname << endl;
      map_bin_vector[current_histname].clear();
      
      // -- making vector<double> for binning
      TString current_binning = ((TObjString *)(tx->At(1)))->String();
      TObjArray *binning_token = current_binning.Tokenize(", ");
      for(int i_bin_part = 0; i_bin_part < binning_token->GetEntries() - 1; i_bin_part++){
	TString current_bin_value = ((TObjString *)(binning_token->At(i_bin_part)))->String();
	if(i_bin_part == 0){
	  current_bin_value.Remove(0,1);
	}
	if(debug) cout << current_bin_value << " | " << endl;
	double current_bin_double = current_bin_value.Atof();
	map_bin_vector[current_histname].push_back(current_bin_double);
      }
      
      draw_histogram(current_histname);
      map_bin_vector.clear();
      
    }//while not end of file
  }//end if file open

  data_file.close();
}
//////////////////////////////////////////////////////////////////////////////

/// Main nction ////////////////////////////////////////////////////////////
void plot(){
  rootlogon();
  open_binning_file("binning_OS.txt");
  //open_binning_file("binning_SS.txt");

  
}
