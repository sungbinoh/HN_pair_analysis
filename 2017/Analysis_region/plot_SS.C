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

bool debug = true;

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
  cout << "[[openfile ]]opening : " << filename << endl;

  TFile *current_file = new TFile ((filename)) ;
  
  cout << "[[openfile ]] Cd : " << dir << endl;
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
void make_histogram(TString nameofhistogram, int N_bin, double binx[]){

  if(debug){
    cout << "[[make_histogram]] checking binning arrary" << endl;
    for(int i = 0; i < N_bin; i++){
      cout << binx[i] << ", " << endl;
    }
  }

  
  TString current_data;
  if(nameofhistogram.Contains("EMu") || nameofhistogram.Contains("DiMu")) current_data = SingleMuon;
  else if(nameofhistogram.Contains("DiEle")) current_data = DoubleEG;
  else return;
  
  TString title_y;
  if(nameofhistogram.Contains("N")) title_y = "Number";
  else title_y = "Events/GeV";
  
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
  gStyle -> SetOptStat(0);
  mapcanvas[canvas] -> SetTopMargin( 0.05 );
  mapcanvas[canvas] -> SetBottomMargin( 0.13 );
  mapcanvas[canvas] -> SetRightMargin( 0.05 );
  mapcanvas[canvas] -> SetLeftMargin( 0.16 );
  
  
  mappad[pad1] = new TPad("", "", 0, 0.25, 1, 1);
  gStyle -> SetOptStat(0);
  mappad[pad1] -> SetTopMargin( 0.07 );
  mappad[pad1] -> SetBottomMargin( 0.02 );
  mappad[pad1] -> SetLeftMargin( 0.15 );
  mappad[pad1] -> SetRightMargin( 0.03 );
  mappad[pad1] -> Draw();
  mappad[pad1] -> cd();
  //mappad[pad1] -> SetLogy();
  
  maplegend[legend] = new TLegend(0.60, 0.60, 0.96, 0.92);
  
  if(debug) cout << "1" << endl;
  
  maphstack[hstack] = new THStack(hstack, "Stacked_" + nameofhistogram);
  gStyle->SetOptTitle(0);
  
  
  int n_kind = 4;
  TString samples_array[] = {VV, WJets, ttbar, DY_high};
  Int_t colour_array[] = {419, 901, 416, 400};
  TString samples_legend[] = {"Other backgrounds", "W+jets", "ttbar", "Z/#gamma + jets"};
  
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
  gStyle->SetOptTitle(0);
  if(debug) cout << "func rebin rebinning" << endl;
  mapfunc[func + "rebin"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "rebin", binx);
  mapfunc[func + "blind_data"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "blind_Data", binx);

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
      
      mapfunc[func + "blind_data"]  -> Add(GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin"));
      // -- divide bin content and error by bin width
      for(int j = 1; j <= N_bin; j++){
	double current_content = GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin") -> GetBinContent(j);
	double current_error = GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin") -> GetBinError(j);
	GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin") -> SetBinContent(j, current_content / (binx[j] - binx[j - 1]) );
	GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin") -> SetBinError(j, current_error / (binx[j] - binx[j - 1]) );
      }

      maphstack[hstack] -> Add(GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin"));
      //maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin"), samples_legend[i], "lf");
      mapfunc[func + "rebin"] -> Add(GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin"));
    }
  }//for loop
  

  if(debug) cout << "4" << endl;

  mapfunc[nameofhistogram + Cycle_name + current_data + overflow] = new TH1F("", "", nx_func, x1_func, x2_func);
  for (Int_t i = 1; i <= nx_func; i++) {
    mapfunc[nameofhistogram + Cycle_name + current_data +overflow] -> SetBinContent(i, GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinContent(i)) ;
    mapfunc[nameofhistogram + Cycle_name + current_data +overflow] -> SetBinError(i, GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinError(i));
  }


  if(blind){
    mapfunc[nameofhistogram + Cycle_name + current_data + "rebin"] = (TH1F*)GetHist(func + "blind_data") -> Clone(clone + "blind");
  }
  else{
    mapfunc[nameofhistogram + Cycle_name + current_data + "rebin"] = dynamic_cast<TH1F*>(GetHist(nameofhistogram + Cycle_name + current_data + overflow) -> Rebin(N_bin, nameofhistogram + Cycle_name + current_data + "rebin", binx));
  }
  //GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> SetMarkerStyle(20);
  //GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> SetMarkerSize(1);
  if(!nameofhistogram.Contains("SR")){
    maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + current_data + "rebin"), "data", "lp");
  }
  //double data_max = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetMaximum();
  
  GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> SetMarkerStyle(20);
  GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> SetMarkerColor(kBlue);

  GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> SetMarkerSize(1);

  if(debug) cout << "5" << endl;
  
  mapfunc[func + "rebin"] -> SetFillColor(416);;
  mapfunc[func + "rebin"] -> SetFillStyle(3013);
  mapfunc[func + "rebin"] -> Draw("hist");
  mapfunc[func + "rebin"]  -> GetYaxis()->SetLabelSize(0.05);;
  mapfunc[func + "rebin"]  -> GetYaxis()->SetTitleSize(0.07);;
  mapfunc[func + "rebin"]  -> GetYaxis()->SetTitleOffset(1.02);;
  mapfunc[func + "rebin"]  -> GetXaxis()->SetLabelSize(0);
  mapfunc[func + "rebin"]  -> GetXaxis() -> SetRangeUser(x1_func, x2_func);
  mapfunc[func + "rebin"] -> GetXaxis() -> SetTitle(nameofhistogram);
  mapfunc[func + "rebin"] -> GetYaxis() -> SetTitle(title_y);
  mapfunc[func + "rebin"] -> Draw("histsame e");


  maplegend[legend] -> AddEntry(mapfunc[func + "rebin"], "Expected SS bkg from MC", "lf"); 
  

  /*
  maphstack[hstack] -> Draw("hist");
  maphstack[hstack] -> GetYaxis()->SetLabelSize(0.05);;
  maphstack[hstack] -> GetYaxis()->SetTitleSize(0.07);;
  maphstack[hstack] -> GetYaxis()->SetTitleOffset(1.02);;
  maphstack[hstack] -> GetXaxis()->SetLabelSize(0);
  maphstack[hstack] -> GetXaxis() -> SetRangeUser(x1_func, x2_func);
  //maphstack[hstack] -> SetMaximum(data_max * 10000.);
  //maphstack[hstack] -> SetMinimum();

  maphstack[hstack] -> GetXaxis() -> SetTitle(nameofhistogram);

  maphstack[hstack] -> GetYaxis() -> SetTitle(title_y);
  maphstack[hstack] -> Draw("histsame");
  */
  vector<double> vx, vy, vexl, vexh, veyl, veyh;
  //GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> Draw("histsamep9");
  
  map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] = new TGraphAsymmErrors(GetHist(nameofhistogram + Cycle_name + current_data + "rebin"));
  //make correct error bars 
  
  for(int i = 0; i < N_bin; i++){
    int N = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetBinContent(i + 1);
    double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
    double U =  (N==0) ? ( ROOT::Math::gamma_quantile_c(alpha,N+1,1) ) : ( ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) );
    if( N!=0 ){
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEYlow(i, (N-L) / (binx[i + 1] - binx[i]) );
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEXlow(i, 0);
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEYhigh(i, (U-N ) / (binx[i + 1] - binx[i]) );
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEXhigh(i, 0);
      
      double current_x = -1., current_y = -1.;
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> GetPoint(i, current_x, current_y);
      double modified_y = -1.;
      modified_y = current_y / ( binx[i + 1] - binx[i] );
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPoint(i, current_x, modified_y);
      
      
      if(debug) cout << "i : " << i << ", current_x : " << current_x << ", current_y : " << current_y << ", modified_y : " << modified_y << endl; 
      
      veyl.push_back( (N-L) / (binx[i + 1] - binx[i]) );
      veyh.push_back( (U-N) / (binx[i + 1] - binx[i]));
    }
    else{
      double zerodata_err_low = 0.1 / (binx[i + 1] - binx[i]);
      double zerodata_err_high = 1.8 / (binx[i + 1] - binx[i]);
      
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
  
  //map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> Draw("p0same");
  
  for (int i = 1; i <= N_bin ; i++) {
    double modified_content = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetBinContent(i) / (binx[i] - binx[i - 1]);
    double modified_error = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetBinError(i) / (binx[i] - binx[i - 1]);;
    GetHist(nameofhistogram + Cycle_name + current_data + "rebin")  -> SetBinContent(i, modified_content );
    GetHist(nameofhistogram + Cycle_name + current_data + "rebin")  -> SetBinError(i, modified_error );
  }
  
  GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetXaxis() -> SetRangeUser(x1_func, x2_func);
  //GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> SetMinimum(1.0);
  
  mappad[pad1] -> Update();
  

  GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> SetMarkerColor(kPink);
  GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> Draw("histsamep");
  
  
  cout << "mapfunc[func + rebin] -> GetBinContent(11) : " << mapfunc[func + "rebin"] -> GetBinContent(11) << endl;
  TH1 *last = (TH1*)maphstack[hstack] ->GetStack()->Last();
  cout << "last->GetBinContent(11) : " << last->GetBinContent(11) << endl;
  
  double data_max = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetMaximum();
  double data_min = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetMinimum();

  mapfunc[func + "rebin"] -> SetMaximum(data_max * 3.);
  mapfunc[func + "rebin"] -> SetMinimum(0.);
  
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kBlack);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(1001);
  maplegend[legend] -> SetShadowColor(0); // 0 = transparent
  
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> Draw("same");

  if(debug) cout << "6" << endl;

  mapcanvas[canvas] -> cd();

  mappad[pad2] = new TPad(pad2, "", 0, 0, 1, 0.25);
  mappad[pad2] -> SetTopMargin( 0.03 );
  mappad[pad2] -> SetBottomMargin( 0.4 );
  mappad[pad2] -> SetLeftMargin( 0.15 );
  mappad[pad2] -> SetRightMargin( 0.03 );
  mappad[pad2] -> Draw();
  mappad[pad2] -> cd();

  mapfunc[clone] = (TH1F*)GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> Clone(clone);

  mapfunc["stat" + nameofhistogram] = (TH1F*)GetHist(func + "rebin") -> Clone(clone);
  
  if(debug) cout << "7" << endl;

  mapline[line] = new TLine(x1_func, 1, x2_func, 1);
  mapline[line] -> Draw();
  mapline[line] -> SetLineStyle(1);
  mapline[line] -> SetLineColor(kRed);

  TString name_x = nameofhistogram;
  mapfunc["stat" + nameofhistogram] -> SetTitle("");
  mapfunc["stat" + nameofhistogram] -> SetYTitle("#frac{Obs.}{Pred.}");
  mapfunc["stat" + nameofhistogram] -> SetFillColor(kOrange);
  mapfunc["stat" + nameofhistogram] -> SetMarkerSize(0);
  mapfunc["stat" + nameofhistogram] -> SetMarkerStyle(0);
  mapfunc["stat" + nameofhistogram] -> SetLineColor(kWhite);
  mapfunc["stat" + nameofhistogram] -> GetXaxis() -> SetTitle(name_x);
  mapfunc["stat" + nameofhistogram] -> GetXaxis() -> SetLabelSize(0.10);
  mapfunc["stat" + nameofhistogram] -> GetYaxis() -> SetLabelSize(0.08);
  mapfunc["stat" + nameofhistogram] -> SetFillStyle(1001);
  mapfunc["stat" + nameofhistogram] -> SetFillColorAlpha(45,0.35);
  mapfunc["stat" + nameofhistogram] -> GetXaxis() -> SetRangeUser(x1_func, x2_func);
  mapfunc["stat" + nameofhistogram] -> GetYaxis() -> SetTitleOffset(0.5);
  mapfunc["stat" + nameofhistogram] -> GetXaxis() -> SetTitleSize(0.15);
  mapfunc["stat" + nameofhistogram] -> GetYaxis() -> SetTitleSize(0.12);
  mapfunc["stat" + nameofhistogram] -> SetMinimum(0.0);
  mapfunc["stat" + nameofhistogram] -> SetMaximum(5.0);
  mapfunc["stat" + nameofhistogram] -> SetStats(0);

  if(debug) cout << "8" << endl;
  
  //mapfunc["stat" + nameofhistogram] -> Divide(mapfunc[func + "rebin"]);
  Int_t ncells = mapfunc["stat" + nameofhistogram] -> GetSize();
  
  for(int i = 0; i < N_bin; i++){
    double current_error = mapfunc["stat" + nameofhistogram] -> GetBinError(i + 1) / mapfunc["stat" + nameofhistogram] -> GetBinContent(i + 1);
    double original_error = mapfunc["stat" + nameofhistogram] -> GetBinError(i + 1);
    double original_content = mapfunc["stat" + nameofhistogram] -> GetBinContent(i + 1);

    cout << "i : " << i << ", current_error : " << current_error << ", original_content : " << original_content << ", original_error : " << original_error << endl;
    bool no_data = false;
    double data_value = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetBinContent(i + 1);;
    if(data_value < 10e-7) current_error = 0.;
    
    mapfunc["stat" + nameofhistogram] -> SetBinContent(i + 1, 1.);
    mapfunc["stat" + nameofhistogram] -> SetBinError(i + 1, current_error);
    
  }
  mapfunc["stat" + nameofhistogram] -> Draw("CE2");
  mapfunc[clone] -> Divide(mapfunc[func + "rebin"]);
  //if(!blind){
    for(int i = 0; i < N_bin; i++){
      double bkg_value = mapfunc[func + "rebin"] -> GetBinContent(i + 1);
      double data_value = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetBinContent(i + 1);
      if(debug) cout << " data_value : " << data_value << endl;
      if(debug) cout << " GetHist(nameofhistogram + Cycle_name + current_data + rebin) -> GetBinContent(" << i + 1 << ") : " << GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetBinContent(i + 1) << endl;
      if(bkg_value < 10e-7){
        mapfunc[clone] -> SetBinContent(i + 1, 999.);
      }
      else if(data_value < 10e-7) mapfunc[clone] -> SetBinContent(i + 1, 0.);
    }
    // -- make correct errorbars
    map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] = new TGraphAsymmErrors(mapfunc[clone]);
    
    for(int i = 0; i < N_bin; i++){
      double bkg_i = mapfunc[func + "rebin"] -> GetBinContent(i+1);
      cout << "veyl.at(i) : " << veyl.at(i) << ", bkg_i : " << bkg_i << ", veyl.at(i) / bkg_i : " << veyl.at(i) / bkg_i << endl;
      
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] -> SetPointEYlow(i, veyl.at(i) / bkg_i);
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] -> SetPointEXlow(i, 0.);
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] -> SetPointEYhigh(i, veyh.at(i) / bkg_i);
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] -> SetPointEXhigh(i, 0.);
      if(bkg_i < 10e-7){
        map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] -> SetPointEYlow(i, 0.);
        map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] -> SetPointEXlow(i, 0.);
        map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] -> SetPointEYhigh(i, 0.);
	map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] -> SetPointEXhigh(i, 0.);
      }
      
    }//end for i on N_bin
    
    //map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] -> Draw("p0same");
    //}
  mapline[line] -> Draw("same");
  
  maplegend["bottom" + legend] = new TLegend(0.2, 0.85, 0.4, 0.95);
  maplegend["bottom" + legend]->SetBorderSize(0);
  maplegend["bottom" + legend]->SetNColumns(2);
  maplegend["bottom" + legend]->AddEntry(mapfunc["stat" + nameofhistogram], "Stat.", "f");
  maplegend["bottom" + legend]->AddEntry(mapfunc[clone] , "Obs./Pred.", "p");
  maplegend["bottom" + legend]->Draw("same");
  
  mapcanvas[canvas] -> cd();
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
  pdfname.Append(".pdf");
  if(debug) cout << "9.1" << endl;

  mapcanvas[canvas] -> SaveAs(pdfname);
  
  if(debug) cout << "10" << endl;
    
}

// -- Setup hist name and open files
void draw_histogram(TString histname){
  
  if(debug) cout << "[[draw_histogram]] histname : " << histname << endl;
  
  unsigned int N_bin = map_bin_vector[histname].size();
  double current_bins[100];
  for(unsigned int i_bin = 0; i_bin < N_bin; i_bin++){
    current_bins[i_bin] = map_bin_vector[histname].at(i_bin);
  }
  
  
  // -- Get which region(CR1, CR2 , ...) and channel (DiEle, DiMu, ...)
  int N_regions = 4;

  TString regions[] = {"CR_Zmass",
		       "CR_ttbar",
		       "CR_inv_mll",
		       "SR"
  };
  /*
  TString regions[] = {                       "SR"
  };
  */
  int N_channels = 6;
  TString channels[] = {"DiEle",
                        "DiMu",
                        "EMu",
			"DiEle_SS",
                        "DiMu_SS",
                        "EMu_SS"
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
  openfile(Cycle_name, VV, current_dir,  histname);
  openfile(Cycle_name, WJets, current_dir,  histname);
  openfile(Cycle_name, ttbar, current_dir, histname);
  
  make_histogram(histname, N_bin, current_bins);
  
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
  //open_binning_file("binning_OS.txt");
  open_binning_file("binning_SS.txt");

  
}
