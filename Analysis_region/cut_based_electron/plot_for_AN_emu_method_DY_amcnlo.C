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
double emu_over_ee = 0.451;
double emu_over_mumu = 0.647;

const double alpha = 1 - 0.6827;

/// Getting Histogram Function ///////////////////////////////////////////////
TH1 * GetHist(TString hname){
  
  TH1 *h = NULL;
  std::map<TString, TH1*>::iterator mapit = maphist.find(hname);
  if(mapit != maphist.end()) return mapit-> second;
  
  return h;
  
}
//////////////////////////////////////////////////////////////////////////////


/// Open ROOT file ///////////////////////////////////////////////////////////
void openfile(TString cyclename, TString samplename){
  
  TString underbar = "_";
  TString version = "dilep_cat_v8-0-7.root";
  TString filename = cyclename + underbar + samplename + underbar + version;
  if(samplename.Contains("data") || samplename.Contains("NNPDF30_13TeV-powheg")) filename = cyclename + underbar + samplename + "_cat_v8-0-7.root";
  
  cout << "opening : " << filename << endl;
  
  mapfile[filename] = new TFile ((filename)) ;
  TString regions[7] = {"CR1", "CR2", "CR3", "CR4", "SR1"};
  TString channels[3] = {"DiEle", "DiMu", "EMu"};
  TString charges[2] = {"OS", "SS"};
  int i_dir = 0;
  TString directories[6];
  for(int j = 0; j < 3; j++){
    for(int k = 0; k < 2; k++){
      directories[i_dir] = charges[k] + "_" + channels[j];
      i_dir++;
    }
  }
  
  int N_regions = 5;
  int N_directories = 6;
  
  for(int i = 0; i < N_regions; i++){
    gDirectory->cd(regions[i]);
    for(int j = 0; j < N_directories; j++){
      gDirectory->cd(regions[i] + "_" + directories[j]);
      TIter next(gDirectory->GetListOfKeys());
      TKey *key;
      vector<TString> histnames;
      while ((key = (TKey*)next())) {
	TClass *cl = gROOT->GetClass(key->GetClassName());
	if (!cl->InheritsFrom("TH1")) continue;
	histnames.push_back(key -> GetName());
      }
      
      for(unsigned int k = 0; k < histnames.size(); k ++){
	maphist[histnames.at(k) + cyclename + samplename] = (TH1*)gDirectory -> Get(histnames.at(k));
      }
      gDirectory->cd("../");
    }//for dir
    if(regions[i].Contains("SR")){
      for(int j = 0; j < N_directories; j++){
	gDirectory->cd(regions[i] + "_" + directories[j] + "_HNcut");
	TIter next(gDirectory->GetListOfKeys());
	TKey *key;
	vector<TString> histnames;
	while ((key = (TKey*)next())) {
	  TClass *cl = gROOT->GetClass(key->GetClassName());
	  if (!cl->InheritsFrom("TH1")) continue;
	  histnames.push_back(key -> GetName());
	}

	for(unsigned int k = 0; k < histnames.size(); k ++){
	  maphist[histnames.at(k) + cyclename + samplename] = (TH1*)gDirectory -> Get(histnames.at(k));
	}
	gDirectory->cd("../");
      }
    }// if SR
    gDirectory->cd("../");
  }//for region
}
//////////////////////////////////////////////////////////////////////////////

/// Open Signal sample root file
void openfile_signal(TString cyclename, TString samplename, TString channel){
  TString directory;
  if(channel.Contains("MuMu")) directory = "MM";
  else if(channel.Contains("ElEl")) directory = "EE";
  else return;

  TString underbar = "_";
  TString version = "cat_v8-0-7.root";
  TString filename = cyclename + underbar + "HNpair_" + channel + "_WR5000_" + samplename +  "_official_cat_v8-0-7.root";

  cout << "opening : " << filename << endl;
  mapfile[filename] = new TFile (("./signal/" + directory + "/" + filename)) ;

  TString regions[7] = {"CR1", "CR2", "CR3", "CR4", "SR1"};
  TString channels[3] = {"DiEle", "DiMu", "EMu"};
  TString charges[2] = {"OS", "SS"};
  int i_dir = 0;
  TString directories[6];
  for(int j = 0; j < 3; j++){
    for(int k = 0; k < 2; k++){
      directories[i_dir] = charges[k] + "_" + channels[j];
      i_dir++;
    }
  }

  int N_regions = 5;
  int N_directories = 6;

  for(int i = 0; i < N_regions; i++){
    gDirectory->cd(regions[i]);
    for(int j = 0; j < N_directories; j++){
      gDirectory->cd(regions[i] + "_" + directories[j]);
      TIter next(gDirectory->GetListOfKeys());
      TKey *key;
      vector<TString> histnames;
      while ((key = (TKey*)next())) {
        TClass *cl = gROOT->GetClass(key->GetClassName());
	if (!cl->InheritsFrom("TH1")) continue;
	histnames.push_back(key -> GetName());
      }

      for(unsigned int k = 0; k < histnames.size(); k ++){
	maphist[histnames.at(k) + channel + samplename] = (TH1*)gDirectory -> Get(histnames.at(k));
      }
      gDirectory->cd("../");
    }
    if(regions[i].Contains("SR")){
      for(int j = 0; j < N_directories; j++){
        gDirectory->cd(regions[i] + "_" + directories[j] + "_HNcut");
        TIter next(gDirectory->GetListOfKeys());
        TKey *key;
        vector<TString> histnames;
        while ((key = (TKey*)next())) {
          TClass *cl = gROOT->GetClass(key->GetClassName());
          if (!cl->InheritsFrom("TH1")) continue;
	  histnames.push_back(key -> GetName());
	}

        for(unsigned int k = 0; k < histnames.size(); k ++){
          maphist[histnames.at(k) + channel + samplename] = (TH1*)gDirectory -> Get(histnames.at(k));
        }
        gDirectory->cd("../");
      }
    }
    gDirectory->cd("../");
  }
  
  gDirectory->cd("Hists");
  TIter next(gDirectory->GetListOfKeys());
  TKey *key;
  vector<TString> histnames;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    histnames.push_back(key -> GetName());
  }

  for(unsigned int i = 0; i < histnames.size(); i ++){
    maphist[histnames.at(i) + channel + samplename] = (TH1*)gDirectory -> Get(histnames.at(i));
  }

  double denom;
  denom = maphist["signal_eff" + channel + samplename] -> GetBinContent(1);
  maphist["denom" + channel + samplename] = new TH1D("denom", "denom", 1, 0., 1.);
  maphist["denom" + channel + samplename] -> Fill(0.5, denom);
  
}
//////////////////////////////////////////////////////////////////////////////


void make_hist_after_emu_method(TString nameofhistogram, TString channel, TString region){
  
  double current_emu_ratio = 1.;
  if(channel.Contains("DiEle")) current_emu_ratio = emu_over_ee;
  else if(channel.Contains("DiMu")) current_emu_ratio = emu_over_mumu;

  TString current_data = SingleMuon;
  TString cyclename = Cycle_name;

  cout << "making data driven plot of " + nameofhistogram + "_" + region + "_" + channel +  cyclename + "emu" << endl;
  
  TString EMu = "_EMu";
  if(channel.Contains("HNcut")) EMu = EMu + "_HNcut";
  
  //clone emu data plot
  mapfunc[nameofhistogram + "_" + region + "_" + channel +  + cyclename + "emu"] = (TH1F*)GetHist(nameofhistogram +"_" + region + EMu + Cycle_name + current_data) -> Clone(nameofhistogram + "_" + region + "_" + channel +  + cyclename + "emu");
  
  TString name_cycle = nameofhistogram + "_" + region + "_EMu" +  + Cycle_name;
  if(channel.Contains("HNcut")) name_cycle = nameofhistogram + "_" + region + "_EMu_HNcut" +  + Cycle_name;
  //subtract non flavour symmetric bkgs
  mapfunc[nameofhistogram + "_" + region + "_" + channel +  + cyclename + "emu"] -> Add(GetHist(name_cycle +  WZ), -1);
  mapfunc[nameofhistogram + "_" + region + "_" + channel +  + cyclename + "emu"] -> Add(GetHist(name_cycle +  ZZ), -1);
  mapfunc[nameofhistogram + "_" + region + "_" + channel +  + cyclename + "emu"] -> Add(GetHist(name_cycle +  SingleTop_s), -1);
  mapfunc[nameofhistogram + "_" + region + "_" + channel +  + cyclename + "emu"] -> Add(GetHist(name_cycle +  SingleTbar_t), -1);
  mapfunc[nameofhistogram + "_" + region + "_" + channel +  + cyclename + "emu"] -> Add(GetHist(name_cycle +  SingleTop_t), -1);
  mapfunc[nameofhistogram + "_" + region + "_" + channel +  + cyclename + "emu"] -> Add(GetHist(name_cycle +  WGtoLNuG), -1);
  mapfunc[nameofhistogram + "_" + region + "_" + channel +  + cyclename + "emu"] -> Add(GetHist(name_cycle +  ZGto2LG), -1);
  mapfunc[nameofhistogram + "_" + region + "_" + channel +  + cyclename + "emu"] -> Add(GetHist(name_cycle +  WWW), -1);
  mapfunc[nameofhistogram + "_" + region + "_" + channel +  + cyclename + "emu"] -> Add(GetHist(name_cycle +  WWZ), -1);
  mapfunc[nameofhistogram + "_" + region + "_" + channel +  + cyclename + "emu"] -> Add(GetHist(name_cycle +  WZZ), -1);
  mapfunc[nameofhistogram + "_" + region + "_" + channel +  + cyclename + "emu"] -> Add(GetHist(name_cycle +  ZZZ), -1);
  
  mapfunc[nameofhistogram + "_" + region + "_" + channel +  + cyclename + "emu"] -> Scale(current_emu_ratio);

  maphist[nameofhistogram + "_" + region + "_" + channel +  + cyclename + "emu"] = mapfunc[nameofhistogram + "_" + region + "_" + channel +  + cyclename + "emu"];
  
}


//////////////////////////////////////////////////////////////////////////////
void makehistogram_variable_bin(TString nameofhistogram, float xmin, float xmax, float ymax, double binx[], int N_bin, TString name_x, bool name_y, TString channel){

  double current_emu_ratio = 1.;
  if(channel.Contains("DiEle")) current_emu_ratio = emu_over_ee;
  else if(channel.Contains("DiMu")) current_emu_ratio = emu_over_mumu;
  
  double additional_weight = 1.;
  //if(channel.Contains("DiEle")) additional_weight = 20150.390 / 35820.052000;

  bool blind = false;
  blind = (channel.Contains("SR")) && (!channel.Contains("EMu"));

  TString current_data;
  if(channel.Contains("EMu") || channel.Contains("DiMu")) current_data = SingleMuon;
  if(channel.Contains("DiEle")) current_data = DoubleEG;
  

  TString title_y;
  if(name_y) title_y = "Events/bin";
  else title_y = "Number";
  
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


  mappad[pad1] = new TPad("", "", 0, 0.25, 1, 1);
  mappad[pad1] -> SetTopMargin( 0.07 );
  mappad[pad1] -> SetBottomMargin( 0.02 );
  mappad[pad1] -> SetLeftMargin( 0.15 );
  mappad[pad1] -> SetRightMargin( 0.03 );
  mappad[pad1] -> Draw();
  mappad[pad1] -> cd();
  mappad[pad1] -> SetLogy();

  maplegend[legend] = new TLegend(0.60, 0.60, 0.96, 0.92);
  
  cout << "1" << endl;

  maphstack[hstack] = new THStack(hstack, "Stacked_" + nameofhistogram);
  gStyle->SetOptTitle(0);
  
  
  int n_kind = 3;
  TString samples_array[] = {WZ, "emu", DY_high};
  Int_t colour_array[] = {419, 416, 400};
  TString samples_legend[] = {"Other backgrounds", "Data Driven Flavour symm. bkg", "Z/#gamma + jets"};
  
  cout << "2" << endl;
  
  TString name_cycle = nameofhistogram + Cycle_name;
  cout << "check1" << endl;
  GetHist(name_cycle +  DY_high) -> Add(GetHist(name_cycle + DY_low));
  
  cout << "2.1" << endl;
  
  GetHist(name_cycle +  SingleTop_s) -> Add(GetHist(name_cycle + SingleTbar_t));
  GetHist(name_cycle +  SingleTop_s) -> Add(GetHist(name_cycle + SingleTop_t));
  GetHist(name_cycle +  WZ) -> Add(GetHist(name_cycle + ZZ));
  GetHist(name_cycle +  WGtoLNuG) -> Add(GetHist(name_cycle + ZGto2LG));
  GetHist(name_cycle +  WWW) -> Add(GetHist(name_cycle + WWZ));
  GetHist(name_cycle +  WWW) -> Add(GetHist(name_cycle + WZZ));
  GetHist(name_cycle +  WWW) -> Add(GetHist(name_cycle + ZZZ));
  maphist[name_cycle + WZ] -> Add(GetHist(name_cycle +  SingleTop_s));
  maphist[name_cycle + WZ] -> Add(GetHist(name_cycle +  WWW));
  maphist[name_cycle + WZ] -> Add(GetHist(name_cycle +  WGtoLNuG));
  
  
  cout << "3" << endl;

  TString overflow = "overflow";
  Int_t nx_func    = GetHist(nameofhistogram + Cycle_name + current_data) -> GetNbinsX()+1;
  Double_t x1_func = GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinLowEdge(1);
  Double_t bw_func = GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinWidth(nx_func);
  Double_t x2_func = GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinLowEdge(nx_func)+bw_func;
  mapfunc[func] = new TH1F("", "", nx_func, x1_func, x2_func);
  mapfunc[func + "rebin"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "rebin", binx);

  for(int i = 0; i < n_kind; i++){
    cout << samples_array[i] << endl;
    if(maphist[nameofhistogram + Cycle_name + samples_array[i]]){
      Int_t nx    = GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetNbinsX()+1;
      Double_t x1 = GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetBinLowEdge(1);
      Double_t bw = GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetBinWidth(nx);
      Double_t x2 = GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetBinLowEdge(nx)+bw;
      //cout << "nx : " << nx << ", x1 : " << x1 << ", bw : " << bw << ", x2 : " << x2 <<  endl;

      TH1F *htmp = new TH1F("", "", nx, x1, x2);
      
      for (Int_t j = 1; j <= nx; j++) {
	//cout << "3.1" << endl;
	//cout << GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetBinContent(j) << endl;
	htmp -> SetBinContent(j, GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetBinContent(j));
	htmp -> SetBinError(j, GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetBinError(j));

	//maphist[nameofhistogram + Cycle_name + samples_array[i] +overflow] -> Fill(maphist[nameofhistogram + Cycle_name + samples_array[i] +overflow] -> GetBinCenter(i), GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> GetBinContent(i));
	//cout << "3.2" << endl;
      }
      
      maphist[nameofhistogram + Cycle_name + samples_array[i] + "rebin"] = htmp -> Rebin(N_bin, nameofhistogram + Cycle_name + samples_array[i] + "rebin", binx);
      GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin") -> SetFillColor(colour_array[i]);
      GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin") -> SetLineColor(colour_array[i]);
      GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin") -> Scale(additional_weight);
      cout << samples_array[i] << endl;
      maphstack[hstack] -> Add(GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin"));
      maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin"), samples_legend[i], "lf");
      mapfunc[func + "rebin"] -> Add(GetHist(nameofhistogram + Cycle_name + samples_array[i] + "rebin"));
    }
  }//for loop                                                                                                                                                                                                                                                                   
  
  cout << "4" << endl;
  
  maphist[nameofhistogram + Cycle_name + current_data + overflow] = new TH1F("", "", nx_func, x1_func, x2_func);
  for (Int_t i = 1; i <= nx_func; i++) {
    maphist[nameofhistogram + Cycle_name + current_data +overflow] -> SetBinContent(i, GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinContent(i));
    maphist[nameofhistogram + Cycle_name + current_data +overflow] -> SetBinError(i, GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinError(i));
  }
   
  maphist[nameofhistogram + Cycle_name + current_data + "rebin"] = GetHist(nameofhistogram + Cycle_name + current_data + overflow) -> Rebin(N_bin, nameofhistogram + Cycle_name + current_data + "rebin", binx);
  GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> SetMarkerStyle(20);
  GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> SetMarkerSize(1);
  if(!channel.Contains("SR")){
    maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + current_data + "rebin"), "data", "lp");
  }
  double data_max = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetMaximum();
  
  cout << "5" << endl;


  maphstack[hstack] -> Draw("hist");
  maphstack[hstack] -> GetYaxis()->SetLabelSize(0.05);;
  maphstack[hstack] -> GetYaxis()->SetTitleSize(0.07);;
  maphstack[hstack] -> GetYaxis()->SetTitleOffset(1.02);;
  maphstack[hstack] -> GetXaxis()->SetLabelSize(0);
  maphstack[hstack] -> GetXaxis() -> SetRangeUser(xmin, xmax);
  maphstack[hstack] -> SetMaximum(data_max * 100.);
  maphstack[hstack] -> SetMinimum(1.0);
  //maphstack[hstack] -> SetMinimum(0.001);                                                                                                                                                                                                                                     
  maphstack[hstack] -> GetXaxis() -> SetTitle(nameofhistogram);
  //maphstack[hstack] -> GetYaxis() -> SetTitle("Entries / 25 GeV");                                                                                                                                                                                                            
  maphstack[hstack] -> GetYaxis() -> SetTitle(title_y);
  maphstack[hstack] -> Draw("histsame");
  
  vector<double> vx, vy, vexl, vexh, veyl, veyh;
  if(!blind){
    GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> Draw("histsamep9");
    
    map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] = new TGraphAsymmErrors(GetHist(nameofhistogram + Cycle_name + current_data + "rebin"));
    //make correct error bars
    
    for(unsigned int i = 0; i < N_bin; i++){
      int N = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetBinContent(i + 1);
      double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
      double U =  (N==0) ? ( ROOT::Math::gamma_quantile_c(alpha,N+1,1) ) : ( ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) );
      if( N!=0 ){
	map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEYlow(i, N-L );
	map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEXlow(i, 0);
	map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEYhigh(i, U-N );
	map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEXhigh(i, 0);
	
	veyl.push_back(N-L);
	veyh.push_back(U-N);
      }
      else{
	double zerodata_err_low = 0.1;
	double zerodata_err_high = 1.8;
	
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
    map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> Draw("p0same");
  }
  GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetXaxis() -> SetRangeUser(xmin, xmax);
  GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> SetMinimum(1.0);
  
  mappad[pad1] -> Update();
  
  
  if(!blind){
    GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> Draw("histsamep");
  }
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kBlack);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(1001);
  maplegend[legend] -> SetShadowColor(0); // 0 = transparent                                                                                                                                                                                                                   
  
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> Draw("same");
  
  //maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + data), "data", "p");                                                                                                                                                                                   
  cout << "6" << endl;

  mapcanvas[canvas] -> cd();

  mappad[pad2] = new TPad(pad2, "", 0, 0, 1, 0.25);
  mappad[pad2] -> SetTopMargin( 0.03 );
  mappad[pad2] -> SetBottomMargin( 0.4 );
  mappad[pad2] -> SetLeftMargin( 0.15 );
  mappad[pad2] -> SetRightMargin( 0.03 );
  mappad[pad2] -> Draw();
  mappad[pad2] -> cd();

  mapfunc[clone] = (TH1F*)GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> Clone(clone);
   
  mapfunc["stat" + nameofhistogram] = (TH1F*)GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> Clone(clone);

  cout << "7" << endl;


  mapline[line] = new TLine(xmin, 1, xmax, 1);
  mapline[line] -> Draw();
  mapline[line] -> SetLineStyle(1);
  mapline[line] -> SetLineColor(kRed);

  mapfunc["stat" + nameofhistogram] -> SetTitle("");
  mapfunc["stat" + nameofhistogram] -> SetYTitle("#frac{Obs.}{Pred.}");
  mapfunc["stat" + nameofhistogram] -> SetFillColor(kOrange);
  mapfunc["stat" + nameofhistogram] -> SetMarkerSize(0);
  mapfunc["stat" + nameofhistogram] -> SetMarkerStyle(0);
  mapfunc["stat" + nameofhistogram] -> SetLineColor(kWhite);
  mapfunc["stat" + nameofhistogram] -> GetXaxis() -> SetTitle(name_x);
  mapfunc["stat" + nameofhistogram] -> GetXaxis() -> SetLabelSize(0.10);
  mapfunc["stat" + nameofhistogram] -> GetYaxis() -> SetLabelSize(0.08);
  mapfunc["stat" + nameofhistogram] -> SetFillColorAlpha(45,0.35);
  mapfunc["stat" + nameofhistogram] -> GetXaxis() -> SetRangeUser(xmin, xmax);
  mapfunc["stat" + nameofhistogram] -> GetYaxis() -> SetTitleOffset(0.5);
  mapfunc["stat" + nameofhistogram] -> GetXaxis() -> SetTitleSize(0.15);
  mapfunc["stat" + nameofhistogram] -> GetYaxis() -> SetTitleSize(0.12);
  mapfunc["stat" + nameofhistogram] -> SetMinimum(0.0);
  mapfunc["stat" + nameofhistogram] -> SetMaximum(5.0);
  mapfunc["stat" + nameofhistogram] -> SetStats(0);
  
  cout << "8" << endl;

  
  
  mapfunc["stat" + nameofhistogram] -> Divide(mapfunc[func + "rebin"]);
  Int_t ncells = mapfunc["stat" + nameofhistogram] -> GetSize();
  for(int i = 0; i < ncells - 1; i++){
    mapfunc["stat" + nameofhistogram] -> SetBinContent(i, 1.);
  }
  mapfunc["stat" + nameofhistogram] -> Draw("CE2");
  mapfunc[clone] -> Divide(mapfunc[func + "rebin"]);
  if(!blind){
    for(unsigned int i = 0; i < N_bin; i++){
      double bkg_value = mapfunc[func + "rebin"] -> GetBinContent(i + 1);
      int data_value = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetBinContent(i + 1);
      //cout << "bkg_value : " << bkg_value << endl;
      if(bkg_value < 10e-7){ 
	//cout << "pass zero" << endl;
	mapfunc[clone] -> SetBinContent(i + 1, 999.);
      }
      else if(data_value < 1) mapfunc[clone] -> SetBinContent(i + 1, 0.);
      //cout << "ratio : " << mapfunc[clone] -> GetBinContent(i + 1) << endl;
    }
    //mapfunc[clone] ->Draw("histp9same");
    
    //make correct error bars
    map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] = new TGraphAsymmErrors(mapfunc[clone]);
    
    for(unsigned int i = 0; i < N_bin; i++){
      double bkg_i = mapfunc[func + "rebin"] -> GetBinContent(i+1);
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
    
    map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] -> Draw("p0same");
  }
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
  latex_Lumi.DrawLatex(0.7, 0.96, "35.9 fb^{-1} (13 TeV)");

  cout << "9" << endl;

  
  TString pdfname;
  if(channel.Contains("EMu")){
    pdfname = "./plots_emu_method/EMu/";

  }
  else if(channel.Contains("DiEle")){
    pdfname = "./plots_emu_method/DiEle/";

  }
  else if(channel.Contains("DiMu")){
    pdfname = "./plots_emu_method/DiMu/";
  }
  else return;

  pdfname.Append(nameofhistogram);
  pdfname.Append("_DY_amcnlo.pdf");
  cout << "9.1" << endl;

  mapcanvas[canvas] -> SaveAs(pdfname);
  
  cout << "10" << endl;


}
//////////////////////////////////////////////////////////////////////////////


void makehistogram_signal_VS_bkg(TString nameofhistogram, double binx[], int N_bin, TString yaxis_name){
  
  TString title_y;
  title_y = "Events/bin";
  
  TString canvas = nameofhistogram;
  TString hstack = nameofhistogram;
  TString legend = nameofhistogram;
  canvas.Insert(0, "c_signal_VS");
  hstack.Insert(0, "hs_signal_VS");
  legend.Insert(0, "legend_");
    
  mapcanvas[canvas] = new TCanvas(canvas,"",800,800);
  gStyle -> SetOptStat(1111);
  mapcanvas[canvas] -> SetTopMargin( 0.05 );
  mapcanvas[canvas] -> SetBottomMargin( 0.13 );
  mapcanvas[canvas] -> SetRightMargin( 0.05 );
  mapcanvas[canvas] -> SetLeftMargin( 0.16 );
  mapcanvas[canvas] -> SetLogy();

  maphstack["hs_" + nameofhistogram] -> GetXaxis()->SetLabelSize(0.02);
  maphstack["hs_" + nameofhistogram] -> GetXaxis()->SetTitle(yaxis_name);
  maphstack["hs_" + nameofhistogram] -> GetXaxis()-> SetRangeUser(0., 5500.);
  
  maphstack["hs_" + nameofhistogram] -> Draw("hist");
  
  TString channel;
  if(nameofhistogram.Contains("DiMu")) channel = "MuMu";
  if(nameofhistogram.Contains("DiEle")) channel = "ElEl";

  int N_signal = 4;
  TString signal_points[] = {"Zp1000_HN200", "Zp2000_HN500", "Zp3000_HN1000", "Zp4000_HN1500"};
  TString signal_legend[] = {"Z' M=1000, N M=200", "Z' M=2000, N M=500", "Z' M=3000, N M=1000", "Z' M=4000, N M=1500"};
  
  Int_t colour_array[] = {623, 900, 432 , 600};
  
  TString overflow = "overflow";
  for(int i = 0; i < N_signal; i++){
    Int_t nx    = GetHist(nameofhistogram + channel + signal_points[i]) -> GetNbinsX()+1;
    Double_t x1 = GetHist(nameofhistogram + channel + signal_points[i]) -> GetBinLowEdge(1);
    Double_t bw = GetHist(nameofhistogram + channel + signal_points[i]) -> GetBinWidth(nx);
    Double_t x2 = GetHist(nameofhistogram + channel + signal_points[i]) -> GetBinLowEdge(nx)+bw;
    //maphist[nameofhistogram + channel + signal_points[i] + overflow] = new TH1F("", "", nx, x1, x2);

    TH1F *htmp = new TH1F("", "", nx, x1, x2);

    for (Int_t j = 1; j <= nx; j++) {
      htmp -> SetBinContent(j, GetHist(nameofhistogram + channel + signal_points[i]) -> GetBinContent(j));
      htmp -> SetBinError(j, GetHist(nameofhistogram + channel + signal_points[i]) -> GetBinError(j));
    }


    maphist[nameofhistogram + channel + signal_points[i] + "rebin"] = htmp -> Rebin(N_bin, nameofhistogram + channel + signal_points[i] + "rebin", binx);
    GetHist(nameofhistogram + channel + signal_points[i] + "rebin") -> SetLineColor(colour_array[i]);
    GetHist(nameofhistogram + channel + signal_points[i] + "rebin") -> SetLineWidth(5);

    GetHist(nameofhistogram + channel + signal_points[i] + "rebin") -> Draw("histsame");
    maplegend[legend] -> AddEntry(GetHist(nameofhistogram + channel + signal_points[i] + "rebin"), signal_legend[i], "l");
  }
  maplegend[legend] -> SetX1NDC(0.65);
  maplegend[legend] -> SetX2NDC(0.95);
  maplegend[legend] -> SetY1NDC(0.7);
  maplegend[legend] -> SetY2NDC(0.95);

  maplegend[legend] -> Draw("same");

  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_Lumi.SetTextSize(0.035);
  latex_Lumi.DrawLatex(0.7, 0.96, "35.9 fb^{-1} (13 TeV)");
  
  TString pdfname;
  pdfname = "./plots_emu_method/Signal_VS_bkg/";

  
  pdfname.Append(nameofhistogram);
  pdfname.Append("_DY_amcnlo.pdf");
  mapcanvas[canvas] -> SaveAs(pdfname);
  
}


//////////////////////////////////////////////////////////////////////////////

void draw_histogram_variable_bin(TString nameofhistogram, float xmin, float xmax, double binx[], int N_bin,  float ymax, TString name_x, bool name_y){

  TString regions[5] = {"CR1", "CR2", "CR3", "CR4", "SR1"};
  TString channels[3] = {"DiEle", "DiMu", "EMu"};
  int i_dir = 0;
  TString directories[15];
  for(int i = 0; i < 5; i++){
    for(int j = 0; j < 3; j++){
      directories[i_dir] = regions[i] + "_OS_" + channels[j];
      i_dir++;
    }
  }

  int N_directories = 15;
  
  //make emu method applied hists
  TString regions_OS[5] = {"CR1_OS", "CR2_OS", "CR3_OS", "CR4_OS", "SR1_OS"};
  TString channels_OS[2] = {"DiMu", "DiEle"};
  for(int i = 0; i < 5; i++){//for regions
    for(int j = 0; j < 2; j++){//for channels
      make_hist_after_emu_method(nameofhistogram, channels_OS[j], regions_OS[i]);
      //if(regions_OS[i].Contains("SR")) make_hist_after_emu_method(nameofhistogram, channels_OS[j] + "_HNcut", regions_OS[i]);
    }
  } 
  
  
  for(int i = 0; i < N_directories; i++){
    if(!directories[i].Contains("EMu"))makehistogram_variable_bin(nameofhistogram + "_" + directories[i], xmin, xmax, ymax, binx, N_bin, name_x, name_y, directories[i]);
    //if(directories[i].Contains("SR")) makehistogram_variable_bin(nameofhistogram + "_" + directories[i] + "_HNcut", xmin, xmax, ymax, binx, N_bin, name_x, name_y, directories[i] + "_HNcut");
  }

}


void draw_signal_VS_bkg(TString nameofhistogram, double binx[], int N_bin, TString yaxis_name){
  TString directories[2] = {"SR1_OS_DiEle", "SR1_OS_DiMu"};
  int N_directories = 2;
  
  for(int i = 0; i < N_directories; i++){
    makehistogram_signal_VS_bkg(nameofhistogram+ "_" + directories[i], binx, N_bin, yaxis_name);
    //makehistogram_signal_VS_bkg(nameofhistogram+ "_" + directories[i] + "_HNcut", binx, N_bin, yaxis_name);
  }
  
}

/// Main nction ////////////////////////////////////////////////////////////
void plot(){
  
  //open data and bkg
  openfile(Cycle_name, SingleMuon);
  openfile(Cycle_name, DoubleEG);
  openfile(Cycle_name, DY_low);
  openfile(Cycle_name, DY_high);
  //openfile(Cycle_name, Wjets);
  openfile(Cycle_name, SingleTbar_tW);
  openfile(Cycle_name, SingleTbar_t);
  openfile(Cycle_name, SingleTop_s);
  openfile(Cycle_name, SingleTop_tW);
  openfile(Cycle_name, SingleTop_t);
  openfile(Cycle_name, WW);
  openfile(Cycle_name, WZ);
  openfile(Cycle_name, ZZ);
  openfile(Cycle_name, WWW);
  openfile(Cycle_name, WWZ);
  openfile(Cycle_name, WZZ);
  openfile(Cycle_name, ZZZ);
  openfile(Cycle_name, WGtoLNuG);
  openfile(Cycle_name, ZGto2LG);
  openfile(Cycle_name, ttbar);
  

  //open signal samples
  
  int Zpmass[8] = {500, 750, 1000, 1500, 2000, 2500, 3000, 4000};
  for(int i = 0; i < 8; i++){
    int HNmass = 100;
    double ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
    while(ratio > 2){
      TString current_name = "Zp" + TString::Itoa(Zpmass[i], 10) + "_HN" + TString::Itoa(HNmass, 10);
      //cout << current_name << endl;                                                                                                                                                                                                                                           
      HNmass += 100;
      ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
      openfile_signal(Cycle_name, current_name, "MuMu");
      openfile_signal(Cycle_name, current_name, "ElEl");
    }
  }
  
  
  
  cout << "open files complete" << endl;

  //make bins for each variables
  double bin_llmass[32], bin_lljjjjmass[15], bin_pt[17], bin_2nd_lep_pt[15];
  int N_bin_llmass = 0, N_bin_lljjjjmass = 0, N_bin_pt = 0, N_bin_2nd_lep_pt = 0;
  for(int i = 0; i < 31; i++){
    N_bin_llmass++;
    bin_llmass[i] = 20. * (i + 0.);
  }
  bin_llmass[31] = 1000;
    
  for(int i = 0; i < 11; i++){
    N_bin_lljjjjmass++;
    bin_lljjjjmass[i] = 200. * (i + 0.);
  }
  bin_lljjjjmass[11] = 2500;
  N_bin_lljjjjmass += 3;
  bin_lljjjjmass[12] = 3500;
  bin_lljjjjmass[13] = 5000;
  bin_lljjjjmass[14] = 5500;//overflow bin
  

  for(int i = 0; i < 8; i++){
    N_bin_2nd_lep_pt++;
    bin_2nd_lep_pt[i] = 30. * (i + 0.);
  }
  bin_2nd_lep_pt[8] = 1000;

  for(int i = 0; i < 8; i++){
    N_bin_pt++;
    bin_pt[i] = 30. * (i + 0.);
  }
  bin_pt[8] = 270;
  bin_pt[9] = 330;
  bin_pt[10] = 390;
  bin_pt[11] = 1000;
  N_bin_pt += 3;

  draw_histogram_variable_bin("h_OS_lljjjjmass", 0., 5500., bin_lljjjjmass, N_bin_lljjjjmass, 100000., "m(Z') (GeV)", true);

  draw_histogram_variable_bin("h_OS_lljjjjmass_AK8_0", 0., 5500., bin_lljjjjmass, N_bin_lljjjjmass, 100000., "m(Z') (GeV)", true);
  draw_histogram_variable_bin("h_OS_lljjjjmass_AK8_1", 0., 5500., bin_lljjjjmass, N_bin_lljjjjmass, 100000., "m(Z') (GeV)", true);
  draw_histogram_variable_bin("h_OS_lljjjjmass_AK8_2", 0., 5500., bin_lljjjjmass, N_bin_lljjjjmass, 100000., "m(Z') (GeV)", true);
  draw_histogram_variable_bin("h_OS_leadingljjmass", 0., 5000., bin_lljjjjmass, N_bin_lljjjjmass, 100000., "m(1st N) (GeV)", true);
  draw_histogram_variable_bin("h_OS_secondljjmass", 0., 5000., bin_lljjjjmass, N_bin_lljjjjmass, 100000., "m(2nd N) (GeV)", true);
  draw_histogram_variable_bin("h_OS_llmass", 50., 1000., bin_llmass, N_bin_llmass, 10000., "m(ll) (GeV)", true);
  draw_histogram_variable_bin("h_OS_leadingLeptonPt", 0., 1000., bin_pt, N_bin_pt, 10000., "Pt(1st lepton) (GeV)", true);
  draw_histogram_variable_bin("h_OS_secondLeptonPt", 0., 1000., bin_2nd_lep_pt, N_bin_2nd_lep_pt, 10000.,"Pt(2nd lepton) (GeV)", true);
  draw_histogram_variable_bin("h_OS_leadingjet_pt", 0., 1000., bin_pt, N_bin_pt, 10000.,"Pt(1st jet) (GeV)", true);
  draw_histogram_variable_bin("h_OS_secondjet_pt", 0., 1000., bin_pt, N_bin_pt, 10000.,"Pt(2nd jet) (GeV)", true);
  
  

  
  draw_signal_VS_bkg("h_OS_lljjjjmass", bin_lljjjjmass, N_bin_lljjjjmass, "m(Z') (GeV)");
  draw_signal_VS_bkg("h_OS_leadingljjmass", bin_lljjjjmass, N_bin_lljjjjmass, "m(1st N) (GeV)");
  draw_signal_VS_bkg("h_OS_secondljjmass", bin_lljjjjmass, N_bin_lljjjjmass, "m(2nd N) (GeV)");
  

  /*
  // -- make file contains m(Z') plots, comment out draw_signal_VS_bkg function above
  TFile *MyFile = new TFile("Bkg_VS_signal_cutbased_e.root","RECREATE");
  gDirectory -> mkdir("MuMu");
  gDirectory -> Cd("MuMu");

  int N_bin = N_bin_lljjjjmass;
  double binx[15];
  for(int i = 0; i < 11; i++){
    binx[i] = 200. * (i + 0.);
  }
  binx[11] = 2500;
  binx[12] = 3500;
  binx[13] = 5000;
  binx[14] = 5500;//overflow bin    
  //cout << "0" << endl;
  
  maphstack["hs_h_OS_lljjjjmass_SR1_OS_DiMu"] -> Write();
  

  //cout << "0.1" << endl;
  
  for(int i = 0; i < 8; i++){
    int HNmass = 100;
    double ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
    TString channel = "MuMu";
    while(ratio > 2){
      TString current_name = "Zp" + TString::Itoa(Zpmass[i], 10) + "_HN" + TString::Itoa(HNmass, 10);
      HNmass += 100;
      ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
      //cout << "1" << endl;
      Int_t nx    = GetHist("h_OS_lljjjjmass_SR1_OS_DiMu" + channel + current_name) -> GetNbinsX()+1;
      Double_t x1 = GetHist("h_OS_lljjjjmass_SR1_OS_DiMu" + channel + current_name) -> GetBinLowEdge(1);
      Double_t bw = GetHist("h_OS_lljjjjmass_SR1_OS_DiMu" + channel + current_name) -> GetBinWidth(nx);
      Double_t x2 = GetHist("h_OS_lljjjjmass_SR1_OS_DiMu" + channel + current_name) -> GetBinLowEdge(nx)+bw;
      
      TH1F *htmp = new TH1F("", "", nx, x1, x2);
      //cout << "2" << endl;

      for (Int_t j = 1; j <= nx; j++) {
	htmp -> SetBinContent(j, GetHist("h_OS_lljjjjmass_SR1_OS_DiMu" + channel + current_name) -> GetBinContent(j));
	htmp -> SetBinError(j, GetHist("h_OS_lljjjjmass_SR1_OS_DiMu" + channel + current_name) -> GetBinError(j));
      }
      //cout << "3" << endl;

      maphist["h_OS_lljjjjmass_SR1_OS_DiMu" + channel + current_name + "rebin"] = htmp -> Rebin(N_bin, "h_OS_lljjjjmass_SR1_OS_DiMu" + channel + current_name + "rebin", binx);

      maphist["h_OS_lljjjjmass_SR1_OS_DiMu" + channel + current_name + "rebin"] -> SetName("h_OS_lljjjjmass_SR1_DiMu_" + current_name);
      maphist["h_OS_lljjjjmass_SR1_OS_DiMu" + channel + current_name + "rebin"] -> Write();
    }
  }

  for(int i = 0; i < 8; i++){
    TString channel = "MuMu";
    int HNmass = 100;
    double ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
    while(ratio > 2){
      TString current_name = "Zp" + TString::Itoa(Zpmass[i], 10) + "_HN" + TString::Itoa(HNmass, 10);
      cout << current_name << endl;
      HNmass += 100;
      ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
      maphist["denom" + channel + current_name] -> SetName("Den_" + current_name + "_h_OS_lljjjjmass_SR1_DiMu");
      maphist["denom" + channel + current_name] -> Write();
    }
  }
  
  gDirectory -> Cd("../");
  gDirectory -> mkdir("MuMu_HNcut");
  gDirectory -> Cd("MuMu_HNcut");
  maphstack["hs_h_OS_lljjjjmass_SR1_OS_DiMu_HNcut"] -> Write();
  
  for(int i = 0; i < 8; i++){
    int HNmass = 100;
    double ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
    TString channel = "MuMu";
    while(ratio > 2){
      TString current_name = "Zp" + TString::Itoa(Zpmass[i], 10) + "_HN" + TString::Itoa(HNmass, 10);
      HNmass += 100;
      ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
      //cout << "1" << endl;                                                                                                                                                                                                                                                    
      Int_t nx    = GetHist("h_OS_lljjjjmass_SR1_OS_DiMu_HNcut" + channel + current_name) -> GetNbinsX()+1;
      Double_t x1 = GetHist("h_OS_lljjjjmass_SR1_OS_DiMu_HNcut" + channel + current_name) -> GetBinLowEdge(1);
      Double_t bw = GetHist("h_OS_lljjjjmass_SR1_OS_DiMu_HNcut" + channel + current_name) -> GetBinWidth(nx);
      Double_t x2 = GetHist("h_OS_lljjjjmass_SR1_OS_DiMu_HNcut" + channel + current_name) -> GetBinLowEdge(nx)+bw;

      TH1F *htmp = new TH1F("", "", nx, x1, x2);
      //cout << "2" << endl;                                                                                                                                                                                                                                                    

      for (Int_t j = 1; j <= nx; j++) {
        htmp -> SetBinContent(j, GetHist("h_OS_lljjjjmass_SR1_OS_DiMu_HNcut" + channel + current_name) -> GetBinContent(j));
        htmp -> SetBinError(j, GetHist("h_OS_lljjjjmass_SR1_OS_DiMu_HNcut" + channel + current_name) -> GetBinError(j));
      }
      //cout << "3" << endl;                                                                                                                                                                                                                                                    

      maphist["h_OS_lljjjjmass_SR1_OS_DiMu_HNcut" + channel + current_name + "rebin"] = htmp -> Rebin(N_bin, "h_OS_lljjjjmass_SR1_OS_DiMu_HNcut" + channel + current_name + "rebin", binx);

      maphist["h_OS_lljjjjmass_SR1_OS_DiMu_HNcut" + channel + current_name + "rebin"] -> SetName("h_OS_lljjjjmass_SR1_DiMu_" + current_name + "_HNcut");
      maphist["h_OS_lljjjjmass_SR1_OS_DiMu_HNcut" + channel + current_name + "rebin"] -> Write();
    }
  }

  for(int i = 0; i < 8; i++){
    TString channel = "MuMu";
    int HNmass = 100;
    double ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
    while(ratio > 2){
      TString current_name = "Zp" + TString::Itoa(Zpmass[i], 10) + "_HN" + TString::Itoa(HNmass, 10);
      cout << current_name << endl;
      HNmass += 100;
      ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
      maphist["denom" + channel + current_name] -> SetName("Den_" + current_name + "_h_OS_lljjjjmass_SR1_DiMu_HNcut");
      maphist["denom" + channel + current_name] -> Write();
    }
  }
  // -- 
  
  gDirectory -> Cd("../");
  gDirectory -> mkdir("ElEl");
  gDirectory -> Cd("ElEl");
  maphstack["hs_h_OS_lljjjjmass_SR1_OS_DiEle"] -> Write();
  
  
  for(int i = 0; i < 8; i++){
    int HNmass = 100;
    double ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
    TString channel = "ElEl";
    while(ratio > 2){
      TString current_name = "Zp" + TString::Itoa(Zpmass[i], 10) + "_HN" + TString::Itoa(HNmass, 10);
      HNmass += 100;
      ratio = (Zpmass[i] + 0.) / (HNmass + 0.);

      Int_t nx    = GetHist("h_OS_lljjjjmass_SR1_OS_DiEle" + channel + current_name) -> GetNbinsX()+1;
      Double_t x1 = GetHist("h_OS_lljjjjmass_SR1_OS_DiEle" + channel + current_name) -> GetBinLowEdge(1);
      Double_t bw = GetHist("h_OS_lljjjjmass_SR1_OS_DiEle" + channel + current_name) -> GetBinWidth(nx);
      Double_t x2 = GetHist("h_OS_lljjjjmass_SR1_OS_DiEle" + channel + current_name) -> GetBinLowEdge(nx)+bw;

      TH1F *htmp = new TH1F("", "", nx, x1, x2);

      for (Int_t j = 1; j <= nx; j++) {
        htmp -> SetBinContent(j, GetHist("h_OS_lljjjjmass_SR1_OS_DiEle" + channel + current_name) -> GetBinContent(j));
        htmp -> SetBinError(j, GetHist("h_OS_lljjjjmass_SR1_OS_DiEle" + channel + current_name) -> GetBinError(j));
      }

      maphist["h_OS_lljjjjmass_SR1_OS_DiEle" + channel + current_name + "rebin"] = htmp -> Rebin(N_bin, "h_OS_lljjjjmass_SR1_OS_DiEle" + channel + current_name + "rebin", binx);

      maphist["h_OS_lljjjjmass_SR1_OS_DiEle" + channel + current_name + "rebin"] -> SetName("h_OS_lljjjjmass_SR1_DiEle_" + current_name);
      maphist["h_OS_lljjjjmass_SR1_OS_DiEle" + channel + current_name + "rebin"] -> Write();
    }
  }

  for(int i = 0; i < 8; i++){
    TString channel = "ElEl";
    int HNmass = 100;
    double ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
    while(ratio > 2){
      TString current_name = "Zp" + TString::Itoa(Zpmass[i], 10) + "_HN" + TString::Itoa(HNmass, 10);
      cout << current_name << endl;
      HNmass += 100;
      ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
      maphist["denom" + channel + current_name] -> SetName("Den_" + current_name + "_h_OS_lljjjjmass_SR1_DiMu");
      maphist["denom" + channel + current_name] -> Write();
    }
  }

  gDirectory -> Cd("../");
  gDirectory -> mkdir("ElEl_HNcut");
  gDirectory -> Cd("ElEl_HNcut");
  maphstack["hs_h_OS_lljjjjmass_SR1_OS_DiEle_HNcut"] -> Write();


  for(int i = 0; i < 8; i++){
    int HNmass = 100;
    double ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
    TString channel = "ElEl";
    while(ratio > 2){
      TString current_name = "Zp" + TString::Itoa(Zpmass[i], 10) + "_HN" + TString::Itoa(HNmass, 10);
      HNmass += 100;
      ratio = (Zpmass[i] + 0.) / (HNmass + 0.);

      Int_t nx    = GetHist("h_OS_lljjjjmass_SR1_OS_DiEle_HNcut" + channel + current_name) -> GetNbinsX()+1;
      Double_t x1 = GetHist("h_OS_lljjjjmass_SR1_OS_DiEle_HNcut" + channel + current_name) -> GetBinLowEdge(1);
      Double_t bw = GetHist("h_OS_lljjjjmass_SR1_OS_DiEle_HNcut" + channel + current_name) -> GetBinWidth(nx);
      Double_t x2 = GetHist("h_OS_lljjjjmass_SR1_OS_DiEle_HNcut" + channel + current_name) -> GetBinLowEdge(nx)+bw;

      TH1F *htmp = new TH1F("", "", nx, x1, x2);

      for (Int_t j = 1; j <= nx; j++) {
        htmp -> SetBinContent(j, GetHist("h_OS_lljjjjmass_SR1_OS_DiEle_HNcut" + channel + current_name) -> GetBinContent(j));
        htmp -> SetBinError(j, GetHist("h_OS_lljjjjmass_SR1_OS_DiEle_HNcut" + channel + current_name) -> GetBinError(j));
      }

      maphist["h_OS_lljjjjmass_SR1_OS_DiEle_HNcut" + channel + current_name + "rebin"] = htmp -> Rebin(N_bin, "h_OS_lljjjjmass_SR1_OS_DiEle_HNcut" + channel + current_name + "rebin", binx);

      maphist["h_OS_lljjjjmass_SR1_OS_DiEle_HNcut" + channel + current_name + "rebin"] -> SetName("h_OS_lljjjjmass_SR1_DiEle_" + current_name + "_HNcut");
      maphist["h_OS_lljjjjmass_SR1_OS_DiEle_HNcut" + channel + current_name + "rebin"] -> Write();
    }
  }

  for(int i = 0; i < 8; i++){
    TString channel = "ElEl";
    int HNmass = 100;
    double ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
    while(ratio > 2){
      TString current_name = "Zp" + TString::Itoa(Zpmass[i], 10) + "_HN" + TString::Itoa(HNmass, 10);
      cout << current_name << endl;
      HNmass += 100;
      ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
      maphist["denom" + channel + current_name] -> SetName("Den_" + current_name + "_h_OS_lljjjjmass_SR1_DiMu_HNcut");
      maphist["denom" + channel + current_name] -> Write();
    }
  }
  








  
  MyFile -> Close();


  cout << "myfile closed" << endl;

  

  for(map<TString, TFile*>::iterator mapit = mapfile.begin(); mapit != mapfile.end(); mapit ++){
    mapit->second->Close();
  }
  
  */



}// End of Main Function ////////////////////////////////////////////////////// 
