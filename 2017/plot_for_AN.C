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

using namespace std;

map<TString, TH1*> maphist;
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
TString DY_low = "SKDYJets_10to50";
TString DY_high = "DYJets";

TString Wjets = "SKWJets";

//single top
TString SingleTbar_tW = "SKSingleTbar_tW";
TString SingleTbar_t = "SKSingleTbar_t";
TString SingleTop_s = "SKSingleTop_s";
TString SingleTop_tW = "SKSingleTop_tW";
TString SingleTop_t = "SKSingleTop_t";

TString WGtoLNuG = "SKWGtoLNuG";

//VV
TString WW = "WWTo2L2Nu_powheg";
TString WZ_2L = "WZTo2L2Q";
TString WZ_3L = "WZTo3LNu";
TString ZZ_2L = "ZZTo2L2Q";
TString ZZ_4L = "ZZTo4L_powheg";

//VVV
TString WWW = "SKWWW";
TString WWZ = "SKWWZ";
TString WZZ = "SKWZZ";
TString ZZZ = "SKZZZ";

//TString ttbar = "SKTTJets_aMC";
TString ttbar = "TTLL_powheg";
//TString ttbar = "SKTTLL_powheg";

//emu ratio results
double emu_over_ee = 0.452;
double emu_over_mumu = 0.642;

bool debug = true;

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
  
  TString filename = cyclename + "_" + samplename + ".root";
  cout << "opening : " << filename << endl;
  
  mapfile[filename] = new TFile ((filename)) ;

  TIter next(gDirectory->GetListOfKeys());
  TKey *key;
  vector<TString> histnames;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    if(debug) cout << "[[openfile]] " << key -> GetName() << endl;
    histnames.push_back(key -> GetName());
  }
  
  for(unsigned int i = 0; i < histnames.size(); i ++){
    maphist[histnames.at(i) + cyclename + samplename] = (TH1*)gDirectory -> Get(histnames.at(i));
  }
  
  
  TString directories[4] = {"CR_Zmass_DiMu", "CR_Zmass_DiEle", "CR_ttbar_DiMu", "CR_ttbar_DiEle"};
  int N_directories = 4;
  
  if(debug) cout << "[[openfile]] go into directories" << endl;
  
  for(int i = 0; i < N_directories; i++){
    gDirectory->cd(directories[i]);
    TIter next_sub(gDirectory->GetListOfKeys());
    TKey *key_sub;
    vector<TString> histnames_sub;
    while ((key_sub = (TKey*)next_sub())) {
      TClass *cl = gROOT->GetClass(key_sub->GetClassName());
      if (!cl->InheritsFrom("TH1")) continue;
      histnames_sub.push_back(key_sub -> GetName());
    }

    for(unsigned int i = 0; i < histnames_sub.size(); i ++){
      maphist[histnames_sub.at(i) + cyclename + samplename] = (TH1*)gDirectory -> Get(histnames_sub.at(i));
    }
    gDirectory->cd("../");
  }
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

  TString directories[18] = {"CR1_DiEle", "CR1_DiMu", "CR1_EMu", "CR2_DiEle", "CR2_DiMu", "CR2_EMu", "CR3_DiEle", "CR3_DiMu", "CR3_EMu", "CR4_DiEle", "CR4_DiMu", "CR4_EMu", "CR5_DiEle", "CR5_DiMu", "CR5_EMu", "SR1_DiEle", "SR1_DiMu", "SR1_EMu"};
  int N_directories = 18;

  for(int i = 0; i < N_directories; i++){
    gDirectory->cd(directories[i]);
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


//////////////////////////////////////////////////////////////////////////////
void makehistogram_variable_bin(TString nameofhistogram, float xmin, float xmax, float ymax, double binx[], int N_bin, TString name_x, bool name_y, TString channel){
  
  cout <<"[[makehistogram_variable_bin]] Draw : " << nameofhistogram << endl;

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

  maplegend[legend] = new TLegend(0.69, 0.60, 0.96, 0.92);
  
  cout << "[[makehistogram_variable_bin]] 1" << endl;

  maphstack[hstack] = new THStack(hstack, "Stacked_" + nameofhistogram);
  gStyle->SetOptTitle(0);
    

  TString WW = "WWTo2L2Nu_powheg";
  TString WZ_2L = "WZTo2L2Q";
  TString WZ_3L = "WZTo3LNu";
  TString ZZ_2L = "ZZTo2L2Q";
  TString ZZ_4L = "ZZTo4L_powheg";

  int n_kind = 3;
  TString samples_array[] = {WW, ttbar, DY_high};
  Int_t colour_array[] = {419, 416, 400};
  TString samples_legend[] = {"Other backgrounds", "t#bar{t}", "Z/#gamma + jets"};

  cout << "[[makehistogram_variable_bin]] 2" << endl;
  
  TString name_cycle = nameofhistogram + Cycle_name;
  cout << "[[makehistogram_variable_bin]] check1" << endl;
  
  
  
  Int_t nx_data    = GetHist(nameofhistogram + Cycle_name + current_data) -> GetNbinsX();
  Double_t x1_data = GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinLowEdge(1);
  Double_t bw_data = GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinWidth(nx_data);
  Double_t x2_data = GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinLowEdge(nx_data) + bw_data;
  for(int i = 0; i < n_kind; i++){
    mapfunc[name_cycle + samples_array[i]] = new TH1F("", "", nx_data, x1_data, x2_data);
  }
  
  
  
  mapfunc[name_cycle + DY_high] -> Add(GetHist(name_cycle + DY_low));
  mapfunc[name_cycle + ttbar] -> Add(GetHist(name_cycle + ttbar));
  mapfunc[name_cycle + WW] -> Add(GetHist(name_cycle + WW));
  mapfunc[name_cycle + WW] -> Add(GetHist(name_cycle + WZ_2L));
  mapfunc[name_cycle + WW] -> Add(GetHist(name_cycle + WZ_3L));
  mapfunc[name_cycle + WW] -> Add(GetHist(name_cycle + ZZ_2L));
  mapfunc[name_cycle + WW] -> Add(GetHist(name_cycle + ZZ_4L));

  cout << "[[makehistogram_variable_bin]] 3" << endl;

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
  GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> SetMarkerSize(1.0);
  if(!channel.Contains("SR")){
    maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + current_data + "rebin"), "data", "lp");
  }
  double data_max = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetMaximum();


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
  
  if(!blind){
    GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> Draw("epsame");
  }
  GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetXaxis() -> SetRangeUser(xmin, xmax);
  GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> SetMinimum(1.0);

  mappad[pad1] -> Update();

  if(!blind){
    GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> Draw("epsame");
  }
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kBlack);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(1001);
  maplegend[legend] -> SetShadowColor(0); // 0 = transparent                                                                                                                                                                                                                    
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> Draw("same");

  //maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + data), "data", "p");                                                                                                                                                                                   

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
  mapfunc["stat" + nameofhistogram] -> SetMaximum(2.0);
  mapfunc["stat" + nameofhistogram] -> SetStats(0);
  
  mapfunc["stat" + nameofhistogram] -> Divide(mapfunc[func + "rebin"]);
  Int_t ncells = mapfunc["stat" + nameofhistogram] -> GetSize();
  for(int i = 0; i < ncells - 1; i++){
    mapfunc["stat" + nameofhistogram] -> SetBinContent(i, 1.);
  }
  mapfunc["stat" + nameofhistogram] -> Draw("CE2");
  mapfunc[clone] -> Divide(mapfunc[func + "rebin"]);
  if(!blind){
    mapfunc[clone] ->Draw("PE1same");
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
  latex_Lumi.DrawLatex(0.7, 0.96, "41.3 fb^{-1} (13 TeV)");

  TString pdfname;
  if(channel.Contains("EMu")){
    pdfname = "./plots/EMu/";

  }
  else if(channel.Contains("DiEle")){
    pdfname = "./plots/DiEle/";

  }
  else if(channel.Contains("DiMu")){
    pdfname = "./plots/DiMu/";
  }
  else return;

  pdfname.Append(nameofhistogram);
  pdfname.Append(".pdf");
  mapcanvas[canvas] -> SaveAs(pdfname);
  

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
  latex_Lumi.DrawLatex(0.7, 0.96, "41.3 fb^{-1} (13 TeV)");
  
  TString pdfname;
  pdfname = "./plots/Signal_VS_bkg/";

  
  pdfname.Append(nameofhistogram);
  pdfname.Append(".pdf");
  mapcanvas[canvas] -> SaveAs(pdfname);
  
}


//////////////////////////////////////////////////////////////////////////////


void draw_histogram_variable_bin(TString nameofhistogram, float xmin, float xmax, double binx[], int N_bin,  float ymax, TString name_x, bool name_y){

  TString directories[4] = {"CR_Zmass_DiMu", "CR_Zmass_DiEle", "CR_ttbar_DiMu", "CR_ttbar_DiEle"};
  int N_directories = 4;

  //void makehistogram_variable_bin(TString nameofhistogram, float xmin, float xmax, float ymax, double binx[], int N_bin, TString name_x, bool name_y, TString channel){
 
  for(int i = 0; i < N_directories; i++){
    makehistogram_variable_bin(nameofhistogram + "_" + directories[i], xmin, xmax, ymax, binx, N_bin, name_x, name_y, directories[i]);
  }
  
}


void draw_signal_VS_bkg(TString nameofhistogram, double binx[], int N_bin, TString yaxis_name){
  TString directories[2] = {"SR1_DiEle", "SR1_DiMu"};
  int N_directories = 2;

  for(int i = 0; i < N_directories; i++){
    makehistogram_signal_VS_bkg(nameofhistogram+ "_" + directories[i], binx, N_bin, yaxis_name);
  }
  
}

/// Main nction ////////////////////////////////////////////////////////////
void plot(){
  
  //open data and bkg
  openfile(Cycle_name, SingleMuon);
  openfile(Cycle_name, DoubleEG);
  openfile(Cycle_name, DY_high);
  openfile(Cycle_name, WW);
  openfile(Cycle_name, WZ_2L);
  openfile(Cycle_name, WZ_3L);
  openfile(Cycle_name, ZZ_2L);
  openfile(Cycle_name, ZZ_4L);
  openfile(Cycle_name, ttbar);

  //open signal samples
  /*
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
  */
  
  
  cout << "open files complete" << endl;
  
  //make bins for each variables
  double bin_llmass[42], bin_lljjjjmass[15], bin_pt[17];
  int N_bin_llmass = 0, N_bin_lljjjjmass = 0, N_bin_pt = 0;
  for(int i = 0; i < 41; i++){
    N_bin_llmass++;
    bin_llmass[i] = i + 70.;
  }
  bin_llmass[41] = 111.;
    
  for(int i = 0; i < 11; i++){
    N_bin_lljjjjmass++;
    bin_lljjjjmass[i] = 200. * (i + 0.);
  }
  bin_lljjjjmass[11] = 2500;
  N_bin_lljjjjmass += 3;
  bin_lljjjjmass[12] = 3500;
  bin_lljjjjmass[13] = 5000;
  bin_lljjjjmass[14] = 5500;//overflow bin
  
  

  for(int i = 0; i < 16; i++){
    N_bin_pt++;
    bin_pt[i] = 30. * (i + 0.);
  }
  bin_pt[16] = 1000;
  
  
  
  draw_histogram_variable_bin("mll", 70., 110., bin_llmass, N_bin_llmass, 10000., "m(ll) (GeV)", true);
  draw_histogram_variable_bin("mZp", 0., 5000., bin_lljjjjmass, N_bin_lljjjjmass, 10000., "m(lljjjj) (GeV)", true);
  


  
  /*
  // -- make file contains m(Z') plots, comment out draw_signal_VS_bkg function above
  TFile *MyFile = new TFile("Bkg_VS_signal.root","RECREATE");
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

  maphstack["hs_h_OS_lljjjjmass_SR1_DiMu"] -> Write();
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
      Int_t nx    = GetHist("h_OS_lljjjjmass_SR1_DiMu" + channel + current_name) -> GetNbinsX()+1;
      Double_t x1 = GetHist("h_OS_lljjjjmass_SR1_DiMu" + channel + current_name) -> GetBinLowEdge(1);
      Double_t bw = GetHist("h_OS_lljjjjmass_SR1_DiMu" + channel + current_name) -> GetBinWidth(nx);
      Double_t x2 = GetHist("h_OS_lljjjjmass_SR1_DiMu" + channel + current_name) -> GetBinLowEdge(nx)+bw;
      
      TH1F *htmp = new TH1F("", "", nx, x1, x2);
      //cout << "2" << endl;

      for (Int_t j = 1; j <= nx; j++) {
	htmp -> SetBinContent(j, GetHist("h_OS_lljjjjmass_SR1_DiMu" + channel + current_name) -> GetBinContent(j));
	htmp -> SetBinError(j, GetHist("h_OS_lljjjjmass_SR1_DiMu" + channel + current_name) -> GetBinError(j));
      }
      //cout << "3" << endl;

      maphist["h_OS_lljjjjmass_SR1_DiMu" + channel + current_name + "rebin"] = htmp -> Rebin(N_bin, "h_OS_lljjjjmass_SR1_DiMu" + channel + current_name + "rebin", binx);

      maphist["h_OS_lljjjjmass_SR1_DiMu" + channel + current_name + "rebin"] -> SetName("h_OS_lljjjjmass_SR1_DiMu_" + current_name);
      maphist["h_OS_lljjjjmass_SR1_DiMu" + channel + current_name + "rebin"] -> Write();
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
  gDirectory -> mkdir("ElEl");
  gDirectory -> Cd("ElEl");
  maphstack["hs_h_OS_lljjjjmass_SR1_DiEle"] -> Write();

  
  for(int i = 0; i < 8; i++){
    int HNmass = 100;
    double ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
    TString channel = "ElEl";
    while(ratio > 2){
      TString current_name = "Zp" + TString::Itoa(Zpmass[i], 10) + "_HN" + TString::Itoa(HNmass, 10);
      HNmass += 100;
      ratio = (Zpmass[i] + 0.) / (HNmass + 0.);

      Int_t nx    = GetHist("h_OS_lljjjjmass_SR1_DiEle" + channel + current_name) -> GetNbinsX()+1;
      Double_t x1 = GetHist("h_OS_lljjjjmass_SR1_DiEle" + channel + current_name) -> GetBinLowEdge(1);
      Double_t bw = GetHist("h_OS_lljjjjmass_SR1_DiEle" + channel + current_name) -> GetBinWidth(nx);
      Double_t x2 = GetHist("h_OS_lljjjjmass_SR1_DiEle" + channel + current_name) -> GetBinLowEdge(nx)+bw;

      TH1F *htmp = new TH1F("", "", nx, x1, x2);

      for (Int_t j = 1; j <= nx; j++) {
        htmp -> SetBinContent(j, GetHist("h_OS_lljjjjmass_SR1_DiEle" + channel + current_name) -> GetBinContent(j));
        htmp -> SetBinError(j, GetHist("h_OS_lljjjjmass_SR1_DiEle" + channel + current_name) -> GetBinError(j));
      }

      maphist["h_OS_lljjjjmass_SR1_DiEle" + channel + current_name + "rebin"] = htmp -> Rebin(N_bin, "h_OS_lljjjjmass_SR1_DiEle" + channel + current_name + "rebin", binx);

      maphist["h_OS_lljjjjmass_SR1_DiEle" + channel + current_name + "rebin"] -> SetName("h_OS_lljjjjmass_SR1_DiEle_" + current_name);
      maphist["h_OS_lljjjjmass_SR1_DiEle" + channel + current_name + "rebin"] -> Write();
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



  MyFile -> Close();






  for(map<TString, TFile*>::iterator mapit = mapfile.begin(); mapit != mapfile.end(); mapit ++){
    mapit->second->Close();
  }

  */






}// End of Main Function ////////////////////////////////////////////////////// 
