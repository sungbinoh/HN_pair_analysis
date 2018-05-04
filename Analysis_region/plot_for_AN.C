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
//#include <stdlib.h>
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
map<TString, TLine*> mapline;
map<TString, TKey*> maphistcheck;
map<TString, TList*> maplist;

//cycle name
TString Cycle_name = "HN_pair_all";
//sample names
//TString data = "data_DoubleMuon";
TString SingleMuon = "data_SingleMuon";
TString DoubleEG = "data_DoubleEG";
//TString periodC = "periodC_SKMuonEG";
//TString periodD = "periodD_SKMuonEG";

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

//TString ttbar = "SKTTJets_aMC";
TString ttbar = "SKTT_powheg";
TString ZGto2LG = "SKZGto2LG";


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
  if(samplename.Contains("data")) filename = cyclename + underbar + samplename + "_cat_v8-0-7.root";
  
  cout << "opening : " << filename << endl;
  
  mapfile[filename] = new TFile ((filename)) ;
  
  //cout << samplename << "'s histograms are" << endl;
  TString directories[12] = {"CR1_DiEle", "CR1_DiMu", "CR1_EMu", "CR2_DiEle", "CR2_DiMu", "CR2_EMu", "CR3_DiEle", "CR3_DiMu", "CR3_EMu", "SR1_DiEle", "SR1_DiMu", "SR1_EMu"}; 
  int N_directories = 12;
  
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

    for(int i = 0; i < histnames.size(); i ++){
      //cout << histnames.at(i) << endl;
      maphist[histnames.at(i) + cyclename + samplename] = (TH1*)gDirectory -> Get(histnames.at(i));
      //cout << histnames.at(i) + cyclename + samplename << endl;
    }
    gDirectory->cd("../");
  }
    
}
//////////////////////////////////////////////////////////////////////////////


/// Make Output Histogram ////////////////////////////////////////////////////
void makehistogram(TString nameofhistogram, float xmin, float xmax, float rebin, float ymax, TString name_x, bool name_y, TString channel){

  double additional_weight = 1.;
  if(channel.Contains("DiEle")) additional_weight = 20150.390 / 35820.052000;
  
  bool blind = false;
  blind = (channel.Contains("SR")) && (!channel.Contains("EMu"));
  
  TString current_data;
  if(channel.Contains("EMu") || channel.Contains("DiMu")) current_data = SingleMuon;
  if(channel.Contains("DiEle")) current_data = DoubleEG; 
  
  
  TString title_y;
  if(name_y) title_y = "Events/" + TString::Itoa(rebin,10) + " (GeV)";
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
  //gStyle -> SetPadTickY(1);
  //gStyle -> SetPadTickX(1);
  
  mappad[pad1] = new TPad("", "", 0, 0.25, 1, 1);
  mappad[pad1] -> SetTopMargin( 0.07 );
  mappad[pad1] -> SetBottomMargin( 0.02 );
  mappad[pad1] -> SetLeftMargin( 0.15 );
  mappad[pad1] -> SetRightMargin( 0.03 );
  //gStyle -> SetPadTickY(1);
  //gStyle -> SetPadTickX(1);
  mappad[pad1] -> Draw();
  mappad[pad1] -> cd();
  mappad[pad1] -> SetLogy();
  
  maplegend[legend] = new TLegend(0.69, 0.50, 0.96, 0.92);
  cout << "0" << endl;
  //cout << nameofhistogram + Cycle_name + data << endl;
  //mapfunc[func] = new TH1F("", "", GetHist(nameofhistogram + Cycle_name + periodC) -> GetNbinsX(),  GetHist(nameofhistogram + Cycle_name + periodC) -> GetXaxis() -> GetXmin(), GetHist(nameofhistogram + Cycle_name + periodC) -> GetXaxis() -> GetXmax());
  mapfunc[func] = new TH1F("", "", GetHist(nameofhistogram + Cycle_name + current_data) -> GetNbinsX(),  GetHist(nameofhistogram + Cycle_name + current_data) -> GetXaxis() -> GetXmin(), GetHist(nameofhistogram + Cycle_name + current_data) -> GetXaxis() -> GetXmax());
  
  cout << "1" << endl;
  
  maphstack[hstack] = new THStack(hstack, "Stacked_" + nameofhistogram);
  gStyle->SetOptTitle(0);
  //maphstack[hstack] -> SetTitle("CMS simulation");
  //gStyle->SetTitleFont(32,"t");
  //gStyle->SetTitleH(0.03);
  //gStyle->SetTitleX(0.83);
  //gStyle->SetTitleY(0.97);
  //gStyle->SetTitleH(0.1);//title perpendicular size  
  //gStyle->SetTitleW(0.9);//title width
  //gStyle->SetTitleX(.78);//title x-axis offset
  
  
  ///To do for each sample
  ///1. add, fill color, markerstyle
  ///2. add to stack
  ///3. set legend
  mapfunc[func] -> Rebin(rebin);
  
  int n_kind = 4;
  //TString samples_array[] = {DY_high, ttbar, Wjets, WGtoLNuG, WZ};
  TString samples_array[] = {WZ, ttbar, DY_high, Wjets};
  
  //kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416, kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
  //kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900
  //Int_t colour_array[] = {906, 416, 867, 600, 800, 432, 632};
  //Int_t colour_array[] = {400, 632, 416, 813, 819};
  Int_t colour_array[] = {819, 632, 400, 416};
  //TString samples_legend[] = {"DY", "top", "Wjets", "X#gamma", "VV"};
  TString samples_legend[] = {"VV", "top", "DY", "Wjets"};
  
  //std::vector<TString> samples(samples_array, samples_array + n_kind);
  //std::vector<Int_t> fill_colour(colour_array, colour_array + n_kind);
  //std::vector<TString> legends(samples_legend, samples_legend + n_kind);
  
  cout << "2" << endl;
  
  
  
  TString name_cycle = nameofhistogram + Cycle_name;
  cout << "check1" << endl;
  if(maphist[name_cycle + DY_high] || maphist[name_cycle + DY_low]){
    GetHist(name_cycle + DY_high) -> Add(GetHist(name_cycle + DY_low));
  }
  cout << "2.1" << endl;
  /*
  if(maphist[name_cycle + SingleTop_s] || maphist[name_cycle + SingleTbar_tW] || maphist[name_cycle + SingleTbar_t] || maphist[name_cycle + SingleTop_tW] || maphist[name_cycle + SingleTop_t]){
    cout << "2.1.1" << endl;
    
    GetHist(name_cycle +  SingleTop_s) -> Add(GetHist(name_cycle + SingleTbar_tW));
    cout << "2.1.2" << endl;
    GetHist(name_cycle +  SingleTop_s) -> Add(GetHist(name_cycle + SingleTbar_t));
    cout << "2.1.3" << endl;
    GetHist(name_cycle +  SingleTop_s) -> Add(GetHist(name_cycle + SingleTop_tW));
    cout << "2.1.4" << endl;
    GetHist(name_cycle +  SingleTop_s) -> Add(GetHist(name_cycle + SingleTop_t));
  }
  cout << "2.2" << endl;
  */
  /*
    if(maphist[name_cycle + WGtoLNuG] || maphist[name_cycle + ZGto2LG]){
    GetHist(name_cycle +  WGtoLNuG) -> Add(GetHist(name_cycle + ZGto2LG));
    }
  */
  if(maphist[name_cycle + WW] || maphist[name_cycle + WZ] || maphist[name_cycle + ZZ]){
    GetHist(name_cycle +  WZ) -> Add(GetHist(name_cycle + WW));
    GetHist(name_cycle +  WZ) -> Add(GetHist(name_cycle + ZZ));
  }
  cout << "2.3" << endl;

  //GetHist(name_cycle +  ttbar) -> Add(GetHist(name_cycle + SingleTop_s));
  
  cout << "3" << endl;
  
  
  cout << "check2" << endl;
  
  
  if(!channel.Contains("SR")){
    maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + current_data), "data", "p");
  }
  
  double data_max = GetHist(nameofhistogram + Cycle_name + current_data) -> GetMaximum();
  
  for(int i = 0; i < n_kind; i++){
    cout << samples_array[i] << endl;
    if(maphist[nameofhistogram + Cycle_name + samples_array[i]]){
      //cout << nameofhistogram + Cycle_name + samples_array[i] << endl;
      GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> SetFillColor(colour_array[i]);
      GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> Rebin(rebin);
      GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> Scale(additional_weight);
      cout << samples_array[i] << endl;
      maphstack[hstack] -> Add(GetHist(nameofhistogram + Cycle_name + samples_array[i]));
      maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + samples_array[n_kind - 1 - i]), samples_legend[n_kind - 1 - i], "f");
      mapfunc[func] -> Add(GetHist(nameofhistogram + Cycle_name + samples_array[i]));
    }
  }//for loop
  
  cout << "4" << endl;
  
  
  /*
    for(int i = 0; i < n_kind; i++){
    maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + samples_array[n_kind - i]), samples_legend[n_kind - i], "f");
    cout << "n_kind - i = " << n_kind - i << endl;
    }
  */
  
  
  ////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  maphstack[hstack] -> Draw("hist");
  maphstack[hstack] -> GetYaxis()->SetLabelSize(0.05);;
  maphstack[hstack] -> GetYaxis()->SetTitleSize(0.07);;
  maphstack[hstack] -> GetYaxis()->SetTitleOffset(1.02);;
  maphstack[hstack] -> GetXaxis()->SetLabelSize(0);;
  maphstack[hstack] -> GetXaxis() -> SetRangeUser(xmin, xmax);
  maphstack[hstack] -> SetMaximum(data_max * 100.);
  maphstack[hstack] -> SetMinimum(1.0);
  //maphstack[hstack] -> SetMinimum(0.001); 
  maphstack[hstack] -> GetXaxis() -> SetTitle(nameofhistogram);
  //maphstack[hstack] -> GetYaxis() -> SetTitle("Entries / 25 GeV");
  maphstack[hstack] -> GetYaxis() -> SetTitle(title_y);  
  maphstack[hstack] -> Draw("histsame");
  
  GetHist(nameofhistogram + Cycle_name + current_data) -> Rebin(rebin);
  GetHist(nameofhistogram + Cycle_name + current_data) -> SetMarkerStyle(20);
  GetHist(nameofhistogram + Cycle_name + current_data) -> SetMarkerSize(1.0);
  if(!blind){
    GetHist(nameofhistogram + Cycle_name + current_data) -> Draw("epsame");
  }
  GetHist(nameofhistogram + Cycle_name + current_data) -> GetXaxis() -> SetRangeUser(xmin, xmax);  
  GetHist(nameofhistogram + Cycle_name + current_data) -> SetMinimum(1.0);
  
  mappad[pad1] -> Update();
  
  if(!blind){
    GetHist(nameofhistogram + Cycle_name + current_data) -> Draw("epsame");       
  }
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kWhite);
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
  
  mapfunc[clone] = (TH1F*)GetHist(nameofhistogram + Cycle_name + current_data) -> Clone(clone);
  //mapfunc[clone] = (TH1F*)GetHist(nameofhistogram + Cycle_name + data) -> Clone(clone);
  
  mapfunc["stat" + nameofhistogram] = (TH1F*)GetHist(nameofhistogram + Cycle_name + current_data) -> Clone(clone);
  
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
  mapfunc["stat" + nameofhistogram] -> SetMinimum(0.5);
  mapfunc["stat" + nameofhistogram] -> SetMaximum(1.5);
  mapfunc["stat" + nameofhistogram] -> SetStats(0);
  //cout << "data integral , MC integral" << endl;
  //cout << mapfunc["stat" + nameofhistogram] -> Integral() << endl;
  //cout << mapfunc[func] -> Integral() << endl;
  
  mapfunc["stat" + nameofhistogram] -> Divide(mapfunc[func]);
  Int_t ncells = mapfunc["stat" + nameofhistogram] -> GetSize();
  for(int i = 0; i < ncells - 2; i++){
    mapfunc["stat" + nameofhistogram] -> SetBinContent(i, 1.);
  }
  mapfunc["stat" + nameofhistogram] -> Draw("CE2");
  mapfunc[clone] -> Divide(mapfunc[func]);
  if(!blind){
    mapfunc[clone] ->Draw("PE1same");
  }
  mapline[line] -> Draw("same");
  
  maplegend["bottom" + legend] = new TLegend(0.2, 0.8, 0.6, 0.9);
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

void draw_histogram(TString nameofhistogram, float xmin, float xmax, float rebin, float ymax, TString name_x, bool name_y){
  
  TString directories[12] = {"CR1_DiEle", "CR1_DiMu", "CR1_EMu", "CR2_DiEle", "CR2_DiMu", "CR2_EMu", "CR3_DiEle", "CR3_DiMu", "CR3_EMu", "SR1_DiEle", "SR1_DiMu", "SR1_EMu"};
  int N_directories = 12;

  //makehistogram(nameofhistogram + "_" + directories[0], xmin, xmax, rebin, ymax, name_x, name_y, directories[0]);

  
  for(int i = 0; i < N_directories; i++){
    makehistogram(nameofhistogram + "_" + directories[i], xmin, xmax, rebin, ymax, name_x, name_y, directories[i]);
  }
  
}


/// Main nction ////////////////////////////////////////////////////////////
void plot(){
  
  openfile(Cycle_name, SingleMuon);
  openfile(Cycle_name, DoubleEG);
  openfile(Cycle_name, DY_low);
  openfile(Cycle_name, DY_high);
  openfile(Cycle_name, Wjets);
  openfile(Cycle_name, SingleTbar_tW);
  openfile(Cycle_name, SingleTbar_t);
  openfile(Cycle_name, SingleTop_s);
  openfile(Cycle_name, SingleTop_tW);
  openfile(Cycle_name, SingleTop_t);
  openfile(Cycle_name, WW);
  openfile(Cycle_name, WZ);
  openfile(Cycle_name, ZZ);
  openfile(Cycle_name, ttbar);
  
  cout << "open files complete" << endl;

  draw_histogram("h_OS_lljjjjmass", 10., 5000., 50., 100000., "M(ll jjjj) (GeV)", true);   
  draw_histogram("h_OS_llmass", 10., 1000., 2., 10000., "M(ll) (GeV)", true);
  draw_histogram("h_OS_leadingLeptonPt", 0., 1000., 2., 10000., "Pt(1st lepton) (GeV)", true);
  draw_histogram("h_OS_secondLeptonPt", 0., 1000., 2.,10000.,"Pt(1st lepton) (GeV)", true);
  draw_histogram("h_OS_leadingjet_pt", 0., 1000., 2.,10000.,"Pt(1st jet) (GeV)", true);
  draw_histogram("h_OS_secondjet_pt", 0., 1000., 2.,10000.,"Pt(1st jet) (GeV)", true);
  



}// End of Main Function ////////////////////////////////////////////////////// 
