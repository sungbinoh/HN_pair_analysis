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

using namespace std;

map<TString, TH1*> maphist;
map<TString, TH2*> maphist_2D;
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
TString Cycle_name = "Electron_FR_cal_all";
//sample names                          
TString data = "data_SinglePhoton";

//DY
TString DY_low = "SKDY10to50";
TString DY_high = "SKDY50plus";

TString Wjets = "SKWJets";

//single top
TString SingleTbar_tW = "SKSingleTbar_tW";
TString SingleTbar_t = "SKSingleTbar_t";
TString SingleTop_s = "SKSingleTop_s";
TString SingleTop_tW = "SKSingleTop_tW";
TString SingleTop_t = "SKSingleTop_t";

TString WGtoLNuG = "SKWGtoLNuG";

//VV
TString WW_lnuqq = "SKWWToLNuQQ";
TString WW_2l2nu = "SKWWTo2L2Nu";
TString WZ = "SKWZTo3LNu_powheg";
TString ZZ = "SKZZTo4L_powheg";

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

TH2 * GetHist_2D(TString hname){

  TH2 *h = NULL;
  std::map<TString, TH2*>::iterator mapit = maphist_2D.find(hname);
  if(mapit != maphist_2D.end()) return mapit-> second;

  return h;

}

//////////////////////////////////////////////////////////////////////////////


/// Open ROOT file ///////////////////////////////////////////////////////////
void openfile(TString cyclename, TString samplename){
  
  TString underbar = "_";
  TString version = "cat_v8-0-8.root";
  TString filename = cyclename + underbar + samplename + underbar + version;
  if(samplename == data) filename = cyclename + underbar + samplename + "_cat_v8-0-8.root";
  
  cout << "opening : " << filename << endl;
  
  mapfile[filename] = new TFile ((filename)) ;
  
  //cout << samplename << "'s histograms are" << endl;
  gDirectory->cd("Hists2D");
  
  TIter next(gDirectory->GetListOfKeys());
  TKey *key;
  vector<TString> histnames;
  vector<TString> histnames_2D;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (cl->InheritsFrom("TH2D")) histnames_2D.push_back(key -> GetName());
    if (!cl->InheritsFrom("TH1")) continue;
    histnames.push_back(key -> GetName());
  }
  
  for(unsigned int i = 0; i < histnames.size(); i ++){
    //cout << histnames.at(i) << endl;
    maphist[histnames.at(i) + cyclename + samplename] = (TH1*)gDirectory -> Get(histnames.at(i));
    //cout << histnames.at(i) + cyclename + samplename << endl;
  }
  for(unsigned int i = 0; i < histnames_2D.size(); i ++){
    //cout << histnames_2D.at(i) << endl; 
    maphist_2D[histnames_2D.at(i) + cyclename + samplename] = (TH2*)gDirectory -> Get(histnames_2D.at(i));
  }
  
  gDirectory -> cd();
  
}
//////////////////////////////////////////////////////////////////////////////


vector<TString> gethistnames(TString cyclename, TString samplename){

  TString underbar = "_";
  TString version = "cat_v8-0-8.root";
  TString filename = cyclename + underbar + samplename + underbar + version;
  if(samplename == data) filename = cyclename + underbar + samplename + "_cat_v8-0-8.root";

  cout << "opening : " << filename << endl;

  mapfile[filename] = new TFile ((filename)) ;
  gDirectory->cd("Hists2D");

  TIter next(gDirectory->GetListOfKeys());
  TKey *key;
  vector<TString> histnames;
  vector<TString> histnames_2D;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (cl->InheritsFrom("TH2D")) histnames_2D.push_back(key -> GetName());
    if (!cl->InheritsFrom("TH1")) continue;
    histnames.push_back(key -> GetName());
  }

  return histnames_2D;

  gDirectory -> cd();
}


/// Make Output Histogram ////////////////////////////////////////////////////
void makehistogram(TString nameofhistogram){
  
  TString name_cycle = nameofhistogram + Cycle_name;
  
  /*
  if(maphist[name_cycle + DY_high] || maphist[name_cycle + DY_low]){
    GetHist(name_cycle + DY_high) -> Add(GetHist(name_cycle + DY_low));
  }
  */
  if(maphist[name_cycle + SingleTop_s] || maphist[name_cycle + SingleTbar_tW] || maphist[name_cycle + SingleTbar_t] || maphist[name_cycle + SingleTop_tW] || maphist[name_cycle + SingleTop_t]){
    GetHist(name_cycle +  SingleTop_s) -> Add(GetHist(name_cycle + SingleTbar_tW));
    GetHist(name_cycle +  SingleTop_s) -> Add(GetHist(name_cycle + SingleTbar_t));
    GetHist(name_cycle +  SingleTop_s) -> Add(GetHist(name_cycle + SingleTop_tW));
    GetHist(name_cycle +  SingleTop_s) -> Add(GetHist(name_cycle + SingleTop_t));
  }
  
  cout << "1" << endl;
  
  maphist_2D[nameofhistogram + "MC_sum"] = (TH2*)GetHist_2D(name_cycle +  DY_high) -> Clone();
  cout << "2" << endl;
  GetHist_2D(nameofhistogram + "MC_sum") -> Add(GetHist_2D(name_cycle +  DY_low));
  cout << "3" << endl;

  //GetHist_2D(nameofhistogram + "MC_sum") -> Add(GetHist(name_cycle +  DY_high));
  GetHist_2D(nameofhistogram + "MC_sum") -> Add(GetHist_2D(name_cycle +  SingleTop_s));
  GetHist_2D(nameofhistogram + "MC_sum") -> Add(GetHist_2D(name_cycle +  ttbar));
  GetHist_2D(nameofhistogram + "MC_sum") -> Add(GetHist_2D(name_cycle +  Wjets));
  //GetHist_2D(nameofhistogram + "MC_sum") -> Add(GetHist_2D(name_cycle +  WGtoLNuG));
  //GetHist_2D(nameofhistogram + "MC_sum") -> Add(GetHist_2D(name_cycle +  ZGto2LG));
  GetHist_2D(nameofhistogram + "MC_sum") -> Add(GetHist_2D(name_cycle +  WW_lnuqq));
  GetHist_2D(nameofhistogram + "MC_sum") -> Add(GetHist_2D(name_cycle +  WW_2l2nu));
  GetHist_2D(nameofhistogram + "MC_sum") -> Add(GetHist_2D(name_cycle +  WZ));
  GetHist_2D(nameofhistogram + "MC_sum") -> Add(GetHist_2D(name_cycle +  ZZ));


  maphist_2D[nameofhistogram + "data_fr"] = (TH2*)GetHist_2D(name_cycle +  data) -> Clone();
  
  GetHist_2D(nameofhistogram + "data_fr") -> Add(GetHist_2D(nameofhistogram + "MC_sum"), -1);
  
}
//////////////////////////////////////////////////////////////////////////////

void draw_colz(TString histname){

  TString canvas = "c1_" + histname;
  
  mapcanvas[canvas] = new TCanvas(canvas,"",800,800);
  gStyle -> SetOptStat(0);

  mapcanvas[canvas] -> SetTopMargin( 0.05 );
  mapcanvas[canvas] -> SetBottomMargin( 0.13 );
  mapcanvas[canvas] -> SetRightMargin( 0.05 );
  mapcanvas[canvas] -> SetLeftMargin( 0.16 );
  
  gStyle->SetPaintTextFormat("4.4f");
  maphist_2D[histname] -> GetXaxis()->SetRangeUser(40., 150.);
  maphist_2D[histname] -> Draw("colztexte1");
  
  mapcanvas[canvas] -> SaveAs("./Fake_Rate/" + histname + ".pdf");
  
}



/// Main Function ////////////////////////////////////////////////////////////
void plot(){
    
  openfile(Cycle_name, data);
  openfile(Cycle_name, DY_low);
  openfile(Cycle_name, DY_high);
  openfile(Cycle_name, Wjets);
  openfile(Cycle_name, SingleTbar_tW);
  openfile(Cycle_name, SingleTbar_t);
  openfile(Cycle_name, SingleTop_s);
  openfile(Cycle_name, SingleTop_tW);
  openfile(Cycle_name, SingleTop_t);
  openfile(Cycle_name, WW_lnuqq);
  openfile(Cycle_name, WW_2l2nu);
  openfile(Cycle_name, WZ);
  openfile(Cycle_name, ZZ);
  openfile(Cycle_name, ttbar);

  /*
  openfile(Cycle_name, ZZTo4L);
  openfile(Cycle_name, WZTo3LNu);
  openfile(Cycle_name, ZGto2LG);
  openfile(Cycle_name, fake);
  */
  
  cout << "open files complete" << endl;
  
  vector<TString> plot_names = gethistnames(Cycle_name, data);
  
  vector<TString> plots_names_F;
  vector<TString> plots_names_F0;
  for(unsigned int i = 0; i < plot_names.size(); i++){
    TString current_name = plot_names.at(i);
    bool ID_syst = current_name.Contains("pt_cone_vs_eta_F0");
    bool ID_syst_else =current_name.Contains("pt_cone_vs_eta");
    //bool Away_jet_syst = current_name.Contains("pt_cone_vs_eta_F0") && current_name.Contains("MUON_HN_LOOSEv7_SIP3");
    //bool Away_jet_syst_else = current_name.Contains("pt_cone_vs_eta") && current_name.Contains("MUON_HN_LOOSEv7_SIP3");
    
    if(ID_syst) plots_names_F0.push_back(current_name);
    else if(ID_syst_else) plots_names_F.push_back(current_name);
  }
  
  cout << "F 2D hists" << endl;
  for(unsigned int i = 0; i < plots_names_F.size(); i++){
    cout << plots_names_F.at(i) << endl;
  }
  cout << "F0 2D hists" << endl;
  for(unsigned int i = 0; i < plots_names_F0.size(); i++){
    cout << plots_names_F0.at(i) << endl;
  }

  
  for(unsigned int i = 0; i < plots_names_F.size(); i++){
    makehistogram(plots_names_F.at(i));
    makehistogram(plots_names_F0.at(i));
  }
  
  for(unsigned int i = 0; i < plots_names_F.size();i++){
    GetHist_2D(plots_names_F.at(i) + "data_fr") -> Divide(GetHist_2D(plots_names_F0.at(i) + "data_fr"));
  }
  
  
  TFile *MyFile = new TFile("Data_driven_FR_syst_DoubleMuon_Nvtx_Reweight.root","RECREATE");
  for(unsigned int i = 0; i < plots_names_F.size(); i++){
    cout << "Writing : " << plots_names_F.at(i) + "data_fr" << endl;;
    maphist_2D[plots_names_F.at(i) + "data_fr"] -> Write();
    //maphist_2D[plots_names_F0.at(i) + "data_fr"] -> Write();
    draw_colz(plots_names_F.at(i) + "data_fr");
  }
  
  MyFile->Close();
  
  
}// End of Main Function ////////////////////////////////////////////////////// 
