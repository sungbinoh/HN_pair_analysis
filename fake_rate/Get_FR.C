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
void openfile(TString lepton){
  TString filename = lepton + "_QCD_result.root";
  
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
  
  for(int i = 0; i < histnames.size(); i ++){
    //cout << histnames.at(i) << endl;
    maphist[histnames.at(i) + lepton] = (TH1*)gDirectory -> Get(histnames.at(i));
    //cout << histnames.at(i) + cyclename + samplename << endl;
  }
  for(int i = 0; i < histnames_2D.size(); i ++){
    //cout << histnames_2D.at(i) << endl; 
    maphist_2D[histnames_2D.at(i) + lepton] = (TH2*)gDirectory -> Get(histnames_2D.at(i));
  }
  
  gDirectory -> cd();
  
}
//////////////////////////////////////////////////////////////////////////////


vector<TString> gethistnames(TString filename){

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



/// Main Function ////////////////////////////////////////////////////////////
void plot(){

  openfile("Electron");
  openfile("Muon");
  
  cout << "open files complete" << endl;
  
  vector<TString> plot_names_e = gethistnames("Electron_QCD_result.root");
  vector<TString> plot_names_mu = gethistnames("Muon_QCD_result.root");
    
  vector<TString> plots_names_L_e;
  vector<TString> plots_names_T_e;
  vector<TString> plots_names_L_mu;
  vector<TString> plots_names_T_mu;
  
  
  for(int i = 0; i < plot_names_e.size(); i++){
    TString current_name = plot_names_e.at(i);
    bool ID_syst = current_name.Contains("pt_vs_eta_L") ;
    bool ID_syst_else =current_name.Contains("pt_vs_eta_T") ;
    
    if(ID_syst) plots_names_L_e.push_back(current_name);
    else if(ID_syst_else) plots_names_T_e.push_back(current_name);
  }

  
  for(int i = 0; i < plot_names_mu.size(); i++){
    TString current_name = plot_names_mu.at(i);
    bool ID_syst = current_name.Contains("pt_vs_eta_L") ;
    bool ID_syst_else =current_name.Contains("pt_vs_eta_T") ;
    
    if(ID_syst) plots_names_L_mu.push_back(current_name);
    else if(ID_syst_else) plots_names_T_mu.push_back(current_name);
  }
  
  
  cout << "Make projections" << endl;
  for(int i = 0; i < plots_names_L_e.size(); i++){
    cout << plots_names_L_e.at(i) << endl;
    maphist[plots_names_L_e.at(i) + "Electron"] = GetHist_2D(plots_names_L_e.at(i) + "Electron") -> ProjectionX(plots_names_L_e.at(i) +"_Electron_1D", 0, 6, "");
  }
  for(int i = 0; i < plots_names_L_mu.size(); i++){
    cout << plots_names_L_mu.at(i) << endl;
    maphist[plots_names_L_mu.at(i) + "Muon"] = GetHist_2D(plots_names_L_mu.at(i) + "Muon") -> ProjectionX(plots_names_L_mu.at(i) +"_Muon_1D", 0, 6, "");
  }
  for(int i = 0; i < plots_names_T_e.size(); i++){
    cout << plots_names_T_e.at(i) << endl;
    maphist[plots_names_T_e.at(i) + "Electron"] = GetHist_2D(plots_names_T_e.at(i) + "Electron") -> ProjectionX(plots_names_T_e.at(i)+"_Electron_1D", 0, 6, "");
  }
  for(int i = 0; i < plots_names_T_mu.size(); i++){
    cout << plots_names_T_mu.at(i) << endl;
    maphist[plots_names_T_mu.at(i) + "Muon"] = GetHist_2D(plots_names_T_mu.at(i) + "Muon") -> ProjectionX(plots_names_T_mu.at(i)+"_Muon_1D", 0, 6, "");
  }
  
  cout << "Dividing histograms" << endl;
  cout << "e" << endl;
  for(int i = 0; i < plots_names_L_e.size();i++){
    cout << plots_names_L_e.at(i) << endl;
    cout << "1D" << endl;
    GetHist(plots_names_T_e.at(i) + "Electron") -> Divide(GetHist(plots_names_L_e.at(i) + "Electron"));
    cout << "2D" << endl;
    GetHist_2D(plots_names_T_e.at(i) + "Electron") -> Divide(GetHist_2D(plots_names_L_e.at(i) + "Electron"));
  }
  cout << "mu" << endl;
  for(int i = 0; i < plots_names_L_mu.size();i++){
    GetHist(plots_names_T_mu.at(i) + "Muon") -> Divide(GetHist(plots_names_L_mu.at(i) + "Muon"));
    GetHist_2D(plots_names_T_mu.at(i) + "Muon") -> Divide(GetHist_2D(plots_names_L_mu.at(i) + "Muon"));
  }
  
  
  cout << "Writing output file" << endl;
  TFile *MyFile = new TFile("FakeRate_QCD.root","RECREATE");
  for(int i = 0; i < plots_names_L_e.size(); i++){
    cout << "Writing : " << plots_names_L_e.at(i) + "data_fr" << endl;
    maphist[plots_names_T_e.at(i) + "Electron"] -> Write();
    maphist_2D[plots_names_T_e.at(i) + "Electron"] -> SetName(plots_names_T_e.at(i) + "Electron");
    maphist_2D[plots_names_T_e.at(i) + "Electron"] -> Write();
  }
  
  for(int i = 0; i < plots_names_L_mu.size(); i++){
    cout << "Writing : " << plots_names_L_mu.at(i) + "data_fr" << endl;
    maphist[plots_names_T_mu.at(i) + "Muon"] -> Write();
    maphist_2D[plots_names_T_mu.at(i) + "Muon"] -> SetName(plots_names_T_mu.at(i) + "Muon");
    maphist_2D[plots_names_T_mu.at(i) + "Muon"] -> Write();
  }
  
  MyFile->Close();
  
  
}// End of Main Function ////////////////////////////////////////////////////// 
