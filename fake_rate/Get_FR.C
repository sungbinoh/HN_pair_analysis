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
  
  vector<TString> plot_names = gethistnames("Electron_QCD_result.root");
  
  vector<TString> plots_names_L;
  vector<TString> plots_names_T;
  
  for(int i = 0; i < plot_names.size(); i++){
    TString current_name = plot_names.at(i);
    bool ID_syst = current_name.Contains("pt_vs_eta_L") ;
    bool ID_syst_else =current_name.Contains("pt_vs_eta_T") ;
    
    if(ID_syst) plots_names_L.push_back(current_name);
    else if(ID_syst_else) plots_names_T.push_back(current_name);
  }
  
  cout << "F 2D hists" << endl;
  for(int i = 0; i < plots_names_L.size(); i++){
    cout << plots_names_L.at(i) << endl;
    maphist[plots_names_L.at(i) + "Electron"] = GetHist_2D(plots_names_L.at(i) + "Electron") -> ProjectionX(plots_names_L.at(i) +"_Electron_1D", 0, 6, "");
    maphist[plots_names_L.at(i) + "Muon"] = GetHist_2D(plots_names_L.at(i) + "Muon") -> ProjectionX(plots_names_L.at(i) +"_Muon_1D", 0, 6, "");
    //maphist[plots_names_F.at(i)] -> Draw();
        
  }

  
  cout << "F0 2D hists" << endl;
  for(int i = 0; i < plots_names_T.size(); i++){
    cout << plots_names_T.at(i) << endl;
    maphist[plots_names_T.at(i) + "Electron"] = GetHist_2D(plots_names_T.at(i) + "Electron") -> ProjectionX(plots_names_T.at(i)+"_Electron_1D", 0, 6, "");
    maphist[plots_names_T.at(i) + "Muon"] = GetHist_2D(plots_names_T.at(i) + "Muon") -> ProjectionX(plots_names_T.at(i)+"_Muon_1D", 0, 6, "");
  }

  
  for(int i = 0; i < plots_names_L.size();i++){
    GetHist(plots_names_T.at(i) + "Electron") -> Divide(GetHist(plots_names_L.at(i) + "Electron"));
    GetHist(plots_names_T.at(i) + "Muon") -> Divide(GetHist(plots_names_L.at(i) + "Muon"));
    GetHist_2D(plots_names_T.at(i) + "Electron") -> Divide(GetHist_2D(plots_names_L.at(i) + "Electron"));
    GetHist_2D(plots_names_T.at(i) + "Muon") -> Divide(GetHist_2D(plots_names_L.at(i) + "Muon"));
  }
  
  
  
  TFile *MyFile = new TFile("FakeRate_QCD.root","RECREATE");
  for(int i = 0; i < plots_names_L.size(); i++){
    /*
    cout << "Writing : " << plots_names_F.at(i) + "data_fr" << endl;
    mapcanvas[plots_names_F.at(i)] = new TCanvas("", "", 800, 600);
    maphist_2D[plots_names_F.at(i)] -> Draw("colztext1e");
    gStyle -> SetPaintTextFormat("0.3f");
    mapcanvas[plots_names_F.at(i)] -> SaveAs("./pdfs/" + plots_names_F.at(i) + ".pdf");
    maphist_2D[plots_names_F.at(i)] -> Write();
    */


    cout << "Writing : " << plots_names_L.at(i) + "data_fr" << endl;
    //mapcanvas[plots_names_F.at(i)] = new TCanvas("", "", 800, 600);
    //maphist[plots_names_F.at(i)] -> Draw("AP");
    //mapcanvas[plots_names_F.at(i)] -> SaveAs("./pdfs/" + plots_names_F.at(i) + ".pdf");
    maphist[plots_names_T.at(i) + "Electron"] -> Write();
    maphist[plots_names_T.at(i) + "Muon"] -> Write();
    maphist_2D[plots_names_T.at(i) + "Electron"] -> SetName(plots_names_T.at(i) + "Electron");
    maphist_2D[plots_names_T.at(i) + "Muon"] -> SetName(plots_names_T.at(i) + "Muon");
    maphist_2D[plots_names_T.at(i) + "Electron"] -> Write();
    maphist_2D[plots_names_T.at(i) + "Muon"] -> Write();
  }
  
  MyFile->Close();

  /*
  draw_histogram("SingleMuonTrigger_Dijet_Awayjet_40_pt_F0", 0., 500., 1., 8000.); 
  draw_histogram("SingleMuonTrigger_Dijet_Awayjet_40_pt_F", 0., 500., 1., 500.);
  draw_histogram("SingleMuonTrigger_Dijet_Awayjet_40_eta_F", -3., 3., 1., 500.);
  draw_histogram("SingleMuonTrigger_Dijet_Awayjet_40_eta_F0", -3., 3., 1., 4000.);
  */
  //draw_histogram("h_tri_Z_tag_dimu_mass", 0., 500., 1., 1000.);
  //draw_histogram("h_tri_non_Z_tag_pt", 0., 500., 5., 1500.);
  

}// End of Main Function ////////////////////////////////////////////////////// 
