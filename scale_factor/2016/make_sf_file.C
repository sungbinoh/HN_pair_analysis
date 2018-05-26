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

void openfile(TString filename){
 
  cout << "opening : " << filename << endl;
  mapfile[filename] = new TFile ((filename)) ;
  
  //cout << samplename << "'s histograms are" << endl;                                                                                                                                                                                                                          
  //gDirectory->cd("Hists2D");

  TIter next(gDirectory->GetListOfKeys());
  TKey *key;
  vector<TString> histnames;
  vector<TString> histnames_2D;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (cl->InheritsFrom("TH2F")) histnames_2D.push_back(key -> GetName());
    if (!cl->InheritsFrom("TH1")) continue;
    histnames.push_back(key -> GetName());
  }

  for(int i = 0; i < histnames.size(); i ++){
    //cout << histnames.at(i) << endl;                                                                                                                                                                                                                                          
    maphist[histnames.at(i)] = (TH1*)gDirectory -> Get(histnames.at(i));
    //cout << histnames.at(i) + cyclename + samplename << endl;                                                                                                                                                                                                                 
  }
  for(int i = 0; i < histnames_2D.size(); i ++){
    //cout << histnames_2D.at(i) << endl;                                                                                                                                                                                                                                       
    maphist_2D[histnames_2D.at(i)] = (TH2*)gDirectory -> Get(histnames_2D.at(i));
  }

  gDirectory -> cd();

}

void openfile_muon(TString filename, TString name){

  cout << "opening : " << filename << endl;
  mapfile[filename] = new TFile ((filename)) ;

  TIter next(gDirectory->GetListOfKeys());
  TKey *key;
  vector<TString> histnames;
  vector<TString> histnames_2D;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (cl->InheritsFrom("TH2F")) histnames_2D.push_back(key -> GetName());
    if (!cl->InheritsFrom("TH1")) continue;
    histnames.push_back(key -> GetName());
  }

  for(int i = 0; i < histnames.size(); i ++){
    maphist[histnames.at(i) + name] = (TH1*)gDirectory -> Get(histnames.at(i));
  }
  for(int i = 0; i < histnames_2D.size(); i ++){
    maphist_2D[histnames_2D.at(i) + name] = (TH2*)gDirectory -> Get(histnames_2D.at(i));
  }

  gDirectory -> cd();
}


void plot(){

  openfile("scaleFactors.root");
  openfile_muon("TnP_NUM_MediumID_DENOM_generalTracks_VAR_map_pt_eta.root", "SB");
  openfile_muon("TnP_NUM_MiniIsoTight_DENOM_MediumID_VAR_map_pt_eta.root", "JS");
  openfile_muon("TnP_NUM_TightIP2D_DENOM_MediumID_VAR_map_pt_eta.root", "SH");

  TFile *MyFile = new TFile("SUSY_IDSF.root","RECREATE");
  /*
  maphist_2D["GsfElectronToMVAVLooseTightIP2D"] -> Multiply(GetHist_2D("MVAVLooseElectronToMini"));
  maphist_2D["GsfElectronToMVAVLooseTightIP2D"] -> Multiply(GetHist_2D("GsfElectronToMVATightTightIP2DSIP3D4"));
  maphist_2D["GsfElectronToMVAVLooseTightIP2D"] -> Multiply(GetHist_2D("MVAVLooseElectronToConvVetoIHit1"));
  maphist_2D["GsfElectronToMVAVLooseTightIP2D"] -> SetName("SF_electron");
  maphist_2D["GsfElectronToMVAVLooseTightIP2D"] -> Write();
  */

 
  maphist_2D["GsfElectronToMVAVLooseTightIP2D"] -> Multiply(GetHist_2D("GsfElectronToCutBasedSpring15L"));
  maphist_2D["GsfElectronToMVAVLooseTightIP2D"] -> Multiply(GetHist_2D("MVAVLooseElectronToMini"));
  maphist_2D["GsfElectronToMVAVLooseTightIP2D"] -> SetName("SF_electron");
  maphist_2D["GsfElectronToMVAVLooseTightIP2D"] -> Write();

  maphist_2D["GsfElectronToMVAVLooseTightIP2D"] -> Write();    
  maphist_2D["SFSB"] -> Multiply(GetHist_2D("SFJS"));
  maphist_2D["SFSB"] -> Multiply(GetHist_2D("SFSH"));
  maphist_2D["SFSB"] -> SetName("SF_muon");
  maphist_2D["SFSB"] -> Write();

  
  MyFile -> Close();
}
