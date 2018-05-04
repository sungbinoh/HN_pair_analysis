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
map<TString, TH2D*> maphist_2D;
map<TString, TFile*> mapfile;
map<TString, TCanvas*> mapcanvas;
map<TString, TPad*> mappad;
map<TString, THStack*> maphstack;
map<TString, TLegend*> maplegend;
map<TString, TH1F*> mapfunc;
map<TString, TLine*> mapline;
map<TString, TKey*> maphistcheck;
map<TString, TList*> maplist;
map<TString, double> map_eff;
map<TString, double> map_eff_err;

//cycle name
TString Cycle_name = "HN_pair_all";
TString acceptancy = "signal_eff";


/// Getting Histogram Function ///////////////////////////////////////////////
TH1 * GetHist(TString hname){
  
  TH1 *h = NULL;
  std::map<TString, TH1*>::iterator mapit = maphist.find(hname);
  if(mapit != maphist.end()) return mapit-> second;
  
  return h;
  
}
//////////////////////////////////////////////////////////////////////////////


void openfile_signal(TString samplename){

  TString filename = "HN_pair_all_HNpair_MuMu_WR5000_" + samplename + "_official_cat_v8-0-7.root";
  cout << "opening : " << filename << endl;
  
  mapfile[filename] = new TFile ((filename)) ;
  gDirectory -> cd("Hists");
  
  maphist[acceptancy + samplename] =  (TH1*)gDirectory -> Get(acceptancy);

  TString content[27] = {"Empty", "after_skim", "MET_filter", "good_PV", "HLT_mu", "HLT_ele", "DiEle_lepton_N_cut", "DiMu_lepton_N_cut", "EMu_lepton_N_cut", "DiEle_lepton_Pt_cut", "DiMu_lepton_Pt_cut", "EMu_leptn_Pt_cut",
		       "DiEle_Mll_cut", "DiMu_Mll_cut", "EMu_Mll_cut", "DiEle_jet_N_cut", "DiMu_jet_N_cut", "EMu jet N cut", "DiEle_Mlljjjj_cut", "DiMu_Mlljjjj_cut", "EMu_Mlljjjj_cut", "DiEle_CR1", "DiMu_CR1", "EMu_CR1", "DiEle_CR2", "DiMu_CR2", "EMu CR2"};
  
  for(int i = 0; i < 27; i++){
    //cout << samplename + content[i] << endl;
    map_eff[samplename + content[i]] = maphist[acceptancy + samplename] -> GetBinContent(i);
    map_eff_err[samplename + content[i]] = maphist[acceptancy + samplename] -> GetBinError(i);

  }
  
  gDirectory -> cd();
}


void make_2D_hist(TString content){
  cout << "making 2D plot for " << content << endl;
  
  Float_t xbins[9] = {500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000};
  Float_t ybins[20];
  for(int i = 0; i < 20; i ++){
    ybins[i] = 100. + 100.0 * i;
  }

  maphist_2D[ content ] = new TH2D(content, content, 8, xbins, 19, ybins);
  gStyle->SetPaintTextFormat("4.4f");
  
  maphist_2D[ content ] -> SetMarkerSize(0.9);
  
  int Zpmass[8] = {500, 750, 1000, 1500, 2000, 2500, 3000, 4000};
  for(int i = 0; i < 8; i++){
    int HNmass = 100;
    double ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
    while(ratio > 2){
      TString current_name = "Zp" + TString::Itoa(Zpmass[i], 10) + "_HN" + TString::Itoa(HNmass, 10);
      //cout << current_name + content << endl;
      maphist_2D[ content ] -> Fill(Zpmass[i] + 1., HNmass + 1., map_eff[current_name + content] / map_eff[current_name + "after_skim"]);
      
      HNmass += 100;
      ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
    }
  }
  
  mapcanvas[ content ] = new TCanvas( content ,"",800,600);
  gStyle -> SetOptStat(0);
  mapcanvas[ content] -> cd();
  
  gStyle -> SetOptStat(0);
  maphist_2D[ content ] -> Draw("COLZTEXT");
  
  mapcanvas[ content ] -> Update();
  mapcanvas[ content ] -> SaveAs("./pdfs/" + content + "_eff_MiniIso.pdf");
}

void make_2D_hist_needed_Nevt(TString content, TString stat_error){
  
  //stat_error : 5_percent, 2p5_percent, 1_percent
  double N_stat;
  if(stat_error.Contains("5_percent")) N_stat = 400.0;
  else if(stat_error.Contains("2_percent")) N_stat = 1600.0;
  else if(stat_error.Contains("1_percent")) N_stat = 10000.0;
  else return;
  
  Float_t xbins[9] = {500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000};
  Float_t ybins[20];
  for(int i = 0; i < 20; i ++){
    ybins[i] = 100. + 100.0 * i;
  }

  maphist_2D[ content ] = new TH2D(content, content, 8, xbins, 19, ybins);
  gStyle->SetPaintTextFormat("4.0f");

  maphist_2D[ content ] -> SetMarkerSize(0.9);

  int Zpmass[8] = {500, 750, 1000, 1500, 2000, 2500, 3000, 4000};
  for(int i = 0; i < 8; i++){
    int HNmass = 100;
    double ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
    while(ratio > 2){
      TString current_name = "Zp" + TString::Itoa(Zpmass[i], 10) + "_HN" + TString::Itoa(HNmass, 10);
      //cout << current_name + content << endl;
      double accep_X_eff =  map_eff[current_name + content] / map_eff[current_name + "after_skim"];
      maphist_2D[ content ] -> Fill(Zpmass[i] + 1., HNmass + 1., N_stat / accep_X_eff);
      
      HNmass += 100;
      ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
    }
  }

  mapcanvas[ content ] = new TCanvas( content ,"",800,600);
  gStyle -> SetOptStat(0);
  mapcanvas[ content] -> cd();

  gStyle -> SetOptStat(0);
  maphist_2D[ content ] -> Draw("COLZTEXT");

  mapcanvas[ content ] -> Update();
  mapcanvas[ content ] -> SaveAs("./pdfs/" + stat_error + "_stat_error.pdf");


}

void make_1D_hist(int Zpmass){
  
  cout << "making 1D plot for Zp mass of " << Zpmass  << endl;
  TString content[11] = {"after_skim", "MET_filter", "good_PV", "HLT_mu", "DiMu_lepton_N_cut", "DiMu_lepton_Pt_cut", "DiMu_Mll_cut", "DiMu_jet_N_cut", "DiMu_Mlljjjj_cut", "DiMu_CR1", "DiMu_CR2"};
  
  int HNmass = 100;
  double ratio = (Zpmass + 0.) / (HNmass + 0.);
  int N_HN = 0;
  while(ratio > 2){
    N_HN++;
    TString current_name = "Zp" + TString::Itoa(Zpmass, 10) + "_HN" + TString::Itoa(HNmass, 10);
    HNmass += 100;
    ratio = (Zpmass + 0.) / (HNmass + 0.);
    maphist["eff" + current_name] = new TH1D("eff" + current_name, "eff" + current_name, 15, 0., 15.);
    for(int i = 0; i < 11; i++){
      maphist["eff" + current_name] -> SetBinContent(i + 2, map_eff[current_name + content[i]] / map_eff[current_name + "after_skim"]);
      //maphist["eff" + current_name] -> SetBinError(i + 2, map_eff_err[current_name + content[i]]);
    }
    for(int i = 0; i < 11; i++){
      maphist["eff" + current_name] -> GetXaxis() -> SetBinLabel(i + 2, content[i]);
    }
    maphist["eff" + current_name] -> SetTitle("");
    //maphist["eff" + current_name] -> Draw();
    gStyle -> SetOptStat(0);
  }
  
  TString map_name = "Zp" + TString::Itoa(Zpmass, 10) + "_eff";
  
  mapcanvas[map_name] = new TCanvas(map_name,"",1500,800);
  gStyle -> SetOptStat(0);
  gStyle->SetPalette(55);//get new palete kDarkBodyRadiator
  mapcanvas[map_name] -> SetTopMargin( 0.05 );
  mapcanvas[map_name] -> SetBottomMargin( 0.13 );
  mapcanvas[map_name] -> SetRightMargin( 0.05 );
  mapcanvas[map_name] -> SetLeftMargin( 0.16 );
  mapcanvas[map_name] -> cd();

  
  maplegend[map_name] = new TLegend(0.80, 0.15, 0.95, 0.95);
  
  HNmass = 100;
  ratio = (Zpmass + 0.) / (HNmass + 0.);
  int color = 20;
  while(ratio > 2){
    TString current_name = "Zp" + TString::Itoa(Zpmass, 10) + "_HN" + TString::Itoa(HNmass, 10);
    
    if(HNmass < 150){
      maphist["eff" + current_name] -> SetMinimum(0.);
      maphist["eff" + current_name] -> SetMaximum(1.2);
      maphist["eff" + current_name] -> SetLineColor(color);
      maphist["eff" + current_name] -> Draw("hist");
    }
    else{
      maphist["eff" + current_name] -> SetLineColor(color);
      maphist["eff" + current_name] -> Draw("histsame");
    }
    color += 76 / N_HN;
    
    maplegend[map_name] -> AddEntry(maphist["eff" + current_name], current_name, "fl");
    
    HNmass += 100;
    ratio = (Zpmass + 0.) / (HNmass + 0.);
  }
  
  //maphist[ map_name ] -> SetMarkerSize(0.9);

  maplegend[map_name] -> Draw("lcsame");

  mapcanvas[ map_name ] -> Update();
  mapcanvas[ map_name ] -> SaveAs("./pdfs/cutflow/cutflow_Zp_" + TString::Itoa(Zpmass, 10) + ".pdf");
}



/// Main Function ////////////////////////////////////////////////////////////
void plot(){
  
  int Zpmass[8] = {500, 750, 1000, 1500, 2000, 2500, 3000, 4000};
  for(int i = 0; i < 8; i++){
    int HNmass = 100;
    double ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
    while(ratio > 2){
      TString current_name = "Zp" + TString::Itoa(Zpmass[i], 10) + "_HN" + TString::Itoa(HNmass, 10);
      //cout << current_name << endl;
      HNmass += 100;
      ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
      openfile_signal(current_name);
    }
  
    make_1D_hist(Zpmass[i]);
  }

  cout << "open files complete" << endl;
  
  
  TString content[27] = {"Empty", "after_skim", "MET_filter", "good_PV", "HLT_mu", "HLT_ele", "DiEle_lepton_N_cut", "DiMu_lepton_N_cut", "EMu_lepton_N_cut", "DiEle_lepton_Pt_cut", "DiMu_lepton_Pt_cut", "EMu_leptn_Pt_cut",
			 "DiEle_Mll_cut", "DiMu_Mll_cut", "EMu_Mll_cut", "DiEle_jet_N_cut", "DiMu_jet_N_cut", "EMu jet N cut", "DiEle_Mlljjjj_cut", "DiMu_Mlljjjj_cut", "EMu_Mlljjjj_cut", "DiEle_CR1", "DiMu_CR1", "EMu_CR1", "DiEle_CR2", "DiMu_CR2", "EMu CR2"};
  /*
  for(int i = 1; i < 27; i++){
    make_2D_hist(content[i]);
  }
  */
  /*
  make_2D_hist_needed_Nevt("DiMu_Mlljjjj_cut", "5_percent");
  make_2D_hist_needed_Nevt("DiMu_Mlljjjj_cut", "2_percent");
  make_2D_hist_needed_Nevt("DiMu_Mlljjjj_cut", "1_percent");
  */
  //make_1D_hist(500);

}// End of Main Function ////////////////////////////////////////////////////// 
