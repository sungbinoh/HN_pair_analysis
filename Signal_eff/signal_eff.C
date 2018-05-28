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


void openfile_signal(TString samplename, TString channel){

  TString filename = "HN_pair_all_HNpair_" + channel + "_WR5000_" + samplename + "_official_cat_v8-0-7.root";
  cout << "opening : " << filename << endl;
  
  mapfile[filename] = new TFile (("./" + channel + "/" + filename)) ;
  gDirectory -> cd("Hists");
  
  maphist[acceptancy + samplename + channel] =  (TH1*)gDirectory -> Get(acceptancy);

  TString content[33] = {"Empty", "after_skim", "MET_filter", "good_PV", "HLT_mu", "HLT_ele", "DiEle_lepton_N_cut", "DiMu_lepton_N_cut", "EMu_lepton_N_cut", "DiEle_lepton_Pt_cut", "DiMu_lepton_Pt_cut", "EMu_leptn_Pt_cut",
			 "DiEle_Mll_cut", "DiMu_Mll_cut", "EMu_Mll_cut", "DiEle_jet_N_cut", "DiMu_jet_N_cut", "EMu jet N cut", "DiEle_Mlljjjj_cut", "DiMu_Mlljjjj_cut", "EMu_Mlljjjj_cut", "DiEle_CR1", "DiMu_CR1", "EMu_CR1", "DiEle_CR2", "DiMu_CR2", "EMu CR2",
			 "DiEle_CR3", "DiMu_CR3", "EMu CR3", "DiMu_Lep_pt_53_20", "DiMu_Lep_pt_53_53", "DiMu_Lep_pt_65_53"};
  
  for(int i = 0; i < 33; i++){
    //cout << samplename + content[i] << endl;
    map_eff[samplename + content[i] + channel] = maphist[acceptancy + samplename + channel] -> GetBinContent(i);
    map_eff_err[samplename + content[i] + channel] = maphist[acceptancy + samplename + channel] -> GetBinError(i);

  }
  
  gDirectory -> cd("../");

  TString directories[] = {"SR1_SS_DiMu", "SR1_SS_DiEle"};
  int n_directories = 2;
  for(int i = 0; i < n_directories; i++){
    gDirectory -> cd(directories[i]);
    TIter next(gDirectory->GetListOfKeys());
    TKey *key;
    vector<TString> histnames;
    while ((key = (TKey*)next())) {
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1")) continue;
      histnames.push_back(key -> GetName());
    }
    
    for(int i = 0; i < histnames.size(); i ++){
      //if(histnames.at(i).Contains("h_OS_lljjjjmass") ) cout << histnames.at(i) + samplename + channel << endl;
      maphist[histnames.at(i) + samplename + channel] = (TH1*)gDirectory -> Get(histnames.at(i));
    }
    gDirectory -> cd("../");
  }  
  
}


void make_2D_hist(TString content, TString channel){
  cout << "making 2D plot for " << content << endl;
  
  Float_t xbins[9] = {500, 750, 1000, 1500, 2000, 2500, 3000, 4000, 5000};
  Float_t ybins[20];
  for(int i = 0; i < 20; i ++){
    ybins[i] = 100. + 100.0 * i;
  }

  maphist_2D[ content + channel ] = new TH2D(content + channel, content + channel, 8, xbins, 19, ybins);
  gStyle->SetPaintTextFormat("4.4f");
  
  maphist_2D[ content + channel] -> SetMarkerSize(0.9);
  
  int Zpmass[8] = {500, 750, 1000, 1500, 2000, 2500, 3000, 4000};
  for(int i = 0; i < 8; i++){
    int HNmass = 100;
    double ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
    while(ratio > 2){
      TString current_name = "Zp" + TString::Itoa(Zpmass[i], 10) + "_HN" + TString::Itoa(HNmass, 10);
      //cout << current_name + content << endl;
      maphist_2D[ content + channel] -> Fill(Zpmass[i] + 1., HNmass + 1., map_eff[current_name + content + channel] / map_eff[current_name + "after_skim" + channel]);
      
      HNmass += 100;
      ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
    }
  }
  
  mapcanvas[ content + channel ] = new TCanvas( content + channel,"",800,600);
  gStyle -> SetOptStat(0);
  mapcanvas[ content + channel ] -> cd();
  
  gStyle -> SetOptStat(0);
  maphist_2D[ content + channel ] -> Draw("COLZTEXT");
  
  mapcanvas[ content + channel ] -> Update();
  mapcanvas[ content + channel ] -> SaveAs("./pdfs/" + channel + "_" + content + "_eff_MiniIso.pdf");
}

void make_2D_hist_needed_Nevt(TString content, TString stat_error, TString channel){
  
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

  maphist_2D[ content  + channel] = new TH2D(content + channel, content + channel, 8, xbins, 19, ybins);
  gStyle->SetPaintTextFormat("4.0f");

  maphist_2D[ content + channel ] -> SetMarkerSize(0.9);

  int Zpmass[8] = {500, 750, 1000, 1500, 2000, 2500, 3000, 4000};
  for(int i = 0; i < 8; i++){
    int HNmass = 100;
    double ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
    while(ratio > 2){
      TString current_name = "Zp" + TString::Itoa(Zpmass[i], 10) + "_HN" + TString::Itoa(HNmass, 10);
      //cout << current_name + content << endl;
      double accep_X_eff =  map_eff[current_name + content + channel] / map_eff[current_name + "after_skim" + channel];
      maphist_2D[ content + channel ] -> Fill(Zpmass[i] + 1., HNmass + 1., N_stat / accep_X_eff);
      
      HNmass += 100;
      ratio = (Zpmass[i] + 0.) / (HNmass + 0.);
    }
  }

  mapcanvas[ content + channel ] = new TCanvas( content + channel,"",800,600);
  gStyle -> SetOptStat(0);
  mapcanvas[ content + channel] -> cd();

  gStyle -> SetOptStat(0);
  maphist_2D[ content + channel ] -> Draw("COLZTEXT");

  mapcanvas[ content + channel ] -> Update();
  mapcanvas[ content + channel ] -> SaveAs("./pdfs/" + channel + "_" + stat_error + "_stat_error.pdf");


}

void make_1D_hist(int Zpmass, TString channel){
  
  cout << "making 1D plot for Zp mass of " << Zpmass << " channel " << channel << endl;
  
  TString content[12] = {"after_skim", "MET_filter", "good_PV", "HLT_mu", "DiMu_lepton_N_cut", "DiMu_lepton_Pt_cut", "DiMu_Mll_cut", "DiMu_jet_N_cut", "DiMu_Mlljjjj_cut", "DiMu_CR1", "DiMu_CR2", "DiMu_CR3"};
  if(channel.Contains("ElEl")){
    TString content_ee[12] = {"after_skim", "MET_filter", "good_PV", "HLT_ele", "DiEle_lepton_N_cut", "DiEle_lepton_Pt_cut", "DiEle_Mll_cut", "DiEle_jet_N_cut", "DiEle_Mlljjjj_cut", "DiEle_CR1", "DiEle_CR2", "DiEle_CR3"};
    for(int i = 0; i < 12; i++){
      content[i] = content_ee[i];
    }
  }
  
  int HNmass = 100;
  double ratio = (Zpmass + 0.) / (HNmass + 0.);
  int N_HN = 0;
  while(ratio > 2){
    N_HN++;
    TString current_name = "Zp" + TString::Itoa(Zpmass, 10) + "_HN" + TString::Itoa(HNmass, 10);
    HNmass += 100;
    ratio = (Zpmass + 0.) / (HNmass + 0.);
    maphist["eff" + current_name + channel] = new TH1D("eff" + current_name + channel, "eff" + current_name + channel, 15, 0., 15.);
    for(int i = 0; i < 12; i++){
      maphist["eff" + current_name + channel] -> SetBinContent(i + 2, map_eff[current_name + content[i] + channel] / map_eff[current_name + "after_skim" + channel]);
      //maphist["eff" + current_name] -> SetBinError(i + 2, map_eff_err[current_name + content[i]]);
    }
    for(int i = 0; i < 12; i++){
      maphist["eff" + current_name + channel] -> GetXaxis() -> SetBinLabel(i + 2, content[i]);
    }
    maphist["eff" + current_name + channel] -> SetTitle("");
    //maphist["eff" + current_name] -> Draw();
    gStyle -> SetOptStat(0);
  }
  
  TString map_name = "Zp" + TString::Itoa(Zpmass, 10) + "_eff" + channel;
  
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
      maphist["eff" + current_name + channel] -> SetMinimum(0.);
      maphist["eff" + current_name + channel] -> SetMaximum(1.2);
      maphist["eff" + current_name + channel] -> SetLineColor(color);
      maphist["eff" + current_name + channel] -> Draw("hist");
    }
    else{
      maphist["eff" + current_name + channel] -> SetLineColor(color);
      maphist["eff" + current_name + channel] -> Draw("histsame");
    }
    color += 76 / N_HN;
    
    maplegend[map_name] -> AddEntry(maphist["eff" + current_name + channel], current_name, "fl");
    
    HNmass += 100;
    ratio = (Zpmass + 0.) / (HNmass + 0.);
  }
  
  //maphist[ map_name ] -> SetMarkerSize(0.9);

  maplegend[map_name] -> Draw("lcsame");

  mapcanvas[ map_name ] -> Update();
  mapcanvas[ map_name ] -> SaveAs("./cutflow/cutflow_Zp_" + TString::Itoa(Zpmass, 10) + "_" + channel + "_doublephoton.pdf");
}


//draw signal shape
void draw_signal_shape(TString nameofhistogram, TString channel, double rebin){
 
  TString canvas = nameofhistogram;
  TString legend = nameofhistogram;
  canvas.Insert(0, "c_");
  legend.Insert(0, "legend_");
  
  mapcanvas[canvas] = new TCanvas(canvas,"",1000,800);
  gStyle -> SetOptStat(0);
  mapcanvas[canvas] -> SetTopMargin( 0.05 );
  mapcanvas[canvas] -> SetBottomMargin( 0.13 );
  mapcanvas[canvas] -> SetRightMargin( 0.05 );
  mapcanvas[canvas] -> SetLeftMargin( 0.16 );
  
  maplegend[legend] = new TLegend(0.65, 0.60, 0.95, 0.95);

  TString masspoints[] = {"Zp2000_HN500", "Zp3000_HN1000", "Zp4000_HN1500"};
  Int_t colour_array[] = {900, 432 , 600};
  int n_points = 3;
  cout << "0"<< endl;

  cout << nameofhistogram + "Zp1000_HN200" + channel << endl;
  maphist[nameofhistogram + "Zp1000_HN200" + channel] -> Rebin(rebin);
  maphist[nameofhistogram + "Zp1000_HN200" + channel] -> SetLineColor(623);
  maphist[nameofhistogram + "Zp1000_HN200" + channel] -> SetTitle("");
  maphist[nameofhistogram + "Zp1000_HN200" + channel] -> Draw("hist");
  maphist[nameofhistogram + "Zp1000_HN200" + channel] -> GetYaxis()->SetLabelSize(0.05);;
  maphist[nameofhistogram + "Zp1000_HN200" + channel] -> GetYaxis()->SetTitleSize(0.07);;
  maphist[nameofhistogram + "Zp1000_HN200" + channel] -> GetYaxis()->SetTitleOffset(1.02);;
  maphist[nameofhistogram + "Zp1000_HN200" + channel] -> SetMaximum(maphist[nameofhistogram + "Zp1000_HN200" + channel] -> GetMaximum() * 1.5);
  maphist[nameofhistogram + "Zp1000_HN200" + channel] -> GetXaxis()->SetLabelSize(0);
  maphist[nameofhistogram + "Zp1000_HN200" + channel] -> GetXaxis() -> SetRangeUser(0., 5000.);
  maphist[nameofhistogram + "Zp1000_HN200" + channel] -> GetXaxis()->SetLabelSize(0.02);
  maphist[nameofhistogram + "Zp1000_HN200" + channel] -> GetXaxis()->SetTitle("m(Z') (GeV)");
  maplegend[legend] -> AddEntry(GetHist(nameofhistogram + "Zp1000_HN200" + channel), "Zp1000_HN200, 0.1 pb", "l");
  cout << "1"<< endl;
  

  //maphist[nameofhistogram + "Zp1000_HN200" + channel] -> SetMaximum(data_max * 100.);
  maphist[nameofhistogram + "Zp1000_HN200" + channel] -> SetMinimum(1.0);
  
  for(int i = 0; i < n_points; i++){
    cout << "2" << endl;
    maphist[nameofhistogram + masspoints[i] + channel] -> Rebin(rebin);
    maphist[nameofhistogram + masspoints[i] + channel] -> SetLineColor(colour_array[i]);
    maphist[nameofhistogram + masspoints[i] + channel] -> Draw("histsame");
    maplegend[legend] -> AddEntry(GetHist(nameofhistogram + masspoints[i] + channel), masspoints[i] + ", 0.1 pb", "l");
  }

  cout << "3" << endl;
  
  maplegend[legend] -> Draw("lcsame");

  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Simulation}}}");
  latex_Lumi.SetTextSize(0.035);
  latex_Lumi.DrawLatex(0.7, 0.96, "35.9 fb^{-1} (13 TeV)");

  
  mapcanvas[ canvas ] -> Update();
  mapcanvas[ canvas ] -> SaveAs("./signal_shape/signal_shape_Zp_" + channel + ".pdf");

}


void draw_signal_shape_normalized(TString nameofhistogram, TString channel, double rebin){

  TString canvas = nameofhistogram;
  TString legend = nameofhistogram;
  canvas.Insert(0, "c_");
  legend.Insert(0, "legend_");

  mapcanvas[canvas] = new TCanvas(canvas,"",1000,800);
  gStyle -> SetOptStat(0);
  mapcanvas[canvas] -> SetTopMargin( 0.05 );
  mapcanvas[canvas] -> SetBottomMargin( 0.13 );
  mapcanvas[canvas] -> SetRightMargin( 0.05 );
  mapcanvas[canvas] -> SetLeftMargin( 0.16 );

  maplegend[legend] = new TLegend(0.16, 0.60, 0.46, 0.95);

  TString masspoints[] = {"Zp4000_HN500", "Zp4000_HN1000", "Zp4000_HN1500"};
  Int_t colour_array[] = {900, 432 , 600};
  int n_points = 3;
  cout << "0"<< endl;

  cout << nameofhistogram + "Zp4000_HN200" + channel << endl;
  maphist[nameofhistogram + "Zp4000_HN200" + channel] -> Rebin(rebin);
  double integ = maphist[nameofhistogram + "Zp4000_HN200" + channel] -> Integral();
  maphist[nameofhistogram + "Zp4000_HN200" + channel] -> Scale(1. / integ);
  maphist[nameofhistogram + "Zp4000_HN200" + channel] -> SetLineColor(623);
  maphist[nameofhistogram + "Zp4000_HN200" + channel] -> SetTitle("");
  maphist[nameofhistogram + "Zp4000_HN200" + channel] -> SetMaximum(0.12);
  maphist[nameofhistogram + "Zp4000_HN200" + channel] -> Draw("hist");
  maphist[nameofhistogram + "Zp4000_HN200" + channel] -> GetYaxis()->SetLabelSize(0.05);
  maphist[nameofhistogram + "Zp4000_HN200" + channel] -> GetYaxis()->SetTitleSize(0.07);
  maphist[nameofhistogram + "Zp4000_HN200" + channel] -> GetYaxis()->SetTitleOffset(1.02);
  maphist[nameofhistogram + "Zp4000_HN200" + channel] -> GetXaxis()->SetLabelSize(0);
  maphist[nameofhistogram + "Zp4000_HN200" + channel] -> GetXaxis() -> SetRangeUser(0., 5000.);
  maphist[nameofhistogram + "Zp4000_HN200" + channel] -> GetXaxis()->SetLabelSize(0.02);
  maphist[nameofhistogram + "Zp4000_HN200" + channel] -> GetXaxis()->SetTitle("m(Z') (GeV)");
  maplegend[legend] -> AddEntry(GetHist(nameofhistogram + "Zp4000_HN200" + channel), "Zp4000_HN200", "l");
  cout << "1"<< endl;


  //maphist[nameofhistogram + "Zp1000_HN200" + channel] -> SetMaximum(data_max * 100.);                                                                                                                                                                                         
  maphist[nameofhistogram + "Zp4000_HN200" + channel] -> SetMinimum(0.);

  for(int i = 0; i < n_points; i++){
    cout << "2" << endl;
    maphist[nameofhistogram + masspoints[i] + channel] -> Rebin(rebin);
    double integ_i =  maphist[nameofhistogram + masspoints[i] + channel] -> Integral();
    maphist[nameofhistogram + masspoints[i] + channel] -> Scale(1. / integ_i);
    maphist[nameofhistogram + masspoints[i] + channel] -> SetLineColor(colour_array[i]);
    maphist[nameofhistogram + masspoints[i] + channel] -> Draw("histsame");
    maplegend[legend] -> AddEntry(GetHist(nameofhistogram + masspoints[i] + channel), masspoints[i], "l");
  }

  cout << "3" << endl;

  maplegend[legend] -> Draw("lcsame");
  
  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Simulation}}}");

  mapcanvas[ canvas ] -> Update();
  mapcanvas[ canvas ] -> SaveAs("./signal_shape/Norm_signal_shape_Zp_" + channel + ".pdf");

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
      openfile_signal(current_name, "MuMu");
      openfile_signal(current_name, "ElEl");
    }
  
    make_1D_hist(Zpmass[i], "MuMu");
    make_1D_hist(Zpmass[i], "ElEl");
  }

  cout << "open files complete" << endl;
  
  
  TString content[27] = {"Empty", "after_skim", "MET_filter", "good_PV", "HLT_mu", "HLT_ele", "DiEle_lepton_N_cut", "DiMu_lepton_N_cut", "EMu_lepton_N_cut", "DiEle_lepton_Pt_cut", "DiMu_lepton_Pt_cut", "EMu_leptn_Pt_cut",
			 "DiEle_Mll_cut", "DiMu_Mll_cut", "EMu_Mll_cut", "DiEle_jet_N_cut", "DiMu_jet_N_cut", "EMu jet N cut", "DiEle_Mlljjjj_cut", "DiMu_Mlljjjj_cut", "EMu_Mlljjjj_cut", "DiEle_CR1", "DiMu_CR1", "EMu_CR1", "DiEle_CR2", "DiMu_CR2", "EMu CR2"};
  
  //draw_signal_shape("h_OS_lljjjjmass_SR1_DiEle", "ElEl", 50.);
  //draw_signal_shape("h_OS_lljjjjmass_SR1_DiMu", "MuMu", 50.);
  
  draw_signal_shape_normalized("h_SS_lljjjjmass_SR1_SS_DiEle", "ElEl", 50.);
  draw_signal_shape_normalized("h_SS_lljjjjmass_SR1_SS_DiMu", "MuMu", 50.);
  
  

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
