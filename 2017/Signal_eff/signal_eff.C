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
#include <stdlib.h>
#include <sstream>
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

  TString filename = "HN_pair_all_HNPairToJJJJ_" + channel + "_" + samplename + "_WR5000.root";
  cout << "[[openfile_signal]] Open " << filename << endl;
  TFile *current_file  = new TFile (("./" + channel + "/" + filename)) ;
  
  maphist["signal_eff" + samplename + channel] =  (TH1*)gDirectory -> Get("signal_eff");
  maphist["sumW" + samplename + channel] =  (TH1*)gDirectory -> Get("sumW");
  
  TString channel_2 = "";
  if(channel.Contains("MuMu")) channel_2 = "DiMu";
  else if(channel.Contains("EE")) channel_2 = "DiEle";
  
  double evt_Medium, evt_Tight, evt_HighPt;

  gDirectory -> cd("SR_" + channel_2);
  maphist["N_event_Medium" + samplename + channel] =  (TH1*)gDirectory -> Get("N_event");
  evt_Medium = maphist["N_event_Medium" + samplename + channel] -> GetBinContent(1);
  
  gDirectory -> cd("../SR_" + channel_2 + "_Tight_Mini");
  maphist["N_event_Tight" + samplename + channel] =  (TH1*)gDirectory -> Get("N_event");
  evt_Tight = maphist["N_event_Tight" + samplename + channel] -> GetBinContent(1);
  
  gDirectory -> cd("../SR_" + channel_2 + "_HighPt_Mini");
  maphist["N_event_HighPt" + samplename + channel] =  (TH1*)gDirectory -> Get("N_event");
  evt_HighPt = maphist["N_event_HighPt" + samplename + channel] -> GetBinContent(1);
  
  gDirectory -> cd("../");
  
  cout << "[[openfile_signal]] Call hist finished" << endl;

  // == Signle eff map
  // == 0.5 : MET filter, 1.5 : HLT_DoublePhoton70_v, 2.5 : HLT_Mu50_v, 3.5 : DiEle N lepton, 4.5 : DiMu N lepton , 5.5 : EMu N lepton
  // == 6.5 : DiEle lepton Pt cut, 7.5 : DiMu lepton Pt cut, 8.5 : EMu lepton Pt cut, 9.5 : DiEle Z CR, 10.5 : DiMu Z mass CR, 11.5 : EMu Z mass CR
  // == 12.5 : DiEle ttbar CR, 13.5 : DiMu ttbar CR, 14.5 : EMu ttbar CR, 15.5 : DiEle inv mll CR, 16.5 : DiMu inv mll CR, 17.5 : EMu inv mll CR
  // == 18.5 : DiEle SR, 19.5 : DiMu SR, 20.5 : EMu SR 
  mapfunc["1D_cutflow" + samplename + channel] = new TH1F("", "", 30, 0., 30.);
  mapfunc["1D_cutflow" + samplename + channel] -> SetDirectory(0);
  TH1::AddDirectory(kFALSE);
  
  double denom = maphist["sumW" + samplename + channel] -> GetBinContent(1) * 41527.5 / (57654 * 100);
  cout << "[[openfile_signal]] denom : " << denom << endl;
  cout << "[[openfile_signal]] evt_Medium : " << evt_Medium << ", evt_Tight : " << evt_Tight << ", evt_HighPt : " << evt_HighPt << endl;

  for(int i = 1; i < 22; i++){
    double current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(i) / denom;
    mapfunc["1D_cutflow" + samplename + channel] -> SetBinContent(i, current_eff);
  }
  
  cout << "[[openfile_signal]] Making 1D Signal Eff. finished " << endl;
  
  TObjArray *tx = samplename.Tokenize("_");
  TString mZp_str = ((TObjString *)(tx->At(0)))->String();
  TString mN_str = ((TObjString *)(tx->At(1)))->String();
  
  cout << "[[openfile_signal]] mZp_str : " << mZp_str << ", mN_str : " << mN_str << endl;
  
  mZp_str.Remove(0,2);
  mN_str.Remove(0,1);
  
  double mZp_double = mZp_str.Atof() + 10.;
  double mN_double = mN_str.Atof() + 10.;

  
  cout << "[[openfile_signal]] Filling 2D eff for each ID" << endl;

  
  maphist_2D["Medium_eff"] -> Fill(mZp_double, mN_double, evt_Medium / denom);
  maphist_2D["Tight_eff"] -> Fill(mZp_double, mN_double, evt_Tight / denom);
  maphist_2D["HighPt_eff"] -> Fill(mZp_double, mN_double, evt_HighPt / denom);

  current_file -> Close();
  delete current_file;
}


void draw_2D_hist(TString content){
  
  mapcanvas[ content ] = new TCanvas( content ,"",800,600);
  gStyle -> SetOptStat(0);
  mapcanvas[ content ] -> cd();
  
  gStyle -> SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3f");
  maphist_2D[ content ] -> Draw("COLZTEXT");
  
  
  
  mapcanvas[ content ] -> Update();
  mapcanvas[ content ] -> SaveAs("./pdfs/" + content + ".pdf");
}


void make_1D_hist(int Zpmass, TString channel){
  
  cout << "making 1D plot for Zp mass of " << Zpmass << " channel " << channel << endl;
  TString map_name;
  map_name = channel;
  mapcanvas[map_name] = new TCanvas(map_name,"",1500,800);
  gStyle -> SetOptStat(0);
  gStyle->SetPalette(55);//get new palete kDarkBodyRadiator
  mapcanvas[map_name] -> SetTopMargin( 0.05 );
  mapcanvas[map_name] -> SetBottomMargin( 0.13 );
  mapcanvas[map_name] -> SetRightMargin( 0.05 );
  mapcanvas[map_name] -> SetLeftMargin( 0.16 );
  mapcanvas[map_name] -> cd();

  
  maplegend[map_name] = new TLegend(0.80, 0.15, 0.95, 0.95);
  
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



/// Main Function ////////////////////////////////////////////////////////////
void plot(){
  
  maphist_2D["Medium_eff"] = new TH2D("", "", 12, 400., 5200., 12, 100., 2500.);
  maphist_2D["Tight_eff"] = new TH2D("", "", 12, 400., 5200., 12, 100., 2500.);
  maphist_2D["HighPt_eff"] = new TH2D("", "", 12, 400., 5200., 12, 100., 2500.);
  
  for(int i_Zp = 0; i_Zp < 12; i_Zp ++){
    int mZp = 400 + 400 * i_Zp; 
    int mN = 100;
    while(mN < mZp / 2){
      TString current_mZp = TString::Itoa(mZp,10);
      TString current_mN = TString::Itoa(mN,10);
      
      openfile_signal("ZP" + current_mZp + "_N" + current_mN, "MuMu");
      
      mN = mN + 200;
    } 
  }
  

  draw_2D_hist("Medium_eff");
  draw_2D_hist("Tight_eff");
  draw_2D_hist("HighPt_eff");

  
}// End of Main Function ////////////////////////////////////////////////////// 
