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
  
  double evt_all, evt_AK8_0, evt_AK8_1, evt_AK8_2;
  cout << "[[openfile_signal]] Go in to " << "SR_" << channel_2 << endl;
  gDirectory -> cd("SR_" + channel_2);
  maphist["N_event_all" + samplename + channel] =  (TH1*)gDirectory -> Get("mZp_SR_" + channel_2);
  if(maphist["N_event_all" + samplename + channel]){
    evt_all = maphist["N_event_all" + samplename + channel] -> Integral();
  }
  else evt_all = 0.;
  
  maphist["N_event_AK8_0" + samplename + channel] =  (TH1*)gDirectory -> Get("mZp_0AK8_SR_" + channel_2);
  if(maphist["N_event_AK8_0" + samplename + channel]){
    evt_AK8_0 = maphist["N_event_AK8_0" + samplename + channel] -> Integral();
  }
  else evt_AK8_0 = 0.;
  
  maphist["N_event_AK8_1" + samplename + channel] =  (TH1*)gDirectory -> Get("mZp_1AK8_SR_" + channel_2);
  if(maphist["N_event_AK8_1" + samplename + channel]){
    evt_AK8_1 = maphist["N_event_AK8_1" + samplename + channel] -> Integral();
  }
  else evt_AK8_1 = 0.;

  maphist["N_event_AK8_2" + samplename + channel] =  (TH1*)gDirectory -> Get("mZp_2AK8_SR_" + channel_2);
  if(maphist["N_event_AK8_2" + samplename + channel]){
    evt_AK8_2 = maphist["N_event_AK8_2" + samplename + channel] -> Integral();
  }
  else evt_AK8_2 = 0.;
  
  gDirectory -> cd("../");
  
  cout << "[[openfile_signal]] Call hist finished" << endl;

  // == Signle eff map
  // == 0.5 : MET filter, 1.5 : HLT_DoublePhoton70_v, 2.5 : HLT_Mu50_v, 3.5 : DiEle N lepton, 4.5 : DiMu N lepton , 5.5 : EMu N lepton
  // == 6.5 : DiEle lepton Pt cut, 7.5 : DiMu lepton Pt cut, 8.5 : EMu lepton Pt cut, 9.5 : DiEle Z CR, 10.5 : DiMu Z mass CR, 11.5 : EMu Z mass CR
  // == 12.5 : DiEle ttbar CR, 13.5 : DiMu ttbar CR, 14.5 : EMu ttbar CR, 15.5 : DiEle inv mll CR, 16.5 : DiMu inv mll CR, 17.5 : EMu inv mll CR
  // == 18.5 : DiEle SR, 19.5 : DiMu SR, 20.5 : EMu SR, 21.5 : DiEle charge, 22.5 : DiMu charge, 23.5 : EMu charge, 24.5 : DiEle mll, 25.5 : DiMu mll
  // == 26.5 : DiEle Njet, 27.5 : DiMu Njet, 28.5 : DiEle N mass, 29.5 : DiMu N mass
  mapfunc["1D_cutflow" + samplename + channel] = new TH1F("", "", 30, 0., 30.);
  mapfunc["1D_cutflow" + samplename + channel] -> SetDirectory(0);
  TH1::AddDirectory(kFALSE);
  
  //double denom = maphist["sumW" + samplename + channel] -> GetBinContent(1) * 41527.5 / (57654 * 100); // times 100 -> percent
  double denom = maphist["signal_eff" + samplename + channel] -> GetBinContent(1); // times 100 -> percent                                                                                                                                                                                                                                              

  cout << "[[openfile_signal]] denom : " << denom << endl;
  //cout << "[[openfile_signal]] evt_Medium : " << evt_Medium << ", evt_Tight : " << evt_Tight << ", evt_HighPt : " << evt_HighPt << endl;
  /*
  for(int i = 1; i < 22; i++){
    double current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(i) / denom;
    mapfunc["1D_cutflow" + samplename + channel] -> SetBinContent(i, current_eff);
  }
  */

  if(channel.Contains("MuMu")){
    
    double current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(1) / denom;//MET filter
    //cout << "MET filter passing evts : " << maphist["signal_eff" + samplename + channel] -> GetBinContent(2) << endl;
    mapfunc["1D_cutflow" + samplename + channel] -> Fill(0.5, current_eff);
    
    current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(3) / denom;//HLT
    mapfunc["1D_cutflow" + samplename + channel] -> Fill(1.5, current_eff);
    
    current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(5) / denom;//two veto mu
    mapfunc["1D_cutflow" + samplename + channel] -> Fill(2.5, current_eff);
    
    current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(23) / denom;//two veto mu charge
    mapfunc["1D_cutflow" + samplename + channel] -> Fill(3.5, current_eff);
    
    current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(8) / denom;//two veto mu pt cut 
    mapfunc["1D_cutflow" + samplename + channel] -> Fill(4.5, current_eff);

    current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(26) / denom;//mll
    mapfunc["1D_cutflow" + samplename + channel] -> Fill(5.5, current_eff);

    current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(28) / denom;//Njet
    mapfunc["1D_cutflow" + samplename + channel] -> Fill(6.5, current_eff);

    current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(30) / denom;//mN
    mapfunc["1D_cutflow" + samplename + channel] -> Fill(7.5, current_eff);

    current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(20) / denom;//SR
    mapfunc["1D_cutflow" + samplename + channel] -> Fill(8.5, current_eff);
    
  }
  
  if(channel.Contains("EE")){
    double current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(1) / denom;//MET 
    mapfunc["1D_cutflow" + samplename + channel] -> SetBinContent(1, current_eff);

    current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(2) / denom;//HLT
    mapfunc["1D_cutflow" + samplename + channel] -> SetBinContent(2, current_eff);

    current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(4) / denom;//two veto el 
    mapfunc["1D_cutflow" + samplename + channel] -> SetBinContent(3, current_eff);

    current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(22 ) / denom;//two veto el charge
    mapfunc["1D_cutflow" + samplename + channel] -> SetBinContent(4, current_eff);
        
    current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(7) / denom;//two veto mu pt 
    mapfunc["1D_cutflow" + samplename + channel] -> SetBinContent(5, current_eff);

    current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(25) / denom;//mll
    mapfunc["1D_cutflow" + samplename + channel] -> SetBinContent(6, current_eff);

    current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(27) / denom;//Njet
    mapfunc["1D_cutflow" + samplename + channel] -> SetBinContent(7, current_eff);

    current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(29) / denom;//mN
    mapfunc["1D_cutflow" + samplename + channel] -> SetBinContent(8, current_eff);

    current_eff = maphist["signal_eff" + samplename + channel] -> GetBinContent(19) / denom;//SR                                                                                                                                                                                                                                                                
    mapfunc["1D_cutflow" + samplename + channel] -> SetBinContent(9, current_eff);
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
  
  maphist_2D["All_eff_" + channel] -> Fill(mZp_double, mN_double,  evt_all / denom);
  maphist_2D["AK8_0_eff_" + channel] -> Fill(mZp_double, mN_double, evt_AK8_0 / denom);
  maphist_2D["AK8_1_eff_" + channel] -> Fill(mZp_double, mN_double, evt_AK8_1 / denom);
  maphist_2D["AK8_2_eff_" + channel] -> Fill(mZp_double, mN_double, evt_AK8_2 / denom);
 

  current_file -> Close();
  delete current_file;
}


void draw_2D_hist(TString content){
  
  mapcanvas[ content ] = new TCanvas( content ,"",800,600);
  gStyle -> SetOptStat(0);
  mapcanvas[ content ] -> cd();
  
  gStyle -> SetOptStat(0);
  gStyle->SetPaintTextFormat("4.3f");
  
  // === Deal with crazy points
  int crazy_bin = maphist_2D[ content ] -> FindBin(4900., 1200.);
  int sidebin_1 = maphist_2D[ content ] -> FindBin(4900., 1000.);
  int sidebin_2 = maphist_2D[ content ] -> FindBin(4900., 1400.);
  int sidebin_3 = maphist_2D[ content ] -> FindBin(4500., 1200.);
  double sidebin_1_content = maphist_2D[ content ] -> GetBinContent(sidebin_1);
  double sidebin_2_content = maphist_2D[ content ] -> GetBinContent(sidebin_2);
  double sidebin_3_content = maphist_2D[ content ] -> GetBinContent(sidebin_3);

  double side_avg = (sidebin_1_content + sidebin_2_content + sidebin_3_content) / 3.0 ;
  maphist_2D[ content ] -> SetBinContent(crazy_bin, side_avg);

  maphist_2D[ content ] -> Draw("COLZTEXT");
  
  
  
  mapcanvas[ content ] -> Update();
  mapcanvas[ content ] -> SaveAs("./pdfs/" + content + "_noBcut.pdf");
}


void make_1D_hist(int Zpmass, TString channel){
  
  TString current_mZp = TString::Itoa(Zpmass,10);
  const int N_label = 9;
  TString label[N_label] = {"MET filter", "HLT", "Two Lep", "Sign(ll)", "Pt(l)", "mll", "Njet", "mN", "SR"};
  
  
  cout << "making 1D plot for Zp mass of " << Zpmass << " channel " << channel << endl;
  TString map_name;
  map_name = channel + current_mZp + "1D";
  mapcanvas[map_name] = new TCanvas(map_name,"",1500,800);
  gStyle -> SetOptStat(0);
  gStyle->SetPalette(55);//get new palete kDarkBodyRadiator
  mapcanvas[map_name] -> SetTopMargin( 0.05 );
  mapcanvas[map_name] -> SetBottomMargin( 0.13 );
  mapcanvas[map_name] -> SetRightMargin( 0.20 );
  mapcanvas[map_name] -> SetLeftMargin( 0.16 );
  mapcanvas[map_name] -> cd();
  //mapcanvas[map_name] -> SetLogy();
  maplegend[map_name] = new TLegend(0.80, 0.13, 0.95, 0.95);

  TString current_file_head = "1D_cutflowZP" + current_mZp + "_N";
  mapfunc[current_file_head + "100" + channel] -> Draw("hist");
  mapfunc[current_file_head + "100" + channel] -> GetYaxis()->SetLabelSize(0.05);
  mapfunc[current_file_head + "100" + channel] -> GetYaxis()->SetTitleSize(0.07);
  mapfunc[current_file_head + "100" + channel] -> GetYaxis()->SetTitleOffset(1.02);
  mapfunc[current_file_head + "100" + channel] -> GetYaxis()->SetRangeUser(0.000, 1.1);
  mapfunc[current_file_head + "100" + channel] -> GetXaxis()->SetLabelSize(0.05);
  mapfunc[current_file_head + "100" + channel] -> GetXaxis()->SetRangeUser(0, N_label);
  
  mapfunc[current_file_head + "100" + channel] -> SetLineColor(632);
  for(int i_label = 0; i_label < N_label; i_label++){
    mapfunc[current_file_head + "100" + channel] -> GetXaxis()-> SetBinLabel(i_label + 1, label[i_label]);
  } 
  
  maplegend[map_name] -> AddEntry(mapfunc[current_file_head + "100" + channel], "ZP" + current_mZp + "N100", "l");
  /*
  Int_t line_colors[11] = {808, 805, 802, 791, 400, 
			   390, 830, 825, 813, 860, 880};
  */
  Int_t line_colors[11] = {800, 400, 416, 600, 880,
			   632, 800, 400, 416, 600, 880};

  int mN = 300;
  int i_color = 0;
  while(mN < Zpmass / 2){
    TString current_mN = TString::Itoa(mN,10);
    mapfunc[current_file_head + current_mN + channel] -> SetLineColor(line_colors[i_color]); 
    if(i_color > 4) mapfunc[current_file_head + current_mN + channel] -> SetLineStyle(5);
    mapfunc[current_file_head + current_mN + channel] -> Draw("histsame");
    maplegend[map_name] -> AddEntry(mapfunc[current_file_head + current_mN + channel], "ZP" + current_mZp + "N" + current_mN, "l");

    i_color++;
    mN = mN + 200;
  }

  
  maplegend[map_name] -> Draw("lcsame");
  
  mapcanvas[ map_name ] -> Update();
  mapcanvas[ map_name ] -> SaveAs("./cutflow/" + channel + "_cutflow_Zp_" + TString::Itoa(Zpmass, 10) + ".pdf");
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
  
  maphist_2D["All_eff_MuMu"] = new TH2D("", "", 12, 400., 5200., 12, 100., 2500.);
  maphist_2D["AK8_0_eff_MuMu"] = new TH2D("", "", 12, 400., 5200., 12, 100., 2500.);
  maphist_2D["AK8_1_eff_MuMu"] = new TH2D("", "", 12, 400., 5200., 12, 100., 2500.);
  maphist_2D["AK8_2_eff_MuMu"] = new TH2D("", "", 12, 400., 5200., 12, 100., 2500.);
  
  maphist_2D["All_eff_EE"] = new TH2D("", "", 12, 400., 5200., 12, 100., 2500.);
  maphist_2D["AK8_0_eff_EE"] = new TH2D("", "", 12, 400., 5200., 12, 100., 2500.);
  maphist_2D["AK8_1_eff_EE"] = new TH2D("", "", 12, 400., 5200., 12, 100., 2500.);
  maphist_2D["AK8_2_eff_EE"] = new TH2D("", "", 12, 400., 5200., 12, 100., 2500.);
  
  /*
  maphist_2D["Medium_eff"] = new TH2D("", "", 12, 400., 5200., 12, 100., 2500.);
  maphist_2D["Tight_eff"] = new TH2D("", "", 12, 400., 5200., 12, 100., 2500.);
  maphist_2D["HighPt_eff"] = new TH2D("", "", 12, 400., 5200., 12, 100., 2500.);
  */
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
  
  for(int i_Zp = 0; i_Zp < 12; i_Zp ++){
    int mZp = 400 + 400 * i_Zp;
    int mN = 100;
    while(mN < mZp / 2){
      TString current_mZp = TString::Itoa(mZp,10);
      TString current_mN = TString::Itoa(mN,10);

      openfile_signal("ZP" + current_mZp + "_N" + current_mN, "EE");

      mN = mN + 200;
    }
  }
  
  for(int i_Zp = 0; i_Zp < 12; i_Zp ++){
    int mZp = 400 + 400 * i_Zp;
    TString current_mZp = TString::Itoa(mZp,10);
    
    make_1D_hist(mZp, "EE");
    make_1D_hist(mZp, "MuMu");

  }

  
  draw_2D_hist("All_eff_MuMu");
  draw_2D_hist("AK8_0_eff_MuMu");
  draw_2D_hist("AK8_1_eff_MuMu");
  draw_2D_hist("AK8_2_eff_MuMu");
  draw_2D_hist("All_eff_EE");
  draw_2D_hist("AK8_0_eff_EE");
  draw_2D_hist("AK8_1_eff_EE");
  draw_2D_hist("AK8_2_eff_EE");
  
  
}// End of Main Function ////////////////////////////////////////////////////// 
