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
void openfile(TString filename){
  
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
    if (!cl->InheritsFrom("TH1")) histnames.push_back(key -> GetName());
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
    if (!cl->InheritsFrom("TH1")) histnames.push_back(key -> GetName());
    histnames.push_back(key -> GetName());
  }

  //return histnames_2D;
  return histnames;
  
  gDirectory -> cd();
}

void draw_dphi_plots(){
  TString sex = "sex";
  mapcanvas[sex] = new TCanvas("", "", 800, 600);
  gStyle -> SetOptStat(1111);
  mapcanvas[sex] -> SetTopMargin( 0.05 );
  mapcanvas[sex] -> SetBottomMargin( 0.13 );
  mapcanvas[sex] -> SetRightMargin( 0.05 );
  mapcanvas[sex] -> SetLeftMargin( 0.16 );
  mapcanvas[sex] -> cd();
  
  
  TString central = "SingleMuonTrigger_Dijet_central_all_MUON_HN_LOOSEv7_SIP3_Nvtx_reweight_events_pt_cone_vs_eta_F_1D";
  TString dphi_1p5 = "SingleMuonTrigger_Dijet_awayjet_dphi_1p5_all_MUON_HN_LOOSEv7_SIP3_Nvtx_reweight_events_pt_cone_vs_eta_F_1D";
  TString dphi_2 = "SingleMuonTrigger_Dijet_awayjet_dphi_2_all_MUON_HN_LOOSEv7_SIP3_Nvtx_reweight_events_pt_cone_vs_eta_F_1D";
  TString dphi_3 = "SingleMuonTrigger_Dijet_awayjet_dphi_3_all_MUON_HN_LOOSEv7_SIP3_Nvtx_reweight_events_pt_cone_vs_eta_F_1D";
  
  maphist[central] -> SetTitle("");
  gStyle -> SetOptStat(0);
  maphist[central] -> GetXaxis() -> SetTitle("P_{T}^{cone} (#mu)");
  maphist[central] -> GetYaxis() -> SetTitle("Fate Rate");
  maphist[central] -> GetYaxis() -> SetRange(10., 60.);
  maphist[central] -> SetMinimum(0.);
  maphist[central] -> SetMaximum(0.35);
  maphist[central] -> SetMarkerSize(10.);
  maphist[central] -> SetLineColor(kBlack);
  maphist[central] -> Draw("E1");
  
  maphist[dphi_1p5] -> SetMarkerColor(kRed);
  maphist[dphi_2] -> SetMarkerColor(kBlue);
  maphist[dphi_3] -> SetMarkerColor(kGreen);
  maphist[dphi_1p5] -> SetLineColor(kRed);
  maphist[dphi_2] -> SetLineColor(kBlue);
  maphist[dphi_3] -> SetLineColor(kGreen);

  maphist[dphi_1p5] -> Draw("E1same");
  maphist[dphi_2] -> Draw("E1same");
  maphist[dphi_3] -> Draw("E1same");
  
  maplegend[sex] = new TLegend(0.65, 0.7, 0.95, 0.95);
  maplegend[sex] -> AddEntry(maphist[central], "Central, #Delta#phi : 2.5", "l");
  maplegend[sex] -> AddEntry(maphist[dphi_1p5], "#Delta#phi : 1.5", "l");
  maplegend[sex] -> AddEntry(maphist[dphi_2], "#Delta#phi : 2", "l");
  maplegend[sex] -> AddEntry(maphist[dphi_3], "#Delta#phi : 3", "l");
  maplegend[sex] -> Draw("same");
  
  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Simulation}}}");
  latex_Lumi.SetTextSize(0.035);
    
  TString pdfname = "./pdfs/";
  pdfname.Append("Dphi_1D_Muon.pdf");
  
  mapcanvas[sex] -> SaveAs(pdfname);
}

void draw_pj_over_pl_plots(){
  TString sex = "sex";
  mapcanvas[sex] = new TCanvas("", "", 800, 600);
  gStyle -> SetOptStat(1111);
  mapcanvas[sex] -> SetTopMargin( 0.05 );
  mapcanvas[sex] -> SetBottomMargin( 0.13 );
  mapcanvas[sex] -> SetRightMargin( 0.05 );
  mapcanvas[sex] -> SetLeftMargin( 0.16 );
  mapcanvas[sex] -> cd();

  TString central = "SingleMuonTrigger_Dijet_central_all_MUON_HN_LOOSEv7_SIP3_Nvtx_reweight_events_pt_cone_vs_eta_F_1D";
  TString pj_over_pl_0p8 = "SingleMuonTrigger_Dijet_pj_over_pl_0p8_all_MUON_HN_LOOSEv7_SIP3_Nvtx_reweight_events_pt_cone_vs_eta_F_1D";
  TString pj_over_pl_1p2 = "SingleMuonTrigger_Dijet_pj_over_pl_1p2_all_MUON_HN_LOOSEv7_SIP3_Nvtx_reweight_events_pt_cone_vs_eta_F_1D";
  TString pj_over_pl_1p5 = "SingleMuonTrigger_Dijet_pj_over_pl_1p5_all_MUON_HN_LOOSEv7_SIP3_Nvtx_reweight_events_pt_cone_vs_eta_F_1D";
  TString pj_over_pl_2 = "SingleMuonTrigger_Dijet_pj_over_pl_2_all_MUON_HN_LOOSEv7_SIP3_Nvtx_reweight_events_pt_cone_vs_eta_F_1D";
  
  
  maphist[central] -> SetTitle("");
  gStyle -> SetOptStat(0);
  maphist[central] -> GetXaxis() -> SetTitle("P_{T}^{cone} (#mu)");
  maphist[central] -> GetYaxis() -> SetTitle("Fate Rate");
  maphist[central] -> GetYaxis() -> SetRange(10., 60.);
  maphist[central] -> SetMinimum(0.);
  maphist[central] -> SetMaximum(0.35);
  maphist[central] -> SetMarkerSize(10.);
  maphist[central] -> SetLineColor(kBlack);
  maphist[central] -> Draw("E1");

  maphist[pj_over_pl_0p8] -> SetLineColor(kRed);
  maphist[pj_over_pl_1p2] -> SetLineColor(kBlue);
  maphist[pj_over_pl_1p5] -> SetLineColor(kYellow);
  maphist[pj_over_pl_2] -> SetLineColor(kGreen);
 
  
  maphist[pj_over_pl_0p8] -> Draw("E1same");
  maphist[pj_over_pl_1p2] -> Draw("E1same");
  maphist[pj_over_pl_1p5] -> Draw("E1same");
  maphist[pj_over_pl_2] -> Draw("E1same");
  
  
  maplegend[sex] = new TLegend(0.65, 0.76, 0.95, 0.95);
  maplegend[sex] -> AddEntry(maphist[central], "Central, P_{T}^{j} / P_{T}^{#mu} : 1", "l");
  maplegend[sex] -> AddEntry(maphist[pj_over_pl_0p8], "P_{T}^{j} / P_{T}^{#mu} : 0.8", "l");
  maplegend[sex] -> AddEntry(maphist[pj_over_pl_1p2], "P_{T}^{j} / P_{T}^{#mu} : 1.2", "l");
  maplegend[sex] -> AddEntry(maphist[pj_over_pl_1p5], "P_{T}^{j} / P_{T}^{#mu} : 1.5", "l");
  maplegend[sex] -> AddEntry(maphist[pj_over_pl_2], "P_{T}^{j} / P_{T}^{#mu} : 2", "l");

  
  maplegend[sex] -> Draw("same");

  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Simulation}}}");
  latex_Lumi.SetTextSize(0.035);

  TString pdfname = "./pdfs/";
  pdfname.Append("Pj_over_Pl_1D_Muon.pdf");

  mapcanvas[sex] -> SaveAs(pdfname);
}
/// Main Function ////////////////////////////////////////////////////////////
void plot(){
    
  //openfile("Data_driven_FR_syst_DoubleMuon_Nvtx_Reweight.root");
  openfile("Data_driven_FR_syst_DoubleMuon_Nvtx_Reweight.root");
  cout << "open files complete" << endl;
  
  draw_dphi_plots();
  draw_pj_over_pl_plots();



}// End of Main Function ////////////////////////////////////////////////////// 
