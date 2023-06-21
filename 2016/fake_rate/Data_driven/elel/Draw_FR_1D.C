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
#include <stdlib.h>
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
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>

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
map<TString, TH2*> maphist_2D;
map<TString, TGraphErrors*> mapTGraphErrors;
map<TString, TGraphAsymmErrors*> mapTGraphAsymmErrors;


/// Getting Histogram Function ///////////////////////////////////////////////
TH1 * GetHist(TString hname){
  
  TH1 *h = NULL;
  std::map<TString, TH1*>::iterator mapit = maphist.find(hname);
  if(mapit != maphist.end()) return mapit-> second;
  
  return h;
  
}
//////////////////////////////////////////////////////////////////////////////


TH2 * GetHist_2D(TString hname){
  
  TH2 *h = NULL;
  std::map<TString, TH2*>::iterator mapit = maphist_2D.find(hname);
  if(mapit != maphist_2D.end()) return mapit-> second;
  
  return h;
  
}


/// Open ROOT file ///////////////////////////////////////////////////////////
void openfile(TString filename){
  
  cout << "opening : " << filename << endl;
  
  mapfile[filename] = new TFile ((filename)) ;
  
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
    cout << histnames.at(i) << endl;
    maphist[histnames.at(i)] = (TH1*)gDirectory -> Get(histnames.at(i));
    //cout << histnames.at(i) + cyclename + samplename << endl;
  }
  for(unsigned int i = 0; i < histnames_2D.size(); i ++){
    cout << histnames_2D.at(i) << endl;
    maphist_2D[histnames_2D.at(i)] = (TH2*)gDirectory -> Get(histnames_2D.at(i));
  }
    
  gDirectory -> cd();
}
//////////////////////////////////////////////////////////////////////////////


/// Make Output Histogram ////////////////////////////////////////////////////
void makehistogram(TString map_string, TString MC_scale){
  
  mapcanvas[map_string] = new TCanvas(map_string,"",800,800);
  gStyle -> SetOptStat(1111);
  mapcanvas[map_string] -> SetTopMargin( 0.05 );
  mapcanvas[map_string] -> SetBottomMargin( 0.13 );
  mapcanvas[map_string] -> SetRightMargin( 0.05 );
  mapcanvas[map_string] -> SetLeftMargin( 0.16 );
  mapcanvas[map_string] -> cd();
  
  mapTGraphErrors["FR_EC" + MC_scale] -> SetTitle("");
  gStyle -> SetOptStat(0);
  mapTGraphErrors["FR_EC" + MC_scale] -> GetXaxis() -> SetTitle("P_{T}^{cone} (#mu)");
  mapTGraphErrors["FR_EC" + MC_scale] -> GetYaxis() -> SetTitle("Fate Rate");
  mapTGraphErrors["FR_EC" + MC_scale] -> SetMinimum(-0.05);
  mapTGraphErrors["FR_EC" + MC_scale] -> SetMaximum(0.2);
  mapTGraphErrors["FR_EC" + MC_scale] -> SetMarkerColor(kGreen);
  mapTGraphErrors["FR_EC" + MC_scale] -> SetLineColor(kGreen);
  mapTGraphErrors["FR_EC" + MC_scale] -> SetFillColorAlpha(kGreen,0.35);
  mapTGraphErrors["FR_EC" + MC_scale] -> Draw("E2AP");
    
  mapTGraphErrors["FR_OB" + MC_scale] -> SetMarkerColor(kRed);
  mapTGraphErrors["FR_OB" + MC_scale] -> SetLineColor(kRed);
  mapTGraphErrors["FR_OB" + MC_scale] -> SetFillColorAlpha(kRed,0.35);
  mapTGraphErrors["FR_OB" + MC_scale] -> Draw("E2same");
  
  mapTGraphErrors["FR_IB" + MC_scale] -> SetMarkerColor(kBlue);
  mapTGraphErrors["FR_IB" + MC_scale] -> SetLineColor(kBlue);
  mapTGraphErrors["FR_IB" + MC_scale] -> Draw("E1same");
  
  maplegend[map_string] = new TLegend(0.5, 0.7, 0.8, 0.95);
  maplegend[map_string] -> AddEntry(mapTGraphErrors["FR_IB" + MC_scale], "|#eta| < 0.8", "lp");
  maplegend[map_string] -> AddEntry(mapTGraphErrors["FR_OB" + MC_scale], "0.8 < |#eta| < 1.479", "f");
  maplegend[map_string] -> AddEntry(mapTGraphErrors["FR_EC" + MC_scale], "1.479 < |#eta| < 2.5", "f");
  maplegend[map_string] -> Draw("same");

  mapcanvas[map_string] -> cd();
  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_Lumi.SetTextSize(0.035);
  latex_Lumi.DrawLatex(0.7, 0.96, "35.9 fb^{-1} (13 TeV)");

    
  TString pdfname = "./Fake_Rate/";
  pdfname.Append("FR_mu" + MC_scale);
  //pdfname.Append(".pdf");
  pdfname.Append("_Trigger_SF.pdf");
  
  mapcanvas[map_string] -> SaveAs(pdfname);
  
}

void same_eta_bin(TString map_string, TString Eta_bin){
  
  mapcanvas[map_string] = new TCanvas(map_string,"",800,800);
  gStyle -> SetOptStat(1111);
  mapcanvas[map_string] -> SetTopMargin( 0.05 );
  mapcanvas[map_string] -> SetBottomMargin( 0.13 );
  mapcanvas[map_string] -> SetRightMargin( 0.05 );
  mapcanvas[map_string] -> SetLeftMargin( 0.16 );
  mapcanvas[map_string] -> cd();
  
  cout << "FR_" + Eta_bin + "_1p15" << endl;

  mapTGraphErrors["FR_" + Eta_bin + "_1p15"] -> SetTitle("");
  gStyle -> SetOptStat(0);
  mapTGraphErrors["FR_" + Eta_bin + "_1p15"] -> GetXaxis() -> SetTitle("P_{T}^{cone} (#mu)");
  mapTGraphErrors["FR_" + Eta_bin + "_1p15"] -> GetYaxis() -> SetTitle("Fate Rate");
  mapTGraphErrors["FR_" + Eta_bin + "_1p15"] -> SetMinimum(-0.05);
  mapTGraphErrors["FR_" + Eta_bin + "_1p15"] -> SetMaximum(0.2);
  mapTGraphErrors["FR_" + Eta_bin + "_1p15"] -> SetMarkerColor(kGreen);
  mapTGraphErrors["FR_" + Eta_bin + "_1p15"] -> SetLineColor(kGreen);
  mapTGraphErrors["FR_" + Eta_bin + "_1p15"] -> SetFillColorAlpha(kGreen,0.35);
  mapTGraphErrors["FR_" + Eta_bin + "_1p15"] -> Draw("E2AP");

  mapTGraphErrors["FR_" + Eta_bin + "_p85"] -> SetMarkerColor(kRed);
  mapTGraphErrors["FR_" + Eta_bin + "_p85"] -> SetLineColor(kRed);
  mapTGraphErrors["FR_" + Eta_bin + "_p85"] -> SetFillColorAlpha(kRed,0.35);
  mapTGraphErrors["FR_" + Eta_bin + "_p85"] -> Draw("E2same");

  mapTGraphErrors["FR_" + Eta_bin] -> SetMarkerColor(kBlue);
  mapTGraphErrors["FR_" + Eta_bin] -> SetLineColor(kBlue);
  mapTGraphErrors["FR_" + Eta_bin] -> Draw("E1same");
  
  maplegend[map_string] = new TLegend(0.5, 0.7, 0.8, 0.95);
  maplegend[map_string] -> AddEntry(mapTGraphErrors["FR_" + Eta_bin], "Central", "lp");
  maplegend[map_string] -> AddEntry(mapTGraphErrors["FR_" + Eta_bin + "_1p15"], "15 % Scale up prompt", "f");
  maplegend[map_string] -> AddEntry(mapTGraphErrors["FR_" + Eta_bin + "_p85"], "15 % Scale down prompt", "f");
  maplegend[map_string] -> Draw("same");

  mapcanvas[map_string] -> cd();
  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_Lumi.SetTextSize(0.035);
  latex_Lumi.DrawLatex(0.7, 0.96, "35.9 fb^{-1} (13 TeV)");


  TString pdfname = "./Fake_Rate/";
  pdfname.Append("FR_mu_" + Eta_bin);
  //pdfname.Append(".pdf");
  pdfname.Append("_Trigger_SF.pdf");
  mapcanvas[map_string] -> SaveAs(pdfname);

}


//////////////////////////////////////////////////////////////////////////////

void draw_FR_plot(TString nameofhistogram){

  mapcanvas[nameofhistogram] = new TCanvas(nameofhistogram,"",800,600);
  gStyle -> SetOptStat(1111);
  mapcanvas[nameofhistogram] -> SetTopMargin( 0.05 );
  mapcanvas[nameofhistogram] -> SetBottomMargin( 0.13 );
  mapcanvas[nameofhistogram] -> SetRightMargin( 0.05 );
  mapcanvas[nameofhistogram] -> SetLeftMargin( 0.16 );
  mapcanvas[nameofhistogram] -> cd();

  maphist_2D[nameofhistogram] -> SetTitle("");
  gStyle -> SetOptStat(0);
  maphist_2D[nameofhistogram] -> GetXaxis() -> SetTitle("P_{T}^{cone} (#mu)");
  maphist_2D[nameofhistogram] -> GetYaxis() -> SetTitle("|#eta_{#mu}|");
  maphist_2D[nameofhistogram] -> Draw("colztext1e");
  gStyle->SetPaintTextFormat("0.4f");

  mapcanvas[nameofhistogram] -> cd();
  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_Lumi.SetTextSize(0.035);
  latex_Lumi.DrawLatex(0.7, 0.96, "35.9 fb^{-1} (13 TeV)");

  TString pdfname = "./Fake_Rate/";
  //pdfname.Append("FR_mu_Nvtx");
  pdfname.Append("FR_mu_Nvtx_trigger_SF");
  pdfname.Append(".pdf");
  mapcanvas[nameofhistogram] -> SaveAs(pdfname);


}

void draw_prompt_syst(TString map_string, TString Eta_bin){
  
  mapcanvas[map_string] = new TCanvas(map_string,"",800,800);
  gStyle -> SetOptStat(1111);
  mapcanvas[map_string] -> SetTopMargin( 0.05 );
  mapcanvas[map_string] -> SetBottomMargin( 0.13 );
  mapcanvas[map_string] -> SetRightMargin( 0.05 );
  mapcanvas[map_string] -> SetLeftMargin( 0.16 );
  mapcanvas[map_string] -> cd();

  cout << "FR_" + Eta_bin + "_syst" << endl;

  mapTGraphAsymmErrors["FR_" + Eta_bin + "_syst"] -> SetTitle("");
  gStyle -> SetOptStat(0);
  mapTGraphAsymmErrors["FR_" + Eta_bin + "_syst"] -> GetXaxis() -> SetTitle("P_{T}^{cone} (#mu)");
  mapTGraphAsymmErrors["FR_" + Eta_bin + "_syst"] -> GetYaxis() -> SetTitle("Fate Rate");
  mapTGraphAsymmErrors["FR_" + Eta_bin + "_syst"] -> GetYaxis() -> SetRange(35., 110.);
  mapTGraphAsymmErrors["FR_" + Eta_bin + "_syst"] -> SetMinimum(0.);
  mapTGraphAsymmErrors["FR_" + Eta_bin + "_syst"] -> SetMaximum(0.2);
  mapTGraphAsymmErrors["FR_" + Eta_bin + "_syst"] -> SetMarkerColor(kRed);
  mapTGraphAsymmErrors["FR_" + Eta_bin + "_syst"] -> SetMarkerSize(10.);
  mapTGraphAsymmErrors["FR_" + Eta_bin + "_syst"] -> SetLineColor(kRed);
  mapTGraphAsymmErrors["FR_" + Eta_bin + "_syst"] -> SetFillColorAlpha(kRed,0.35);
  mapTGraphAsymmErrors["FR_" + Eta_bin + "_syst"] -> Draw("E1AP");
  mapTGraphAsymmErrors["FR_" + Eta_bin + "_syst"] -> Draw("E2same");

  mapTGraphErrors["FR_" + Eta_bin] -> SetMarkerColor(kBlue);
  mapTGraphErrors["FR_" + Eta_bin] -> SetLineColor(kBlue);
  //mapTGraphErrors["FR_" + Eta_bin] -> Draw("E1same");
  
  maplegend[map_string] = new TLegend(0.5, 0.8, 0.8, 0.95);
  //maplegend[map_string] -> AddEntry(mapTGraphErrors["FR_" + Eta_bin], "Central + stat.", "lp");
  maplegend[map_string] -> AddEntry(mapTGraphAsymmErrors["FR_" + Eta_bin + "_syst"], "13 % Prompt Syst.", "f");
  maplegend[map_string] -> Draw("same");

  mapcanvas[map_string] -> cd();
  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_Lumi.SetTextSize(0.035);
  latex_Lumi.DrawLatex(0.7, 0.96, "35.9 fb^{-1} (13 TeV)");


  TString pdfname = "./Fake_Rate/";
  pdfname.Append("FR_mu_" + Eta_bin);
  //pdfname.Append(".pdf");                                                                                                                                                        
  pdfname.Append("_syst.pdf");
  mapcanvas[map_string] -> SaveAs(pdfname);

}

void make_root_file(TString nameofhistogram, double FR_IB[], double FR_OB[], double FR_EC[], double FR_IB_syst_bellow[], double FR_OB_syst_bellow[], double FR_EC_syst_bellow[]){
  
  double pt_bins[6] = {15., 25., 35., 45., 55., 65.};
  TString data_F0 = nameofhistogram + "_F0_data";
  TString data_F = nameofhistogram + "_F_data";

  
  for(int i = 0; i < 6; i++){
    int bin_IB = maphist_2D[data_F]->FindBin(pt_bins[i], 0.5);
    int bin_OB = maphist_2D[data_F]->FindBin(pt_bins[i], 1.0);
    int bin_EC = maphist_2D[data_F]->FindBin(pt_bins[i], 2.0);
    
    maphist_2D[data_F0 + "FR_syst"] -> SetBinContent(bin_IB, FR_IB[i]);
    maphist_2D[data_F0 + "FR_syst"] -> SetBinError(bin_IB, FR_IB_syst_bellow[i]);
    maphist_2D[data_F0 + "FR_syst"] -> SetBinContent(bin_OB, FR_OB[i]);
    maphist_2D[data_F0 + "FR_syst"] -> SetBinError(bin_OB, FR_OB_syst_bellow[i]);
    maphist_2D[data_F0 + "FR_syst"] -> SetBinContent(bin_EC, FR_EC[i]);
    maphist_2D[data_F0 + "FR_syst"] -> SetBinError(bin_EC, FR_EC_syst_bellow[i]);
  }
  maphist_2D[data_F0 + "FR_syst"] -> GetXaxis() -> SetRangeUser(0., 60.);
  
  maphist_2D[data_F0 + "FR_syst"] -> SetName("FR_ptcone_central");
  
  TFile *MyFile = new TFile("FR_muon_central.root","RECREATE");
  maphist_2D[data_F0 + "FR_syst"] -> Write();
  MyFile->Close();
    
  
}

void draw_all_eta(TString map_string){
  
  mapcanvas[map_string] = new TCanvas(map_string,"",800,800);
  gStyle -> SetOptStat(1111);
  mapcanvas[map_string] -> SetTopMargin( 0.05 );
  mapcanvas[map_string] -> SetBottomMargin( 0.13 );
  mapcanvas[map_string] -> SetRightMargin( 0.05 );
  mapcanvas[map_string] -> SetLeftMargin( 0.16 );
  mapcanvas[map_string] -> cd();
  
  mapTGraphAsymmErrors["FR_IB_syst"] -> SetTitle("");
  gStyle -> SetOptStat(0);
  mapTGraphAsymmErrors["FR_IB_syst"] -> GetXaxis() -> SetTitle("P_{T}^{cone} (#mu)");
  mapTGraphAsymmErrors["FR_IB_syst"] -> GetYaxis() -> SetTitle("Fate Rate");
  mapTGraphAsymmErrors["FR_IB_syst"] -> GetXaxis() -> SetRangeUser(40., 85.);
  mapTGraphAsymmErrors["FR_IB_syst"] -> SetMinimum(0.);
  mapTGraphAsymmErrors["FR_IB_syst"] -> SetMaximum(0.3);
  mapTGraphAsymmErrors["FR_IB_syst"] -> SetMarkerColor(kRed);
  mapTGraphAsymmErrors["FR_IB_syst"] -> SetMarkerSize(10.);
  mapTGraphAsymmErrors["FR_IB_syst"] -> SetLineColor(kRed);
  mapTGraphAsymmErrors["FR_IB_syst"] -> SetFillColorAlpha(kRed,0.35);
  mapTGraphAsymmErrors["FR_IB_syst"] -> Draw("E1AP");
  
  mapTGraphAsymmErrors["FR_OB_syst"] -> SetLineColor(kGreen);
  mapTGraphAsymmErrors["FR_OB_syst"] -> Draw("E1same");
  
  mapTGraphAsymmErrors["FR_EC_syst"] -> SetLineColor(kBlue);
  mapTGraphAsymmErrors["FR_EC_syst"] -> Draw("E1same");

  
  maplegend[map_string] = new TLegend(0.5, 0.8, 0.8, 0.95);
  maplegend[map_string] -> AddEntry(mapTGraphAsymmErrors["FR_IB_syst"], "Fake Rate, Inner Barrel", "lp");
  maplegend[map_string] -> AddEntry(mapTGraphAsymmErrors["FR_OB_syst"], "Fake Rate, Outer Barrel", "lp");
  maplegend[map_string] -> AddEntry(mapTGraphAsymmErrors["FR_EC_syst"], "Fake Rate, Endcap", "lp");
  maplegend[map_string] -> Draw("same");
  
  mapcanvas[map_string] -> cd();
  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_Lumi.SetTextSize(0.035);
  latex_Lumi.DrawLatex(0.7, 0.96, "35.9 fb^{-1} (13 TeV)");
  
  TString pdfname = "./Fake_Rate/";
  pdfname.Append("FR_mu_");
  pdfname.Append("all_eta.pdf");
  mapcanvas[map_string] -> SaveAs(pdfname);
  
}


void draw_histogram(TString nameofhistogram){

  TString data_F = nameofhistogram + "_F_data";
  TString data_F0 = nameofhistogram + "_F0_data";
  TString MC_F = nameofhistogram + "_F_MC";
  TString MC_F0 = nameofhistogram + "_F0_MC";
  
  maphist_2D[data_F + "1p15"] = (TH2*)GetHist_2D(data_F) -> Clone();
  maphist_2D[data_F + "p85"] =(TH2*)GetHist_2D(data_F) -> Clone();
  maphist_2D[data_F0 + "1p15"] =(TH2*)GetHist_2D(data_F0) -> Clone();
  maphist_2D[data_F0 + "p85"] =(TH2*)GetHist_2D(data_F0) -> Clone();

  maphist_2D[data_F0 + "FR_syst"] =(TH2*)GetHist_2D(data_F0) -> Clone();
  
  maphist_2D[MC_F + "1p15"] = (TH2*)GetHist_2D(MC_F) -> Clone();// increase MC 15%
  maphist_2D[MC_F + "p85"] = (TH2*)GetHist_2D(MC_F) -> Clone();// decrease MC 15%
  maphist_2D[MC_F + "1p15"] -> Scale(1.13);
  maphist_2D[MC_F + "p85"] -> Scale(0.87);
  maphist_2D[MC_F0 + "1p15"] = (TH2*)GetHist_2D(MC_F0) -> Clone();// increase MC 15%
  maphist_2D[MC_F0 + "p85"] = (TH2*)GetHist_2D(MC_F0) -> Clone();// decrease MC 15%
  maphist_2D[MC_F0 + "1p15"] -> Scale(1.13);
  maphist_2D[MC_F0 + "p85"] -> Scale(0.87);
  
  cout << GetHist_2D(MC_F) -> Integral() << ", " << maphist_2D[MC_F + "p85"] -> Integral() << ", " << maphist_2D[MC_F + "1p15"] -> Integral() << endl;
  
  TString syst_string[3] = {"", "p85", "1p15"};
  for(int i = 0; i < 3; i++){
    cout << MC_F + syst_string[i] << "-------" << endl;
    
    maphist_2D[data_F + syst_string[i]] ->  Add(GetHist_2D(MC_F + syst_string[i]), -1);
    maphist_2D[data_F0 + syst_string[i]] ->  Add(GetHist_2D(MC_F0 + syst_string[i]), -1);
    GetHist_2D(data_F + syst_string[i]) -> Divide(GetHist_2D(data_F0 + syst_string[i]));
  }
  
  double FR_IB[6], FR_OB[6], FR_EC[6];
  double FR_IB_1p15[6], FR_OB_1p15[6], FR_EC_1p15[6]; //increase prompt 15%
  double FR_IB_p85[6], FR_OB_p85[6], FR_EC_p85[6]; // decrease prompt 15%

  double FR_IB_e[6], FR_OB_e[6], FR_EC_e[6];
  double FR_IB_1p15_e[6], FR_OB_1p15_e[6], FR_EC_1p15_e[6];
  double FR_IB_p85_e[6], FR_OB_p85_e[6], FR_EC_p85_e[6];
  
  
  double pt_bins[6] = {47.5, 60., 75., 102.5, 135., 172.5};
  //double pt_bins_e[6] = {2.5, 2.5, 2.5, 2.5, 2.5, 2.5}; 
  double pt_bins_e[6] = {7.5, 5., 10., 17.5, 15., 22.5};
  

  cout << "for for TGraphError" << endl;
  for(int i = 0; i < 6; i++){
    int bin_IB = maphist_2D[data_F]->FindBin(pt_bins[i], 0.5);
    int bin_OB = maphist_2D[data_F]->FindBin(pt_bins[i], 1.0);
    int bin_EC = maphist_2D[data_F]->FindBin(pt_bins[i], 2.0);
    
    FR_IB[i] =  double(maphist_2D[data_F]->GetBinContent(bin_IB));
    FR_OB[i] =  double(maphist_2D[data_F]->GetBinContent(bin_OB));
    FR_EC[i] =  double(maphist_2D[data_F]->GetBinContent(bin_EC));
    FR_IB_e[i] =  double(maphist_2D[data_F]->GetBinError(bin_IB));
    FR_OB_e[i] =  double(maphist_2D[data_F]->GetBinError(bin_OB));
    FR_EC_e[i] =  double(maphist_2D[data_F]->GetBinError(bin_EC));
    
    FR_IB_p85[i] =  double(maphist_2D[data_F + "p85"]->GetBinContent(bin_IB));
    FR_OB_p85[i] =  double(maphist_2D[data_F + "p85"]->GetBinContent(bin_OB));
    FR_EC_p85[i] =  double(maphist_2D[data_F + "p85"]->GetBinContent(bin_EC));
    FR_IB_p85_e[i] =  double(maphist_2D[data_F + "p85"]->GetBinError(bin_IB));
    FR_OB_p85_e[i] =  double(maphist_2D[data_F + "p85"]->GetBinError(bin_OB));
    FR_EC_p85_e[i] =  double(maphist_2D[data_F + "p85"]->GetBinError(bin_EC));

    FR_IB_1p15[i] =  double(maphist_2D[data_F + "1p15"]->GetBinContent(bin_IB));
    FR_OB_1p15[i] =  double(maphist_2D[data_F + "1p15"]->GetBinContent(bin_OB));
    FR_EC_1p15[i] =  double(maphist_2D[data_F + "1p15"]->GetBinContent(bin_EC));
    FR_IB_1p15_e[i] =  double(maphist_2D[data_F + "1p15"]->GetBinError(bin_IB));
    FR_OB_1p15_e[i] =  double(maphist_2D[data_F + "1p15"]->GetBinError(bin_OB));
    FR_EC_1p15_e[i] =  double(maphist_2D[data_F + "1p15"]->GetBinError(bin_EC));
    
    cout << "FR_IB[i] : " << FR_IB[i] << ", FR_IB_e[i] : " << FR_IB_e[i] << endl;
  }
  
  mapTGraphErrors["FR_IB"] = new TGraphErrors(6, pt_bins, FR_IB, pt_bins_e, FR_IB_e);
  mapTGraphErrors["FR_OB"] = new TGraphErrors(6, pt_bins, FR_OB, pt_bins_e, FR_OB_e);
  mapTGraphErrors["FR_EC"] = new TGraphErrors(6, pt_bins, FR_EC, pt_bins_e, FR_EC_e);

  mapTGraphErrors["FR_IB"]  -> Draw();
  
  mapTGraphErrors["FR_IB_p85"] = new TGraphErrors(6, pt_bins, FR_IB_p85, pt_bins_e, FR_IB_p85_e);
  mapTGraphErrors["FR_OB_p85"] = new TGraphErrors(6, pt_bins, FR_OB_p85, pt_bins_e, FR_OB_p85_e);
  mapTGraphErrors["FR_EC_p85"] = new TGraphErrors(6, pt_bins, FR_EC_p85, pt_bins_e, FR_EC_p85_e);

  mapTGraphErrors["FR_IB_1p15"] = new TGraphErrors(6, pt_bins, FR_IB_1p15, pt_bins_e, FR_IB_1p15_e);
  mapTGraphErrors["FR_OB_1p15"] = new TGraphErrors(6, pt_bins, FR_OB_1p15, pt_bins_e, FR_OB_1p15_e);
  mapTGraphErrors["FR_EC_1p15"] = new TGraphErrors(6, pt_bins, FR_EC_1p15, pt_bins_e, FR_EC_1p15_e);

  double FR_IB_syst_upper[6];
  double FR_IB_syst_bellow[6];
  double FR_OB_syst_upper[6];
  double FR_OB_syst_bellow[6];
  double FR_EC_syst_upper[6];
  double FR_EC_syst_bellow[6];

  for(int i = 0; i < 6; i ++){
    FR_IB_syst_upper[i] = FR_IB_p85[i] - FR_IB[i];
    FR_IB_syst_bellow[i] = FR_IB[i] - FR_IB_1p15[i];

    FR_OB_syst_upper[i] = FR_OB_p85[i] - FR_OB[i];
    FR_OB_syst_bellow[i] = FR_OB[i] - FR_OB_1p15[i];

    FR_EC_syst_upper[i] = FR_EC_p85[i] - FR_EC[i];
    FR_EC_syst_bellow[i] = FR_EC[i] - FR_EC_1p15[i];
  }
  
  //FR_IB[3] = (FR_IB[2] / FR_IB_syst_bellow[2] + FR_IB[4] / FR_IB_syst_bellow[4]) / ( 1./FR_IB_syst_bellow[2] + 1./FR_IB_syst_bellow[4]);
  //FR_IB_syst_bellow[3] = sqrt(2) * (FR_IB_syst_bellow[2] * FR_IB_syst_bellow[4]) / ( FR_IB_syst_bellow[2] + FR_IB_syst_bellow[4] );
  //FR_IB_syst_upper[3] = FR_IB_syst_bellow[3];
  
  mapTGraphAsymmErrors["FR_IB_syst"] = new TGraphAsymmErrors(6, pt_bins, FR_IB, pt_bins_e, pt_bins_e, FR_IB_syst_bellow, FR_IB_syst_upper);
  mapTGraphAsymmErrors["FR_OB_syst"] = new TGraphAsymmErrors(6, pt_bins, FR_OB, pt_bins_e, pt_bins_e, FR_OB_syst_bellow, FR_OB_syst_upper);
  mapTGraphAsymmErrors["FR_EC_syst"] = new TGraphAsymmErrors(6, pt_bins, FR_EC, pt_bins_e, pt_bins_e, FR_EC_syst_bellow, FR_EC_syst_upper);
  
  /*
  makehistogram("sex", "");
  makehistogram("fuck", "_p85");
  makehistogram("suck", "_1p15");
  */
  /*
    same_eta_bin("sex", "IB");
    same_eta_bin("fuck", "OB");
    same_eta_bin("suck", "EC");
  */
  //draw_FR_plot(data_F);
  /*
  draw_prompt_syst("sex", "IB");
  draw_prompt_syst("unki", "OB");
  draw_prompt_syst("fuck", "EC");
  */
  draw_all_eta("sibal");
  //make_root_file(nameofhistogram, FR_IB, FR_OB, FR_EC, FR_IB_syst_bellow, FR_OB_syst_bellow, FR_EC_syst_bellow);
  
}


/// Main Function ////////////////////////////////////////////////////////////
void plot(){
    
  openfile("FR_electron_plots.root");
  cout << "file open complete" << endl;
  
  //draw_histogram("SingleMuonTrigger_Dijet_Awayjet_40_all_MUON_HN_LOOSEv7_SIP3_events_pt_cone_vs_eta");
  //draw_histogram("SingleMuonTrigger_Dijet_Awayjet_40_all_MUON_HN_LOOSEv7_SIP3_Nvtx_reweight_events_pt_cone_vs_eta");
  //draw_histogram("SingleMuonTrigger_Dijet_central_all_MUON_HN_LOOSEv7_SIP3_Nvtx_reweight_events_pt_cone_vs_eta");
  draw_histogram("SingleElectronTrigger_Dijet_central_all_ELECTRON_SUSY_HNPAIR_LOOSE_central_Nvtx_reweight_events_pt_cone_vs_eta");
  //SingleMuonTrigger_Dijet_Awayjet_40_all_MUON_HN_LOOSEv7_SIP3_Nvtx_reweight_events_pt_cone_vs_eta_F
		    /*
  draw_histogram("SingleMuonTrigger_Dijet_Awayjet_40_all_MUON_HN_LOOSESTv5_SIP_10_pt_cone_F", 0., 500., 1., 1500.); 
  draw_histogram("SingleMuonTrigger_Dijet_Awayjet_40_all_MUON_HN_LOOSESTv5_SIP_10_pt_cone_F0", 0., 500., 1., 8000.);
  draw_histogram("SingleMuonTrigger_Dijet_Awayjet_40_all_MUON_HN_LOOSESTv5_SIP_10_eta_F", -3., 3., 1., 500.);
  draw_histogram("SingleMuonTrigger_Dijet_Awayjet_40_all_MUON_HN_LOOSESTv5_SIP_10_eta_F0", -3., 3., 1., 5000.);
  */
  
  
  //draw_histogram("h_tri_Z_tag_dimu_mass", 0., 500., 1., 1000.);
  //draw_histogram("h_tri_non_Z_tag_pt", 0., 500., 5., 1500.);
  

}// End of Main Function ////////////////////////////////////////////////////// 
