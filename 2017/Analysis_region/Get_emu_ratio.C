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
#include <TF1.h>
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
#include <TLatex.h>


using namespace std;

map<TString, TH1*> maphist;
map<TString, TH2*> maphist_2D;
map<TString, TF1*> map_fitfunc;
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
  
  TString samplename;
  if(filename.Contains("TTLL")) samplename = "TTLL";
  else samplename = "TT";
  
  TIter next_main(gDirectory->GetListOfKeys());
  TKey *key_main;
  vector<TString> histnames_main;
  while ((key_main = (TKey*)next_main())) {
    TClass *cl = gROOT->GetClass(key_main->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    histnames_main.push_back(key_main -> GetName());
  }

  for(int k = 0; k < histnames_main.size(); k ++){
    maphist[histnames_main.at(k) + samplename] = (TH1*)gDirectory -> Get(histnames_main.at(k));
  }


  //cout << samplename << "'s histograms are" << endl;
  TString regions[3] = {"CR_Zmass", "CR_ttbar", "SR"};
  TString channels[3] = {"DiEle", "DiMu", "EMu"};
  //TString charges[2] = {"OS", "SS"};
  int i_dir = 0;
  TString directories[9];
  for(int j = 0; j < 3; j++){
    for(int k = 0; k < 3; k++){
      directories[i_dir] = regions[j] + "_" + channels[k];
      i_dir++;
    }
  }
  
  int N_directories = 9;
  
  for(int i = 0; i < i_dir; i++){
    gDirectory->cd(directories[i]);
    TIter next(gDirectory->GetListOfKeys());
    TKey *key;
    vector<TString> histnames;
    while ((key = (TKey*)next())) {
      TClass *cl = gROOT->GetClass(key->GetClassName());
      if (!cl->InheritsFrom("TH1")) continue;
      histnames.push_back(key -> GetName());
    }
    
    for(int k = 0; k < histnames.size(); k ++){
      maphist[histnames.at(k) + samplename] = (TH1*)gDirectory -> Get(histnames.at(k));
    }
    gDirectory->cd("../");
  }//for region
    
}
//////////////////////////////////////////////////////////////////////////////


void make_hist_with_overflow(TString nameofhistogram, TString sample, int N_bin, double binx[]){
  
  Int_t nx    = GetHist(nameofhistogram + sample) -> GetNbinsX()+1;
  Double_t x1 = GetHist(nameofhistogram + sample) -> GetBinLowEdge(1);
  Double_t bw = GetHist(nameofhistogram + sample) -> GetBinWidth(nx);
  Double_t x2 = GetHist(nameofhistogram + sample) -> GetBinLowEdge(nx)+bw;
  //cout << "nx : " << nx << ", x1 : " << x1 << ", bw : " << bw << ", x2 : " << x2 <<  endl;                                                                                                                                                                                

  TH1F *htmp = new TH1F("", "", nx, x1, x2);

  for (Int_t j = 1; j <= nx; j++) {
    htmp -> SetBinContent(j, GetHist(nameofhistogram + sample) -> GetBinContent(j));
    htmp -> SetBinError(j, GetHist(nameofhistogram + sample) -> GetBinError(j));
  }
  maphist[nameofhistogram + sample + "rebin"] = htmp -> Rebin(N_bin, nameofhistogram + sample + "rebin", binx);
}

void make_ratio_plot(TString nameofhistogram, TString sample, TString region, TString channel){
  
  mapfunc[nameofhistogram + channel + sample + region] = (TH1F*)GetHist(nameofhistogram + "_" + channel + sample + "rebin") -> Clone(nameofhistogram + channel + sample + region);
  mapfunc[nameofhistogram + channel + sample + region] -> Divide(GetHist(nameofhistogram + "_" + region + "_OS_EMu" + sample + "rebin"));

  //mapfunc[nameofhistogram + channel + sample + region] -> Draw();
}



void save_ratio_plot(TString nameofhistogram, TString sample, TString region, TString channel){
  
  TString title_y = "/m(e#mu jjjj)";
  TString title_x = "m(Z') (GeV)";
  TString directory;
  
  if(channel.Contains("DiEle")){
    title_y = "m(eejjjj)" + title_y;
    directory = "ee";
  }
  else if(channel.Contains("DiMu")){
    title_y= "m(#mu#mu jjjj)" + title_y;
    directory = "mumu";
  }  
  else return;
  
  TString canvas = nameofhistogram + channel + sample;
  TString legend = nameofhistogram + channel + sample;
  canvas.Insert(0, "c_signal_VS");
  legend.Insert(0, "legend_");
  
  mapcanvas[canvas] = new TCanvas(canvas,"",800,800);
  gStyle -> SetOptStat(0);
  mapcanvas[canvas] -> SetTopMargin( 0.05 );
  mapcanvas[canvas] -> SetBottomMargin( 0.13 );
  mapcanvas[canvas] -> SetRightMargin( 0.05 );
  mapcanvas[canvas] -> SetLeftMargin( 0.16 );
  
  mapfunc[nameofhistogram + channel + sample + region] -> SetLineColor(kRed);
  mapfunc[nameofhistogram + channel + sample + region] -> GetYaxis()->SetLabelSize(0.05);
  mapfunc[nameofhistogram + channel + sample + region] -> GetYaxis()->SetTitle(title_y);
  mapfunc[nameofhistogram + channel + sample + region] -> GetYaxis()->SetRangeUser(0., 100.0);
  mapfunc[nameofhistogram + channel + sample + region] -> GetXaxis()->SetLabelSize(0.05);
  mapfunc[nameofhistogram + channel + sample + region] -> GetXaxis()->SetTitle(title_x);
  
  
  //  TF1 *func_eta = new TF1("func", "[0]*TMath::Gaus(x, [1], [2], 1)", -0.05, 0.3);
  map_fitfunc[nameofhistogram + channel + sample + region] = new TF1(nameofhistogram + channel + sample + region, "[1]", 0., 5500.);
  map_fitfunc[nameofhistogram + channel + sample + region]  -> SetLineColor(kBlue);

  
  mapfunc[nameofhistogram + channel + sample + region]  -> Fit(map_fitfunc[nameofhistogram + channel + sample + region], "R");
  double chi2, chi2_norm, Ndf, p0, p0_err;
  int ndf;
  
  
  chi2 = map_fitfunc[nameofhistogram + channel + sample + region] -> GetChisquare();
  ndf = map_fitfunc[nameofhistogram + channel + sample + region] -> GetNDF();
  Ndf = 0. + ndf;
  chi2_norm = chi2/ Ndf;
  p0 = map_fitfunc[nameofhistogram + channel + sample + region] -> GetParameter(1);
  p0_err = map_fitfunc[nameofhistogram + channel + sample + region] -> GetParError(1);
  
  TString chi2_Str, chi2_norm_str, Ndf_str, p0_str, p0_err_str;
  chi2_Str = Form("%4.2f", chi2);
  chi2_norm_str = Form("%4.3f", chi2_norm);
  Ndf_str = Form("%4.0f", Ndf);
  p0_str = Form("%4.3f", p0);
  p0_err_str = Form("%4.3f", p0_err);
  
  mapfunc[nameofhistogram + channel + sample + region] -> Draw();
  
  
  TLatex latex_CMSPriliminary, latex_Lumi, chi2_result, p0_result;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  chi2_result.SetNDC();
  p0_result.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Simulation}}}");
  latex_Lumi.SetTextSize(0.035);
  latex_Lumi.DrawLatex(0.7, 0.96, "35.9 fb^{-1} (13 TeV)");
  chi2_result.SetTextSize(0.035);
  chi2_result.DrawLatex(0.25, 0.7, "#chi^{2} / ndof = " + chi2_Str + "/" + Ndf_str + " = " + chi2_norm_str);
  p0_result.SetTextSize(0.035);
  p0_result.DrawLatex(0.25, 0.7 - 0.035, "ratio = " + p0_str + " #pm " + p0_err_str);

  TString pdfname;
  pdfname = "./emu_ratio/" + directory + "/"+ channel + "_emu_ratio_" + sample + ".pdf";
  
  mapcanvas[canvas] -> SaveAs(pdfname);

}

void run_plots_with_overflow(TString nameofhistogram, int N_bin, double binx[]){
  
  int N_regions = 3;
  int N_channels = 3;
  
  TString regions[3] = {"CR_Zmass", "CR_ttbar", "SR"};
  TString channels[3] = {"DiMu", "DiEle", "EMu"};
  
  
  for(int i = 0; i < N_regions; i++){
    for(int j = 0; j < N_channels; j++){
      make_hist_with_overflow(nameofhistogram + "_" + regions[i] + "_" + channels[j], "TTLL", N_bin, binx);
    }
  }
  
  
  for(int i = 0; i < N_regions; i++){
    make_ratio_plot(nameofhistogram, "TTLL", regions[i], regions[i] + "_DiMu");
    make_ratio_plot(nameofhistogram, "TTLL", regions[i], regions[i] + "_DiEle");
  }
  
  for(int i = 0; i < N_regions; i++){
    save_ratio_plot(nameofhistogram, "TTLL", regions[i], regions[i] + "_DiMu");
    save_ratio_plot(nameofhistogram, "TTLL", regions[i], regions[i] + "_DiEle");
  }

}




/// Main Function ////////////////////////////////////////////////////////////
void plot(){
  
  openfile("HN_pair_all_TTLL_powheg.root");
  
  double bin_lljjjjmass[15];
  int N_bin_lljjjjmass = 0;
  for(int i = 0; i < 11; i++){
    N_bin_lljjjjmass++;
    bin_lljjjjmass[i] = 200. * (i + 0.);
  }
  bin_lljjjjmass[11] = 2500;
  N_bin_lljjjjmass += 3;
  bin_lljjjjmass[12] = 3500;
  bin_lljjjjmass[13] = 5000;
  bin_lljjjjmass[14] = 5500;//overflow bin
  
  run_plots_with_overflow("mZp", N_bin_lljjjjmass, bin_lljjjjmass);
  
  

}// End of Main Function ////////////////////////////////////////////////////// 
