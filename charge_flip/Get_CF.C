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

TString cycle_name = "CF_MC_";
TString file_end = "_cat_v8-0-7.root";

TString dy_low = "SKDYJets_10to50";
TString dy_high = "SKDYJets";
TString dy_binned_50_120 = "SKZToEE_NNPDF30_13TeV-powheg_M_50_120";
TString dy_binned_120_200 = "SKZToEE_NNPDF30_13TeV-powheg_M_120_200";
TString dy_binned_200_400 = "SKZToEE_NNPDF30_13TeV-powheg_M_200_400";
TString dy_binned_400_800 = "SKZToEE_NNPDF30_13TeV-powheg_M_400_800";
TString dy_binned_800_1400 = "SKZToEE_NNPDF30_13TeV-powheg_M_800_1400";
TString dy_binned_1400_2300 = "SKZToEE_NNPDF30_13TeV-powheg_M_1400_2300";
TString dy_binned_2300_3500 = "SKZToEE_NNPDF30_13TeV-powheg_M_2300_3500";
TString dy_binned_3500_4500 = "SKZToEE_NNPDF30_13TeV-powheg_M_3500_4500";
TString dy_binned_4500_6000 = "SKZToEE_NNPDF30_13TeV-powheg_M_4500_6000";
TString dy_binned_6000_Inf= "SKZToEE_NNPDF30_13TeV-powheg_M_6000_Inf";


TString TTLL = "SKTTLL_powheg";
TString TT = "SKTT_powheg";

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
void openfile(TString sample){
  
  TString filename = cycle_name + sample + file_end;
  
  cout << "opening : " << filename << endl;
  mapfile[filename] = new TFile ((filename)) ;
  
  //cout << samplename << "'s histograms are" << endl;
  gDirectory -> Cd("Hists");
  
  TIter next(gDirectory->GetListOfKeys());
  TKey *key;
  vector<TString> histnames;
  while ((key = (TKey*)next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl->InheritsFrom("TH1")) continue;
    histnames.push_back(key -> GetName());
  }
  
  for(int i = 0; i < histnames.size(); i ++){
    maphist[histnames.at(i) + sample] = (TH1*)gDirectory -> Get(histnames.at(i));
  }
      
}
//////////////////////////////////////////////////////////////////////////////


void save_ratio_plot(TString nameofhistogram, TString sample){
  
  TString title_y = "Charge mis. measurement prob.";
  TString title_x = "1/P_{T} (1/GeV)";
  TString directory;
  
  TString canvas = nameofhistogram + sample;
  TString legend = nameofhistogram + sample;
  canvas.Insert(0, "c_signal_VS");
  legend.Insert(0, "legend_");
  
  mapcanvas[canvas] = new TCanvas(canvas,"",800,800);
  gStyle -> SetOptStat(0);
  mapcanvas[canvas] -> SetTopMargin( 0.05 );
  mapcanvas[canvas] -> SetBottomMargin( 0.13 );
  mapcanvas[canvas] -> SetRightMargin( 0.05 );
  mapcanvas[canvas] -> SetLeftMargin( 0.16 );
  
  mapfunc[nameofhistogram + "_num" + sample] = (TH1F*)GetHist(nameofhistogram + "_num" + sample) -> Clone(nameofhistogram + "_num" + sample);
  mapfunc[nameofhistogram + "_num" + sample] -> SetTitle("");
  mapfunc[nameofhistogram + "_num" + sample] -> Divide(maphist[nameofhistogram + "_den" + sample]);

  mapfunc[nameofhistogram + "_num" + sample] -> SetLineColor(kRed);
  mapfunc[nameofhistogram + "_num" + sample] -> GetYaxis()->SetLabelSize(0.02);
  mapfunc[nameofhistogram + "_num" + sample] -> GetYaxis()->SetTitle(title_y);
  mapfunc[nameofhistogram + "_num" + sample] -> GetXaxis()->SetLabelSize(0.02);
  mapfunc[nameofhistogram + "_num" + sample] -> GetXaxis()->SetRangeUser(0., 0.02);
  mapfunc[nameofhistogram + "_num" + sample] -> GetXaxis()->SetTitle(title_x);
  
  mapfunc[nameofhistogram + "_num" + sample] -> Draw();
  
  TLatex latex_CMSPriliminary, latex_Lumi, chi2_result, p0_result;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  chi2_result.SetNDC();
  p0_result.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.25, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Simulation}}}");
  latex_Lumi.SetTextSize(0.035);
  latex_Lumi.DrawLatex(0.7, 0.96, "35.9 fb^{-1} (13 TeV)");
  

  TString pdfname;
  pdfname = "./CF/CF_" + sample + "_" + nameofhistogram + ".pdf";
  
  mapcanvas[canvas] -> SaveAs(pdfname);

}

void add_all_MC(TString nameofhistogram){
  
  TString title_y = "Charge miss measure. prob.";
  TString title_x = "1/P_{T} (1/GeV)";
  TString directory;

  TString canvas = nameofhistogram;
  TString legend = nameofhistogram;
  canvas.Insert(0, "c_signal_VS");
  legend.Insert(0, "legend_");

  mapcanvas[canvas] = new TCanvas(canvas,"",800,800);
  gStyle -> SetOptStat(0);
  mapcanvas[canvas] -> SetTopMargin( 0.05 );
  mapcanvas[canvas] -> SetBottomMargin( 0.13 );
  mapcanvas[canvas] -> SetRightMargin( 0.05 );
  mapcanvas[canvas] -> SetLeftMargin( 0.16 );

  mapfunc[nameofhistogram + "_num"] = (TH1F*)GetHist(nameofhistogram + "_num" + dy_low) -> Clone(nameofhistogram + "_num");
  mapfunc[nameofhistogram + "_num"] -> Add(GetHist(nameofhistogram + "_num" + dy_high));
  mapfunc[nameofhistogram + "_num"] -> Add(GetHist(nameofhistogram + "_num" + dy_binned_50_120));
  mapfunc[nameofhistogram + "_num"] -> Add(GetHist(nameofhistogram + "_num" + dy_binned_120_200));
  mapfunc[nameofhistogram + "_num"] -> Add(GetHist(nameofhistogram + "_num" + dy_binned_200_400));
  mapfunc[nameofhistogram + "_num"] -> Add(GetHist(nameofhistogram + "_num" + dy_binned_400_800));
  mapfunc[nameofhistogram + "_num"] -> Add(GetHist(nameofhistogram + "_num" + dy_binned_800_1400));
  mapfunc[nameofhistogram + "_num"] -> Add(GetHist(nameofhistogram + "_num" + dy_binned_1400_2300));
  mapfunc[nameofhistogram + "_num"] -> Add(GetHist(nameofhistogram + "_num" + dy_binned_2300_3500));
  mapfunc[nameofhistogram + "_num"] -> Add(GetHist(nameofhistogram + "_num" + dy_binned_3500_4500));
  mapfunc[nameofhistogram + "_num"] -> Add(GetHist(nameofhistogram + "_num" + dy_binned_4500_6000));
  mapfunc[nameofhistogram + "_num"] -> Add(GetHist(nameofhistogram + "_num" + dy_binned_6000_Inf));
  mapfunc[nameofhistogram + "_num"] -> Add(GetHist(nameofhistogram + "_num" + TT));
  mapfunc[nameofhistogram + "_num"] -> Add(GetHist(nameofhistogram + "_num" + TTLL));

  mapfunc[nameofhistogram + "_num"] -> SetTitle("");

  mapfunc[nameofhistogram + "_den"] = (TH1F*)GetHist(nameofhistogram + "_den" + dy_low) -> Clone(nameofhistogram + "_den");
  mapfunc[nameofhistogram + "_den"] -> Add(GetHist(nameofhistogram + "_den" + dy_high));
  mapfunc[nameofhistogram + "_den"] -> Add(GetHist(nameofhistogram + "_den" + dy_binned_50_120));
  mapfunc[nameofhistogram + "_den"] -> Add(GetHist(nameofhistogram + "_den" + dy_binned_120_200));
  mapfunc[nameofhistogram + "_den"] -> Add(GetHist(nameofhistogram + "_den" + dy_binned_200_400));
  mapfunc[nameofhistogram + "_den"] -> Add(GetHist(nameofhistogram + "_den" + dy_binned_400_800));
  mapfunc[nameofhistogram + "_den"] -> Add(GetHist(nameofhistogram + "_den" + dy_binned_800_1400));
  mapfunc[nameofhistogram + "_den"] -> Add(GetHist(nameofhistogram + "_den" + dy_binned_1400_2300));
  mapfunc[nameofhistogram + "_den"] -> Add(GetHist(nameofhistogram + "_den" + dy_binned_2300_3500));
  mapfunc[nameofhistogram + "_den"] -> Add(GetHist(nameofhistogram + "_den" + dy_binned_3500_4500));
  mapfunc[nameofhistogram + "_den"] -> Add(GetHist(nameofhistogram + "_den" + dy_binned_4500_6000));
  mapfunc[nameofhistogram + "_den"] -> Add(GetHist(nameofhistogram + "_den" + dy_binned_6000_Inf));
  mapfunc[nameofhistogram + "_den"] -> Add(GetHist(nameofhistogram + "_den" + TT));
  mapfunc[nameofhistogram + "_den"] -> Add(GetHist(nameofhistogram + "_den" + TTLL));
  
  mapfunc[nameofhistogram + "_num"] -> Divide(mapfunc[nameofhistogram + "_den"]);
  
  map_fitfunc[nameofhistogram] = new TF1(nameofhistogram, "[0] + [1]*x", 0., 0.02);
  map_fitfunc[nameofhistogram]  -> SetLineColor(kBlue);
  
  mapfunc[nameofhistogram + "_num"]  -> Fit(map_fitfunc[nameofhistogram], "R");

  double chi2, chi2_norm, Ndf, p0, p0_err, p1, p1_err;
  int ndf;
  
  chi2 = map_fitfunc[nameofhistogram] -> GetChisquare();
  ndf = map_fitfunc[nameofhistogram] -> GetNDF();
  Ndf = 0. + ndf;
  chi2_norm = chi2/ Ndf;
  p0 = map_fitfunc[nameofhistogram] -> GetParameter(0);
  p0_err = map_fitfunc[nameofhistogram] -> GetParError(0);
  p1 = map_fitfunc[nameofhistogram] -> GetParameter(1);
  p1_err = map_fitfunc[nameofhistogram] -> GetParError(1);

  TString chi2_Str, chi2_norm_str, Ndf_str, p0_str, p0_err_str, p1_str, p1_err_str;
  chi2_Str = Form("%4.2f", chi2);
  chi2_norm_str = Form("%4.3f", chi2_norm);
  Ndf_str = Form("%4.0f", Ndf);
  p0_str = Form("%4.5f", p0);
  p0_err_str = Form("%4.5f", p0_err);
  p1_str = Form("%4.5f", p1);
  p1_err_str = Form("%4.5f", p1_err);
  

  mapfunc[nameofhistogram + "_num"] -> SetLineColor(kRed);
  mapfunc[nameofhistogram + "_num"] -> GetYaxis()->SetLabelSize(0.02);
  mapfunc[nameofhistogram + "_num"] -> GetYaxis()->SetTitle(title_y);
  //mapfunc[nameofhistogram + "_num" + sample] -> GetYaxis()->SetRangeUser(0., 2.0);                                                                                                                                                                                            
  mapfunc[nameofhistogram + "_num"] -> GetXaxis()->SetLabelSize(0.02);
  mapfunc[nameofhistogram + "_num"] -> GetXaxis()->SetRangeUser(0., 0.02);
  mapfunc[nameofhistogram + "_num"] -> GetXaxis()->SetTitle(title_x);


  mapfunc[nameofhistogram + "_num"] -> Draw();

  TLatex latex_CMSPriliminary, latex_Lumi, chi2_result, p0_result, p1_result;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  chi2_result.SetNDC();
  p0_result.SetNDC();
  p1_result.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.25, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Simulation}}}");
  latex_Lumi.SetTextSize(0.035);
  latex_Lumi.DrawLatex(0.7, 0.96, "35.9 fb^{-1} (13 TeV)");
  chi2_result.SetTextSize(0.025);
  chi2_result.DrawLatex(0.35, 0.7, "#chi^{2} / ndof = " + chi2_Str + "/" + Ndf_str + " = " + chi2_norm_str);
  p0_result.SetTextSize(0.025);
  p0_result.DrawLatex(0.35, 0.7 - 0.025, "y = a + bx, where");
  p1_result.SetTextSize(0.025);
  p1_result.DrawLatex(0.35, 0.7 - 0.025 - 0.025, "a = " + p0_str + " #pm " + p0_err_str + ", b = " + p1_str + " #pm " + p1_err_str );
  

  TString pdfname;
  pdfname = "./CF/CF_all_" + nameofhistogram + ".pdf";

  mapcanvas[canvas] -> SaveAs(pdfname);
}

void run_for_sample(TString sample){
  save_ratio_plot("CF_RECOpT_EC", sample);
  save_ratio_plot("CF_RECOpT_IB", sample);
  save_ratio_plot("CF_RECOpT_OB", sample);
  
  save_ratio_plot("CF_TRUEpT_EC", sample);
  save_ratio_plot("CF_TRUEpT_IB", sample);
  save_ratio_plot("CF_TRUEpT_OB", sample);
}



/// Main Function ////////////////////////////////////////////////////////////
void plot(){
   
  openfile(dy_low);
  openfile(dy_high);
  openfile(dy_binned_50_120);
  openfile(dy_binned_120_200);
  openfile(dy_binned_200_400);
  openfile(dy_binned_400_800);
  openfile(dy_binned_800_1400);
  openfile(dy_binned_1400_2300);
  openfile(dy_binned_2300_3500);
  openfile(dy_binned_3500_4500);
  openfile(dy_binned_4500_6000);
  openfile(dy_binned_6000_Inf);
  openfile(TTLL);
  openfile(TT);

  run_for_sample(dy_low);
  run_for_sample(dy_high);
  run_for_sample(dy_binned_50_120);
  run_for_sample(dy_binned_120_200);
  run_for_sample(dy_binned_200_400);
  run_for_sample(dy_binned_400_800);
  run_for_sample(dy_binned_800_1400);
  run_for_sample(dy_binned_1400_2300);
  run_for_sample(dy_binned_2300_3500);
  run_for_sample(dy_binned_3500_4500);
  run_for_sample(dy_binned_4500_6000);
  run_for_sample(dy_binned_6000_Inf);
  run_for_sample(TTLL);
  run_for_sample(TT);
  
  add_all_MC("CF_RECOpT_EC");
  add_all_MC("CF_RECOpT_IB");
  add_all_MC("CF_RECOpT_OB");
  add_all_MC("CF_TRUEpT_EC");
  add_all_MC("CF_TRUEpT_IB");
  add_all_MC("CF_TRUEpT_OB");

  

}// End of Main Function ////////////////////////////////////////////////////// 
