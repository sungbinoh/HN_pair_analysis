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
#include <TMath.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
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
#include <Math/DistFunc.h>
#include "setTDRStyle.C"

using namespace std;

map<TString, TH1*> maphist;
map<TString, TGraphAsymmErrors*> map_asym_gr;
map<TString, TFile*> mapfile;
map<TString, TCanvas*> mapcanvas;
map<TString, TPad*> mappad;
map<TString, THStack*> maphstack;
map<TString, TLegend*> maplegend;
map<TString, TH1F*> mapfunc;
map<TString, double> map_overflow;
map<TString, TLine*> mapline;
map<TString, TKey*> maphistcheck;
map<TString, TList*> maplist;
map<TString, std::vector<double>> map_bin_vector;

//cycle name
TString Cycle_name = "HN_pair_all";
//sample names
//TString data = "data_DoubleMuon";
TString bkgs = "bkgs";
TString prompt = "prompt";

TString DY_high = "DYJets";

const double alpha = 1 - 0.6827;

bool debug = true;

double DY_scale = 1.0798;

// -- Set TDR Style
void rootlogon()
{
  cout << "[rootlogon] get tdrstyle.C" << endl;
  // Load CMS style
  gROOT->LoadMacro("tdrstyle.C");
  setTDRStyle();
}

/// Getting Histogram Function ///////////////////////////////////////////////
TH1F * GetHist(TString hname){
  
  TH1F *h = NULL;
  std::map<TString, TH1F*>::iterator mapit = mapfunc.find(hname);
  if(mapit != mapfunc.end()) return mapit-> second;
  
  return h;
  
}
//////////////////////////////////////////////////////////////////////////////


/// Open ROOT file ///////////////////////////////////////////////////////////
void openfile(TString cyclename, TString samplename, TString dir, TString histname){
  
  TString filename = cyclename + "_" + samplename + ".root";
  //cout << "[[openfile ]]opening : " << filename << endl;

  TFile *current_file = new TFile ((filename)) ;
  
  //cout << "[[openfile ]] Cd : " << dir << endl;
  gDirectory->cd(dir);
      
  TH1F * current_hist = (TH1F*)gDirectory -> Get(histname);
  if(current_hist){
    current_hist -> SetDirectory(0);
  }
  TH1::AddDirectory(kFALSE);
  
  if(samplename.Contains("DYJets") && current_hist) current_hist->Scale(DY_scale);
  mapfunc[histname + cyclename + samplename] = current_hist;
  
  current_file -> Close();
  delete current_file;
}
//////////////////////////////////////////////////////////////////////////////


// -- Maing drawing function
void make_histogram(TString nameofhistogram, TString current_histname, int N_bin, double binx[]){

  if(debug){
    cout << "[[make_histogram]] checking binning arrary" << endl;
    for(int i = 0; i < N_bin; i++){
      cout << binx[i] << ", " << endl;
    }
  }

  TString title_y;
  if(nameofhistogram.Contains("N")) title_y = "Number";
  else title_y = "Events/bin";
  
  bool blind = false;
  blind = (nameofhistogram.Contains("SR")) && (!nameofhistogram.Contains("EMu"));
  
  TString pad1 = nameofhistogram;
  TString pad2 = nameofhistogram;
  TString canvas = nameofhistogram;
  TString hstack = nameofhistogram;
  TString legend = nameofhistogram;
  TString func = nameofhistogram;
  TString clone = nameofhistogram;
  TString line = nameofhistogram;
  pad1.Insert(0, "pad1_");
  pad2.Insert(0, "pad2_");
  canvas.Insert(0, "c_");
  hstack.Insert(0, "hs_");
  legend.Insert(0, "legend_");
  func.Insert(0, "ratio_");
  clone.Insert(0, "h3_");
  line.Insert(0, "l_");
  
  mapcanvas[canvas] = new TCanvas(canvas,"",800,800);
  gStyle -> SetOptStat(0);
  mapcanvas[canvas] -> SetTopMargin( 0.05 );
  mapcanvas[canvas] -> SetBottomMargin( 0.13 );
  mapcanvas[canvas] -> SetRightMargin( 0.05 );
  mapcanvas[canvas] -> SetLeftMargin( 0.16 );
  
  
  mappad[pad1] = new TPad("", "", 0, 0.25, 1, 1);
  mappad[pad1] -> SetTopMargin( 0.07 );
  mappad[pad1] -> SetBottomMargin( 0.02 );
  mappad[pad1] -> SetLeftMargin( 0.15 );
  mappad[pad1] -> SetRightMargin( 0.03 );
  mappad[pad1] -> Draw();
  mappad[pad1] -> cd();
  mappad[pad1] -> SetLogy();
  
  maplegend[legend] = new TLegend(0.60, 0.60, 0.96, 0.92);
  
  if(debug) cout << "1" << endl;
  
  int n_kind = 2;
  TString samples_array[] = {bkgs, prompt};
  Int_t colour_array[] = {419, 901};
  TString samples_legend[] = {"All bkgs", "Prompt Leptons (no CF)"};
  
  TString overflow = "overflow";
  Int_t nx_func    = GetHist(nameofhistogram + Cycle_name + bkgs) -> GetNbinsX()+1;
  Double_t x1_func = GetHist(nameofhistogram + Cycle_name + bkgs) -> GetBinLowEdge(1);
  Double_t bw_func = GetHist(nameofhistogram + Cycle_name + bkgs) -> GetBinWidth(nx_func);
  Double_t x2_func = GetHist(nameofhistogram + Cycle_name + bkgs) -> GetBinLowEdge(nx_func)+bw_func;
  binx[N_bin] = binx[N_bin - 1] + bw_func;
  if(debug) cout << "[[make_histogram]] binx[N_bin] ; " << binx[N_bin]  << endl;
  mapfunc[func] = new TH1F("", "", nx_func, x1_func, x2_func);
  if(debug) cout << "func rebin rebinning" << endl;
  mapfunc[func + "rebin"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "rebin", binx);
  mapfunc[func + "rebin"] -> Draw("hist");
  mapfunc[func + "rebin"] -> GetYaxis()->SetLabelSize(0.05);;
  mapfunc[func + "rebin"] -> GetYaxis()->SetTitleSize(0.07);;
  mapfunc[func + "rebin"] -> GetYaxis()->SetTitleOffset(1.02);;
  mapfunc[func + "rebin"] -> GetXaxis()->SetLabelSize(0);
  mapfunc[func + "rebin"] -> GetXaxis() -> SetRangeUser(x1_func, x2_func);
  mapfunc[func + "rebin"] -> GetXaxis() -> SetTitle(nameofhistogram);
  mapfunc[func + "rebin"] -> GetYaxis() -> SetTitle(title_y);

  mapfunc[func + "rebin" + bkgs] = (TH1F*)GetHist(nameofhistogram + Cycle_name + bkgs) -> Rebin(N_bin, func + "rebin", binx);
  mapfunc[func + "rebin" + bkgs] -> SetLineColor(kRed);
  //mapfunc[func + "rebin" + bkgs] -> Draw("histsame");
  
  double data_max = mapfunc[func + "rebin" + bkgs] -> GetMaximum();
  mapfunc[func + "rebin"] -> SetMaximum(data_max * 100000.);
  mapfunc[func + "rebin"] -> SetMinimum(data_max * 0.0001);

  mapfunc[func + "rebin"] -> Draw("hist");
  
  mapfunc[func + "rebin" + bkgs] -> Draw("histsame");

  
  mapfunc[func + "rebin" + prompt] = (TH1F*)GetHist(nameofhistogram + Cycle_name + prompt) -> Rebin(N_bin, func + "rebin", binx);
  mapfunc[func + "rebin" + prompt] -> SetLineColor(kBlue);
  mapfunc[func + "rebin" + prompt] -> Draw("histsame");
  
  mappad[pad1] -> Update();
  
  maplegend[legend] -> AddEntry(GetHist(func + "rebin" + bkgs), "All bkgs", "lf");
  maplegend[legend] -> AddEntry(GetHist(func + "rebin" + prompt), "Prompt Leptons (no CF)", "lf");
  
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kBlack);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(1001);
  maplegend[legend] -> SetShadowColor(0); // 0 = transparent
  
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> Draw("same");

  if(debug) cout << "6" << endl;

  mapcanvas[canvas] -> cd();

  mappad[pad2] = new TPad(pad2, "", 0, 0, 1, 0.25);
  mappad[pad2] -> SetTopMargin( 0.03 );
  mappad[pad2] -> SetBottomMargin( 0.4 );
  mappad[pad2] -> SetLeftMargin( 0.15 );
  mappad[pad2] -> SetRightMargin( 0.03 );
  mappad[pad2] -> Draw();
  mappad[pad2] -> cd();
  //mappad[pad2] -> SetLogy();
  
  mapfunc[clone] = (TH1F*)GetHist(func + "rebin" + bkgs) -> Clone(clone);
  mapfunc[clone] -> Add(GetHist(func + "rebin" + prompt), -1.);
  //mapfunc[clone] = (TH1F*)GetHist(func + "rebin" + prompt) -> Clone(clone);
  mapfunc[clone] -> Divide(mapfunc[func + "rebin" + bkgs] );
  mapfunc[clone] -> Draw("hist");
  TString name_x = nameofhistogram;
  mapfunc[clone] -> SetTitle("");
  //mapfunc[clone] -> SetYTitle("#frac{prompt.}{All bkgs.}");
  mapfunc[clone] -> SetYTitle("#frac{All - prompt}{All}");
  mapfunc[clone] -> SetFillColor(kOrange);
  mapfunc[clone] -> SetMarkerSize(20);
  mapfunc[clone] -> SetMarkerStyle(0);
  mapfunc[clone] -> SetLineColor(kBlack);
  mapfunc[clone] -> GetXaxis() -> SetTitle(name_x);
  mapfunc[clone] -> GetXaxis() -> SetLabelSize(0.10);
  mapfunc[clone] -> GetYaxis() -> SetLabelSize(0.08);
  mapfunc[clone] -> SetFillStyle(1001);
  mapfunc[clone] -> SetFillColorAlpha(45,0.35);
  mapfunc[clone] -> GetXaxis() -> SetRangeUser(x1_func, x2_func);
  mapfunc[clone] -> GetYaxis() -> SetTitleOffset(0.5);
  mapfunc[clone] -> GetXaxis() -> SetTitleSize(0.15);
  mapfunc[clone] -> GetYaxis() -> SetTitleSize(0.12);
  mapfunc[clone] -> SetMinimum(-2.0);
  mapfunc[clone] -> SetMaximum(2.0);
  mapfunc[clone] -> SetStats(0);
  mapfunc[clone] -> Draw("epsame");
  if(debug) cout << "8" << endl;
  
  mapcanvas[canvas] -> cd();
  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Simulation}}}");
  latex_Lumi.SetTextSize(0.035);
  latex_Lumi.DrawLatex(0.7, 0.96, "41.3 fb^{-1} (13 TeV)");
  
  if(debug) cout << "9" << endl;
  
  
  TString pdfname;
  pdfname = "./plots/ratio/";
  
  pdfname.Append(nameofhistogram);
  pdfname.Append("_ratio.pdf");
  if(debug) cout << "9.1" << endl;
  
  mapcanvas[canvas] -> SaveAs(pdfname);
  
  if(debug) cout << "10" << endl;
  
}

// -- Setup hist name and open files
void draw_histogram(TString histname){
  
  TString current_histname = "";
  
  TObjArray *token = histname.Tokenize("_");
  int N_token = token -> GetEntries();
  current_histname = ((TObjString *)(token->At(0)))->String();
  for(int i_token = 1; i_token < N_token - 1; i_token++){
    current_histname = current_histname + "_" + ((TObjString *)(token->At(i_token)))->String();
  }
  
  if(debug) cout << "[[draw_histogram]] histname : " << current_histname << endl;
  
  unsigned int N_bin = map_bin_vector[histname].size();
  double current_bins[100];
  for(unsigned int i_bin = 0; i_bin < N_bin; i_bin++){
    current_bins[i_bin] = map_bin_vector[histname].at(i_bin);
  }
  
  
  // -- Get which region(CR1, CR2 , ...) and channel (DiEle, DiMu, ...)
  int N_regions = 6;

  TString regions[] = {"CR_Zmass",
		       "CR_ttbar",
		       "CR_inv_mll",
		       "tight_CR_Zmass",
		       "tight_CR_ttbar",
		       "SR"
  };
  

  int N_channels = 3;
  TString channels[] = {"DiEle",
                        "DiMu",
                        "EMu"
  };
  
  const int N_syst = 19;
  TString systematics[N_syst] = {"central",
				 "ElectronScaleDown", "ElectronScaleUp",
				 "ElectronSmearDown", "ElectronSmearUp",
				 "JetsResDown", "JetsResUp",
				 "JetsScaleDown", "JetsScaleUp",
				 "SD_JMR_Down", "SD_JMR_Up",
                                 "SD_JMS_Down", "SD_JMS_Up",
				 "PUReweight_Down", "PUReweight_Up",
				 "MuonScaleDown", "MuonScaleUp",
                                 "MuonSmearDown", "MuonSmearUp",
				 //"PDFNormDown", "PDFNormUp"
  };
  
  TString current_dir = "empty";
  TString current_channel = "empty";
  
  for(int i_region = 0; i_region < N_regions; i_region++){
    if(histname.Contains(regions[i_region])) current_dir = regions[i_region]; 
  }
  for(int i_channel = 0; i_channel < N_channels; i_channel++){
    if(histname.Contains(channels[i_channel])) current_channel = channels[i_channel];
  }
  
  current_dir = current_dir + "_" + current_channel;
  
  if(debug) cout << "current_dir : " << current_dir << "/" << current_histname << "+systflags" << endl;
  if(current_dir.Contains("empty") || current_channel.Contains("empty")) return;

  for(int i_syst = 0; i_syst < N_syst; i_syst++){
    TString current_dir_syst = current_dir + "_" + systematics[i_syst];
    TString current_hist_syst = current_histname + "_" + systematics[i_syst];

    openfile(Cycle_name, bkgs, current_dir_syst, current_hist_syst);
    openfile(Cycle_name, prompt, current_dir_syst, current_hist_syst);
  }
  
  make_histogram(histname, current_histname, N_bin, current_bins);
  
}

// -- Open binning txt file
void open_binning_file(TString filename){

  cout << "[open_binning_file] start to open binngin file : " << filename << endl;
  
  ifstream data_file;
  data_file.open(filename);
  
  char line[500];
  if(data_file.is_open()){
    while(!data_file.eof()){
      data_file.getline(line, 500);
      if(debug) cout << line << endl;
      TString this_line = line;
      if(this_line.Contains("resign")) continue;
      
      TObjArray *tx = this_line.Tokenize("\t");
      int N_line_part = tx -> GetEntries();
      if(N_line_part < 1) continue;
      if(debug){
	cout << "[[open_binning_file]] histname : " << ((TObjString *)(tx->At(0)))->String() << endl;
	cout << "[[open_binning_file]] binning : " << ((TObjString *)(tx->At(1)))->String() << endl;
      }
      
      // -- saving current histname for mapping : TString && clear the map vector
      TString current_histname = ((TObjString *)(tx->At(0)))->String();
      if(debug) cout << "[[open_binning_file]] current_histname : " << current_histname << endl;
      map_bin_vector[current_histname].clear();
      
      // -- making vector<double> for binning
      TString current_binning = ((TObjString *)(tx->At(1)))->String();
      TObjArray *binning_token = current_binning.Tokenize(", ");
      for(int i_bin_part = 0; i_bin_part < binning_token->GetEntries() - 1; i_bin_part++){
	TString current_bin_value = ((TObjString *)(binning_token->At(i_bin_part)))->String();
	if(i_bin_part == 0){
	  current_bin_value.Remove(0,1);
	}
	if(debug) cout << current_bin_value << " | " << endl;
	double current_bin_double = current_bin_value.Atof();
	map_bin_vector[current_histname].push_back(current_bin_double);
      }
      
      draw_histogram(current_histname);
      map_bin_vector.clear();
      
    }//while not end of file
  }//end if file open

  data_file.close();
}
//////////////////////////////////////////////////////////////////////////////

/// Main nction ////////////////////////////////////////////////////////////
void plot(){
  rootlogon();
  open_binning_file("binning_SR.txt");
  //open_binning_file("binning_SS.txt");
  //open_binning_file("binning_Nvtx.txt");

}
