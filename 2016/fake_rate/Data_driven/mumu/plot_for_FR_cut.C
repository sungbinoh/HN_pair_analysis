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
map<TString, TFile*> mapfile;
map<TString, TCanvas*> mapcanvas;
map<TString, TPad*> mappad;
map<TString, THStack*> maphstack;
map<TString, TLegend*> maplegend;
map<TString, TH1F*> mapfunc;
map<TString, TLine*> mapline;
map<TString, TKey*> maphistcheck;
map<TString, TList*> maplist;

//cycle name
//TString Cycle_name = "ExampleAnalyzerElectronMuon";
//TString Cycle_name = "ExampleAnalyzerDiMuon";
TString Cycle_name = "Muon_FR_cal_all";
//sample names                          
//TString data = "data_DoubleMuon";
TString data = "data_SingleMuon";

//TString periodC = "periodC_SKMuonEG";
//TString periodD = "periodD_SKMuonEG";

//DY
TString DY_low = "SKDYJets_10to50";
TString DY_high = "SKDYJets";

TString Wjets = "SKWJets";

//single top
TString SingleTbar_tW = "SKSingleTbar_tW";
TString SingleTbar_t = "SKSingleTbar_t";
TString SingleTop_s = "SKSingleTop_s";
TString SingleTop_tW = "SKSingleTop_tW";
TString SingleTop_t = "SKSingleTop_t";

TString WGtoLNuG = "SKWGtoLNuG";

//VV
TString WW = "SKWW";
TString WZ = "SKWZ";
TString ZZ = "SKZZ";

//TString ttbar = "SKTTJets_aMC";
TString ttbar = "SKTT_powheg";
TString ZGto2LG = "SKZGto2LG";


 /// Getting Histogram Function ///////////////////////////////////////////////
 TH1 * GetHist(TString hname){

   TH1 *h = NULL;
   std::map<TString, TH1*>::iterator mapit = maphist.find(hname);
   if(mapit != maphist.end()) return mapit-> second;

   return h;

 }
 //////////////////////////////////////////////////////////////////////////////


 /// Open ROOT file ///////////////////////////////////////////////////////////
 void openfile(TString cyclename, TString samplename){
  
   TString underbar = "_";
   TString version = "cat_v8-0-7.root";
   TString filename = cyclename + underbar + samplename + underbar + version;
   if(samplename == data) filename = cyclename + underbar + samplename + "_cat_v8-0-7.root";

   cout << "opening : " << filename << endl;
   
   mapfile[filename] = new TFile ((filename)) ;
   
   //cout << samplename << "'s histograms are" << endl;
   gDirectory->cd("Hists");
   TIter next(gDirectory->GetListOfKeys());
   TKey *key;
   vector<TString> histnames;
   while ((key = (TKey*)next())) {
     TClass *cl = gROOT->GetClass(key->GetClassName());
     if (!cl->InheritsFrom("TH1")) continue;
     histnames.push_back(key -> GetName());
   }
   
   for(int i = 0; i < histnames.size(); i ++){
     //cout << histnames.at(i) << endl;
     maphist[histnames.at(i) + cyclename + samplename] = (TH1*)gDirectory -> Get(histnames.at(i));
     //cout << histnames.at(i) + cyclename + samplename << endl;
   }
   
   gDirectory -> cd();
   
 }
//////////////////////////////////////////////////////////////////////////////


 /// Make Output Histogram ////////////////////////////////////////////////////
void makehistogram(TString nameofhistogram, float xmin, float xmax, float rebin, float ymax, TString name_x, bool name_y){
  
  TString title_y;
  if(name_y) title_y = "Events/" + TString::Itoa(rebin,10) + " (GeV)";
  else title_y = "Number";
  
  
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
   gStyle -> SetOptStat(1111);
   mapcanvas[canvas] -> SetTopMargin( 0.05 );
   mapcanvas[canvas] -> SetBottomMargin( 0.13 );
   mapcanvas[canvas] -> SetRightMargin( 0.05 );
   mapcanvas[canvas] -> SetLeftMargin( 0.16 );
   //gStyle -> SetPadTickY(1);
   //gStyle -> SetPadTickX(1);
      
   mappad[pad1] = new TPad("", "", 0, 0.25, 1, 1);
   mappad[pad1] -> SetTopMargin( 0.07 );
   mappad[pad1] -> SetBottomMargin( 0.02 );
   mappad[pad1] -> SetLeftMargin( 0.15 );
   mappad[pad1] -> SetRightMargin( 0.03 );
   //gStyle -> SetPadTickY(1);
   //gStyle -> SetPadTickX(1);
   mappad[pad1] -> Draw();
   mappad[pad1] -> cd();
   mappad[pad1] -> SetLogy();

   maplegend[legend] = new TLegend(0.54, 0.6, 0.81, 0.92);
   cout << "0" << endl;
   //cout << nameofhistogram + Cycle_name + data << endl;
   //mapfunc[func] = new TH1F("", "", GetHist(nameofhistogram + Cycle_name + periodC) -> GetNbinsX(),  GetHist(nameofhistogram + Cycle_name + periodC) -> GetXaxis() -> GetXmin(), GetHist(nameofhistogram + Cycle_name + periodC) -> GetXaxis() -> GetXmax());
   mapfunc[func] = new TH1F("", "", GetHist(nameofhistogram + Cycle_name + data) -> GetNbinsX(),  GetHist(nameofhistogram + Cycle_name + data) -> GetXaxis() -> GetXmin(), GetHist(nameofhistogram + Cycle_name + data) -> GetXaxis() -> GetXmax());

   cout << "1" << endl;
   
   maphstack[hstack] = new THStack(hstack, "Stacked_" + nameofhistogram);
   gStyle->SetOptTitle(0);
   //maphstack[hstack] -> SetTitle("CMS simulation");
   //gStyle->SetTitleFont(32,"t");
   //gStyle->SetTitleH(0.03);
   //gStyle->SetTitleX(0.83);
   //gStyle->SetTitleY(0.97);
   //gStyle->SetTitleH(0.1);//title perpendicular size  
   //gStyle->SetTitleW(0.9);//title width
   //gStyle->SetTitleX(.78);//title x-axis offset
   
   
   ///To do for each sample
   ///1. add, fill color, markerstyle
   ///2. add to stack
   ///3. set legend
   mapfunc[func] -> Rebin(rebin);
   
   int n_kind = 4;
   //TString samples_array[] = {DY_high, ttbar, Wjets, WGtoLNuG, WZ};
   TString samples_array[] = {WZ, ttbar, DY_high, Wjets};
   
   //kWhite  = 0,   kBlack  = 1,   kGray    = 920,  kRed    = 632,  kGreen  = 416, kBlue   = 600, kYellow = 400, kMagenta = 616,  kCyan   = 432,  kOrange = 800,
   //kSpring = 820, kTeal   = 840, kAzure   =  860, kViolet = 880,  kPink   = 900
   //Int_t colour_array[] = {906, 416, 867, 600, 800, 432, 632};
   //Int_t colour_array[] = {400, 632, 416, 813, 819};
   Int_t colour_array[] = {819, 632, 400, 416};
   //TString samples_legend[] = {"DY", "top", "Wjets", "X#gamma", "VV"};
   TString samples_legend[] = {"VV", "top", "DY", "Wjets"};

   //std::vector<TString> samples(samples_array, samples_array + n_kind);
   //std::vector<Int_t> fill_colour(colour_array, colour_array + n_kind);
   //std::vector<TString> legends(samples_legend, samples_legend + n_kind);
   
   cout << "2" << endl;
   
   
   
   TString name_cycle = nameofhistogram + Cycle_name;
   cout << "check1" << endl;
   if(maphist[name_cycle + DY_high] || maphist[name_cycle + DY_low]){
     GetHist(name_cycle + DY_high) -> Add(GetHist(name_cycle + DY_low));
   }
   if(maphist[name_cycle + SingleTop_s] || maphist[name_cycle + SingleTbar_tW] || maphist[name_cycle + SingleTbar_t] || maphist[name_cycle + SingleTop_tW] || maphist[name_cycle + SingleTop_t]){
     GetHist(name_cycle +  SingleTop_s) -> Add(GetHist(name_cycle + SingleTbar_tW));
     GetHist(name_cycle +  SingleTop_s) -> Add(GetHist(name_cycle + SingleTbar_t));
     GetHist(name_cycle +  SingleTop_s) -> Add(GetHist(name_cycle + SingleTop_tW));
     GetHist(name_cycle +  SingleTop_s) -> Add(GetHist(name_cycle + SingleTop_t));
   }
   /*
   if(maphist[name_cycle + WGtoLNuG] || maphist[name_cycle + ZGto2LG]){
     GetHist(name_cycle +  WGtoLNuG) -> Add(GetHist(name_cycle + ZGto2LG));
   }
   */
   if(maphist[name_cycle + WW] || maphist[name_cycle + WZ] || maphist[name_cycle + ZZ]){
     GetHist(name_cycle +  WZ) -> Add(GetHist(name_cycle + WW));
     GetHist(name_cycle +  WZ) -> Add(GetHist(name_cycle + ZZ));
   }
   GetHist(name_cycle +  ttbar) -> Add(GetHist(name_cycle + SingleTop_s));
   
   cout << "3" << endl;


   cout << "check2" << endl;

   
   maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + data), "data", "p");
   
   for(int i = 0; i < n_kind; i++){
     cout << samples_array[i] << endl;
     if(maphist[nameofhistogram + Cycle_name + samples_array[i]]){
       GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> SetFillColor(colour_array[i]);
       GetHist(nameofhistogram + Cycle_name + samples_array[i]) -> Rebin(rebin);
       cout << samples_array[i] << endl;
       maphstack[hstack] -> Add(GetHist(nameofhistogram + Cycle_name + samples_array[i]));
       maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + samples_array[n_kind - 1 - i]), samples_legend[n_kind - 1 - i], "f");
       mapfunc[func] -> Add(GetHist(nameofhistogram + Cycle_name + samples_array[i]));
     }
   }//for loop

   cout << "4" << endl;


   /*
   for(int i = 0; i < n_kind; i++){
     maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + samples_array[n_kind - i]), samples_legend[n_kind - i], "f");
     cout << "n_kind - i = " << n_kind - i << endl;
   }
   */


   ////////////////////////////////////////////////////////////////////////////////////////////////////////////

   maphstack[hstack] -> Draw("hist");
   maphstack[hstack] -> GetYaxis()->SetLabelSize(0.05);;
   maphstack[hstack] -> GetYaxis()->SetTitleSize(0.07);;
   maphstack[hstack] -> GetYaxis()->SetTitleOffset(1.02);;
   maphstack[hstack] -> GetXaxis()->SetLabelSize(0);;
   maphstack[hstack] -> GetXaxis() -> SetRangeUser(xmin, xmax);
   //maphstack[hstack] -> SetMaximum(ymax);
   double ymax_data = GetHist(nameofhistogram + Cycle_name + data) ->GetBinContent(GetHist(nameofhistogram + Cycle_name + data) -> GetMaximumBin());
   maphstack[hstack] -> SetMaximum(ymax_data * rebin * 1000.);
   maphstack[hstack] -> SetMinimum(1.0);
   //maphstack[hstack] -> SetMinimum(0.001); 
   maphstack[hstack] -> GetXaxis() -> SetTitle(nameofhistogram);
   //maphstack[hstack] -> GetYaxis() -> SetTitle("Entries / 25 GeV");
   maphstack[hstack] -> GetYaxis() -> SetTitle(title_y);  
   maphstack[hstack] -> Draw("histsame");

   GetHist(nameofhistogram + Cycle_name + data) -> Rebin(rebin);
   GetHist(nameofhistogram + Cycle_name + data) -> SetMarkerStyle(20);
   GetHist(nameofhistogram + Cycle_name + data) -> SetMarkerSize(1.0);
   GetHist(nameofhistogram + Cycle_name + data) -> Draw("epsame");
   GetHist(nameofhistogram + Cycle_name + data) -> GetXaxis() -> SetRangeUser(xmin, xmax);  
   GetHist(nameofhistogram + Cycle_name + data) -> SetMinimum(1.0);
   
   mappad[pad1] -> Update();
   
   GetHist(nameofhistogram + Cycle_name + data) -> Draw("epsame");       

   maplegend[legend] -> SetFillColor(kWhite);
   maplegend[legend] -> SetLineColor(kWhite);
   maplegend[legend] -> SetBorderSize(1);
   maplegend[legend] -> SetFillStyle(1001);
   maplegend[legend] -> SetShadowColor(0); // 0 = transparent
   maplegend[legend] -> SetEntrySeparation(0.3);
   maplegend[legend] -> Draw("same");
   
   //maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + data), "data", "p");
   
   mapcanvas[canvas] -> cd();

   mappad[pad2] = new TPad(pad2, "", 0, 0, 1, 0.25);
   mappad[pad2] -> SetTopMargin( 0.03 );
   mappad[pad2] -> SetBottomMargin( 0.4 );
   mappad[pad2] -> SetLeftMargin( 0.15 );
   mappad[pad2] -> SetRightMargin( 0.03 );
   mappad[pad2] -> Draw();
   mappad[pad2] -> cd();

   mapfunc[clone] = (TH1F*)GetHist(nameofhistogram + Cycle_name + data) -> Clone(clone);
   //mapfunc[clone] = (TH1F*)GetHist(nameofhistogram + Cycle_name + data) -> Clone(clone);

   mapfunc["stat" + nameofhistogram] = (TH1F*)GetHist(nameofhistogram + Cycle_name + data) -> Clone(clone);

   mapline[line] = new TLine(xmin, 1, xmax, 1);
   mapline[line] -> Draw();
   mapline[line] -> SetLineStyle(1);
   mapline[line] -> SetLineColor(kRed);

   mapfunc["stat" + nameofhistogram] -> SetTitle("");
   mapfunc["stat" + nameofhistogram] -> SetYTitle("#frac{Obs.}{Pred.}");
   mapfunc["stat" + nameofhistogram] -> SetFillColor(kOrange);
   mapfunc["stat" + nameofhistogram] -> SetMarkerSize(0);
   mapfunc["stat" + nameofhistogram] -> SetMarkerStyle(0);
   mapfunc["stat" + nameofhistogram] -> SetLineColor(kWhite);
   mapfunc["stat" + nameofhistogram] -> GetXaxis() -> SetTitle(name_x);
   mapfunc["stat" + nameofhistogram] -> GetXaxis() -> SetLabelSize(0.10);
   mapfunc["stat" + nameofhistogram] -> GetYaxis() -> SetLabelSize(0.08);
   mapfunc["stat" + nameofhistogram] -> SetFillColorAlpha(45,0.35);
   mapfunc["stat" + nameofhistogram] -> GetXaxis() -> SetRangeUser(xmin, xmax);
   mapfunc["stat" + nameofhistogram] -> GetYaxis() -> SetTitleOffset(0.5);
   mapfunc["stat" + nameofhistogram] -> GetXaxis() -> SetTitleSize(0.15);
   mapfunc["stat" + nameofhistogram] -> GetYaxis() -> SetTitleSize(0.12);
   mapfunc["stat" + nameofhistogram] -> SetMinimum(0.6);
   mapfunc["stat" + nameofhistogram] -> SetMaximum(5.0);
   if(nameofhistogram.Contains("loose")) mapfunc["stat" + nameofhistogram] -> SetMaximum(10.0);
   
   mapfunc["stat" + nameofhistogram] -> SetStats(0);
   cout << "data integral , MC integral" << endl;
   cout << mapfunc["stat" + nameofhistogram] -> Integral() << endl;
   cout << mapfunc[func] -> Integral() << endl;

   mapfunc["stat" + nameofhistogram] -> Divide(mapfunc[func]);
   Int_t ncells = mapfunc["stat" + nameofhistogram] -> GetSize();
   for(int i = 0; i < ncells - 2; i++){
     mapfunc["stat" + nameofhistogram] -> SetBinContent(i, 1.);
     mapfunc["stat" + nameofhistogram] -> SetBinError(i, 0.13);
   }
   mapfunc["stat" + nameofhistogram] -> Draw("CE2");
   mapfunc[clone] -> Divide(mapfunc[func]);
   mapfunc[clone] ->Draw("PE1same");
   mapline[line] -> Draw("same");

   maplegend["bottom" + legend] = new TLegend(0.2, 0.8, 0.6, 0.9);
   maplegend["bottom" + legend]->SetBorderSize(0);
   maplegend["bottom" + legend]->SetNColumns(2);
   maplegend["bottom" + legend]->AddEntry(mapfunc["stat" + nameofhistogram], "13% Band", "f");
   maplegend["bottom" + legend]->AddEntry(mapfunc[clone] , "Obs./Pred.", "p");
   maplegend["bottom" + legend]->Draw("same");
   
   mapline[line] -> Draw("same");

   
   mapcanvas[canvas] -> cd();
   TLatex latex_CMSPriliminary, latex_Lumi;
   latex_CMSPriliminary.SetNDC();
   latex_Lumi.SetNDC();
   latex_CMSPriliminary.SetTextSize(0.035);
   latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
   latex_Lumi.SetTextSize(0.035);
   latex_Lumi.DrawLatex(0.7, 0.96, "35.9 fb^{-1} (13 TeV)");
   
   TString pdfname = "./pdfs/";
   pdfname.Append(nameofhistogram);
   pdfname.Append(".pdf");
   mapcanvas[canvas] -> SaveAs(pdfname);

 }
 //////////////////////////////////////////////////////////////////////////////

void draw_histogram(TString nameofhistogram, float xmin, float xmax, float rebin, float ymax, TString name_x, bool name_y){

  vector<TString> directories;
  //directories.push_back("basic_PU");
  //directories.push_back("basic_tempPU");
  //directories.push_back("Fake_CR_1");

  makehistogram(nameofhistogram, xmin, xmax, rebin, ymax, name_x, name_y);
  //makehistogram(nameofhistogram + "_Nvtx_reweight_central", xmin, xmax, rebin, ymax, name_x, name_y);
  //makehistogram(nameofhistogram + "_Nvtx_reweight", xmin, xmax, rebin, ymax, name_x, name_y);
  
}


/// Main Function ////////////////////////////////////////////////////////////
void plot(){
    
  openfile(Cycle_name, data);
  openfile(Cycle_name, DY_low);
  openfile(Cycle_name, DY_high);
  openfile(Cycle_name, Wjets);
  openfile(Cycle_name, SingleTbar_tW);
  openfile(Cycle_name, SingleTbar_t);
  openfile(Cycle_name, SingleTop_s);
  openfile(Cycle_name, SingleTop_tW);
  openfile(Cycle_name, SingleTop_t);
  openfile(Cycle_name, WW);
  openfile(Cycle_name, WZ);
  openfile(Cycle_name, ZZ);
  openfile(Cycle_name, ttbar);

  cout << "open files complete" << endl;

  //draw_histogram("SingleMuonTrigger_Dijet_Awayjet_40_all_MUON_HN_LOOSEv7_SIP3_pt_cone_F0", 0., 200., 1., 5000., "P_{T}^{cone}(1st #mu) (GeV)", true);
  //draw_histogram("SingleMuonTrigger_Dijet_Awayjet_40_all_MUON_HN_LOOSEv7_SIP3_Nvtx_reweight_pt_cone_F0", 0., 200., 1., 5000., "P_{T}^{cone}(1st #mu) (GeV)", true);
  
  draw_histogram("Ptcone_UseEvent_loose", 0., 200., 5., 5000., "P_{T}^{cone}(1st #mu) (GeV)", true);
  draw_histogram("Ptcone_UseEvent_tight", 0., 200., 5., 200., "P_{T}^{cone}(1st #mu) (GeV)", true);
  draw_histogram("Eta_UseEvent_loose", -3., 3., 10., 5000., "#eta (#mu)", false);
  draw_histogram("Eta_UseEvent_tight", -3., 3., 10., 200., "#eta (#mu)", false);
  
  /*
  draw_histogram("HLT_Mu17_TrkIsoVVL_v_ZPeak_mll", 70., 110., 1., 20000., "M(#mu#mu) (GeV)", true);
  draw_histogram("HLT_Mu8_TrkIsoVVL_v_ZPeak_mll", 70., 110., 1., 1000., "M(#mu#mu) (GeV)", true);
  draw_histogram("HLT_Mu3_PFJet40_v_ZPeak_mll", 70., 110., 1., 60., "M(#mu#mu) (GeV)", true);
  
  draw_histogram("HLT_Mu3_PFJet40_v_ZPeak_leadpt", 0., 200., 1., 20., "P_{T}(1st #mu) (GeV)", true);
  draw_histogram("HLT_Mu8_TrkIsoVVL_v_ZPeak_leadpt", 0., 200., 1., 400., "P_{T}(1st #mu) (GeV)", true);
  draw_histogram("HLT_Mu17_TrkIsoVVL_v_ZPeak_leadpt", 0., 200., 1., 8000., "P_{T}(1st #mu) (GeV)", true);

  draw_histogram("HLT_Mu17_TrkIsoVVL_v_ZPeak_PFMET", 0., 200., 5., 30000., "MET (GeV)", true);
  draw_histogram("HLT_Mu8_TrkIsoVVL_v_ZPeak_PFMET", 0., 200., 5., 2000., "MET (GeV)", true);
  draw_histogram("HLT_Mu3_PFJet40_v_ZPeak_PFMET", 0., 200., 5., 60., "MET (GeV)", true);
  
  draw_histogram("HLT_Mu17_TrkIsoVVL_v_ZPeak_Nvtx", 0., 50., 1., 8000., "N(vertice)", false);
  draw_histogram("HLT_Mu8_TrkIsoVVL_v_ZPeak_Nvtx", 0., 50., 1., 4000., "N(vertice)", false);
  draw_histogram("HLT_Mu3_PFJet40_v_ZPeak_Nvtx", 0., 50., 1., 1000., "N(vertice)", false);
    
  
  //John's W dominant region
  draw_histogram("HLT_Mu17_TrkIsoVVL_v_W_John_MT", 0., 200., 1., 30000., "M_{T}(#mu, METvec) (GeV)", true);
  draw_histogram("HLT_Mu8_TrkIsoVVL_v_W_John_MT", 0., 200., 1., 2000., "M_{T}(#mu, METvec) (GeV)", true);
  draw_histogram("HLT_Mu3_PFJet40_v_W_John_MT", 0., 200., 1., 60., "M_{T}(#mu, METvec) (GeV)", true);
  draw_histogram("HLT_Mu17_TrkIsoVVL_v_W_John_PFMET", 0., 200., 5., 200000., "MET (GeV)", true);
  draw_histogram("HLT_Mu8_TrkIsoVVL_v_W_John_PFMET", 0., 200., 5., 10000., "MET (GeV)", true);
  draw_histogram("HLT_Mu3_PFJet40_v_W_John_PFMET", 0., 200., 5., 1000., "MET (GeV)", true);
  draw_histogram("HLT_Mu17_TrkIsoVVL_v_W_John_leadpt", 0., 200., 1., 50000., "P_{T}(1st #mu) (GeV)", true);
  draw_histogram("HLT_Mu8_TrkIsoVVL_v_W_John_leadpt", 0., 200., 1., 2000., "P_{T}(1st #mu) (GeV)", true);
  draw_histogram("HLT_Mu3_PFJet40_v_W_John_leadpt", 0., 200., 1., 200., "P_{T}(1st #mu) (GeV)", true);
  */
  //draw_histogram("HLT_Mu17_TrkIsoVVL_v_W_John_Nvtx", 0., 50., 1., 80000., "N(vertice)", false);
  //draw_histogram("HLT_Mu8_TrkIsoVVL_v_W_John_Nvtx", 0., 50., 1., 4000., "N(vertice)", false);
  //draw_histogram("HLT_Mu3_PFJet40_v_W_John_Nvtx", 0., 50., 1., 1000., "N(vertice)", false);

  /*
  draw_histogram("dPhi_jet_lepton_loose", 0., 3.5, 10., 10000., "#Delta #phi (jet, lepton)", false);
  draw_histogram("dPhi_jet_lepton_tight", 0., 3.5, 10., 2000., "#Delta #phi (jet, lepton)", false);
  
  draw_histogram("Pt_jet_over_Pt_lepton_loose", 0., 5., 10., 10000., "P_{T}^{jet} / P_{T}^{lepton}", false);
  draw_histogram("Pt_jet_over_Pt_lepton_tight", 0., 5., 10., 10000., "P_{T}^{jet} / P_{T}^{lepton}", false);
  */
  /*
  draw_histogram("HLT_Mu17_TrkIsoVVL_v_Ptcone_MET_inclusive_loose", 0., 300., 5, 100., "P_{T}^{cone} (#mu)", true);
  draw_histogram("HLT_Mu17_TrkIsoVVL_v_Ptcone_MET_inclusive_tight", 0., 200., 5, 100., "P_{T}^{cone} (#mu)", true);
  draw_histogram("HLT_Mu17_TrkIsoVVL_v_Eta_MET_inclusive_loose", -3., 3., 10, 100., "|#eta|(#mu)", false);
  draw_histogram("HLT_Mu17_TrkIsoVVL_v_Eta_MET_inclusive_tight", -3., 3., 10, 100., "|#eta|(#mu)", false);
  */
  /*  
  draw_histogram("HLT_Mu8_TrkIsoVVL_vMET_MET_MT_inclusive_tight", 0., 200., 5., 100000., "MET (GeV)", true);
  draw_histogram("HLT_Mu3_PFJet40_vMET_MET_MT_inclusive_tight",0., 200., 5., 1000., "MET (GeV)", true);
  draw_histogram("HLT_Mu17_TrkIsoVVL_vMET_MET_MT_inclusive_tight",0., 200., 5., 400000., "MET (GeV)", true);
  
  draw_histogram("HLT_Mu8_TrkIsoVVL_vMET_MET_MT_inclusive_tight_loose", 0., 200., 5., 1000000., "MET (GeV)", true);
  draw_histogram("HLT_Mu3_PFJet40_vMET_MET_MT_inclusive_tight_loose",0., 200., 5., 10000., "MET (GeV)", true);
  draw_histogram("HLT_Mu17_TrkIsoVVL_vMET_MET_MT_inclusive_tight_loose",0., 200., 5., 60000000., "MET (GeV)", true);
  
  draw_histogram("HLT_Mu3_PFJet40_vMT_MET_MT_inclusive_tight", 0., 200., 5., 800., "M_{T}(#mu, METvec) (GeV)", true);
  draw_histogram("HLT_Mu8_TrkIsoVVL_vMT_MET_MT_inclusive_tight", 0., 200., 5., 100000., "M_{T}(#mu, METvec) (GeV)", true);
  draw_histogram("HLT_Mu17_TrkIsoVVL_vMT_MET_MT_inclusive_tight", 0., 200., 5., 400000., "M_{T}(#mu, METvec) (GeV)", true);
  
  draw_histogram("HLT_Mu3_PFJet40_vMT_MET_MT_inclusive_tight_loose", 0., 200., 5., 8000., "M_{T}(#mu, METvec) (GeV)", true);
  draw_histogram("HLT_Mu8_TrkIsoVVL_vMT_MET_MT_inclusive_tight_loose", 0., 200., 5., 1000000., "M_{T}(#mu, METvec) (GeV)", true);
  draw_histogram("HLT_Mu17_TrkIsoVVL_vMT_MET_MT_inclusive_tight_loose", 0., 200., 5., 40000000., "M_{T}(#mu, METvec) (GeV)", true);
  
  draw_histogram("HLT_Mu3_PFJet40_vPtcone_MET_MT_inclusive_tight", 0., 200., 5., 4000., "P_{T}^{cone} (#mu)", true);
  draw_histogram("HLT_Mu8_TrkIsoVVL_vPtcone_MET_MT_inclusive_tight", 0., 200., 5., 100000., "P_{T}^{cone} (#mu)",true);
  draw_histogram("HLT_Mu17_TrkIsoVVL_vPtcone_MET_MT_inclusive_tight", 0., 200., 5., 600000., "P_{T}^{cone} (#mu)",true);
  
  draw_histogram("HLT_Mu3_PFJet40_vPtcone_MET_MT_inclusive_loose", 0., 200., 5., 30000., "P_{T}^{cone} (#mu)",true);
  draw_histogram("HLT_Mu8_TrkIsoVVL_vPtcone_MET_MT_inclusive_loose", 0., 200., 5., 4000000., "P_{T}^{cone} (#mu)",true);
  draw_histogram("HLT_Mu17_TrkIsoVVL_vPtcone_MET_MT_inclusive_loose", 0., 200., 5., 200000000., "P_{T}^{cone} (#mu)",true);


  /*
  draw_histogram("SingleMuonTrigger_Dijet_Awayjet_40_all_MUON_HN_LOOSESTv5_SIP_10_pt_cone_F", 0., 500., 1., 1500.); 
  draw_histogram("SingleMuonTrigger_Dijet_Awayjet_40_all_MUON_HN_LOOSESTv5_SIP_10_pt_cone_F0", 0., 500., 1., 8000.);
  draw_histogram("SingleMuonTrigger_Dijet_Awayjet_40_all_MUON_HN_LOOSESTv5_SIP_10_eta_F", -3., 3., 1., 500.);
  draw_histogram("SingleMuonTrigger_Dijet_Awayjet_40_all_MUON_HN_LOOSESTv5_SIP_10_eta_F0", -3., 3., 1., 5000.);
  */
  
  
  //draw_histogram("h_tri_Z_tag_dimu_mass", 0., 500., 1., 1000.);
  //draw_histogram("h_tri_non_Z_tag_pt", 0., 500., 5., 1500.);
  

}// End of Main Function ////////////////////////////////////////////////////// 
