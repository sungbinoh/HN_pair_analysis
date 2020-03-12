#include "canvas_margin.h"
#include "mylib.h"

void openfile_background(TString cyclename, TString samplename, TString dir, TString histname, int year){
  
  cout << "[openfile_background] " << samplename << ", " << histname << ", " << year << endl;

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/rootfiles/" +TString::Itoa(year,10) + "/Background/";
  TString filename = root_file_path + cyclename + "_" + samplename + ".root";
  TFile *current_file = new TFile ((filename)) ;
  gDirectory->cd(dir);
  
  TH1F * current_hist = (TH1F*)gDirectory -> Get(histname);
  if(current_hist){
    current_hist -> SetDirectory(0);
  }
  TH1::AddDirectory(kFALSE);
  
  //if(samplename.Contains("DYJets") && current_hist) current_hist->Scale(DY_scale);
  mapfunc[histname + cyclename + samplename + TString::Itoa(year,10)] = current_hist;
  
  if(mapfunc[histname + cyclename + samplename + TString::Itoa(year,10)]) cout << "integral of [" << histname + cyclename + samplename + TString::Itoa(year,10) << "] : " << mapfunc[histname + cyclename + samplename + TString::Itoa(year,10)] -> Integral() << endl;
  current_file -> Close();
  delete current_file;
}

void openfile_DATA(TString cyclename, TString samplename, TString dir, TString histname, int year){

  cout << "[openfile_DATA] " << samplename << endl;
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/rootfiles/" +TString::Itoa(year,10) + "/DATA/";
  TString filename = root_file_path + cyclename + "_" + samplename + ".root";
  TFile *current_file = new TFile ((filename)) ;
  gDirectory->cd(dir);
  
  TH1F * current_hist = (TH1F*)gDirectory -> Get(histname);
  if(current_hist){
    current_hist -> SetDirectory(0);
  }
  TH1::AddDirectory(kFALSE);

  mapfunc[histname + cyclename + samplename + TString::Itoa(year,10)] = current_hist;

  current_file -> Close();
  delete current_file;
}
    


void make_histogram(TString nameofhistogram, TString current_histname, int N_bin, double binx[]){

  //cout << "[make_histogram] nameofhistogram : " << nameofhistogram << ", current_histname : " << current_histname << endl;
 
  if(debug){
    cout << "[[make_histogram]] checking binning arrary" << endl;
    for(int i = 0; i < N_bin; i++){
      cout << binx[i] << ", " << endl;
    }
  }
  
  TString current_data;
  if(nameofhistogram.Contains("EMu") || nameofhistogram.Contains("DiMu")) current_data = SingleMuon;
  else if(nameofhistogram.Contains("DiEle")) current_data = DoubleEG;
  else return;


  double N_DY = 0.;
  double N_MC = 0.;
  double N_Data = 0.;
  double DY_norm_2016 = 1., DY_norm_2017 = 1., DY_norm_2018 = 1.;;
  double DY_norm = 1.;
  if( GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2016"].at(0) + "2016") && GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2017"].at(0) + "2017") && GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2018"].at(0) + "2018")){
      
    N_DY += GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2016"].at(0) + "2016") -> Integral();
    N_DY += GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2017"].at(0) + "2017") -> Integral();
    N_DY += GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2018"].at(0) + "2018") -> Integral();

    cout << "N_DY : " << N_DY << endl;
    
    if(GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2016"].at(0) + "2016"))     N_MC += GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2016"].at(0) + "2016") -> Integral();
    if(GetHist(nameofhistogram + Cycle_name + map_sample_names["ttbar2016"].at(0) + "2016") ) N_MC += GetHist(nameofhistogram + Cycle_name + map_sample_names["ttbar2016"].at(0) + "2016") -> Integral();
    if(GetHist(nameofhistogram + Cycle_name + map_sample_names["WJets2016"].at(0) + "2016") )  N_MC += GetHist(nameofhistogram + Cycle_name + map_sample_names["WJets2016"].at(0) + "2016") -> Integral();
    if(GetHist(nameofhistogram + Cycle_name + map_sample_names["Other2016"].at(0) + "2016") )        N_MC += GetHist(nameofhistogram + Cycle_name + map_sample_names["Other2016"].at(0) + "2016") -> Integral();
    
    if(GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2017"].at(0) + "2017"))     N_MC += GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2017"].at(0) + "2017") -> Integral();
    if(GetHist(nameofhistogram + Cycle_name + map_sample_names["ttbar2017"].at(0) + "2017") ) N_MC += GetHist(nameofhistogram + Cycle_name + map_sample_names["ttbar2017"].at(0) + "2017") -> Integral();
    if(GetHist(nameofhistogram + Cycle_name + map_sample_names["WJets2017"].at(0) + "2017") )  N_MC += GetHist(nameofhistogram + Cycle_name + map_sample_names["WJets2017"].at(0) + "2017") -> Integral();
    if(GetHist(nameofhistogram + Cycle_name + map_sample_names["Other2017"].at(0) + "2017") )        N_MC += GetHist(nameofhistogram + Cycle_name + map_sample_names["Other2017"].at(0) + "2017") -> Integral();

    if(GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2018"].at(0) + "2018"))     N_MC += GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2018"].at(0) + "2018") -> Integral();
    if(GetHist(nameofhistogram + Cycle_name + map_sample_names["ttbar2018"].at(0) + "2018") ) N_MC += GetHist(nameofhistogram + Cycle_name + map_sample_names["ttbar2018"].at(0) + "2018") -> Integral();
    if(GetHist(nameofhistogram + Cycle_name + map_sample_names["WJets2018"].at(0) + "2018") )  N_MC += GetHist(nameofhistogram + Cycle_name + map_sample_names["WJets2018"].at(0) + "2018") -> Integral();
    if(GetHist(nameofhistogram + Cycle_name + map_sample_names["Other2018"].at(0) + "2018") )        N_MC += GetHist(nameofhistogram + Cycle_name + map_sample_names["Other2018"].at(0) + "2018") -> Integral();

    if(GetHist(nameofhistogram + Cycle_name + current_data + "2016") )N_Data += GetHist(nameofhistogram + Cycle_name + current_data + "2016") -> Integral();
    if(GetHist(nameofhistogram + Cycle_name + current_data + "2017") )N_Data += GetHist(nameofhistogram + Cycle_name + current_data + "2017") -> Integral();
    if(GetHist(nameofhistogram + Cycle_name + current_data + "2018") )N_Data += GetHist(nameofhistogram + Cycle_name + current_data + "2018") -> Integral();

    DY_norm = (N_DY + N_Data - N_MC) / N_DY;
   
   
    if(nameofhistogram.Contains("DiMu")){
      DY_norm_2016 = 0.801;
      DY_norm_2017 = 1.084;
      DY_norm_2018 = 1.06;
    }
    if(nameofhistogram.Contains("DiEle")){
      DY_norm_2016 = 0.907;
      DY_norm_2017 = 1.078;
      DY_norm_2018 = 1.031;
    }

    GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2016"].at(0) + "2016") -> Scale(DY_norm_2016); 
    GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2017"].at(0) + "2017") -> Scale(DY_norm_2017);
    GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2018"].at(0) + "2018") -> Scale(DY_norm_2018);

    
    outfile << nameofhistogram << "'s DY_norm : " << DY_norm << endl;
    
    for(int i_syst = 1; i_syst < N_syst; i_syst++){
      cout << current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["DYJets2016"].at(0) + "2016" << endl;
      if(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["DYJets2016"].at(0) + "2016")){
	GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["DYJets2016"].at(0) + "2016") -> Scale(DY_norm_2016);
	cout << "DY interal 2016 syst " <<  systematics[i_syst] << " : " << GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["DYJets2016"].at(0) + "2016") -> Integral() << endl;
      }
      if(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["DYJets2017"].at(0) + "2017")){
        GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["DYJets2017"].at(0) + "2017") -> Scale(DY_norm_2017);
      }
      if(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["DYJets2018"].at(0) + "2018")){
        GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["DYJets2018"].at(0) + "2018") -> Scale(DY_norm_2018);
      }
    }
  }
  
  // -- Add three years histograms
  Int_t nx_ttbar    = GetHist(nameofhistogram + Cycle_name + map_sample_names["ttbar2016"].at(0) + "2016") -> GetNbinsX();
  Double_t x1_ttbar = GetHist(nameofhistogram + Cycle_name + map_sample_names["ttbar2016"].at(0) + "2016") -> GetBinLowEdge(1);
  Double_t x2_ttbar = GetHist(nameofhistogram + Cycle_name + map_sample_names["ttbar2016"].at(0) + "2016") -> GetBinLowEdge(nx_ttbar+1);
  TH1F *htmp_add = new TH1F("", "", nx_ttbar, x1_ttbar, x2_ttbar);
  
  mapfunc[nameofhistogram + Cycle_name + map_sample_names["DYJets"].at(0)] = (TH1F*)htmp_add -> Clone("DYJets");
  if(GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2016"].at(0) + "2016")) mapfunc[nameofhistogram + Cycle_name + map_sample_names["DYJets"].at(0)] -> Add(GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2016"].at(0) + "2016"));
  if(GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2017"].at(0) + "2017")) mapfunc[nameofhistogram + Cycle_name + map_sample_names["DYJets"].at(0)] -> Add(GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2017"].at(0) + "2017"));
  if(GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2018"].at(0) + "2018")) mapfunc[nameofhistogram + Cycle_name + map_sample_names["DYJets"].at(0)] -> Add(GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets2018"].at(0) + "2018"));
  cout << "[make_histogram] intgeral DYJet : " << mapfunc[nameofhistogram + Cycle_name + map_sample_names["DYJets"].at(0)] -> Integral() << endl;
  
  mapfunc[nameofhistogram + Cycle_name + map_sample_names["ttbar"].at(0)] = (TH1F*)htmp_add -> Clone("ttbar");
  if(GetHist(nameofhistogram + Cycle_name + map_sample_names["ttbar2016"].at(0) + "2016")) mapfunc[nameofhistogram + Cycle_name + map_sample_names["ttbar"].at(0)] -> Add(GetHist(nameofhistogram + Cycle_name + map_sample_names["ttbar2016"].at(0) + "2016"));
  if(GetHist(nameofhistogram + Cycle_name + map_sample_names["ttbar2017"].at(0) + "2017")) mapfunc[nameofhistogram + Cycle_name + map_sample_names["ttbar"].at(0)] -> Add(GetHist(nameofhistogram + Cycle_name + map_sample_names["ttbar2017"].at(0) + "2017"));
  if(GetHist(nameofhistogram + Cycle_name + map_sample_names["ttbar2018"].at(0) + "2018")) mapfunc[nameofhistogram + Cycle_name + map_sample_names["ttbar"].at(0)] -> Add(GetHist(nameofhistogram + Cycle_name + map_sample_names["ttbar2018"].at(0) + "2018"));
  
  mapfunc[nameofhistogram + Cycle_name + map_sample_names["WJets"].at(0)] = (TH1F*)htmp_add -> Clone("WJets");
  if(GetHist(nameofhistogram + Cycle_name + map_sample_names["WJets2016"].at(0) + "2016")) mapfunc[nameofhistogram + Cycle_name + map_sample_names["WJets"].at(0)] -> Add( GetHist(nameofhistogram + Cycle_name + map_sample_names["WJets2016"].at(0) + "2016"));
  if(GetHist(nameofhistogram + Cycle_name + map_sample_names["WJets2017"].at(0) + "2017")) mapfunc[nameofhistogram + Cycle_name + map_sample_names["WJets"].at(0)] -> Add(GetHist(nameofhistogram + Cycle_name + map_sample_names["WJets2017"].at(0) + "2017"));
  if(GetHist(nameofhistogram + Cycle_name + map_sample_names["WJets2018"].at(0) + "2018")) mapfunc[nameofhistogram + Cycle_name + map_sample_names["WJets"].at(0)] -> Add(GetHist(nameofhistogram + Cycle_name + map_sample_names["WJets2018"].at(0) + "2018"));
  
  mapfunc[nameofhistogram + Cycle_name + map_sample_names["Other"].at(0)] = (TH1F*)htmp_add -> Clone("Other");
  if(GetHist(nameofhistogram + Cycle_name + map_sample_names["Other2016"].at(0) + "2016")) mapfunc[nameofhistogram + Cycle_name + map_sample_names["Other"].at(0)] -> Add( GetHist(nameofhistogram + Cycle_name + map_sample_names["Other2016"].at(0) + "2016"));
  if(GetHist(nameofhistogram + Cycle_name + map_sample_names["Other2017"].at(0) + "2017")) mapfunc[nameofhistogram + Cycle_name + map_sample_names["Other"].at(0)] -> Add(GetHist(nameofhistogram + Cycle_name + map_sample_names["Other2017"].at(0) + "2017"));
  if(GetHist(nameofhistogram + Cycle_name + map_sample_names["Other2018"].at(0) + "2018")) mapfunc[nameofhistogram + Cycle_name + map_sample_names["Other"].at(0)] -> Add(GetHist(nameofhistogram + Cycle_name + map_sample_names["Other2018"].at(0) + "2018"));
  
  mapfunc[nameofhistogram + Cycle_name + current_data] = GetHist(nameofhistogram + Cycle_name + current_data + "2016");
  mapfunc[nameofhistogram + Cycle_name + current_data] -> Add(GetHist(nameofhistogram + Cycle_name + current_data + "2017"));
  mapfunc[nameofhistogram + Cycle_name + current_data] -> Add(GetHist(nameofhistogram + Cycle_name + current_data + "2018"));
  

  /*
  for(int i_syst = 1; i_syst < N_syst; i_syst++){
    mapfunc[current_histname + systematics[i_syst] + Cycle_name + map_sample_names["DYJets"].at(0)] = (TH1F*)htmp_add -> Clone("DYJets" + systematics[i_syst]);
    if(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["DYJets2016"].at(0) + "2016")) mapfunc[current_histname + systematics[i_syst] + Cycle_name + map_sample_names["DYJets"].at(0)] -> Add(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["DYJets2016"].at(0) + "2016"));
    if(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["DYJets2017"].at(0) + "2017")) mapfunc[current_histname + systematics[i_syst] + Cycle_name + map_sample_names["DYJets"].at(0)] -> Add(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["DYJets2017"].at(0) + "2017"));
    if(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["DYJets2018"].at(0) + "2018")) mapfunc[current_histname + systematics[i_syst] + Cycle_name + map_sample_names["DYJets"].at(0)] -> Add(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["DYJets2018"].at(0) + "2018"));
    
    cout << "DY interal all syst " << systematics[i_syst] << " : " << mapfunc[current_histname + systematics[i_syst] + Cycle_name + map_sample_names["DYJets"].at(0)] -> Integral() << endl;

    
    mapfunc[current_histname + systematics[i_syst] + Cycle_name + map_sample_names["ttbar"].at(0)] = (TH1F*)htmp_add -> Clone("ttbar" + systematics[i_syst]);
    if(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["ttbar2016"].at(0) + "2016")) mapfunc[current_histname + systematics[i_syst] + Cycle_name + map_sample_names["ttbar"].at(0)] -> Add(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["ttbar2016"].at(0) + "2016"));
    if(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["ttbar2017"].at(0) + "2017")) mapfunc[current_histname + systematics[i_syst] + Cycle_name + map_sample_names["ttbar"].at(0)] -> Add(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["ttbar2017"].at(0) + "2017"));
    if(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["ttbar2018"].at(0) + "2018")) mapfunc[current_histname + systematics[i_syst] + Cycle_name + map_sample_names["ttbar"].at(0)] -> Add(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["ttbar2018"].at(0) + "2018"));
    
    mapfunc[current_histname + systematics[i_syst] + Cycle_name + map_sample_names["WJets"].at(0)] = (TH1F*)htmp_add -> Clone("WJets" + systematics[i_syst]);
    if(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["WJets2016"].at(0) + "2016")) mapfunc[current_histname + systematics[i_syst] + Cycle_name + map_sample_names["WJets"].at(0)] -> Add(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["WJets2016"].at(0) + "2016"));
    if(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["WJets2017"].at(0) + "2017")) mapfunc[current_histname + systematics[i_syst] + Cycle_name + map_sample_names["WJets"].at(0)] -> Add(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["WJets2017"].at(0) + "2017"));
    if(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["WJets2018"].at(0) + "2018")) mapfunc[current_histname + systematics[i_syst] + Cycle_name + map_sample_names["WJets"].at(0)] -> Add(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["WJets2018"].at(0) + "2018"));

    mapfunc[current_histname + systematics[i_syst] + Cycle_name + map_sample_names["Other"].at(0)] = (TH1F*)htmp_add -> Clone("Other" + systematics[i_syst]);
    if(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["Other2016"].at(0) + "2016")) mapfunc[current_histname + systematics[i_syst] + Cycle_name + map_sample_names["Other"].at(0)] -> Add(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["Other2016"].at(0) + "2016"));
    if(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["Other2017"].at(0) + "2017")) mapfunc[current_histname + systematics[i_syst] + Cycle_name + map_sample_names["Other"].at(0)] -> Add(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["Other2017"].at(0) + "2017"));
    if(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["Other2018"].at(0) + "2018")) mapfunc[current_histname + systematics[i_syst] + Cycle_name + map_sample_names["Other"].at(0)] -> Add(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["Other2018"].at(0) + "2018"));
  }
  */
  
  
  if(debug){
    cout << "[[make_histogram]] Added three years" << endl;
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
  canvas_margin(mapcanvas[canvas]);
  gStyle -> SetOptStat(1111);

  ////////////////////////////////////
  /// Pad 1
  ////////////////////////////////////
  mappad[pad1] = new TPad("", "", 0, 0.25, 1, 1);
  mappad[pad1] -> SetTopMargin( 0.07 );
  mappad[pad1] -> SetBottomMargin( 0.02 );
  mappad[pad1] -> SetLeftMargin( 0.15 );
  mappad[pad1] -> SetRightMargin( 0.03 );
  mappad[pad1] -> Draw();
  mappad[pad1] -> cd();
  //mappad[pad1] -> SetLogy();
  
  maplegend[legend] = new TLegend(0.60, 0.60, 0.96, 0.92);

  if(debug) cout << "1" << endl;
  
  maphstack[hstack] = new THStack(hstack, "Stacked_" + nameofhistogram);
  gStyle->SetOptTitle(0);
  
  TString legend_list[4] = {"Other", "WJets", "ttbar", "DYJets"};
  Int_t colour_array[] = {419, 901, 416, 400};

  TString name_cycle = nameofhistogram + Cycle_name;
  if(debug) cout << "check1" << endl;
  
  // -- Add Overflow bin
  TString overflow = "overflow";
  //cout << nameofhistogram + Cycle_name + current_data << endl;
  //cout << "nameofhistogram + Cycle_name + current_data : " << nameofhistogram + Cycle_name + current_data << endl;
  Int_t nx_template = GetHist(nameofhistogram + Cycle_name + current_data) -> GetNbinsX()+1;
  Double_t x1_template = GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinLowEdge(1);
  Double_t bw_template = binx[N_bin - 1] - binx[N_bin- 2];
  Double_t x2_template = GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinLowEdge(nx_template)+bw_template;
  binx[N_bin] = binx[N_bin - 1] + bw_template;
  if(debug) cout << "[[make_histogram]] binx[N_bin] ; " << binx[N_bin]  << endl;
  mapfunc[func] = new TH1F("", "", nx_template, x1_template, x2_template);
  if(debug) cout << "func rebin rebinning" << endl;
  mapfunc[func + "rebin"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "rebin", binx);
  mapfunc[func + "rebin_stat_err"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "rebin", binx);
  mapfunc[func + "blind_data"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "blind_Data", binx);
  // -- For each Bkg Legend : Add overflow, Rebin, Syst error bar, Add on legend, and set fill & line colors
  for(int i_legend = 0; i_legend < 4; i_legend++){
    if(debug) cout << legend_list[i_legend] << endl;
    TString current_sample = map_sample_names[legend_list[i_legend]].at(0);
    if(mapfunc[nameofhistogram + Cycle_name + current_sample]){
      Int_t nx    = GetHist(nameofhistogram + Cycle_name + current_sample) -> GetNbinsX()+1;
      Double_t x1 = GetHist(nameofhistogram + Cycle_name + current_sample) -> GetBinLowEdge(1);
      Double_t bw = bw_template;
      Double_t x2 = GetHist(nameofhistogram + Cycle_name + current_sample) -> GetBinLowEdge(nx)+bw_template;
      
      TH1F *htmp = new TH1F("", "", nx, x1, x2);
      TH1F *htmp_stat_err = new TH1F("", "", nx, x1, x2);
      
      // -- Make rebinned central hists
      for(Int_t j = 1; j <= nx; j++){
        double current_bin_content = GetHist(nameofhistogram + Cycle_name + current_sample) -> GetBinContent(j);
	double current_bin_error = GetHist(nameofhistogram + Cycle_name + current_sample) -> GetBinError(j);
        htmp -> SetBinContent(j, current_bin_content);
        htmp -> SetBinError(j, current_bin_error);
        htmp_stat_err -> SetBinContent(j, current_bin_content);
        htmp_stat_err -> SetBinError(j, current_bin_error);
      }
      
      if(debug) cout << current_sample + " : called all syst hists" << endl;
      
      mapfunc[nameofhistogram + Cycle_name + current_sample + "rebin"] = dynamic_cast<TH1F*>(htmp -> Rebin(N_bin, nameofhistogram + Cycle_name + current_sample + "rebin", binx));
      mapfunc[nameofhistogram + Cycle_name + current_sample + "rebin_stat_err"] = dynamic_cast<TH1F*>(htmp_stat_err -> Rebin(N_bin, nameofhistogram + Cycle_name + current_sample + "rebin", binx));
      
      
      // -- Make rebinned hists for all systematic categories
      for(int i_syst = 1; i_syst < N_syst; i_syst++){
        //if(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + current_sample)){
	if(GetHist(current_histname + systematics[i_syst] + Cycle_name + current_sample)){ // -- FIXME

	  if(debug) cout << "Filling " << systematics[i_syst] << endl;
          TH1F *htmp_syst = new TH1F("", "", nx, x1, x2);
          for(Int_t j = 1; j <= nx; j++){
            //double current_bin_content = GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + current_sample) -> GetBinContent(j);
            double current_bin_content = GetHist(current_histname + systematics[i_syst] + Cycle_name + current_sample) -> GetBinContent(j); // -- FIXME
	    htmp_syst -> SetBinContent(j, current_bin_content);
          }
          mapfunc[nameofhistogram + systematics[i_syst] + current_sample + "rebin"] = dynamic_cast<TH1F*>(htmp_syst -> Rebin(N_bin, nameofhistogram + Cycle_name + current_sample + "rebin", binx));
        }
      }

      vector<double> error_vector = Get_Syst_Error(nameofhistogram, current_sample);

      for(unsigned int j_syst = 0; j_syst < error_vector.size(); j_syst++){
	mapfunc[nameofhistogram + Cycle_name + current_sample + "rebin"] -> SetBinError(j_syst, error_vector.at(j_syst));
      }
      mapfunc[nameofhistogram + Cycle_name + current_sample + "rebin"] ->  SetFillColor(colour_array[i_legend]);
      mapfunc[nameofhistogram + Cycle_name + current_sample + "rebin"] ->  SetLineColor(colour_array[i_legend]);
      
      // -- Prepare blind data point (= bkg)
      mapfunc[func + "blind_data"]  -> Add(GetHist(nameofhistogram + Cycle_name + current_sample + "rebin"));
      
      maphstack[hstack] -> Add(GetHist(nameofhistogram + Cycle_name + current_sample + "rebin"));
      //maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + current_sample + "rebin"), legend_list[i_legend], "lf");
      mapfunc[func + "rebin"] -> Add(GetHist(nameofhistogram + Cycle_name + current_sample + "rebin"));
      mapfunc[func + "rebin_stat_err"] -> Add(GetHist(nameofhistogram + Cycle_name + current_sample + "rebin_stat_err"));
    }
  }
  
  for(int i_legend = 0; i_legend < 4; i_legend++){
    TString current_sample = map_sample_names[legend_list[4 - i_legend - 1]].at(0);
    if(GetHist(nameofhistogram + Cycle_name + current_sample + "rebin")) maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + current_sample + "rebin"), legend_list[4 - i_legend - 1], "lf");
  }
  
  if(debug) cout << "4" << endl;

  mapfunc[nameofhistogram + Cycle_name + current_data + overflow] = new TH1F("", "", nx_template, x1_template, x2_template);
  for (Int_t i = 1; i <= nx_template; i++) {
    mapfunc[nameofhistogram + Cycle_name + current_data +overflow] -> SetBinContent(i, GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinContent(i)) ;
    mapfunc[nameofhistogram + Cycle_name + current_data +overflow] -> SetBinError(i, GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinError(i));
  }
  if(blind){
    mapfunc[nameofhistogram + Cycle_name + current_data + "rebin"] = (TH1F*)GetHist(func + "blind_data") -> Clone(clone + "blind");
  }
  else{
    mapfunc[nameofhistogram + Cycle_name + current_data + "rebin"] = dynamic_cast<TH1F*>(GetHist(nameofhistogram + Cycle_name + current_data + overflow) -> Rebin(N_bin, nameofhistogram + Cycle_name + current_data + "rebin", binx));
  }
  
  // -- Add DATA on Legend
  if(!nameofhistogram.Contains("SR")){
    maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + current_data + "rebin"), "data", "lp");
  }
  GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> SetMarkerStyle(20);
  GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> SetMarkerColor(kBlack);
  GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> SetMarkerSize(1);

  if(debug) cout << "5" << endl;
  
  maphstack[hstack] -> Draw("hist");
  maphstack[hstack] -> GetYaxis()->SetLabelSize(0.05);;
  maphstack[hstack] -> GetYaxis()->SetTitleSize(0.07);;
  maphstack[hstack] -> GetYaxis()->SetTitleOffset(1.02);;
  maphstack[hstack] -> GetXaxis()->SetLabelSize(0);
  maphstack[hstack] -> GetXaxis() -> SetRangeUser(x1_template, x2_template);
  
  maphstack[hstack] -> GetXaxis() -> SetTitle(nameofhistogram);

  maphstack[hstack] -> GetYaxis() -> SetTitle(title_y);
  maphstack[hstack] -> Draw("histsame");
  
  
  // -- Proper error for data : 
  Proper_error_data(nameofhistogram, current_data, N_bin, binx);
  map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetMarkerColor(kBlack);
  map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> Draw("p0same");

  // -- Draw data
  GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetXaxis() -> SetRangeUser(x1_template, x2_template);
  mappad[pad1] -> Update();
  //GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> SetMarkerColor(kRed);
  GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> Draw("histsamep");
  
    
  
  // -- Draw Bkg Error bar
  mapfunc[func + "rebin"] -> SetFillColor(kBlack);;
  mapfunc[func + "rebin"] -> SetFillStyle(3013);
  mapfunc[func + "rebin"] -> SetMarkerSize(0);
  mapfunc[func + "rebin"] -> Draw("e2same");
  maplegend[legend] -> AddEntry(mapfunc[func + "rebin"], "Syst. + Stat.", "pf");
  
  
  // -- Set y-axis range
  double data_max = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetMaximum();
  double data_min = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetMinimum();
  //maphstack[hstack] -> SetMaximum(data_max * 100000.);//logy
  maphstack[hstack] -> SetMaximum(data_max * 1.5);
  maphstack[hstack] -> SetMinimum(0.0001);
  
  // -- Draw Legend
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kBlack);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(1001);
  maplegend[legend] -> SetShadowColor(0); // 0 = transparent 
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> Draw("same");
  
  if(debug) cout << "6" << endl;


  ////////////////////////////////////
  /// Pad 2
  ////////////////////////////////////
  mapcanvas[canvas] -> cd();
  mappad[pad2] = new TPad(pad2, "", 0, 0, 1, 0.25);
  mappad[pad2] -> SetTopMargin( 0.03 );
  mappad[pad2] -> SetBottomMargin( 0.4 );
  mappad[pad2] -> SetLeftMargin( 0.15 );
  mappad[pad2] -> SetRightMargin( 0.03 );
  mappad[pad2] -> Draw();
  mappad[pad2] -> cd();
  
  // -- Make Clone hist of Data (For ratio)
  mapfunc[clone] = (TH1F*)GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> Clone(clone);
  
  // -- Make Clone hist for Bkg. One for total error, and one for stat error
  mapfunc["stat" + nameofhistogram] = (TH1F*)GetHist(func + "rebin") -> Clone(clone);
  mapfunc["stat_err" + nameofhistogram] = (TH1F*)GetHist(func + "rebin_stat_err") -> Clone(clone);
  
  if(debug) cout << "7" << endl;

  // -- Draw Line at Data/MC = 1
  mapline[line] = new TLine(binx[0], 1, binx[N_bin - 1], 1);
  mapline[line] -> Draw();
  mapline[line] -> SetLineStyle(1);
  mapline[line] -> SetLineColor(kBlue);
  
  // -- Draw syst box hist with all contents as 1
  TString name_x = nameofhistogram;
  mapfunc["stat" + nameofhistogram] -> SetTitle("");
  mapfunc["stat" + nameofhistogram] -> SetFillColor(kOrange);
  mapfunc["stat" + nameofhistogram] -> SetMarkerSize(0);
  mapfunc["stat" + nameofhistogram] -> SetMarkerStyle(0);
  mapfunc["stat" + nameofhistogram] -> SetLineColor(kWhite);
  mapfunc["stat" + nameofhistogram] -> GetXaxis() -> SetTitle(name_x);
  mapfunc["stat" + nameofhistogram] -> GetXaxis() -> SetTitleSize(0.15);
  mapfunc["stat" + nameofhistogram] -> GetXaxis() -> SetLabelSize(0.10);
  mapfunc["stat" + nameofhistogram] -> GetXaxis() -> SetRangeUser(x1_template, x2_template);
  
  mapfunc["stat" + nameofhistogram] -> SetYTitle("#frac{Obs.}{Pred.}");
  mapfunc["stat" + nameofhistogram] -> GetYaxis() -> SetTitleSize(0.12);
  mapfunc["stat" + nameofhistogram] -> GetYaxis() -> SetTitleOffset(0.5);
  mapfunc["stat" + nameofhistogram] -> GetYaxis() -> SetLabelSize(0.08);
  mapfunc["stat" + nameofhistogram] -> GetYaxis() -> SetNdivisions(505);
  
  mapfunc["stat" + nameofhistogram] -> SetFillStyle(1001);
  mapfunc["stat" + nameofhistogram] -> SetFillColorAlpha(45,0.35);
  mapfunc["stat" + nameofhistogram] -> SetMinimum(0.5);
  mapfunc["stat" + nameofhistogram] -> SetMaximum(1.5);
  mapfunc["stat" + nameofhistogram] -> SetStats(0);
  
  if(debug) cout << "8" << endl;

  Int_t ncells = mapfunc["stat" + nameofhistogram] -> GetSize();
  for(int i = 0; i < N_bin; i++){
    double current_error = mapfunc["stat" + nameofhistogram] -> GetBinError(i + 1) / mapfunc["stat" + nameofhistogram] -> GetBinContent(i + 1);
    double original_error = mapfunc["stat" + nameofhistogram] -> GetBinError(i + 1);
    double original_content = mapfunc["stat" + nameofhistogram] -> GetBinContent(i + 1);

    bool no_data = false;
    double data_value = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetBinContent(i + 1);;
    if(data_value < 10e-7) current_error = 0.;

    mapfunc["stat" + nameofhistogram] -> SetBinContent(i + 1, 1.);
    mapfunc["stat" + nameofhistogram] -> SetBinError(i + 1, current_error);

  }
  mapfunc["stat" + nameofhistogram] -> Draw("CE2");
  

  // -- Draw stat error box hist with all contents as 1
  mapfunc["stat_err" + nameofhistogram] -> SetMarkerSize(0);
  mapfunc["stat_err" + nameofhistogram] -> SetMarkerStyle(0);
  mapfunc["stat_err" + nameofhistogram] -> SetLineColor(kRed);
  mapfunc["stat_err" + nameofhistogram] -> SetFillColor(kRed);
  mapfunc["stat_err" + nameofhistogram] -> Draw("E2same");

  for(int i = 0; i < N_bin; i++){
    double current_error = mapfunc["stat_err" + nameofhistogram] -> GetBinError(i + 1) / mapfunc["stat_err" + nameofhistogram] -> GetBinContent(i + 1);
    double original_error = mapfunc["stat_err" + nameofhistogram] -> GetBinError(i + 1);
    double original_content = mapfunc["stat_err" + nameofhistogram] -> GetBinContent(i + 1);

    bool no_data = false;
    double data_value = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetBinContent(i + 1);;
    if(data_value < 10e-7) current_error = 0.;

    mapfunc["stat_err" + nameofhistogram] -> SetBinContent(i + 1, 1.);
    mapfunc["stat_err" + nameofhistogram] -> SetBinError(i + 1, current_error);

  }

  // -- Make Ratio hist
  mapfunc[clone] -> Divide(mapfunc[func + "rebin"]);
  for(int i = 0; i < N_bin; i++){
    double bkg_value = mapfunc[func + "rebin"] -> GetBinContent(i + 1);
    double data_value = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetBinContent(i + 1);
    if(debug) cout << " data_value : " << data_value << endl;
    if(debug) cout << " GetHist(nameofhistogram + Cycle_name + current_data + rebin) -> GetBinContent(" << i + 1 << ") : " << GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetBinContent(i + 1) << endl;
    if(bkg_value < 10e-7){
      mapfunc[clone] -> SetBinContent(i + 1, 999.);
    }
    else if(data_value < 10e-7) mapfunc[clone] -> SetBinContent(i + 1, 0.);
  }
  

  // -- Correct Data Stat Error bar
  map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] = new TGraphAsymmErrors(mapfunc[clone]);
  for(int i = 0; i < N_bin; i++){
    double bkg_i = mapfunc[func + "rebin"] -> GetBinContent(i+1);
    //cout << "veyl.at(i) : " << veyl.at(i) << ", bkg_i : " << bkg_i << ", veyl.at(i) / bkg_i : " << veyl.at(i) / bkg_i << endl;
    
    map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] -> SetPointEYlow(i, veyl.at(i) / bkg_i);
    map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] -> SetPointEXlow(i, 0.);
    map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] -> SetPointEYhigh(i, veyh.at(i) / bkg_i);
    map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] -> SetPointEXhigh(i, 0.);
    if(bkg_i < 10e-7){
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] -> SetPointEYlow(i, 0.);
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] -> SetPointEXlow(i, 0.);
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] -> SetPointEYhigh(i, 0.);
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] -> SetPointEXhigh(i, 0.);
    }

  }//end for i on N_bin
  
  map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error_ratio"] -> Draw("p0same");
  mapline[line] -> Draw("same");

  maplegend["bottom" + legend] = new TLegend(0.2, 0.85, 0.4, 0.95);
  maplegend["bottom" + legend]->SetBorderSize(0);
  maplegend["bottom" + legend]->SetNColumns(3);
  maplegend["bottom" + legend]->AddEntry(mapfunc["stat_err" + nameofhistogram], "Stat.", "f");
  maplegend["bottom" + legend]->AddEntry(mapfunc["stat" + nameofhistogram], "Stat. + Syst.", "f");
  maplegend["bottom" + legend]->AddEntry(mapfunc[clone] , "Obs./Pred.", "p");
  maplegend["bottom" + legend]->Draw("same");

  
  ////////////////////////////////////
  /// -- Latex
  ////////////////////////////////////
  mapcanvas[canvas] -> cd();
  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_Lumi.SetTextSize(0.035);
  TString total_lumi = "137.4";
  latex_Lumi.DrawLatex(0.7, 0.96, total_lumi + " fb^{-1} (13 TeV)");

  if(debug) cout << "9" << endl;


  TString pdfname;
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  pdfname = WORKING_DIR + "/plots/" + TString::Itoa(tag_year,10) + "/";
  if(nameofhistogram.Contains("EMu")){
    pdfname = pdfname + "EMu/";

  }
  else if(nameofhistogram.Contains("DiEle")){
    pdfname = pdfname + "DiEle/";

  }
  else if(nameofhistogram.Contains("DiMu")){
    pdfname = pdfname + "DiMu/";
  }
  else return;
  
  pdfname.Append(nameofhistogram);
  pdfname.Append(".pdf");
  if(debug) cout << "9.1" << endl;

  mapcanvas[canvas] -> SaveAs(pdfname);

  if(debug) cout << "10" << endl;
  
  vx.clear();
  vy.clear(); 
  vexl.clear();
  vexh.clear();
  veyl.clear();
  veyh.clear();
  
}


void open_files(TString histname){
  
  if(debug) cout << "[[draw_histogram]] histname : " << histname << endl;  
  
  unsigned int N_bin = map_bin_vector[histname].size();
  double current_bins[100];
  for(unsigned int i_bin = 0; i_bin < N_bin; i_bin++){
    current_bins[i_bin] = map_bin_vector[histname].at(i_bin);
  }


  // -- Get which region(CR1, CR2 , ...) and channel (DiEle, DiMu, ...)
  //int N_regions = 6;
  int N_regions = 4;
  
  TString regions[] = {"CR_Zmass",
                       //"CR_ttbar",
                       //"CR_inv_mll",
                       "tight_CR_Zmass",
                       //"tight_CR_ttbar",
                       "SR",
		       "SR_b_veto",
  };


  int N_channels = 6;
  TString channels[] = {"DiEle",
                        "DiMu",
                        "EMu",
			"DiEle_DYreweight",
                        "DiMu_DYreweight",
                        "EMu_DYreweight",
  };
  

  TString current_dir = "empty";
  TString current_channel = "empty";
  
  // -- Fill current_dir and current_channel strings
  for(int i_region = 0; i_region < N_regions; i_region++){
    if(histname.Contains(regions[i_region])) current_dir = regions[i_region];
  }
  for(int i_channel = 0; i_channel < N_channels; i_channel++){
    // -- FIXME
    if(histname.Contains(channels[i_channel])){
      current_channel = channels[i_channel];
    }
  }

  current_dir = current_dir + "_" + current_channel;

  if(debug) cout << "current_dir : " << current_dir << "/" << histname << "+systflags" << endl;
  if(current_dir.Contains("empty") || current_channel.Contains("empty")) return;
  
  
  for(int i_syst = 0; i_syst < N_syst; i_syst++){
    // -- FIXME
    //TString current_dir_syst = current_dir + "_" + systematics[i_syst];
    TString current_hist_syst = histname + "_" + systematics[i_syst];
    TString current_dir_syst = current_dir + "_" + systematics[i_syst];
    //TString current_hist_syst = histname + "_DYreweight_" + systematics[i_syst]; // -- FIXME
    TString legend_list[4] = {"DYJets", "ttbar", "WJets", "Other"};
    for(int i_legend_list = 0; i_legend_list < 4; i_legend_list++){
      
      unsigned int vector_size = map_sample_names[legend_list[i_legend_list]].size();
      //cout << "legend_list[i_legend_list] : " << legend_list[i_legend_list] << endl;
      //cout << "vector_size : " << vector_size << endl;
      
      for(int i_vector = 0; i_vector < vector_size; i_vector++){
	openfile_background(Cycle_name, map_sample_names[legend_list[i_legend_list] + "2016"].at(i_vector), current_dir_syst, current_hist_syst, 2016);
	openfile_background(Cycle_name, map_sample_names[legend_list[i_legend_list] + "2017"].at(i_vector), current_dir_syst, current_hist_syst, 2017);
        openfile_background(Cycle_name, map_sample_names[legend_list[i_legend_list] + "2018"].at(i_vector), current_dir_syst, current_hist_syst, 2018);

      }
    }
    
    //cout << "current_dir_syst : " << current_dir_syst << endl;
    //cout << "current_hist_syst : " << current_hist_syst << endl;
    openfile_DATA(Cycle_name, DoubleEG, current_dir_syst, current_hist_syst, 2016);
    openfile_DATA(Cycle_name, DoubleEG, current_dir_syst, current_hist_syst, 2017);
    openfile_DATA(Cycle_name, DoubleEG, current_dir_syst, current_hist_syst, 2018);

    openfile_DATA(Cycle_name, SingleMuon, current_dir_syst, current_hist_syst, 2016);
    openfile_DATA(Cycle_name, SingleMuon, current_dir_syst, current_hist_syst, 2017);
    openfile_DATA(Cycle_name, SingleMuon, current_dir_syst, current_hist_syst, 2018);

  }
  // -- FIXME
  //histname = histname + "_central";
  //cout << "histname : " << histname << ", current_histname : " << current_histname << endl;
  make_histogram(histname + "_central", histname, N_bin, current_bins);
  
}

void open_binning_file(TString filename){
  
  cout << "[open_binning_file] start to open binngin file : " << filename << endl;
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString binning_file_path=WORKING_DIR+"/script/"+filename;
  ifstream data_file;
  data_file.open(binning_file_path);
 
  char line[500];
  if(data_file.is_open()){
    while(!data_file.eof()){
      data_file.getline(line, 500);
      if(debug) cout << line << endl;
      TString this_line = line;
      if(this_line.Contains("#")) continue;
      
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
      map_bin_vector[current_histname + "_DYreweight"].clear();

      
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
	map_bin_vector[current_histname + "_DYreweight"].push_back(current_bin_double);
      }

      cout << "[[open_binning_file]] current_histname : " << current_histname << endl;
      
      open_files(current_histname);
      open_files(current_histname + "_DYreweight");
      map_bin_vector.clear();
      
    }//while not end of file
  }//end if file open
  data_file.close();
  
}

void QuickPlotFullRun2(int year=2019){
  setTDRStyle();
  
  // == Set Input Sample List
  tag_year = year;
  
  //outfile.open("output_" + TString::Itoa(tag_year,10) + ".txt");
  
  //map_sample_names["DYJets"] = {"DYJets_MG"};
  map_sample_names["DYJets2018"] = {"DYJets_MG_HT"};
  map_sample_names["ttbar2018"] = {"TT_powheg"};
  map_sample_names["WJets2018"] = {"WJets_MG"};
  map_sample_names["Other2018"] = {"VV"};
  map_sample_names["Muon2018"] = {"data_SingleMuon"};
  map_sample_names["EGamma2018"] = {"data_DoubleEG"};
    
  SingleMuon = "data_SingleMuon";
  DoubleEG = "data_DoubleEG";
  

  //map_sample_names["DYJets"] = {"DYJets_MG"};
  map_sample_names["DYJets2017"] = {"DYJets_MG_HT"};
  map_sample_names["ttbar2017"] = {"TT_powheg"};
  map_sample_names["WJets2017"] = {"WJets_MG"};
  map_sample_names["Other2017"] = {"VV"};
  map_sample_names["Muon2017"] = {"data_SingleMuon"};
  map_sample_names["EGamma2017"] = {"data_DoubleEG"};
  
  //SingleMuon = "data_SingleMuon";
  //DoubleEG = "data_DoubleEG";


  //map_sample_names["DYJets"] = {"DYJets_MG"};
  map_sample_names["DYJets2016"] = {"DYJets_MG_HT"};
  map_sample_names["ttbar2016"] = {"TT_powheg"};
  map_sample_names["WJets2016"] = {"WJets_MG"};
  map_sample_names["Other2016"] = {"VV"};
  map_sample_names["Muon2016"] = {"data_SingleMuon"};
  map_sample_names["EGamma2016"] = {"data_DoubleEG"};


  map_sample_names["DYJets"] = {"DYJets_MG_HT"};
  map_sample_names["ttbar"] = {"TT_powheg"};
  map_sample_names["WJets"] = {"WJets_MG"};
  map_sample_names["Other"] = {"VV"};
  map_sample_names["Muon"] = {"data_SingleMuon"};
  map_sample_names["EGamma"] = {"data_DoubleEG"};
  
  
  //SingleMuon = "data_SingleMuon";
  //DoubleEG = "data_DoubleEG";

  Cycle_name = "HN_pair_all_SkimTree_LRSMHighPt";
  
  //open_binning_file("binning_uniform_test.txt");
  open_binning_file("binning_Zp.txt");
  //open_binning_file("binning_test.txt");

  outfile.close();
  
}









