#include "canvas_margin.h"
#include "mylib.h"

void openfile_background(TString cyclename, TString samplename, TString dir, TString histname){
  
  cout << "[openfile_background] " << samplename << endl;

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/rootfiles/" +TString::Itoa(tag_year,10) + "/Background/";
  TString filename = root_file_path + cyclename + "_" + samplename + ".root";
  TFile *current_file = new TFile ((filename)) ;
  gDirectory->cd(dir);

  TH1F * current_hist = (TH1F*)gDirectory -> Get(histname);
  if(current_hist){
    current_hist -> SetDirectory(0);
  }
  TH1::AddDirectory(kFALSE);
  
  //if(samplename.Contains("DYJets") && current_hist) current_hist->Scale(DY_scale);
  mapfunc[histname + cyclename + samplename] = current_hist;
  
  current_file -> Close();
  delete current_file;
}

void openfile_DATA(TString cyclename, TString samplename, TString dir, TString histname){

  cout << "[openfile_DATA] " << samplename << endl;
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/rootfiles/" +TString::Itoa(tag_year,10) + "/DATA/";
  TString filename = root_file_path + cyclename + "_" + samplename + ".root";
  TFile *current_file = new TFile ((filename)) ;
  gDirectory->cd(dir);
  
  TH1F * current_hist = (TH1F*)gDirectory -> Get(histname);
  if(current_hist){
    current_hist -> SetDirectory(0);
  }
  TH1::AddDirectory(kFALSE);

  mapfunc[histname + cyclename + samplename] = current_hist;

  cout << "[openfile_DATA] " << histname + cyclename + samplename << endl;
  current_file -> Close();
  delete current_file;
}
 
void openfile_signal(TString samplename, TString channel, TString dir, TString histname){

  TString filename = "SR_ZpNN_HNPairToJJJJ_" + channel + "_" + samplename + "_WR5000.root";
  cout << "[[openfile_signal]] Open " << filename << endl;

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/rootfiles/" +TString::Itoa(tag_year,10) + "/Signal/";
  filename = root_file_path + filename;
  TFile *current_file = new TFile ((filename)) ;
  gDirectory->cd(dir);

  TH1F * current_hist = (TH1F*)gDirectory -> Get(histname);
  if(current_hist){
    current_hist -> SetDirectory(0);
  }
  TH1::AddDirectory(kFALSE);

  mapfunc[histname + samplename] = current_hist;

  current_file -> Close();
  delete current_file;
}
   
void draw_syst_lines(TString current_histname, int N_bin, double x1_template, double x2_template){
  
  TString legend_list[4] = {"Other", "Non-prompt", "ttbar", "DYJets"};
  TString func = current_histname + "_central";
  TString nameofhistogram = current_histname + "_central";
  func.Insert(0, "ratio_");
  std::vector<double> central_content;
  central_content.clear();
  for(int j = 1; j < N_bin; j++){
    double current_content = mapfunc[func+ "rebin"] -> GetBinContent(j);
    central_content.push_back(current_content);
  }
  
  
  for(int i_legend = 0; i_legend < 4; i_legend++){
    TString current_sample = map_sample_names[legend_list[4 - i_legend - 1]].at(0);
    TString cycle_and_sample = Cycle_name + current_sample;
    for(int i_syst = 0; i_syst < N_syst_comparison; i_syst++){
      TString current_syst = systematics_comparison[i_syst];
      map_gr[current_histname + current_sample + current_syst + "line_Up"] = new TGraph();
      map_gr[current_histname + current_sample + current_syst + "line_Up"] = new TGraph();
      map_gr[current_histname + current_sample + current_syst + "line_Down"] = new TGraph();
      
      for(int j = 1; j < N_bin - 1; j++){
	cout << "[draw_syst_lines] getting syst arrays : " << current_sample << ", " << current_syst << ", " << j << endl;
	
	if(map_syst_array[current_histname + cycle_and_sample + current_syst + "Up"].size() > 0){
	  double current_content = mapfunc[current_histname + "_central" + cycle_and_sample + "rebin"] -> GetBinContent(j+1);
	  double current_up = fabs(map_syst_array[current_histname + cycle_and_sample + current_syst + "Up"].at(j) - current_content)/ central_content.at(j);
	  double current_down = (-1.) * fabs(map_syst_array[current_histname + cycle_and_sample + current_syst + "Down"].at(j) - current_content)/ central_content.at(j);
	  cout << "current_up : " << current_up << ", current_down : " << current_down << endl;
	  double current_x = mapfunc[func + "rebin"] -> GetBinCenter(j + 1);
	  map_gr[current_histname + current_sample + current_syst + "line_Up"]  -> SetPoint(j, current_x, current_up);
	  map_gr[current_histname + current_sample + current_syst + "line_Down"]  -> SetPoint(j, current_x, current_down);
	}
      }
    }
  }
  TString canvas = current_histname + "syst_line";
  mapcanvas[canvas] = new TCanvas(canvas,"",800,800);
  canvas_margin(mapcanvas[canvas]);
  gStyle -> SetOptStat(1111);
  
  mapfunc["template_syst_line" + nameofhistogram] = (TH1F*)GetHist(func + "rebin") -> Clone("template_syst_line" + nameofhistogram);
  for(int i = 1; i < N_bin; i++){
    mapfunc["template_syst_line" + nameofhistogram] -> SetBinContent(i, 1.);
    mapfunc["template_syst_line" + nameofhistogram] -> SetBinError(i, 0.);
  }

  TString name_x = nameofhistogram;

  mapfunc["template_syst_line" + nameofhistogram] -> SetTitle("");
  mapfunc["template_syst_line" + nameofhistogram] -> SetLineColor(kWhite);
  mapfunc["template_syst_line" + nameofhistogram] -> GetXaxis() -> SetTitle(name_x);
  mapfunc["template_syst_line" + nameofhistogram] -> GetXaxis() -> SetTitleSize(0.05);
  mapfunc["template_syst_line" + nameofhistogram] -> GetXaxis() -> SetLabelSize(0.05);
  mapfunc["template_syst_line" + nameofhistogram] -> GetXaxis() -> SetRangeUser(x1_template, x2_template);
  mapfunc["template_syst_line" + nameofhistogram] -> GetYaxis() -> SetTitle("#frac{Syst.}{Pred.}");
  mapfunc["template_syst_line" + nameofhistogram] -> GetYaxis() -> SetTitleSize(0.05);
  mapfunc["template_syst_line" + nameofhistogram] -> GetYaxis() -> SetTitleOffset(1.0);
  mapfunc["template_syst_line" + nameofhistogram] -> GetYaxis() -> SetLabelSize(0.05);
  mapfunc["template_syst_line" + nameofhistogram] -> GetYaxis() -> SetNdivisions(505);
  mapfunc["template_syst_line" + nameofhistogram] -> SetMinimum(-0.5);
  mapfunc["template_syst_line" + nameofhistogram] -> SetMaximum(0.5);
  mapfunc["template_syst_line" + nameofhistogram] -> SetStats(0);
  mapfunc["template_syst_line" + nameofhistogram] -> Draw("hist");

  TString legend = current_histname + "syst_line";
  maplegend[legend] = new TLegend(0.30, 0.75, 0.90, 0.92);
  
  Int_t line_array[] = {1, 2, 3, 4};
  Int_t color_array[] = {632, 901, 800, 400, 416, 432, 600, 880, 883,
			 879, 878, 877, 876, 875, 874, 873, 872, 871, 870,
  };
  int N_legend_columns = 4;
  bool is_Wjet = true;
  for(int i_legend = 0; i_legend < 4; i_legend++){
    TString current_sample = map_sample_names[legend_list[4 - i_legend - 1]].at(0);
    TString cycle_and_sample = Cycle_name + current_sample;
    for(int i_syst = 0; i_syst < N_syst_comparison; i_syst++){
      TString current_syst = systematics_comparison[i_syst];

      
      if(map_syst_array[current_histname + cycle_and_sample + current_syst + "Up"].size() > 0){
	map_gr[current_histname + current_sample + current_syst + "line_Up"] -> SetLineColor(color_array[i_syst]);
	map_gr[current_histname + current_sample + current_syst + "line_Up"] -> SetLineStyle(line_array[i_legend]);
	map_gr[current_histname + current_sample + current_syst + "line_Up"] -> Draw("lsame");
	
	map_gr[current_histname + current_sample + current_syst + "line_Down"] -> SetLineColor(color_array[i_syst]);
	map_gr[current_histname + current_sample + current_syst + "line_Down"] -> SetLineStyle(line_array[i_legend]);
	map_gr[current_histname + current_sample + current_syst + "line_Down"] -> Draw("lsame");
	
	maplegend[legend] -> AddEntry(map_gr[current_histname + current_sample + current_syst + "line_Up"], legend_list[4 - i_legend - 1] + "_" + current_syst, "l");
      }
      else is_Wjet = false;
    }
  }
  if(is_Wjet) maplegend[legend] -> SetNColumns(4);
  else maplegend[legend] -> SetNColumns(3);

  maplegend[legend] -> Draw("same");
  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_Lumi.SetTextSize(0.035);
  TString total_lumi = "41.3";
  if(tag_year == 2016) total_lumi = "35.9";
  if(tag_year == 2018) total_lumi = "59.7";
  latex_Lumi.DrawLatex(0.7, 0.96, total_lumi + " fb^{-1} (13 TeV)");

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
  pdfname.Append("_syst_lines.pdf");

  mapcanvas[canvas] -> SaveAs(pdfname);

}
void make_histogram(TString nameofhistogram, TString current_histname, int N_bin, double binx[]){

  //cout << "[make_histogram] nameofhistogram : " << nameofhistogram << ", current_histname : " << current_histname << endl;
 
  if( GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets"].at(0))){

     
    /*
    Int_t idx0 = nameofhistogram.Index("__DYreweight");
    TString new_histname = "";
    for(int i_index = 0; i_index < idx0; i_index++){
      new_histname += nameofhistogram[i_index];
    }
    new_histname += "_central";
   
    if(GetHist(new_histname + Cycle_name + map_sample_names["DYJets"].at(0)) ){
      double integ_before = GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets"].at(0) ) -> Integral();
      double integ_after  = GetHist(new_histname + Cycle_name + map_sample_names["DYJets"].at(0)) -> Integral();
      double DY_norm = integ_before / integ_after;
      //cout << "integ_before : " << integ_before << endl;
      //cout << "integ_after : " << integ_after << endl;
      //cout << "DY_norm : " << DY_norm << endl;
      GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets"].at(0)) -> Scale(DY_norm);
    }
    */
    
  }
  
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
  double DY_norm = 1.;
  if( GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets"].at(0))){
      
    N_DY = GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets"].at(0)) -> Integral();
  
    cout << "N_DY " << map_sample_names["DYJets"].at(0) << " : " << N_DY << endl;
    
    if(GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets"].at(0)))     N_MC += GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets"].at(0)) -> Integral();
    if(GetHist(nameofhistogram + Cycle_name + map_sample_names["ttbar"].at(0)) ) N_MC += GetHist(nameofhistogram + Cycle_name + map_sample_names["ttbar"].at(0)) -> Integral();
    if(GetHist(nameofhistogram + Cycle_name + map_sample_names["Non-prompt"].at(0)) )  N_MC += GetHist(nameofhistogram + Cycle_name + map_sample_names["Non-prompt"].at(0)) -> Integral();
    if(GetHist(nameofhistogram + Cycle_name + map_sample_names["Other"].at(0)) )        N_MC += GetHist(nameofhistogram + Cycle_name + map_sample_names["Other"].at(0)) -> Integral();
    
    if(GetHist(nameofhistogram + Cycle_name + current_data) )N_Data = GetHist(nameofhistogram + Cycle_name + current_data) -> Integral();
    
    DY_norm = (N_DY + N_Data - N_MC) / N_DY;
    if(tag_year == 2016){
      if(nameofhistogram.Contains("0AK8")){
        if(nameofhistogram.Contains("DiMu")){
          DY_norm = 0.812;
        }
        if(nameofhistogram.Contains("DiEle")){
          DY_norm = 0.962;
        }
      }
      if(nameofhistogram.Contains("1AK8")){
        if(nameofhistogram.Contains("DiMu")){
          DY_norm = 0.709;
        }
        if(nameofhistogram.Contains("DiEle")){
          DY_norm = 0.689;
        }
      }
      if(nameofhistogram.Contains("2AK8")){
        if(nameofhistogram.Contains("DiMu")){
          DY_norm = 0.604;
        }
        if(nameofhistogram.Contains("DiEle")){
          DY_norm = 0.511;
        }
      }
    }
    if(tag_year == 2017){
      if(nameofhistogram.Contains("0AK8")){
        if(nameofhistogram.Contains("DiMu")){
          DY_norm = 1.095;
        }
        if(nameofhistogram.Contains("DiEle")){
          DY_norm = 1.093;
        }
      }
      if(nameofhistogram.Contains("1AK8")){
        if(nameofhistogram.Contains("DiMu")){
          DY_norm = 1.056;
        }
        if(nameofhistogram.Contains("DiEle")){
          DY_norm = 1.062;
        }
      }
      if(nameofhistogram.Contains("2AK8")){
        if(nameofhistogram.Contains("DiMu")){
          DY_norm = 0.747;
	}
        if(nameofhistogram.Contains("DiEle")){
          DY_norm = 0.881;
	}
      }
    }
    if(tag_year == 2018){
      if(nameofhistogram.Contains("0AK8")){
        if(nameofhistogram.Contains("DiMu")){
          DY_norm = 1.127;
        }
        if(nameofhistogram.Contains("DiEle")){
          DY_norm = 1.087;
        }
      }
      if(nameofhistogram.Contains("1AK8")){
        if(nameofhistogram.Contains("DiMu")){
          DY_norm = 0.947;
        }
        if(nameofhistogram.Contains("DiEle")){
          DY_norm = 0.952;
        }
      }
      if(nameofhistogram.Contains("2AK8")){
        if(nameofhistogram.Contains("DiMu")){
          DY_norm = 0.997;
        }
        if(nameofhistogram.Contains("DiEle")){
          DY_norm = 0.957;
        }
      }
    }

    GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets"].at(0)) -> Scale(DY_norm); 
    
    outfile << nameofhistogram << "'s DY_norm : " << DY_norm << endl;
    
    for(int i_syst = 1; i_syst < N_syst; i_syst++){
      if(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["DYJets"].at(0))){
	GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["DYJets"].at(0)) -> Scale(DY_norm);
      }
    }
  }
  

  TString title_y;
  if(nameofhistogram.Contains("N")) title_y = "Number";
  else title_y = "Events/bin";

  bool blind = false;
  if((nameofhistogram.Contains("SR")) && (!nameofhistogram.Contains("EMu"))) blind = blind_SR;

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
  mappad[pad1] -> SetLogy();
  
  maplegend[legend] = new TLegend(0.60, 0.60, 0.96, 0.92);

  if(debug) cout << "1" << endl;
  
  maphstack[hstack] = new THStack(hstack, "Stacked_" + nameofhistogram);
  gStyle->SetOptTitle(0);

  TString legend_list[4] = {"Other", "Non-prompt", "ttbar", "DYJets"};
  //TString legend_list[5] = {"Other", "Non-prompt", "ttbar", "DYJets", "QCD"}; 
  Int_t colour_array[] = {419, 901, 416, 400, 600};

  TString name_cycle = nameofhistogram + Cycle_name;
  if(debug) cout << "check1" << endl;
  
  // -- Add Overflow bin
  TString overflow = "overflow";
  mapfunc[func] = new TH1F("", "", N_bin - 1, binx);
  mapfunc[func + "rebin"] = new TH1F("", "", N_bin - 1, binx);
  mapfunc[func + "rebin_stat_err"] = new TH1F("", "", N_bin - 1, binx);
  mapfunc[func + "blind_data"] = new TH1F("", "", N_bin - 1, binx);
  double x1_template = binx[0];
  double x2_template = binx[N_bin - 1];
  // -- For each Bkg Legend : Add overflow, Rebin, Syst error bar, Add on legend, and set fill & line colors
  for(int i_legend = 0; i_legend < 4; i_legend++){
    if(debug) cout << legend_list[i_legend] << endl;
    TString current_sample = map_sample_names[legend_list[i_legend]].at(0);
    if(mapfunc[nameofhistogram + Cycle_name + current_sample]){
      
      Rebin_with_overflow(nameofhistogram + Cycle_name + current_sample, N_bin, binx);
      Remove_Negative_Bins(nameofhistogram + Cycle_name + current_sample + "rebin");
      //cout << "GetHist  integral " << current_sample << " : " << GetHist(nameofhistogram + Cycle_name + current_sample)->Integral() << endl;
      //cout << "temp integral " << current_sample << " : " << htmp->Integral() << endl;
      
      
      if(debug) cout << current_sample + " : calling all syst hists" << endl;
      // -- Make rebinned hists for all systematic categories
      for(int i_syst = 0; i_syst < N_syst_comparison; i_syst++){
	Save_syst_array(current_histname, systematics_comparison[i_syst], Cycle_name + current_sample, N_bin, binx);
      }
      
      mapfunc[nameofhistogram + Cycle_name + current_sample + "rebin"] ->  SetFillColor(colour_array[i_legend]);
      mapfunc[nameofhistogram + Cycle_name + current_sample + "rebin"] ->  SetLineColor(colour_array[i_legend]);
      
      // -- Prepare blind data point (= bkg)
      mapfunc[func + "blind_data"]  -> Add(GetHist(nameofhistogram + Cycle_name + current_sample + "rebin"));
      
      maphstack[hstack] -> Add(GetHist(nameofhistogram + Cycle_name + current_sample + "rebin"));
      //maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + current_sample + "rebin"), legend_list[i_legend], "lf");
      mapfunc[func + "rebin"] -> Add(GetHist(nameofhistogram + Cycle_name + current_sample + "rebin"));
      mapfunc[func + "rebin_stat_err"] -> Add(GetHist(nameofhistogram + Cycle_name + current_sample + "rebin"));
    }
  }

  for(int i_legend = 0; i_legend < 4; i_legend++){
    TString current_sample = map_sample_names[legend_list[4 - i_legend - 1]].at(0);
    if(GetHist(nameofhistogram + Cycle_name + current_sample + "rebin")) {
      maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + current_sample + "rebin"), legend_list[4 - i_legend - 1], "lf");
      sum_syst_error(current_histname, Cycle_name + current_sample, N_bin);
    }  
  }

  vector<double> bkg_central;
  vector<double> bkg_up;
  vector<double> bkg_down;
  bkg_central.clear();
  bkg_up.clear();
  bkg_down.clear();
  for(int j = 1; j < N_bin ; j++){
    double current_bkg = mapfunc[func + "rebin"] -> GetBinContent(j);
    bkg_central.push_back(current_bkg);
    
    double current_up_sum = 0.;
    double current_down_sum = 0.;
    
    for(int i_legend = 0; i_legend < 4; i_legend++){
      TString current_sample = map_sample_names[legend_list[4 - i_legend - 1]].at(0);
      if(GetHist(nameofhistogram + Cycle_name + current_sample + "rebin")) {
	cout << "[make_histogram] current_sample : " << current_sample << endl;
	cout << "[make_histogram] Up size : " << map_syst_array[current_histname + Cycle_name + current_sample + "Up"].size() << endl;
	cout << "[make_histogram] Down size : " << map_syst_array[current_histname + Cycle_name + current_sample + "Down"].size() << endl;
	
	current_up_sum   = current_up_sum   + map_syst_array[current_histname + Cycle_name + current_sample + "Up"].at(j-1);
	current_down_sum = current_down_sum + map_syst_array[current_histname + Cycle_name + current_sample + "Down"].at(j-1);
      }
    }
    bkg_up.push_back(current_up_sum);
    bkg_down.push_back(current_down_sum);
  }


  map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error"] = new TGraphAsymmErrors(mapfunc[func + "rebin"]);
  map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error_ratio"] = new TGraphAsymmErrors(mapfunc[func + "rebin"]);
  
  for(int i = 0; i < N_bin - 1; i++){
    map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error"] -> SetPointEYlow(i,  bkg_down.at(i));
    map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error"] -> SetPointEYhigh(i, bkg_up.at(i));
    
    double current_x = mapfunc[func + "rebin"] -> GetBinCenter(i + 1);
    double current_y = mapfunc[func + "rebin"] -> GetBinContent(i + 1);
    map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error_ratio"] -> SetPoint(i, current_x, 1.);
    map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error_ratio"] -> SetPointEYlow(i, bkg_down.at(i) / current_y);
    map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error_ratio"] -> SetPointEYhigh(i, bkg_up.at(i) / current_y);
  }
    
  //draw_syst_lines(current_histname, N_bin, x1_template, x2_template);
  mappad[pad1] -> cd();

  if(debug) cout << "4" << endl;
  cout << "nameofhistogram + Cycle_name + current_data : " << nameofhistogram + Cycle_name + current_data << endl;
  Rebin_with_overflow(nameofhistogram + Cycle_name + current_data, N_bin, binx);

  if(blind){
    mapfunc[nameofhistogram + Cycle_name + current_data + "rebin"] = (TH1F*)GetHist(func + "blind_data") -> Clone(clone + "blind");
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
  map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error"] -> SetFillColor(kBlack);;
  map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error"] -> SetFillStyle(3013);
  map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error"] -> SetMarkerSize(0);
  map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error"] -> Draw("e2same");
  maplegend[legend] -> AddEntry(map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error"], "Syst. + Stat.", "pf");
  
  // -- Draw Signal
  Int_t signal_colour_array[] = {622, 632, 432, 600};
  /*
  for(unsigned int i_signal = 0; i_signal < map_sample_names["Signal"].size(); i_signal++){
    Rebin_with_overflow(nameofhistogram + map_sample_names["Signal"].at(i_signal), N_bin, binx);
    mapfunc[nameofhistogram + map_sample_names["Signal"].at(i_signal) + "rebin"] -> SetLineColor(signal_colour_array[i_signal]);
    mapfunc[nameofhistogram + map_sample_names["Signal"].at(i_signal) + "rebin"] -> SetLineWidth(2);
    maplegend[legend] -> AddEntry(mapfunc[nameofhistogram + map_sample_names["Signal"].at(i_signal) + "rebin"],  map_sample_names["Signal"].at(i_signal), "l");
    mapfunc[nameofhistogram + map_sample_names["Signal"].at(i_signal) + "rebin"] -> Draw("histsame");
  }
  */
  // -- Set y-axis range
  double data_max = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetMaximum();
  double data_min = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetMinimum();
  maphstack[hstack] -> SetMaximum(data_max * 1000.);//logy
  //maphstack[hstack] -> SetMaximum(data_max * 1.5);
  maphstack[hstack] -> SetMinimum(0.01);
  
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
  mapfunc["pad2_template" + nameofhistogram] = (TH1F*)GetHist(func + "rebin") -> Clone(clone);
  mapfunc["stat_err" + nameofhistogram] = (TH1F*)GetHist(func + "rebin_stat_err") -> Clone(clone);
  
  for(int i = 1; i < N_bin; i++){
    mapfunc["pad2_template" + nameofhistogram] -> SetBinContent(i, 1.);
    mapfunc["pad2_template" + nameofhistogram] -> SetBinError(i, 0.);
  }
  if(debug) cout << "7" << endl;

  // -- Draw Line at Data/MC = 1
  mapline[line] = new TLine(binx[0], 1, binx[N_bin - 1], 1);
  mapline[line] -> Draw();
  mapline[line] -> SetLineStyle(1);
  mapline[line] -> SetLineColor(kBlue);
  
  // -- Draw syst box hist with all contents as 1
  TString name_x = nameofhistogram;
    
  mapfunc["pad2_template" + nameofhistogram] -> SetTitle("");
  mapfunc["pad2_template" + nameofhistogram] -> SetLineColor(kWhite);
  mapfunc["pad2_template" + nameofhistogram] -> GetXaxis() -> SetTitle(name_x);
  mapfunc["pad2_template" + nameofhistogram] -> GetXaxis() -> SetTitleSize(0.15);
  mapfunc["pad2_template" + nameofhistogram] -> GetXaxis() -> SetLabelSize(0.10);
  mapfunc["pad2_template" + nameofhistogram] -> GetXaxis() -> SetRangeUser(x1_template, x2_template);
  mapfunc["pad2_template" + nameofhistogram] -> GetYaxis() -> SetTitle("#frac{Obs.}{Pred.}");
  mapfunc["pad2_template" + nameofhistogram] -> GetYaxis() -> SetTitleSize(0.15);
  mapfunc["pad2_template" + nameofhistogram] -> GetYaxis() -> SetTitleOffset(0.5);
  mapfunc["pad2_template" + nameofhistogram] -> GetYaxis() -> SetLabelSize(0.08);
  mapfunc["pad2_template" + nameofhistogram] -> GetYaxis() -> SetNdivisions(505);
  mapfunc["pad2_template" + nameofhistogram] -> SetMinimum(0.5);
  mapfunc["pad2_template" + nameofhistogram] -> SetMaximum(1.5);
  mapfunc["pad2_template" + nameofhistogram] -> SetStats(0);
  mapfunc["pad2_template" + nameofhistogram] -> Draw("hist");

  //map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error"] -> SetFillColor(kBlack);;
  //map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error"] -> SetFillStyle(3013);
  //map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error"] -> SetMarkerSize(0);
  //map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error"] -> Draw("e2same");
  
  map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error_ratio"] -> SetFillColor(kOrange);
  map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error_ratio"] -> SetMarkerSize(0);
  map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error_ratio"] -> SetLineColor(kWhite);
  //map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error_ratio"] -> SetFillStyle(3005);
  //map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error_ratio"] -> SetFillColorAlpha(45,0.35);
  
  if(debug) cout << "8" << endl;

  for(int i = 0; i < N_bin; i++){
    bool no_data = false;
    double data_value = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetBinContent(i + 1);;
    if(data_value < 10e-7) no_data = true;
    if(no_data){
      map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error_ratio"] -> SetPointEYlow(i, 0.);
      map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error_ratio"] -> SetPointEYhigh(i, 0.);
    }
    cout << "bkg_error_ration " << i << " error : " << map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error_ratio"] -> GetErrorYhigh(i) << ", " << map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error_ratio"] -> GetErrorYlow(i) << endl;
  }
  map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error_ratio"] -> Draw("e2same");
  

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
  maplegend["bottom" + legend]->AddEntry(map_asym_gr[nameofhistogram + Cycle_name + "Bkg_Error_ratio"], "Stat. + Syst.", "f");
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
  TString total_lumi = "41.3";
  if(tag_year == 2016) total_lumi = "35.9";
  if(tag_year == 2018) total_lumi = "59.7";
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
  pdfname.Append("_signal.pdf");
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
  
  maphist.clear();
  map_syst_array.clear();
  mapfunc.clear();
  if(debug) cout << "[[draw_histogram]] histname : " << histname << endl;  
  
  unsigned int N_bin = map_bin_vector[histname].size();
  double current_bins[100];
  for(unsigned int i_bin = 0; i_bin < N_bin; i_bin++){
    current_bins[i_bin] = map_bin_vector[histname].at(i_bin);
  }


  // -- Get which region(CR1, CR2 , ...) and channel (DiEle, DiMu, ...)
  //int N_regions = 6;
  int N_regions = 7;
  
  TString regions[] = {"CR_Zmass",
		       "CR_Zmass_nPV40bigger",
		       "CR_Zmass_nPV40smaller",
                       //"CR_ttbar",
                       //"CR_inv_mll",
                       "tight_CR_Zmass",
                       //"tight_CR_ttbar",
                       "SR",
		       "SR_b_veto",
		       "SR_1b",
  };


  int N_channels = 6;
  TString channels[] = {"DiEle",
                        "DiMu",
                        "EMu",
			"DiEle_DYreweight",
                        "DiMu_DYreweight",
                        "EMu_DYreweight",
  };

  int N_AK8_bins = 3;
  TString AK8_bins[] = {"0AK8",
			"1AK8",
			"2AK8"
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
  if(!histname.Contains("mZp")){
    for(int i_AK8_bin = 0; i_AK8_bin < N_AK8_bins; i_AK8_bin++){
      if(histname.Contains(AK8_bins[i_AK8_bin])){
	current_dir = AK8_bins[i_AK8_bin] + "_" + current_dir;
      } 
    }
  }
  
  if(histname.Contains("Problematic")){
    current_dir = "Problematic_" + current_dir;
  }

  if(debug) cout << "current_dir : " << current_dir << "/" << histname << "+systflags" << endl;
  if(current_dir.Contains("empty") || current_channel.Contains("empty")) return;
  
  
  for(int i_syst = 0; i_syst < N_syst; i_syst++){
    // -- FIXME
    //TString current_dir_syst = current_dir + "_" + systematics[i_syst];
    TString current_hist_syst = histname + "_" + systematics[i_syst];
    TString current_dir_syst = current_dir + "_" + systematics[i_syst];
    //TString current_hist_syst = histname + "_DYreweight_" + systematics[i_syst]; // -- FIXME
    TString legend_list[4] = {"DYJets", "ttbar", "Non-prompt", "Other"};
    //TString legend_list[5] = {"QCD", "DYJets", "ttbar", "Non-prompt", "Other"}; 
    for(int i_legend_list = 0; i_legend_list < 4; i_legend_list++){
      
      unsigned int vector_size = map_sample_names[legend_list[i_legend_list]].size();
      //cout << "legend_list[i_legend_list] : " << legend_list[i_legend_list] << endl;
      //cout << "vector_size : " << vector_size << endl;
      
      for(int i_vector = 0; i_vector < vector_size; i_vector++){
	openfile_background(Cycle_name, map_sample_names[legend_list[i_legend_list]].at(i_vector), current_dir_syst, current_hist_syst);
      }
    }

    //cout << "current_dir_syst : " << current_dir_syst << endl;
    //cout << "current_hist_syst : " << current_hist_syst << endl;
    openfile_DATA(Cycle_name, DoubleEG, current_dir_syst, current_hist_syst);
    openfile_DATA(Cycle_name, SingleMuon, current_dir_syst, current_hist_syst);
    
    if(systematics[i_syst] == "central"){
      TString current_signal_channel = "";
      if(current_hist_syst.Contains("DiMu")) current_signal_channel = "MuMu";
      if(current_hist_syst.Contains("DiEle")) current_signal_channel = "EE";

      for(unsigned int i_signal = 0; i_signal < map_sample_names["Signal"].size(); i_signal++){
	openfile_signal(map_sample_names["Signal"].at(i_signal), current_signal_channel, current_dir_syst, current_hist_syst);
      }
    }

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
      
      //open_files(current_histname);
      open_files(current_histname + "_DYreweight");
      map_bin_vector.clear();
      
    }//while not end of file
  }//end if file open
  data_file.close();
  
}

void QuickPlot_CR_signal(int year=2018){
  setTDRStyle();
  
  // == Set Input Sample List
  tag_year = year;
  
  outfile.open("output_" + TString::Itoa(tag_year,10) + ".txt");
  
  if(tag_year == 2018){
    //map_sample_names["DYJets"] = {"DYJets_MG"};
    map_sample_names["DYJets"] = {"DYJets_MG_HT"};
    map_sample_names["ttbar"] = {"TT_powheg"};
    //map_sample_names["WJets"] = {"WJets_MG_HT"};
    map_sample_names["Non-prompt"] = {"fake"};
    map_sample_names["Other"] = {"VV"};
    map_sample_names["QCD"] = {"QCD"};
    map_sample_names["Muon"] = {"data_SingleMuon"};
    map_sample_names["EGamma"] = {"data_EGamma"};
    
    SingleMuon = "data_SingleMuon";
    DoubleEG = "data_EGamma";
  }

  if(tag_year == 2017){
    //map_sample_names["DYJets"] = {"DYJets_MG"};
    map_sample_names["DYJets"] = {"DYJets_MG_HT"};
    map_sample_names["ttbar"] = {"TT_powheg"};
    //map_sample_names["WJets"] = {"WJets_MG_HT"};
    map_sample_names["Non-prompt"] = {"fake"};
    map_sample_names["Other"] = {"VV"};
    map_sample_names["QCD"] = {"QCD"};
    map_sample_names["Muon"] = {"data_SingleMuon"};
    map_sample_names["EGamma"] = {"data_DoubleEG"};

    SingleMuon = "data_SingleMuon";
    DoubleEG = "data_DoubleEG";
  }

  if(tag_year == 2016){
    //map_sample_names["DYJets"] = {"DYJets_MG"};
    map_sample_names["DYJets"] = {"DYJets_MG_HT"};
    map_sample_names["ttbar"] = {"TT_powheg"};
    //map_sample_names["WJets"] = {"WJets_MG_HT"};
    map_sample_names["Non-prompt"] = {"fake"};
    map_sample_names["Other"] = {"VV"};
    map_sample_names["QCD"] = {"QCD"};
    map_sample_names["Muon"] = {"data_SingleMuon"};
    map_sample_names["EGamma"] = {"data_DoubleEG"};

    SingleMuon = "data_SingleMuon";
    DoubleEG = "data_DoubleEG";
  }
  
  //map_sample_names["Signal"] = {"ZP1000_N200", "ZP2000_N500", "ZP3000_N1000", "ZP4000_N1500"};
  map_sample_names["Signal"] = {""};

  Cycle_name = "HN_pair_all_SkimTree_LRSMHighPt";

  
  //open_binning_file("binning_limit_postfit.txt");
  
  /*
  open_binning_file("binning_limit_merged_mZp_2AK8_SR_DiEle.txt");
  open_binning_file("binning_limit_merged_mZp_1AK8_SR_DiEle.txt");
  open_binning_file("binning_limit_merged_mZp_0AK8_SR_DiEle.txt");
  */
  //open_binning_file("binning_limit_merged_mZp_2AK8_SR_DiMu.txt");
  //open_binning_file("binning_limit_merged_mZp_1AK8_SR_DiMu.txt");
  //open_binning_file("binning_limit_merged_mZp_0AK8_SR_DiMu.txt");
  
  //open_binning_file("binning_2017_comparion.txt");
  open_binning_file("binning_problematic.txt");
  //open_binning_file("binning_signal_vs_data.txt");
  
  outfile.close();
  
}
