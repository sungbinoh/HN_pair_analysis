#include "canvas_margin.h"
#include "mylib.h"


void openfile(TString channel, TString which_fit, TString histname, TString mass_point){
  //if(debug) cout << "[openfile_background] " << samplename << endl;

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString filename = WORKING_DIR+"/rootfiles/fitDiagnostics/fitDiagnostics_fit_result_" + mass_point + ".root";
  //TString filename = WORKING_DIR+"/rootfiles/fitDiagnostics/fitDiagnostics_fit_result.root";

  TFile *current_file = new TFile ((filename)) ;
  gDirectory->cd(which_fit);
  gDirectory->cd(histname + "_" + TString::Itoa(tag_year,10));

  //TString histograms[4] = {"DYJets_MG_HT", "TT_powheg", "VV", "WJets_MG_HT"};
  TString histograms[4] = {"DYJets_MG_HT", "TT_powheg", "VV", "fake"};


  for(int i_hist = 0; i_hist < 4; i_hist++){
    TH1F * current_hist = (TH1F*)gDirectory -> Get(histname + "_DYreweight_" + histograms[i_hist]);
    if(current_hist){
      current_hist -> SetDirectory(0);
    }
    TH1::AddDirectory(kFALSE);

    mapfunc[histname + histograms[i_hist] + which_fit] = current_hist;
  }
  
  TGraphAsymmErrors * current_data = (TGraphAsymmErrors*) gDirectory -> Get("data");
  if(current_data){
    //current_data -> SetDirectory(0);
  }
  //TGraphAsymmErrors::AddDirectory(kFALSE);
  map_asym_gr[histname + "data" + which_fit] = current_data;

  current_file -> Close();
  delete current_file;
}

void make_histogram(TString histname, TString which_fit, int N_bin, double binx[], TString mass_point){

  if(debug){
    cout << "[[make_histogram]] checking binning arrary" << endl;
    for(int i = 0; i < N_bin; i++){
      cout << binx[i] << ", " << endl;
    }
  }
  
  TString nameofhistogram = histname + "_" + which_fit;
  
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
  mappad[pad1] -> SetLogy();
  
  maplegend[legend] = new TLegend(0.60, 0.60, 0.96, 0.92);

  if(debug) cout << "1" << endl;
  

  TH1F * pad1_template = new TH1F("", "", N_bin-1, binx);
  gStyle->SetOptTitle(0);

  // -- Bkgs
  //TString histograms[4] = {"VV", "WJets_MG_HT", "TT_powheg", "DYJets_MG_HT"};
  //TString legend_list[4] = {"Other", "WJets", "ttbar", "DYJets"};
  TString histograms[4] = {"VV", "fake", "TT_powheg", "DYJets_MG_HT"};
  TString legend_list[4] = {"Other", "Non-prompt", "ttbar", "DYJets"};


  Int_t colour_array[] = {419, 901, 416, 400};

  maphstack[hstack] = new THStack(hstack, "Stacked_" + nameofhistogram);
  TH1F * bkg_sum = new TH1F("", "", N_bin-1, binx);
  
  for(int i_hist = 0; i_hist < 4; i_hist++){
    TString current_hist = histname + histograms[i_hist] + which_fit;
    if(mapfunc[current_hist]){
      mapfunc[current_hist + "template"] = new TH1F(current_hist + "template", current_hist + "template", N_bin-1, binx); 
      for(int i_bin = 1; i_bin < N_bin; i_bin++){
	double current_content = mapfunc[current_hist] -> GetBinContent(i_bin);
	double current_error = mapfunc[current_hist] -> GetBinError(i_bin);
	mapfunc[current_hist + "template"]  -> SetBinContent(i_bin, current_content);
	
	
	double bkg_sum_content = bkg_sum -> GetBinContent(i_bin);
	double bkg_sum_error = bkg_sum -> GetBinError(i_bin);
	bkg_sum -> SetBinContent(i_bin, bkg_sum_content + current_content);
	bkg_sum -> SetBinError(i_bin, bkg_sum_error + current_error);
      }
  
      mapfunc[current_hist + "template"] -> SetFillColor(colour_array[i_hist]);
      mapfunc[current_hist + "template"] -> SetLineColor(colour_array[i_hist]);
      maphstack[hstack] -> Add(mapfunc[current_hist + "template"]);
     
    }  
  }

  for(int i_legend = 0; i_legend < 4; i_legend++){
    TString current_hist = histname + histograms[4 - i_legend - 1] + which_fit;
    if(mapfunc[current_hist]){
      maplegend[legend] -> AddEntry(mapfunc[current_hist + "template"] , legend_list[4 - i_legend - 1], "lf");
    }
  }
  
  // -- Data
  TH1F * data_hist = new TH1F("", "", N_bin-1, binx);
  double max_data = 0.;
  const int N_points = N_bin - 1;
  double data_y[N_points], data_x[N_points];
  double data_ey_high[N_points], data_ey_low[N_points];
  double data_ex_high[N_points], data_ex_low[N_points];
  for(int i_point = 0; i_point < N_points; i_point++){
    double current_x, current_y;
    double current_ey_high, current_ey_low;
    map_asym_gr[histname + "data" + which_fit] -> GetPoint(i_point, current_x, current_y);
    current_ey_high = map_asym_gr[histname + "data" + which_fit] -> GetErrorYhigh(i_point);
    current_ey_low = map_asym_gr[histname + "data" + which_fit] -> GetErrorYlow(i_point);
    //cout << "current_y : " << current_y << ", +" << current_ey_high << ", -" << current_ey_low << endl;
    data_y[i_point] = current_y;
    data_ey_high[i_point] = current_ey_high;
    data_ey_low[i_point] = current_ey_low;
    data_hist -> SetBinContent(i_point + 1, current_y); 
    if(current_y > max_data) max_data = current_y;
  }
  for(int i_x = 0; i_x < N_points; i_x++){
    double current_ex = (binx[i_x + 1] - binx[i_x]) / 2.0;
    double current_x = binx[i_x] + current_ex;
    
    data_x[i_x] = current_x;
    data_ex_high[i_x] = current_ex; 
    data_ex_low[i_x] =current_ex;
  }
   
  TGraphAsymmErrors * data_gr = new TGraphAsymmErrors(N_points, data_x, data_y, data_ex_low, data_ex_high, data_ey_low, data_ey_high);
  
  // -- Draw Pad
  pad1_template -> Draw("hist");
  pad1_template -> GetYaxis()->SetLabelSize(0.05);;
  pad1_template -> GetYaxis()->SetTitleSize(0.07);;
  pad1_template -> GetYaxis()->SetTitleOffset(1.02);;
  pad1_template -> GetXaxis()->SetLabelSize(0);
  pad1_template -> GetXaxis() -> SetTitle(nameofhistogram);
  //pad1_template -> GetYaxis() -> SetRangeUser(0., max_data * 1.8);
  pad1_template -> GetYaxis() -> SetRangeUser(0.1, max_data * 1000.); // logy
  pad1_template -> GetYaxis() -> SetTitle(title_y);
  pad1_template -> SetStats(0);
  pad1_template -> Draw("histsame");
  
  maphstack[hstack] -> Draw("histsame");
  
  bkg_sum -> SetFillColor(kBlack);
  bkg_sum -> SetFillStyle(3013);
  bkg_sum -> SetMarkerSize(0);;
  bkg_sum -> Draw("e2same");;

  data_gr -> Draw("epsame");
  
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kBlack);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(1001);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> Draw("same");

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

  TString name_x = nameofhistogram;
  
  TH1F * pad2_template = new TH1F("", "", N_bin-1, binx);
  pad2_template -> Draw("hist");

  pad2_template -> SetTitle("");
  pad2_template -> SetLineColor(kWhite);
  pad2_template -> GetXaxis() -> SetTitle(name_x);
  pad2_template -> GetXaxis() -> SetTitleSize(0.15);
  pad2_template -> GetXaxis() -> SetLabelSize(0.10);
  pad2_template -> GetYaxis() -> SetTitle("#frac{Obs.}{Pred.}");
  pad2_template -> GetYaxis() -> SetTitleSize(0.15);
  pad2_template -> GetYaxis() -> SetTitleOffset(0.5);
  pad2_template -> GetYaxis() -> SetLabelSize(0.08);
  pad2_template -> GetYaxis() -> SetNdivisions(505);
  pad2_template -> GetYaxis() -> SetRangeUser(0.5, 1.5);
  pad2_template -> SetStats(0);
  pad2_template -> Draw("histsame");
  
  TH1F * pad2_bkg_err = new TH1F("", "", N_bin-1, binx);

  double ratio_data_y[N_points], ratio_data_ey_high[N_points], ratio_data_ey_low[N_points];
  double ratio_bkg_ey[N_points];
  for(int i_bin = 1; i_bin < N_bin; i_bin++){
    double current_data = data_hist -> GetBinContent(i_bin);
    double current_bkg = bkg_sum -> GetBinContent(i_bin);
    double current_data_ey_high = data_ey_high[i_bin - 1] / current_data;
    double current_data_ey_low = data_ey_low[i_bin - 1] / current_data;
    double current_bkg_ey = bkg_sum -> GetBinError(i_bin) / current_bkg;
    
    data_hist -> SetBinContent(i_bin, current_data/current_bkg);
    data_hist -> SetBinError(i_bin, 0.);
    ratio_data_y[i_bin - 1] = current_data/current_bkg;
    ratio_data_ey_high[i_bin - 1] = current_data_ey_high;
    ratio_data_ey_low[i_bin - 1] = current_data_ey_low;
    ratio_bkg_ey[i_bin -1] = current_bkg_ey;
    
    pad2_bkg_err -> SetBinContent(i_bin, 1.0);
    pad2_bkg_err -> SetBinError(i_bin, current_bkg_ey);
  }
  
  pad2_bkg_err -> SetFillColor(kOrange);
  pad2_bkg_err -> SetMarkerSize(0);
  pad2_bkg_err -> SetLineColor(kWhite);
  pad2_bkg_err -> Draw("E2same");
  
  data_hist -> Draw("epsame");
  
  TGraphAsymmErrors * ratio_data_gr = new TGraphAsymmErrors(N_points, data_x, ratio_data_y, data_ex_low, data_ex_high, ratio_data_ey_low, ratio_data_ey_high);
  ratio_data_gr -> Draw("epsame");
  
  maplegend["bottom" + legend] = new TLegend(0.2, 0.85, 0.4, 0.95);
  maplegend["bottom" + legend]->SetBorderSize(0);
  maplegend["bottom" + legend]->SetNColumns(3);
  maplegend["bottom" + legend]->AddEntry(pad2_bkg_err, "Stat. + Syst.", "f");
  maplegend["bottom" + legend]->AddEntry(data_hist , "Obs./Pred.", "p");
  maplegend["bottom" + legend]->Draw("same");

  TLine *pad2_line = new TLine(binx[0], 1, binx[N_bin - 1], 1);
  pad2_line -> SetLineStyle(1);
  pad2_line -> SetLineColor(kBlue);
  pad2_line -> Draw("same");

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
  pdfname.Append("_" + mass_point);
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


void open_files(TString histname, TString mass_point){
  
  maphist.clear();
  map_syst_array.clear();
  mapfunc.clear();
  if(debug) cout << "[[draw_histogram]] histname : " << histname << endl;  
  
  unsigned int N_bin = map_bin_vector[histname].size();
  double current_bins[100];
  for(unsigned int i_bin = 0; i_bin < N_bin; i_bin++){
    current_bins[i_bin] = map_bin_vector[histname].at(i_bin);
  }
  
  TString channel = "";
  if(histname.Contains("DiEle")) channel = "DiEle";
  else if(histname.Contains("DiMu") || histname.Contains("EMu") ) channel = "DiMu";
  else channel = "NULL";

  openfile(channel, "shapes_prefit", histname, mass_point);
  openfile(channel, "shapes_fit_b", histname, mass_point);
  make_histogram(histname, "shapes_prefit", N_bin, current_bins, mass_point);
  make_histogram(histname, "shapes_fit_b", N_bin, current_bins, mass_point);

}

void open_binning_file(TString filename, TString mass_point){
 
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
      }

      cout << "[[open_binning_file]] current_histname : " << current_histname << endl;
      
      open_files(current_histname, mass_point);
      map_bin_vector.clear();
      
    }//while not end of file
  }//end if file open
  data_file.close();
  
}

void Plot_Post_fit_signal(int int_tag_year=2018){
  setTDRStyle();
  
  tag_year = int_tag_year;
  TString binning_file = "binning_limit_postfit.txt";

  cout << "binning_file : " << binning_file << endl;
  
  open_binning_file(binning_file);
  
  open_binning_file("binning_limit_merged_mZp_2AK8_SR_DiEle.txt", "HNPairToJJJJ_EE_ZP4000_N1200_WR5000");
  /*
  open_binning_file("binning_limit_merged_mZp_1AK8_SR_DiEle.txt");
  open_binning_file("binning_limit_merged_mZp_0AK8_SR_DiEle.txt");
  open_binning_file("binning_limit_merged_mZp_2AK8_SR_DiMu.txt");
  open_binning_file("binning_limit_merged_mZp_1AK8_SR_DiMu.txt");
  open_binning_file("binning_limit_merged_mZp_0AK8_SR_DiMu.txt");
  */

}
