#include "canvas_margin.h"
#include "mylib.h"

void openfile_signal(TString samplename, TString channel, TString dir, TString histname, int year){

  TString filename = "SR_ZpNN_HNPairToJJJJ_" + channel + "_" + samplename + "_WR5000.root";
  cout << "[[openfile_signal]] Open " << filename << endl;

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/rootfiles/" +TString::Itoa(year,10) + "/Signal/";
  filename = root_file_path + filename;
  TFile *current_file = new TFile ((filename)) ;
  gDirectory->cd(dir);
  
  TH1F * current_hist = (TH1F*)gDirectory -> Get(histname);
  if(current_hist){
    current_hist -> SetDirectory(0);
    if(debug) cout << "[[openfile_signal]] current_hist TRUE" << endl;
  }
  TH1::AddDirectory(kFALSE);

  mapfunc[histname + samplename + TString::Itoa(year,10)] = current_hist;
  if(debug) cout << "[[openfile_signal]] Called hist : " << histname + samplename + TString::Itoa(year,10) << endl;
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
  
  int year = tag_year;

  TString current_data;
  if(nameofhistogram.Contains("EMu") || nameofhistogram.Contains("DiMu")) current_data = SingleMuon;
  else if(nameofhistogram.Contains("DiEle")) current_data = DoubleEG;
  else return;


  TString title_y;
  if(nameofhistogram.Contains("N")) title_y = "Number";
  else title_y = "Events/bin";

  bool blind = false;
  blind = (nameofhistogram.Contains("SR")) && (!nameofhistogram.Contains("EMu"));

  TString canvas = nameofhistogram;
  TString legend = nameofhistogram;
  TString func = nameofhistogram;
  TString clone = nameofhistogram;
  TString line = nameofhistogram;
  canvas.Insert(0, "c_");
  legend.Insert(0, "legend_");
  func.Insert(0, "ratio_");
  clone.Insert(0, "h3_");
  line.Insert(0, "l_");
  
  mapcanvas[canvas] = new TCanvas(canvas,"",800,800);
  canvas_margin(mapcanvas[canvas]);
  gStyle -> SetOptStat(1111);
  mapcanvas[canvas] -> SetLogy();
  
  maplegend[legend] = new TLegend(0.60, 0.60, 0.96, 0.92);

  if(debug) cout << "[make_histogram] Check Point 1" << endl;
    
  // -- Rebin Histograms
  Int_t signal_colour_array[] = {622, 632, 432, 600, 416};
  
  for(unsigned int i_signal = 0; i_signal < map_sample_names["Signal"].size(); i_signal++){
    Rebin_with_overflow_limit(nameofhistogram + map_sample_names["Signal"].at(i_signal) + TString::Itoa(year,10), N_bin, binx);
  }

  // -- Make Template Histogram and Set Y-max
  mapfunc[func + "template"] = (TH1F*)mapfunc[nameofhistogram + map_sample_names["Signal"].at(3) + TString::Itoa(year,10) + "rebin"] -> Clone(func + "template");
  mapfunc[func + "template"]  -> SetStats(0);
  double hist_max = mapfunc[func + "template"] -> GetMaximum();
  mapfunc[func + "template"] -> SetMaximum(hist_max * 10000);
  mapfunc[func + "template"] -> SetMinimum(0.1);
  mapfunc[func + "template"] -> GetXaxis() -> SetTitle(nameofhistogram);
  mapfunc[func + "template"] -> GetYaxis()->SetLabelSize(0.05);;
  mapfunc[func + "template"] -> GetYaxis()->SetTitleSize(0.07);;
  mapfunc[func + "template"] -> GetYaxis()->SetTitleOffset(1.02);;
  mapfunc[func + "template"] -> GetYaxis()->SetTitle("Events/bin");
  //mapfunc[func + "template"] -> GetXaxis()->SetLabelSize(0);
  mapfunc[func + "template"] -> Draw("hist");
  
  // -- Draw histograms
  for(unsigned int i_signal = 0; i_signal < map_sample_names["Signal"].size(); i_signal++){
    mapfunc[nameofhistogram + map_sample_names["Signal"].at(i_signal) + TString::Itoa(year,10) + "rebin"] -> SetLineColor(signal_colour_array[i_signal]);
    mapfunc[nameofhistogram + map_sample_names["Signal"].at(i_signal) + TString::Itoa(year,10) + "rebin"] -> SetLineWidth(2);
    maplegend[legend] -> AddEntry(mapfunc[nameofhistogram + map_sample_names["Signal"].at(i_signal) + TString::Itoa(year,10) + "rebin"],  map_sample_names["Signal"].at(i_signal), "l");
    mapfunc[nameofhistogram + map_sample_names["Signal"].at(i_signal) + TString::Itoa(year,10) + "rebin"] -> Draw("histsame");
  }

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
  /// -- Latex
  ////////////////////////////////////
  mapcanvas[canvas] -> cd();
  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Simulation}}}");
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
  pdfname.Append("_signals.pdf");
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
  int N_regions = 5;
  
  TString regions[] = {"CR_Zmass",
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
  
  if(debug) cout << "current_dir : " << current_dir << "/" << histname << "+systflags" << endl;
  if(current_dir.Contains("empty") || current_channel.Contains("empty")) return;
  
  
  for(int i_syst = 0; i_syst < N_syst; i_syst++){
    // -- FIXME
    //TString current_dir_syst = current_dir + "_" + systematics[i_syst];
    TString current_hist_syst = histname + "_" + systematics[i_syst];
    TString current_dir_syst = current_dir + "_" + systematics[i_syst];
    //TString current_hist_syst = histname + "_DYreweight_" + systematics[i_syst]; // -- FIXME

    if(systematics[i_syst] == "central"){
      TString current_signal_channel = "";
      if(current_hist_syst.Contains("DiMu")) current_signal_channel = "MuMu";
      if(current_hist_syst.Contains("DiEle")) current_signal_channel = "EE";

      for(unsigned int i_signal = 0; i_signal < map_sample_names["Signal"].size(); i_signal++){
	openfile_signal(map_sample_names["Signal"].at(i_signal), current_signal_channel, current_dir_syst, current_hist_syst, tag_year);
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

void QuickPlot_signal(int year=2018){
  setTDRStyle();
  
  // == Set Input Sample List
  tag_year = year;
  
  map_sample_names["Signal"] = {"ZP1000_N200", "ZP2000_N500", "ZP3000_N1000", "ZP4000_N1500", "ZP4000_N100"};
  
  Cycle_name = "HN_pair_all_SkimTree_LRSMHighPt";
  
  //open_binning_file("binning_uniform_test.txt");
  open_binning_file("binning_signal_vs_data.txt");
  //open_binning_file("binning_test.txt");

  outfile.close();
  
}









