#include "canvas_margin.h"
#include "mylib.h"


void openfile_signal(TString cyclename, TString samplename, TString histname, TString systematics){

  cout << "[openfile_background] " << samplename << endl;
  //cout << "samplename : " << samplename << endl;
  //cout << "histname : " << histname << endl;
  //cout << "cyclename : " << cyclename << endl;
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/output/LimitTemplate/" +TString::Itoa(tag_year,10) + "/";
  TString filename = root_file_path + "shape_" + histname + "_HNPairToJJJJ_" + samplename + "_" + cyclename + "_WR5000.root";
    

  TFile *current_file = new TFile ((filename)) ;
  TString current_hist_name = histname + "_HNPairToJJJJ_" + samplename + "_" + cyclename + "_WR5000_" + systematics;
  if(systematics == "") current_hist_name = histname + "_HNPairToJJJJ_" + samplename + "_" + cyclename + "_WR5000";
  TH1F * current_hist = (TH1F*)gDirectory -> Get(current_hist_name);
  if(current_hist){
    current_hist -> SetDirectory(0);
  }
  TH1::AddDirectory(kFALSE);

  cout << "current_hist_name : " <<current_hist_name << endl;
  mapfunc[histname + cyclename + samplename + systematics] = current_hist;
  current_file -> Close();
  delete current_file;

}

void make_histogram(TString cyclename, TString samplename, TString histname, int N_bin, double binx[]){
  
  TString systematics_PDF[3] = {"", "PDFUp", "PDFDown"};
  TString nameofhistogram = histname;
  int year = tag_year;
  
  TString title_y;
  if(nameofhistogram.Contains("N")) title_y = "Number";
  else title_y = "Events/bin";

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
  //mapcanvas[canvas] -> SetLogy();

  maplegend[legend] = new TLegend(0.60, 0.80, 0.96, 0.95);

  if(debug) cout << "[make_histogram] Check Point 1" << endl;

  Int_t signal_colour_array[] = {622, 632, 432, 600, 416};
  TString current_hist = histname + cyclename + samplename;
  TString current_hist_up = current_hist + "PDFUp";
  TString current_hist_down = current_hist + "PDFDown";

  cout <<"current_hist : " <<current_hist << endl;
  Rebin_with_overflow(current_hist,N_bin, binx);
  Rebin_with_overflow(current_hist_up,N_bin, binx);
  Rebin_with_overflow (current_hist_down,N_bin, binx);
  
  mapfunc[func + "template"] = (TH1F*)mapfunc[current_hist + "rebin"] -> Clone(func + "template");
  mapfunc[func + "template"]  -> SetStats(0);
  double hist_max = mapfunc[func + "template"] -> GetMaximum();
  mapfunc[func + "template"] -> SetMaximum(hist_max * 2.5);
  //mapfunc[func + "template"] -> SetMaximum(hist_max * 100);//logy
  mapfunc[func + "template"] -> SetMinimum(0.1);
  mapfunc[func + "template"] -> GetXaxis() -> SetTitle(nameofhistogram);
  mapfunc[func + "template"] -> GetYaxis()->SetLabelSize(0.05);;
  mapfunc[func + "template"] -> GetYaxis()->SetTitleSize(0.07);;
  mapfunc[func + "template"] -> GetYaxis()->SetTitleOffset(1.02);;
  mapfunc[func + "template"] -> GetYaxis()->SetTitle("Events/bin");
  mapfunc[func + "template"] -> Draw("hist");
  
  mapfunc[current_hist + "rebin"] -> SetLineColor(kBlack);
  mapfunc[current_hist_up + "rebin"] -> SetLineColor(kBlue);
  mapfunc[current_hist_down + "rebin"] -> SetLineColor(kRed);
  mapfunc[current_hist + "rebin"] -> SetLineWidth(2);
  mapfunc[current_hist_up + "rebin"] -> SetLineWidth(2);
  mapfunc[current_hist_down + "rebin"] -> SetLineWidth(2);
  
  mapfunc[current_hist + "rebin"] -> Draw("histsame");
  mapfunc[current_hist_up + "rebin"] -> Draw("histsame");
  mapfunc[current_hist_down + "rebin"] -> Draw("histsame");

  maplegend[legend] -> AddEntry(mapfunc[current_hist + "rebin"], cyclename + " " + samplename + " central", "l");
  maplegend[legend] -> AddEntry(mapfunc[current_hist_up + "rebin"], cyclename + " " + samplename + " PDF Up", "l");
  maplegend[legend] -> AddEntry(mapfunc[current_hist_down + "rebin"], cyclename + " " + samplename + " PDF Down", "l");

  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kBlack);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(1001);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> Draw("same");

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

  pdfname.Append(cyclename);
  pdfname.Append("_");
  pdfname.Append(nameofhistogram);
  pdfname.Append("_PDF.pdf");
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

  int N_regions = 7;

  TString regions[] = {"CR_Zmass",
                       "CR_Zmass_nPV40bigger",
                       "CR_Zmass_nPV40smaller",
		       "tight_CR_Zmass",
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

  TString systematics_PDF[3] = {"", "PDFUp", "PDFDown"};
  for(int i_syst = 0; i_syst < 3; i_syst++){
    TString current_hist_syst = histname + "_" + systematics_PDF[i_syst];
    TString current_dir_syst = current_dir + "_" + systematics_PDF[i_syst];

    TString current_signal_channel = "";
    if(current_hist_syst.Contains("DiMu")) current_signal_channel = "MuMu";
    if(current_hist_syst.Contains("DiEle")) current_signal_channel = "EE";
    
    for(unsigned int i_signal = 0; i_signal < map_sample_names["Signal"].size(); i_signal++){
      openfile_signal(map_sample_names["Signal"].at(i_signal), current_signal_channel, histname, systematics_PDF[i_syst]);
    }
  }
  for(unsigned int i_signal = 0; i_signal < map_sample_names["Signal"].size(); i_signal++){
    TString current_signal_channel = "";
    if(histname.Contains("DiMu")) current_signal_channel = "MuMu";
    if(histname.Contains("DiEle")) current_signal_channel = "EE";
    make_histogram(map_sample_names["Signal"].at(i_signal), current_signal_channel, histname, N_bin, current_bins);
  }

}

void open_binning_file(TString filename){
  
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

      open_files(current_histname + "_DYreweight");
      map_bin_vector.clear();
    }
  }

  data_file.close();

}

void Plot_signal_syst_shape(int year=2018){
  setTDRStyle();

  tag_year = year;

  map_sample_names["Signal"] = {"ZP2000_N200", "ZP3200_N600", "ZP4400_N600"};
  map_sample_names["syst"] = {"PDF", ""};
  Cycle_name = "HN_pair_all_SkimTree_LRSMHighPt";

  open_binning_file("binning_limit_merged_mZp_2AK8_SR_DiEle.txt");
  open_binning_file("binning_limit_merged_mZp_1AK8_SR_DiEle.txt");
  open_binning_file("binning_limit_merged_mZp_0AK8_SR_DiEle.txt");
  open_binning_file("binning_limit_merged_mZp_2AK8_SR_DiMu.txt");
  open_binning_file("binning_limit_merged_mZp_1AK8_SR_DiMu.txt");
  open_binning_file("binning_limit_merged_mZp_0AK8_SR_DiMu.txt");


}
