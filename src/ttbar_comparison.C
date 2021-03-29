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
    cout << "[openfile_background] " << histname + cyclename + samplename + TString::Itoa(year,10) << " exists" << endl;
    current_hist -> SetDirectory(0);
  }
  TH1::AddDirectory(kFALSE);

  mapfunc[histname + cyclename + samplename + TString::Itoa(year,10)] = current_hist;
  //mapfunc[histname + cyclename + samplename + TString::Itoa(year,10)] -> SetLineColor(kRed);
  current_file -> Close();
  delete current_file;
}


void open_files(TString histname){

  if(debug) cout << "[[draw_histogram]] histname : " << histname << endl;

  unsigned int N_bin = map_bin_vector[histname].size();
  double current_bins[100];
  for(unsigned int i_bin = 0; i_bin < N_bin; i_bin++){
    current_bins[i_bin] = map_bin_vector[histname].at(i_bin);
  }

  int N_regions = 5;

  TString regions[] = {"CR_Zmass",
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

  for(int i_syst = 0; i_syst < 1; i_syst++){

    TString current_hist_syst = histname + "_" + systematics[i_syst];
    TString current_dir_syst = current_dir + "_" + systematics[i_syst];

    TString legend_list[1] = {"ttbar"};

    for(int i_legend_list = 0; i_legend_list < 1; i_legend_list++){

      unsigned int vector_size = map_sample_names[legend_list[i_legend_list]].size();
      for(int i_vector = 0; i_vector < vector_size; i_vector++){
        openfile_background(Cycle_name, map_sample_names[legend_list[i_legend_list] + "2016"].at(i_vector), current_dir_syst, current_hist_syst, 2016);
        openfile_background(Cycle_name, map_sample_names[legend_list[i_legend_list] + "2017"].at(i_vector), current_dir_syst, current_hist_syst, 2017);
        openfile_background(Cycle_name, map_sample_names[legend_list[i_legend_list] + "2018"].at(i_vector), current_dir_syst, current_hist_syst, 2018);
      }
      TString histname_2016 = current_hist_syst + Cycle_name + map_sample_names[legend_list[i_legend_list] + "2016"].at(0) + "2016";
      TString histname_2017 = current_hist_syst + Cycle_name + map_sample_names[legend_list[i_legend_list] + "2017"].at(0) + "2017";
      TString histname_2018 = current_hist_syst + Cycle_name + map_sample_names[legend_list[i_legend_list] + "2018"].at(0) + "2018";
      
      cout << "histname_2016 : " << histname_2016 << endl;
      
      Rebin_with_overflow(histname_2016, N_bin, current_bins);
      Rebin_with_overflow(histname_2017, N_bin, current_bins);
      Rebin_with_overflow(histname_2018, N_bin, current_bins);
    }
  }

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

      TString current_histname = ((TObjString *)(tx->At(0)))->String();
      if(debug) cout << "[[open_binning_file]] current_histname : " << current_histname << endl;
      map_bin_vector[current_histname].clear();
      map_bin_vector[current_histname + "_DYreweight"].clear();

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


void overlap_ttbar(TString NAK8){
  //mZp_0AK8_SR_DiMu_DYreweight_centralHN_pair_all_SkimTree_LRSMHighPtTT_powheg2016
  TString nameofhistogram = "mZp_" + NAK8;
  TString title_y;
  if(nameofhistogram.Contains("N")) title_y = "Number";
  else title_y = "Events/bin";

  TString canvas = nameofhistogram;
  TString legend = nameofhistogram;
  TString func = nameofhistogram;
  TString clone = nameofhistogram;
  TString line = nameofhistogram;
  TString pad1 = nameofhistogram;
  TString pad2 = nameofhistogram;
  canvas.Insert(0, "c_");
  pad1.Insert(0, "pad1_");
  pad2.Insert(0, "pad2_");
  legend.Insert(0, "legend_");
  func.Insert(0, "ratio_");
  clone.Insert(0, "h3_");
  line.Insert(0, "l_");

  mapcanvas[canvas] = new TCanvas(canvas,"",800,800);
  canvas_margin(mapcanvas[canvas]);
  gStyle -> SetOptStat(1111);
  //mapcanvas[canvas] -> SetLogy();

  ////////////////////////////////////
  /// -- Pad 1
  ////////////////////////////////////
  mappad[pad1] = new TPad("", "", 0, 0.25, 1, 1);
  mappad[pad1] -> SetTopMargin( 0.07 );
  mappad[pad1] -> SetBottomMargin( 0.02 );
  mappad[pad1] -> SetLeftMargin( 0.15 );
  mappad[pad1] -> SetRightMargin( 0.03 );
  mappad[pad1] -> Draw();
  mappad[pad1] -> cd();
  mappad[pad1] -> SetLogy();
  
  maplegend[legend] = new TLegend(0.70, 0.60, 0.95, 0.95);

  TString histname_ee = "mZp_" + NAK8 + "_SR_DiEle_DYreweight_centralHN_pair_all_SkimTree_LRSMHighPtTT_powheg";
  TString histname_mumu = "mZp_" + NAK8 + "_SR_DiMu_DYreweight_centralHN_pair_all_SkimTree_LRSMHighPtTT_powheg";
  TString histname_emu = "mZp_" + NAK8 + "_SR_EMu_DYreweight_centralHN_pair_all_SkimTree_LRSMHighPtTT_powheg";
  
  mapfunc[histname_ee] = (TH1F*)mapfunc[histname_ee + "2016rebin"] -> Clone(histname_ee);
  mapfunc[histname_ee] -> Add(mapfunc[histname_ee + "2017rebin"]);
  mapfunc[histname_ee] -> Add(mapfunc[histname_ee + "2018rebin"]);

  mapfunc[histname_mumu] = (TH1F*)mapfunc[histname_mumu + "2016rebin"] -> Clone(histname_mumu);
  mapfunc[histname_mumu] -> Add(mapfunc[histname_mumu + "2017rebin"]);
  mapfunc[histname_mumu] -> Add(mapfunc[histname_mumu + "2018rebin"]);

  mapfunc[histname_emu] = (TH1F*)mapfunc[histname_emu + "2016rebin"] -> Clone(histname_emu);
  mapfunc[histname_emu] -> Add(mapfunc[histname_emu + "2017rebin"]);
  mapfunc[histname_emu] -> Add(mapfunc[histname_emu + "2018rebin"]);
  mapfunc[histname_emu] -> Scale(0.5);
  
  mapfunc[histname_ee + "ratio"] = (TH1F*)mapfunc[histname_ee] -> Clone(histname_ee + "ratio"); 
  mapfunc[histname_mumu + "ratio"] = (TH1F*)mapfunc[histname_mumu] -> Clone(histname_mumu + "ratio");
  mapfunc[histname_emu + "ratio"] = (TH1F*)mapfunc[histname_emu] -> Clone(histname_emu + "ratio");

  mapfunc[func + "template"] = (TH1F*)mapfunc[histname_emu] -> Clone(func + "template");
  mapfunc[func + "template"]  -> SetStats(0);
  double hist_max = mapfunc[func + "template"] -> GetMaximum();
  //mapfunc[func + "template"] -> SetMaximum(hist_max * 1.5);
  mapfunc[func + "template"] -> SetMaximum(hist_max * 10);//logy
  mapfunc[func + "template"] -> SetMinimum(0.1);
  mapfunc[func + "template"] -> GetXaxis() -> SetTitle(nameofhistogram);
  mapfunc[func + "template"] -> GetXaxis() -> SetLabelSize(0);
  mapfunc[func + "template"] -> GetYaxis()->SetLabelSize(0.05);;
  mapfunc[func + "template"] -> GetYaxis()->SetTitleSize(0.07);;
  mapfunc[func + "template"] -> GetYaxis()->SetTitleOffset(1.02);;
  mapfunc[func + "template"] -> GetYaxis()->SetTitle("Events/bin");
  mapfunc[func + "template"] -> SetFillColor(0);
  mapfunc[func + "template"] -> Draw("ep");

  mapfunc[histname_ee] -> SetLineColor(kBlack);
  mapfunc[histname_mumu] -> SetLineColor(kBlue);
  mapfunc[histname_emu] -> SetLineColor(kRed);

  mapfunc[histname_ee] -> SetMarkerColor(kBlack);
  mapfunc[histname_mumu] -> SetMarkerColor(kBlue);
  mapfunc[histname_emu] -> SetMarkerColor(kRed);


  mapfunc[histname_ee] -> SetLineWidth(2);
  mapfunc[histname_mumu] -> SetLineWidth(2);
  mapfunc[histname_emu] -> SetLineWidth(2);

  mapfunc[histname_ee] -> Draw("epsame");
  mapfunc[histname_mumu] -> Draw("epsame");
  mapfunc[histname_emu]-> Draw("epsame");

  maplegend[legend] -> AddEntry(mapfunc[histname_ee], "ee", "l");
  maplegend[legend] -> AddEntry(mapfunc[histname_mumu], "mumu", "l");
  maplegend[legend] -> AddEntry(mapfunc[histname_emu], "emu", "l");

  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kBlack);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(1001);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> Draw("same");

  ////////////////////////////////////
  // --- Pad 2
  ////////////////////////////////////
  mapcanvas[canvas] -> cd();
  mapcanvas[canvas] -> cd();
  mappad[pad2] = new TPad(pad2, "", 0, 0, 1, 0.25);
  mappad[pad2] -> SetTopMargin( 0.03 );
  mappad[pad2] -> SetBottomMargin( 0.4 );
  mappad[pad2] -> SetLeftMargin( 0.15 );
  mappad[pad2] -> SetRightMargin( 0.03 );
   mappad[pad2] -> Draw();
  mappad[pad2] -> cd();

  mapfunc[func + "ratio"] = (TH1F*)mapfunc[histname_emu] -> Clone(func + "ratio");
  mapfunc[func + "ratio"] -> SetTitle("");
  mapfunc[func + "ratio"] -> SetLineColor(kWhite);
  //mapfunc[func + "ratio"] -> GetXaxis() -> SetTitle(name_x);
  mapfunc[func + "ratio"] -> GetXaxis() -> SetTitleSize(0.15);
  mapfunc[func + "ratio"] -> GetXaxis() -> SetLabelSize(0.10);
  mapfunc[func + "ratio"] -> GetXaxis() -> SetRangeUser(0., 6000.);
  mapfunc[func + "ratio"] -> GetYaxis() -> SetTitle("#frac{Obs.}{Pred.}");
  mapfunc[func + "ratio"] -> GetYaxis() -> SetTitleSize(0.15);
  mapfunc[func + "ratio"] -> GetYaxis() -> SetTitleOffset(0.5);
  mapfunc[func + "ratio"] -> GetYaxis() -> SetLabelSize(0.08);
  mapfunc[func + "ratio"] -> GetYaxis() -> SetNdivisions(505);
  mapfunc[func + "ratio"] -> SetMinimum(0.5);
  mapfunc[func + "ratio"] -> SetMaximum(1.5);
  mapfunc[func + "ratio"] -> SetStats(0);
  mapfunc[func + "ratio"] -> Draw("hist");

  double norm_emu = mapfunc[histname_emu + "ratio"]  -> Integral();
  double norm_ee = mapfunc[histname_ee + "ratio"]  -> Integral();
  double norm_mumu = mapfunc[histname_mumu + "ratio"]  -> Integral();

  cout << "norm_emu : " <<norm_emu << endl;

  mapfunc[histname_ee + "ratio"] -> Scale(1./norm_ee);
  mapfunc[histname_emu + "ratio"] -> Scale(1./norm_emu);
  mapfunc[histname_mumu + "ratio"] -> Scale(1./norm_mumu);

  mapfunc[histname_ee + "ratio"] -> Divide(mapfunc[histname_emu + "ratio"] );
  mapfunc[histname_mumu + "ratio"] -> Divide(mapfunc[histname_emu + "ratio"] );
  
  cout << "mapfunc[histname_ee ]-> GetBinContent(3) : " << mapfunc[histname_ee + "ratio"] -> GetBinContent(3) << endl;

  mapfunc[histname_ee + "ratio"] -> SetLineColor(kBlack);
  mapfunc[histname_mumu + "ratio"] -> SetLineColor(kBlue);

  mapfunc[histname_ee + "ratio"] -> Draw("epsame");
  mapfunc[histname_mumu + "ratio"] -> Draw("epsame");

  mapline[line] = new TLine(0., 1., 6000., 1.);
  mapline[line] -> Draw();
  mapline[line] -> SetLineStyle(1);
  mapline[line] -> SetLineColor(kRed);
  
  maplegend["bottom" + legend] = new TLegend(0.2, 0.85, 0.6, 0.95);
  maplegend["bottom" + legend]->SetBorderSize(0);
  maplegend["bottom" + legend]->SetNColumns(3);
  maplegend["bottom" + legend]->AddEntry(mapfunc[histname_ee + "ratio"], "ee/emu (normalized)", "l");
  maplegend["bottom" + legend]->AddEntry(mapfunc[histname_mumu + "ratio"], "mumu/emu (normalized)", "l");
  maplegend["bottom" + legend]->Draw("same");

  
  mapcanvas[canvas] -> cd();
  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_Lumi.SetTextSize(0.035);
  latex_Lumi.DrawLatex(0.7, 0.96, "137.4 fb^{-1} (13 TeV)");

  TString pdfname;
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  pdfname = WORKING_DIR + "/plots/Signal/ttbar_comparison_" + NAK8 + ".pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);



}



void ttbar_comparison(int year=2019){
  
  setTDRStyle();
  map_sample_names["ttbar2018"] = {"TT_powheg"};
  map_sample_names["ttbar2017"] = {"TT_powheg"};
  map_sample_names["ttbar2016"] = {"TT_powheg"};
  map_sample_names["ttbar"] = {"TT_powheg"};
  
  Cycle_name = "HN_pair_all_SkimTree_LRSMHighPt";

  TString binning_file = "binning_limit_postfit.txt";
  cout << "binning_file : " << binning_file << endl;

  open_binning_file(binning_file);
  open_binning_file("binning_limit_merged_mZp_2AK8_SR_DiEle.txt");
  open_binning_file("binning_limit_merged_mZp_1AK8_SR_DiEle.txt");
  open_binning_file("binning_limit_merged_mZp_0AK8_SR_DiEle.txt");
  open_binning_file("binning_limit_merged_mZp_2AK8_SR_DiMu.txt");
  open_binning_file("binning_limit_merged_mZp_1AK8_SR_DiMu.txt");
  open_binning_file("binning_limit_merged_mZp_0AK8_SR_DiMu.txt");

  overlap_ttbar("0AK8");
  overlap_ttbar("1AK8");
  overlap_ttbar("2AK8");

  outfile.close();


}
