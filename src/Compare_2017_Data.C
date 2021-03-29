#include "canvas_margin.h"
#include "mylib.h"

void openfile_DATA(TString cyclename, TString samplename, TString dir, TString histname){

  cout << "[openfile_DATA] " << samplename << endl;
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/rootfiles/" +TString::Itoa(tag_year,10) + "/DATA/";
  TString filename = root_file_path + cyclename + "_data_" + samplename + ".root";
  TFile *current_file = new TFile ((filename)) ;
  gDirectory->cd(dir);

  TH1F * current_hist = (TH1F*)gDirectory -> Get(histname);
  if(current_hist){
    current_hist -> SetDirectory(0);
  }
  TH1::AddDirectory(kFALSE);

  mapfunc[histname + cyclename + samplename] = current_hist;
  cout << "[openfile_DATA] histname + cyclename + samplename : " << histname + cyclename + samplename << endl;
  current_file -> Close();
  delete current_file;
}

void openfile_single_period(TString cyclename, TString samplename, TString dir, TString histname){

  cout << "[openfile_single_period] " << samplename << endl;
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
  cout << "[openfile_DATA] histname + cyclename + samplename : " << histname + cyclename + samplename << endl;
  current_file -> Close();
  delete current_file;
}

void print_contents(TString cyclename, TString samplename, TString histname){
  TString hist_all = histname + cyclename + samplename + "rebin";
  TString hist_B = histname + cyclename + samplename + "_Brebin";
  TString hist_C = histname + cyclename + samplename + "_Crebin";
  TString hist_D = histname + cyclename + samplename + "_Drebin";
  TString hist_E = histname + cyclename + samplename + "_Erebin";
  TString hist_F = histname + cyclename + samplename + "_Frebin";

  cout << "all : " << mapfunc[hist_all] -> GetBinContent(5) << endl;
  cout << "B : " << mapfunc[hist_B] -> GetBinContent(5) << ",\t" << (mapfunc[hist_B] -> GetBinContent(5)) / 4.823 << " per fb-1" << endl;
  cout << "C : " << mapfunc[hist_C] -> GetBinContent(5) << ",\t" << (mapfunc[hist_C] -> GetBinContent(5)) / 9.664 << " per fb-1" << endl;
  cout << "D : " << mapfunc[hist_D] -> GetBinContent(5) << ",\t" << (mapfunc[hist_D] -> GetBinContent(5)) / 4.252 << " per fb-1" << endl;
  cout << "E : " << mapfunc[hist_E] -> GetBinContent(5) << ",\t" << (mapfunc[hist_E] -> GetBinContent(5)) / 9.278 << " per fb-1" << endl;
  cout << "F : " << mapfunc[hist_F] -> GetBinContent(5) << ",\t" << (mapfunc[hist_F] -> GetBinContent(5)) / 13.540 << " per fb-1" << endl;
  //4.823 + 9.664 + 4.252 + 9.278 + 13.540
  

}

void overlap_data(TString cyclename, TString samplename, TString histname){

  TString hist_BCD = histname + cyclename + samplename + "_BCDrebin";
  TString hist_EF = histname + cyclename + samplename + "_EFrebin";
  cout << "hist_BCD : " <<hist_BCD << endl;
  
  TString nameofhistogram = histname;
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

  ////////////////////////////////////
  // -- Pad 1
  ////////////////////////////////////
  mappad[pad1] = new TPad("", "", 0, 0.25, 1, 1);
  mappad[pad1] -> SetTopMargin( 0.07 );
  mappad[pad1] -> SetBottomMargin( 0.02 );
  mappad[pad1] -> SetLeftMargin( 0.15 );
  mappad[pad1] -> SetRightMargin( 0.03 );
  mappad[pad1] -> Draw();
  mappad[pad1] -> cd();
  //mappad[pad1] -> SetLogy();

  //double norm_BCD = mapfunc[hist_BCD] -> Integral();
  //double norm_EF = mapfunc[hist_EF] -> Integral();
  double norm_BCD = 18.739000;
  double norm_EF = 22.818000;
  mapfunc[hist_BCD] -> Scale(1./norm_BCD);
  mapfunc[hist_EF] -> Scale(1./norm_EF);

  maplegend[legend] = new TLegend(0.70, 0.60, 0.95, 0.95);
  mapfunc[func + "template"] = (TH1F*)mapfunc[hist_BCD] -> Clone(func + "template");
  mapfunc[func + "template"]  -> SetStats(0);
  double hist_max = mapfunc[func + "template"] -> GetMaximum();
  mapfunc[func + "template"] -> SetMaximum(hist_max * 1.5); 
  //mapfunc[func + "template"] -> SetMaximum(hist_max * 10);//logy
  mapfunc[func + "template"] -> SetMinimum(0.001);
  mapfunc[func + "template"] -> GetXaxis() -> SetTitle(nameofhistogram);
  mapfunc[func + "template"] -> GetXaxis() -> SetLabelSize(0);
  mapfunc[func + "template"] -> GetYaxis()->SetLabelSize(0.05);;
  mapfunc[func + "template"] -> GetYaxis()->SetTitleSize(0.07);;
  mapfunc[func + "template"] -> GetYaxis()->SetTitleOffset(1.02);;
  mapfunc[func + "template"] -> GetYaxis()->SetTitle("Events/bin");
  mapfunc[func + "template"] -> SetFillColor(0);
  mapfunc[func + "template"] -> Draw("ep");
  
  mapfunc[hist_BCD] -> SetLineColor(kBlue);
  mapfunc[hist_EF] -> SetLineColor(kRed);

  mapfunc[hist_BCD] -> Draw("epsame");
  mapfunc[hist_EF] -> Draw("epsame");

  maplegend[legend] -> AddEntry(mapfunc[hist_BCD], "Period BCD", "l");
  maplegend[legend] -> AddEntry(mapfunc[hist_EF], "Period EF", "l");
  
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kBlack);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(1001);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> Draw("same");

  ////////////////////////////////////
  // -- Pad 2
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

  mapfunc[func + "ratio"] = (TH1F*)mapfunc[hist_BCD] -> Clone(func + "ratio");
  mapfunc[func + "ratio"] -> SetTitle("");
  mapfunc[func + "ratio"] -> SetLineColor(kWhite);
  mapfunc[func + "ratio"] -> GetXaxis() -> SetTitleSize(0.15);
  mapfunc[func + "ratio"] -> GetXaxis() -> SetLabelSize(0.10);
  mapfunc[func + "ratio"] -> GetXaxis() -> SetRangeUser(0., 6000.);
  mapfunc[func + "ratio"] -> GetYaxis() -> SetTitle("#frac{Obs.}{Pred.}");
  mapfunc[func + "ratio"] -> GetYaxis() -> SetTitleSize(0.15);
  mapfunc[func + "ratio"] -> GetYaxis() -> SetTitleOffset(0.5);
  mapfunc[func + "ratio"] -> GetYaxis() -> SetLabelSize(0.08);
  mapfunc[func + "ratio"] -> GetYaxis() -> SetNdivisions(505);
  mapfunc[func + "ratio"] -> SetMinimum(0.0);
  mapfunc[func + "ratio"] -> SetMaximum(3.0);
  mapfunc[func + "ratio"] -> SetStats(0);
  mapfunc[func + "ratio"] -> Draw("hist");

  mapfunc[hist_BCD + "ratio"] = (TH1F*)mapfunc[hist_BCD] -> Clone(hist_BCD + "ratio");
  mapfunc[hist_BCD + "ratio"] -> Divide(mapfunc[hist_EF] );
  mapfunc[hist_BCD + "ratio"] -> Draw("epsame");

  mapline[line] = new TLine(0., 1., 6000., 1.);
  mapline[line] -> Draw();
  mapline[line] -> SetLineStyle(1);
  mapline[line] -> SetLineColor(kRed);
  
  maplegend["bottom" + legend] = new TLegend(0.2, 0.85, 0.6, 0.95);
  maplegend["bottom" + legend]->SetBorderSize(0);
  maplegend["bottom" + legend]->SetNColumns(3);
  maplegend["bottom" + legend]->AddEntry(mapfunc[hist_BCD + "ratio"], "BCD/EF", "l");
  maplegend["bottom" + legend]->Draw("same");
  
  mapcanvas[canvas] -> cd();
  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_Lumi.SetTextSize(0.035);
  latex_Lumi.DrawLatex(0.7, 0.96, "41.3 fb^{-1} (13 TeV)");

  TString pdfname;
  TString NAK8 = "";
  if(histname.Contains("0AK8")) NAK8 = "0AK8";
  if(histname.Contains("1AK8")) NAK8 = "1AK8";
  if(histname.Contains("2AK8")) NAK8 = "2AK8";

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  pdfname = WORKING_DIR + "/plots/Signal/Data_comparison_" + NAK8 + "_" + histname + ".pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);


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

  if(histname.Contains("Problematic")){
    current_dir = "Problematic_" + current_dir;
  }

  for(int i_syst = 0; i_syst < 1; i_syst++){

    TString current_hist_syst = histname + "_" + systematics[i_syst];
    TString current_dir_syst = current_dir + "_" + systematics[i_syst];

    openfile_DATA(Cycle_name, "SingleMuon", current_dir_syst, current_hist_syst);
    openfile_DATA(Cycle_name, "SingleMuon_BCD", current_dir_syst, current_hist_syst);
    openfile_DATA(Cycle_name, "SingleMuon_EF", current_dir_syst, current_hist_syst);
    openfile_single_period(Cycle_name, "SingleMuon_B", current_dir_syst, current_hist_syst);
    openfile_single_period(Cycle_name, "SingleMuon_C", current_dir_syst, current_hist_syst);
    openfile_single_period(Cycle_name, "SingleMuon_D", current_dir_syst, current_hist_syst);
    openfile_single_period(Cycle_name, "SingleMuon_E", current_dir_syst, current_hist_syst);
    openfile_single_period(Cycle_name, "SingleMuon_F", current_dir_syst, current_hist_syst);

    TString hist_all = current_hist_syst  + Cycle_name + "SingleMuon";
    TString hist_BCD = current_hist_syst  + Cycle_name + "SingleMuon_BCD";
    TString hist_EF = current_hist_syst  + Cycle_name + "SingleMuon_EF";
    TString hist_B = current_hist_syst  + Cycle_name + "SingleMuon_B";
    TString hist_C = current_hist_syst  + Cycle_name + "SingleMuon_C";
    TString hist_D = current_hist_syst  + Cycle_name + "SingleMuon_D";
    TString hist_E = current_hist_syst  + Cycle_name + "SingleMuon_E";
    TString hist_F = current_hist_syst  + Cycle_name + "SingleMuon_F";

    Rebin_with_overflow(hist_all, N_bin, current_bins);
    Rebin_with_overflow(hist_BCD, N_bin, current_bins);
    Rebin_with_overflow(hist_EF, N_bin, current_bins);
    Rebin_with_overflow(hist_B, N_bin, current_bins);
    Rebin_with_overflow(hist_C, N_bin, current_bins);
    Rebin_with_overflow(hist_D, N_bin, current_bins);
    Rebin_with_overflow(hist_E, N_bin, current_bins);
    Rebin_with_overflow(hist_F, N_bin, current_bins);
  }

  for(int i_syst = 0; i_syst < 1; i_syst++){

    TString current_hist_syst = histname + "_" + systematics[i_syst];
    TString current_dir_syst = current_dir + "_" + systematics[i_syst];
    print_contents(Cycle_name, "SingleMuon", current_hist_syst);
    overlap_data(Cycle_name, "SingleMuon", current_hist_syst);
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


void Compare_2017_Data(int year=2017){
  tag_year = year;
  
  setTDRStyle();
  Cycle_name = "HN_pair_all_SkimTree_LRSMHighPt";

  TString binning_file = "binning_limit_postfit.txt";

  //open_binning_file(binning_file);
  open_binning_file("binning_limit_merged_mZp_2AK8_SR_DiMu.txt");
  open_binning_file("binning_limit_merged_mZp_1AK8_SR_DiMu.txt");
  open_binning_file("binning_limit_merged_mZp_0AK8_SR_DiMu.txt");
  open_binning_file("binning_problematic.txt");

}
