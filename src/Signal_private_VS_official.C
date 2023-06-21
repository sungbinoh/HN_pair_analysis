#include "canvas_margin.h"
#include "mylib.h"

void openfile(TString samplename, TString dir, TString histname){
  
  cout << "[openfile] " << samplename << endl;
  
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/rootfiles/" +TString::Itoa(tag_year,10) + "/Signal/";
  TString filename = root_file_path + Cycle_name + "_" + samplename + ".root";
  TFile *current_file = new TFile ((filename)) ;
  gDirectory->cd(dir);

  TH1F * current_hist = (TH1F*)gDirectory -> Get(histname);
  if(current_hist){
    current_hist -> SetDirectory(0);
  }
  TH1::AddDirectory(kFALSE);

  mapfunc[histname + Cycle_name + samplename] = current_hist;

  current_file -> Close();
  delete current_file;

}

void make_histogram(TString nameofhistogram, TString current_histname, int N_bin, double binx[], TString mass_point){

  TString title_y;
  if(nameofhistogram.Contains("N")) title_y = "Number";
  else title_y = "Events/bin";
  
  cout << "N_bin : " << N_bin << endl;
  
  TString canvas = nameofhistogram;
  TString hstack = nameofhistogram;
  TString legend = nameofhistogram;
  TString func = nameofhistogram;
  TString clone = nameofhistogram;
  TString line = nameofhistogram;
  canvas.Insert(0, "c_");
  hstack.Insert(0, "hs_");
  legend.Insert(0, "legend_");
  func.Insert(0, "ratio_");
  clone.Insert(0, "h3_");
  line.Insert(0, "l_");

  mapcanvas[canvas] = new TCanvas(canvas,"",800,800);
  canvas_margin(mapcanvas[canvas]);
  gStyle -> SetOptStat(1111);
  //mapcanvas[canvas] -> SetLogy();
  
  TH1F *official_signal = mapfunc[nameofhistogram + Cycle_name + map_sample_names["Official_" + mass_point].at(0)];
  TH1F *private_signal = mapfunc[nameofhistogram + Cycle_name + map_sample_names["Private_" + mass_point].at(0)];
  official_signal -> Rebin(200);
  private_signal -> Rebin(200);
 
  official_signal -> SetLineColor(kBlue);
  private_signal -> SetLineColor(kRed);
  maplegend[legend] = new TLegend(0.20, 0.80, 0.46, 0.92);
  maplegend[legend] -> AddEntry(official_signal, "official", "l");
  maplegend[legend] -> AddEntry(private_signal, "private", "l");
  
  TH1F * pad1_template = (TH1F*)official_signal -> Clone("template");
  gStyle->SetOptTitle(0);
  int N_bin_template = pad1_template -> GetNbinsX();
  for(int i = 1; i < N_bin_template; i++){
    pad1_template -> SetBinContent(i, 0.);
    pad1_template -> SetBinError(i, 0.);
  }

  double max_official = 0.;
  for(int i = 1; i < N_bin_template; i++){
    double current_content = official_signal -> GetBinContent(i);
    if(current_content > max_official) max_official = current_content;
  }



  pad1_template -> Draw("hist");
  pad1_template -> GetYaxis()->SetLabelSize(0.05);;
  pad1_template -> GetYaxis()->SetTitleSize(0.07);;
  pad1_template -> GetYaxis()->SetTitleOffset(1.02);;
  pad1_template -> GetXaxis()->SetLabelSize(0.05);
  pad1_template -> GetXaxis() -> SetTitle(nameofhistogram);
  pad1_template -> GetYaxis() -> SetRangeUser(0., max_official * 1.8);
  //pad1_template -> GetYaxis() -> SetRangeUser(0.01, max_official * 1000.); // logy 
  pad1_template -> GetYaxis() -> SetTitle(title_y);
  pad1_template -> SetStats(0);
  pad1_template -> Draw("histsame");
  
  official_signal -> Draw("histsame");
  private_signal -> Draw("histsame");
  maplegend[legend] -> Draw("same");

  
  // -- Latex
  TLatex latex_CMSPriliminary, latex_Lumi, latex_masspoint;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_masspoint.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Simulation}}}");
  latex_Lumi.SetTextSize(0.035);
  TString total_lumi = "41.3";
  if(tag_year == 2016) total_lumi = "35.9";
  if(tag_year == 2018) total_lumi = "59.7";
  latex_Lumi.DrawLatex(0.7, 0.96, total_lumi + " fb^{-1} (13 TeV)");
  latex_masspoint.SetTextSize(0.035);
  TString latex_mass_point = "EE_" + mass_point + " of " + TString::Itoa(tag_year,10);
  latex_masspoint.DrawLatex(0.5, 0.90, latex_mass_point);
  
  TString pdfname;
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  pdfname = WORKING_DIR + "/plots/Signal/";
  pdfname.Append(mass_point + "_" + nameofhistogram);
  pdfname.Append("_" + TString::Itoa(tag_year,10) + ".pdf");

  mapcanvas[canvas] -> SaveAs(pdfname);
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

  TString legend_list[1] = {"Zp4000_N500"};

  for(int i_syst = 0; i_syst < 1; i_syst++){
    TString current_hist_syst = histname + "_" + systematics[i_syst];
    TString current_dir_syst = current_dir + "_" + systematics[i_syst];

    for(int i_legend_list = 0; i_legend_list < 1; i_legend_list++){
      openfile(map_sample_names["Official_" + legend_list[i_legend_list]].at(0), current_dir_syst, current_hist_syst);
      openfile(map_sample_names["Private_" + legend_list[i_legend_list]].at(0), current_dir_syst, current_hist_syst);
    }
  }
  
  for(int i_legend_list = 0; i_legend_list < 1; i_legend_list++){
    make_histogram(histname + "_central", histname, N_bin, current_bins, legend_list[i_legend_list]);
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



void Signal_private_VS_official(int year=2018){
  setTDRStyle();
  tag_year = year;
  Cycle_name = "SR_ZpNN";
  map_sample_names["Official_Zp4000_N500"] = {"EE_Zp4000_N500"};
  map_sample_names["Private_Zp4000_N500"] = {"HNPairToJJJJ_EE_ZP4000_N500_WR5000"};

  open_binning_file("binning_private_VS_official.txt");

}
