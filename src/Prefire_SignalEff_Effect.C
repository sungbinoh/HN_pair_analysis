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

  current_file -> Close();
  delete current_file;
}

void openfile_signal(TString samplename, TString channel, TString dir, TString histname){

  TString filename = "SR_ZpNN_HNPairToJJJJ_" + channel + "_" + samplename + "_WR5000.root";
  cout << "[[openfile_signal]] Open " << filename << endl;

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/rootfiles/" +TString::Itoa(tag_year,10) + "/Signal_Study/";
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

void make_histogram(TString nameofhistogram, TString current_histname){

  // -- Draw Signal
  Int_t signal_colour_array[] = {622, 632, 432, 600};
  
  for(unsigned int i_signal = 0; i_signal < map_sample_names["Signal"].size(); i_signal++){
    int N_bin = mapfunc[current_histname + "_central" + map_sample_names["Signal"].at(i_signal)] -> GetNbinsX();
    
    double PF_applied = 0., err_PF_applied = 0.;
    double PF_no = 0., err_PF_no = 0.;
    PF_applied = mapfunc[current_histname + "_central" + map_sample_names["Signal"].at(i_signal)] -> Integral();
    PF_no = mapfunc[current_histname + "_NoPFweightcentral" + map_sample_names["Signal"].at(i_signal)] -> Integral();
    for(int i = 1; i < N_bin + 1; i++){
      err_PF_applied = err_PF_applied + pow(mapfunc[current_histname + "_central" + map_sample_names["Signal"].at(i_signal)] -> GetBinError(i), 2);
      err_PF_no = err_PF_no + pow(mapfunc[current_histname + "_NoPFweightcentral" + map_sample_names["Signal"].at(i_signal)] -> GetBinError(i), 2);
    }
    
    err_PF_applied = pow(err_PF_applied, 0.5);
    err_PF_no = pow(err_PF_no, 0.5);
        
    double ratio = PF_applied / PF_no;
    double error = ratio * (err_PF_applied / PF_applied + err_PF_no / PF_no);
    cout << map_sample_names["Signal"].at(i_signal) << ", ratio : " << ratio << " pm " << error << endl;
    
  }
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


  int N_channels = 6;
  TString channels[] = {"DiEle",
                        "DiMu",
                        "EMu",
                        "DiEle_DYreweight",
                        "DiMu_DYreweight",
                        "EMu_DYreweight",
  };

  TString current_dir = "SR";
  TString current_channel = "empty";

  for(int i_channel = 0; i_channel < N_channels; i_channel++){
    if(histname.Contains(channels[i_channel])){
      current_channel = channels[i_channel];
    }
  }
  current_dir = current_dir + "_" + current_channel;

  if(current_dir.Contains("empty") || current_channel.Contains("empty")) return;
  
  TString current_signal_channel = "";
  TString current_hist_syst = histname;;

  if(current_hist_syst.Contains("DiMu")) current_signal_channel = "MuMu";
  if(current_hist_syst.Contains("DiEle")) current_signal_channel = "EE";
  
  for(unsigned int i_signal = 0; i_signal < map_sample_names["Signal"].size(); i_signal++){
    openfile_signal(map_sample_names["Signal"].at(i_signal), current_signal_channel, current_dir + "_central" , histname + "_central");
    openfile_signal(map_sample_names["Signal"].at(i_signal), current_signal_channel, current_dir + "_NoPFweightcentral" , histname + "_NoPFweightcentral");

  }
  make_histogram(histname + "_central", histname);
  
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



void Prefire_SignalEff_Effect(int year=2018){
  setTDRStyle();
  
  // == Set Input Sample List
  tag_year = year;
  
  Cycle_name = "SR_ZpNN";

  map_sample_names["Signal"] = {"ZP1000_N200", "ZP2000_N900", "ZP3000_N200", "ZP3000_N1000", "ZP4000_N300", "ZP4000_N1500"};
  
  open_binning_file("binning_limit_merged.txt");

}
