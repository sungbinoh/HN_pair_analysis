#include "mylib.h"

void open_files(TString histname, int year){
  
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  //TString signal_list_file_path=WORKING_DIR+"/script/Signal_list/MC_signal_" + TString::Itoa(tag_year,10) + ".txt";
  TString signal_list_file_path=WORKING_DIR+"/script/Signal_list/MC_signal_2018.txt";
  TString root_file_path = WORKING_DIR+"/output/LimitTemplate/" + TString::Itoa(year,10) + "/";
  char line[500];
  ifstream signal_file;
  signal_file.open(signal_list_file_path);
 
  //TString bkgs[4] = {"VV", "TT_powheg", "DYJets_MG_HT", "WJets_MG_HT"};
  //TString bkgs[4] = {"DYJets_MG_HT", "TT_powheg", "WJets_MG_HT", "VV"};
  TString bkgs[4] = {"DYJets_MG_HT", "TT_powheg", "fake", "VV"};

  double lumi_error = 1.;
  if(year == 2016) lumi_error = 1.025;
  if(year == 2017) lumi_error = 1.023;
  if(year == 2018) lumi_error = 1.025;
  
  if(signal_file.is_open()){

    int abort = 0;
    while(!signal_file.eof()){
      signal_file.getline(line, 500);
      
      cout << line << endl;
      TString this_line = line;
      TString channel = "";
      TString channel2 = "";
      if(histname.Contains("DiEle")){
	channel = "DiEle";
	channel2 = "EE";
      }
      if(histname.Contains("DiMu")){
	channel= "DiMu";
	channel2 = "MuMu";
      }
      if(!this_line.Contains(channel2)) continue;
      
      TString current_input = root_file_path + "shape_" + histname + "_" + this_line + ".root";
      //cout << "current_input : " << current_input << endl;
      TFile *current_input_file = new TFile ((current_input)) ;
      //TH1F *obs_hist = (TH1F *)current_input_file->Get(histname);
      //double obs = obs_hist -> Integral();
      bool bkg_bool[4] = {false, false, false, false};
      //double signal_rate = 0.;
      
      double denom_bkgs[4] = {1., 1., 1., 1.};
      double denom_signal = 1.;
      int N_bkg = 0;
      for(int i= 0; i < 4; i++){
	TH1F *current_hist = (TH1F *)current_input_file->Get(histname + "_" + bkgs[i]);
	if(current_hist){
	  //cout << "current_hist : " << histname + "_" + bkgs[i] << endl;
	  bkg_bool[i] = true;
	  N_bkg++;
	  denom_bkgs[i] = current_hist -> Integral();
	}
      }
      TH1F *signal_hist = (TH1F *)current_input_file->Get(histname + "_" + this_line);
      denom_signal = signal_hist -> Integral();

      for(int i = 0; i < N_bkg; i++){
	for(int j = 0; j < N_syst_comparison; j++){
	  //cout << systematics_comparison[j] << endl;
	  TString str_current_hist_Up = histname + "_" + bkgs[i] + "_" + systematics_comparison[j] + "Up";
	  TString str_current_hist_Down = histname + "_" + bkgs[i] + "_" + systematics_comparison[j] + "Down";
	  //if(systematics_comparison[j].Contains("Jets") || systematics_comparison[j].Contains("SD") || systematics_comparison[j].Contains("ZPtRw")){
	  if(systematics_comparison[j].Contains("TriggerSF")){
	    str_current_hist_Up = histname + "_" + bkgs[i] + "_" + systematics_comparison[j] + TString::Itoa(year,10) + "Up";
	    str_current_hist_Down = histname + "_" + bkgs[i] + "_" + systematics_comparison[j] + TString::Itoa(year,10) + "Down";
	    //str_current_hist_Up = histname + "_" + bkgs[i] + "_" + systematics_comparison[j] + "Up";
	    //str_current_hist_Down = histname + "_" + bkgs[i] + "_" + systematics_comparison[j] + "Down";
	  }
	  
	  cout << "[open_files] str_current_hist_Up : " << str_current_hist_Up << endl;
	  TH1F *current_hist_Up = (TH1F *)current_input_file->Get(str_current_hist_Up);
	  TH1F *current_hist_Down = (TH1F *)current_input_file->Get(str_current_hist_Down);
	  double current_Up = current_hist_Up -> Integral();
	  double current_Down = current_hist_Down -> Integral();
	  current_Up = (current_Up - denom_bkgs[i]) / denom_bkgs[i];
	  current_Down = (current_Down - denom_bkgs[i]) / denom_bkgs[i];
	  double current_syst = pow(current_Up, 2) + pow(current_Down, 2);
	  current_syst = pow(current_syst, 0.5);
	  current_syst = current_syst * 100.;
	  
	  double new_min, new_max;
	  if(map_syst_table[bkgs[i] + systematics_comparison[j] + channel2].at(0) > current_syst) new_min = current_syst;
	  else new_min = map_syst_table[bkgs[i] + systematics_comparison[j] + channel2 + TString::Itoa(year,10) ].at(0);
	  
	  if(map_syst_table[bkgs[i] + systematics_comparison[j] + channel2].at(1) < current_syst) new_max = current_syst;
	  else new_max = map_syst_table[bkgs[i] + systematics_comparison[j] + channel2 + TString::Itoa(year,10)].at(1);

	  map_syst_table[bkgs[i] + systematics_comparison[j] + channel2 + TString::Itoa(year,10)] = {new_min, new_max};
	}
      }


      for(int j = 0; j < N_syst_comparison; j++){
	TString str_current_hist_Up = histname + "_" + this_line + "_" + systematics_comparison[j] + "Up";
	TString str_current_hist_Down = histname + "_" + this_line + "_" + systematics_comparison[j] + "Down";
	if(systematics_comparison[j].Contains("TriggerSF")){
	  str_current_hist_Up = histname + "_" + this_line + "_" + systematics_comparison[j] + TString::Itoa(year,10) + "Up";
	  str_current_hist_Down = histname + "_" + this_line + "_" + systematics_comparison[j] + TString::Itoa(year,10) + "Down";
	}
	
	TH1F *current_hist_Up = (TH1F *)current_input_file->Get(str_current_hist_Up);
	TH1F *current_hist_Down = (TH1F *)current_input_file->Get(str_current_hist_Down);
	double current_Up = current_hist_Up -> Integral();
	double current_Down = current_hist_Down -> Integral();
	current_Up = (current_Up - denom_signal) / denom_signal;
	current_Down = (current_Down - denom_signal) / denom_signal;
	double current_syst = pow(current_Up, 2) + pow(current_Down, 2);
	current_syst = pow(current_syst, 0.5);
	current_syst = current_syst * 100.;

	double new_min, new_max;
	if(map_syst_table[systematics_comparison[j] + channel2 + TString::Itoa(year,10)].at(0) > current_syst) new_min = current_syst;
	else new_min = map_syst_table[systematics_comparison[j] + channel2 + TString::Itoa(year,10)].at(0);

	if(map_syst_table[systematics_comparison[j] + channel2 + TString::Itoa(year,10)].at(1) < current_syst) new_max = current_syst;
	else new_max = map_syst_table[systematics_comparison[j] + channel2 + TString::Itoa(year,10)].at(1);

	map_syst_table[systematics_comparison[j] + channel2 + TString::Itoa(year,10)] = {new_min, new_max};
      }

      TString PDF_uncertainty[2] = {"PDF", "Scale"};
      for(int j = 0; j < 2; j++){

        TString str_current_hist_Up = histname + "_" + this_line + "_" + PDF_uncertainty[j] + "Up";
        TString str_current_hist_Down = histname + "_" + this_line + "_" + PDF_uncertainty[j] + "Down";

        TH1F *current_hist_Up = (TH1F *)current_input_file->Get(str_current_hist_Up);
        TH1F *current_hist_Down = (TH1F *)current_input_file->Get(str_current_hist_Down);
        double current_Up = current_hist_Up -> Integral();
        double current_Down = current_hist_Down -> Integral();
        current_Up = (current_Up - denom_signal) / denom_signal;
        current_Down = (current_Down - denom_signal) / denom_signal;
	double current_syst = pow(current_Up, 2) + pow(current_Down, 2);
        current_syst = pow(current_syst, 0.5);
        current_syst = current_syst * 100.;

        double new_min, new_max;
        if(map_syst_table[PDF_uncertainty[j] + channel2 + TString::Itoa(year,10)].at(0) > current_syst) new_min = current_syst;
        else new_min = map_syst_table[PDF_uncertainty[j] + channel2 + TString::Itoa(year,10)].at(0);

        if(map_syst_table[PDF_uncertainty[j] + channel2 + TString::Itoa(year,10)].at(1) < current_syst) new_max = current_syst;
        else new_max = map_syst_table[PDF_uncertainty[j] + channel2 + TString::Itoa(year,10)].at(1);

        map_syst_table[PDF_uncertainty[j] + channel2 + TString::Itoa(year,10)] = {new_min, new_max};
      }


      
      abort++;
      if(abort>5) break;
    }
  }
}

void open_binning_file(TString filename, int year){

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
      
      cout << "[[open_binning_file]] current_histname : " << current_histname << endl;

      open_files(current_histname + "_DYreweight", year);
    }
  }
  data_file.close();
}

void MakeSystematicTable_FullRun2(int ahahah){

  //TString bkgs[4] = {"VV", "TT_powheg", "DYJets_MG_HT", "WJets_MG_HT"};
  //TString bkgs[4] = {"DYJets_MG_HT", "TT_powheg", "WJets_MG_HT", "VV"};
  TString bkgs[4] = {"DYJets_MG_HT", "TT_powheg", "fake", "VV"};

  for(int i = 0; i < 4; i++){
    for(int j = 0; j < N_syst_comparison; j++){
      map_syst_table[bkgs[i] + systematics_comparison[j] + "EE2016"]   = {100., 0.};
      map_syst_table[bkgs[i] + systematics_comparison[j] + "MuMu2016"] = {100., 0.};
      map_syst_table[systematics_comparison[j] + "EE2016"] = {100., 0.};
      map_syst_table[systematics_comparison[j] + "MuMu2016"] = {100., 0.};

      map_syst_table[bkgs[i] + systematics_comparison[j] + "EE2017"]   = {100., 0.};
      map_syst_table[bkgs[i] + systematics_comparison[j] + "MuMu2017"] = {100., 0.};
      map_syst_table[systematics_comparison[j] + "EE2017"] = {100., 0.};
      map_syst_table[systematics_comparison[j] + "MuMu2017"] = {100., 0.};

      map_syst_table[bkgs[i] + systematics_comparison[j] + "EE2018"]   = {100., 0.};
      map_syst_table[bkgs[i] + systematics_comparison[j] + "MuMu2018"] = {100., 0.};
      map_syst_table[systematics_comparison[j] + "EE2018"] = {100., 0.};
      map_syst_table[systematics_comparison[j] + "MuMu2018"] = {100., 0.};

      map_syst_table[bkgs[i] + systematics_comparison[j] + "EE"]   = {100., 0.};
      map_syst_table[bkgs[i] + systematics_comparison[j] + "MuMu"] = {100., 0.};
      map_syst_table[systematics_comparison[j] + "EE"] = {100., 0.};
      map_syst_table[systematics_comparison[j] + "MuMu"] = {100., 0.};
    }
  }

  TString PDF_uncertainty[2] = {"PDF", "Scale"};
  for(int j = 0; j < 2; j++){
    map_syst_table[PDF_uncertainty[j] + "EE2016"] = {100., 0.};
    map_syst_table[PDF_uncertainty[j] + "MuMu2016"] = {100., 0.};

    map_syst_table[PDF_uncertainty[j] + "EE2017"] = {100., 0.};
    map_syst_table[PDF_uncertainty[j] + "MuMu2017"] = {100., 0.};

    map_syst_table[PDF_uncertainty[j] + "EE2018"] = {100., 0.};
    map_syst_table[PDF_uncertainty[j] + "MuMu2018"] = {100., 0.};

    map_syst_table[PDF_uncertainty[j] + "EE"] = {100., 0.};
    map_syst_table[PDF_uncertainty[j] + "MuMu"] = {100., 0.};

  }
   
  open_binning_file("binning_limit_merged.txt", 2016);
  open_binning_file("binning_limit_merged.txt", 2017);
  open_binning_file("binning_limit_merged.txt", 2018);

  for(int i = 0; i < 4; i++){
    for(int j = 0; j < N_syst_comparison; j++){
      map_syst_table[bkgs[i] + systematics_comparison[j] + "EE"] = { 
	std::min(std::min(map_syst_table[bkgs[i] + systematics_comparison[j] + "EE2016"].at(0), map_syst_table[bkgs[i] + systematics_comparison[j] + "EE2017"].at(0)), map_syst_table[bkgs[i] + systematics_comparison[j] + "EE2018"].at(0)),
	std::max(std::max(map_syst_table[bkgs[i] + systematics_comparison[j] + "EE2016"].at(1), map_syst_table[bkgs[i] + systematics_comparison[j] + "EE2017"].at(1)), map_syst_table[bkgs[i] + systematics_comparison[j] + "EE2018"].at(1))
      };
      
      map_syst_table[bkgs[i] + systematics_comparison[j] + "MuMu"] = {
	std::min(std::min(map_syst_table[bkgs[i] + systematics_comparison[j] + "MuMu2016"].at(0), map_syst_table[bkgs[i] + systematics_comparison[j] + "MuMu2017"].at(0)), map_syst_table[bkgs[i] + systematics_comparison[j] + "MuMu2018"].at(0)),
        std::max(std::max(map_syst_table[bkgs[i] + systematics_comparison[j] + "MuMu2016"].at(1), map_syst_table[bkgs[i] + systematics_comparison[j] + "MuMu2017"].at(1)), map_syst_table[bkgs[i] + systematics_comparison[j] + "MuMu2018"].at(1))
      };
      
      map_syst_table[systematics_comparison[j] + "EE"] = {
	std::min(std::min(map_syst_table[systematics_comparison[j] + "EE2016"].at(0), map_syst_table[systematics_comparison[j] + "EE2017"].at(0)), map_syst_table[systematics_comparison[j] + "EE2018"].at(0)),
        std::max(std::max(map_syst_table[systematics_comparison[j] + "EE2016"].at(1), map_syst_table[systematics_comparison[j] + "EE2017"].at(1)), map_syst_table[systematics_comparison[j] + "EE2018"].at(1))
      };
      
      map_syst_table[systematics_comparison[j] + "MuMu"] = {
	std::min(std::min(map_syst_table[systematics_comparison[j] + "MuMu2016"].at(0), map_syst_table[systematics_comparison[j] + "MuMu2017"].at(0)), map_syst_table[systematics_comparison[j] + "MuMu2018"].at(0)),
        std::max(std::max(map_syst_table[systematics_comparison[j] + "MuMu2016"].at(1), map_syst_table[systematics_comparison[j] + "MuMu2017"].at(1)), map_syst_table[systematics_comparison[j] + "MuMu2018"].at(1))
      };
    }
  }

  for(int j = 0; j < N_syst_comparison; j++){
    map_syst_table["bkgs" + systematics_comparison[j] + "EE"] = {
      std::min(std::min(std::min(map_syst_table[bkgs[0] + systematics_comparison[j] + "EE"].at(0),
				 map_syst_table[bkgs[1] + systematics_comparison[j] + "EE"].at(0)),
			map_syst_table[bkgs[2] + systematics_comparison[j] + "EE"].at(0) ),
	       map_syst_table[bkgs[3] + systematics_comparison[j] + "EE"].at(0)),
      std::max(std::max(std::max(map_syst_table[bkgs[0] + systematics_comparison[j] + "EE"].at(1),
                                 map_syst_table[bkgs[1] + systematics_comparison[j] + "EE"].at(1)),
			map_syst_table[bkgs[2] + systematics_comparison[j] + "EE"].at(1) ),
	       map_syst_table[bkgs[3] + systematics_comparison[j] + "EE"].at(1)),
    };
    
    map_syst_table["bkgs" + systematics_comparison[j] + "MuMu"] = {
      std::min(std::min(std::min(map_syst_table[bkgs[0] + systematics_comparison[j] + "MuMu"].at(0),
                                 map_syst_table[bkgs[1] + systematics_comparison[j] + "MuMu"].at(0)),
			map_syst_table[bkgs[2] + systematics_comparison[j] + "MuMu"].at(0) ),
	       map_syst_table[bkgs[3] + systematics_comparison[j] + "MuMu"].at(0)),
      std::max(std::max(std::max(map_syst_table[bkgs[0] + systematics_comparison[j] + "MuMu"].at(1),
				 map_syst_table[bkgs[1] + systematics_comparison[j] + "MuMu"].at(1)),
			map_syst_table[bkgs[2] + systematics_comparison[j] + "MuMu"].at(1) ),
               map_syst_table[bkgs[3] + systematics_comparison[j] + "MuMu"].at(1)),
    };
  }
  
  for(int j = 0; j < 2; j++){
    map_syst_table[PDF_uncertainty[j] + "EE"] = {
      std::min(std::min(map_syst_table[PDF_uncertainty[j] + "EE2016"].at(0), map_syst_table[PDF_uncertainty[j] + "EE2017"].at(0)), map_syst_table[PDF_uncertainty[j] + "EE2018"].at(0)),
      std::max(std::max(map_syst_table[PDF_uncertainty[j] + "EE2016"].at(1), map_syst_table[PDF_uncertainty[j] + "EE2017"].at(1)), map_syst_table[PDF_uncertainty[j] + "EE2018"].at(1))
    };
    
    map_syst_table[PDF_uncertainty[j] + "MuMu"] = {
      std::min(std::min(map_syst_table[PDF_uncertainty[j] + "MuMu2016"].at(0), map_syst_table[PDF_uncertainty[j] + "MuMu2017"].at(0)), map_syst_table[PDF_uncertainty[j] + "MuMu2018"].at(0)),
      std::max(std::max(map_syst_table[PDF_uncertainty[j] + "MuMu2016"].at(1), map_syst_table[PDF_uncertainty[j] + "MuMu2017"].at(1)), map_syst_table[PDF_uncertainty[j] + "MuMu2018"].at(1))
    };
  }

  TString systematics_table[N_syst] = {"ElectronScale",
				       "ElectronSmear",
				       "JetsRes",
				       "JetsScale",
				       //"SD Mass Res",
				       "SD Mass Scale",
				       "PU Reweight",
				       "MuonScale",
				       //"MuonSmear",
				       "L1 Prefire",
				       "MuonRecoSF",
				       "MuonIDSF",
				       "MuonISOSF",
				       "MuonTriggerSF",
				       "ElectronRecoSF",
				       "ElectronIDSF",
				       "ElectronTriggerSF",
				       "ZPtRw",
  };  
  
  TString syst_correlation[N_syst] = {"Correlated",
				      "Correlated",
				      "Correlated",
				      "Correlated",
				      //"Correlated",
				      "Correlated",
				      "Correlated",
				      "Correlated",
				      //"Not Done",
				      "Correlated",
				      "Correlated",
				      "Correlated",
				      "Correlated",
				      "Uncorrelated",
				      "Correlated",
				      "Correlated",
				      "Uncorrelated",
				      "Correlated",
  };

  TString syst_bkg_signal[N_syst] = {"All bkgd./Signal",
				     "All bkgd./Signal",
				     "All bkgd./Signal",
				     "All bkgd./Signal",
				     "All bkgd./Signal",
				     "All bkgd./Signal",
				     "All bkgd./Signal",
				     "All bkgd./Signal",
				     "All bkgd./Signal",
				     "All bkgd./Signal",
				     "All bkgd./Signal",
				     "All bkgd./Signal",
				     "All bkgd./Signal",
				     "All bkgd./Signal",
				     "All bkgd./Signal",
				     "DY+jets"
  };
  
  ofstream file_syst_table("./output/latex/syst_table.tex", ios::trunc);
  
  file_syst_table << "\\begin{table}"<< endl;
  file_syst_table << "\\centering" << endl;
  file_syst_table << "\\topcaption{The list of systematics and their impacts on total number of events at siganl region.}" << endl;
  file_syst_table << "\\label{tab:syst_ee}" << endl;
  file_syst_table << "\\resizebox{\\textwidth}{!}{" << endl;
  file_syst_table << "\\begin{tabular}{l l l l l l l}" << endl;
  file_syst_table << "\\hline\\hline" << endl;
  file_syst_table << "\\multirow{2}{*}{Source} & \\multirow{2}{*}{Bkgd./Signal process} & \\multirow{2}{*}{Year-to-year treatment} & $\\Pe\\Pe$ bkgd. & $\\Pe\\Pe$ signal & $\\mu\\mu$ bkgd. & $\\mu\\mu$ signal \\\\" << endl;
  file_syst_table << "&                                       &     &  (\\%)           & (\\%)            & (\\%)           & (\\%) \\\\" << endl;
  file_syst_table << "\\hline" << endl;

  file_syst_table << "Integrated luminosity & All bkgd./Signal & Uncorrelated & 2.3 -- 2.5 & 2.3 -- 2.5 & 2.3 -- 2.5 & 2.3 -- 2.5\\\\ " << endl;
  for(int j = 0; j < N_syst_comparison; j++){
    file_syst_table << systematics_table[j] << " & " << syst_bkg_signal[j] << " & " << syst_correlation[j];
    file_syst_table << fixed <<setprecision(1) << " & " ;
    if(map_syst_table["bkgs" + systematics_comparison[j] + "EE"].at(1) > 0.1){
      file_syst_table << map_syst_table["bkgs" + systematics_comparison[j] + "EE"].at(0) << " -- " << map_syst_table["bkgs" + systematics_comparison[j] + "EE"].at(1);
    }
    else file_syst_table << " $< 0.1$ ";
    file_syst_table << " & " ;
    
    if( systematics_table[j] == "ZPtRw"){
      file_syst_table << "\\NA";
    }
    else{
      if(map_syst_table[systematics_comparison[j] + "EE"].at(1) > 0.1){
	file_syst_table << map_syst_table[systematics_comparison[j] + "EE"].at(0) << " -- " << map_syst_table[systematics_comparison[j] + "EE"].at(1);
      }
      else file_syst_table << " $< 0.1$ ";
    }

    file_syst_table << " & " ;
    if(map_syst_table["bkgs" + systematics_comparison[j] + "MuMu"].at(1) > 0.1){
      file_syst_table << map_syst_table["bkgs" + systematics_comparison[j] + "MuMu"].at(0) << " -- " << map_syst_table["bkgs" + systematics_comparison[j] + "MuMu"].at(1);
    }    
    else file_syst_table << " $< 0.1$ ";
    file_syst_table << " & " ;

    if( systematics_table[j] == "ZPtRw"){
      file_syst_table << "\\NA";
    }
    else{
      if(map_syst_table[systematics_comparison[j] + "MuMu"].at(1) > 0.1){
	file_syst_table << map_syst_table[systematics_comparison[j] + "MuMu"].at(0) << " -- " << map_syst_table[systematics_comparison[j] + "MuMu"].at(1);
      }
      else file_syst_table << " $< 0.1$ ";
    }
    file_syst_table << "\\\\" << endl;
  }
  for(int j = 0; j < 2; j++){
    file_syst_table << PDF_uncertainty[j];
    file_syst_table << " & Signal & Correlated & \\NA & " << map_syst_table[PDF_uncertainty[j] + "EE"].at(0) << " -- " << map_syst_table[PDF_uncertainty[j] + "EE"].at(1) << " & \\NA & ";
    file_syst_table << map_syst_table[PDF_uncertainty[j] + "MuMu"].at(0) << " -- " << map_syst_table[PDF_uncertainty[j] + "MuMu"].at(1) << "\\\\" << endl;
  }
  file_syst_table << "Nonprompt norm. & Nonprompt & Correlated & 100 & \\NA & 100 & \\NA  \\\\" << endl;
  file_syst_table << "Rare SM norm. & Others & Correlated & 50 & \\NA & 50 & \\NA \\\\" << endl;
  file_syst_table << "\\hline\\hline" << endl;
  file_syst_table << "\\end{tabular}" << endl;
  file_syst_table << "}" << endl;
  file_syst_table << "\\end{table}" << endl;
  file_syst_table << "\n\n";
  

  /*
  for(int i = 0; i < 4; i++){
    for(int j = 0; j < N_syst_comparison; j++){
      cout << bkgs[i] << ", EE, " << systematics_comparison[j] << ", " << map_syst_table[bkgs[i] + systematics_comparison[j] + "EE"].at(0) << " - " << map_syst_table[bkgs[i] + systematics_comparison[j] + "EE"].at(1) << endl; 
    }
  }
  for(int j = 0; j < N_syst_comparison; j++){
    cout << "signal , EE, " << systematics_comparison[j] << ", " << map_syst_table[systematics_comparison[j] + "EE"].at(0) << " - " << map_syst_table[systematics_comparison[j] + "EE"].at(1) << endl;
  }
  for(int j = 0; j < 2; j++){
    cout << "signal , EE, " << PDF_uncertainty[j] << ", " << map_syst_table[PDF_uncertainty[j] + "EE"].at(0) << " - " << map_syst_table[PDF_uncertainty[j] + "EE"].at(1) << endl;
  }

  for(int i = 0; i < 4; i++){
    for(int j = 0; j < N_syst_comparison; j++){
      cout << bkgs[i] << ", MuMu, " << systematics_comparison[j] << ", " << map_syst_table[bkgs[i] + systematics_comparison[j] + "MuMu"].at(0) << " - " << map_syst_table[bkgs[i] + systematics_comparison[j] + "EE"].at(1) << endl;
    }
  }
  for(int j = 0; j < N_syst_comparison; j++){
    cout << "signal , MuMu, " << systematics_comparison[j] << ", " << map_syst_table[systematics_comparison[j] + "MuMu"].at(0) << " - " << map_syst_table[systematics_comparison[j] + "MuMu"].at(1) << endl;
  }
  for(int j = 0; j < 2; j++){
    cout << "signal , MuMu, " << PDF_uncertainty[j] << ", " << map_syst_table[PDF_uncertainty[j] + "MuMu"].at(0) << " - " << map_syst_table[PDF_uncertainty[j] + "MuMu"].at(1) << endl;
  }
  */

}
