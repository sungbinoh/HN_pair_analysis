#include "mylib.h"

void open_files(TString histname){
  
  bool exo_datacard = false;
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  //TString signal_list_file_path=WORKING_DIR+"/script/Signal_list/MC_signal_" + TString::Itoa(tag_year,10) + ".txt";
  TString signal_list_file_path=WORKING_DIR+"/script/Signal_list/MC_signal_2016.txt";
  TString root_file_path = WORKING_DIR+"/output/LimitTemplate/" + TString::Itoa(tag_year,10) + "/";
  char line[500];
  ifstream signal_file;
  signal_file.open(signal_list_file_path);
 
  //TString bkgs[4] = {"VV", "TT_powheg", "DYJets_MG_HT", "WJets_MG_HT"};
  TString bkgs[4] = {"VV", "TT_powheg", "DYJets_MG_HT", "fake"};

  double lumi_error = 1.;
  /*
  if(tag_year == 2016) lumi_error = 1.025;
  if(tag_year == 2017) lumi_error = 1.023;
  if(tag_year == 2018) lumi_error = 1.025;
  */  
  // -- Correlated Run2 Lumi error = 1.8 https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiLUM
  if(tag_year == 2016) lumi_error = 1.016;
  if(tag_year == 2017) lumi_error = 1.016;
  if(tag_year == 2018) lumi_error = 1.016;
  
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
      TFile *current_input_file = new TFile ((current_input)) ;
      //TH1F *obs_hist = (TH1F *)current_input_file->Get(histname);
      //double obs = obs_hist -> Integral();
      bool bkg_bool[4] = {false, false, false, false};
      //double signal_rate = 0.;
      
      int N_bkg = 0;
      for(int i= 0; i < 4; i++){
	TH1F *current_hist = (TH1F *)current_input_file->Get(histname + "_" + bkgs[i]);
	if(current_hist){
	  bkg_bool[i] = true;
	  N_bkg++;
	}
      }
      //TH1F *current_signal_hist = (TH1F *)current_input_file->Get(histname + "_" + this_line + "_central");
      //cout << histname + "_" + this_line + "_central" << endl;
      //signal_rate = current_signal_hist -> Integral();

      ofstream file_shape("./output/DataCards/shape_" + histname + "_" + this_line + "_" + TString::Itoa(tag_year,10) + ".txt", ios::trunc);
      file_shape << "imax 1" << endl;
      file_shape << "jmax " << N_bkg << endl;
      file_shape << "kmax *" << endl;
      file_shape << "---------------" << endl;
      if(exo_datacard) 
	file_shape << "shapes * * shapes/" + TString::Itoa(tag_year,10) + "/shape_" << histname << "_" << this_line << ".root $PROCESS $PROCESS_$SYSTEMATIC" << endl;
      else file_shape << "shapes * * /data6/Users/suoh/Limit/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/ZpNN/LimitTemplate/" + TString::Itoa(tag_year,10) + "/shape_"<< histname << "_" << this_line << ".root $PROCESS $PROCESS_$SYSTEMATIC" << endl;
      file_shape << "---------------" << endl;
      file_shape << "bin\tbin1" << endl;
      //file_shape << "observation " << obs << endl;
      file_shape << "observation\t-1" << endl;
      file_shape << "------------------------------" << endl;
      file_shape << "bin\tbin1";//\tbin1\tbin1\tbin1\tbin1" << endl;
      for(int i= 0; i < 4; i++){
	if(bkg_bool[i]) file_shape << "\tbin1";
      }
      file_shape << "\n";
      file_shape << "process";
      file_shape << "\t" << histname + "_" + this_line;;
      for(int i = 0; i < 4; i++){
	if(bkg_bool[i]) file_shape << "\t" << histname + "_" + bkgs[i];
      }
      file_shape << "\n";
      file_shape << "process\t0";//\t1\t2\t3\t4" << endl;
      int current_process = 1;
      for(int i = 0; i < 4; i++){
	if(bkg_bool[i]){
	  file_shape << "\t" << current_process;
	  current_process++;
	}
      }
      file_shape << "\n";
      file_shape << "rate\t-1";//\t-1\t-1\t-1\t-1" << endl;;
      for(int i = 0; i < 4; i++){
        if(bkg_bool[i]) file_shape << "\t-1";
      }
      file_shape << "\n";
      /*
      file_shape << "rate";
      for(int i= 0; i < 4; i++){
	file_shape << "\t" << rate[i];
      }
      file_shape << "\t" << signal_rate << endl;
      */
      TString adding_bkgs = "\t1";
      TString adding_bkgs_signal_only = "\t1";
      for(int i = 0; i < 4; i++){
        if(bkg_bool[i]){
	  adding_bkgs = adding_bkgs + "\t1";
	  adding_bkgs_signal_only = adding_bkgs_signal_only + "\t-";
	}	
      }
      
      TString ttbar_xsec = "\t1";
      for(int i = 0; i < 4; i++){
	if(bkg_bool[i]){
	  if(i == 1) ttbar_xsec = ttbar_xsec + "\t1.05";
	  else ttbar_xsec = ttbar_xsec + "\t1";
	}
      }

      TString fake_norm = "\t1";
      for(int i = 0; i < 4; i++){
	if(bkg_bool[i]){
	  if(i == 3) fake_norm = fake_norm + "\t2.0";
	  else fake_norm = fake_norm + "\t1";
	}
      }
      
      TString Minor_MC_norm = "\t1";
      for(int i = 0; i < 4; i++){
        if(bkg_bool[i]){
          if(i == 0) Minor_MC_norm = Minor_MC_norm + "\t1.50";
          else Minor_MC_norm = Minor_MC_norm + "\t1";
        }
      }

      TString NAK8 = "";
      if(histname.Contains("0AK8")){
	NAK8 = "0AK8";
      }
      if(histname.Contains("1AK8")){
	NAK8 = "1AK8";
      }
      if(histname.Contains("2AK8")){
	NAK8 = "2AK8";
      }

      file_shape << "--------------------------------" << endl;
      //file_shape << "lumi" + TString::Itoa(tag_year,10) + "\tlnN\t" << lumi_error;// << "\t"<< lumi_error << "\t"<< lumi_error << "\t" << lumi_error << "\t" << lumi_error << endl; // -- lumi UnCo.
      file_shape << "lumi" << "\tlnN\t" << lumi_error;// << "\t"<< lumi_error << "\t"<< lumi_error << "\t" << lumi_error << "\t" << lumi_error << endl; // -- lumi correlated
      for(int i = 0; i < 4; i++){
	if(bkg_bool[i]) file_shape << "\t" << lumi_error;
      }
      file_shape << "\n";
      file_shape << "fake_norm" + TString::Itoa(tag_year,10) + "\tlnN" + fake_norm << endl;
      file_shape << "Minor_MC_norm\tlnN" + Minor_MC_norm << endl;
      file_shape << "ElectronScale\tshapeN2" + adding_bkgs << endl;
      file_shape << "ElectronSmear\tshapeN2" + adding_bkgs << endl;
      file_shape << "ElectronRecoSF\tshapeN2" + adding_bkgs << endl;
      file_shape << "ElectronIDSF\tshapeN2" + adding_bkgs << endl;
      file_shape << "ElectronTriggerSF" + TString::Itoa(tag_year,10) + "\tshapeN2" + adding_bkgs << endl;
      file_shape << "MuonScale\tshapeN2" + adding_bkgs << endl;
      //file_shape << "MuonSmear\tshapeN2" + adding_bkgs << endl;
      file_shape << "MuonRecoSF\tshapeN2" + adding_bkgs << endl;
      file_shape << "MuonIDSF\tshapeN2" + adding_bkgs << endl;
      file_shape << "MuonISOSF\tshapeN2" + adding_bkgs << endl;
      file_shape << "MuonTriggerSF" + TString::Itoa(tag_year,10) + "\tshapeN2" + adding_bkgs << endl;
      file_shape << "JetsRes" + TString::Itoa(tag_year,10) + "\tshapeN2" + adding_bkgs << endl;
      file_shape << "JetsScale\tshapeN2" + adding_bkgs << endl;
      //file_shape << "SD_JMR_\tshapeN2" + adding_bkgs << endl;
      file_shape << "SD_JMS_\tshapeN2" + adding_bkgs << endl;
      file_shape << "PUReweight_\tshapeN2" + adding_bkgs << endl;
      file_shape << "Prefire_" + TString::Itoa(tag_year,10) + "\tshapeN2" + adding_bkgs << endl;
      file_shape << "ZPtRw\tshapeN2" + adding_bkgs << endl;
      file_shape << "PDF\tshapeN2" + adding_bkgs_signal_only << endl;
      file_shape << "Scale\tshapeN2" + adding_bkgs_signal_only << endl;
      file_shape << "R_ttbar_" + NAK8 + "_" + TString::Itoa(tag_year,10) + "\trateParam\tbin1\t" << histname + "_" + bkgs[1] << "\t1" << endl; // -- ttbar
      file_shape << "R_DY_" + NAK8 + "_" + TString::Itoa(tag_year,10) + "\trateParam\tbin1\t" << histname + "_" + bkgs[2] << "\t1" << endl; // -- DY
      file_shape << "* autoMCStats 0 0 1" << endl;
      file_shape.close();
      current_input_file -> Close();

      abort++;
      //if(abort>5) break;
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
      
      cout << "[[open_binning_file]] current_histname : " << current_histname << endl;

      open_files(current_histname + "_DYreweight");
    }
  }
  data_file.close();
}

void MakeDataCards(int year=2018){

  tag_year = year;
  
  open_binning_file("binning_limit_merged.txt");

}
