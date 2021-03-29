#include "mylib.h"

void open_files(TString histname){
  
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/output/LimitTemplate/" + TString::Itoa(tag_year,10) + "/";
   
  //TString bkgs[4] = {"VV", "TT_powheg", "DYJets_MG_HT", "WJets_MG_HT"};
  TString bkgs[4] = {"VV", "TT_powheg", "DYJets_MG_HT", "fake"};
 
  double lumi_error = 1.;
  /*
  if(tag_year == 2016) lumi_error = 1.025;
  if(tag_year == 2017) lumi_error = 1.023;
  if(tag_year == 2018) lumi_error = 1.025;
  */
  // -- Correlated Run2 Lumi error = 1.8 https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiLUM
  if(tag_year == 2016) lumi_error = 1.018;
  if(tag_year == 2017) lumi_error = 1.018;
  if(tag_year == 2018) lumi_error = 1.018;

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
  
  //TString current_input = root_file_path + "shape_" + histname + "_fake_signal.root";
  TString current_input = root_file_path + "shape_" + histname + ".root";

  TFile *current_input_file = new TFile ((current_input)) ;
  //TH1F *obs_hist = (TH1F *)current_input_file->Get(histname);
  //double obs = obs_hist -> Integral();
  bool bkg_bool[4] = {false, false, false, false};
  //double signal_rate = 0.;
  
  int N_bkg = 0;
  for(int i= 0; i < 4; i++){
    TH1F *current_hist = (TH1F *)current_input_file->Get(histname + "_" + bkgs[i]);
    if(current_hist){
      double current_integ = current_hist -> Integral();
      if(current_integ > 0){
	bkg_bool[i] = true;
	N_bkg++;
      }
    }
  }
  //TH1F *current_signal_hist = (TH1F *)current_input_file->Get(histname + "_" + this_line + "_central");
  //cout << histname + "_" + this_line + "_central" << endl;
  //signal_rate = current_signal_hist -> Integral();

  
  //ofstream file_shape("./output/DataCards/shape_" + histname + "_" + TString::Itoa(tag_year,10) + "_fake_signal.txt", ios::trunc);
  ofstream file_shape("./output/DataCards/shape_" + histname + "_" + TString::Itoa(tag_year,10) + ".txt", ios::trunc);

  //file_shape << "imax 1" << endl;
  //file_shape << "jmax " << N_bkg << endl;
  file_shape << "imax *" << endl;
  file_shape << "jmax *" << endl;
  file_shape << "kmax *" << endl;
  file_shape << "---------------" << endl;
  //file_shape << "shapes * * /data6/Users/suoh/Limit/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/ZpNN/LimitTemplate/" + TString::Itoa(tag_year,10) + "/shape_"<< histname << "_fake_signal.root $PROCESS $PROCESS_$SYSTEMATIC" << endl;
  file_shape << "shapes * * /data6/Users/suoh/Limit/CMSSW_10_2_13/src/HiggsAnalysis/CombinedLimit/ZpNN/LimitTemplate/" + TString::Itoa(tag_year,10) + "/shape_"<< histname << ".root $PROCESS $PROCESS_$SYSTEMATIC" << endl;
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
  file_shape << "\tsignal";
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
  for(int i = 0; i < 4; i++){
    if(bkg_bool[i]){
      adding_bkgs = adding_bkgs + "\t1";
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
  file_shape << "lumi" << "\tlnN\t" << lumi_error; // -- lumi correlated
  for(int i = 0; i < 4; i++){
    if(bkg_bool[i]) file_shape << "\t" << lumi_error;
  }
  file_shape << "\n";
  //file_shape << "ttbar_xsec\tlnN" + ttbar_xsec << endl;
  file_shape << "fake_norm\tlnN" + fake_norm << endl;
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
  file_shape << "JetsRes\tshapeN2" + adding_bkgs << endl;
  file_shape << "JetsScale\tshapeN2" + adding_bkgs << endl;
  file_shape << "SD_JMS_\tshapeN2" + adding_bkgs << endl;
  file_shape << "PUReweight_\tshapeN2" + adding_bkgs << endl;
  file_shape << "Prefire_\tshapeN2" + adding_bkgs << endl;
  file_shape << "ZPtRw\tshapeN2" + adding_bkgs << endl;
  if(histname.Contains("EMu")) file_shape << "R_ttbar_" + NAK8 + "_" + TString::Itoa(tag_year,10) + "\trateParam\tbin1\t" << histname + "_" + bkgs[1] << "\t1" << endl; // -- ttbar in emu sideband CR
  else file_shape << "R_DY_" + NAK8 + "_" + TString::Itoa(tag_year,10) + "\trateParam\tbin1\t" << histname + "_" + bkgs[2] << "\t1" << endl; // -- DY in ee & mumu CR
  file_shape << "* autoMCStats 0 0 1" << endl;
  file_shape.close();
  current_input_file -> Close();

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

void MakeDataCards_CR_fake_signal(int year=2018){

  tag_year = year;
  
  open_binning_file("binning_limit_CR.txt");

}
