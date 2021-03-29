#include "canvas_margin.h"
#include "mylib.h"

void openfile_background(TString cyclename, TString samplename, TString dir, TString histname){
  
  if(debug) cout << "[openfile_background] " << samplename << endl;

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

  if(debug) cout << "[openfile_DATA] " << samplename << endl;
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
    
void openfile_Signal(TString cyclename, TString samplename, TString dir, TString histname){

  if(debug) cout << "[openfile_Signal] " << samplename << endl;
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/rootfiles/" +TString::Itoa(tag_year,10) + "/Signal/";
  TString filename = root_file_path + cyclename + "_" + samplename + ".root";
  TFile *current_file = new TFile ((filename)) ;
  gDirectory->cd(dir);

  TH1F * current_hist = (TH1F*)gDirectory -> Get(histname);
  if(current_hist){
    current_hist -> SetDirectory(0);
  }
  TH1::AddDirectory(kFALSE);

  mapfunc[histname + cyclename + samplename] = current_hist;
  mapfunc[histname + cyclename + samplename] -> Scale(signal_scale);

  // -- Call Hists for PDF, alpha_s, and Scale errors
  if(histname.Contains("central")){
    
    // -- Call Hessian Hists
    for(int i_err = 1; i_err < 101; i_err++){
      TString current_histname = histname + "_Hessian_" + TString::Itoa(i_err, 10);
      //cout << "[openfile_Signal] current_histname : " <<  current_histname << endl;
      TH1F * current_hessian =  (TH1F*)gDirectory -> Get(current_histname);
      mapfunc[current_histname + cyclename + samplename] = current_hessian;
      mapfunc[current_histname + cyclename + samplename] -> Scale(signal_scale);
      if(current_hessian){
	current_hessian -> SetDirectory(0);
      }
      TH1::AddDirectory(kFALSE);
    }
    
    // -- Call Alpha_s Hists
    for(int i_err = 0; i_err < 2; i_err++){
      TString current_histname = histname + "_AlphaS_" + TString::Itoa(i_err, 10);
      //cout << "[openfile_Signal] current_histname : " <<  current_histname << endl;
      TH1F * current_alpha_s = (TH1F*)gDirectory -> Get(current_histname);
      mapfunc[current_histname + cyclename + samplename] = current_alpha_s;
      mapfunc[current_histname + cyclename + samplename] -> Scale(signal_scale);
      if(current_alpha_s){
	current_alpha_s -> SetDirectory(0);
      }
      TH1::AddDirectory(kFALSE);
    }
    
    // -- Call Scale Hists
    for(int i_err = 1; i_err < 9; i_err++){
      if(i_err == 5 || i_err == 7) continue;
      TString current_histname = histname + "_Scale_" + TString::Itoa(i_err, 10);
      cout << "[openfile_Signal] current_histname : " <<  current_histname << endl;
      TH1F * current_scale = (TH1F*)gDirectory -> Get(current_histname);
      mapfunc[current_histname + cyclename + samplename] = current_scale;
      mapfunc[current_histname + cyclename + samplename] -> Scale(signal_scale);
      if(current_scale){
	current_scale -> SetDirectory(0);
      }
      TH1::AddDirectory(kFALSE);
    }
  }// -- if hist name contains "central"
  
  current_file -> Close();
  delete current_file;
}

void make_histogram(TString nameofhistogram, TString current_histname, int N_bin, double binx[]){

  //cout << "[make_histogram] nameofhistogram : " << nameofhistogram << ", current_histname : " << current_histname << endl;
 
  if( GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets"].at(0))){

     
    /*
    Int_t idx0 = nameofhistogram.Index("__DYreweight");
    TString new_histname = "";
    for(int i_index = 0; i_index < idx0; i_index++){
      new_histname += nameofhistogram[i_index];
    }
    new_histname += "_central";
   
    if(GetHist(new_histname + Cycle_name + map_sample_names["DYJets"].at(0)) ){
      double integ_before = GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets"].at(0) ) -> Integral();
      double integ_after  = GetHist(new_histname + Cycle_name + map_sample_names["DYJets"].at(0)) -> Integral();
      double DY_norm = integ_before / integ_after;
      //cout << "integ_before : " << integ_before << endl;
      //cout << "integ_after : " << integ_after << endl;
      //cout << "DY_norm : " << DY_norm << endl;
      GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets"].at(0)) -> Scale(DY_norm);
    }
    */
    
  }
  
  if(debug){
    cout << "[[make_histogram]] checking binning arrary" << endl;
    for(int i = 0; i < N_bin; i++){
      cout << binx[i] << ", " << endl;
    }
  }
  
  TString current_data;
  if(nameofhistogram.Contains("EMu") || nameofhistogram.Contains("DiMu")) current_data = SingleMuon;
  else if(nameofhistogram.Contains("DiEle")) current_data = DoubleEG;
  else return;

  //gDirectory -> mkdir(current_data);
  //gDirectory -> cd(current_data);

  double N_DY = 0.;
  double N_MC = 0.;
  double N_Data = 0.;
  double DY_norm = 1.;
  if( GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets"].at(0))){
      
    if(debug) cout << "N_DY " << map_sample_names["DYJets"].at(0) << " : " << N_DY << endl;
    
    if(tag_year == 2016){
      if(nameofhistogram.Contains("DiMu")){
	if(nameofhistogram.Contains("0AK8")) DY_norm = 0.812047; 
	if(nameofhistogram.Contains("1AK8")) DY_norm = 0.708608;
        if(nameofhistogram.Contains("2AK8")) DY_norm = 0.603515;
      }
      if(nameofhistogram.Contains("DiEle")){
	if(nameofhistogram.Contains("0AK8")) DY_norm = 0.961725;
	if(nameofhistogram.Contains("1AK8")) DY_norm = 0.688508;
        if(nameofhistogram.Contains("2AK8")) DY_norm = 0.511236;
      }
    }
    if(tag_year == 2017){
      if(nameofhistogram.Contains("DiMu")){
	if(nameofhistogram.Contains("0AK8")) DY_norm = 1.09545;
	if(nameofhistogram.Contains("1AK8")) DY_norm = 1.05588;
        if(nameofhistogram.Contains("2AK8")) DY_norm = 0.74721;
      }
      if(nameofhistogram.Contains("DiEle")){
	if(nameofhistogram.Contains("0AK8")) DY_norm = 1.09268;
	if(nameofhistogram.Contains("1AK8")) DY_norm = 1.06172;
        if(nameofhistogram.Contains("2AK8")) DY_norm = 0.880778;
      }
    }if(tag_year == 2018){
      if(nameofhistogram.Contains("DiMu")){
	if(nameofhistogram.Contains("0AK8")) DY_norm = 1.12716;
	if(nameofhistogram.Contains("1AK8")) DY_norm = 0.946549;
        if(nameofhistogram.Contains("2AK8")) DY_norm = 0.997479;
      }
      if(nameofhistogram.Contains("DiEle")){
	if(nameofhistogram.Contains("0AK8")) DY_norm = 1.08663;
	if(nameofhistogram.Contains("1AK8")) DY_norm = 0.951804;
        if(nameofhistogram.Contains("2AK8")) DY_norm = 0.957299;
      }
    }
    
    GetHist(nameofhistogram + Cycle_name + map_sample_names["DYJets"].at(0)) -> Scale(DY_norm); 
    
    outfile << nameofhistogram << "'s DY_norm : " << DY_norm << endl;
    
    for(int i_syst = 1; i_syst < N_syst; i_syst++){
      if(GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["DYJets"].at(0))){
	GetHist(current_histname + "_" + systematics[i_syst] + Cycle_name + map_sample_names["DYJets"].at(0)) -> Scale(DY_norm);
      }
    }
  }
  
  bool blind = false;
  blind = (nameofhistogram.Contains("SR")) && (!nameofhistogram.Contains("EMu"));

  TString func = nameofhistogram;
  TString clone = nameofhistogram;
  func.Insert(0, "ratio_");
  clone.Insert(0, "h3_");
  
  TString legend_list[4] = {"Other", "WJets", "ttbar", "DYJets"};
  
  TString name_cycle = nameofhistogram + Cycle_name;
  if(debug) cout << "check1" << endl;
  
  // -- Add Overflow bin
  TString overflow = "overflow";
  Int_t nx_template = GetHist(nameofhistogram + Cycle_name + current_data) -> GetNbinsX();
  Double_t x1_template = GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinLowEdge(1);
  double bw = GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinWidth(nx_template);
  Double_t x2_template = GetHist(nameofhistogram + Cycle_name + current_data) -> GetBinLowEdge(nx_template) + bw;
  binx[N_bin] = binx[N_bin - 1] + 500.;
  //if(debug) cout << "[[make_histogram]] binx[N_bin] ; " << binx[N_bin]  << endl;
  //cout << "nx_template : " << nx_template << endl;
  //cout << "x1_template : " << x1_template << endl;
  //cout << "x2_template : " << x2_template << endl;
  mapfunc[func] = new TH1F("", "", nx_template, x1_template, x2_template);
  if(debug) cout << "func rebin rebinning" << endl;
  mapfunc[func + "rebin"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "rebin", binx);
  mapfunc[func + "rebin_stat_err"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "rebin", binx);
  mapfunc[func + "blind_data"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "blind_Data", binx);
  // -- For each Bkg Legend : Add overflow, Rebin, Syst error bar, Add on legend, and set fill & line colors
  for(int i_legend = 0; i_legend < 4; i_legend++){
    if(debug) cout << legend_list[i_legend] << endl;
    TString current_sample = map_sample_names[legend_list[i_legend]].at(0);
    if(mapfunc[nameofhistogram + Cycle_name + current_sample]){
      
      //Rebin_with_overflow_limit(nameofhistogram + Cycle_name + current_sample, N_bin, binx);
      Rebin_with_overflow(nameofhistogram + Cycle_name + current_sample, N_bin, binx);
      //cout << "GetHist  integral " << current_sample << " : " << GetHist(nameofhistogram + Cycle_name + current_sample)->Integral() << endl;
      //cout << "temp integral " << current_sample << " : " << htmp->Integral() << endl;
      
      
      if(debug) cout << current_sample + " : calling all syst hists" << endl;
      // -- Make rebinned hists for all systematic categories
      for(int i_syst = 0; i_syst < N_syst_comparison; i_syst++){
	//Write_syst_error_limit(current_histname, systematics_comparison[i_syst], Cycle_name, current_sample, N_bin, binx);
	Write_syst_error(current_histname, systematics_comparison[i_syst], Cycle_name, current_sample, N_bin, binx);
      }
      
      // -- Prepare blind data point (= bkg)
      mapfunc[func + "blind_data"]  -> Add(GetHist(nameofhistogram + Cycle_name + current_sample + "rebin"));
      //cout << nameofhistogram + Cycle_name + current_sample + "rebin" << endl;
      //cout << "data Nbin : " << mapfunc[func + "blind_data"] -> GetNbinsX() << ", add Nbin : " << GetHist(nameofhistogram + Cycle_name + current_sample + "rebin") -> GetNbinsX() << endl;
      mapfunc[nameofhistogram + Cycle_name + current_sample + "rebin"] -> SetName(current_histname + "_" + current_sample);// + "_central");
      //mapfunc[nameofhistogram + Cycle_name + current_sample + "rebin"] -> Write();
      
      //maplegend[legend] -> AddEntry(GetHist(nameofhistogram + Cycle_name + current_sample + "rebin"), legend_list[i_legend], "lf");
      mapfunc[func + "rebin"] -> Add(GetHist(nameofhistogram + Cycle_name + current_sample + "rebin"));
      mapfunc[func + "rebin_stat_err"] -> Add(GetHist(nameofhistogram + Cycle_name + current_sample + "rebin"));
    }
  }
  
  //Rebin_with_overflow_limit(nameofhistogram + Cycle_name + current_data, N_bin, binx);
  Rebin_with_overflow(nameofhistogram + Cycle_name + current_data, N_bin, binx);

  if(blind){
    mapfunc[nameofhistogram + Cycle_name + current_data + "rebin"] = (TH1F*)GetHist(func + "blind_data") -> Clone(clone + "blind");
  }
  //mapfunc[nameofhistogram + Cycle_name + current_data + "rebin"] -> SetName(current_histname);
  mapfunc[nameofhistogram + Cycle_name + current_data + "rebin"] -> SetName("data_obs");

  //mapfunc[nameofhistogram + Cycle_name + current_data + "rebin"] -> Write();
  
}

void Write_data_bkg(TString current_histname){
  
  TString nameofhistogram = current_histname + "_central";
  TString legend_list[4] = {"Other", "WJets", "ttbar", "DYJets"};
  for(int i_legend = 0; i_legend < 4; i_legend++){
    if(debug) cout << legend_list[i_legend] << endl;
    TString current_sample = map_sample_names[legend_list[i_legend]].at(0);
    if(mapfunc[nameofhistogram + Cycle_name + current_sample + "rebin"]){
      Remove_Negative_Bins(nameofhistogram + Cycle_name + current_sample + "rebin");
      mapfunc[nameofhistogram + Cycle_name + current_sample + "rebin"] -> Write();
    }
    for(int i_syst = 0; i_syst < N_syst_comparison; i_syst++){
      TString systematics = systematics_comparison[i_syst];
      TString cycle_and_sample = Cycle_name + current_sample;
      TString histname_Up   = current_histname + "_" + systematics + "Up"   + cycle_and_sample;
      TString histname_Down = current_histname + "_" + systematics + "Down" + cycle_and_sample;
      
      if(mapfunc[histname_Up + "rebin"]){
	Remove_Negative_Bins(histname_Up + "rebin");
	mapfunc[histname_Up + "rebin"] -> Write();
      }      
      if(mapfunc[histname_Down + "rebin"]){
	Remove_Negative_Bins(histname_Down + "rebin");
	mapfunc[histname_Down + "rebin"] -> Write();
      }    
    }
  }
  
  TH1F *fake_signal = (TH1F*)mapfunc[nameofhistogram + Cycle_name + map_sample_names["Other"].at(0) + "rebin"] -> Clone();
  fake_signal -> Scale(0.1);
  fake_signal -> SetName("signal");
  fake_signal -> Write();
  for(int i_syst = 0; i_syst < N_syst_comparison; i_syst++){
    TString current_sample = map_sample_names["Other"].at(0);
    TString systematics = systematics_comparison[i_syst];
    TString cycle_and_sample = Cycle_name + current_sample;
    TString histname_Up   = current_histname + "_" + systematics + "Up"   + cycle_and_sample;
    TString histname_Down = current_histname + "_" + systematics + "Down" + cycle_and_sample;
    
    TString name_up = "signal_" + systematics + "Up";
    TString name_down = "signal_" + systematics + "Down";
    // -- Uncorrelated syst
    if(systematics.Contains("TriggerSF")){
      name_up = "signal_" + systematics + TString::Itoa(tag_year,10) + "Up";
      name_down = "signal_" + systematics +  TString::Itoa(tag_year,10) + "Down";
    }
 
    TH1F *current_syst_up = (TH1F*)mapfunc[histname_Up + "rebin"] -> Clone();
    TH1F *current_syst_down = (TH1F*)mapfunc[histname_Down + "rebin"] -> Clone();
    current_syst_up -> Scale(0.1);
    current_syst_down -> Scale(0.1);
    current_syst_up -> SetName(name_up);
    current_syst_down -> SetName(name_down);
    current_syst_up -> Write();
    current_syst_down -> Write();
  }
  

  TString current_data;
  if(nameofhistogram.Contains("EMu") || nameofhistogram.Contains("DiMu")) current_data = SingleMuon;
  else if(nameofhistogram.Contains("DiEle")) current_data = DoubleEG;
  else return;
  mapfunc[nameofhistogram + Cycle_name + current_data + "rebin"] -> Write();

}

void make_histogram_signal(TString nameofhistogram, TString current_histname, TString mass_point, int N_bin, double binx[]){

  TString overflow = "overflow";
  /*
  Int_t nx_template = GetHist(nameofhistogram + Cycle_name + mass_point) -> GetNbinsX()+1;
  Double_t x1_template = GetHist(nameofhistogram + Cycle_name + mass_point) -> GetBinLowEdge(1);
  Double_t bw_template = (binx[N_bin - 1] - binx[0]) / 20.;
  Double_t x2_template = GetHist(nameofhistogram + Cycle_name + mass_point) -> GetBinLowEdge(nx_template)+bw_template;
  binx[N_bin] = binx[N_bin - 1] + bw_template;
  if(debug) cout << "[[make_histogram]] binx[N_bin] ; " << binx[N_bin]  << endl;
  mapfunc[func] = new TH1F("", "", nx_template, x1_template, x2_template);
  if(debug) cout << "func rebin rebinning" << endl;
  mapfunc[func + "rebin"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "rebin", binx);
  mapfunc[func + "rebin_stat_err"] = (TH1F*)mapfunc[func] -> Rebin(N_bin, func + "rebin", binx);
  */

  TString current_sample = mass_point;
  if(mapfunc[nameofhistogram + Cycle_name_signal + current_sample]){
    
    //Rebin_with_overflow_limit(nameofhistogram + Cycle_name_signal + current_sample, N_bin, binx);
    Rebin_with_overflow(nameofhistogram + Cycle_name_signal + current_sample, N_bin, binx);
    
    if(debug) cout << current_sample + " : calling all syst hists" << endl;
    for(int i_syst = 0; i_syst < N_syst_comparison; i_syst++){
      //Write_syst_error_limit(current_histname, systematics_comparison[i_syst], Cycle_name_signal, current_sample, N_bin, binx);
      Write_syst_error(current_histname, systematics_comparison[i_syst], Cycle_name_signal, current_sample, N_bin, binx);
    }
    mapfunc[nameofhistogram + Cycle_name_signal + current_sample + "rebin"] -> SetName(current_histname + "_" + current_sample);// + "_central");
    mapfunc[nameofhistogram + Cycle_name_signal + current_sample + "rebin"] -> Write();
    

    for(int i_syst = 0; i_syst < N_syst_comparison; i_syst++){
      TString systematics = systematics_comparison[i_syst];
      TString cycle_and_sample = Cycle_name_signal + current_sample;
      TString histname_Up   = current_histname + "_" + systematics + "Up"   + cycle_and_sample;
      TString histname_Down = current_histname + "_" + systematics + "Down" + cycle_and_sample;

      if(debug) cout << "[make_histogram_signal] histname_Up : " << histname_Up <<endl;

      mapfunc[histname_Up + "rebin"] -> Write();
      mapfunc[histname_Down + "rebin"] -> Write();
    }

    Estimate_PDF_Error(current_histname, Cycle_name_signal, current_sample, N_bin, binx);
    mapfunc[current_histname + "_" + current_sample + "_PDFUp"] -> Write();
    mapfunc[current_histname + "_" + current_sample + "_PDFDown"] -> Write();

    Estimate_Scale_Error(current_histname, Cycle_name_signal, current_sample, N_bin, binx);
    mapfunc[current_histname + "_" + current_sample + "_ScaleUp"] -> Write();
    mapfunc[current_histname + "_" + current_sample + "_ScaleDown"] -> Write();

  }
}


void open_files(TString histname){
  
  maphist.clear();
  map_syst_array.clear();
  mapfunc.clear();
  mapcanvas.clear();
  map_asym_gr.clear();
  
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
    TString legend_list[4] = {"DYJets", "ttbar", "WJets", "Other"};
    for(int i_legend_list = 0; i_legend_list < 4; i_legend_list++){
      
      unsigned int vector_size = map_sample_names[legend_list[i_legend_list]].size();
      //cout << "legend_list[i_legend_list] : " << legend_list[i_legend_list] << endl;
      //cout << "vector_size : " << vector_size << endl;
      
      for(int i_vector = 0; i_vector < vector_size; i_vector++){
	openfile_background(Cycle_name, map_sample_names[legend_list[i_legend_list]].at(i_vector), current_dir_syst, current_hist_syst);
      }
    }

    //cout << "current_dir_syst : " << current_dir_syst << endl;
    //cout << "current_hist_syst : " << current_hist_syst << endl;
    openfile_DATA(Cycle_name, DoubleEG, current_dir_syst, current_hist_syst);
    openfile_DATA(Cycle_name, SingleMuon, current_dir_syst, current_hist_syst);
    
  }
  // -- FIXME
  //histname = histname + "_central";
  //cout << "histname : " << histname << ", current_histname : " << current_histname << endl;
  //TFile *MyFile = new TFile("/output/Bkg_VS_signal_" + histname + ".root","RECREATE");
  //TFile *MyFile = new TFile("Bkg_VS_signal_" + histname + ".root","RECREATE");
  make_histogram(histname + "_central", histname, N_bin, current_bins);
  
  
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");

  TString channel2 = "";
  if(histname.Contains("DiEle")){
    channel2 = "EE";
  }
  if(histname.Contains("DiMu")){
    channel2 = "MuMu";
  }
  
  TFile *current_shape = new TFile("shape_" + histname + "_fake_signal.root","RECREATE");
  
  Write_data_bkg(histname);
  current_shape -> Close();
    
  //MyFile -> Close();
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

void MakeLimitTemplate_CR_fake_signal(int year=2018){
  setTDRStyle();
  
  // == Set Input Sample List
  tag_year = year;
  
  outfile.open("output_" + TString::Itoa(tag_year,10) + ".txt");
  
  if(tag_year == 2018){
    //map_sample_names["DYJets"] = {"DYJets_MG"};
    map_sample_names["DYJets"] = {"DYJets_MG_HT"};
    map_sample_names["ttbar"] = {"TT_powheg"};
    //map_sample_names["WJets"] = {"WJets_MG_HT"};
    map_sample_names["WJets"] = {"fake"};
    map_sample_names["Other"] = {"VV"};
    map_sample_names["Muon"] = {"data_SingleMuon"};
    map_sample_names["EGamma"] = {"data_EGamma"};
    
    SingleMuon = "data_SingleMuon";
    DoubleEG = "data_EGamma";
  }

  if(tag_year == 2017){
    //map_sample_names["DYJets"] = {"DYJets_MG"};
    map_sample_names["DYJets"] = {"DYJets_MG_HT"};
    map_sample_names["ttbar"] = {"TT_powheg"};
    //map_sample_names["WJets"] = {"WJets_MG_HT"};
    map_sample_names["WJets"] = {"fake"};
    map_sample_names["Other"] = {"VV"};
    map_sample_names["Muon"] = {"data_SingleMuon"};
    map_sample_names["EGamma"] = {"data_DoubleEG"};

    SingleMuon = "data_SingleMuon";
    DoubleEG = "data_DoubleEG";
  }

  if(tag_year == 2016){
    //map_sample_names["DYJets"] = {"DYJets_MG"};
    map_sample_names["DYJets"] = {"DYJets_MG_HT"};
    map_sample_names["ttbar"] = {"TT_powheg"};
    //map_sample_names["WJets"] = {"WJets_MG_HT"};
    map_sample_names["WJets"] = {"fake"};
    map_sample_names["Other"] = {"VV"};
    map_sample_names["Muon"] = {"data_SingleMuon"};
    map_sample_names["EGamma"] = {"data_DoubleEG"};

    SingleMuon = "data_SingleMuon";
    DoubleEG = "data_DoubleEG";
  }
  
  

  Cycle_name = "HN_pair_all_SkimTree_LRSMHighPt";
  open_binning_file("binning_limit_CR.txt");

  outfile.close();
  
}
