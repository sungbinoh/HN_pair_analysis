#include "canvas_margin.h"
#include "mylib.h"

void openfile_signal(TString samplename, TString channel){
  
  TString filename = "HN_pair_all_HNPairToJJJJ_" + channel + "_" + samplename + "_WR5000.root";
  cout << "[[openfile_signal]] Open " << filename << endl;

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/rootfiles/" +TString::Itoa(tag_year,10) + "/Signal/";
  filename = root_file_path + filename;
  TFile *current_file = new TFile ((filename)) ;
  
  gDirectory->cd("SR_DiMu_DYreweight_central");
  maphist[samplename + channel] =  (TH1*)gDirectory -> Get("mZp_SR_DiMu_DYreweight_central");
  maphist[samplename + channel + "b_veto"] =  (TH1*)gDirectory -> Get("mZp_SR_b_veto_DiMu_DYreweight_central");
  
  double N_origianl = maphist[samplename + channel] -> Integral();
  double N_b_veto = maphist[samplename + channel + "b_veto"] -> Integral();
  double ratio = N_b_veto / N_origianl;
  cout << "ratio : " << ratio << endl;
  
  current_file -> Close();
  
}    

void openfile_signal_hist(TString samplename, TString channel, TString dir, TString histname){
  
  TString filename = "HN_pair_all_HNPairToJJJJ_" + channel + "_" + samplename + "_WR5000.root";
  cout << "[[openfile_signal]] Open " << filename << endl;

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/rootfiles/" +TString::Itoa(tag_year,10) + "/Signal/";
  filename = root_file_path + filename;
  TFile *current_file = new TFile ((filename)) ;

  gDirectory->cd(dir);

  TH1F * current_hist = (TH1F*)gDirectory -> Get(histname);
  if(current_hist){
    current_hist -> SetDirectory(0);
  }
  TH1::AddDirectory(kFALSE);
  
  cout << "[openfile_signal_hist] " << samplename +  channel + histname << endl;
  mapfunc[samplename +  channel + histname] = current_hist;
  
  current_file -> Close();
  delete current_file;
  
}

void plot_N_bjet(TString channel){
  
  TString channel_dir;
  if(channel.Contains("MuMu")) channel_dir = "DiMu";
  TString current_hist_name = "Nbjet_SR_" + channel_dir + "_DYreweight_central";
  mapcanvas[current_hist_name] = new TCanvas(current_hist_name,"",800,800);
  
  Int_t color_array[] = {632, 901, 800, 400, 416, 432, 600, 880, 883};
  Int_t line_array[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

  cout << "ZP400_N100" + channel + current_hist_name << endl;
  mapfunc["ZP400_N100" + channel + current_hist_name] -> Draw("hist");
  double integral = mapfunc["ZP400_N100" + channel + current_hist_name] -> Integral();
  mapfunc["ZP400_N100" + channel + current_hist_name] -> Scale(1./integral);
  mapfunc["ZP400_N100" + channel + current_hist_name] -> GetYaxis()->SetLabelSize(0.05);
  mapfunc["ZP400_N100" + channel + current_hist_name] -> GetYaxis()->SetTitleSize(0.07);
  mapfunc["ZP400_N100" + channel + current_hist_name] -> GetYaxis()->SetTitleOffset(1.02);
  //mapfunc["ZP400_N100" + channel + current_hist_name] -> GetXaxis()->SetLabelSize(0);
  int maximum_bin = mapfunc["ZP400_N100" + channel + current_hist_name] -> GetMaximumBin();
  double maximum = mapfunc["ZP400_N100" + channel + current_hist_name] -> GetBinContent(maximum_bin);
  mapfunc["ZP400_N100" + channel + current_hist_name] -> GetYaxis() -> SetRangeUser(0., maximum * 1.5);
  mapfunc["ZP400_N100" + channel + current_hist_name] -> GetXaxis() -> SetTitle("N(bjet)");
  mapfunc["ZP400_N100" + channel + current_hist_name] -> SetLineColor(632);
  mapfunc["ZP400_N100" + channel + current_hist_name] -> SetLineStyle(1);

  mapfunc["ZP400_N100" + channel + current_hist_name] -> GetYaxis() -> SetTitle("Arbitary Unit");
  mapfunc["ZP400_N100" + channel + current_hist_name] -> Draw("histsame");
  
  TString map_name = current_hist_name;
  maplegend[map_name] = new TLegend(0.70, 0.15, 0.95, 0.95);
  maplegend[map_name] -> AddEntry(mapfunc["ZP400_N100" + channel + current_hist_name], "ZP400_N100", "fl");
  int i_line = 0;
  for(int i_Zp = 0; i_Zp < 12; i_Zp ++){
    
    int mZp = 400 + 400 * i_Zp;
    int i_color = 0;
    if(mZp == 400 || mZp == 4800 || mZp == 800) continue;
    i_line ++;
    int mN = 100;
    
    while(mN < mZp / 2){
      TString current_mZp = TString::Itoa(mZp,10);
      TString current_mN = TString::Itoa(mN,10);

      TString channel ="MuMu";
      TString current_dir = "SR_" + channel + "_DYreweight_central";
      
      double current_integral = mapfunc["ZP" + current_mZp + "_N" + current_mN + channel + current_hist_name] -> Integral();
      mapfunc["ZP" + current_mZp + "_N" + current_mN + channel + current_hist_name] -> Scale(1./current_integral);
      mapfunc["ZP" + current_mZp + "_N" + current_mN + channel + current_hist_name] -> SetLineColor(color_array[i_color]);
      i_color++;
      mapfunc["ZP" + current_mZp + "_N" + current_mN + channel + current_hist_name] -> SetLineStyle(line_array[i_line]);

      mapfunc["ZP" + current_mZp + "_N" + current_mN + channel + current_hist_name] -> Draw("histsame");
      maplegend[map_name] -> AddEntry(mapfunc["ZP" + current_mZp + "_N" + current_mN + channel + current_hist_name], "ZP" + current_mZp + "_N" + current_mN, "fl");
      mN = mN + 200;
    }
  }

  maplegend[map_name] -> Draw("lcsame");

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString pdfname = WORKING_DIR + "/plots/Signal/";
  mapcanvas[current_hist_name] -> SaveAs(pdfname + "Nbjet_signal.pdf");
  
}

void QuickSignal(int year=2018){
  setTDRStyle();
  
  // == Set Input Sample List
  tag_year = year;

  maphist2D["b_veto_effect_MuMu"] = new TH2D("", "", 12, 400., 5200., 12, 100., 2500.);
  
  for(int i_Zp = 0; i_Zp < 12; i_Zp ++){
    int mZp = 400 + 400 * i_Zp;
    int mN = 100;
    while(mN < mZp / 2){
      TString current_mZp = TString::Itoa(mZp,10);
      TString current_mN = TString::Itoa(mN,10);
      
      //openfile_signal("ZP" + current_mZp + "_N" + current_mN, "MuMu");
      TString channel = "MuMu";
      TString channel_dir = "DiMu";

      TString current_dir = "SR_" + channel_dir + "_DYreweight_central";
	
      openfile_signal_hist("ZP" + current_mZp + "_N" + current_mN, channel, current_dir, "Nbjet_" + current_dir);
      mN = mN + 200;
    }
  }
  
  plot_N_bjet("MuMu");
  

}









