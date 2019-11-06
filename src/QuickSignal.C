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

      openfile_signal("ZP" + current_mZp + "_N" + current_mN, "MuMu");

      mN = mN + 200;
    }
  }
  


}









