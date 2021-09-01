#include "canvas_margin.h"
#include "mylib.h"

map<TString, std::vector<double> > mean;
map<TString, std::vector<double> > std_dev;
map<TString, std::vector<double> > eff;
map<TString, std::vector<double> > eff_err;

void open_files(TString histname, int year){

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  //TString signal_list_file_path=WORKING_DIR+"/script/Signal_list/MC_signal_2018.txt";
  TString signal_list_file_path=WORKING_DIR+"/script/Signal_list/list_match.txt";
  //TString root_file_path = WORKING_DIR+"/rootfiles/" + TString::Itoa(year,10) + "/Signal/";
  TString root_file_path = "/data6/Users/suoh/SKFlatOutput/Run2Legacy_v4/SR_ZpNN/" + TString::Itoa(year,10) + "/";
  char line[500];
  ifstream signal_file;
  signal_file.open(signal_list_file_path);

  if(signal_file.is_open()){
    int abort = 0;
    while(!signal_file.eof()){
      signal_file.getline(line, 500);

      cout << line << endl;
      TString this_line = line;
      TString channel = "";
      TString channel2 = "";
  
      if(this_line.Contains("EE")){
        channel = "DiEle";
        channel2 = "EE";
      }
      if(this_line.Contains("MuMu")){
        channel= "DiMu";
        channel2 = "MuMu";
      }
      if(!this_line.Contains("HNPairToJJJJ")) continue;
      // -- Open signal root file
      TString current_input = root_file_path + "SR_ZpNN_" + line + ".root";
      TFile *current_input_file = new TFile ((current_input)) ;
      TString dir = "SR_" + channel + "_DYreweight_central";
      gDirectory->cd(dir);
            
      TH1F * current_hist = (TH1F*)gDirectory -> Get(histname + "_SR_" + channel + "_DYreweight_central");
      //cout << dir << endl;
      //cout << histname + "_SR_" + channel + "_DYreweight_central" << endl;
      if(current_hist){
      	cout << current_hist ->GetMean() << "\t" << current_hist -> GetStdDev() << endl;
	cout << current_hist -> GetBinContent(0) / current_hist -> Integral() << "\t" << 
	  current_hist ->GetBinContent(1) / current_hist-> Integral() << "\t" <<
	  current_hist ->GetBinContent(2) / current_hist-> Integral() << "\t" <<
	  current_hist ->GetBinContent(3) / current_hist-> Integral() << endl;
	mean[this_line + histname].push_back(current_hist ->GetMean());
	std_dev[this_line + histname].push_back(current_hist -> GetStdDev());
	eff[this_line + histname].push_back(current_hist ->GetBinContent(2) / current_hist-> Integral());
	eff_err[this_line + histname].push_back(current_hist ->GetBinError(2) / current_hist-> Integral());
      }
      else cout << "NULL Hist" << endl;
      current_input_file->Close();

      abort++;
      //if(abort > 5) break; 

    }// -- End while signal txt open (eof)
  }// -- End if signal txt open


}

void make_plots(TString histname){

  TCanvas *c1 = new TCanvas("c1", "", 750, 600);
  canvas_margin(c1);
  c1->cd();

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString signal_list_file_path=WORKING_DIR+"/script/Signal_list/list_match.txt";  
  char line[500];
  ifstream signal_file;
  signal_file.open(signal_list_file_path);
  
  TString mass_points[8] = {"DiMu_Zp4000_N1500", "DiEle_Zp4000_N1500",
			   "DiMu_Zp4000_N1000", "DiEle_Zp4000_N1000",
			   "DiMu_Zp4000_N500", "DiEle_Zp4000_N500",
			   "DiMu_Zp4000_N100", "DiEle_Zp4000_N100",
  };

  TH1D *template_h = new TH1D("", "", 8, 0., 8.);
  template_h -> SetTitle(histname);
  template_h -> SetStats(0);
  double range_bottom = 0.;
  double range_top = 1.3;
  if(histname.Contains("ratio")){
    range_top = 2.0;
  }
  template_h -> GetYaxis() -> SetRangeUser(range_bottom, range_top);
  template_h -> GetYaxis() -> SetLabelSize(0.03);
  for(int i = 0; i < 8; i++){
    template_h -> GetXaxis() -> SetBinLabel( i + 1, mass_points[i]);
  }
  template_h -> GetXaxis() -> SetLabelSize(0.02);
  template_h -> LabelsOption("v");
  template_h -> Draw();

  double x[8]     = {0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
  double x_err[8] = { 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.};
  double y[8] =     { 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.};
  double y_err[8] = { 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.};
  if(signal_file.is_open()){
    int i_content = 0;
    while(!signal_file.eof()){
      signal_file.getline(line, 500);

      cout << line << endl;
      TString this_line = line;
      if(!this_line.Contains("HNPairToJJJJ")) continue;  
      if(histname.Contains("ratio")){
	y[i_content] = mean[this_line + histname].at(0);
	y_err[i_content] = std_dev[this_line + histname].at(0);
      }
      else{
	y[i_content] = eff[this_line + histname].at(0);
	y_err[i_content] = eff_err[this_line + histname].at(0);
      }
      i_content ++;
    }// -- End while signal txt open (eof)
  }// -- End if signal txt open

  cout << "contents filled" << endl;
  TGraphErrors *gr= new TGraphErrors(8, x, y, x_err, y_err);
  gr -> SetMarkerColor(4);
  gr -> SetLineColor(4);
  gr -> Draw("sameP");

  for(int i = 0; i < 7; i++){
    TString current_string = "Line" + TString::Itoa(i, 10);
    double current_x = i + 1;
    mapline[current_string] = new TLine(current_x, range_bottom, current_x, range_top);
    mapline[current_string] -> Draw();
    mapline[current_string] -> SetLineStyle(1);
    mapline[current_string] -> SetLineColor(kBlack);
  }
  TLegend *lg = new TLegend(0.65, 0.85, 0.95, 0.95);
  lg->SetBorderSize(1);
  if(histname.Contains("ratio")) lg->AddEntry(gr, "Ratio_Gen./Reco.", "lp");
  else lg->AddEntry(gr, "Eff #Delta R(Gen., Reco.) < 0.8", "lp");
  lg -> Draw("same");

  TString ENV_PLOT_PATH = getenv("PLOT_PATH");
  ENV_PLOT_PATH =  ENV_PLOT_PATH + "/Signal/";
  cout << "ENV_PLOT_PATH : " <<ENV_PLOT_PATH << endl;
  c1 -> SaveAs(ENV_PLOT_PATH + "/" + histname + ".pdf");
  
}

void HN_matching(){

  TString plots[15] = {"N_matched_Ns_0AK8", "N_matched_Ns_1AK8", "N_matched_Ns_2AK8",
		     "Pt_ratio_Ns0_0AK8", "Pt_ratio_Ns0_1AK8", "Pt_ratio_Ns0_2AK8", 
		     "Pt_ratio_Ns1_0AK8", "Pt_ratio_Ns1_1AK8", "Pt_ratio_Ns1_2AK8",
		     "Mass_ratio_Ns0_0AK8", "Mass_ratio_Ns0_1AK8", "Mass_ratio_Ns0_2AK8",
                     "Mass_ratio_Ns1_0AK8", "Mass_ratio_Ns1_1AK8", "Mass_ratio_Ns1_2AK8",
  };

  for(int i_plots = 0; i_plots < 15; i_plots ++){
    cout << plots[i_plots] << endl;
    open_files(plots[i_plots], 2016);
    make_plots(plots[i_plots]);
  }

  
  //open_files("Lepton_0_Pt", 2016);

}
