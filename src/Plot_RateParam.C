#include "canvas_margin.h"
#include "mylib.h"

void Plot_RateParam(){

  setTDRStyle();

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2(true);
  TH1::AddDirectory(kFALSE);

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString ENV_PLOT_PATH = getenv("PLOT_PATH");
  WORKING_DIR = WORKING_DIR + "/rootfiles/fitDiagnostics/";
  ENV_PLOT_PATH = ENV_PLOT_PATH + "/limit/Rate_param/";

   
  TString channel = "MuMu";
  TString binning = "";
  
  TString file_DiMu = WORKING_DIR + "test_DiMu.txt";
  TString file_DiEle = WORKING_DIR + "test_DiEle.txt";

  const int N_Rate_Param = 18;
  TString Rate_param_names[N_Rate_Param] = {};
  double Rate_param_DiMu[N_Rate_Param] = {};
  double Rate_param_DiEle[N_Rate_Param] = {};
  double Rate_param_DiMu_err[N_Rate_Param] = {};
  double Rate_param_DiEle_err[N_Rate_Param] = {};
  
  TString channels[2] = {"DiMu", "DiEle"};
  for(int i = 0; i < 2; i++){
    TString current_channel = channels[i];
    string elline;
    ifstream in(WORKING_DIR+"/test_"+ current_channel + ".txt");
    
    int j = 0;
    while(getline(in,elline)){
      
      TString current_line = elline;
      if(!current_line.Contains("R_DY") && !current_line.Contains("R_ttbar")) continue;
      //cout << "current_line : " << current_line << endl;
      std::istringstream is( elline );
      
      // 2000  900 0.0042160                                                                                                                                                                              
      TString Param_name, a, b, d, f, g;
      double c, e;
      
      is >> Param_name;
      is >> a;
      is >> b;
      is >> c; // -- b-only fit Param
      is >> d;
      is >> e; // -- b-only fit Param Err
      is >> f;
      /*
      cout << "Param_name : " << Param_name << endl;
      cout << "a : " << a << endl;
      cout << "b : " << b << endl;
      cout << "c : " << c << endl;
      cout << "d : " << d << endl;
      cout << "e : " << e << endl;
      cout << "f : " << f << endl;
      cout << "g : " << g << endl;
      */
      Rate_param_names[j] = Param_name;
      if(current_channel.Contains("DiMu")){
	Rate_param_DiMu[j] = c;
	Rate_param_DiMu_err[j] = e;
      }
      if(current_channel.Contains("DiEle")){
	Rate_param_DiEle[j] = c;
	Rate_param_DiEle_err[j] = e;
      }

      j ++;
    }
  }

  // -- Check arraies
  for(int i = 0; i < N_Rate_Param; i++){
    
    cout << Rate_param_names[i] << "\t" << Rate_param_DiMu[i] << " \pm " << Rate_param_DiMu_err[i] << "\t" << Rate_param_DiEle[i] << " \pm " << Rate_param_DiEle_err[i] << endl;
    
  }

  TCanvas *c1 = new TCanvas("c1", "", 600, 600);
  gStyle -> SetOptStat(1111);
  canvas_margin(c1);
  c1->cd();

  TH1D* template_h = new TH1D("", "", 18, 0., 18.);
  template_h -> SetTitle("Limit ee Asym. VS HybridNew");
  template_h -> SetStats(0);
  template_h -> GetYaxis() -> SetRangeUser(0.5, 1.5);
  template_h -> GetYaxis() -> SetLabelSize(0.03);

  for(int i = 0; i < N_Rate_Param ; i++){
    template_h -> GetXaxis() -> SetBinLabel( i + 1, Rate_param_names[i]);
  }
  template_h -> GetXaxis() -> SetLabelSize(0.02);
  template_h -> LabelsOption("v");
  template_h -> Draw();

  double x[N_Rate_Param]  = {0.4, 1.4, 2.4, 3.4, 4.4, 5.4, 6.4, 7.4, 8.4, 9.4, 10.4, 11.4, 12.4, 13.4, 14.4, 15.4, 16.4, 17.4};
  double x2[N_Rate_Param] = {0.6, 1.6, 2.6, 3.6, 4.6, 5.6, 6.6, 7.6, 8.6, 9.6, 10.6, 11.6, 12.6, 13.6, 14.6, 15.6, 16.6, 17.6};
  double x_err[N_Rate_Param] = {0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.};
  TGraphErrors *gr_Rate_param_DiMu = new TGraphErrors(18, x, Rate_param_DiMu, x_err, Rate_param_DiMu_err);
  TGraphErrors *gr_Rate_param_DiEle = new TGraphErrors(18, x2, Rate_param_DiEle, x_err, Rate_param_DiEle_err);
  gr_Rate_param_DiMu -> SetMarkerColor(4);
  gr_Rate_param_DiMu -> SetLineColor(4);
  gr_Rate_param_DiMu -> Draw("sameP");
  
  gr_Rate_param_DiEle -> SetMarkerColor(kRed);
  gr_Rate_param_DiEle -> SetLineColor(kRed);
  gr_Rate_param_DiEle -> Draw("sameP");
  
  // -- Draw vertical lines
  for(int i = 0; i < 5; i++){
    TString current_string = "Line" + TString::Itoa(i, 10);
    double current_x = 3 * (i + 1);
    mapline[current_string] = new TLine(current_x, 0.5, current_x, 1.5);
    mapline[current_string] -> Draw();
    mapline[current_string] -> SetLineStyle(1);
    mapline[current_string] -> SetLineColor(kBlack);
  }
  TLegend *lg = new TLegend(0.65, 0.85, 0.95, 0.95);
  lg->SetBorderSize(1);
  //lg->SetFillStyle(0);
  lg->AddEntry(gr_Rate_param_DiMu, "Rate_Param_DiMu", "lp");
  lg->AddEntry(gr_Rate_param_DiEle, "Rate_Param_DiEle", "lp");
  lg -> Draw("same");

  c1 -> SaveAs(ENV_PLOT_PATH + "/Rate_Param_comparison.pdf");
  


}
 





