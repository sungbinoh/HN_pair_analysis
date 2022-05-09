#include "canvas_margin.h"
#include "mylib.h"

using namespace RooStats;
//using namespace RooFit ;

void QuickLimit_ratio(int xxx=0){
  
  double pvalue_1sig = SignificanceToPValue(1.);
  double pvalue_2sig = SignificanceToPValue(2.);
  double pvalue_3sig = SignificanceToPValue(3.);

  setTDRStyle();

  gStyle->SetOptStat(0);
  gStyle->SetLineWidth(2);

  TH1::SetDefaultSumw2(true);
  TH1::AddDirectory(kFALSE);

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString ENV_PLOT_PATH = getenv("PLOT_PATH");
  ENV_PLOT_PATH = ENV_PLOT_PATH + "/limit/";
  TString dirname = "limit_result";

  TString channel = "MuMu";
  TString binning = "3AK8_CRs";
  if(xxx== 1) channel = "EE";

  TString base_filepath_default = WORKING_DIR+"/output/Result_3AK8_CRs_" + channel + "_Zpt_big_syst.txt";
  TString base_filepath_new_syst = WORKING_DIR+"/output/Result_3AK8_CRs_" + channel + ".txt";
  TString base_plotpath = ENV_PLOT_PATH+"/Significance/";

  cout << "Open " << base_filepath_default << endl;
  if( !gSystem->mkdir(base_plotpath, kTRUE) ){
    cout
      << "###################################################" << endl
      << "Directoy " << base_plotpath << " is created" << endl
      << "###################################################" << endl
      << endl;
  }
  
  // -- Get Maximum Significance
  double Maximum_significance = 0.;  int mZP_max_sig = 0; int mN_max_sig = 0;
  map<TString, double> map_significance;
  vector<int> ZPmasses = {400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400};
  vector<int> ZPmasses_quarter = {400, 800, 1200, 1600, 2000, 2400, 2800, 3200, 3600, 4000, 4400};
  const int n_ZPmass = ZPmasses.size();
  const int n_ZPmass_quarter = ZPmasses_quarter.size();
  double x[n_ZPmass], x_quarter[n_ZPmass_quarter], exp_limit_default_boosted[n_ZPmass], exp_limit_default_resolved[n_ZPmass], exp_limit_new_syst_boosted[n_ZPmass], exp_limit_new_syst_resolved[n_ZPmass];
  
  for(unsigned int i = 0; i < n_ZPmass_quarter; i++){
    int mZP = ZPmasses_quarter.at(i);
    x_quarter[i] = mZP;
    vector<int> hnmasses = GetHNMassRange(mZP,true);
    const int n_Nmass = hnmasses.size();
    for(unsigned int j = 0; j < hnmasses.size(); j++){
      int mN = hnmasses.at(j);
    
      string line_str;
      ifstream in_default(base_filepath_default);
      ifstream in_new_syst(base_filepath_new_syst);
      double current_exp_limit_default = 0.;
      double current_exp_limit_new_syst = 0.;
      while(getline(in_default,line_str)){
	std::istringstream is( line_str );
	int a, b;
	double c_95_down;
	double c_68_down;
	double c;
	double c_68_up;
	double c_95_up;
	double c_obs;
	
	is >> a;
	is >> b;
	is >> c_95_down;
	is >> c_68_down;
	is >> c;
	is >> c_68_up;
	is >> c_95_up;
	is >> c_obs;
	
	if(a==mZP && b==mN){
	  current_exp_limit_default = c;
	}
      }
      while(getline(in_new_syst,line_str)){
	std::istringstream is( line_str );
	int a, b;
	double c_95_down;
	double c_68_down;
	double c;
	double c_68_up;
	double c_95_up;
	double c_obs;
	
	is >> a;
	is >> b;
	is >> c_95_down;
	is >> c_68_down;
	is >> c;
	is >> c_68_up;
	is >> c_95_up;
	is >> c_obs;
	
	if(a==mZP && b==mN){
	  current_exp_limit_new_syst = c;
	}
      }
      
      // -- Save Exp Limit for mN = mZP/4
      if(mN == mZP/4){
	exp_limit_default_resolved[i] = current_exp_limit_default;
	exp_limit_new_syst_resolved[i] = current_exp_limit_new_syst;
      }
    }
  }

  for(unsigned int i = 0; i < n_ZPmass; i++){
    int mZP = ZPmasses.at(i);
    x[i] = mZP;
    
    vector<int> hnmasses = GetHNMassRange(mZP,true);
    const int n_Nmass = hnmasses.size();

    for(unsigned int j = 0; j < hnmasses.size(); j++){
      int mN = hnmasses.at(j);
      
      string line_str;
      ifstream in_default(base_filepath_default);
      ifstream in_new_syst(base_filepath_new_syst);
      double current_exp_limit_default = 0.;
      double current_exp_limit_new_syst = 0.;
      while(getline(in_default,line_str)){
	std::istringstream is( line_str );
	int a, b;
        double c_95_down;
        double c_68_down;
        double c;
        double c_68_up;
        double c_95_up;
	double c_obs;

	is >> a;
        is >> b;
	is >> c_95_down;
        is >> c_68_down;
	is >> c;
        is >> c_68_up;
	is >> c_95_up;
        is >> c_obs;

	if(a==mZP && b==mN){
	  current_exp_limit_default = c;
	  //cout << "current_exp_limit_new_syst : " << current_exp_limit_new_syst << endl;
	}
      }
      while(getline(in_new_syst,line_str)){
	std::istringstream is( line_str );
        int a, b;
        double c_95_down;
        double c_68_down;
        double c;
        double c_68_up;
        double c_95_up;
        double c_obs;

        is >> a;
        is >> b;
        is >> c_95_down;
        is >> c_68_down;
        is >> c;
	is >> c_68_up;
        is >> c_95_up;
        is >> c_obs;

	if(a==mZP && b==mN){
          current_exp_limit_new_syst = c;
        }
      }

      TString mass_point = "mZP" + TString::Itoa(mZP,10) + "_mN" + TString::Itoa(mN,10);
      
      // -- Save Exp Limit for mN = 100 GeV
      if(mN == 100){
	exp_limit_default_boosted[i] = current_exp_limit_default;
	exp_limit_new_syst_boosted[i] = current_exp_limit_new_syst;

      }
    }
  }


  // -- Estimate ratio (new syst)/(default syst)
  for(unsigned int i = 0; i < n_ZPmass; i++){
    exp_limit_new_syst_boosted[i] = exp_limit_new_syst_boosted[i] / exp_limit_default_boosted[i];
    exp_limit_new_syst_resolved[i] = exp_limit_new_syst_resolved[i] / exp_limit_default_resolved[i];
    cout << "exp_limit_new_syst_boosted[i] : " << exp_limit_new_syst_boosted[i] << endl;
 }

  TCanvas *c1 = new TCanvas("c1", "", 600, 600);
  canvas_margin(c1);
  c1->cd();
  TH1D *hist_dummy = new TH1D("hist_dummy", "", 5000, 0., 5000.);
  hist_axis(hist_dummy);
  hist_dummy->Draw("hist");
  hist_dummy->GetYaxis()->SetRangeUser(0,2);
  hist_dummy->GetYaxis()->SetTitle("#sigma (pb)");
  hist_dummy->GetYaxis()->SetTitleSize(0.04);
  hist_dummy->GetYaxis()->SetTitleOffset(1.7);
  hist_dummy->GetYaxis()->SetLabelSize(0.045);
  hist_dummy->GetXaxis()->SetLabelSize(0.045);
  hist_dummy->GetYaxis()->SetTitle("Expected Limit Ratio (Half/Stat)");
  
  hist_dummy->GetXaxis()->SetRangeUser(400., 4400.);
  hist_dummy->GetXaxis()->SetTitle("m_{Z'} (GeV)");
    
  TLine *one_sigma = new TLine(400., 1., 4400., 1.);
  one_sigma -> SetLineStyle(9);
  one_sigma -> Draw("same");
  /*
  TLine *two_sigma = new TLine(400., pvalue_2sig, 4400., pvalue_2sig);
  two_sigma -> SetLineStyle(9);
  two_sigma -> Draw("same");
  TLine *three_sigma = new TLine(400., pvalue_3sig, 4400., pvalue_3sig);
  three_sigma -> SetLineStyle(9);
  three_sigma -> Draw("same");
  */
  TGraph *significance = new TGraph(n_ZPmass,x,exp_limit_new_syst_boosted);
  significance -> Draw("lsame");

  TLatex one_sigma_latex;
  TLatex two_sigma_latex;
  TLatex three_sigma_latex;
  one_sigma_latex.SetTextSize(0.035);
  two_sigma_latex.SetTextSize(0.035);
  three_sigma_latex.SetTextSize(0.035);
  one_sigma_latex.SetNDC();
  two_sigma_latex.SetNDC();
  three_sigma_latex.SetNDC();
  //one_sigma_latex.DrawLatex(0.20, 0.75 ,"1#sigma");
  //two_sigma_latex.DrawLatex(0.20, 0.57,"2#sigma");
  //three_sigma_latex.DrawLatex(0.20, 0.33,"3#sigma");

  c1->SaveAs(base_plotpath+"/Limit_ratio_N100_" + channel + ".pdf");

  c1->GetListOfPrimitives()->Remove(significance); c1->Modified();

  TGraph *significance_resolved = new TGraph(n_ZPmass_quarter,x_quarter,exp_limit_new_syst_resolved);
  significance_resolved -> Draw("lsame");

  c1->SaveAs(base_plotpath+"/Limit_ratio_quater_" + channel + ".pdf");

}
