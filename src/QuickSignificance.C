#include "canvas_margin.h"
#include "mylib.h"

using namespace RooStats;
//using namespace RooFit ;

void QuickSignificance(int xxx=0){
  
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

  TString base_filepath = WORKING_DIR+"/output/Significance_3AK8_CRs_significance_" + channel + ".txt";
  TString base_plotpath = ENV_PLOT_PATH+"/Significance/";

  cout << "Open " << base_filepath << endl;
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
  const int n_ZPmass = ZPmasses.size();
  double x[n_ZPmass], y_pvalue[n_ZPmass];
  

  for(unsigned int i = 0; i < n_ZPmass; i++){
    int mZP = ZPmasses.at(i);
    x[i] = mZP;
    
    vector<int> hnmasses = GetHNMassRange(mZP,true);
    const int n_Nmass = hnmasses.size();

    for(unsigned int j = 0; j < hnmasses.size(); j++){
      int mN = hnmasses.at(j);
      
      string line_str;
      ifstream in(base_filepath);

      double current_sig = 0.;
      while(getline(in,line_str)){
	std::istringstream is( line_str );
	int a, b;
	double sig;

	is >> a;
        is >> b;
        is >> sig;

	if(a==mZP && b==mN){
	  current_sig = sig;
	}
      }
      TString mass_point = "mZP" + TString::Itoa(mZP,10) + "_mN" + TString::Itoa(mN,10);
      
      // --- Save Maximum Significance
      if(current_sig > Maximum_significance){
	Maximum_significance = current_sig;
	mZP_max_sig = mZP;
	mN_max_sig = mN;
      }
      // -- Save P-value for mN = 100 GeV
      if(mN == 100){
	double current_pvalue = SignificanceToPValue(current_sig);
	y_pvalue[i] = current_pvalue;
      }
    }
  }

  cout << "Maximum_significance : " << Maximum_significance << ", p-value : " << SignificanceToPValue(Maximum_significance) << " , (mZ', mN) = (" << mZP_max_sig << ", " << mN_max_sig << ")" << endl;
  TCanvas *c1 = new TCanvas("c1", "", 600, 600);
  canvas_margin(c1);
  c1->cd();
  TH1D *hist_dummy = new TH1D("hist_dummy", "", 5000, 0., 5000.);
  hist_axis(hist_dummy);
  hist_dummy->Draw("hist");
  hist_dummy->GetYaxis()->SetRangeUser(1E-4,1);
  hist_dummy->GetYaxis()->SetTitle("#sigma (pb)");
  hist_dummy->GetYaxis()->SetTitleSize(0.04);
  hist_dummy->GetYaxis()->SetTitleOffset(1.7);
  hist_dummy->GetYaxis()->SetLabelSize(0.045);
  hist_dummy->GetXaxis()->SetLabelSize(0.045);
  hist_dummy->GetYaxis()->SetTitle("Local Significance (p-value)");
  
  hist_dummy->GetXaxis()->SetRangeUser(400., 4400.);
  hist_dummy->GetXaxis()->SetTitle("m_{Z'} (GeV)");
  c1->SetLogy();
  
  /*
  for(int i = 0; i < n_ZPmass; i++){
    cout << "x[" << i << "], y[" << i << "] = " << x[i] << ", " << y_pvalue[i] << endl;
  }
  */
  TLine *one_sigma = new TLine(400., pvalue_1sig, 4400., pvalue_1sig);
  one_sigma -> SetLineStyle(9);
  one_sigma -> Draw("same");
  TLine *two_sigma = new TLine(400., pvalue_2sig, 4400., pvalue_2sig);
  two_sigma -> SetLineStyle(9);
  two_sigma -> Draw("same");
  TLine *three_sigma = new TLine(400., pvalue_3sig, 4400., pvalue_3sig);
  three_sigma -> SetLineStyle(9);
  three_sigma -> Draw("same");
  TGraph *significance = new TGraph(n_ZPmass,x,y_pvalue);
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
  one_sigma_latex.DrawLatex(0.20, 0.75 ,"1#sigma");
  two_sigma_latex.DrawLatex(0.20, 0.57,"2#sigma");
  three_sigma_latex.DrawLatex(0.20, 0.33,"3#sigma");

  c1->SaveAs(base_plotpath+"/Significance_N100_" + channel + ".pdf");
  
  // -- Look elsewhere effect
  int N_1 = 0;
  // -- i = 13 : 3000 GeV
  for(int i = 21; i < n_ZPmass; i++){
    cout << x[i] << ", " << y_pvalue[i] << endl;
    if(y_pvalue[i] > 0) N_1++;
  }
  double N_1_err = pow(N_1 + 0., 0.5);
  double minimum_pvalue = SignificanceToPValue(Maximum_significance);
  double global_pvalue = minimum_pvalue + (N_1  + 0.) * exp((-pow(Maximum_significance,2))/2.);
  double global_pvalue_err = N_1_err * exp((-pow(Maximum_significance,2))/2.);
  double global_significance = PValueToSignificance(global_pvalue);
  double global_significance_up = PValueToSignificance(global_pvalue - global_pvalue_err);
  double global_significance_down = PValueToSignificance(global_pvalue + global_pvalue_err);
  
  cout << "global p-value : " << global_pvalue << " +- " << global_pvalue_err << 
    " , global significance : " << global_significance << " + " << global_significance_up - global_significance << " - " << global_significance - global_significance_down << endl;

}
