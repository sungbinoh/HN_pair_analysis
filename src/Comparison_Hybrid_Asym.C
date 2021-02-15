#include "canvas_margin.h"
#include "mylib.h"

void Comparison_Hybrid_Asym(){

  TString ENV_PLOT_PATH = getenv("PLOT_PATH");
  ENV_PLOT_PATH = ENV_PLOT_PATH + "/limit/";
  TString dirname = "limit_result";
  TString base_plotpath = ENV_PLOT_PATH+"/3AK8_CRs/";

  TCanvas *c1 = new TCanvas("c1", "", 600, 600);
  gStyle -> SetOptStat(1111);
  canvas_margin(c1);
  c1->cd();

  double x[6] = {1., 2., 3., 4., 5., 6.};
  double x_asym[6] = {1.2, 2.2, 3.2, 4.2, 5.2, 6.2};
  double x_err[6] = {0., 0., 0., 0., 0., 0.};
  double asym_ee[6] = {6.875, 16.0625, 4.2031, 3.7969, 18.8125, 21.5};
  double asym_mumu[6] = {7.5938, 14.1875, 2.4141, 3.8281, 7.0, 22.375};
  double hybrid_ee[6] = {7.21048, 16.406, 4.16018, 4.03914, 15.625, 23.4211};
  double hybrid_mumu[6] = {7.80224, 15.4939, 2.43982, 3.95724, 7.55617, 25.022};
  double hybrid_ee_err[6] = {0.163571, 0.5, 0.066129, 0.125, 0.625, 1.07967};
  double hybrid_mumu_err[6] = {0.173548, 0.25, 0.0236424, 0.125, 0.275661, 1.04122};
  TString label[6] = {"(2000, 200)", "(2000, 700)", "(3200, 200)", "(3200, 600)", "(4400, 200)", "(4400, 1800)"};
  
  TH1D* template_h = new TH1D("", "", 7, 0., 7);
  template_h -> SetTitle("Limit ee Asym. VS HybridNew");
  template_h -> SetStats(0);
  template_h -> GetYaxis() -> SetRangeUser(0., 30.);
  for(int i = 0; i < 6;i++){
    template_h -> GetXaxis() -> SetBinLabel( i + 1, label[i]);
  }
  template_h -> Draw();
  TGraphErrors *Hybrid_ee = new TGraphErrors(6, x, hybrid_ee, x_err, hybrid_ee_err);
  TGraphErrors *Hybrid_mumu = new TGraphErrors(6, x, hybrid_mumu, x_err, hybrid_mumu_err);
  
  TGraph *Asym_ee = new TGraph(6, x_asym, asym_ee);
  TGraph *Asym_mumu = new TGraph(6, x_asym, asym_mumu);

  // -- For EE
  //Hybrid_ee -> SetTitle("Limit ee Asym. VS HybridNew");
  Hybrid_ee -> SetMarkerColor(4);
  Hybrid_ee -> SetMarkerStyle(21);
  Hybrid_ee -> Draw("sameP");
  Asym_ee -> SetMarkerStyle(21);
  Asym_ee -> SetMarkerColor(kRed);
  Asym_ee -> Draw("sameP");

  TLegend *lg = new TLegend(0.2, 0.8, 0.3, 0.95);
  lg->SetBorderSize(0);
  lg->SetFillStyle(0);
  lg->AddEntry(Hybrid_ee, "HybridNew", "lp");
  lg->AddEntry(Asym_ee, "Asymptotic", "lp");
  lg -> Draw("same");
  c1 -> SaveAs(base_plotpath + "/EE/Asym_vs_Hybridnew.pdf");

  // -- For MuMu
  template_h -> SetTitle("Limit mumu Asym. VS HybridNew");
  template_h ->Draw();
  Hybrid_mumu -> SetMarkerColor(4);
  Hybrid_mumu -> SetMarkerStyle(21);
  for(int i = 0; i < 6;i++){
    Hybrid_mumu -> GetXaxis()->SetBinLabel( i + 1, label[i]);
  }
  Hybrid_mumu -> Draw("sameP");
  Asym_mumu -> SetMarkerStyle(21);
  Asym_mumu -> SetMarkerColor(kRed);
  Asym_mumu -> Draw("sameP");

  TLegend *lg_mumu = new TLegend(0.2, 0.8, 0.3, 0.95);
  lg_mumu->SetBorderSize(0);
  lg_mumu->SetFillStyle(0);
  lg_mumu->AddEntry(Hybrid_mumu, "HybridNew", "lp");
  lg_mumu->AddEntry(Asym_mumu, "Asymptotic", "lp");
  lg_mumu -> Draw("same");
  c1 -> SaveAs(base_plotpath + "/MuMu/Asym_vs_Hybridnew.pdf");

}
