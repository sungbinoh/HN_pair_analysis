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
  //c1 -> SetLogy();

  double x[6] = {1., 2., 3., 4., 5., 6.};
  double x_asym[6] = {1.2, 2.2, 3.2, 4.2, 5.2, 6.2};
  double x_err[6] = {0., 0., 0., 0., 0., 0.};
  double asym_ee[6] = {1.4148, 2.8714, 0.8331, 0.8104, 1.8439, 0.8086 };
  double asym_mumu[6] = {0.8364, 1.1380, 0.1541, 0.2896, 0.1786, 0.3631};
  double asym_ee_err[6] = {0.0062318, 0.0126248, 0.0032709, 0.0032728, 0.0064874, 0.0018440};
  double asym_mumu_err[6] = {0.0031431, 0.0027777, 0.0005559, 0.0013467, 0.0007584, 0.0016108};
  double hybrid_ee[6] = {1.43112, 2.97378, 0.8239, 0.824599, 1.93009, 0.865366 };
  double hybrid_mumu[6] = {0.866629, 1.18734, 0.158411, 0.29194, 0.178933, 0.306487};
  double hybrid_ee_err[6] = {0.0134352, 0.0310753, 0.012052, 0.0114564, 0.0376278, 0.0148956 };
  double hybrid_mumu_err[6] = {0.0900444, 0.0220893, 0.00535674, 0.00536097, 0.00503208, 0.0304598};
  TString label[6] = {"(2000, 200)", "(2000, 700)", "(3200, 200)", "(3200, 600)", "(4000, 200)", "(4000, 1600)"};
  
  TH1D* template_h = new TH1D("", "", 7, 0., 7);
  template_h -> SetTitle("Limit ee Asym. VS HybridNew");
  template_h -> SetStats(0);
  //template_h -> GetYaxis() -> SetRangeUser(0.01, 100.); // logy
  template_h -> GetYaxis() -> SetRangeUser(0.1, 4.);
  for(int i = 0; i < 6;i++){
    template_h -> GetXaxis() -> SetBinLabel( i + 1, label[i]);
  }
  template_h -> Draw();
  TGraphErrors *Hybrid_ee = new TGraphErrors(6, x, hybrid_ee, x_err, hybrid_ee_err);
  TGraphErrors *Hybrid_mumu = new TGraphErrors(6, x, hybrid_mumu, x_err, hybrid_mumu_err);
  
  TGraphErrors *Asym_ee = new TGraphErrors(6, x_asym, asym_ee, x_err, asym_ee_err);
  TGraphErrors *Asym_mumu = new TGraphErrors(6, x_asym, asym_mumu, x_err, asym_mumu_err);

  // -- For EE
  //Hybrid_ee -> SetTitle("Limit ee Asym. VS HybridNew");
  Hybrid_ee -> SetMarkerColor(4);
  Hybrid_ee -> SetLineColor(4);
  //Hybrid_ee -> SetMarkerStyle(24);
  Hybrid_ee -> Draw("sameP");
  //Asym_ee -> SetMarkerStyle(24);
  Asym_ee -> SetMarkerColor(kRed);
  Asym_ee -> SetLineColor(kRed);
  Asym_ee -> Draw("sameP");

  TLegend *lg = new TLegend(0.2, 0.8, 0.3, 0.95);
  lg->SetBorderSize(0);
  lg->SetFillStyle(0);
  lg->AddEntry(Hybrid_ee, "HybridNew", "l");
  lg->AddEntry(Asym_ee, "Asymptotic", "l");
  lg -> Draw("same");
  c1 -> SaveAs(base_plotpath + "/EE/Asym_vs_Hybridnew.pdf");

  // -- For MuMu
  template_h -> SetTitle("Limit mumu Asym. VS HybridNew");
  template_h -> GetYaxis() -> SetRangeUser(0.1, 1.5);
  template_h ->Draw();
  Hybrid_mumu -> SetMarkerColor(4);
  Hybrid_mumu -> SetLineColor(4);
  //Hybrid_mumu -> SetMarkerStyle(21);
  for(int i = 0; i < 6;i++){
    Hybrid_mumu -> GetXaxis()->SetBinLabel( i + 1, label[i]);
  }
  Hybrid_mumu -> Draw("sameP");
  //Asym_mumu -> SetMarkerStyle(21);
  Asym_mumu -> SetMarkerColor(kRed);
  Asym_mumu -> SetLineColor(kRed);
  Asym_mumu -> Draw("sameP");

  TLegend *lg_mumu = new TLegend(0.2, 0.8, 0.3, 0.95);
  lg_mumu->SetBorderSize(0);
  lg_mumu->SetFillStyle(0);
  lg_mumu->AddEntry(Hybrid_mumu, "HybridNew", "l ");
  lg_mumu->AddEntry(Asym_mumu, "Asymptotic", "l");
  lg_mumu -> Draw("same");
  c1 -> SaveAs(base_plotpath + "/MuMu/Asym_vs_Hybridnew.pdf");

}
