#include "canvas_margin.h"
#include "mylib.h"

void openfile(TString binning, TString channel){
  
  cout << "[openfile] " << binning + "_" + channel << endl;
  gStyle->SetLineWidth(2);

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/plots/limit/" + binning + "/" + channel + "/";
  TString filename = root_file_path + "Limits.root";
  TFile *current_file = new TFile ((filename)) ;
  TString gr_name = "Contour";
  TGraph * current_gr = (TGraph*)gDirectory -> Get(gr_name);
  TGraph * current_gr_up = (TGraph*)gDirectory -> Get(gr_name + "_up");
  TGraph * current_gr_down = (TGraph*)gDirectory -> Get(gr_name + "_down");
  TGraph * current_gr_obs = (TGraph*)gDirectory -> Get(gr_name + "_obs");
  TGraph * current_gr_atlas = (TGraph*)gDirectory -> Get(gr_name + "_ATLAS");

  TH1::AddDirectory(kFALSE);
  
  map_gr[binning + channel] = current_gr;
  map_gr[binning + channel + "_up"] = current_gr_up;
  map_gr[binning + channel + "_down"] = current_gr_down;
  map_gr[binning + channel + "_obs"] = current_gr_obs;
  map_gr[binning + channel + "_atlas"] = current_gr_atlas;

  current_file -> Close();
  delete current_file;
}

void make_histogram(TString channel){
  TCanvas *c_ct = new TCanvas("c_ct", "", 750, 600);
  canvas_margin(c_ct);
  c_ct->cd();

  TH1D *dummy_ct_ForEachZP = new TH1D("dummy_ct_ForEachZP", "", 6000, 0., 6000.);
  hist_axis(dummy_ct_ForEachZP);
  dummy_ct_ForEachZP->GetYaxis()->SetLabelSize(0.05);
  dummy_ct_ForEachZP->GetYaxis()->SetTitleSize(0.06);
  dummy_ct_ForEachZP->GetYaxis()->SetTitleOffset(1.25);
  dummy_ct_ForEachZP->Draw("hist");
  dummy_ct_ForEachZP->GetYaxis()->SetRangeUser(100., 2500.);
  dummy_ct_ForEachZP->GetYaxis()->SetTitle("m_{N} (GeV)");
  dummy_ct_ForEachZP->GetXaxis()->SetLabelSize(0.045);
  dummy_ct_ForEachZP->GetXaxis()->SetRangeUser(400., 5000.);
  dummy_ct_ForEachZP->GetXaxis()->SetTitleSize(0.06);
  dummy_ct_ForEachZP->GetXaxis()->SetTitle("m_{Z'} (GeV)");


  double x_Z2N[2] = {0,10000};
  double y_Z2N[2] = {0,5000};
  TGraph *gr_Z2N = new TGraph(2,x_Z2N,y_Z2N);
  gr_Z2N->SetLineStyle(4);
  gr_Z2N->SetLineWidth(5);
  gr_Z2N->Draw("lsame");

  // -- 868 = kAzure + 8, 632 = kRed, 416 = kGreen, 400 = kYellow, 800 = kOrange
  Int_t signal_colour_array[] = {868, 400, 416, 800};

  TLegend *lg = new TLegend(0.2, 0.52, 0.7, 0.93);
  lg->SetBorderSize(0);
  lg->SetFillStyle(0);
  //lg->AddEntry(gr_Z2N, "m_{N} = m_{Z'} / 2 ", "l");

  TGraph *gr_central_fill = (TGraph*)map_gr[map_sample_names["limit"].at(0) + channel]->Clone();
  gr_central_fill -> SetFillColor(kAzure + 8);
  gr_central_fill -> Draw("fsame");

  TString binnings = "";
  map_gr[map_sample_names["limit"].at(0) + channel] -> SetLineColor(kBlue);
  map_gr[map_sample_names["limit"].at(0) + channel] -> SetLineStyle(8);
  map_gr[map_sample_names["limit"].at(0) + channel] -> SetLineWidth(3);
  map_gr[map_sample_names["limit"].at(0) + channel] -> Draw("lsame");

  map_gr[map_sample_names["limit"].at(0) + channel + "_up"] -> SetLineColor(kBlue);
  map_gr[map_sample_names["limit"].at(0) + channel + "_up"] -> SetLineStyle(2);
  map_gr[map_sample_names["limit"].at(0) + channel + "_up"] -> SetLineWidth(1);
  map_gr[map_sample_names["limit"].at(0) + channel + "_up"] -> Draw("lsame");

  map_gr[map_sample_names["limit"].at(0) + channel + "_down"] -> SetLineColor(kBlue);
  map_gr[map_sample_names["limit"].at(0) + channel + "_down"] -> SetLineStyle(2);
  map_gr[map_sample_names["limit"].at(0) + channel + "_down"] -> SetLineWidth(1);
  map_gr[map_sample_names["limit"].at(0) + channel + "_down"] -> Draw("lsame");

  lg -> AddEntry(map_gr[map_sample_names["limit"].at(0) + channel], "Combined (exp.)", "l");
  lg -> AddEntry(map_gr[map_sample_names["limit"].at(0) + channel + "_up"], "Combined (exp. #pm s.d.)", "l");
  lg -> AddEntry(map_gr[map_sample_names["limit"].at(0) + channel + "_obs"], "Combined (obs.)", "l");
  
  // -- Draw 0AK8, 1AK8, 2AK8 observed
  for(unsigned int i = 1; i < map_sample_names["limit"].size(); i++){
    cout << map_sample_names["limit"].at(i) + channel + "_obs" << endl;
    map_gr[map_sample_names["limit"].at(i) + channel + "_obs"] -> SetLineColor(signal_colour_array[i]);
    map_gr[map_sample_names["limit"].at(i) + channel + "_obs"] -> SetLineStyle(1);
    map_gr[map_sample_names["limit"].at(i) + channel + "_obs"] -> Draw("lsame");
    
    lg -> AddEntry(map_gr[map_sample_names["limit"].at(i) + channel + "_obs"], map_sample_names["limit"].at(i) + " (obs.)", "l");
    binnings = binnings + map_sample_names["limit"].at(i);
  }

  map_gr[map_sample_names["limit"].at(0) + channel + "_atlas"] -> SetLineStyle(8);
  map_gr[map_sample_names["limit"].at(0) + channel + "_atlas"] -> SetLineWidth(3);
  map_gr[map_sample_names["limit"].at(0) + channel + "_atlas"] -> SetFillColor(kAzure + 8);
  map_gr[map_sample_names["limit"].at(0) + channel + "_atlas"] -> SetLineColor(kRed);
  //map_gr[map_sample_names["limit"].at(0) + channel + "_atlas"] -> Draw("lsame");
  //lg -> AddEntry(map_gr[map_sample_names["limit"].at(0) + channel + "_atlas"], "ATLAS 8 TeV", "l");
   
  map_gr[map_sample_names["limit"].at(0) + channel + "_obs"] -> SetLineWidth(3);
  map_gr[map_sample_names["limit"].at(0) + channel + "_obs"] -> SetLineColor(kBlack);
  map_gr[map_sample_names["limit"].at(0) + channel + "_obs"] -> Draw("lsame");


  lg->Draw("same");

  TLatex latex_CMSPriliminary, latex_Lumi, latex_channel, latex_xy;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  //latex_CMSPriliminary.DrawLatex(0.16, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_CMSPriliminary.DrawLatex(0.16, 0.96, "#font[62]{CMS}");
  latex_Lumi.SetTextSize(0.035 );
  latex_Lumi.DrawLatex(0.60, 0.96, "138 fb^{-1} (13 TeV, Full Run2)");
  latex_channel.SetNDC();
  latex_channel.SetTextSize(0.05);
  if(channel.Contains("MuMu")) latex_channel.DrawLatex(0.22, 0.45, "#mu#mu OS + SS");
  else latex_channel.DrawLatex(0.22, 0.45, "#font[42]{ee} OS + SS");
  dummy_ct_ForEachZP -> Draw("axissame");
  latex_xy.SetNDC();
  latex_xy.SetTextSize(0.05);
  latex_xy.DrawLatex(0.70, 0.65, "m_{N} = m_{Z'} / 2");

  TString pdfname;
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  pdfname = WORKING_DIR + "/plots/limit/comparison/";
  pdfname = pdfname + binnings + "_" + channel + ".pdf";
  c_ct->SaveAs(pdfname);
    
}


void open_files(){
  
  for(unsigned int i = 0; i < map_sample_names["limit"].size(); i++){
    openfile(map_sample_names["limit"].at(i), "EE");
    openfile(map_sample_names["limit"].at(i), "MuMu");
  }

  make_histogram("EE");
  make_histogram("MuMu");

}

void Quick_Contour_overlap_obs(int year=2018){
  setTDRStyle();
  
  // == Set Input Sample List
  tag_year = year;
  
  map_sample_names["limit"] = {"3AK8_CRs", "0AK8", "1AK8", "2AK8"};
  //map_sample_names["limit"] = {"3AK8", "4AK8"};

  open_files();
}
