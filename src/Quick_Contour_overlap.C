#include "canvas_margin.h"
#include "mylib.h"

void openfile(TString binning, TString channel){
  
  cout << "[openfile_background] " << binning + "_" + channel << endl;

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/plots/limit/" + binning + "/" + channel + "/";
  TString filename = root_file_path + "Limits.root";
  TFile *current_file = new TFile ((filename)) ;
  TString gr_name = "Cotour";
  TGraph * current_gr = (TGraph*)gDirectory -> Get(gr_name);
  TH1::AddDirectory(kFALSE);
  
  //if(samplename.Contains("DYJets") && current_hist) current_hist->Scale(DY_scale);
  map_gr[binning + channel] = current_gr;
  
  current_file -> Close();
  delete current_file;
}

void make_histogram(TString channel){
  TCanvas *c_ct = new TCanvas("c_ct", "", 750, 600);
  canvas_margin(c_ct);
  c_ct->cd();

  TH1D *dummy_ct_ForEachZP = new TH1D("dummy_ct_ForEachZP", "", 6000, 0., 6000.);
  hist_axis(dummy_ct_ForEachZP);
  dummy_ct_ForEachZP->GetYaxis()->SetLabelSize(0.04);
  dummy_ct_ForEachZP->GetYaxis()->SetTitleSize(0.06);
  dummy_ct_ForEachZP->GetYaxis()->SetTitleOffset(1.25);
  dummy_ct_ForEachZP->Draw("hist");
  dummy_ct_ForEachZP->GetYaxis()->SetRangeUser(100., 2500.);
  dummy_ct_ForEachZP->GetYaxis()->SetTitle("m_{N} (GeV)");
  dummy_ct_ForEachZP->GetXaxis()->SetRangeUser(400., 5000.);
  dummy_ct_ForEachZP->GetXaxis()->SetTitle("m_{Z'} (GeV)");


  double x_Z2N[2] = {0,10000};
  double y_Z2N[2] = {0,5000};
  TGraph *gr_Z2N = new TGraph(2,x_Z2N,y_Z2N);
  gr_Z2N->SetLineStyle(4);
  gr_Z2N->SetLineWidth(5);
  gr_Z2N->Draw("lsame");

  Int_t signal_colour_array[] = {868, 632, 416, 400};

  TLegend *lg = new TLegend(0.2, 0.7, 0.6, 0.90);
  lg->SetBorderSize(0);
  lg->SetFillStyle(0);
  lg->AddEntry(gr_Z2N, "m_{N} = m_{Z'} / 2 ", "l");

  TString binnings = "";
  for(unsigned int i = 0; i < map_sample_names["limit"].size(); i++){
    map_gr[map_sample_names["limit"].at(i) + channel] -> SetLineColor(signal_colour_array[i]);
    map_gr[map_sample_names["limit"].at(i) + channel] -> SetLineStyle(1);
    map_gr[map_sample_names["limit"].at(i) + channel] -> Draw("lsame");
    if(i == 0){
      map_gr[map_sample_names["limit"].at(i) + channel] -> SetFillColor(signal_colour_array[i]);
      map_gr[map_sample_names["limit"].at(i) + channel] -> Draw("fsame");
    }
    lg -> AddEntry(map_gr[map_sample_names["limit"].at(i) + channel], map_sample_names["limit"].at(i), "l");

    binnings = binnings + map_sample_names["limit"].at(i);
  }
  
  lg->Draw("same");

  TLatex latex_CMSPriliminary, latex_Lumi, latex_channel;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.16, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_Lumi.SetTextSize(0.035 );
  latex_Lumi.DrawLatex(0.60, 0.96, "137.4 fb^{-1} (13 TeV, Full Run2)");
  latex_channel.SetNDC();
  latex_channel.SetTextSize(0.08);
  if(channel.Contains("MuMu")) latex_channel.DrawLatex(0.2, 0.6, "#mu#mu OS + SS");
  else latex_channel.DrawLatex(0.2, 0.6, "ee OS + SS");
  dummy_ct_ForEachZP -> Draw("axissame");
  
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

void Quick_Contour_overlap(int year=2018){
  setTDRStyle();
  
  // == Set Input Sample List
  tag_year = year;
  
  map_sample_names["limit"] = {"3AK8_CRs", "0AK8", "1AK8", "2AK8"};
  //map_sample_names["limit"] = {"3AK8", "4AK8"};

  open_files();
}
