#include "canvas_margin.h"
#include "mylib.h"


void openfile_background(TString cyclename, TString samplename, TString dir, TString histname, int tag_year){

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/rootfiles/" +TString::Itoa(tag_year,10) + "/Background/";
  TString filename = root_file_path + cyclename + "_" + samplename + ".root";
  TFile *current_file = new TFile ((filename)) ;
  gDirectory->cd(dir);

  TH1F * current_hist = (TH1F*)gDirectory -> Get(histname);
  if(current_hist){
    current_hist -> SetDirectory(0);
  }
  TH1::AddDirectory(kFALSE);

  //if(samplename.Contains("DYJets") && current_hist) current_hist->Scale(DY_scale);                                                                                                                                                                                            
  mapfunc[histname + cyclename + samplename] = current_hist;

  current_file -> Close();
  delete current_file;
}

void openfile_DATA(TString cyclename, TString samplename, TString dir, TString histname, int tag_year){

  //cout << "[openfile_DATA] " << samplename << endl;                                                                                                                                                                                                                           
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/rootfiles/" +TString::Itoa(tag_year,10) + "/DATA/";
  TString filename = root_file_path + cyclename + "_" + samplename + ".root";
  TFile *current_file = new TFile ((filename)) ;
  gDirectory->cd(dir);

  TH1F * current_hist = (TH1F*)gDirectory -> Get(histname);
  if(current_hist){
    current_hist -> SetDirectory(0);
  }
  TH1::AddDirectory(kFALSE);
  
  mapfunc[histname + cyclename + samplename + "_" + TString::Itoa(tag_year,10)] = current_hist;
  
  current_file -> Close();
  delete current_file;
}

void openfile_background_2D(TString cyclename, TString samplename, TString dir, TString histname, int tag_year){

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/rootfiles/" +TString::Itoa(tag_year,10) + "/Background/";
  TString filename = root_file_path + cyclename + "_" + samplename + ".root";
  TFile *current_file = new TFile ((filename)) ;
  gDirectory->cd(dir);

  TH2D * current_hist = (TH2D*)gDirectory -> Get(histname);
  if(current_hist){
    current_hist -> SetDirectory(0);
  }
  TH1::AddDirectory(kFALSE);

  maphist2D[histname + cyclename + samplename] = current_hist;

  current_file -> Close();
  delete current_file;
}

void openfile_DATA_2D(TString cyclename, TString samplename, TString dir, TString histname, int tag_year){

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString root_file_path = WORKING_DIR+"/rootfiles/" +TString::Itoa(tag_year,10) + "/DATA/";
  TString filename = root_file_path + cyclename + "_" + samplename + ".root";
  TFile *current_file = new TFile ((filename)) ;
  gDirectory->cd(dir);

  TH2D * current_hist = (TH2D*)gDirectory -> Get(histname);
  if(current_hist){
    current_hist -> SetDirectory(0);
  }
  TH1::AddDirectory(kFALSE);
  cout << "[openfile_DATA_2D] " + histname + cyclename + samplename << endl;
  maphist2D[histname + cyclename + samplename + "_" + TString::Itoa(tag_year,10)] = current_hist;
  
  current_file -> Close();
  delete current_file;
}

void draw_2D_DATA(TString nameofhistogram, TString current_data){
  
  cout << "[draw_2D_DATA] nameofhistogram :" << nameofhistogram << endl;
  double total_lumi_2016 = 35.9;
  double total_lumi_2017 = 41.3;
  double total_lumi_2018 = 59.7;

  TCanvas *canvas = new TCanvas(nameofhistogram + "canvas","",800,800);
  canvas_margin(canvas);
  gStyle -> SetOptStat(1111);
  
  cout << "[draw_2D_DATA] nameofhistogram + Cycle_name + current_data : " << nameofhistogram + Cycle_name + current_data << endl;
  if(!maphist2D[nameofhistogram + Cycle_name + current_data + "_2016"]) return;
  if(!maphist2D[nameofhistogram + Cycle_name + current_data + "_2017"]) return;
  if(!maphist2D[nameofhistogram + Cycle_name + current_data + "_2018"]) return;
  
  /*
  maphist2D[nameofhistogram + Cycle_name + current_data + "_2016"] -> Scale(1. / total_lumi_2016);
  maphist2D[nameofhistogram + Cycle_name + current_data + "_2017"] -> Scale(1. / total_lumi_2017);
  maphist2D[nameofhistogram + Cycle_name + current_data + "_2018"] -> Scale(1. / total_lumi_2018);
  */
  //maphist2D[nameofhistogram + Cycle_name + current_data] -> SetStats(0);
  //maphist2D[nameofhistogram + Cycle_name + current_data] -> GetXaxis() -> SetTitle(nameofhistogram);
  //maphist2D[nameofhistogram + Cycle_name + current_data] -> GetYaxis() -> SetTitle(nameofhistogram);
  
  //maphist2D[nameofhistogram + Cycle_name + current_data] -> Draw("COLZTEXT");
  int bin_2016[4];
  bin_2016[0] = maphist2D[nameofhistogram + Cycle_name + current_data + "_2016"] -> FindBin(1.5, 1.5);
  bin_2016[1] = maphist2D[nameofhistogram + Cycle_name + current_data + "_2016"] -> FindBin(-1.5, 1.5);
  bin_2016[2] = maphist2D[nameofhistogram + Cycle_name + current_data + "_2016"] -> FindBin(-1.5, -1.5);
  bin_2016[3] = maphist2D[nameofhistogram + Cycle_name + current_data + "_2016"] -> FindBin(1.5, -1.5);
  
  int bin_2017[4];
  bin_2017[0] = maphist2D[nameofhistogram + Cycle_name + current_data + "_2017"] -> FindBin(1.5, 1.5);
  bin_2017[1] = maphist2D[nameofhistogram + Cycle_name + current_data + "_2017"] -> FindBin(-1.5, 1.5);
  bin_2017[2] = maphist2D[nameofhistogram + Cycle_name + current_data + "_2017"] -> FindBin(-1.5, -1.5);
  bin_2017[3] = maphist2D[nameofhistogram + Cycle_name + current_data + "_2017"] -> FindBin(1.5, -1.5);

  int bin_2018[4];
  bin_2018[0] = maphist2D[nameofhistogram + Cycle_name + current_data + "_2018"] -> FindBin(1.5, 1.5);
  bin_2018[1] = maphist2D[nameofhistogram + Cycle_name + current_data + "_2018"] -> FindBin(-1.5, 1.5);
  bin_2018[2] = maphist2D[nameofhistogram + Cycle_name + current_data + "_2018"] -> FindBin(-1.5, -1.5);
  bin_2018[3] = maphist2D[nameofhistogram + Cycle_name + current_data + "_2018"] -> FindBin(1.5, -1.5);
  
  
  int N = 4;
  double bin_content_2016[4], bin_error_2016[4];
  double bin_content_2017[4], bin_error_2017[4];
  double bin_content_2018[4], bin_error_2018[4];

  double x[4] = {0.5, 1.5, 2.5, 3.5};
  double ex[4] = {0.5, 0.5, 0.5, 0.5};
  for(int i_2D = 0; i_2D < 4; i_2D++){
    //cout << i_2D + 1 << " quadrature : " << maphist2D[nameofhistogram + Cycle_name + current_data] -> GetBinContent(bin[i_2D]) << endl;;
    bin_content_2016[i_2D] = maphist2D[nameofhistogram + Cycle_name + current_data + "_2016"] -> GetBinContent(bin_2016[i_2D]);
    bin_content_2017[i_2D] = maphist2D[nameofhistogram + Cycle_name + current_data + "_2017"] -> GetBinContent(bin_2017[i_2D]);
    bin_content_2018[i_2D] = maphist2D[nameofhistogram + Cycle_name + current_data + "_2018"] -> GetBinContent(bin_2018[i_2D]);

    bin_error_2016[i_2D] = maphist2D[nameofhistogram + Cycle_name + current_data + "_2016"] -> GetBinError(bin_2016[i_2D]);
    bin_error_2017[i_2D] = maphist2D[nameofhistogram + Cycle_name + current_data + "_2017"] -> GetBinError(bin_2017[i_2D]);
    bin_error_2018[i_2D] = maphist2D[nameofhistogram + Cycle_name + current_data + "_2018"] -> GetBinError(bin_2018[i_2D]);
    
  }
  
  double denom_2016 = bin_content_2016[0];
  double denom_2017 = bin_content_2017[0];
  double denom_2018 = bin_content_2018[0];

  for(int i_2D = 0; i_2D < 4; i_2D++){
    cout << i_2D + 1 <<" quadrature : " << bin_content_2016[i_2D] << endl;
    bin_content_2016[i_2D] = bin_content_2016[i_2D] / denom_2016;
    bin_content_2017[i_2D] = bin_content_2017[i_2D] / denom_2017;
    bin_content_2018[i_2D] = bin_content_2018[i_2D] / denom_2018;
    
    bin_error_2016[i_2D] = bin_error_2016[i_2D] / denom_2016;
    bin_error_2017[i_2D] = bin_error_2017[i_2D] / denom_2017;
    bin_error_2018[i_2D] = bin_error_2018[i_2D] / denom_2018;
    
    cout << i_2D + 1 << " quadrature : " << bin_content_2016[i_2D] << endl;
  }
  
  TGraphErrors *gr_2016 = new TGraphErrors(4,x,bin_content_2016,ex,bin_error_2016);
  gr_2016->SetTitle("");
  gr_2016->SetMarkerColor(4);
  gr_2016->SetMarkerStyle(21);
  
  std::vector<TString> x_labels = {"1st", "2nd", "3rd", "4th"};
  for(int i = 0; i < 4; i++){
    gr_2016->GetXaxis()->SetBinLabel(gr_2016->GetXaxis()->FindBin(i + 0.5), x_labels[i]);
  }
  
  double max_y = 0.;
  double min_y = 99999999.;
  for(int i = 0; i < 4; i++){
    if(max_y < bin_content_2016[i]) max_y = bin_content_2016[i];
    if(min_y > bin_content_2016[i]) min_y = bin_content_2016[i];
  }
  for(int i = 0; i < 4;i++){
    if(max_y < bin_content_2017[i]) max_y = bin_content_2017[i];
    if(min_y > bin_content_2017[i]) min_y = bin_content_2017[i];
  }for(int i = 0; i < 4;i++){
    if(max_y < bin_content_2018[i]) max_y = bin_content_2018[i];
    if(min_y > bin_content_2018[i]) min_y = bin_content_2018[i];
  }
  
  double y_range_bottom = 0.;
  double y_range_top = max_y * 2.;

  gr_2016 -> SetMarkerColor(kBlue);
  gr_2016 -> SetLineColor(kBlue);
  gr_2016 -> GetXaxis() -> SetTitleSize(0.15);
  gr_2016 -> GetXaxis() -> SetTitleOffset(0.5);
  gr_2016 -> GetXaxis() -> SetTitle("");
  gr_2016 -> GetXaxis() -> LabelsOption("h"); 
  gr_2016 -> GetYaxis() -> SetTitle("Count / 1st Quadrature");
  gr_2016 -> GetYaxis() -> SetTitleSize(0.07);
  gr_2016 -> GetYaxis() -> SetTitleOffset(1.02);
  gr_2016 -> GetYaxis() -> SetLabelSize(0.05);
  gr_2016 -> GetYaxis() -> SetRangeUser(y_range_bottom, y_range_top);


  gr_2016->Draw("ALP");


  TGraphErrors *gr_2017 = new TGraphErrors(4,x,bin_content_2017,ex,bin_error_2017);
  gr_2017 -> SetMarkerStyle(21);
  gr_2017 -> SetMarkerColor(kRed);
  gr_2017 -> SetLineColor(kRed);
  gr_2017 -> Draw("LPsame");

  TGraphErrors *gr_2018 = new TGraphErrors(4,x,bin_content_2018,ex,bin_error_2018);
  gr_2018 -> SetMarkerStyle(21);
  gr_2018 -> SetMarkerColor(kBlack);
  gr_2018 -> SetLineColor(kBlack);
  gr_2018 -> Draw("LPsame");
  
  TLegend *legend = new TLegend(0.60, 0.75, 0.95, 0.95);
  legend -> AddEntry(gr_2016, "2016", "lp");
  legend -> AddEntry(gr_2017, "2017", "lp");
  legend -> AddEntry(gr_2018, "2018", "lp");

  legend -> SetFillColor(kWhite);
  legend -> SetLineColor(kBlack);
  legend -> SetBorderSize(1);
  legend -> SetFillStyle(1001);
  legend -> SetShadowColor(0); 
  legend -> SetEntrySeparation(0.3);
  legend -> Draw("same");
  
  ////////////////////////////////////
  // -- Latex
  ////////////////////////////////////
  canvas -> cd();
  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_Lumi.SetTextSize(0.035);
  TString total_lumi_str = "41.3";
  if(tag_year == 2016) total_lumi_str = "35.9";
  if(tag_year == 2018) total_lumi_str = "59.7";
  latex_Lumi.DrawLatex(0.7, 0.96, "(13 TeV)");
  
  TString pdfname;
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  pdfname = WORKING_DIR + "/plots/HEM_issue/";
  pdfname.Append(nameofhistogram);
  pdfname.Append("_1D.pdf");
  
  canvas -> SaveAs(pdfname);
  
}

void draw_1D_DATA(TString nameofhistogram, TString current_data){

  

}

void TwoD_Yield(){
  
  int N_2D_plst = 6;
  TString list_2D_lots[] = {"FatJet_0phi_VS_eta", "FatJet_1phi_VS_eta",
			    "Jet_0phi_VS_eta", "Jet_1phi_VS_eta",
			    "Lepton_0phi_VS_eta", "Lepton_1phi_VS_eta",
  };

  int N_regions = 3;
  TString regions[] = {"CR_Zmass",
                       "CR_ttbar",
                       //"CR_inv_mll",
		       "SR"
  };


  int N_channels = 3;
  TString channels[] = {"DiEle",
                        "DiMu",
                        "EMu"
  };
  
  for(int i_list_2D = 0; i_list_2D < N_2D_plst; i_list_2D++){
    for(int i_regions = 0; i_regions < N_regions; i_regions++){
      for(int i_channels = 0; i_channels < N_channels; i_channels++){

	TString current_2D_plot = list_2D_lots[i_list_2D] + "_" + regions[i_regions] + "_" + channels[i_channels] + "_central";
	cout << current_2D_plot << endl;
	
	TString current_dir = "empty";
	TString current_channel = "empty";
	
	current_dir = regions[i_regions];
	current_channel = channels[i_channels];
	current_dir = current_dir + "_" + current_channel + "_central";
	
	SingleMuon = "data_SingleMuon";
        DoubleEG = "data_DoubleEG";

	openfile_DATA_2D(Cycle_name, DoubleEG, current_dir, current_2D_plot, 2016);
	openfile_DATA_2D(Cycle_name, SingleMuon, current_dir, current_2D_plot, 2016);

	openfile_DATA_2D(Cycle_name, DoubleEG, current_dir, current_2D_plot, 2017);
	openfile_DATA_2D(Cycle_name, SingleMuon, current_dir, current_2D_plot, 2017);
	
	openfile_DATA_2D(Cycle_name, DoubleEG, current_dir, current_2D_plot, 2018);
        openfile_DATA_2D(Cycle_name, SingleMuon, current_dir, current_2D_plot, 2018);
	
	if(channels[i_channels].Contains("DiEle")) draw_2D_DATA(current_2D_plot, DoubleEG);
	else draw_2D_DATA(current_2D_plot, SingleMuon);
	
      }
    }
  }
  
  return;
  
}


void OneD_Kinematics(){

  int N_1D_plst = 8;
  TString list_1D_lots[] = {"FatJet_0_Pt_1st_Quadrant", 
			    "FatJet_0_Pt_2nd_Quadrant",
			    "FatJet_0_Pt_3rd_Quadrant",
			    "FatJet_0_Pt_4th_Quadrant",
			    "Jet_0_Pt_1st_Quadrant",
			    "Jet_0_Pt_2nd_Quadrant",
			    "Jet_0_Pt_3rd_Quadrant",
			    "Jet_0_Pt_4th_Quadrant",
  };

  int N_regions = 3;
  TString regions[] = {"CR_Zmass",
		       "CR_ttbar",
		       //"CR_inv_mll",
                       "SR"
  };


  int N_channels = 3;
  TString channels[] = {"DiEle",
                        "DiMu",
                        "EMu"
  };
  
  for(int i_list_1D = 0; i_list_1D < N_1D_plst; i_list_1D++){
    for(int i_regions = 0; i_regions < N_regions; i_regions++){
      for(int i_channels = 0; i_channels < N_channels; i_channels++){

	TString current_1D_plot = list_1D_lots[i_list_1D] + "_" + regions[i_regions] + "_" + channels[i_channels] + "_central";
        cout << current_1D_plot << endl;

	TString current_dir = "empty";
	TString current_channel = "empty";

	current_dir = regions[i_regions];
	current_channel = channels[i_channels];
	current_dir = current_dir + "_" + current_channel + "_central";

	SingleMuon = "data_SingleMuon";
	DoubleEG = "data_DoubleEG";
	
	openfile_DATA(Cycle_name, DoubleEG, current_dir, current_1D_plot, 2016);
	openfile_DATA(Cycle_name, SingleMuon, current_dir, current_1D_plot, 2016);

	openfile_DATA(Cycle_name, DoubleEG, current_dir, current_1D_plot, 2017);
	openfile_DATA(Cycle_name, SingleMuon, current_dir, current_1D_plot, 2017);

	openfile_DATA(Cycle_name, DoubleEG, current_dir, current_1D_plot, 2018);
	openfile_DATA(Cycle_name, SingleMuon, current_dir, current_1D_plot, 2018);
	
        if(channels[i_channels].Contains("DiEle")) draw_1D_DATA(current_1D_plot, DoubleEG);
        else draw_1D_DATA(current_1D_plot, SingleMuon);
	
      }
    }
  }

  return;
  


}

void HEM_issue(){
  
  Cycle_name = "HN_pair_all_SkimTree_LRSMHighPt";

  TwoD_Yield();
  OneD_Kinematics();

}
