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

  mapfunc[histname + cyclename + samplename] = current_hist;

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
  maphist2D[histname + cyclename + samplename] = current_hist;
  
  current_file -> Close();
  delete current_file;
}

void draw_2D_DATA(TString nameofhistogram, TString current_data, int tag_year){
  
  cout << "[draw_2D_DATA] nameofhistogram :" << nameofhistogram << endl;
  double total_lumi = 35.9;
  if(tag_year == 2017) total_lumi = 41.3;
  if(tag_year == 2018) total_lumi = 59.7;

  TString pad1 = nameofhistogram;
  TString pad2 = nameofhistogram;
  TString canvas = nameofhistogram;
  TString hstack = nameofhistogram;
  TString legend = nameofhistogram;
  TString func = nameofhistogram;
  TString clone = nameofhistogram;
  TString line = nameofhistogram;
  pad1.Insert(0, "pad1_");
  pad2.Insert(0, "pad2_");
  canvas.Insert(0, "c_");
  hstack.Insert(0, "hs_");
  legend.Insert(0, "legend_");
  func.Insert(0, "ratio_");
  clone.Insert(0, "h3_");
  line.Insert(0, "l_");

  mapcanvas[canvas] = new TCanvas(canvas,"",800,800);
  canvas_margin(mapcanvas[canvas]);
  gStyle -> SetOptStat(1111);
  
  cout << "[draw_2D_DATA] nameofhistogram + Cycle_name + current_data : " << nameofhistogram + Cycle_name + current_data << endl;
  if(!maphist2D[nameofhistogram + Cycle_name + current_data]) return;
  maphist2D[nameofhistogram + Cycle_name + current_data] -> Scale(1. / total_lumi);
  //maphist2D[nameofhistogram + Cycle_name + current_data] -> SetStats(0);
  //maphist2D[nameofhistogram + Cycle_name + current_data] -> GetXaxis() -> SetTitle(nameofhistogram);
  //maphist2D[nameofhistogram + Cycle_name + current_data] -> GetYaxis() -> SetTitle(nameofhistogram);
  
  //maphist2D[nameofhistogram + Cycle_name + current_data] -> Draw("COLZTEXT");
  int bin[4];
  bin[0] = maphist2D[nameofhistogram + Cycle_name + current_data] -> FindBin(1.5, 1.5);
  bin[1] = maphist2D[nameofhistogram + Cycle_name + current_data] -> FindBin(-1.5, 1.5);
  bin[2] = maphist2D[nameofhistogram + Cycle_name + current_data] -> FindBin(-1.5, -1.5);
  bin[3] = maphist2D[nameofhistogram + Cycle_name + current_data] -> FindBin(1.5, -1.5);
  
  int N = 4;
  double bin_content[4], bin_error[4];
  double x[4] = {0.5, 1.5, 2.5, 3.5};
  double ex[4] = {0.5, 0.5, 0.5, 0.5};
  for(int i_2D = 0; i_2D < 4; i_2D++){
    cout << i_2D + 1 << " quadrature : " << maphist2D[nameofhistogram + Cycle_name + current_data] -> GetBinContent(bin[i_2D]) << endl;;
    bin_content[i_2D] = maphist2D[nameofhistogram + Cycle_name + current_data] -> GetBinContent(bin[i_2D]);
    bin_error[i_2D] = maphist2D[nameofhistogram + Cycle_name + current_data] -> GetBinError(bin[i_2D]);
  }
  TGraphErrors *gr = new TGraphErrors(4,x,bin_content,ex,bin_error);
  gr->SetTitle("");
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  
  std::vector<TString> x_labels = {"1st", "2nd", "3rd", "4th"};
  for(int i = 0; i < 4; i++){
    gr->GetXaxis()->SetBinLabel(gr->GetXaxis()->FindBin(i + 0.5), x_labels[i]);
  }
  
  double y_mean = 0., ye_mean = 0.;
  for(int i = 0; i < 4; i++){
    y_mean += bin_content[i];
    ye_mean += bin_error[i];
  }
  y_mean = y_mean / 4.;
  ye_mean = ye_mean / 4.;
  double y_range_bottom = y_mean - 10. * ye_mean;
  double y_range_top = y_mean + 10. * ye_mean;



  gr -> GetXaxis() -> SetTitleSize(0.15);
  gr -> GetXaxis() -> SetTitleOffset(0.5);
  gr -> GetXaxis() -> SetTitle("");
  gr -> GetXaxis() -> LabelsOption("h"); 
  gr -> GetYaxis() -> SetTitle("Count per 1 fb^{-1}");
  gr -> GetYaxis() -> SetTitleSize(0.07);
  gr -> GetYaxis() -> SetTitleOffset(1.02);
  gr -> GetYaxis() -> SetLabelSize(0.05);
  gr -> GetYaxis() -> SetRangeUser(y_range_bottom, y_range_top);


  gr->Draw("ALP");
  
  
  ////////////////////////////////////
  // -- Latex
  ////////////////////////////////////
  mapcanvas[canvas] -> cd();
  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_Lumi.SetTextSize(0.035);
  TString total_lumi_str = "41.3";
  if(tag_year == 2016) total_lumi_str = "35.9";
  if(tag_year == 2018) total_lumi_str = "59.7";
  latex_Lumi.DrawLatex(0.7, 0.96, total_lumi_str + " fb^{-1} (13 TeV)");

  

  TString pdfname;
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  pdfname = WORKING_DIR + "/plots/HEM_issue/";
  pdfname.Append(nameofhistogram);
  pdfname.Append("_1D.pdf");
  
  mapcanvas[canvas] -> SaveAs(pdfname);
  
}

void TwoD_Yield(int tag_year){
  
  SingleMuon = "data_SingleMuon";
  DoubleEG = "data_DoubleEG";
  if(tag_year == 2018) DoubleEG = "data_EGamma";
  
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
	
	openfile_DATA_2D(Cycle_name, DoubleEG, current_dir, current_2D_plot, tag_year);
	openfile_DATA_2D(Cycle_name, SingleMuon, current_dir, current_2D_plot, tag_year);
	
	if(channels[i_channels].Contains("DiEle")) draw_2D_DATA(current_2D_plot, DoubleEG, tag_year);
	else draw_2D_DATA(current_2D_plot, SingleMuon, tag_year);
	
	
      }
    }
  }
  
  return;
  
}

void HEM_issue(){
  
  Cycle_name = "HN_pair_all_SkimTree_LRSMHighPt";

  TwoD_Yield(2018);
  

}
