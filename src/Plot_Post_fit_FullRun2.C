#include "canvas_margin.h"
#include "mylib.h"


void openfile(TString plot_name, TString which_fit, TString mass_point){
  // -- plot_name = 2AK8_SR_EMu ..., which_fit = prefit, postfit
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString NAK8 = "";
  if(plot_name.Contains("0AK8")){
    NAK8 = "0AK8";
  }
  if(plot_name.Contains("1AK8")){
    NAK8 = "1AK8";
  }
  if(plot_name.Contains("2AK8")){
    NAK8 = "2AK8";
  }
  
  TString filename = WORKING_DIR+"/rootfiles/fitDiagnostics/output_" + plot_name + "_" + mass_point + ".root";
  if(plot_name.Contains("SR")){
    if(plot_name.Contains("DiEle") || plot_name.Contains("DiMu")){
      filename = WORKING_DIR+"/rootfiles/fitDiagnostics/output_mZp_"+ NAK8+ "_SR_"  + mass_point + ".root";
    }
  }
  cout << "[openfile] " << filename << endl;
  //TString filename = WORKING_DIR+"/rootfiles/fitDiagnostics/fitDiagnostics_fit_result_100fake.root";

  TFile *current_file = new TFile ((filename)) ;
  //gDirectory->cd(which_fit);

  TString histograms[4] = {"DYJets_MG_HT", "TT_powheg", "VV", "fake"};
  TString years[3] = {"2016", "2017", "2018"};
  
  
  // -- Get MC plots
  for(int i_year = 0; i_year < 3; i_year++){
    TString current_year = years[i_year];
    
    gDirectory->cd(plot_name + "_" + current_year + "_" + which_fit);
    
    for(int i_hist = 0; i_hist < 4; i_hist++){
      TH1F * current_hist = (TH1F*)gDirectory -> Get(plot_name + "_DYreweight_" + histograms[i_hist]);
      if(current_hist){
	current_hist -> SetDirectory(0);
      }
      TH1::AddDirectory(kFALSE);
      
      mapfunc[plot_name + histograms[i_hist] + which_fit + current_year] = current_hist;
    }  
    gDirectory -> cd("../");
  }

  // -- Get total bkg/data
  gDirectory -> cd(which_fit);
  TH1F * current_bkg = (TH1F*)gDirectory -> Get("TotalBkg");
  TH1F * current_data = (TH1F*)gDirectory -> Get("data_obs");
  if(current_bkg){
    current_bkg -> SetDirectory(0);
  }
  TH1::AddDirectory(kFALSE);
  mapfunc[plot_name + "bkg" + which_fit] = current_bkg;
  if(current_data){
    current_data -> SetDirectory(0);
  }
  TH1::AddDirectory(kFALSE);
  mapfunc[plot_name + "data" + which_fit] = current_data;

  current_file -> Close();
  delete current_file;
}

void make_histogram(TString plot_name, TString which_fit, int N_bin, double binx[]){

  TString nameofhistogram = plot_name + "_" + which_fit;
  
  TString title_y;
  if(nameofhistogram.Contains("N")) title_y = "Number";
  else title_y = "Events";

  bool blind = false;
  blind = (nameofhistogram.Contains("SR")) && (!nameofhistogram.Contains("EMu"));
  
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

  ////////////////////////////////////
  /// Pad 1
  ////////////////////////////////////
  mappad[pad1] = new TPad("", "", 0, 0.25, 1, 1);
  mappad[pad1] -> SetTopMargin( 0.07 );
  mappad[pad1] -> SetBottomMargin( 0.05 );
  mappad[pad1] -> SetLeftMargin( 0.15 );
  mappad[pad1] -> SetRightMargin( 0.04 );
  mappad[pad1] -> Draw();
  mappad[pad1] -> cd();
  mappad[pad1] -> SetLogy();
  
  maplegend[legend] = new TLegend(0.55, 0.45, 0.95, 0.90);

  if(debug) cout << "1" << endl;
  

  TH1F * pad1_template = new TH1F("", "", N_bin-1, binx);
  gStyle->SetOptTitle(0);
  gStyle->SetLineWidth(2);
  // -- Bkgs
  //TString histograms[4] = {"VV", "WJets_MG_HT", "TT_powheg", "DYJets_MG_HT"};
  //TString legend_list[4] = {"Other", "WJets", "ttbar", "DYJets"};
  TString histograms[4] = {"VV", "fake", "TT_powheg", "DYJets_MG_HT"};
  TString legend_list[4] = {"Other", "Nonprompt", "ttbar", "DYJets"};
  TString years[3] = {"2016", "2017", "2018"};

  //Int_t colour_array[] = {870, 632, 416, 400}; // -- Original
  //Int_t colour_array[] = {870, 901, 416, 400}; // -- For John's comment
  //Int_t colour_array[] = {869, 901, 416, 400}; // -- Test for graylish blue
  Int_t colour_array[] = {867, 901, 416, 393}; // -- Test for, more pastel Blue, darker Yellow and Green
  maphstack[hstack] = new THStack(hstack, "Stacked_" + nameofhistogram);
  
  for(int i_hist = 0; i_hist < 4; i_hist++){
    TString current_hist = plot_name + histograms[i_hist] + which_fit;
    mapfunc[current_hist + "sum"] = new TH1F(current_hist + "sum", current_hist + "sum", N_bin-1, binx);
    for(int i_year = 0;i_year < 3; i_year++){
      if(mapfunc[current_hist + years[i_year]]){
	for(int i_bin = 1; i_bin < N_bin; i_bin++){
	  double current_content = mapfunc[current_hist + years[i_year]] -> GetBinContent(i_bin) + mapfunc[current_hist + "sum"]  -> GetBinContent(i_bin);
	  mapfunc[current_hist + "sum"]  -> SetBinContent(i_bin, current_content);
	}
      }  
    }
    mapfunc[current_hist + "sum"] -> SetFillColor(colour_array[i_hist]);
    mapfunc[current_hist + "sum"] -> SetLineColor(colour_array[i_hist]);
    maphstack[hstack] -> Add(mapfunc[current_hist + "sum"]);
  }
  
  for(int i_legend = 0; i_legend < 4; i_legend++){
    TString current_hist = plot_name + histograms[4 - i_legend - 1] + which_fit;
    if(mapfunc[current_hist + "sum"]){
      maplegend[legend] -> AddEntry(mapfunc[current_hist + "sum"] , legend_list[4 - i_legend - 1], "lf");
    }
  }

  // -- Bkg sums
  TString current_bkg_sums = plot_name + "bkg" + which_fit;
  TH1F* bkg_sum = new TH1F(current_bkg_sums + "rebin", current_bkg_sums + "rebin", N_bin-1, binx);
  for(int i_bin = 1; i_bin < N_bin; i_bin++){
    double current_content = mapfunc[current_bkg_sums] -> GetBinContent(i_bin);
    double current_error = mapfunc[current_bkg_sums] -> GetBinError(i_bin);
    bkg_sum -> SetBinContent(i_bin, current_content);
    bkg_sum -> SetBinError(i_bin, current_error);
  }
  
  // -- Data
  double max_data = 0.;
  double last_bin_data = 0.;
  TString current_data = plot_name + "data" + which_fit;
  TH1F* data_hist = new TH1F(current_data + "rebin", current_data + "rebin", N_bin-1, binx);
  for(int i_bin = 1; i_bin < N_bin; i_bin++){
    double current_content = mapfunc[current_data] -> GetBinContent(i_bin);
    double current_error = mapfunc[current_data] -> GetBinError(i_bin);
    data_hist -> SetBinContent(i_bin, current_content);
    data_hist -> SetBinError(i_bin, current_error);
    if(current_content > max_data) max_data = current_content;
  }

  map_asym_gr[current_data] = new TGraphAsymmErrors(data_hist);
  for(int i = 0; i < N_bin; i++){
    int N = data_hist -> GetBinContent(i + 1);
    double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
    double U =  (N==0) ? ( ROOT::Math::gamma_quantile_c(alpha,N+1,1) ) : ( ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) );
    if( N!=0 ){
      map_asym_gr[current_data] -> SetPointEYlow(i, (N-L) );
      map_asym_gr[current_data] -> SetPointEXlow(i, 0);
      map_asym_gr[current_data] -> SetPointEYhigh(i, (U-N) );
      map_asym_gr[current_data] -> SetPointEXhigh(i, 0);
      
      double current_x = -1., current_y = -1.;
      map_asym_gr[current_data] -> GetPoint(i, current_x, current_y);
      double modified_y = -1.;
      modified_y = current_y;
      map_asym_gr[current_data] -> SetPoint(i, current_x, modified_y);
      
      veyl.push_back( (N-L));
      veyh.push_back( (U-N));
    }
    else{
      double zerodata_err_low = 0.1 ;
      double zerodata_err_high = 1.8 ;

      veyl.push_back(zerodata_err_low);
      veyh.push_back(zerodata_err_high);

      double current_x = data_hist -> GetBinCenter(i + 1);
      double current_ex = data_hist -> GetBinWidth(i + 1) /2.;

      map_asym_gr[current_data] -> SetPointEYlow(i, zerodata_err_low);
      map_asym_gr[current_data] -> SetPointEXlow(i, 0.);
      map_asym_gr[current_data] -> SetPointEYhigh(i, zerodata_err_high);
      map_asym_gr[current_data] -> SetPointEXhigh(i, 0.);

      vx.push_back(current_x);
      vexl.push_back(current_ex);
      vexh.push_back(current_ex);
    }
  }
  
  const int N_points = N_bin - 1;
  double data_y[N_points], data_x[N_points];
  double data_ey_high[N_points], data_ey_low[N_points];
  double data_ex_high[N_points], data_ex_low[N_points];
  for(int i_bin = 1; i_bin < N_bin; i_bin++){
    double current_data = data_hist -> GetBinContent(i_bin);
    double current_data_ey_high = veyh.at(i_bin-1);
    double current_data_ey_low = veyl.at(i_bin-1);

    data_y[i_bin -1] = current_data;
    data_ey_high[i_bin - 1] = current_data_ey_high;
    data_ey_low[i_bin - 1] = current_data_ey_low;
  }
  for(int i_x = 0; i_x < N_points; i_x++){
    double current_ex = (binx[i_x + 1] - binx[i_x]) / 2.0;
    double current_x = binx[i_x] + current_ex;

    data_x[i_x] = current_x;
    data_ex_high[i_x] = 0.;
    data_ex_low[i_x] = 0.;
  }
  last_bin_data = data_y[N_points - 1];
  TGraphAsymmErrors * data_gr = new TGraphAsymmErrors(N_points, data_x, data_y, data_ex_low, data_ex_high, data_ey_low, data_ey_high);

  maplegend[legend] -> AddEntry(bkg_sum, "Stat.+syst. uncert.", "f");
  maplegend[legend] -> AddEntry(data_gr, "Observed", "lp");

  // -- Draw Pad
  pad1_template -> Draw("hist");
  pad1_template -> GetYaxis()->SetLabelSize(0.07);
  pad1_template -> GetYaxis()->SetTitleSize(0.07);
  pad1_template -> GetYaxis()->SetTitleOffset(1.02);
  pad1_template -> GetXaxis()->SetLabelSize(0);
  pad1_template -> GetXaxis() -> SetTitle(nameofhistogram);
  //pad1_template -> GetYaxis() -> SetRangeUser(0., max_data * 1.8);
  //pad1_template -> GetYaxis() -> SetRangeUser(0.1, max_data * 10.); // logy
  pad1_template -> GetYaxis() -> SetRangeUser(0.1, last_bin_data * 5000.); // last bin data logy
  pad1_template -> GetYaxis() -> SetTitle(title_y);
  pad1_template -> SetStats(0);
  pad1_template -> Draw("histsame");
  
  maphstack[hstack] -> Draw("histsame");
  
  bkg_sum -> SetFillColor(kBlack);
  bkg_sum -> SetFillStyle(3013);
  bkg_sum -> SetMarkerSize(0);;
  bkg_sum -> SetLineWidth(0);
  bkg_sum -> Draw("e2same");

  //data_hist -> Draw("psame");
  data_gr -> SetMarkerSize(1.3);
  data_gr -> SetLineWidth(2);
  data_gr -> Draw("ezpsame");





  maplegend[legend] -> SetFillColor(kWhite);
  //maplegend[legend] -> SetLineColor(kBlack);
  //maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(1001);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> Draw("same");

  gPad->RedrawAxis();

  

  ////////////////////////////////////
  /// Pad 2
  ////////////////////////////////////
  mapcanvas[canvas] -> cd();
  mappad[pad2] = new TPad(pad2, "", 0, 0, 1, 0.25);
  mappad[pad2] -> SetTopMargin( 0.05 );
  mappad[pad2] -> SetBottomMargin( 0.4 );
  mappad[pad2] -> SetLeftMargin( 0.15 );
  mappad[pad2] -> SetRightMargin( 0.04 );
  mappad[pad2] -> Draw();
  mappad[pad2] -> cd();

  //TString name_x = nameofhistogram;
  TString name_x = "Reco. m_{Z'} (GeV)";

  TH1F * pad2_template = new TH1F("", "", N_bin-1, binx);
  pad2_template -> Draw("hist");

  pad2_template -> SetTitle("");
  pad2_template -> SetLineColor(kWhite);
  pad2_template -> GetXaxis() -> SetTitle(name_x);
  pad2_template -> GetXaxis() -> SetTitleSize(0.20);
  //pad2_template -> GetXaxis() -> SetLabelSize(0.125);
  pad2_template -> GetXaxis() -> SetLabelSize(0.15);
  pad2_template -> GetYaxis() -> SetTitle("#frac{Obs.}{Pred.}");
  pad2_template -> GetYaxis() -> SetTitleSize(0.15);
  pad2_template -> GetYaxis() -> SetTitleOffset(0.4);
  //pad2_template -> GetYaxis() -> SetLabelSize(0.10);
  pad2_template -> GetYaxis() -> SetLabelSize(0.125);

  pad2_template -> GetYaxis() -> SetNdivisions(505);
  pad2_template -> GetYaxis() -> SetRangeUser(0.0, 2.0);
  pad2_template -> SetStats(0);
  pad2_template -> Draw("histsame");
  
  TH1F * pad2_bkg_err = new TH1F("", "", N_bin-1, binx);
  TH1F * pad2_ratio = new TH1F("", "", N_bin-1, binx);
  //const int N_points = N_bin - 1;
  //double data_y[N_points], data_x[N_points];
  double ratio_data_y[N_points], ratio_data_ey_high[N_points], ratio_data_ey_low[N_points];
  //double data_ex_high[N_points], data_ex_low[N_points];
  for(int i_bin = 1; i_bin < N_bin; i_bin++){
    double current_data = data_hist -> GetBinContent(i_bin);
    double current_bkg = bkg_sum -> GetBinContent(i_bin);
    //double current_data_ey_high = veyh.at(i_bin-1) / current_data;
    //double current_data_ey_low = veyl.at(i_bin-1) / current_data;
    double current_data_ey_high = veyh.at(i_bin-1) / current_bkg;
    double current_data_ey_low = veyl.at(i_bin-1) / current_bkg;
    double current_bkg_ey = bkg_sum -> GetBinError(i_bin) / current_bkg;
    
    ratio_data_y[i_bin -1] = current_data/current_bkg;
    ratio_data_ey_high[i_bin - 1] = current_data_ey_high;
    ratio_data_ey_low[i_bin - 1] = current_data_ey_low;
    
    pad2_ratio -> SetBinContent(i_bin, current_data/current_bkg);
    pad2_ratio -> SetBinError(i_bin, 0.);
    pad2_bkg_err -> SetBinContent(i_bin, 1.0);
    pad2_bkg_err -> SetBinError(i_bin, current_bkg_ey);
  }
  for(int i_x = 0; i_x < N_points; i_x++){
    double current_ex = (binx[i_x + 1] - binx[i_x]) / 2.0;
    double current_x = binx[i_x] + current_ex;

    data_x[i_x] = current_x;
    data_ex_high[i_x] = current_ex;
    data_ex_low[i_x] =current_ex;
  }
  
  pad2_bkg_err -> SetFillColor(kBlack);
  pad2_bkg_err -> SetMarkerSize(0);
  //pad2_bkg_err -> SetLineColor(kWhite);
  pad2_bkg_err -> SetFillStyle(3013);
  pad2_bkg_err -> Draw("E2same");
  
  pad2_ratio -> Draw("psame");
  
  TGraphAsymmErrors * ratio_data_gr = new TGraphAsymmErrors(N_points, data_x, ratio_data_y, data_ex_low, data_ex_high, ratio_data_ey_low, ratio_data_ey_high);
  ratio_data_gr -> SetMarkerSize(1.3);
  ratio_data_gr -> SetLineWidth(2);
  ratio_data_gr -> Draw("ezpsame");

  maplegend["bottom" + legend] = new TLegend(0.2, 0.85, 0.4, 0.95);
  maplegend["bottom" + legend]->SetBorderSize(0);
  maplegend["bottom" + legend]->SetNColumns(3);
  maplegend["bottom" + legend]->AddEntry(pad2_bkg_err, "Stat. + Syst.", "f");
  maplegend["bottom" + legend]->AddEntry(ratio_data_gr, "Obs./Pred.", "p");
  //maplegend["bottom" + legend]->Draw("same");

  TLine *pad2_line = new TLine(binx[0], 1, binx[N_bin - 1], 1);
  pad2_line -> SetLineStyle(1);
  pad2_line -> SetLineColor(kBlue);
  pad2_line -> Draw("same");

  ratio_data_gr -> Draw("ezpsame");

  gPad->RedrawAxis();

  ////////////////////////////////////
  /// -- Latex
  ////////////////////////////////////
  mapcanvas[canvas] -> cd();
  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  //latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS}");
  latex_Lumi.SetTextSize(0.035);
  latex_Lumi.DrawLatex(0.7, 0.96, "138 fb^{-1} (13 TeV)");

  TLatex latex_channel;
  latex_channel.SetTextSize(0.035);
  if(nameofhistogram.Contains("EMu")){
    latex_channel.DrawLatex(0.20, 0.90, "#font[42]{e#mu}");
  }
  else if(nameofhistogram.Contains("DiEle")){
    latex_channel.DrawLatex(0.20, 0.90, "#font[42]{ee}");
  }
  else if(nameofhistogram.Contains("DiMu")){
    latex_channel.DrawLatex(0.20, 0.90, "#mu#mu");

  }
  else return;

  TLatex latex_region;
  latex_region.SetTextSize(0.035);
  if(nameofhistogram.Contains("0AK8")){
    latex_region.DrawLatex(0.20, 0.85, "0AK8");
  }
  else if(nameofhistogram.Contains("1AK8")){
    latex_region.DrawLatex(0.20, 0.85, "1AK8");
  }
  else if(nameofhistogram.Contains("2AK8")){
    latex_region.DrawLatex(0.20, 0.85, "2AK8");

  }
  else return;

  TLatex latex_fit;
  latex_fit.SetTextSize(0.035);
  if(nameofhistogram.Contains("prefit")){
    latex_fit.DrawLatex(0.20, 0.80, "Pre-fit");
  }
  else if(nameofhistogram.Contains("postfit")){
    latex_fit.DrawLatex(0.20, 0.80, "Post-fit");
  }
  else return;


  if(debug) cout << "9" << endl;

  ////////////////////////////////////
  /// -- Save PDF file
  ////////////////////////////////////
  TString pdfname;
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  pdfname = WORKING_DIR + "/plots/" + TString::Itoa(tag_year,10) + "/";
  if(nameofhistogram.Contains("EMu")){
    pdfname = pdfname + "EMu/";

  }
  else if(nameofhistogram.Contains("DiEle")){
    pdfname = pdfname + "DiEle/";

  }
  else if(nameofhistogram.Contains("DiMu")){
    pdfname = pdfname + "DiMu/";
  }
  else return;
  
  pdfname.Append(nameofhistogram);
  pdfname.Append(".pdf");
  if(debug) cout << "9.1" << endl;

  mapcanvas[canvas] -> SaveAs(pdfname);

  TString out_root_name = WORKING_DIR + "/output/root/" + TString::Itoa(tag_year,10) + "/" + nameofhistogram + ".root";
  TFile *outfile = new TFile(out_root_name, "RECREATE");
  outfile -> cd();
  bkg_sum -> SetName("top_bkg");
  data_gr -> SetName("top_data");
  bkg_sum -> Write();
  data_gr -> Write();
  for(int i_legend = 0; i_legend < 4; i_legend++){
    TString current_hist = plot_name + histograms[4 - i_legend - 1] + which_fit;
    if(mapfunc[current_hist + "sum"]){
      mapfunc[current_hist + "sum"] -> SetName(legend_list[4 - i_legend - 1]);
      mapfunc[current_hist + "sum"] -> Write();
    }
  }
  pad2_bkg_err -> SetName("bottom_unc");
  ratio_data_gr -> SetName("bottom_ratio");
  pad2_bkg_err -> Write();
  ratio_data_gr -> Write();
  outfile -> Close();

  if(debug) cout << "10" << endl;
  
  vx.clear();
  vy.clear(); 
  vexl.clear();
  vexh.clear();
  veyl.clear();
  veyh.clear();
  
}


void open_files(TString plot_name){
  
  //maphist.clear();
  //mapfunc.clear();
  if(debug) cout << "[[open_files]] plot_name : " << plot_name << endl;  
  unsigned int N_bin = map_bin_vector[plot_name].size();
  double current_bins[100];
  for(unsigned int i_bin = 0; i_bin < N_bin; i_bin++){
    current_bins[i_bin] = map_bin_vector[plot_name].at(i_bin);
  }
  
  TString which_fits[2] = {"prefit", "postfit"};
  for(int i_which_fit = 0; i_which_fit < 2; i_which_fit++){
    TString which_fit = which_fits[i_which_fit];
    if(plot_name.Contains("DiEle")){
      openfile(plot_name, which_fit, "HNPairToJJJJ_EE_ZP4000_N1200_WR5000");
    }
    else{
      openfile(plot_name, which_fit, "HNPairToJJJJ_MuMu_ZP4000_N1200_WR5000");
    }   

    make_histogram(plot_name, which_fit, N_bin, current_bins);
  }
}

void open_binning_file(TString filename){

  cout << "[open_binning_file] start to open binngin file : " << filename << endl;
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString binning_file_path=WORKING_DIR+"/script/"+filename;
  ifstream data_file;
  data_file.open(binning_file_path);

  char line[500];
  if(data_file.is_open()){
    while(!data_file.eof()){
      data_file.getline(line, 500);
      if(debug) cout << line << endl;
      TString this_line = line;
      if(this_line.Contains("#")) continue;

      TObjArray *tx = this_line.Tokenize("\t");
      int N_line_part = tx -> GetEntries();
      if(N_line_part < 1) continue;
      if(debug){
        cout << "[[open_binning_file]] histname : " << ((TObjString *)(tx->At(0)))->String() << endl;
        cout << "[[open_binning_file]] binning : " << ((TObjString *)(tx->At(1)))->String() << endl;
      }

      TString current_histname = ((TObjString *)(tx->At(0)))->String();
      if(debug) cout << "[[open_binning_file]] current_histname : " << current_histname << endl;
      map_bin_vector[current_histname].clear();
      TString current_binning = ((TObjString *)(tx->At(1)))->String();
      TObjArray *binning_token = current_binning.Tokenize(", ");
      for(int i_bin_part = 0; i_bin_part < binning_token->GetEntries() - 1; i_bin_part++){
        TString current_bin_value = ((TObjString *)(binning_token->At(i_bin_part)))->String();
        if(i_bin_part == 0){
          current_bin_value.Remove(0,1);
        }
        if(debug) cout << current_bin_value << " | " << endl;
        double current_bin_double = current_bin_value.Atof();
        map_bin_vector[current_histname].push_back(current_bin_double);
      }

      cout << "[[open_binning_file]] current_histname : " << current_histname << endl;

      open_files(current_histname);
      map_bin_vector.clear();
    }
  }
  data_file.close();
}


void overlap_ttbar(TString NAK8){
  
  TString nameofhistogram = "mZp_" + NAK8;
  TString title_y;
  if(nameofhistogram.Contains("N")) title_y = "Number";
  else title_y = "Events";

  TString canvas = nameofhistogram;
  TString legend = nameofhistogram;
  TString func = nameofhistogram;
  TString clone = nameofhistogram;
  TString line = nameofhistogram;
  canvas.Insert(0, "c_");
  legend.Insert(0, "legend_");
  func.Insert(0, "ratio_");
  clone.Insert(0, "h3_");
  line.Insert(0, "l_");

  mapcanvas[canvas] = new TCanvas(canvas,"",800,800);
  canvas_margin(mapcanvas[canvas]);
  gStyle -> SetOptStat(1111);
  mapcanvas[canvas] -> SetLogy();

  maplegend[legend] = new TLegend(0.60, 0.80, 0.96, 0.95);

  TString histname_ee = "mZp_" + NAK8 + "_SR_DiEleTT_powhegprefitsum";
  TString histname_mumu = "mZp_" + NAK8 + "_SR_DiMuTT_powhegprefitsum";
  TString histname_emu = "mZp_" + NAK8 + "_SR_EMuTT_powhegprefitsum";
  
  cout << "histname_emu : " << histname_emu << endl;
  
  mapfunc[func + "template"] = (TH1F*)mapfunc[histname_emu] -> Clone(func + "template");
  mapfunc[func + "template"]  -> SetStats(0);
  double hist_max = mapfunc[func + "template"] -> GetMaximum();
  //mapfunc[func + "template"] -> SetMaximum(hist_max * 1.5);
  mapfunc[func + "template"] -> SetMaximum(hist_max * 100);//logy 
  mapfunc[func + "template"] -> SetMinimum(0.1);
  mapfunc[func + "template"] -> GetXaxis() -> SetTitle(nameofhistogram);
  mapfunc[func + "template"] -> GetYaxis()->SetLabelSize(0.05);;
  mapfunc[func + "template"] -> GetYaxis()->SetTitleSize(0.07);;
  mapfunc[func + "template"] -> GetYaxis()->SetTitleOffset(1.02);;
  mapfunc[func + "template"] -> GetYaxis()->SetTitle("Events/bin");
  mapfunc[func + "template"] -> SetFillColor(0);
  mapfunc[func + "template"] -> Draw("l");

  mapfunc[histname_ee] -> SetLineColor(kBlack);
  mapfunc[histname_mumu] -> SetLineColor(kBlue);
  mapfunc[histname_emu] -> SetLineColor(kRed);
  
  mapfunc[histname_ee] -> SetFillColor(0);
  mapfunc[histname_mumu] -> SetFillColor(0);
  mapfunc[histname_emu] -> SetFillColor(0);

  mapfunc[histname_ee] -> SetLineWidth(2);
  mapfunc[histname_mumu] -> SetLineWidth(2);
  mapfunc[histname_emu] -> SetLineWidth(2);
  
  mapfunc[histname_ee] -> Draw("histsame");
  mapfunc[histname_mumu] -> Draw("histsame");
  mapfunc[histname_emu]-> Draw("histsame");

  maplegend[legend] -> AddEntry(mapfunc[histname_ee], "#font[42]{ee}", "l");
  maplegend[legend] -> AddEntry(mapfunc[histname_mumu], "mumu", "l");
  maplegend[legend] -> AddEntry(mapfunc[histname_emu], "#font[42]{e#mu}", "l");
  
  maplegend[legend] -> SetFillColor(kWhite);
  maplegend[legend] -> SetLineColor(kBlack);
  maplegend[legend] -> SetBorderSize(1);
  maplegend[legend] -> SetFillStyle(1001);
  maplegend[legend] -> SetShadowColor(0);
  maplegend[legend] -> SetEntrySeparation(0.3);
  maplegend[legend] -> Draw("same");

  mapcanvas[canvas] -> cd();

  TLatex latex_CMSPriliminary, latex_Lumi;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  //latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_CMSPriliminary.DrawLatex(0.15, 0.96, "#font[62]{CMS}");
  latex_Lumi.SetTextSize(0.035);
  latex_Lumi.DrawLatex(0.7, 0.96, "138 fb^{-1} (13 TeV)");
  
  TString pdfname;
  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  pdfname = WORKING_DIR + "/plots/Signal/ttbar_comparison_" + NAK8 + ".pdf";
  mapcanvas[canvas] -> SaveAs(pdfname);

  
}




void Plot_Post_fit_FullRun2(int year=2019){
  setTDRStyle();
  tag_year = year;

  TString binning_file = "binning_limit_postfit.txt";
  cout << "binning_file : " << binning_file << endl;

  open_binning_file(binning_file);  

  /*
  open_binning_file("binning_limit_merged_mZp_2AK8_SR_DiEle.txt");
  open_binning_file("binning_limit_merged_mZp_1AK8_SR_DiEle.txt");
  open_binning_file("binning_limit_merged_mZp_0AK8_SR_DiEle.txt");
  open_binning_file("binning_limit_merged_mZp_2AK8_SR_DiMu.txt");
  open_binning_file("binning_limit_merged_mZp_1AK8_SR_DiMu.txt");
  open_binning_file("binning_limit_merged_mZp_0AK8_SR_DiMu.txt");
  */
  /*
  overlap_ttbar("0AK8");
  overlap_ttbar("1AK8");
  overlap_ttbar("2AK8");
  */

}
