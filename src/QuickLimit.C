#include "canvas_margin.h"
#include "mylib.h"

void QuickLimit(int xxx=0){

  double r_value_time = signal_scale * 0.1; // restore r value by weighting to pb from fb
  
  setTDRStyle();

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2(true);
  TH1::AddDirectory(kFALSE);

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString ENV_PLOT_PATH = getenv("PLOT_PATH");
  ENV_PLOT_PATH = ENV_PLOT_PATH + "/limit/";
  TString dirname = "limit_result";
  
  TString channel = "MuMu";
  TString binning = "";
  // --- xxx : 0-3 MuMu, 4-7 EE, 1 5 Inclusive, 2 6 3AK8, 3 7 3AK8_CRs
  if(xxx>3){
    channel = "EE";
  }
  if(xxx== 1 || xxx == 5){
    binning = "Inclusive";
  }
  if(xxx== 2 || xxx == 6){
    binning = "3AK8";
  }
  if(xxx== 3 || xxx == 7){
    binning = "3AK8_CRs";
  }
  // --- xxx additional definitions
  if(xxx == 8){
    channel = "EE";
    binning = "0AK8";
  }
  if(xxx == 9){
    channel = "EE";
    binning = "1AK8";
  }
  if(xxx == 10){
    channel = "EE";
    binning = "2AK8";
  }
  if(xxx == 11){
    channel = "MuMu";
    binning = "0AK8";
  }
  if(xxx == 12){
    channel = "MuMu";
    binning = "1AK8";
  }
  if(xxx == 13){
    channel = "MuMu";
    binning = "2AK8";
  }
  if(xxx == 14){
    channel = "EE";
    binning = "3AK8_2016";
  }
  if(xxx == 15){
    channel = "EE";
    binning = "3AK8_2017";
  }
  if(xxx == 16){
    channel = "EE";
    binning = "3AK8_2018";
  }
  if(xxx == 17){
    channel = "MuMu";
    binning = "3AK8_2016";
  }
  if(xxx == 18){
    channel = "MuMu";
    binning = "3AK8_2017";
  }
  if(xxx == 19){
    channel = "MuMu";
    binning = "3AK8_2018";
  }
  
  
  
  TString ATLAS_result_file;
  if(channel.Contains("MuMu")) ATLAS_result_file = "script/ATLAS_Obs_Limit/ATLAS_mumu_limit.txt";
  if(channel.Contains("EE")) ATLAS_result_file = "script/ATLAS_Obs_Limit/ATLAS_elel_limit.txt";
  
  TString base_filepath = WORKING_DIR+"/output/Result_" + binning + "_" + channel + ".txt";
  TString base_plotpath = ENV_PLOT_PATH+"/" + binning + "/" +  channel;
  cout << "Open " << base_filepath << endl;
  if( !gSystem->mkdir(base_plotpath, kTRUE) ){
    cout
    << "###################################################" << endl
    << "Directoy " << base_plotpath << " is created" << endl
    << "###################################################" << endl
    << endl;
  }
  TFile *out = new TFile(base_plotpath+"/Limits.root","RECREATE");

  TString XSEC_file = "script/cross_section/ZprimeToNN.dat";
  

  //=================================
  //==== 1) xsec vs mN, for each ZP
  //=================================
 
  cout <<" #### 1) xsec vs mN, for each ZP ####" << endl;
  
  vector<int> ZPmasses = {400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000};
  vector<double> nN_crosses;
  
  for(unsigned int i=0; i<ZPmasses.size(); i++){
    
    int mZP = ZPmasses.at(i);
    vector<int> hnmasses = GetHNMassRange(mZP,true);
    
    const int n_Nmass = hnmasses.size();
    
    double x[n_Nmass], y_theory[n_Nmass], y_exp[n_Nmass]; 
    double y_95_down[n_Nmass], y_68_down[n_Nmass], y_68_up[n_Nmass], y_95_up[n_Nmass];
    
    double prev_theory = 10;
    double prev_limit = 1;
    double prev_mass = 1;
    
    cout << "#### " << mZP << " ####" << endl;
    
    vector<double> mN_Cross;
    
    bool HasExclusion = false;
    
    for(unsigned int j=0; j<hnmasses.size(); j++){
      
      int mN = hnmasses.at(j);
      
      x[j] = mN;

      TString this_string = "Zp"+TString::Itoa(mZP,10)+"_HN"+TString::Itoa(mN,10);
      cout << this_string << endl;

      double this_xsec = 1.;

      //==== Find Theory XSEC
      string elline;
      ifstream in(WORKING_DIR+"/"+XSEC_file);
      while(getline(in,elline)){
        std::istringstream is( elline );
	
        // 2000  900 0.0042160
        int a,b;
        double xsec;
        is >> a;
        is >> b;
        is >> xsec;
        if(a==mZP&&mN==b) this_xsec = xsec;
      }
      cout << this_string << "\t" << this_xsec << endl;

      //==== Find Limit Value
      string elline2;
      ifstream in2(base_filepath);
      
      double this_limit(1.), this_limit_95_down(1.), this_limit_68_down(1.), this_limit_68_up(1.), this_limit_95_up(1.);
      
      while(getline(in2,elline2)){
        std::istringstream is( elline2 );

        // 2000  900 0.0042160
        int a, b;
        double c_95_down;
	double c_68_down;
	double c;
        double c_68_up;
        double c_95_up;

        is >> a;
	is >> b;
        is >> c_95_down;
        is >> c_68_down;
	is >> c; 
	is >> c_68_up;
	is >> c_95_up;
	
	
	if(a==mZP && b==mN){
	  cout << "a : " << a << ", b : " << b << ", c : " << c << endl;

	  this_limit_95_down = c_95_down * r_value_time;
	  this_limit_68_down = c_68_down * r_value_time;
	  this_limit = c*r_value_time;
	  this_limit_68_up = c_68_up * r_value_time;
          this_limit_95_up = c_95_up * r_value_time;
	}     
      }

      y_theory[j] = this_xsec;
      
      y_95_down[j] = this_limit_95_down;
      y_68_down[j] = this_limit_68_down;
      y_exp[j]= this_limit;
      y_68_up[j] = this_limit_68_up;
      y_95_up[j] = this_limit_95_up;

      
      if(this_limit < this_xsec) HasExclusion = true;

      if(j>0){
        if(  (prev_theory-prev_limit)*(y_theory[j]-y_exp[j]) < 0 ){

          double N1 = prev_mass;
          double dN = mN-N1;

          double L1 = prev_limit;
          double L2 = y_exp[j];
          double T1 = prev_theory;
          double T2 = y_theory[j];

          double dL = L2-L1;
          double dT = T2-T1;

          double x_meet = (L1-T1)*dN/(dT-dL)+N1;
          mN_Cross.push_back(x_meet);
          cout << x_meet << endl;

        }
      }

      prev_theory = this_xsec;
      prev_limit = this_limit;
      prev_mass = mN;

      cout << this_string << "\t" << this_xsec << "\t" << this_limit << endl;

    } // END loop N mass

    if(mN_Cross.size()==0 && HasExclusion){
      nN_crosses.push_back(hnmasses.at(0));
      nN_crosses.push_back(hnmasses.at(n_Nmass-1));
    }
    else if(mN_Cross.size()==1){
      nN_crosses.push_back(hnmasses.at(0));
      nN_crosses.push_back(mN_Cross.at(0));
    }
    else if(mN_Cross.size()==2){
      nN_crosses.push_back(mN_Cross.at(0));
      nN_crosses.push_back(mN_Cross.at(1));
    }
    else{

    }

    // -- make arrays for each 95% and 68% band                                                                                                                                                                                                                               
    double e_95_down[n_Nmass], e_68_down[n_Nmass], e_68_up[n_Nmass], e_95_up[n_Nmass], e_x_down[n_Nmass], e_x_up[n_Nmass];
    for(int i_band = 0; i_band < n_Nmass; i_band ++){
      e_95_down[i_band] = fabs(y_exp[i_band] - y_95_down[i_band]);
      e_68_down[i_band] = fabs(y_exp[i_band] - y_68_down[i_band]);
      e_68_up[i_band] = fabs(y_exp[i_band] - y_68_up[i_band]);
      e_95_up[i_band] = fabs(y_exp[i_band] - y_95_up[i_band]);
      e_x_down[i_band] = 100.;
      e_x_up[i_band] = 100.;
    }
    
    TGraph *Theory_for_this_mZP = new TGraph(n_Nmass, x, y_theory);
    TGraph *Limit_for_this_mZP = new TGraph(n_Nmass, x, y_exp);
    TGraphAsymmErrors *Limit_for_this_mZP_95_band = new TGraphAsymmErrors(n_Nmass, x, y_exp, e_x_down, e_x_up, e_95_down, e_95_up);
    TGraphAsymmErrors *Limit_for_this_mZP_68_band = new TGraphAsymmErrors(n_Nmass, x, y_exp, e_x_down, e_x_up, e_68_down, e_68_up);
    
    

    cout << "#### Z' = " << mZP << " ####" << endl;
    cout << "---- Theory ----" << endl;
    Theory_for_this_mZP->Print();
    cout << "---- Limit ----" << endl;
    Limit_for_this_mZP->Print();


    TCanvas *c1 = new TCanvas("c1", "", 600, 600);
    canvas_margin(c1);
    c1->cd();

    TH1D *hist_dummy = new TH1D("hist_dummy", "", 5000, 0., 5000.);
    hist_axis(hist_dummy);
    hist_dummy->Draw("hist");
    hist_dummy->GetYaxis()->SetRangeUser(1E-6,100);
    //hist_dummy->GetYaxis()->SetTitle("#sigma(pp #rightarrow Z')BR(Z' #rightarrow NN #rightarrow ee) (pb)");
    hist_dummy->GetYaxis()->SetTitleSize(0.04);
    hist_dummy->GetYaxis()->SetTitleOffset(1.5);
    if(channel == "EE"){
      hist_dummy->GetYaxis()->SetTitle("#sigma(pp #rightarrow Z')BR(Z' #rightarrow NN #rightarrow eejjjj) (pb)");
    }
    if(channel == "MuMu"){
      hist_dummy->GetYaxis()->SetTitle("#sigma(pp #rightarrow Z')BR(Z' #rightarrow NN #rightarrow #mu#mujjjj) (pb)");
    }
    hist_dummy->GetXaxis()->SetRangeUser(hnmasses.at(0),hnmasses.at(hnmasses.size()-1));
    hist_dummy->GetXaxis()->SetTitle("m_{N} (GeV)");
    c1->SetLogy();
    
    Limit_for_this_mZP_95_band -> SetFillColor(kYellow);
    Limit_for_this_mZP_95_band -> SetLineStyle(3);
    if(n_Nmass > 1) Limit_for_this_mZP_95_band -> Draw("3same");
    else Limit_for_this_mZP_95_band -> Draw("e2same");

    Limit_for_this_mZP_68_band -> SetFillColor(kGreen);
    Limit_for_this_mZP_68_band -> SetLineStyle(3);
    if(n_Nmass > 1) Limit_for_this_mZP_68_band -> Draw("3same");
    else Limit_for_this_mZP_68_band -> Draw("e2same");
    
    Theory_for_this_mZP->SetLineColor(kRed);
    Theory_for_this_mZP->SetLineWidth(3);
    Theory_for_this_mZP->SetMarkerStyle(15);
    Theory_for_this_mZP->SetMarkerColor(kRed);
    Theory_for_this_mZP->Draw("lpsame");

    Limit_for_this_mZP->SetLineColor(kBlack);
    Limit_for_this_mZP->SetLineStyle(3);
    Limit_for_this_mZP->SetLineWidth(3);
    Limit_for_this_mZP->SetMarkerStyle(15);
    Limit_for_this_mZP->SetMarkerColor(kBlack);
    Limit_for_this_mZP->Draw("lpsame");

    TLegend *lg = new TLegend(0.7, 0.7, 0.95, 0.95);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(Theory_for_this_mZP, "Theory", "lp");
    lg->AddEntry(Limit_for_this_mZP, "Expected", "lp");
    lg->AddEntry(Limit_for_this_mZP_68_band, "Expected #pm 1 #sigma", "lf");
    lg->AddEntry(Limit_for_this_mZP_95_band, "Expected #pm 2 #sigma", "lf");


    lg->Draw();

    TLatex latex_zpmass;
    latex_zpmass.SetTextSize(0.05);
    latex_zpmass.SetNDC();
    latex_zpmass.DrawLatex(0.2, 0.9,"m_{Z'} = "+TString::Itoa(mZP,10)+" GeV");

    c1->SaveAs(base_plotpath+"/Limit_Zp"+TString::Itoa(mZP,10)+".pdf");
    //c1->SaveAs(base_plotpath+"/Limit_Zp"+TString::Itoa(mZP,10)+".png");
    
    out->cd();
    Theory_for_this_mZP->SetName("Zp"+TString::Itoa(mZP,10)+"_Theory");
    Limit_for_this_mZP->SetName("Zp"+TString::Itoa(mZP,10)+"_Limit");
    Theory_for_this_mZP->Write();
    Limit_for_this_mZP->Write();

    c1->Close();

  }

  //==================================
  //==== 2) xsec vs mZP, for each nN
  //==================================

  cout <<" #### 2) xsec vs mZP, for each N ####" << endl;

  vector<int> Nmasses = {100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400};
  vector<double> mZP_crosses;

  for(unsigned int i=0; i<Nmasses.size(); i++){

    int mN = Nmasses.at(i);
    vector<int> zpmasses = GetZPMassRange(mN);

    const int n_ZPmass = zpmasses.size();

    double x[n_ZPmass], y_theory[n_ZPmass], y_exp[n_ZPmass];
    double y_95_down[n_ZPmass], y_68_down[n_ZPmass], y_68_up[n_ZPmass], y_95_up[n_ZPmass];

    
    double prev_theory = 10;
    double prev_limit = 1;
    double prev_mass = 1;

    cout << "#### " << mN << " ####" << endl;

    vector<double> mZP_Cross;

    bool HasExclusion = false;

    for(unsigned int j=0; j<zpmasses.size(); j++){

      int mZP = zpmasses.at(j);

      x[j] = mZP;

      TString this_string = "Zp"+TString::Itoa(mZP,10)+"_HN"+TString::Itoa(mN,10);
      //cout << this_string << endl;

      double this_xsec = 1.;

      //==== Find Theory XSEC
      string elline;
      ifstream in(WORKING_DIR+"/"+XSEC_file);
      while(getline(in,elline)){
        std::istringstream is( elline );

        // 2000  900 0.0042160 
        int a,b;
        double xsec;
        is >> a;
        is >> b;
        is >> xsec;
        if(a==mZP&&mN==b) this_xsec = xsec;

      }
      //cout << this_string << "\t" << this_xsec << endl;

      //==== Find Limit Value
      string elline2;
      ifstream in2(base_filepath);
      double this_limit(1.), this_limit_95_down(1.), this_limit_68_down(1.), this_limit_68_up(1.), this_limit_95_up(1.);
      while(getline(in2,elline2)){
        std::istringstream is( elline2 );
	
        // 2000  900 0.0042160
        int a, b;
        double c_95_down;
        double c_68_down;
        double c;
        double c_68_up;
        double c_95_up;

	is >> a;
	is >> b;
        is >> c_95_down;
	is >> c_68_down;
        is >> c;
        is >> c_68_up;
        is >> c_95_up;

        if(a==mZP && b==mN){
          this_limit_95_down = c_95_down * r_value_time;
          this_limit_68_down = c_68_down * r_value_time;
          this_limit = c*r_value_time;
          this_limit_68_up = c_68_up * r_value_time;
          this_limit_95_up = c_95_up * r_value_time;
	}
      }

      y_theory[j] = this_xsec;

      y_95_down[j] = this_limit_95_down;
      y_68_down[j] = this_limit_68_down;
      y_exp[j]= this_limit;
      y_68_up[j] = this_limit_68_up;
      y_95_up[j] = this_limit_95_up;


      if(this_limit < this_xsec) HasExclusion = true;

      if(j>0){
        if(  (prev_theory-prev_limit)*(y_theory[j]-y_exp[j]) < 0 ){

          double ZP1 = prev_mass;
          double dZP = mZP-ZP1;

          double L1 = prev_limit;
          double L2 = y_exp[j];
          double T1 = prev_theory;
          double T2 = y_theory[j];

/*
          cout << "prev_mass = " << prev_mass << endl;
          cout << "mZP = " << mZP << endl;
          cout << "dZP = " << dZP << endl;
          cout << "L1 = " << L1 << endl;
          cout << "L2 = " << L2 << endl;
          cout << "T1 = " << T1 << endl;
          cout << "T2 = " << T2 << endl;
*/

          double dL = L2-L1;
          double dT = T2-T1;

          double x_meet = (L1-T1)*dZP/(dT-dL)+ZP1;
          mZP_Cross.push_back(x_meet);
          cout << x_meet << endl;

        }
      }

      prev_theory = this_xsec;
      prev_limit = this_limit;
      prev_mass = mZP;

      //cout << this_string << "\t" << this_xsec << "\t" << this_limit << endl;

    } // END loop N mass

    if(mZP_Cross.size()==0 && HasExclusion){
      mZP_crosses.push_back(zpmasses.at(0));
      mZP_crosses.push_back(zpmasses.at(n_ZPmass-1));
    }
    else if(mZP_Cross.size()==1){
      mZP_crosses.push_back(zpmasses.at(0));
      mZP_crosses.push_back(mZP_Cross.at(0));
    }
    else if(mZP_Cross.size()==2){
      mZP_crosses.push_back(mZP_Cross.at(0));
      mZP_crosses.push_back(mZP_Cross.at(1));
    }
    else{

    }


    double e_95_down[n_ZPmass], e_68_down[n_ZPmass], e_68_up[n_ZPmass], e_95_up[n_ZPmass], e_x_down[n_ZPmass], e_x_up[n_ZPmass];
    for(int i_band = 0; i_band < n_ZPmass; i_band ++){
      e_95_down[i_band] = fabs(y_exp[i_band] - y_95_down[i_band]);
      e_68_down[i_band] = fabs(y_exp[i_band] - y_68_down[i_band]);
      e_68_up[i_band] = fabs(y_exp[i_band] - y_68_up[i_band]);
      e_95_up[i_band] = fabs(y_exp[i_band] - y_95_up[i_band]);
      e_x_down[i_band] = 100.;
      e_x_up[i_band] = 100.;
    }

    TGraph *Theory_for_this_mN = new TGraph(n_ZPmass,x,y_theory);
    TGraph *Limit_for_this_mN = new TGraph(n_ZPmass,x,y_exp);
    TGraphAsymmErrors *Limit_for_this_mZP_95_band = new TGraphAsymmErrors(n_ZPmass, x, y_exp, e_x_down, e_x_up, e_95_down, e_95_up);
    TGraphAsymmErrors *Limit_for_this_mZP_68_band = new TGraphAsymmErrors(n_ZPmass, x, y_exp, e_x_down, e_x_up, e_68_down, e_68_up);
/*
    cout << "#### Z' = " << mN << " ####" << endl;
    cout << "---- Theory ----" << endl;
    Theory_for_this_mN->Print();
    cout << "---- Limit ----" << endl;
    Limit_for_this_mN->Print();
*/

    TCanvas *c1 = new TCanvas("c1", "", 600, 600);
    canvas_margin(c1);
    c1->cd();

    TH1D *hist_dummy = new TH1D("hist_dummy", "", 5000, 0., 5000.);
    hist_axis(hist_dummy);
    hist_dummy->Draw("hist");
    hist_dummy->GetYaxis()->SetRangeUser(1E-6,100);
    hist_dummy->GetYaxis()->SetTitle("#sigma (pb)");
    hist_dummy->GetYaxis()->SetTitleSize(0.04);
    hist_dummy->GetYaxis()->SetTitleOffset(1.5);
    if(channel == "EE"){
      hist_dummy->GetYaxis()->SetTitle("#sigma(pp #rightarrow Z')BR(Z' #rightarrow NN #rightarrow eejjjj) (pb)");
    }
    if(channel == "MuMu"){
      hist_dummy->GetYaxis()->SetTitle("#sigma(pp #rightarrow Z')BR(Z' #rightarrow NN #rightarrow #mu#mujjjj) (pb)");
    }

    //hist_dummy->GetXaxis()->SetRangeUser(zpmasses.at(0),zpmasses.at(zpmasses.size()-1));
    hist_dummy->GetXaxis()->SetRangeUser(400., 4400.); // -- FIXME
    hist_dummy->GetXaxis()->SetTitle("m_{Z'} (GeV)");
    c1->SetLogy();

    Limit_for_this_mZP_95_band -> SetFillColor(kYellow);
    Limit_for_this_mZP_95_band -> SetLineStyle(3);
    if(n_ZPmass > 1) Limit_for_this_mZP_95_band -> Draw("3same");
    else Limit_for_this_mZP_95_band -> Draw("e2same");

    Limit_for_this_mZP_68_band -> SetFillColor(kGreen);
    Limit_for_this_mZP_68_band -> SetLineStyle(3);
    if(n_ZPmass > 1) Limit_for_this_mZP_68_band -> Draw("3same");
    else Limit_for_this_mZP_68_band -> Draw("e2same");



    Theory_for_this_mN->SetLineColor(kRed);
    Theory_for_this_mN->SetLineWidth(3);
    Theory_for_this_mN->SetMarkerStyle(15);
    Theory_for_this_mN->SetMarkerColor(kRed);
    Theory_for_this_mN->Draw("lpsame");

    Limit_for_this_mN->SetLineColor(kBlack);
    Limit_for_this_mN->SetLineStyle(3);
    Limit_for_this_mN->SetLineWidth(3);
    Limit_for_this_mN->SetMarkerStyle(15);
    Limit_for_this_mN->SetMarkerColor(kBlack);
    Limit_for_this_mN->Draw("lpsame");

    TLegend *lg = new TLegend(0.7, 0.7, 0.95, 0.95);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(Theory_for_this_mN, "Theory", "lp");
    lg->AddEntry(Limit_for_this_mN, "Expected", "lp");
    lg->AddEntry(Limit_for_this_mZP_68_band, "Expected #pm 1 #sigma", "lf");
    lg->AddEntry(Limit_for_this_mZP_95_band, "Expected #pm 2 #sigma", "lf");
    lg->Draw();

    TLatex latex_zpmass;
    latex_zpmass.SetTextSize(0.035);
    latex_zpmass.SetNDC();
    latex_zpmass.DrawLatex(0.25, 0.9,"m_{N} = "+TString::Itoa(mN,10)+" GeV");

    TLatex latex_channel;
    latex_channel.SetTextSize(0.035);
    latex_channel.SetNDC();
    if(channel == "EE"){
      latex_channel.DrawLatex(0.25, 0.85,"Combined ee channel");
    }
    if(channel == "MuMu"){
      latex_channel.DrawLatex(0.25, 0.85,"Combined #mu#mu channel");
    }

    c1->SaveAs(base_plotpath+"/Limit_N"+TString::Itoa(mN,10)+".pdf");
    //c1->SaveAs(base_plotpath+"/Limit_N"+TString::Itoa(mN,10)+".png");

    out->cd();
    Theory_for_this_mN->SetName("Zp"+TString::Itoa(mN,10)+"_Theory");
    Limit_for_this_mN->SetName("Zp"+TString::Itoa(mN,10)+"_Limit");
    Theory_for_this_mN->Write();
    Limit_for_this_mN->Write();

    c1->Close();

  }

  //========================================
  //==== 3) Save Limits on Map
  //========================================
  map<TString, double> limit_exp;
  map<TString, double> limit_exp_95_down;
  map<TString, double> limit_exp_68_down;
  map<TString, double> limit_exp_68_up;
  map<TString, double> limit_exp_95_up;
  map<TString, double> limit_theory;

  string elline;
  ifstream in(WORKING_DIR+"/"+XSEC_file);
  while(getline(in,elline)){
    std::istringstream is( elline );
    int a,b;
    double xsec;
    is >> a;
    is >> b;
    is >> xsec;
    TString mZP = TString::Itoa(a, 10);
    TString mN = TString::Itoa(b, 10);
    limit_theory[mZP + "_" + mN] = xsec;
  }

  string elline2;
  ifstream in2(base_filepath);
  while(getline(in2,elline2)){
    std::istringstream is( elline2 );
    TString a, b;
    double c_95_down;
    double c_68_down;
    double c;
    double c_68_up;
    double c_95_up;

    is >> a;
    is >> b;
    is >> c_95_down;
    is >> c_68_down;
    is >> c;
    is >> c_68_up;
    is >> c_95_up;

    TString mZp = a;
    TString mN = b;

    if(a.Atoi() > 0 && b.Atoi() > 0){
      limit_exp_95_down[mZp + "_" + mN] = c_95_down * r_value_time;
      limit_exp_68_down[mZp + "_" + mN] = c_68_down * r_value_time;
      limit_exp[mZp + "_" + mN] = c*r_value_time;
      limit_exp_68_up[mZp + "_" + mN] = c_68_up * r_value_time;
      limit_exp_95_up[mZp + "_" + mN] = c_95_up * r_value_time;

      double mZp_double = mZp.Atof() + 10.;
      double mN_double = mN.Atof() + 10.;
      //exp_limit_2D -> Fill(mZp_double, mN_double, c*r_value_time);
    }
  }



  //==================================
  //==== 4) Cross section limit for  mN = 0.25 mZ'
  //==================================
  TCanvas *c_quarter = new TCanvas("c_quarter", "", 600, 600);
  canvas_margin(c_quarter);
  c_quarter->cd();
  
  double a_quarter = 0.25; // y = 0.25 x for mN = 0.25 mZ'
  int N_mZ_quarter = 21; // 21 = to 4400
  double current_mN_quarter[N_mZ_quarter];
  int mZp_array_quarter[] = {400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400};
  vector<double> Cross_95_down;
  vector<double> Cross_68_down;
  vector<double> Cross;
  vector<double> Cross_68_up;
  vector<double> Cross_95_up;
  
  double y_theory_quarter[N_mZ_quarter];
  double y_exp_quarter[N_mZ_quarter], y_exp_quarter_95_down[N_mZ_quarter], y_exp_quarter_68_down[N_mZ_quarter], y_exp_quarter_68_up[N_mZ_quarter], y_exp_quarter_95_up[N_mZ_quarter];
  double x_quarter[N_mZ_quarter];
  
  //Get mN values needed
  for(int i_N = 0; i_N < N_mZ_quarter; i_N++){
    current_mN_quarter[i_N] = a_quarter * (mZp_array_quarter[i_N] + 0.); // mN = 0.25 mZ'
    x_quarter[i_N] = mZp_array_quarter[i_N];
  }
  
  for(int i_N = 0; i_N < N_mZ_quarter; i_N++){
    double mN = current_mN_quarter[i_N];
    int mN_low;
    int mN_high;
    mN_low = GetLowMN(mN);
    mN_high = mN_low + 100;
    TString current_mZp_str = TString::Itoa(mZp_array_quarter[i_N], 10);
    TString current_mN_low_str = TString::Itoa(mN_low, 10);
    TString current_mN_high_str = TString::Itoa(mN_high, 10);
    
    //Linear extrapolate theory cross sections                                                                                                                                                                                                                                
    double mN_low_xsec = limit_theory[current_mZp_str + "_" +current_mN_low_str];
    double mN_high_xsec = limit_theory[current_mZp_str + "_" +current_mN_high_str];
    y_theory_quarter[i_N] = mN_low_xsec - 0.005 * (mN_low_xsec - mN_high_xsec) * (mN - mN_low + 0.);
   
    // -- Linear extrapolate cross section limits
    mN_low_xsec = limit_exp[current_mZp_str + "_" +current_mN_low_str];
    mN_high_xsec = limit_exp[current_mZp_str + "_" +current_mN_high_str];
    y_exp_quarter[i_N] = mN_low_xsec - 0.005 * (mN_low_xsec - mN_high_xsec) * (mN - mN_low + 0.);

    mN_low_xsec = limit_exp_95_down[current_mZp_str + "_" +current_mN_low_str];
    mN_high_xsec = limit_exp_95_down[current_mZp_str + "_" +current_mN_high_str];
    y_exp_quarter_95_down[i_N] = mN_low_xsec - 0.005 * (mN_low_xsec - mN_high_xsec) * (mN - mN_low + 0.);

    mN_low_xsec = limit_exp_68_down[current_mZp_str + "_" +current_mN_low_str];
    mN_high_xsec = limit_exp_68_down[current_mZp_str + "_" +current_mN_high_str];
    y_exp_quarter_68_down[i_N] = mN_low_xsec - 0.005 * (mN_low_xsec - mN_high_xsec) * (mN - mN_low + 0.);

    mN_low_xsec = limit_exp_68_up[current_mZp_str + "_" +current_mN_low_str];
    mN_high_xsec = limit_exp_68_up[current_mZp_str + "_" +current_mN_high_str];
    y_exp_quarter_68_up[i_N] = mN_low_xsec - 0.005 * (mN_low_xsec - mN_high_xsec) * (mN - mN_low + 0.);

    mN_low_xsec = limit_exp_95_up[current_mZp_str + "_" +current_mN_low_str];
    mN_high_xsec = limit_exp_95_up[current_mZp_str + "_" +current_mN_high_str];
    y_exp_quarter_95_up[i_N] = mN_low_xsec - 0.005 * (mN_low_xsec - mN_high_xsec) * (mN - mN_low + 0.);
    
  }

  double e_95_down[N_mZ_quarter], e_68_down[N_mZ_quarter], e_68_up[N_mZ_quarter], e_95_up[N_mZ_quarter], e_x_down[N_mZ_quarter], e_x_up[N_mZ_quarter];
  for(int i_band = 0; i_band < N_mZ_quarter; i_band ++){
    e_95_down[i_band] = fabs(y_exp_quarter[i_band] - y_exp_quarter_95_down[i_band]);
    e_68_down[i_band] = fabs(y_exp_quarter[i_band] - y_exp_quarter_68_down[i_band]);
    e_68_up[i_band] = fabs(y_exp_quarter[i_band] - y_exp_quarter_68_up[i_band]);
    e_95_up[i_band] = fabs(y_exp_quarter[i_band] - y_exp_quarter_95_up[i_band]);
    e_x_down[i_band] = 100.;
    e_x_up[i_band] = 100.;
  }

  TGraph *Theory_quarter = new TGraph(N_mZ_quarter, x_quarter , y_theory_quarter);
  TGraph *Limit_quarter = new TGraph(N_mZ_quarter, x_quarter , y_exp_quarter);
  TGraphAsymmErrors *Limit_quarter_95_band = new TGraphAsymmErrors(N_mZ_quarter, x_quarter, y_exp_quarter, e_x_down, e_x_up, e_95_down, e_95_up);
  TGraphAsymmErrors *Limit_quarter_68_band = new TGraphAsymmErrors(N_mZ_quarter, x_quarter, y_exp_quarter, e_x_down, e_x_up, e_68_down, e_68_up);
  
  TH1D *hist_dummy = new TH1D("hist_dummy", "", 4000, 400., 4400.);
  hist_axis(hist_dummy);
  hist_dummy->Draw("hist");
  hist_dummy->GetYaxis()->SetRangeUser(1E-6,100);
  hist_dummy->GetYaxis()->SetTitle("#sigma (pb)");
  hist_dummy->GetYaxis()->SetTitleSize(0.04);
  hist_dummy->GetYaxis()->SetTitleOffset(1.5);
  if(channel == "EE"){
    hist_dummy->GetYaxis()->SetTitle("#sigma(pp #rightarrow Z')BR(Z' #rightarrow NN #rightarrow eejjjj) (pb)");
  }
  if(channel == "MuMu"){
    hist_dummy->GetYaxis()->SetTitle("#sigma(pp #rightarrow Z')BR(Z' #rightarrow NN #rightarrow #mu#mujjjj) (pb)");
  }
  hist_dummy->GetXaxis()->SetTitle("m_{Z'} (GeV)");
  c_quarter->SetLogy();

  Limit_quarter_95_band -> SetFillColor(kYellow);
  Limit_quarter_95_band -> SetLineStyle(3);
  if(N_mZ_quarter > 1) Limit_quarter_95_band -> Draw("3same");
  else Limit_quarter_95_band -> Draw("e2same");

  Limit_quarter_68_band -> SetFillColor(kGreen);
  Limit_quarter_68_band -> SetLineStyle(3);
  if(N_mZ_quarter > 1) Limit_quarter_68_band -> Draw("3same");
  else Limit_quarter_68_band -> Draw("e2same");

  Theory_quarter->SetLineColor(kRed);
  Theory_quarter->SetLineWidth(3);
  Theory_quarter->SetMarkerStyle(15);
  Theory_quarter->SetMarkerColor(kRed);
  Theory_quarter->Draw("lpsame");

  Limit_quarter->SetLineColor(kBlack);
  Limit_quarter->SetLineStyle(3);
  Limit_quarter->SetLineWidth(3);
  Limit_quarter->SetMarkerStyle(15);
  Limit_quarter->SetMarkerColor(kBlack);
  Limit_quarter->Draw("lpsame");

  TLegend *lg_quarter = new TLegend(0.7, 0.7, 0.95, 0.95);
  lg_quarter->SetBorderSize(0);
  lg_quarter->SetFillStyle(0);
  lg_quarter->AddEntry(Theory_quarter, "Theory", "lp");
  lg_quarter->AddEntry(Limit_quarter, "Expected", "lp");
  lg_quarter->AddEntry(Limit_quarter_68_band, "Expected #pm 1 #sigma", "lf");
  lg_quarter->AddEntry(Limit_quarter_95_band, "Expected #pm 2 #sigma", "lf");
  lg_quarter->Draw();

  TLatex latex_quarter;
  latex_quarter.SetTextSize(0.035);
  latex_quarter.SetNDC();
  latex_quarter.DrawLatex(0.25, 0.9,"m_{N} = m_{Z'}/4");

  TLatex latex_channel_quarter;
  latex_channel_quarter.SetTextSize(0.035);
  latex_channel_quarter.SetNDC();
  if(channel == "EE"){
    latex_channel_quarter.DrawLatex(0.25, 0.85,"Combined ee channel");
  }
  if(channel == "MuMu"){
    latex_channel_quarter.DrawLatex(0.25, 0.85,"Combined #mu#mu channel");
  }
  
  
  c_quarter->SaveAs(base_plotpath+"/Limit_quarter.pdf");
  
  c_quarter->Close();

  //=============================
  //==== Okay, Now Draw Contour
  //=============================
  
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
  
  TH2D *exp_limit_2D = new TH2D("", "", 12 ,400., 5200., 16, 100., 3300); //for 2D xsec limit plot
  
  //========================================
  //==== 2) For loop over y = ax + 100(1 - 4a)
  //========================================
  
  int N_a = 50; // how many different a's
  int N_mZp = 24; // how many Zp mass points
  int mZp_array[] = {400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000};
  mZP_crosses.clear();
  vector<double> mZP_crosses_95_down, mZP_crosses_68_down, mZP_crosses_68_up, mZP_crosses_95_up;
  vector<double> mN_crosses_95_down, mN_crosses_68_down, mN_crosses, mN_crosses_68_up, mN_crosses_95_up;
  mZP_crosses_95_down.push_back(400.);
  mZP_crosses_68_down.push_back(400.);
  mZP_crosses.push_back(400.);
  mZP_crosses_68_up.push_back(400.);
  mZP_crosses_95_up.push_back(400.);

  mN_crosses_95_down.push_back(100.);
  mN_crosses_68_down.push_back(100.);
  mN_crosses.push_back(100.);
  mN_crosses_68_up.push_back(100.);
  mN_crosses_95_up.push_back(100.);
  
  
  
  for(int i_a = 0; i_a < N_a; i_a++){
    
    double current_a = 0.5 * (i_a + 0.) / (N_a + 0.);
    cout << "current_a : " << current_a << endl;
    double current_mN[N_mZp];
    vector<double> mZP_Cross_95_down;
    vector<double> mZP_Cross_68_down;
    vector<double> mZP_Cross;
    vector<double> mZP_Cross_68_up;
    vector<double> mZP_Cross_95_up;

    vector<double> mN_Cross_95_down;
    vector<double> mN_Cross_68_down;
    vector<double> mN_Cross;
    vector<double> mN_Cross_68_up;
    vector<double> mN_Cross_95_up;
    
    //Get mN values needed
    for(int i_N = 0; i_N < N_mZp; i_N++){
      current_mN[i_N] = current_a * (mZp_array[i_N] + 0.) + 100.0 * (1.0 - 4.0 * current_a); // mN = a*mZp + 100(1 - 4a)
    }
    
    //Get linearly extrapolated limits
    double previous_diff[5] = {999., 999., 999., 999., 999.};
    bool HasExclusion[5] = {false, false, false, false, false};
    
    int n_ZPmass = N_mZp;
    double x[n_ZPmass], y_theory[n_ZPmass], y_exp[n_ZPmass];
    double y_exp_95_down[n_ZPmass], y_exp_68_down[n_ZPmass], y_exp_68_up[n_ZPmass], y_exp_95_up[n_ZPmass];
    
    for(int i_N = 0; i_N < N_mZp; i_N++){
      double mN = current_mN[i_N];
      int mN_low;
      int mN_high;
      mN_low = GetLowMN(mN);
      mN_high = mN_low + 100;
      TString current_mZp_str = TString::Itoa(mZp_array[i_N], 10);
      TString current_mN_low_str = TString::Itoa(mN_low, 10);
      TString current_mN_high_str = TString::Itoa(mN_high, 10);
      
      //Linear extrapolate theory cross sections
      double mN_low_xsec = limit_theory[current_mZp_str + "_" +current_mN_low_str];
      double mN_high_xsec = limit_theory[current_mZp_str + "_" +current_mN_high_str];
      
      //y_theory[i_N] = ((mN - mN_low) * mN_low_xsec + (mN_high - mN) * mN_high) / 200.;
      y_theory[i_N] = mN_low_xsec - 0.005 * (mN_low_xsec - mN_high_xsec) * (mN - mN_low + 0.); 
      
      
      //Linear extrapolate cross section limits
      mN_low_xsec = limit_exp[current_mZp_str + "_" +current_mN_low_str];
      mN_high_xsec = limit_exp[current_mZp_str + "_" +current_mN_high_str];
      //y_exp[i_N] = ((mN - mN_low) * mN_low_xsec + (mN_high - mN) * mN_high) / 200.;
      y_exp[i_N] = mN_low_xsec - 0.005 * (mN_low_xsec - mN_high_xsec) * (mN - mN_low + 0.);
      
      mN_low_xsec = limit_exp_95_down[current_mZp_str + "_" +current_mN_low_str];
      mN_high_xsec = limit_exp_95_down[current_mZp_str + "_" +current_mN_high_str];
      //y_exp_95_down[i_N] = ((mN - mN_low) * mN_low_xsec + (mN_high - mN) * mN_high) / 200.;
      y_exp_95_down[i_N] = mN_low_xsec - 0.005 * (mN_low_xsec - mN_high_xsec) * (mN - mN_low + 0.);

      mN_low_xsec = limit_exp_68_down[current_mZp_str + "_" +current_mN_low_str];
      mN_high_xsec = limit_exp_68_down[current_mZp_str + "_" +current_mN_high_str];
      //y_exp_68_down[i_N] = ((mN - mN_low) * mN_low_xsec + (mN_high - mN) * mN_high) / 200.;
      y_exp_68_down[i_N] = mN_low_xsec - 0.005 * (mN_low_xsec - mN_high_xsec) * (mN - mN_low + 0.);

      mN_low_xsec = limit_exp_68_up[current_mZp_str + "_" +current_mN_low_str];
      mN_high_xsec = limit_exp_68_up[current_mZp_str + "_" +current_mN_high_str];
      //y_exp_68_up[i_N] = ((mN - mN_low) * mN_low_xsec + (mN_high - mN) * mN_high) / 200.;
      y_exp_68_up[i_N] = mN_low_xsec - 0.005 * (mN_low_xsec - mN_high_xsec) * (mN - mN_low + 0.);

      mN_low_xsec = limit_exp_95_up[current_mZp_str + "_" +current_mN_low_str];
      mN_high_xsec = limit_exp_95_up[current_mZp_str + "_" +current_mN_high_str];
      //y_exp_95_up[i_N] = ((mN - mN_low) * mN_low_xsec + (mN_high - mN) * mN_high) / 200.;
      y_exp_95_up[i_N] = mN_low_xsec - 0.005 * (mN_low_xsec - mN_high_xsec) * (mN - mN_low + 0.);

      
      if(y_exp_95_down[i_N] < y_theory[i_N]) HasExclusion[0] = true;
      if(y_exp_68_down[i_N] < y_theory[i_N]) HasExclusion[1] = true;
      if(y_exp[i_N] < y_theory[i_N]) HasExclusion[2] = true;
      if(y_exp_68_up[i_N] < y_theory[i_N]) HasExclusion[3] = true;
      if(y_exp_95_up[i_N] < y_theory[i_N]) HasExclusion[4] = true;
      
      double current_diff[5];
      current_diff[0] = y_theory[i_N] - y_exp_95_down[i_N];
      current_diff[1] = y_theory[i_N] - y_exp_68_down[i_N];
      current_diff[2] = y_theory[i_N] - y_exp[i_N];
      current_diff[3] = y_theory[i_N] - y_exp_68_up[i_N];
      current_diff[4] = y_theory[i_N] - y_exp_95_up[i_N];

      
      if(i_N>0){
	double ZP1 = mZp_array[i_N - 1] + 0.;
	double dZP = (mZp_array[i_N] + 0.) - (mZp_array[i_N - 1] + 0.);
	
	double T1 = y_theory[i_N - 1];
	double T2 = y_theory[i_N];
	double dT = T2-T1;
	
	//95 down
	if(current_diff[0] * previous_diff[0] < 0){
	  double L1 = y_exp_95_down[i_N - 1];
          double L2 = y_exp_95_down[i_N];
	  double dL = L2-L1;
          
          double x_meet = (L1-T1)*dZP/(dT-dL)+ZP1;
          mZP_Cross_95_down.push_back(x_meet);
	  mN_Cross_95_down.push_back(current_a * (x_meet + 0.) + 100.0 * (1.0 - 4.0 * current_a));
	}

	//68 down
	if(current_diff[1] * previous_diff[1] < 0){
          double L1 = y_exp_68_down[i_N - 1];
          double L2 = y_exp_68_down[i_N];
          double dL = L2-L1;

          double x_meet = (L1-T1)*dZP/(dT-dL)+ZP1;
          mZP_Cross_68_down.push_back(x_meet);
	  mN_Cross_68_down.push_back(current_a * (x_meet + 0.) + 100.0 * (1.0 - 4.0 * current_a));

	}

	//medium
	if(current_diff[2] * previous_diff[2] < 0){
          double L1 = y_exp[i_N - 1];
          double L2 = y_exp[i_N];
          double dL = L2-L1;

          double x_meet = (L1-T1)*dZP/(dT-dL)+ZP1;
          mZP_Cross.push_back(x_meet);
	  mN_Cross.push_back(current_a * (x_meet + 0.) + 100.0 * (1.0 - 4.0 * current_a));
	}
	
	//68 up
	if(current_diff[3] * previous_diff[3] < 0){
          double L1 = y_exp_68_up[i_N - 1];
          double L2 = y_exp_68_up[i_N];
          double dL = L2-L1;

          double x_meet = (L1-T1)*dZP/(dT-dL)+ZP1;
          mZP_Cross_68_up.push_back(x_meet);
	  mN_Cross_68_up.push_back(current_a * (x_meet + 0.) + 100.0 * (1.0 - 4.0 * current_a));
	}
	
	//95 up
	if(current_diff[4] * previous_diff[4] < 0){
          double L1 = y_exp_95_up[i_N - 1];
          double L2 = y_exp_95_up[i_N];
          double dL = L2-L1;

          double x_meet = (L1-T1)*dZP/(dT-dL)+ZP1;
          mZP_Cross_95_up.push_back(x_meet);
	  mN_Cross_95_up.push_back(current_a * (x_meet + 0.) + 100.0 * (1.0 - 4.0 * current_a));
	}
      }
      
      previous_diff[0] = y_theory[i_N] - y_exp_95_down[i_N];
      previous_diff[1] = y_theory[i_N] - y_exp_68_down[i_N];
      previous_diff[2] = y_theory[i_N] - y_exp[i_N];
      previous_diff[3] = y_theory[i_N] - y_exp_68_up[i_N];
      previous_diff[4] = y_theory[i_N] - y_exp_95_up[i_N];

    }//end for loop over mZp's (i_N)
    
    
    //Fill in crossing points
    if(mZP_Cross_95_down.size()>0){
      mZP_crosses_95_down.push_back(mZP_Cross_95_down.at(0));
      mN_crosses_95_down.push_back(mN_Cross_95_down.at(0));
    }
    if(mZP_Cross_68_down.size()>0){
      mZP_crosses_68_down.push_back(mZP_Cross_68_down.at(0));
      mN_crosses_68_down.push_back(mN_Cross_68_down.at(0));
    }
    if(mZP_Cross.size()>0){
      mZP_crosses.push_back(mZP_Cross.at(0));
      mN_crosses.push_back(mN_Cross.at(0));
    }
    if(mZP_Cross_68_up.size()>0){
      mZP_crosses_68_up.push_back(mZP_Cross_68_up.at(0));
      mN_crosses_68_up.push_back(mN_Cross_68_up.at(0));
    }
    if(mZP_Cross_95_up.size()>0){
      mZP_crosses_95_up.push_back(mZP_Cross_95_up.at(0));
      mN_crosses_95_up.push_back(mN_Cross_95_up.at(0));
    }
    
  }//end for loop over a's

  ////////////////////////////////////////////////////
  // === one more time on the line of y = x/2 - 50
  ///////////////////////////////////////////////////
  vector<double> theory_cross_50;
  vector<double> limit_cross_50;
  vector<double> limit_cross_95_down_50;
  vector<double> limit_cross_68_down_50;
  vector<double> limit_cross_68_up_50;
  vector<double> limit_cross_95_up_50;
 
  vector<double> mZP_Cross_95_down;
  vector<double> mZP_Cross_68_down;
  vector<double> mZP_Cross;
  vector<double> mZP_Cross_68_up;
  vector<double> mZP_Cross_95_up;

  vector<double> mN_Cross_95_down;
  vector<double> mN_Cross_68_down;
  vector<double> mN_Cross;
  vector<double> mN_Cross_68_up;
  vector<double> mN_Cross_95_up;
 
  theory_cross_50.clear();
  theory_cross_50.push_back(limit_theory["400_150"]);
  theory_cross_50.push_back(limit_theory["1600_750"]);
  theory_cross_50.push_back(limit_theory["2600_1250"]);
  theory_cross_50.push_back(limit_theory["2800_1350"]);
  
  limit_cross_50.clear();
  limit_cross_50.push_back(limit_exp["400_150"]);
  limit_cross_50.push_back(limit_exp["1600_750"]);
  limit_cross_50.push_back(limit_exp["2600_1250"]);
  limit_cross_50.push_back(limit_exp["2800_1350"]);

  limit_cross_95_down_50.clear();
  limit_cross_95_down_50.push_back(limit_exp_95_down["400_150"]);
  limit_cross_95_down_50.push_back(limit_exp_95_down["1600_750"]);
  limit_cross_95_down_50.push_back(limit_exp_95_down["2600_1250"]);
  limit_cross_95_down_50.push_back(limit_exp_95_down["2800_1350"]);

  limit_cross_68_down_50.clear();
  limit_cross_68_down_50.push_back(limit_exp_68_down["400_150"]);
  limit_cross_68_down_50.push_back(limit_exp_68_down["1600_750"]);
  limit_cross_68_down_50.push_back(limit_exp_68_down["2600_1250"]);
  limit_cross_68_down_50.push_back(limit_exp_68_down["2800_1350"]);

  limit_cross_68_up_50.clear();
  limit_cross_68_up_50.push_back(limit_exp_68_up["400_150"]);
  limit_cross_68_up_50.push_back(limit_exp_68_up["1600_750"]);
  limit_cross_68_up_50.push_back(limit_exp_68_up["2600_1250"]);
  limit_cross_68_up_50.push_back(limit_exp_68_up["2800_1350"]);

  limit_cross_95_up_50.clear();
  limit_cross_95_up_50.push_back(limit_exp_95_up["400_150"]);
  limit_cross_95_up_50.push_back(limit_exp_95_up["1600_750"]);
  limit_cross_95_up_50.push_back(limit_exp_95_up["2600_1250"]);
  limit_cross_95_up_50.push_back(limit_exp_95_up["2800_1350"]);
  bool HasExclusion[5] = {false, false, false, false, false};
  double previous_diff[5] = {999., 999., 999., 999., 999.};


  const int N_Zp_50 = 4;
  double mZp_50_arr[N_Zp_50] = {400., 1600., 2600., 2800.};
  double mN_50_arr[N_Zp_50] = {150., 750., 1250., 1350.};
  
  for(unsigned int i_50 = 0; i_50 < N_Zp_50; i_50++){
    cout << "i_50 : " << i_50 << endl;
    if(limit_cross_95_down_50.at(i_50) < theory_cross_50.at(i_50)) HasExclusion[0] = true;
    if(limit_cross_68_down_50.at(i_50) < theory_cross_50.at(i_50)) HasExclusion[1] = true;
    if(theory_cross_50.at(i_50) < theory_cross_50.at(i_50)) HasExclusion[2] = true;
    if(limit_cross_68_up_50.at(i_50) < theory_cross_50.at(i_50)) HasExclusion[3] = true;
    if(limit_cross_95_up_50.at(i_50) < theory_cross_50.at(i_50)) HasExclusion[4] = true;
    

    double current_diff[5];
    current_diff[0] = theory_cross_50.at(i_50) - limit_cross_95_down_50.at(i_50);
    current_diff[1] = theory_cross_50.at(i_50) - limit_cross_68_down_50.at(i_50);
    current_diff[2] = theory_cross_50.at(i_50) - limit_cross_50.at(i_50);
    current_diff[3] = theory_cross_50.at(i_50) - limit_cross_68_up_50.at(i_50);
    current_diff[4] = theory_cross_50.at(i_50) - limit_cross_95_up_50.at(i_50);
  
    
    if(i_50>0){
      double ZP1 = mZp_50_arr[i_50 - 1] + 0.;
      double dZP = (mZp_50_arr[i_50] + 0.) - (mZp_50_arr[i_50 - 1] + 0.);

      double T1 = theory_cross_50.at(i_50 - 1);
      double T2 = theory_cross_50.at(i_50);
      double dT = T2-T1;

      if(current_diff[0] * previous_diff[0] < 0){
	double L1 = limit_cross_95_down_50.at(i_50 - 1);
	double L2 = limit_cross_95_down_50.at(i_50);
	double dL = L2-L1;

	double x_meet = (L1-T1)*dZP/(dT-dL)+ZP1;
	mZP_Cross_95_down.push_back(x_meet);
	mN_Cross_95_down.push_back(0.5 * (x_meet + 0.) -50.0);
      }

      if(current_diff[1] * previous_diff[1] < 0){
        double L1 = limit_cross_68_down_50.at(i_50 - 1);
        double L2 = limit_cross_68_down_50.at(i_50);
        double dL = L2-L1;

        double x_meet = (L1-T1)*dZP/(dT-dL)+ZP1;
        mZP_Cross_68_down.push_back(x_meet);
        mN_Cross_68_down.push_back(0.5 * (x_meet + 0.) -50.0);
      }
      
      if(current_diff[2] * previous_diff[2] < 0){
	double L1 = limit_cross_50.at(i_50 - 1);
	double L2 = limit_cross_50.at(i_50);
	double dL = L2-L1;

	double x_meet = (L1-T1)*dZP/(dT-dL)+ZP1;
	mZP_Cross.push_back(x_meet);
	mN_Cross.push_back(0.5 * (x_meet + 0.) -50.0);
      }


      if(current_diff[3] * previous_diff[3] < 0){
	double L1 = limit_cross_68_up_50.at(i_50 - 1);
	double L2 = limit_cross_68_up_50.at(i_50);
	double dL = L2-L1;

	double x_meet = (L1-T1)*dZP/(dT-dL)+ZP1;
	mZP_Cross_68_up.push_back(x_meet);
	mN_Cross_68_up.push_back(0.5 * (x_meet + 0.) -50.0);
      }

      if(current_diff[4] * previous_diff[4] < 0){
        double L1 = limit_cross_95_up_50.at(i_50 - 1);
        double L2 = limit_cross_95_up_50.at(i_50);
        double dL = L2-L1;

        double x_meet = (L1-T1)*dZP/(dT-dL)+ZP1;
        mZP_Cross_95_up.push_back(x_meet);
        mN_Cross_95_up.push_back(0.5 * (x_meet + 0.) -50.0);
      }
    }

    previous_diff[0] = theory_cross_50.at(i_50) - limit_cross_95_down_50.at(i_50);
    previous_diff[1] = theory_cross_50.at(i_50) - limit_cross_68_down_50.at(i_50);
    previous_diff[2] = theory_cross_50.at(i_50) - limit_cross_50.at(i_50);
    previous_diff[3] = theory_cross_50.at(i_50) - limit_cross_68_up_50.at(i_50);
    previous_diff[4] = theory_cross_50.at(i_50) - limit_cross_95_up_50.at(i_50);

  }

  cout << "end 4 -50" << endl;
  
  if(mZP_Cross_95_down.size()>0){
    mZP_crosses_95_down.push_back(mZP_Cross_95_down.at(0));
    mN_crosses_95_down.push_back(mN_Cross_95_down.at(0));
  }
  else{
    mZP_crosses_95_down.push_back(1600.);
    mN_crosses_95_down.push_back(750.);
  }

  if(mZP_Cross_68_down.size()>0){
    mZP_crosses_68_down.push_back(mZP_Cross_68_down.at(0));
    mN_crosses_68_down.push_back(mN_Cross_68_down.at(0));
  }
  else{
    mZP_crosses_68_down.push_back(1600.);
    mN_crosses_68_down.push_back(750.);
  }
  
  if(mZP_Cross.size()>0){
    mZP_crosses.push_back(mZP_Cross.at(0));
    mN_crosses.push_back(mN_Cross.at(0));
  }
  else{
    mZP_crosses.push_back(1600.);
    mN_crosses.push_back(750.);
  }

  if(mZP_Cross_68_up.size()>0){
    mZP_crosses_68_up.push_back(mZP_Cross_68_up.at(0));
    mN_crosses_68_up.push_back(mN_Cross_68_up.at(0));
  }
  else{
    mZP_crosses_68_up.push_back(1600.);
    mN_crosses_68_up.push_back(750.);
  }
  
  if(mZP_Cross_95_up.size()>0){
    mZP_crosses_95_up.push_back(mZP_Cross_95_up.at(0));
    mN_crosses_95_up.push_back(mN_Cross_95_up.at(0));
  }
  else{
    mZP_crosses_95_up.push_back(1600.);
    mN_crosses_95_up.push_back(750.);
  }

  
  //========================================
  //==== 3) Let's draw found points 
  //========================================
  
  c_ct -> SetLogz();
  //exp_limit_2D -> Draw("colzsame");
    
  const int n_point_95_down = mZP_crosses_95_down.size()+1;
  const int n_point_68_down = mZP_crosses_68_down.size()+1;
  const int n_point = mZP_crosses.size()+1;
  const int n_point_68_up = mZP_crosses_68_up.size()+1;
  const int n_point_95_up = mZP_crosses_95_up.size()+1;
  
  double x_mZp_95_down[n_point_95_down];
  double x_mZp_68_down[n_point_68_down];
  double x_mZp[n_point];
  double x_mZp_68_up[n_point_68_up];
  double x_mZp_95_up[n_point_95_up];

  double y_mZp_95_down[n_point_95_down];
  double y_mZp_68_down[n_point_68_down];
  double y_mZp[n_point];
  double y_mZp_68_up[n_point_68_up];
  double y_mZp_95_up[n_point_95_up];
  for(int i_fill = 0; i_fill < n_point_95_down-2; i_fill++){
    x_mZp_95_down[i_fill] = mZP_crosses_95_down.at(i_fill+1);
    y_mZp_95_down[i_fill] = mN_crosses_95_down.at(i_fill+1);
  }
  for(int i_fill = 0; i_fill < n_point_68_down-2; i_fill++){
    x_mZp_68_down[i_fill] = mZP_crosses_68_down.at(i_fill+1);
    y_mZp_68_down[i_fill] = mN_crosses_68_down.at(i_fill+1);
  }
  for(int i_fill = 0; i_fill < n_point-2; i_fill++){
    x_mZp[i_fill] = mZP_crosses.at(i_fill+1);
    y_mZp[i_fill] = mN_crosses.at(i_fill+1);
  }
  for(int i_fill = 0; i_fill < n_point_68_up-2; i_fill++){
    x_mZp_68_up[i_fill] = mZP_crosses_68_up.at(i_fill+1);
    y_mZp_68_up[i_fill] = mN_crosses_68_up.at(i_fill+1);
  }
  for(int i_fill = 0; i_fill < n_point_95_up-2; i_fill++){
    x_mZp_95_up[i_fill] = mZP_crosses_95_up.at(i_fill+1);
    y_mZp_95_up[i_fill] = mN_crosses_95_up.at(i_fill+1);
  }
  x_mZp_95_down[n_point_95_down-2] = 400.;
  y_mZp_95_down[n_point_95_down-2] = 150.;
  x_mZp_68_down[n_point_68_down-2] = 400.;
  y_mZp_68_down[n_point_68_down-2] = 150.;
  x_mZp[n_point-2] = 400.;
  y_mZp[n_point-2] = 150.;
  x_mZp_68_up[n_point_68_up-2] = 400.;
  y_mZp_68_up[n_point_68_up-2] = 150.;
  x_mZp_95_up[n_point_95_up-2] = 400.;
  y_mZp_95_up[n_point_95_up-2] = 150.;

  x_mZp_95_down[n_point_95_down-1] = mZP_crosses_95_down.at(0);
  y_mZp_95_down[n_point_95_down-1] = mN_crosses_95_down.at(0);
  x_mZp_68_down[n_point_68_down-1] = mZP_crosses_68_down.at(0);
  y_mZp_68_down[n_point_68_down-1] = mN_crosses_68_down.at(0);
  x_mZp[n_point-1] = mZP_crosses.at(0);
  y_mZp[n_point-1] = mN_crosses.at(0);
  x_mZp_68_up[n_point_68_up-1] = mZP_crosses_68_up.at(0);
  y_mZp_68_up[n_point_68_up-1] = mN_crosses_68_up.at(0);
  x_mZp_95_up[n_point_95_up-1] = mZP_crosses_95_up.at(0);
  y_mZp_95_up[n_point_95_up-1] = mN_crosses_95_up.at(0);
  

 
  TGraph *gr_central = new TGraph(n_point, x_mZp, y_mZp);
  TGraph *gr_central_fill = (TGraph*)gr_central->Clone();
  gr_central_fill -> SetFillColor(kAzure + 8);
  gr_central_fill -> Draw("fsame");
  
  gr_central->Print();
  gr_central-> SetLineStyle(8);
  gr_central->SetLineWidth(3);
  gr_central->SetLineColor(kBlue);
  gr_central->Draw("lsame");
  
  TGraph *gr_68_down = new TGraph(n_point_68_down, x_mZp_68_down, y_mZp_68_down);
  gr_68_down->Print();
  gr_68_down->SetLineStyle(2);
  gr_68_down->SetLineWidth(1);
  gr_68_down->SetLineColor(kBlue);
  gr_68_down->Draw("lsame");
  
  TGraph *gr_68_up = new TGraph(n_point_68_up, x_mZp_68_up, y_mZp_68_up);
  gr_68_up->Print();
  gr_68_up->SetLineStyle(2);
  gr_68_up->SetLineWidth(1);
  gr_68_up->SetLineColor(kBlue);
  gr_68_up->Draw("lsame");
  
 
  //========================================
  //==== 4) Add ATLAS 8 TeV results 
  //========================================
  
  string atlas_line;
  vector<double> atlas_mZp, atlas_mN;
  ifstream atlas_in(WORKING_DIR+"/"+ATLAS_result_file);
  while(getline(atlas_in,atlas_line)){
    std::istringstream is( atlas_line );
    
    //atlas_line
    TString mark;
    double a,b;
    is >> mark;
    is >> a;
    is >> b;
    
    if(mark.Contains("TRUE")){
      atlas_mZp.push_back(a);
      atlas_mN.push_back(b);
    }
  }
  
  const int N_atlas = atlas_mZp.size();
  double mZp_atlas[N_atlas+2], mN_atlas[N_atlas+2];
  for(int i_atlas = 0; i_atlas < N_atlas; i_atlas++){
    mZp_atlas[i_atlas] = atlas_mZp.at(i_atlas);
    mN_atlas[i_atlas] = atlas_mN.at(i_atlas);
  }
  
  mZp_atlas[N_atlas] = 400.;
  mN_atlas[N_atlas] = 50.;
  mZp_atlas[N_atlas+1] = atlas_mZp.at(0);
  mN_atlas[N_atlas+1] = atlas_mN.at(0);

  TGraph *gr_atlas = new TGraph(N_atlas + 2, mZp_atlas, mN_atlas);
  gr_atlas->Print();
  gr_atlas-> SetLineStyle(8);
  gr_atlas->SetLineWidth(3);
  gr_atlas->SetLineColor(kRed);
  gr_atlas->Draw("lsame");
  

  //========================================
  //==== 5) Draw Legend and Latex 
  //========================================

  TLegend *lg = new TLegend(0.2, 0.7, 0.6, 0.90);
  lg->SetBorderSize(0);
  lg->SetFillStyle(0);
  lg->AddEntry(gr_Z2N, "m_{N} = m_{Z'} / 2 ", "l");
  lg->AddEntry(gr_central, "95% CL Expected Limit #pm 1 #sigma", "l");
  lg->AddEntry(gr_atlas, "#splitline{95% CL Observed Limit}{ATLAS @ 8 TeV}", "l");
  lg->Draw();

  TLatex latex_CMSPriliminary, latex_Lumi, latex_channel;
  latex_CMSPriliminary.SetNDC();
  latex_Lumi.SetNDC();
  latex_CMSPriliminary.SetTextSize(0.035);
  latex_CMSPriliminary.DrawLatex(0.16, 0.96, "#font[62]{CMS} #font[42]{#it{#scale[0.8]{Preliminary}}}");
  latex_Lumi.SetTextSize(0.035 );
  if(binning.Contains("2016"))latex_Lumi.DrawLatex(0.60, 0.96, "35.9 fb^{-1} (13 TeV, 2016)");
  else if(binning.Contains("2017"))latex_Lumi.DrawLatex(0.60, 0.96, "41.3 fb^{-1} (13 TeV, 2017)");
  else if(binning.Contains("2018")) latex_Lumi.DrawLatex(0.60, 0.96, "59.7 fb^{-1} (13 TeV, 2018)");
  else latex_Lumi.DrawLatex(0.60, 0.96, "137.4 fb^{-1} (13 TeV, Full Run2)");
  latex_channel.SetNDC();
  latex_channel.SetTextSize(0.08);
  if(channel.Contains("MuMu")) latex_channel.DrawLatex(0.2, 0.6, "#mu#mu OS + SS");
  else latex_channel.DrawLatex(0.2, 0.6, "ee OS + SS");
  
  
  dummy_ct_ForEachZP -> Draw("axissame");
  
  //==== Save

  c_ct->SaveAs(base_plotpath+"/Limit_2D.pdf");
  //c_ct->SaveAs(base_plotpath+"/Limit_2D.png");

  gr_central->SetName("Cotour");
  gr_central->Write();
  out->Close();


}
 





