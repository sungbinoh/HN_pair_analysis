#include "canvas_margin.h"
#include "mylib.h"

void QuickLimit(int xxx=0){

  setTDRStyle();

  gStyle->SetOptStat(0);

  TH1::SetDefaultSumw2(true);
  TH1::AddDirectory(kFALSE);

  TString WORKING_DIR = getenv("PLOTTER_WORKING_DIR");
  TString dataset = getenv("CATANVERSION");
  TString ENV_PLOT_PATH = getenv("PLOT_PATH");

  TString dirname = "180930_OS_bkgd30";

  TString channel = "MuMu";
  TString method = "shape";

  if(xxx==1){
    channel = "MuMu";
    method = "count";
  }
  if(xxx==2){
    channel = "ElEl";
    method = "shape";
  }
  if(xxx==3){
    channel = "ElEl";
    method = "count";
  }

  TString base_filepath = WORKING_DIR+"/rootfiles/"+dataset+"/Limit/"+dirname+"/";
  TString base_plotpath = ENV_PLOT_PATH+"/"+dataset+"/QuickLimit/"+dirname+"/"+channel+"/"+method;

  if( !gSystem->mkdir(base_plotpath, kTRUE) ){
    cout
    << "###################################################" << endl
    << "Directoy " << base_plotpath << " is created" << endl
    << "###################################################" << endl
    << endl;
  }
  TFile *out = new TFile(base_plotpath+"/Limits.root","RECREATE");

  TString XSEC_file = "script/2017SampleXsec/values_MG_NLO.txt";


  //=================================
  //==== 1) xsec vs mN, for each ZP
  //=================================

  cout <<" #### 1) xsec vs mN, for each ZP ####" << endl;

  vector<int> ZPmasses = {400, 800, 1200, 1600, 2000, 2400, 2800, 3200, 3600, 4000, 4400, 4800};
  vector<double> nN_crosses;

  for(unsigned int i=0; i<ZPmasses.size(); i++){

    int mZP = ZPmasses.at(i);
    vector<int> hnmasses = GetHNMassRange(mZP,true);

    const int n_Nmass = hnmasses.size();

    double x[n_Nmass], y_theory[n_Nmass], y_exp[n_Nmass];

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
      ifstream in2(base_filepath+"/results_"+method+"_"+channel+".txt");
      double this_limit(1.);
      while(getline(in2,elline2)){
        std::istringstream is( elline2 );

        // 2000  900 0.0042160
        TString a;
        double b;
        is >> a;
        is >> b;
        if(a==this_string) this_limit = b*0.001;
      }

      y_theory[j] = this_xsec;
      y_exp[j]= this_limit;

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

      //cout << this_string << "\t" << this_xsec << "\t" << this_limit << endl;

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

    TGraph *Theory_for_this_mZP = new TGraph(n_Nmass,x,y_theory);
    TGraph *Limit_for_this_mZP = new TGraph(n_Nmass,x,y_exp);

/*
    cout << "#### Z' = " << mZP << " ####" << endl;
    cout << "---- Theory ----" << endl;
    Theory_for_this_mZP->Print();
    cout << "---- Limit ----" << endl;
    Limit_for_this_mZP->Print();
*/

    TCanvas *c1 = new TCanvas("c1", "", 600, 600);
    canvas_margin(c1);
    c1->cd();

    TH1D *hist_dummy = new TH1D("hist_dummy", "", 5000, 0., 5000.);
    hist_axis(hist_dummy);
    hist_dummy->Draw("hist");
    hist_dummy->GetYaxis()->SetRangeUser(1E-6,100);
    hist_dummy->GetYaxis()->SetTitle("#sigma (pb)");
    hist_dummy->GetXaxis()->SetRangeUser(hnmasses.at(0),hnmasses.at(hnmasses.size()-1));
    hist_dummy->GetXaxis()->SetTitle("m_{N} (GeV)");
    c1->SetLogy();

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

    TLegend *lg = new TLegend(0.2, 0.2, 0.5, 0.4);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(Theory_for_this_mZP, "Theory", "lp");
    lg->AddEntry(Limit_for_this_mZP, "Expected", "lp");
    lg->Draw();

    TLatex latex_zpmass;
    latex_zpmass.SetTextSize(0.05);
    latex_zpmass.SetNDC();
    latex_zpmass.DrawLatex(0.6, 0.3,"m_{Z'} = "+TString::Itoa(mZP,10)+" GeV");

    c1->SaveAs(base_plotpath+"/Limit_Zp"+TString::Itoa(mZP,10)+".pdf");
    c1->SaveAs(base_plotpath+"/Limit_Zp"+TString::Itoa(mZP,10)+".png");

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

  vector<int> Nmasses = {100, 300, 500, 700, 900, 1100, 1300, 1500, 1700, 1900, 2100, 2300};
  vector<double> mZP_crosses;

  for(unsigned int i=0; i<Nmasses.size(); i++){

    int mN = Nmasses.at(i);
    vector<int> zpmasses = GetZPMassRange(mN);

    const int n_ZPmass = zpmasses.size();

    double x[n_ZPmass], y_theory[n_ZPmass], y_exp[n_ZPmass];

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
      ifstream in2(base_filepath+"/results_"+method+"_"+channel+".txt");
      double this_limit(1.);
      while(getline(in2,elline2)){
        std::istringstream is( elline2 );

        // 2000  900 0.0042160
        TString a;
        double b;
        is >> a;
        is >> b;
        if(a==this_string) this_limit = b*0.001;
      }

      y_theory[j] = this_xsec;
      y_exp[j]= this_limit;

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

    TGraph *Theory_for_this_mN = new TGraph(n_ZPmass,x,y_theory);
    TGraph *Limit_for_this_mN = new TGraph(n_ZPmass,x,y_exp);

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
    hist_dummy->GetXaxis()->SetRangeUser(zpmasses.at(0),zpmasses.at(zpmasses.size()-1));
    hist_dummy->GetXaxis()->SetTitle("m_{Z'} (GeV)");
    c1->SetLogy();

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

    TLegend *lg = new TLegend(0.2, 0.2, 0.5, 0.4);
    lg->SetBorderSize(0);
    lg->SetFillStyle(0);
    lg->AddEntry(Theory_for_this_mN, "Theory", "lp");
    lg->AddEntry(Limit_for_this_mN, "Expected", "lp");
    lg->Draw();

    TLatex latex_zpmass;
    latex_zpmass.SetTextSize(0.05);
    latex_zpmass.SetNDC();
    latex_zpmass.DrawLatex(0.6, 0.3,"m_{N} = "+TString::Itoa(mN,10)+" GeV");

    c1->SaveAs(base_plotpath+"/Limit_N"+TString::Itoa(mN,10)+".pdf");
    c1->SaveAs(base_plotpath+"/Limit_N"+TString::Itoa(mN,10)+".png");

    out->cd();
    Theory_for_this_mN->SetName("Zp"+TString::Itoa(mN,10)+"_Theory");
    Limit_for_this_mN->SetName("Zp"+TString::Itoa(mN,10)+"_Limit");
    Theory_for_this_mN->Write();
    Limit_for_this_mN->Write();

    c1->Close();

  }


  //=============================
  //==== Okay, Now Draw Contour
  //=============================

  TCanvas *c_ct = new TCanvas("c_ct", "", 600, 600);
  canvas_margin(c_ct);
  c_ct->cd();

  TH1D *dummy_ct_ForEachZP = new TH1D("dummy_ct_ForEachZP", "", 4000, 0., 4000.);
  hist_axis(dummy_ct_ForEachZP);
  dummy_ct_ForEachZP->GetYaxis()->SetLabelSize(0.04);
  dummy_ct_ForEachZP->GetYaxis()->SetTitleSize(0.06);
  dummy_ct_ForEachZP->GetYaxis()->SetTitleOffset(1.25);
  dummy_ct_ForEachZP->Draw("hist");
  dummy_ct_ForEachZP->GetYaxis()->SetRangeUser(50., 2000.);
  dummy_ct_ForEachZP->GetYaxis()->SetTitle("m_{N} (GeV)");
  dummy_ct_ForEachZP->GetXaxis()->SetRangeUser(500., 4000.);
  dummy_ct_ForEachZP->GetXaxis()->SetTitle("m_{Z'} (GeV)");

  double x_Z2N[2] = {0,10000};
  double y_Z2N[2] = {0,5000};
  TGraph *gr_Z2N = new TGraph(2,x_Z2N,y_Z2N);
  gr_Z2N->SetLineStyle(3);
  gr_Z2N->SetLineWidth(3);
  gr_Z2N->Draw("lsame");

  //========================================
  //==== 1) For each ZP mass, find N range
  //========================================

  const int n_ZP_has_exlclusion = nN_crosses.size()/2;
  double x_ct_ForEachZP[1+2*n_ZP_has_exlclusion], y_ct_ForEachZP[1+2*n_ZP_has_exlclusion];
  for(int i=0; i<n_ZP_has_exlclusion; i++){
    x_ct_ForEachZP[i] = ZPmasses.at(i);
    y_ct_ForEachZP[i] = nN_crosses.at(1+2*i);
  }
  for(int i=0; i<n_ZP_has_exlclusion; i++){
    int j = n_ZP_has_exlclusion-1-i;
    x_ct_ForEachZP[n_ZP_has_exlclusion+i] = ZPmasses.at(j);
    y_ct_ForEachZP[n_ZP_has_exlclusion+i] = nN_crosses.at(2*j);
  }
  x_ct_ForEachZP[2*n_ZP_has_exlclusion] = x_ct_ForEachZP[0];
  y_ct_ForEachZP[2*n_ZP_has_exlclusion] = y_ct_ForEachZP[0];


  TGraph *gr_ForEachZP = new TGraph(1+2*n_ZP_has_exlclusion,x_ct_ForEachZP,y_ct_ForEachZP);
  gr_ForEachZP->Print();
  gr_ForEachZP->SetLineColor(kBlue);
  gr_ForEachZP->SetLineWidth(2);
  gr_ForEachZP->SetMarkerColor(kBlue);
  gr_ForEachZP->SetMarkerStyle(15);
  gr_ForEachZP->Draw("lpsame");

  //========================================
  //==== 2) For each N mass, find ZP range
  //========================================

  const int n_N_has_exlclusion = mZP_crosses.size()/2;
  double x_ct_ForEachN[1+2*n_N_has_exlclusion], y_ct_ForEachN[1+2*n_N_has_exlclusion];
  for(int i=0; i<n_N_has_exlclusion; i++){
    x_ct_ForEachN[i] = mZP_crosses.at(1+2*i);
    y_ct_ForEachN[i] = Nmasses.at(i);

  }
  for(int i=0; i<n_N_has_exlclusion; i++){
    int j = n_N_has_exlclusion-1-i;
    x_ct_ForEachN[n_N_has_exlclusion+i] = mZP_crosses.at(2*j);
    y_ct_ForEachN[n_N_has_exlclusion+i] = Nmasses.at(j);

  }
  x_ct_ForEachN[2*n_N_has_exlclusion] = x_ct_ForEachN[0];
  y_ct_ForEachN[2*n_N_has_exlclusion] = y_ct_ForEachN[0];


  TGraph *gr_ForEachN = new TGraph(1+2*n_N_has_exlclusion,x_ct_ForEachN,y_ct_ForEachN);
  gr_ForEachN->Print();
  gr_ForEachN->SetLineWidth(2);
  gr_ForEachN->SetLineColor(kRed);
  gr_ForEachN->SetMarkerColor(kRed);
  gr_ForEachN->SetMarkerStyle(15);
  gr_ForEachN->Draw("lpsame");

  //==== Save

  c_ct->SaveAs(base_plotpath+"/Limit_2D.pdf");
  c_ct->SaveAs(base_plotpath+"/Limit_2D.png");

  gr_ForEachZP->SetName("Cotour");
  gr_ForEachZP->Write();
  out->Close();


}









