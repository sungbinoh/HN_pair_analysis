#ifndef mylib_h
#define mylib_h

TGraphAsymmErrors* hist_to_graph(TH1D* hist, bool YErrorZero=false){

  TH1::SetDefaultSumw2(true);

  int Nbins = hist->GetXaxis()->GetNbins();
  double x[Nbins], y[Nbins], xlow[Nbins], xup[Nbins], ylow[Nbins], yup[Nbins];
  TAxis *xaxis = hist->GetXaxis();
  for(Int_t i=0; i<Nbins; i++){
    x[i] = xaxis->GetBinCenter(i+1);
    y[i] = hist->GetBinContent(i+1);
    xlow[i] = xaxis->GetBinCenter(i+1)-xaxis->GetBinLowEdge(i+1);
    xup[i] = xaxis->GetBinUpEdge(i+1)-xaxis->GetBinCenter(i+1);
    ylow[i] = hist->GetBinError(i+1);
    yup[i] = hist->GetBinError(i+1);
    if(YErrorZero){
      ylow[i] = 0;
      yup[i] = 0;
    }
    //ylow[i] = 0;
    //yup[i] = 0;
    //cout << "x = " << x[i] << ", y = " << y[i] << ", x_low = " << xlow[i] << ", xup = " << xup[i] << ", ylow = " << ylow[i] << ", yup = " << yup[i] << endl;
  }
  TGraphAsymmErrors *out = new TGraphAsymmErrors(Nbins, x, y, xlow, xup, ylow, yup);
  out->SetTitle("");
  return out;

}

TGraphAsymmErrors* hist_to_graph(TH1D* hist, int n_skip_x_left){

  TH1::SetDefaultSumw2(true);

  int Nbins = hist->GetXaxis()->GetNbins()-n_skip_x_left;
  double x[Nbins], y[Nbins], xlow[Nbins], xup[Nbins], ylow[Nbins], yup[Nbins];
  TAxis *xaxis = hist->GetXaxis();
  for(Int_t i=1; i<=Nbins; i++){
    x[i-1] = xaxis->GetBinCenter(i+n_skip_x_left);
    y[i-1] = hist->GetBinContent(i+n_skip_x_left);
    xlow[i-1] = xaxis->GetBinCenter(i+n_skip_x_left)-xaxis->GetBinLowEdge(i+n_skip_x_left);
    xup[i-1] = xaxis->GetBinUpEdge(i+n_skip_x_left)-xaxis->GetBinCenter(i+n_skip_x_left);
    ylow[i-1] = hist->GetBinError(i+n_skip_x_left);
    yup[i-1] = hist->GetBinError(i+n_skip_x_left);
    //ylow[i-1] = 0;
    //yup[i-1] = 0;
    //cout << "x = " << x[i-1] << ", y = " << y[i-1] << ", x_low = " << xlow[i-1] << ", xup = " << xup[i-1] << ", ylow = " << ylow[i-1] << ", yup = " << yup[i-1] << endl;
  }
  TGraphAsymmErrors *out = new TGraphAsymmErrors(Nbins, x, y, xlow, xup, ylow, yup);
  out->SetTitle("");
  return out;

}

TGraphAsymmErrors* hist_to_graph(TH1D* hist, int n_skip_x_left, int n_x_shift, int i_x_shift){

  TH1::SetDefaultSumw2(true);

  int Nbins = hist->GetXaxis()->GetNbins()-n_skip_x_left;
  double x[Nbins], y[Nbins], xlow[Nbins], xup[Nbins], ylow[Nbins], yup[Nbins];
  TAxis *xaxis = hist->GetXaxis();
  for(Int_t i=1; i<=Nbins; i++){
    x[i-1] = xaxis->GetBinCenter(i+n_skip_x_left);
    y[i-1] = hist->GetBinContent(i+n_skip_x_left);
    xlow[i-1] = xaxis->GetBinCenter(i+n_skip_x_left)-xaxis->GetBinLowEdge(i+n_skip_x_left);
    xup[i-1] = xaxis->GetBinUpEdge(i+n_skip_x_left)-xaxis->GetBinCenter(i+n_skip_x_left);
    ylow[i-1] = hist->GetBinError(i+n_skip_x_left);
    yup[i-1] = hist->GetBinError(i+n_skip_x_left);

    double dx = (xaxis->GetBinUpEdge(i+n_skip_x_left)-xaxis->GetBinLowEdge(i+n_skip_x_left))/double(n_x_shift+1);
    x[i-1] = xaxis->GetBinLowEdge(i+n_skip_x_left) + double(i_x_shift+1) * dx;
    xlow[i-1] = double(i_x_shift+1) * dx;
    xup[i-1] = xaxis->GetBinUpEdge(i+n_skip_x_left)-x[i-1];

    //ylow[i-1] = 0;
    //yup[i-1] = 0;
    //cout << "x = " << x[i-1] << ", y = " << y[i-1] << ", x_low = " << xlow[i-1] << ", xup = " << xup[i-1] << ", ylow = " << ylow[i-1] << ", yup = " << yup[i-1] << endl;
  }
  TGraphAsymmErrors *out = new TGraphAsymmErrors(Nbins, x, y, xlow, xup, ylow, yup);
  out->SetTitle("");
  return out;

}


TGraphAsymmErrors* GraphSubtract(TGraphAsymmErrors *a, TGraphAsymmErrors *b, bool Rel){

  //==== do a-b

  int NX = a->GetN();
  TGraphAsymmErrors* gr_out = (TGraphAsymmErrors*)a->Clone();

  for(int i=0; i<NX; i++){

    double a_x, a_y, b_x, b_y;

    a->GetPoint(i, a_x, a_y);
    b->GetPoint(i, b_x, b_y);

    if(Rel==true){
      gr_out->SetPoint( i, a_x, (a_y-b_y)/b_y );
    }
    else{
      gr_out->SetPoint( i, a_x, a_y-b_y );
    }

  }

  return gr_out;

}

void RemoveLargeError(TGraphAsymmErrors *a){

  int NX = a->GetN();

  for(int i=0; i<NX; i++){

    double x, y, yerr_low, yerr_high;

    a->GetPoint(i, x, y);
    yerr_low  = a->GetErrorYlow(i);
    yerr_high = a->GetErrorYhigh(i);

    if(y+yerr_high==1.){
      a->SetPointEYhigh( i, yerr_low );
    }

  }

}

void ScaleGraph(TGraphAsymmErrors *a, double c){

  int NX = a->GetN();

  for(int i=0; i<NX; i++){

    double x, y, yerr_low, yerr_high;

    a->GetPoint(i, x, y);
    yerr_low  = a->GetErrorYlow(i);
    yerr_high = a->GetErrorYhigh(i);

    a->SetPoint(i, x, c*y);
    a->SetPointEYlow(i, c*yerr_low);
    a->SetPointEYhigh(i, c*yerr_high);

  }

}



double GetMaximum(TH1D* hist){

  TAxis *xaxis = hist->GetXaxis();

  double maxval(-1.);
  for(int i=1; i<=xaxis->GetNbins(); i++){
    if( hist->GetBinContent(i) + hist->GetBinError(i) > maxval ){
      maxval = hist->GetBinContent(i) + hist->GetBinError(i);
    }
  }

  return maxval;

}

double GetMaximum(TGraphAsymmErrors *a){

  int NX = a->GetN();

  double maxval(-9999.);
  for(int i=0; i<NX; i++){

    double x, y, yerr_low, yerr_high;

    a->GetPoint(i, x, y);
    yerr_low  = a->GetErrorYlow(i);
    yerr_high = a->GetErrorYhigh(i);

    if( y+yerr_high > maxval ){
      maxval = y+yerr_high;
    }

  }

  return maxval;

}

double GetYieldSystematics(TH1D *hist){

  int n_syst = hist->GetXaxis()->GetNbins();
  int n_source = (n_syst-1)/2;

  //==== Bin index
  //==== i=1 : central
  //==== i=2 : _MuonEn_up
  //==== i=3 : _MuonEn_down
  //==== -> n_syst = 3
  //==== -> n_source = (n_syst-1)/2 = (3-1)/2 = 1

  double y_central = hist->GetBinContent(1);

  double sum_syst = 0.;
  for(int i=1; i<=n_source; i++){
    double yield_up = hist->GetBinContent(i*2);
    double yield_down = hist->GetBinContent(i*2+1);

    double syst_up = fabs(yield_up-y_central);
    double syst_down = fabs(yield_down-y_central);

    sum_syst += 0.5*(syst_up*syst_up+syst_down*syst_down);

  }

  return sqrt(sum_syst);

}

TDirectory *MakeTemporaryDirectory(){

  gROOT->cd();
  TDirectory* tempDir = 0;
  int counter = 0;
  while (not tempDir) {
    // First, let's find a directory name that doesn't exist yet:
    std::stringstream dirname;
    dirname << "HNCommonLeptonFakes_%i" << counter;
    if (gROOT->GetDirectory((dirname.str()).c_str())) {
      ++counter;
      continue;
    }
    // Let's try to make this directory:
    tempDir = gROOT->mkdir((dirname.str()).c_str());

  }

  return tempDir;

}

void AddPhantomZero(double a, TString align, int digit_int, int digit_frac){

  if(align=="r"){
    int number_maxdigit = 0;
    for(int i=10; i>=0; i--){
      if(a/pow(10,i)>=1.){
        number_maxdigit = i;
        break;
      }
    }
    //cout << "number_maxdigit = " << number_maxdigit << endl;
    for(int i=0; i<digit_int-(number_maxdigit+1); i++) cout << "\\phantom{0}";
    cout << std::fixed<<std::setprecision(digit_frac) << a;
  }

  else if(align=="l"){
    int target_total_digit = digit_int+digit_frac;
    int number_maxdigit = 0;
    for(int i=10; i>=0; i--){
      if(a/pow(10,i)>=1.){
        number_maxdigit = i;
        break;
      }
    }
    //cout << "target_total_digit = " << target_total_digit << endl;
    //cout << "number_maxdigit = " << number_maxdigit << endl;
    //cout << "--> # of \\phantom{0} = " << target_total_digit-digit_int-(number_maxdigit+1) << endl;

    cout << std::fixed<<std::setprecision(digit_frac) << a;
    for(int i=0; i<target_total_digit-(number_maxdigit+1)-digit_frac; i++) cout << "\\phantom{0}";
  }

}

vector<int> GetZPMassRange(int mN){

  vector<int> this_masses;

  if(mN==100){
    this_masses = {
      400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 
    };
    return this_masses;
  }
  if(mN==200){
    this_masses = {
      600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==300){
    this_masses = {
      800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==400){
    this_masses = {
      1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==500){
    this_masses = {
      1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==600){
    this_masses = {
      1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==700){
    this_masses = {
      1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==800){
    this_masses = {
      1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==900){
    this_masses = {
      2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1000){
    this_masses = {
      2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1100){
    this_masses = {
      2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1200){
    this_masses = {
      2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1300){
    this_masses = {
      2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1400){
    this_masses = {
      3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1500){
    this_masses = {
      3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1600){
    this_masses = {
      3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1700){
    this_masses = {
      3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1800){
    this_masses = {
      3800, 4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==1900){
    this_masses = {
      4000, 4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==2000){
    this_masses = {
      4200, 4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==2100){
    this_masses = {
      4400, 4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==2200){
    this_masses = {
      4600, 4800, 5000,
    };
    return this_masses;
  }
  if(mN==2300){
    this_masses = {
      4800, 5000,
    };
    return this_masses;
  }
  if(mN==2400){
    this_masses = {
      5000,
    };
    return this_masses;
  }
  
  return this_masses;
  
}

vector<int> GetHNMassRange(int mZP, bool useallmass=false){

  vector<int> this_masses;

  if(useallmass){

    if(mZP==400){
      this_masses = {
	100, 
      };
      return this_masses;
    }
    if(mZP==600){
      this_masses = {
        100, 200,
      };
      return this_masses;
    }
    if(mZP==800){
      this_masses = {
	100, 200, 300,
      };
      return this_masses;
    }
    if(mZP==1000){
      this_masses = {
        100, 200, 300, 400,
      };
      return this_masses;
    }
    if(mZP==1200){
      this_masses = {
	100, 200, 300, 400, 500, 
      };
      return this_masses;
    }
    if(mZP==1400){
      this_masses = {
        100, 200, 300, 400, 500, 600,
      };
      return this_masses;
    }
    if(mZP==1600){
      this_masses = {
	100, 200, 300, 400, 500, 600, 700, 
      };
      return this_masses;
    }
    if(mZP==1800){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800,
      };
      return this_masses;
    }
    if(mZP==2000){
      this_masses = {
	100, 200, 300, 400, 500, 600, 700, 800, 900, 
      };
      return this_masses;
    }
    if(mZP==2200){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
      };
      return this_masses;
    }
    if(mZP==2400){
      this_masses = {
	100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 
      };
      return this_masses;
    }
    if(mZP==2600){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200,
      };
      return this_masses;
    }
    if(mZP==2800){
      this_masses = {
	100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 
      };
      return this_masses;
    }
    if(mZP==3000){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400,
      };
      return this_masses;
    }
    if(mZP==3200){
      this_masses = {
	100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 
      };
      return this_masses;
    }
    if(mZP==3400){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600,
      };
      return this_masses;
    }
    if(mZP==3600){
      this_masses = {
	100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700,
      };
      return this_masses;
    }
    if(mZP==3800){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800,
      };
      return this_masses;
    }
    if(mZP==4000){
      this_masses = {
	100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900,
      };
      return this_masses;
    }
    if(mZP==4200){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000,
      };
      return this_masses;
    }
    if(mZP==4400){
      this_masses = {
	100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100,
      };
      return this_masses;
    }
    if(mZP==4600){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200,
      };
      return this_masses;
    }
    if(mZP==4800){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300,
      };
      return this_masses;
    }
    if(mZP==5000){
      this_masses = {
        100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400,
      };
      return this_masses;
    }
    
  }
  else{
    
    if(mZP==400){
      this_masses = {
        100,
      };
      return this_masses;
    }
    if(mZP==800){
      this_masses = {
	100, 300,
      };
      return this_masses;
    }
    if(mZP==1200){
      this_masses = {
	//100, 300, 500,
	100, 500,
      };
      return this_masses;
    }
    if(mZP==1600){
      this_masses = {
	//100, 300, 500, 700,
	100, 700,
      };
      return this_masses;
    }
    if(mZP==2000){
      this_masses = {
	//100, 300, 500, 700, 900,
	100, 900,
      };
      return this_masses;
    }
    if(mZP==2400){
      this_masses = {
	//100, 300, 500, 700, 900, 1100,
	100, 1100,
      };
      return this_masses;
    }
    if(mZP==2800){
      this_masses = {
	//100, 300, 500, 700, 900, 1100, 1300,
	100, 1300,
      };
      return this_masses;
    }
    if(mZP==3200){
      this_masses = {
	//100, 300, 500, 700, 900, 1100, 1300, 1500,
	100, 1500,
      };
      return this_masses;
    }
    if(mZP==3600){
      this_masses = {
	//100, 300, 500, 700, 900, 1100, 1300, 1500, 1700,
	100, 700, 1700,
      };
      return this_masses;
    }
    if(mZP==4000){
      this_masses = {
	//100, 300, 500, 700, 900, 1100, 1300, 1500, 1700, 1900,
	100, 900, 1900,
      };
      return this_masses;
    }
    if(mZP==4400){
      this_masses = {
	//100, 300, 500, 700, 900, 1100, 1300, 1500, 1700, 1900, 2100,
	100, 1100, 2100,
      };
      return this_masses;
    }
    if(mZP==4800){
      this_masses = {
        //100, 300, 500, 700, 900, 1100, 1300, 1500, 1700, 1900, 2100, 2300,
	100, 1300, 2300,
      };
      return this_masses;
    }


  }

  cout << "[GetHNMassRange(int mZP)] Wrong mZP!!" << endl;
  return this_masses;




}

bool IsSignalCATSamaple(TString samplename){

  if(samplename.Contains("HNpair")) return true;

  return false;

}


int GetLowMN(double current_mN){
  current_mN = current_mN / 100.;
  int low_mN = current_mN;
  low_mN = low_mN * 100;
  
  return low_mN;
}









#endif
