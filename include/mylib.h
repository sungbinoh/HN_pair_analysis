#ifndef mylib_h
#define mylib_h

ofstream outfile;

using namespace std;
const double alpha = 1 - 0.6827;
// == Debugging Mode
bool debug = true;
int tag_year = 0;
TString Cycle_name;
TString SingleMuon;
TString DoubleEG;
// == Call all needed maps
map<TString, TH1*> maphist;
map<TString, TH2D*> maphist2D;
map<TString, TGraphAsymmErrors*> map_asym_gr;
map<TString, TFile*> mapfile;
map<TString, TCanvas*> mapcanvas;
map<TString, TPad*> mappad;
map<TString, THStack*> maphstack;
map<TString, TLegend*> maplegend;
map<TString, TH1F*> mapfunc;
map<TString, double> map_overflow;
map<TString, TLine*> mapline;
map<TString, TKey*> maphistcheck;
map<TString, TList*> maplist;
map<TString, std::vector<double> > map_bin_vector;
map<TString, std::vector<TString> > map_sample_names;
map<TString, std::vector<double> > map_syst_array;

const int N_syst = 19;
TString systematics[N_syst] = {"central",
			       "ElectronScaleDown", "ElectronScaleUp",
			       "ElectronSmearDown", "ElectronSmearUp",
			       "JetsResDown", "JetsResUp",
			       "JetsScaleDown", "JetsScaleUp",
			       "SD_JMR_Down", "SD_JMR_Up",
			       "SD_JMS_Down", "SD_JMS_Up",
			       "PUReweight_Down", "PUReweight_Up",
			       "MuonScaleDown", "MuonScaleUp",
			       "MuonSmearDown", "MuonSmearUp",
			       //"PDFNormDown", "PDFNormUp" 
};

const int N_syst_comparison = 9;
TString systematics_comparison[N_syst] = {"ElectronScale",
					  "ElectronSmear",
					  "JetsRes",
					  "JetsScale",
					  "SD_JMR_",
					  "SD_JMS_",
					  "PUReweight_",
					  "MuonScale",
					  "MuonSmear"
};

vector<double> vx, vy, vexl, vexh, veyl, veyh;

TH1F * GetHist(TString hname){

  TH1F *h = NULL;
  std::map<TString, TH1F*>::iterator mapit = mapfunc.find(hname);
  if(mapit != mapfunc.end()) return mapit-> second;

  return h;

}

vector<double> Get_Syst_Error(TString nameofhistogram, TString sample){
  
  vector<double> error_vector;
  error_vector.clear();
  if(debug) cout << nameofhistogram << ", " << sample << endl;  
  Int_t N_rebinned_x = mapfunc[nameofhistogram + Cycle_name + sample + "rebin"] -> GetNbinsX();
  
  for(Int_t j = 1; j < N_rebinned_x; j++){
    double current_bin_error =  mapfunc[nameofhistogram + Cycle_name + sample + "rebin"] -> GetBinError(j);
    for(int i_syst = 1; i_syst < N_syst; i_syst++){
      if(GetHist(nameofhistogram + systematics[i_syst] + sample + "rebin")){
	double diff = fabs(mapfunc[nameofhistogram + Cycle_name + sample + "rebin"] -> GetBinContent(j) - mapfunc[nameofhistogram + systematics[i_syst] + sample + "rebin"] -> GetBinContent(j));
	if(debug) cout << systematics[i_syst] << " : " << diff << endl;;
	current_bin_error = pow(current_bin_error, 2);
	current_bin_error += pow(diff, 2);
	current_bin_error = pow(current_bin_error, 0.5);
      }
    }
    error_vector.push_back(current_bin_error);
    if(debug) cout << mapfunc[nameofhistogram + Cycle_name + sample + "rebin"] -> GetBinContent(j) << ", syst : " << current_bin_error << endl;
  }

  return error_vector;
}

void Rebin_with_overflow(TString histname, int N_bin, double binx[]){

  if(debug) cout << "[Rebin_with_overflow]" << endl;

  Int_t nx    = mapfunc[histname] -> GetNbinsX()+1;
  Double_t x1 = mapfunc[histname] -> GetBinLowEdge(1);
  Double_t bw = (binx[N_bin - 1] - binx[0]) / 20.;
  Double_t x2 = mapfunc[histname] -> GetBinLowEdge(nx)+bw;

  TH1F *htmp = new TH1F("", "", nx, x1, x2);
  for(Int_t j = 1; j <= nx; j++){
    double current_bin_content = mapfunc[histname] -> GetBinContent(j);
    double current_bin_error = mapfunc[histname] -> GetBinError(j);
    htmp -> SetBinContent(j, current_bin_content);
    htmp -> SetBinError(j, current_bin_error);
  }
  mapfunc[histname + "rebin"] = dynamic_cast<TH1F*>(htmp -> Rebin(N_bin, histname + "rebin", binx));

}

void Save_syst_array(TString current_histname, TString systematics, TString cycle_and_sample, int N_bin, double binx[]){
  
  cout << "[Save_syst_array]" << endl;
  
  TString histname_Up   = current_histname + "_" + systematics + "Up"   + cycle_and_sample;
  TString histname_Down = current_histname + "_" + systematics + "Down" + cycle_and_sample;
  TString histname_central = current_histname + "_central" + cycle_and_sample + "rebin"; // -- central hists are already rebinned

  // -- booleans to check if there are up/down hists
  bool up_hist = false;
  bool down_hist = false;
  
  if(!mapfunc[histname_Up]) return;

  Int_t nx    = mapfunc[histname_Up] -> GetNbinsX()+1;
  Double_t x1 = mapfunc[histname_Up] -> GetBinLowEdge(1);
  Double_t bw = (binx[N_bin - 1] - binx[0]) / 20.;
  Double_t x2 = mapfunc[histname_Up] -> GetBinLowEdge(nx)+bw;

  TH1F *htmp_syst_Up = new TH1F("", "", nx, x1, x2);
  TH1F *htmp_syst_Down = new TH1F("", "", nx, x1, x2);

  for(Int_t j = 1; j <= nx; j++){
    double current_bin_content_Up = mapfunc[histname_central] -> GetBinContent(j);
    double current_bin_error_Up = mapfunc[histname_central] -> GetBinError(j);
    if(mapfunc[histname_Up]){
      current_bin_content_Up = mapfunc[histname_Up] -> GetBinContent(j);
      current_bin_error_Up = mapfunc[histname_Up] -> GetBinError(j);
      up_hist =true;
    }
    htmp_syst_Up -> SetBinContent(j, current_bin_content_Up);
    htmp_syst_Up -> SetBinError(j, current_bin_error_Up);
    
    double current_bin_content_Down = mapfunc[histname_central] -> GetBinContent(j);
    double current_bin_error_Down = mapfunc[histname_central] -> GetBinError(j);
    if(mapfunc[histname_Down]){
      current_bin_content_Down = mapfunc[histname_Down] -> GetBinContent(j);
      current_bin_error_Down = mapfunc[histname_Down] -> GetBinError(j);
      down_hist = true;
    }
    htmp_syst_Down -> SetBinContent(j, current_bin_content_Down);
    htmp_syst_Down -> SetBinError(j, current_bin_error_Down);
  }
  
  mapfunc[histname_Up + "rebin"] = dynamic_cast<TH1F*>(htmp_syst_Up -> Rebin(N_bin, histname_Up + "rebin", binx));
  mapfunc[histname_Down + "rebin"] = dynamic_cast<TH1F*>(htmp_syst_Down -> Rebin(N_bin, histname_Down + "rebin", binx));

  map_syst_array[current_histname + systematics + "Up"].clear();
  map_syst_array[current_histname + systematics + "Down"].clear();
  //map_syst_array
  Int_t new_nx = mapfunc[histname_Up + "rebin"] -> GetNbinsX() + 1;
  for(Int_t j = 1; j < new_nx; j++){
    double current_entry_Up = mapfunc[histname_Up + "rebin"] -> GetBinContent(j);
    double current_entry_Down =  mapfunc[histname_Down + "rebin"] ->GetBinContent(j);

    if(current_entry_Up > current_entry_Down){
      map_syst_array[current_histname + cycle_and_sample + systematics + "Up"].push_back(current_entry_Up);
      map_syst_array[current_histname + cycle_and_sample + systematics + "Down"].push_back(current_entry_Down);
    }
    else{
      map_syst_array[current_histname + cycle_and_sample + systematics + "Up"].push_back(current_entry_Down);
      map_syst_array[current_histname + cycle_and_sample + systematics + "Down"].push_back(current_entry_Up);
    }
  }

  cout << "[Save_syst_array] " << current_histname + cycle_and_sample + systematics << "Up size : " << map_syst_array[current_histname + cycle_and_sample + systematics + "Up"].size() << endl;
  cout << "[Save_syst_array] " << current_histname + cycle_and_sample + systematics << "Down size : " << map_syst_array[current_histname + cycle_and_sample + systematics + "Down"].size() << endl;


}

void sum_syst_error(TString current_histname, TString cycle_and_sample, int N_bin){
  
  // -- get central entries
  std::vector<double> central_content;
  central_content.clear();
  
  for(int j = 1; j < N_bin + 1; j++){
    double current_content = mapfunc[current_histname + "_central" + cycle_and_sample + "rebin"] -> GetBinContent(j);
    central_content.push_back(current_content);
  }
  
  cout << "[sum_syst_error] " << endl;
  map_syst_array[current_histname + cycle_and_sample + "Up"].clear();
  map_syst_array[current_histname + cycle_and_sample + "Down"].clear();

  for(int j = 1; j < N_bin + 1; j++){
    double current_error_up = 0.;
    double current_error_down = 0.;

    double current_stat_error = mapfunc[current_histname + "_central" + cycle_and_sample + "rebin"] -> GetBinError(j);
    current_error_up = current_error_up + pow(current_stat_error, 2);
    current_error_down = current_error_down+ pow(current_stat_error, 2);
    
    for(int i_syst = 0; i_syst < N_syst_comparison; i_syst++){
      cout << "[sum_syst_error] for j = " << j << ", i_syst = " << i_syst << endl;
      
      TString current_syst = systematics_comparison[i_syst];
      double current_up = fabs(central_content.at(j-1) - map_syst_array[current_histname + cycle_and_sample + current_syst + "Up"].at(j-1) );
      double current_down = fabs(central_content.at(j-1) - map_syst_array[current_histname + cycle_and_sample + current_syst + "Down"].at(j-1) );
      
      cout << "[sum_syst_error] current_content : " << central_content.at(j-1) << ", up : " << map_syst_array[current_histname + cycle_and_sample + current_syst + "Up"].at(j-1) << ", down : " << map_syst_array[current_histname + cycle_and_sample + current_syst + "Down"].at(j-1) << endl;
      
      current_error_up = current_error_up + pow(current_up, 2);
      current_error_down = current_error_down + pow(current_down, 2);
    }
    
    current_error_up = pow(current_error_up, 0.5);
    current_error_down = pow(current_error_down, 0.5);
    
    map_syst_array[current_histname + cycle_and_sample + "Up"].push_back(current_error_up);
    map_syst_array[current_histname + cycle_and_sample + "Down"].push_back(current_error_down);
  }

  cout << "[sum_syst_error] Up size : " << map_syst_array[current_histname + cycle_and_sample + "Up"].size() << endl;
  cout << "[sum_syst_error] Down size : " << map_syst_array[current_histname + cycle_and_sample + "Down"].size() << endl;

}

void Proper_error_data(TString nameofhistogram, TString current_data, int N_bin, double binx[]){
  
  //vector<double> vx, vy, vexl, vexh, veyl, veyh;
  
  map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] = new TGraphAsymmErrors(GetHist(nameofhistogram + Cycle_name + current_data + "rebin"));
  
  for(int i = 0; i < N_bin; i++){
    int N = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetBinContent(i + 1);
    double L =  (N==0) ? 0  : (ROOT::Math::gamma_quantile(alpha/2,N,1.));
    double U =  (N==0) ? ( ROOT::Math::gamma_quantile_c(alpha,N+1,1) ) : ( ROOT::Math::gamma_quantile_c(alpha/2,N+1,1) );
    if( N!=0 ){
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEYlow(i, (N-L) ); // / (binx[i + 1] - binx[i]) );
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEXlow(i, 0);
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEYhigh(i, (U-N) );// / (binx[i + 1] - binx[i]) );
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEXhigh(i, 0);
      
      double current_x = -1., current_y = -1.;
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> GetPoint(i, current_x, current_y);
      double modified_y = -1.;
      modified_y = current_y; // ( binx[i + 1] - binx[i] );
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPoint(i, current_x, modified_y);

      
      //if(debug) cout << "i : " << i << ", current_x : " << current_x << ", current_y : " << current_y << ", modified_y : " << modified_y << endl;

      veyl.push_back( (N-L)); // / (binx[i + 1] - binx[i]) );
      veyh.push_back( (U-N)); // / (binx[i + 1] - binx[i]));
    }
    else{
      double zerodata_err_low = 0.1 ;// / (binx[i + 1] - binx[i]);
      double zerodata_err_high = 1.8 ;// / (binx[i + 1] - binx[i]);

      veyl.push_back(zerodata_err_low);
      veyh.push_back(zerodata_err_high);
      
      double current_x = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetBinCenter(i + 1);
      double current_ex = GetHist(nameofhistogram + Cycle_name + current_data + "rebin") -> GetBinWidth(i + 1) /2.;

      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEYlow(i, zerodata_err_low);
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEXlow(i, 0.);
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEYhigh(i, zerodata_err_high);
      map_asym_gr[nameofhistogram + Cycle_name + current_data + "correct_error"] -> SetPointEXhigh(i, 0.);

      vx.push_back(current_x);
      vexl.push_back(current_ex);
      vexh.push_back(current_ex);
    }
  }//end for i on N_bin
    
}


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
      400, 800, 1200, 1600, 2000, 2400, 2800, 3200, 3600, 4000, 4400, 4800, 
    };
    return this_masses;
  }
  if(mN==300){
    this_masses = {
      800, 1200, 1600, 2000, 2400, 2800, 3200, 3600, 4000, 4400, 4800,
    };
    return this_masses;
  }
  if(mN==500){
    this_masses = {
      1200, 1600, 2000, 2400, 2800, 3200, 3600, 4000, 4400, 4800,
    };
    return this_masses;
  }
  if(mN==700){
    this_masses = {
      1600, 2000, 2400, 2800, 3200, 3600, 4000, 4400, 4800,
    };
    return this_masses;
  }
  if(mN==900){
    this_masses = {
      2000, 2400, 2800, 3200, 3600, 4000, 4400, 4800,
    };
    return this_masses;
  }
  if(mN==1100){
    this_masses = {
      2400, 2800, 3200, 3600, 4000, 4400, 4800,
    };
    return this_masses;
  }
  if(mN==1300){
    this_masses = {
      2800, 3200, 3600, 4000, 4400, 4800,
    };
    return this_masses;
  }
  if(mN==1500){
    this_masses = {
      3200, 3600, 4000, 4400, 4800,
    };
    return this_masses;
  }
  if(mN==1700){
    this_masses = {
      3600, 4000, 4400, 4800,
    };
    return this_masses;
  }
  if(mN==1900){
    this_masses = {
      4000, 4400, 4800,
    };
    return this_masses;
  }
  if(mN==2100){
    this_masses = {
      4400, 4800,
    };
    return this_masses;
  }
  if(mN==2300){
    this_masses = {
      4800,
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
    if(mZP==800){
      this_masses = {
	100, 300,
      };
      return this_masses;
    }
    if(mZP==1200){
      this_masses = {
	100, 300, 500, 
      };
      return this_masses;
    }
    if(mZP==1600){
      this_masses = {
      100, 300, 500, 700, 
      };
      return this_masses;
    }
    if(mZP==2000){
      this_masses = {
      100, 300, 500, 700, 900, 
      };
      return this_masses;
    }
    if(mZP==2400){
      this_masses = {
      100, 300, 500, 700, 900, 1100, 
      };
      return this_masses;
    }
    if(mZP==2800){
      this_masses = {
      100, 300, 500, 700, 900, 1100, 1300, 
      };
      return this_masses;
    }
    if(mZP==3200){
      this_masses = {
      100, 300, 500, 700, 900, 1100, 1300, 1500, 
      };
      return this_masses;
    }
    if(mZP==3600){
      this_masses = {
	100, 300, 500, 700, 900, 1100, 1300, 1500, 1700,
      };
      return this_masses;
    }
    if(mZP==4000){
      this_masses = {
	100, 300, 500, 700, 900, 1100, 1300, 1500, 1700, 1900,
      };
      return this_masses;
    }
    if(mZP==4400){
      this_masses = {
	100, 300, 500, 700, 900, 1100, 1300, 1500, 1700, 1900, 2100,
      };
      return this_masses;
    }
    if(mZP==4800){
      this_masses = {
        100, 300, 500, 700, 900, 1100, 1300, 1500, 1700, 1900, 2100, 2300,
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



#endif
