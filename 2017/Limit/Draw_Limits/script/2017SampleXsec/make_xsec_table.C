#include <TH2.h>
#include <TH1F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TFile.h>
#include <TChain.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TKey.h>
#include <TList.h>
#include <TObject.h>
#include <TCollection.h>
#include <TLatex.h>
#include <Math/DistFunc.h>

using namespace std;

map<TString, TH1*> maphist;
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


///////////////////
// -- Main Function
///////////////////
void plot(){
  
  ofstream Xsec_output;
  Xsec_output.open("values_MG_NLO.txt");
  
  ifstream Xsec_input;
  Xsec_input.open("Xsec.txt");
  
  char line[500];
  if(Xsec_input.is_open()){
    while(!Xsec_input.eof()){
      Xsec_input.getline(line, 500);
      
      TString this_line = line;
      if(this_line.Contains("ZprimetoNN")){
        TObjArray *tx = this_line.Tokenize("\t");
        TString Sample_name = ((TObjString *)(tx->At(0)))->String();
        TString Xsec = ((TObjString *)(tx->At(1)))->String();

	cout << "Sample_name : " << Sample_name << ", Xsec : " << Xsec << endl;
	
	TObjArray *tx_sample = Sample_name.Tokenize("_");
	TString mZp = ((TObjString *)(tx_sample->At(2)))->String();
	TString mN = ((TObjString *)(tx_sample->At(3)))->String();
	
	mZp.Remove(0,6);
	mN.Remove(0,1);
	
	cout << "mZp : " << mZp << ", mN : " << mN << endl;
	
	Xsec_output << mZp << "\t" << mN << "\t" << Xsec << endl;
	
      }
    }
  }
  
  
  Xsec_input.close();
  Xsec_output.close(); 
  
}
