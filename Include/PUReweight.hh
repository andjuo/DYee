#ifndef PUReweight_HH
#define PUReweight_HH

#include <TROOT.h>
#include <TString.h>
#include <TFile.h>
#include <TH1D.h>
#include <iostream>
#include <assert.h>
#include "../Include/DYTools.hh"
//#include "../Include/MyTools.hh"


class PUReweight_t {
public:
  typedef enum { maxPVs=65 } TConst_t;
  typedef enum { _none, _Hildreth, 
		 _Hildreth_plus5percent, _Hildreth_minus5percent,
		 _TwoHistos } 
    TReweightMethod_t;
protected:
  TString FName; // file name
  TFile *FFile; // pointer to a file
  TH1D *hRef;   // reference histogram -- what we want to have
  TH1D *hActive; // active histogram (source distribution that needs reshaping)
  TH1D *hWeight; // histogram of weights
  TH1D *hWeightHildreth; // histogram of weights according to the Hildreth's method
  int FCreate; // whether a file is being created (1) or updated (2), otherwise - reading (0)
  TReweightMethod_t FActiveMethod;
public:
  PUReweight_t(TReweightMethod_t method=_Hildreth);
  PUReweight_t(DYTools::TSystematicsStudy_t systMode, TReweightMethod_t method=_Hildreth);
  ~PUReweight_t() { this->clear(); }

  void clear() {
    if (hRef) { delete hRef; hRef=0; }
    if ((FCreate!=0) && FFile && hActive) { FFile->cd(); hActive->Write(); }
    if (hActive) { delete hActive; hActive=0; }
    if (hWeight) { delete hWeight; hWeight=0; }
    if (FFile) { delete FFile; FFile=0; }
    FCreate=0;
  }

  // access
  const TString& fileName() const { return FName; }
  const TH1D* getHRef() const { return hRef; }
  const TH1D* getHActive() const { return hActive; }
  const TH1D* getHWeigth() const { return hWeight; }
  int getCreate() const { return FCreate; }

  int setActiveMethod(TReweightMethod_t method) {
    int res=1;
    switch(method) {
    case _none: ; break;
    case _Hildreth: 
    case _Hildreth_plus5percent:
    case _Hildreth_minus5percent:
      if (!hWeightHildreth) assert(initializeHildrethWeights(method));
      break;
    case _TwoHistos: ; break; // RecoLevel initialization cannot be done here
    default:
      std::cout << "PUReweight::setActiveMethod is not ready for method=<"
		<< method
		<< ">\n";
      res=0;
    }
    FActiveMethod=method;
    return res;
  }

  // The generic getWeight is commented out to force the users
  // make explicit choice of which exactly method is called.
  //
  double getWeight(int nPV) const {
    double weight=0.;
    switch(FActiveMethod) {
    case _none: assert(0); break;
    case _Hildreth_plus5percent:
    case _Hildreth_minus5percent:
    case _Hildreth: weight=getWeightHildreth(nPV); break;
    case _TwoHistos: weight=getWeightTwoHistos(nPV); break;
    default:
      std::cout << "PUReweight::getWeight is not ready for method=<"
 		<< FActiveMethod
 		<< ">\n";
    }
    return weight;
  }

  double getWeightTwoHistos(int nGoodPV) const {
    double w=0;
    if (!hWeight) {
      std::cout << "PUReweight::getWeightTwoHistos: call setActiveSample first\n";
    }
    else {
      if (nGoodPV> PUReweight_t::maxPVs+1) nGoodPV=PUReweight_t::maxPVs+1;
      int idx=hWeight->FindBin( nGoodPV );
      w=hWeight->GetBinContent( idx );
    }
    return w;
  }

  double getWeightHildreth(int nPU) const {
    double w=0;
    if (!hWeightHildreth) {
      std::cout << " The weights for Hildreth's method not available\n";
      std::cout << " A problem during intialization of PUReweight?\n";
    }
    else {
      if (nPU > hWeightHildreth->GetNbinsX() ) 
	nPU = hWeightHildreth->GetNbinsX(); 

      int idx = hWeightHildreth->FindBin( nPU );
      w = hWeightHildreth->GetBinContent( idx );
    }
    return w;
  }

  int setHildrethWeights(TReweightMethod_t method=_Hildreth) { return initializeHildrethWeights(method); }

  // setup weights from two histograms: weights=targetHisto/sourceHisto
  int setSimpleWeights(const TString &targetFile, 
		       const TString &targetHistoName,
		       const TString &sourceFile, 
		       const TString &sourceHistoName,
		       int loadTH1F=1);
  
  int setDefaultFile(const TString &use_dirTag, const TString &analysisTag, 
		     int create=0) {
    TString fname=TString("../root_files/selected_events/") + use_dirTag + 
      TString("/npv") + analysisTag + TString(".root");
    int res=setFile(fname,create);
    if (!res) std::cout << "Error in PUReweight::setDefaultFile\n";
    return res;
  }

  int setFile(const TString &fname, int create=0);
  int setReference(const TString &setName);
  int setReference(const TH1D *new_hRef) {
    hRef=(TH1D*)new_hRef->Clone(new_hRef->GetName()+TString("_clone"));
    hRef->SetDirectory(0);
    return 1; 
  }
  int setActiveSample(const TString &setName); // set active sample and calculate weights
  int prepareWeights(int save_weights); // needs to be called if setReference(TH1D) was called

  int Fill(UInt_t nGoodPV, double weight) {
    if ((FCreate==0) || !hActive) {
      std::cout << "cannot fill the histogram:\n";
      if (FCreate==0) std::cout << " - file is not opened for creation\n";
      if (!hActive) std::cout << " - active histogram is not set\n";
      return 0;
    }
    if (nGoodPV> PUReweight_t::maxPVs+1) nGoodPV=PUReweight_t::maxPVs+1;
    hActive->Fill(nGoodPV,weight);
    return 1;
  }

  int printActiveDistr_and_Weights(std::ostream& out=std::cout) const;
  void print(std::ostream& out=std::cout) const;

  int printWeights(std::ostream &out) const { 
    int res=this->printHisto(out,hWeight,"hWeight");
    if (!res) out << "in printWeights\n";
    return res;
  }

  int printActive(std::ostream &out) const { 
    int res=this->printHisto(out,hActive,"hActive");
    if (!res) out << "in printActive\n";
    return res;
  }

protected:
  TH1D *newHisto(const TString &name) const {
    // histogram: PUs from 0 to maxPVs with an overflow bin (maxPVs+1) 
    TH1D *h= new TH1D(name,name,PUReweight_t::maxPVs+2,-0.5,Double_t(PUReweight_t::maxPVs+1.5));
    h->Sumw2();
    h->GetXaxis()->SetTitle("nGoodPVs"); h->GetYaxis()->SetTitle("weight (a.u.)");
    return h;
  }

  int printHisto(std::ostream& out, const TH1D* histo, const TString &name) const;

  int initializeHildrethWeights(TReweightMethod_t method=_Hildreth);

  // weights=target/source
  int initializeTwoHistoWeights(TH1D* hTarget, TH1D* hSource);
};


#endif
