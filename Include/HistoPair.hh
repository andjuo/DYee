#ifndef HistoPair_HH
#define HistoPair_HH

#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/BaseClass.hh"


// ---------------------------------------------

class HistoPair2D_t : public BaseClass_t {
protected:
  TH2D *fHisto,*fHistoSystErr;

public:

  HistoPair2D_t() : BaseClass_t("HistoPair2D_t"),
		  fHisto(NULL), fHistoSystErr(NULL) {}

  // ----------------

  HistoPair2D_t(const TString &nameBase) :
    BaseClass_t("HistoPair2D_t"), fHisto(NULL), fHistoSystErr(NULL) {
    fHisto=createBaseH2(nameBase);
    fHistoSystErr=createBaseH2(nameBase + TString("Syst"));
  }

  // ----------------

  TH2D *histo() const { return fHisto; }
  TH2D *histoSystErr() const { return fHistoSystErr; }
  TH2D *editHisto() { return fHisto; }
  TH2D *editHistoSystErr() { return fHistoSystErr; }

  // ----------------

  int cloneHisto(TH2D* h, int setTitle=1) {
    TString hname=fHisto->GetName();
    if (fHisto) delete fHisto;
    fHisto=Clone(h,hname,setTitle);
    return (fHisto) ? 1:0;
  }

  // ----------------

  int isInitialized() const { return (fHisto && fHistoSystErr) ? 1:0; }

  // ----------------

  Bool_t add(TH2D* h, double weight=1.) { return fHisto->Add(h,weight); }
  Bool_t addSystErr(TH2D* h, double weight=1.) { return fHistoSystErr->Add(h,weight); }

  // ----------------

  Bool_t add(HistoPair2D_t &pair, double weight=1.) {
    return (fHisto->Add(pair.fHisto, weight) &&
	    fHistoSystErr->Add(pair.fHistoSystErr, weight));
  }

  // ----------------

  double ZpeakCount(double *err=NULL) { return ::ZpeakCount(fHisto,err); }

  // ----------------

  int correctNegativeValues() {
    int count=0;
    for (int i=1; i<=fHisto->GetNbinsX(); ++i) {
      for (int j=1; j<=fHisto->GetNbinsY(); ++j) {
	if (fHisto->GetBinContent(i,j)<0) {
	  count++;
	  fHisto->SetBinContent(i,j, 0.);
	}
      }
    }
    return count;
  }

  // ----------------

  int loadThreeMatrices(const TString &fname,
			const TString &field, const TString &errField, 
			const TString &systErrField,
			int checkBinning=1, int absoluteRapidity=1) {
    if (fHisto) delete fHisto;
    if (fHistoSystErr) delete fHistoSystErr;
    int res=LoadThreeMatrices(fname,
			      &this->fHisto, &this->fHistoSystErr,
			      field,errField,systErrField,
			      checkBinning,absoluteRapidity);
    if (!res) this->printError("in LoadThreeMatrices");
    return res;
  }

  // ----------------
  
  // assume that a file is open
  int Write() const {
    if (!this->isInitialized()) return 0;
    fHisto->Write();
    fHistoSystErr->Write();
    return 1;
  }

  // ----------------
  
  // assume that a file is open
  int Read() const {
    if (!this->isInitialized()) return 0;
    fHisto->Read(fHisto->GetName());
    fHistoSystErr->Read(fHistoSystErr->GetName());
    if (!this->isInitialized()) return 0;
    fHisto->SetDirectory(0);
    fHistoSystErr->SetDirectory(0);
    return 1;
  }

  // ----------------

  void print() const {
    printHistoErr(this->fHisto, this->fHistoSystErr);
  }

  // ----------------

};

// ---------------------------------------------


#endif
