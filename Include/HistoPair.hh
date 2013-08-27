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

  int getNrows() const { return fHisto->GetNbinsX(); }
  int getNcols() const { return fHisto->GetNbinsY(); }

  // ----------------

  //double get(int idx1, int idx2) const { return fHisto->GetBinContent(idx1+1,idx2+1); }
  //double getErr(int idx1, int idx2) const { return fHisto->GetBinError(idx1+1,idx2+1); }
  //double getSystErr(int idx1, int idx2) const { return fHistoSystErr->GetBinError(idx1+1,idx2+1); }
  double getBinContent(int bin1, int bin2) const { return fHisto->GetBinContent(bin1,bin2); }
  double getBinError(int bin1, int bin2) const { return fHisto->GetBinError(bin1,bin2); }
  double getBinSystError(int bin1, int bin2) const { return fHistoSystErr->GetBinError(bin1,bin2); }
  
  void setBinContent(int ibin1, int ibin2, double value) { fHisto->SetBinContent(ibin1,ibin2, value); }
  void setBinError(int ibin1, int ibin2, double value) { fHisto->SetBinError(ibin1,ibin2, value); }
  void setBinSystError(int ibin1, int ibin2, double value) { fHistoSystErr->SetBinError(ibin1,ibin2, value); }

  // ----------------

  int cloneHisto(TH2D* h, int setTitle=1) {
    TString hname=fHisto->GetName();
    if (fHisto) delete fHisto;
    fHisto=Clone(h,hname,setTitle);
    return (fHisto) ? 1:0;
  }

  // ----------------

  int isInitialized() const { return (fHisto && fHistoSystErr) ? 1:0; }

  // ---------------

  TMatrixD *histoAsM() const {
    if (!isInitialized()) { printError("getContentsAsM: object is not initialized"); return NULL; }
    TMatrixD *M=new TMatrixD(fHisto->GetNbinsX(),fHisto->GetNbinsY());
    for (int i=0; i<fHisto->GetNbinsX(); ++i) {
      for (int j=0; j<fHisto->GetNbinsY(); ++j) {
	(*M)(i,j)=fHisto->GetBinContent(i+1,j+1);
      }
    }
    return M;
  }

  // ---------------

  TMatrixD *histoErrAsM() const {
    if (!isInitialized()) { printError("getContentsAsM: object is not initialized"); return NULL; }
    TMatrixD *M=new TMatrixD(fHisto->GetNbinsX(),fHisto->GetNbinsY());
    for (int i=0; i<fHisto->GetNbinsX(); ++i) {
      for (int j=0; j<fHisto->GetNbinsY(); ++j) {
	(*M)(i,j)=fHisto->GetBinError(i+1,j+1);
      }
    }
    return M;
  }

  // ---------------

  TMatrixD *histoSystErrAsM() const {
    if (!isInitialized()) { printError("getContentsAsM: object is not initialized"); return NULL; }
    TMatrixD *M=new TMatrixD(fHistoSystErr->GetNbinsX(),fHistoSystErr->GetNbinsY());
    for (int i=0; i<fHistoSystErr->GetNbinsX(); ++i) {
      for (int j=0; j<fHistoSystErr->GetNbinsY(); ++j) {
	(*M)(i,j)=fHistoSystErr->GetBinError(i+1,j+1);
      }
    }
    return M;
  }

  // ----------------

  int assign(const TMatrixD &val, const TMatrixD &valErr, const TMatrixD &valSystErr) {
    if (!isInitialized()) { return reportError("assign(TMatrixD): object is not initialized"); }
    if ((val.GetNrows() != fHisto->GetNbinsX()) ||
	(val.GetNcols() != fHisto->GetNbinsY())) {
      return reportError("assign(TMatrixD): dim error");
    }
    for (int ibin=1; ibin<=val.GetNrows(); ++ibin) {
      for (int jbin=1; jbin<=val.GetNcols(); ++jbin) {
	fHisto->SetBinContent(ibin,jbin, val(ibin-1,jbin-1));
	fHisto->SetBinError  (ibin,jbin, valErr(ibin-1,jbin-1));
	fHistoSystErr->SetBinContent(ibin,jbin, 0.);
	fHistoSystErr->SetBinError  (ibin,jbin, valSystErr(ibin-1,jbin-1));
      }
    }
    return 1;
  }

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
  int Read() {
    if (!this->isInitialized()) return 0;
    fHisto->Read(fHisto->GetName());
    fHistoSystErr->Read(fHistoSystErr->GetName());
    if (!this->isInitialized()) return 0;
    fHisto->SetDirectory(0);
    fHistoSystErr->SetDirectory(0);
    return 1;
  }

  // ----------------

  int Load(const TString &fname, const TString& subdir="") {
    if (!isInitialized()) return this->reportError("Load(%s): object is not initialized",fname);
    TFile file(fname,"read");
    if (!file.IsOpen()) return this->reportError("Load(%s): failed to open file");
    if (subdir.Length()) file.cd(subdir);
    int res=this->Read();
    file.Close();
    return (res) ? 1 : this->reportError("Load(%s): failed to load <%s>",fname,fHisto->GetName());
  }

  // ----------------

  HistoPair2D_t& operator+= (HistoPair2D_t &h) { this->add(h); return *this; }
  HistoPair2D_t& operator-= (HistoPair2D_t &h) { this->add(h,-1); return *this; }
    
  // ----------------

  int unfold(const TMatrixD &U, const HistoPair2D_t &iniHP,
	     int locUnfoldingBins=DYTools::nUnfoldingBins,
	     const int *locYBins= DYTools::nYBins) {
    int res=1;
    if ((fHisto->GetNbinsX() != iniHP.fHisto->GetNbinsX() ) ||
	(fHisto->GetNbinsY() != iniHP.fHisto->GetNbinsY())) {
      res=reportError("Dim error in unfold(HistoPair2D): (this)[%s], iniHP[%s]",Form("%d,%d", fHisto->GetNbinsX(),fHisto->GetNbinsY()),Form("%d,%d",iniHP.fHisto->GetNbinsX(),iniHP.fHisto->GetNbinsY()));
    }
    if (fHisto->GetNbinsX()*fHisto->GetNbinsY() < locUnfoldingBins) {
      res=reportError("Dim error in unfold(HistoPair2D): (this)[%d,%d], locUnfoldingBins=%d", fHisto->GetNbinsX(),fHisto->GetNbinsY(),locUnfoldingBins);
    }
    if ( (locUnfoldingBins != U.GetNrows()) ||
	 (locUnfoldingBins != U.GetNcols()) ) {
      res=reportError("Dim error in unfold(HistoPair2D): locUnfoldingBins=%d, UnfM[%d,%d]",locUnfoldingBins, U.GetNrows(), U.GetNcols());
    }
    if (res) {
      for (int ir=0, ini=0; ir<fHisto->GetNbinsX(); ir++) {
	if (fHisto->GetNbinsY() < locYBins[ir]) {
	  res=reportError("Detected dim problem in unfold(HistoPair2D): (this)->GetNbinsY=%d, nYBins[ir]=%d",fHisto->GetNbinsY(),locYBins[ir]);
	}
	for (int ic=0; 
	     //(ic<finHPM.GetNcols()) && (ini<locUnfoldingBins); 
	     (ic<fHisto->GetNbinsY()) &&
	       (ic<locYBins[ir]) && (ini<locUnfoldingBins);
	     ++ic, ++ini) {
	  double sum=0;
	  double sumErr2=0;
	  double sumSystErr2=0;
	  for (int jr=0, fin=0; jr<iniHP.fHisto->GetNbinsX(); jr++) {
	    for (int jc=0;
		 //(jc<iniHPM.GetNcols()) && (fin<locUnfoldingBins);
		 (jc<iniHP.fHisto->GetNbinsY()) &&
		   (jc<locYBins[jr]) && (fin<locUnfoldingBins);
		 ++jc, ++fin) {
	      sum += U(fin,ini) * iniHP.fHisto->GetBinContent(jr+1,jc+1);
	      const double tmpErr= U(fin,ini) * iniHP.fHisto->GetBinError(jr+1,jc+1);
	      sumErr2 += tmpErr*tmpErr;
	      const double tmpSystErr= U(fin,ini) * iniHP.fHistoSystErr->GetBinError(jr+1,jc+1);
	      sumSystErr2 += tmpSystErr*tmpSystErr;
	    }
	  }
	  fHisto->SetBinContent(ir+1,ic+1, sum);
	  fHisto->SetBinError(ir+1, ic+1, sqrt(sumErr2));
	  fHistoSystErr->SetBinContent(ir+1, ic+1, 0.);
	  fHistoSystErr->SetBinError(ir+1, ic+1, sqrt(sumSystErr2));
	}
      }
    }
    if (!res) printError("in unfold for %s",fHisto->GetName());
    return res;
  }

  // ----------------

  void print() const {
    printHistoErr(this->fHisto, this->fHistoSystErr);
  }

  // ----------------

};

// ---------------------------------------------


#endif
