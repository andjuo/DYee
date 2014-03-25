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
    if (!this->isInitialized()) printError("HistoPair2D_t(TString)");
  }

  // ----------------

  HistoPair2D_t(const TString &nameBase, const HistoPair2D_t &orig, int setTitle=1) :
    BaseClass_t("HistoPair2D_t"), fHisto(NULL), fHistoSystErr(NULL) {
    fHisto= Clone(orig.fHisto,nameBase,setTitle);
    fHistoSystErr= Clone(orig.fHistoSystErr,nameBase + TString("Syst"),setTitle);
    if (!this->isInitialized()) printError("HistoPair2D_t(HistoPair2D)");
  }

  // ----------------

  void clear() {
    if (fHisto) delete fHisto;
    if (fHistoSystErr) delete fHistoSystErr;
  }

  // ----------------

  TH2D *histo() const { return fHisto; }
  TH2D *histoSystErr() const { return fHistoSystErr; }
  TH2D *editHisto() { return fHisto; }
  TH2D *editHistoSystErr() { return fHistoSystErr; }

  int getNrows() const { return fHisto->GetNbinsX(); }
  int getNcols() const { return fHisto->GetNbinsY(); }

  TString GetName() const { return fHisto->GetName(); }


  // ----------------

  void swapPtrs(HistoPair2D_t &hp) {
    TH2D* htmp=fHisto;  fHisto=hp.fHisto;  hp.fHisto=htmp;
    htmp=fHistoSystErr; fHistoSystErr=hp.fHistoSystErr; hp.fHistoSystErr=htmp;
  }

  // ----------------

  int changeName(TString newName, TString newTitle="") {
    int res=1;
    if (fHisto) { 
      fHisto->SetName(newName);
      if (newTitle.Length()) fHisto->SetTitle(newTitle);
    }
    else res=0;
    if (fHistoSystErr) {
      fHistoSystErr->SetName(newName + TString("Syst"));
      if (newTitle.Length()) fHistoSystErr->SetTitle(newTitle + TString("Syst"));
    }
    return (res) ? 1 : reportError("changeName(%s) - fHisto is null",newName);
  }

  // ----------------

  TH2D* createHistoWithFullError(const TString histoName) const {
    TH2D* h=Clone(fHisto,histoName);
    for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin) {
      for (int jbin=1; jbin<=h->GetNbinsY(); ++jbin) {
	double err=h->GetBinError(ibin,jbin);
	double systErr=fHistoSystErr->GetBinError(ibin,jbin);
	h->SetBinError(ibin,jbin,sqrt(err*err+systErr*systErr));
      }
    }
    return h;
  }

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

  // ----------------

  int chkSameDim(const TH2D* h) const {
    if (!isInitialized()) return reportError("chkSameDim: object is not initialized");
    if (!h) return reportError("chkSameDim: provided histo is NULL");
    if ((fHisto->GetNbinsX() != h->GetNbinsX()) ||
	(fHisto->GetNbinsY() != h->GetNbinsY())) {
      return reportError(Form("chkSameDim: provided histoDims are different: %s[%d,%d], supplied %s[%d,%d]",fHisto->GetName(),fHisto->GetNbinsX(),fHisto->GetNbinsY(),h->GetName(),h->GetNbinsX(),h->GetNbinsY()));
    }
    return 1;
  }

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

  int assign(const HistoPair2D_t &hp) {
    if (&hp == this) return 1; // no self-assignment
    TString hname=fHisto->GetName();
    TString hnameSystErr=fHistoSystErr->GetName();
    delete fHisto;
    delete fHistoSystErr;
    fHisto=Clone(hp.fHisto, hname);
    fHistoSystErr=Clone(hp.fHistoSystErr, hnameSystErr);
    return this->isInitialized();
  }

  // ----------------

  int assign(const TH2D *hBase, const TH2D *hSystError, int resetSystErrCentralVals=1) {
    if (!hBase && !hSystError) return this->reportError("assign(TH2D*,TH2D*): both pointers cannot be 0 at the same time");
    const TH2D* setHBase=(hBase) ? hBase : hSystError;
    const TH2D* setHSyst=(hSystError) ? hSystError : hBase;
    if ((setHBase->GetNbinsX() != setHSyst->GetNbinsX()) ||
	(setHBase->GetNbinsY() != setHSyst->GetNbinsY())) {
      return this->reportError("assign: %s","dim mismatch in the supplied hBase and hSystError");
    }
    TString hname=fHisto->GetName();
    TString hnameSystErr=fHistoSystErr->GetName();
    fHisto=Clone(setHBase, hname);
    fHistoSystErr=Clone(setHSyst, hnameSystErr);

    if (!hBase) fHisto->Reset();

    if (!hSystError) fHistoSystErr->Reset();
    else if (resetSystErrCentralVals) {
      for (int ibin=1; ibin<=fHistoSystErr->GetNbinsX(); ++ibin) {
	for (int jbin=1; jbin<=fHistoSystErr->GetNbinsY(); ++jbin) {
	  fHistoSystErr->SetBinContent(ibin,jbin, 0.);
	}
      }
    }
    return 1;
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

  Bool_t add(const TH2D* h, double weight=1.) { 
    //std::cout << "\nWARN: add(TH2D*) behavior is untested\n";
    return fHisto->Add(h,weight); 
  }

  // ----------------

  Bool_t addSystErr(const TH2D* h, double weight=1.) { 
    /*
    for (int ibin=1; ibin<=fHistoSystErr->GetNbinsX(); ++ibin) {
      for (int jbin=1; jbin<=fHistoSystErr->GetNbinsY(); ++jbin) {
	double err=fHistoSystErr->GetBinError(ibin,jbin);
	double addErr= weight * h->GetBinError(ibin,jbin);
	fHistoSystErr->SetBinError(ibin,jbin, sqrt(err*err + addErr*addErr));
      }
    }
    */
    //std::cout << "add syst err\n";
    return fHistoSystErr->Add(h,weight);
  }

  // ----------------

  Bool_t addSystErrPercent(const TH2D* h, double extra_weight=1.) {
    if (!chkSameDim(h)) return Bool_t(reportError("in function addSystErrPercent"));
    for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin) {
      for (int jbin=1; jbin<=h->GetNbinsY(); ++jbin) {
	double errA= fHistoSystErr->GetBinError(ibin,jbin);
	double errB= fHisto->GetBinContent(ibin,jbin) * h->GetBinError(ibin,jbin) * extra_weight;
	double err=sqrt(errA*errA + errB*errB);
	fHistoSystErr->SetBinError(ibin,jbin, err);
      }
    }
    return true;
  }

  // ----------------

  Bool_t add(const HistoPair2D_t &pair, double weight=1.) {
    return (fHisto->Add(pair.fHisto, weight) &&
	    fHistoSystErr->Add(pair.fHistoSystErr, weight));
  }

  // ----------------

  // divide uncorrelated quantities
  // A = B/C
  // full error:
  // (dA)^2 = 1/C^2 (dB)^2 + B^2/C^4 (dC)^2
  // here we use a decomposition:
  // (dA_stat)^2 = 1/C^2 (dB_stat)^2
  // (dA_syst)^2 = 1/C^2 (dB_syst)^2 + B^2/C^4 (dC_tot)^2
  //             = ( (dB_syst)^2 + A^2 (dC_tot)^2 ) /C^2
  //

  int divide(const HistoPair2D_t &p, const TH2D *h2) {
    //if (!h2) std::cout << "h2 is null" << std::endl; else std::cout << "h2 ok" << std::endl;
    //printHisto(h2);
    TH2D *h2tmp=Clone(h2,h2->GetName() + TString("tmp"),"");
    removeError(h2tmp);
    fHisto->Divide(p.fHisto,h2tmp);
    // evaluate systErr
    for (int ibin=1; ibin<=fHistoSystErr->GetNbinsX(); ++ibin) {
      for (int jbin=1; jbin<=fHistoSystErr->GetNbinsY(); ++jbin) {
	const double term1= p.fHistoSystErr->GetBinError(ibin,jbin);
	const double term2= this->fHisto->GetBinContent(ibin,jbin) * h2->GetBinError(ibin,jbin);
	const double division=h2->GetBinContent(ibin,jbin);
	fHistoSystErr->SetBinContent(ibin,jbin, 0.);
	fHistoSystErr->SetBinError(ibin,jbin, sqrt( term1*term1 + term2*term2 )/division );
      }
    }
    return 1;
  }


  // ----------------

  // if r=a/b,  (dr)^2/r^2 = (da)^2/a^2 + (db)^2/b^2
  // if (dr)^2=(dr)_stat^2 + (dr)_syst^2
  // Statistical and systematical errors are kept separate
  // see also divide_allErrSyst

  int divide(double norm, double normErr, double normSystErr) {
    if (norm==double(0)) {
      return reportError("divide(norm,normErr,normSystErr): norm=0");
    }
    fHisto->Scale(1./norm);
    fHistoSystErr->Scale(1./norm);
    double eSqrStat=normErr*normErr;
    double eSqrSyst=normSystErr*normSystErr;
    if (eSqrStat!=double(0)) {
      // add the extra systematic error
      eSqrStat/=(norm*norm);
      for (int ibin=1; ibin<=fHisto->GetNbinsX(); ++ibin) {
	for (int jbin=1; jbin<=fHisto->GetNbinsY(); ++jbin) {
	  double val = fHisto->GetBinContent(ibin,jbin);
	  double baseStat= fHisto->GetBinError(ibin,jbin);
	  double newErrStatSqr=baseStat*baseStat + val*val*eSqrStat;
	  fHisto->SetBinError(ibin,jbin, sqrt(newErrStatSqr));
	}
      }
    }
    if (eSqrSyst!=double(0)) {
      // add the extra systematic error
      eSqrSyst/=(norm*norm);
      for (int ibin=1; ibin<=fHisto->GetNbinsX(); ++ibin) {
	for (int jbin=1; jbin<=fHisto->GetNbinsY(); ++jbin) {
	  double val = fHisto->GetBinContent(ibin,jbin);
	  double baseSyst= fHistoSystErr->GetBinError(ibin,jbin);
	  double newErrSystSqr=baseSyst*baseSyst + val*val*eSqrSyst;
	  fHistoSystErr->SetBinError(ibin,jbin, sqrt(newErrSystSqr));
	}
      }
    }
    return 1;
  }

  // ----------------

  // if r=a/b,  (dr)^2/r^2 = (da)^2/a^2 + (db)^2/b^2
  // if (dr)^2=(dr)_stat^2 + (dr)_syst^2
  // and (dr)_stat= (da)_stat/b,
  // then (dr)_syst^2 = (da)_syst^2/b^2 + (r/b)^2 (db)^2

  int divide_allErrSyst(double norm, double normErr, double normSystErr) {
    if (norm==double(0)) {
      return reportError("divide_allErrSyst(norm,normErr,normSystErr): norm=0");
    }
    fHisto->Scale(1./norm);
    fHistoSystErr->Scale(1./norm);
    double eSqr=normErr*normErr + normSystErr*normSystErr;
    if (eSqr!=double(0)) {
      // add the extra systematic error
      eSqr/=(norm*norm);
      for (int ibin=1; ibin<=fHisto->GetNbinsX(); ++ibin) {
	for (int jbin=1; jbin<=fHisto->GetNbinsY(); ++jbin) {
	  double val = fHisto->GetBinContent(ibin,jbin);
	  double base= fHistoSystErr->GetBinError(ibin,jbin);
	  double newErrSqr=base*base + val*val*eSqr;
	  fHistoSystErr->SetBinError(ibin,jbin, sqrt(newErrSqr));
	}
      }
    }
    return 1;
  }

  // ----------------

  double ZpeakCount(double *err=NULL) const { return ::ZpeakCount(fHisto,err); }

  double ZpeakCountSystErr() const { 
    double systErr=0.;
    ::ZpeakCount(fHistoSystErr,&systErr);
    return systErr;
  }

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
    if (systErrField && this->fHistoSystErr) {
      fHistoSystErr->SetName(fHisto->GetName() + TString("Syst"));
    }
    if (!res) this->printError("in LoadThreeMatrices");
    return res;
  }

  // ----------------
  
  // assume that a file is open
  int Write(TString saveAsName="") const {
    if (!this->isInitialized()) return 0;
    if (saveAsName.Length()==0) saveAsName=fHisto->GetName();
    fHisto->Write(saveAsName);
    saveAsName.Append("Syst");
    fHistoSystErr->Write(saveAsName);
    return 1;
  }

  // ----------------
  
  // assume that a file is open
  int Read(TString fieldName="", TString fieldErrName="") {
    if (!this->isInitialized()) return 0;
    if (fieldName.Length()==0) fieldName=fHisto->GetName();
    if (fieldErrName.Length()==0) fieldErrName=fHistoSystErr->GetName();
    bool ok=fHisto->Read(fieldName);
    if (ok) ok=fHistoSystErr->Read(fieldErrName);
    if (!ok) { 
      std::cout << "failed to load <" << fieldName << "> and <" 
		<< fieldErrName << "> from "
		<< "current directory <"
		<< gDirectory->GetName() << ">\n";
      return 0;
    }
    if (!this->isInitialized()) return 0;
    fHisto->SetDirectory(0);
    fHistoSystErr->SetDirectory(0);
    return 1;
  }

  // ----------------
  // ----------------

  int Write(TFile &fout, TString subDir="", TString saveWithName="") const {
    if (!this->isInitialized()) {
      return this->reportError("Write(fout=<%s>) object is not initialized",fout.GetName());
    }
    int res=saveHisto(fout,fHisto,subDir,saveWithName);
    if (saveWithName.Length()) saveWithName.Append("_Syst");
    if (res) res=saveHisto(fout,fHistoSystErr,subDir,saveWithName);
    return (res) ? 1 : reportError("Write(fout=<%s>",fout.GetName());
  }

  // ----------------

  int Read(TFile &fin, TString subDir="", TString loadWithName="", int checkBinning=1)  {
    if (!this->isInitialized()) {
      return this->reportError("Read(fin=<%s>) object is not initialized",fin.GetName());      
    }
    TString loadHName=(loadWithName.Length()) ? loadWithName : TString(fHisto->GetName());
    TH2D *h2    =LoadHisto2D(fin,loadHName,subDir,checkBinning);
    TString loadHSystName=(loadWithName.Length()) ? (loadWithName + TString("_Syst")) : TString(fHistoSystErr->GetName());
    TH2D *h2syst=LoadHisto2D(fin,loadHSystName,subDir,checkBinning);
    int res=1;
    if (h2 && h2syst) {
      delete fHisto;
      delete fHistoSystErr;
      fHisto=h2;
      fHistoSystErr=h2syst;
    }
    else res=0;
    return (res) ? 1 : reportError("Read(fin=<%s>",fin.GetName());
  }

  // ----------------
  // ----------------

  int Load(const TString &fname, int checkBinning, TString subDir="") {
    if (!isInitialized()) return this->reportError("Load(%s): object is not initialized",fname);
    TFile file(fname,"read");
    if (!file.IsOpen()) return this->reportError("Load(%s): failed to open file");
    int res=1;
    if (res && checkBinning) res=checkBinningArrays(file);
    if (!res) this->printError("Load(%s): binning error");
    TString fieldName,fieldErrName;
    if (res && subDir.Length()) { 
      if (subDir[subDir.Length()-1]!='/') subDir.Append("/");
      fieldName=subDir + TString(fHisto->GetName());
      fieldErrName=subDir + TString(fHistoSystErr->GetName());
    }
    if (res) res=this->Read(fieldName,fieldErrName);
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
// ---------------------------------------------

HistoPair2D_t* getDiff(const TString &nameBase, 
		       const HistoPair2D_t &var1, 
		       const HistoPair2D_t &var2,
		       int resetSystError) {
  HistoPair2D_t *diff=new HistoPair2D_t(nameBase,var1);
  diff->add(var2,-1);
  if (resetSystError) diff->editHistoSystErr()->Reset();
  return diff;
}

// ---------------------------------------------
// ---------------------------------------------

inline 
TH2D* getRelDifference(const std::vector<HistoPair2D_t*> &hpV, 
		       TString newName, int includeVariants) {
  std::vector<const TH2D*> vec;
  for (unsigned int i=0; i<hpV.size(); ++i) {
    if (hpV[i]) vec.push_back(hpV[i]->histo());
    else std::cout << "getRelDifference(vec<HistoPair*>: null ptr\n";
  }
  TH2D* h2Diff=getRelDifference(vec,newName,includeVariants);
  if (!h2Diff) std::cout << "error in getRelDifference(vec<HistoPair*>)\n";
  return h2Diff;
}

// ---------------------------------------------
// ---------------------------------------------


#endif
