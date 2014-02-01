#ifndef CovariantMatrix_HH
#define CovariantMatrix_HH

#include "../Include/BaseClass.hh"
#include "CovariantTensor.hh"

// ------------------------------------------------------------
// ------------------------------------------------------------

struct SaveFigParams_t {
  int canvWidth,canvHeight;
  TColorRange_t centerRange;
  int massBins, nColorBins;
  double maxValUser;
  int do_delete;
  
public:
  SaveFigParams_t() : canvWidth(600),canvHeight(600),
		      centerRange(_colrange_center),
		      massBins(1),
		      nColorBins(51),
		      maxValUser(0.),
		      do_delete(1) 
  { }

  SaveFigParams_t(TColorRange_t set_range, int set_massBins=1,
		  int set_nColorBins=21, double set_maxValUser=0.,
		  int set_do_delete=1) : 
    canvWidth(600),canvHeight(600),
    centerRange(set_range),
    massBins(set_massBins),
    nColorBins(set_nColorBins),
    maxValUser(set_maxValUser),
    do_delete(set_do_delete) 
  { }

  void SetCanvasSize(int width, int height) { canvWidth=width; canvHeight=height; }
  void SetCenterRange(TColorRange_t range) { centerRange=range; }
  void SetMassBins(int val) { massBins=val; }
  void SetNColorBins(int val) { nColorBins=val; }
  void SetMaxValUser(double val) { maxValUser=val; }
  void SetDeletion(int val) { do_delete=val; }
};

// ------------------------------------------------------------
// ------------------------------------------------------------

class CovarianceMatrix_t : public BaseClass_t {
protected:
  TMatrixD CovM;
  TString FName;
public:
  CovarianceMatrix_t(const TString &name, int dim) : 
    BaseClass_t("CovarianceMatrix_t"),
    CovM(dim,dim),
    FName(name)
  {}

  CovarianceMatrix_t(const TString &name, const CovarianceMatrix_t &M) :
    BaseClass_t("CovarianceMatrix_t"),
    CovM(M.CovM),
    FName(name)
  {}

  CovarianceMatrix_t(const TString &name, const TMatrixD& M) :
    BaseClass_t("CovarianceMatrix_t"),
    CovM(M),
    FName(name)
  {}

  CovarianceMatrix_t(const TString &name,
		     TMatrixT<double>::EMatrixCreatorsOp1 op,
		     const CovarianceMatrix_t &prototype) :
    BaseClass_t("CovarianceMatrix_t"),

    CovM(op,prototype.CovM),
    FName(name) 
  {}

  template<class MatrixType_t>
  CovarianceMatrix_t(const TString &name,
		     TMatrixT<double>::EMatrixCreatorsOp1 op,
		     const MatrixType_t &prototype) :
    BaseClass_t("CovarianceMatrix_t"),
    CovM(op,prototype),
    FName(name) 
  {}

  // --------------------------------

  const TMatrixD& GetCovM() const { return CovM; }
  TMatrixD& EditCovM() { return CovM; }
  const TString& GetName() const { return FName; }
  TString &EditName() { return FName; }
  void SetName(const TString &new_name) { FName=new_name; }
  double operator()(int i, int j) const { return CovM(i,j); }
  
  void Zero() { CovM.Zero(); }

  // ---------------

  void SetValue(int i, int j, double val) { CovM(i,j) =  val; }
  void AddValue(int i, int j, double val) { CovM(i,j) += val; }

  // ---------------

  int Assign(const CovarianceMatrix_t &CM) {
    CovM=CM.CovM;
    return 1;
  }

  // ---------------

  double ResetElement(int i, int j, double newVal=0.) { 
    double val=CovM(i,j); 
    CovM(i,j)=newVal; 
    return val; 
  }

  // ---------------

  int SetDiagonalElements(const TVectorD &vec) {
    CovM.Zero();
    TMatrixDDiag diag(CovM);
    for (int i=0; i<vec.GetNoElements(); ++i) {
      diag(i)=vec[i];
    }
    return 1;
  }
  // ---------------

  int SetDiagonalElementsSqr(const TVectorD &vec) {
    CovM.Zero();
    TMatrixDDiag diag(CovM);
    for (int i=0; i<vec.GetNoElements(); ++i) {
      diag(i)=vec[i]*vec[i];
    }
    return 1;
  }

  // ---------------

  double checkDiagonal(const TMatrixD &referenceMatrix, int *storeIdx=NULL) const {
    double max=0., maxAbs=-1.;
    int idx=-1;
    for (int ir=0; ir<CovM.GetNrows(); ++ir) {
      double diff= CovM(ir,ir) - referenceMatrix(ir,ir);
      if (fabs(diff)>maxAbs) {
	max=diff;
	maxAbs=fabs(diff);
	idx=ir;
      }
    }
    if (storeIdx) *storeIdx=idx;
    return max;
  }

  // ---------------

  double checkDiagonal(const CovarianceMatrix_t &CM, int *storeIdx=NULL) const {
    return checkDiagonal(CM.CovM,storeIdx);
  }

  // ---------------

  double checkDiagonal(const TVectorD &referenceVector, int *storeIdx=NULL) const {
    double max=0., maxAbs=-1.;
    int idx=-1;
    for (int ir=0; ir<CovM.GetNrows(); ++ir) {
      double diff= CovM(ir,ir) - referenceVector(ir);
      if (fabs(diff)>maxAbs) {
	max=diff;
	maxAbs=fabs(diff);
	idx=ir;
      }
    }
    if (storeIdx) *storeIdx=idx;
    return max;
  }

  // ---------------

  double checkDiagonalSqrt(const TMatrixD &referenceMatrix, int *storeIdx=NULL) const {
    double max=0., maxAbs=-1.;
    int idx=-1;
    for (int ir=0; ir<CovM.GetNrows(); ++ir) {
      double diff= sqrt(CovM(ir,ir)) - referenceMatrix(ir,ir);
      if (fabs(diff)>maxAbs) {
	max=diff;
	maxAbs=fabs(diff);
	idx=ir;
      }
    }
    if (storeIdx) *storeIdx=idx;
    return max;
  }

  // ---------------

  double checkDiagonalSqrt(const TVectorD &referenceVector, int *storeIdx=NULL) const {
    double max=0., maxAbs=-1.;
    int idx=-1;
    for (int ir=0; ir<CovM.GetNrows(); ++ir) {
      double diff= sqrt(CovM(ir,ir)) - referenceVector(ir);
      if (fabs(diff)>maxAbs) {
	max=diff;
	maxAbs=fabs(diff);
	idx=ir;
      }
    }
    if (storeIdx) *storeIdx=idx;
    return max;
  }

  // --------------------------------

  int Unsquare() {
    int ok=1;
    for (int ir=0; ok && (ir<CovM.GetNrows()); ++ir) {
      TMatrixDRow_const row(CovM,ir);
      for (int ic=0; ok && (ic<CovM.GetNcols()); ++ic) {
	if (row[ic]<0) ok=0;
      }
    }
    if (ok) {
      for (int ir=0; ok && (ir<CovM.GetNrows()); ++ir) {
	TMatrixDRow row(CovM,ir);
	for (int ic=0; ok && (ic<CovM.GetNcols()); ++ic) {
	  row(ic) = sqrt(row[ic]);
	}
      }
    }
    return ok;
  }

  // --------------------------------

  int Add(const CovarianceMatrix_t &CM, double scale=1.) {
    if (scale==1) CovM += CM.CovM;
    else CovM += scale *CM.CovM;
    return 1;
  }

  // --------------------------------

  int Add(const TMatrixD &M, double scale=1.) {
    if (scale!=1.0) CovM += scale*M;
    else CovM += M;
    return 1;
  }

  // --------------------------------

  int AddDiagonalElements(const CovarianceMatrix_t &CM) {
    //std::cout << "CovM:\n"; CovM.Print(); 
    //std::cout << "CM.CovM:\n"; CM.CovM.Print();
    for (int ir=0; ir<CovM.GetNrows(); ++ir) {
      TMatrixDRow(CovM,ir)[ir] += TMatrixDRow_const(CM.CovM,ir)[ir];
    }
    //std::cout << "CovM += diag(CM.CovM):\n"; CovM.Print();
    return 1;
  }

  // --------------------------------

  int SetDiagonalElements(const CovarianceMatrix_t &CM) {
    //std::cout << "CovM:\n"; CovM.Print(); 
    //std::cout << "CM.CovM:\n"; CM.CovM.Print();
    for (int ir=0; ir<CovM.GetNrows(); ++ir) {
      TMatrixDRow(CovM,ir)[ir] = TMatrixDRow_const(CM.CovM,ir)[ir];
    }
    //std::cout << "CovM += diag(CM.CovM):\n"; CovM.Print();
    return 1;
  }

  // --------------------------------

  int GetDiagonalElements(TVectorD &Vsqr) const {
    if (CovM.GetNrows()!=Vsqr.GetNoElements()) 
      return reportError("[%s]GetDiagonalElements(TVectorD): dims %d vs %d",FName.Data(),CovM.GetNrows(),Vsqr.GetNoElements());

    for (int ir=0; ir<CovM.GetNrows(); ++ir) {
      Vsqr[ir]=TMatrixDDiag_const(CovM)[ir];
    }
    return 1;
  }

  // --------------------------------

  int GetDiagonalError(TVectorD &Verr) const {
    if (!GetDiagonalElements(Verr)) return reportError(" in GetDiagonalError");
    for (int i=0; i<Verr.GetNoElements(); ++i) {
      Verr[i] = sqrt(Verr[i]);
    }
    return 1;
  }

  // --------------------------------

  int KeepOnlyDiagonalElements() {
    CovarianceMatrix_t tmp("tmp",*this);
    this->CovM.Zero();
    this->SetDiagonalElements(tmp);
    return 1;
  }

  // --------------------------------

  // propagate:  M' = U M U^T

  int Propagate(const TMatrixD &U, int transpose=0) {
    TMatrixD Utrans(TMatrixD::kTransposed,U);
    if (!transpose) { // U M U^T
      TMatrixD tmp(U);
      tmp *= CovM;
      CovM .Mult(tmp,Utrans);
    }
    else { // U^T M U
      // following UnfoldingTools.C prescription
      TMatrixD tmp(Utrans);
      tmp *= CovM;
      CovM .Mult(tmp,U);
    }
    return 1;
  }

  // --------------------------------

  void Scale(double scale) {
    CovM *= scale;
  }

  // --------------------------------

  double Invert(double extraScale=1.) {
    double det=-1;
    if (extraScale!=1.) CovM *= extraScale;
    CovM.Invert(&det);
    std::cout << this->GetName() << ".Invert  det=" << det << ", extraScale=" << extraScale << "\n";
    if (extraScale!=1.) CovM *= extraScale;
    return det;
  }

  // --------------------------------

  // calculate element_ij = V_i V_j covM_ij

  int ScaleByVec(const TVectorD &multV) {
    int res=1;
    for (int ir=0; ir<CovM.GetNrows(); ++ir) {
      TMatrixDRow row(CovM,ir);
      row *= multV[ir];
      TMatrixDColumn col(CovM,ir);
      col *= multV[ir];
    }
    if (!res) this->printError("ScaleByVec");
    return res;
  }
  
  // --------------------------------

  // calculate element_ij = 1/(V_i V_j) covM_ij

  int InvScaleByVec(const TVectorD &multV) {
    int res=1;
    for (int ir=0; ir<CovM.GetNrows(); ++ir) {
      TMatrixDRow row(CovM,ir);
      row *= 1./multV[ir];
      TMatrixDColumn col(CovM,ir);
      col *= 1./multV[ir];
    }
    if (!res) this->printError("InvScaleByVec");
    return res;
  }
  
  // --------------------------------

  // calculate element_ij = Sum(m,n) V_m V_n covT(im,jn)

  int Contract(const TVectorD &multV,
	       const CovarianceTensor_t &covT,
	       int onlyDiagonalTerms=0) {
    CovM.Zero();
    if (covT.GetDim()==0) {
      //this->printWarning(Form("Contract [Name=%s]: covT.GetDim=0\n",FName.Data()));
      return 1;
    }
    int res=1;
    for (int ir=0; ir<CovM.GetNrows(); ++ir) {
      TMatrixDRow row(CovM,ir);
      for (int ic=0; ic<CovM.GetNcols(); ++ic) {
	double sum=0.;
	for (int m=0; m<covT.GetNrows(); ++m) {
	  const double w1=multV[m];
	  const int idx1=covT.FlatIdx(ir,m);
	  TMatrixDRow_const covTrow(covT.GetCovTensor(), idx1);
	  for (int n=0; n<covT.GetNcols(); ++n) {
	    if (onlyDiagonalTerms && ((m!=n) || (ir!=ic))) continue; // only diagonal covariance
	    const int idx2=covT.FlatIdx(ic,n);
	    //sum+= w1*multV[n] * covT(idx1,idx2);
	    sum+= w1*multV[n] * covTrow[idx2];
	  }
	}
	row(ic)=sum;
      }
    }
    if (!res) this->printError("contract");
    return res;
  }
  
  // --------------------------------
  
  // calculate element_ij = 1/(W_i W_j) Sum(m,n) V_m V_n covT(im,jn)

  int ContractAndInvWeight(const TVectorD &multV, 
			   const CovarianceTensor_t &covT,
			   const TVectorD &weightV) {
    int res=
      Contract(multV,covT) &&
      InvScaleByVec(weightV);
    if (!res) this->printError("contractAndInvWeight");
    return res;
  }

  // --------------------------------

  // calculate correlation matrix
  int Correlations(const CovarianceMatrix_t &CM) {
    for (int ir=0; ir<CovM.GetNrows(); ++ir) {
      for (int ic=0; ic<CovM.GetNcols(); ++ic) {
	this->CovM(ir,ic) = CM.CovM(ir,ic)/sqrt(CM.CovM(ir,ir)*CM.CovM(ic,ic));
      }
    }
    return 1;
  }

  // --------------------------------

  TH2D* DrawSubpad(TCanvas *c, int subPad, TColorRange_t centerRange, int massBins, int nColorBins=51, double maxValUser=0.) const {
    TString name=TString("h2D_") + FName;
    TString title=name;
    TH2D* h2=createHisto2D(CovM, NULL, name.Data(), title.Data(), centerRange, massBins,maxValUser);
    if (c && h2) drawHistoSubpad(c,subPad,h2,centerRange,massBins,nColorBins);
    else std::cout << "CovariantMatrix(" << FName << ")::DrawSubpad: either canvas or histogram is null. Nothing drawn\n";
    return h2;
  }

  // --------------------------------

  TH2D* Draw(TCanvas *c, TColorRange_t centerRange, int massBins, int nColorBins=51, double maxValUser=0.) const {
    return DrawSubpad(c,0,centerRange,massBins,nColorBins,maxValUser);
  }

  // --------------------------------

  TH2D* DrawCorrSubpad(TCanvas *c, int subPad, TColorRange_t centerRange, int massBins, int nColorBins=51, double maxValUser=0.) const {
    TString newName=this->FName + TString("_corr");
    CovarianceMatrix_t Mcorr(newName,*this);
    Mcorr.Correlations(*this);
    return Mcorr.DrawSubpad(c,subPad,centerRange,massBins,nColorBins,maxValUser);
  }

  // --------------------------------

  TH2D* DrawCorr(TCanvas *c, TColorRange_t centerRange, int massBins, int nColorBins=51, double maxValUser=0.) const {
    return DrawCorrSubpad(c,0,centerRange,massBins,nColorBins,maxValUser);
  }

  // --------------------------------

  int SaveFig(const TString &figFileName, int canvWidth=600, int canvHeight=600,
	      TColorRange_t centerRange=_colrange_center, 
	      int massBins=1, int nColorBins=51,
	      double maxValUser=0., int do_delete=1) const {
    int index=figFileName.Index('.');
    TString canvName=TString("canv_") + 
      ((index>0) ? figFileName(0,index) : figFileName);
    TCanvas *c=new TCanvas(canvName,canvName, canvWidth,canvHeight);
    AdjustFor2DplotWithHeight(c);
    TH2D *h=this->Draw(c,centerRange,massBins,nColorBins,maxValUser);
    h->SetDirectory(0);
    h->Draw("colz");
    c->Update();
    c->SaveAs(figFileName);
    if (do_delete) {
      delete c;
      delete h;
    }
    return 1;
  }

  // --------------------------------

  int SaveFig(const TString &figFileName, const SaveFigParams_t &p) const {
    return SaveFig(figFileName,
		   p.canvWidth,p.canvHeight,
		   p.centerRange,
		   p.massBins,
		   p.nColorBins,
		   p.maxValUser,
		   p.do_delete);
  }

  // --------------------------------

  int SaveFig(const TString &figSavePath, const TString &figFileNameTag, const SaveFigParams_t &p) const {
    TString figFileName= figSavePath + 
      TString("fig_") + this->FName + figFileNameTag + TString(".png");
    return SaveFig(figFileName,
		   p.canvWidth,p.canvHeight,
		   p.centerRange,
		   p.massBins,
		   p.nColorBins,
		   p.maxValUser,
		   p.do_delete);
  }

  // --------------------------------

  int SaveCorrFig(const TString &figSavePath, const TString &figFileNameTag, const SaveFigParams_t &p) const {
    CovarianceMatrix_t C(this->FName + TString("_corr"),this->CovM);
    C.Correlations(*this);
    TString saveName= figSavePath + 
      TString("fig_") + C.GetName() + figFileNameTag + TString(".png");
    return C.SaveFig(saveName,
		     p.canvWidth,p.canvHeight,
		     p.centerRange,
		     p.massBins,
		     p.nColorBins,
		     p.maxValUser,
		     p.do_delete);
  }

  // --------------------------------

  template<class Histo1D_t>
  int FillError(Histo1D_t *h) const {
    if (!h) return reportError("FillError: null histo");
    if (CovM.GetNrows()!=h->GetNbinsX()) return reportError("FillError: CovM.GetNrows=%d, h->GetNbinsX=%d",CovM.GetNrows(),h->GetNbinsX());
    for (int i=0; i<CovM.GetNrows(); ++i) {
      h->SetBinContent(i+1, sqrt(TMatrixDDiag_const(CovM)[i]));
      h->SetBinError  (i+1, 0.);
    }
    return 1;
  }

  // --------------------------------

  int Save(TFile &fout, const char *saveName=NULL) const {
    if (!fout.IsOpen()) return 0;
    fout.cd();
    TString name=(saveName) ? TString(saveName) : FName;
    CovM.Write(name);
    return 1;
  }

  // --------------------------------

  int Load(TFile &fin, const char *newName=NULL) {
    if (!fin.IsOpen()) return 0;
    if (newName) FName=newName;
    TMatrixD* tmp=(TMatrixD*) fin.Get(FName);
    if (!tmp) return 0;
    CovM = *tmp;
    return 1;
  }


  // --------------------------------

  int Write(TFile &fout, const char *saveName=NULL) const { return this->Save(fout,saveName); }
  int Read(TFile &fin, const char *newName=NULL) { return this->Load(fin,newName); }

  // --------------------------------

  int Load(const TString &fname, const char *newName=NULL) {
    TFile fin(fname,"read");
    int res=this->Load(fin,newName);
    fin.Close();
    return res;
  }

  // --------------------------------

  int Save(const TString &fname, const char *oldName=NULL) const {
    TFile fout(fname,"recreate");
    int res=this->Save(fout,oldName);
    fout.Close();
    return res;
  }

  // --------------------------------

  void Print() const {
    std::cout << "CovarianceMatrix_t(" << FName << "): \n";
    CovM.Print();
    std::cout << std::endl;
  }

  // --------------------------------

  void PrintCorr() const {
    CovarianceMatrix_t C(this->FName + TString("_corr"),this->CovM);
    C.Correlations(*this);
    std::cout << "Correlations in ";
    C.Print();
  }

  // --------------------------------

  void PrintDiag(const char *extra_msg, const char *number_format=NULL, int take_sqrt=0) const {
    TString format= 
      Form(" %c2d %s\n", '%', 
	   (number_format==NULL) ? "%6.2g" : number_format);
    if (extra_msg) std::cout << extra_msg << "\n";
    std::cout << "diagonal elements of CovarianceMatrix_t(" << FName << "): \n";
    if (!take_sqrt) {
      for (int i=0; i<CovM.GetNrows(); ++i) {
	std::cout << Form(format.Data(), i,CovM(i,i));
      }
    }
    else {
      for (int i=0; i<CovM.GetNrows(); ++i) {
	std::cout << Form(format.Data(), i,sqrt(CovM(i,i)));
      }
    }
    std::cout << std::endl;
  }

  // --------------------------------
  // --------------------------------
};

// ------------------------------------------------------------
// ------------------------------------------------------------

inline
double calculateZPeakCS(const TVectorD &xsec, int idxMin, int idxMax) {
  double sum=0.;
  for (int i=idxMin; i<idxMax; ++i) {
    sum+= xsec[i];
  }
  return sum;
}

// ------------------------------------------------------------

inline
double calculateZPeakCSErr(const TVectorD &xsecTotErr, int idxMin, int idxMax) {
  double sum=0.;
  for (int i=idxMin; i<idxMax; ++i) {
    sum+= pow(xsecTotErr[i],2);
  }
  return sqrt(sum);
}

// ------------------------------------------------------------

inline
int calculateZPeakCov(const CovarianceMatrix_t &covM, int idxMin, int idxMax,
		      double &covSigmaZ_val, TVectorD &covSigmaZSigma_vec,
		      int printDebug=0) {

  covSigmaZ_val=0.;
  double diagSum=0.0;
  for (int i=idxMin; i<idxMax; ++i) {
    for (int j=idxMin; j<idxMax; ++j) {
      covSigmaZ_val  += covM(i,j);
      if (i==j) diagSum += covM(i,j);
    }
  }
  if (printDebug) std::cout << "covSigmaZ_val=" << covSigmaZ_val << ", diag contributions=" << diagSum << "=(" << sqrt(diagSum) << ")^2\n";
  covSigmaZSigma_vec.Zero();
  for (int i=0; i<covSigmaZSigma_vec.GetNoElements(); ++i) {
    double sum=0;
    for (int j=idxMin; j<idxMax; ++j) {
      sum+= covM(i,j);
    }
    covSigmaZSigma_vec(i)=sum;
  }
  return 1;
}


// ------------------------------------------------------------
// ------------------------------------------------------------

inline 
void testMaxDiff(const TString &msg, 
		 const CovarianceMatrix_t &a, 
		 const CovarianceMatrix_t &b) {
  return testMaxDiff(msg,a.GetCovM(),b.GetCovM());
}

// ------------------------------------------------------------



#endif
