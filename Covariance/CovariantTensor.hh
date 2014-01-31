#ifndef CovariantTensor_HH
#define CovariantTensor_HH

#include <TROOT.h>
#include <TMatrixD.h>
//#include <TRandom3.h>
#include <TH2D.h>
#include <TRandom.h>
#include <TString.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TPaletteAxis.h>
#include <TStyle.h>
#include <TColor.h>
#include <TText.h>

#define DYee

//for getting matrix condition number
#include <TDecompLU.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string>

#ifdef DYee
#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#else
#include "CrossSection.hh"
#endif

#include "../Include/colorPalettes.hh"
#include "../Include/BaseClass.hh"


// ------------------------------------------------------------

#ifndef ColorPalettes_HH
typedef enum { _colrange_none=0, _colrange_center=1, _colrange_positive=2, 
	       _colrange_three, _colrange_nice, 
	       _colrange_last } TColorRange_t;
#endif

// ------------------------------------------------------------

double maxOnDiagonal(const TMatrixD &M) {
  double max=0., absMax=-1.;
  for (int ir=0; ir<M.GetNrows(); ++ir) {
    for (int ic=0; ic<M.GetNcols(); ++ic) {
      const double x=M(ir,ic);
      if (fabs(x)>absMax) {
	max=x;
	absMax=fabs(x);
      }
    }
  }
  return max;
}


// ------------------------------------------------------------

double maxElement(const TMatrixD &M, int *storeRow=NULL, int *storeCol=NULL) {
  double max=0., absMax=-1.;
  int sr=-1, sc=-1;
  for (int ir=0; ir<M.GetNrows(); ++ir) {
    for (int ic=0; ic<M.GetNcols(); ++ic) {
      const double x=M(ir,ic);
      if (fabs(x)>absMax) {
	max=x;
	absMax=fabs(x);
	sr=ir; sc=ic;
      }
    }
  }
  if (storeRow) *storeRow=sr;
  if (storeCol) *storeCol=sc;
  return max;
}

// ------------------------------------------------------------

double maxElement(const TVectorD &V, int *storeRow=NULL, int *storeCol=NULL) {
  double max=0., absMax=-1.;
  int sr=-1;
  for (int ir=0; ir<V.GetNoElements(); ++ir) {
      const double x=V[ir];
      if (fabs(x)>absMax) {
	max=x;
	absMax=fabs(x);
	sr=ir;
      }
  }
  if (storeRow) *storeRow=sr;
  if (storeCol) *storeCol=-1;
  return max;
}

// ------------------------------------------------------------

double maxElement(const TH2D *h, int *storeRow=NULL, int *storeCol=NULL) {
  double max=0., absMax=-1.;
  int sr=-1, sc=-1;
  for (int ir=1; ir<=h->GetNbinsX(); ++ir) {
    for (int ic=1; ic<=h->GetNbinsY(); ++ic) {
      const double x=h->GetBinContent(ir,ic);
      if (fabs(x)>absMax) {
	max=x;
	absMax=fabs(x);
	sr=ir; sc=ic;
      }
    }
  }
  if (storeRow) *storeRow=sr;
  if (storeCol) *storeCol=sc;
  return max;
}

// ------------------------------------------------------------

template<class Container_t>
double maxDiff(const Container_t &A, const Container_t &B, int *storeRow=NULL, int *storeCol=NULL) {
  Container_t diff=A; diff-=B;
  double max= maxElement(diff,storeRow,storeCol);
  return max;
}
  
// ------------------------------------------------------------

void testMaxDiff(const TString &msg, const TVectorD &a, const TVectorD &b) {
  int i=-1;
  std::cout << msg << ". MaxDiff=" << maxDiff(a,b,&i) << "\n";
  std::cout << " values there (at " << i << "): " << a[i] << " and " << b[i] << "\n";
  return;
}

// ------------------------------------------------------------

void testMaxDiff(const TString &msg, const TMatrixD &a, const TMatrixD &b) {
  int ir=-1, ic=-1;
  std::cout << msg << ". MaxDiff=" << maxDiff(a,b,&ir,&ic) << "\n";
  std::cout << " values there, at (" << ir << "," << ic << "): " 
	    << a(ir,ic) << " and " << b(ir,ic) << "\n";
  return;
}

// ------------------------------------------------------------

double maxValue(const TVectorD &V) {
  double max=0., absMax=-1;
  for (int i=0; i<V.GetNoElements(); ++i) {
    const double x=V[i];
    if (fabs(x)>absMax) {
      max=x;
      absMax=fabs(x);
    }
  }
  return max;
}


// ------------------------------------------------------------
#ifndef ColorPalettes_HH

void set_nice_style(int nb=51) {
  const Int_t NRGBs = 5;
  const Int_t NCont = nb;
   Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
   Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
   Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
   Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
   TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
   gStyle->SetNumberContours(NCont);
}


// ------------------------------------------------------------

void set_three_style(int nb=51) {
  const Int_t NRGBs = 3;
  const Int_t NCont = nb;
   Double_t stops[NRGBs] = { 0.00, 0.50, 1.00 };
   Double_t red[NRGBs]   = { 0.10, 0.10, 1.00 };
   Double_t green[NRGBs] = { 0.10, 0.80, 0.00 };
   Double_t blue[NRGBs]  = { 1.00, 0.10, 0.00 };
   TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
   gStyle->SetNumberContours(NCont);
}


// ------------------------------------------------------------

void set_bottom_white_style(int nb=51) {
   const UInt_t Number = 3;
   Double_t Red[Number]    = { 1.00, 1.00, 1.00};
   Double_t Green[Number]  = { 1.00, 0.50, 0.00};
   Double_t Blue[Number]   = { 1.00, 0.50, 0.00};
   Double_t Length[Number] = { 0.00, 0.50, 1.00 };
   TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
   gStyle->SetNumberContours(nb);
}

// ------------------------------------------------------------

void set_center_white_style(int nb=51) {
   const UInt_t Number = 3;
   Double_t Red[Number]    = { 0.00, 1.00, 1.00};
   Double_t Green[Number]  = { 0.00, 1.00, 0.00};
   Double_t Blue[Number]   = { 0.50, 1.00, 1.00};
   Double_t Length[Number] = { 0.00, 0.50, 1.00 };
   TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
   gStyle->SetNumberContours(nb);
}

// ------------------------------------------------------------

void set_white_style_palette(TColorRange_t colRange, int nColorBins=51) {
  switch(colRange) {
  case _colrange_none: ; break;
  case _colrange_center: set_center_white_style(nColorBins); break;
  case _colrange_positive: set_bottom_white_style(nColorBins); break;
  case _colrange_three: set_three_style(nColorBins); break;
  case _colrange_nice: set_nice_style(nColorBins); break;
  default: 
    std::cout << "\n\n\tError in set_white_style_palette: not ready for int(colRange)=" << int(colRange) << "\n\n";
  }
  return;
}

// ------------------------------------------------------------

int centerHistoZAxis(TH2D *h, TColorRange_t centerRange, double maxValUser=0.) {
  if (centerRange <= _colrange_none) return 1;

  double zmin=1e9, zmax=-1e9;
  for (int i=1; i<=h->GetNbinsX(); ++i) {
    //std::cout << "i=" << i << "\n";
    for (int j=1; j<=h->GetNbinsY(); ++j) {
      double v=h->GetBinContent(i,j);
      if (v<zmin) zmin=v;
      if (v>zmax) zmax=v;
    }
  }

  if (centerRange > _colrange_none) {
    std::cout << h->GetName() << ": zmin=" << zmin << ", zmax=" << zmax << "\n";
    if ((centerRange==_colrange_center) ||
	(centerRange==_colrange_three)) {
      if (zmin*zmax >= 0) {
	std::cout << " warning in centerHistoZAxis: cannot center (zmin=" 
		  << zmin << ", zmax=" << zmax << ")\n";
      }

      {
	double z=zmax;
	if (-zmin > z) z= -zmin;
	std::cout << "histoName=<" << h->GetName() << ">, max|z|= " << z << "\n";
	if (maxValUser>0) {
	  z=maxValUser;
	  std::cout << " maxValUser=" << maxValUser << "\n";
	}
	h->GetZaxis()->SetRangeUser(-z,z);
      }
    }
    else if (centerRange==_colrange_positive) {
      if (maxValUser>0.) {
	if (maxValUser < zmax) {
	  std::cout << "histoName=<" << h->GetName() << ">, warning: zmax=" << zmax << ", maxValUser=" << maxValUser << "\n"; 
	}
	h->GetZaxis()->SetRangeUser(0,maxValUser);
      }
    }
  }

  return 1;
}

#endif

// ------------------------------------------------------------

TH2D* createHisto2D(const TMatrixD &M, const TMatrixD *Merr, const char *histoName, const char *histoTitle, TColorRange_t centerRange, int massBins=0, double maxValUser=0.) {
  int nRows=M.GetNrows();
  int nCols=M.GetNcols();
  TH2D *h=NULL;

  if (massBins) {
#ifndef DYTools_HH
    std::cout << "createHisto2D. Warning: massBins=1, but DYTools is not available\n";
#endif
    if ((nRows!=DYTools::nMassBins) || (nCols!=DYTools::nMassBins)) {
      std::cout << "createHisto2D. Warning: massBins=1, DYTools is available, but the matrix has incorrect number of bins (" << nRows << "x" << nCols << "), instead of a square matrix of dim=" << DYTools::nMassBins << "\n";
    }
    else {
      double *mbins=new double[DYTools::nMassBins+1];
      for (int i=0; i<=DYTools::nMassBins; ++i) 
	mbins[i]=DYTools::massBinLimits[i];
      mbins[DYTools::nMassBins] += 5e-1;
      h= new TH2D(histoName,histoTitle,
		  DYTools::nMassBins, mbins,
		  DYTools::nMassBins, mbins);
      delete mbins;
    }
    //h->GetXaxis()->SetTitle("mass_{ee;1}");
    //h->GetYaxis()->SetTitle("mass_{ee;2}");
  }

  if (!h) {
    std::cout << "createHisto2D: nRows=" << nRows << ", nCols=" << nCols << "\n";
    h= new TH2D(histoName,histoTitle,
		nRows,0.5,nRows+0.5, 
		nCols,0.5,nCols+0.5);
    //h->GetXaxis()->SetTitle("idx_{1}");
    //h->GetYaxis()->SetTitle("idx_{2}");
  }
  h->SetDirectory(0);
  for (int ir=0; ir<nRows; ir++) {
    for (int ic=0; ic<nCols; ic++) {
      h->SetBinContent(ir+1,ic+1, M(ir,ic));
      //std::cout << "set (" << (ir+1) << "," << (ic+1) << ")=" << M(ir,ic) << "\n";
      if (Merr) h->SetBinError(ir+1,ic+1, (*Merr)(ir,ic));
      else h->SetBinError(ir+1,ic+1, 0.);
    }
  }
  
  if (centerRange > _colrange_none) {
    if (!centerHistoZAxis(h,centerRange,maxValUser)) {
      std::cout << "error in createHisto2D\n";
    }
  }

  return h;
}

 
// ------------------------------------------------------------
#ifndef ColorPalettes_HH

int drawHistoSubpad(TCanvas *c, int subPad, TH2D *h2, TColorRange_t centerRange, int useMassBins, int nColorBins=51) {
  c->cd(subPad);
  TPad *pad=(TPad*)c->GetPad(subPad);

  if (useMassBins) {
    pad->SetLogy(1);
    pad->SetLogx(1);

    h2->GetXaxis()->SetMoreLogLabels();
    h2->GetXaxis()->SetNoExponent();
    h2->GetYaxis()->SetMoreLogLabels();
    h2->GetYaxis()->SetNoExponent();
    int ndiv=703;
    h2->GetXaxis()->SetNdivisions(ndiv,true);
    h2->GetYaxis()->SetNdivisions(ndiv,true);
  }
  switch (centerRange) {
  case _colrange_center: set_center_white_style(nColorBins); break;
  case _colrange_positive: set_bottom_white_style(nColorBins); break;
  default: ; 
  }
  h2->Draw("colz");
  //if (useMassBins) {
  //  h2->GetXaxis()->SetRangeUser(15,4000.);
  //  h2->GetYaxis()->SetRangeUser(15,4000.);
  //}
  c->Modified();
  return 1;
}

// ------------------------------------------------------------

int drawHisto(TCanvas *c, TH2D *h2, TColorRange_t centerRange, int useMassBins, int nColorBins=51) {
  int res=drawHistoSubpad(c,0,h2,centerRange,useMassBins,nColorBins);
  c->Update();
  return res;
}
 
#endif
// ------------------------------------------------------------

template<class Histo_t>
void MovePalette (TCanvas *c, Histo_t *h, 
		  double x1NDC, double x2NDC, double y1NDC, double y2NDC) {
  TPaletteAxis *palette = 
    (TPaletteAxis*)h->GetListOfFunctions()->FindObject("palette");
  palette->SetX1NDC(x1NDC);
  palette->SetX2NDC(x2NDC);
  palette->SetY1NDC(y1NDC);
  palette->SetY2NDC(y2NDC);
  c->Modified();
}


// ------------------------------------------------------------

void Transpose(TMatrixD &M) {
  TMatrixD tmp(M);
  M.Transpose(tmp);
}

// ------------------------------------------------------------
// ------------------------------------------------------------

  // (nRows,nCols) are the dimensions of a matrix. 
  // The CovarianceTensor contains covariances of each element (i,j) with (k,l)
  // thus the size is (nR x nC) x (nR x nC).

class CovarianceTensor_t : public BaseClass_t {
  int FDim,FNrows,FNcols;
  TMatrixD CovT;
  TMatrixD TinvAverage, TinvError;
  TMatrixD CorrT;
  TString FName;
public:
  // if controlDim=0, the tensor is empty
  CovarianceTensor_t(const TString &name, int controlDim, int nDim) : 
    BaseClass_t("CovarianceTensor_t"), 
    FDim(controlDim*nDim), FNrows(nDim), FNcols(nDim),
    CovT(FDim,FDim), 
    TinvAverage(nDim,nDim), TinvError(nDim,nDim),
    CorrT(FDim,FDim),
    FName(name)
    {
      std::cout << "CovariantTensor_t(" << name << ", controlDim=" 
		<< controlDim << ", nDim=" << nDim << "); FDim=" 
		<< FDim << "\n";
    }

  void SetName(const TString &name) { FName=name; }
  int GetDim() const { return FDim; }
  int GetNrows() const { return FNrows; }
  int GetNcols() const { return FNcols; }
  const TString& GetName() const { return FName; }
  const TMatrixD& GetCovTensor() const { return CovT; }
  //const TMatrixD* GetCovTensorPtr() const { return CovT; }
  const TMatrixD& GetCorrTensor() const { return CorrT; }
  //const TMatrixD* GetCorrTensorPtr() const { return CorrT; }
  const TMatrixD& GetInvAverage() const { return TinvAverage; }
  const TMatrixD& GetInvAverageErr() const { return TinvError; }

  // hack access
  TMatrixD &EditInvAverage() { return TinvAverage; }
  TMatrixD &EditAverageErr() { return TinvError; }

  void Zero() { 
    FDim=0; FNcols=0; FNrows=0;
    CovT.Zero(); 
    TinvAverage.Zero();
    TinvError.Zero();
    CorrT.Zero();
  }

  int FlatIdx(int row, int col) const { return row*FNcols + col; }
  double operator()(int idx1, int idx2) { return CovT(idx1,idx2); }

  // edits

  int SetDiagonal(const TMatrixD &M, const TMatrixD &Merr) {
    TinvAverage=M;
    TinvError=Merr;
    if (FDim) {
      for (int ir=0, idx=0; ir<FNrows; ir++) {
	TMatrixDRow_const rowErr(Merr,ir);
	for (int ic=0; ic<FNcols; ++ic, ++idx) {
	  CovT(idx,idx)= rowErr[ic]*rowErr[ic];
	}
      }
    }
    return 1;
  }

  
  // -------- Main method
  int calculateInvertedMatrixCovariance(
      const TMatrixD &T,
      const TMatrixD &TErrPos, const TMatrixD &TErrNeg,
      int nExps=10000, int debug=0);

  // --------------- Tests

  int checkSymmetry() const {
    int ok=1;
    std::cout << "CovariantTensor::checkSymmetry -- \n";
    if (FDim) {
      for (int ir=0; ok && (ir<FNrows); ++ir) {
	TMatrixDRow_const row(CovT,ir);
	TMatrixDColumn_const col(CovT,ir);
	for (int ic=0; ok && (ic<FNcols); ++ic) {
	  if (row[ic] != col[ic]) {
	    ok=0;
	    std::cout << "  failing at (ir,ic)=(" << ir << "," << ic << "): "
		      << row[ic] << " vs " << col[ic] << "\n";
	  }
	}
      }
    }
    if (ok) std::cout << "  test passed\n";
    return ok;
  }

  // ---------------

  double checkDiagonal() const {
    double max=0., maxAbs=0.;
    if (FDim) {
      for (int ir=0, idx=0; ir<FNrows; ++ir) {
	for (int ic=0; ic<FNcols; ++ic, ++idx) {
	  double diff= CovT(idx,idx) - TinvError(ir,ic);
	  if (fabs(diff)>maxAbs) {
	    max=diff;
	    maxAbs=fabs(diff);
	  }
	}
      }
    }
    return max;
  }

  // ---------------

  double checkDiagonalSqrt() const {
    double max=0., maxAbs=0.;
    if (FDim) {
      for (int ir=0, idx=0; ir<FNrows; ++ir) {
	for (int ic=0; ic<FNcols; ++ic, ++idx) {
	  double diff= sqrt( CovT(idx,idx) )  - TinvError(ir,ic);
	  if (fabs(diff)>maxAbs) {
	    max=diff;
	    maxAbs=fabs(diff);
	  }
	}
      }
    }
    return max;
  }

  // ---------------

  double checkDiagonal(const TMatrixD &referenceMatrix) const {
    double max=0., maxAbs=0.;
    if (FDim) {
      for (int ir=0, idx=0; ir<FNrows; ++ir) {
	for (int ic=0; ic<FNcols; ++ic, ++idx) {
	  double diff= CovT(idx,idx) - referenceMatrix(ir,ic);
	  if (fabs(diff)>maxAbs) {
	    max=diff;
	    maxAbs=fabs(diff);
	  }
	}
      }
    }
    return max;
  }

  // ---------------

  int SaveObj(TFile &fout) {
    fout.cd();
    CovT.Write("covarianceTensor");
    CorrT.Write("correlationTensor");
    TinvAverage.Write("invRespAverage");
    TinvError.Write("invRespError");
    return 1;
  }

  // ---------------

  int Save(TFile &fout, const char *origMName=NULL, TMatrixD *origM=NULL, TMatrixD *origMErr=NULL) {
    if (!fout.IsOpen()) return 0;
    if (!SaveObj(fout)) return 0;
    if (origMName && origM) { 
      origM->Write(origMName);
      if (origMErr) origMErr->Write(Form("%sErr",origMName));
    }
    return 1;
  }

  // ---------------

  int Save(const TString &fname, const char *origMName=NULL, TMatrixD *origM=NULL, TMatrixD *origMErr=NULL) {
    TFile fout(fname,"recreate");
    if (!fout.IsOpen()) return 0;
    if (!Save(fout,origMName,origM,origMErr)) return 0;
    fout.Close();
    return 1;
  }

  // ---------------

  int Load(TFile &fin) {
    fin.cd();
    TMatrixD* covTptr=(TMatrixD*)fin.Get("covarianceTensor");
    TMatrixD* corrTptr=(TMatrixD*)fin.Get("correlationTensor");
    TMatrixD* tinvAveragePtr=(TMatrixD*)fin.Get("invRespAverage");
    TMatrixD* tinvErrorPtr=(TMatrixD*)fin.Get("invRespError");
    if (!covTptr || !corrTptr || !tinvAveragePtr || !tinvErrorPtr) return 0;
    CovT= *covTptr;
    CorrT= *corrTptr;
    TinvAverage= *tinvAveragePtr;
    TinvError= *tinvErrorPtr;
    delete covTptr;
    delete corrTptr;
    delete tinvAveragePtr;
    delete tinvErrorPtr;
    FNrows= TinvAverage.GetNrows();
    FNcols= TinvAverage.GetNcols();
    FDim= CovT.GetNrows();
    return 1;
  }

  // ---------------

  int Load(const TString &fname, const char *origMName=NULL, TMatrixD *origM=NULL, TMatrixD *origMErr=NULL) {
    TFile fin(fname,"read");
    if (!fin.IsOpen()) return 0;
    if (!Load(fin)) return 0;
    if (origMName && origM) { 
      TMatrixD *tmpM= (TMatrixD*)fin.Get(origMName);
      TMatrixD *tmpMErr= (TMatrixD*)fin.Get(Form("%sErr",origMName));
      if (!tmpM || (origMErr && !tmpMErr)) return 0;
      *origM = *tmpM;
      if (origMErr) *origMErr = *tmpMErr;
      if (tmpM) delete tmpM;
      if (tmpMErr) delete tmpMErr;
    }
    fin.Close();
    return 1;
  }

  // ---------------

  TH2D* Draw(TCanvas *c, TColorRange_t centerRange=_colrange_center, int massBins=0, int nColorBins=51, double maxValUser=0.) const {
    TString name=TString("h2D_") + FName + TString("_covT");
    TString title=name;
    TH2D* histo=createHisto2D(CovT,NULL, name.Data(),title.Data(), 
			      centerRange, (FDim) ? massBins : 0,
			      maxValUser);
    if (massBins==0) {
      int ndiv=504;
      bool optimize=true;
      histo->GetXaxis()->SetRangeUser(0,CovT.GetNrows());
      histo->GetYaxis()->SetRangeUser(0,CovT.GetNcols());
      histo->GetXaxis()->SetNdivisions(ndiv,optimize);
      histo->GetYaxis()->SetNdivisions(ndiv,optimize);
    }

    drawHisto(c,histo,centerRange,massBins,nColorBins);
    return histo;
  }

  // ---------------
  /*
  TH2D* DrawLines(TCanvas *c, TColorRange_t centerRange=_colrange_center, int massBins=0, int nColorBins=15, double maxValUser=0.) const {
    TString name=TString("h2D_") + FName + TString("_covT_lines");
    TString title=name;
    TH2D* histo=createHisto2D(CovT,NULL, name.Data(),title.Data(), 
			      centerRange, massBins,maxValUser);
    if (massBins==0) {
      int ndiv=504;
      bool optimize=true;
      histo->GetXaxis()->SetRangeUser(0,CovT.GetNrows());
      histo->GetYaxis()->SetRangeUser(0,CovT.GetNcols());
      histo->GetXaxis()->SetNdivisions(ndiv,optimize);
      histo->GetYaxis()->SetNdivisions(ndiv,optimize);
    }

    drawHisto(c,histo,centerRange,massBins,nColorBins);
    for (int ir=0; ir<CovT.GetNrows(); ++ir) {
      TMatrixDRow_const row(CovT,ir);
      double maxVal=0;
      for (int ic=0; ic<CovT.GetNcols(); ++ic) {
	
      }
    }
    return histo;
  }
  */

  // ---------------
  // ---------------

  void Print() const {
    std::cout << "CovarianceTensor(" << FName << ")\n";
    std::cout << "  -- TinvAverage:\n";
    TinvAverage.Print();
    std::cout << "  -- TinvError:\n";
    TinvError.Print();
    std::cout << "  -- CovT (covariance tensor):\n";
    CovT.Print();
    std::cout << "  -- CorrT (correlation tensor):\n";
    CorrT.Print();
    std::cout << std::endl;
  }
};

// ------------------------------------------------------------
// ------------------------------------------------------------

// ------------------------------------------------------------

int CovarianceTensor_t::calculateInvertedMatrixCovariance(
   const TMatrixD &T, 
   const TMatrixD &TErrPos, const TMatrixD &TErrNeg,
   int nExps, int debug) {

  // Calculate covariances of the inverted matrix by the Monte Carlo
  // method

  Double_t det = 0;
  FNrows = T.GetNrows();
  FNcols = T.GetNcols();

  TMatrixD TinvSum(FNrows,FNcols);
  TMatrixD TinvSumSquares(FNrows,FNcols);
  CovarianceTensor_t crossTerms("crossTerms",FNrows,FNcols);

  // Reset Matrix where we will be accumulating RMS/sigma:
  TinvSum        .Zero();
  TinvSumSquares .Zero();
  //assert(this->SetSize(FNrows,FNcols)); // also sets to 0

  double detSqSum=0, detAbsSum=0;

  // Do many tries, accumulate RMS
  for(int iTry = 0; iTry<nExps; iTry++){
    HERE("iTry=%d/%d",iTry,nExps-1);
    // Find the smeared matrix
    TMatrixD Tsmeared = T;
    for(int i = 0; i<FNrows; i++){
      for(int j = 0; j<FNcols; j++){
	double central = T(i,j);
	double sigPos = TErrPos(i,j);
	double sigNeg = TErrNeg(i,j);
 	// Switch to symmetric errors: approximation, but much simpler
	double sig = (sigPos+sigNeg)/2.0;
	Tsmeared(i,j) = gRandom->Gaus(central,sig);
      }
    }

    if ((debug==10) && (iTry==0)) {
      if (1) {
	Tsmeared.Draw("colz");
	TDecompLU lu(Tsmeared);
	std::cout << "lu.Condition=" << lu.Condition() << ", condNumber=" << (-lu.Condition() * Tsmeared.Norm1()) << "\n";
      }
      else if (0) {
	TMatrixD diff=T;
	diff-=Tsmeared;
	diff.Draw("colz");
      }
    }

    //std::cout << "Tsmeared:\n"; Tsmeared.Print();

    // Find the inverted to smeared matrix
    TMatrixD TinvSmeared = Tsmeared;
    TinvSmeared.Invert(&det);
    HERE("det=%6.4e",det);
    detAbsSum+=fabs(det);
    detSqSum+= det*det;

    // Accumulate sum and sum of squares for each element
    for(int i = 0; i<FNrows; i++){
      for(int j = 0; j<FNcols; j++){
	TinvSum       (i,j) += TinvSmeared(i,j);
	TinvSumSquares(i,j) += TinvSmeared(i,j)*TinvSmeared(i,j);
      }
    }

    // Accumulate the terms R[r1,c1]*R[r2,c2]
    TMatrixD tmpM(FNrows,FNcols);

    for (int r1 = 0; r1<FNrows; r1++) {
      //HERE("r1=%d",r1);
      TMatrixDRow TinvRow(TinvSmeared,r1);
      for (int c1 = 0; c1<FNcols; c1++) {
	// prepare TinvRow[r1,c1] * TinvRow;
	tmpM = TinvSmeared;
	tmpM *= TinvRow[c1];
	// prepare data access -- flat indexing
	int idx1=this->FlatIdx(r1,c1);
	TMatrixDRow ctRow( crossTerms.CovT, idx1);
	// add values to this row
	for (int r2 = 0, idx2=0; r2<FNrows; r2++) {
	  TMatrixDRow_const tmpMrow(tmpM,r2);
	  for (int c2 = 0; c2<FNcols; c2++, idx2++) {
	    ctRow(idx2) += tmpMrow[c2];
	  }
	}
      }
    }
    //HERE("accumulation done");
  }

  std::cout << "checkDiagonal(TinvSumSquares) = " << crossTerms.checkDiagonal(TinvSumSquares) << "\n";

  // Calculate the error matrix
  for(int i = 0; i<FNrows; i++){
    for(int j = 0; j<FNcols; j++){
      (this->TinvAverage)(i,j) = TinvSum(i,j)/double(nExps);
      (this->TinvError)(i,j) = 
	sqrt( 
	       TinvSumSquares(i,j)/double(nExps) 
		- (TinvSum(i,j)/double(nExps))*(TinvSum(i,j)/double(nExps)) 
	       );
    }
  }

  // Calculate the covariance matrix

  crossTerms.CovT *= 1/double(nExps);
  //std::cout << "checkDiagonal(TinAverage) = " << crossTerms.checkDiagonal(this->GetInvAverage()) << "\n";

  if (0) {
    for (int r1=0; r1<FNrows; ++r1) {
      for (int c1=0; c1<FNcols; ++c1) {
	const int idx1=this->FlatIdx(r1,c1);
	TMatrixDRow covRow( this->CovT, idx1 );
	TMatrixDRow_const ctRow( crossTerms.CovT, idx1 );
	const double aAvg= (this->TinvAverage)(r1,c1);
	for (int r2=0; r2<FNrows; ++r2) {
	  for (int c2=0; c2<FNcols; ++c2) {
	    const int idx2=this->FlatIdx(r2,c2);
	    const double bAvg= (this->TinvAverage)(r2,c2);
	    covRow(idx2) =  ctRow(idx2) - aAvg * bAvg;
	  }
	}
      }
    }
  }
  else {

    for (int r1=0, idx1=0; r1<FNrows; ++r1) {
      for (int c1=0; c1<FNcols; ++c1, ++idx1) {
	TMatrixDRow covRow( this->CovT, idx1 );
	TMatrixDRow_const ctRow( crossTerms.CovT, idx1 );
	const double aAvg= (this->TinvAverage)(r1,c1);
	for (int r2=0, idx2=0; r2<FNrows; ++r2) {
	  for (int c2=0; c2<FNcols; ++c2, ++idx2) {
	    const double bAvg= (this->TinvAverage)(r2,c2);
	    covRow(idx2) =  ctRow(idx2) - aAvg * bAvg;
	  }
	}
      }
    }

  }

  // fill correlation matrix
  TMatrixDDiag_const covDiag( this->CovT );
  for (int r1=0, idx1=0; r1<FNrows; ++r1) {
    for (int c1=0; c1<FNcols; ++c1, ++idx1) {
      TMatrixDRow corrRow( this->CorrT, idx1 );
      TMatrixDRow_const covRow( this->CovT, idx1 );
      for (int r2=0, idx2=0; r2<FNrows; ++r2) {
	for (int c2=0; c2<FNcols; ++c2, ++idx2) {
	  corrRow(idx2) = covRow(idx2)/sqrt(covDiag(idx1)*covDiag(idx2));
	}
      }
    }
  }

  return 1;
}

// ------------------------------------------------------------
// ------------------------------------------------------------



#endif
