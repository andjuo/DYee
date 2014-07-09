#include "../Include/MyTools.hh"
#include <fstream>
#include <TBenchmark.h>

//--------------------------------------------------
//--------------------------------------------------

void printMeanAndRMS(TH1D *histo, const char *extra_msg) {
  if (extra_msg) std::cout << extra_msg;
  if (!histo) { std::cout << "histo is null\n"; return; }
  std::cout << histo->GetName()
	    << " mean=" << histo->GetMean()
	    << " +- " << histo->GetRMS()
	    << "\n";
}

//--------------------------------------------------
//--------------------------------------------------

int PrintHisto2Dvec(const char *msg, const std::vector<TH2D*> &vec,
		     int exponent, int maxLines) {
  const char *format_range= " %5.2f-%5.2f  %5.2f-%5.2f ";
  const char *format_char= (exponent) ?  "    %e    %e" : "    %f    %f";

  if (msg) std::cout << msg;
  std::cout << "values of histos:";
  for (unsigned int i=0; i<vec.size(); ++i) {
    std::cout << vec[i]->GetName() << " ";
  }
  std::cout << "\n";

  TH2D* histo=vec[0];
  int truncated=0;
  int imax= histo->GetNbinsX() * histo->GetNbinsY();
  if ((maxLines>0) && (imax>maxLines)) { imax=maxLines; truncated=1; }

  int shown=0;
  for(int i=1; i<=histo->GetNbinsX(); i++) {
    double x=histo->GetBinLowEdge(i);
    double w=histo->GetBinWidth(i);
    for (int j=1; j<=histo->GetNbinsY(); ++j, shown++) {
      double y=histo->GetYaxis()->GetBinLowEdge(j);
      double h=histo->GetYaxis()->GetBinWidth(j);
      std::cout << Form(format_range,x,x+w,y,y+h);
      for (unsigned int iv=0; iv<vec.size(); ++iv) {
	std::cout << Form(format_char,
			  vec[iv]->GetBinContent(i,j),
			  vec[iv]->GetBinError(i,j));
      }
      if (shown>=imax) break;
    }
    std::cout << "\n";
    if (shown>=imax) break;
  }
  if (truncated) std::cout << "... output truncated to " << maxLines << " lines\n";
  return 1;
}

//--------------------------------------------------

int PrintTwoHistos(const char *msg, TH2D *h2a, TH2D *h2b,
		   int exponent, int maxLines) {
  std::vector<TH2D*> hV;
  hV.push_back(h2a);
  hV.push_back(h2b);
  return PrintHisto2Dvec(msg,hV,exponent,maxLines);
}

//--------------------------------------------------

int printHistoAssumeMassBinning(const TH1D* histo,
				int exponent, int maxLines) {
  TH1D *h1= convertBinNo2Mass(histo);
  if (!h1 ||
      !printHisto(h1,exponent,maxLines)) {
    std::cout << "printHistoAssumeMassBinning (orig histoName="
	      << histo->GetName() << "): failed to prepare or print\n";
    return 0;
  }
  delete h1;
  return 1;
}


//--------------------------------------------------
//--------------------------------------------------

std::vector<TString>* createMassRangeVec(TString prependStr) {
  std::vector<TString>* vec=new std::vector<TString>();
  vec->reserve(DYTools::nMassBins);
  for (int im=0; im<DYTools::nMassBins; ++im) {
    TString mStr=Form("M_%1.0lf_%1.0lf",
		  DYTools::massBinLimits[im],DYTools::massBinLimits[im+1]);
    if (prependStr.Length()) mStr.Prepend(prependStr);
    vec->push_back(mStr);
  }
  return vec;
}

//--------------------------------------------------

int replaceAll(std::vector<TString*> &vec, TString oldText, TString newText) {
  int count=0;
  for (unsigned int i=0; i<vec.size(); ++i) {
    TString *s=vec[i];
    if (s->Index(oldText)!=-1) {
      count++;
      s->ReplaceAll(oldText,newText);
    }
  }
  return count;
}

//--------------------------------------------------

int replaceAll(std::vector<TString> &vec, TString oldText, TString newText) {
  int count=0;
  for (unsigned int i=0; i<vec.size(); ++i) {
    TString s=vec[i];
    if (s.Index(oldText)!=-1) {
      count++;
      vec[i].ReplaceAll(oldText,newText);
    }
  }
  return count;
}

//--------------------------------------------------

TString niceNumber(int iVal, int iValMax) {
  if (iValMax<0) iValMax=-iValMax;
  int cnt=0;
  while (iValMax>0) {
    cnt++;
    iValMax/=10;
  }
  const char *format=Form("%c%dd",'%',cnt);
  TString str=Form(format,iVal);
  str.ReplaceAll(" ","0");
  return str;
}

//--------------------------------------------------
//--------------------------------------------------

void printHisto(const std::vector<TH2D*> hV, int exponent, int maxLines, int maxEntries) {
  unsigned int imax=hV.size();
  if ((maxEntries>0) && (imax>(unsigned int)(maxEntries))) {
    imax=(unsigned int)(maxEntries);
  }
  std::cout << "printing " << imax << " entries of a vector hV[" << hV.size() << "]\n";
  for (unsigned int i=0; i<imax; ++i) {
    printHisto(hV[i],exponent,maxLines);
  }
  return ;
}

//---------------------------------------------------------------

TH2D* removeUnderflow(TH2D* h2orig, TString newName) {
  if ((DYTools::massBinningSet!=DYTools::_MassBins_2011) &&
      (DYTools::massBinningSet!=DYTools::_MassBins_2012)) {
    std::cout << dashline;
    std::cout << "removeUnderflow is not ready for "
	      << DYTools::massBinningSet << "\n";
    return NULL;
  }
  if (DYTools::study2D==0) {
    TH2D* h2=Clone(h2orig,newName);
    return h2;
  }

  double *newMassBinLimits= new double[DYTools::nMassBins];
  for (int i=0; i<DYTools::nMassBins; ++i) {
    newMassBinLimits[i]= DYTools::massBinLimits[i+1];
  }

  TH2D* h2=new TH2D(newName,newName,
		    DYTools::nMassBins-1,newMassBinLimits,
		    DYTools::nYBinsMax,DYTools::yRangeMin,DYTools::yRangeMax);
  delete newMassBinLimits;

  h2->SetDirectory(0);
  h2->SetStats(0);
  h2->Sumw2();
  h2->GetXaxis()->SetTitle("M_{ee} [GeV]");
  h2->GetYaxis()->SetTitle("|y|");
  for (int ibin=1; ibin<=h2->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=h2->GetNbinsY(); ++jbin) {
      h2->SetBinContent( ibin,jbin, h2orig->GetBinContent(ibin+1,jbin) );
      h2->SetBinError  ( ibin,jbin, h2orig->GetBinError  (ibin+1,jbin) );
    }
  }
  return h2;
}


//---------------------------------------------------------------
//---------------------------------------------------------------

TMatrixD* deriveCovMFromRndStudies(const std::vector<TH2D*> &rndV,
				   int unbiasedEstimate,
				   TH2D *avgDistr,
				   TMatrixD *covError) {
  //const int unbiasedEstimate=1;
  if (rndV.size()<(unsigned int)(1+unbiasedEstimate)) {
    HERE("deriveCovMFromRndStudies(vec(TH2D*)) vec is empty or has 1 entry");
    return NULL;
  }
  int dim=DYTools::nUnfoldingBins;
  //dim=2; std::cout << "deriveCovMFromRndStudies : dimension Hack\n";
  // container for sum(x)
  TVectorD VSumX(dim);
  VSumX.Zero();
  // container for sum(x y)
  TMatrixD *MSumXY=new TMatrixD(dim,dim);
  if (!MSumXY) {
    HERE("deriveCovMFromRndStudies: failed to create the container");
    return NULL;
  }
  MSumXY->Zero();

  // 1. Accumulate sums
  // iterate over 1st histo
  int rMax=rndV[0]->GetNbinsX();
  int cMax=rndV[0]->GetNbinsY();
  for (unsigned int i=0; i<rndV.size(); i++) {
    const TH2D *hi=rndV[i];
    for (int ibr1=1; ibr1<=rMax; ++ibr1) {
      for (int ibc1=1; ibc1<=cMax; ++ibc1) {
	int idxFlat1=DYTools::findIndexFlat(ibr1-1,ibc1-1);
	if (idxFlat1<0) continue; // outside of considered space
	double val1=hi->GetBinContent(ibr1,ibc1);
	VSumX(idxFlat1) += val1;
      }
    }
  }

  for (unsigned int i=0; i<rndV.size(); i++) {
    const TH2D *hi=rndV[i];

    for (int ibr1=1; ibr1<=rMax; ++ibr1) {
      for (int ibc1=1; ibc1<=cMax; ++ibc1) {
	int idxFlat1=DYTools::findIndexFlat(ibr1-1,ibc1-1);
	if (idxFlat1<0) continue; // outside of considered space
	double val1=hi->GetBinContent(ibr1,ibc1);

	for (int ibr2=ibr1; ibr2<=rMax; ++ibr2) {
	  for (int ibc2=(ibr1==ibr2) ? ibc1 : 1; ibc2<=cMax; ++ibc2) {
	    int idxFlat2=DYTools::findIndexFlat(ibr2-1,ibc2-1);
	    if (idxFlat2<0) continue; // outside of considered space
	    double val2=hi->GetBinContent(ibr2,ibc2);
	    (*MSumXY)(idxFlat1,idxFlat2) += val1*val2;
	    if (idxFlat1 != idxFlat2) {
	      // due to the optimization ibr2=ibr1, we have to
	      // fill the mirrored term
	      (*MSumXY)(idxFlat2,idxFlat1) += val1*val2;
	    }
	  }
	}
	/*
	for (int ibr2=ibr1; ibr2<=rMax; ++ibr2) {
	  for (int ibc2=ibc1; ibc2<=cMax; ++ibc2) {
	    int idxFlat2=DYTools::findIndexFlat(ibr2-1,ibc2-1);
	    if (idxFlat2<0) continue; // outside of considered space
	    double val2=hi->GetBinContent(ibr2,ibc2);
	    (*MSumXY)(idxFlat1,idxFlat2) += val1*val2;
	    if (idxFlat1 != idxFlat2) {
	      // due to the optimization ibr2=ibr1, we have to
	      // fill the mirrored term
	      (*MSumXY)(idxFlat2,idxFlat1) += val1*val2;
	    }
	  }
	}
	*/
      }
    }
  }

  // 2. Derive the averages and the covariances
  double nExps=double(rndV.size());
  double factor=(unbiasedEstimate) ? (nExps-1) : nExps;
  double correctionFactor = nExps/factor;

  VSumX    *= (1/nExps);

  //std::cout << "nExps=" << nExps << ", VSumX="; VSumX.Print();

  (*MSumXY)*= (1/factor);
  for (int idxFlat1=0; idxFlat1<dim; ++idxFlat1) {
    for (int idxFlat2=0; idxFlat2<dim; ++idxFlat2) {
      (*MSumXY)(idxFlat1,idxFlat2) -= correctionFactor*(VSumX(idxFlat1) * VSumX(idxFlat2));
    }
  }

  // 3. If needed fill the average container
  if (avgDistr) {
    avgDistr->Reset();
    for (int ibr=1; ibr<=avgDistr->GetNbinsX(); ++ibr) {
      for (int ibc=1; ibc<=avgDistr->GetNbinsY(); ++ibc) {
	int idxFlat= DYTools::findIndexFlat(ibr-1,ibc-1);
	if (idxFlat<0) continue; // outside of considered space
	double errSqr=(*MSumXY)(idxFlat,idxFlat);
	if (errSqr<0) {
	  std::cout << "deriveCovMFromRndStudies: negative errSqr at "
		    << "idx=" << idxFlat
		    << ", centralValue=" << VSumX(idxFlat)
		    << ", err^2=" << errSqr << "\n";
	  const double tol=1e-6;
	  if (-errSqr/VSumX(idxFlat) < tol) {
	    std::cout << " - the relative error is below tolerance "
		      << Form("(%5.2e)",tol) << ", reverting the sign\n";
	    errSqr=-errSqr;
	    (*MSumXY)(idxFlat,idxFlat)= errSqr;
	  }
	  else {
	    std::cout << errdash;
	  }
	}
	avgDistr->SetBinContent(ibr,ibc, VSumX(idxFlat));
	avgDistr->SetBinError  (ibr,ibc, sqrt(errSqr));
      }
    }
  }

  if (covError) {
    covError->Zero();
    for (unsigned int i=0; i<rndV.size(); i++) {
      const TH2D *hi=rndV[i];
      for (int ibr1=1; ibr1<=rMax; ++ibr1) {
	for (int ibc1=1; ibc1<=cMax; ++ibc1) {
	  int idxFlat1=DYTools::findIndexFlat(ibr1-1,ibc1-1);
	  if (idxFlat1<0) continue; // outside of considered space
	  double val1=hi->GetBinContent(ibr1,ibc1);

	  for (int ibr2=ibr1; ibr2<=rMax; ++ibr2) {
	    for (int ibc2=(ibr1==ibr2) ? ibc1 : 1; ibc2<=cMax; ++ibc2) {
	      int idxFlat2=DYTools::findIndexFlat(ibr2-1,ibc2-1);
	      if (idxFlat2<0) continue; // outside of considered space
	      double val2=hi->GetBinContent(ibr2,ibc2);
	      double covEstimate= (val1-VSumX(idxFlat1))*(val2-VSumX(idxFlat2));
	      double dev= covEstimate - (*MSumXY)(idxFlat1,idxFlat2);
	      (*covError)(idxFlat1,idxFlat2) += dev*dev;
	      if (idxFlat1 != idxFlat2) {
		(*covError)(idxFlat2,idxFlat1) += dev*dev;
	      }
	    }
	  }
	}
      }
    }
    (*covError) *= (1/nExps);
    //unsquareElements(*covError);
  }

  return MSumXY;
}

//--------------------------------------------------
//--------------------------------------------------

TMatrixD* corrFromCov(const TMatrixD &cov) {
  TMatrixD *corr= new TMatrixD(cov);
  if (!corr) {
    std::cout << "corrFromCov: failed to create a new matrix\n";
    return NULL;
  }
  corr->Zero();
  for (int ir=0; ir<cov.GetNrows(); ++ir) {
    double eR=sqrt(cov(ir,ir));
    for (int ic=0; ic<cov.GetNcols(); ++ic) {
      double eC=sqrt(cov(ic,ic));
      (*corr)(ir,ic) = cov(ir,ic)/(eR*eC);
    }
  }
  return corr;
}

//--------------------------------------------------

TH2D* errorFromCov(const TMatrixD &cov, TString newName) {
  TH2D* h2=createBaseH2(newName,newName,1);
  int iflat=0;
  for (int xbin=1; xbin<=h2->GetNbinsX(); ++xbin) {
    for (int ybin=1; ybin<=h2->GetNbinsY(); ++ybin, ++iflat) {
      if (iflat >= cov.GetNrows()) break;
      h2->SetBinContent(xbin,ybin, sqrt(cov(iflat,iflat)));
      h2->SetBinError(xbin,ybin, 0.);
    }
  }
  return h2;
}

//--------------------------------------------------

TMatrixD* covFromCorr(const TMatrixD &corr, const TVectorD &errs) {
  if ((corr.GetNrows() != errs.GetNoElements()) ||
      (corr.GetNcols() != errs.GetNoElements())) {
    std::cout << "covFromCorr dimension mismatch:\n";
    std::cout << "  corr[" << corr.GetNrows() << " x "
	      << corr.GetNcols() << "]\n";
    std::cout << "  errs[" << errs.GetNoElements() << "]\n";
    return NULL;
  }
  TMatrixD *cov= new TMatrixD(corr);
  if (!cov) {
    std::cout << "covFromCorr: failed to create a new matrix\n";
    return NULL;
  }
  cov->Zero();
  for (int ir=0; ir<corr.GetNrows(); ++ir) {
    double eR=errs[ir];
    for (int ic=0; ic<corr.GetNcols(); ++ic) {
      double eC=errs[ic];
      (*cov)(ir,ic) = corr(ir,ic)*(eR*eC);
    }
  }
  return cov;
}

//--------------------------------------------------

TMatrixD* covFromCorr(const TMatrixD &corr, const TH2D* h2Errs) {
  int nUnfBins=DYTools::nUnfoldingBins;
  if ((corr.GetNrows() != nUnfBins) ||
      (corr.GetNcols() != nUnfBins)) {
    std::cout << "covFromCorr(TMatrixD,TH2D*) dimension mismatch:\n";
    std::cout << "  corr[" << corr.GetNrows() << " x "
	      << corr.GetNcols() << "]\n";
    std::cout << "Expected DYTools::nUnfoldingBins=" << nUnfBins << "\n";
    return NULL;
  }
  if ((h2Errs->GetNbinsX() != DYTools::nMassBins) ||
      (h2Errs->GetNbinsY() != DYTools::nYBinsMax)) {
    std::cout << "covFromCorr(TMatrixD,TH2D*) dimension mismatch:\n";
    std::cout << "Expected h2Errs to be " << DYTools::nMassBins << "x"
	      << DYTools::nYBinsMax << ", got "
	      << h2Errs->GetNbinsX() << "x" << h2Errs->GetNbinsY() << "\n";
    return NULL;
  }

  TMatrixD *cov= new TMatrixD(corr);
  if (!cov) {
    std::cout << "covFromCorr: failed to create a new matrix\n";
    return NULL;
  }
  int iUnf=0;
  for (int ibin=1; ibin<=h2Errs->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=h2Errs->GetNbinsY(); ++jbin, ++iUnf) {
      double err= h2Errs->GetBinContent(ibin,jbin);
      if (iUnf>=nUnfBins) break;
      for (int i=0; i<nUnfBins; ++i) {
	(*cov)(iUnf,i) *= err;
	(*cov)(i,iUnf) *= err;
      }
    }
  }

  return cov;
}

//--------------------------------------------------

TMatrixD* partialCorrFromCov(const TMatrixD &totCov, const TMatrixD &cov) {
  TMatrixD *corr= new TMatrixD(cov);
  if (!corr) {
    std::cout << "partialCorrFromCov: failed to create a new matrix\n";
    return NULL;
  }
  corr->Zero();
  for (int ir=0; ir<cov.GetNrows(); ++ir) {
    double eR=sqrt(totCov(ir,ir));
    for (int ic=0; ic<cov.GetNcols(); ++ic) {
      double eC=sqrt(totCov(ic,ic));
      (*corr)(ir,ic) = cov(ir,ic)/(eR*eC);
    }
  }
  return corr;
}

//--------------------------------------------------

TMatrixD *relativeCov(const TVectorD &centralValue, const TMatrixD &cov) {
  TMatrixD *M = new TMatrixD(cov);
  if (!M) {
    std::cout << "relativeCov: failed to create a new matrix\n";
    return NULL;
  }
  M->Zero();
  for (int ir=0; ir<cov.GetNrows(); ++ir) {
    double cR=centralValue(ir);
    for (int ic=0; ic<cov.GetNcols(); ++ic) {
      double cC=centralValue(ic);
      (*M)(ir,ic) = cov(ir,ic)/(cR*cC);
    }
  }
  return M;
}

//--------------------------------------------------

TH2D* createHisto2D(const TMatrixD &M, const TMatrixD *Merr, const char *histoName, const char *histoTitle, TColorRange_t centerRange, int massBins, double maxValUser) {
  int nRows=M.GetNrows();
  int nCols=M.GetNcols();
  TH2D *h=NULL;

  if (massBins) {
#ifndef DYTools_HH
    std::cout << "createHisto2D. Warning: massBins=1, but DYTools is not available\n";
#endif
    if ((nRows==DYTools::nMassBins) && (nCols==1)) {
      double *mbins=new double[DYTools::nMassBins+1];
      for (int i=0; i<=DYTools::nMassBins; ++i)
	mbins[i]=DYTools::massBinLimits[i];
      //mbins[DYTools::nMassBins] += 5e-1;
      h= new TH2D(histoName,histoTitle,
		  DYTools::nMassBins, mbins,
		  1, DYTools::yRangeMin, DYTools::yRangeMax);
      delete mbins;
    }
    else if ((nRows!=DYTools::nMassBins) || (nCols!=DYTools::nMassBins)) {
      std::cout << "createHisto2D. Warning: massBins=1, DYTools is available, but the matrix has incorrect number of bins (" << nRows << "x" << nCols << "), instead of a square matrix of dim=" << DYTools::nMassBins << "\n";
    }
    else {
      double *mbins=new double[DYTools::nMassBins+1];
      for (int i=0; i<=DYTools::nMassBins; ++i)
	mbins[i]=DYTools::massBinLimits[i];
      //mbins[DYTools::nMassBins] += 5e-1;
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

//--------------------------------------------------

TMatrixD* createMatrixD(const TH2D* h2, int useErr) {
  int nRows=h2->GetNbinsX();
  int nCols=h2->GetNbinsY();
  TMatrixD *M=new TMatrixD(nRows,nCols);
  if (!M) {
    std::cout << "createMatrixD: failed to create a new matrix\n";
    return NULL;
  }
  M->Zero();

  for (int ir=0; ir<nRows; ir++) {
    for (int ic=0; ic<nCols; ic++) {
      double val=(useErr) ?
	h2->GetBinError(ir+1,ic+1) : h2->GetBinContent(ir+1,ic+1);
      (*M)(ir,ic)=val;
    }
  }
  return M;
}

//--------------------------------------------------

int divideMatrix(TMatrixD &nom, const TMatrixD &denom) {
  if (!sameNumBins(nom,denom)) {
    std::cout << "divideMatrix: different sizes\n";
    return 0;
  }
  for (int ir=0; ir<nom.GetNrows(); ++ir) {
    for (int ic=0; ic<nom.GetNcols(); ++ic) {
      double x= nom(ir,ic);
      double y= denom(ir,ic);
      double r= x/y;
      if ((y==double(0)) && (x==double(0))) r=0.;
      nom(ir,ic) = r;
    }
  }
  return 1;
}

//--------------------------------------------------

TH1D* createHisto1D(const TVectorD &vec, const TVectorD* vecErr,
		    const char *histoName, const char* histoTitle,
		    const TString &xAxisLabel, const TString &yAxisLabel) {
  TString name= (histoName) ? histoName : "histo1D";
  TString title=(histoTitle) ? histoTitle : name.Data();
  TH1D* h= createAnyTH1D(name,title,
			 vec.GetNoElements(),0.,vec.GetNoElements(),
			 xAxisLabel,yAxisLabel);
  for (int i=0; i<vec.GetNoElements(); i++) {
    h->SetBinContent(i+1, vec[i]);
    h->SetBinError  (i+1, (vecErr) ? (*vecErr)[i] : 0.);
  }
  return h;
}

//--------------------------------------------------
//--------------------------------------------------
/*
TH2D* getRelDifferenceVA(const TH2D *baseValue, TString newName, int nVariations, TH2D* hVar1, ...) {
  HERE("entered getRelDifference");
  TString newTitle=baseValue->GetTitle() + TString(" rel.Diff");
  TH2D *hRes=Clone(baseValue,newName,newTitle);
  HERE("ok");

  // remove error
  for (int ibin=1; ibin<=hRes->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=hRes->GetNbinsY(); ++jbin) {
      hRes->SetBinError(ibin,jbin, 0.);
    }
  }
  HERE("start variations");

  if (nVariations>=1) {
    va_list vl;
    va_start(vl,hVar1);
    for (int i=0; i<nVariations; ++i) {
      HERE("i=%d",i);
      TH2D* hv_inp=(TH2D*)va_arg(vl,TH2D*);
      if (!hv_inp) HERE("hv_inp is null");
      HERE("got hv_inp");
      //printHisto(hv_inp);
      hv_inp->Print("range");  // crashed here!
      TH2D* hv=(TH2D*)hv_inp->Clone(Form("var_%d",i));
      HERE("cloned ok %d",((hv==NULL)?0:1));
      hv->Add(baseValue,-1.);
      HERE("added base val");
      for (int ibin=1; ibin<=hRes->GetNbinsX(); ++ibin) {
	for (int jbin=1; jbin<=hRes->GetNbinsY(); ++jbin) {
	  double err=hRes->GetBinError(ibin,jbin);
	  double var=fabs(hv->GetBinError(ibin,jbin));
	  if (var>err) hRes->SetBinError(ibin,jbin, var);
	}
      }
      delete hv;
    }
    va_end(vl);
  }
  return hRes;
}
*/

//--------------------------------------------------

TH2D* getRelDifference(const TH2D *baseValue, TString newName, int includeVariants, const TH2D* hVar1, const TH2D *hVar2, const TH2D *hVar3, const TH2D* hVar4) {
  if (!baseValue) {
    std::cout << "getRelDifference: base histogram is NULL" << std::endl;
    return NULL;
  }

  std::vector<const TH2D*> hV;
  hV.reserve(4);
  hV.push_back(hVar1);
  if (hVar2!=NULL) hV.push_back(hVar2);
  if (hVar3!=NULL) hV.push_back(hVar3);
  if (hVar4!=NULL) hV.push_back(hVar4);

  TH2D *h2Diff=getRelDifference(hV,newName,includeVariants);
  if (!h2Diff) std::cout << "error in getRelDifference(TH2D*,etc)\n";
  return h2Diff;
}

//--------------------------------------------------

TH2D* getRelDifference(const std::vector<const TH2D*> &var, TString newName, int includeVariants) {
  if (var.size()==0) {
    std::cout << "getRelDifference(vec): vector is empty" << std::endl;
    return NULL;
  }

  TH2D *hRes=Clone(var[0],newName,newName);
  if (!hRes) {
    std::cout << "getRelDifference(vec): failed to copy the 1st histogram" << std::endl;
    return NULL;
  }

  // remove error
  removeError2D(hRes);

  // If includeVariants is 1, then differences between all histograms are
  // considered. The j index will run from 0 (baseValue) to
  // hV.size(). If j>0, the reference is the variant histogram
  unsigned int imax=(includeVariants) ? (var.size()+1) : 1;

  for (unsigned int i=0; i<imax; ++i) {
    for (unsigned int j=i+1; j<var.size(); ++j) {
      for (int ibin=1; ibin<=hRes->GetNbinsX(); ++ibin) {
	for (int jbin=1; jbin<=hRes->GetNbinsY(); ++jbin) {
	  // current max difference
	  double err=hRes->GetBinError(ibin,jbin);
	  // get diff
	  double diff=fabs(var[j]->GetBinContent(ibin,jbin) -
			   var[i]->GetBinContent(ibin,jbin));
	  if (diff>err) hRes->SetBinError(ibin,jbin, diff);
	}
      }
    }
  }
  return hRes;
}

//--------------------------------------------------

TH2D* getRelDifference(const std::vector<TH2D*> &var, TString newName, int includeVariants) {
  std::vector<const TH2D*> vec;
  for (unsigned int i=0; i<var.size(); ++i) {
    vec.push_back((const TH2D*)var[i]);
  }
  TH2D *diff=getRelDifference(vec,newName,includeVariants);
  return diff;
}

//--------------------------------------------------

TH2D* getRelError(const TH2D* h2, TString newName, int storeAsError) {
  TH2D* h2r=Clone(h2,newName);
  h2r->Reset();
  for (int ibin=1; ibin<=h2->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=h2->GetNbinsY(); ++jbin) {
      double v=h2->GetBinContent(ibin,jbin);
      double e=h2->GetBinError  (ibin,jbin);
      if (storeAsError) h2r->SetBinError(ibin,jbin, e/v);
      else h2r->SetBinContent(ibin,jbin, e/v);
    }
  }
  return h2r;
}

//--------------------------------------------------
//--------------------------------------------------

// write flags to a file.
// Assumes that the file is open
int writeIntFlagValues(const TString &fieldName, int nFlags, int flag1, ...) {
  TVectorD *v= new TVectorD(nFlags);
  (*v)(0)=flag1;
  if (nFlags>1) {
    va_list vl;
    va_start(vl,flag1);
    // one flag is already specified
    for (int i=0; i<nFlags-1; i++) {
    (*v)(i+1)=va_arg(vl,int);
    }
    va_end(vl);
  }
  v->Write(fieldName);
  delete v;
  return 1;
}

//--------------------------------------------------

// write flags to a file.
// Assumes that the file is open
int writeFlagValues(const TString &fieldName, int nFlags, double flag1, ...) {
  TVectorD *v= new TVectorD(nFlags);
  (*v)(0)=flag1;
  if (nFlags>1) {
    va_list vl;
    va_start(vl,flag1);
    // one flag is already specified
    for (int i=0; i<nFlags-1; i++) {
    (*v)(i+1)=va_arg(vl,double);
    }
    va_end(vl);
  }
  v->Write(fieldName);
  delete v;
  return 1;
}

//---------------------------------------------------------------

TVectorD* readFlagValues(TFile &fin, const TString &fieldName, int nFlags) {
  int doWarn=(nFlags>0) ? 1:0;
  TVectorD* v=(TVectorD*)fin.Get(fieldName);
  int got=(v) ? v->GetNoElements() : 0;
  if (v && (got!=abs(nFlags))) { delete v; v=NULL; }
  if (!v && doWarn) {
    std::cout << "readFlagValues: couldn't get full field=<"
	      << fieldName << ">. Got " << got << "/" << nFlags
	      << " elements\n";
  }
  return v;
}

//---------------------------------------------------------------

int writeFlagValues(const TString &fieldName, int count, const double *flags) {
  if (!flags) {
    std::cout  << "Error writeFlagValues(fieldName=" << fieldName << ","
	       << count << ",double*): null ptr\n";
    return 0;
  }
  TVectorD v(count);
  for (int i=0; i<count; ++i) {
    v(i)=flags[i];
  }
  v.Write(fieldName);
  return 1;
}

//---------------------------------------------------------------

int writeFlagValues(const TString &fieldName, const std::vector<int> &flagsV) {
  TVectorD v(flagsV.size());
  for (unsigned int i=0; i<flagsV.size(); ++i) {
    v(i)=flagsV[i];
  }
  v.Write(fieldName);
  return 1;
}

//---------------------------------------------------------------

int writeFlagValues(const TString &fieldName,
		    const std::vector<double> &flagsV) {
  TVectorD v(flagsV.size());
  for (unsigned int i=0; i<flagsV.size(); ++i) {
    v(i)=flagsV[i];
  }
  v.Write(fieldName);
  return 1;
}

//------------------------------------------------------------------------------------------------------------------------

//--------------------------------------------------
//--------------------------------------------------

TString getTimeStrForPrint(Float_t tf) {
  int t=int(tf);
  int mins=t/60;
  int hours=mins/60;
  TString str;
  if (hours>0) {
    mins= t%60;
    str=Form("(%2d hour(s) %2d min)",hours,mins);
  }
  else str=Form("(%2d min)",mins);
  return str;
}

//-----------------------

void ShowBenchmarkTime(const char *clock_name) {
  gBenchmark->Show(clock_name);
  Float_t realT=gBenchmark->GetRealTime(clock_name);
  Float_t cpuT =gBenchmark->GetCpuTime(clock_name);
  printf("  realTime= %5.2lf sec %s\n",realT,getTimeStrForPrint(realT).Data());
  printf("  cpuTime = %5.2lf sec %s\n",cpuT ,getTimeStrForPrint(cpuT ).Data());
  std::cout.flush();
}

//--------------------------------------------------
//--------------------------------------------------


//  convert m[nMassBins][ybins] -> v[flat_idx]
TH1D* flattenHisto(const TH2D *h2, TString newName) {
  TH1D* h=createAnyTH1D(newName,newName,DYTools::nUnfoldingBins,
			0.,DYTools::nUnfoldingBins,"bin","count");
  int fi=0;
  for (int i=0;
       (i<DYTools::nMassBins) && (fi<DYTools::nUnfoldingBins);
       ++i) {
    for (int yi=0;
	 (yi<DYTools::nYBins[i]) && (fi<DYTools::nUnfoldingBins);
	 ++yi, ++fi) {
      h->SetBinContent(fi+1, h2->GetBinContent(i+1,yi+1));
      h->SetBinError  (fi+1, h2->GetBinError  (i+1,yi+1));
    }
  }
  return h;
}

// -------------------------------------------

//  convert v[flat_idx] -> m[nMassBins][ybins]
TH2D* deflattenHisto(const TH1D *h, TString newName) {
  TH2D* h2= createBaseH2(newName,newName,1);
  int fi=0;
  for (int i=0;
       (i<DYTools::nMassBins) && (fi<DYTools::nUnfoldingBins);
       ++i) {
    for (int yi=0;
	 (yi<DYTools::nYBins[i]) && (fi<DYTools::nUnfoldingBins);
	 ++yi, ++fi) {
      h2->SetBinContent(i+1,yi+1, h->GetBinContent(fi+1));
      h2->SetBinError  (i+1,yi+1, h->GetBinError  (fi+1));
    }
  }
  return h2;
}

//--------------------------------------------------
//--------------------------------------------------

TH2D* convertBaseH2actual(const TH2D* h2, TString newHistoName, int setTitle) {
  TH2D* newH=Clone(h2,newHistoName,newHistoName);
  if (!setTitle) newH->SetTitle("");
  for (int ibin=1; ibin<=h2->GetNbinsX(); ++ibin) {
    // adjustment is needed only if the bin count is different
    if (DYTools::nYBins[ibin-1]!=DYTools::nYBinsMax) {
      // ready only for double spacing
      if (2*DYTools::nYBins[ibin-1]!=DYTools::nYBinsMax) {
	std::cout << "convertBaseH2actual is ready only for double spacing\n";
	delete newH;
	return NULL;
      }
      for (int jbin=1; jbin<=h2->GetNbinsY(); jbin+=2) {
	double x=
	  h2->GetBinContent(ibin,jbin  ) +
	  h2->GetBinContent(ibin,jbin+1);
	double e2=
	  pow(h2->GetBinError(ibin,jbin  ),2) +
	  pow(h2->GetBinError(ibin,jbin+1),2);
	double e=sqrt(e2);
	int newJbin=(jbin+1)/2;
	newH->SetBinContent(ibin, newJbin, x);
	newH->SetBinError  (ibin, newJbin, e);
      }
      for (int jbin=DYTools::nYBins[ibin-1]+1; jbin<=DYTools::nYBinsMax; ++jbin) {
	newH->SetBinContent(ibin,jbin,  0.);
	newH->SetBinError  (ibin,jbin,  1.e9);
      }
    }
  }
  return newH;
}

//--------------------------------------------------

TH1D* convertBinNo2Mass(const TH1D* histo) {
  if (histo->GetNbinsX() != DYTools::nMassBins) {
    std::cout << "convertBinNo2Mass: histo->GetName="
	      << histo->GetName() << " has " << histo->GetNbinsX()
	      << " bins instead of "
	      << DYTools::nMassBins << "\n";
    return NULL;
  }
  TH1D *h1=createBaseH1(TString(histo->GetName()) + TString("_massClone"));
  if (!h1 ||
      !Copy(histo,h1)) {
    std::cout << "convertBinNo2Mass: failed to prepare\n";
    return NULL;
  }

  return h1;
}

//--------------------------------------------------

TH1D* convertRange2BinNo(const TH1D* histo) {
  TString name=histo->GetName() + TString("_binNo");
  int nBins= histo->GetNbinsX();
  TH1D *h1=new TH1D(name,name, nBins, -0.5, nBins-0.5);
  if (!h1) {
    std::cout << "convertRange2BinNo(TH1D): failed to create histo\n";
    return NULL;
  }
  h1->SetDirectory(0);
  h1->SetStats(0);
  if (!Copy(histo,h1)) {
    std::cout << "convertRange2BinNo: failed to copy\n";
    return NULL;
  }
  return h1;
}

//--------------------------------------------------

TH2D* convertRange2BinNo(const TH2D* histo) {
  TString name=histo->GetName() + TString("_binNo");
  int nBinsX= histo->GetNbinsX();
  int nBinsY= histo->GetNbinsY();
  TH2D *h2=new TH2D(name,name, nBinsX, -0.5, nBinsX-0.5, nBinsY, -0.5, nBinsY-0.5);
  if (!h2) {
    std::cout << "convertRange2BinNo(TH2D): failed to create histo\n";
    return NULL;
  }
  h2->SetDirectory(0);
  h2->SetStats(0);
  if (!Copy(histo,h2)) {
    std::cout << "convertRange2BinNo(TH2D): failed to copy\n";
    return NULL;
  }
  return h2;
}

//--------------------------------------------------

TH1D* createProfileX(TH2D *h2, int iyBin, const TString &name, int setTitle, const char *title) {
  if ((iyBin<=0) || (iyBin>h2->GetNbinsY())) {
    std::cout << "\n\n\tcreateProfileX(" << h2->GetName() << ", iyBin=" << iyBin << "(bad value!!), name=" << name << ")\n\n";
  }

  // prepare the range info
  int nxBins=h2->GetNbinsX();
  double *xv=new double[nxBins+1];
  TAxis *ax=h2->GetXaxis();
  for (int i=1; i<=nxBins; i++) {
    xv[i-1] = ax->GetBinLowEdge(i);
  }
  xv[nxBins]=ax->GetBinLowEdge(nxBins) + ax->GetBinWidth(nxBins);
  TH1D *h=new TH1D(name,"",h2->GetNbinsX(),xv);
  delete [] xv;

  // copy the profile
  h->SetDirectory(0);
  h->SetStats(0);
  if (setTitle) {
    if (title) h->SetTitle(title);
    else h->SetTitle(name);
  }
  h->GetXaxis()->SetTitle( h2->GetXaxis()->GetTitle() );
  for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
    h->SetBinContent(ibin,h2->GetBinContent(ibin,iyBin));
    h->SetBinError(ibin,h2->GetBinError(ibin,iyBin));
  }
  return h;
}
// -------------------------------------------
/*
inline
TH1D* createProfileY(TH2D *h2, int ixBin, const TString &name, int setTitle=0, const char *title=NULL) {
  if ((ixBin<=0) || (ixBin>h2->GetNbinsX())) {
    std::cout << "\n\n\tcreateProfileY(" << h2->GetName() << ", ixBin=" << ixBin << "(bad value!!), name=" << name << ")\n\n";
  }
  TH1D *h=new TH1D(name,"",h2->GetNbinsY(),h2->GetYaxis()->GetXbins()->GetArray());
  h->SetDirectory(0);
  h->SetStats(0);
  if (setTitle) {
    if (title) h->SetTitle(title);
    else h->SetTitle(name);
  }
  h->GetXaxis()->SetTitle( h2->GetYaxis()->GetTitle() );
  for (int ibin=1; ibin<=h2->GetNbinsY(); ibin++) {
    h->SetBinContent(ibin,h2->GetBinContent(ixBin,ibin));
    h->SetBinError(ibin,h2->GetBinError(ixBin,ibin));
  }
  return h;
}
*/

// -------------------------------------------

TH1D* createProfileY(TH2D *h2, int ixBin, const TString &name, int setTitle, const char *title, int set_nYbins, double set_ymin, double set_ymax) {
  if ((ixBin<=0) || (ixBin>h2->GetNbinsX())) {
    std::cout << "\n\n\tcreateProfileY(" << h2->GetName() << ", ixBin=" << ixBin << "(bad value!!), name=" << name << ")\n\n";
  }
  TH1D *h=NULL;
  if (set_nYbins==-1) {
    int nyBins=h2->GetNbinsY();
    double *yv=new double[nyBins+1];
    TAxis *ay=h2->GetYaxis();
    for (int i=1; i<=nyBins; i++) {
      yv[i-1] = ay->GetBinLowEdge(i);
      //std::cout << "yv[" << i-1 << "]=" << yv[i-1] << "\n";
    }
    yv[nyBins]=ay->GetBinLowEdge(nyBins) + ay->GetBinWidth(nyBins);
    std::cout << "yv[nYBins=" << nyBins << "]=" << yv[nyBins] << "\n";
    h=new TH1D(name,"",h2->GetNbinsY(),yv);
    delete [] yv;
    set_nYbins=h2->GetNbinsY();
  }
  else h=new TH1D(name,"",set_nYbins,set_ymin,set_ymax);

  h->SetDirectory(0);
  h->SetStats(0);
  if (setTitle) {
    if (title) h->SetTitle(title);
    else h->SetTitle(name);
  }
  h->GetXaxis()->SetTitle( h2->GetYaxis()->GetTitle() );
  for (int ibin=1; (ibin<=h2->GetNbinsY()) && (ibin<=set_nYbins); ibin++) {
    h->SetBinContent(ibin,h2->GetBinContent(ixBin,ibin));
    h->SetBinError(ibin,h2->GetBinError(ixBin,ibin));
  }
  return h;
}

//--------------------------------------------------

TH1D* createProfileAuto(TH2D* h2, int iBin, const TString &name, int setTitle, const char *title) {
  TH1D *hProf= (DYTools::study2D) ?
    createProfileY(h2,iBin,name,setTitle,title) :
    createProfileX(h2,iBin,name,setTitle,title);
  if (!hProf) std::cout << "error from createProfileAuto\n";
  return hProf;
}

//--------------------------------------------------

int createRapidityProfileVec(const std::vector<TH2D*> &h2SrcV, std::vector<std::vector<TH1D*>*> &hProfV, const std::vector<TString> &labelsV, int markerStyle, double markerSize) {
  if (h2SrcV.size()==0) {
    std::cout << "createRapidityProfileVec: h2SrcV.size=0\n";
    return 0;
  }

  int res=1;
  for (int im=0; res && (im<DYTools::nMassBins); im++) {
    std::vector<TH1D*> *hPVec = new std::vector<TH1D*>();
    hProfV.push_back(hPVec);

    TString mStr=Form("_m%1.0f_%1.0f",DYTools::massBinLimits[im],
		      DYTools::massBinLimits[im+1]);

    for (unsigned int i=0; res && (i<h2SrcV.size()); ++i) {
      TH1D* hProf=NULL;
      if (h2SrcV[i]) {
	TString name=Form("hProf_%s_%s",labelsV[i].Data(),mStr.Data());
	//std::cout << "Creating hProfRaw=" << name << "\n";
	hProf= createProfileY(h2SrcV[i],im+1,name,1,name, DYTools::nYBins[im],0.,DYTools::yRangeMax+1e-3);
	if (!hProf) res=0;
	else {
	  hProf->SetMarkerStyle(markerStyle);
	  hProf->SetMarkerSize(markerSize);
	}
      }
      hPVec->push_back(hProf);
    }
  }
  if (!res) std::cout << "error in createRapidityProfileVec\n";
  return res;
}

//--------------------------------------------------

int createMassProfileVec(const std::vector<TH2D*> &h2SrcV, std::vector<std::vector<TH1D*>*> &hProfV, const std::vector<TString> &labelsV, int markerStyle, double markerSize) {
  if (h2SrcV.size()==0) {
    std::cout << "createMassProfileVec: h2SrcV.size=0\n";
    return 0;
  }

  if (DYTools::study2D==1) {
    int yCount=DYTools::nYBins[0];
    for (int im=1; im<DYTools::nMassBins; ++im) {
      if (yCount!=DYTools::nYBins[im]) {
	std::cout << "createMassProfileVec: the rapidity bin count should be the same in all mass bins\n";
	return 0;
      }
    }
  }


  TMatrixD *ybinLimits=DYTools::getYBinLimits();
  int res=(ybinLimits) ? 1:0;

  for (int iy=0; res && (iy<DYTools::nYBins[0]); iy++) {
    std::vector<TH1D*> *hPVec = new std::vector<TH1D*>();
    hProfV.push_back(hPVec);

    TString yStr=Form("_y%1.0f_%1.0f",(*ybinLimits)(0,iy),(*ybinLimits)(0,iy+1));

    for (unsigned int i=0; res && (i<h2SrcV.size()); ++i) {
      TH1D* hProf=NULL;
      if (h2SrcV[i]) {
	TString name=Form("hProf_%s_%s",labelsV[i].Data(),yStr.Data());
	std::cout << "Creating hProfRaw=" << name << "\n";
	hProf= createProfileX(h2SrcV[i],iy+1,name,1,name.Data());
	if (!hProf) res=0;
	else {
	  hProf->SetMarkerStyle(markerStyle);
	  hProf->SetMarkerSize(markerSize);
	}
      }
      hPVec->push_back(hProf);
    }
  }

  if (ybinLimits) delete ybinLimits;
  if (!res) std::cout << "error in createMassProfileVec\n";
  return res;
}

//--------------------------------------------------
//--------------------------------------------------

TH1D* removeLastBin(const TH1D* hOrig, TString newName, int setTitle, const char *newTitle) {
  // prepare the range info
  int nxBins=hOrig->GetNbinsX()-1; // we will remove last bin
  double *xv=new double[nxBins+1];
  TAxis *ax=hOrig->GetXaxis();
  for (int i=1; i<=nxBins; i++) {
    xv[i-1] = ax->GetBinLowEdge(i);
  }
  xv[nxBins]=ax->GetBinLowEdge(nxBins) + ax->GetBinWidth(nxBins); // last bin  right edge

  TH1D *h=new TH1D(newName,newName,nxBins,xv);
  delete [] xv;

  // copy the profile
  h->SetDirectory(0);
  h->SetStats(0);
  if (setTitle) {
    if (newTitle) h->SetTitle(newTitle);
    else h->SetTitle(newName);
  }
  h->GetXaxis()->SetTitle( hOrig->GetXaxis()->GetTitle() );
  for (int ibin=1; ibin<=nxBins; ibin++) {
    h->SetBinContent(ibin,hOrig->GetBinContent(ibin));
    h->SetBinError(ibin,hOrig->GetBinError(ibin));
  }
  return h;
}

//--------------------------------------------------

void linearApprox(TH1D *h, int ibin1, int ibin2,
		  double valOnEdge1, double valOnEdge2) {
  double y1=h->GetBinContent(ibin1);
  double y2=h->GetBinContent(ibin2);
  if (valOnEdge1!=double(-1.)) y1=valOnEdge1;
  if (valOnEdge2!=double(-1.)) y2=valOnEdge2;
  double dy=(ibin1==ibin2) ? 0. : (y2-y1)/(ibin2-ibin1);
  for (int ib=ibin1; ib<=ibin2; ++ib) {
    h->SetBinContent(ib, y1 + dy*(ib-ibin1));
  }
}

//--------------------------------------------------

void linearApprox(TH1D *h, double mass1, double mass2,
		  double valOnEdge1, double valOnEdge2) {
  linearApprox(h,h->FindBin(mass1),h->FindBin(mass2),
	       valOnEdge1,valOnEdge2);
}

//--------------------------------------------------
//--------------------------------------------------

int convertBaseH2actualVec(const std::vector<TH2D*> &baseV, std::vector<TH2D*> &actualV, const TString histoNameBase, const std::vector<TString> &sample_labels, int setHistoTitle) {
  int res=1;
  if (baseV.size()==0) return 0;
  for (unsigned int i=0; i<baseV.size(); ++i) {
    TString histoName= histoNameBase + sample_labels[i];
    TH2D* h=convertBaseH2actual(baseV[i],histoName,setHistoTitle);
    if (!h) { res=0; break; }
    actualV.push_back(h);
  }
  if (!res) std::cout << "error in convertBaseH2actualVec\n";
  return res;
}

//--------------------------------------------------
//--------------------------------------------------

// error of denom is ignored
int scaleHisto(TH1D *histoNom, const TH1D *histoDenom, int mult) {
  if (histoNom->GetNbinsX() != histoDenom->GetNbinsX()) {
    std::cout << "scaleHisto: different number of bins\n";
    assert(0);
    return 0;
  }
  for (int ibin=1; ibin<=histoNom->GetNbinsX(); ++ibin) {
    double v=histoNom->GetBinContent(ibin);
    double e=histoNom->GetBinError(ibin);
    double x=histoDenom->GetBinContent(ibin);
    if (mult) {
      histoNom->SetBinContent(ibin, v*x);
      histoNom->SetBinError(ibin, e*x);
    }
    else {
      histoNom->SetBinContent(ibin, v/x);
      histoNom->SetBinError(ibin, e/x);
    }
  }
  return 1;
}

//--------------------------------------------------

// error of denom is ignored
int scaleHisto(TH2D *histoNom, const TH2D *histoDenom, int mult) {
  if ((histoNom->GetNbinsX() != histoDenom->GetNbinsX()) ||
      (histoNom->GetNbinsY() != histoDenom->GetNbinsY())) {
    std::cout << "scaleHisto(2D): different number of bins\n";
    return 0;
  }
  for (int ibin=1; ibin<=histoNom->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=histoNom->GetNbinsY(); ++jbin) {
      double v=histoNom->GetBinContent(ibin,jbin);
      double e=histoNom->GetBinError(ibin,jbin);
      double x=histoDenom->GetBinContent(ibin,jbin);
      if (mult) {
	histoNom->SetBinContent(ibin,jbin, v*x);
	histoNom->SetBinError(ibin,jbin, e*x);
      }
      else {
	histoNom->SetBinContent(ibin,jbin, v/x);
	histoNom->SetBinError(ibin,jbin, e/x);
      }
    }
  }
  return 1;
}

//--------------------------------------------------

// error of denom is NOT ignored
int multiplyHisto(TH1D *histoNom, const TH1D *histoDenom, int mult) {
  if (histoNom->GetNbinsX() != histoDenom->GetNbinsX()) {
    std::cout << "scaleHisto: different number of bins\n";
    assert(0);
    return 0;
  }
  for (int ibin=1; ibin<=histoNom->GetNbinsX(); ++ibin) {
    double v=histoNom->GetBinContent(ibin);
    double ve=histoNom->GetBinError(ibin);
    double x=histoDenom->GetBinContent(ibin);
    double xe=histoDenom->GetBinError(ibin);
    double relErr=sqrt((ve/v)*(ve/v) + (xe/x)*(xe/x));
    double result=(mult) ? v*x : v/x;
    histoNom->SetBinContent(ibin, result);
    histoNom->SetBinError(ibin, result*relErr);
  }
  return 1;
}

//--------------------------------------------------

// error of denom is NOT ignored
int multiplyHisto(TH2D *histoNom, const TH2D *histoDenom, int mult) {
  if ((histoNom->GetNbinsX() != histoDenom->GetNbinsX()) ||
      (histoNom->GetNbinsY() != histoDenom->GetNbinsY())) {
    std::cout << "scaleHisto(2D): different number of bins\n";
    return 0;
  }
  for (int ibin=1; ibin<=histoNom->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=histoNom->GetNbinsY(); ++jbin) {
      double v=histoNom->GetBinContent(ibin,jbin);
      double ve=histoNom->GetBinError(ibin,jbin);
      double x=histoDenom->GetBinContent(ibin,jbin);
      double xe=histoDenom->GetBinError(ibin,jbin);
      double relErr=sqrt((ve/v)*(ve/v) + (xe/x)*(xe/x));
      double result=(mult) ? v*x : v/x;
      histoNom->SetBinContent(ibin,jbin, result);
      histoNom->SetBinError(ibin,jbin, result*relErr);
    }
  }
  return 1;
}

//--------------------------------------------------
//--------------------------------------------------

template<class histo_t>
double* getXrange(const histo_t *h) {
  double *arr=new double[h->GetNbinsX()+1];
  for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin) {
    double x1=h->GetBinLowEdge(ibin);
    arr[ibin-1]=x1;
  }
  int i=h->GetNbinsX();
  double x=h->GetBinLowEdge(i);
  double w=h->GetBinWidth(i);
  arr[i]=x+w;
  return arr;
}

// ----------------------------------------------

template<class histo_t>
double* getYrange(const histo_t *h) {
  double *arr=new double[h->GetNbinsY()+1];
  TAxis *ay=h->GetYaxis();
  for (int ibin=1; ibin<=h->GetNbinsY(); ++ibin) {
    double x1=ay->GetBinLowEdge(ibin);
    arr[ibin-1]=x1;
  }
  int i=h->GetNbinsY();
  double x=ay->GetBinLowEdge(i);
  double w=ay->GetBinWidth(i);
  arr[i]=x+w;
  return arr;
}

//--------------------------------------------------

TH1D* convert_TH1F_to_TH1D(const TH1F *h, TString newName) {
  double *arr=getXrange(h);
  int size=h->GetNbinsX();
  TH1D *hd=new TH1D(newName,newName,size,arr);
  hd->SetDirectory(0);
  delete arr;
  for (int ibin=1; ibin<=size; ++ibin) {
    hd->SetBinContent(ibin, h->GetBinContent(ibin));
    hd->SetBinError(ibin, h->GetBinError(ibin));
  }
  return hd;
}

//--------------------------------------------------

TH1F* convert_TH1D_to_TH1F(const TH1D *h, TString newName) {
  double *arr=getXrange(h);
  TH1F *hf=new TH1F(newName,newName,h->GetNbinsX(),arr);
  hf->SetDirectory(0);
  delete arr;
  for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin) {
    hf->SetBinContent(ibin, h->GetBinContent(ibin));
    hf->SetBinError(ibin, h->GetBinError(ibin));
  }
  return hf;
}

//--------------------------------------------------

TH2D* convert_TH2F_to_TH2D(const TH2F *h, TString newName) {
  double *arrX=getXrange(h);
  double *arrY=getYrange(h);
  int sizeX=h->GetNbinsX();
  int sizeY=h->GetNbinsY();
  TH2D *hd=new TH2D(newName,newName,sizeX,arrX,sizeY,arrY);
  hd->SetDirectory(0);
  delete arrX;
  delete arrY;
  for (int ibin=1; ibin<=sizeX; ++ibin) {
    for (int jbin=1; jbin<=sizeY; ++jbin) {
      hd->SetBinContent(ibin,jbin, h->GetBinContent(ibin,jbin));
      hd->SetBinError(ibin,jbin, h->GetBinError(ibin,jbin));
    }
  }
  return hd;
}

//--------------------------------------------------

TH2F* convert_TH2D_to_TH2F(const TH2D *h, TString newName) {
  double *arrX=getXrange(h);
  double *arrY=getYrange(h);
  int sizeX=h->GetNbinsX();
  int sizeY=h->GetNbinsY();
  TH2F *hf=new TH2F(newName,newName,sizeX,arrX,sizeY,arrY);
  hf->SetDirectory(0);
  delete arrX;
  delete arrY;
  for (int ibin=1; ibin<=sizeX; ++ibin) {
    for (int jbin=1; jbin<=sizeY; ++jbin) {
      hf->SetBinContent(ibin,jbin, h->GetBinContent(ibin,jbin));
      hf->SetBinError(ibin,jbin, h->GetBinError(ibin,jbin));
    }
  }
  return hf;
}

//--------------------------------------------------
//--------------------------------------------------

TH1D* convert_absolute2rshape(const TH1D *h, double norm) {
  if (DYTools::study2D) {
    std::cout << "convert_absolute2rshape: 1D study is expected\n";
    return NULL;
  }

  TString name=h->GetName() + TString("_rshape");
  TH1D *hR=Clone(h,name,name);

  for (int ibin=1; ibin<=DYTools::nMassBins; ++ibin) {
    double dM=DYTools::massBinLimits[ibin] - DYTools::massBinLimits[ibin-1];
    double v= h->GetBinContent(ibin);
    double e= h->GetBinError  (ibin);
    hR->SetBinContent(ibin, v/(norm*dM));
    hR->SetBinError  (ibin, e/(norm*dM));
  }
  return hR;
}

//--------------------------------------------------

TH1D* convert_rshape2absolute(const TH1D *h, double norm) {
  if (DYTools::study2D) {
    std::cout << "convert_rshape2absolute: 1D study is expected\n";
    return NULL;
  }

  TString name=h->GetName() + TString("_absolute");
  TH1D *hR=Clone(h,name,name);

  for (int ibin=1; ibin<=DYTools::nMassBins; ++ibin) {
    double dM=DYTools::massBinLimits[ibin] - DYTools::massBinLimits[ibin-1];
    double v= h->GetBinContent(ibin);
    double e= h->GetBinError  (ibin);
    hR->SetBinContent(ibin, v*dM*norm);
    hR->SetBinError  (ibin, e*dM*norm);
  }
  return hR;
}

//--------------------------------------------------
//--------------------------------------------------

TH1D* LoadHisto1D(TString histoName, const TString &fname, TString subDir, int checkBinning) {
  TString theCall=TString("LoadHisto1D(<") + histoName + TString(">,<") + fname + TString(">,<") + subDir + TString(Form(">, checkBinning=%d)",checkBinning));

  TFile fin(fname,"read");
  if (!fin.IsOpen()) {
    std::cout << theCall << ": failed to open the file\n";
    return NULL;
  }
  if (checkBinning && !checkBinningArrays(fin)) {
    std::cout << theCall << ": binning test failed\n";
    return NULL;
  }

  if (histoName.Index("/")!=-1) {
    Ssiz_t idx=histoName.Index("/");
    subDir.Append(histoName(0,idx+1));
    histoName.Remove(0,idx+1);
  }

  if (subDir.Length()) { fin.cd(subDir); }
  TH1D *h1=(TH1D*)fin.Get(histoName);
  if (h1) h1->SetDirectory(0);
  fin.Close();

  if (!h1) std::cout << "failed to get histo " << histoName
		     << ", from subDir=<" << subDir << ">\n";
  return h1;
}

//--------------------------------------------------

TH2D* LoadHisto2D(TString histoName, const TString &fname, TString subDir, int checkBinning) {
  TString theCall=TString("LoadHisto2D(<") + histoName + TString(">,<") + fname + TString(">,<") + subDir + TString(Form(">, checkBinning=%d)",checkBinning));

  TFile fin(fname,"read");
  if (!fin.IsOpen()) {
    std::cout << theCall << ": failed to open the file\n";
    return NULL;
  }
  if (checkBinning && !checkBinningArrays(fin)) {
    std::cout << theCall << ": binning test failed\n";
    return NULL;
  }

  if (histoName.Index("/")!=-1) {
    Ssiz_t idx=histoName.Index("/");
    subDir.Append(histoName(0,idx+1));
    histoName.Remove(0,idx+1);
  }

  //std::cout << "calling LoadHisto2D(fin," << histoName << "," << subDir << ")\n";
  TH2D* h2=LoadHisto2D(fin,histoName,subDir,0);
  if (!h2) std::cout << theCall << ": failed to load the histo\n";
  else {
    h2->SetDirectory(0);
  }
  fin.Close();
  return h2;
}

//--------------------------------------------------

TH2D* LoadHisto2D(TFile &fin, TString histoName, TString subDir, int checkBinning) {
  TString theCall=TString("LoadHisto2D(file=<") + fin.GetName() + TString(">,<") + histoName + TString(">,<") + subDir + TString(Form(">, checkBinning=%d)",checkBinning));

  if (checkBinning && !checkBinningArrays(fin)) {
    std::cout << theCall << ": binning test failed\n";
    return NULL;
  }

  TString loadHistoName;
  if (subDir.Length()) {
    if (subDir[subDir.Length()-1]!='/') subDir.Append("/");
    loadHistoName=subDir + histoName;
  }
  else loadHistoName=histoName;
  //std::cout << "loadHistoName=<" << loadHistoName << ">\n";

  TH2D* h2=(TH2D*)fin.Get(loadHistoName);
  if (!h2) {
    std::cout << theCall << ": failed to load the histo <"
	      << loadHistoName << ">\n";
  }
  else {
    h2->SetDirectory(0);
  }
  return h2;
}

//--------------------------------------------------
//--------------------------------------------------

void writeBinningArrays(TFile &fout, TString producedBy, int allInfo) {
  if (!fout.IsOpen()) {
    std::cout << "writeBinningArrays error: file is not open\n";
    return;
  }
  fout.cd();
  int correction= (DYTools::study2D==-1) ? 1:0;
  if (allInfo || (DYTools::study2D)) {
    writeIntFlagValues("DYTools_study2D",1,DYTools::study2D);
  }
  else {
    writeIntFlagValues("DYTools_study1D",1,1);
  }
  TVectorD mass(DYTools::nMassBins+1);
  TVectorD rapidityCounts(DYTools::nMassBins + correction);
  if (correction) {
    mass[0]=0;
    rapidityCounts[0]=0;
  }
  else {
    for (int i=0; i<DYTools::nMassBins+1; i++)
      mass[i]=DYTools::massBinLimits[i];
    for (int i=0; i<DYTools::nMassBins  ; i++)
      rapidityCounts[i]=DYTools::nYBins[i];
  }
  mass.Write("massBinLimits");
  rapidityCounts.Write("rapidityCounts");

  // meta data
  TObjString info((producedBy.Length()) ?
		  producedBy : TString("Non-specified macro"));
  TObjString timeTag(DayAndTimeTag(0));
  TObjString explain="productionTime";
  if (allInfo) {
    info.Write("producedBy");
    timeTag.Write("timeTag");
  }
  info.Write(TString("infoProducedBy: ") + producedBy);
  explain.Write(timeTag.String());
}

//--------------------------------------------------

int checkBinningArrays(TFile &fin, int printMetaData) {
  const char *fncname="checkBinningArrays: ";
  TString fileInfo=TString("on file <") + fin.GetName() + ">";
  fin.cd();
  if (printMetaData) {
    std::cout << "Meta data " << fileInfo << "\n";
    TObjString *info=(TObjString*)fin.Get("producedBy");
    TObjString *timeTag=(TObjString*)fin.Get("timeTag");
    if (info) {
      std::cout << " - producedBy <" << info->String() << ">\n";
      delete info;
    }
    if (timeTag) {
      std::cout << " - timeTag <" << timeTag->String() << ">\n";
      delete timeTag;
    }
  }
  TVectorD *initOk=readFlagValues(fin,"DYTools_study2D",-1); // -1 : don't warn
  //if (!initOk) std::cout << " .. called from checkBinningArrays\n";
  TVectorD* mass= (TVectorD*)fin.FindObjectAny("massBinLimits");
  TVectorD* rapidityCounts= (TVectorD*)fin.FindObjectAny("rapidityCounts");
  if (!mass || !rapidityCounts) {
    std::cout << fncname << "file <" << fin.GetName() << "> does not contain keys massBinning and/or rapidityCounts\n";
    return false;
  }
  int comparisonOk=1;
  int noInitOnFile=0;
  if (initOk && ((*initOk)[0]!=DYTools::study2D)) {
    std::cout << " file recorded study2D=" << (*initOk)[0]
	      << ", current DYTools::study2D=" << DYTools::study2D << "\n";
    noInitOnFile=((*initOk)[0]==-1) ? 1:0;
    comparisonOk=0;
    delete initOk;
  }
  if (noInitOnFile) {
    std::cout << "binning ranges on file are not checked\n"
	      << "( file was produced w/o DYTools::setup() )\n";
  }
  else {
    int comparisonOk_ranges=checkBinningRanges(*mass,*rapidityCounts,
					       fin.GetName());
    if (comparisonOk) comparisonOk=comparisonOk_ranges;
  }
  if (mass) delete mass;
  if (rapidityCounts) delete rapidityCounts;
  return comparisonOk;
}

//--------------------------------------------------

int checkBinningRanges(const TVectorD &mass, const TVectorD &rapidityCounts, const TString &fname) {
  const char *fncname="checkBinningRanges: ";
  TString fileInfo=TString("on file <") + fname + ">";

  bool massOk=true, rapidityOk=true;
  int first_err=1;
  if (mass.GetNoElements() != DYTools::nMassBins+1) {
    std::cout << errdash; first_err=0;
    std::cout << "\n" << fncname
	      << " number of mass bins " << fileInfo
	      << " is " << mass.GetNoElements()
	      << " while " << (DYTools::nMassBins+1) << " is expected\n";
    massOk=false;
  }

  if (rapidityCounts.GetNoElements() != DYTools::nMassBins ) {
    if (first_err) { first_err=0; std::cout << errdash; }
    std::cout << "\n" << fncname
	      << "number of mass bins in rapidityCounts " << fileInfo
		<< " is " << rapidityCounts.GetNoElements()
	      << " while " << DYTools::nMassBins << " is expected\n";
    rapidityOk=false;
  }

  if (massOk) {
    for(int i=0; i<DYTools::nMassBins+1; i++){
      if( DYTools::massBinLimits[i] != mass[i] ) {
	if (first_err) { first_err=0; std::cout << errdash; }
	std::cout << fncname
		  << " mass limit " << fileInfo
		  << " at i=" << i
		  << " is " << mass[i]
		  << " instead of expected " << DYTools::massBinLimits[i] << "\n";
	massOk=false;
      }
    }
  }
  if (rapidityOk) {
    for (int i=0; i<DYTools::nMassBins; i++) {
      if ( DYTools::nYBins[i] != rapidityCounts[i] ) {
	if (first_err) { first_err=0; std::cout << errdash; }
	std::cout << fncname
		  << "y bin count " << fileInfo
		  << " at i=" << i
		  << " is " << rapidityCounts[i]
		  << " instead of expected " << DYTools::nYBins[i] << "\n";
	rapidityOk=false;
      }
    }
  }
  if (!massOk || !rapidityOk) {
    if (first_err) { first_err=0; std::cout << errdash; }
    std::cout << "file info: mass[" << mass.GetNoElements() << "]: ";
    for (int i=0; i<mass.GetNoElements(); ++i) {
      std::cout << " " << mass[i];
    }
    std::cout << "\n";
    std::cout << "file info: rapidityCounts["
	      << rapidityCounts.GetNoElements() << "]: ";
    for (int i=0; i<rapidityCounts.GetNoElements(); ++i) {
      std::cout << " " << rapidityCounts[i];
    }
    std::cout << "\n" << std::endl;
  }
  if (!first_err) std::cout << errdash;
  return (massOk && rapidityOk) ? 1:0;
}

//--------------------------------------------------

int checkMatrixSize(const TMatrixD &m, const TString &name) {
  int res= ((m.GetNrows() == DYTools::nMassBins) &&
	    (m.GetNcols() == DYTools::nYBinsMax)) ? 1:0;
  if (!res) {
    std::cout << "matrix (name=" << name << ") failed range check\n";
    std::cout << "expected: " << DYTools::nMassBins << "x" << DYTools::nYBinsMax << ", got " << m.GetNrows() << "x" << m.GetNcols() << std::endl;
  }
  return res;
}

//--------------------------------------------------
//--------------------------------------------------

void printMatrix(const TString &name, const TMatrixD &M, int exponent) {
  std::cout << "Matrix " << name << "\n";
  if (!exponent) M.Print();
  else {
    const char *format="%8.4e";
    for (int ir=0; ir<M.GetNrows(); ++ir) {
      for (int ic=0; ic<M.GetNcols(); ++ic) {
	std::cout << " " << Form(format,M(ir,ic));
      }
      std::cout << "\n";
    }
  }
}

//--------------------------------------------------

TMatrixD* loadMatrix(const TString &fname, const TString &fieldName, int expect_nRows, int expect_nCols, int reportFieldError) {
  TFile f(fname,"read");
  TMatrixD *M=NULL;
  int ok=1;
  if (!f.IsOpen()) ok=0;
  if (ok==1) {
    M=(TMatrixD*)f.Get(fieldName);
    f.Close();
    if (!M) ok=-1;
    else {
      if ((M->GetNrows()!=expect_nRows) ||
	  (M->GetNcols()!=expect_nCols)) {
	ok=-2;
      }
    }
  }
  if (ok!=1) {
    int report=1;
    if ((ok==-1) && !reportFieldError) report=0;
    if (report) {
      std::cout << "Error in loadMatrix(fname=<" << fname << ">, fieldName=<" << fieldName << ">, nRows=" << expect_nRows << ", nCols=" << expect_nCols << "):\n";
      if (ok==0) std::cout << " - failed to open the file\n";
      else if (ok==-1) std::cout << " - failed to load the field\n";
      else if (ok==-2) {
	std::cout << " - size mistmatch. Expect " << expect_nRows << "x" << expect_nCols << ", got " << M->GetNrows() << "x" << M->GetNcols() << "\n";
	delete M;
	M=NULL;
      }
    }
  }
  return M;
}

//--------------------------------------------------

TH2D* LoadMatrixFields(TFile &fin, const TString &field, const TString &fieldErr, int loadErr, int absoluteRapidity) {
  TMatrixD* val= NULL;
  if (field.Length()) val=(TMatrixD*)fin.Get(field);
  TMatrixD* err=(loadErr) ? (TMatrixD*)fin.Get(fieldErr) : NULL;
  if (   (field.Length() && !val)
	 || (loadErr && !err)
	 || (val && !checkMatrixSize(*val,field))
	 || (err && !checkMatrixSize(*err,fieldErr))
      //|| ((loadErr) && err && !checkMatrixSize(*err,fieldErr))
      ) {
    std::cout << "error in LoadMatrixFields(fname=" << fin.GetName() << ",""" << field << """, """ << fieldErr << """, loadErr=" << loadErr << "):\n";
    if (!val) std::cout << "\t failed to get values\n";
    if ((loadErr==1) && !err) std::cout << "\t failed to get errors\n";
    return NULL;
  }
  //std::cout << "loaded field=<" << field << "> from file <" << fin.GetName() << ">\n";
  //std::cout << "values are "; val->Print();
  TH2D *h2=createBaseH2(field,field,absoluteRapidity);
  for (int ibin=1; ibin<=h2->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=h2->GetNbinsY(); ++jbin) {
      double value=(val) ? (*val)(ibin-1,jbin-1) : 0.;
      double e=(loadErr==1) ? (*err)(ibin-1,jbin-1) : 0.;
      if (loadErr==2) e= ( (*err)(ibin-1,jbin-1) * (*err)(ibin-1,jbin-1) );
      if (loadErr==3) { value=0; e=(*err)(ibin-1,jbin-1); }
      //std::cout << "putting ibin=" << ibin << ", jbin=" << jbin << ", value=" << value << "\n";
      h2->SetBinContent(ibin,jbin, value);
      h2->SetBinError(ibin,jbin, e);
    }
  }
  return h2;
}

//--------------------------------------------------

TH2D* LoadMatrixFields(const TString &fname, int checkBinning, const TString &field, const TString &fieldErr, int loadErr, int absoluteRapidity) {
  TFile fin(fname,"read");
  if (!fin.IsOpen()) {
    std::cout << "LoadMatrixField(fname): failed to open file <" << fname << ">\n";
    return NULL;
  }
  int res=1;
  TH2D *h2=NULL;
  if (res && checkBinning) res=checkBinningArrays(fin);
  if (res) h2=LoadMatrixFields(fin,field,fieldErr,loadErr,absoluteRapidity);
  fin.Close();
  if (!res || !h2) std::cout << "error in LoadMatrixField(fname)\n";
  return h2;
}

//--------------------------------------------------

int LoadThreeMatrices(TFile &fin, TH2D **h2, TH2D **h2syst, const TString &field, const TString &fieldErr, const TString &fieldSystErr, int absoluteRapidity) {
  int res=1;
  (*h2)=LoadMatrixFields(fin,field,fieldErr,1,absoluteRapidity);
  if (!(*h2)) res=0;
  if (res) (*h2syst)=LoadMatrixFields(fin,"",fieldSystErr,3,absoluteRapidity);
  if (!(*h2syst)) res=0;
  if (!res) {
    std::cout << "error in LoadThreeMatrices(fname=" << fin.GetName() << ", """ << field << """, """ << fieldErr << """, """ << fieldSystErr << """)\n";
  }
  return res;
}

//--------------------------------------------------

int LoadThreeMatrices(const TString &fileName, TH2D **h2, TH2D **h2syst, const TString &field, const TString &fieldErr, const TString &fieldSystErr, int checkBinning, int absoluteRapidity) {
  TFile file(fileName,"read");
  if (!file.IsOpen()) {
    std::cout << "failed to open <" << fileName << ">\n";
    return 0;
  }
  int res=1;
  if (res && checkBinning) res=checkBinningArrays(file);
  if (res) {
    res=LoadThreeMatrices(file,h2,h2syst,
			  field,fieldErr,fieldSystErr,
			  absoluteRapidity);
  }
  file.Close();
  if (!res) {
    std::cout << "failed to load data from <" << fileName << ">\n";
  }
  return res;
}

//--------------------------------------------------
//--------------------------------------------------


inline
TH2D* extractSubArea(const TH2D *histo,
		     int xbin1, int xbin2, int ybin1, int ybin2,
		     const TString &newName, int setTitle, int resetAxis) {
  int nXBins=histo->GetNbinsX();
  int nYBins=histo->GetNbinsY();
  if ((xbin1==0) || (ybin1==0) ||
      (xbin1>xbin2) || (ybin1>ybin2)) {
    std::cout << "extractSubArea: problem with the requested area ";
    printf("(xbin1=%d, xbin2=%d, ybin1=%d, ybin2=%d)\n",xbin1,xbin2,ybin1,ybin2);
    std::cout << std::endl;
    return NULL;
  }
  if ((xbin2>nXBins) || (ybin2>nYBins)) {
    std::cout << "extractSubArea: \n";
    printf("requested last bin x=%d (available %d), y=%d (%d)",xbin2,nXBins,ybin2,nYBins);
    std::cout << std::endl;
    return NULL;
  }

  const TAxis *ax=histo->GetXaxis();
  const TAxis *ay=histo->GetYaxis();

  double *xnew=new double[xbin2-xbin1+2];
  double *ynew=new double[ybin2-ybin1+2];
  if (resetAxis) {
    // ignore the actual shift of the axes
    int lastXbin=xbin2-xbin1+1;
    int lastYbin=ybin2-ybin1+1;
    for (int i=0; i<lastXbin; ++i) xnew[i]=ax->GetBinLowEdge(i+1);
    xnew[lastXbin]= ax->GetBinLowEdge(lastXbin) + ax->GetBinWidth(lastXbin);
    for (int i=0; i<lastYbin; ++i) ynew[i]=ay->GetBinLowEdge(i+1);
    ynew[lastYbin]= ay->GetBinLowEdge(lastYbin) + ay->GetBinWidth(lastYbin);
  }
  else {
    // keep original labels of the axis
    for (int i=0; i<xbin2-xbin1+1; ++i) xnew[i]= ax->GetBinLowEdge(i+xbin1);
    xnew[xbin2-xbin1+1]= ax->GetBinLowEdge(xbin2) + ax->GetBinWidth(xbin2);
    for (int i=0; i<ybin2-ybin1+1; ++i) ynew[i]= ay->GetBinLowEdge(i+ybin1);
    ynew[ybin2-ybin1+1]= ay->GetBinLowEdge(ybin2) + ay->GetBinWidth(ybin2);
  }

  TH2D *h2=new TH2D(newName,"",xbin2-xbin1+1,xnew,ybin2-ybin1+1,ynew);
  if (setTitle) h2->SetTitle(newName);
  h2->SetDirectory(0);
  h2->SetStats(0);
  h2->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
  h2->GetYaxis()->SetTitle(histo->GetYaxis()->GetTitle());

  for (int ibin=xbin1; ibin<=xbin2; ++ibin) {
    for (int jbin=ybin1; jbin<=ybin2; ++jbin) {
      const int new_ibin= ibin-xbin1+1;
      const int new_jbin= jbin-ybin1+1;
      h2->SetBinContent(new_ibin,new_jbin,histo->GetBinContent(ibin,jbin));
      h2->SetBinError  (new_ibin,new_jbin,histo->GetBinError  (ibin,jbin));
    }
  }
  delete [] xnew;
  delete [] ynew;
  return h2;
}


//--------------------------------------------------
//--------------------------------------------------

void prepare(int count,
	     std::vector<TString> &pathV,
	     std::vector<TString> &fnameV,
	     std::vector<TString> &fieldV,
	     std::vector<TString> &labelV,
	     int clear,
	     int addEmptyElements) {
  if (clear) {
    pathV.clear();
    fnameV.clear();
    fieldV.clear();
    labelV.clear();
  }
  else count+=int(pathV.size());

  pathV.reserve(count);
  fnameV.reserve(count);
  fieldV.reserve(count);
  labelV.reserve(count);

  if (addEmptyElements) {
    for (int i=0; i<count; ++i) {
      TString empty=Form("empty_%d",i);
      pathV.push_back(empty);
      fnameV.push_back(empty);
      fieldV.push_back(empty);
      labelV.push_back(empty);
    }
  }
}

//--------------------------------------------------

void prepare(int count,
	     std::vector<TString> &pathV,
	     std::vector<TString> &fieldV,
	     std::vector<TString> &labelV,
	     int clear,
	     int addEmptyElements) {
  if (clear) {
    pathV.clear();
    fieldV.clear();
    labelV.clear();
  }
  else count+=int(pathV.size());

  pathV.reserve(count);
  fieldV.reserve(count);
  labelV.reserve(count);

  if (addEmptyElements) {
    for (int i=0; i<count; ++i) {
      TString empty=Form("empty_%d",i);
      pathV.push_back(empty);
      fieldV.push_back(empty);
      labelV.push_back(empty);
    }
  }
}

//--------------------------------------------------
//--------------------------------------------------

int saveLatexTable(TString fileTag,
		   const std::vector<TH2D*> &histosV,
		   const std::vector<TString> &labelsV_inp,
		   const char *format,
		   int printErrors,
		   int savePlainFormat) {
  TString fileName=Form("table_%dD_%s.tex",DYTools::study2D+1,
			fileTag.Data());
  std::ofstream fout(fileName);
  if (!fout.is_open()) {
    std::cout << "failed to create a file <" << fileName << ">\n";
    return 0;
  }

  int count = int(histosV.size());

  std::vector<TString> labelsV,subLabelsV;
  int hasSubLabels=0;
  labelsV.reserve(labelsV_inp.size());
  subLabelsV.reserve(labelsV_inp.size());
  for (unsigned int i=0; i<labelsV_inp.size(); i++) {
    TString txt=labelsV_inp[i];
    int p=txt.Index("//");
    if (p!=-1) {
      hasSubLabels=1;
      labelsV.push_back(txt(0,p));
      subLabelsV.push_back(txt(p+2,txt.Length()));
    }
    else {
      labelsV.push_back(txt);
      subLabelsV.push_back("");
    }
  }


  fout << "%\\documentclass{article}\n";
  fout << "%\\begin{document}\n";
  fout << "%\\usepackage{rotating}\n";
  fout << "%\\begin{document}\n";
  fout << "%\\begin{sidewaystable}[tbhp]\n";
  fout << "\\begin{table}[tbhp]\n";
  fout << "\\caption{\\label{tbl-" << fileTag << "} " << fileTag << "}\n";

  fout << "%\\hspace{-3cm}\n";
  fout << "\\begin{tabular}{|";
  for (int i=0; i<count+1+DYTools::study2D; ++i) fout << "c|";
  fout << "}\n";
  fout << "\\hline\n";

  if (DYTools::study2D==0) fout << " $m$ &";
  else fout << " $m$ & rapidity &";

  for (int i=0; i<count; ++i) {
    fout << " " << labelsV[i];
    if (i!=count-1) fout << " &";
  }
  fout << "\\\\\n";
  if (DYTools::study2D==0) fout << " (GeV) &";
  else fout << " (GeV) &       &";
  for (int i=0; i<count; ++i) {
    fout << " " << subLabelsV[i];
    if (i!=count-1) fout << " &";
  }
  fout << "\\\\\n";
  fout << "\\hline\n";

  TMatrixD *yBinLimits= DYTools::getYBinLimits();

  int ibinMin=(DYTools::study2D==0) ? 1:2;
  for (int ibin=ibinMin; ibin<=DYTools::nMassBins; ++ibin) {
    int yRangeCount= int(histosV[0]->GetNbinsY());
    if (DYTools::study2D && (ibin==DYTools::nMassBins)) yRangeCount=12;
    for (int jbin=1; jbin<=yRangeCount; ++jbin) {

      fout << " " << DYTools::massBinLimits[ibin-1] << " - "
	   << DYTools::massBinLimits[ibin] << " & ";
      if (DYTools::study2D==1) {
	fout << " " << Form("%3.1lf",(*yBinLimits)(ibin-1,jbin-1)) << " - "
	     << Form("%3.1lf",(*yBinLimits)(ibin-1,jbin)) << " & ";
      }

      for (int ih=0; ih<count; ++ih) {
	if (printErrors==0) { // only central values
	  fout << " " << Form(format,histosV[ih]->GetBinContent(ibin,jbin));
	  if (ih!=count-1) fout << " &";
	}
	else if (printErrors==1) { // central value and the error
	  fout << " " << Form(format,histosV[ih]->GetBinContent(ibin,jbin));
	  fout << " \\pm ";
	  fout << " " << Form(format,histosV[ih]->GetBinError(ibin,jbin));
	  if (ih!=count-1) fout << " &";
	}
	else if (printErrors==2) { // only errors
	  fout << " " << Form(format,histosV[ih]->GetBinError(ibin,jbin));
	  if (ih!=count-1) fout << " &";
	}
      }
      fout << "\\\\\\hline\n";
    }
  }

  fout << "\\hline\n";
  fout << "\\end{tabular}\n";
  fout << "\\end{table}\n";
  fout << "%\\end{sidewaystable}\n";
  fout << "%\\end{document}\n";
  fout.close();

  delete yBinLimits;
  std::cout << " saved file " << fileName << "\n";

  if (savePlainFormat) {
    TString fname=fileName;
    fname.ReplaceAll(".tex",".txt");
    std::ofstream foutTxt(fname);
    if (!foutTxt.is_open()) {
      std::cout << "failed to create text file <" << fname << ">\n";
      return 0;
    }
    foutTxt << count  << " " << histosV[0]->GetNbinsX() << " "
	    << histosV[0]->GetNbinsY()  << "\n";
    for (unsigned int i=0; i<histosV.size(); ++i) {
      foutTxt << "<<" << labelsV[i] << ">>\n";
      const TH2D* h=histosV[i];
      for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin) {
	for (int jbin=1; jbin<=h->GetNbinsY(); ++jbin) {
	  double x=h->GetBinContent(ibin,jbin);
	  if (x!=x) foutTxt << " " << double(0);
	  else foutTxt << " " << x;
	}
	foutTxt << "\n";
      }
    }
    foutTxt.close();
  }
  return 1;
}

//--------------------------------------------------

int loadLatexTableTextFile(TString fileTag,
			   std::vector<TH2D*> &histosV,
			   std::vector<TString> &labelsV,
			   int printErrors) {
  if (printErrors) {
    std::cout << "loadLatexTableTextFile is not ready for printError=" << printErrors << "\n";
    return 0;
  }

  TString fileName=Form("table_%dD_%s.txt",DYTools::study2D+1,
			fileTag.Data());
  std::ifstream fin(fileName);
  if (!fin.is_open()) {
    std::cout << "failed to open a file <" << fileName << ">\n";
    return 0;
  }

  int count=0, xsize=0, ysize=0;
  fin >> count >> xsize >> ysize;
  std::cout << "count=" << count << ", xsize=" << xsize << ", ysize=" << ysize << "\n";


  labelsV.reserve(count);
  for (int i=0; i<count; i++) labelsV.push_back(niceNumber(i+1,count+1));
  createBaseH2Vec(histosV,"histo_",labelsV,1,1);
  labelsV.clear();

  std::string s;
  for (int i=0; i<count; ++i) {
    getline(fin,s);
    if (s.size()==0) getline(fin,s);
    std::cout << "s.size=" << s.size() << "\n";
    size_t pos1=s.find_last_of('<');
    size_t pos2=s.find_first_of('>');
    std::string label=s.substr(pos1+1,pos2-pos1-1);
    std::cout << "label=<" << label << ">\n";
    HERE("aa");
    labelsV.push_back(label.c_str());
    TH2D* h=histosV[i];
    double val;
    for (int ibin=1; ibin<=xsize; ++ibin) {
      for (int jbin=1; jbin<=ysize; ++jbin) {
	fin >> val;
	h->SetBinContent(ibin,jbin, val);
      }
    }
  }
  fin.close();
  return 1;
}

//--------------------------------------------------
//--------------------------------------------------

TCanvas *compareHistos(TH2D* hA, TString labelA,
		       TH2D* hB, TString labelB,
		       int removeError, TString yAxisLabel) {
  std::vector<TH2D*> hV;
  std::vector<TString> labelV;
  hV.push_back(hA);
  hV.push_back(hB);
  labelV.push_back(labelA);
  labelV.push_back(labelB);
  TCanvas *cx= plotProfiles("cx",hV,labelV,NULL,
			    removeError,yAxisLabel);
  return cx;
}

//--------------------------------------------------

TCanvas *plotProfiles(TString canvName, TH2D* histo, TString label,
		      int do_removeError, TString yAxisLabel) {
  if (!histo) {
    std::cout << "plotProfiles: null histo" << std::endl;
    return NULL;
  }
  if (DYTools::study2D==-1) {
    std::cout << "plotProfiles: setup the analysis. "
	      << "Call DYTools::setup(DIM-1)\n";
    return NULL;
  }
  std::vector<TH2D*> hV;
  std::vector<TString> labelV;
  hV.push_back(histo);
  labelV.push_back(label);
  TCanvas *cx= plotProfiles(canvName,hV,labelV,NULL,
			    do_removeError,yAxisLabel);
  return cx;
}

//--------------------------------------------------
//--------------------------------------------------

TCanvas* plotProfiles(TString canvName,
		      const std::vector<TH2D*> &histosV,
		      const std::vector<TString> &labelsV,
		      std::vector<int> *colorsV,
		      int do_removeError,
		      TString yAxisLabel,
		      std::vector<std::vector<TH1D*>*> *hProfV,
		      std::vector<ComparisonPlot_t*> *cpV,
		      int delayDraw) {
  if (!hProfV) hProfV= new std::vector<std::vector<TH1D*>*>();

  const unsigned int colorCount=6;
  const int autoColor[colorCount] = { kBlack, kBlue, kGreen+1, kOrange+1, kRed+1, kViolet };
  int ourColors=0;
  if (colorsV==NULL) {
    ourColors=1;
    colorsV=new std::vector<int>();
    colorsV->reserve(histosV.size());
    for (unsigned int i=0; i<histosV.size(); ++i) {
      colorsV->push_back(autoColor[i%colorCount]);
    }
  }

  int oneHisto=(histosV.size()==1) ? 1 : 0;
  int canvWidth=(DYTools::study2D==1) ? 1200 : ((oneHisto) ? 800 : 700);
  int canvHeight=800;
  TCanvas *c1=new TCanvas(canvName,canvName, canvWidth,canvHeight);

  if (DYTools::study2D==1) {

    if (!createRapidityProfileVec(histosV,*hProfV,labelsV)) {
      std::cout << "failed to create profiles\n";
      return NULL;
    }
    std::cout << "there are " << hProfV->size() << " profiles\n";

    for (int im=1; im<7; ++im) {

      TString mStr=Form("M_%2.0lf_%2.0lf",DYTools::massBinLimits[im],DYTools::massBinLimits[im+1]);
      TString cpName="cp_" + mStr;
      TString cpTitle=mStr;
      ComparisonPlot_t *cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,cpName,cpTitle,"|y|",yAxisLabel,"ratio");
      if (cpV) cpV->push_back(cp);
      if (oneHisto) {
	cp->SetRefIdx(-111); // no ratio
	cp->SetNoLegend(1);
      }
      if (im==1) {
	if (!oneHisto) cp->Prepare6Pads(c1,1);
	else c1->Divide(3,2);
      }

      for (unsigned int ih=0; ih<(*hProfV)[im]->size(); ++ih) {
	TH1D* h=(*(*hProfV)[im])[ih];
	if (do_removeError) removeError1D(h);
	if (ih==0) {
	  h->SetMarkerStyle(20);
	}
	if ( ourColors && (ih/colorCount > 0) ) {
	  //std::cout << "changing the marker\n";
	  h->SetMarkerStyle(5);
	}
	cp->AddHist1D(h,labelsV[ih],"LPE1",(*colorsV)[ih],1,0,1);
      }
      if (!delayDraw) {
	if (!oneHisto) cp->Draw6(c1,1,im);
	else { cp->Draw(c1,false,"png",im); }
      }
    }
  }
  else {
    // 1D

    if (!createMassProfileVec(histosV,*hProfV,labelsV)) {
      std::cout << "failed to create profiles\n";
      return NULL;
    }
    std::cout << "there are " << hProfV->size() << " profiles\n";

    for (int iy=0; iy<DYTools::nYBinsMax; ++iy) {

      TString yStr=Form("iy_%d",iy);
      TString cpName=TString("cp_") + yStr;
      TString cpTitle; //=yStr;
      ComparisonPlot_t *cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,cpName,cpTitle,"#it{M}_{ee} [GeV]",yAxisLabel,"ratio");
      if (cpV) cpV->push_back(cp);
      if (oneHisto) {
	cp->SetRefIdx(-111); // no ratio
	cp->SetNoLegend(1);
      }
      cp->SetLogx(1);
      if ((iy==0) && !oneHisto) cp->Prepare2Pads(c1);

      for (unsigned int ih=0; ih<(*hProfV)[iy]->size(); ++ih) {
	TH1D* h=(*(*hProfV)[iy])[ih];
	if (oneHisto && (ih==0)) {
	  h->GetXaxis()->SetMoreLogLabels();
	  h->GetXaxis()->SetNoExponent();
	}
	if (do_removeError) removeError1D(h);
	if (ih==0) {
	  h->SetMarkerStyle(20);
	}
	if ( ourColors && (ih/colorCount > 0) ) {
	  //std::cout << "changing the marker\n";
	  h->SetMarkerStyle(5);
	}
	cp->AddHist1D(h,labelsV[ih],"LPE1",(*colorsV)[ih],(ih+1)%3,0,1);
      }
      if (!delayDraw) {
	if (!oneHisto) cp->Draw(c1);
	else cp->Draw(c1,false,"png",0);
      }

    }
  }
  c1->Update();
  if (ourColors) delete colorsV;
  return c1;
}

//--------------------------------------------------
//--------------------------------------------------

