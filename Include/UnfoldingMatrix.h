#ifndef UnfoldingMatrix_H
#define UnfoldingMatrix_H

#include <TFile.h>
#include <TRandom.h>
#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/MitStyleRemix.hh"
#include "../Include/CPlot.hh"
#include "../Include/FlatIndex.h"

// a pair of histograms
#include "../Include/HistoPair.hh"

//for getting matrix condition number
#include <TDecompLU.h>


//=== FUNCTION DECLARATIONS ======================================================================================

/*
int nUnfoldingBins = DYTools::getTotalNumberOfBins();

void computeNormalizedBinContent(double subset, double subsetErr,
				 double total, double totalErr,
				 double& ratio, double& ratioErr);
void calculateInvertedMatrixErrors(const TMatrixD &T,
          const TMatrixD &TErrPos, const TMatrixD &TErrNeg,
				   TMatrixD &TinvErr, int NsmearingExps);
*/

/*
inline 
int validFlatIndices(const FlatIndex_t &fi_a, const FlatIndex_t &fi_b) {
  return (DYTools::validFlatIndex(fi_a.idx()) && DYTools::validFlatIndex(fi_b.idx()));
}
*/

//=== FUNCTION DEFINITIONS ======================================================================================

inline
int computeNormalizedBinContent(double subset, double subsetErr,
				 double total, double totalErr,
				 double& ratio, double& ratioErr){

  ratio=0.; ratioErr=0.;
  if(total == 0) {
    printf("makeUnfoldingMatrix::Possible problem\n");
    printf("     empty column in the response matrix\n");
    return 0;
  }
  
  ratio = subset/total;

  // The formula for the ratio = subset/total is obtained by error
  // propagation. The subsetErr and totalErr are NOT assumed ot be
  // the sqrt(subset) and sqrt(total). (If one assume that, the formula
  // below reduces to the familiar sqrt(ratio*(1-ratio)/total) ).
  // The subset and subsetErr are part of the total and totalErr.
  // The formula is easiest to derive if we take "A +- dA" and
  // "B +- dB" as independent numbers, with total = A+B and
  // totalErr^2 = dA^2 + dB^2. One then does error propagation of
  // the expression ratio = A/(A+B).
  // The outcome of it is found below (the absolute error on the ratio)
  ratioErr = (1/total)*sqrt( subsetErr*subsetErr*(1-2*ratio)
			     + totalErr*totalErr*ratio*ratio );

  return 1;
}

// -----------------------------------------------

inline
void calculateInvertedMatrixErrors(const TMatrixD &T, 
	      const TMatrixD &TErrPos, const TMatrixD &TErrNeg,
		    TMatrixD &TinvErr, int Nsmearings=1000){

  // Calculate errors on the inverted matrix by the Monte Carlo
  // method

  Double_t det;
  int nRow = T.GetNrows();
  int nCol = T.GetNcols();
  TMatrixD TinvSum(nRow,nCol);
  TMatrixD TinvSumSquares(nRow,nCol);

  // Reset Matrix where we will be accumulating RMS/sigma:
  TinvSum        = 0;
  TinvSumSquares = 0;
  TinvErr        = 0;

  // Do many tries, accumulate RMS
  //int N = 10000;
  //N=1000; std::cout << "\ncalculateInvertedMatrixErrors: setting Ntries=" << N << "\n\n";
  int N=Nsmearings;
  std::cout << "\tcalculateInvertedMatrixErrors: Nsmearings=" << Nsmearings
	    << "\n";
  for(int iTry = 0; iTry<N; iTry++){
    // Find the smeared matrix
    TMatrixD Tsmeared = T;
    for(int i = 0; i<nRow; i++){
      for(int j = 0; j<nCol; j++){
	double central = T(i,j);
	double sigPos = TErrPos(i,j);
	double sigNeg = TErrNeg(i,j);
 	// Switch to symmetric errors: approximation, but much simpler
	double sig = (sigPos+sigNeg)/2.0;
	Tsmeared(i,j) = gRandom->Gaus(central,sig);
      }
    }
    // Find the inverted to smeared matrix
    TMatrixD TinvSmeared = Tsmeared;
    TinvSmeared.Invert(&det);
    // Accumulate sum and sum of squares for each element
    for(int i2 = 0; i2<nRow; i2++){
      for(int j2 = 0; j2<nCol; j2++){
	TinvSum       (i2,j2) += TinvSmeared(i2,j2);
	TinvSumSquares(i2,j2) += TinvSmeared(i2,j2)*TinvSmeared(i2,j2);
      }
    }
  }

  // Calculate the error matrix
  //TMatrixD TinvAverage = TinvSum; // average values
  for(int i = 0; i<nRow; i++){
    for(int j = 0; j<nCol; j++){
      TinvErr(i,j) = sqrt( TinvSumSquares(i,j)/double(N) 
			   - (TinvSum(i,j)/double(N))*(TinvSum(i,j)/double(N)) );
      //TinvAverage(i,j) = TinvSum(i,j)/double(N);
    }
  }

  return;
}

// -----------------------------------------------

// ------------------------------------------------------------

inline
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

inline
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

template<class Container_t>
double maxDiff(const Container_t &A, const Container_t &B, int *storeRow=NULL, int *storeCol=NULL, int printDiff=0) {
  Container_t diff=A; diff-=B;
  if (printDiff) { std::cout << "diff="; diff.Print(); }
  double max= maxElement(diff,storeRow,storeCol);
  return max;
}
  
// ------------------------------------------------------------

inline
void testMaxDiff(const TString &msg, const TVectorD &a, const TVectorD &b, int printDiff=0) {
  int i=-1;
  double md=maxDiff(a,b,&i,NULL,printDiff);
  std::cout << msg << ". MaxDiff=" << md << "\n";
  std::cout << " values there (at " << i << "): " << a[i] << " and " << b[i] << "\n";
  return;
}

// ------------------------------------------------------------

inline
void testMaxDiff(const TString &msg, const TMatrixD &a, const TMatrixD &b, int printDiff=0) {
  int ir=-1, ic=-1;
  double md=maxDiff(a,b,&ir,&ic,printDiff);
  std::cout << msg << ". MaxDiff=" << md << "\n";
  std::cout << " values there, at (" << ir << "," << ic << "): " 
	    << a(ir,ic) << " and " << b(ir,ic) << "\n";
  return;
}

// ------------------------------------------------------------

inline
void testMaxRelDiff(const TString &msg, const TVectorD &a, const TVectorD &b) {
  int i=-1;
  TVectorD diff=a-b;
  TVectorD avg=a+b;  avg*=0.5;
  for (int ii=0; ii<diff.GetNoElements(); ++ii) diff(ii) /= avg[ii];
  double max= maxElement(diff,&i);
  std::cout << msg << ". MaxRelDiff=" << max << "\n";
  std::cout << " values there (at " << i << "): " << a[i] << " and " << b[i] << "\n";
  return;
}

// ------------------------------------------------------------

inline
void testMaxRelDiff(const TString &msg, const TMatrixD &a, const TMatrixD &b) {
  int ir=-1, ic=-1;
  TMatrixD diff=a-b;
  TMatrixD avg=a+b;  avg*=0.5;
  for (int i=0; i<diff.GetNrows(); ++i) {
    for (int j=0; j<diff.GetNcols(); ++j) {
      diff(i,j) /= avg(i,j);
    }
  }
  double max= maxElement(diff,&ir,&ic);
  std::cout << msg << ". MaxRelDiff=" << max << "\n";
  std::cout << " values there, at (" << ir << "," << ic << "): " 
	    << a(ir,ic) << " and " << b(ir,ic) << "\n";
  return;
}

// ------------------------------------------------------------


// ---------------------------------------------------------------------

//=== Typedef =================================================================================================

namespace UnfoldingMatrix {
  typedef enum { _cDET_Response, _cFSR, _cFSR_DET, _cFSR_DETcorrFactors } TUnfoldingMatrixType_t;

  inline
  TString kindName(TUnfoldingMatrixType_t aKind) {
    TString s="unknown";
    switch(aKind) {
    case _cDET_Response: s="DetResponse"; break;
    case _cFSR: s="FSR"; break;
    case _cFSR_DET: s="FSR_DET"; break;
    case _cFSR_DETcorrFactors: s="FSR_DETcorrFactors"; break;
    }
    return s;
  }


  inline
  void getYieldNames(TUnfoldingMatrixType_t theKind, TString &iniName, TString &finName) {
    switch(theKind) {
    case _cDET_Response: 
      iniName="yieldsMcPostFsrGen";
      finName="yieldsMcPostFsrRec";
      break;
    case _cFSR:
      iniName="yieldsMcPreFsrGen";
      finName="yieldsMcPostFsrGen";
      break;
    case _cFSR_DET:
      iniName="yieldsMcPreFsrGenDET";
      finName="yieldsMcPostFsrGenDET";
      break;
    case _cFSR_DETcorrFactors:
      iniName="fsrDETcorrFactorsGen";
      finName="fsrDETcorrFactorsReco";
      break;
    default:
      std::cout << "getYieldNames cannot handle this 'kind'=" << theKind << "\n";
      assert(0);
    }
  }

};

//=== Class DECLARATIONS =================================================================================================

class UnfoldingMatrix_t {
public:
  // DET_Response: ini --> PostFsrGen, fin --> PostFsrReco
  // FSR, FSR_DET: ini --> PreFsrGen,  fin -->PostFsrGen
  // FSR_DETcorrFactors: factors, correcting preFSR and postFSR acceptance
  // yieldsIni: gen, yieldsFin: reco
public:
  UnfoldingMatrix::TUnfoldingMatrixType_t kind;
  TString name, iniYieldsName, finYieldsName;
  TMatrixD *yieldsIni; //(DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD *yieldsFin; //(DYTools::nMassBins,DYTools::nYBinsMax);

  // Matrices for unfolding
  TMatrixD *DetMigration; //(DYTools::nUnfoldingBins, nUnfoldingBins);
  TMatrixD *DetMigrationErr; //(nUnfoldingBins, nUnfoldingBins);
  TMatrixD *DetResponse; //(nUnfoldingBins, nUnfoldingBins);
  TMatrixD *DetResponseErrPos; //(nUnfoldingBins, nUnfoldingBins);
  TMatrixD *DetResponseErrNeg; //(nUnfoldingBins, nUnfoldingBins);
  TMatrixD *DetInvertedResponse; //(nUnfoldingBins, nUnfoldingBins);
  TMatrixD *DetInvertedResponseErr; //(nUnfoldingBins, nUnfoldingBins);

  TVectorD *DetResponseArr; //(nUnfoldingBins);
  TVectorD *DetInvertedResponseArr; //(nUnfoldingBins);
  TVectorD *DetInvertedResponseErrArr; //(nUnfoldingBins);
  TVectorD *yieldsIniArr; //(nUnfoldingBins)
  TVectorD *yieldsFinArr; //(nUnfoldingBins);

public:
  TString ourKindName() const { return UnfoldingMatrix::kindName(this->kind); }

public:
  UnfoldingMatrix_t(UnfoldingMatrix::TUnfoldingMatrixType_t set_kind, const TString &set_name) : 
    kind(set_kind),
    name(set_name), iniYieldsName(), finYieldsName(),
    yieldsIni(0), yieldsFin(0),
    DetMigration(0), DetMigrationErr(0),
    DetResponse(0), DetResponseErrPos(0), DetResponseErrNeg(0),
    DetResponseArr(0), DetInvertedResponseArr(0), DetInvertedResponseErrArr(0),
    yieldsIniArr(0), yieldsFinArr(0)
  {
    TMatrixD my(DYTools::nMassBins,DYTools::nYBinsMax);
    TMatrixD unf(DYTools::nUnfoldingBins, DYTools::nUnfoldingBins);
    TVectorD arr(DYTools::nUnfoldingBins);
    my=0; unf=0; arr=0;
    yieldsIni=new TMatrixD(my); yieldsFin=new TMatrixD(my);
    DetMigration=new TMatrixD(unf); DetMigrationErr=new TMatrixD(unf);
    DetResponse=new TMatrixD(unf); 
    DetResponseErrPos=new TMatrixD(unf); DetResponseErrNeg=new TMatrixD(unf);
    DetInvertedResponse=new TMatrixD(unf); DetInvertedResponseErr=new TMatrixD(unf);
    DetResponseArr=new TVectorD(arr);
    DetInvertedResponseArr=new TVectorD(arr); DetInvertedResponseErrArr=new TVectorD(arr);
    yieldsIniArr=new TVectorD(arr); yieldsFinArr=new TVectorD(arr);
    UnfoldingMatrix_t::getYieldNames(set_kind, iniYieldsName, finYieldsName);
  }

  UnfoldingMatrix_t(const UnfoldingMatrix_t &U, const TString &set_name) :
    kind(U.kind),
    name(set_name),
    iniYieldsName(), finYieldsName(),
    yieldsIni(0), yieldsFin(0),
    DetMigration(0), DetMigrationErr(0),
    DetResponse(0), DetResponseErrPos(0), DetResponseErrNeg(0),
    DetResponseArr(0), DetInvertedResponseArr(0), DetInvertedResponseErrArr(0),
    yieldsIniArr(0), yieldsFinArr(0)
  {
    yieldsIni=new TMatrixD(*U.yieldsIni); 
    yieldsFin=new TMatrixD(*U.yieldsFin);
    DetMigration=new TMatrixD(*U.DetMigration); 
    DetMigrationErr=new TMatrixD(*U.DetMigrationErr);
    DetResponse=new TMatrixD(*U.DetResponse); 
    DetResponseErrPos=new TMatrixD(*U.DetResponseErrPos); 
    DetResponseErrNeg=new TMatrixD(*U.DetResponseErrNeg);
    DetInvertedResponse=new TMatrixD(*U.DetInvertedResponse); 
    DetInvertedResponseErr=new TMatrixD(*U.DetInvertedResponseErr);
    DetResponseArr=new TVectorD(*U.DetResponseArr);
    DetInvertedResponseArr=new TVectorD(*U.DetInvertedResponseArr);
    DetInvertedResponseErrArr=new TVectorD(*U.DetInvertedResponseErrArr);
    yieldsIniArr=new TVectorD(*U.yieldsIniArr); 
    yieldsFinArr=new TVectorD(*U.yieldsFinArr);
    UnfoldingMatrix_t::getYieldNames(kind, iniYieldsName, finYieldsName);
  }

  ~UnfoldingMatrix_t() {
    if (yieldsIni) delete yieldsIni;
    if (yieldsFin) delete yieldsFin;
    if (DetMigration) delete DetMigration;
    if (DetMigrationErr) delete DetMigrationErr;
    if (DetResponse) delete DetResponse;
    if (DetResponseErrPos) delete DetResponseErrPos;
    if (DetResponseErrNeg) delete DetResponseErrNeg;
    if (DetInvertedResponse) delete DetInvertedResponse;
    if (DetInvertedResponseErr) delete DetInvertedResponseErr;
    if (DetResponseArr) delete DetResponseArr;
    if (DetInvertedResponseArr) delete DetInvertedResponseArr;
    if (DetInvertedResponseErrArr) delete DetInvertedResponseErrArr;
    if (yieldsIniArr) delete yieldsIniArr;
    if (yieldsFinArr) delete yieldsFinArr;
  }

  const TString& getName() const { return name; }
  const TMatrixD* getDetResponse() const { return DetResponse; }
  const TMatrixD* getDetInvResponse() const { return DetInvertedResponse; }
  const TMatrixD* getIniM() const { return yieldsIni; }
  const TMatrixD* getFinM() const { return yieldsFin; }

  // Note that they are constructed from TMatrixD after the call of
  // prepareFIArrays
  const TVectorD* getIniVec() const { return yieldsIniArr; }
  const TVectorD* getFinVec() const { return yieldsFinArr; }
  // Read the note above.

  int sizesMatch(const UnfoldingMatrix_t &U) const {
    int ok=((yieldsIni->GetNrows() == U.yieldsIni->GetNrows()) &&
	    (yieldsIni->GetNcols() == U.yieldsIni->GetNcols()) &&
	    (yieldsIniArr->GetNoElements() == U.yieldsIniArr->GetNoElements())) ? 1:0;
    return ok;
  }

  void fillIni(int iMassBinGen, int iYBinGen, double fullWeight) {
    using namespace DYTools;
    if ((iMassBinGen==-1) || (iYBinGen==-1)) {
      return;
    }
    if ((iMassBinGen >= nMassBins) ||
	(iYBinGen >= nYBins[iMassBinGen])) {
      std::cout << "UnfoldingMatrix_t::fillIni(" << iMassBinGen << "," << iYBinGen << "): incorrect indices. Max values (" << nMassBins << "," << (((iMassBinGen>=0) && (iMassBinGen<nMassBins)) ? nYBins[iMassBinGen] : nYBinsMax) << "; matrixName=<" << name << ">\n";
      assert(0);
    }
    (*yieldsIni)(iMassBinGen,iYBinGen) += fullWeight;
  }

  void fillFin(int iMassBinReco, int iYBinReco, double fullWeight) {
    using namespace DYTools;
    if ((iMassBinReco==-1) || (iYBinReco==-1)) {
      return;
    }
    if ((iMassBinReco >= nMassBins) ||
	(iYBinReco >= nYBins[iMassBinReco])) {
      std::cout << "UnfoldingMatrix_t::fillPostFin(" << iMassBinReco << "," << iYBinReco << "): incorrect indices. Max values (" << nMassBins << "," << (((iMassBinReco>=0) && (iMassBinReco<nMassBins)) ? nYBins[iMassBinReco] : nYBinsMax) << "; matrixName=<" << name << ">\n";
      assert(0);
    }
    (*yieldsFin)(iMassBinReco,iYBinReco) += fullWeight;
  }

  void fillMigration(int idx1, int idx2, double weight) {
    if ( !DYTools::validFlatIndices(idx1,idx2) ) {
      std::cout << "UnfoldingMatrix_t::fillMigration: idx1=" << idx1 << ", idx2=" << idx2 << ". Max allowed values: " << DYTools::nUnfoldingBins << "; matrixName=<" << name << ">\n";
      return;
    }
    (*DetMigration)(idx1,idx2) += weight;
    (*DetMigrationErr)(idx1,idx2) += weight * weight;
  }

  void fillIni(const FlatIndex_t &fi, double weight) { return fillIni(fi.iM(),fi.iY(),weight); }
  void fillFin(const FlatIndex_t &fi, double weight) { return fillFin(fi.iM(),fi.iY(),weight); }
  void fillMigration(const FlatIndex_t &fi_ini, const FlatIndex_t &fi_fin, double weight) { fillMigration(fi_ini.idx(),fi_fin.idx(),weight); }


  void finalizeDetMigrationErr() {
    for(int i=0; i < (*DetMigration).GetNrows(); i++)
      for(int j=0; j < (*DetMigration).GetNcols(); j++)
	if( (*DetMigrationErr)(i,j) >=0 )
	  (*DetMigrationErr)(i,j) = sqrt( (*DetMigrationErr)(i,j) );
	else {
	  printf("UnfoldingMatrix_t::finalizeDetMigrationErr Error: negative weights in DetMigrationErr\n");
	  std::cout << "matrixName=<" << name << ">\n";
	  assert(0);
	}
  }

  void squareDetMigrationErr() {
    for (int i=0; i<(*DetMigration).GetNrows(); ++i) {
      for (int j=0; j<(*DetMigration).GetNcols(); ++j) {
	double x=(*DetMigrationErr)(i,j);
	(*DetMigrationErr)(i,j)=x*x;
      }
    }
  }

  int randomizeMigrationMatrix(const UnfoldingMatrix_t &U,
			       const UnfoldingMatrix_t *Uexact=NULL,
			       int NsmearingsForError=1000) {
    if (!this->sizesMatch(U)) {
      std::cout << "UnfoldingMatrix_t::randomizeMigrationMatrix(U) sizes do not match\n";
      return 0;
    }

    if (Uexact) {
      if (!this->sizesMatch(*Uexact)) {
	std::cout << "UnfoldingMatrix_t::randomizeMigrationMatrix(U,*Uexact) sizes do not match\n";
	return 0;
      }
    }

    for (int ir=0; ir<U.DetMigration->GetNrows(); ++ir) {
      for (int ic=0; ic<U.DetMigration->GetNcols(); ++ic) {
	(*this->DetMigration)(ir,ic) = (*U.DetMigration)(ir,ic) + 
	  (*U.DetMigrationErr)(ir,ic) * gRandom->Gaus(0,1);
      }
    }

    *DetMigrationErr= (*U.DetMigrationErr);
    this->computeResponseMatrix();
    this->invertResponseMatrix(NsmearingsForError);
    if (Uexact) {
      // this is a FSR matrix
      *yieldsIni= (*U.yieldsIni);
      *yieldsFin= (*U.yieldsFin);
      this->modifyDETResponseMatrices(*Uexact);
    }
    this->prepareFIArrays();
    return 1;
  }

  // errors should be squared!
  int setErrorOnMigrationMatrix(const UnfoldingMatrix_t &U1, const UnfoldingMatrix_t &U2) {
    if (!this->sizesMatch(U1) || !this->sizesMatch(U2)) {
      std::cout << "UnfoldingMatrix_t::setErrorOnMigrationMatrix sized do not match\n";
      return 0;
    }
    for (int ir=0; ir<this->DetMigration->GetNrows(); ++ir) {
      for (int ic=0; ic<this->DetMigration->GetNcols(); ++ic) {
	double ref=(*DetMigration)(ir,ic);
	double v1 =(*U1.DetMigration)(ir,ic);
	double v2 =(*U2.DetMigration)(ir,ic);
	double diff=fabs(v1-ref);
	if (fabs(v2-ref)>diff) diff=fabs(v2-ref);
	if (fabs(v1-v2 )>diff) diff=fabs(v1-v2);
	(*DetMigrationErr)(ir,ic) = diff;
      }
    }
    return 1;
  }

  // errors should be squared!
  int addMigration(const UnfoldingMatrix_t &U, double weight) {
    TMatrixD Ytmp= *U.yieldsIni;
    Ytmp *= weight;
    *this->yieldsIni += Ytmp;
    Ytmp= *U.yieldsFin;
    Ytmp *= weight;
    *this->yieldsFin += Ytmp;
    TMatrixD Mtmp= *U.DetMigration;
    Mtmp *= weight;
    *this->DetMigration += Mtmp;
    Mtmp= *U.DetMigrationErr;
    Mtmp *= (weight*weight);
    *this->DetMigrationErr += Mtmp;
    return 1;
  }

  void computeResponseMatrix() {
    double tCentral, tErr;
    for(int igen = 0; igen < (*DetMigration).GetNrows(); igen++){
      // First find the normalization for the given generator level slice
      double nEventsInGenBin = 0;
      double nEventsInGenBinErr = 0;
      for(int ireco = 0; ireco < (*DetMigration).GetNcols(); ireco++){
	nEventsInGenBin += (*DetMigration)(igen,ireco);
	nEventsInGenBinErr += ((*DetMigrationErr)(igen,ireco)*
			       (*DetMigrationErr)(igen,ireco));
      }
      nEventsInGenBinErr = sqrt(nEventsInGenBinErr);
      
      // Now normalize each element and find errors
      for(int ireco = 0; ireco < (*DetMigration).GetNcols(); ireco++){
	tCentral = 0;
	tErr     = 0;
	int res=computeNormalizedBinContent((*DetMigration)(igen,ireco),
					    (*DetMigrationErr)(igen,ireco),
					    nEventsInGenBin,
					    nEventsInGenBinErr,
					    tCentral, tErr);
	if (!res) std::cout << name << ": igen=" << igen << ", ireco=" << ireco << "\n";
	(*DetResponse)      (igen,ireco) = tCentral;
	(*DetResponseErrPos)(igen,ireco) = tErr;
	(*DetResponseErrNeg)(igen,ireco) = tErr;
      }
    }
  }

  /* // A.J. test on Sept 20, 2012. This method is not accurate 
  void computeResponseMatrix_MdfBeforeNormalization(const UnfoldingMatrix_t &exact) {
    //this->computeResponseMatrix();
    std::cout << "computeResponseMatrix_Mdf\n";
    double tCentral, tErr;
    TVectorD iniV(DYTools::nUnfoldingBins), finV(DYTools::nUnfoldingBins);
    TVectorD iniVexact(DYTools::nUnfoldingBins), finVexact(DYTools::nUnfoldingBins);
    int resFlatten= 
      (
       (flattenMatrix(*yieldsIni, iniV) == 1 ) &&
       (flattenMatrix(*yieldsFin, finV) == 1 ) &&
       (flattenMatrix(*exact.yieldsIni, iniVexact) == 1 ) &&
       (flattenMatrix(*exact.yieldsFin, finVexact) == 1 ) 
       ) ? 1:0;
    assert(resFlatten);

    for(int igen = 0; igen < (*DetMigration).GetNrows(); igen++){
      // First find the normalization for the given generator level slice
      double nEventsInGenBin = 0;
      double nEventsInGenBinErr = 0;
      for(int ireco = 0; ireco < (*DetMigration).GetNcols(); ireco++){
	nEventsInGenBin += (*DetMigration)(igen,ireco);
	nEventsInGenBinErr += ((*DetMigrationErr)(igen,ireco)*
			       (*DetMigrationErr)(igen,ireco));
      }
      nEventsInGenBinErr = sqrt(nEventsInGenBinErr);
      
      // rescale event number in the generated bin
      printf("igen=%d, nEventsInGenBin=%9.4lf, ini/fin=%9.4lf/%9.4lf, exact ini/fin=%9.4lf/%9.4lf\n",igen,nEventsInGenBin,iniV[igen],finV[igen],iniVexact[igen],finVexact[igen]);
      double scaleIni=
	( iniVexact[igen] == 0 ) ? 0. : (iniV[igen]/iniVexact[igen]);
      nEventsInGenBin *= scaleIni;
      nEventsInGenBinErr *= scaleIni;
      printf(" scaled number of events: %9.4lf\n",nEventsInGenBin);

      // Now normalize each element and find errors
      for(int ireco = 0; ireco < (*DetMigration).GetNcols(); ireco++){
      // prepare scale in the reconstructed bin
	double scaleFin=
	  ( finVexact[ireco] == 0 ) ? 0. : (finV[ireco]/finVexact[ireco] );

	tCentral = 0;
	tErr     = 0;
	computeNormalizedBinContent((*DetMigration)(igen,ireco) * scaleFin,
				    (*DetMigrationErr)(igen,ireco) * scaleFin,
				    nEventsInGenBin,
				    nEventsInGenBinErr,
				    tCentral, tErr);
	(*DetResponse)      (igen,ireco) = tCentral;
	(*DetResponseErrPos)(igen,ireco) = tErr;
	(*DetResponseErrNeg)(igen,ireco) = tErr;
      }
    }
  }
  */

  /*
  // This is Ilya's suggestion: modify the response matrix
  void computeResponseMatrix_Mdf(const UnfoldingMatrix_t &exact) {
    //this->computeResponseMatrix();
    std::cout << "computeResponseMatrix_Mdf\n";
    double tCentral, tErr;
    TVectorD iniV(DYTools::nUnfoldingBins), finV(DYTools::nUnfoldingBins);
    TVectorD iniVexact(DYTools::nUnfoldingBins), finVexact(DYTools::nUnfoldingBins);
    int resFlatten= 
      (
       (flattenMatrix(*yieldsIni, iniV) == 1 ) &&
       (flattenMatrix(*yieldsFin, finV) == 1 ) &&
       (flattenMatrix(*exact.yieldsIni, iniVexact) == 1 ) &&
       (flattenMatrix(*exact.yieldsFin, finVexact) == 1 ) 
       ) ? 1:0;
    assert(resFlatten);

    for(int igen = 0; igen < (*DetMigration).GetNrows(); igen++){
      // First find the normalization for the given generator level slice
      double nEventsInGenBin = 0;
      double nEventsInGenBinErr = 0;
      for(int ireco = 0; ireco < (*DetMigration).GetNcols(); ireco++){
	nEventsInGenBin += (*DetMigration)(igen,ireco);
	nEventsInGenBinErr += ((*DetMigrationErr)(igen,ireco)*
			       (*DetMigrationErr)(igen,ireco));
      }
      nEventsInGenBinErr = sqrt(nEventsInGenBinErr);
      
      //printf("igen=%d, nEventsInGenBin=%9.4lf, ini/fin=%9.4lf/%9.4lf, exact ini/fin=%9.4lf/%9.4lf\n",igen,nEventsInGenBin,iniV[igen],finV[igen],iniVexact[igen],finVexact[igen]);
      double scaleIniInv=
	( iniV[igen] == 0 ) ? 0. : (iniVexact[igen]/iniV[igen]);

      // Now normalize each element and find errors
      for(int ireco = 0; ireco < (*DetMigration).GetNcols(); ireco++){
      // prepare scale in the reconstructed bin
	double scaleFin=
	  ( finVexact[ireco] == 0 ) ? 0. : (finV[ireco]/finVexact[ireco] );

	tCentral = 0;
	tErr     = 0;
	computeNormalizedBinContent((*DetMigration)(igen,ireco),
				    (*DetMigrationErr)(igen,ireco),
				    nEventsInGenBin,
				    nEventsInGenBinErr,
				    tCentral, tErr);
	if ((igen<2) && (ireco<2)) {
	  printf("igen=%d, ireco=%d, tCentral=%8.4lf, tErr=%8.4lf, scaleIniInv=%6.4lf, scaleFin=%6.4lf\n",igen,ireco,tCentral,tErr,scaleIniInv,scaleFin);
	}
	(*DetResponse)      (igen,ireco) = tCentral * scaleIniInv*scaleFin;
	(*DetResponseErrPos)(igen,ireco) = tErr  *scaleIniInv*scaleFin;
	(*DetResponseErrNeg)(igen,ireco) = tErr *scaleIniInv*scaleFin;
      }
    }
  }
  */


  // DET correction factors defined as all.yields/restricted.yields, where
  // all.yields contain events in the acceptance which pass only the relevant
  // pre-FSR or post-FSR cuts, while restricted.yields contain events in the
  // acceptance that passed both pre-FSR and post-FSR cuts (but their masses
  // are either pre-FSR /ini/ or post-FSR /fin/)
  void prepareFsrDETcorrFactors(const UnfoldingMatrix_t &all, const UnfoldingMatrix_t &restricted) {
    for (int iVec=0; iVec<2; ++iVec) {
      TMatrixD *nomV=(iVec==0) ? all.yieldsIni : all.yieldsFin;
      TMatrixD *denomV=(iVec==0) ? restricted.yieldsIni : restricted.yieldsFin;
      TMatrixD *ratioV=(iVec==0) ? yieldsIni : yieldsFin;
      for (int iM=0; iM<nomV->GetNrows(); ++iM) {
	for (int iY=0; iY<nomV->GetNcols(); ++iY) {
	  double d=(*denomV)[iM][iY];
	  double r=(d==0.) ? 0. : ((*nomV)[iM][iY] / d);
	  //std::cout << "divide " << (*nomV)[iM][iY] << " by " << d << " = " << r << "\n";
	  (*ratioV)[iM][iY] = r;
	}
      }
    }
    if (0) {
      std::cout << std::string(80,'-') << "\n";
      std::cout << "prepareFsrDETcorrFactors\n";
      std::cout << std::string(80,'-') << "\n";
      all.printYields();
      restricted.printYields();
      this->printYields();
      std::cout << std::string(80,'-') << "\n";
    }
  }

  void invertResponseMatrix(int NsmearingExperiments=1000) {
  // Find inverted response matrix
    (*DetInvertedResponse) = (*DetResponse);
    Double_t det;
    int info=DYTools::study2D; // slower calculation, so keep the user informed
    (*DetInvertedResponse).Invert(&det);
    if (info) std::cout << name << " inverted" << std::endl;
    //calculateInvertedMatrixErrors(*DetResponse, *DetResponseErrPos, *DetResponseErrNeg, *DetInvertedResponseErr);
    calculateInvertedMatrixErrors(*DetResponse,
			 *DetResponseErrPos, *DetResponseErrPos,
			 *DetInvertedResponseErr, NsmearingExperiments);
    if (info) std::cout << name << " inverted errs calculated" << std::endl;
  }


  // apply factors to the response and inv.response matrices to account
  // for the event loss in DET
  void modifyDETResponseMatrices(const UnfoldingMatrix_t &exact) {
    HERE("modifyDETResponseMatrices");
    TVectorD iniV(DYTools::nUnfoldingBins), finV(DYTools::nUnfoldingBins);
    TVectorD iniVexact(DYTools::nUnfoldingBins), finVexact(DYTools::nUnfoldingBins);
    int resFlatten= 
      (
       (flattenMatrix(*yieldsIni, iniV) == 1 ) &&
       (flattenMatrix(*yieldsFin, finV) == 1 ) &&
       (flattenMatrix(*exact.yieldsIni, iniVexact) == 1 ) &&
       (flattenMatrix(*exact.yieldsFin, finVexact) == 1 ) 
       ) ? 1:0;
    assert(resFlatten);

    if (0) {
      std::cout << "iniV="; iniV.Print();
      std::cout << "finV="; finV.Print();
      std::cout << "iniVexact="; iniVexact.Print();
      std::cout << "finVexact="; finVexact.Print();
    }

    for(int igen = 0; igen < (*DetMigration).GetNrows(); igen++){
      for(int ireco = 0; ireco < (*DetMigration).GetNcols(); ireco++){
	double scaleIni=
	  ( iniVexact[igen] == 0 ) ? 0. : (iniV[igen]/iniVexact[igen]);
	double scaleIniInv=
	  ( iniV[igen] == 0 ) ? 0. : (iniVexact[igen]/iniV[igen]);
	double scaleFin=
	  ( finVexact[ireco] == 0 ) ? 0. : (finV[ireco]/finVexact[ireco] );
	double scaleFinInv=
	  ( finV[ireco] == 0 ) ? 0. : (finVexact[ireco]/finV[ireco] );

	(*DetResponse)      (igen,ireco) *= scaleIniInv*scaleFin;
	(*DetResponseErrPos)(igen,ireco) *= scaleIniInv*scaleFin;
	(*DetResponseErrNeg)(igen,ireco) *= scaleIniInv*scaleFin;
	(*DetInvertedResponse)   (ireco,igen) *= scaleIni*scaleFinInv;
	(*DetInvertedResponseErr)(ireco,igen) *= scaleIni*scaleFinInv;
      }
    }
  }


  void prepareFIArrays() {
    int resFlatten=
      (flattenMatrix(*DetResponse, *DetResponseArr) == 1) &&
      (flattenMatrix(*DetInvertedResponse, *DetInvertedResponseArr) == 1) &&
      (flattenMatrix(*DetInvertedResponseErr, *DetInvertedResponseErrArr) == 1) &&
      (flattenMatrix(*yieldsIni, *yieldsIniArr) == 1) &&
      (flattenMatrix(*yieldsFin, *yieldsFinArr) == 1);
    if (!resFlatten) {
      std::cout << "Error : failed to flatten the arrays\n";
      assert(0);
    }
  }


  static void getYieldNames(UnfoldingMatrix::TUnfoldingMatrixType_t theKind, TString &iniName, TString &finName) { return UnfoldingMatrix::getYieldNames(theKind,iniName,finName); }

  void getFileNames(const TString &outputDir,
		    const TString &fileTag,
		    TString &matrixFName, TString &yieldsFName) const {
    matrixFName=outputDir + this->name + TString("_unfolding_constants") + 
      fileTag + TString(".root");
    yieldsFName=outputDir + TString("yields_") + this->name + 
      fileTag + TString(".root");
  }

  // - ------------------

  static TString generateFNameTag(DYTools::TSystematicsStudy_t systMode,
				  int seed) {
    TString u="_";
    TString fnameTag="UNKNOWN";
    switch(systMode) {
    case DYTools::NO_SYST: 
      fnameTag=DYTools::analysisTag; 
      break;
    /*
    case DYTools::SYST_RND: 
      fnameTag=TString("_replica_") + DYTools::analysisTag; 
      //saveIdxMin=0;
      //fnameTag+=seed;
      break;
    */
    case DYTools::RESOLUTION_STUDY: 
      fnameTag=TString("_seed_") + DYTools::analysisTag;
      //fnameTag+=seed;
      break;
    case DYTools::FSR_RND_STUDY:
    case DYTools::PU_RND_STUDY:
      fnameTag=TString(Form("seed%d_",seed)) + DYTools::analysisTag;
      break;
    case DYTools::FSR_STUDY:
    case DYTools::FSR_5plus:
    case DYTools::FSR_5minus:
      fnameTag=TString("_fsrStudy_") + DYTools::analysisTag;
      //fnameTag=TString("_reweight_") + DYTools::analysisTag;
      //fnameTag+= int(100*reweightFsr);
      break;
    case DYTools::PU_STUDY:
    case DYTools::PILEUP_5plus:
    case DYTools::PILEUP_5minus:
      fnameTag=TString("_puStudy_") + DYTools::analysisTag;
      break;
    /*
    case DYTools::ESCALE_STUDY:
      fnameTag=DYTools::analysisTag+TString("_escale") + u;
      break;
    */
    case DYTools::ESCALE_RESIDUAL:
      fnameTag=DYTools::analysisTag+TString("_escaleResidual");
      break;
    default:
      std::cout<<"\n\tERRORrequested mode not recognized when determining fnameTag"<<std::endl;
    }
    return fnameTag;
  }

  void autoSaveToFile(const TString &outputDir, const TString &fileTag) const {
    TString matrixFName,yieldsFName;
    this->getFileNames(outputDir,fileTag, matrixFName,yieldsFName);
    std::cout << "saving to files <" << matrixFName << "> and <" << yieldsFName << ">\n";
    this->saveToFile(matrixFName,yieldsFName);
  }
  
  int autoLoadFromFile(const TString &outputDir, const TString &fileTag) {
    TString matrixFName,yieldsFName;
    this->getFileNames(outputDir,fileTag, matrixFName,yieldsFName);
    std::cout << "loading from files <" << matrixFName << "> and <" << yieldsFName << ">\n";
    return this->loadFromFile(matrixFName,yieldsFName);
  }
  

  void saveToFile(const TString &fileName, const TString &refFileName) const {
    std::cout << "UnfoldingMatrix_t::saveToFile(\n  <" << fileName << ">\n  <" << refFileName << ">) for name=" << this->name << "\n";
    if (kind!=UnfoldingMatrix::_cFSR_DETcorrFactors) {
      TFile fConst(fileName, "recreate" );
      //name.Write("matrixName");
      (*DetMigration)            .Write("DetMigration");
      (*DetMigrationErr)         .Write("DetMigrationErr");
      (*DetResponse)             .Write("DetResponse");
      (*DetResponseErrPos)       .Write("DetResponseErrPos");
      (*DetResponseErrNeg)       .Write("DetResponseErrNeg");
      (*DetInvertedResponse)     .Write("DetInvertedResponse");
      (*DetInvertedResponseErr)  .Write("DetInvertedResponseErr");
      (*DetResponseArr)          .Write("DetResponseFIArray");
      (*DetInvertedResponseArr)  .Write("DetInvertedResponseFIArray");
      (*DetInvertedResponseErrArr).Write("DetInvertedResponseErrFIArray");
      (*yieldsIni).Write(iniYieldsName);
      (*yieldsFin).Write(finYieldsName);
      (*yieldsIniArr).Write(iniYieldsName + TString("FIArray"));
      (*yieldsFinArr).Write(finYieldsName + TString("FIArray"));
      writeBinningArrays(fConst);
      fConst.Close();
    }

    // Store reference MC arrays in a file
    TFile fRef(refFileName, "recreate" );
    (*yieldsIni).Write(iniYieldsName);
    (*yieldsFin).Write(finYieldsName);
    (*yieldsIniArr).Write(iniYieldsName + TString("FIArray"));
    (*yieldsFinArr).Write(finYieldsName + TString("FIArray"));
    writeBinningArrays(fRef);
    fRef.Close();
  }

  // ------------------------------------------

  int loadFromFile(const TString &fileName, const TString &refFileName) {
    std::cout << "UnfoldingMatrix_t::loadFromFile(\n  <" << fileName << ">\n  <" << refFileName << ">) for name=" << this->name << "\n";
    if (kind!=UnfoldingMatrix::_cFSR_DETcorrFactors) {
      TFile fConst(fileName);
      if (!fConst.IsOpen()) {
	std::cout << "failed to open the file <" << fileName << ">\n";
	return 0;
      }
      if (!checkBinningArrays(fConst)) {
	fConst.Close();
	return 0;
      }
      (*DetMigration)            .Read("DetMigration");
      (*DetMigrationErr)         .Read("DetMigrationErr");
      (*DetResponse)             .Read("DetResponse");
      (*DetResponseErrPos)       .Read("DetResponseErrPos");
      (*DetResponseErrNeg)       .Read("DetResponseErrNeg");
      (*DetInvertedResponse)     .Read("DetInvertedResponse");
      (*DetInvertedResponseErr)  .Read("DetInvertedResponseErr");
      (*DetResponseArr)          .Read("DetResponseFIArray");
      (*DetInvertedResponseArr)  .Read("DetInvertedResponseFIArray");
      (*DetInvertedResponseErrArr).Read("DetInvertedResponseErrFIArray");
      (*yieldsIni).Read(iniYieldsName);
      (*yieldsFin).Read(finYieldsName);
      (*yieldsIniArr).Read(iniYieldsName + TString("FIArray"));
      (*yieldsFinArr).Read(finYieldsName + TString("FIArray"));
      fConst.Close();
    }

    if (kind==UnfoldingMatrix::_cFSR_DETcorrFactors) {
      // Retrieve reference MC arrays in a file
      TFile fRef(refFileName);
      if (!fRef.IsOpen()) {
	std::cout << "failed to open the file <" << refFileName << ">\n";
	return 0;
      }
      if (!checkBinningArrays(fRef)) {
	fRef.Close();
	return 0;
      }
      (*yieldsIni).Read(iniYieldsName);
      (*yieldsFin).Read(finYieldsName);
      (*yieldsIniArr).Read(iniYieldsName + TString("FIArray"));
      (*yieldsFinArr).Read(finYieldsName + TString("FIArray"));
      fRef.Close();
    }
    return 1;
  }

  // ------------------------------------------

  int loadFromFileAny(const TString &fileName, const TString &refFileName,
		      TVectorD **massBinning, TVectorD **nRapidityBins) {
    std::cout << "UnfoldingMatrix_t::loadFromFileAny(\n  <" << fileName << ">\n  <" << refFileName << ">) for name=" << this->name << "\n";
    if (!massBinning || !nRapidityBins) {
      std::cout << "loadFromFileAny: the provided massBinning and nRapidityBins pointers are null\n";
      return 0;
    }
    if (kind!=UnfoldingMatrix::_cFSR_DETcorrFactors) {
      TFile fConst(fileName);
      if (!fConst.IsOpen()) {
	std::cout << "failed to open the file <" << fileName << ">\n";
	return 0;
      }
      if (!checkBinningArrays(fConst)) {
	fConst.Close();
	return 0;
      }
      (*massBinning) = (TVectorD*)fConst.Get("massBinLimits");
      (*nRapidityBins) = (TVectorD*)fConst.Get("rapidityCounts");

      (*DetMigration)            .Read("DetMigration");
      (*DetMigrationErr)         .Read("DetMigrationErr");
      (*DetResponse)             .Read("DetResponse");
      (*DetResponseErrPos)       .Read("DetResponseErrPos");
      (*DetResponseErrNeg)       .Read("DetResponseErrNeg");
      (*DetInvertedResponse)     .Read("DetInvertedResponse");
      (*DetInvertedResponseErr)  .Read("DetInvertedResponseErr");
      (*DetResponseArr)          .Read("DetResponseFIArray");
      (*DetInvertedResponseArr)  .Read("DetInvertedResponseFIArray");
      (*DetInvertedResponseErrArr).Read("DetInvertedResponseErrFIArray");
      (*yieldsIni).Read(iniYieldsName);
      (*yieldsFin).Read(finYieldsName);
      (*yieldsIniArr).Read(iniYieldsName + TString("FIArray"));
      (*yieldsFinArr).Read(finYieldsName + TString("FIArray"));
      fConst.Close();
    }

    if (kind==UnfoldingMatrix::_cFSR_DETcorrFactors) {
      // Retrieve reference MC arrays in a file
      TFile fRef(refFileName);
      if (!fRef.IsOpen()) {
	std::cout << "failed to open the file <" << refFileName << ">\n";
	return 0;
      }
      if (!checkBinningArrays(fRef)) {
	fRef.Close();
	return 0;
      }
      (*massBinning) = (TVectorD*)fRef.Get("massBinLimits");
      (*nRapidityBins) = (TVectorD*)fRef.Get("rapidityCounts");

      (*yieldsIni).Read(iniYieldsName);
      (*yieldsFin).Read(finYieldsName);
      (*yieldsIniArr).Read(iniYieldsName + TString("FIArray"));
      (*yieldsFinArr).Read(finYieldsName + TString("FIArray"));
      fRef.Close();
    }
    return 1;
  }

  // ------------------------------------------
  // ------------------------------------------

  static TMatrixD* LoadUnfM(//UnfoldingMatrix::TUnfoldingMatrixType_t set_kind, 
			    const TString &set_name,
			    const TString &outputDir,
			    const TString &fileTag,
			    int inverse,
			    int checkBinArrays=1
			    ) {
    // construct file name
    // emulate this->getFileNames()
    TString matrixFName=
      outputDir + set_name + TString("_unfolding_constants") +
      fileTag + TString(".root");
    TFile fin(matrixFName,"read");
    if (!fin.IsOpen()) {
      std::cout << "loadUnfM: failed to open file <" << fin.GetName() << ">\n";
      return NULL;
    }
    if (checkBinArrays && !checkBinningArrays(fin)) {
      std::cout << "loadUnfM(" << fin.GetName() << ") does not match the binning (checkBinArrays is ON)\n";
      return NULL;
    }
    TString fieldName= (inverse) ? "DetInvertedResponse" : "DetResponse";
    TMatrixD *unf=(TMatrixD*) fin.Get(fieldName);
    fin.Close();
    if (!unf) std::cout << "loadUnfM(" << fin.GetName() << "): failed to load field<" << fieldName << ">\n";
    return unf;
  }

  // ------------------------------------------
  // ------------------------------------------

  void printYields() const {
    std::cout << "Yields of matrix=" << name << " (" 
	      << iniYieldsName << " and " << finYieldsName << ")\n";
    for (int ir=0; ir<(*yieldsIni).GetNrows(); ++ir) {
      std::cout << "ir=" << ir << "\n";
      for (int ic=0; ic<(*yieldsIni).GetNcols(); ++ic) {
	printf(" % 9.6lf  % 9.6lf\n",(*yieldsIni)[ir][ic],(*yieldsFin)[ir][ic]);
      }
      if (DYTools::study2D) printf("\n");
    }
  }

  void printMigration() const {
    std::cout << "DetMigration of <" << name << ">:\n";
    int printSystErr=0;
    TMatrixD zeroErr=*DetMigrationErr;
    zeroErr=0;
    printCSMatrixValues("DetMigration",*DetMigration,*DetMigrationErr,zeroErr,printSystErr);
  }

  void printResponse() const {
    std::cout << "DetResponse,ErrPos,ErrNeg of <" << name << ">:\n";
    int printSystErr=1;
    printCSMatrixValues("DetResponse",*DetResponse,*DetResponseErrPos,*DetResponseErrNeg,printSystErr);
  }

  void printInvResponse() const {
    std::cout << "DetInvertedResponse of <" << name << ">:\n";
    int printSystErr=0;
    TMatrixD zeroErr=*DetInvertedResponseErr;
    zeroErr=0;
    printCSMatrixValues("DetInvertedResponse",*DetInvertedResponse,*DetInvertedResponseErr,zeroErr,printSystErr);
  }

  void printMatrices() const {
    std::string line(80,'-');
    std::cout << "\n" << line << "\n";
    printMigration();
    printResponse();
    printInvResponse();
    std::cout << line << "\n";
  }

  void printConditionNumber() const {
    //matrix condition number
    TDecompLU lu(*DetResponse);
    double condLU=lu.Condition();
    std::cout << "Matrix=" << name << "\n";
    std::cout << " condition number from TDecompLU condLU= " << condLU << std::endl;
    std::cout << " condition number ||DetResponse||*||DetResponseInv||=" << DetResponse->Norm1()*DetInvertedResponse->Norm1() << std::endl;
    std::cout << " chk ROOT bug: -condLU*||DetResponse||=" << (-condLU*DetResponse->Norm1()) << "\n" << std::endl;
  }

  void prepareHResponse(TH2D **hResponse_out=NULL,
			TH2D **hInvResponse_out=NULL,
			TCanvas **canv=NULL,
			CPlot **plotResponse_out=NULL,
			CPlot **plotInvResponse_out=NULL
			) {
    // Plot response and inverted response matrices
    TString kName=this->ourKindName();
    kName.Append("_"); kName.Append(this->name);
    TH2D *hResponse = new TH2D(TString("hResponse_") + kName,"",
			       DYTools::nUnfoldingBins, -0.5, DYTools::nUnfoldingBins-0.5,
			       DYTools::nUnfoldingBins, -0.5, DYTools::nUnfoldingBins-0.5);
    TH2D *hInvResponse = new TH2D(TString("hInvResponse") + kName,"",
				  DYTools::nUnfoldingBins, -0.5, DYTools::nUnfoldingBins-0.5,
				  DYTools::nUnfoldingBins, -0.5, DYTools::nUnfoldingBins-0.5);
    hResponse->SetDirectory(0);
    hInvResponse->SetDirectory(0);
    for(int i=0; i<(*DetResponse).GetNrows(); i++){
      for(int j=0; j<(*DetResponse).GetNcols(); j++){
	hResponse->SetBinContent(i,j, (*DetResponse)(i,j));
	hInvResponse->SetBinContent(i,j, (*DetInvertedResponse)(i,j));
      }
    }
    hResponse->GetYaxis()->SetTitleOffset(1.1);
    hInvResponse->GetYaxis()->SetTitleOffset(1.1);


    TString canvName=TString("canvResponse") + kName;
    TCanvas *e1 = MakeCanvas(canvName,canvName,1200,600);
    e1->Divide(2,1);
    AdjustFor2DplotWithHeight(e1);
    CPlot *plotResponse=
      new CPlot(TString("response") + kName,"",
		"flat index gen",
		"flat index reco");
    plotResponse->AddHist2D(hResponse,"COLZ");
    plotResponse->Draw(e1,false,"png",1);

    CPlot *plotInvResponse=
      new CPlot(TString("invResponse") + kName,"",
		"flat index reco",
		"flat index gen");
    plotInvResponse->AddHist2D(hInvResponse,"COLZ");
    plotInvResponse->Draw(e1,false,"png",2);
    e1->Update();
    SaveCanvas(e1,Form("hResponse_%s_",DYTools::analysisTag.Data()) + kName);
  
    if (hResponse_out) *hResponse_out=hResponse;
    if (hInvResponse_out) *hInvResponse_out=hInvResponse;
    if (canv) *canv=e1;
    if (plotResponse_out) *plotResponse_out=plotResponse;
    if (plotInvResponse_out) *plotInvResponse_out=plotInvResponse;
  }

  // ------------------------------------------------------

  TMatrixD* getReconstructionEffect(const UnfoldingMatrix_t &inexact) const {
    TMatrixD *res=new TMatrixD(*yieldsIni);
    *res=0;
    for(int i=0; i < res->GetNrows(); i++){
      for(int j=0; j < res->GetNcols(); j++){
	double nexact = (*yieldsIni)(i,j);
	double nactual = (*inexact.yieldsIni)(i,j);
	if( nexact != 0 )
	  (*res)(i,j) = (nexact-nactual)/nexact;
      }
    }
    return res;
  }

  // ------------------------------------------------------

  // ------------------------------------------------------

};


// -----------------------------------------------
//               More functions: unfolding
// -----------------------------------------------

inline
int unfold(TVectorD &finVec, const TMatrixD &U, const TVectorD &iniVec) {
  int res=1;
  if (finVec.GetNoElements() != iniVec.GetNoElements() ) {
    std::cout << "Dim error in unfold: finVec[" << finVec.GetNoElements() << "], iniVec[" << iniVec.GetNoElements() << "]\n";
    res=0;
  }
  finVec.Zero();
  if ( (finVec.GetNoElements() != U.GetNrows()) ||
       (finVec.GetNoElements() != U.GetNcols()) ) {
    std::cout << "Dim error in unfold: finVec[" << finVec.GetNoElements() << "], U[" << U.GetNrows() << ", " << U.GetNcols() << "]\n";
    res=0;
  }
  if (res) {
    for (int i=0; i<finVec.GetNoElements(); i++) {
      double sum=0;
      for (int j=0; j<finVec.GetNoElements(); j++) {
	sum += U(j,i) * iniVec(j);
      }
      finVec(i) = sum;
    }
  }
  return res;
}

// -----------------------------------------------

inline
int unfold_true2reco(TVectorD &vecReco, const UnfoldingMatrix_t &U, const TVectorD &vecTrue) {
  return unfold(vecReco, *U.getDetResponse(), vecTrue);
}
// -----------------------------------------------

inline
int unfold_reco2true(TVectorD &vecTrue, const UnfoldingMatrix_t &U, const TVectorD &vecReco) {
  return unfold(vecTrue, *U.getDetInvResponse(), vecReco);
}

// -----------------------------------------------
// -----------------------------------------------

inline
int unfold(TMatrixD &finM, const TMatrixD &U, const TMatrixD &iniM) {
  int res=1;
  if ((finM.GetNrows() != iniM.GetNrows() ) ||
      (finM.GetNcols() != iniM.GetNcols()) ) {
    std::cout << "Dim error in unfold(Matrix): finM[" << finM.GetNrows() << ", " << finM.GetNcols() << "], iniM[" << iniM.GetNrows() << ", " << iniM.GetNcols() << "]\n";
    res=0;
  }
  if (finM.GetNrows()*finM.GetNcols() < DYTools::nUnfoldingBins) {
    std::cout << "Dim error in unfold(Matrix): finM[" << finM.GetNrows() << ", " << finM.GetNcols() << "], nUnfoldingBins=" << DYTools::nUnfoldingBins << "\n";
    res=0;
  }
  if ( (DYTools::nUnfoldingBins != U.GetNrows()) ||
       (DYTools::nUnfoldingBins != U.GetNcols()) ) {
    std::cout << "Dim error in unfold(Matrix): nUnfoldingBins=" << DYTools::nUnfoldingBins << ", U[" << U.GetNrows() << ", " << U.GetNcols() << "]\n";
    res=0;
  }
  finM.Zero();
  if (res) {
    for (int ir=0, ini=0; ir<finM.GetNrows(); ir++) {
      if (finM.GetNcols()<DYTools::nYBins[ir]) {
	std::cout << "Dim problem in unfold(Matrix) at ir=" << ir << ", finM.GetNcols=" << finM.GetNcols() << ", nYBins=" << DYTools::nYBins[ir] << "\n";
      }
      for (int ic=0; 
	   //(ic<finM.GetNcols()) && (ini<DYTools::nUnfoldingBins); 
	   (ic<finM.GetNcols()) &&
	   (ic<DYTools::nYBins[ir]) && (ini<DYTools::nUnfoldingBins);
	   ++ic, ++ini) {
	double sum=0;
	for (int jr=0, fin=0; jr<iniM.GetNrows(); jr++) {
	  for (int jc=0;
	       (jc<iniM.GetNcols()) &&
	       (jc<DYTools::nYBins[jr]) && (fin<DYTools::nUnfoldingBins);
	       ++jc, ++fin) {
	    sum += U(fin,ini) * iniM(jr,jc);
	  }
	}
	finM(ir,ic) = sum;
      }
    }
  }
  return res;
}

// -----------------------------------------------

inline
int unfold_true2reco(TMatrixD &matReco, const UnfoldingMatrix_t &U, const TMatrixD &matTrue) {
  return unfold(matReco, *U.getDetResponse(), matTrue);
}
// -----------------------------------------------

inline
int unfold_reco2true(TMatrixD &matTrue, const UnfoldingMatrix_t &U, const TMatrixD &matReco) {
  return unfold(matTrue, *U.getDetInvResponse(), matReco);
}

// -----------------------------------------------

// -----------------------------------------------

// -----------------------------------------------
// -----------------------------------------------

#ifdef HistoPair_HH
inline
int unfold(HistoPair2D_t &finHP, const TMatrixD &U, const HistoPair2D_t &iniHP) {
  int res=1;
  if ((finHP.getNrows() != iniHP.getNrows() ) ||
      (finHP.getNcols() != iniHP.getNcols()) ) {
    std::cout << "Dim error in unfold(HistoPair2D): finHP[" << finHP.getNrows() << ", " << finHP.getNcols() << "], iniHP[" << iniHP.getNrows() << ", " << iniHP.getNcols() << "]\n";
    res=0;
  }
  if (finHP.getNrows()*finHP.getNcols() < DYTools::nUnfoldingBins) {
    std::cout << "Dim error in unfold(HistoPair2D): finHP[" << finHP.getNrows() << ", " << finHP.getNcols() << "], nUnfoldingBins=" << DYTools::nUnfoldingBins << "\n";
    res=0;
  }
  if ( (DYTools::nUnfoldingBins != U.GetNrows()) ||
       (DYTools::nUnfoldingBins != U.GetNcols()) ) {
    std::cout << "Dim error in unfold(HistoPair2D): nUnfoldingBins=" << DYTools::nUnfoldingBins << ", U[" << U.GetNrows() << ", " << U.GetNcols() << "]\n";
    res=0;
  }
  finHP.Reset();
  if (res) {
    for (int ir=0, ini=0; ir<finHP.getNrows(); ir++) {
      if (finHP.getNcols()<DYTools::nYBins[ir]) {
	std::cout << "Dim problem in unfold(HistoPair2D) at ir=" << ir << ", finM.GetNcols=" << finHP.getNcols() << ", nYBins=" << DYTools::nYBins[ir] << "\n";
      }
      for (int ic=0; 
	   //(ic<finHPM.GetNcols()) && (ini<DYTools::nUnfoldingBins); 
	   (ic<DYTools::nYBins[ir]) && (ini<DYTools::nUnfoldingBins);
	   ++ic, ++ini) {
	double sum=0;
	double sumErr2=0;
	double sumSystErr2=0;
	for (int jr=0, fin=0; jr<iniHP.getNrows(); jr++) {
	  for (int jc=0;
	       //(jc<iniM.GetNcols()) && (fin<DYTools::nUnfoldingBins);
	       (jc<DYTools::nYBins[jr]) && (fin<DYTools::nUnfoldingBins);
	       ++jc, ++fin) {
	    sum += U(fin,ini) * iniHP.getBinContent(jr+1,jc+1);
	    const double tmpErr= U(fin,ini) * iniHP.getBinError(jr+1,jc+1);
	    sumErr2 += tmpErr*tmpErr;
	    const double tmpSystErr= U(fin,ini) * iniHP.getBinSystError(jr+1,jc+1);
	    sumSystErr2 += tmpSystErr*tmpSystErr;
	  }
	}
	finHP.setBinContent(ir+1,ic+1, sum);
	finHP.setBinError(ir+1, ic+1, sqrt(sumErr2));
	finHP.setBinSystError(ir+1, ic+1, sqrt(sumSystErr2));
      }
    }
  }
  return res;
}
#endif

// -----------------------------------------------

#ifdef HistoPair_HH
inline
int unfold_true2reco(HistoPair2D_t &hpReco, const UnfoldingMatrix_t &U, const HistoPair2D_t &hpTrue) {
  return unfold(hpReco, *U.getDetResponse(), hpTrue);
}
#endif

// -----------------------------------------------

#ifdef HistoPair_HH
inline
int unfold_reco2true(HistoPair2D_t &hpTrue, const UnfoldingMatrix_t &U, const HistoPair2D_t &hpReco) {
  return unfold(hpTrue, *U.getDetInvResponse(), hpReco);
}
#endif

// -----------------------------------------------
// -----------------------------------------------
// -----------------------------------------------

#endif
