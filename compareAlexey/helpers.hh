#ifndef helpers_HH
#define helpers_HH

#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include <fstream>
#include <sstream>

// -----------------------------------------------------------------------

typedef enum { _cmp_RawYield, _cmp_FakeBkg, _cmp_TrueBkg,
	       _cmp_Eff, _cmp_MCeff, _cmp_EffRho, _cmp_Acc,
	       _cmp_UnfYields } TCompareCase_t;

typedef enum { _cmp_fsrGood, _cmp_fsrExact } TFsrUnfCompareCase_t;

// -----------------------------------------------------------------------

TH2D *loadTextFile(TString fname, TString histoName,
		   int hassError=1, int hasDivider=1);
TH2D *loadXYZTextFile(TString fname, TString histoName, int hasError,
		      int transposed=0);

TH2D *loadHisto1D_convert_TH2D(TString fname, TString histoName);

void printProfileSums(const TH2D *h2, TVectorD *sumX=NULL,TVectorD *sumY=NULL);
void compareProfileSums(const TH2D *h2A, const TH2D *h2B);

// -----------------------------------------------------------------------
// -----------------------------------------------------------------------

#endif
