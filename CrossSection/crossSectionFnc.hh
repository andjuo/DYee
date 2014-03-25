#ifndef crossSectionFnc_HH
#define crossSectionFnc_HH

#include <TBenchmark.h>

#include "../Include/DYTools.hh"
#include "../Include/DYToolsUI.hh"
#include "../Include/MyTools.hh"

#include "../Include/EventSelector.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/HistoPair.hh"
#include "../Include/UnfoldingMatrix.h"

extern TString codeDebugFilePath;


// ----------------------------------------------
// ----------------------------------------------

struct InputArgs_t {
  InputFileMgr_t *fInpMgr;
  DYTools::TSystematicsStudy_t fSystMode;
  TString fResNameBase;
  int fNeedsDETUnfolding;
  int fAllNormErrorIsSyst;
public:
  InputArgs_t(InputFileMgr_t *set_InpMgr, DYTools::TSystematicsStudy_t set_systMode,
	      TString set_resNameBase="",
	      int set_needsDetUnfolding=1,
	      int set_allNormErrorIsSyst=0) :
    fInpMgr(set_InpMgr), 
    fSystMode(set_systMode),
    fResNameBase(set_resNameBase), 
    fNeedsDETUnfolding(set_needsDetUnfolding),
    fAllNormErrorIsSyst(set_allNormErrorIsSyst)
  {}

  InputArgs_t(const InputArgs_t &ia, 
	      TString set_resNameBase, int set_needsDetUnfolding=-1,
	      int set_allNormErrorIsSyst=-1) :
    fInpMgr(ia.fInpMgr), 
    fSystMode(ia.fSystMode),
    fResNameBase(set_resNameBase),
    fNeedsDETUnfolding(ia.fNeedsDETUnfolding),
    fAllNormErrorIsSyst(ia.fAllNormErrorIsSyst)
  {
    if (set_needsDetUnfolding!=-1) fNeedsDETUnfolding=set_needsDetUnfolding;
    if (set_allNormErrorIsSyst!=-1) fAllNormErrorIsSyst=set_allNormErrorIsSyst;
  }
  
  const InputFileMgr_t* inpMgr() const { return fInpMgr; }
  DYTools::TSystematicsStudy_t systMode() const { return fSystMode; }
  TString resNameBase() const { return fResNameBase; }

  void needsDetUnfolding(int yes) { fNeedsDETUnfolding=yes; }
  int  needsDetUnfolding() const { return fNeedsDETUnfolding; }
  void allNormErrorIsSyst(int yes) { fAllNormErrorIsSyst=yes; }
  int  allNormErrorIsSyst() const { return fAllNormErrorIsSyst; }

};

// ----------------------------------------------
// ----------------------------------------------

struct CSResults_t {
  double fZXSec,fZXSecErr,fZXSecSystErr;
public:
  CSResults_t() : fZXSec(0.), fZXSecErr(0.), fZXSecSystErr(0.) {}
  
  CSResults_t(const CSResults_t &r) :
    fZXSec(r.fZXSec), fZXSecErr(r.fZXSecErr), fZXSecSystErr(r.fZXSecSystErr) {}

  double xs() const { return fZXSec; }
  double xsErr() const { return fZXSecErr; }
  double xsErrSyst() const { return fZXSecSystErr; }
  double xsSystErr() const { return fZXSecSystErr; }

  void assignCS(double set_xs, double set_xsErr, double set_xsSystErr) {
    fZXSec=set_xs; fZXSecErr=set_xsErr; fZXSecSystErr=set_xsSystErr;
  }

  friend std::ostream& operator<<(std::ostream& out, const CSResults_t &r) {
    out << "CSResults(" << Form("%7.4lf +- %7.4lf +- %7.4lf",r.fZXSec,r.fZXSecErr,r.fZXSecSystErr) << ")";
    return out;
  }
};

// ----------------------------------------------
// ----------------------------------------------


//=== FUNCTION DECLARATIONS ======================================================================================

//int getNormalizationMBinRange(int &firstNormBin, int &lastNormBin);

//=== CS Calculation FUNCTION DECLARATIONS ======================================================================================

// step 1
int unfoldDetResolution(const InputArgs_t &ia, const HistoPair2D_t &ini, HistoPair2D_t &fin);
// step 2
int efficiencyCorrection(const InputArgs_t &ia, const HistoPair2D_t &ini, HistoPair2D_t &fin);
// step 3
int efficiencyScaleCorrection(const InputArgs_t &ia, const HistoPair2D_t &ini, HistoPair2D_t &fin);
// step 4
int acceptanceCorrection(const InputArgs_t &ia, const HistoPair2D_t &ini, HistoPair2D_t &fin);
// step 5
int fsrCorrection_det(const InputArgs_t &ia, const HistoPair2D_t &ini, HistoPair2D_t &fin);
int fsrCorrection_fullSpace(const InputArgs_t &ia, const HistoPair2D_t &ini, HistoPair2D_t &fin);

// systematics
// of step 1
int addSystError_DetResUnf_unfold(const InputArgs_t &ia, HistoPair2D_t &hp, HistoPair2D_t *res=NULL);
int addSystError_DetResUnf_escale(const InputArgs_t &ia, HistoPair2D_t &hp, HistoPair2D_t *res=NULL);


// saving the result
int saveResult(const InputArgs_t &ia, const HistoPair2D_t &hp, const TString &extraTag);



// cross section calculation
int calculateCSdistribution(const InputArgs_t &ia, const HistoPair2D_t &hp_ini, 
			    DYTools::TCrossSectionKind_t csKind,
			    HistoPair2D_t &hp_fin);
int calculateCS(const InputArgs_t &ia, const HistoPair2D_t &hp_ini, 
		DYTools::TCrossSectionKind_t csKind,
		HistoPair2D_t &hp_fin, CSResults_t &res);

#endif
