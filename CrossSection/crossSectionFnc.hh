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



// ----------------------------------------------
// ----------------------------------------------

struct InputArgs_t {
  InputFileMgr_t *inpMgr;
  DYTools::TSystematicsStudy_t systMode;
  TString resNameBase;
  int needsDETUnfolding;
public:
  InputArgs_t(InputFileMgr_t *set_InpMgr, DYTools::TSystematicsStudy_t set_systMode,
	      TString set_resNameBase="",
	      int set_needsDetUnfolding=1) :
    inpMgr(set_InpMgr), systMode(set_systMode),
    resNameBase(set_resNameBase), 
    needsDETUnfolding(set_needsDetUnfolding)
  {}

  InputArgs_t(const InputArgs_t &ia, 
	      TString set_resNameBase, int set_needsDetUnfolding=-1) :
    inpMgr(ia.inpMgr), systMode(ia.systMode),
    resNameBase(set_resNameBase),
    needsDETUnfolding(ia.needsDETUnfolding) 
  {
    if (set_needsDetUnfolding!=-1) needsDETUnfolding=set_needsDetUnfolding;
  }
    
};

// ----------------------------------------------
// ----------------------------------------------

struct CSResults_t {
  double ZXSec,ZXSecErr,ZXSecSystErr;
public:
  CSResults_t() : ZXSec(0.), ZXSecErr(0.), ZXSecSystErr(0.) {}
  
  CSResults_t(const CSResults_t &r) :
    ZXSec(r.ZXSec), ZXSecErr(r.ZXSecErr), ZXSecSystErr(r.ZXSecSystErr) {}

  void assignCS(double xs, double xsErr, double xsSystErr) {
    ZXSec=xs; ZXSecErr=xsErr; ZXSecSystErr=xsSystErr;
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
