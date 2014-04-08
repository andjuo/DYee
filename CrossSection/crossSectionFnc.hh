#ifndef crossSectionFnc_HH
#define crossSectionFnc_HH

#include <TBenchmark.h>

#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"

#include "../Include/EventSelector.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/HistoPair.hh"
#include "../Include/UnfoldingMatrix.h"

extern TString codeDebugFilePath;


// ----------------------------------------------
// ----------------------------------------------

struct InputArgs_t {
  TString fName;
  InputFileMgr_t *fInpMgr;
  DYTools::TSystematicsStudy_t fSystMode;
  DYTools::TCrossSectionKind_t fCSKind;
  TString fResNameBase;
  int fNeedsDETUnfolding;
  int fNeedsEffCorr;
  int fNeedsEffScaleCorr;
  int fNeedsAccCorr;
  int fNeedsFsrCorr;
  int fAllNormErrorIsSyst;
  int fIncludeCorrError; // whether add the error of each correction
  int fSilentMode;
  int fNoSave;
public:
  InputArgs_t(TString set_name,
	      InputFileMgr_t *set_InpMgr,
	      DYTools::TSystematicsStudy_t set_systMode,
	      DYTools::TCrossSectionKind_t set_csKind,
	      TString set_resNameBase="",
	      int set_needsDetUnfolding=1,
	      int set_allNormErrorIsSyst=0,
	      int set_includeCorrErr=1) :
    fName(set_name),
    fInpMgr(set_InpMgr), 
    fSystMode(set_systMode),
    fCSKind(set_csKind),
    fResNameBase(set_resNameBase), 
    fNeedsDETUnfolding(set_needsDetUnfolding),
    fNeedsEffCorr(1),
    fNeedsEffScaleCorr(1),
    fNeedsAccCorr(1),
    fNeedsFsrCorr(1),
    fAllNormErrorIsSyst(set_allNormErrorIsSyst),
    fIncludeCorrError(set_includeCorrErr),
    fSilentMode(0),
    fNoSave(0)
  {}

  InputArgs_t(TString set_name,
	      const InputArgs_t &ia,
	      TString set_resNameBase, int set_needsDetUnfolding=-1,
	      int set_allNormErrorIsSyst=-1,
	      int set_includeCorrError=-1) :
    fName(set_name),
    fInpMgr(ia.fInpMgr), 
    fSystMode(ia.fSystMode),
    fCSKind(ia.fCSKind),
    fResNameBase(set_resNameBase),
    fNeedsDETUnfolding(ia.fNeedsDETUnfolding),
    fNeedsEffCorr(ia.fNeedsEffCorr),
    fNeedsEffScaleCorr(ia.fNeedsEffScaleCorr),
    fNeedsAccCorr(ia.fNeedsAccCorr),
    fNeedsFsrCorr(ia.fNeedsFsrCorr),
    fAllNormErrorIsSyst(ia.fAllNormErrorIsSyst),
    fIncludeCorrError(ia.fIncludeCorrError),
    fSilentMode(ia.fSilentMode),
    fNoSave(ia.fNoSave)
  {
    if (set_needsDetUnfolding!=-1) fNeedsDETUnfolding=set_needsDetUnfolding;
    if (set_allNormErrorIsSyst!=-1) fAllNormErrorIsSyst=set_allNormErrorIsSyst;
    if (set_includeCorrError!=-1) fIncludeCorrError=set_includeCorrError;
  }
  
  TString name() const { return fName; }
  TString GetName() const { return fName; }
  void name(const TString &new_name) { fName=new_name; }
  void SetName(const TString &new_name) { fName=new_name; }

  const InputFileMgr_t* inpMgr() const { return fInpMgr; }
  DYTools::TSystematicsStudy_t systMode() const { return fSystMode; }
  DYTools::TCrossSectionKind_t csKind() const { return fCSKind; }
  void csKind(DYTools::TCrossSectionKind_t set_csKind) { fCSKind=set_csKind; }
  TString resNameBase() const { return fResNameBase; }

  void needsDetUnfolding(int yes) { fNeedsDETUnfolding=yes; }
  int  needsDetUnfolding() const { return fNeedsDETUnfolding; }
  void needsEffCorr(int yes) { fNeedsEffCorr=yes; }
  int  needsEffCorr() const { return fNeedsEffCorr; }
  void needsEffScaleCorr(int yes) { fNeedsEffScaleCorr=yes; }
  int  needsEffScaleCorr() const { return fNeedsEffScaleCorr; }
  void needsAccCorr(int yes) { fNeedsAccCorr=yes; }
  int  needsAccCorr() const { return fNeedsAccCorr; }
  void needsFsrCorr(int yes) { fNeedsFsrCorr=yes; }
  int  needsFsrCorr() const { return fNeedsFsrCorr; }

  void allNormErrorIsSyst(int yes) { fAllNormErrorIsSyst=yes; }
  int  allNormErrorIsSyst() const { return fAllNormErrorIsSyst; }
  void includeCorrError(int yes) { fIncludeCorrError=yes; }
  int  includeCorrError() const { return fIncludeCorrError; }

  void silentMode(int yes) { fSilentMode=yes; }
  int  silentMode() const { return fSilentMode; }
  void noSave(int yes) { fNoSave=yes; }
  int  noSave() const { return fNoSave; }

  friend std::ostream& operator<<(std::ostream &out, const InputArgs_t &ia) {
    out << "InputArgs(name=" << ia.fName << "):\n";
    out << "  inpMgr.loadedFileName=<" << ia.fInpMgr->loadedFileName() << ">\n";
    out << "  systMode=" << SystematicsStudyName(ia.fSystMode)
	<< ", csKind=" << CrossSectionKindName(ia.fCSKind) << "\n";
    out << "  resNameBase=" << ia.fResNameBase << "\n";
    out << "  - needsUnfolding=" << ia.fNeedsDETUnfolding << "\n";
    out << "  - needsEffCorr  =" << ia.fNeedsEffCorr << "\n";
    out << "  - needsRhoCorr  =" << ia.fNeedsEffScaleCorr << "\n";
    out << "  - needsAccCorr  =" << ia.fNeedsAccCorr << "\n";
    out << "  - needsFsrCorr  =" << ia.fNeedsFsrCorr << "\n";
    out << "  allNormErrorIsSyst=" << ia.fAllNormErrorIsSyst << "\n";
    out << "  includeCorrError=" << ia.fIncludeCorrError << "\n";
    out << "  silentMode=" << ia.fSilentMode << "\n";
    out << "  noSave    =" << ia.fNoSave << "\n";
    return out;
  }

};

// ----------------------------------------------
// ----------------------------------------------

struct CSResults_t {
  double fZXSec,fZXSecErr,fZXSecSystErr;
public:
  CSResults_t() : fZXSec(0.), fZXSecErr(0.), fZXSecSystErr(0.) {}
  
  CSResults_t(const CSResults_t &r) :
    fZXSec(r.fZXSec), fZXSecErr(r.fZXSecErr), fZXSecSystErr(r.fZXSecSystErr) {}

  void clear() {
    fZXSec=0; fZXSecErr=0.; fZXSecSystErr=0.;
  }

  double xs() const { return fZXSec; }
  double xsErr() const { return fZXSecErr; }
  double xsErrSyst() const { return fZXSecSystErr; }
  double xsSystErr() const { return fZXSecSystErr; }

  void assignCS(double set_xs, double set_xsErr, double set_xsSystErr) {
    fZXSec=set_xs; fZXSecErr=set_xsErr; fZXSecSystErr=set_xsSystErr;
  }

  void assignZpeakCS(const HistoPair2D_t &hp) {
    fZXSec = hp.ZpeakCount(&fZXSecErr);
    fZXSecSystErr= hp.ZpeakCountSystErr();
  }

  int doNormalize(HistoPair2D_t &hp, int allErrIsSyst) const {
    int res=1;
    if (allErrIsSyst) res=hp.divide_allErrSyst(fZXSec,fZXSecErr,fZXSecSystErr);
    else res=hp.divide(fZXSec,fZXSecErr,fZXSecSystErr);
    if (!res) std::cout << "error in CSResults::useToNormalize\n";
    return res;
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
int saveResult(const InputArgs_t &ia, const HistoPair2D_t &hp,
	       TString extraTag);



// cross section calculation
int calculateCSdistribution(const InputArgs_t &ia, const HistoPair2D_t &hp_ini, 
			    DYTools::TCrossSectionKind_t csKind,
			    HistoPair2D_t &hp_fin);
int calculateCS(const InputArgs_t &ia, const HistoPair2D_t &hp_ini, 
		DYTools::TCrossSectionKind_t csKind,
		HistoPair2D_t &hp_fin, CSResults_t &res);

int calcVecOfCSdistributions(const InputArgs_t &ia,
			     const std::vector<TH2D*> &yieldIniV,
			     DYTools::TCrossSectionKind_t csKind,
			     std::vector<TH2D*> &yieldFinV);

#endif
