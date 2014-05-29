#ifndef CSCovWorkFlags_HH
#define CSCovWorkFlags_HH

#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"

// ----------------------------------------------------------------
// ----------------------------------------------------------------

struct CSCovCalcFlags_t;  // what to calculate
struct WorkFlags_t; // auxiliary container of tags. Contains CSCovCalcFlags_t
struct TCovData_t; // the covariance matrices

// -----------------------------------------------------------
// -----------------------------------------------------------

typedef enum { _corrNone=0,
	       _yield,
	       _corrUnf,
	       _corrEff,
	       _corrESF,
	       _corrAcc,
	       _corrFSR,
	       _corrGlobalFSR, _corrGlobalPU,
	       _corrGlobalFEWZ,
	       _corrLast } TCorrCase_t;

// -----------------------------------------------------------
// -----------------------------------------------------------

TString corrCaseName(TCorrCase_t cs) {
  TString name;
  switch(cs) {
  case _corrNone: name="NONE"; break;
  case _yield: name="yield"; break;
  case _corrUnf: name="corrUnf"; break;
  case _corrEff: name="corrEff"; break;
  case _corrESF: name="corrESF"; break;
  case _corrAcc: name="corrAcc"; break;
  case _corrFSR: name="corrFSR"; break;
  case _corrGlobalFSR: name="corrGlobalFSR"; break;
  case _corrGlobalPU: name="corrGlobalPU"; break;
  case _corrGlobalFEWZ: name="corrGlobalFEWZ"; break;
  case _corrLast: name="LAST"; break;
  default: name="UNKNOWN";
  }
  return name;
}

// ----------------------------------------------------------------
// ----------------------------------------------------------------

//=== Functions =================================================================================================

int loadYieldCovMatrices(const TString &fnameBase,std::vector<TMatrixD*> &covs,
			 std::vector<TString> &labels, const WorkFlags_t &wf);
int loadUnfCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs,
		       std::vector<TString> &labels, const WorkFlags_t &wf);
int loadEffCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs,
		       std::vector<TString> &labels, const WorkFlags_t &wf);
int loadEsfCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs,
		       std::vector<TString> &labels, const WorkFlags_t &wf);
int loadAccCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs,
		       std::vector<TString> &labels, const WorkFlags_t &wf);
int loadFsrCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs,
		       std::vector<TString> &labels, const WorkFlags_t &wf);

int loadGlobalCovMatrices(const TString &fnameBase,
			  std::vector<TMatrixD*> &covs,
			  std::vector<TString> &labels,
			  const WorkFlags_t &wf);

TH2D *loadMainCSResult(int crossSection=0, TH2D** h2SystErr=NULL);

// ----------------------------------------------------------------
// ----------------------------------------------------------------


struct CSCovCalcFlags_t {
  int fYieldStat, fYieldSyst; // ignore the way stat error was obtained
  int fYieldStatDetailed; // take only observed yield stat error
  int fYieldSystDetailed; // all background error is syst (fake/true)
  int fYieldEScale;
  int fUnfPU, fUnfFSR, fUnfRnd, fUnfEScale;
  int fEffPU, fEffFSR, fEffRnd;
  int fESFtot;
  int fESFtotCheck;
  int fAccPU, fAccFSR, fAccRnd;
  int fFsrPU, fFsrFSR, fFsrRnd;
  int fGlobalFSR, fGlobalPU, fGlobalFEWZ;
public:
  CSCovCalcFlags_t() :
    fYieldStat(0), fYieldSyst(0),
    fYieldStatDetailed(0), fYieldSystDetailed(0),
    fYieldEScale(0),
    fUnfPU(0), fUnfFSR(0), fUnfRnd(0), fUnfEScale(0),
    fEffPU(0), fEffFSR(0), fEffRnd(0),
    fESFtot(0), fESFtotCheck(0),
    fAccPU(0), fAccFSR(0), fAccRnd(0),
    fFsrPU(0), fFsrFSR(0), fFsrRnd(0),
    fGlobalFSR(0), fGlobalPU(0), fGlobalFEWZ(0)
  {}


  CSCovCalcFlags_t(const CSCovCalcFlags_t &flags) :
    fYieldStat(0), fYieldSyst(0),
    fYieldStatDetailed(0), fYieldSystDetailed(0),
    fYieldEScale(0),
    fUnfPU(0), fUnfFSR(0), fUnfRnd(0), fUnfEScale(0),
    fEffPU(0), fEffFSR(0), fEffRnd(0),
    fESFtot(0), fESFtotCheck(0),
    fAccPU(0), fAccFSR(0), fAccRnd(0),
    fFsrPU(0), fFsrFSR(0), fFsrRnd(0),
    fGlobalFSR(0), fGlobalPU(0), fGlobalFEWZ(0)
  {
    this->assign(flags);
  }

  int calc_YieldStat() const { return fYieldStat; }
  void calc_YieldStat(int yes) { fYieldStat=yes; }
  int calc_YieldSyst() const { return fYieldSyst; }
  void calc_YieldSyst(int yes) { fYieldSyst=yes; }
  int calc_YieldStatDetailed() const { return fYieldStatDetailed; }
  void calc_YieldStatDetailed(int yes) { fYieldStatDetailed=yes; }
  int calc_YieldSystDetailed() const { return fYieldSystDetailed; }
  void calc_YieldSystDetailed(int yes) { fYieldSystDetailed=yes; }
  int calc_YieldEscale() const { return fYieldEScale; }
  void calc_YieldEscale(int yes) { fYieldEScale=yes; }
  int calc_UnfPU() const { return fUnfPU; }
  void calc_UnfPU(int yes) { fUnfPU=yes; }
  int calc_UnfFSR() const { return fUnfFSR; }
  void calc_UnfFSR(int yes) { fUnfFSR=yes; }
  int calc_UnfRnd() const { return fUnfRnd; }
  void calc_UnfRnd(int yes) { fUnfRnd=yes; }
  int calc_UnfEScale() const { return fUnfEScale; }
  void calc_UnfEScale(int yes) { fUnfEScale=yes; }
  int calc_EffPU() const { return fEffPU; }
  void calc_EffPU(int yes) { fEffPU=yes; }
  int calc_EffFSR() const { return fEffFSR; }
  void calc_EffFSR(int yes) { fEffFSR=yes; }
  int calc_EffRnd() const { return fEffRnd; }
  void calc_EffRnd(int yes) { fEffRnd=yes; }
  int calc_ESFtot() const { return fESFtot; }
  void calc_ESFtot(int yes) { fESFtot=yes; }
  int calc_ESFtotCheck() const { return fESFtotCheck; }
  void calc_ESFtotCheck(int yes) { fESFtotCheck=yes; }
  int calc_AccPU() const { return fAccPU; }
  void calc_AccPU(int yes) { fAccPU=yes; }
  int calc_AccFSR() const { return fAccFSR; }
  void calc_AccFSR(int yes) { fAccFSR=yes; }
  int calc_AccRnd() const { return fAccRnd; }
  void calc_AccRnd(int yes) { fAccRnd=yes; }
  int calc_FsrPU() const { return fFsrPU; }
  void calc_FsrPU(int yes) { fFsrPU=yes; }
  int calc_FsrFSR() const { return fFsrFSR; }
  void calc_FsrFSR(int yes) { fFsrFSR=yes; }
  int calc_FsrRnd() const { return fFsrRnd; }
  void calc_FsrRnd(int yes) { fFsrRnd=yes; }
  int calc_globalFSR() const { return fGlobalFSR; }
  void calc_globalFSR(int yes) { fGlobalFSR=yes; }
  int calc_globalPU() const { return fGlobalPU; }
  void calc_globalPU(int yes) { fGlobalPU=yes; }
  int calc_globalFEWZ() const { return fGlobalFEWZ; }
  void calc_globalFEWZ(int yes) { fGlobalFEWZ=yes; }

  int doCalcYieldCov() const {
    return (fYieldStat + fYieldStatDetailed +
	    fYieldSyst + fYieldSystDetailed +
	    fYieldEScale);
  }

  int doCalcUnfCov() const {
    return (fUnfPU + fUnfFSR + fUnfRnd + fUnfEScale);
  }

  int doCalcEffCov() const {
    return (fEffPU + fEffFSR + fEffRnd);
  }

  int doCalcESFCov() const {
    return (fESFtot + fESFtotCheck);
  }

  int doCalcAccCov() const {
    return (fAccPU + fAccFSR + fAccRnd);
  }

  int doCalcFSRCov() const {
    return (fFsrPU + fFsrFSR + fFsrRnd);
  }

  int doCalcGlobalCov() const {
    return (fGlobalFSR + fGlobalPU + fGlobalFEWZ);
  }

  // ---------------------------

  int finalizeFlags();

  // ---------------------------

  void assign(const CSCovCalcFlags_t &f) {
    if (this == &f) return;
    fYieldStat=f.fYieldStat;
    fYieldSyst=f.fYieldSyst;
    fYieldStatDetailed=f.fYieldStatDetailed;
    fYieldSystDetailed=f.fYieldSystDetailed;
    fYieldEScale=f.fYieldEScale;
    fUnfPU=f.fUnfPU;
    fUnfFSR=f.fUnfFSR;
    fUnfRnd=f.fUnfRnd;
    fUnfEScale=f.fUnfEScale;
    fEffPU=f.fEffPU;
    fEffFSR=f.fEffFSR;
    fEffRnd=f.fEffRnd;
    fESFtot=f.fESFtot;
    fESFtotCheck=f.fESFtotCheck;
    fAccPU =f.fAccPU;
    fAccFSR=f.fAccFSR;
    fAccRnd=f.fAccRnd;
    fFsrPU =f.fFsrPU;
    fFsrFSR=f.fFsrFSR;
    fFsrRnd=f.fFsrRnd;
    fGlobalFSR=f.fGlobalFSR;
    fGlobalPU=f.fGlobalPU;
    fGlobalFEWZ=f.fGlobalFEWZ;
  }
};



// ----------------------------------------------------------------
// ----------------------------------------------------------------


struct WorkFlags_t {
  int fCase;
  int fCSCov;
  int fSaveTotCovDetails;
  TString fExtraTag;
  std::vector<TString> fExtraTagV;
  CSCovCalcFlags_t fCalcFlags;
public:
  WorkFlags_t(int the_case=0, int set_showCSCov=1, TString set_extra_tag="") :
    fCase(the_case), fCSCov(set_showCSCov),
    fSaveTotCovDetails(0),
    fExtraTag(set_extra_tag),
    fExtraTagV(),
    fCalcFlags()
  {
    init_ExtraTagV(1);
  }

  WorkFlags_t(const WorkFlags_t &w) :
    fCase(w.fCase), fCSCov(w.fCSCov),
    fSaveTotCovDetails(w.fSaveTotCovDetails),
    fExtraTag(w.fExtraTag),
    fExtraTagV(w.fExtraTagV),
    fCalcFlags(w.fCalcFlags)
  {}

  int theCase() const { return fCase; }
  void theCase(int the_case) { fCase=the_case; }
  int showCSCov() const { return fCSCov; }
  void showCSCov(int show) { fCSCov=show; }
  int saveTotCovDetails() const { return fSaveTotCovDetails; }
  void saveTotCovDetails(int yes) { fSaveTotCovDetails=yes; }
  int hasExtraTag() const { return (fExtraTag.Length()>0) ? 1:0; }
  TString extraFileTag() const { return fExtraTag; }
  void extraFileTag(TString setTag) { fExtraTag=setTag; }

  const CSCovCalcFlags_t & calcFlags() const { return fCalcFlags; }
  const CSCovCalcFlags_t & cf() const { return fCalcFlags; }
  CSCovCalcFlags_t & editCalcFlags() { return fCalcFlags; }
  //int doCalcYieldCov() const { return fCalcFlags.doCalcYieldCof; }

  int finalizeFlags() { return fCalcFlags.finalizeFlags(); }

  TString extraFileTag(TCorrCase_t idx1) const {
    int idx=int(idx1);
    return fExtraTagV[idx];
  }

  void extraFileTag(TCorrCase_t idx1, TString tag) {
    int idx=int(idx1);
    fExtraTagV[idx]=tag;
  }

  // the user should call with only_global=0, if different extraTags are needed
  void init_ExtraTagV(int only_global) {
    if (fExtraTagV.size()>0) fExtraTagV.clear();
    if (!only_global) {
      fExtraTagV.reserve(int(_corrLast));
      for (int idx=0; idx<int(_corrLast); ++idx) {
	fExtraTagV.push_back(TString());
      }
    }
  }

  TString fieldName(TString tag) const {
    TString field=(fCSCov) ? "covCS_" : "cov_";
    field.Append(tag);
    return field;
  }

  void adjustFName(TString &fname) const {
    if (fExtraTag.Length()) {
      if (fname.Index(".root")>0) {
	fname.ReplaceAll(".root",fExtraTag + TString(".root"));
      }
      else {
	fname.Append(fExtraTag);
      }
      std::cout << "WorkFlags_t::adjustFName: fname=<" << fname << ">\n";
    }
  }

  void adjustFName(TString &fname, TCorrCase_t cs) const {
    int idx=int(cs);
    TString eTag;
    std::cout << "idx=" << idx << ", fExtraTagV.size=" << fExtraTagV.size() << "\n";
    if (idx<int(fExtraTagV.size())) eTag=fExtraTagV[idx];
    else eTag=fExtraTag;
    std::cout << "eTag=<" << eTag << ">\n";
    if (eTag.Length()) {
      if (fname.Index(".root")>0) {
	fname.ReplaceAll(".root",eTag + TString(".root"));
      }
      else {
	fname.Append(eTag);
      }
      std::cout << "WorkFlags_t::adjustFName(corrCase="
		<< corrCaseName(cs) << "): fname=<" << fname << ">\n";
    }
  }


};

// ----------------------------------------------------------------
// ----------------------------------------------------------------

struct TCovData_t {
  std::vector<int> isActive;
  std::vector<TMatrixD*> covYieldV, covUnfV, covEffV, covEsfV;
  std::vector<TMatrixD*> covAccV, covFsrV, covGlobalV;
  std::vector<TString> labelYieldV, labelUnfV, labelEffV, labelEsfV;
  std::vector<TString> labelAccV, labelFsrV, labelGlobalV;
public:

  // ----------------

  TCovData_t() : isActive(0),
		 covYieldV(), covUnfV(), covEffV(), covEsfV(),
		 covAccV(), covFsrV(), covGlobalV(),
		 labelYieldV(), labelUnfV(), labelEffV(), labelEsfV(),
		 labelAccV(), labelFsrV(), labelGlobalV()
  {
    isActive.reserve(7);
    for (int i=0; i<7; ++i) isActive.push_back(0);
  }


  // ----------------

  const std::vector<TMatrixD*>* getCovV(int i) const {
    const std::vector<TMatrixD*>* ptr=NULL;
    switch(i) {
    case 0: ptr=&covYieldV; break;
    case 1: ptr=&covUnfV; break;
    case 2: ptr=&covEffV; break;
    case 3: ptr=&covEsfV; break;
    case 4: ptr=&covAccV; break;
    case 5: ptr=&covFsrV; break;
    case 6: ptr=&covGlobalV; break;
    default:
      std::cout << "index error in getCovV\n";
    }
    return ptr;
  }

  // ----------------

  const std::vector<TString>* getLabelV(int i) const {
    const std::vector<TString>* ptr=NULL;
    switch(i) {
    case 0: ptr=&labelYieldV; break;
    case 1: ptr=&labelUnfV; break;
    case 2: ptr=&labelEffV; break;
    case 3: ptr=&labelEsfV; break;
    case 4: ptr=&labelAccV; break;
    case 5: ptr=&labelFsrV; break;
    case 6: ptr=&labelGlobalV; break;
    default:
      std::cout << "index error in getLabelV\n";
    }
    return ptr;
  }

  // ----------------

  TMatrixD* calcTotalCov() const {
    TMatrixD* totcov=NULL;
    for (unsigned int idx=0; idx<isActive.size(); ++idx) {
      if (!isActive[idx]) continue;
      const std::vector<TMatrixD*> *covs=this->getCovV(idx);
      for (unsigned int i=0; i<covs->size(); ++i) {
	TMatrixD *cov=(*covs)[i];
	if (totcov==NULL) totcov = new TMatrixD(*cov);
	else (*totcov) += (*cov);
      }
    }
    return totcov;
  }

  // ----------------

  // increase the error on the efficiency scale factors
  int addESFsyst(const TString ver="20140525");

  // ----------------

  int Write(TString subDir) const;

  // ----------------

};


// ----------------------------------------------------------------
// ----------------------------------------------------------------

struct NormCS_t {
  DYTools::TCrossSectionKind_t fcsKind;
  double fcs,fcsErrNoLumi;
public:
  NormCS_t(DYTools::TCrossSectionKind_t set_csKind=DYTools::_cs_None,
	   const TString fname="dir-CSInfo/inpCS_sigmaZ.dat");

  DYTools::TCrossSectionKind_t csKind() const { return fcsKind; }
  double cs() const { return fcs; }
  double csErrNoLumi() const { return fcsErrNoLumi; }

  int loadValues(const TString fname);

  int loadValues(const TString fname,
		 DYTools::TCrossSectionKind_t set_csKind) {
    fcsKind=set_csKind;
    return loadValues(fname);
  }

  friend
  std::ostream& operator<<(std::ostream& out, const NormCS_t &zcs) {
    out << "NormCS(" << zcs.fcsKind << ", " << zcs.fcs << ", "
	<< zcs.fcsErrNoLumi << ")";
    return out;
  }
};

// ----------------------------------------------------------------
// ----------------------------------------------------------------


#endif

