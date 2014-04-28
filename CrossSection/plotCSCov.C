#include "../Include/DYTools.hh"
#include "../CrossSection/crossSectionFnc.hh"
#include "../Include/colorPalettes.hh"
#include <TBenchmark.h>

//=== Global flags =================================================================================================

int calc_YieldStat=1;
int calc_YieldSyst=1;
int calc_YieldUnregEn=0; // EScale includes all escale uncertainty
int calc_YieldEScale=0;
int calc_YieldApplyEScale=0; // EScale include all escale uncertainty

int doCalcYieldCov=(calc_YieldStat + calc_YieldSyst + calc_YieldUnregEn +
		    calc_YieldEScale + calc_YieldApplyEScale) ? 1:0;

int calc_UnfPU=0;
int calc_UnfFSR=0;
int calc_UnfRnd=1;
int calc_UnfEScale=1;

int doCalcUnfCov=(calc_UnfPU + calc_UnfFSR +
		  calc_UnfRnd + calc_UnfEScale) ? 1:0;

int calc_EffPU=0;
int calc_EffFSR=0;
int calc_EffRnd=1;

int doCalcEffCov=(calc_EffPU + calc_EffFSR + calc_EffRnd) ? 1:0;

int calc_ESFtot=1;
int calc_ESFtotCheck=0;

int doCalcESFCov=(calc_ESFtot+calc_ESFtotCheck) ? 1:0;

// acceptance correction is later adjusted once 1D/2D measurement is known
int calc_AccFSR=0;
int calc_AccRnd=1;

int doCalcAccCov=(calc_AccFSR + calc_AccRnd) ? 1:0;

int calc_FsrFSR=0; // not ready!
int calc_FsrRnd=1;

int doCalcFSRCov=(calc_FsrFSR + calc_FsrRnd) ? 1:0;

int calc_globalFSR=1;
int calc_globalPU=1;
int calc_globalFEWZ=0;

int doCalcGlobalCov=(calc_globalFSR + calc_globalPU + calc_globalFEWZ) ? 1:0;

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

// -----------------------------------------------------------

struct WorkFlags_t {
  int fCase;
  int fCSCov;
  TString fExtraTag;
  std::vector<TString> fExtraTagV;
public:
  WorkFlags_t(int the_case=0, int set_showCSCov=1, TString set_extra_tag="") :
    fCase(the_case), fCSCov(set_showCSCov),
    fExtraTag(set_extra_tag),
    fExtraTagV()
  {
    init_ExtraTagV(1);
  }

  WorkFlags_t(const WorkFlags_t &w) :
    fCase(w.fCase), fCSCov(w.fCSCov),
    fExtraTag(w.fExtraTag),
    fExtraTagV(w.fExtraTagV)
  {}

  int theCase() const { return fCase; }
  void theCase(int the_case) { fCase=the_case; }
  int showCSCov() const { return fCSCov; }
  void showCSCov(int show) { fCSCov=show; }
  int hasExtraTag() const { return (fExtraTag.Length()>0) ? 1:0; }
  TString extraFileTag() const { return fExtraTag; }
  void extraFileTag(TString setTag) { fExtraTag=setTag; }

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

// -----------------------------------------------------------

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

};


//=== Functions =================================================================================================

int loadYieldCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf);
int loadUnfCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf);
int loadEffCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf);
int loadEsfCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf);
int loadAccCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf);
int loadFsrCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf);
int loadGlobalCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf);

TH2D *loadMainCSResult();
int workWithData(TCovData_t &dt, const WorkFlags_t &wf);


//=== MAIN MACRO =================================================================================================


int plotCSCov(int analysisIs2D, TString conf, int the_case, int showCSCov=1,
	      TString outFileExtraTag_UserInput="")
{

  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return retCodeError;
  }

  calc_AccFSR*=(1-analysisIs2D);
  calc_AccRnd*=(1-analysisIs2D);
  doCalcAccCov*=(1-analysisIs2D);

  // Settings 
  //==============================================================================================================

  DYTools::TCrossSectionKind_t csKind=(DYTools::study2D) ? DYTools::_cs_preFsrDet : DYTools::_cs_preFsr;

  DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN;
  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;

  TCovData_t dt;
  WorkFlags_t work(the_case,showCSCov,outFileExtraTag_UserInput);

  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(conf)) return retCodeError;

  // Construct eventSelector, update mgr and plot directory
  TString extraTag; // empty
  TString plotExtraTag;

  EventSelector_t evtSelector(inpMgr,runMode,systMode,
			      extraTag,plotExtraTag,
			      EventSelector::_selectDefault);

  int systFileFlag=1;
  TString fnameBase=inpMgr.crossSectionFullFileName(systMode,
					       csKind,0,systFileFlag);
  std::cout << "fnameBase=<" << fnameBase << ">\n";

  if ((the_case==2) || (the_case==3) || (the_case==4)) {
    // global study I
    work.init_ExtraTagV(0);
    work.extraFileTag(_yield, "-yieldOnly_nExps1000");
    work.extraFileTag(_corrUnf, "-unfOnly");
    work.extraFileTag(_corrEff, "-effRndOnly");
    work.extraFileTag(_corrESF, "-esfOnly");
    work.extraFileTag(_corrAcc, "-accRndOnly");
    work.extraFileTag(_corrFSR, "-fsrRndOnly");
    if (the_case==3) work.extraFileTag(_corrFSR, "-fsrRndOnly_nExps1000");
    //work.extraFileTag(_corrGlobalFSR, "-globalFSROnly_nExps20");
    //work.extraFileTag(_corrGlobalPU, "-globalPUOnly_nExps20");
    work.extraFileTag(_corrGlobalFSR, "-globalFSROnly");
    work.extraFileTag(_corrGlobalPU, "-globalPUOnly");
    if (the_case==4) {
      std::cout << "\n\n\t EXCLUDING RNDs beyond rho\n";
      calc_EffRnd=0;
      calc_AccRnd=0;
      calc_FsrRnd=0;
    }
  }

  int res=1;
  if (res && doCalcYieldCov) { 
    if (!loadYieldCovMatrices(fnameBase,dt.covYieldV,dt.labelYieldV,work)) return 0;
    dt.isActive[0]=1;
  }
  if (res && doCalcUnfCov) { 
    if (!loadUnfCovMatrices(fnameBase,dt.covUnfV,dt.labelUnfV,work)) return 0;
    dt.isActive[1]=1;
  }
  if (res && doCalcEffCov) { 
    if (!loadEffCovMatrices(fnameBase,dt.covEffV,dt.labelEffV,work)) return 0;
    dt.isActive[2]=1;
  }
  if (res && doCalcESFCov) {
    if (!loadEsfCovMatrices(fnameBase,dt.covEsfV,dt.labelEsfV,work)) return 0;
    dt.isActive[3]=1;
  }
  if (res && doCalcAccCov) { 
    if (!loadAccCovMatrices(fnameBase,dt.covAccV,dt.labelAccV,work)) return 0;
    dt.isActive[4]=1;
  }
  if (res && doCalcFSRCov) { 
    if (!loadFsrCovMatrices(fnameBase,dt.covFsrV,dt.labelFsrV,work)) return 0;
    dt.isActive[5]=1;
  }

  if (res && doCalcGlobalCov) {
    if (calc_globalFSR &&
	(calc_FsrFSR || calc_AccFSR || calc_EffFSR || calc_UnfFSR)) {
      std::cout << "since calc_globalFSR is on, "
		<< "individual FSR studies have to be switched off\n";
      return 0;
    }
    if (calc_globalPU && (calc_EffPU || calc_UnfPU)) {
      std::cout << "since calc_globalPU is on, "
		<< "individual PU studies have to be switched off\n";
      return 0;
    }
    if(!loadGlobalCovMatrices(fnameBase,dt.covGlobalV,dt.labelGlobalV,work))
      return 0;
    dt.isActive[6]=1;
  }

  if (!workWithData(dt,work)) return 0;

  return retCodeOk;
}

// ---------------------------------------------------------------------------
  // Implementations
  //==============================================================================================================


int loadYieldCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf) {

  TString inpFileName=fnameBase;
  //if (wf.extraFileTag().Index("-yield")==-1) {
  //  inpFileName.ReplaceAll(".root","-yieldOnly.root");
  //}
  wf.adjustFName(inpFileName,_yield);
  TFile fin(inpFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << inpFileName << ">\n";
    return 0;
  }
  TMatrixD *ptr;
  if (calc_YieldStat) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("YieldStat"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("signal stat");
    }
  }
  if (calc_YieldSyst) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("YieldSyst"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("signal syst");
    }
  }
  if (calc_YieldUnregEn) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("YieldUnregEn"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("signal unreg.en.");
    }
  }
  if (calc_YieldEScale) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("YieldEScale"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("signal EScale uncert.");
    }
  }
  if (calc_YieldApplyEScale) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("YieldApplyEScale"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("signal adhoc escale");
      }
  }
  std::cout << "loaded " << covs.size() << " entries from file <" << fin.GetName() << ">\n";
  fin.Close();
  return 1;
}


// -----------------------------------------------------------


int loadUnfCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf) {

  TString inpFileName=fnameBase;
  //if (wf.extraFileTag().Index("-unf")==-1) {
  //  inpFileName.ReplaceAll(".root","-unfOnly.root");
  //}
  wf.adjustFName(inpFileName,_corrUnf);
  TFile fin(inpFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << inpFileName << ">\n";
    return 0;
  }
  TMatrixD *ptr;
  if (calc_UnfPU) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("UnfPU"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("unf pile-up");
    }
  }
  if (calc_UnfFSR) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("UnfFSR"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("unf FSR");
    }
  }
  if (calc_UnfRnd) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("UnfRnd"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("unf stat");
    }
  }
  if (calc_UnfEScale) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("UnfEScale"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("unf e-scale");
    }
  }

  std::cout << "loaded " << covs.size() << " entries from file <" << fin.GetName() << ">\n";
  fin.Close();
  return 1;
}


// -----------------------------------------------------------


int loadEffCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf) {

  TString inpFileName=fnameBase;
  //if (wf.extraFileTag().Index("-eff")==-1) {
  //  inpFileName.ReplaceAll(".root","-effOnly.root");
  //}
  wf.adjustFName(inpFileName,_corrEff);
  TFile fin(inpFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << inpFileName << ">\n";
    return 0;
  }
  TMatrixD *ptr;
  if (calc_EffPU) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("EffPU"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("eff pile-up");
    }
  }
  if (calc_EffFSR) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("EffFSR"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("eff FSR");
    }
  }
  if (calc_EffRnd) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("EffRnd"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("eff stat");
    }
  }

  std::cout << "loaded " << covs.size() << " entries from file <" << fin.GetName() << ">\n";
  fin.Close();
  return 1;
}


// -----------------------------------------------------------


int loadEsfCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf) {

  TString inpFileName=fnameBase;
  //if (wf.extraFileTag().Index("-esf")==-1) {
  //  inpFileName.ReplaceAll(".root","-esfOnly.root");
  //}
  wf.adjustFName(inpFileName,_corrESF);
  TFile fin(inpFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << inpFileName << ">\n";
    return 0;
  }
  TMatrixD *ptr;
  if (calc_ESFtot) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("ESFtot"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("ESF tot");
    }
  }
  else if (calc_ESFtotCheck) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("ESFtotCheck"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("ESF tot");
    }
  }

  std::cout << "loaded " << covs.size() << " entries from file <" << fin.GetName() << ">\n";
  fin.Close();
  return 1;
}



// -----------------------------------------------------------


int loadAccCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf) {

  TString inpFileName=fnameBase;
  //inpFileName.ReplaceAll(".root","-accOnly.root");
  wf.adjustFName(inpFileName,_corrAcc);
  TFile fin(inpFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << inpFileName << ">\n";
    return 0;
  }
  TMatrixD *ptr;
  if (calc_AccFSR) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("AccFSR"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("acc FSR");
    }
  }
  if (calc_AccRnd) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("AccRnd"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("acc stat");
    }
  }

  std::cout << "loaded " << covs.size() << " entries from file <" << fin.GetName() << ">\n";
  fin.Close();
  return 1;
}



// -----------------------------------------------------------


int loadFsrCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf) {

  TString inpFileName=fnameBase;
  //inpFileName.ReplaceAll(".root","-fsrUnfOnly.root");
  wf.adjustFName(inpFileName,_corrFSR);
  TFile fin(inpFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << inpFileName << ">\n";
    return 0;
  }
  TMatrixD *ptr;
  if (calc_FsrFSR) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("FsrFSR"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("FSR FSR");
    }
  }
  if (calc_FsrRnd) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("FsrRnd"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("FSR stat");
    }
  }

  std::cout << "loaded " << covs.size() << " entries from file <" << fin.GetName() << ">\n";
  fin.Close();
  return 1;
}



// -----------------------------------------------------------


int loadGlobalCovMatrices(const TString &fnameBase,
			  std::vector<TMatrixD*> &covs,
			  std::vector<TString> &labels,
			  const WorkFlags_t &wf) {

  if (fnameBase.Length()==0) {
    std::cout << "loadGlobalCovMatrices warning: fnameBase.Length=0\n";
  }

  if (wf.showCSCov()==0) {
    std::cout << "loadGlobalCovMatrices: results are available only for "
	      << "the final cross section\n";
    return 0;
  }

  for (int i=0; i<3; ++i) {
    TString tag;
    int calc=0;
    TCorrCase_t corrCase=_corrNone;
    switch(i) {
    case 0: tag="puRndStudy"; calc=calc_globalPU;
            corrCase=_corrGlobalPU;
	    break;
    case 1: tag="fsrRndStudy"; calc=calc_globalFSR;
            corrCase=_corrGlobalFSR;
	    break;
    case 2: tag="fewzRndStudy"; calc=calc_globalFEWZ;
            corrCase=_corrGlobalFEWZ;
	    break;
    default:
      std::cout << "loadGlobalCovMatrices: the case i=" << i
		<< " is not ready\n";
      return 0;
    }
    if (!calc) continue;

    TString inpFileName;
    if (0) {
      // local file
      inpFileName=Form("csSyst-%s-%dD.root",tag.Data(),
		       DYTools::study2D+1);
    }
    else {
      inpFileName=fnameBase;
      wf.adjustFName(inpFileName,corrCase);
    }

    TFile fin(inpFileName,"read");
    if (!fin.IsOpen()) {
      std::cout << "failed to open a file <" << inpFileName << ">\n";
      return 0;
    }
    TString field= TString("covCS_") + tag;
    std::cout << "loading <" << field << "> from <" << inpFileName << ">\n";
    TMatrixD *ptr = (TMatrixD*)fin.Get(field);
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back(tag);
    }
    else {
      std::cout << "failed to get <" << field << "> from <"
		<< fin.GetName() << ">\n";
      return 0;
    }
    fin.Close();
  }

  std::cout << "loaded " << covs.size() << " entries from global files\n";
  return 1;
}


// -----------------------------------------------------------
// -----------------------------------------------------------

TH2D *loadMainCSResult() {
  TString csFileName="../../Results-DYee/root_files_reg/xsec/DY_j22_19712pb/xSec_preFsr_1DpreFsrFullSp.root";
  TString fieldName="hpPreFsrFullSp";
  if (DYTools::study2D) {
    csFileName="../../Results-DYee/root_files_reg/xsec/DY_j22_19712pb/xSec_preFsrDet_2DpreFsrDet.root";
    fieldName="hpPreFsrDet";
  }
  TFile fin(csFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open the cross-section file <"
	      << fin.GetName() << ">\n";
    return NULL;
  }
  TH2D *h2=LoadHisto2D(fin,fieldName,"",1);
  fin.Close();
  if (!h2) {
    std::cout << "loadMainCSResult error\n";
  }
  return h2;
}

// -----------------------------------------------------------
// -----------------------------------------------------------

void plotAllCovs(TCovData_t &dt, const WorkFlags_t &wf) {
  //gStyle->SetPalette(1);

  TMatrixD *totalCov=dt.calcTotalCov();

  // plot covariances, correlations, and partial correlations
  if (1)
  for (int iCorr=0; iCorr<3; ++iCorr) {
    //if (iCorr!=2) continue;
    TString covStr;
    switch(iCorr) {
    case 0: covStr="Cov_"; break;
    case 1: covStr="Corr_"; break;
    case 2: covStr="partCorr_"; break;
    default:
      std::cout << "plotAllCovs: unknown iCorr=" << iCorr << "\n";
      return;
    }
    if (wf.showCSCov()) covStr.Prepend("CS");

    for (unsigned int idx=0; idx<dt.isActive.size(); ++idx) {
      //if (idx!=3) continue;

      if (!dt.isActive[idx]) continue;
      const std::vector<TMatrixD*> * covV= dt.getCovV(idx);
      const std::vector<TString>* labelV= dt.getLabelV(idx);
      for (unsigned int i=0; i<covV->size(); ++i) {
	TString cName=TString("canv") + covStr + (*labelV)[i];
	TCanvas *cx= new TCanvas(cName,cName, 750,700);
	AdjustFor2DplotWithHeight(cx);
	TString explain;
	TString histoTag=covStr + (*labelV)[i];
	if (wf.hasExtraTag()) {
	  histoTag.Append(" ");
	  wf.adjustFName(histoTag);
	}
	eliminateSeparationSigns(histoTag);
	TH2D* h2raw=NULL;
	TString histoName;
	TString histoTitle;

	if (iCorr==0) {
	  //set_nice_style(51);
	  gStyle->SetPalette(1);
	  //(*covV)[i]->Draw("COLZ");
	  h2raw=createHisto2D(*(*covV)[i], NULL,
			      "h2Covariance_base",
			      "h2Covariance_base",
			      _colrange_none,0,1.0001);
	  //explain="Covariance ";
	  //wf.adjustFName(explain);
	  histoName="h2Covariance";
	  histoTitle="Covariance ";
	}
	else if (iCorr==1) {
	  gStyle->SetPalette(1);
	  set_center_white_style5(11);
	  TMatrixD *corr= corrFromCov( *(*covV)[i] );
	  TString histoNameBase=TString("hCorr_base") + histoTag;
	  histoName="h2Corr";
	  histoTitle=TString("Correlations ");
	  h2raw=createHisto2D(*corr,NULL,histoNameBase,histoTitle,_colrange_center,0,1.0001);

	  delete corr;
	  //explain="Correlation ";
	}
	else if (iCorr==2) {
	  set_center_white_style(21);
	  //if (DYTools::study2D)
	  //else
	  set_center_white_style5(11);
	  //gStyle->SetPalette(1);
	  TMatrixD *corr= partialCorrFromCov( *totalCov, *(*covV)[i] );
	  TString histoNameBase=TString("hPartCorrBase") + histoTag;
	  histoName=TString("hPartCorr");
	  histoTitle=TString("Partial correlations ");
	  h2raw=createHisto2D(*corr,NULL,histoNameBase,histoTitle,_colrange_center,0,1.0001);
	  delete corr;
	  //explain="Correlation ";
	}

	histoName.Append(histoTag);
	histoTitle.Append(histoTag);
	TH2D* h2save=clipToAnalysisUnfBins(h2raw,histoName,histoTitle,1); // reset axis!
	if (iCorr) h2save->GetZaxis()->SetRangeUser(-1.0001,1.0001);
	h2save->Draw("COLZ");

	if (explain.Length()) {
	  TText *txt=new TText();
	  txt->DrawTextNDC(0.25,0.93,explain + (*labelV)[i]);
	}
	cx->Update();

	if (1) {
	  TString figName=TString("fig-") + DYTools::analysisTag + TString("--") + histoTag;
	  //eliminateSeparationSigns(figName);
	  std::cout << "figName=<" << figName << ">\n";
	  SaveCanvas(cx,figName);
	}
      }
    }
  }

  if (wf.showCSCov()) {
    // plot error profile
    for (int iCorr=0; iCorr<2; ++iCorr) {
      TString covStr;
      TH2D* h2Main=NULL;
      TString yAxisLabel="uncertainty from cov";
      switch(iCorr) {
      case 0:  covStr="CSCov_"; break;
      case 1:
	covStr="CSCov_";
	h2Main=loadMainCSResult();
	if (!h2Main) return;
	removeError(h2Main);
	h2Main->Scale(0.01);
	yAxisLabel="relative uncertainty from cov (%)";
	break;
      default:
	std::cout << "plotAllCovs: unknown iCorr=" << iCorr << " /2nd loop/\n";
	return;
      }

      std::vector<TH2D*> errFromCovV;
      std::vector<TString> errFromCovLabelV;

      for (unsigned int idx=0; idx<dt.isActive.size(); ++idx) {
	if (!dt.isActive[idx]) continue;
	const std::vector<TMatrixD*> * covV= dt.getCovV(idx);
	const std::vector<TString>* labelV= dt.getLabelV(idx);

	errFromCovV.reserve(errFromCovV.size() + covV->size());
	errFromCovLabelV.reserve(errFromCovLabelV.size() + covV->size());

	for (unsigned int i=0; i<covV->size(); ++i) {
	  TString histoLabel=TString("histoErr_") + (*labelV)[i];
	  TH2D* h2=errorFromCov(*(*covV)[i],histoLabel);
	  if (!h2) {
	    std::cout << "failed to create the error histogram "
		      << histoLabel << "\n";
	    return;
	  }
	  if (iCorr==1) {
	    if (!scaleHisto(h2,h2Main)) return;
	  }
	  errFromCovV.push_back(h2);
	  errFromCovLabelV.push_back((*labelV)[i]);
	}
      }

      TString canvName=Form("canvErr_%d",iCorr);
      std::vector<std::vector<TH1D*>*> hProfV;
      std::vector<ComparisonPlot_t*> cpV;
      int delayDraw=0;

      TCanvas *cx=plotProfiles(canvName,
			       errFromCovV, errFromCovLabelV,
			       NULL,1, yAxisLabel,
			       &hProfV, &cpV,
			       delayDraw);
      if (delayDraw) cx->Update();

      if (1) {
	TString figName=TString("fig-") + DYTools::analysisTag +
	  TString("--") + "errorProfiles";
	if (wf.hasExtraTag()) {
	  figName.Append("-");
	  figName.Append(wf.extraFileTag());
	}
	//eliminateSeparationSigns(figName);
	std::cout << "figName=<" << figName << ">\n";
	SaveCanvas(cx,figName);
      }
    }
  }

  if (totalCov) delete totalCov;
  return;
}


// -----------------------------------------------------------

void plotTotCov(TCovData_t &dt, const WorkFlags_t &wf) {
  TMatrixD *totalCov=dt.calcTotalCov();
  for (int iCorr=0; iCorr<4; ++iCorr) {
    if (iCorr==2) continue; // not ready
    TString covStr;
    switch(iCorr) {
    case 0: covStr="Cov_"; break;
    case 1: covStr="Corr_"; break;
    case 2: covStr="partCorr_"; break;
    case 3: covStr="RelCov_"; break;
    default:
      std::cout << "plotAllCovs: unknown iCorr=" << iCorr << "\n";
      return;
    }
    if (wf.showCSCov()) covStr.Prepend("CS");

    TString cName=TString("canvTot") + covStr;
    TCanvas *cx= new TCanvas(cName,cName, 750,700);
    AdjustFor2DplotWithHeight(cx);
    TString explain;
    TMatrixD *plotMatrix=NULL;
    int removePlotMatrix=0;
    TString histoName,histoTitle;

    if (iCorr==0) {
      //set_nice_style(51);
      gStyle->SetPalette(1);
      plotMatrix=totalCov; removePlotMatrix=0;
      histoName="hCov";
      histoTitle="Total covariance";
    }
    else if (iCorr==1) {
      gStyle->SetPalette(1);
      plotMatrix= corrFromCov( *totalCov );
      removePlotMatrix=1;
      histoName=TString("hCorr");
      histoTitle=TString("Total correlations");
    }
    else if (iCorr==3) {
      TH2D *h2=loadMainCSResult();
      TMatrixD *csValAsM=createMatrixD(h2,0);
      if (!csValAsM) return;
      TVectorD csV(DYTools::nUnfoldingBins);
      if (!flattenMatrix(*csValAsM,csV)) return;
      //if (DYTools::study2D==0) {
      //	std::cout << "changing last value\n";
      //	csV(39)*=10;
      //	csV(40)=70;
      //      }
      plotMatrix=relativeCov(csV,*totalCov);
      removePlotMatrix=1;
      if (!plotMatrix) return;
      if (0) { // check
	std::cout << "check the numbers\n";
	printHisto(h2);
	std::cout << "same histo as matrix\n";
	printMatrix("cs",*csValAsM,0);
	std::cout << "same histo as vector\n";
	csV.Print();
	std::cout << "total covariance\n";
	printMatrix("totalCov",*totalCov,1);
	std::cout << "relative covariance\n";
	printMatrix("covDivCS",*plotMatrix,1);
      }
      delete csValAsM;
      delete h2;
      histoName="hPartCov";
      histoTitle="Partial covariance";
    }

    TH2D* h2=createHisto2D(*plotMatrix,NULL,
			   histoName+TString("_base"),
			   histoTitle+TString("_base"),
			   _colrange_none,0,1.0001);

    TH2D* h2save=clipToAnalysisUnfBins(h2,histoName,histoTitle,1); // reset axis!
    /*
    int rangeMin=(DYTools::study2D) ? 25  : 1;
    int rangeMax=(DYTools::study2D) ? 156 : DYTools::nMassBins;
    TH2D *h2save=extractSubArea(h2,rangeMin,rangeMax,rangeMin,rangeMax,histoName,0,1); // reset axis!
    h2save->SetTitle(histoTitle);
    */
    h2save->Draw("COLZ");

    if (explain.Length()) {
      TText *txt=new TText();
      txt->DrawTextNDC(0.25,0.93,explain);
    }
    cx->Update();

    if (1) {
      TString figName=TString("fig-") + DYTools::analysisTag + TString("--total-") + covStr;
      //eliminateSeparationSigns(figName);
      std::cout << "figName=<" << figName << ">\n";
      SaveCanvas(cx,figName);
    }
  }
}

// -----------------------------------------------------------

int workWithData(TCovData_t &dt, const WorkFlags_t &wf) {

  if ((wf.theCase()==0) || (wf.theCase()==2) || (wf.theCase()==3)
      || (wf.theCase()==4)) {
    plotAllCovs(dt,wf);
  }
  else if (wf.theCase()==1) {
    plotTotCov(dt,wf);
  }

  return 1;
}




// -----------------------------------------------------------
// -----------------------------------------------------------
