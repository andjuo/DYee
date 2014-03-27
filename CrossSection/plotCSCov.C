#include "../Include/DYTools.hh"
#include "../CrossSection/crossSectionFnc.hh"
#include <TBenchmark.h>

//=== Global flags =================================================================================================

int calc_YieldStat=1;
int calc_YieldSyst=1;
int calc_YieldUnregEn=1;
int calc_YieldEScale=1;
int calc_YieldApplyEScale=1;

int doCalcYieldCov=(calc_YieldStat + calc_YieldSyst + calc_YieldUnregEn +
		    calc_YieldEScale + calc_YieldApplyEScale) ? 1:0;

int calc_UnfPU=1;
int calc_UnfFSR=1;
int calc_UnfRnd=1;

int doCalcUnfCov=(calc_UnfPU + calc_UnfFSR + calc_UnfRnd) ? 1:0;

int calc_EffPU=1;
int calc_EffFSR=1;
int calc_EffRnd=1;

int doCalcEffCov=(calc_EffPU + calc_EffFSR + calc_EffRnd) ? 1:0;

int calc_AccFSR=1;
int calc_AccRnd=1;

int doCalcAccCov=(calc_AccFSR + calc_AccRnd) ? 1:0;

int calc_FsrFSR=1; // not ready!
int calc_FsrRnd=1;

int doCalcFSRCov=(calc_FsrFSR + calc_FsrRnd) ? 1:0;


struct TCovData_t {
  std::vector<TMatrixD*> covYieldV, covUnfV, covEffV, covAccV, covFsrV;
  std::vector<TString> labelYieldV, labelUnfV, labelEffV, labelAccV, labelFsrV;
};


//=== Functions =================================================================================================

int loadYieldCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels);
int loadUnfCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels);
int loadEffCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels);
int loadAccCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels);
int loadFsrCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels);


//=== MAIN MACRO =================================================================================================


int plotCSCov(TString conf, int the_case ) 
{


  // Settings 
  //==============================================================================================================

  if (conf==TString("default")) {
    conf="../config_files/data_vilnius8TeV_regSSD.conf.py";
  }

  DYTools::TCrossSectionKind_t csKind=(DYTools::study2D) ? DYTools::_cs_preFsrDet : DYTools::_cs_preFsr;

  DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN;
  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;

  TCovData_t dt;

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

  int res=1;
  if (res && doCalcYieldCov) { 
    if (!loadYieldCovMatrices(fnameBase,dt.covYieldV,dt.labelYieldV)) return 0;
  }
  if (res && doCalcUnfCov) { 
    if (!loadUnfCovMatrices(fnameBase,dt.covUnfV,dt.labelUnfV)) return 0;
  }
  if (res && doCalcEffCov) { 
    if (!loadUnfCovMatrices(fnameBase,dt.covEffV,dt.labelEffV)) return 0;
  }
  if (res && doCalcAccCov) { 
    if (!loadUnfCovMatrices(fnameBase,dt.covAccV,dt.labelAccV)) return 0;
  }
  if (res && doCalcFSRCov) { 
    if (!loadFsrCovMatrices(fnameBase,dt.covFsrV,dt.labelFsrV)) return 0;
  }


  return retCodeOk;
}

// ---------------------------------------------------------------------------
  // Implementations
  //==============================================================================================================


int loadYieldCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels) {

  TString inpFileName=fnameBase;
  inpFileName.ReplaceAll(".root","-yieldOnly.root");
  TFile fin(inpFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << inpFileName << ">\n";
    return 0;
  }
  TMatrixD *ptr;
  if (calc_YieldStat) {
    ptr=(TMatrixD*)fin.Get("covCS_YieldStat");
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("signal stat");
    }
  }
  if (calc_YieldSyst) {
    ptr=(TMatrixD*)fin.Get("covCS_YieldSyst");
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("signal syst");
    }
  }
  if (calc_YieldUnregEn) {
    ptr=(TMatrixD*)fin.Get("covCS_YieldUnregEn");
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("signal unreg.en.");
    }
  }
  if (calc_YieldEScale) {
    ptr=(TMatrixD*)fin.Get("covCS_YieldEScale");
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("signal EScale uncert.");
    }
  }
  if (calc_YieldApplyEScale) {
    ptr=(TMatrixD*)fin.Get("covCS_YieldApplyEScale");
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


int loadUnfCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels) {

  TString inpFileName=fnameBase;
  inpFileName.ReplaceAll(".root","-unfOnly.root");
  TFile fin(inpFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << inpFileName << ">\n";
    return 0;
  }
  TMatrixD *ptr;
  if (calc_UnfPU) {
    ptr=(TMatrixD*)fin.Get("covCS_UnfPU");
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("unf pile-up");
    }
  }
  if (calc_UnfFSR) {
    ptr=(TMatrixD*)fin.Get("covCS_UnfFSR");
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("unf FSR");
    }
  }
  if (calc_UnfRnd) {
    ptr=(TMatrixD*)fin.Get("covCS_UnfRnd");
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("unf stat");
    }
  }

  std::cout << "loaded " << covs.size() << " entries from file <" << fin.GetName() << ">\n";
  fin.Close();
  return 1;
}


// -----------------------------------------------------------


int loadEffCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels) {

  TString inpFileName=fnameBase;
  inpFileName.ReplaceAll(".root","-effOnly.root");
  TFile fin(inpFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << inpFileName << ">\n";
    return 0;
  }
  TMatrixD *ptr;
  if (calc_EffPU) {
    ptr=(TMatrixD*)fin.Get("covCS_EffPU");
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("eff pile-up");
    }
  }
  if (calc_EffFSR) {
    ptr=(TMatrixD*)fin.Get("covCS_EffFSR");
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("eff FSR");
    }
  }
  if (calc_EffRnd) {
    ptr=(TMatrixD*)fin.Get("covCS_EffRnd");
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


int loadAccCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels) {

  TString inpFileName=fnameBase;
  inpFileName.ReplaceAll(".root","-accOnly.root");
  TFile fin(inpFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << inpFileName << ">\n";
    return 0;
  }
  TMatrixD *ptr;
  if (calc_AccFSR) {
    ptr=(TMatrixD*)fin.Get("covCS_AccFSR");
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("acc FSR");
    }
  }
  if (calc_AccRnd) {
    ptr=(TMatrixD*)fin.Get("covCS_AccRnd");
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


int loadFsrCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels) {

  TString inpFileName=fnameBase;
  inpFileName.ReplaceAll(".root","-fsrOnly.root");
  TFile fin(inpFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << inpFileName << ">\n";
    return 0;
  }
  TMatrixD *ptr;
  if (calc_FsrFSR) {
    ptr=(TMatrixD*)fin.Get("covCS_FsrFSR");
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("FSR FSR");
    }
  }
  if (calc_FsrRnd) {
    ptr=(TMatrixD*)fin.Get("covCS_FsrRnd");
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
