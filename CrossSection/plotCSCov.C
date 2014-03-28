#include "../Include/DYTools.hh"
#include "../CrossSection/crossSectionFnc.hh"
#include "../Include/colorPalettes.hh"
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

int calc_ESFtot=0;  // wrong correlations
int calc_ESFtotCheck=1;

int doCalcESFCov=(calc_ESFtot+calc_ESFtotCheck) ? 1:0;

int calc_AccFSR=1 * (1-DYTools::study2D);
int calc_AccRnd=1 * (1-DYTools::study2D);

int doCalcAccCov=(calc_AccFSR + calc_AccRnd) ? 1:0;

int calc_FsrFSR=1; // not ready!
int calc_FsrRnd=1;

int doCalcFSRCov=(calc_FsrFSR + calc_FsrRnd) ? 1:0;

// -----------------------------------------------------------

struct TCovData_t {
  std::vector<int> isActive;
  std::vector<TMatrixD*> covYieldV, covUnfV, covEffV, covEsfV;
  std::vector<TMatrixD*> covAccV, covFsrV;
  std::vector<TString> labelYieldV, labelUnfV, labelEffV, labelEsfV;
  std::vector<TString> labelAccV, labelFsrV;
public:

  // ----------------

  TCovData_t() : isActive(0),
		 covYieldV(), covUnfV(), covEffV(), covEsfV(),
		 covAccV(), covFsrV(),
		 labelYieldV(), labelUnfV(), labelEffV(), labelEsfV(),
		 labelAccV(), labelFsrV()
  {
    isActive.reserve(6);
    for (int i=0; i<6; ++i) isActive.push_back(0);
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

int loadYieldCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels);
int loadUnfCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels);
int loadEffCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels);
int loadEsfCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels);
int loadAccCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels);
int loadFsrCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels);

int workWithData(TCovData_t &dt, int the_case);


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
    dt.isActive[0]=1;
  }
  if (res && doCalcUnfCov) { 
    if (!loadUnfCovMatrices(fnameBase,dt.covUnfV,dt.labelUnfV)) return 0;
    dt.isActive[1]=1;
  }
  if (res && doCalcEffCov) { 
    if (!loadEffCovMatrices(fnameBase,dt.covEffV,dt.labelEffV)) return 0;
    dt.isActive[2]=1;
  }
  if (res && doCalcESFCov) {
    if (!loadEsfCovMatrices(fnameBase,dt.covEsfV,dt.labelEsfV)) return 0;
    dt.isActive[3]=1;
  }
  if (res && doCalcAccCov) { 
    if (!loadAccCovMatrices(fnameBase,dt.covAccV,dt.labelAccV)) return 0;
    dt.isActive[4]=1;
  }
  if (res && doCalcFSRCov) { 
    if (!loadFsrCovMatrices(fnameBase,dt.covFsrV,dt.labelFsrV)) return 0;
    dt.isActive[5]=1;
  }

  if (!workWithData(dt,the_case)) return 0;

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


int loadEsfCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels) {

  TString inpFileName=fnameBase;
  inpFileName.ReplaceAll(".root","-esfOnly.root");
  TFile fin(inpFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << inpFileName << ">\n";
    return 0;
  }
  TMatrixD *ptr;
  if (calc_ESFtot) {
    ptr=(TMatrixD*)fin.Get("covCS_ESFtot");
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("ESF tot");
    }
  }
  else if (calc_ESFtotCheck) {
    ptr=(TMatrixD*)fin.Get("covCS_ESFtotCheck");
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
  inpFileName.ReplaceAll(".root","-fsrUnfOnly.root");
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
// -----------------------------------------------------------

void plotAllCovs(TCovData_t &dt) {
  //gStyle->SetPalette(1);

  TMatrixD *totalCov=dt.calcTotalCov();

  for (int iCorr=0; iCorr<3; ++iCorr) {
    if (iCorr!=2) continue;
    TString covStr;
    switch(iCorr) {
    case 0: covStr="Cov_"; break;
    case 1: covStr="Corr_"; break;
    case 2: covStr="partCorr_"; break;
    default:
      std::cout << "plotAllCovs: unknown iCorr=" << iCorr << "\n";
      return;
    }

    for (int idx=0; idx<6; ++idx) {
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
	eliminateSeparationSigns(histoTag);
	if (iCorr==0) {
	  set_nice_style(51);
	  (*covV)[i]->Draw("COLZ");
	  explain="Covariance ";
	}
	else if (iCorr==1) {
	  gStyle->SetPalette(1);
	  TMatrixD *corr= corrFromCov( *(*covV)[i] );
	  TString histoName=TString("hCorr") + histoTag;
	  TString histoTitle=TString("Correlations ") + histoTag;
	  TH2D* h2=createHisto2D(*corr,NULL,histoName,histoTitle,_colrange_center,0,1.0001);
	  h2->Draw("COLZ");
	  delete corr;
	  //explain="Correlation ";
	}
	else if (iCorr==2) {
	  if (DYTools::study2D) set_center_white_style(21);
	  else set_center_white_style5(21);
	  TMatrixD *corr= partialCorrFromCov( *totalCov, *(*covV)[i] );
	  TString histoName=TString("hPartCorr") + histoTag;
	  TString histoTitle=TString("Partial correlations ") + histoTag;
	  TH2D* h2=createHisto2D(*corr,NULL,histoName,histoTitle,_colrange_center,0,1.0001);
	  h2->Draw("COLZ");
	  delete corr;
	  //explain="Correlation ";
	}
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
  if (totalCov) delete totalCov;
  return;
}


// -----------------------------------------------------------

void plotTotCov(TCovData_t &dt) {
  TMatrixD *totalCov=dt.calcTotalCov();
  for (int iCorr=0; iCorr<2; ++iCorr) {
    TString covStr;
    switch(iCorr) {
    case 0: covStr="Cov_"; break;
    case 1: covStr="Corr_"; break;
    case 2: covStr="partCorr_"; break;
    default:
      std::cout << "plotAllCovs: unknown iCorr=" << iCorr << "\n";
      return;
    }
    TString cName=TString("canvTot") + covStr;
    TCanvas *cx= new TCanvas(cName,cName, 750,700);
    AdjustFor2DplotWithHeight(cx);
    TString explain;

    if (iCorr==0) {
      set_nice_style(51);
      totalCov->Draw("COLZ");
      explain="Total covariance ";
    }
    else if (iCorr==1) {
      gStyle->SetPalette(1);
      TMatrixD *corr= corrFromCov( *totalCov );
      TString histoName=TString("hCorr");
      TString histoTitle=TString("Total correlations");
      TH2D* h2=createHisto2D(*corr,NULL,histoName,histoTitle,_colrange_center,0,1.0001);
      h2->Draw("COLZ");
      delete corr;
      //explain="Correlation ";
    }
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

int workWithData(TCovData_t &dt, int the_case) {

  if (the_case==0) {
    plotAllCovs(dt);
  }
  else if (the_case==1) {
    plotTotCov(dt);
  }

  return 1;
}




// -----------------------------------------------------------
// -----------------------------------------------------------
