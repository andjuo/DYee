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
int calc_UnfRnd=0;

int doCalcUnfCov=(calc_UnfPU + calc_UnfFSR + calc_UnfRnd) ? 1:0;

int calc_EffPU=0;
int calc_EffFSR=0;
int calc_EffRnd=0;

int doCalcEffCov=(calc_EffPU + calc_EffFSR + calc_EffRnd) ? 1:0;

int calc_ESFtot=0;  // wrong correlations
int calc_ESFtotCheck=0;

int doCalcESFCov=(calc_ESFtot+calc_ESFtotCheck) ? 1:0;

int calc_AccFSR=0 * (1-DYTools::study2D);
int calc_AccRnd=0 * (1-DYTools::study2D);

int doCalcAccCov=(calc_AccFSR + calc_AccRnd) ? 1:0;

int calc_FsrFSR=0; // not ready!
int calc_FsrRnd=0;

int doCalcFSRCov=(calc_FsrFSR + calc_FsrRnd) ? 1:0;

// -----------------------------------------------------------

struct WorkFlags_t {
  int fCase;
  int fCSCov;
  TString fExtraTag;
public:
  WorkFlags_t(int the_case=0, int set_showCSCov=1, TString set_extra_tag="") :
    fCase(the_case), fCSCov(set_showCSCov),
    fExtraTag(set_extra_tag)
  {}

  WorkFlags_t(const WorkFlags_t &w) :
    fCase(w.fCase), fCSCov(w.fCSCov),
    fExtraTag(w.fExtraTag)
  {}

  int theCase() const { return fCase; }
  void theCase(int the_case) { fCase=the_case; }
  int showCSCov() const { return fCSCov; }
  void showCSCov(int show) { fCSCov=show; }
  int hasExtraTag() const { return (fExtraTag.Length()>0) ? 1:0; }
  TString extraFileTag() const { return fExtraTag; }
  void extraFileTag(TString setTag) { fExtraTag=setTag; }

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
};

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

int loadYieldCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf);
int loadUnfCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf);
int loadEffCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf);
int loadEsfCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf);
int loadAccCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf);
int loadFsrCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf);

int workWithData(TCovData_t &dt, const WorkFlags_t &wf);


//=== MAIN MACRO =================================================================================================


int plotCSCov(int analysisIs2D, TString conf, int the_case, int showCSCov=1,
	      TString outFileExtraTag_UserInput="")
{

  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return retCodeError;
  }

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

  if (!workWithData(dt,work)) return 0;

  return retCodeOk;
}

// ---------------------------------------------------------------------------
  // Implementations
  //==============================================================================================================


int loadYieldCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf) {

  TString inpFileName=fnameBase;
  inpFileName.ReplaceAll(".root","-yieldOnly.root");
  wf.adjustFName(inpFileName);
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
  inpFileName.ReplaceAll(".root","-unfOnly.root");
  wf.adjustFName(inpFileName);
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

  std::cout << "loaded " << covs.size() << " entries from file <" << fin.GetName() << ">\n";
  fin.Close();
  return 1;
}


// -----------------------------------------------------------


int loadEffCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs, std::vector<TString> &labels, const WorkFlags_t &wf) {

  TString inpFileName=fnameBase;
  inpFileName.ReplaceAll(".root","-effOnly.root");
  wf.adjustFName(inpFileName);
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
  inpFileName.ReplaceAll(".root","-esfOnly.root");
  wf.adjustFName(inpFileName);
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
  inpFileName.ReplaceAll(".root","-accOnly.root");
  wf.adjustFName(inpFileName);
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
  inpFileName.ReplaceAll(".root","-fsrUnfOnly.root");
  wf.adjustFName(inpFileName);
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
// -----------------------------------------------------------

void plotAllCovs(TCovData_t &dt, const WorkFlags_t &wf) {
  //gStyle->SetPalette(1);

  TMatrixD *totalCov=dt.calcTotalCov();

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
      TString csFileName="../root_files_reg/xsec/DY_j22_19712pb/xSec_preFsr_1DpreFsrFullSp.root";
      TString fieldName="hpPreFsrFullSp";
      if (DYTools::study2D) {
	csFileName="../root_files_reg/xsec/DY_j22_19712pb/xSec_preFsrDet_2DpreFsrDet.root";
	fieldName="hpPreFsrDet";
      }
      TFile fin(csFileName,"read");
      if (!fin.IsOpen()) {
	std::cout << "failed to open the cross-section file <"
		  << fin.GetName() << ">\n";
	return;
      }
      TH2D *h2=LoadHisto2D(fin,fieldName,"",1);
      fin.Close();
      if (!h2) return;
      TMatrixD *csValAsM=createMatrixD(h2,0);
      if (!csValAsM) return;
      TVectorD csV(DYTools::nUnfoldingBins);
      if (!flattenMatrix(*csValAsM,csV)) return;
      if (DYTools::study2D==0) {
	std::cout << "changing last value\n";
	csV(39)*=10;
	csV(40)=70;
      }
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

  if (wf.theCase()==0) {
    plotAllCovs(dt,wf);
  }
  else if (wf.theCase()==1) {
    plotTotCov(dt,wf);
  }

  return 1;
}




// -----------------------------------------------------------
// -----------------------------------------------------------
