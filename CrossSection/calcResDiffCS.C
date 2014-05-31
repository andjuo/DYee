#include <TBenchmark.h>

#include "../Include/DYTools.hh"
#include "../CrossSection/crossSectionFnc.hh"
#include "../Include/UnfoldingMatrix.h"


//=== MAIN MACRO =================================================================================================

int calcResDiffCS(int analysisIs2D,
		  TString conf,
		  DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN,
		  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST,
		  DYTools::TCrossSectionKind_t csKind=DYTools::_cs_None
		  ) {

  gBenchmark->Start("calcCrossSection");

  {
    DYTools::printExecMode(runMode,systMode);
    const int debug_print=1;
    if (!DYTools::checkSystMode(systMode,debug_print,1, DYTools::NO_SYST))
      {
	return retCodeError;
      }
  }

  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return retCodeError;
  }

  if (csKind==DYTools::_cs_None) {
    csKind=(DYTools::study2D) ? DYTools::_cs_preFsrDet : DYTools::_cs_preFsr;
    std::cout << "default csKind " << CrossSectionKindName(csKind) << "\n";
  }

  //--------------------------------------------------------------------------------------------------------------
  // Settings
  //==============================================================================================================

  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(conf)) return retCodeError;
  InputFileMgr_t inpMgrEScale(inpMgr);

  // Construct eventSelector, update mgr and plot directory
  TString extraTag=""; // empty
  TString plotExtraTag;

  EventSelector_t evtSelectorEScale(inpMgrEScale,runMode,DYTools::APPLY_ESCALE,
				    extraTag,plotExtraTag,
				    EventSelector::_selectDefault);
  EventSelector_t evtSelector(inpMgr,runMode,systMode,
			      extraTag,plotExtraTag,
			      EventSelector::_selectDefault);

  // Prepare output directory
  inpMgr.crossSectionDir(systMode,1);

  InputArgs_t inpArgs("default",&inpMgr,systMode,csKind);
  inpArgs.resNameBaseAppend("-test");

  //codeDebugFilePath="/home/andriusj/cms/DYee8TeV-20140403/root_files/constants/DY_j22_19712pb/";
  inpArgs.noSave(codeDebugFilePath.Length());

  InputFileMgr_t inpMgrResDiff(inpMgr);
  inpMgrResDiff.rootFileBaseDir("/media/sdb/andriusj/Results-DYee-unfResidual/root_files_reg/");
  std::cout << dashline;
  std::cout << "rootFileBaseDir=<"
	    << inpMgrResDiff.rootFileBaseDir() << ">\n";
  std::cout << dashline;

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code
  //==============================================================================================================

  TString yieldName=(1 && (inpMgr.userKeyValueAsInt("DDBKG")==1)) ?
    "signalYieldDDbkg" : "signalYieldMCbkg";
  HistoPair2D_t hpSignalYield(yieldName);

  TString resCSName=yieldName;
  resCSName.ReplaceAll("signal","cs");
  HistoPair2D_t hpCS(resCSName);

  TString resCSNameResDiff=resCSName;
  resCSNameResDiff.Append("resDiff");
  HistoPair2D_t hpCSResDiff(resCSNameResDiff);

  CSResults_t csResult;
  CSResults_t csResultResDiff;

  // --------------------- Start work

  // load yields
  int res=1;
  if (res) {
    const int loadNormalRunSelection=1;
    DYTools::TSystematicsStudy_t yieldSystMode=systMode;
    TString fnameBgSubtracted;
    yieldSystMode=DYTools::APPLY_ESCALE;
    //inpArgs.resNameBaseAppend("-peakPosCorr");
    fnameBgSubtracted=inpMgrEScale.signalYieldFullFileName(yieldSystMode,
						     loadNormalRunSelection);
    std::cout << "fnameBgSubtracted=<" << fnameBgSubtracted << ">\n";
    res=hpSignalYield.Load(fnameBgSubtracted,1);
  }
  if (res) res=saveResult(inpArgs,hpSignalYield,"raw");

  // unfolding correction
  HistoPair2D_t hpUnfoldedYield("unfYield");
  if (res) res=unfoldDetResolution(inpArgs, hpSignalYield, hpUnfoldedYield);

  if (0 && res) {
    UnfoldingMatrix_t UnfM(UnfoldingMatrix::_cDET_Response,"detResponse");
    TString outputDir=inpMgr.constDir(systMode,0);
    TString fnameTag=UnfoldingMatrix_t::generateFNameTag(systMode,-1);
    if (!UnfM.autoLoadFromFile(outputDir,fnameTag)) {
      return retCodeError;
    }
    HistoPair2D_t hpUnfChk("unfYieldChk");
    if (!unfold_reco2true(hpUnfChk,UnfM,hpSignalYield)) {
      return retCodeError;
    }
    printHisto(hpUnfoldedYield);
    printHisto(hpUnfChk);
    return retCodeStop;
  }

  HistoPair2D_t hpUnfoldedYieldResDiff("unfYieldResDiff");
  if (res) {
    UnfoldingMatrix_t UnfM(UnfoldingMatrix::_cDET_Response,"detResponse_0_nonRnd");
    DYTools::TSystematicsStudy_t systModeRes=DYTools::ESCALE_RESIDUAL;
    TString outputDir=inpMgrResDiff.constDir(systModeRes,0);
    TString fnameTag=UnfoldingMatrix_t::generateFNameTag(systModeRes,-1);
    if (!UnfM.autoLoadFromFile(outputDir,fnameTag)) {
      return retCodeError;
    }
    if (!unfold_reco2true(hpUnfoldedYieldResDiff,UnfM,hpSignalYield)) {
      return retCodeError;
    }
  }

  //if (res) res=saveResult(inpArgs,hpUnfoldedYield,"unf");

  inpArgs.needsDetUnfolding(0);
  //inpArgsNET.allNormErrorIsSyst(1);
  if (res) res=calculateCS(inpArgs,hpUnfoldedYield,csKind,hpCS,csResult);
  if (res) res=calculateCS(inpArgs,hpUnfoldedYieldResDiff,csKind,hpCSResDiff,csResultResDiff);

  if (res) res=saveResult(inpArgs,hpCS,"");
  if (res) {
    inpArgs.resNameBase("-testResDiff");
    res=saveResult(inpArgs,hpCSResDiff,"");
  }


  if (1 && res) {
    HistoPair2D_t hpRatio("hpRatio");
    hpRatio.divide(hpCS,hpCSResDiff.histo());
    std::vector<TH2D*> histosV;
    std::vector<TString> labelsV;
    if (0) {
      histosV.push_back(hpRatio.histo());
      labelsV.push_back("data");
    }
    else {
      histosV.push_back(hpCS.histo());
      histosV.push_back(hpCSResDiff.histo());
      labelsV.push_back("regular xsec");
      labelsV.push_back("xsec with res.diff");

      PrintHisto2Dvec("\n\nCompare the distributions",histosV,0,-1);
    }
    TCanvas *cx= plotProfiles("cx",histosV,labelsV,NULL,1,
			      "count",NULL,NULL,0);
    cx->Update();
    return retCodeStop;
  }

  //gBenchmark->Show("calcCrossSection");
  ShowBenchmarkTime("calcCrossSection");
  return (res) ? retCodeOk : retCodeError;
}


  //--------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------
