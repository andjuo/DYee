#include "../Include/DYTools.hh"
#include "../CrossSection/crossSectionFnc.hh"
#include <TBenchmark.h>

//=== MAIN MACRO =================================================================================================

int calcCSCov(TString conf, int nExps=100,
	      DYTools::TCrossSectionKind_t csKind=DYTools::_cs_preFsr) {
  gBenchmark->Start("calcCSCov");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================

  if (conf==TString("default")) {
    conf="../config_files/data_vilnius8TeV_regSSD.conf.py";
  }

  DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN;
  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;

  int storeDetails=1;

  int calc_YieldStat=1;
  int calc_YieldSyst=1;
  int calc_YieldEScale=1;


  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(conf)) return retCodeError;

  // Construct eventSelector, update mgr and plot directory
  TString extraTag; // empty
  TString plotExtraTag;

  EventSelector_t evtSelector(inpMgr,runMode,systMode,
			      extraTag,plotExtraTag,
			      EventSelector::_selectDefault);

  // Prepare output directory
  inpMgr.crossSectionDir(systMode,1);

  InputArgs_t inpArgs(&inpMgr,systMode);

  int systFileFlag=1;
  TString outFileName=inpMgr.crossSectionFullFileName(systMode,
					       csKind,0,systFileFlag);
  std::cout << "outFileName=<" << outFileName << ">\n";
  TFile outFile(outFileName,"recreate");
  if (!outFile.IsOpen()) {
    std::cout << "failed to (re)create a file <" << outFile.GetName() << ">\n";
    return retCodeError;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================

  int useDDBkg=(inpMgr.userKeyValueAsInt("DDBKG")==1) ? 1:0;
  TString yieldName=(useDDBkg) ? "signalYieldDDbkg" : "signalYieldMCbkg";
  HistoPair2D_t hpSignalYield(yieldName);
  
  TString resCSName=yieldName;
  resCSName.ReplaceAll("signal","cs");
  HistoPair2D_t hpCS(resCSName);

  CSResults_t csResult;

  TMatrixD *covYieldStat=NULL, *covYieldSyst=NULL;
  TMatrixD *covYieldEScale=NULL;

  TMatrixD *covCS_YieldStat=NULL, *covCS_YieldSyst=NULL;
  TMatrixD *covCS_YieldEScale=NULL;

  int res=1;

  // --------------------- Start work

  // Load basic distribution
  const int loadNormalRunSelection=1;
  TString fnameBgSubtracted=inpMgr.signalYieldFullFileName(systMode,loadNormalRunSelection);
  std::cout << "fnameBgSubtracted=<" << fnameBgSubtracted << ">\n";
  res=hpSignalYield.Load(fnameBgSubtracted,1);

  if (res) res= calculateCS(inpArgs,hpSignalYield,csKind,hpCS,csResult);

  //if (storeDetails && res) {
  //  
  //}

  ////////////////////////////
  // Uncertainties in yields
  ////////////////////////////

  TString yieldField;
  TString yieldFieldExtraSyst=(useDDBkg) ? "h2RelDiffDDbkg" : "h2RelDiffMCbkg";

  if (calc_YieldStat || calc_YieldSyst) {
    int saveSilentMode=inpArgs.silentMode();
    inpArgs.silentMode(2);

    for (int iSyst=0; res && (iSyst<2); ++iSyst) {
      if (!iSyst && !calc_YieldStat) continue;
      if ( iSyst && !calc_YieldSyst) continue;
      
      std::vector<TH2D*> vecRnd;
      if (!createRandomizedVec(hpSignalYield,iSyst,nExps,"hRnd_yield_",vecRnd)) {
	std::cout << "failed to create randomized esemble for iSyst=" << iSyst << "\n";
	return retCodeError;
      }

      //printHisto(vecRnd,0,5,2);

      int unbiasedEstimate=1;
      TH2D* avgDistr=createBaseH2(Form("hAvgDistr_%d",iSyst));
      TH2D* csAvgDistr=createBaseH2(Form("hCSAvgDistr_%d",iSyst));
      TMatrixD* cov= deriveCovMFromRndStudies(vecRnd,unbiasedEstimate,avgDistr);
      TMatrixD* csCov=NULL;
      if (iSyst) covYieldSyst=cov; else covYieldStat=cov;

      if (!cov) res=0;
      else {
	std::vector<TH2D*> csRndV;
	res=calcVecOfCSdistributions(inpArgs,vecRnd,csKind,csRndV);
	csCov= deriveCovMFromRndStudies(csRndV,unbiasedEstimate,csAvgDistr);
	if (iSyst) covCS_YieldSyst=csCov; else covCS_YieldStat=csCov;
	if (!csCov) res=0;

	//printHisto(csRndV,0,5,2);

	ClearVec(csRndV);
      }
      ClearVec(vecRnd);
    }
    
    inpArgs.silentMode(saveSilentMode);
  }

  outFile.Close();
  gBenchmark->Show("calcCSCov");
  return retCodeOk;

}
