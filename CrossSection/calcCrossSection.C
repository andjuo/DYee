#include <TBenchmark.h>

#include "../Include/DYTools.hh"
#include "../CrossSection/crossSectionFnc.hh"


//=== MAIN MACRO =================================================================================================

int calcCrossSection(int analysisIs2D,
		     TString conf,
		     DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN,
		     DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST,
		     DYTools::TCrossSectionKind_t csKind=DYTools::_cs_None,
		     int special_case=-1) {

  gBenchmark->Start("calcCrossSection");

  {
    DYTools::printExecMode(runMode,systMode);
    const int debug_print=1;
    if (!DYTools::checkSystMode(systMode,debug_print,5, DYTools::NO_SYST,
				DYTools::FSR_5plus, DYTools::FSR_5minus,
			DYTools::PILEUP_5plus, DYTools::PILEUP_5minus)) {
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

  //codeDebugFilePath="/home/andriusj/cms/DYee8TeV-20140403/root_files/constants/DY_j22_19712pb/";
  inpArgs.noSave(codeDebugFilePath.Length());

  // additional manipulation
  if (special_case!=-1) {
    if (special_case<10) {
      //covRhoFileSF_nMB41_asymHLT_FSR_5minus-allSyst_100.root
      //covRhoFileSF_nMB41_asymHLT_FSR_5plus-allSyst_100.root
      //covRhoFileSF_nMB41_asymHLT_Pileup5minus-allSyst_100.root
      //covRhoFileSF_nMB41_asymHLT_Pileup5plus-allSyst_100.root
      //covRhoFileSF_nMB41_asymHLT_regEn-allSyst_100.root
      //covRhoFileSF_nMB41_asymHLT_Unregressed_energy-allSyst_100.root
      TString tag;
      TString id;
      switch(special_case) {
      case 0: tag="Unregressed_energy-allSyst_100"; id="-ESFUnregEn"; break;
      case 1: tag="regEn-allSyst_100"; id="-ESFregEn"; break;
      case 2: tag="FSR_5minus-allSyst_100"; id="-ESFFSR5m"; break;
      case 3: tag="FSR_5plus-allSyst_100"; id="-ESFFSR5p"; break;
      case 4: tag="Pileup5minus-allSyst_100"; id="-ESFPU5m"; break;
      case 5: tag="Pileup5plus-allSyst_100"; id="-ESFPU5p"; break;
      default: ;
      }
      if (!tag.Length()) {
	std::cout << "special_case for event efficiency scale factors failed\n";
	return retCodeError;
      }
      inpArgs.resNameBaseAppend(id);
      std::cout << "event efficiency scale factor tag=<" << tag << ">\n";
      inpMgr.addUserKey("SpecialESFTag",tag);
    }
  }

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================

  TString yieldName=(1 && (inpMgr.userKeyValueAsInt("DDBKG")==1)) ?
    "signalYieldDDbkg" : "signalYieldMCbkg";
  HistoPair2D_t hpSignalYield(yieldName);

  TString resCSName=yieldName;
  resCSName.ReplaceAll("signal","cs");
  HistoPair2D_t hpCS(resCSName);

  CSResults_t csResult;
  
  // --------------------- Start work

  // load yields 
  int res=1;
  if (res) {
    const int loadNormalRunSelection=1;
    DYTools::TSystematicsStudy_t yieldSystMode=systMode;
    TString fnameBgSubtracted;
    if (0) {
      yieldSystMode=DYTools::ESCALE_DIFF_0000;
      inpArgs.resNameBaseAppend("-regEn");
      fnameBgSubtracted=
          inpMgr.signalYieldFullFileName(yieldSystMode,loadNormalRunSelection);
    }
    else {
      yieldSystMode=DYTools::APPLY_ESCALE;
      //inpArgs.resNameBaseAppend("-peakPosCorr");
      fnameBgSubtracted=inpMgrEScale.signalYieldFullFileName(yieldSystMode,
						     loadNormalRunSelection);
    }
    const int load_debug_file=(codeDebugFilePath.Length()) ? 1:0;
    if ( ! load_debug_file ) {
      std::cout << "fnameBgSubtracted=<" << fnameBgSubtracted << ">\n";
      res=hpSignalYield.Load(fnameBgSubtracted,1);
    }
    else {
      TString tmpFName=codeDebugFilePath +
	    TString(Form("yields_bg-subtracted%dD.root",DYTools::study2D+1));
      tmpFName.ReplaceAll("constants","yields");
      std::cout << "debug fnameBgSubtracted=<" << tmpFName << ">\n";
      TString field="YieldsSignal";
      TString fieldErr="YieldsSignalErr";
      TString fieldSystErr="YieldsSignalSystErr";
      res=hpSignalYield.loadThreeMatrices(tmpFName,field,fieldErr,fieldSystErr,
					  1,1);
    }
  }
  if (res) res=saveResult(inpArgs,hpSignalYield,"raw");

  // unfolding correction
  HistoPair2D_t hpUnfoldedYield("unfYield");
  if (res) res=unfoldDetResolution(inpArgs, hpSignalYield, hpUnfoldedYield);
  //if (res) res=addSystError_DetResUnf_unfold(inpArgs, hpUnfoldedYield);
  //if (res) res=addSystError_DetResUnf_escale(inpArgs, hpUnfoldedYield);
  if (res) res=saveResult(inpArgs,hpUnfoldedYield,"unf");


  if (0) {
    // ratio of the unfolded to pre-unfolded yield
    HistoPair2D_t hpRatio("hpRatio");
    hpRatio.divide(hpSignalYield,hpUnfoldedYield.histo());
    std::vector<TH2D*> histosV;
    std::vector<TString> labelsV;
    if (0) {
      histosV.push_back(hpRatio.histo());
      labelsV.push_back("data");
    }
    else {
      histosV.push_back(hpUnfoldedYield.histo());
      histosV.push_back(hpSignalYield.histo());
      labelsV.push_back("unfolded yield");
      labelsV.push_back("signal yield");
    }
    TCanvas *cx= plotProfiles("cx",histosV,labelsV,NULL,1,
			      "count",NULL,NULL,0);
    cx->Update();
    return retCodeStop;
  }

  if (0) {
    HistoPair2D_t hpUnf_FSR5minus("unf_FSR5minus");
    HistoPair2D_t hpUnf_FSR5plus("unf_FSR5plus");
    InputArgs_t iaFSR5minus("iaFSR5minus",inpArgs,"-fsr5minus");
    InputArgs_t iaFSR5plus("iaFSR5plus",inpArgs,"-fsr5plus");
    iaFSR5minus.systMode(DYTools::FSR_5minus);
    iaFSR5plus.systMode(DYTools::FSR_5plus);
    unfoldDetResolution(iaFSR5minus,hpSignalYield,hpUnf_FSR5minus);
    unfoldDetResolution(iaFSR5plus,hpSignalYield,hpUnf_FSR5plus);
    std::vector<TH2D*> histosV;
    std::vector<TString> labelsV;

    histosV.push_back(hpSignalYield.histo());
    histosV.push_back(hpUnfoldedYield.histo());
    histosV.push_back(hpUnf_FSR5minus.histo());
    histosV.push_back(hpUnf_FSR5plus.histo());
    labelsV.push_back("signal yield");
    labelsV.push_back("unfolded yield");
    labelsV.push_back("unfolded yield (FSR weight-5%)");
    labelsV.push_back("unfolded yield (FSR weight+5%)");

    TCanvas *cx= plotProfiles("cx",histosV,labelsV,NULL,1,
			      "count",NULL,NULL,0);
    cx->Update();
    return retCodeStop;
  }

  // --------------------------------------
  // Normal calculation
  // --------------------------------------
  inpArgs.needsDetUnfolding(0);
  //inpArgsNET.allNormErrorIsSyst(1);
  if (res) res=calculateCS(inpArgs,hpUnfoldedYield,csKind,hpCS,csResult);

  if (res) res=saveResult(inpArgs,hpCS,"");

  // -----------------------------------------------
  // Get the error propagation estimate
  // -----------------------------------------------

  if (0 && res && (systMode==DYTools::NO_SYST)) {
    std::cout << dashline << "Working on error propagation estimate\n"
	      << dashline;
    InputArgs_t iaNoExtraErr("iaNoExtraErr",inpArgs,"noExtraErr");
    iaNoExtraErr.silentMode(0);
    iaNoExtraErr.needsDetUnfolding(0);
    iaNoExtraErr.includeCorrError(0);

    InputArgs_t iaWithExtraErr("iaWithExtraErr",iaNoExtraErr,"wExtraErr");
    iaWithExtraErr.includeCorrError(1);

    // Result containers
    HistoPair2D_t hpCS_yieldErr("hpCS_yieldErr");
    HistoPair2D_t hpCS_effErr("hpCS_effErr");
    HistoPair2D_t hpCS_accErr("hpCS_accErr");
    CSResults_t csResult_yieldErr;
    CSResults_t csResult_effErr;
    CSResults_t csResult_accErr;

    // 1. calculate propagated yield error
    if (res) res= calculateCS(iaNoExtraErr,hpUnfoldedYield,csKind,
			      hpCS_yieldErr,csResult_yieldErr);

    // 2. calculate propagated efficiency error
    // First construct unf distribution without error.
    // Then obtain the efficiency correction with error, and finish
    // without additional errors
    HistoPair2D_t hpUnfYieldNoErr("hpUnfYieldNoErr",hpUnfoldedYield);
    removeError(hpUnfYieldNoErr.editHisto());
    hpUnfYieldNoErr.editHistoSystErr()->Reset();

    iaWithExtraErr.needsEffCorr(1);
    iaWithExtraErr.needsEffScaleCorr(0);
    HistoPair2D_t hpEffCorrYield("hpEffCorrYield");
    if (res) res= calculateCSdistribution(iaWithExtraErr,hpUnfYieldNoErr,
				    DYTools::_cs_postFsrDet,hpEffCorrYield);
    iaNoExtraErr.needsEffCorr(0);
    if (res) res= calculateCS(iaNoExtraErr,hpEffCorrYield,csKind,
			      hpCS_effErr,csResult_effErr);

    // 3. calculate propagated acceptance error
    // First produce post-FSR cross section, ignoring the error effects
    // (hpEffCorrYield does not contain the efficiency scale correction)
    // FSR Unfolding does not add extra error
    HistoPair2D_t hpPostFsrCS("hpPostFsrCS");
    iaNoExtraErr.needsDetUnfolding(0);
    iaNoExtraErr.needsEffCorr(1);
    iaNoExtraErr.needsEffScaleCorr(1);
    iaNoExtraErr.needsAccCorr(0);
    iaNoExtraErr.needsFsrCorr(0);
    if (res) res= calculateCSdistribution(iaNoExtraErr,hpUnfYieldNoErr,
					  DYTools::_cs_postFsrDet,hpPostFsrCS);

    iaWithExtraErr.needsDetUnfolding(0);
    iaWithExtraErr.needsEffCorr(0);
    iaWithExtraErr.needsEffScaleCorr(0);
    iaWithExtraErr.needsAccCorr(1);
    iaWithExtraErr.needsFsrCorr(1);
    if (res) res= calculateCS(iaWithExtraErr,hpPostFsrCS,csKind,
			      hpCS_accErr,csResult_accErr);

    // check the consistency
    if (0) {
      HERE(dashline.c_str());
      HERE("the central values have to match perfectly");
      int truncX=-1;
      int truncY=3;
      printHisto(hpCS,truncX,truncY);
      printHisto(hpCS_yieldErr,truncX,truncY);
      printHisto(hpCS_effErr,truncX,truncY);
      printHisto(hpCS_accErr,truncX,truncY);
    }

    // Get the corrections
    TString effCorrFName= inpMgr.correctionFullFileName("efficiency",systMode,0);
    TH2D *hEff=LoadHisto2D("hEfficiency",effCorrFName,"",1);
    TString accCorrFName= inpMgr.correctionFullFileName("acceptance",systMode,0);
    TH2D *hAcc=(DYTools::study2D) ?
      NULL : LoadHisto2D("hAcceptance",accCorrFName,"",1);


    // Save the results
    TString fname=inpMgr.crossSectionFullFileName(systMode,csKind,1,0);
    TString extraTag_loc="_errProp.root";
    fname.ReplaceAll(".root",extraTag_loc);
    std::cout << "output file for error propagation is <" << fname << ">\n";
    TFile fout(fname,"recreate");
    if (!fout.IsOpen()) return retCodeError;
    if (res) res=hpCS.Write(fout,"","mainCS");
    if (res) res=hpCS_yieldErr.Write(fout,"","mainCS_yieldErr");
    if (res) res=hpCS_effErr.Write(fout,"","mainCS_effErr");
    if (res) res=hpCS_accErr.Write(fout,"","mainCS_accErr");
    if (res) res=saveHisto(fout,hEff,"","hEfficiency");
    if (res && hAcc) res=saveHisto(fout,hAcc,"","hAcceptance");
    if (res) writeBinningArrays(fout,"calcCrossSection");
    fout.Close();
    if (!res) std::cout << "failed to save file <" << fout.GetName() << ">\n";
    else std::cout << "file <" << fout.GetName() << "> saved\n";
  }

  //gBenchmark->Show("calcCrossSection");
  ShowBenchmarkTime("calcCrossSection");
  return (res) ? retCodeOk : retCodeError;
}


  //--------------------------------------------------------------------------------------------------------------
  //--------------------------------------------------------------------------------------------------------------
