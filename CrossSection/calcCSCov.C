#include "../Include/DYTools.hh"
#include "../CrossSection/crossSectionFnc.hh"
#include <TBenchmark.h>

//=== MAIN MACRO =================================================================================================

int calcCSCov(TString conf, int nExps=100,
	      DYTools::TCrossSectionKind_t csKind=DYTools::_cs_preFsr,
	      TString outFileExtraTag_UserInput="") {
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

  int calc_YieldStat=0;
  int calc_YieldSyst=0;
  int calc_YieldUnregEn=0;
  int calc_YieldEScale=0;
  int calc_YieldApplyEScale=0;

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

  int calc_ESFtot=0; // wrong correlations
  int calc_ESFtotCheck=0;

  int doCalcESFCov=(calc_ESFtot+calc_ESFtotCheck) ? 1:0;

  int calc_AccFSR=0;
  int calc_AccRnd=0;

  int doCalcAccCov=(calc_AccFSR + calc_AccRnd) ? 1:0;

  int calc_FsrFSR=0; // not ready
  int calc_FsrRnd=1;

  int doCalcFSRCov=(calc_FsrFSR + calc_FsrRnd) ? 1:0;

  TMatrixD* covCS_fromYield=NULL;
  TMatrixD* covCS_fromUnf=NULL;
  TMatrixD* covCS_fromEff=NULL;
  TMatrixD* covCS_fromESF=NULL;
  TMatrixD* covCS_fromAcc=NULL;
  TMatrixD* covCS_fromFSR=NULL;

  int needsDetResUnfM= (doCalcUnfCov + doCalcEffCov) ? 1:0;
  int needsFSRUnfM   = ((doCalcAccCov + doCalcFSRCov) && DYTools::isFullSpaceCS(csKind)) ? 1:0;
  int needsFSRUnfM_det= (doCalcFSRCov && !DYTools::isFullSpaceCS(csKind)) ? 1:0;

  // The "default" unfolding matrix
  UnfoldingMatrix_t *detResponse=NULL; //(UnfoldingMatrix::_cDET_Response,"detResponse");
  UnfoldingMatrix_t *fsrResponse=NULL; //(UnfoldingMatrix::_cFSR, "fsrGood");
  UnfoldingMatrix_t *fsrResponseDet=NULL; //(UnfoldingMatrix::_cFSR_DET, "fsrDETgood");
  UnfoldingMatrix_t *fsrResponseExact=NULL; //(UnfoldingMatrix::_cFSR, "fsrExact");
  UnfoldingMatrix_t *fsrResponseDetExact=NULL; //(UnfoldingMatrix::_cFSR_DET, "fsrDETexact");

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
  if (outFileExtraTag_UserInput.Length()) {
    std::cout << "applying extra tag from input= " << outFileExtraTag_UserInput << "\n";
    outFileName.ReplaceAll(".root",outFileExtraTag_UserInput);
    if (outFileName.Index(".root")==-1) outFileName.Append(".root");
  }
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

  int res=1;

  // --------------------- Start work

  // Load basic distribution
  const int loadNormalRunSelection=1;
  TString fnameBgSubtracted=inpMgr.signalYieldFullFileName(systMode,loadNormalRunSelection);
  std::cout << "fnameBgSubtracted=<" << fnameBgSubtracted << ">\n";
  res=hpSignalYield.Load(fnameBgSubtracted,1);
  printHisto(hpSignalYield,6);

  if (res) res= calculateCS(inpArgs,hpSignalYield,csKind,hpCS,csResult);


  // load the needed unfolding matrices
  if (res && (needsDetResUnfM || needsFSRUnfM || needsFSRUnfM_det)) {
    TString constDirDef=inpMgr.constDir(DYTools::NO_SYST,0);
    TString fnameTagDef=UnfoldingMatrix_t::generateFNameTag(DYTools::NO_SYST);
    if (res && needsDetResUnfM) {
      detResponse= new UnfoldingMatrix_t(UnfoldingMatrix::_cDET_Response,"detResponse");
      if (!detResponse) res=0;
      if (res) res=detResponse->autoLoadFromFile(constDirDef,fnameTagDef);
    }
    if (res && needsFSRUnfM) {
      fsrResponse= new UnfoldingMatrix_t(UnfoldingMatrix::_cFSR, "fsrGood");
      fsrResponseExact= new UnfoldingMatrix_t(UnfoldingMatrix::_cFSR, "fsrExact");
      if (!fsrResponse || !fsrResponseExact) res=0;
      if (res) res=fsrResponse->autoLoadFromFile(constDirDef,fnameTagDef);
      if (res) res=fsrResponseExact->autoLoadFromFile(constDirDef,fnameTagDef);
    }
    if (res && needsFSRUnfM_det) {
      fsrResponseDet= new UnfoldingMatrix_t(UnfoldingMatrix::_cFSR_DET, "fsrDETgood");
      fsrResponseDetExact= new UnfoldingMatrix_t(UnfoldingMatrix::_cFSR_DET, "fsrDETexact");
      if (!fsrResponseDet || !fsrResponseDetExact) res=0;
      if (res) res=fsrResponseDet->autoLoadFromFile(constDirDef,fnameTagDef);
      if (res) res=fsrResponseDetExact->autoLoadFromFile(constDirDef,fnameTagDef);
    }
    if (!res) {
      std::cout << "failed to load the unfolding matrix\n";
      return retCodeError;
    }
  }

  ////////////////////////////
  // Uncertainties in yields
  ////////////////////////////

  TString yieldFieldExtraSyst=(useDDBkg) ? "h2RelDiffDDbkg" : "h2RelDiffMCbkg";

  if (res && doCalcYieldCov) {
    int saveSilentMode=inpArgs.silentMode();
    inpArgs.silentMode(2);
    TMatrixD *covYieldTot=NULL;

    for (int iSyst=0; res && (iSyst<5); ++iSyst) {
      int run=0;
      HistoPair2D_t *hpIni=NULL;
      DYTools::TSystematicsStudy_t runSystMode=DYTools::NO_SYST;
      TString csCovName;
      switch(iSyst) {
      case 0:
	run=calc_YieldStat; hpIni=&hpSignalYield;
	csCovName="covCS_YieldStat";
	break;
      case 1:
	run=calc_YieldSyst; hpIni=&hpSignalYield;
	csCovName="covCS_YieldSyst";
	break;
      case 2:
	run=calc_YieldUnregEn; runSystMode=DYTools::UNREGRESSED_ENERGY;
	csCovName="covCS_YieldUnregEn";
	break;
      case 3:
	run=calc_YieldEScale; runSystMode=DYTools::ESCALE_STUDY;
	csCovName="covCS_YieldEScale";
	break;
      case 4:
	run=calc_YieldApplyEScale; runSystMode=DYTools::APPLY_ESCALE;
	csCovName="covCS_YieldApplyEScale";
	break;
      default:
	std::cout << "not ready for yields iSyst=" << iSyst << "\n";
	return retCodeError;
      }

      if (!run) continue;
      std::cout << " - will produce " << csCovName << "\n";

      if (!hpIni) {
	int checkBinning=0;
	int ignoreDebugFlag=1;
	TString loadFileName= inpMgr.signalYieldFullFileName(runSystMode,ignoreDebugFlag,0,1);
	TH2D *h2=LoadHisto2D(yieldFieldExtraSyst,loadFileName,"",checkBinning);
	if (!h2) {
	  std::cout << "failed to load systematics\n";
	  return retCodeError;
	}
	h2->SetName("h2tmp");
	hpIni=new HistoPair2D_t(yieldFieldExtraSyst);
	hpIni->add(h2,1.);
	delete h2;
	std::cout << "Loaded "; printHisto(*hpIni,6);
      }
      
      std::vector<TH2D*> vecRnd;
      if (!createRandomizedVec(hpSignalYield,iSyst,nExps,"hRnd_yield_",vecRnd)) {
	std::cout << "failed to create randomized esemble for iSyst=" << iSyst << "\n";
	return retCodeError;
      }

      //printHisto(vecRnd,0,5,2);

      int unbiasedEstimate=1;
      TH2D* avgDistr=createBaseH2(Form("hYieldAvgDistr_%d",iSyst));
      TH2D* csAvgDistr=createBaseH2(Form("hCSYieldAvgDistr_%d",iSyst));
      TMatrixD* cov= deriveCovMFromRndStudies(vecRnd,unbiasedEstimate,avgDistr);
      TMatrixD* csCov=NULL;

      if (cov) {
	if (!covYieldTot) covYieldTot=new TMatrixD(*cov);
	else (*covYieldTot) += (*cov);
      }

      if (!cov) res=0;
      else {
	std::vector<TH2D*> csRndV;
	res=calcVecOfCSdistributions(inpArgs,vecRnd,csKind,csRndV);
	csCov= deriveCovMFromRndStudies(csRndV,unbiasedEstimate,csAvgDistr);
	if (!csCov) res=0;

	if (csCov) {
	  if (!covCS_fromYield) covCS_fromYield=new TMatrixD(*csCov);
	  else (*covCS_fromYield) += (*csCov);
	}

	//printHisto(csRndV,0,5,2);

	ClearVec(csRndV);
      }
      ClearVec(vecRnd);

      if (res) {
	outFile.cd();
	TString covName=csCovName;
	covName.ReplaceAll("covCS","cov");
	if (cov  ) cov->Write(covName);
	if (csCov) csCov->Write(csCovName);
      }

      // release memory
      if (csCov) delete csCov;
      if (cov  ) delete cov;
      if (csAvgDistr) delete csAvgDistr;
      if (avgDistr  ) delete avgDistr;
      if (runSystMode!=DYTools::NO_SYST) {
	hpIni->clear();
	delete hpIni;
      }
    }

    if (res) {
      outFile.cd();
      if (covYieldTot) covYieldTot->Write("covYieldTot");
    }

    // clean-up
    if (covYieldTot) delete covYieldTot;
    
    inpArgs.silentMode(saveSilentMode);
  }

  HERE("after calculating yield cov res=%d",res);

  //////////////////////////////
  // Uncertainties in unfolding
  //////////////////////////////

  if (res && doCalcUnfCov) {
    int saveSilentMode=inpArgs.silentMode();
    inpArgs.silentMode(2);
    TMatrixD *covUnfTot=NULL;

    for (int iSyst=0; res && (iSyst<3); ++iSyst) {
      int run=0;
      DYTools::TSystematicsStudy_t runSystMode=DYTools::NO_SYST;
      TString csCovName;
      // unf matrices for bounds
      TString uNamePlus, uNameMinus;
      DYTools::TSystematicsStudy_t smPlus=DYTools::NO_SYST;
      DYTools::TSystematicsStudy_t smMinus=smPlus;
      UnfoldingMatrix_t *Uplus=NULL, *Uminus=NULL;
      // unf matrix for the reference
      UnfoldingMatrix_t *Uref=NULL;
      switch(iSyst) {
      case 0:
	run=calc_UnfRnd; Uref=detResponse;
	csCovName="covCS_UnfRnd";
	break;
      case 1:
	run=calc_UnfPU;
	runSystMode=DYTools::PU_STUDY;
	uNamePlus ="detResponse_PU5plus" ; smPlus =DYTools::PILEUP_5plus;
	uNameMinus="detResponse_PU5minus"; smMinus=DYTools::PILEUP_5minus;
	csCovName="covCS_UnfPU";
	break;
      case 2:
	run=calc_UnfFSR;
	runSystMode=DYTools::FSR_STUDY;;
	uNamePlus ="detResponse_105"; smPlus =DYTools::FSR_5plus;
	uNameMinus="detResponse_095"; smMinus=DYTools::FSR_5minus;
	csCovName="covCS_UnfFSR";
	break;
      default:
	std::cout << "not ready for unf iSyst=" << iSyst << "\n";
	return retCodeError;
      }

      if (!run) continue;
      std::cout << " - will produce " << csCovName << "\n";

      if (smPlus!=DYTools::NO_SYST) {
	Uplus =new UnfoldingMatrix_t(UnfoldingMatrix::_cDET_Response,uNamePlus);
	Uminus=new UnfoldingMatrix_t(UnfoldingMatrix::_cDET_Response,uNameMinus);
	if (!Uplus || !Uminus) res=0;
	if (res) res= Uplus->autoLoadFromFile(inpMgr.constDir(runSystMode,0),
				 UnfoldingMatrix_t::generateFNameTag(smPlus));
	if (res) res= Uminus->autoLoadFromFile(inpMgr.constDir(runSystMode,0),
				 UnfoldingMatrix_t::generateFNameTag(smMinus));
	if (res) {
	  Uref=new UnfoldingMatrix_t(*detResponse,"detResponseRef");
	  if (!Uref) res=0;
	}
	if (res) res=Uref->setErrorOnMigrationMatrix(*Uplus,*Uminus);
	if (res) {
	  if (Uplus) { delete Uplus; Uplus=NULL; }
	  if (Uminus) { delete Uminus; Uminus=NULL; }
	}
      }

      if (!res) {
	std::cout << "error in preparation the unfolding matrix\n";
	return retCodeError;
      }

      std::vector<TH2D*> vecUnfRnd;
      if (res) {
	UnfoldingMatrix_t Urnd(UnfoldingMatrix::_cDET_Response,"detResponseRnd");
	HistoPair2D_t hpUnf("hpUnf");

	for (int iexp=0; res && (iexp<nExps); ++iexp) {
	  res=Urnd.randomizeMigrationMatrix(*Uref);
	  if (res) res=unfold_reco2true(hpUnf,Urnd,hpSignalYield);
	  if (res) {
	    TString rndVec=Form("h2unfRnd_%d",iexp);
	    TH2D* h2=Clone(hpUnf.histo(),rndVec);
	    if (!h2) res=0;
	    vecUnfRnd.push_back(h2);
	  }
	}
	hpUnf.clear(); // delete histos
      }

      //printHisto(vecUnfRnd,0,5,2);

      int unbiasedEstimate=1;
      TH2D* avgDistr=createBaseH2(Form("hUnfAvgDistr_%d",iSyst));
      TH2D* csAvgDistr=createBaseH2(Form("hCSUnfAvgDistr_%d",iSyst));
      TMatrixD* cov= deriveCovMFromRndStudies(vecUnfRnd,unbiasedEstimate,avgDistr);
      TMatrixD* csCov=NULL;

      if (cov) {
	if (!covUnfTot) covUnfTot=new TMatrixD(*cov);
	else (*covUnfTot) += (*cov);
      }

      if (!cov) res=0;
      else {
	std::vector<TH2D*> csRndV;
	int needsUnf=0;
	InputArgs_t iaUnf(inpArgs,"unfSyst",needsUnf);
	res=calcVecOfCSdistributions(iaUnf,vecUnfRnd,csKind,csRndV);
	csCov= deriveCovMFromRndStudies(csRndV,unbiasedEstimate,csAvgDistr);
	if (!csCov) res=0;

	if (csCov) {
	  if (!covCS_fromUnf) covCS_fromUnf=new TMatrixD(*csCov);
	  else (*covCS_fromUnf) += (*csCov);
	}

	//printHisto(csRndV,0,5,2);

	ClearVec(csRndV);
      }
      ClearVec(vecUnfRnd);

      if (res) {
	outFile.cd();
	TString covName=csCovName;
	covName.ReplaceAll("covCS","cov");
	if (cov  ) cov->Write(covName);
	if (csCov) csCov->Write(csCovName);
      }

      // release memory
      if (csCov) delete csCov;
      if (cov  ) delete cov;
      if (csAvgDistr) delete csAvgDistr;
      if (avgDistr  ) delete avgDistr;
      if ((smPlus!=DYTools::NO_SYST) && Uref) {
	delete Uref;
      }
    }

    // clean-up
    if (covUnfTot) delete covUnfTot;

    inpArgs.silentMode(saveSilentMode);
  }

  HERE("after calculating unfolding cov res=%d",res);

  ////////////////////////////////
  // Uncertainties in MC efficiency
  ////////////////////////////////

  if (res && doCalcEffCov) {
    int saveSilentMode=inpArgs.silentMode();
    inpArgs.silentMode(2);
    TMatrixD *covEffTot=NULL;

    // "default" unfolded yields
    HistoPair2D_t hpUnf("hpUnf");
    res= unfold_reco2true(hpUnf,*detResponse,hpSignalYield);
    printHisto(hpUnf,5);

    for (int iSyst=0; res && (iSyst<3); ++iSyst) {
      int run=0;
      DYTools::TSystematicsStudy_t runSystMode=DYTools::NO_SYST;
      TString csCovName;
      TString systTag;
      switch(iSyst) {
      case 0:
	run=calc_EffRnd;
	csCovName="covCS_EffRnd";
	break;
      case 1:
	run=calc_EffPU;
	runSystMode=DYTools::PU_STUDY; systTag="_PileUpSyst";
	csCovName="covCS_EffPU";
	break;
      case 2:
	run=calc_EffFSR;
	runSystMode=DYTools::FSR_STUDY; systTag="_FSRsyst";
	csCovName="covCS_EffFSR";
	break;
      default:
	std::cout << "not ready for eff iSyst=" << iSyst << "\n";
	return retCodeError;
      }

      if (!run) continue;
      std::cout << " - will produce " << csCovName << "\n";

      std::vector<TH2D*> vecRnd;
      if (res) {
	TString effField=TString("hefficiency") + systTag;
	TString hpCorrName=TString("hpEff") + systTag;
	HistoPair2D_t hpEffCorr(effField);
	int checkBinning=0;
	int applyExtraTag=1;
	TString loadFileName= inpMgr.correctionSystFullFileName("efficiency",DYTools::NO_SYST,applyExtraTag);

	if (res) {
	  TH2D *h2corr=LoadHisto2D(effField,loadFileName,"",checkBinning);
	  if (!h2corr) {
	    std::cout << "failed to load systematics\n";
	    return retCodeError;
	  }
	  h2corr->SetName("h2tmp");
	  hpEffCorr.add(h2corr,1.);
	  delete h2corr;
	}
	std::cout << "Loaded "; printHisto(hpEffCorr,6);


	int useSystErr=0; // randomize within statistical error
	if (!createRandomizedVec(hpEffCorr,useSystErr,nExps,TString("hRnd_eff")+systTag,vecRnd)) {
	  std::cout << "failed to create randomized esemble for iSyst=" << iSyst << "\n";
	  return retCodeError;
	}
	hpEffCorr.clear();
      }

      if (nExps==2) printHisto(vecRnd,0,5,2);

      int unbiasedEstimate=1;
      TH2D* avgDistr=createBaseH2(Form("hEffAvgDistr_%d",iSyst));
      TH2D* csAvgDistr=createBaseH2(Form("hCSEffAvgDistr_%d",iSyst));
      TMatrixD* cov= deriveCovMFromRndStudies(vecRnd,unbiasedEstimate,avgDistr);
      TMatrixD* csCov=NULL;

      if (cov) {
	if (!covEffTot) covEffTot=new TMatrixD(*cov);
	else (*covEffTot) += (*cov);
      }
      if (!cov) res=0;

      // The randomized efficiency values need to be applied
      // on the unfolded vector
      if (res) {
	HistoPair2D_t hpUnfEffYield("hpUnfEffYield");
	for (unsigned int i=0; res && (i<vecRnd.size()); ++i) {
	  TH2D *rndEff=vecRnd[i];
	  removeError(rndEff);
	  if (res) {
	    res=hpUnfEffYield.divide(hpUnf,rndEff);
	    //printHisto(hpUnfEffYield,6);
	  }
	  if (res) {
	    // swap pointers to histograms, since
	    // vecRnd has to contain eff-corrected yields
	    hpUnfEffYield.swapHistoPtr(&vecRnd[i]);
	    // By default, the 1st entry will have name hpUnfEffYield,
	    // which is not correct.
	    // To avoid confusion, rename the histo
	    TString newName=Form("h2UnfEffYield_%d",i);
	    newName+=systTag;
	    vecRnd[i]->SetName(newName);
	    vecRnd[i]->SetTitle(newName);
	  }
	}
	hpUnfEffYield.clear();
      }

      if (nExps==2) printHisto(vecRnd,0,5,2);

      if (res) {
	std::vector<TH2D*> csRndV;
	InputArgs_t inpArgsEff(inpArgs,"effSyst");
	inpArgsEff.needsDetUnfolding(0);
	inpArgsEff.needsEffCorr(0);
	res=calcVecOfCSdistributions(inpArgsEff,vecRnd,csKind,csRndV);
	csCov= deriveCovMFromRndStudies(csRndV,unbiasedEstimate,csAvgDistr);
	if (!csCov) res=0;

	if (csCov) {
	  if (!covCS_fromEff) covCS_fromEff=new TMatrixD(*csCov);
	  else (*covCS_fromEff) += (*csCov);
	}

	//printHisto(csRndV,0,5,2);

	ClearVec(csRndV);
      }
      ClearVec(vecRnd);

      if (res) {
	outFile.cd();
	TString covName=csCovName;
	covName.ReplaceAll("covCS","cov");
	if (cov  ) cov->Write(covName);
	if (csCov) csCov->Write(csCovName);
      }

      // release memory
      if (csCov) delete csCov;
      if (cov  ) delete cov;
      if (csAvgDistr) delete csAvgDistr;
      if (avgDistr  ) delete avgDistr;
    }

    if (res) {
      outFile.cd();
      if (covEffTot) covEffTot->Write("covEffTot");
    }

    // clean-up
    hpUnf.clear();
    if (covEffTot) delete covEffTot;

    inpArgs.silentMode(saveSilentMode);
  }

  HERE("after calculating eff cov res=%d",res);

  //////////////////////////////////////////////
  // Uncertainties in efficiency scale factors
  //////////////////////////////////////////////

  if (res && doCalcESFCov) {
    int saveSilentMode=inpArgs.silentMode();
    inpArgs.silentMode(2);
    TMatrixD *covESFTot=NULL;

    // "default": eff corrected yield
    // this is almost postFSR cross section in the acceptance
    HistoPair2D_t hpUnfEff("hpUnfEff");
    if (res) {
      InputArgs_t iaNoRho(inpArgs,"rhoSyst");
      iaNoRho.needsEffScaleCorr(0);
      res=calculateCSdistribution(iaNoRho,hpSignalYield,
				  DYTools::_cs_postFsrDet,
				  hpUnfEff);
    }

    printHisto(hpUnfEff,5);

    for (int iSyst=0; res && (iSyst<2); ++iSyst) {
      int run=0;
      //DYTools::TSystematicsStudy_t runSystMode=DYTools::NO_SYST;
      TString csCovName;
      TString systTag;
      switch(iSyst) {
      case 0:
	run=calc_ESFtot;
	csCovName="covCS_ESFtot";
	break;
      case 1:
	run=calc_ESFtotCheck;
	csCovName="covCS_ESFtotCheck";
	break;
      default:
	std::cout << "not ready for eff iSyst=" << iSyst << "\n";
	return retCodeError;
      }

      if (!run) continue;
      std::cout << " - will produce " << csCovName << "\n";

      std::vector<TH2D*> vecRnd;

      if (iSyst==0) {
	// load the scale factors and their total errors
	int checkBinning=0;
	TString rhoCorrFName=inpMgr.correctionFullFileName("scale_factors",DYTools::NO_SYST,0);
	TH2D* hRho=LoadMatrixFields(rhoCorrFName,checkBinning,"scaleFactor","scaleFactorErr",1);
	if (!hRho) res=0;

	HistoPair2D_t hpRho("hpRho");
	if (res) {
	  hpRho.add(hRho,1.);
	  delete hRho;
	}
	else continue;

	std::cout << "Loaded "; printHisto(hpRho,6);

	int useSystErr=0; // randomize within statistical error
	if (!createRandomizedVec(hpRho,useSystErr,nExps,TString("hRnd_esf")+systTag,vecRnd)) {
	  std::cout << "failed to create randomized esemble for iSyst=" << iSyst << "\n";
	  return retCodeError;
	}
	hpRho.clear();
      }
      else if (iSyst==1) {
	// load the randomized vectors
	std::vector<TString> namesV; namesV.reserve(nExps);
	for (int iexp=0; res && (iexp<nExps); ++iexp) {
	  TString histoName=Form("h2RndScaleFactor_%dD_%d",DYTools::study2D+1,iexp);
	  namesV.push_back(histoName);
	  TH2D* h2=createBaseH2(histoName,histoName,1); // absRapidity
	  if (!h2) res=0;
	  vecRnd.push_back(h2);
	}
	if (res) {
	  TString fname=inpMgr.correctionSystFullFileName("scale_factors",DYTools::NO_SYST,0);
	  TFile inpF(fname,"read");
	  if (!inpF.IsOpen()) {
	    std::cout << "failed to open a file <" << inpF.GetName() << ">\n";
	    res=0;
	  }
	  if (res) res=loadVec(inpF,vecRnd,"rndSF");
	  inpF.Close();
	}
	if (vecRnd.size()<21) TCanvas *cx=plotProfiles("cx",vecRnd,namesV);
	//cx->
      }
      else {
	std::cout << "not ready for the ESF systematics iSyst=" << iSyst << "\n";
	return retCodeError;
      }
      if (!res) continue;

      if (nExps==2) printHisto(vecRnd,0,5,2);

      int unbiasedEstimate=1;
      TH2D* avgDistr=createBaseH2(Form("hESFAvgDistr_%d",iSyst));
      TH2D* csAvgDistr=createBaseH2(Form("hCSEsfAvgDistr_%d",iSyst));
      TMatrixD* cov= deriveCovMFromRndStudies(vecRnd,unbiasedEstimate,avgDistr);
      TMatrixD* csCov=NULL;

      if (cov) {
	if (!covESFTot) covESFTot=new TMatrixD(*cov);
	else (*covESFTot) += (*cov);
      }
      if (!cov) res=0;

      // The randomized efficiency scale values need to be applied
      // on the unfolded vector
      if (res) {
	HistoPair2D_t hpUnfEffSFYield("hpUnfEffSFYield");
	for (unsigned int i=0; res && (i<vecRnd.size()); ++i) {
	  TH2D *rndESF=vecRnd[i];
	  removeError(rndESF);
	  if (res) {
	    res=hpUnfEffSFYield.divide(hpUnfEff,rndESF);
	    //printHisto(hpUnfEffSFYield,6);
	  }
	  if (res) {
	    // swap pointers to histograms, since
	    // vecRnd has to contain ESF-corrected yields
	    hpUnfEffSFYield.swapHistoPtr(&vecRnd[i]);
	    // By default, the 1st entry will have name hpUnfEffSFYield,
	    // which is not correct.
	    // To avoid confusion, rename the histo
	    TString newName=Form("h2PostFsrYield_%d",i);
	    newName+=systTag;
	    vecRnd[i]->SetName(newName);
	    vecRnd[i]->SetTitle(newName);
	  }
	}
	hpUnfEffSFYield.clear();
      }

      //if (nExps==2) printHisto(vecRnd,0,5,2);

      if (res) {
	std::vector<TH2D*> csRndV;
	InputArgs_t inpArgsRho(inpArgs,"rhoSyst");
	inpArgsRho.needsDetUnfolding(0);
	inpArgsRho.needsEffCorr(0);
	inpArgsRho.needsEffScaleCorr(0);
	res=calcVecOfCSdistributions(inpArgsRho,vecRnd,csKind,csRndV);
	csCov= deriveCovMFromRndStudies(csRndV,unbiasedEstimate,csAvgDistr);
	if (!csCov) res=0;

	if (nExps==2) printHisto(csRndV,0,5,2);

	if (csCov) {
	  if (!covCS_fromESF) covCS_fromESF=new TMatrixD(*csCov);
	  else (*covCS_fromESF) += (*csCov);

	  outFile.cd();
	  csCov->Write(Form("covESF_iSyst_%d",iSyst));
	}

	//printHisto(csRndV,0,5,2);

	ClearVec(csRndV);
      }
      ClearVec(vecRnd);

      if (res) {
	outFile.cd();
	TString covName=csCovName;
	covName.ReplaceAll("covCS","cov");
	if (cov  ) cov->Write(covName);
	if (csCov) csCov->Write(csCovName);
      }

      // release memory
      if (csCov) delete csCov;
      if (cov  ) delete cov;
      if (csAvgDistr) delete csAvgDistr;
      if (avgDistr  ) delete avgDistr;
    }

    if (res) {
      outFile.cd();
      if (covESFTot) covESFTot->Write("covESFTot");
    }

    // clean-up
    hpUnfEff.clear();
    if (covESFTot) delete covESFTot;

    inpArgs.silentMode(saveSilentMode);
  }

  HERE("after calculating ESF cov res=%d",res);

  ////////////////////////////////
  // Uncertainties in MC acceptance
  ////////////////////////////////

  if (res && doCalcAccCov && DYTools::isFullSpaceCS(csKind)) {
    int saveSilentMode=inpArgs.silentMode();
    inpArgs.silentMode(2);
    TMatrixD *covAccTot=NULL;

    // "default" eff+effScale corrected yield
    // this is post-FSR cross section in the acceptance
    HistoPair2D_t hpPostFsrCS("hpPostFsrCS");
    if (res) {
      res=calculateCSdistribution(inpArgs,hpSignalYield,
				  DYTools::_cs_postFsrDet,
				  hpPostFsrCS);
    }

    printHisto(hpPostFsrCS,5);

    for (int iSyst=0; res && (iSyst<2); ++iSyst) {
      int run=0;
      DYTools::TSystematicsStudy_t runSystMode=DYTools::NO_SYST;
      TString csCovName;
      TString systTag;
      switch(iSyst) {
      case 0:
	run=calc_AccRnd;
	csCovName="covCS_AccRnd";
	break;
      case 1:
	run=calc_AccFSR;
	runSystMode=DYTools::FSR_STUDY; systTag="_FSRsyst";
	csCovName="covCS_AccFSR";
	break;
      default:
	std::cout << "not ready for acc iSyst=" << iSyst << "\n";
	return retCodeError;
      }

      if (!run) continue;
      std::cout << " - will produce " << csCovName << "\n";

      std::vector<TH2D*> vecRnd;
      if (res) {
	TString accField=TString("hacceptance") + systTag;
	TString hpCorrName=TString("hpAcc") + systTag;
	HistoPair2D_t hpAccCorr(accField);
	int checkBinning=0;
	int applyExtraTag=1;
	TString loadFileName= inpMgr.correctionSystFullFileName("acceptance",DYTools::NO_SYST,applyExtraTag);

	if (res) {
	  TH2D *h2corr=LoadHisto2D(accField,loadFileName,"",checkBinning);
	  if (!h2corr) {
	    std::cout << "failed to load systematics\n";
	    return retCodeError;
	  }
	  h2corr->SetName("h2tmp");
	  hpAccCorr.add(h2corr,1.);
	  delete h2corr;
	}
	std::cout << "Loaded "; printHisto(hpAccCorr,6);


	int useSystErr=0; // randomize within statistical error
	if (!createRandomizedVec(hpAccCorr,useSystErr,nExps,TString("hRnd_acc")+systTag,vecRnd)) {
	  std::cout << "failed to create randomized esemble for iSyst=" << iSyst << "\n";
	  return retCodeError;
	}
	hpAccCorr.clear();
      }

      printHisto(vecRnd,0,5,2);

      int unbiasedEstimate=1;
      TH2D* avgDistr=createBaseH2(Form("hAccAvgDistr_%d",iSyst));
      TH2D* csAvgDistr=createBaseH2(Form("hCSAccAvgDistr_%d",iSyst));
      TMatrixD* cov= deriveCovMFromRndStudies(vecRnd,unbiasedEstimate,avgDistr);
      TMatrixD* csCov=NULL;

      if (cov) {
	if (!covAccTot) covAccTot=new TMatrixD(*cov);
	else (*covAccTot) += (*cov);
      }
      if (!cov) res=0;

      // The randomized acceptance values need to be applied
      // on the post-FSR CS
      if (res) {
	HistoPair2D_t hpFullSpPostFsr("hpFullSpPostFsr");
	HistoPair2D_t hpFullSpPreFsr ("hpFullSpPreFsr");
	HistoPair2D_t *hpFinal=(csKind==DYTools::_cs_preFsr) ? &hpFullSpPreFsr : &hpFullSpPostFsr;
	TString newName;
	for (unsigned int i=0; res && (i<vecRnd.size()); ++i) {
	  TH2D *rndEff=vecRnd[i];
	  removeError(rndEff);
	  if (res) {
	    res=hpFullSpPostFsr.divide(hpPostFsrCS,rndEff);
	    //printHisto(hpFullSp,6);
	  }
	  if (csKind == DYTools::_cs_preFsr) {
	    // apply FSR unfolding locally
	    if (res) res=unfold_reco2true(hpFullSpPreFsr,*fsrResponse,hpFullSpPostFsr);
	    newName=Form("h2FullSpPreFsr_%d",i);
	  }
	  else {
	    newName=Form("h2FullSpPostFsr_%d",i);
	  }
	  if (res) {
	    // swap pointers to histograms, since
	    // vecRnd has to contain final yields
	    hpFinal->swapHistoPtr(&vecRnd[i]);

	    // To avoid confusion, rename the histo
	    newName+=systTag;
	    vecRnd[i]->SetName(newName);
	    vecRnd[i]->SetTitle(newName);
	  }
	}
	hpFullSpPostFsr.clear();
	hpFullSpPreFsr.clear();
      }

      printHisto(vecRnd,0,5,2);

      if (res) {
	// Cross section is already calculated, so we calculate the
	// covariance between the items of vecRnd
	csCov= deriveCovMFromRndStudies(vecRnd,unbiasedEstimate,csAvgDistr);
	if (!csCov) res=0;

	if (csCov) {
	  if (!covCS_fromAcc) covCS_fromAcc=new TMatrixD(*csCov);
	  else (*covCS_fromAcc) += (*csCov);
	}
      }
      ClearVec(vecRnd);

      if (res) {
	outFile.cd();
	TString covName=csCovName;
	covName.ReplaceAll("covCS","cov");
	if (cov  ) cov->Write(covName);
	if (csCov) csCov->Write(csCovName);
      }

      // release memory
      if (csCov) delete csCov;
      if (cov  ) delete cov;
      if (csAvgDistr) delete csAvgDistr;
      if (avgDistr  ) delete avgDistr;
    }

    if (res) {
      outFile.cd();
      if (covAccTot) covAccTot->Write("covAccTot");
    }

    // clean-up
    hpPostFsrCS.clear();
    if (covAccTot) delete covAccTot;

    inpArgs.silentMode(saveSilentMode);
  }

  HERE("after calculating acc cov res=%d",res);


  /////////////////////////////////////
  // Uncertainties in MC FSR unfolding
  /////////////////////////////////////

  if (doCalcFSRCov && DYTools::isPreFsrCS(csKind)) {
    int saveSilentMode=inpArgs.silentMode();
    inpArgs.silentMode(2);

    UnfoldingMatrix_t *fsrU=(DYTools::isFullSpaceCS(csKind)) ? fsrResponse : fsrResponseDet;
    UnfoldingMatrix_t *fsrUexact=(DYTools::isFullSpaceCS(csKind)) ? fsrResponseExact : fsrResponseDetExact;

    if (0) {
      fsrU->yieldsIni->Print();
      fsrU->yieldsFin->Print();
      //fsrUexact->yieldsIni->Print();
      //fsrUexact->yieldsFin->Print();
      return retCodeStop;
    }

    UnfoldingMatrix::TUnfoldingMatrixType_t fsrUnfKind=(DYTools::isFullSpaceCS(csKind)) ? UnfoldingMatrix::_cFSR : UnfoldingMatrix::_cFSR_DET;

    TMatrixD *covFsrTot=NULL;

    // this is post-FSR cross section (in acceptance or full space)
    HistoPair2D_t hpPostFsrCS("hpPostFsrCS");
    if (res) {
      DYTools::TCrossSectionKind_t iniCS=
	(DYTools::isFullSpaceCS(csKind)) ?
	             DYTools::_cs_postFsr : DYTools::_cs_postFsrDet;
      res=calculateCSdistribution(inpArgs,hpSignalYield,
				  iniCS,
				  hpPostFsrCS);
    }

    printHisto(hpPostFsrCS,5);

    for (int iSyst=0; res && (iSyst<2); ++iSyst) {
      int run=0;
      DYTools::TSystematicsStudy_t runSystMode=DYTools::NO_SYST;
      TString csCovName;
      TString systTag;
      // unf matrices for bounds
      TString uNamePlus, uNameMinus;
      DYTools::TSystematicsStudy_t smPlus=DYTools::NO_SYST;
      DYTools::TSystematicsStudy_t smMinus=smPlus;
      UnfoldingMatrix_t *Uplus=NULL, *Uminus=NULL;
      // unf matrix for the reference
      UnfoldingMatrix_t *Uref=NULL;
      switch(iSyst) {
      case 0:
	run=calc_FsrRnd;
	csCovName="covCS_FsrRnd";
	Uref=fsrU;
	break;
      case 1:
	run=calc_FsrFSR;
	runSystMode=DYTools::FSR_STUDY; systTag="_FSRsyst";
	uNamePlus ="fsrResponse_105"; smPlus =DYTools::FSR_5plus;
	uNameMinus="fsrResponse_095"; smMinus=DYTools::FSR_5minus;
	csCovName="covCS_FsrFSR";
	break;
      default:
	std::cout << "not ready for FSR iSyst=" << iSyst << "\n";
	return retCodeError;
      }

      if (!run) continue;
      std::cout << " - will produce " << csCovName << "\n";


      if (smPlus!=DYTools::NO_SYST) {
	Uplus =new UnfoldingMatrix_t(fsrUnfKind,uNamePlus);
	Uminus=new UnfoldingMatrix_t(fsrUnfKind,uNameMinus);
	if (!Uplus || !Uminus) res=0;
	if (res) res= Uplus->autoLoadFromFile(inpMgr.constDir(runSystMode,0),
				 UnfoldingMatrix_t::generateFNameTag(smPlus));
	if (res) res= Uminus->autoLoadFromFile(inpMgr.constDir(runSystMode,0),
				 UnfoldingMatrix_t::generateFNameTag(smMinus));
	if (res) {
	  Uref=new UnfoldingMatrix_t(*fsrU,"fsrResponseRef");
	  if (!Uref) res=0;
	  if (res && 0) {
	    Uref->yieldsIni->Print();
	    Uref->yieldsFin->Print();
	    return retCodeStop;
	  }
	}
	if (res) res=Uref->setErrorOnMigrationMatrix(*Uplus,*Uminus);
	if (res) {
	  if (Uplus) { delete Uplus; Uplus=NULL; }
	  if (Uminus) { delete Uminus; Uminus=NULL; }
	}
      }

      if (!res) {
	std::cout << "error in preparation the unfolding matrix\n";
	return retCodeError;
      }

      // The randomized unfolding matrices need to be applied
      // on the post-FSR CS
      std::vector<TH2D*> vecUnfRnd;
      if (res) {
	UnfoldingMatrix_t Urnd(fsrUnfKind,"fsrResponseRnd");
	HistoPair2D_t hpUnf("hpFsrUnf");

	for (int iexp=0; res && (iexp<nExps); ++iexp) {

	  res=Urnd.randomizeMigrationMatrix(*Uref, fsrUexact);
	  if (res) {
	    res=unfold_reco2true(hpUnf,Urnd,hpPostFsrCS);
	    printHisto(hpUnf,5);
	  }
	  if (res) {
	    TString rndVec=Form("h2fsrUnfRnd_%d",iexp);
	    TH2D* h2=Clone(hpUnf.histo(),rndVec);
	    if (!h2) res=0;
	    vecUnfRnd.push_back(h2);
	  }
	}
	hpUnf.clear(); // delete histos
      }

      if (nExps==2) printHisto(vecUnfRnd,0,5,2);

      int unbiasedEstimate=1;
      TH2D* csAvgDistr=createBaseH2(Form("hCSFsrAvgDistr_%d",iSyst));
      TMatrixD* csCov= deriveCovMFromRndStudies(vecUnfRnd,unbiasedEstimate,csAvgDistr);

      if (csCov) {
	if (!covFsrTot) covFsrTot=new TMatrixD(*csCov);
	else (*covFsrTot) += (*csCov);
	if (!covCS_fromFSR) covCS_fromFSR=new TMatrixD(*csCov);
	else (*covCS_fromFSR) += (*csCov);
      }
      if (!csCov) res=0;

      ClearVec(vecUnfRnd);

      if (res) {
	outFile.cd();
	if (csCov) csCov->Write(csCovName);
      }

      // release memory
      if (csCov) delete csCov;
      if (csAvgDistr) delete csAvgDistr;
    }

    if (res) {
      outFile.cd();
      if (covFsrTot) covFsrTot->Write("covFsrTot");
    }

    // clean-up
    hpPostFsrCS.clear();
    if (covFsrTot) delete covFsrTot;

    inpArgs.silentMode(saveSilentMode);
  }

  HERE("after calculating FSR cov res=%d",res);


  ////////////////////////////
  // Finalize
  ////////////////////////////

  if (res) {
    outFile.cd();
    if (covCS_fromYield) covCS_fromYield->Write("covCStot_fromYield");
    if (covCS_fromUnf  ) covCS_fromUnf  ->Write("covCStot_fromUnf");
    if (covCS_fromEff  ) covCS_fromEff  ->Write("covCStot_fromEff");
    if (covCS_fromESF  ) covCS_fromESF  ->Write("covCStot_fromESF");
    if (covCS_fromAcc  ) covCS_fromAcc  ->Write("covCStot_fromAcc");
    if (covCS_fromFSR  ) covCS_fromFSR  ->Write("covCStot_fromFSR");
  }
  else {
    std::cout << "ERROR: res=" << res << "\n";
  }


  outFile.Close();
  std::cout << "output file <" << outFile.GetName() << "> created\n";
  gBenchmark->Show("calcCSCov");
  return retCodeOk;

}
