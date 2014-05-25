// Check the results of the unfolding test
// The unfolding matrix is produced by splitting the sample
// into 2 pieces (SYST_RND)
// Created on May 19, 2014

#include "../Include/DYTools.hh"
#include "../Include/UnfoldingMatrix.h"
#include "../Include/HistoPair.hh"
#include "../Include/MyTools.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/EventSelector.hh"
#include <TRandom3.h>

//=== MAIN MACRO =================================================================================================

int checkClosure_systRnd(int analysisIs2D,
			 TString conf="defaultAdHoc",
			 int iSeed=1001,
			 int iBr=0,
			 int iSave=0,
			 int applyOnSigYield=0)
{
  DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN;
  DYTools::TSystematicsStudy_t systMode=DYTools::SYST_RND;

  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return retCodeError;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Settings
  //==============================================================================================================

  TString plotSavePath="plots-chk-systRnd";
  TString seedStr=Form("-seed%d",iSeed);

  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(conf)) return retCodeError;
  InputFileMgr_t *inpMgrForYield=NULL;
  if (applyOnSigYield) inpMgrForYield= new InputFileMgr_t(inpMgr);

  // Construct eventSelector, update mgr and plot directory
  TString extraTag,plotExtraTag;
  EventSelector_t *evtSelectorForYield=NULL;
  DYTools::TSystematicsStudy_t systModeForYield=DYTools::APPLY_ESCALE;
  if (inpMgrForYield) {
    evtSelectorForYield=new EventSelector_t(*inpMgrForYield,runMode,
					    systModeForYield,
					    extraTag,plotExtraTag,
					    EventSelector::_selectDefault);
  }
  EventSelector_t evtSelector(inpMgr,runMode,systMode,
			      extraTag, "", EventSelector::_selectDefault);
  evtSelector.setTriggerActsOnData(false);

  TString yieldName=(1 && (inpMgr.userKeyValueAsInt("DDBKG")==1)) ?
    "signalYieldDDbkg" : "signalYieldMCbkg";
  // signal yield
  HistoPair2D_t hpSignalYield(yieldName);

  // load signal yield
  const int loadNormalRunSelection=1;
  if (inpMgrForYield) {
    TString fnameBgSubtracted=
      inpMgrForYield->signalYieldFullFileName(systModeForYield,
					      loadNormalRunSelection);
    if (!hpSignalYield.Load(fnameBgSubtracted,1)) {
      std::cout << "failed to load signal yield\n";
      return retCodeError;
    }
  }

  // Prepare our unfolding matrix
  TString constDir= inpMgr.constDir(systMode,0);
  TString fnameTag= UnfoldingMatrix_t::generateFNameTag(systMode,iSeed);
  // global unf.M
  UnfoldingMatrix_t detResponse(UnfoldingMatrix::_cDET_Response,"detResponse");
  // first unf.M
  UnfoldingMatrix_t detResponseA(UnfoldingMatrix::_cDET_Response,
				 Form("detResponse_seed%d_replica%d",iSeed,0));
  // second unf.M, independent of the 1st one
  UnfoldingMatrix_t detResponseB(UnfoldingMatrix::_cDET_Response,
				 Form("detResponse_seed%d_replica%d",iSeed,1));
  if (!detResponse .autoLoadFromFile(constDir,fnameTag) ||
      !detResponseA.autoLoadFromFile(constDir,fnameTag) ||
      !detResponseB.autoLoadFromFile(constDir,fnameTag)) {
    std::cout << "failed to get the unfolding matrices\n";
    return retCodeError;
  }

  /*
  // unfolded yield
  HistoPair2D_t hpUnfoldedYield("unfYield");
  if (!hpUnfoldedYield.unfold(*detResponse.getDetInvResponse(),hpSignalYield)){
    std::cout << "failed to perform unfolding\n";
    return retCodeError;
    }*/

  //printHisto(hpUnfoldedYield);

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code
  //==============================================================================================================

  std::cout << mainpart;

  // compare initial and final MC yields
  if ((iBr==0) || (iBr==1) || (iBr==2) || (iBr==3)) {
    TColorRange_t colR=_colrange_default;
    int useMB=1; // use mass bins

    TH2D *h2Ini=createHisto2D(*detResponse.getIniM(),detResponse.getIniMerr(),
			      "h2Ini","h2Ini",colR,useMB);
    TH2D *h2IniA=
      createHisto2D(*detResponseA.getIniM(),detResponseA.getIniMerr(),
		    "h2IniA","h2IniA",colR,useMB);
    TH2D *h2IniB=
      createHisto2D(*detResponseB.getIniM(),detResponseB.getIniMerr(),
		    "h2IniB","h2IniB",colR,useMB);

    TH2D *h2Fin=createHisto2D(*detResponse.getFinM(),detResponse.getFinMerr(),
			      "h2Fin","h2Fin",colR,useMB);
    TH2D *h2FinA=
      createHisto2D(*detResponseA.getFinM(),detResponseA.getFinMerr(),
		    "h2FinA","h2FinA",colR,useMB);
    TH2D *h2FinB=
      createHisto2D(*detResponseB.getFinM(),detResponseB.getFinMerr(),
		    "h2FinB","h2FinB",colR,useMB);

    // scale full distributions
    double sumABini=h2IniA->Integral() + h2IniB->Integral();
    h2Ini->Scale(0.5*sumABini/h2Ini->Integral());
    double sumABfin=h2FinA->Integral() + h2FinB->Integral();
    h2Fin->Scale(0.5*sumABfin/h2Fin->Integral());

    if (1) {
      HERE("xx");
      //h2FinA->SetBinContent(41,1,1.452241);
      //h2FinA->SetBinContent(41,1,1.366610);
      printHisto(h2FinA);
      printHisto(h2FinB);
      printHisto(h2Fin);
    }

    std::vector<TH2D*> hRecoV;
    std::vector<TString> labelsRecoV;
    std::vector<ComparisonPlot_t*> cpRecoV;

    if (iBr<2) {
      hRecoV.push_back(h2Fin);
      labelsRecoV.push_back("reco (scaled)");
    }
    if (iBr!=3) { hRecoV.push_back(h2FinA); labelsRecoV.push_back("recoA"); }
    if (iBr!=2) { hRecoV.push_back(h2FinB); labelsRecoV.push_back("recoB"); }

    if (iBr==1) {
      hRecoV.push_back(unfold_true2reco(detResponseA,h2IniA,"mcA fold A"));
      labelsRecoV.push_back("mcA fold A");
      hRecoV.push_back(unfold_true2reco(detResponseB,h2IniB,"mcB fold B"));
      labelsRecoV.push_back("mcB fold B");
    }
    else if (iBr==2) {
      hRecoV.push_back(unfold_true2reco(detResponseB,h2IniA,"mcB fold A"));
      labelsRecoV.push_back("mcB fold A");
    }
    else if (iBr==3) {
      hRecoV.push_back(unfold_true2reco(detResponseA,h2IniB,"mcA fold B"));
      labelsRecoV.push_back("mcA fold B");
    }

    TCanvas *cxReco=plotProfiles("cxReco",hRecoV,labelsRecoV,NULL,0,
				 "yield",NULL,&cpRecoV,1);
    for (unsigned int i=0; i<cpRecoV.size(); ++i) {
      TString title=cpRecoV[i]->GetTitle();
      title.Append(Form(" seed=%d",iSeed));
      cpRecoV[i]->SetTitle(title);
      if (DYTools::study2D==0) {
	if (1) {
	  ComparisonPlot_t *cpp=cpRecoV[i];
	  cpp->SetPrintValues(1);
	  cpp->SetPrintRatios(1);
	}
	cpRecoV[i]->Draw(cxReco);
      }
      else {
	cpRecoV[i]->Draw6(cxReco,1,i+1);
      }
    }
    cxReco->Update();


    std::vector<TH2D*> hGenV;
    std::vector<TString> labelsGenV;
    std::vector<ComparisonPlot_t*> cpGenV;
    if (iBr<2) {
      hGenV.push_back(h2Ini); labelsGenV.push_back("gen (scaled)");
    }
    if (iBr!=3) { hGenV.push_back(h2IniA); labelsGenV.push_back("genA"); }
    if (iBr!=2) { hGenV.push_back(h2IniB); labelsGenV.push_back("genB"); }

    if (iBr==1) {
      hGenV.push_back(unfold_reco2true(detResponseA,h2FinA,"mcA unf A"));
      labelsGenV.push_back("mcA unf A");
      hGenV.push_back(unfold_reco2true(detResponseB,h2FinB,"mcB unf B"));
      labelsGenV.push_back("mcB unf B");
    }
    else if (iBr==2) {
      hGenV.push_back(unfold_reco2true(detResponseB,h2FinA,"mcB unf A"));
      labelsGenV.push_back("mcB unf A");
    }
    else if (iBr==3) {
      hGenV.push_back(unfold_reco2true(detResponseA,h2FinB,"mcA unf B"));
      labelsGenV.push_back("mcA unf B");
    }

   TCanvas *cxGen=plotProfiles("cxGen",hGenV,labelsGenV,NULL,0,
				"yield",NULL,&cpGenV,1);

    for (unsigned int i=0; i<cpGenV.size(); ++i) {
      TString title=cpGenV[i]->GetTitle();
      title.Append(Form(" seed=%d",iSeed));
      cpGenV[i]->SetTitle(title);
      if (DYTools::study2D==0) {
	cpGenV[i]->Draw(cxGen);
      }
      else {
	cpGenV[i]->Draw6(cxGen,1,i+1);
      }
    }
    cxGen->Update();

    TString fnameBase=Form("fig-%s-mcYields",DYTools::analysisTag.Data());
    fnameBase.Append(seedStr);
    fnameBase.Append("-");
    TString fnameReco=fnameBase + TString("Reco");
    TString fnameGen =fnameBase + TString("Gen");
    if (iBr==1) {
      fnameReco.Append("Closure");
      fnameGen.Append("Closure");
    }
    else if (iBr==2) {
      fnameReco.Append("ChkA");
      fnameGen.Append("ChkA");
    }
    else if (iBr==3) {
      fnameReco.Append("ChkB");
      fnameGen.Append("ChkB");
    }
    std::cout << "FName=<" << fnameReco << "> and <" << fnameGen << ">\n";
    std::cout << " at path=<" << plotSavePath << ">\n";
    if (iSave) {
      SaveCanvas(cxReco,fnameReco,plotSavePath);
      SaveCanvas(cxGen ,fnameGen ,plotSavePath);
    }
    else {
      std::cout << "Not saved as requested\n";
    }
  }


  // ---------------------------------------
  //  create pull plots
  // ---------------------------------------

  if ((iBr==10)) {
    TColorRange_t colR=_colrange_default;
    int useMB=1; // use mass bins
    gRandom->SetSeed(iSeed);

    // gen
    TH2D *h2Ini=createHisto2D(*detResponse.getIniM(),detResponse.getIniMerr(),
			      "h2Ini","h2Ini",colR,useMB);

    // reco
    TH2D *h2Fin=createHisto2D(*detResponse.getFinM(),detResponse.getFinMerr(),
			      "h2Fin","h2Fin",colR,useMB);

     HistoPair2D_t hpReco("hpReco");
    if (!hpReco.assign(h2Fin,NULL)) {
      std::cout << "failed assigning to hpReco\n";
      return retCodeError;
    }

    int nExps=1000;
    unsigned int nExps100=100;
    std::vector<TH2D*> histoRecoV, histoUnfV;
    if (!createRandomizedVec(hpReco,0,nExps,"hRecoRnd_",histoRecoV)) {
      std::cout << "failed randomization\n";
      return retCodeError;
    }

    TH2D *h2Sum=Clone(h2Ini,"h2Sum","h2Sum");
    h2Sum->Reset();
    TH2D *h2Sum100=Clone(h2Ini,"h2Sum100","h2Sum100");
    h2Sum100->Reset();
    histoUnfV.reserve(histoRecoV.size());

    for (unsigned int i=0; i<histoRecoV.size(); ++i) {
      TString name=Form("hUnf_%d",i);
      TH2D *hUnf=unfold_reco2true(detResponse,histoRecoV[i],name);
      histoUnfV.push_back(hUnf);
      accumulateForRndStudies(h2Sum,hUnf);
      if (i<nExps100) accumulateForRndStudies(h2Sum100,hUnf);
    }
    accumulateForRndStudies_finalize(h2Sum,nExps,0);
    accumulateForRndStudies_finalize(h2Sum100,nExps100,0);
    printHisto(h2Sum);

    TH2D* h2Pull=Clone(h2Ini,"h2Pull","h2Pull");
    h2Pull->Reset();
    TH2D* h2Pull100=Clone(h2Ini,"h2Pull100","h2Pull100");
    h2Pull100->Reset();
    swapContentAndError2D(h2Sum);
    swapContentAndError2D(h2Sum100);
    for (unsigned int i=0; i<histoUnfV.size(); ++i) {
      TString name=Form("hDiff_%d",i);
      TH2D *hTmp=Clone(histoUnfV[i],name);
      hTmp->Add(h2Ini,-1);
      if (i<nExps100) {
	TH2D* hTmp100=Clone(hTmp,"hTmp100");
	scaleHisto(hTmp100,h2Sum100);
	accumulateForRndStudies(h2Pull100,hTmp100);
	delete hTmp100;
      }
      scaleHisto(hTmp,h2Sum);
      accumulateForRndStudies(h2Pull,hTmp);
      delete hTmp;
    }
    accumulateForRndStudies_finalize(h2Pull,nExps,0);
    accumulateForRndStudies_finalize(h2Pull100,nExps100,0);
    swapContentAndError2D(h2Sum); // restore
    swapContentAndError2D(h2Sum100); // restore

    std::vector<TH2D*> histoV;
    std::vector<TString> labelsV;
    std::vector<ComparisonPlot_t*> cpV;

    histoV.push_back(h2Sum); labelsV.push_back("unfolded avg");
    histoV.push_back(h2Ini); labelsV.push_back("gen");

    TCanvas *cx=plotProfiles("cx",histoV,labelsV,NULL,0,
			     "yield",NULL,&cpV,1);

    for (unsigned int i=0; i<cpV.size(); ++i) {
      TString title=cpV[i]->GetTitle();
      title.Append(Form(" seed=%d",iSeed));
      cpV[i]->SetTitle(title);
      if (DYTools::study2D==0) {
	if (0) {
	  ComparisonPlot_t *cpp=cpV[i];
	  cpp->SetPrintValues(1);
	  cpp->SetPrintRatios(1);
	}
	cpV[i]->Draw(cx);
      }
      else {
	cpV[i]->Draw6(cx,1,i+1);
      }
    }
    cx->Update();


    std::vector<TH2D*> hPullV;
    std::vector<TString> labelsPullV;
    std::vector<ComparisonPlot_t*> cpPullV;

    hPullV.push_back(h2Pull);
    labelsPullV.push_back(Form("pull (nExps=%d)",nExps));
    hPullV.push_back(h2Pull100);
    labelsPullV.push_back(Form("pull (nExps=%d)",nExps100));

    TCanvas *cxPull=plotProfiles("cxPull",hPullV,labelsPullV,NULL,0,
				 "pull",NULL,&cpPullV,1);

    for (unsigned int i=0; i<cpPullV.size(); ++i) {
      TString title=cpPullV[i]->GetTitle();
      title.Append(Form(" seed=%d",iSeed));
      cpPullV[i]->SetTitle(title);
      cpPullV[i]->ErrorsOnRatios(0);
      if (DYTools::study2D==0) {
	cpPullV[i]->Draw(cxPull);
      }
      else {
	cpPullV[i]->Draw6(cxPull,1,i+1);
      }
    }
    cxPull->Update();

    std::vector<TH2D*> hPullErrV;
    std::vector<TString> labelsPullErrV;
    std::vector<ComparisonPlot_t*> cpPullErrV;

    TH2D *h2PullErr=Clone(h2Pull,"h2PullError");
    swapContentAndError(h2PullErr);
    removeError(h2PullErr);
    TH2D *h2Pull100Err=Clone(h2Pull100,"h2Pull100Error");
    swapContentAndError(h2Pull100Err);
    removeError(h2Pull100Err);

    hPullErrV.push_back(h2PullErr);
    labelsPullErrV.push_back(Form("nExps=%d",nExps));
    hPullErrV.push_back(h2Pull100Err);
    labelsPullErrV.push_back(Form("nExps=%d",nExps100));

    TCanvas *cxPullErr=plotProfiles("cxPullErr",hPullErrV,labelsPullErrV,
				    NULL,0,
				    "pull error",NULL,&cpPullErrV,1);

    for (unsigned int i=0; i<cpPullErrV.size(); ++i) {
      TString title=cpPullErrV[i]->GetTitle();
      title.Append(Form(" seed=%d",iSeed));
      cpPullErrV[i]->SetTitle(title);
      cpPullErrV[i]->ErrorsOnRatios(0);
      if (DYTools::study2D==0) {
	cpPullErrV[i]->Draw(cxPullErr);
      }
      else {
	cpPullErrV[i]->Draw6(cxPullErr,1,i+1);
      }
    }
    cxPullErr->Update();

    TString fnameBase=Form("fig-%s-mcYields",DYTools::analysisTag.Data());
    fnameBase.Append(seedStr);
    fnameBase.Append("-");
    TString fname    =fnameBase + TString("AvgUnf");
    TString fnamePull =fnameBase + TString("Pull");
    TString fnamePullErr =fnameBase + TString("PullErr");
    std::cout << "FName=<" << fname << ">, <" << fnamePull
	      << "> and <" << fnamePullErr << ">\n";
    std::cout << " at path=<" << plotSavePath << ">\n";
    if (iSave) {
      SaveCanvas(cx    ,fname    ,plotSavePath);
      SaveCanvas(cxPull,fnamePull,plotSavePath);
      SaveCanvas(cxPullErr,fnamePullErr,plotSavePath);
    }
    else {
      std::cout << "Not saved as requested\n";
    }
  }

  // ---------------------------------------
  // check Alexey
  // ---------------------------------------

  if ((iBr==20)) {
    TColorRange_t colR=_colrange_default;
    int useMB=1; // use mass bins
    gRandom->SetSeed(iSeed);

    // gen
    TH2D *h2Ini=createHisto2D(*detResponse.getIniM(),detResponse.getIniMerr(),
			      "h2Ini","h2Ini",colR,useMB);

    // reco
    TH2D *h2Fin=createHisto2D(*detResponse.getFinM(),detResponse.getFinMerr(),
			      "h2Fin","h2Fin",colR,useMB);

     HistoPair2D_t hpReco("hpReco");
     HistoPair2D_t hpGen("hpGen");
    if (!hpReco.assign(h2Fin,NULL) ||
	!hpGen .assign(h2Ini,NULL)) {
      std::cout << "failed assigning to hpReco or hpGen\n";
      return retCodeError;
    }

    int nExps=1000;
    unsigned int nExps100=100;
    std::vector<TH2D*> histoRecoV, histoUnfV;
    std::vector<TH2D*> histoGenV;
    std::vector<TString> histoGenLabelsV;

    if (!createRandomizedVec(hpReco,0,nExps,"hRecoRnd_",histoRecoV) ||
	!createRandomizedVec(hpGen ,0,nExps,"hGenRnd_" ,histoGenV,&histoGenLabelsV)) {
      std::cout << "failed randomization\n";
      return retCodeError;
    }

    if (0) {
      TCanvas *cy=plotProfiles("cy",histoGenV,histoGenLabelsV,NULL,0,
			       "gen yield",NULL,NULL,0);
      cy->Update();
      return retCodeStop;
    }

    TH2D *h2RecoSum=Clone(h2Ini,"h2RecoSum","h2RecoSum");
    h2RecoSum->Reset();
    TH2D *h2Sum=Clone(h2Ini,"h2Sum","h2Sum");
    h2Sum->Reset();
    TH2D *h2Sum100=Clone(h2Ini,"h2Sum100","h2Sum100");
    h2Sum100->Reset();
    histoUnfV.reserve(histoRecoV.size());

    for (unsigned int i=0; i<histoRecoV.size(); ++i) {
      accumulateForRndStudies(h2RecoSum,histoRecoV[i]);
      TString name=Form("hUnf_%d",i);
      TH2D *hUnf=unfold_reco2true(detResponse,histoRecoV[i],name);
      histoUnfV.push_back(hUnf);
      accumulateForRndStudies(h2Sum,hUnf);
      if (i<nExps100) accumulateForRndStudies(h2Sum100,hUnf);
    }
    accumulateForRndStudies_finalize(h2RecoSum,nExps,0);
    accumulateForRndStudies_finalize(h2Sum,nExps,0);
    accumulateForRndStudies_finalize(h2Sum100,nExps100,0);
    printHisto(h2Sum);

    if (1) {
      TH2D *h2RecoSumErr=Clone(h2RecoSum,"h2RecoSumErr");
      swapContentAndError(h2RecoSumErr);
      removeError(h2RecoSumErr);
      TH2D *h2RecoErr=Clone(h2Fin,"h2RecoErr");
      swapContentAndError(h2RecoErr);
      removeError(h2RecoErr);
      std::vector<TH2D*> hVchk;
      std::vector<TString> labelsVchk;
      hVchk.push_back(h2RecoSumErr); labelsVchk.push_back("ensemble err");
      hVchk.push_back(h2RecoErr); labelsVchk.push_back("recoErr");
      TCanvas *cz= plotProfiles("cz",hVchk,labelsVchk,NULL,0,
				"error",NULL,NULL,0);
      cz->Update();
      return retCodeStop;
    }
    if (1) {
      TH2D *h2SumErr=Clone(h2Sum,"h2SumErr");
      swapContentAndError(h2SumErr);
      removeError(h2SumErr);
      TH2D *h2GenErr=Clone(h2Ini,"h2GenErr");
      swapContentAndError(h2GenErr);
      removeError(h2GenErr);
      std::vector<TH2D*> hVchk;
      std::vector<TString> labelsVchk;
      hVchk.push_back(h2SumErr); labelsVchk.push_back("unf.ensemble err");
      hVchk.push_back(h2GenErr); labelsVchk.push_back("genErr");
      TCanvas *cz= plotProfiles("cz",hVchk,labelsVchk,NULL,0,
				"error",NULL,NULL,0);
      cz->Update();
      return retCodeStop;
    }

    //delete h2Sum;
    //h2Sum=Clone(h2Ini,"h2Sum_chk");

    TH2D* h2Pull=Clone(h2Ini,"h2Pull","h2Pull");
    h2Pull->Reset();
    TH2D* h2Pull100=Clone(h2Ini,"h2Pull100","h2Pull100");
    h2Pull100->Reset();
    swapContentAndError2D(h2Sum);
    swapContentAndError2D(h2Sum100);
    for (unsigned int i=0; i<histoUnfV.size(); ++i) {
      TString name=Form("hDiff_%d",i);
      TH2D *hTmp=Clone(histoUnfV[i],name);
      hTmp->Add(histoGenV[i],-1);
      if (i<nExps100) {
	TH2D* hTmp100=Clone(hTmp,"hTmp100");
	scaleHisto(hTmp100,h2Sum100);
	accumulateForRndStudies(h2Pull100,hTmp100);
	delete hTmp100;
      }
      scaleHisto(hTmp,h2Sum);
      accumulateForRndStudies(h2Pull,hTmp);
      delete hTmp;
    }
    accumulateForRndStudies_finalize(h2Pull,nExps,0);
    accumulateForRndStudies_finalize(h2Pull100,nExps100,0);
    swapContentAndError2D(h2Sum); // restore
    swapContentAndError2D(h2Sum100); // restore

    std::vector<TH2D*> histoV;
    std::vector<TString> labelsV;
    std::vector<ComparisonPlot_t*> cpV;

    histoV.push_back(h2Sum); labelsV.push_back("unfolded avg");
    histoV.push_back(h2Ini); labelsV.push_back("gen");

    TCanvas *cx=plotProfiles("cx",histoV,labelsV,NULL,0,
			     "yield",NULL,&cpV,1);

    for (unsigned int i=0; i<cpV.size(); ++i) {
      TString title=cpV[i]->GetTitle();
      title.Append(Form(" seed=%d",iSeed));
      cpV[i]->SetTitle(title);
      if (DYTools::study2D==0) {
	if (0) {
	  ComparisonPlot_t *cpp=cpV[i];
	  cpp->SetPrintValues(1);
	  cpp->SetPrintRatios(1);
	}
	cpV[i]->Draw(cx);
      }
      else {
	cpV[i]->Draw6(cx,1,i+1);
      }
    }
    cx->Update();


    std::vector<TH2D*> hPullV;
    std::vector<TString> labelsPullV;
    std::vector<ComparisonPlot_t*> cpPullV;

    hPullV.push_back(h2Pull);
    labelsPullV.push_back(Form("pull (nExps=%d)",nExps));
    hPullV.push_back(h2Pull100);
    labelsPullV.push_back(Form("pull (nExps=%d)",nExps100));

    TCanvas *cxPull=plotProfiles("cxPull",hPullV,labelsPullV,NULL,0,
				 "pull",NULL,&cpPullV,1);

    for (unsigned int i=0; i<cpPullV.size(); ++i) {
      TString title=cpPullV[i]->GetTitle();
      title.Append(Form(" seed=%d",iSeed));
      cpPullV[i]->SetTitle(title);
      cpPullV[i]->ErrorsOnRatios(0);
      if (DYTools::study2D==0) {
	cpPullV[i]->Draw(cxPull);
      }
      else {
	cpPullV[i]->Draw6(cxPull,1,i+1);
      }
    }
    cxPull->Update();

    std::vector<TH2D*> hPullErrV;
    std::vector<TString> labelsPullErrV;
    std::vector<ComparisonPlot_t*> cpPullErrV;

    TH2D *h2PullErr=Clone(h2Pull,"h2PullError");
    swapContentAndError(h2PullErr);
    removeError(h2PullErr);
    TH2D *h2Pull100Err=Clone(h2Pull100,"h2Pull100Error");
    swapContentAndError(h2Pull100Err);
    removeError(h2Pull100Err);

    hPullErrV.push_back(h2PullErr);
    labelsPullErrV.push_back(Form("nExps=%d",nExps));
    hPullErrV.push_back(h2Pull100Err);
    labelsPullErrV.push_back(Form("nExps=%d",nExps100));

    TCanvas *cxPullErr=plotProfiles("cxPullErr",hPullErrV,labelsPullErrV,
				    NULL,0,
				    "pull error",NULL,&cpPullErrV,1);

    for (unsigned int i=0; i<cpPullErrV.size(); ++i) {
      TString title=cpPullErrV[i]->GetTitle();
      title.Append(Form(" seed=%d",iSeed));
      cpPullErrV[i]->SetTitle(title);
      cpPullErrV[i]->ErrorsOnRatios(0);
      if (DYTools::study2D==0) {
	cpPullErrV[i]->Draw(cxPullErr);
      }
      else {
	cpPullErrV[i]->Draw6(cxPullErr,1,i+1);
      }
    }
    cxPullErr->Update();

    TString fnameBase=Form("fig-%s-mcYields",DYTools::analysisTag.Data());
    fnameBase.Append(seedStr);
    fnameBase.Append("-");
    TString fname    =fnameBase + TString("AvgUnf");
    TString fnamePull =fnameBase + TString("Pull");
    TString fnamePullErr =fnameBase + TString("PullErr");
    std::cout << "FName=<" << fname << ">, <" << fnamePull
	      << "> and <" << fnamePullErr << ">\n";
    std::cout << " at path=<" << plotSavePath << ">\n";
    if (iSave) {
      SaveCanvas(cx    ,fname    ,plotSavePath);
      SaveCanvas(cxPull,fnamePull,plotSavePath);
      SaveCanvas(cxPullErr,fnamePullErr,plotSavePath);
    }
    else {
      std::cout << "Not saved as requested\n";
    }
  }

  // ---------------------------------------
  //
  // ---------------------------------------




  return retCodeOk;
}


