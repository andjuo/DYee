#include "../Include/DYTools.hh"
#include "../CrossSection/crossSectionFnc.hh"
#include "../Include/colorPalettes.hh"
#include <TBenchmark.h>


//=== Global flags =================================================================================================


// Function plots randomized yields
//   iSyst=0, statError
//   iSyst=1, systError
//   iSyst=2, statError and systError
//   iSyst=3, statError randomized vectors
//   iSyst=4, systError randomized vectors
// For iSyst=3 and 4, the error plots have no meaning
int plotSignalYieldRndResults(TString fname, TString histoDirName,
			      TString yieldName, int iSyst, int saveCanvas,
			      TString extraTagForSaving,
			      TH2D **hp=NULL);

// Plot errors from signal yields (propagated and from cov)
int plotSignalYieldCSResults(TString fname, TString csFieldName,
			     int saveCanvas, TString extraTagForSaving);


//=== MAIN MACRO =================================================================================================

int plotCSCovStudy(int analysisIs2D,
		   TString conf,
		   int iBr=0,
		   int saveCanvas=0,
		   DYTools::TCrossSectionKind_t csKind=DYTools::_cs_None,
		   TString outFileExtraTag_UserInput="") {

  gBenchmark->Start("plotCSCovStudy");

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

  DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN;
  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;
  //DYTools::TSystematicsStudy_t yieldSystMode=DYTools::APPLY_ESCALE;

  //InputFileMgr_t inpMgrEScale;
  //if (!inpMgrEScale.Load(conf)) return retCodeError;

  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(conf)) return retCodeError;

  // Construct eventSelector, update mgr and plot directory
  TString extraTag; // empty
  TString plotExtraTag;

  //EventSelector_t evtSelectorEScale(inpMgrEScale,runMode,yieldSystMode,
  //			      extraTag,plotExtraTag,
  //			      EventSelector::_selectDefault);
  EventSelector_t evtSelector(inpMgr,runMode,systMode,
			      extraTag,plotExtraTag,
			      EventSelector::_selectDefault);

  int systFileFlag=1;
  TString outFileName=inpMgr.crossSectionFullFileName(systMode,
					       csKind,0,systFileFlag);
  if (outFileExtraTag_UserInput.Length()) {
    std::cout << "applying extra tag from input= " << outFileExtraTag_UserInput << "\n";
    outFileName.ReplaceAll(".root",outFileExtraTag_UserInput);
    if (outFileName.Index(".root")==-1) outFileName.Append(".root");
  }
  std::cout << "outFileName=<" << outFileName << ">\n";

  // -----------------------
  //   Additional settings
  // -----------------------

  int useDDBkg=(inpMgr.userKeyValueAsInt("DDBKG")==1) ? 1:0;
  TString yieldName=(useDDBkg) ? "signalYieldDDbkg" : "signalYieldMCbkg";
  HistoPair2D_t hpSignalYield(yieldName);

  TString resCSName=yieldName;
  resCSName.ReplaceAll("signal","cs");
  HistoPair2D_t hpCS(resCSName);

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code
  //==============================================================================================================

  int res=1;
  if (iBr==0) {
    res=plotSignalYieldRndResults(outFileName,"covCS_YieldStat_details",
				  yieldName,0,saveCanvas,
				  outFileExtraTag_UserInput);
  }
  else if (iBr==1) {
    res=plotSignalYieldRndResults(outFileName,"covCS_YieldSyst_details",
				  yieldName,1,saveCanvas,
				  outFileExtraTag_UserInput);
  }
  else if (iBr==2) {
    TH2D* h2YieldFromStatErr=NULL;
    res=plotSignalYieldRndResults(outFileName,"covCS_YieldStat_details",
				  yieldName,0,0,
				  outFileExtraTag_UserInput,
				  &h2YieldFromStatErr);
    res=plotSignalYieldRndResults(outFileName,"covCS_YieldSyst_details",
				  yieldName,2,saveCanvas,
				  outFileExtraTag_UserInput,
				  &h2YieldFromStatErr);
  }
  else if (iBr==3) {
    res=plotSignalYieldRndResults(outFileName,"covCS_YieldStat_details",
				  yieldName,3,saveCanvas,
				  outFileExtraTag_UserInput);
  }
  else if (iBr==4) {
    res=plotSignalYieldRndResults(outFileName,"covCS_YieldSyst_details",
				  yieldName,4,saveCanvas,
				  outFileExtraTag_UserInput);
  }
  else if (iBr==5) {
    res=plotSignalYieldCSResults(outFileName,"hpFinal_yieldErrOnly",
				  saveCanvas,
				  outFileExtraTag_UserInput);
  }

  ShowBenchmarkTime("plotCSCovStudy");

  return retCodeOk;
}


// ----------------------------------------------------------------------------
//     Function definitions
//============================================================================

int plotSignalYieldRndResults(TString fname, TString histoDirName,
			      TString yieldName, int iSyst,
			      int saveCanvas, TString extraTagForSaving,
			      TH2D **hp) {

  const char *fncName="plotSignalYieldRndResults";
  HistoPair2D_t hpYield(yieldName);
  TString hAvgName=Form("hYieldAvgDistr_%d",(iSyst<2) ? iSyst : 1);
  if (iSyst>2) hAvgName=Form("hYieldAvgDistr_%d", iSyst-3);
  TH2D* hAvg=createBaseH2(hAvgName);

  std::cout << "plotSignalYieldRndResults extraTagForSaving_inp=<"
	    << extraTagForSaving << ">, ";
  extraTagForSaving.ReplaceAll("-yieldOnly","");
  std::cout << " extraTag_final=<" << extraTagForSaving << "\n";

  int res=1;
  TFile fin(fname,"read");
  res=fin.IsOpen();
  if (!res) {
    std::cout << fncName << ": failed to open the file <" << fname << ">\n";
    return 0;
  }
  res = hpYield.Read(fin,histoDirName,"");
  if (res) {
    hAvg=LoadHisto2D(fin,hAvg->GetName(),histoDirName,1);
    if (!hAvg) res=0;
    if (res && hp && (iSyst!=2)) {
      (*hp)=Clone(hAvg,TString(hAvg->GetName()) + TString("_clone"),"clone");
      if (!(*hp)) res=0;
    }
  }
  fin.Close();

  if (!res) {
    std::cout << "failed to get the average distribution\n";
    return 0;
  }

  TH2D* h2YieldSystErr=NULL;

  std::vector<TH2D*> histoV;
  std::vector<TH2D*> histoErrV;
  std::vector<TString> labelV;
  std::vector<TString> labelErrV;
  TString fileTagBase, fileTag;
  histoV.reserve(2);
  labelV.reserve(2);

  if ((iSyst==0) || (iSyst==3)) {
    histoV.push_back(hpYield.histo());
    labelV.push_back("signal yield w/stat.err");
    labelErrV.push_back("stat.err. of signal yield");
    fileTagBase="statErr";
  }
  else if ((iSyst==1) || (iSyst==2) || (iSyst==4)) {
    h2YieldSystErr=hpYield.createHistoWithSystError(Form("hSystErr_%d",iSyst));
    histoV.push_back(h2YieldSystErr);
    labelV.push_back("signal yield w/syst.err");
    labelErrV.push_back("syst.err. of signal yield");
    fileTagBase="systErr";
  }

  if (iSyst!=2) {
    histoV.push_back(hAvg);
    labelV.push_back("randomization average");
    labelErrV.push_back("error of randomization");
  }
  else {
    if (!hp) { std::cout << "iSyst=2 requires hp!=NULL\n"; return 0; }
    histoV.push_back((*hp));
    labelV.push_back("stat.err randomization");
    labelErrV.push_back("stat.err from randomization");
    histoV.push_back(hAvg);
    labelV.push_back("syst.err randomization");
    labelErrV.push_back("syst.err from randomization");
  }

  if ((iSyst==3) || (iSyst==4)) {
    TFile fin2(fin.GetName(),"read");
    if (!fin2.IsOpen()) {
      std::cout << "failed to re-open the file\n";
      return 0;
    }

    std::vector<TH2D*> histoTmpV;
    std::vector<TString> sampleLabelsV;
    int nExps=100;
    sampleLabelsV.reserve(nExps);
    for (int i=0; i<nExps; ++i) {
      sampleLabelsV.push_back(Form("%d",i));
    }

    TString histoNameBase=(iSyst==3) ? "hRnd_yield_stat_" : "hRnd_yield_syst_";
    TString histoSubDir=(iSyst==3) ? "covCS_YieldStat_details" : "covCS_YieldSyst_details";
    res=(createBaseH2Vec(histoTmpV,histoNameBase,sampleLabelsV,1,1) &&
	 loadVec(fin2,histoTmpV,histoSubDir)) ? 1:0;
    fin2.Close();
    if (!res) {
      std::cout << "error preparing the histos\n";
      return 0;
    }

    if (res) {
      histoV.reserve(histoV.size() + nExps);
      labelV.reserve(labelV.size() + nExps);
      labelErrV.reserve(labelErrV.size() + nExps);
      for (unsigned int i=0; i<histoTmpV.size(); ++i) {
	histoV.push_back(histoTmpV[i]);
	labelV.push_back(histoTmpV[i]->GetName());
	labelErrV.push_back(histoTmpV[i]->GetName());
      }
    }
  }

  for (int plotErr=0; plotErr<2; ++plotErr) {
    TString canvName=Form("cx%d",plotErr);
    std::vector<std::vector<TH1D*>*> hProfV;
    std::vector<ComparisonPlot_t*> cpV;
    TCanvas *cx=NULL;


    if (plotErr==0) {
      fileTag = fileTagBase + TString("-yields");
      cx= plotProfiles(canvName,histoV,labelV,NULL,0,"yield counts",
		       &hProfV,&cpV,1);
    }
    else if (plotErr==1) {
      fileTag = fileTagBase + TString("-errOnly");
      histoErrV.reserve(histoV.size());
      for (unsigned int i=0; i<histoV.size(); ++i) {
	TH2D* h2src=histoV[i];
	TString newName=TString(h2src->GetName()) + TString("_err");
	TH2D* h2=Clone(h2src,newName,newName);
	swapContentAndError2D(h2);
	removeError2D(h2);
	histoErrV.push_back(h2);
      }
      cx= plotProfiles(canvName,histoErrV,labelErrV,
		       NULL,0,"yield error",&hProfV,&cpV,1);
    }

    std::cout << "cpV.size()=" << cpV.size() << "\n";

    if ((iSyst==3) || (iSyst==4)) {
      for (unsigned int i=0; i<cpV.size(); ++i) {
	ComparisonPlot_t *cp=cpV[i];
	for (int ii=0; ii<2; ii++) {
	  TH1D* hSrc=cp->GetHisto(ii);
	  TH1D* h=(TH1D*)hSrc->Clone(TString(hSrc->GetName()) + TString("_clone"));
	  h->SetDirectory(0);
	  h->SetMarkerStyle(33);
	  int color=(ii==0) ? kWhite : (kRed+1);
	  cp->AddHist1D(h,labelV[ii],"LP",color,0,0,-1);
	}
      }
    }

    if (cpV.size()==1) {
      cpV[0]->ChangeLegendPos(-0.15,0.0,-0.08,0.);
      if (1) {
	cpV[0]->SetLogy(1);
	cpV[0]->TransLegend(-0.2,-0.65);
      }
      cpV[0]->Draw(cx);
    }
    else {
      for (unsigned int i=0; i<cpV.size(); ++i) {
	cpV[i]->ChangeLegendPos(-0.2,0.0,-0.08,0.);
	if (1) {
	  cpV[i]->TransLegend(-0.2,-0.65);
	}
	cpV[i]->Draw6(cx,1,i+1);
      }
    }

    if (cx) cx->Update();

    if (iSyst==2) fileTag.Append("CombiPlot");
    TString str2D=Form("_%dD",DYTools::study2D+1);
    TString dirName=Form("plots_csCovStudy%s",str2D.Data());
    TString figName=Form("fig-%s-",yieldName.Data()) +
                         fileTag + str2D + extraTagForSaving;
    std::cout << "figName=<" << figName << ">\n";
    std::cout << "dirName=<" << dirName << ">\n";
    if (saveCanvas) {
      SaveCanvas(cx,figName,dirName);
    }
    else std::cout << " .. not saved as requested\n";
  }

  return 1;
}

// ------------------------------------------------------------------

int plotSignalYieldCSResults(TString fname, TString csFieldName,
			     int saveCanvas, TString extraTagForSaving) {

  const char *fncName="plotSignalYieldCSResults";
  HistoPair2D_t hpFinCS(csFieldName);

  std::cout << "plotSignalYieldRndResults extraTagForSaving_inp=<"
	    << extraTagForSaving << ">, ";
  extraTagForSaving.ReplaceAll("-yieldOnly","");
  std::cout << " extraTag_final=<" << extraTagForSaving << "\n";

  int res=1;
  TFile fin(fname,"read");
  res=fin.IsOpen();
  if (!res) {
    std::cout << fncName << ": failed to open the file <" << fname << ">\n";
    return 0;
  }
  res = hpFinCS.Read(fin,"","");
  TMatrixD *covCS_YieldStat=(TMatrixD*)fin.Get("covCS_YieldStat");
  TMatrixD *covCS_YieldSyst=(TMatrixD*)fin.Get("covCS_YieldSyst");
  if (!covCS_YieldStat || !covCS_YieldSyst) res=0;
  fin.Close();

  if (res) {
    if (!covCS_YieldStat) res=0;
    if (!covCS_YieldSyst) res=0;
  }

  if (!res) {
    std::cout << "failed to get the final cross section and/or covariance\n";
    return 0;
  }

  TH2D *h2_CSStatErr=hpFinCS.createHistoClone("h2CS_StatErr");
  TH2D *h2_CSSystErr=hpFinCS.createHistoWithSystError("h2CS_SystErr");
  if (res) res=setErrorAsContent(h2_CSStatErr);
  if (res) res=setErrorAsContent(h2_CSSystErr);

  TH2D *h2_covCS_YieldStat=errorFromCov(*covCS_YieldStat,"h2_covCS_YieldStat");
  TH2D *h2_covCS_YieldSyst=errorFromCov(*covCS_YieldSyst,"h2_covCS_YieldSyst");
  delete covCS_YieldStat;
  delete covCS_YieldSyst;
  if (!h2_covCS_YieldStat || !h2_covCS_YieldSyst) res=0;

  if (!res) {
    std::cout << "failed to prepare error histograms\n";
    return 0;
  }

  for (int iSyst=0; iSyst<3; ++iSyst) {

    std::vector<TH2D*> histoV;
    std::vector<TString> labelV;

    if ((iSyst==0) || (iSyst==2)) {
      histoV.push_back(h2_CSStatErr);
      labelV.push_back("propagated yield stat err");
      histoV.push_back(h2_covCS_YieldStat);
      labelV.push_back("stat err from cov");
    }
    if ((iSyst==1) || (iSyst==2)) {
      histoV.push_back(h2_CSSystErr);
      labelV.push_back("propagated yield syst err");
      histoV.push_back(h2_covCS_YieldSyst);
      labelV.push_back("syst err from cov");
    }

    std::vector<std::vector<TH1D*>*> hProfV;
    std::vector<ComparisonPlot_t*> cpV;
    TString canvName=Form("cx%d",iSyst);
    TCanvas *cx=NULL;

    cx= plotProfiles(canvName,histoV,labelV,NULL,0,"cross section error",
		     &hProfV,&cpV,1);

    std::cout << "cpV.size()=" << cpV.size() << "\n";

    if (cpV.size()==1) {
      cpV[0]->ChangeLegendPos(-0.15,0.0,-0.08,0.);
      if (1) {
	//cpV[0]->SetLogy(1);
	//cpV[0]->TransLegend(-0.2,-0.65);
      }
      cpV[0]->Draw(cx);
    }
    else {
      for (unsigned int i=0; i<cpV.size(); ++i) {
	cpV[i]->ChangeLegendPos(-0.2,0.0,-0.08,0.);
	if (1) {
	  cpV[i]->TransLegend(-0.2,-0.65);
	}
	cpV[i]->Draw6(cx,1,i+1);
      }
    }

    if (cx) cx->Update();

    TString fileTag=(iSyst==0) ? "StatErr" : "SystErr";
    if (iSyst==2) fileTag="BothErr";
    TString str2D=Form("_%dD",DYTools::study2D+1);
    TString dirName=Form("plots_csCovStudy%s",str2D.Data());
    TString figName=Form("fig-signalYieldProp") +
      fileTag + str2D + extraTagForSaving;
    std::cout << "figName=<" << figName << ">\n";
    std::cout << "dirName=<" << dirName << ">\n";
    if (saveCanvas) {
      SaveCanvas(cx,figName,dirName);
    }
    else std::cout << " .. not saved as requested\n";
  }

  return 1;
}

// ------------------------------------------------------------------
