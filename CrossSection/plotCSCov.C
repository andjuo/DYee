#include "../Include/DYTools.hh"
#include "../CrossSection/crossSectionFnc.hh"
#include "../Include/colorPalettes.hh"
#include "CSCovWorkFlags.hh"
#include <TBenchmark.h>

//================================================================================================================

int workWithData(TCovData_t &dt, const WorkFlags_t &wf);


//=== MAIN MACRO =================================================================================================

// the_case: 1,5 -- plotTotCov, 0,2,3,4 -- plotAllCovs
// actually useful: 2,3,4,5


int plotCSCov(int analysisIs2D, TString conf, int the_case, int workBranch,
	      int showCSCov=1,
	      TString outFileExtraTag_UserInput="",
	      int saveTotCovDetails_user=0)
{

  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return retCodeError;
  }

  TString extraESFsystematics;
  extraESFsystematics="20140525";

  // Settings 
  //==============================================================================================================

  DYTools::TCrossSectionKind_t csKind=(DYTools::study2D) ? DYTools::_cs_preFsrDet : DYTools::_cs_preFsr;

  DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN;
  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;

  TCovData_t dt;
  WorkFlags_t work(the_case,showCSCov,outFileExtraTag_UserInput);
  work.saveTotCovDetails(saveTotCovDetails_user);

  CSCovCalcFlags_t *cf= & work.editCalcFlags();

  switch(abs(workBranch)) {
  case 0:
    cf->calc_YieldStatDetailed(1);
    cf->calc_YieldSystDetailed(1);
    break;
  case 1:
    cf->calc_YieldStatDetailed(1);
    cf->calc_YieldSystDetailed(1);
    cf->calc_YieldEscale(1);
    break;
  case 2:
    cf->calc_UnfRnd(1);
    break;
  case 3:
    cf->calc_UnfRnd(1);
    cf->calc_UnfEScale(1);
    break;
  case 4:
    cf->calc_UnfEScale(1);
    break;
  case 5:
    cf->calc_ESFtot(1);
    break;
  case 6:
    cf->calc_YieldStatDetailed(1);
    cf->calc_YieldSystDetailed(1);
    cf->calc_YieldEscale(1);  // added on May 27, 2014
    cf->calc_UnfRnd(1);
    //cf->calc_UnfEScale(1); // commented out on May 27, 2014
    cf->calc_ESFtot(1);
    cf->calc_EffRnd(1);
    cf->calc_AccRnd(1); // not active in 2D!
    cf->calc_FsrRnd(1);
    cf->calc_globalFSR(1);
    cf->calc_globalPU(1);
    break;
  case 7:
    cf->calc_globalPU(1);
    break;
  case 8:
    cf->calc_globalFSR(1);
    break;
  case 9:
    cf->calc_YieldStatDetailed(1);
    break;
  default:
    std::cout << "workBranch=" << workBranch << " is not ready\n";
    return retCodeError;
  }

  if (!work.finalizeFlags()) {
    std::cout << "error from finalizeFlags\n";
    return retCodeError;
  }

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

  if ((the_case==2) || (the_case==3) || (the_case==4)
      || (the_case==5)) {
    // global study I
    work.init_ExtraTagV(0);
    work.extraFileTag(_yield, "-yieldOnly_nExps1000");
    //work.extraFileTag(_yield, "-yieldOnly");
    work.extraFileTag(_corrUnf, "-unfOnly");
    if ((workBranch==2) ||
	(workBranch==5)) work.extraFileTag(_corrUnf, "-unfRndOnly_nExps1000");
    //if ((workBranch==6)) work.extraFileTag(_corrUnf,"-unfOnly_nExps1000");
    if (workBranch==4) work.extraFileTag(_corrUnf, "-unfOnly_nExps20");
    work.extraFileTag(_corrEff, "-effRndOnly_nExps1000");
    work.extraFileTag(_corrESF, "-esfOnly");
    work.extraFileTag(_corrAcc, "-accRndOnly");
    work.extraFileTag(_corrFSR, "-fsrRndOnly");
    //if (the_case==3) 
    //work.extraFileTag(_corrFSR, "-fsrRndOnly_nExps1000");
    //work.extraFileTag(_corrGlobalFSR, "-globalFSROnly_nExps20");
    //work.extraFileTag(_corrGlobalPU, "-globalPUOnly_nExps20");
    work.extraFileTag(_corrGlobalFSR, "-globalFSROnly");
    work.extraFileTag(_corrGlobalPU, "-globalPUOnly");
    /*
    if (the_case==4) {
      std::cout << "\n\n\t EXCLUDING RNDs beyond rho\n";
      cf->calc_EffRnd(0);
      cf->calc_AccRnd(0);
      cf->calc_FsrRnd(0);
    }
    */
  }

  int res=1;
  if (res && cf->doCalcYieldCov()) {
    if (!loadYieldCovMatrices(fnameBase,dt.covYieldV,dt.labelYieldV,work)) return 0;
    dt.isActive[0]=1;
  }
  if (res && cf->doCalcUnfCov()) {
    if (!loadUnfCovMatrices(fnameBase,dt.covUnfV,dt.labelUnfV,work)) return 0;
    dt.isActive[1]=1;
  }
  if (res && cf->doCalcEffCov()) {
    if (!loadEffCovMatrices(fnameBase,dt.covEffV,dt.labelEffV,work)) return 0;
    dt.isActive[2]=1;
  }
  if (res && cf->doCalcESFCov()) {
    if (!loadEsfCovMatrices(fnameBase,dt.covEsfV,dt.labelEsfV,work)) return 0;
    if (extraESFsystematics.Length()) {
      if (!dt.addESFsyst(extraESFsystematics)) {
	std::cout << "error adding extra ESF systematics\n";
	return 0;
      }
    }
    dt.isActive[3]=1;
  }
  if (res && cf->doCalcAccCov()) {
    if (!loadAccCovMatrices(fnameBase,dt.covAccV,dt.labelAccV,work)) return 0;
    dt.isActive[4]=1;
  }
  if (res && cf->doCalcFSRCov()) {
    if (!loadFsrCovMatrices(fnameBase,dt.covFsrV,dt.labelFsrV,work)) return 0;
    dt.isActive[5]=1;
  }
  if (res && cf->doCalcGlobalCov()) {
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
      //if (iCorr) continue;
      TString covStr;
      TH2D* h2Main=NULL;
      TString yAxisLabel="uncertainty from cov";
      TString figTag;
      switch(iCorr) {
      case 0:  covStr="CSCov_"; figTag.Clear(); break;
      case 1:
	covStr="CSCov_";
	h2Main=loadMainCSResult(1);  // load CS
	if (!h2Main) return;
	removeError(h2Main);
	h2Main->Scale(0.01);
	yAxisLabel="relative uncertainty from cov (%)";
	figTag="frac";
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
	  TString histoLabel=
	    TString(Form("histoErr_%d_",iCorr)) + (*labelV)[i];
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
      int delayDraw=1;

      TCanvas *cx=plotProfiles(canvName,
			       errFromCovV, errFromCovLabelV,
			       NULL,1, yAxisLabel,
			       &hProfV, &cpV,
			       delayDraw);
      for (unsigned int i=0; i<cpV.size(); ++i) {
	if (DYTools::study2D) {
	  double dy=(iCorr==0) ? -0.65 : -0.2;
	  cpV[i]->TransLegend(-0.4,dy);
	  cpV[i]->Draw6(cx,1,i+1);
	}
	else {
	  double dx=(iCorr==0) ? -0.1 : -0.4;
	  cpV[i]->TransLegend(dx,0.);
	  cpV[i]->Draw(cx);
	}
      }
      if (delayDraw) cx->Update();

      if (1) {
	TString figName=TString("fig-") + DYTools::analysisTag +
	  TString("--") + "errorProfiles";
	if (figTag.Length()) figName.Append(figTag);
	figName.Append(wf.extraFileTag());
	//eliminateSeparationSigns(figName);
	std::cout << "figName=<" << figName << ">\n";
	SaveCanvas(cx,figName);
      }

      // Save table
      if (0 && (iCorr==1)) {
	TH2D *totErrorFromCov= errorFromCov( *totalCov, "totalErr" );
	if (!scaleHisto(totErrorFromCov,h2Main)) return;
	errFromCovV.push_back(totErrorFromCov);
	errFromCovLabelV.push_back("total error");
	TString tableTag=figTag + wf.extraFileTag();
	if (!saveLatexTable(tableTag,errFromCovV,errFromCovLabelV,
			    "%5.2lf",0,1)) {
	  std::cout << "failed to save table\n";
	  return;
	}
      }
    }
  }

  if (totalCov) delete totalCov;
  return;
}


// -----------------------------------------------------------

void plotTotCov(TCovData_t &dt, const WorkFlags_t &wf) {
  if (wf.showCSCov()==0) {
    std::cout << "plotTotCov needs CS covariance\n";
    return;
  }

  TMatrixD *totalCov=dt.calcTotalCov();

  // save total covariance
  if (1) {
    TString fname=TString(Form("finalCov-%dD",DYTools::study2D+1));
    fname.Append(wf.extraFileTag() + TString(".root"));
    TFile fout(fname,"recreate");
    totalCov->Write("totalCov");
    TMatrixD *corr= corrFromCov(*totalCov);
    corr->Write("totalCorr");
    delete corr;
    if (wf.saveTotCovDetails()) {
      if (!dt.Write("details")) {
	std::cout << "failed to save the individual covariances\n";
	return;
      }
    }
    else { std::cout << " details not saved, as requested\n"; }
    writeBinningArrays(fout,"plotCSCov");
    fout.Close();
    std::cout << "\n\tfile <" << fout.GetName() << "> created\n\n";
  }

  for (int iCorr=0; iCorr<4; ++iCorr) {
    if (iCorr==2) continue; // not ready
    if (iCorr==3) continue; // not ready
    TString covStr;
    switch(iCorr) {
    case 0: covStr="Cov"; break;
    case 1: covStr="Corr"; break;
    case 2: covStr="partCorr"; break;
    case 3: covStr="RelCov"; break;
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
      TH2D *h2=loadMainCSResult(1);
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
      TString figName=TString("fig-") + DYTools::analysisTag +
	TString("--total-") + covStr;
      figName.Append(wf.extraFileTag());
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
  else if ((wf.theCase()==1) || (wf.theCase()==5)) {
    plotTotCov(dt,wf);
  }

  return 1;
}




// -----------------------------------------------------------
// -----------------------------------------------------------
