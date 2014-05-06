#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include <TRandom.h>
#include "../Include/HistoPair.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/EventSelector.hh"

// --------------------------------------------------------

//int load(const std::vector<TString> &pathV,
//	 const std::vector<TString> &fnameV,
//	 const std::vector<TString> &fieldV,
//	 std::vector<HistoPair2D_t*> &csV);

// --------------------------------------------------------
// --------------------------------------------------------
// Load data from bg-subtracted files


void compareRndSigYields(int analysisIs2D=1,
			 TString conf="defaultAdHoc",
			 int iBr=0,
			 int doSave=0,
			 TString *figName=NULL,
			 TString *dirName=NULL) {

  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Settings
  //==============================================================================================================

  int totalErr=1;
  int loadSyst=1;

  //int the_set=0;
  std::vector<TString> pathV;
  std::vector<TString> fnameV;
  std::vector<TString> fieldV;
  std::vector<TString> labelV;
  TString canvasSaveName, canvasSaveDir;

  std::vector<HistoPair2D_t*> csV;

  double set_ratio_y[2];
  double transLegendX=(DYTools::study2D==1) ? -0.42 : -0.2;
  double transLegendY=0.;

  set_ratio_y[0]=1.00;
  set_ratio_y[1]=1.00;


  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(conf)) return;

  DYTools::TRunMode_t runMode= DYTools::NORMAL_RUN;
  DYTools::TSystematicsStudy_t systModeRef, systMode1, systMode2, systModeV;
  systModeRef=DYTools::APPLY_ESCALE;
  systMode1  =DYTools::SYST_MODE_FAILURE;
  systMode2  =DYTools::SYST_MODE_FAILURE;
  systModeV  =DYTools::SYST_MODE_FAILURE;

  int seedMin=inpMgr.userKeyValueAsInt("SEEDMIN");
  int seedMax=inpMgr.userKeyValueAsInt("SEEDMAX");
  int dSeed=1;
  unsigned int idxRndVec=(unsigned int)(-1);

  //--------------------------------------------------------------------------------------------------------------
  // Define branches
  //==============================================================================================================

  TString extraTag;
  TString plotExtraTag;

  if ((iBr==0) || (iBr==1)) { // added on 2014.04.12

    if (iBr==1) {
      seedMin=-111;
      seedMax= 111;
      dSeed= 222;
      canvasSaveName="fig-CmpEscaleExtremes";
      canvasSaveDir="dir-plots-EScaleStudyExtreme";
    }

    loadSyst=0;
    //prepare(0,pathV,fnameV,fieldV,labelV);
    // Construct eventSelector, update mgr and plot directory
    systModeRef=DYTools::APPLY_ESCALE;
    EventSelector_t evtSelector1(inpMgr,runMode,systModeRef,
		       extraTag, plotExtraTag, EventSelector::_selectDefault);
    if (iBr!=1) {
      pathV .push_back("");
      fnameV.push_back(inpMgr.signalYieldFullFileName(systModeRef,0));
      fieldV.push_back("Input/hpMCsignal");
      labelV.push_back("Zee MC signal");
    }

    pathV .push_back("");
    fnameV.push_back(inpMgr.signalYieldFullFileName(systModeRef,0));
    fieldV.push_back("signalYieldDDbkg");
    labelV.push_back("Data with peak corr.");

    if (1) {
      systMode1=DYTools::ESCALE_DIFF_0000;
      InputFileMgr_t inpMgrUncorr(inpMgr);
      EventSelector_t evtSelector2(inpMgrUncorr,runMode,systMode1,
			extraTag, plotExtraTag, EventSelector::_selectDefault);
      pathV .push_back("");
      fnameV.push_back(inpMgrUncorr.signalYieldFullFileName(systMode1,0));
      fieldV.push_back("signalYieldDDbkg");
      labelV.push_back("Data (regressed)");
    }

    prepare(int((seedMax-seedMin)/dSeed),pathV,fnameV,fieldV,labelV,0,0);
    systModeV=DYTools::ESCALE_STUDY_RND;
    idxRndVec=pathV.size();
    for (int iseed=seedMin; iseed<=seedMax; iseed+=dSeed) {
      //if ((0 || (dSeed==1)) && (iseed-seedMin>79)) break;
      //if (iseed-seedMin>2) break;
      InputFileMgr_t inpMgrRnd(inpMgr);
      inpMgrRnd.editEnergyScaleTag().Append(Form("_RANDOMIZED%d",iseed));
      EventSelector_t evtSelector3(inpMgrRnd,runMode,systModeV,
		       extraTag, plotExtraTag, EventSelector::_selectDefault);
      pathV.push_back("");
      fnameV.push_back(inpMgrRnd.signalYieldFullFileName(systModeV,0));
      fieldV.push_back("signalYieldDDbkg");
      labelV.push_back(Form("Data rnd%d",iseed));
    }

    if (iBr==1) {
      pathV .push_back("");
      fnameV.push_back(inpMgr.signalYieldFullFileName(systModeRef,0));
      fieldV.push_back("Input/hpMCsignal");
      labelV.push_back("Zee MC signal");
    }

    transLegendX=(DYTools::study2D==1) ? -0.42 : -0.1;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code
  //==============================================================================================================

  HERE("main analysis code");

  /*
  {
    canvasSaveName="fig-puStudy-";
    canvasSaveDir="plots-puStudy";
    transLegendX=(DYTools::study2D==1) ? -0.35 : -0.1;
    transLegendY=(DYTools::study2D==1) ? -0.55 : -0.0;
    set_ratio_y[0]=(DYTools::study2D==1) ? 0.9 : 0.96;
    set_ratio_y[1]=(DYTools::study2D==1) ? 1.1 : 1.04;
  }
  */

  if (DYTools::study2D) {
    for (unsigned int i=0; i<fnameV.size(); ++i) {
      fnameV[i].ReplaceAll("preFsr_1D","preFsrDet_2D");
    }
  }

  HERE("loadHistoPairV");

  if (!loadHistoPairV(pathV,fnameV,fieldV,labelV, csV, loadSyst)) {
    std::cout << "failed to load data\n";
    return;
  }

  HERE("loaded");

  std::vector<TH2D*> histoV;
  if (!convertHistoPairVec2HistoVec(csV, histoV, totalErr)) {
    std::cout << "failed to prepare histos\n";
    return;
  }

  std::vector<ComparisonPlot_t*> cpV;
  std::vector<std::vector<TH1D*>*> hProfV;
  int delayDraw=1;
  TCanvas *cx=plotProfiles("cx",histoV,labelV,NULL,0,"signal yield counts",
			   &hProfV, &cpV,delayDraw);
  if (!cx) {
    std::cout << "failed to create canvas with profiles\n";
    return;
  }

  if (0) { // default plottting
    for (unsigned int ic=0; ic<cpV.size(); ++ic) {
      ComparisonPlot_t *cp=cpV[ic];
      cp->TransLegend(transLegendX,transLegendY);
      cp->SetRatioYRange(set_ratio_y[0], set_ratio_y[1]);
      TH1D* hTarget=cp->GetHisto(2); hTarget->SetMarkerStyle(kMultiply);
      if (DYTools::study2D) cp->Draw6(cx,1,ic+1);
      else cp->Draw(cx);
    }
    cx->Update();
  }
  else {
    delete cx;

    if ((iBr==1)) {
      if (DYTools::study2D==0)  cpV[0]->SetLogy(0);

      if (1) { // added on May 05, 2015
	HERE("shift the notation\n");
	for (unsigned int ic=0; ic<cpV.size(); ++ic) {
	  //if (DYTools::study2D==0) cpV[ic]->SetLogy(1);
	  TH1D *hMC=cpV[ic]->GetHisto(4);
	  cpV[ic]->SkipInRatioPlots(hMC);
	  //if (DYTools::study2D==1)
	    cpV[ic]->ErrorsOnRatios(0);

	  TH1D* h1=cpV[ic]->GetHisto(0);
	  h1->SetMarkerStyle(24);
	  h1->SetMarkerColor(kBlue);
	  h1->SetLineColor(kBlue);
	  h1->SetLineStyle(2);
	  TH1D* h2=cpV[ic]->GetHisto(1);
	  h2->SetMarkerStyle(5);
	  h2->SetMarkerColor(30);
	  h2->SetMarkerSize(1.5);
	  h2->SetLineStyle(1);
	  h2->SetLineColor(30);
	  TH1D* h3=cpV[ic]->GetHisto(2);
	  h3->SetMarkerStyle(20);
	  h3->SetMarkerColor(39);
	  h3->SetMarkerSize(1.);
	  h3->SetLineStyle(3);
	  h3->SetLineColor(39);
	  TH1D* h4=cpV[ic]->GetHisto(3);
	  h4->SetMarkerStyle(2);
	  h4->SetMarkerColor(kRed);
	  h4->SetLineStyle(2);
	  h4->SetLineColor(kRed);
	  // signal MC
	  hMC->Scale(0.968);
	  hMC->SetMarkerStyle(20);
	  hMC->SetMarkerColor(kBlack);
	  hMC->SetMarkerSize(0.8);
	  hMC->SetLineStyle(0);
	  hMC->SetLineColor(kBlack);

	}
      }
    }

    for (unsigned int ic=0; ic<cpV.size(); ++ic) {
      TString canvName=Form("cx%d",ic);
      cx=new TCanvas(canvName,canvName,700,800);
      cx->Divide(2,1);
      double div=(DYTools::study2D==0) ? 0.7 : 0.5;
      div=(DYTools::study2D==0) ? 0.29 : 0.5;
      cpV[0]->PreparePads(cx,1,2,"comp",div);

      ComparisonPlot_t *cp=cpV[ic];
      if (DYTools::study2D==1) {
	cp->TransLegend(transLegendX+0.03,-0.5);
	cp->SetXAxisTextSizes(0.1,1.2, 0.1);
	cp->SetYAxisTextSizes(0.09,1., 0.09);
	cp->SetRatioYTitleSize(0.1, 0.9);
	cp->SetRatioLabelSize(0.09);
      }

      cp->SetRatioYRange(set_ratio_y[0], set_ratio_y[1]);
      TH1D* hTarget=cp->GetHisto(2); hTarget->SetMarkerStyle(kMultiply);
      cp->Draw(cx);
      cx->Update();
    }
  }


  // plot count dirstributions of a mass bin
  if (0 && (DYTools::study2D==0)) {
    std::cout << "hProfV.size()=" << hProfV.size() << "\n";
    std::vector<TH1D*> hMassV;
    std::vector<TString> massStrV;
    massStrV.reserve(DYTools::nMassBins);
    for (int im=0; im<DYTools::nMassBins; ++im) {
      massStrV.push_back(Form("M_%1.0lf_%1.0lf",
		  DYTools::massBinLimits[im],DYTools::massBinLimits[im+1]));
    }
    if (!createAnyH1Vec(hMassV,"hMassProf_",massStrV,200,-1e6,1e6,
			"yield count","count",1)) return;
    std::vector<TH1D*> *histos=hProfV[0];
    for (int im=14; im<20; im++) {
      double avg=0.;
      int count=0;
      for (unsigned int ii=0; ii<histos->size(); ++ii) {
	avg+=(*histos)[ii]->GetBinContent(im);
	count++;
      }
      avg/=double(count);
      for (unsigned int ii=0; ii<histos->size(); ++ii) {
	hMassV[im]->Fill((*histos)[ii]->GetBinContent(im) - avg);
      }

      TString canvName=TString("canvProfAt_") + massStrV[im];
      TCanvas *cm=new TCanvas(canvName,canvName,600,600);
      hMassV[im]->Draw();
      cm->Update();
    }
  }

  TString outName=canvasSaveName + DYTools::analysisTag;
  if (doSave) {
    SaveCanvas(cx,outName,canvasSaveDir);
  }
  else {
    std::cout << "would save to <" << outName << "> in <" << canvasSaveDir << ">\n";
  }
  if (figName) *figName= outName;
  if (dirName) *dirName= canvasSaveDir;

  // Covariance study
  if (1) {
    std::vector<TH2D*> hRndVec;
    hRndVec.reserve(histoV.size());
    for (unsigned int i=idxRndVec; i<histoV.size(); ++i) {
      hRndVec.push_back(histoV[i]);
    }

    int unbiasedEstimate=1;
    TH2D* avgDistr=createBaseH2("hSignalYieldAvgDistr");
    TMatrixD* cov= deriveCovMFromRndStudies(hRndVec,unbiasedEstimate,avgDistr);
    TMatrixD* corr=corrFromCov(*cov);

    TH2D* h2Cov=createHisto2D(*cov,NULL,"h2Cov","",_colrange_default);
    TH2D *h2Corr=createHisto2D(*corr,NULL,"h2Corr","",_colrange_center,0,1.0001);
    TH2D* h2CovSave=clipToAnalysisUnfBins(h2Cov,"h2CovSave","",1);
    TH2D* h2CorrSave=clipToAnalysisUnfBins(h2Corr,"h2CorrSave","",1);
    h2CorrSave->GetZaxis()->SetRangeUser(-1.0001,1.0001);

    TCanvas *canvCov= new TCanvas("ccov","ccov",900,900);
    AdjustFor2DplotWithHeight(canvCov);
    //cov->Draw("COLZ");
    h2CovSave->Draw("COLZ");
    canvCov->Update();
    TCanvas *canvCorr= new TCanvas("ccorr","ccorr",900,900);
    AdjustFor2DplotWithHeight(canvCorr);
    //corr->Draw("COLZ");
    h2CorrSave->Draw("COLZ");
    canvCorr->Update();
  }

  return;
}

// --------------------------------------------------------
// --------------------------------------------------------
// --------------------------------------------------------
// --------------------------------------------------------
