#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include <TRandom.h>
#include "../Include/HistoPair.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/EventSelector.hh"
#include "../Include/UnfoldingMatrix.h"

// ----------------------------------------

void swapStr(TString &a, TString &b) {
  TString tmp=a; a=b; b=tmp;
}

// --------------------------------------------------------
// --------------------------------------------------------

void compareRndYields(TString conf="defaultAdHoc",
		      int iBr=0,
		      int analysisIs2D=1,
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
  double unfYieldsExtraFactor=1/0.148613;

  int the_set=0;
  std::vector<TString> pathV, m_pathV;
  std::vector<TString> fnameV, m_fnameV;
  std::vector<TString> fieldV, m_fieldV;
  std::vector<TString> labelV, m_labelV;
  TString canvasSaveName, canvasSaveDir;

  std::vector<HistoPair2D_t*> csV;

  double set_ratio_y[2];
  double transLegendX=(DYTools::study2D==1) ? -0.42 : -0.2;
  double transLegendY=0.;

  set_ratio_y[0]=1.00;
  set_ratio_y[1]=1.00;


  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(conf)) return;
  InputFileMgr_t inpMgrRemote;
  if (!inpMgrRemote.Load(conf + TString("Remote"))) return;

  DYTools::TRunMode_t runMode= DYTools::NORMAL_RUN;
  DYTools::TSystematicsStudy_t systModeRef, systMode1, systMode2, systMode3;
  DYTools::TSystematicsStudy_t systModeV;
  systModeRef=DYTools::APPLY_ESCALE;
  systMode1  =DYTools::SYST_MODE_FAILURE;
  systMode2  =DYTools::SYST_MODE_FAILURE;
  systModeV  =DYTools::SYST_MODE_FAILURE;

  int seedMin=inpMgr.userKeyValueAsInt("SEEDMIN");
  int seedMax=inpMgr.userKeyValueAsInt("SEEDMAX");
  int dSeed=1;

  //--------------------------------------------------------------------------------------------------------------
  // Define branches
  //==============================================================================================================

  TString extraTag;
  TString plotExtraTag;

  // Check that MC reco yields are similar to data yields
  if (iBr==0) { // added on 2014.04.12
    loadSyst=0;
    if (1) {
    prepare(2,pathV,fnameV,fieldV,labelV);
    // Construct eventSelector, update mgr and plot directory
    systModeRef=DYTools::APPLY_ESCALE;
    EventSelector_t evtSelector1(inpMgr,runMode,systModeRef,
		       extraTag, plotExtraTag, EventSelector::_selectDefault);
    pathV [0]="";
    fnameV[0]=inpMgr.yieldFullFileName(-1,systModeRef,0);
    fieldV[0]="yields/hYield_data";
    labelV[0]="Raw data with peak corr.";

    systMode2=DYTools::ESCALE_DIFF_0000;
    EventSelector_t evtSelector2(inpMgr,runMode,systMode2,
		       extraTag, plotExtraTag, EventSelector::_selectDefault);
    pathV [1]="";
    fnameV[1]=inpMgr.yieldFullFileName(-1,systMode2,0);
    fieldV[1]="yields/hYield_data";
    labelV[1]="Raw data (regressed)";
    }

    // vectors for matrices
    prepare(1,m_pathV,m_fnameV,m_fieldV,m_labelV);
   
    systMode3=DYTools::NO_SYST;
    EventSelector_t evtSelector3(inpMgr,runMode,systMode3,
		       extraTag, plotExtraTag, EventSelector::_selectDefault);
    UnfoldingMatrix_t detResponse(UnfoldingMatrix::_cDET_Response,"detResponse");
    TString matrixFName,yieldsFName;
    TString yieldGenName, yieldRecName;
    detResponse.getFileNames(inpMgrRemote.constDir(systMode3,0),
			     UnfoldingMatrix_t::generateFNameTag(systMode3),
			     matrixFName,yieldsFName);
    detResponse.getYieldNames(UnfoldingMatrix::_cDET_Response,
			      yieldGenName,yieldRecName);
    //swapStr(yieldGenName,yieldRecName);
    m_pathV [0]="";
    m_fnameV[0]=yieldsFName;
    m_fieldV[0]=yieldRecName;
    m_labelV[0]="MC reco vec";
    std::cout << "yieldsFName=" << yieldsFName << ", yieldRecName=" << yieldRecName << "\n";
    //return;

    seedMin=-111;
    seedMax= 111;
    dSeed= 222;
    // add space
    prepare(int((seedMax-seedMin)/dSeed),m_pathV,m_fnameV,m_fieldV,m_labelV,0,0);
    systModeV=DYTools::RESOLUTION_STUDY;
    for (int iseed=seedMin; iseed<=seedMax; iseed+=dSeed) {
      //if (iseed<0) continue;
      //if (iseed-seedMin>2) break;
      InputFileMgr_t inpMgrRnd(inpMgrRemote);
      inpMgrRnd.editEnergyScaleTag().Append(Form("_INVERTED_RANDOMIZED%d",iseed));
      EventSelector_t evtSelectorRnd(inpMgrRnd,runMode,systModeV,
		       extraTag, plotExtraTag, EventSelector::_selectDefault);
      TString nameRnd=Form("detResponse_seed%d",iseed);
      UnfoldingMatrix_t detResponseRnd(UnfoldingMatrix::_cDET_Response,nameRnd);
      TString matrixFNameRnd,yieldsFNameRnd;
      TString yieldGenNameRnd, yieldRecNameRnd;
      detResponseRnd.getFileNames(inpMgrRnd.constDir(systModeV,0),
			    UnfoldingMatrix_t::generateFNameTag(systModeV),
				  matrixFNameRnd,yieldsFNameRnd);
      detResponseRnd.getYieldNames(UnfoldingMatrix::_cDET_Response,
				   yieldGenNameRnd,yieldRecNameRnd);
      //swapStr(yieldGenNameRnd,yieldRecNameRnd);
      m_pathV.push_back("");
      m_fnameV.push_back(yieldsFNameRnd);
      m_fieldV.push_back(yieldRecNameRnd);
      m_labelV.push_back(Form("MC reco rnd%d",iseed));
      std::cout << "yieldsFName=" << yieldsFName << ", yieldRecName=" << yieldRecName << "\n";
    }

    transLegendX=(DYTools::study2D==1) ? -0.42 : -0.1;
  } // (iBr==0)

  // ----------------------------------------------
  // Check that MC reco yields are similar to data _signal_ yields
  if (iBr==1) { // added on 2014.04.14
    loadSyst=0;
    if (0) {
    prepare(2,pathV,fnameV,fieldV,labelV);
    // Construct eventSelector, update mgr and plot directory
    systModeRef=DYTools::APPLY_ESCALE;
    EventSelector_t evtSelector1(inpMgr,runMode,systModeRef,
		       extraTag, plotExtraTag, EventSelector::_selectDefault);
    pathV [0]="";
    fnameV[0]=inpMgr.signalYieldFullFileName(systModeRef,0);
    fieldV[0]="signalYieldDDbkg";
    labelV[0]="Data signal with peak corr.";

    systMode2=DYTools::ESCALE_DIFF_0000;
    EventSelector_t evtSelector2(inpMgr,runMode,systMode2,
		       extraTag, plotExtraTag, EventSelector::_selectDefault);
    pathV [1]="";
    fnameV[1]=inpMgr.signalYieldFullFileName(systMode2,0);
    fieldV[1]="signalYieldDDbkg";
    labelV[1]="Data signal (regressed)";
    }

    // vectors for matrices
    prepare(1,m_pathV,m_fnameV,m_fieldV,m_labelV);
   
    systMode3=DYTools::NO_SYST;
    EventSelector_t evtSelector3(inpMgr,runMode,systMode3,
		       extraTag, plotExtraTag, EventSelector::_selectDefault);
    UnfoldingMatrix_t detResponse(UnfoldingMatrix::_cDET_Response,"detResponse");
    TString matrixFName,yieldsFName;
    TString yieldGenName, yieldRecName;
    detResponse.getFileNames(inpMgrRemote.constDir(systMode3,0),
			     UnfoldingMatrix_t::generateFNameTag(systMode3),
			     matrixFName,yieldsFName);
    detResponse.getYieldNames(UnfoldingMatrix::_cDET_Response,
			      yieldGenName,yieldRecName);
    //swapStr(yieldGenName,yieldRecName);
    m_pathV [0]="";
    m_fnameV[0]=yieldsFName;
    m_fieldV[0]=yieldRecName;
    m_labelV[0]="MC reco vec";
    std::cout << "yieldsFName=" << yieldsFName << ", yieldRecName=" << yieldRecName << "\n";
    //return;

    seedMin=-111;
    seedMax= 111;
    dSeed= 222;
    // add space
    prepare(int((seedMax-seedMin)/dSeed),m_pathV,m_fnameV,m_fieldV,m_labelV,0,0);
    systModeV=DYTools::RESOLUTION_STUDY;
    for (int iseed=seedMin; iseed<=seedMax; iseed+=dSeed) {
      //if (iseed<0) continue;
      //if (iseed-seedMin>2) break;
      InputFileMgr_t inpMgrRnd(inpMgrRemote);
      inpMgrRnd.editEnergyScaleTag().Append(Form("_INVERTED_RANDOMIZED%d",iseed));
      EventSelector_t evtSelectorRnd(inpMgrRnd,runMode,systModeV,
		       extraTag, plotExtraTag, EventSelector::_selectDefault);
      TString nameRnd=Form("detResponse_seed%d",iseed);
      UnfoldingMatrix_t detResponseRnd(UnfoldingMatrix::_cDET_Response,nameRnd);
      TString matrixFNameRnd,yieldsFNameRnd;
      TString yieldGenNameRnd, yieldRecNameRnd;
      detResponseRnd.getFileNames(inpMgrRnd.constDir(systModeV,0),
			    UnfoldingMatrix_t::generateFNameTag(systModeV),
				  matrixFNameRnd,yieldsFNameRnd);
      detResponseRnd.getYieldNames(UnfoldingMatrix::_cDET_Response,
				   yieldGenNameRnd,yieldRecNameRnd);
      //swapStr(yieldGenNameRnd,yieldRecNameRnd);
      m_pathV.push_back("");
      m_fnameV.push_back(yieldsFNameRnd);
      m_fieldV.push_back(yieldRecNameRnd);
      m_labelV.push_back(Form("MC reco rnd%d",iseed));
      std::cout << "yieldsFName=" << yieldsFName << ", yieldRecName=" << yieldRecName << "\n";
    }

    transLegendX=(DYTools::study2D==1) ? -0.42 : -0.1;
  } // (iBr==1)


  if (iBr==100) { // added on 2014.04.12
    loadSyst=0;
    prepare(2,pathV,fnameV,fieldV,labelV);
    // Construct eventSelector, update mgr and plot directory
    systModeRef=DYTools::APPLY_ESCALE;
    EventSelector_t evtSelector1(inpMgr,runMode,systModeRef,
		       extraTag, plotExtraTag, EventSelector::_selectDefault);
    pathV [0]="";
    fnameV[0]=inpMgr.yieldFullFileName(-1,systModeRef,0);
    fieldV[0]="yields/hYield_data";
    labelV[0]="Data with peak corr.";

    systMode1=DYTools::ESCALE_DIFF_0000;
    EventSelector_t evtSelector2(inpMgr,runMode,systMode1,
		       extraTag, plotExtraTag, EventSelector::_selectDefault);
    pathV [1]="";
    fnameV[1]=inpMgr.yieldFullFileName(-1,systMode1,0);
    fieldV[1]="yields/hYield_data";
    labelV[1]="Data (regressed)";

    transLegendX=(DYTools::study2D==1) ? -0.42 : -0.1;
    systModeV=DYTools::ESCALE_STUDY_RND;
    for (int iseed=seedMin; iseed<=seedMax; ++iseed) {
      //if (iseed-seedMin>2) break;
      InputFileMgr_t inpMgrRnd(inpMgr);
      inpMgrRnd.editEnergyScaleTag().Append(Form("_RANDOMIZED%d",iseed));
      EventSelector_t evtSelector3(inpMgrRnd,runMode,systModeV,
		       extraTag, plotExtraTag, EventSelector::_selectDefault);
      pathV.push_back("");
      fnameV.push_back(inpMgrRnd.yieldFullFileName(-1,systModeV,0));
      fieldV.push_back("yields/hYield_data");
      labelV.push_back(Form("Data rnd%d",iseed));
    }
    
  }

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //============================================================================================================== 


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

  if (!loadHistoPairV(pathV,fnameV,fieldV,labelV, csV, loadSyst)) {
    std::cout << "failed to load data\n";
    return;
  }

  

  std::vector<TH2D*> histoV;
  if (!convertHistoPairVec2HistoVec(csV, histoV, totalErr)) {
    std::cout << "failed to prepare histos\n";
    return;
  }

  for (unsigned int im=0; im < m_fnameV.size(); ++im) {
    TString fullName= m_pathV[im] + m_fnameV[im];
    HERE("load %s from %s",m_fieldV[im],fullName);
    TH2D* h2= LoadMatrixFields(fullName,1,m_fieldV[im],"",0,1);
    if (!h2) {
      std::cout << "error\n";
      return;
    }
    h2->Scale(unfYieldsExtraFactor);
    histoV.push_back(h2);
    labelV.push_back(m_labelV[im]);
  }

  std::vector<ComparisonPlot_t*> cpV;
  int delayDraw=1;
  TCanvas *cx=plotProfiles("cx",histoV,labelV,NULL,0,"yield counts",
			   NULL, &cpV,delayDraw);
  if (!cx) {
    std::cout << "failed to create canvas with profiles\n";
    return;
  }

  for (unsigned int ic=0; ic<cpV.size(); ++ic) {
    ComparisonPlot_t *cp=cpV[ic];
    cp->TransLegend(transLegendX,transLegendY);
    cp->SetRatioYRange(set_ratio_y[0], set_ratio_y[1]);
    if (DYTools::study2D) cp->Draw6(cx,1,ic+1);
    else cp->Draw(cx);
  }
  cx->Update();


  TString outName=canvasSaveName + DYTools::analysisTag;
  if (doSave) {
    SaveCanvas(cx,outName,canvasSaveDir);
  }
  else {
    std::cout << "would save to <" << outName << "> in <" << canvasSaveDir << ">\n";
  }
  if (figName) *figName= outName;
  if (dirName) *dirName= canvasSaveDir;

  return;
}

// --------------------------------------------------------
// --------------------------------------------------------
// --------------------------------------------------------
// --------------------------------------------------------
