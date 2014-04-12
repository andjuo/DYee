#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include <TRandom.h>
#include "../Include/HistoPair.hh"

// --------------------------------------------------------

void prepare(int count,
	     std::vector<TString> &pathV,
	     std::vector<TString> &fnameV,
	     std::vector<TString> &fieldV,
	     std::vector<TString> &labelV);

//int load(const std::vector<TString> &pathV,
//	 const std::vector<TString> &fnameV,
//	 const std::vector<TString> &fieldV,
//	 std::vector<HistoPair2D_t*> &csV);

// --------------------------------------------------------
// --------------------------------------------------------

void compareCS(int iBr=0,
	       int analysisIs2D=1,
	       int doSave=0,
	       TString *figName=NULL,
	       TString *dirName=NULL) {

  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return;
  }

  int the_set=0;
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

  TString gen_path="../../Results-DYee/root_files_reg/xsec/DY_j22_19712pb";
  TString csDDBkg="csYieldDDbkg";
  TString csMCBkg="csYieldMCbkg";

  if (0) { // added on 2014.04.11
    prepare(2,pathV,fnameV,fieldV,labelV);
    pathV[0]=gen_path + TString("-ESFUnregEn/");
    pathV[1]=gen_path + TString("-ESFregEn/");
    fnameV[0]="xSec_preFsr_1D.root";
    fnameV[1]="xSec_preFsr_1D.root";
    fieldV[0]=csDDBkg;
    fieldV[1]=csDDBkg;
    labelV[0]="ESF unreg.en";
    labelV[1]="ESF reg.en";
    canvasSaveName="fig-regVsUnreg-";
    canvasSaveDir="plots-regVsUnreg";
    transLegendX=(DYTools::study2D==1) ? -0.35 : -0.2;
    transLegendY=(DYTools::study2D==1) ? -0.55 : -0.2;
    set_ratio_y[0]=0.98;
    set_ratio_y[1]=1.02;
  }

  if (0) { // added on 2014.04.11
    prepare(3,pathV,fnameV,fieldV,labelV);
    pathV[0]=gen_path + TString("-ESFregEn/");
    pathV[1]=gen_path + TString("_FSR_5plus-ESFFSR5p/");
    pathV[2]=gen_path + TString("_FSR_5minus-ESFFSR5m/");
    fnameV[0]="xSec_preFsr_1D.root";
    fnameV[1]="xSec_preFsr_1D.root";
    fnameV[2]="xSec_preFsr_1D.root";
    fieldV[0]=csDDBkg;
    fieldV[1]=csDDBkg;
    fieldV[2]=csDDBkg;
    labelV[0]="ESF reg.en";
    labelV[1]="FSR 5plus";
    labelV[2]="FSR 5minus";
    canvasSaveName="fig-fsrStudy-";
    canvasSaveDir="plots-fsrStudy";
    transLegendX=(DYTools::study2D==1) ? -0.35 : -0.1;
    transLegendY=(DYTools::study2D==1) ? -0.55 : -0.0;
    set_ratio_y[0]=(DYTools::study2D==1) ? 0.98 : 0.97;
    set_ratio_y[1]=(DYTools::study2D==1) ? 1.02 : 1.03;
  }

  if (1) { // added on 2014.04.11
    prepare(3,pathV,fnameV,fieldV,labelV);
    pathV[0]=gen_path + TString("-ESFregEn/");
    pathV[1]=gen_path + TString("_Pileup5plus-ESFPU5p/");
    pathV[2]=gen_path + TString("_Pileup5minus-ESFPU5m/");
    fnameV[0]="xSec_preFsr_1D.root";
    fnameV[1]="xSec_preFsr_1D.root";
    fnameV[2]="xSec_preFsr_1D.root";
    fieldV[0]=csDDBkg;
    fieldV[1]=csDDBkg;
    fieldV[2]=csDDBkg;
    labelV[0]="ESF reg.en";
    labelV[1]="PU 5plus";
    labelV[2]="PU 5minus";
    canvasSaveName="fig-puStudy-";
    canvasSaveDir="plots-puStudy";
    transLegendX=(DYTools::study2D==1) ? -0.35 : -0.1;
    transLegendY=(DYTools::study2D==1) ? -0.55 : -0.0;
    set_ratio_y[0]=(DYTools::study2D==1) ? 0.9 : 0.96;
    set_ratio_y[1]=(DYTools::study2D==1) ? 1.1 : 1.04;
  }

  if (DYTools::study2D) {
    for (unsigned int i=0; i<fnameV.size(); ++i) {
      fnameV[i].ReplaceAll("preFsr_1D","preFsrDet_2D");
    }
  }


  if (!loadHistoPairV(pathV,fnameV,fieldV,labelV, csV)) {
    std::cout << "failed to load data\n";
    return;
  }

  int totalErr=1;
  std::vector<TH2D*> histoV;
  if (!convertHistoPairVec2HistoVec(csV, histoV, totalErr)) {
    std::cout << "failed to prepare histos\n";
    return;
  }

  std::vector<ComparisonPlot_t*> cpV;
  int delayDraw=1;
  TCanvas *cx=plotProfiles("cx",histoV,labelV,NULL,0,"cross section",
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

void prepare(int count,
		    std::vector<TString> &pathV,
		    std::vector<TString> &fnameV,
		    std::vector<TString> &fieldV,
		    std::vector<TString> &labelV) {
  pathV.clear();
  fnameV.clear();
  fieldV.clear();
  labelV.clear();
  pathV.reserve(count);
  fnameV.reserve(count);
  fieldV.reserve(count);
  labelV.reserve(count);
  for (int i=0; i<count; ++i) {
    TString empty=Form("empty_%d",i);
    pathV.push_back(empty);
    fnameV.push_back(empty);
    fieldV.push_back(empty);
    labelV.push_back(empty);
  }
}

// --------------------------------------------------------
// --------------------------------------------------------

// --------------------------------------------------------
// --------------------------------------------------------
