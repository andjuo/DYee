#include <TBenchmark.h>
#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/CPlot.hh"
#include "../Include/MitStyleRemix.hh"
#include "../Include/ComparisonPlot.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/EventSelector.hh"
#include "../Include/HistoPair.hh"

//Int_t minMutualMultiple();
//Int_t minMutualMultipleTwo(Int_t n1, Int_t n2);
//Bool_t checkMatrixSize(TMatrixD m);

//void unsquareMatrixElements(TMatrixD &m);

//template<class T>
//inline T SQR(const T &x) { return (x)*(x); }

/*
void bkgTablesToLatex(TMatrixD true2eBackground, TMatrixD true2eBackgroundError, TMatrixD true2eBackgroundErrorSyst, 
                      TMatrixD fakeEleBackground, TMatrixD  fakeEleBackgroundError, TMatrixD fakeEleBackgroundErrorSyst, 
                      TMatrixD true2eBackgroundFromData, TMatrixD true2eBackgroundFromDataError, TMatrixD true2eBackgroundFromDataErrorSyst, 
                      TMatrixD fakeEleBackgroundFromData, TMatrixD  fakeEleBackgroundFromDataError, TMatrixD fakeEleBackgroundFromDataErrorSyst, 
                      TMatrixD totalBackground, TMatrixD totalBackgroundError, TMatrixD totalBackgroundErrorSyst);
*/

// -----------------------------------------------------------------------------

int subtractBackground(int analysisIs2D,
		       const TString conf = "default",
		       DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN,
		       DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST,
		       int iSeed=-1) {


  gBenchmark->Start("subtractBackground");

  {
    using namespace DYTools;
    DYTools::printExecMode(runMode,systMode);
    const int debug_print=1;
    if (!DYTools::checkSystMode(systMode,debug_print,9,
				DYTools::NO_SYST, //DYTools::ESCALE_STUDY, 
				DYTools::ESCALE_STUDY_RND,
	        DYTools::UNREGRESSED_ENERGY,DYTools::APPLY_ESCALE,
		ESCALE_DIFF_0000, ESCALE_DIFF_0005, ESCALE_DIFF_0010, ESCALE_DIFF_0015, ESCALE_DIFF_0020
				)) 
      return retCodeError;
  }

  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return retCodeError;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================

  // Read from configuration file only the location of the root files
  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(conf)) 
    return retCodeError;
  if (systMode==DYTools::ESCALE_STUDY_RND) {
    inpMgr.editEnergyScaleTag().Append(Form("_RANDOMIZED%d",iSeed));
  }
  //inpMgr.Print();

  // Construct eventSelector, update inpMgr and plot directory
  TString extraTag;
  TString plotExtraTag;
  EventSelector_t evtSelector(inpMgr,runMode,systMode,
		       extraTag, plotExtraTag, EventSelector::_selectDefault);

  //std::cout <<" " << inpMgr.yieldFullName(-1,systMode,0) << "\n";
  //std::cout <<" " << inpMgr.signalYieldFullName(systMode) << "\n";

  // Save background-subtracted signal yields
  // inputDir+TString("/yields_bg-subtracted.root")
  const int ignoreDebugRunFlag=0;
  TString outFileName= inpMgr.signalYieldFullFileName(systMode,ignoreDebugRunFlag);
  std::cout << "generated outFileName=<" << outFileName << ">\n";

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //============================================================================================================== 


  vector<TH2D*> yields;
  createBaseH2Vec(yields,"hYield_",inpMgr.sampleNames());

  {
    int createDir=0;
    TString yieldFullName= inpMgr.yieldFullFileName(-1,systMode,createDir);
    std::cout << "loading from <" << yieldFullName << ">\n";
    TFile file(yieldFullName,"read");
    int res=1;
    if (res) res=checkBinningArrays(file);
    if (res) res=loadVec(file,yields,"yields");
    file.Close();
    if (!res) {
      std::cout << "error occurred during load of file <" << yieldFullName << ">\n";
      return retCodeError;
    }
    std::cout << dashline;
  }

  if (0) {
    TFile file("test_new_2D.root","recreate");
    for (unsigned int i=0; i<yields.size(); ++i) {
      yields[i]->Write(inpMgr.sampleName(i));
    }
    file.Close();
  }

  if (0) {
    for (unsigned int i=0; i<inpMgr.sampleCount(); ++i) {
      //std::cout << "sample " << inpMgr.sampleName(i) << "\n";
      //yields[i]->Print("ranges");
      printHisto(yields[i]);
    }
  }
  //return retCodeStop;
  
  bool useTrue2eBkgDataDriven = true;
  bool useFakeBkgDataDriven = true;

  HistoPair2D_t ddbkgTrue2e;
  HistoPair2D_t ddbkgFake;

  const int checkBinning=1;
  if (useTrue2eBkgDataDriven) {
    TString fname=inpMgr.GetTrue2eDDBkgFName();
    std::cout << "fname=" << fname << "\n";
    if (!ddbkgTrue2e.loadThreeMatrices(fname,
			   "true2eBackgroundFromData",
			   "true2eBackgroundFromDataError",
			   "true2eBackgroundFromDataErrorSyst",
			   checkBinning,1)) {
       std::cout << "failed to load data from <" << fname << ">\n";
       return retCodeError;
    }
    ddbkgTrue2e.print();
  }

  if (useFakeBkgDataDriven) {
    TString fname=inpMgr.GetFakeDDBkgFName();
    if (!ddbkgFake.loadThreeMatrices(fname,
			   "fakeBackgroundFromData",
			   "fakeBackgroundFromDataError",
			   "fakeBackgroundFromDataErrorSyst",
			   checkBinning,1)) {
       std::cout << "failed to load data from <" << fname << ">\n";
       return retCodeError;
    }
    ddbkgFake.print();
  }

  if (0) {
    TCanvas *cx=new TCanvas("cx","cx",600,600);
    if (ddbkgTrue2e.histo()) ddbkgTrue2e.histo()->Draw("COLZ");
    cx->Update();
  }

  // background
  HistoPair2D_t mcbkgTrue2e("mcBkgTrue2e");
  HistoPair2D_t mcbkgFake("mcBkgFake");

  TH2D *hWZZZ=NULL;
  int add_WZZZ_error_old_way=1;

  // Calculate true dielectron background, which includes
  // WW,WZ,ZZ, ttbar, Wt, and DY->tautau
  for (unsigned int i=0; i<yields.size(); ++i) {
    TString sname=inpMgr.sampleName(i);
    double systWeight=0.;
    int add=0;
    if (sname == "ztt") { systWeight=0.; add=1; }
    else if (sname == "ttbar") { systWeight=0.5; add=1; }
    else if (sname == "ww") { systWeight=1.0; add=1; }
    else if (sname == "wtop") { systWeight=1.0; add=1; }
    else if ((sname == "zz") || (sname == "wz")) { add=2; }
    else if ((sname == "qcd") || (sname == "wjets")) { systWeight=0.5; add=-1; }
    else if ((sname != "data") && (sname != "zee")) {
      std::cout << "the sample " << sname << " is not considered in the systematics\n";
      return retCodeError;
    }

    if (add==1) { // weighted error
      mcbkgTrue2e.add(yields[i]);
      if (systWeight > 0.) {
	TH2D *weights=(TH2D*)yields[i]->Clone(yields[i]->GetName() + TString("_clone"));
	removeError(weights);
	swapContentAndError(weights);
	mcbkgTrue2e.addSystErr(weights,systWeight);
	delete weights;
      }
    }
    else if (add==2) { 
      mcbkgTrue2e.add(yields[i]);
      // 100% syst error
      if (add_WZZZ_error_old_way) {
	// Previous version of DYee package (for 2011 data) 
	// had the syst error on zz and wz as
	// a linear sum of yields. Here it is added in quadrature
	if (!hWZZZ) {
	  hWZZZ=Clone(yields[i],"hWZZZ","hWZZZ");
	}
	else {
	  hWZZZ->Add(yields[i]);
	}
      }
      else {
	// Here it is added in quadrature
	TH2D *weights=(TH2D*)yields[i]->Clone(yields[i]->GetName() + TString("_clone"));
	removeError(weights);
	swapContentAndError(weights);
	mcbkgTrue2e.addSystErr(weights,1.);
	delete weights;
      }
    }
    else if (add==-1) { 
      mcbkgFake.add(yields[i]);
      if (systWeight > 0.) {
	TH2D *weights=(TH2D*)yields[i]->Clone(yields[i]->GetName() + TString("_clone"));
	removeError(weights);
	swapContentAndError(weights);
	mcbkgFake.addSystErr(weights,systWeight);
	delete weights;
      }
    }
  }

  if (add_WZZZ_error_old_way) {
    TH2D *weights=(TH2D*)hWZZZ->Clone(hWZZZ->GetName() + TString("_clone"));
    removeError(weights);
    swapContentAndError(weights);
    mcbkgTrue2e.addSystErr(weights,1.);
    delete weights;
  }

  if (0) {
    mcbkgTrue2e.print();
    mcbkgFake.print();
  }

  HistoPair2D_t observedYield("observedYield");
  HistoPair2D_t mcbkgTotal("mcBkgTotal");
  HistoPair2D_t ddbkgTotal("ddBkgTotal");

  observedYield.add(yields[0]);

  mcbkgTotal.add(mcbkgTrue2e);
  mcbkgTotal.add(mcbkgFake);

  if (useTrue2eBkgDataDriven && useFakeBkgDataDriven) {
    ddbkgTotal.add(ddbkgTrue2e);
    ddbkgTotal.add(ddbkgFake);
  }

  // Construct signalYield by first assigning the observedYield
  HistoPair2D_t signalYieldMCbkg("signalYieldMCbkg",observedYield);
  HistoPair2D_t signalYieldDDbkg("signalYieldDDbkg",observedYield);

  signalYieldMCbkg.add_allErrorIsSyst(mcbkgTotal, -1.);

  int negValSigMCbkg=signalYieldMCbkg.correctNegativeValues();
  std::cout << "corrected " << negValSigMCbkg << " negative entries in the signal after MCBkg\n";

  if (useTrue2eBkgDataDriven && useFakeBkgDataDriven) {
    signalYieldDDbkg.add_allErrorIsSyst(ddbkgTotal,  -1.);
    int negValSigDDbkg=signalYieldDDbkg.correctNegativeValues();
    std::cout << "corrected " << negValSigDDbkg << " negative entries in the signal after DDBkg\n";
  }
  else {
    // remove assigned observedYieldValues
    signalYieldDDbkg.Reset();
  }


  /*
  TMatrixD bkgRatesUsual(DYTools::nMassBins,nYBinsMax);
  for (int i=0; i<DYTools::nMassBins; i++) { 
    for (int j=0; j<DYTools::nYBins[i]; j++) { 
      double denom=observedYields(i,j);
      if (denom==0.0) bkgRatesUsual(i,j)=0;
      bkgRatesUsual(i,j)= 100.0*totalBackground(i,j)/denom;
    } 
  } 
  */


  //
  // data-MC residual difference: 
  // create weight factors to make MC prediction signal-like
  //
  //HistoPair2D_t mcSignal("mcSignal");
  //mcSignal.cloneHisto(yields.back(),1);
  TH2D *mcSignal=Clone(yields.back(),"mcSignal",1);

  // 1. determine counts in Z-peak area in data and MC
  std::cout << "Event counts in Z-peak region:\n";
  double countDataNoMCbkg=signalYieldMCbkg.ZpeakCount(NULL);
  std::cout << "\t countDataNoMCbkg=" << countDataNoMCbkg << "\n";
  double countDataNoDDbkg=signalYieldDDbkg.ZpeakCount(NULL);
  std::cout << "\t countDataNoDDbkg=" << countDataNoDDbkg << "\n";
  double countMCsignal=ZpeakCount(mcSignal,NULL);
  std::cout << "\t countMCsignal=" << countMCsignal << "\n";

  // 2. create weight maps
  TH2D* zeeMCShapeReweight_mcBkg=Clone(signalYieldMCbkg.histo(),"zeeMCShapeReweight_mcBkg",1);
  zeeMCShapeReweight_mcBkg->Divide(mcSignal);
  zeeMCShapeReweight_mcBkg->Scale(countMCsignal/countDataNoMCbkg);
  eliminateNaNs(zeeMCShapeReweight_mcBkg,1.,0.);
  std::cout << " Z peak count from reweight (MCbkg)=" << ZpeakCount(zeeMCShapeReweight_mcBkg) << "\n";

  TH2D* zeeMCShapeReweight_ddBkg=Clone(signalYieldDDbkg.histo(),"zeeMCShapeReweight_ddBkg",1);
  zeeMCShapeReweight_ddBkg->Divide(mcSignal);
  zeeMCShapeReweight_ddBkg->Scale(countMCsignal/countDataNoDDbkg);
  double setRWValue_tmp=(useTrue2eBkgDataDriven && useFakeBkgDataDriven) ? 1.:0.;
  eliminateNaNs(zeeMCShapeReweight_ddBkg,setRWValue_tmp,0.);
  std::cout << " Z peak count from reweight (DDbkg)=" << ZpeakCount(zeeMCShapeReweight_ddBkg) << "\n";


  signalYieldMCbkg.print();
  signalYieldDDbkg.print();
  zeeMCShapeReweight_mcBkg->Print("range");

  HistoPair2D_t hpMCSignal("hpMCsignal");
  if (!hpMCSignal.cloneHisto(mcSignal)) {
    std::cout << "failed to create hpMCSignal\n";
    return retCodeError;
  }

  double observedInZpeak = observedYield.ZpeakCount(NULL);
  double mcBkgInZpeak    = mcbkgTotal.ZpeakCount(NULL);
  double ddBkgInZpeak    = ddbkgTotal.ZpeakCount(NULL);
  double zeeSignalInZpeak= ZpeakCount(mcSignal,NULL);
  double ratio_data2mcAll= observedInZpeak/(mcBkgInZpeak + zeeSignalInZpeak);
  double ratio_data2ddZee= observedInZpeak/(ddBkgInZpeak + zeeSignalInZpeak);

  std::cout << "outFileName=<" << outFileName << ">\n";

  TFile fileOut(outFileName,"recreate");
  // main results 
  signalYieldMCbkg.Write();
  signalYieldDDbkg.Write();
  writeIntFlagValues("ddbkgUsed",1,int(useTrue2eBkgDataDriven && useFakeBkgDataDriven));
  writeFlagValues("ZpeakCounts",4,observedInZpeak,mcBkgInZpeak,ddBkgInZpeak,
		  zeeSignalInZpeak);
  writeFlagValues("ZpeakRatios",2,ratio_data2mcAll,ratio_data2ddZee);
  // systematics studies
  fileOut.cd();
  fileOut.mkdir("ShapeReweight");
  fileOut.cd("ShapeReweight");
  zeeMCShapeReweight_mcBkg->Write();
  zeeMCShapeReweight_ddBkg->Write();
  // input arrays
  fileOut.mkdir("Input");
  fileOut.cd("Input");
  observedYield.Write();
  mcbkgTrue2e.Write();
  mcbkgFake.Write();
  mcbkgTotal.Write();
  ddbkgTrue2e.Write();
  ddbkgFake.Write();
  ddbkgTotal.Write();
  hpMCSignal.Write();
  writeBinningArrays(fileOut);
  fileOut.Close();
  std::cout << "file <" << fileOut.GetName() << "> saved\n";
  

  /*
  TString outFileNamePlots=outFileName;
  outFileNamePlots.Replace(outFileNamePlots.Index(".root"),sizeof(".root"),"-plots.root");
  TFile *foutPlots=new TFile(outFileNamePlots,"recreate");
  PlotMatrixVariousBinning(bkgRatesUsual,"bkgRatesPercent","COLZ",foutPlots);
  PlotMatrixVariousBinning(signalYields,"signalYields","COLZ",foutPlots);
  PlotMatrixVariousBinning(signalYieldsError,"signalYieldsError","COLZ",foutPlots);
  foutPlots->Close();
  */

  /*
  if (0) {
    for(int i=0; i<DYTools::nMassBins; i++){
      // Print tables of yields and background
      if ( (DYTools::study2D==1) ||
	   ((DYTools::study2D==0) && (i==0)) ) {
	   
	printf("\n\n\t\tTables for iMass=%d\n\n",i);
	
	// Table 1: split background into true, wz/zz, and qcd
	printf(" Note: stat error in signal yield contain stat error on background,\n");
	printf("   and syst error on signal yield contains syst error on background\n");
	printf("mass range   rapidity range   observed       true2e-bg         wz-zz-bg               fake-bg                 total-bg            signal\n");
      }
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	double absYmin=0, absYmax=0;
	DYTools::findAbsYValueRange(i,yi, absYmin,absYmax);
	printf("%5.1f-%5.1f GeV ", DYTools::massBinLimits[i], DYTools::massBinLimits[i+1]);
	printf("%4.2f-%4.2f ", absYmin,absYmax);
	printf(" %7.0f+-%3.0f ", observedYield.getBinContent(i,yi), observedYield.getBinError(i,yi));
	printf(" %5.1f+-%4.1f+-%4.1f ", true2eBackground[i][yi], true2eBackgroundError[i][yi], true2eBackgroundErrorSyst[i][yi]);
	printf(" %6.2f+-%4.2f+-%4.2f ", wzzz[i][yi], wzzzError[i][yi], wzzzErrorSyst[i][yi]);
	printf(" %5.1f+-%5.1f+-%5.1f ", fakeEleBackground[i][yi], fakeEleBackgroundError[i][yi], fakeEleBackgroundErrorSyst[i][yi]);
	printf("    %5.1f+-%4.1f+-%4.1f ", totalBackground[i][yi], totalBackgroundError[i][yi], totalBackgroundErrorSyst[i][yi]);
	printf("    %8.1f+-%5.1f+-%5.1f ", signalYields[i][yi], signalYieldsError[i][yi], signalYieldsErrorSyst[i][yi]);
	printf("\n");
      }
    }
    std::cout << std::endl;
  }
  */

  if (0) {
    for (int ibackground=0; ibackground<2; ++ibackground) {
      const HistoPair2D_t *bkgTrue2e=(ibackground==0) ? &mcbkgTrue2e : &ddbkgTrue2e;
      const HistoPair2D_t *bkgFake=(ibackground==0) ? &mcbkgFake : &ddbkgFake;
      const HistoPair2D_t *totBkg= (ibackground==0) ? &mcbkgTotal : &ddbkgTotal;
      const HistoPair2D_t *sigYield=(ibackground==0) ? &signalYieldMCbkg : &signalYieldDDbkg;

      if ((ibackground==1) && (!useTrue2eBkgDataDriven || !useFakeBkgDataDriven)) {
	std::cout << "\nData-driven background is not used";
	if (( useTrue2eBkgDataDriven && !useFakeBkgDataDriven) ||
	    (!useTrue2eBkgDataDriven &&  useFakeBkgDataDriven)) std::cout << "fully";
	std::cout << "\n";
	continue;
      }

      for(int iidx=0; iidx<DYTools::nMassBins; iidx++){
	// Print tables of yields and background
	if ( (DYTools::study2D==1) ||
	     ((DYTools::study2D==0) && (iidx==0)) ) {
	   
	  printf("\n\n\t\tTables for iMass=%d, %s background\n\n",iidx,(ibackground) ? "data-driven" : "MC");
	
	  // Table 1: split background into true and fake
	  printf(" Note: stat error in signal yield contain stat error on background,\n");
	  printf("   and syst error on signal yield contains syst error on background\n");
	  printf("mass range   rapidity range   observed       true2e-bg         fake-bg                 total-bg            signal\n");
	}
	for (int yiidx=0; yiidx<DYTools::nYBins[iidx]; ++yiidx) {
	  double absYmin=0, absYmax=0;
	  DYTools::findAbsYValueRange(iidx,yiidx, absYmin,absYmax);
	  printf("%5.1f-%5.1f GeV ", DYTools::massBinLimits[iidx], DYTools::massBinLimits[iidx+1]);
	  printf("%4.2f-%4.2f ", absYmin,absYmax);
	  const int i=iidx+1;
	  const int yi=yiidx+1;
	  printf(" %7.0f+-%3.0f ", observedYield.getBinContent(i,yi), observedYield.getBinError(i,yi));
	  printf(" %5.1f+-%4.1f+-%4.1f ", bkgTrue2e->getBinContent(i,yi),bkgTrue2e->getBinError(i,yi),bkgTrue2e->getBinSystError(i,yi));
	  printf(" %5.1f+-%5.1f+-%5.1f ", bkgFake->getBinContent(i,yi),bkgFake->getBinError(i,yi),bkgFake->getBinSystError(i,yi));
	  printf("    %5.1f+-%4.1f+-%4.1f ",totBkg->getBinContent(i,yi),totBkg->getBinError(i,yi),totBkg->getBinSystError(i,yi));
	  printf("    %8.1f+-%5.1f+-%5.1f ", sigYield->getBinContent(i,yi),sigYield->getBinError(i,yi),sigYield->getBinSystError(i,yi));
	  printf("\n");
	}
      }
      std::cout << std::endl;
    }
  }

  /*
  if (0) {
    int yi=0;
    // Table 2: combined true2e and WZ/ZZ backgrounds only
    printf("\n  only true2e-bg + ww-wz\n");
    printf("mass range      true2e, includingwz/zz\n");
    for(int i=0; i<DYTools::nMassBins; i++){
      printf("%5.1f-%5.1f GeV: ", DYTools::massBinLimits[i], DYTools::massBinLimits[i+1]);
      double val = true2eBackground[i][yi] + wzzz[i][yi];
      double err=0, sys=0;
      err = sqrt(SQR(true2eBackgroundError[i][yi])
		 + SQR(wzzzError[i][yi]));
      sys = sqrt(SQR(true2eBackgroundErrorSyst[i][yi])
		 + SQR(wzzzErrorSyst[i][yi]) );
      printf(" %5.1f+-%4.1f+-%4.1f ", val,err, sys);
      printf("\n");
    }

    // Table 3: Systematic error on signal yields assuming that it includes
    // only the syst. error on the background.
    printf("\n  Systematics, %% relative to background subtracted yields\n");
    printf("mass range            subtr-signal    total-bg      syst-from-bg-frac      syst-from-bg-percent\n");
    for(int i=0; i<DYTools::nMassBins; i++){
      printf("%5.1f-%5.1f GeV: ", DYTools::massBinLimits[i], DYTools::massBinLimits[i+1]);
      printf("    %8.1f+-%5.1f+-%4.1f ", signalYields[i][yi], signalYieldsError[i][yi],signalYieldsErrorSyst[i][yi]);
      printf("    %5.1f+-%4.1f+-%4.1f ", totalBackground[i][yi], totalBackgroundError[i][yi], totalBackgroundErrorSyst[i][yi]);
      printf("    %6.4f ", totalBackgroundErrorSyst[i][yi]/signalYields[i][yi]);
      printf("    %6.1f ", totalBackgroundErrorSyst[i][yi]*100.0/signalYields[i][yi]);
      printf("\n");
    }
  }

//Latex printout
  if (DYTools::study2D==1)
     latexPrintoutBackgroundRates2D(observedYields, observedYieldsErr, 
                                  totalBackground, totalBackgroundError, 
                                  totalBackgroundErrorSyst, bkgRatesUsual, 
                                    "YieldsAndBackgrounds/subtractBackground.C");
   else if (DYTools::study2D==0)
     {
       latexPrintoutBackgroundRates1D(observedYields, observedYieldsErr, 
                                  totalBackground, totalBackgroundError, 
                                  totalBackgroundErrorSyst, bkgRatesUsual, 
                                    "YieldsAndBackgrounds/subtractBackground.C");
       if (doDDvsMCcomparisonTable) bkgTablesToLatex(true2eBackground,true2eBackgroundError,true2eBackgroundErrorSyst, 
                        fakeEleBackground, fakeEleBackgroundError, fakeEleBackgroundErrorSyst, 
                        true2eBackgroundFromData,true2eBackgroundFromDataError,true2eBackgroundFromDataErrorSyst, 
                        fakeEleBackgroundFromData, fakeEleBackgroundFromDataError, fakeEleBackgroundFromDataErrorSyst, 
                        totalBackground, totalBackgroundError, totalBackgroundErrorSyst);
     }


  //
  // debug plots
  //
  const int makeDebugPlots=1;
  if (makeDebugPlots) {
    if (0) {
      const int plot_idx=0;
      const int fnc_of_rapidity=0;
      TMatrixD tmpErr=zeeMCShapeReweight;
      tmpErr=0;      
      PlotMatrixMYSlices(plot_idx,fnc_of_rapidity,zeeMCShapeReweight,  "zeeMCShapeReweight");
    }

    if (1) {
      std::vector<int> indices;
      indices.push_back(0);
      std::vector<TMatrixD> matrV;
      std::vector<TMatrixD> matrErrV;
      std::vector<TString> labelV;
      TMatrixD tmpErr=zeeMCShapeReweight;  tmpErr=0;
      TMatrixD mcRew=zeePredictedYield;
      for (int i=0; i<zeePredictedYield.GetNrows(); ++i) {
	for (int j=0; j<zeePredictedYield.GetNcols(); ++j) {
	  mcRew(i,j) = zeePredictedYield(i,j) * zeeMCShapeReweight(i,j);
	}
      }
      matrV.push_back(signalYields); matrV.push_back(zeePredictedYield);
      matrV.push_back(mcRew);
      matrErrV.push_back(tmpErr); matrErrV.push_back(tmpErr);
      matrErrV.push_back(tmpErr);
      labelV.reserve(matrV.size());
      labelV.push_back("signalYields (data)");
      labelV.push_back("zeeYields (MC)");
      labelV.push_back("reweighted MC");
      PlotMatrixMYSlices(indices,0,matrV, matrErrV, labelV, "dataVsMC",
		       "hist", NULL, "dataVsMC");
    }

    if (1) {
      const int iYBin=0;
      const int perMassBinWidth=0;
      const int perRapidityBinWidth=0;
      TCanvas *canvZpeak=MakeCanvas("canvZpeak","canvZpeak",800,900);
      TH1F* hDataNoBkg=extractMassDependence("hDataNoBkg","", 
					     signalYields,signalYieldsError,
					     iYBin,
					     perMassBinWidth,perRapidityBinWidth);
      TH1F *hZee=extractMassDependence("hZee","",
				       zeePredictedYield,zeePredictedYieldErr,
				       iYBin,
				       perMassBinWidth,perRapidityBinWidth);
      if (!hZee) {
	std::cout << "\n\n\tError: failed to locate Zee sample\n\n";
      }
      else {
	ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,"compPlot","",
			    "mass [GeV]", "counts", "MC/data");
	std::cout << "hZee normalization factor=" << (hDataNoBkg->Integral()/hZee->Integral()) << "\n";
	hZee->Scale(hDataNoBkg->Integral()/hZee->Integral());
	hZee->SetMarkerStyle(24);
	double dy=0.2;
	cp.SetRatioYRange(1-dy,1+dy);
#ifndef _check_Zpeak
	cp.SetLogx();
#endif
	removeError(hZee);
	canvZpeak->Divide(1,2);
	cp.PreparePads(canvZpeak);
	cp.AddHist1D(hDataNoBkg, "data signal", "LPE", kBlack, 1,0,1);
	cp.AddHist1D(hZee, "MC (normalized)", "LP same", kBlue, 1,0,1);
	cp.Draw(canvZpeak,false,"png");
	SaveCanvas(canvZpeak,canvZpeak->GetName());
      }
    }

  }

*/

  // Make yield distributions
  if (1) {
    std::cout << dashline;
    std::cout << dashline;
    std::cout << "Make yield distributions\n";
    std::cout << dashline;

    int scale=1;
    for (int useDDBkg=0; useDDBkg<2; ++useDDBkg) {

      std::vector<TH2D*> histosV;
      std::vector<TString> labelsV;
      std::vector<int> colorsV;
      if (useDDBkg) {
	histosV.reserve(4);
	labelsV.reserve(4);
	colorsV.reserve(4);
	histosV.push_back(yields[0]); // data
	histosV.push_back(ddbkgTrue2e.histo());
	histosV.push_back(ddbkgFake.histo());
	histosV.push_back(yields.back()); // zeeMC
	labelsV.push_back("data");
	labelsV.push_back("ddbkg true2e");
	labelsV.push_back("ddbkg fake");
	labelsV.push_back("Zee MC");
	colorsV.push_back(kBlack);
	colorsV.push_back(kViolet);
	colorsV.push_back(kMagenta);
	colorsV.push_back(inpMgr.sampleInfos().back()->colors[0]);
      }
      else {
	histosV.reserve(yields.size());
	labelsV=inpMgr.sampleNames();
	colorsV.reserve(inpMgr.sampleCount());
	for (unsigned int i=0; i<inpMgr.sampleCount(); ++i) {
	  histosV.push_back(yields[i]);
	  colorsV.push_back(inpMgr.sampleInfo(i)->colors[0]);
	}
      }

      std::vector<std::vector<TH1D*>*> hProfV;
      TString canvName=(useDDBkg) ? "cDD" : "cMC";
      int canvWidth=(DYTools::study2D==1) ? 1100 : 700;
      TCanvas *c1=new TCanvas(canvName,canvName, canvWidth,800);

      if (DYTools::study2D==1) {

	if (!createRapidityProfileVec(histosV,hProfV,labelsV)) {
	  std::cout << "failed to create profiles for useDDBkg=" << useDDBkg << "\n";
	  return retCodeError;
	}
	std::cout << "there are " << hProfV.size() << " profiles\n";
      
	for (int im=1; im<7; ++im) {
	  
	  TString mStr=Form("M_%2.0lf_%2.0lf",DYTools::massBinLimits[im],DYTools::massBinLimits[im+1]);
	  TString cpName="cp_" + mStr;
	  TString cpTitle=mStr;
	  ComparisonPlot_t *cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,cpName,cpTitle,"|y|","uncorrected yield","ratio");
	  if (im==1) cp->Prepare6Pads(c1,1);
	  
	  TString hName=Form("hSumMC_%d",im);
	  TH1D *hSum=(TH1D*)(*hProfV[im])[0]->Clone(hName);
	  hSum->Reset();
	  hSum->SetDirectory(0);
	  hSum->SetTitle(hName);
	  
	  for (unsigned int ih=0; ih<hProfV[im]->size(); ++ih) {
	    TH1D* h=(*hProfV[im])[ih];
	    if (scale && (ih>0)) {
	      const double scaleVal=(useDDBkg) ?
		ratio_data2ddZee : ratio_data2mcAll;
	      h->Scale(scaleVal);
	    }
	    if (ih==0) {
	      h->SetMarkerStyle(20);
	      cp->AddHist1D(h,labelsV[0],"LP",colorsV[ih]);
	    }
	    else {
	      cp->AddToStack(h,labelsV[ih],colorsV[ih]);
	      hSum->Add(h,1.);
	    }
	}
	  cp->AddHist1D(hSum,"simulation","LP skip",kRed+1,1,1,-1);
	  cp->Draw6(c1,1,im);
	  if (!useDDBkg) cp->ChangeLegendPos(0.2,0.,0.,0.);
	  c1->Update();
	  
	}
      }
      else {
	// 1D

	if (!createMassProfileVec(histosV,hProfV,labelsV)) {
	  std::cout << "failed to create profiles for useDDBkg=" << useDDBkg << "\n";
	  return retCodeError;
	}
	std::cout << "there are " << hProfV.size() << " profiles\n";
      
	for (int iy=0; iy<DYTools::nYBinsMax; ++iy) {
	  
	  TString yStr=Form("iy_%d",iy);
	  TString cpName=TString("cp_") + yStr;
	  TString cpTitle; //=yStr;
	  ComparisonPlot_t *cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,cpName,cpTitle,"#it{M}_{ee} [GeV]","uncorrected yield","ratio");
	  cp->SetLogx(1);
	  //cp->SetLogy(1);
	  if (iy==0) cp->Prepare2Pads(c1);
	  
	  TString hName=Form("hSumMC_%d",iy);
	  TH1D *hSum=(TH1D*)(*hProfV[iy])[0]->Clone(hName);
	  hSum->Reset();
	  hSum->SetDirectory(0);
	  hSum->SetTitle(hName);
	  
	  for (unsigned int ih=0; ih<hProfV[iy]->size(); ++ih) {
	    TH1D* h=(*hProfV[iy])[ih];
	    if (cp->logY()) h->SetMinimum(1e-7);
	    std::cout << "h->GetMinimum()=" << h->GetMinimum() << "\n";
	    if (scale && (ih>0)) {
	      const double scaleVal=(useDDBkg) ?
		ratio_data2ddZee : ratio_data2mcAll;
	      h->Scale(scaleVal);
	    }
	    if (ih==0) {
	      h->SetMarkerStyle(20);
	      cp->AddHist1D(h,labelsV[0],"LP",colorsV[ih]);
	    }
	    else {
	      cp->AddToStack(h,labelsV[ih],colorsV[ih]);
	      hSum->Add(h,1.);
	    }
	  }
	  if (cp->logY()) hSum->SetMinimum(1e-7);
	  cp->SkipInRatioPlots((*hProfV[iy])[0]);
	  cp->AddHist1D(hSum,"simulation","LP skip",kRed+1,0,1,-1);
	  cp->SetRefIdx(hSum);
	  TH1D* hDataClone=Clone((*hProfV[iy])[0],"hDataClone","");
	  cp->AddHist1D(hDataClone,"data","LP skip",TAttMarker(kBlack,20,1),0,1,-1);
	  cp->SetLogy(0);
	  cp->Draw(c1);
	  if (!useDDBkg) cp->ChangeLegendPos(0.2,0.,0.,0.);
	  c1->Update();
	  
	}
      }
      
      TString fname="fig-Yields-";
      fname.Append((useDDBkg) ? "DDBkg" : "MCBkg");
      fname.Append(DYTools::analysisTag);
      SaveCanvas(c1,fname);
    }

  }


  //gBenchmark->Show("subtractBackground");
  ShowBenchmarkTime("subtractBackground");
  return retCodeOk;

}

/*
Bool_t checkMatrixSize(TMatrixD m)
{  
  if (m.GetNrows()==DYTools::nMassBins) {
    if ((m.GetNcols()==DYTools::findMaxYBins()) ||
	(m.GetNcols()==DYTools::nYBinsMax) ) {
      return 1;
    }
  }

  std::cout << "m.dims (" << m.GetNrows() << " x " << m.GetNcols() << ") instead of expected (" << DYTools::nMassBins << " x " << DYTools::findMaxYBins() << ") or (" << DYTools::nMassBins << " x " << DYTools::nYBinsMax << ")" << std::endl;
  return 0;
}


void unsquareMatrixElements(TMatrixD &m) {
  for (int i=0; i<m.GetNrows(); ++i) {
    for (int j=0; j<m.GetNcols(); ++j) {
      m(i,j) = sqrt(m(i,j));
    }
  }
}

void bkgTablesToLatex(TMatrixD true2eBackground, TMatrixD true2eBackgroundError, TMatrixD true2eBackgroundErrorSyst, 
                      TMatrixD fakeEleBackground, TMatrixD  fakeEleBackgroundError, TMatrixD fakeEleBackgroundErrorSyst, 
                      TMatrixD true2eBackgroundFromData, TMatrixD true2eBackgroundFromDataError, TMatrixD true2eBackgroundFromDataErrorSyst, 
                      TMatrixD fakeEleBackgroundFromData, TMatrixD  fakeEleBackgroundFromDataError, TMatrixD fakeEleBackgroundFromDataErrorSyst, 
                      TMatrixD totalBackground, TMatrixD totalBackgroundError, TMatrixD totalBackgroundErrorSyst)
  {
   TString path;
   if (DYTools::study2D==0) path="tables1D/";
   else path="tables2D/";
   gSystem->mkdir(path,kTRUE);

    TString txtFileName=path+"trueFakeTotal.txt";
    FILE* txtFile = fopen(txtFileName,"w");
   TString str=DayAndTimeTag();
   fprintf(txtFile,"%50s",str.Data());
    for(int i=0; i<DYTools::nMassBins; i++){
      // Print tables of yields and background
      if ( (DYTools::study2D==1) ||
	   ((DYTools::study2D==0) && (i==0)) ) {
	   
	fprintf(txtFile,"\n\n\t\tTables for iMass=%d\n\n",i);
	
	// Table 1: split background into true, wz/zz, and qcd
	fprintf(txtFile,"mass range     true2e-bg         fake-bg             total-bg\n");
      }
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	double absYmin=0, absYmax=0;
	DYTools::findAbsYValueRange(i,yi, absYmin,absYmax);
	fprintf(txtFile,"%5.0f-%5.0f & ", DYTools::massBinLimits[i], DYTools::massBinLimits[i+1]);
	//printf("%4.2f-%4.2f ", absYmin,absYmax);
	fprintf(txtFile," $%5.0f\\pm%4.0f\\pm%4.0f$& ", true2eBackground[i][yi], true2eBackgroundError[i][yi], true2eBackgroundErrorSyst[i][yi]);
	fprintf(txtFile," $%5.0f\\pm%4.0f\\pm%4.0f$& ", fakeEleBackground[i][yi], fakeEleBackgroundError[i][yi], fakeEleBackgroundErrorSyst[i][yi]);
	fprintf(txtFile," $%5.0f\\pm%4.0f\\pm%4.0f$", totalBackground[i][yi], totalBackgroundError[i][yi], totalBackgroundErrorSyst[i][yi]);
	fprintf(txtFile,"\\\\\n");
      }
    }
   fclose (txtFile);



    txtFileName=path+"trueFakeMCvsDD.txt";
    txtFile = fopen(txtFileName,"w");
   str=DayAndTimeTag();
   fprintf(txtFile,"%50s",str.Data());
    for(int i=0; i<DYTools::nMassBins; i++){
      // Print tables of yields and background
      if ( (DYTools::study2D==1) ||
	   ((DYTools::study2D==0) && (i==0)) ) {
	   
	fprintf(txtFile,"\n\n\t\tTables for iMass=%d\n\n",i);
	
	// Table 1: split background into true, wz/zz, and qcd
	fprintf(txtFile," Note: stat error in signal yield contain stat error on background,\n");
	fprintf(txtFile,"   and syst error on signal yield contains syst error on background\n");
	fprintf(txtFile,"mass range         MM-true2e          DD-true2e          MC-fake            DD-fake\n");
      }
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	double absYmin=0, absYmax=0;
	DYTools::findAbsYValueRange(i,yi, absYmin,absYmax);
	fprintf(txtFile,"%5.0f-%5.0f & ", DYTools::massBinLimits[i], DYTools::massBinLimits[i+1]);
	fprintf(txtFile," $%5.0f\\pm%4.0f\\pm%4.0f$& ", true2eBackground[i][yi], true2eBackgroundError[i][yi], true2eBackgroundErrorSyst[i][yi]);
	fprintf(txtFile," $%5.0f\\pm%4.0f\\pm%4.0f$", true2eBackgroundFromData[i][yi], true2eBackgroundFromDataError[i][yi], true2eBackgroundFromDataErrorSyst[i][yi]);
	fprintf(txtFile," $%5.0f\\pm%4.0f\\pm%4.0f$& ", fakeEleBackground[i][yi], fakeEleBackgroundError[i][yi], fakeEleBackgroundErrorSyst[i][yi]);
	fprintf(txtFile," $%5.0f\\pm%4.0f\\pm%4.0f$", fakeEleBackgroundFromData[i][yi], fakeEleBackgroundFromDataError[i][yi], fakeEleBackgroundFromDataErrorSyst[i][yi]);
	fprintf(txtFile,"\\\\\n");
      }
    }


   fclose (txtFile);
  }

*/
