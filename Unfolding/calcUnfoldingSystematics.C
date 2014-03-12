#include <TFile.h>
#include "../Include/DYTools.hh"
#include "../Include/DYToolsUI.hh"
#include "../Include/HistoPair.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/EventSelector.hh"
#include "../Include/UnfoldingMatrix.h"
#include "../Include/ComparisonPlot.hh"


int unfoldDET_local(const HistoPair2D_t &iniYields, 
		    HistoPair2D_t &finYields,
		    const InputFileMgr_t &inpMgr,
		    DYTools::TSystematicsStudy_t systMode);

//TCanvas* plotHistos(


// --------------------------------------------------
//  Main code
// --------------------------------------------------

int calcUnfoldingSystematics(const TString conf){

  // check whether it is a calculation
  if (conf.Contains("_DebugRun_")) {
    std::cout << "calcUnfoldingSystematics: _DebugRun_ detected. Terminating the script\n";
    return retCodeOk;
  }

  //----------------------------------------
  // Settings 
  //========================================
  
  DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN;
  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;

  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(conf)) return retCodeError;
  // plotDetResponse uses escale!
  // no energy correction for this evaluation
  inpMgr.clearEnergyScaleTag();

  // Construct eventSelector, update mgr and plot directory
  EventSelector_t evtSelector(inpMgr,runMode,systMode,
			      "", "", EventSelector::_selectDefault);
  evtSelector.setTriggerActsOnData(false);

  //const int seedMin=inpMgr.userKeyValueAsInt("SEEDMIN");
  //const int seedMax=inpMgr.userKeyValueAsInt("SEEDMAX");
  //int seedFirst=0;
  //int seedLast=-1;
  //const int seedDiff=(systMode==DYTools::FSR_STUDY) ? 3 : 0; //(seedMax-seedMin+1);

  //std::cout << "seedMin..seedMax=" << seedMin << ".." << seedMax << "\n";

  //return retCodeOk;

  TString debugFName=Form("debug_unfSyst_%s.root",DYTools::analysisTag.Data());
  TFile *debugFile=NULL;
  if (1) debugFile=new TFile(debugFName,"recreate");

  //---------------------------------------
  // Main analysis code 
  //=======================================


  std::cout << mainpart;

  HistoPair2D_t signalYieldMCbkg("signalYieldMCbkg");
  HistoPair2D_t signalYieldDDbkg("signalYieldDDbkg");
  const int ignoreDebugRunFlag=0;
  TString yieldsFileName= inpMgr.signalYieldFullFileName(systMode,ignoreDebugRunFlag);
  int loadYields=1;
  if (loadYields) {
    TFile fileYields(yieldsFileName,"read");
    if (!fileYields.IsOpen()) {
      std::cout << "Failed to load yields from <" << fileYields.GetName() << ">\n";
      return retCodeError;
    }
    signalYieldMCbkg.Read();
    signalYieldDDbkg.Read();
    fileYields.Close();
  }

  int useMCbkg=1;
  HistoPair2D_t *signal=(useMCbkg) ? &signalYieldMCbkg : &signalYieldDDbkg;

  int res=1;

  if (debugFile) { debugFile->cd(); signal->Write(); }

/////////////////////////////////
//calculate central values
/////////////////////////////////

  HistoPair2D_t unfYields("unfYields");
  if (res) res=unfoldDET_local(*signal,unfYields,inpMgr,DYTools::NO_SYST);

  if (debugFile) { debugFile->cd(); unfYields.Write(); }

/////////////////////////////////
//calculate smearing systematics 
/////////////////////////////////

  /*
  int nseeds = 0;

  for(int i=seedFirst; i<=seedLast; i++){
    nseeds++;
    applyUnfoldingLocal(signalYields, unfoldedYields, 1, i, 100);
    for(int idx = 0; idx < nUnfoldingBins; idx++){
      unfoldedYieldsMean[idx] += unfoldedYields[idx];
      unfoldedYieldsSquaredMean[idx] += unfoldedYields[idx]*unfoldedYields[idx];
    }
  }


  // Final calculation of the mean and RMS for Smearing
  TVectorD unfoldingSystPercentSmear(nUnfoldingBins);
  for(int idx = 0; idx < nUnfoldingBins; idx++){
    unfoldedYieldsMean[idx] = unfoldedYieldsMean[idx]/double(nseeds);
    unfoldedYieldsSquaredMean[idx] = 
      unfoldedYieldsSquaredMean[idx]/double(nseeds);
    unfoldedYieldsRMS[idx] = 
      sqrt(unfoldedYieldsSquaredMean[idx] - 
	   unfoldedYieldsMean[idx]*unfoldedYieldsMean[idx]);
    unfoldingSystPercentSmear[idx] = 
      unfoldedYieldsRMS[idx]*100.0/unfoldedYieldsMean[idx];
  }
  */

  
/////////////////////////////////
//calculate Fsr systematics 
/////////////////////////////////

  HistoPair2D_t unfYieldsFsr5plus("unfYieldsFsr5plus");
  HistoPair2D_t unfYieldsFsr5minus("unfYieldsFsr5minus");

  if (res) res=unfoldDET_local(*signal,unfYieldsFsr5plus,inpMgr,DYTools::FSR_5plus);
  if (res) res=unfoldDET_local(*signal,unfYieldsFsr5minus,inpMgr,DYTools::FSR_5minus);

  HistoPair2D_t *unfYields_dev_Fsr5plus = 
    getDiff("unfYields_dev_Fsr5plus",unfYields,unfYieldsFsr5plus,1);
  HistoPair2D_t *unfYields_dev_Fsr5minus =
    getDiff("unfYields_dev_Fsr5minus",unfYields,unfYieldsFsr5minus,1);
  HistoPair2D_t *unfYields_diff_FSR=
    getDiff("unfYields_diff_FSR",unfYieldsFsr5plus,unfYieldsFsr5minus,1);

  if (debugFile) {
    debugFile->cd();
    unfYieldsFsr5plus.Write();
    unfYieldsFsr5minus.Write();
    unfYields_dev_Fsr5plus->Write();
    unfYields_dev_Fsr5minus->Write();
    unfYields_diff_FSR->Write();
  }

  /*
  TVectorD unfoldedYieldsFsrMax(nUnfoldingBins);
  TVectorD unfoldedYieldsFsrMin(nUnfoldingBins);
  TVectorD unfoldedYieldsFsrErr(nUnfoldingBins);
  applyUnfoldingLocal(signalYields, unfoldedYieldsFsrMax, 0, 1000, 105);
  applyUnfoldingLocal(signalYields, unfoldedYieldsFsrMin, 0, 1000, 95);

  TVectorD unfoldingSystPercentFsr(nUnfoldingBins); 

  for(int idx = 0; idx < nUnfoldingBins; idx++) {
    unfoldedYieldsFsrErr[idx] = 
      fabs(unfoldedYieldsFsrMax[idx]-unfoldedYieldsFsrMin[idx]) /
      (unfoldedYieldsFsrMax[idx]+unfoldedYieldsFsrMin[idx]);
    unfoldingSystPercentFsr[idx] = unfoldedYieldsFsrErr[idx]*100.0;
  }
  */

/////////////////////////////////
//calculate pile-up systematics 
/////////////////////////////////

  HistoPair2D_t unfYieldsPU5plus("unfYieldsPU5plus");
  HistoPair2D_t unfYieldsPU5minus("unfYieldsPU5minus");

  if (res) res=unfoldDET_local(*signal,unfYieldsPU5plus,inpMgr,DYTools::FSR_5plus);
  if (res) res=unfoldDET_local(*signal,unfYieldsPU5minus,inpMgr,DYTools::FSR_5minus);

  HistoPair2D_t *unfYields_dev_PU5plus = 
    getDiff("unfYields_dev_PU5plus",unfYields,unfYieldsPU5plus,1);
  HistoPair2D_t *unfYields_dev_PU5minus =
    getDiff("unfYields_dev_PU5minus",unfYields,unfYieldsPU5minus,1);
  HistoPair2D_t *unfYields_diff_PU=
    getDiff("unfYields_diff_PU",unfYieldsPU5plus,unfYieldsPU5minus,1);

  if (debugFile) {
    debugFile->cd();
    unfYieldsPU5plus.Write();
    unfYieldsPU5minus.Write();
    unfYields_dev_PU5plus->Write();
    unfYields_dev_PU5minus->Write();
    unfYields_diff_PU->Write();
  }


/////////////////////////////////
//calculate total systematics 
/////////////////////////////////
  /*
    TVectorD unfoldingSystPercent(nUnfoldingBins); 
    for(int idx = 0; idx < nUnfoldingBins; idx++) {
      unfoldingSystPercent[idx] = 
	sqrt(
	     unfoldingSystPercentSmear[idx]*unfoldingSystPercentSmear[idx] +
	     unfoldingSystPercentFsr[idx]*unfoldingSystPercentFsr[idx] 
	     );
    
  }
  */

//printing out to the screen

  /*
   printf("mass     rapidity   mean-unfolded   RMS-unfolded   rel-error     rel-err-percent-Smear     rel-err-percent-Fsr      rel-err-percent-total \n");
   for(int i=0, idx=0; i<DYTools::nMassBins; i++){
    double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi, ++idx) {
      printf("%4.0f-%4.0f  %4.2lf-%4.2lf  %7.1f      %7.1f          %6.4f             %6.1f                 %6.3f                       %6.1f\n",
	     DYTools::massBinLimits[i],DYTools::massBinLimits[i+1],
	     rapidityBinLimits[yi],rapidityBinLimits[yi+1],
	     unfoldedYieldsMean[idx], unfoldedYieldsRMS[idx],
	     unfoldedYieldsRMS[idx]/unfoldedYieldsMean[idx],
	     unfoldingSystPercentSmear[idx], unfoldingSystPercentFsr[idx], 
	     unfoldingSystPercent[idx]);
      //unfoldedYieldsRMS[idx]*100.0/unfoldedYieldsMean[idx]);
    }
  }
  


  // Store constants in the file
  TString outputDir(TString("../root_files/systematics/")+tagDirConstants);
  gSystem->mkdir(outputDir,kTRUE);
  TString unfoldingSystFileName(outputDir+TString("/unfolding_systematics") 
				+ DYTools::analysisTag + TString(".root"));

  TFile fa(unfoldingSystFileName,"recreate");
  unfolding::writeBinningArrays(fa);
  unfoldedYieldsMean.Write("unfoldedYieldsMeanFI");
  unfoldedYieldsRMS.Write("unfoldedYieldsRMSFI");
  unfoldingSystPercentSmear.Write("unfoldingSystPercentSmearFI");
  unfoldingSystPercentFsr.Write("unfoldingSystPercentFsrFI");
  unfoldingSystPercent.Write("unfoldingSystPercentFI");
  fa.Close();

  */
  if (debugFile) { 
    debugFile->Close();
    std::cout << "created a debug file: <" << debugFile->GetName() << ">\n";
    delete debugFile; 
  }

  return retCodeOk;
}

//-----------------------------------------------------------------
// Unfold
//-----------------------------------------------------------------

int unfoldDET_local(const HistoPair2D_t &iniYields, 
		 HistoPair2D_t &finYields,
		 const InputFileMgr_t &inpMgr,
		 DYTools::TSystematicsStudy_t systMode) 
{

  TString wStr;
  const TString u="_";
  switch (systMode) {
  case DYTools::FSR_5plus: wStr="_105"; break;
  case DYTools::FSR_5minus: wStr="_095"; break;
  case DYTools::PILEUP_5plus: wStr="_PU5plus"; break;
  case DYTools::PILEUP_5minus: wStr="_PU5minus"; break;
  default: ;
  }
  TString umName=TString("detResponse") + wStr;

  DYTools::TSystematicsStudy_t dirSystMode=DYTools::NO_SYST;
  switch(systMode) {
  case DYTools::FSR_5plus:
  case DYTools::FSR_5minus:
    dirSystMode= DYTools::FSR_STUDY;
    break;
  case DYTools::PILEUP_5plus:
  case DYTools::PILEUP_5minus:
    dirSystMode= DYTools::PU_STUDY;
    break;
  default: ;
  }
  
  UnfoldingMatrix_t detResponse(UnfoldingMatrix::_cDET_Response,umName);
  TString constDir=inpMgr.constDir(dirSystMode,0);
  TString fnameTag=UnfoldingMatrix_t::generateFNameTag(dirSystMode);
  int res=detResponse.autoLoadFromFile(constDir,fnameTag);
  if (!res) res=-1;
  if (res==1) res=finYields.unfold(*detResponse.getDetInvResponse(), iniYields);

  if (res!=1) {
    std::cout << "ERROR in unfoldDET_local:\n  ";
    if (res==-1) std::cout << "failed to load unfolding matrix\n";
    else std::cout << "error during unfolding\n";
    std::cout << "for systMode=" 
	      << SystematicsStudyName(systMode)
	      << ", derived dirSystMode="
	      << SystematicsStudyName(dirSystMode)
	      << "\n";
    return 0;
  }
  return res;
}
