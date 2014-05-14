#include <TBenchmark.h>
#include <TUnfold.h>
#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/EventSelector.hh"
#include "../Include/EventWeight.hh"
#include "../Include/UnfoldingMatrix.h"

//=== MAIN MACRO =================================================================================================

int useTUnfold(int analysisIs2D,
	       TString conf="defaultAdHoc",
	       DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN,
	       DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST) 
{
  // normal calculation
  gBenchmark->Start("useTUnfold");
  
  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return retCodeError;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Settings
  //==============================================================================================================

  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(conf)) return retCodeError;
  InputFileMgr_t inpMgrForYield(inpMgr);

  // Construct eventSelector, update mgr and plot directory
  TString extraTag,plotExtraTag;
  EventSelector_t evtSelectorForYield(inpMgrForYield,runMode,
				      DYTools::APPLY_ESCALE,
				      extraTag,plotExtraTag,
				      EventSelector::_selectDefault);
  EventSelector_t evtSelector(inpMgr,runMode,systMode,
			      extraTag, "", EventSelector::_selectDefault);
  evtSelector.setTriggerActsOnData(false);

  TString yieldName=(1 && (inpMgr.userKeyValueAsInt("DDBKG")==1)) ?
    "signalYieldDDbkg" : "signalYieldMCbkg";
  // signal yield
  HistoPair2D_t hpSignalYield(yieldName);

  // load signal yield
  DYTools::TSystematicsStudy_t yieldSystMode=DYTools::APPLY_ESCALE;
  const int loadNormalRunSelection=1;
  TString fnameBgSubtracted=
       inpMgrForYield.signalYieldFullFileName(yieldSystMode,
					      loadNormalRunSelection);
  if (!hpSignalYield.Load(fnameBgSubtracted,1)) {
    std::cout << "failed to load signal yield\n";
    return retCodeError;
  }

  // Prepare our unfolding matrix
  TString constDir= inpMgr.constDir(systMode,0);
  TString fnameTag= UnfoldingMatrix_t::generateFNameTag(systMode,0);
  UnfoldingMatrix_t detResponse(UnfoldingMatrix::_cDET_Response,"detResponse");
  if (!detResponse.autoLoadFromFile(constDir,fnameTag)) {
    std::cout << "failed to get the unfolding matrix\n";
    return retCodeError;
  }

  // unfolded yield
  HistoPair2D_t hpUnfoldedYield("unfYield");
  if (!hpUnfoldedYield.unfold(*detResponse.getDetInvResponse(),hpSignalYield)){
    std::cout << "failed to perform unfolding\n";
    return retCodeError;
  }

  //printHisto(hpUnfoldedYield);

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code
  //==============================================================================================================

  std::cout << mainpart;

  TH1D* h1Sig=createProfileX(hpSignalYield.histo(),1,"h1Sig",1);
  TH1D* h1UnfOur=createProfileX(hpUnfoldedYield.histo(),1,"h1UnfOur",1);

  TH1D *h1UnfChk= (TH1D*)h1UnfOur->Clone("h1UnfChk");
  if (!unfold_reco2true(h1UnfChk,detResponse,h1Sig)) {
    std::cout << "failed to produce h1UnfChk\n";
    return retCodeError;
  }

  //printHisto(h1UnfOur);
  printHisto(h1UnfChk);

  // For test create histogram, ignoring the mass bins,etc.
  TH2D* h2Mig= createHisto2D(*detResponse.getMigration(),
			     detResponse.getMigrationErr(),
			     "h2migration","detMigration");

  // initialize the unfolding object
  TUnfold unfold(h2Mig, TUnfold::kHistMapOutputHoriz, TUnfold::kRegModeNone);

  // set input distribution and bias scale (=0)
  if(unfold.SetInput(h1Sig,0.0)>=10000) {
    std::cout<<"Unfolding result may be wrong\n";
  }

  // set up a bin map, excluding underflow and overflow bins
  // the binMap relates the output of the unfolding to the final
  // histogram bins
  Int_t *binMap=new Int_t[h1UnfOur->GetNbinsX()+2];
  for(Int_t i=0; i<h1UnfOur->GetNbinsX()+2; i++) binMap[i]=i;
  //binMap[0]=-1; // underflow

  TH1D *h1Unf=(TH1D*)h1UnfOur->Clone("h1Unf");
  h1Unf->Reset();
  unfold.GetOutput(h1Unf,binMap);

  HERE("aa");

  printHisto(h1Unf);

  // clean-up
  delete [] binMap;

  ShowBenchmarkTime("useTUnfold");
  return retCodeOk;
}
