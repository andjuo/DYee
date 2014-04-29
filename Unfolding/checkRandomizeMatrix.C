#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include <TRandom.h>
#include "../Include/HistoPair.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/EventSelector.hh"
#include "../Include/UnfoldingMatrix.h"


// --------------------------------------------------------
// --------------------------------------------------------

int checkRandomizeMatrix(int analysisIs2D=0,
			 TString conf="defaultAdHocRemote",
			 int fsr=1) {

  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return retCodeError;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================

  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(conf)) return retCodeError;

  DYTools::TRunMode_t runMode= DYTools::NORMAL_RUN;
  DYTools::TSystematicsStudy_t systMode= DYTools::NO_SYST;

  TString extraTag;
  TString plotExtraTag;

  EventSelector_t evtSelector(inpMgr,runMode,systMode,
		       extraTag, plotExtraTag, EventSelector::_selectDefault);

  TString constDirDef=inpMgr.constDir(systMode,0);
  TString fnameTagDef=UnfoldingMatrix_t::generateFNameTag(systMode);

  UnfoldingMatrix::TUnfoldingMatrixType_t kind;
  TString uName, uNameExact;
  switch(fsr) {
  case 0:
    kind= UnfoldingMatrix::_cDET_Response;
    uName="detResponse";
    uNameExact="detResponseExact";
    break;
  case 1:
    kind= UnfoldingMatrix::_cFSR;
    uName="fsrGood";
    uNameExact="fsrExact";
    break;
  case 2:
    kind= UnfoldingMatrix::_cFSR_DET;
    uName="fsrDET";
    uNameExact="fsrDETexact";
    break;
  default:
    std::cout << "unknown fsr flag = " << fsr << "\n";
    return retCodeError;
  }

  UnfoldingMatrix_t U(kind,uName);
  UnfoldingMatrix_t Uexact(kind,uNameExact);

  if (!U.autoLoadFromFile(constDirDef,fnameTagDef) ||
      !Uexact.autoLoadFromFile(constDirDef,fnameTagDef)) {
    std::cout << "failed to get the response matrices\n";
    return retCodeError;
  }


  return retCodeOk;
}

// --------------------------------------------------------
// --------------------------------------------------------
// --------------------------------------------------------
// --------------------------------------------------------
