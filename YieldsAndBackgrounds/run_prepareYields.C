#include "prepareYields.C"

int run_prepareYields(int analysisIs2D,
		      int the_case, int systModeFlag=0, int loadData=1) {
  TString confName="default";
  DYTools::TRunMode_t runMode=(loadData) ? 
    DYTools::LOAD_DATA : DYTools::NORMAL_RUN;

  // Default systematics will be ESCALE_DIFF_0000 to apply the selection cuts
  //DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;
  DYTools::TSystematicsStudy_t systMode=DYTools::ESCALE_DIFF_0000;

  switch(systModeFlag) {
  case 0: ; break;
  case 1: systMode=DYTools::UNREGRESSED_ENERGY; break;
  }

  switch(the_case) {
  case 0: ; break;
  case 1: 
    confName="defaultAdHoc";
    systMode=DYTools::APPLY_ESCALE;
    break;
  default:
    std::cout << "the_case=" << the_case << " is not ready\n";
    return retCodeError;
  }

  // Setup the mass-rapidity setting
  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return retCodeError;
  }

  analysisIs2D=-111; // further calls to DYTools::setup do nothing
  int res=prepareYields(analysisIs2D,confName,runMode,systMode);
  return res;
}
