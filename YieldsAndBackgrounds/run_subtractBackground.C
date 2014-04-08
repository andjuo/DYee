#include "subtractBackground.C"

int run_subtractBackground(int analysisIs2D,
			   int the_case, int systModeFlag=0, int loadData=1) {
  TString confName="default";

  DYTools::TRunMode_t runMode=(loadData) ? DYTools::LOAD_DATA : DYTools::NORMAL_RUN;
  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;

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

  int res=subtractBackgroundR9(confName,runMode,systMode);
  return res;
}
