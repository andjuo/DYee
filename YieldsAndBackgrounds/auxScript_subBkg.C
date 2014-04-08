#include <TROOT.h>
#include "../Include/DYTools.hh"
#include "subtractBackground.C"

int auxScript_subBkg(int analysisIs2D) {
  //gROOT->ProcessLine(".L prepareYields.C+");
  //prepareYields("../config_files/data_vilnius8TeV_regSSD.conf.py",DYTools::NORMAL_RUN,DYTools::NO_SYST);
  //prepareYields("../config_files/data_vilnius8TeV_reg.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);

  //gROOT->ProcessLine(".L subtractBackground.C+");
  //int res=subtractBackground(analysisIs2D,"default",DYTools::NORMAL_RUN,DYTools::NO_SYST);
  int res=subtractBackground(analysisIs2D,"default",DYTools::NORMAL_RUN,DYTools::ESCALE_DIFF_0000);
  return res;
}
