#include "CovariantEff.h"

void plotEffs(int debugMode) {
  TString confFileName="../config_files/data_vilnius8TeV_regSSD.conf.py";

  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;
  CovariantEffMgr_t mgr;
  int nExps=1;

  HERE("calling setup");
  assert(mgr.Setup(confFileName,nExps));
  assert(mgr.initOk());

  std::cout << "\n\nok. Start studies\n";

  TString plotFName="effs.root";
  TFile faPlots(plotFName,"recreate");
  drawEfficiencies(&faPlots);
  faPlots.Close();
}

