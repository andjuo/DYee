#include "CovariantEff.h"

int extractOldSF=0; // from DMDY package

void plotEffs(int debugMode, TString confFileName="../config_files/data_vilnius8TeV_regSSD.conf.py", DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST) {
  std::cout << "debugMode=" << debugMode << " is ignored\n";

  CovariantEffMgr_t mgr;
  int nExps=1;

  if (extractOldSF) {
    confFileName="../config_files/data_DMDY8TeV_oldNtuples.conf.py";
    mgr.editMgr().rootFileBaseDir("../../DMDY-Ilya-20130808-my-copy/root_files/");
  }

  HERE("calling setup");
  assert(mgr.Setup(confFileName,nExps,systMode));
  assert(mgr.initOk());

  std::cout << "\n\nok. Start studies\n";

  TString tnpString=TString(mgr.mgr().getTNP_etetaBinningString().c_str());
  tnpString.ReplaceAll(" ","");
  TString plotFName=Form("el-effs-%s.root",tnpString.Data());
  if (extractOldSF) plotFName.ReplaceAll(".root","-summer2013.root");
  TFile faPlots(plotFName,"recreate");
  drawEfficiencies(&faPlots);
  drawScaleFactors(&faPlots,1); // save arrays
  faPlots.Close();
  std::cout << "file " << faPlots.GetName() << " saved.\n";
}

