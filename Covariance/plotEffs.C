#include "CovariantEff.h"

int extractOldSF=0;

void plotEffs(int debugMode) {
  TString confFileName="../config_files/data_vilnius8TeV_regSSD.conf.py";
  std::cout << "debugMode=" << debugMode << " is ignored\n";

  if (extractOldSF) {
    confFileName="../config_files/data_DMDY8TeV_oldNtuples.conf.py";
  }

  //DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;
  CovariantEffMgr_t mgr;
  int nExps=1;

  if (extractOldSF) {
    mgr.editMgr().rootFileBaseDir("../../DMDY-Ilya-20130808-my-copy/root_files/");
  }

  HERE("calling setup");
  assert(mgr.Setup(confFileName,nExps));
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

