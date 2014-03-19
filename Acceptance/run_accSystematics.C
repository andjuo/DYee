#include "plotDYAcceptance.C"
#include "../Include/calcCorrectionSyst.h"

int run_accSystematics(int debug, TString flagStr="10011") {
  TString conf="../config_files/data_vilnius8TeV_regSSD.conf.py";

  DYTools::TRunMode_t runMode=DebugInt2RunMode(debug);

  ApplySystFlags_t flagsPlus,flagsMinus;
  if (!flagsPlus.assignFlagsOdd(flagStr) ||
      !flagsMinus.assignFlagsEven(flagStr)) return retCodeError;

  flagsPlus.adjustForAcceptance();
  flagsMinus.adjustForAcceptance();

  int ok=retCodeOk;
  if (!DYTools::loadData(runMode)) {
    std::cout << "creating distributions\n";
    if ((ok==retCodeOk) && flagsPlus.noSyst()   ) ok=plotDYAcceptance(conf,runMode,DYTools::NO_SYST);
    if ((ok==retCodeOk) && flagsPlus.fsr() ) ok=plotDYAcceptance(conf,runMode,DYTools::FSR_5plus);
    if ((ok==retCodeOk) && flagsMinus.fsr()) ok=plotDYAcceptance(conf,runMode,DYTools::FSR_5minus);
    if (ok!=retCodeOk) std::cout << "error in run_accSystematics\n";
  }
  else {
    std::vector<TH2D*> resHistos;
    int printTable=1;
    int save=1;

    ok=calcCorrectionSyst(debug,conf,"acceptance","hAcceptance",flagsPlus,flagsMinus,printTable,&resHistos,&save);
  }

  if (ok!=retCodeOk) std::cout << "error in run_accSystematics\n";
  return ok;
}
