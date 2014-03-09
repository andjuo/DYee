#include "plotDYAcceptance.C"

int run_accSystematics(int debug, TString flags="-") {
  TString conf="../config_files/data_vilnius8TeV_regSSD.conf.py";

  DYTools::TRunMode_t runMode=DebugInt2RunMode(debug);

  int manyFlags=(flags.Length()==5) ? 1:0;
  int noSyst=1;
  if (manyFlags) noSyst=(flags[0]=='1') ? 1:0;
  int fsr5plus= (manyFlags && (flags[1]=='1')) ? 1:0;
  int fsr5minus=(manyFlags && (flags[2]=='1')) ? 1:0;
  int pu5plus=  (manyFlags && (flags[3]=='1')) ? 1:0;
  int pu5minus= (manyFlags && (flags[4]=='1')) ? 1:0;

  int ok=retCodeOk;
  if ((ok==retCodeOk) && noSyst   ) ok=plotDYAcceptance(conf,runMode,DYTools::NO_SYST);
  if ((ok==retCodeOk) && fsr5plus ) ok=plotDYAcceptance(conf,runMode,DYTools::FSR_5plus);
  if ((ok==retCodeOk) && fsr5minus) ok=plotDYAcceptance(conf,runMode,DYTools::FSR_5minus);
  if ((ok==retCodeOk) && pu5plus  ) ok=plotDYAcceptance(conf,runMode,DYTools::PILEUP_5plus);
  if ((ok==retCodeOk) && pu5minus ) ok=plotDYAcceptance(conf,runMode,DYTools::PILEUP_5minus);
  if (ok!=retCodeOk) std::cout << "error in run_accSystematics\n";
  return ok;
}
