#define InputFileMgr_defaultFile_HH


TString getDefaultFileName(TString index) {
  TString fname;
  if (index=="default") {
    fname="../config_files/data_vilnius8TeV_regSSD.conf.py";
  }
  else if (index=="defaultOld") {
    fname="../config_files/old-data_vilnius8TeV_regSSD.conf.py";
  }
  return fname;
}

