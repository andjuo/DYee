#define InputFileMgr_defaultFile_HH


TString getDefaultFileName(TString index) {
  TString fname;
  if (index=="defaultRemote") {
    fname="../config_files/data_vilnius8TeVremote_regSSD.conf.py";
  }
  else if (index=="default") {
    fname="../config_files/data_vilnius8TeV_regSSD.conf.py";
  }
  if (index=="defaultAdHocRemote") {
    fname="../config_files/data_vilnius8TeVremote_regSSD_adHoc.conf.py";
  }
  else if (index=="defaultAdHoc") {
    fname="../config_files/data_vilnius8TeV_regSSD_adHoc.conf.py";
  }
  else if (index=="defaultOld") {
    fname="../config_files/old-data_vilnius8TeV_regSSD.conf.py";
  }
  else fname=index;
  return fname;
}

