#define InputFileMgr_defaultFile_HH


TString getDefaultFileName(TString index) {
  TString fname;
  int replaced=1;
  if (index=="defaultRemote") {
    fname="../config_files/data_vilnius8TeVremote_regSSD.conf.py";
  }
  else if (index=="default") {
    fname="../config_files/data_vilnius8TeV_regSSD.conf.py";
  }
  else if (index=="defaultEgamma") {
    fname="../config_files/data_vilnius8TeV_regSSD_egamma.conf.py";
  }
  else if (index=="defaultAdHocRemote") {
    fname="../config_files/data_vilnius8TeVremote_regSSD_adHoc.conf.py";
  }
  else if ((index=="cernRemoteAdHoc") || (index=="cernAdHocRemote")) {
    fname="../config_files/data_cern8TeVremote_adHoc.conf.py";
  }
  else if (index=="defaultAdHoc") {
    fname="../config_files/data_vilnius8TeV_regSSD_adHoc.conf.py";
  }
  else if (index=="cernDefault") {
    fname="../config_files/data_cern8TeV.conf.py";
  }
  else if (index=="cernRnd") {
    fname="../config_files/data_cern8TeV_rndStudy.conf.py";
  }
  else if (index=="cernRemote") {
    fname="../config_files/data_cern8TeVremote.conf.py";
  }
  else if (index=="defaultOld") {
    fname="../config_files/old-data_vilnius8TeV_regSSD.conf.py";
  }
  else if (index=="dressed") {
    fname="../config_files/data_vilnius8TeV_dressed.conf.py";
  }
  else { fname=index; replaced=0; }
  if (replaced) std::cout << "getDefaultFileName: "
			  << "converted to inputFileName=<" << fname << ">\n";
  return fname;
}

