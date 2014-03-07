{
  gROOT->ProcessLine(".L selectEventsR9.C+");
  //selectEvents("../config_files/data_vilnius7TeV.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);
  //selectEvents("../config_files/data_vilnius7TeV.conf.py",DYTools::NORMAL,DYTools::NO_SYST);
  //selectEventsR9("../config_files/data_vilnius8TeV_regSSD.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);
  selectEventsR9("../config_files/data_vilnius8TeV_regSSD.conf.py",DYTools::NORMAL_RUN,DYTools::NO_SYST);
}
