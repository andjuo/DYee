{
  gROOT->ProcessLine(".L selectEvents.C+");
  //selectEvents("../config_files/data_vilnius7TeV.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);
  //selectEvents("../config_files/data_vilnius7TeV.conf.py",DYTools::NORMAL,DYTools::NO_SYST);
  selectEvents("../config_files/data_vilnius7TeV_test.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);

}
