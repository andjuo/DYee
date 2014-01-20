{
  gROOT->ProcessLine(".L selectEvents.C+");
  //selectEvents("../config_files/data_vilnius7TeV.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);
  selectEvents("../config_files/data_vilnius8TeV_reg.conf.py",DYTools::NORMAL_RUN,DYTools::NO_SYST);
  //selectEvents("../config_files/data_vilnius8TeV_reg.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);

}
