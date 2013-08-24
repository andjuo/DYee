{
  gROOT->ProcessLine(".L plotDYEfficiency.C+");
  plotDYEfficiency("../config_files/data_vilnius7TeV_test.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);
  //plotDYEfficiency("../config_files/data_vilnius7TeV_test.conf.py",DYTools::DEBUG_LOAD,DYTools::NO_SYST);

  //gROOT->ProcessLine(".L plotDYEfficiencyChk.C+");
  //plotDYEfficiencyChk("../config_files/data_vilnius7TeV_test.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);
}
