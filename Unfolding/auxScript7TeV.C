{
  gROOT->ProcessLine(".L plotUnfoldingMatrix.C+");
  plotUnfoldingMatrix("../config_files/data_vilnius7TeV_test.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);
  //plotUnfoldingMatrix("../config_files/data_vilnius7TeV_test.conf.py",DYTools::NORMAL_RUN,DYTools::NO_SYST);
  //plotUnfoldingMatrix("../config_files/data_vilnius7TeV_test.conf.py",DYTools::LOAD_DATA,DYTools::NO_SYST);

  //plotUnfoldingMatrix("../config_files/data_vilnius7TeV_test.conf.py",DYTools::DEBUG_RUN,DYTools::FSR_STUDY);
  //plotUnfoldingMatrix("../config_files/data_vilnius7TeV_test.conf.py",DYTools::DEBUG_RUN,DYTools::FSR_STUDY,1.05,1.);

}
