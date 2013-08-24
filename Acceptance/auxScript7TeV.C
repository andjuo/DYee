{
  gROOT->ProcessLine(".L plotDYAcceptance.C+");
  //plotDYAcceptance("../config_files/data_vilnius7TeV_test.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);
  plotDYAcceptance("../config_files/data_vilnius7TeV_test.conf.py",DYTools::DEBUG_RUN,DYTools::FSR_STUDY,1.05,1.);

}
