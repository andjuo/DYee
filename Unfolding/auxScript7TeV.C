{
  gROOT->ProcessLine(".L plotDetResponse.C+");
  //calcDetResponse("../config_files/data_vilnius7TeV_test.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);
  plotDetResponse("../config_files/data_vilnius7TeV_test.conf.py",DYTools::DEBUG_RUN,DYTools::FSR_STUDY);
  //plotDYAcceptance("../config_files/data_vilnius7TeV_test.conf.py",DYTools::DEBUG_RUN,DYTools::FSR_STUDY,1.05,1.);

}
