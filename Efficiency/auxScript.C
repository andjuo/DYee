{
  gROOT->ProcessLine(".L plotDYEfficiency.C+");
  //plotDYEfficiency(1,"../config_files/data_vilnius8TeV.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);
  plotDYEfficiency(1,"default",DYTools::NORMAL_RUN,DYTools::NO_SYST);
  //plotDYEfficiency(1,"../config_files/data_vilnius8TeV.conf.py",DYTools::DEBUG_RUN,DYTools::FSR_STUDY,1.05,1.);

}
