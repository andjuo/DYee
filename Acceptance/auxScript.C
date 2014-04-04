{
  gROOT->ProcessLine(".L plotDYAcceptance.C+");
  //plotDYAcceptance(1,"../config_files/data_vilnius8TeV.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);
  plotDYAcceptance(1,"default",DYTools::NORMAL_RUN,DYTools::NO_SYST);
  //plotDYAcceptance(1,"default",DYTools::DEBUG_RUN,DYTools::FSR_STUDY,1.05,1.);
}
