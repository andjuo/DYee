{
  gROOT->ProcessLine(".L plotUnfoldingMatrix.C+");
  //plotUnfoldingMatrix("../config_files/data_vilnius8TeV.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);
  plotUnfoldingMatrix("../config_files/data_vilnius8TeV.conf.py",DYTools::NORMAL_RUN,DYTools::NO_SYST);
  //plotUnfoldingMatrix("../config_files/data_vilnius8TeV.conf.py",DYTools::DEBUG_RUN,DYTools::FSR_STUDY,1.05,1.);

}
