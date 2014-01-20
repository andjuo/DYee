{
  gROOT->ProcessLine(".L prepareYields.C+");
  prepareYields("../config_files/data_vilnius8TeV_reg.conf.py",DYTools::NORMAL_RUN,DYTools::NO_SYST);
  //prepareYields("../config_files/data_vilnius8TeV_reg.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);

  gROOT->ProcessLine(".L subtractBackground.C+");
  subtractBackground("../config_files/data_vilnius8TeV_reg.conf.py",DYTools::NORMAL_RUN,DYTools::NO_SYST);
}
