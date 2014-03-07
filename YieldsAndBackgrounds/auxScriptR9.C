{
  gROOT->ProcessLine(".L prepareYieldsR9.C+");
  //prepareYieldsR9("../config_files/data_vilnius8TeV_regSSD.conf.py",DYTools::NORMAL_RUN,DYTools::NO_SYST);
  prepareYieldsR9("../config_files/data_vilnius8TeV_regSSD.conf.py",DYTools::LOAD_DATA,DYTools::NO_SYST);
  prepareYieldsR9("../config_files/data_vilnius8TeV_regSSD.conf.py",DYTools::LOAD_DATA,DYTools::UNREGRESSED_ENERGY);
}
