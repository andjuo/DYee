{
  gROOT->ProcessLine(".L calcCrossSection.C+");
  calcCrossSection("../config_files/data_vilnius7TeV_test.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);
}
