{
  gROOT->ProcessLine(".L calcCrossSection.C+");
  calcCrossSection("../config_files/data_vilnius8TeV.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);
}
