{
  gROOT->ProcessLine(".L calcEventEff.C+");
  //calcEventEff("../config_files/data_vilnius8TeV.conf.py",1,DYTools::NORMAL_RUN);
  calcEventEff("../config_files/data_vilnius8TeV.conf.py",0,DYTools::NORMAL_RUN);
  //calcEventEff("../config_files/data_vilnius8TeV.conf.py",0,DYTools::DEBUG_RUN);

}
