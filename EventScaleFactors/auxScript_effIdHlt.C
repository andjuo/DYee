{
  gROOT->ProcessLine(".L eff_IdHlt.C+");
  //eff_IdHlt("../config_files/data_vilnius8TeV_regSSD-tagSyst.conf.py","HLTleg1",0,DYTools::NORMAL_RUN,DYTools::PILEUP_5minus);
  eff_IdHlt("../config_files/data_vilnius8TeV_regSSD-tagSyst.conf.py","HLTleg1",0,DYTools::DEBUG_RUN,DYTools::UNREGRESSED_ENERGY);

}
