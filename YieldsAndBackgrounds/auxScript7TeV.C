{
  gROOT->ProcessLine(".L prepareYields.C+");
  //prepareYields("../config_files/data_vilnius7TeV.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);
  //prepareYields("../config_files/data_vilnius7TeV.conf.py",DYTools::NORMAL_RUN,DYTools::NO_SYST);
  //prepareYields("../config_files/data_vilnius7TeV.conf.py",DYTools::LOAD_DATA,DYTools::NO_SYST);
  
  //prepareYields("../config_files/data_vilnius7TeV_test.conf.py",DYTools::NORMAL_RUN,DYTools::NO_SYST);

  gROOT->ProcessLine(".L subtractBackground.C+");
  subtractBackground("../config_files/data_vilnius7TeV_test.conf.py",DYTools::NORMAL_RUN,DYTools::NO_SYST);
		    
}
