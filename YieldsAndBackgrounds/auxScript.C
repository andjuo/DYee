{
  gROOT->ProcessLine(".L prepareYields.C+");
  //prepareYields("../config_files/data_vilnius8TeV_regSSD.conf.py",DYTools::NORMAL_RUN,DYTools::NO_SYST);
  //prepareYields(1,"default",DYTools::LOAD_DATA,DYTools::NO_SYST);
  //prepareYields(1,"default",DYTools::LOAD_DATA,DYTools::UNREGRESSED_ENERGY);

  // assuming selectEvents was run with systMode=LOWER_ET_CUT, apply cuts here
  prepareYields(1,"default",DYTools::NORMAL_RUN,DYTools::ESCALE_DIFF_0000);
}
