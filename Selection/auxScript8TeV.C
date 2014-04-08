{
  gROOT->ProcessLine(".L selectEvents.C+");
  // The first number tells if the analysis is 2D
  //selectEvents(1,"default",DYTools::DEBUG_RUN,DYTools::NO_SYST);
  //selectEvents(1,"default",DYTools::NORMAL_RUN,DYTools::NO_SYST);
  selectEvents(1,"default",DYTools::NORMAL_RUN,DYTools::LOWER_ET_CUT);
}
