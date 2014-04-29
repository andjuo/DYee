{
  gROOT->ProcessLine(".L calcEventEff.C+");
  //calcEventEff(0,"default",1,DYTools::NORMAL_RUN);
  calcEventEff(0,"default",0,DYTools::NORMAL_RUN);
  //calcEventEff(0,"default",0,DYTools::DEBUG_RUN);
}
