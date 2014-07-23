/*
{
  gROOT->ProcessLine(".L plotUnfoldingMatrix.C+");
  //plotUnfoldingMatrix(1,"../config_files/data_vilnius8TeV.conf.py",DYTools::DEBUG_RUN,DYTools::NO_SYST);
  //plotUnfoldingMatrix(1,"default",DYTools::NORMAL_RUN,DYTools::NO_SYST);
  plotUnfoldingMatrix(1,"default",DYTools::DEBUG_RUN,DYTools::FSR_STUDY);

}
*/

{
  gROOT->ProcessLine(".L plotUnfoldingMatrixDressed.C+");
  plotUnfoldingMatrixDressed(0,"../config_files/data_vilnius8TeV_dressed.conf.py",DYTools::NORMAL_RUN,DYTools::NO_SYST,20.);
}
