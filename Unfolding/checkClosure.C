#include "../Include/DYTools.hh"
#include "../Include/UnfoldingMatrix.h"
#include "../Include/HistoPair.hh"


void checkClosure() {
 TString dir="../root_files7TeV/constants/DY_20130801_test_DebugRun/";
#ifdef DYee8TeV
 dir="../root_files/constants/DY_j22_19789pb/";
 dir="../root_files/constants/DY_j22_19789pb_DebugRun/";
#endif
 UnfoldingMatrix_t *M=NULL;

 UnfoldingMatrix::TUnfoldingMatrixType_t kind=UnfoldingMatrix::_cFSR;
 M=new UnfoldingMatrix_t(kind, "fsrExact");
 //M=new UnfoldingMatrix_t(kind, "fsrGood");
 kind=UnfoldingMatrix::_cDET_Response;
 M=new UnfoldingMatrix_t(kind, "detResponseExact");
 //M=new UnfoldingMatrix_t(kind, "detResponse");

 M->autoLoadFromFile(dir,DYTools::analysisTag);
 
 TVectorD ini= * (M->getIniVec());
 TVectorD fin= * (M->getFinVec());
 TMatrixD iniM = * (M->getIniM());
 TMatrixD finM = * (M->getFinM());

 TVectorD chk(fin);
 //TVectorD diff(fin);
 chk.Zero();

 std::cout << "ini="; ini.Print();

 //unfold_true2reco(chk,*M,fin);
 //std::cout << "true2reco  "; chk.Print();

 unfold_reco2true(chk,*M,fin);
 //std::cout << "reco2true "; chk.Print();
 testMaxDiff("\n reco2true ",ini,chk);
 testMaxRelDiff("\n reco2true ",ini,chk);
 //diff=ini-chk;
 //std::cout << "diff "; diff.Print();

 std::cout << dashline;

 std::cout << "fin="; fin.Print();
 chk.Zero();
 unfold_true2reco(chk,*M,ini);
 //std::cout << "true2reco "; chk.Print();
 testMaxDiff("\n true2reco ",fin,chk);
 testMaxRelDiff("\n true2reco ",fin,chk);

 std::cout << "\n\n";

 {
   TMatrixD chkM(iniM);
   chkM.Zero();

   std::cout << dashline <<dashline << "\n";
   //std::cout << "iniM="; iniM.Print();
   
   unfold_reco2true(chkM,*M,finM);
   //std::cout << "reco2true "; chkM.Print();
   testMaxDiff("\n reco2true(matrix) ",iniM,chkM);
   testMaxRelDiff("\n reco2true(matrix) " ,iniM,chkM);

   std::cout << dashline;
   
   //std::cout << "finM="; finM.Print();
   chkM.Zero();
   unfold_true2reco(chkM,*M,iniM);
   //std::cout << "true2reco "; chkM.Print();
   testMaxDiff("\n true2reco(matrix) ",finM,chkM);
   testMaxRelDiff("\n true2reco(matrix) ",finM,chkM);
 }


 std::cout << "\n\n";

 if (0) {
   TMatrixD chkM(iniM);
   std::cout << "\n\n" << dashline << dashline << "\n";

   TMatrixD *unfTrue2Reco=UnfoldingMatrix_t::loadUnfM(M->getName(),dir,DYTools::analysisTag,0,1);
   TMatrixD *unfReco2True=UnfoldingMatrix_t::loadUnfM(M->getName(),dir,DYTools::analysisTag,1,1);

   unfold(chkM, *unfReco2True, finM);
   testMaxDiff("\n reco2true(matrix|load) ",iniM,chkM);

   unfold(chkM, *unfTrue2Reco, iniM);
   testMaxDiff("\n true2reco(matrix|load) ",finM,chkM);
   
   delete unfTrue2Reco;
   delete unfReco2True;

   std::cout << "\n\n" << dashline << dashline << "\n";
 }

 if (0) {
   HistoPair2D_t iniHP("iniHP"), finHP("finHP"), chkHP("chkHP");
   std::cout << "\n\n" << dashline << "\tCheck histoPair\n" << dashline << "\n";

   TMatrixD chkM(iniM);
   chkM.Zero();
   iniHP.assign(iniM,chkM,chkM);
   finHP.assign(finM,chkM,chkM);

   TMatrixD *unfTrue2Reco=UnfoldingMatrix_t::loadUnfM(M->getName(),dir,DYTools::analysisTag,0,1);
   TMatrixD *unfReco2True=UnfoldingMatrix_t::loadUnfM(M->getName(),dir,DYTools::analysisTag,1,1);

   unfold(chkHP, *unfReco2True, finHP);
   TMatrixD *resMtrue=chkHP.histoAsM();
   testMaxDiff("\n reco2true(histoPair|load) ",iniM,*resMtrue,1);

   unfold(chkHP, *unfTrue2Reco, iniHP);
   TMatrixD *resMreco=chkHP.histoAsM();
   testMaxDiff("\n true2reco(histoPair|load) ",finM,*resMreco,1);
   
   delete resMtrue;
   delete resMreco;
   delete unfTrue2Reco;
   delete unfReco2True;

   std::cout << "\n\n" << dashline << dashline << "\n";
 }

 if (1) { // check HistoPair2D_t::unfold
   HistoPair2D_t iniHP("iniHP"), finHP("finHP"), chkHP("chkHP");
   std::cout << "\n\n" << dashline << "\tCheck histoPair\n" << dashline << "\n";

   TMatrixD chkM(iniM);
   chkM.Zero();
   iniHP.assign(iniM,chkM,chkM);
   finHP.assign(finM,chkM,chkM);

   TMatrixD *unfTrue2Reco=UnfoldingMatrix_t::loadUnfM(M->getName(),dir,DYTools::analysisTag,0,1);
   TMatrixD *unfReco2True=UnfoldingMatrix_t::loadUnfM(M->getName(),dir,DYTools::analysisTag,1,1);

   chkHP.unfold(*unfReco2True, finHP);
   TMatrixD *resMtrue=chkHP.histoAsM();
   testMaxDiff("\n reco2true(histoPair|load) ",iniM,*resMtrue,1);

   chkHP.unfold(*unfTrue2Reco, iniHP);
   TMatrixD *resMreco=chkHP.histoAsM();
   testMaxDiff("\n true2reco(histoPair|load) ",finM,*resMreco,1);
   
   delete resMtrue;
   delete resMreco;
   delete unfTrue2Reco;
   delete unfReco2True;

   std::cout << "\n\n" << dashline << dashline << "\n";
 }

 delete M;
}
