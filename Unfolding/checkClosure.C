#include "../Include/DYTools.hh"
#include "../Include/UnfoldingMatrix.h"


void checkClosure() {
 TString dir="../root_files7TeV/constants/DY_20130801_test_DebugRun/";
 UnfoldingMatrix_t *M=NULL;

 UnfoldingMatrix::TUnfoldingMatrixType_t kind=UnfoldingMatrix::_cFSR;
 //M=new UnfoldingMatrix_t(kind, "fsrExact");
 //M=new UnfoldingMatrix_t(kind, "fsrGood");
 kind=UnfoldingMatrix::_cDET_Response;
 M=new UnfoldingMatrix_t(kind, "detResponseExact");
 //M=new UnfoldingMatrix_t(kind, "detResponse");

 M->autoLoadFromFile(dir,DYTools::analysisTag);
 
 TVectorD ini= * (M->getIniVec());
 TVectorD fin= * (M->getFinVec());

 TVectorD chk(fin);
 //TVectorD diff(fin);
 chk.Zero();

 std::cout << "ini="; ini.Print();

 //unfold_true2reco(chk,*M,fin);
 //std::cout << "true2reco  "; chk.Print();

 unfold_reco2true(chk,*M,fin);
 std::cout << "reco2true "; chk.Print();
 testMaxDiff("\n reco2true ",ini,chk);
 testMaxRelDiff("\n reco2true ",ini,chk);
 //diff=ini-chk;
 //std::cout << "diff "; diff.Print();

 std::cout << dashline;

 std::cout << "fin="; fin.Print();
 chk.Zero();
 unfold_true2reco(chk,*M,ini);
 std::cout << "true2reco "; chk.Print();
 testMaxDiff("\n true2reco ",fin,chk);
 testMaxRelDiff("\n true2reco ",fin,chk);

 std::cout << "\n\n";
 delete M;
}
