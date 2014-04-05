#include "../Include/DYTools.hh"
#include "../Include/UnfoldingMatrix.h"

void printUnfMatrix(int analysisIs2D, int printCase=-1, int dataIdx=1) {

  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return;
  }

  TString dir1="../root_files7TeV/constants/DY_20130801_test/";
  TString dir2="/home/andriusj/cms/CMSSW_3_8_4/src/DrellYanDMDY-20130131/root_files/constants_debug/DY_m10+pr+a05+o03+pr_4839pb/";
  TString ending="_withFEWZ1D.root";

  if (1) { // full run
    //dir2="/home/andriusj/cms/CMSSW_3_8_4/src/DrellYanDMDY-20130324-closure20121217/root_files/constants/DY_m10+pr+a05+o03+pr_4839pb/";
    dir2=" /media/spektras/DrellYanDMDY-20121217/root_files/constants/DY_m10+pr+a05+o03+pr_4839pb/";
    ending="1D_PU.root";
  }
  if (1) {
    dir2="/home/andriusj/cms/CMSSW_3_8_4/src/DrellYanDMDY-20130131/root_files/constants_tmp/DY_m10+pr+a05+o03+pr_4839pb/";
    ending="_withFEWZ1D.root";
  }

  TString fname="detResponseExact_unfolding_constants1D.root";
  //fname="detResponse_unfolding_constants1D.root";
  //fname="efficiency.root";
  //fname="fsrDETexact_unfolding_constants1D.root";
  //fname="fsrDET_good_unfolding_constants1D.root";
  //fname="fsrDET_unfolding_constants1D.root";
  //fname="fsrExact_unfolding_constants1D.root";
  //fname="fsrGood_unfolding_constants1D.root";

  UnfoldingMatrix_t *M1=NULL, *M2=NULL;

  if (0) {
    M1 = new UnfoldingMatrix_t(UnfoldingMatrix::_cDET_Response,"detResponse");
    M1->autoLoadFromFile(dir1,DYTools::analysisTag);

    M2 = new UnfoldingMatrix_t(UnfoldingMatrix::_cDET_Response,"detResponse");
    M2->loadFromFile(dir2 + TString("detResponse_unfolding_constants") + ending,"");
  }
  else if (0) {
    M1 = new UnfoldingMatrix_t(UnfoldingMatrix::_cDET_Response,"detResponseExact");
    M1->autoLoadFromFile(dir1,DYTools::analysisTag);

    M2 = new UnfoldingMatrix_t(UnfoldingMatrix::_cDET_Response,"detResponseExact");
    M2->loadFromFile(dir2 + TString("detResponseExact_unfolding_constants") + ending,"");
  }
  else if (0) {
    M1 = new UnfoldingMatrix_t(UnfoldingMatrix::_cFSR,"fsrExact");
    M1->autoLoadFromFile(dir1,DYTools::analysisTag);

    M2 = new UnfoldingMatrix_t(UnfoldingMatrix::_cFSR,"fsrExact");
    M2->loadFromFile(dir2 + TString("fsrExact_unfolding_constants") + ending,"");
  }
  else if (0) {
    M1 = new UnfoldingMatrix_t(UnfoldingMatrix::_cFSR,"fsrGood");
    M1->autoLoadFromFile(dir1,DYTools::analysisTag);

    M2 = new UnfoldingMatrix_t(UnfoldingMatrix::_cFSR,"fsrGood");
    M2->loadFromFile(dir2 + TString("fsrGood_unfolding_constants") + ending,"");
  }
  else if (0) {
    M1 = new UnfoldingMatrix_t(UnfoldingMatrix::_cFSR_DET,"fsrDET");
    M1->autoLoadFromFile(dir1,DYTools::analysisTag);

    M2 = new UnfoldingMatrix_t(UnfoldingMatrix::_cFSR_DET,"fsrDET");
    M2->loadFromFile(dir2 + TString("fsrDET_unfolding_constants") + ending,"");
  }
  else if (0) {
    M1 = new UnfoldingMatrix_t(UnfoldingMatrix::_cFSR_DET,"fsrDETexact");
    M1->autoLoadFromFile(dir1,DYTools::analysisTag);

    M2 = new UnfoldingMatrix_t(UnfoldingMatrix::_cFSR_DET,"fsrDETexact");
    M2->loadFromFile(dir2 + TString("fsrDETexact_unfolding_constants") + ending,"");
  }
  else if (1) {
    M1 = new UnfoldingMatrix_t(UnfoldingMatrix::_cFSR_DET,"fsrDET_good");
    M1->autoLoadFromFile(dir1,DYTools::analysisTag);

    M2 = new UnfoldingMatrix_t(UnfoldingMatrix::_cFSR_DET,"fsrDET_good");
    M2->loadFromFile(dir2 + TString("fsrDET_good_unfolding_constants") + ending,"");
  }

  int print1=((printCase==-1) || (printCase==1)) ? 1:0;
  int print2=((printCase==-1) || (printCase==2)) ? 1:0;

  if (dataIdx==1) {
    std::cout << "\n\nYields\n\n";
    if (print1) M1->printYields();
    if (print2) M2->printYields();
  }

  if (dataIdx==2) {
    std::cout << "\n\nMigration\n\n";
    if (print1) M1->printMigration();
    if (print2) M2->printMigration();
  }

  if (dataIdx==3) {
    std::cout << "\n\nResponse (error is symmetric!)\n\n";
    if (print1) M1->printResponse();
    if (print2) M2->printResponse();
  }
 
 if (dataIdx==4) {
   std::cout << "\n\nInvResponse\n\n";
    if (print1) M1->printInvResponse();
    if (print2) M2->printInvResponse();
  }

 if (M1) delete M1;
 if (M2) delete M2;
}
