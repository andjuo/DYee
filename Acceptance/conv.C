#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include <TFile.h>
#include <TMatrixD.h>
#include <iostream>
#include <fstream>

int conv() {
  std::ifstream fin("alexey_preFSR_acc.txt");
  int ibin;
  double a,ae;
  TMatrixD A(7,24);
  TMatrixD Aerr(7,24);
  A.Zero();
  Aerr.Zero();

  TH2D *hAcc=createBaseH2("hAcc","hAcc (Alexey)",1);

  int im=1;
  int iy=0;
  while (!fin.eof()) {
    fin >> ibin >> a >> ae;
    if (fin.eof()) break;
    A(im,iy)=a;
    Aerr(im,iy)=ae;
    hAcc->SetBinContent(im+1,iy+1,a);
    hAcc->SetBinError  (im+1,iy+1,ae);
    iy++;
    if (iy>=24) {
      im++;
      iy=0;
    }
  }
  fin.close();

  A.Print();

  TFile fout("Acc_preFSR.root","recreate");
  A.Write("Acc");
  Aerr.Write("AccErr");
  hAcc->Write("hAcc");
  fout.Close();
  return 1;
}
