#include "../BLUE/resultCombiner.C"

void runBLUE_selfTest1D(int normalized=0) {
  //TString xs_fname="xsecForBlue-1D-yieldStatOnly-fromCS.dat";
  //TString cov_fname="covForBlue-1D-yieldStatOnly.dat";
  TString xs_fname="dir-forBlue/xsecForBlue-1D-fromCov.dat";
  TString cov_fname="dir-forBlue/covForBlue-1D.dat";

  TString tag;
  switch(normalized) {
  case 0: tag="1D-absolute"; break;
  case 1: tag="1D-normalized"; break;
  case 2: tag="1D-rshape"; break;
  default:
    std::cout << "wrong value for normalized=" << normalized << "\n";
    return;
  }
  xs_fname.ReplaceAll("1D",tag);
  cov_fname.ReplaceAll("1D",tag);

  resultCombiner(xs_fname,cov_fname, "41",
		 xs_fname,cov_fname, "41");
  return;
}
