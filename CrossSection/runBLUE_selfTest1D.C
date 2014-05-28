#include "../BLUE/resultCombiner.C"

void runBLUE_selfTest1D() {
  //TString xs_fname="xsecForBlue-1D-yieldStatOnly-fromCS.dat";
  //TString cov_fname="covForBlue-1D-yieldStatOnly.dat";
  TString xs_fname="xsecForBlue-1D-fromCov.dat";
  TString cov_fname="covForBlue-1D.dat";

  resultCombiner(xs_fname,cov_fname, "41",
		 xs_fname,cov_fname, "41");
  return;
}
