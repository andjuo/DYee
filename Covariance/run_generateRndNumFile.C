#include "generateRndNumFile.C"

int run_generateRndNumFile(int iSeedMin, int iSeedMax,
			   int nNumbers=11940000,
	TString outFNameBase="../root_files_reg/theory/rndSequence_")
{
  int ok=1;
  for (int iSeed=iSeedMin; (ok==1) && (iSeed<=iSeedMax); ++iSeed) {
    TString fname=outFNameBase + TString(Form("%d",iSeed));
    ok=generateRndNumFile(iSeed,nNumbers,fname);
  }
  return ok;
}
