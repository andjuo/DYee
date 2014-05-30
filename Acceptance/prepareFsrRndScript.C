#include <TROOT.h>
#include <TRandom.h>
#include <TString.h>
#include <iostream>
#include <fstream>


int prepareFsrRndScript(int seed, int nExps, TString outName,
			TString macroName, TString shortName) {

  TString systModeFormat="FSR_RND_STUDYid%dfixed%6.5lf";

  std::ofstream fout(outName.Data());
  fout << "#!/bin/bash\n";
  fout << "\n";
  fout << "analysisIs2D=$1\n";
  fout << "config=$2\n";
  fout << "runMode=$3\n";
  fout << "\n";
  fout << "if [ ${#runMode} -eq 0 ] ; then\n";
  fout << "  echo \"this_script  analysisIs2D  config  runMode\"\n";
  fout << "  exit\n";
  fout << "fi\n";
  fout << "\n";

  gRandom->SetSeed(seed);

  for (int i=0; i<nExps; i++) {
    int id=i+1000;
    fout << "root -l -q -b " << macroName
	 <<"+\\(${analysisIs2D},\\\"${config}\\\",";
    fout << "${runMode},";
    fout << "DYTools::FSR_RND_STUDY,";
    fout << "\\\"" << Form(systModeFormat.Data(),id,gRandom->Gaus(0,1.));
    fout << "\\\"\\) > log-" << shortName << "-$((${analysisIs2D}+1))D-"
	 << "${runModeStr}-${rndStudyStr}" << id << ".out\n";
  }
  fout.close();
  return 1;
}
