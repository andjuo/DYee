#include "subtractBackground.C"

int run_subtractBackground(int analysisIs2D,
			   int the_case, int systModeFlag=0, int debug=-1) {
  TString confName="default";
  DYTools::TRunMode_t runMode=DebugInt2RunMode(debug);

  //DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;
  DYTools::TSystematicsStudy_t systMode=DYTools::ESCALE_DIFF_0000;

  switch(systModeFlag) {
  case 0: ; break;
  case 1: systMode=DYTools::UNREGRESSED_ENERGY; break;
  }

  switch(the_case) {
  case 0: ; break;
  case 1: 
    confName="defaultAdHoc";
    systMode=DYTools::APPLY_ESCALE;
    break;
  case 2:
  case 3:
    confName="defaultAdHoc";
    systMode=DYTools::ESCALE_STUDY_RND;
    break;
  default:
    std::cout << "the_case=" << the_case << " is not ready\n";
    return retCodeError;
  }

  // Setup the mass-rapidity setting
  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return retCodeError;
  }

  analysisIs2D=-111; // further calls to DYTools::setup do nothing
  int res=1;

  if (systMode!=DYTools::ESCALE_STUDY_RND) {
    res=subtractBackground(analysisIs2D,confName,runMode,systMode);
  }
  else {
    gBenchmark->Start("run_subtractBackground");
    int iSeedMin=-1;
    int iSeedMax=-1;
    int dSeed=1;
    if (the_case!=3) {
      InputFileMgr_t inpMgr;
      if (!inpMgr.Load(confName)) return retCodeError;
      iSeedMin=inpMgr.userKeyValueAsInt("SEEDMIN");
      iSeedMax=inpMgr.userKeyValueAsInt("SEEDMAX");
    }
    else {
      iSeedMin=-111;
      iSeedMax= 111;
      dSeed=iSeedMax-iSeedMin;
    }
    if (res) res=retCodeOk;
    for (int iSeed=iSeedMin; (res==retCodeOk) && (iSeed<=iSeedMax);
	 iSeed+=dSeed) {
      std::cout << "\n\n\tstarting iSeed=" << iSeed << "\n\n";
      res=subtractBackground(analysisIs2D,confName,runMode,systMode,iSeed);
    }
    ShowBenchmarkTime("run_subtractBackground");
  }
  return res;
}
