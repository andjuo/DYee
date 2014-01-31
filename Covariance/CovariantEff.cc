#include "CovariantEff.h"

//const TString tnpTag="et6_eta5";
//const TString tnpTag="et6spec_eta5";

//const TString tnpMCFileName=TString("../config_files/sf_mc_${tnpTag}_vilnius.conf").ReplaceAll("${tnpTag}",tnpTag);
//const TString tnpDataFileName=TString("../config_files/sf_data_${tnpTag}_vilnius.conf").ReplaceAll("${tnpTag}",tnpTag);
//const TString filename_mc="../config_files/fall11mc_vilnius.input";

const TString triggerSetString="Full2011_hltEffOld";

//EffArray_t *ro_Data_loc=NULL;
//EffArray_t *ro_MC_loc=NULL;

// ----------------------------------------------------------------
// ----------------------------------------------------------------

CovariantEffMgr_t::CovariantEffMgr_t() :  
  BaseClass_t("CovariantEffMgr_t"), 
  FInpMgr(),
  FEvtSelector(),
  FPUStr(),
  FTriggers(TrigSet_UNDEFINED,false,0),
  FInitOk(0)
{
}


// ----------------------------------------------------------------

int CovariantEffMgr_t::Setup(const TString &confFileName, int nExps) {

  int res=1;
  int createDestDir=0;
  if (!initManagers(confFileName,DYTools::NORMAL_RUN,FInpMgr,FEvtSelector,FPUStr,createDestDir)) {
    res=reportError("Setup","failed to initialize managers");
    return res;
  }

  //ro_Data_loc=new EffArray_t[nExps];
  //ro_MC_loc= new EffArray_t[nExps];
  //assert(ro_Data_loc && ro_MC_loc);

  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;

  // Read efficiency constants from ROOT files
  // This has to be done AFTER configuration file is parsed
  if (!fillEfficiencyConstants( FInpMgr, systMode ) ) {
    return reportError("Setup");
  }

  int debug_pseudo_exps=0;
  if (!preparePseudoExps(nExps,debug_pseudo_exps)) 
    return reportError("Setup");

  FInitOk=res;

  return 1;
}

// ----------------------------------------------------------------



// ----------------------------------------------------------------
// ----------------------------------------------------------------

StudyArr_t::StudyArr_t(const TString &baseName, int nMassBins, int nExps,
		       int binCount, double minVal, double maxVal) :
  BaseClass_t(Form("StudyArr_t(%s)",baseName.Data())),
  FVec(),
  FBaseName(baseName),
  FReady(0) {
  this->Init(nMassBins,nExps, binCount,minVal,maxVal);
}


// ----------------------------------------------------------------

int StudyArr_t::Init(int nMassBins, int nExps, 
		     int binCount, double minVal, double maxVal) {
  FReady=0;
  if (binCount==0) return 0;
  if (FVec.size()) { BaseClass_t::printWarning("::Init will clear FVec"); }
  FVec.clear();
  FVec.reserve(nMassBins);
  for (int i=0; i<nMassBins; ++i) {
    std::vector<TH1D*> *vec=new std::vector<TH1D*>();
    FVec.push_back(vec);
    vec->reserve(nExps);
    for (int j=0; j<nExps; ++j) {
      TString name=Form("%s_%d_exp%d",FBaseName.Data(),i,j);
      TH1D* h=new TH1D(name,name, binCount,minVal,maxVal);
      h->SetDirectory(0);
      h->Sumw2();
      vec->push_back(h);
    }
  }
  FReady=1;
  return 1;
}

// ----------------------------------------------------------------
// ----------------------------------------------------------------


// ----------------------------------------------------------------
// ----------------------------------------------------------------

