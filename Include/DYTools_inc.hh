#ifndef DYToolsUI_HH
#define DYToolsUI_HH

//#include "../Include/DYTools.hh"
//#include <stdarg.h>

// ------------------------------------------------------------------
//   Names and recognition
// ------------------------------------------------------------------

//inline
TString MassBinningName(DYTools::TMassBinning_t set) {
  using namespace DYTools;
  TString name;
  switch(set) {
  case _MassBins_Undefined: name="UNDEFINED"; break;
  case _MassBins_default: name="MassBinsDefault"; break;
  case _MassBins_2012: name="MassBins2012"; break;
  case _MassBins_2011: name="MassBins2011"; break;
  case _MassBins_test4: name="MassBinsTest4"; break;
  case _MassBins_Zpeak: name="MassBinsZpeak"; break;
  case _MassBins_noUnderflow: name="MassBins_noUnderflow"; break;
  case _MassBins_withFullOverflow: name="MassBins_withFullOverflow"; break;
  case _MassBins_bins100GeV: name="MassBins_bins100GeV"; break;
  case _MassBins_bins5GeV: name="MassBins_bins5GeV"; break;
  case _MassBins_finerMassRange: name="MassBins_finerMassRange"; break;
  default: name="UNKNOWN_MASS_BINNING";
  }
  return name;
}

// ------------------------------------------------------------------

//inline
TString CurrentMassBinningName() {
  return MassBinningName(DYTools::massBinningSet);
}

// ------------------------------------------------------------------
// ------------------------------------------------------------------

//inline
TString SystematicsStudyName(DYTools::TSystematicsStudy_t study) {
  using namespace DYTools;
  TString name;
  switch(study) {
  case SYST_MODE_FAILURE: name="Syst_Mode_Failure(error in the code)"; break;
  case NO_SYST: name="NormalRun_noSyst"; break;
  case RESOLUTION_STUDY: name="ResolutionStudy"; break;
  case FSR_STUDY: name="FSRStudy"; break;
  case FSR_RND_STUDY: name="FSRRndStudy"; break;
  case FSR_5plus: name="FSR_5plus"; break;
  case FSR_5minus: name="FSR_5minus"; break;
  case APPLY_ESCALE: name="ApplyEscale"; break;
  case ESCALE_RESIDUAL: name="EScale_residual"; break;
  case ESCALE_STUDY: name="EScale_study"; break;
  case ESCALE_STUDY_RND: name="EScale_study_randomized"; break;
  case ESCALE_DIFF_0000: name="EScaleDiff0000"; break;
  case ESCALE_DIFF_0005: name="EScaleDiff0005"; break;
  case ESCALE_DIFF_0010: name="EScaleDiff0010"; break;
  case ESCALE_DIFF_0015: name="EScaleDiff0015"; break;
  case ESCALE_DIFF_0020: name="EScaleDiff0020"; break;
  case UNREGRESSED_ENERGY: name="Unregressed_energy"; break;
  case LOWER_ET_CUT: name="LowerEtCut"; break;
  case NO_REWEIGHT: name="NoReweight"; break;
  case NO_REWEIGHT_PU: name="NoPUReweight"; break;
  case NO_REWEIGHT_FEWZ: name="NoFewzReweight"; break;
  case SYST_RND: name="Syst_random"; break;
  case TAG_ID: name="TagID"; break;
  case TAG_PT: name="TagPt"; break;
  case PU_STUDY: name="PUStudy"; break;
  case PU_RND_STUDY: name="PURndStudy"; break;
  case PILEUP_5plus: name="Pileup5plus"; break;
  case PILEUP_5minus: name="Pileup5minus"; break;
  case UNREG_PU5plus: name="UnregEnPU5plus"; break;
  case UNREG_PU5minus: name="UnregEnPU5minus"; break;
  case UNREG_FSR5plus: name="UnregEnFSR5plus"; break;
  case UNREG_FSR5minus: name="UnregEnFSR5minus"; break;
  case UNREG_TagID: name="UnregTagID"; break;
  case UNREG_TagPt: name="UnregTagPt"; break;
  default: name="UNKNOWN_SYSTEMATICS_NAME";
  }
  return name;
}

// ------------------------------------------------------------------

//inline
TString generateSystTag(DYTools::TSystematicsStudy_t systMode) {
  TString tag;
  if (systMode!=DYTools::NO_SYST) {
    tag=TString("_") + SystematicsStudyName(systMode);
  }
  return tag;
}

// ------------------------------------------------------------------

//inline
TString RunModeName(DYTools::TRunMode_t study) {
  using namespace DYTools;
  TString name;
  switch(study) {
  case NORMAL_RUN: name="NormalRun"; break;
  case LOAD_DATA: name="LoadData"; break;
  case DEBUG_RUN: name="DebugRun"; break;
  case DEBUG_LOAD: name="DebugLoad"; break;
  default: name="UNKNOWN_RunMode_NAME";
  }
  return name;
}

// ------------------------------------------------------------------

//inline
DYTools::TRunMode_t DebugInt2RunMode(int debug) {
  using namespace DYTools;
  TRunMode_t runMode=NO_RUN;
  if (debug==1) runMode=DEBUG_RUN;
  else if (debug==0) runMode=NORMAL_RUN;
  else if ((debug==-1) || (debug==-2)) runMode=LOAD_DATA;
  else if (debug== 2) runMode=DEBUG_LOAD;
  return runMode;
}

// ------------------------------------------------------------------

//inline
TString generateRunTag(DYTools::TRunMode_t runMode) {
  TString tag= (runMode==DYTools::NORMAL_RUN) ? "" : RunModeName(runMode);
  return tag;
}

// ------------------------------------------------------------------

//inline
TString TnPMethodName(DYTools::TTnPMethod_t method) {
  using namespace DYTools;
  TString name;
  switch(method) {
  case COUNTnCOUNT: name="Count&Count"; break;
  case COUNTnFIT: name="Count&Fit"; break;
  case FITnFIT: name="Fit&Fit"; break;
  default: name="Unknown_TnP_method";
  }
  return name;
}

// ------------------------------------------------------------------

//inline
TString EfficiencyKindName(DYTools::TEfficiencyKind_t kind) {
  using namespace DYTools;
  TString name;
  switch(kind) {
  case RECO: name="RECO"; break;
  case ID: name="ID"; break;
  case HLT: name="HLT"; break;
  case HLT_leg1: name="HLTleg1"; break;
  case HLT_leg2: name="HLTleg2"; break;
  case effNONE: name="effNONE"; break;
  default: name="Unknown_Efficiency_Kind";
  }
  return name;
}

// ------------------------------------------------------------------

//inline
TString EtBinSetName(DYTools::TEtBinSet_t set) {
  using namespace DYTools;
  TString name;
  switch(set) {
  case ETBINS_UNDEFINED: name="EtBinsUndefined"; break;
  case ETBINS1: name="EtBins1"; break;
  case ETBINS5: name="EtBins5"; break;
  case ETBINS6: name="EtBins6"; break;
  case ETBINS6x: name="EtBins6x"; break;
  case ETBINS6short: name="EtBins6short"; break;
  case ETBINS6syst: name="EtBins6syst"; break;
  case ETBINS6alt: name="EtBins6alt"; break;
  case ETBINS6altB: name="EtBins6altB"; break;
  case ETBINS7: name="EtBins7"; break;
  case ETBINS7alt: name="EtBins7alt"; break;
  case ETBINS7altB: name="EtBins7altB"; break;
  case ETBINS8: name="EtBins8"; break;
  case ETBINS9: name="EtBins9"; break;
  default: name="Unknown_Et_BinSet";
    std::cout << "EtBinSetName=<" << name << ">\n";
    assert(0);
  }
  return name;
}

// ------------------------------------------------------------------

//inline
TString EtaBinSetName(DYTools::TEtaBinSet_t set) {
  using namespace DYTools;
  TString name;
  switch(set) {
  case ETABINS_UNDEFINED: name="EtaBinsUndefined"; break;
  case ETABINS1: name="EtaBins1"; break;
  case ETABINS2: name="EtaBins2"; break;
  case ETABINS2Negs: name="EtaBins2Negs"; break;
  case ETABINS3: name="EtaBins3"; break;
  case ETABINS3Negs: name="EtaBins3Negs"; break;
  case ETABINS5: name="EtaBins5"; break;
  case ETABINS5nomerge: name="EtaBins5nomerge"; break;
  case ETABINS5egamma: name="EtaBins5egamma"; break;
  case ETABINS5_max25: name="EtaBins5_max25"; break;
  case ETABINS5Negs: name="EtaBins5Negs"; break;
  case ETABINS4test: name="EtaBins4test"; break;
  case ETABINS4testNegs: name="EtaBins4testNegs"; break;
  case ETABINS4alt: name="EtaBins4alt"; break;
  case ETABINS4altNegs: name="EtaBins4altNegs"; break;
  case ETABINS5alt: name="EtaBins5alt"; break;
  case ETABINS5altNegs: name="EtaBins5altNegs"; break;
  case ETABINS7: name="EtaBins7"; break;
  case ETABINS8alt: name="EtaBins8alt"; break;
  case ETABINS8altNegs: name="EtaBins8altNegs"; break;
  case ETABINS9: name="EtaBins9"; break;
  case ETABINS14: name="EtaBins14"; break;
  default: name="Unknown_Eta_BinSet";
    std::cout << "EtaBinSetName=<" << name << ">\n";
    assert(0);
  }
  return name;
}

// ------------------------------------------------------------------

//inline
TString DataKindName(DYTools::TDataKind_t kind) {
  using namespace DYTools;
  TString name;
  switch(kind) {
  case DATA: name="Data"; break;
  case MC: name="MC"; break;
  default: name="Unknown_DataKind";
  }
  return name;
}

// ------------------------------------------------------------------
// ------------------------------------------------------------------

//inline
DYTools::TSystematicsStudy_t DetermineSystematicsStudy(const TString &str) {
  using namespace DYTools;
  DYTools::TSystematicsStudy_t study=NO_SYST;
  if (str.Contains("NORMALRUN_NOSYST") || str.Contains("NormalRun_noSyst") ||
      str.Contains("NO_SYST") || str.Contains("NoSyst")) { study=NO_SYST; }
  else if (str.Contains("RESOLUTION_STUDY") || str.Contains("ResolutionStudy")) study=RESOLUTION_STUDY;
  else if (str.Contains("FSR_STUDY") || str.Contains("FSRStudy")) study=FSR_STUDY;
  else if (str.Contains("FSR_RND_STUDY") || str.Contains("FSRRndStudy")) study=FSR_RND_STUDY;
  else if (str.Contains("UNREG_FSR5PLUS") || str.Contains("UnregEnFSR5plus")) study=UNREG_FSR5plus;
  else if (str.Contains("UNREG_FSR5MINUS") || str.Contains("UnregEnFSR5minus")) study=UNREG_FSR5minus;
  else if (str.Contains("FSR_5PLUS") || str.Contains("FSR_5plus")) study=FSR_5plus;
  else if (str.Contains("FSR_5MINUS") || str.Contains("FSR_5minus")) study=FSR_5minus;
  else if (str.Contains("APPLY_ESCALE") || str.Contains("ApplyEscale")) study=APPLY_ESCALE;
  else if (str.Contains("ESCALE_RESIDUAL") || str.Contains("EScale_residual")) study=ESCALE_RESIDUAL;
  else if (str.Contains("ESCALE_STUDY_RND") || str.Contains("EScale_study_randomized")) study=ESCALE_STUDY_RND;
  else if (str.Contains("ESCALE_STUDY") || str.Contains("EScale_study")) study=ESCALE_STUDY;
  else if (str.Contains("ESCALE_DIFF_0000") || str.Contains("EScale0000")) study=ESCALE_DIFF_0000;
  else if (str.Contains("ESCALE_DIFF_0005") || str.Contains("EScale0005")) study=ESCALE_DIFF_0005;
  else if (str.Contains("ESCALE_DIFF_0010") || str.Contains("EScale0010")) study=ESCALE_DIFF_0010;
  else if (str.Contains("ESCALE_DIFF_0015") || str.Contains("EScale0015")) study=ESCALE_DIFF_0015;
  else if (str.Contains("ESCALE_DIFF_0020") || str.Contains("EScale0020")) study=ESCALE_DIFF_0020;
  else if (str.Contains("UNREGRESSED_ENERGY") || str.Contains("Unregressed_energy")) study=UNREGRESSED_ENERGY;
  else if (str.Contains("LOWER_ET_CUT") || str.Contains("LowerEtCut")) study=LOWER_ET_CUT;
  else if (str.Contains("NO_REWEIGHT_PU") || str.Contains("NoPUReweight")) study=NO_REWEIGHT_PU;
  else if (str.Contains("NO_REWEIGHT_FEWZ") || str.Contains("NoFewzReweight")) study=NO_REWEIGHT_FEWZ;
  else if (str.Contains("NO_REWEIGHT") || str.Contains("NoReweight")) study=NO_REWEIGHT;
  else if (str.Contains("SYST_RND") || str.Contains("SYST_RANDOM") || str.Contains("Syst_random")) study=SYST_RND;
  else if (str.Contains("PU_STUDY") || str.Contains("PUStudy")) study=PU_STUDY;
  else if (str.Contains("PU_RND_STUDY") || str.Contains("PURndStudy")) study=PU_RND_STUDY;
  else if (str.Contains("PILEUP_5plus") || str.Contains("PILEUP_5PLUS")) study=PILEUP_5plus;
  else if (str.Contains("PILEUP_5minus") || str.Contains("PILEUP_5MINUS")) study=PILEUP_5minus;
  else if (str.Contains("UNREG_PU5plus") || str.Contains("UnregEnPU5plus")) study=UNREG_PU5plus;
  else if (str.Contains("UNREG_PU5minus") || str.Contains("UnregEnPU5minus")) study=UNREG_PU5minus;
  else if (str.Contains("UNREG_TagID") || str.Contains("UnregTagID")) study=UNREG_TagID;
  else if (str.Contains("UNREG_TagPt") || str.Contains("UnregTagPt")) study=UNREG_TagPt;
  else if (str.Contains("TAG_ID") || str.Contains("TagID")) study=TAG_ID;
  else if (str.Contains("TAG_PT") || str.Contains("TagPt")) study=TAG_PT;
  else {
    std::cout << "DetermineSystematicsStudy: failure at <" << str << ">\n";
    assert(0);
  }
  return study;
}

// ------------------------------------------------------------------

int usesUnregEnergy(DYTools::TSystematicsStudy_t systMode) {
  using namespace DYTools;
  int yes=-1;
  switch(systMode) {
  case NO_SYST:
  case RESOLUTION_STUDY:
  case FSR_STUDY:
  case FSR_RND_STUDY:
  case FSR_5plus:
  case FSR_5minus:
  case ESCALE_RESIDUAL:
  case LOWER_ET_CUT:
  case NO_REWEIGHT:
  case NO_REWEIGHT_PU:
  case NO_REWEIGHT_FEWZ:
  case SYST_RND:
  case TAG_ID:
  case TAG_PT:
  case PU_STUDY:
  case PU_RND_STUDY:
  case PILEUP_5plus:
  case PILEUP_5minus:
    yes=0;
    break;
  case UNREGRESSED_ENERGY:
  case UNREG_PU5plus:
  case UNREG_PU5minus:
  case UNREG_FSR5plus:
  case UNREG_FSR5minus:
  case UNREG_TagID:
  case UNREG_TagPt:
    yes=1;
    break;
  case ESCALE_STUDY:
  case ESCALE_STUDY_RND:
  case APPLY_ESCALE:
  case ESCALE_DIFF_0000:
  case ESCALE_DIFF_0005:
  case ESCALE_DIFF_0010:
  case ESCALE_DIFF_0015:
  case ESCALE_DIFF_0020:
    yes=2;
    break;
  default:
    std::cout << "usesUnregEn: unhandled case\n";
    assert(0);
  }
  return yes;
}

// ------------------------------------------------------------------

int isGlobalRndStudy(DYTools::TSystematicsStudy_t systMode) {
  using namespace DYTools;
  int yes=0;
  switch(systMode) {
  case NO_SYST:
  case RESOLUTION_STUDY:
  case FSR_STUDY:
  case FSR_5plus:
  case FSR_5minus:
  case ESCALE_RESIDUAL:
  case LOWER_ET_CUT:
  case NO_REWEIGHT:
  case NO_REWEIGHT_PU:
  case NO_REWEIGHT_FEWZ:
  case SYST_RND:
  case TAG_ID:
  case TAG_PT:
  case PU_STUDY:
  case PILEUP_5plus:
  case PILEUP_5minus:
  case UNREGRESSED_ENERGY:
  case UNREG_PU5plus:
  case UNREG_PU5minus:
  case UNREG_FSR5plus:
  case UNREG_FSR5minus:
  case UNREG_TagID:
  case UNREG_TagPt:
  case ESCALE_STUDY:
  case ESCALE_STUDY_RND:
  case APPLY_ESCALE:
  case ESCALE_DIFF_0000:
  case ESCALE_DIFF_0005:
  case ESCALE_DIFF_0010:
  case ESCALE_DIFF_0015:
  case ESCALE_DIFF_0020:
    yes=0;
    break;
  case FSR_RND_STUDY:
  case PU_RND_STUDY:
    yes=1;
    break;
  default:
    std::cout << "isGlobalRndStudy: unhandled case\n";
    assert(0);
  }
  return yes;
}

// ------------------------------------------------------------------

//inline
DYTools::TTnPMethod_t DetermineTnPMethod(const TString &str) {
  using namespace DYTools;
  DYTools::TTnPMethod_t method=COUNTnCOUNT;
  if (str.Contains("COUNTnCOUNT") || str.Contains("Count&Count")) method=COUNTnCOUNT;
  else if (str.Contains("COUNTnFIT") || str.Contains("Count&Fit")) method=COUNTnFIT;
  else if (str.Contains("FITnFIT") || str.Contains("Fit&Fit")) method=FITnFIT;
  else {
    std::cout << "DetermineTnPMethod failed at <" << str << ">\n";
    assert(0);
  }
  return method;
}

// ------------------------------------------------------------------

//inline
DYTools::TEfficiencyKind_t DetermineEfficiencyKind(const TString &str) {
  using namespace DYTools;
  DYTools::TEfficiencyKind_t kind=RECO;
  if (str.Contains("GSF") || str.Contains("Reco") || str.Contains("RECO")) kind=RECO;
  else if (str.Contains("ID")) kind=ID;
  else if (str.Contains("HLTleg1")) kind=HLT_leg1;
  else if (str.Contains("HLT_leg1")) kind=HLT_leg1;
  else if (str.Contains("HLTleg2")) kind=HLT_leg2;
  else if (str.Contains("HLT_leg2")) kind=HLT_leg2;
  else if (str.Contains("HLT")) kind=HLT;
  else if (str.Contains("effNONE")) kind=effNONE;
  else {
    std::cout << "DetermineEfficiencyKind failed at <" << str << ">\n";
    assert(0);
  }
  return kind;
}

// ------------------------------------------------------------------

//inline
DYTools::TEtBinSet_t DetermineEtBinSet(const TString& str) {
  using namespace DYTools;
  DYTools::TEtBinSet_t kind=ETBINS1;
  if (str.Contains("ETBINS1") || str.Contains("EtBins1")) kind=ETBINS1;
  else if (str.Contains("ETBINS5") || str.Contains("EtBins5")) kind=ETBINS5;
  else if (str.Contains("ETBINS6altB") || str.Contains("EtBins6altB")) kind=ETBINS6altB;
  else if (str.Contains("ETBINS6alt") || str.Contains("EtBins6alt")) kind=ETBINS6alt;
  else if (str.Contains("ETBINS6syst") || str.Contains("EtBins6syst")) kind=ETBINS6syst;
  else if (str.Contains("ETBINS6short") || str.Contains("EtBins6short")) kind=ETBINS6short;
  else if (str.Contains("ETBINS6x") || str.Contains("EtBins6x")) kind=ETBINS6x;
  else if (str.Contains("ETBINS6") || str.Contains("EtBins6")) kind=ETBINS6;
  else if (str.Contains("ETBINS7altB") || str.Contains("EtBins7altB")) kind=ETBINS7altB;
  else if (str.Contains("ETBINS7alt") || str.Contains("EtBins7alt")) kind=ETBINS7alt;
  else if (str.Contains("ETBINS7") || str.Contains("EtBins7")) kind=ETBINS7;
  else if (str.Contains("ETBINS8") || str.Contains("EtBins8")) kind=ETBINS8;
  else if (str.Contains("ETBINS9") || str.Contains("EtBins9")) kind=ETBINS9;
  else {
    std::cout << "DetermineEtBinSet failed at <" << str  << ">\n";
    assert(0);
  }
  return kind;
}

// ------------------------------------------------------------------

//inline
DYTools::TEtaBinSet_t DetermineEtaBinSet(const TString& str) {
  using namespace DYTools;
  DYTools::TEtaBinSet_t kind=ETABINS1;
  if (str.Contains("ETABINS14") || str.Contains("EtaBins14")) kind=ETABINS14;
  else if (str.Contains("ETABINS1") || str.Contains("EtaBins1")) kind=ETABINS1;
  else if (str.Contains("ETABINS2Negs") || str.Contains("EtaBins2Negs")) kind=ETABINS2Negs;
  else if (str.Contains("ETABINS2") || str.Contains("EtaBins2")) kind=ETABINS2;
  else if (str.Contains("ETABINS3Negs") || str.Contains("EtaBins3Negs")) kind=ETABINS3Negs;
  else if (str.Contains("ETABINS3") || str.Contains("EtaBins3")) kind=ETABINS3;
  else if (str.Contains("ETABINS5altNegs") || str.Contains("EtaBins5altNegs")) kind=ETABINS5altNegs;
  else if (str.Contains("ETABINS5alt") || str.Contains("EtaBins5alt")) kind=ETABINS5alt;
  else if (str.Contains("ETABINS5_max25") || str.Contains("EtaBins5_max25")) kind=ETABINS5_max25;
  else if (str.Contains("ETABINS5egamma") || str.Contains("EtaBins5egamma")) kind=ETABINS5egamma;
  else if (str.Contains("ETABINS5nomerge") || str.Contains("EtaBins5nomerge")) kind=ETABINS5nomerge;
  else if (str.Contains("ETABINS5") || str.Contains("EtaBins5")) kind=ETABINS5;
  else if (str.Contains("ETABINS4altNegs") || str.Contains("EtaBins4altNegs")) kind=ETABINS4altNegs;
  else if (str.Contains("ETABINS4alt") || str.Contains("EtaBins4alt")) kind=ETABINS4alt;
  else if (str.Contains("ETABINS4testNegs") || str.Contains("EtaBins4testNegs")) kind=ETABINS4testNegs;
  else if (str.Contains("ETABINS4test") || str.Contains("EtaBins4test")) kind=ETABINS4test;
  else if (str.Contains("ETABINS7") || str.Contains("EtaBins7")) kind=ETABINS7;
  else if (str.Contains("ETABINS8altNegs") || str.Contains("EtaBins8altNegs")) kind=ETABINS8altNegs;
  else if (str.Contains("ETABINS8alt") || str.Contains("EtaBins8alt")) kind=ETABINS8alt;
  else if (str.Contains("ETABINS9") || str.Contains("EtaBins9")) kind=ETABINS9;
  else {
    std::cout << "DetermineEtaBinSet failed at <" << str  << ">\n";
    assert(0);
  }
  return kind;
}

// ------------------------------------------------------------------

//inline
DYTools::TDataKind_t DetermineDataKind(const TString &str) {
  using namespace DYTools;
  DYTools::TDataKind_t kind=DATA;
  if ((str.Contains("DATA") || str.Contains("data")) &&
      !(str.Contains("MC") || str.Contains("mc"))) kind=DATA;
  else if ((str.Contains("MC") || str.Contains("mc")) &&
	   !(str.Contains("DATA") || str.Contains("data"))) kind=MC;
  else {
    std::cout << "DetermineDataKind failed at <" << str << ">\n";
    assert(0);
  }
  return kind;
}

// ------------------------------------------------------------------
// ------------------------------------------------------------------

//inline
TString CrossSectionKindName(DYTools::TCrossSectionKind_t kind) {
  using namespace DYTools;
  TString name;
  switch(kind) {
  case _cs_None: name="none"; break;
  case _cs_spec: name="specCS"; break;
  case _cs_preFsr: name="preFsr"; break;
  case _cs_preFsrNorm: name="preFsrNorm"; break;
  case _cs_preFsrDet: name="preFsrDet"; break;
  case _cs_preFsrDetNorm: name="preFsrDetNorm"; break;
  case _cs_preFsrDetErr: name="preFsrDetErr"; break;
  case _cs_preFsrDetNormErr: name="preFsrDetNormErr"; break;
  case _cs_preFsrDetSystErr: name="preFsrDetSystErr"; break;
  case _cs_preFsrDetNormSystErr: name="preFsrDetNormSystErr"; break;
  case _cs_postFsr: name="postFsr"; break;
  case _cs_postFsrNorm: name="postFsrNorm"; break;
  case _cs_postFsrDet: name="postFsrDet"; break;
  case _cs_postFsrDetNorm: name="postFsrDetNorm"; break;
  default:
    std::cout << "CrossSectionKindName: cannot determine the name\n";
    assert(0);
  }
  return name;
}

// ------------------------------------------------------------------

//inline
DYTools::TCrossSectionKind_t DetermineCrossSectionKind(const TString &str) {
  using namespace DYTools;
  DYTools::TCrossSectionKind_t kind = _cs_None;
  if (str.Contains("preFsrNorm")) kind=_cs_preFsrNorm;
  else if (str.Contains("preFsrDetNormSystErr")) kind=_cs_preFsrDetNormSystErr;
  else if (str.Contains("preFsrDetNormErr")) kind=_cs_preFsrDetNormErr;
  else if (str.Contains("preFsrDetNorm")) kind=_cs_preFsrDetNorm;
  else if (str.Contains("preFsrDetSystErr")) kind=_cs_preFsrDetSystErr;
  else if (str.Contains("preFsrDetErr")) kind=_cs_preFsrDetErr;
  else if (str.Contains("preFsrDet")) kind=_cs_preFsrDet;
  else if (str.Contains("preFsr")) kind=_cs_preFsr;
  else if (str.Contains("postFsrNorm")) kind=_cs_postFsrNorm;
  else if (str.Contains("postFsrDetNorm")) kind=_cs_postFsrDetNorm;
  else if (str.Contains("postFsrDet")) kind=_cs_postFsrDet;
  else if (str.Contains("postFsr")) kind=_cs_postFsr;
  else if (str.Contains("specCS")) kind=_cs_spec;
  else if (str.Contains("none") || str.Contains("None")) kind=_cs_None;
  else {
    std::cout << "DetermineCrossSectionKind: cannod identify <" << str << ">\n";
    assert(0);
  }
  return kind;
}

// ------------------------------------------------------------------

//inline
TString CrossSectionKindLongName(DYTools::TCrossSectionKind_t kind) {
  using namespace DYTools;
  TString name;
  switch(kind) {
  case _cs_None: name="none"; break;
  case _cs_spec: name="special cross-section (user-defined)"; break;
  case _cs_preFsr: name="counts in full space"; break;
  case _cs_preFsrNorm: name="normalized cross section"; break;
  case _cs_preFsrDet: name="counts in DET space"; break;
  case _cs_preFsrDetNorm: name="normalized cross section (DET)"; break;
  case _cs_preFsrDetErr: name="error in DET space"; break;
  case _cs_preFsrDetNormErr: name="error of norm.cross section (DET)"; break;
  case _cs_preFsrDetSystErr: name="syst.error in DET space"; break;
  case _cs_preFsrDetNormSystErr: name="syst.error of norm.cross section (DET)"; break;
  case _cs_postFsr: name="postFsr counts in full space"; break;
  case _cs_postFsrNorm: name="normalized postFsr cross section"; break;
  case _cs_postFsrDet: name="postFsr counts in DET space"; break;
  case _cs_postFsrDetNorm: name="normalized postFsr cross section (DET)"; break;
  default:
    std::cout << "CrossSectionKindLongName: cannot determine the name\n";
    assert(0);
  }
  return name;
}

// ------------------------------------------------------------------

void AdjustFileNameEnding(TString &fname,DYTools::TSystematicsStudy_t systMode,
			  int iSeed) {
  TString ending;
  if (systMode==DYTools::FSR_RND_STUDY) {
    ending=Form("__FSR_RND_STUDY%d.root",iSeed);
  }
  else if (systMode==DYTools::PU_RND_STUDY) {
    ending=Form("__PU_RND_STUDY%d.root",iSeed);
  }
  if (ending.Length()) fname.ReplaceAll(".root",ending);
}

// ------------------------------------------------------------------
//     Printing
// ------------------------------------------------------------------
/*
inline std::ostream& operator<<(std::ostream& out, DYTools::TMassBinning_t massBinning) { out << MassBinningName(massBinning); return out; }
inline std::ostream& operator<<(std::ostream& out, DYTools::TSystematicsStudy_t study) { out << SystematicsStudyName(study); return out; }
inline std::ostream& operator<<(std::ostream& out, DYTools::TRunMode_t mode) { out << RunModeName(mode); return out; }
inline std::ostream& operator<<(std::ostream& out, DYTools::TTnPMethod_t method) { out << TnPMethodName(method); return out; }
inline std::ostream& operator<<(std::ostream& out, DYTools::TEfficiencyKind_t kind) { out << EfficiencyKindName(kind); return out; }
inline std::ostream& operator<<(std::ostream& out, DYTools::TEtBinSet_t set) { out << EtBinSetName(set); return out; }
inline std::ostream& operator<<(std::ostream& out, DYTools::TEtaBinSet_t set) { out << EtaBinSetName(set); return out; }
inline std::ostream& operator<<(std::ostream& out, DYTools::TDataKind_t kind) { out << DataKindName(kind); return out; }
inline std::ostream& operator<<(std::ostream& out, DYTools::TCrossSectionKind_t kind) { out << CrossSectionKindName(kind); return out; }
*/

// ------------------------------------------------------------------
// ------------------------------------------------------------------

namespace DYTools {
  //inline
  void printRunMode(TRunMode_t mode) {
    std::cout << "\nRun Mode: ";
    switch(mode) {
    case NORMAL_RUN: std::cout << "Normal execution\n"; break;
    case LOAD_DATA: std::cout << "Selected data will be loaded\n"; break;
    case DEBUG_RUN: std::cout << "\tDEBUG MODE is ON\n\n"; break;
    case DEBUG_LOAD: std::cout << "\tDEBUG_LOAD mode is ON\n\n"; break;
    default:
      std::cout << "\n\nprintRunMode(" << mode << ")\n";
    }
    return;
  }
};

// ------------------------------------------------------------------

namespace DYTools {
  //inline
  void printSystMode(TSystematicsStudy_t study) {
    std::cout << "Sytematics: " << SystematicsStudyName(study) << "\n";
    return;
  }
};

// ------------------------------------------------------------------

namespace DYTools {
  //inline
  void printExecMode(TRunMode_t mode, TSystematicsStudy_t study) {
    printRunMode(mode);
    printSystMode(study);
    return;
  }
};

// ------------------------------------------------------------------

namespace DYTools {
  //inline
  bool checkSystMode(DYTools::TSystematicsStudy_t mode, int debug_print, int nEntries, ...) {
    if (debug_print) std::cout << "\n\ncheckSystMode: Syst mode = " << SystematicsStudyName(mode) << "\n";
    bool ok=false;
    va_list vl;
    va_start(vl,nEntries);
    for (int i=0; !ok && (i<nEntries); ++i) {
      DYTools::TSystematicsStudy_t allowed= DYTools::TSystematicsStudy_t(va_arg(vl,int));
      if (debug_print) std::cout << " allowed#" << (i+1) << " is " << SystematicsStudyName(allowed) << "\n";
      ok=(allowed==mode) ? true : false;
    }
    va_end(vl);
    if (debug_print) std::cout << "answer is " << ((ok) ? "true":"false") << "\n";
    return ok;
  }
};

// ------------------------------------------------------------------

namespace DYTools {
  //inline
  bool checkCSKind(DYTools::TCrossSectionKind_t kind, int debug_print, int nEntries, ...) {
    if (debug_print) std::cout << "\n\ncheckCSKind: Syst mode = " << CrossSectionKindName(kind) << "\n";
    bool ok=false;
    va_list vl;
    va_start(vl,nEntries);
    for (int i=0; !ok && (i<nEntries); ++i) {
      DYTools::TCrossSectionKind_t allowed= DYTools::TCrossSectionKind_t(va_arg(vl,int));
      if (debug_print) std::cout << " allowed#" << (i+1) << " is " << CrossSectionKindName(allowed) << "\n";
      ok=(allowed==kind) ? true : false;
    }
    va_end(vl);
    if (debug_print) std::cout << "answer is " << ((ok) ? "true":"false") << "\n";
    return ok;
  }
};

// ------------------------------------------------------------------
/*
namespace DYTools {
  //inline
  bool checkSystMode(DYTools::TSystematicsStudy_t mode, int nEntries, const DYTools::TSystematicsStudy_t *allowedModes, int print) {
    if (print) std::cout << "\n\ncheckSystMode: Syst mode = " << SystematicsStudyName(mode) << "\n";
    bool ok=false;
    for (int i=0; !ok && (i<nEntries); ++i)
      ok=(allowedModes[i]==mode) ? true : false;
    return ok;
  }
};

// ------------------------------------------------------------------

namespace DYTools {
  //inline
  bool checkSystMode(DYTools::TSystematicsStudy_t mode, const std::vector<DYTools::TSystematicsStudy_t> &allowedModes, int print) {
    if (print) std::cout << "\n\ncheckSystMode: Syst mode = " << SystematicsStudyName(mode) << "\n";
    bool ok=false;
    for (unsigned int i=0; !ok && (i<allowedModes.size()); ++i)
      ok=(allowedModes[i]==mode) ? true : false;
    return ok;
  }
};
*/
// ------------------------------------------------------------------
// ------------------------------------------------------------------


#endif
