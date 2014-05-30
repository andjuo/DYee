#ifndef EventSelector_HH
#define EventSelector_HH

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <iostream>
#include <ostream>
#include <string>

#include "../Include/DYTools.hh"
#include "../Include/ElectronEnergyScale.hh"
#include "../Include/TDielectron.hh"
#include "../Include/TVertex.hh"
#include "../Include/EleIDCuts.hh"
#include "../Include/EventCounter.h"
#include "../Include/TriggerSelection.hh"
#include "../Include/BaseClass.hh"
#include "../Include/JsonParser.hh"
#include "../Include/CPlot.hh"

#include "../Include/AccessOrigNtuples.hh"
#include "../Include/InputFileMgr.hh"


// -------------------------------------------------
// -------------------------------------------------

namespace EventSelector {
  typedef enum { _selectNone, _selectDefault, _selectHLTOrdered, 
		 _selectLowerEtCut,
#ifdef DYee8TeV_reg
		 _selectHLTOrdered_nonRegressed, // use non-regressed values
#endif
		 _selectZeeSignalGen, _selectZeeSignalReco } TSelectionType_t;
  typedef enum { _escaleNone, _escaleUncorrected, _escaleData, _escaleDataRnd, _escaleMC } TEScaleCorrection_t;

  std::string selectionName(TSelectionType_t selection) {
    std::string name;
    switch(selection) {
    case _selectNone: name="None"; break;
    case _selectHLTOrdered: name="ordered"; break;
    case _selectZeeSignalGen: name="zeeSignalGen"; break;
    case _selectZeeSignalReco: name="zeeSignalReco"; break;
    case _selectDefault: name="default"; break;
    default: name="UNKNOWN";
    }
    return name;
  }
};

// -------------------------------------------------
// -------------------------------------------------

class EventSelector_t : BaseClass_t {
public:
protected:
  EventSelector::TSelectionType_t fSelection;
  EventSelector::TEScaleCorrection_t fEScaleCorrType;
  ElectronEnergyScale *fEScale;
  TriggerSelection_t fTrigger;
  EventCounter_t fEC;
  int fEScaleOwner;
public:
  EventSelector_t(); // use together with init(....)

  // most complete constructor
  EventSelector_t(InputFileMgr_t &mgr,
		  DYTools::TRunMode_t runMode,
		  DYTools::TSystematicsStudy_t systMode,
		  const TString& extraTag,
		  const TString& plotsExtraTag,
		  EventSelector::TSelectionType_t set_selection); // e.g. EventSelector::_selectDefault 

  EventSelector_t(EventSelector::TSelectionType_t set_selection, 
		  const TriggerSelection_t &set_trigger,
		  ElectronEnergyScale *set_escale=NULL);

  EventSelector_t(EventSelector::TSelectionType_t set_selection, 
		  const TString &trigger_string,
		  ElectronEnergyScale *set_escale=NULL);

  EventSelector_t(const EventSelector_t &es, ElectronEnergyScale *set_escale=NULL);

  ~EventSelector_t();

  // initialization to be used with an empty constructor
  int init(InputFileMgr_t &mgr, 
	   DYTools::TRunMode_t runMode, DYTools::TSystematicsStudy_t systMode, 
	   const TString &extraTag="", 
	   const TString &plotsExtraTag="",
	   EventSelector::TSelectionType_t set_selection= EventSelector::_selectDefault);

  int initEScale(const TString &escaleTag, int printEScale);

  void noEScaleCorrection() { fEScaleCorrType=EventSelector::_escaleNone; }
  int setEScaleCorrectionType(DYTools::TDataKind_t dataKind, DYTools::TSystematicsStudy_t systMode);

  int setEScaleCorrectionType(bool isData, DYTools::TSystematicsStudy_t systMode) {
    return this->setEScaleCorrectionType((isData) ? DYTools::DATA : DYTools::MC, systMode);
  }

  const TriggerSelection_t& trigger() const { return fTrigger; }
  void setTriggerActsOnData(bool yes) { fTrigger.actOnData(yes); }
  void setTriggerActsOnData(DYTools::TDataKind_t dataKind) { fTrigger.actOnData((dataKind==DYTools::DATA) ? true : false); }
  bool triggerActsOnData() const { return fTrigger.actOnData(); }

  bool eventTriggerOk(const AccessOrigNtuples_t &accessInfo) const {
    return fTrigger.matchEventTriggerBit(accessInfo.evtInfoPtr()->triggerBits,accessInfo.evtInfoPtr()->runNum);
  }

  TString &editECName() { return fEC.editName(); }

  bool inAcceptance(const AccessOrigNtuples_t &accessInfo) const {
    const mithep::TGenInfo *gen=accessInfo.genPtr();
    return DYTools::goodEtEtaPair( gen->pt_1, gen->eta_1, gen->pt_2, gen->eta_2 );
  }

  bool inAcceptancePreFsr(const AccessOrigNtuples_t &accessInfo) const {
    const mithep::TGenInfo *gen=accessInfo.genPtr();
    return DYTools::goodEtEtaPair( gen->vpt_1, gen->veta_1, gen->vpt_2, gen->veta_2 );
  }

  const ElectronEnergyScale *escale() const { return fEScale; }

  // ----------------------------
  // auxiliary functions

  TString generateFullTag(DYTools::TRunMode_t runMode,
			  DYTools::TSystematicsStudy_t systMode,
			  TString externalExtraTag="",
			  int includeSystMode=0) const {
    //std::cout << "runMode=<" << RunModeName(runMode) << ">\n";
    TString tag= (runMode!=DYTools::LOAD_DATA) ? generateRunTag(runMode) : "";
    if (includeSystMode) {
      TString tag2= generateSystTag(systMode);
      //if (tag2.Length()) tag.Append(TString("_") + tag2);
      if (tag2.Length()) tag.Append(tag2);
    }

    if ((systMode==DYTools::ESCALE_STUDY) ||
	(systMode==DYTools::ESCALE_STUDY_RND) ||
	(systMode==DYTools::APPLY_ESCALE)) {
      if (!fEScale) return this->reportError("generateBaseTag: fEScale is not set");
      tag.Append("_");
      tag.Append(fEScale->calibrationSetShortName());
    }

    if (systMode==DYTools::FSR_RND_STUDY) {
      int idx=externalExtraTag.Index("fixed");
      if (idx>0) externalExtraTag.Remove(idx,externalExtraTag.Length());
    }
    if (externalExtraTag.Length()) tag.Append(TString("_") + externalExtraTag);
    return tag; 
  }

  // -----

  void SetPlotOutDir(DYTools::TRunMode_t runMode,
		     DYTools::TSystematicsStudy_t systMode,
		     const TString &externalExtraTag="",
		     int printOutDirName=1) const {
    int includeSystMode=1;
    TString tag= this->generateFullTag(runMode,systMode,externalExtraTag,includeSystMode);
    CPlot::sOutDir= TString("plots");
    if (tag.Length()) {
      if ((tag[0]!='_') && (tag[0]!='-')) tag.Prepend(TString("_"));
      CPlot::sOutDir.Append(tag);
    }
    if (CPlot::sOutDir.Index(DYTools::analysisTag)==-1) {
      int len=CPlot::sOutDir.Length();
      if (CPlot::sOutDir[len-1]!='_') CPlot::sOutDir.Append("_");
      CPlot::sOutDir.Append(DYTools::analysisTag);
    }
    if (printOutDirName) std::cout << "CPlot::sOutDir=<" << CPlot::sOutDir << ">\n";
    gSystem->mkdir(CPlot::sOutDir,true);
  }

  // -----

  // dielectron may be modified by escale corrections
  bool testDielectron_default(mithep::TDielectron *dielectron, 
			      const mithep::TEventInfo *evtInfo,
			      EventCounter_t *ec=NULL); // evtInfo is for eleID

  // dielectron may be modified by escale corrections
  // The Et cuts are lower by 5GeV
  bool testDielectron_lowerEtCut(mithep::TDielectron *dielectron, 
				 const mithep::TEventInfo *evtInfo,
				 EventCounter_t *ec=NULL); // evtInfo is for eleID

  // HLTordered= default, but:
  //   - the trigger matching is strict (lead,trail)
  //   - mass cut is not applied
  bool testDielectron_HLT_ordered(mithep::TDielectron *dielectron, 
				  const mithep::TEventInfo *evtInfo,
				  EventCounter_t *ec=NULL); // evtInfo is for eleID
#ifdef DYee8TeV_reg
  bool testDielectron_HLT_ordered_nonRegressed
                     (mithep::TDielectron *dielectron, 
		      const mithep::TEventInfo *evtInfo,
		      EventCounter_t *ec=NULL); // evtInfo is for eleID
#endif

  // dielectron may be modified by escale corrections
  //bool testDielectron_zeeSignalGen(mithep::TDielectron *dielectron, 
				//ULong_t leadingTriggerObjectBit, ULong_t trailingTriggerObjectBit, 
  //				const mithep::TEventInfo *evtInfo,
  //			   EventCounter_t *ec=NULL) { return false; }// evtInfo is for eleID

  // dielectron may be modified by escale corrections
  bool testDielectron(mithep::TDielectron *dielectron, 
		      const mithep::TEventInfo *evtInfo,
		      EventCounter_t *ec=NULL) {// evtInfo is for eleID
    bool ok=kFALSE;
    switch(fSelection) {
    case EventSelector::_selectNone: break;
    case EventSelector::_selectDefault: ok=testDielectron_default(dielectron,evtInfo,ec); break;
    case EventSelector::_selectLowerEtCut: ok=testDielectron_lowerEtCut(dielectron,evtInfo,ec); break;
    case EventSelector::_selectHLTOrdered: ok=testDielectron_HLT_ordered(dielectron,evtInfo,ec); break;
#ifdef DYee8TeV_reg
    case EventSelector::_selectHLTOrdered_nonRegressed: ok=testDielectron_HLT_ordered_nonRegressed(dielectron,evtInfo,ec); break;
#endif
    //case EventSelector::_selectZeeSignalGen: ok=testDielectron_zeeSignalGen(dielectron,evtInfo,ec); break;
    default:
      ok=kFALSE;
    }
    return ok;
  }

  // dielectron may be modified by escale corrections
  /*
  bool testDielectron(mithep::TDielectron *dielectron, TEScaleCorrection_t applyEScale, 
		      ULong_t leadingTriggerObjectBit, ULong_t trailingTriggerObjectBit,
		      double rho, 
		      eventCounter_t *ec=NULL) {// the last argument is rho for PU correction of the isolation
    bool ok=kFALSE;
    switch(fSelection) {
    case _selectNone: break;
    case _selectDefault: ok=testDielectron_default(dielectron,applyEScale,leadingTriggerObjectBit,trailingTriggerObjectBit,rho,ec); break;
    default:
      ok=kFALSE;
    }
    return ok;
  }
  */

  /*
  bool operator()(mithep::TDielectron *dielectron, TEScaleCorrection_t applyEScale, 
		  ULong_t leadingTriggerObjectBit, ULong_t trailingTriggerObjectBit,
		  double rho, eventCounter_t *ec=NULL) {
    return testDielectron(dielectron,applyEScale,
			  leadingTriggerObjectBit,trailingTriggerObjectBit,
			  rho,ec);
  }
  */

  // -----

  template<class tData_t>
  int smearEvent(tData_t *data) const {
    double smearCorr= fEScale->generateMCSmear(data->scEta_1,data->scEta_2);
    data->mass += smearCorr;
    return 1;
  }

  // -----

  template<class tData_t>
  int applyEScale(tData_t *data, int isData=0, int randomized=0) const {
    return fEScale->applyEnergyScale(data,isData,randomized);
  }
  // -----


  std::ostream& printCounts(std::ostream&);
  void printCounts() { printCounts(std::cout); }
};

// -------------------------------------------------
// -------------------------------------------------


// for compatibility with previous code

//pvArr: TClonesArray("mithep::TVertex");
inline
UInt_t countGoodVertices(const TClonesArray *pvArr) {
  UInt_t nGoodPV=0;
  for(Int_t ipv=0; ipv<pvArr->GetEntriesFast(); ipv++) {
    const mithep::TVertex *pv = (mithep::TVertex*)((*pvArr)[ipv]);                  
    if(pv->nTracksFit                        < 1)  continue;
    if(pv->ndof                              < 4)  continue;
    if(fabs(pv->z)                           > 24) continue;
    if(sqrt((pv->x)*(pv->x)+(pv->y)*(pv->y)) > 2)  continue;
    nGoodPV++;
  }
  return nGoodPV;
}

// -------------------------------------------------
// -------------------------------------------------


#endif
