#include "EventSelector.hh"
#include <TLorentzVector.h>

// ---------------------------------------------------------------

EventSelector_t::EventSelector_t(InputFileMgr_t &mgr,
				 DYTools::TRunMode_t runMode,
				 DYTools::TSystematicsStudy_t systMode,
				 const TString& extraTag,
				 EventSelector::TSelectionType_t set_selection) :
  BaseClass_t("EventSelector_t"),
  fSelection(set_selection),
  fEScaleCorrType(EventSelector::_escaleNone),
  fEScale(NULL),
  fTrigger(mgr.triggerTag(),true),
  fEC("eventSelector"),
  fEScaleOwner(0) 
{
  const int printEScale=1;
  int res= this->initEScale(mgr.energyScaleTag(),printEScale);
  if (res) this->SetPlotOutDir(runMode,systMode,extraTag);

  DYTools::TRunMode_t runModeLocal=runMode;
  if (runMode==DYTools::DEBUG_LOAD) runModeLocal=DYTools::DEBUG_RUN;
  else if (runMode==DYTools::LOAD_DATA) runModeLocal=DYTools::NORMAL_RUN;
  
  if (res) {
    TString auto_tag=this->generateFullTag(runModeLocal,systMode,extraTag);
    if (auto_tag.Length()) auto_tag.Prepend("_");
    mgr.setNtupleNameExtraTag(auto_tag);
    mgr.setDirNameExtraTag(auto_tag);
  }
  if (!res) this->reportError("EventSelector_t::EventSelector_t(mgr)");
}

// ---------------------------------------------------------------

EventSelector_t::EventSelector_t(EventSelector::TSelectionType_t set_selection, 
				 const TriggerSelection_t &set_trigger,
				 ElectronEnergyScale *set_escale) :
  BaseClass_t("EventSelector_t"),
  fSelection(set_selection),
  fEScaleCorrType(EventSelector::_escaleNone),
  fEScale(set_escale),
  fTrigger(set_trigger),
  fEC("eventSelector"),
  fEScaleOwner(0)
{
}

// ---------------------------------------------------------------

EventSelector_t::EventSelector_t(EventSelector::TSelectionType_t set_selection, 
				 const TString &trigger_string_tag,
				 ElectronEnergyScale *set_escale) :
  BaseClass_t("EventSelector_t"),
  fSelection(set_selection),
  fEScaleCorrType(EventSelector::_escaleNone),
  fEScale(set_escale),
  fTrigger(trigger_string_tag,true),
  fEC("eventSelector"),
  fEScaleOwner(0)
{
}

// ---------------------------------------------------------------

EventSelector_t::EventSelector_t(const EventSelector_t &es,
				 ElectronEnergyScale *set_escale) :
  BaseClass_t("EventSelector_t"),
  fSelection(es.fSelection),
  fEScaleCorrType(es.fEScaleCorrType),
  fEScale(set_escale),
  fTrigger(fTrigger),
  fEC("eventSelector"),
  fEScaleOwner(0)
{
}

// ---------------------------------------------------------------

EventSelector_t::~EventSelector_t() {
  if (fEScaleOwner && fEScale) delete fEScale;
}

// ---------------------------------------------------------------

int EventSelector_t::initEScale(const TString &escaleTag, int print) {
  if (fEScaleOwner && fEScale) delete fEScale;
  fEScaleOwner=1;
  fEScale=new ElectronEnergyScale(escaleTag);
  int ok= (fEScale && fEScale->isInitialized()) ? 1:0;
  if (ok && print) fEScale->print();
  return ok;
}

// ---------------------------------------------------------------

int EventSelector_t::setEScaleCorrectionType(DYTools::TDataKind_t dataKind, DYTools::TSystematicsStudy_t systMode) {
  if (dataKind==DYTools::DATA) {
    switch (systMode) {
    case DYTools::NO_SYST: 
    case DYTools::ESCALE_STUDY:
      fEScaleCorrType=EventSelector::_escaleData;
      break;
    case DYTools::ESCALE_STUDY_RND:
      fEScaleCorrType=EventSelector::_escaleDataRnd;
      break;
    default:
      return this->reportError("setEScaleCorrectionType(%s,%s)",DataKindName(dataKind),SystematicsStudyName(systMode));
    }
  }
  else {
    fEScaleCorrType=EventSelector::_escaleNone;
    //return this->reportError("setEScaleCorrectionType(%s,%s)",DataKindName(dataKind),SystematicsStudyName(systMode));
  }
  return 1;
}

// ---------------------------------------------------------------

// dielectron may be modified by escale corrections
bool EventSelector_t::testDielectron_default(mithep::TDielectron *dielectron, 
					     //DYTools::TDataKind_t dataKind,
					     const mithep::TEventInfo *evtInfo,
					     EventCounter_t *ec) {// evtInfo is for eleID

  fEC.numDielectrons_inc();

  //
  // Energy scale corrections for data
  // NOTE: the electrons and dielectron 4-vectors are updated, the supercluster quantities are not
  //
  Double_t scEt1 = dielectron->scEt_1;
  Double_t scEt2 = dielectron->scEt_2;
  // Electron energy scale correction
  if((fEScaleCorrType==EventSelector::_escaleData) || 
     (fEScaleCorrType==EventSelector::_escaleDataRnd)) {
    if (!fEScale) {
      return this->reportError("testDielectron_default: escaleData is requested, but the pointer is null");
    }

    double corr1 = 1, corr2= 1;
    if (fEScaleCorrType==EventSelector::_escaleData) {
      corr1=fEScale->getEnergyScaleCorrection(dielectron->scEta_1);
      corr2=fEScale->getEnergyScaleCorrection(dielectron->scEta_2);
    }
    else {
      corr1=fEScale->getEnergyScaleCorrectionRandomized(dielectron->scEta_1);
      corr2=fEScale->getEnergyScaleCorrectionRandomized(dielectron->scEta_2);
    }
    scEt1 = dielectron->scEt_1 * corr1;
    scEt2 = dielectron->scEt_2 * corr2;
    
    TLorentzVector ele1; 
    ele1.SetPtEtaPhiM(dielectron->pt_1,dielectron->eta_1,dielectron->phi_1,0.000511);
    ele1 *= corr1;
    dielectron->pt_1  = ele1.Pt();
    dielectron->eta_1 = ele1.Eta();
    dielectron->phi_1 = ele1.Phi();
    
    TLorentzVector ele2; 
    ele2.SetPtEtaPhiM(dielectron->pt_2,dielectron->eta_2,dielectron->phi_2,0.000511);
    ele2 *= corr2;
    dielectron->pt_2  = ele2.Pt();
    dielectron->eta_2 = ele2.Eta();
    dielectron->phi_2 = ele2.Phi();
      
    TLorentzVector vDiEle = ele1+ele2;            
    dielectron->mass = vDiEle.M();
    dielectron->pt   = vDiEle.Pt();
    dielectron->y    = vDiEle.Rapidity();
    dielectron->phi  = vDiEle.Phi(); 
  }
       	  
  // requirements on BOTH electrons
  // For DY ET cuts are asymmetric.
  if( !DYTools::goodEtEtaPair(scEt1, dielectron->scEta_1,
			      scEt2, dielectron->scEta_2) ) return false;

  if (ec) ec->numDielectronsGoodEtEta_inc();
  fEC.numDielectronsGoodEtEta_inc();
  
  // Both electrons must match trigger objects. At least one ordering
  // must match
  if (! fTrigger.matchTriggerObjectBitAnyOrder(dielectron->hltMatchBits_1, dielectron->hltMatchBits_2, evtInfo->runNum) ) 
    return false;

  if (ec) ec->numDielectronsHLTmatched_inc();
  fEC.numDielectronsHLTmatched_inc();

  // Other cuts to both electrons
  if (!passEleID(dielectron,evtInfo)) return false;

  if (ec) ec->numDielectronsIDpassed_inc();
  fEC.numDielectronsIDpassed_inc();

  // loose mass window 
  double minMass= DYTools::massBinLimits[0];
  if (minMass>=10) {
    minMass= (fEScaleCorrType==EventSelector::_escaleUncorrected) ? 5 : 10;
  }
  if( dielectron->mass < minMass ) return false;

  if (ec) ec->numDielectronsGoodMass_inc();
  fEC.numDielectronsGoodMass_inc();

  // selection PASSED
  return true;
}


// ---------------------------------------------------------------

std::ostream& EventSelector_t::printCounts(std::ostream &out) {
  /*
  char buf[20];
  out << "DielectronSelector (selection=" << EventSelector_t::selectionName(fSelection) << ")\n";
  const char *format="%10u";
  sprintf(buf,format,fTotalCandidates);
  out << "Total number of dielectron candidates  " << buf << "\n";
  sprintf(buf,format,fCandidatesGoodEta);
  out << "Total candidates with good eta         " << buf << "\n";
  sprintf(buf,format,fCandidatesGoodEt);
  out << "Total candidates with good Et          " << buf << "\n";
  sprintf(buf,format,fCandidatesHLTMatched);
  out << "Total candidates HLT matched           " << buf << "\n";
  sprintf(buf,format,fCandidatesIDPassed);
  out << "Total candidates ID passed             " << buf << "\n";
  sprintf(buf,format,fCandidatesMassAboveMinLimit);
  out << "Total candidates with mass above limit " << buf << "\n";
  */
  return out;
}

// ---------------------------------------------------------------

// ---------------------------------------------------------------
