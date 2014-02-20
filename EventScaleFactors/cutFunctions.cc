#if !defined(__CINT__) || defined(__MAKECINT__)
#include "../EventScaleFactors/cutFunctions.hh"
#include "../Include/DYToolsUI.hh"
#endif

Bool_t dielectronMatchedToGeneratorLevel(const mithep::TGenInfo *gen, const mithep::TDielectron *dielectron){

  Bool_t result = kTRUE;

  // Matching is done as prescribed in AN-12-116
  //   - dR < 0.2
  //   - match to status 1 (post-FSR)
  //   - no charge matching

  // In the generator branch of this ntuple, first particle is always
  // negative, and second always positive. In the Dielectron block
  // of the ntuple, the first particle is always the one with larger Pt.
  double dR1_1=999, dR2_1=999;
  double dR1_2=999, dR2_2=999;
  TLorentzVector v1reco, v2reco, v1gen, v2gen;
  v1reco.SetPtEtaPhiM(dielectron->pt_1, dielectron->eta_1, dielectron->phi_1, 0.000511);
  v2reco.SetPtEtaPhiM(dielectron->pt_2, dielectron->eta_2, dielectron->phi_2, 0.000511);
  v1gen .SetPtEtaPhiM(gen->pt_1, gen->eta_1, gen->phi_1, 0.000511);
  v2gen .SetPtEtaPhiM(gen->pt_2, gen->eta_2, gen->phi_2, 0.000511);
  // Try both assignments, at least one assignment should match
  dR1_1 = v1reco.DeltaR(v1gen);
  dR2_1 = v2reco.DeltaR(v2gen);
  dR1_2 = v1reco.DeltaR(v2gen);
  dR2_2 = v2reco.DeltaR(v1gen);
  // Require that both are within the required cone
  bool matchAssignment1 = (fabs(dR1_1) < 0.2 && fabs(dR2_1) < 0.2 );
  bool matchAssignment2 = (fabs(dR1_2) < 0.2 && fabs(dR2_2) < 0.2 );
  if( ! (matchAssignment1 || matchAssignment2) ) result = kFALSE; 
  
  return result;
}

Bool_t electronMatchedToGeneratorLevel(const mithep::TGenInfo *gen, const mithep::TElectron *electron){
  
  Bool_t result = kTRUE;
  
  // Matching is done as prescribed in AN-12-116
  //   - dR < 0.2
  //   - match to status 1 (post-FSR)
  //   - no charge matching

  // In the generator branch of this ntuple, first particle is always
  // negative, and second always positive (but this is not used at present
  // as there is no charge matching).
  double dR1=999;
  double dR2=999;
  TLorentzVector vreco, v1gen, v2gen;
  vreco.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, 0.000511);
  v1gen .SetPtEtaPhiM(gen->pt_1, gen->eta_1, gen->phi_1, 0.000511);
  v2gen .SetPtEtaPhiM(gen->pt_2, gen->eta_2, gen->phi_2, 0.000511);
  
  dR1 = vreco.DeltaR(v1gen);
  dR2 = vreco.DeltaR(v2gen);

  if( !( fabs(dR1) < 0.2 || fabs(dR2) < 0.2 ) ) result = kFALSE; 
  
  return result;
}

Bool_t scMatchedToGeneratorLevel(const mithep::TGenInfo *gen, const mithep::TPhoton *sc){

  Bool_t result = kTRUE;
  // We do not know which of the gen electrons possibly
  // produced this supercluster, so we check both.
  double dR1=999, dR2=999;
  TLorentzVector vreco, v1gen, v2gen;
  vreco.SetPtEtaPhiM(sc->pt, sc->eta, sc->phi, 0.000511);
  v1gen .SetPtEtaPhiM(gen->pt_1, gen->eta_1, gen->phi_1, 0.000511);
  v2gen .SetPtEtaPhiM(gen->pt_2, gen->eta_2, gen->phi_2, 0.000511);
  dR1 = vreco.DeltaR(v1gen);
  dR2 = vreco.DeltaR(v2gen);
  // Require that at least one is  within dR of 0.2, otherwise bail out
  if( fabs(dR1) > 0.2 && fabs(dR2) > 0.2 ) result = kFALSE; 
  
  return result;
}

bool passID(const mithep::TElectron *electron, double rho){

  bool result = false;
#ifdef DYee8TeV
  result = passEGMID2012(electron, WP_MEDIUM, rho);
#endif
#ifdef DYee8TeV_reg
  result = passEGMID2012(electron, WP_MEDIUM, rho);
#endif
#ifdef DYee7TeV
  result = passEGMID2011(electron, WP_MEDIUM, rho);
#endif
  return result;
}

bool passIDTag(const mithep::TElectron *electron, double rho){

  bool result = false;
#ifdef DYee8TeV
  result = passEGMID2012(electron, WP_TIGHT, rho);
#endif
#ifdef DYee8TeV_reg
  result = passEGMID2012(electron, WP_TIGHT, rho);
#endif
#ifdef DYee7TeV
  result = passEGMID2011(electron, WP_TIGHT, rho);
#endif
  return result;
}

bool isTag(const mithep::TElectron *electron, ULong_t trigger, double rho){

  bool elePassID  = passIDTag(electron, rho);
  bool elePassHLT =  (electron ->hltMatchBits & trigger);
  bool notInGap = ! DYTools::isEcalGap( electron->scEta );  
  bool result = ( elePassID && elePassHLT && notInGap && (electron->pt > 25) );

  return result;
}

// -----------------------------------------------------
// to select for tag-systematics studies
// reduce the requirements for the tag when selecting the events
// systMode=RESOLUTION_STUDY : lower tag pt cut
// systMODE=FSR_STUDY : mediumID instead of tightID

bool isTag_systStudy(const mithep::TElectron *electron, ULong_t trigger, double rho, DYTools::TSystematicsStudy_t systMode){
  double elePtCut= (systMode==DYTools::RESOLUTION_STUDY) ? 20. : 25.;
  if (electron->pt <= elePtCut) return false;

  bool elePassHLT =  (electron ->hltMatchBits & trigger);
  bool notInGap = ! DYTools::isEcalGap( electron->scEta );  
  if ( !elePassHLT || !notInGap) return false;

  bool elePassID=false;
  if (systMode!=DYTools::FSR_STUDY) elePassID= passIDTag(electron, rho);
  else {
    // lower the ID requirement  
#ifdef DYee8TeV
    elePassID = passEGMID2012(electron, WP_MEDIUM, rho);
#endif
#ifdef DYee8TeV_reg
    elePassID = passEGMID2012(electron, WP_MEDIUM, rho);
#endif
#ifdef DYee7TeV
    elePassID = passEGMID2011(electron, WP_MEDIUM, rho);
#endif
  }

  bool result = ( elePassID && elePassHLT && notInGap );

  return result;
}

// -------------------------------------------------------------------

TString getLabel(int sample, DYTools::TEfficiencyKind_t effType, int method,  DYTools::TEtBinSet_t etBinning, DYTools::TEtaBinSet_t etaBinning, const TriggerSelection_t &trigSet){
  using namespace DYTools;

  TString label = analysisTag;
  if (analysisTag.Length()>0) label.Append("_");

  assert ( trigSet.isDefined() );
  if (sample != -1111) {
    label+= trigSet.triggerSetName();
  }

  if(sample == DATA)
    label += "_data";
  else if((sample == MC) || (sample == -1111))
    label += "_mc";
  else
    assert(0);

  label+= EfficiencyKindName( effType );
  //if(effType == HLT ) {
  //  if ((sample==DATA) && trigSet.hltEffMethodIs2011New()) label += "_hlt2011new";
  //}

  if (sample != -1111) {
    if(method == COUNTnCOUNT)
      label += "_count-count";
    else if( method == COUNTnFIT ) 
      label += "_count-fit";
    else if( method == FITnFIT ) 
      label += "_fit-fit";
    else
      assert(0);
  }

  const int old_style=0;
  if (old_style) {
    label += "_bins-et";
    label += getNEtBins(etBinning);
    label += "-eta";
    label += getNEtaBins(etaBinning);
  }
  else {
    label += EtBinSetName(etBinning);
    label += EtaBinSetName(etaBinning);
  }

  return label;
}

// -------------------------------------------------------------------
