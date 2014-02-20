#if !defined(__CINT__) || defined(__MAKECINT__)

#include <TLorentzVector.h>         // class for Lorentz vector computations

// define classes and constants to read in ntuple
#include "../Include/DYTools.hh"
#include "../Include/TGenInfo.hh"
#include "../Include/TDielectron.hh"
#include "../Include/TElectron.hh"
#include "../Include/TPhoton.hh"
#include "../Include/DYToolsUI.hh"
#include "../Include/EleIDCuts.hh"
#include "../Include/TriggerSelection.hh"

#endif

Bool_t dielectronMatchedToGeneratorLevel(const mithep::TGenInfo *gen, const mithep::TDielectron *dielectron);

Bool_t electronMatchedToGeneratorLevel(const mithep::TGenInfo *gen, const mithep::TElectron *electron);

Bool_t scMatchedToGeneratorLevel(const mithep::TGenInfo *gen, const mithep::TPhoton *sc);

bool passID(const mithep::TElectron *electron, double rho);

bool isTag(const mithep::TElectron *electron, ULong_t trigger, double rho);

// to select for tag-systematics studies
// reduce the requirements for the tag when selecting the events
// systMode=RESOLUTION_STUDY : lower tag pt cut
// systMODE=FSR_STUDY : mediumID instead of tightID
bool isTag_systStudy(const mithep::TElectron *electron, ULong_t trigger, double rho, DYTools::TSystematicsStudy_t systMode);

TString getLabel(int sample, DYTools::TEfficiencyKind_t effType, int method, 
		 DYTools::TEtBinSet_t etBinning, 
		 DYTools::TEtaBinSet_t etaBinning, 
		 const TriggerSelection_t &trigSet);

