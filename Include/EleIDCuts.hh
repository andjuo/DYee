#ifndef EleIDCuts_HH
#define EleIDCuts_HH

#include "../Include/DYTools.hh"

// --------------- 7 TeV analysis -------------------
#ifdef DYee7TeV
#include "Inc_7TeV/EleIDCuts.hh"
#endif
// --------------- 8 TeV analysis -------------------
#ifdef DYee8TeV
#include "Inc_8TeV/EleIDCuts.hh"
#endif
// --------------- 8 TeV analysis (regressed) -------------------
#ifdef DYee8TeV_reg
#include "Inc_8TeV_reg/EleIDCuts.hh"
#endif
// --------------------------------------------------


//
// Define electron IDs
//

typedef enum { _EleID_Smurf2011,
	       _EleID_EGM2011_Medium,
	       _EleID_EGM2012_Medium} TEleID_t;

//const TEleID_t _electronID=_EleID_Smurf2011;
// Started working on this but not finished.

#ifdef DYee7TeV
const TEleID_t _electronID=_EleID_EGM2011_Medium;
#endif
#if (defined DYee8TeV) || (defined DYee8TeV_reg)
const TEleID_t _electronID=_EleID_EGM2012_Medium;
#endif


//
// Generic passEleID function
//

template<class EleObj_t>
Bool_t passEleID(const EleObj_t *electron, const mithep::TEventInfo *info=NULL) {
  Bool_t pass=kFALSE;
  switch(_electronID) {
#ifdef DYee7TeV
  case _EleID_Smurf2011: pass=passSmurf(electron); break;
  case _EleID_EGM2011_Medium: 
    assert(info);
    //pass=passEGMID2011(electron,_EleID_EGM2011_Medium,info->rhoLowEta);
    pass=passEGMID2011(electron,WP_MEDIUM,info->rhoLowEta);
    break;
#endif
#if (defined DYee8TeV) || (defined DYee8TeV_reg)
  case _EleID_EGM2012_Medium: 
    assert(info);
    pass=passEGMID2012(electron,WP_MEDIUM,info->rhoLowEta);
    break;
#endif
  default:
    std::cout << "passEleID: ElectronID is not prepared\n";
  }
  return pass;
}


#endif
