#ifndef TGenPhoton_HH
#define TGenPhoton_HH


#include "../Include/DYTools.hh"
#include <TObject.h>

// --------------- 8 TeV analysis (regressed) -------------------

#ifdef DYee8TeV_reg

namespace mithep
{
  class TGenPhoton : public TObject
  {
  public:
    TGenPhoton():
      status(0),parentId(0),pt(0),eta(0),phi(0)
    {}
    ~TGenPhoton(){}

    Int_t   status;                    // particle status
    Int_t   parentId;                  // PDG ID of parent
    Float_t pt, eta, phi;              // kinematics

    ClassDef(TGenPhoton,1)
 };
};

#endif

// ------------------------------------------------------

#endif
