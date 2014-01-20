#ifndef TVertex_HH
#define TVertex_HH

#include "DYTools.hh"
#include <TObject.h>

// --------------- 7 TeV analysis -------------------

#ifdef DYee7TeV

namespace mithep 
{
  class TVertex : public TObject
  {
    public:
      TVertex():
      nTracksFit(0), ndof(0), chi2(0), sumPt(0), x(0), y(0), z(0)
      {}
      ~TVertex(){}

      UInt_t  nTracksFit;
      Float_t ndof;      
      Float_t chi2;
      Float_t sumPt;
      Float_t x,y,z;					 

    ClassDef(TVertex,2)
  };
}
#endif

// --------------- 8 TeV analysis -------------------

#ifdef DYee8TeV

namespace mithep 
{
  class TVertex : public TObject
  {
    public:
      TVertex():
      nTracksFit(0), ndof(0), chi2(0), sumPt(0), x(0), y(0), z(0)
      {}
      ~TVertex(){}

      UInt_t  nTracksFit;
      Float_t ndof;      
      Float_t chi2;
      Float_t sumPt;
      Float_t x,y,z;					 

    ClassDef(TVertex,2)
  };
}

#endif

// --------------- 8 TeV analysis (regressed) -------------------

#ifdef DYee8TeV_reg

namespace mithep 
{
  class TVertex : public TObject
  {
    public:
      TVertex():
      nTracksFit(0), ndof(0), chi2(0), sumPt(0), x(0), y(0), z(0)
      {}
      ~TVertex(){}

      UInt_t  nTracksFit;
      Float_t ndof;      
      Float_t chi2;
      Float_t sumPt;
      Float_t x,y,z;					 

    ClassDef(TVertex,2)
  };
}
#endif


// --------------------------------------------------



#endif
