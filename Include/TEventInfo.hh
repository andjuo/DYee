#ifndef TEventInfo_HH
#define TEventInfo_HH

#include "DYTools.hh"
#include <TObject.h>

// --------------- 7 TeV analysis -------------------

#ifdef DYee7TeV

namespace mithep 
{
  class TEventInfo : public TObject
  {
    public:
      TEventInfo():
      runNum(0), evtNum(0), lumiSec(0), nPU(0), nPUminus(0), nPUplus(0),
      pvx(0), pvy(0), pvz(0), bsx(0), bsy(0), bsz(0), 
      pfMET(0), pfMETphi(0), pfSumET(0), trkMET(0), trkMETphi(0), trkSumET(0),
      rhoLowEta(0), rhoHighEta(0), triggerBits(0), hasGoodPV(0)		    
      {}
      ~TEventInfo(){}

      UInt_t  runNum; 			             // run number in data
      UInt_t  evtNum; 			             // event number in data
      UInt_t  lumiSec;			             // lumi section
      UInt_t  nPU, nPUminus, nPUplus;                // number of reconstructed pile up vertices in event (MC only)       
      Float_t pvx, pvy, pvz;		             // best primary vertex
      Float_t bsx, bsy, bsz;		             // beamspot
      Float_t pfMET, pfMETphi, pfSumET;	             // particle flow MET
      Float_t trkMET, trkMETphi, trkSumET;           // track MET
      Float_t rhoLowEta, rhoHighEta;                 // average energy density for isolation correction
      ULong_t triggerBits;		             // HLT trigger bits
      Bool_t  hasGoodPV;                             // event has a good PV?

    ClassDef(TEventInfo,2)
  };
}

#endif

// --------------- 8 TeV analysis -------------------

#ifdef DYee8TeV
//#include "Inc_8TeV/TEventInfo.hh"

namespace mithep 
{
  class TEventInfo : public TObject
  {
    public:
      TEventInfo():
      runNum(0), evtNum(0), lumiSec(0), nPU(0), nPUminus(0), nPUplus(0),
      pvx(0), pvy(0), pvz(0), bsx(0), bsy(0), bsz(0), 
      pfMET(0), pfMETphi(0), pfSumET(0), trkMET(0), trkMETphi(0), trkSumET(0),
      rhoLowEta(0), rhoHighEta(0), triggerBits(0), hasGoodPV(0)		    
      {}
      ~TEventInfo(){}

      UInt_t  runNum; 			             // run number in data
      UInt_t  evtNum; 			             // event number in data
      UInt_t  lumiSec;			             // lumi section
      UInt_t  nPU, nPUminus, nPUplus;                // number of reconstructed pile up vertices in event (MC only)       
      Float_t pvx, pvy, pvz;		             // best primary vertex
      Float_t bsx, bsy, bsz;		             // beamspot
      Float_t pfMET, pfMETphi, pfSumET;	             // particle flow MET
      Float_t trkMET, trkMETphi, trkSumET;           // track MET
      Float_t rhoLowEta, rhoHighEta;                 // average energy density for isolation correction
      ULong_t triggerBits;		             // HLT trigger bits
      Bool_t  hasGoodPV;                             // event has a good PV?

    ClassDef(TEventInfo,2)
  };
}

#endif

// --------------- 8 TeV analysis (regressed) -------------------

#ifdef DYee8TeV_reg

namespace mithep 
{
  class TEventInfo : public TObject
  {
    public:
      TEventInfo():
      runNum(0), evtNum(0), lumiSec(0),
      nPU(0), nPUminus(0), nPUplus(0),
      nPUmean(0), nPUmeanminus(0), nPUmeanplus(0),
      pvx(0), pvy(0), pvz(0), bsx(0), bsy(0), bsz(0), 
      pfMET(0), pfMETphi(0), pfSumET(0), trkMET(0), trkMETphi(0), trkSumET(0),
      rhoLowEta(0), rhoHighEta(0), triggerBits(0), hasGoodPV(0)		    
      {}
      ~TEventInfo(){}

      UInt_t  runNum; 			             // run number in data
      UInt_t  evtNum; 			             // event number in data
      UInt_t  lumiSec;			             // lumi section
      UInt_t  nPU, nPUminus, nPUplus;                // number of generated pile up vertices in event (MC only)       
      Float_t nPUmean, nPUmeanminus, nPUmeanplus;    // number of expected pile up vertices in event (MC only)
      Float_t pvx, pvy, pvz;		             // best primary vertex
      Float_t bsx, bsy, bsz;		             // beamspot
      Float_t pfMET, pfMETphi, pfSumET;	             // particle flow MET
      Float_t trkMET, trkMETphi, trkSumET;           // track MET
      Float_t rhoLowEta, rhoHighEta;                 // average energy density for isolation correction, jet energy correction
      ULong_t triggerBits;		             // HLT trigger bits
      Bool_t  hasGoodPV;                             // event has a good PV?

    ClassDef(TEventInfo,3)
  };
}
#endif


// ------------------------------------------------------


#endif
