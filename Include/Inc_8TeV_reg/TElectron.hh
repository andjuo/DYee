#ifndef EWKANA_NTUPLER_TELECTRON_HH
#define EWKANA_NTUPLER_TELECTRON_HH

#include <TObject.h>

namespace mithep
{
  class TElectron : public TObject
  {
    public:
      TElectron():
      pt(0), ptUncorr(0), eta(0), phi(0), trkIso03(0), emIso03(0), hadIso03(0),	     
      chIso_00_01(0), chIso_01_02(0), chIso_02_03(0), chIso_03_04(0), chIso_04_05(0),
      gammaIso_00_01(0), gammaIso_01_02(0), gammaIso_02_03(0), gammaIso_03_04(0), gammaIso_04_05(0),
      neuHadIso_00_01(0), neuHadIso_01_02(0), neuHadIso_02_03(0), neuHadIso_03_04(0), neuHadIso_04_05(0),
      pfPt(0), pfEta(0), pfPhi(0), d0(0), dz(0), scEt(0), scEtUncorr(0), scEta(0), scPhi(0),
      ecalE(0), HoverE(0), EoverP(0), fBrem(0), deltaEtaIn(0), deltaPhiIn(0), sigiEtaiEta(0),
      partnerDeltaCot(0), partnerDist(0), mva(0), q(0), nExpHitsInner(0), scID(0), trkID(0),
      typeBits(0), hltMatchBits(0), isConv(0)
      {}
      ~TElectron(){}
    
    

    void replace2UncorrEn() {
      pt=ptUncorr;
      scEt=scEtUncorr;
    }

    void assign(const TElectron *e) {
      pt= e->pt; 
      ptUncorr= e->ptUncorr;
      eta= e->eta;
      phi= e->phi;
      trkIso03= e->trkIso03;
      emIso03= e->emIso03;
      hadIso03= e->hadIso03;
      chIso_00_01= e->chIso_00_01;
      chIso_01_02= e->chIso_01_02;
      chIso_02_03= e->chIso_02_03;
      chIso_03_04= e->chIso_03_04;
      chIso_04_05= e->chIso_04_05;
      gammaIso_00_01= e->gammaIso_00_01;
      gammaIso_01_02= e->gammaIso_01_02;
      gammaIso_02_03= e->gammaIso_02_03;
      gammaIso_03_04= e->gammaIso_03_04;
      gammaIso_04_05= e->gammaIso_04_05;
      neuHadIso_00_01= e->neuHadIso_00_01;
      neuHadIso_01_02= e->neuHadIso_01_02;
      neuHadIso_02_03= e->neuHadIso_02_03;
      neuHadIso_03_04= e->neuHadIso_03_04;
      neuHadIso_04_05= e->neuHadIso_04_05;
      pfPt= e->pfPt;
      pfEta= e->pfEta;
      pfPhi= e->pfPhi;
      d0= e->d0;
      dz= e->dz;
      scEt= e->scEt;
      scEtUncorr= e->scEtUncorr;
      scEta= e->scEta;
      scPhi= e->scPhi;
      ecalE= e->ecalE;
      HoverE= e->HoverE;
      EoverP= e->EoverP;
      fBrem= e->fBrem;
      deltaEtaIn= e->deltaEtaIn;
      deltaPhiIn= e->deltaPhiIn;
      sigiEtaiEta= e->sigiEtaiEta;
      partnerDeltaCot= e->partnerDeltaCot;
      partnerDist= e->partnerDist;
      mva= e->mva;
      q= e->q;
      nExpHitsInner= e->nExpHitsInner;
      scID= e->scID;
      trkID= e->trkID;
      typeBits= e->typeBits;
      hltMatchBits= e->hltMatchBits;
      isConv= e->isConv;
    }


    // ---------- Data fields -------------

    Float_t pt, ptUncorr, eta, phi;          // kinematics
      Float_t trkIso03;                        // track isolation
      Float_t emIso03;                         // ECAL-based isolation
      Float_t hadIso03;                        // HCAL-based isolation
      Float_t chIso_00_01;                     // Particle Flow charged isolation
      Float_t chIso_01_02;
      Float_t chIso_02_03;
      Float_t chIso_03_04;
      Float_t chIso_04_05;
      Float_t gammaIso_00_01;                  // Particle Flow gamma isolation
      Float_t gammaIso_01_02;
      Float_t gammaIso_02_03;
      Float_t gammaIso_03_04;
      Float_t gammaIso_04_05;
      Float_t neuHadIso_00_01;                 // Particle Flow neutral hadron isolation
      Float_t neuHadIso_01_02;
      Float_t neuHadIso_02_03;
      Float_t neuHadIso_03_04;
      Float_t neuHadIso_04_05;
      Float_t pfPt, pfEta, pfPhi;              // Matching Particle Flow candidate kinematics
      Float_t d0, dz;                          // impact parameter
      Float_t scEt, scEtUncorr, scEta, scPhi;  // supercluster
      Float_t ecalE;                           // ECAL energy
      Float_t HoverE;                          // H / E
      Float_t EoverP;                          // E / p
      Float_t fBrem;                           // brem fraction  
      Float_t deltaEtaIn;                      // eta difference between track (at vertex) and SC
      Float_t deltaPhiIn;                      // phi difference between track (at vertex) and SC
      Float_t sigiEtaiEta;                     // eta-width of shower in number of crystals
      Float_t partnerDeltaCot;                 // cot(theta) difference with conversion partner track	
      Float_t partnerDist;                     // distance in x-y plane to nearest conversion partner track
      Float_t mva;                             // MVA electron ID
      Int_t   q;                               // charge
      UInt_t  nExpHitsInner;                   // number of hits expected before first hit      	             
      UInt_t  scID;                            // supercluster ID (for matching to photon superclusters)
      UInt_t  trkID;                           // tracker track ID (for matching to muons)
      UInt_t  typeBits;                        // bits for electron type
      ULong_t hltMatchBits;                    // bits for matching with HLT primitives
      Bool_t  isConv;                          // is conversion? (vertexing method)
          
    ClassDef(TElectron,3)
  };
}
#endif
