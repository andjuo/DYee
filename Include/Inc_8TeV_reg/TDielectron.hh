#ifndef EWKANA_NTUPLER_TDIELECTRON_HH
#define EWKANA_NTUPLER_TDIELECTRON_HH

#include <TObject.h>
#include <TLorentzVector.h>

namespace mithep
{
  class TDielectron : public TObject
  {
    public:
      TDielectron():
      mass(0), pt(0), y(0), phi(0),
      pt_1(0), ptUncorr_1(0), eta_1(0), phi_1(0), trkIso03_1(0), emIso03_1(0), hadIso03_1(0),
      chIso_00_01_1(0), chIso_01_02_1(0), chIso_02_03_1(0), chIso_03_04_1(0), chIso_04_05_1(0),
      gammaIso_00_01_1(0), gammaIso_01_02_1(0), gammaIso_02_03_1(0), gammaIso_03_04_1(0), gammaIso_04_05_1(0),
      neuHadIso_00_01_1(0), neuHadIso_01_02_1(0), neuHadIso_02_03_1(0), neuHadIso_03_04_1(0), neuHadIso_04_05_1(0),
      pfPt_1(0), pfEta_1(0), pfPhi_1(0), d0_1(0), dz_1(0), scEt_1(0), scEtUncorr_1(0), scEta_1(0), scPhi_1(0),
      ecalE_1(0), HoverE_1(0), EoverP_1(0), fBrem_1(0), deltaEtaIn_1(0), deltaPhiIn_1(0), sigiEtaiEta_1(0),
      partnerDeltaCot_1(0), partnerDist_1(0), mva_1(0), q_1(0), nExpHitsInner_1(0), scID_1(0), trkID_1(0),
      typeBits_1(0), hltMatchBits_1(0), isConv_1(0),
      pt_2(0), ptUncorr_2(0), eta_2(0), phi_2(0), trkIso03_2(0), emIso03_2(0), hadIso03_2(0),
      chIso_00_01_2(0), chIso_01_02_2(0), chIso_02_03_2(0), chIso_03_04_2(0), chIso_04_05_2(0),
      gammaIso_00_01_2(0), gammaIso_01_02_2(0), gammaIso_02_03_2(0), gammaIso_03_04_2(0), gammaIso_04_05_2(0),
      neuHadIso_00_01_2(0), neuHadIso_01_02_2(0), neuHadIso_02_03_2(0), neuHadIso_03_04_2(0), neuHadIso_04_05_2(0),
      pfPt_2(0), pfEta_2(0), pfPhi_2(0), d0_2(0), dz_2(0), scEt_2(0), scEtUncorr_2(0), scEta_2(0), scPhi_2(0),
      ecalE_2(0), HoverE_2(0), EoverP_2(0), fBrem_2(0), deltaEtaIn_2(0), deltaPhiIn_2(0), sigiEtaiEta_2(0),
      partnerDeltaCot_2(0), partnerDist_2(0), mva_2(0), q_2(0), nExpHitsInner_2(0), scID_2(0), trkID_2(0),
      typeBits_2(0), hltMatchBits_2(0), isConv_2(0)
      {}

    TDielectron(const TDielectron &d)  :
      TObject(d),
      mass(d.mass), pt(d.pt), y(d.y), phi(d.phi),
      // 1st ele
      pt_1(d.pt_1), ptUncorr_1(d.ptUncorr_1), eta_1(d.eta_1), phi_1(d.phi_1), 
      trkIso03_1(d.trkIso03_1), emIso03_1(d.emIso03_1), hadIso03_1(d.hadIso03_1),
      chIso_00_01_1(d.chIso_00_01_1), chIso_01_02_1(d.chIso_01_02_1),
      chIso_02_03_1(d.chIso_02_03_1), chIso_03_04_1(d.chIso_03_04_1),
      chIso_04_05_1(d.chIso_04_05_1),
      gammaIso_00_01_1(d.gammaIso_00_01_1), gammaIso_01_02_1(d.gammaIso_01_02_1),
      gammaIso_02_03_1(d.gammaIso_02_03_1), gammaIso_03_04_1(d.gammaIso_03_04_1),
      gammaIso_04_05_1(d.gammaIso_04_05_1),
      neuHadIso_00_01_1(d.neuHadIso_00_01_1), neuHadIso_01_02_1(d.neuHadIso_01_02_1),
      neuHadIso_02_03_1(d.neuHadIso_02_03_1), neuHadIso_03_04_1(d.neuHadIso_02_03_1),
      neuHadIso_04_05_1(d.neuHadIso_04_05_1),
      pfPt_1(d.pfPt_1), pfEta_1(d.pfEta_1), pfPhi_1(d.pfPhi_1),
      d0_1(d.d0_1), dz_1(d.dz_1), 
      scEt_1(d.scEt_1), scEtUncorr_1(d.scEtUncorr_1), 
      scEta_1(d.scEta_1), scPhi_1(d.scPhi_1),
      ecalE_1(d.ecalE_1), HoverE_1(d.HoverE_1), EoverP_1(d.EoverP_1),
      fBrem_1(d.fBrem_1), deltaEtaIn_1(d.deltaEtaIn_1),
      deltaPhiIn_1(d.deltaPhiIn_1), sigiEtaiEta_1(d.sigiEtaiEta_1),
      partnerDeltaCot_1(d.partnerDeltaCot_1),
      partnerDist_1(d.partnerDist_1),
      mva_1(d.mva_1), q_1(d.q_1),
      nExpHitsInner_1(d.nExpHitsInner_1), scID_1(d.scID_1), trkID_1(d.trkID_1),
      typeBits_1(d.typeBits_1), hltMatchBits_1(d.hltMatchBits_1),
      isConv_1(d.isConv_1),
      // 2nd ele
      pt_2(d.pt_2), ptUncorr_2(d.ptUncorr_2), eta_2(d.eta_2), phi_2(d.phi_2), 
      trkIso03_2(d.trkIso03_2), emIso03_2(d.emIso03_2), hadIso03_2(d.hadIso03_2),
      chIso_00_01_2(d.chIso_00_01_2), chIso_01_02_2(d.chIso_01_02_2),
      chIso_02_03_2(d.chIso_02_03_2), chIso_03_04_2(d.chIso_03_04_2),
      chIso_04_05_2(d.chIso_04_05_2),
      gammaIso_00_01_2(d.gammaIso_00_01_2), gammaIso_01_02_2(d.gammaIso_01_02_2),
      gammaIso_02_03_2(d.gammaIso_02_03_2), gammaIso_03_04_2(d.gammaIso_03_04_2),
      gammaIso_04_05_2(d.gammaIso_04_05_2),
      neuHadIso_00_01_2(d.neuHadIso_00_01_2), neuHadIso_01_02_2(d.neuHadIso_01_02_2),
      neuHadIso_02_03_2(d.neuHadIso_02_03_2), neuHadIso_03_04_2(d.neuHadIso_02_03_2),
      neuHadIso_04_05_2(d.neuHadIso_04_05_2),
      pfPt_2(d.pfPt_2), pfEta_2(d.pfEta_2), pfPhi_2(d.pfPhi_2),
      d0_2(d.d0_2), dz_2(d.dz_2), 
      scEt_2(d.scEt_2), scEtUncorr_2(d.scEtUncorr_2), 
      scEta_2(d.scEta_2), scPhi_2(d.scPhi_2),
      ecalE_2(d.ecalE_2), HoverE_2(d.HoverE_2), EoverP_2(d.EoverP_2),
      fBrem_2(d.fBrem_2), deltaEtaIn_2(d.deltaEtaIn_2),
      deltaPhiIn_2(d.deltaPhiIn_2), sigiEtaiEta_2(d.sigiEtaiEta_2),
      partnerDeltaCot_2(d.partnerDeltaCot_2),
      partnerDist_2(d.partnerDist_2),
      mva_2(d.mva_2), q_2(d.q_2),
      nExpHitsInner_2(d.nExpHitsInner_2), scID_2(d.scID_2), trkID_2(d.trkID_2),
      typeBits_2(d.typeBits_2), hltMatchBits_2(d.hltMatchBits_2),
      isConv_2(d.isConv_2)
    {}

      ~TDielectron(){} 

    void replace2UncorrEn(int check=0) {
      if (!check) {
	pt_1=ptUncorr_1;
	scEt_1=scEtUncorr_1;
	pt_2=ptUncorr_2;
	scEt_2=scEtUncorr_2;
      }
      TLorentzVector ele1,ele2;
      ele1.SetPtEtaPhiM(pt_1,eta_1,phi_1,0.000511);
      ele2.SetPtEtaPhiM(pt_2,eta_2,phi_2,0.000511);
      TLorentzVector diEle = ele1+ele2;
      mass= diEle.M();
      pt  = diEle.Pt();
      y   = diEle.Rapidity();
      phi = diEle.Phi();
   }

   
      Float_t mass, pt, y, phi;  // dielectron kinematics
      
      // leading electron
      Float_t pt_1, ptUncorr_1, eta_1, phi_1;      
      Float_t trkIso03_1;          
      Float_t emIso03_1;           
      Float_t hadIso03_1;          
      Float_t chIso_00_01_1;
      Float_t chIso_01_02_1;
      Float_t chIso_02_03_1;
      Float_t chIso_03_04_1;
      Float_t chIso_04_05_1;
      Float_t gammaIso_00_01_1;
      Float_t gammaIso_01_02_1;
      Float_t gammaIso_02_03_1;
      Float_t gammaIso_03_04_1;
      Float_t gammaIso_04_05_1;
      Float_t neuHadIso_00_01_1;
      Float_t neuHadIso_01_02_1;
      Float_t neuHadIso_02_03_1;
      Float_t neuHadIso_03_04_1;
      Float_t neuHadIso_04_05_1;
      Float_t pfPt_1, pfEta_1, pfPhi_1;        
      Float_t d0_1, dz_1;            
      Float_t scEt_1, scEtUncorr_1, scEta_1, scPhi_1;
      Float_t ecalE_1;
      Float_t HoverE_1;            
      Float_t EoverP_1;            
      Float_t fBrem_1;             
      Float_t deltaEtaIn_1;        
      Float_t deltaPhiIn_1;        
      Float_t sigiEtaiEta_1;       
      Float_t partnerDeltaCot_1;   
      Float_t partnerDist_1;       
      Float_t mva_1;            
      Int_t   q_1;                 
      UInt_t  nExpHitsInner_1;     
      UInt_t  scID_1;              
      UInt_t  trkID_1;        
      UInt_t  typeBits_1;     
      ULong_t hltMatchBits_1;      
      Bool_t  isConv_1;
            
      // lagging electron
      Float_t pt_2, ptUncorr_2, eta_2, phi_2;      
      Float_t trkIso03_2;          
      Float_t emIso03_2;           
      Float_t hadIso03_2;          
      Float_t chIso_00_01_2;
      Float_t chIso_01_02_2;
      Float_t chIso_02_03_2;
      Float_t chIso_03_04_2;
      Float_t chIso_04_05_2;
      Float_t gammaIso_00_01_2;
      Float_t gammaIso_01_02_2;
      Float_t gammaIso_02_03_2;
      Float_t gammaIso_03_04_2;
      Float_t gammaIso_04_05_2;
      Float_t neuHadIso_00_01_2;
      Float_t neuHadIso_01_02_2;
      Float_t neuHadIso_02_03_2;
      Float_t neuHadIso_03_04_2;
      Float_t neuHadIso_04_05_2;
      Float_t pfPt_2, pfEta_2, pfPhi_2;        
      Float_t d0_2, dz_2;            
      Float_t scEt_2, scEtUncorr_2, scEta_2, scPhi_2;
      Float_t ecalE_2;
      Float_t HoverE_2;            
      Float_t EoverP_2;            
      Float_t fBrem_2;             
      Float_t deltaEtaIn_2;        
      Float_t deltaPhiIn_2;        
      Float_t sigiEtaiEta_2;       
      Float_t partnerDeltaCot_2;   
      Float_t partnerDist_2;       
      Float_t mva_2;            
      Int_t   q_2;                 
      UInt_t  nExpHitsInner_2;     
      UInt_t  scID_2;              
      UInt_t  trkID_2;       
      UInt_t  typeBits_2;      
      ULong_t hltMatchBits_2;    
      Bool_t  isConv_2;

    ClassDef(TDielectron,3)
  };
}
#endif
