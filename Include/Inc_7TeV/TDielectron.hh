#ifndef EWKANA_NTUPLER_TDIELECTRON_HH
#define EWKANA_NTUPLER_TDIELECTRON_HH

#include <TObject.h>
#include "../Include/TElectron.hh"

namespace mithep
{
  class TDielectron : public TObject
  {
    public:
      TDielectron():
      mass(0), pt(0), y(0), phi(0),
      pt_1(0), eta_1(0), phi_1(0), trkIso03_1(0), emIso03_1(0), hadIso03_1(0),
      chIso_00_01_1(0), chIso_01_02_1(0), chIso_02_03_1(0), chIso_03_04_1(0), chIso_04_05_1(0),
      gammaIso_00_01_1(0), gammaIso_01_02_1(0), gammaIso_02_03_1(0), gammaIso_03_04_1(0), gammaIso_04_05_1(0),
      neuHadIso_00_01_1(0), neuHadIso_01_02_1(0), neuHadIso_02_03_1(0), neuHadIso_03_04_1(0), neuHadIso_04_05_1(0),
      pfPt_1(0), pfEta_1(0), pfPhi_1(0), d0_1(0), dz_1(0), scEt_1(0), scEta_1(0), scPhi_1(0),
      ecalE_1(0), HoverE_1(0), EoverP_1(0), fBrem_1(0), deltaEtaIn_1(0), deltaPhiIn_1(0), sigiEtaiEta_1(0),
      partnerDeltaCot_1(0), partnerDist_1(0), mva_1(0), q_1(0), nExpHitsInner_1(0), scID_1(0), trkID_1(0),
      typeBits_1(0), hltMatchBits_1(0), isConv_1(0),
      pt_2(0), eta_2(0), phi_2(0), trkIso03_2(0), emIso03_2(0), hadIso03_2(0),
      chIso_00_01_2(0), chIso_01_02_2(0), chIso_02_03_2(0), chIso_03_04_2(0), chIso_04_05_2(0),
      gammaIso_00_01_2(0), gammaIso_01_02_2(0), gammaIso_02_03_2(0), gammaIso_03_04_2(0), gammaIso_04_05_2(0),
      neuHadIso_00_01_2(0), neuHadIso_01_02_2(0), neuHadIso_02_03_2(0), neuHadIso_03_04_2(0), neuHadIso_04_05_2(0),
      pfPt_2(0), pfEta_2(0), pfPhi_2(0), d0_2(0), dz_2(0), scEt_2(0), scEta_2(0), scPhi_2(0),
      ecalE_2(0), HoverE_2(0), EoverP_2(0), fBrem_2(0), deltaEtaIn_2(0), deltaPhiIn_2(0), sigiEtaiEta_2(0),
      partnerDeltaCot_2(0), partnerDist_2(0), mva_2(0), q_2(0), nExpHitsInner_2(0), scID_2(0), trkID_2(0),
      typeBits_2(0), hltMatchBits_2(0), isConv_2(0)
      {}
      ~TDielectron(){} 
   
      Float_t mass, pt, y, phi;  // dielectron kinematics
      
      // leading electron
      Float_t pt_1, eta_1, phi_1;      
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
      Float_t scEt_1, scEta_1, scPhi_1;
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
      Float_t pt_2, eta_2, phi_2;      
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
      Float_t scEt_2, scEta_2, scPhi_2;
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

  public:

    void extractElectron(int index, mithep::TElectron &ele) {
      mithep::TDielectron *dielectron = this;
    
      if(index == 1){
	ele. pt                  = dielectron-> pt_1                 ;
	ele. eta                 = dielectron-> eta_1                ;
	ele. phi                 = dielectron-> phi_1                ;      
	
	ele. trkIso03            = dielectron-> trkIso03_1           ;
	ele. emIso03             = dielectron-> emIso03_1            ;
	ele. hadIso03            = dielectron-> hadIso03_1           ;
	ele. chIso_00_01         = dielectron-> chIso_00_01_1        ;
	ele. chIso_01_02         = dielectron-> chIso_01_02_1        ;
	ele. chIso_02_03         = dielectron-> chIso_02_03_1        ;
	ele. chIso_03_04         = dielectron-> chIso_03_04_1        ;
	ele. chIso_04_05         = dielectron-> chIso_04_05_1        ;

	ele. gammaIso_00_01      = dielectron-> gammaIso_00_01_1        ;
	ele. gammaIso_01_02      = dielectron-> gammaIso_01_02_1        ;
	ele. gammaIso_02_03      = dielectron-> gammaIso_02_03_1        ;
	ele. gammaIso_03_04      = dielectron-> gammaIso_03_04_1        ;
	ele. gammaIso_04_05      = dielectron-> gammaIso_04_05_1        ;
	
	ele. neuHadIso_00_01     = dielectron-> neuHadIso_00_01_1        ;
	ele. neuHadIso_01_02     = dielectron-> neuHadIso_01_02_1        ;
	ele. neuHadIso_02_03     = dielectron-> neuHadIso_02_03_1        ;
	ele. neuHadIso_03_04     = dielectron-> neuHadIso_03_04_1        ;
	ele. neuHadIso_04_05     = dielectron-> neuHadIso_04_05_1        ;
	
	ele.pfPt                 = dielectron->pfPt_1 ;
	ele.pfEta                = dielectron->pfEta_1 ;
	ele.pfPhi                = dielectron->pfPhi_1 ;

	ele. d0                  = dielectron-> d0_1                 ;
	ele. dz                  = dielectron-> dz_1                 ;
	ele. scEt                = dielectron-> scEt_1               ;
	ele. scEta               = dielectron-> scEta_1              ;
	ele. scPhi               = dielectron-> scPhi_1              ;
	ele. ecalE               = dielectron-> ecalE_1              ;
	ele. HoverE              = dielectron-> HoverE_1             ;
	ele. EoverP              = dielectron-> EoverP_1             ;
	ele. fBrem               = dielectron-> fBrem_1              ;
	ele. deltaEtaIn          = dielectron-> deltaEtaIn_1         ;
	ele. deltaPhiIn          = dielectron-> deltaPhiIn_1         ;
	ele. sigiEtaiEta         = dielectron-> sigiEtaiEta_1        ;
	ele. partnerDeltaCot     = dielectron-> partnerDeltaCot_1    ;
	ele. partnerDist         = dielectron-> partnerDist_1        ;
	ele.mva                  = dielectron->mva_1                 ;
	
	ele. q                   = dielectron-> q_1                  ;
	ele. nExpHitsInner       = dielectron-> nExpHitsInner_1      ;
	ele.  scID               = dielectron-> scID_1               ;
	ele.  trkID              = dielectron-> trkID_1              ;
	ele.  typeBits           = dielectron-> typeBits_1           ;
	
	ele. hltMatchBits        = dielectron-> hltMatchBits_1       ;
	ele. isConv              = dielectron-> isConv_1             ;
      }
      else{
	ele. pt                  = dielectron-> pt_2                 ;
	ele. eta                 = dielectron-> eta_2                ;
	ele. phi                 = dielectron-> phi_2                ;
	ele. trkIso03            = dielectron-> trkIso03_2           ;
	ele. emIso03             = dielectron-> emIso03_2            ;
	ele. hadIso03            = dielectron-> hadIso03_2           ;
	ele. chIso_00_01         = dielectron-> chIso_00_01_2        ;
	ele. chIso_01_02         = dielectron-> chIso_01_02_2        ;
	ele. chIso_02_03         = dielectron-> chIso_02_03_2        ;
	ele. chIso_03_04         = dielectron-> chIso_03_04_2        ;
	ele. chIso_04_05         = dielectron-> chIso_04_05_2        ;

	ele. gammaIso_00_01      = dielectron-> gammaIso_00_01_2        ;
	ele. gammaIso_01_02      = dielectron-> gammaIso_01_02_2        ;
	ele. gammaIso_02_03      = dielectron-> gammaIso_02_03_2        ;
	ele. gammaIso_03_04      = dielectron-> gammaIso_03_04_2        ;
	ele. gammaIso_04_05      = dielectron-> gammaIso_04_05_2        ;

	ele. neuHadIso_00_01     = dielectron-> neuHadIso_00_01_2        ;
	ele. neuHadIso_01_02     = dielectron-> neuHadIso_01_02_2        ;
	ele. neuHadIso_02_03     = dielectron-> neuHadIso_02_03_2        ;
	ele. neuHadIso_03_04     = dielectron-> neuHadIso_03_04_2        ;
	ele. neuHadIso_04_05     = dielectron-> neuHadIso_04_05_2        ;

	ele.pfPt                 = dielectron->pfPt_2 ;
	ele.pfEta                = dielectron->pfEta_2 ;
	ele.pfPhi                = dielectron->pfPhi_2 ;
	
	ele. d0                  = dielectron-> d0_2                 ;
	ele. dz                  = dielectron-> dz_2                 ;
	
	ele. scEt                = dielectron-> scEt_2               ;
	ele. scEta               = dielectron-> scEta_2              ;
	ele. scPhi               = dielectron-> scPhi_2              ;
	ele. ecalE               = dielectron-> ecalE_2              ;
	ele. HoverE              = dielectron-> HoverE_2             ;
	ele. EoverP              = dielectron-> EoverP_2             ;
	ele. fBrem               = dielectron-> fBrem_2              ;
	
	ele. deltaEtaIn          = dielectron-> deltaEtaIn_2         ;
	ele. deltaPhiIn          = dielectron-> deltaPhiIn_2         ;
	ele. sigiEtaiEta         = dielectron-> sigiEtaiEta_2        ;
	ele. partnerDeltaCot     = dielectron-> partnerDeltaCot_2    ;
	ele. partnerDist         = dielectron-> partnerDist_2        ;
	ele.mva                  = dielectron->mva_2                 ;
	
	ele. q                   = dielectron-> q_2                  ;
	ele. nExpHitsInner       = dielectron-> nExpHitsInner_2      ;
	ele.  scID               = dielectron-> scID_2               ;
	ele.  trkID              = dielectron-> trkID_2              ;
	ele.  typeBits           = dielectron-> typeBits_2           ;
	
	ele. hltMatchBits        = dielectron-> hltMatchBits_2       ;
	ele. isConv              = dielectron-> isConv_2             ;
      }

    }

    // ---------------------------------------

    ClassDef(TDielectron,2)

  };
}
#endif
