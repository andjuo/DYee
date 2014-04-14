#ifndef TDielectron_HH
#define TDielectron_HH

#include "../Include/DYTools.hh"
#include <TObject.h>
#include "../Include/TElectron.hh"
#include <TLorentzVector.h>

// --------------- 7 TeV analysis -------------------


#ifdef DYee7TeV

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

    void restoreEScaleModifiedValues(const mithep::TDielectron &di) {
      pt_1 = di.pt_1; eta_1=di.eta_1; phi_1=di.phi_1;
      pt_2 = di.pt_2; eta_2=di.eta_2; phi_2=di.phi_2;
      scEt_1 = di.scEt_1; scEta_1=di.scEta_1; scPhi_1=di.scPhi_1;
      scEt_2 = di.scEt_2; scEta_2=di.scEta_2; scPhi_2=di.scPhi_2;
      mass = di.mass;
      pt= di.pt;
      y = di.y;
      phi=di.phi;
    }

    void extractElectron(int index, mithep::TElectron &ele) const {
      const mithep::TDielectron *dielectron = this;
    
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

// --------------- 8 TeV analysis -------------------

#ifdef DYee8TeV

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

    void restoreEScaleModifiedValues(const mithep::TDielectron &di) {
      pt_1 = di.pt_1; eta_1=di.eta_1; phi_1=di.phi_1;
      pt_2 = di.pt_2; eta_2=di.eta_2; phi_2=di.phi_2;
      scEt_1 = di.scEt_1; scEta_1=di.scEta_1; scPhi_1=di.scPhi_1;
      scEt_2 = di.scEt_2; scEta_2=di.scEta_2; scPhi_2=di.scPhi_2;
      mass = di.mass;
      pt= di.pt;
      y = di.y;
      phi=di.phi;
    }

    void extractElectron(int index, mithep::TElectron &ele) const {
      const mithep::TDielectron *dielectron = this;
    
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
	ele. ecalE               = dielectron-> ecalE_1             ;
	ele. HoverE              = dielectron-> HoverE_1             ;
	ele. EoverP              = dielectron-> EoverP_1             ;
	ele. fBrem               = dielectron-> fBrem_1             ;

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
      }else{
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
	ele. ecalE               = dielectron-> ecalE_2             ;
	ele. HoverE              = dielectron-> HoverE_2             ;
	ele. EoverP              = dielectron-> EoverP_2             ;
	ele. fBrem               = dielectron-> fBrem_2             ;
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


    ClassDef(TDielectron,2)
  };
}

#endif

// --------------------------------------------------

// --------------- 8 TeV analysis (regressed) -------------------

#ifdef DYee8TeV_reg

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

    void replace2UncorrEn(int check=0, double scale=1.0) {
      if (!check) {
	pt_1  =scale*(ptUncorr_1 - pt_1) + pt_1;
	scEt_1=scale*(scEtUncorr_1 - scEt_1 ) + scEt_1;
	pt_2  =scale*(ptUncorr_2 - pt_2) + pt_2;
	scEt_2=scale*(scEtUncorr_2 - scEt_2 ) + scEt_2;
      }
      TLorentzVector ele1, ele2;
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

  public:

    void restoreEScaleModifiedValues(const mithep::TDielectron &di) {
      pt_1 = di.pt_1; eta_1=di.eta_1; phi_1=di.phi_1;
      pt_2 = di.pt_2; eta_2=di.eta_2; phi_2=di.phi_2;
      scEt_1 = di.scEt_1; scEta_1=di.scEta_1; scPhi_1=di.scPhi_1;
      scEt_2 = di.scEt_2; scEta_2=di.scEta_2; scPhi_2=di.scPhi_2;
      mass = di.mass;
      pt= di.pt;
      y = di.y;
      phi=di.phi;
    }

    void extractElectron(int index, mithep::TElectron &ele) const {
      const mithep::TDielectron *dielectron = this;
    
    if(index == 1){
      ele. pt                  = dielectron-> pt_1                 ;
      ele. ptUncorr            = dielectron-> ptUncorr_1           ;
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
      ele. scEtUncorr          = dielectron-> scEtUncorr_1         ;
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
    else {
      ele. pt                  = dielectron-> pt_2                 ;
      ele. ptUncorr            = dielectron-> ptUncorr_2           ;
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
      ele. scEtUncorr          = dielectron-> scEtUncorr_2         ;
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
    
    return;
  }

    void assign(const TDielectron *d) {
      mass = d->mass; pt= d->pt; y= d->y; phi= d->phi;
      // 1st ele
      pt_1= d->pt_1; ptUncorr_1= d->ptUncorr_1;
      eta_1= d->eta_1; phi_1= d->phi_1;
      trkIso03_1= d->trkIso03_1; 
      emIso03_1= d->emIso03_1;
      hadIso03_1= d->hadIso03_1;
      chIso_00_01_1= d->chIso_00_01_1;
      chIso_01_02_1= d->chIso_01_02_1;
      chIso_02_03_1= d->chIso_02_03_1;
      chIso_03_04_1= d->chIso_03_04_1;
      chIso_04_05_1= d->chIso_04_05_1;
      gammaIso_00_01_1= d->gammaIso_00_01_1;
      gammaIso_01_02_1= d->gammaIso_01_02_1;
      gammaIso_02_03_1= d->gammaIso_02_03_1;
      gammaIso_03_04_1= d->gammaIso_03_04_1;
      gammaIso_04_05_1= d->gammaIso_04_05_1;
      neuHadIso_00_01_1= d->neuHadIso_00_01_1;
      neuHadIso_01_02_1= d->neuHadIso_01_02_1;
      neuHadIso_02_03_1= d->neuHadIso_02_03_1;
      neuHadIso_03_04_1= d->neuHadIso_03_04_1;
      neuHadIso_04_05_1= d->neuHadIso_04_05_1;
      pfPt_1= d->pfPt_1;
      pfEta_1= d->pfEta_1;
      pfPhi_1= d->pfPhi_1;
      d0_1= d->d0_1;
      dz_1= d->dz_1;
      scEt_1= d->scEt_1;
      scEtUncorr_1= d->scEtUncorr_1;
      scEta_1= d->scEta_1;
      scPhi_1= d->scPhi_1;
      ecalE_1= d->ecalE_1;
      HoverE_1= d->HoverE_1;
      EoverP_1= d->EoverP_1;
      fBrem_1= d->fBrem_1;
      deltaEtaIn_1= d->deltaEtaIn_1;
      deltaPhiIn_1= d->deltaPhiIn_1;
      sigiEtaiEta_1= d->sigiEtaiEta_1;
      partnerDeltaCot_1= d->partnerDeltaCot_1;
      partnerDist_1= d->partnerDist_1;
      mva_1= d->mva_1;
      q_1= d->q_1;
      nExpHitsInner_1= d->nExpHitsInner_1;
      scID_1= d->scID_1;
      trkID_1= d->trkID_1;
      typeBits_1= d->typeBits_1;
      hltMatchBits_1= d->hltMatchBits_1;
      isConv_1= d->isConv_1;
      // 2nd ele
      pt_2= d->pt_2; ptUncorr_2= d->ptUncorr_2;
      eta_2= d->eta_2; phi_2= d->phi_2;
      trkIso03_2= d->trkIso03_2; 
      emIso03_2= d->emIso03_2;
      hadIso03_2= d->hadIso03_2;
      chIso_00_01_2= d->chIso_00_01_2;
      chIso_01_02_2= d->chIso_01_02_2;
      chIso_02_03_2= d->chIso_02_03_2;
      chIso_03_04_2= d->chIso_03_04_2;
      chIso_04_05_2= d->chIso_04_05_2;
      gammaIso_00_01_2= d->gammaIso_00_01_2;
      gammaIso_01_02_2= d->gammaIso_01_02_2;
      gammaIso_02_03_2= d->gammaIso_02_03_2;
      gammaIso_03_04_2= d->gammaIso_03_04_2;
      gammaIso_04_05_2= d->gammaIso_04_05_2;
      neuHadIso_00_01_2= d->neuHadIso_00_01_2;
      neuHadIso_01_02_2= d->neuHadIso_01_02_2;
      neuHadIso_02_03_2= d->neuHadIso_02_03_2;
      neuHadIso_03_04_2= d->neuHadIso_03_04_2;
      neuHadIso_04_05_2= d->neuHadIso_04_05_2;
      pfPt_2= d->pfPt_2;
      pfEta_2= d->pfEta_2;
      pfPhi_2= d->pfPhi_2;
      d0_2= d->d0_2;
      dz_2= d->dz_2;
      scEt_2= d->scEt_2;
      scEtUncorr_2= d->scEtUncorr_2;
      scEta_2= d->scEta_2;
      scPhi_2= d->scPhi_2;
      ecalE_2= d->ecalE_2;
      HoverE_2= d->HoverE_2;
      EoverP_2= d->EoverP_2;
      fBrem_2= d->fBrem_2;
      deltaEtaIn_2= d->deltaEtaIn_2;
      deltaPhiIn_2= d->deltaPhiIn_2;
      sigiEtaiEta_2= d->sigiEtaiEta_2;
      partnerDeltaCot_2= d->partnerDeltaCot_2;
      partnerDist_2= d->partnerDist_2;
      mva_2= d->mva_2;
      q_2= d->q_2;
      nExpHitsInner_2= d->nExpHitsInner_2;
      scID_2= d->scID_2;
      trkID_2= d->trkID_2;
      typeBits_2= d->typeBits_2;
      hltMatchBits_2= d->hltMatchBits_2;
      isConv_2= d->isConv_2;
    }

    ClassDef(TDielectron,3)
  };
}
#endif


// --------------------------------------------------


#endif
