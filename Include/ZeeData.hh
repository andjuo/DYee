#ifndef ZEE_DATA_HH
#define ZEE_DATA_HH

#include <TObject.h>
#include "../Include/DYTools.hh"

#include "../Include/TEventInfo.hh"
#include "../Include/TDielectron.hh"
#include "../Include/EventWeight.hh"

#ifdef DYee8TeV_reg
#include <TLorentzVector.h>
#include "../Include/TPhoton.hh"
#endif

// define what extra fields to store

#define ZeeData_storeGoldenFlag 

#ifdef DYee8TeV_reg
# define ZeeData_storeUnregEn
#endif



class ZeeData_t : public TObject
{
public:
  UInt_t  runNum;                          // run number in data
  UInt_t  evtNum;                          // event number in data
  UInt_t  lumiSec;                         // lumi section      
  //UInt_t  nTracks0;                        // number of reconstructed tracks in event
  //UInt_t  nCaloTowers0;                    // number of reconstructed calorimeter towers in event
  UInt_t  nPV;                             // number of valid reconstructed primary vertices in event                                          
  //UInt_t  nGoodPV;                         // number of good PVs
  //UInt_t  nJets;                           // number of jets (with some requirements)
  //Float_t caloMEx, caloMEy, caloSumET;     // calorimeter MET 
  //Float_t tcMEx, tcMEy, tcSumET;           // track-corrected MET
  //Float_t pfMEx, pfMEy, pfSumET;           // particle flow MET

  Float_t mass, pt, y, phi;                // dielectron kinematics , y is |y|
  Float_t y_signed;        // y_signed includes the sign
  
  Float_t pt_1, eta_1, phi_1;              // leading electron
  Float_t scEt_1, scEta_1, scPhi_1;  
  ULong_t  hltMatchBits_1;
  Int_t   q_1;
  
  Float_t pt_2, eta_2, phi_2;              // lagging electron
  Float_t scEt_2, scEta_2, scPhi_2;
  ULong_t  hltMatchBits_2;
  Int_t   q_2;

  Double_t weight;                          // total event weight  
  Double_t gen_weight;                          // generated weight  
  Double_t pu_weight;                 // pu weight
  Double_t fewz_weight;              // fewz_weight
#ifdef ZeeData_storeGoldenFlag
  Int_t golden_1, golden_2;    // whether electron R9>0.94 (golden electron), or <0.94 (showering electron)
  Float_t R9_1, R9_2;
#endif
#ifdef ZeeData_storeUnregEn
  Float_t ptUncorr_1, ptUncorr_2;
  Float_t scEtUncorr_1, scEtUncorr_2;
#endif


//"runNum/i:evtNum:lumiSec:nTracks0:nCaloTowers0:nPV:nGoodPV:nJets:caloMEx/F:caloMEy:caloSumET:tcMEx:tcMEy:tcSumET:pfMEx:pfMEy:pfSumET:mass:pt:y:phi:pt_1:eta_1:phi_1:scEt_1:scEta_1:scPhi_1:hltMatchBits_1/i:q_1/I:pt_2/F:eta_2:phi_2:scEt_2:scEta_2:scPhi_2:hltMatchBits_2/i:q_2/I:weight/F"
//"runNum/i:evtNum:lumiSec:nTracks0:nCaloTowers0:nPV:nJets:caloMEx/F:caloMEy:caloSumET:tcMEx:tcMEy:tcSumET:pfMEx:pfMEy:pfSumET:mass:pt:y:phi:pt_1:eta_1:phi_1:scEt_1:scEta_1:scPhi_1:hltMatchBits_1/i:q_1/I:pt_2/F:eta_2:phi_2:scEt_2:scEta_2:scPhi_2:hltMatchBits_2/i:q_2/I:weight/F"

  void assign(const ZeeData_t &a) {
    runNum=a.runNum;
    evtNum=a.evtNum;
    lumiSec=a.lumiSec;
    //nTracks0=a.nTracks0;
    //nCaloTowers0=a.nCaloTowers0;
    nPV=a.nPV;
    //nGoodPV=a.nGoodPV;
    //nJets=a.nJets;
    //caloMEx=a.caloMEx;   caloMEy=a.caloMEy;   caloSumET=a.caloSumET;
    //tcMEx=a.tcMEx;   tcMEy=a.tcMEy;   tcSumET=a.tcSumET;
    //pfMEx=a.tcMEx;   pfMEy=a.pfMEy;   pfSumET=a.pfSumET;
    mass=a.mass;
    pt=a.pt;
    y=a.y;
    phi=a.phi;
    y_signed=a.y_signed;

    pt_1=a.pt_1; eta_1=a.eta_1; phi_1=a.phi_1;
    scEt_1=a.scEt_1; scEta_1=a.scEta_1; scPhi_1=a.scPhi_1;
    hltMatchBits_1=a.hltMatchBits_1;
    q_1=a.q_1;
    pt_2=a.pt_2; eta_2=a.eta_2; phi_2=a.phi_2;
    scEt_2=a.scEt_2; scEta_2=a.scEta_2; scPhi_2=a.scPhi_2;
    hltMatchBits_2=a.hltMatchBits_2;
    q_2=a.q_2;
    weight=a.weight;
    gen_weight=a.gen_weight;
    pu_weight=a.pu_weight;
    fewz_weight=a.fewz_weight;
#ifdef ZeeData_storeGoldenFlag
    golden_1=a.golden_1; golden_2=a.golden_2;
    R9_1=a.R9_1; R9_2=a.R9_2;
#endif
#ifdef ZeeData_storeUnregEn
    ptUncorr_1=a.ptUncorr_1;
    ptUncorr_2=a.ptUncorr_2;
    scEtUncorr_1=a.scEtUncorr_1;
    scEtUncorr_2=a.scEtUncorr_2;
#endif
  }


  void Assign(const mithep::TEventInfo *info, const mithep::TDielectron *dielectron, 
		const UInt_t npv, 
		//const UInt_t njets, 
		const Double_t set_gen_weight, 
		const Double_t set_pu_weight,
		const Double_t set_fewz_weight
#ifdef ZeeData_storeGoldenFlag
	      , const Float_t setR9val_1, const Float_t setR9val_2
#endif
		) {
    ZeeData_t *data= this;
  data->runNum         = info->runNum;
  data->evtNum         = info->evtNum;
  data->lumiSec        = info->lumiSec;
  //data->nTracks0       = 0;
  //data->nCaloTowers0   = 0;
  data->nPV            = npv;
  //data->nJets          = njets;
  //data->caloMEx        = 0;
  //data->caloMEy        = 0;
  //data->caloSumET      = 0;
  //data->tcMEx          = 0; //info->trkMET * cos(info->trkMETphi);
  //data->tcMEy          = 0; //info->trkMET * sin(info->trkMETphi);
  //data->tcSumET        = info->trkSumET;
  //data->pfMEx          = 0; //info->pfMET * cos(info->pfMETphi);
  //data->pfMEy          = 0; //info->pfMET * sin(info->pfMETphi);
  //data->pfSumET        = info->pfSumET;
  data->mass           = dielectron->mass;
  data->pt             = dielectron->pt;
  data->y              = fabs(dielectron->y);
  data->y_signed       = dielectron->y;
  data->phi            = dielectron->phi; 
  // no endorsement that the 1st electron is the leading one
  data->pt_1           = dielectron->pt_1;
  data->eta_1          = dielectron->eta_1;
  data->phi_1          = dielectron->phi_1;
  data->scEt_1         = dielectron->scEt_1;
  data->scEta_1        = dielectron->scEta_1;
  data->scPhi_1        = dielectron->scPhi_1;
  data->hltMatchBits_1 = dielectron->hltMatchBits_1;
  data->q_1            = dielectron->q_1;
  data->pt_2           = dielectron->pt_2;
  data->eta_2          = dielectron->eta_2;
  data->phi_2          = dielectron->phi_2;
  data->scEt_2         = dielectron->scEt_2;
  data->scEta_2        = dielectron->scEta_2;
  data->scPhi_2        = dielectron->scPhi_2;
  data->hltMatchBits_2 = dielectron->hltMatchBits_2;
  data->q_2            = dielectron->q_2;
  data->weight         = set_gen_weight * set_pu_weight * set_fewz_weight;
  data->gen_weight     = set_gen_weight;
  data->pu_weight      = set_pu_weight;
  data->fewz_weight    = set_fewz_weight;
#ifdef ZeeData_storeGoldenFlag
  data->golden_1 = (setR9val_1 > 0.94) ? 1:0;
  if (setR9val_1<0) data->golden_1=-1;
  data->golden_2 = (setR9val_2 > 0.94) ? 1:0;
  if (setR9val_2<0) data->golden_2=-1;
  R9_1=setR9val_1;
  R9_2=setR9val_2;
#endif
#ifdef ZeeData_storeUnregEn
  ptUncorr_1= dielectron->ptUncorr_1;
  ptUncorr_2= dielectron->ptUncorr_2;
  scEtUncorr_1= dielectron->scEtUncorr_1;
  scEtUncorr_2= dielectron->scEtUncorr_2;
#endif
  }

  // ----------------------

  void Assign(const mithep::TEventInfo *info, 
	      const mithep::TDielectron *dielectron, 
#ifdef ZeeData_storeGoldenFlag
	      const mithep::TPhoton *sc_1,
	      const mithep::TPhoton *sc_2,
#endif
	      const UInt_t npv, 
	      //const UInt_t njets, 
	      const EventWeight_t &ew
	      ) {

#ifdef ZeeData_storeGoldenFlag

    int scIDok_1=((sc_1!=NULL) && (sc_1->scID == dielectron->scID_1)) ? 1:0;
    int scIDok_2=((sc_2!=NULL) && (sc_2->scID == dielectron->scID_2)) ? 1:0;
    Float_t setR9_1= (scIDok_1) ? sc_1->R9 : -1.;
    Float_t setR9_2= (scIDok_2) ? sc_2->R9 : -1.;
    Assign(info,dielectron,
	   npv,
	   ew.baseWeight(),
	   ew.puWeight(),
	   ew.fewzWeight(),
	   setR9_1,
	   setR9_2
	   );

#else

    Assign(info,dielectron,
	   npv,
	   ew.baseWeight(),
	   ew.puWeight(),
	   ew.fewzWeight()
	   );

#endif
  }

  // ----------------------

#ifdef ZeeData_storeUnregEn
  void replace2UncorrEn(int check=0, double scale=static_cast<double>(1.0)) {
    if (!check) {
      pt_1  =scale*(ptUncorr_1 - pt_1) + pt_1;
      scEt_1=scale*(scEtUncorr_1 - scEt_1 ) + scEt_1;
      pt_2  =scale*(ptUncorr_2 - pt_2) + pt_2;
      scEt_2=scale*(scEtUncorr_2 - scEt_2 ) + scEt_2;
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
#endif

  // ----------------------


#ifdef ZeeData_storeGoldenFlag
  ClassDef(ZeeData_t,3)
#else
  ClassDef(ZeeData_t,2)
#endif

};


// -----------------------------------------------------


#endif
