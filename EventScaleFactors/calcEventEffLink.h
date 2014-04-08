#ifndef calcEventEffLink_H
#define calcEventEffLink_H

#include <TROOT.h>
#include <TString.h>
#include <TMatrixD.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include "../Include/DYTools.hh"
#include "../Include/TriggerSelection.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/EventSelector.hh"
#include "../EventScaleFactors/tnpSelectEvents.hh"



//extern TString dirTag;
extern const int nonUniversalHLT; // if nonUniversalHLT=1, HLT_leg1 and HLT_leg2 are used

extern vector<TMatrixD*> dataEff,mcEff;
//extern vector<TMatrixD*> dataEffErrLo,mcEffErrLo;
//extern vector<TMatrixD*> dataEffErrHi,mcEffErrHi;
extern vector<TMatrixD*> dataEffAvgErr,mcEffAvgErr;

extern DYTools::TEtBinSet_t etBinning;
extern int etBinCount;
extern double *etBinLimits;

extern DYTools::TEtaBinSet_t etaBinning;
extern int etaBinCount;
extern double *etaBinLimits;

int initManagers(const TString &confFileName, DYTools::TRunMode_t runMode,
		 InputFileMgr_t &inpMgr, EventSelector_t &evtSelector,
		 TString &puStr, int createDestinationDir,
		 DYTools::TSystematicsStudy_t systMode);

int preparePseudoExps(int nExps, int debug_pseudo_exps);

double getHLTefficiency(DYTools::TDataKind_t dataKind,
			int etBin1, int etaBin1, double et1,
			int etBin2, int etaBin2, double et2);
double getHLTefficiency_IK(DYTools::TDataKind_t dataKind,
			   int etBin1, int etaBin1, double et1,
			   int etBin2, int etaBin2, double et2);
double getHLTefficiencyErr(DYTools::TDataKind_t dataKind,
			   int etBin1, int etaBin1, double et1,
			   int etBin2, int etaBin2, double et2);

int fillEfficiencyConstants( const InputFileMgr_t &inpMgr,
		  DYTools::TSystematicsStudy_t systMode = DYTools::NO_SYST );


void drawEfficiencies(TFile *fRootOutput, 
		      int ignoreAsymHLT=0,
		      std::vector<TGraphErrors*> *grDataVec=NULL,
		      std::vector<TGraphErrors*> *grMCVec=NULL);

// access scale factors
double findEventScaleFactor(int kind, const esfSelectEvent_t &data); // -1 for total
double findEventScaleFactorSmeared(int kind, const esfSelectEvent_t &data,
				   int iexp); // -1 for total

//
// Original functions
//

TGraphAsymmErrors* getAsymGraph_vsEt(DYTools::TEtBinSet_t etBinning_inp,
				     DYTools::TEtaBinSet_t etaBinning_inp,
				     int iEta,
				     const TMatrixD &Meff,
				     const TMatrixD &MeffLo,
				     const TMatrixD &MeffHi,
				     TH1D  **histo=NULL,
				     const char *histoName=NULL);



TGraphAsymmErrors* getAsymGraph_vsEta(DYTools::TEtBinSet_t etBinning_inp,
				      DYTools::TEtaBinSet_t etaBinning_inp,
				      int iEt,
				      const TMatrixD &Meff,
				      const TMatrixD &MeffLo,
				      const TMatrixD &MeffHi,
				      TH1D  **histo=NULL,
				      const char *histoName=NULL);

inline
TGraphAsymmErrors* getAsymGraph(int vsEt,
				DYTools::TEtBinSet_t etBinning_inp,
				DYTools::TEtaBinSet_t etaBinning_inp,
				int iBin,
				const TMatrixD &Meff,
				const TMatrixD &MeffLo,
				const TMatrixD &MeffHi,
				TH1D  **histo=NULL,
				const char *histoName=NULL) {
  return (vsEt) ? 
    getAsymGraph_vsEt(etBinning_inp,etaBinning_inp,iBin,Meff,MeffLo,MeffHi,histo,histoName) :
    getAsymGraph_vsEta(etBinning_inp,etaBinning_inp,iBin,Meff,MeffLo,MeffHi,histo,histoName);
}

// ---------------------------------------

// load efficiency as matrix
// if the method is count-count, load it as weighted
int loadEff(const TString &fname, int weighted, TMatrixD **eff, TMatrixD **effLo, TMatrixD **effHi);

// ---------------------------------------

#endif
