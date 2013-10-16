#ifndef calcEventEffLink_H
#define calcEventEffLink_H

#include <TROOT.h>
#include <TString.h>
#include <TMatrixD.h>
#include <TGraphErrors.h>
#include "../Include/DYTools.hh"
#include "../Include/TriggerSelection.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/EventSelector.hh"



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
		 TString &puStr, int createDestinationDir);

double getHLTefficiency(DYTools::TDataKind_t dataKind,
			int etBin1, int etaBin1, double et1,
			int etBin2, int etaBin2, double et2);
double getHLTefficiency_IK(DYTools::TDataKind_t dataKind,
			   int etBin1, int etaBin1, double et1,
			   int etBin2, int etaBin2, double et2);
double getHLTefficiencyErr(DYTools::TDataKind_t dataKind,
			   int etBin1, int etaBin1, double et1,
			   int etBin2, int etaBin2, double et2);

int fillEfficiencyConstants( const InputFileMgr_t &inpMgr );


void drawEfficiencies(TFile *fRootOutput, 
		      int ignoreAsymHLT=0,
		      std::vector<TGraphErrors*> *grDataVec=NULL,
		      std::vector<TGraphErrors*> *grMCVec=NULL);





#endif
