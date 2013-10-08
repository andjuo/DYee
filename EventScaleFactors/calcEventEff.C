#include <TFile.h>
#include <TChain.h>
#include <TString.h>
#include <TBenchmark.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
//#include <TH1D.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TRandom.h>
#include <TVectorD.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TTimeStamp.h>
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <vector>                   // STL vector class
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

//#define debug_systScaleArrs

#include <vector>

#include "../Include/CPlot.hh"
#include "../Include/MitStyleRemix.hh"

#include "../Include/DYTools.hh"
//#include "../Include/MyTools.hh"
//#include "../Include/EleIDCuts.hh"

//#include "../Include/EWKAnaDefs.hh"
//#include "../Include/TGenInfo.hh"
//#include "../Include/TEventInfo.hh"
//#include "../Include/TDielectron.hh"
//#include "../Include/TElectron.hh"
//#include "../Include/TVertex.hh"
//#include "../Include/TriggerSelection.hh"
#include "../EventScaleFactors/cutFunctions.hh"
#include "../EventScaleFactors/tnpSelectEvents.hh"

#include "../Include/AccessOrigNtuples.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/EventSelector.hh"
#include "../Include/EventWeight.hh"
//#include "../Include/PUReweight.hh"
//#include "../Include/UnfoldingTools.hh"
//#include "../Include/FEWZ.hh"

using namespace mithep;
using namespace std;

// nonUniversalHLT=1 causes to use the double-trigger efficiencies
// EffHLT = eff(1,trig1)*eff(2,trig2) + eff(1,trig2)*eff(2,trig1) - eff(1,trig1)*eff(2,trig1)
// where the eff(i,trig1) is the efficiency i'th electron in HLT_leg1,
// eff(i,trig2) is the efficiency of i'th electron in HLT_leg2,
// the threshold for HLT_leg1 is higher than for HLT_leg2
const int nonUniversalHLT=1;// if nonUniversalHLT=1, HLT_leg1 and HLT_leg2 are used
const int NEffTypes=3 + 2*nonUniversalHLT;
const int allowToIgnoreAnalysisTag=1; // efficiencies and selected events...
// .... are the same for 1D and 2D

// Declaration of arrays for systematic studies
typedef double EffArray_t[NEffTypes][DYTools::nEtBinsMax][DYTools::nEtaBinsMax]; // largest storage

const int nexp=100;
typedef TH1D* SystTH1DArray_t[DYTools::nMassBins][nexp];   // mass index
typedef TH1D* SystTH1DArrayFI_t[DYTools::nUnfoldingBinsMax][nexp]; // flat(mass,y) index

template<class T> T SQR(const T& x) { return x*x; }


//=== FUNCTION DECLARATIONS ======================================================================================

int fillEfficiencyConstants( const InputFileMgr_t &inpMgr );

int fillOneEfficiency(const TString &dirTag, const TString filename, 
  UInt_t kindIdx, vector<TMatrixD*> &data, vector<TMatrixD*> &dataErrLo, 
  vector<TMatrixD*> &dataErrHi, vector<TMatrixD*> &dataAvgErr, int weightedCnC);

// moved to AccessOrigNtuples.hh
//Bool_t matchedToGeneratorLevel(const TGenInfo *gen, 
//  const TDielectron *dielectron);

int createSelectionFile(const InputFileMgr_t &mcMgr, 
			EventSelector_t &evtSelector,
			const TString &outSkimFName, DYTools::TRunMode_t runMode);

double findEventScaleFactor(int kind, const esfSelectEvent_t &data); // -1 for total

double findEventScaleFactor(int kind, 
			    double Et1, double eta1,
			    double Et2, double eta2);

double findScaleFactor(DYTools::TEfficiencyKind_t kind, int etBin, int etaBin);

// nonUniversalHLT scale factor
double findScaleFactorHLT(int etBin1, int etaBin1, double et1,
			  int etBin2, int etaBin2, double et2);
// nonUniversalHLT (double-trigger) efficiency
double getHLTefficiency(DYTools::TDataKind_t dataKind,
			int etBin1, int etaBin1, double et1,
			int etBin2, int etaBin2, double et2);
double getHLTefficiencyErr(DYTools::TDataKind_t dataKind,
			   int etBin1, int etaBin1, double et1,
			   int etBin2, int etaBin2, double et2);

double findEventScaleFactorSmeared(int kind, const esfSelectEvent_t &data, 
				   int iexp); // -1 for total

double findEventScaleFactorSmeared(int kind, 
				   double Et1, double eta1,
				   double Et2, double eta2,
				   int iexp);

double findScaleFactorSmeared(
       DYTools::TEfficiencyKind_t kind, 
       int scEtBin, int scEtaBin, 
       int iexp);

double findScaleFactorSmeared(
       DYTools::TEfficiencyKind_t kind, 
       int scEtBin, int scEtaBin, 
       const EffArray_t &dataRndWeight, const EffArray_t &mcRndWeight);

// nonUniversalHLT scale factor
double findScaleFactorHLTSmeared(int scEtBin1, int scEtaBin1, double scEt1,
				 int scEtBin2, int scEtaBin2, double scEt2,
				 int iexp);

// nonUniversalHLT scale factor
double findScaleFactorHLTSmeared(int scEtBin1, int scEtaBin1, double scEt1,
				 int scEtBin2, int scEtaBin2, double scEt2,
				 const EffArray_t &dataRndWeight, 
				 const EffArray_t &mcRndWeight);

// nonUniversalHLT (double-trigger) efficiency
double getHLTefficiencySmeared(DYTools::TDataKind_t dataKind,
			       int etBin1, int etaBin1, double Et1,
			       int etBin2, int etaBin2, double Et2,
			       const EffArray_t &rndWeight);

void drawEfficiencies(TFile *fRoot);

void drawEfficiencyGraphs(
     TGraphErrors *grData, TGraphErrors *grMc,
     TString yAxisTitle, TString text, TString plotName,
     TFile *fRoot);

//template<class Graph_t>
//void drawEfficiencyGraphsAsymmErrors(Graph_t *grData, Graph_t *grMc,
//			  TString yAxisTitle, TString text, TString plotName);

void drawScaleFactors(TFile *fRoot);

void drawScaleFactorGraphs(TGraphErrors *gr, TString yAxisTitle, TString text, 
			   TString plotName, TFile *fRoot);

void drawEventScaleFactors(
     const TVectorD &scaleRecoV, const TVectorD &scaleRecoErrV,
     const TVectorD &scaleIdV , const TVectorD &scaleIdErrV ,
     const TVectorD &scaleHltV, const TVectorD &scaleHltErrV,
     const TVectorD &scaleV   , const TVectorD &scaleErrV   ,
     TFile *fRoot);

// nonUniversalHLT
void drawEventScaleFactorsHLT(
     const TVectorD &scaleHLTleg1V, const TVectorD &scaleHLTleg1ErrV,
     const TVectorD &scaleHLTleg2V, const TVectorD &scaleHLTleg2ErrV ,
     TFile *fRoot);

// nonUniversalHLT
void drawEventScaleFactorsFI(
     const TVectorD &scaleRecoFIV, const TVectorD &scaleRecoErrFIV,
     const TVectorD &scaleIdFIV , const TVectorD &scaleIdErrFIV ,
     const TVectorD &scaleHltFIV, const TVectorD &scaleHltErrFIV,
     const TVectorD &scaleV   , const TVectorD &scaleErrFIV   ,
     int rapidityBinIndex0,
     TFile *fRoot,
     std::vector<CPlot*> *cplotV=NULL);

void drawEventScaleFactorsHltFI(
     const TVectorD &scaleHltLeg1FIV, const TVectorD &scaleHltLeg1ErrFIV,
     const TVectorD &scaleHltLeg2FIV , const TVectorD &scaleHltLeg2ErrFIV ,
     int rapidityBinIndex0,
     TFile *fRoot,
     std::vector<CPlot*> *cplotV=NULL);

void drawEventScaleFactorGraphs(
     TGraphErrors *gr, TString yAxisTitle, 
     TString plotName, TFile *fRoot);

double errOnRatio(double a, double da, double b, double db);


// Global variables
//const int nexp = 100;

template<class TSystTH1DArray_t>
void deriveScaleMeanAndErr(const int binCount, const int nexpCount, 
			   TSystTH1DArray_t &systScale, 
			   TVectorD &scaleMean, TVectorD &scaleMeanErr) {

  if ((scaleMean.GetNoElements() != binCount) ||
      (scaleMeanErr.GetNoElements() != binCount)) {
    std::cout << "Error in derive: binCount=" << binCount 
	      << ", scaleMean[" << scaleMean.GetNoElements() 
	      << "], scaleMeanErr[" << scaleMeanErr.GetNoElements() << "\n";
    assert(0);
  }
  for(int ibin = 0; ibin < binCount; ibin++){
    scaleMean[ibin] = 0;
    scaleMeanErr[ibin] = 0;
    for(int iexp = 0; iexp < nexpCount; iexp++){
      scaleMean[ibin] += systScale[ibin][iexp]->GetMean();
      scaleMeanErr[ibin] += SQR(systScale[ibin][iexp]->GetMean());
    }
    scaleMean[ibin] = scaleMean[ibin]/double(nexpCount);
    scaleMeanErr[ibin] = sqrt( scaleMeanErr[ibin] / double(nexpCount) 
			       - SQR(scaleMean[ibin]) );
    //std::cout << "scaleMean[ibin=" << ibin << "]=" << scaleMean[ibin] << ", err=" << scaleMeanErr[ibin] << "\n";
  }
  return;
}


//=== Constants ==========================

const bool savePlots = true;
const bool correlationStudy = false;

// File names for efficiency measurements from tag and probe
//TString          dirTag;

vector<TMatrixD*> dataEff,mcEff;
vector<TMatrixD*> dataEffErrLo,mcEffErrLo;
vector<TMatrixD*> dataEffErrHi,mcEffErrHi;
vector<TMatrixD*> dataEffAvgErr,mcEffAvgErr;

// Global variables
//const int nexp = 100;

EffArray_t *ro_Data=NULL;
EffArray_t *ro_MC=NULL;


DYTools::TEtBinSet_t etBinning= DYTools::ETBINS_UNDEFINED;
int etBinCount=0;
double *etBinLimits=NULL;

DYTools::TEtaBinSet_t etaBinning= DYTools::ETABINS_UNDEFINED;
int etaBinCount=0;
double *etaBinLimits=NULL;

// ---------- Setup main variables ---------------

int initManagers(const TString &confFileName, DYTools::TRunMode_t runMode,
		 InputFileMgr_t &inpMgr, EventSelector_t **evtSelector,
		 TString &puStr, int createDestinationDir) {
  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;

  //InputFileMgr_t inpMgr;
  if (!inpMgr.LoadTnP(confFileName)) return retCodeError;
  // no energy correction for this evaluation
  inpMgr.clearEnergyScaleTag();

  const int puReweight=inpMgr.puReweightFlag();
  const int useFewzWeights=inpMgr.fewzFlag();
  //TString 
  puStr = (puReweight) ? "PU" : "";
  if (!useFewzWeights) puStr.Append("_noFEWZ");
  if (nonUniversalHLT) puStr.Append("_HLTlegs");
  //CPlot::sOutDir = TString("plots") + DYTools::analysisTag + puStr;
  //std::cout << "changed CPlot::sOutDir=<" << CPlot::sOutDir << ">\n";
  //gSystem->mkdir(CPlot::sOutDir,true);

  // Construct eventSelector, update mgr and plot directory
  TString plotsExtraTag=puStr;
  //EventSelector_t 
  *evtSelector = new EventSelector_t(inpMgr,runMode,systMode,
			  "",plotsExtraTag, EventSelector::_selectOrdered);
  if (!(*evtSelector)) {
    std::cout << "failed to create EventSelector_t\n";
    return retCodeError;
  }
  (*evtSelector)->setTriggerActsOnData(false);

  // Prepare output directory
  inpMgr.constDir(systMode,createDestinationDir);

  //
  // prepare global variables
  // 

  ro_Data= new EffArray_t[nexp];
  ro_MC  = new EffArray_t[nexp];
  assert(ro_Data && ro_MC);


  //dirTag=inpMgr.tnpTag();

  etBinning=inpMgr.etBinsKind();
  etBinCount=DYTools::getNEtBins(etBinning);
  etBinLimits=DYTools::getEtBinLimits(etBinning);
  
  etaBinning=inpMgr.etaBinsKind();
  etaBinCount=DYTools::getNEtaBins(etaBinning);
  etaBinLimits=DYTools::getEtaBinLimits(etaBinning);

  if (etBinCount>DYTools::nEtBinsMax) {
    printf("ERROR in the code: etBinCount=%d, hard-coded DYTools::nEtBinsMax=%d\n",etBinCount,DYTools::nEtBinsMax);
    std::cout << std::endl;
    assert(0);
  }
  if (etaBinCount>DYTools::nEtaBinsMax) {
    printf("ERROR in the code: etaBinCount=%d, hard-coded DYTools::nEtaBinsMax=%d\n",etaBinCount,DYTools::nEtaBinsMax);
    std::cout << std::endl;
    assert(0);
  }

  return retCodeOk;
}


//=== MAIN MACRO =================================================================================================

int calcEventEff(const TString confFileName,
		  int selectEvents, 
		  DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN)
{

  gBenchmark->Start("calcEventEff");

//  ---------------------------------
//       Preliminary checks
//  ---------------------------------

  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;

  // verify whether it was a compilation check
  if (confFileName.Contains("_DebugRun_")) {
    std::cout << "calcEventEff: _DebugRun_ detected. Terminating the script\n";
    return retCodeOk;
  }

  {
    DYTools::printExecMode(runMode,systMode);
    const int debug_print=1;
    if (!DYTools::checkSystMode(systMode,debug_print,1, DYTools::NO_SYST)) 
      return retCodeError;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  InputFileMgr_t inpMgr;
  if (!inpMgr.LoadTnP(confFileName)) return retCodeError;
  // no energy correction for this evaluation
  inpMgr.clearEnergyScaleTag();

  const int puReweight=inpMgr.puReweightFlag();
  const int useFewzWeights=inpMgr.fewzFlag();
  TString puStr = (puReweight) ? "PU" : "";
  if (!useFewzWeights) puStr.Append("_noFEWZ");
  if (nonUniversalHLT) puStr.Append("_HLTlegs");
  //CPlot::sOutDir = TString("plots") + DYTools::analysisTag + puStr;
  //std::cout << "changed CPlot::sOutDir=<" << CPlot::sOutDir << ">\n";
  //gSystem->mkdir(CPlot::sOutDir,true);

  // Construct eventSelector, update mgr and plot directory
  TString plotsExtraTag=puStr;
  EventSelector_t evtSelector(inpMgr,runMode,systMode,
			      "",plotsExtraTag, EventSelector::_selectOrdered);
  evtSelector.setTriggerActsOnData(false);

  // Prepare output directory
  inpMgr.constDir(systMode,1);

//  ---------------------------------
//         Normal execution
//  ---------------------------------

  ro_Data= new EffArray_t[nexp];
  ro_MC  = new EffArray_t[nexp];
  assert(ro_Data && ro_MC);


  //dirTag=inpMgr.tnpTag();

  etBinning=inpMgr.etBinsKind();
  etBinCount=DYTools::getNEtBins(etBinning);
  etBinLimits=DYTools::getEtBinLimits(etBinning);
  
  etaBinning=inpMgr.etaBinsKind();
  etaBinCount=DYTools::getNEtaBins(etaBinning);
  etaBinLimits=DYTools::getEtaBinLimits(etaBinning);

  if (etBinCount>DYTools::nEtBinsMax) {
    printf("ERROR in the code: etBinCount=%d, hard-coded DYTools::nEtBinsMax=%d\n",etBinCount,DYTools::nEtBinsMax);
    std::cout << std::endl;
    assert(0);
  }
  if (etaBinCount>DYTools::nEtaBinsMax) {
    printf("ERROR in the code: etaBinCount=%d, hard-coded DYTools::nEtaBinsMax=%d\n",etaBinCount,DYTools::nEtaBinsMax);
    std::cout << std::endl;
    assert(0);
  }


  TString selectEventsFName=inpMgr.tnpSelEventsFullName(systMode,0);
  
  // remove HLTlegs tag from the file name, if needed
  if (nonUniversalHLT) selectEventsFName.ReplaceAll("_HLTlegs","");
  //selectEventsFName.ReplaceAll(".root","-chk.root");
  std::cout << "selectEventsFName=<" << selectEventsFName << ">\n";
  
  if (selectEvents) {
    if (!createSelectionFile(inpMgr, evtSelector,selectEventsFName, runMode)) {
      std::cout << "failed to create selection file <" 
		<< selectEventsFName << ">\n";
    }
    else {
      std::cout << "selection file <" << selectEventsFName << "> created\n";
      gBenchmark->Show("calcEventEff");
    }
    return retCodeError;
  }

  int applyNtupleExtraTag=1;
  TString corrName="scale_factors";
  if (nonUniversalHLT) corrName.Append("_asymHLT");
  TString sfConstFileName=inpMgr.correctionFullFileName(corrName,systMode,applyNtupleExtraTag);
  std::cout << "\n\tsfConstFileName=<" << sfConstFileName << ">\n";
  if (nonUniversalHLT) std::cout << "\t\tnote: HLT legs are considered!\n\n";


  // Read efficiency constants from ROOT files
  // This has to be done AFTER configuration file is parsed
  if (!fillEfficiencyConstants( inpMgr ) ) {
    return retCodeError;
  }

  //std::cout << "stopping\n";
  //return retCodeStop;


  TH1D *hScale = new TH1D("hScale", "", 150, 0.0, 1.5);
  TH1D *hScaleReco = new TH1D("hScaleReco", "", 150, 0.0, 1.5);
  TH1D *hScaleId  = new TH1D("hScaleId" , "", 150, 0.0, 1.5);
  TH1D *hScaleHlt = new TH1D("hScaleHlt", "", 150, 0.0, 1.5);
  TH1D *hScaleHltLeg1 = (nonUniversalHLT) ? new TH1D("hScaleHltLeg1", "", 150, 0.0, 1.5) : NULL;
  TH1D *hScaleHltLeg2 = (nonUniversalHLT) ? new TH1D("hScaleHltLeg2", "", 150, 0.0, 1.5) : NULL;
  std::vector<TH1D*> hScaleEffV;
  hScaleEffV.reserve(NEffTypes);
  hScaleEffV.push_back(hScaleReco);
  hScaleEffV.push_back(hScaleId);
  hScaleEffV.push_back(hScaleHlt);
  if (nonUniversalHLT) {
    hScaleEffV.push_back(hScaleHltLeg1);
    hScaleEffV.push_back(hScaleHltLeg2);
  }

  TH1D *hZpeakEt = new TH1D("hZpeakEt", "", etBinCount, etBinLimits);
  vector<TH1D*> hLeadingEtV;
  vector<TH1D*> hTrailingEtV;
  vector<TH1D*> hElectronEtV;

  vector<TH1D*> hScaleV;     // mass indexing
  vector<TH1D*> hScaleRecoV;
  vector<TH1D*> hScaleIdV;
  vector<TH1D*> hScaleHltV;
  vector<TH1D*> hScaleHltLeg1V , hScaleHltLeg2V;
  vector<TH1D*> hScaleFIV; // flat indexing
  vector<TH1D*> hScaleRecoFIV;
  vector<TH1D*> hScaleIdFIV;
  vector<TH1D*> hScaleHltFIV;
  vector<TH1D*> hScaleHltLeg1FIV, hScaleHltLeg2FIV;

  hLeadingEtV.reserve(DYTools::nMassBins);
  hTrailingEtV.reserve(DYTools::nMassBins);
  hElectronEtV.reserve(DYTools::nMassBins);

  hScaleV.reserve(DYTools::nMassBins);
  hScaleRecoV.reserve(DYTools::nMassBins);
  hScaleIdV.reserve(DYTools::nMassBins);
  hScaleHltV.reserve(DYTools::nMassBins);
  if (nonUniversalHLT) {
    hScaleHltLeg1V.reserve(DYTools::nMassBins);
    hScaleHltLeg2V.reserve(DYTools::nMassBins);
  }

  int nUnfoldingBins = DYTools::getTotalNumberOfBins();
  hScaleFIV.reserve(nUnfoldingBins);
  hScaleRecoFIV.reserve(nUnfoldingBins);
  hScaleIdFIV.reserve(nUnfoldingBins);
  hScaleHltFIV.reserve(nUnfoldingBins);
  if (nonUniversalHLT) {
    hScaleHltLeg1FIV.reserve(nUnfoldingBins);
    hScaleHltLeg2FIV.reserve(nUnfoldingBins);
  }

  for(int i=0; i<DYTools::nMassBins; i++){
    double *rapidityBinLimits=DYTools::getYBinArray(i);

    TString base = TString("hScaleV_massBin");
    base += i;
    hScaleV.push_back(new TH1D(base,base,150,0.0,1.5));
    hScaleRecoV.push_back(new TH1D(base+TString("_reco"),
				   base+TString("_reco"),150,0.0,1.5));
    hScaleIdV.push_back(new TH1D(base+TString("_id" ),
				 base+TString("_id" ),150,0.0,1.5));
    hScaleHltV.push_back(new TH1D(base+TString("_hlt" ),
				  base+TString("_hlt" ),150,0.0,1.5));
    if (nonUniversalHLT) {
      hScaleHltLeg1V.push_back(new TH1D(base+TString("_hltLeg1"),
					base+TString("_hltLeg1"),150,0.0,1.5));
      hScaleHltLeg2V.push_back(new TH1D(base+TString("_hltLeg2"),
					base+TString("_hltLeg2"),150,0.0,1.5));
    }
    base = "hLeadingEt_";
    base += i;
    hLeadingEtV.push_back(new TH1D(base,base,etBinCount, etBinLimits));
    base = "hTrailingEt_";
    base += i;
    hTrailingEtV.push_back(new TH1D(base,base,etBinCount, etBinLimits));
    base = "hElectronEt_";
    base += i;
    hElectronEtV.push_back(new TH1D(base,base,etBinCount, etBinLimits));
    
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
      char buf[50];
      sprintf(buf,"mIdx%d_y%4.2lf-%4.2lf",
	      i,rapidityBinLimits[yi],rapidityBinLimits[yi+1]);
      TString idxStr=buf;
      idxStr.ReplaceAll(".","_");
      base = TString("hScaleV_")+idxStr;
      hScaleFIV.push_back(new TH1D(base,base,150,0.0,1.5));
      hScaleRecoFIV.push_back(new TH1D(base+TString("_reco"),
				       base+TString("_reco"),150,0.0,1.5));
      hScaleIdFIV.push_back(new TH1D(base+TString("_id" ),
				     base+TString("_id" ),150,0.0,1.5));
      hScaleHltFIV.push_back(new TH1D(base+TString("_hlt"),
				      base+TString("_hlt"),150,0.0,1.5));
      if (nonUniversalHLT) {
	hScaleHltLeg1FIV.push_back(new TH1D(base+TString("_hltLeg1"),
					    base+TString("_hltLeg1"),150,0.0,1.5));
	hScaleHltLeg2FIV.push_back(new TH1D(base+TString("_hltLeg2"),
					    base+TString("_hltLeg2"),150,0.0,1.5));
      }
    }
    delete rapidityBinLimits;
  }

  // Create Gaussian-distributed random offsets for each pseudo-experiment
  /*
  for(int i=0; i<nexp; i++){
    ro_D_B_reco[i] = gRandom->Gaus(0.0,1.0);
    ro_D_E_reco[i] = gRandom->Gaus(0.0,1.0);
    ro_M_B_reco[i] = gRandom->Gaus(0.0,1.0);
    ro_M_E_reco[i] = gRandom->Gaus(0.0,1.0);
    
    ro_D_B_id[i] = gRandom->Gaus(0.0,1.0);
    ro_D_E_id[i] = gRandom->Gaus(0.0,1.0);
    ro_M_B_id[i] = gRandom->Gaus(0.0,1.0);
    ro_M_E_id[i] = gRandom->Gaus(0.0,1.0);
    
    ro_D_B_hlt[i] = gRandom->Gaus(0.0,1.0);
    ro_D_E_hlt[i] = gRandom->Gaus(0.0,1.0);
    ro_M_B_hlt[i] = gRandom->Gaus(0.0,1.0);
    ro_M_E_hlt[i] = gRandom->Gaus(0.0,1.0);
  }
  */
  int debug_pseudo_exps=0;
  for (int i=0; i<nexp; i++) {
    EffArray_t *arr= & ro_Data[i];
    for (int kind=0; kind<NEffTypes; ++kind) {
      for (int iEt=0; iEt<DYTools::nEtBinsMax; ++iEt) {
	for (int iEta=0; iEta<DYTools::nEtaBinsMax; ++iEta) {
	  // In the special case of the RECO efficiency for low Et
	  // electrons, some eta bins are MERGED in tag and probe.
	  // However the binning is kept standard, so the values
	  // for the efficiencies in the merged bins are the same,
	  // and the errors are 100% correlated. Take this into
	  // account and make the smearing 100% correlated as well.
	  if( kind == 0 && etaBinning == DYTools::ETABINS5 
	      && (getEtBinLimits(etBinning))[iEt+1] <= 20.0 
	      && (iEta == 1 || iEta == 4)    ) {
	    // We are dealing with merged bins of the RECO efficiency with the right binning
	    // For iEta == 1 or 4, fall back to the values for iEta == 0 or 3.
	    (*arr)[kind][iEt][iEta]= (*arr)[kind][iEt][iEta-1];
	  }else{
	    // The default case, all other efficiencies and bins
	    (*arr)[kind][iEt][iEta]=
	      (debug_pseudo_exps) ? ((kind+1)*100 + (iEt+1)*10 + iEta+1) :
	      gRandom->Gaus(0.0,1.0);
	  }
	}
      }
    }
    arr= & ro_MC[i];
    for (int kind=0; kind<NEffTypes; ++kind) {
      for (int iEt=0; iEt<DYTools::nEtBinsMax; ++iEt) {
	for (int iEta=0; iEta<DYTools::nEtaBinsMax; ++iEta) {
	  // In the case of RECO efficiency with MERGED bins for tag and
	  // probe and ETABINS5, the randomization is different. See 
	  // comments above. 
	  if( kind == 0 && etaBinning == DYTools::ETABINS5 
	      && (getEtBinLimits(etBinning))[iEt+1] <= 20.0 
	      && (iEta == 1 || iEta == 4)    ) {
	    // We are dealing with merged bins of the RECO efficiency with the right binning
	    // For iEta == 1 or 4, fall back to the values for iEta == 0 or 3.
	    (*arr)[kind][iEt][iEta]= (*arr)[kind][iEt][iEta-1];
	  }else{
	    // The default case, all other efficiencies and bins
	    (*arr)[kind][iEt][iEta]=
	      (debug_pseudo_exps) ? -((kind+1)*100 + (iEt+1)*10 + iEta+1) :
	      gRandom->Gaus(0.0,1.0);
	  }
	}
      }
    }
  }
  
  // Create container for data for error estimates based on pseudo-experiments
  //TH1D *systScale[DYTools::nMassBins][nexp];

#ifdef debug_systScaleArrs
  SystTH1DArray_t systScale;
  SystTH1DArray_t systScaleReco;
  SystTH1DArray_t systScaleId;
  SystTH1DArray_t systScaleHlt;
  SystTH1DArray_t systScaleHltLeg1;
  SystTH1DArray_t systScaleHltLeg2;
#endif

  //TH1D *systScaleFI[nUnfoldingBins][nexp];
  SystTH1DArrayFI_t systScaleFI;
  SystTH1DArrayFI_t systScaleRecoFI;
  SystTH1DArrayFI_t systScaleIdFI;
  SystTH1DArrayFI_t systScaleHltFI;
  SystTH1DArrayFI_t systScaleHltLeg1FI;
  SystTH1DArrayFI_t systScaleHltLeg2FI;

  for(int i=0; i<nUnfoldingBins; i++) {
    for(int j=0; j<nexp; j++){
      TString base;
#ifdef debug_systScaleArrs
      if (i<DYTools::nMassBins) {
	base = "hScaleM_massBin";
	base += i;
	base += "_exp";
	base += j;
	systScale[i][j] = new TH1D(base,base,150,0.0,1.5);
	systScaleReco[i][j] = new TH1D(base+TString("_reco"),
				       base+TString("_reco"),150,0.0,1.5);
	systScaleId [i][j] = new TH1D(base+TString("_id" ),
				      base+TString("_id" ),150,0.0,1.5);
	systScaleHlt[i][j] = new TH1D(base+TString("_hlt"),
				      base+TString("_hlt"),150,0.0,1.5);
	systScaleHltLeg1[i][j] = (nonUniversalHLT) ?
	                 new TH1D(base+TString("_hltLeg1"),
				  base+TString("_hltLeg1"),150,0.0,1.5) : NULL;
	systScaleHltLeg2[i][j] = (nonUniversalHLT) ?
	                 new TH1D(base+TString("_hltLeg2"),
				  base+TString("_hltLeg2"),150,0.0,1.5) : NULL;
      }
#endif
      base = "hScaleM_flatIdx";
      base += i;
      base += "_exp";
      base += j;
      systScaleFI[i][j] = new TH1D(base,base,150,0.0,1.5);
      systScaleRecoFI[i][j] = new TH1D(base+TString("_reco"),
				       base+TString("_reco"),150,0.0,1.5);
      systScaleIdFI [i][j] = new TH1D(base+TString("_id" ),
				      base+TString("_id" ),150,0.0,1.5);
      systScaleHltFI[i][j] = new TH1D(base+TString("_hlt"),
				      base+TString("_hlt"),150,0.0,1.5);
      systScaleHltLeg1FI[i][j] = (nonUniversalHLT) ?
			new TH1D(base+TString("_hltLeg1"),
				 base+TString("_hltLeg1"),150,0.0,1.5) : NULL;
      systScaleHltLeg2FI[i][j] = (nonUniversalHLT) ?
			new TH1D(base+TString("_hltLeg2"),
				 base+TString("_hltLeg2"),150,0.0,1.5) : NULL;
    }
  }

  // for correlation studies
  TH1D *hEvtW=new TH1D("hEvtW","hEvtW",nUnfoldingBins,0.,double(nUnfoldingBins));
  TH1D *hEsfEvtW=new TH1D("hEsfEvtW","hEsfEvtW",nUnfoldingBins,0.,double(nUnfoldingBins));
  double sumEvtW_Zpeak=0., sumEsfEvtW_Zpeak=0.;
  std::vector<double> systSumEsfEvtW_ZpeakV(nexp);
  std::vector<TH1D*> hSystEsfEvtWV;
  hSystEsfEvtWV.reserve(nexp);
  for (int i=0; i<nexp; i++) {
    TString base = Form("hSystEsfEvtW_exp%d",i);
    TH1D* h=new TH1D(base,base,nUnfoldingBins,0.,double(nUnfoldingBins));
    h->SetDirectory(0);
    hSystEsfEvtWV.push_back(h);
  }

  PUReweight_t PUReweight(PUReweight_t::_Hildreth);
  // For Hildreth method of PU reweighting, the lines below are not needed
//   if (puReweight) {
//     assert(PUReweight.setDefaultFile(mcMgr.dirTag(),DYTools::analysisTag_USER,0));
//     assert(PUReweight.setReference("hNGoodPV_data"));
//     assert(PUReweight.setActiveSample("hNGoodPV_zee"));
//   }

  //HERE("break the macro"); return;


  if (selectEventsFName.Index("_DebugRun")!=-1) selectEventsFName.ReplaceAll("_DebugRun","");
  HERE("opening selectEventsFile=<%s>",selectEventsFName.Data());

  TFile *skimFile=new TFile(selectEventsFName);
  if (!skimFile || !skimFile->IsOpen()) {
    if (allowToIgnoreAnalysisTag) {
      std::cout << ".... changing analysis tag in selectEventsFName\n";
      if (DYTools::analysisTag=="2D") selectEventsFName.ReplaceAll("2D","1D");
      else selectEventsFName.ReplaceAll("1D","2D");
      if (skimFile) delete skimFile;
      skimFile=new TFile(selectEventsFName);
    }
    if (!skimFile || !skimFile->IsOpen()) {
      std::cout << "failed to open file <" << selectEventsFName << ">\n";
      assert(0);
    }
  }
  TTree *skimTree = (TTree*)skimFile->Get("Events");
  assert(skimTree);
  esfSelectEvent_t selData;
  selData.setBranchAddress(skimTree);

  //std::cout << "DYTools::nMassBins=" << DYTools::nMassBins << ", nUnfoldingBins=" << nUnfoldingBins << std::endl;
  ULong_t maxEvents=skimTree->GetEntries();
  std::cout << "there are " << skimTree->GetEntries() 
	    << " entries in the <" << selectEventsFName << "> file\n";
  for (UInt_t ientry=0; ientry < maxEvents; ++ientry) {
    if (DYTools::isDebugMode(runMode) && (ientry>10000)) break;
    printProgress(100000," ientry=",ientry,maxEvents);

    skimTree->GetEntry(ientry);

    double scaleFactor = findEventScaleFactor(-1, selData); // HLT formula is inside
    double scaleFactorReco = sqrt(findEventScaleFactor(DYTools::RECO,selData));
    double scaleFactorId  = sqrt(findEventScaleFactor(DYTools::ID,selData));
    double scaleFactorHlt = sqrt(findEventScaleFactor(DYTools::HLT,selData)); // HLT formula is inside
    double scaleFactorHltLeg1 = (nonUniversalHLT) ? 
      sqrt(findEventScaleFactor(DYTools::HLT_leg1,selData)) : 0.;
    double scaleFactorHltLeg2 = (nonUniversalHLT) ? 
      sqrt(findEventScaleFactor(DYTools::HLT_leg2,selData)) : 0;
    double weight=selData.weight;
    //if (puReweight) weight *= PUReweight.getWeightHildreth(selData.nGoodPV);
    if ( 0 || ( ientry%20000 == 0 )) std::cout << "ientry=" << ientry << ", weight=" << weight << ", scaleFactor=" << scaleFactor << "\n";

    hScale->Fill(scaleFactor, weight);
    hScaleReco->Fill( scaleFactorReco, weight);
    hScaleId ->Fill( scaleFactorId, weight);
    hScaleHlt->Fill( scaleFactorHlt, weight);
    if (nonUniversalHLT) {
      hScaleHltLeg1->Fill( scaleFactorHltLeg1, weight);
      hScaleHltLeg2->Fill( scaleFactorHltLeg2, weight);
    }

    // Use generator-level post-FSR mass, y
    int ibin= DYTools::findMassBin(selData.genMass);
    int idx = DYTools::findIndexFlat(selData.genMass,selData.genY);
    if ((ibin>=0) && (ibin<DYTools::nMassBins)) {
      hLeadingEtV [ibin]->Fill( selData.et_1, weight);
      hTrailingEtV[ibin]->Fill( selData.et_2, weight);
      hElectronEtV[ibin]->Fill( selData.et_1, weight);
      hElectronEtV[ibin]->Fill( selData.et_2, weight);
      if( selData.insideMassWindow(60,120) ) {
	hZpeakEt->Fill(selData.et_1, weight);
	hZpeakEt->Fill(selData.et_2, weight);
	if ((idx>=0) && (idx<nUnfoldingBins)) {
	  sumEvtW_Zpeak+=weight;
	  sumEsfEvtW_Zpeak+=weight*scaleFactor;
	}
      }

      hScaleRecoV[ibin]->Fill( scaleFactorReco, weight);
      hScaleIdV [ibin]->Fill( scaleFactorId, weight);
      hScaleHltV[ibin]->Fill( scaleFactorHlt, weight);
      if (nonUniversalHLT) {
	hScaleHltLeg1V[ibin]->Fill( scaleFactorHltLeg1, weight);
	hScaleHltLeg2V[ibin]->Fill( scaleFactorHltLeg2, weight);
      }
      hScaleV   [ibin]->Fill( scaleFactor, weight);
      //if (ibin==39) std::cout << " sf=" << scaleFactor << " * " << weight << "\n";

      if ((idx>=0) && (idx<nUnfoldingBins)) {
	hScaleRecoFIV[idx]->Fill( scaleFactorReco, weight);
	hScaleIdFIV [idx]->Fill( scaleFactorId, weight);
	hScaleHltFIV[idx]->Fill( scaleFactorHlt, weight);
	if (nonUniversalHLT) {
	  hScaleHltLeg1FIV[idx]->Fill( scaleFactorHltLeg1, weight);
	  hScaleHltLeg2FIV[idx]->Fill( scaleFactorHltLeg2, weight);
	}
	hScaleFIV   [idx]->Fill( scaleFactor, weight);
	hEvtW->Fill(idx,weight);
	hEsfEvtW->Fill(idx,scaleFactor*weight);
      }
	
      // Acumulate pseudo-experiments for error estimate
      for(int iexp = 0; iexp<nexp; iexp++){
	scaleFactor = findEventScaleFactorSmeared(-1, selData, iexp); // HLT formula is inside
	//if (ibin==39) std::cout << " rndSf= " << scaleFactor << " * " << weight << "\n";
	//std::cout << "findEventScaleFactor(selData)=" << findEventScaleFactor(selData) << ", smeared scaleFactor=" << scaleFactor << "\n";
	scaleFactorReco = sqrt(findEventScaleFactorSmeared(DYTools::RECO,selData,iexp));
	scaleFactorId  = sqrt(findEventScaleFactorSmeared(DYTools::ID,selData,iexp));
	scaleFactorHlt = sqrt(findEventScaleFactorSmeared(DYTools::HLT,selData,iexp)); // HLT formula is inside
	if (nonUniversalHLT) {
	  scaleFactorHltLeg1 = sqrt(findEventScaleFactorSmeared(DYTools::HLT_leg1,selData,iexp));
	  scaleFactorHltLeg2 = sqrt(findEventScaleFactorSmeared(DYTools::HLT_leg2,selData,iexp));
	}
#ifdef debug_systScaleArrs
	systScale    [ibin][iexp]->Fill(scaleFactor, weight);
	systScaleReco[ibin][iexp]->Fill(scaleFactorReco, weight);
	systScaleId [ibin][iexp]->Fill(scaleFactorId, weight);
	systScaleHlt[ibin][iexp]->Fill(scaleFactorHlt, weight);
	if (nonUniversalHLT) {
	  systScaleHltLeg1[ibin][iexp]->Fill(scaleFactorHltLeg1, weight);
	  systScaleHltLeg2[ibin][iexp]->Fill(scaleFactorHltLeg2, weight);
	}
#endif
	if ((idx>=0) && (idx<nUnfoldingBins)) {
	  systScaleFI    [idx][iexp]->Fill(scaleFactor, weight);
	  systScaleRecoFI[idx][iexp]->Fill(scaleFactorReco, weight);
	  systScaleIdFI [idx][iexp]->Fill(scaleFactorId, weight);
	  systScaleHltFI[idx][iexp]->Fill(scaleFactorHlt, weight);
	  if (nonUniversalHLT) {
	    systScaleHltLeg1FI[idx][iexp]->Fill(scaleFactorHltLeg1, weight);
	    systScaleHltLeg2FI[idx][iexp]->Fill(scaleFactorHltLeg2, weight);
	  }
	  hSystEsfEvtWV[iexp]->Fill(idx,scaleFactor*weight);
	  if( selData.insideMassWindow(60,120) ) {
	    systSumEsfEvtW_ZpeakV[iexp]+=weight*scaleFactor;
	  }
	}
      }
    } // if (ibin is ok)
    
    // 	if(scaleFactor>1.3)
    // 	  printf("  leading:   %f    %f      trailing:   %f   %f     mass: %f\n",
// 		 leading->scEt, leading->scEta, trailing->scEt, trailing->scEta, dielectron->mass);
    
  } // end loop over selected events
   
 
  delete skimTree;
  delete skimFile;

  std::cout << "loop over selected events done" << std::endl;
  //std::cout << "forced stop\n";
  //return retCodeStop;
  
  // Calculate errors on the scale factors
  // The "Mean" are the mean among all pseudo-experiments, very close to the primary scale factor values
#ifdef debug_systScaleArrs
  TVectorD scaleMeanV(DYTools::nMassBins);
  TVectorD scaleMeanErrV(DYTools::nMassBins);
  TVectorD scaleMeanRecoV(DYTools::nMassBins);
  TVectorD scaleMeanRecoErrV(DYTools::nMassBins);
  TVectorD scaleMeanIdV(DYTools::nMassBins);
  TVectorD scaleMeanIdErrV(DYTools::nMassBins);
  TVectorD scaleMeanHltV(DYTools::nMassBins);
  TVectorD scaleMeanHltErrV(DYTools::nMassBins);
  TVectorD scaleMeanHltLeg1V(DYTools::nMassBins);
  TVectorD scaleMeanHltLeg1ErrV(DYTools::nMassBins);
  TVectorD scaleMeanHltLeg2V(DYTools::nMassBins);
  TVectorD scaleMeanHltLeg2ErrV(DYTools::nMassBins);
  // Put into these vectors the content of the mean of the primary scale factor distributions
  TVectorD scaleV(DYTools::nMassBins);
  TVectorD scaleRecoV(DYTools::nMassBins);
  TVectorD scaleIdV(DYTools::nMassBins);
  TVectorD scaleHltV(DYTools::nMassBins);
  TVectorD scaleHltLeg1V(DYTools::nMassBins);
  TVectorD scaleHltLeg2V(DYTools::nMassBins);
#endif

  TVectorD scaleMeanFIV(nUnfoldingBins);
  TVectorD scaleMeanErrFIV(nUnfoldingBins);
  TVectorD scaleMeanRecoFIV(nUnfoldingBins);
  TVectorD scaleMeanRecoErrFIV(nUnfoldingBins);
  TVectorD scaleMeanIdFIV(nUnfoldingBins);
  TVectorD scaleMeanIdErrFIV(nUnfoldingBins);
  TVectorD scaleMeanHltFIV(nUnfoldingBins);
  TVectorD scaleMeanHltErrFIV(nUnfoldingBins);
  TVectorD scaleMeanHltLeg1FIV(nUnfoldingBins);
  TVectorD scaleMeanHltLeg1ErrFIV(nUnfoldingBins);
  TVectorD scaleMeanHltLeg2FIV(nUnfoldingBins);
  TVectorD scaleMeanHltLeg2ErrFIV(nUnfoldingBins);
  // Put into these vectors the content of the mean of the primary scale factor distributions
  TVectorD scaleFIV(nUnfoldingBins);
  TVectorD scaleRecoFIV(nUnfoldingBins);
  TVectorD scaleIdFIV(nUnfoldingBins);
  TVectorD scaleHltFIV(nUnfoldingBins);
  TVectorD scaleHltLeg1FIV(nUnfoldingBins);
  TVectorD scaleHltLeg2FIV(nUnfoldingBins);

#ifdef debug_systScaleArrs
  deriveScaleMeanAndErr(DYTools::nMassBins,nexp, 
			systScale, scaleMeanV,scaleMeanErrV);
  deriveScaleMeanAndErr(DYTools::nMassBins,nexp, 
			systScaleReco, scaleMeanRecoV,scaleMeanRecoErrV);
  deriveScaleMeanAndErr(DYTools::nMassBins,nexp, 
			systScaleId  , scaleMeanIdV  ,scaleMeanIdErrV  );
  deriveScaleMeanAndErr(DYTools::nMassBins,nexp, 
			systScaleHlt , scaleMeanHltV ,scaleMeanHltErrV );
  if (nonUniversalHLT) {
    deriveScaleMeanAndErr(DYTools::nMassBins,nexp, 
			  systScaleHltLeg1 , scaleMeanHltLeg1V ,scaleMeanHltLeg1ErrV );
    deriveScaleMeanAndErr(DYTools::nMassBins,nexp, 
			  systScaleHltLeg2 , scaleMeanHltLeg2V ,scaleMeanHltLeg2ErrV );
  }
#endif

  deriveScaleMeanAndErr(nUnfoldingBins,nexp, 
			systScaleFI, scaleMeanFIV,scaleMeanErrFIV);
  deriveScaleMeanAndErr(nUnfoldingBins,nexp, 
			systScaleRecoFI, scaleMeanRecoFIV,scaleMeanRecoErrFIV);
  if (1) {
  deriveScaleMeanAndErr(nUnfoldingBins,nexp, 
			systScaleIdFI  , scaleMeanIdFIV  ,scaleMeanIdErrFIV  );
  deriveScaleMeanAndErr(nUnfoldingBins,nexp, 
			systScaleHltFI , scaleMeanHltFIV ,scaleMeanHltErrFIV );
  if (nonUniversalHLT) {
    deriveScaleMeanAndErr(nUnfoldingBins,nexp, 
			  systScaleHltLeg1FI , scaleMeanHltLeg1FIV ,scaleMeanHltLeg1ErrFIV );
    deriveScaleMeanAndErr(nUnfoldingBins,nexp, 
			  systScaleHltLeg2FI , scaleMeanHltLeg2FIV ,scaleMeanHltLeg2ErrFIV );
  }
  }


#ifdef debug_systScaleArrs
  for(int ibin = 0; ibin < DYTools::nMassBins; ibin++){
    scaleRecoV[ibin] = hScaleRecoV[ibin]->GetMean();
    scaleIdV [ibin] = hScaleIdV [ibin] ->GetMean();
    scaleHltV[ibin] = hScaleHltV[ibin]->GetMean();
    scaleV   [ibin] = hScaleV   [ibin]->GetMean();
  }
#endif

  for(int idx = 0; idx < nUnfoldingBins; idx++){
    scaleRecoFIV[idx] = hScaleRecoFIV[idx]->GetMean();
    scaleIdFIV [idx] = hScaleIdFIV [idx] ->GetMean();
    scaleHltFIV[idx] = hScaleHltFIV[idx]->GetMean();
    scaleFIV   [idx] = hScaleFIV   [idx]->GetMean();
  }
  if (nonUniversalHLT) {
#ifdef debug_systScaleArrs
    for(int ibin = 0; ibin < DYTools::nMassBins; ibin++){
      scaleHltLeg1V[ibin] = hScaleHltLeg1V[ibin]->GetMean();
      scaleHltLeg2V[ibin] = hScaleHltLeg2V[ibin]->GetMean();
    }
#endif
    for(int idx = 0; idx < nUnfoldingBins; idx++){
      scaleHltLeg1FIV[idx] = hScaleHltLeg1FIV[idx]->GetMean();
      scaleHltLeg2FIV[idx] = hScaleHltLeg2FIV[idx]->GetMean();
    }
  }

  /* superceded
  for(int ibin = 0; ibin < DYTools::nMassBins; ibin++){
    scaleMeanV[ibin] = 0;
    scaleMeanErrV[ibin] = 0;
    scaleMeanRecoV[ibin] = 0;
    scaleMeanRecoErrV[ibin] = 0;
    scaleMeanIdV[ibin] = 0;
    scaleMeanIdErrV[ibin] = 0;
    scaleMeanHltV[ibin] = 0;
    scaleMeanHltErrV[ibin] = 0;
    for(int iexp = 0; iexp < nexp; iexp++){
      scaleMeanV[ibin] += systScale[ibin][iexp]->GetMean();
      scaleMeanErrV[ibin] += SQR(systScale[ibin][iexp]->GetMean());

      scaleMeanRecoV[ibin] += systScaleReco[ibin][iexp]->GetMean();
      scaleMeanRecoErrV[ibin] += SQR(systScaleReco[ibin][iexp]->GetMean());

      scaleMeanIdV[ibin] += systScaleId[ibin][iexp]->GetMean();
      scaleMeanIdErrV[ibin] += SQR(systScaleId[ibin][iexp]->GetMean());

      scaleMeanHltV[ibin] += systScaleHlt[ibin][iexp]->GetMean();
      scaleMeanHltErrV[ibin] += SQR(systScaleHlt[ibin][iexp]->GetMean());
    }
    scaleRecoV[ibin] = hScaleRecoV[ibin]->GetMean();
    scaleIdV [ibin] = hScaleIdV [ibin] ->GetMean();
    scaleHltV[ibin] = hScaleHltV[ibin]->GetMean();
    scaleV   [ibin] = hScaleV   [ibin]->GetMean();

    scaleMeanV[ibin] = scaleMeanV[ibin]/double(nexp);
    scaleMeanErrV[ibin] = sqrt( scaleMeanErrV[ibin] / double(nexp) 
				- scaleMeanV[ibin]*scaleMeanV[ibin] ); 
				
    scaleMeanRecoV[ibin] = scaleMeanRecoV[ibin]/double(nexp);
    scaleMeanRecoErrV[ibin] = sqrt( scaleMeanRecoErrV[ibin] / double(nexp)
				- scaleMeanRecoV[ibin]*scaleMeanRecoV[ibin] ); 
				
    scaleMeanIdV[ibin] = scaleMeanIdV[ibin]/double(nexp);
    scaleMeanIdErrV[ibin] = sqrt( scaleMeanIdErrV[ibin] / double(nexp) 
				- scaleMeanIdV[ibin]*scaleMeanIdV[ibin] ); 
				
    scaleMeanHltV[ibin] = scaleMeanHltV[ibin]/double(nexp);
    scaleMeanHltErrV[ibin] = sqrt( scaleMeanHltErrV[ibin] / double(nexp) 
				- scaleMeanHltV[ibin]*scaleMeanHltV[ibin] ); 
				
  }
  for(int idx = 0; idx < nUnfoldingBins; idx++){
    scaleMeanFIV[idx] = 0;
    scaleMeanErrFIV[idx] = 0;
    scaleMeanRecoFIV[idx] = 0;
    scaleMeanRecoErrFIV[idx] = 0;
    scaleMeanIdFIV[idx] = 0;
    scaleMeanIdErrFIV[idx] = 0;
    scaleMeanHltFIV[idx] = 0;
    scaleMeanHltErrFIV[idx] = 0;
    for(int iexp = 0; iexp < nexp; iexp++){
      scaleMeanFIV[idx] += systScaleFI[idx][iexp]->GetMean();
      scaleMeanErrFIV[idx] += SQR(systScaleFI[idx][iexp]->GetMean());

      scaleMeanRecoFIV[idx] += systScaleRecoFI[idx][iexp]->GetMean();
      scaleMeanRecoErrFIV[idx] += SQR(systScaleRecoFI[idx][iexp]->GetMean());

      scaleMeanIdFIV[idx] += systScaleIdFI[idx][iexp]->GetMean();
      scaleMeanIdErrFIV[idx] += SQR(systScaleIdFI[idx][iexp]->GetMean());

      scaleMeanHltFIV[idx] += systScaleHltFI[idx][iexp]->GetMean();
      scaleMeanHltErrFIV[idx] += SQR(systScaleHltFI[idx][iexp]->GetMean());
    }
    scaleRecoFIV[idx] = hScaleRecoFIV[idx]->GetMean();
    scaleIdFIV [idx] = hScaleIdFIV [idx] ->GetMean();
    scaleHltFIV[idx] = hScaleHltFIV[idx]->GetMean();
    scaleFIV   [idx] = hScaleFIV   [idx]->GetMean();

    scaleMeanFIV[idx] = scaleMeanFIV[idx]/double(nexp);
    scaleMeanErrFIV[idx] = sqrt( scaleMeanErrFIV[idx] / double(nexp) 
				- scaleMeanFIV[idx]*scaleMeanFIV[idx] ); 
				
    scaleMeanRecoFIV[idx] = scaleMeanRecoFIV[idx]/double(nexp);
    scaleMeanRecoErrFIV[idx] = sqrt( scaleMeanRecoErrFIV[idx] / double(nexp)
				- scaleMeanRecoFIV[idx]*scaleMeanRecoFIV[idx] ); 
				
    scaleMeanIdFIV[idx] = scaleMeanIdFIV[idx]/double(nexp);
    scaleMeanIdErrFIV[idx] = sqrt( scaleMeanIdErrFIV[idx] / double(nexp) 
				- scaleMeanIdFIV[idx]*scaleMeanIdFIV[idx] ); 
				
    scaleMeanHltFIV[idx] = scaleMeanHltFIV[idx]/double(nexp);
    scaleMeanHltErrFIV[idx] = sqrt( scaleMeanHltErrFIV[idx] / double(nexp) 
				- scaleMeanHltFIV[idx]*scaleMeanHltFIV[idx] ); 
				
  }
  */

  // Store constants in the file

  TMatrixD scaleMatrix(DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD scaleMatrixErr(DYTools::nMassBins,DYTools::nYBinsMax);
  deflattenMatrix(scaleFIV, scaleMatrix);
  deflattenMatrix(scaleMeanErrFIV, scaleMatrixErr);

  TVectorD vecEtBins(etBinCount+1);
  for (int i=0; i<etBinCount+1; ++i) vecEtBins[i]=etBinLimits[i];
  TVectorD vecEtaBins(etaBinCount+1);
  for (int i=0; i<etaBinCount+1; ++i) vecEtaBins[i]=etaBinLimits[i];

  TFile fa(sfConstFileName, "recreate");
#ifdef debug_systScaleArrs
  scaleV.Write("scaleFactorMassIdxArray");
  scaleMeanErrV.Write("scaleFactorErrMassIdxArray");
#endif
  scaleFIV.Write("scaleFactorFlatIdxArray");
  scaleMeanErrFIV.Write("scaleFactorErrFlatIdxArray");
  scaleMatrix.Write("scaleFactor");
  scaleMatrixErr.Write("scaleFactorErr");
  vecEtBins.Write("etBinLimits");
  vecEtaBins.Write("etaBinLimits");

  scaleRecoFIV.Write("scaleRecoFactor_FIArray");
  scaleMeanRecoErrFIV.Write("scaleRecoFactorErr_FIArray");
  scaleIdFIV.Write("scaleIdFactor_FIArray");
  scaleMeanIdErrFIV.Write("scaleIdFactorErr_FIArray");
  scaleHltFIV.Write("scaleHltFactor_FIArray");
  scaleMeanHltErrFIV.Write("scaleHltFactorErr_FIArray");
  if (nonUniversalHLT) {
    scaleHltLeg1FIV.Write("scaleHltLeg1Factor_FIArray");
    scaleMeanHltLeg1ErrFIV.Write("scaleHltLeg1FactorErr_FIArray");
    scaleHltLeg2FIV.Write("scaleHltLeg2Factor_FIArray");
    scaleMeanHltLeg2ErrFIV.Write("scaleHltLeg2FactorErr_FIArray");
  }
  writeBinningArrays(fa);
  fa.Close();

  //std::cout << "forced stop 2\n";   return retCodeStop;

  //
  // Correlation study
  //
  
  if (correlationStudy) {
    TString corrFileNameDebug=sfConstFileName;
    corrFileNameDebug.ReplaceAll("scale_factors","esf_correlation_debug");
    TFile fCorrDebug(corrFileNameDebug,"recreate");
    if (!fCorrDebug.IsOpen()) {
      std::cout << "failed to create a file <" << corrFileNameDebug << ">" << std::endl;
      assert(0);
    }
    
    const double c_esfLimit_low = 0.0;
    const double c_esfLimit_high = 1.5;
    const int c_theEsfBins = int((c_esfLimit_high-c_esfLimit_low)/5e-3 + 1e-3);

    // omit the underflow mass bin in 2D case
    const int iMmin=(DYTools::study2D==0) ? 0 : 1;
    const int idxMin= iMmin * DYTools::nYBins[0];
    const int nCorrelationBins=nUnfoldingBins-idxMin;

    TH2D *hCorrelation=new TH2D("hCorrelation","hCorrelation",nCorrelationBins,0.5,nCorrelationBins+0.5,nCorrelationBins,0.5,nCorrelationBins+0.5);
    hCorrelation->SetDirectory(0);
    TH2D *hCorrelationNorm=(TH2D*)hCorrelation->Clone("hCorrelation_Norm");
    hCorrelationNorm->SetDirectory(0);

    std::vector<TH1D*> hRatioV, hRatio_NormV;
    hRatioV.reserve(nexp);
    hRatio_NormV.reserve(nexp);
    for (int iexp=0; iexp<nexp; ++iexp) {
      TString nameHRatio= Form("hRatio_iexp%d",iexp);
      TH1D *hRatio=(TH1D*)hSystEsfEvtWV[iexp]->Clone(nameHRatio);
      hRatio->SetDirectory(0);
      hRatio->Divide(hSystEsfEvtWV[iexp],hEvtW);
      hRatioV.push_back(hRatio);
      
      TString nameHRatioNorm= Form("hRatio_Norm_iexp%d",iexp);
      TH1D *hRatio_Norm=(TH1D*)hSystEsfEvtWV[iexp]->Clone(nameHRatioNorm);
      hRatio_Norm->SetDirectory(0);
      hRatio_Norm->Divide(hSystEsfEvtWV[iexp],hEvtW,sumEvtW_Zpeak/systSumEsfEvtW_ZpeakV[iexp],1.);
      hRatio_NormV.push_back(hRatio_Norm);
    }
    
    for (int iM=iMmin, i=idxMin; iM<DYTools::nMassBins; ++iM) {
      for (int iY=0; (iY<DYTools::nYBins[iM]) && (i<nUnfoldingBins); ++iY, ++i) {
	for (int jM=iMmin, j=idxMin; jM<DYTools::nMassBins; ++jM) {
	  for (int jY=0; (jY<DYTools::nYBins[jM]) && (j<nUnfoldingBins); ++jY, ++j) {
	    TString nameHESF= Form("hESF_%d_%d__iM%d_iY%d__jM%d_jY%d",i-idxMin,j-idxMin,iM-iMmin,iY,jM-iMmin,jY);
	    TString nameHESF_Norm= Form("hESFNorm_%d_%d__iM%d_iY%d__jM%d_jY%d",i-idxMin,j-idxMin,iM-iMmin,iY,jM-iMmin,jY);
	    if (DYTools::isDebugMode(runMode)) HERE(nameHESF.Data());
	    TH2D *hESF=new TH2D(nameHESF,nameHESF,c_theEsfBins,c_esfLimit_low,c_esfLimit_high, c_theEsfBins,c_esfLimit_low,c_esfLimit_high);
	    TH2D *hESF_Norm=(TH2D*)hESF->Clone(nameHESF_Norm);
	    hESF->SetDirectory(0); hESF_Norm->SetDirectory(0);
	    
	    for (int iexp=0; iexp<nexp; ++iexp) {
	      TH1D *hRatio=hRatioV[iexp];
	      hESF->Fill(hRatio->GetBinContent(i+1),hRatio->GetBinContent(j+1));
	      TH1D *hRatio_Norm=hRatio_NormV[iexp];
	      hESF_Norm->Fill(hRatio_Norm->GetBinContent(i+1),hRatio_Norm->GetBinContent(j+1));
	    }
	    
	    if ((jM==0) && (jY==0)) {
	      TString tmpDirName=Form("%dDcorrel_histos_McorrBin%d_Ybin%d",(DYTools::study2D) ? 2:1,iM+1-iMmin,iY+1);
	      fCorrDebug.cd();
	      fCorrDebug.mkdir(tmpDirName);
	      fCorrDebug.cd(tmpDirName);
	    }
	    
	    hESF->Write();
	    hESF_Norm->Write();
	    
	    hCorrelation->Fill(i+1-idxMin,j+1-idxMin, hESF->GetCorrelationFactor(1,2) );
	    hCorrelationNorm->Fill(i+1-idxMin,j+1-idxMin, hESF_Norm->GetCorrelationFactor(1,2) );
	    delete hESF;
	    delete hESF_Norm;
	  }
	}
      }
    }
    fCorrDebug.cd();
    hCorrelation->Write();
    hCorrelationNorm->Write();
    fCorrDebug.Close();
    std::cout << "file <" << corrFileNameDebug << "> created" << std::endl;

    TString corrFileName=corrFileNameDebug;
    corrFileName.ReplaceAll("_debug","");
    TFile fCorr(corrFileName,"recreate");
    if (!fCorr.IsOpen()) {
      std::cout << "failed to create a file <" << corrFileName << ">" << std::endl;
      assert(0);
    }

    fCorr.cd();
    hCorrelation->Write();
    hCorrelationNorm->Write();

    TCanvas *canvCorr = MakeCanvas("canvEffCorr","canvEffCorr",600,600);
    canvCorr->SetRightMargin(0.12);
    CPlot plotEffCorr("efficiencyCorrelations","",
		       "idx_{1}",
		       "idx_{2}");
    gStyle->SetPalette(1);
    plotEffCorr.AddHist2D(hCorrelation,"COLZ");
    plotEffCorr.Draw(canvCorr);
    SaveCanvas(canvCorr,"figEffCorrelations");
    canvCorr->Write();

    TCanvas *canvCorrNorm = MakeCanvas("canvEffCorrNorm","canvEffCorrNorm",600,600);
    canvCorrNorm->SetRightMargin(0.12);
    CPlot plotEffCorrNorm("efficiencyCorrelations_Norm","",
                       "idx_{1}",
                       "idx_{2}");
    gStyle->SetPalette(1);
    plotEffCorrNorm.AddHist2D(hCorrelationNorm,"COLZ");
    plotEffCorrNorm.Draw(canvCorrNorm);
    SaveCanvas(canvCorrNorm,"figEffCorrelations_Norm");
    canvCorrNorm->Write();

    fCorr.Close();
    std::cout << "file <" << corrFileName << "> created" << std::endl;

    for (unsigned int i=0; i<hRatioV.size(); ++i) delete hRatioV[i];
    for (unsigned int i=0; i<hRatio_NormV.size(); ++i) delete hRatio_NormV[i];
  }

  //
  // Some plots
  //

  TString plotsFName=sfConstFileName;
  plotsFName.Replace(plotsFName.Index(".root"),sizeof(".root"),"-plots.root");
  TFile *faPlots=new TFile(plotsFName,"recreate");
  assert(faPlots);

  int c1Height=500+nonUniversalHLT*300;
  TCanvas *c1 = 
    new TCanvas("canvScaleFactors","canvScaleFactors",10,10,500,c1Height);
  c1->Divide(2,2+nonUniversalHLT);
  c1->cd(1);   hScale->Draw();
  c1->cd(2);   hScaleReco->Draw();
  c1->cd(3);   hScaleId->Draw();
  c1->cd(4);   hScaleHlt->Draw();
  if (nonUniversalHLT) {
    c1->cd(5); hScaleHltLeg1->Draw();
    c1->cd(6); hScaleHltLeg2->Draw();
  }
  c1->Write(); // save to file

#ifdef debug_systScaleArrs
  if (1) {
    std::cout << "\nScale factors as a function of mass bin\n";
    std::cout << "    mass          rho_reco          rho_id"
	      << "       rho_hlt        rho_total\n";
    std::string format1=
      std::string("   %3.0f - %3.0f     %5.3f +- %5.3f    %5.3f +- %5.3f")+
      std::string("     %5.3f +- %5.3f     %5.3f +- %5.3f\n");
    for(int i=0; i<DYTools::nMassBins; i++){
      printf(format1.c_str(),
	     DYTools::massBinLimits[i], DYTools::massBinLimits[i+1],
	     hScaleRecoV[i]->GetMean(), scaleMeanRecoErrV[i],
	     hScaleIdV[i]->GetMean(),  scaleMeanIdErrV[i],
	     hScaleHltV[i]->GetMean(), scaleMeanHltErrV[i],
	     hScaleV[i]->GetMean()   , scaleMeanErrV[i]);
    }
  }

  if (1 && nonUniversalHLT) {
    std::cout << "\nHLT scale factors as a function of mass bin\n";
    std::cout << "    mass          rho_hlt_leg1          rho_hlt_leg2"
	      << "       rho_hlt\n";
    std::string format1=
      std::string("   %3.0f - %3.0f     %5.3f +- %5.3f    %5.3f +- %5.3f")+
      std::string("     %5.3f +- %5.3f\n");
    for(int i=0; i<DYTools::nMassBins; i++){
      printf(format1.c_str(),
	     DYTools::massBinLimits[i], DYTools::massBinLimits[i+1],
	     hScaleHltLeg1V[i]->GetMean(), scaleMeanHltLeg1ErrV[i],
	     hScaleHltLeg2V[i]->GetMean(), scaleMeanHltLeg2ErrV[i],
	     hScaleHltV[i]->GetMean(), scaleMeanHltErrV[i]);
    }
  }
#endif

  if (1) {
    std::cout << "\nScale factors as a function of mass and rapidity bin\n";
    std::cout << "    mass     rapidity     rho_reco          rho_id"
	      << "       rho_hlt        rho_total\n";
    std::string format2=
      std::string("   %3.0f - %3.0f  %4.2f-%4.2f     %5.3f +- %5.3f    %5.3f +- %5.3f")+
      std::string("     %5.3f +- %5.3f     %5.3f +- %5.3f\n");
    for(int i=0; i<DYTools::nMassBins; i++){
      double *rapidityBinLimits=DYTools::getYBinArray(i);
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	int idx=DYTools::findIndexFlat(i,yi);
	std::cout << "idx=" << idx << "\n";
	printf(format2.c_str(),
	       DYTools::massBinLimits[i], DYTools::massBinLimits[i+1],
	       rapidityBinLimits[yi], rapidityBinLimits[yi+1],
	       hScaleRecoFIV[idx]->GetMean(), scaleMeanRecoErrFIV[idx],
	       hScaleIdFIV[idx]->GetMean(),  scaleMeanIdErrFIV[idx],
	       hScaleHltFIV[idx]->GetMean(), scaleMeanHltErrFIV[idx],
	       hScaleFIV[idx]->GetMean()   , scaleMeanErrFIV[idx]);
      }
    }
  }

  if (1 && nonUniversalHLT) {
    std::cout << "\nHLT scale factors as a function of mass and rapidity bin\n";
    std::cout << "    mass     rapidity     rho_hlt_leg1          rho_hlt_leg2"
	      << "       rho_hlt\n";
    std::string format2=
      std::string("   %3.0f - %3.0f  %4.2f-%4.2f     %5.3f +- %5.3f    %5.3f +- %5.3f")+
      std::string("     %5.3f +- %5.3f\n");
    for(int i=0; i<DYTools::nMassBins; i++){
      double *rapidityBinLimits=DYTools::getYBinArray(i);
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	int idx=DYTools::findIndexFlat(i,yi);
	std::cout << "idx=" << idx << "\n";
	printf(format2.c_str(),
	       DYTools::massBinLimits[i], DYTools::massBinLimits[i+1],
	       rapidityBinLimits[yi], rapidityBinLimits[yi+1],
	       hScaleHltLeg1FIV[idx]->GetMean(), scaleMeanHltLeg1ErrFIV[idx],
	       hScaleHltLeg2FIV[idx]->GetMean(), scaleMeanHltLeg2ErrFIV[idx],
	       hScaleHltFIV[idx]->GetMean(), scaleMeanHltErrFIV[idx]);
      }
    }
  }

  drawEfficiencies(faPlots);
  HERE("next: drawScaleFactors");
  drawScaleFactors(faPlots);

#ifdef debug_systScaleArrs
  if (1)
  drawEventScaleFactors(scaleRecoV, scaleMeanRecoErrV,
			scaleIdV , scaleMeanIdErrV ,
			scaleHltV, scaleMeanHltErrV,
			scaleV   , scaleMeanErrV   ,
			faPlots);

  if (1 && nonUniversalHLT) {
    drawEventScaleFactorsHLT(scaleHltLeg1V, scaleMeanHltLeg1ErrV,
			     scaleHltLeg2V, scaleMeanHltLeg2ErrV ,
			     faPlots);
  }
#endif


  const int rapidityIdxCount= (DYTools::study2D==1) ? 3 : 1;
  const int rapidityIdx_2D[] = { 0, 10, 20 };
  const int rapidityIdx_1D[] = {0};
  const int *rapidityIdx = (DYTools::study2D==1) ? rapidityIdx_2D : rapidityIdx_1D;
  if (1) {
    HERE("\n\tmake yDepSF plots");
    std::vector<CPlot*> cplotsV;
    if (1) {
      CPlot *yDepSF_Full= new CPlot("yDepSF_Full","", "m(e^{+}e^{-}) [GeV]", "scale factor");
      CPlot *yDepSF_Reco= new CPlot("yDepSF_Reco","", "m(e^{+}e^{-}) [GeV]", "sqrt( RECO scale factor )");
      CPlot *yDepSF_ID  = new CPlot("yDepSF_Id"  ,"", "m(e^{+}e^{-}) [GeV]", "sqrt( ID scale factor )");
      CPlot *yDepSF_HLT = new CPlot("yDepSF_Hlt" ,"", "m(e^{+}e^{-}) [GeV]", "sqrt( HLT scale factor )");
      CPlot *yDepSF_HLTleg1 = (nonUniversalHLT) ? new CPlot("yDepSF_Hlt_leg1" ,"", "m(e^{+}e^{-}) [GeV]", "sqrt( HLT_leg1 scale factor )") : NULL;
      CPlot *yDepSF_HLTleg2 = (nonUniversalHLT) ? new CPlot("yDepSF_Hlt_leg2" ,"", "m(e^{+}e^{-}) [GeV]", "sqrt( HLT_leg2 scale factor )") : NULL;
      cplotsV.reserve(NEffTypes+1);
      cplotsV.push_back(yDepSF_Full); cplotsV.push_back(yDepSF_Reco);
      cplotsV.push_back(yDepSF_ID);   cplotsV.push_back(yDepSF_HLT);
      if (nonUniversalHLT) {
	cplotsV.push_back(yDepSF_HLTleg1);
	cplotsV.push_back(yDepSF_HLTleg2);
      }
      for (unsigned int i=0; i<cplotsV.size(); ++i) {
	cplotsV[i]->SetLogx();
	cplotsV[i]->SetYRange(0.5, 1.5);
	cplotsV[i]->AddLine(0,1.0, 1500,1.0, kBlack, kDashed);
      }
    }
    HERE("calling drawEventScaleFactorsFI");
    for (int ri=0; ri<rapidityIdxCount; ++ri) {
      HERE("ri=%d",ri);
      drawEventScaleFactorsFI(scaleRecoFIV, scaleMeanRecoErrFIV,
			      scaleIdFIV , scaleMeanIdErrFIV ,
			      scaleHltFIV, scaleMeanHltErrFIV,
			      scaleFIV   , scaleMeanErrFIV   ,
			      rapidityIdx[ri],
			      faPlots, &cplotsV
			      );
    }

    if (1 && nonUniversalHLT) {
      for (int ri=0; ri<rapidityIdxCount; ++ri) {
	HERE("nonUniHLT ri=%d",ri);
	drawEventScaleFactorsHltFI(scaleHltLeg1FIV, scaleMeanHltLeg1ErrFIV,
				   scaleHltLeg2FIV , scaleMeanHltLeg2ErrFIV ,
				   rapidityIdx[ri],
				   faPlots, &cplotsV
				   );
      }
    }


    // This part does not work: combined plots are empty for some reason
    /*
    if (cplotsV.size()>=4) {
      HERE("cplotsV.size>=4 (it is %d)",cplotsV.size());
      TCanvas *c4 = MakeCanvas("canvYDepScaleFactors","canvYDepScaleFactors",800,800);
      c4->Divide(2,2);
      cplotsV[0]->Draw(c4,false,"png",1);
      cplotsV[1]->Draw(c4,false,"png",2);
      cplotsV[2]->Draw(c4,false,"png",3);
      cplotsV[3]->Draw(c4,savePlots,"png",4);
      c4->Update();
      TString plotName=TString("figYDepScaleFactors") + DYTools::analysisTag + puStr
	+ TString(".png");
      c4->SaveAs(plotName);
      c4->Write();
    }

    if (cplotsV.size()==6) {
      HERE("\ncplotsV.size==6");
      TCanvas *c3 = MakeCanvas("canvYDepHltScaleFactors","canvYDepHltScaleFactors",900,300);
      c3->Divide(3,1);
      cplotsV[4]->Draw(c3,false,"png",1);
      cplotsV[5]->Draw(c3,false,"png",2);
      cplotsV[3]->Draw(c3,savePlots,"png",3);
      c3->Update();
      TString plotName=TString("figYDepHltScaleFactors") + DYTools::analysisTag + puStr
	+ TString(".png");
      c3->SaveAs(plotName);
      c3->Write();
    }
    */
  }

  //
  // Make plots of Et spectra
  //
  HERE("plot Et spectra");

  // Normalize first
  for(int i=0; i<DYTools::nMassBins; i++){
    printf("Total events in mass bin %3d     %10.0f\n", 
	   i, hLeadingEtV[i]->GetSumOfWeights());
    hLeadingEtV  [i]->Sumw2();
    hTrailingEtV [i]->Sumw2();
    hElectronEtV [i]->Sumw2();
    hLeadingEtV [i]->Scale(1.0/hLeadingEtV[i]->GetSumOfWeights());
    hTrailingEtV[i]->Scale(1.0/hTrailingEtV[i]->GetSumOfWeights());
    hElectronEtV[i]->Scale(1.0/hElectronEtV[i]->GetSumOfWeights());
  }
  printf("Total events around Z peak    %10.0f\n", 
	 hZpeakEt->GetSumOfWeights()/2.0);
  hZpeakEt->Sumw2();
  hZpeakEt->Scale(1.0/hZpeakEt->GetSumOfWeights());


  const int plotMassBinIdxCount=5;
  const int plotMassBinIdx[plotMassBinIdxCount] = { 0, 3, 4, 5, DYTools::nMassBins-1 };
  std::vector<TCanvas*> canvV;
  canvV.reserve(plotMassBinIdxCount);

  for (int plot_i=0; plot_i<plotMassBinIdxCount; ++plot_i) {
    char buf[30];
    sprintf(buf,"canvET_%d",plot_i+1);
    TCanvas *c3 = MakeCanvas(buf,buf);
    canvV.push_back(c3);
    sprintf(buf,"etplot_bin%d",plotMassBinIdx[plot_i]);
    CPlot etplot1(buf, "","E_{T} [GeV]", "N_{ele}, normalized");
    const int cbin=plotMassBinIdx[plot_i];
    if (cbin<=DYTools::nMassBins) {
      sprintf(buf,"%1.0f-%1.0f mass range",DYTools::massBinLimits[cbin],DYTools::massBinLimits[cbin+1]);
    }
    else {
      std::cout << "ERROR: cbin=" << cbin << " is greater than DYTools::nMassBins=" << DYTools::nMassBins << "\n";
      break;
    }
    TString label = buf;
    etplot1.SetLogx();
    etplot1.AddHist1D(hElectronEtV[cbin] , label , "PE", kRed);
    etplot1.AddHist1D(hZpeakEt       , "60-120 mass range", "hist,f", kBlack);
    if (cbin==DYTools::nMassBins-1) etplot1.SetYRange(0.0,1.0);
    else etplot1.SetYRange(0.0,0.6);
    hElectronEtV[cbin]->GetXaxis()->SetMoreLogLabels();
    hElectronEtV[cbin]->GetXaxis()->SetNoExponent();
    hZpeakEt->SetFillStyle(3001);
    hElectronEtV[cbin]->SetMarkerColor(kRed);
    etplot1.Draw(c3,savePlots,"png");
    c3->Write();
  }

  saveVec(*faPlots,hLeadingEtV,"leadingEt");
  saveVec(*faPlots,hTrailingEtV,"trailingEt");
  saveVec(*faPlots,hElectronEtV,"electronEt");
  
  faPlots->Close();
  gBenchmark->Show("calcEventEff");
  return retCodeOk;
}

// ------------------------------------------------------------
/*
Bool_t matchedToGeneratorLevel(const TGenInfo *gen, 
			       const TDielectron *dielectron){

  Bool_t result = kTRUE;
  // In the generator branch of this ntuple, first particle is always
  // negative, and second always positive. In the Dielectron block
  // of the ntuple, the first particle is always the one with larger Pt.
  double dR1=999, dR2=999;
  TLorentzVector v1reco, v2reco, v1gen, v2gen;
  v1reco.SetPtEtaPhiM(dielectron->pt_1, dielectron->eta_1, 
		      dielectron->phi_1, 0.000511);
  v2reco.SetPtEtaPhiM(dielectron->pt_2, dielectron->eta_2, 
		      dielectron->phi_2, 0.000511);
  v1gen .SetPtEtaPhiM(gen->pt_1, gen->eta_1, gen->phi_1, 0.000511);
  v2gen .SetPtEtaPhiM(gen->pt_2, gen->eta_2, gen->phi_2, 0.000511);
  if( dielectron->q_1 < 0 ){
    dR1 = v1reco.DeltaR(v1gen);
    dR2 = v2reco.DeltaR(v2gen);
  }else{
    dR1 = v1reco.DeltaR(v2gen);
    dR2 = v2reco.DeltaR(v1gen);
  }
  // Require that both are within loose dR of 0.4, otherwise bail out
  if( fabs(dR1) > 0.4 || fabs(dR2) > 0.4 ) result = kFALSE; 
  
  return result;
}
*/

// ------------------------------------------------------------

int createSelectionFile(const InputFileMgr_t &inpMgr, 
			EventSelector_t &evtSelector,
			const TString &outSkimFName, 
			DYTools::TRunMode_t runMode) {

  // Event weight handler
  EventWeight_t evWeight;
  evWeight.init(inpMgr.puReweightFlag(),inpMgr.fewzFlag());

#ifdef esfSelectEventsIsObject
  esfSelectEvent_t::Class()->IgnoreTObjectStreamer();
#endif
  esfSelectEvent_t selData;
  //outSkimFName.ReplaceAll(".root","-chk.root");
  TFile *skimFile=new TFile(outSkimFName,"recreate");
  cout << "createSelectionFileName=<" << outSkimFName << ">\n";
  if (!skimFile || !skimFile->IsOpen()) {
    cout << "createSelectionFile: failed to create a file <" 
	 << outSkimFName << ">\n";
    return 0;
  }
  TTree *skimTree= new TTree("Events","Events");
  assert(skimTree);
  selData.createBranches(skimTree);

  // access original n-tuple
  AccessOrigNtuples_t accessInfo;

  double extraWeightFactor=1.0;
  EventCounterExt_t ecTotal("total");

  // Loop over files
  for (unsigned int isample=0; isample<inpMgr.mcSampleCount(); ++isample) {
    //if (isample>0) break;
    const CSample_t *mcSample=inpMgr.mcSampleInfo(isample);
    std::cout << "Processing " << mcSample->getLabel() << "..." << std::endl;
    std::cout << " of size " << mcSample->size() << "\n";
    if (mcSample->size()!=1) {
      std::cout << "mcSample->size is expected to be 1\n";
      return retCodeError;
    }

    for (unsigned int ifile=0; ifile<mcSample->size(); ++ifile) {
      //if (ifile>0) break;

      // Read input file
      TFile infile(mcSample->getFName(ifile),"read");
      assert(infile.IsOpen());
      std::cout << " Reading file <" << mcSample->getFName(ifile) << ">\n";

      // Get the TTrees
      if (!accessInfo.setTree(infile,"Events",true)) {
	return retCodeError;
      }

      // Find weight for events for this file
      // The first file in the list comes with weight 1*extraWeightFactor,
      // all subsequent ones are normalized to xsection and luminosity
      ULong_t maxEvents = accessInfo.getEntries();
      // to match old version package (DYee 7TeV paper), 
      if ((inpMgr.userKeyValueAsInt("USE7TEVMCWEIGHT")==1) && 
	  (isample==0) && (ifile==0)) {
	extraWeightFactor=maxEvents / (inpMgr.totalLumi() * inpMgr.mcSampleInfo(0)->getXsec(ifile));
      }
      //std::cout << "extraWeightFactor=" << extraWeightFactor << ", chk=" << (maxEvents0/inpMgr.mcSampleInfo(0)->getXsec(ifile)) << "\n";
      //const double extraWeightFactor=1.0;
      if (! evWeight.setWeight_and_adjustMaxEvents(maxEvents, inpMgr.totalLumi(), mcSample->getXsec(ifile), 
						   extraWeightFactor, inpMgr.selectEventsFlag())) {
	std::cout << "adjustMaxEvents failed\n";
	return retCodeError;
      }

      std::cout << "       -> sample base weight is " << evWeight.baseWeight() << "\n";
 

      // loop through events
      EventCounterExt_t ec(Form("%s_file%d",mcSample->name.Data(),ifile));
      ec.setIgnoreScale(0); // 1 - count events, 0 - take weight in account
      // adjust the scale in the counter
      // if FEWZ weight should be considered, use evWeight.totalWeight() after
      // the FEWZ weight has been identified (see a line below)
      ec.setScale(evWeight.baseWeight());

      std::cout << "numEntries = " << accessInfo.getEntriesFast() 
		<< ", " << maxEvents << " events will be used" << std::endl;


      for(ULong_t ientry=0; ientry<maxEvents; ientry++) {
	//if (ientry<=290) continue;
	//if (ientry>294) break;
	//if (ientry>100) break;
	if (DYTools::isDebugMode(runMode) && (ientry>1000000)) break; // debug option
	//if (DYTools::isDebugMode(runMode) && (ientry>100)) break; // debug option
	printProgress(500000," ientry=",ientry,maxEvents);
	ec.numEvents_inc();
	
	// Load generator level info
	accessInfo.GetGen(ientry);
	// If the Z->ll leptons are not electrons, discard this event.
	// This is needed for signal MC samples such as Madgraph Z->ll
	// where all 3 lepton flavors are possible
	if (!accessInfo.genLeptonsAreElectrons()) continue;

 	// Load event info
	accessInfo.GetInfoEntry(ientry);

	// Adjust event weight
	// .. here "false" = "not data"
	if (0) {
	  // main analysis way to set weight
	  // nPVs=info->nPU
	  evWeight.set_PU_and_FEWZ_weights(accessInfo,false);
	}
	else {
	  // a different way to calculate nPVs
	  accessInfo.GetPVs(ientry);
	  evWeight.setFewzWeight(accessInfo.genPtr());
	  evWeight.setPUWeight(countGoodVertices(accessInfo.getPVArr()));
	}

	//std::cout << "ientry=" << ientry << ", totalWeight=" << evWeight.totalWeight() << "\n";
	/*
	if (0 || (ientry%20000==0)) {
	  const mithep::TGenInfo *gen=accessInfo.genPtr();
	  std::cout << "ientry=" << ientry << ", gen->vmass=" << gen->vmass << ", gen->vpt=" <<gen->vpt << ", gen->vy=" << gen->vy << "; reco gen->mass=" << gen->mass << ", gen->pt=" << gen->pt << ", gen->y=" << gen->y << "\n";
	  std::cout << "evWeight: " << evWeight << "\n";
	}
	*/


	// adjust the scale in the counter to include FEWZ 
	// (and possibly PU) weight
	//ec.setScale(evWeight.totalWeight());

	// check event trigger
	if (!evtSelector.eventTriggerOk(accessInfo)) {
	  //hFail[ifile]->Fill(gen->mass, fabs(gen->y), evWeight.totalWeight());
	  continue; // no trigger accept? Skip to next event...	
	}
	ec.numEventsPassedEvtTrigger_inc();

	// load dielectron array
	accessInfo.GetDielectrons(ientry);

	const int trace=0; //(ientry==293) ? 1:0;
	if (trace) std::cout << "tracing\n";

	// loop through dielectrons
	//int pass=0;
	for(Int_t i=0; i<accessInfo.dielectronCount(); i++) {
	  mithep::TDielectron *dielectron = accessInfo.editDielectronPtr(i);
	  ec.numDielectrons_inc();
	  if (trace) std::cout << "i=" << i << "\n";

	  // For MC-only, do generator level matching
	  if( ! accessInfo.dielectronMatchedToGenLevel(i) ) {
	    if (trace) printf("failed gen match\n");
	    continue;
	  }
	  ec.numDielectronsGenMatched_inc();

	  // Consider only events in the mass range of interest
	  // Use generator level post-FSR mass.
	  if( accessInfo.genPtr()->mass < DYTools::massBinLimits[0] || 
	      accessInfo.genPtr()->mass > DYTools::massBinLimits[DYTools::nMassBins]) {
	    if (trace) printf("failed mass window\n");
	    continue;
	  }
	  ec.numDielectronsGoodMass_inc();

	  // escale may modify dielectron! But it should not be applied here
	  if (!evtSelector.testDielectron(dielectron,accessInfo.evtInfoPtr(),&ec)) {
	    if (trace) printf("failed selection\n");
	    continue;
	  }

	  if (0) {
	    const mithep::TGenInfo *gen=accessInfo.genPtr();
	    std::cout << "ientry=" << ientry << ", iD=" << i << ", gen->vec(" << gen->vmass << "," << gen->vpt << "," << gen->vy << "), weight=" << evWeight << "\n";
	  }

	  if (int(ec.numDielectronsPass[0]+1e-3)%20000 == 0) {
	    const mithep::TGenInfo *gen=accessInfo.genPtr();
	    std::cout << "ientry=" << ientry << ", iDielectron=" << i << ", gen->mass=" << gen->vmass << ", gen->vy=" <<gen->vy << ", gen->vpt=" << gen->vpt << "; gen->mass=" << gen->mass << ", gen->y=" << gen->y << "\n";
	    std::cout << "evWeight: " << evWeight << "\n";
	  }

	  ec.numDielectronsPass_inc();

	  const int isData=0;
	  selData.assign(accessInfo,isData,i,evWeight.totalWeight());

	  skimTree->Fill();
	} // end loop over dielectrons

	//if (ec.numDielectronsPass[0]>21000) break;

      } // end loop over events
      ec.print();  // print info about file
      //ecSample.add(ec); // accumulate event counts
      ecTotal.add(ec);

      //if (ecTotal.numDielectronsPass[0]>20000) break;
    } // end loop over files
    //ecSample.print(); // print info about sample
    evtSelector.printCounts();
    //if (ecTotal.numDielectronsPass[0]>20000) break;
  } // end loop over samples
  ecTotal.print();
    
  skimFile->cd();
  skimTree->Write();
  skimFile->Close();
  delete skimFile;

  return 1;
}


// ------------------------------------------------------------

// ------------------------------------------------------------

double findEventScaleFactor(int kind, const esfSelectEvent_t &data) {

  return findEventScaleFactor(kind, 
			      data.et_1, data.eta_1,
			      data.et_2, data.eta_2);
}

// ------------------------------------------------------------

double findEventScaleFactor(int kind, 
			    double Et1, double eta1,
			    double Et2, double eta2) {

  const int etBin1 = DYTools::findEtBin(Et1, etBinning);
  const int etaBin1 = DYTools::findEtaBin(eta1, etaBinning);
  const int etBin2 = DYTools::findEtBin(Et2, etBinning);
  const int etaBin2 = DYTools::findEtaBin(eta2, etaBinning);

  const int bin1ok= ((etBin1==-1) || (etBin1==-1)) ? 0:1;
  const int bin2ok= ((etBin2==-1) || (etBin2==-1)) ? 0:1;

  double esf1=1.0;
  double esf2=1.0;
  double esfHLT=1.0; // nonUniversalHLT (double-trigger)

  const int getRECOsf=((kind==DYTools::RECO) || (kind==-1)) ? 1:0;
  const int getIDsf  =((kind==DYTools::ID  ) || (kind==-1)) ? 1:0;
  const int getHLTsf =((kind==DYTools::HLT ) || (kind==-1)) ? 1:0;
  // HLTlegs are not called for nonUniversalHLT directly, when kind==-1
  const int getHLTleg1sf =(kind==DYTools::HLT_leg1) ? 1:0;
  const int getHLTleg2sf =(kind==DYTools::HLT_leg2) ? 1:0;

  if (bin1ok) {
    if (getRECOsf) {
      esf1 *= findScaleFactor(DYTools::RECO, etBin1, etaBin1);
    }
    if (getIDsf) {
      esf1 *= findScaleFactor(DYTools::ID, etBin1, etaBin1);
    }
  }

  if (bin2ok) {
    if (getRECOsf) {
      esf2 *= findScaleFactor(DYTools::RECO, etBin2, etaBin2);
    }
    if (getIDsf) {
      esf2 *= findScaleFactor(DYTools::ID, etBin2, etaBin2);
    }
  }
  
  if ((kind!=-1) && bin1ok) {
    if (getHLTleg1sf) {
      esf1 *= findScaleFactor(DYTools::HLT_leg1, etBin1, etaBin1);
    }
    if (getHLTleg2sf) {
      esf1 *= findScaleFactor(DYTools::HLT_leg2, etBin1, etaBin1);
    }
  }
  if ((kind!=-1) && bin2ok) {
    if (getHLTleg1sf) {
      esf2 *= findScaleFactor(DYTools::HLT_leg1, etBin2, etaBin2);
    }
    if (getHLTleg2sf) {
      esf2 *= findScaleFactor(DYTools::HLT_leg2, etBin2, etaBin2);
    }
  }
  

  if (getHLTsf) {
    if (nonUniversalHLT) {
      if (bin1ok && bin2ok) {
	esfHLT = findScaleFactorHLT(etBin1,etaBin1, Et1,
				    etBin2,etaBin2, Et2);
      }
    }
    else {
      if (bin1ok) {
	esf1 *= findScaleFactor(DYTools::HLT, etBin1, etaBin1);
      }
      if (bin2ok) {
	esf2 *= findScaleFactor(DYTools::HLT, etBin2, etaBin2);
      }
    }
  }

  return esf1*esf2*esfHLT;
}

// --------------------------------------

double findScaleFactor(DYTools::TEfficiencyKind_t kind, int etBin, int etaBin) {

  double result = 0;
  if( (etBin == -1) || (etaBin == -1)) {
    // Found bin outside of calibrated range, return 1.0
    result = 1.0;
    return result;
  }

  double effData=(*dataEff[kind])[etBin][etaBin];
  double effMC  =(*mcEff  [kind])[etBin][etaBin];

  result = effData/effMC;
  //std::cout << "scaleFactor[kind=" << kind << "][etBin=" << etBin << "][etaBin=" << etaBin << "]=" << result << "\n";

  return result;
}

// --------------------------------------

// nonUniversalHLT scale factor
double findScaleFactorHLT(int etBin1, int etaBin1, double et1,
			  int etBin2, int etaBin2, double et2) {

  double result = 0;
  if( (etBin1 == -1) || (etaBin1 == -1) ||
      (etBin2 == -1) || (etaBin2 == -1)) {
    // Found bin outside of calibrated range, return 1.0
    result = 1.0;
    return result;
  }

  double effData=getHLTefficiency(DYTools::DATA,
				  etBin1,etaBin1, et1,
				  etBin2,etaBin2, et2);
  double effMC  =getHLTefficiency(DYTools::MC,
				  etBin1,etaBin1, et1,
				  etBin2,etaBin2, et2);

  result = effData/effMC;
  //std::cout << "nonUniversal scaleFactorHLT[etBin=" << etBin << "][etaBin=" << etaBin << "]=" << result << "\n";

  return result;
}

// --------------------------------------

// nonUniversalHLT efficiency
double getHLTefficiency(DYTools::TDataKind_t dataKind,
			int etBin1, int etaBin1, double et1,
			int etBin2, int etaBin2, double et2) {

  double result = 0;
  if( (etBin1 == -1) || (etaBin1 == -1) ||
      (etBin2 == -1) || (etaBin2 == -1)) {
    // Found bin outside of calibrated range, return 1.0
    result = 1.0;
    return result;
  }

  const int iLeg1=DYTools::HLT_leg1;
  const int iLeg2=DYTools::HLT_leg2;
  TMatrixD *effLeg1=
    (dataKind==DYTools::DATA) ? dataEff[iLeg1] : mcEff[iLeg1];
  TMatrixD *effLeg2=
    (dataKind==DYTools::DATA) ? dataEff[iLeg2] : mcEff[iLeg2];
    
  double eff11=(et1>17.) ? (*effLeg1)[etBin1][etaBin1] : 0.0;
  double eff12=(et1> 8.) ? (*effLeg2)[etBin1][etaBin1] : 0.0;
  double eff21=(et2>17.) ? (*effLeg1)[etBin2][etaBin2] : 0.0;
  double eff22=(et2> 8.) ? (*effLeg2)[etBin2][etaBin2] : 0.0;

  result= eff11*eff22 + eff12*eff21 - eff11*eff21;

  /* // debug branch
  if (0) {
    double oldEff= (*effLeg2)[etBin1][etaBin1] * (*effLeg2)[etBin2][etaBin2];
    double diff=result-oldEff;
    int studyIt=1; //((et1<17.) || (et2<17.)) ? 1:0;
    if (studyIt || (fabs(diff)>0.1)) {
      printf("\nEt1=%4.1lf, etaBin1=%d;  Et2=%4.1lf, etaBin2=%d; hltEff=%5.2lf, hltEff_old=%5.2lf, diff=%5.2lf\n",et1,etaBin1,et2,etaBin2,result,oldEff,result-oldEff);
      printf("  eff11=%6.4lf, eff12=%6.4lf;  eff21=%6.4lf, eff22=%6.4lf\n",eff11,eff12,eff21,eff22);
    }
  }
  */

  return result;
}

// --------------------------------------

// nonUniversalHLT efficiency: What if we ignore HLT leg-dependence?
double getHLTefficiency_IK(DYTools::TDataKind_t dataKind,
			   int etBin1, int etaBin1, double et1,
			   int etBin2, int etaBin2, double et2) {

  double result = 0;
  if( (etBin1 == -1) || (etaBin1 == -1) ||
      (etBin2 == -1) || (etaBin2 == -1)) {
    // Found bin outside of calibrated range, return 1.0
    result = 1.0;
    return result;
  }

  const int iLeg1=DYTools::HLT;  // not leg
  const int iLeg2=DYTools::HLT;  // not leg
  TMatrixD *effLeg1=
    (dataKind==DYTools::DATA) ? dataEff[iLeg1] : mcEff[iLeg1];
  TMatrixD *effLeg2=
    (dataKind==DYTools::DATA) ? dataEff[iLeg2] : mcEff[iLeg2];
    
  double eff11=(et1>17.) ? (*effLeg1)[etBin1][etaBin1] : 0.0;
  double eff12=(et1> 8.) ? (*effLeg2)[etBin1][etaBin1] : 0.0;
  double eff21=(et2>17.) ? (*effLeg1)[etBin2][etaBin2] : 0.0;
  double eff22=(et2> 8.) ? (*effLeg2)[etBin2][etaBin2] : 0.0;

  result= eff11*eff22 + eff12*eff21 - eff11*eff21;

  /* // debug branch
  if (0) {
    double oldEff= (*effLeg2)[etBin1][etaBin1] * (*effLeg2)[etBin2][etaBin2];
    double diff=result-oldEff;
    int studyIt=1; //((et1<17.) || (et2<17.)) ? 1:0;
    if (studyIt || (fabs(diff)>0.1)) {
      printf("\nEt1=%4.1lf, etaBin1=%d;  Et2=%4.1lf, etaBin2=%d; hltEff=%5.2lf, hltEff_old=%5.2lf, diff=%5.2lf\n",et1,etaBin1,et2,etaBin2,result,oldEff,result-oldEff);
      printf("  eff11=%6.4lf, eff12=%6.4lf;  eff21=%6.4lf, eff22=%6.4lf\n",eff11,eff12,eff21,eff22);
    }
  }
  */

  return result;
}

// --------------------------------------

// nonUniversalHLT efficiency
double getHLTefficiencyErr(DYTools::TDataKind_t dataKind,
			   int etBin1, int etaBin1, double et1,
			   int etBin2, int etaBin2, double et2) {

  double result = 0;
  if( (etBin1 == -1) || (etaBin1 == -1) ||
      (etBin2 == -1) || (etaBin2 == -1)) {
    // Found bin outside of calibrated range, return 1.0
    result = 1.0;
    return result;
  }

  const int iLeg1=DYTools::HLT_leg1;
  const int iLeg2=DYTools::HLT_leg2;
    
  TMatrixD *effLeg1=
    (dataKind==DYTools::DATA) ? dataEff[iLeg1] : mcEff[iLeg1];
  TMatrixD *effLeg2=
    (dataKind==DYTools::DATA) ? dataEff[iLeg2] : mcEff[iLeg2];
  TMatrixD *effLeg1Err=
    (dataKind==DYTools::DATA) ? dataEffAvgErr[iLeg1] : mcEffAvgErr[iLeg1];
  TMatrixD *effLeg2Err=
    (dataKind==DYTools::DATA) ? dataEffAvgErr[iLeg2] : mcEffAvgErr[iLeg2];
    
  double eff11=(et1>17.) ? (*effLeg1)[etBin1][etaBin1] : 0.0;
  double eff12=(et1> 8.) ? (*effLeg2)[etBin1][etaBin1] : 0.0;
  double eff21=(et2>17.) ? (*effLeg1)[etBin2][etaBin2] : 0.0;
  double eff22=(et2> 8.) ? (*effLeg2)[etBin2][etaBin2] : 0.0;

  double eff11err=(et1>17.) ? (*effLeg1Err)[etBin1][etaBin1] : 0.0;
  double eff12err=(et1> 8.) ? (*effLeg2Err)[etBin1][etaBin1] : 0.0;
  double eff21err=(et2>17.) ? (*effLeg1Err)[etBin2][etaBin2] : 0.0;
  double eff22err=(et2> 8.) ? (*effLeg2Err)[etBin2][etaBin2] : 0.0;

  // efficiency //result= eff11*eff22 + eff12*eff21 - eff11*eff21;
  result=sqrt(
	      pow( (eff22-eff21)*eff11err, 2) +
	      pow( eff21*eff12err, 2) +
	      pow( (eff12-eff11)*eff21err, 2) +
	      pow( eff11*eff22err, 2)
	      );

  return result;
}

// ---------------------- all scale factors smeared ------------------------------

// ------------------------------------------------------------

double findEventScaleFactorSmeared(int kind, const esfSelectEvent_t &data,
				   int iexp) {
  return findEventScaleFactorSmeared(kind,
				     data.et_1, data.eta_1,
				     data.et_2, data.eta_2,
				     iexp);
}

// ------------------------------------------------------------

double findEventScaleFactorSmeared(int kind, 
				   double Et1, double eta1,
				   double Et2, double eta2,
				   int iexp) {

  const int debug=0; //(iexp==0) ? 1 : 0;
  if (debug) {
    std::cout << "findEventScaleFactorSmeared("
	      << EfficiencyKindName(DYTools::TEfficiencyKind_t(kind)) 
	      << ", " << Et1 << ", " << eta1 << "; "
	      << Et2 << ", " << eta2 << "; iexp=" << iexp << ")\n";
  }

  double esf1=1.0;
  double esf2=1.0;
  double esfHLT=1.0;
  int etBin1 = DYTools::findEtBin(Et1, etBinning);
  int etaBin1 = DYTools::findEtaBin(eta1, etaBinning);
  int etBin2 = DYTools::findEtBin(Et2, etBinning);
  int etaBin2 = DYTools::findEtaBin(eta2, etaBinning);
  if (debug) std::cout << Form("bins: (%d,%d; %d,%d)\n",etBin1,etaBin1,etBin2,etaBin2) << "\n";

  int bin1ok=0;
  int bin2ok=0;

  int getRECOsf=((kind==DYTools::RECO) || (kind==-1)) ? 1:0;
  int getIDsf  =((kind==DYTools::ID  ) || (kind==-1)) ? 1:0;
  int getHLTsf =((kind==DYTools::HLT ) || (kind==-1)) ? 1:0;
  // HLTlegs are not called for nonUniversalHLT directly, when kind==-1
  const int getHLTleg1sf =(kind==DYTools::HLT_leg1) ? 1:0;
  const int getHLTleg2sf =(kind==DYTools::HLT_leg2) ? 1:0;


  if ((etBin1!=-1) && (etaBin1!=-1)) {
    bin1ok=1;
    if (getRECOsf) {
      if (debug) std::cout << "getRECOsf(1)\n";
      esf1*= findScaleFactorSmeared(DYTools::RECO, etBin1, etaBin1, iexp);
    }
    if (getIDsf) {
      if (debug) std::cout << "getIDsf(1)\n";
      esf1*= findScaleFactorSmeared(DYTools::ID, etBin1, etaBin1, iexp);
    }
  }

  if ((etBin2!=-1) && (etaBin2!=-1)) {
    bin2ok=1;
    if (getRECOsf) {
      if (debug) std::cout << "getRECOsf(2)\n";
      esf2*= findScaleFactorSmeared(DYTools::RECO, etBin2, etaBin2, iexp);
    }
    if (getIDsf) {
      if (debug) std::cout << "getIDsf(2)\n";
      esf2*= findScaleFactorSmeared(DYTools::ID, etBin2, etaBin2, iexp);
    }
  }

  if ((kind!=-1) && bin1ok) {
    if (getHLTleg1sf) {
      if (debug) std::cout << "getHLTleg1sf(1)\n";
      esf1 *= findScaleFactorSmeared(DYTools::HLT_leg1, etBin1, etaBin1, iexp);
    }
    if (getHLTleg2sf) {
      if (debug) std::cout << "getHLTleg2sf(1)\n";
      esf1 *= findScaleFactorSmeared(DYTools::HLT_leg2, etBin1, etaBin1, iexp);
    }
  }
  if ((kind!=-1) && bin2ok) {
    if (getHLTleg1sf) {
      if (debug) std::cout << "getHLTleg1sf(2)\n";
      esf2 *= findScaleFactorSmeared(DYTools::HLT_leg1, etBin2, etaBin2, iexp);
    }
    if (getHLTleg2sf) {
      if (debug) std::cout << "getHLTleg2sf(2)\n";
      esf2 *= findScaleFactorSmeared(DYTools::HLT_leg2, etBin2, etaBin2, iexp);
    }
  }

  if (getHLTsf) {
    if (nonUniversalHLT) {
      if (bin1ok && bin2ok) {
	if (debug) std::cout << "getHLTsf nonUni\n";
	esfHLT = findScaleFactorHLTSmeared(etBin1,etaBin1, Et1,
					   etBin2,etaBin2, Et2, iexp);
      }
    }
    else {
      if (bin1ok) {
	if (debug) std::cout << "getHLT sf(1)\n";
	esf1 *= findScaleFactorSmeared(DYTools::HLT, etBin1, etaBin1, iexp);
      }
      if (bin2ok) {
	if (debug) std::cout << "getHLT sf(2)\n";
	esf2 *= findScaleFactorSmeared(DYTools::HLT, etBin2, etaBin2, iexp);
      }
    }
  }

  return esf1*esf2*esfHLT;
}

// --------------------------------------

double findScaleFactorSmeared(
       DYTools::TEfficiencyKind_t kind, 
       int scEtBin, int scEtaBin, 
       int iexp) {

  return findScaleFactorSmeared(kind,scEtBin,scEtaBin,ro_Data[iexp],ro_MC[iexp]);

}

// --------------------------------------

double findScaleFactorSmeared(
       DYTools::TEfficiencyKind_t kind, 
       int etBin, int etaBin, 
       const EffArray_t &dataRndWeight, const EffArray_t &mcRndWeight) {

  double result = 0;
  if( (etBin == -1) || (etaBin == -1)) {
    std::cout << "err: findScaleFactorSmeared(etBin=" << etBin << ", etaBin=" << etaBin << ")\n";
    // Found bin outside of calibrated range, return 1.0
    result = 1.0;
    return result;
  }

  double effData=
    (*dataEff[kind])[etBin][etaBin] + 
    dataRndWeight[kind][etBin][etaBin] * (*dataEffAvgErr[kind])[etBin][etaBin];
  //if (effData>100.) effData=100.;
  double effMC=
    (*mcEff[kind])[etBin][etaBin] + 
    mcRndWeight[kind][etBin][etaBin] * (*mcEffAvgErr[kind])[etBin][etaBin];
  //if (effMC>100.) effMC=100.;
  //std::cout << "scaleFactorSmeared[kind=" << kind << "][etBin=" << etBin << "][etaBin=" << etaBin << "]=" << effData << '/' << effMC << "=" << (effData/effMC) << "\n";
  return (effData/effMC);
}

// --------------------------------------

// nonUniversalHLT
double findScaleFactorHLTSmeared(int scEtBin1, int scEtaBin1, double scEt1,
				 int scEtBin2, int scEtaBin2, double scEt2,
				 int iexp) {

  return findScaleFactorHLTSmeared(scEtBin1,scEtaBin1, scEt1,
				   scEtBin2,scEtaBin2, scEt2,
				   ro_Data[iexp],ro_MC[iexp]);

}

// --------------------------------------

double findScaleFactorHLTSmeared(int etBin1, int etaBin1, double Et1,
				 int etBin2, int etaBin2, double Et2,
	   const EffArray_t &dataRndWeight, const EffArray_t &mcRndWeight) {

  double result = 0;
  if( (etBin1 == -1) || (etaBin1 == -1) ||
      (etBin2 == -1) || (etaBin2 == -1)) {
    // Found bin outside of calibrated range, return 1.0
    std::cout << "err: findScaleFactorHLTSmeared(etBin1=" << etBin1 << ", etaBin1=" << etaBin1 << "; etBin2=" << etBin1 << ", etaBin2=" << etaBin2 << ")\n";

    result = 1.0;
    return result;
  }

  double effData=getHLTefficiencySmeared(DYTools::DATA,
					 etBin1,etaBin1, Et1,
					 etBin2,etaBin2, Et2,
					 dataRndWeight);
  double effMC  =getHLTefficiencySmeared(DYTools::MC,
					 etBin1,etaBin1, Et1,
					 etBin2,etaBin2, Et2,
					 mcRndWeight);

  if (effMC==double(0.)) {
    std::cout << "err: findScaleFactorHLTSmeared  got effMC=0 !\n";
  }
  result = effData/effMC;
  //std::cout << "scaleFactorHLTsmeared[etBin=" << etBin << "][etaBin=" << etaBin << "]=" << result << "\n";

  return result;
}

// --------------------------------------

// nonUniversalHLT efficiency
double getHLTefficiencySmeared(DYTools::TDataKind_t dataKind,
			       int etBin1, int etaBin1, double Et1,
			       int etBin2, int etaBin2, double Et2,
			       const EffArray_t &rndWeight) {

  double result = 0;
  if( (etBin1 == -1) || (etaBin1 == -1) ||
      (etBin2 == -1) || (etaBin2 == -1)) {
    // Found bin outside of calibrated range, return 1.0
    result = 1.0;
    return result;
  }

  const int iLeg1=DYTools::HLT_leg1;
  const int iLeg2=DYTools::HLT_leg2;
  TMatrixD *effLeg1=
    (dataKind==DYTools::DATA) ? dataEff[iLeg1] : mcEff[iLeg1];
  TMatrixD *effLeg1Err=
    (dataKind==DYTools::DATA) ? dataEffAvgErr[iLeg1] : mcEffAvgErr[iLeg1];
  TMatrixD *effLeg2=
    (dataKind==DYTools::DATA) ? dataEff[iLeg2] : mcEff[iLeg2];
  TMatrixD *effLeg2Err=
    (dataKind==DYTools::DATA) ? dataEffAvgErr[iLeg2] : mcEffAvgErr[iLeg2];

  double eff11=0., eff12=0., eff21=0., eff22=0.;
 
  if (Et1>17.) {
    eff11 = (*effLeg1)[etBin1][etaBin1] +
      rndWeight[iLeg1][etBin1][etaBin1] * (*effLeg1Err)[etBin1][etaBin1];
  }
  if (Et1> 8.) {
    eff12 = (*effLeg2)[etBin1][etaBin1] +
      rndWeight[iLeg2][etBin1][etaBin1] * (*effLeg2Err)[etBin1][etaBin1];
  }
  if (Et2>17.) {
    eff21 = (*effLeg1)[etBin2][etaBin2] +
      rndWeight[iLeg1][etBin2][etaBin2] * (*effLeg1Err)[etBin2][etaBin2];
  }
  if (Et2> 8.) {
    eff22 = (*effLeg2)[etBin2][etaBin2] +
      rndWeight[iLeg2][etBin2][etaBin2] * (*effLeg2Err)[etBin2][etaBin2];
  }

  result= eff11*eff22 + eff12*eff21 - eff11*eff21;

  return result;

}

// --------------------------------------
// --------------------------------------


 void drawEfficiencyGraphs(TGraphErrors *grData, TGraphErrors *grMc,
			   TString yAxisTitle, TString text, TString plotName,
			   TFile *fRoot){
   
   // Generate "random" canvas name
//    TTimeStamp time;
//    TString c = "c";
//    c += time.GetNanoSec();
   TString c = plotName;
//    printf("Canvas name %s\n", c.Data());

   TCanvas *c2 = MakeCanvas(c,c);
   CPlot plot1(c,"","E_{T} [GeV]", yAxisTitle);
   plot1.SetLogx(); 
   plot1.AddGraph(grData,"data","PE2", kBlue);
   plot1.AddGraph(grMc  ,"MC"  ,"PE", kBlack);
   plot1.Draw(c2);
   plot1.SetYRange(0.0,1.1);
   grData->GetXaxis()->SetTitle("E_{T} [GeV]");
   grData->GetXaxis()->SetMoreLogLabels();
   grData->GetXaxis()->SetNoExponent();
   grData->SetFillStyle(3001);
   grData->SetFillColor(kBlue);
   grMc->SetMarkerStyle(24);
   plot1.TransLegend(0.0, -0.5);
   plot1.AddTextBox(text, 0.6,0.4,0.87,0.55, 0);

   TLine *line = new TLine(10,1.0,500,1.0);
   line->SetLineStyle(kDashed);
   line->Draw("same");

   plot1.Draw(c2, savePlots, "png");
   if (fRoot) c2->Write();

   return;
 }

// -------------------------------------------------------------------------
/*
template<class Graph_t>
void drawEfficiencyGraphsAsymmErrors(Graph_t *grData, Graph_t *grMc,
                     TString yAxisTitle, TString text, TString plotName){
  
   // Generate "random" canvas name
//    TTimeStamp time;
//    TString c = "c";
//    c += time.GetNanoSec();
   TString c = plotName;
//    printf("Canvas name %s\n", c.Data());

   TCanvas *c2 = MakeCanvas(c,c);
   CPlot plot1(c,"","E_{T} [GeV]", yAxisTitle);
   plot1.SetLogx(); 
   plot1.AddGraph(grData,"data","PE2", kBlue);
   plot1.AddGraph(grMc  ,"MC"  ,"PE", kBlack);
   plot1.Draw(c2);
   plot1.SetYRange(0.2,1.1);
   grData->GetXaxis()->SetTitle("E_{T} [GeV]");
   grData->GetXaxis()->SetMoreLogLabels();
   grData->GetXaxis()->SetNoExponent();
   grData->SetFillStyle(3001);
   grData->SetFillColor(kBlue);
   grMc->SetMarkerStyle(24);
   plot1.TransLegend(0.0, -0.5);
   plot1.AddTextBox(text, 0.6,0.4,0.8,0.5, 0);
   plot1.Draw(c2, savePlots, "png");

   TLine *line = new TLine(10,1.0,500,1.0);
   line->SetLineStyle(kDashed);
   line->Draw("same");

   return;
 }

// -------------------------------------------------------------------------

template<class Graph_t>
void drawEfficiencyGraphsAsymmErrorsPU(Graph_t *grData, Graph_t *grMc,
				     const TString &yAxisTitle, 
				     const TString &text, 
				     const TString &plotName){
  
   // Generate "random" canvas name
//    TTimeStamp time;
//    TString c = "c";
//    c += time.GetNanoSec();
   TString c = plotName;
//    printf("Canvas name %s\n", c.Data());

   TString xAxisTitle = "nGoodPV";
   TCanvas *c2 = MakeCanvas(c,c);
   CPlot plot1(c,"",xAxisTitle, yAxisTitle);
   //plot1.SetLogx(); 
   plot1.AddGraph(grData,"data","PE2", kBlue);
   plot1.AddGraph(grMc  ,"MC"  ,"PE", kBlack);
   plot1.Draw(c2);
   plot1.SetYRange(0.2,1.1);
   grData->GetXaxis()->SetTitle(xAxisTitle);
   //grData->GetXaxis()->SetMoreLogLabels();
   //grData->GetXaxis()->SetNoExponent();
   grData->SetFillStyle(3001);
   grData->SetFillColor(kBlue);
   grMc->SetMarkerStyle(24);
   plot1.TransLegend(0.0, -0.5);
   plot1.AddTextBox(text, 0.6,0.4,0.8,0.5, 0);
   plot1.Draw(c2, savePlots, "png");

   TLine *line = new TLine(0,1.0,DYTools::nPVLimits[DYTools::nPVBinCount],1.0);
   line->SetLineStyle(kDashed);
   line->Draw("same");

   return;
 }
*/

// -------------------------------------------------------------------------

void drawScaleFactorGraphs(TGraphErrors *gr, TString yAxisTitle, TString text,
			   TString plotName, TFile *fRoot){
   
   // Generate "random" canvas name
//    TTimeStamp time;
//    TString c = "c";
//    c += time.GetNanoSec();
  TString c = plotName;
//   printf("Canvas name %s\n", c.Data());

   TCanvas *c2 = MakeCanvas(c,c);
   CPlot plot1(c,"","E_{T} [GeV]", yAxisTitle);
   plot1.SetLogx(); 
   plot1.AddGraph(gr,"present E_{T}/#eta dep.","PE", kBlack);
   plot1.Draw(c2);
   plot1.SetYRange(0.5,1.5);
   gr->GetXaxis()->SetTitle("E_{T} [GeV]");
   gr->GetXaxis()->SetMoreLogLabels();
   gr->GetXaxis()->SetNoExponent();
   plot1.AddTextBox(text, 0.6,0.35,0.87,0.5, 0);
   
   TLine *line = new TLine(10,1.0,500,1.0);
   line->SetLineStyle(kDashed);
   line->Draw("same");

   plot1.Draw(c2,savePlots, "png");
   if (fRoot) c2->Write();
   return;
 }

// -------------------------------------------------------------------------

void drawEventScaleFactorGraphs(TGraphErrors *gr, TString yAxisTitle, 
				TString plotName, TFile *fRoot) {
  
  // Generate "random" canvas name
//   TTimeStamp time;
//   TString c = "c";
//   c += time.GetNanoSec();
  TString c = plotName;
//   printf("Canvas name %s\n", c.Data());
  
  TCanvas *c2 = MakeCanvas(c,c);
  CPlot plot1(c,"","m(e^{+}e^{-}) [GeV]", yAxisTitle);
  plot1.SetLogx(); 
  plot1.AddGraph(gr,"present E_{T}/#eta dep.","PE", kBlack);
  plot1.Draw(c2);
  plot1.SetYRange(0.0,1.5);
  gr->GetXaxis()->SetTitle("m(e^{+}e^{-}) [GeV]");
  gr->GetXaxis()->SetMoreLogLabels();
  gr->GetXaxis()->SetNoExponent();
//   cout << "From main progam CPlot::sOutDir " << CPlot::sOutDir << endl;
  //CPlot::sOutDir = "plots";
  
  TLine *line = new TLine(15,1.0,500,1.0);
  line->SetLineStyle(kDashed);
  line->Draw("same");

  plot1.Draw(c2, savePlots, "png");
  if (fRoot) c2->Write();
  return;
 }

// -------------------------------------------------------------------------

void drawEfficiencies(TFile *fRoot){
  // Make graphs
  double x[etBinCount];
  double dx[etBinCount];
  for(int i=0; i<etBinCount; i++){
    x[i]  = 0.5*(etBinLimits[i] + etBinLimits[i+1]);
    dx[i] = 0.5*(etBinLimits[i+1] - etBinLimits[i]);
  }


  double effData[etBinCount],effDataErr[etBinCount];
  double effMC[etBinCount],effMCErr[etBinCount];
  const int bufsize=30;
  char plotLabel[bufsize];

  int signedEta=DYTools::signedEtaBinning(etaBinning);
  const char *etaSignStr=(signedEta) ? "_eta" : "_abs_eta";
  const char *etaLabelStr=(signedEta) ? "#eta" : "|#eta|";

  for (int kind=0; kind<NEffTypes; ++kind) {
    if ((nonUniversalHLT) && (kind==DYTools::HLT)) continue; // skip for now
    for (int iEta=0; iEta<etaBinCount; ++iEta) {
      TString etaStr=
	Form("%s_%5.3lf__%5.3lf", etaSignStr,
	     etaBinLimits[iEta],etaBinLimits[iEta+1]);
      etaStr.ReplaceAll(".","_");
      snprintf(plotLabel,bufsize,"%5.3lf < %s < %5.3lf",
	       etaBinLimits[iEta],etaLabelStr,etaBinLimits[iEta+1]);

      for (int iEt=0; iEt<etBinCount; ++iEt) {
	effData[iEt]= (*dataEff[kind])[iEt][iEta];
	effDataErr[iEt]= (*dataEffAvgErr[kind])[iEt][iEta];
      }
      for (int iEt=0; iEt<etBinCount; ++iEt) {
	effMC[iEt]= (*mcEff[kind])[iEt][iEta];
	effMCErr[iEt]= (*mcEffAvgErr[kind])[iEt][iEta];
      }
      TGraphErrors *grDataEff 
	= new TGraphErrors(etBinCount, x,  effData, dx, effDataErr);
  
      TGraphErrors *grMcEff 
	= new TGraphErrors(etBinCount, x,  effMC, dx, effMCErr);
  
      // Draw all graphs
      TString effName=EfficiencyKindName(DYTools::TEfficiencyKind_t(kind));
      TString ylabel=TString("efficiency_{") + effName + TString("}");
      TString plotName = TString("plot_eff_") + DYTools::analysisTag + TString("_") +
	effName + etaStr;
      
      drawEfficiencyGraphs(grDataEff, grMcEff,
			   ylabel, plotLabel, plotName, fRoot);

      //delete grMcEff;
      //delete grDataEff;
    }
  }

  if (nonUniversalHLT) {
    // create a special range variable
    // First, collect bin edges
    std::vector<double> tmpEtBinLimits;
    tmpEtBinLimits.reserve(20);
    tmpEtBinLimits.push_back(etBinLimits[0]);
    for (int iEt=1; iEt<etBinCount; ++iEt) {
      if ((etBinLimits[iEt-1]<8.) && (etBinLimits[iEt]>8.)) {
	tmpEtBinLimits.push_back(8.);
      }
      else if ((etBinLimits[iEt-1]<17.) && (etBinLimits[iEt]>17.)) {
	tmpEtBinLimits.push_back(17.);
      }
      tmpEtBinLimits.push_back(etBinLimits[iEt]);
    }
    tmpEtBinLimits.push_back(etBinLimits[etBinCount]);

    // Second, create a special range variable
    const unsigned int specEtBinCount=tmpEtBinLimits.size()-1;
    double *tmpBinCs = new double[specEtBinCount];
    double *dTmpBins= new double[specEtBinCount];
    for (unsigned int i=0; i<specEtBinCount; ++i) {
      double tmin=tmpEtBinLimits[i];
      double tmax=tmpEtBinLimits[i+1];
      tmpBinCs [i]=0.5*(tmin+tmax);
      dTmpBins[i]=0.5*(tmax-tmin);
    }

    if (0) {
      std::cout << "tmpEtBinLimits[" << tmpEtBinLimits.size() << "]: ";
      for (unsigned int i=0; i<tmpEtBinLimits.size(); ++i) {
	std::cout << " " << tmpEtBinLimits[i];
      }
      std::cout << "\n";
      std::cout << "tmpBinCs[" << specEtBinCount << "]: ";
      for (unsigned int i=0; i<specEtBinCount; ++i) {
	std::cout << " " << tmpBinCs[i];
      }
      std::cout << "\n";
      HERE("test");
    }

    // Now, create the plots
    for (int iEta=0; iEta<etaBinCount; ++iEta) {
      TString etaStr=
	Form("%s_%5.3lf__%5.3lf", etaSignStr,
	     etaBinLimits[iEta],etaBinLimits[iEta+1]);
      etaStr.ReplaceAll(".","_");
      snprintf(plotLabel,bufsize,"%5.3lf < %s < %5.3lf",
	       etaBinLimits[iEta],etaLabelStr,etaBinLimits[iEta+1]);

      for (unsigned int iEt=0; iEt<specEtBinCount; ++iEt) {
	int iEt_tmp= DYTools::findEtBin( tmpBinCs[iEt], etBinning);
	effData[iEt]= getHLTefficiency(DYTools::DATA,
				       iEt_tmp, iEta, tmpBinCs[iEt],
				       iEt_tmp, iEta, tmpBinCs[iEt]);
	effDataErr[iEt]= getHLTefficiencyErr(DYTools::DATA,
					     iEt_tmp, iEta, tmpBinCs[iEt],
					     iEt_tmp, iEta, tmpBinCs[iEt]);
      }

      for (unsigned int iEt=0; iEt<specEtBinCount; ++iEt) {
	int iEt_tmp= DYTools::findEtBin( tmpBinCs[iEt], etBinning);
	effMC[iEt]= getHLTefficiency(DYTools::MC,
				     iEt_tmp, iEta, tmpBinCs[iEt],
				     iEt_tmp, iEta, tmpBinCs[iEt]);
	effMCErr[iEt]= getHLTefficiencyErr(DYTools::MC,
					   iEt_tmp, iEta, tmpBinCs[iEt],
					   iEt_tmp, iEta, tmpBinCs[iEt]);
      }
      TGraphErrors *grDataEff 
	= new TGraphErrors(specEtBinCount, tmpBinCs,  effData, dTmpBins, effDataErr);
  
      TGraphErrors *grMcEff 
	= new TGraphErrors(specEtBinCount, tmpBinCs,  effMC, dTmpBins, effMCErr);
  
      // Draw all graphs
      TString effName=EfficiencyKindName(DYTools::HLT);
      TString ylabel=TString("efficiency_{") + effName + TString(";asym,diag}");
      TString plotName = TString("plot_eff_") + DYTools::analysisTag + TString("_") +
	effName + etaStr;
      
      drawEfficiencyGraphs(grDataEff, grMcEff,
			   ylabel, plotLabel, plotName, fRoot);

   
    }
    if (tmpBinCs) delete tmpBinCs;
    if (dTmpBins) delete dTmpBins;
  }

  return;
}

// -------------------------------------------------------------------------

void drawScaleFactors(TFile *fRoot){

  double x[etBinCount];
  double dx[etBinCount];

  for(int i=0; i<etBinCount; i++){
    x[i]  = (etBinLimits[i] + etBinLimits[i+1])/2.0;
    dx[i] = (etBinLimits[i+1] - etBinLimits[i])/2.0;
  }


  double scale[etBinCount];
  double scaleErr[etBinCount];
  const int bufsize=30;
  char plotLabel[bufsize];

  int signedEta=DYTools::signedEtaBinning(etaBinning);
  const char *etaSignStr=(signedEta) ? "_eta" : "_abs_eta";
  const char *etaLabelStr=(signedEta) ? "#eta" : "|#eta|";

  for (int kind=0; kind<NEffTypes; ++kind) {
    for (int iEta=0; iEta<etaBinCount; ++iEta) {
      TString etaStr=
	Form("%s_%5.3lf__%5.3lf", etaSignStr,
	     etaBinLimits[iEta],etaBinLimits[iEta+1]);
      etaStr.ReplaceAll(".","_");
      snprintf(plotLabel,bufsize,"%5.3lf < %s < %5.3lf",
	       etaBinLimits[iEta],etaLabelStr,etaBinLimits[iEta+1]);

      for (int iEt=0; iEt<etBinCount; ++iEt) {
	scale[iEt]= (*dataEff[kind])[iEt][iEta] / (*mcEff[kind])[iEt][iEta];
	scaleErr[iEt]= errOnRatio( (*dataEff[kind])[iEt][iEta], 
				   (*dataEffAvgErr[kind])[iEt][iEta],
				   (*mcEff[kind])[iEt][iEta], 
				   (*mcEffAvgErr[kind])[iEt][iEta] );
      }

      TGraphErrors *grScaleFactor
	= new TGraphErrors(etBinCount, x, scale, dx, scaleErr);
      TString effName=EfficiencyKindName(DYTools::TEfficiencyKind_t(kind));
      TString ylabel=TString("scale factor ") + effName;
      TString plotName = TString("plot_scale_") + DYTools::analysisTag + TString("_") +
	effName + etaStr;

      drawScaleFactorGraphs(grScaleFactor, ylabel, plotLabel, plotName, fRoot);
      //delete grScaleFactor;
    }
  }
}

// -------------------------------------------------------------------------

double errOnRatio(double a, double da, double b, double db){

  double result = 0;
  if(a == 0 || b == 0)
    return result;
 
  result = (a/b)*sqrt( (da/a)*(da/a) + (db/b)*(db/b) );
  return result;
}

// -------------------------------------------------------------------------

TGraphErrors* createGraph_vsMass(const TVectorD &v, const TVectorD &vErr) {
  if (v.GetNoElements()!=DYTools::nMassBins) {
    printf("createGraph_vsMass works with vectors of size DYTools::nMassBins=%d, not with size %d\n",DYTools::nMassBins,v.GetNoElements());
    assert(0);
  }

  // repackage into arrays
  double x[DYTools::nMassBins];
  double dx[DYTools::nMassBins];
  double y   [DYTools::nMassBins];
  double yErr [DYTools::nMassBins];

  for(int i=0; i<DYTools::nMassBins; i++){
    x[i] = (DYTools::massBinLimits[i] + DYTools::massBinLimits[i+1])/2.0;
    dx[i]= (DYTools::massBinLimits[i+1] - DYTools::massBinLimits[i])/2.0;
    y[i] = v[i];
    yErr[i]= vErr[i];
  }

  TGraphErrors *gr = 
    new TGraphErrors(DYTools::nMassBins, x, y, dx, yErr);
  return gr;
}

// -------------------------------------------------------------------------

void drawEventScaleFactors(const TVectorD &scaleRecoV, const TVectorD &scaleRecoErrV,
			   const TVectorD &scaleIdV , const TVectorD &scaleIdErrV ,
			   const TVectorD &scaleHltV, const TVectorD &scaleHltErrV,
			   const TVectorD &scaleV   , const TVectorD &scaleErrV   ,
			   TFile *fRoot)
{

  if (scaleRecoV.GetNoElements()!=DYTools::nMassBins) {
    std::cout << "\n\nERROR: drawEventScaleFactors should be called for "
	      << "mass-bin indexed vectors\n\n";
    return;
  }

  TGraphErrors *grScale = createGraph_vsMass(scaleV,scaleErrV);
  TGraphErrors *grScaleReco = createGraph_vsMass(scaleRecoV,scaleRecoErrV);
  TGraphErrors *grScaleId  = createGraph_vsMass(scaleIdV,scaleIdErrV);
  TGraphErrors *grScaleHlt = createGraph_vsMass(scaleHltV,scaleHltErrV);

  TString plotName;
  TString plotNameBase = 
    TString("plot_event_scale_") + DYTools::analysisTag + TString("_");
  plotName = plotNameBase + TString("reco");
  drawEventScaleFactorGraphs(grScaleReco, "RECO scale factor", plotName,fRoot);
  //plotName = "plot_event_scale_id";
  plotName = plotNameBase + TString("id");
  drawEventScaleFactorGraphs(grScaleId , "ID scale factor", plotName, fRoot);
  plotName = plotNameBase + TString("hlt");
  drawEventScaleFactorGraphs(grScaleHlt, "HLT scale factor", plotName, fRoot);
  plotName = plotNameBase + TString("full");
  drawEventScaleFactorGraphs(grScale   , "event scale factor", plotName,fRoot);

}

// -------------------------------------------------------------------------

void drawEventScaleFactorsHLT(const TVectorD &scaleHltLeg1V, const TVectorD &scaleHltLeg1ErrV,
			   const TVectorD &scaleHltLeg2V , const TVectorD &scaleHltLeg2ErrV ,
			   TFile *fRoot)
{

  if (scaleHltLeg1V.GetNoElements()!=DYTools::nMassBins) {
    std::cout << "\n\nERROR: drawEventScaleFactorsHLT should be called for "
	      << "mass-bin indexed vectors\n\n";
    return;
  }

  TGraphErrors *grScaleHltLeg1 = createGraph_vsMass(scaleHltLeg1V,scaleHltLeg1ErrV);
  TGraphErrors *grScaleHltLeg2 = createGraph_vsMass(scaleHltLeg2V,scaleHltLeg2ErrV);
  if (!grScaleHltLeg1 || !grScaleHltLeg2) {
    std::cout << "drawEventScaleFactorsHLT: got null ptr(s) to graphs\n";
    return;
  }

  TString plotName;
  TString plotNameBase = 
    TString("plot_event_scale_") + DYTools::analysisTag + TString("_");
  plotName = plotNameBase + TString("hltLeg1");
  drawEventScaleFactorGraphs(grScaleHltLeg1, "HLT_leg1 scale factor", plotName,fRoot);
  plotName = plotNameBase + TString("hltLeg2");
  drawEventScaleFactorGraphs(grScaleHltLeg2 , "HLT_leg2 scale factor", plotName, fRoot);

}

// -------------------------------------------------------------------------

TGraphErrors* createGraph_vsMassFI(const TVectorD &v, const TVectorD &vErr,
				 int rapidityIndex) {
  int nUnfoldingBins = DYTools::getTotalNumberOfBins();
  if (v.GetNoElements()!=nUnfoldingBins) {
    std::cout << "\n\nERROR: createGraph_vsMassFI should be called for "
	      << "flat-indexed vectors\n\n";
    assert(0);
  }

  // repackage into arrays
  double x[DYTools::nMassBins];
  double dx[DYTools::nMassBins];
  double y   [DYTools::nMassBins];
  double yErr [DYTools::nMassBins];

  double rapidity=DYTools::findAbsYValue(0,rapidityIndex);

  for(int i=0; i<DYTools::nMassBins; i++){
    std::cout << "i=" << i << ", rapidityIndex=" << rapidityIndex << "\n";
    std::cout << "nYBins[massBin]=" << DYTools::nYBins[i] << "\n";
    int idx=DYTools::findIndexFlat(i,rapidityIndex);
    if (i==DYTools::nMassBins-1) {
      // last mass bin is special
      std::cout << "rapidity=" << rapidity << ", findAbsYBin=" << DYTools::findAbsYBin(i,rapidity) << "\n";
      idx=DYTools::findIndexFlat(i, DYTools::findAbsYBin(i,rapidity));
    }
    if (idx<0) {
      std::cout <<"WARNING: createGraph_vsMassFI: massBin=" << i << ", rapidityIndex=" << rapidityIndex << "(supplied to subroutine), idx=" << idx << std::endl;
      //assert(0);
      return NULL;
    }

    x[i] = (DYTools::massBinLimits[i] + DYTools::massBinLimits[i+1])/2.0;
    dx[i]= (DYTools::massBinLimits[i+1] - DYTools::massBinLimits[i])/2.0;
    y[i] = v[idx];
    yErr[i]= vErr[idx];
  }

  TGraphErrors *gr = 
    new TGraphErrors(DYTools::nMassBins, x, y, dx, yErr);
  return gr;
}

// -------------------------------------------------------------------------

void drawEventScaleFactorsFI(const TVectorD &scaleRecoFIV, const TVectorD &scaleRecoErrFIV,
			     const TVectorD &scaleIdFIV , const TVectorD &scaleIdErrFIV ,
			     const TVectorD &scaleHltFIV, const TVectorD &scaleHltErrFIV,
			     const TVectorD &scaleFIV   , const TVectorD &scaleErrFIV   ,
			     int rapidityIndex,
			     TFile *fRoot,
			     std::vector<CPlot*> *cplotV)
{
  int nUnfoldingBins = DYTools::getTotalNumberOfBins();
  if (scaleRecoFIV.GetNoElements()!=nUnfoldingBins) {
    std::cout << "\n\nERROR: drawEventScaleFactorsFI should be called for "
	      << "flat-indexed vectors\n\n";
    return;
  }


  TGraphErrors *grScale = createGraph_vsMassFI(scaleFIV,scaleErrFIV,rapidityIndex);
  TGraphErrors *grScaleReco = createGraph_vsMassFI(scaleRecoFIV,scaleRecoErrFIV,rapidityIndex);
  TGraphErrors *grScaleId  = createGraph_vsMassFI(scaleIdFIV,scaleIdErrFIV,rapidityIndex);
  TGraphErrors *grScaleHlt = createGraph_vsMassFI(scaleHltFIV,scaleHltErrFIV,rapidityIndex);

  char buf[20];
  double rapidity=DYTools::findAbsYValue(0,rapidityIndex);
  sprintf(buf,"_y=%4.2lf_",rapidity); // for file name
  TString yStr=buf;
  yStr.ReplaceAll(".","_"); yStr.ReplaceAll("=","_");
  sprintf(buf,"|y|=%4.2lf",rapidity); // for label
  TString plotName;
  TString plotNameBase = 
    TString("plot_event_scale_") + DYTools::analysisTag + TString("_") + yStr;
  plotName = plotNameBase + TString("reco");
  drawEventScaleFactorGraphs(grScaleReco, 
			     TString("RECO scale factor ") + TString(buf),
			     plotName, fRoot);
  plotName = plotNameBase + TString("id");
  drawEventScaleFactorGraphs(grScaleId , 
			     TString("ID scale factor ") + TString(buf),
			     plotName, fRoot);
  plotName = plotNameBase + TString("hlt");
  drawEventScaleFactorGraphs(grScaleHlt, 
			     TString("HLT scale factor ") + TString(buf), 
			     plotName, fRoot);
  plotName = plotNameBase + TString("full");
  drawEventScaleFactorGraphs(grScale   , 
			     TString("event scale factor") + TString(buf),
			     plotName, fRoot);

  if (cplotV && (cplotV->size()==4)) {
    const int colorCount=5;
    const int colors[colorCount] = { kBlack, kRed+1, kBlue+1, kGreen+2, kViolet  };
    CPlot* cpFull = (*cplotV)[0];
    CPlot* cpReco = (*cplotV)[1];
    CPlot* cpId = (*cplotV)[2];
    CPlot* cpHlt = (*cplotV)[3];
    int currCollIdx= cpFull->getItemCount() % colorCount;
    int currCol= colors[currCollIdx];
    cpFull->AddGraph(grScale, buf, "PE", currCol);
    cpReco->AddGraph(grScaleReco, buf, "PE", currCol);
    cpId->AddGraph(grScaleId, buf, "PE", currCol);
    cpHlt->AddGraph(grScaleHlt, buf, "PE", currCol);
  }
}

// -------------------------------------------------------------------------

void drawEventScaleFactorsHltFI
   (const TVectorD &scaleHltLeg1FIV, const TVectorD &scaleHltLeg1ErrFIV,
    const TVectorD &scaleHltLeg2FIV, const TVectorD &scaleHltLeg2ErrFIV,
    int rapidityIndex,
    TFile *fRoot,
    std::vector<CPlot*> *cplotV)
{
  int nUnfoldingBins = DYTools::getTotalNumberOfBins();
  if (scaleHltLeg1FIV.GetNoElements()!=nUnfoldingBins) {
    std::cout << "\n\nERROR: drawEventScaleFactorsHltFI should be called for "
	      << "flat-indexed vectors\n\n";
    return;
  }

  
  TGraphErrors *grScaleHltLeg1 = createGraph_vsMassFI(scaleHltLeg1FIV,scaleHltLeg1ErrFIV,rapidityIndex);
  TGraphErrors *grScaleHltLeg2  = createGraph_vsMassFI(scaleHltLeg2FIV,scaleHltLeg2ErrFIV,rapidityIndex);

  char buf[20];
  double rapidity=DYTools::findAbsYValue(0,rapidityIndex);
  sprintf(buf,"_y=%4.2lf_",rapidity); // for file name
  TString yStr=buf;
  yStr.ReplaceAll(".","_"); yStr.ReplaceAll("=","_");
  sprintf(buf,"|y|=%4.2lf",rapidity); // for label
  TString plotName;
  TString plotNameBase = 
    TString("plot_event_scale_") + DYTools::analysisTag + TString("_") + yStr;
  plotName = plotNameBase + TString("hltLeg1");
  drawEventScaleFactorGraphs(grScaleHltLeg1, 
			     TString("HLT_leg1 scale factor ") + TString(buf),
			     plotName, fRoot);
  plotName = plotNameBase + TString("hltLeg2");
  drawEventScaleFactorGraphs(grScaleHltLeg2 , 
			     TString("HLT_leg2 scale factor ") + TString(buf),
			     plotName, fRoot);

  if (cplotV && (cplotV->size()==6)) {
    const int colorCount=5;
    const int colors[colorCount] = { kBlack, kRed+1, kBlue+1, kGreen+2, kViolet  };
    CPlot* cpHltLeg1 = (*cplotV)[4];
    CPlot* cpHltLeg2 = (*cplotV)[5];
    int currCollIdx= cpHltLeg1->getItemCount() % colorCount;
    int currCol= colors[currCollIdx];
    cpHltLeg1->AddGraph(grScaleHltLeg1, buf, "PE", currCol);
    cpHltLeg2->AddGraph(grScaleHltLeg2, buf, "PE", currCol);
  }
}

// -------------------------------------------------------------------------

// This method reads all ROOT files that have efficiencies from
// tag and probe in TMatrixD form and converts the matrices into 
// more simple arrays.
int fillEfficiencyConstants(  const InputFileMgr_t &inpMgr ) {

  int puReweight = inpMgr.puReweightFlag();
  TriggerSelection_t triggers(inpMgr.triggerTag(),true);

  TString fnStart="efficiency_TnP_"; //+ DYTools::analysisTag;
  TString fnEnd=".root";
  if (puReweight) fnEnd= TString("_PU.root");

  if (dataEff.size()) dataEff.clear();
  if (dataEffErrLo.size()) dataEffErrLo.clear();
  if (dataEffErrHi.size()) dataEffErrHi.clear();
  if (dataEffAvgErr.size()) dataEffAvgErr.clear();
  if (mcEff.size()) mcEff.clear();
  if (mcEffErrLo.size()) mcEffErrLo.clear();
  if (mcEffErrHi.size()) mcEffErrHi.clear();
  if (mcEffAvgErr.size()) mcEffAvgErr.clear();
  dataEff.reserve(NEffTypes); 
  dataEffErrLo.reserve(NEffTypes); dataEffErrHi.reserve(NEffTypes);
  dataEffAvgErr.reserve(NEffTypes);
  mcEff.reserve(NEffTypes); 
  mcEffErrLo.reserve(NEffTypes); mcEffErrHi.reserve(NEffTypes);
  mcEffAvgErr.reserve(NEffTypes);

  int res=1;
  for (int kind=0; res && (kind<NEffTypes); ++kind) {    
    DYTools::TEtBinSet_t localEtBinning= etBinning;
    DYTools::TEfficiencyKind_t effKind=DYTools::TEfficiencyKind_t(kind);
      //(etBinning == DYTools::ETBINS6spec) ? DYTools::ETBINS6 : etBinning;
    triggers.actOnData(true);
    TString dataFName=fnStart + 
      getLabel(DYTools::DATA, effKind,
	       inpMgr.tnpEffCalcMethod(DYTools::DATA,effKind),
	       localEtBinning, etaBinning, triggers)
      + fnEnd;
    int weightedCnC= 
      (DYTools::efficiencyIsHLT(effKind))  ?
      		puReweight : 0;
    res=fillOneEfficiency(inpMgr.tnpTag(), dataFName, kind, 
			  dataEff, dataEffErrLo, dataEffErrHi, dataEffAvgErr,
			  weightedCnC);
  }
  for (int kind=0; res && (kind<NEffTypes); ++kind) {
    DYTools::TEtBinSet_t localEtBinning= etBinning;
    DYTools::TEfficiencyKind_t effKind=DYTools::TEfficiencyKind_t(kind);
      // (etBinning == DYTools::ETBINS6spec) ? DYTools::ETBINS6 : etBinning;
    triggers.actOnData(false);
    TString mcFName=fnStart + 
      getLabel(DYTools::MC, effKind,
	       inpMgr.tnpEffCalcMethod(DYTools::MC,effKind),
	       localEtBinning, etaBinning, triggers)
      + fnEnd;
    res=fillOneEfficiency(inpMgr.tnpTag(), mcFName, kind, 
			  mcEff, mcEffErrLo, mcEffErrHi, mcEffAvgErr, puReweight);
  }
  if (res!=1) std::cout << "Error in fillEfficiencyConstants\n"; 
  else std::cout << "fillEfficiencyConstants ok\n";
  return res;
}

// -------------------------------------------------------------------------

int fillOneEfficiency(const TString &dirTag, const TString filename, 
   UInt_t kindIdx, vector<TMatrixD*> &effV, vector<TMatrixD*> &errLoV, 
   vector<TMatrixD*> &errHiV, vector<TMatrixD*> &avgErrV, 
   int weightedCnC) {

  TString fullFName=TString("../root_files/tag_and_probe/")+dirTag+TString("/")+filename;
  TFile *f=new TFile(fullFName);
  if(!f->IsOpen()) {
    if (allowToIgnoreAnalysisTag) {
      std::cout << "... changing analysis tag in efficiency file name=<" 
		<< fullFName << ">\n";
      if (DYTools::analysisTag=="2D") fullFName.ReplaceAll("2D","1D");
      else fullFName.ReplaceAll("1D","2D");
      if(f) delete f;
      f = new TFile(fullFName);
    }
    if (!f || !f->IsOpen()) {
      std::cout << "failed to open a file <" << fullFName << ">" << std::endl;
      assert(0);
    }
  }
  std::cout << "reading <" << filename << ">\nfullFName=<" << fullFName << ">" << std::endl;
  
  TMatrixD *effMatrix        = NULL;
  TMatrixD *effMatrixErrLow  = NULL;
  TMatrixD *effMatrixErrHigh = NULL;
  if (weightedCnC) {
    effMatrix        = (TMatrixD*)f->Get("effArray2DWeighted");
    effMatrixErrLow  = (TMatrixD*)f->Get("effArrayErrLow2DWeighted");
    effMatrixErrHigh = (TMatrixD*)f->Get("effArrayErrHigh2DWeighted");
  }
  else {
    effMatrix        = (TMatrixD*)f->Get("effArray2D");
    effMatrixErrLow  = (TMatrixD*)f->Get("effArrayErrLow2D");
    effMatrixErrHigh = (TMatrixD*)f->Get("effArrayErrHigh2D");
  }
  f->Close();
  delete f;

  // Make sure that the objects are present
  if( !(effMatrix && effMatrixErrLow && effMatrixErrHigh) ) {
    std::cout << "file <" << fullFName << "> does not contain effMatrix,effMatrixErrLow or effMatrixErrHigh" << std::endl;
    assert(0);
  }

  // Make sure that there are only two eta bins and appropriate number of ET bins
  if( effMatrix->GetNcols() != etaBinCount ) {
    std::cout << "The number of eta bins stored in constants file ("
	      << effMatrix->GetNcols() << ") is not "
	      << etaBinCount << "\n";
    return 0;
  }

  /*
  if (etBinning== DYTools::ETBINS6spec) {
    TMatrixD *tmpEff=effMatrix;
    effMatrix=new TMatrixD(etBinCount,etaBinCount);
    TMatrixD *tmpEffErrLow=effMatrixErrLow;
    effMatrixErrLow=new TMatrixD(etBinCount,etaBinCount);
    TMatrixD *tmpEffErrHigh=effMatrixErrHigh;
    effMatrixErrHigh=new TMatrixD(etBinCount,etaBinCount);

    for (int ir=0, idxr=0; ir<effMatrix->GetNrows(); ++ir,++idxr) {
      if (ir==2) idxr--;
      else if (ir==7) idxr--;
      for (int ic=0; ic<effMatrix->GetNcols(); ++ic) {
	(*effMatrix)(ir,ic) = (*tmpEff)(idxr,ic);
	(*effMatrixErrLow)(ir,ic) = (*tmpEffErrLow)(idxr,ic);
	(*effMatrixErrHigh)(ir,ic) = (*tmpEffErrHigh)(idxr,ic);
      }
    }

    if (1) {
      std::cout << "tmpEff= " << tmpEff->GetNrows() << " x " << tmpEff->GetNcols();
      std::cout << "; EtBins= "; 
      for (int i=0; i<DYTools::nEtBins6; i++) std::cout << " " << DYTools::etBinLimits6[i];
      std::cout << "\n";
      std::cout << "effMatrix= " << effMatrix->GetNrows() << " x " << effMatrix->GetNcols();
      std::cout << "; EtBins= "; 
      for (int i=0; i<DYTools::nEtBins6spec; i++) std::cout << " " << DYTools::etBinLimits6spec[i];
      std::cout << "\n";
      std::cout << "\ntmpEff\n"; tmpEff->Print();
      std::cout << "effMatrix\n"; effMatrix->Print();
    }
    delete tmpEff;
    delete tmpEffErrLow;
    delete tmpEffErrHigh;
  }
  */

  if( effMatrix->GetNrows() != etBinCount ) {
    std::cout << "The number of ET bins stored in constants file ("
	      << effMatrix->GetNrows() << ") is different form expected ("
	      << etBinCount << ")\n";
    return 0;
  }

  if ((effV.size()!=kindIdx) ||
      (errLoV.size()!=kindIdx) ||
      (errHiV.size()!=kindIdx) ||
      (avgErrV.size()!=kindIdx)) {
    cout << "Error: effV.size=" << effV.size() 
	 << ", errLoV.size=" << errLoV.size() 
	 << ", errHiV.size=" << errHiV.size() 
	 << ", avgErrV.size=" << avgErrV.size() 
	 << ", kindIdx=" <<  kindIdx << "\n";
    return 0;
  }

  effV.push_back(effMatrix);
  errLoV.push_back(effMatrixErrLow);
  errHiV.push_back(effMatrixErrHigh);
  TMatrixD* avgErr=(TMatrixD*)effMatrixErrLow->Clone("avgErr");
  for (int col=0; col<effMatrixErrLow->GetNcols(); ++col) {
    for (int row=0; row<effMatrixErrLow->GetNrows(); ++row) {
      (*avgErr)[row][col]=
	0.5*((*effMatrixErrLow)[row][col] + (*effMatrixErrHigh)[row][col]);
    }
  }
  avgErrV.push_back(avgErr);

  HERE("leaving fillOneEfficiency");
  return 1;
}


// -------------------------------------------------------------------------

void PrintEffInfoLines(const char *msg, int effKind, int effMethod, 
		       int binCount, const double *eff, const double *effErr) {
  std::cout << "PrintEffInfoLines(" << ((msg) ? msg : "<null>") 
	    << ", effKind=" << effKind << ", effMethod=" << effMethod 
	    << ", binCount=" << binCount << "\n";
  for (int i=0; i<binCount; i++) {
    printf("  i=%d  eff=%6.4lf effErr=%8.6e\n",i,eff[i],effErr[i]);
  }
  std::cout << std::endl;
  return;
}


// -------------------------------------------------------------------------

