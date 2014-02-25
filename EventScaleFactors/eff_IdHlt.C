#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TChain.h>
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TGraphErrors.h>           // graphs
#include <TProfile.h>               // profile histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for Lorentz vector computations
#include <TGraphAsymmErrors.h>
#include <TClonesArray.h>
#include <TMatrixD.h>
#include <TRandom.h>
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <vector>                   // STL vector class
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooBreitWigner.h"
#include "RooExponential.h"
#include "RooCBShape.h"
#include "RooGenericPdf.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFFTConvPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooGaussian.h"

#include "../Include/CPlot.hh"          // helper class for plots
#include "../Include/MitStyleRemix.hh"  // style settings for drawing

// define classes and constants to read in ntuple
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TGenInfo.hh"
#include "../Include/TEventInfo.hh"
#include "../Include/TDielectron.hh"
#include "../Include/TElectron.hh"
#include "../Include/TVertex.hh"
#include "../Include/DYTools.hh"
#include "../Include/DYToolsUI.hh"
#include "../Include/EleIDCuts.hh"
#include "../Include/TriggerSelection.hh"

#include "../EventScaleFactors/cutFunctions.hh"
#include "../EventScaleFactors/fitFunctions.hh"
#include "../EventScaleFactors/fitFunctionsCore.hh"

#include "../EventScaleFactors/tnpSelectEvents.hh"

#include "../Include/EventSelector.hh"

// lumi section selection with JSON files
#include "../Include/JsonParser.hh"

#endif

using namespace mithep;


//=== COMMON CONSTANTS ===========================================================================================

const int evaluate_efficiencies=0;
//const int performPUReweight=0;
const int performOppositeSignTest=1;

//=== FUNCTION DECLARATIONS ======================================================================================

//=== MAIN MACRO =================================================================================================

int eff_IdHlt(const TString configFile, 
	      const TString effTypeString, 
	      int runOnData,
	      DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN,
	      DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST)
{

  gBenchmark->Start("eff_IdHlt");
  

  //  ---------------------------------
  //       Preliminary checks
  //  ---------------------------------

  if (!effTypeString.Contains("ID") &&
      !effTypeString.Contains("HLT")) {
    std::cout << "eff_IdHlt: effTypeString should be either \"ID\" or \"HLT\"\n";
    //return retCodeError; // will abort later, after an additional check
  }

  DYTools::printExecMode(runMode,systMode);

  if (configFile.Contains("_check_")) {
    return retCodeStop;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================

  TDescriptiveInfo_t tnpSection;
  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(configFile,&tnpSection) ||
      !inpMgr.KeepFirstAndLastSample()
      //|| !inpMgr.SetSkimsToNtuples()
      ) return retCodeError;

  // Construct eventSelector, update inpMgr and plot directory
  EventSelector_t evtSelector(inpMgr,runMode,systMode,
			      "","", EventSelector::_selectDefault);

  // Event weight handler
  EventWeight_t evWeight;
  evWeight.init(inpMgr.puReweightFlag(),inpMgr.fewzFlag(),systMode);
  TriggerSelection_t triggers(evtSelector.trigger());

  // Prepare output directory
  TString tagAndProbeDir=inpMgr.tnpDir(systMode,1);
  _TnP_fitPlots_dir = tagAndProbeDir;  // defined in fitFunctionsCore.hh

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  const tnpSelectEvent_t::TCreateBranchesOption_t weightBranch1stStep=
    (inpMgr.puReweightFlag()) ?
            tnpSelectEvent_t::_skipWeight :
	    tnpSelectEvent_t::_dontSkipWeight;
  //TString puStr = (inpMgr.puReweightFlag()) ? "_PU" : "";
 
  Double_t massLow  = 60;
  Double_t massHigh = 120;

  DYTools::TEfficiencyKind_t effType = DetermineEfficiencyKind(effTypeString);
  printf("Efficiency type to measure: %s\n", EfficiencyKindName(effType).Data());
  if ((effType!=DYTools::ID) && !DYTools::efficiencyIsHLT(effType)) {
    std::cout << "eff_IdHlt does not work with <" << EfficiencyKindName(effType) << "> efficiency\n";
    return retCodeError;
  }

  // Read in the configuration file
  DYTools::TDataKind_t dataKind = (runOnData) ? DYTools::DATA : DYTools::MC;
  TString sampleTypeString = (runOnData) ? "DATA" : "MC";
  TString calcMethodString = inpMgr.getTNP_calcMethod(tnpSection,dataKind,effType);
  TString etBinningString  = inpMgr.getTNP_etBinningString(tnpSection);
  TString etaBinningString = inpMgr.getTNP_etaBinningString(tnpSection);
  TString dirTag= inpMgr.selectionTag();

  vector<TString> ntupleFileNames;
  vector<TString> jsonFileNames;
  inpMgr.getTNP_ntuples(tnpSection,runOnData,ntupleFileNames,jsonFileNames);
  if (1) {
    for (unsigned int i=0; i<ntupleFileNames.size(); ++i) {
      std::cout << " i=" << i << ": " << ntupleFileNames[i];
      if (jsonFileNames.size()>i) std::cout << "; json " << jsonFileNames[i];
      std::cout << "\n";
    }
  }

  printf("Efficiency calculation method: %s\n", calcMethodString.Data());
  int calcMethod= DetermineTnPMethod(calcMethodString);

  DYTools::TEtBinSet_t etBinning = DetermineEtBinSet(etBinningString);
  printf("SC ET binning: %s\n", EtBinSetName(etBinning).Data());

  DYTools::TEtaBinSet_t etaBinning = DetermineEtaBinSet(etaBinningString);
  printf("SC eta binning: %s\n", EtaBinSetName(etaBinning).Data());

  printf("Sample: %s\n", sampleTypeString.Data());
  int sample=DetermineDataKind(sampleTypeString);

  // Correct the trigger object
  triggers.actOnData((sample==DYTools::DATA)?true:false);
  evtSelector.setTriggerActsOnData((sample==DYTools::DATA)?true:false);

  if (effType==DYTools::HLT) {
    std::cout //<< "\tHLT efficiency calculation method " 
	      //<< triggers.hltEffCalcName() << 
	      << " triggerSet=" 
	      << triggers.triggerSetName() << "\n";
  }
  else {
    //triggers.hltEffCalcMethod(HLTEffCalc_2011Old);
    std::cout << "WARNING: not setting HLTEffCalc_2011Old\n";
  }

  TRandom *rnd= new TRandom();
  rnd->SetSeed(0); 

 // The label is a string that contains the fields that are passed to
  // the function below, to be used to name files with the output later.
  TString label = getLabel(sample, effType, calcMethod, etBinning, etaBinning, triggers);

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================
  
  
#ifdef tnpSelectEventsIsObject
  tnpSelectEvent_t::Class()->IgnoreTObjectStreamer();
#endif

  //  
  // Set up histograms
  //
  TH1F* hMass           = new TH1F("hMass"     ,"",30,massLow,massHigh);
  TH1F* hMassTotal      = new TH1F("hMassTotal","",30,massLow,massHigh);
  TH1F* hMassPass       = new TH1F("hMassPass" ,"",30,massLow,massHigh);
  TH1F* hMassFail       = new TH1F("hMassFail" ,"",30,massLow,massHigh);

#if (defined DYee7TeV)
  // The source/target histograms are unbiased pile-up distributions
  // for data(target) and MC(source) prepared according to Hildreth's instructions.
  // Note: for 7 TeV, the data histogram for TnP takes into account TnP trigger prescales
  // and the corresponding file is non-standard, so the name is set below.
  // For 8 TeV, these names are blank and the default Hildreth method is used.
  TString puTargetFName = "";
  TString puSourceFName = "";
  if( DYTools::energy8TeV ){
    puTargetFName = "undefined";
  }else{
    if (effType==DYTools::HLT) {
      // HLT efficiency source file
      puTargetFName = "../root_files/pileup/dataPileupHildreth_full2011_TnP_HLT_20121118_repacked.root";
    } else {
      // ID efficiency source file
      puTargetFName = "../root_files/pileup/dataPileupHildreth_full2011_TnP_ID_20121118_repacked.root";
    }
    puSourceFName = "../root_files/pileup/mcPileupHildreth_full2011_20121110_repacked.root";
  }
  
  if (inpMgr.puReweightFlag()) {
    if( DYTools::energy8TeV ){
      // Do not do the check for 8 TeV where the custom files are not used
    }else{
      // Check file existance for 7 TeV
      TFile tmpFile(puTargetFName);
      int npvOk=tmpFile.IsOpen();
      tmpFile.Close();
      if (!npvOk) {
	std::cout << "the file needed of PV-reweighting, <" << puTargetFName << "> does not exist. Run selectEvents.C first\n";
	assert(0);
      }
    }
  }
#endif // end of block for 7TeV analysis


  TFile *templatesFile = 0;
  vector<TH1F*> hPassTemplateV;
  vector<TH1F*> hFailTemplateV;
  if( sample != DYTools::DATA) {
    // For simulation, we will be saving templates
    TString labelMC = 
      getLabel(-1111, effType, calcMethod, etBinning, etaBinning, triggers);
    TString templatesLabel = tagAndProbeDir + 
      TString("/mass_templates_")+labelMC+TString(".root");
    templatesFile = new TFile(templatesLabel,"recreate");
    for(int i=0; i<DYTools::getNEtBins(etBinning); i++){
      for(int j=0; j<DYTools::getNEtaBins(etaBinning); j++){
	TString hname = "hMassTemplate_Et";
	hname += i;
	hname += "_eta";
	hname += j;
	hPassTemplateV.push_back(new TH1F(hname+TString("_pass"),"",60,massLow,massHigh));
	hFailTemplateV.push_back(new TH1F(hname+TString("_fail"),"",60,massLow,massHigh));
      }
    }
  } else {
    // For data, we will be using templates,
    // however, if the request is COUNTnCOUNT, do nothing
    if( calcMethod != DYTools::COUNTnCOUNT ){
      TString labelMC = 
	getLabel(-1111, effType, calcMethod, etBinning, etaBinning, triggers);
      TString templatesLabel = 
	tagAndProbeDir+TString("/mass_templates_")+labelMC+TString(".root");
      templatesFile = new TFile(templatesLabel);
      if( ! templatesFile->IsOpen() )
	assert(0);
    }
  }

  // This file can be utilized in the future, but for now
  // opening it just removes complaints about memory resident
  // trees. No events are actually written.
  TString selectEventsFName=inpMgr.tnpSelectEventsFName(systMode,sampleTypeString,effTypeString,triggers.triggerSetName());
  if (runMode==DYTools::DEBUG_RUN) {
    selectEventsFName.ReplaceAll("/selectEvents","_DEBUG/selectEvents");
  }
  /*
  TString uScore="_";
  TString selectEventsFName=tagAndProbeDir + TString("/selectEvents_") 
    + DYTools::analysisTag + uScore+
    + sampleTypeString + uScore +
    + effTypeString + uScore +  triggers.triggerSetName() + puStr;
  selectEventsFName.Append(".root");
  */
  std::cout << "selectEventsFName=<" << selectEventsFName << ">\n"; 
  TFile *selectedEventsFile = new TFile(selectEventsFName,"recreate");
  if (!selectedEventsFile) {
    assert(0);
  }

  tnpSelectEvent_t storeData;
  const int new_store_data_code=1;
  TTree *passTree = new TTree("passTree","passTree");
  Double_t storeMass, storeEt, storeEta;
  UInt_t storeNGoodPV;
  if (new_store_data_code) {
    storeData.createBranches(passTree, weightBranch1stStep);
  }
  else {
    passTree->Branch("mass",&storeMass,"mass/D");
    passTree->Branch("et",&storeEt  ,"et/D");
    passTree->Branch("eta",&storeEta ,"eta/D");
    passTree->Branch("nGoodPV",&storeNGoodPV,"nGoodPV/i");
  }

  TTree *failTree = new TTree("failTree","failTree");
  if (new_store_data_code) {
    storeData.createBranches(failTree, weightBranch1stStep);
  }
  else {
    failTree->Branch("mass",&storeMass,"mass/D");
    failTree->Branch("et",&storeEt  ,"et/D");
    failTree->Branch("eta",&storeEta ,"eta/D");
    failTree->Branch("nGoodPV",&storeNGoodPV,"nGoodPV/i");
  }

  int nDivisions = DYTools::getNEtBins(etBinning)*DYTools::getNEtaBins(etaBinning);
  double ymax = 800;
  if(nDivisions <4 )
    ymax = nDivisions * 200;
  else if (nDivisions>DYTools::maxTnPCanvasDivisions) {
    nDivisions=DYTools::maxTnPCanvasDivisions;
  }

  TCanvas *c1 = MakeCanvas("c1","c1", 600, int(ymax));
  c1->Divide(2,nDivisions);

  int eventsInNtuple = 0;
  int eventsAfterJson = 0;
  int eventsAfterTrigger = 0;
  int totalCand = 0;
  int totalCandInMassWindow = 0;
  int totalCandInEtaAcceptance = 0;
  int totalCandEtAbove10GeV = 0;
  int totalCandMatchedToGen = 0;
  int totalCandOppositeSign = 0;
  int totalTagProbePairs = 0;

  // Loop over files
  for(UInt_t ifile=0; ifile<ntupleFileNames.size(); ifile++){

    //
    // Access samples and fill histograms
    //  
    TFile *infile = 0;
    TTree *eventTree = 0;
        
    // Data structures to store info from TTrees
    mithep::TEventInfo *info = new mithep::TEventInfo();
    mithep::TGenInfo   *gen  = new mithep::TGenInfo();
    TClonesArray *dielectronArr = new TClonesArray("mithep::TDielectron");
    TClonesArray *pvArr   = new TClonesArray("mithep::TVertex");
     
    // Read input file
    cout << "Processing " << ntupleFileNames[ifile] << "..." << endl;
    infile = new TFile(ntupleFileNames[ifile]); 
    assert(infile);
    
    // Set up JSON for data
    Bool_t hasJSON = kFALSE;
    // mithep::RunLumiRangeMap rlrm;
    JsonParser jsonParser;
    if((jsonFileNames.size()>0) && (jsonFileNames[ifile].CompareTo("NONE")!=0)) { 
      hasJSON = kTRUE;
      // rlrm.AddJSONFile(samp->jsonv[ifile].Data());
      std::cout << "JSON file " << jsonFileNames[ifile] << "\n";
      jsonParser.Initialize(jsonFileNames[ifile].Data()); 
    }

    // Get the TTrees
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",&info);                
    TBranch *infoBr       = eventTree->GetBranch("Info");
    assert(infoBr);

    // check whether the file is suitable for the requested run range
    UInt_t runNumMin = UInt_t(eventTree->GetMinimum("runNum"));
    UInt_t runNumMax = UInt_t(eventTree->GetMaximum("runNum"));
    std::cout << "runNumMin=" << runNumMin << ", runNumMax=" << runNumMax << "\n";
    if (!triggers.validRunRange(runNumMin,runNumMax)) {
      std::cout << "... file contains uninteresting run range\n";
      continue;
    }

    // Define other branches
    eventTree->SetBranchAddress("Dielectron",&dielectronArr); 
    TBranch *dielectronBr = eventTree->GetBranch("Dielectron");
    eventTree->SetBranchAddress("PV", &pvArr); 
    TBranch *pvBr         = eventTree->GetBranch("PV");
    assert(dielectronBr); assert(pvBr);

    TBranch *genBr = 0;
    if(sample != DYTools::DATA){
      eventTree->SetBranchAddress("Gen",&gen);
      genBr = eventTree->GetBranch("Gen");
    }

    TElectron *ele1= new TElectron();
    TElectron *ele2= new TElectron();

    // loop over events    
    eventsInNtuple += eventTree->GetEntries();
     for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
       if (DYTools::isDebugMode(runMode) && (ientry>100000)) break;
       
       if(sample != DYTools::DATA){
	genBr->GetEntry(ientry);
	// If the Z->ll leptons are not electrons, discard this event.
	// This is needed for signal MC samples such as Madgraph Z->ll
	// where all 3 lepton flavors are possible
	if(abs(gen->lid_1) != 11 || abs(gen->lid_2) != 11)
	  continue;
       }
      // Check that the whole event has fired the appropriate trigger
      infoBr->GetEntry(ientry);
      
      // Apply JSON file selection to data
      if(hasJSON && !jsonParser.HasRunLumi(info->runNum, info->lumiSec)) continue;  // not certified run? Skip to next event...
      eventsAfterJson++;

      // Event level trigger cut
      bool idEffTrigger = (effType==DYTools::ID) ? true:false;
      ULong_t eventTriggerBit= triggers.getEventTriggerBit_TagProbe(info->runNum, idEffTrigger);

      if(!(info->triggerBits & eventTriggerBit)) continue;  // no trigger accept? Skip to next event... 
      eventsAfterTrigger++;

      ULong_t tagTriggerObjectBit= triggers.getTagTriggerObjBit(info->runNum,idEffTrigger);
      ULong_t probeTriggerObjectBit_Tight= triggers.getProbeTriggerObjBit_Tight(info->runNum,idEffTrigger);
      ULong_t probeTriggerObjectBit_Loose= triggers.getProbeTriggerObjBit_Loose(info->runNum,idEffTrigger);
//       ULong_t probeTriggerObjectBit= probeTriggerObjectBit_Tight | probeTriggerObjectBit_Loose;

      // loop through dielectrons
      dielectronArr->Clear();
      dielectronBr->GetEntry(ientry);

      mithep::TDielectron *dielectron=NULL;
#ifdef DYee8TeV_reg
      mithep::TDielectron unregDielectron;
#endif

      for(Int_t i=0; i<dielectronArr->GetEntriesFast(); i++) {
	totalCand++;
	dielectron = (mithep::TDielectron*)((*dielectronArr)[i]);

#ifdef DYee8TeV_reg
	if ((systMode==DYTools::UNREGRESSED_ENERGY) ||
	    (systMode==DYTools::UNREG_PU5plus) ||
	    (systMode==DYTools::UNREG_PU5minus) ||
	    (systMode==DYTools::UNREG_TagID) ||
	    (systMode==DYTools::UNREG_TagPt))
	  {
	  unregDielectron.assign(dielectron);
	  //std::cout << "dielectron info    : " << dielectron->mass << ", " << dielectron->pt << ", " << dielectron->y << ", " << dielectron->phi << "\n";
	  unregDielectron.replace2UncorrEn(0); // 0 - do replace, 1 - don't replace
	  dielectron= &unregDielectron;
	  //std::cout << "dielectron info (2): " << dielectron->mass << ", " << dielectron->pt << ", " << dielectron->y << ", " << dielectron->phi << "\n";
	}
#endif

	// Tag and probe is done around the Z peak
	if((dielectron->mass < massLow) || (dielectron->mass > massHigh)) continue;
	totalCandInMassWindow++;

	// Check whether each electron is in the rapidity gap, but not cut yet
	bool isGap1 = DYTools::isEcalGap( dielectron->scEta_1);
	bool isGap2 = DYTools::isEcalGap( dielectron->scEta_2);
	// The bollean for the cut on the probe. We remove the rapidity
	// gap for probes only for a special eta binning
	bool passGapCut1 = true;
	bool passGapCut2 = true;
	if(etaBinning == DYTools::ETABINS2){
	  if(isGap1)
	    passGapCut1 = false;
	  if(isGap2)
	    passGapCut2 = false;
	}

	// ECAL acceptance cut on supercluster Et
	// Note: the main analysis may or may not allow eta go to kECAL_MAX_ETA.
	// Here, we simply require max possible for ecal.
	// For the tag, this is always ok. For the probe, the tag and probe
	// procedure is done in eta bins anyways.
	if((fabs(dielectron->scEta_1) > DYTools::kECAL_MAX_ETA)       
	   || (fabs(dielectron->scEta_2) > DYTools::kECAL_MAX_ETA)) continue;  // outside eta range? Skip to next event...
	totalCandInEtaAcceptance++;
	// None of the electrons should be below 10 GeV
	if((dielectron->pt_1 < 10)               || (dielectron->pt_2 < 10))	      continue;  // below supercluster ET cut? Skip to next event...
	totalCandEtAbove10GeV++;
	
	// Next, we will do a loose kinematic matching to generator level
	// info. 
	// For the data, this is not needed and not done. We take all
	// candidates, and take care of background by fitting.
	// For MC, however, we do not fit, but count pass/fail events.
	// So we need to make sure there is no background. However, even
	// in the signal Z->ee MC sample there jets and therefore fake
	// electrons. So we drop all candidates that do not have both leptons
	// matched.
	// 
	if( sample != DYTools::DATA )
	  if( ! dielectronMatchedToGeneratorLevel(gen, dielectron) ) continue;
	totalCandMatchedToGen++;

	if (performOppositeSignTest && ( dielectron->q_1 == dielectron->q_2 )) continue;
	totalCandOppositeSign++;

	// ECAL driven: this condition is NOT applied	

	// Preliminary selection is complete. Now work on tags and probes.
	
	dielectron->extractElectron(1,*ele1);
	dielectron->extractElectron(2,*ele2);
	bool isTag1 = false;
	bool isTag2 = false;
	
	// Syst mode for tag considers two cases:
	//   Resolution_study or TAG_PT: lower elePt
	//   FSR_study or TAG_ID: mediumID
	// other systematics do not affect the selection, thus
	// systMode branching for isTag() might be removed.
	if (systMode == DYTools::NO_SYST) {
	  isTag1=isTag(ele1, tagTriggerObjectBit, info->rhoLowEta);
	  isTag2=isTag(ele2, tagTriggerObjectBit, info->rhoLowEta);
	}
	else {
	  isTag1=isTag_systStudy(ele1, tagTriggerObjectBit, info->rhoLowEta, systMode);
	  isTag2=isTag_systStudy(ele2, tagTriggerObjectBit, info->rhoLowEta, systMode);
	}
	
	// Any electron that made it here is eligible to be a probe
	// for ID cuts.
	bool isIDProbe1     = true && passGapCut1;
	bool isIDProbe2     = true && passGapCut2;
	bool isIDProbePass1 = passID(ele1, info->rhoLowEta) && passGapCut1;
	bool isIDProbePass2 = passID(ele2, info->rhoLowEta) && passGapCut2;
	
	// Probes for HLT cuts:
	// For a univeral pass/fail of either leading or trailing electron, 
	// take into account the trigger thresholds with respect to the Et of the
	// probe.
	//    The "Tight", or the leading trigger leg, is presently Ele17,
	// so match to it only electrons above 20 GeV to avoid turn-on.
	//    WARNING: this is a bit dangerous because some day we will
	// switch to higher-Et triggers. To be improved.
	ULong_t probeTriggerObjectBit_probe1 = probeTriggerObjectBit_Loose;
	if( ele1->scEt > 20 ) 
	  probeTriggerObjectBit_probe1 |= probeTriggerObjectBit_Tight;
	ULong_t probeTriggerObjectBit_probe2 = probeTriggerObjectBit_Loose;
	if( ele2->scEt > 20 ) 
	  probeTriggerObjectBit_probe2 |= probeTriggerObjectBit_Tight;
	
	bool isHLTProbe1     = passID(ele1, info->rhoLowEta) && passGapCut1;
	bool isHLTProbe2     = passID(ele2, info->rhoLowEta) && passGapCut2;
	bool isHLTProbePass1 = ( isHLTProbe1 && (ele1->hltMatchBits & probeTriggerObjectBit_probe1) && passGapCut1 ) ;
	bool isHLTProbePass2 = ( isHLTProbe2 && (ele2->hltMatchBits & probeTriggerObjectBit_probe2) && passGapCut2);
	bool isHLTProbePass1tight = ( isHLTProbe1 && (ele1 ->hltMatchBits & probeTriggerObjectBit_Tight) && passGapCut1);
	bool isHLTProbePass2tight = ( isHLTProbe2 && (ele2 ->hltMatchBits & probeTriggerObjectBit_Tight) && passGapCut2);
	bool isHLTProbePass1loose = ( isHLTProbe1 && (ele1 ->hltMatchBits & probeTriggerObjectBit_Loose) && passGapCut1);
	bool isHLTProbePass2loose = ( isHLTProbe2 && (ele2 ->hltMatchBits & probeTriggerObjectBit_Loose) && passGapCut2);

	// 
	//  Apply tag and probe, and accumulate counters or histograms
	//       
	
	bool isProbe1     = false;
	bool isProbe2     = false;
	bool isProbePass1 = false;
	bool isProbePass2 = false;
	switch( effType ) {
	case DYTools::ID:
	  isProbe1     = isIDProbe1;
	  isProbe2     = isIDProbe2;
	  isProbePass1 = isIDProbePass1;
	  isProbePass2 = isIDProbePass2;
	  break;
	case DYTools::HLT:
	 //case DYTools::HLT_rndTag:
	  isProbe1     = isHLTProbe1;
	  isProbe2     = isHLTProbe2;
	  isProbePass1 = isHLTProbePass1;
	  isProbePass2 = isHLTProbePass2;
	  //if ((effType==DYTools::HLT_rndTag) 
	  //    && triggers.useRandomTagTnPMethod(info->runNum)) {
	  //  std::cout << "random tag\n";
	  //  if (rnd->Uniform() <= 0.5) {
	  //    // tag is 1st electron
	  //    if (!isTag1) continue;
	  //    isTag2=0; // ignore whether ele2 can be a tag
	  //  }
	  //  else {
	  //    if (!isTag2) continue;
	  //    isTag1=0; // ignore whether ele1 can be a tag
	  //  }
	  //}
	  break;
	case DYTools::HLT_leg1:
	  isProbe1 = isHLTProbe1;
	  isProbe2 = isHLTProbe2;
	  isProbePass1 = isHLTProbePass1tight;
	  isProbePass2 = isHLTProbePass2tight;
	  break;
	case DYTools::HLT_leg2:
	  isProbe1 = isHLTProbe1;
	  isProbe2 = isHLTProbe2;
	  isProbePass1 = isHLTProbePass1loose;
	  isProbePass2 = isHLTProbePass2loose;
	  break;
	default:
	  printf("ERROR: unknown efficiency type requested\n");
	}

	// get the number of goodPVs
	pvBr->GetEntry(ientry);
	storeNGoodPV=0;
	// For data, we use the count of good reconstructed vertices
	// but for MC, since this will be used for PU reweighting,
	// we are using the gen-level number of simulated PU.
	if( (sample==DYTools::DATA) ) {
	  storeNGoodPV = countGoodVertices(pvArr);
	}else{
#ifdef DYee8TeV_reg
	  storeNGoodPV = int(info->nPUmean);
#else
	  storeNGoodPV = info->nPU;
#endif
	}
	
	storeMass = dielectron->mass;
	double event_weight=1.0;

	if(isTag1)
	  totalTagProbePairs++;
	if(isTag2)
	  totalTagProbePairs++;

	// First electron is the tag, second is the probe
	if( isTag1 && isProbe2){
	  // total probes
	  hMassTotal->Fill(dielectron->mass);
	  storeEt   = dielectron->scEt_2;
	  storeEta  = dielectron->scEta_2;
	  if (new_store_data_code) {
#ifdef tnpStoreTag
	    storeData.assignTag(dielectron->scEt_1,dielectron->scEta_1);
#endif
	    storeData.assign(dielectron->mass,dielectron->y,
			     dielectron->scEt_2,dielectron->scEta_2,
			     storeNGoodPV,event_weight,1.);
	  }
	  int templateBin = 
	    getTemplateBin( DYTools::findEtBin(storeEt,etBinning),
			    DYTools::findEtaBin(storeEta,etaBinning),
			    etaBinning);
	  if( isProbePass2 ){
	    // passed
	    hMassPass->Fill(dielectron->mass);
	    passTree->Fill();
	    if(sample != DYTools::DATA && templateBin != -1)
	      hPassTemplateV[templateBin]->Fill(dielectron->mass);
	  }else{
	    // fail
	    hMassFail->Fill(dielectron->mass);
	    failTree->Fill();
	    if(sample != DYTools::DATA && templateBin != -1)
	      hFailTemplateV[templateBin]->Fill(dielectron->mass);
	  }
	}
	// Second electron is the tag, first is the probe
	if( isTag2 && isProbe1 ){
	  // total probes
	  hMassTotal->Fill(dielectron->mass);
	  storeEt   = dielectron->scEt_1;
	  storeEta  = dielectron->scEta_1;
	  if (new_store_data_code) {
#ifdef tnpStoreTag
	    storeData.assignTag(dielectron->scEt_2,dielectron->scEta_2);
#endif
	    storeData.assign(dielectron->mass,dielectron->y,
			     dielectron->scEt_1,dielectron->scEta_1,
			     storeNGoodPV,event_weight,1.);
	  }
	  int templateBin = 
	    getTemplateBin( DYTools::findEtBin(storeEt,etBinning),
			    DYTools::findEtaBin(storeEta,etaBinning),
			    etaBinning);
	  if( isProbePass1 ){
	    // passed
	    hMassPass->Fill(dielectron->mass);
	    passTree->Fill();
	    if(sample != DYTools::DATA && templateBin != -1)
	      hPassTemplateV[templateBin]->Fill(dielectron->mass);
	  }else{
	    // fail
	    hMassFail->Fill(dielectron->mass);
	    failTree->Fill();
	    if(sample != DYTools::DATA && templateBin != -1)
	      hFailTemplateV[templateBin]->Fill(dielectron->mass);
	  }
	}
	
	// In case the full selection is applied:
	//       if( !(isTag1 && ele2_passID) && !(isTag2 && ele1_passID) ) continue;
	if( !(isTag1 && isIDProbePass2) && !(isTag2 && isIDProbePass1) ) continue;
	//       if( !(isTag1) && !(isTag2) ) continue;
	// Fill histogram
	hMass->Fill(dielectron->mass);
	
      } // end loop over dielectron candidates
    } // end loop over events
  
    delete infile;
    infile=0;
    eventTree=0;
    
    delete gen;
    delete info;
    delete dielectronArr;
  } // end loop over files

  // save the selected trees
  selectedEventsFile->cd();
  passTree->Write();
  failTree->Write();
  selectedEventsFile->Write();

  if (inpMgr.puReweightFlag()) {
    selectedEventsFile->Close();
    delete selectedEventsFile;

    TString outFNamePV = tagAndProbeDir + 
      TString("/npv_tnp") + effTypeString + TString("_") + 
      sampleTypeString + DYTools::analysisTag + TString(".root");
    //
    // The names below are actually used for 7 TeV data/MC.
    // For 8 TeV, they are not used.
    TString puTargetDistrName="pileup_lumibased_data";
    TString puSourceDistrName="pileup_simulevel_mc";
#if !(defined DYee7TeV)
    TString puTargetFName="", puSourceFName=""; // not used
#endif
//     TString refDistribution="hNGoodPV_data";

    TString sampleNameBase= effTypeString + TString("_") + 
      sampleTypeString + DYTools::analysisTag;
    bool isMC = (sample==DYTools::MC);
    int res=CreatePUWeightedBranch(systMode,
				   selectEventsFName,
				   puTargetFName, puTargetDistrName,
				   puSourceFName, puSourceDistrName,
// 				   outFNamePV, sampleNameBase, 
				   isMC);
    assert(res);
    selectedEventsFile=new TFile(selectEventsFName);
    assert(selectedEventsFile);
    passTree= (TTree*)selectedEventsFile->Get("passTree");
    failTree= (TTree*)selectedEventsFile->Get("failTree");
    assert(passTree); assert(failTree);
  }

  //
  // Efficiency analysis
  //
  
//   printf("Number of regular candidates:      %15.0f\n", hMass->GetSumOfWeights());
  printf("Total events in ntuple                                       %15d\n",eventsInNtuple);
  printf("    events after JSON selection (data)                       %15d\n",eventsAfterJson);
  printf("    events after event level trigger cut                     %15d\n",eventsAfterTrigger);
  printf("\nTotal candidates (no cuts)                                   %15d\n",totalCand);
  printf("        candidates in 60-120 mass window                     %15d\n",totalCandInMassWindow);
  printf("        candidates in pseudorapidity acceptance              %15d\n",totalCandInEtaAcceptance);
  printf("        candidates, both electrons above 10 GeV              %15d\n",totalCandEtAbove10GeV);
  printf("        candidates matched to GEN level (if MC)              %15d\n",totalCandMatchedToGen);
  if (performOppositeSignTest) 
    printf("        candidates opposite sign                             %15d\n",totalCandOppositeSign);
  printf("Number of tag-probe pairs                                    %15d\n", totalTagProbePairs);
  printf("\nNumber of probes, total                                      %15.0f\n", hMassTotal->GetSumOfWeights());
  printf("Number of probes, passed                                     %15.0f\n", hMassPass->GetSumOfWeights());
  printf("Number of probes, failed                                     %15.0f\n", hMassFail->GetSumOfWeights());

  if (evaluate_efficiencies) {
  // Human-readbale text file to store measured efficiencies
  TString reslog = tagAndProbeDir+
    TString("/efficiency_TnP_")+label+TString(".txt");
  ofstream effOutput;
  effOutput.open(reslog);
  // Print into the results file the header.
  effOutput << "Efficiency calculation method: " << calcMethodString.Data() << endl;
  effOutput << "Efficiency type to measure: " << effTypeString.Data() << endl;
  effOutput << "SC ET binning: " << etBinningString.Data() << endl;
  effOutput << "SC eta binning: " << etaBinningString.Data() << endl;
  effOutput << "Sample: " << sampleTypeString.Data() << endl;
  effOutput << "Files processed: " << endl;
  for(UInt_t i=0; i<ntupleFileNames.size(); i++)
    effOutput << "   " << ntupleFileNames[i].Data() << endl;
  
  // ROOT file to store measured efficiencies in ROOT format
  TString resrootBase = tagAndProbeDir+
    TString("/efficiency_TnP_")+label;
  //TString resroot = tagAndProbeDir+
  //  TString("/efficiency_TnP_")+label+TString(".root");
  //TFile *resultsRootFile = new TFile(resroot,"recreate");

  // Fit log 
  TString fitlogname = tagAndProbeDir+
    TString("/efficiency_TnP_")+label+TString("_fitlog.dat");
  ofstream fitLog;
  fitLog.open(fitlogname);

  //
  //  Find efficiency
  //
  bool useTemplates = false;
  if(sample == DYTools::DATA && effType == DYTools::ID &&
     (calcMethod == DYTools::COUNTnFIT || DYTools::FITnFIT) )
    useTemplates = true;

  int NsetBins=120;
//   bool isRECO=0;
  const char* setBinsType="fft";

  measureEfficiencyPU(passTree, failTree,
		    calcMethod, etBinning, etaBinning, c1, effOutput, fitLog,
		      useTemplates, templatesFile, 
		      resrootBase,
		      //resultsRootFile,
		      NsetBins, effType, setBinsType, 
		      dirTag, triggers.triggerSetName(),0);

  effOutput.close();
  fitLog.close();
  TString command = "cat ";
  command += reslog;
  system(command.Data());

  TString fitpicname = tagAndProbeDir+
    TString("/efficiency_TnP_")+label;
  if (calcMethod==DYTools::COUNTnCOUNT) fitpicname.Append(".png"); else fitpicname.Append("_fit.png");
  c1->SaveAs(fitpicname);

  // Save MC templates
  if(sample != DYTools::DATA){
    templatesFile->cd();
    for(int i=0; i<DYTools::getNEtBins(etBinning); i++){
      for(int j=0; j<DYTools::getNEtaBins(etaBinning); j++){
	int templateBin = getTemplateBin( i, j, etaBinning);
	hPassTemplateV[templateBin]->Write();
	hFailTemplateV[templateBin]->Write();
      }
    }
    templatesFile->Close();
  }
  }

  selectedEventsFile->Close();
  std::cout << "selectedEventsFile <" << selectEventsFName << "> saved\n";
  gBenchmark->Show("eff_IdHlt");
  
  return retCodeOk;
}


