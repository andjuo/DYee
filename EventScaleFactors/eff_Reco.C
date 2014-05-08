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
#include "RooGaussian.h"
#include "RooGenericPdf.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooFFTConvPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooDataSet.h"
#include "RooFitResult.h"

#include "../Include/CPlot.hh"          // helper class for plots
#include "../Include/MitStyleRemix.hh"  // style settings for drawing

// define classes and constants to read in ntuple
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TGenInfo.hh"
#include "../Include/TEventInfo.hh"
#include "../Include/TDielectron.hh"
#include "../Include/TElectron.hh"
#include "../Include/TPhoton.hh"
#include "../Include/TVertex.hh"
#include "../Include/DYTools.hh"
#include "../Include/EleIDCuts.hh"
#include "../Include/TriggerSelection.hh"

#include "../EventScaleFactors/cutFunctions.hh"
#include "../EventScaleFactors/fitFunctions.hh"
#include "../EventScaleFactors/fitFunctionsCore.hh"

#include "../Include/EventSelector.hh"
#include "../EventScaleFactors/tnpSelectEvents.hh"


#endif


using namespace mithep;


//=== COMMON CONSTANTS ===========================================================================================

const int evaluate_efficiencies=0;
//const int performPUReweight=0;

//=== FUNCTION DECLARATIONS ======================================================================================

//=== MAIN MACRO =================================================================================================

int eff_Reco(int analysisIs2D,
	     const TString configFile, 
	     const TString effTypeString, 
	     int runOnData,
	     DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN,
	     DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST)
{
  using namespace mithep; 
  gBenchmark->Start("eff_Reco");

  //  ---------------------------------
  //       Preliminary checks
  //  ---------------------------------

  if (!effTypeString.Contains("RECO")) {
    std::cout << "eff_Reco: effTypeString should be \"RECO\"\n";
    return retCodeError;
  }

  DYTools::printExecMode(runMode,systMode);

  if (configFile.Contains("_check_")) {
    return retCodeStop;
  }

  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return retCodeError;
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
			      "", "", EventSelector::_selectDefault);
  TriggerSelection_t triggers(evtSelector.trigger());

  // Event weight handler
  EventWeight_t evWeight;
  evWeight.init(inpMgr.puReweightFlag(),inpMgr.fewzFlag(),systMode);

  // Prepare output directory
  TString tagAndProbeDir=inpMgr.tnpDir(systMode,1);
  _TnP_fitPlots_dir = tagAndProbeDir;  // defined in fitFunctionsCore.hh

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  const tnpSelectEvent_t::TCreateBranchesOption_t weightBranch1stStep=
    (inpMgr.puReweightFlag()) ?
      tnpSelectEvent_t::_skipWeight :
      tnpSelectEvent_t::_dontSkipWeight;
  //TString puStr=(inpMgr.puReweightFlag()) ? "_PU" : "";
  
  Double_t massLow  = 60;
  Double_t massHigh = 120;

  // Read in the configuration file
  DYTools::TDataKind_t dataKind = (runOnData) ? DYTools::DATA : DYTools::MC;
  TString sampleTypeString = (runOnData) ? "DATA" : "MC";
  TString calcMethodString = inpMgr.getTNP_calcMethod(tnpSection,dataKind,DYTools::RECO);
  std::cout << "calcMethodString=" << calcMethodString << "\n";

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

  DYTools::TEfficiencyKind_t effType = DetermineEfficiencyKind(effTypeString);
  printf("Efficiency type to measure: %s\n", EfficiencyKindName(effType).Data());
  if ( effType != DYTools::RECO ) {
    std::cout << "effReco works with RECO efficiency only\n";
    assert(0);
  }

  DYTools::TEtBinSet_t etBinning = DetermineEtBinSet(etBinningString);
  printf("SC ET binning: %s\n", EtBinSetName(etBinning).Data());

  DYTools::TEtaBinSet_t etaBinning = DetermineEtaBinSet(etaBinningString);
  printf("SC eta binning: %s\n", EtaBinSetName(etaBinning).Data());

  printf("Sample: %s\n", sampleTypeString.Data());
  int sample=DetermineDataKind(sampleTypeString);

  // Correct the trigger object
  triggers.actOnData((sample==DYTools::DATA)?true:false);
  evtSelector.setTriggerActsOnData((sample==DYTools::DATA)?true:false);

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
  //   TH1F* hMass           = new TH1F("hMass"     ,"",30,massLow,massHigh);
  TH1F* hMassTotal      = new TH1F("hMassTotal","",30,massLow,massHigh);
  TH1F* hMassPass       = new TH1F("hMassPass" ,"",30,massLow,massHigh);
  TH1F* hMassFail       = new TH1F("hMassFail" ,"",30,massLow,massHigh);

#if (defined DYee7TeV)
  // Here, for 7 TeV we set the names of the source/target histograms with 
  // unbiased pile-up distributions that are  prepared according to Hildreth's instructions.
  // The data histogram for TnP takes into account TnP trigger prescales.
  // For 8 TeV, this is not needed and the default ones will be used
  TString puTargetFName = "";
  TString puSourceFName = "";
  if( DYTools::energy8TeV ){
    puTargetFName = "undefined";
    puSourceFName = "undefined";
  }else{
    puTargetFName = "../root_files/pileup/dataPileupHildreth_full2011_TnP_RECO_20121118_repacked.root";
    puSourceFName = "../root_files/pileup/mcPileupHildreth_full2011_20121110_repacked.root";
  }

  if (inpMgr.puReweightFlag()) {
    if( DYTools::energy8TeV ){
      // Do not do the check for 8 TeV where the custom files are not used
    }else{
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
    TString templatesLabel = 
      tagAndProbeDir + TString("/mass_templates_")+ labelMC+TString(".root");
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

    if (evaluate_efficiencies) {
      // For data, we will be using templates
      // however, if the request is COUNTnCOUNT, do nothing
      if( calcMethod != DYTools::COUNTnCOUNT ){
	TString labelMC =
	  getLabel(-1111,effType, calcMethod, etBinning, etaBinning, triggers);
	TString templatesLabel =
	  tagAndProbeDir+TString("/mass_templates_")+labelMC+TString(".root");
	templatesFile = new TFile(templatesLabel);
	if( ! templatesFile->IsOpen() ) {
	  std::cout << "templatesFile name " << templatesLabel << "\n";
	  assert(0);
	}
      }
    }
  }

  // This file is utilized by fit_EffReco
  TString selectEventsFName=inpMgr.tnpSelectEventsFName(systMode,sampleTypeString,effTypeString,triggers.triggerSetName());
    /*
  TString uScore="_";
  TString selectEventsFName=tagAndProbeDir + TString("/selectEvents_") 
    + DYTools::analysisTag + uScore
    + sampleTypeString + uScore +
    + effTypeString + uScore +  triggers.triggerSetName();
  if (inpMgr.puReweightFlag()) selectEventsFName.Append(puStr);
  selectEventsFName.Append(".root");
    */
  std::cout << "selectEventsFName=<" << selectEventsFName << ">\n"; 
  TFile *selectedEventsFile = new TFile(selectEventsFName,"recreate");
  if(!selectedEventsFile) 
    assert(0);

  tnpSelectEvent_t storeData;
  const int new_store_data_code=1;
  TTree *passTree = new TTree("passTree","passTree");
  assert(passTree);
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
  assert(failTree);
  if (new_store_data_code) {
    storeData.createBranches(failTree, weightBranch1stStep);
  }
  else {
    failTree->Branch("mass",&storeMass,"mass/D");
    failTree->Branch("et",&storeEt  ,"et/D");
    failTree->Branch("eta",&storeEta ,"eta/D");
    failTree->Branch("nGoodPV",&storeNGoodPV,"nGoodPV/i");
  }


  int eventsInNtuple = 0;
  int eventsAfterTrigger = 0;
  int eventsAfterJson = 0;
  int eventsAfterMET  = 0;
  int tagCand = 0;
  int tagCandPassEt = 0;
  int tagCandPassEta = 0;
  int tagCandGenMatched = 0;
  int tagCandEcalDriven = 0;
  int tagCandFinalCount = 0;
  int numTagProbePairs = 0;
  int numTagProbePairsPassEt = 0;
  int numTagProbePairsPassEta = 0;
  int numTagProbePairsGenMatched = 0;
  int numTagProbePairsInMassWindow = 0;
  int numTagProbePairsPassSCIso = 0;
  
  // Loop over files
  for(UInt_t ifile=0; ifile<ntupleFileNames.size(); ifile++){

    //
    // Access samples and fill histograms
    //  
    TFile *infile = 0;
    TTree *eventTree = 0;
        
    // Data structures to store info from TTrees
    mithep::TEventInfo *info = new mithep::TEventInfo();
    mithep::TGenInfo      *gen  = new mithep::TGenInfo();
    TClonesArray *scArr   = new TClonesArray("mithep::TPhoton");
    TClonesArray *eleArr  = new TClonesArray("mithep::TElectron");
    TClonesArray *pvArr   = new TClonesArray("mithep::TVertex");
    
    // Read input file
    cout << "Processing " << ntupleFileNames[ifile] << "..." << endl;
    if (ntupleFileNames[ifile].Index("tight-loose_skim")!=-1) {
      std::cout << "eff_Reco cannot work with 'tight-loose_skim' file: TElectron and TPhoton branches are needed" << endl;
      assert(0);
    }
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

    // check whether the file is suitable for the requested run range
    UInt_t runNumMin = UInt_t(eventTree->GetMinimum("runNum"));
    UInt_t runNumMax = UInt_t(eventTree->GetMaximum("runNum"));
    std::cout << "runNumMin=" << runNumMin << ", runNumMax=" << runNumMax << "\n";
    if (!triggers.validRunRange(runNumMin,runNumMax)) {
      std::cout << "... file contains uninteresting run range\n";
      continue;
    }
    
    // Define other branches
    eventTree->SetBranchAddress("Photon"  ,&scArr); 
    eventTree->SetBranchAddress("Electron",&eleArr); 
    eventTree->SetBranchAddress("PV", &pvArr); 
    TBranch *electronBr   = eventTree->GetBranch("Electron");
    TBranch *photonBr     = eventTree->GetBranch("Photon");
    TBranch *pvBr         = eventTree->GetBranch("PV");
    assert(electronBr); assert(photonBr);
    assert(pvBr);
    TBranch *genBr = 0;
    if(sample != DYTools::DATA){
      eventTree->SetBranchAddress("Gen",&gen);
      genBr = eventTree->GetBranch("Gen");
    }

    // loop over events    
    eventsInNtuple += eventTree->GetEntries();
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      //for(UInt_t ientry=0; ientry<1000; ientry++) { 
      if (DYTools::isDebugMode(runMode) && (ientry>100000)) break;  // This is for faster turn-around in testing
      
      // Check that the whole event has fired the appropriate trigger
      infoBr->GetEntry(ientry);

      // Apply JSON file selection to data
      if(hasJSON && !jsonParser.HasRunLumi(info->runNum, info->lumiSec)) continue;  // not certified run? Skip to next event...
      eventsAfterJson++;

      // Event level trigger cut
      ULong_t eventTriggerBit= triggers.getEventTriggerBit_SCtoGSF(info->runNum);
      ULong_t tagTriggerObjectBit= triggers.getLeadingTriggerObjectBit_SCtoGSF(info->runNum);

      if(!(info->triggerBits & eventTriggerBit)) continue;  // no trigger accept? Skip to next event... 
      eventsAfterTrigger++;
      
      // Apply event-level pfMET
      if( !(info->pfMET < 20 ) ) continue;
      eventsAfterMET++;

      // get analysis object arrays
      if(sample != DYTools::DATA)  
	genBr->GetEntry(ientry);

      eleArr->Clear();
      electronBr->GetEntry(ientry);
      scArr->Clear();
      photonBr->GetEntry(ientry);


      mithep::TElectron *electron=NULL;
#ifdef DYee8TeV_reg
      mithep::TElectron unregElectron;
#endif

      // Loop over the tag electrons
      for(int iele = 0; iele < eleArr->GetEntriesFast(); iele++){
	
	tagCand++;
	electron = (mithep::TElectron*)((*eleArr)[iele]);

#ifdef DYee8TeV_reg
	if ((systMode==DYTools::UNREGRESSED_ENERGY) ||
	    (systMode==DYTools::UNREG_PU5plus) ||
	    (systMode==DYTools::UNREG_PU5minus) ||
	    (systMode==DYTools::UNREG_TagID) ||
	    (systMode==DYTools::UNREG_TagPt))
	  {
	  unregElectron.assign(electron);
	  unregElectron.replace2UncorrEn();
	  electron= &unregElectron;
	}
#endif

	// All cuts for the tag electron should be applied here
	if(electron->scEt<20) continue;
	tagCandPassEt++;

	// For the tag, always exclude rapidity gap
	bool isInGap = DYTools::isEcalGap(electron->scEta);
	if ( isInGap ) continue;

	if( fabs(electron->scEta) > DYTools::kECAL_MAX_ETA) continue;
	tagCandPassEta++;
	
	if( sample != DYTools::DATA)
	  if( ! electronMatchedToGeneratorLevel(gen, electron) ) continue;
	
	tagCandGenMatched++;

	// ECAL driven: this condition is NOT applied	
	
	if( !isTag( electron, tagTriggerObjectBit, info->rhoLowEta) ) continue;
	// Syst mode for tag considers two cases:
	//   Resolution_study or TAG_PT: lower elePt
	//   FSR_study or TAG_ID: mediumID
	// other systematics do not affect the selection, thus
	// systMode branching for isTag() might be removed.
	bool tagOk=false;
	if (systMode == DYTools::NO_SYST) {
	  tagOk=isTag( electron, tagTriggerObjectBit, info->rhoLowEta);
	}
	else {
	  tagOk=isTag_systStudy( electron, tagTriggerObjectBit, info->rhoLowEta, systMode);
	}
	if (!tagOk) continue;

	tagCandFinalCount++;
	
	
	// Loop over superclusters in this event: the probes
	// Note: each supercluster has a number assigned: scID,
	// and each electron object from TElectron collection has
	// the field scID that tells from which supercluster this electron
	// comes from. That allows to make the match between the 
	// object in TPhoton collection and TElectron collection to
	// find superclusters reconstructed as electrons.
	
	for(int isc = 0; isc < scArr->GetEntriesFast(); isc++){
	  
	  const TPhoton *sc = (TPhoton*)((*scArr)[isc]);
	  // Avoid probe that is same as tag
	  if( sc->scID == electron->scID ) continue;

	  numTagProbePairs++;
	  // Apply probe cuts
	  if(sc->scEt < 10) continue;
	  numTagProbePairsPassEt++;

	  // For the probe, exclude eta gap only for one specific eta 
	  // binning, barrel/endcap split
	  if(etaBinning == DYTools::ETABINS2){
	    if( DYTools::isEcalGap( sc->scEta ) ) continue;
	  }

	  if( fabs(sc->scEta) > DYTools::kECAL_MAX_ETA) continue;
	  numTagProbePairsPassEta++;

	  if( sample != DYTools::DATA)
	    if( ! scMatchedToGeneratorLevel(gen, sc) ) continue;
	  numTagProbePairsGenMatched++;

	  // Tracker isolation cut helps to clean up
	  // the probes, but introduces a small bias
	  // ~3% below 20 GeV, and <1% above 20 GeV
 	  if( ientry<1) printf("\n\n WARNING! Cut on probe isolation is applied! Beware of a small bias\n");
 	  if( fabs(sc->trkIso04)/sc->pt > 0.15 ) continue;
	  numTagProbePairsPassSCIso++;
	  
	  // Find mass of the electron-supercluster pair
	  TLorentzVector ele4V, sc4V, dycand4V;
	  ele4V.SetPtEtaPhiM(electron->pt, electron->eta, electron->phi, 0.000511);
	  sc4V .SetPtEtaPhiM(sc->pt, sc->eta, sc->phi, 0.000511);
	  dycand4V = ele4V + sc4V;
	  double mass = dycand4V.M();	  
	  // Tag and probe is done around the Z peak
	  if((mass < massLow) || (mass > massHigh)) continue;
	  numTagProbePairsInMassWindow++;

	  // The probes are fully selected at this point.

	  // Loop over electron collection again to find a match to this supercluster
	  // Match only to ECAL-driven GSF electrons
	  const TElectron *electronMatch = 0;
	  for(int iele2 = 0; iele2 < eleArr->GetEntriesFast(); iele2++){
	    const TElectron *electron2 = (TElectron*)((*eleArr)[iele2]);
	    if( sc->scID == electron2->scID ){
	      // Check ecal driven bits
	      if( electron2->typeBits & kEcalDriven )
		electronMatch = electron2;
	    }
	  } // end loop over electrons searching for SC match
	  
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

	  // total probes
	  double event_weight=1.0;
	  double ee_rapidity=0.;
	  hMassTotal->Fill(mass);
	  storeMass = mass;
	  storeEt   = sc->scEt;
	  storeEta  = sc->scEta;
	  if (new_store_data_code) {
#ifdef tnpStoreTag
	    storeData.assignTag(electron->scEt,electron->scEta);
#endif
	    storeData.assign(mass,ee_rapidity,sc->scEt,sc->scEta, storeNGoodPV,

			     event_weight, 1.);
	  }
	  int templateBin = 
	    getTemplateBin( DYTools::findEtBin(sc->scEt,etBinning),
			    DYTools::findEtaBin(sc->scEta,etaBinning),
			    etaBinning);
	  if( electronMatch != 0 ){
	    // supercluster has match in reconstructed electrons: "pass"
	    hMassPass->Fill(mass);
	    passTree->Fill();
	    if(sample != DYTools::DATA && templateBin != -1)
	      hPassTemplateV[templateBin]->Fill(mass);
	  }else{
	    // supercluster is not reconstructed as an electron
	    hMassFail->Fill(mass);
	    failTree->Fill();
	    if(sample != DYTools::DATA && templateBin != -1)
	      hFailTemplateV[templateBin]->Fill(mass);
	  }
	  
	  } // end loop over superclusters - probes	  
      } // end loop over electrons - tags      
    } // end loop over events
    
    delete infile;
    infile=0;
    eventTree=0;
    
    delete gen;
    delete info;
    delete eleArr;
    delete scArr;
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
      TString("/npv_tnp") + effTypeString + TString("_") + sampleTypeString +
      DYTools::analysisTag + TString(".root");
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
  printf("    events after event level MET cut                         %15d\n",eventsAfterMET);
  printf("\nTotal electron tag candidates (no cuts)                      %15d\n",tagCand);
  printf("                 tag candidates Et>20                        %15d\n",tagCandPassEt);
  printf("                 tag candidates, eta in acceptance           %15d\n",tagCandPassEta);
  printf("                 tag candidates, matched to GEN (if MC)      %15d\n",tagCandGenMatched);
  printf("                 tag candidates, ECAL driven                 %15d\n",tagCandEcalDriven);
  printf("                 tag candidates, full selection(ID,HLT)      %15d\n",tagCandFinalCount);

  printf("\nTotal tag(electron)-probe(supercluster) pairs                %15d\n",numTagProbePairs);
  printf("               probe Et>10                                   %15d\n",numTagProbePairsPassEt);
  printf("               probe eta in acceptance                       %15d\n",numTagProbePairsPassEta);
  printf("               probe matched to GEN (if MC)                  %15d\n",numTagProbePairsGenMatched);
  printf("               tag-probe mass in 60-120 GeV window           %15d\n",numTagProbePairsInMassWindow);
  printf("               probe passes SC trk isolation                 %15d\n",numTagProbePairsPassSCIso);

  printf("\nNumber of probes, total                                      %15.0f\n", hMassTotal->GetSumOfWeights());
  printf("Number of probes, passed                                     %15.0f\n", hMassPass->GetSumOfWeights());
  printf("Number of probes, failed                                     %15.0f\n", hMassFail->GetSumOfWeights());

  if (evaluate_efficiencies) {
  // Human-readable text file to store measured efficiencies
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

  // Fit log 
  TString fitlogname = 
    TString("results_unsorted/efficiency_TnP_")+label+TString("_fitlog.dat");
  ofstream fitLog;
  fitLog.open(fitlogname);

  //
  //  Find efficiency
  //
  bool useTemplates = false;
  if(sample == DYTools::DATA)
    useTemplates = true;

  int NsetBins=30;
//   bool isRECO=1;
  const char* setBinsType="cache";
  
  int nDivisions = 
    DYTools::getNEtBins(etBinning)*DYTools::getNEtaBins(etaBinning);
  double ymax = 800;
  if(nDivisions <4 )
    ymax = nDivisions * 200;
  else if (nDivisions>DYTools::maxTnPCanvasDivisions) {
    nDivisions=DYTools::maxTnPCanvasDivisions;
  }

  TCanvas *c1 = MakeCanvas("c1","c1", 600, (int)ymax);
  c1->Divide(2,nDivisions);
  measureEfficiencyPU(passTree, failTree,
		    calcMethod, etBinning, etaBinning, c1, effOutput, fitLog,
		      useTemplates, templatesFile, 
		      resrootBase,
		      //resultsRootFile, //resultsRootFilePlots,
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
  //c1->Update();
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
 
  gBenchmark->Show("eff_Reco");

  return retCodeOk;
}


