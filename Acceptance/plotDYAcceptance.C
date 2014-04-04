#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TGraphErrors.h>           // graphs
#include <TProfile.h>               // profile histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for Lorentz vector computations
#include <TVectorD.h>
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <vector>                   // STL vector class
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "../Include/CPlot.hh"          // helper class for plots
#include "../Include/MitStyleRemix.hh"  // style settings for drawing
#include "../Include/MyTools.hh"        // miscellaneous helper functions
#include "../Include/DYTools.hh"
//#include "../Include/plotFunctions.hh"
//#include "../Include/latexPrintouts.hh"

// define classes and constants to read in ntuple
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TEventInfo.hh"
#include "../Include/TGenInfo.hh"
#include "../Include/TDielectron.hh"
#include "../Include/TriggerSelection.hh"
#include "../Include/TVertex.hh"

// Helper functions for Electron ID selection
#include "../Include/EleIDCuts.hh"

//#include "../Include/FEWZ.hh"
//#include "../Include/PUReweight.hh"

#include "../Include/AccessOrigNtuples.hh"
#include "../Include/EventSelector.hh"
#include "../Include/EventWeight.hh"
#include "../Include/InputFileMgr.hh"

#endif


#define calculate_PreFsrAcc


//=== FUNCTION DECLARATIONS ======================================================================================

//=== MAIN MACRO =================================================================================================

int plotDYAcceptance(int analysisIs2D,
		     const TString conf,
		     DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN,
		     DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST)
{
  gBenchmark->Start("plotDYAcceptance");

  {
    DYTools::printExecMode(runMode,systMode);
    const int debug_print=1;
    if (!DYTools::checkSystMode(systMode,debug_print,5, DYTools::NO_SYST,
				DYTools::FSR_5plus, DYTools::FSR_5minus,
				DYTools::NO_REWEIGHT, DYTools::NO_REWEIGHT_FEWZ))
      return retCodeError;
  }
  //
  // A note on systematics mode
  // - FSR_5plus, FSR_5minus should be taken care here

  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return retCodeError;
  }


  //--------------------------------------------------------------------------------------------------------------
  // Settings
  //==============================================================================================================

  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(conf)) return retCodeError;
  // no energy correction for this evaluation
  inpMgr.clearEnergyScaleTag();

  // Construct eventSelector, update mgr and plot directory
  EventSelector_t evtSelector(inpMgr,runMode,systMode,
			      "","",EventSelector::_selectDefault);
  evtSelector.setTriggerActsOnData(false);

  // Acceptance is generator-level quantity and should not
  // depend on the pile-up. The flag is disabled.
  EventWeight_t evWeight;
  evWeight.init(0*inpMgr.puReweightFlag(),inpMgr.fewzFlag(),systMode);
  std::cout << "evWeight: "; evWeight.PrintDetails();

  int useSpecWeight=1;
  double specWeight=1.0;
  const double FSRmassDiff=1.0;
  if (systMode==DYTools::FSR_5plus) { specWeight=1.05; useSpecWeight=1; }
  else if (systMode==DYTools::FSR_5minus) { specWeight=0.95; useSpecWeight=1; }

  // Prepare output directory
  inpMgr.constDir(systMode,1);
  TString outFileName=inpMgr.correctionFullFileName("acceptance",systMode,0);
  std::cout << "generated outFileName=<" << outFileName << ">\n";

#ifdef calculate_PreFsrAcc
  std::cout << dashline;
  std::cout << "\t\tcalculate_PreFsrAcc is defined\n";
  std::cout << dashline;
#endif

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code
  //==============================================================================================================

  std::cout << mainpart;

  //
  // Set up histograms
  //

  // containers to accumulate the events
  std::vector<TH2D*> hvPass, hvTotal;
  std::vector<TH2D*> hvFail;

#ifdef calculate_PreFsrAcc
  std::vector<TH2D*> hvPassPreFsr, hvTotalPreFsr;
#endif

  // debug containters
  std::vector<TH1D*> hMassv;
  std::vector<TH1D*> hMassBinsv;
  std::vector<TH1D*> hZpeakv;
  TH1D *hSelEvents=NULL;

  // the main result of the macro
  createBaseH2Vec(hvPass ,"hvPass_" ,inpMgr.mcSampleNames());
  createBaseH2Vec(hvTotal,"hvTotal_",inpMgr.mcSampleNames());
  createBaseH2Vec(hvFail,"hvFail_",inpMgr.mcSampleNames());

#ifdef calculate_PreFsrAcc
  createBaseH2Vec(hvPassPreFsr ,"hvPassPreFsr_" ,inpMgr.mcSampleNames());
  createBaseH2Vec(hvTotalPreFsr,"hvTotalPreFsr_",inpMgr.mcSampleNames());
#endif

  // debug distributions: 1GeV bins
  createAnyH1Vec(hMassv,"hGenMass_",inpMgr.mcSampleNames(),1990,10.,2000.,"#it{M}_{ee} [GeV]","counts/1GeV");
  // debug distributions for current mass bin
  createBaseH1Vec(hMassBinsv,"hGenMassBins_",inpMgr.mcSampleNames());
  // debug: accumulate info about the selected events in the samples
  hSelEvents=createAnyTH1D("hSelEvents","hSelEvents",inpMgr.mcSampleCount(),0,inpMgr.mcSampleCount(),"sampleId","event count");
  // collect number of events in the Z-peak
  createAnyH1Vec(hZpeakv,"hZpeak_",inpMgr.mcSampleNames(),60,60.,120.,"#it{M}_{ee} [GeV]","counts/1GeV");

  //
  // Access samples and fill histograms
  //
  AccessOrigNtuples_t accessInfo;

  //
  // loop over samples
  //
  if (DYTools::processData(runMode)) {

  double extraWeightFactor=1.0;
  EventCounterExt_t ecTotal("total");
  for (unsigned int isample=0; isample<inpMgr.mcSampleCount(); ++isample) {
    const CSample_t *mcSample=inpMgr.mcSampleInfo(isample);
    std::cout << "Processing " << mcSample->getLabel() << "..." << std::endl;
    std::cout << " of size " << mcSample->size() << "\n";
    if (mcSample->size()!=1) {
      std::cout << "mcSample->size is expected to be 1\n";
      return retCodeError;
    }

    // accumulate info about processed files
    EventCounterExt_t ecSample(mcSample->name);

    for (unsigned int ifile=0; ifile<mcSample->size(); ++ifile) {
      // Read input file
      TFile *infile=new TFile(mcSample->getFName(ifile),"read");
      if (!infile || !infile->IsOpen()) {
	TString skimName=inpMgr.convertSkim2Ntuple(mcSample->getFName(ifile));
	std::cout <<  "  .. failed. Trying <" << skimName << ">" << std::endl;
	infile= new TFile(skimName,"read");
      }
      assert(infile->IsOpen());
      std::cout << " Reading file <" << mcSample->getFName(ifile) << ">\n";

      // Get the TTrees
      if (!accessInfo.setTree(*infile,"Events",true)) {
	return retCodeError;
      }

    // Find weight for events for this file
    // The first file in the list comes with weight 1*extraWeightFactor,
    // all subsequent ones are normalized to xsection and luminosity
      ULong_t maxEvents = accessInfo.getEntries();
      // to match old version package (DYee 7TeV paper),
      if ((inpMgr.userKeyValueAsInt("USE7TEVMCWEIGHT")==1) && (isample==0) && (ifile==0)) {
	extraWeightFactor=maxEvents / (inpMgr.totalLumi() * inpMgr.mcSampleInfo(0)->getXsec(ifile));
      }
      //std::cout << "extraWeightFactor=" << extraWeightFactor << ", chk=" << (maxEvents0/inpMgr.mcSampleInfo(0)->getXsec(ifile)) << "\n";
      //const double extraWeightFactor=1.0;
      if (! evWeight.setWeight_and_adjustMaxEvents(maxEvents, inpMgr.totalLumi(), mcSample->getXsec(ifile),
						   extraWeightFactor, inpMgr.selectEventsFlag())) {
	std::cout << "adjustMaxEvents failed\n";
	return retCodeError;
      }
      std::cout << "mcSample xsec=" << mcSample->getXsec(ifile) << ", nEntries=" << maxEvents << "\n";


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
	if (DYTools::isDebugMode(runMode) && (ientry>1000000)) break; // debug option
	//if (DYTools::isDebugMode(runMode) && (ientry>100)) break; // debug option
	printProgress(250000," ientry=",ientry,maxEvents);
	ec.numEvents_inc();
	
	// Load generator level info
	accessInfo.GetGen(ientry);
	// If the Z->ll leptons are not electrons, discard this event.
	// This is needed for signal MC samples such as Madgraph Z->ll
	// where all 3 lepton flavors are possible
	if (!accessInfo.genLeptonsAreElectrons()) continue;

	// Load event info to get nPU
	accessInfo.GetInfoEntry(ientry);

	// FSR study correction for weight
	if (useSpecWeight) {
	  evWeight.setSpecWeightValue(accessInfo,FSRmassDiff,specWeight);
	}

	// Adjust event weight
	// .. here "false" = "not data"
	evWeight.set_PU_and_FEWZ_weights(accessInfo,false);

	// adjust the scale in the counter to include FEWZ
	// (and possibly PU) weight
	//ec.setScale(evWeight.totalWeight());

	// accumulate denominator
	const mithep::TGenInfo *gen= accessInfo.genPtr();
	hMassv[isample]->Fill(gen->vmass, evWeight.totalWeight());
	hMassBinsv[isample]->Fill(gen->vmass, evWeight.totalWeight());

	hvTotal[isample]->Fill(gen->mass, fabs(gen->y), evWeight.totalWeight());

#ifdef calculate_PreFsrAcc
	hvTotalPreFsr[isample]->Fill(gen->vmass, fabs(gen->vy), evWeight.totalWeight());
#endif

	int failAcc=0;

	// check acceptance at Gen PostFSR level
	if (!evtSelector.inAcceptance(accessInfo)) {
	  hvFail[isample]->Fill(gen->mass, fabs(gen->y), evWeight.totalWeight());
	  failAcc=1;
	}

#ifdef calculate_PreFsrAcc
	if (evtSelector.inAcceptancePreFsr(accessInfo)) {
	  hvPassPreFsr[isample]->Fill(gen->vmass, fabs(gen->vy), evWeight.totalWeight());
	}
#endif

	if (failAcc) continue;

	ec.numEventsPassedAcceptance_inc();

	// accumulate denominator
	hvPass[isample]->Fill(gen->mass, fabs(gen->y), evWeight.totalWeight());
	hSelEvents->Fill(isample,evWeight.totalWeight());

      } // loop over events
      ec.print(2);  // print info about file
      ecSample.add(ec); // accumulate event counts
      ecTotal.add(ec);

      infile->Close();
      delete infile;
    }
    ecSample.print(2); // print info about sample
    evtSelector.printCounts();
  }
  ecTotal.print(2);
  } // if (processData)


  std::cout << "outFileName=<" << outFileName << ">\n";
  if (DYTools::processData(runMode)) {
    TFile file(outFileName,"recreate");
    int res=file.IsOpen();
    if (res) res=saveVec(file,hvPass,"accPassDir");
    if (res) res=saveVec(file,hvTotal,"accTotalDir");
    if (res) res=saveVec(file,hvFail,"accFailDir");
#ifdef calculate_PreFsrAcc
    if (res) res=saveVec(file,hvPassPreFsr,"accPassPreFsrDir");
    if (res) res=saveVec(file,hvTotalPreFsr,"accTotalPreFsrDir");
#endif
    if (res) res=saveVec(file,hMassv,"mass_1GeV_bins");
    if (res) res=saveVec(file,hMassBinsv,"mass_analysis_bins");
    if (res) res=saveVec(file,hZpeakv,"mass_Zpeak_1GeV");
    if (res) res=saveHisto(file,hSelEvents,"procFileInfo");
    if (res) writeBinningArrays(file);
    file.Close();
    if (!res) {
      std::cout << "error occurred during save to file <" << outFileName << ">\n";
      return retCodeError;
    }
  }
  else {
    TFile file(outFileName,"read");
    int res=file.IsOpen();
    if (res) res=checkBinningArrays(file);
    if (res) res=loadVec(file,hvPass,"accPassDir");
    if (res) res=loadVec(file,hvTotal,"accTotalDir");
    if (res) res=loadVec(file,hvFail,"accFailDir");
#ifdef calculate_PreFsrAcc
    if (res) res=saveVec(file,hvPassPreFsr,"accPassPreFsrDir");
    if (res) res=saveVec(file,hvTotalPreFsr,"accTotalPreFsrDir");
#endif
    if (res) res=loadVec(file,hMassv,"mass_1GeV_bins");
    if (res) res=loadVec(file,hMassBinsv,"mass_analysis_bins");
    if (res) res=loadVec(file,hZpeakv,"mass_Zpeak_1GeV");
    if (res) res=loadHisto(file,&hSelEvents,"procFileInfo");
    file.Close();
    if (!res) {
      std::cout << "error occurred during save to file <" << outFileName << ">\n";
      return retCodeError;
    }
  }

  TH2D *hSumPass_BaseH2=createBaseH2("hSumPass_baseH2","hSumPass_baseH2",1);
  TH2D *hSumTotal_BaseH2=createBaseH2("hSumTotal_baseH2","hSumTotal_baseH2",1);
  TH2D *hSumFail_BaseH2=createBaseH2("hSumFail_baseH2","hSumFail_baseH2",1);
  addHistos(hSumPass_BaseH2,hvPass);
  addHistos(hSumTotal_BaseH2,hvTotal);
  addHistos(hSumFail_BaseH2,hvFail);

  TH2D *hAcc_BaseH2 = createBaseH2("hAcceptance_baseH2","hAcc_baseH2",1);
  // We need the binomial error for the acceptance
  // eff=Pass/Tot,
  // (dEff)^2= (1-eff)^2/T^2 (dPass)^2 + eff^2/T^2 (dFail)^2
  // (dFail)^2 = (dTot)^2 - (dPass)^2
  hAcc_BaseH2->Divide(hSumPass_BaseH2,hSumTotal_BaseH2,1,1,"b");

  TH2D *hSumPass=convertBaseH2actual(hSumPass_BaseH2,"hSumPass",1);
  TH2D *hSumFail=convertBaseH2actual(hSumFail_BaseH2,"hSumFail",1);
  TH2D *hSumTotal=convertBaseH2actual(hSumTotal_BaseH2,"hSumTotal",1);
  TH2D *hAcc=Clone(hAcc_BaseH2,"hAcceptance","hAcc");
  hAcc->Divide(hSumPass,hSumTotal,1,1,"b");

#ifdef calculate_PreFsrAcc
  TH2D *hSumPassPreFsr_BaseH2=createBaseH2("hSumPassPreFsr_baseH2","hSumPassPreFsr_baseH2",1);
  TH2D *hSumTotalPreFsr_BaseH2=createBaseH2("hSumTotalPreFsr_baseH2","hSumTotalPreFsr_baseH2",1);
  addHistos(hSumPassPreFsr_BaseH2,hvPassPreFsr);
  addHistos(hSumTotalPreFsr_BaseH2,hvTotalPreFsr);

  TH2D *hSumPassPreFsr=convertBaseH2actual(hSumPassPreFsr_BaseH2,"hSumPassPreFsr",1);
  TH2D *hSumTotalPreFsr=convertBaseH2actual(hSumTotalPreFsr_BaseH2,"hSumTotalPreFsr",1);
  TH2D *hAccPreFsr=Clone(hSumPassPreFsr,"hAccPreFsr","hAccPreFsr");
  hAccPreFsr->Divide(hSumPassPreFsr,hSumTotalPreFsr,1,1,"b");
#endif

  std::cout << dashline;
  std::cout << dashline;

  printHisto(hSumPass_BaseH2);
  printHisto(hSumTotal_BaseH2);
  printHisto(hAcc_BaseH2);

  printHisto(hSumPass);
  printHisto(hSumTotal);
  printHisto(hAcc);

  //printHisto(hSumFail);
  std::cout << dashline;

#ifdef calculate_PreFsrAcc
  std::cout << dashline;
  printHisto(hSumPassPreFsr);
  printHisto(hSumTotalPreFsr);
  printHisto(hAccPreFsr);
  std::cout << dashline;
#endif

  if (DYTools::processData(runMode)) {
    TFile file(outFileName,"update");
    int res=file.IsOpen();
    std::cout << "res=" << res << "\n";
    if (res) res=saveHisto(file,hAcc,"");
    if (res) res=saveHisto(file,hSumPass,"");
    if (res) res=saveHisto(file,hSumTotal,"");
    if (res) res=saveHisto(file,hSumFail,"");
#ifdef calculate_PreFsrAcc
    if (res) res=saveHisto(file,hAccPreFsr,"");
    if (res) res=saveHisto(file,hSumPassPreFsr,"");
    if (res) res=saveHisto(file,hSumTotalPreFsr,"");
#endif
    file.Close();
    if (!res) {
      std::cout << "error occurred during additional save to file <" << outFileName << ">\n";
      return retCodeError;
    }
  }


  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  /*
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
    */

  //gBenchmark->Show("plotDYAcceptance");
  ShowBenchmarkTime("plotDYAcceptance");
  return retCodeOk;
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
