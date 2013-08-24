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



//=== FUNCTION DECLARATIONS ======================================================================================

//=== MAIN MACRO =================================================================================================

int plotDYAcceptance(const TString conf,
		     DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN,
		     DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST,
		     double FSRreweight=1.0, double FSRmassDiff=1.)
{
  gBenchmark->Start("plotDYAcceptance");

  {
    DYTools::printExecMode(runMode,systMode);
    const int debug_print=1;
    if (!DYTools::checkSystMode(systMode,debug_print,2, DYTools::NO_SYST, DYTools::FSR_STUDY)) 
      return retCodeError;
  }

  TString extraTag;
  if (systMode==DYTools::FSR_STUDY) {
    extraTag=Form("_FSR%3.2lf",FSRreweight);
    extraTag.ReplaceAll(".","");
    std::cout << "FSR study with extra factor=" << FSRreweight << ", if mass diff is below " << FSRmassDiff << "; extraTag=<" << extraTag << ">\n";
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
			      extraTag, EventSelector::_selectDefault);
  evtSelector.setTriggerActsOnData(false);

  // Event weight handler
  EventWeight_t evWeight;
  //evWeight.init(inpMgr.puReweightFlag(),inpMgr.fewzFlag());
  // ignore PU-reweight flag
  evWeight.init(0,inpMgr.fewzFlag());

  // Prepare output directory
  inpMgr.constDir(systMode,1);

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================
  
  //  
  // Set up histograms
  //

  // containers to accumulate the events
  std::vector<TH2D*> hvPass, hvTotal;
  std::vector<TH2D*> hvFail;
  std::vector<TH2D*> hvPassPreFsr, hvTotalPreFsr;

  // debug containters
  std::vector<TH1D*> hMassv;//, hZMass2v, hZPtv, hZPt2v, hZyv, hZPhiv;
  std::vector<TH1D*> hMassBinsv;
  std::vector<TH1D*> hZpeakv;
  TH1D *hSelEvents=NULL;
  
  // the main result of the macro
  createBaseH2Vec(hvPass ,"hvPass_" ,inpMgr.mcSampleNames());
  createBaseH2Vec(hvTotal,"hvTotal_",inpMgr.mcSampleNames());
  createBaseH2Vec(hvFail,"hvFail_",inpMgr.mcSampleNames());

  createBaseH2Vec(hvPassPreFsr ,"hvPassPreFsr_" ,inpMgr.mcSampleNames());
  createBaseH2Vec(hvTotalPreFsr,"hvTotalPreFsr_",inpMgr.mcSampleNames());

  // debug distributions: 1GeV bins
  //createAnyH1Vec(hMassv,"hMass_",inpMgr.sampleNames(),2500,0.,2500.,"M_{ee} [GeV]","counts/1GeV");
  createAnyH1Vec(hMassv,"hMass_",inpMgr.mcSampleNames(),1490,10.,1500.,"M_{ee} [GeV]","counts/1GeV");
  // debug distributions for current mass bin
  createBaseH1Vec(hMassBinsv,"hMassBins_",inpMgr.mcSampleNames());
  // debug: accumulate info about the selected events in the samples
  hSelEvents=createAnyTH1D("hSelEvents","hSelEvents",inpMgr.mcSampleCount(),0,inpMgr.mcSampleCount(),"sampleId","event count");
  // collect number of events in the Z-peak
  createAnyH1Vec(hZpeakv,"hZpeak_",inpMgr.mcSampleNames(),60,60.,120.,"M_{ee} [GeV]","counts/1GeV");

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
      TFile infile(mcSample->getFName(ifile),"read");
      assert(infile.IsOpen());

      // Get the TTrees
      if (!accessInfo.setTree(infile,"Events",true)) {
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
	ec.numEvents_inc();
	if (DYTools::isDebugMode(runMode) && (ientry>1000000)) break; // debug option
	//if (DYTools::isDebugMode(runMode) && (ientry>100)) break; // debug option
	printProgress(100000," ientry=",ientry,maxEvents);
	
	// Load generator level info
	accessInfo.GetGen(ientry);
	// If the Z->ll leptons are not electrons, discard this event.
	// This is needed for signal MC samples such as Madgraph Z->ll
	// where all 3 lepton flavors are possible
	if (!accessInfo.genLeptonsAreElectrons()) continue;

	// Load event info to get nPU
	accessInfo.GetInfoEntry(ientry);

	// FSR study correction for weight
	if (systMode==DYTools::FSR_STUDY) {
	  evWeight.setSpecWeightValue(accessInfo,FSRmassDiff,FSRreweight);
	}

	// Adjust event weight
	// .. here "false" = "not data"
	evWeight.set_PU_and_FEWZ_weights(accessInfo,false);

	// adjust the scale in the counter to include FEWZ 
	// (and possibly PU) weight
	//ec.setScale(evWeight.totalWeight());

	// accumulate denominator
	const mithep::TGenInfo *gen= accessInfo.genPtr();
	hvTotal[ifile]->Fill(gen->mass, fabs(gen->y), evWeight.totalWeight());
	hvTotalPreFsr[ifile]->Fill(gen->vmass, fabs(gen->vy), evWeight.totalWeight());

	int failAcc=0;

	// check acceptance at Gen level
	if (!evtSelector.inAcceptance(accessInfo)) {
	  hvFail[ifile]->Fill(gen->mass, fabs(gen->y), evWeight.totalWeight());
	  failAcc=1;
	}

	if (evtSelector.inAcceptancePreFsr(accessInfo)) {
	  hvPassPreFsr[ifile]->Fill(gen->vmass, fabs(gen->vy), evWeight.totalWeight());
	}
	if (failAcc) continue;

	ec.numEventsPassedAcceptance_inc();

	// accumulate denominator
	hvPass[ifile]->Fill(gen->mass, fabs(gen->y), evWeight.totalWeight());

      } // loop over events
      ec.print(2);  // print info about file
      ecSample.add(ec); // accumulate event counts
      ecTotal.add(ec);
      
      //delete infile;
      infile.Close();
      //infile=0; 
      //eventTree=0;
    }
    ecSample.print(2); // print info about sample
    evtSelector.printCounts();
  }
  ecTotal.print(2);
  } // if (processData)


  TString outFileName=inpMgr.correctionFullFileName("acceptance",systMode,0);
  std::cout << "outFileName=<" << outFileName << ">\n";
  if (DYTools::processData(runMode)) {
    TFile file(outFileName,"recreate");
    int res=file.IsOpen();
    if (res) res=saveVec(file,hvPass,"accPassDir");
    if (res) res=saveVec(file,hvTotal,"accTotalDir");
    if (res) res=saveVec(file,hvFail,"accFailDir");
    if (res) res=saveVec(file,hvPassPreFsr,"accPassPreFsrDir");
    if (res) res=saveVec(file,hvTotalPreFsr,"accTotalPreFsrDir");
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
    if (res) res=saveVec(file,hvPassPreFsr,"accPassPreFsrDir");
    if (res) res=saveVec(file,hvTotalPreFsr,"accTotalPreFsrDir");
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

  TH2D *hSumPass=createBaseH2("hSumPass","hSumPass",1);
  TH2D *hSumTotal=createBaseH2("hSumTotal","hSumTotal",1);
  TH2D *hSumFail=createBaseH2("hSumFail","hSumFail",1);
  addHistos(hSumPass,hvPass);
  addHistos(hSumTotal,hvTotal);
  addHistos(hSumFail,hvFail);

  TH2D *hAcc = createBaseH2("hAcc","hAcc",1);
  // We need the binomial error for the acceptance
  // eff=Pass/Tot, 
  // (dEff)^2= (1-eff)^2/T^2 (dPass)^2 + eff^2/T^2 (dFail)^2
  // (dFail)^2 = (dTot)^2 - (dPass)^2
  hAcc->Divide(hSumPass,hSumTotal,1,1,"b");

  TH2D *hSumPassPreFsr=createBaseH2("hSumPassPreFsr","hSumPassPreFsr",1);
  TH2D *hSumTotalPreFsr=createBaseH2("hSumTotalPreFsr","hSumTotalPreFsr",1);
  TH2D *hAccPreFsr = createBaseH2("hAccPreFsr","hAccPreFsr",1);
  addHistos(hSumPassPreFsr,hvPassPreFsr);
  addHistos(hSumTotalPreFsr,hvTotalPreFsr);
  hAccPreFsr->Divide(hSumPassPreFsr,hSumTotalPreFsr,1,1,"b");

  std::cout << dashline;
  std::cout << dashline;

  printHisto(hSumPass);
  printHisto(hSumTotal);
  printHisto(hAcc);
  //printHisto(hSumFail);
  std::cout << dashline;
  printHisto(hSumPassPreFsr);
  printHisto(hSumTotalPreFsr);
  printHisto(hAccPreFsr);
  std::cout << dashline;

  if (DYTools::processData(runMode)) {
    TFile file(outFileName,"update");
    int res=file.IsOpen();
    std::cout << "res=" << res << "\n";
    if (res) res=saveHisto(file,hAcc,"");
    if (res) res=saveHisto(file,hSumPass,"");
    if (res) res=saveHisto(file,hSumTotal,"");
    if (res) res=saveHisto(file,hSumFail,"");
    if (res) res=saveHisto(file,hSumPassPreFsr,"");
    if (res) res=saveHisto(file,hSumTotalPreFsr,"");
    file.Close();
    if (!res) {
      std::cout << "error occurred during additional save to file <" << outFileName << ">\n";
      return retCodeError;
    }
  }

    /*



  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  

     
  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl; 
    */  

  gBenchmark->Show("plotDYAcceptance");
  return retCodeOk;
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
