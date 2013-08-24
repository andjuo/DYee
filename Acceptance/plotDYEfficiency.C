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

int plotDYEfficiency(const TString conf,
		  DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN,
		  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST)
{
  gBenchmark->Start("plotDYEfficiency");

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
  if (!inpMgr.Load(conf)) return retCodeError;
  // no energy correction for this evaluation
  inpMgr.clearEnergyScaleTag();

  // Construct eventSelector, update mgr and plot directory
  EventSelector_t evtSelector(inpMgr,runMode,systMode,
			      "", EventSelector::_selectDefault);
  evtSelector.setTriggerActsOnData(false);

  // Event weight handler
  EventWeight_t evWeight;
  evWeight.init(inpMgr.puReweightFlag(),inpMgr.fewzFlag());

  // Prepare output directory
  inpMgr.constDir(systMode,1);

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================
  
  //  
  // Set up histograms
  //

  // containers to accumulate the events
  std::vector<TH2D*> hPass, hTotal;
  //std::vector<TH2D*> hFail;

  // debug containters
  std::vector<TH1D*> hMassv;//, hZMass2v, hZPtv, hZPt2v, hZyv, hZPhiv;
  std::vector<TH1D*> hMassBinsv;
  std::vector<TH1D*> hZpeakv;
  TH1D *hSelEvents=NULL;
  
  // the main result of the macro
  createBaseH2Vec(hPass ,"hPass_" ,inpMgr.mcSampleNames());
  createBaseH2Vec(hTotal,"hTotal_",inpMgr.mcSampleNames());
  //createBaseH2Vec(hFail,"hFail_",inpMgr.mcSampleNames());

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
      if ((isample==0) && (ifile==0)) {
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

	// check acceptance at Gen level
	if (!evtSelector.inAcceptance(accessInfo)) continue;
	ec.numEventsPassedAcceptance_inc();

	// Load event info
	accessInfo.GetInfoEntry(ientry);

	// Adjust event weight
	// .. here "false" = "not data"
	evWeight.set_PU_and_FEWZ_weights(accessInfo,false);

	// adjust the scale in the counter to include FEWZ 
	// (and possibly PU) weight
	//ec.setScale(evWeight.totalWeight());

	// accumulate denominator
	const mithep::TGenInfo *gen= accessInfo.genPtr();
	hTotal[ifile]->Fill(gen->mass, fabs(gen->y), evWeight.totalWeight());

	// check event trigger
	if (!evtSelector.eventTriggerOk(accessInfo)) {
	  //hFail[ifile]->Fill(gen->mass, fabs(gen->y), evWeight.totalWeight());
	  continue; // no trigger accept? Skip to next event...	
	}
	ec.numEventsPassedEvtTrigger_inc();

	// load dielectron array
	accessInfo.GetDielectrons(ientry);

	// loop through dielectrons
	//int pass=0;
	for(Int_t i=0; i<accessInfo.dielectronCount(); i++) {
	  mithep::TDielectron *dielectron = accessInfo.editDielectronPtr(i);
	  ec.numDielectrons_inc();
	  
	  // escale may modify dielectron! But it should not here
	  if (!evtSelector.testDielectron(dielectron,accessInfo.evtInfoPtr(),&ec)) continue;
	  //pass=1;


          /******** We have a Z candidate! HURRAY! ********/
	  
	  ec.numDielectronsPass_inc();
	  if (ec.numDielectronsOkSameSign_inc(dielectron->q_1,dielectron->q_2)) {
	    // same sign event
	  }

	  hPass[ifile]->Fill(gen->mass, fabs(gen->y), evWeight.totalWeight());
	  hSelEvents->Fill(ifile,evWeight.totalWeight());

	} // loop over dielectrons
	//if (!pass) {
	//  hFail[ifile]->Fill(gen->mass, fabs(gen->y), evWeight.totalWeight());
	//}
      } // loop over events
      ec.print();  // print info about file
      ecSample.add(ec); // accumulate event counts
      ecTotal.add(ec);
      
      //delete infile;
      infile.Close();
      //infile=0; 
      //eventTree=0;
    }
    ecSample.print(); // print info about sample
    evtSelector.printCounts();
  }
  ecTotal.print();
  } // if (processData)


  TString outFileName=inpMgr.correctionFullFileName("efficiency",systMode,0);
  if (DYTools::processData(runMode)) {
    TFile file(outFileName,"recreate");
    int res=1;
    if (res) res=saveVec(file,hPass,"effPassDir");
    if (res) res=saveVec(file,hTotal,"effTotalDir");
    //if (res) res=saveVec(file,hFail,"effFailDir");
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
    if (res) res=loadVec(file,hPass,"effPassDir");
    if (res) res=loadVec(file,hTotal,"effTotalDir");
    //if (res) res=loadVec(file,hFail,"effFailDir");
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
  //TH2D *hSumFail=createBaseH2("hSumFail","hSumFail",1);
  addHistos(hSumPass,hPass);
  addHistos(hSumTotal,hTotal);
  //addHistos(hSumFail,hFail);

  TH2D *hEff = createBaseH2("hEff","hEff",1);
  // We need the binomial error for the efficiency:
  // eff=Pass/Tot, 
  // (dEff)^2= (1-eff)^2/T^2 (dPass)^2 + eff^2/T^2 (dFail)^2
  // (dFail)^2 = (dTot)^2 - (dPass)^2
  hEff->Divide(hSumPass,hSumTotal,1,1,"b");

  //TH2D *hEffChk=Clone(hEff,

  printHisto(hSumPass);
  printHisto(hSumTotal);
  printHisto(hEff);
  //printHisto(hSumFail);

    /*

      // The A_FSR starting point: gen level quantities
      // The FSR subscript means we work with post-FSR generator level quantities
      double et1 = gen->pt_1;
      double et2 = gen->pt_2;
      double eta1 = gen->eta_1;
      double eta2 = gen->eta_2;

      // Apply acceptance requirement
      if( ! DYTools::goodEtEtaPair( et1, eta1, et2, eta2 ) ) continue;

      // These events are in acceptance, use them for efficiency denominator
      Bool_t isBGen1 = DYTools::isBarrel(eta1);
      Bool_t isBGen2 = DYTools::isBarrel(eta2);
      // determine number of good vertices
      pvBr->GetEntry(ientry);
      int iPUBin=-1;
      int nGoodPV=-1;
      double puWeight=1.0;
	nGoodPV=0;
	const int new_pv_count_code=1;
	if (new_pv_count_code) {
	  nGoodPV=countGoodVertices(pvArr);
	}
	else {
	for (Int_t ipv=0; ipv<pvArr->GetEntriesFast(); ipv++) {
	  const mithep::TVertex *pv = (mithep::TVertex*)((*pvArr)[ipv]);
	  if(pv->nTracksFit                        < 1)  continue;
	  if(pv->ndof                              < 4)  continue;
	  if(fabs(pv->z)                           > 24) continue;
	  if(sqrt((pv->x)*(pv->x)+(pv->y)*(pv->y)) > 2)  continue;
	  nGoodPV++;
	}
	}
      if ((gen->mass>=60) && (gen->mass<=120)) {
	if (nGoodPV>0) {
           if (nGoodPV<=nEventsZPeakPURaw.GetNoElements()) 
	        nEventsZPeakPURaw[nGoodPV] += scale * gen->weight;
	   iPUBin=DYTools::findPUBin(nGoodPV);
	   //std::cout << "iPUBin=" << iPUBin << ", nGoodPV=" << nGoodPV << "\n";
	   if ((iPUBin!=-1) && (iPUBin < nEventsZPeakPU.GetNoElements())) {
	     nEventsZPeakPU[iPUBin] += scale * gen->weight;
	   }
	   //else {
	   //  std::cout << "error in PU bin indexing iPUBin=" << iPUBin << ", nGoodPV=" << nGoodPV << "\n";
	   //}
	}
	//else std::cout << "nGoodPV=" << nGoodPV << "\n";
      }

#ifdef usePUReweight
      puWeight = puReweight.getWeightHildreth(info->nPU);
      nZv_puUnweighted += scale * gen->weight;
      nZv_puWeighted += scale * gen->weight * puWeight;
#endif
    
      // Use post-FSR generator level mass for binning
      int ibinGenM = DYTools::findMassBin(gen->mass);
      int ibinGenY = DYTools::findAbsYBin(ibinGenM,gen->y);
      double totalWeight= scale * gen->weight * puWeight;
      if (useFewzWeights) totalWeight *= fewz.getWeight(gen->vmass,gen->vpt,gen->vy);

      // Accumulate denominator for efficiency calculations
      if(ibinGenM != -1 && ibinGenY != -1 && ibinGenM < DYTools::nMassBins && ibinGenY < DYTools::nYBins[ibinGenM]){
	nEventsv(ibinGenM,ibinGenY) += totalWeight;
        sumWeightsTotaSq(ibinGenM,ibinGenY) += totalWeight*totalWeight;
	// Split events barrel/endcap using matched supercluster or particle eta
	if(isBGen1 && isBGen2)                                  { nEventsBBv(ibinGenM,ibinGenY) += totalWeight; } 
	else if(!isBGen1 && !isBGen2)                           { nEventsEEv(ibinGenM,ibinGenY) += totalWeight; } 
	else if((isBGen1 && !isBGen2) || (!isBGen1 && isBGen2)) { nEventsBEv(ibinGenM,ibinGenY) += totalWeight; }
      }else
        binProblem++;

      // The line below is replaced by the superseeding method, can be cleaned up
      // if(!(info->triggerBits & eventTriggerBit)) continue;  // no trigger accept? Skip to next event...                                   

      if( !(requiredTriggers.matchEventTriggerBit(info->triggerBits, 
						  info->runNum))) 
	continue;
      
      // loop through dielectrons
      dielectronArr->Clear();
      dielectronBr->GetEntry(ientry);

      for(Int_t i=0; i<dielectronArr->GetEntriesFast(); i++) {
        const mithep::TDielectron *dielectron = (mithep::TDielectron*)((*dielectronArr)[i]);
	
	// Apply selection
	// Eta cuts and Et cuts
	if( ! DYTools::goodEtEtaPair( dielectron->scEt_1, dielectron->scEta_1,
				      dielectron->scEt_2, dielectron->scEta_2 ) ) continue;
        Bool_t isB1 = DYTools::isBarrel(dielectron->scEta_1);
        Bool_t isB2 = DYTools::isBarrel(dielectron->scEta_2);
	
	if( !( (isB1 == isBGen1 && isB2 == isBGen2 ) 
	       || (isB1 == isBGen2 && isB2 == isBGen1 ) ) )
	  countMismatch++;

	// Both electrons must match trigger objects. At least one ordering
	// must match
	if( ! requiredTriggers.matchTwoTriggerObjectsAnyOrder( dielectron->hltMatchBits_1,
							       dielectron->hltMatchBits_2,
							       info->runNum) ) continue;

	// The clause below can be deleted, it is superseeded by new methods
//  	if( ! ( 
//  	       (dielectron->hltMatchBits_1 & leadingTriggerObjectBit && 
//  		dielectron->hltMatchBits_2 & trailingTriggerObjectBit )
//  	       ||
//  	       (dielectron->hltMatchBits_1 & trailingTriggerObjectBit && 
//  		dielectron->hltMatchBits_2 & leadingTriggerObjectBit ) ) ) continue;

	// *** Smurf ID is superseeded by new selection ***
// 	// The Smurf electron ID package is the same as used in HWW analysis
// 	// and contains cuts like VBTF WP80 for pt>20, VBTF WP70 for pt<10
// 	// with some customization, plus impact parameter cuts dz and dxy
//  	if(!passSmurf(dielectron)) continue;

	// The selection below is for the EGM working points from spring 2012
	// recommended for both 2011 and 2012 data
	if( DYTools::energy8TeV == 1){
	  if(!passEGMID2012(dielectron, WP_MEDIUM, info->rhoLowEta)) continue;
	}else{
	  if(!passEGMID2011(dielectron, WP_MEDIUM, info->rhoLowEta)) continue;
	}

        // ******** We have a Z candidate! HURRAY! ******** /

	hZMassv[ifile]->Fill(gen->mass,totalWeight);

	// DEBUG
// 	if(ibinGen == 12)
// 	  printf("Gen mass %f  reco mass %f  scEt_1= %f  scEt_2= %f  scEta_1= %f  scEta_2= %f\n",
// 		 gen->mass, dielectron->mass, dielectron->scEt_1, dielectron->scEt_2,
// 		 dielectron->scEta_1, dielectron->scEta_2);
	
	// Accumulate numerator for efficiency calculations
	if ((nGoodPV>0) && (iPUBin!=-1)) { // -1 may also indicate that the mass was not in Z-peak range
	  if ((nGoodPV>=0) && (nGoodPV<=nEventsZPeakPURaw.GetNoElements())) nPassZPeakPURaw[nGoodPV] += totalWeight;
	  if (iPUBin < nPassZPeakPU.GetNoElements()) {
	    nPassZPeakPU[iPUBin] += totalWeight;
	  }
	  //else {
	  //  std::cout << "error in PU bin indexing\n";
	  //}
	}
	if(ibinGenM != -1 && ibinGenY != -1 && ibinGenM < DYTools::nMassBins && ibinGenY < DYTools::nYBins[ibinGenM]){
	  nPassv(ibinGenM,ibinGenY) += totalWeight;
          sumWeightsPassSq(ibinGenM,ibinGenY) += totalWeight*totalWeight;
	  if(isB1 && isB2)                            { nPassBBv(ibinGenM,ibinGenY) += totalWeight; } 
	  else if(!isB1 && !isB2)                     { nPassEEv(ibinGenM,ibinGenY) += totalWeight; } 
	  else if((isB1 && !isB2) || (!isB1 && isB2)) { nPassBEv(ibinGenM,ibinGenY) += totalWeight; }
	}

      } // end loop over dielectrons
    } // end loop over events 
    delete infile;
    infile=0, eventTree=0;
  } // end loop over files
  delete gen;

  cout << "ERROR: binning problem (" << binProblem <<" events in ECAL gap)"<<endl;

  effv      = 0;
  effErrv   = 0;
  effBBv    = 0;
  effErrBBv = 0;
  effBEv    = 0;
  effErrBEv = 0;
  effEEv    = 0;
  effErrEEv = 0;
  for(int i=0; i<DYTools::nMassBins; i++)
    for(int j=0; j<DYTools::nYBins[i]; j++){
      if(nEventsv(i,j) != 0){
        double nPass, nFail, nPassErr, nFailErr;
        nPass=nPassv(i,j);
        nFail=nEventsv(i,j)-nPassv(i,j); 
        nPassErr=sqrt(sumWeightsPassSq(i,j));
        nFailErr=sqrt(sumWeightsTotaSq(i,j)-sumWeightsPassSq(i,j));
        effv(i,j) = nPassv(i,j)/nEventsv(i,j);
        //effErrv(i,j) = sqrt(effv(i,j)*(1-effv(i,j))/nEventsv(i,j));
        effErrv(i,j) = sqrt(( nFail*nFail * nPassErr*nPassErr + nPass*nPass * nFailErr*nFailErr)) / (nEventsv(i,j)*nEventsv(i,j));
      }
    
      if (nEventsBBv(i,j) != 0) {
        effBBv(i,j) = nPassBBv(i,j)/nEventsBBv(i,j);
        effErrBBv(i,j) = sqrt(effBBv(i,j)*(1-effBBv(i,j))/nEventsBBv(i,j));
      }

      if (nEventsBEv(i,j) != 0) {
        effBEv(i,j) = nPassBEv(i,j)/nEventsBEv(i,j);
        effErrBEv(i,j) = sqrt(effBEv(i,j)*(1-effBEv(i,j))/nEventsBEv(i,j));
      }
      
      if (nEventsEEv(i,j) != 0) {
        effEEv(i,j) = nPassEEv(i,j)/nEventsEEv(i,j);
        effErrEEv(i,j) = sqrt(effEEv(i,j)*(1-effEEv(i,j))/nEventsEEv(i,j));
      }
    };

  effZPeakPU=0; effErrZPeakPU=0;
  for (int i=0; i<DYTools::nPVBinCount; ++i) {
    effZPeakPU[i]= nPassZPeakPU[i]/nEventsZPeakPU[i];
    effErrZPeakPU[i]= sqrt(effZPeakPU[i]*(1-effZPeakPU[i])/nEventsZPeakPU[i]);
  }

  printf("Sanity check: gen vs reco barrel-endcap assignment mismatches: %d\n",countMismatch);

  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  

  // destination dir
  TString extraTag= DYTools::analysisTag;
  if (!useFewzWeights) extraTag.Append("_noFEWZ");
  CPlot::sOutDir="plots" + extraTag;
  TString outputDir(TString("../root_files/constants/")+dirTag);
  gSystem->mkdir(outputDir,kTRUE);
  TString fnamePlots=outputDir + TString("/event_efficiency_plots") + extraTag + TString(".root");
  TFile *filePlots=new TFile(fnamePlots,"recreate");
  if (!filePlots || !filePlots->IsOpen()) {
    std::cout << "failed to create a file <" << fnamePlots << ">\n";
    throw 2;
  }

  TCanvas *c = MakeCanvas("canvEfficiency","canvEfficiency",800,600);

  // string buffers
  char ylabel[50];   // y-axis label

  // Z mass
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hZMassv[0]->GetBinWidth(1));
  CPlot plotZMass1("zmass1","","m(Z) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) { 
    plotZMass1.SetYRange(1, 10000000);
    plotZMass1.AddHist1D(hZMassv[i],labelv[i],"hist",colorv[i],linev[i]); 
  }
  plotZMass1.SetLogy();
  plotZMass1.Draw(c);
  SaveCanvas(c, "zmass1");

  PlotMatrixVariousBinning(effv, "efficiency", "LEGO2", filePlots);
  filePlots->Close();
  if (DYTools::study2D==0)
    Plot1D(effv,effErrv,"efficiency1D","efficiency");

  // Store constants in the file
  //TString effConstFileName(outputDir+TString("/event_efficiency_constants.root"));
  TString effConstFileName(outputDir+TString("/event_efficiency_constants") + extraTag + TString(".root"));

  TMatrixD nEventsvErr=nEventsv;
  TMatrixD nEventsBBvErr=nEventsBBv;
  TMatrixD nEventsBEvErr=nEventsBEv;
  TMatrixD nEventsEEvErr=nEventsEEv;
  TMatrixD nPassBBvErr=nPassBBv;
  TMatrixD nPassBEvErr=nPassBEv;
  TMatrixD nPassEEvErr=nPassEEv;

  TMatrixD nPassvErr=nPassv;
  for (int i=0; i<DYTools::nMassBins; ++i) {
    for (int j=0; j<DYTools::nYBins[i]; ++j) {
      nEventsvErr(i,j)=sqrt(sumWeightsTotaSq(i,j));
    }
  }
  for (int i=0; i<DYTools::nMassBins; ++i) {
    for (int j=0; j<DYTools::nYBins[i]; ++j) {
      nPassvErr(i,j)=sqrt(sumWeightsPassSq(i,j));
    }
  }
  nEventsBBvErr=0;
  nEventsBEvErr=0;
  nEventsEEvErr=0;
  nPassBBvErr=0;
  nPassBEvErr=0;
  nPassEEvErr=0;

   TFile fa(effConstFileName,"recreate");
   effv.Write("efficiencyArray");
   effErrv.Write("efficiencyErrArray");

   /
   nPassv.Write("effEval_nPass");
   nPassvErr.Write("effEval_nPassErr");
   nEventsv.Write("effEval_nTotal");
   nEventsvErr.Write("effEval_nTotalErr");

   nEventsBBv.Write("effEval_nTotalBB");
   nEventsBBvErr.Write("effEval_nTotalBBErr");
   nEventsBEv.Write("effEval_nTotalBE");
   nEventsBEvErr.Write("effEval_nTotalBEErr");
   nEventsEEv.Write("effEval_nTotalEE");
   nEventsEEvErr.Write("effEval_nTotalEEErr");
   nPassBBv.Write("effEval_nPassBB");
   nPassBBvErr.Write("effEval_nPassBBErr");
   nPassBEv.Write("effEval_nPassBE");
   nPassBEvErr.Write("effEval_nPassBEErr");
   nPassEEv.Write("effEval_nPassEE");
   nPassEEvErr.Write("effEval_nPassEEErr");

   effBBv.Write("efficiencyBB");
   effErrBBv.Write("efficiencyBBErr");
   effBEv.Write("efficiencyBE");
   effErrBEv.Write("efficiencyBEErr");
   effEEv.Write("efficiencyEE");
   effErrEEv.Write("efficiencyEEErr");
   /
   //effZPeakPU.Write("efficiencyZPeakPUArray");
   //effErrZPeakPU.Write("efficiencyErrZPeakPUArray");
   //nEventsZPeakPU.Write("nEventsZPeakPUArray");
   //nPassZPeakPU.Write("nPassZPeakPUArray");
   //nEventsZPeakPURaw.Write("nEventsZPeakPURawArray");
   //nPassZPeakPURaw.Write("nPassZPeakPURawArray");
   //unfolding::writeBinningArrays(fa);
   fa.Close();
     
  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl; 
  
  cout << labelv[0] << " file: " << fnamev[0] << endl;
  printf("     Number of generated events: %8.1lf\n",nZv);
#ifdef usePUReweight
  printf("     Number of selected events (puUnweighted): %8.1lf\n",nZv_puUnweighted);
  printf("     Number of selected events (puWeighted)  : %8.1lf\n",nZv_puWeighted);
#endif

  const char *yRangeStr=(DYTools::study2D) ? "rapidity range" : "";
  printf(" mass range  %s preselected      passed     total_Eff        BB-BB_Eff        EB-BB_Eff        EB-EB_Eff\n",yRangeStr);
  for(int i=0; i<DYTools::nMassBins; i++){
    double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
      printf(" %4.0f-%4.0f ", DYTools::massBinLimits[i], DYTools::massBinLimits[i+1]);
      if (DYTools::study2D!=0) printf(" %4.2f-%4.2f ", rapidityBinLimits[yi], rapidityBinLimits[yi+1]);
      printf("    %10.0f   %10.0f   %7.4f+-%6.4f  %7.4f+-%6.4f  %7.4f+-%6.4f  %7.4f+-%6.4f \n",
	     nEventsv(i,yi), nPassv(i,yi),
	     effv(i,yi), effErrv(i,yi),
	     effBBv(i,yi), effErrBBv(i,yi),
	     effBEv(i,yi), effErrBEv(i,yi),
	     effEEv(i,yi), effErrEEv(i,yi));
    }
    delete rapidityBinLimits;
  }

  if (1) {
    printf(" mass_range  %s total_preselected  total_passed     BB-BB_preselected  passed    EB-BB_preselected  passed        EB-EB_preselected  passed\n",yRangeStr);
    for(int i=0; i<DYTools::nMassBins; i++){
      double *rapidityBinLimits=DYTools::getYBinLimits(i);
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	printf(" %4.0f-%4.0f ", DYTools::massBinLimits[i], DYTools::massBinLimits[i+1]);
	if (DYTools::study2D!=0) printf(" %4.2f-%4.2f ", rapidityBinLimits[yi], rapidityBinLimits[yi+1]);
	printf("    %8.1f %8.1f   %8.1f %8.1f  %8.1f %8.1f  %8.1f %8.1f\n",
	       nEventsv(i,yi), nPassv(i,yi),
	       nEventsBBv(i,yi), nPassBBv(i,yi),
	       nEventsBEv(i,yi), nPassBEv(i,yi),
	       nEventsEEv(i,yi), nPassEEv(i,yi));
      }
      delete rapidityBinLimits;
    }
  }

  //printout to the txtFile in latex format
  if (DYTools::study2D)
    {
      latexPrintoutEfficiency2D(effv,effErrv,"Efficiency/plotDYEfficiency.C");
    }
  else if (DYTools::study2D==0)
    {
      latexPrintoutEfficiency1D(effv,effErrv,"Efficiency/plotDYEfficiency.C");
    }

  cout << endl;
  std::cout<<"printout in the Latex format is saved to the text file"<<std::endl;

  printf("\n\nZ-peak efficiency\n");
  printf(" PU bin    preselected      passed     total_Eff\n");
  for(int i=0; i<DYTools::nPVBinCount; i++){
    printf(" %4.0f-%4.0f   %10.0f   %10.0f   %7.4f+-%6.4f\n",
	   DYTools::nPVLimits[i], DYTools::nPVLimits[i+1],
	   nEventsZPeakPU[i], nPassZPeakPU[i],
	   effZPeakPU[i], effErrZPeakPU[i]);
  }

  //sanity check printout
  printSanityCheck(effv, effErrv, "eff");

  cout << endl;
*/  
  gBenchmark->Show("plotDYEfficiency");
  return retCodeOk;
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
