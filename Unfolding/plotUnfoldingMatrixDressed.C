// Derived on July 22, 2014 to do FSR unfolding to dressed-leptons

#include <TBenchmark.h>
#include "../Include/DYTools.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/EventSelector.hh"
#include "../Include/EventWeight.hh"
#include "../Include/UnfoldingMatrix.h"



//=== MAIN MACRO =================================================================================================

int plotUnfoldingMatrixDressed(int analysisIs2D,
			const TString conf,
			DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN,
			DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST,
			double dR_thr=0.1,
			TString rndStudyStr=""
			) {

  const double FSRmassDiff=1.; // largest energy of FSR photon to consider
  std::cout << "dR_thr=" << dR_thr << "\n";


  // check whether it is a calculation
  if (conf.Contains("_DebugRun_")) {
    std::cout << "plotUnfoldingMatrix: _DebugRun_ detected. Terminating the script\n";
    return retCodeOk;
  }

  // normal calculation
  gBenchmark->Start("makeUnfoldingMatrix");

   {
    DYTools::printExecMode(runMode,systMode);
    const int debug_print=1;
    if (!DYTools::checkSystMode(systMode,debug_print,10,
				DYTools::NO_SYST, DYTools::SYST_RND,
				DYTools::FSR_STUDY,
				DYTools::PU_STUDY,
				DYTools::FSR_5plus, DYTools::FSR_5minus,
				DYTools::PILEUP_5plus, DYTools::PILEUP_5minus,
				//DYTools::ESCALE_STUDY,
				DYTools::FSR_RND_STUDY, DYTools::PU_RND_STUDY))
      return retCodeError;
  }

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
  TString extraTag=rndStudyStr;
  EventSelector_t evtSelector(inpMgr,runMode,systMode,
			      extraTag, "", EventSelector::_selectDefault);
  evtSelector.setTriggerActsOnData(false);

  // PU and FSR RND studies have to provide the seed externally
  int globalSeed=-1;
  for (int i=0; (globalSeed<=0) && (i<rndStudyStr.Length()); ++i) {
    globalSeed=atoi(rndStudyStr.Data() + i);
  }

  // Event weight handler
  EventWeight_t evWeight;
  int res=evWeight.init(inpMgr.puReweightFlag(),inpMgr.fewzFlag(),
			systMode,rndStudyStr);
  // May 01, 2014: PU weights have to be applied at all steps
  //EventWeight_t evWeightNoPU; // for FSR unfolding weights
  //if (res) res=evWeightNoPU.init(0,inpMgr.fewzFlag(),systMode,rndStudyStr);
  if (!res) {
    std::cout << "failed to prepare weights\n";
    return retCodeError;
  }

  // Prepare output directory
  inpMgr.constDir(systMode,1);

  //return retCodeOk;

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code
  //==============================================================================================================


  std::cout << mainpart;

  TRandom random;

  std::vector<EventWeight_t*> specEWeightsV;
  std::vector<double> specReweightsV;
  std::vector<EventSelector_t*> evtSelectorV;
  std::vector<TH2D*> specTH2DWeightV; // used for ESCALE_RESIDUAL

  double specWeight=1.;
  int useSpecWeight=0;
  if (systMode==DYTools::FSR_5plus) { specWeight=1.05; useSpecWeight=1; }
  else if (systMode==DYTools::FSR_5minus) { specWeight=0.95; useSpecWeight=1; }
  else if (systMode==DYTools::FSR_RND_STUDY) useSpecWeight=1;

  // check random seed. Special weights use their own,
  // built-in dependencies on seed
  {
    int startSeed=-1;
    if (systMode==DYTools::SYST_RND) {
      std::cout << "setting startSeed=" << globalSeed << "\n";
      startSeed= globalSeed;
    }
    random.SetSeed(startSeed);
    gRandom->SetSeed(startSeed);
  }

  //
  // Set up histograms
  //
  std::vector<TH1D*> hMassv;
  std::vector<TH1D*> hMassBinsv;
  //TH1D *hSelEvents=NULL;

  // debug distributions: 1GeV bins
  //createAnyH1Vec(hMassv,"hMass_",inpMgr.sampleNames(),2500,0.,2500.,"M_{ee} [GeV]","counts/1GeV");
  createAnyH1Vec(hMassv,"hMass_",inpMgr.mcSampleNames(),1490,10.,1500.,"M_{ee} [GeV]","counts/1GeV");
  // debug distributions for current mass bin
  createBaseH1Vec(hMassBinsv,"hMassBins_",inpMgr.mcSampleNames());
  // debug: accumulate info about the selected events in the samples
  //hSelEvents=createAnyTH1D("hSelEvents","hSelEvents",inpMgr.mcSampleCount(),0,inpMgr.mcSampleCount(),"sampleId","event count");


  UnfoldingMatrix_t fsrGood(UnfoldingMatrix::_cFSR, "fsrGood");
  UnfoldingMatrix_t fsrExact(UnfoldingMatrix::_cFSR, "fsrExact");
  UnfoldingMatrix_t fsrDET(UnfoldingMatrix::_cFSR_DET,"fsrDET"); // only relevant indices are checked for ini,fin
  UnfoldingMatrix_t fsrDETexact(UnfoldingMatrix::_cFSR_DET,"fsrDETexact"); // all indices are checked
  // a good working version: response matrix and invResponse are modified after the inversion
  UnfoldingMatrix_t fsrDET_good(UnfoldingMatrix::_cFSR_DET,"fsrDETgood");

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

    for (unsigned int ifile=0; ifile<mcSample->size(); ++ifile) {
      // Read input file
      TFile *infile= new TFile(mcSample->getFName(ifile),"read");
      if (!infile || !infile->IsOpen()) {
	TString skimName=inpMgr.convertSkim2Ntuple(mcSample->getFName(ifile));
	std::cout <<  "  .. failed. Trying <" << skimName << ">" << std::endl;
	infile= new TFile(skimName,"read");
      }
      assert(infile->IsOpen());

      // Get the TTrees
      if (!accessInfo.setTree_withGenPhoton(*infile,"Events")) {
	return retCodeError;
      }

      // Find weight for events for this file
      // The first file in the list comes with weight 1*extraWeightFactor,
      // all subsequent ones are normalized to xsection and luminosity
      ULong_t maxEvents = accessInfo.getEntries();
      // to match old version package (DYee 7TeV paper),
      if (inpMgr.userKeyValueAsInt("USE7TEVMCWEIGHT") &&
	  (isample==0) && (ifile==0)) {
	extraWeightFactor=maxEvents / (inpMgr.totalLumi() * inpMgr.mcSampleInfo(0)->getXsec(ifile));
	//extraWeightFactor=maxEvents / inpMgr.mcSampleInfo(0)->getXsec(ifile);
      }
      //std::cout << "extraWeightFactor=" << extraWeightFactor << ", chk=" << (maxEvents0/inpMgr.mcSampleInfo(0)->getXsec(ifile)) << "\n";
      //const double extraWeightFactor=1.0;
      if (! evWeight.setWeight_and_adjustMaxEvents(maxEvents,
						   inpMgr.totalLumi(),
						   mcSample->getXsec(ifile),
						   extraWeightFactor, inpMgr.selectEventsFlag())) {
	std::cout << "adjustMaxEvents failed\n";
	return retCodeError;
      }
      std::cout << "mcSample xsec=" << mcSample->getXsec(ifile) << ", nEntries=" << maxEvents << "\n";

      std::cout << "       -> sample base weight is " << evWeight.baseWeight() << "\n";
      for (unsigned int iSt=0; iSt<specEWeightsV.size(); ++iSt) {
	specEWeightsV[iSt]->setBaseWeight(evWeight);
      }

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
	if (DYTools::isDebugMode(runMode) &&
	    (ientry>ULong_t(1000000)+DYTools::study2D*ULong_t(2000000))) break; // debug option
	//if (DYTools::isDebugMode(runMode) && (ientry>100)) break; // debug option
	printProgress(250000," ientry=",ientry,maxEvents);
	ec.numEvents_inc();

	// Load generator level info
	accessInfo.GetGen(ientry);
	// If the Z->ll leptons are not electrons, discard this event.
	// This is needed for signal MC samples such as Madgraph Z->ll
	// where all 3 lepton flavors are possible
	if (!accessInfo.genLeptonsAreElectrons()) continue;

	// Load event info
	accessInfo.GetInfoEntry(ientry);

	// Load genPhoton array
	if (accessInfo.GetGenPhotons(ientry)==-1) {
	  std::cout << "error loading photons\n";
	  return retCodeError;
	}

	// Adjust event weight
	// .. here "false" = "not data"
	evWeight.set_PU_and_FEWZ_weights(accessInfo,false);
	//evWeightNoPU.set_PU_and_FEWZ_weights(accessInfo,false);
	if (useSpecWeight) {
	  evWeight.setSpecWeightValue(accessInfo,FSRmassDiff,specWeight);
	  //evWeightNoPU.setSpecWeightValue(accessInfo,FSRmassDiff,specWeight);
	}

	// FSR study correction for weight
	if (systMode==DYTools::FSR_STUDY) {
	  for (unsigned int iSt=0; iSt<specEWeightsV.size(); ++iSt) {
	    specEWeightsV[iSt]->setSpecWeightValue(accessInfo,FSRmassDiff,specReweightsV[iSt]);
	  }
	}

	// setup spec weights
	// .. here "false" = "not data"
	for (unsigned int iSt=0; iSt<specEWeightsV.size(); ++iSt) {
	  specEWeightsV[iSt]->set_PU_and_FEWZ_weights(accessInfo,false);
	}

	if (ientry<20) {
	  std::cout << "ientry=" << ientry << ", "; evWeight.Print(0);
	  //printf("reweight=%4.2lf, fewz_weight=%4.2lf,dE_fsr=%+6.4lf\n",reweight,fewz_weight,(gen->mass-gen->vmass));
	  if (systMode!=DYTools::RESOLUTION_STUDY) {
	    for (unsigned int iSt=0; iSt<specEWeightsV.size(); ++iSt) {
	      std::cout << " specEWeight[" << iSt << "] = "; specEWeightsV[iSt]->Print(0); //std::cout << "\n";
	    }
	  }
	  std::cout << "\n";
	}

	// adjust the scale in the counter to include FEWZ
	// (and possibly PU) weight
	//ec.setScale(evWeight.totalWeight());

	FlatIndex_t fiGenPreFsr, fiGenDressed;
	fiGenPreFsr.setGenPreFsrIdx(accessInfo);

	TLorentzVector dressedE1,dressedE2;
	TLorentzVector *dressedEE = accessInfo.getDressedGenDielectron(dR_thr,
							&dressedE1,&dressedE2);
	if (!dressedEE) {
	  std::cout << "failed to get dressedEE\n";
	  return retCodeError;
	}
	fiGenDressed.setIdx(dressedEE->M(),dressedEE->Rapidity());

	// begin FSR unfolding block
	fsrGood.fillIni(fiGenPreFsr , evWeight.totalWeight());
	fsrGood.fillFin(fiGenDressed, evWeight.totalWeight());
	if (fiGenPreFsr.isValid() && fiGenDressed.isValid()) {
	  fsrGood.fillMigration(fiGenPreFsr, fiGenDressed, evWeight.totalWeight());
	  fsrExact.fillIni(fiGenPreFsr , evWeight.totalWeight());
	  fsrExact.fillFin(fiGenDressed, evWeight.totalWeight());
	  fsrExact.fillMigration(fiGenPreFsr, fiGenDressed, evWeight.totalWeight());
	}

	int preFsrOk=0, dressedOk=0;
	if (evtSelector.inAcceptancePreFsr(accessInfo) &&
	    fiGenPreFsr.isValid()) {
	  preFsrOk=1;
	  fsrDET     .fillIni(fiGenPreFsr, evWeight.totalWeight());
	  fsrDET_good.fillIni(fiGenPreFsr, evWeight.totalWeight());
	}

	if (fiGenDressed.isValid() &&
	    DYTools::goodEtEtaPair(dressedE1.Pt(),dressedE1.Eta(),
				   dressedE2.Pt(),dressedE2.Eta())) {
	  dressedOk=1;
	  fsrDET     .fillFin(fiGenDressed, evWeight.totalWeight());
	  fsrDET_good.fillFin(fiGenDressed, evWeight.totalWeight());
	}

	if (preFsrOk && dressedOk) {
	  fsrDET.fillMigration(fiGenPreFsr, fiGenDressed, evWeight.totalWeight());
	  fsrDET_good.fillMigration(fiGenPreFsr, fiGenDressed, evWeight.totalWeight());
	  fsrDETexact.fillIni(fiGenPreFsr , evWeight.totalWeight());
	  fsrDETexact.fillFin(fiGenDressed, evWeight.totalWeight());
	  fsrDETexact.fillMigration(fiGenPreFsr, fiGenDressed, evWeight.totalWeight());
	}

	delete dressedEE;
	// end of FSR unfolding block

      } // end loop over events
      infile->Close();
      delete infile;
      std::cout << ec << "\n";
      ecTotal.add(ec);
    } // end loop over files
    std::cout << "total counts : " << ecTotal << "\n";

  } // loop over iSample
  } // runMode


  UnfoldingMatrix_t fsrDETcorrections(UnfoldingMatrix::_cFSR_DETcorrFactors,"fsrCorrFactors");

  if (DYTools::processData(runMode)) {
    // Compute the errors on the elements of migration matrix
    fsrGood.finalizeDetMigrationErr();
    fsrExact.finalizeDetMigrationErr();
    fsrDET.finalizeDetMigrationErr();
    fsrDETexact.finalizeDetMigrationErr();
    fsrDET_good.finalizeDetMigrationErr();

  // Find response matrix, which is simply the normalized migration matrix
    std::cout << "find response matrix" << std::endl;
    fsrGood.computeResponseMatrix();
    fsrExact.computeResponseMatrix();
    fsrDET.computeResponseMatrix();
    fsrDETexact.computeResponseMatrix();
    fsrDET_good.computeResponseMatrix();

    std::cout << "find inverted response matrix" << std::endl;
    fsrGood.invertResponseMatrix();
    fsrExact.invertResponseMatrix();
    fsrDET.invertResponseMatrix();
    fsrDETexact.invertResponseMatrix();
    fsrDET_good.invertResponseMatrix();

    fsrDETcorrections.prepareFsrDETcorrFactors(fsrDET,fsrDETexact);
    //fsrDETcorrections.printYields();

    std::cout << "finalize fsrDET_good" << std::endl;
    fsrGood.modifyDETResponseMatrices(fsrExact);
    fsrDET_good.modifyDETResponseMatrices(fsrDETexact);

    std::cout << "prepare flat-index arrays" << std::endl;
    fsrGood.prepareFIArrays();
    fsrExact.prepareFIArrays();
    fsrDET.prepareFIArrays();
    fsrDETexact.prepareFIArrays();
    fsrDET_good.prepareFIArrays();
    fsrDETcorrections.prepareFIArrays();
  }

  //
  // Store constants and reference arrays in files
  //
  if (DYTools::processData(runMode))  std::cout << "store constants in a file" << std::endl;

  //TString outFile=inpMgr.correctionFullFileName("unfolding",systMode,0);
  TString outputDir=inpMgr.constDir(systMode,0);

  //int saveIdxMin=-1;
  TString fnameTag=UnfoldingMatrix_t::generateFNameTag(systMode,globalSeed);
  fnameTag.Prepend("dressed");

  if (DYTools::isDebugMode(runMode)) fnameTag.Prepend("_DebugRun_");
  std::cout << "fnameTag=<" << fnameTag << ">\n";
  CPlot::sOutDir.Append(fnameTag);
  CPlot::sOutDir.ReplaceAll(DYTools::analysisTag,"");
  CPlot::sOutDir.Append(DYTools::analysisTag);

  TString callingMacro="plotUnfoldingMatrix.systMode=";
  callingMacro.Append(SystematicsStudyName(systMode));

  if (DYTools::processData(runMode)) {
    if (//(systMode!=DYTools::NORMAL_RND) &&
	(systMode!=DYTools::RESOLUTION_STUDY) &&
	//(systMode!=DYTools::FSR_STUDY) &&
	(systMode!=DYTools::ESCALE_STUDY)) {
      fsrGood.autoSaveToFile(outputDir,fnameTag,callingMacro);
      fsrExact.autoSaveToFile(outputDir,fnameTag,callingMacro);
      fsrDET.autoSaveToFile(outputDir,fnameTag,callingMacro);
      fsrDETexact.autoSaveToFile(outputDir,fnameTag,callingMacro);
      fsrDET_good.autoSaveToFile(outputDir,fnameTag,callingMacro);
      fsrDETcorrections.autoSaveToFile(outputDir,fnameTag,callingMacro);
    }
  }
  else {
    if (//(systMode!=DYTools::NORMAL_RND) &&
	(systMode!=DYTools::RESOLUTION_STUDY) &&
	//(systMode!=DYTools::FSR_STUDY) &&
	(systMode!=DYTools::ESCALE_STUDY)) {
      if (!fsrGood.autoLoadFromFile(outputDir,fnameTag) ||
	  !fsrExact.autoLoadFromFile(outputDir,fnameTag) ||
	  !fsrDET.autoLoadFromFile(outputDir,fnameTag) ||
	  !fsrDETexact.autoLoadFromFile(outputDir,fnameTag) ||
	  !fsrDET_good.autoLoadFromFile(outputDir,fnameTag) ||
	  !fsrDETcorrections.autoLoadFromFile(outputDir,fnameTag)) {
	std::cout << "loading failed\n";
	return retCodeError;
      }
    }
  }


  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================

  // Plot response and inverted response matrices
  fsrGood.prepareHResponse();
  fsrExact.prepareHResponse();
  fsrDET.prepareHResponse();
  fsrDETexact.prepareHResponse();
  fsrDET_good.prepareHResponse();

  ShowBenchmarkTime("makeUnfoldingMatrixDressed");
  return retCodeOk;
}

