#include <TBenchmark.h>
#include "../Include/DYTools.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/EventSelector.hh"
#include "../Include/EventWeight.hh"
#include "../Include/UnfoldingMatrix.h"



//=== MAIN MACRO =================================================================================================

int plotUnfoldingMatrixRho(int analysisIs2D,
			const TString conf,
			DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN,
			DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST,
			TString rndStudyStr=""
			) {

  const double FSRmassDiff=1.; // largest energy of FSR photon to consider


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
    if (!DYTools::checkSystMode(systMode,debug_print,12,
				DYTools::NO_SYST, DYTools::SYST_RND,
				DYTools::RESOLUTION_STUDY, DYTools::FSR_STUDY,
				DYTools::PU_STUDY,
				DYTools::FSR_5plus, DYTools::FSR_5minus,
				DYTools::PILEUP_5plus, DYTools::PILEUP_5minus,
				//DYTools::ESCALE_STUDY,
				DYTools::ESCALE_RESIDUAL,
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
  InputFileMgr_t *yieldInpMgr=NULL; // needed for ESCALE_RESIDUAL
  if (!inpMgr.Load(conf)) return retCodeError;

  // plotDetResponse uses escale!
  if (systMode==DYTools::ESCALE_RESIDUAL) {
    yieldInpMgr= new InputFileMgr_t(inpMgr);
    // create a temporary object to set proper directories
    EventSelector_t tmpEventSelector(*yieldInpMgr,runMode,
		   DYTools::APPLY_ESCALE,"","",EventSelector::_selectDefault);
  }
  else if (systMode!=DYTools::RESOLUTION_STUDY) {
    // no energy correction for this evaluation
    inpMgr.clearEnergyScaleTag();
  }
  else {
    if (inpMgr.energyScaleTag() == "UNCORRECTED") {
      std::cout << "RESOLUTION_STUDY needs energy scale correction\n";
      return retCodeError;
    }
  }

  // Construct eventSelector, update mgr and plot directory
  TString extraTag=TString("_rho_") + rndStudyStr;
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

  if (inpMgr.userKeyValueExists("SpecFile_EffScaleFactor")) {
    inpMgr.addUserKey(TString("SpecFile_EffScaleFactor"),"../../Results-DYee/root_files_reg/constants/DY_j22_19712pb_egamma_Unregressed_energy/covRhoFileSF_nMB41_asymHLT_Unregressed_energy-allSyst_100_v2.root");
  }
  TString rhoCorrFName=inpMgr.correctionFullFileName("scale_factors_asymHLT",systMode,0);
  // permitted setting of the special file
  if (inpMgr.correctionSpecFileName("SpecFile_EffScaleFactor",
					      rhoCorrFName)) {
    std::cout << "inpMgr contained spec file definition <" <<
      rhoCorrFName << ">\n";
  }
  TH2D *hRho=LoadMatrixFields(rhoCorrFName,1,
			      "scaleFactor","scaleFactorErr",1);
  if (!hRho) return retCodeError;


  // Prepare output directory
  inpMgr.rootFileBaseDir("root_files_reg_Rho/");
  std::cout << "\n\trootFileBaseDir=<" << inpMgr.rootFileBaseDir() << ">\n";
  inpMgr.constDir(systMode,1);


  int seedMin=inpMgr.userKeyValueAsInt("SEEDMIN");
  int seedMax=inpMgr.userKeyValueAsInt("SEEDMAX");
  int dSeed=1;
  int seedDiff=(systMode==DYTools::FSR_STUDY) ? 3 : (seedMax-seedMin+1);

  //std::cout << "seedMin..seedMax=" << seedMin << ".." << seedMax << "\n";

  //return retCodeOk;

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code
  //==============================================================================================================


  std::cout << mainpart;

  TRandom random;
  std::vector<ElectronEnergyScale*> escaleV;
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

  // The random seeds are needed only if we are running this script in systematics mode

  if (systMode==DYTools::FSR_STUDY) {
    specReweightsV.reserve(seedDiff);
    specEWeightsV.reserve(seedDiff);
    for (int i=0; i<seedDiff; ++i) {
      double specW= 1 + 0.05*(i-1);
      specReweightsV.push_back(specW);
      specEWeightsV.push_back(new EventWeight_t(evWeight));
    }
  }
  else if (systMode==DYTools::PU_STUDY) {
    if (inpMgr.puReweightFlag()==0) {
      std::cout << "systMode=PU_STUDY needs puReweightFlag=1 in the input file\n";
      return retCodeError;
    }
    specEWeightsV.reserve(2);
    for (int i=0; i<2; ++i) {
      DYTools::TSystematicsStudy_t study=(i==0) ? DYTools::PILEUP_5minus : DYTools::PILEUP_5plus;
      EventWeight_t *ew=new EventWeight_t();
      if (!ew->init(inpMgr.puReweightFlag(),inpMgr.fewzFlag(),study,rndStudyStr)) {
	std::cout << "in plotUnfoldingMatrix.C\n";
	return retCodeError;
      }
      specEWeightsV.push_back(ew);
    }
  }
  else if (systMode==DYTools::SYST_RND) {
    // nothing special about weights
  }
  else if (systMode==DYTools::RESOLUTION_STUDY) {
    if (seedMax==-111) {
      seedMin=-111;
      seedMax= 111;
      dSeed=seedMax-seedMin;
      seedDiff=2;
    }
    if (seedMax < seedMin) {
      printf("error: randomSeedMax=%d, seedMin=%d\n",seedMax,seedMin);
      return retCodeError;
    }
    specEWeightsV.reserve(seedDiff); // not used, but needed as a check
    escaleV.reserve(seedDiff);
    for (int i=seedMin; i<=seedMax; i+=dSeed) {
      TString escaleTag=inpMgr.energyScaleTag() +
	TString(Form("_MIRROR_RANDOMIZED%d",i));
      ElectronEnergyScale *ees= new ElectronEnergyScale(escaleTag);
      if (1) {
	std::cout << "randomSeed=" << i << ". EScale=";
	ees->print();
	std::cout<<"\n";
      }
      escaleV.push_back(ees);
      specEWeightsV.push_back(new EventWeight_t(evWeight));
      EventSelector_t *evtSel=new EventSelector_t(evtSelector,ees);
      // correction acts like on data!
      evtSel->setEScaleCorrectionType(DYTools::DATA,DYTools::ESCALE_STUDY_RND);
      //evtSel->editECName().Append(Form("_idx%d",i+seedMin));
      evtSelectorV.push_back(evtSel);
    }
  }

  // prepare tools for ESCALE_RESIDUAL
  TH2D* h2ShapeWeights=NULL;
  if (systMode==DYTools::ESCALE_RESIDUAL) {
    if (!yieldInpMgr) {
      std::cout << "yieldInpMgr had to be created\n";
      return retCodeError;
    }
    DYTools::TSystematicsStudy_t yieldSystMode=DYTools::APPLY_ESCALE;
    TString shapeFName=yieldInpMgr->signalYieldFullFileName(yieldSystMode,1);
    delete yieldInpMgr; // no longer needed
    if (rndStudyStr.Length()) {
      shapeFName.ReplaceAll(TString("__") + rndStudyStr,"");
    }
    TString subdir="ShapeReweight";
    TString field="zeeMCShapeReweight_";
    TString ddBkg=(inpMgr.userKeyValueAsInt("DDBKG")==1) ? "ddBkg" : "mcBkg";
    field.Append(ddBkg);
    std::cout << "Obtaining shape weights from <"
	      << shapeFName << ">"
	      << "(use" << ddBkg << ")\n";
    h2ShapeWeights=LoadHisto2D(field,shapeFName,subdir,1);
    if (!h2ShapeWeights) {
      std::cout << "failed to find histo \"ZeeMCShapeReweight\"\n";
      return retCodeError;
    }
    std::cout << "shapeWeights:\n"; printHisto(h2ShapeWeights);
    std::cout << "adjusting weights by rho\n";
    if (!multiplyHisto(h2ShapeWeights,hRho,0)) return retCodeError;
    std::cout << "adjusted weights:\n"; printHisto(h2ShapeWeights);

    int ensembleSize= inpMgr.userKeyValueAsInt("RESIDUAL_STUDY_SIZE");
    if (ensembleSize<=0) ensembleSize=100;
    ensembleSize++;
    std::cout << "EScale_residual ensemble size=" << ensembleSize
	      << " (one added for non-randomized entry)\n";

    std::vector<TString> tmpLabelV; // local variable for testing
    specTH2DWeightV.reserve(ensembleSize);
    tmpLabelV.reserve(ensembleSize);

    specTH2DWeightV.push_back(Clone(h2ShapeWeights,
				    "h2NonRndShapeW","h2NonRndShapeW"));
    tmpLabelV.push_back("NonRndShape");

    // prepare histo for randomization. Assume 10% error on the deviation
    for (int ibin=1; ibin<=h2ShapeWeights->GetNbinsX(); ++ibin) {
      for (int jbin=1; jbin<=h2ShapeWeights->GetNbinsY(); ++jbin) {
	double dev=h2ShapeWeights->GetBinContent(ibin,jbin);
	h2ShapeWeights->SetBinError(ibin,jbin, 0.1*dev);
      }
    }

    HistoPair2D_t hpRnd("hpRnd",h2ShapeWeights);

    for (int i=1; i<ensembleSize; ++i) {
      TString name=Form("rndShapeWeight_%d",i);
      TH2D* h2Rnd=hpRnd.randomizedWithinErr(0,name);
      specTH2DWeightV.push_back(h2Rnd);
      tmpLabelV.push_back(name);
    }

    if (0) {
      TCanvas *cx= plotProfiles("cx",specTH2DWeightV,tmpLabelV,NULL,1,
				"MC/data shape reweight");
      cx->Update();
      return retCodeStop;
    }
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


  /*
  TH1F *hMassDiff   = new TH1F("hMassDiff","", 100, -30, 30);
  TH1F *hMassDiffBB = new TH1F("hMassDiffBB","", 100, -30, 30);
  TH1F *hMassDiffEB = new TH1F("hMassDiffEB","", 100, -30, 30);
  TH1F *hMassDiffEE = new TH1F("hMassDiffEE","", 100, -30, 30);

  // These histograms will contain (gen-reco) difference
  // for each (mass, Y) bin in a flattened format
  TH2F *hMassDiffV = new TH2F("hMassDiffV","",
			      nUnfoldingBins, -0.5, nUnfoldingBins-0.5,
			      100, -50.0, 50.0);
  TH2F *hYDiffV = new TH2F("hYDiffV","",
			   nUnfoldingBins, -0.5, nUnfoldingBins-0.5,
			   100, -5.0, 5.0);
  */

//   TH1F *hMassDiffV[nUnfoldingBins];
//   for(int i=0; i<nUnfoldingBins; i++){
//     sprintf(hname,"hMassDiffV_%d",i);
//     hMassDiffV[i] = new TH1F(hname,"",100,-50,50);
//   }

  UnfoldingMatrix_t detResponse(UnfoldingMatrix::_cDET_Response,"detResponse");
  UnfoldingMatrix_t detResponseExact(UnfoldingMatrix::_cDET_Response,"detResponseExact");
  UnfoldingMatrix_t detResponseReversed(UnfoldingMatrix::_cDET_Response,"detResponseReversed");
  UnfoldingMatrix_t detResponsePostFsrDet(UnfoldingMatrix::_cDET_Response,"detResponsePostFsrDet");
  UnfoldingMatrix_t detResponsePreFsrDet(UnfoldingMatrix::_cDET_Response,"detResponsePreFsrDet");

  UnfoldingMatrix_t fsrGood(UnfoldingMatrix::_cFSR, "fsrGood");
  UnfoldingMatrix_t fsrExact(UnfoldingMatrix::_cFSR, "fsrExact");
  UnfoldingMatrix_t fsrDET(UnfoldingMatrix::_cFSR_DET,"fsrDET"); // only relevant indices are checked for ini,fin
  UnfoldingMatrix_t fsrDETexact(UnfoldingMatrix::_cFSR_DET,"fsrDETexact"); // all indices are checked
  // a good working version: response matrix and invResponse are modified after the inversion
  UnfoldingMatrix_t fsrDET_good(UnfoldingMatrix::_cFSR_DET,"fsrDETgood");

  std::vector<UnfoldingMatrix_t*> detRespV;

  if (systMode==DYTools::NO_SYST) {}
  else if (systMode==DYTools::SYST_RND) {
    detRespV.reserve(2);
    for (int ir=0; ir<2; ++ir) {
      TString name=Form("detResponse_seed%d_replica%d",globalSeed,ir);
      detRespV.push_back(new UnfoldingMatrix_t(UnfoldingMatrix::_cDET_Response,name));
    }
  }
  else if (systMode==DYTools::RESOLUTION_STUDY) {
    detRespV.reserve(escaleV.size());
    for (int i=seedMin; i<=seedMax; i+=dSeed) {
      TString name=Form("detResponse_seed%d",i);
      detRespV.push_back(new UnfoldingMatrix_t(UnfoldingMatrix::_cDET_Response,name));
    }
  }
  else if (systMode==DYTools::FSR_STUDY) {
    detRespV.reserve(specReweightsV.size());
    for (unsigned int i=0; i<specReweightsV.size(); i++) {
      TString wStr=(i==0) ? Form("0%2.0f",specReweightsV[i]*100.) : Form("%3.0f",specReweightsV[i]*100.);
      TString name=TString("detResponse_") + wStr;
      detRespV.push_back(new UnfoldingMatrix_t(UnfoldingMatrix::_cDET_Response,name));
    }
  }
  else if (systMode==DYTools::PU_STUDY) {
    if (specEWeightsV.size()!=2) { std::cout << "expected specEWeights.size=2\n"; return retCodeError; }
    detRespV.reserve(specEWeightsV.size());
    for (unsigned int i=0; i<specEWeightsV.size(); i++) {
      TString wStr=(i==0) ? "PU5minus" : "PU5plus";
      TString name=TString("detResponse_") + wStr;
      detRespV.push_back(new UnfoldingMatrix_t(UnfoldingMatrix::_cDET_Response,name));
    }
  }
  else if (systMode==DYTools::ESCALE_RESIDUAL) {
    unsigned int count=specTH2DWeightV.size();
    detRespV.reserve(count);
    for (unsigned int i=0; i<count; ++i) {
      TString name=Form("detResponse_%s",niceNumber(i,count).Data());
      if (i==0) name="detResponse_0_nonRnd";
      detRespV.push_back(new UnfoldingMatrix_t(UnfoldingMatrix::_cDET_Response,name));
    }
  }
  /*
  else if (systMode==DYTools::ESCALE_STUDY) {
    detRespV.reserve(escaleV.size());
    for (unsigned int i=0; i<escaleV.size(); ++i) {
      TString name=TString("detResponse_") + escaleV[i]->calibrationSetShortName();
      detRespV.push_back(new UnfoldingMatrix_t(UnfoldingMatrix_t::_cDET_Response,name));
    }
  }
  */
  if (detRespV.size()) {
    std::cout << "names in detRespV:\n";
    for (unsigned int i=0; i<detRespV.size(); ++i) {
      std::cout << "  - " << detRespV[i]->getName() << "\n";
    }
  }

  // check
  if ((systMode==DYTools::RESOLUTION_STUDY) ||
      (systMode==DYTools::FSR_STUDY) ||
      (systMode==DYTools::PU_STUDY)
      //|| (systMode==DYTools::ESCALE_STUDY)
      ) {
    if (//(detRespV.size() != escaleV.size()) ||
	(detRespV.size() != specEWeightsV.size())) {
      std::cout << "error: detRespV.size=" << detRespV.size()
	//<< ", escaleV.size=" << escaleV.size()
	//<< ", specReweightsV.size=" << specReweightsV.size()
		<< ", specEWeightsV.size=" << specEWeightsV.size()
		<< "\n";
      assert(0);
    }
  }

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
      if (!accessInfo.setTree(*infile,"Events",true)) {
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

	FlatIndex_t fiGenPreFsr, fiGenPostFsr;
	fiGenPreFsr.setGenPreFsrIdx(accessInfo);
	fiGenPostFsr.setGenPostFsrIdx(accessInfo);

	// begin FSR unfolding block
	fsrGood.fillIni(fiGenPreFsr , evWeight.totalWeight());
	fsrGood.fillFin(fiGenPostFsr, evWeight.totalWeight());
	if (fiGenPreFsr.isValid() && fiGenPostFsr.isValid()) {
	  fsrGood.fillMigration(fiGenPreFsr, fiGenPostFsr, evWeight.totalWeight());
	  fsrExact.fillIni(fiGenPreFsr , evWeight.totalWeight());
	  fsrExact.fillFin(fiGenPostFsr, evWeight.totalWeight());
	  fsrExact.fillMigration(fiGenPreFsr, fiGenPostFsr, evWeight.totalWeight());
	}

	int preFsrOk=0, postFsrOk=0;
	if (evtSelector.inAcceptancePreFsr(accessInfo) &&
	    fiGenPreFsr.isValid()) {
	  preFsrOk=1;
	  fsrDET     .fillIni(fiGenPreFsr, evWeight.totalWeight());
	  fsrDET_good.fillIni(fiGenPreFsr, evWeight.totalWeight());
	}

	if (evtSelector.inAcceptance(accessInfo) &&
	    fiGenPostFsr.isValid()) {
	  postFsrOk=1;
	  fsrDET     .fillFin(fiGenPostFsr, evWeight.totalWeight());
	  fsrDET_good.fillFin(fiGenPostFsr, evWeight.totalWeight());
	}

	if (preFsrOk && postFsrOk) {
	  fsrDET.fillMigration(fiGenPreFsr, fiGenPostFsr, evWeight.totalWeight());
	  fsrDET_good.fillMigration(fiGenPreFsr, fiGenPostFsr, evWeight.totalWeight());
	  fsrDETexact.fillIni(fiGenPreFsr , evWeight.totalWeight());
	  fsrDETexact.fillFin(fiGenPostFsr, evWeight.totalWeight());
	  fsrDETexact.fillMigration(fiGenPreFsr, fiGenPostFsr, evWeight.totalWeight());
	}
	// end of FSR unfolding block


	// check event trigger
	if (!evtSelector.eventTriggerOk(accessInfo)) {
	  continue; // no trigger accept? Skip to next event...	
	}
	ec.numEventsPassedEvtTrigger_inc();

	// load dielectron array
	accessInfo.GetDielectrons(ientry);

	// loop through dielectrons
	//int pass=0;
	int candCount=0;
	mithep::TDielectron uncorrDielectron;
	for(Int_t i=0; i<accessInfo.dielectronCount(); i++) {
	  mithep::TDielectron *dielectron = accessInfo.editDielectronPtr(i);
	  ec.numDielectrons_inc();

	  // keep unmodified dielectron
	  if (escaleV.size()) uncorrDielectron.restoreEScaleModifiedValues(*dielectron);

	  // escale may modify dielectron! But it should not here
	  if (!evtSelector.testDielectron(dielectron,accessInfo.evtInfoPtr(),&ec)) continue;
	  //pass=1;

          /******** We have a Z candidate! HURRAY! ********/

	  candCount++;
	  ec.numDielectronsPass_inc();
	  if (ec.numDielectronsOkSameSign_inc(dielectron->q_1,dielectron->q_2)) {
	    // same sign event
	  }

	  //
	  // Fill structures for response matrix

	  FlatIndex_t fiReco;
	  fiReco.setRecoIdx(dielectron);

	  // Get the rho factor
	  double rhoVal=hRho->GetBinContent(fiReco.iM()+1,fiReco.iY()+1);
	  if (ientry<20) {
	    std::cout << "dielectron (m,y)=("
		      << dielectron->mass << ", " << dielectron->y
		      << "), rho=" << rhoVal << "\n";
	  }

	  // Fill the matrix of post-FSR generator level invariant mass and rapidity
	  double diWeight=evWeight.totalWeight() * rhoVal;
	  detResponse.fillIni(fiGenPostFsr, diWeight);
	  detResponse.fillFin(fiReco      , diWeight);

	  int bothFIValid=fiGenPostFsr.isValid() && fiReco.isValid();
	  if (bothFIValid) {
	    ec.numDielectronsGoodMass_inc();
	    detResponse.fillMigration(fiGenPostFsr, fiReco, diWeight);

	    detResponseExact.fillIni(fiGenPostFsr, diWeight);
	    detResponseExact.fillFin(fiReco      , diWeight);
	    detResponseExact.fillMigration(fiGenPostFsr, fiReco, diWeight);
	  }

	  if (postFsrOk) {
	    detResponsePostFsrDet.fillIni(fiGenPostFsr, diWeight);
	    detResponsePostFsrDet.fillFin(fiReco      , diWeight);
	    if (bothFIValid) {
	      detResponsePostFsrDet.fillMigration(fiGenPostFsr,fiReco, diWeight);
	    }
	  }
	  if (preFsrOk) {
	    detResponsePreFsrDet.fillIni(fiGenPreFsr, diWeight);
	    detResponsePreFsrDet.fillFin(fiReco     , diWeight);
	    if (preFsrOk && fiReco.isValid()) {
	      detResponsePreFsrDet.fillMigration(fiGenPreFsr, fiReco, diWeight);
	    }
	  }


	  detResponseReversed.fillIni(fiReco,       diWeight);
	  detResponseReversed.fillFin(fiGenPostFsr, diWeight);
	  if (bothFIValid) {
	    detResponseReversed.fillMigration(fiReco,fiGenPostFsr, diWeight);
	  }

	  if (systMode != DYTools::RESOLUTION_STUDY) {
	    switch(systMode) {
	    case DYTools::SYST_RND: {
	      double rnd=gRandom->Gaus(0,1.);
	      if (rnd==double(0.)) rnd=gRandom->Gaus(0,1.);
	      int idx=(rnd<double(0.)) ? 0:1;
	      detRespV[idx]->fillIni( fiGenPostFsr, diWeight );
	      detRespV[idx]->fillFin( fiReco      , diWeight );
	      if (bothFIValid) {
		detRespV[idx]->fillMigration( fiGenPostFsr, fiReco, diWeight);
	      }
	    }
	      break;

	    case DYTools::ESCALE_RESIDUAL:
	      for (unsigned int iSt=0; iSt<detRespV.size(); ++iSt) {
		const TH2D *h2Rnd= specTH2DWeightV[iSt];
		double w=1.;
		if (fiReco.isValid()) {
		  w=h2Rnd->GetBinContent(fiReco.iM()+1,fiReco.iY()+1);
		  if ((iSt==0) && (ientry<20)) {
		    std::cout << "dielectron(M,Y)=" << dielectron->mass
			      << "," << dielectron->y << ", fiReco="
			      << fiReco << ", specWeight=" << w << "\n";
		  }
		}
		double studyWeight= diWeight * w;
		detRespV[iSt]->fillIni( fiGenPostFsr, studyWeight );
		detRespV[iSt]->fillFin( fiReco      , studyWeight );
		if (bothFIValid) {
		  detRespV[iSt]->fillMigration( fiGenPostFsr, fiReco,
						studyWeight);
		}
	      }
	      break;

	    default:
	      for (unsigned int iSt=0; iSt<detRespV.size(); ++iSt) {
		double studyWeight=specEWeightsV[iSt]->totalWeight();
		detRespV[iSt]->fillIni( fiGenPostFsr, studyWeight );
		detRespV[iSt]->fillFin( fiReco      , studyWeight );
		if (bothFIValid) {
		  detRespV[iSt]->fillMigration( fiGenPostFsr, fiReco,
						studyWeight );
		}
	      }
	    }
	  }

	  if (escaleV.size() && (systMode==DYTools::RESOLUTION_STUDY)) {
	    for (unsigned int iESc=0; iESc<escaleV.size(); ++iESc) {
	      dielectron->restoreEScaleModifiedValues(uncorrDielectron);
	      if (evtSelectorV[iESc]->testDielectron(dielectron,
						accessInfo.evtInfoPtr())) {
		FlatIndex_t fiRecoMdf;
		fiRecoMdf.setRecoIdx(dielectron);
		detRespV[iESc]->fillIni(fiGenPostFsr, diWeight);
		detRespV[iESc]->fillFin(fiRecoMdf   , diWeight);
		if (fiGenPostFsr.isValid() && fiRecoMdf.isValid()) {
		  detRespV[iESc]->fillMigration(fiGenPostFsr,fiRecoMdf,
						diWeight);
		}
	      }
	    }
	  }

	  /*
	  Bool_t isB1 = DYTools::isBarrel(dielectron->scEta_1);
	  Bool_t isB2 = DYTools::isBarrel(dielectron->scEta_2);

	  hMassDiff->Fill(massResmeared - gen->mass);
	  if( isB1 && isB2 )
	    hMassDiffBB->Fill(massResmeared - gen->mass);
	  if( (isB1 && !isB2) || (!isB1 && isB2) )
	    hMassDiffEB->Fill(massResmeared - gen->mass);
	  if( !isB1 && !isB2 )
	    hMassDiffEE->Fill(massResmeared - gen->mass);

	  hMassDiffV->Fill(iIndexFlatGen, massResmeared - gen->mass);
	  hYDiffV   ->Fill(iIndexFlatGen, dielectron->y - gen->y);
	  // 	if(iIndexFlatGen != -1){
	  // 	  hMassDiffV[iIndexFlatGen]->Fill(massResmeared - gen->mass);
	  // 	}
	  */

	} // end loop over dielectrons
	if (candCount>1) ec.numMultiDielectronsOk_inc();

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
    detResponse.finalizeDetMigrationErr();
    detResponseExact.finalizeDetMigrationErr();
    detResponseReversed.finalizeDetMigrationErr();
    detResponsePostFsrDet.finalizeDetMigrationErr();
    detResponsePreFsrDet.finalizeDetMigrationErr();
    fsrGood.finalizeDetMigrationErr();
    fsrExact.finalizeDetMigrationErr();
    fsrDET.finalizeDetMigrationErr();
    fsrDETexact.finalizeDetMigrationErr();
    fsrDET_good.finalizeDetMigrationErr();
    for (unsigned int i=0; i<detRespV.size(); ++i) detRespV[i]->finalizeDetMigrationErr();

  // Find response matrix, which is simply the normalized migration matrix
    std::cout << "find response matrix" << std::endl;
    detResponse.computeResponseMatrix();
    detResponseExact.computeResponseMatrix();
    detResponseReversed.computeResponseMatrix();
    detResponsePostFsrDet.computeResponseMatrix();
    detResponsePreFsrDet.computeResponseMatrix();
    fsrGood.computeResponseMatrix();
    fsrExact.computeResponseMatrix();
    fsrDET.computeResponseMatrix();
    fsrDETexact.computeResponseMatrix();
    fsrDET_good.computeResponseMatrix();
    for (unsigned int i=0; i<detRespV.size(); ++i) detRespV[i]->computeResponseMatrix();

    std::cout << "find inverted response matrix" << std::endl;
    detResponse.invertResponseMatrix();
    detResponseExact.invertResponseMatrix();
    detResponseReversed.invertResponseMatrix();
    detResponsePostFsrDet.invertResponseMatrix();
    detResponsePreFsrDet.invertResponseMatrix();
    fsrGood.invertResponseMatrix();
    fsrExact.invertResponseMatrix();
    fsrDET.invertResponseMatrix();
    fsrDETexact.invertResponseMatrix();
    fsrDET_good.invertResponseMatrix();
    for (unsigned int i=0; i<detRespV.size(); ++i) detRespV[i]->invertResponseMatrix();

    fsrDETcorrections.prepareFsrDETcorrFactors(fsrDET,fsrDETexact);
    //fsrDETcorrections.printYields();

    std::cout << "finalize fsrDET_good" << std::endl;
    fsrGood.modifyDETResponseMatrices(fsrExact);
    fsrDET_good.modifyDETResponseMatrices(fsrDETexact);

    std::cout << "prepare flat-index arrays" << std::endl;
    detResponse.prepareFIArrays();
    detResponseExact.prepareFIArrays();
    detResponseReversed.prepareFIArrays();
    detResponsePostFsrDet.prepareFIArrays();
    detResponsePreFsrDet.prepareFIArrays();
    fsrGood.prepareFIArrays();
    fsrExact.prepareFIArrays();
    fsrDET.prepareFIArrays();
    fsrDETexact.prepareFIArrays();
    fsrDET_good.prepareFIArrays();
    fsrDETcorrections.prepareFIArrays();
    for (unsigned int i=0; i<detRespV.size(); ++i) detRespV[i]->prepareFIArrays();
  }

  //
  // Store constants and reference arrays in files
  //
  if (DYTools::processData(runMode))  std::cout << "store constants in a file" << std::endl;

  //TString outFile=inpMgr.correctionFullFileName("unfolding",systMode,0);
  TString outputDir=inpMgr.constDir(systMode,0);

  //int saveIdxMin=-1;
  TString fnameTag=UnfoldingMatrix_t::generateFNameTag(systMode,globalSeed);
  /*
  {
    TString u="_";
    switch(systMode) {
    case DYTools::NO_SYST:
      fnameTag=DYTools::analysisTag;
      break;
    case DYTools::SYST_RND:
      fnameTag=TString("_replica_") + DYTools::analysisTag;
      //saveIdxMin=0;
     //fnameTag+=seed;
       break;
    case DYTools::RESOLUTION_STUDY:
      fnameTag=TString("_seed_") + DYTools::analysisTag;
      //fnameTag+=seed;
      break;
    case DYTools::FSR_STUDY:
      fnameTag=TString("_fsrStudy_") + DYTools::analysisTag;
      //fnameTag=TString("_reweight_") + DYTools::analysisTag;
      //fnameTag+= int(100*reweightFsr);
      break;
    case DYTools::PU_STUDY:
      fnameTag=TString("_puStudy_") + DYTools::analysisTag;
      break;
    case DYTools::ESCALE_STUDY:
      fnameTag=DYTools::analysisTag+TString("_escale") + u;
      break;
    case DYTools::ESCALE_RESIDUAL:
      fnameTag=DYTools::analysisTag+TString("_escaleResidual");
      break;
    default:
      std::cout<<"requested mode not recognized when determining fnameTag"<<std::endl;
      assert(0);
    }
  }
  */
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
      detResponse.autoSaveToFile(outputDir,fnameTag,callingMacro);  // detResponse, reference mc arrays
      detResponseExact.autoSaveToFile(outputDir,fnameTag,callingMacro);
      detResponseReversed.autoSaveToFile(outputDir,fnameTag,callingMacro);
      detResponsePostFsrDet.autoSaveToFile(outputDir,fnameTag,callingMacro);
      detResponsePreFsrDet.autoSaveToFile(outputDir,fnameTag,callingMacro);
      fsrGood.autoSaveToFile(outputDir,fnameTag,callingMacro);
      fsrExact.autoSaveToFile(outputDir,fnameTag,callingMacro);
      fsrDET.autoSaveToFile(outputDir,fnameTag,callingMacro);
      fsrDETexact.autoSaveToFile(outputDir,fnameTag,callingMacro);
      fsrDET_good.autoSaveToFile(outputDir,fnameTag,callingMacro);
      fsrDETcorrections.autoSaveToFile(outputDir,fnameTag,callingMacro);
    }
    for (unsigned int i=0; i<detRespV.size(); i++)
      detRespV[i]->autoSaveToFile(outputDir,fnameTag,callingMacro);

    // additional saving for systematics
    if (systMode==DYTools::FSR_STUDY) {
      detRespV[0]->autoSaveToFile(inpMgr.constDir(DYTools::FSR_5minus,0),
	  UnfoldingMatrix_t::generateFNameTag(DYTools::FSR_5minus,globalSeed),
				  callingMacro);
      detRespV[2]->autoSaveToFile(inpMgr.constDir(DYTools::FSR_5plus,0),
	  UnfoldingMatrix_t::generateFNameTag(DYTools::FSR_5plus,globalSeed),
				  callingMacro);
    }
    else if (systMode==DYTools::PU_STUDY) {
      TString dir0=inpMgr.constDir(DYTools::PILEUP_5minus,0);
      TString tag0=UnfoldingMatrix_t::generateFNameTag(DYTools::PILEUP_5minus,
						       globalSeed);
      TString dir1=inpMgr.constDir(DYTools::PILEUP_5plus,0);
      TString tag1=UnfoldingMatrix_t::generateFNameTag(DYTools::PILEUP_5plus,
						       globalSeed);
      detRespV[0]->autoSaveToFile(dir0,tag0,callingMacro);
      detRespV[1]->autoSaveToFile(dir1,tag1,callingMacro);
    }
  }
  else {
    if (//(systMode!=DYTools::NORMAL_RND) &&
	(systMode!=DYTools::RESOLUTION_STUDY) &&
	//(systMode!=DYTools::FSR_STUDY) &&
	(systMode!=DYTools::ESCALE_STUDY)) {
      if (!detResponse.autoLoadFromFile(outputDir,fnameTag) ||
	  !detResponseExact.autoLoadFromFile(outputDir,fnameTag) ||
	  !detResponseReversed.autoLoadFromFile(outputDir,fnameTag) ||
	  !detResponsePostFsrDet.autoLoadFromFile(outputDir,fnameTag) ||
	  !detResponsePreFsrDet.autoLoadFromFile(outputDir,fnameTag) ||
	  !fsrGood.autoLoadFromFile(outputDir,fnameTag) ||
	  !fsrExact.autoLoadFromFile(outputDir,fnameTag) ||
	  !fsrDET.autoLoadFromFile(outputDir,fnameTag) ||
	  !fsrDETexact.autoLoadFromFile(outputDir,fnameTag) ||
	  !fsrDET_good.autoLoadFromFile(outputDir,fnameTag) ||
	  !fsrDETcorrections.autoLoadFromFile(outputDir,fnameTag)) {
	std::cout << "loading failed\n";
	return retCodeError;
      }
    }
    for (unsigned int i=0; i<detRespV.size(); i++) detRespV[i]->autoLoadFromFile(outputDir,fnameTag);
  }


  UnfoldingMatrix_t detRespAvg(detResponse.kind, "detResponseAvg");
  if (1 && detRespV.size()) { //computeAverage
    const double weight=1/double(detRespV.size());
    for (unsigned int i=0; i<detRespV.size(); ++i) {
      TString name=Form("tmp_%d",i);
      UnfoldingMatrix_t tmp(*detRespV[i],name);
      tmp.squareDetMigrationErr();
      detRespAvg.addMigration(tmp,weight);
    }
    detRespAvg.finalizeDetMigrationErr();
    detRespAvg.computeResponseMatrix();
    detRespAvg.invertResponseMatrix();
    detRespAvg.prepareFIArrays();
    detRespAvg.autoSaveToFile(outputDir,fnameTag,callingMacro); // detResponse, reference mc arrays
  }


  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================

  /*
  std::cout << "making plots" << std::endl;

  TString unfoldingConstFileName, yieldsFileName;
  detResponse.getFileNames(outputDir,fnameTag, unfoldingConstFileName, yieldsFileName);
  TString unfoldingConstantsPlotFName=unfoldingConstFileName;
  unfoldingConstantsPlotFName.Replace(unfoldingConstantsPlotFName.Index(".root"),
				      sizeof(".root"),
				      "_plots.root");
  TFile *fPlots=new TFile(unfoldingConstantsPlotFName,"recreate");
  if (!fPlots) {
    std::cout << "failed to create a file <" << unfoldingConstantsPlotFName << ">\n";
  }

#ifdef CrossSection_HH
  if (1) {  // study Ini and Fin vecs
    std::vector<VXSectD_t*> dataV;
    TString canvName="canvChk";
    TString fewzTag=(useFewzWeights) ? "_withFEWZ" : "_noFEWZ";
    if (useFewzWeights && regularizeFEWZ) fewzTag="_withMdfFEWZ";
    canvName.Append(fewzTag);
    const int twoPads=0;
    TCanvas *c=new TCanvas(canvName,canvName,600*(1+twoPads),600);
    if (twoPads) {
      c->Divide(2,1);
      c->GetPad(1)->SetLogx(1);
      c->GetPad(2)->SetLogx(1);
      c->GetPad(1)->SetLogy(1);
      c->GetPad(2)->SetLogy(1);
    }
    else { c->SetLogx(1); c->SetLogy(1); }
    int ok=1;

    TH1F *hGenAvg=new TH1F("hGenAvg","hGenAvg",DYTools::nMassBins,DYTools::massBinLimits);
    TH1F *hRecAvg=new TH1F("hRecAvg","hRecAvg",DYTools::nMassBins,DYTools::massBinLimits);
    hGenAvg->SetDirectory(0);
    hRecAvg->SetDirectory(0);
    hGenAvg->Sumw2();
    hRecAvg->Sumw2();

    TH1F *hGenAvg_chk=NULL; //new TH1F("hGenAvg_chk","hGenAvg_chk",DYTools::nMassBins,DYTools::massBinLimits);
    TH1F *hRecAvg_chk=NULL; //new TH1F("hRecAvg_chk","hRecAvg_chk",DYTools::nMassBins,DYTools::massBinLimits);
    if (ok && hGenAvg_chk && hRecAvg_chk) {
      hGenAvg_chk->SetDirectory(0);
      hRecAvg_chk->SetDirectory(0);
      VXSectD_t d("yields_repAvg",DYTools::nUnfoldingBinsMax);
      TString matrixFName,yieldsFName;
      detRespAvg.getFileNames(outputDir,fnameTag, matrixFName,yieldsFName);
      ok=d.Load(matrixFName,"yieldsMcPostFsrGenFIArray","yieldsMcPostFsrRecFIArray","");
      if (ok) ok= (d.FillHisto(hGenAvg_chk,1) && d.FillHistoWithError(hRecAvg_chk));
    }

    for (unsigned int i=0; ok && (i<detRespV.size()); ++i) {
      TString name=Form("yields_rep%d",i);
      VXSectD_t *d=new VXSectD_t(name,DYTools::nUnfoldingBinsMax);
      dataV.push_back(d);
      TString matrixFName,yieldsFName;
      detRespV[i]->getFileNames(outputDir,fnameTag, matrixFName,yieldsFName);
      ok=d->Load(matrixFName,"yieldsMcPostFsrGenFIArray","yieldsMcPostFsrRecFIArray","");
      if (!ok) continue;
      TString hGenName=Form("hGen_rep%d",i);
      TString hRecName=Form("hRec_rep%d",i);
      TH1F *hGen=new TH1F(hGenName,hGenName,DYTools::nMassBins,DYTools::massBinLimits);
      TH1F *hRec=new TH1F(hRecName,hRecName,DYTools::nMassBins,DYTools::massBinLimits);
      hGen->SetDirectory(0);
      hRec->SetDirectory(0);
      if (i==0) {
	hGen->SetTitle(hGen->GetTitle() + fewzTag);
	hRec->SetTitle(hRec->GetTitle() + fewzTag);
      }
      ok= (d->FillHisto(hGen,1) && d->FillHistoWithError(hRec));
      if (!ok) continue;
      hGenAvg->Add(hGen,1/double(detRespV.size()));
      hRecAvg->Add(hRec,1/double(detRespV.size()));
      TString opt="L";
      if (i>0) opt.Append("same");
      int color = i%(50-20) + 20;
      hGen->SetLineColor(color);
      hRec->SetLineColor(color);
      hGen->GetYaxis()->SetRangeUser(5.,3.e6);
      hRec->GetYaxis()->SetRangeUser(5.,3.e6);
      if (twoPads) c->cd(1);
      hGen->Draw(opt);
      if (twoPads) c->cd(2);
      else { if (i==0) opt.Append("same"); }
      hRec->Draw(opt);
    }
    hGenAvg->SetLineColor(kAzure+1);
    hRecAvg->SetLineColor(kRed+1);
    if (twoPads) c->cd(1);
    hGenAvg->Draw("L same");
    if (hGenAvg_chk) hGenAvg_chk->Draw("L same");
    if (twoPads) c->cd(2);
    hRecAvg->Draw("L same");
    if (hRecAvg_chk) hRecAvg_chk->Draw("L same");
    c->Update();
    return;
  }
#endif


  TCanvas *c = MakeCanvas("canvZmass1","canvZmass1",800,600);

  // string buffers
  char ylabel[50];   // y-axis label

  //
  // Draw DY candidate mass at the reconstruction level. Extra
  // smearing is applied. This figure allows one to judge the
  // correctness of the weights aplied to different samples from the
  // smoothness of the combined result.
  //
  sprintf(ylabel,"a.u. / %.1f GeV/c^{2}",hZMassv[0]->GetBinWidth(1));
  CPlot plotZMass1("zmass1","","m(ee) [GeV/c^{2}]",ylabel);
  for(UInt_t i=0; i<fnamev.size(); i++) {
    plotZMass1.AddHist1D(hZMassv[i],labelv[i],"hist",colorv[i],linev[i]);
  }
  plotZMass1.SetLogy();
  plotZMass1.Draw(c);
  SaveCanvas(c,"zmass1");
//   plotZMass1.Draw(c,doSave,format);
//   if (fPlots) { fPlots->cd(); c->Write(); }

  //
  // Draw a plot that illustrates the detector resolution effects.
  // We plot (gen-rec)/gen as a function of mass and rapidity.
  //
  TMatrixD resolutionEffect(DYTools::nMassBins,DYTools::nYBinsMax);
  resolutionEffect = 0;
  for(int i=0; i < resolutionEffect.GetNrows(); i++){
    for(int j=0; j < resolutionEffect.GetNcols(); j++){
      double ngen = (*detResponse.yieldsIni)(i,j);
      double nrec = (*detResponse.yieldsFin)(i,j);
      if( ngen != 0 )
	resolutionEffect(i,j) = (ngen-nrec)/ngen;
    }
  }
  resolutionEffect.Print();
  PlotMatrixVariousBinning(resolutionEffect, "resolution_effect", "LEGO2", NULL);

  //
  // Draw a plot that illustrates the losses due to reconstruction
  // We plot (preFsrExact-preFsr)/preFsrExact as a
  // function of mass and rapidity.
  //
  TMatrixD *unfRecoEffect=detResponseExact.getReconstructionEffect(detResponse);
  unfRecoEffect->Print();
  PlotMatrixVariousBinning(*unfRecoEffect, "reconstruction_effect", "LEGO2", NULL);
  delete unfRecoEffect;

  TMatrixD *unfFsrDETRecoEffect=fsrDETexact.getReconstructionEffect(fsrDET);

  PlotMatrixVariousBinning(*unfFsrDETRecoEffect, "reconstruction_effect_fsrDET", "LEGO2", NULL);
  delete unfFsrDETRecoEffect;

  */

  // Plot response and inverted response matrices
  //std::vector<TH2F*> hResponseV, hInvResponseV;
  //std::vector<TCanvas*> canvV;
  //std::vector<CPlot*> cpResponseV;

  //TH2F *hR, *hIR;
  //TCanvas *e2;
  //CPlot *cpR, *cpIR;
  detResponse.prepareHResponse();
  fsrGood.prepareHResponse();
  fsrExact.prepareHResponse();
  fsrDET.prepareHResponse();
  fsrDETexact.prepareHResponse();
  fsrDET_good.prepareHResponse();

  for (unsigned int iESc=0; iESc<detRespV.size(); ++iESc) {
    detRespV[iESc]->prepareHResponse();
  }

  /*
  // Create a plot of detector resolution without mass binning
  TCanvas *g = MakeCanvas("canvMassDiff","canvMassDiff",600,600);
  CPlot plotMassDiff("massDiff","","reco mass - gen post-FSR mass [GeV/c^{2}]","a.u.");
  hMassDiffBB->Scale(1.0/hMassDiffBB->GetSumOfWeights());
  hMassDiffEB->Scale(1.0/hMassDiffEB->GetSumOfWeights());
  hMassDiffEE->Scale(1.0/hMassDiffEE->GetSumOfWeights());
  plotMassDiff.AddHist1D(hMassDiffBB,"EB-EB","hist",kBlack);
  plotMassDiff.AddHist1D(hMassDiffEB,"EE-EB","hist",kBlue);
  plotMassDiff.AddHist1D(hMassDiffEE,"EE-EE","hist",kRed);
  plotMassDiff.Draw(g);
  SaveCanvas(g,"massDiff");
//   if (fPlots) g->Write();

  // Create a plot of reco - gen post-FSR mass and rapidity difference
  TCanvas *h1 = MakeCanvas("canvMassDiffV","canvMassDiffV",600,600);
  CPlot plotMassDiffV("massDiffV","",
		      "flat index",
		      "reco mass - gen post-FSR mass [GeV/c^{2}]");
  plotMassDiffV.AddHist2D(hMassDiffV,"LEGO");
  plotMassDiffV.Draw(h1);
  SaveCanvas(h1,"hMassDiffV");

  // Create a plot of reco - gen post-FSR mass and rapidity difference
  TCanvas *h2 = MakeCanvas("canvYDiffV","canvYDiffV",600,600);
  CPlot plotYDiffV("massDiffV","",
		      "flat index",
		      "reco Y - gen post-FSR Y");
  plotYDiffV.AddHist2D(hYDiffV,"LEGO");
  plotYDiffV.Draw(h2);
  SaveCanvas(h2,"hYDiffV");

  if (fPlots) {
    fPlots->Close();
    delete fPlots;
    std::cout << "plots saved to a file <" << unfoldingConstantsPlotFName << ">\n";
  }

  //draw errors of Unfolding matrix
  TCanvas *cErrorsResp = MakeCanvas("cErrorsResp","detResponse.DetInvertedResponseErr", 600,600);
  detResponse.DetInvertedResponseErr->Draw("LEGO2");
  cErrorsResp->Update();
  SaveCanvas(cErrorsResp,"cErrorsResp");

  TCanvas *cFsrErrorsResp = MakeCanvas("cErrorsFsr","fsr__.DetInvertedResponseErr", 1200, 600);
  cFsrErrorsResp->Divide(2,1);
  cFsrErrorsResp->cd(1);
  fsrExact.DetInvertedResponseErr->Draw("LEGO2");
  cFsrErrorsResp->cd(2);
  fsrDET.DetInvertedResponseErr->Draw("LEGO2");
  SaveCanvas(cFsrErrorsResp,"cErrorsFsr");



  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;

  detResponse.printConditionNumber();
  fsrExact.printConditionNumber();
  fsrDET.printConditionNumber();
  fsrDETexact.printConditionNumber();

  if (0) {
    //detResponse.printMatrices();
    //fsr.printMatrices();
    //fsrDET.printMatrices();
    fsrDET_Mdf.printMatrices();
    fsrDET_good.printMatrices();
  }

  //Print errors of the Unfolding matrix when they exceed 0.1
  /
  for (int iM=0; iM<DYTools::nMassBins; iM++)
    for (int iY=0; iY<DYTools::nYBins[iM]; iY++)
      for (int jM=0; jM<DYTools::nMassBins; jM++)
        for (int jY=0; jY<DYTools::nYBins[jM]; jY++)
          {
	    int i=DYTools::findIndexFlat(iM,iY);
	    int j=DYTools::findIndexFlat(jM,jY);
             if (DetInvertedResponseErr(i,j)>0.1)
                {
                   std::cout<<"DetInvertedResponseErr("<<i<<","<<j<<")="<<DetInvertedResponseErr(i,j);
                   std::cout<<", DetInvertedResponse("<<i<<","<<j<<")="<<DetInvertedResponse(i,j)<<std::endl;
                   std::cout<<"(iM="<<iM<<", iY="<<iY<<", jM="<<jM<<", jY="<<jY<<")"<<std::endl<<std::endl;
                }
             if (DetInvertedResponseErr2(i,j)>0.1)
                {
                   std::cout<<"DetInvertedResponseErr2("<<i<<","<<j<<")="<<DetInvertedResponseErr2(i,j);
                   std::cout<<", DetInvertedResponse("<<i<<","<<j<<")="<<DetInvertedResponse(i,j)<<std::endl;
                   std::cout<<"(iM="<<iM<<", iY="<<iY<<", jM="<<jM<<", jY="<<jY<<")"<<std::endl<<std::endl;
                }
          }
  /

  /
  if (0) {
    // Printout of all constants, uncomment if needed
    //printf("DetCorrFactor:\n"); DetCorrFactor.Print();
    printf("DetMigration:\n"); DetMigration.Print();
    printf("DetResponse:\n"); DetResponse.Print();

    printf("DetInvertedResponse:\n"); DetInvertedResponse.Print();
    //printf("DetInvertedResponseErr:\n"); DetInvertedResponseErr.Print();
    //printf("DetResponseArr:\n"); DetResponseArr.Print();
    //printf("DetInvertedResponseArr:\n"); DetInvertedResponseArr.Print();
    //printf("DetInvertedResonseErrArr:\n"); DetInvertedResponseErrArr.Print();

    //   printf("Detector corr factor numerator:\n");
    //   DetCorrFactorNumerator.Print();

    printf("yieldsMcPostFsrGen:\n");
    yieldsMcPostFsrGen.Print();

    printf("yieldsMcPostFsrRec:\n");
    yieldsMcPostFsrRec.Print();


    //   printf("Detector corr factor denominator:\n");
    //   DetCorrFactorDenominator.Print();
    //   printf("yieldsMcPostFsrRecArr:\n");
    //   yieldsMcPostFsrRecArr.Print();

    //printf("yieldsMcGen:\n");
    //yieldsMcGen.Print();
  }
  /

  if (0) {
    detResponse.printYields();
    fsrExact.printYields();
    fsrDET.printYields();
  }
*/
  //gBenchmark->Show("makeUnfoldingMatrix");
  ShowBenchmarkTime("makeUnfoldingMatrix");
  return retCodeOk;
}

