#include <TBenchmark.h>
#include "../Include/DYTools.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/EventSelector.hh"
#include "../Include/EventWeight.hh"
#include "../Include/UnfoldingMatrix.h"



//=== MAIN MACRO =================================================================================================

int plotDetResponse(const TString conf,
		    DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN,
		    DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST) {
//systematicsMode 0 (NORMAL) - no systematic calc
//1 (RESOLUTION_STUDY) - systematic due to smearing, 2 (FSR_STUDY) - systematics due to FSR, reweighting
//check mass spectra with reweightFsr = 0.95; 1.00; 1.05  
//mass value until which do reweighting


  // check whether it is a calculation
  if (conf.Contains("_DebugRun_")) {
    std::cout << "plotDetResponse: _DebugRun_ detected. Terminating the script\n";
    return retCodeOk;
  }
 
  // normal calculation
  gBenchmark->Start("makeUnfoldingMatrix");

   {
    DYTools::printExecMode(runMode,systMode);
    const int debug_print=1;
    if (!DYTools::checkSystMode(systMode,debug_print,6, 
				DYTools::NO_SYST, DYTools::SYST_RND,
				DYTools::RESOLUTION_STUDY, DYTools::FSR_STUDY,
				DYTools::ESCALE_STUDY,DYTools::ESCALE_RESIDUAL))
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

 

  const int seedMin=inpMgr.userKeyValueAsInt("SEEDMIN");
  const int seedMax=inpMgr.userKeyValueAsInt("SEEDMAX");
  const int seedDiff=(systMode==DYTools::FSR_STUDY) ? 3 : (seedMax-seedMin+1);

  std::cout << "seedMin..seedMax=" << seedMin << ".." << seedMax << "\n";

  return retCodeOk;

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================

  TRandom random;
  std::vector<ElectronEnergyScale*> escaleV;
  std::vector<double> specReweightsV;
  std::vector<EventSelector_t*> evtSelectorV;

  random.SetSeed(-1);
  gRandom->SetSeed(-1);

  // The random seeds are needed only if we are running this script in systematics mode
  /*
  if ((systMode==DYTools::SYST_RND) || (systMode==DYTools::FSR_STUDY)) {
    specReweightsV.reserve(seedDiff);
    escaleV.reserve(seedDiff);
    evtSelectorV.reserve(seedDiff);
    for (int i=0; i<seedDiff; ++i) {
      double specW= (systMode==DYTools::FSR_STUDY) ?  (1 + 0.05*(i-1)) : 1.0;
      specReweightsV.push_back(specW);
      ElectronEnergyScale *ees= new ElectronEnergyScale(*evtSelector.escale());
      escaleV.push_back(ees);
      EventSelector_t *evtSel=new EventSelector_t(evtSelector,ees);
      evtSel->editECName().Append(Form("_idx%d",i+seedMin));
      evtSelectorV.push_back(evtSel);
    }
  }
  */
  /*
  else if (systematicsMode==DYTools::RESOLUTION_STUDY) {
    if (randomSeedMax < seedMin) {
      printf("error: randomSeedMax=%d, seedMin=%d\n",randomSeedMax,seedMin);
      return;
    }
    specWeightsV.reserve(randomSeedMax-seedMin+1);
    escaleV.reserve(randomSeedMax-seedMin+1);
    for (int i=seedMin; i<randomSeedMax; ++i) {
      ElectronEnergyScale *ees= new ElectronEnergyScale(escale);
      random.SetSeed(i);
      gRandom->SetSeed(i);
      ees->randomizeSmearingWidth(i);
      std::cout << "randomSeed=" << i << ". EScale="; ees->print(); std::cout<<"\n";
      escaleV.push_back(ees);
      specWeightsV.push_back(1.);
    }
  }
  else if (systematicsMode==DYTools::ESCALE_STUDY) {
    EScaleTagFileMgr_t mgr;
    assert(mgr.Load(escaleStudyDefs_FileName));
    specWeightsV.reserve(mgr.size());
    escaleV.reserve(mgr.size());
    for (unsigned int i=0; i<mgr.size(); ++i) {
      ElectronEnergyScale *ees= new ElectronEnergyScale(mgr.escaleTag(i));
      if (!ees ||
	  !ees->isInitialized()) {
	std::cout << "failed to identify escale calibration from tag: >>" 
		  << mgr.escaleTag(i) 
		  << "<<\n";
	assert(0);
      }
      escaleV.push_back(ees);
      specWeightsV.push_back(1.0);
    }
  }
  */

  // prepare tools for ESCALE_RESIDUAL
  /*
  TMatrixD *shapeWeights=NULL;
  if (systematicsMode==DYTools::ESCALE_RESIDUAL) {
    TString shapeFName=TString("../root_files/yields/") + dirTag + 
      TString("/yields_bg-subtracted") + DYTools::analysisTag + TString(".root");
    std::cout << "Obtaining shape weights from <" << shapeFName << ">\n";
    TFile fshape(shapeFName);
    if (!fshape.IsOpen()) {
      std::cout << "failed to open a file <" << shapeFName << ">\n";
      throw 2;
    }
    shapeWeights = (TMatrixD*)fshape.Get("ZeeMCShapeReweight");
    if (!shapeWeights) {
      std::cout << "failed to find object \"ZeeMCShapeReweight\"\n";
      throw 2;
    }
    dirTag += TString("_escale_residual");
    std::cout << "changing dirTag to <" << dirTag << ">\n";
    (*shapeWeights)(0,0)=1; (*shapeWeights)(1,0)=1; (*shapeWeights)(2,0)=1;
    std::cout << "shapeWeights:\n"; shapeWeights->Print(); // return;
  }
  */

  //  
  // Set up histograms
  //
  std::vector<TH1D*> hMassv;
  std::vector<TH1D*> hMassBinsv;
  TH1D *hSelEvents=NULL;

  // debug distributions: 1GeV bins
  //createAnyH1Vec(hMassv,"hMass_",inpMgr.sampleNames(),2500,0.,2500.,"M_{ee} [GeV]","counts/1GeV");
  createAnyH1Vec(hMassv,"hMass_",inpMgr.mcSampleNames(),1490,10.,1500.,"M_{ee} [GeV]","counts/1GeV");
  // debug distributions for current mass bin
  createBaseH1Vec(hMassBinsv,"hMassBins_",inpMgr.mcSampleNames());
  // debug: accumulate info about the selected events in the samples
  hSelEvents=createAnyTH1D("hSelEvents","hSelEvents",inpMgr.mcSampleCount(),0,inpMgr.mcSampleCount(),"sampleId","event count");


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

  UnfoldingMatrix_t detResponse(UnfoldingMatrix_t::_cDET_Response,"detResponse");
  UnfoldingMatrix_t detResponseExact(UnfoldingMatrix_t::_cDET_Response,"detResponseExact");

  UnfoldingMatrix_t fsrGood(UnfoldingMatrix_t::_cFSR, "fsrGood");
  UnfoldingMatrix_t fsrExact(UnfoldingMatrix_t::_cFSR, "fsrExact");
  UnfoldingMatrix_t fsrDET(UnfoldingMatrix_t::_cFSR_DET,"fsrDET"); // only relevant indices are checked for ini,fin
  UnfoldingMatrix_t fsrDETexact(UnfoldingMatrix_t::_cFSR_DET,"fsrDETexact"); // all indices are checked
  // a good working version: response matrix and invResponse are modified after the inversion
  UnfoldingMatrix_t fsrDET_good(UnfoldingMatrix_t::_cFSR_DET,"fsrDET_good"); 

  std::vector<UnfoldingMatrix_t*> detRespV;
  
  /*
  if (systematicsMode==DYTools::SYST_RND) {
    detRespV.reserve(escaleV.size());
    for (unsigned int i=0; i<escaleV.size(); ++i) {
      TString name=Form("detResponse_replica%d",i);
      detRespV.push_back(new UnfoldingMatrix_t(UnfoldingMatrix_t::_cDET_Response,name));
    }
  }
  else if (systematicsMode==DYTools::RESOLUTION_STUDY) {
    detRespV.reserve(escaleV.size());
    for (int i=seedMin; i<randomSeedMax; ++i) {
      TString name=Form("detResponse_seed%d",i);
      detRespV.push_back(new UnfoldingMatrix_t(UnfoldingMatrix_t::_cDET_Response,name));
    }
  }
  else if (systematicsMode==DYTools::FSR_STUDY) {
    detRespV.reserve(escaleV.size());
    for (int i=0; i<3; i++) {
      TString wStr=(i==0) ? Form("0%2.0f",specWeightsV[i]*100.) : Form("%3.0f",specWeightsV[i]*100.);
      TString name=TString("detResponse_") + wStr;
      detRespV.push_back(new UnfoldingMatrix_t(UnfoldingMatrix_t ::_cDET_Response,name));
    }
  }
  else if (systematicsMode==DYTools::ESCALE_STUDY) {
    detRespV.reserve(escaleV.size());
    for (unsigned int i=0; i<escaleV.size(); ++i) {
      TString name=TString("detResponse_") + escaleV[i]->calibrationSetShortName();
      detRespV.push_back(new UnfoldingMatrix_t(UnfoldingMatrix_t::_cDET_Response,name));
    }
  }
  if (detRespV.size()) {
    std::cout << "names in detRespV:\n";
    for (unsigned int i=0; i<detRespV.size(); ++i) {
      std::cout << "  - " << detRespV[i]->getName() << "\n";
    }
  }
  */

  // check
  if ((systMode==DYTools::SYST_RND) ||
      (systMode==DYTools::RESOLUTION_STUDY) ||
      (systMode==DYTools::FSR_STUDY) ||
      (systMode==DYTools::ESCALE_STUDY)) {
    if ((detRespV.size() != escaleV.size()) ||
	(detRespV.size() != specReweightsV.size())) {
      std::cout << "error: detRespV.size=" << detRespV.size() 
		<< ", escaleV.size=" << escaleV.size()
		<< ", specReweightsV.size=" << specReweightsV.size()
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

  } // unfinished!!
  } // runMode

    /*
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
      

    eventCounter_t ec;

    // Read input file
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]); 
    assert(infile);

    // Get the TTrees
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

    // Find weight for events for this file
    // The first file in the list comes with weight 1,
    // all subsequent ones are normalized to xsection and luminosity
    double xsec=xsecv[ifile];
    AdjustXSectionForSkim(infile,xsec,eventTree->GetEntries(),1);
    lumiv[ifile] = eventTree->GetEntries()/xsec;
    double extraScale=1.; // 4839*1666/27166257.; // MC Zee scale in selectEvents
    double scale = extraScale*lumiv[0]/lumiv[ifile];
    cout << "       -> sample weight is " << scale << endl;

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",&info);                TBranch *infoBr       = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Gen",&gen);                  TBranch *genBr = eventTree->GetBranch("Gen");
    eventTree->SetBranchAddress("Dielectron",&dielectronArr); TBranch *dielectronBr = eventTree->GetBranch("Dielectron");
  
    // loop over events    
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      if (debugMode && (ientry>1000000)) break;
      if (ientry%1000000==0) { printProgress("ientry=",ientry,eventTree->GetEntriesFast()); }
      else if (ientry%100000==0) { printProgress("ientry=",ientry,eventTree->GetEntriesFast()); }
      ec.numEvents++;

      genBr->GetEntry(ientry);
      infoBr->GetEntry(ientry);

      double wPU=1.0;
      if (performPUReweight) {
	wPU = puWeight.getWeightHildreth(info->nPU);
      }

      double reweight=1.;
      if (systematicsMode!=DYTools::FSR_STUDY) reweight=1.0;
      else if (((gen->mass)-(gen->vmass))>massLimit) reweight=1.0;
      //else reweight=reweightFsr; // should be taken care by extraWeight

      double fewz_weight = 1.0;
      if (useFewzWeights) fewz_weight=fewz.getWeight(gen->vmass,gen->vpt,gen->vy);
 
      if (ientry<20) {
	printf("reweight=%4.2lf, fewz_weight=%4.2lf,dE_fsr=%+6.4lf\n",reweight,fewz_weight,(gen->mass-gen->vmass));
      }

      int iMassBinGenPreFsr = DYTools::findMassBin(gen->vmass);
      int iYBinGenPreFsr = DYTools::findAbsYBin(iMassBinGenPreFsr, gen->vy);
      int iMassBinGenPostFsr = DYTools::findMassBin(gen->mass);
      int iYBinGenPostFsr = DYTools::findAbsYBin(iMassBinGenPostFsr, gen->y);
      int idxGenPreFsr = DYTools::findIndexFlat(iMassBinGenPreFsr, iYBinGenPreFsr);
      int idxGenPostFsr = DYTools::findIndexFlat(iMassBinGenPostFsr, iYBinGenPostFsr);

      // full fullGenWeight is not affected by reweighting
      double fullGenWeight_tmp = reweight * scale * gen->weight * fewz_weight;
      double fullGenWeightPU = fullGenWeight_tmp * wPU;
      if (ientry<20) std::cout << "fullGenWeightPU= (rew=" << reweight << ")*(scale=" << scale << ")*(gen.w=" << gen->weight << ")*(fewz=" << fewz_weight << ")*(wPU=" << wPU << ") = " << fullGenWeightPU << "\n";

      { // a block for debug purposes
	double fullGenWeight=fullGenWeightPU;

      fsrGood.fillIni(iMassBinGenPreFsr,iYBinGenPreFsr, fullGenWeight);
      fsrGood.fillFin(iMassBinGenPostFsr,iYBinGenPostFsr, fullGenWeight);
      if (validFlatIndices(idxGenPreFsr, idxGenPostFsr)) {
	fsrGood.fillMigration(idxGenPreFsr,idxGenPostFsr, fullGenWeight);
	fsrExact.fillIni(iMassBinGenPreFsr,iYBinGenPreFsr, fullGenWeight);
	fsrExact.fillFin(iMassBinGenPostFsr,iYBinGenPostFsr, fullGenWeight);
	fsrExact.fillMigration(idxGenPreFsr,idxGenPostFsr, fullGenWeight);
      }
 
      int preFsrOk=0, postFsrOk=0;
      if( DYTools::goodEtEtaPair(gen->vpt_1, gen->veta_1,
				 gen->vpt_2, gen->veta_2) ) {
	if (validFlatIndex(idxGenPreFsr)) {
	  preFsrOk=1;
	  fsrDET    .fillIni(iMassBinGenPreFsr,iYBinGenPreFsr,fullGenWeight);
	  fsrDET_Mdf.fillIni(iMassBinGenPreFsr,iYBinGenPreFsr,fullGenWeight);
	  fsrDET_good.fillIni(iMassBinGenPreFsr,iYBinGenPreFsr,fullGenWeight);
	}
      }
      
      if( DYTools::goodEtEtaPair(gen->pt_1, gen->eta_1,
				 gen->pt_2, gen->eta_2 ) ) {
	if (validFlatIndex(idxGenPostFsr)) {
	  postFsrOk=1;
	  fsrDET    .fillFin(iMassBinGenPostFsr,iYBinGenPostFsr,fullGenWeight);
	  fsrDET_Mdf.fillFin(iMassBinGenPostFsr,iYBinGenPostFsr,fullGenWeight);
	  fsrDET_good.fillFin(iMassBinGenPostFsr,iYBinGenPostFsr,fullGenWeight);
	}
      }

      if (preFsrOk && postFsrOk) {
	fsrDET.fillMigration(idxGenPreFsr,idxGenPostFsr, fullGenWeight);
	fsrDET_Mdf.fillMigration(idxGenPreFsr,idxGenPostFsr, fullGenWeight);
	fsrDET_good.fillMigration(idxGenPreFsr,idxGenPostFsr, fullGenWeight);
      	fsrDETexact.fillIni(iMassBinGenPreFsr,iYBinGenPreFsr,fullGenWeight);
	fsrDETexact.fillFin(iMassBinGenPostFsr,iYBinGenPostFsr,fullGenWeight);
	fsrDETexact.fillMigration(idxGenPreFsr,idxGenPostFsr, fullGenWeight);
      }
      } // a block for debug purposes

	
      if( !(requiredTriggers.matchEventTriggerBit(info->triggerBits, 
						  info->runNum))) 
	continue;
      ec.numEventsPassedEvtTrigger++;

      // possible optimization
      // do not consider the event, if reweighting factor is 0.
      //if (wPU==double(0.0)) continue;

     // loop through dielectrons
      dielectronArr->Clear();
      dielectronBr->GetEntry(ientry);    
      for(Int_t i=0; i<dielectronArr->GetEntriesFast(); i++) {
	ec.numDielectronsUnweighted++;
	ec.numDielectrons_inc();

        const mithep::TDielectron *dielectron = (mithep::TDielectron*)((*dielectronArr)[i]);
	
	// Apply selection
	// Eta cuts
	// Asymmetric SC Et cuts
	if (! DYTools::goodEtEtaPair(dielectron->scEt_1, dielectron->scEta_1,
				     dielectron->scEt_2, dielectron->scEta_2) ) {
	  continue;
	}
	ec.numDielectronsGoodEta_inc();
	ec.numDielectronsGoodEt_inc();
   	
	// Both electrons must match trigger objects. At least one ordering
	// must match
	if( ! requiredTriggers.matchTwoTriggerObjectsAnyOrder( dielectron->hltMatchBits_1,
							       dielectron->hltMatchBits_2,
							       info->runNum) ) continue;
	ec.numDielectronsHLTmatched_inc();
	
	// *** Smurf ID is superseeded by new selection ***
// 	// The Smurf electron ID package is the same as used in HWW analysis
// 	// and contains cuts like VBTF WP80 for pt>20, VBTF WP70 for pt<10
// 	// with some customization, plus impact parameter cuts dz and dxy
// 	if(!passSmurf(dielectron)) continue;  

	// The selection below is for the EGM working points from spring 2012
	// recommended for both 2011 and 2012 data
	if(!passEGM2011(dielectron, WP_MEDIUM, info->rhoLowEta)) continue;  
	ec.numDielectronsIDpassed_inc();

        // We have a Z candidate! HURRAY! 

// 	// Apply extra smearing to MC reconstructed dielectron mass
// 	// to better resemble the data
// 	// In systematics mode, use randomized MC smear factors
	double smearingCorrection = (systematicsMode == DYTools::RESOLUTION_STUDY) ?
          escale.generateMCSmearRandomized(dielectron->scEta_1,dielectron->scEta_2) :
          escale.generateMCSmear(dielectron->scEta_1,dielectron->scEta_2);
	double massResmeared = dielectron->mass + smearingCorrection;

	hZMassv[ifile]->Fill(massResmeared,scale * gen->weight * wPU);

	//
	// Fill structures for response matrix and bin by bin corrections
	// Note: there is no handling of overflow, underflow at present,
	// those entries are just dropped. This can be improved.
	// The only possible cases are: underflow in mass and overflow in Y.


	// Fill the matrix of the reconstruction level mass and rapidity
	int iMassReco = DYTools::findMassBin(massResmeared);
	int iYReco = DYTools::findAbsYBin(iMassReco, dielectron->y);

	double shape_weight = 1.0;
	if( shapeWeights && iMassReco != -1 && iYReco != -1) {
	    shape_weight = (*shapeWeights)[iMassReco][iYReco];
	    //std::cout << "massResmeared=" << massResmeared << ", iMassReco=" << iMassReco << ", shapeWeight=" << shape_weight << "\n";
	}
	double fullWeightPU = fullGenWeightPU * shape_weight;

	// Fill the matrix of post-FSR generator level invariant mass and rapidity
	detResponse.fillIni( iMassBinGenPostFsr, iYBinGenPostFsr, fullWeightPU );
	detResponse.fillFin( iMassReco, iYReco, fullWeightPU );

	
        // Unlike the mass vs Y reference yields matrices, to prepare the
	// migration matrix we flatten (mass,Y) into a 1D array, and then
	// store (mass,Y in 1D)_gen vs (mass,Y in 1D)_rec
	int iIndexFlatGen  = DYTools::findIndexFlat(iMassBinGenPostFsr, iYBinGenPostFsr);
 	int iIndexFlatReco = DYTools::findIndexFlat(iMassReco, iYReco);
	if ( validFlatIndices(iIndexFlatGen, iIndexFlatReco) ) {
	  ec.numDielectronsGoodMass_inc();
	  //std::cout << "adding DetMig(" << iIndexFlatGen << "," << iIndexFlatReco << ") = " << reweight << "*" << scale << "*" << gen->weight << "*" << shape_weight << " = "  << (reweight * scale * gen->weight * shape_weight) << "\n";
	  detResponse.fillMigration(iIndexFlatGen, iIndexFlatReco, fullWeightPU );
	  detResponseExact.fillIni( iMassBinGenPostFsr, iYBinGenPostFsr, fullGenWeightPU );
	  detResponseExact.fillFin( iMassReco, iYReco, fullGenWeightPU );
	  detResponseExact.fillMigration(iIndexFlatGen, iIndexFlatReco, fullWeightPU );
	}

	for (unsigned int iESc=0; iESc<escaleV.size(); ++iESc) {
	  double massCorr= (systematicsMode == DYTools::RESOLUTION_STUDY) ?
	    escaleV[iESc]->generateMCSmearRandomized(dielectron->scEta_1,dielectron->scEta_2) :
	    escaleV[iESc]->generateMCSmear(dielectron->scEta_1,dielectron->scEta_2);
	  double mRes= dielectron->mass + massCorr;
	  int iM = DYTools::findMassBin(mRes);
	  int iY = DYTools::findAbsYBin(iM, dielectron->y);
	  detRespV[iESc]->fillIni( iMassBinGenPostFsr, iYBinGenPostFsr, fullWeightPU * specWeightsV[iESc] );
	  detRespV[iESc]->fillFin( iM, iY, fullWeightPU * specWeightsV[iESc] );
	  int idxFReco= DYTools::findIndexFlat(iM, iY);
	  if ( validFlatIndices(iIndexFlatGen, idxFReco) ) {
	    detRespV[iESc]->fillMigration(iIndexFlatGen, idxFReco, fullWeightPU * specWeightsV[iESc] );
	  }
	}

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

      } // end loop over dielectrons

    } // end loop over events 
    delete infile;
    infile=0, eventTree=0;
    std::cout << ec << "\n";
    totEC.add(ec);
  } // end loop over files
  std::cout << "total counts : " << totEC << "\n";
  } 
  delete gen;

  //return;

  if (debugMode==1) return;

  UnfoldingMatrix_t fsrDETcorrections(UnfoldingMatrix_t::_cFSR_DETcorrFactors,"fsrCorrFactors");

  if (debugMode!=-1) {
  // Compute the errors on the elements of migration matrix
  detResponse.finalizeDetMigrationErr();
  detResponseExact.finalizeDetMigrationErr();
  fsrGood.finalizeDetMigrationErr();
  fsrExact.finalizeDetMigrationErr();
  fsrDET.finalizeDetMigrationErr();
  fsrDETexact.finalizeDetMigrationErr();
  fsrDET_Mdf.finalizeDetMigrationErr();
  fsrDET_good.finalizeDetMigrationErr();
  for (unsigned int i=0; i<detRespV.size(); ++i) detRespV[i]->finalizeDetMigrationErr();

  // Find response matrix, which is simply the normalized migration matrix
  std::cout << "find response matrix" << std::endl;
  detResponse.computeResponseMatrix();
  detResponseExact.computeResponseMatrix();
  fsrGood.computeResponseMatrix();
  fsrExact.computeResponseMatrix();
  fsrDET.computeResponseMatrix();
  fsrDETexact.computeResponseMatrix();
  //fsrDET_Mdf.computeResponseMatrix_MdfBeforeNormalization(fsrDETexact);
  fsrDET_Mdf.computeResponseMatrix_Mdf(fsrDETexact);
  fsrDET_good.computeResponseMatrix();
  for (unsigned int i=0; i<detRespV.size(); ++i) detRespV[i]->computeResponseMatrix();

  std::cout << "find inverted response matrix" << std::endl;
  detResponse.invertResponseMatrix();
  detResponseExact.invertResponseMatrix();
  fsrGood.invertResponseMatrix();
  fsrExact.invertResponseMatrix();
  fsrDET.invertResponseMatrix();
  fsrDETexact.invertResponseMatrix();
  fsrDET_Mdf.invertResponseMatrix();
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
  fsrGood.prepareFIArrays();
  fsrExact.prepareFIArrays();
  fsrDET.prepareFIArrays();
  fsrDETexact.prepareFIArrays();
  fsrDET_Mdf.prepareFIArrays();
  fsrDET_good.prepareFIArrays();
  fsrDETcorrections.prepareFIArrays();
  for (unsigned int i=0; i<detRespV.size(); ++i) detRespV[i]->prepareFIArrays();  
  }

  //
  // Store constants and reference arrays in files
  //
  if (debugMode!=-1) std::cout << "store constants in a file" << std::endl;

  TString outputDir(TString("../root_files/constants_tmp/")+dirTag);
  if((systematicsMode==DYTools::NORMAL_RND) || 
     (systematicsMode==DYTools::RESOLUTION_STUDY) || 
     (systematicsMode==DYTools::FSR_STUDY) ||
     (systematicsMode==DYTools::ESCALE_STUDY))
    outputDir = TString("../root_files/systematics_tmp/")+dirTag;
  gSystem->mkdir(outputDir,kTRUE);
  outputDir.Append("/");

  //int saveIdxMin=-1;

  TString fnameTag="";
  {
    TString u="_";
    switch(systematicsMode) {
    case DYTools::NORMAL: 
      fnameTag=DYTools::analysisTag; 
      break;
    case DYTools::NORMAL_RND: 
      fnameTag=TString("_replica_") + DYTools::analysisTag; 
      //saveIdxMin=0;
     //fnameTag+=seed;
       break;
    case DYTools::RESOLUTION_STUDY: 
      fnameTag=TString("_seed_") + DYTools::analysisTag;
      //fnameTag+=seed;
      break;
    case DYTools::FSR_STUDY:
      fnameTag=TString("_reweight_") + DYTools::analysisTag;
      //fnameTag+= int(100*reweightFsr);
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
  //if (!useFewzWeights) { fnameTag=TString("_noFEWZ") + fnameTag; }
  if (useFewzWeights) { 
    TString extra=(regularizeFEWZ) ? "_withMdfFEWZ" : "_withFEWZ";
    fnameTag=extra + fnameTag; 
  }
  std::cout << "fnameTag=<" << fnameTag << ">\n";
  CPlot::sOutDir=TString("plots") + fnameTag;

  if (debugMode!=-1) {
    if (//(systematicsMode!=DYTools::NORMAL_RND) && 
	(systematicsMode!=DYTools::RESOLUTION_STUDY) && 
	(systematicsMode!=DYTools::FSR_STUDY) && 
	(systematicsMode!=DYTools::ESCALE_STUDY)) {
      detResponse.autoSaveToFile(outputDir,fnameTag);  // detResponse, reference mc arrays
      detResponseExact.autoSaveToFile(outputDir,fnameTag);
      fsrGood.autoSaveToFile(outputDir,fnameTag);
      fsrExact.autoSaveToFile(outputDir,fnameTag);
      fsrDET.autoSaveToFile(outputDir,fnameTag);
      fsrDETexact.autoSaveToFile(outputDir,fnameTag);
      fsrDET_Mdf.autoSaveToFile(outputDir,fnameTag);
      fsrDET_good.autoSaveToFile(outputDir,fnameTag);
      fsrDETcorrections.autoSaveToFile(outputDir,fnameTag);
    }
    for (unsigned int i=0; i<detRespV.size(); i++) detRespV[i]->autoSaveToFile(outputDir,fnameTag);
  }
  else {
    if (//(systematicsMode!=DYTools::NORMAL_RND) && 
	(systematicsMode!=DYTools::RESOLUTION_STUDY) && 
	(systematicsMode!=DYTools::FSR_STUDY) &&
	(systematicsMode!=DYTools::ESCALE_STUDY)) {
      if (!detResponse.autoLoadFromFile(outputDir,fnameTag) ||
	  !detResponseExact.autoLoadFromFile(outputDir,fnameTag) ||
	  !fsrGood.autoLoadFromFile(outputDir,fnameTag) ||
	  !fsrExact.autoLoadFromFile(outputDir,fnameTag) ||
	  !fsrDET.autoLoadFromFile(outputDir,fnameTag) ||
	  !fsrDETexact.autoLoadFromFile(outputDir,fnameTag) ||
	  !fsrDET_Mdf.autoLoadFromFile(outputDir,fnameTag) ||
	  !fsrDET_good.autoLoadFromFile(outputDir,fnameTag) ||
	  !fsrDETcorrections.autoLoadFromFile(outputDir,fnameTag)) {
	std::cout << "loading failed\n";
	return;
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
    detRespAvg.autoSaveToFile(outputDir,fnameTag);  // detResponse, reference mc arrays
  }


  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  

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
  gBenchmark->Show("makeUnfoldingMatrix");
  return retCodeOk;
}

