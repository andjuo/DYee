#include <TBenchmark.h>
#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/AccessOrigNtuples.hh"
#include "../Include/EventSelector.hh"
#include "../Include/EventWeight.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/ComparisonPlot.hh"
#include "../Include/MitStyleRemix.hh"
#include "../Include/RangeUser.h"


//=== MAIN MACRO ========================================

int derivePostFsrCS(int analysisIs2D,
		    const TString conf,
		    DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN,
		    DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST,
		    int inAcceptance=0,
		    TString extraTag="",
		    RangeUser_t *rangeUser=NULL) 
{

  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return retCodeError;
  }

  if (DYTools::study2D==0) {
    std::cout << "call this macro for 2D case\n";
    return retCodeError;
  }

  if (inAcceptance==1) { std::cout << "cut on pre-FSR quantities\n"; }
  else if (inAcceptance==2) { std::cout << "cut on post-FSR quantities\n"; }

  gBenchmark->Start("derivePostFsrCS");
  {
    DYTools::printExecMode(runMode,systMode);
    const int debug_print=1;
    if (!DYTools::checkSystMode(systMode,debug_print,2, DYTools::NO_SYST, DYTools::NO_REWEIGHT))
      return retCodeError;
  }

  //-----------------------------------------------------
  // Settings 
  //=====================================================
  
  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(conf)) return retCodeError;
  // no energy correction for this evaluation
  inpMgr.clearEnergyScaleTag();

  // Construct eventSelector, update mgr and plot directory
  EventSelector_t evtSelector(inpMgr,runMode,systMode,
			      extraTag,"",EventSelector::_selectDefault);
  evtSelector.setTriggerActsOnData(false);

  // Acceptance is generator-level quantity and should not 
  // depend on the pile-up. The flag is disabled.
  EventWeight_t evWeight;
  evWeight.init(0*inpMgr.puReweightFlag(),inpMgr.fewzFlag(),systMode);

  // Prepare output directory
  inpMgr.theoryDir(systMode,1);

  TString outFileName=inpMgr.theoryFullFileName("postFsrCS",systMode,1);
  if (inAcceptance==1) outFileName.ReplaceAll(".root","-preFsrDet.root");
  else if (inAcceptance==2) outFileName.ReplaceAll(".root","-postFsrDet.root");
  std::cout << " constructed outFileName=<" << outFileName << ">\n";
  //return retCodeOk;

  //-----------------------------------------------------
  // Main analysis code 
  //=====================================================

  std::cout << mainpart;
  
  //  
  // Set up histograms
  //

  // Accumulate distributions
  std::vector<TH1D*> hMass1Dv;
  std::vector<TH2D*> hMass2Dv;
  std::vector<TH2D*> hMass2DasymV;
  
  // distributions in 1GeV bins
  createAnyH1Vec(hMass1Dv,"hMass1D_",inpMgr.mcSampleNames(),2990,10.,3000.,"#it{M}_{ee} [GeV]","counts/1GeV");
  createBaseH2Vec(hMass2Dv,"hMass2D_",inpMgr.mcSampleNames(),1,1);
  createBaseH2Vec(hMass2DasymV,"hMass2Dasym_",inpMgr.mcSampleNames(),0,1);

  for (unsigned int i=0; i<hMass1Dv.size(); ++i) {
    std::cout << i << "  " << hMass1Dv[i]->GetName() << "\n";
  }
  //return retCodeStop;

  //
  // Access samples and fill histograms
  //
  AccessOrigNtuples_t accessInfo;
  
  //
  // loop over samples
  //
  if (DYTools::processData(runMode)) {

  EventCounterExt_t ecTotal("total");
  double extraWeightFactor=1.;

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
      if (!infile) {
	std::cout <<  "  .. failed" << std::endl;
      }
      assert(infile->IsOpen());
      std::cout << " Reading file <" << mcSample->getFName(ifile) << ">\n";

      // Get the TTrees
      if (!accessInfo.setTree_withGenBranch(*infile,"Events")) {
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
	if (DYTools::isDebugMode(runMode) && (ientry>100000)) break; // debug option
	//if (DYTools::isDebugMode(runMode) && (ientry>100)) break; // debug option
	printProgress(100000," ientry=",ientry,maxEvents);
	ec.numEvents_inc();
	
	// Load generator level info
	accessInfo.GetGen(ientry);
	// If the Z->ll leptons are not electrons, discard this event.
	// This is needed for signal MC samples such as Madgraph Z->ll
	// where all 3 lepton flavors are possible
	if (!accessInfo.genLeptonsAreElectrons()) {
	  std::cout << "genLeptons aren't electrons\n";
	  continue;
	}

	// Load event info to get nPU
	accessInfo.GetInfoEntry(ientry);

	// Adjust event weight
	// .. here "false" = "not data"
	evWeight.set_PU_and_FEWZ_weights(accessInfo,false);

	// adjust the scale in the counter to include FEWZ 
	// (and possibly PU) weight
	//ec.setScale(evWeight.totalWeight());

	// accumulate denominator
	const mithep::TGenInfo *gen= accessInfo.genPtr();
	//if (ientry<10) std::cout << "totWeight=" << evWeight.totalWeight() << "\n";

	if (inAcceptance==1) {
	  if ((gen->vpt_1<10) || (gen->vpt_2<10) ||
	      (fabs(gen->veta_1)>2.5) || (fabs(gen->veta_2)>2.5)) continue;
	}
	if (inAcceptance==2) {
	  if ((gen->pt_1<10) || (gen->pt_2<10) ||
	      (fabs(gen->eta_1)>2.5) || (fabs(gen->eta_2)>2.5)) continue;
	}

	hMass1Dv[isample]->Fill(gen->mass, evWeight.totalWeight());
	hMass2Dv[isample]->Fill(gen->mass, fabs(gen->y), evWeight.totalWeight());
	hMass2DasymV[isample]->Fill(gen->mass, gen->y, evWeight.totalWeight());

	ec.numEventsPassedAcceptance_inc();

      } // loop over events
      ec.print(2);  // print info about file
      ecSample.add(ec); // accumulate event counts
      ecTotal.add(ec);
      
      infile->Close();
      
      delete infile;
    }
    ecSample.print(2); // print info about sample
    //std::cout << "check isample=" << isample << ": " << hMass1Dv[isample]->Integral() << "\n";
    //printHisto(hMass1Dv[isample]);
  }
  ecTotal.print(2);
  } // if (processData)


  std::cout << "outFileName=<" << outFileName << ">\n";
  if (DYTools::processData(runMode)) {
    TFile file(outFileName,"recreate");
    int res=file.IsOpen();
    if (res) res=saveVec(file,hMass1Dv,"mass_1D_1GeV_bins");
    if (res) res=saveVec(file,hMass2Dv,"mass_2D_base");
    if (res) res=saveVec(file,hMass2DasymV,"mass_2D_asym_base");
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
    if (res) res=loadVec(file,hMass1Dv,"mass_1D_1GeV_bins");
    if (res) res=loadVec(file,hMass2Dv,"mass_2D_base");
    if (res) res=loadVec(file,hMass2DasymV,"mass_2D_asym_base");
    file.Close();
    if (!res) {
      std::cout << "error occurred during load from file <" << outFileName << ">\n";
      return retCodeError;
    }
  }

  //for (unsigned int i=0; i<hMass1Dv.size(); ++i) {
  //  std::cout << i << "  " << hMass1Dv[i]->GetName() << "\n";
  //  printHisto(hMass1Dv[i]);
  //}
  //return retCodeStop;


  //-----------------------------------------------------
  // Make plots 
  //=====================================================

  int make_1D_plots=1;
  int make_2D_plots=0;

  if (make_1D_plots) {
    int plot_case=1;
    TString plotTag=TString("1D-") + extraTag;
    TString cpName=TString("cp_") + plotTag;
    TString cpTitle=extraTag;
    ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,cpName,cpTitle,
			"#it{M}_{ee} [GeV]","counts/GeV","ratio");
    cp.SetRefIdx(-111);
    cp.SetYTitleSize(0.06,1.2);
    cp.SetYRange(1e-5,1e7);
    cp.SetLogx(1);
    if (rangeUser) rangeUser->apply(cp);
    //cp.SetLogy(1);
    //cp.SetLogx(0); cp.SetXRange(59.99,120.01);
    //cp.SetYRange(ymin,ymax);

    TString hName=Form("hSumMC_%d",plot_case);
    TH1D *hSum=(TH1D*)hMass1Dv[0]->Clone(hName);
    hSum->Reset();
    hSum->SetDirectory(0);
    hSum->SetTitle(hName);
    
    for (unsigned int ih=0; ih<hMass1Dv.size(); ++ih) {
      hSum->Add(hMass1Dv[ih],1.);
      int color=inpMgr.mcInfo(ih)->colors[0];
      std::cout << "color=" << color << "\n";
      cp.AddToStack(hMass1Dv[ih],inpMgr.mcSampleName(ih),color);
    }
    hSum->SetMarkerStyle(24);
    //cp.AddHist1D(hSum,"all MC","LP",kBlack,1,1,1);

    TString canvName=Form("cx");
    TString canvTitle=Form("cx");
    TCanvas *cx=new TCanvas(canvName,canvTitle,850,700);
    //cp.Prepare2Pads(cx);
    SetSideSpaces(cx,0.05,0.15,0.,0.01);
    cp.Draw(cx,false,"png",0);
    cp.ChangeLegendPos(0.1,0.,0.1,0.);
    cx->Update();
    
    TString figName= TString("fig-") + plotTag;
    SaveCanvas(cx,figName);
  }
     
  if (make_2D_plots) {
    int plot_case=1;
    for (int im=0; im<DYTools::nMassBins; ++im) {
      if (im>0) break;
      TString mStr=Form("M_%2.0lf_%2.0lf-",DYTools::massBinLimits[im],DYTools::massBinLimits[im+1]);
      TString plotTag=TString("2D-") + mStr + extraTag;
      TString cpName=TString("cp_") + plotTag;
      TString cpTitle=extraTag + TString("  ") + mStr;
      ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,cpName,cpTitle,
			  "|y|","counts","ratio");
      cp.SetRefIdx(-111);
      cp.SetYTitleSize(0.06,1.2);
      //cp.SetYRange(1e-5,1e7);
      //cp.SetLogx(1);
      if (rangeUser) rangeUser->apply(cp);
      //cp.SetLogy(1);
      //cp.SetLogx(0); cp.SetXRange(59.99,120.01);
      //cp.SetYRange(ymin,ymax);
      
      std::vector<TH1D*> h1V;
      for (unsigned int ih=0; ih<hMass2Dv.size(); ++ih) {
	//int set_nYBins=DYTools::nYBins[im];
	int set_nYBins=DYTools::nYBinsMax;
	TString hName=Form("h%d_im%d",ih,im);
	TH1D *h= createProfileY(hMass2Dv[ih],im+1,hName,1,NULL, set_nYBins,0.,DYTools::yRangeMax+1e-4);
	h1V.push_back(h);
	if (1 && (ih==0) && (im==1)) {
	  std::cout << mStr;
	  std::cout << " checking :";
	  printHisto(hMass2Dv[ih]);
	  printHisto(h);
	}
      }

      TString hName=Form("hSumMC_%d",plot_case);
      TH1D *hSum=(TH1D*)h1V[0]->Clone(hName);
      hSum->Reset();
      hSum->SetDirectory(0);
      hSum->SetTitle(hName);
    
      std::vector<unsigned int> index;
      index.reserve(h1V.size());
      if (h1V.size()==4) {
	index.push_back(3); index.push_back(2); index.push_back(1);
	index.push_back(0);
      }
      else {
	for (unsigned int ih=0; ih<h1V.size(); ++ih) index.push_back(ih);
      }

      for (unsigned int idx=0; idx<h1V.size(); ++idx) {
	unsigned int ih=index[idx];
	hSum->Add(h1V[ih],1.);
	int color=inpMgr.mcInfo(ih)->colors[0];
	std::cout << "color=" << color << "\n";
	cp.AddToStack(h1V[ih],inpMgr.mcSampleName(ih),color);
	if (DYTools::massBinningSet == DYTools::_MassBins_bins100GeV) {
	  TH1D *hc=(TH1D*)h1V[ih]->Clone(Form("tmp_%d_%d",im,ih));
	  hc->SetMarkerStyle(24);
	  cp.AddHist1D(hc,inpMgr.mcSampleName(ih),"LPE1",color+2,1,1,-1);
	}
      }
      hSum->SetMarkerStyle(24);

      TString canvName=Form("cx%d",im);
      TString canvTitle=Form("cx%d",im);
      TCanvas *cx=new TCanvas(canvName,canvTitle,900,700);
      //cp.Prepare2Pads(cx);
      SetSideSpaces(cx,0.0,0.25,0.,0.02);
      cp.Draw(cx,false,"png",0);
      cp.ChangeLegendPos(0.05,0.,0.05,-0.05);
      if (DYTools::massBinningSet == DYTools::_MassBins_bins100GeV) {
	cp.ChangeLegendPos(0.,0.,0.,-0.15);
      }
      if (h1V.size()!=4) cp.TransLegend(0.0,0.3);
      cx->Update();
    
      TString figName="fig-" + plotTag;
      SaveCanvas(cx,figName);
    }
  }
     
  //-----------------------------------------------------
  // Summary print out
  //=====================================================
  std::cout << std::endl;

  //gBenchmark->Show("derivePostFsrCS");
  ShowBenchmarkTime("derivePostFsrCS");
  return retCodeOk;
}



//=== FUNCTION DEFINITIONS ==============================

//-------------------------------------------------------


  

