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

int derivePreFsrCS_PU(int analysisIs2D,
		   const TString conf,
		   DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN,
		   DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST,
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

  gBenchmark->Start("derivePreFsrCS");
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
  EventWeight_t evWeightFewz;
  evWeightFewz.init(0,1,systMode);
  EventWeight_t evWeightPUFewz;
  evWeightPUFewz.init(1,1,systMode);
  EventWeight_t evWeightNoPUNoFewz;
  evWeightNoPUNoFewz.init(0,0,systMode);
  EventWeight_t evWeightPUNoFewz;
  evWeightPUNoFewz.init(1,0,systMode);

  // Prepare output directory
  inpMgr.theoryDir(systMode,1);

  TString outFileName=inpMgr.theoryFullFileName("preFsrCS_PU",systMode,1);
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
  std::vector<TH2D*> hMass2DnoPUnoFEWZv;
  std::vector<TH2D*> hMass2DnoPUwFEWZv;
  std::vector<TH2D*> hMass2DwPUnoFEWZv;
  std::vector<TH2D*> hMass2DwPUwFEWZv;
  std::vector<TH2D*> hMass2DasymV;
  PUReweight_t puNoPUNoFewz(PUReweight_t::_none);
  PUReweight_t puNoPUwFewz(PUReweight_t::_none);
  PUReweight_t puWPUNoFewz(PUReweight_t::_none);
  PUReweight_t puWPUwFewz(PUReweight_t::_none);

  TString fname_puNoPUNoFewz=Form("pu_noPU_noFEWZ.root");
  TString fname_puNoPUwFewz=Form("pu_noPU_wFEWZ.root");
  TString fname_puWPUNoFewz=Form("pu_wPU_noFEWZ.root");
  TString fname_puWPUwFewz=Form("pu_wPU_wFEWZ.root");

  int createPUFile=(DYTools::processData(runMode)) ? 1:0;
  if (!puNoPUNoFewz.setFile(fname_puNoPUNoFewz,createPUFile) ||
      !puNoPUwFewz.setFile(fname_puNoPUwFewz,createPUFile) ||
      !puWPUNoFewz.setFile(fname_puWPUNoFewz,createPUFile) ||
      !puWPUwFewz.setFile(fname_puWPUwFewz,createPUFile)) {
    std::cout << "failed to create the PU distribution files\n";
    return retCodeError;
  }
  if (!puNoPUNoFewz.setActiveSample("npv_noPU_noFewz") ||
      !puNoPUwFewz.setActiveSample("npv_noPU_wFewz") ||
      !puWPUNoFewz.setActiveSample("npv_wPU_noFewz") ||
      !puWPUwFewz.setActiveSample("npv_wPU_wFewz")) {
    std::cout << "failed to set active sample\n";
    return retCodeError;
  }

  // distributions in 1GeV bins
  createAnyH1Vec(hMass1Dv,"hMass1D_",inpMgr.mcSampleNames(),2990,10.,3000.,"#it{M}_{ee} [GeV]","counts/1GeV");
  createBaseH2Vec(hMass2Dv,"hMass2D_",inpMgr.mcSampleNames(),1,1);
  createBaseH2Vec(hMass2DnoPUnoFEWZv,"hMass_noPU_noFewz_2D_",inpMgr.mcSampleNames(),1,1);
  createBaseH2Vec(hMass2DnoPUwFEWZv,"hMass_noPU_wFewz_2D_",inpMgr.mcSampleNames(),1,1);
  createBaseH2Vec(hMass2DwPUnoFEWZv,"hMass_wPU_noFewz_2D_",inpMgr.mcSampleNames(),1,1);
  createBaseH2Vec(hMass2DwPUwFEWZv,"hMass_wPU_wFewz_2D_",inpMgr.mcSampleNames(),1,1);
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
      if (! evWeightFewz.setWeight_and_adjustMaxEvents(maxEvents, inpMgr.totalLumi(), mcSample->getXsec(ifile), 
						   extraWeightFactor, inpMgr.selectEventsFlag())) {
	std::cout << "adjustMaxEvents failed\n";
	return retCodeError;
      }
      evWeightPUFewz.setBaseWeight(evWeightFewz);
      evWeightNoPUNoFewz.setBaseWeight(evWeightFewz);
      evWeightPUNoFewz.setBaseWeight(evWeightFewz);
      std::cout << "mcSample xsec=" << mcSample->getXsec(ifile) << ", nEntries=" << maxEvents << "\n";
      
      std::cout << "       -> sample base weight is " << evWeightFewz.baseWeight() << "\n";
    
      // loop through events
      EventCounterExt_t ec(Form("%s_file%d",mcSample->name.Data(),ifile));
      ec.setIgnoreScale(0); // 1 - count events, 0 - take weight in account
      // adjust the scale in the counter
      // if FEWZ weight should be considered, use evWeight.totalWeight() after
      // the FEWZ weight has been identified (see a line below)
      ec.setScale(evWeightFewz.baseWeight());

      std::cout << "numEntries = " << accessInfo.getEntriesFast() 
		<< ", " << maxEvents << " events will be used" << std::endl;


      //if (isample<2) continue;

      for(ULong_t ientry=0; ientry<maxEvents; ientry++) {
	if (DYTools::isDebugMode(runMode) && (ientry>1000)) break; // debug option
	//if (DYTools::isDebugMode(runMode) && (ientry>100)) break; // debug option
	printProgress(250000," ientry=",ientry,maxEvents);
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
	evWeightFewz.set_PU_and_FEWZ_weights(accessInfo,false);
	evWeightPUFewz.set_PU_and_FEWZ_weights(accessInfo,false);
	evWeightNoPUNoFewz.set_PU_and_FEWZ_weights(accessInfo,false);
	evWeightPUNoFewz.set_PU_and_FEWZ_weights(accessInfo,false);

	// adjust the scale in the counter to include FEWZ 
	// (and possibly PU) weight
	//ec.setScale(evWeight.totalWeight());

	// accumulate denominator
	const mithep::TGenInfo *gen= accessInfo.genPtr();

	if (DYTools::isDebugMode(runMode)) {
	  std::cout << "ientry=" << ientry << "\n";
	  std::cout << " nn " << evWeightNoPUNoFewz << "\n";
	  std::cout << " nw " << evWeightFewz << "\n";
	  std::cout << " wn " << evWeightPUNoFewz << "\n";
	  std::cout << " ww " << evWeightPUFewz << "\n";
	}


	int nPVs= accessInfo.getNPV(false);
	puNoPUNoFewz.Fill(nPVs, evWeightNoPUNoFewz.totalWeight());
	puNoPUwFewz.Fill(nPVs, evWeightFewz.totalWeight());
	puWPUNoFewz.Fill(nPVs, evWeightPUNoFewz.totalWeight());
	puWPUwFewz.Fill(nPVs, evWeightPUFewz.totalWeight());

	//if (ientry<10) std::cout << "totWeight=" << evWeight.totalWeight() << "\n";
	hMass1Dv[isample]->Fill(gen->vmass, evWeightNoPUNoFewz.totalWeight());
	hMass2Dv[isample]->Fill(gen->vmass, fabs(gen->vy), evWeightNoPUNoFewz.totalWeight());
	hMass2DasymV[isample]->Fill(gen->vmass, gen->vy, evWeightNoPUNoFewz.totalWeight());

	hMass2DnoPUnoFEWZv[isample]->Fill(gen->vmass, fabs(gen->vy), evWeightNoPUNoFewz.totalWeight());
	hMass2DnoPUwFEWZv[isample]->Fill(gen->vmass, fabs(gen->vy), evWeightFewz.totalWeight());
	hMass2DwPUnoFEWZv[isample]->Fill(gen->vmass, fabs(gen->vy), evWeightPUNoFewz.totalWeight());
	hMass2DwPUwFEWZv[isample]->Fill(gen->vmass, fabs(gen->vy), evWeightPUFewz.totalWeight());

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
    if (res) res=saveVec(file,hMass2DnoPUnoFEWZv,"mass_noPU_noFEWZ_2D_base");
    if (res) res=saveVec(file,hMass2DnoPUwFEWZv,"mass_noPU_wFEWZ_2D_base");
    if (res) res=saveVec(file,hMass2DwPUnoFEWZv,"mass_wPU_noFEWZ_2D_base");
    if (res) res=saveVec(file,hMass2DwPUwFEWZv,"mass_wPU_wFEWZ_2D_base");
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
    if (res) res=loadVec(file,hMass2DnoPUnoFEWZv,"mass_noPU_noFEWZ_2D_base");
    if (res) res=loadVec(file,hMass2DnoPUwFEWZv,"mass_noPU_wFEWZ_2D_base");
    if (res) res=loadVec(file,hMass2DwPUnoFEWZv,"mass_wPU_noFEWZ_2D_base");
    if (res) res=loadVec(file,hMass2DwPUwFEWZv,"mass_wPU_wFEWZ_2D_base");
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

  int make_PU_plots=1;
  int make_1D_plots=0;
  int make_2D_plots=0;

  if (make_PU_plots) {
    //TH1D* h_noPU_noFEWZ= convert_TH1F_to_TH1D(puNoPUNoFewz.getHActive(),"npv_noPU_noFewz");
    //TH1D* h_noPU_wFEWZ = convert_TH1F_to_TH1D(puNoPUwFewz.getHActive() ,"npv_noPU_wFewz");
    //TH1D *h_wPU_noFEWZ = convert_TH1F_to_TH1D(puWPUNoFewz.getHActive() ,"npv_wPU_noFewz");
    //TH1D *h_wPU_wFEWZ  = convert_TH1F_to_TH1D(puWPUwFewz.getHActive()  ,"npv_wPU_wFewz");
    TH1D* h_noPU_noFEWZ= Clone(puNoPUNoFewz.getHActive(),"npv_nPU_noFewz","");
    TH1D* h_noPU_wFEWZ = Clone(puNoPUwFewz.getHActive() ,"npv_noPU_wFewz","");
    TH1D *h_wPU_noFEWZ = Clone(puWPUNoFewz.getHActive() ,"npv_wPU_noFewz","");
    TH1D *h_wPU_wFEWZ  = Clone(puWPUwFewz.getHActive()  ,"npv_wPU_wFewz" ,"");

    removeError(h_noPU_noFEWZ);
    removeError(h_noPU_wFEWZ);
    removeError(h_wPU_noFEWZ);
    removeError(h_wPU_wFEWZ);

    TString cpName="compPlot";
    TString cpTitle="PU distributions in Zee MC";
    ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,cpName,cpTitle,
			"#PU_{mean}","count","ratio");
    cp.AddHist1D(h_noPU_noFEWZ,"no corr.","LP",TAttMarker(kBlack,24,0.8),1,1,1);
    cp.AddHist1D(h_noPU_wFEWZ, "w/FEWZ","LP",TAttMarker(kBlue,20,0.8),2,1,1);
    cp.AddHist1D(h_wPU_noFEWZ, "w/PU","LP",TAttMarker(kGreen+1,5,0.8),3,1,1);
    cp.AddHist1D(h_wPU_wFEWZ, "w/PU w/FEWZ","LP",TAttMarker(kRed+1,3,0.8),4,1,1);
    
    TCanvas *cx=new TCanvas("cPU","cPU",800,800);
    cp.Prepare2Pads(cx);
    cp.Draw(cx);
    cx->Update();
  }

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

  //gBenchmark->Show("derivePreFsrCS");
  ShowBenchmarkTime("derivePreFsrCS");
  return retCodeOk;
}



//=== FUNCTION DEFINITIONS ==============================

//-------------------------------------------------------


  

