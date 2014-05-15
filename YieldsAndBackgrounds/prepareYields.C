
// This is extra macro Jul 22, 2013


//================================================================================================
//
// Prepare binned histograms with signal and background events for further analysis.
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1D.h>                   // 1D histograms
#include <THStack.h>
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TRandom.h>
#include <TDatime.h>                // time stamp

#include <TLatex.h> 

#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "../Include/CPlot.hh"          // helper class for plots
#include "../Include/MitStyleRemix.hh"  // style settings for drawing
#include "../Include/MyTools.hh"        // miscellaneous helper functions

// define structures to read in ntuple
#include "../Include/ZeeData.hh"

#include "../Include/ElectronEnergyScale.hh"        // energy scale correction
#include "../Include/DYTools.hh"
#include "../Include/ComparisonPlot.hh"
#include "../Include/EventSelector.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/DielectronChargeCounter.hh"

#endif

#ifndef ZeeData_storeGoldenFlag
#error prepareYieldsR9 needs event selected with "ZeeData_storeGoldenFlag" defined
#endif

// Forward declarations
/*
void DrawMassPeak(vector<TH1F*> hMassv, vector<CSample*> samplev, vector<TString> snamev, TH1F* hMassDibosons, bool hasData, 
                   bool mergeDibosons, TString labelDibosons, Int_t colorDibosons, Double_t lumi, char* lumitext,  bool actualBinning);

void DrawFlattened(vector<TMatrixD*> yields, vector<TMatrixD*> yieldsSumw2, vector<CSample*> samplev, vector<TString> snamev, bool hasData, 
                   bool mergeDibosons, TString labelDibosons, Int_t colorDibosons, Double_t lumi, char* lumitext);

void Draw6Canvases(vector<TMatrixD*> yields, vector<TMatrixD*> yieldsSumw2,
                    vector<CSample*> samplev, vector<TString> snamev, 
                    bool hasData, double dataOverMc, double* dataOverMcEachBin, bool normEachBin=1, bool singleCanvas=0);
void SomeHistAttributes (TH1F* hist, TString samplename);
*/
//void SaveCanvas(TCanvas* canv, TString canvName);

//=== MAIN MACRO =================================================================================================

int prepareYields(int analysisIs2D,
		  const TString conf  = "default",
		  DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN,
		  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST,
		  int iSeed=-1)
{  
  gBenchmark->Start("prepareYields");

  {
    using namespace DYTools;
    DYTools::printExecMode(runMode,systMode);
    const int debug_print=1;
    int systModeCount=3;
#ifdef ZeeData_storeUnregEn
    systModeCount+=6;
#endif
    if (!DYTools::checkSystMode(systMode,debug_print,systModeCount,
				DYTools::NO_SYST, //DYTools::ESCALE_STUDY, 
				DYTools::ESCALE_STUDY_RND,
				DYTools::APPLY_ESCALE,
#ifdef ZeeData_storeUnregEn
	        DYTools::UNREGRESSED_ENERGY,
		ESCALE_DIFF_0000, ESCALE_DIFF_0005, ESCALE_DIFF_0010, ESCALE_DIFF_0015, ESCALE_DIFF_0020
#endif
				)) 
      return retCodeError;
  }

  DYTools::TSystematicsStudy_t outputSystMode=systMode;
  if (systMode==DYTools::UNREGRESSED_ENERGY) {
#ifndef ZeeData_storeUnregEn
    std::cout << "systMode=" << SystematicsStudyName(systMode) << " requires ZeeData_storeUnregEn to be defined\n";
    return retCodeErr;
#endif
    systMode=DYTools::NO_SYST;  // input files have no sytematics
  }
 
  int applyEScale=0;
  int applyEtEtaCut=0;
  double mdfEnFactor=0.;
  switch(outputSystMode) {
  case DYTools::ESCALE_DIFF_0000:
  case DYTools::ESCALE_DIFF_0005:
  case DYTools::ESCALE_DIFF_0010:
  case DYTools::ESCALE_DIFF_0015:
  case DYTools::ESCALE_DIFF_0020:
    systMode=DYTools::LOWER_ET_CUT;
    applyEtEtaCut=1;
    break;
  case DYTools::APPLY_ESCALE:
    systMode=DYTools::LOWER_ET_CUT;
    applyEtEtaCut=1;
    applyEScale=1;
    break;
  case DYTools::ESCALE_STUDY_RND:
    systMode=DYTools::LOWER_ET_CUT;
    applyEtEtaCut=1;
    applyEScale=2; // randomized scale factors
    break;
  default: ;
  }
  switch(outputSystMode) {
  case DYTools::ESCALE_DIFF_0000: mdfEnFactor=0.000; break;
  case DYTools::ESCALE_DIFF_0005: mdfEnFactor=0.005; break;
  case DYTools::ESCALE_DIFF_0010: mdfEnFactor=0.010; break;
  case DYTools::ESCALE_DIFF_0015: mdfEnFactor=0.015; break;
  case DYTools::ESCALE_DIFF_0020: mdfEnFactor=0.020; break;
  default: ;
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

  InputFileMgr_t inpMgr_forOutput(inpMgr);

  //mgr.Print();

  // no energy correction for this evaluation
  //inpMgr.clearEnergyScaleTag();

  // Construct eventSelector, update mgr and plot directory
  TString extraTag;
  TString plotExtraTag;
  EventSelector_t evtSelector(inpMgr,runMode,systMode,
		      extraTag, plotExtraTag, EventSelector::_selectDefault);

  if (outputSystMode==DYTools::ESCALE_STUDY_RND) {
    inpMgr_forOutput.editEnergyScaleTag().Append(Form("_RANDOMIZED%d",iSeed));
  }
  EventSelector_t evtSelector_forOutput(inpMgr_forOutput,runMode,
					outputSystMode,
		      extraTag, plotExtraTag, EventSelector::_selectDefault);
  // However, the plots should be saved according to the outputSystMode
  // note: this is taken case by evtSelector_forOutput
  //evtSelector.SetPlotOutDir(runMode,outputSystMode,plotExtraTag,1);

  int createDir=DYTools::processData(runMode);
  TString yieldFullName= inpMgr_forOutput.yieldFullFileName(-1,outputSystMode,createDir);
  std::cout << "will (work with)/(produce) yieldFile=<" 
	    << yieldFullName << ">\n";

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //============================================================================================================== 



  vector<TH2D*> yieldsBaseH2;
  vector<TH2D*> yields; // non-uniform rapidity grid taken into account
  vector<TH2D*> yieldsSameEventLead; // multi-dielectron events
  vector<TH2D*> yieldsSameEventTail; // multi-dielectron events

  vector<TH1D*> hMassv; // 1GeV bins
  vector<TH1D*> hMassR9ggBBv, hMassR9ggEEv, hMassR9ggBEv; // 1GeV bins
  vector<TH1D*> hMassR9gnBBv, hMassR9gnEEv, hMassR9gnBEv, hMassR9ngBEv; // 1GeV bins
  vector<TH1D*> hMassR9nnBBv, hMassR9nnEEv, hMassR9nnBEv; // 1GeV bins
  vector<TH1D*> hR9MissEventsv;

  vector<TH1D*> hMassBinsv;
  TH1D* hSelEvents = NULL;

  vector<TH1D*> hZpeakv;
  vector<TH1D*> hZpeakBBv, hZpeakBEv, hZpeakEEv;
  vector<DielectronChargeCounter_t*> dccV; 

  // the main result of the macro
  createBaseH2Vec(yieldsBaseH2,"hYield_BaseH2_",inpMgr.sampleNames());
  // create distribution of the same-event candidates
  createBaseH2Vec(yieldsSameEventLead,"hYieldsSameEventLead_BaseH2_",inpMgr.sampleNames());
  createBaseH2Vec(yieldsSameEventTail,"hYieldsSameEventTail_BaseH2_",inpMgr.sampleNames());
  // debug distributions: 1GeV bins
  //createAnyH1Vec(hMassv,"hMass_",inpMgr.sampleNames(),2500,0.,2500.,"M_{ee} [GeV]","counts/1GeV");
  createAnyH1Vec(hMassv,"hMass_",inpMgr.sampleNames(),1490,10.,1500.,"M_{ee} [GeV]","counts/1GeV");

  // mass distributions selecting by golden R9 value
  int nBins_R9=1990;
  double massMin_R9=10.;
  double massMax_R9=2000.;

  createAnyH1Vec(hMassR9ggBBv,"hMassR9ggBB_",inpMgr.sampleNames(),nBins_R9,massMin_R9,massMax_R9,"M_{ee} [GeV]","ggBB counts/1GeV");
  createAnyH1Vec(hMassR9ggEEv,"hMassR9ggEE_",inpMgr.sampleNames(),nBins_R9,massMin_R9,massMax_R9,"M_{ee} [GeV]","ggEE counts/1GeV");
  createAnyH1Vec(hMassR9ggBEv,"hMassR9ggBE_",inpMgr.sampleNames(),nBins_R9,massMin_R9,massMax_R9,"M_{ee} [GeV]","ggBE counts/1GeV");
  createAnyH1Vec(hMassR9nnBBv,"hMassR9nnBB_",inpMgr.sampleNames(),nBins_R9,massMin_R9,massMax_R9,"M_{ee} [GeV]","nnBB counts/1GeV");
  createAnyH1Vec(hMassR9nnEEv,"hMassR9nnEE_",inpMgr.sampleNames(),nBins_R9,massMin_R9,massMax_R9,"M_{ee} [GeV]","nnEE counts/1GeV");
  createAnyH1Vec(hMassR9nnBEv,"hMassR9nnBE_",inpMgr.sampleNames(),nBins_R9,massMin_R9,massMax_R9,"M_{ee} [GeV]","nnBE counts/1GeV");
  createAnyH1Vec(hMassR9gnBBv,"hMassR9gnBB_",inpMgr.sampleNames(),nBins_R9,massMin_R9,massMax_R9,"M_{ee} [GeV]","gnBB counts/1GeV");
  createAnyH1Vec(hMassR9gnEEv,"hMassR9gnEE_",inpMgr.sampleNames(),nBins_R9,massMin_R9,massMax_R9,"M_{ee} [GeV]","gnEE counts/1GeV");
  createAnyH1Vec(hMassR9gnBEv,"hMassR9gnBE_",inpMgr.sampleNames(),nBins_R9,massMin_R9,massMax_R9,"M_{ee} [GeV]","gBnE counts/1GeV");
  createAnyH1Vec(hMassR9ngBEv,"hMassR9ngBE_",inpMgr.sampleNames(),nBins_R9,massMin_R9,massMax_R9,"M_{ee} [GeV]","nBgE counts/1GeV");

  createAnyH1Vec(hR9MissEventsv,"hR9MissEvents_",inpMgr.sampleNames(), 3, 0.5, 3.5, "gg/ng/nn","unassigned barrel/endcap counts");

  // debug distributions for current mass bin
  createBaseH1Vec(hMassBinsv,"hMassBins_",inpMgr.sampleNames());
  // debug: accumulate info about the selected events in the samples
  hSelEvents=createAnyTH1D("hSelEvents","hSelEvents",inpMgr.sampleCount(),0,inpMgr.sampleCount(),"sampleId","event count");
  // collect number of events in the Z-peak
  createAnyH1Vec(hZpeakv,"hZpeak_",inpMgr.sampleNames(),60,60.,120.,"M_{ee} [GeV]","counts/1GeV");
  createAnyH1Vec(hZpeakBBv,"hZpeakBB_",inpMgr.sampleNames(),60,60.,120.,"M_{ee} [GeV]","counts (BB)/1GeV");
  createAnyH1Vec(hZpeakBEv,"hZpeakBE_",inpMgr.sampleNames(),60,60.,120.,"M_{ee} [GeV]","counts (BE)/1GeV");
  createAnyH1Vec(hZpeakEEv,"hZpeakEE_",inpMgr.sampleNames(),60,60.,120.,"M_{ee} [GeV]","counts (EE)/1GeV");

  dccV.reserve(inpMgr.sampleCount());
  for (unsigned int isam=0; isam<inpMgr.sampleCount(); ++isam) {
    TString dccName="dcc_" + inpMgr.sampleName(isam);
    dccV.push_back(new DielectronChargeCounter_t(dccName));
  }
  
  ZeeData_t *data = new ZeeData_t();
  TRandom random;

  //
  //  Diboson backgrounds need to be saved separately, but plotted
  //  together. Detect whether there are separate ww/wz/zz contribution,
  //  and whether merging is needed later. Of course, this relies on
  // the fact that the file data.conf has names ww, wz, zz for those
  // contributions.
  //bool mergeDibosons = (inpMgr.CountDibosonSamples()==3) ? true : false;

  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0; 
  
  if (DYTools::processData(runMode)) {
  for(UInt_t isam=0; isam<inpMgr.sampleCount(); ++isam) {

    TString fname = inpMgr.nTupleFullFileName(isam,systMode);
    if (inpMgr.userKeyValueAsInt("IGNOREDEBUGRUNFORYIELDS")==1) {
      fname.ReplaceAll("_DebugRun","");
    }
    cout << "Processing " << fname << "..." << endl;
    bool isData=(inpMgr.sampleName(isam) == "data") ? 1:0;

    infile = new TFile(fname);
    assert(infile); 

    // Get the TTree and set branch address
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree); 
    eventTree->SetBranchAddress("Events",&data);

    std::cout << "here are " << eventTree->GetEntries() << " entries in " << inpMgr.sampleName(isam) << " sample\n";

    // keep track of multi-dielectron events
    UInt_t mdePrevRunNum=(UInt_t)(-1);
    UInt_t mdePrevEvtNum=mdePrevRunNum;
    double mdeOldMass=-1;
    double mdeOldY=-1;
    double mdeOldWeight=-1;
    int mdeFirstEncounter=0;
    std::vector<double> multiEEmass, multiEEy, multiEEweight;

    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      if ((runMode==DYTools::DEBUG_RUN) && (ientry>100000)) break;
      printProgress(1000000," ientry=",ientry,eventTree->GetEntriesFast());

      eventTree->GetEntry(ientry);
      if (applyEtEtaCut || (outputSystMode==DYTools::UNREGRESSED_ENERGY)) {
#ifdef ZeeData_storeUnregEn
	// study variations due to the regression
	/*
	if (0 && !((data->ptUncorr_1 > data->pt_1) &&
		   (data->scEtUncorr_1 > data->scEt_1) &&
		   (data->ptUncorr_2 > data->pt_2) &&
		   (data->scEtUncorr_2 > data->scEt_2))) {
	  std::cout << "data: uncorr set " << data->ptUncorr_1 << ", " << data->scEtUncorr_1 << "; " << data->ptUncorr_2 << ", " << data->scEtUncorr_2 << "; corr set " << data->pt_1 << ", " << data->scEt_1 << "; " << data->pt_2 << ", " << data->scEt_2 << "\n";
	}
	*/

	data->replace2UncorrEn(0,mdfEnFactor);
#else
	std::cout << "outputSystMode=UNREGRESSED_ENERGY needs that the values are stored in ZeeData\n";
#endif
      }

      // If This is MC, add extra smearing to the mass
      // We apply extra smearing to all MC samples: it is may be
      // not quite right for fake electron backgrounds, but these
      // are not dominant, and in any case we do not have corrections
      // for fake electrons.
      if (applyEScale) {
	int randomized=(isData) ? (applyEScale-1) : 0;
	if (!evtSelector_forOutput.applyEScale(data,isData,randomized)) {
	  std::cout << "error when smearing/scaling an event\n";
	}
      }

      if (applyEtEtaCut) {
	if (!DYTools::goodEtEtaPair(data->scEt_1, data->scEta_1,
				    data->scEt_2, data->scEta_2)) continue;
      }
      Double_t weight = data->weight;


      // Find the 2D bin for this event:
      int massBin = DYTools::findMassBin(data->mass);
      int yBin    = DYTools::findAbsYBin(massBin, data->y);

      if ((massBin==-1) || (yBin==-1)) // out of range
	continue;

      if ((mdePrevRunNum==data->runNum) &&
	  (mdePrevEvtNum==data->evtNum)) {
	if (mdeFirstEncounter) {
	  //yieldsSameEvent[isam]->Fill(mdeOldMass,mdeOldY, mdeOldWeight);
	  mdeFirstEncounter=0;
	  multiEEmass.push_back(mdeOldMass);
	  multiEEy.push_back(mdeOldY);
	  multiEEweight.push_back(mdeOldWeight);
	}
	//yieldsSameEvent[isam]->Fill(data->mass,fabs(data->y), weight);
	multiEEmass.push_back(data->mass);
	multiEEy.push_back(fabs(data->y));
	multiEEweight.push_back(weight);
      }
      else {
	mdeFirstEncounter=1;
	if (multiEEmass.size()) {
	  double maxE=*std::max_element(multiEEmass.begin(),multiEEmass.end());
	  for (unsigned int ii=0; ii<multiEEmass.size(); ++ii) {
	    if (multiEEmass[ii]==maxE) {
	      yieldsSameEventLead[isam]->Fill(multiEEmass[ii],
					      multiEEy[ii], multiEEweight[ii]);
	    }
	    else {
	      yieldsSameEventTail[isam]->Fill(multiEEmass[ii],
					      multiEEy[ii], multiEEweight[ii]);
	    }
	  }
	  multiEEmass.clear();
	  multiEEy.clear();
	  multiEEweight.clear();
	}
      }
      mdePrevRunNum=data->runNum;
      mdePrevEvtNum=data->evtNum;
      mdeOldMass=data->mass;
      mdeOldY=fabs(data->y);
      mdeOldWeight=weight;

      yieldsBaseH2[isam]->Fill(data->mass,fabs(data->y), weight);
      //if (data->q_1==data->q_2) std::cout << "same-sign event weight=" << weight << "\n";
      dccV[isam]->Fill(data,weight);

      hMassv[isam]->Fill(data->mass,weight);
      hMassBinsv[isam]->Fill(data->mass,weight);
      hZpeakv[isam]->Fill(data->mass,weight);

      int isB1=DYTools::isBarrel(data->eta_1);
      int isB2=DYTools::isBarrel(data->eta_2);
      int isE1=DYTools::isEndcap(data->eta_1);
      int isE2=DYTools::isEndcap(data->eta_2);

      if (isB1 && isB2) hZpeakBBv[isam]->Fill(data->mass,weight);
      else if (isE1 && isE2) hZpeakEEv[isam]->Fill(data->mass,weight);
      else if ((isB1 && isE2) || (isE1 && isB2)) hZpeakBEv[isam]->Fill(data->mass,weight);

      int isG1=(data->golden_1==1) ? 1:0;
      int isG2=(data->golden_2==1) ? 1:0;
      if (isG1 && isG2) {
	if (isB1 && isB2) hMassR9ggBBv[isam]->Fill(data->mass,weight);
	else if (isE1 && isE2) hMassR9ggEEv[isam]->Fill(data->mass,weight);
	else if ((isB1 && isE2) || (isE1 && isB2)) hMassR9ggBEv[isam]->Fill(data->mass,weight);
	else hR9MissEventsv[isam]->Fill(1,1.);
      }
      else if ((isG1 && !isG2) || (!isG1 && isG2)) {
	if (isB1 && isB2) hMassR9gnBBv[isam]->Fill(data->mass,weight);
	else if (isE1 && isE2) hMassR9gnEEv[isam]->Fill(data->mass,weight);
	else if ((isB1 &&  isG1 && isE2 && !isG2) ||
		 (isE1 && !isG1 && isB2 &&  isG2)) {
	  hMassR9gnBEv[isam]->Fill(data->mass,weight);
	}
	else if ((isB1 && !isG1 && isE2 &&  isG2) ||
		 (isE1 &&  isG1 && isB2 && !isG2)) {
	  hMassR9ngBEv[isam]->Fill(data->mass,weight);
	}
	else hR9MissEventsv[isam]->Fill(2,1.);
      }
      else if (!isG1 && !isG2) {
	if (isB1 && isB2) hMassR9nnBBv[isam]->Fill(data->mass,weight);
	else if (isE1 && isE2) hMassR9nnEEv[isam]->Fill(data->mass,weight);
	else if ((isB1 && isE2) ||
		 (isB2 && isE1)) {
	  hMassR9nnBEv[isam]->Fill(data->mass,weight);
	}
	else hR9MissEventsv[isam]->Fill(3,1.);
      }

      if ( DYTools::validMass(data->mass) ) {
	hSelEvents->Fill(isam, weight);
      }
    }
    delete infile;
    infile=0, eventTree=0;
  }

  if (0) {
    TFile file("test_new.root","recreate");
    for (unsigned int i=0; i<hMassv.size(); ++i) {
      //hMassv[i]->Write(inpMgr.sampleName(i));
      hMassBinsv[i]->Write(inpMgr.sampleName(i));
    }
    file.Close();
  }

  // Prepare the correct distribution
  if (!convertBaseH2actualVec(yieldsBaseH2, yields,"hYield_",inpMgr.sampleNames(),1)) return retCodeError;

  }
  else {
    // to be able to load data
    createBaseH2Vec(yields,"hYield_",inpMgr.sampleNames());
  }



  if (DYTools::processData(runMode)) {
    std::cout << "saving to <" << yieldFullName << ">\n";
    TFile file(yieldFullName,"recreate");
    int res=1;
    if (res) res=saveVec(file,yieldsBaseH2,"yieldsBaseH2");
    if (res) res=saveVec(file,yields,"yields");
    if (res) res=saveVec(file,yieldsSameEventLead,"yieldsSameEventLead");
    if (res) res=saveVec(file,yieldsSameEventTail,"yieldsSameEventTail");
    if (res) res=saveVec(file,hMassv,"mass_1GeV_bins");
    if (res) res=saveVec(file,hMassBinsv,"mass_analysis_bins");
    if (res) res=saveHisto(file,hSelEvents,"");
    if (res) res=saveVec(file,hZpeakv,"mass_Zpeak_1GeV");
    if (res) res=saveVec(file,hZpeakBBv,"mass_ZpeakBB_1GeV");
    if (res) res=saveVec(file,hZpeakBEv,"mass_ZpeakBE_1GeV");
    if (res) res=saveVec(file,hZpeakEEv,"mass_ZpeakEE_1GeV");
    if (res) res=saveVec(file,dccV,"ee_charge_counter");
    if (res) res=saveVec(file,hMassR9ggBBv,"mass_R9ggBB");
    if (res) res=saveVec(file,hMassR9ggEEv,"mass_R9ggEE");
    if (res) res=saveVec(file,hMassR9ggBEv,"mass_R9ggBE");
    if (res) res=saveVec(file,hMassR9nnBBv,"mass_R9nnBB");
    if (res) res=saveVec(file,hMassR9nnEEv,"mass_R9nnEE");
    if (res) res=saveVec(file,hMassR9nnBEv,"mass_R9nnBE");
    if (res) res=saveVec(file,hMassR9gnBBv,"mass_R9gnBB");
    if (res) res=saveVec(file,hMassR9gnEEv,"mass_R9gnEE");
    if (res) res=saveVec(file,hMassR9gnBEv,"mass_R9gnBE");
    if (res) res=saveVec(file,hMassR9ngBEv,"mass_R9ngBE");
    if (res) res=saveVec(file,hR9MissEventsv,"R9_misses");
    if (res) writeBinningArrays(file,"prepareYields");
    file.Close();
    if (!res) {
      std::cout << "error occurred during save to file <" << yieldFullName << ">\n";
      return retCodeError;
    }
  }
  else {
    std::cout << "runMode= loadData" << std::endl;
    std::cout << "loading from <" << yieldFullName << ">\n";
    TFile file(yieldFullName,"read");
    int res=1;
    if (res) res=checkBinningArrays(file);
    if (res) res=loadVec(file,yieldsBaseH2,"yieldsBaseH2");
    if (res) res=loadVec(file,yields,"yields");
    if (res) res=loadVec(file,yieldsSameEventLead,"yieldsSameEventLead");
    if (res) res=loadVec(file,yieldsSameEventTail,"yieldsSameEventTail");
    if (res) res=loadVec(file,hMassv,"mass_1GeV_bins");
    if (res) res=loadVec(file,hMassBinsv,"mass_analysis_bins");
    if (res) res=loadHisto(file,&hSelEvents,"");
    if (res) res=loadVec(file,hZpeakv,"mass_Zpeak_1GeV");
    if (res) res=loadVec(file,hZpeakBBv,"mass_ZpeakBB_1GeV");
    if (res) res=loadVec(file,hZpeakBEv,"mass_ZpeakBE_1GeV");
    if (res) res=loadVec(file,hZpeakEEv,"mass_ZpeakEE_1GeV");
    if (res) res=loadVec(file,dccV,"ee_charge_counter");
    if (res) res=loadVec(file,hMassR9ggBBv,"mass_R9ggBB");
    if (res) res=loadVec(file,hMassR9ggEEv,"mass_R9ggEE");
    if (res) res=loadVec(file,hMassR9ggBEv,"mass_R9ggBE");
    if (res) res=loadVec(file,hMassR9nnBBv,"mass_R9nnBB");
    if (res) res=loadVec(file,hMassR9nnEEv,"mass_R9nnEE");
    if (res) res=loadVec(file,hMassR9nnBEv,"mass_R9nnBE");
    if (res) res=loadVec(file,hMassR9gnBBv,"mass_R9gnBB");
    if (res) res=loadVec(file,hMassR9gnEEv,"mass_R9gnEE");
    if (res) res=loadVec(file,hMassR9gnBEv,"mass_R9gnBE");
    if (res) res=loadVec(file,hMassR9ngBEv,"mass_R9ngBE");
    if (res) res=loadVec(file,hR9MissEventsv,"R9_misses");
    file.Close();
    if (!res) {
      std::cout << "error occurred during load of file <" << yieldFullName << ">\n";
      return retCodeError;
    }
  }

  std::cout << "className=<" << hSelEvents->Class()->ClassName() << ">\n";

  if (0) {
    HERE("create canvas");
    TCanvas *cx=new TCanvas("cx","cx",600,600);
    //hSelEvents->Draw();
    //hMassBinsv[8]->Draw();
    //TH1D *hh=dccV[0]->calcH2TotProfile(NULL,1,1,NULL);
    //hh->Draw("LP");
    //hh->Print("ranges");
    
    dccV[dccV.size()-1]->calcH2SSFracProfile(NULL,NULL,1,1,NULL)->DrawClone("LP");
    //dccV[9]->calcH2SSFrac(NULL,1,NULL)->DrawClone("COLZ");
    cx->Update();
    HERE("cx->Update() called");
  }

  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================
  /*

  // Buffers for labels and comments
  //char ylabel[100];   // y-axis label
  char lumitext[50];
  if(lumi>0) {
    if(lumi<1) { sprintf(lumitext,"#int#font[12]{L}dt = %.0f nb^{-1}",1000.*lumi); }
    if(lumi<1000) { sprintf(lumitext,"#int#font[12]{L}dt = %.0f nb^{-1}",1000.*lumi); }
    else       { sprintf(lumitext,"#int#font[12]{L}dt = %.3g fb^{-1}",lumi/1000.0); }
  }


  if (1) {
    TCanvas *cx=new TCanvas("cx","cx",900,900);
    /
    TH2D *h2Tot_data=dccV[0]->calcH2Tot(1);
    h2Tot_data->Print();
    TH2D *h2SSFrac_data=dccV[0]->calcH2SSFrac(h2Tot_data,1);
    TH2D *h2Tot_Zee=dccV[dccV.size()-1]->calcH2Tot(1);
    TH2D *h2SSFrac_Zee=dccV[dccV.size()-1]->calcH2SSFrac(h2Tot_Zee,1);

    std::cout << "h2Tot_data->Integral =" << h2Tot_data->Integral() << "\n";
    std::cout << "h2Tot_Zee->Integral =" << h2Tot_Zee->Integral() << "\n";
    std::cout << "hMass_data->Integral = " << hMassv[0]->Integral() << "\n";
    std::cout << "hMass_Zee->Integral = " << hMassv[hMassv.size()-1]->Integral() << "\n";
    /

    std::cout << "sumDataSSTot=" << sumDataSSTot << ", sumDataTot=" << sumDataTot << ", ratio=" << sumDataSSTot/sumDataTot << "\n";
    std::cout << "sumDataSSTot1=" << sumDataSSTot1 << ", sumDataTot1=" << sumDataTot1 << ", ratio=" << sumDataSSTot1/sumDataTot1 << "\n";
    std::cout << "sumDataSSTot2=" << sumDataSSTot2 << ", sumDataTot2=" << sumDataTot2 << ", ratio=" << sumDataSSTot2/sumDataTot2 << "\n";

    cx->Divide(3,3);
    if (!DYTools::study2D) {
      int ipad=1;
      for (unsigned int isam=0; isam<dccV.size(); ++isam) {
	if (isam==dccV.size()-2) continue;
	cx->cd(ipad);
	TPad *pad=(TPad*)cx->GetPad(ipad); pad->SetLogx();
	ipad++;

	TH2D *h2SSFrac = dccV[isam]->calcH2SSFrac(NULL,1);
	TH1D *h1SSFrac=getProfileX(h2SSFrac,1,h2SSFrac->GetName() + TString("_profile"),1);

	//if (0) { // get average fraction
	  TH2D *h2Tot=dccV[isam]->calcH2Tot(1);
	  TH1D *h1Tot=getProfileX(h2Tot,1,h2Tot->GetName() + TString("_profile"));
	  //double tot=h1Tot->Integral();
	  TH2D *h2ss=dccV[isam]->getH2SS();
	  TH1D *h1ss=getProfileX(h2ss,1,h2ss->GetName() + TString("_profile"));
	  std::cout << " sample=" << snamev[isam] << ", sameSign=" << h1ss->Integral() << ", total=" << h1Tot->Integral() << ", avgFracSameSign=" << h1ss->Integral()/h1Tot->Integral() << "\n";
	  //}
	for (int the_case=-1; the_case<=2; ++the_case) {
	  std::cout << "  the_case=" << the_case << ", integral=" << dccV[isam]->calcIntegral(1,the_case) << "\n";
	}


	TString yAxisLabel;
	TH1D *histo=h1SSFrac; yAxisLabel="same sign ee fraction";
	histo=h1ss; yAxisLabel="same sign ee count";

	removeError(histo);
	///TString title=histo->GetTitle();
	//title.ReplaceAll("dcc_",""); 
	//title.ReplaceAll("_SSfrac_profile","");
	histo->SetTitle(snamev[isam]);
	histo->GetYaxis()->SetTitle(yAxisLabel);
	histo->GetYaxis()->SetTitleOffset(1.4);
	histo->Draw();
      }
    }
    else { // 2D
      AdjustFor2DplotWithHeight(cx);
      int ipad=1;
      for (unsigned int isam=0; isam<dccV.size(); ++isam) {
	if (isam==dccV.size()-2) continue;
	cx->cd(ipad);
	ipad++;
	TH2D *h2SSFrac = dccV[isam]->calcH2SSFrac(NULL,1);
	h2SSFrac->Draw("COLZ");
      }
    }
  }
  return retCodeStop;

  // ---------- END Of TEST ------------------


  // Merge diboson histograms if needed
  int clone_idx=(hasData) ? 1:0;
  TH1F *hMassBinsDibosons = (TH1F*)hMassBinsv[clone_idx]->Clone("hMassBinsDibosons");
  TH1F *hMassDibosons = (TH1F*)hMassv[clone_idx]->Clone("hMassDibosons");
  hMassBinsDibosons->Reset();
  hMassDibosons->Reset();
  Int_t colorDibosons = 1;
  TString labelDibosons = "WW/WZ/ZZ";
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    if( snamev[isam] == "ww" || snamev[isam] == "wz" || snamev[isam] == "zz"){
      hMassDibosons->Add(hMassv[isam]);
      hMassBinsDibosons->Add(hMassBinsv[isam]);
      // Use color of the last diboson entry
      colorDibosons = samplev[isam]->color;
    }
  }

  // Additional normalization for MC:
  //
  // Ideally, we would normalize all MC samples to data luminosity
  // In practice, however, it is not easy because of two reasons:
  //  - the luminosity is known with an error (systematic shift of 6% is
  //       for example suspected in mid-2011)
  //  - data/MC scale factors for efficiency to select events may
  //       move normalization off by another 5%
  // Therefore, we normalize total MC to the data Z peak. This gives us
  // the scale factor that is applied to all samples. 
  // In the following calculation it is assumed that the first histogram
  // is data and the last is signal MC.  
  TH1F *totalMCMass = (TH1F*)hMassv[0]->Clone("totalMCMass");
  totalMCMass->Reset();
  for(UInt_t isam=(hasData)?1:0; isam<samplev.size(); isam++) {
    totalMCMass->Add(hMassv[isam]);
  }
  double massNormMin = 60.0;
  double massNormMax = 120.0;
  double dataOverMc = hMassv[0]->Integral(hMassv[0]->FindBin(massNormMin+0.001),
					  hMassv[0]->FindBin(massNormMax-0.001)) /
    totalMCMass->Integral(totalMCMass->FindBin(massNormMin+0.001),
			  totalMCMass->FindBin(massNormMax-0.001));
  printf("data to MC extra correction from Z peak normalization: %f\n",dataOverMc);

  double dataOverMcEachBin[DYTools::nMassBins+1];
  for (int i=0; i<DYTools::nMassBins; i++)
    {
      dataOverMcEachBin[i] = hMassv[0]->Integral(hMassv[0]->FindBin(DYTools::massBinLimits[i]+0.001),hMassv[0]->FindBin(DYTools::massBinLimits[i+1]-0.001)) /
	totalMCMass->Integral(totalMCMass->FindBin(DYTools::massBinLimits[i]+0.001),totalMCMass->FindBin(DYTools::massBinLimits[i+1]-0.001));
      printf("data to MC %i bin norm: %f\n",i,dataOverMcEachBin[i]);
    }

  //std::cout << "\n\nsetting dataOverMc=1\n"; dataOverMc=1; // for 1_to_1 comparison with 1D
  
  // Rescale all MC samples. This is not totally proper for fake lepton
  // backgrounds, but ok for backgrounds with true leptons, and those are dominant
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hMassv[isam]->Scale(dataOverMc);
    hMassBinsv[isam]->Scale(dataOverMc);
    printf("  MC %s IS RESCALED for plotting\n", snamev[isam].Data());
  }

  hMassDibosons->Scale(dataOverMc); printf(" compound MC %s IS RESCALED for plotting\n", hMassDibosons->GetName());
  hMassBinsDibosons->Scale(dataOverMc); printf("  compound MC %s IS RESCALED for plotting\n", hMassBinsDibosons->GetName());

  //
  // Prepare outputDir and the plot file
  //

  TString outputDirYields(outputDir.Data());
  if (performPUReweight) outputDirYields.ReplaceAll("selected_events","yields");
  else outputDirYields.ReplaceAll("selected_events","yields_noPU");
  gSystem->mkdir(outputDirYields,kTRUE);
  TString fNameOutYieldPlots(outputDirYields+TString("/yield_plots") + DYTools::analysisTag);
  fNameOutYieldPlots += ".root";
  TFile *fYieldPlots = new TFile( fNameOutYieldPlots, "recreate" );
  if (!fYieldPlots) {
    std::cout << "Failed to create a file <" << fNameOutYieldPlots << ">\n";
    throw 2;
  }

  //
  // Draw mass spectrum without rapidity binning
  //


  // First, draw the mass histograms with fine mass binning
  DrawMassPeak(hMassv, samplev, snamev, hMassDibosons, hasData, mergeDibosons, labelDibosons, colorDibosons, lumi, lumitext, 0, fYieldPlots);

  // Second, draw the mass histograms with the mass binning used in the analysis
  DrawMassPeak(hMassBinsv, samplev, snamev, hMassBinsDibosons, hasData, mergeDibosons, labelDibosons, colorDibosons, lumi, lumitext, 1, fYieldPlots);

  // Draw the flattened figure (Y histograms for different mass regions)
  DrawFlattened(yields, yieldsSumw2, samplev, snamev, hasData, mergeDibosons, labelDibosons, colorDibosons, lumi, lumitext, fYieldPlots);

  // Draw rapidity in mass slices 
  if (DYTools::study2D==1)
    {
       Draw6Canvases(yields, yieldsSumw2, samplev, snamev, hasData, dataOverMc, dataOverMcEachBin, 1, 0, fYieldPlots);
       Draw6Canvases(yields, yieldsSumw2, samplev, snamev, hasData, dataOverMc, dataOverMcEachBin, 1, 1, fYieldPlots);
       Draw6Canvases(yields, yieldsSumw2, samplev, snamev, hasData, dataOverMc, dataOverMcEachBin, 1,-1, fYieldPlots);
       Draw6Canvases(yields, yieldsSumw2, samplev, snamev, hasData, dataOverMc, dataOverMcEachBin, 0, 0, fYieldPlots);
       Draw6Canvases(yields, yieldsSumw2, samplev, snamev, hasData, dataOverMc, dataOverMcEachBin, 0, 1, fYieldPlots);
       Draw6Canvases(yields, yieldsSumw2, samplev, snamev, hasData, dataOverMc, dataOverMcEachBin, 0,-1, fYieldPlots);
    }

  fYieldPlots->Close();
*/

  //--------------------------------------------------------------------------------------------------------------
  // Save the results
  //==============================================================================================================

  // Pack info into writable objects
  //TVectorD massBinLimits(nMassBins+1);
  //TVectorD rapidityBinning(nMassBins+1);
  //for(int i=0; i <= nMassBins; i++){
  //  massBinLimits(i) = DYTools::massBinLimits[i];
  //  rapidityBinning(i) = DYTools::nYBins[i];
  //}

  // This dummy object is only needed to convey the number
  // of samples considered. The method below is rather convoluted,
  // but I do not know a better one. Ideally, we would just store
  // a list of strings, each string containing the sample name.
  /*
  TVectorD dummySampleCount(sampleTags.size());
  dummySampleCount = 0;

  gSystem->mkdir(outputDirYields,kTRUE);
  TString fNameOutYields(outputDirYields+TString("/yields") + DYTools::analysisTag);
  fNameOutYields += ".root";
  TFile fYields( fNameOutYields, "recreate" );
  //massBinLimits      .Write("massBinLimits");
  //rapidityBinning    .Write("rapidityBinning");
  writeBinningArrays(fYields,"prepareYields");
  dummySampleCount   .Write("dummySampleCount");
  writeBinningArrays(fYields,"prepareYields");
  char objName[100];
  for(UInt_t isam = 0; isam < yields.size(); isam++){
    sprintf(objName,"sample_name_%i",isam);
    TObjString *sampleNameStorable = new TObjString( sampleTags.at(isam) );
    sampleNameStorable->Write(objName);
    sprintf(objName,"yields_%s",sampleTags.at(isam).Data());
    yields[isam]       ->Write(objName);
    sprintf(objName,"yieldsSumw2_%s",sampleTags.at(isam).Data());
    yieldsSumw2[isam]  ->Write(objName);
  }
  fYields.Close();
  */

  /*
  // Save mass histograms into a separate file for a direct comparison 
  // with DrellYan(1D) package
  TString fNameOutHists(outputDirYields+"/massHist");
  fNameOutHists += DYTools::analysisTag + TString(".root");
  TFile fMassHists(fNameOutHists,"recreate");
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    hMassBinsv[isam]->Write(snamev[isam]);
  }
  std::cout << "file <" << fNameOutHists << "> created\n";
  fMassHists.Close();
  */


  //--------------------------------------------------------------------------------------------------------------
  // Print R9 plots in categories
  //==============================================================================================================

  if (0) {
    double ymin=1.,ymax=1e6;
    TString titleBase="Reg.En., ";
    TString fileBase="fig-RegEn--";
    if (outputSystMode==DYTools::UNREGRESSED_ENERGY) {
      titleBase="UnReg.En, ";  fileBase="fig-UnRegEn--"; 
    }
    else if (outputSystMode==DYTools::APPLY_ESCALE) {
      titleBase="Reg.En.Mdf., "; fileBase="fig-RegEnMdf--";
    }
    for (int plot_case=0; plot_case<10; ++plot_case) {
      TString plotTag;
      TString titleTag;
      switch(plot_case) {
      case 0: plotTag="R9ggBB"; titleTag="barrel-barrel, both golden"; break;
      case 1: plotTag="R9ggEE"; titleTag="endcap-endcap, both golden"; break;
      case 2: plotTag="R9ggBE"; titleTag="barrel-endcap, both golden"; break;
      case 3: plotTag="R9nnBB"; titleTag="barrel-barrel, both showering"; break;
      case 4: plotTag="R9nnEE"; titleTag="endcap-endcap, both showering"; break;
      case 5: plotTag="R9nnBE"; titleTag="barrel-endcap, both showering"; break;
      case 6: plotTag="R9gnBB"; titleTag="barrel-barrel, golden & showering"; break;
      case 7: plotTag="R9gnEE"; titleTag="endcap-endcap, golden & showering"; break;
      case 8: plotTag="R9gnBE"; titleTag="golden in barrel, showering in endcap"; break;
      case 9: plotTag="R9ngBE"; titleTag="showering in barrel, golden in endcap"; break;
      default:
	std::cout << "plot_case=" << plot_case << " is not ready for R9xxYY studies\n";
	return retCodeError;
      }
      vector<TH1D*> *histos=NULL;
      switch(plot_case) {
      case 0: histos=&hMassR9ggBBv; ymax=1e5; break;
      case 1: histos=&hMassR9ggEEv; ymax=1e5; break;
      case 2: histos=&hMassR9ggBEv; ymax=1e5; break;
      case 3: histos=&hMassR9nnBBv; ymax=1e6; break;
      case 4: histos=&hMassR9nnEEv; ymax=1e5; break;
      case 5: histos=&hMassR9nnBEv; ymax=1e6; break;
      case 6: histos=&hMassR9gnBBv; ymax=1e6; break;
      case 7: histos=&hMassR9gnEEv; ymax=1e5; break;
      case 8: histos=&hMassR9gnBEv; ymax=1e5; break;
      case 9: histos=&hMassR9ngBEv; ymax=1e5; break;
      default:
	std::cout << "plot_case=" << plot_case << " is not ready for R9xxYY studies (2)\n";
	return retCodeError;
      }
      ymax=1e6;


      TString cpName=TString("cp_") + plotTag;
      TString cpTitle=titleBase;
      cpTitle.Append(titleTag);
      ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,cpName,cpTitle,
			  "#it{M}_{ee} [GeV]","counts/GeV","ratio");
      cp.SetYTitleSize(0.06,1.2);
      cp.SetLogx(1);
      cp.SetLogy(1);
      cp.SetLogx(0); cp.SetXRange(59.99,120.01);
      cp.SetYRange(ymin,ymax);

      TString hName=Form("hSumMC_%d",plot_case);
      TH1D *h=(TH1D*)(*histos)[0]->Clone(hName);
      h->Reset();
      h->SetDirectory(0);
      h->SetTitle(hName);

      for (unsigned int ih=0; ih<histos->size(); ++ih) {
	if (ih==0) cp.AddHist1D((*histos)[ih],inpMgr.sampleName(ih),"LP",kBlack);
	else {
	  h->Add((*histos)[ih],1.);
	  int color=inpMgr.sampleInfo(ih)->colors[0];
	  cp.AddToStack((*histos)[ih],inpMgr.sampleName(ih),color);
	}
      }
      h->SetMarkerStyle(24);
      cp.AddHist1D(h,"all MC","LP skip",kRed,1,1,1);

      TString canvName=Form("cx%d",plot_case);
      TString canvTitle=Form("canv_%d",plot_case);
      TCanvas *cx=new TCanvas(canvName,canvTitle,800,700);
      cp.Prepare2Pads(cx);
      SetSideSpaces(cx,0.0,0.2,0.,0.01);
      cp.Draw(cx);
      cp.ChangeLegendPos(0.2,0.,0.1,0.);
      cx->Update();

      TString figName=fileBase + plotTag;
      SaveCanvas(cx,figName);
    }
  }


  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================

  std::cout << std::endl;
  std::cout << "*" << std::endl;
  std::cout << "* SUMMARY" << std::endl;
  std::cout << "*--------------------------------------------------" << std::endl;
  std::cout << std::endl;
  

  double sumData_Zpeak=0.;
  double sumZee_Zpeak=0.;
  double dataOverMC=0.;

  if (1) {
    std::cout << "   Data: " << setprecision(1) << fixed << hSelEvents->GetBinContent(1) << " Drell-Yan events!\n";
    const CSample_t *sample = inpMgr.sampleInfo(0);
    for(UInt_t ifile=0; ifile<sample->size(); ifile++) {
      std::cout << "     " << sample->getFName(ifile) << "\n";
    }
  }
  std::cout << std::endl;
  if(inpMgr.sampleCount()>1) {
    std::cout << "   MC:\n";
    for(UInt_t isam=1; isam<inpMgr.sampleCount(); isam++) {
      const CSample_t *sample = inpMgr.sampleInfo(isam);
      std::cout << "      " << inpMgr.sampleName(isam) << " " << setw(8) << setprecision(3) << fixed << hSelEvents->GetBinContent(isam+1) << " +/- ";
      std::cout << setw(4) << setprecision(3) << fixed << hSelEvents->GetBinError(isam+1) << "\n";
      for(UInt_t ifile=0; ifile<sample->size(); ifile++) {
	std::cout << "     " << sample->getFName(ifile) << "\n";
      }
      std::cout << std::endl;
    }
  }
  std::cout << std::endl;


  //
  // Summary printout in mass bins, integrated over rapidity
  //
  // Add yields over rapidity bins
  double totalData            [DYTools::nMassBins];
  double totalSignalMC        [DYTools::nMassBins];
  double totalSignalMCError   [DYTools::nMassBins];
  double totalBg     [DYTools::nMassBins];
  double totalBgError[DYTools::nMassBins];

  for( int im=0; im<DYTools::nMassBins; im++){
    totalData         [im] = 0;
    totalSignalMC     [im] = 0;
    totalSignalMCError[im] = 0;
    totalBg           [im] = 0;
    totalBgError      [im] = 0;
    int isZpeak=((DYTools::massBinLimits[im]>=60.) && 
		 (DYTools::massBinLimits[im]<120.)) ? 1:0;
    for(int iy = 0; iy < DYTools::nYBins[im]; iy++){
      for( UInt_t isam = 0; isam < yields.size(); isam++){
	double yield_val=yields[isam]->GetBinContent(im+1,iy+1);
	double yield_err=yields[isam]->GetBinError  (im+1,iy+1);
	if( inpMgr.sampleName(isam) == TString("data") ){
	  totalData[im] += yield_val;
	  if (isZpeak) sumData_Zpeak += yield_val;
	}else if ( inpMgr.sampleName(isam) == TString("zee") ){
	  totalSignalMC[im] += yield_val;
	  totalSignalMCError[im] += yield_err*yield_err;
	  if (isZpeak) sumZee_Zpeak+= yield_val;
	}else{
	  // what remains are background samples
	  totalBg[im] += yield_val;
	  totalBgError[im] += yield_err*yield_err;
	}
      } // end loop over samples
    } // end loop over rapidity bins
    totalBgError[im] = sqrt( totalBgError[im] );
    totalSignalMCError[im] = sqrt( totalSignalMCError[im] );
  } // end loop over mass bins

  if (sumZee_Zpeak>double(0.)) dataOverMC= sumData_Zpeak/sumZee_Zpeak;

  printf("Printout of the data, MC signal and MC backgrounds integrated over rapidity\n");
  for (int scale=0; scale<2; ++scale) {
    double factor=(scale) ? dataOverMC : 1.;
    printf("     mass bin        data      MC-signal     MC-backgr\n");
    for(int im = 0; im < DYTools::nMassBins; im++){
      printf("%5.0f-%5.0f GeV: ", DYTools::massBinLimits[im],
	     DYTools::massBinLimits[im+1]);
      printf(" %7.0f", totalData[im]);
      printf(" %9.1f+-%6.1f", totalSignalMC[im]*factor, 
	                      totalSignalMCError[im]*factor);
      printf(" %7.2f+-%5.2f", totalBg[im]*factor, 
	                      totalBgError[im]*factor);
      printf("\n");
    }
    if (scale==0) printf("Note: these MC numbers are not rescaled!\n");
    else printf("Note: these MC numbers are scaled by %8.5lf\n",dataOverMC);
    printf("\n");
  }

  if (1) {
    for (int scale=0; scale<2; ++scale) {
      double factor=(scale) ? dataOverMC : 1.;
      if (factor==double(0)) {
	std::cout << "dataOverMC scaling factor is 0\n";
	continue;
      }
      // A different view of background table
      printf("\n\nPrintout of the backgrounds for all mass bins, view II\n");
      if (scale) printf(" (scaled by %8.5lf)\n",factor);
      else printf(" (not scaled)\n");

      printf("            ");
      for(UInt_t isam=0; isam<inpMgr.sampleCount(); isam++) {
	printf(" %14s ",inpMgr.sampleName(isam).Data());
      }
      printf("           total          fraction\n");
      for(int ibin=0; ibin<DYTools::nMassBins; ibin++){
	printf("%5.1f-%5.1f GeV: ",
	       hMassBinsv[0]->GetXaxis()->GetBinLowEdge(ibin+1),
	       hMassBinsv[0]->GetXaxis()->GetBinUpEdge(ibin+1));
	// Data:
	printf(" %7.0f+-%5.0f ",hMassBinsv[0]->GetBinContent(ibin+1),hMassBinsv[0]->GetBinError(ibin+1) );
	// Individual MC samples
	double total=0., totalError=0.;
	for(UInt_t isam=1; isam<inpMgr.sampleCount(); isam++) {
	  double thisContent = factor*hMassBinsv[isam]->GetBinContent(ibin+1);
	  double thisError = factor*hMassBinsv[isam]->GetBinError(ibin+1);
	  printf(" %7.2f+-%5.2f ",thisContent, thisError);
	  if ( (isam!=0) && (inpMgr.sampleName(isam)!=TString("zee"))) {
	    total+= thisContent;
	    totalError+=thisError*thisError;
	  }
	}
	totalError = sqrt(totalError);
	// Total
	printf("  %8.2f+-%6.2f",total, totalError);
	printf("  %5.1f\n",100*total/hMassBinsv[0]->GetBinContent(ibin+1));
      }
    }
  }

  /*
  if (1 && hZpeakv.size()) {
    std::cout << "Zpeak region average\n";
    unsigned int idx=(hasData) ? 1 : 0;
    TH1F *hTotMC=  (TH1F*)hZpeakv[idx]->Clone("hZpeak_all");
    hTotMC->Reset();
    hTotMC->SetDirectory(0);
    TH1F *hZee=NULL;
    for (unsigned int i=0; i<hZpeakv.size(); ++i) {
      printf(" mean in <%s> is %6.4lf pm %6.4lf\n",hZpeakv[i]->GetName(),hZpeakv[i]->GetMean(),hZpeakv[i]->GetMeanError());
      if (i>=idx) hTotMC->Add(hZpeakv[i]);
      if (snamev[i] == "zee") {
	hZee=(TH1F*)hZpeakv[i]->Clone("hZpeak_zee_clone");
	hZee->SetDirectory(0);
      }
    }
    printf(" mean in <%s> is %6.4lf pm %6.4lf\n",hTotMC->GetName(),hTotMC->GetMean(),hTotMC->GetMeanError());

    if (!hZee) std::cout << "failed to identify Zee\n";
    else if (!hasData) std::cout << "there no experimental data points available\n";
    else {
      TH1F* hData=(TH1F*) hZpeakv[0]->Clone("hZpeak_data_clone");
      hData->SetDirectory(0);

      TCanvas *cZp=MakeCanvas("cZpeak","", 600,700);
      ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,"compPlot","",
			  "mass [GeV]","counts", "MC/data");
      std::cout << "hZee normalization factor=" << (hData->Integral()/hZee->Integral()) << "\n";
      std::cout << "hTotMC normalization factor=" << (hData->Integral()/hTotMC->Integral()) << "\n";
      hZee->Scale( hData->Integral()/hZee->Integral() );
      hTotMC->Scale( hData->Integral()/hTotMC->Integral() );
      cp.Prepare2Pads(cZp);
      //cp.SetLogy();
      cp.SetRatioYRange(0.8,1.1);
      cp.AddHist1D(hData,"data", "LPE", kBlack, 1,0,1);
      cp.AddHist1D(hZee,"MC (norm)", "hist", kBlue, 1,0,1);
      cp.AddHist1D(hTotMC,"tot MC (norm)","hist",kRed, 1,0,1);
      cp.Draw(cZp,false,"png");
      SaveCanvas(cZp,cZp->GetName());
    }
  }
    */

  std::cout << "produced yieldFile=<" << yieldFullName << ">\n";

  //gBenchmark->Show("prepareYields");
  ShowBenchmarkTime("prepareYields");

  // memory clean-up
  ClearVec(yieldsBaseH2);
  ClearVec(yields);
  ClearVec(yieldsSameEventLead);
  ClearVec(yieldsSameEventTail);
  ClearVec(hMassv);
  ClearVec(hMassR9ggBBv);
  ClearVec(hMassR9ggEEv);
  ClearVec(hMassR9ggBEv);
  ClearVec(hMassR9gnBBv);
  ClearVec(hMassR9gnEEv);
  ClearVec(hMassR9gnBEv);
  ClearVec(hMassR9ngBEv);
  ClearVec(hMassR9nnBBv);
  ClearVec(hMassR9nnEEv);
  ClearVec(hMassR9nnBEv);
  ClearVec(hR9MissEventsv);
  ClearVec(hMassBinsv);
  ClearVec(hZpeakv);
  ClearVec(dccV);
  if (data) delete data;


  return retCodeOk;
}

