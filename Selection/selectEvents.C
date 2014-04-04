// 2014.03.04 Based on selectEvents.C
// The name of the created file has extra 'R9'

//================================================================================================
//
// Z->e e selection macro
//
//  * plots distributions associated with selected events
//  * prints list of selected events from data
//  * outputs ROOT files of events passing selection for each sample, 
//    which can be processed by plotSelect.C
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2D.h>
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

using namespace std;

#include "../Include/CPlot.hh"          // helper class for plots
#include "../Include/MitStyleRemix.hh"  // style settings for drawing
#include "../Include/MyTools.hh"        // miscellaneous helper functions
#include "../Include/DYTools.hh"
#include "../Include/TriggerSelection.hh"

// define structures to read in ntuple
#include "../Include/EWKAnaDefs.hh"
#include "../Include/TEventInfo.hh"
#include "../Include/TGenInfo.hh"
#include "../Include/TDielectron.hh"
#include "../Include/TVertex.hh"

// Helper functions for Electron ID selection
#include "../Include/EleIDCuts.hh" 

#include "../Include/AccessOrigNtuples.hh"
#include "../Include/EventSelector.hh"
#include "../Include/EventWeight.hh"
#include "../Include/InputFileMgr.hh"

#endif

// define structure for output ntuple
#include "../Include/ZeeData.hh"

//=== FUNCTION DECLARATIONS ======================================================================================

//=== MAIN MACRO =================================================================================================

int selectEventsR9(const TString conf, 
		   DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN,
		   DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST)
{  
  gBenchmark->Start("selectEvents");


  {
    DYTools::printExecMode(runMode,systMode);
    const int debug_print=1;
    if (!DYTools::checkSystMode(systMode,debug_print,4, DYTools::NO_SYST, DYTools::ESCALE_STUDY, DYTools::ESCALE_STUDY_RND,DYTools::LOWER_ET_CUT)) 
      return retCodeError;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================

  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(conf)) return retCodeError;

  int keepFirstAndLast=0;
  unsigned int keepSample=-1;
  if (keepFirstAndLast){ 
    if (!inpMgr.KeepFirstAndLastSample()) {
      std::cout << "failed to eliminate the middle samples" << std::endl;
      return retCodeError;
    }
    std::cout << "kept only 1st and last sample:\n";
    inpMgr.Print();
  }
  else if (keepSample<inpMgr.sampleCount()) {
    if (!inpMgr.KeepOnlySample(keepSample)) {
      std::cout << "failed to eliminate unneeded samples" << std::endl;
      return retCodeError;
    }
    std::cout << "kept only sample #" << keepSample << "\n";
    inpMgr.Print();
  }
  
  //return retCodeStop;


  // Construct eventSelector, update inpMgr and plot directory
  TString extraTag="R9";
  EventSelector_t evtSelector(inpMgr,runMode,systMode,
			      extraTag,"", EventSelector::_selectDefault);

  // Event weight handler
  EventWeight_t evWeight;
  evWeight.init(inpMgr.puReweightFlag(),inpMgr.fewzFlag(),systMode);

  // Prepare output directory
  gSystem->mkdir(inpMgr.nTupleDir(systMode),true);

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //EventCounter_t ec;

  //
  // Access samples and fill histograms
  //
  const int scBrIsActive=1;
  AccessOrigNtuples_t accessInfo(scBrIsActive);
  
  // Data structures to store info from TTrees
  //mithep::TEventInfo *info    = new mithep::TEventInfo();
  //mithep::TGenInfo *gen       = new mithep::TGenInfo();
  //TClonesArray *dielectronArr = new TClonesArray("mithep::TDielectron");
  //TClonesArray *pvArr         = new TClonesArray("mithep::TVertex");


  // Reduce output file
  ZeeData_t *data=new ZeeData_t();
  ZeeData_t::Class()->IgnoreTObjectStreamer();
  
  //
  // loop over samples
  //
  if (DYTools::processData(runMode)) {
  for(UInt_t isam=0; isam<inpMgr.sampleCount(); isam++) {        

    // Configure the object for trigger matching	

    const CSample_t* samp = inpMgr.sampleInfo(isam);
    // prepare event selector
    bool isData = (samp->name=="data") ? true : false;
    //requiredTriggers.actOnData(isData);
    evtSelector.setTriggerActsOnData(isData);
    if (!evtSelector.setEScaleCorrectionType(isData,systMode)) 
      return retCodeError;

    //
    // Prepare ntuple file name
    //
    TString outName = inpMgr.nTupleFullFileName(isam,systMode);
    if ((systMode!=DYTools::NO_SYST) && 
	(systMode!=DYTools::LOWER_ET_CUT) &&
	(isam!=0)) {
      std::cout << "... systMode=<" << SystematicsStudyName(systMode) 
		<< ">, skipping the non-data files\n";
      break;
    }

    //
    // Set up output ntuple file for the sample
    //
    std::cout << "recreating file <" << outName << ">\n";
    TFile outFile(outName,"RECREATE");
    TTree *outTree = new TTree("Events","Events");
    outTree->Branch("Events","ZeeData_t",&data);

    //
    // loop through files
    //
    EventCounterExt_t ecSample(samp->name);
    for(UInt_t ifile=0; ifile<samp->size(); ifile++) {
      cout << "Processing " << samp->getFName(ifile) << "... " << std::endl;
      TFile *infile=new TFile(samp->getFName(ifile),"read");
      if (!infile || !infile->IsOpen()) {
	TString skimName=inpMgr.convertSkim2Ntuple(samp->getFName(ifile));
	std::cout <<  "  .. failed. Trying <" << skimName << ">" << std::endl;
	infile= new TFile(skimName,"read");
      }
      assert(infile->IsOpen());
    
      if (!accessInfo.prepareJson(samp->getJsonFName(ifile))) {
	std::cout << "\nfailed at json file <" << samp->getJsonFName(ifile) << ">\n";
	return retCodeError;
      }
      
      // Get the TTree
      int isSignalMC=(samp->name=="zee") ? 1:0;
      // Gen branch is activated only for signal MC:
      //int setupGenBranch=isSignalMC;
      // Photon branch activation is set in the constructor
      if (!accessInfo.setTree(*infile,"Events", isSignalMC)) return retCodeError;
      
      // Determine maximum number of events to consider
      // *** CASES ***
      // <> lumi < 0                             => use all events in the sample
      // <> xsec = 0                             => for data (use all events)
      // <> lumi > 0, xsec > 0, doWeight = true  => use all events and scale to lumi
      // <> lumi > 0, xsec > 0, doWeight = false => compute expected number of events
      ULong_t maxEvents = accessInfo.getEntries();
      const double extraWeightFactor=1.0;
      if (! evWeight.setWeight_and_adjustMaxEvents(maxEvents, inpMgr.totalLumi(), samp->getXsec(ifile), 
						   extraWeightFactor, inpMgr.selectEventsFlag())) {
	std::cout << "adjustMaxEvents failed\n";
	return retCodeError;
      }


      // loop through events
      EventCounterExt_t ec(Form("%s_file%d",samp->name.Data(),ifile));
      ec.setIgnoreScale(0); // 1 - count events, 0 - take weight in account
      // adjust the scale in the counter
      // if FEWZ weight should be considered, use evWeight.totalWeight() after
      // the FEWZ weight has been identified (see a line below)
      ec.setScale(evWeight.baseWeight());

      std::cout << "numEntries = " << accessInfo.getEntriesFast() 
		<< ", " << maxEvents << " events will be used" << std::endl;

      for(ULong_t ientry=0; ientry<maxEvents; ientry++) {
	ec.numEvents_inc();
	//if (DYTools::isDebugMode(runMode) && (ientry>10000)) break; // debug option
	if (DYTools::isDebugMode(runMode) && (ientry>100)) break; // debug option
	printProgress(100000," ientry=",ientry,maxEvents);
	
	if( isSignalMC ) {
	  // Load generator level info
	  accessInfo.GetGen(ientry);
	  // If the Z->ll leptons are not electrons, discard this event.
	  // This is needed for signal MC samples such as Madgraph Z->ll
	  // where all 3 lepton flavors are possible
	  if (!accessInfo.genLeptonsAreElectrons()) continue;
	}

	accessInfo.GetInfoEntry(ientry);
	if (!accessInfo.eventInJson()) continue;   // not certified run? Skip to next event...

	//if (!accessInfo.eventTriggerOk(requiredTriggers)) continue; // no trigger accept? Skip to next event...	
	if (!evtSelector.eventTriggerOk(accessInfo)) continue; // no trigger accept? Skip to next event...	
	ec.numEventsPassedEvtTrigger_inc();

	// for data we need the PVArr
	if (isData) accessInfo.GetPVs(ientry);

	// init fewz weight for the event
	// access info should return a pointer for a non-Zee sample
	evWeight.setFewzWeight(accessInfo.genPtr());

	// if PU weight is needed for statistics, uncomment it here
	//if (inpMgr.puReweightFlag() && !isData) {
	//  int nPV=accessInfo.getNPV(isData);
	//  evWeight.setPUWeight( nPV );
	//}

	// adjust the scale in the counter to include FEWZ 
	// (and possibly PU) weight
	//ec.setScale(evWeight.totalWeight());
	
	// load dielectron array
	accessInfo.GetDielectrons(ientry);

	// loop through dielectrons
	UInt_t dielCount=0;
	int photonBrLoaded=0;

	for(Int_t i=0; i<accessInfo.dielectronCount(); i++) {
	  mithep::TDielectron *dielectron = accessInfo.editDielectronPtr(i);
	  ec.numDielectrons_inc();
	  
	  // escale may modify dielectron!
	  if (!evtSelector.testDielectron(dielectron,accessInfo.evtInfoPtr(),&ec)) continue;

	  //hMass2v[isam]->Fill(dielectron->mass,weight);
	  //hMass3v[isam]->Fill(dielectron->mass,weight);


          /******** We have a Z candidate! HURRAY! ********/
	  
	  dielCount++;
	  ec.numDielectronsPass_inc();
	  if (ec.numDielectronsOkSameSign_inc(dielectron->q_1,dielectron->q_2)) {
	    // same sign event
	  }

	  // event printout
	  /*
          if((isam==0) && evtfile.is_open())
            eventDump(evtfile, dielectron, info->runNum, info->lumiSec, info->evtNum, 
		      leadingTriggerObjectBit, trailingTriggerObjectBit);
	  */
	  
	  //
	  // Fill histograms
	  // 
	  //hMassv[isam]->Fill(dielectron->mass,weight);

	  //if ((isam==0) && (ientry<100)) std::cout << "ientry=" << ientry << ", weight: " << evWeight << "\n";

	  // fill ntuple data
	  // For PU reweighting following the
	  // Hildreth scheme, we need the simulation level number of PU interactions

	  int nPVs=accessInfo.getNPV(isData);
	  if (inpMgr.puReweightFlag() && !isData) {
	    evWeight.setPUWeight( nPVs );
	  }
	  else evWeight.setPUWeightValue(1.0); // redundant

	  if (!photonBrLoaded) {
	    photonBrLoaded=1;
	    accessInfo.GetPhotons(ientry);
	  }
	  int phoIdx_1= accessInfo.locateSCID(dielectron->scID_1);
	  int phoIdx_2= accessInfo.locateSCID(dielectron->scID_2);
	  const mithep::TPhoton *sc_1=(phoIdx_1>=0) ? accessInfo.photonPtr(phoIdx_1) : NULL;
	  const mithep::TPhoton *sc_2=(phoIdx_2>=0) ? accessInfo.photonPtr(phoIdx_2) : NULL;

	  data->Assign(accessInfo.evtInfoPtr(), dielectron,
		       sc_1, sc_2,
		       nPVs,
		       evWeight
		       );

	  if (1 && isDebugMode(runMode)) {
	    //std::cout << "ientry=" << ientry << ", nPU=" << accessInfo.getNPV(isData) << ", mass=" << dielectron->mass << ", weight=" << evWeight.baseWeight() << ", weightSave=" << evWeight.totalWeight() << "\n";
	    std::cout << "ientry=" << ientry << ", nPU=" << accessInfo.getNPV(isData) << ", mass=" << dielectron->mass << ", weightSave=" << evWeight.totalWeight() << "\n";
	  }
	  outTree->Fill();
	}

	if (dielCount>1) ec.numMultiDielectronsOk_inc();
      }
      ec.print();  // print info about file
      ecSample.add(ec); // accumulate event counts

      infile->Close();
      delete infile;
      //eventTree=0;
    }
    ecSample.print(); // print info about sample

    std::cout << "next file" << std::endl;
    outFile.Write();
    delete outTree;
    outFile.Close();        
    //delete outFile;

    evtSelector.printCounts();
  }
  }
  else {
    std::cout << " DYTools::LOAD_DATA detected\n";
  }

  //if (evtfile.is_open()) evtfile.close();

  if (systMode!=DYTools::NO_SYST) {
    std::cout << "\n\tsystMode=" << SystematicsStudyName(systMode) << ". Terminating the macro\n";
    return retCodeOk;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================

  /*
  TString outFNameHistos = outputDir + TString("/selectEvents") + DYTools::analysisTag_USER + TString("-plots.root");
  TFile *outFileHistos = new TFile(outFNameHistos,"RECREATE");
  outFileHistos->cd();
  for(UInt_t isam=0; isam<samplev.size(); isam++) hNGoodPVv[isam]->Write();

  TString canvasName="selectEvents" + DYTools::study2Dstr;
  TCanvas *c = MakeCanvas(canvasName,canvasName,canw,canh);

  printf("Make plots\n");fflush(stdout);
  // string buffers
  char ylabel[100];   // y-axis label
  char lumitext[50];
  if(lumi>0) {
    if(lumi<1) { sprintf(lumitext,"#int#font[12]{L}dt = %.0f nb^{-1}",1000.*lumi); }
    else       { sprintf(lumitext,"#int#font[12]{L}dt = %.3g pb^{-1}",lumi); }
  }
      
  // scale factor for yield in MC to equal yield in data
  Double_t mcscale=1;
  if(hasData) {
    Double_t numer = nSelv[0];
    Double_t denom = 0;
    for(UInt_t isam=1; isam<samplev.size(); isam++)
      denom += nSelv[isam];
    mcscale = (denom>0) ? numer/denom : 1.0;
  }
  
  printf("Plot dielectron mass\n");fflush(stdout);
  // dielectron mass
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMassv[0]->GetBinWidth(1));
  CPlot plotMass("mass","","m(e^{+}e^{-}) [GeV/c^{2}]",ylabel);
  if(hasData) { plotMass.AddHist1D(hMassv[0],samplev[0]->label,"E"); }
  printf("  debug1\n");fflush(stdout);
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    hMassv[isam]->Scale(mcscale);
    plotMass.AddToStack(hMassv[isam],samplev[isam]->label,samplev[isam]->color);
  }
  printf("  debug2\n");fflush(stdout);
  if(samplev.size()>5)
    plotMass.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotMass.TransLegend(0.1,0);
  if(lumi>0) plotMass.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotMass.Draw(c,kFALSE,format);
  
  plotMass.SetName("masslog");
  plotMass.SetLogy();
  if( plotMass.GetStack() != NULL)
    plotMass.SetYRange((1e-4)*(plotMass.GetStack()->GetMaximum()),10.*(plotMass.GetStack()->GetMaximum()));  
  plotMass.Draw(c,kFALSE,format);

  c->Write();
  //outFileHistos->Close();
  */

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //============================================================================================================== 
  /*
  ofstream txtfile;
  char txtfname[100];    
  sprintf(txtfname,"%s/summary.txt",outputDir.Data());
  txtfile.open(txtfname);
  assert(txtfile.is_open());
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;
  txtfile << endl;

  txtfile << "  L_int = " << lumi << "/pb" << endl;
  txtfile << endl;
  
  if(hasData) {
    txtfile << "   Data: " << setprecision(1) << fixed << nProcessedEvents << " events processed!" << endl;
    txtfile << "         " << setprecision(1) << fixed << nSelv[0] << " Z events!" << endl;
    txtfile << "         " << setprecision(1) << fixed << nPosSSv[0] << " SS (+) events!" << endl;
    txtfile << "         " << setprecision(1) << fixed << nNegSSv[0] << " SS (-) events!" << endl;
    for(UInt_t ifile=0; ifile<samplev[0]->fnamev.size(); ifile++)
      txtfile << "     " << samplev[0]->fnamev[ifile] << endl;
      txtfile << endl;
  } 
  
  if(samplev.size()>1) {
    txtfile << "   MC:" << endl;
    for(UInt_t isam=1; isam<samplev.size(); isam++) {      
      for(UInt_t ifile=0; ifile<samplev[isam]->fnamev.size(); ifile++) {
        if(ifile==0) {
          txtfile << setw(10) << snamev[isam];
          txtfile << setw(10) << setprecision(2) << fixed << nSelv[isam] << " +/- ";
          txtfile << setw(5) << setprecision(2) << fixed << sqrt(nSelVarv[isam]);
          txtfile << "   " << "SS (+) = " << setw(5) << setprecision(3) << nPosSSv[isam];
	  txtfile << "   " << "SS (-) = " << setw(5) << setprecision(3) << nNegSSv[isam];
          txtfile << "   " << samplev[isam]->fnamev[ifile] << endl;
        } else {
          txtfile << setw(48) << "" << "   " << samplev[isam]->fnamev[ifile] << endl;
        }
      }
      txtfile << endl;
    }
  }
  txtfile.close();
  cout << "file <" << txtfname << "> created\n";

  cout << endl;
  cout << " <> Output saved in " << outputDir << "/" << endl;
  cout << endl;
  */
        
  gBenchmark->Show("selectEvents");
  std::cout << "exit" << std::endl;
  return retCodeOk;
} 


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------

/*
void eventDump(ofstream &ofs, const mithep::TDielectron *dielectron, 
               const UInt_t runNum, const UInt_t lumiSec, const UInt_t evtNum, 
	       const UInt_t triggerObj1, const UInt_t triggerObj2)
{
  ofs << endl;
  ofs << "Run:" << runNum;
  ofs << "  Lumi:" << lumiSec;
  ofs << "  Event:" << evtNum;
  ofs << "  mass: " << dielectron->mass;
  ofs << "  pt: " << dielectron->pt << endl;
  
  ofs << "----------+-----------+-----------+-----------+-----------+-----------+-------------+------------+------------+-----------+------" << endl;
  ofs << "  SC ET   |  SC eta   |   SC phi  | trkiso/pt | emiso/pt  | hadiso/pt | sigiEtaiEta |    deta    |    dphi    |    H/E    | HLT" << endl;
  ofs << "----------+-----------+-----------+-----------+-----------+-----------+-------------+------------+------------+-----------+------" << endl;
      
  ofs << setw(9) << dielectron->scEt_1 << " |";
  ofs << setw(10) << dielectron->scEta_1 << " |";
  ofs << setw(10) << dielectron->scPhi_1 << " |";
  ofs << setw(10) << dielectron->trkIso03_1/dielectron->pt_1 << " |";
  ofs << setw(10) << dielectron->emIso03_1/dielectron->pt_1 << " |";
  ofs << setw(10) << dielectron->hadIso03_1/dielectron->pt_1 << " |";
  ofs << setw(12) << dielectron->sigiEtaiEta_1 << " |";
  ofs << setw(12) << dielectron->deltaEtaIn_1 << "|";
  ofs << setw(12) << dielectron->deltaPhiIn_1 << "|";
  ofs << setw(10) << dielectron->HoverE_1 << " |";
  if(dielectron->hltMatchBits_1 & triggerObj1)
    ofs << " LEAD" << endl; 
  else if(dielectron->hltMatchBits_1 & triggerObj2)
    ofs << " TRAIL" << endl;
  else
    ofs << " NOMAT" << endl;
    
  ofs << setw(9) << dielectron->scEt_2 << " |";
  ofs << setw(10) << dielectron->scEta_2 << " |";
  ofs << setw(10) << dielectron->scPhi_2 << " |";
  ofs << setw(10) << dielectron->trkIso03_2/dielectron->pt_2 << " |";
  ofs << setw(10) << dielectron->emIso03_2/dielectron->pt_2 << " |";
  ofs << setw(10) << dielectron->hadIso03_2/dielectron->pt_2 << " |";
  ofs << setw(12) << dielectron->sigiEtaiEta_2 << " |";
  ofs << setw(12) << dielectron->deltaEtaIn_2 << "|";
  ofs << setw(12) << dielectron->deltaPhiIn_2 << "|";
  ofs << setw(10) << dielectron->HoverE_2 << " |";
  if(dielectron->hltMatchBits_2 & triggerObj1)
    ofs << " LEAD" << endl; 
  else if(dielectron->hltMatchBits_2 & triggerObj2)
    ofs << " TRAIL" << endl;
  else
    ofs << " NOMAT" << endl;
}    

*/
