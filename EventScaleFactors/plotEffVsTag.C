#include "../Include/DYTools.hh"
#include "../Include/EventSelector.hh"
#include "../Include/InputFileMgr.hh"
#include "../EventScaleFactors/effCalc.hh"
#include "../Include/colorPalettes.hh"

#include "tnpSelectEvents.hh"

#ifndef tnpStoreTag
#error  tnpStoreTag should be enabled and files should contain tagEt,tagEta fields
#endif

// ------------------------------------------------------




// ------------------------------------------------------

int plotEffVsTag(const TString configFile, const TString effTypeString="RECO", int runOnData=0 ) {

  // -------
  // setup -
  // -------
  
  if (!effTypeString.Contains("RECO") &&
      !effTypeString.Contains("ID") &&
      !effTypeString.Contains("HLT")) {
    std::cout << "calcEff: effTypeString should be \"RECO\",\"ID\" or \"HLT*\"\n";
    return retCodeError;
  }

  const DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;
  const DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN;

  TDescriptiveInfo_t tnpSection;
  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(configFile,&tnpSection) ||
      !inpMgr.KeepFirstAndLastSample()
      //|| !inpMgr.SetSkimsToNtuples()
      ) return retCodeError;

  // Construct eventSelector, update inpMgr and plot directory
  EventSelector_t evtSelector(inpMgr,runMode,systMode,
			      "", "", EventSelector::_selectDefault);
  TriggerSelection_t triggers(evtSelector.trigger());

  // Prepare output directory. No creation
  TString tagAndProbeDir=inpMgr.tnpDir(systMode,0);

 //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================

  //Double_t massLow  = 60;
  //Double_t massHigh = 120;

  DYTools::TEfficiencyKind_t effType = DetermineEfficiencyKind(effTypeString);
  printf("Efficiency type to measure: %s\n", EfficiencyKindName(effType).Data());
  // Read in the configuration file
  DYTools::TDataKind_t dataKind = (runOnData) ? DYTools::DATA : DYTools::MC;
  TString sampleTypeString = (runOnData) ? "DATA" : "MC";
  //TString calcMethodString = inpMgr.getTNP_calcMethod(tnpSection,dataKind,effType);
  TString etBinningString  = inpMgr.getTNP_etBinningString(tnpSection);
  //TString etaBinningString = inpMgr.getTNP_etaBinningString(tnpSection);
  TString etaBinningString= "DYTools::ETABINS14";
  TString dirTag= inpMgr.selectionTag();

  vector<TString> ntupleFileNames;
  vector<TString> jsonFileNames;
  inpMgr.getTNP_ntuples(tnpSection,runOnData,ntupleFileNames,jsonFileNames);
  if (1) {
    for (unsigned int i=0; i<ntupleFileNames.size(); ++i) {
      std::cout << " i=" << i << ": " << ntupleFileNames[i];
      if (jsonFileNames.size()>i) std::cout << "; json " << jsonFileNames[i];
      std::cout << "\n";
    }
  }

  //printf("Efficiency calculation method: %s\n", calcMethodString.Data());
  //int calcMethod= DetermineTnPMethod(calcMethodString);

  DYTools::TEtBinSet_t etBinning = DetermineEtBinSet(etBinningString);
  printf("SC ET binning: %s\n", EtBinSetName(etBinning).Data());

  DYTools::TEtaBinSet_t etaBinning = DetermineEtaBinSet(etaBinningString);
  printf("SC eta binning: %s\n", EtaBinSetName(etaBinning).Data());

  printf("Sample: %s\n", sampleTypeString.Data());
  int sample=DetermineDataKind(sampleTypeString);

  // Correct the trigger object
  triggers.actOnData((sample==DYTools::DATA)?true:false);
  evtSelector.setTriggerActsOnData((sample==DYTools::DATA)?true:false);

  // The label is a string that contains the fields that are passed to
  // the function below, to be used to name files with the output later.
  //TString label = getLabel(sample, effType, calcMethod, etBinning, etaBinning, triggers);

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================
  
  // Load selected events
  TString selectEventsFName=inpMgr.tnpSelectEventsFName(systMode,sampleTypeString,effTypeString,triggers.triggerSetName());
  std::cout << "selectEventsFName=<" << selectEventsFName << ">\n"; 
  TFile *selectedEventsFile = new TFile(selectEventsFName);
  if(!selectedEventsFile || !selectedEventsFile->IsOpen()) {
    std::cout << "failed to open file <" << selectEventsFName << ">\n";
    assert(0);
  }

  TTree *passTree = (TTree*)selectedEventsFile->Get("passTree");
  assert(passTree);
  TTree *failTree = (TTree*)selectedEventsFile->Get("failTree");
  assert(failTree);

  //tnpSelectEvent_t storeData;
  //storeData.setBranchAddress(passTree);
  //storeData.setBranchAddress(failTree);

  int nEtaProbe                = DYTools::getNEtaBins(etaBinning);
  const double *limitsEtaProbe = DYTools::getEtaBinLimits(etaBinning);

  TString etaCutFormat= DYTools::signedEtaBinning(etaBinning) ?
    " ( eta >= %5.3f && eta < %5.3f ) " :
    " ( abs(eta) >= %5.3f && abs(eta) < %5.3f ) ";
 
  TString tagEtaCutFormat= " ( abs(tagEta) >= %5.3f && abs(tagEta) < %5.3f ) ";
  
  //double tagEtaMin=0.;
  //double tagEtaMax=2.5;
  //int tagEtaCount=25;
  //double tagH=(tagEtaMax-tagEtaMin)/double(tagEtaCount);
  
  TString histoName="histo";
  TString histoTitle="efficiency";
  //TH2D *h2=new TH2D(histoName,histoTitle, 
  //		    nEtaProbe,limitsEtaProbe,
  //		    tagEtaCount,tagEtaMin,tagEtaMax);
  TH2D *h2=new TH2D(histoName,histoTitle,
		    nEtaProbe,limitsEtaProbe,
		    nEtaProbe,limitsEtaProbe);
  h2->SetDirectory(0);
  h2->GetXaxis()->SetTitle("probe #eta");
  h2->GetYaxis()->SetTitle("tag #eta");

  TH2D *h2pass= (TH2D*)h2->Clone("h2pass");
  h2pass->SetTitle("h2pass");
  h2pass->SetDirectory(0);
  TH2D *h2fail= (TH2D*)h2->Clone("h2fail");
  h2fail->SetTitle("h2fail");
  h2fail->SetDirectory(0);

  TCanvas *canvW=new TCanvas("canvW","pass & fail distributions",400,800);
  const int nCWrow=4;
  const int nCWcol=2;
  canvW->Divide(nCWcol,nCWrow);

  int ipad=1;
  TFile ftmp("ftmp.root","recreate");

  // calculate the efficiencies
  for (int ibin=1; ibin<=h2->GetNbinsX(); ++ibin) {
    double etaPmin=h2->GetXaxis()->GetBinLowEdge(ibin);
    double etaPmax=etaPmin + h2->GetXaxis()->GetBinWidth(ibin);
    TString etaPcut=TString::Format( etaCutFormat, etaPmin,etaPmax);
    TTree *passSubTree= passTree->CopyTree(etaPcut);
    TTree *failSubTree= failTree->CopyTree(etaPcut);
    for (int jbin=1; jbin<=h2->GetNbinsY(); ++jbin) {
      double etaTmin=h2->GetYaxis()->GetBinLowEdge(jbin);
      double etaTmax=etaTmin + h2->GetYaxis()->GetBinWidth(jbin);
      TString etaTcut=TString::Format( tagEtaCutFormat, etaTmin,etaTmax);

      TString cut = etaPcut + TString(" && ") + etaTcut;
      std::cout << "cut=" << cut << "\n";

      canvW->cd(ipad);
      //passTree->Draw("weight >> hwPass", cut);
      passSubTree->Draw("weight >> hwPass", cut);
      canvW->cd(ipad+1);
      //failTree->Draw("weight >> hwFail", cut);
      failSubTree->Draw("weight >> hwFail", cut);
      ipad=(ipad+2)%(nCWrow*nCWcol);
      TH1 *hwPass=(TH1*)gDirectory->Get("hwPass");
      TH1 *hwFail=(TH1*)gDirectory->Get("hwFail");
      double probesPassWeighted = hwPass->GetMean() * hwPass->GetEntries();
      double probesFailWeighted = hwFail->GetMean() * hwFail->GetEntries();
      double effCountWeighted, effErrLowCountWeighted, effErrHighCountWeighted;

      DYTools::calcEfficiency( probesPassWeighted, 
			       probesPassWeighted+probesFailWeighted,
			       DYTools::EFF_CLOPPER_PEARSON,
			       effCountWeighted, 
			       effErrLowCountWeighted,effErrHighCountWeighted);
      double avgErr=0.5*(effErrLowCountWeighted+effErrHighCountWeighted);
      h2->SetBinContent(ibin,jbin,effCountWeighted);
      h2->SetBinError  (ibin,jbin,avgErr);

      double probesPassErr= hwPass->GetMeanError();
      double probesFailErr= hwFail->GetMeanError();
      h2pass->SetBinContent(ibin,jbin, probesPassWeighted);
      h2pass->SetBinError  (ibin,jbin, probesPassErr);
      h2fail->SetBinContent(ibin,jbin, probesFailWeighted);
      h2fail->SetBinError  (ibin,jbin, probesFailErr);

      std::cout << "ibin,jbin=" << ibin << "," << jbin << " " << effCountWeighted << " +/- " << avgErr << "\n";
      canvW->Update();
    }
    delete passSubTree;
    delete failSubTree;
  }
  ftmp.Close();

  gStyle->SetPalette(1);
  TCanvas *cx=new TCanvas("cx","cx",1500,500);
  cx->Divide(3,1);
  AdjustFor2DplotWithHeight(cx);

  cx->cd(1);
  h2pass->Draw("COLZ");
  cx->cd(2);
  h2fail->Draw("COLZ");

  cx->cd(3);
  h2->Draw("COLZ");
  cx->Update();

  SaveCanvas(cx,TString("canvVer4_passFail_") + effTypeString,"plots_tagEta");

  return retCodeOk;
}



// ------------------------------------------------------




