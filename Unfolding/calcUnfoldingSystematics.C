#include <TFile.h>
#include "../Include/DYTools.hh"
#include <TRandom.h>
#include "../Include/HistoPair.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/EventSelector.hh"
#include "../Include/UnfoldingMatrix.h"
#include "../Include/ComparisonPlot.hh"


int unfoldDET_local(const HistoPair2D_t &iniYields, 
		    HistoPair2D_t &finYields,
		    const InputFileMgr_t &inpMgr,
		    DYTools::TSystematicsStudy_t systMode,
		    UnfoldingMatrix_t **ptrDetResponse=NULL);

//TCanvas* plotHistos(


struct TInfoLoc_t {
public:
  const HistoPair2D_t* hp;
  TString slabel, label;
public:
  TInfoLoc_t(const HistoPair2D_t *setHP, TString setShortLabel, 
	     TString setLongLabel) :
    hp(setHP), slabel(setShortLabel), label(setLongLabel)
  {}
};


// --------------------------------------------------
//  Main code
// --------------------------------------------------

int calcUnfoldingSystematics(int analysisIs2D,
			     TString conf, int debug,
	     int needInfo=0, std::vector<const TInfoLoc_t*> *infoV=NULL) {

  // check whether it is a calculation
  if (conf.Contains("_DebugRun_")) {
    std::cout << "calcUnfoldingSystematics: _DebugRun_ detected. Terminating the script\n";
    return retCodeOk;
  }

  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return retCodeError;
  }

  //----------------------------------------
  // Settings 
  //========================================
  
  int calcFSR=1;
  int calcPU=1;
  int calcRnd=1;

  DYTools::TRunMode_t runMode=DebugInt2RunMode(debug);
  int doCalc=processData(runMode);
  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;
  TString correctionKind="unfYields";

  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(conf)) return retCodeError;
  // plotDetResponse uses escale!
  // no energy correction for this evaluation
  inpMgr.clearEnergyScaleTag();

  // Construct eventSelector, update mgr and plot directory
  TString extraTag;
  extraTag="";
  EventSelector_t evtSelector(inpMgr,runMode,systMode,
			      extraTag, "", EventSelector::_selectDefault);
  evtSelector.setTriggerActsOnData(false);

  const int seedMin=inpMgr.userKeyValueAsInt("SEEDMIN");
  const int seedMax=inpMgr.userKeyValueAsInt("SEEDMAX");
  int seedFirst=0;
  int seedLast=(calcRnd) ? (seedMax-seedMin+1) : -1;
  std::cout << "seedMin..seedMax=" << seedMin << ".." << seedMax << "; seedLast=" << seedLast << "\n";

  //return retCodeOk;

  //TString ioFName=Form("debug_unfSyst_%s.root",DYTools::analysisTag.Data());
  TString ioFName=inpMgr.correctionSystFullFileName(correctionKind,DYTools::NO_SYST,0);
  TFile *ioFile=NULL;

  if (doCalc) ioFile= new TFile(ioFName,"recreate");
  else ioFile= new TFile(ioFName,"read");

  UnfoldingMatrix_t *DetResponse=NULL;

  //---------------------------------------
  // Main analysis code 
  //=======================================


  std::cout << mainpart;
  int res=1;

  HistoPair2D_t signalYieldMCbkg("signalYieldMCbkg");
  HistoPair2D_t signalYieldDDbkg("signalYieldDDbkg");
  const int ignoreDebugRunFlag=0;
  TString yieldsFileName= inpMgr.signalYieldFullFileName(systMode,ignoreDebugRunFlag);
  int loadYields=1;
  if (loadYields && doCalc) {
    TFile fileYields(yieldsFileName,"read");
    if (!fileYields.IsOpen()) {
      std::cout << "Failed to load yields from <" << fileYields.GetName() << ">\n";
      return retCodeError;
    }
    if (res) res=signalYieldMCbkg.Read();
    if (res) res=signalYieldDDbkg.Read();
    fileYields.Close();
  }

  int useMCbkg=1;
  HistoPair2D_t *signal=(useMCbkg) ? &signalYieldMCbkg : &signalYieldDDbkg;

  if (doCalc) { 
    if (res && ioFile) { 
      ioFile->cd(); 
      signalYieldMCbkg.Write(*ioFile,"rawInput");
      signalYieldDDbkg.Write(*ioFile,"rawInput");
      ioFile->cd();
      signal->Write(); 
    } 
  }
  else { res=signal->Read(); }

  HERE("res=%d",res);

/////////////////////////////////
//calculate central values
/////////////////////////////////

  HistoPair2D_t unfYields("hunfYields");
  if (doCalc) {
    if (res) res=unfoldDET_local(*signal,unfYields,inpMgr,DYTools::NO_SYST,
				 (calcRnd) ? &DetResponse : NULL);

    if (res && ioFile) { ioFile->cd(); res=unfYields.Write(); }
  }
  else {
    if (res) res=unfYields.Read();
  }

  HERE("res=%d after central values",res);

/////////////////////////////////
//calculate smearing systematics 
/////////////////////////////////

  TH2D *h2RndSyst=NULL;

  if (calcRnd) {
    std::vector<TH2D*> rndUnfVec;
    rndUnfVec.reserve(seedLast+1);
    rndUnfVec.push_back(Clone(unfYields.histo(),"hNonRndUnfYield","hNonRndUnfYield"));
    h2RndSyst=Clone(unfYields.histo(),"h2RndSyst","h2RndSyst");
    
    int nseeds=0;

    TString UrndName="detResponse_rnd";
    UnfoldingMatrix_t Urnd(UnfoldingMatrix::_cDET_Response,UrndName);
    HistoPair2D_t rndUnfYield("rndYield");
    TH2D *hSum = Clone(unfYields.histo(),"hRndSystSum","hRndSystSum");
    hSum->Reset();
    std::vector<TString> labelsV;

    //seedLast=seedFirst+1;
    for(int i=seedFirst; i<=seedLast; i++){
      nseeds++;
      gRandom->SetSeed(i+seedMin);
      if (!Urnd.randomizeMigrationMatrix(*DetResponse)) {
	std::cout << "randomization failed for iseed=" << i << "\n";
	return retCodeError;
      }
      if (!rndUnfYield.unfold(*Urnd.getDetInvResponse(), *signal)) {
	std::cout << "unfold failed for iseed=" << i << "\n";
	return retCodeError;
      }
      TString newName=Form("hRndUnf_seedNo%d",i);
      TH2D* hres=Clone(rndUnfYield.histo(),newName,newName);
      labelsV.push_back(Form("seedNo%d",i));
      rndUnfVec.push_back(hres);
      
      TH2D *hAdd=Clone(rndUnfYield.histo(),"hAdd","hAdd");
      if (!setErrorAsContent(hAdd,rndUnfYield.histo())) {
	std::cout << "error at iseed=" << i << "\n";
	return retCodeError;
      }
      //printHisto(hAdd);
      
      hSum->Add(hAdd,1.);
      delete hAdd;
    }
    //printHisto(hSum);

    hSum->Scale(1/double(nseeds));
    // Final calculation of the mean and RMS for randomization
    for (int ibin=1; ibin<=hSum->GetNbinsX(); ++ibin) {
      for (int jbin=1; jbin<=hSum->GetNbinsY(); ++jbin) {
	double avg=hSum->GetBinContent(ibin,jbin);
	double sqrAvg=hSum->GetBinError(ibin,jbin);
	double varSqr=sqrAvg*sqrAvg*nseeds - avg*avg;
	h2RndSyst->SetBinError(ibin,jbin, sqrt(varSqr));
      }
    }

    TCanvas *cx=plotProfiles("cRnd",rndUnfVec,labelsV,NULL,1);
    //cx->cd(1);
    //removeError(hSum);
    //hSum->Draw("same");
    SaveCanvas(cx,"fig-chkRnd");

    printHisto(hSum);
    printHisto(h2RndSyst);

    if (res && ioFile) {
      res=saveVec(*ioFile,rndUnfVec,"rndUnf");
      if (res) res=saveHisto(*ioFile,hSum,"rndUnfSum");
      if (res) res=saveHisto(*ioFile,h2RndSyst,"");
    }
  }

/////////////////////////////////
//calculate Fsr systematics 
/////////////////////////////////

  HistoPair2D_t unfYieldsFsr5plus("FSR5plus");
  HistoPair2D_t unfYieldsFsr5minus("FSR5minus");

  if (calcFSR) {
    if (doCalc) {
      if (res) res=unfoldDET_local(*signal,unfYieldsFsr5plus,inpMgr,DYTools::FSR_5plus);
      if (res) res=unfoldDET_local(*signal,unfYieldsFsr5minus,inpMgr,DYTools::FSR_5minus);
      if (res && ioFile) {
	unfYieldsFsr5plus.print();
	ioFile->cd();
	if (res) res= unfYieldsFsr5plus.Write(*ioFile,"details");
	if (res) res= unfYieldsFsr5minus.Write(*ioFile,"details");
      }
    }
    else {
      if (res) res= unfYieldsFsr5plus.Read(*ioFile,"details");
      if (res) res= unfYieldsFsr5minus.Read(*ioFile,"details");
    }
  }

  HERE("res=%d after load FSR samples",res);

  TString hFSRsystName=Form("h%s_FSRsyst",correctionKind.Data());
  TH2D* h2FSRsyst= getRelDifference(unfYields.histo(),hFSRsystName,1,unfYieldsFsr5plus.histo(),unfYieldsFsr5minus.histo());

  if (ioFile && doCalc) {
    if (res) { if (!h2FSRsyst->Write()) res=0; }
  }

  HERE("res=%d after FSR syst",res);

/////////////////////////////////
//calculate pile-up systematics 
/////////////////////////////////

  HistoPair2D_t unfYieldsPU5plus("PU5plus");
  HistoPair2D_t unfYieldsPU5minus("PU5minus");

  if (calcPU) {
    if (doCalc) {
      std::cout << "res=" << res << "\n";
      if (res) res=unfoldDET_local(*signal,unfYieldsPU5plus,inpMgr,DYTools::PILEUP_5plus);
      unfYieldsPU5plus.print();
      if (res) res=unfoldDET_local(*signal,unfYieldsPU5minus,inpMgr,DYTools::PILEUP_5minus);
      if (res && ioFile) {
	if (res) res= unfYieldsPU5plus.Write(*ioFile,"details");
	if (res) res= unfYieldsPU5minus.Write(*ioFile,"details");
      }
    }
    else {
      if (res) res= unfYieldsPU5plus.Read(*ioFile,"details");
      if (res) res= unfYieldsPU5minus.Read(*ioFile,"details");
    }
  }
    
  TString hPUsystName=Form("h%s_PileUpSyst",correctionKind.Data());
  TH2D* h2PUsyst= getRelDifference(unfYields.histo(),hPUsystName,1,
				   unfYieldsPU5plus.histo(),unfYieldsPU5minus.histo());

  if (ioFile && doCalc) {
    if (res) { if (!h2PUsyst->Write()) res=0; }
  }

  HERE("res=%d after pile-up syst",res);

/////////////////////////////////
//calculate total systematics 
/////////////////////////////////
  /*
    TVectorD unfoldingSystPercent(nUnfoldingBins); 
    for(int idx = 0; idx < nUnfoldingBins; idx++) {
      unfoldingSystPercent[idx] = 
	sqrt(
	     unfoldingSystPercentSmear[idx]*unfoldingSystPercentSmear[idx] +
	     unfoldingSystPercentFsr[idx]*unfoldingSystPercentFsr[idx] 
	     );
    
  }
  */

//printing out to the screen

  /*
   printf("mass     rapidity   mean-unfolded   RMS-unfolded   rel-error     rel-err-percent-Smear     rel-err-percent-Fsr      rel-err-percent-total \n");
   for(int i=0, idx=0; i<DYTools::nMassBins; i++){
    double *rapidityBinLimits=DYTools::getYBinLimits(i);
    for (int yi=0; yi<DYTools::nYBins[i]; ++yi, ++idx) {
      printf("%4.0f-%4.0f  %4.2lf-%4.2lf  %7.1f      %7.1f          %6.4f             %6.1f                 %6.3f                       %6.1f\n",
	     DYTools::massBinLimits[i],DYTools::massBinLimits[i+1],
	     rapidityBinLimits[yi],rapidityBinLimits[yi+1],
	     unfoldedYieldsMean[idx], unfoldedYieldsRMS[idx],
	     unfoldedYieldsRMS[idx]/unfoldedYieldsMean[idx],
	     unfoldingSystPercentSmear[idx], unfoldingSystPercentFsr[idx], 
	     unfoldingSystPercent[idx]);
      //unfoldedYieldsRMS[idx]*100.0/unfoldedYieldsMean[idx]);
    }
  }
  


  // Store constants in the file
  TString outputDir(TString("../root_files/systematics/")+tagDirConstants);
  gSystem->mkdir(outputDir,kTRUE);
  TString unfoldingSystFileName(outputDir+TString("/unfolding_systematics") 
				+ DYTools::analysisTag + TString(".root"));

  TFile fa(unfoldingSystFileName,"recreate");
  unfolding::writeBinningArrays(fa);
  unfoldedYieldsMean.Write("unfoldedYieldsMeanFI");
  unfoldedYieldsRMS.Write("unfoldedYieldsRMSFI");
  unfoldingSystPercentSmear.Write("unfoldingSystPercentSmearFI");
  unfoldingSystPercentFsr.Write("unfoldingSystPercentFsrFI");
  unfoldingSystPercent.Write("unfoldingSystPercentFI");
  fa.Close();

  */

  HERE("res=%d at finishing up",res);

  if (ioFile) { 
    if (doCalc) writeBinningArrays(*ioFile);
    ioFile->Close();
    if (doCalc) {
      std::cout << "created a debug file: <" << ioFile->GetName() << ">\n";
    }
    else {
      std::cout << "used a debug file: <" << ioFile->GetName() << ">\n";
    }
    delete ioFile; 
  }


  if (res && infoV) {
    infoV->push_back(new TInfoLoc_t(signal,"signal",(useMCbkg) ? "signal (MC bkg)" : "signal (DD bkg)"));
    infoV->push_back(new TInfoLoc_t(&unfYields,"unfYields","unfolded signal"));
    if (needInfo==1) {
      infoV->push_back(new TInfoLoc_t(&unfYieldsFsr5plus,"FSR5plus","unf.signal FSR5plus"));
      infoV->push_back(new TInfoLoc_t(&unfYieldsFsr5minus,"FSR5minus","unf.signal FSR5minus"));
      //infoV->push_back(new TInfoLoc_t( unfYields_dev_Fsr5plus,"devFSR5plus","diff.FSR5plus"));
      //infoV->push_back(new TInfoLoc_t( unfYields_dev_Fsr5minus,"devFSR5minus","diff.FSR5minus"));
      //infoV->push_back(new TInfoLoc_t( unfYields_diff_FSR,"diffFSR","diff.FSR"));
      infoV->push_back(new TInfoLoc_t(&unfYieldsPU5plus,"PU5plus","unf.signal PU5plus"));
      infoV->push_back(new TInfoLoc_t(&unfYieldsPU5minus,"PU5minus","unf.signal PU5minus"));
      //infoV->push_back(new TInfoLoc_t( unfYields_dev_PU5plus,"devPU5plus","diff.PU5plus"));
      //infoV->push_back(new TInfoLoc_t( unfYields_dev_PU5minus,"devPU5minus","diff.PU5minus"));
      //infoV->push_back(new TInfoLoc_t( unfYields_diff_PU,"diffPU","diff.PU"));
    }
  }

  return (res) ? retCodeOk : retCodeError;
}

//-----------------------------------------------------------------
// Unfold
//-----------------------------------------------------------------

int unfoldDET_local(const HistoPair2D_t &iniYields, 
		    HistoPair2D_t &finYields,
		    const InputFileMgr_t &inpMgr,
		    DYTools::TSystematicsStudy_t systMode,
		    UnfoldingMatrix_t **ptrDetResponse) 
{

  TString wStr;
  const TString u="_";
  switch (systMode) {
  case DYTools::FSR_5plus: wStr="_105"; break;
  case DYTools::FSR_5minus: wStr="_095"; break;
  case DYTools::PILEUP_5plus: wStr="_PU5plus"; break;
  case DYTools::PILEUP_5minus: wStr="_PU5minus"; break;
  default: ;
  }
  TString umName=TString("detResponse") + wStr;

  DYTools::TSystematicsStudy_t dirSystMode=DYTools::NO_SYST;
  switch(systMode) {
  case DYTools::FSR_5plus:
  case DYTools::FSR_5minus:
    dirSystMode= DYTools::FSR_STUDY;
    break;
  case DYTools::PILEUP_5plus:
  case DYTools::PILEUP_5minus:
    dirSystMode= DYTools::PU_STUDY;
    break;
  default: ;
  }
  
  UnfoldingMatrix_t detResponse(UnfoldingMatrix::_cDET_Response,umName);
  TString constDir=inpMgr.constDir(dirSystMode,0);
  TString fnameTag=UnfoldingMatrix_t::generateFNameTag(dirSystMode);
  int res=detResponse.autoLoadFromFile(constDir,fnameTag);
  //std::cout << "autoLoad of " << umName << " from " << constDir << " with fnameTag=" << fnameTag << " has res=" << res << "\n";
  //std::cout << "iniYields"; iniYields.print();
  //std::cout << "detResponse.getDetInvResponse(): "; detResponse.getDetInvResponse()->Print();
  if (!res) res=-1;
  if ((res==1) && ptrDetResponse) {
    TString setName=TString("detResponse_") + UnfoldingMatrix_t::generateFNameTag(systMode);
    (*ptrDetResponse)=new UnfoldingMatrix_t(detResponse,setName);
  }

  if (res==1) res=finYields.unfold(*detResponse.getDetInvResponse(), iniYields);

  if (res!=1) {
    std::cout << "ERROR in unfoldDET_local:\n  ";
    if (res==-1) std::cout << "failed to load unfolding matrix\n";
    else std::cout << "error during unfolding\n";
    std::cout << "for systMode=" 
	      << SystematicsStudyName(systMode)
	      << ", derived dirSystMode="
	      << SystematicsStudyName(dirSystMode)
	      << "\n";
    return 0;
  }
  return res;
}
