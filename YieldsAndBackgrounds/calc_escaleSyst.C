// Enhanced macro of calc_escaleDiffSyst.C
// works with a greater number of systematics

#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/EventSelector.hh"
#include "../Include/HistoPair.hh"
#include "../Include/ComparisonPlot.hh"

int calc_escaleSyst(int debug, 
		    TString conf,
		    DYTools::TSystematicsStudy_t systMode=DYTools::ESCALE_STUDY) 

{

  std::cout << "debug=" << debug << " (sets 'checkTheCode')\n";
  DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN;

  if (conf==TString("default")) {
    conf="../config_files/data_vilnius8TeV_regSSD.conf.py";
  }

  {
    using namespace DYTools;
    DYTools::printExecMode(runMode,systMode);

    const int debug_print=1;
    if (!DYTools::checkSystMode(systMode,debug_print,3,
				DYTools::ESCALE_STUDY,
				DYTools::UNREGRESSED_ENERGY,
				DYTools::APPLY_ESCALE))
      return retCodeError;
  }
 
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================

  // Read from configuration file only the location of the root files
  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(conf)) 
    return retCodeError;
  // no energy correction for this evaluation
  inpMgr.clearEnergyScaleTag();

  //inpMgr.Print();


  std::vector<TString> fnames;

  // Construct eventSelector, update inpMgr and plot directory
  TString extraTag="R9";
  TString plotExtraTag;
  EventSelector_t evtSelector(inpMgr,runMode,DYTools::NO_SYST,
			      extraTag, plotExtraTag, EventSelector::_selectDefault);

  std::vector<DYTools::TSystematicsStudy_t> systModeVec;
  std::vector<TString> systModeNameV;

  switch(systMode) {
  case DYTools::ESCALE_STUDY: {
    systModeVec.push_back(DYTools::ESCALE_DIFF_0000);
    systModeVec.push_back(DYTools::ESCALE_DIFF_0005);
    systModeVec.push_back(DYTools::ESCALE_DIFF_0010);
    systModeVec.push_back(DYTools::ESCALE_DIFF_0015);
    systModeVec.push_back(DYTools::ESCALE_DIFF_0020);
    systModeVec.push_back(DYTools::NO_SYST );
  }
    break;
  case DYTools::UNREGRESSED_ENERGY: {
    systModeVec.push_back(DYTools::UNREGRESSED_ENERGY);
    systModeVec.push_back(DYTools::NO_SYST);
  }
    break;
  case DYTools::APPLY_ESCALE: {
    systModeVec.push_back(DYTools::APPLY_ESCALE);
    systModeVec.push_back(DYTools::NO_SYST);
  }
    break;
  default:
    std::cout << "calc_escaleSyst: not ready for " << SystematicsStudyName(systMode) << "\n";
    return retCodeError;
  }
  int checkTheCode=debug;
  int skipNoSystInPlot=(!checkTheCode && (systMode==DYTools::ESCALE_STUDY)) ? 1:0;

  for (unsigned int i=0; i<systModeVec.size(); i++) {
    DYTools::TSystematicsStudy_t iSystMode=(DYTools::TSystematicsStudy_t)(systModeVec[i]);
    if (checkTheCode) {
      if (i==1) iSystMode=DYTools::NO_SYST;
      else if (i>1) break;;
    }
    TString systModeName=SystematicsStudyName(iSystMode);
    systModeNameV.push_back(systModeName);
    std::cout << "systMode=" << systModeName << "\n";
    const int ignoreDebugRunFlag=0;
    TString resF= inpMgr.signalYieldFullFileName(iSystMode,ignoreDebugRunFlag);
    std::cout << "resF=" << resF << "\n";
    TFile fin(resF,"read");
    if (!fin.IsOpen()) std::cout << " ... failed to open\n";
    fin.Close();
    fnames.push_back(resF);
  }

  // However, the plots should be saved according to the systMode
  // Restore the correct value
  evtSelector.SetPlotOutDir(runMode,systMode,plotExtraTag,1);

  // inputDir+TString("/yields_bg-subtracted.root")
  const int ignoreDebugRunFlag=0;
  const int createDir=1;
  const int systematics=1;
  TString outFileName= inpMgr.signalYieldFullFileName(systMode,ignoreDebugRunFlag,createDir,systematics);
  std::cout << "generated outFileName=<" << outFileName << ">\n";

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //============================================================================================================== 

  std::vector<HistoPair2D_t*> signalYieldMCbkgV, signalYieldDDbkgV;

  for (unsigned int i=0; i<fnames.size(); ++i) {
    HistoPair2D_t sMC("signalYieldMCbkg");
    HistoPair2D_t sDD("signalYieldDDbkg");
    if (!sMC.Load(fnames[i],1,"") ||
	!sDD.Load(fnames[i],0,"")) {
      std::cout << "failed to load signalYield\n";
      return retCodeError;
    }
    //printHisto(sMC.histo());
    TString newNameMC=TString("hSigMCbkg_") + systModeNameV[i];
    TString newNameDD=TString("hSigDDbkg_") + systModeNameV[i];
    signalYieldMCbkgV.push_back(new HistoPair2D_t(newNameMC,sMC,1));
    signalYieldDDbkgV.push_back(new HistoPair2D_t(newNameDD,sDD,1));
    sMC.clear();
    sDD.clear();
    //printHisto(signalYieldMCbkgV.back()->histo());
  }

  if (0 && (systModeVec[0]==DYTools::NO_SYST)) {
    TH2D* h2RelDiffMCbkg_chk= getRelDifference(signalYieldMCbkgV,
					       "h2RelDiffMCbkg_chk", 1);
    printHisto(h2RelDiffMCbkg_chk);
    
    TH2D* h2RelDiffDDbkg_chk= getRelDifference(signalYieldDDbkgV,
					       "h2RelDiffDDbkg_chk", 1);
    printHisto(h2RelDiffDDbkg_chk);
    return retCodeStop;
  }

  TH2D* h2RelDiffMCbkg= getRelDifference(signalYieldMCbkgV,
					       "h2RelDiffMCbkg", 1);
  TH2D* h2RelDiffDDbkg= getRelDifference(signalYieldDDbkgV,
					 "h2RelDiffDDbkg", 1);

  printHisto(h2RelDiffMCbkg);
  printHisto(h2RelDiffMCbkg);
  
  if (1) {
    int res1=1;
    TFile fout(outFileName,"recreate");
    res1=fout.IsOpen();
    if (res1) {
      h2RelDiffMCbkg->Write();
      h2RelDiffDDbkg->Write();
    }
    if (res1) res1=saveVec(fout,signalYieldMCbkgV,"mcBkg");
    if (res1) res1=saveVec(fout,signalYieldDDbkgV,"ddBkg");
    writeBinningArrays(fout);
    fout.Close();
    if (!res1) {
      std::cout << "\n\tError saving file <" << fout.GetName() << ">\n\n";
    }
  }

  // Prepare the plot
  if (1) {
    
    for (int useDDBkg=0; useDDBkg<2; ++useDDBkg) {

      std::vector<HistoPair2D_t*> *src=(useDDBkg) ? &signalYieldDDbkgV : &signalYieldMCbkgV;
      std::vector<TH2D*> histosV;
      std::vector<TString> labelsV=systModeNameV;
      std::vector<int> colorsV;
      for (unsigned int i=0; i<src->size(); ++i) {
	histosV.push_back((*src)[i]->editHisto());
      }
      colorsV.push_back(kBlack);
      colorsV.push_back(kBlue);
      colorsV.push_back(kGreen+1);
      colorsV.push_back(kOrange);
      colorsV.push_back(kRed+1);
      colorsV.push_back(kViolet);

      std::vector<std::vector<TH1D*>*> hProfV;
      TString canvName=(useDDBkg) ? "cDD" : "cMC";
      int canvWidth=(DYTools::study2D==1) ? 1100 : 700;
      TCanvas *c1=new TCanvas(canvName,canvName, canvWidth,900);

      if (DYTools::study2D==1) {

	if (!createRapidityProfileVec(histosV,hProfV,labelsV)) {
	  std::cout << "failed to create profiles for useDDBkg=" << useDDBkg << "\n";
	  return retCodeError;
	}
	std::cout << "there are " << hProfV.size() << " profiles\n";
      
	for (int im=1; im<7; ++im) {
	  
	  TString mStr=Form("M_%2.0lf_%2.0lf",DYTools::massBinLimits[im],DYTools::massBinLimits[im+1]);
	  TString cpName="cp_" + mStr;
	  TString cpTitle=mStr;
	  ComparisonPlot_t *cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,cpName,cpTitle,"|y|","signal yield","ratio");
	  if (im==1) cp->Prepare6Pads(c1,1);
	  
	  for (unsigned int ih=0; ih<hProfV[im]->size(); ++ih) {
	    if (skipNoSystInPlot && (systModeVec[ih]==DYTools::NO_SYST)) continue;
	    TH1D* h=(*hProfV[im])[ih];
	    removeError1D(h);
	    if (ih==0) {
	      h->SetMarkerStyle(20);
	    }
	    cp->AddHist1D(h,labelsV[ih],"LP",colorsV[ih]);
	  }
	  cp->Draw6(c1,1,im);
	  //if (!useDDBkg) cp->ChangeLegendPos(0.2,0.,0.,0.);
	  c1->Update();
	  
	}
      }
      else {
	// 1D

	if (!createMassProfileVec(histosV,hProfV,labelsV)) {
	  std::cout << "failed to create profiles for useDDBkg=" << useDDBkg << "\n";
	  return retCodeError;
	}
	std::cout << "there are " << hProfV.size() << " profiles\n";
      
	for (int iy=0; iy<DYTools::nYBinsMax; ++iy) {
	  
	  TString yStr=Form("iy_%d",iy);
	  TString cpName=TString("cp_") + yStr;
	  TString cpTitle; //=yStr;
	  ComparisonPlot_t *cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,cpName,cpTitle,"#it{M}_{ee} [GeV]","signale yield","ratio");
	  cp->SetLogx(1);
	  if (iy==0) cp->Prepare2Pads(c1);
	  
	  for (unsigned int ih=0; ih<hProfV[iy]->size(); ++ih) {
	    if (skipNoSystInPlot && (systModeVec[ih]==DYTools::NO_SYST)) continue;
	    TH1D* h=(*hProfV[iy])[ih];
	    removeError1D(h);
	    if (ih==0) {
	      h->SetMarkerStyle(20);
	    }
	    cp->AddHist1D(h,labelsV[ih],"LP",colorsV[ih]);
	  }
	  cp->Draw(c1);
	  //if (!useDDBkg) cp->ChangeLegendPos(0.2,0.,0.,0.);
	  c1->Update();
	  
	}
      }
      
      TString fname="fig-";
      fname.Append(SystematicsStudyName(systMode));
      fname.Append((useDDBkg) ? "-DDBkg-" : "-MCBkg-");
      fname.Append(DYTools::analysisTag);
      TString outDir="plots-escaleSyst";
      SaveCanvas(c1,fname,outDir);
    }
      
  }

  return retCodeOk;
}
