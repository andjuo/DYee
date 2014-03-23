#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/EventSelector.hh"
#include "../Include/HistoPair.hh"
#include "../Include/ComparisonPlot.hh"

int calc_escaleDiffSyst(int debug, 
			TString conf,
			DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN) 

{

  std::cout << "debug=" << debug << " (ignored)\n";

  if (conf==TString("default")) {
    conf="../config_files/data_vilnius8TeV_regSSD.conf.py";
  }

  DYTools::TSystematicsStudy_t systMode=DYTools::ESCALE_STUDY;
  {
    using namespace DYTools;
    DYTools::printExecMode(runMode,systMode);
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

  const int systModeCount=5;
  const int systModeVec[systModeCount]= { DYTools::ESCALE_DIFF_0000, DYTools::ESCALE_DIFF_0005, DYTools::ESCALE_DIFF_0010, DYTools::ESCALE_DIFF_0015, DYTools::ESCALE_DIFF_0020 };
  std::vector<TString> systModeNameV;

  for (int i=0; i<systModeCount; i++) {
    DYTools::TSystematicsStudy_t iSystMode=(DYTools::TSystematicsStudy_t)(systModeVec[i]);
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
  TString outFileName= inpMgr.signalYieldFullFileName(systMode,ignoreDebugRunFlag,createDir);
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

  TH2D* h2RelDiffMCbkg= 
    getRelDifference(signalYieldMCbkgV[0]->histo(),
		     "h2RelDiffMCbkg",
		     1,
		     signalYieldMCbkgV[1]->histo(), 
		     signalYieldMCbkgV[2]->histo(),
		     signalYieldMCbkgV[3]->histo(), 
		     signalYieldMCbkgV[4]->histo());
  //printHisto(h2RelDiffMCbkg);

  TH2D* h2RelDiffDDbkg= 
    getRelDifference(signalYieldDDbkgV[0]->histo(),
		     "h2RelDiffDDbkg",
		     1,
		     signalYieldDDbkgV[1]->histo(), 
		     signalYieldDDbkgV[2]->histo(),
		     signalYieldDDbkgV[3]->histo(), 
		     signalYieldDDbkgV[4]->histo());
  //printHisto(h2RelDiffDDbkg);
  //printHisto(signalYieldDDbkgV[0]->histo());

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
	    TH1D* h=(*hProfV[im])[ih];
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
	    TH1D* h=(*hProfV[iy])[ih];
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
      
      TString fname="fig-escaleDiffSyst-";
      fname.Append((useDDBkg) ? "DDBkg-" : "MCBkg-");
      fname.Append(DYTools::analysisTag);
      SaveCanvas(c1,fname);
    }
      
  }

  return retCodeOk;
}
