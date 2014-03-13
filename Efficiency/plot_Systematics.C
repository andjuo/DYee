#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/calcCorrectionSyst.h"
#include "../Include/ComparisonPlot.hh"
#include "../Include/HistoPair.hh"
#include "../Include/MitStyleRemix.hh" // SetSideSpaces

void plot_Systematics(TString correctionKind, TString flags="111", int saveCanvas=0) {
  std::cout << "correctionKind=" << correctionKind << "\n";

  int totSyst=1;
  int FSRsyst=1;
  int PUsyst=1;
  if (flags.Length()==3) {
    totSyst=(flags[0]=='1') ? 1:0;
    FSRsyst=(flags[1]=='1') ? 1:0;
    PUsyst =(flags[2]=='1') ? 1:0;
  }

  TString foutName;
  if (correctionKind=="efficiency") {
    foutName=Form("../root_files_reg/constants/DY_j22_19712pb/efficiencySyst_%s.root",DYTools::analysisTag.Data());
  }
  else if (correctionKind=="acceptance") {
    foutName=Form("../root_files_reg/constants/DY_j22_19712pb/acceptanceSyst_%s.root",DYTools::analysisTag.Data());
  }
  else if (correctionKind=="unfYields") {
    foutName=Form("../root_files_reg/constants/DY_j22_19712pb/unfYieldsSyst_%s.root",DYTools::analysisTag.Data());
  }
  else {
    std::cout << "not ready for correctionKind=<" << correctionKind << ">\n";
    return;
  }

  std::vector<TH2D*> histos;
  std::vector<std::string> labelsV;
  if (!loadCorrectionSystFromFile(foutName,correctionKind,histos,labelsV)) return;

  TCanvas *c1=NULL, *cFSR=NULL, *cPU=NULL;

  std::vector<TH2D*> hCombiV;
  std::vector<TString> combiLabels;

  HistoPair2D_t hpBase("hpBase");
  hpBase.assign(histos[0],NULL);
  HistoPair2D_t hpFSR("hpFSR",hpBase);
  hpFSR.addSystErr(histos[1]);
  HistoPair2D_t hpPU("hpPU",hpBase);
  hpPU.addSystErr(histos[2]);
  HistoPair2D_t hpTotErr("hpTotErr",hpFSR);
  hpTotErr.addSystErr(histos[2]);
    
  hCombiV.push_back(histos[0]);
  hCombiV.push_back(hpFSR.createHistoWithFullError("hWithFSRsyst"));
  hCombiV.push_back(hpPU.createHistoWithFullError("hWithPileUpSyst"));
  hCombiV.push_back(hpTotErr.createHistoWithFullError("hWithTotSyst"));
  combiLabels.push_back(correctionKind + TString(" +stat.err"));
  combiLabels.push_back(correctionKind + TString(" +stat.err +FSR.syst."));
  combiLabels.push_back(correctionKind + TString(" +stat.err +Pile-up.syst."));
  combiLabels.push_back(correctionKind + TString(" +tot.err"));

  if (DYTools::study2D==1) {

    // -----------  2D case ----------------------

    std::vector<std::vector<TH1D*>*> hProfRawW;
    std::vector<std::vector<TH1D*>*> hProfW;

    for (int mbin=2; mbin<=7; ++mbin) {
      std::vector<TH1D*>* hProfRawV = new std::vector<TH1D*>();
      hProfRawW.push_back(hProfRawV);

      TString mStr=Form("_m%1.0f_%1.0f",DYTools::massBinLimits[mbin-1],
			DYTools::massBinLimits[mbin]);

      for (unsigned int i=0; i<histos.size(); ++i) {
	TH1D* hProf=NULL;
	if (histos[i]) {
	  TString name=Form("hProfRaw_%s_%s",labelsV[i].c_str(),mStr.Data());
	  hProf= createProfileY(histos[i],mbin,name,1,name, DYTools::nYBins[mbin-1],0.,DYTools::yRangeMax+1e-3);
	  
	  //hProf->GetXaxis()->SetMoreLogLabels();
	  //hProf->GetXaxis()->SetNoExponent();
	  hProf->SetMarkerStyle(20);
	  hProf->SetMarkerSize(0.8);
	}
	hProfRawV->push_back(hProf);
      }
    }

    for (int mbin=2; mbin<=7; ++mbin) {
      std::vector<TH1D*>* hProfV= new std::vector<TH1D*>();
      hProfW.push_back(hProfV);

      TString mStr=Form("_m%1.0f_%1.0f",DYTools::massBinLimits[mbin-1],
			DYTools::massBinLimits[mbin]);

      for (unsigned int i=0; i<hCombiV.size(); ++i) {
	TString name=Form("hProf_%s_%s",combiLabels[i].Data(),mStr.Data());
	if (!hCombiV[i]) { HERE("combiV[%u] is null",i); return ; }
	TH1D *h=createProfileY(hCombiV[i],mbin,name,1,name, DYTools::nYBins[mbin-1],0.,DYTools::yRangeMax+1e-3);
	//h->GetXaxis()->SetMoreLogLabels();
	//h->GetXaxis()->SetNoExponent();
	hProfV->push_back( h );
      }
    }

    // -----------------
    // Partial systematics
    // ------------------
    if (PUsyst+FSRsyst) {
      std::vector<ComparisonPlot_t*> cpFSRsyst,cpPUsyst;
    
      for (int isPU=0; isPU<2; ++isPU) {
	if ( isPU && !PUsyst) continue;
	if (!isPU && !FSRsyst) continue;

	TString canvName=(isPU) ? "cPU" : "cFSR";
	c1 = new TCanvas(canvName,canvName, 1100,900);
	if (isPU) cPU=c1; else cFSR=c1;

	for (int im=0; im<6; ++im) {
	  TString title=(isPU) ? "Pile-up systematics " : "FSR systematics ";
	  TString mRange=Form("M_%2.0lf_%2.0lf",DYTools::massBinLimits[im+1],DYTools::massBinLimits[im+2]);
	  std::vector<TH1D*> *hProfRawV= hProfRawW[im];
	  if (!hProfRawV) {
	    std::cout << "hProfRawV for " << title << mRange << " is null\n";
	    continue;
	  }
	  
	  ComparisonPlot_t *cp = new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,"cpSyst",title+mRange,"|y|",correctionKind,"ratio");
	  if (isPU) cpPUsyst.push_back(cp); else cpFSRsyst.push_back(cp);

	  //cp->SetLogx();
	  cp->SetPrintValues(0);
	  cp->SetYTitleSize(0.07,1.1);
	  cp->SetRatioYTitleSize(0.17,0.47);
      
	  int imin=(isPU) ? 5 : 3;
	  if ((*hProfRawV)[0]) {
	    removeError((*hProfRawV)[0]);
	    cp->AddHist1D((*hProfRawV)[0],labelsV[0],"LP",kBlack,1,0);
	  }
	  if ((*hProfRawV)[imin]) {
	    removeError((*hProfRawV)[imin]);
	    cp->AddHist1D((*hProfRawV)[imin],labelsV[imin],"LP",kBlue,1,0);
	  }
	  if ((*hProfRawV)[imin+1]) {
	    removeError((*hProfRawV)[imin+1]);
	    cp->AddHist1D((*hProfRawV)[imin+1],labelsV[imin+1],"LP",kGreen+1,1,0);
	  }
	  
	  if (im==0) cp->Prepare6Pads(c1, 1);
	  cp->Draw6(c1,1,im+1);

	  cp->TransLegend(-0.4,-0.55);
	  cp->WidenLegend(0.17,0.);
	  c1->Update();
	}
      }
    }

    // -----------------
    // Total systematics
    // ------------------
    if (totSyst) {

      std::vector<ComparisonPlot_t*> cpSystV;
    
      int noComparison=1;
      
      c1 = new TCanvas("c1","c1", 1200,900-noComparison*150);
      if (noComparison) {
	c1->Divide(3,2);
	SetSideSpaces(c1,0.05,0.,0.,0.02);
      }

      for (int im=0; im<6; ++im) {
	TString title="Total systematics ";
	TString mRange=Form("M_%2.0lf_%2.0lf",DYTools::massBinLimits[im+1],DYTools::massBinLimits[im+2]);
	std::vector<TH1D*> *hProfV= hProfW[im];
	if (!hProfV) {
	  std::cout << "hProfV for " << title << mRange << " is null\n";
	  continue;
	}
	
	ComparisonPlot_t *cpSyst=NULL;
	cpSyst = new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,title+mRange,"","|y|",correctionKind,"ratio");
	cpSystV.push_back(cpSyst);
	
	//cpSyst->SetLogx();
	cpSyst->SetPrintValues(0);
	if (noComparison) cpSyst->SetRefIdx(-111);
	cpSyst->SetYTitleSize(0.07,1.4);
	cpSyst->SetRatioYTitleSize(0.17,0.47);
      
	if ((*hProfV)[3]) cpSyst->AddHist1D((*hProfV)[3],combiLabels[3],"LPE1",kRed+1,1,0);
	if (0 && (*hProfV)[2]) cpSyst->AddHist1D((*hProfV)[2],combiLabels[2],"LPE3",kGreen+1,1,0);
	if (0 && (*hProfV)[1]) cpSyst->AddHist1D((*hProfV)[1],combiLabels[1],"LPE3",kBlue,1,0);
	if ((*hProfV)[0]) cpSyst->AddHist1D((*hProfV)[0],combiLabels[0],"LPE1",kBlack,1,0);
	//if ((*hProfV)[0]) cpSyst->AddHist1D((*hProfV)[0],combiLabels[0],"LPE3",kBlack,2,0,-1);
      
	if (noComparison) {
	  cpSyst->Draw(c1,false,"png",im+1);
	}
	else {
	  if (im==0) cpSyst->Prepare6Pads(c1,1);
	  cpSyst->Draw6(c1,1,im+1);
	}
	cpSyst->TransLegend(-0.35,-0.57);
	cpSyst->WidenLegend(0.17,0.);
	c1->Update();
      }
    }
 


  }
  else {
    // -----------  1D case ----------------------

    std::vector<TH1D*> hProfRawV;
    for (unsigned int i=0; i<histos.size(); ++i) {
      TH1D* hProf=NULL;
      if (histos[i]) {
	TString name=Form("hProfRaw_%s",labelsV[i].c_str());
	hProf= createProfileX(histos[i],1,name,1,name);
	hProf->GetXaxis()->SetMoreLogLabels();
	hProf->GetXaxis()->SetNoExponent();
	hProf->SetMarkerStyle(20);
	hProf->SetMarkerSize(0.8);
      }
      hProfRawV.push_back(hProf);
    }

    std::vector<TH1D*> hProfV;
    for (unsigned int i=0; i<hCombiV.size(); ++i) {
      TString name=Form("hProf_%s",combiLabels[i].Data());
      if (!hCombiV[i]) { HERE("combiV[%u] is null",i); return ; }
      TH1D *h=createProfileX(hCombiV[i],1,name,1,name);
      h->GetXaxis()->SetMoreLogLabels();
      h->GetXaxis()->SetNoExponent();
      hProfV.push_back( h );
    }

    // -----------------
    // Partial systematics
    // ------------------
    if (PUsyst+FSRsyst) {
      ComparisonPlot_t *cpFSRsyst=NULL, *cpPUsyst=NULL;
    
      for (int isPU=0; isPU<2; ++isPU) {
	if ( isPU && !PUsyst) continue;
	if (!isPU && !FSRsyst) continue;
	TString title=(isPU) ? "Pile-up systematics" : "FSR systematics";
	ComparisonPlot_t *cp = new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,"cpSyst",title,"#it{M}_{ee} [GeV]",correctionKind,"ratio");
	if (isPU) cpPUsyst=cp; else cpFSRsyst=cp;

	cp->SetLogx();
	cp->SetPrintValues(0);
	cp->SetYTitleSize(0.07,1.);
	cp->SetRatioYTitleSize(0.17,0.47);
      
	int imin=(isPU) ? 5 : 3;
	if (hProfRawV[0]) {
	  removeError(hProfRawV[0]);
	  cp->AddHist1D(hProfRawV[0],labelsV[0],"LP",kBlack,1,0);
	}
	if (hProfRawV[imin]) {
	  removeError(hProfRawV[imin]);
	  cp->AddHist1D(hProfRawV[imin],labelsV[imin],"LP",kBlue,1,0);
	}
	if (hProfRawV[imin+1]) {
	  removeError(hProfRawV[imin+1]);
	  cp->AddHist1D(hProfRawV[imin+1],labelsV[imin+1],"LP",kGreen+1,1,0);
	}
      
	TString canvName=(isPU) ? "cPU" : "cFSR";
	c1 = new TCanvas(canvName,canvName, 600,700);
	if (isPU) cPU=c1; else cFSR=c1;
	cp->Prepare2Pads(c1);
	cp->Draw(c1);

	cp->TransLegend(-0.1,-0.5);
	cp->WidenLegend(0.15,0.);
	c1->Update();
      }
    }

    // -----------------
    // Total systematics
    // ------------------
    if (totSyst) {
      ComparisonPlot_t *cpSyst=NULL;
    
      int noComparison=1;

      cpSyst = new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,"cpSyst","","#it{M}_{ee} [GeV]",correctionKind,"ratio");
      cpSyst->SetLogx();
      cpSyst->SetPrintValues(1);
      if (noComparison) cpSyst->SetRefIdx(-111);
      cpSyst->SetYTitleSize(0.07,1.);
      cpSyst->SetRatioYTitleSize(0.17,0.47);
      
      if (hProfV[3]) cpSyst->AddHist1D(hProfV[3],combiLabels[3],"LPE1",kRed+1,1,0);
      if (0 && hProfV[2]) cpSyst->AddHist1D(hProfV[2],combiLabels[2],"LPE3",kGreen+1,1,0);
      if (0 && hProfV[1]) cpSyst->AddHist1D(hProfV[1],combiLabels[1],"LPE3",kBlue,1,0);
      if (hProfV[0]) cpSyst->AddHist1D(hProfV[0],combiLabels[0],"LPE1",kBlack,1,0);
      //if (hProfV[0]) cpSyst->AddHist1D(hProfV[0],combiLabels[0],"LPE3",kBlack,2,0,-1);
      
      c1 = new TCanvas("c1","c1", 600,700-noComparison*100);
      if (noComparison) {
	cpSyst->Draw(c1,false,"png",0);
      }
      else {
	cpSyst->Prepare2Pads(c1);
	cpSyst->Draw(c1);
      }
      cpSyst->TransLegend(-0.1,-0.55);
      cpSyst->WidenLegend(0.15,0.05);
      c1->Update();
    }
  }

  TString dirName="plots-syst";
  TString fnameBase=Form("fig-%s-%s-",correctionKind.Data(),DYTools::analysisTag.Data());
  TString fnamePU= fnameBase + TString("PileUpSyst");
  TString fnameFSR=fnameBase + TString("FSRsyst");
  TString fnameTot=fnameBase + TString("total");
  std::cout << "plot file names:\n";
  if (totSyst) std::cout << " - " << fnameTot << "\n";
  if (FSRsyst) std::cout << " - " << fnameFSR << "\n";
  if (PUsyst) std::cout << " - " << fnamePU << "\n";
  if (!saveCanvas) { std::cout << ".. not saved (as requested)\n"; }
  else {
    if (c1) SaveCanvas(c1,fnameTot,dirName);
    if (cFSR) SaveCanvas(cFSR,fnameFSR,dirName);
    if (cPU) SaveCanvas(cPU,fnamePU,dirName);
  }
  
  return;
}
