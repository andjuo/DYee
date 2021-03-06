#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/ComparisonPlot.hh"
#include "../CrossSection/crossSectionFnc.hh"

//=== MAIN MACRO =================================================================================================

int compare2theory(int analysisIs2D,
		   TString confFile,
		   TString csFNameExtraTag="",
		   TString loadDataField="auto",
		   int savePlots=0,
		   DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST,
		   DYTools::TCrossSectionKind_t csKind=DYTools::_cs_None,
		   int special_case=-1
		   )
{

  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return retCodeError;
  }

  if (csKind==DYTools::_cs_None) {
    csKind=(analysisIs2D) ? DYTools::_cs_preFsrDet : DYTools::_cs_preFsr;
    std::cout << "default csKind " << CrossSectionKindName(csKind) << "\n";
  }

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================

  DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN;

  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(confFile)) return retCodeError;

  // Construct eventSelector, update mgr and plot directory
  TString extraTag=""; // empty
  TString plotExtraTag;

  EventSelector_t evtSelector(inpMgr,runMode,systMode,
			      extraTag,plotExtraTag,
			      EventSelector::_selectDefault);


  TString figNameTag= CrossSectionKindName(csKind);
  if (csFNameExtraTag.Length()) {
    figNameTag.Append(TString("-") + csFNameExtraTag);
  }
  if (loadDataField!="auto") {
    figNameTag.Append(TString("-") + loadDataField);
  }

  // ------------------------------------------------

  std::vector<TString> theoryFilesV;
  TString theory_path="../root_files_reg/theory/FEWZ3_prediction/";
  TString theory_field;
  if (DYTools::study2D) {
    theoryFilesV.reserve(6);
    TString path=theory_path;
    theory_field="xsec";
    theoryFilesV.push_back(path + TString("2Dabsxsec20to30_NNLO_cuts_CTEQ12NNLO-mdf.root"));
    theoryFilesV.push_back(path + TString("2Dabsxsec30to45_NNLO_cuts_CTEQ12NNLO-mdf.root"));
    theoryFilesV.push_back(path + TString("2Dabsxsec45to60_NNLO_cuts_CTEQ12NNLO-mdf.root"));
    theoryFilesV.push_back(path + TString("2Dabsxsec60to120_NNLO_cuts_CTEQ12NNLO-mdf.root"));
    theoryFilesV.push_back(path + TString("2Dabsxsec120to200_NNLO_cuts_CTEQ12NNLO-mdf.root"));
    theoryFilesV.push_back(path + TString("2Dabsxsec200to1500_NNLO_cuts_CTEQ12NNLO-mdf.root"));
  }
  else {
    theoryFilesV.push_back(theory_path +
			   TString("1Dabsxsec_NNLO_CTEQ12NNLO_41.root"));
    theory_field="invm_FEWZ_41";
  }

  // ------------------------------------------------

  TString inpFName=inpMgr.crossSectionFullFileName(systMode,csKind,0,0);
  if (csFNameExtraTag.Length()) {
    inpFName.ReplaceAll(".root",csFNameExtraTag + TString(".root"));
  }
  int nBayesIters=-1;
  if ((special_case>10) && (special_case<20)) {
    nBayesIters= special_case - 10;
    inpFName.ReplaceAll("DY_j22_19712pb",Form("DY_j22_19712pb-BayesUnf%d",nBayesIters));
  }
  std::cout << "input fname is <" << inpFName << ">\n";

  // ------------------------------------------------

  std::vector<TH1D*> theoryV;

  std::vector<HistoPair2D_t*> xsecHPV;
  std::vector<std::vector<TH1D*>*> xsecStatErrVV;
  std::vector<std::vector<TH1D*>*> xsecTotErrVV;


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================

  // ------------------------------------------------

  for (unsigned int i=0; i<theoryFilesV.size(); ++i) {
    TFile finTh(theoryFilesV[i],"read");
    TH1D* h=(TH1D*)finTh.Get(theory_field);
    if (!h) {
      finTh.Close();
      std::cout << "failed to get " << theory_field << " from file "
		<< finTh.GetName() << "\n";
      return retCodeError;
    }
    h->SetDirectory(0);
    h->SetMarkerStyle(24);
    h->SetMarkerSize(0.7);
    theoryV.push_back(h);
    finTh.Close();
      
    TString hTitle=finTh.GetName();
    Ssiz_t pos=hTitle.Last('/');
    hTitle.Remove(0,pos+1);
    pos=hTitle.Index(".root");
    hTitle.Remove(pos,hTitle.Length());
    std::cout << "hTitle=" << hTitle << "\n";
    h->SetTitle(hTitle);

    if (DYTools::study2D==0) h->GetYaxis()->SetRangeUser(1e-8,1e4);
  }

  // ------------------------------------------------

  int res=1;
  int dataMarkerStyle=20;
  double dataMarkerSize=0.8;

  if (res) {
    TFile fin(inpFName,"read");
    if (!fin.IsOpen()) {
      std::cout << "failed to open a file <" << inpFName << ">\n";
      return retCodeError;
    }
    
    // load xsec
    HistoPair2D_t *hp= new HistoPair2D_t("xsec");
    if (!hp) res=0; else xsecHPV.push_back(hp);
    TString fieldName=loadDataField;
    if (loadDataField=="auto") {
      fieldName=(csFNameExtraTag.Length()) ? "xsec" : "count";
      std::cout << "loading field <" << fieldName << ">\n";
      if (res) res=hp->Read(fin,"auto",fieldName,1);
    }
    else {
      std::cout << "loading field <" << fieldName << ">\n";
      if (res) res=hp->Read(fin,"",fieldName,1);
    }
    fin.Close();
    if (!res) {
      std::cout << "failed to get field <" << fieldName << ">\n";
      return retCodeError;
    }

    std::vector<TH2D*> xsStatErrH2V, xsTotErrH2V;
    for (unsigned int i=0; i<xsecHPV.size(); ++i) {
      TH2D* h2=xsecHPV[i]->createHistoWithFullError(Form("h2_fullErr_%d",i));
      xsTotErrH2V.push_back(h2);
      xsStatErrH2V.push_back(xsecHPV[i]->histo());
    }
    std::vector<TString> labelsTotErrV, labelsStatErrV;
    if (DYTools::study2D) {
      for (int im=0; im<DYTools::nMassBins; ++im) {
	//TString mStr=Form("M_%1.0lf_%1.0lf",
	//		  DYTools::massBinLimits[im],
	//		  DYTools::massBinLimits[im]);
	labelsTotErrV.push_back(TString("wTotErr_"));
	labelsStatErrV.push_back(TString("wStatErr_"));
      }
      if (!createRapidityProfileVec(xsTotErrH2V,xsecTotErrVV,labelsTotErrV,
				    dataMarkerStyle,dataMarkerSize) ||
	  !createRapidityProfileVec(xsStatErrH2V,xsecStatErrVV,labelsStatErrV,
				    dataMarkerStyle,dataMarkerSize)){
	std::cout << "failed to create rapidity profiles\n";
	return retCodeError;
      }

      std::cout << "xsecTotErrVV.size=" << xsecTotErrVV.size() << "\n";
      std::cout << "xsecTotErrVV[0]->size()=" << xsecTotErrVV[0]->size() << "\n";
      int firstBinRemoved=0;
      if (DYTools::massBinLimits[0]==double(0)) {
	firstBinRemoved=1;
	std::cout << "\ndeleting 1st mass bin ("
		  << DYTools::massBinLimits[0] << " .. "
		  << DYTools::massBinLimits[1] << ")\n";
	std::vector<TH1D*> *vec1= xsecTotErrVV[0];
	std::vector<TH1D*> *vec2= xsecStatErrVV[0];
	xsecTotErrVV.erase( xsecTotErrVV.begin() );
	xsecStatErrVV.erase( xsecStatErrVV.begin() );
	if (1) {
	  for (unsigned int i=0; i<vec1->size(); ++i) {
	    std::cout << " - removing " << (*vec1)[i]->GetName() << "\n";
	  }
	  for (unsigned int i=0; i<vec2->size(); ++i) {
	    std::cout << " - removing " << (*vec2)[i]->GetName() << "\n";
	  }
	}
	std::cout << "\n";
	ClearVec(*vec1);
	ClearVec(*vec2);
	delete vec1;
	delete vec2;
      }

      for (int ivv=0; ivv<2; ++ivv) {
	std::vector<std::vector<TH1D*>*>* vv=NULL;
	vv= (ivv==0) ? &xsecTotErrVV : &xsecStatErrVV;
	for (unsigned int ivec=0; ivec < vv->size(); ++ivec) {
	  double factor=
	    (DYTools::nYBins[ivec+firstBinRemoved]==24) ? 10. : 5.;
	  std::vector<TH1D*> *vec= (*vv)[ivec];
	  for (unsigned int ii=0; ii < vec->size(); ++ii) {
	    std::cout << " " << (*vec)[ii]->GetName()
		      << " scaled by factor 1/" << (1./factor) << "\n";
	    (*vec)[ii]->Scale(factor);
	  }
	}
      }
    }
    else {
      labelsTotErrV.push_back("totErr_");
      labelsStatErrV.push_back("statErr_");
      if (!createMassProfileVec(xsTotErrH2V,xsecTotErrVV,labelsTotErrV,
				dataMarkerStyle,dataMarkerSize) ||
	  !createMassProfileVec(xsStatErrH2V,xsecStatErrVV,labelsStatErrV,
				dataMarkerStyle,dataMarkerSize)) {
	std::cout << "failed to create mass profiles\n";
	return retCodeError;
      }
      std::cout << "xsecTotErrVV.size=" << xsecTotErrVV.size() << "\n";
      std::cout << "xsecTotErrVV[0]->size()=" << xsecTotErrVV[0]->size() << "\n";
    }
    // clear total error vector, but not stat error vector
    ClearVec(xsTotErrH2V);
  }

  // ------------------------------------------------

  if (DYTools::study2D) {
    for (unsigned int im=0; im<theoryV.size(); ++im) {
      TString mStr=Form("M_%1.0lf_%1.0lf", DYTools::massBinLimits[im+1],
			DYTools::massBinLimits[im+2]);
      std::cout << "plotting " << mStr << "\n";
      ComparisonPlot_t *cp=NULL;
      cp= new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,
			       TString("cp") + mStr,
			       csFNameExtraTag + TString(" ") + mStr,
			       "|y|",TString("cross section"),
			       "ratio");
      cp->SetPrintValues(0);
      cp->SetRatioLabelSize(0.11);

      // theory
      cp->AddHist1D(theoryV[im], "CTEQ12NNLO", "LP", kBlue, 1,0);

      std::cout << "theoryV[im]->GetNbinsX()=" << theoryV[im]->GetNbinsX() << "\n";
      std::cout << "(*xsecTotErrVV[im])[0]->GetNbinsX()=" << (*xsecTotErrVV[im])[0]->GetNbinsX() << "\n";

      // calculation
      int color1=kBlack;
      int color2=46; // red-brown
      cp->AddHist1D((*xsecTotErrVV[im])[0], "tot err", "LPE1", color2, 1,0);
      cp->AddHist1D((*xsecStatErrVV[im])[0],"stat err", "LPE1", color1, 1,0);

      TString canvName=Form("cx%d",im);
      TCanvas *cx=new TCanvas(canvName,canvName,700,800);
      cp->Prepare2Pads(cx);
      cp->Draw(cx);
      cp->TransLegend(-0.05,0);
      //cp->WidenLegend(0,0);
      cx->Update();

      if (savePlots) {
	TString figName=TString("fig-2D-") + figNameTag +
	  TString("-") + mStr;
	eliminateSeparationSigns(figName);
	SaveCanvas(cx,figName);
      }
    }
  }
  else {
    ComparisonPlot_t *cp=NULL;
    cp= new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,"cpXS",
			     csFNameExtraTag,
			     "M_{ee} [GeV]",TString("cross section"),
			     "ratio");
    cp->SetPrintValues(1);
    cp->SetRatioLabelSize(0.11);
    cp->SetLogx(1);
    cp->SetLogy(1);

    // theory
    cp->AddHist1D(theoryV[0], "CTEQ12NNLO", "L", kBlue, 1,0);

    // calculation
    int color1=kBlack;
    int color2=46; // red-brown
    cp->AddHist1D((*xsecTotErrVV[0])[0], "tot err", "LPE1", color2, 1,0);
    cp->AddHist1D((*xsecStatErrVV[0])[0],"stat err", "LPE1", color1, 1,0);

    TCanvas *cx=new TCanvas("cx","cx",700,800);
    cp->Prepare2Pads(cx);
    cp->Draw(cx);
    cp->TransLegend(-0.05,0);
    //cp->WidenLegend(0,0);
    cx->Update();

    if (savePlots) {
      TString figName=TString("fig-1D-") + figNameTag;
      eliminateSeparationSigns(figName);
      SaveCanvas(cx,figName);
    }
  }

  return retCodeOk;
}
