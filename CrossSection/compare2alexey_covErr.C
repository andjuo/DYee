#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/ComparisonPlot.hh"
#include "../CrossSection/crossSectionFnc.hh"
#include "CSCovWorkFlags.hh"
#include <fstream>
#include <sstream>

// --------------------------------------------

TH1D* loadAlexey1D(int rshape=0);


//=== MAIN MACRO =================================================================================================

int compare2alexey_covErr(int analysisIs2D=0,
			  TString csFNameExtraTag="",
			  int savePlots=0,
		 DYTools::TCrossSectionKind_t csKind=DYTools::_cs_None
		   )
{
  if (analysisIs2D) {
    std::cout << "not ready for 2D\n";
    return retCodeError;
  }

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

  TString inpFName=Form("table_%dD_frac.root",DYTools::study2D+1);
  TString figNameTag=csFNameExtraTag;

  if (csFNameExtraTag.Length()) {
    inpFName.ReplaceAll(".root",csFNameExtraTag + TString(".root"));
    
  }
  std::cout << "input fname is <" << inpFName << ">\n";

  // ------------------------------------------------

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

  std::vector<TH1D*> theoryV;

  //std::vector<HistoPair2D_t*> xsecHPV;
  std::vector<TH1D*> xsecStatErrV;
  std::vector<TH1D*> xsecTotErrV;


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

  TString fieldXsec="xsec";
  TString fieldStatErr="signal_stat_err";
  TString fieldTotalErr="total_err";

  if (res) {
    TFile fin(inpFName,"read");
    if (!fin.IsOpen()) {
      std::cout << "failed to open a file <" << inpFName << ">\n";
      return retCodeError;
    }

    TH1D* hXS=createBaseH1(fieldXsec,"Xsec","#sigma");
    TH1D* hStatErr=createBaseH1(fieldStatErr,"Stat.err","stat err");
    TH1D *hTotalErr=createBaseH1(fieldTotalErr,"Tot.err","total err");

    res= (loadHisto(fin,&hXS,"") &&
	  loadHisto(fin,&hStatErr,"") &&
	  loadHisto(fin,&hTotalErr,"")) ? 1:0;
    fin.Close();
    if (!res) {
      std::cout << "failed to load histos\n";
      return retCodeError;
    }

    TH1D *hXStotErr=Clone(hXS,"hXStotErr","hXStotErr");
    for (int ibin=1; ibin<=hXS->GetNbinsX(); ++ibin) {
      double xs=hXS->GetBinContent(ibin);
      double extraErr=(ibin==15) ? 0.25 : 0.;
      extraErr=0;
      hXS->SetBinError(ibin, xs*hStatErr->GetBinContent(ibin));
      hXStotErr->SetBinError(ibin, xs*(hTotalErr->GetBinContent(ibin)+extraErr));
    }

    

    delete hStatErr;
    delete hTotalErr;

    xsecStatErrV.push_back(hXS);
    xsecTotErrV.push_back(hXStotErr);
  }

  TH1D* hAlexey=loadAlexey1D(0);
  if (!hAlexey) return retCodeError;

  //printHisto(hAlexey); return retCodeStop;

  // ------------------------------------------------

  /*
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
  */

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
    //removeError(theoryV[0]);
    cp->AddHist1D(theoryV[0], "CTEQ12NNLO", "L", kBlue, 1,0);

    // calculation
    int color1=kBlack;
    int color2=kGreen+1;//46; // red-brown
    //cp->AddHist1D(xsecTotErrV[0], "Andrius tot err", "LPE1", color2, 1,0);
    //cp->AddHist1D(xsecStatErrV[0],"Andrius stat err", "LPE1", TAttMarker(color1,24,1), 1,0);
    cp->AddHist1D(hAlexey,"Alexey","LPE1",TAttMarker(9,kMultiply,1.2),1,0,1);

    TCanvas *cx=new TCanvas("cx","cx",700,800);
    cp->Prepare2Pads(cx);
    cp->Draw(cx);
    cp->TransLegend(-0.05,0);
    //cp->WidenLegend(0,0);
    cx->Update();

    if (savePlots) {
      TString figName=TString("figErrFromCov-1D-") + figNameTag;
      eliminateSeparationSigns(figName);
      SaveCanvas(cx,figName);
    }
    //}

  return retCodeOk;
}

// ----------------------------------------------------------
// ----------------------------------------------------------

TH1D* loadAlexey1D(int rshape) {

  TString fname="preFSRfullAcc1D_EE.txt";
  std::ifstream fin(fname);
  if (!fin.is_open()) {
    std::cout << "failed to open the file <" << fname << ">\n";
    return NULL;
  }

  TH1D *h1= createBaseH1("hAlexey","hAlexey","#sigma");

  std::string s, tmp;
  int ibin=1;
  double r,rerr;
  while (!fin.eof() && getline(fin,s)) {
    std::stringstream ss(s);
    ss >> tmp >> ibin >> r >> rerr;
    std::cout << "loaded <s>\n"
	      << "  decoded tmp=<" << tmp << ">, "
	      << r << " +- " << rerr << "\n";
    h1->SetBinContent(ibin, r   );
    h1->SetBinError  (ibin, rerr);
  }

  fin.close();

  if (!rshape) {
    NormCS_t sigmaZ(DYTools::_cs_None,"");
    if (!sigmaZ.loadValues("dir-CSInfo/inpCS_sigmaZ.dat") ||
	(sigmaZ.cs()==double(0.))) {
      std::cout << "failed to load the cross section\n";
      return NULL;
    }

    double norm=(0) ?  sigmaZ.cs() : 1.;
    TH1D *hAbs=convert_rshape2absolute(h1, norm);

    if (0) {
      TH1D *hRchk=convert_absolute2rshape(hAbs, norm);
      TH1D *hDiff=Clone(h1,"hdiff","hdiff");
      if (0) {
	hDiff->Add(hRchk,-1);
      }
      else {
	for (ibin=1; ibin<=hDiff->GetNbinsX(); ++ibin) {
	  hDiff->SetBinContent(ibin,
		       h1->GetBinContent(ibin) - hRchk->GetBinContent(ibin));
	  hDiff->SetBinError (ibin,
			    h1->GetBinError(ibin) - hRchk->GetBinError(ibin));
	}
      }
      printHisto(hDiff);
      delete hRchk;
      delete hDiff;
      return NULL;
    }

    swapPtrs(&hAbs,&h1);
    delete hAbs;
  }

  return h1;
}

// ----------------------------------------------------------
// ----------------------------------------------------------
