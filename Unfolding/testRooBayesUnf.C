#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/UnfoldingMatrix.h"

int testRooBayesUnf(int analysisIs2D) {
  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return retCodeError;
  }

  const int markerCount=5;
  const TAttMarker markerV[markerCount] = {
    TAttMarker(kBlack,24,1.),
    TAttMarker(kBlue,20,0.9),
    TAttMarker(kRed+1,22,0.9),
    TAttMarker(kGreen+1,4,0.9),
    TAttMarker(kViolet,3,0.9)
  };


  TString path="../../Results-DYee/root_files_reg/";
  TString yieldFName=Form("yield/DY_j22_19712pb_ApplyEscale/bg-subtracted_yield_%dD__peak20140220.root",analysisIs2D+1);
  TString yieldField="signalYieldDDbkg";

  TH2D* h2SigYield=LoadHisto2D(yieldField,path+yieldFName,"",1);
  if (!h2SigYield) return retCodeError;

  //TString detRespFName=Form("constants/DY_j22_19712pb/detResponse_unfolding_constants%dD.root",analysisIs2D+1);
  TString detRespName="detResponse";
  //detRespName="detResponseExact";
  UnfoldingMatrix_t *detResponse= new UnfoldingMatrix_t(UnfoldingMatrix::_cDET_Response,detRespName);
  if (!detResponse) return retCodeError;
  TString constDirDef=path + TString("constants/DY_j22_19712pb/");
  TString fnameTagDef=UnfoldingMatrix_t::generateFNameTag(DYTools::NO_SYST,0);
  int res=detResponse->autoLoadFromFile(constDirDef,fnameTagDef);
  //res=detResponse->autoLoadFromFile_forRooUnfold(constDirDef,fnameTagDef);
  if (!res) return retCodeError;

  TH2D* h2Unf= Clone(h2SigYield,"h2Unf");
  if (!unfold_reco2true(h2Unf,*detResponse,h2SigYield)) return retCodeError;
  TH1D* hSigYield_flat= flattenHisto(h2SigYield,"hSigYield_flat");
  TH1D* hUnfYield_flat= flattenHisto(h2Unf, "hUnfYield_flat");

  // ---------------------------------------------------------------------
  //  Study part
  // ---------------------------------------------------------------------

  RooUnfBayes_t rooUnf(*detResponse);

  // Closure test
  if (0) {
    UnfoldingMatrix_t *U= detResponse;
    TH2D* h2MCReco= createHisto2D(*U->getFinM(),U->getFinMerr(),
				  "h2MCReco","h2MCReco",0,0,0.);
    TH1D* hMCReco_flat= flattenHisto(h2MCReco,"hMCReco_flat");
    if (!hMCReco_flat) return retCodeError;
    TH2D* h2MCRecoUnf= Clone(h2MCReco,"h2MCRecoUnf");
    if (!unfold_reco2true(h2MCRecoUnf,*detResponse,h2MCReco))
      return retCodeError;

    TH2D* h2MCTrue= createHisto2D(*U->getIniM(),U->getIniMerr(),
				  "h2MCTrue","h2MCTrue",0,0,0.);
    TH1D* hMCTrue_flat= flattenHisto(h2MCTrue,"hMCTrue_flat");
    TH1D* hMCRecoUnf_flat= flattenHisto(h2MCRecoUnf,"hMCRecoUnf_flat");
    ComparisonPlot_t *cp= new ComparisonPlot_t("cp","closure",
					       "bin","count","ratio");

    printZpeakCount(h2MCReco);
    printZpeakCount(h2MCTrue);
    printZpeakCount(h2MCRecoUnf);

    cp->AddHist1D(hMCReco_flat,"MC reco","LP",TAttMarker(kRed,5,1),1,0,1);
    cp->AddHist1D(hMCTrue_flat,"MC true (unf) distr",
		  "LP",markerV[0],1,0,1);
    cp->AddHist1D(hMCRecoUnf_flat,"unf.yield (matrix inv.)",
		  "LP",TAttMarker(kOrange,20,1.),1,0,1);
    if (1)
    for (int iStep=1; iStep<=4; ++iStep) {
      //if (iStep>1) break;
      if (!rooUnf.doUnfoldBayes(hMCReco_flat,iStep)) return retCodeError;
      TH1D* h1= rooUnf.unfYield();
      printZpeakCount(h1);
      cp->AddHist1D(h1, Form("rooUnf. (%d steps)",iStep),
		    "LP",markerV[iStep%markerCount],iStep/markerCount+1,0,1);
    }

    HERE("draw");
    TCanvas *cx=new TCanvas("cx","cx",800,800);
    cp->Prepare2Pads(cx);
    cp->Draw(cx);
    cx->Update();
    return retCodeStop;
  }

  // Unfolding test
  if (1) {
    UnfoldingMatrix_t *U= detResponse;
    TH2D* h2MCTrue= createHisto2D(*U->getIniM(),U->getIniMerr(),
				  "h2MCTrue","h2MCTrue",0,0,0.);
    TH1D* hMCTrue_flat= flattenHisto(h2MCTrue,"hMCTrue_flat");

    ComparisonPlot_t *cp= new ComparisonPlot_t("cp","unfolding",
					       "bin","count","ratio");
    cp->SetLogy(0);

    ComparisonPlot_t *cpErr= NULL;
    if (1) {
      cpErr = new ComparisonPlot_t("cpErr","unfolding error",
				   "bin","count error","ratio");
    }

    //cp->AddHist1D(hSigYield_flat,"sig.yield","LP",TAttMarker(kRed,5,1),1,0,1);

    hMCTrue_flat->Scale(6.8);
    if (0) {
      cp->AddHist1D(hMCTrue_flat,"MC gen distr (scaled)",
		    "LP",markerV[0],1,0,1);
      if (cpErr) {
	TH1D *hMCTrue_flatErr= Clone(hMCTrue_flat,"hMCTrue_flatErr");
	setErrorAsContent(hMCTrue_flatErr);
	cpErr->AddHist1D(hMCTrue_flatErr,"MC gen distr (scaled)",
			 "LP",markerV[0],1,0,1);
      }
    }

    if (1) {
      cp->AddHist1D(hUnfYield_flat,"unf.sig.yield (matrix inv.)",
		    "LP",TAttMarker(kOrange,20,1.),1,0,1);
      if (cpErr) {
	TH1D* hUnfYield_flatErr= Clone(hUnfYield_flat,"hUnfYield_flatErr");
	setErrorAsContent(hUnfYield_flatErr);
	cpErr->AddHist1D(hUnfYield_flatErr,"unf.sig.yield (matrix inv.)",
			 "LP",TAttMarker(kOrange,20,1.),1,0,1);
      }
    }

    printZpeakCount(h2SigYield);
    printZpeakCount(hSigYield_flat);
    printZpeakCount(hUnfYield_flat);

    if (1)
    for (int iStep=1; iStep<=4; ++iStep) {
      //if (iStep!=4) continue;
      //if (iStep!=10) continue;
      //if (iStep%2==0) continue;
      //if (iStep>1) break;
      if (!rooUnf.doUnfoldBayes(hSigYield_flat,iStep)) return retCodeError;
      TH1D* h1= rooUnf.unfYield();
      printZpeakCount(h1);
      //std::cout << dashline << "iStep=" << iStep << "\n\n"; printHisto(h1);
      TString loc_label=Form("rooUnf. (%d steps)",iStep);
      cp->AddHist1D(h1, loc_label,
		    "LP",markerV[iStep%markerCount],iStep/markerCount+1,0,1);
      if (cpErr) {
	TH1D *h1err= Clone(h1,Form("hErr_iStep%d",iStep));
	setErrorAsContent(h1err);
	cpErr->AddHist1D(h1err, loc_label,
		 "LP",markerV[iStep%markerCount],iStep/markerCount+1,0,1);
      }
    }

    HERE("draw");
    TCanvas *cx=new TCanvas("cx","cx",800,800);
    cp->Prepare2Pads(cx);
    cp->Draw(cx);
    cx->Update();

    if (cpErr) {
      TCanvas *cxErr=new TCanvas("cxErr","cxErr",800,800);
      cpErr->Prepare2Pads(cxErr);
      cpErr->Draw(cxErr);
      cxErr->Update();
    }

    return retCodeStop;
  }

  return retCodeOk;
}
