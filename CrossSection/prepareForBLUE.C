#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "CSCovWorkFlags.hh"
#include <fstream>

// ------------------------------------------------------------

// if h2CSsystErr=NULL, the statistical error is plotted

int compareError(const TH2D* h2CS, const TMatrixD &cov,
		 int relativeError=1,
		 TH2D* h2CSsystErr=NULL, int totalErr=1);


// ------------------------------------------------------------
//    Main macro
// ------------------------------------------------------------

int prepareForBLUE(int analysisIs2D) {
  if (analysisIs2D) {
    std::cout << "code is not ready of 2D case\n";
  }
  if (!DYTools::setup(analysisIs2D)) return retCodeError;

  TH2D *h2CS=NULL, *h2CSsystErr=NULL;
  h2CS=loadMainCSResult(1,&h2CSsystErr);
  if (!h2CS) return retCodeError;

  if (1) {
    printHisto(h2CS);
    printHisto(h2CSsystErr);
  }

  int test=1;
  TString covFName=(!analysisIs2D) ? "finalCov-1D.root" : "finalCov-2D.root";
  if (test) covFName.ReplaceAll(".root","-yieldStatOnly.root");
  //TFile fCov(covFName,"read");
  //if (!fCov.
  int nBins=DYTools::nUnfoldingBins;
  //if (analysisIs2D) nBins-=24;
  TMatrixD* covPtr=loadMatrix(covFName,"totalCov",nBins,nBins,1);
  if (!covPtr) return retCodeError;
  TMatrixD cov(*covPtr);
  delete covPtr;

  if (0) {
    TCanvas *ct=new TCanvas("ct","ct",700,700);
    AdjustFor2DplotWithHeight(ct);
    gStyle->SetPalette(1);
    cov.Draw("COLZ");
    ct->Update();
    return retCodeStop;
  }

  // Compare the statistical error
  if (1) {
    if (1 && test) {
      int relativeErr=1;
      if (!compareError(h2CS,cov,0,NULL,1)) return retCodeError;
      if (!compareError(h2CS,cov,relativeErr,NULL,1)) return retCodeError;
      return retCodeStop;
    }
    // Compare the systematic error
    if (1 && !test) {
      int relativeErr=1;
      if (!compareError(h2CS,cov,0,h2CSsystErr,1)) return retCodeError;
      if (!compareError(h2CS,cov,relativeErr,h2CSsystErr,1)) return retCodeError;
      return retCodeStop;
    }
  }

  // error from cross section calculation
  TH2D *h2errFromCS=Clone(h2CS,"h2errFromCS");
  swapContentAndError(h2errFromCS);
  removeError(h2errFromCS);
  // error from covariance
  TH2D *h2errFromCov= errorFromCov(cov,"h2errFromCov");

  // save xsec error
  if (1) {
    TString outFNameBase="xsecForBlue-" + DYTools::analysisTag;
    if (test) outFNameBase.Append("-yieldStatOnly");

    for (int iSrc=0; iSrc<2; ++iSrc) {
      TH2D *h2err=(iSrc==0) ? h2errFromCS : h2errFromCov;
      TString outFName=outFNameBase;
      outFName.Append((iSrc==0) ? "-fromCS" : "-fromCov");
      outFName.Append(".dat");

      std::ofstream fout(outFName);
      TMatrixD *yV= DYTools::getYBinLimits();
      //yV->Print();
      for (int ibin=1; ibin<=h2err->GetNbinsX(); ++ibin) {
	for (int jbin=1; jbin<=h2err->GetNbinsY(); ++jbin) {
	  if (analysisIs2D) {
	    if (ibin==1) continue;
	    if ((ibin==7) && (jbin>12)) break;
	    fout << Form("%3d  %3d  %6.1lf  %6.1lf  %5.2lf %5.2lf  %.16e   %.16e\n",
			 ibin,jbin,
			 h2CS->GetBinLowEdge(ibin),
			 h2CS->GetBinLowEdge(ibin) + h2CS->GetBinWidth(ibin),
			 (*yV)(ibin-1,jbin-1),(*yV)(ibin-1,jbin),
			 h2CS->GetBinContent(ibin,jbin),
			 h2err->GetBinContent(ibin,jbin));
	  }
	  else {
	    fout << Form("%2d   %6.1lf   %6.1lf   %.16e   %.16e\n",
			 ibin,
			 h2CS->GetBinLowEdge(ibin),
			 h2CS->GetBinLowEdge(ibin) + h2CS->GetBinWidth(ibin),
			 h2CS->GetBinContent(ibin,jbin),
			 h2err->GetBinContent(ibin,jbin));
	  }
	}
      }
      delete yV;
      fout.close();
      std::cout << "file <" << outFName << "> created\n";
    }
  }

  if (1) {
    TString covOutFName="covForBlue-" + DYTools::analysisTag;
    if (test) covOutFName.Append("-yieldStatOnly");
    covOutFName.Append(".dat");
    std::ofstream foutCov(covOutFName);
    if (analysisIs2D) {
      for (int i=24; i<cov.GetNrows(); ++i) {
	for (int j=24; j<cov.GetNcols(); ++j) {
	  foutCov << Form("%3d %3d % 12.8e\n",i-23,j-23,cov(i,j));
	}
      }
    }
    else {
      for (int i=0; i<cov.GetNrows(); ++i) {
	for (int j=0; j<cov.GetNcols(); ++j) {
	  foutCov << Form("%3d %3d % 12.8e\n",i+1,j+1,cov(i,j));
	}
      }
    }
    foutCov.close();
    std::cout << "file <" << covOutFName << "> created\n";
  }

  return retCodeOk;
}


// ------------------------------------------------------------
// ------------------------------------------------------------

int compareError(const TH2D* h2CS, const TMatrixD &cov,
		 int relativeError,
		 TH2D* h2CSsystErr, int totalErr) {

  if (h2CSsystErr) {
    std::cout << "compareError: comparison of systematic error\n";
    if (!totalErr) std::cout << "NOTE: are you sure it is possible to compare"
			     << " incomplete (not the total) error\n";
  }
  std::cout << "compareError: relativeError=" << relativeError
	    << ", totalError=" << totalErr << "\n";

  // error from cross section calculation
  TH2D *h2errFromCS=NULL;
  TString errorKind;
  if (h2CSsystErr) {
    h2errFromCS=Clone(h2CSsystErr,"h2errFromCS_loc");
    if (totalErr) {
      h2errFromCS->Add(h2CS,1.);
      errorKind="tot.err ";
    }
    else errorKind="syst.err ";
  }
  else {
    h2errFromCS=Clone(h2CS,"h2errFromCS_loc");
    errorKind="stat.err ";
  }

  swapContentAndError(h2errFromCS);
  removeError(h2errFromCS);
  // error from covariance
  TH2D *h2errFromCov= errorFromCov(cov,"h2errFromCov_loc");

  TString yAxisLabel="propagated " + errorKind;

  if (relativeError) {
    if (!scaleHisto( h2errFromCS, h2CS, 0) ||
	!scaleHisto( h2errFromCov,h2CS, 0)) {
      return retCodeError;
    }
    yAxisLabel.Append(" (relative)");
  }

  std::vector<TH2D*> histoV;
  std::vector<TString> labelV;
  histoV.push_back(h2errFromCS);
  labelV.push_back(errorKind + TString("from CS"));
  histoV.push_back(h2errFromCov);
  labelV.push_back(errorKind + TString("from Cov"));

  std::vector<ComparisonPlot_t*> cpV;
  int delayDraw=0;

  TString canvName=Form("cx_r%d_t%d",relativeError,totalErr);
  TCanvas *cx= plotProfiles(canvName,histoV,labelV,NULL,1,
			    yAxisLabel,NULL,&cpV,delayDraw);
  if (!cx) return 0;
  return 1;
}


// ------------------------------------------------------------
// ------------------------------------------------------------
