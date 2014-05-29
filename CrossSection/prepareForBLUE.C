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

int prepareForBLUE(int analysisIs2D, int normalized=0) {
  if (analysisIs2D) {
    std::cout << "code is not ready of 2D case\n";
  }
  if (!DYTools::setup(analysisIs2D)) return retCodeError;

  // load normalized cross section
  NormCS_t sigmaZ(DYTools::_cs_None,"");
  if (normalized) {
    if (!sigmaZ.loadValues("dir-CSInfo/inpCS_sigmaZ.dat") ||
	(sigmaZ.cs()==double(0.))) {
      std::cout << "failed to load the cross section\n";
      return retCodeError;
    }
    if (!DYTools::checkCSKind(sigmaZ.csKind(),1,2,
			      DYTools::_cs_preFsr,//DYTools::_cs_postFsr,
			      DYTools::_cs_preFsrDet//,DYTools::_cs_postFsrDet
			      )) {
      std::cout << "prepareForBLUE: code is not ready for "
		<< sigmaZ.csKind() << "\n";
      return retCodeError;
    }
    std::cout << "sigmaZ=" << sigmaZ << "\n";
  }

  TString outputTag;
  switch (normalized) {
  case 0: outputTag="-absolute"; break;
  case 1: outputTag="-normalized"; break;
  case 2: outputTag="-rshape"; break;
  default:
    std::cout << "cannot interpret normalized=" << normalized << "\n";
    return retCodeError;
  }

  TH2D *h2CS=NULL, *h2CSsystErr=NULL;
  h2CS=loadMainCSResult(1,&h2CSsystErr);
  if (!h2CS) return retCodeError;

  if (0) {
    printHisto(h2CS);
    printHisto(h2CSsystErr);
  }

  int test=0;
  TString covFName=(!analysisIs2D) ? "finalCov-1D.root" : "finalCov-2D.root";
  if (test) covFName.ReplaceAll(".root","-yieldStatOnly.root");

  int nBins=DYTools::nUnfoldingBins;
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
  if (0) {
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

  TH2D *xsecFinal=Clone(h2CS,"xsecFinal");
  //TH2D *xsecFinalSystErr=Clone(h2CSsystErr,"xsecFinalSystErr");

  if (normalized) {
    xsecFinal->Reset();
    //xsecFinalSystErr->Reset();

    TMatrixD *yV= DYTools::getYBinLimits();
    //yV->Print();
    int idx=0;

    // keep the value of the uncertainty due to the norm.factor error
    double normRelErrSqr=pow(sigmaZ.csErrNoLumi()/sigmaZ.cs(),2);
    TMatrixD UncNorm(cov);
    UncNorm=normRelErrSqr;
    //UncNorm.Print();

    //std::cout << "sigmaZ.cs()=" << sigmaZ
    for (int ibin=1; ibin<=h2CS->GetNbinsX(); ++ibin) {
      for (int jbin=1; jbin<=h2CS->GetNbinsY(); ++jbin, idx++) {
	double massMin=h2CS->GetBinLowEdge(ibin);
	double massMax=h2CS->GetBinLowEdge(ibin) + h2CS->GetBinWidth(ibin);
	double yMin= (*yV)(ibin-1,jbin-1);
	double yMax= (*yV)(ibin-1,jbin);
	double dMdY= (analysisIs2D) ? (yMax-yMin) : (massMax-massMin);
	double factor=1./sigmaZ.cs();

	if (normalized==2) {
	  factor /= dMdY;
	}

	xsecFinal->SetBinContent(ibin,jbin,
				 factor * h2CS->GetBinContent(ibin,jbin));
	xsecFinal->SetBinError  (ibin,jbin,
				 factor * h2CS->GetBinError  (ibin,jbin));
	//xsecFinalSystErr->SetBinContent(ibin,jbin,
	//			 factor*h2CSsystErr->GetBinContent(ibin,jbin));
	//xsecFinalSystErr->SetBinError  (ibin,jbin,
	//			 factor*h2CSsystErr->GetBinError  (ibin,jbin));

	// normalize the covariance
	// the uncertainty due to the normalization will added later
	if (idx < cov.GetNrows()) {
	  double xsFactor= factor * h2CS->GetBinContent(ibin,jbin);
	  for (int ii=0; ii<cov.GetNrows(); ++ii) {
	    cov(ii,idx) *= factor;
	    cov(idx,ii) *= factor;

	    UncNorm(ii,idx) *= xsFactor;
	    UncNorm(idx,ii) *= xsFactor;
	  }
	}
      }
    }

    for (int ir=0; ir<UncNorm.GetNrows(); ++ir) {
      for (int ic=ir; ic<UncNorm.GetNcols(); ++ic) {
	double tmp=fabs(cov(ir,ic));
	if (tmp==double(0)) continue;
	if (UncNorm(ir,ic)/tmp > 0.01) {
	  std::cout << "ir=" << ir << ", ic=" << ic << ", UncNorm="
		    << UncNorm(ir,ic) << ", cov=" << cov(ir,ic) << "\n";
	}
      }
    }

    cov+= UncNorm;
    delete yV;
  }

  // error from cross section calculation
  TH2D *h2errFromCS=Clone(xsecFinal,"h2errFromCS");
  swapContentAndError(h2errFromCS);
  removeError(h2errFromCS);
  // error from covariance
  TH2D *h2errFromCov= errorFromCov(cov,"h2errFromCov");


  TString outPath="dir-forBlue/";
  gSystem->mkdir(outPath,1);

  // save xsec error
  if (1) {
    TString outFNameBase=outPath;
    outFNameBase.Append("xsecForBlue-" + DYTools::analysisTag + outputTag);
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
	  double xsecVal=xsecFinal->GetBinContent(ibin,jbin);
	  double xsecErr=h2err->GetBinContent(ibin,jbin);
	  double massMin=xsecFinal->GetBinLowEdge(ibin);
	  double massMax= massMin + xsecFinal->GetBinWidth(ibin);
	  double yMin= (*yV)(ibin-1,jbin-1);
	  double yMax= (*yV)(ibin-1,jbin);

	  if (analysisIs2D) {
	    if (ibin==1) continue;
	    if ((ibin==7) && (jbin>12)) break;
	    fout << Form("%3d  %3d  %6.1lf  %6.1lf  %5.2lf %5.2lf  %.16e   %.16e\n",
			 ibin-1,jbin,
			 massMin,massMax,
			 yMin,yMax,
			 xsecVal,xsecErr);
	  }
	  else {
	    fout << Form("%2d   %6.1lf   %6.1lf   %.16e   %.16e\n",
			 ibin,
			 massMin,massMax,
			 xsecVal,xsecErr);
	  }
	}
      }
      delete yV;
      fout.close();
      std::cout << "file <" << outFName << "> created\n";
    }
  }

  if (1) {
    TString covOutFName=outPath;
    covOutFName.Append("covForBlue-" + DYTools::analysisTag + outputTag);
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
