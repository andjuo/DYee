#ifndef covStudy_H
#define covStudy_H

inline
void performCovStudy(const TMatrixD &covInput)
{

  int dim=covInput.GetNrows();
  TMatrixD cov(covInput);
  TMatrixD corr(cov);
  for (int i=0; i<dim; i++) {
    for (int j=0; j<dim; j++) {
      corr(i,j) = cov(i,j)/sqrt(cov(i,i)*cov(j,j));
    }
  }

  if (1) {
    TCanvas *canvCov= new TCanvas("cCov","cCov",600,600);
    cov.Draw("COLZ");
    canvCov->Update();

    TCanvas *canvCorr= new TCanvas("cCorr","cCorr",600,600);
    corr.Draw("COLZ");
    canvCorr->Update();
  }

  //
  // Try to reconstruct the errors
  //

  int transpose=1;
  TDecompChol chol(cov);
  chol.Decompose();
  TMatrixD genM(chol.GetU());
  if (transpose) genM.Transpose(chol.GetU());

  std::cout << "check : 3x3 elements\n";
  std::cout << genM(0,0) << " " << genM(0,1) << " " << genM(0,2) << "\n";
  std::cout << genM(1,0) << " " << genM(1,1) << " " << genM(1,2) << "\n";
  std::cout << genM(2,0) << " " << genM(2,1) << " " << genM(2,2) << "\n";

  if (0) {
    TMatrixD U(genM);
    TMatrixD Ut(genM);
    Ut.Transpose(genM);

    TMatrixD testCov(U);
    if (transpose) {
      testCov = U*Ut;
    }
    else {
      testCov = Ut*U;
    }
    TCanvas *ct= new TCanvas("ct","ct",600,600);
    testCov.Draw("COLZ");
    ct->Update();
    return;
  }

  TCanvas *cTest= new TCanvas("cTest","cTest",600,600);

  // create the ensemble
  int nExps=1000;
  TVectorD avgX(dim);
  avgX.Zero();

  std::vector<TH1D*> h1RndVec;
  h1RndVec.reserve(nExps);

  double sumAt0=0;

  for (int iexp=0; iexp<nExps; ++iexp) {
    // create random uncorrelated distribution
    TVectorD rndX(dim), rndY(dim);
    for (int i=0; i<dim; i++) rndX[i]= gRandom->Gaus(0,1);

    rndY = genM * rndX;
    if (cTest) {
      TString name=Form("h_%d",iexp);
      TH1D* h= new TH1D(name,name,dim,0.,dim);
      h1RndVec.push_back(h);
      for (int i=0; i<dim; ++i) {
	h->SetBinContent(i+1, rndY[i]);
      }
      TString opt="hist";
      if (iexp) opt.Append("same");
      h->Draw(opt);
      cTest->Update();
    }

    rndY *= (1/double(nExps));
    avgX += rndY;
  }

  TH1D* hModel= new TH1D("hErrModel","hErrModel",dim,0.,dim);
  hModel->Reset();
  for (int i=0; i<dim; ++i) {
    hModel->SetBinContent(i+1, avgX[i]);
  }
  hModel->SetLineColor(kGreen+1);
  hModel->SetMarkerColor(kGreen+1);
  printHisto(hModel);

  if (cTest) {
    cTest->cd();
    hModel->Draw("hist same");
    cTest->Update();
  }

  // check the covariance and correlations
  if (1) {
    TMatrixD chkCov(dim,dim);
    chkCov.Zero();

    // Calculate
    // Cov(x,y) = avg(sum x*y) - avg(x) avg(y)
    // in two steps:
    // 1) the (sum x*y)
    for (unsigned int i=0; i<h1RndVec.size(); ++i) {
      const TH1D *h= h1RndVec[i];
      for (int ibin=1; ibin<=dim; ++ibin) {
	for (int jbin=1; jbin<=dim; ++jbin) {
	  double mult= h->GetBinContent(ibin) * h->GetBinContent(jbin);
	  chkCov(ibin-1,jbin-1) += mult;
	}
      }
    }
    chkCov*= 1/double(nExps);

    // 2) subtract the averages
    for (int i=0; i<dim; ++i) {
      for (int j=0; j<dim; ++j) {
	chkCov(i,j) -= hModel->GetBinContent(i+1)*hModel->GetBinContent(j+1);
      }
    }

    if (1) {
      std::cout << " i   cov(i,i)  chkCov(i,i)\n";
      for (int i=0; i<dim; ++i) {
	std::cout << i << "  " << cov(i,i) << " " << chkCov(i,i) << "\n";
      }
    }


    // Extract the errors
    for (int ibin=1; ibin<=dim; ++ibin) {
      hModel->SetBinError(ibin, sqrt(chkCov(ibin-1,ibin-1)));
    }


    TMatrixD chkCorr(chkCov);
    for (int i=0; i<dim; i++) {
      for (int j=0; j<dim; j++) {
	chkCorr(i,j) = chkCov(i,j)/sqrt(chkCov(i,i)*chkCov(j,j));
      }
    }

    TCanvas *canvChkCov= new TCanvas("cChkCov","cChkCov",600,600);
    chkCov.Draw("COLZ");
    canvChkCov->Update();

    TCanvas *canvChkCorr= new TCanvas("cChkCorr","cChkCorr",600,600);
    chkCorr.Draw("COLZ");
    canvChkCorr->Update();
  }
}

#endif
