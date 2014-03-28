#include <TROOT.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TString.h>
#include <TFile.h>
#include <TH2D.h>
#include <iostream>
#include <vector>
#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"

void extractRndESF(int fileIdx) {
  TString fname=(fileIdx==0) ? "rhoFileSF_nMB41_egamma_asymHLT_Unregressed_energy-allSyst_100.root" : "rhoFileSF_nMB7_egamma_asymHLT_Unregressed_energy-allSyst_100.root";
  TString fname2=fname;
  fname2.ReplaceAll("rhoFileSF_","covRhoFileSF_");

  int nExps=100;
  int nTotBins=(fileIdx==0) ? 41 : 156;

  int imMax=(fileIdx==0) ? DYTools::_nMassBins2012 : DYTools::_nMassBins2D;
  int iyMax=(fileIdx==0) ? 1 : 24;
  const double *massEdges=(fileIdx==0) ? DYTools::_massBinLimits2012 : DYTools::_massBinLimits2D;
  double yMax=(fileIdx==0) ? 8.9 : 2.4;

  TVectorD sumWeight(nTotBins);
  TMatrixD sumWeightRho_Rnd(nExps,nTotBins);
  TMatrixD sumWeightRhoSqr_Rnd(nExps,nTotBins);
  TMatrixD scaleFactor(imMax,iyMax);
  TMatrixD scaleFactorErr(imMax,iyMax);

  std::vector<TH2D*> rndSF;  
  

  int res=1;

  // Load randomized scale factor matrices
  if (res) {
    TFile fin(fname,"read");
    if (!fin.IsOpen()) {
      std::cout << "failed to open file <" << fin.GetName() << ">\n";
      return;
    }
    TVectorD* V=(TVectorD*)fin.Get("sumWeight");
    if (!V) {
      std::cout << "failed to get sumWeight\n";
      return;
    }
    sumWeight = *V;
    TMatrixD* M=(TMatrixD*)fin.Get("sumWeightRho_Rnd");
    if (!M) {
      std::cout << "failed to get sumWeightRho_Rnd\n";
      return;
    }
    sumWeightRho_Rnd = *M;
    TMatrixD* M2=(TMatrixD*)fin.Get("sumWeightRhoSqr_Rnd");
    if (!M2) {
      std::cout << "failed to get sumWeightRhoSqr_Rnd\n";
      return;
    }
    sumWeightRhoSqr_Rnd = *M2;

    delete V;
    delete M;
    delete M2;

    fin.Close();
    std::cout << "data from <" << fin.GetName() << "> loaded\n";
  }

  // Load scale factors
  if (res) {
    TFile fin(fname2,"read");
    if (!fin.IsOpen()) {
      std::cout << "failed to open file <" << fin.GetName() << ">\n";
      return;
    }
    TMatrixD* M=(TMatrixD*)fin.Get("scaleFactor");
    if (!M) {
      std::cout << "failed to get scaleFactor\n";
      return;
    }
    scaleFactor = *M;
    TMatrixD* M2=(TMatrixD*)fin.Get("scaleFactorErr");
    if (!M2) {
      std::cout << "failed to get scaleFactorErr\n";
      return;
    }
    scaleFactorErr = *M2;

    delete M;
    delete M2;

    fin.Close();
    std::cout << "data from <" << fin.GetName() << "> loaded\n";
  }


  // convert rnd matrix to the TH2D distributions

  for (int iexp=0; iexp<nExps; ++iexp) {
    TH2D* h2=NULL;
    TString histoName=Form("h2RndScaleFactor_%dD_%d",fileIdx+1,iexp);
    TString xTitle;
    h2= new TH2D(histoName,histoName,
		 imMax, massEdges,
		 iyMax, 0., yMax);
    if (fileIdx==0) {
      xTitle="#it{M}_{ee} [GeV]";
    }
    else {
      xTitle="|y|";
    }
    rndSF.push_back(h2);
    h2->SetDirectory(0);
    h2->SetStats(0);
    h2->Sumw2();
    h2->GetXaxis()->SetTitle(xTitle);
    h2->GetYaxis()->SetTitle("rnd scale factor");

    
    int idx=0;
    for (int im=0; im<imMax; ++im) {
      for (int iy=0; iy<iyMax; ++iy, idx++) {
	if (idx>=nTotBins) break;
	double val=sumWeightRho_Rnd(iexp,idx) / sumWeight(idx);
	double valErrSqr= sumWeightRhoSqr_Rnd(iexp,idx)/sumWeight(idx) 
	  - val*val;
	h2->SetBinContent(im+1,iy+1, val);
	h2->SetBinError  (im+1,iy+1, sqrt(valErrSqr));
      }
    }
  }

  TVectorD mass(imMax+1);
  for (int im=0; im<=imMax; ++im) mass[im]=massEdges[im];
  TVectorD rapidityCounts(imMax);
  if (fileIdx==0) for (int im=0; im<imMax; ++im) rapidityCounts[im]=iyMax;
  else {
    for (int im=0; im<imMax-1; ++im) rapidityCounts[im]=24;
    rapidityCounts[imMax-1]=12;
  }

  TString outFName=fname;
  outFName.ReplaceAll(".root","-histos.root");
  TFile outF(outFName,"recreate");
  scaleFactor.Write("scaleFactor");
  scaleFactorErr.Write("scaleFactorErr");
  saveVec(outF,rndSF,"rndSF");
  outF.cd();
  mass.Write("massBinLimits");
  rapidityCounts.Write("rapidityCounts");
  outF.Close();
  std::cout << "file <" << outF.GetName() << "> created\n";
  

  return;
}
