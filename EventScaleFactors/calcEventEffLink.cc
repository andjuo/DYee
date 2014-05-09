#include "calcEventEffLink.h"

#include "calcEventEff.C"

// ------------------------------------------------------------

//int fillEfficiencyConstants( const InputFileMgr_t &inpMgr,
//			     DYTools::TSystematicsStudy_t systMode ) {
//  return fillEfficiencyConstants(inpMgr,systMode);
//}

// ------------------------------------------------------------

TGraphAsymmErrors* getAsymGraph_vsEt(DYTools::TEtBinSet_t etBinning_inp,
				     DYTools::TEtaBinSet_t etaBinning_inp,
				     int iEta,
				     const TMatrixD &Meff,
				     const TMatrixD &MeffLo,
				     const TMatrixD &MeffHi,
				     TH1D  **histo,
				     const char *histoName) {
  int loc_etBinCount=DYTools::getNEtBins(etBinning_inp);
  double *loc_etBinLimits=DYTools::getEtBinLimits(etBinning_inp);
  //int loc_etaBinCount=DYTools::getNEtaBins(etaBinning_inp);
  double *loc_etaBinLimits=DYTools::getEtaBinLimits(etaBinning_inp);

  double x[loc_etBinCount], dx[loc_etBinCount];
  double eff[loc_etBinCount], effLo[loc_etBinCount], effHi[loc_etBinCount];
  for (int i=0; i<loc_etBinCount; ++i) {
    x[i] = 0.5*(loc_etBinLimits[i  ] + loc_etBinLimits[i+1]);
    dx[i]= 0.5*(loc_etBinLimits[i+1] - loc_etBinLimits[i  ]);
    eff[i] = Meff[i][iEta];
    effLo[i] = MeffLo[i][iEta];
    effHi[i] = MeffHi[i][iEta];
  }

  if ((histo!=NULL) && (histoName!=NULL)) {
    TH1D *h=new TH1D(histoName,histoName,loc_etBinCount,loc_etBinLimits);
    h->SetDirectory(0);
    h->SetStats(0);
    h->GetXaxis()->SetTitle("E_{T}");
    h->GetYaxis()->SetTitle("efficiency");
    for (int i=0; i<loc_etBinCount; ++i) {
      h->SetBinContent(i+1, eff[i]);
      h->SetBinError  (i+1, 0.5*(effLo[i]+effHi[i]));
    }
    *histo=h;
  }
  
  int signedEta=DYTools::signedEtaBinning(etaBinning_inp);
  TString etaStr=Form("%s_%5.3lf_%5.3lf",(signedEta)?"#eta":"abs(#eta)",loc_etaBinLimits[iEta],loc_etaBinLimits[iEta+1]);
  TGraphAsymmErrors *gr=new TGraphAsymmErrors(loc_etBinCount,x,eff,dx,dx,effLo,effHi);
  gr->SetTitle(etaStr);
  return gr;
}

// ------------------------------------------------------------

TGraphAsymmErrors* getAsymGraph_vsEta(DYTools::TEtBinSet_t etBinning_inp,
				      DYTools::TEtaBinSet_t etaBinning_inp,
				      int iEt,
				      const TMatrixD &Meff,
				      const TMatrixD &MeffLo,
				      const TMatrixD &MeffHi,
				      TH1D  **histo,
				      const char *histoName) {
  //int loc_etBinCount=DYTools::getNEtBins(etBinning_inp);
  double *loc_etBinLimits=DYTools::getEtBinLimits(etBinning_inp);
  int loc_etaBinCount=DYTools::getNEtaBins(etaBinning_inp);
  double *loc_etaBinLimits=DYTools::getEtaBinLimits(etaBinning_inp);

  double x[loc_etaBinCount], dx[loc_etaBinCount];
  double eff[loc_etaBinCount], effLo[loc_etaBinCount], effHi[loc_etaBinCount];
  for (int i=0; i<loc_etaBinCount; ++i) {
    x[i] = 0.5*(loc_etaBinLimits[i  ] + loc_etaBinLimits[i+1]);
    dx[i]= 0.5*(loc_etaBinLimits[i+1] - loc_etaBinLimits[i  ]);
    eff[i]   = Meff[iEt][i];
    effLo[i] = MeffLo[iEt][i];
    effHi[i] = MeffHi[iEt][i];
  }

  int signedEta=DYTools::signedEtaBinning(etaBinning_inp);
  
  if ((histo!=NULL) && (histoName!=NULL)) {
    TH1D *h=new TH1D(histoName,histoName,loc_etaBinCount,loc_etaBinLimits);
    h->SetDirectory(0);
    h->SetStats(0);
    TString xaxisTitle= (signedEta) ? "|#eta|" : "#eta";
    h->GetXaxis()->SetTitle(xaxisTitle);
    h->GetYaxis()->SetTitle("efficiency");
    for (int i=0; i<loc_etaBinCount; ++i) {
      h->SetBinContent(i+1, eff[i]);
      h->SetBinError  (i+1, 0.5*(effLo[i]+effHi[i]));
    }
    *histo=h;
  }
  
  TString titleStr=Form("Et_%1.0lf_%1.0lf",loc_etBinLimits[iEt],loc_etBinLimits[iEt+1]);
  TGraphAsymmErrors *gr=new TGraphAsymmErrors(loc_etaBinCount,x,eff,dx,dx,effLo,effHi);
  gr->SetTitle(titleStr);
  return gr;
}


// ------------------------------------------------------------

TGraphAsymmErrors* addErrors(const TGraphAsymmErrors *gr, const TVectorD &errs,
			     TString newName, int relative) {
  if (errs.GetNoElements()!=gr->GetN()) {
    std::cout << "addErrors: size mismatch gr->getN=" << gr->GetN()
	      << ", errs.GetNoElements=" << errs.GetNoElements() << "\n";
    return NULL;
  }
  TGraphAsymmErrors* gr2=(TGraphAsymmErrors*)gr->Clone(newName);
  for (int ibin=0; ibin<gr->GetN(); ++ibin) {
    double ylo=gr->GetErrorYlow(ibin);
    double yhi=gr->GetErrorYhigh(ibin);
    double addErr=errs(ibin);
    if (relative) addErr *= gr->GetY()[ibin];
    double newLo=sqrt(pow(ylo,2) + pow(addErr,2));
    double newHi=sqrt(pow(yhi,2) + pow(addErr,2));
    gr2->SetPointEYlow(ibin,newLo);
    gr2->SetPointEYhigh(ibin,newHi);
    //std::cout << " ibin=" << ibin << ": " << gr->GetY()[ibin] << " + "
    //      << yhi << " - " << ylo << "\n";
  }
  return gr2;
}

// ------------------------------------------------------------
// ------------------------------------------------------------

int loadEff(const TString &fname, int weighted, TMatrixD **eff, TMatrixD **effLo, TMatrixD **effHi) {
  std::cout  << "loading <" << fname << ">\n";
  TFile file1(fname,"read");
  if (!file1.IsOpen()) { std::cout << "failed to open file <" << file1.GetName() << ">\n"; return 0; }
  TString field="effArray2D";
  TString fieldErrLo="effArrayErrLow2D";
  TString fieldErrHi="effArrayErrHigh2D";
  if (weighted) {
    field.Append("Weighted");
    fieldErrLo.Append("Weighted");
    fieldErrHi.Append("Weighted");
  }
  TMatrixD *eff1= (TMatrixD*)file1.Get(field);
  TMatrixD *eff1ErrLo= (TMatrixD*)file1.Get(fieldErrLo);
  TMatrixD *eff1ErrHi= (TMatrixD*)file1.Get(fieldErrHi);
  file1.Close();
  if (!eff1 || !eff1ErrLo || !eff1ErrHi) {
    std::cout << "failed to get fields from <" << file1.GetName() << ">\n"; 
    return 0;
  }
  *eff=eff1;
  *effLo=eff1ErrLo;
  *effHi=eff1ErrHi;
  return 1;
}

// ------------------------------------------------------------

int loadEffVsPU(const TString &fnameBase, int weighted,
		std::vector<TMatrixD*> &effV,
		std::vector<TMatrixD*> &effLoV, std::vector<TMatrixD*> &effHiV,
		std::vector<TString> &labelsV)
{
  ClearVec(effV); ClearVec(effLoV); ClearVec(effHiV);
  labelsV.clear();


  effV.reserve(DYTools::nPVBinCount);
  effLoV.reserve(DYTools::nPVBinCount);
  effHiV.reserve(DYTools::nPVBinCount);
  labelsV.reserve(DYTools::nPVBinCount);

  TMatrixD *eff=NULL, *effLo=NULL, *effHi=NULL;
  for (int pu_i=0; pu_i<DYTools::nPVBinCount; ++pu_i) {
    UInt_t pvMin=UInt_t(DYTools::nPVLimits[pu_i  ]+0.6);
    UInt_t pvMax=UInt_t(DYTools::nPVLimits[pu_i+1]-0.4);
    TString puRangeStr=Form("_%u_%u",pvMin,pvMax);
    labelsV.push_back(TString("eff") + puRangeStr);
    TString resRootFName= fnameBase + puRangeStr + TString(".root");
    //std::cout << "loading <" << resRootFName << ">\n";
    if (!loadEff(resRootFName,weighted,&eff,&effLo,&effHi)) {
      std::cout << " .. error in loadEffVsPU\n";
      return 0;
    }
    effV.push_back(eff);
    effLoV.push_back(effLo);
    effHiV.push_back(effHi);
  }
  return 1;
}


// ------------------------------------------------------------

TGraphAsymmErrors* effVsPU_asGraph(int iEtBin, int iEtaBin,
				   const std::vector<TMatrixD*> &effV,
				   const std::vector<TMatrixD*> &effLoV,
				   const std::vector<TMatrixD*> &effHiV)
{
  if ((effV.size()!=int(DYTools::nPVBinCount)) ||
      (effLoV.size()!=effV.size()) ||
      (effHiV.size()!=effV.size())) {
    std::cout << "effVsPU_asGraph: vec size mismatch:\n";
    std::cout << "  - Expected size: " << DYTools::nPVBinCount << "\n";
    std::cout << "  - effV[" << effV.size() << "], effLoV["
	      << effLoV.size() << "], effHiV[" << effHiV.size() << "]\n";
    return NULL;
  }

  double x[DYTools::nPVBinCount], dx[DYTools::nPVBinCount];
  double eff[DYTools::nPVBinCount];
  double effLo[DYTools::nPVBinCount], effHi[DYTools::nPVBinCount];

  for (int pu_i=0; pu_i<DYTools::nPVBinCount; ++pu_i) {
    x [pu_i]= 0.5*( DYTools::nPVLimits[pu_i] + DYTools::nPVLimits[pu_i+1] );
    dx[pu_i]= 0.5*( DYTools::nPVLimits[pu_i+1] - DYTools::nPVLimits[pu_i] );
    eff  [pu_i] = (*effV[pu_i])(iEtBin,iEtaBin);
    effLo[pu_i] = (*effLoV[pu_i])(iEtBin,iEtaBin);
    effHi[pu_i] = (*effHiV[pu_i])(iEtBin,iEtaBin);
  }

  TGraphAsymmErrors *gr= new TGraphAsymmErrors(DYTools::nPVBinCount,
					       x,eff,dx,dx,effLo,effHi);
  if (!gr) {
    std::cout << "effVsPU_asGraph: failed to create TGraphAsymmErrors\n";
  }
  return gr;
}


// ------------------------------------------------------------
// ------------------------------------------------------------
// ------------------------------------------------------------
