#include "calcEventEffLink.h"

#include "calcEventEff.C"

// ------------------------------------------------------------

int fillEfficiencyConstants( const InputFileMgr_t &inpMgr ) {
  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;
  return fillEfficiencyConstants(inpMgr,systMode);
}

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
