#include "../Include/DYTools.hh"
#include "calcEventEffLink.h"
#include <TGraphAsymmErrors.h>
#include "../Include/ComparisonPlot.hh"

// ------------------------------------------------------------


// ------------------------------------------------------------

TGraphAsymmErrors* getAsymGraph(DYTools::TEtBinSet_t etBinning_inp,
				DYTools::TEtaBinSet_t etaBinning_inp,
				int iEta,
				const TMatrixD &Meff,
				const TMatrixD &MeffLo,
				const TMatrixD &MeffHi,
				TH1D  **histo=NULL,
				const char *histoName=NULL) {
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
  
  //const int bufsize=30;
  //char plotLabel[bufsize];
  int signedEta=DYTools::signedEtaBinning(etaBinning_inp);
  TString etaStr=Form("%s_%5.3lf_%5.3lf",(signedEta)?"#eta":"abs(#eta)",loc_etaBinLimits[iEta],loc_etaBinLimits[iEta+1]);
  TGraphAsymmErrors *gr=new TGraphAsymmErrors(loc_etBinCount,x,eff,dx,dx,effLo,effHi);
  gr->SetTitle(etaStr);
  return gr;
}



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

int loadEGammaEff(TString kindStr, TMatrixD **eff, TMatrixD **effLo, TMatrixD **effHi) {
  TString field=(kindStr.Index("mc")!=-1) ? "eff_mc" : "eff_data";
  TString fname="mediumID.root";
  TFile fin(fname,"read");
  if (!fin.IsOpen()) {
    std::cout << "Failed to open a file <" << fin.GetName() << ">\n";
    return 0;
  }
  TMatrixD *Meff=(TMatrixD*)fin.Get(field);
  TMatrixD *MeffLo=(TMatrixD*)fin.Get(field + TString("_errLo"));
  TMatrixD *MeffHi=(TMatrixD*)fin.Get(field + TString("_errHi"));
  fin.Close();
  *eff=Meff;
  *effLo=MeffLo;
  *effHi=MeffHi;
  std::cout << "eff="; Meff->Print();
  std::cout << "effLo="; MeffLo->Print();
  return 1;
}

// ------------------------------------------------------------
// does not work
/*
TGraphAsymmErrors* divideEffs(const TGraphAsymmErrors *gr1, const TGraphAsymmErrors *gr2) {
  TH1F *h1=gr1->GetHistogram();
  TH1F *h2=gr2->GetHistogram();
  std::cout << "h1="; h1->Print("range");
  std::cout << "h2="; h2->Print("range");
  TGraphAsymmErrors *div=(TGraphAsymmErrors*)gr1->Clone("temp");
  div->Divide(h1,h2,"cl=0.683 b(1,1) mode");
  return div;
}
*/

// ------------------------------------------------------------

TString effDataKindString(const TString str) {
  TString effKind="xx", dataKind="xx";
  if (str.Index("RECO")!=-1) effKind="RECO";
  else if (str.Index("ID")!=-1) effKind="ID";
  else if (str.Index("HLTleg1")!=-1) effKind="HLTleg1";
  else if (str.Index("HLTleg2")!=-1) effKind="HLTleg2";
  else if (str.Index("HLT")!=-1) effKind="HLT";
  if (str.Index("data")!=-1) dataKind="data";
  else if (str.Index("mc")!=-1) dataKind="mc";
  TString final=dataKind+TString(" ")+effKind;
  return final;
}

// ------------------------------------------------------------


void compareEff(int iBr=0, int iEta=0, 
		double transLegendY_user=0.,
		double ratioTitleOffset=0.58) {
  TString path1, path2;
  TString effKindLongStr1,effKindLongStr2;
  TString fnameBase="efficiency_TnP_1D_Full2012_";

  TString label1="new n-tuples";
  TString label2="old n-tuples";
  TString fnameTag;

  TString path3="";
  TString effKindLongStr3="";
  TString label3="";

  int HLTcomparison=0;
  //int divideBy1st=0;
  int relRatio=0;
  int HLTsystematics=0;
  double transLegendX=-0.2;
  double transLegendY=-0.4;


  if (0) {
    path1="/home/andriusj/cms/DYee8TeV-20140118/root_files/tag_and_probe/DY_j22_19712pb/";
    path2="/home/andriusj/cms/CMSSW_3_8_4/src/DYee8TeV-20130801/DrellYanDMDY/root_files/tag_and_probe/DY_j22_19789pb/";
    label1="new n-tuples";
    label2="old n-tuples";
    fnameTag="-new_vs_old--";
    effKindLongStr1="dataRECO_fit-fitEtBins6EtaBins5_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins6EtaBins5_PU";
    transLegendX=-0.05;
    transLegendY=-0.6;
  }

  if (0) {
    path1="/home/andriusj/cms/DYee-20131024/root_files_reg/tag_and_probe/DY_j22_19712pb/"; 
    path2="/home/andriusj/cms/DYee8TeV-20140118/root_files/tag_and_probe/DY_j22_19712pb/";
    effKindLongStr1="dataRECO_fit-fitEtBins6EtaBins5_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins6EtaBins5_PU";
    label1="DYee |#eta|<2.5";
    label2="DMDY |#eta|<2.5";
    fnameTag="-cmpPkg-Eta25--";
  }

  if (0) {
    //etaBinSet2=DYTools::ETABINS5corr;
    path1="/home/andriusj/cms/DYee8TeV-20140118/root_files/tag_and_probe/DY_j22_19712pb/";
    path2="/home/andriusj/cms/DYee8TeV-20140118-maxEta24/root_files/tag_and_probe/DY_j22_19712pb/";
    effKindLongStr1="dataRECO_fit-fitEtBins6EtaBins5_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins6EtaBins5_PU";
    label1="etaMax = 2.5";
    label2="etaMax = 2.4";
    fnameTag="-diffEtaMaxDMDY--";
  }

  if (0) {
    path1="/home/andriusj/cms/DYee-20131024/root_files_reg/tag_and_probe/DY_j22_19712pb/"; 
    path2=path1;

    effKindLongStr1="dataRECO_fit-fitEtBins6EtaBins5_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins6EtaBins5corr_PU";
    label1="DYee |#eta|<2.5";
    label2="DYee |#eta|<2.4";
    fnameTag="-diffEtaMax--";
  }

  if (0) {
    path1="/home/andriusj/cms/DYee-20131024/root_files_reg/tag_and_probe/DY_j22_19712pb/"; 
    path2="/home/andriusj/cms/DYee8TeV-20140118-maxEta24/root_files/tag_and_probe/DY_j22_19712pb/";
    effKindLongStr1="dataRECO_fit-fitEtBins6EtaBins5corr_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins6EtaBins5_PU";
    label1="DYee |#eta|<2.4";
    label2="DMDY |#eta|<2.4";
    fnameTag="-cmpPkg--";
  }

  if (0) {
    HLTcomparison=1;
    path1="/home/andriusj/cms/DYee-20131024/root_files_reg/tag_and_probe/DY_j22_19712pb/"; 
    path2=path1;
    path3=path1;
    effKindLongStr1="dataHLT_count-countEtBins6EtaBins5corr_PU";
    effKindLongStr2="dataHLTleg1_count-countEtBins6EtaBins5corr_PU";
    effKindLongStr3="dataHLTleg2_count-countEtBins6EtaBins5corr_PU";
    if (iBr==1) {
      effKindLongStr1.ReplaceAll("data","mc");
      effKindLongStr2.ReplaceAll("data","mc");
      effKindLongStr3.ReplaceAll("data","mc");
    }
    label1="HLT (7TeV analysis style)";
    label2="HLT leg1";
    label3="HLT leg2";
    fnameTag="-cmpHLT-etaMax24--";
  }

  if (0) {
    path1="/home/andriusj/cms/CMSSW_3_8_4/src/DYee8TeV-20130801/DrellYanDMDY/root_files/tag_and_probe/DY_j22_19789pb/";
    //path2="/media/spektras/DrellYanDMDY-Ilya-20130808/root_files/constants/DY_j22_19789pb/";
    path2="/home/andriusj/cms/DMDY-Ilya-20130808-my-copy/root_files/tag_and_probe/DY_j22_19789pb/";
    effKindLongStr1="dataRECO_fit-fitEtBins6EtaBins5_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins6EtaBins5_PU";
    label1="DYDM (Summer2013,""my"")";
    label2="DYDM (Ilya)";
    fnameTag="-checkSummer2013--";
  }

  if (0) {
    path1="/home/andriusj/cms/DYee-20131024/root_files_reg/tag_and_probe/DY_j22_19712pb/"; 
    path2="./";
    effKindLongStr1="dataID_fit-fitEtBins6EtaBins5corr_PU";
    effKindLongStr2="EGAMMA_dataID_EtBins6EtaBins5";
    label1="DYee |#eta|<2.4";
    label2="EGamma";
    if (1) {
      path3="/home/andriusj/cms/DYee-20131024/root_files_reg/tag_and_probe/DY_j22_19712pb/"; 
      effKindLongStr3="dataID_fit-fitEtBins6EtaBins5_PU";
      label3="DYee |#eta|<2.5";
    }
    fnameTag="-cmpEGamma--";
    transLegendX=-0.1;
  }

  if (0) { // tag-related HLT systematics
    HLTcomparison=1;
    HLTsystematics=1;
    //divideBy1st=1;
    relRatio=1;
    path1="../root_files_reg/tag_and_probe/DY_j22_19712pb/";
    path2="../root_files_reg/tag_and_probe/DY_j22_19712pb___tagLowerPt_tagLowerPt/";
    path3="../root_files_reg/tag_and_probe/DY_j22_19712pb___tagMediumID_tagMediumID/";
    effKindLongStr1="dataHLTleg1_count-countEtBins6systEtaBins5_PU";
    effKindLongStr2="dataHLTleg1_count-countEtBins6systEtaBins5_PU";
    effKindLongStr3="dataHLTleg1_count-countEtBins6systEtaBins5_PU";
    label1="Tight, pT>25GeV";
    label2="Tight, pT>20GeV";
    label3="Medium, pT>25GeV";
    fnameTag="-HLTsyst-tagLeg1--";
    if (1) {
      path1="../root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy/";
      path2="../root_files_reg/tag_and_probe/DY_j22_19712pb_UnregTagID/";
      path3="../root_files_reg/tag_and_probe/DY_j22_19712pb_UnregTagPt/";
      fnameTag="-unregHLTsyst-tagLeg1--";
    }
    transLegendX=-0.2;
  }

  if (1) { // PU-related HLT systematics
    HLTcomparison=1;
    HLTsystematics=1;
    //divideBy1st=1;
    relRatio=1;
    path1="../root_files_reg/tag_and_probe/DY_j22_19712pb/";
    path2="../root_files_reg/tag_and_probe/DY_j22_19712pb_Pileup5minus/";
    path3="../root_files_reg/tag_and_probe/DY_j22_19712pb_Pileup5plus/";
    effKindLongStr1="dataHLTleg1_count-countEtBins6systEtaBins5_PU";
    effKindLongStr2="dataHLTleg1_count-countEtBins6systEtaBins5_PU";
    effKindLongStr3="dataHLTleg1_count-countEtBins6systEtaBins5_PU";
    label1="reference";
    label2="PU -5%";
    label3="PU +5%";
    fnameTag="-HLTsyst-puLeg1--";
    if (1) { 
      path1="../root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy/";
      path2="../root_files_reg/tag_and_probe/DY_j22_19712pb_UnregEnPU5minus/";
      path3="../root_files_reg/tag_and_probe/DY_j22_19712pb_UnregEnPU5plus/";
      fnameTag="-unregHLTsyst-puLeg1--";
    }
    transLegendX=-0.02;
  }

  if (0) { // enReg-related HLT systematics
    HLTcomparison=1;
    HLTsystematics=1;
    //divideBy1st=1;
    relRatio=1;
    path1="../root_files_reg/tag_and_probe/DY_j22_19712pb/";
    path2="../root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy/";
    effKindLongStr1="dataHLTleg1_count-countEtBins6systEtaBins5_PU";
    effKindLongStr2="dataHLTleg1_count-countEtBins6systEtaBins5_PU";
    label1="reference";
    label2="unregr.energy";
    fnameTag="-HLTsyst-enRegLeg1--";
    transLegendX=-0.02;
  }

  if (0) { // RECO syst max|eta|<2.5
    HLTcomparison=0;
    //divideBy1st=1;
    relRatio=1;
    int loc_unregressed=1;
    path1="../root_files_reg/tag_and_probe/DY_j22_19712pb/";
    if (loc_unregressed) path1="../root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy/";
    path2=path1;
    effKindLongStr1="dataRECO_fit-fitEtBins6EtaBins5_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins6EtaBins5_max25_PU";
    label1="max|#eta|<2.4";
    label2="max|#eta|<2.5";

    if (0) {
      path3="/home/andriusj/cms/DYee8TeV-20140118/root_files/tag_and_probe/DY_j22_19712pb/";
      effKindLongStr3="dataRECO_fit-fitEtBins6EtaBins5_PU";
      label3="DMDY etaMax = 2.5";
    }
    if (0) {
      path3="/home/andriusj/cms/DYee8TeV-20140118-maxEta24/root_files/tag_and_probe/DY_j22_19712pb/";
     effKindLongStr3="dataRECO_fit-fitEtBins6EtaBins5_PU";
     label3="DMDY etaMax = 2.4";
   }
    fnameTag="-recoMaxEtaSyst--";
    if (loc_unregressed) fnameTag="-unReg-recoMaxEtaSyst--";
    transLegendX=0.0;
  }


  // -------------------------------
  // processing
  // -------------------------------

  if (transLegendY_user!=0.) transLegendY=transLegendY_user;

  if (label2 == TString("EGamma")) {
    if (iBr==0) {
    }
    else if (iBr==1) {
      effKindLongStr1.ReplaceAll("dataID_fit-fit","mcID_count-count");
      effKindLongStr2.ReplaceAll("dataID","mcID");
    }
    else {
      std::cout << "iBr error\n";
      return;
    }
  }
  else if (HLTsystematics) {
    if (iBr==0) {}
    else if (iBr==1) {
      effKindLongStr1.ReplaceAll("leg1","leg2");
      effKindLongStr2.ReplaceAll("leg1","leg2");
      effKindLongStr3.ReplaceAll("leg1","leg2");
      fnameTag.ReplaceAll("Leg1","Leg2");
    }
    else if (iBr==2) {
      effKindLongStr1.ReplaceAll("data","mc");
      effKindLongStr2.ReplaceAll("data","mc");
      effKindLongStr3.ReplaceAll("data","mc");
    }
    else if (iBr==3) { 
      effKindLongStr1.ReplaceAll("dataHLTleg1","mcHLTleg2");
      effKindLongStr2.ReplaceAll("dataHLTleg1","mcHLTleg2");
      effKindLongStr3.ReplaceAll("dataHLTleg1","mcHLTleg2");
      fnameTag.ReplaceAll("Leg1","Leg2");
   }
  }
  else {
    if (iBr==0) {
    }
    else if (iBr==1) {
      effKindLongStr1.ReplaceAll("RECO","ID");
      effKindLongStr2.ReplaceAll("RECO","ID");
    }
    else if (iBr==2) {
      effKindLongStr1.ReplaceAll("RECO_fit-fit","HLT_count-count");
      effKindLongStr2.ReplaceAll("RECO_fit-fit","HLT_count-count");
    }
    else if (iBr==3) {
      effKindLongStr1.ReplaceAll("dataRECO_fit-fit","mcRECO_count-count");
      effKindLongStr2.ReplaceAll("dataRECO_fit-fit","mcRECO_count-count");
      effKindLongStr3.ReplaceAll("dataRECO_fit-fit","mcRECO_count-count");
    }
    else if (iBr==4) {
      effKindLongStr1.ReplaceAll("dataRECO_fit-fit","mcID_count-count");
      effKindLongStr2.ReplaceAll("dataRECO_fit-fit","mcID_count-count");
    }
    else if (iBr==5) {
      effKindLongStr1.ReplaceAll("dataRECO_fit-fit","mcHLT_count-count");
      effKindLongStr2.ReplaceAll("dataRECO_fit-fit","mcHLT_count-count");
    }
    else if (iBr==6) {
      effKindLongStr1.ReplaceAll("dataRECO_fit-fit","dataHLTleg1_count-count");
      effKindLongStr2.ReplaceAll("dataRECO_fit-fit","dataHLTleg1_count-count");
    }
    else if (iBr==7) {
      effKindLongStr1.ReplaceAll("dataRECO_fit-fit","dataHLTleg2_count-count");
      effKindLongStr2.ReplaceAll("dataRECO_fit-fit","dataHLTleg2_count-count");
    }
    else if (iBr==8) {
      effKindLongStr1.ReplaceAll("dataRECO_fit-fit","mcHLTleg1_count-count");
      effKindLongStr2.ReplaceAll("dataRECO_fit-fit","mcHLTleg1_count-count");
    }
    else if (iBr==9) {
      effKindLongStr1.ReplaceAll("dataRECO_fit-fit","mcHLTleg2_count-count");
      effKindLongStr2.ReplaceAll("dataRECO_fit-fit","mcHLTleg2_count-count");
    }
    else {
      std::cout << "iBr error\n";
      return;
    }
  }


  TString fname1=path1 + fnameBase + effKindLongStr1 + TString(".root");
  TString fname2=path2 + fnameBase + effKindLongStr2 + TString(".root");

  DYTools::TEtBinSet_t etBinSet1=DetermineEtBinSet(effKindLongStr1);
  DYTools::TEtaBinSet_t etaBinSet1=DetermineEtaBinSet(effKindLongStr1);
  DYTools::TEtBinSet_t etBinSet2=DetermineEtBinSet(effKindLongStr2);
  DYTools::TEtaBinSet_t etaBinSet2=DetermineEtaBinSet(effKindLongStr2);
  std::cout << "sets: "<< EtBinSetName(etBinSet1) << "," << EtaBinSetName(etaBinSet1) << "  " << EtBinSetName(etBinSet2) << "," << EtaBinSetName(etaBinSet2) << "\n";

  TString effKind =effDataKindString(effKindLongStr1);
  TString effKind2=effDataKindString(effKindLongStr2);
  if (effKind != effKind2) {
    if ( !efficiencyIsHLT(DetermineEfficiencyKind(effKind )) ||
	 !efficiencyIsHLT(DetermineEfficiencyKind(effKind2)) ) {
      std::cout << "effKind1=<" << effKind << ">\n";
      std::cout << "effKind2=<" << effDataKindString(effKindLongStr2) << ">\n";
      return;
    }
  }

  TString dataKind=effKind + TString(" ");
  int weighted1=(effKindLongStr1.Index("count-count")!=-1) ? 1 : 0;
  int weighted2=(effKindLongStr2.Index("count-count")!=-1) ? 1 : 0;
  //int iEta=0;

  TMatrixD *eff1=NULL, *eff1ErrLo=NULL, *eff1ErrHi=NULL;
  TMatrixD *eff2=NULL, *eff2ErrLo=NULL, *eff2ErrHi=NULL;
  TH1D *histo1=NULL, *histo2=NULL;
  const char *histo1Name="histo1";
  const char *histo2Name="histo2";


  if (!loadEff(fname1,weighted1,&eff1,&eff1ErrLo,&eff1ErrHi)) {
    std::cout << "failed to get fields from <" << fname1 << "> (1)\n"; 
    return ;
  }

  if (label2 == TString("EGamma")) {
    if (!loadEGammaEff(effKindLongStr2,&eff2,&eff2ErrLo,&eff2ErrHi)) {
      std::cout << "failed to load EGammaEff\n";
      return;
    }
  }
  else {
    if (!loadEff(fname2,weighted2,&eff2,&eff2ErrLo,&eff2ErrHi)) {
      std::cout << "failed to get fields from <" << fname2 << "> (2)\n"; 
      return ;
    }
  }

  TGraphAsymmErrors* gr1=getAsymGraph(etBinSet1,etaBinSet1,iEta,*eff1,*eff1ErrLo,*eff1ErrHi,&histo1,histo1Name);
  //std::cout << gr1->GetTitle() << ": "; gr1->Print("range");

  TGraphAsymmErrors* gr2=getAsymGraph(etBinSet2,etaBinSet2,iEta,*eff2,*eff2ErrLo,*eff2ErrHi,&histo2,histo2Name);
  //std::cout << gr2->GetTitle() << ": "; gr2->Print("range");

  TGraphAsymmErrors* div=(TGraphAsymmErrors*)gr1->Clone("div");
  //TH1D *div=(TH1D*)histo1->Clone("div");
  //div->Divide(histo1,histo2,1.,1.,"b");
  if (relRatio) {
    TH1D *diff12=(TH1D*)histo1->Clone("diff12");
    diff12->Add(histo2,-1.);
    div->Divide(diff12,histo1,"pois");
  }
  else div->Divide(histo1,histo2,"pois");
  div->Print("range");

  TH1D *histo3=NULL;
  TGraphAsymmErrors* gr3=NULL;
  //TH1D* div31=NULL;
  TGraphAsymmErrors* div31=NULL;

  if (effKindLongStr3.Length() && label3.Length()) {
    TString fname3= path3 + fnameBase + effKindLongStr3 + TString(".root");
    int weighted3=(effKindLongStr3.Index("count-count")!=-1) ? 1:0;
    TMatrixD *eff3=NULL, *eff3ErrLo=NULL, *eff3ErrHi=NULL;
    if (!loadEff(fname3,weighted3,&eff3,&eff3ErrLo,&eff3ErrHi)) {
      std::cout << "failed to get field from <" << fname3 << "> (3)\n";
      return ;
    }
    DYTools::TEtBinSet_t etBinSet3=DetermineEtBinSet(effKindLongStr3);
    DYTools::TEtaBinSet_t etaBinSet3=DetermineEtaBinSet(effKindLongStr3);
   
    TString histo3Name="histo3";
    gr3=getAsymGraph(etBinSet3,etaBinSet3,iEta,*eff3,*eff3ErrLo,*eff3ErrHi, &histo3,histo3Name);
    //std::cout << gr3->GetTitle() << ": "; gr3->Print("range");
    delete eff3;
    delete eff3ErrLo;
    delete eff3ErrHi;

    //div31=(TH1D*)histo1->Clone("div31");
    div31=(TGraphAsymmErrors*)gr1->Clone("div31");
    div31->SetTitle("div31");
    //if (HLTcomparison) div31->Divide(histo1,histo3,1.,1.,"b");
    //else div31->Divide(histo2,histo3,1.,1.,"b");
    if (HLTcomparison) {
      if (relRatio) {
	TH1D *diff13=(TH1D*)histo1->Clone("diff13");
	diff13->Add(histo3,-1.0);
	div31->Divide(diff13,histo1,"pois");
      }
      else div31->Divide(histo1,histo3,"pois");
    }
    else div31->Divide(histo2,histo3,"pois");
    div31->Print("range");
    div31->SetLineColor(kGreen+1);
    div31->SetMarkerColor(kGreen+1);

    //div->Divide(histo2,histo1,1.,1.,"b");
    if (!relRatio) div->Divide(histo2,histo1,"pois");
    div->SetLineColor(kBlue);
    div->SetMarkerColor(kBlue);
  }

  double *loc_etaBinLimits=DYTools::getEtaBinLimits(etaBinSet1);
  int signedEta=DYTools::signedEtaBinning(etaBinSet1);
  TString cpTitle=dataKind+ TString(Form(" %5.3lf #leq %s #leq %5.3lf",loc_etaBinLimits[iEta],(signedEta)?"#eta":"abs(#eta)",loc_etaBinLimits[iEta+1]));

  //CPlot cp("comp",cpTitle,"E_{T} [GeV]","efficiency");
  ComparisonPlot_t cp(ComparisonPlot_t::_ratioRel,"comp",cpTitle,
			  "E_{T} [GeV]","efficiency","ratio");
  cp.SetRefIdx(-111); // do not plot lower panel
  cp.SetLogx();

  if (gr3 && HLTcomparison) { // for HLT efficiency
    cp.SetYRange(0.0,1.02);
    /*
    if (DetermineDataKind(effKind)==DYTools::DATA) {
      if (iEta==2) cp.SetYRange(0.3,1.01);
      else if (iEta==1) cp.SetYRange(0.3,1.01);
      else if (iEta==0) cp.SetYRange(0.3,1.01);
    }
    else {
      if (iEta==0) cp.SetYRange(0.5,1.01);
      else if (iEta==1) cp.SetYRange(0.5,1.01);
      else if (iEta==2) cp.SetYRange(0.5,1.01);
      else if (iEta==3) cp.SetYRange(0.5,1.01);
      else if (iEta==4) cp.SetYRange(0.5,1.01);
    }
    */
  }

  TCanvas *cx=new TCanvas("cx","cx",600,700);
  cp.Prepare2Pads(cx);

  gr1->GetYaxis()->SetTitleOffset(1.4);

  if (gr3 && !HLTcomparison) {
    std::cout << "\n\tInverted plotting order 2,1\n";
    gr2->GetYaxis()->SetTitleOffset(1.2);
    //gr1->SetMarkerStyle(24);
    div->SetMarkerStyle(24);
    cp.AddGraph(gr2,label2,"LPE1",kBlue);
    cp.AddGraph(gr1,label1," PE1",kBlack,24);
  }
  else {
    cp.AddGraph(gr1,label1,"LPE1",kBlack);
    cp.AddGraph(gr2,label2,"LPE1",kBlue,24);
  }
  if (gr3) {
    //gr3->SetMarkerStyle(27);
    div31->SetMarkerStyle(27);
    cp.AddGraph(gr3,label3," PE1",kGreen+1,27);
  }
  cp.Draw(cx,0,"png",1);
  cp.TransLegend(transLegendX, transLegendY);
  cp.WidenLegend(0.2,0.);
  cx->cd(2);
  cx->GetPad(2)->SetLogx(cp.fLogx);
  div->SetTitle("");
  div->GetXaxis()->SetTitle("E_{T}");
  div->GetXaxis()->SetTitleSize(0.17);
  div->GetXaxis()->SetLabelSize(0.15);
  div->GetXaxis()->SetNoExponent();
  div->GetXaxis()->SetMoreLogLabels();

  div->GetYaxis()->SetTitle("ratio");
  if (relRatio) div->GetYaxis()->SetTitle("([1]-[i])/[1]");
  div->GetYaxis()->SetTitleOffset(0.5);
  if (label2 == TString("EGamma")) {
    div->GetYaxis()->SetTitle("EG/our");
    div->GetYaxis()->SetTitleOffset(0.45);
  }
  if (ratioTitleOffset>0) div->GetYaxis()->SetTitleOffset(ratioTitleOffset);
  div->GetYaxis()->SetTitleSize(0.13);
  div->GetYaxis()->SetLabelSize(0.13);
  div->GetYaxis()->SetNdivisions(805);  
  if (div31) {
    if (HLTcomparison) {
      if (relRatio) {
	div->GetYaxis()->SetRangeUser(-0.1,0.1);
      }
      else {
	div->GetYaxis()->SetRangeUser(0.99,1.01);
	if (iEta==2) div->GetYaxis()->SetRangeUser(0.9,1.1);
      }
    }
  }
  div->Draw("ALPE");
  if (div31) {
    div31->Draw("LPE1 same");
  }

  double one=(relRatio) ? 0. : 1.;
  TLine *lineAtOne =   new TLine(10,one, 500,one);
  lineAtOne->SetLineStyle(kDashed);
  lineAtOne->SetLineWidth(1);
  lineAtOne->SetLineColor(kBlack);
  lineAtOne->Draw();

  cx->Update();

  if (fnameTag.Length()) {
    TString fname=TString("fig-eff-") + fnameTag + cpTitle;
    fname.ReplaceAll(" #leq "," ");
    fname.ReplaceAll(" ","_");
    fname.ReplaceAll("(#eta)","Eta");
    fname.ReplaceAll("#eta","eta");
    fname.ReplaceAll(".","_");
    //fname.Append(".png");
    std::cout << "fname=" << fname << "\n";

    TString locOutDir=TString("plots") + fnameTag;
    locOutDir.ReplaceAll("--","");
    SaveCanvas(cx,fname,locOutDir);
  }

  return ;
}
