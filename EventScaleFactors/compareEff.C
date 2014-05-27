//
// Made on March 14, 2014
// Enhanced version of compareEff.C. Allows plotting vs eta
//

#include "../Include/DYTools.hh"
#include "calcEventEffLink.h"
#include <TGraphAsymmErrors.h>
#include "../Include/ComparisonPlot.hh"

// ------------------------------------------------------------


// ------------------------------------------------------------
// ------------------------------------------------------------
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


void compareEff(int iBr=0, int iBin=0, int vsEt=1,
		int doSave=0,
		double transLegendY_user=0.,
		double transLegendX_user=0.,
		TString *outFileName_ptr=NULL,
		TString *outDir_ptr=NULL) {

  double ratioTitleOffset=0.58;

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
  int relRatio=0;   // if relRatio==-1, the ratio plot is not displayed
  int HLTsystematics=0;
  double transLegendX=-0.2;
  double transLegendY=-0.4;
  int allowInvert12=1; // added on Mar 19, 2014
  int setId=0;

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
    effKindLongStr1="dataID_fit-fitEtBins6EtaBins5_PU";
    effKindLongStr2="EGAMMA_dataID_EtBins6EtaBins5";
    label1="DYee |#eta|<2.4";
    label2="EGamma";
    if (0) {
      path3="/home/andriusj/cms/DYee-20131024/root_files_reg/tag_and_probe/DY_j22_19712pb/"; 
      effKindLongStr3="dataID_fit-fitEtBins6EtaBins5_maxEta25_PU";
      label3="DYee |#eta|<2.5";
    }
    fnameTag="-cmpEGamma2--";
    transLegendX=(vsEt) ? -0.1 : -0.4;
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

  if (0) { // PU-related HLT systematics
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

  if (0) { // 2014.03.15
    HLTcomparison=0;
    relRatio=-1;
    path1="dir-Rami/";
    path2=path1; path3=path1;
    effKindLongStr1="dataRECO_fit-fitEtBins6EtaBins5_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins6EtaBins7_PU";
    effKindLongStr3="dataRECO_fit-fitEtBins6EtaBins9_PU";
    label1="Rami: Et6Eta5";
    label2="Rami: Et6Eta7";
    label3="Rami: Et6Eta9";
    fnameTag="-Rami-vars--";
    transLegendX=-0.4;
  }

  if (0) { // 2014.05.10
    if (vsEt) {
      std::cout << "Eta bins are different\n";
      return;
    }
    HLTcomparison=0;
    relRatio=-1;
    path1="Results-DYee-effBinStudy/root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy/";
    path2=path1; path3=path1;
    effKindLongStr1="dataRECO_fit-fitEtBins6systEtaBins5_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins6systEtaBins7_PU";
    effKindLongStr3="dataRECO_fit-fitEtBins6systEtaBins9_PU";
    label1="Et6Eta5";
    label2="Et6Eta7";
    label3="Et6Eta9";
    fnameTag="-et6-vars--";
    transLegendX=-0.4;
    allowInvert12=0;
  }

  if (0) { // 2014.05.13
    HLTcomparison=0;
    relRatio=-1;
    path1="Results-DYee-effBinStudy/root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy/";
    path2=path1; path3=path1;
    effKindLongStr1="dataRECO_fit-fitEtBins6systEtaBins5_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins6systEtaBins5nomerge_PU";
    //effKindLongStr3="dataRECO_fit-fitEtBins6systEtaBins9_PU";
    label1="Et6Eta5";
    label2="Et6Eta5nomerge";
    label3="Et6Eta9";
    fnameTag="-et6eta5nomerge-vars--";
    transLegendX=-0.4;
    allowInvert12=0;
  }

  if (0) { // 2014.05.13
    if (vsEt) {
      std::cout << "Eta bins are different\n";
      return;
    }
    HLTcomparison=0;
    relRatio=-1;
    path1="Results-DYee-effBinStudy/root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy/";
    path2=path1; path3=path1;
    effKindLongStr1="dataRECO_fit-fitEtBins6systEtaBins5nomerge_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins6systEtaBins7_PU";
    effKindLongStr3="dataRECO_fit-fitEtBins6systEtaBins9_PU";
    label1="Et6Eta5nomerge";
    label2="Et6Eta7";
    label3="Et6Eta9";
    fnameTag="-et6eta5nomerge-vars--";
    transLegendX=-0.4;
    allowInvert12=0;
  }

  if (0) { // 2014.05.10
    HLTcomparison=0;
    relRatio=-1;
    path1="Results-DYee-effBinStudy/root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy/";
    path2=path1;
    effKindLongStr1="dataRECO_fit-fitEtBins6systEtaBins5_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins6systEtaBins5_max25_PU";
    label1="Et6Eta5";
    label2="Et6Eta5_max25";
    fnameTag="-et6eta5-vars--";
    transLegendX=-0.4;
    allowInvert12=0;
  }

  if (0) { // 2014.05.10
    if (!vsEt) {
      std::cout << "Et bins are different\n";
      return;
    }
    HLTcomparison=0;
    relRatio=-1;
    path1="Results-DYee-effBinStudy/root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy/";
    path2=path1;
    effKindLongStr1="dataRECO_fit-fitEtBins6systEtaBins5_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins7altEtaBins5_PU";
    label1="Et6Eta5";
    label2="Et7altEta5";
    fnameTag="-eta5-vars--";
    //transLegendX=-0.4;
    transLegendX=-0.07;
    allowInvert12=0;
  }

  if (0) { // 2014.05.10
    HLTcomparison=0;
    relRatio=-1;
    path1="Results-DYee-effBinStudy/root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy/";
    path2=path1;
    effKindLongStr1="dataRECO_fit-fitEtBins6systEtaBins7_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins7altEtaBins7_PU";
    label1="Et6Eta7";
    label2="Et7altEta7";
    fnameTag="-eta7-vars--";
    //transLegendX=-0.4;
    transLegendX=-0.07;
    allowInvert12=0;
  }

  if (0) { // 2014.03.19
    HLTcomparison=0;
    allowInvert12=0;
    relRatio=-1;
    path1="../Covariance/dir-toys/";
    path2=path1;
    path3=path1;
    effKindLongStr1="dataRECO_fit-fitEtBins6EtaBins2_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins6EtaBins5_PU";
    effKindLongStr3="dataRECO_fit-fitEtBins6EtaBins9_PU";
    fnameBase="efficiency_TnP_var0-2D_Full2012_";
    label1="Et6Eta2 v.0";
    label2="Et6Eta5 v.0";
    label3="Et6Eta9 v.0";
    fnameTag="-toy-varEt6Et-var0";
    transLegendX=-0.4;
  }

  if (0) { // 2014.03.19
    HLTcomparison=0;
    allowInvert12=0;
    relRatio=-1;
    path1="../Covariance/dir-toys/";
    path2=path1;
    path3=path1;
    effKindLongStr1="dataRECO_fit-fitEtBins6EtaBins2_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins6EtaBins5_PU";
    effKindLongStr3="dataRECO_fit-fitEtBins6EtaBins9_PU";
    fnameBase="efficiency_TnP_var1-2D_Full2012_";
    label1="Et6Eta2 v.1";
    label2="Et6Eta5 v.1";
    label3="Et6Eta9 v.1";
    fnameTag="-toy-varEt6Et-var1";
    transLegendX=-0.4;
  }

  if (0) { // 2014.03.19
    HLTcomparison=0;
    allowInvert12=0;
    relRatio=-1;
    path1="../Covariance/dir-toys/";
    path2=path1;
    path3=path1;
    effKindLongStr1="dataRECO_fit-fitEtBins6EtaBins5_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins6EtaBins5_PU";
    fnameBase="efficiency_TnP_var2-2D_Full2012_";
    label1="Et6Eta5 v.2";
    label2="";
    fnameTag="-toy-varEt6Et-var2";
    transLegendX=-0.4;
  }

  if (0) { // 2014.03.19
    HLTcomparison=0;
    allowInvert12=0;
    relRatio=-1;
    path1="../Covariance/dir-toys/";
    path2=path1;
    path3=path1;
    effKindLongStr1="dataRECO_fit-fitEtBins6EtaBins2_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins6EtaBins5_PU";
    effKindLongStr3="dataRECO_fit-fitEtBins6EtaBins9_PU";
    fnameBase="efficiency_TnP_var3-2D_Full2012_";
    label1="Et6Eta2 v.3";
    label2="Et6Eta5 v.3";
    label3="Et6Eta9 v.3";
    fnameTag="-toy-varEt6Et-var3";
    transLegendX=-0.4;
  }

  if (0) { // 2014.05.23
    HLTcomparison=0;
    allowInvert12=0;
    relRatio=-1;
    path1="Results-DYee-effBinStudy/root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy/";
    path2="dir-Rami-20140511/dir-Effs/";
    effKindLongStr1="dataRECO_fit-fitEtBins6systEtaBins5_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins6EtaBins5_PU";
    label1="Et6Eta5";
    label2="Et6Eta5 (Rami)";
    fnameTag="-eff-cmpRami-vars--";
    transLegendX=-0.4;
    transLegendY=-0.06;
  }

  if (0) { // 2014.05.15
    HLTcomparison=0;
    relRatio=-1;
    path1="Results-DYee-effBinStudy/root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy/";
    path2="Results-DYee-idEffStudy/root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy/";
    path3=path2;
    effKindLongStr1="dataID_fit-fitEtBins6systEtaBins5_PU";
    effKindLongStr2="dataID_fit-fitEtBins6systEtaBins5-idSystdxyz_095_PU";
    effKindLongStr3="dataID_fit-fitEtBins6systEtaBins5-idSystdxyz_105_PU";
    label1="Et6Eta5 base";
    label2="Et6Eta5 dxyz#times0.95";
    label3="Et6Eta5 dxyz#times1.05";
    fnameTag="-idEffStudy-dxyz";
    transLegendX=-0.4;
    transLegendY=-0.06;
    allowInvert12=0;
    setId=1;
  }


  // -------------------------------
  // processing
  // -------------------------------

  if (transLegendY_user!=0.) transLegendY=transLegendY_user;
  if (transLegendX_user!=0.) transLegendX=transLegendX_user;

  if (setId==0) {
  if (label2 == TString("EGamma")) {
    if (iBr==0) {
    }
    else if (iBr==1) {
      effKindLongStr1.ReplaceAll("dataID_fit-fit","mcID_count-count");
      effKindLongStr2.ReplaceAll("dataID","mcID");
    }
    else if (iBr==2) {
      effKindLongStr1.ReplaceAll("dataID_fit-fit","dataRECO_fit-fit");
      effKindLongStr2.ReplaceAll("dataID","dataRECO");
    }
    else if (iBr==3) {
      effKindLongStr1.ReplaceAll("dataRECO_fit-fit","mcRECO_count-count");
      effKindLongStr2.ReplaceAll("dataRECO","mcRECO");
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
      effKindLongStr3.ReplaceAll("RECO","ID");
    }
    else if (iBr==2) {
      effKindLongStr1.ReplaceAll("RECO_fit-fit","HLT_count-count");
      effKindLongStr2.ReplaceAll("RECO_fit-fit","HLT_count-count");
      effKindLongStr3.ReplaceAll("RECO_fit-fit","HLT_count-count");
    }
    else if (iBr==3) {
      effKindLongStr1.ReplaceAll("dataRECO_fit-fit","mcRECO_count-count");
      effKindLongStr2.ReplaceAll("dataRECO_fit-fit","mcRECO_count-count");
      effKindLongStr3.ReplaceAll("dataRECO_fit-fit","mcRECO_count-count");
    }
    else if (iBr==4) {
      effKindLongStr1.ReplaceAll("dataRECO_fit-fit","mcID_count-count");
      effKindLongStr2.ReplaceAll("dataRECO_fit-fit","mcID_count-count");
      effKindLongStr3.ReplaceAll("dataRECO_fit-fit","mcID_count-count");
    }
    else if (iBr==5) {
      effKindLongStr1.ReplaceAll("dataRECO_fit-fit","mcHLT_count-count");
      effKindLongStr2.ReplaceAll("dataRECO_fit-fit","mcHLT_count-count");
      effKindLongStr3.ReplaceAll("dataRECO_fit-fit","mcHLT_count-count");
    }
    else if (iBr==6) {
      effKindLongStr1.ReplaceAll("dataRECO_fit-fit","dataHLTleg1_count-count");
      effKindLongStr2.ReplaceAll("dataRECO_fit-fit","dataHLTleg1_count-count");
      effKindLongStr3.ReplaceAll("dataRECO_fit-fit","dataHLTleg1_count-count");
    }
    else if (iBr==7) {
      effKindLongStr1.ReplaceAll("dataRECO_fit-fit","dataHLTleg2_count-count");
      effKindLongStr2.ReplaceAll("dataRECO_fit-fit","dataHLTleg2_count-count");
      effKindLongStr3.ReplaceAll("dataRECO_fit-fit","dataHLTleg2_count-count");
    }
    else if (iBr==8) {
      effKindLongStr1.ReplaceAll("dataRECO_fit-fit","mcHLTleg1_count-count");
      effKindLongStr2.ReplaceAll("dataRECO_fit-fit","mcHLTleg1_count-count");
      effKindLongStr3.ReplaceAll("dataRECO_fit-fit","mcHLTleg1_count-count");
    }
    else if (iBr==9) {
      effKindLongStr1.ReplaceAll("dataRECO_fit-fit","mcHLTleg2_count-count");
      effKindLongStr2.ReplaceAll("dataRECO_fit-fit","mcHLTleg2_count-count");
      effKindLongStr3.ReplaceAll("dataRECO_fit-fit","mcHLTleg2_count-count");
    }
    else {
      std::cout << "iBr error\n";
      return;
    }
  }
  }
  else if (setId==1) { // idEffStudy
    std::vector<TString*> tsv;
    tsv.reserve(10);
    tsv.push_back(&effKindLongStr1); // to change data->mc
    tsv.push_back(&effKindLongStr2);
    tsv.push_back(&effKindLongStr3);
    tsv.push_back(&label2);
    tsv.push_back(&label3);
    tsv.push_back(&fnameTag);

    TString oldText="dxyz";
    TString newText;

    if (iBr>=10) {
      if (!replaceAll(tsv,"dataID_fit-fit","mcID_count-count")) {
	std::cout << "failed to change 'data' to 'mc'\n";
	return;
      }
    }


    switch(iBr%10) {
    case 0: break;
    case 1: newText="invEminusInvP"; break;
    case 2: newText="Aeff"; break;
    case 3: newText="relPFIso"; break;
    case 4: newText="dEta"; break;
    case 5: newText="dPhi"; break;
    case 6: newText="sigmaIEtaIEta"; break;
    case 7: newText="HoverE"; break;
    default:
      std::cout << "setId=" << setId << " not ready for iBr=" << iBr << "\n";
      return;
    }
    if (newText.Length()) {
      if (!replaceAll(tsv,oldText,newText)) {
	std::cout << "could not find the replacements\n";
	return;
      }
    }
  }
  else {
    std::cout << "setId=" << setId << " is not ready\n";
    return;
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

  TGraphAsymmErrors* gr1=getAsymGraph(vsEt,etBinSet1,etaBinSet1,iBin,*eff1,*eff1ErrLo,*eff1ErrHi,&histo1,histo1Name);
  //std::cout << gr1->GetTitle() << ": "; gr1->Print("range");

  TGraphAsymmErrors* gr2=getAsymGraph(vsEt,etBinSet2,etaBinSet2,iBin,*eff2,*eff2ErrLo,*eff2ErrHi,&histo2,histo2Name);
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
    gr3=getAsymGraph(vsEt,etBinSet3,etaBinSet3,iBin,*eff3,*eff3ErrLo,*eff3ErrHi, &histo3,histo3Name);
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

  double *loc_etBinLimits=DYTools::getEtBinLimits(etBinSet1);
  double *loc_etaBinLimits=DYTools::getEtaBinLimits(etaBinSet1);
  int signedEta=DYTools::signedEtaBinning(etaBinSet1);
  TString cpTitle;
  if (vsEt) cpTitle= dataKind+ TString(Form(" %5.3lf #leq %s #leq %5.3lf",loc_etaBinLimits[iBin],(signedEta)?"#eta":"abs(#eta)",loc_etaBinLimits[iBin+1]));
  else cpTitle= dataKind+ TString(Form(" %2.0lf #leq #it{E}_{T} #leq %2.0lf GeV",loc_etBinLimits[iBin],loc_etBinLimits[iBin+1]));
  TString xaxisTitle=(vsEt) ? "#it{E}_{T}" : ((signedEta) ? "#eta" : "|#eta|");

  ComparisonPlot_t cp(ComparisonPlot_t::_ratioRel,"comp",cpTitle,
			  xaxisTitle,"efficiency","ratio");
  cp.SetRefIdx(-111); // do not plot lower panel
  if (vsEt) cp.SetLogx();

  if (gr3 && HLTcomparison) { // for HLT efficiency
    cp.SetYRange(0.0,1.02);
  }
  cp.SetYRange(0,1.2);

  // Square canvas if ratio is not plotted
  int cxHeight=(relRatio==-1) ? 600 : 700;
  TCanvas *cx=new TCanvas("cx","cx",600,cxHeight);
  if (relRatio!=-1) cp.Prepare2Pads(cx);

  gr1->GetYaxis()->SetNdivisions(506);
  gr1->GetYaxis()->SetTitleOffset(1.4);

  int color1=kBlack;
  int color2=kBlue;
  int color3=kGreen+1;

  if (1) { // for comparison with Rami
    color1=kBlue;
    color2=kRed;
    color3=kGreen;
  }

  if (gr3 && !HLTcomparison && allowInvert12) {
    std::cout << "\n\tInverted plotting order 2,1\n";
    gr2->GetYaxis()->SetTitleOffset(1.2);
    //gr1->SetMarkerStyle(24);
    div->SetMarkerStyle(24);
    cp.AddGraph(gr2,label2,"LPE1",color2);
    cp.AddGraph(gr1,label1," PE1",color1,24);
  }
  else {
    cp.AddGraph(gr1,label1,"PE1",color1);
    if (label2.Length()) cp.AddGraph(gr2,label2,"PE1",color2,24);
  }
  if (gr3) {
    //gr3->SetMarkerStyle(27);
    div31->SetMarkerStyle(27);
    cp.AddGraph(gr3,label3," PE1",color3,27);
  }

  int targetPad=(relRatio==-1) ? 0 : 1;
  cp.Draw(cx,0,"png",targetPad);
  cp.TransLegend(transLegendX, transLegendY);
  cp.WidenLegend(0.2,0.);
  if (setId==1) cp.WidenLegend(0.,-0.03);

  if (relRatio!=-1) {
  cx->cd(2);
  cx->GetPad(2)->SetLogx(cp.fLogx);
  div->SetTitle("");
  div->GetXaxis()->SetTitle(xaxisTitle);
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
	if (iBin==2) div->GetYaxis()->SetRangeUser(0.9,1.1);
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
  }

  cx->Update();

  if (fnameTag.Length()) {
    TString fname=TString("fig-eff-") + fnameTag + cpTitle;
    fname.ReplaceAll(" #leq "," ");
    fname.ReplaceAll(" ","_");
    fname.ReplaceAll("(#eta)","Eta");
    fname.ReplaceAll("#eta","eta");
    fname.ReplaceAll(".","_");
    fname.ReplaceAll("#it{E}_{T}","Et");
    //fname.Append(".png");
    TString locOutDir=TString("plots") + fnameTag;
    locOutDir.ReplaceAll("--","");

    std::cout << "fnameBase=<" << fname << "> in <" << locOutDir << ">\n";
    if (outFileName_ptr) *outFileName_ptr=fname;
    if (outDir_ptr) *outDir_ptr=locOutDir;

    if (doSave) {
      SaveCanvas(cx,fname,locOutDir);
    }
    else {
      std::cout <<  " ... not saved (as requested)\n";
    }
  }

  return ;
}
