#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
//#include "../Include/MitStyleRemix.hh"
#include "../Include/colorPalettes.hh"
#include "../Include/ComparisonPlot.hh"

// ------------------------------------------------------------

TH2D* loadESF(TString fname, TString label) {
  int checkBinning=0;
  TH2D* h2=LoadMatrixFields(fname, checkBinning, "scaleFactor", "scaleFactorErr", 1, 1);
  h2->SetTitle(label);
  h2->GetXaxis()->SetMoreLogLabels();
  h2->GetXaxis()->SetNoExponent();
  std::cout << "loading from <" << fname << ">\n";
  std::cout << "got "; h2->Print("range");
  return h2;
}

// ------------------------------------------------------------

TH2D* loadESF_2D_7TeV(TString fname, TString label) {
  TFile fin(fname,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << fname << ">\n";
    return NULL;
  }
  TMatrixD* sf=(TMatrixD*) fin.Get("scaleFactor");
  TMatrixD* sfErr=(TMatrixD*) fin.Get("scaleFactorErr");
  fin.Close();

  const int nMBins=7;
  const double mbLims2D[nMBins+1]= {0, // first bin is underflow
			  20, 30, 45, 60, 120, 200, 1500 };
  const int yBinCount=24;
  const double yMax=2.5;

  if ((sf->GetNrows()!=nMBins) || (sf->GetNcols()!=yBinCount)) {
    std::cout << "loadESF_2D_7TeV:\n";
    std::cout << "\tfailed the range check\n";
    assert(0);
  }

  TString hname=TString("h2_ScaleFactor_") + label;
  TH2D* h2=new TH2D(hname,label, nMBins, mbLims2D, yBinCount,0.,yMax);
  h2->SetDirectory(0);
  h2->SetStats(0);

  for (int im=0; im<nMBins; ++im) {
    for (int iy=0; iy<yBinCount; ++iy) {
      h2->SetBinContent(im+1,iy+1, (*sf)(im,iy));
      h2->SetBinError  (im+1,iy+1, (*sfErr)(im,iy));
    }
  }
  delete sf;
  delete sfErr;

  h2->SetTitle(label);
  h2->GetXaxis()->SetMoreLogLabels();
  h2->GetXaxis()->SetNoExponent();
  return h2;
}

// ------------------------------------------------------------

TMatrixD* enforceYRanges(int the_case) {
  TMatrixD* r=new TMatrixD(6,2);
  if (the_case==1) {
    (*r)(0,0)=0.85; (*r)(0,1)=1.15;
    (*r)(1,0)=0.85; (*r)(1,1)=1.10;
    (*r)(2,0)=0.89; (*r)(2,1)=1.11;
    (*r)(3,0)=0.94; (*r)(3,1)=1.10;
    (*r)(4,0)=0.95; (*r)(4,1)=1.12;
    (*r)(5,0)=0.95; (*r)(5,1)=1.12;
  }
  return r;
}

// ------------------------------------------------------------


void compareESF(int iBr=0,
		int nDim=1,
		int doSave=0,
		TString *figName=NULL,
		TString *dirName=NULL) {

  TString esfLongStr1="1D_Full2012_hltEffOld_PU";
  TString esfLongStr2="1D";
  TString esfLongStr3,esfLongStr4;
  TString path1="/home/andriusj/cms/DYee8TeV-20140118/root_files/constants/DY_j22_19712pb/";
  TString path2="/home/andriusj/cms/DYee-20131024/root_files_reg/constants/DY_j22_19712pb/";
  TString path3,path4;
  //path2="/home/andriusj/cms/CMSSW_3_8_4/src/DYee8TeV-20130801/DrellYanDMDY/root_files/tag_and_probe/DY_j22_19789pb/";
  //path2="/home/andriusj/cms/DYee8TeV-20140118-maxEta24/root_files/tag_and_probe/DY_j22_19712pb/";

  if (((nDim==1) &&  DYTools::study2D) ||
      ((nDim==2) && !DYTools::study2D)) {
    std::cout << "the macro uses basic implementations that require nDim info to match study2D in DYTools.hh\n";
    return;
  }


  TString fnameBase1="scale_factors_";
  TString fnameBase2="scale_factors_";
  TString fnameBase3, fnameBase4;

  TString dimStr;
  if (nDim==1) dimStr.Append("1D");
  else if (nDim==2) dimStr.Append("2D");
  else {
    std::cout << "error: nDim=" << nDim << "\n";
    return;
  }


  TString saveDirTag,fnameTag;
  TString label1="DrellYanDMDY";
  TString label2="DYee";
  TString label3,label4;
  int secondIs7TeV=0;
  double set_ratio_y_min=0.96;
  double set_ratio_y_max=1.04;
  double transLegendX=(nDim==1) ? -0.2 : -0.42;

  int compSet=-1;
  int swapColors=0;

  TMatrixD *setYRanges2D=NULL;

  if (0) {
    path1="/home/andriusj/cms/CMSSW_3_8_4/src/DYee8TeV-20130801/DrellYanDMDY/root_files/constants/DY_j22_19789pb/";
    path2="/home/andriusj/cms/DMDY-Ilya-20130808-my-copy/root_files/constants/DY_j22_19789pb/";
    esfLongStr1="1D_Full2012_hltEffOld_PU";
    esfLongStr2="1D_Full2012_hltEffOld_PU";
    label1="DYDM |#eta|<2.5 (Summer2013,my)";
    label2="DYDM |#eta|<2.5 (Ilya)";
    //secondIs7TeV=1;
    fnameTag="-DMDY";
    saveDirTag="-checkSummer2013-esf--";
  }

  if (0) {
    path1="/home/andriusj/cms/CMSSW_3_8_4/src/DYee8TeV-20130801/DrellYanDMDY/root_files/constants/DY_j22_19789pb/";
    path1="/home/andriusj/cms/DMDY-Ilya-20130808-my-copy/root_files/constants/DY_j22_19789pb/";
    path2="/home/andriusj/cms/DYee8TeV-20140118/root_files/constants/DY_j22_19712pb/"; // new n-tuples
    esfLongStr1="1D_Full2012_hltEffOld_PU";
    esfLongStr2="1D_Full2012_hltEffOld_PU";
    label1="DYDM |#eta|<2.5 (Summer 2013)";
    label2="DYDM |#eta|<2.5 (new n-tuples)";
    if (nDim==2) {  set_ratio_y_min=0.96; set_ratio_y_max=1.06; }
    else { set_ratio_y_min=0.98; set_ratio_y_max=1.08; }
    fnameTag="-DMDY";
    saveDirTag="-new_vs_old-esf--";
  }

  if (0) {
    path1="/home/andriusj/cms/DYee-20131024/root_files_reg/constants/DY_j22_19712pb/"; 
    path2="/home/andriusj/cms/DYee8TeV-20140118/root_files/constants/DY_j22_19712pb/";
    esfLongStr1="1D";
    esfLongStr2="1D_Full2012_hltEffOld_PU";
    label1="DYee |#eta|<2.5";
    label2="DMDY |#eta|<2.5";
    set_ratio_y_min=0.96; 
    set_ratio_y_max=1.04;
    fnameTag="-cmpPkg-Eta25";
    saveDirTag="-cmpPkg-Eta25-esf--";
  }

  if (0) {
    path1="/home/andriusj/cms/DYee-20131024/root_files_reg/constants/DY_j22_19712pb/"; 
    path2="/home/andriusj/cms/DYee8TeV-20140118-maxEta24/root_files/constants/DY_j22_19712pb/";
    esfLongStr1="etaMax24_1D";
    esfLongStr2="1D_Full2012_hltEffOld_PU";
    label1="DYee |#eta|<2.4";
    label2="DMDY |#eta|<2.4";
    set_ratio_y_min=0.99; 
    set_ratio_y_max=1.01;
    fnameTag="-cmpPkg-Eta24";
    saveDirTag="-cmpPkg-esf--";
  }

  if (0) {
    path1="/home/andriusj/cms/DYee-20131024/root_files_reg/constants/DY_j22_19712pb/"; 
    path2=path1;
    esfLongStr1="1D";
    esfLongStr2="etaMax24_1D";
    label1="DYee |#eta|<2.5";
    label2="DYee |#eta|<2.4";
    if (nDim==2) { set_ratio_y_min=0.95; set_ratio_y_max=1.05; }
    else { set_ratio_y_min=0.99; set_ratio_y_max=1.01; }
    fnameTag="-DYee";
    saveDirTag="-diffEtaMax-esf--";
  }

  if (0) {
    path1="/home/andriusj/cms/DYee-20131024/root_files_reg/constants/DY_j22_19712pb/"; 
    path2=path1;
    esfLongStr1="etaMax24_1D";
    esfLongStr2="etaMax24_asymHLT_1D";
    label1="DYee |#eta|<2.4";
    label2="DYee |#eta|<2.4 asym.HLT";
    set_ratio_y_min=0.99; 
    set_ratio_y_max=1.01;
    fnameTag="-DYee-asymHLT";
    saveDirTag="-diffEtaMax-esf--";
  }

  if (0) {
    path1="/home/andriusj/cms/CMSSW_3_8_4/src/DYee8TeV-20130801/DrellYanDMDY/root_files/constants/DY_j22_19789pb/";
    //path1="/home/andriusj/cms/DYee8TeV-20140118/root_files/constants/DY_j22_19712pb/"; // new n-tuples
    path2="/home/andriusj/cms/DYee-20131024/root_files_reg/constants/DY_j22_19712pb/"; 
    esfLongStr1="1D_Full2012_hltEffOld_PU";
    esfLongStr2="etaMax24_asymHLT_1D";
    label1="Summer 2013";
    label2="DYee |#eta|<2.4 asym.HLT";
    if (nDim==2) { set_ratio_y_min=0.95; set_ratio_y_max=1.05; }
    else { set_ratio_y_min=0.98; set_ratio_y_max=1.05; }
    fnameTag="-totChange";
    saveDirTag="-diffEtaMax-esf--";
  }

  if (0) { // added 2014.02.19
    path1="/home/andriusj/cms/DYee-20140217-EESF/root_files_reg/constants/DY_j22_19712pb/";
    path2="/home/andriusj/cms/DYee-20140217-EESF/root_files_reg/constants/DY_j22_19712pb_egamma/"; 
    esfLongStr1="ourSF_asymHLT_1D";
    esfLongStr2="egamma_asymHLT_1D";
    label1="our SF";
    label2="EGamma";
    if (nDim==2) { set_ratio_y_min=0.95; set_ratio_y_max=1.05; }
    else { set_ratio_y_min=0.98; set_ratio_y_max=1.05; }
    fnameTag="-egamma";
    saveDirTag="-egamma--";
  }

  if (0) { // added 2014.02.26
    path1="../Covariance/";
    path2="../Covariance/";
    fnameBase1="covRhoFileSF"; fnameBase2=fnameBase1;
    esfLongStr1="_nMB41_egamma_asymHLT_Unregressed_energy-allSyst_100";
    esfLongStr2="_nMB41_egamma_asymHLT_Unregressed_energy-noSyst_100";
    label1="EGamma (unreg.en.+all syst.)";
    label2="EGamma (unreg.en.)";
    fnameTag="-egamma-unregEn-allSyst";
    saveDirTag="-egamma-unregEn--";
    if (nDim==2) { set_ratio_y_min=0.95; set_ratio_y_max=1.05; }
    else { set_ratio_y_min=0.94; set_ratio_y_max=1.06; }
    compSet=10;
    swapColors=0;
  }

  if (0) { // added 2014.03.01
    path1="../Covariance/";
    path2="../Covariance/";
    fnameBase1="covRhoFileSF"; fnameBase2=fnameBase1;
    esfLongStr1="_nMB41_egamma_asymHLT_Unregressed_energy-allSyst_100";
    esfLongStr2="_nMB41_egamma_asymHLT_Unregressed_energy-allSyst-egammaSystOnly_100";
    label1="EGamma (unreg.en.+all syst.)";
    label2="EGamma (unreg.en.+EG(reco,id) +HLT)";
    fnameTag="-egamma-unregEn-allSyst-egammaSystOnly";
    saveDirTag="-egamma-unregEn--";
    if (nDim==2) { set_ratio_y_min=0.95; set_ratio_y_max=1.05; }
    else { set_ratio_y_min=0.93; set_ratio_y_max=1.07; }
    compSet=13;
    swapColors=(nDim==1) ? 2 : 0;
    setYRanges2D=enforceYRanges(1);
  }

  if (0) { // added 2014.03.01
    path1="../Covariance/";
    path2="../Covariance/";
    fnameBase1="covRhoFileSF"; fnameBase2=fnameBase1;
    esfLongStr1="_nMB41_egamma_asymHLT_Unregressed_energy-allSyst_100";
    esfLongStr2="_nMB41_egamma_asymHLT_Unregressed_energy-allSyst_1000";
    label1="EGamma (unreg.en.+all syst. 100)";
    label2="EGamma (unreg.en.+all.syst. 1000)";
    fnameTag="-egamma-unregEn-allSyst-1000";
    saveDirTag="-egamma-unregEn--";
    if (nDim==2) { set_ratio_y_min=0.95; set_ratio_y_max=1.05; }
    else { set_ratio_y_min=0.93; set_ratio_y_max=1.07; }
    compSet=13;
    swapColors=0;
    setYRanges2D=enforceYRanges(1);
  }

  if (0) { // added 2014.03.03
    path1="../Covariance/";
    path2="../Covariance/";
    fnameBase1="covRhoFileSF"; fnameBase2=fnameBase1;
    esfLongStr1="_nMB41_egamma_asymHLT_Unregressed_energy-recoSyst_100";
    esfLongStr2="_nMB41_egamma_asymHLT_Unregressed_energy-recoSyst-egammaSystOnly_100";
    label1="EGamma (unreg.en.+reco syst.)";
    label2="EGamma (unreg.en.+EG(reco))";
    fnameTag="-egamma-unregEn-recoSyst-egammaSystOnly";
    saveDirTag="-egamma-unregEn--";
    if (nDim==2) { set_ratio_y_min=0.95; set_ratio_y_max=1.05; }
    else { set_ratio_y_min=0.93; set_ratio_y_max=1.07; }
    compSet=12;
    swapColors=0;
    setYRanges2D=enforceYRanges(1);
  }

  if (0) { // added 2014.02.27
    path1="../Covariance/";
    path2="../Covariance/";
    fnameBase1="covRhoFileSF"; fnameBase2=fnameBase1;
    esfLongStr1="_nMB41_egamma_asymHLT_regEn-allSyst_100";
    esfLongStr2="_nMB41_egamma_asymHLT_regEn_100";
    label1="EGamma (reg.en.+all syst.)";
    label2="EGamma (reg.en.)";
    fnameTag="-egamma-regEn-allSyst";
    saveDirTag="-egamma-regEn--";
    if (nDim==2) { set_ratio_y_min=0.95; set_ratio_y_max=1.05; }
    else { set_ratio_y_min=0.94; set_ratio_y_max=1.06; }
    compSet=10;
  }

  if (0) { // added 2014.02.28
    path1="../Covariance/";
    path2="../Covariance/";
    fnameBase1="covRhoFileSF"; fnameBase2=fnameBase1;
    esfLongStr1="_nMB41_egamma_asymHLT_regEn-allSyst_100";
    esfLongStr2="_nMB41_egamma_asymHLT_Unregressed_energy-allSyst_100";
    label1="EGamma (reg.en.+all syst.)";
    label2="EGamma (unreg.en.+all syst.)";
    fnameTag="-egamma-cmpEn-allSyst";
    saveDirTag="-egamma-cmpEn--";
    if (nDim==2) { set_ratio_y_min=0.95; set_ratio_y_max=1.05; }
    else { set_ratio_y_min=0.94; set_ratio_y_max=1.06; }
    compSet=11;
    swapColors=1;
  }

  if (0) { // added 2014.03.08
    path1="../Covariance/";
    path2="../Covariance/";
    path3="../Covariance/";
    path4="../Covariance/";
    fnameBase1="covRhoFileSF";
    fnameBase2="covRhoFileSF_RECO"; 
    fnameBase3="covRhoFileSF_ID";
    fnameBase4="covRhoFileSF_dtHLT";
    esfLongStr1="_nMB41_egamma_asymHLT_Unregressed_energy-allSyst_100";
    esfLongStr2="_nMB41_egamma_asymHLT_Unregressed_energy-recoSyst_100";
    esfLongStr3="_nMB41_egamma_asymHLT_Unregressed_energy-idSyst_100";
    esfLongStr4="_nMB41_egamma_asymHLT_Unregressed_energy-hltSyst_100";
    label1="full SF (unreg.en.+all syst.)";
    label2="RECO SF (unreg.en.+RECO syst.)";
    label3="ID SF   (unreg.en.+ID syst)";
    label4="HLT SF  (unreg.en.+HLT syst)";
    fnameTag="-egamma-unregEn-diffSF";
    saveDirTag="-egamma-unregEn--";
    if (nDim==2) { set_ratio_y_min=0.9; set_ratio_y_max=1.15; }
    else { set_ratio_y_min=0.95; set_ratio_y_max=1.15; }
    compSet=-10;
    swapColors=0;
    //setYRanges2D=enforceYRanges(1);
  }

  if (1) { // added 2014.03.19
    path1="../Covariance/";
    path2="../Covariance/";
    path3="../Covariance/";
    path4="../Covariance/";
    fnameBase1="covRhoFileSF_RECO";
    fnameBase2="covRhoFileSF_RECO"; 
    fnameBase3="covRhoFileSF_RECO";
    fnameBase4="covRhoFileSF_RECO";
    esfLongStr1="_nMB41_etaMax24_asymHLT_Unregressed_energy_100-EtaBins5_var0";
    esfLongStr2="_nMB41_etaMax24_asymHLT_Unregressed_energy_100-EtaBins5_var1";
    esfLongStr3="_nMB41_etaMax24_asymHLT_Unregressed_energy_100-EtaBins5_var2";
    esfLongStr4="_nMB41_etaMax24_asymHLT_Unregressed_energy_100-EtaBins5_var3";
    label1="EtaBins5 var.0";
    label2="EtaBins5 var.1";
    label3="EtaBins5 var.2";
    label4="EtaBins5 var.3";
    fnameTag="-etaBins5-vars";
    saveDirTag="-toyStudy--";
    //if (nDim==2) { set_ratio_y_min=0.9; set_ratio_y_max=1.15; }
    //else { set_ratio_y_min=0.95; set_ratio_y_max=1.15; }
    set_ratio_y_max=set_ratio_y_min; // force auto-setup of range
    transLegendX=-0.2;
    compSet=15;
    swapColors=0;
    //setYRanges2D=enforceYRanges(1);
  }

  if (iBr!=0) {
    if (compSet==10) {
      if (iBr==1) {
	esfLongStr1.ReplaceAll("allSyst","recoSyst");
	label1.ReplaceAll("all syst","reco syst");
	fnameTag.ReplaceAll("allSyst","recoSyst");
      }
      else if (iBr==2) {
	esfLongStr1.ReplaceAll("allSyst","idSyst");
	label1.ReplaceAll("all syst","ID syst");
	fnameTag.ReplaceAll("allSyst","idSyst");
      }
      else if (iBr==3) {
	esfLongStr1.ReplaceAll("allSyst","hltSyst");
	label1.ReplaceAll("all syst","HLT syst");
	fnameTag.ReplaceAll("allSyst","hltSyst");
      }
      else if (iBr==-1) {
	esfLongStr2.ReplaceAll("-noSyst_100","-recoSyst_100");
	label2.ReplaceAll("en.","en.+reco syst");
	fnameTag.ReplaceAll("allSyst","allSyst_vsRecoSyst");
      }
      else if (iBr==-2) {
	esfLongStr2.ReplaceAll("-noSyst_100","-idSyst_100");
	label2.ReplaceAll("en.","en.+ID syst");
	fnameTag.ReplaceAll("allSyst","allSyst_vsIDSyst");
      }
      else if (iBr==-3) {
	esfLongStr2.ReplaceAll("-noSyst_100","-hltSyst_100");
	label2.ReplaceAll("en.","en.+HLT syst");
	fnameTag.ReplaceAll("allSyst","allSyst_vsHLTSyst");
      }
    }
    else if (compSet==12) {
      if (iBr==0) {
      }
      else if (iBr==1) {
	esfLongStr1.ReplaceAll("reco","id");
	esfLongStr2.ReplaceAll("reco","id");
	label1.ReplaceAll("reco","id");
	label2.ReplaceAll("reco","id");
	fnameTag.ReplaceAll("reco","id");
      }
    }
    else if (compSet==15) {
      if (iBr==0) {
	// no change
      }
      else {
	fnameBase4.Clear();
	if (iBr==1) {
	  transLegendX=-0.4;
	  esfLongStr1="_nMB41_etaMax24_asymHLT_Unregressed_energy_100-EtaBins2_var0";
	  esfLongStr2="_nMB41_etaMax24_asymHLT_Unregressed_energy_100-EtaBins5_var0";
	  esfLongStr3="_nMB41_etaMax24_asymHLT_Unregressed_energy_100-EtaBins9_var0";
	  label1="EtaBins2 var.0";
	  label2="EtaBins5 var.0";
	  label3="EtaBins9 var.0";
	  fnameTag="-varEtaBins-var0";
	}
	else if (iBr==2) {
	  transLegendX=-0.4;
	  esfLongStr1="_nMB41_etaMax24_asymHLT_Unregressed_energy_100-EtaBins2_var1";
	  esfLongStr2="_nMB41_etaMax24_asymHLT_Unregressed_energy_100-EtaBins5_var1";
	  esfLongStr3="_nMB41_etaMax24_asymHLT_Unregressed_energy_100-EtaBins9_var1";
	  label1="EtaBins2 var.1";
	  label2="EtaBins5 var.1";
	  label3="EtaBins9 var.1";
	  fnameTag="-varEtaBins-var1";
	}
	else if (iBr==3) {
	  transLegendX=-0.2;
	  esfLongStr1="_nMB41_etaMax24_asymHLT_Unregressed_energy_100-EtaBins2_var3";
	  esfLongStr2="_nMB41_etaMax24_asymHLT_Unregressed_energy_100-EtaBins5_var3";
	  esfLongStr3="_nMB41_etaMax24_asymHLT_Unregressed_energy_100-EtaBins9_var3";
	  label1="EtaBins2 var.3";
	  label2="EtaBins5 var.3";
	  label3="EtaBins9 var.3";
	  fnameTag="-varEtaBins-var3";
	}
      }
    }
  }


  if (DYTools::study2D) {
    esfLongStr1.ReplaceAll("1D","2D");
    esfLongStr2.ReplaceAll("1D","2D");
    esfLongStr3.ReplaceAll("1D","2D");
    esfLongStr4.ReplaceAll("1D","2D");
    esfLongStr1.ReplaceAll("nMB41","nMB7");
    esfLongStr2.ReplaceAll("nMB41","nMB7");
    esfLongStr3.ReplaceAll("nMB41","nMB7");
    esfLongStr4.ReplaceAll("nMB41","nMB7");
  }


  TString fname1=path1 + fnameBase1 + esfLongStr1 + TString(".root");
  TString fname2=path2 + fnameBase2 + esfLongStr2 + TString(".root");

  TH2D *h2esf1= loadESF(fname1,label1);
  h2esf1->Print("range");
  TH2D *h2esf2= (!secondIs7TeV) ? loadESF(fname2,label2) : loadESF_2D_7TeV(fname2,label2);

  TH2D *h2esf3=NULL, *h2esf4=NULL;

  if (fnameBase3.Length()) {
    TString fname3= path3 + fnameBase3 + esfLongStr3 + TString(".root");
    h2esf3= loadESF(fname3,label3);
  }
  if (fnameBase4.Length()) {
    TString fname4= path4 + fnameBase4 + esfLongStr4 + TString(".root");
    h2esf4= loadESF(fname4,label4);
  }

  TH2D *h2diff= (TH2D*)h2esf1->Clone("diff");
  if (!secondIs7TeV) {
    h2diff->Add(h2esf2,-1);
  }
  h2diff->SetTitle("diff");

  gStyle->SetPalette(1);
  TCanvas *cx=NULL;

  if (0) {
    cx=new TCanvas("cx","cx",1200,400);

    cx->Divide(3,1);
    AdjustFor2DplotWithHeight(cx);
    
    cx->GetPad(1)->SetLogx();
    cx->GetPad(2)->SetLogx();
    cx->GetPad(3)->SetLogx();

    cx->cd(1);
    h2esf1->Draw("COLZ");
    cx->cd(2);
    h2esf2->Draw("COLZ");
    cx->cd(3);
    h2diff->Draw("COLZ");
    cx->Update();
  }

  TCanvas *cy=NULL;

  if (DYTools::study2D) {
    h2esf1->Print("range");
    h2esf2->Print("range");
    if (h2esf3) h2esf3->Print("range");
    if (h2esf4) h2esf4->Print("range");

    std::vector<TH1D*> hProfEsf1,hProfEsf2;
    std::vector<TH1D*> hProfEsf3,hProfEsf4;

    for (int mbin=2; mbin<=7; mbin++) {
      for (int iset=0; iset<4; ++iset) {
	TH2D *AllEsf=NULL;
	std::vector<TH1D*> *profV;
	switch(iset) {
	case 0: AllEsf=h2esf1; profV=&hProfEsf1; break;
	case 1: AllEsf=h2esf2; profV=&hProfEsf2; break;
	case 2: AllEsf=h2esf3; profV=&hProfEsf3; break;
	case 3: AllEsf=h2esf4; profV=&hProfEsf4; break;
	}
	if (AllEsf==NULL) continue;

	TString hname1=Form("hProf%d_m%d",iset+1,mbin);
	int set_nYBins=DYTools::nYBins[mbin-1]; //(mbin==7) ? 12:24;
	TH1D *h1= createProfileY(AllEsf,mbin,hname1,1,hname1, set_nYBins,0.,DYTools::yRangeMax+1e-4);

	h1->GetYaxis()->SetTitleSize(0.07);
	h1->GetYaxis()->SetTitleOffset(1.08);
	h1->GetYaxis()->SetLabelSize(0.06);
	//if (iset>1) h1->SetMarkerStyle(24);
	if (h2esf4) h1->SetMarkerSize(0.7);
	
	profV->push_back(h1);
      }
    }

    std::vector<ComparisonPlot_t*> cpV;
    for (unsigned int i=0; i<6; ++i) {
      TString mRange=Form("M_%2.0lf_%2.0lf",DYTools::massBinLimits[i+1],DYTools::massBinLimits[i+2]);
      ComparisonPlot_t *cp=NULL;
      //cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioRel,Form("cp_%s",mRange.Data()),mRange,"|y|","event scale factor","rel.ratio");
      cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,Form("cp_%s",mRange.Data()),mRange,"|y|","event scale factor","ratio");
      cp->SetRatioNdivisions(404);
      cp->SetYTitleSize(0.08,0.97);
      cp->SetYLabelSize(0.06);
      cp->SetRatioYTitleSize(0.18,0.40);

      if (setYRanges2D) {
	cp->SetYRange((*setYRanges2D)(i,0),(*setYRanges2D)(i,1));
      }
      //cp->SetRatioYRange(0.97,1.01);
      cp->SetRatioYRange(set_ratio_y_min,set_ratio_y_max);
      //cp->SetRatioYRange(-0.01,0.01);

      int showLegend=1;
      if ((compSet==-10) && (i!=4)) showLegend=-1;

      cp->AddHist1D(hProfEsf1[i],label1,"LPE1",kBlack,1,0,showLegend);
      cp->AddHist1D(hProfEsf2[i],label2,"PE3",kBlue ,2,0,showLegend);
      if (i<hProfEsf3.size()) cp->AddHist1D(hProfEsf3[i],label3,"PE3",kRed+1, 3,0,showLegend);
      if (i<hProfEsf4.size()) cp->AddHist1D(hProfEsf4[i],label4,"PE3",kGreen+2, 3,0,showLegend);
      //cp->AddHist1D(hProfEsf2[i],label2,"LPE3",kBlue ,2,0);
      cp->AddLine(0.,1.0, 2.4,1.0, kBlue+1,kDashed);
      cpV.push_back(cp);
    }

    cy=new TCanvas("cy","cy",1100,800);
    cpV[0]->Prepare6Pads(cy, 1);

    int printRatio=0;
    
    for (int i=0; i<6; ++i) {
      if (printRatio) {
	std::cout << "---- " << cpV[i]->GetTitle() << "\n";
	cpV[i]->SetPrintRatio(1);
      }
      cpV[i]->Draw6(cy,1,i+1);
      cpV[i]->TransLegend(transLegendX,-0.1);
      cpV[i]->WidenLegend(0.1,0.1);
      cpV[i]->WidenLegend(0.2,0.0);
      //cpV[i]->
    }
    cy->Update();
  }
  else {
    TH1D *h1= createProfileX(h2esf1,1,"hProf1",1,"hProf1");
    TH1D *h2= createProfileX(h2esf2,1,"hProf2",1,"hProf2");
    TH1D *h3= (h2esf3) ? createProfileX(h2esf3,1,"hProf3",1,"hProf3") : NULL;
    TH1D *h4= (h2esf4) ? createProfileX(h2esf4,1,"hProf4",1,"hProf4") : NULL;

    //std::cout << "\n\ncreated profile from " << h2esf2->GetName() << "\n";
    //h2esf2->Print("range");
    //h2->Print("range");

    int relRatio=0;
    ComparisonPlot_t *cp=NULL;
    if (!relRatio) cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,"cpESF","","M_{ee} [GeV]", "event scale factor","ratio");
    else cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioRel,"cpESF","","M_{ee} [GeV]", "event scale factor","rel.ratio");
    cp->SetYTitleSize(0.07,1.16);
    cp->SetRatioYTitleSize(0.17,0.47);
    cp->SetPrintRatio(1);

    int color1=kBlack;
    int color2=46; // red-brown
    if (h3) color2=kBlue;
    if (swapColors) { int icol=color1; color1=color2; color2=icol; }
    if (swapColors==2) color1=kRed;

    //cp->SetRatioYRangeC(1,0.04);
    cp->SetRatioYRange(set_ratio_y_min, set_ratio_y_max);
    if (relRatio) cp->SetRatioYRangeC(0.,0.01);
    cp->SetLogx();
    //h1->GetYaxis()->SetTitleOffset(1.47);
    cp->AddHist1D(h1, label1,"LPE1",color1,1,0);
    TString opt2=(h3) ? "LPE3" : "LPE1";
    cp->AddHist1D(h2, label2,opt2,color2,1,0);
    if (h3) cp->AddHist1D(h3, label3, "LPE3",kRed+1,1,0);
    if (h4) cp->AddHist1D(h4, label4, "LPE3",kGreen+2,1,0);

    cp->AddLine(DYTools::massBinLimits[0],1.,DYTools::massBinLimits[DYTools::nMassBins],1., kBlack,kDashed);

    cy=new TCanvas("cy","cy",600,700);
    cp->Prepare2Pads(cy);
    cp->Draw(cy);
    cp->TransLegend(transLegendX,-0.6);
    cp->WidenLegend(0.2,0.);
    cy->Update();
  }

  if (cy && fnameTag.Length() && saveDirTag.Length()) {
    TString fname=TString("fig") + saveDirTag;
    fname.Append(Form("%dD",(DYTools::study2D) ? 2:1));
    fname.Append(fnameTag);
    TString locOutDir=TString("plots") + saveDirTag;
    locOutDir.ReplaceAll("--","");
    std::cout << "save <" << fname << "> in <" << locOutDir << ">\n";
    if (figName) *figName=fname;
    if (dirName) *dirName=locOutDir;
    if (doSave) {
      SaveCanvas(cy,fname,locOutDir);
    }
    else {
      std::cout << "not saved due to a request\n";
    }
  }

  //h2diff->Print("range");
  //h2esf1->Print("range");

}
