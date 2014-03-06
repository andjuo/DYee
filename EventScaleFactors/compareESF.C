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
		int doSave=0) {

  TString esfLongStr1="1D_Full2012_hltEffOld_PU";
  TString esfLongStr2="1D";
  TString path1="/home/andriusj/cms/DYee8TeV-20140118/root_files/constants/DY_j22_19712pb/";
  TString path2="/home/andriusj/cms/DYee-20131024/root_files_reg/constants/DY_j22_19712pb/";
  //path2="/home/andriusj/cms/CMSSW_3_8_4/src/DYee8TeV-20130801/DrellYanDMDY/root_files/tag_and_probe/DY_j22_19789pb/";
  //path2="/home/andriusj/cms/DYee8TeV-20140118-maxEta24/root_files/tag_and_probe/DY_j22_19712pb/";

  if (((nDim==1) &&  DYTools::study2D) ||
      ((nDim==2) && !DYTools::study2D)) {
    std::cout << "the macro uses basic implementations that require nDim info to match study2D in DYTools.hh\n";
    return;
  }


  TString fnameBase="scale_factors_";

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
  int secondIs7TeV=0;
  double set_ratio_y_min=0.96;
  double set_ratio_y_max=1.04;

  int compSet=-1;
  int swapColors=0;

  TMatrixD *setYRanges2D=NULL;

  /*
  if (0) {
    path1="/home/andriusj/cms/DYee8TeV-20140118/root_files/constants/DY_j22_19712pb/";
    path2="/home/andriusj/cms/DYee8TeV-20140118-maxEta24/root_files/constants/DY_j22_19712pb/";
    esfLongStr1="1D_Full2012_hltEffOld_PU";
    esfLongStr2="1D_Full2012_hltEffOld_PU";
    label1="DMDY |#eta|<2.5";
    label2="DMDY |#eta|<2.4";
    fnameTag="-DMDY";
    saveDirTag="-diffEtaMax-esf--";
  }

  if (0) {
    path1="/home/andriusj/cms/DYee-20131024/root_files_reg/constants/DY_j22_19712pb/"; 
    path2=path1;
    esfLongStr1="1D";
    esfLongStr2="etaMax24_asymHLT_1D";
    label1="DYee |#eta|<#color[3]{2.5}";
    label2="DYee |#eta|<2.4 asym.HLT";
    //fnameTag="-cmpPkg--";
  }
  */


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
    fnameBase="covRhoFileSF";
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
    fnameBase="covRhoFileSF";
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
    fnameBase="covRhoFileSF";
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

  if (1) { // added 2014.03.03
    path1="../Covariance/";
    path2="../Covariance/";
    fnameBase="covRhoFileSF";
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
    fnameBase="covRhoFileSF";
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
    fnameBase="covRhoFileSF";
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


  if (DYTools::study2D) {
    esfLongStr1.ReplaceAll("1D","2D");
    esfLongStr2.ReplaceAll("1D","2D");
    esfLongStr1.ReplaceAll("nMB41","nMB7");
    esfLongStr2.ReplaceAll("nMB41","nMB7");
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
  }

  TString fname1=path1 + fnameBase + esfLongStr1 + TString(".root");
  TString fname2=path2 + fnameBase + esfLongStr2 + TString(".root");

  TH2D *h2esf1= loadESF(fname1,label1);
  h2esf1->Print("range");
  TH2D *h2esf2= (!secondIs7TeV) ? loadESF(fname2,label2) : loadESF_2D_7TeV(fname2,label2);

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

    std::vector<TH1D*> hProfEsf1,hProfEsf2;
    for (int mbin=2; mbin<=7; mbin++) {
      TString hname1=Form("hProf1_m%d",mbin);
      TString hname2=Form("hProf2_m%d",mbin);
      int set_nYBins=DYTools::nYBins[mbin-1]; //(mbin==7) ? 12:24;
      TH1D *h1= createProfileY(h2esf1,mbin,hname1,1,hname1, set_nYBins,0.,DYTools::yRangeMax+1e-4);
      TH1D *h2= createProfileY(h2esf2,mbin,hname2,1,hname2, set_nYBins,0.,DYTools::yRangeMax+1e-4);

      h1->GetYaxis()->SetTitleSize(0.07);
      h1->GetYaxis()->SetTitleOffset(1.08);
      h1->GetYaxis()->SetLabelSize(0.06);

      hProfEsf1.push_back(h1);
      hProfEsf2.push_back(h2);
    }

    std::vector<ComparisonPlot_t*> cpV;
    for (int i=0; i<6; ++i) {
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
      cp->AddHist1D(hProfEsf1[i],label1,"LPE1",kBlack,1,0);
      cp->AddHist1D(hProfEsf2[i],label2,"PE3",kBlue ,2,0);
      //cp->AddHist1D(hProfEsf2[i],label2,"LPE3",kBlue ,2,0);
      cp->AddLine(0.,1.0, 2.4,1.0, kBlue+1,kDashed);
      cpV.push_back(cp);
    }

    cy=new TCanvas("cy","cy",1200,800);
    cpV[0]->Prepare6Pads(cy, 1);

    int printRatio=0;
    
    for (int i=0; i<6; ++i) {
      if (printRatio) {
	std::cout << "---- " << cpV[i]->GetTitle() << "\n";
	cpV[i]->SetPrintRatio(1);
      }
      cpV[i]->Draw6(cy,1,i+1);
      cpV[i]->TransLegend(-0.42,-0.1);
      cpV[i]->WidenLegend(0.1,0.1);
      cpV[i]->WidenLegend(0.2,0.0);
      //cpV[i]->
    }
    cy->Update();
  }
  else {
    TH1D *h1= createProfileX(h2esf1,1,"hProf1",1,"hProf1");
    TH1D *h2= createProfileX(h2esf2,1,"hProf2",1,"hProf2");

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
    if (swapColors) { int icol=color1; color1=color2; color2=icol; }
    if (swapColors==2) color1=kRed;

    //cp->SetRatioYRangeC(1,0.04);
    cp->SetRatioYRange(set_ratio_y_min, set_ratio_y_max);
    if (relRatio) cp->SetRatioYRangeC(0.,0.01);
    cp->SetLogx();
    //h1->GetYaxis()->SetTitleOffset(1.47);
    cp->AddHist1D(h1, label1,"LPE1",color1,1,0);
    cp->AddHist1D(h2, label2,"LPE1",color2,1,0);

    cp->AddLine(DYTools::massBinLimits[0],1.,DYTools::massBinLimits[DYTools::nMassBins],1., kBlack,kDashed);

    cy=new TCanvas("cy","cy",600,700);
    cp->Prepare2Pads(cy);
    cp->Draw(cy);
    cp->TransLegend(-0.2,-0.6);
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
