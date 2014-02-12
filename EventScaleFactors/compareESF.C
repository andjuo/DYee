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


void compareESF(TString esfLongStr1="1D_Full2012_hltEffOld_PU",
		TString esfLongStr2="1D",
		int nDim=1,
		int doSave=1) {
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

  if (1) {
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



  if (DYTools::study2D) {
    esfLongStr1.ReplaceAll("1D","2D");
    esfLongStr2.ReplaceAll("1D","2D");
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
      ComparisonPlot_t *cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,Form("cp_%s",mRange.Data()),mRange,"|y|","event scale factor","ratio");
      cp->SetRatioNdivisions(404);
      cp->SetYTitleSize(0.08,0.97);
      cp->SetYLabelSize(0.06);
      cp->SetRatioYTitleSize(0.18,0.40);

      //cp->SetRatioYRange(0.97,1.01);
      cp->SetRatioYRange(set_ratio_y_min,set_ratio_y_max);
      cp->AddHist1D(hProfEsf1[i],label1,"LPE1",kBlack,1,0);
      cp->AddHist1D(hProfEsf2[i],label2,"PE3",kBlue ,2,0);
      //cp->AddHist1D(hProfEsf2[i],label2,"LPE3",kBlue ,2,0);
      cpV.push_back(cp);
    }

    cy=new TCanvas("cy","cy",1200,800);
    cpV[0]->Prepare6Pads(cy, 1);
    
    for (int i=0; i<6; ++i) {
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

    ComparisonPlot_t *cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,"cpESF","","M_{ee} [GeV]", "event scale factor","ratio");
    cp->SetYTitleSize(0.07,1.16);
    cp->SetRatioYTitleSize(0.17,0.47);
    cp->SetPrintRatio(1);

    //cp->SetRatioYRangeC(1,0.04);
    cp->SetRatioYRange(set_ratio_y_min, set_ratio_y_max);
    cp->SetLogx();
    //h1->GetYaxis()->SetTitleOffset(1.47);
    cp->AddHist1D(h1, label1,"LPE1",kBlack,1,0);
    cp->AddHist1D(h2, label2,"LPE1",kBlue ,2,0);

    cy=new TCanvas("cy","cy",600,700);
    cp->Prepare2Pads(cy);
    cp->Draw(cy);
    cp->TransLegend(-0.2,-0.6);
    cp->WidenLegend(0.2,0.);
    cy->Update();
  }

  if (cy && fnameTag.Length() && saveDirTag.Length()) {
    if (doSave) {
      TString fname=TString("fig") + saveDirTag;
      fname.Append(Form("%dD",(DYTools::study2D) ? 2:1));
      fname.Append(fnameTag);
      TString locOutDir=TString("plots") + saveDirTag;
      locOutDir.ReplaceAll("--","");
      SaveCanvas(cy,fname,locOutDir);
    }
    else {
      std::cout << "not saved due to a request\n";
    }
  }

  //h2diff->Print("range");
  //h2esf1->Print("range");

}
