// Created on Feb 12, 2014. Adapted compareEff.C

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

int loadSF(TString fname, TString label, 
	   TMatrixD **sf_out, 
	   TMatrixD **sfErrLo_out, TMatrixD **sfErrHi_out) {
  TFile fin(fname,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << fname << ">\n";
    return 0;
  }
  TString field=Form("scale_factor_%s",label.Data());
  TString fieldErr=field + TString("_err");
  TMatrixD* sf=(TMatrixD*) fin.Get(field);
  TMatrixD* sfErr=(TMatrixD*) fin.Get(fieldErr);
  fin.Close();

  if (!sf || !sfErr) {
    std::cout << "loadSF(" << fname << ", label=" << label << ") error\n";
    if (!sf) std::cout << " - failed to load <" << field << ">\n";
    if (!sfErr) std::cout << " - failed to load <" << fieldErr << ">\n";
    return 0;
  }

  *sf_out   = sf;
  *sfErrLo_out= sfErr;
  *sfErrHi_out= new TMatrixD(*sfErr);

  return 1;
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

int loadEGammaEff(TString fname, TString kindStr, TMatrixD **eff, TMatrixD **effLo, TMatrixD **effHi) {
  TString field=(kindStr.Index("mc")!=-1) ? "eff_mc" : "eff_data";
  if (kindStr==TString("sf")) field="sf";
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

TString effDataKindString(const TString str) {
  TString effKind="xx", dataKind="xx";
  if (str.Index("RECO")!=-1) effKind="RECO";
  else if (str.Index("ID")!=-1) effKind="ID";
  else if (str.Index("HLTleg1")!=-1) effKind="HLTleg1";
  else if (str.Index("HLTleg2")!=-1) effKind="HLTleg2";
  else if (str.Index("HLT")!=-1) effKind="HLT";
  if (str.Index("data")!=-1) dataKind="data";
  else if (str.Index("mc")!=-1) dataKind="mc";
  else if (str.Index("sf")!=-1) dataKind="";
  TString final=dataKind+TString(" ")+effKind;
  return final;
}

// ------------------------------------------------------------
// ------------------------------------------------------------

void compareSF(int iBr=0, int iEta=0, 
	       int doSave=0,
	       double ratioTitleOffset=0.58,
	       double transLegendY_user=0.) {


  TString path1, path2; 
  TString fnameBase="scale_factors_";

  TString label1,label2;
  TString fnameTag;
  TString sfKindLongStr1,sfKindLongStr2;
  TString fileTag1, fileTag2;
  TString sfKind;

  TString path3, sfKindLongStr3,label3,fileTag3;

  TString egammaFName;

  int HLTcomparison=0;
  double transLegendX=-0.2;
  double transLegendY=-0.4;
  int exchange12=0;

  if (1) { // compare to EGamma
    path1="/home/andriusj/cms/DYee-20131024/root_files_reg/constants/DY_j22_19712pb/"; 
    path1="./";
    path2="./";
    fnameBase="el-effs-";
    egammaFName="mediumID.root";
    fileTag1="ETBINS6ETABINS5corr2";
    //fileTag2="EGAMMA";
    sfKindLongStr1="sf_ID_ETBINS6ETABINS5corr2";
    sfKindLongStr2="sf_ID_EtBins6EtaBins5_EGamma";
    label1="DYee |#eta|<2.4";
    label2="EGamma";
    if (1) {
      path3="./";
      fileTag3="ETBINS6ETABINS5";
      sfKindLongStr3="sf_ID_ETBINS6ETABINS5";
      label3="DYee |#eta|<2.5";
      exchange12=1;
    }
    fnameTag="-cmpEGamma-LABEL";
    transLegendX=-0.1;
  }

  if (0) {
    path1="./";
    path2="./";
    fnameBase="el-effs-";
    fileTag1="ETBINS6ETABINS5-summer2013";
    fileTag2="ETBINS6ETABINS5corr2";
    sfKindLongStr1="sf_ID_ETBINS6ETABINS5";
    sfKindLongStr2="sf_ID_ETBINS6ETABINS5corr2";
    label1="Summer2013";
    label2="DYee |#eta|<2.4";
    if (1) {
      path3="./";
      fileTag3="ETBINS6ETABINS5";
      sfKindLongStr3="sf_ID_ETBINS6ETABINS5";
      label3="DYee |#eta|<2.5";
      exchange12=0;
    }
    fnameTag="-cmpSummer2013-LABEL";
    transLegendX=-0.1;
    //transLegendY=-0.2;
  }

  if (transLegendY_user!=0.) transLegendY=transLegendY_user;

  if (iBr==0) sfKind="ID";
  else if (iBr==1) sfKind="RECO";
  else if (iBr==2) sfKind="HLT";
  else {
    std::cout << "iBr error\n";
    return;
  }

  fnameTag.ReplaceAll("LABEL",sfKind);
  if (iBr!=0) {
    sfKindLongStr1.ReplaceAll("ID",sfKind);
    sfKindLongStr2.ReplaceAll("ID",sfKind);
    sfKindLongStr3.ReplaceAll("ID",sfKind);
  }

  DYTools::TEtBinSet_t etBinSet1=DetermineEtBinSet(sfKindLongStr1);
  DYTools::TEtaBinSet_t etaBinSet1=DetermineEtaBinSet(sfKindLongStr1);
  DYTools::TEtBinSet_t etBinSet2=DetermineEtBinSet(sfKindLongStr2);
  DYTools::TEtaBinSet_t etaBinSet2=DetermineEtaBinSet(sfKindLongStr2);
  std::cout << "sets: "<< EtBinSetName(etBinSet1) << "," << EtaBinSetName(etaBinSet1) << "  " << EtBinSetName(etBinSet2) << "," << EtaBinSetName(etaBinSet2) << "\n";

  TString effKind =effDataKindString(sfKindLongStr1);
  TString effKind2=effDataKindString(sfKindLongStr2);
  if (effKind != effKind2) {
    if ( !efficiencyIsHLT(DetermineEfficiencyKind(effKind )) ||
	 !efficiencyIsHLT(DetermineEfficiencyKind(effKind2)) ) {
      std::cout << "effKind1=<" << effKind << ">\n";
      std::cout << "effKind2=<" << effDataKindString(sfKindLongStr2) << ">\n";
      return;
    }
  }
  TString dataKind=effKind;

  TMatrixD *sf1=NULL, *sf1ErrLo=NULL, *sf1ErrHi=NULL;
  TMatrixD *sf2=NULL, *sf2ErrLo=NULL, *sf2ErrHi=NULL;
  TH1D *histo1=NULL, *histo2=NULL;
  const char *histo1Name="histo1";
  const char *histo2Name="histo2";

  TString fname1=path1 + fnameBase + fileTag1 + TString(".root");
  TString fname2=path2 + fnameBase + fileTag2 + TString(".root");

  if (!loadSF(fname1,sfKind,&sf1,&sf1ErrLo,&sf1ErrHi)) {
    std::cout << "failed to get fields from <" << fname1 << "> (1)\n"; 
    return ;
  }
  HERE("loadSF ok");

  if (label2 == TString("EGamma")) {
    fname2=egammaFName;
    if (!loadEGammaEff(fname2,"sf",&sf2,&sf2ErrLo,&sf2ErrHi)) {
      std::cout << "failed to load EGammaSf\n";
      return;
    }
    HERE("load egamma ok");
  }
  else {
    if (!loadSF(fname2,sfKind,&sf2,&sf2ErrLo,&sf2ErrHi)) {
      std::cout << "failed to get fields from <" << fname2 << "> (2)\n"; 
      return ;
    }
  }

  HERE("create graphs");

  TGraphAsymmErrors* gr1=getAsymGraph(etBinSet1,etaBinSet1,iEta,*sf1,*sf1ErrLo,*sf1ErrHi,&histo1,histo1Name);
  gr1->Print("range");

  TGraphAsymmErrors* gr2=getAsymGraph(etBinSet2,etaBinSet2,iEta,*sf2,*sf2ErrLo,*sf2ErrHi,&histo2,histo2Name);
  gr2->Print("range");

  //TGraphAsymmErrors* div=(TGraphAsymmErrors*)gr1->Clone("div");
  TH1D *div=(TH1D*)histo1->Clone("div");
  div->Divide(histo1,histo2,1.,1.,"b");
  div->Print("range");

  TH1D *histo3=NULL;
  TGraphAsymmErrors* gr3=NULL;
  TH1D* div31=NULL;

  if (sfKindLongStr3.Length() && fileTag3.Length() && label3.Length()) {
    TString fname3= path3 + fnameBase + fileTag3 + TString(".root");
    TMatrixD *sf3=NULL, *sf3ErrLo=NULL, *sf3ErrHi=NULL;
    if (!loadSF(fname3,sfKind,&sf3,&sf3ErrLo,&sf3ErrHi)) {
      std::cout << "failed to get field from <" << fname3 << "> (3)\n";
      return ;
    }
    DYTools::TEtBinSet_t etBinSet3=DetermineEtBinSet(sfKindLongStr3);
    DYTools::TEtaBinSet_t etaBinSet3=DetermineEtaBinSet(sfKindLongStr3);
   
    TString histo3Name="histo3";
    gr3=getAsymGraph(etBinSet3,etaBinSet3,iEta,*sf3,*sf3ErrLo,*sf3ErrHi, &histo3,histo3Name);
    gr3->Print("range");
    delete sf3;
    delete sf3ErrLo;
    delete sf3ErrHi;

    div31=(TH1D*)histo1->Clone("div31");
    div31->SetTitle("div31");
    if (HLTcomparison) div31->Divide(histo1,histo3,1.,1.,"b");
    else {
      if (!exchange12) div31->Divide(histo3,histo1,1.,1.,"b");
      else div31->Divide(histo2,histo3,1.,1.,"b");
    }
    div31->Print("range");
    div31->SetLineColor(kGreen+1);
    div31->SetMarkerColor(kGreen+1);

    div->Divide(histo2,histo1,1.,1.,"b");
    div->SetLineColor(kBlue);
    div->SetMarkerColor(kBlue);
  }

  double *loc_etaBinLimits=DYTools::getEtaBinLimits(etaBinSet1);
  int signedEta=DYTools::signedEtaBinning(etaBinSet1);
  TString cpTitle=dataKind+ TString(Form(" %5.3lf #leq %s #leq %5.3lf",loc_etaBinLimits[iEta],(signedEta)?"#eta":"abs(#eta)",loc_etaBinLimits[iEta+1]));

  ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,"comp",cpTitle,
		      "E_{T}",effKind + TString(" scale factor"),"ratio");
  cp.SetLogx();
  cp.AddLine(10.,1.,500.,1.,kBlack,2);

  if (gr3 && HLTcomparison) { // for HLT sficiency
    cp.SetYRange(0.0,1.02);
    /*
    if (DetermineDataKind(sfKind)==DYTools::DATA) {
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

  if (gr3 && !HLTcomparison && exchange12) {
    std::cout << "\n\tInverted plotting order 2,1\n";
    gr2->GetYaxis()->SetTitleOffset(1.4);
    //gr1->SetMarkerStyle(24);
    div->SetMarkerStyle(24);
    cp.AddGraph(gr2,label2,"LPE1",kBlue);
    cp.AddGraph(gr1,label1," PE1",kBlack,24);
  }
  else {
    cp.AddGraph(gr1,label1,"LPE1",kBlack);
    cp.AddGraph(gr2,label2," PE1",kBlue,24);
  }
  if (gr3) {
    //gr3->SetMarkerStyle(27);
    div31->SetMarkerStyle(27);
    cp.AddGraph(gr3,label3," PE1",kGreen+1,27);
  }

  cp.SetRefIdx(-1); // not not plot ratios
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
      div->GetYaxis()->SetRangeUser(0.99,1.01);
      if (iEta==2) div->GetYaxis()->SetRangeUser(0.9,1.1);
    }
  }
  div->Draw("LP");
  if (div31) {
    div31->Draw("LP same");
  }

  TLine *lineAtOne =   new TLine(10,1., 500,1.);
  lineAtOne->SetLineStyle(kDashed);
  lineAtOne->SetLineWidth(1);
  lineAtOne->SetLineColor(kBlack);
  lineAtOne->Draw();

  cx->Update();

  // Save file
  if (fnameTag.Length()) {
    TString fname=TString("fig-sf") + cpTitle;
    fname.ReplaceAll(" #leq "," ");
    fname.ReplaceAll(" ","_");
    fname.ReplaceAll("(#eta)","Eta");
    fname.ReplaceAll("#eta","eta");
    fname.ReplaceAll(".","_");
    //fname.Append(".png");
    std::cout << "fname=" << fname << "\n";

    TString locOutDir=TString("plots") + fnameTag;
    if (doSave) {
      locOutDir.ReplaceAll("--","");
      SaveCanvas(cx,fname,locOutDir);
    }
    else {
      std::cout << "... canvas not saved, as requested\n";
      std::cout << "   locOutDir=" << locOutDir << "\n";
    }
  }

  return ;
}
