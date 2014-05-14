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

TGraphAsymmErrors* calculateSF(TString dataFName, int vsEt, double targetVal,
			       TH1D** histo, TString histoName) {
  int weighted=(dataFName.Index("count-count")!=-1) ? -1 : 0;

  TString mcFName=dataFName;
  mcFName.ReplaceAll("data","mc");
  mcFName.ReplaceAll("fit-fit","count-count");

  TMatrixD *eff1=NULL, *eff1ErrLo=NULL, *eff1ErrHi=NULL;
  TMatrixD *eff2=NULL, *eff2ErrLo=NULL, *eff2ErrHi=NULL;

  if (!loadEff(dataFName,weighted,&eff1,&eff1ErrLo,&eff1ErrHi)) {
    std::cout << "failed to get fields from <" << dataFName << "> (1)\n";
    return NULL;
  }

  if (!loadEff(mcFName,1,&eff2,&eff2ErrLo,&eff2ErrHi)) {
    std::cout << "failed to get fields from <" << mcFName << "> (2)\n";
    return NULL;
  }

  TString effKindLongStr1=dataFName;
  DYTools::TEtBinSet_t etBinSet1=DetermineEtBinSet(effKindLongStr1);
  DYTools::TEtaBinSet_t etaBinSet1=DetermineEtaBinSet(effKindLongStr1);
  int loc_etBinCount=DYTools::getNEtBins(etBinSet1);
  int loc_etaBinCount=DYTools::getNEtaBins(etaBinSet1);
  double *loc_etBinLimits=DYTools::getEtBinLimits(etBinSet1);
  double *loc_etaBinLimits=DYTools::getEtaBinLimits(etaBinSet1);

  int use_binCount=(vsEt) ? loc_etBinCount : loc_etaBinCount;
  double *use_binLimits=(vsEt) ? loc_etBinLimits : loc_etaBinLimits;

  int etMin=0, etMax=loc_etBinCount;
  int etaMin=0, etaMax=loc_etaBinCount;
  TString grTitle;
  if (vsEt) {
    grTitle=Form("eta_%2.1lf",targetVal);
    etaMin=DYTools::_findMassBin(targetVal,loc_etaBinCount,loc_etaBinLimits);
    if (etaMin==-1) {
      std::cout << "failed to get etaMin value\n";
      return NULL;
    }
    etaMax=etaMin+1;
  }
  else {
    grTitle=Form("Et_%2.1lf",targetVal);
    etMin=DYTools::_findMassBin(targetVal,loc_etBinCount,loc_etBinLimits);
    if (etMin==-1) {
      std::cout << "failed to get etMin value\n";
      return NULL;
    }
    etMax=etMin+1;
  }

  eliminateSeparationSigns(histoName);
  (*histo) = new TH1D(histoName, histoName, use_binCount,use_binLimits);
  (*histo)->SetDirectory(0);

  double x[use_binCount], dx[use_binCount];
  double sf[use_binCount], sfLo[use_binCount], sfHi[use_binCount];

  int iPt=0;
  for (int iEt=etMin; iEt<etMax; ++iEt) {
    for (int iEta=etaMin; iEta<etaMax; iEta++, iPt++) {
      x[iPt] = 0.5*(use_binLimits[iPt] + use_binLimits[iPt+1]);
      dx[iPt]= 0.5*(use_binLimits[iPt+1] - use_binLimits[iPt]);

      double effData=(*eff1)(iEt,iEta);
      double effDataErr= 0.5*((*eff1ErrLo)(iEt,iEta) + (*eff1ErrHi)(iEt,iEta));
      double effMC=(*eff2)(iEt,iEta);
      double effMCErr= 0.5*((*eff2ErrLo)(iEt,iEta) + (*eff2ErrHi)(iEt,iEta));
      sf[iPt] = effData/effMC;
      sfLo[iPt]=errOnRatio(effData,effDataErr, effMC,effMCErr);
      sfHi[iPt]=sfLo[iPt];

      (*histo)->SetBinContent(iPt+1, sf[iPt]);
      (*histo)->SetBinError  (iPt+1, sfLo[iPt]);
    }
  }

  delete loc_etBinLimits;
  delete loc_etaBinLimits;

  TGraphAsymmErrors* gr= new TGraphAsymmErrors(use_binCount,
					       x,sf,dx,dx,sfLo,sfHi);
  gr->SetTitle(grTitle);

  return gr;
}

// ------------------------------------------------------------
// ------------------------------------------------------------

// Divide two histograms, assuming that they have overlapping ranges

TGraphAsymmErrors* sfRatioGraph(const TH1D* sf1, const TH1D *sf2,
				TString effString1, TString effString2,
				int vsEt)
{

  DYTools::TEtBinSet_t etBinSet1=DetermineEtBinSet(effString1);
  DYTools::TEtaBinSet_t etaBinSet1=DetermineEtaBinSet(effString1);
  DYTools::TEtBinSet_t etBinSet2=DetermineEtBinSet(effString2);
  DYTools::TEtaBinSet_t etaBinSet2=DetermineEtaBinSet(effString2);

  int useBCount1= (vsEt) ?
    DYTools::getNEtBins(etBinSet1) : DYTools::getNEtaBins(etaBinSet1);
  int useBCount2= (vsEt) ?
    DYTools::getNEtBins(etBinSet2) : DYTools::getNEtaBins(etaBinSet2);
  double *loc_bins1= (vsEt) ?
    DYTools::getEtBinLimits(etBinSet1) : DYTools::getEtaBinLimits(etaBinSet1);
  double *loc_bins2= (vsEt) ?
    DYTools::getEtBinLimits(etBinSet2) : DYTools::getEtaBinLimits(etaBinSet2);

  std::vector<double> xVals;
  std::vector<int> yIdx1, yIdx2;
  xVals.reserve(useBCount1+useBCount2);
  yIdx1.reserve(useBCount1+useBCount2);
  yIdx2.reserve(useBCount1+useBCount2);

  double x=loc_bins1[0];
  if (x<loc_bins2[0]) {
    std::cout << "myDivide: assumption about lower edge is not correct\n";
    return NULL;
  }
  if (loc_bins1[useBCount1] != loc_bins2[useBCount2]) {
    std::cout << "sfRatioGraph: assumption about last edge is not correct\n";
    for (int i=0; i<=useBCount1; ++i) std::cout << " " <<  loc_bins1[i];
    std::cout << "\n";
    for (int i=0; i<=useBCount2; ++i) std::cout << " " <<  loc_bins2[i];
    std::cout << "\n";
    if (fabs(loc_bins1[useBCount1]-loc_bins2[useBCount2])<1e-3) {
      std::cout << "autocorrection: changing loc_bins1 upper limit\n";
      loc_bins1[useBCount1]=loc_bins2[useBCount2];
    }
    else return NULL;
  }

  xVals.push_back(x);
  yIdx1.push_back(1);
  yIdx2.push_back(1);

  int finished=0;
  int i1=0, i2=0;

  //std::cout << "x=" << x << ", i1=" << i1 << ", i2=" << i2 << "\n";
  while (!finished) {
    //std::cout << "x=" << x << ", i1=" << i1 << ", i2=" << i2 << "\n";
    //std::cout << "x=" << x << ", next x=" << loc_bins1[i1+1] << " " << loc_bins2[i2+1] << "\n";
    if (loc_bins1[i1+1] > loc_bins2[i2+1]) {
      x=loc_bins2[i2+1];
      i2++;
      yIdx1.push_back(i1+1); // +1 to convert to bin number
      yIdx2.push_back(i2+1);
    }
    else {
      x=loc_bins1[i1+1];
      if (loc_bins2[i2+1]==loc_bins1[i1+1]) i2++;
      i1++;
      yIdx1.push_back(i1+1);
      yIdx2.push_back(i2+1);
    }
    xVals.push_back(x);
    finished = ((i1>=useBCount1) && (i2>=useBCount2)) ? 1:0;
  }
  //xVals.push_back(x);

  delete loc_bins1;
  delete loc_bins2;

  int ptCount=int(xVals.size()-1);
  double xv[ptCount], dx[ptCount];
  double rv[ptCount], dr[ptCount];

  //printHisto(sf1);
  //printHisto(sf2);

  for (int i=0; i<ptCount; ++i) {
    xv[i]= 0.5*( xVals[i] + xVals[i+1]);
    dx[i]= 0.5*( xVals[i+1] - xVals[i]);
    double y1=sf1->GetBinContent(yIdx1[i]);
    double y1err=sf1->GetBinError(yIdx1[i]);
    double y2=sf2->GetBinContent(yIdx2[i]);
    double y2err=sf2->GetBinError(yIdx2[i]);
    //std::cout << "i=" << i << ", y1=" << y1 << ", y2=" << y2 << "\n";
    rv[i]= y1/y2;
    dr[i]= errOnRatio(y1,y1err, y2,y2err);
  }

  TGraphAsymmErrors* gr= new TGraphAsymmErrors(ptCount,
					       xv,rv,dx,dx,dr,dr);
  TString grTitle=TString("gr_") + effString1 + TString("_div_") + effString2;
  gr->SetTitle(grTitle);

  return gr;
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


void compareSFv2(int iBr=0, int iBin=0, int vsEt=1,
		 int doSave=0,
		 double transLegendY_user=0.,
		 //double ratioTitleOffset=0.58,
		 TString *outFileName_ptr=NULL,
		 TString *outDir_ptr=NULL,
		 double val_user=0.) {
  TString path1, path2;
  TString effKindLongStr1,effKindLongStr2;
  TString fnameBase="efficiency_TnP_1D_Full2012_";

  TString label1="new n-tuples";
  TString label2="old n-tuples";
  TString fnameTag;

  TString path3="";
  TString effKindLongStr3="";
  TString label3="";

  double transLegendX=-0.2;
  double transLegendY=-0.4;


  if (0) { // 2014.05.10
    path1="Results-DYee-effBinStudy/root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy/";
    path2=path1; path3=path1;
    effKindLongStr1="dataRECO_fit-fitEtBins6systEtaBins5_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins6systEtaBins7_PU";
    effKindLongStr3="dataRECO_fit-fitEtBins6systEtaBins9_PU";
    label1="Et6Eta5";
    label2="Et6Eta7";
    label3="Et6Eta9";
    fnameTag="-sf-et6-vars--";
    transLegendX=-0.4;
  }

  if (0) { // 2014.05.10
    path1="Results-DYee-effBinStudy/root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy/";
    path2=path1;
    effKindLongStr1="dataRECO_fit-fitEtBins6systEtaBins5_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins6systEtaBins5_max25_PU";
    label1="Et6Eta5";
    label2="Et6Eta5_max25";
    fnameTag="-sf-et6eta5-vars--";
    transLegendX=-0.4;
  }

  if (0) { // 2014.05.10
    path1="Results-DYee-effBinStudy/root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy/";
    path2=path1;
    effKindLongStr1="dataRECO_fit-fitEtBins6systEtaBins5_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins7altEtaBins5_PU";
    label1="Et6Eta5";
    label2="Et7altEta5";
    fnameTag="-sf-eta5-vars--";
    //transLegendX=-0.4;
    transLegendX=-0.07;
  }

  if (0) { // 2014.05.10
    path1="Results-DYee-effBinStudy/root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy/";
    path2=path1;
    effKindLongStr1="dataRECO_fit-fitEtBins6systEtaBins7_PU";
    effKindLongStr2="dataRECO_fit-fitEtBins7altEtaBins7_PU";
    label1="Et6Eta7";
    label2="Et7altEta7";
    fnameTag="-sf-eta7-vars--";
    //transLegendX=-0.4;
    transLegendX=-0.07;
  }

  if (1) { // 2014.05.10
    path1="Results-DYee-effBinStudy/root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy/";
    path2=path1; path3=path1;
    effKindLongStr2="dataRECO_fit-fitEtBins7altEtaBins7_PU";
    effKindLongStr1="dataRECO_fit-fitEtBins6systEtaBins5_PU";
    effKindLongStr3="dataRECO_fit-fitEtBins6systEtaBins9_PU";
    label1="Et7Eta7";
    label2="Et6Eta5";
    label3="Et6Eta9";
    fnameTag="-sf-cmpBase77-vars--";
    transLegendX=-0.4;
  }


  // -------------------------------
  // processing
  // -------------------------------

  if (transLegendY_user!=0.) transLegendY=transLegendY_user;

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
    effKindLongStr1.ReplaceAll("dataRECO_fit-fit","dataHLTleg1_count-count");
    effKindLongStr2.ReplaceAll("dataRECO_fit-fit","dataHLTleg1_count-count");
    effKindLongStr3.ReplaceAll("dataRECO_fit-fit","dataHLTleg1_count-count");
  }
  else if (iBr==4) {
    effKindLongStr1.ReplaceAll("dataRECO_fit-fit","dataHLTleg2_count-count");
    effKindLongStr2.ReplaceAll("dataRECO_fit-fit","dataHLTleg2_count-count");
    effKindLongStr3.ReplaceAll("dataRECO_fit-fit","dataHLTleg2_count-count");
  }
  else {
    std::cout << "iBr error\n";
    return;
  }

  DYTools::TEtBinSet_t etBinSet1=DetermineEtBinSet(effKindLongStr1);
  DYTools::TEtaBinSet_t etaBinSet1=DetermineEtaBinSet(effKindLongStr1);
  //int loc_etBinCount=DYTools::getEtBinCount(etBinSet1);
  //int loc_etaBinCount=DYTools::getEtaBinCount(etaBinSet1);
  double *loc_etBinLimits=DYTools::getEtBinLimits(etBinSet1);
  double *loc_etaBinLimits=DYTools::getEtaBinLimits(etaBinSet1);

  double *target_BinLimits=(vsEt) ? loc_etaBinLimits : loc_etBinLimits;
  double targetVal= 0.5 * ( target_BinLimits[iBin] + target_BinLimits[iBin+1] );
  if (val_user!=double(0.)) {
    std::cout << "user supplied value = " << val_user << "\n";
    targetVal=val_user;
  }
  std::cout << "iBin=" << iBin << ", targetVal=" << targetVal << "\n";


  TString fname1=path1 + fnameBase + effKindLongStr1 + TString(".root");
  TString fname2=path2 + fnameBase + effKindLongStr2 + TString(".root");
  TString fname3;
  if (effKindLongStr3.Length()) {
    fname3= path3 + fnameBase + effKindLongStr3 + TString(".root");
  }

  TGraphAsymmErrors *gr1=NULL, *gr2=NULL, *gr3=NULL;
  TH1D *h1=NULL, *h2=NULL, *h3=NULL;

  gr1=calculateSF(fname1, vsEt, targetVal, &h1, TString("histo_") + label1);
  gr2=calculateSF(fname2, vsEt, targetVal, &h2, TString("histo_") + label2);
  if (fname3.Length()) {
    gr3=calculateSF(fname3, vsEt, targetVal, &h3, TString("histo_") + label3);
  }

  TGraphAsymmErrors* div21=sfRatioGraph(h2,h1,effKindLongStr2,effKindLongStr1,vsEt);
  if (!div21) { std::cout << "error detected\n"; return; }
  div21->Print("range");

  TGraphAsymmErrors* div31=NULL;
  if (h3) {
    div31=sfRatioGraph(h3,h1,effKindLongStr3,effKindLongStr1,vsEt);
    if (!div31) { std::cout << "error detected\n"; return; }
    div31->Print("range");
  }

  TString effKindName=EfficiencyKindName(DetermineEfficiencyKind(effKindLongStr1));
  TString targetTitle=(vsEt) ? "|#eta|" : "ET";
  TString xaxisTitle =(vsEt) ? "#it{E}_{T} [GeV]" : "|eta|";
  TString cpTitle= effKindName +
    TString(Form(" %s = %4.2lf",targetTitle.Data(), targetVal));

  ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,"comp",cpTitle,
		      xaxisTitle,"scale factor","ratio");
  cp.SetRefIdx(-111); // do not plot lower panel
  if (vsEt) cp.SetLogx();

  TCanvas *cx=new TCanvas("cx","cx",800,800);
  cp.Prepare2Pads(cx);

  gr1->GetYaxis()->SetTitleOffset(1.4);
  gr1->GetXaxis()->SetLabelOffset(2.);

  cp.AddGraph(gr1,label1,"LPE1",kBlue);
  cp.AddGraph(gr2,label2," PE1",kBlack,24);
  if (gr3) cp.AddGraph(gr3,label3," PE1",kGreen+2,27);

  cp.Draw(cx,0,"png",1);
  cp.TransLegend(transLegendX, transLegendY);
  cp.WidenLegend(0.2,0.);


  ComparisonPlot_t cpRatio(ComparisonPlot_t::_ratioPlain,"compRatio","",
			   xaxisTitle,"ratio","nothing");
  cpRatio.NoLegend(1);
  cpRatio.SetRefIdx(-111); // no lower panel
  if (vsEt) cpRatio.SetLogx();
  cpRatio.SetXAxisTextSizes(0.15,1.1, 0.11);
  cpRatio.SetYAxisTextSizes(0.15,0.5, 0.11);

  cpRatio.AddGraph(div21,label2,"PE1",kBlack,24);
  if (div31) cpRatio.AddGraph(div31,label3,"PE1",kGreen+2,27);

  TAxis *ax=div21->GetXaxis();
  int lineColor=(div31) ? kRed+1 : kBlue;
  cpRatio.AddLine(ax->GetBinLowEdge(1),1.,
		  ax->GetBinUpEdge(ax->GetNbins()),1.,lineColor,2);

  cpRatio.Draw(cx,0,"png",2);

  if (0) {
    double one=1;
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
    fname.ReplaceAll("|#eta|","absEta");
    fname.ReplaceAll("(#eta)","Eta");
    fname.ReplaceAll("#eta","eta");
    fname.ReplaceAll(".","_");
    fname.ReplaceAll("#it{E}_{T}","Et");
    fname.ReplaceAll("_=_","_");
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
