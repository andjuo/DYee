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
// does not work

TGraphAsymmErrors* divideEffs(const TGraphAsymmErrors *gr1, const TGraphAsymmErrors *gr2) {
  TH1F *h1=gr1->GetHistogram();
  TH1F *h2=gr2->GetHistogram();
  std::cout << "h1="; h1->Print("range");
  std::cout << "h2="; h2->Print("range");
  TGraphAsymmErrors *div=(TGraphAsymmErrors*)gr1->Clone("temp");
  div->Divide(h1,h2,"cl=0.683 b(1,1) mode");
  return div;
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
  TString final=dataKind+TString(" ")+effKind;
  return final;
}

// ------------------------------------------------------------


void compareEff(TString effKindLongStr1="mcRECO_count-countEtBins6EtaBins5_PU", 
		TString effKindLongStr2="mcRECO_count-countEtBins6EtaBins5_PU",
		int iEta=0) {
  TString path1="/home/andriusj/cms/DYee8TeV-20140118/root_files/tag_and_probe/DY_j22_19712pb/";
  TString path2="/home/andriusj/cms/DYee-20131024/root_files_reg/tag_and_probe/DY_j22_19712pb/";
  //path2="/home/andriusj/cms/CMSSW_3_8_4/src/DYee8TeV-20130801/DrellYanDMDY/root_files/tag_and_probe/DY_j22_19789pb/";
  path2="/home/andriusj/cms/DYee8TeV-20140118-maxEta24/root_files/tag_and_probe/DY_j22_19712pb/";
  TString fnameBase="efficiency_TnP_1D_Full2012_";

  TString fname1=path1 + fnameBase + effKindLongStr1 + TString(".root");
  TString fname2=path2 + fnameBase + effKindLongStr2 + TString(".root");

  DYTools::TEtBinSet_t etBinSet1=DYTools::ETBINS6;
  DYTools::TEtaBinSet_t etaBinSet1=DYTools::ETABINS5;
  DYTools::TEtBinSet_t etBinSet2=DYTools::ETBINS6;
  DYTools::TEtaBinSet_t etaBinSet2=DYTools::ETABINS5;

  TString label1="new n-tuples";
  TString label2="old n-tuples";
  TString fnameTag="-new_vs_old--";

  etaBinSet2=DYTools::ETABINS5corr;
  label1="etaMax = 2.5";
  label2="etaMax = 2.4";
  fnameTag="-diffEtaMax--";


  TString effKind=effDataKindString(effKindLongStr1);
  if (effKind != effDataKindString(effKindLongStr2)) {
    std::cout << "effKind1=<" << effKind << ">\n";
    std::cout << "effKind2=<" << effDataKindString(effKindLongStr2) << ">\n";
    return;
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

  if (!loadEff(fname2,weighted2,&eff2,&eff2ErrLo,&eff2ErrHi)) {
    std::cout << "failed to get fields from <" << fname2 << "> (2)\n"; 
    return ;
  }

  TGraphAsymmErrors* gr1=getAsymGraph(etBinSet1,etaBinSet1,iEta,*eff1,*eff1ErrLo,*eff1ErrHi,&histo1,histo1Name);
  gr1->Print("range");

  TGraphAsymmErrors* gr2=getAsymGraph(etBinSet2,etaBinSet2,iEta,*eff2,*eff2ErrLo,*eff2ErrHi,&histo2,histo2Name);
  gr2->Print("range");

  //TGraphAsymmErrors* div=(TGraphAsymmErrors*)gr1->Clone("div");
  TH1D *div=(TH1D*)histo1->Clone("div");
  div->Divide(histo1,histo2,1.,1.,"b");
  div->Print("range");


  ComparisonPlot_t cpTemp(ComparisonPlot_t::_ratioPlain,"comp","comp",
			  "E_{T}","eff","ratio");

  double *loc_etaBinLimits=DYTools::getEtaBinLimits(etaBinSet1);
  int signedEta=DYTools::signedEtaBinning(etaBinSet1);
  TString cpTitle=dataKind+ TString(Form(" %5.3lf #leq %s #leq %5.3lf",loc_etaBinLimits[iEta],(signedEta)?"#eta":"abs(#eta)",loc_etaBinLimits[iEta+1]));

  CPlot cp("comp",cpTitle,"E_{T}","efficiency");
  cp.SetLogx();
  TCanvas *cx=new TCanvas("cx","cx",600,700);
  cpTemp.Prepare2Pads(cx);

  gr1->GetYaxis()->SetTitleOffset(1.4);

  cp.AddGraph(gr1,label1,"LP",kBlack);
  cp.AddGraph(gr2,label2,"LP",kBlue);
  cp.Draw(cx,0,"png",1);
  cp.TransLegend(0, -0.4);
  cx->cd(2);
  cx->GetPad(2)->SetLogx(cp.fLogx);
  div->SetTitle("");
  div->GetXaxis()->SetTitle("E_{T}");
  div->GetXaxis()->SetTitleSize(0.17);
  div->GetXaxis()->SetLabelSize(0.17);
  div->GetXaxis()->SetNoExponent();
  div->GetXaxis()->SetMoreLogLabels();

  div->GetYaxis()->SetTitle("ratio");
  div->GetYaxis()->SetTitleOffset(0.7);
  div->GetYaxis()->SetTitleSize(0.1);
  div->Draw("LP");

  cx->Update();

  TString fname=TString("fig-eff-") + fnameTag + cpTitle;
  fname.ReplaceAll(" #leq "," ");
  fname.ReplaceAll(" ","_");
  fname.ReplaceAll("(#eta)","Eta");
  fname.ReplaceAll("#eta","eta");
  fname.ReplaceAll(".","_");
  //fname.Append(".png");
  std::cout << "fname=" << fname << "\n";

  SaveCanvas(cx,fname);

  return ;
}
