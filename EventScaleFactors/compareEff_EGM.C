//
// Made on May 07, 2014 from compareEff.C. Make a plot to match EG-13-001 paper
// Made on March 14, 2014
// Enhanced version of compareEff.C. Allows plotting vs eta
//

#include "../Include/DYTools.hh"
#include "calcEventEffLink.h"
#include <TGraphAsymmErrors.h>
#include "../Include/ComparisonPlot.hh"
#include "../Include/colorPalettes.hh"

// ------------------------------------------------------------
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


void compareEff_EGM(int iBr=0, int iBin=0, int vsEt=1,
		     int doSave=0,
		     double transLegendY_user=0.,
		     double ratioTitleOffset=0.58,
		     TString *outFileName_ptr=NULL,
		     TString *outDir_ptr=NULL) {
  TString path1, path2;
  TString effKindLongStr1,effKindLongStr2;
  TString fnameBase="efficiency_TnP_1D_Full2012_";

  TString label1;
  TString label2;
  TString fnameTag;

  TString path3="";
  TString effKindLongStr3="";
  TString label3="";

  double transLegendX=-0.2;
  double transLegendY=-0.6;


  if (1) { // 2014.05.07
    path1="../root_files_reg/tag_and_probe/DY_j22_19712pb_egamma_Unregressed_energy/";
    path2= path1;
    effKindLongStr1="dataRECO_fit-fitEtBins6EtaBins5egamma_PU";
    effKindLongStr2="mcRECO_count-countEtBins6EtaBins5egamma_PU";
    label1="data";
    label2="simulation";
    fnameTag="-reco";
    transLegendX=0.1;
  }



  // -------------------------------
  // processing
  // -------------------------------

  if (transLegendY_user!=0.) transLegendY=transLegendY_user;

  if ((iBr==0) || (iBr==2)) {
    if (iBr==2) fnameTag.Append("extraSyst");
  }
  else if ((iBr==1) || (iBr==3)) {
    effKindLongStr1.ReplaceAll("RECO","ID");
    effKindLongStr2.ReplaceAll("RECO","ID");
    fnameTag="-id";
    if (iBr==3) fnameTag.Append("extraSyst");
  }
  else {
    std::cout << "unknown branch = " << iBr << "\n";
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
  if (effKind == effKind2) {
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

  if (!loadEff(fname2,weighted2,&eff2,&eff2ErrLo,&eff2ErrHi)) {
    std::cout << "failed to get fields from <" << fname2 << "> (2)\n";
    return ;
  }

  if (etaBinSet1==DYTools::ETABINS5egamma) etaBinSet1=DYTools::ETABINS5_max25;
  if (etaBinSet2==DYTools::ETABINS5egamma) etaBinSet2=DYTools::ETABINS5_max25;
  if (etBinSet1==DYTools::ETBINS6) etBinSet1=DYTools::ETBINS6short;
  if (etBinSet2==DYTools::ETBINS6) etBinSet2=DYTools::ETBINS6short;

  TGraphAsymmErrors* gr1=getAsymGraph(vsEt,etBinSet1,etaBinSet1,iBin,*eff1,*eff1ErrLo,*eff1ErrHi,&histo1,histo1Name);
  //std::cout << gr1->GetTitle() << ": "; gr1->Print("range");

  TGraphAsymmErrors* gr2=getAsymGraph(vsEt,etBinSet2,etaBinSet2,iBin,*eff2,*eff2ErrLo,*eff2ErrHi,&histo2,histo2Name);
  //std::cout << gr2->GetTitle() << ": "; gr2->Print("range");

  if (1) {
    gr1->Print("range");
    gr2->Print("range");
  }


  TGraphAsymmErrors* div=(TGraphAsymmErrors*)gr1->Clone("div");
  //TH1D *div=(TH1D*)histo1->Clone("div");
  //div->Divide(histo1,histo2,1.,1.,"b");
  div->Divide(histo1,histo2,"pois");
  div->Print("range");

  TMatrixD *sfSystErrEgamma=loadMatrix(fname1,"sf_syst_rel_error_egamma",
				       6,5,1);
  TMatrixD *sfSystErr25=loadMatrix(fname1,"sf_syst_rel_error_maxEta25",
				   6,5,1);

  sfSystErrEgamma->Print();
  sfSystErr25->Print();

  TVectorD egmSystErr(6), ourSystErr(6);
  for (int ic=0; ic<6; ++ic) {
    egmSystErr(ic) = (*sfSystErrEgamma)(ic,iBin);
    ourSystErr(ic) = (*sfSystErr25)(ic,iBin);
  }

  TGraphAsymmErrors* sfEG=addErrors(div,egmSystErr,"sfEG",1);
  TGraphAsymmErrors* sfOur=addErrors(sfEG,ourSystErr,"sfOur",1);

  sfEG->Print();

  double *loc_etBinLimits=DYTools::getEtBinLimits(etBinSet1);
  double *loc_etaBinLimits=DYTools::getEtaBinLimits(etaBinSet1);
  int signedEta=DYTools::signedEtaBinning(etaBinSet1);
  TString cpTitle;
  if (vsEt) cpTitle= dataKind+ TString(Form(" %5.3lf #leq %s #leq %5.3lf",loc_etaBinLimits[iBin],(signedEta)?"#eta":"abs(#eta)",loc_etaBinLimits[iBin+1]));
  else cpTitle= dataKind+ TString(Form(" %2.0lf #leq #it{E}_{T} #leq %2.0lf GeV",loc_etBinLimits[iBin],loc_etBinLimits[iBin+1]));
  TString xaxisTitle=(vsEt) ? "#it{E}_{T}" : ((signedEta) ? "#eta" : "|#eta|");


  ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,"comp",cpTitle,
			  xaxisTitle,"efficiency","ratio");
  cp.SetRefIdx(-111); // do not plot lower panel
  //if (vsEt) cp.SetLogx();
  cp.SetLogx(0);
  cp.SetGridx(1);
  cp.SetGridy(1);
  cp.SetYRange(0.0, 1.1);
  cp.SetXRange(0,100);

  // Square canvas if ratio is not plotted
  TCanvas *cx=new TCanvas("cx","cx",700,700);
  if (0) {
    cp.Prepare2Pads(cx);
  }
  else {
    cx->Divide(2,1);
    cp.PreparePads(cx,1,2,"comp",0.29,0.001);
  }

  SetSideSpaces(cx,0,-0.1, 0.12, 0.0,1);
  SetSideSpaces(cx,0,-0.1,-0.08,-0.18,2);

  gr1->GetYaxis()->SetTitleOffset(1.4);
  gr1->GetXaxis()->SetNdivisions(506);

  gr2->SetFillStyle(3001);
  gr2->SetFillColor(38);

  gr1->SetMarkerSize(0.8);
  gr2->SetMarkerSize(0.8);

  cp.AddGraph(gr1,label1,"LPE1",kBlack,20, 1,1,1.);
  cp.AddGraph(gr2,label2,"LPE2",38,22, 1,1,1.);

  cp.Draw(cx,0,"png",1);
  cp.TransLegend(transLegendX, transLegendY);
  //cp.WidenLegend(0.2,0.);

  ComparisonPlot_t cpSF(ComparisonPlot_t::_ratioPlain,"compSF","",
			  xaxisTitle,"scale factor","ratio");
  cpSF.SetRefIdx(-111);
  cpSF.NoLegend(1);
  cpSF.SetGridx(1);
  cpSF.SetGridy(1);
  //cpSF.SetXAxisTextSizes(0.08,1.1,0.08);
  //cpSF.SetYAxisTextSizes(0.08,1.0,0.08);
  cpSF.SetXAxisTextSizes(0.11,1.10,0.11);
  cpSF.SetYAxisTextSizes(0.14,0.57,0.11);

  cpSF.SetYRange(0.8,1.2);
  if ((iBr==1) || (iBr==3)) cpSF.SetYRange(0.7,1.3);
  cpSF.SetXRange(0,100);

  int color=38;
  sfEG->SetLineColor(color);
  sfEG->SetMarkerColor(color);
  sfEG->SetMarkerSize(0.8);
  sfEG->SetFillColor(color);
  sfEG->SetFillStyle(3001);
  sfOur->SetFillStyle(3002);

  if ((iBr==2) || (iBr==3)) {
    cpSF.AddGraph(sfOur,"","LPE2",kRed,1, 1,1,1.);
    TGraphAsymmErrors *sfEG3=(TGraphAsymmErrors*)sfEG->Clone("sfEG_clone3");
    sfEG3->SetFillStyle(3002);
    cpSF.AddGraph(sfEG3,"","LPE2", 9,1, 1,1,1.);
  }
  cpSF.AddGraph(sfEG,"","LPE2",color,20, 1,1,1.);
  TGraphAsymmErrors *sfEG2=(TGraphAsymmErrors*)sfEG->Clone("sfEG_clone");
  cpSF.AddGraph(sfEG2,"","LPE1", 9,20, 1,1,1.);

  cpSF.Draw(cx,0,"png",2);

  // ---------- line at 1

  if (0) {
    double one= 1.;
    TLine *lineAtOne =   new TLine(0,one, 100,one);
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
