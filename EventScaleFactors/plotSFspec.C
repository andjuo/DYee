// Created on Feb 26, 2014. Adapted compareSF.C
// Created on Feb 12, 2014. Adapted compareEff.C

#include "../Include/DYTools.hh"
#include "calcEventEffLink.h"
#include <TGraphAsymmErrors.h>
#include "../Include/ComparisonPlot.hh"
#include "../Include/MitStyleRemix.hh" // SetSideSpaces

// ------------------------------------------------------------

int addInQuadrature(const TMatrixD &sf, const TMatrixD &sfErr, const TMatrixD &sfSystRelErr, TMatrixD &totErr) {
  for (int ir=0; ir<sf.GetNrows(); ++ir) {
    for (int ic=0; ic<sf.GetNcols(); ++ic) {
      double err2=pow(sfErr(ir,ic),2) + pow(sf(ir,ic)*sfSystRelErr(ir,ic),2);
      totErr(ir,ic)=sqrt(err2);
    }
  }
  return 1;
}

// ----------------------------------------------------------------

TMatrixD* loadMatrix(const TString &fname, const TString &fieldName, int expect_nRows, int expect_nCols, int reportFieldError) {
  TFile f(fname,"read");
  TMatrixD *M=NULL;
  int ok=1;
  if (!f.IsOpen()) ok=0;
  if (ok==1) {
    M=(TMatrixD*)f.Get(fieldName);
    f.Close();
    if (!M) ok=-1;
    else {
      if ((M->GetNrows()!=expect_nRows) ||
	  (M->GetNcols()!=expect_nCols)) {
	ok=-2;
      }
    }
  }
  if (ok!=1) {
    int report=1;
    if ((ok==-1) && !reportFieldError) report=0;
    if (report) {
      std::cout << "Error in loadMatrix(fname=<" << fname << ">, fieldName=<" << fieldName << ">, nRows=" << expect_nRows << ", nCols=" << expect_nCols << "):\n";
      if (ok==0) std::cout << " - failed to open the file\n";
      else if (ok==-1) std::cout << " - failed to load the field\n";
      else if (ok==-2) {
	std::cout << " - size mistmatch. Expect " << expect_nRows << "x" << expect_nCols << ", got " << M->GetNrows() << "x" << M->GetNcols() << "\n";
	delete M;
	M=NULL;
      }
    }
  }
  return M;
}

// ----------------------------------------------------------------



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
/*
// Moved to calcEventEffLink.h

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
*/

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

void plotSFspec(int iBr=0, int iBin=0, 
		int doSave=0,
		double transLegendY_user=0.,
		int vsEt=1) {


  TString fname;

  TString label1,label2,label3;
  TString sfKindLongStr;
  //TString sfKind;

  double transLegendX=-0.2;
  double transLegendY=-0.4;

  TString saveFileTag;

  label1="stat.error";
  label2="stat.+syst.err.";
  label3="stat.+syst.err.+add.syst.err";

  if (iBr==0) {
    fname="./egammaRECO.root";
    sfKindLongStr="sf_RECO_ETBINS6ETABINS5";

    saveFileTag="-cmpEgammaRECO";
    transLegendX=-0.23;
    //transLegendY=-0.2;
  }

  if (iBr==1) {
    fname="./mediumID.root";
    sfKindLongStr="sf_mediumID_ETBINS6ETABINS5";

    saveFileTag="-cmpEgammaMediumID";
    transLegendX=-0.23;
    //transLegendY=-0.2;
  }

  //if (!vsEt) saveFileTag.Append("-vsEta");

  if (transLegendY_user!=0.) transLegendY=transLegendY_user;

  /*
  if (iBr==0) sfKind="ID";
  else if (iBr==1) sfKind="RECO";
  else if (iBr==2) sfKind="HLT";
  else {
    std::cout << "iBr error\n";
    return;
  }
  */

  DYTools::TEtBinSet_t etBinSet=DetermineEtBinSet(sfKindLongStr);
  DYTools::TEtaBinSet_t etaBinSet=DetermineEtaBinSet(sfKindLongStr);
  int nEtBins=DYTools::getNEtBins(etBinSet);
  int nEtaBins=DYTools::getNEtaBins(etaBinSet);

  std::cout << "sets: "<< EtBinSetName(etBinSet) << "," << EtaBinSetName(etaBinSet) << "\n";

  TString effKind =effDataKindString(sfKindLongStr);
  TString dataKind=effKind;

  TMatrixD *sf=NULL, *sf1ErrLo=NULL, *sf1ErrHi=NULL;
  TMatrixD *systRelErr=NULL, sf2ErrLo(nEtBins,nEtaBins), sf2ErrHi(nEtBins,nEtaBins);
  TMatrixD *systRelErrTot=NULL, sf3ErrLo(nEtBins,nEtaBins), sf3ErrHi(nEtBins,nEtaBins);

  // load the scale factors
  if (!loadEGammaEff(fname,"sf",&sf,&sf1ErrLo,&sf1ErrHi)) {
    std::cout << "failed to load EGammaSf\n";
    return;
  }
  HERE("load egamma ok");

  systRelErr=loadMatrix(fname,"sf_syst_rel_error_egamma",nEtBins,nEtaBins,1);
  if (!systRelErr) return;
  systRelErrTot=loadMatrix(fname,"sf_syst_rel_error",nEtBins,nEtaBins,1);
  if (!systRelErrTot) return;

  HERE("add errors");
  addInQuadrature(*sf,*sf1ErrLo, *systRelErr, sf2ErrLo);
  addInQuadrature(*sf,*sf1ErrHi, *systRelErr, sf2ErrHi);
  addInQuadrature(*sf,*sf1ErrLo, *systRelErrTot, sf3ErrLo);
  addInQuadrature(*sf,*sf1ErrHi, *systRelErrTot, sf3ErrHi);

  HERE("create graphs");

  TGraphAsymmErrors* gr1=getAsymGraph(vsEt, etBinSet,etaBinSet,iBin,*sf,*sf1ErrLo,*sf1ErrHi);
  gr1->GetXaxis()->SetMoreLogLabels();
  gr1->GetXaxis()->SetNoExponent();
  gr1->Print("range");

  TGraphAsymmErrors* gr2=getAsymGraph(vsEt, etBinSet,etaBinSet,iBin,*sf,sf2ErrLo,sf2ErrHi);
  gr2->GetXaxis()->SetMoreLogLabels();
  gr2->GetXaxis()->SetNoExponent();
  gr2->Print("range");

  TGraphAsymmErrors* gr3=getAsymGraph(vsEt, etBinSet,etaBinSet,iBin,*sf,sf3ErrLo,sf3ErrHi);
  gr3->GetXaxis()->SetMoreLogLabels();
  gr3->GetXaxis()->SetNoExponent();
  gr3->Print("range");


  double *loc_etBinLimits=DYTools::getEtBinLimits(etBinSet);
  double *loc_etaBinLimits=DYTools::getEtaBinLimits(etaBinSet);
  int signedEta=DYTools::signedEtaBinning(etaBinSet);
  TString binStrForTitle=(vsEt) ? TString(Form(" %5.3lf #leq %s #leq %5.3lf",loc_etaBinLimits[iBin],(signedEta)?"#eta":"abs(#eta)",loc_etaBinLimits[iBin+1])) :
    TString(Form(" %2.0lf #leq #it{E}_{T} #leq %2.0lf GeV",loc_etBinLimits[iBin],loc_etBinLimits[iBin+1]));
  TString cpTitle=dataKind+ binStrForTitle;
  TString xAxisTitle="#it{E}_{T} [GeV]";
  if (!vsEt) xAxisTitle=(signedEta) ? "#eta" : "|#eta|";

  ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,"comp",cpTitle,
		      xAxisTitle,effKind + TString(" scale factor"),"ratio");
  cp.SetRefIdx(-111); // no ratio plot
  if (vsEt) cp.SetLogx();
  cp.AddLine(10.,1.,500.,1.,kBlack,2);

  TCanvas *cx=new TCanvas("cx","cx",600,600);
  SetSideSpaces(cx,0.05,0.,0.,0.02);

  gr1->GetYaxis()->SetTitleOffset(1.4);

  //cp.AddGraph(gr3,label3,"LPE",43); //kRed+3);
  cp.AddGraph(gr3,label3,"LPE",kRed); //kRed+3);
  cp.AddGraph(gr2,label2,"LPE",38); //kBlue+2);
  cp.AddGraph(gr1,label1," PE",kBlack,20); //24

  cp.Draw(cx,0,"png",0);

  cp.TransLegend(transLegendX, transLegendY);
  cp.WidenLegend(0.25,0.);

  cx->Update();

  // Save file
  if (saveFileTag.Length()) {
    TString outfname=TString("fig-sf-egamma-") + cpTitle;
    outfname.ReplaceAll(" #leq "," ");
    outfname.ReplaceAll(" ","_");
    outfname.ReplaceAll("(#eta)","Eta");
    outfname.ReplaceAll("#eta","eta");
    outfname.ReplaceAll(".","_");
    outfname.ReplaceAll("#it{E}_{T}","Et");
    //fname.Append(".png");
    std::cout << "outfname=" << outfname << "\n";

    TString locOutDir=TString("plots") + saveFileTag;
    if (doSave) {
      locOutDir.ReplaceAll("--","");
      SaveCanvas(cx,outfname,locOutDir);
    }
    else {
      std::cout << "... canvas not saved, as requested\n";
      std::cout << "   locOutDir=" << locOutDir << "\n";
    }
  }

  return ;
}
