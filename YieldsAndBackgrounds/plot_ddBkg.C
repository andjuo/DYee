#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
//#include "../Include/MitStyleRemix.hh"
#include "../Include/colorPalettes.hh"
#include "../Include/ComparisonPlot.hh"

// ------------------------------------------------------------

TH2D* loadESF(TString fname, TString label, TString fieldName, TString fieldErrName) {
  int checkBinning=0;
  TH2D* h2=LoadMatrixFields(fname, checkBinning, fieldName, fieldErrName, 1, 1);
  h2->SetTitle(label);
  h2->GetXaxis()->SetMoreLogLabels();
  h2->GetXaxis()->SetNoExponent();
  std::cout << "loading from <" << fname << ">\n";
  std::cout << "got "; h2->Print("range");
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


void plot_ddBkg(int iBr=0,
		int nDim=1,
		int doSave=0) {

  TString esfLongStr1, esfLongStr2;
  TString esfLongStr3,esfLongStr4;
  TString path1,path2;
  TString path3,path4;

  if (((nDim==1) &&  DYTools::study2D) ||
      ((nDim==2) && !DYTools::study2D)) {
    std::cout << "the macro uses basic implementations that require nDim info to match study2D in DYTools.hh\n";
    return;
  }


  TString fnameBase1="fakeBkgDataPoints";
  TString fnameBase2=fnameBase1;

  TString fnameBase3, fnameBase4;
  TString loadFieldName="fakeBackgroundFromData";

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
  double set_ratio_y_min=0.96;
  double set_ratio_y_max=1.04;

  int compSet=-1;
  int swapColors=0;

  TMatrixD *setYRanges2D=NULL;

  if (1) { // added 2014.03.14
    path1="../root_files_reg/ddbkgYield/DY_j22_19712pb/";
    path2=path1;
    esfLongStr1="_20131231_1D";
    esfLongStr2="_20140312_1D";
    label1="2013.12.31";
    label2="2014.03.12";
    fnameTag="-cmpDDBkg";
    saveDirTag="-cmpDDBkg";
    if (nDim==2) { set_ratio_y_min=0.9; set_ratio_y_max=1.15; }
    else { set_ratio_y_min=0.95; set_ratio_y_max=1.15; }
    compSet=-10;
    swapColors=0;
    //setYRanges2D=enforceYRanges(1);
  }


  if (DYTools::study2D) {
    esfLongStr1.ReplaceAll("1D","2D");
    esfLongStr2.ReplaceAll("1D","2D");
    esfLongStr3.ReplaceAll("1D","2D");
    esfLongStr4.ReplaceAll("1D","2D");
  }

  if (iBr==1) {
    fnameBase1="true2eBkgDataPoints";
    loadFieldName="true2eBackgroundFromData";
    fnameBase2=fnameBase1;
  }

  TString fname1=path1 + fnameBase1 + esfLongStr1 + TString(".root");
  TString fname2=path2 + fnameBase2 + esfLongStr2 + TString(".root");

  TString loadFieldErrName=loadFieldName + TString("Error");
  TH2D *h2esf1= loadESF(fname1,label1, loadFieldName,loadFieldErrName);
  h2esf1->Print("range");
  TH2D *h2esf2= loadESF(fname2,label2, loadFieldName,loadFieldErrName);

  TH2D *h2esf3=NULL, *h2esf4=NULL;

  if (fnameBase3.Length()) {
    TString fname3= path3 + fnameBase3 + esfLongStr3 + TString(".root");
    h2esf3= loadESF(fname3,label3, loadFieldName,loadFieldErrName);
  }
  if (fnameBase4.Length()) {
    TString fname4= path4 + fnameBase4 + esfLongStr4 + TString(".root");
    h2esf4= loadESF(fname4,label4, loadFieldName,loadFieldErrName);
  }

  TH2D *h2diff= (TH2D*)h2esf1->Clone("diff");
  h2diff->Add(h2esf2,-1);
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
      cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,Form("cp_%s",mRange.Data()),mRange,"|y|",loadFieldName,"ratio");
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
    TH1D *h3= (h2esf3) ? createProfileX(h2esf3,1,"hProf3",1,"hProf3") : NULL;
    TH1D *h4= (h2esf4) ? createProfileX(h2esf4,1,"hProf4",1,"hProf4") : NULL;

    //std::cout << "\n\ncreated profile from " << h2esf2->GetName() << "\n";
    //h2esf2->Print("range");
    //h2->Print("range");

    int relRatio=0;
    ComparisonPlot_t *cp=NULL;
    if (!relRatio) cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,"cpESF","","M_{ee} [GeV]", loadFieldName,"ratio");
    else cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioRel,"cpESF","","M_{ee} [GeV]", loadFieldName,"rel.ratio");
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
