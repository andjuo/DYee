// adapted from EventScaleFactors/compareESF.C

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
  return h2;
}

// ------------------------------------------------------------


void compareESF(TString esfFile1="scale_factors_etaMax24_asymHLT_1D.root",
		TString esfFile2="covRhoFile_nMB41_etaMax24_asymHLT_100.root",
		int nDim=1) {
  TString path1="../root_files_reg/constants/DY_j22_19712pb/";
  TString path2="./";

  if (((nDim==1) &&  DYTools::study2D) ||
      ((nDim==2) && !DYTools::study2D)) {
    std::cout << "the macro uses basic implementations that require nDim info to match study2D in DYTools.hh\n";
    return;
  }

  TString dimStr;
  if (nDim==1) dimStr.Append("1D");
  else if (nDim==2) dimStr.Append("2D");
  else {
    std::cout << "error: nDim=" << nDim << "\n";
    return;
  }


  TString label1="calcEventEff";
  TString label2="studyEffCov";

  /*
  if (DYTools::study2D) {
    esfLongStr1.ReplaceAll("1D","2D");
    esfLongStr2.ReplaceAll("1D","2D");
  }
  */

  TString fname1=path1 + esfFile1;
  TString fname2=path2 + esfFile2;

  TH2D *h2esf1= loadESF(fname1,label1);
  TH2D *h2esf2= loadESF(fname2,label2);

  TH2D *h2diff= (TH2D*)h2esf1->Clone("diff");
  h2diff->Add(h2esf2,-1);
  h2diff->SetTitle("diff");

  gStyle->SetPalette(1);
  TCanvas *cx=new TCanvas("cx","cx",1200,400);
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

  if (DYTools::study2D) {
    std::vector<TH1D*> hProfEsf1,hProfEsf2;
    for (int mbin=2; mbin<=7; mbin++) {
      TString hname1=Form("hProf1_m%d",mbin);
      TString hname2=Form("hProf2_m%d",mbin);
      int set_nYBins=(mbin==7) ? 12:24;
      TH1D *h1= createProfileY(h2esf1,mbin,hname1,1,hname1, set_nYBins,0.,2.401);
      TH1D *h2= createProfileY(h2esf2,mbin,hname2,1,hname2, set_nYBins,0.,2.401);
      hProfEsf1.push_back(h1);
      hProfEsf2.push_back(h2);
    }

    std::vector<ComparisonPlot_t*> cpV;
    for (int i=0; i<6; ++i) {
      TString mRange=Form("M_%2.0lf_%2.0lf",DYTools::massBinLimits[i+1],DYTools::massBinLimits[i+2]);
      ComparisonPlot_t *cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,Form("cp_%s",mRange.Data()),mRange,"|y|","event scale factor","ratio");
      cp->AddHist1D(hProfEsf1[i],label1,"LPE1",kBlack,1,0);
      cp->AddHist1D(hProfEsf2[i],label2,"LPE1",kBlue ,2,0);
      cpV.push_back(cp);
    }

    TCanvas *cy=new TCanvas("cy","cy",1200,800);
    cpV[0]->Prepare6Pads(cy, 1);
    
    for (int i=0; i<6; ++i) {
      cpV[i]->Draw6(cy,1,i+1);
    }
    cy->Update();
  }
  else {
    TH1D *h1= createProfileX(h2esf1,1,"hProf1",1,"hProf1");
    TH1D *h2= createProfileX(h2esf2,1,"hProf2",1,"hProf2");

    ComparisonPlot_t *cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,"cpESF","","M_{ee} [GeV]", "event scale factor","ratio");
    cp->SetRatioYRangeC(1,0.002);
    cp->SetRatioYRange(0.96, 1.00);
    cp->SetLogx();
    cp->AddHist1D(h1, label1,"LPE1",kBlack,1,0);
    cp->AddHist1D(h2, label2,"LPE1",kBlue ,2,0);

    TCanvas *cy=new TCanvas("cy","cy",600,700);
    cp->Prepare2Pads(cy);
    cp->Draw(cy);
    cp->TransLegend(-0.1,-0.6);
    cy->Update();
  }

  //h2diff->Print("range");
  h2esf1->Print("range");

}
