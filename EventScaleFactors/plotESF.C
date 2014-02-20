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


void plotESF(int doSave=0) {

  TString fnameBase="scale_factors_";

  TString dimStr;
  if (DYTools::study2D==0) dimStr.Append("1D");
  else dimStr.Append("2D");

  TString inpFName;
  TString fnameTag, saveDirTag="-plot-esf";
  TString label1;

  if (1) { 
    inpFName="/home/andriusj/cms/DYee-20140217-EESF/root_files_reg/constants/DY_j22_19712pb/scale_factors_ourSF_asymHLT_1D.root";
    label1="our SF";
    fnameTag="-ourSF";
  }

  if (1) {
    inpFName="/home/andriusj/cms/DYee-20140217-EESF/root_files_reg/constants/DY_j22_19712pb_egamma/scale_factors_egamma_asymHLT_1D.root"; 
    label1="EGamma";
    fnameTag="-egamma";
  }

  if (0) {
    inpFName="../Covariance/covRhoFile_nMB41_egamma_asymHLT_100.root";
    label1="EGamma";
    fnameTag="-egammaV2";
  }

  if (DYTools::study2D) {
    inpFName.ReplaceAll("1D","2D");
  }

  TH2D *h2esf1= loadESF(inpFName,label1);
  //h2esf1->Print("range");

  gStyle->SetPalette(1);
  TCanvas *cx=NULL;

  if (0) {
    cx=new TCanvas("cx","cx",500,500);
    AdjustFor2DplotWithHeight(cx);
    
    cx->SetLogx();
    h2esf1->Draw("COLZ");
    cx->Update();
  }

  TCanvas *cy=NULL;

  if (DYTools::study2D) {
    h2esf1->Print("range");

    std::vector<TH1D*> hProfEsf1;
    for (int mbin=2; mbin<=7; mbin++) {
      TString hname1=Form("hProf1_m%d",mbin);
      TString hname2=Form("hProf2_m%d",mbin);
      int set_nYBins=DYTools::nYBins[mbin-1]; //(mbin==7) ? 12:24;
      TH1D *h1= createProfileY(h2esf1,mbin,hname1,1,hname1, set_nYBins,0.,DYTools::yRangeMax+1e-4);

      hProfEsf1.push_back(h1);
    }

    std::vector<ComparisonPlot_t*> cpV;
    for (int i=0; i<6; ++i) {
      TString mRange=Form("M_%2.0lf_%2.0lf",DYTools::massBinLimits[i+1],DYTools::massBinLimits[i+2]);
      ComparisonPlot_t *cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,Form("cp_%s",mRange.Data()),mRange,"|y|","#it{#rho}_{Data/MC}","ratio");
      cp->SetRefIdx(-111); // no ratio
      cp->SetRatioNdivisions(404);
      cp->SetYAxisTextSizes(0.08,1.35, 0.06);
      cp->SetXAxisTextSizes(0.07,0.83, 0.06);
      cp->SetRatioYTitleSize(0.18,0.40);

      //cp->SetRatioYRange(0.97,1.01);
      //cp->SetRatioYRange(set_ratio_y_min,set_ratio_y_max);
      cp->AddHist1D(hProfEsf1[i],label1,"LPE1",kBlack,1,0,  -1);

      cp->AddLine(0.,1.0, 2.4,1.0, kBlue+1,kDashed);
      cpV.push_back(cp);
    }

    cy=new TCanvas("cy","cy",1200,700);
    cy->Divide(3,2);
    SetSideSpaces(cy,0.08,-0.04,-0.02,0.0);
    
    for (int i=0; i<6; ++i) {
      //cpV[i]->Draw6(cy,1,i+1);
      cpV[i]->Draw(cy,false,"png",i+1);
      cpV[i]->TransLegend(-0.42,-0.1);
      cpV[i]->WidenLegend(0.1,0.1);
      cpV[i]->WidenLegend(0.2,0.0);
      //cpV[i]->
    }
    cy->Update();
  }
  else {
    TH1D *h1= createProfileX(h2esf1,1,"hProf1",1,"hProf1");
    h1->GetXaxis()->SetNoExponent();
    h1->GetXaxis()->SetMoreLogLabels();

    ComparisonPlot_t *cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,"cpESF","","M_{ee} [GeV]", "#it{#rho}_{Data/MC}","ratio");
    cp->SetRefIdx(-111); // no ratio
    cp->SetYTitleSize(0.07,1.35);
    cp->SetRatioYTitleSize(0.17,0.47);
    cp->SetPrintRatio(1);

    //cp->SetRatioYRangeC(1,0.04);
    //cp->SetRatioYRange(set_ratio_y_min, set_ratio_y_max);
    cp->SetLogx();
    //h1->GetYaxis()->SetTitleOffset(1.47);
    cp->AddHist1D(h1, label1,"LPE1",kBlack,1,0,  -1);

    //cp->AddLine(15,1.0, 1500.,1.0, kBlue+1,kDashed);

    cy=new TCanvas("cy","cy",600,580);
    SetSideSpaces(cy,0.05,0.02,-0.00,0.02);


   //cp->Prepare2Pads(cy);
    cp->Draw(cy,false,"png",0);
    //cp->TransLegend(-0.2,-0.6);
    //cp->WidenLegend(0.2,0.);
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


}
