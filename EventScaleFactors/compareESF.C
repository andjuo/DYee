#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
//#include "../Include/MitStyleRemix.hh"
#include "../Include/colorPalettes.hh"

// ------------------------------------------------------------

TH2D* loadESF(TString fname, TString label) {
  TH2D* h2=LoadMatrixFields(fname, 1, "scaleFactor", "scaleFactorErr", 1, 1);
  h2->SetTitle(label);
  return h2;
}

// ------------------------------------------------------------


void compareESF(TString esfLongStr1="1D_Full2012_hltEffOld_PU",
		TString esfLongStr2="1D",
		int nDim=1) {
  TString path1="/home/andriusj/cms/DYee8TeV-20140118/root_files/constants/DY_j22_19712pb/";
  TString path2="/home/andriusj/cms/DYee-20131024/root_files_reg/constants/DY_j22_19712pb/";
  //path2="/home/andriusj/cms/CMSSW_3_8_4/src/DYee8TeV-20130801/DrellYanDMDY/root_files/tag_and_probe/DY_j22_19789pb/";
  //path2="/home/andriusj/cms/DYee8TeV-20140118-maxEta24/root_files/tag_and_probe/DY_j22_19712pb/";


  TString fnameBase="scale_factors_";

  TString dimStr;
  if (nDim==1) dimStr.Append("1D");
  else if (nDim==2) dimStr.Append("2D");
  else {
    std::cout << "error: nDim=" << nDim << "\n";
    return;
  }

  TString fname1=path1 + fnameBase + esfLongStr1 + TString(".root");
  TString fname2=path2 + fnameBase + esfLongStr2 + TString(".root");


  TString label1="DrellYanDMDY";
  TString label2="DYee";

  TH2D *h2esf1= loadESF(fname1,label1);
  TH2D *h2esf2= loadESF(fname2,label2);

  TH2D *h2diff= (TH2D*)h2esf1->Clone("diff");
  h2diff->Add(h2esf2,-1);
  h2diff->SetTitle("diff");

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

}
