#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/ComparisonPlot.hh"

void compareAtEdge(int iBr, int iSet=0, double scaleHisto1=1.) {
  TString fileName;
  TString histo1Name, histo2Name;
  TString label1,label2;
  TString cpTitle;
  double xMin=30., xMax=70.;
  double yMin=0., yMax=1.e5;

  if ((iBr==0) || (iBr==1)) {
    fileName="postFsrCS_2D_mb5GeV-preFsrDet.root";
    histo1Name="mass_1D_1GeV_bins/hMass1D_Madgraph_Zll_10_50";
    histo2Name="mass_1D_1GeV_bins/hMass1D_Madgraph_Zll__50";
    label1="Zll_10_50";
    label2="Zll_50";
    cpTitle="preFsrDet acceptance";
  }

  if (iBr==1) {
    fileName.ReplaceAll("preFsrDet","postFsrDet");
    cpTitle.ReplaceAll("preFsrDet","postFsrDet");
  }

  if (iSet==1) {
    std::vector<TString*> tsv;
    tsv.push_back(&histo1Name);
    tsv.push_back(&histo2Name);
    replaceAll(tsv,"mass_","massPreFsr_");
    replaceAll(tsv,"hMass","hMassPreFsr");
  }

  TFile fin(fileName,"read");
  if (!fin.IsOpen()) return;
  TH1D* h1=(TH1D*)fin.Get(histo1Name);
  h1->SetDirectory(0);
  TH1D* h2=(TH1D*)fin.Get(histo2Name);
  h2->SetDirectory(0);
  fin.Close();

  h1->Scale(scaleHisto1);
  TH1D* hSum=(TH1D*)h1->Clone("hSum");
  hSum->Add(h2);

  TH1D* hDer=(TH1D*)hSum->Clone("hDer");
  hDer->Reset();
  for (int ibin=1; ibin<hDer->GetNbinsX(); ++ibin) {
    hDer->SetBinContent(ibin,(hSum->GetBinContent(ibin+1)-hSum->GetBinContent(ibin-1))/2.);
  }

  TString yAxisTitle=(iSet==0) ? "post-FSR counts" : "pre-FSR counts";

  ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,"cp",cpTitle,
		      "M_{ee}",yAxisTitle,"ratio");
  cp.SetRefIdx(-111);
  cp.SetXRange(xMin,xMax);
  cp.SetYRange(yMin,yMax);
  cp.AddHist1D(h1,label1,"LP",kBlack);
  cp.AddHist1D(h2,label2,"LP",kBlue);
  cp.AddHist1D(hSum,"sum","LP",kRed+1);

  ComparisonPlot_t cpDer(ComparisonPlot_t::_ratioPlain,"cpDer","",
			 "M_{ee}","dSum/dx","ratio");
  cpDer.NoLegend(1);
  cpDer.SetRefIdx(-111);
  cpDer.SetXRange(xMin,xMax);
  cpDer.SetYRange(-2000,1000);
  cpDer.AddHist1D(hDer,"","LP",kRed+1);

  TCanvas *cx=new TCanvas("cx","cx",700,700);
  cp.Prepare2Pads(cx);
  cp.Draw(cx,false,"png",1);

  cpDer.Draw(cx,false,"png",2);

  cx->Update();
}
