#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/ComparisonPlot.hh"
#include "../Include/ElectronEnergyScale.hh"
#include <TRandom.h>

void plotRndEScale(int theCase) {
  ElectronEnergyScale esBase(ElectronEnergyScale::Date20140220_2012_j22_peak_position);
  TString escaleStr=esBase.calibrationSetName();

  int seedMin=1000;
  int seedMax=1100;
  int dSeed=1;

  gRandom->SetSeed(seedMin);

  std::vector<TH1D*> hV;

  TH1D* hEScale=esBase.createScaleHisto("hEScaleBase");
  HERE("xx");
  printHisto(hEScale);

  hV.reserve(esBase.numberOfEtaBins());
  for (int ieta=0; ieta<esBase.numberOfEtaBins(); ++ieta) {
    const TAxis* ax=hEScale->GetXaxis();
    TString hName=Form("hEScale_eta_%2.1lf_%2.1lf",
		       esBase.getEtaLowEdge(ieta),
		       esBase.getEtaHighEdge(ieta));
    TH1D *h=new TH1D(hName,hName,200,0.98,1.02);
    h->SetDirectory(0);
    std::cout << " - histo " << hName << "\n";
    hV.push_back(h);
  }

  for (int iSeed=seedMin; iSeed<=seedMax; iSeed+=dSeed) {
    TString rndEScale=escaleStr + TString(Form("_RANDOMIZED%d",iSeed));
    ElectronEnergyScale es(rndEScale);
    es.print();
    for (int ieta=0; ieta<es.numberOfEtaBins(); ++ieta) {
      double eta=es.getEtaBinCenter(ieta+1);
      double factor=es.getEnergyScaleCorrectionRandomized(eta);
      if (theCase==1) factor/=0.99841;
      std::cout << " ieta=" << ieta << ", eta=" << eta << ", factor=" << factor << "\n";
      hV[ieta]->Fill(factor);
    }
  }

  for (unsigned int i=0; i<hV.size(); ++i) {
    TString canvName=hV[i]->GetName();
    canvName.ReplaceAll("hEScaleEta","canvEScale_");
    TCanvas *cx=new TCanvas(canvName,canvName,600,600);
    hV[i]->Draw();
    cx->Update();
  }
  return;
}
