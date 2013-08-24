#include <TROOT.h>
#include <TH1F.h>

#include "../Include/MyTools.hh"

void tryIt() {
  TH1F *h1=new TH1F("h1","h1",2,0,1);
  TH1F *h2=new TH1F("h2","h2",2,0,1);
  h1->SetBinContent(1,0.5);
  h1->SetBinError(1,0.5);
  h2->SetBinContent(1,0.5);
  h2->SetBinError(1,0.5);
  h2->Add(h1,2);
  h1->Print("range");
  h2->Print("range");


  if (1) {
    TH2D *h2d=new TH2D("h2d","",2,0.,2.,4,10.,20.);
    TH2D *h2dErr=(TH2D*)h2d->Clone("h2dErr");

    h2d->SetBinContent(1,1, 0.1);
    h2d->SetBinError(1,1,0.5);
    h2dErr->SetBinContent(1,2,100.);
    h2dErr->SetBinError(1,2,0.2);
    printHistoErr(h2d,h2dErr);
  }
}
