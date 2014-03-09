#include "../DYTools.hh"
#include "../MyTools.hh"

void chkTH2DRelDiff() {
  TH2D *h2=new TH2D("h2","h2",5,0.5,5.5,2,0.5,2.5);
  TH2D *h2p1=(TH2D*)h2->Clone("h2p1");
  TH2D *h2m2=(TH2D*)h2->Clone("h2m2");

  for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
    for (int jbin=1; jbin<=h2->GetNbinsY(); jbin++) {
      h2->SetBinContent(ibin,jbin, ibin);
      h2p1->SetBinContent(ibin,jbin, ibin+1+jbin*0.1);
      h2m2->SetBinContent(ibin,jbin, ibin-2+jbin+0.2);
    }
  }

  std::cout << h2->GetName() << "\n"; h2->Print("range");
  std::cout << h2p1->GetName() << "\n"; h2p1->Print("range");
  std::cout << h2m2->GetName() << "\n"; h2m2->Print("range");

  TH2D *h2NoDiff=getRelDifference(h2,"nodiff",h2);
  printHisto(h2NoDiff);
  
  TH2D *h2Diff1=getRelDifference(h2,"diff1",h2m2);
  printHisto(h2Diff1);
  
  TH2D *h2Diff2=getRelDifference(h2,"diff2",h2p1,h2m2);
  printHisto(h2Diff2);
  
}
