#include "../DYTools.hh"
#include "../MyTools.hh"

void chkExtractHisto() {
  int dimX=10;
  int dimY=5;

  TH2D *hSrc=new TH2D("hSrc","",dimX,0.,100., dimY,0.,500.);
  hSrc->SetDirectory(0);
  hSrc->GetXaxis()->SetTitle("myX");
  hSrc->GetYaxis()->SetTitle("myY");

  for (int ibin=1; ibin<=hSrc->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=hSrc->GetNbinsY(); ++jbin) {
      hSrc->SetBinContent(ibin,jbin, ibin + 0.01*jbin);
      hSrc->SetBinError  (ibin,jbin, 0.01*(ibin + 0.01*jbin));
    }
  }

  if (0) { // full range
    if (1) {
      TH2D *h2=extractSubArea(hSrc, 1,10,1,5, "h2",1);
      h2->Print("range");
    }
    
    if (1) {
      TH2D *h2reset=extractSubArea(hSrc, 1,10,1,5, "h2_reset",1,1);
      h2reset->Print("range");
    }
  }

  if (1) { // partial area
    if (1) {
      TH2D *h2=extractSubArea(hSrc, 2,5,3,4, "h2",1);
      h2->Print("range");
    }
    
    if (1) {
      TH2D *h2reset=extractSubArea(hSrc, 2,5,3,4, "h2_reset",1,1);
      h2reset->Print("range");
    }
  }

  if (0) {
    TH1D *h1=createProfileX(hSrc, 2, "hProfile2x",1);
    h1->Print("range");
  }

  if (0) {
    TH1D *h1=createProfileY(hSrc, 2, "hProfile2y",1);
    h1->Print("range");
  }
  return;
}

