#include "../DYTools.hh"
#include "../MyTools.hh"

void chkAddHistos(int chk_study2D) {
  if (!DYTools::setup(chk_study2D)) {
    return;
  }
 
  std::vector<TH1D*> h1V;
  std::vector<TString> labelsV;
  labelsV.push_back("1st");
  labelsV.push_back("2nd");

  createBaseH1Vec(h1V,"hChk_",labelsV,"counts",1);
  for (unsigned int ih=0; ih<h1V.size(); ++ih) {
    TH1D* h=h1V[ih];
    double factor=(ih==0) ? 0.01:100;
    for (int imb=1; imb<=h->GetNbinsX(); ++imb) {
      h->SetBinContent(imb, factor*imb);
    }
  }
  TH1D* hSum=addHistos("hSum",h1V);
  printHisto(hSum);
  return;
}
