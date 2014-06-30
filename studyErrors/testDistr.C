#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/ComparisonPlot.hh"
#include <TMath.h>

void testDistr(int the_case=0, double par1=0., double par2=0.) {
  TH1D *h1=new TH1D("h1","h1",1000,0.,5.);
  h1->SetDirectory(0);
  h1->SetStats(true);
  gStyle->SetOptStat("mr");
  //  h1->SetOptStat("mr"); // To display the mean and RMS:   SetOptStat("mr");
 

  TString label;
  switch(the_case) {
  case 0: label="log-normal"; break;
  case 1: label="gamma"; break;
  default:
    std::cout << "not ready\n";
    return;
  }

  double x,y;
  for (int ibin=1; ibin<=h1->GetNbinsX(); ++ibin) {
    x= h1->GetBinCenter(ibin);
    switch(the_case) {
    case 0: y= TMath::LogNormal(x,par1,0,par2); break;
    case 1: y= TMath::GammaDist(x,par1,0,par2); break;
    default:
      std::cout << "fnc not ready\n";
      return;
    }
    h1->SetBinContent(ibin,y);
  }

  std::cout << " mean=" << h1->GetMean() << " +- " << h1->GetRMS() << "\n";

  TCanvas *cx=new TCanvas("cx",label,600,600);
  h1->Draw();
  cx->Update();
}
