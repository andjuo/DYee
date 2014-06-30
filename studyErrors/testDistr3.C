#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/ComparisonPlot.hh"
#include <TMath.h>

void testDistr3(int the_case=0, double parMean=0., double par1a=0., double par1step=0.) {
  TH1D *h1a=new TH1D("h1","h1",2000,0.,2.);
  h1a->SetDirectory(0);
  h1a->SetStats(true);
  gStyle->SetOptStat("mr");
  //  h1->SetOptStat("mr"); // To display the mean and RMS:   SetOptStat("mr");
  TH1D* h1b=Clone(h1a,"h1b","h1b");
  TH1D* h1c=Clone(h1a,"h1c","h1c");
  std::vector<TH1D*> hV;
  hV.push_back(h1a);
  hV.push_back(h1b);
  hV.push_back(h1c);

  h1a->SetLineColor(kBlack);
  h1b->SetLineColor(kBlue);
  h1c->SetLineColor(kGreen+1);


  TString label;
  switch(the_case) {
  case 0: label="log-normal"; break;
  case 1: label="gamma"; break;
  default:
    std::cout << "not ready\n";
    return;
  }

  double x,y;
  for (int ibin=1; ibin<=h1a->GetNbinsX(); ++ibin) {
    x= h1a->GetBinCenter(ibin);
    for (int ipar=0; ipar<3; ++ipar) {
      double par1=par1a+par1step*ipar;
      switch(the_case) {
      case 0: y= TMath::LogNormal(x,par1,0,parMean); break;
      case 1: y= TMath::GammaDist(x,par1,0,parMean); break;
      default:
	std::cout << "fnc not ready\n";
	return;
      }
      hV[ipar]->SetBinContent(ibin,y);
    }
  }

  for (int ipar=0; ipar<3; ++ipar) {
    double par1=par1a+par1step*ipar;
    std::cout << "par=" << par1 << ",";
    std::cout << " mean=" << hV[ipar]->GetMean()
	      << " +- " << hV[ipar]->GetRMS() << "\n";
  }

  TCanvas *cx=new TCanvas("cx",label,600,600);
  h1a->Draw();
  h1b->Draw("L same");
  h1c->Draw("L same");
  cx->Update();
}
