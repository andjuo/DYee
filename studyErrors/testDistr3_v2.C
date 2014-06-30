#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/ComparisonPlot.hh"
#include "../Include/CrossSection.hh"
#include <TMath.h>
#include <TRandom3.h>

void testDistr3_v2(int the_case=0, double parMean=0., double par1a=0., double par1step=0.) {
  TH1D *h1a=new TH1D("h1","h1",20000,0.,10.);
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

  VXSectD_t vx("vx",3);
  VXSectD_t vtest("vtest",3);
  vx.Zero();
  vtest.Zero();

  TString label;
  switch(the_case) {
  case 0: label="log-normal"; break;
  case 1: label="gamma"; break;
  default:
    std::cout << "not ready\n";
    return;
  }

  gRandom->SetSeed(10007);

  double y, xtest;
  const int nExps=100000;
  for (int ipar=0; ipar<3; ++ipar) {
    for (int iexp=0; iexp<nExps; ++iexp) {
      double par1=par1a+par1step*ipar;
      switch(the_case) {
      case 0:
	xtest= parMean + par1 * gRandom->Gaus(double(0.),double(1.));
	y= exp(xtest);
	//y=xtest;
	break;
	//case 1: y= TMath::GammaDist(x,par1,0,parMean); break;
      default:
	std::cout << "fnc not ready\n";
	return;
      }
      //if (ipar==0) std::cout << "iexp=" << iexp << ", ipar=" << ipar << ", y=" << y << "\n";
      hV[ipar]->Fill(y);
      vx.Accumulate(ipar,y);
      vtest.Accumulate(ipar,xtest);
    }
  }
  vx.Accumulate_end(nExps);
  vtest.Accumulate_end(nExps);

  for (int ipar=0; ipar<3; ++ipar) {
    double par1=par1a+par1step*ipar;
    std::cout << "par=" << par1 << ",";
    std::cout << " mean=" << hV[ipar]->GetMean()
	      << " +- " << hV[ipar]->GetRMS() << "\n";
  }
  vx.Print();
  //vtest.Print();

  TCanvas *cx=new TCanvas("cx",label,600,600);
  h1a->Draw();
  h1b->Draw("L same");
  h1c->Draw("L same");
  cx->Update();
}
