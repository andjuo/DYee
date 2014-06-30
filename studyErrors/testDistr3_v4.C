#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/ComparisonPlot.hh"
#include "../Include/CrossSection.hh"
#include <TMath.h>
#include <TRandom3.h>
#include <TF1.h>

void testDistr3_v4(int the_case=0, double parMean=0., double par1a=0., double par1step=0.) {
  TH1D *h1a=new TH1D("h1","h1",100,0.,10.);
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

  std::vector<TF1*> fncV;

  h1a->SetLineColor(kBlack);
  h1b->SetLineColor(kBlue);
  h1c->SetLineColor(kGreen+1);

  VXSectD_t vx("vx",3);
  VXSectD_t vtest("vtest",3);
  vx.Zero();
  vtest.Zero();

  TString label;
  switch(the_case) {
  case -1: label="gaussian"; break;
  case 0: label="log-normal"; break;
  case 1: label="gamma"; break;
  case 2: label="polynomial"; break;
  default:
    std::cout << "not ready\n";
    return;
  }

  gRandom->SetSeed(10007);

  TString fncString;
  int nPars=3;
  double xRangeMin=0, xRangeMax=20.;
  double scale=0.5;
  switch(the_case) {
  case -1: fncString="gaus"; xRangeMin=-20;
    break;
  case 0: fncString= "[0]*TMath::LogNormal(x, [1], [2], [3])"; break;
  case 1: fncString= "[0]*TMath::GammaDist(x, [1], [2], [3])"; break;
  case 2: {
    double a=-1.4;
    double b= 13.2;
    double c=-7.8;
    double det=sqrt(b*b-4*a*c);
    if (det!=det) {
      std::cout << "negative det^2\n";
      return;
    }
    xRangeMin= (-b-det)/(2*a);
    xRangeMax= (-b+det)/(2*a);
    fncString=Form("%6.4lf*x*x+%6.4lf*x+%6.4lf",a,b,c);
  }
    break;
  default:
    std::cout << "fncNotReady\n";
    return;
  }

  for (int ipar=0; ipar<3; ++ipar) {
    double par1=par1a+par1step*ipar;
    
    TF1 *fnc= new TF1(Form("flog_%d",ipar),fncString, xRangeMin, xRangeMax);
    fnc->SetNpx(5000);
    if (nPars==3) fnc->SetParameters(scale,par1,0,parMean);
    else fnc->SetParameters(scale,parMean,par1);
    fncV.push_back(fnc);
  }

  double y=0, xtest=0;
  const int nExps=1000;
  for (int ipar=0; ipar<3; ++ipar) {
    for (int iexp=0; iexp<nExps; ++iexp) {
      //double par1=par1a+par1step*ipar;
      switch(the_case) {
      case -1:
      case 0:
      case 1:
      case 2:
	y= 1.0*fncV[ipar]->GetRandom();
	//y=xtest;
	break;
	//case 1: y= TMath::GammaDist(x,par1,0,parMean); break;
      default:
	std::cout << "fnc not ready (application)\n";
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

  //printHisto(hV[0]);

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
