#include "../DYTools.hh"
#include "../MyTools.hh"

void chkCov() {
  TH2D *hBase=new TH2D("hBase","hBase",2,0.5,2.5, 1, 0., 1.);
  TH2D *h2a= Clone(hBase,"h2a");
  TH2D *h2b= Clone(hBase,"h2b");
  TH2D *h2c= Clone(hBase,"h2c");

  h2a->SetBinContent(1,1, -0.07);
  h2a->SetBinError  (1,1,  0.18);
  h2a->SetBinContent(2,1,  0.17);
  h2a->SetBinError  (2,1,  0.1);
  
  h2b->SetBinContent(1,1,  0.12);
  h2b->SetBinError  (1,1,  0.01);
  h2b->SetBinContent(2,1,  0.07);
  h2b->SetBinError  (2,1,  0.);

  h2c->SetBinContent(1,1,  0.28);
  h2c->SetBinError  (1,1,  0.17);
  h2c->SetBinContent(2,1, -0.03);
  h2c->SetBinError  (2,1,  0.1);

  std::vector<TH2D*> rndV;
  rndV.reserve(3);
  rndV.push_back(h2a);
  rndV.push_back(h2b);
  rndV.push_back(h2c);

  TH2D* avg=Clone(hBase,"hAvg");
  TMatrixD covErrM(2,2);

  TMatrixD* cov= deriveCovMFromRndStudies(rndV,0,avg,&covErrM);
  TMatrixD* corr= corrFromCov(*cov);

  printHisto(avg);
  std::cout << "covariance (off diagonal should be -0.01167)\n"; 
  cov->Print();
  std::cout << "correlation (off diagonal should be -0.999)\n";
  corr->Print();
  std::cout << "covariance error\n";
  covErrM.Print();

  std::cout << "fractional error\n";
  TMatrixD frErr(covErrM);
  divideMatrix(frErr,*cov);
  frErr.Print();
}
