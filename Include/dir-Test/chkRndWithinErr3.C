#include "../DYTools.hh"
#include "../MyTools.hh"
#include <TRandom3.h>
#include "../HistoPair.hh"

void chkRndWithinErr3(int nExps) {
  HistoPair2D_t hSrc("hSrc");
  for (int ir=1; ir<=hSrc.getNrows(); ir++) {
    for (int ic=1; ic<=hSrc.getNcols(); ic++) {
      hSrc.setBinContent(ir,ic, ir*10+0.1*ic);
      hSrc.setBinError(ir,ic, 0.2*ir + 0.1*ic);
      hSrc.setBinSystError(ir,ic, 0.3*ir + 0.2*ic);
    }
  }

  printHisto(hSrc);

  gRandom->SetSeed(-1);

  HistoPair2D_t hRnd("hRnd");
  TH2D *hSumStat=createBaseH2("hSumStat");
  //TH2D *hSumSyst=createBaseH2("hSumSyst");

  std::vector<HistoPair2D_t*> hRndVec;
  hRndVec.reserve(nExps);

  for (int iexp=0; iexp<nExps; ++iexp) {
    hRnd.randomizeWithinErr(hSrc,0);
    if (nExps==1) printHisto(hRnd);
    accumulateForRndStudies(hSumStat,hRnd.histo());
    if (nExps==1) printHisto(hSumStat);

    TString tmpName=Form("hpRnd_%d",iexp);
    hRndVec.push_back(new HistoPair2D_t(tmpName,hRnd,1));

    //hRnd.randomizeWithinErr(hSrc,1);
    //accumulateForRndStudies(hSumSyst,hRnd.histo());
  }

  int unbiasedEstimate=0;
  accumulateForRndStudies_finalize(hSumStat,nExps,unbiasedEstimate);
  //accumulateForRndStudies_finalize(hSumSyst,nExps);

  std::cout << "hSumStat error should be close to the stat error of hSrc\n";
  printHisto(hSumStat);

  // derive the covariance matrix
  HistoPair2D_t hAvg("hAvg");
  TMatrixD *covM=deriveCovMFromRndStudies(hRndVec,unbiasedEstimate,&hAvg);
  if (!covM) {
    std::cout << "derivation of covM failed\n";
    return ;
  }
  else {
    std::cout << "hAvg central value should be the same as that of hSumStat\n";
    std::cout << "hAvg error should be close to that of the stat error of hSrc\n";
    printHisto(hAvg);
  }

  if (1 && (hRndVec.size()==2)) {
    int cnt=2;
    for (int i=0; i<cnt; ++i) {
      for (int j=0; j<cnt; ++j) {
	for (unsigned int ii=0; ii<hRndVec.size(); ++ii) {
	  std::cout << Form("vec[%u](%d,%d)= ",ii,i,j) << hRndVec[ii]->GetBinContent(i+1,j+1) << "\n";
	}
      }
    }
    for (int i=0; i<cnt; ++i) {
      for (int j=0; j<cnt; ++j) {
	std::cout << Form("cov(%d,%d)=",i,j) << " " << (*covM)(i,j) << "\n";
      }
    }
  }

  //std::cout << "hSumSyst error should be close to the syst error of hSrc\n";
  //printHisto(hSumSyst);
}
