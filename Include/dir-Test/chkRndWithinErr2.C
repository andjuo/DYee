#include "../DYTools.hh"
#include "../MyTools.hh"
#include <TRandom3.h>
#include "../HistoPair.hh"
#include "../CrossSection.hh"

void chkRndWithinErr2(int nExps) {
  HistoPair2D_t hSrc("hSrc");
  for (int ir=1; ir<=hSrc.getNrows(); ir++) {
    for (int ic=1; ic<=hSrc.getNcols(); ic++) {
      hSrc.setBinContent(ir,ic, ir*10+0.1*ic);
      hSrc.setBinError(ir,ic, 0.2*ir + 0.1*ic);
      hSrc.setBinSystError(ir,ic, 0.3*ir + 0.2*ic);
    }
  }

  printHisto(hSrc);

  HistoPair2D_t hRnd("hRnd");
  TH2D *hSumStat=createBaseH2("hSumStat");
  TH2D *hSumSyst=createBaseH2("hSumSyst");

  MXSectD_t MSumStat("mSumStat",hRnd.getNrows(),hRnd.getNcols());

  for (int iexp=0; iexp<nExps; ++iexp) {
    hRnd.randomizeWithinErr(hSrc,0);
    if (nExps==1) printHisto(hRnd);
    accumulateForRndStudies(hSumStat,hRnd.histo());

    for (int ir=1; ir<=hRnd.getNrows(); ++ir) {
      for (int ic=1; ic<=hRnd.getNcols(); ++ic) {
	MSumStat.AccumulateForAvg(ir-1,ic-1,hRnd.getBinContent(ir,ic),1/double(nExps));
      }
    }

    if (nExps==1) printHisto(hSumStat);
    hRnd.randomizeWithinErr(hSrc,1);
    accumulateForRndStudies(hSumSyst,hRnd.histo());
  }
  accumulateForRndStudies_finalize(hSumStat,nExps);
  accumulateForRndStudies_finalize(hSumSyst,nExps);

  MSumStat.GetXS().Print();

  MXSectD_t MresStat("MresStat",MSumStat);
  MSumStat.CalcAvgAndRMS(MresStat);
  std::cout << "MresStat error should be close to the stat error of hSrc\n";
  std::cout << "Is it the same as hSumStat?\n";
  MresStat.GetXS().Print();
  MresStat.GetError().Print();

  std::cout << "hSumStat error should be close to the stat error of hSrc\n";
  printHisto(hSumStat);
  std::cout << "hSumSyst error should be close to the syst error of hSrc\n";
  printHisto(hSumSyst);
}
