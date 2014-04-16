#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/EventWeight.hh"

void chkEventWeight(int the_case=0) {
  EventWeight_t ev;
  TString rndStudyStr;

  if (the_case==0) {
    std::cout << "default case\n";
  }
  else if (the_case==999) {
    std::cout << "probe default case\n";
    rndStudyStr="PU_RND_0";
  }
  else {
    rndStudyStr=Form("PU_RND_%d",the_case);
    std::cout << "test <" << rndStudyStr << ">\n";
  }

  if (!ev.init(1,0,DYTools::NO_SYST,rndStudyStr)) {
    std::cout << " failure\n";
    return;
  }

  for (int nPV=10; nPV<25; nPV++) {
    ev.setPUWeight(nPV);
    printf("nPV=%2d w=%7.4lf\n",nPV,ev.puWeight());
  }
}
