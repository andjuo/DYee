#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/EventWeight.hh"
#include "../Include/AccessOrigNtuples.hh"

// -------------------------------------------------------------

class AccessOrigNtuplesDebug_t : public AccessOrigNtuples_t {
public:
  AccessOrigNtuplesDebug_t() : AccessOrigNtuples_t() {
    genBrIsActive=1;
  }

  mithep::TGenInfo* editGenPtr() { return gen; }
};


// -------------------------------------------------------------

class EventWeightDebug_t : public EventWeight_t {
public:
  EventWeightDebug_t() : EventWeight_t() {}

  TFSRSystematics_t *getFSR() { return fFSR; }

  double nextRndValue() {
    if (!fFSR) return 0.;
    return fFSR->nextValue();
  }

  UInt_t maxEntries() const { return (fFSR) ? fFSR->maxEntries() : 0; }

};

// -----------------------------------

void chkEventWeight(int pu_study=1, int the_case=0) {
  EventWeightDebug_t ev;
  TString rndStudyStr;

  if (pu_study==1) {
    std::cout << "Pile-up reweight study\n";

    if (the_case==0) {
      std::cout << "default case\n";
    }
    else if (the_case==999) {
      std::cout << "probe default case\n";
      rndStudyStr="PU_RND_STUDY0";
    }
    else {
      rndStudyStr=Form("PU_RND_STUDY%d",the_case);
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

  else {
    std::cout << "FSR study\n";
    rndStudyStr=Form("FSR_RND_STUDY%d",the_case);
    std::cout << "test <" << rndStudyStr << ">\n";

    if (!ev.init(1,0,DYTools::NO_SYST,rndStudyStr)) {
      std::cout << " failure\n";
      return;
    }

    ev.PrintDetails();

    if (0) {
      std::cout << " next=" << ev.nextRndValue() << "\n";
      ev.PrintDetails();
    }

    AccessOrigNtuplesDebug_t accessInfo;
    accessInfo.editGenPtr()->vmass = 20;
    accessInfo.editGenPtr()->mass  = 10;

    if (0) {
      for (int i=0; i<5; ++i) {
	ev.setSpecWeightValue(accessInfo,1,1.11111);
	ev.PrintDetails();
      }
    }
    else {
      TH1D *h= new TH1D("hRnd",Form("file %s",rndStudyStr.Data()),
			400,0.,2.);
      h->SetDirectory(0);
      for (UInt_t i=0; i<ev.maxEntries(); ++i) {
	ev.setSpecWeightValue(accessInfo,1,1.11111);
	h->Fill(ev.specWeight());
      }

      TCanvas *cx=new TCanvas("cx","cx",800,800);
      h->SetMarkerStyle(24);
      h->Draw("LPE");
      cx->Update();
    }
  }
}
