#include "../Include/EventWeight.hh"


// --------------------------------------------------------

EventWeight_t::~EventWeight_t() {
  if (fPUReweight) delete fPUReweight;
  if (fFEWZ) delete fFEWZ;
}

// --------------------------------------------------------

int EventWeight_t::init(int do_puReweight, int do_fewzCorr, DYTools::TSystematicsStudy_t systMode) {
  fDoPUReweight=0;
  fFewzCorr=0;

  if (fPUReweight) { delete fPUReweight; fPUReweight=NULL; }
  if (fFEWZ) { delete fFEWZ; fFEWZ=NULL; }
  
  if (do_puReweight &&
      ( (systMode==DYTools::NO_REWEIGHT) || (systMode==DYTools::NO_REWEIGHT_PU) )) {
	std::cout << "EventWeight::init: puReweight requested but the SystMode is " << systMode << "\n";
	do_puReweight=0;
  }
  if (do_fewzCorr &&
      ( (systMode==DYTools::NO_REWEIGHT) || (systMode==DYTools::NO_REWEIGHT_FEWZ) )) {
    std::cout << "EventWeight::init: FEWZ correction requested but the SystMode is " << systMode << "\n";
    do_fewzCorr=0;
  }

  int ok=1;
  if (do_puReweight) {

    fPUReweight=new PUReweight_t(systMode);
    if (!fPUReweight) ok=0; else fDoPUReweight=do_puReweight;

    // hard-coded check
#ifndef DYee8TeV_reg
    double expectWeight_at_11=(DYTools::energy8TeV) ? 1.283627 : 1.32768;
    if (fabs(fPUReweight->getWeight(11) - expectWeight_at_11)>1e-4) {
      std::cout << "EventWeight::init failed hard-coded check\n";
      assert(0);
    }
#endif
  }

  if (do_fewzCorr) {
    bool cutZPt100= ((do_fewzCorr && 2)!=0) ? true : false;
    fFEWZ = new FEWZ_t(true, cutZPt100);
    if (!fFEWZ || !fFEWZ->isInitialized()) ok=0;
    else fFewzCorr=do_fewzCorr;
  }
  return ok;
}


// --------------------------------------------------------

// --------------------------------------------------------


