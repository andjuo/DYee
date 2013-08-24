#include "../Include/EventWeight.hh"


// --------------------------------------------------------

EventWeight_t::~EventWeight_t() {
  if (fPUReweight) delete fPUReweight;
  if (fFEWZ) delete fFEWZ;
}

// --------------------------------------------------------

int EventWeight_t::init(int do_puReweight, int do_fewzCorr) {
  fDoPUReweight=0;
  fFewzCorr=0;

  if (fPUReweight) { delete fPUReweight; fPUReweight=NULL; }
  if (fFEWZ) { delete fFEWZ; fFEWZ=NULL; }
  
  int ok=1;
  if (do_puReweight) {
    fPUReweight=new PUReweight_t();
    if (!fPUReweight) ok=0; else fDoPUReweight=do_puReweight;

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


