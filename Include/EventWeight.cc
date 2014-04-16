#include "../Include/EventWeight.hh"


// --------------------------------------------------------

EventWeight_t::EventWeight_t(const EventWeight_t &ev) :
  fDoPUReweight(ev.fDoPUReweight),
  fFewzCorr(ev.fFewzCorr),
  fPUReweight(ev.fPUReweight),
  fFEWZ(ev.fFEWZ),
  fPtrOwner(0), // uses ptrs of another object
  baseW(ev.baseW),
  puW(ev.puW),
  fewzW(ev.fewzW),
  specW(ev.specW)
{}

// --------------------------------------------------------

EventWeight_t::~EventWeight_t() {
  if (fPtrOwner && fPUReweight) delete fPUReweight;
  if (fPtrOwner && fFEWZ) delete fFEWZ;
}

// --------------------------------------------------------

int EventWeight_t::init(int do_puReweight, int do_fewzCorr,
	         DYTools::TSystematicsStudy_t systMode, TString rndStudyStr) {
  fDoPUReweight=0;
  fFewzCorr=0;

  if (fPtrOwner && fPUReweight) { delete fPUReweight; }
  if (fPtrOwner && fFEWZ) { delete fFEWZ; }
  fPtrOwner=1;
  fPUReweight=NULL;
  fFEWZ=NULL;
  
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
  if (ok && do_puReweight) {

    int regularPileup=1;
    int idxPU=-1;
    if (rndStudyStr.Length()!=0) {
      int idx=rndStudyStr.Index("PU_RND_");
      if (idx>=0) {
	idxPU=atoi(rndStudyStr.Data() + idx + strlen("PU_RND_"));
	std::cout << "PU_RND detected. idx=" << idxPU << "\n";
	regularPileup=0;
      }
    }

    if (regularPileup) {
      fPUReweight=new PUReweight_t(systMode);
    }
    else {
#ifndef DYee8TeV_reg
      std::cout << "non-regular pileup reweight is ready for DYee_8TeV_reg\n";
      ok=0;
#else
      fPUReweight=new PUReweight_t(PUReweight_t::_none);
      if (fPUReweight) {
	TString targetFile="../root_files_reg/pileup/8TeV_reg/randomized_pileup_20140415.root";
	TString targetField=Form("hRnd_lumibased_data_%d",idxPU);
	if (idxPU==0) targetField="pileup_lumibased_data_base";
	else if (idxPU== 111) targetField="pileup_lumibased_data_111";
	else if (idxPU==-111) targetField="pileup_lumibased_data_-111";
	TString sourceFile="../root_files_reg/pileup/8TeV_reg/mcPileupHildreth_mean_full2012_20131106_repacked.root";
	ok=fPUReweight->setSimpleWeights(targetFile,targetField,
					 sourceFile,"pileup_simulevel_mc");
	if (!ok) {
	  std::cout << "error setting simple weights\n";
	}
      }
#endif
    }

    if (ok && !fPUReweight) ok=0; else fDoPUReweight=do_puReweight;

    // hard-coded check
#ifndef DYee8TeV_reg
    double expectWeight_at_11=(DYTools::energy8TeV) ? 1.283627 : 1.32768;
    if (fabs(fPUReweight->getWeight(11) - expectWeight_at_11)>1e-4) {
      std::cout << "EventWeight::init failed hard-coded check\n";
      assert(0);
    }
#endif
  }

  if (ok && do_fewzCorr) {
    bool cutZPt100= ((do_fewzCorr && 2)!=0) ? true : false;
    fFEWZ = new FEWZ_t(true, cutZPt100);
    if (!fFEWZ || !fFEWZ->isInitialized()) ok=0;
    else fFewzCorr=do_fewzCorr;
  }
  return ok;
}


// --------------------------------------------------------

// --------------------------------------------------------


