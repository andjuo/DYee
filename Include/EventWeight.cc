#include "../Include/EventWeight.hh"

#include <TRandom3.h>
#warning EventWeight.cc included TRandom3


// --------------------------------------------------------
// --------------------------------------------------------

int TFSRSystematics_t::init() {
  fReady=0;
  if (DetermineSystematicsStudy(fRndStudyStr)==DYTools::FSR_RND_STUDY) {
    // determine seed
    fSeed=0;
    // Special seed FSR_RND_STUDYid#fixed#. The value after
    // 'fixed' is the random number
    int specSeed= (fRndStudyStr.Index("fixed")!=-1) ? 1:0;
    // try to determine the seed
    for (int i=0; (fSeed==0) && (i<fRndStudyStr.Length()); i++) {
      fSeed=atoi(fRndStudyStr.Data() + i);
      //std::cout << "chk fSeed=" << fSeed << "\n";
    }
    if (fSeed==0) {
      std::cout << "TFSRSystematics::init - failed to determine seed\n";
      return 0;
    }
    // got seed
    if (!specSeed && (abs(fSeed)!=111)) {
      TString fname=Form("../root_files_reg/theory/rndSequence_%d.root",fSeed);
      fFile = new TFile(fname,"read");
      if (!fFile || !fFile->IsOpen()) {
	std::cout << "failed to open the file <" << fname << ">\n";
	return 0;
      }
      fTree= (TTree*)fFile->Get("rndSequence");
      if (!fTree) {
	std::cout << "failed to get the tree\n";
	return 0;
      }
      fTree->SetBranchAddress("rnd",&fRnd);
      fBr= fTree->GetBranch("rnd");
      if (!fBr) {
	std::cout << "failed to get the branch\n";
	return 0;
      }

      fAvailableEntries= fTree->GetEntries();
      fEntry=0;
    }
    else {
      fAvailableEntries=-1;
      fEntry=0;
      if (specSeed) {
	if ((fRndStudyStr.Index("fixed")!=-1) && (fRndStudyStr.Index("id")!=-1)) {
	  std::cout << "fRndStudyStr=<" << fRndStudyStr << ">\n";
	  int idx=fRndStudyStr.Index("fixed");
	  fRnd=atof(fRndStudyStr.Data() + idx + 5);
	}
	else {
	  gRandom->SetSeed(fSeed);
	  fRnd=gRandom->Gaus(0,1.);
	}
      }
      else fRnd=(fSeed<0) ? -1 : 1;
    }
    fReady=1;
  }
  return fReady;
}


// --------------------------------------------------------

void TFSRSystematics_t::printDetails() const {
  std::cout << "FSRSystematics(str=" << fRndStudyStr << "):\n";
  if (!fReady) std::cout << "  NOT READY\n";
  else {
    std::cout << " - seed       = " << fSeed << "\n";
    std::cout << " - entry(max) = " << fEntry << "("
	      << fAvailableEntries << ")\n";
    std::cout << " - current rnd= " << fRnd << "\n";
  }
}

// --------------------------------------------------------
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
  specW(ev.specW),
  fFSR(NULL)
{
  if (ev.fFSR) fFSR=new TFSRSystematics_t(*ev.fFSR);
}

// --------------------------------------------------------

EventWeight_t::~EventWeight_t() {
  if (fPtrOwner && fPUReweight) delete fPUReweight;
  if (fPtrOwner && fFEWZ) delete fFEWZ;
  if (fFSR) delete fFSR;
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
    if ((rndStudyStr.Length()!=0) &&
	(DetermineSystematicsStudy(rndStudyStr)==DYTools::PU_RND_STUDY)) {
      for (int i=0; (idxPU<=0) && (i<rndStudyStr.Length()); ++i) {
	idxPU=atoi(rndStudyStr.Data() + i);
      }
      if (idxPU>=0) {
	std::cout << "PU_RND_STUDY detected. idx=" << idxPU << "\n";
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
	TString targetFile="../root_files_reg/pileup/8TeV_reg/randomized_pileup_20140424.root";
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

  if (ok && rndStudyStr.Length()
      && (DetermineSystematicsStudy(rndStudyStr)==DYTools::FSR_RND_STUDY)) {
    fFSR= new TFSRSystematics_t(rndStudyStr);
    if (!fFSR) ok=0;
    if (ok && !fFSR->ready()) ok=0;
  }

  return ok;
}


// --------------------------------------------------------

// --------------------------------------------------------


