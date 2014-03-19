#ifndef EventWeight_HH
#define EventWeight_HH

// ----------------------------------------------------

#include "../Include/DYTools.hh"
#include "../Include/PUReweight.hh"
#include "../Include/FEWZ.hh"
#include "../Include/TGenInfo.hh"
#include "../Include/AccessOrigNtuples.hh"


// The class EventWeigh_t has four fields:
//  baseW -- determined from input of the file
//  puW, fewzW -- correction factors are set automatically
//  specW -- user can control this factor (e.g. FSR reweight)

// ----------------------------------------------------
// ----------------------------------------------------


class EventWeight_t  {
protected:
  int fDoPUReweight;
  int fFewzCorr;
  
  PUReweight_t *fPUReweight;
  FEWZ_t *fFEWZ;

  double baseW;
  double puW;
  double fewzW;
  double specW;
public:
  EventWeight_t() :
    fDoPUReweight(1),
    fFewzCorr(0),
    fPUReweight(NULL),
    fFEWZ(NULL),
    baseW(1.),
    puW(1.),
    fewzW(1.),
    specW(1.)
  {}

  ~EventWeight_t();

  int init(int do_puReweight, int do_fewzCorr, DYTools::TSystematicsStudy_t systMode);

  void setFewzWeight(const mithep::TGenInfo* gen) {
    fewzW=(gen && fFEWZ) ? fFEWZ->getWeight(gen->vmass,gen->vpt,gen->vy) : 1.0;
  }

  void setPUWeight(int nPV) {
    puW= (fPUReweight) ? fPUReweight->getWeight(nPV) : 1.0;
    //std::cout << "setPUWeight(nPV=" << nPV << ") = " << puW << "\n";
  }

  void setPUWeightValue(double val) { puW=val; }
  void setSpecWeightValue(double val) { specW=val; }

  void setSpecWeightValue(const AccessOrigNtuples_t &accessInfo, double massDiff, double specWeight_val) {
    const mithep::TGenInfo* gen=accessInfo.genPtr();
    specW = (gen->vmass - gen->mass <= massDiff) ? 1.0 : specWeight_val;
  }

  void setFewzWeight(const AccessOrigNtuples_t &accessInfo) {
    this->setFewzWeight(accessInfo.genPtr());
  }

  int set_PU_and_FEWZ_weights(const AccessOrigNtuples_t &accessInfo, bool isData) {
    const mithep::TGenInfo* gen=accessInfo.genPtr();
    this->setFewzWeight(gen);
    if (fPUReweight) {
      int nPVs=accessInfo.getNPV(isData);
      this->setPUWeight(nPVs);
    }
    else puW=1.0;
    return 1;
  }


  double baseWeight() const { return baseW; }
  double puWeight() const { return puW; }
  double fewzWeight() const { return fewzW; }
  double specWeight() const { return specW; }
  double totalWeight() const { return baseW * puW * fewzW * specW; }

  void setBaseWeight(double w) { baseW=w; }
  void setBaseWeight(const EventWeight_t &w) { baseW=w.baseW; }

  template<class TLong_t>
  int setWeight_and_adjustMaxEvents(TLong_t &maxEvents, double lumi, double xsec, double extraFactor, int doWeight) {
    // Determine maximum number of events to consider
    // *** CASES ***
    // <> lumi < 0                             => use all events in the sample
    // <> xsec = 0                             => for data (use all events)
    // <> lumi > 0, xsec > 0, doWeight = true  => use all events and scale to lumi
    // <> lumi > 0, xsec > 0, doWeight = false => compute expected number of events
    //std::cout << "setWeight_and_adjustMaxEvents(maxEvents=" << maxEvents << ", lumi=" << lumi << ", xsec=" << xsec << ", extraFactor=" << extraFactor << ", doWeight=" << doWeight << ")\n";

    baseW=1.0;
    if ((lumi>0) && (xsec!=0.)) {
      if (doWeight) { 
	baseW = lumi * xsec/double(maxEvents); 
	std::cout << "setting baseW=" << lumi << " * " << xsec << " / " << maxEvents << " = " << baseW << " (1/baseW=" << 1/baseW << ")\n";
      }
      else { 
	TLong_t needsEvents=(TLong_t)(lumi*xsec); 
	if (needsEvents>maxEvents) {
	  std::cout << "Not enough events for " << lumi << " pb^-1 in file\n";
	  return 0;
	}
	maxEvents=needsEvents;
      }
    }
    baseW*=extraFactor;
    std::cout << " weight corrected by extraFactor=" << extraFactor << ". final=" << baseW << "\n";
    return 1;
  }

  void PrintDetails() const {
    std::cout << " baseW=" << baseW << ", specW=" << specW << ", flag_puReweight=" << fDoPUReweight << ", flag_FewzReweight=" << fFewzCorr << "\n";
  }


  void Print(int short_output=1) const {
    if (short_output) {
      std::cout << "w(" << baseW << "b," << fewzW << "fewz)";
    }
    else std::cout << *this << "\n";
  }

  friend
  std::ostream& operator<<(std::ostream &out, const EventWeight_t &ew) {
    out << "w(" << ew.baseW << "b," << ew.puW << "pu, " << ew.fewzW << "f, " << ew.specW << "s)=" << ew.totalWeight();
    return out;
  }

};





// ----------------------------------------------------



#endif
