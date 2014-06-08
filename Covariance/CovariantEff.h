#ifndef CovariantEff_H
#define CovariantEff_H

//#include "../EventScaleFactors/calcEventEff.C"
#include "../EventScaleFactors/calcEventEffLink.h"
#include "CovariantMatrix.hh"

//extern EffArray_t *ro_Data_loc;
//extern EffArray_t *ro_MC_loc;

// ---------------------------------

class EtEtaIndexer_t;

// ---------------------------------
// ---------------------------------
// ---------------------------------

class CovariantEffMgr_t : public BaseClass_t {
 protected:
  InputFileMgr_t FInpMgr;
  EventSelector_t FEvtSelector;
  TString FPUStr;
  TriggerSelection_t FTriggers;
  int FInitOk;
  std::vector<TMatrixD*> FRhoRelSystErrs;  // what was provided by user
  std::vector<TMatrixD*> FRhoExtraFactors;  // pseudo-exps
 public:
  CovariantEffMgr_t(); 

  const InputFileMgr_t& mgr() const { return FInpMgr; }
  InputFileMgr_t& editMgr() { return FInpMgr; }

  int initOk() const { return FInitOk; }
  int Setup(const TString &confFileName, int nExps, DYTools::TSystematicsStudy_t systMode);
  int SetupSFsyst(const TString &confFileName, const TString &recoSystFName, const TString &idSystFName, const TString &hltSystFName, int nExps, DYTools::TSystematicsStudy_t systMode, int egammaSystOnly); 

  template<class idx_t>
  double getExtraSF(idx_t iExp, const EtEtaIndexer_t &fidx1, const EtEtaIndexer_t &fidx2) const;
};

// ---------------------------------

class EtEtaIndexer_t //: public BaseClass_t
{
 protected:
  int FEtBinCount,FEtaBinCount;
  int FEtIdx,FEtaIdx;
  int FIdx;
 public:
  EtEtaIndexer_t(int set_etBinCount, int set_etaBinCount) :
    //BaseClass_t("EtEtaIndexer_t"), 
    FEtBinCount(set_etBinCount), FEtaBinCount(set_etaBinCount),
    FEtIdx(-1), FEtaIdx(-1), FIdx(-1)
    {}

  EtEtaIndexer_t(const EtEtaIndexer_t &a) :
    //BaseClass("EtEtaIndexer_t"),
    FEtBinCount(a.FEtBinCount), FEtaBinCount(a.FEtaBinCount),
    FEtIdx(a.FEtIdx), FEtaIdx(a.FEtaIdx), FIdx(a.FIdx)
    {}
    
  int maxFlatIdx() const { return flatEtEtaIdx(FEtBinCount,FEtaBinCount-1); }
  int isValid() const { return ((FEtIdx>=0) && (FEtIdx<FEtBinCount) && (FEtaIdx>=0) && (FEtaIdx<FEtaBinCount)); }


  // Serve as converter
  // bin indexing starts from 0
  int flatEtEtaIdx(int etBin, int etaBin) const {
    return (etBin+etaBin*FEtBinCount);
  }
  int getEtBin(int flatIdx) const { return flatIdx%FEtBinCount; }
  int getEtaBin(int flatIdx) const { return flatIdx/FEtBinCount; }

  // Access to member data
  int getIdx() const { return FIdx; }
  int getEtBin() const { return FEtIdx; }
  int getEtaBin() const { return FEtaIdx; }

  void setEtBin(int ibin) { FEtIdx=ibin; }
  void setEtaBin(int ibin) { FEtaIdx=ibin; }

  void setEtEta(int etBin, int etaBin) {
    FEtIdx=etBin;
    FEtaIdx=etaBin;
    FIdx=flatEtEtaIdx(FEtIdx,FEtaIdx);
  }

  void setEtEtaBin(int flatIdx) {
    FEtIdx=getEtBin(flatIdx);
    FEtaIdx=getEtaBin(flatIdx);
    FIdx=flatEtEtaIdx(FEtIdx,FEtaIdx);
  }

  int setEtEta(double et, double eta) {
    FEtIdx=DYTools::findEtBin(et, etBinning);
    FEtaIdx=DYTools::findEtaBin(eta, etaBinning);
    FIdx=flatEtEtaIdx(FEtIdx,FEtaIdx);
    return this->isValid();
  }

  int Start() { FEtIdx=0; FEtaIdx=0; FIdx=0; return 0; }

  int Next() { 
    FEtIdx++;
    FIdx++;
    if (FEtIdx >= FEtBinCount) {
      FEtIdx=0;
      FEtaIdx++;
    }
    if (FEtaIdx >= FEtaBinCount) return -1;
    return FIdx;
  }

  friend
  inline int operator==(const EtEtaIndexer_t &a, const EtEtaIndexer_t &b) {
    return ((a.FEtBinCount==b.FEtBinCount) &&
	    (a.FEtaBinCount==b.FEtaBinCount) &&
	    (a.FEtIdx==b.FEtIdx) &&
	    (a.FEtaIdx==b.FEtaIdx) &&
	    (a.FIdx==b.FIdx)) ? 1:0;
  }

  friend
  inline int operator!=(const EtEtaIndexer_t &a, const EtEtaIndexer_t &b) {
    return (a==b) ? 0:1;
  }

  friend
  inline std::ostream& operator<< (std::ostream &out, const EtEtaIndexer_t &a) {
    out << "(" << a.FIdx << " : " << a.FEtIdx << "," << a.FEtaIdx << ")";
    return out;
  }

  void PrintEtEtaVals(const double *loc_etBinLimits, const double *loc_etaBinLimits, std::ostream &out=std::cout) {
    out << "(" << loc_etBinLimits[FEtIdx] << "," << loc_etaBinLimits[FEtaIdx] << ")";
  }

};

// ---------------------------------
// split the Et bin containing 17GeV

class EtEtaIndexer17_t //: public BaseClass_t
{
 protected:
  int FEtBinCount,FEtaBinCount;
  int FEtIdx,FEtaIdx;
  int FIdx;
 public:
  EtEtaIndexer17_t(int set_etBinCount, int set_etaBinCount) :
    //BaseClass_t("EtEtaIndexer_t"),
    FEtBinCount(set_etBinCount+1), FEtaBinCount(set_etaBinCount),
    FEtIdx(-1), FEtaIdx(-1), FIdx(-1)
    {}

  EtEtaIndexer17_t(const EtEtaIndexer17_t &a) :
    //BaseClass("EtEtaIndexer_t"),
    FEtBinCount(a.FEtBinCount), FEtaBinCount(a.FEtaBinCount),
    FEtIdx(a.FEtIdx), FEtaIdx(a.FEtaIdx), FIdx(a.FIdx)
    {}

  int maxFlatIdx() const { return flatEtEtaIdx(FEtBinCount,FEtaBinCount-1); }
  int isValid() const { return ((FEtIdx>=0) && (FEtIdx<FEtBinCount) && (FEtaIdx>=0) && (FEtaIdx<FEtaBinCount)); }


  // Serve as converter
  // bin indexing starts from 0
  int flatEtEtaIdx(int etBin, int etaBin) const {
    return (etBin+etaBin*FEtBinCount);
  }
  int getEtBin(int flatIdx) const { return flatIdx%FEtBinCount; }
  int getEtaBin(int flatIdx) const { return flatIdx/FEtBinCount; }

  // Access to member data
  int getIdx() const { return FIdx; }
  int getEtBin() const { return FEtIdx; }
  int getEtaBin() const { return FEtaIdx; }

  void setEtBin(int ibin) { FEtIdx=ibin; }
  void setEtaBin(int ibin) { FEtaIdx=ibin; }

  void setEtEta(int etBin, int etaBin) {
    FEtIdx=etBin;
    FEtaIdx=etaBin;
    FIdx=flatEtEtaIdx(FEtIdx,FEtaIdx);
  }

  void setEtEtaBin(int flatIdx) {
    FEtIdx=getEtBin(flatIdx);
    FEtaIdx=getEtaBin(flatIdx);
    FIdx=flatEtEtaIdx(FEtIdx,FEtaIdx);
  }

  int setEtEta(double et, double eta) {
    FEtIdx=DYTools::findEtBin(et, etBinning);
    if (et>17.) { if (et>499.) FEtIdx=FEtBinCount-1; else FEtIdx++; }
    FEtaIdx=DYTools::findEtaBin(eta, etaBinning);
    FIdx=flatEtEtaIdx(FEtIdx,FEtaIdx);
    return this->isValid();
  }

  int Start() { FEtIdx=0; FEtaIdx=0; FIdx=0; return 0; }

  int Next() {
    FEtIdx++;
    FIdx++;
    if (FEtIdx >= FEtBinCount) {
      FEtIdx=0;
      FEtaIdx++;
    }
    if (FEtaIdx >= FEtaBinCount) return -1;
    return FIdx;
  }

  friend
    inline int operator==(const EtEtaIndexer17_t &a, const EtEtaIndexer17_t &b) {
    return ((a.FEtBinCount==b.FEtBinCount) &&
	    (a.FEtaBinCount==b.FEtaBinCount) &&
	    (a.FEtIdx==b.FEtIdx) &&
	    (a.FEtaIdx==b.FEtaIdx) &&
	    (a.FIdx==b.FIdx)) ? 1:0;
  }

  friend
  inline int operator!=(const EtEtaIndexer17_t &a, const EtEtaIndexer17_t &b) {
    return (a==b) ? 0:1;
  }

  friend
  inline std::ostream& operator<< (std::ostream &out, const EtEtaIndexer17_t &a) {
    out << "(" << a.FIdx << " : " << a.FEtIdx << "," << a.FEtaIdx << ")";
    return out;
  }

  void PrintEtEtaVals(const double *loc_etBinLimits, const double *loc_etaBinLimits, std::ostream &out=std::cout) {
    if ((loc_etBinLimits[FEtIdx+1]<17.) || (loc_etBinLimits[FEtIdx]>17.)) {
      int useEtIdx=FEtIdx;
      if (loc_etBinLimits[FEtIdx]>17.) useEtIdx--;
      out << "(" << loc_etBinLimits[useEtIdx] << "-"
	  << loc_etBinLimits[useEtIdx+1] << ",";
    }
    else {
      // split bin
      if (FEtIdx==DYTools::findEtBin(16.99, etBinning)) {
	out << "(" << loc_etBinLimits[FEtIdx] << "-17,";
      }
      else {
	out << "(17-" << loc_etBinLimits[FEtIdx+1] << ",";
      }
    }
    out << loc_etaBinLimits[FEtaIdx] << "-"
	<< loc_etaBinLimits[FEtaIdx+1] << ")";
  }
};

// ---------------------------------
// ---------------------------------

template<class idx_t>
double CovariantEffMgr_t::getExtraSF(idx_t iExp, const EtEtaIndexer_t &fidx1, const EtEtaIndexer_t &fidx2) const {
  if (FRhoExtraFactors.size()==0) return 1.;
  if ((unsigned int)(iExp) >= FRhoExtraFactors.size()) {
    reportError("getExtraSF(%d): there are %d experiments prepared",(int)(iExp),(int)(FRhoExtraFactors.size()));
    return 0.;
  }
  double factor1= (*FRhoExtraFactors[iExp])(fidx1.getEtBin(),fidx1.getEtaBin());
  double factor2= (fidx2.isValid()) ? (*FRhoExtraFactors[iExp])(fidx2.getEtBin(),fidx2.getEtaBin()) : 1.;
  return factor1*factor2;
}

// ---------------------------------
// ---------------------------------

inline
void AddToCovMElement(CovarianceMatrix_t &M, const EtEtaIndexer_t &idx1, const EtEtaIndexer_t &idx2, double val) {
  M.AddValue(idx1.getIdx(),idx2.getIdx(), val);
}

// ---------------------------------
// ---------------------------------

class StudyArr_t : public BaseClass_t {
  std::vector<std::vector<TH1D*> *> FVec;
  TString FBaseName;
  int FReady;
 public:
  StudyArr_t(const TString &baseName, int nMassBins, int nExps,
	     int binCount=0, double minVal=0., double maxVal=0.);

  int Init(int nMassBins, int nExps,
	   int binCount, double minVal, double maxVal);

  TH1D* operator()(int i, int j) { return (*FVec[i])[j]; }
};

//const int cSystArrayBinCount=DYTools::nMassBins;
//typedef TH1D** SystTH1DArray_t[cSystArrayBinCount];   // mass index

//int PrepareStudyArray(const TString &name_base, int nExps, SystTH1DArray_t &arr, int binCount=150, double minVal=0.0, double maxVal=1.5);

// ---------------------------------
// ---------------------------------

class ScaleFactor_t : public BaseClass_t {
 protected:
  EtEtaIndexer17_t fTmpIdx;
  int fKind;
  TMatrixD fSF;
  std::vector<TMatrixD*> fSFrndV;
 public:
  ScaleFactor_t(int set_kind,
		int use_etBinCount, int use_etaBinCount);
  ~ScaleFactor_t() { ClearVec(fSFrndV); }

  int kind() const { return fKind; }
  const TMatrixD& scaleFactorM() const { return fSF; }
  const std::vector<TMatrixD*>& rndScaleFactorMvec() const { return fSFrndV; }
  const TMatrixD* rndScaleFactorM(int i) const { return fSFrndV[i]; }

  // no optimization. Only interface to calcEventEff.C
  double eventScaleFactor(const esfSelectEvent_t &selData) const;
  double rndEventScaleFactor(int idx, const esfSelectEvent_t &selData) const;

  // optimization
  double eventScaleFactor(const EtEtaIndexer17_t &fidx1,
			  const EtEtaIndexer17_t &fidx2) const {
    return fSF(fidx1.getIdx(),fidx2.getIdx());
  }

  double rndEventScaleFactor(int iexp,
			     const EtEtaIndexer17_t &fidx1,
			     const EtEtaIndexer17_t &fidx2) const {
    return (*fSFrndV[iexp])(fidx1.getIdx(),fidx2.getIdx());
  }

  int init(int set_kind,
	   DYTools::TEtBinSet_t set_etBinning,
	   DYTools::TEtaBinSet_t set_etaBinning,
	   int nExps);

};

// ---------------------------------
// ---------------------------------
/*

double GetEffVal(DYTools::TDataKind_t data_or_mc, DYTools::TEfficiencyKind_t kind, int scEtBin, int scEtaBin) {
  const std::vector<TMatrixD*> *eff= (data_or_mc==DYTools::DATA) ?  &dataEff : &mcEff;
  double val=(*(*eff)[int(kind)])[scEtBin][scEtaBin];
  return val;
}

// ---------------------------------

double GetEffVal(DYTools::TDataKind_t data_or_mc, DYTools::TEfficiencyKind_t kind, const EtEtaIndexer_t &fidx) {
  return GetEffVal(data_or_mc,kind,fidx.getEtBin(),fidx.getEtaBin());
}

// ---------------------------------

int GetEffValAndErr(DYTools::TDataKind_t data_or_mc, DYTools::TEfficiencyKind_t kind, int scEtBin, int scEtaBin, double &val, double &avgErr) {
  const std::vector<TMatrixD*> *eff= (data_or_mc==DYTools::DATA) ?  &dataEff : &mcEff;
  const std::vector<TMatrixD*> *effErrLo= (data_or_mc==DYTools::DATA) ?  &dataEffErrLo : &mcEffErrLo;
  const std::vector<TMatrixD*> *effErrHi= (data_or_mc==DYTools::DATA) ?  &dataEffErrHi : &mcEffErrHi;

  val=(*(*eff)[int(kind)])[scEtBin][scEtaBin];
  avgErr=0.5*( (*(*effErrLo)[int(kind)])[scEtBin][scEtaBin] +
	       (*(*effErrHi)[int(kind)])[scEtBin][scEtaBin]);
  return 1;
}

// ---------------------------------

double GetEffAvgErr(DYTools::TDataKind_t data_or_mc, DYTools::TEfficiencyKind_t kind, int scEtBin, int scEtaBin) {
  const std::vector<TMatrixD*> *effAvgErr= (data_or_mc==DYTools::DATA) ?  &dataEffAvgErr : &mcEffAvgErr;

  double avgErr=(*(*effAvgErr)[int(kind)])[scEtBin][scEtaBin];
  return avgErr;
}

// ---------------------------------

double GetEffAvgErr(DYTools::TDataKind_t data_or_mc, DYTools::TEfficiencyKind_t kind, EtEtaIndexer_t fidx) {
  return GetEffAvgErr(data_or_mc,kind,fidx.getEtBin(),fidx.getEtaBin());
}

// ---------------------------------

double GetEffValSmeared(DYTools::TDataKind_t data_or_mc, DYTools::TEfficiencyKind_t kind, int scEtBin, int scEtaBin, int iexp) {
  //HERE("GetEffValSmeared: iexp=%d",iexp);
  const std::vector<TMatrixD*> *eff= (data_or_mc==DYTools::DATA) ?  &dataEff : &mcEff;
  const std::vector<TMatrixD*> *effErrAvg= (data_or_mc==DYTools::DATA) ?  &dataEffAvgErr : &mcEffAvgErr;
  //const EffArray_t *rnd= (data_or_mc==DYTools::DATA) ? &(ro_Data_loc[iexp]) : &(ro_MC_loc[iexp]);
  const EffArray_t *rnd= (data_or_mc==DYTools::DATA) ? &(ro_Data[iexp]) : &(ro_MC[iexp]);

  double rndSigma= (*rnd)[kind][scEtBin][scEtaBin];

  //HERE("GetEffValSmeared");
  double val=(*(*eff)[int(kind)])[scEtBin][scEtaBin]
    + rndSigma * (*(*effErrAvg)[int(kind)])[scEtBin][scEtaBin];
  return val;
}

// ---------------------------------

double GetEffValSmeared(DYTools::TDataKind_t data_or_mc, DYTools::TEfficiencyKind_t kind, const EtEtaIndexer_t &idx, int iexp) {
  return GetEffValSmeared(data_or_mc,kind,idx.getEtBin(),idx.getEtaBin(),iexp);
}
 
// ---------------------------------

int GetEffValAndErrSmeared(DYTools::TDataKind_t data_or_mc, DYTools::TEfficiencyKind_t kind, int scEtBin, int scEtaBin, int iexp, double &val, double &valErr) {
  val=GetEffValSmeared(data_or_mc,kind,scEtBin,scEtaBin,iexp);
  valErr=0.;
  return 1;
}

// ---------------------------------


int FillEffConstantValue(TH1D *histo, DYTools::TDataKind_t data_or_mc, DYTools::TEfficiencyKind_t kind, int scEtBin, int scEtaBin, int iexp=-1, int print=0) {
  double val=0.,valErr=0.;
  int res=(iexp<0) ?
    GetEffValAndErr(data_or_mc,kind,scEtBin,scEtaBin, val,valErr) :
    GetEffValAndErrSmeared(data_or_mc,kind,scEtBin,scEtaBin,iexp, val,valErr);
  if (!res) return 0;
  if (print) std::cout << "FillEffConstantValue(" << histo->GetName() 
		       << ", " << data_or_mc
		       << ", " << kind 
		       << ", EtBin=" << scEtBin 
		       << ", EtaBin=" << scEtaBin
		       << ", iexp=" << iexp
		       << ") : " << val << "  " << valErr 
		       << std::endl;
  histo->Reset();
  for (int i=1; i<=histo->GetNbinsX(); ++i) {
    histo->SetBinContent(i, val);
    histo->SetBinError(i, valErr);
  }
  return 1;
}

// ---------------------------------
*/

#endif
