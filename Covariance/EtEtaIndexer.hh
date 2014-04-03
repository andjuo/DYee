#ifndef EtEtaIndexer_HH
#define EtEtaIndexer_HH

#include <iostream>

// ---------------------------------

class EtEtaIndexer_t //: public BaseClass_t 
{
 protected:
  int FEtBinCount,FEtaBinCount;
  int FEtIdx,FEtaIdx;
 public:
  EtEtaIndexer_t(int set_etBinCount, int set_etaBinCount) :
    //BaseClass_t("EtEtaIndexer_t"), 
    FEtBinCount(set_etBinCount), FEtaBinCount(set_etaBinCount),
    FEtIdx(-1), FEtaIdx(-1)
    {}

  EtEtaIndexer_t(const EtEtaIndexer_t &a) :
    //BaseClass("EtEtaIndexer_t"),
    FEtBinCount(a.FEtBinCount), FEtaBinCount(a.FEtaBinCount),
    FEtIdx(a.FEtIdx), FEtaIdx(a.FEtaIdx)
    {}
    
  int maxFlatIdx() const { return flatEtEtaIdx(FEtBinCount,FEtaBinCount-1); }
  int isValid() const { return ((FEtIdx>=0) && (FEtIdx<FEtBinCount) && (FEtaIdx>=0) && (FEtaIdx<FEtaBinCount)); }

    // bin indexing starts from 0
  int flatEtEtaIdx(int etBin, int etaBin) const {
    return (etBin+etaBin*FEtBinCount);
  }

  int flatEtEtaIdx() const { return (FEtIdx + FEtaIdx*FEtBinCount); }
  int getIdx() const { return flatEtEtaIdx(); }
  int getEtBin() const { return FEtIdx; }
  int getEtaBin() const { return FEtaIdx; }

  int getEtBin(int flatIdx) const { return flatIdx%FEtBinCount; }
  int getEtaBin(int flatIdx) const { return flatIdx/FEtBinCount; }

  void setEtBin(int ibin) { FEtIdx=ibin; }
  void setEtaBin(int ibin) { FEtaIdx=ibin; }
  void setEtEta(int etBin, int etaBin) { FEtIdx=etBin; FEtaIdx=etaBin; }
  void setEtEtaBin(int flatIdx) { FEtIdx=getEtBin(flatIdx); FEtaIdx=getEtaBin(flatIdx); }

  /*
#ifdef DYTools_HH
  int setEtEta(double et, double eta) {
    FEtIdx=DYTools::findEtBin(et, etBinning);
    FEtaIdx=DYTools::findEtaBin(eta, etaBinning);
    return this->isValid();
  }
#endif
  */

  int Start() { FEtIdx=0; FEtaIdx=0; return 0; }

  int Next() { 
    FEtaIdx++;
    if (FEtaIdx >= FEtaBinCount) {
      FEtaIdx=0;
      FEtIdx++;
    }
    if (FEtIdx >= FEtBinCount) return -1;
    return this->flatEtEtaIdx();
  }

  EtEtaIndexer_t& operator=(int flatIdx) {
    FEtIdx=this->getEtBin(flatIdx);
    FEtaIdx=this->getEtaBin(flatIdx);
    return *this;
  }

  friend
  inline int operator==(const EtEtaIndexer_t &a, const EtEtaIndexer_t &b) {
    return ((a.FEtBinCount==b.FEtBinCount) &&
	    (a.FEtaBinCount==b.FEtaBinCount) &&
	    (a.FEtIdx==b.FEtIdx) &&
	    (a.FEtaIdx==b.FEtaIdx)) ? 1:0;
  }

  friend
  inline int operator!=(const EtEtaIndexer_t &a, const EtEtaIndexer_t &b) {
    return (a==b) ? 0:1;
  }

  friend
  inline std::ostream& operator<< (std::ostream &out, const EtEtaIndexer_t &a) {
    out << "(" << a.FEtIdx << "," << a.FEtaIdx << ")";
    return out;
  }

  void PrintEtEtaVals(const double *loc_etBinLimits, const double *loc_etaBinLimits, std::ostream &out=std::cout) {
    out << "(" << loc_etBinLimits[FEtIdx] << "," << loc_etaBinLimits[FEtaIdx] << ")";
  }

};


#endif
