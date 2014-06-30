#ifndef studyHelpers_HH
#define studyHelpers_HH

#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"

// -------------------------------------------------------------
// -------------------------------------------------------------

typedef enum { _valTBYieldGen=0, _valTBYieldTot, _valYieldTot, 
	       _valUnfYield, _valCS,
	       _valLast }
  TValueType_t;

typedef enum { _syval_sigYieldOrig=0, 
	       _syval_sigYieldGen, _syval_sigYieldActual,
	       _syval_unfYieldOrig, _syval_unfYield,
	       _syval_csVal,
	       _syvalLast }
  TSigYieldValueType_t;


// -------------------------------------------------------------

inline
int next(TValueType_t &it) {
  if (it!=_valLast) {
    it=TValueType_t(int(it)+1);
  }
  return (it!=_valLast) ? 1:0;
}

// -------------------------------------------------------------

TString valTypeName(TValueType_t valType);

inline
std::ostream& operator<<(std::ostream& out, TValueType_t valType) {
  out << valTypeName(valType);
  return out;
}

// -------------------------------------------------------------

inline
int next(TSigYieldValueType_t &it) {
  if (it!=_syvalLast) {
    it=TSigYieldValueType_t(int(it)+1);
  }
  return (it!=_syvalLast) ? 1:0;
}

// -------------------------------------------------------------

TString sigYieldValTypeName(TSigYieldValueType_t valType);

inline
std::ostream& operator<<(std::ostream& out, TSigYieldValueType_t valType) {
  out << sigYieldValTypeName(valType);
  return out;
}

// -------------------------------------------------------------
// -------------------------------------------------------------

int calculateSidedErrors(const TH1D* hVals, double centralValue,
			 double &errPos, double &errNeg, double &errRMS);

// -------------------------------------------------------------

void printSidedErrors(const TH1D* hVals, double centralValue=-9999.);

// -------------------------------------------------------------
// -------------------------------------------------------------

struct StudyInfo_t {
  std::vector<TH1D*> fRawHistos40,fRawHistos41;
  double fErr40pos[_valLast], fErr40neg[_valLast], fErr40[_valLast];
  double fErr41pos[_valLast], fErr41neg[_valLast], fErr41[_valLast];
  TH2D *fH2SigYield, *fH2UnfYield, *fH2TrueBkg, *fH2MainCS;
  TH2D *fH2AvgCS;
public:
  StudyInfo_t() :
    fRawHistos40(), fRawHistos41(),
    fErr40pos(), fErr40neg(), fErr40(),
    fErr41pos(), fErr41neg(), fErr41(),
    fH2SigYield(NULL), fH2UnfYield(NULL),
    fH2TrueBkg(NULL), fH2MainCS(NULL),
    fH2AvgCS(NULL)
  {
    this->Zero();
  }

  const TH1D* getRawHisto40(int i) const { return fRawHistos40[i]; }
  const TH1D* getRawHisto41(int i) const { return fRawHistos41[i]; }

  const TH1D* getRawHisto(int ibin, int i) const {
    return (ibin==40) ? fRawHistos40[i] : fRawHistos41[i];
  }

  double getCenter(int ibin, TValueType_t valType, int &ierr) const;

  double getErrPos(int ibin, int i) const {
    return (ibin==40) ? fErr40pos[i] : fErr41pos[i];
  }

  double getErrNeg(int ibin, int i) const {
    return (ibin==40) ? fErr40neg[i] : fErr41neg[i];
  }

  double getErrRMS(int ibin, int i) const {
    return (ibin==40) ? fErr40[i] : fErr41[i];
  }

  const TH2D* getSigYield() const { return fH2SigYield; }
  const TH2D* getUnfYield() const { return fH2UnfYield; }
  const TH2D* getTrueBkg() const { return fH2TrueBkg; }
  const TH2D* getMainCS() const { return fH2MainCS; }
  const TH2D* getAvgCS() const { return fH2AvgCS; }


  void Clear();
  void Zero();

  // init is called by getSidedErrors
  int init(TString fname, TString histoTag);
  // prepare the object
  int getSidedErrors(TString fname, TString histoTag);

  void PrintValues(int ibin, TValueType_t valType) const;
};


// -------------------------------------------------------------
// -------------------------------------------------------------

struct SigErrStudyInfo_t {
  std::vector<TH1D*> fRawHistos41;
  double fErr41[_syvalLast];
  TH2D *fH2SigYield, *fH2UnfYield, *fH2MainCS, *fH2AvgCS;
public:
  SigErrStudyInfo_t() :
    fRawHistos41(),
    fErr41(),
    fH2SigYield(NULL), fH2UnfYield(NULL),
    fH2MainCS(NULL),
    fH2AvgCS(NULL)
  {
    this->Zero();
  }

  const TH1D* getRawHisto41(int i) const { return fRawHistos41[i]; }
  TH1D* rawHisto41(int i) { return fRawHistos41[i]; }

  double getCenter(TSigYieldValueType_t valType) const
         { return fRawHistos41[valType]->GetMean(); }
  double getCenter(int i) const
         { return fRawHistos41[i]->GetMean(); }
  double getErrRMS(int i) const { return fRawHistos41[i]->GetRMS(); }


  const TH2D* getSigYield() const { return fH2SigYield; }
  const TH2D* getUnfYield() const { return fH2UnfYield; }
  const TH2D* getMainCS() const { return fH2MainCS; }
  const TH2D* getAvgCS() const { return fH2AvgCS; }


  void Clear();
  void Zero();

  // init is called by getSidedErrors
  int init(TString fname, TString histoTag);

  void PrintValues(TSigYieldValueType_t valType) const;
  void PrintValuesAll() const;
};


// -------------------------------------------------------------
// -------------------------------------------------------------







#endif
