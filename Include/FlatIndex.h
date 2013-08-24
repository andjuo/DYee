#ifndef FlatIndex_H
#define FlatIndex_H

#include "../Include/DYTools.hh"
#include "../Include/AccessOrigNtuples.hh"
//#include 


struct FlatIndex_t {
  int FIdxMass, FIdxY, FFIdx;
public:
  FlatIndex_t () : FIdxMass(-1), FIdxY(-1), FFIdx(-1) {}
  FlatIndex_t (const FlatIndex_t &fi) : FIdxMass(fi.FIdxMass), FIdxY(fi.FIdxY), FFIdx(fi.FFIdx) {}

  // --------------

  int iM() const { return FIdxMass; }
  int iY() const { return FIdxY; }
  int idx() const { return FFIdx; }
  int isValid() const { return (FFIdx!=-1) ? 1:0; }


  // --------------

  // --------------

  int setGenPreFsrIdx(const mithep::TGenInfo *gen) {
    FIdxMass=DYTools::findMassBin(gen->vmass);
    FIdxY   =DYTools::findAbsYBin(FIdxMass,gen->vy);
    FFIdx   =DYTools::findIndexFlat(FIdxMass,FIdxY);
    return FFIdx;
  }

  // --------------

  int setGenPostFsrIdx(const mithep::TGenInfo *gen) {
    FIdxMass=DYTools::findMassBin(gen->mass);
    FIdxY   =DYTools::findAbsYBin(FIdxMass,gen->y);
    FFIdx   =DYTools::findIndexFlat(FIdxMass,FIdxY);
    return FFIdx;
  }

  // --------------

  int setGenPreFsrIdx(const AccessOrigNtuples_t &ai) { return this->setGenPreFsrIdx(ai.genPtr()); }
  int setGenPostFsrIdx(const AccessOrigNtuples_t &ai) { return this->setGenPostFsrIdx(ai.genPtr()); }

  // --------------

  int setRecoIdx(const mithep::TDielectron *dielectron) {
    FIdxMass=DYTools::findMassBin(dielectron->mass);
    FIdxY   =DYTools::findAbsYBin(FIdxMass, dielectron->y);
    FFIdx   =DYTools::findIndexFlat(FIdxMass,FIdxY);
    return FFIdx;
  }

  // --------------

  void print() const {
    std::cout << *this;
  }

  friend std::ostream& operator<<(std::ostream& out, const FlatIndex_t &fi) {
    out << "FlatIdx(iM=" << fi.FIdxMass << ", iY=" << fi.FIdxY << ", idx=" << fi.FFIdx << ")";
    return out;
  }
};



#endif

