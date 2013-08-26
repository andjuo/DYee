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

  int setIdx(double mass, double y) {
    FIdxMass=DYTools::findMassBin(mass);
    FIdxY   =DYTools::findAbsYBin(FIdxMass,y);
    FFIdx   =DYTools::findIndexFlat(FIdxMass,FIdxY);
    return FFIdx;
  }

  // --------------

  int setGenPreFsrIdx(const mithep::TGenInfo *gen) {
    return this->setIdx(gen->vmass, gen->vy);
  }

  // --------------

  int setGenPostFsrIdx(const mithep::TGenInfo *gen) {
    return this->setIdx(gen->mass, gen->y);
  }

  // --------------

  int setGenPreFsrIdx(const AccessOrigNtuples_t &ai) { return this->setGenPreFsrIdx(ai.genPtr()); }
  int setGenPostFsrIdx(const AccessOrigNtuples_t &ai) { return this->setGenPostFsrIdx(ai.genPtr()); }

  // --------------

  int setRecoIdx(const mithep::TDielectron *dielectron) {
    return this->setIdx(dielectron->mass, dielectron->y);
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

