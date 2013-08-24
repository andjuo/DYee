#include <TVectorD.h>
#include <TMatrixD.h>

#include <fstream>
#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/UnfoldingTools.hh"

typedef TH1F Histo_t; // match the type returned in MyTools.hh

// -----------------------------------------

void AdjustDim(std::string &line) {
  const std::string expect =( DYTools::study2D) ? "2D" : "1D";
  const std::string canHave=(!DYTools::study2D) ? "2D" : "1D";
  if (!PosOk(line,expect)) {
    size_t pos=line.find(canHave);
    if (PosOk(pos)) {
      std::cout << "changing <" << line << ">\n";
      line[pos]=expect[0];
      std::cout << "      to <" << line << ">\n";
    }
  }
  return;
}

// -----------------------------------------

int LoadMatrix(const TString &fname, TMatrixD &M, TMatrixD &Merr, const TString &field, const TString &fieldErr) {
  TFile f(fname,"read");
  TMatrixD *Mptr=(TMatrixD*)f.Get(field);
  TMatrixD *MerrPtr=(TMatrixD*)f.Get(fieldErr);
  f.Close();
  if (!Mptr || !MerrPtr) {
    std::cout << "error in loading from <" << fname << ">\n";
    if (!Mptr) std::cout << " failed to load <" << field << ">\n";
    if (!MerrPtr) std::cout << " failed to load <" << fieldErr << ">\n";
    return 0;
  }
  if ((Mptr->GetNrows() != M.GetNrows()) ||
      (Mptr->GetNcols() != M.GetNcols()) ) {
    std::cout << "dim mismatch in <" << fname << ">\n";
    return 0;
  }
  M = *Mptr;
  Merr= *MerrPtr;
  delete Mptr;
  delete MerrPtr;
  return 1;
}

// -----------------------------------------

int GetMassProfile(const TMatrixD &m, int yIdx, TVectorD &v) {
  v=0;
  for (int i=0; i<m.GetNrows(); ++i) {
    v[i] = m(i,yIdx);
  }
  return 1;
}

// -----------------------------------------

int GetMassProfile1D(const TMatrixD &m, const TMatrixD &mErr,
		     TVectorD &v, TVectorD &vErr) {
  int res=
    GetMassProfile(m,0,v) &&
    GetMassProfile(mErr,0,vErr);
  return res;
}

// -----------------------------------------

// -----------------------------------------
// -----------------------------------------

void gatherA_spec() {
  if (DYTools::study2D) {
    std::cout << "The code is for 1D\n";
    return;
  }

  TMatrixD baseM(DYTools::nMassBins,DYTools::nYBinsMax);
  TMatrixD M(baseM), MErr(baseM);
  //TVectorD baseV(DYTools::nMassBin);
  //TVectorD V(baseV), VErr(baseV);

  TString outFName = "acc_check_20130821_1D.root";

  TFile fout(outFName,"recreate");

  // pre FSR without FEWZ
  
  TString inpName1="../root_files/constants_preFsr/DY_j22_19789pb/acceptance_constants_1D.root";
  if (!LoadMatrix(inpName1, M,MErr, "acceptanceMatrix_noFEWZ","acceptanceMatrix_noFEWZ_err")) return;

  Histo_t *h1=extractMassDependence("hAcc_postFSR_noFEWZ","",M,MErr,0,0,0);
  fout.cd();
  h1->Write();

  if (!LoadMatrix(inpName1, M,MErr, "acceptanceMatrixPreFsr_noFEWZ","acceptanceMatrixPreFsrErr_noFEWZ")) return;
  Histo_t *h2=extractMassDependence("hAcc_preFSR_noFEWZ","",M,MErr,0,0,0);
  fout.cd();
  h2->Write();

  if (!LoadMatrix(inpName1, M,MErr, "nEventsPreFsr_noFEWZ","nEventsPreFsr_noFEWZ_err")) return;
  Histo_t *h3=extractMassDependence("hNEvents_preFSR_noFEWZ","",M,MErr,0,0,0);
  fout.cd();
  h3->Write();

  fout.Close();

  std::cout << "file <" << outFName << "> created\n";

  return;
}
