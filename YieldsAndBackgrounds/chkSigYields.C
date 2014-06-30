#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"

void chkSigYields() {
  if (!DYTools::setup(0)) return;

  TString path="/media/sdb/andriusj/Results-DYee-escaleRndFlat-20140601/root_files_reg/yield/DY_j22_19712pb_EScale_study_randomized/";
  TString fname1="bg-subtracted_yield_1D__peak20140220flat_DtRND1018.root";
  TString fname2="bg-subtracted_yield_1D__peak20140220flat_DtRND1014.root";
  TString fname3="bg-subtracted_yield_1D__peak20140220flat_DtRND1021.root";
  TString field="signalYieldDDbkg";
  TH2D* h2a=LoadHisto2D(field,path+fname1,"",1);
  TH2D* h2b=LoadHisto2D(field,path+fname2,"",1);
  TH2D* h2c=LoadHisto2D(field,path+fname3,"",1);

  std::vector<TH2D*> hV;
  std::vector<TString> lV;

  hV.push_back(h2a); lV.push_back("a");
  hV.push_back(h2b); lV.push_back("b");
  hV.push_back(h2c); lV.push_back("c");

  TCanvas *cx= plotProfiles("cx",hV,lV);
  cx->Update();
}
