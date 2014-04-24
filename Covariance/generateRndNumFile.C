#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TString.h>
#include <TObjString.h>
#include <TRandom.h>
#include <iostream>

void generateRndNumFile(int iSeed, int nNumbers, TString outFName) {
  gRandom->SetSeed(iSeed);


  TFile fout(outFName,"recreate");
  if (!fout.IsOpen()) {
    std::cout << "failed to create the file <" << outFName << ">\n";
    return;
  }

  double r;
  TTree *tree =new TTree("rndSequence","rndSequence");
  tree->Branch("rnd",&r,"rnd/D");
  for (int i=0; i<nNumbers; ++i) {
    r=gRandom->Gaus(0.,1.);
    tree->Fill();
  }
  fout.cd();
  tree->Write();
  TObjString info(Form("seed_%d__%d",iSeed,nNumbers));
  info.Write();
  fout.Close();
  return;
}

