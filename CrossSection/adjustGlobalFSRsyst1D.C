#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/ComparisonPlot.hh"
#include "CSCovWorkFlags.hh" // for comparison with main result


// ---------------------------------------------------

void adjustGlobalFSRsyst1D(int iSave=0) {

  if (!DYTools::setup(0)) return;

  TString inpFile="inp_globalFSRsyst_1D.root";
  TFile fin(inpFile,"read");
  TH2D *h2Err=createBaseH2("totErr","",1);
  TH2D *h2CS =createBaseH2("mainCS_div100","",1);
  if (!loadHisto(fin,&h2Err,"") ||
      !loadHisto(fin,&h2CS,"")) {
    std::cout << "failed to load histos\n";
    return;
  }
  fin.Close();

  // ---------------------
  // Manual adjustment
  // ---------------------

  TH2D* h2ErrMdf=Clone(h2Err,"h2ErrMdf");

  printHisto(h2Err);

  int ibin=0;
  int jbin=1;
  //ibin=7;
  //h2ErrMdf->SetBinContent(ibin,jbin, 0.5*h2CS->GetBinContent(ibin,jbin));
  //ibin=8;
  //h2ErrMdf->SetBinContent(ibin,jbin, 0.5*h2CS->GetBinContent(ibin,jbin));
  //ibin=9;
  //h2ErrMdf->SetBinContent(ibin,jbin, 0.25*h2CS->GetBinContent(ibin,jbin));
  //ibin=16;
  //h2ErrMdf->SetBinContent(ibin,jbin, 1.2*h2CS->GetBinContent(ibin,jbin));
  ibin=18;
  h2ErrMdf->SetBinContent(ibin,jbin, 0.8*h2CS->GetBinContent(ibin,jbin));
  ibin=20;
  h2ErrMdf->SetBinContent(ibin,jbin, 0.6*h2CS->GetBinContent(ibin,jbin));
  ibin=22;
  h2ErrMdf->SetBinContent(ibin,jbin, 0.48*h2CS->GetBinContent(ibin,jbin));
  ibin=24;
  h2ErrMdf->SetBinContent(ibin,jbin, 0.44*h2CS->GetBinContent(ibin,jbin));
  ibin=26;
  h2ErrMdf->SetBinContent(ibin,jbin, 0.46*h2CS->GetBinContent(ibin,jbin));
  ibin=28;
  h2ErrMdf->SetBinContent(ibin,jbin, 0.48*h2CS->GetBinContent(ibin,jbin));
  ibin=31;
  h2ErrMdf->SetBinContent(ibin,jbin, 0.42*h2CS->GetBinContent(ibin,jbin));
  ibin=33;
  h2ErrMdf->SetBinContent(ibin,jbin, 0.35*h2CS->GetBinContent(ibin,jbin));
   //linearApprox(hErrMdf, 2.11, 2.31);

  std::cout << "absolute additional errors errors of the scale factors\n";
  printHisto(h2ErrMdf);


  // ----------------------
  // Display the result
  // ----------------------

  TH2D* h2RelErr=Clone(h2Err,"h2RelErr");
  divide(h2RelErr,h2CS);
  TH2D* h2RelErrMdf=Clone(h2ErrMdf,"h2RelErrMdf");
  divide(h2RelErrMdf,h2CS);

  std::vector<TH2D*> hAbsV, hRelV;
  std::vector<TString> labelV;

  hAbsV.push_back(h2Err);
  hRelV.push_back(h2RelErr);
  labelV.push_back("orig");
  hAbsV.push_back(h2ErrMdf);
  hRelV.push_back(h2RelErrMdf);
  labelV.push_back("assigned");

  TCanvas *cx= plotProfiles("cx",hAbsV,labelV,NULL,1,"abs error");
  TCanvas *cy= plotProfiles("cy",hRelV,labelV,NULL,1,"rel error");

  cx->Update();
  cy->Update();

  /*
  if (iSave) {
    TString newFName=fname;
    newFName.ReplaceAll(".root","-mdf.root");
    TFile fout(newFName,"recreate");
    for (unsigned int i=0; i<hVar.size(); ++i) {
      hVar[i]->Write(fieldNamesV[i]);
    }
    hErrMdf->Write("rhoAbsSyst");
    hRelErrMdf->Write("rhoRelSyst");
    writeBinningArrays(fout);
    fout.Close();
    std::cout << "file <" << fout.GetName() << "> saved\n";
  }
  */

  return;
}

// ------------------------------------------------------------------
// ------------------------------------------------------------------

