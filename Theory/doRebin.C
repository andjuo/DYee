// A macro to rebin 2D theory prediction from Alexey

#include <TROOT.h>
#include <TH1D.h>
#include <TString.h>
#include <TFile.h>
#include <iostream>

int doRebin(TString fname) {
  TFile file(fname,"read");
  if (!file.IsOpen()) {
    std::cout << "failed to open the file <" << fname << ">\n";
    return 0;
  }
  TH1D *h1=(TH1D*)file.Get("xsec");
  h1->SetDirectory(0);
  file.Close();
  h1->SetName("xsec_inp");

  int lastMassBin=(fname.Index("200to1500")>0) ? 1:0;
  if (lastMassBin) std::cout << "lastMassBin determined\n";

  TString htitle=fname;
  htitle.ReplaceAll(".root","");
  int nYbins=(lastMassBin) ? 12 : 24;
  TH1D *hout=new TH1D("xsec",htitle,nYbins,0.,2.4);
  hout->SetDirectory(0);
  hout->GetXaxis()->SetTitle( h1->GetXaxis()->GetTitle() );
  hout->GetXaxis()->SetNdivisions( h1->GetXaxis()->GetNdivisions() );
  hout->GetYaxis()->SetTitle( h1->GetYaxis()->GetTitle() );

  if (!lastMassBin) {
    double factor=20.;  // 2/0.1
    for (int ibin=1; ibin<=h1->GetNbinsX(); ++ibin) {
      hout->SetBinContent(ibin, factor * h1->GetBinContent(ibin));
      hout->SetBinError  (ibin, factor * h1->GetBinError(ibin));
    }
  }
  else {
    double factor=10.;  // 2/0.2
    for (int ibin=1; ibin<=hout->GetNbinsX(); ++ibin) {
      double v1=h1->GetBinContent(2*ibin-1);
      double v2=h1->GetBinContent(2*ibin );
      double e1=h1->GetBinError (2*ibin-1);
      double e2=h1->GetBinError (2*ibin  );
      hout->SetBinContent(ibin, factor * (v1+v2) );
      hout->SetBinError  (ibin, factor * sqrt(e1*e1+e2*e2));
    }
  }

  TString foutName=fname;
  foutName.ReplaceAll(".root","-mdf.root");
  TFile fout(foutName,"recreate");
  hout->Write();
  fout.Close();
  std::cout << "file <" << fout.GetName() << "> created\n";
  return 1;
}
