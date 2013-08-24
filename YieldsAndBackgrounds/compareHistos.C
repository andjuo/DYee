#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>

#include "../Include/DYTools.hh"
#include "../Include/InputFileMgr.hh"

void compareHistos(TString conf, TString fName1, TString fName2) {
  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(conf)) return;

  TFile file1(fName1,"read");
  TFile file2(fName2,"read");

  for (unsigned int i=0; i<inpMgr.sampleCount(); ++i) {
    //if (i) break;
    TH1F *h1=(TH1F*)file1.Get(inpMgr.sampleName(i));
    TH1D *h2=(TH1D*)file2.Get(inpMgr.sampleName(i));

    if (!h1 || !h2) std::cout << "failed to get " << inpMgr.sampleName(i) << "\n";

    std::cout << "sample " << inpMgr.sampleName(i) << "\n";
    int ok=1;
    if (h1->GetNbinsX() != h2->GetNbinsX()) {
      std::cout << "\tbin count is different\n";
      ok=0;
    }
    else {
      for (int ibin=1; ibin<=h1->GetNbinsX(); ++ibin) {
	double v1=h1->GetBinContent(ibin);
	double v2=h2->GetBinContent(ibin);
	
	double sum=v1+v2;
	if (sum>0.) {
	  double rdiff=fabs(v1-v2)/sum;
	  if (rdiff > 1e-2) {
	    ok=0;
	    std::cout << "ibin=" << ibin << ", v1=" << v1 << ", v2=" << v2 << "(rdiff=" << rdiff << ")\n";
	  }
	}
      }
    }
    delete h1;
    delete h2;
    if (ok) std::cout << "\tsample ok\n";
  }
  file1.Close();
  file2.Close();
}

