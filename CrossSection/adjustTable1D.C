#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "CSCovWorkFlags.hh"

void adjustTable1D(int saveLatex=0, int saveStatErr=0) {
  if (!DYTools::setup(0)) return;
  TString fileTag="frac";
  std::vector<TH2D*> histosV_old,histosV;
  std::vector<TString> labelsV_old,labelsV;

  if (0) {
    // check
    if (!loadLatexTableTextFile(fileTag,histosV,labelsV,0)) return;
    return;
  }

  // modify
  if (!loadLatexTableTextFile(fileTag,histosV_old,labelsV_old,0)) return;

  HERE("loaded ok");

  histosV.reserve(histosV_old.size());
  labelsV.reserve(labelsV_old.size());

  unsigned int imin=(saveStatErr) ? 0 : 1;
  for (unsigned int i=imin; i<histosV_old.size()-1; ++i) {
    histosV.push_back(Clone(histosV_old[i],Form("clone_%d",i)));
    labelsV.push_back(labelsV_old[i]);
    std::cout << "storing " << labelsV.back() << "\n";
  }

  HERE("calculating the total");

  TH2D* h2=Clone(histosV[0],"h2Sum");
  h2->Reset();
  for (int ibin=1; ibin<=h2->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=h2->GetNbinsY(); ++jbin) {
      double sum=0;
      for (unsigned int i=0; i<histosV.size(); ++i) {
	sum+= pow(histosV[i]->GetBinContent(ibin,jbin),2);
      }
      h2->SetBinContent(ibin,jbin, sqrt(sum));
    }
  }
  histosV.push_back(h2);
  labelsV.push_back("total error");


  if (saveLatex) {
    HERE("replacing for latex");

    replaceAll(labelsV,"Mass (GeV)","$m$//(GeV)");
    replaceAll(labelsV,"signal stat","signal//stat");
    replaceAll(labelsV,"signal syst","Bkgr.est.//($\\%$)");
    replaceAll(labelsV,"signal EScale uncert.","E-scale//($\\%$)");
    replaceAll(labelsV,"unf stat","Det.resol.//($\\%$)");
    replaceAll(labelsV,"eff stat","Eff.//($\\%$)");
    replaceAll(labelsV,"ESF tot","$\\rho$//($\\%$)");
    replaceAll(labelsV,"acc stat","Acc.stat.//($\\%$)");
    replaceAll(labelsV,"FSR stat","FSR.unf.//($\\%$)");
    replaceAll(labelsV,"puRndStudy","Coll.CS//($\\%$)");
    replaceAll(labelsV,"fsrRndStudy","FSR model//($\\%$)");
    replaceAll(labelsV,"total error","Total//($\\%$)");

    for (unsigned int i=0; i<labelsV.size(); ++i) {
      std::cout << " - " << i << "  " << labelsV[i] << "\n";
    }

    HERE("saving");

    if (!saveLatexTable(fileTag + TString("-noStat"),histosV,labelsV,"%5.2lf",0,0)) return;

  }
  else {
    HERE(dashline);
    HERE("replacing");

    //replaceAll(labelsV,"Mass (GeV)","mass");
    replaceAll(labelsV,"signal stat","signal_stat_err");
    replaceAll(labelsV,"signal syst","bkgr_est_err");
    replaceAll(labelsV,"signal EScale uncert.","escale_err");
    replaceAll(labelsV,"unf stat","det_resolution_err");
    replaceAll(labelsV,"eff stat","eff_rnd_err");
    replaceAll(labelsV,"ESF tot","rho_err");
    replaceAll(labelsV,"acc stat","acc_rnd_err");
    replaceAll(labelsV,"FSR stat","fsr_rnd_err");
    replaceAll(labelsV,"puRndStudy","pileup_err");
    replaceAll(labelsV,"fsrRndStudy","fsr_model_err");
    replaceAll(labelsV,"total error","total_err");

    for (unsigned int i=0; i<labelsV.size(); ++i) {
      std::cout << " - " << i << "  " << labelsV[i] << "\n";
    }

    TH2D* hCS=loadMainCSResult(1);
    if (!hCS) return;
    TH1D* h1CS=createProfileAuto(hCS,1,"xsec");
    if (!h1CS) return;

    TString fname="table_1D_frac.root";
    TFile fout(fname,"recreate");
    if (!saveHisto(fout,h1CS,"","xsec")) return;
    for (unsigned int i=0; i<histosV.size(); ++i) {
      TH1D* h1=createProfileAuto(histosV[i],1,labelsV[i]);
      h1->Scale(0.01);
      if (!saveHisto(fout,h1,"",labelsV[i])) return;
    }
    writeBinningArrays(fout,"adjustTable1D",0);
    fout.Close();
    std::cout << "file <" << fout.GetName() << "> closed\n";
  }

  return;
}
