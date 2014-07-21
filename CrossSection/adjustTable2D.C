#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include <ostream>
#include <fstream>
#include "CSCovWorkFlags.hh"


void adjustTable2D(int saveLatex=0, int saveStatErr=0,
		   int nBayesIters=-1) {
  if (!DYTools::setup(1)) return;
  TString fileTag="frac-total";
  if (nBayesIters!=-1) fileTag=Form("frac_nBayes%d",nBayesIters);
  std::vector<TH2D*> histosV_old,histosV;
  std::vector<TString> labelsV_old,labelsV;

  if (0) {
    // check
    if (!loadLatexTableTextFile(fileTag,histosV,labelsV,0)) return;

    HERE("saving (check)");

    TString extraTag=(saveStatErr) ? "-wStat" : "-noStat";
    if (!saveLatexTable(fileTag + extraTag,histosV,labelsV,"%5.2lf",0,0))
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
  TH2D* h2Syst=Clone(histosV[0],"h2SumSystOnly");
  h2->Reset();
  h2Syst->Reset();
  for (int ibin=1; ibin<=h2->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=h2->GetNbinsY(); ++jbin) {
      double sum=0;
      double sumSyst=0;
      for (unsigned int i=0; i<histosV.size(); ++i) {
	double term2= pow(histosV[i]->GetBinContent(ibin,jbin),2);
	sum+= term2;
	if (!saveStatErr || (saveStatErr && (i>0))) {
	  sumSyst += term2;
	}
      }
      h2->SetBinContent(ibin,jbin, sqrt(sum));
      h2Syst->SetBinContent(ibin,jbin, sqrt(sumSyst));
    }
  }
  if (saveStatErr) {
    histosV.push_back(h2Syst);
    labelsV.push_back("total syst error");
  }
  histosV.push_back(h2);
  labelsV.push_back("total error");


  if (!saveLatex) {
    HERE(dashline);
    HERE("replacing for root file");

    //replaceAll(labelsV,"Mass (GeV)","mass");
    replaceAll(labelsV,"signal stat","signal_stat_err");
    replaceAll(labelsV,"signal syst","bkgr_est_err");
    replaceAll(labelsV,"signal EScale uncert.","escale_err");
    replaceAll(labelsV,"unf stat","det_resolution_err");
    replaceAll(labelsV,"unf e-scale residual","unf_escale_res");
    replaceAll(labelsV,"eff stat","eff_rnd_err");
    replaceAll(labelsV,"ESF tot","rho_err");
    //replaceAll(labelsV,"acc stat","acc_rnd_err");
    replaceAll(labelsV,"FSR stat","fsr_rnd_err");
    replaceAll(labelsV,"puRndStudy","pileup_err");
    replaceAll(labelsV,"fsrRndStudy","fsr_model_err");
    replaceAll(labelsV,"total syst error","total_syst_err");
    replaceAll(labelsV,"total error","total_err");

    for (unsigned int i=0; i<labelsV.size(); ++i) {
      std::cout << " - " << i << "  " << labelsV[i] << "\n";
    }

    TH2D* hCS=loadMainCSResult(1,NULL,nBayesIters);
    if (!hCS) return;

    TString fname="table_2D_frac.root";
    if (nBayesIters!=-1) {
      fname.ReplaceAll("_frac",Form("_frac_nBayes%d",nBayesIters));
    }
    TFile fout(fname,"recreate");

    TH2D* hCSok=removeUnderflow(hCS,"hCSok");
    printHisto(hCSok);

    if (!saveHisto(fout,hCSok,"","xsec")) return;
    for (unsigned int i=0; i<histosV.size(); ++i) {
      TH2D* h2a=removeUnderflow(histosV[i],labelsV[i]);
      h2a->Scale(0.01);
      printHisto(h2a);
      if (!saveHisto(fout,h2a,"",labelsV[i])) return;
    }
    writeBinningArrays(fout,"adjustTable2D",0);
    fout.Close();
    std::cout << "file <" << fout.GetName() << "> closed\n";

    return;
  }

  // -------------------------
  // Save latex table
  // -------------------------

  HERE("replacing");

  replaceAll(labelsV,"Mass (GeV)","$m$//(GeV)");
  replaceAll(labelsV,"signal stat","signal//stat");
  replaceAll(labelsV,"signal syst","Bkgr.est.//($\\%$)");
  replaceAll(labelsV,"signal EScale uncert.","E-scale//($\\%$)");
  replaceAll(labelsV,"unf stat","Det.resol.//($\\%$)");
  replaceAll(labelsV,"unf e-scale residual","E-Scale//res.($\\%$)");
  replaceAll(labelsV,"eff stat","Eff.//($\\%$)");
  replaceAll(labelsV,"ESF tot","$\\rho$//($\\%$)");
  replaceAll(labelsV,"acc stat","Acc.stat.//($\\%$)");
  replaceAll(labelsV,"FSR stat","FSR.//unf.($\\%$)");
  replaceAll(labelsV,"puRndStudy","Coll.//CS($\\%$)");
  replaceAll(labelsV,"fsrRndStudy","FSR model//($\\%$)");
  replaceAll(labelsV,"total syst error","Total//syst. ($\\%$)");
  replaceAll(labelsV,"total error","Total//($\\%$)");

  for (unsigned int i=0; i<labelsV.size(); ++i) {
    std::cout << " - " << i << "  " << labelsV[i] << "\n";
  }

  HERE("saving");

  if (0) {
    if (!saveLatexTable(fileTag + TString("-noStat"),histosV,labelsV,"%5.2lf",0,0)) return;
  }

  TString extraTag=(saveStatErr) ? "-wStat" : "-noStat";
  TString fname="table" + DYTools::analysisTag +
    fileTag + extraTag + TString(".tex");
  std::ofstream fout(fname);
  TMatrixD *yBinLimits= DYTools::getYBinLimits();

  std::vector<TString> tlabelsV,subLabelsV;
  int hasSubLabels=0;
  tlabelsV.reserve(labelsV.size());
  subLabelsV.reserve(labelsV.size());
  for (unsigned int i=0; i<labelsV.size(); i++) {
    TString txt=labelsV[i];
    int p=txt.Index("//");
    if (p!=-1) {
      hasSubLabels=1;
      tlabelsV.push_back(txt(0,p));
      subLabelsV.push_back(txt(p+2,txt.Length()));
    }
    else {
      tlabelsV.push_back(txt);
      subLabelsV.push_back("");
    }
  }
  const char *comment=(1) ? "%" : "";

  fout << comment << "\\documentclass{article}\n";
  fout << comment << "\\usepackage{rotating}\n";
  fout << comment << "\\newcommand{\\GeV}{\\ GeV\\ }\n";
  fout << comment << "\\begin{document}\n";
  const char *format="%5.2lf ";

  if (1) {
    // Rotation is needed
    for (int itable=0; itable<6; ++itable) {
      fout << "\\begin{sidewaystable}[tbhp]\n";
      fout << "%\\begin{table}[tbhp]\n";
      fout << "\\caption{\\label{tbl-" << itable << "} ";
      fout << "Summary of systematic uncertainties in the dielectron channel for ";
      fout << "$" << DYTools::massBinLimits[itable+1] << "<m<"
	   << DYTools::massBinLimits[itable+2] << "$\\GeV bin as a function"
	   << " of $|y|$. The \"Total\" is a quadratic sum of all sources.}\n";

      if (strlen(comment)==0) fout << "\\hspace{-1cm}\n";
      fout << "\\begin{tabular}{|";

      int count=int(histosV.size());
      for (int i=0; i<count+1; ++i) fout << "c|";
      fout << "}\n";
      fout << "\\hline\n";

      fout << " $|y|$ & ";

      for (int i=0; i<count; ++i) {
	fout << " " << tlabelsV[i];
	if (i!=count-1) fout << " &";
      }
      fout << "\\\\\n";
      fout << "       &";
      for (int i=0; i<count; ++i) {
	fout << " " << subLabelsV[i];
	if (i!=count-1) fout << " &";
      }
      fout << "\\\\\n";
      fout << "\\hline\n";

      int ibinMin=itable+1;
      for (int ibin=ibinMin+1; ibin<ibinMin+2; ++ibin) {
	int yRangeCount= int(histosV[0]->GetNbinsY());
	if (DYTools::study2D && (ibin==DYTools::nMassBins)) yRangeCount=12;

	for (int jbin=1; jbin<=yRangeCount; ++jbin) {
	  fout << " " << Form("%3.1lf",(*yBinLimits)(ibin-1,jbin-1)) << "-"
	       << Form("%3.1lf",(*yBinLimits)(ibin-1,jbin)) << " & ";

	  for (int ih=0; ih<count; ++ih) {
	    fout << " " << Form(format,histosV[ih]->GetBinContent(ibin,jbin));
	    if (ih!=count-1) fout << " &";
	  }
	  fout << "\\\\\\hline\n";
	}
      }
      fout << "\\hline\n";
      fout << "\\end{tabular}\n";
      fout << "%\\end{table}\n";
      fout << "\\end{sidewaystable}\n";
      fout << "\n\n";
    }

  }
  else {
    // rotation is not needed
  for (int itable=0; itable<3; ++itable) {
    fout << "\\begin{sidewaystable}[tbhp]\n";
    fout << "%\\begin{table}[tbhp]\n";
    fout << "\\caption{\\label{tbl-" << itable << "} ";
    fout << "Summary of systematic uncertainties in the dielectron channel for ";
    fout << "$" << DYTools::massBinLimits[2*itable+1] << "<m<"
	 << DYTools::massBinLimits[2*itable+2] << "$\\GeV and $"
	 << DYTools::massBinLimits[2*itable+2] << "<m<"
	 << DYTools::massBinLimits[2*itable+3] << "$\\GeV bins as a function"
	 << " of $|y|$. The \"Total\" is a quadratic sum of all sources.}\n";

    if (strlen(comment)==0) fout << "\\hspace{-1cm}\n";
    fout << "\\begin{tabular}{|";

    int count=int(histosV.size());
    for (int i=0; i<count+1; ++i) fout << "c|";
    fout << "}\n";
    fout << "\\hline\n";

    fout << " $|y|$ & ";

    for (int i=0; i<count; ++i) {
      fout << " " << tlabelsV[i];
      if (i!=count-1) fout << " &";
    }
    fout << "\\\\\n";
    fout << "       &";
    for (int i=0; i<count; ++i) {
      fout << " " << subLabelsV[i];
      if (i!=count-1) fout << " &";
    }
    fout << "\\\\\n";
    fout << "\\hline\n";

    int ibinMin=2*itable+1;
    for (int ibin=ibinMin+1; ibin<=ibinMin+2; ++ibin) {
      int yRangeCount= int(histosV[0]->GetNbinsY());
      if (DYTools::study2D && (ibin==DYTools::nMassBins)) yRangeCount=12;
      fout << Form("   & \\multicolumn{%lu}{c|}{$",histosV.size())
	   << DYTools::massBinLimits[ibin-1] << "<m<"
	   << DYTools::massBinLimits[ibin] << "$\\GeV} \\\\\\hline\n";

      for (int jbin=1; jbin<=yRangeCount; ++jbin) {
	fout << " " << Form("%3.1lf",(*yBinLimits)(ibin-1,jbin-1)) << "-"
	     << Form("%3.1lf",(*yBinLimits)(ibin-1,jbin)) << " & ";

	for (int ih=0; ih<count; ++ih) {
	  fout << " " << Form(format,histosV[ih]->GetBinContent(ibin,jbin));
	  if (ih!=count-1) fout << " &";
	}
	fout << "\\\\\\hline\n";
      }
    }
    fout << "\\hline\n";
    fout << "\\end{tabular}\n";
    fout << "%\\end{table}\n";
    fout << "\\end{sidewaystable}\n";
    fout << "\n\n";
  }
  }
  fout << comment << "\\end{document}\n";
  fout.close();

  delete yBinLimits;
  std::cout << " saved file " << fname << "\n";

}
