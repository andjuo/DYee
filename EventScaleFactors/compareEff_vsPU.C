//
// Made on May 08, 2014
// to plot pu-dependence of the efficiency
//

#include "../Include/DYTools.hh"
#include "calcEventEffLink.h"
#include <TGraphAsymmErrors.h>
#include "../Include/ComparisonPlot.hh"
#include "../Include/colorPalettes.hh"

// ------------------------------------------------------------


// ------------------------------------------------------------
// ------------------------------------------------------------
// ------------------------------------------------------------

TString effDataKindString(const TString str) {
  TString effKind="xx", dataKind="xx";
  if (str.Index("RECO")!=-1) effKind="RECO";
  else if (str.Index("ID")!=-1) effKind="ID";
  else if (str.Index("HLTleg1")!=-1) effKind="HLTleg1";
  else if (str.Index("HLTleg2")!=-1) effKind="HLTleg2";
  else if (str.Index("HLT")!=-1) effKind="HLT";
  if (str.Index("data")!=-1) dataKind="data";
  else if (str.Index("mc")!=-1) dataKind="mc";
  TString final=dataKind+TString(" ")+effKind;
  return final;
}

// ------------------------------------------------------------
// ------------------------------------------------------------

int loadData(TString etEtaStr,
	     int iBr, int iEt, int iEta, int iSelect,
	     TString &explainChoice,
	     std::vector<TGraphAsymmErrors*> &dataV,
	     std::vector<TString> &grLabelsV);

// ------------------------------------------------------------
// ------------------------------------------------------------


void compareEff_vsPU(int iBr=0,
		     int iEt=0, int iEta=0,
		     int doSave=0,
		     double transLegendY_user=0.,
		     TString *outFileName_ptr=NULL,
		     TString *outDir_ptr=NULL,
		     int iSelect=0) {

  TString etEtaStr="EtBins1EtaBins1";
  int iEtMax=1, iEtaMax=1;
  if (iBr>=10) {
    etEtaStr="EtBins1EtaBins5";
    iEtMax=1; iEtaMax=5;
  }

  double transLegendX=0.;
  double transLegendY=transLegendY_user;

  std::vector<TGraphAsymmErrors*> dataV;
  std::vector<TString> grLabelsV;
  TString theCaseStr;

  if (!loadData(etEtaStr,iBr,iEt,iEta, iSelect,
		theCaseStr,dataV,grLabelsV)) return;

  if (iSelect==1) {
    std::cout << "data selected\n";
    return;
  }

  std::cout << theCaseStr << ": got " << dataV.size() << " data graphs with "
	    << grLabelsV.size() << " explanations\n";
  if (dataV.size()==0) return;

  int production=0;
  std::vector<TString> labelsForPlotV;
  labelsForPlotV.reserve(grLabelsV.size());
  for (unsigned int i=0; i<grLabelsV.size(); ++i) {
    TString s=grLabelsV[i];
    TString isMC= (s.Index("data")!=-1) ? "Data " : "MC ";
    int idx=s.Index("_run");
    TString runStr= s(idx+1,s.Length());
    TString final=isMC + runStr;
    if (production && (isMC.Index("Data")==-1)) final="MC";
    labelsForPlotV.push_back(final);
  }

  if (production) {
    transLegendX=-0.1;
    if (transLegendY_user==0.) transLegendY=-0.06;
  }


  TString cpTitle;
  if (iBr>=10) {
    //DYTools::TEtBinSet_t etBinSet=DetermineEtBinSet(etEtaStr);
    DYTools::TEtaBinSet_t etaBinSet=DetermineEtaBinSet(etEtaStr);
    //double *loc_etBinLimits=DYTools::getEtBinLimits(etBinSet);
    double *loc_etaBinLimits=DYTools::getEtaBinLimits(etaBinSet);
    cpTitle= TString(Form("%5.3lf #leq |#eta| #leq %5.3lf",
			  loc_etaBinLimits[iEta],loc_etaBinLimits[iEta+1]));
  }

  //if (vsEt) cpTitle= dataKind+ TString(Form(" %5.3lf #leq %s #leq %5.3lf",loc_etaBinLimits[iBin],(signedEta)?"#eta":"abs(#eta)",loc_etaBinLimits[iBin+1]));
  //else cpTitle= dataKind+ TString(Form(" %2.0lf #leq #it{E}_{T} #leq %2.0lf GeV",loc_etBinLimits[iBin],loc_etBinLimits[iBin+1]));
  TString xaxisTitle="nPV";
  TString yaxisTitle=" efficiency";
  if (iBr%10==0) yaxisTitle.Prepend("RECO");
  else if (iBr%10==1) yaxisTitle.Prepend("ID");

  ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,"comp",cpTitle,
			  xaxisTitle,yaxisTitle,"ratio");
  cp.SetYRange(0.9,1.02);
  if (iBr==1) cp.SetYRange(0.73,0.90);
  if (iBr>10) cp.SetYRange(0.6,0.9);

  cp.SetRefIdx(-111); // do not plot lower panel
  cp.SetGridx(true);
  cp.SetGridy(true);

  // Square canvas if ratio is not plotted
  TCanvas *cx=new TCanvas("cx","cx",800,800);

  const int nColors=6;
  const int colors[nColors] = { kBlack, kGreen+2, kBlue, kRed+1, 46, 38 };

  if (production) {
    unsigned int igr=1;
    cp.AddGraph(dataV[igr],labelsForPlotV[igr],"LPE",colors[igr%nColors],1, 1,1,1.);
  }

  for (unsigned int igr=0; igr<dataV.size(); ++igr) {
    if ((igr%2==1) && (production || !production && (igr!=1))) continue;
    cp.AddGraph(dataV[igr],labelsForPlotV[igr],"LPE",colors[igr%nColors],1, 1,1,1.);
  }

  int targetPad=0;
  cp.Draw(cx,0,"png",targetPad);
  cp.TransLegend(transLegendX, transLegendY);
  cp.WidenLegend(0.2,0.);

  /*
  if (fnameTag.Length()) {
    TString fname=TString("fig-eff-") + fnameTag + cpTitle;
    fname.ReplaceAll(" #leq "," ");
    fname.ReplaceAll(" ","_");
    fname.ReplaceAll("(#eta)","Eta");
    fname.ReplaceAll("#eta","eta");
    fname.ReplaceAll(".","_");
    fname.ReplaceAll("#it{E}_{T}","Et");
    //fname.Append(".png");
    TString locOutDir=TString("plots") + fnameTag;
    locOutDir.ReplaceAll("--","");

    std::cout << "fnameBase=<" << fname << "> in <" << locOutDir << ">\n";
    if (outFileName_ptr) *outFileName_ptr=fname;
    if (outDir_ptr) *outDir_ptr=locOutDir;

    if (doSave) {
      SaveCanvas(cx,fname,locOutDir);
    }
    else {
      std::cout <<  " ... not saved (as requested)\n";
    }
  }
  */

  return ;
}


// -------------------------------------------------------------------------
// -------------------------------------------------------------------------

int loadData(TString etEtaStr,
	     int iBr, int iEt, int iEta, int iSelect,
	     TString &explainChoice,
	     std::vector<TGraphAsymmErrors*> &dataV,
	     std::vector<TString> &grLabelsV) {


  TString outFileName=TString("res_effVsPU_") + etEtaStr + TString(".root");
  TString fileOpt=(iSelect) ? "recreate" : "read";
  TFile fout(outFileName,fileOpt);
  if (!fout.IsOpen()) {
    std::cout << "failed to " << fileOpt << " the file <"
	      << fout.GetName() << ">\n";
    return 0;
  }
  if (iSelect) {
    if (!writeFlagValues("nPVLimits",
			 DYTools::nPVBinCount,DYTools::nPVLimits)) return 0;
  }


  for (int iRun=0; iRun<3; ++iRun) {
    if (iBr>=10) iRun=10;
    TString runTag;
    switch(iRun) {
    case 0: runTag="_runAB"; break;
    case 1: runTag="_runC"; break;
    case 2: runTag="_runD"; break;
    case 10: break; // no specialRunTag
    default:
      std::cout << "not ready for iRun=" << iRun << "\n";
      return 0;
    }
    for (int iCase=0; iCase<4; ++iCase) {
      int choice=(iSelect==0) ? (iBr%10) : -1;
      if ((choice!=-1) && (choice%2 != iCase%2)) continue;
      TString dataEffKind;

      switch(iCase) {
      case 0: dataEffKind="dataRECO_fit-fit"; break;
      case 1: dataEffKind="dataID_fit-fit"; break;
      case 2: dataEffKind="mcRECO_count-count"; break;
      case 3: dataEffKind="mcID_count-count"; break;
      default:
	std::cout << "not ready for iCase=" << iCase << "\n";
	return 0;
      }
      explainChoice=dataEffKind;

      int idx=dataEffKind.Index("_");
      TString resDir=dataEffKind(0,idx) + runTag;

      std::vector<TMatrixD*> effV,effLoV,effHiV;
      std::vector<TString> labelsV;

      if (iSelect) {
	//TString fnameBase="Results-DYee/root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy_runC/efficiency_TnP_1D_Full2012_mcID_count-countEtBins1EtaBins1_PU_varPU";
	TString path= TString("Results-DYee-") + etEtaStr;
	path.Append("/root_files_reg/tag_and_probe/DY_j22_19712pb_Unregressed_energy");
	path.Append(runTag); path.Append("/");
	TString fileNameBase="efficiency_TnP_1D_Full2012_";
	fileNameBase.Append(dataEffKind);
	fileNameBase.Append(etEtaStr);
	fileNameBase.Append("_PU_varPU");

	int weighted=(dataEffKind.Index("count-count")!=-1) ? 1:0;
	if (!loadEffVsPU(path + fileNameBase,
			 weighted,effV,effLoV,effHiV,labelsV)) {
	  std::cout << "failed to load the effs\n";
	  return 0;
	}
	if (0) {
	  for (unsigned int i=0; i<effV.size(); i++) {
	    std::cout << "i= " << i << "  "; effV[i]->Print();
	    std::cout << "\n";
	  }
	}

	fout.mkdir(resDir);
	fout.cd(resDir);
	for (unsigned int i=0; i<effV.size(); ++i) {
	  if (!effV[i] || !effLoV[i] || !effHiV[i]) {
	    HERE("null ptr");
	    return 0;
	  }
	  //HERE(labelsV[i]);
	  effV  [i]->Write(labelsV[i]);
	  effLoV[i]->Write(labelsV[i] + TString("_errLo"));
	  effHiV[i]->Write(labelsV[i] + TString("_errHi"));
	}
	fout.cd();
      }
      else {
	int cnt=DYTools::nPVBinCount;
	labelsV.reserve(cnt);
	effV.reserve(cnt); effLoV.reserve(cnt); effHiV.reserve(cnt);
	for (int pu_i=0; pu_i<DYTools::nPVBinCount; ++pu_i) {
	  TString label=Form("eff_%u_%u",
			     UInt_t(DYTools::nPVLimits[pu_i  ]+0.6),
			     UInt_t(DYTools::nPVLimits[pu_i+1]-0.4));
	  labelsV.push_back(label);
	  TString fieldName=resDir+TString("/") + label;
	  TMatrixD *M = (TMatrixD*)fout.Get(fieldName);
	  TMatrixD *Mlo=(TMatrixD*)fout.Get(fieldName + TString("_errLo"));
	  TMatrixD *Mhi=(TMatrixD*)fout.Get(fieldName + TString("_errHi"));
	  if (!M || !Mlo || !Mhi) {
	    std::cout << "failed to load data for field <" <<
	      fieldName << ">\n";
	    return 0;
	  }
	  effV.push_back(M);
	  effLoV.push_back(Mlo);
	  effHiV.push_back(Mhi);
	}
	TGraphAsymmErrors* gr= effVsPU_asGraph(iEt,iEta,effV,effLoV,effHiV);
	dataV.push_back(gr);
	grLabelsV.push_back(dataEffKind + runTag);
      }
      ClearVec(effV); ClearVec(effLoV); ClearVec(effHiV);
      labelsV.clear();
    }
  }

  if (iSelect) {
    fout.cd();
    writeBinningArrays(fout,"compareEff_vsPU");
  }
  fout.Close();
  HERE("\n\tfile closed\n");
  return 1;
}

// -------------------------------------------------------------------------
