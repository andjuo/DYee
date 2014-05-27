#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/ComparisonPlot.hh"
#include "CSCovWorkFlags.hh" // for comparison with main result


// ---------------------------------------------------

void adjustESFextraSyst(int analysisIs2D, int iM=0, int iSave=0,
			int compareToMainErr=0) {

  if (!DYTools::setup(analysisIs2D)) return;

  TString path="dir-ESFsyst-20140525/";
  TString fname=path + TString("esf_syst_varEt5_") + DYTools::analysisTag;
  if (analysisIs2D) fname.Append(Form("_iM%d",iM));
  fname.Append(".root");
  std::cout << "working with file=<" << fname << ">\n";

  double transLegendX=0., transLegendY=0.;

  std::vector<TString> fieldNamesV;
  TString refField;

  refField="rho";
  fieldNamesV.push_back("Et6Eta7");
  fieldNamesV.push_back("Et6Eta9");

  TFile file(fname,"read");
  if (!file.IsOpen()) {
    std::cout << "failed to open a file <"
	      << file.GetName() << ">\n";
    return;
  }

  TH1D *hRef=new TH1D(refField,refField,1,0.,1.);
  std::vector<TH1D*> hVar;
  hVar.reserve(fieldNamesV.size());

  // the values in the histograms are relative errors: (val-ref)/ref
  // we need to convert to absolute errors
  int res=loadHisto(file,&hRef,"");
  if (res) {
    std::cout << "rho[2]=" << hRef->GetBinContent(2) << "\n";
    for (unsigned int i=0; res && (i<fieldNamesV.size()); ++i) {
      TH1D* h=new TH1D(fieldNamesV[i],"",1,0.,1.);
      res=loadHisto(file,&h,"");
      std::cout << "relDiff[2]=" << h->GetBinContent(2) << "\n";
      if (res) {
	res=scaleHisto(h,hRef,1); // multiply
	std::cout << "diff[2]=" << h->GetBinContent(2) << ", i.e. (ref*relDiff)\n";
	hVar.push_back(h);
      }
    }
  }
  file.Close();

  std::cout << "File <" << file.GetName() << "> loaded\n";

  TH1D *hErr=(TH1D*)hRef->Clone("hErr");
  hErr->Reset();
  hErr->SetDirectory(0);
  hErr->SetTitle("hErr");

  for (int ibin=1; ibin<=hRef->GetNbinsX(); ++ibin) {
    double dRho=0.;
    for (unsigned int i=0; i<hVar.size(); ++i) {
      double yval=hVar[i]->GetBinContent(ibin);
      double x=fabs(yval);
      if (yval<0) hVar[i]->SetBinContent(ibin, x);
      if (dRho<x) dRho=x;
    }
    hErr->SetBinContent(ibin,dRho);
  }

  // ---------------------
  // Manual adjustment
  // ---------------------

  TH1D *hErrMdf=(TH1D*)hErr->Clone("hErrMdf");
  hErrMdf->SetDirectory(0);
  hErrMdf->SetTitle("hErrMdf");

  if (!analysisIs2D) {
    transLegendY=-0.1;
    hErrMdf->SetBinContent(2, 0.5*(hErrMdf->GetBinContent(1)+hErrMdf->GetBinContent(3)));
    linearApprox(hErrMdf, 66., 77.);
    linearApprox(hErrMdf, 97., 190.);
    linearApprox(hErrMdf, 190., 1999.);
  }
  else {
    transLegendX=-0.35;
    transLegendY=-0.5;

    if (iM==1) {
      linearApprox(hErrMdf, 0.01, 0.2, 0.019, -1.);
      linearApprox(hErrMdf, 1.81, 2.39, -1., -1.);
    }
    if (iM==2) {
      transLegendX=-0.05;
      transLegendY=-0.07;
      //linearApprox(hErrMdf, 0.00, 1.41);
      linearApprox(hErrMdf, 1.81, 2.21);
    }
    if (iM==3) {
      transLegendX=-0.35;
      transLegendY=-0.07;
      //linearApprox(hErrMdf, 0.00, 1.21);
      //linearApprox(hErrMdf, 1.41, 1.91);
      linearApprox(hErrMdf, 2.11, 2.31);
    }
    if (iM==4) {
      transLegendX=-0.35;
      transLegendY=-0.07;
      //linearApprox(hErrMdf, 0.00, 0.81);
      //linearApprox(hErrMdf, 1.01, 1.21);
      //linearApprox(hErrMdf, 1.31, 1.71);
      //linearApprox(hErrMdf, 1.71, 2.01);
      linearApprox(hErrMdf, 2.11, 2.31);
    }
    if (iM==5) {
      transLegendX=-0.35;
      transLegendY=-0.07;
      //linearApprox(hErrMdf, 0.00, 1.31);
      linearApprox(hErrMdf, 1.31, 1.61);
      //linearApprox(hErrMdf, 1.71, 2.01);
      linearApprox(hErrMdf, 2.11, 2.31);
    }
    if (iM==6) {
      transLegendX=-0.35;
      transLegendY=-0.07;
      linearApprox(hErrMdf, 0.00, 0.41, 0.000355);
      //linearApprox(hErrMdf, 1.31, 1.71);
      linearApprox(hErrMdf, 2.11, 2.31, -1.,  0.006624);
    }
  }

  printHisto(hErrMdf);

  TH1D* hMainErr=NULL;
  if (compareToMainErr) {
    HERE("prepare comparison to main error");
    TString fnameBase="../../Results-DYee/root_files_reg/xsec/DY_j22_19712pb/xSec_preFsrSyst_1D.root";
    if (analysisIs2D) fnameBase="../../Results-DYee/root_files_reg/xsec/DY_j22_19712pb/xSec_preFsrDetSyst_2D.root";
    std::vector<TMatrixD*> covs;
    std::vector<TString> labelEffV;
    WorkFlags_t wf;
    wf.editCalcFlags().calc_ESFtot(1);
    wf.finalizeFlags();
    wf.init_ExtraTagV(0);
    wf.extraFileTag(_corrESF,"-esfOnly");
    HERE("load cov matrix");
    if (!loadEsfCovMatrices(fnameBase,covs,labelEffV,wf)) {
      std::cout << "failed to load main err\n";
      return;
    }
    HERE("got cov matrix");
    int idxMin=0;
    int idxMax=hErr->GetNbinsX();
    if (analysisIs2D) {
      idxMin=24*iM;
      idxMax=idxMin + ((iM==6) ? 12 : 24);
    }
    if (idxMax-idxMin != hErr->GetNbinsX()) {
      std::cout << "size mismatch\n";
      return;
    }
    hMainErr=(TH1D*)hErr->Clone("hMainErr");
    hMainErr->Reset();
    hMainErr->SetDirectory(0);

    for (int ibin=1; ibin<=hMainErr->GetNbinsX(); ++ibin) {
      int ii= idxMin+ibin-1;
      hMainErr->SetBinContent( ibin, sqrt((*covs[0])(ii,ii)) );
    }
  }

  TString cpTitle=DYTools::analysisTag;
  if (analysisIs2D) cpTitle.Append(Form(" iM=%d",iM));
  TString xAxisLabel=(analysisIs2D) ? "|y|" : "M_{ee}";
  ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,"cp",cpTitle,
		      xAxisLabel,"#Delta#rho","ratio");
  cp.SetRefIdx(-111);
  if (!analysisIs2D) cp.SetLogx();
  if (!analysisIs2D && hMainErr) cp.SetLogy();
  cp.SetXAxisTextSizes(0.07,0.9, 0.06);
  cp.SetYAxisTextSizes(0.06,1.9, 0.06);

  const int ncolors=3;
  const int colors[ncolors] = { kBlack, kBlue+1, kGreen+1 };
  for (unsigned int i=0; i<fieldNamesV.size(); ++i) {
    hVar[i]->GetXaxis()->SetLabelOffset(0.008);
    removeError(hVar[i]);
    cp.AddHist1D(hVar[i],fieldNamesV[i],"LP",
		 TAttMarker(colors[i%ncolors],20,0.8),i+1,0,1);
  }
  //cp.AddHist1D(hErr,"max","LP",TAttMarker(kOrange,24,1.),0,0,1);
  cp.AddHist1D(hErrMdf,"assigned","LP",TAttMarker(kRed,5,1),1,0,1);

  if (hMainErr) {
    cp.AddHist1D(hMainErr,"err from cov","LP",TAttMarker(38,30,1.),0,0,1);
  }

  TCanvas *cx=new TCanvas("cx","cx",700,700);
  SetSideSpaces(cx,0.05,-0.05,0.,0.02,0);
  cp.Draw(cx,false,"png",0);
  cp.TransLegend(transLegendX,transLegendY);
  cx->Update();

  if (iSave) {
    TString newFName=fname;
    newFName.ReplaceAll(".root","-mdf.root");
    TFile fout(newFName,"recreate");
    for (unsigned int i=0; i<hVar.size(); ++i) {
      hVar[i]->Write(fieldNamesV[i]);
    }
    hErrMdf->Write("rhoSyst");
    writeBinningArrays(fout);
    fout.Close();
    std::cout << "file <" << fout.GetName() << "> saved\n";
  }

  return;
}

// ------------------------------------------------------------------
// ------------------------------------------------------------------

