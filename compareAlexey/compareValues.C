#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "helpers.hh"

int compareValues(int analysisIs2D, TCompareCase_t theCase) {

  if (!DYTools::setup(analysisIs2D)) return retCodeError;

  // ----------------------------------
  // Main part

  TH2D* h1=NULL, *h2=NULL, *h3=NULL;
  TH2D *h1Syst=NULL;
  TString fname1,field1,field1err;
  TString fname2,field2,field2err;
  TString fname3,field3,field3err;
  TString label1="DYee", label2="Alexey", label3="Alexey(v2)";

  TString path1="../../Results-DYee/root_files_reg/";
  TString path2="dir-Alexey/";
  TString path3="dir-Alexey/";
  int loadText2=0;
  int is1Dhisto2=0, is1Dhisto3=0;

  TString yAxisLabel="y";

  switch(theCase) {
  case _cmp_RawYield:
    fname1="yield/DY_j22_19712pb_ApplyEscale/bg-subtracted_yield_1D__peak20140220.root";
    field1="Input/observedYield";
    fname2="raw_yield1D_EE.txt";
    loadText2=1;
    yAxisLabel="raw yield";
    break;

  case _cmp_FakeBkg:
    fname1="yield/DY_j22_19712pb_ApplyEscale/bg-subtracted_yield_1D__peak20140220.root";
    field1="Input/fakeBackgroundFromData";
    field1err="Input/fakeBackgroundFromDataSyst";
    fname2="fakeBkg1D_EE.txt";
    loadText2=1;
    yAxisLabel="fake bkg";
    break;

  case _cmp_TrueBkg:
    fname1="yield/DY_j22_19712pb_ApplyEscale/bg-subtracted_yield_1D__peak20140220.root";
    field1="Input/true2eBackgroundFromData";
    field1err="Input/true2eBackgroundFromDataSyst";
    fname2="trueBkg1D_EE.txt";
    loadText2=1;
    yAxisLabel="true bkg";
    break;

  case _cmp_Eff:
    fname1="constants/DY_j22_19712pb/efficiency_1D.root";
    field1="hEfficiency";
    fname2="efficiencyTotal1D_EE.txt";
    loadText2=1;
    fname3="acceff.root";
    field3="eff_postFSRcorr";
    is1Dhisto3=1;
    label3="Alexey (acceff.root)";
    yAxisLabel="efficiency";
    break;

  case _cmp_MCeff:
    fname1="constants/DY_j22_19712pb/efficiency_1D.root";
    field1="hEfficiency";
    fname2="acceff.root";
    field2="eff_postFSRcorr";
    is1Dhisto2=1;
    label2="Alexey (acceff.root)";
    yAxisLabel="MC efficiency";
    break;

  case _cmp_EffRho:
    fname1="constants/DY_j22_19712pb/efficiency_1D.root";
    field1="hEfficiency";
    fname2="efficiencyTotal1D_EE.txt";
    loadText2=1;
    yAxisLabel="efficiency #times #rho";
    break;

  case _cmp_Acc:
    fname1="constants/DY_j22_19712pb/acceptance_1D.root";
    field1="hAcceptance";
    fname2="acceptance1D_EE.txt";
    loadText2=1;
    fname3="acceff.root";
    field3="acc_postFSRcorr";
    is1Dhisto3=1;
    label3="Alexey (acceff.root)";
    yAxisLabel="acceptance";
    break;

  case _cmp_UnfYields:
    path1="./";
    fname1="cmp_Alexey_unfYields-20140604.root";
    field1="unfYields_DYee";
    path2="./";
    fname2="dir-Alexey-20140603-unf/unfData.txt";
    loadText2=3;
    if (1) {
      path3="./";
      fname3=fname1;
      field3="unfYields_Alexey";
      label3="Alexey's unfMatrix";
    }
    yAxisLabel="unfolded yield";
    break;

  default:
    std::cout << "Not ready for the case\n";
    return retCodeError;
  }

  // ----------------------------------
  // Load data

  h1=LoadHisto2D(field1,path1+fname1,"",1);
  if (!h1) return retCodeError;
  if (field1err.Length()) {
    h1Syst=LoadHisto2D(field1err,path1+fname1,"",1);
    if (!h1Syst) return retCodeError;
    h1->Add(h1Syst,1.);
  }
  
  if (loadText2) {
    int hasDivider=((loadText2==1) || (loadText2!=3)) ? 1:0;
    int hasError=  ((loadText2==1) || (loadText2!=3)) ? 1:0;
    h2=loadTextFile(path2+fname2,"h2",hasDivider,hasError);
    if (!h2) return retCodeError;

    if (0 && (loadText2==3)) {
      std::cout << "check the results (loadText2==3)\n";
      printHisto(h2);
      return retCodeStop;
    }
  }
  else if (is1Dhisto2) h2= loadHisto1D_convert_TH2D(path2+fname2,field2);

  if (fname3.Length()) {
    if (is1Dhisto3) h3= loadHisto1D_convert_TH2D(path3+fname3,field3);
    else h3=LoadHisto2D(field3,path3+fname3,"",1);
  }

  // ----------------------------------
  // Special adjustments

  TH2D* rhoErrChk=NULL;
  TH2D* effRhoErrChk=NULL;

  if (theCase==_cmp_EffRho) {
    // load the scale factors
    TString rhoCorrFName="../../Results-DYee/root_files_reg/constants/DY_j22_19712pb_egamma_Unregressed_energy/scale_factors_asymHLT_1D.root";
    rhoCorrFName="../../Results-DYee/root_files_reg/constants/DY_j22_19712pb_egamma_Unregressed_energy/covRhoFileSF_nMB41_asymHLT_Unregressed_energy-allSyst_100_v2.root";
    TH2D* hRho=LoadMatrixFields(rhoCorrFName,1,"scaleFactor","scaleFactorErr",1);
    if (!hRho) return retCodeError;
    int check=1;
    if (check) {
      printHisto(h1);
      printHisto(hRho);
    }
    TH2D *hEff_loc=Clone(h1,"hEff_loc");
    if (!multiplyHisto(h1,hRho,1)) return retCodeError;
    h1->SetName("hEffRho");
    if (check) printHisto(h1);

    TH2D* rhoRelErr_loc=getRelError(hRho,"rhoRelErr_loc",0);
    printHisto(rhoRelErr_loc);

    if (0) {
      rhoErrChk=Clone(hRho,"rhoErrChk");
      swapContentAndError(rhoErrChk);
      removeError(rhoErrChk);
    }

    if (1) { // check if the rror is added in quadrature
      effRhoErrChk=Clone(hRho,"effRhoErrChk");
      for (int ibin=1; ibin<=h1->GetNbinsX(); ++ibin) {
	for (int jbin=1; jbin<=h1->GetNbinsY(); ++jbin) {
	  double e1= hRho->GetBinError(ibin,jbin);
	  double e2= hEff_loc->GetBinError(ibin,jbin);
	  effRhoErrChk->SetBinContent(ibin,jbin,sqrt(e1*e1+e2*e2));
	}
      }
      removeError(effRhoErrChk);
    }
  }

  // ----------------------------------
  // Plot data

  std::vector<TH2D*> hV;
  std::vector<TString> labelV;

  hV.push_back(h1); labelV.push_back(label1);
  hV.push_back(h2); labelV.push_back(label2);
  if (h3) { hV.push_back(h3); labelV.push_back(label3); }

  TCanvas *cx= plotProfiles("cx",hV,labelV,NULL,1,
			    yAxisLabel);
  cx->Update();

  // ----------------------------------
  // Plot error

  if (theCase!=_cmp_UnfYields) {

  std::vector<TH2D*> hErrV;

  for (unsigned int i=0; i<hV.size(); ++i) {
    TH2D* hErr=Clone(hV[i],hV[i]->GetName() + TString("_err"));
    swapContentAndError(hErr);
    removeError(hErr);
    hErrV.push_back(hErr);
  }
  if (rhoErrChk) {
    hErrV.push_back(rhoErrChk);
    labelV.push_back("#rho errors");
  }
  if (effRhoErrChk) {
    hErrV.push_back(effRhoErrChk);
    labelV.push_back("errors added in quadrature");
  }

  TCanvas *cy= plotProfiles("cy",hErrV,labelV,NULL,1,
			    yAxisLabel + TString(" error"));
  cy->Update();
  }

  return retCodeOk;
}

// ----------------------------------------------------------------
