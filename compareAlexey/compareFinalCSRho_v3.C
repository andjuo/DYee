// An updated file was provided

#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/UnfoldingMatrix.h"
#include "helpers.hh"

int compareFinalCSRho_v3(int analysisIs2D=0, int workCase=0,
			 int stopAt=0, int saveCS=0) {

  if (!DYTools::setup(analysisIs2D)) return retCodeError;

  // ----------------------------------------
  // Main part

  TH2D *h2SigYield_global=NULL;
  TH2D *h2UnfYield1_global=NULL;
  TH2D *h2UnfYield2_global=NULL;

  TH2D *hEffRho1_global=NULL;
  TH2D *hEffRho2_global=NULL;
  TH2D *hAcc1_global=NULL;
  TH2D *hAcc2_global=NULL;

  TMatrixD *fsrInvUnf1_global=NULL;
  TMatrixD *fsrInvUnf2_global=NULL;

  // ----------------------------------------------------------
  //  Unfolding
  // ============================================

  // ----------------------------------------
  // load signal yield and unfold it

  if (1) {
  TString path="../../Results-DYee/root_files_reg/";
  TString yieldFName="yield/DY_j22_19712pb_ApplyEscale/bg-subtracted_yield_1D__peak20140220.root";
  TString yieldField="signalYieldDDbkg";

  TH2D* h2SigYield=LoadHisto2D(yieldField,path+yieldFName,"",1);
  if (!h2SigYield) return retCodeError;

  TH2D* h2SigYield2=NULL;
  if (1) {
    std::cout << "\n\tusing Zee signal\n";
    h2SigYield2=LoadHisto2D("ShapeReweight/hYield_zee",path+yieldFName,"",1);
    if (!h2SigYield2) return retCodeError;
  }
  else {
    h2SigYield2=Clone(h2SigYield,"h2SigYield2");
    std::cout << "\n\tusing the same signal yield\n";
  }


  h2SigYield_global= Clone(h2SigYield,"h2SigYield_ini");
  if (!h2SigYield_global) return retCodeError;

  TMatrixD* UnfM1= loadMatrix(path+TString("constants/DY_j22_19712pb/detResponse_unfolding_constants1D.root"),"DetResponse",41,41,1);
  if (!UnfM1) return retCodeError;
  TH2D* detResp1= createHisto2D(*UnfM1,NULL,
				"detResponse_DYee","detResponse_DYee",
				_colrange_default,1,0.);
  if (!detResp1) return retCodeError;

  TString usePath=TString("Results-DYee-rhoUnf/root_files_reg/");
  usePath=TString("root_files_reg_Rho/");
  TString detResponseName="detResponse";
  switch(workCase) {
  case 0: break;
  case 1: detResponseName="detResponsePostFsrDet"; break;
  default:
    std::cout << "workCase=" << workCase << " is not ready\n";
    return retCodeError;
  }
  TString unfMFileName=usePath + TString(Form("constants/DY_j22_19712pb/%s_unfolding_constants1D.root",detResponseName.Data()));
  TMatrixD* UnfM2= loadMatrix(unfMFileName,"DetResponse",41,41,1);
  if (!UnfM2) return retCodeError;
  TH2D* detResp2= createHisto2D(*UnfM2,NULL,
				"detResponseRho_DYee","detResponseRho_DYee",
				_colrange_default,1,0.);
  if (!detResp2) return retCodeError;

  if (stopAt==1) {
    // compare the response matrices
    // ------------- begin inset
    TH2D* detRespDiff=Clone(detResp1,"detRespDiff");
    detRespDiff->Add(detResp2,-1);
    detRespDiff->SetTitle("difference");

    detResp1->GetZaxis()->SetRangeUser(0,1.);
    detResp2->GetZaxis()->SetRangeUser(0,1.);

    printProfileSums(detResp1);
    printProfileSums(detResp2);
    compareProfileSums(detResp1,detResp2);

    TCanvas *cdiff=new TCanvas("cdiff","cdiff",1200,400);
    cdiff->Divide(3,1);
    for (int i=1; i<3; ++i) {
      TPad *pad=(TPad*)cdiff->GetPad(i);
      pad->SetLogx();
      pad->SetLogy();
    }
    AdjustFor2DplotWithHeight(cdiff);
    cdiff->cd(1);
    detResp1->Draw("COLZ");
    cdiff->cd(2);
    detResp2->Draw("COLZ");
    cdiff->cd(3);
    //detRespDiff->GetZaxis()->SetRangeUser(-0.01,0.);
    detRespDiff->Draw("COLZ");
    //UnfM->Draw("COLZ");
    cdiff->Update();
    return retCodeStop;
    // ------------- end inset
  }

  TMatrixD invUnf1(*UnfM1);
  TMatrixD invUnf2(*UnfM2);

  double det;
  invUnf1.Invert(&det);
  invUnf2.Invert(&det);

  if (stopAt==2) {
    // compare the inverted response matrices
    // ------------- begin inset
    TH2D* detInvResp1= createHisto2D(invUnf1,NULL,
				   "detInvResponse_DYee","detInvResponse_DYee",
				   _colrange_default,0,0.);
    TH2D* detInvResp2= createHisto2D(invUnf2,NULL,
				   "detInvResponse","detInvResponse",
				   _colrange_default,0,0.);
    
    TH2D* detInvRespDiff=Clone(detInvResp1,"detInvRespDiff");
    detInvRespDiff->Add(detInvResp2,-1);
    detInvRespDiff->SetTitle("difference");

    TCanvas *cdiff=new TCanvas("cdiff","cdiff",1200,400);
    cdiff->Divide(3,1);
    for (int i=1; i<3; ++i) {
      TPad *pad=(TPad*)cdiff->GetPad(i);
      pad->SetLogx();
      pad->SetLogy();
    }
    AdjustFor2DplotWithHeight(cdiff);
    cdiff->cd(1);
    detInvResp1->Draw("COLZ");
    cdiff->cd(2);
    detInvResp2->Draw("COLZ");
    cdiff->cd(3);
    detInvRespDiff->Draw("COLZ");
    //UnfM->Draw("COLZ");
    cdiff->Update();
    return retCodeStop;
    // ------------- end inset
  }

  TH2D *h2Unf1=Clone(h2SigYield,"h2Unf_DYee");
  TH2D *h2Unf2=Clone(h2SigYield,"h2Unf2");

  if ( !unfold(h2Unf1, invUnf1, h2SigYield) ||
       !unfold(h2Unf2, invUnf2, h2SigYield2) ) return retCodeError;
  h2UnfYield1_global=Clone(h2Unf1, "h2UnfYield1");
  h2UnfYield2_global=Clone(h2Unf2, "h2UnfYield2");
  if (!h2UnfYield1_global || !h2UnfYield2_global) return retCodeError;

  if (stopAt==3) {
    // ------------ begin inset
    TString label1="DYee";
    TString label2="DYee (#rho corr)";

    std::vector<int> colors;
    colors.push_back(kRed+1);
    colors.push_back(kBlue);

    // ----------------------------------
    // Plot data

    TString yAxisLabel="unfolded yield";
    std::vector<TH2D*> hV;
    std::vector<TString> labelV;

    hV.push_back(h2Unf1); labelV.push_back(label1);
    hV.push_back(h2Unf2); labelV.push_back(label2);

    TCanvas *cx= plotProfiles("cxUnf",hV,labelV,&colors,1,
			      yAxisLabel);
    cx->Update();
      
    // ----------------------------------
    // Plot error

    std::vector<TH2D*> hErrV;

    for (unsigned int i=0; i<hV.size(); ++i) {
      TH2D* hErr=Clone(hV[i],hV[i]->GetName() + TString("_err"));
      swapContentAndError(hErr);
      removeError(hErr);
      hErrV.push_back(hErr);
    }

    TCanvas *cy= plotProfiles("cerrUnf",hErrV,labelV,&colors,1,
			      yAxisLabel + TString(" error"));
    cy->Update();
      // ------------------------ end inset
  }
  }

  // ----------------------------------------------------------
  //  Load corrective factors: EffRho, Acc
  // ============================================

  if (1) {

  for (int iter=0; iter<2; ++iter) {
  TCompareCase_t theCase=(iter==0) ? _cmp_EffRho : _cmp_Acc;

  TH2D* h1=NULL, *h2=NULL, *h3=NULL;
  TH2D *h1Syst=NULL;
  TString fname1,field1,field1err;
  TString fname2,field2,field2err;
  TString fname3,field3,field3err;
  TString label1="DYee", label2="DYeeRho", label3="unknown";

  TString path1="../../Results-DYee/root_files_reg/";
  TString path2=path1;
  TString path3=path1;
  int loadText2=0;
  int is1Dhisto2=0, is1Dhisto3=0;

  TString yAxisLabel="y";

  switch(theCase) {
  case _cmp_RawYield:
    fname1="yield/DY_j22_19712pb_ApplyEscale/bg-subtracted_yield_1D__peak20140220.root";
    field1="Input/observedYield";
    //fname2="raw_yield1D_EE.txt";
    //loadText2=1;
    yAxisLabel="raw yield";
    break;

  case _cmp_FakeBkg:
    fname1="yield/DY_j22_19712pb_ApplyEscale/bg-subtracted_yield_1D__peak20140220.root";
    field1="Input/fakeBackgroundFromData";
    field1err="Input/fakeBackgroundFromDataSyst";
    //fname2="fakeBkg1D_EE.txt";
    //loadText2=1;
    yAxisLabel="fake bkg";
    break;

  case _cmp_TrueBkg:
    fname1="yield/DY_j22_19712pb_ApplyEscale/bg-subtracted_yield_1D__peak20140220.root";
    field1="Input/true2eBackgroundFromData";
    field1err="Input/true2eBackgroundFromDataSyst";
    //fname2="trueBkg1D_EE.txt";
    //loadText2=1;
    yAxisLabel="true bkg";
    break;

  case _cmp_Eff:
    fname1="constants/DY_j22_19712pb/efficiency_1D.root";
    field1="hEfficiency";
    fname2=fname1;
    field2=field2;
    //fname2="efficiencyTotal1D_EE.txt";
    //loadText2=1;
    //fname3="acceff.root";
    //field3="eff_postFSRcorr";
    //is1Dhisto3=1;
    //label3="Alexey (acceff.root)";
    yAxisLabel="efficiency";
    break;

  case _cmp_MCeff:
    fname1="constants/DY_j22_19712pb/efficiency_1D.root";
    field1="hEfficiency";
    //fname2="acceff.root";
    //field2="eff_postFSRcorr";
    //is1Dhisto2=1;
    //label2="Alexey (acceff.root)";
    yAxisLabel="MC efficiency";
    break;

  case _cmp_EffRho:
    fname1="constants/DY_j22_19712pb/efficiency_1D.root";
    field1="hEfficiency";
    fname2=fname1;
    field2=field1;
    //fname2="efficiencyTotal1D_EE.txt";
    //loadText2=1;
    yAxisLabel="efficiency #times #rho";
    break;

  case _cmp_Acc:
    fname1="constants/DY_j22_19712pb/acceptance_1D.root";
    field1="hAcceptance";
    fname2=fname1;
    field2=field1;
    //fname2="acceptance1D_EE.txt";
    //loadText2=1;
    //fname3="acceff.root";
    //field3="acc_postFSRcorr";
    //is1Dhisto3=1;
    //label3="Alexey (acceff.root)";
    yAxisLabel="acceptance";
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
    h2=loadTextFile(path2+fname2,"h2");
    if (!h2) return retCodeError;
  }
  else if (is1Dhisto2) h2= loadHisto1D_convert_TH2D(path2+fname2,field2);
  else h2=LoadHisto2D(field2,path2+fname2,"",1);

  if (fname3.Length()) {
    if (is1Dhisto3) h3= loadHisto1D_convert_TH2D(path3+fname3,field3);
  }

  // ----------------------------------
  // Special adjustments

  if (theCase==_cmp_EffRho) {
    // load the scale factors
    TString rhoCorrFName="../../Results-DYee/root_files_reg/constants/DY_j22_19712pb_egamma_Unregressed_energy/covRhoFileSF_nMB41_asymHLT_Unregressed_energy-allSyst_100_v2.root";
    TH2D* hRho=LoadMatrixFields(rhoCorrFName,1,"scaleFactor","scaleFactorErr",1);
    if (!hRho) return retCodeError;

    int check=0;
    if (check) {
      printHisto(h1);
      printHisto(hRho);
    }
    if (!multiplyHisto(h1,hRho,1)) return retCodeError;
    if (!multiplyHisto(h2,hRho,1)) return retCodeError;
    if (check) printHisto(h1);

    if (0) {
    TH2D *rhoRelErr=getRelError(hRho,"rhoRelErr",0);
    printHisto(rhoRelErr);
    }
  }

  if (theCase==_cmp_EffRho) {
    hEffRho1_global=Clone(h1,"hEffRho_DYee");
    hEffRho2_global=Clone(h2,"hEffRho_DYeeRho");
    if (!hEffRho1_global || !hEffRho2_global) return retCodeError;
  }
  else if (theCase==_cmp_Acc) {
    hAcc1_global=Clone(h1,"hAcc_DYee");
    hAcc2_global=Clone(h2,"hAcc_DYeeRho");
  }

  if (stopAt==4) {
  // ----------------------------------
  // Plot data

  std::vector<TH2D*> hV;
  std::vector<TString> labelV;

  hV.push_back(h1); labelV.push_back(label1);
  hV.push_back(h2); labelV.push_back(label2);
  if (h3) { hV.push_back(h3); labelV.push_back(label3); }

  TString cxName=Form("cxCorr_%d",iter);
  TCanvas *cx= plotProfiles(cxName,hV,labelV,NULL,1,
			    yAxisLabel);
  cx->Update();

  // ----------------------------------
  // Plot error

  std::vector<TH2D*> hErrV;

  for (unsigned int i=0; i<hV.size(); ++i) {
    TH2D* hErr=Clone(hV[i],hV[i]->GetName() + TString("_err"));
    swapContentAndError(hErr);
    removeError(hErr);
    hErrV.push_back(hErr);
  }

  TString cyName=Form("cyCorr_%d",iter);
  TCanvas *cy= plotProfiles(cyName,hErrV,labelV,NULL,1,
			    yAxisLabel + TString(" error"));
  cy->Update();
  }
  }
  } // load corrective factors


  // ----------------------------------------------------------
  //  FSR unfolding
  // ============================================

  if (1) {
    TFsrUnfCompareCase_t theCase=_cmp_fsrGood;

  // ----------------------------------
  // Main part

  TString yAxisLabel="unfolded yield";

  TString path="../../Results-DYee/root_files_reg/";
  TString yieldFName="../../Results-DYee/root_files_reg/xsec/DY_j22_19712pb/xSec_preFsr_1DpostFsrFullSp.root";
  TString yieldField="hpPostFsrFullSp_divLumi";

  TH2D* h2PostFsrCS=LoadHisto2D(yieldField,path+yieldFName,"",1);
  if (!h2PostFsrCS) return retCodeError;

  TString fName=TString("constants/DY_j22_19712pb/detResponse_unfolding_constants1D.root");
  TString tag;
  switch(theCase) {
  case _cmp_fsrGood:
    fName.ReplaceAll("detResponse","fsrGood");
    tag="good";
    break;
  case _cmp_fsrExact:
    fName.ReplaceAll("detResponse","fsrExact");
    tag="exact";
    break;
  default:
    std::cout << "macro is not ready for this case\n";
    return retCodeError;
  }

  TMatrixD* UnfM1= loadMatrix(path+fName,"DetResponse",41,41,1);
  if (!UnfM1) return retCodeError;
  TString histoName1="detFSRResponse_DYee" + tag;
  TH2D* detResp1= createHisto2D(*UnfM1,NULL,histoName1,histoName1,
				_colrange_default,1,0.);
  detResp1->GetZaxis()->SetRangeUser(0,1);
  if (!detResp1) return retCodeError;

  TH2D* detResp2=detResp1;
  TMatrixD* UnfM2=new TMatrixD(*UnfM1);
  if (!UnfM2) return retCodeError;

  if (stopAt==10) {
    // compare the response matrices
    // ------------- begin inset
    TH2D* detRespDiff=Clone(detResp1,"detRespDiff");
    detRespDiff->Add(detResp2,-1);
    detRespDiff->SetTitle("difference");

    TCanvas *cdiff=new TCanvas("cdiff","cdiff",1200,400);
    cdiff->Divide(3,1);
    AdjustFor2DplotWithHeight(cdiff);
    cdiff->cd(1);
    detResp1->Draw("COLZ");
    cdiff->cd(2);
    detResp2->Draw("COLZ");
    cdiff->cd(3);
    detRespDiff->Draw("COLZ");
    //UnfM->Draw("COLZ");
    cdiff->Update();
    printHisto(detRespDiff);
    return retCodeStop;
    // ------------- end inset
  }

  TMatrixD invUnf1(*UnfM1);
  TMatrixD invUnf2(*UnfM2);

  double det;
  invUnf1.Invert(&det);
  invUnf2.Invert(&det);

  if (stopAt==11) {
    // compare the inverted response matrices
    // ------------- begin inset
    TH2D* detInvResp1= createHisto2D(invUnf1,NULL,
				   "detInvResponse_DYee","detInvResponse_DYee",
				   _colrange_default,0,0.);
    TH2D* detInvResp2= createHisto2D(invUnf2,NULL,
				   "detInvResponse","detInvResponse",
				   _colrange_default,0,0.);

    TH2D* detInvRespDiff=Clone(detInvResp1,"detInvRespDiff");
    detInvRespDiff->Add(detInvResp2,-1);
    detInvRespDiff->SetTitle("difference");

    TCanvas *cdiff=new TCanvas("cdiff","cdiff",1200,400);
    cdiff->Divide(3,1);
    AdjustFor2DplotWithHeight(cdiff);
    cdiff->cd(1);
    detInvResp1->Draw("COLZ");
    cdiff->cd(2);
    detInvResp2->Draw("COLZ");
    cdiff->cd(3);
    detInvRespDiff->Draw("COLZ");
    //UnfM->Draw("COLZ");
    cdiff->Update();
    return retCodeStop;
    // ------------- end inset
  }

  TH2D *h2Unf1=Clone(h2PostFsrCS,"h2Unf_DYee");
  TH2D *h2Unf2=Clone(h2PostFsrCS,"h2Unf2");

  if ( !unfold(h2Unf1, invUnf1, h2PostFsrCS) ||
       !unfold(h2Unf2, invUnf2, h2PostFsrCS) ) return retCodeError;

  fsrInvUnf1_global= new TMatrixD(invUnf1);
  fsrInvUnf2_global= new TMatrixD(invUnf2);
  if (!fsrInvUnf1_global || !fsrInvUnf2_global) return retCodeError;

  if (0) { // plot unfolded cross section

  TString label1="DYee";
  TString label2="DYeeRho";

  // ----------------------------------
  // Plot data

  std::vector<TH2D*> hV;
  std::vector<TString> labelV;

  hV.push_back(h2Unf1); labelV.push_back(label1);
  hV.push_back(h2Unf2); labelV.push_back(label2);

  TCanvas *cx= plotProfiles("cxFsr",hV,labelV,NULL,1,
			    yAxisLabel);
  cx->Update();

  // ----------------------------------
  // Plot error

  std::vector<TH2D*> hErrV;

  for (unsigned int i=0; i<hV.size(); ++i) {
    TH2D* hErr=Clone(hV[i],hV[i]->GetName() + TString("_err"));
    swapContentAndError(hErr);
    removeError(hErr);
    hErrV.push_back(hErr);
  }

  TCanvas *cy= plotProfiles("cyFsr",hErrV,labelV,NULL,1,
			    yAxisLabel + TString(" error"));
  cy->Update();
  }
  } // load FSR unfolding

  // ----------------------------------------------------------------
  // ----------------------------------------------------------
  //  Final calculation
  // ============================================

  if (1) {
    TH2D* h2UnfEffRhoYield1= Clone(h2UnfYield1_global,"h2UnfEffRhoYield_DYee");
    TH2D* h2UnfEffRhoYield2= Clone(h2UnfYield2_global,"h2UnfEffRhoYield_DYeeRho");
    int multiply=0; // 1 - multiply, 0 - divide
    if (!multiplyHisto(h2UnfEffRhoYield1, hEffRho1_global, multiply) ||
	!multiplyHisto(h2UnfEffRhoYield2, hEffRho2_global, multiply)) {
      return retCodeError;
    }

    TH2D* h2PostFsrYield1= Clone(h2UnfEffRhoYield1, "h2PostFsrYield_DYee");
    TH2D* h2PostFsrYield2= Clone(h2UnfEffRhoYield2, "h2PostFsrYield_DYeeRho");
    if (!multiplyHisto(h2PostFsrYield1, hAcc1_global, multiply) ||
	!multiplyHisto(h2PostFsrYield2, hAcc2_global, multiply)) {
      return retCodeError;
    }

    TH2D* h2PostFsrCS1= Clone(h2PostFsrYield1, "h2PostFsrCS_DYee");
    TH2D* h2PostFsrCS2= Clone(h2PostFsrYield2, "h2PostFsrCS_DYeeRho");
    h2PostFsrCS1->Scale(1/DYTools::lumiAtECMS);
    h2PostFsrCS2->Scale(1/DYTools::lumiAtECMS);

    TH2D* h2PreFsrCS1= Clone(h2PostFsrCS1, "h2PreFsrCS_DYee");
    TH2D* h2PreFsrCS2= Clone(h2PostFsrCS2, "h2PreFsrCS_DYeeRho");
    if ( !unfold(h2PreFsrCS1, *fsrInvUnf1_global, h2PostFsrCS1) ||
	 !unfold(h2PreFsrCS2, *fsrInvUnf2_global, h2PostFsrCS2) ) {
      return retCodeError;
    }

    if (saveCS) {
      TString fname="cmp_UnfRho-20140604.root";
      TFile fout(fname,"recreate");
      if (!saveHisto(fout,h2PreFsrCS1,"","") ||
	  !saveHisto(fout,h2PreFsrCS2,"","")) {
	return retCodeError;
      }
      writeBinningArrays(fout,"compareAlexey/compareFinalCSRho.C");
      fout.Close();
      std::cout << "file <" << fout.GetName() << "> created\n";
    }

    const int plotSteps[4]= { 1, 1, 1, 1 };

    for (int iCase=0; iCase<4; ++iCase) {
      if (!plotSteps[iCase]) continue;

      TH2D *h1=NULL, *h2=NULL;
      TString canvTitle="cx";
      TString yAxisLabel="y";
      switch(iCase) {
      case 0:
	h1=h2UnfEffRhoYield1;
	h2=h2UnfEffRhoYield2;
	canvTitle.Append("_unfEffRhoYield");
	yAxisLabel="N_{u}/(#epsilon#rho)";
	break;
      case 1:
	h1=h2PostFsrYield1;
	h2=h2PostFsrYield2;
	canvTitle.Append("_postFsrYield");
	yAxisLabel="N_{u}/(A#epsilon#rho)";
	break;
      case 2:
	h1=h2PostFsrCS1;
	h2=h2PostFsrCS2;
	canvTitle.Append("_postFsrCS");
	yAxisLabel="#sigma_{postFsr}=N_{u}/(A#epsilon#rhoL) [pb]";
	break;
      case 3:
	h1=h2PreFsrCS1;
	h2=h2PreFsrCS2;
	canvTitle.Append("_preFsrCS");
	yAxisLabel="#sigma_{preFsr} [pb]";
	break;
      default:
	std::cout << "not ready for iCase=" << iCase << "\n";
	return retCodeError;
      }
      if (!h1 || !h2) {
	std::cout << "histos were not assigned\n";
	return retCodeError;
      }

      TString label1="DYee";
      TString label2="DYeeRho";

      // ----------------------------------
      // Plot data

      std::vector<TH2D*> hV;
      std::vector<TString> labelV;

      hV.push_back(h1); labelV.push_back(label1);
      hV.push_back(h2); labelV.push_back(label2);

      TCanvas *cx= plotProfiles(canvTitle,hV,labelV,NULL,1,
				yAxisLabel);
      cx->Update();

      // ----------------------------------
      // Plot error

      if (0) {
      std::vector<TH2D*> hErrV;

      for (unsigned int i=0; i<hV.size(); ++i) {
	TH2D* hErr=Clone(hV[i],hV[i]->GetName() + TString("_err"));
	swapContentAndError(hErr);
	removeError(hErr);
	hErrV.push_back(hErr);
      }

      TCanvas *cy= plotProfiles(canvTitle+TString("_err"),hErrV,labelV,NULL,1,
				yAxisLabel + TString(" error"));
      cy->Update();
      }
    }
  }

  return retCodeOk;
}
