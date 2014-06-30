#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/UnfoldingMatrix.h"
#include "helpers.hh"

int compareFsrUnfValues(int analysisIs2D=0,
			TFsrUnfCompareCase_t theCase=_cmp_fsrGood) {

  if (!DYTools::setup(analysisIs2D)) return retCodeError;

  // ----------------------------------
  // Main part

  TString yAxisLabel="unfolded yield";

  TString path="../../Results-DYee/root_files_reg/";
  TString yieldFName="../../Results-DYee/root_files_reg/xsec/DY_j22_19712pb/xSec_preFsr_1DpostFsrFullSp.root";
  TString yieldField="hpPostFsrFullSp_divLumi";

  TH2D* h2PostFsrCS=LoadHisto2D(yieldField,path+yieldFName,"",1);
  if (!h2PostFsrCS) return retCodeError;

  TH2D* detResp2=LoadHisto2D("DetResponse","dir-Alexey/FSRresMatrixProd.root","",0);
  if (!detResp2) return retCodeError;
  //detResp2->SetTitle("DetFSRResponse");
  detResp2->GetZaxis()->SetRangeUser(0.,1.);

  TMatrixD* UnfM2=createMatrixD(detResp2,0);

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


  if (0) {
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

  if (0) {
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

  TString label1="DYee";
  TString label2="Alexey";

  // ----------------------------------
  // Plot data

  std::vector<TH2D*> hV;
  std::vector<TString> labelV;

  hV.push_back(h2Unf1); labelV.push_back(label1);
  hV.push_back(h2Unf2); labelV.push_back(label2);

  TCanvas *cx= plotProfiles("cx",hV,labelV,NULL,1,
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

  TCanvas *cy= plotProfiles("cy",hErrV,labelV,NULL,1,
			    yAxisLabel + TString(" error"));
  cy->Update();

  return retCodeOk;
}

// ----------------------------------------------------------------
