// An updated file was provided

#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/UnfoldingMatrix.h"
#include "helpers.hh"

int compareUnfValues2(int analysisIs2D=0) {

  if (!DYTools::setup(analysisIs2D)) return retCodeError;

  // ----------------------------------
  // Main part

  TString yAxisLabel="unfolded yield";

  TString path="../../Results-DYee/root_files_reg/";
  TString yieldFName="yield/DY_j22_19712pb_ApplyEscale/bg-subtracted_yield_1D__peak20140220.root";
  TString yieldField="signalYieldDDbkg";

  TH2D* h2SigYield=LoadHisto2D(yieldField,path+yieldFName,"",1);
  if (!h2SigYield) return retCodeError;

  TH2D* detResp2raw=loadXYZTextFile("dir-Alexey-20140602/updated-response-matrix.txt","DetResponse_raw",0,1);
  if (!detResp2raw) return retCodeError;

  TMatrixD* UnfM2=createMatrixD(detResp2raw,0);
  if (!UnfM2) return retCodeError;
  TH2D *detResp2= createHisto2D(*UnfM2,NULL,"detResponse_Alexey","DetResponse_Alexey",_colrange_positive, 1,1.);
  if (!detResp2) return retCodeError;

  if (0) {
    printHisto(detResp2raw);
    UnfM2->Print();
    printHisto(detResp2);
    return retCodeStop;
  }

  TMatrixD* UnfM1= loadMatrix(path+TString("constants/DY_j22_19712pb/detResponse_unfolding_constants1D.root"),"DetResponse",41,41,1);
  if (!UnfM1) return retCodeError;
  TH2D* detResp1= createHisto2D(*UnfM1,NULL,
				"detResponse_DYee","detResponse_DYee",
				_colrange_default,1,0.);
  if (!detResp1) return retCodeError;


  if (0) {
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
       !unfold(h2Unf2, invUnf2, h2SigYield) ) return retCodeError;

  if (1) {
    // ------------ begin inset
    TString label1="DYee";
    TString label2="Alexey";

    std::vector<int> colors;
    colors.push_back(kRed+1);
    colors.push_back(kBlue);

    // ----------------------------------
    // Plot data

    std::vector<TH2D*> hV;
    std::vector<TString> labelV;

    hV.push_back(h2Unf1); labelV.push_back(label1);
    hV.push_back(h2Unf2); labelV.push_back(label2);

    TCanvas *cx= plotProfiles("cx",hV,labelV,&colors,1,
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

    TCanvas *cy= plotProfiles("cy",hErrV,labelV,&colors,1,
			      yAxisLabel + TString(" error"));
    cy->Update();
  }

  return retCodeOk;
}

// ----------------------------------------------------------------
