#include "../Include/DYTools.hh"
#include "../Include/DYToolsUI.hh"
#include "../Include/ComparisonPlot.hh"
#include "../Include/colorPalettes.hh"
#include "EtEtaIndexer.hh"


void plotEtEtaPairs(int idx) {
  TString fname="rhoFileSF_nMB41_RamiUnregEn_100.root";
  TString etEtaString="EtBins6EtaBins5";
  TString mbIdx;

  if (DYTools::study2D==1) {
    fname.ReplaceAll("nMB41","nMB7");
    int im=idx/24;
    int iy=idx-im*24;
    double dy=(im==6) ? 0.2 : 0.1;
    mbIdx=Form("etEtaPairs_Mlo_%1.0lf_ylo_%3.1lf",
	       DYTools::_massBinLimits2D[im], iy*dy);
  }
  else {
    mbIdx=Form("etEtaPairs_M_%1.0lf_%1.0lf",DYTools::_massBinLimits2012[idx],
	       DYTools::_massBinLimits2012[idx+1]);
  }
  std::cout << "idx=" << idx << ", mbIdx=" << mbIdx << "\n";


  TFile fin(fname,"read");
  TString field=Form("etEtaPairs_mfidx_%d",idx);
  TMatrixD *M=(TMatrixD*)fin.Get(field);
  fin.Close();
  if (!M) {
    std::cout << "Failed to load " << mbIdx << "\n";
    return;
  }
  //M->Print();
  
  DYTools::TEtBinSet_t etBinSet=DetermineEtBinSet(etEtaString);
  DYTools::TEtaBinSet_t etaBinSet=DetermineEtaBinSet(etEtaString);
  int nEtBins= DYTools::getNEtBins(etBinSet);
  int nEtaBins= DYTools::getNEtaBins(etaBinSet);
  double *loc_etBinLimits=DYTools::getEtBinLimits(etBinSet);
  double *loc_etaBinLimits=DYTools::getEtaBinLimits(etaBinSet);

  //  TMatrixD 
  TH2D *h2=new TH2D(Form("h_%s",mbIdx.Data()),mbIdx,
		    nEtaBins,loc_etaBinLimits,
		    nEtBins, loc_etBinLimits);
  h2->GetYaxis()->SetNoExponent();
  h2->GetYaxis()->SetMoreLogLabels();
  h2->GetYaxis()->SetTitle("#it{E}_{T}");
  h2->GetXaxis()->SetTitle("|#it{#eta}|");

  TH2D *h2lead=Clone(h2,Form("h2Lead_%s",mbIdx.Data()),Form("lead %s",mbIdx.Data()));
  TH2D *h2trail=Clone(h2,Form("h2Trail_%s",mbIdx.Data()),Form("trail %s",mbIdx.Data()));

  EtEtaIndexer_t fi1(nEtBins,nEtaBins);
  EtEtaIndexer_t fi2(nEtBins,nEtaBins);

  for (int ir=0; ir<M->GetNrows(); ir++) {
    fi1.setEtEtaBin(ir);
    double et1=loc_etBinLimits[fi1.getEtBin()]+1e-2;
    double eta1=loc_etaBinLimits[fi1.getEtaBin()]+1e-2;
    std::cout << "ir=" <<ir << ", et1=" << et1 << ", eta1=" << eta1 << "\n";
    for (int ic=ir; ic<M->GetNcols(); ++ic) {
      fi2.setEtEtaBin(ic);
      double et2=loc_etBinLimits[fi2.getEtBin()]+1e-2;
      double eta2=loc_etaBinLimits[fi2.getEtaBin()]+1e-2;
      h2->Fill(eta1,et1, (*M)(ir,ic));
      h2->Fill(eta2,et2, (*M)(ir,ic));
      h2lead->Fill(eta2,et2, (*M)(ir,ic));
      h2trail->Fill(eta1,et1, (*M)(ir,ic));
    }
  }

  set_nice_style(21);

  TCanvas *cx=new TCanvas("cx","cx",700,700);
  cx->SetLogy(1);
  AdjustFor2DplotWithHeight(cx);
  h2->Draw("COLZ");
  cx->Update();

  TCanvas *cx1=new TCanvas("cx1","cx1",700,700);
  cx1->SetLogy(1);
  AdjustFor2DplotWithHeight(cx1);
  h2lead->Draw("COLZ");
  cx1->Update();

  TCanvas *cx2=new TCanvas("cx2","cx2",700,700);
  cx2->SetLogy(1);
  AdjustFor2DplotWithHeight(cx2);
  h2trail->Draw("COLZ");
  cx2->Update();

  eliminateSeparationSigns(mbIdx,1);
  TString cfname =TString("fig-sum-") + mbIdx;
  TString cfname1=TString("fig-lead-") + mbIdx;
  TString cfname2=TString("fig-trail-") + mbIdx;
  TString path=TString("dir-pairPlot-") + DYTools::analysisTag;
  SaveCanvas(cx ,cfname ,path);
  SaveCanvas(cx1,cfname1,path);
  SaveCanvas(cx2,cfname2,path);

}
