#include <TROOT.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include "../Include/MitStyleRemix.hh"
#include "../Include/ComparisonPlot.hh"
#include "../Include/MyTools.hh"
#include "../Include/UnfoldingTools.hh"
#include "calcEventEffLink.h"
#include "../Covariance/colorPalettes.hh"

// --------------------------------------------------------------

TH1F* getHistoFromFile(const char *fname, 
		       const char *varName, const char *varErrName,
		       int isMatrix,
		       const char *newHistoName) {
   TFile miofile(fname,"read");
   const int nMBins=DYTools::nMassBins;
   const int nYBins=DYTools::nYBinsMax;
   TMatrixD var(nMBins,nYBins), varErr(nMBins,nYBins);

   if (isMatrix) {
     TMatrixD *x=(TMatrixD*) miofile.Get(varName);
     TMatrixD *xErr=(TMatrixD*) miofile.Get(varErrName);
     if (!x || !xErr) {
       std::cout << "failed to get <" << varName << "> or <" << varErrName << ">\n";
       return NULL;
     }
     var= *x; 
     varErr= *xErr;
     delete x;
     delete xErr;
   }
   else {
     TVectorD *x=(TVectorD*) miofile.Get(varName);
     TVectorD *xErr=(TVectorD*) miofile.Get(varErrName);
     if (!x || !xErr) {
       std::cout << "failed to get <" << varName << "> or <" << varErrName << ">\n";
       return NULL;
     }
     unfolding::deflattenMatrix(*x, var);
     unfolding::deflattenMatrix(*xErr, varErr);
     delete x;
     delete xErr;
   }

   TH1F* histo=
     extractMassDependence(newHistoName,"",
			   var,varErr,
			   0,0,0);
   miofile.Close();
   return histo;
}


// --------------------------------------------------------------

TH1F* getHistoFromFile(const std::string &fname,
		       const std::string &varName,
		       const std::string &varErrName,
		       int isMatrix,
		       const std::string &newHistoName) {
  return getHistoFromFile(fname.c_str(), 
			  varName.c_str(), 
			  varErrName.c_str(),
			  isMatrix,
			  newHistoName.c_str());
}

// --------------------------------------------------------------

void calcHltEffTrunc(DYTools::TDataKind_t dataKind,
		 int iEt1, int iEta1,
		 int iEt2, int iEta2,
		 double &hlt_eff, double &hlt_effErr) {
  const vector<TMatrixD*> *eff=(dataKind==DYTools::DATA) ? &dataEff : &mcEff;
  const vector<TMatrixD*> *effErr=(dataKind==DYTools::DATA) ? &dataEffAvgErr : &mcEffAvgErr;

  double eff1=(*(*eff)[DYTools::HLT])[iEt1][iEta1];
  double eff2=(*(*eff)[DYTools::HLT])[iEt2][iEta2];
  double eff1err=(*(*effErr)[DYTools::HLT])[iEt1][iEta1];
  double eff2err=(*(*effErr)[DYTools::HLT])[iEt2][iEta2];
  hlt_eff= eff1*eff2;
  hlt_effErr=sqrt( pow(eff1err*eff2,2) + pow(eff1*eff2err,2));
  return;
}

// --------------------------------------------------------------

void calcHltEff(DYTools::TDataKind_t dataKind,
		int iEtBin1, int iEtaBin1, double et1,
		int iEtBin2, int iEtaBin2, double et2,
		double &hlt_eff, double &hlt_effErr) {

  hlt_eff=getHLTefficiency(dataKind,
			   iEtBin1,iEtaBin1, et1,
			   iEtBin2,iEtaBin2, et2);
  hlt_effErr=getHLTefficiencyErr(dataKind,
				 iEtBin1,iEtaBin1, et1,
				 iEtBin2,iEtaBin2, et2);
  //if ((dataKind==DYTools::MC) && (et1>100) && (et2>100)) {
  //  std::cout << "MC HLT eff et1,et2>100 " << hlt_eff << " +/- " << hlt_effErr << "\n";
  //}
  return;
}


// --------------------------------------------------------------



// --------------------------------------------------------------

void compareDoubleHLTeff() {

  const TString mcInputFile="../config_files/fall11mc_vilnius.input";
  const TString tnpMCInputFile="../config_files/sf_mc_et6_eta5_vilnius.conf";
  const TString tnpDataInputFile="../config_files/sf_data_et6_eta5_vilnius.conf";
  const TString triggerSetString="Full2011_hltEffOld";
  int puReweight=1;

  if (nonUniversalHLT==0) {
    std::cout << "nonUniversalHLT has to be 1\n";
    return;
  }

  // Construct the trigger object
  TriggerSelection triggers(triggerSetString, false, 0); 
  assert ( triggers.isDefined() );

  MCInputFileMgr_t mcMgr;
  TnPInputFileMgr_t tnpDataMgr,tnpMCMgr;

  if (!mcMgr.Load(mcInputFile) ||
      !tnpMCMgr.Load(tnpMCInputFile) ||
      !tnpDataMgr.Load(tnpDataInputFile)) {
    return;
  }
  if (!tnpMCMgr.hasSameBinCounts(tnpDataMgr)) {
    cout << "Files tnpMCInputFile=<" << tnpMCInputFile 
	 << ">, tnpDataInputFile=<" << tnpDataInputFile 
	 << "> have different bin counts:\n";
    cout << "MC   input: " << tnpMCMgr;
    cout << "Data input: " << tnpDataMgr;
    return;
  }
  dirTag=tnpMCMgr.dirTag();
 
  etBinning=tnpMCMgr.etBinsKind();
  etBinCount=DYTools::getNEtBins(etBinning);
  etBinLimits=DYTools::getEtBinLimits(etBinning);
  
  etaBinning=tnpMCMgr.etaBinsKind();
  etaBinCount=DYTools::getNEtaBins(etaBinning);
  etaBinLimits=DYTools::getEtaBinLimits(etaBinning);

  if (!fillEfficiencyConstants( tnpMCMgr, tnpDataMgr, triggers, puReweight )) {
    return;
  }

  std::vector<TGraphErrors*> grDataV, grMCV;
  std::vector<TString> etaStrV, plotLabelV;
  std::vector<TH2D*> effData2DV, effMC2DV;

  if (1) {
    double etTrunc[etBinCount+1], dEtTrunc[etBinCount+1];
    double etBinLimitsExt[etBinCount+2];
    double etExt[etBinCount+2], dEtExt[etBinCount+2];
    int iEtExtra=2;
    for (int i=0; i<etBinCount+1; ++i) {
      etTrunc[i] = 0.5*(etBinLimits[i] + etBinLimits[i+1]);
      dEtTrunc[i]= 0.5*(etBinLimits[i+1] - etBinLimits[i]);
    }
    for (int i=0, idx=0; i<etBinCount+1; ++i,++idx) {
      etBinLimitsExt[idx]=etBinLimits[i];
      if (idx==0) etBinLimitsExt[idx]-=1e-3;
      if ((etBinLimits[i]<17.) && (etBinLimits[i+1]>17.)) {
	etBinLimitsExt[idx+1]=17.;
	iEtExtra=idx+1;
	std::cout << "setting iEtExtra=" << idx+1 << "\n";
	idx++;
      }
      std::cout << "i=" << i << ", idx=" << idx << std::endl;
    }
    for (int i=0; i<etBinCount+2; ++i) {
      etExt[i] = 0.5*(etBinLimitsExt[i] + etBinLimitsExt[i+1]);
      dEtExt[i]= 0.5*(etBinLimitsExt[i+1] - etBinLimitsExt[i]);
    }
    if (1) {
      std::cout << "etBinLimitsExt: ";
      for (int i=0; i<etBinCount+2; ++i) {
	std::cout << " " << etBinLimitsExt[i];
      }
      std::cout << "\n";
    }
    
    double effDataTrunc[etBinCount], effDataTruncErr[etBinCount];
    double effMCTrunc[etBinCount], effMCTruncErr[etBinCount];
    double effData[etBinCount+1], effDataErr[etBinCount+1];
    double effMC[etBinCount+1], effMCErr[etBinCount+1];
    TString plotLabel;
    
    TH2D *effData2Dtrunc=new TH2D("effDataTrunc","effDataTrunc",
				  etBinCount+1,etBinLimitsExt,
				  etBinCount+1,etBinLimitsExt);
    TH2D *effMC2Dtrunc=new TH2D("effMCTrunc","effMCTrunc",
				etBinCount+1,etBinLimitsExt,
				etBinCount+1,etBinLimitsExt);
    TH2D *effData2D=new TH2D("effData","effData",
			     etBinCount+1,etBinLimitsExt,
			     etBinCount+1,etBinLimitsExt);
    TH2D *effMC2D=new TH2D("effMC","effMC",
			   etBinCount+1,etBinLimitsExt,
			   etBinCount+1,etBinLimitsExt);
    effData2Dtrunc->SetDirectory(0);
    effMC2Dtrunc->SetDirectory(0);
    effData2D->SetDirectory(0);
    effMC2D->SetDirectory(0);
    
    
    for (int iEta1=0; iEta1<etaBinCount; ++iEta1) {
      for (int iEta2=iEta1; iEta2<etaBinCount; ++iEta2) {
	// skip the gap
	if ((etaBinning==DYTools::ETABINS5) && ((iEta1==2) || (iEta2==2))) continue;

	TString etaStr=Form("_abs_eta_%1.1lf_%1.1lf__%1.1lf_%1.1lf",
			    etaBinLimits[iEta1],etaBinLimits[iEta1+1],
			    etaBinLimits[iEta2],etaBinLimits[iEta2+1]);
	etaStrV.push_back(etaStr);
	plotLabel=Form("#splitline{%3.1lf < |#eta_{1}| < %3.1lf}{%3.1lf < |#eta_{2}| < %3.1lf}",
		       etaBinLimits[iEta1],etaBinLimits[iEta1+1],
		       etaBinLimits[iEta2],etaBinLimits[iEta2+1]);
	plotLabelV.push_back(plotLabel);
	
	effData2Dtrunc->Reset();
	effMC2Dtrunc->Reset();
	effData2D->Reset();
	effMC2D->Reset();

	double eff,effErr;
	for (int iEt=0; iEt<etBinCount; ++iEt) {
	  calcHltEffTrunc(DYTools::DATA,iEt,iEta1,iEt,iEta2, eff,effErr);
	  effDataTrunc[iEt]=eff;
	  effDataTruncErr[iEt]=effErr;
	  calcHltEffTrunc(DYTools::MC,iEt,iEta1,iEt,iEta2, eff,effErr);
	  effMCTrunc[iEt]=eff;
	  effMCTruncErr[iEt]=effErr;
	  
	  int idxEt=(iEt>=iEtExtra) ? iEt+1 : iEt;
	  double effMCVal,effMCerrVal;
	  int print=1;
	  for (int iEt2=0; iEt2<etBinCount; ++iEt2) {
	    int idxEt2=(iEt2>=iEtExtra) ? iEt2+1 : iEt2;
	    calcHltEffTrunc(DYTools::DATA,iEt,iEta1,iEt2,iEta2, eff,effErr);
	    effData2Dtrunc->SetBinContent(idxEt+1,idxEt2+1, eff);
	    effData2Dtrunc->SetBinError(idxEt+1,idxEt2+1, effErr);
	    calcHltEffTrunc(DYTools::MC,iEt,iEta1,iEt2,iEta2, effMCVal,effMCerrVal);
	    effMC2Dtrunc->SetBinContent(idxEt+1,idxEt2+1, effMCVal);
	    effMC2Dtrunc->SetBinError(idxEt+1,idxEt2+1, effMCerrVal);
	    if (print) std::cout << "effMC   (" << idxEt+1 << "," << idxEt2+1 << ")=" << effMCVal << " +/- " << effMCerrVal << "\n";
	    if ((iEt==iEtExtra-1) && (iEt2==iEtExtra-1)) {
	      effData2Dtrunc->SetBinContent(iEt+2,iEt2+2, eff);
	      effData2Dtrunc->SetBinError(iEt+2,iEt2+2, effErr);
	      effMC2Dtrunc->SetBinContent(iEt+2,iEt2+2, effMCVal);
	      effMC2Dtrunc->SetBinError(iEt+2,iEt2+2, effMCerrVal);
	      if (print) std::cout << "effMC*a (" << iEt+2 << "," << iEt2+2 << ")=" << effMCVal << " +/- " << effMCerrVal << "\n";
	    }
	    if (iEt==iEtExtra-1) {
	      effData2Dtrunc->SetBinContent(iEt+2,idxEt2+1, eff);
	      effData2Dtrunc->SetBinError(iEt+2,idxEt2+1, effErr);
	      effMC2Dtrunc->SetBinContent(iEt+2,idxEt2+1, effMCVal);
	      effMC2Dtrunc->SetBinError(iEt+2,idxEt2+1, effMCerrVal);
	      if (print) std::cout << "effMC*b (" << iEt+2 << "," << idxEt2+1 << ")=" << effMCVal << " +/- " << effMCerrVal << "\n";
	    }
	    if (iEt2==iEtExtra-1) {
	      effData2Dtrunc->SetBinContent(idxEt+1,iEt2+2, eff);
	      effData2Dtrunc->SetBinError(idxEt+1,iEt2+2, effErr);
	      effMC2Dtrunc->SetBinContent(idxEt+1,iEt2+2, effMCVal);
	      effMC2Dtrunc->SetBinError(idxEt+1,iEt2+2, effMCerrVal);
	      if (print) std::cout << "effMC*c (" << idxEt+1 << "," << iEt2+2 << ")=" << effMCVal << " +/- " << effMCerrVal << "\n";
	    }
	  }
	  
	}
	
	for (int iEt=0; iEt<etBinCount+1; ++iEt) {
	  double etCenter=0.5*(etBinLimitsExt[iEt] + etBinLimitsExt[iEt+1]);
	  int idxEt=(iEt>=iEtExtra) ? iEt-1 : iEt;
	  if ((iEta1==0) && (iEta2==0)) {
	    std::cout << "iEt=" << iEt << ", idxEt=" << idxEt << ", etCenter=" << etCenter << "\n";
	  }
	  calcHltEff(DYTools::DATA,
		     idxEt,iEta1,etCenter,
		     idxEt,iEta2,etCenter,  eff,effErr);
	  effData[iEt]=eff;
	  effDataErr[iEt]=effErr;
	  std::cout << " data eff iEt=" << iEt  << ", iEta1=" << iEta1 << ", iEta2=" << iEta2 << ", eff=" << eff << " +/- " << effErr << "\n";
	  calcHltEff(DYTools::MC,
		     idxEt,iEta1,etCenter,
		     idxEt,iEta2,etCenter,  eff,effErr);
	  effMC[iEt]=eff;
	  effMCErr[iEt]=effErr;
	  //std::cout << "MC eff iEt=" << iEt << ", eff=" << eff << " +/- " << effErr << std::endl;
	  
	  for (int iEt2=0; iEt2<etBinCount+1; ++iEt2) {
	    //std::cout << "iEt2=" << iEt2 << ", iEtExtra=" << iEtExtra << "\n";
	    double etCenter2=0.5*(etBinLimitsExt[iEt2] + etBinLimitsExt[iEt2+1]);
	    
	    int idxEt2=(iEt2>=iEtExtra) ? (iEt2-1) : iEt2;
	    if ((iEta1==0) && (iEta2==0)) {
	      std::cout << "iEt2=" << iEt2 << ", idxEt2=" << idxEt2 << ", etCenter2=" << etCenter2 << " (iEtExtra=" << iEtExtra << ")\n";
	    }
	    calcHltEff(DYTools::DATA,
		       idxEt ,iEta1,etCenter,
		       idxEt2,iEta2,etCenter2,  eff,effErr);
	    effData2D->SetBinContent(iEt+1,iEt2+1, eff);
	    effData2D->SetBinError(iEt+1,iEt2+1, effErr);
	    calcHltEff(DYTools::MC,
		       idxEt ,iEta1,etCenter,
		       idxEt2,iEta2,etCenter2,  eff,effErr);
	    effMC2D->SetBinContent(iEt+1,iEt2+1, eff);
	    effMC2D->SetBinError(iEt+1,iEt2+1, effErr);
	    //HERE("xx");
	  }
	}
	HERE("creating graphs");

	TGraphErrors *grDataEffTrunc=
	  new TGraphErrors(etBinCount, etTrunc, effDataTrunc, dEtTrunc, effDataTruncErr);
	TGraphErrors *grMCEffTrunc=
	  new TGraphErrors(etBinCount, etTrunc, effMCTrunc, dEtTrunc, effMCTruncErr);
	TGraphErrors *grDataEff=
	  new TGraphErrors(etBinCount+1, etExt, effData, dEtExt, effDataErr);
	TGraphErrors *grMCEff=
	  new TGraphErrors(etBinCount+1, etExt, effMC, dEtExt, effMCErr);
	grDataV.push_back(grDataEffTrunc);
	grMCV.push_back(grMCEffTrunc);
	grDataV.push_back(grDataEff);
	grMCV.push_back(grMCEff);

	TH2D *h2=NULL;
	h2=(TH2D*)effData2Dtrunc->Clone(TString("effDataTrunc_")+etaStr);
	h2->SetDirectory(0);
	effData2DV.push_back(h2);
	h2=(TH2D*)effData2D->Clone(TString("effData_")+etaStr);
	h2->SetDirectory(0);
	effData2DV.push_back(h2);
	h2=(TH2D*)effMC2Dtrunc->Clone(TString("effMCTrunc_")+etaStr);
	h2->SetDirectory(0);
	effMC2DV.push_back(h2);
	h2=(TH2D*)effMC2D->Clone(TString("effMC_")+etaStr);
	h2->SetDirectory(0);
	effMC2DV.push_back(h2);
      }
    }
  }

  // make plots
  if (1) { // Et-diagonal efficiencies
    for (unsigned int i=0; i<grDataV.size(); i+=2) {
      if (i!=4) continue;
      unsigned int idx=i/2;
      TString canvName=TString("canv_") + etaStrV[idx];
      TString cpName=TString("cp_") + etaStrV[idx];
      TCanvas *c=new TCanvas(canvName,canvName, 600,600);
      c->SetLeftMargin(0.2);
      CPlot *plot1=new CPlot(cpName,"","E_{T} [GeV]","e_{1}e_{2} HLT efficiency");

      plot1->SetYRange(0.6,1.1);
      
      plot1->SetLogx();
      grDataV[i]->GetXaxis()->SetMoreLogLabels();
      grDataV[i]->GetXaxis()->SetNoExponent();
      grDataV[i]->SetFillStyle(3005);
      grDataV[i]->SetFillColor(kBlue);
      grDataV[i+1]->SetFillStyle(3004);
      grDataV[i+1]->SetFillColor(kRed+1);
      plot1->AddGraph(grDataV[i],"data eff (trunc)","PE2",kBlue);
      plot1->AddGraph(grMCV[i],"MC eff (trunc)","PE",kBlack,24);
      plot1->AddGraph(grDataV[i+1],"data eff (full)","PE2",kRed+1);
      plot1->AddGraph(grMCV[i+1],"MC eff (full)","PE",kGreen+2,24);

      plot1->AddTextBox(plotLabelV[idx], 0.6,0.45,0.87,0.6, 0);
      plot1->TransLegend(0.02,-0.5);
      
      plot1->Draw(c);
      
      TLine *line = new TLine(10,1.0,500,1.0);
      line->SetLineStyle(kDashed);
      line->Draw("same");
      c->Update();
      //break;
    }
  }


  if (1) { // 2D distributions
    int studyData=0;
    for (unsigned int i=0; i<effData2DV.size(); i+=2) {
      if (i!=2) continue;
      TH2D *h2trunc=(studyData) ? effData2DV[i  ] : effMC2DV[i  ];
      h2trunc->GetXaxis()->SetTitle("E_{T;1}");
      h2trunc->GetYaxis()->SetTitle("E_{T;2}");
      TH2D *h2full =(studyData) ? effData2DV[i+1] : effMC2DV[i+1];
      h2full->GetXaxis()->SetTitle("E_{T;1}");
      h2full->GetYaxis()->SetTitle("E_{T;2}");
      TString diffName=h2trunc->GetTitle() + TString("_minus_") + h2full->GetTitle();
      TH2D *h2Diff=(TH2D*)h2trunc->Clone(diffName);
      h2Diff->Add(h2full,-1.);
      h2Diff->SetTitle(diffName);
      unsigned int idx=i/2;
      TString canvName=TString("canv_") + etaStrV[idx];
      TString cpName=TString("cp_") + etaStrV[idx];
      TCanvas *c=new TCanvas(canvName,canvName, 1200,400);
      c->Divide(3,1);
      AdjustFor2DplotWithHeight(c,0.12);
      for (int ipad=1; ipad<4; ++ipad) c->GetPad(ipad)->SetLeftMargin(0.17);
      //for (int ipad=1; ipad<4; ++ipad) c->GetPad(ipad)->SetRightMargin(0.13);
      //c->cd(1);
      //h2trunc->Draw("colz");
      //c->cd(2);
      //h2full->Draw("colz");
      //c->cd(3);
      //h2Diff->Draw("colz");
      if (0) {
	h2trunc->GetXaxis()->SetRangeUser(0,30);
	h2full->GetXaxis()->SetRangeUser(0,30);
	h2Diff->GetXaxis()->SetRangeUser(0,30);
      }
      if (0) {
 	h2trunc->GetYaxis()->SetRangeUser(0,30);
	h2full->GetYaxis()->SetRangeUser(0,30);
	h2Diff->GetYaxis()->SetRangeUser(0,30);
      }
      if (0) {
	h2trunc->Print("ranges");
      }
      TColorRange_t colors1=_colrange_positive;
      TColorRange_t colors2=_colrange_center;
      colors2=_colrange_positive;
      colors1=_colrange_none; colors2=_colrange_none; gStyle->SetPalette(1);
      colors1=_colrange_default; colors2=_colrange_default;
      drawHistoSubpadAdjustZ(c,1,h2trunc,colors1, 1,1.001,51);
      drawHistoSubpadAdjustZ(c,2,h2full ,colors1, 1,1.001,51);
      drawHistoSubpadAdjustZ(c,3,h2Diff ,colors2, 1,1.001,51);

      TLatex *txt=new TLatex();
      h2Diff->GetZaxis()->SetRangeUser(1e-3,1.0);
      txt->SetTextSize(0.035);
      txt->DrawLatex(700.,9.5, "10^{-3}");


      c->cd(2);
      txt->SetTextSize(0.06);
      txt->DrawLatex(40.,300.,plotLabelV[i]);

      c->Update();
      //break;
    }
  }

  return;
}
