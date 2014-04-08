// from 2014.02.25 ver. studyCovEff.C
// from 2013.08.01 ver. studyEffCovD

#include <TBenchmark.h>
#include "CovariantEff.h"
//#include "../EventScaleFactors/tnpSelectEvents.hh"

const int nExps=100;

// ------------------------------------------------------
// ------------------------------------------------------

int PrepareEtEtaIdx(const esfSelectEvent_t &selData, EtEtaIndexer_t &fidx1, EtEtaIndexer_t &fidx2) {
  int ok= (fidx1.setEtEta(selData.et_1,selData.eta_1) && 
	   fidx2.setEtEta(selData.et_2,selData.eta_2)) ? 1:0;
  if (!ok) {
    //std::cout << "PrepareEtEtaIdx: failing at " << selData.et_1 << ", " << selData.eta_1 << "; " << selData.et_2 << ", " << selData.eta_2 << "\n";
    ok=1;
    if (selData.et_1 > 499.) ok=fidx1.setEtEta(499.,selData.eta_1);
    if (selData.et_2 > 499.) ok= ok && fidx2.setEtEta(499.,selData.eta_2);
    //std::cout << "Et=500GeV correction " << ((ok) ? "suceeded":"failed") << "\n";
  }
  return ok;
}

// ------------------------------------------------------
// ------------------------------------------------------

// If the scaleFactorKind is not full (-1), then syst mode should
// be manually adjusted to reflect this fact. Namely, it is best to 
// include only RECO systematics if, RECO SF is needed.
// or only ID syst, for ID SF.


int studyEffCov_SFsyst(int analysisIs2D,
		       int debugMode,
		       TString systCode="111",
		       int scaleFactorKind=-1, TString variant="") {
  gBenchmark->Start("studyEffCov");

  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return retCodeError;
  }

  if (( DYTools::study2D && (DYTools::nYBins[0]==1)) ||
      (!DYTools::study2D && (DYTools::nYBins[0] >1))) {
    std::cout << "\n\nlinking error:\n";
    std::cout << " either study2D=1 but DYTools::nYBins[0]=1\n";
    std::cout << " or study2D=0 but DYTools::nYBins[0] >1\n";
    std::cout << "Try to remove *.d *.so and recompile\n\n";
    return retCodeError;
  }

  int nTotBins=DYTools::getTotalNumberOfBins(); //same as DYTools::nMassBins in 1D;
  std::cout << "nTotBins=" << nTotBins << "\n";

  //TString confFileName="../config_files/data_vilnius8TeV_regSSD.conf.py";
  TString confFileName="../config_files/data_vilnius8TeV_regSSD_egamma.conf.py";
  //TString confFileName="toy_EtaBins2_var1.conf.py";
  if (variant.Length()) confFileName=Form("toy_%s.conf.py",variant.Data());

  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;
  systMode=DYTools::UNREGRESSED_ENERGY;
  CovariantEffMgr_t mgr;

  TString recoSystFName, idSystFName, hltSystFName;
  TString includedSyst;
  int egammaSystOnly=0;
  if (systCode.Length()>=3) {
    int recoOn=(systCode[0]=='1') ? 1:0;
    int idOn  =(systCode[1]=='1') ? 1:0;
    int hltOn =(systCode[2]=='1') ? 1:0;
    if (recoOn) recoSystFName="../../Results-DYee/root_files_reg/extra_systematics/DY_j22_19712pb_egamma_Unregressed_energy/efficiency_TnP_1D_Full2012_dataRECO_fit-fitEtBins6EtaBins5egamma_PU.root";
    if (idOn) idSystFName="../../Results-DYee/root_files_reg/extra_systematics/DY_j22_19712pb_egamma_Unregressed_energy/efficiency_TnP_1D_Full2012_dataID_fit-fitEtBins6EtaBins5egamma_PU.root";
    if (hltOn) {
      hltSystFName="../../Results-DYee/root_files_reg/extra_systematics/DY_j22_19712pb_egamma_Unregressed_energy/unregHLTSystematics20140226.root";
    }
    if (recoOn+idOn+hltOn==3) includedSyst="-allSyst";
    else if (recoOn+idOn+hltOn==0) includedSyst="-noSyst";
    else {
      if (recoOn) includedSyst.Append("-recoSyst");
      if (idOn)   includedSyst.Append("-idSyst");
      if (hltOn)  includedSyst.Append("-hltSyst");
    }
    if ((recoOn+idOn!=0) && egammaSystOnly) {
      includedSyst.Append("-egammaSystOnly");
    }
  }
  if (1) {
    std::cout << "with systematics: \n";
    std::cout << "  - (reco) " << recoSystFName << "\n";
    std::cout << "  - (id  ) " << idSystFName << "\n";
    std::cout << "  - (hlt ) " << hltSystFName << "\n";
  }

  HERE("calling setup");
  assert(mgr.SetupSFsyst(confFileName,recoSystFName,idSystFName,hltSystFName,nExps,systMode,egammaSystOnly));
  assert(mgr.initOk());

  TString name_extraTag=Form("nMB%d",DYTools::nMassBins);
  name_extraTag.Append(mgr.mgr().userKeyValueAsTString("T&P_ESF_extra"));
  if (nonUniversalHLT) name_extraTag.Append("_asymHLT");
  if (systMode==DYTools::NO_SYST) name_extraTag.Append("_regEn");
  else name_extraTag.Append(generateSystTag(systMode));
  name_extraTag.Append(includedSyst);
  
  TString fileKindStr="SF";
  if (scaleFactorKind!=-1) {
    switch(scaleFactorKind) {
      //case 0:
    case int(DYTools::RECO): fileKindStr.Append("_RECO"); break;
      //case 1:
    case int(DYTools::ID)  : fileKindStr.Append("_ID"); break;
      //case 2:
    case int(DYTools::HLT) : fileKindStr.Append("_HLT"); break;
      //case 3:
    case int(DYTools::HLT_leg1): fileKindStr.Append("_HLTleg1"); break;
      //case 4:
    case int(DYTools::HLT_leg2): fileKindStr.Append("_HLTleg2"); break;
    case 11: 
      fileKindStr.Append("_dtHLT"); 
      if (nonUniversalHLT==0) {
	std::cout << "double-trigger scale factor is requested, but the calcEventEff.C macro contains nonUniversalHLT=0\n";
	return retCodeError;
      }
      else {
	scaleFactorKind=int(DYTools::HLT);
      }
      break;
    default:
      std::cout << "cannot determine scaleFactorKind=" << scaleFactorKind << "\n";
      return retCodeError;
    }
  }

  TString name_covRhoRho=Form("covRhoRho%s_%s_%d",fileKindStr.Data(),name_extraTag.Data(),nExps);
  
  std::cout << "\n\nok. Start studies\n";
  std::cout << "The resulting cov matrix name is " << name_covRhoRho << "\n";

  TString rhoFileName=Form("rhoFile%s_%s_%d.root",fileKindStr.Data(),name_extraTag.Data(),nExps);
  TString covFileName=Form("covRhoFile%s_%s_%d.root",fileKindStr.Data(),name_extraTag.Data(),nExps);

  TString covFileNamePublic=
    TString(Form("covRhoFile%s",fileKindStr.Data())) +
    TString(Form("_el%dD_%dexps_",1+DYTools::study2D,nExps)) +
    //mgr.mgr().userKeyValueAsTString("T&P_ESF_extra") +
    //TString((nonUniversalHLT) ? "_asymHLT" : "") +
    name_extraTag +
    TString(".root");

  if (variant.Length()) {
    TString newEnd=Form("-%s.root",variant.Data());
    rhoFileName.ReplaceAll(".root",newEnd);
    covFileName.ReplaceAll(".root",newEnd);
    covFileNamePublic.ReplaceAll(".root",newEnd);
  }

  std::cout << "\nNames of files to be produced:\n";
  std::cout << " - " << rhoFileName << "\n";
  std::cout << " - " << covFileName << "\n";
  std::cout << " - " << covFileNamePublic << "\n";
  std::cout << std::endl;

  // there is a 'nMB' tag
  //if (DYTools::study2D) {
  //  rhoFileName.ReplaceAll(".root","_2D.root");
  //  covFileName.ReplaceAll(".root","_2D.root");
  //}

  if (debugMode==1) {
    covFileName.ReplaceAll(".root","_Debug.root");
    rhoFileName.ReplaceAll(".root","_Debug.root");
  }

  TVectorD sumWeight(nTotBins);
  //std::cout << sumWeight.GetNoElements() << " vs " << DYTools::getTotalNumberOfBins() << std::endl;
  //return;
  TVectorD sumWeightRho(nTotBins);
  TVectorD sumWeightRhoSqr(nTotBins);
  sumWeight.Zero();
  sumWeightRho.Zero();
  sumWeightRhoSqr.Zero();
    
  TMatrixD sumWeightRho_Rnd(nExps,nTotBins);
  TMatrixD sumWeightRhoSqr_Rnd(nExps,nTotBins);
  sumWeightRho_Rnd.Zero();
  sumWeightRhoSqr_Rnd.Zero();

  std::vector<TMatrixD*> etEtaPairsV;
  etEtaPairsV.reserve(DYTools::nUnfoldingBinsMax);
  EtEtaIndexer_t fidx1(etBinCount,etaBinCount);
  EtEtaIndexer_t fidx2(etBinCount,etaBinCount);

  std::vector<TH1D*> hScaleFIV_150; // flat indexing. 150 bins, like in calcEventEff.C
  std::vector<TH1D*> hScaleFIV_1500;

  std::vector<TString> sample_labels;
  sample_labels.reserve(nTotBins);
  for (int i=0; i<nTotBins; ++i) { sample_labels.push_back(TString(Form("_fib%d",i))); }
  if (!createAnyH1Vec(hScaleFIV_150,"hScaleFIV_150",sample_labels,150,0.,1.5,"scale factor","counts",1) ||
      !createAnyH1Vec(hScaleFIV_1500,"hScaleFIV_1500",sample_labels,1500,0,1.5,"scale factor","counts",1)) {
    std::cout << "failed to prepare scale factor histo-arrays\n";
    return retCodeError;
  }

  if (debugMode!=-1) {


    for (int i=0; i<nTotBins; ++i) {
      TMatrixD *Mptr= new TMatrixD(fidx1.maxFlatIdx(),fidx1.maxFlatIdx());
      if (!Mptr) std::cout << "null Mptr at i=" << i << std::endl;
      (*Mptr)=0;
      etEtaPairsV.push_back(Mptr);
    }

    // PU-weight is present in the weighted of the selected events
    //PUReweight_t PUReweight(PUReweight_t::_Hildreth);
  
    TString puStr = (mgr.mgr().puReweightFlag()) ? "_PU" : "";
    if (!mgr.mgr().fewzFlag()) puStr.Append("_noFEWZ");
    //if (nonUniversalHLT) puStr.Append("_HLTlegs");
    TString selectEventsFName=mgr.mgr().tnpSFSelEventsFullName(systMode,0);
    // remove HLTlegs tag from the file name, if needed
    if (nonUniversalHLT) selectEventsFName.ReplaceAll("_HLTlegs","");

    HERE("opening selectEventsFile=<%s>",selectEventsFName.Data());
    TFile *skimFile=new TFile(selectEventsFName);
    if (!skimFile || !skimFile->IsOpen()) {
      std::cout << "failed to open file <" << selectEventsFName << ">\n";
      assert(0);
    }
    TTree *skimTree = (TTree*)skimFile->Get("Events");
    assert(skimTree);
    esfSelectEvent_t selData;
    selData.setBranchAddress(skimTree);
        

    std::cout << "there are " << skimTree->GetEntries() 
	      << " entries in the <" << selectEventsFName << "> file\n";
    const UInt_t maxEvents=skimTree->GetEntries();
    for (UInt_t ientry=0; ientry<maxEvents; ++ientry) {
      //if (debugMode && (ientry>100000)) break;
      if (debugMode && (ientry%10000!=0)) continue;
      printProgress(200000," ientry=",ientry,maxEvents);
      
      skimTree->GetEntry(ientry);

      // Use generator-level post-FSR mass, y
      //int massBin= DYTools::findMassBin(selData.genMass);
      int massFIdx = DYTools::findIndexFlat(selData.genMass,selData.genY);
      int massBin=massFIdx;
      if (massBin<0) continue;

      double weight=selData.weight;
           
      int fidxOk=PrepareEtEtaIdx(selData, fidx1, fidx2);
      if (!fidxOk) {
	std::cout << "failed to prepare EtEtaIdx\n"; 
	//return;
	continue;
      }

      (*etEtaPairsV[massFIdx])(fidx1.getIdx(),fidx2.getIdx()) += weight;
      if (fidx1!=fidx2) (*etEtaPairsV[massFIdx])(fidx2.getIdx(),fidx1.getIdx()) += weight;

      double scaleFactor = findEventScaleFactor(scaleFactorKind, selData); // HLT formula is inside
      //if (puReweight) weight *= PUReweight.getWeightHildreth(selData.nGoodPV);
      if ( 0 || ( ientry%20000 == 0 )) std::cout << "ientry=" << ientry << ", weight=" << weight << ", scaleFactor=" << scaleFactor << "\n";
          
      //if (massBin==39) std::cout << " sf=" << scaleFactor << " * " << weight << "\n";
      if (0) {
	double tmpEff1=(*dataEff[DYTools::HLT])(fidx1.getEtBin(),fidx1.getEtaBin());
	double tmpEff2=(*dataEff[DYTools::HLT])(fidx2.getEtBin(),fidx2.getEtaBin());
	double tmpEffx=getHLTefficiency(DYTools::DATA,
					fidx1.getEtBin(),fidx1.getEtaBin(),selData.et_1,
					fidx2.getEtBin(),fidx2.getEtaBin(),selData.et_2);
	std::cout << Form("(%3.1lf,%3.1lf; %3.1lf,%3.1lf)  ",selData.et_1,selData.eta_1,selData.et_2,selData.eta_2)
		  << " eff=" << tmpEff1 << "*" << tmpEff2 << "=" << tmpEff1*tmpEff2 << "; vs " << tmpEffx << "\n";
	return retCodeStop;
      }

      hScaleFIV_150[massBin]->Fill(scaleFactor,weight);
      hScaleFIV_1500[massBin]->Fill(scaleFactor,weight);

      sumWeight(massBin) += weight;
      sumWeightRho(massBin) += weight*scaleFactor;
      sumWeightRhoSqr(massBin) += weight*pow(scaleFactor, 2.);
      
      for (int iexp=0; iexp < nExps; ++iexp) {
	double scaleFactorSmeared = findEventScaleFactorSmeared(-1, selData, iexp); // HLT formula is inside
	double extraFactor= mgr.getExtraSF(iexp, fidx1,fidx2);
	//if ((ientry<10)&&(iexp<10)) std::cout << "ientry=" << ientry << ", iexp=" << iexp << ", extraFactor=" << extraFactor << "\n";
	sumWeightRho_Rnd(iexp,massBin) += weight*scaleFactorSmeared*extraFactor;
	sumWeightRhoSqr_Rnd(iexp,massBin) += weight *pow(scaleFactorSmeared*extraFactor, 2.);
      }
    }
    
    delete skimTree;
    delete skimFile;

    std::cout << "loop over selected events done" << std::endl;

    TFile rhoFile(rhoFileName,"recreate");
    sumWeight.Write("sumWeight");
    sumWeightRho.Write("sumWeightRho");
    sumWeightRhoSqr.Write("sumWeightRhoSqr");
    sumWeightRho_Rnd.Write("sumWeightRho_Rnd");
    sumWeightRhoSqr_Rnd.Write("sumWeightRhoSqr_Rnd");
    for (unsigned int i=0; i<etEtaPairsV.size(); ++i) {
      etEtaPairsV[i]->Write(Form("etEtaPairs_mfidx_%u",i));
    }
    saveVec(rhoFile,hScaleFIV_150,"esf_histos_150");
    saveVec(rhoFile,hScaleFIV_1500,"esf_histos_1500");

    writeBinningArrays(rhoFile);
    rhoFile.Close();
  }
  else {
    TFile rhoFile(rhoFileName,"read");
    sumWeight= * (TVectorD*)rhoFile.Get("sumWeight");
    sumWeightRho= * (TVectorD*)rhoFile.Get("sumWeightRho");
    sumWeightRhoSqr= * (TVectorD*)rhoFile.Get("sumWeightRhoSqr");
    sumWeightRho_Rnd = * (TMatrixD*)rhoFile.Get("sumWeightRho_Rnd");
    sumWeightRhoSqr_Rnd = * (TMatrixD*)rhoFile.Get("sumWeightRhoSqr_Rnd");
    for (int i=0; i<DYTools::nUnfoldingBinsMax; ++i) {
      TMatrixD *Mptr= (TMatrixD*)rhoFile.Get(Form("etEtaPairs_mfidx_%d",i));
      etEtaPairsV.push_back(Mptr);
    }
    loadVec(rhoFile,hScaleFIV_150,"esf_histos_150");
    loadVec(rhoFile,hScaleFIV_1500,"esf_histos_1500");

    rhoFile.Close();
  }


  // -=----   calculate further

  TVectorD avgRho(nTotBins), avgRhoErr(nTotBins);
  for (int i=0; i<nTotBins; ++i) {
    double mean=sumWeightRho[i]/sumWeight[i];
    avgRho(i) = mean;
    avgRhoErr(i) = sumWeightRhoSqr[i]/sumWeight[i] - mean*mean;
  }


  TMatrixD avgRho_Rnd(nExps,nTotBins);
  TMatrixD avgRhoSqr_Rnd(nExps,nTotBins);
  std::vector<TMatrixD*> avgRhoRho_RndV;
  avgRhoRho_RndV.reserve(nExps);
  for (int iexp=0; iexp<nExps; iexp++) {
    for (int j=0; j<nTotBins; ++j) {
      avgRho_Rnd(iexp,j) = sumWeightRho_Rnd(iexp,j)/sumWeight(j);
      avgRhoSqr_Rnd(iexp,j) = sumWeightRhoSqr_Rnd(iexp,j)/sumWeight(j);
    }
    TMatrixD* M=new TMatrixD(nTotBins,nTotBins);
    avgRhoRho_RndV.push_back(M);
    for (int i=0; i<nTotBins; ++i) {
      for (int j=0; j<nTotBins; ++j) {
	(*M)(i,j) = avgRho_Rnd(iexp,i) * avgRho_Rnd(iexp,j);
      }
    }
  }

  const int chkMassBin_i=12;
  const int chkMassBin_j=32;
  if (0) {
    double xmin=0.85;
    double xmax=1.05;
    TH2D *hChkMassBin= new TH2D("hChkMassBin",
      Form("Bin(%d,%d) nExps=%d",chkMassBin_i,chkMassBin_j,nExps),
			      200, xmin, xmax, 200, xmin, xmax);
    hChkMassBin->SetDirectory(0);
    for (int iexp=0; iexp<nExps; ++iexp) {
      hChkMassBin->Fill(avgRho_Rnd(iexp,chkMassBin_i),
			avgRho_Rnd(iexp,chkMassBin_j));
    }
    std::cout << Form("covariance in (%d,%d) = %6.4e, corr=%6.4lf\n", chkMassBin_i,chkMassBin_j, hChkMassBin->GetCovariance(),hChkMassBin->GetCorrelationFactor());
    double cf=hChkMassBin->GetCorrelationFactor();
    TString cTitle=
      Form("Bin(%d,%d)   nExps=%d   corr=%4.2lf",chkMassBin_i,chkMassBin_j,nExps,cf);
    TCanvas *cx=new TCanvas("cx",cTitle,500,500);
    AdjustFor2DplotWithHeight(cx,0.1);
    cx->SetLeftMargin(0.17);
    //gStyle->SetPalette(1);
    hChkMassBin->GetXaxis()->SetTitleOffset(0.95);
    hChkMassBin->GetXaxis()->SetTitle("<#rho_{1}>");
    hChkMassBin->GetYaxis()->SetTitleOffset(1.5);
    hChkMassBin->GetYaxis()->SetTitle("<#rho_{2}>");
    hChkMassBin->SetTitle(cTitle);
    hChkMassBin->SetMarkerColor(kBlue+1);
    hChkMassBin->SetMarkerSize(0.95);
    
    set_bottom_white_style(7);
    hChkMassBin->Draw("COLZ");
    cx->Update();
    TString fname=Form("fig_MassBin_%d_%d__nExps%d.png",chkMassBin_i,chkMassBin_j,nExps);
    cx->SaveAs(fname);
    return retCodeStop;
  }

  TVectorD avgRhoMean(nTotBins);
  TVectorD avgRhoRMS(nTotBins);
  avgRhoMean.Zero();
  avgRhoRMS.Zero();
  for (int i=0; i<nTotBins; ++i) {
    double sum=0;
    double sumw2=0;
    for (int iexp=0; iexp<nExps; iexp++) {
      sum+= avgRho_Rnd(iexp,i);
      sumw2 += pow(avgRho_Rnd(iexp,i), 2);
    }
    sum /= double(nExps);
    avgRhoMean(i)=sum;
    //double sumw2=0;
    //for (int iexp=0; iexp<nExps; ++iexp) {
    //  sumw2+= avgRhoSqr_Rnd(iexp,i);
    //}
    sumw2 /= double(nExps);
    sumw2 -= sum*sum;
    avgRhoRMS(i)=(sumw2>0.) ? sqrt(sumw2) : -sqrt(-sumw2);
  }
  
  if (1) {
    int idx=-1;
    double max=maxDiff(avgRho,avgRhoMean,&idx);
    std::cout << "maxDiff(avgRho,avgRhoMean)=" << max << " @ idx=" << idx << "\n";
    std::cout << " - values there: " << avgRho[idx] << " and " << avgRhoMean[idx] << "\n";
  }
  if (0) {
    std::cout << "avgRho="; avgRho.Print();
    std::cout << "avgRhoMean="; avgRhoMean.Print();
    TVectorD diff=avgRho-avgRhoMean;
    //std::cout << "diff="; diff.Print();
    //std::cout << "\nweights"; sumWeight.Print();
  }
  if (0) {
    std::cout << "avgRhoMean and RMS\n";
    for (int i=0; i<nTotBins; ++i) {
      printf(" %2d %10.6lf +- % 8.6lf\n",i,avgRhoMean[i],avgRhoRMS[i]);
    }
    std::cout << "\n";
    //return retCodeStop;
  }

  if (0) { // study (et,eta; et,eta) distributions
    for (int i=0; i<DYTools::nUnfoldingBinsMax; ++i) {
      TMatrixD *Mptr= etEtaPairsV[i];
      TString canvName=Form("canv_FlatIdx_%d",i);
      if (!Mptr) { std::cout << "null Mptr for " << canvName << "\n"; continue; }
      TCanvas *cx=new TCanvas(canvName,canvName,600,600);
      AdjustFor2DplotWithHeight(cx);
      TString hname=Form("h2_etEta_massFI%d_%1.0lf_%1.0lf",i,DYTools::massBinLimits[i],DYTools::massBinLimits[i+1]);
      TH2D *h2=createHisto2D(*Mptr,NULL,hname,hname,_colrange_positive,0,0);
      h2->Draw("colz");
      cx->Update();
    }
  }

  if (0) { // study (et,eta; et,eta) distributions
    //double *loc_etBinLimits=DYTools::getEtBinLimits(etBinning);
    //double *loc_etaBinLimits=DYTools::getEtaBinLimits(etaBinning);
    for (int i=0; i<etBinCount;  i++) std::cout << " " << etBinLimits[i] << "\n";
    TVectorD sumEtLt20(DYTools::nUnfoldingBinsMax);
    TVectorD sumEtLt20Both(DYTools::nUnfoldingBinsMax);
    TVectorD allPairs(DYTools::nUnfoldingBinsMax);
    for (int i=0; i<DYTools::nUnfoldingBinsMax; ++i) {
      TMatrixD *Mptr= etEtaPairsV[i];
      if (!Mptr) { std::cout << "null Mptr for i=" << i << "\n"; continue; }

      int ir=0, ic=0;
      double maxVal=maxElement(*Mptr,&ir,&ic);
      fidx1.setEtEtaBin(ir);
      fidx2.setEtEtaBin(ic);
      std::cout << Form("Mass range %3.0lf -- %4.0lf",DYTools::massBinLimits[i],DYTools::massBinLimits[i+1]) << ", maxVal=" << maxVal << " at (" << ir << "," << ic << ") = {";
      std::cout << fidx1;
      //fidx1.PrintEtEtaVals(etBinLimits,etaBinLimits);
      std::cout << ";";
      std::cout << fidx2;
      //fidx2.PrintEtEtaVals(etBinLimits,etaBinLimits);
      std::cout << "}\n";

      double sumEt_lt_20=0., sumEt_lt_20both=0;
      double sumPairs=0.;
      for (int iEt1=0; iEt1<etBinCount; ++iEt1) {
	const int et1_lt_20 = (etBinLimits[iEt1]<20.) ? 1:0;
	//std::cout << "iEt1=" << iEt1 << ", etBinLimits[iEt1]=" << etBinLimits[iEt1] << ", et1_lt_20=" << et1_lt_20 << "\n";
	for (int iEt2=iEt1; iEt2<etBinCount; ++iEt2) {
	  const int et2_lt_20 = (etBinLimits[iEt2]<20.) ? 1:0;
	  const int sum_et_lt_20=(et1_lt_20 + et2_lt_20 > 0) ? 1:0;
	  const int sum_et_lt_20_both=(et1_lt_20 + et2_lt_20 == 2) ? 1:0;
	  for (int iEta1=0; iEta1<etaBinCount; ++iEta1) {
	    fidx1.setEtEta(iEt1,iEta1);
	    for (int iEta2=iEta1; iEta2<etaBinCount; ++iEta2) {
	      fidx2.setEtEta(iEt2,iEta2);
	      double cnt=(*Mptr)(fidx1.flatEtEtaIdx(),fidx2.flatEtEtaIdx());
	      if (sum_et_lt_20) sumEt_lt_20 += cnt;
	      if (sum_et_lt_20_both) sumEt_lt_20both += cnt;
	      sumPairs += cnt;
	    }
	  }
	}
      }
      sumEtLt20[i]=sumEt_lt_20;
      sumEtLt20Both[i]=sumEt_lt_20both;
      allPairs[i]=sumPairs;
    }

    for (int i=0; i<DYTools::nUnfoldingBinsMax; ++i) {
      std::cout << Form("Mass range %3.0lf -- %4.0lf",DYTools::massBinLimits[i],DYTools::massBinLimits[i+1]) 
		<< ", sum{Et_i<20}=" << sumEtLt20[i] 
		<< ", sum{all}=" << allPairs[i]
		<< ", sum{Et_i<20}/sum{all}=" << sumEtLt20[i]*100/allPairs[i] << "\%\n";
    }

    for (int i=0; i<DYTools::nUnfoldingBinsMax; ++i) {
      std::cout << Form("Mass range %3.0lf -- %4.0lf",DYTools::massBinLimits[i],DYTools::massBinLimits[i+1]) 
		<< ", sum{two Et_i<20}=" << sumEtLt20Both[i] 
		<< ", sum{all}=" << allPairs[i]
		<< ", sum{two Et_i<20}/sum{all}=" << sumEtLt20Both[i]*100/allPairs[i] << "\%\n";
    }

    //delete etBinLimits;
    //delete etaBinLimits;
  }


  if (1) { // save covariance
    HERE("save covariance");

    CovarianceMatrix_t covRhoRho(name_covRhoRho,nTotBins);
    for (int i=0; i<nTotBins; i++) {
      for (int j=0; j<nTotBins; ++j) {
	double sum=0;
	for (int iexp=0; iexp<nExps; ++iexp) {
	  sum+= ( * avgRhoRho_RndV[iexp] )(i,j) / double(nExps);
	}
	covRhoRho.SetValue(i,j, sum - avgRhoMean(i)*avgRhoMean(j));
      }
    }
    //covRhoRho.Print();

    TString name_corrRhoRho=name_covRhoRho;
    name_corrRhoRho.ReplaceAll("cov","corr");
    CovarianceMatrix_t corrRhoRho(name_corrRhoRho,nTotBins);
    assert(corrRhoRho.Correlations(covRhoRho));
    
    TCanvas *cx=new TCanvas("canvRho","cx",1200,600);
    cx->Divide(2,1);
    AdjustFor2DplotWithHeight(cx);
    //TH2D* h2cov =covRhoRho .DrawSubpad(cx,1,_colrange_positive,0,20,1e-3);
    TH2D* h2cov =covRhoRho .DrawSubpad(cx,1,_colrange_positive,0,20);
    TH2D* h2corr=corrRhoRho.DrawSubpad(cx,2,_colrange_positive,0,20,1.);
    std::cout << "histos: " << h2cov->GetName() << " and " << h2corr->GetName() << "\n";
  //set_bottom_white_style(20);
    gStyle->SetPalette(1);
    cx->Update();

    TMatrixD esfM(DYTools::nMassBins,DYTools::nYBinsMax);
    TMatrixD esfMerr(DYTools::nMassBins,DYTools::nYBinsMax);
    TMatrixD esfMpseudo(DYTools::nMassBins,DYTools::nYBinsMax);
    TMatrixD esfMpseudoErr(DYTools::nMassBins,DYTools::nYBinsMax);
    if (!deflattenMatrix(avgRho,esfM) ||
	!deflattenMatrix(avgRhoErr,esfMerr) ||
	!deflattenMatrix(avgRhoMean,esfMpseudo) ||
	!deflattenMatrix(avgRhoRMS,esfMpseudoErr)) {
      std::cout << "failed to deflatten scale factors\n";
      return retCodeError;
    }

    TVectorD esfFromHisto150(nTotBins), esfFromHisto150err(nTotBins);
    TVectorD esfFromHisto1500(nTotBins), esfFromHisto1500err(nTotBins);
    for (int i=0; i<nTotBins; ++i) {
      esfFromHisto150(i)= hScaleFIV_150[i]->GetMean();
      esfFromHisto150err(i)= hScaleFIV_150[i]->GetRMS();
      esfFromHisto1500(i)= hScaleFIV_1500[i]->GetMean();
      esfFromHisto1500err(i)= hScaleFIV_1500[i]->GetRMS();
    }
    TMatrixD esfMFromHisto150(DYTools::nMassBins,DYTools::nYBinsMax);
    TMatrixD esfMFromHisto150err(esfMFromHisto150);
    TMatrixD esfMFromHisto1500(esfMFromHisto150);
    TMatrixD esfMFromHisto1500err(esfMFromHisto150);
    if (!deflattenMatrix(esfFromHisto150,esfMFromHisto150) ||
	!deflattenMatrix(esfFromHisto150err,esfMFromHisto150err) ||
	!deflattenMatrix(esfFromHisto1500,esfMFromHisto1500) ||
	!deflattenMatrix(esfFromHisto1500err,esfMFromHisto1500err)) {
      std::cout << "failed to deflatten scale factors from histos\n";
      return retCodeError;
    }

    TFile fCov(covFileName,"recreate");
    covRhoRho.Write(fCov,"covRhoRho");
    corrRhoRho.Write(fCov,"corrRhoRho");
    cx->Write("canvRhoRho");
    avgRhoMean.Write("scaleFactorFlatIdxArray");
    avgRhoRMS.Write("scaleFactorErrFlatIdxArray");
    esfM.Write("scaleFactor");
    esfMerr.Write("scaleFactorErr_stat");
    esfMpseudo.Write("scaleFactor_pseudo");
    esfMpseudoErr.Write("scaleFactor_pseudoErr");
    esfMpseudoErr.Write("scaleFactorErr");
    esfMFromHisto150.Write("scaleFactor_hb150");
    esfMFromHisto150err.Write("scaleFactor_hb150err");
    esfMFromHisto1500.Write("scaleFactor_hb1500");
    esfMFromHisto1500err.Write("scaleFactor_hb1500err");
    writeBinningArrays(fCov);
    fCov.Close();
    std::cout << "file <" << covFileName << "> recreated\n";


    TFile fResult(covFileNamePublic,"recreate");
    int rangeMin=(DYTools::study2D) ? 25  : 1;
    int rangeMax=(DYTools::study2D) ? 156 : DYTools::nMassBins;
    TH2D *h2covSave=extractSubArea(h2cov,rangeMin,rangeMax,rangeMin,rangeMax,Form("covRhoRho_%d",nExps),1,1); // reset axis!
    TH2D *h2corrSave=extractSubArea(h2corr,rangeMin,rangeMax,rangeMin,rangeMax,Form("corrRhoRho_%d",nExps),1,1); // reset axis!

    TH2D *h2ScaleFactorsFull=createHisto2D(esfM,&esfMpseudoErr,"hScaleFactorFull","scale factor (all bins)",_colrange_positive,0);
    int esfXRangeMin=(DYTools::study2D) ? 2 : 1;
    int esfXRangeMax=(DYTools::study2D) ? 7 : DYTools::nMassBins;
    int esfYRangeMin=(DYTools::study2D) ?  1 : 1;
    int esfYRangeMax=(DYTools::study2D) ? 24 : 1;
    TH2D *h2ScaleFactors=extractSubArea(h2ScaleFactorsFull,esfXRangeMin,esfXRangeMax,esfYRangeMin,esfYRangeMax,"scaleFactor",1, 1); // reset axis!
    h2covSave->Write("covRhoRho");
    h2corrSave->Write("corrRhoRho");
    h2ScaleFactors->Write("scaleFactors");
    writeBinningArrays(fResult);
    fResult.Close();
    std::cout << "file <" << fResult.GetName() << "> saved\n";

    if (0) {
      h2covSave->Print("range");
    }

    if (0) {
      TCanvas *ctest=new TCanvas("ctest","ctest",1200,600);
      ctest->Divide(2,1);
      AdjustFor2DplotWithHeight(ctest);
      ctest->cd(1); h2covSave->Draw("COLZ");
      ctest->cd(2); h2corrSave->Draw("COLZ");
      ctest->Update();
    }

    if (1) {  // for Alexey
      gStyle->SetPalette(1);
      TCanvas *ccov=new TCanvas("ccov","ccov",800,800);
      AdjustFor2DplotWithHeight(ccov);
      h2covSave->Draw("COLZ");
      ccov->Update();
      TCanvas *ccorr=new TCanvas("ccorr","ccorr",800,800);
      AdjustFor2DplotWithHeight(ccorr);
      h2corrSave->Draw("COLZ");
      ccorr->Update();
    }
  }

  //gBenchmark->Show("studyEffCov");
  ShowBenchmarkTime("selectEvents");
  return retCodeOk;
}
