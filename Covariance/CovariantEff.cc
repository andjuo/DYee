#include "CovariantEff.h"

//const TString tnpTag="et6_eta5";
//const TString tnpTag="et6spec_eta5";

//const TString tnpMCFileName=TString("../config_files/sf_mc_${tnpTag}_vilnius.conf").ReplaceAll("${tnpTag}",tnpTag);
//const TString tnpDataFileName=TString("../config_files/sf_data_${tnpTag}_vilnius.conf").ReplaceAll("${tnpTag}",tnpTag);
//const TString filename_mc="../config_files/fall11mc_vilnius.input";

const TString triggerSetString="Full2011_hltEffOld";

//EffArray_t *ro_Data_loc=NULL;
//EffArray_t *ro_MC_loc=NULL;

// ----------------------------------------------------------------
// ----------------------------------------------------------------

// ----------------------------------------------------------------
// ----------------------------------------------------------------

CovariantEffMgr_t::CovariantEffMgr_t() :  
  BaseClass_t("CovariantEffMgr_t"), 
  FInpMgr(),
  FEvtSelector(),
  FPUStr(),
  FTriggers(TrigSet_UNDEFINED,false,0),
  FInitOk(0),
  FRhoRelSystErrs(),
  FRhoExtraFactors()
{
}


// ----------------------------------------------------------------

int CovariantEffMgr_t::Setup(const TString &confFileName, int nExps,
			     DYTools::TSystematicsStudy_t systMode) {

  int res=1;
  int createDestDir=0;
  if (!initManagers(confFileName,DYTools::NORMAL_RUN,FInpMgr,FEvtSelector,FPUStr,createDestDir,systMode)) {
    res=reportError("Setup","failed to initialize managers");
    return res;
  }

  //ro_Data_loc=new EffArray_t[nExps];
  //ro_MC_loc= new EffArray_t[nExps];
  //assert(ro_Data_loc && ro_MC_loc);

  if (nExps>0) {

    // Read efficiency constants from ROOT files
    // This has to be done AFTER configuration file is parsed
    if (!fillEfficiencyConstants( FInpMgr, systMode ) ) {
      return reportError("Setup");
    }

    if (nExps>1) {
      int debug_pseudo_exps=0;
      if (!preparePseudoExps(nExps,debug_pseudo_exps)) 
	return reportError("Setup");
    }
    else std::cout << "CovariantEffMgr::Setup:  nExps=1\n";
  }
  else {
    std::cout << "CovariantEffMgr::Setup:  nExps=0\n";
  }

  FInitOk=res;

  return 1;
}

// ----------------------------------------------------------------

int CovariantEffMgr_t::SetupSFsyst(const TString &confFileName, 
				   const TString &recoSystFName, 
				   const TString &idSystFName, 
				   const TString &hltSystFName, 
				   int nExps,
				   DYTools::TSystematicsStudy_t systMode,
				   int egammaSystOnly) {

  int res=this->Setup(confFileName,nExps,systMode);
  if (!res) {
    res=reportError("SetupSFsyst: failed initial Setup");
    return res;
  }
  
  int rhoSystCount=((recoSystFName.Length()) ? 1:0) +  ((idSystFName.Length()) ? 1:0) + ((hltSystFName.Length()) ? 1:0);

  FInitOk=0; // assume failure

  // Load systematics
  if (rhoSystCount && (nExps>1)) {
    
    FRhoRelSystErrs.reserve(rhoSystCount);
    FRhoExtraFactors.reserve(nExps);
    std::vector<int> err_ids; // we need to track RECO syst due to merging of bins
    err_ids.reserve(rhoSystCount);
    
    TString etStr=FInpMgr.getTNP_etBinningString();
    std::cout << "etStr=" << etStr << "\n";
    DYTools::TEtBinSet_t loc_etBinning=DetermineEtBinSet(etStr);
    int loc_etBinCount=DYTools::getNEtBins(loc_etBinning);
    //etBinLimits=DYTools::getEtBinLimits(etBinning);
    
    TString etaStr=FInpMgr.getTNP_etaBinningString();
    std::cout << "etaStr=" << etaStr << "\n";
    DYTools::TEtaBinSet_t loc_etaBinning=DetermineEtaBinSet(etaStr);
    int loc_etaBinCount=DYTools::getNEtaBins(loc_etaBinning);
    //etaBinLimits=DYTools::getEtaBinLimits(etaBinning);
    
    TString systFieldName=(egammaSystOnly) ? "sf_syst_rel_error_egamma" : "sf_syst_rel_error";

    if (recoSystFName.Length()) {
      TMatrixD *Mreco= loadMatrix(recoSystFName,systFieldName,loc_etBinCount,loc_etaBinCount);
      if (!Mreco) return 0;
      FRhoRelSystErrs.push_back(Mreco);
      err_ids.push_back(int(DYTools::RECO));
    }
    if (idSystFName.Length()) {
      TMatrixD *Mid= loadMatrix(idSystFName,systFieldName,loc_etBinCount,loc_etaBinCount);
      if (!Mid) return 0;
      FRhoRelSystErrs.push_back(Mid);
      err_ids.push_back(int(DYTools::ID));
    }
    if (hltSystFName.Length()) {
      TMatrixD* Mhlt= loadMatrix(hltSystFName,systFieldName,loc_etBinCount,loc_etaBinCount,0);
      if (Mhlt) {
	FRhoRelSystErrs.push_back(Mhlt);
	err_ids.push_back(int(DYTools::HLT));
      }
      else if (nonUniversalHLT) {
	// assume that the field was not found
	TMatrixD* MhltLeg1mc= loadMatrix(hltSystFName,"hltLeg1_mc_syst_rel_error",loc_etBinCount,loc_etaBinCount,1);
	if (!MhltLeg1mc) return 0;
	TMatrixD* MhltLeg2mc= loadMatrix(hltSystFName,"hltLeg2_mc_syst_rel_error",loc_etBinCount,loc_etaBinCount,1);
	if (!MhltLeg2mc) return 0;
	TMatrixD* MhltLeg1data= loadMatrix(hltSystFName,"hltLeg1_data_syst_rel_error",loc_etBinCount,loc_etaBinCount,1);
	if (!MhltLeg1data) return 0;
	TMatrixD* MhltLeg2data= loadMatrix(hltSystFName,"hltLeg2_data_syst_rel_error",loc_etBinCount,loc_etaBinCount,1);
	if (!MhltLeg2data) return 0;
	// Add this systematic error to the one already employed
	// The containers are in calcEventEff.C
	for (int isData=0; isData<2; ++isData) {
	  const int iLeg1=DYTools::HLT_leg1;
	  const int iLeg2=DYTools::HLT_leg2;
	  TMatrixD *effLeg1= (isData) ? dataEff[iLeg1] : mcEff[iLeg1];
	  TMatrixD *effLeg1Err= (isData) ? dataEffAvgErr[iLeg1] : mcEffAvgErr[iLeg1];
	  TMatrixD *effLeg2= (isData) ? dataEff[iLeg2] : mcEff[iLeg2];
	  TMatrixD *effLeg2Err= (isData) ? dataEffAvgErr[iLeg2] : mcEffAvgErr[iLeg2];
	  TMatrixD *eSystLeg1= (isData) ? MhltLeg1data : MhltLeg1mc;
	  TMatrixD *eSystLeg2= (isData) ? MhltLeg2data : MhltLeg2mc;
	  for (int iEt=0; iEt<loc_etBinCount; ++iEt) {
	    for (int iEta=0; iEta<loc_etaBinCount; ++iEta) {
	      double e1=pow( (*effLeg1Err)[iEt][iEta], 2 ) + 
		pow( (*effLeg1)[iEt][iEta] * (*eSystLeg1)[iEt][iEta], 2 );
	      (*effLeg1Err)[iEt][iEta] = sqrt(e1);
	      double e2=pow( (*effLeg2Err)[iEt][iEta], 2 ) + 
		pow( (*effLeg2)[iEt][iEta] * (*eSystLeg2)[iEt][iEta], 2 );
	      (*effLeg2Err)[iEt][iEta] = sqrt(e2);
	    }
	  }
	}
	delete MhltLeg2data;
	delete MhltLeg1data;
	delete MhltLeg2mc;
	delete MhltLeg1mc;
      }
      else {
	HERE("\nsystematics for universal HLT is not known\n");
	return 0;
      }
    }

    for (int i=0; i<nExps; i++) {
      
      TMatrixD extra(loc_etBinCount,loc_etaBinCount);
      for (int iEt=0; iEt<loc_etBinCount; ++iEt) {
	for (int iEta=0; iEta<loc_etaBinCount; ++iEta) {
	  extra(iEt,iEta) = 1.;
	}
      }
      
      for (unsigned int kind=0; kind<FRhoRelSystErrs.size(); ++kind) {
	// kind has to be external loop, since we might merge some bins (and force factors to be identical)
	for (int iEt=0; iEt<loc_etBinCount; ++iEt) {
	  for (int iEta=0; iEta<loc_etaBinCount; ++iEta) {
	    // In the special case of the RECO efficiency for low Et
	    // electrons, some eta bins are MERGED in tag and probe.
	    // However the binning is kept standard, so the values
	    // for the efficiencies in the merged bins are the same,
	    // and the errors are 100% correlated. Take this into
	    // account and make the smearing 100% correlated as well.
	    if ((err_ids[kind]==int(DYTools::RECO))
		&& DYTools::mergeEtBins(etaBinning)
		&& (getEtBinLimits(etBinning))[iEt+1] <= 20.0 
		&& (iEta == 1 || iEta == 4)    ) {
	      // We are dealing with merged bins of the RECO efficiency with the right binning
	      // For iEta == 1 or 4, fall back to the values for iEta == 0 or 3.
	      extra(iEt,iEta) = extra(iEt,iEta);
	    }
	    else {
	      // The default case, all other efficiencies and bins
	      const double relErr=(*FRhoRelSystErrs[kind])(iEt,iEta);
	      double factor=( 1. + gRandom->Gaus(0.0,1.0) * relErr );
	      extra(iEt,iEta) *= factor;
	    }
	  }
	}
      }
      //std::cout << "iexp=" << i << " rhoExtraFactors: "; extra.Print();
      FRhoExtraFactors.push_back(new TMatrixD(extra));
    }
  }
  else {
    if (nExps==1) std::cout << "CovariantEffMgr::SetupSFsyst: nExps=1\n";
  }

  FInitOk=1;

  return 1;
}


// ----------------------------------------------------------------



// ----------------------------------------------------------------
// ----------------------------------------------------------------

StudyArr_t::StudyArr_t(const TString &baseName, int nMassBins, int nExps,
		       int binCount, double minVal, double maxVal) :
  BaseClass_t(Form("StudyArr_t(%s)",baseName.Data())),
  FVec(),
  FBaseName(baseName),
  FReady(0) {
  this->Init(nMassBins,nExps, binCount,minVal,maxVal);
}


// ----------------------------------------------------------------

int StudyArr_t::Init(int nMassBins, int nExps, 
		     int binCount, double minVal, double maxVal) {
  FReady=0;
  if (binCount==0) return 0;
  if (FVec.size()) { BaseClass_t::printWarning("::Init will clear FVec"); }
  FVec.clear();
  FVec.reserve(nMassBins);
  for (int i=0; i<nMassBins; ++i) {
    std::vector<TH1D*> *vec=new std::vector<TH1D*>();
    FVec.push_back(vec);
    vec->reserve(nExps);
    for (int j=0; j<nExps; ++j) {
      TString name=Form("%s_%d_exp%d",FBaseName.Data(),i,j);
      TH1D* h=new TH1D(name,name, binCount,minVal,maxVal);
      h->SetDirectory(0);
      h->Sumw2();
      vec->push_back(h);
    }
  }
  FReady=1;
  return 1;
}

// ----------------------------------------------------------------
// ----------------------------------------------------------------

ScaleFactor_t::ScaleFactor_t(int set_kind,
			     int use_etBinCount, int use_etaBinCount) :
  BaseClass_t("ScaleFactor_t"),
  fTmpIdx(use_etBinCount, use_etaBinCount),
  fKind(set_kind),
  fSF(fTmpIdx.maxFlatIdx(),fTmpIdx.maxFlatIdx()),
  fSFrndV()
{
  fSF.Zero();
}

// ----------------------------------------------------------------

double ScaleFactor_t::eventScaleFactor(const esfSelectEvent_t &selData) const {
  return findEventScaleFactor(fKind,selData);
}

// ----------------------------------------------------------------

double ScaleFactor_t::rndEventScaleFactor(int iexp,
				   const esfSelectEvent_t &selData) const {
  return findEventScaleFactorSmeared(fKind,selData,iexp);
}

// ----------------------------------------------------------------

int ScaleFactor_t::init (int set_kind,
			 DYTools::TEtBinSet_t set_etBinning,
			 DYTools::TEtaBinSet_t set_etaBinning,
			 int nExps)
{
  ClearVec(fSFrndV);
  fKind=set_kind;

  EtEtaIndexer17_t fi1(fTmpIdx), fi2(fTmpIdx);
  fSF.Zero();

  int loc_etBinCount=  DYTools::getNEtBins(set_etBinning) + 1;
  int loc_etaBinCount= DYTools::getNEtaBins(set_etaBinning);
  double *loc_etBinLimits_tmp = DYTools::getEtBinLimits(set_etBinning);
  double *loc_etBinLimits = new double[loc_etBinCount+1];
  double *loc_etaBinLimits= DYTools::getEtaBinLimits(set_etaBinning);
  double *loc_etBinCenter= new double[loc_etBinCount];
  double *loc_etaBinCenter=new double[loc_etaBinCount];

  if ((loc_etBinCount*loc_etaBinCount != fSF.GetNrows()) ||
      (loc_etBinCount*loc_etaBinCount != fSF.GetNcols())) {
    reportError("init: cannot change dimension after object creation");
    std::cout << " current dim " << fSF.GetNrows() << " x " << fSF.GetNcols() << "\n";
    std::cout << " new dim " << loc_etBinCount << " x " << loc_etaBinCount << "\n";
    return 0;
  }

  if (!loc_etBinLimits || !loc_etBinLimits_tmp || !loc_etaBinLimits ||
      !loc_etBinCenter || !loc_etaBinCenter) {
    this->reportError("init: Failed to allocate temporary arrays");
    return 0;
  }

  int shift=0;
  for (int ibin=0; ibin<loc_etBinCount+1; ++ibin) {
    if ((loc_etBinLimits_tmp[ibin+shift]  <17.) &&
	(loc_etBinLimits_tmp[ibin+shift+1]>17.)) {
      shift=-1;
      loc_etBinLimits[ibin]=  loc_etBinLimits_tmp[ibin];
      loc_etBinLimits[ibin+1]=17;
      ibin++;
    }
    else {
      loc_etBinLimits[ibin] = loc_etBinLimits_tmp[ibin+shift];
    }
  }

  for (int i=0; i<loc_etBinCount+1; ++i) {
    std::cout << " i=" << i << ", val=" << loc_etBinLimits[i] << "\n";
  }

  for (int ibin=0; ibin<loc_etBinCount; ++ibin) {
    loc_etBinCenter[ibin] =
      0.5*(loc_etBinLimits[ibin] + loc_etBinLimits[ibin+1]);
  }

  for (int i=0; i<loc_etBinCount; ++i) {
    std::cout << " i=" << i << ", center=" << loc_etBinCenter[i] << "\n";
  }

  for (int ibin=0; ibin<loc_etaBinCount; ++ibin) {
    std::cout << "iEtaBin=" << ibin << ", val=" << loc_etaBinLimits[ibin] << "\n";
    loc_etaBinCenter[ibin] =
      0.5*(loc_etaBinLimits[ibin] + loc_etaBinLimits[ibin+1]);
  }

  fSFrndV.reserve(nExps);
  for (int iexp=0; iexp<nExps; ++iexp) {
    TMatrixD *m= new TMatrixD(fSF);
    fSFrndV.push_back(m);
  }

  //HERE("fill the matrix\n");
  fi1.Start();
  do {
    double et1 = loc_etBinCenter [fi1.getEtBin()];
    double eta1= loc_etaBinCenter[fi1.getEtaBin()];
    fi2.Start();
    do {
      double et2 = loc_etBinCenter [fi2.getEtBin()];
      double eta2= loc_etaBinCenter[fi2.getEtaBin()];

      //HERE("zx");
      //std::cout << "fi1 =" << fi1 << ", fi2=" << fi2 << "\n";
      //std::cout << "(et,eta)=( " << et1 << "," << eta1 << "), (" << et2
      //<< "," << eta2 << ")\n";
      //std::cout << Form("(%4.1lf-%4.1lf,%4.2lf-%4.2lf) ",loc_etBinLimits[fi1.getEtBin()],loc_etBinLimits[fi1.getEtBin()+1],loc_etaBinLimits[fi1.getEtaBin()],loc_etaBinLimits[fi1.getEtaBin()+1]);
      //std::cout << Form("(%4.1lf-%4.1lf,%4.2lf-%4.2lf) ",loc_etBinLimits[fi2 .getEtBin()],loc_etBinLimits[fi2.getEtBin()+1],loc_etaBinLimits[fi2.getEtaBin()],loc_etaBinLimits[fi2.getEtaBin()+1]);
      //HERE("get the scale factor\n");

      // call function from calcEventEff.C
      double value= findEventScaleFactor(fKind,et1,eta1,et2,eta2);
      if (value!=value) value=1e9;
      //HERE(" value=%lf",value);
      fSF(fi1.getIdx(),fi2.getIdx()) = value;
      //fSF(fi2.getIdx(),fi1.getIdx()) = value;

      for (int iexp=0; iexp<nExps; ++iexp) {
	value= findEventScaleFactorSmeared(fKind,et1,eta1,et2,eta2,iexp);
	if (value!=value) value=1e9;
	(*fSFrndV[iexp])(fi1.getIdx(),fi2.getIdx()) = value;
	//(*fSFrndV[iexp])(fi2.getIdx(),fi1.getIdx()) = value;
      }

    } while(fi2.Next()!=-1);
  } while(fi1.Next()!=-1);
  //HERE("done");

  delete loc_etBinLimits_tmp;
  delete loc_etBinLimits;
  delete loc_etaBinLimits;
  delete loc_etBinCenter;
  delete loc_etaBinCenter;

  return 1;
}

// ----------------------------------------------------------------



// ----------------------------------------------------------------
// ----------------------------------------------------------------

