#include "../Include/DYTools.hh"
#include "../CrossSection/crossSectionFnc.hh"

const double ZMass_low = 60.;
const double ZMass_high=120.;

// the variable codeDebugFilePath indicates the location from which
// files of another compatible package 
// (another version of DYee or DrellYanDMDY) should be loaded.
// If no debug path is assigned, the functions will use the files
// of the current analysis.
TString codeDebugFilePath;
//TString codeDebugFilePath= "/home/andriusj/cms/DYee8TeV-20140118-2D/root_files/constants/DY_j22_19712pb/";

//=== FUNCTION DECLARATIONS ======================================================================================

// ----------------------------------------------
/*
int getNormalizationMBinRange(int &firstNormBin, int &lastNormBin) {
  firstNormBin=-1;
  lastNormBin=-1;

  for (int i=0; i<=DYTools::nMassBins; i++) {
    if (DYTools::massBinLimits[i] == ZMass_low ) firstNormBin=i;
    if (DYTools::massBinLimits[i] == ZMass_high) lastNormBin=i;
  }

  if(firstNormBin == -1 || lastNormBin == -1){
    printf("\nERROR: normalization limits are misaligned with mass binning!\n\n");
    return 0;
  }
  printf("\nCross section normalization is to the bins %d - %d from %5.1f to %5.1f GeV\n", 
	 firstNormBin, lastNormBin, 
	 DYTools::massBinLimits[firstNormBin], DYTools::massBinLimits[lastNormBin+1]);

  return 1;
}
*/

// ----------------------------------------------

//=== Cross-section calculation FUNCTION IMPLEMENTATIONS ======================================================================================


// ----------------------------------------------

int unfoldDetResolution(const InputArgs_t &inpArg, const HistoPair2D_t &ini, HistoPair2D_t &fin) {
  if (inpArg.silentMode()<2) HERE(" -- unfoldDetResolution");
  // To load the unfolding matrix we have to create variables used
  // for the construction of
  // the class instance as well as autoLoadFromFile
  TString constDir= inpArg.inpMgr()->constDir(inpArg.systMode(),0);
  TString fnameTag= UnfoldingMatrix_t::generateFNameTag(inpArg.systMode(),
							inpArg.externalSeed());
  TMatrixD *UnfM=NULL;
  int inverse=1;
  const int load_debug_file=(codeDebugFilePath.Length()) ? 1:0;
  if ( load_debug_file ) {
    constDir=codeDebugFilePath;
    fnameTag=DYTools::analysisTag + TString("_PU");
  }
  UnfM=UnfoldingMatrix_t::LoadUnfM("detResponse", constDir, fnameTag, inverse);
  if (!UnfM) return 0;
  // unfolding does not include additional error of the unfolding matrix
  int res= fin.unfold(*UnfM, ini);
  if (!res) return 0;
  delete UnfM;
  return 1;
}

// ----------------------------------------------

int fsrCorrection_det(const InputArgs_t &inpArg, const HistoPair2D_t &ini, HistoPair2D_t &fin) {
  if (inpArg.silentMode()<2) HERE("fsrCorrection_det");
  // To load the unfolding matrix we have to create variables used
  // for the construction of
  // the class instance as well as autoLoadFromFile
  TString constDir= inpArg.inpMgr()->constDir(inpArg.systMode(),0);
  TString fnameTag= UnfoldingMatrix_t::generateFNameTag(inpArg.systMode(),
							inpArg.externalSeed());
  TMatrixD *UnfM=NULL;
  int inverse=1;
  const int load_debug_file=(codeDebugFilePath.Length()) ? 1:0;
  if ( load_debug_file ) {
    constDir=codeDebugFilePath;
    fnameTag=DYTools::analysisTag + TString("_PU");
  }
  UnfM=UnfoldingMatrix_t::LoadUnfM("fsrDETgood", constDir, fnameTag, inverse);
  if (!UnfM) return 0;
  // unfolding does not include additional error of the unfolding matrix
  int res= fin.unfold(*UnfM, ini);
  if (!res) return 0;
  delete UnfM;
  return 1;
}

// ----------------------------------------------

int fsrCorrection_fullSpace(const InputArgs_t &inpArg, const HistoPair2D_t &ini, HistoPair2D_t &fin) {
  if (inpArg.silentMode()<2) HERE("fsrCorrection_fullSpace");
  // To load the unfolding matrix we have to create variables used
  // for the construction of
  // the class instance as well as autoLoadFromFile
  TString constDir= inpArg.inpMgr()->constDir(inpArg.systMode(),0);
  TString fnameTag= UnfoldingMatrix_t::generateFNameTag(inpArg.systMode(),
							inpArg.externalSeed());
  TMatrixD *UnfM=NULL;
  int inverse=1;
  const int load_debug_file=(codeDebugFilePath.Length()) ? 1:0;
  if (load_debug_file) {
    constDir=codeDebugFilePath;
    fnameTag=DYTools::analysisTag + TString("_PU");
  }
  UnfM=UnfoldingMatrix_t::LoadUnfM("fsrGood", constDir, fnameTag, inverse);
  if (!UnfM) return 0;
  // unfolding does not include additional error of the unfolding matrix
  int res= fin.unfold(*UnfM, ini);
  if (!res) return 0;
  delete UnfM;
  return 1;
}

// ----------------------------------------------

int efficiencyCorrection(const InputArgs_t &inpArg, const HistoPair2D_t &ini, HistoPair2D_t &fin) {
  if (inpArg.silentMode()<2) HERE(" -- efficiencyCorrection");
  TString effCorrFName=inpArg.inpMgr()->correctionFullFileName("efficiency",inpArg.systMode(),0);
  TH2D *hEff=NULL;
  const int load_debug_file=(codeDebugFilePath.Length()) ? 1:0;
  if ( ! load_debug_file ) {
    hEff=LoadHisto2D("hEfficiency",effCorrFName,"",1);
  }
  else {
    TString constDir=codeDebugFilePath;
    effCorrFName=constDir + TString(Form("event_efficiency_constants%dD.root",DYTools::study2D+1));
    std::cout << "loading from <" << effCorrFName << ">\n";
    hEff=LoadMatrixFields(effCorrFName,0,"efficiencyArray","efficiencyErrArray",1);
  }
  //std::cout << "effCorrFName=<" << effCorrFName << ">\n";
  int res=(hEff!=NULL) ? 1:0;
  if (res) res=fin.divide(ini,hEff,!inpArg.includeCorrError());
  if (!res) std::cout << "error in efficiencyCorrection\n";
  return res;
}
// ----------------------------------------------

int efficiencyScaleCorrection(const InputArgs_t &inpArg, const HistoPair2D_t &ini, HistoPair2D_t &fin) {
  if (inpArg.silentMode()<2) HERE(" -- efficiencyScaleCorrection");
  DYTools::TSystematicsStudy_t systMode=inpArg.systMode();
  if (1) {
    if (!inpArg.silentMode()) {
      std::cout << "Efficiency scale factors have built-in systMode\n";
    }
    systMode=DYTools::NO_SYST;
  }
  TString rhoCorrFName=inpArg.inpMgr()->correctionFullFileName("scale_factors_asymHLT",systMode,0);
  TH2D *hRho=NULL;
  const int load_debug_file=(codeDebugFilePath.Length()) ? 1:0;
  if ( ! load_debug_file ) {
    TString sfTag=inpArg.inpMgr()->userKeyValueAsTString("SpecialESFTag");
    if (!inpArg.silentMode())  std::cout << "\n\n";
    std::cout << "\tsfTag=<" << sfTag << ">\n";
    if (!inpArg.silentMode()) std::cout << "\n";
    if (sfTag.Length()) {
      // special
      rhoCorrFName="/home/andriusj/cms/DYee-20140407-ESF-systVars/Covariance/";
      //covRhoFileSF_nMB7_asymHLT_FSR_5plus-allSyst_100.root
      rhoCorrFName.Append(Form("covRhoFileSF_nMB%d_asymHLT_%s.root",
			       DYTools::nMassBins,sfTag.Data()));
      std::cout << "special name of scale factor file <" << rhoCorrFName
		<< ">\n";
    }
    else {
      // default manipulation
      sfTag=inpArg.inpMgr()->userKeyValueAsTString("ScaleFactorTag");
      if (sfTag.Length()) {
	std::cout << "special scale factor tag=<" << sfTag << ">\n";
	TString constTag=inpArg.inpMgr()->constTag();
	rhoCorrFName.ReplaceAll(constTag,sfTag);
      }
    }
    //hRho=LoadHisto2D("hEffScaleFactor",rhoCorrFName,"",1);
    int checkBinning=0;
    hRho=LoadMatrixFields(rhoCorrFName,checkBinning,"scaleFactor","scaleFactorErr",1);
  }
  else {
    TString constDir=codeDebugFilePath;
    rhoCorrFName=constDir + TString(Form("scale_factors_%dD_Full2012_hltEffOld_PU.root",DYTools::study2D+1));
    hRho=LoadMatrixFields(rhoCorrFName,1,"scaleFactor","scaleFactorErr",1);
  }
  //std::cout << "rhoCorrFName=<" << rhoCorrFName << ">\n";
  int res=(hRho!=NULL) ? 1:0;
  if (res) res=fin.divide(ini,hRho,!inpArg.includeCorrError());
  if (!res) std::cout << "error in efficiencyScaleCorrection\n";
  return res;
}

// ----------------------------------------------

int acceptanceCorrection(const InputArgs_t &inpArg, const HistoPair2D_t &ini, HistoPair2D_t &fin) {
  if (inpArg.silentMode()<2) HERE(" -- acceptanceCorrection");
  DYTools::TSystematicsStudy_t systMode=inpArg.systMode();
  if ((systMode==DYTools::PILEUP_5plus) ||
      (systMode==DYTools::PILEUP_5minus)) {
    std::cout << "correction! acceptance does not depend on pile-up\n";
    systMode=DYTools::NO_SYST;
  }
  TString accCorrFName=inpArg.inpMgr()->correctionFullFileName("acceptance",systMode,0);
  TH2D* hAcc=NULL;
  const int load_debug_file=(codeDebugFilePath.Length()) ? 1:0;
  if ( ! load_debug_file ) {
    hAcc=LoadHisto2D("hAcceptance",accCorrFName,"",1);
  }
  else {
    TString constDir=codeDebugFilePath;
    accCorrFName=constDir + TString(Form("acceptance_constants%dD.root",DYTools::study2D+1));
    std::cout << "loading from <" << accCorrFName << ">\n";
    hAcc=LoadMatrixFields(accCorrFName,1,"acceptanceMatrix","acceptanceErrMatrix",1);
  }
  //std::cout << "accCorrFName=<" << accCorrFName << ">\n";

  int res=(hAcc!=NULL) ? 1:0;
  if (res) res=fin.divide(ini,hAcc,!inpArg.includeCorrError());
  if (!res) std::cout << "error in acceptanceCorrection\n";
  return res;
}


// ----------------------------------------------
// Adding systematic errors
// ----------------------------------------------

int addSystError_DetResUnf_unfold(const InputArgs_t &ia, HistoPair2D_t &hp, HistoPair2D_t *resHP) {
  return 0;
  /*
  HERE(" --- addSystError_DetResUnf_unfold");
  if (resHP==NULL) resHP=&hp; else resHP->assign(hp);
  TString systErrDir= inpArg.inpMgr()->systErrDir(inpArg.systMode(),0);
  TString fnameTag= DYTools::analysisTag;
  if (1) {
    systErrDir="/home/andriusj/cms/CMSSW_3_8_4/src/DrellYanDMDY-20130131/root_files/systematics/DY_m10+pr+a05+o03+pr_4839pb/";
  }
  */
  return 1;
}

// ----------------------------------------------

int addSystError_DetResUnf_escale(const InputArgs_t &ia, HistoPair2D_t &hp, HistoPair2D_t *resHP) {
  return 0;
  if (ia.silentMode()<2) HERE(" --- addSystError_DetResUnf_escale");
  if (resHP==NULL) resHP=&hp; else resHP->assign(hp);
  int res=(resHP==NULL) ? 0:1;
  //TString systErrDir= inpArg.inpMgr()->systErrDir(inpArg.systMode(),0);
  //TString fnameTag= DYTools::analysisTag;
  TString fname;//=inpArg.inpMgr()->systematicsFullFileName("escale",inpArg.systMode(),0);
  if (1) {
    //systErrDir="/home/andriusj/cms/CMSSW_3_8_4/src/DrellYanDMDY-20130131/root_files/systematics/DY_m10+pr+a05+o03+pr_4839pb/";
    fname=TString("/home/andriusj/cms/CMSSW_3_8_4/src/DrellYanDMDY-20130131/root_files/systematics/DY_m10+pr+a05+o03+pr_4839pb/escale_systematics1D_tmp.root");
  }
  TH2D *escaleSyst=NULL;
  if (res) { 
    escaleSyst=LoadMatrixFields(fname,1,"","escaleSystError",1);
    res=(escaleSyst==NULL) ? 0:1;
  }
  if (res) res=resHP->addSystErrPercent(escaleSyst);
  if (escaleSyst) delete escaleSyst;
  if (!res) std::cout << "error in addSystError_DetResUnf_escale\n";
  return res;
}


// ----------------------------------------------
// ----------------------------------------------

int saveResult(const InputArgs_t &ia, const HistoPair2D_t &hp, 
	       TString extraTag) {
  if (!ia.silentMode()) {
    std::cout << "would save result with extraTag_inp=<" << extraTag << ">\n";
    TString fname=ia.inpMgr()->crossSectionFullFileName(ia.systMode(),
							ia.csKind(),
							1,0);
    if (!ia.includeCorrError()) extraTag.Append("_noCorrErr");
    if (extraTag.Length()) fname.ReplaceAll(".root",extraTag+TString(".root"));
    if (ia.resNameBase().Length()) {
      Ssiz_t index=fname.Last('/');
      if (index==-1) fname.Prepend(ia.resNameBase());
      else fname.Insert(index,ia.resNameBase());
    }
    std::cout << "output fname would be <" << fname << ">\n";

    hp.print();

    HistoPair2D_t hpPerLumi(hp.GetName() + TString("_divLumi"),hp);
    hpPerLumi.scale(1./ia.inpMgr()->totalLumi());
    //hpPerLumi.print();
    //std::cout << "\n\n";
    CSResults_t zpResult;
    zpResult.assignZpeakCS(hpPerLumi);
    std::cout << "normalization factor: " << zpResult << "\n";
    HistoPair2D_t hpNorm(hp.GetName() + TString("_norm"),hp);
    hpNorm.scale(1./ia.inpMgr()->totalLumi());
    zpResult.doNormalize(hpNorm,ia.allNormErrorIsSyst());

    if (!ia.noSave()) {
      int res=1;
      CreateDir(fname,1);
      TFile fout(fname,"recreate");
      if (res) res=hp.Write(fout);
      if (res) res=hpPerLumi.Write(fout);
      // for plotting
      if (res) {
	res= (hp.Write(fout,"auto","count") &&
	      hpPerLumi.Write(fout,"auto","xsec") &&
	      hpNorm.Write(fout,"auto","xsecNorm") &&
	      zpResult.write("auto/norm")) ? 1:0;
      }
      if (res) writeBinningArrays(fout);
      fout.Close();
      if (!res) std::cout << "error while saving -- ";
      std::cout << "saved to file <" << fout.GetName() << ">\n";
      if (!res) return 0;
    }
  }
  return 1;
}


// ----------------------------------------------
// ----------------------------------------------
//=== MAIN MACRO =================================================================================================

// cross section calculation
int calculateCSdistribution(const InputArgs_t &ia, const HistoPair2D_t &hp_ini, 
			    DYTools::TCrossSectionKind_t csKind,
			    HistoPair2D_t &hp_fin) {

  {
    using namespace DYTools;
    int debugPrint=0;
    if (!checkCSKind(csKind,debugPrint,
		     4,
		     _cs_preFsr, _cs_preFsrDet,
		     _cs_postFsr,_cs_postFsrDet)) {
      std::cout << "in calculateCSdistribution\n";
      return 0;
    }
  }

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //============================================================================================================== 

  int res=1;

  TString hpInpName= ia.resNameBase() + TString("_inp");
  TString hpOutName= ia.resNameBase() + TString("_out");
  HistoPair2D_t hpInp(hpInpName,hp_ini);
  HistoPair2D_t hpOut(hpOutName);

  HistoPair2D_t hpUnfoldedYield(ia.resNameBase() + TString("unfYield"));
  HistoPair2D_t hpEffCorrYield(ia.resNameBase() + TString("effCorrYield"));
  HistoPair2D_t hpEffRhoCorrYield(ia.resNameBase() + TString("effRhoCorrYield"));

  if (ia.needsDetUnfolding()) {
    // unfolding correction
    if (res) res=unfoldDetResolution(ia, hpInp, hpUnfoldedYield);
    //if (res) res=addSystError_DetResUnf_unfold(ia, hpUnfoldedYield);
    //if (res) res=addSystError_DetResUnf_escale(ia, hpUnfoldedYield);
    if (res) res=saveResult(ia,hpUnfoldedYield,ia.resNameBase() + TString("unf"));
  }
  else {
    res=hpUnfoldedYield.assign(hpInp);
  }

  // efficiency correction
  if (ia.needsEffCorr()) {
    if (res) res=efficiencyCorrection(ia, hpUnfoldedYield,hpEffCorrYield);
    if (res) res=saveResult(ia,hpEffCorrYield,
			    ia.resNameBase() + TString("eff"));
  }
  else {
    res=hpEffCorrYield.assign(hpUnfoldedYield);
  }

  // efficiency scale correction
  if (ia.needsEffScaleCorr()) {
    if (res) res=efficiencyScaleCorrection(ia, hpEffCorrYield,hpEffRhoCorrYield);
    if (res) res=saveResult(ia,hpEffRhoCorrYield,
			    ia.resNameBase() + TString("effRho"));
  }
  else {
    res=hpEffRhoCorrYield.assign(hpEffCorrYield);
  }


  // post-FSR in detector acceptance
  if (csKind==DYTools::_cs_postFsrDet) {
    if (res) res=hp_fin.assign(hpEffRhoCorrYield);
    if (!res) {
      std::cout << "failed to produce " << CrossSectionKindName(csKind) << "\n";
    }
    return res;
  }

  // pre-FSR in acceptance result
  if (csKind==DYTools::_cs_preFsrDet) {
    if (!ia.needsFsrCorr()) {
      std::cout << "result is preFsrDet, but the FSR correction is switched off in inpArgs\n";
      res=0;
    }
    if (res) {
      HistoPair2D_t hpPreFsrDet(ia.resNameBase() + TString("hpPreFsrDet"));
      if (res) res=fsrCorrection_det(ia, hpEffRhoCorrYield, hpPreFsrDet);
      if (res) res=saveResult(ia,hpPreFsrDet,
			      ia.resNameBase() + TString("preFsrDet"));
      if (res) res=hp_fin.assign(hpPreFsrDet);
      if (!res) {
	std::cout << "failed to produce " << CrossSectionKindName(csKind) << "\n";
      }
    }
    return res;
  }
    
  // full space post-FSR result
  HistoPair2D_t hpPostFsrFullSp(ia.resNameBase() + TString("hpPostFsrFullSp"));
  if (ia.needsAccCorr()) {
    if (res) res=acceptanceCorrection(ia, hpEffRhoCorrYield, hpPostFsrFullSp);
    if (res) res=saveResult(ia,hpPostFsrFullSp,
			    ia.resNameBase() + TString("postFsrFullSp"));
  }
  else {
    std::cout << "result is in full space (not detector acceptance), but the acceptance correction is switched off in inpArgs\n";
    res=0;
  }

  if (csKind==DYTools::_cs_postFsr) {
    if (res) res=hp_fin.assign(hpPostFsrFullSp);
    if (!res) {
      std::cout << "failed to produce " << CrossSectionKindName(csKind) << "\n";
    }
    return res;
  }
  
  // pre-FSR in full space result
  HistoPair2D_t hpPreFsrFullSp(ia.resNameBase() + TString("hpPreFsrFullSp"));
  if (ia.needsFsrCorr()) {
    if (res) res=fsrCorrection_fullSpace(ia, hpPostFsrFullSp, hpPreFsrFullSp);
    if (res) res=saveResult(ia,hpPreFsrFullSp,
			    ia.resNameBase() + TString("preFsrFullSp"));
  }
  else {
    std::cout << "result is preFsr in full space, but the FSR correction is switched off in inpArgs\n";
    res=0;
  }

  if (csKind==DYTools::_cs_preFsr) {
    if (res) res=hp_fin.assign(hpPreFsrFullSp);
    if (!res) {
      std::cout << "failed to produce " << CrossSectionKindName(csKind) << "\n";
    }
    return res;
  }

  // all cross-section should be produced before this point
  std::cout << "CODE ERROR: failed to produce " << CrossSectionKindName(csKind) << "\n";
  return 0;
}

// ----------------------------------------------
// ----------------------------------------------

int calculateCS(const InputArgs_t &ia, const HistoPair2D_t &hp_ini, 
		DYTools::TCrossSectionKind_t csKind,
		HistoPair2D_t &hp_fin, CSResults_t &zpResult)
{
  zpResult.clear();

  // initial test
  {
    using namespace DYTools;
    int debugPrint=0;
    if (!checkCSKind(csKind,debugPrint,
		     8,
		     _cs_preFsr, _cs_preFsrDet,
		     _cs_postFsr,_cs_postFsrDet,
		     _cs_preFsrNorm, _cs_preFsrDetNorm,
		     _cs_postFsrNorm,_cs_postFsrDetNorm
		     )) {
      std::cout << "in calculateCSdistribution\n";
      return 0;
    }
  }

  // do the calculation

  DYTools::TCrossSectionKind_t csAbsKind=DYTools::getAbsCSKind(csKind);
  int res=calculateCSdistribution(ia,hp_ini,csAbsKind,hp_fin);
  if (!res) {
    std::cout << "calculateCS could not obtain the distribution\n";
    return 0;
  }

  if (DYTools::study2D) {
    hp_fin.eliminateNaNs();
  }

  // Divide by the luminosity
  hp_fin.scale(1/ia.inpMgr()->totalLumi());

  zpResult.assignZpeakCS(hp_fin);
  if (ia.silentMode()<2) std::cout << "normalization factor: " << zpResult << "\n";

  if (res && DYTools::isNormalizedCS(csKind)) {
    res=zpResult.doNormalize(hp_fin, ia.allNormErrorIsSyst());
    if (!res) {
      std::cout << "in calculateCS\n";
      return 0;
    }
    res=saveResult(ia,hp_fin,
		   ia.resNameBase() + TString("_final"));
  }

  return 1;
}

// ----------------------------------------------
// ----------------------------------------------

int calcVecOfCSdistributions(const InputArgs_t &ia,
			     const std::vector<TH2D*> &yieldIniV,
			     DYTools::TCrossSectionKind_t csKind,
			     std::vector<TH2D*> &yieldFinV)
{

  if (yieldIniV.size()==0) {
    std::cout << "calcVecOfCSdistributions: vector is empty\n";
    return 0;
  }
  yieldFinV.clear();
  yieldFinV.reserve(yieldIniV.size());

  int res=1;
  HistoPair2D_t hpIni("hp"), hpFin("hpres");
  CSResults_t csRes;
  for (unsigned int i=0; res && (i<yieldIniV.size()); ++i) {
    TH2D* hSrc=yieldIniV[i];
    res=hpIni.changeName(TString("hpIni_") + TString(hSrc->GetName()));
    hpIni.Reset();
    if (res) res=hpIni.add(hSrc,1);
    if (res) res=hpFin.changeName(TString("hpCS_") + TString(hSrc->GetName()));
    if (res) res=calculateCS(ia,hpIni,csKind,hpFin,csRes);
    TH2D *hRes=NULL;
    if (res) {
      hRes=Clone(hpFin.histo(),TString("hCS_") + TString(hSrc->GetName()));
    }
    if (res && !hRes) res=0;
    if (res && hRes) yieldFinV.push_back(hRes);
  }
  hpIni.clear(); hpFin.clear();
  if (!res) std::cout << "error in calcVecOfCSdistributions\n";
  return res;
}

// ----------------------------------------------
// ----------------------------------------------

