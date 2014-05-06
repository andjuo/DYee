#include "CSCovWorkFlags.hh"


// -----------------------------------------------------------
// -----------------------------------------------------------

int CSCovCalcFlags_t::finalizeFlags() {
  // Acc correction only in 1D
  fAccFSR *= (1-DYTools::study2D);
  fAccRnd *= (1-DYTools::study2D);

  if (this->calc_globalFSR() &&
      (fFsrFSR || fAccFSR || fEffFSR || fUnfFSR)) {
    std::cout << "since calc_globalFSR is on, "
	      << "individual FSR studies have to be switched off\n";
    return 0;
  }
  if (this->calc_globalPU() &&
      (fAccPU || fEffPU || fUnfPU || fFsrPU)) {
    std::cout << "since calc_globalPU is on, "
	      << "individual PU studies have to be switched off\n";
    return 0;
  }
  return 1;
}

// -----------------------------------------------------------
// -----------------------------------------------------------


int loadYieldCovMatrices(const TString &fnameBase,std::vector<TMatrixD*> &covs,
			 std::vector<TString> &labels, const WorkFlags_t &wf) {

  if (!wf.cf().doCalcYieldCov()) return 1;
  TString inpFileName=fnameBase;
  //if (wf.extraFileTag().Index("-yield")==-1) {
  //  inpFileName.ReplaceAll(".root","-yieldOnly.root");
  //}
  wf.adjustFName(inpFileName,_yield);
  TFile fin(inpFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << inpFileName << ">\n";
    return 0;
  }
  TMatrixD *ptr;
  if (wf.cf().calc_YieldStat()) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("YieldStat"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("signal stat (w/bkg)");
    }
  }
  if (wf.cf().calc_YieldSyst()) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("YieldSyst"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("signal syst (bkg)");
    }
  }
  if (wf.cf().calc_YieldStatDetailed()) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("YieldStatDetailed"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("signal stat");
    }
  }
  if (wf.cf().calc_YieldSystDetailed()) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("YieldSystDetailed"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("signal syst");
    }
  }
  if (wf.cf().calc_YieldEscale()) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("YieldEScale"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("signal EScale uncert.");
    }
  }

  fin.Close();

  std::cout << "loaded " << covs.size() << " entries from file <"
	    << fin.GetName() << ">\n";
  if (wf.cf().doCalcYieldCov() != int(covs.size())) {
    std::cout << " .. had to load " << wf.cf().doCalcYieldCov()
	      << " entries\n";
  }
  return 1;
}


// -----------------------------------------------------------


int loadUnfCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs,
		       std::vector<TString> &labels, const WorkFlags_t &wf) {

  if (!wf.cf().doCalcUnfCov()) return 1;

  TString inpFileName=fnameBase;
  //if (wf.extraFileTag().Index("-unf")==-1) {
  //  inpFileName.ReplaceAll(".root","-unfOnly.root");
  //}
  wf.adjustFName(inpFileName,_corrUnf);
  TFile fin(inpFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << inpFileName << ">\n";
    return 0;
  }
  TMatrixD *ptr;
  if (wf.cf().calc_UnfPU()) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("UnfPU"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("unf pile-up");
    }
  }
  if (wf.cf().calc_UnfFSR()) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("UnfFSR"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("unf FSR");
    }
  }
  if (wf.cf().calc_UnfRnd()) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("UnfRnd"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("unf stat");
    }
  }
  if (wf.cf().calc_UnfEScale()) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("UnfEScale"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("unf e-scale");
    }
  }
  fin.Close();

  std::cout << "loaded " << covs.size() << " entries from file <"
	    << fin.GetName() << ">\n";
  if (wf.cf().doCalcUnfCov() != int(covs.size())) {
    std::cout << "error the request was " << wf.cf().doCalcUnfCov()
	      << " matrices\n";
    for (unsigned int i=0; i<labels.size(); ++i) {
      std::cout << " - " << labels[i] << "\n";
    }
    return 0;
  }

  return 1;
}


// -----------------------------------------------------------


int loadEffCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs,
		       std::vector<TString> &labels, const WorkFlags_t &wf) {

  if (!wf.cf().doCalcEffCov()) return 1;

  TString inpFileName=fnameBase;
  //if (wf.extraFileTag().Index("-eff")==-1) {
  //  inpFileName.ReplaceAll(".root","-effOnly.root");
  //}
  wf.adjustFName(inpFileName,_corrEff);
  TFile fin(inpFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << inpFileName << ">\n";
    return 0;
  }
  TMatrixD *ptr;
  if (wf.cf().calc_EffPU()) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("EffPU"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("eff pile-up");
    }
  }
  if (wf.cf().calc_EffFSR()) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("EffFSR"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("eff FSR");
    }
  }
  if (wf.cf().calc_EffRnd()) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("EffRnd"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("eff stat");
    }
  }

  fin.Close();

  std::cout << "loaded " << covs.size() << " entries from file <"
	    << fin.GetName() << ">\n";
  if (wf.cf().doCalcEffCov() != int(covs.size())) {
    std::cout << "error: expected to load " << wf.cf().doCalcEffCov()
	      << " matrices\n";
    return 0;
  }
  return 1;
}


// -----------------------------------------------------------


int loadEsfCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs,
		       std::vector<TString> &labels, const WorkFlags_t &wf) {

  if (!wf.cf().doCalcESFCov()) return 1;

  TString inpFileName=fnameBase;
  //if (wf.extraFileTag().Index("-esf")==-1) {
  //  inpFileName.ReplaceAll(".root","-esfOnly.root");
  //}
  wf.adjustFName(inpFileName,_corrESF);
  TFile fin(inpFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << inpFileName << ">\n";
    return 0;
  }
  TMatrixD *ptr;
  if (wf.cf().calc_ESFtot()) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("ESFtot"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("ESF tot");
    }
  }
  else if (wf.cf().calc_ESFtotCheck()) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("ESFtotCheck"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("ESF tot");
    }
  }

  fin.Close();

  std::cout << "loaded " << covs.size() << " entries from file <"
	    << fin.GetName() << ">\n";
  if (wf.cf().doCalcESFCov() != int(covs.size())) {
    std::cout << "error: expected to load " << wf.cf().doCalcESFCov()
	      << " matrices\n";
    return 0;
  }
  return 1;
}


// -----------------------------------------------------------


int loadAccCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs,
		       std::vector<TString> &labels, const WorkFlags_t &wf) {

  if (!wf.cf().doCalcAccCov()) return 1;

  TString inpFileName=fnameBase;
  //inpFileName.ReplaceAll(".root","-accOnly.root");
  wf.adjustFName(inpFileName,_corrAcc);
  TFile fin(inpFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << inpFileName << ">\n";
    return 0;
  }
  TMatrixD *ptr;
  if (wf.cf().calc_AccFSR()) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("AccFSR"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("acc FSR");
    }
  }
  if (wf.cf().calc_AccRnd()) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("AccRnd"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("acc stat");
    }
  }

  fin.Close();

  std::cout << "loaded " << covs.size() << " entries from file <"
	    << fin.GetName() << ">\n";
  if (wf.cf().doCalcAccCov() != int(covs.size())) {
    std::cout << "error: expected to load " << wf.cf().doCalcAccCov()
	      << " matrices\n";
  }
  return 1;
}



// -----------------------------------------------------------


int loadFsrCovMatrices(const TString &fnameBase, std::vector<TMatrixD*> &covs,
		       std::vector<TString> &labels, const WorkFlags_t &wf) {

  if (!wf.cf().doCalcFSRCov()) return 1;

  TString inpFileName=fnameBase;
  //inpFileName.ReplaceAll(".root","-fsrUnfOnly.root");
  wf.adjustFName(inpFileName,_corrFSR);
  TFile fin(inpFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << inpFileName << ">\n";
    return 0;
  }
  TMatrixD *ptr;
  if (wf.cf().calc_FsrFSR()) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("FsrFSR"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("FSR FSR");
    }
  }
  if (wf.cf().calc_FsrRnd()) {
    ptr=(TMatrixD*)fin.Get(wf.fieldName("FsrRnd"));
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back("FSR stat");
    }
  }

  fin.Close();

  std::cout << "loaded " << covs.size() << " entries from file <"
	    << fin.GetName() << ">\n";
  if (wf.cf().doCalcFSRCov() != int(covs.size())) {
    std::cout << "error: expected to load " << wf.cf().doCalcFSRCov()
	      << " matrices\n";
  }
  return 1;
}



// -----------------------------------------------------------


int loadGlobalCovMatrices(const TString &fnameBase,
			  std::vector<TMatrixD*> &covs,
			  std::vector<TString> &labels,
			  const WorkFlags_t &wf) {

  if (!wf.cf().doCalcGlobalCov()) return 1;

  if (fnameBase.Length()==0) {
    std::cout << "loadGlobalCovMatrices warning: fnameBase.Length=0\n";
  }

  if (wf.showCSCov()==0) {
    std::cout << "loadGlobalCovMatrices: results are available only for "
	      << "the final cross section\n";
    return 0;
  }

  for (int i=0; i<3; ++i) {
    TString tag;
    int calc=0;
    TCorrCase_t corrCase=_corrNone;
    switch(i) {
    case 0: tag="puRndStudy"; calc=wf.cf().calc_globalPU();
            corrCase=_corrGlobalPU;
	    break;
    case 1: tag="fsrRndStudy"; calc=wf.cf().calc_globalFSR();
            corrCase=_corrGlobalFSR;
	    break;
    case 2: tag="fewzRndStudy"; calc=wf.cf().calc_globalFEWZ();
            corrCase=_corrGlobalFEWZ;
	    break;
    default:
      std::cout << "loadGlobalCovMatrices: the case i=" << i
		<< " is not ready\n";
      return 0;
    }
    if (!calc) continue;

    TString inpFileName;
    if (0) {
      // local file
      inpFileName=Form("csSyst-%s-%dD.root",tag.Data(),
		       DYTools::study2D+1);
    }
    else {
      inpFileName=fnameBase;
      wf.adjustFName(inpFileName,corrCase);
    }

    TFile fin(inpFileName,"read");
    if (!fin.IsOpen()) {
      std::cout << "failed to open a file <" << inpFileName << ">\n";
      return 0;
    }
    TString field= TString("covCS_") + tag;
    std::cout << "loading <" << field << "> from <" << inpFileName << ">\n";
    TMatrixD *ptr = (TMatrixD*)fin.Get(field);
    if (ptr) {
      covs.push_back(ptr);
      labels.push_back(tag);
    }
    else {
      std::cout << "failed to get <" << field << "> from <"
		<< fin.GetName() << ">\n";
      return 0;
    }
    fin.Close();
  }

  std::cout << "loaded " << covs.size() << " entries from global files\n";
  return 1;
}


// -----------------------------------------------------------
// -----------------------------------------------------------

TH2D *loadMainCSResult(int crossSection) {
  TString csFileName1D_CS="../../Results-DYee/root_files_reg/xsec/DY_j22_19712pb/xSec_preFsr_1D.root";
  TString csFileName1D_Count="../../Results-DYee/root_files_reg/xsec/DY_j22_19712pb/xSec_preFsr_1DpreFsr.root";
  TString csFileName2D_CS="../../Results-DYee/root_files_reg/xsec/DY_j22_19712pb/xSec_preFsrDet_2D.root";
  TString csFileName2D_Count="../../Results-DYee/root_files_reg/xsec/DY_j22_19712pb/xSec_preFsrDet_2DpreFsrDet.root";
  TString csFileName;
  TString fieldName="csYieldDDbkg";

  switch ( DYTools::study2D * 2 + ((crossSection) ? 1:0) ) {
  case 0: csFileName= csFileName1D_Count; fieldName="hpPreFsrFullSp";
    break;
  case 1: csFileName= csFileName1D_CS; break;
  case 2: csFileName= csFileName2D_Count; fieldName="hpPreFsrDet";
    break;
  case 3: csFileName= csFileName2D_CS; break;
  default:
   std::cout << "code error in loadMainCSResult\n";
   return NULL;
 }

  TFile fin(csFileName,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open the cross-section file <"
	      << fin.GetName() << ">\n";
    return NULL;
  }
  TH2D *h2=LoadHisto2D(fin,fieldName,"",1);
  fin.Close();
  if (!h2) {
    std::cout << "loadMainCSResult error\n";
  }
  printHisto(h2);
  return h2;
}

// -----------------------------------------------------------
// -----------------------------------------------------------

