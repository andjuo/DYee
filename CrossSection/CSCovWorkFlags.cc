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

// increase the error on the efficiency scale factors
int TCovData_t::addESFsyst(const TString ver) {
  TString fileBase;
  TString refField="rhoSyst";
  if (ver==TString("20140525")) {
    fileBase="dir-ESFsyst-20140525/esf_syst_varEt5_" + DYTools::analysisTag;
  }
  else {
    std::cout << "TCovData_t::addESFsyst -- not ready for ver="
	      << ver << "\n";
    return 0;
  }

  if (covEsfV.size()!=1) {
    std::cout << "TCovData_t::addESFsyst: code assumes covEsfV.size=1\n";
    return 0;
  }

  TMatrixD* cov=covEsfV[0];
  TVectorD errs(cov->GetNrows());
  errs.Zero();

  int iUnfBin=(DYTools::study2D==0) ? 0 : DYTools::nYBinsMax;
  int sliceCount=(DYTools::study2D==0) ? 1 : (DYTools::nMassBins-1);

  for (int iSlice=0; iSlice<sliceCount; ++iSlice) {
    TString fname= fileBase;
    if (sliceCount==1) fname.Append("-mdf.root");
    else fname.Append(Form("_iM%d-mdf.root",iSlice+1));
    TFile file(fname,"read");
    if (!file.IsOpen()) {
      std::cout << "TCovData_t::addESFsyst - failed to open a file <"
		<< file.GetName() << ">\n";
      return 0;
    }
    TH1D *hRef=new TH1D(refField,refField,1,0.,1.);
    int res=loadHisto(file,&hRef,"");
    file.Close();
    if (!res) {
      std::cout << "error loading hRef\n";
      return 0;
    }

    std::cout << "File <" << file.GetName() << "> loaded\n";
    for (int ibin=1; ibin<=hRef->GetNbinsX(); ++ibin) {
      if (iUnfBin>=errs.GetNoElements()) {
	std::cout << "iUnfBin=" << iUnfBin << ", errs["
		  << errs.GetNoElements() << "]\n";
	std::cout << "Error in code during uncertainty accumulation\n";
	return 0;
      }
      errs[iUnfBin]=hRef->GetBinContent(ibin);
      iUnfBin++; // next entry
    }
    delete hRef;
  }
  std::cout << "iUnfBin=" << iUnfBin << "\n";
  if (iUnfBin!=errs.GetNoElements()) {
    std::cout << "accumulated iUnfBin=" << iUnfBin << " from needed "
	      << errs.GetNoElements() << " elements\n";
    return 0;
  }
  //std::cout << "errs="; errs.Print();

  TMatrixD* corr = corrFromCov(*cov);
  TVectorD newErr(errs);

  for (int i=0; i<errs.GetNoElements(); ++i) {
    double e1   =errs[i];
    double e2sqr=(*cov)(i,i);
    double eTot=sqrt( e1*e1 + e2sqr );
    std::cout << "i=" << i << ", ";
    std::cout << "old error=" << sqrt(e2sqr) << ", extraErr=" << e1;
    std::cout << ((e1>sqrt(e2sqr)) ? "(gt)" : "(le)");
    std::cout << ", totErr=" << eTot << "\n";
    newErr[i] = eTot;
  }

  TMatrixD* newCov= covFromCorr(*corr,newErr);
  if (!newCov) {
    std::cout << "failed to create new covariance\n";
    return 0;
  }

  // update the contents in the vector
  (*cov) = (*newCov);

  // clean-up. Do not remove cov, since it is still used
  delete corr;
  delete newCov;

  return 1;
}

// -----------------------------------------------------------

int TCovData_t::Write(TString subdir) const {
  if (subdir.Length()) {
    gDirectory->cd();
    gDirectory->mkdir(subdir);
    gDirectory->cd(subdir);
  }
  TString fieldName= TString("active");
  int res=writeFlagValues(fieldName,this->isActive);
  if (res) {
    for (unsigned int i=0; i < isActive.size(); ++i) {
      if (!isActive[i]) continue;
      const std::vector<TMatrixD*> *Cov=NULL;
      const std::vector<TString> *labels=NULL;
      switch(i) {
      case 0: Cov=&covYieldV; labels=&labelYieldV; break;
      case 1: Cov=&covUnfV; labels=&labelUnfV; break;
      case 2: Cov=&covEffV; labels=&labelEffV; break;
      case 3: Cov=&covEsfV; labels=&labelEsfV; break;
      case 4: Cov=&covAccV; labels=&labelAccV; break;
      case 5: Cov=&covFsrV; labels=&labelFsrV; break;
      case 6: Cov=&covGlobalV; labels=&labelGlobalV; break;
      default:
	std::cout << "TCovData::Write is not ready for set=" << i << "\n";
	return 0;
      }
      for (unsigned int ii=0; ii<Cov->size(); ++ii) {
	fieldName=(*labels)[ii];
	eliminateSeparationSigns(fieldName);
	(*Cov)[ii]->Write(fieldName);
      }
    }
  }
  if (subdir.Length()) gDirectory->cd();
  return 1;
}

// ----------------------------------------------------------------
// ----------------------------------------------------------------
