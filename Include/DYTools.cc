#include "../Include/DYTools.hh"


// -----------------------------------------------------------

namespace DYTools {
  int study2D=-1;
  int extendYRangeFor1D=1; // whether |ymax|=14 for 1D study

  TString analysisTag_USER;
  TString analysisTag="anTag_UNKNOWN";
  //TString study2Dstr;

  // For 8 TeV, we do unconventional 2.4 to match the muon channel
  const double electronEtaMax = (energy8TeV) ? 2.4 : 2.5;

  // Et cuts for the electrons
  double etMinLead  = 20;
  double etMinTrail = 10;

  // Luminosity definitions
  TString strLumiAtECMS="noLumiString";
  double lumiAtECMS=0.;

  DYTools::TMassBinning_t massBinningSet= _MassBins_Undefined;
  TString analysisTag_binning="mbUndefined";
  int nMassBins=0;
  double *massBinLimits=NULL;
  int *nYBins=NULL;
  int nYBinsMax=0;
  double yRangeEdge=0;
  double yRangeMax=0;

  int nUnfoldingBinsMax=0;
  int nUnfoldingBins=0;
};

// -----------------------------------------------------------
// -----------------------------------------------------------

#ifdef useDefaultBinSet  // defined in DYTools.hh
//
// Define mass and rapidity ranges
//

namespace DYTools {

// Constants that define binning in mass and rapidity
// Note: bin zero is underflow, overflow is neglected
const int _nMassBins2D = 7;
const double _massBinLimits2D[_nMassBins2D+1] =
  {0, // first bin is underflow
   20, 30, 45, 60, 120, 200, 1500
  }; // overflow is very unlikely, do not account for it
const int _nYBins2D[_nMassBins2D] =
  { _nBinsYLowMass,// underflow, binned like first mass bin
    _nBinsYLowMass, _nBinsYLowMass, _nBinsYLowMass, _nBinsYLowMass,
    _nBinsYLowMass,
    _nBinsYHighMass
    }; // overflow is neglected
const int _nYBinsMax2D=_nBinsYLowMass; // the largest division into Y bins

//
// Define default mass binning for 1D
//

// 2011 mass binning
const int _nMassBins2011 = 40;
const double _massBinLimits2011[_nMassBins2011+1] =
  {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76,
   81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141,
   150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440,
   510, 600, 1000, 1500}; // 40 bins
const int _nYBinsMax2011=1; // the largest division into Y bins
const int _nYBins2011[_nMassBins2011] = {
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};

// 2012 mass binning
const int _nMassBins2012 = 41;
const double _massBinLimits2012[_nMassBins2012+1] =
  {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76,
   81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141,
   150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440,
   510, 600, 1000, 1500, 2000 }; // 41 bin
const int _nYBinsMax2012=1; // the largest division into Y bins
const int _nYBins2012[_nMassBins2012] = {
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};

  const DYTools::TMassBinning_t _defined_massBinningSet= _MassBins_2012;
  const TString _defined_analysisTag_binning="";
  const double _extend_Y = 6.5; //additional space for the rapidity range in 1D

  const int _nMassBins1D= _nMassBins2012;
  const double * _massBinLimits1D= _massBinLimits2012;
  const int * _nYBins1D = _nYBins2012;
  const int _nYBinsMax1D= _nYBinsMax2012;


}; // namespace DYTools

// ------------------------------

#else

  // non standard ranges
  //#include "rangeDef_noUnderflow_inc.h"
  //#include "rangeDef_withFullOverflow_inc.h"
  //#include "rangeDef_massBinTest4_inc.h"
  //#include "rangeDef_ZpeakRegion_inc.h"
  //#include "rangeDef_bins100GeV_inc.h"
  #include "rangeDef_finerMassRange_inc.h"
#endif


// -----------------------------------------------------------
// -----------------------------------------------------------

namespace DYTools {

int assignMassRapidityValues() {
  massBinningSet= _defined_massBinningSet;
  analysisTag_binning= _defined_analysisTag_binning;
  const double *mblim=NULL;
  const int *ylim=NULL;
  if (DYTools::study2D==1) {
    nMassBins= _nMassBins2D;
    mblim= _massBinLimits2D;
    ylim = _nYBins2D;
    nYBinsMax= _nYBinsMax2D;
  }
  else if (DYTools::study2D==0) {
    nMassBins= _nMassBins1D;
    mblim= _massBinLimits1D;
    ylim = _nYBins1D;
    nYBinsMax= _nYBinsMax1D;
  }
  else {
    std::cout << "ERROR detected in assignMassRapidityValues: "
	      << "DYTools::study2D=" << DYTools::study2D << "\n";
    return 0;
  }

  double addY=(!study2D && extendYRangeFor1D) ? _extend_Y : 0;
  yRangeEdge= electronEtaMax + addY;
  yRangeMax = yRangeEdge;

  massBinLimits= new double[nMassBins+1];
  nYBins= new int[nMassBins];
  for (int i=0; i<=nMassBins; ++i) {
    massBinLimits[i] = mblim[i];
  }
  for (int i=0; i< nMassBins; ++i) {
    nYBins[i] = ylim[i];
  }
  return 1;
}

}; // namespace DYTools

// -----------------------------------------------------------
// -----------------------------------------------------------

namespace DYTools {
  int reset() {
    study2D=0;
    extendYRangeFor1D=1;

    analysisTag_USER="";

    strLumiAtECMS="noLumiString";
    lumiAtECMS=0.;
    return 1;
  }
};

// -----------------------------------------------------------

namespace DYTools {
  int setup(int analysisIs2D) {
    if (analysisIs2D==-111) {
      std::cout << "DYTools::setup(-111) is called. Nothing changed\n";
      return 1;
    }
    study2D=analysisIs2D;
    if (!assignMassRapidityValues()) {
      std::cout << "error in DYTools::setup(analysisIs2D=" << analysisIs2D << ")\n";
      return 0;
    }

#ifdef unlimited_Y_in_1D
    if (!study2D) analysisTag_USER.Append("-unlimY");
#endif

#ifdef DYee7TeV
    strLumiAtECMS= "4.8 fb^{-1} at #sqrt{s} = 7 TeV";
    lumiAtECMS= 4839.;
#endif
#ifdef DYee8TeV // not -regressed n-tuples
    strLumiAtECMS= "19.8 fb^{-1} at #sqrt(s) = 8 TeV";
    lumiAtECMS= 19789.;
#endif
#ifdef DYee8TeV_reg // energy-regressed files
    strLumiAtECMS= "19.7 fb^{-1} at #sqrt(s) = 8 TeV";
    lumiAtECMS= 19712.;
#endif

    TString study2Dstr=TString((study2D) ? "2D" : "1D") + analysisTag_binning;
    analysisTag= study2Dstr + analysisTag_USER;
    nUnfoldingBinsMax= nMassBins * nYBinsMax;
    nUnfoldingBins= nMassBins;
    if (study2D) {
      nUnfoldingBins=0;
      for (int im=0; im<nMassBins; ++im) {
	nUnfoldingBins += nYBins[im];
      }
    }

    return 1;
  }
};

// -----------------------------------------------------------
// -----------------------------------------------------------

// Bin limits in rapidity for a particular mass slice

namespace DYTools {

  double *getYBinArray(int massBin){
    double *result = 0;
    if( massBin < 0 || massBin >= nMassBins ) {
      return result;
    }
    int nYBinsThisSlice = nYBins[massBin];
    result = new double[nYBinsThisSlice+1];

    double delta = (yRangeMax - yRangeMin)/double(nYBinsThisSlice);
    if (massBinningSet== _MassBins_withFullOverflow) {
      // yRangeEdge <> yRangeMax in this case
      delta = (yRangeEdge - yRangeMin)/double (nYBinsThisSlice-1);
    }
    for(int i=0; i<nYBinsThisSlice; i++){
      result[i] = yRangeMin + i * delta;
    }

    result[nYBinsThisSlice] = yRangeMax;
    return result;
  }

  // ----------------------------

  // Rapidity bin limits
  // yBins should have dimensions (nMassBins, nYBinsMax+1)
  // The range is given by (im,iy)-(im,iy+1)
  TMatrixD* getYBinLimits() {
    TMatrixD *yBins=new TMatrixD(nMassBins, nYBinsMax+1);
    (*yBins).Zero();
    for (int im=0; im<nMassBins; ++im) {
      double delta = (yRangeMax - yRangeMin)/double(nYBins[im]);
      if (massBinningSet== _MassBins_withFullOverflow) {
	// yRangeEdge <> yRangeMax in this case
	delta = (yRangeEdge - yRangeMin)/double (nYBins[im]-1);
      }
      for (int iy=0; iy<nYBins[im]; ++iy) {
	(*yBins)(im,iy) = yRangeMin + iy * delta;
      }
      (*yBins)(im,nYBins[im]) = yRangeMax;
    }
    return yBins;
  }

};

// -----------------------------------------------------------
// -----------------------------------------------------------

namespace DYTools {

  // ----------------------------

  int getNEtBins(int binning){
    int n=0;
    switch(binning) {
    case ETBINS_UNDEFINED: n = 0; break;
    case ETBINS1: n = nEtBins1; break;
    case ETBINS5: n = nEtBins5; break;
    case ETBINS6syst: // same as ETBINS6, avoid built-in systematics
    case ETBINS6: n = nEtBins6; break;
    case ETBINS6x: n = nEtBins6x; break;
    case ETBINS6short: n = nEtBins6; break;
    case ETBINS6alt: n = nEtBins6alt; break;
    case ETBINS6altB: n = nEtBins6altB; break;
    case ETBINS7: n = nEtBins7; break;
    case ETBINS7alt: n = nEtBins7alt; break;
    case ETBINS7altB: n = nEtBins7altB; break;
    case ETBINS8: n = nEtBins8; break;
    case ETBINS9: n = nEtBins9; break;
    default:
      printf("ERROR: unknown binning requested\n");
      n=0;
    }
    return n;
  }

  // ----------------------------

  double *getEtBinLimits(int binning){
    int n = getNEtBins(binning);
    double *limitsOut = new double[n+1];
    const double *limits = NULL;
    switch(binning) {
    case ETBINS1: limits=etBinLimits1; break;
    case ETBINS5: limits=etBinLimits5; break;
    case ETBINS6syst: // same as ETBINS6, avoid built-in systematics
    case ETBINS6: limits=etBinLimits6; break;
    case ETBINS6x: limits=etBinLimits6x; break;
    case ETBINS6short: limits=etBinLimits6short; break;
    case ETBINS6alt: limits=etBinLimits6alt; break;
    case ETBINS6altB: limits=etBinLimits6altB; break;
    case ETBINS7: limits=etBinLimits7; break;
    case ETBINS7alt: limits=etBinLimits7alt; break;
    case ETBINS7altB: limits=etBinLimits7altB; break;
    case ETBINS8: limits=etBinLimits8; break;
    case ETBINS9: limits=etBinLimits9; break;
    default:
      printf("ERROR: unknown/undefined binning requested\n");
      assert(0);
    }
    for(int i=0; i<=n; i++)
      limitsOut[i] = limits[i];
   
    return limitsOut;
  }

  // ----------------------------
  // ----------------------------

  int getNEtaBins(int binning){
    int n=0;
    switch(binning) {
    case ETABINS_UNDEFINED: n = 0; break;
    case ETABINS1: n = nEtaBins1; break;
    case ETABINS2: n = nEtaBins2; break;
    case ETABINS2Negs: n = nEtaBins2Negs; break;
    case ETABINS3: n = nEtaBins3; break;
    case ETABINS3Negs: n = nEtaBins3Negs; break;
    case ETABINS5egamma: // identical, although eff is with |eta|<2.5
    case ETABINS5nomerge:
    case ETABINS5: n = nEtaBins5; break;
    case ETABINS5_max25: n = nEtaBins5_max25; break;
    case ETABINS4test: n = nEtaBins4test; break;
    case ETABINS4testNegs: n = nEtaBins4testNegs; break;
    case ETABINS4alt: n = nEtaBins4alt; break;
    case ETABINS4altNegs: n = nEtaBins4altNegs; break;
    case ETABINS5alt: n = nEtaBins5alt; break;
    case ETABINS5altNegs: n = nEtaBins5altNegs; break;
    case ETABINS7: n = nEtaBins7; break;
    case ETABINS8alt: n = nEtaBins8alt; break;
    case ETABINS8altNegs: n = nEtaBins8altNegs; break;
    case ETABINS9: n = nEtaBins9; break;
    case ETABINS14: n = nEtaBins14; break;
    default:
      printf("ERROR: unknown binning requested\n");
      assert(0);
      n=0;
    }
    return n;
  }

  // ----------------------------

  double *getEtaBinLimits(int binning){
    int n = getNEtaBins(binning);
    double *limitsOut = new double[n+1];
    const double *limits = NULL;
    switch(binning) {
    case ETABINS1: limits = etaBinLimits1; break;
    case ETABINS2: limits = etaBinLimits2; break;
    case ETABINS2Negs: limits = etaBinLimits2Negs; break;
    case ETABINS3: limits = etaBinLimits3; break;
    case ETABINS3Negs: limits = etaBinLimits3Negs; break;
    case ETABINS5egamma: // identical, although eff is with |eta|<2.5
    case ETABINS5nomerge:
    case ETABINS5: limits = etaBinLimits5; break;
    case ETABINS5_max25: limits = etaBinLimits5_max25; break;
    case ETABINS4test: limits = etaBinLimits4test; break;
    case ETABINS4testNegs: limits = etaBinLimits4testNegs; break;
    case ETABINS4alt: limits = etaBinLimits4alt; break;
    case ETABINS4altNegs: limits = etaBinLimits4altNegs; break;
    case ETABINS5alt: limits = etaBinLimits5alt; break;
    case ETABINS5altNegs: limits = etaBinLimits5altNegs; break;
    case ETABINS7: limits = etaBinLimits7; break;
    case ETABINS8alt: limits = etaBinLimits8alt; break;
    case ETABINS8altNegs: limits = etaBinLimits8altNegs; break;
    case ETABINS9: limits = etaBinLimits9; break;
    case ETABINS14: limits = etaBinLimits14; break;
    default:
      printf("ERROR: unknown/undefined binning requested\n");
      assert(0);
      n=0;
    }
    for (int i=0; i<=n; ++i) {
      limitsOut[i] = limits[i];
    }
    return limitsOut;
  }

  // ----------------------------

  int signedEtaBinning(int binning) {
    int yes=0;
    switch(binning) {
    case ETABINS1:
    case ETABINS2:
    case ETABINS3:
    case ETABINS5egamma: // identical, although eff is with |eta|<2.5
    case ETABINS5nomerge:
    case ETABINS5:
    case ETABINS5_max25:
    case ETABINS4test:
    case ETABINS4alt:
    case ETABINS5alt:
    case ETABINS7:
    case ETABINS8alt:
    case ETABINS9:
    case ETABINS14:
      yes=0;
      break;
    case ETABINS2Negs:
    case ETABINS3Negs:
    case ETABINS4testNegs:
    case ETABINS4altNegs:
    case ETABINS5altNegs:
    case ETABINS8altNegs:
      yes=1;
      break;
    default:
      printf("ERROR: unknown/undefined binning requested\n");
      assert(0);
    }
    return yes;
  }

  // ----------------------------

};

// -----------------------------------------------------------
// -----------------------------------------------------------

namespace DYTools {

  // ---------------------

  TCrossSectionKind_t getAbsCSKind(TCrossSectionKind_t normCS) {
    TCrossSectionKind_t kind = _cs_None;
    switch(normCS) {
    case _cs_spec: kind=_cs_spec; break;
    case _cs_preFsr:
    case _cs_preFsrNorm: kind=_cs_preFsr; break;
    case _cs_preFsrDet:
    case _cs_preFsrDetNorm:
    case _cs_preFsrDetErr:
    case _cs_preFsrDetNormErr:
    case _cs_preFsrDetSystErr:
    case _cs_preFsrDetNormSystErr: kind=_cs_preFsrDet; break;
    case _cs_postFsr:
    case _cs_postFsrNorm: kind=_cs_postFsr; break;
    case _cs_postFsrDet:
    case _cs_postFsrDetNorm: kind=_cs_postFsrDet; break;
    default:
      std::cout << "unhandled case in getAbsCSKind\n";
    }
    return kind;
  }

  // ---------------------

  int isAbsoluteCS(TCrossSectionKind_t cs) {
    int yes=0;
    switch(cs) {
    case _cs_spec:
    case _cs_preFsr:
    case _cs_preFsrDet:
    case _cs_postFsr:
    case _cs_postFsrDet:
      yes=1;
      break;
    default: ;
    }
    return yes;
  }

  // ---------------------

  int isNormalizedCS(TCrossSectionKind_t cs) {
    int yes=0;
    switch(cs) {
    case _cs_preFsrNorm:
    case _cs_preFsrDetNorm:
    case _cs_postFsrNorm:
    case _cs_postFsrDetNorm:
      yes=1;
      break;
    default: ;
    }
    return yes;
  }

  // ---------------------

  int isFullSpaceCS(TCrossSectionKind_t cs) {
    int yes=0;
    switch(cs) {
    case _cs_preFsr:
    case _cs_preFsrNorm:
    case _cs_postFsr:
    case _cs_postFsrNorm:
      yes=1;
      break;
    default: ;
    }
    return yes;
  }

  // ---------------------

  int isPreFsrCS(TCrossSectionKind_t cs) {
    int yes=0;
    switch(cs) {
    case _cs_spec: yes=-1; break;
    case _cs_preFsr:
    case _cs_preFsrNorm:
    case _cs_preFsrDet:
    case _cs_preFsrDetNorm:
    case _cs_preFsrDetErr:
    case _cs_preFsrDetNormErr:
    case _cs_preFsrDetSystErr:
    case _cs_preFsrDetNormSystErr: yes=1; break;
    case _cs_postFsr:
    case _cs_postFsrNorm:
    case _cs_postFsrDet:
    case _cs_postFsrDetNorm: yes=0; break;
    default:
      std::cout << "unhandled case in isPreFsrCS\n";
    }
    return yes;
  }

  // ---------------------


};

// -----------------------------------------------------------
// -----------------------------------------------------------

// --------------------------------------------------------

namespace DYTools {

  int checkTotalLumi(double lumi) {
    if (fabs(lumi-DYTools::lumiAtECMS)>lumiAtECMS_accuracy) {
      std::cout << "checkTotalLumi detected mismatch: lumi=" << lumi
		<< ", DYTools::lumiAtECMS=" << DYTools::lumiAtECMS
		<< " +/- " << DYTools::lumiAtECMS_accuracy << std::endl;
      return 0;
    }
    return 1;
  }

};

// -----------------------------------------------------------
// -----------------------------------------------------------

#include "DYTools_inc.hh"

// -----------------------------------------------------------
// -----------------------------------------------------------
