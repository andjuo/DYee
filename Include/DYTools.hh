#ifndef DYTools_HH
#define DYTools_HH

//#define DYee7TeV
//#define DYee8TeV
#define DYee8TeV_reg


// Mass-Rapidity binning
// Non-default rapidity binnings are defined in separate files,
// included in DYTools.cc

#define useDefaultBinSet


// Check that the analysis is clearly defined
#if (defined DYee7TeV && (defined DYee8TeV || defined DYee8TeV_reg)) || (defined DYee8TeV && defined DYee8TeV_reg)
#   error define only 7TeV, 8TeV or 8TeV_reg analysis
#elif !(defined DYee7TeV || defined DYee8TeV || defined DYee8TeV_reg)
#   error analysis kind is not specified
#endif


// to ignore the rapidity range in 1D for indexing
// define preprocessor variable unlimited_Y_in_1D
// Note: this variable is not compatible with _MassBins_withFullOverflows
#ifdef useDefaultBinSet
//#  define unlimited_Y_in_1D
#endif

#include <iostream>
#include <math.h>
#include <assert.h>
#include <TMatrixD.h>
#include <TString.h>
#include <TH2D.h>


namespace DYTools {

  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  // ----------------------------

  // Declare mass binnings
  typedef enum { _MassBins_Undefined, _MassBins_default,
		 _MassBins_2011, _MassBins_2012,
		 _MassBins_test4, _MassBins_Zpeak,
		 _MassBins_noUnderflow, _MassBins_withFullOverflow,
		 _MassBins_bins5GeV,
		 _MassBins_bins100GeV,
		 _MassBins_finerMassRange
  }     TMassBinning_t;

  // ----------------------------

#ifdef DYee7TeV
  const int energy8TeV=0;
#endif
#ifdef DYee8TeV
  const int energy8TeV=1;
#endif
#ifdef DYee8TeV_reg
  const int energy8TeV=2;
#endif

  // ----------------------------

  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  // Global variable: analysis is either 1D or 2D
  // 1D analysis collapses the rapidity binning
  // The default choice of binnings for 1D and 2D
  //   assign correct values to
  //             nMassBins, massBinLimits, nYBins
  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  extern int study2D;
  extern int extendYRangeFor1D;
 
  extern TString analysisTag_USER;
  extern TString analysisTag; // final, to be used in the macros
  extern TString analysisTag_binning; // string describing the mass-y binning

  // For 8 TeV, we do unconventional 2.4 to match the muon channel
  extern const double electronEtaMax; // = (energy8TeV) ? 2.4 : 2.5;
  // Et cuts for the electrons
  extern double etMinLead; //  = 20; // GeV
  extern double etMinTrail; // = 10; // GeV

  // Luminosity definitions
  extern TString strLumiAtECMS;
  extern double lumiAtECMS;

  // Rapidity binning is different for different mass bins
  // Note: this implementation neglects underflow and overflow
  // in rapidity.
  const double yRangeMin =  0.0;
  const int _nBinsYLowMass  = (energy8TeV) ? 24 : 25;
  const int _nBinsYHighMass = (energy8TeV) ? 12 : 10;

  //
  // Define mass and rapidity ranges
  //
  extern DYTools::TMassBinning_t massBinningSet;
  extern int nMassBins;
  extern double *massBinLimits;
  extern int *nYBins;
  extern int nYBinsMax;
  extern double yRangeEdge;
  extern double yRangeMax;

  extern int nUnfoldingBinsMax;
  extern int nUnfoldingBins;

  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  int reset();
  int setup(int analysisIs2D);
  int checkTotalLumi(double lumi);

  // -----------------------------------------------
  // -----------------------------------------------

  // Global parameters, kinematics, etc
  const bool excludeEcalGap   = (energy8TeV) ? false : true; // For 8 TeV we are not excluding the transition region
  // ECAL constants, analysis-invariant
  const Double_t kECAL_MAX_ETA  = 2.5;
  const Double_t kECAL_GAP_MIDDLE = 1.479; // This is where the split barrel/endcap happens
  const Double_t kECAL_GAP_LOW  = 1.4442;
  const Double_t kECAL_GAP_HIGH = 1.566;

  // maximum number of canvas divisions in EventScaleFactors macros
  const int maxTnPCanvasDivisions=10;

  // luminosity
  const double lumiAtECMS_accuracy=1.;


  // ----------------------------
  // ----------------------------

  // Systematics modes for unfolding and acceptance
  typedef enum { SYST_MODE_FAILURE=0, NO_SYST, RESOLUTION_STUDY,
		 FSR_STUDY, FSR_RND_STUDY, FSR_5plus, FSR_5minus,
		 APPLY_ESCALE, ESCALE_RESIDUAL,
		 ESCALE_STUDY, ESCALE_STUDY_RND, UNREGRESSED_ENERGY,
		 LOWER_ET_CUT,
		 ESCALE_DIFF_0000, ESCALE_DIFF_0005, ESCALE_DIFF_0010, ESCALE_DIFF_0015, ESCALE_DIFF_0020,
		 NO_REWEIGHT, NO_REWEIGHT_PU, NO_REWEIGHT_FEWZ,
		 SYST_RND,
		 TAG_ID, TAG_PT,
		 PU_STUDY, PU_RND_STUDY, PILEUP_5plus, PILEUP_5minus,
		 UNREG_PU5plus, UNREG_PU5minus,
		 UNREG_FSR5plus, UNREG_FSR5minus,
		 UNREG_TagID, UNREG_TagPt } TSystematicsStudy_t;

  typedef enum { NORMAL_RUN, LOAD_DATA, DEBUG_RUN, DEBUG_LOAD, NO_RUN } TRunMode_t;

  inline int processData(TRunMode_t mode) { return ((mode!=LOAD_DATA) && (mode!=DEBUG_LOAD)) ? 1:0; }
  inline int loadData(TRunMode_t mode) { return ((mode==LOAD_DATA) || (mode==DEBUG_LOAD)) ? 1:0; }
  inline int isDebugMode(TRunMode_t mode) { return ((mode==DEBUG_RUN) || (mode==DEBUG_LOAD)) ? 1:0; }

  // Tag and probe fitting constants
  typedef enum {COUNTnCOUNT, COUNTnFIT, FITnFIT} TTnPMethod_t;
  typedef enum {RECO=0, ID=1, HLT=2, HLT_leg1=3, HLT_leg2=4, effNONE=99 } TEfficiencyKind_t;

  inline
  bool efficiencyIsHLT(TEfficiencyKind_t eff) {
    bool yes=false;
    switch(eff) {
    case HLT: case HLT_leg1: case HLT_leg2: yes=true; break;
    default: yes=false;
    }
    return yes;
  }

  // ----------------------------

  // asymmetric et cuts
  inline
  bool goodEtPair(double et1, double et2) {
    return ( (et1>etMinTrail) && (et2>etMinTrail) &&
	     ((et1>etMinLead) || (et2>etMinLead)) ) ? true : false;
  }

  // ----------------------------

  // generic bin idx finder
  inline
  int _findMassBin(double mass, int nMassBinsLoc, const double *massBinLimitsLoc){ 

    int result =-1;
    for(int ibin=0; ibin < nMassBinsLoc; ibin++){
      if( (mass >= massBinLimitsLoc[ibin]) && (mass < massBinLimitsLoc[ibin+1])) {
	result = ibin;
	break;
      }
    }
    return result;

  };

  // ----------------------------

  //
  // Define functions which should be used in the code
  //

  inline int validMass(double mass) {
    if (massBinningSet == _MassBins_withFullOverflow) return 1;
    return ((mass>=massBinLimits[0]) && (mass<=massBinLimits[nMassBins]));
  }

  // ----------------------------

  // find mass bin idx
  inline int findMassBin(double mass) {
    int result=-1;
    if (mass >= massBinLimits[nMassBins]) {
      if (massBinningSet == _MassBins_withFullOverflow) {
	result=nMassBins-1;
      }
    }
    else result=_findMassBin(mass,nMassBins,massBinLimits);
    return result;
  }

  // ----------------------------

  template<class Idx_t>
  inline double findMassBinCenter(Idx_t idx, int warn_on_error=1) {
    if ((idx<Idx_t(0)) || (idx>=Idx_t(nMassBins))) {
      if (warn_on_error) std::cout << "\n\tDYTools::findMassBinCenter(" << idx << ") - index error\n" << std::endl;
      return 0;
    }
    return 0.5*(massBinLimits[idx] + massBinLimits[idx+1]);
  }

  // ----------------------------

  inline
  int findYBin(int massBin, double y){
 
    int result = -1;
    if( massBin < 0 || massBin >= nMassBins) return result;

#ifdef unlimited_Y_in_1D
    if (massBinningSet == _MassBins_withFullOverflow) {
      std::cout << "findYBin: massBinningSet=_MassBins_withFullOverflows. unlimited_Y_in_1D cannot be set\n";
      return -1;
    }
    if (study2D==0) return 0; // for 1D case, include everything
#endif

    int nYBinsThisMassRange = nYBins[massBin];
    if (massBinningSet == _MassBins_withFullOverflow) {
      if (study2D==0) return 0;
      if (y >= yRangeEdge) return (nYBinsThisMassRange-1);
      nYBinsThisMassRange--;
    }
    if ( y < yRangeMin  ||  y > yRangeEdge ) return result;

    result = int( ( y - yRangeMin )* 1.000001*nYBinsThisMassRange/( yRangeEdge - yRangeMin ) );
    double dy=((yRangeEdge-yRangeMin)/double(nYBinsThisMassRange));
    //std::cout << "y=" << y << ", dy=" << dy  << ", binIdx=" << result << ", yRange=" << (result*dy) << ".." << (result*dy+dy) << "\n";
    if (result < 0 || result >= nYBinsThisMassRange) {
      std::cout << "y=" << y << ", dy=" << dy  << ", binIdx=" << result << ", yRange=" << (result*dy) << ".." << (result*dy+dy) << "\n";
      result=-1;
    }
   
    return result;
  }
 
  // ----------------------------

  inline int findAbsYBin(int massBin, double y){ return findYBin(massBin,fabs(y)); }
  inline int findYBin(double mass, double y) { return findYBin(findMassBin(mass),y); }
  inline int findAbsYBin(double mass, double y) { return findYBin(findMassBin(mass),fabs(y)); }

  // ----------------------------

  // return rapidity value at the Ybin center
  inline
  double findAbsYValue(int massBin, int yBin) {
    if (massBin>=nMassBins) {
      std::cout << "ERROR in findAbsYValue: massBin=" << massBin << ", nMassBins=" << nMassBins << "\n";
      return 0;
    }
    int ybinCount=nYBins[massBin];
    if (yBin>ybinCount) {
      std::cout << "ERROR in findAbsYValue: massBin=" << massBin << ", yBin=" << yBin << ", nYBins[massBin]=" << ybinCount << "\n";
      return 0;
    }
    double val=-1;
    if (massBinningSet!=_MassBins_withFullOverflow) {
      double yh=(yRangeMax-yRangeMin)/double(ybinCount);
      val= (yBin+0.5)*yh;
    }
    else {
      if (yBin==ybinCount-1) val=0.5*(yRangeMax+yRangeEdge);
      else {
	double yh=(yRangeEdge-yRangeMin)/double(ybinCount-1);
	val= (yBin+0.5)*yh;
      }
    }
    return val;
  }
 
  // ----------------------------

  // return rapidity range of Ybin
  inline
  void findAbsYValueRange(int massBin, int yBin, double &absYMin, double &absYMax) {
    absYMin=0.; absYMax=0.;
    if (massBin>=nMassBins) {
      std::cout << "ERROR in findAbsYValueRange: massBin=" << massBin << ", nMassBins=" << nMassBins << "\n";
      return;
    }
    int ybinCount=nYBins[massBin];
    if (yBin>ybinCount) {
      std::cout << "ERROR in findAbsYValueRange: massBin=" << massBin << ", yBin=" << yBin << ", nYBins[massBin]=" << ybinCount << "\n";
      return;
    }
    if (massBinningSet!=_MassBins_withFullOverflow) {
      double yh=(yRangeMax-yRangeMin)/double(ybinCount);
      absYMin=yBin*yh;
      absYMax=yBin*yh+yh;
    }
    else {
      if (yBin==ybinCount-1) {
	absYMin=yRangeEdge;
	absYMax=yRangeMax;
      }
      else {
	double yh=(yRangeEdge-yRangeMin)/double(ybinCount-1);
	absYMin=yBin*yh;
	absYMax=yBin*yh+yh;
      }
    }
    return;
  }
 
  // ----------------------------


  // This function finds a unique 1D index for (index_m, index_y) pair
  inline
  int findIndexFlat(int massBin, int yBin){
   
    int result = -1;
    if( massBin < 0 || massBin >= nMassBins || yBin < 0 || yBin >= nYBins[massBin] )
      return result;
   
    result = 0;
    for(int i=0; i< massBin; i++)
      result += nYBins[i];
   
    result += yBin;
   
    return result;
  }

  // ----------------------------

  // This function finds a unique 1D index for (index_m, index_y) pair
  inline
  int findIndexFlat(double mass, double y){
    int massBin=findMassBin(mass);
    int yBin=findAbsYBin(massBin,y);
    return findIndexFlat(massBin,yBin);
  }

  // ----------------------------

  inline bool validFlatIndex(int idx) {
    return ( idx != -1 && idx < nUnfoldingBins ) ? true : false;
  }

  // ----------------------------

  inline bool validFlatIndices(int idx1, int idx2) {
    return (validFlatIndex(idx1) && validFlatIndex(idx2)) ? true : false;
  }

  // ----------------------------

  //
  // Unfolding matrix binning
  // info: was getNumberOf2DBins
  inline
  int getTotalNumberOfBins(){
    int result = 0;
    for(int i=0; i<nMassBins; i++)
      result += nYBins[i];
    return result;
  }

  // ----------------------------

  // Largest number of Ybins
  inline
  int findMaxYBins(){
    int nYBMax=nYBins[0];
    for (int i=1; i<nMassBins; i++)
      if (nYBins[i]>nYBMax) nYBMax=nYBins[i];
    return nYBMax;
  }

  // ----------------------------

  // Bin limits in rapidity for a particular mass slice
  double *getYBinArray(int massBin);
  // ----------------------------

  // Rapidity bin limits
  // yBins should have dimensions (nMassBins, nYBinsMax+1)
  // The range is given by (im,iy)-(im,iy+1)
  TMatrixD* getYBinLimits();

  // ----------------------------
  // ----------------------------

  //
  // Define Et and Eta binning
  //
  typedef enum {ETBINS_UNDEFINED=-1, ETBINS1=1, ETBINS5, ETBINS6, ETBINS6syst,
		ETBINS6x, ETBINS6short,
		ETBINS6alt, ETBINS6altB, ETBINS7, ETBINS7alt, ETBINS7altB,
		ETBINS8, ETBINS9} TEtBinSet_t;
  const int nEtBins1 = 1;
  const double etBinLimits1[nEtBins1 + 1] =
    {10, 500};
  const int nEtBins5 = 5;
  const double etBinLimits5[nEtBins5 + 1] =
    {10, 20, 30, 40, 50, 500};
  const int nEtBins6 = 6;
  const double etBinLimits6[nEtBins6 + 1] =
    {10, 15, 20, 30, 40, 50, 500};
  const int nEtBins6x = 6;
  const double etBinLimits6x[nEtBins6x + 1] =
    {10, 15, 20, 30, 40, 50, 900};
  const double etBinLimits6short[nEtBins6 + 1] =
    {10, 15, 20, 30, 40, 50, 100};
  const int nEtBins6alt = 6;
  const double etBinLimits6alt[nEtBins6alt + 1] =
    {10, 20, 30, 40, 50, 75, 500};
  const int nEtBins6altB = 6;
  const double etBinLimits6altB[nEtBins6altB + 1] =
    {10, 20, 30, 40, 50, 60, 500};
  const int nEtBins7 = 7;
  const double etBinLimits7[nEtBins7 + 1] =
    {10, 15, 20, 30, 40, 50, 100, 500};
  const int nEtBins7alt = 7;
  const double etBinLimits7alt[nEtBins7alt + 1] =
    {10, 15, 20, 30, 40, 50, 75, 500};
  const int nEtBins7altB = 7;
  const double etBinLimits7altB[nEtBins7altB + 1] =
    {10, 20, 30, 40, 50, 75, 100, 500};
  const int nEtBins8 = 8;
  const double etBinLimits8[nEtBins8 + 1] =
    {10, 15, 20, 30, 40, 50, 75, 100, 500};
  const int nEtBins9 = 9;
  const double etBinLimits9[nEtBins9 + 1] =
    {10, 15, 20, 30, 40, 50, 75, 100, 125, 500};

  const int nEtBinsMax = nEtBins9;

  // ----------------------------

  // Get information about Et binning

  int getNEtBins(int binning);
  double *getEtBinLimits(int binning);

  // ----------------------------

  inline
  int findEtBin(double et, int binning){
   
    int result =-1;
    int n = getNEtBins(binning);
    double *limits = getEtBinLimits(binning);
    for(int ibin=0; ibin < n; ibin++){
      if( et >= limits[ibin] && et < limits[ibin+1]) {
	result = ibin;
	break;
      }
    }
    delete limits;
    return result;
  };

  // ----------------------------

  typedef enum {ETABINS_UNDEFINED=-1, ETABINS1=1, ETABINS2, ETABINS2Negs,
		ETABINS3, ETABINS3Negs, ETABINS5, ETABINS5egamma,
		ETABINS5nomerge,
		ETABINS5Negs, ETABINS5_max25, ETABINS4test, ETABINS4testNegs,
		ETABINS4alt, ETABINS4altNegs, ETABINS5alt, ETABINS5altNegs,
		ETABINS7, ETABINS8alt, ETABINS8altNegs, ETABINS9,
		ETABINS14} TEtaBinSet_t;

  const int nEtaBins1 = 1;
  const double etaBinLimits1[nEtBins1 + 1] =
    {0, 2.4000001};
  const int nEtaBins2 = 2;
  const double etaBinLimits2[nEtaBins2 + 1] =
    {0, 1.479, 2.4000001};
  const int nEtaBins3 = 3;
  const double etaBinLimits3[nEtaBins3 + 1] =
    {0, 1.479, 2.0, 2.4000001};
  const int nEtaBins3Negs = 6;
  const double etaBinLimits3Negs[nEtaBins3Negs + 1] =
    {-2.4000001, -2.0, -1.479, 0., 1.479, 2.0, 2.4000001};
  const int nEtaBins2Negs = 4;
  const double etaBinLimits2Negs[nEtaBins2Negs + 1] =
    {-2.400001, -1.479, 0, 1.479, 2.4000001};
  const int nEtaBins5 = 5;
  const double etaBinLimits5[nEtaBins5 + 1] =
    {0, 0.8, 1.4442, 1.566, 2.0, 2.400001 };
  const int nEtaBins5_max25 = 5;
  const double etaBinLimits5_max25[nEtaBins5_max25 + 1 ] =
    {0, 0.8, 1.4442, 1.566, 2.0, 2.500001 };
  const int nEtaBins4test = 4;
  const double etaBinLimits4test[nEtaBins4test + 1] =
    {0, 0.8, 1.479, 2.0, 2.400001 };
  const int nEtaBins4testNegs = 8;
  const double etaBinLimits4testNegs[nEtaBins4testNegs + 1] =
    {-2.40001, -2.0, -1.479, -0.8, 0, 0.8, 1.479, 2.0, 2.400001 };
  const int nEtaBins4alt = 4;
  const double etaBinLimits4alt[nEtaBins4alt + 1 ] =
    {0, 0.8, 1.479, 2.2, 2.400001 };
  const int nEtaBins4altNegs = 8;
  const double etaBinLimits4altNegs[nEtaBins4altNegs + 1 ] =
    {-2.400001, -2.2, -1.479, -0.8, 0, 0.8, 1.479, 2.2, 2.400001 };
  const int nEtaBins5alt = 5;
  const double etaBinLimits5alt[nEtaBins5alt + 1 ] =
    {0, 0.8, 1.4442, 1.566, 2.2, 2.400001 };
  const int nEtaBins5altNegs = 10;
  const double etaBinLimits5altNegs[nEtaBins5altNegs + 1 ] =
    {-2.400001, -2.2, -1.566, -1.4442, -0.8, 0, 0.8, 1.4442, 1.566, 2.2, 2.400001 };
  const int nEtaBins7 = 7;
  const double etaBinLimits7[nEtaBins7 + 1 ] =
    {0, 0.5, 1.0, 1.4442, 1.566, 1.9, 2.2, 2.40001 };
  const int nEtaBins8alt = 8;
  const double etaBinLimits8alt[nEtaBins8alt + 1 ] =
    {0, 0.4, 0.8, 1.0, 1.479, 1.8, 2.0, 2.2, 2.40001 };
  const int nEtaBins8altNegs = 16;
  const double etaBinLimits8altNegs[nEtaBins8altNegs + 1 ] =
    {-2.40001, -2.2, -2.0, -1.8, -1.479, -1.0, -0.8, -0.4, 0, 0.4, 0.8, 1.0, 1.479, 1.8, 2.0, 2.2, 2.4 };
  const int nEtaBins9 = 9;
  const double etaBinLimits9[nEtaBins9 + 1 ] =
    {0, 0.4, 0.8, 1.2, 1.4442, 1.566, 1.8, 2.0, 2.2, 2.40001 };
  const int nEtaBins14 = 14;
  const double etaBinLimits14[nEtaBins14 + 1 ] =
    {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4442, 1.566, 1.8, 2.0, 2.2, 2.3, 2.4, 2.500001 };

  const int nEtaBinsMax= nEtaBins5altNegs;

  // ----------------------------

  inline
  int mergeEtBins(int binning) {
    int res=0;
    switch(binning) {
    case ETABINS5:
    case ETABINS5_max25:
    case ETABINS5egamma:
      res=1;
      break;
    default: ;
    }
    return res;
  }

  // ----------------------------

  int getNEtaBins(int binning);
  double *getEtaBinLimits(int binning);
  int signedEtaBinning(int binning);

  // ----------------------------

  inline
  int findEtaBin(double eta, int binning){
    int result =-1;
    int n = getNEtaBins(binning);
    const double *limits = getEtaBinLimits(binning);
    if (!signedEtaBinning(binning) && (eta<0)) eta=-eta;
    for(int ibin=0; ibin < n; ibin++){
      if( eta >= limits[ibin] && eta < limits[ibin+1]) {
	result = ibin;
	break;
      }
    }
    delete limits;
    return result;
  }

  // ----------------------------

  // Primary vertices

  //const int nPVBinCount=7;
  //const double nPVLimits[nPVBinCount+1] = { 0., 5., 10., 15., 20., 25., 30., 100. };
  const int nPVBinCount=11;
  const double nPVLimits[nPVBinCount+1] = { 0.5, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5, 16.5, 20.5, 24.5, 40.5 };
  //  1., 3., 5., 7., 9., 11., 13., 15., 17., 21., 25., 40. };

  inline int findPUBin(int nPV) { return _findMassBin(double(nPV),nPVBinCount,nPVLimits); }
 
  // ---------------------
  //
  // Cross section types
  //
  typedef enum { _cs_None=0, _cs_spec,
		 _cs_preFsr, _cs_preFsrNorm,
		 _cs_preFsrDet, _cs_preFsrDetNorm,
		 _cs_preFsrDetErr, _cs_preFsrDetNormErr,
		 _cs_preFsrDetSystErr, _cs_preFsrDetNormSystErr,
		 _cs_postFsr, _cs_postFsrNorm,
		 _cs_postFsrDet, _cs_postFsrDetNorm } TCrossSectionKind_t;


  // ---------------------

  TCrossSectionKind_t getAbsCSKind(TCrossSectionKind_t normCS);
  int isAbsoluteCS(TCrossSectionKind_t cs);
  int isNormalizedCS(TCrossSectionKind_t cs);
  int isFullSpaceCS(TCrossSectionKind_t cs);
  int isPreFsrCS(TCrossSectionKind_t cs);

  // ---------------------

  //
  // Triggers vs run numbers
  //
  //  enum { UNDEF, REL38X, REL39X};
  typedef enum { DATA, MC} TDataKind_t;

 // ---------------------

  //
  // Barrel or endcap
  //

  inline
  bool isBarrel(double eta){
    bool result = false;
    const double etaBarrel=(excludeEcalGap) ? kECAL_GAP_LOW : kECAL_GAP_MIDDLE;
    if(fabs(eta) <= etaBarrel)
      result = true;
    return result;
  }

 // ---------------------

  inline
  bool isEndcap(double eta){
    bool result = false;
    const double etaEndcap=(excludeEcalGap) ? kECAL_GAP_HIGH : kECAL_GAP_MIDDLE;
    if(
#ifdef DYee7TeV
       fabs(eta) >= etaEndcap
#else
       fabs(eta) > etaEndcap   // kECAL_GAP_LOW=kECAL_GAP_HIGH for 8TeV, if gap is excluded
#endif
       && fabs(eta)<electronEtaMax )
      result = true;
    return result;
  }

 // ---------------------

  inline
  bool isEcalGap(double eta) {
    return ( fabs(eta) > kECAL_GAP_LOW && fabs(eta) < kECAL_GAP_HIGH );
  }

 // ---------------------

  inline
  bool goodEta(double eta) {
    bool result = true;
    // General requirement to be in range
    if (fabs(eta) >= electronEtaMax)
      result=false;

    if (excludeEcalGap && isEcalGap(eta))
      result=false;

    return result;
  }

 // ---------------------

  inline
  bool goodEtaPair(double eta1, double eta2) {
    return (goodEta(eta1) && goodEta(eta2));
  }

 // ---------------------

  inline int goodEtEtaPair(double et1, double eta1, double et2, double eta2) {
    return (goodEtPair(et1,et2) && goodEtaPair(eta1,eta2));
  }

 // ---------------------

};


// ------------------------------------------------------------------

const int retCodeOk=7;
const int retCodeError=13;
const int retCodeStop=11;

// ------------------------------------------------------------------
// ------------------------------------------------------------------

TString MassBinningName(DYTools::TMassBinning_t set);
TString CurrentMassBinningName();
TString SystematicsStudyName(DYTools::TSystematicsStudy_t study);
TString generateSystTag(DYTools::TSystematicsStudy_t systMode);
int usesUnregEnergy(DYTools::TSystematicsStudy_t systMode);
int isGlobalRndStudy(DYTools::TSystematicsStudy_t systMode);
TString RunModeName(DYTools::TRunMode_t study);
DYTools::TRunMode_t DebugInt2RunMode(int debug);
TString generateRunTag(DYTools::TRunMode_t runMode);
TString TnPMethodName(DYTools::TTnPMethod_t method);
TString EfficiencyKindName(DYTools::TEfficiencyKind_t kind);
TString EtBinSetName(DYTools::TEtBinSet_t set);
TString EtaBinSetName(DYTools::TEtaBinSet_t set);
TString DataKindName(DYTools::TDataKind_t kind);
DYTools::TSystematicsStudy_t DetermineSystematicsStudy(const TString &str);
DYTools::TTnPMethod_t DetermineTnPMethod(const TString &str);
DYTools::TEfficiencyKind_t DetermineEfficiencyKind(const TString &str);
DYTools::TEtBinSet_t DetermineEtBinSet(const TString& str);
DYTools::TEtaBinSet_t DetermineEtaBinSet(const TString& str);
DYTools::TDataKind_t DetermineDataKind(const TString &str);
TString CrossSectionKindName(DYTools::TCrossSectionKind_t kind);
DYTools::TCrossSectionKind_t DetermineCrossSectionKind(const TString &str);
TString CrossSectionKindLongName(DYTools::TCrossSectionKind_t kind);
void AdjustFileNameEnding(TString &fname,DYTools::TSystematicsStudy_t systMode,
			  int iSeed);

// ------------------------------------------------------------------
//     Printing
// ------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& out, DYTools::TMassBinning_t massBinning) { out << MassBinningName(massBinning); return out; }
inline std::ostream& operator<<(std::ostream& out, DYTools::TSystematicsStudy_t study) { out << SystematicsStudyName(study); return out; }
inline std::ostream& operator<<(std::ostream& out, DYTools::TRunMode_t mode) { out << RunModeName(mode); return out; }
inline std::ostream& operator<<(std::ostream& out, DYTools::TTnPMethod_t method) { out << TnPMethodName(method); return out; }
inline std::ostream& operator<<(std::ostream& out, DYTools::TEfficiencyKind_t kind) { out << EfficiencyKindName(kind); return out; }
inline std::ostream& operator<<(std::ostream& out, DYTools::TEtBinSet_t set) { out << EtBinSetName(set); return out; }
inline std::ostream& operator<<(std::ostream& out, DYTools::TEtaBinSet_t set) { out << EtaBinSetName(set); return out; }
inline std::ostream& operator<<(std::ostream& out, DYTools::TDataKind_t kind) { out << DataKindName(kind); return out; }
inline std::ostream& operator<<(std::ostream& out, DYTools::TCrossSectionKind_t kind) { out << CrossSectionKindName(kind); return out; }

// ------------------------------------------------------------------
// ------------------------------------------------------------------

namespace DYTools {

  void printRunMode(TRunMode_t mode);
  void printSystMode(TSystematicsStudy_t study);
  void printExecMode(TRunMode_t mode, TSystematicsStudy_t study);
  bool checkSystMode(DYTools::TSystematicsStudy_t mode, int debug_print, int nEntries, ...);
  bool checkCSKind(DYTools::TCrossSectionKind_t kind, int debug_print, int nEntries, ...);

};

// ------------------------------------------------------------------
// ------------------------------------------------------------------



#endif
