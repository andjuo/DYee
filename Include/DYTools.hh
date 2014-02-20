#ifndef DYTools_HH
#define DYTools_HH

//#define DYee7TeV
//#define DYee8TeV
#define DYee8TeV_reg

#if (defined DYee7TeV && (defined DYee8TeV || defined DYee8TeV_reg)) || (defined DYee8TeV && defined DYee8TeV_reg)
#   error define only 7TeV, 8TeV or 8TeV_reg analysis
#elif !(defined DYee7TeV || defined DYee8TeV || defined DYee8TeV_reg) 
#   error analysis kind is not specified
#endif


#define useDefaultBinSet

// to ignore the rapidity range in 1D for indexing
// define preprocessor variable unlimited_Y_in_1D
// Note: this variable is not compatible with _MassBins_withFullOverflows
#ifdef useDefaultBinSet
//#  define unlimited_Y_in_1D
#endif

#include <iostream>
#include <math.h>
#include <assert.h>
//#include <TEfficiency.h>
#include <TMatrixD.h>
#include <TString.h>
#include <TH2D.h>

#include "../Include/EWKAnaDefs.hh"


namespace DYTools {

  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  // Global variable: analysis is either 1D or 2D
  // 1D analysis collapses the rapidity binning
  // The default choice of binnings for 1D and 2D
  //   assign correct values to 
  //             nMassBins, massBinLimits, nYBins
  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

  const int study2D=0;
  const int extendYRangeFor1D=1; // whether |ymax|=14 for 1D study
#ifdef unlimited_Y_in_1D
  const TString analysisTag_USER=(!study2D) ? "-unlimY" : "";
#else
  const TString analysisTag_USER=""; // extra name to differentiate the analysis files
#endif

  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
  // ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !

#ifdef DYee7TeV
  const int energy8TeV=0;
#endif
#ifdef DYee8TeV 
  const int energy8TeV=1;
#endif
#ifdef DYee8TeV_reg
  const int energy8TeV=2;
#endif

  // Global parameters, kinematics, etc
  const bool excludeEcalGap   = (energy8TeV) ? false : true; // For 8 TeV we are not excluding the transition region
  // ECAL constants, analysis-invariant
  const Double_t kECAL_MAX_ETA  = 2.5;
  const Double_t kECAL_GAP_MIDDLE = 1.479; // This is where the split barrel/endcap happens
  const Double_t kECAL_GAP_LOW  = 1.4442;
  const Double_t kECAL_GAP_HIGH = 1.566;
  // For 8 TeV, we do unconventional 2.4 to match the muon channel
  const double electronEtaMax = (energy8TeV) ? 2.4 : 2.5;
  // Et cuts for the electrons
  const double etMinLead  = 20;
  const double etMinTrail = 10;

  const int maxTnPCanvasDivisions=10; // maximum number of canvas divisions in EventScaleFactors macros

#ifndef DYee8TeV_reg // luminosity for not-regressed files or 7TeV analysis
  const TString strLumiAtECMS=(energy8TeV) ? "19.8 fb^{-1} at #sqrt(s) = 8 TeV" : "4.8 fb^{-1} at #sqrt{s} = 7 TeV";
  const double lumiAtECMS=(energy8TeV) ? 19789. : 4839.;
#else // luminosity for regressed files or 7TeV analysis
  const TString strLumiAtECMS=(energy8TeV==2) ? "19.7 fb^{-1} at #sqrt(s) = 8 TeV" : "4.8 fb^{-1} at #sqrt{s} = 7 TeV";
  const double lumiAtECMS=(energy8TeV==2) ? 19712. : 4839.;
#endif

  const double lumiAtECMS_accuracy=1.;


  // asymmetric et cuts
  inline
  bool goodEtPair(double et1, double et2) {
    return ( (et1>etMinTrail) && (et2>etMinTrail) &&
	     ((et1>etMinLead) || (et2>etMinLead)) ) ? true : false;
  }

  // Rapidity binning is different for different mass bins
  // Note: this implementation neglects underflow and overflow
  // in rapidity.
  const double yRangeMin =  0.0;
  const double _yRangeMax_base =  electronEtaMax + ((!study2D && extendYRangeFor1D) ? ((energy8TeV) ? 6.5:6.5) : 0);
  const double _yRangeEdge_base = _yRangeMax_base;
  const int _nBinsYLowMass  = (energy8TeV) ? 24 : 25;
  const int _nBinsYHighMass = (energy8TeV) ? 12 : 10;

  // ----------------------------

  // Declare mass binnings
  typedef enum { _MassBins_Undefined, _MassBins_default, 
		 _MassBins_2011, _MassBins_2012,
		 _MassBins_test4, _MassBins_Zpeak,
		 _MassBins_noUnderflow, _MassBins_withFullOverflow 
  }     TMassBinning_t;

  // ----------------------------

  //
  // Define mass and rapidity ranges
  //

#ifdef useDefaultBinSet

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

  // 2011 mass binning
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

  const TString analysisTag_binning="";
  const DYTools::TMassBinning_t massBinningSet= _MassBins_2012; //_MassBins_default;
  const int _nMassBins1D=_nMassBins2012;
  const double *_massBinLimits1D=_massBinLimits2012;
  const int *_nYBins1D=_nYBins2012;
  const int _nYBinsMax1D=_nYBinsMax2012;
  const double _yRangeEdge=_yRangeEdge_base;
  const double _yRangeMax=_yRangeMax_base; 

#else
  // non standard ranges
  //#include "rangeDef_noUnderflow_inc.h"
  #include "rangeDef_withFullOverflow_inc.h"
  //#include "rangeDef_massBinTest4_inc.h"
  //#include "rangeDef_ZpeakRegion_inc.h"

#endif

  // variables to be used in the analysis

  const int nMassBins=(study2D) ? _nMassBins2D : _nMassBins1D;
  const double *massBinLimits=(study2D) ? _massBinLimits2D : _massBinLimits1D;
  const int *nYBins=(study2D) ? _nYBins2D : _nYBins1D;
  const int nYBinsMax=(study2D) ? _nYBinsMax2D : _nYBinsMax1D;
  const double yRangeEdge=_yRangeEdge;
  const double yRangeMax=_yRangeMax; 

  // ----------------------------

  // Systematics modes for unfolding and acceptance 
  typedef enum { SYST_MODE_FAILURE=0, NO_SYST, RESOLUTION_STUDY, FSR_STUDY, ESCALE_RESIDUAL, ESCALE_STUDY, ESCALE_STUDY_RND, SYST_RND } TSystematicsStudy_t;

  typedef enum { NORMAL_RUN, LOAD_DATA, DEBUG_RUN, DEBUG_LOAD } TRunMode_t;

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

  const TString study2Dstr=TString((study2D) ? "2D" : "1D") + analysisTag_binning;
  const TString analysisTag=study2Dstr + analysisTag_USER;
  const int nUnfoldingBinsMax= nMassBins * nYBinsMax;
  const int nUnfoldingBins= (study2D) ? 
    ((nMassBins-2) * nYBinsMax + 
     nYBins[nMassBins-2] + nYBins[nMassBins-1]) : nMassBins;


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

  // 
  // Define functions which should be used in the code
  //

  inline int validMass(double mass) { 
    if (massBinningSet == _MassBins_withFullOverflow) return 1;
    return ((mass>=massBinLimits[0]) && (mass<=massBinLimits[nMassBins])); 
  }

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


  template<class Idx_t>
  inline double findMassBinCenter(Idx_t idx, int warn_on_error=1) { 
    if ((idx<Idx_t(0)) || (idx>=Idx_t(nMassBins))) {
      if (warn_on_error) std::cout << "\n\tDYTools::findMassBinCenter(" << idx << ") - index error\n" << std::endl;
      return 0;
    }
    return 0.5*(massBinLimits[idx] + massBinLimits[idx+1]);
  }

  inline 
  int findYBin(int massBin, double y){
  
    int result = -1;
    if( massBin < 0 || massBin > nMassBins) return result;

#ifdef unlimited_Y_in_1D
    if (massBinningSet == _MassBins_withFullOverflow) {
      std::cout << "findYBin: massBinningSet=_MassBins_withFullOverflows. unlimited_Y_in_1D cannot be set\n";
      return -1;
    }
    if (study2D==0) return 0; // for 1D case, include everything
#endif

    int nYBinsThisMassRange = nYBins[massBin];
    if (massBinningSet == _MassBins_withFullOverflow) {
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
  
  inline int findAbsYBin(int massBin, double y){ return findYBin(massBin,fabs(y)); }
  inline int findYBin(double mass, double y) { return findYBin(findMassBin(mass),y); }
  inline int findAbsYBin(double mass, double y) { return findYBin(findMassBin(mass),fabs(y)); }

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
  

  // This function finds a unique 1D index for (index_m, index_y) pair
  inline 
  int findIndexFlat(int massBin, int yBin){
    
    int result = -1;
    if( massBin < 0 || massBin > nMassBins || yBin < 0 || yBin > nYBins[massBin] )
      return result;
    
    result = 0;
    for(int i=0; i< massBin; i++)
      result += nYBins[i];
    
    result += yBin;
    
    return result;
  }

  // This function finds a unique 1D index for (index_m, index_y) pair
  inline 
  int findIndexFlat(double mass, double y){
    int massBin=findMassBin(mass);
    int yBin=findAbsYBin(massBin,y);
    return findIndexFlat(massBin,yBin);
  }


  inline bool validFlatIndex(int idx) {
    return ( idx != -1 && idx < nUnfoldingBins ) ? true : false;
  }

  inline bool validFlatIndices(int idx1, int idx2) {
    return (validFlatIndex(idx1) && validFlatIndex(idx2)) ? true : false;
  }



  // 
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

  // Largest number of Ybins
  inline 
  int findMaxYBins(){
    int nYBMax=nYBins[0];
    for (int i=1; i<nMassBins; i++)
      if (nYBins[i]>nYBMax) nYBMax=nYBins[i];
    return nYBMax;
  }

  // Bin limits in rapidity for a particular mass slice
  inline 
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


  // Rapidity bin limits
  // yBins should have dimensions (nMassBins, nYBinsMax+1)
  // The range is given by (im,iy)-(im,iy+1)
  inline
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



  // NEEDS clean up
//   // Note: this HAS TO BE the total number of all 2D bins, that is
//   // the sum of the contents of the nYBins array above
//   const int nUnfoldingBins = 160;
  


  //
  // Define single electron Pt binning
  //
  /*
  const int nPtBins = 5;
  const double ptBinLimits[nPtBins+1] = 
    {10, 20, 30, 40, 50, 500};
  
  int findPtBin(double pt){
    
    int result =-1;
    for(int ibin=0; ibin < nPtBins; ibin++){
      if( pt >= ptBinLimits[ibin] && pt < ptBinLimits[ibin+1]) {
	result = ibin;
	break;
      }
    }
    
    return result;
  };
  */

  //
  // Define Et and Eta binning
  //
  typedef enum {ETBINS_UNDEFINED=-1, ETBINS1=1, ETBINS5, ETBINS6, ETBINS6syst, ETBINS6alt, ETBINS6altB, ETBINS7, ETBINS7alt, ETBINS7altB, ETBINS8, ETBINS9} TEtBinSet_t;
  const int nEtBins1 = 1;
  const double etBinLimits1[nEtBins1 + 1] = 
    {10, 500};
  const int nEtBins5 = 5;
  const double etBinLimits5[nEtBins5 + 1] = 
    {10, 20, 30, 40, 50, 500};
  const int nEtBins6 = 6;
  const double etBinLimits6[nEtBins6 + 1] = 
    {10, 15, 20, 30, 40, 50, 500};
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


  inline 
  int getNEtBins(int binning){
    int n=0;
    switch(binning) {
    case ETBINS_UNDEFINED: n = 0; break;
    case ETBINS1: n = nEtBins1; break;
    case ETBINS5: n = nEtBins5; break;
    case ETBINS6syst: // same as ETBINS6, avoid built-in systematics
    case ETBINS6: n = nEtBins6; break;
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

  inline 
  double *getEtBinLimits(int binning){
    int n = getNEtBins(binning);
    double *limitsOut = new double[n+1];
    const double *limits = NULL;
    switch(binning) {
    case ETBINS1: limits=etBinLimits1; break;
    case ETBINS5: limits=etBinLimits5; break;
    case ETBINS6syst: // same as ETBINS6, avoid built-in systematics
    case ETBINS6: limits=etBinLimits6; break;
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

  typedef enum {ETABINS_UNDEFINED=-1, ETABINS1=1, ETABINS2, ETABINS2Negs, ETABINS3, ETABINS3Negs, ETABINS5, ETABINS5egamma, ETABINS5Negs, ETABINS5_max25, ETABINS4test, ETABINS4testNegs,  ETABINS4alt, ETABINS4altNegs, ETABINS5alt, ETABINS5altNegs, ETABINS8alt, ETABINS8altNegs, ETABINS14} TEtaBinSet_t;
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
  const int nEtaBins8alt = 8;
  const double etaBinLimits8alt[nEtaBins8alt + 1 ] =
    {0, 0.4, 0.8, 1.0, 1.479, 1.8, 2.0, 2.2, 2.40001 };
  const int nEtaBins8altNegs = 16;
  const double etaBinLimits8altNegs[nEtaBins8altNegs + 1 ] =
    {-2.40001, -2.2, -2.0, -1.8, -1.479, -1.0, -0.8, -0.4, 0, 0.4, 0.8, 1.0, 1.479, 1.8, 2.0, 2.2, 2.4 };
  const int nEtaBins14 = 14;
  const double etaBinLimits14[nEtaBins14 + 1 ] = 
    {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4442, 1.566, 1.8, 2.0, 2.2, 2.3, 2.4, 2.500001 };

  const int nEtaBinsMax= nEtaBins5altNegs;

  inline
  int mergeEtBins(int binning) {
    return ((binning==ETABINS5) || (binning==ETABINS5_max25)
	    || (binning==ETABINS5egamma)
	    ) ? 1:0;
  }


  inline 
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
    case ETABINS5: n = nEtaBins5; break;
    case ETABINS5_max25: n = nEtaBins5_max25; break;
    case ETABINS4test: n = nEtaBins4test; break;
    case ETABINS4testNegs: n = nEtaBins4testNegs; break;
    case ETABINS4alt: n = nEtaBins4alt; break;
    case ETABINS4altNegs: n = nEtaBins4altNegs; break;
    case ETABINS5alt: n = nEtaBins5alt; break;
    case ETABINS5altNegs: n = nEtaBins5altNegs; break;
    case ETABINS8alt: n = nEtaBins8alt; break;
    case ETABINS8altNegs: n = nEtaBins8altNegs; break;
    case ETABINS14: n = nEtaBins14; break;
    default:
      printf("ERROR: unknown binning requested\n");
      assert(0);
      n=0;
    }
    return n;
  }

  inline 
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
    case ETABINS5: limits = etaBinLimits5; break;
    case ETABINS5_max25: limits = etaBinLimits5_max25; break;
    case ETABINS4test: limits = etaBinLimits4test; break;
    case ETABINS4testNegs: limits = etaBinLimits4testNegs; break;
    case ETABINS4alt: limits = etaBinLimits4alt; break;
    case ETABINS4altNegs: limits = etaBinLimits4altNegs; break;
    case ETABINS5alt: limits = etaBinLimits5alt; break;
    case ETABINS5altNegs: limits = etaBinLimits5altNegs; break;
    case ETABINS8alt: limits = etaBinLimits8alt; break;
    case ETABINS8altNegs: limits = etaBinLimits8altNegs; break;
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

  inline
  int signedEtaBinning(int binning) {
    int yes=0;
    switch(binning) {
    case ETABINS1: 
    case ETABINS2: 
    case ETABINS3: 
    case ETABINS5egamma: // identical, although eff is with |eta|<2.5
    case ETABINS5:
    case ETABINS5_max25:
    case ETABINS4test:
    case ETABINS4alt:
    case ETABINS5alt:
    case ETABINS8alt:
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


  // Primary vertices 

  //const int nPVBinCount=7;
  //const double nPVLimits[nPVBinCount+1] = { 0., 5., 10., 15., 20., 25., 30., 100. };
  const int nPVBinCount=11;
  const double nPVLimits[nPVBinCount+1] = { 0.5, 2.5, 4.5, 6.5, 8.5, 10.5, 12.5, 14.5, 16.5, 20.5, 24.5, 40.5 };
  //  1., 3., 5., 7., 9., 11., 13., 15., 17., 21., 25., 40. };

  inline int findPUBin(int nPV) { return _findMassBin(double(nPV),nPVBinCount,nPVLimits); }
  
  // 
  // Cross section types
  //
  typedef enum { _cs_None=0, _cs_preFsr, _cs_preFsrNorm, 
		 _cs_preFsrDet, _cs_preFsrDetNorm,
		 _cs_preFsrDetErr, _cs_preFsrDetNormErr,
		 _cs_preFsrDetSystErr, _cs_preFsrDetNormSystErr,
		 _cs_postFsr, _cs_postFsrNorm, 
		 _cs_postFsrDet, _cs_postFsrDetNorm } TCrossSectionKind_t;

  // 
  // Triggers vs run numbers
  //
  //  enum { UNDEF, REL38X, REL39X};
  typedef enum { DATA, MC} TDataKind_t;


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

  inline
  bool isEcalGap(double eta) {
    return ( fabs(eta) > kECAL_GAP_LOW && fabs(eta) < kECAL_GAP_HIGH );
  }

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


  inline 
  bool goodEtaPair(double eta1, double eta2) {
    return (goodEta(eta1) && goodEta(eta2));
  }

  inline int goodEtEtaPair(double et1, double eta1, double et2, double eta2) {
    return (goodEtPair(et1,et2) && goodEtaPair(eta1,eta2));
  }

  //
  // Several functions to calculate efficiency or ratio
  //
  /*
  enum {EFF_POISSON, EFF_BINOMIAL, EFF_CLOPPER_PEARSON};

  inline 
  void calcEfficiency(double nPass, double nTotal, int method, 
		      double &eff, double &effErrLow, double &effErrHigh)
  {

    eff        = 0;
    effErrLow  = 0;
    effErrHigh = 0;

    if(nTotal == 0) return;
      
    // The efficiency is always just the ratio, it is in the error
    // calculation that the methods differ
    eff = nPass/nTotal;

    if( nPass > nTotal && method != EFF_POISSON){
      printf("calcEfficiency::WARNING requested efficiency calculation for nPass > nTotal,\n");
      printf("     switching to EFF_SIMPLE_DIVIDE method\n");
      method = EFF_POISSON;
    }

    if(method == EFF_POISSON){
      // Numerator and denominator are assumed uncorrelated
      if( nPass != 0 )
	effErrLow =  eff*sqrt( 1/nPass + 1/nTotal );
      else
	effErrLow = 0;
      effErrHigh = effErrLow;
    } else if (method == EFF_BINOMIAL){
      // Numerator is a subset of denominator
      effErrLow = sqrt(eff*(1-eff)/nTotal);
      effErrHigh = effErrLow;
    } else if (method == EFF_CLOPPER_PEARSON) {
      // Use Clopper-Pearson method implemented in ROOT to get +-1 sigma
      // asymmertic errors.
      // a) Clopper-Pearson method requires integer pass/total.
      // b) C++ does not have a built-in rounding function
      int nPassInt = int((nPass-floor(nPass)<0.5)?floor(nPass):ceil(nPass)+1e-3);
      int nTotalInt = int((nTotal-floor(nTotal)<0.5)?floor(nTotal):ceil(nTotal)+1e-3);
      effErrLow  = eff - TEfficiency::ClopperPearson(nTotalInt, nPassInt, 0.68269, kFALSE);
      effErrHigh = TEfficiency::ClopperPearson(nTotalInt, nPassInt, 0.68269, kTRUE) - eff;
    } else{
      printf("CalcEfficiency::ERROR: unknown method requested\n");
    }

    return;
  }
  */
};

// --------------------------------------------------------

namespace DYTools {

  inline int checkTotalLumi(double lumi) {
    if (fabs(lumi-DYTools::lumiAtECMS)>lumiAtECMS_accuracy) {
      std::cout << "checkTotalLumi detected mismatch: lumi=" << lumi 
		<< ", DYTools::lumiAtECMS=" << DYTools::lumiAtECMS
		<< " +/- " << DYTools::lumiAtECMS_accuracy << std::endl;
      return 0;
    }
    return 1;
  }
 
};

// ------------------------------------------------------------------

const int retCodeOk=7;
const int retCodeError=13;
const int retCodeStop=11;

// ------------------------------------------------------------------



#endif
