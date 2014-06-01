#ifndef ElectronEnergyScale_HH
#define ElectronEnergyScale_HH

#include <TString.h>
#include <TF1.h>
#include <TRandom.h>
#include <TH1D.h>

//#define UseEEM
#define UseZeeData

#ifdef UseEEM
#include "../Include/EtaEtaMass.hh"
#endif

#ifdef UseZeeData
#include "../Include/ZeeData.hh"
#endif


class ElectronEnergyScale {

public:

  enum CalibrationSet {
    UNDEFINED=0,
    UNCORRECTED,
    Date20110901_EPS11_default,
    Date20120101_default,
    Date20120802_default,
    Date20121003FEWZ_default, // FEWZ correction and MC backgroud were considered in derivation
    Date20121025FEWZPU_default, // FEWZ correction and PU reweighting were considered in derivation
    Date20130529_2012_j22_adhoc,  // a guess to what the calib. should be for Jan22 8 TeV data
    Date20140220_2012_j22_peak_position_IK, // scaling derived so that the peak in data and mc matches, another try. Suggested by Ilya K.
    Date20140220_2012_j22_peak_position, // fewer eta ranges to closer match the EGamma binning
                                         // for reco with momentum regression. No MC smearing.
    Date20140220_2012_j22_peak_position_flat, // one coefficient
    CalSet_File_Gauss,
    CalSet_File_BreitWigner,
    CalSet_File_Voigt
  };

  // Constructor
  ElectronEnergyScale(CalibrationSet calibrationSet);
  ElectronEnergyScale(const TString &escaleTagName);
  // Destructor
  //~ElectronEnergyScale() { this->clear(); }
  // Clean-up
  void clear();

  // Initialization
  void init(CalibrationSet calSet, int debug=0);
  void init(const TString &stringWithEScaleTagName, int debug=0);
  bool initializeAllConstants(int debug=0);
  bool initializeExtraSmearingFunction(int normalize=0);
  bool isInitialized() const;
  bool randomizedStudy() const { return _randomizedStudy; }
  void randomizedStudy(bool doRandomizedStudy) { _randomizedStudy=doRandomizedStudy; }
  bool isScaleRandomized() const { return _energyScaleCorrectionRandomizationDone; }
  bool isSmearRandomized() const { return _smearingWidthRandomizationDone; }

  bool loadInputFile(const TString &fileName, int debug=0);
  bool assignConstants(const std::vector<std::string> &lines, int debug=0);
  static bool AssignConstants(const std::vector<std::string> &lines, int count, double *eta, double *scale, double *scaleErr, double *smear, double *smearErr, int debug=0);
  static bool AssignConstants2(const std::vector<std::string> &lines, int count, const char *parameterNameStart, double *par, double *parErr, int debug=0);

  // Access
  CalibrationSet getCalibrationSet() const { return _calibrationSet; }
  bool setCalibrationSet(CalibrationSet calSet); // in specific cases allow to change the calibration set
  int numberOfEtaBins() const { return _nEtaBins; }
  int numberOfEtaEtaBins() const { return _nEtaBins*(_nEtaBins+1)/2; }

  // Eta bin index in the same way as TH1D bin index
  int getEtaBinIdx(double eta) const {
    int idx=-1;
    for (int i=0; (idx<0) && (i<_nEtaBins); ++i) {
      if ((eta>_etaBinLimits[i]) && (eta<=_etaBinLimits[i+1])) idx=i+1;
    }
    return idx;
  }

  double getEtaLowEdge(int idx) const { return _etaBinLimits[idx]; }
  double getEtaHighEdge(int idx) const { return _etaBinLimits[idx+1]; }

  // triangular (i,j) index
  int getEtaEtaIdx(double eta1, double eta2) const {
    int idx1=getEtaBinIdx(eta1)-1;
    int idx2=getEtaBinIdx(eta2)-1;
    int idx=-1;
    if ((idx1>=0) && (idx2>=0)) {
      if (idx1>idx2) { idx=idx1; idx1=idx2; idx2=idx; }
      idx=(2*_nEtaBins-1-idx1)*idx1/2 + idx2;
    }
    return idx;
  }

  double getEtaBinCenter(int i) const { return 0.5*(_etaBinLimits[i-1]+_etaBinLimits[i]); }

  double getEnergyScaleCorrection(double eta) const;
  // old-style smear (event shift)
  double generateMCSmear(double eta1, double eta2) const;
  // updated smear (distribution) : one event
  bool addSmearedWeight(TH1D *hMassDestination, int eta1Bin, int eta2Bin, double mass, double weight) const {
    assert(this->isInitialized());
    return addSmearedWeightAny(hMassDestination,eta1Bin,eta2Bin,mass,weight,kFALSE);
  }
  // updated smear (distribution) : smear collection
  void smearDistribution(TH1D *destination, int eta1Bin, int eta2Bin, const TH1D *source) const {
    assert(this->isInitialized());
    smearDistributionAny(destination,eta1Bin,eta2Bin,source,kFALSE);
  }

  // Several functions for systematic studies. These
  // randomize constants for energy scale corrections
  // within their errors.
  int   randomizeEnergyScaleCorrections(int seed);
  double getEnergyScaleCorrectionRandomized(double eta) const;
  // for systematics studies with 2012 data, MC has to be scaled by the
  // inverse factors as compared to data
  int   invertRandomizedEnergyScaleCorrections();
  int   mirrorRandomizedEnergyScaleCorrections();

  void   randomizeSmearingWidth(int seed);

  // old-style smear (event shift)
  double generateMCSmearRandomized(double eta1, double eta2) const;
  // updated smear (distribution) : one event
  bool addSmearedWeightRandomized(TH1D *hMassDestination, int eta1Bin, int eta2Bin, double mass, double weight) const {
    assert(this->isInitialized()); assert(this->isSmearRandomized());
    return addSmearedWeightAny(hMassDestination,eta1Bin,eta2Bin,mass,weight,kTRUE);
  }
  // updated smear (distribution) : smear collection
  void smearDistributionRandomized(TH1D *destination, int eta1Bin, int eta2Bin, const TH1D *source) const {
    assert(this->isInitialized()); assert(this->isSmearRandomized());
    smearDistributionAny(destination,eta1Bin,eta2Bin,source,kTRUE);
  }

  // Several functions for working with smearing single MC electrons
  double generateMCSmearSingleEle(double eta) const;
  double generateMCSmearSingleEleRandomized(double eta) const;
  double generateMCSmearAnySingleEle(double eta, bool randomize) const;

  // Encapsulated electron energy scaling

  template<class tData_t>
  int applyEnergyScale(tData_t *data, int isData, int randomized=0) const {
    int res=1;
    if (!isData) {
      double smearCorr= this->generateMCSmearAny(data->scEta_1,data->scEta_2,randomized);
      data->mass += smearCorr;
    }
    else {
      res= this->applyDataEnergyScale(data,randomized);
    }
    return res;
  }

  template<class tData_t>
  int applyDataEnergyScale(tData_t *data, int randomized) const {
    double corr1= getEnergyScaleCorrectionAny(data->scEta_1,randomized);
    double corr2= getEnergyScaleCorrectionAny(data->scEta_2,randomized);
    TLorentzVector ele1,ele2,ee;
    ele1.SetPtEtaPhiM(data->pt_1,data->eta_1,data->phi_1,0.000511);
    ele2.SetPtEtaPhiM(data->pt_2,data->eta_2,data->phi_2,0.000511);
    ele1*=corr1;
    ele2*=corr2;
    data->scEt_1 *= (ele1.Pt()/data->pt_1);
    data->pt_1 = ele1.Pt();
    data->scEta_1 *= (ele1.Eta()/data->eta_1);
    data->eta_1= ele1.Eta();
    data->scPhi_1 *= (ele1.Phi()/data->scPhi_1);
    data->phi_1= ele1.Phi();
    data->scEt_2 *= (ele2.Pt()/data->pt_2);
    data->pt_2 = ele2.Pt();
    data->scEta_2 *= (ele2.Eta()/data->eta_2);
    data->eta_2= ele2.Eta();
    data->scPhi_2 *= (ele2.Phi()/data->scPhi_2);
    data->phi_2= ele2.Phi();
    ee = ele1 + ele2;
    data->mass = ee.M();
    data->pt   = ee.Pt();
    data->y    = ee.Rapidity();
    data->phi  = ee.Phi();
    return 1;
  }

#ifdef UseZeeData
  int applyEnergyScale_instantiate(ZeeData_t *data) const {
    return this->applyEnergyScale(data,1,1);
  }
#endif

  // print
  void print() const;
  void printAsTexTable(const TString &fname) const;

  // Useful functions
  static CalibrationSet DetermineCalibrationSet(const TString &escaleTagName, TString *inputFileName=NULL);  // TString -> CalibrationSet
  static TString ExtractFileName(const TString &escaleTagName);
  static TString CalibrationSetName(CalibrationSet escaleTag, const TString *fileName);   // CalibrationSet -> TString
  static TString CalibrationSetFunctionName(CalibrationSet escaleTag); // name of the calibrating function

  TString calibrationSetName() const { return ElectronEnergyScale::CalibrationSetName(this->_calibrationSet, &this->_inpFileName); }
  TString calibrationSetFunctionName() const { return ElectronEnergyScale::CalibrationSetFunctionName(this->_calibrationSet); }
  TString calibrationSetShortName() const;

  TH1D* createScaleHisto(const TString &namebase) const;
  TH1D* createSmearHisto(const TString &namebase, int parameterNo) const;

#ifdef UseEEM
  int loadEEMFile(const TString &eemFileName, std::vector<std::vector<double>*> &eemData) const;
  int loadEEMFile(const TString &eemFileName, std::vector<std::vector<double>*> &eemData, double massMin, double massMax) const;
  int loadEEMFile(const TString &eemFileName, std::vector<std::vector<EtaEtaMassData_t>*> &eemData, double massMin=-1e9, double massMax=1e9) const;
#endif

#ifdef ZEE_DATA_HH
  int ProcessZeeDataFile(const char *mc_file_name, const char *data_file_name, std::vector<std::vector<double>*> &mcData, std::vector<std::vector<double>*> &expData, unsigned int debugLoad=0);
  int ProcessZeeDataFile(const char *mc_file_name, const char *data_file_name, std::vector<std::vector<ZeeData_t*>*> &mcData, std::vector<std::vector<ZeeData_t*>*> &expData, double mass_min, double mass_max, int load_only_idx=-1, unsigned int debugLoad=0);
  int ProcessZeeDataFiles(const std::vector<std::string> &mc_file_names, const char *data_file_name, std::vector<std::vector<double>*> &mcData, std::vector<std::vector<double>*> &expData, unsigned int debugLoad=0);
  int ProcessZeeDataFileApproximateMCWeight(const char *mc_file_name, const char *data_file_name, std::vector<std::vector<double>*> &mcData, std::vector<std::vector<double>*> &expData, double unitWeight=0., int include_fewz_weigth=1, int puReweight=1, unsigned int debugLoad=0);
  int ProcessZeeDataFileApproximateMCWeight(const std::vector<std::string> &mc_file_names, const char *data_file_name, std::vector<std::vector<double>*> &mcData, std::vector<std::vector<double>*> &expData, double unitWeight=0., int include_fewz_weight=1, int puReweight=1, unsigned int debugLoad=0);
    //int ProcessEEMFileApproximateMCWeights(const char *mc_file_name, const char *data_file_name, std::vector<std::vector<double>*> &mcData, std::vector<std::vector<double>*> &expData);
  //int ProcessEEMFile(const char *mc_file_name, const char *data_file_name, std::vector<std::vector<double>*> &mcData, std::vector<std::vector<double>*> &mcWeights, std::vector<std::vector<double>*> &expData);
#endif


  //protected: 
  // Internal functions, not for general use
public: // made public to be able to check whether the randomized value is different from non-randomized

  double getEnergyScaleCorrectionAny(double eta, bool randomize) const;
  double getEnergySmearCorrection(double eta, int the_set=0) const;

  // old-style smear (event shift)
  double generateMCSmearAny(double eta1, double eta2, bool randomize) const;
  // updated smear (distribution) : one event
  bool addSmearedWeightAny(TH1D *hMassDestination, int eta1Bin, int eta2Bin, double mass, double weight, bool randomize) const;
  // updated smear (distribution) : smear collection
  void smearDistributionAny(TH1D *destination, int eta1Bin, int eta2Bin, const TH1D *source, bool randomize) const;

  double* editDataConst() {
    if (!_randomizedStudy) {
      for (int i=0; i<_nEtaBins; i++) _dataConstRandomized[i]=_dataConst[i];
      _randomizedStudy=true;
    }
    return _dataConstRandomized; 
  }

protected:
  TH1D* createParamHisto(const TString &namebase, const TString &nameTag, const double *params, const double *paramErrs) const;

private:

  CalibrationSet     _calibrationSet;
  TString            _inpFileName;
  bool               _isInitialized;
  bool               _randomizedStudy;

  // The data structure assumes the following:
  //  - the binning is in eta only
  //  - data constants contain a single scale factor 
  //  - MC constants contain several variable describing the shape
  //         of the smearing function
  //  - each constant, data or MC, comes with an error

  // Data constants
  int                    _nEtaBins;
  double *               _etaBinLimits;
  double *               _dataConst;
  double *               _dataConstErr;
  // randomized for systematics studies
  double *               _dataConstRandomized;
  bool                   _energyScaleCorrectionRandomizationDone;
  int                    _dataSeed;
  
  // MC constants. The number of constants and the names depend
  // on particular calibration set. For now, maximum possible is four.
  int                    _nMCConstants;
  TString                _mcConst1Name;
  TString                _mcConst2Name;
  TString                _mcConst3Name;
  TString                _mcConst4Name;
  double *               _mcConst1;
  double *               _mcConst2;
  double *               _mcConst3;
  double *               _mcConst4;
  double *               _mcConst1Err;
  double *               _mcConst2Err;
  double *               _mcConst3Err;
  double *               _mcConst4Err;
  bool                   _smearingWidthRandomizationDone;
  int                    _mcSeed;

protected:
  // Functions to be used for extra smearing
  static const int nMaxFunctions = 50;
  TF1 *smearingFunctionGrid[nMaxFunctions][nMaxFunctions];
  TF1 *smearingFunctionGridRandomized[nMaxFunctions][nMaxFunctions];

};


#endif
