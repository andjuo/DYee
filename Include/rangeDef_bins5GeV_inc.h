// Based on 100GeV
// May 20, 2014

#warning included rangeDef_bins5GeV_inc.h

namespace DYTools {

const TString _defined_analysisTag_binning="_mb5GeV";
const DYTools::TMassBinning_t _defined_massBinningSet= _MassBins_bins5GeV;

// Constants that define binning in mass and rapidity
const int _nMassBins2D = 26;
const double _massBinLimits2D[_nMassBins2D+1] =
  {
    0, // underflow
    10, 15, 20, 25, 35, 40, 45, 50, 55, 60, 65, 70,
    120, 200,
    300, 400, 500, 600, 700, 800, 900, 1000,
    1200, 1500, 2000,
    3000 // overflow
  };
// Rapidity binning is different for different mass bins
// yRangeMax has no meaning. There is no limit in findYAbsIndex
const double _extend_Y = 6.5;
const double _yRangeEdge= electronEtaMax;
const double _yRangeMax = electronEtaMax;
const int _nYBinsMax2D = _nBinsYLowMass;
const int _nYBins2D[_nMassBins2D] =
  {
    _nBinsYLowMass, _nBinsYLowMass, _nBinsYLowMass, _nBinsYLowMass,
    _nBinsYLowMass, _nBinsYLowMass, _nBinsYLowMass, _nBinsYLowMass,
    _nBinsYLowMass, _nBinsYLowMass, _nBinsYLowMass, _nBinsYLowMass,
    _nBinsYLowMass, _nBinsYLowMass, _nBinsYLowMass,
    _nBinsYLowMass, _nBinsYLowMass
    //_nBinsYHighMass,
    //_nBinsYHighMass
  }; // overflow is neglected

const int _nMassBins1D = _nMassBins2D;
const double *_massBinLimits1D= _massBinLimits2D;
const int _nYBinsMax1D= 1;

const int _nYBins1D[_nMassBins1D] = {
  1, 1, 1, 1, 1,
  1, 1, 1, 1, 1,
  1, 1, 1, 1, 1,
  1, 1, 1, 1, 1,
  1, 1, 1, 1, 1,
  1
};

}; // namespace DYTools
