#warning  included rangeDef_ZpeakRegion_inc.h

namespace DYTools {

  // Z-peak region 1GeV bins
const TString _defined_analysisTag_binning="_mbZpeak";
const DYTools::TMassBinning_t _defined_massBinningSet= _MassBins_Zpeak;

const int _nMassBins1D = 60;
const double _massBinLimits1D[_nMassBins1D+1] = 
  { 60, 61, 62, 63, 64, 65, 66, 67, 68, 69,
    70, 71, 72, 73, 74, 75, 76, 77, 78, 79,
    80, 81, 82, 83, 84, 85, 86, 87, 88, 89,
    90, 91, 92, 93, 94, 95, 96, 97, 98, 99,
    100, 101, 102, 103, 104, 105, 106, 107, 108, 109,
    110, 111, 112, 113, 114, 115, 116, 117, 118, 119,
    120 };
const int _nYBinsMax1D=1; // the largest division into Y bins
const int _nYBins1D[_nMassBins1D] = { 
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1
};

const int _nMassBins2D = _nMassBins1D;
const double *_massBinLimits2D= _massBinLimits1D;
const int _nYBinsMax2D= _nYBinsMax1D;
const int *_nYBins2D= _nYBins1D;

const double _extend_Y = 6.5;
const double _yRangeEdge= electronEtaMax;
const double _yRangeMax = electronEtaMax;

}; // namespace DYTools
