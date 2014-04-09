#warning included rangeDef_massBinTest4_inc.h

namespace DYTools {

const TString _defined_analysisTag_binning="_mbTest4";
const DYTools::TMassBinning_t _defined_massBinningSet= _MassBins_test4;

const double _extend_Y = 6.5;
const double _yRangeEdge= electronEtaMax;
const double _yRangeMax = electronEtaMax;

const int _nMassBins1D = 4;
const double _massBinLimits1D[_nMassBins1D+1] =
  {15,45,60,120,1500}; // 4 bins with Z-peak region singled-out
const int _nYBinsMax1D=1; // the largest division into Y bins
const int _nYBins1D[_nMassBins1D] = { 1, 1, 1, 1 };

const int _nMassBins2D = _nMassBins1D;
const double *_massBinLimits2D= _massBinLimits1D;
const int _nYBinsMax2D=2;
const int _nYBins2D[_nMassBins1D] = { 2, 2, 2, 1 };

}; // namespace DYTools
