
#warning included rangeDef_withFullOverflow_inc.h

const TString analysisTag_binning="_mbFullOverflow";
const DYTools::TMassBinning_t massBinningSet= _MassBins_withFullOverflow;

// Constants that define binning in mass and rapidity
const int _nMassBins2D = 8;
const double _massBinLimits2D[_nMassBins2D+1] = 
  {
    0, // underflow
    20, 30, 45, 60, 120, 200, 1500,
    3000 // overflow
  }; 
// Rapidity binning is different for different mass bins
// yRangeMax has no meaning. There is no limit in findYAbsIndex
const double _yRangeEdge= _yRangeEdge_base;
const double _yRangeMax = _yRangeEdge + 0.1;
const int _nBinsYLowMass1  = _nBinsYLowMass + 1;
const int _nBinsYHighMass1 = _nBinsYHighMass + 1;
const int _nYBinsMax2D = _nBinsYLowMass1;
const int _nYBins2D[_nMassBins2D] = 
  { 
    _nBinsYLowMass1, _nBinsYLowMass1,
    _nBinsYLowMass1, _nBinsYLowMass1, 
    _nBinsYLowMass1, _nBinsYLowMass1,
    _nBinsYHighMass1,
    _nBinsYHighMass1
  }; // overflow is neglected

const int _nMassBins1D = _nMassBins2D;
const double *_massBinLimits1D= _massBinLimits2D;
const int _nYBinsMax1D= 1;

const int _nYBins1D[_nMassBins1D] = { 1, 1, 1, 1, 1, 1, 1, 1 };

