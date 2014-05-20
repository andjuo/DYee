
#warning included rangeDef_withFullOverflow_inc.h

const TString _defined_analysisTag_binning="_mbFullOverflow";
const DYTools::TMassBinning_t _defined_massBinningSet= DYTools::_MassBins_withFullOverflow;

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
const double _extend_Y = 6.5;
const double _yRangeEdge= 2.4;
const double _yRangeMax = _yRangeEdge + 0.1;
const int _nBinsYLowMass1  = 25;
const int _nBinsYHighMass1 = 12 + 1;
const int _nYBinsMax2D = _nBinsYLowMass1;
const int _nYBins2D[_nMassBins2D] = 
  { 
    _nBinsYLowMass1, _nBinsYLowMass1,
    _nBinsYLowMass1, _nBinsYLowMass1, 
    _nBinsYLowMass1, _nBinsYLowMass1,
    _nBinsYHighMass1,
    _nBinsYHighMass1
  }; // overflow is neglected

// 2012 mass binning
const int _nMassBins1D = 43;
const double _massBinLimits1D[_nMassBins1D+1] =
  {10,
   15, 20, 25, 30, 35, 40, 45, 50, 55, 60,
   64, 68, 72, 76, 81, 86, 91, 96, 101,106,
   110, 115, 120, 126, 133, 141, 150, 160, 171, 185,
   200, 220, 243, 273, 320, 380, 440, 510, 600, 1000,
   1500, 2000, 3000}; // 43 bin
const int _nYBinsMax1D=1; // the largest division into Y bins
const int _nYBins1D[_nMassBins1D] = {
  1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1
};

