// Define ranges without underflow in 2D

#warning included rangeDef_noUnderflow_inc.h


const TString analysisTag_binning="_mbNoUnderflow";
const DYTools::TMassBinning_t massBinningSet= _MassBins_noUnderflow;

// Constants that define binning in mass and rapidity
const int _nMassBins2D = 6;
const double _massBinLimits2D[_nMassBins2D+1] = 
  {
    20, 30, 45, 60, 120, 200, 1500
  }; // overflow is very unlikely, do not account for it
const int _nYBinsMax2D=_nBinsYLowMass;
const int _nYBins2D[_nMassBins2D] = 
  { 
    _nBinsYLowMass, 
    _nBinsYLowMass, _nBinsYLowMass, 
    _nBinsYLowMass, _nBinsYLowMass,
    _nBinsYHighMass
  }; // overflow is neglected


const int _nMassBins1D = _nMassBins2D;
const double *_massBinLimits1D= _massBinLimits2D;
const int _nYBinsMax1D= 1;

const int _nYBins1D[_nMassBins1D] = {
  1, 1, 1, 1, 1, 1
};

const double _yRangeEdge=_yRangeEdge_base;
const double _yRangeMax=_yRangeMax_base; 
