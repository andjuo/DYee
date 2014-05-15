// Define finer mass range

#warning included rangeDef_finerMassRange_inc.h

namespace DYTools {

const TString _defined_analysisTag_binning="_mbFinerM";
const DYTools::TMassBinning_t _defined_massBinningSet=_MassBins_finerMassRange;

// Constants that define binning in mass and rapidity
const int _nMassBins2D = 14;
const double _massBinLimits2D[_nMassBins2D+1] = 
  { 0,
    20, 25, 30, 35, 45, 50, 60, 90, 120, 160, 200, 1000, 1500,
    2000
  }; // overflow is very unlikely, do not account for it
const int _nYBinsMax2D=24;
const int _nYBins2D[_nMassBins2D] =
  {
    24,
    24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 12, 12,
    12
  }; // overflow is neglected


const int _nMassBins1D = 81;
const double _massBinLimits1D[_nMassBins1D+1] =
  {
   15, 17.5, 20, 22.5, 25, 27.5, 30, 32.5, 35, 37.5,
   40, 42.5, 45, 47.5, 50, 52.5, 55, 57.5, 60, 62.0,
   64, 66.0, 68, 70.0, 72, 74.0, 76, 78.5, 81, 83.5,
    86,  88.5,  91,  93.5,  96,  98.5, 101,  103.5, 106, 108.0,
   110, 112.5, 115, 117.5, 120, 123.0, 126, 129.5, 133, 137.0,
   141, 145.5, 150, 155.0, 160, 165.5, 171, 178.0, 185, 192.5,
   200, 210.0, 220, 231.5, 243, 258.0, 273, 296.5, 320, 350.0,
   380, 410.0, 440, 475.0, 510, 555.0, 600, 800.0, 1000, 1250.,
   1500, 2000 }; // 82 bin  // no data between 1750 and 2000
const int _nYBinsMax1D=1; // the largest division into Y bins
const int _nYBins1D[_nMassBins1D] = {
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
  1
};

 const double _extend_Y = 6.5;
 const double _yRangeEdge=electronEtaMax;
 const double _yRangeMax=electronEtaMax;

};
