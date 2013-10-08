#ifndef EffCalc_HH
#define EffCalc_HH

#include "../Include/DYTools.hh"
#include <TEfficiency.h>


namespace DYTools {
  //
  // Several functions to calculate efficiency or ratio
  //

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
};


#endif
