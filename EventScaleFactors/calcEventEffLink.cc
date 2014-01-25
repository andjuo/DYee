#include "calcEventEffLink.h"

#include "calcEventEff.C"

int fillEfficiencyConstants( const InputFileMgr_t &inpMgr ) {
  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;
  return fillEfficiencyConstants(inpMgr,systMode);
}
