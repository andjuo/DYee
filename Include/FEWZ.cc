#ifndef FEWZ_CC
#define FEWZ_CC

#include "DYTools.hh"
#include "FEWZ.hh"

#ifdef DYee7TeV
#include "Inc_7TeV/FEWZ.cc"
#endif
#ifdef DYee8TeV
#include "Inc_8TeV/FEWZ.cc"
#endif
#ifdef DYee8TeV_reg
#include "Inc_8TeV_reg/FEWZ.cc"
#endif


#endif
