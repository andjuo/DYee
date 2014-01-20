#ifndef FEWZ___HH
#define FEWZ___HH

#include "DYTools.hh"

// --------------- 7 TeV analysis -------------------
#ifdef DYee7TeV
#include "Inc_7TeV/FEWZ.hh"
#endif
// --------------- 8 TeV analysis -------------------
#ifdef DYee8TeV
#include "Inc_8TeV/FEWZ.hh"
#endif
// --------------- 8 TeV analysis (regressed) -------------------
#ifdef DYee8TeV_reg
#include "Inc_8TeV_reg/FEWZ.hh"
#endif
// --------------------------------------------------


#endif
