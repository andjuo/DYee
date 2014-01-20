#ifndef EWKAnaDefs_HH
#define EWKAnaDefs_HH

#include "DYTools.hh"

// --------------- 7 TeV analysis -------------------
#ifdef DYee7TeV
#include "Inc_7TeV/EWKAnaDefs.hh"
#endif
// --------------- 8 TeV analysis -------------------
#ifdef DYee8TeV
#include "Inc_8TeV/EWKAnaDefs.hh"
#endif
// --------------- 8 TeV analysis (regressed) -------------------
#ifdef DYee8TeV_reg
#include "Inc_8TeV_reg/EWKAnaDefs.hh"
#endif
// --------------------------------------------------


#endif
