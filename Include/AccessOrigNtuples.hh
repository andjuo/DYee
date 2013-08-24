#ifndef AccessOrigNtuples_HH
#define AccessOrigNtuples_HH


#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TClonesArray.h>
#include <iostream>
#include <ostream>
#include <string>
#include <math.h>

#include "../Include/DYTools.hh"

/*
#include "../Include/ElectronEnergyScale.hh"
#include "../Include/TDielectron.hh"
#include "../Include/TVertex.hh"
#include "../Include/EleIDCuts.hh"
#include "../Include/eventCounter.h"
#include "../Include/BaseClass.hh"
*/

#include "../Include/TEventInfo.hh"
#include "../Include/TGenInfo.hh"
#include "../Include/TDielectron.hh"
#include "../Include/TVertex.hh"
#include "../Include/JsonParser.hh"

#include "../Include/BaseClass.hh"
#include "../Include/TriggerSelection.hh"


// -------------------------------------------------

class AccessOrigNtuples_t : BaseClass_t {
protected:
  int hasJson;
  JsonParser json;
  TTree *tree;
  mithep::TEventInfo *info;
  mithep::TGenInfo *gen;
  TClonesArray *dielectronArr;
  TClonesArray *pvArr;

  TBranch *infoBr;
  TBranch *genBr;
  TBranch *dielectronBr;
  TBranch *pvBr;
  int genBrIsActive;
public:
  AccessOrigNtuples_t() : BaseClass_t("AccessOrigNtuples_t"),
    hasJson(0),
    json(),
    tree(NULL),
    info(new mithep::TEventInfo()),
    gen(new mithep::TGenInfo()),
    dielectronArr(new TClonesArray("mithep::TDielectron")),
    pvArr(new TClonesArray("mithep::TVertex")),
    infoBr(NULL),
    genBr(NULL),
    dielectronBr(NULL),
    pvBr(NULL),
    genBrIsActive(0)
  {}

  /* this dealoccation causes ROOT crash
  ~AccessOrigNtuples_t() {
    if (tree) delete tree;
    if (info) delete info;
    if (gen) delete gen;
    if (dielectronArr) delete dielectronArr;
    if (pvArr) delete pvArr;
    if (infoBr) delete infoBr;
    if (genBr) delete genBr;
    if (dielectronBr) delete dielectronBr;
    if (pvBr) delete pvBr;
  }
  */

  const mithep::TGenInfo* genPtr() const { return (genBrIsActive) ? gen : NULL; }
  //const mithep::TEventInfo* infoPtr() const { return info; }
  const mithep::TEventInfo* evtInfoPtr() const { return info; }
  int dielectronCount() const { return dielectronArr->GetEntriesFast(); }
  const mithep::TDielectron* dielectronPtr(int i) const { return (const mithep::TDielectron*)((*dielectronArr)[i]); }
  mithep::TDielectron* editDielectronPtr(int i) { return (mithep::TDielectron*)((*dielectronArr)[i]); }

  int getNPV(int isData) const {
    int npv=0;
    if (isData) npv=pvArr->GetEntriesFast();
    else npv=info->nPU;
    return npv;
  }

  int setTree(TFile &infile, const TString &treeName, int setupGenBr) {
    tree = (TTree*)infile.Get(treeName);
    infoBr=NULL;
    genBr=NULL;
    genBrIsActive=0;
    dielectronBr=NULL;
    pvBr=NULL;
    dielectronArr->Clear();
    pvArr->Clear();
    int ok=(tree) ? 1:0;
    if (ok) {
      tree->SetBranchAddress("Info",&info);
      infoBr=tree->GetBranch("Info");
      if (!infoBr) ok=0;
    }
    if (ok) {
      tree->SetBranchAddress("Dielectron",&dielectronArr);
      dielectronBr=tree->GetBranch("Dielectron");
      if (!dielectronBr) ok=0;
    }
    if (ok) {
      tree->SetBranchAddress("PV",&pvArr);
      pvBr=tree->GetBranch("PV");
      if (!pvBr) ok=0;
    }
    if (ok && setupGenBr) {
      tree->SetBranchAddress("Gen",&gen);
      genBr=tree->GetBranch("Gen");
      if (!genBr) ok=0; else genBrIsActive=1;
    }
    return ok;
  }

  ULong_t getEntries() { 
    if (!tree) { reportError("getEntries: tree is null\n"); return 0; }
    return tree->GetEntries();
  }

  ULong_t getEntriesFast() {
    if (!tree) { reportError("getEntriesFast: tree is null\n"); return 0; }
    return tree->GetEntriesFast();
  }

  template<class UInt_type>
  void GetInfoEntry(UInt_type ientry) { 
    dielectronArr->Clear();
    pvArr->Clear();
    infoBr->GetEntry(ientry); 
  }

  template<class UInt_type>
  void GetGen(UInt_type ientry) { 
    genBr->GetEntry(ientry);
  }

  template<class UInt_type>
  void GetDielectrons(UInt_type ientry) { 
    dielectronArr->Clear();
    dielectronBr->GetEntry(ientry); 
  }

  template<class UInt_type>
  void GetPVs(UInt_type ientry) { 
    pvArr->Clear();
    pvBr->GetEntry(ientry); 
  }

  int genLeptonsAreElectrons() const {
    return ((abs(gen->lid_1)==11) && (abs(gen->lid_2)==11)) ? 1:0;
  }

  int prepareJson(const TString &jsonFname) {
    hasJson=0;
    if (jsonFname.Length() &&
	(jsonFname.CompareTo("NONE")!=0)) {
      if (!json.Initialize(jsonFname)) return 0;
      hasJson=1;
    }
    return 1;
  }

  int eventInJson() const {
    if (!hasJson) return 1;
    return json.HasRunLumi(info->runNum,info->lumiSec);
  }
  
  //int eventTriggerOk(const TriggerSelection_t &trigger) const {
  //  return (trigger.matchEventTriggerBit(info->triggerBits,info->runNum)) ? 1:0;
  //}

};


// -------------------------------------------------
// -------------------------------------------------


#endif
