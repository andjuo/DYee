#ifndef tnpSelectEvents_HH
#define tnpSelectEvents_HH

#define tnpSelectEventsIsObject
#define esfSelectEventsIsObject

#define tnpStoreTag

#if defined(tnpSelectEventsIsObject) || defined(esfSelectEventsIsObject)
#include <TObject.h>
#endif


#include <TTree.h>
#include <TBranch.h>

#include "../Include/PUReweight.hh"
#include "../Include/DYTools.hh"
#include "../Include/AccessOrigNtuples.hh"



// ------------------------------------------------------

#ifdef tnpSelectEventsIsObject
class tnpSelectEvent_t : public TObject 
#else
struct tnpSelectEvent_t
#endif
{
public:
  typedef enum { _skipWeight, _dontSkipWeight } TCreateBranchesOption_t;
public:

  Double_t mass, yZ; 
  Double_t et,eta; // probe
  UInt_t nGoodPV;
  Double_t eventWeight,weight;
#ifdef tnpStoreTag
  Double_t tagEt,tagEta;
#endif

#ifdef tnpStoreTag
  void assignTag(Double_t tag_et, Double_t tag_eta) {
    tagEt=tag_et; tagEta=tag_eta;
  }
#endif

  void assign(Double_t mass_1, Double_t yZee1, Double_t et_1, Double_t eta_1, 
	      UInt_t nGoodPV_1, Double_t event_weight1, Double_t weight1) {
    mass=mass_1; yZ=yZee1; et=et_1; eta=eta_1; 
    nGoodPV=nGoodPV_1;
    eventWeight=event_weight1; weight=weight1;
  }

  void createBranches(TTree *tree, TCreateBranchesOption_t opt) {
    tree->Branch("mass",&this->mass,"mass/D");
    tree->Branch("yZee",&this->yZ,"yZee/D");
    tree->Branch("et",&this->et,"et/D");
    tree->Branch("eta",&this->eta,"eta/D");
    tree->Branch("nGoodPV",&this->nGoodPV,"nGoodPV/i");
    tree->Branch("eventWeight",&this->eventWeight,"eventWeight/D");
    if (opt!=_skipWeight) tree->Branch("weight",&this->weight,"weight/D");
#ifdef tnpStoreTag
    tree->Branch("tagEt",&this->tagEt,"tagEt/D");
    tree->Branch("tagEta",&this->tagEta,"tagEta/D");
#endif
  }

  void setBranchAddress(TTree *tree) {
    tree->SetBranchAddress("mass",&this->mass);
    tree->SetBranchAddress("yZee",&this->yZ);
    tree->SetBranchAddress("et",&this->et);
    tree->SetBranchAddress("eta",&this->eta);
    tree->SetBranchAddress("nGoodPV",&this->nGoodPV);
    tree->SetBranchAddress("eventWeight",&this->eventWeight);
    tree->SetBranchAddress("weight",&this->weight);
#ifdef tnpStoreTag
    tree->SetBranchAddress("tagEt",&this->tagEt);
    tree->SetBranchAddress("tagEta",&this->tagEta);
#endif
  }

  bool insideMassWindow(double mass_low, double mass_high) const {
    return ((mass>=mass_low) && (mass<mass_high));
  }

  bool insideRapidityWindow(double rapidity_min, double rapidity_max) const {
    return ((yZ>=rapidity_min) && (yZ<rapidity_max));
  }

#ifdef tnpSelectEventsIsObject
  ClassDef(tnpSelectEvent_t,1)
#endif
};

// ------------------------------------------------------

#ifdef esfSelectEventsIsObject
class esfSelectEvent_t : public TObject 
#else
struct esfSelectEvent_t
#endif
{
public:

  Double_t genMass,genY,mass,y;
  Double_t et_1,eta_1,et_2,eta_2;
  Double_t weight;
  UInt_t nGoodPV;
  UInt_t largeFSR;

  esfSelectEvent_t() : TObject(),
		       genMass(0.), genY(0.), mass(0.), y(0.),
		       et_1(0.), eta_1(0.), et_2(0.), eta_2(0.),
		       weight(0.),
		       nGoodPV(0), largeFSR(0)
  {}

  /*
  void assign(double _genMass, double _genY,
	      double _mass, double _y1, 
	      double _et1, double _eta1,
	      double _et2, double _eta2, double _weight, UInt_t _nGoodPV) {
    genMass=_genMass; genY=_genY;
    mass=_mass; y=_y1;
    et_1=_et1; eta_1=_eta1;
    et_2=_et2; eta_2=_eta2;
    weight=_weight;
    nGoodPV=_nGoodPV;
  }
  */

  void assign(const mithep::TDielectron *dielectron, const mithep::TGenInfo *gen, UInt_t _nGoodPV, double _weight, UInt_t _largeFSR=0) {
    genMass= gen->mass;
    genY= gen->y;
    mass= dielectron->mass;
    y= dielectron->y;
    int leadingIs1st= ( dielectron->scEta_1 > dielectron->scEta_2 ) ? 1:0;
    if (leadingIs1st) {
      et_1= dielectron->scEt_1;
      eta_1= dielectron->scEta_1;
      et_2= dielectron->scEt_2;
      eta_2= dielectron->scEta_2;
    }
    else {
      et_1= dielectron->scEt_2;
      eta_1= dielectron->scEta_2;
      et_2= dielectron->scEt_1;
      eta_2= dielectron->scEta_1;
    }
    nGoodPV=_nGoodPV;
    weight=_weight;
    largeFSR=_largeFSR;
  }

  // 
  //void assign(const AccessOrigNtuples_t &accessInfo, int isData, int dielectronIdx, double _weight) {
  //  this->assign(accessInfo.dielectronPtr(dielectronIdx),
  //		 accessInfo.genPtr(),
  //		 accessInfo.getNPV(isData), 
  //		 _weight);
  //}


  void createBranches(TTree *tree) {
    tree->Branch("genMass",&this->genMass,"genMass/D");
    tree->Branch("genY",&this->genY,"genY/D");
    tree->Branch("mass",&this->mass,"mass/D");
    tree->Branch("y",&this->y,"y/D");
    tree->Branch("et_1",&this->et_1,"et_1/D");
    tree->Branch("eta_1",&this->eta_1,"eta_1/D");
    tree->Branch("et_2",&this->et_2,"et_2/D");
    tree->Branch("eta_2",&this->eta_2,"eta_2/D");
    tree->Branch("weight",&this->weight,"weight/D");
    tree->Branch("nGoodPV",&this->nGoodPV,"nGoodPV/i");
    tree->Branch("largeFSR",&this->largeFSR,"largeFSR/i");
  }

  void setBranchAddress(TTree *tree) {
    tree->SetBranchAddress("genMass",&this->genMass);
    tree->SetBranchAddress("genY",&this->genY);
    tree->SetBranchAddress("mass",&this->mass);
    tree->SetBranchAddress("y",&this->y);
    tree->SetBranchAddress("et_1",&this->et_1);
    tree->SetBranchAddress("eta_1",&this->eta_1);
    tree->SetBranchAddress("et_2",&this->et_2);
    tree->SetBranchAddress("eta_2",&this->eta_2);
    tree->SetBranchAddress("weight",&this->weight);
    tree->SetBranchAddress("nGoodPV",&this->nGoodPV);
    tree->SetBranchAddress("largeFSR",&this->largeFSR);
  }

  bool insideMassWindow(double mass_low, double mass_high) const {
    return ((mass>=mass_low) && (mass<=mass_high));
  }

#ifdef esfSelectEventsIsObject
  ClassDef(esfSelectEvent_t,2)
#endif
};

// ------------------------------------------------------
// ------------------------------------------------------

// THIS COMMENT NEEDS REVISION
// Load default PU distribution (data) from file <puReferenceFName> and use it to create
// weight branch in tag-and-probe selected events file <fname>.
// The selected events file is read to determine PU distribution (which is saved to file
// <savePUFName> in two histograms {savePUHistoNameBase}_pass and {savePUHistoNameBase}_fail
			  
int CreatePUWeightedBranch(DYTools::TSystematicsStudy_t systMode,
			   const TString &fName,
			   const TString &puTargetFName, const TString &puTargetDistrName,
			   const TString &puSourceFName, const TString &puSourceDistrName,			   
// 			   const TString &savePUFName, 
// 			   const TString &savePUHistoNameBase,
			   bool isMC 
			   ) {
  std::cout << "entered CreatePUWeightedBranch (" << fName << ")" << std::endl;

  // Set up pile-up reweighting
  PUReweight_t puReweight(systMode);
  // Only MC will be reweighted. 
  if( isMC ){
    if( DYTools::energy8TeV ){
      // For 8 TeV, standard Hildreth weights are used.
      // The method does not need to be set because above
      // it is set in the constructor by default.
    }else{
      // For 7 TeV, tt is reweighted to the 
      // tag and probe data sample. The method is not the plain
      // Hildreth method for 2011 because in 2011 tag and probe
      // triggers are prescaled, so a specially prepared 
      // "target" histogram is used.
      puReweight.setActiveMethod(PUReweight_t::_TwoHistos);
      puReweight.setSimpleWeights( puTargetFName, puTargetDistrName,
				   puSourceFName, puSourceDistrName);
    }
  }

  // COMMENTED OUT BELOW: THE PREVIOIS VERSION WHERE BOTH T&P
  // MC AND DATA ARE REWEIGHTED WITH RESPECT TO THE SIGNAL SAMPLE
  //
//   // We will reweight both TnP data and MC to the signal data
//   // using the Hildreth method. 
//   // For MC, the Hildreth's weights are set up already in the 
//   // constructor. For the data, we set up the source and the target
//   // histogram manually.
//   if( isMC ){
//     // Really, this is already done in the constructor, but just to 
//     // make it clear to a code reader
//     puReweight.setActiveMethod(PUReweight_t::_Hildreth);    
//     // There is no need to set anything else because
//     // the set method above also sets the weight values
//   } else {
//     puReweight.setActiveMethod(PUReweight_t::_TwoHistos);
//     puReweight.setSimpleWeights( puTargetFName, puTargetDistrName,
// 				 puSourceFName, puSourceDistrName);
//   }

  // Print everything
  if( isMC ){
    if( DYTools::energy8TeV ){
      printf("PU reweight info: this is MC. The weights are standard Hildreth weights.\n");
    }else{
      printf("PU reweight info: this is MC. The weights are prepared based on:\n");
      printf("  target (reference) histo= %s   from file= %s\n", 
	     puTargetDistrName.Data(),puTargetFName.Data());
      printf("  source (active)    histo= %s   from file= %s\n", 
	     puSourceDistrName.Data(),puSourceFName.Data());
    }
    for(int i=1; i<=45; i++){
      double ww = 0;
      if( DYTools::energy8TeV ){
	ww = puReweight.getWeightHildreth(i);
      }else{
	ww = puReweight.getWeightTwoHistos(i);
      }
      printf("   PU=%3d     weight= %f\n", i, ww);
    }
  }else{
    printf("PU reweight info: this is data. The PU weights are 1.0\n");
  }

  TFile *selectedEventsFile= new TFile(fName,"UPDATE");
  assert(selectedEventsFile && selectedEventsFile->IsOpen());
  //std::cout << "selectedEventsFile tree entry counts: " << ((TTree*)selectedEventsFile->Get("passTree"))->GetEntries() << ", " << ((TTree*)selectedEventsFile->Get("failTree"))->GetEntries() << "\n";

  // first accumulate the PU distribution
  // loop over the trees (pass/fail)
  for (int fail=0; (fail<2); ++fail) {
    TString pass_fail_str= (fail) ? "_fail" : "_pass";
    TString treeName= (fail) ? "failTree" : "passTree";
    // get the tree
    TTree *tree = (TTree*)selectedEventsFile->Get(treeName); assert(tree);
    // set addresses 
    double evWeight;
    UInt_t nGoodPV;
    tree->SetBranchAddress("eventWeight",&evWeight);   
    tree->SetBranchAddress("nGoodPV",&nGoodPV);
    TBranch *evWeightBr= tree->GetBranch("eventWeight");
    TBranch *pvCountBr= tree->GetBranch("nGoodPV");
    assert(evWeightBr && pvCountBr);

  }

  // set weights
  for (int fail=0; (fail<2); ++fail) {
    TString pass_fail_str= (fail) ? "_fail" : "_pass";
    TString treeName= (fail) ? "failTree" : "passTree";
    // get the tree
    TTree *tree = (TTree*)selectedEventsFile->Get(treeName); assert(tree);
    // set addresses and create a new branch
    double evWeight,pvWeight;
    UInt_t nGoodPV;
    tree->SetBranchAddress("eventWeight",&evWeight);   
    tree->SetBranchAddress("nGoodPV",&nGoodPV);
    tree->Branch("weight",&pvWeight,"weight/D");
    TBranch *evWeightBr= tree->GetBranch("eventWeight");
    TBranch *pvCountBr= tree->GetBranch("nGoodPV");
    TBranch *weightBr= tree->GetBranch("weight");
    assert(evWeightBr && pvCountBr && weightBr);

    // fill the weight branch
    for (UInt_t i=0; i<tree->GetEntries(); ++i) {
      pvCountBr->GetEntry(i);
      evWeightBr->GetEntry(i);
      // For MC, apply PU weights
      if(isMC){
	if( DYTools::energy8TeV ){
	  pvWeight=evWeight * puReweight.getWeightHildreth(nGoodPV);
	}else{
	  pvWeight=evWeight * puReweight.getWeightTwoHistos(nGoodPV);
	}
      }else{
	pvWeight=evWeight;
      }
      weightBr->Fill();
    }
    selectedEventsFile->cd();
    tree->Write();
  }
  selectedEventsFile->Write();
  //std::cout << "selectedEventsFile tree entry counts: " << ((TTree*)selectedEventsFile->Get("passTree"))->GetEntries() << ", " << ((TTree*)selectedEventsFile->Get("failTree"))->GetEntries() << "\n";
  delete selectedEventsFile;
  return 1;
}

// ------------------------------------------------------

#endif
