#include "../Include/AccessOrigNtuples.hh"


// ---------------------
// ---------------------

void print(const char *msg, const TLorentzVector &v) {
  if (msg) std::cout << msg;
  std::cout << " (Pt,eta,phi)="
	    << Form("(%f,%f,%f)",v.Pt(),v.Eta(),v.Phi()) << "\n";
}

// ---------------------
// ---------------------

TLorentzVector* AccessOrigNtuples_t::getDressedGenDielectron(double dRthr,
	 TLorentzVector *dressedE1_ext, TLorentzVector *dressedE2_ext) const
{
  if (!genPhBrIsActive) {
    std::cout <<"error detected in AccessOrigNtuples::getDressedGenDielectron:"
	      <<" genPhBrIsActive=0\n";
    return NULL;
  }
  // E1 will be an electron (PdgId=11), E2 will be a positron (PdgId=-11)
  TLorentzVector *dressedE1= dressedE1_ext;
  TLorentzVector *dressedE2= dressedE2_ext;
  if (!dressedE1) dressedE1= new TLorentzVector();
  if (!dressedE2) dressedE2= new TLorentzVector();
  if ((gen->lid_1 + gen->lid_2!=0) || (abs(gen->lid_1)!=11)) {
    std::cout << "getDressedGenDielectron: wrong lid_x info\n";
    return NULL;
  }
  dressedE1->SetPtEtaPhiM(gen->pt_1, gen->eta_1, gen->phi_1, 0.000511);
  dressedE2->SetPtEtaPhiM(gen->pt_2, gen->eta_2, gen->phi_2, 0.000511);
  if (gen->lid_1!=11) {
    std::cout << "getDressedGenDielectron: swapping values\n";
    TLorentzVector tmp=*dressedE1; *dressedE1=*dressedE2; *dressedE2=tmp;
  }

  const int debug=0;
  if (debug) std::cout << "\n";

  int nGenPhotons= genPhArr->GetEntries();
  TLorentzVector ph;
  for (Int_t i=0; i<nGenPhotons; i++) {
    const mithep::TGenPhoton *gph=
                           (const mithep::TGenPhoton*)((*genPhArr)[i]);
    TLorentzVector *mother= NULL;
    if (gph->parentId == gen->lid_1) mother=dressedE1;
    else if (gph->parentId == gen->lid_2) mother=dressedE2;
    //else std::cout << "gph->parentId=" << gph->parentId << "\n";
    if (mother) {
      ph.SetPtEtaPhiM(gph->pt, gph->eta, gph->phi, 0.);
      double dr=mother->DeltaR(ph);
      if (debug) {
	int idx= (gph->parentId == gen->lid_1) ? 1:2;
	std::cout << "i=" << i << ", dr from mother" << idx
		  << " = " << dr << "\n";
      }
      if (dr<dRthr) {
	(*mother) += ph;
      }
    }
  }

  TLorentzVector *dressedDiE=
                    new TLorentzVector( (*dressedE1) + (*dressedE2) );

  if (debug && (fabs(gen->vpt_1-dressedE1->Pt())>0.01)) {
    std::cout << "postFsr e1: (Pt,eta,phi)="
	      << Form("(%f,%f,%f)",gen->pt_1,gen->eta_1,gen->phi_1) << "\n";
    std::cout << "preFsr e1: (Pt,eta,phi)="
	      << Form("(%f,%f,%f)",gen->vpt_1,gen->veta_1,gen->vphi_1) << "\n";
    print("dressedE1 ",*dressedE1);
    std::cout << "postFsr e2: (Pt,eta,phi)="
	      << Form("(%f,%f,%f)",gen->pt_2,gen->eta_2,gen->phi_2) << "\n";
    std::cout << "preFsr e2: (Pt,eta,phi)="
	      << Form("(%f,%f,%f)",gen->vpt_2,gen->veta_2,gen->vphi_2) << "\n";
    print("dressedE2 ",*dressedE2);
  }

  if (!dressedE1_ext) delete dressedE1;
  if (!dressedE2_ext) delete dressedE2;

  return dressedDiE;
}

// ---------------------
