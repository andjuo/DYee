#include "../Inc_7TeV/FEWZ.hh"
#include <TFile.h>
#include <TString.h>

// ----------------------------------------------------------

FEWZ_t::FEWZ_t(bool loadWeights, bool do_cutZPT100, bool regularize) : 
  fInitialized(kFALSE), fCutZPT100(do_cutZPT100),
  fRegularized(regularize)
{
  if (loadWeights) {
    TFile fweights("../root_files7TeV/fewz/weights_stepwise_prec10-5_fine12.root");
    if(do_cutZPT100)
      cout << "FEWZ_t NOTE: in MC, for Z/gamma* PT>100 the FEWZ weights for 80<PT<100 GeV are used!" << endl;
    if( !fweights.IsOpen() ) assert(0);
    bool ok=kTRUE;
    for(int i=0; ok && (i<_nMassBinsFEWZ); i++){
      TString hnames = TString::Format("weight_%02d",i+1);
      weights[i] = (TH2D*)fweights.Get(hnames);
      hnames = TString::Format("h_weighterror_%02d",i+1);
      weightErrors[i] = (TH2D*)fweights.Get(hnames);
      if (!weights[i] || !weightErrors[i]) ok=kFALSE;
      else {
	weights[i]->SetDirectory(0);
	weightErrors[i]->SetDirectory(0);
      }
    }
    fInitialized=ok;

    if (ok) {
      if (fRegularized) {
	std::cout << "Regularization weights[7] (" << weights[7]->GetTitle() << ") scale=" << (1.021/0.986) << ", etc.\n";
	weights     [7]->Scale(1.021/0.986);
	weightErrors[7]->Scale(1.021/0.986);
	weights     [8]->Scale(1.024/1.05 );
	weightErrors[8]->Scale(1.024/1.05 );
	weights     [9]->Scale(1.027/1.003);
	weightErrors[9]->Scale(1.027/1.003);
	weights     [10]->Scale(1.030/0.985);
	weightErrors[10]->Scale(1.030/0.985);
	weights     [14]->Scale(1.062/1.052);
	weightErrors[14]->Scale(1.062/1.052);
      }
    }
  }
}

// ----------------------------------------------------------
