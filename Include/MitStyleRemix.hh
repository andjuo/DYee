#ifndef MitStyleRemix_HH
#define MitStyleRemix_HH

#include <TCanvas.h>
#include <TH1.h>
#include <TPad.h>

void     MitStyleRemix();
TCanvas* MakeCanvas   (const char* name, const char *title, int dX = 500, int dY = 500);
void     InitSubPad   (TPad* pad, int i);
void     InitHist     (TH1 *hist, const char *xtit, const char *ytit  = "Number of Entries",
		       EColor color = kBlack);
void     SetStyle     ();



inline
void AdjustFor2DplotWithHeight(TCanvas *c, double rmargin=0.18) {
  int count=0;
  for (int i=0; i<50; i++) {
    TPad *pad=(TPad*)c->GetPad(i);
    if (pad) {
      count++;
      pad->SetRightMargin(rmargin);
    }
  }
  if (count==1) c->SetRightMargin(rmargin);
}

#endif
