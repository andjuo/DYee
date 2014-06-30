#ifndef MitStyleRemix_HH
#define MitStyleRemix_HH

#include <TCanvas.h>
#include <TH1.h>
#include <TPad.h>

void     MitStyleRemix();
TCanvas* MakeCanvas   (const char* name, const char *title,
		       int dX = 500, int dY = 500);
void     InitSubPad   (TPad* pad, int i);
void     InitHist     (TH1 *hist, const char *xtit,
		       const char *ytit  = "Number of Entries",
		       EColor color = kBlack);
void     SetStyle     ();
void     SetTdrStyle  ();  // probably the "official" CMS style



#ifndef ColorPalettes_HH
inline
void AdjustFor2DplotWithHeight(TCanvas *c, double rmargin=0.18,
			       int logx=-1, int logy=-1) {
  int count=0;
  for (int i=0; i<50; i++) {
    TPad *pad=(TPad*)c->GetPad(i);
    if (pad) {
      count++;
      pad->SetRightMargin(rmargin);
      if (logx!=-1) pad->SetLogx(logx);
      if (logy!=-1) pad->SetLogy(logy);
    }
  }
  if (count==1) {
    c->SetRightMargin(rmargin);
    if (logx!=-1) c->SetLogx(logx);
    if (logy!=-1) c->SetLogy(logy);
  }
}
#endif


// ------------------------------------------------------------

// A function to control space that histogram takes.
// E.g. if the y-labels+title need more space call SetSideSpaces(c,0.2,0,0,0);

#ifndef ColorPalettes_HH
inline
void SetSideSpaces(TCanvas *c, double dxLeft=0.05, double dxRight=0.02,
		   double dyTop=0., double dyBottom=0.02,
		   int padIdx=-1) {
  int count=0;
  for (int i=0; i<50; i++) {
    if ((padIdx!=-1) && (i!=padIdx)) continue;
    TPad *pad=(TPad*)c->GetPad(i);
    if (pad) {
      count++;
      pad->Range(pad->GetX1()-dxLeft,pad->GetY1()-dyBottom,
		 pad->GetX2()-dxRight,pad->GetY2()-dyTop);
      pad->SetLeftMargin(pad->GetLeftMargin() + dxLeft);
      pad->SetRightMargin(pad->GetRightMargin() + dxRight);
      pad->SetTopMargin(pad->GetTopMargin() + dyTop);
      pad->SetBottomMargin(pad->GetBottomMargin() + dyBottom);
    }
  }
}
#endif


#endif
