#include <TCanvas.h>
#include <TPad.h>
#include <TH1.h>
#include <TStyle.h>
#include "MitStyleRemix.hh"

void MitStyleRemix() {
  /*
  const char* author   = "$Author: ikrav $$";
  const char* modified = "$Modified by ksung $$";
  printf(" MIT root style REMIX(%s,%s).\n",author,modified);
  printf("\n");
  printf(" Use: MakeCanvas(name,title)\n");
  printf("      InitSubPad(pad,nPad)\n");
  printf("      InitHist(hist,xTitle,yTitle,color)\n");
  printf("\n");
  */
  if (0) SetStyle();
  else SetTdrStyle(); // "official" CMS style

}

TCanvas* MakeCanvas(const char* name, const char *title, int dX, int dY)
{
  // Start with a canvas
  TCanvas *canvas = new TCanvas(name,title,0,0,dX,dY);
  // General overall stuff
  canvas->SetFillColor      (0);
  canvas->SetBorderMode     (0);
  canvas->SetBorderSize     (10);
  // Set margins to reasonable defaults
  canvas->SetLeftMargin     (0.18);
  canvas->SetRightMargin    (0.05);
  canvas->SetTopMargin      (0.08);
  canvas->SetBottomMargin   (0.15);
  // Setup a frame which makes sense
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);

  return canvas;
}

void InitSubPad(TPad* pad, int i)
{
  //printf("Pad: %p, index: %d\n",pad,i);

  pad->cd(i);
  TPad *tmpPad = (TPad*) pad->GetPad(i);
  tmpPad->SetLeftMargin  (0.18);
  tmpPad->SetTopMargin   (0.05);
  tmpPad->SetRightMargin (0.07);
  tmpPad->SetBottomMargin(0.15);
  return;
}

void InitHist(TH1 *hist, const char *xtit, const char *ytit, EColor color)
{
  hist->SetXTitle(xtit);
  hist->SetYTitle(ytit);
  hist->SetLineColor(color);
  hist->SetTitleSize  (0.055,"Y");
  hist->SetTitleOffset(1.600,"Y");
  hist->SetLabelOffset(0.014,"Y");
  hist->SetLabelSize  (0.050,"Y");
  hist->SetLabelFont  (42   ,"Y");
  hist->SetTitleSize  (0.055,"X");
  hist->SetTitleOffset(1.300,"X");
  hist->SetLabelOffset(0.014,"X");
  hist->SetLabelSize  (0.050,"X");
  hist->SetLabelFont  (42   ,"X");
  hist->SetMarkerStyle(20);
  hist->SetMarkerColor(color);
  hist->SetMarkerSize (0.6);
  // Strangely enough this cannot be set anywhere else??
  hist->GetYaxis()->SetTitleFont(42);
  hist->GetXaxis()->SetTitleFont(42);
  hist->SetTitle("");
  return;
}

void SetStyle()
{
  TStyle *MITStyle = new TStyle("MIT-Style","The Perfect Style for Plots ;-)");
  gStyle = MITStyle;

  // Canvas
  MITStyle->SetCanvasColor     (0);
  MITStyle->SetCanvasBorderSize(10);
  MITStyle->SetCanvasBorderMode(0);
  MITStyle->SetCanvasDefH      (700);
  MITStyle->SetCanvasDefW      (700);
  MITStyle->SetCanvasDefX      (100);
  MITStyle->SetCanvasDefY      (100);

  // Pads
  MITStyle->SetPadColor       (0);
  MITStyle->SetPadBorderSize  (10);
  MITStyle->SetPadBorderMode  (0);
  MITStyle->SetPadBottomMargin(0.13);
  MITStyle->SetPadTopMargin   (0.08);
  MITStyle->SetPadLeftMargin  (0.18);
  MITStyle->SetPadRightMargin (0.05);
  MITStyle->SetPadGridX       (0);
  MITStyle->SetPadGridY       (0);
  MITStyle->SetPadTickX       (0);
  MITStyle->SetPadTickY       (0);

  // Frames
  MITStyle->SetFrameFillStyle ( 0);
  MITStyle->SetFrameFillColor ( 0);
  MITStyle->SetFrameLineColor ( 1);
  MITStyle->SetFrameLineStyle ( 0);
  MITStyle->SetFrameLineWidth ( 1);
  MITStyle->SetFrameBorderSize(10);
  MITStyle->SetFrameBorderMode( 0);

  // Histograms
  MITStyle->SetHistFillColor(2);
  MITStyle->SetHistFillStyle(0);
  MITStyle->SetHistLineColor(1);
  MITStyle->SetHistLineStyle(0);
  MITStyle->SetHistLineWidth(2);
  MITStyle->SetNdivisions(505);

  // Functions
  MITStyle->SetFuncColor(1);
  MITStyle->SetFuncStyle(0);
  MITStyle->SetFuncWidth(2);

  // Various
  MITStyle->SetMarkerStyle(20);
  MITStyle->SetMarkerColor(kBlack);
  MITStyle->SetMarkerSize (1.2);

  MITStyle->SetTitleBorderSize(0);
  MITStyle->SetTitleFillColor (0);
  MITStyle->SetTitleX         (0.2);

  MITStyle->SetTitleSize  (0.055,"X");
  MITStyle->SetTitleOffset(1.200,"X");
  MITStyle->SetLabelOffset(0.005,"X");
  MITStyle->SetLabelSize  (0.050,"X");
  MITStyle->SetLabelFont  (42   ,"X");

  MITStyle->SetStripDecimals(kFALSE);

  MITStyle->SetTitleSize  (0.055,"Y");
  MITStyle->SetTitleOffset(1.600,"Y");
  MITStyle->SetLabelOffset(0.010,"Y");
  MITStyle->SetLabelSize  (0.050,"Y");
  MITStyle->SetLabelFont  (42   ,"Y");

  MITStyle->SetTextSize   (0.055);
  MITStyle->SetTextFont   (42);

  MITStyle->SetStatFont   (42);
  MITStyle->SetTitleFont  (42);
  MITStyle->SetTitleFont  (42,"X");
  MITStyle->SetTitleFont  (42,"Y");

  MITStyle->SetOptStat    (0);
  return;
}


// Probably "official" CMS style
void SetTdrStyle() {

  TStyle *tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

  printf("TDR Style initialized");

// For the canvas:
  tdrStyle->SetCanvasBorderMode(0);
  tdrStyle->SetCanvasColor(kWhite);
  tdrStyle->SetCanvasDefH(600); //Height of canvas
  tdrStyle->SetCanvasDefW(600); //Width of canvas
  tdrStyle->SetCanvasDefX(0);   //POsition on screen
  tdrStyle->SetCanvasDefY(0);

// For the Pad:
  tdrStyle->SetPadBorderMode(0);
  // tdrStyle->SetPadBorderSize(Width_t size = 1);
  tdrStyle->SetPadColor(kWhite);
  tdrStyle->SetPadGridX(false);
  tdrStyle->SetPadGridY(false);
  tdrStyle->SetGridColor(0);
  tdrStyle->SetGridStyle(3);
  tdrStyle->SetGridWidth(1);

// For the frame:
  tdrStyle->SetFrameBorderMode(0);
  tdrStyle->SetFrameBorderSize(1);
  tdrStyle->SetFrameFillColor(0);
  tdrStyle->SetFrameFillStyle(0);
  tdrStyle->SetFrameLineColor(1);
  tdrStyle->SetFrameLineStyle(1);
  tdrStyle->SetFrameLineWidth(1);

// For the histo:
  // tdrStyle->SetHistFillColor(1);
  // tdrStyle->SetHistFillStyle(0);
  tdrStyle->SetHistLineColor(1);
  tdrStyle->SetHistLineStyle(0);
  tdrStyle->SetHistLineWidth(1);
  // tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
  // tdrStyle->SetNumberContours(Int_t number = 20);

  tdrStyle->SetEndErrorSize(2);
  //tdrStyle->SetErrorMarker(20);
  tdrStyle->SetErrorX(0.);

  tdrStyle->SetMarkerStyle(20);

//For the fit/function:
  tdrStyle->SetOptFit(1);
  tdrStyle->SetFitFormat("5.4g");
  tdrStyle->SetFuncColor(2);
  tdrStyle->SetFuncStyle(1);
  tdrStyle->SetFuncWidth(1);

//For the date:
  tdrStyle->SetOptDate(0);
  // tdrStyle->SetDateX(Float_t x = 0.01);
  // tdrStyle->SetDateY(Float_t y = 0.01);

// For the statistics box:
  tdrStyle->SetOptFile(0);
  tdrStyle->SetOptStat(0); // To display the mean and RMS:   SetOptStat("mr");
  tdrStyle->SetStatColor(kWhite);
  tdrStyle->SetStatFont(42);
  tdrStyle->SetStatFontSize(0.025);
  tdrStyle->SetStatTextColor(1);
  tdrStyle->SetStatFormat("6.4g");
  tdrStyle->SetStatBorderSize(1);
  tdrStyle->SetStatH(0.1);
  tdrStyle->SetStatW(0.15);
  // tdrStyle->SetStatStyle(Style_t style = 1001);
  // tdrStyle->SetStatX(Float_t x = 0);
  // tdrStyle->SetStatY(Float_t y = 0);

// Margins:
  tdrStyle->SetPadTopMargin(0.15);
  tdrStyle->SetPadBottomMargin(0.13);
  tdrStyle->SetPadLeftMargin(0.18);
  tdrStyle->SetPadRightMargin(0.15);

// For the Global title:

//  tdrStyle->SetOptTitle(0);
  tdrStyle->SetTitleFont(42);
  tdrStyle->SetTitleColor(1);
  tdrStyle->SetTitleTextColor(1);
  tdrStyle->SetTitleFillColor(10);
  tdrStyle->SetTitleFontSize(0.05);
  // tdrStyle->SetTitleH(0); // Set the height of the title box
  // tdrStyle->SetTitleW(0); // Set the width of the title box
  tdrStyle->SetTitleX(0.255); // Set the position of the title box
  tdrStyle->SetTitleY(0.985); // Set the position of the title box
  // tdrStyle->SetTitleStyle(Style_t style = 1001);
  tdrStyle->SetTitleBorderSize(0); // no border around the title

// For the axis titles:

  tdrStyle->SetTitleColor(1, "XYZ");
  tdrStyle->SetTitleFont(42, "XYZ");
  tdrStyle->SetTitleSize(0.06, "XYZ");
  // tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
  // tdrStyle->SetTitleYSize(Float_t size = 0.02);
  tdrStyle->SetTitleXOffset(1.0);
//   tdrStyle->SetTitleYOffset(1.5);
  tdrStyle->SetTitleOffset(1.5, "Y"); // Another way to set the Offset

// For the axis labels:

  tdrStyle->SetLabelColor(1, "XYZ");
  tdrStyle->SetLabelFont(42, "XYZ");
  tdrStyle->SetLabelOffset(0.007, "XYZ");
  tdrStyle->SetLabelSize(0.05, "XYZ");

// For the axis:

  tdrStyle->SetAxisColor(1, "XYZ");
  tdrStyle->SetStripDecimals(false);
  tdrStyle->SetTickLength(0.03, "XYZ");
  //tdrStyle->SetNdivisions(510, "XYZ");
  tdrStyle->SetNdivisions(505, "XYZ");
  tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
  tdrStyle->SetPadTickY(1);

// Change for log plots:
  tdrStyle->SetOptLogx(0);
  tdrStyle->SetOptLogy(0);
  tdrStyle->SetOptLogz(0);

  tdrStyle->SetPalette(1,0);

// Postscript options:
  // tdrStyle->SetPaperSize(15.,15.);
  // tdrStyle->SetLineScalePS(Float_t scale = 3);
  // tdrStyle->SetLineStyleString(Int_t i, const char* text);
  // tdrStyle->SetHeaderPS(const char* header);
  // tdrStyle->SetTitlePS(const char* pstitle);

  // tdrStyle->SetBarOffset(Float_t baroff = 0.5);
  // tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
  // tdrStyle->SetPaintTextFormat(const char* format = "g");
  // tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
  // tdrStyle->SetTimeOffset(Double_t toffset);
  // tdrStyle->SetHistMinimumZero(kTRUE);

  tdrStyle->cd();
}
