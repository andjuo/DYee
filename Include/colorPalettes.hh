#ifndef ColorPalettes_HH
#define ColorPalettes_HH

#include <TROOT.h>
#include <TColor.h>
#include <TStyle.h>

// ------------------------------------------------------------
// ------------------------------------------------------------


typedef enum { _colrange_none=0, _colrange_center=1, _colrange_positive=2, 
	       _colrange_positive_bw,
	       _colrange_three, _colrange_nice, 
	       _colrange_default,
	       _colrange_last } TColorRange_t;

// ------------------------------------------------------------
// ------------------------------------------------------------

void set_nice_style(int nb=51) {
  const Int_t NRGBs = 5;
  const Int_t NCont = nb;
   Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
   Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
   Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
   Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
   TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
   gStyle->SetNumberContours(NCont);
}


// ------------------------------------------------------------

void set_three_style(int nb=51) {
  const Int_t NRGBs = 3;
  const Int_t NCont = nb;
   Double_t stops[NRGBs] = { 0.00, 0.50, 1.00 };
   Double_t red[NRGBs]   = { 0.10, 0.10, 1.00 };
   Double_t green[NRGBs] = { 0.10, 0.80, 0.00 };
   Double_t blue[NRGBs]  = { 1.00, 0.10, 0.00 };
   TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
   gStyle->SetNumberContours(NCont);
}


// ------------------------------------------------------------

void set_bottom_white_style(int nb=51) {
   const UInt_t Number = 3;
   Double_t Red[Number]    = { 1.00, 1.00, 1.00};
   Double_t Green[Number]  = { 1.00, 0.50, 0.00};
   Double_t Blue[Number]   = { 1.00, 0.50, 0.00};
   Double_t Length[Number] = { 0.00, 0.50, 1.00 };
   TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
   gStyle->SetNumberContours(nb);
}

// ------------------------------------------------------------

void set_bottom_white_style_bw(int nb=51) {
   const UInt_t Number = 3;
   Double_t Red[Number]    = { 0.50, 0.750, 1.00};
   Double_t Green[Number]  = { 0.50, 0.750, 1.00};
   Double_t Blue[Number]   = { 0.50, 0.750, 1.00};
   Double_t Length[Number] = { 0.00, 0.750, 1.00 };
   TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
   gStyle->SetNumberContours(nb);
}

// ------------------------------------------------------------

void set_wbr_style(int nb=51, double blue_pos=0.25) {
  const UInt_t Number = 3;
   Double_t Red[Number]    = { 1.00, 0.000, 1.00};
   Double_t Green[Number]  = { 1.00, 0.000, 0.00};
   Double_t Blue[Number]   = { 1.00, 1.000, 0.00};
   Double_t Length[Number] = { 0.00, blue_pos, 1.00 };
   TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
   gStyle->SetNumberContours(nb);
}

// ------------------------------------------------------------

void set_center_white_style(int nb=51) {
   const UInt_t Number = 3;
   Double_t Red[Number]    = { 0.00, 1.00, 0.90};
   Double_t Green[Number]  = { 0.00, 1.00, 0.00};
   Double_t Blue[Number]   = { 0.50, 1.00, 0.00};
   Double_t Length[Number] = { 0.00, 0.50, 1.00 };
   TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
   gStyle->SetNumberContours(nb);
}

// ------------------------------------------------------------

void set_center_white_style5(int nb=51) {
   const UInt_t Number = 5;
   Double_t Red[Number]    = { 0.00, 0.00, 1.00, 0.90, 0.90};
   Double_t Green[Number]  = { 0.50, 0.00, 1.00, 0.00, 0.90};
   Double_t Blue[Number]   = { 0.30, 0.60, 1.00, 0.00, 0.30};
   Double_t Length[Number] = { 0.00, 0.25, 0.50, 0.75, 1.00 };
   TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb,0.9);
   gStyle->SetNumberContours(nb);
}

// ------------------------------------------------------------

void set_center_white_style2(int nb=51) {
   const UInt_t Number = 3;
   Double_t Red[Number]    = { 0.00, 1.00, 1.00};
   Double_t Green[Number]  = { 0.00, 1.00, 0.00};
   Double_t Blue[Number]   = { 0.50, 1.00, 1.00};
   Double_t Length[Number] = { 0.00, 0.50, 1.00 };
   TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
   gStyle->SetNumberContours(nb);
}

// ------------------------------------------------------------
// blue is less dark

void set_center_white_style_mod(int nb=51, double rMin=0.40, double gMin=0.20, double bMin=0.8) {
   const UInt_t Number = 3;
   Double_t Red[Number]    = { rMin, 1.00, 1.00};
   Double_t Green[Number]  = { gMin, 1.00, 0.00};
   Double_t Blue[Number]   = { bMin, 1.00, 1.00};
   Double_t Length[Number] = { 0.00, 0.50, 1.00 };
   TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
   gStyle->SetNumberContours(nb);
}

// ------------------------------------------------------------

void set_center_white_style_half(int nb=51) {
   const UInt_t Number = 2;
   Double_t Red[Number]    = { 1.00, 1.00};
   Double_t Green[Number]  = { 1.00, 0.00};
   Double_t Blue[Number]   = { 1.00, 1.00};
   Double_t Length[Number] = { 0.00, 1.00 };
   TColor::CreateGradientColorTable(Number,Length,Red,Green,Blue,nb);
   gStyle->SetNumberContours(nb);
}

// ------------------------------------------------------------

void set_white_style_palette(TColorRange_t colRange, int nColorBins=51) {
  switch(colRange) {
  case _colrange_none: ; break;
  case _colrange_center: set_center_white_style(nColorBins); break;
  case _colrange_positive: set_bottom_white_style(nColorBins); break;
  case _colrange_three: set_three_style(nColorBins); break;
  case _colrange_nice: set_nice_style(nColorBins); break;
  case _colrange_default: gStyle->SetPalette(1); break;
  default: 
    std::cout << "\n\n\tError in set_white_style_palette: not ready for int(colRange)=" << int(colRange) << "\n\n";
  }
  return;
}

// ------------------------------------------------------------

template<class Histo2D_t>
void getZrange(const Histo2D_t *h, double &zmin, double &zmax) {

  zmin=1e9; zmax=-1e9;
  for (int i=1; i<=h->GetNbinsX(); ++i) {
    //std::cout << "i=" << i << "\n";
    for (int j=1; j<=h->GetNbinsY(); ++j) {
      double v=h->GetBinContent(i,j);
      if (v<zmin) zmin=v;
      if (v>zmax) {
	//std::cout << "i=" << i << ", j=" << j << ", value=" << v << "\n";
	zmax=v;
      }
    }
  }
  return;
}

// ------------------------------------------------------------

template<class Histo2D_t>
void printZrange(const Histo2D_t *h) {
  double zmin=0, zmax=0.;
  getZrange(h,zmin,zmax);
  std::cout << "histo(" << h->GetName() << ")  zrange=" 
	    << zmin << " .. " << zmax << "\n";
}

// ------------------------------------------------------------

int centerHistoZAxis(TH2D *h, TColorRange_t centerRange, double maxValUser=0.) {
  if (centerRange <= _colrange_none) return 1;

  double zmin=1e9, zmax=-1e9;
  for (int i=1; i<=h->GetNbinsX(); ++i) {
    //std::cout << "i=" << i << "\n";
    for (int j=1; j<=h->GetNbinsY(); ++j) {
      double v=h->GetBinContent(i,j);
      if (v<zmin) zmin=v;
      if (v>zmax) zmax=v;
    }
  }

  if (centerRange > _colrange_none) {
    std::cout << h->GetName() << ": zmin=" << zmin << ", zmax=" << zmax << "\n";
    if ((centerRange==_colrange_center) ||
	(centerRange==_colrange_three)) {
      if (zmin*zmax >= 0) {
	std::cout << " warning in centerHistoZAxis: cannot center (zmin=" 
		  << zmin << ", zmax=" << zmax << ")\n";
      }

      {
	double z=zmax;
	if (-zmin > z) z= -zmin;
	std::cout << "histoName=<" << h->GetName() << ">, max|z|= " << z << "\n";
	if (maxValUser>0) {
	  z=maxValUser;
	  std::cout << " maxValUser=" << maxValUser << "\n";
	}
	h->GetZaxis()->SetRangeUser(-z,z);
      }
    }
    else if ((centerRange==_colrange_positive) ||
	     (centerRange==_colrange_positive_bw)) {
      if (maxValUser>0.) {
	if (maxValUser < zmax) {
	  std::cout << "histoName=<" << h->GetName() << ">, warning: zmax=" << zmax << ", maxValUser=" << maxValUser << "\n"; 
	}
	h->GetZaxis()->SetRangeUser(0,maxValUser);
      }
    }
  }

  return 1;
}

// ------------------------------------------------------------
// ------------------------------------------------------------

int drawHistoSubpad(TCanvas *c, int subPad, TH2D *h2, TColorRange_t centerRange, int useMassBins, int nColorBins=51) {
  c->cd(subPad);
  TPad *pad=(TPad*)c->GetPad(subPad);

  if (useMassBins) {
    pad->SetLogy(1);
    pad->SetLogx(1);

    h2->GetXaxis()->SetMoreLogLabels();
    h2->GetXaxis()->SetNoExponent();
    h2->GetYaxis()->SetMoreLogLabels();
    h2->GetYaxis()->SetNoExponent();
    int ndiv=703;
    h2->GetXaxis()->SetNdivisions(ndiv,true);
    h2->GetYaxis()->SetNdivisions(ndiv,true);
  }
  switch (centerRange) {
  case _colrange_center: set_center_white_style(nColorBins); break;
  case _colrange_positive: set_bottom_white_style(nColorBins); break;
  case _colrange_positive_bw: set_bottom_white_style_bw(nColorBins); break;
  default: ; 
  }
  h2->Draw("colz");
  //if (useMassBins) {
  //  h2->GetXaxis()->SetRangeUser(15,4000.);
  //  h2->GetYaxis()->SetRangeUser(15,4000.);
  //}
  c->Modified();
  return 1;
}

// ------------------------------------------------------------

int drawHisto(TCanvas *c, TH2D *h2, TColorRange_t centerRange, int useMassBins, int nColorBins=51) {
  int res=drawHistoSubpad(c,0,h2,centerRange,useMassBins,nColorBins);
  c->Update();
  return res;
}

// ------------------------------------------------------------
// ------------------------------------------------------------

int drawHistoSubpadAdjustZ(TCanvas *c, int subPad, TH2D *h2, TColorRange_t centerRange, int useMassBins, double maxValUser=0., int nColorBins=51) {
  int res=(centerHistoZAxis(h2,centerRange,maxValUser) &&
	   drawHistoSubpad(c,subPad,h2,centerRange,useMassBins,nColorBins)) ? 1:0;
  return res;
}

// ------------------------------------------------------------

int drawHistoAdjustZ(TCanvas *c, TH2D *h2, TColorRange_t centerRange, int useMassBins, double maxValUser=0., int nColorBins=51) {
  return drawHistoSubpadAdjustZ(c,0,h2,centerRange,useMassBins,maxValUser,nColorBins);
}

// ------------------------------------------------------------
// ------------------------------------------------------------

#ifndef MitStyleRemix_HH
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

#ifndef MitStyleRemix_HH
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

// ------------------------------------------------------------


#endif
