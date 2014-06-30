#include "helpers.hh"


// -------------------------------------------------------

TH2D *loadTextFile(TString fname, TString histoName, int hasError,
		   int hasDivider) {

  std::ifstream fin(fname);
  if (!fin.is_open()) {
    std::cout << "failed to open the file <" << fname << ">\n";
    return NULL;
  }

  TH2D *h2=createBaseH2(histoName,histoName,1);

  std::string s, tmp;
  double val=0., err=0.;

  int ibin=1;
  int jbin=1;
  while (!fin.eof() && getline(fin,s)) {
    std::stringstream ss(s);
    ss >> val;
    if (hasDivider) ss >> tmp;
    if (hasError)   ss>> err;
    h2->SetBinContent(ibin,jbin, val);
    h2->SetBinError  (ibin,jbin, err);
    ibin++;
    if (ibin>DYTools::nMassBins) {
      ibin=1;
      jbin++;
    }
  }
  fin.close();
  return h2;
}

// -------------------------------------------------------

TH2D *loadXYZTextFile(TString fname, TString histoName, int hasError,
		      int transposed) {

  std::vector<std::string> lines;

  std::ifstream fin(fname);
  if (!fin.is_open()) {
    std::cout << "failed to open the file <" << fname << ">\n";
    return NULL;
  }
  lines.reserve(1000);
  std::string s;
  int xMax=0, yMax=0;
  int x,y;
  while (!fin.eof() && getline(fin,s)) {
    if (lines.size()==lines.capacity()) {
      std::cout << "increate the capacity\n";
      lines.reserve(lines.size() + 1000);
    }
    lines.push_back(s);
    std::stringstream ss(s);
    ss >> x >> y;
    if (x>xMax) xMax=x;
    if (y>yMax) yMax=y;
  }
  fin.close();
  if (transposed) { int tmp=xMax; xMax=yMax; yMax=tmp; }

  std::cout << "xMax=" << xMax << ", yMax=" << yMax << "\n";

  TH2D *h2= new TH2D(histoName,histoName, xMax,0.,xMax,yMax,0.,yMax);
  h2->SetDirectory(0);
  h2->SetStats(0);
  h2->Sumw2();

  double val=0., err=0.;

  int ibin=1;
  int jbin=1;
  for (unsigned int i=0; i<lines.size(); ++i) {
    std::stringstream ss(lines[i]);
    ss >> ibin >> jbin >> val;
    if (hasError) ss >> err;
    if (transposed) { int tmp=ibin; ibin=jbin; jbin=tmp; }
    h2->SetBinContent(ibin,jbin, val);
    h2->SetBinError  (ibin,jbin, err);
  }
  return h2;
}

// --------------------------------------------------------

TH2D *loadHisto1D_convert_TH2D(TString fname, TString fieldName) {
  TH1D* htmp=LoadHisto1D(fieldName,fname,"",0);
  if (!htmp) return NULL;
  htmp->SetName("htmp");
  TH2D *h2=createBaseH2(fieldName,fieldName,1);
  if (!h2) return NULL;
  for (int ibin=1; ibin<=htmp->GetNbinsX(); ++ibin) {
    h2->SetBinContent(ibin,1, htmp->GetBinContent(ibin));
    h2->SetBinError  (ibin,1, htmp->GetBinError(ibin));
  }
  return h2;
}


// --------------------------------------------------------

void printProfileSums(const TH2D* h2,TVectorD *sumX_user,TVectorD *sumY_user) {

  TVectorD sumX(h2->GetNbinsX());
  TVectorD sumY(h2->GetNbinsY());
  sumX.Zero(); sumY.Zero();

  for (int ibin=1; ibin<=h2->GetNbinsX(); ++ibin) {
    double sumI=0.;
    for (int jbin=1; jbin<=h2->GetNbinsY(); ++jbin) {
      double x=h2->GetBinContent(ibin,jbin);
      sumI+=x;
      sumY[jbin-1] += x;
    }
    sumX[ibin-1] =sumI;
  }
  std::cout << "profile sums for " << h2->GetName() << "\n";
  sumX.Print();
  sumY.Print();

  if (sumX_user) (*sumX_user) = sumX;
  if (sumY_user) (*sumY_user) = sumY;
  return;
}


// --------------------------------------------------------

void compareProfileSums(const TH2D *h2A, const TH2D *h2B) {
  int dimX=h2A->GetNbinsX();
  int dimY=h2A->GetNbinsY();
  TVectorD sumXa(dimX), sumXb(dimX);
  TVectorD sumYa(dimY), sumYb(dimY);
  printProfileSums(h2A,&sumXa,&sumYa);
  printProfileSums(h2B,&sumXb,&sumYb);

  TVectorD diffX(sumXa);
  diffX -= sumXb;
  TVectorD diffY(sumYa);
  diffY -= sumYb;
  std::cout << "differences in profiles of " << h2A->GetName()
	    << " and " << h2B->GetName() << "\n";
  diffX.Print();
  diffY.Print();
}


// --------------------------------------------------------
