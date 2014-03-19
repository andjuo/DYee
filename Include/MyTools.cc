#include "../Include/MyTools.hh"

//--------------------------------------------------
//--------------------------------------------------
/*
TH2D* getRelDifferenceVA(const TH2D *baseValue, TString newName, int nVariations, TH2D* hVar1, ...) {
  HERE("entered getRelDifference");
  TString newTitle=baseValue->GetTitle() + TString(" rel.Diff");
  TH2D *hRes=Clone(baseValue,newName,newTitle);
  HERE("ok");

  // remove error
  for (int ibin=1; ibin<=hRes->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=hRes->GetNbinsY(); ++jbin) {
      hRes->SetBinError(ibin,jbin, 0.);
    }
  }
  HERE("start variations");

  if (nVariations>=1) {
    va_list vl;
    va_start(vl,hVar1);
    for (int i=0; i<nVariations; ++i) {
      HERE("i=%d",i);
      TH2D* hv_inp=(TH2D*)va_arg(vl,TH2D*);
      if (!hv_inp) HERE("hv_inp is null");
      HERE("got hv_inp");
      //printHisto(hv_inp);
      hv_inp->Print("range");  // crashed here!
      TH2D* hv=(TH2D*)hv_inp->Clone(Form("var_%d",i));
      HERE("cloned ok %d",((hv==NULL)?0:1));
      hv->Add(baseValue,-1.);
      HERE("added base val");
      for (int ibin=1; ibin<=hRes->GetNbinsX(); ++ibin) {
	for (int jbin=1; jbin<=hRes->GetNbinsY(); ++jbin) {
	  double err=hRes->GetBinError(ibin,jbin);
	  double var=fabs(hv->GetBinError(ibin,jbin));
	  if (var>err) hRes->SetBinError(ibin,jbin, var);
	}
      }
      delete hv;
    }
    va_end(vl);
  }
  return hRes;
}
*/

//--------------------------------------------------

TH2D* getRelDifference(const TH2D *baseValue, TString newName, int includeVariants, const TH2D* hVar1, const TH2D *hVar2, const TH2D *hVar3, const TH2D* hVar4) {
  if (!baseValue) {
    std::cout << "getRelDifference: base histogram is NULL" << std::endl;
    return NULL;
  }

  // TString newTitle=baseValue->GetTitle() + TString(" rel.Diff");
  TH2D *hRes=Clone(baseValue,newName,newName);
  if (!hRes) {
    std::cout << "getRelDifference: failed to copy the histogram" << std::endl;
    return NULL;
  }

  // remove error
  for (int ibin=1; ibin<=hRes->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=hRes->GetNbinsY(); ++jbin) {
      hRes->SetBinError(ibin,jbin, 0.);
    }
  }

  if (!hVar1) {
    std::cout << "\n\n\t\tERROR getRelDifference: hVar1 is null\n" << std::endl;
    return hRes;
  }

  std::vector<const TH2D*> hV;
  hV.reserve(4);
  hV.push_back(hVar1);
  if (hVar2!=NULL) hV.push_back(hVar2);
  if (hVar3!=NULL) hV.push_back(hVar3);
  if (hVar4!=NULL) hV.push_back(hVar4);

  // If includeVariants is 1, then differences between all histograms are
  // considered. The j index will run from 0 (baseValue) to 
  // hV.size(). If j>0, the reference is the variant histogram
  unsigned int jmax=(includeVariants) ? (hV.size()+1) : 1;

  for (unsigned int i=0; i<hV.size(); ++i) {
    for (unsigned int j=0; j<jmax; ++j) {
      if ((j>0) && (i+1==j)) continue; // avoid comparing distribution to itself
      TH2D* hd=(TH2D*)hV[i]->Clone(Form("var_%d",i));
      const TH2D* hRef= (j==0) ? baseValue : hV[j-1];
      hd->Add(hRef,-1.);
      for (int ibin=1; ibin<=hRes->GetNbinsX(); ++ibin) {
	for (int jbin=1; jbin<=hRes->GetNbinsY(); ++jbin) {
	  double err=hRes->GetBinError(ibin,jbin);
	  double var=fabs(hd->GetBinContent(ibin,jbin)); // !content
	  if (var>err) hRes->SetBinError(ibin,jbin, var);
	}
      }
      delete hd;
    }
  }
  return hRes;
}

//--------------------------------------------------
//--------------------------------------------------

TH2D* convertBaseH2actual(const TH2D* h2, TString newHistoName, int setTitle) {
  TH2D* newH=Clone(h2,newHistoName,newHistoName);
  if (!setTitle) newH->SetTitle("");
  for (int ibin=1; ibin<=h2->GetNbinsX(); ++ibin) {
    // adjustment is needed only if the bin count is different
    if (DYTools::nYBins[ibin-1]!=DYTools::nYBinsMax) {
      // ready only for double spacing
      if (2*DYTools::nYBins[ibin-1]!=DYTools::nYBinsMax) {
	std::cout << "convertBaseH2actual is ready only for double spacing\n";
	delete newH;
	return NULL;
      }
      for (int jbin=1; jbin<=h2->GetNbinsY(); jbin+=2) {
	double x=
	  h2->GetBinContent(ibin,jbin  ) +
	  h2->GetBinContent(ibin,jbin+1);
	double e2=
	  pow(h2->GetBinError(ibin,jbin  ),2) +
	  pow(h2->GetBinError(ibin,jbin+1),2);
	double e=sqrt(e2);
	int newJbin=(jbin+1)/2;
	newH->SetBinContent(ibin, newJbin, x);
	newH->SetBinError  (ibin, newJbin, e);
      }
      for (int jbin=DYTools::nYBins[ibin-1]+1; jbin<=DYTools::nYBinsMax; ++jbin) {
	newH->SetBinContent(ibin,jbin,  0.);
	newH->SetBinError  (ibin,jbin,  1.e9);
      }
    }
  }
  return newH;
}

//--------------------------------------------------

TH1D* createProfileX(TH2D *h2, int iyBin, const TString &name, int setTitle, const char *title) {
  if ((iyBin<=0) || (iyBin>h2->GetNbinsY())) {
    std::cout << "\n\n\tcreateProfileX(" << h2->GetName() << ", iyBin=" << iyBin << "(bad value!!), name=" << name << ")\n\n";
  }

  // prepare the range info
  int nxBins=h2->GetNbinsX();
  double *xv=new double[nxBins+1];
  TAxis *ax=h2->GetXaxis();
  for (int i=1; i<=nxBins; i++) {
    xv[i-1] = ax->GetBinLowEdge(i);
  }
  xv[nxBins]=ax->GetBinLowEdge(nxBins) + ax->GetBinWidth(nxBins);
  TH1D *h=new TH1D(name,"",h2->GetNbinsX(),xv);
  delete [] xv;

  // copy the profile
  h->SetDirectory(0);
  h->SetStats(0);
  if (setTitle) {
    if (title) h->SetTitle(title);
    else h->SetTitle(name);
  }
  h->GetXaxis()->SetTitle( h2->GetXaxis()->GetTitle() );
  for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
    h->SetBinContent(ibin,h2->GetBinContent(ibin,iyBin));
    h->SetBinError(ibin,h2->GetBinError(ibin,iyBin));
  }
  return h;
}
// -------------------------------------------
/*
inline
TH1D* createProfileY(TH2D *h2, int ixBin, const TString &name, int setTitle=0, const char *title=NULL) {
  if ((ixBin<=0) || (ixBin>h2->GetNbinsX())) {
    std::cout << "\n\n\tcreateProfileY(" << h2->GetName() << ", ixBin=" << ixBin << "(bad value!!), name=" << name << ")\n\n";
  }
  TH1D *h=new TH1D(name,"",h2->GetNbinsY(),h2->GetYaxis()->GetXbins()->GetArray());
  h->SetDirectory(0);
  h->SetStats(0);
  if (setTitle) {
    if (title) h->SetTitle(title);
    else h->SetTitle(name);
  }
  h->GetXaxis()->SetTitle( h2->GetYaxis()->GetTitle() );
  for (int ibin=1; ibin<=h2->GetNbinsY(); ibin++) {
    h->SetBinContent(ibin,h2->GetBinContent(ixBin,ibin));
    h->SetBinError(ibin,h2->GetBinError(ixBin,ibin));
  }
  return h;
}
*/

// -------------------------------------------

TH1D* createProfileY(TH2D *h2, int ixBin, const TString &name, int setTitle, const char *title, int set_nYbins, double set_ymin, double set_ymax) {
  if ((ixBin<=0) || (ixBin>h2->GetNbinsX())) {
    std::cout << "\n\n\tcreateProfileY(" << h2->GetName() << ", ixBin=" << ixBin << "(bad value!!), name=" << name << ")\n\n";
  }
  TH1D *h=NULL;
  if (set_nYbins==-1) {
    int nyBins=h2->GetNbinsY();
    double *yv=new double[nyBins+1];
    TAxis *ay=h2->GetYaxis();
    for (int i=1; i<=nyBins; i++) {
      yv[i-1] = ay->GetBinLowEdge(i);
      //std::cout << "yv[" << i-1 << "]=" << yv[i-1] << "\n";
    }
    yv[nyBins]=ay->GetBinLowEdge(nyBins) + ay->GetBinWidth(nyBins);
    std::cout << "yv[nYBins=" << nyBins << "]=" << yv[nyBins] << "\n";
    h=new TH1D(name,"",h2->GetNbinsY(),yv);
    delete [] yv;
    set_nYbins=h2->GetNbinsY();
  }
  else h=new TH1D(name,"",set_nYbins,set_ymin,set_ymax);

  h->SetDirectory(0);
  h->SetStats(0);
  if (setTitle) {
    if (title) h->SetTitle(title);
    else h->SetTitle(name);
  }
  h->GetXaxis()->SetTitle( h2->GetYaxis()->GetTitle() );
  for (int ibin=1; (ibin<=h2->GetNbinsY()) && (ibin<=set_nYbins); ibin++) {
    h->SetBinContent(ibin,h2->GetBinContent(ixBin,ibin));
    h->SetBinError(ibin,h2->GetBinError(ixBin,ibin));
  }
  return h;
}

//--------------------------------------------------
//--------------------------------------------------

TH1D* removeLastBin(const TH1D* hOrig, TString newName, int setTitle, const char *newTitle) {
  // prepare the range info
  int nxBins=hOrig->GetNbinsX()-1; // we will remove last bin
  double *xv=new double[nxBins+1]; 
  TAxis *ax=hOrig->GetXaxis();
  for (int i=1; i<=nxBins; i++) {
    xv[i-1] = ax->GetBinLowEdge(i);
  }
  xv[nxBins]=ax->GetBinLowEdge(nxBins) + ax->GetBinWidth(nxBins); // last bin  right edge

  TH1D *h=new TH1D(newName,newName,nxBins,xv);
  delete [] xv;

  // copy the profile
  h->SetDirectory(0);
  h->SetStats(0);
  if (setTitle) {
    if (newTitle) h->SetTitle(newTitle);
    else h->SetTitle(newName);
  }
  h->GetXaxis()->SetTitle( hOrig->GetXaxis()->GetTitle() );
  for (int ibin=1; ibin<=nxBins; ibin++) {
    h->SetBinContent(ibin,hOrig->GetBinContent(ibin));
    h->SetBinError(ibin,hOrig->GetBinError(ibin));
  }
  return h;
}

//--------------------------------------------------
//--------------------------------------------------

TH2D* LoadHisto2D(TString histoName, const TString &fname, TString subDir, int checkBinning) {
  TString theCall=TString("LoadHisto2D(<") + histoName + TString(">,<") + fname + TString(">,<") + subDir + TString(Form(">, checkBinning=%d)",checkBinning));

  TFile fin(fname,"read");
  if (!fin.IsOpen()) {
    std::cout << theCall << ": failed to open the file\n";
    return NULL;
  }
  if (checkBinning && !checkBinningArrays(fin)) {
    std::cout << theCall << ": binning test failed\n";
    return NULL;
  }

  TH2D* h2=LoadHisto2D(fin,histoName,subDir,0);
  if (!h2) std::cout << theCall << ": failed to load the histo\n";
  else {
    h2->SetDirectory(0);
  }
  fin.Close();
  return h2;
}

//--------------------------------------------------

TH2D* LoadHisto2D(TFile &fin, TString histoName, TString subDir, int checkBinning) {
  TString theCall=TString("LoadHisto2D(file=<") + fin.GetName() + TString(">,<") + histoName + TString(">,<") + subDir + TString(Form(">, checkBinning=%d)",checkBinning));
  
  if (checkBinning && !checkBinningArrays(fin)) {
    std::cout << theCall << ": binning test failed\n";
    return NULL;
  }

  TString loadHistoName;
  if (subDir.Length()) {
    if (subDir[subDir.Length()-1]!='/') subDir.Append("/");
    loadHistoName=subDir + histoName; 
  }
  else loadHistoName=histoName;
  //std::cout << "loadHistoName=<" << loadHistoName << ">\n";

  TH2D* h2=(TH2D*)fin.Get(loadHistoName);
  if (!h2) std::cout << theCall << ": failed to load the histo\n";
  else {
    h2->SetDirectory(0);
  }
  return h2;
}

//--------------------------------------------------
//--------------------------------------------------

void writeBinningArrays(TFile &fout) {
  fout.cd();
  TVectorD mass(DYTools::nMassBins+1);
  TVectorD rapidityCounts(DYTools::nMassBins);
  for (int i=0; i<DYTools::nMassBins+1; i++) mass[i]=DYTools::massBinLimits[i];
  for (int i=0; i<DYTools::nMassBins  ; i++) rapidityCounts[i]=DYTools::nYBins[i];
  mass.Write("massBinLimits"); // was 'massBinning'
  rapidityCounts.Write("rapidityCounts");
}

//--------------------------------------------------

int checkBinningArrays(TFile &fin) {
  const char *fncname="unfolding::checkBinningArrays: ";
  TString fileInfo=TString("on file <") + fin.GetName() + ">";
  fin.cd();
  TVectorD* mass= (TVectorD*)fin.FindObjectAny("massBinLimits");
  TVectorD* rapidityCounts= (TVectorD*)fin.FindObjectAny("rapidityCounts");
  if (!mass || !rapidityCounts) {
    std::cout << fncname << "file <" << fin.GetName() << "> does not contain keys massBinning and/or rapidityCounts\n";
    return false;
  }
  int comparisonOk=checkBinningRanges(*mass,*rapidityCounts, fin.GetName());
  if (mass) delete mass;
  if (rapidityCounts) delete rapidityCounts;
  return comparisonOk;
}

//--------------------------------------------------

int checkBinningRanges(const TVectorD &mass, const TVectorD &rapidityCounts, const TString &fname) {
  const char *fncname="unfolding::checkBinningRanges: ";
  TString fileInfo=TString("on file <") + fname + ">";
  
  bool massOk=true, rapidityOk=true;
  int first_err=1;
  if (mass.GetNoElements() != DYTools::nMassBins+1) {
    std::cout << errdash; first_err=0;
    std::cout << "\n" << fncname 
	      << " number of mass bins " << fileInfo
	      << " is " << mass.GetNoElements() 
	      << " while " << (DYTools::nMassBins+1) << " is expected\n";
    massOk=false;
  }
  
  if (rapidityCounts.GetNoElements() != DYTools::nMassBins ) {
    if (first_err) { first_err=0; std::cout << errdash; }
    std::cout << "\n" << fncname
	      << "number of mass bins in rapidityCounts " << fileInfo
		<< " is " << rapidityCounts.GetNoElements() 
	      << " while " << DYTools::nMassBins << " is expected\n";
    rapidityOk=false;
  }
  
  if (massOk) {
    for(int i=0; i<DYTools::nMassBins+1; i++){
      if( DYTools::massBinLimits[i] != mass[i] ) {
	if (first_err) { first_err=0; std::cout << errdash; }
	std::cout << fncname
		  << " mass limit " << fileInfo
		  << " at i=" << i 
		  << " is " << mass[i] 
		  << " instead of expected " << DYTools::massBinLimits[i] << "\n";
	massOk=false;
      }
    }
  }
  if (rapidityOk) {
    for (int i=0; i<DYTools::nMassBins; i++) {
      if ( DYTools::nYBins[i] != rapidityCounts[i] ) {
	if (first_err) { first_err=0; std::cout << errdash; }
	std::cout << fncname 
		  << "y bin count " << fileInfo
		  << " at i=" << i 
		  << " is " << rapidityCounts[i] 
		  << " instead of expected " << DYTools::nYBins[i] << "\n";
	rapidityOk=false;
      }
    }
  }
  if (!massOk || !rapidityOk) {
    if (first_err) { first_err=0; std::cout << errdash; }
    std::cout << "file info: mass[" << mass.GetNoElements() << "]: ";
    for (int i=0; i<mass.GetNoElements(); ++i) {
      std::cout << " " << mass[i];
    }
    std::cout << "\n";
    std::cout << "file info: rapidityCounts[" 
	      << rapidityCounts.GetNoElements() << "]: ";
    for (int i=0; i<rapidityCounts.GetNoElements(); ++i) {
      std::cout << " " << rapidityCounts[i];
    }
    std::cout << "\n" << std::endl;
  }
  if (!first_err) std::cout << errdash;
  return (massOk && rapidityOk) ? 1:0;
}

//--------------------------------------------------

int checkMatrixSize(const TMatrixD &m, const TString &name) {
  int res= ((m.GetNrows() == DYTools::nMassBins) &&
	    (m.GetNcols() == DYTools::nYBinsMax)) ? 1:0;
  if (!res) {
    std::cout << "matrix (name=" << name << ") failed range check\n";
    std::cout << "expected: " << DYTools::nMassBins << "x" << DYTools::nYBinsMax << ", got " << m.GetNrows() << "x" << m.GetNcols() << std::endl;
  }
  return res;
}

//--------------------------------------------------
//--------------------------------------------------

void printMatrix(const TString &name, const TMatrixD &M, int exponent) {
  std::cout << "Matrix " << name << "\n";
  if (!exponent) M.Print();
  else {
    const char *format="%8.4e";
    for (int ir=0; ir<M.GetNrows(); ++ir) {
      for (int ic=0; ic<M.GetNcols(); ++ic) {
	std::cout << " " << Form(format,M(ir,ic));
      }
      std::cout << "\n";
    }
  }
}

//--------------------------------------------------

TMatrixD* loadMatrix(const TString &fname, const TString &fieldName, int expect_nRows, int expect_nCols, int reportFieldError) {
  TFile f(fname,"read");
  TMatrixD *M=NULL;
  int ok=1;
  if (!f.IsOpen()) ok=0;
  if (ok==1) {
    M=(TMatrixD*)f.Get(fieldName);
    f.Close();
    if (!M) ok=-1;
    else {
      if ((M->GetNrows()!=expect_nRows) ||
	  (M->GetNcols()!=expect_nCols)) {
	ok=-2;
      }
    }
  }
  if (ok!=1) {
    int report=1;
    if ((ok==-1) && !reportFieldError) report=0;
    if (report) {
      std::cout << "Error in loadMatrix(fname=<" << fname << ">, fieldName=<" << fieldName << ">, nRows=" << expect_nRows << ", nCols=" << expect_nCols << "):\n";
      if (ok==0) std::cout << " - failed to open the file\n";
      else if (ok==-1) std::cout << " - failed to load the field\n";
      else if (ok==-2) {
	std::cout << " - size mistmatch. Expect " << expect_nRows << "x" << expect_nCols << ", got " << M->GetNrows() << "x" << M->GetNcols() << "\n";
	delete M;
	M=NULL;
      }
    }
  }
  return M;
}

//--------------------------------------------------

TH2D* LoadMatrixFields(TFile &fin, const TString &field, const TString &fieldErr, int loadErr, int absoluteRapidity) {
  TMatrixD* val= NULL;
  if (field.Length()) val=(TMatrixD*)fin.Get(field);
  TMatrixD* err=(loadErr) ? (TMatrixD*)fin.Get(fieldErr) : NULL;
  if (   (field.Length() && !val) 
	 || (loadErr && !err)
	 || (val && !checkMatrixSize(*val,field))
	 || (err && !checkMatrixSize(*err,fieldErr))
      //|| ((loadErr) && err && !checkMatrixSize(*err,fieldErr))
      ) {
    std::cout << "error in LoadMatrixFields(fname=" << fin.GetName() << ",""" << field << """, """ << fieldErr << """, loadErr=" << loadErr << "):\n";
    if (!val) std::cout << "\t failed to get values\n";
    if ((loadErr==1) && !err) std::cout << "\t failed to get errors\n";
    return NULL;
  }
  //std::cout << "loaded field=<" << field << "> from file <" << fin.GetName() << ">\n";
  //std::cout << "values are "; val->Print();
  TH2D *h2=createBaseH2(field,field,absoluteRapidity);
  for (int ibin=1; ibin<=h2->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=h2->GetNbinsY(); ++jbin) {
      double value=(val) ? (*val)(ibin-1,jbin-1) : 0.;
      double e=(loadErr==1) ? (*err)(ibin-1,jbin-1) : 0.;
      if (loadErr==2) e= ( (*err)(ibin-1,jbin-1) * (*err)(ibin-1,jbin-1) );
      if (loadErr==3) { value=0; e=(*err)(ibin-1,jbin-1); }
      //std::cout << "putting ibin=" << ibin << ", jbin=" << jbin << ", value=" << value << "\n";
      h2->SetBinContent(ibin,jbin, value);
      h2->SetBinError(ibin,jbin, e);
    }
  }
  return h2;
}

//--------------------------------------------------

TH2D* LoadMatrixFields(const TString &fname, int checkBinning, const TString &field, const TString &fieldErr, int loadErr, int absoluteRapidity) {
  TFile fin(fname,"read");
  if (!fin.IsOpen()) {
    std::cout << "LoadMatrixField(fname): failed to open file <" << fname << ">\n";
    return NULL;
  }
  int res=1;
  TH2D *h2=NULL;
  if (res && checkBinning) res=checkBinningArrays(fin);
  if (res) h2=LoadMatrixFields(fin,field,fieldErr,loadErr,absoluteRapidity);
  fin.Close();
  if (!res || !h2) std::cout << "error in LoadMatrixField(fname)\n";
  return h2;
}

//--------------------------------------------------

int LoadThreeMatrices(TFile &fin, TH2D **h2, TH2D **h2syst, const TString &field, const TString &fieldErr, const TString &fieldSystErr, int absoluteRapidity) {
  int res=1;
  (*h2)=LoadMatrixFields(fin,field,fieldErr,1,absoluteRapidity);
  if (!(*h2)) res=0;
  if (res) (*h2syst)=LoadMatrixFields(fin,"",fieldSystErr,3,absoluteRapidity);
  if (!(*h2syst)) res=0;
  if (!res) {
    std::cout << "error in LoadThreeMatrices(fname=" << fin.GetName() << ", """ << field << """, """ << fieldErr << """, """ << fieldSystErr << """)\n";
  }
  return res;
}

//--------------------------------------------------

int LoadThreeMatrices(const TString &fileName, TH2D **h2, TH2D **h2syst, const TString &field, const TString &fieldErr, const TString &fieldSystErr, int checkBinning, int absoluteRapidity) {
  TFile file(fileName,"read");
  if (!file.IsOpen()) { 
    std::cout << "failed to open <" << fileName << ">\n"; 
    return 0;
  }
  int res=1;
  if (res && checkBinning) res=checkBinningArrays(file);
  if (res) {
    res=LoadThreeMatrices(file,h2,h2syst,
			  field,fieldErr,fieldSystErr,
			  absoluteRapidity);
  }
  file.Close();
  if (!res) {
    std::cout << "failed to load data from <" << fileName << ">\n";
  }
  return res;
}

//--------------------------------------------------
//--------------------------------------------------


inline
TH2D* extractSubArea(TH2D *histo, int xbin1, int xbin2, int ybin1, int ybin2, const TString &newName, int setTitle, int resetAxis) {
  int nXBins=histo->GetNbinsX();
  int nYBins=histo->GetNbinsY();
  if ((xbin1==0) || (ybin1==0) ||
      (xbin1>xbin2) || (ybin1>ybin2)) {
    std::cout << "extractSubArea: problem with the requested area ";
    printf("(xbin1=%d, xbin2=%d, ybin1=%d, ybin2=%d)\n",xbin1,xbin2,ybin1,ybin2);
    std::cout << std::endl;
    return NULL;
  }
  if ((xbin2>nXBins) || (ybin2>nYBins)) {
    std::cout << "extractSubArea: \n";
    printf("requested last bin x=%d (available %d), y=%d (%d)",xbin2,nXBins,ybin2,nYBins);
    std::cout << std::endl;
    return NULL;
  }
  
  const TAxis *ax=histo->GetXaxis();
  const TAxis *ay=histo->GetYaxis();
  
  double *xnew=new double[xbin2-xbin1+2];
  double *ynew=new double[ybin2-ybin1+2];
  if (resetAxis) {
    // ignore the actual shift of the axes
    int lastXbin=xbin2-xbin1+1;
    int lastYbin=ybin2-ybin1+1;
    for (int i=0; i<lastXbin; ++i) xnew[i]=ax->GetBinLowEdge(i+1);
    xnew[lastXbin]= ax->GetBinLowEdge(lastXbin) + ax->GetBinWidth(lastXbin);
    for (int i=0; i<lastYbin; ++i) ynew[i]=ay->GetBinLowEdge(i+1);
    ynew[lastYbin]= ay->GetBinLowEdge(lastYbin) + ay->GetBinWidth(lastYbin);
  }
  else { 
    // keep original labels of the axis
    for (int i=0; i<xbin2-xbin1+1; ++i) xnew[i]= ax->GetBinLowEdge(i+xbin1);
    xnew[xbin2-xbin1+1]= ax->GetBinLowEdge(xbin2) + ax->GetBinWidth(xbin2);
    for (int i=0; i<ybin2-ybin1+1; ++i) ynew[i]= ay->GetBinLowEdge(i+ybin1);
    ynew[ybin2-ybin1+1]= ay->GetBinLowEdge(ybin2) + ay->GetBinWidth(ybin2);
  }

  TH2D *h2=new TH2D(newName,"",xbin2-xbin1+1,xnew,ybin2-ybin1+1,ynew);
  if (setTitle) h2->SetTitle(newName);
  h2->SetDirectory(0);
  h2->SetStats(0);
  h2->GetXaxis()->SetTitle(histo->GetXaxis()->GetTitle());
  h2->GetYaxis()->SetTitle(histo->GetYaxis()->GetTitle());

  for (int ibin=xbin1; ibin<=xbin2; ++ibin) {
    for (int jbin=ybin1; jbin<=ybin2; ++jbin) {
      const int new_ibin= ibin-xbin1+1;
      const int new_jbin= jbin-ybin1+1;
      h2->SetBinContent(new_ibin,new_jbin,histo->GetBinContent(ibin,jbin));
      h2->SetBinError  (new_ibin,new_jbin,histo->GetBinError  (ibin,jbin));
    }
  }
  delete [] xnew;
  delete [] ynew;
  return h2;
}


//--------------------------------------------------
//--------------------------------------------------

