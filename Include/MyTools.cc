#include "../Include/MyTools.hh"

//--------------------------------------------------
//--------------------------------------------------

TH2D* LoadHisto2D(const TString &histoName, const TString &fname, const TString &subDir, int checkBinning) {
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

  //HERE("LoadHisto2D get histo %s",histoName.Data());
  
  TH2D* h2=(TH2D*)fin.Get(histoName);
  h2->SetDirectory(0);
  fin.Close();
  if (!h2) std::cout << theCall << ": failed to load the histo\n";
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
  TH2D *h2=createBaseH2(field,field,absoluteRapidity);
  for (int ibin=1; ibin<=h2->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=h2->GetNbinsY(); ++jbin) {
      double value=(val) ? (*val)(ibin-1,jbin-1) : 0.;
      double e=(loadErr==1) ? (*err)(ibin-1,jbin-1) : 0.;
      if (loadErr==2) e= ( (*err)(ibin-1,jbin-1) * (*err)(ibin-1,jbin-1) );
      if (loadErr==3) { value=0; e=(*err)(ibin-1,jbin-1); }
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

