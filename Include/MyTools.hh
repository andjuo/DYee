#ifndef MYTOOLS_HH
#define MYTOOLS_HH

#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TVector.h>
#include <TDirectory.h>
//#include <TKey.h>
#include <iostream>
#include <math.h>
#include <stdlib.h> // atoi, atof
#include "../Include/CPlot.hh"
#include "../Include/DYTools.hh"
#include "../Include/colorPalettes.hh"
#include "../Include/ComparisonPlot.hh"

//----------------------------------------

const std::string dashline=std::string(65,'-') + std::string("\n");
const std::string mainpart=std::string("\n") + dashline + std::string("\t\t\t--- Main Part ---\n") + dashline + std::string("\n");
const std::string errdash=std::string(30,'=') + std::string("\n\t\tError\n") + std::string(30,'=') + std::string("\n");

//----------------------------------------
/*
namespace toolbox {

inline
Double_t deltaR(Double_t eta1, Double_t phi1, Double_t eta2, Double_t phi2)
{
  const Double_t pi = 3.14159265358979;
  Double_t dphi = fabs(phi1-phi2);
  while (dphi>pi)
    dphi = fabs(dphi - 2.0*pi);

  Double_t deta = eta1-eta2;

  return sqrt(dphi*dphi + deta*deta);
}

}
//------------------------------------------------------------------------------------------------------------------------
namespace toolbox {

inline
Double_t deltaPhi(Double_t phi1, Double_t phi2){
  // Compute dPhi between two given angles. Results is in [0,pi].
  const Double_t pi = 3.14159265358979;
  Double_t dphi = fabs(phi1-phi2);
  while (dphi>pi)
    dphi = fabs(dphi - 2.0*pi);

  return dphi;
}

}
//------------------------------------------------------------------------------------------------------------------------
namespace toolbox {

inline
Int_t roundToInt(Double_t x)
{
  if(x>0)
    return ((x-floor(x)) < (ceil(x)-x)) ? (Int_t)floor(x) : (Int_t)ceil(x);
  else
    return ((x-floor(x)) < (ceil(x)-x)) ? (Int_t)ceil(x) : (Int_t)floor(x);
}

}
*/


//------------------------------------------------------------------------------------------------------------------------
// error on ratio when counts are pass and fail (i.e. independent)

inline
double errOnRatio(double a, double da, double b, double db){

  const double zero=static_cast<double>(0);
  double result = zero;
  if(a == zero || b == zero)
    return result;
  result = (a/b)*sqrt( (da/a)*(da/a) + (db/b)*(db/b) );
  return result;
}


//----------------------------------------------------------------
// error on ratio when counts are pass and total

inline
double errOnRatioPT(double p, double dp, double tot, double dtot){
  double f=tot-p;
  double df=sqrt(dtot*dtot-dp*dp);
  return errOnRatio(p,dp,f,df);
}


//----------------------------------------------------------------
// print error on ratio when counts are pass and total

inline
void printEff_PT(double p, double dp, double tot, double dtot) {
  double eff=p/tot;
  double deff=errOnRatioPT(p,dp,tot,dtot);
  std::cout << "pass=" << p << " +- " << dp << ", tot="
	    << tot << " +- " << dtot << ", eff="
	    << eff << " +- " << deff << "\n";
}

//----------------------------------------------------------------
// print error on ratio when counts are pass and total

inline
void printError_ab(double a, double da, double b, double db) {
  double r=a*b;
  double relErr=sqrt( (da/a)*(da/a) + (db/b)*(db/b) );
  std::cout << "(" << a << " +- " << da << ") x ("
	    << b << " +- " << db << ") = ("
	    << r << " +- " << r*relErr << "); RelErr=" << relErr << "\n";
}

//------------------------------------------------------------------------------------------------------------------------

template<class T>
inline
void PrintVec(const char *msg, const std::vector<T>& vec, int prneol=0) {
  if (msg) std::cout << msg;
  std::cout << "vec[" << vec.size() << "]: ";
  for (unsigned int i=0; i<vec.size(); ++i) {
    if (prneol) std::cout << "\n" << i << ") ";
    std::cout << " " << vec[i];
  }
  if (prneol) std::cout << "\n";
}

//------------------------------------------------------------------------------------------------------------------------

void printMeanAndRMS(TH1D *histo, const char *extra_msg=NULL);

int PrintHisto2Dvec(const char *msg, const std::vector<TH2D*> &vec, int exponent=0, int maxLines=-1);
int PrintTwoHistos(const char *msg, TH2D *h2a, TH2D *h2b, int exponent=0,
		   int maxLines=-1);

//------------------------------------------------------------------------------------------------------------------------

inline
void SaveCanvas(TCanvas* canv, const TString &canvName, TString destDir=CPlot::sOutDir){
  gSystem->mkdir(destDir,kTRUE);
  gSystem->mkdir(destDir+TString("/png"),kTRUE);
  gSystem->mkdir(destDir+TString("/pdf"),kTRUE);
  gSystem->mkdir(destDir+TString("/root"),kTRUE);

  TString saveName=destDir+TString("/png/");
  saveName+=canvName;
  saveName+=".png";
  canv->SaveAs(saveName);
  saveName.ReplaceAll("png","pdf");
  canv->SaveAs(saveName);
  saveName.ReplaceAll("pdf","root");
  canv->SaveAs(saveName);
  return;
}
// ----------------------------------------------------------

template<class T>
inline
void ClearVec(std::vector<T*> &vec) {
  for (unsigned int i=0; i<vec.size(); ++i) if (vec[i]) delete vec[i];
  vec.clear();
}

//-----------------------------------------------------------

std::vector<TString>* createMassRangeVec(TString prependStr="");
int replaceAll(std::vector<TString*> &vec, TString oldText, TString newText);
int replaceAll(std::vector<TString> &vec, TString oldText, TString newText);
TString niceNumber(int iVal, int iValMax);

//------------------------------------------------------------------------------------------------------------------------

#ifndef BaseClass_HH

inline
void HERE(const char *msg) {
  std::cout << ((msg) ? msg : "HERE") << std::endl;
}

//--------------------------------------------------------------

inline
void HERE(const std::string &msg) {
  std::cout << msg << std::endl;
}

//--------------------------------------------------------------

template<class Type_t>
inline
void HERE(const char *format, const Type_t &a) {
  std::cout << Form(format,a) << std::endl;
}

//--------------------------------------------------------------

template<class Type1_t, class Type2_t>
inline
void HERE(const char *format, const Type1_t &a, const Type2_t &b) {
  std::cout << Form(format,a,b) << std::endl;
}

//--------------------------------------------------------------

template<>
inline
void HERE(const char *format, const TString &a) {
  const char *aStr=(a.Length()) ? a.Data() : " ";
  std::cout << Form(format,aStr) << std::endl;
}
//--------------------------------------------------------------

template<>
inline
void HERE(const char *format, const TString &a, const TString &b) {
  const char *aStr=(a.Length()) ? a.Data() : " ";
  const char *bStr=(b.Length()) ? b.Data() : " ";
  std::cout << Form(format,aStr,bStr) << std::endl;
}

//--------------------------------------------------------------

template<class Type_t>
inline
void HERE(const char *format, const TString &a, const TString &b, const Type_t &x ) {
  const char *aStr=(a.Length()) ? a.Data() : " ";
  const char *bStr=(b.Length()) ? b.Data() : " ";
  std::cout << Form(format,aStr,bStr,x) << std::endl;
}
#endif

//------------------------------------------------------------------------------------------------------------------------

inline bool PosOk(size_t pos) { return (pos==std::string::npos) ? 0 : 1; }

//----------------------------------------------------------------------

template<class T>
inline bool PosOk(const std::string &s, const T& substr) {
  return (s.find(substr)==std::string::npos) ? 0 : 1;
}

//------------------------------------------------------------------------------------------------------------------------

inline
int printHisto(std::ostream& out, const TH1D* histo, int exponent=0, int maxLines=-1) {
  if (!histo) {
    out << "printHisto: histo is null\n";
    return 0;
  }
  char buf[100];
  const char *format= (exponent) ?
    " %5.2f-%5.2f    %e    %e\n" :
    " %5.2f-%5.2f    %f    %f\n";

  out << "values of " << histo->GetName() << "\n";
  int imax=histo->GetNbinsX();
  int truncated=0;
  if ((maxLines>0) && (imax>maxLines)) { imax=maxLines; truncated=1; }
  for(int i=1; i<=imax; i++) {
    double x=histo->GetBinLowEdge(i);
    double w=histo->GetBinWidth(i);
    sprintf(buf,format,
	    x,x+w,histo->GetBinContent(i),histo->GetBinError(i));
    out << buf;
  }
  if (truncated) out << "... output truncated to " << maxLines << " lines\n";
  return 1;
}

//------------------------------------------------------------------------------------------------------------------------

inline int printHisto(const TH1D* histo, int exponent=0, int maxLines=-1) { return printHisto(std::cout, histo, exponent, maxLines); }

int printHistoAssumeMassBinning(const TH1D* histo,
				int exponent=0, int maxLines=-1);

//------------------------------------------------------------------------------------------------------------------------

inline
int printHistoErr(std::ostream& out, const TH1D* histo, const TH1D* histoSystErr, int exponent=0, int maxLines=-1) {
  if (!histo || !histoSystErr) {
    out << "printHisto: histo is null\n";
    return 0;
  }
  char buf[100];
  const char *format= (exponent) ?
    " %5.2f-%5.2f    %e    %e    %e\n" :
    " %5.2f-%5.2f    %f    %f    %f\n";

  out << "values of " << histo->GetName() << "\n";
  int imax=histo->GetNbinsX();
  int truncated=0;
  if ((maxLines>0) && (imax>maxLines)) { imax=maxLines; truncated=1; }
  for(int i=1; i<=imax; i++) {
    double x=histo->GetBinLowEdge(i);
    double w=histo->GetBinWidth(i);
    sprintf(buf,format,
	    x,x+w,histo->GetBinContent(i),histo->GetBinError(i),histoSystErr->GetBinError(i));
    out << buf;
  }
  if (truncated) out << "... output truncated to " << maxLines << " lines\n";
  return 1;
}

//------------------------------------------------------------------------------------------------------------------------

inline int printHistoErr(const TH1D* histo, const TH1D* histoSystErr, int exponent=0, int maxLines=-1) { return printHistoErr(std::cout, histo, histoSystErr, exponent,maxLines); }

//------------------------------------------------------------------------------------------------------------------------

inline
int printHistoErr(std::ostream& out, const TH2D* histo, const TH2D* histoSystErr, int exponent=0, int maxLinesX=-1, int maxLinesY=-1) {
  if (!histo) {
    out << "printHistoErr: histo is null\n";
    return 0;
  }
  char buf[200];
  const char *format_char= (exponent) ?
    " %5.2f-%5.2f  %5.2f-%5.2f    %e    %e    %e\n" :
    " %5.2f-%5.2f  %5.2f-%5.2f    %f    %f    %f\n";
  std::string format=format_char;
  if (!histoSystErr) {
    size_t p1=format.find_last_of(' ');
    size_t p2=format.find_last_not_of(' ',p1);
    format.erase(p2+1,format.size()-p2-2);
  }

  out << "values of " << histo->GetName() << "\n";
  int imax=histo->GetNbinsX();
  int jmax=histo->GetNbinsY();
  int truncated=0;
  if ((maxLinesX>0) && (imax>maxLinesX)) { imax=maxLinesX; truncated=1; }
  if ((maxLinesY>0) && (jmax>maxLinesY)) { jmax=maxLinesY; truncated=1; }
  for(int i=1; i<=imax; i++) {
    double x=histo->GetBinLowEdge(i);
    double w=histo->GetBinWidth(i);
    for (int j=1; j<=jmax; ++j) {
      double y=histo->GetYaxis()->GetBinLowEdge(j);
      double h=histo->GetYaxis()->GetBinWidth(j);
      double err2=(histoSystErr) ? histoSystErr->GetBinError(i,j) : 0.0;
      sprintf(buf,format.c_str(),
	      x,x+w,y,y+h,
	      histo->GetBinContent(i,j),
	      histo->GetBinError(i,j),
	      err2);
      out << buf;
    }
    if (jmax!=histo->GetNbinsY()) out << "  ...\n";
  }
  if (truncated) out << "... output truncated to (" << maxLinesX
		     << ", " << maxLinesY << ") lines\n";
  return 1;
}

//------------------------------------------------------------------------------------------------------------------------

inline int printHistoErr(const TH2D* histo, const TH2D* histoSystErr, int exponent=0, int maxLinesX=-1, int maxLinesY=-1) { return printHistoErr(std::cout, histo, histoSystErr, exponent, maxLinesX, maxLinesY); }

inline int printHisto(const TH2D* histo, int exponent=0, int maxLines=-1) { return printHistoErr(std::cout, histo, NULL, exponent, maxLines); }

void printHisto(const std::vector<TH2D*> hV, int exponent=0, int maxLines=-1, int maxEntries=-1);

//------------------------------------------------------------------------------------------------------------------------


inline
TH1D *extractRapidityDependence(const TString &name, const TString &title,
				const TMatrixD &m, const TMatrixD &mErr,
				int iMassBin, int perMassBinWidth=0) {
  TString hName= name + Form("_massBin%d",iMassBin);
  TString hTitle= title + Form("_massBin%d",iMassBin);
  TH1D *h=new TH1D(name,title,DYTools::nYBins[iMassBin],DYTools::yRangeMin,DYTools::yRangeMax);
  h->SetDirectory(0);
  h->Sumw2();
  for (int iY=0; iY<DYTools::nYBins[iMassBin]; ++iY) {
    double factor= (perMassBinWidth==0) ?
      1 : 1/(DYTools::massBinLimits[iMassBin+1] - DYTools::massBinLimits[iMassBin]);
    h->SetBinContent(iY+1, m[iMassBin][iY]*factor);
    h->SetBinError(iY+1, fabs(mErr[iMassBin][iY])*factor);
  }
  return h;
}

// -----------------------------------------------------------------------------

inline
TH1D *extractMassDependence(const TString &name, const TString &title,
			    const TMatrixD &m, const TMatrixD &mErr,
			    int iYBin,
			    int perMassBinWidth=1, int perRapidityBinWidth=0) {
  if (perRapidityBinWidth) std::cout << "\n\tWARNING: extractMassDependence: perRapidityBinWidth=1 is experimental\n\n";
  std::cout << "dims  m[" << m.GetNrows() << ',' << m.GetNcols() << "], "
	    << " mErr[" << mErr.GetNrows() << ',' << mErr.GetNcols() << "]"
	    << std::endl;
  if ((m.GetNrows()!=DYTools::nMassBins) ||
      (m.GetNcols()!=DYTools::nYBinsMax) ||
      (mErr.GetNrows()!=DYTools::nMassBins) ||
      (mErr.GetNcols()!=DYTools::nYBinsMax)) {
    std::cout << "extractMassDependence: expecting matrices "
	      << DYTools::nMassBins << "x" << DYTools::nYBinsMax
	      << "  instead of  "
	      << m.GetNrows() << "x" << m.GetNcols() << "\n";
    assert(0);
  }
  TString hName= name + Form("_yBin%d",iYBin);
  TString hTitle= title + Form("_yBin%d",iYBin);
  TH1D *h=new TH1D(name,title,DYTools::nMassBins,DYTools::massBinLimits);
  h->SetDirectory(0);
  h->Sumw2();
  for (int iM=0; iM<DYTools::nMassBins; ++iM) {
    bool lastBin= (iM==DYTools::nMassBins-1) ? true : false;
    if (lastBin && (iYBin!=0)) {
      double yc=DYTools::findAbsYValue(0,iYBin);
      iYBin=DYTools::findAbsYBin(iM,yc);
    }
    double val=m[iM][iYBin];
    double factor=1.;
    if ((iM==DYTools::nMassBins-1) && DYTools::study2D) {
      // this correction may be needed for the 2D case
     factor *= DYTools::nYBins[iM]/double(DYTools::nYBins[0]);
     std::cout << "extractMassDependence: change factor: 1 --> " << factor << "\n";
    }
    if (perRapidityBinWidth) {
      factor *= (DYTools::yRangeMax-DYTools::yRangeMin)/double(DYTools::nYBins[iM]);
    }
    if (perMassBinWidth) factor/=(DYTools::massBinLimits[iM+1]-DYTools::massBinLimits[iM]);
    //std::cout << "factor=" << factor << "\n";
    h->SetBinContent(iM+1, val*factor);
    h->SetBinError(iM+1, fabs(mErr[iM][iYBin])*factor);
  }
  return h;
}

//-----------------------------------------------------------------

/*
inline
TH1D *extractMassDependenceSpec(const TString &name, const TString &title,
			    const TMatrixD &m, const TMatrixD &mErr,
			    int iYBin,
			const TVectorD &massGrid, const TVectorD &yBinCount,
			    int perMassBinWidth=1, int perRapidityBinWidth=0) {
  if (perRapidityBinWidth) std::cout << "\n\tWARNING: extractMassDependenceSpec: perRapidityBinWidth=1 is experimental\n\n";
  std::cout << "dims  m[" << m.GetNrows() << ',' << m.GetNcols() << "], "
	    << " mErr[" << mErr.GetNrows() << ',' << mErr.GetNcols() << "]"
	    << ", massGrid[" << massGrid.GetNoElements() << "]"
	    << ", yBinCount[" << yBinCount.GetNoElements() << "]"
	    << std::endl;
  TVectorD yBinCountCheck(massGrid.GetNoElements()-1);
  yBinCountCheck=1;
  if (iYBin!=0) {
    std::cout << "extractMassDependenceSpec: iYBin should be 0\n";
    assert(0);
  }
  if ((m.GetNrows()!=massGrid.GetNoElements()-1) ||
      (m.GetNcols()!=1) ||
      (mErr.GetNrows()!=massGrid.GetNoElements()-1) ||
      (mErr.GetNcols()!=1)) {
    std::cout << "extractMassDependenceSpec: expecting matrices "
	      << massGrid.GetNoElements() << "x" << 1
	      << "  instead of  "
	      << m.GetNrows() << "x" << m.GetNcols() << "\n";
    assert(0);
  }
  TString hName= name + Form("_yBin%d",iYBin);
  TString hTitle= title + Form("_yBin%d",iYBin);
  double *massBins=new double[massGrid.GetNoElements()];
  for (int i=0; i<massGrid.GetNoElements(); ++i) {
    massBins[i]=massGrid[i];
    //std::cout << "massBins[i=" << i << "]=" << massBins[i] << "\n";
  }
  TH1D *h=new TH1D(name,title,massGrid.GetNoElements()-1,massBins);
  delete massBins;

  h->SetDirectory(0);
  h->Sumw2();
  for (int iM=0; iM<m.GetNrows(); ++iM) {
    double val=m[iM][iYBin];
    double factor=1.;
    if (perRapidityBinWidth) {
      factor *= (DYTools::yRangeMax-DYTools::yRangeMin)/double(yBinCount[iM]);
    }
    if (perMassBinWidth) factor/=(massGrid[iM+1]-massGrid[iM]);
    h->SetBinContent(iM+1, val*factor);
    h->SetBinError(iM+1, fabs(mErr[iM][iYBin])*factor);
  }

  return h;
}
*/

//-----------------------------------------------------------------
//-----------------------------------------------------------------

void printMatrix(const TString &name, const TMatrixD &M, int exponent=0);

//-----------------------------------------------------------------

inline
void printYields(const TString &name, const TMatrixD &cs, const TMatrixD &csErr, const TMatrixD &csErrSyst, int printSystError=1) {
  std::cout << "\nprintYields for name=<" << name << ">\n";
  if ((cs.GetNrows() != csErr.GetNrows()) ||
      (cs.GetNrows() != csErrSyst.GetNrows())) {
    printf(" -- numbers of rows is different: %d, %d, %d\n",cs.GetNrows(),csErr.GetNrows(),csErrSyst.GetNrows());
    assert(0);
  }
  if ((cs.GetNcols() != csErr.GetNcols()) ||
      (cs.GetNcols() != csErrSyst.GetNcols())) {
    printf(" -- numbers of rows is different: %d, %d, %d\n",cs.GetNcols(),csErr.GetNcols(),csErrSyst.GetNcols());
    assert(0);
  }
  std::cout << "cs.GetNcols()=" << cs.GetNcols() << ", cs.GetNrows="<< cs.GetNrows() << "\n";
  if (printSystError) {
    for (int ir=0; ir<cs.GetNrows(); ++ir) {
      for (int ic=0; ic<cs.GetNcols(); ++ic) {
	//printf(" %6.4lf %8.4e %8.4e\n",cs[ir][ic],csErr[ir][ic],csErrSyst[ir][ic]);
	printf(" %8.4e %8.4e %8.4e\n",cs[ir][ic],csErr[ir][ic],csErrSyst[ir][ic]);
      }
      printf("\n");
    }
  }
  else {
    for (int ir=0; ir<cs.GetNrows(); ++ir) {
      for (int ic=0; ic<cs.GetNcols(); ++ic) {
	//printf(" %6.4lf %8.4e\n",cs[ir][ic],csErr[ir][ic]);
	printf(" %8.4e %8.4e\n",cs[ir][ic],csErr[ir][ic]);
      }
      printf("\n");
    }
  }
  std::cout << "done" << std::endl;
}

//------------------------------------------------------------------------------------------------------------------------

inline void printCSMatrixValues(const TString &name, const TMatrixD &cs, const TMatrixD &csErr, const TMatrixD &csSystErr, int printSystError=1) {
  printYields(name,cs,csErr,csSystErr,printSystError);
}

//------------------------------------------------------------------------------------------------------------------------

inline void printProgress(int printIfDivisiveBy, const char *msg, int idx, int idxMax) {
  if (idx%printIfDivisiveBy == 0) {
    double r=trunc(idx/double(idxMax)*1000)*0.1;
    std::cout << msg << idx << "/" << idxMax << " (" << r << "%)" << std::endl;
  }
}

//------------------------------------------------------------------------------------------------------------------------
/*
inline void printProgress(const char *msg, int idx, int idxMax) {
  double r=trunc(idx/double(idxMax)*1000)*0.1;
  std::cout << msg << idx << "/" << idxMax << " (" << r << "%)" << std::endl;
}
*/

//------------------------------------------------------------------------------------------------------------------------
/*
inline
void printSanityCheck(const TMatrixD &val, const TMatrixD &err, const TString &name)
{
  using namespace DYTools;
  std::cout<<"Sanity check printout"<<std::endl;
  if (val.GetNrows()!=nMassBins || val.GetNcols()!=findMaxYBins())
    {
      std::cout<<name<<" matrix has wrong size"<<std::cout;
      return;
    }
  if (err.GetNrows()!=nMassBins || err.GetNcols()!=findMaxYBins())
    {
      std::cout<<name<<"Err matrix has wrong size"<<std::cout;
      return;
    }
  std::cout<<"Nan values of "<<name<<" or/and "<<name<<"Err:"<<std::endl;
  for(int i=0; i<DYTools::nMassBins; i++)
    for (int yi=0; yi<nYBins[i]; ++yi)
      {
        if ( (val(i,yi)!=val(i,yi)) || (err(i,yi)!=err(i,yi)) )
           std::cout<<name<<"("<<i<<","<<yi<<")="<<val(i,yi)<<", "<<name<<"Err("<<i<<","<<yi<<")="<<err(i,yi)<<std::endl;
      }
  std::cout<<"Large errors ("<<name<<"Errv>0.1*"<<name<<" or "<<name<<"Err>0.1) :"<<std::endl;
  for(int i=0; i<DYTools::nMassBins; i++)
    for (int yi=0; yi<nYBins[i]; ++yi)
      {
        if ( fabs(err(i,yi))>0.1*fabs(val(i,yi)) || fabs(err(i,yi))>0.1)
           std::cout<<name<<"("<<i<<","<<yi<<")="<<val(i,yi)<<", "<<name<<"Err("<<i<<","<<yi<<")="<<err(i,yi)<<std::endl;
      }
}

//------------------------------------------------------------------------------------------------------------------------
*/

template<class Histo_t>
inline
int AppendToHistoName(Histo_t* h, TString add, int setAsTitle=0) {
  if (!h) {
    std::cout << "AppendToHistoName got null histo ptr\n";
    return 0;
  }
  TString newName=h->GetName() + add;
  h->SetName(newName);
  if (setAsTitle) h->SetTitle(newName);
  return 1;
}

//------------------------------------------------------------------------------------------------------------------------

template<class histo_t>
inline
void removeError1D(histo_t* h) {
  //std::cout << "nulifying error for " << h->GetName() << "\n";
  for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin) {
    h->SetBinError(ibin,0);
  }
}

//---------------------------------------------------------------

template<class histo_t>
inline
void removeError2D(histo_t* h) {
  //std::cout << "nulifying error for " << h->GetName() << "\n";
  for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=h->GetNbinsY(); ++jbin) {
      h->SetBinError(ibin,jbin,0);
    }
  }
}

//---------------------------------------------------------------

inline void removeError(TH1F *h) { removeError1D(h); }
inline void removeError(TH1D *h) { removeError1D(h); }
inline void removeError(TH2F *h) { removeError2D(h); }
inline void removeError(TH2D *h) { removeError2D(h); }

TH2D* removeUnderflow(TH2D* h, TString newName);

//------------------------------------------------------------------------------------------------------------------------

template<class histo_t>
inline int eliminateNaNs2D(histo_t *h, double setValue, double setError) {
  int count=0;
  for (int i=1; i<=h->GetNbinsX(); ++i) {
    for (int j=1; j<=h->GetNbinsY(); ++j) {
      if (h->GetBinContent(i,j)!=h->GetBinContent(i,j)) {
	h->SetBinContent(i,j, setValue);
	count++;
      }
      if (h->GetBinError(i,j)!=h->GetBinError(i,j)) {
	h->SetBinError(i,j, setError);
	count++;
      }
    }
  }
  return count;
}

//---------------------------------------------------------------

inline int eliminateNaNs(TH2D* h, double setValue, double setError=0.) { return eliminateNaNs2D(h,setValue,setError); }

//---------------------------------------------------------------

inline
int correctNegativeValues(TH2D* h) {
  int count=0;
  for (int i=1; i<=h->GetNbinsX(); ++i) {
      for (int j=1; j<=h->GetNbinsY(); ++j) {
	if (h->GetBinContent(i,j)<0) {
	  count++;
	  h->SetBinContent(i,j, 0.);
	}
      }
  }
  return count;
}

//------------------------------------------------------------------------------------------------------------------------

template<class histo_t>
inline
void squareBinContent1D(histo_t* h) {
  for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin) {
    double v=h->GetBinContent(ibin);
    h->SetBinContent(ibin, v*v);
  }
}

//---------------------------------------------------------------

template<class histo_t>
inline
void squareBinContent2D(histo_t* h) {
  for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=h->GetNbinsY(); ++jbin) {
      double v=h->GetBinContent(ibin,jbin);
      h->SetBinContent(ibin,jbin, v*v);
    }
  }
}

//---------------------------------------------------------------

inline void squareBinContent(TH2F *h) { squareBinContent2D(h); }
inline void squareBinContent(TH2D *h) { squareBinContent2D(h); }

//------------------------------------------------------------------------------------------------------------------------

template<class histo_t>
inline
void unsquareBinContent1D(histo_t* h) {
  for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin) {
    double v=h->GetBinContent(ibin);
    h->SetBinContent(ibin, sqrt(v));
  }
}

//---------------------------------------------------------------

template<class histo_t>
inline
void unsquareBinContent2D(histo_t* h) {
  for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=h->GetNbinsY(); ++jbin) {
      double v=h->GetBinContent(ibin,jbin);
      h->SetBinContent(ibin,jbin, sqrt(v));
    }
  }
}

//---------------------------------------------------------------

inline void unsquareBinContent(TH2F *h, int removeErr) {
  if (removeErr) removeError2D(h);
  unsquareBinContent2D(h);
}
inline void unsquareBinContent(TH2D *h, int removeErr) {
  if (removeErr) removeError2D(h);
  unsquareBinContent2D(h);
}

//---------------------------------------------------------------

inline void squareElements(TMatrixD &m) {
  for (int ir=0; ir<m.GetNrows(); ++ir) {
    for (int ic=0; ic<m.GetNcols(); ++ic) {
      const double x = m(ir,ic);
      m(ir,ic) = x*x;
    }
  }
}

//---------------------------------------------------------------

inline void unsquareElements(TMatrixD &m) {
  for (int ir=0; ir<m.GetNrows(); ++ir) {
    for (int ic=0; ic<m.GetNcols(); ++ic) {
      m(ir,ic) = sqrt(m(ir,ic));
    }
  }
}

//------------------------------------------------------------------------------------------------------------------------

template<class ptr_t>
inline void swapPtrs(ptr_t **ptr1, ptr_t **ptr2) {
  ptr_t* tmp=(*ptr1); (*ptr1)=(*ptr2); (*ptr2)=tmp;
}

//---------------------------------------------------------------

inline void swapHistoPtrs(TH1D **ptr1, TH1D **ptr2) {
  TH1D* tmp=(*ptr1); (*ptr1)=(*ptr2); (*ptr2)=tmp;
}
inline void swapHistoPtrs(TH2D **ptr1, TH2D **ptr2) {
  TH2D* tmp=(*ptr1); (*ptr1)=(*ptr2); (*ptr2)=tmp;
}

//---------------------------------------------------------------

template<class histo_t>
inline
void swapContentAndError1D(histo_t* h) {
  for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin) {
    double v=h->GetBinContent(ibin);
    h->SetBinContent(ibin, h->GetBinError(ibin));
    h->SetBinError(ibin, v);
  }
}

//---------------------------------------------------------------

template<class histo_t>
inline
void swapContentAndError2D(histo_t* h) {
  for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=h->GetNbinsY(); ++jbin) {
      double v=h->GetBinContent(ibin,jbin);
      h->SetBinContent(ibin,jbin, h->GetBinError(ibin,jbin));
      h->SetBinError(ibin,jbin, v);
    }
  }
}

//---------------------------------------------------------------

inline void swapContentAndError(TH2F *h) { swapContentAndError2D(h); }
inline void swapContentAndError(TH2D *h) { swapContentAndError2D(h); }

//---------------------------------------------------------------

inline int setErrorAsContent(TH2D *h) {
  if (!h) return 0;
  swapContentAndError2D(h);
  removeError2D(h);
  return 1;
}

//---------------------------------------------------------------

inline int setErrorAsContent(TH1D *h) {
  if (!h) return 0;
  swapContentAndError1D(h);
  removeError(h);
  return 1;
}

//---------------------------------------------------------------

inline
int setErrorAsContent(TH2D* hDest, const TH2D* hSrc) {
  if (!hDest || !hSrc ||
      (hDest->GetNbinsX() != hSrc->GetNbinsX()) ||
      (hDest->GetNbinsY() != hSrc->GetNbinsY())) {
    std::cout << "error in setErrorAsContent\n";
    return 0;
  }
  for (int ibin=1; ibin<=hDest->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=hDest->GetNbinsY(); ++jbin) {
      hDest->SetBinError(ibin,jbin, hSrc->GetBinContent(ibin,jbin));
    }
  }
  return 1;
}

//---------------------------------------------------------------

inline
int accumulateForRndStudies(TH2D* hSum, const TH2D* hAdd) {
  if (!hSum || !hAdd ||
      (hSum->GetNbinsX() != hAdd->GetNbinsX()) ||
      (hSum->GetNbinsY() != hAdd->GetNbinsY())) {
    std::cout << "error in accumulateForRndStudies\n";
    return 0;
  }
  for (int ibin=1; ibin<=hSum->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=hSum->GetNbinsY(); ++jbin) {
      double v=hSum->GetBinContent(ibin,jbin);
      double verr=hSum->GetBinError(ibin,jbin);
      double a=hAdd->GetBinContent(ibin,jbin);
      hSum->SetBinContent(ibin,jbin,v+a);
      hSum->SetBinError(ibin,jbin, sqrt(verr*verr + a*a));
    }
  }
  return 1;
}

//---------------------------------------------------------------

inline
int accumulateForRndStudies_finalize(TH2D* hSum, int nexps, int unbiasedEstimate=1) {
  //const int unbiasedEstimate=0;
  if (!hSum || (nexps<1+unbiasedEstimate)) {
    std::cout << "error in accumulateForRndStudies_finalize\n";
    return 0;
  }
  double factor=(unbiasedEstimate) ? double(nexps-1) : double(nexps);
  hSum->Scale(1/double(nexps));
  for (int ibin=1; ibin<=hSum->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=hSum->GetNbinsY(); ++jbin) {
      double avg=hSum->GetBinContent(ibin,jbin);
      double sqrtX2=hSum->GetBinError(ibin,jbin);
      // mean(x^2) - mean(x)^2
      double varSqr= nexps/factor*(nexps*sqrtX2*sqrtX2 - avg*avg);
      hSum->SetBinError(ibin,jbin, sqrt(varSqr));
    }
  }
  return 1;
}

//---------------------------------------------------------------
//---------------------------------------------------------------

TMatrixD* deriveCovMFromRndStudies(const std::vector<TH2D*> &rndVinp,
				   int unbiasedEstimate=1,
				   TH2D *avgDistr=NULL,
				   TMatrixD* covError=NULL);

// Calculate correlation matrix
TMatrixD* corrFromCov(const TMatrixD &cov);

// Extract error from matrix diagonal
TH2D* errorFromCov(const TMatrixD &cov, TString newName);

// Create covariance matrix from correlation matrix and errors
TMatrixD* covFromCorr(const TMatrixD &corr, const TVectorD &errs);
TMatrixD* covFromCorr(const TMatrixD &corr, const TH2D *h2Errs);

// Calculate partial correlation matrix
// assuming that cov is a part of totCov
TMatrixD* partialCorrFromCov(const TMatrixD &totCov, const TMatrixD &cov);

// Calculate relative covariance matrix
TMatrixD* relativeCov(const TVectorD &centralValue, const TMatrixD &cov);

TH2D* createHisto2D(const TMatrixD &M, const TMatrixD *Merr,
		    const char *histoName, const char *histoTitle,
		    TColorRange_t centerRange=_colrange_default,
		    int massBins=0,
		    double maxValUser=0.);

inline
TH2D* createHisto2D(const TMatrixD &M, const TMatrixD *Merr,
		    const char *histoName, const char *histoTitle,
		    int centerRange,
		    int massBins=0,
		    double maxValUser=0.) {
  return createHisto2D(M,Merr,histoName,histoTitle,
		       TColorRange_t(centerRange),massBins,maxValUser);
}

TMatrixD* createMatrixD(const TH2D *h2, int useErr=0);

// nom -> nom(i,j)/denom(i,j)
int divideMatrix(TMatrixD &nom, const TMatrixD &denom);

//----------------------------------------------------------------------------

TH1D* createHisto1D(const TVectorD &vec, const TVectorD* vecErr,
		    const char *histoName, const char* histoTitle=NULL,
		const TString &xAxisLabel="x", const TString &yAxisLabel="y");

//------------------------------------------------------------------------------------------------------------------------

//
// Produce a histogram that contains central values of the 'baseline' histogram
// and put the differences in the central values as an error of this central value
//

//TH2D *getRelDifferenceVA(const TH2D *baseValue, TString newName,
//		       int nVariations, TH2D *hVar1, ...);

// if includeVariants<>0, the differences among the variants are considered
TH2D *getRelDifference(const TH2D *baseValue, TString newName,
		       int includeVariants,
		       const TH2D *hVar1, const TH2D *hVar2=NULL,
		       const TH2D *hVar3=NULL, const TH2D *hVar4=NULL);

TH2D *getRelDifference(const std::vector<TH2D*> &vec, TString newName,
		       int includeVariants);
TH2D *getRelDifference(const std::vector<const TH2D*> &vec, TString newName,
		       int includeVariants);

TH2D* getRelError(const TH2D* h2, TString newName, int storeAsError=0);

//------------------------------------------------------------------------------------------------------------------------

// write flags to a file.
// Assumes that the file is open
int writeIntFlagValues(const TString &fieldName, int nFlags, int flag1, ...);
int writeFlagValues(const TString &fieldName, int nFlags, double flag1, ...);
int writeFlagValues(const TString &fieldName, int nFlags, const double *flags);
TVectorD* readFlagValues(TFile &fin, const TString &fieldName, int nFlags);

int writeFlagValues(const TString &fieldName, const std::vector<int> &flagsV);
int writeFlagValues(const TString &fieldName, const std::vector<double> &flagsV);

//------------------------------------------------------------------------------------------------------------------------

/*
inline
void subdivideBinWeightByLinearApprox(
   double m_1, double wsigma_1to0,
   double m0, double wsigma0to1, double wsigma0to1Err,
   double m1,
   double wsigma1to2, // needed if m2>0
   double m2,
   double mStar,
   double &wSigmaStar1, double &wSigmaStarErr1,
   double &wSigmaStar2, double &wSigmaStarErr2
   ) {
  // linear approximation
  // the cross section is per bin!
  // (m_1,m0, wsigma_1to0); (m0,m1, wsigma0to1) ; (m1,m2, wsigma1to2)
  // We want to divide (m0,m1) to (m0,mStar), (mStar,m2)
  // The cross section at the center of a bin
  // ( mc_05, wsigma_1to0/w0 ); ( mc05, wsigma0to1/w1 ) ; ( mc15, wsigma1to2/w2 )
  // or
  // ( mc_05, sigma_1to0 ) ; ( mc05, sigma0to1 ); ( mc15, sigma1to2 )
  // defines linear dependencies
  //  sigma= (sigma0to1 - sigma_1to0)/(mc05-mc_05) * ( m - mc_05 ) + sigma_1to0
  // and
  //  sigma= (sigma1to2 - sigma0to1)/(mc15-mc05) * ( m - mc05 ) + sigma0to1

  const int debug=0;
  const int warn_on_branches=0;
  const int pure_linear_branch=0; // recommended is 0 -- preserves total count in all cases
  const double thr=-9e8;  // if m2<thr, ignore sigma1to2

  if (debug) {
    std::cout << "subdivideBinWeightByLinearApprox: \n";
    std::cout << "  m_1=" << m_1 << ", m0=" << m0 << ", wsigma_1to0=" << wsigma_1to0 << "\n";
    std::cout << "  m0=" << m0 << ", m1=" << m1 << ", wsigma0to1=( " << wsigma0to1 << " pm " << wsigma0to1Err << " )\n";
    if (m2<thr) std::cout << " m2, wsigma1to2 are ignored\n";
    else std::cout << "  m1=" << m1 << ", m2=" << m2 << ", wsigma1to2=" << wsigma1to2 << "\n";
    std::cout << "  mStar=" << mStar << "\n";
    std::cout << "\n";
  }

  double mc_05=0.5*(m_1 + m0);
  double w0=m0-m_1;
  double w1=m1-m0;
  double sigma_1to0 = wsigma_1to0 / w0;
  double sigma0to1 =  wsigma0to1 / w1;
  double slope=(sigma0to1 - sigma_1to0)*2./(m1 - m_1);
  if (debug) std::cout << "slope=" << slope << " (1)\n";
  if ((m2>thr) && pure_linear_branch) {
    // Try to improve using slope with respect to the next bin
    double w2=m2-m1;
    double sigma1to2 =  wsigma1to2 / w2;
    double slopePlus= (sigma1to2 - sigma0to1) * 2./(m2 - m0);
    slope = 0.5*( slope + slopePlus );
    if (debug) std::cout << "slope=" << slope << " (2)\n";
  }

  if ((mStar>=m0) && (mStar<=m1)) {     // interpolation

    // center of the bin to consider, and the width
    double mcStar1=0.5*(m0 + mStar);
    double wStar1= mStar-m0;
    if (debug) std::cout << " mcStar1=" << mcStar1 << ", wStar1=" << wStar1 << "\n";
    // the cross section there
    double sigma1_div_width= slope*( mcStar1 - mc_05 ) + sigma_1to0;
    if (debug) std::cout << " sigma1_div_width=" << sigma1_div_width << "\n";
    // 1st answer
    wSigmaStar1 = sigma1_div_width * wStar1;

    double wStar2= m1-mStar;
    if (pure_linear_branch) {  // Pure linear approximation - use formula
      if (warn_on_branches) std::cout << " branch: pure linear approximation\n";
      // center to the other part of the bin, and the width
      double mcStar2=0.5*(mStar + m1);
      double sigma2_div_width= slope*( mcStar2 - mc_05 ) + sigma_1to0;
      // 2nd answer
      wSigmaStar2 = sigma2_div_width* wStar2;
    }
    else {  // Take the remaining counts
      if (warn_on_branches) std::cout << " branch: taking remaining area\n";
      wSigmaStar2= wsigma0to1 - wSigmaStar1;
    }

    // Error calculation
    double err=wsigma0to1Err/w1;
    double frac1=  wStar1 / sqrt( wStar1*wStar1 + wStar2*wStar2 );
    double frac2=  wStar2 / sqrt( wStar1*wStar1 + wStar2*wStar2 );
    wSigmaStarErr1= frac1*err *wStar1;
    wSigmaStarErr2= frac2*err *wStar2;

  }
  else { // extrapolation
    wSigmaStar1=0; wSigmaStarErr1=0;
    wSigmaStar2=0; wSigmaStarErr2=0;
    if (mStar>m1) {
      if (warn_on_branches) std::cout << " branch:  Extrapolation to the right\n";
      double mc01=0.5*(m0 + m1);
      double sigmaAtM1= slope*( m1 - mc01 ) + sigma0to1;
      double sigmaAtMStar= slope*( mStar - mc01 ) + sigma0to1;
      double w= mStar-m1;
      if (debug) std::cout << "sigmaAtM1=" << sigmaAtM1 << ", sigmaAtMStar=" << sigmaAtMStar << ", w=" << w << "\n";
      wSigmaStar1 = 0.5*(sigmaAtM1 + sigmaAtMStar) * w;
      wSigmaStarErr1= (wsigma0to1Err/w1) *w;
    }
    else {
      // this is not a branch -- this is an important warning
      std::cout << "\n\tExtrapolation to the left is not possible by construction\n";
    }
  }

  return;
}
*/
// --------------------------------
// --------------------------------
/*
inline
TH1D* createZpeakHisto(const char *hname="hZpeak_mass", const char *htitle="Z-peak region") {
  const int mbCount=13;
  const double mbins[mbCount+1]={  60,  64,  68,  72, 76,
				   81,  86,  91,  96,101,
				  106, 110, 115, 120 };
  TH1D *h= new TH1D(hname,htitle,mbCount,mbins);
  h->Sumw2();
  h->SetDirectory(0);
  std::cout << "mcCount=13\n";
  return h;
}

// --------------------------------

inline
TH1D* createZpeakHisto1(const char *hname="hZpeak_mass", const char *htitle="Z-peak region", double mass_min=60., double mass_max=120., int binCount=60) {
  TH1D *h= new TH1D(hname,htitle,binCount,mass_min,mass_max);
  h->Sumw2();
  h->SetDirectory(0);
  return h;
}

// --------------------------------

template<class TH1D_t>
inline
int printZpeakInfo(TH1D_t *h) {
  if (!h) return 0;
  std::cout << Form("(mean,mean_err)=(%6.4lf,%6.4lf)", h->GetMean(),h->GetMeanError())
            << " in " << h->GetName() << "\n";
  return 1;
}
*/

// --------------------------------

TString getTimeStrForPrint(Float_t tf);
void ShowBenchmarkTime(const char *clock_name);

// --------------------------------

inline
TString DayAndTimeTag(int eliminateSigns=1)
{
   time_t ltime;
   ltime=time(NULL);
   TString str = TString(asctime(localtime(&ltime)));
   if (str[str.Length()-1]=='\n') str.Remove(str.Length()-1,1);
   if (eliminateSigns) {
     str.ReplaceAll(" ","_");
     str.ReplaceAll(":","");
   }
   return str;
}

//------------------------------------------------------------------------------------------------------------------------

  // -------------------------------------------
    //  convert m[nMassBins][ybins] -> v[flat_idx]
inline
  int flattenMatrix(const TMatrixD &m, TVectorD &v) {
    for (int i=0; i<DYTools::nMassBins; ++i) {
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	int flatIdx=DYTools::findIndexFlat(i,yi);
	v[flatIdx]=m[i][yi];
      }
    }
    return 1;
  }

  // -------------------------------------------

  //  convert v[flat_idx] -> m[nMassBins][ybins]
inline
  int deflattenMatrix(const TVectorD &v, TMatrixD &m) {
    for (int i=0; i<DYTools::nMassBins; ++i) {
      for (int yi=0; yi<DYTools::nYBins[i]; ++yi) {
	int flatIdx=DYTools::findIndexFlat(i,yi);
	m[i][yi]=v[flatIdx];
      }
    }
    return 1;
  }


// -------------------------------------------
//  convert m[nMassBins][ybins] -> v[flat_idx]
TH1D* flattenHisto(const TH2D *h2, TString newName);
//  convert v[flat_idx] -> m[nMassBins][ybins]
TH2D* deflattenHisto(const TH1D *h, TString newName);


//------------------------------------------------------------------------------------------------------------------------


// -------------------------------------------

inline
TH1D* createBaseH1(const TString &histoName, const TString &histoTitle="", const TString &yAxisLabel="counts") {
  TH1D* h=new TH1D(histoName,histoTitle,
		   DYTools::nMassBins,DYTools::massBinLimits);
  h->SetDirectory(0);
  h->SetStats(0);
  h->Sumw2();
  h->GetXaxis()->SetTitle("M_{ee} [GeV]");
  h->GetYaxis()->SetTitle(yAxisLabel);
  return h;
}

// -------------------------------------------

inline
TH2D* createBaseH2(const TString &histoName, const TString &histoTitle="", int absRapidity=1) {
  double yMin=(absRapidity) ? DYTools::yRangeMin : -DYTools::yRangeMax;
  int yDivCount=(absRapidity) ? DYTools::nYBinsMax : 2*DYTools::nYBinsMax;
  TH2D* h2=new TH2D(histoName,histoTitle,
		    DYTools::nMassBins,DYTools::massBinLimits,
		    yDivCount,yMin, DYTools::yRangeMax);
  h2->SetDirectory(0);
  h2->SetStats(0);
  h2->Sumw2();
  h2->GetXaxis()->SetTitle("M_{ee} [GeV]");
  h2->GetYaxis()->SetTitle("y");
  return h2;
}

// -------------------------------------------
// -------------------------------------------

TH2D* convertBaseH2actual(const TH2D* h2, TString newHistoName, int setTitle=0);
TH1D* convertBinNo2Mass(const TH1D* hSrc);

TH1D* convertRange2BinNo(const TH1D* hSrc);
TH2D* convertRange2BinNo(const TH2D* hSrc);

TH1D* createProfileY(TH2D *h2, int ixBin, const TString &name, int setTitle=0, const char *title=NULL, int set_nYbins=-1, double set_ymin=0., double set_ymax=1.);
TH1D* createProfileX(TH2D *h2, int iyBin, const TString &name, int setTitle=0, const char *title=NULL);

TH1D* createProfileAuto(TH2D* h2, int iBin, const TString &name, int setTitle=0, const char *title=NULL);

int createRapidityProfileVec(const std::vector<TH2D*> &h2V, std::vector<std::vector<TH1D*>*> &hProfV, const std::vector<TString> &labels, int markerStyle=24, double markerSize=0.8);

// works for 1D or uniform binning
int createMassProfileVec(const std::vector<TH2D*> &h2V, std::vector<std::vector<TH1D*>*> &hProfV, const std::vector<TString> &labels, int markerStyle=24, double markerSize=0.8);

//int createAllProfileXVec(const std::vector<TH2D*> &h2V, const std::vector<TString> &names, int markerStyle=24, double markerSize=0.8);

/*
inline
TH1D* createProfileX(TH2D *h2, int iyBin, const TString &name, int setTitle=0, const char *title=NULL) {
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
*/

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
/*
inline
TH1D* createProfileY(TH2D *h2, int ixBin, const TString &name, int setTitle=0, const char *title=NULL, int set_nYbins=-1, double set_ymin=0., double set_ymax=1.) {
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
      std::cout << "yv[" << i-1 << "]=" << yv[i-1] << "\n";
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
*/

// -------------------------------------------
// -------------------------------------------

// For a fast switch from 41 to 40 mass bins

 TH1D* removeLastBin(const TH1D* hIni, TString newName, int setTitle=0, const char *newTitle=NULL);


// Modify histogram contents. Useful adjusting the assigned systematics
void linearApprox(TH1D *h, int ibin1, int ibin2,
		  double valOnEdge1=-1., double valOnEdge2=-1.);
void linearApprox(TH1D *h, double mass1, double mass2,
		  double valOnEdge1=-1., double valOnEdge2=-1.);

// -------------------------------------------

inline
void eliminateSeparationSigns(TString &name, int noDashForDot=0,
			      int dontChangeMinus=0) {
  name.ReplaceAll(" ","_");
  if (!dontChangeMinus) name.ReplaceAll("-","_");
  name.ReplaceAll("+","_");
  if (noDashForDot) name.ReplaceAll(".",""); else name.ReplaceAll(".","_");
  name.ReplaceAll(",","_");
  name.ReplaceAll(";","_");
  name.ReplaceAll("<","_");
  name.ReplaceAll(">","_");
  name.ReplaceAll("=","_");
  name.ReplaceAll("(","_");
  name.ReplaceAll(")","_");
}

// -------------------------------------------
// -------------------------------------------

inline
int createBaseH1Vec(std::vector<TH1D*> &histosV, const TString &histoNameBase, const std::vector<TString> &sample_labels, const TString &yAxisLabel="counts", int setHistoTitle=1) {
  int res=1;
  if (sample_labels.size()==0) return 0;
  histosV.reserve(sample_labels.size());
  for (unsigned int i=0; i<sample_labels.size(); ++i) {
    TString histoName= histoNameBase + sample_labels[i];
    TString histoTitle=(setHistoTitle) ? histoName : "";
    eliminateSeparationSigns(histoName);
    TH1D *h=createBaseH1(histoName,histoTitle,yAxisLabel);
    if (!h) { res=0; break; }
    histosV.push_back(h);
  }
  return res;
}

// -------------------------------------------

//
// Note that BaseH2 contains the division of (nMassBins,nYBinsMax)
//
inline
int createBaseH2Vec(std::vector<TH2D*> &histosV, const TString &histoNameBase, const std::vector<TString> &sample_labels, int absRapidity=1, int setHistoTitle=1) {
  int res=1;
  if (sample_labels.size()==0) return 0;
  histosV.reserve(sample_labels.size());
  for (unsigned int i=0; i<sample_labels.size(); ++i) {
    TString histoName= histoNameBase + sample_labels[i];
    TString histoTitle=(setHistoTitle) ? histoName : "";
    eliminateSeparationSigns(histoName);
    TH2D *h=createBaseH2(histoName,histoTitle,absRapidity);
    if (!h) { res=0; break; }
    histosV.push_back(h);
  }
  return res;
}

// -------------------------------------------

int convertBaseH2actualVec(const std::vector<TH2D*> &baseV, std::vector<TH2D*> &actualV, const TString histoNameBase, const std::vector<TString> &sample_labels, int setHistoTitle=1);

// -------------------------------------------

inline
TH1D* createAnyTH1D(const TString &hname, const TString &htitle, int nBins, double xMin, double xMax,
		 const TString &xAxisLabel="x", const TString &yAxisLabel="y") {
  TH1D *h=new TH1D(hname,htitle,nBins,xMin,xMax);
  h->SetDirectory(0);
  h->SetStats(0);
  h->Sumw2();
  h->GetXaxis()->SetTitle(xAxisLabel);
  h->GetYaxis()->SetTitle(yAxisLabel);
  return h;
}

// -------------------------------------------

inline
int createAnyH1Vec(std::vector<TH1D*> &histosV, const TString &histoNameBase,
		   const std::vector<TString> &sample_labels,
		   int nBins, const double xMin, const double xMax,
		   const TString &xAxisLabel="x",
		   const TString &yAxisLabel="counts", int setHistoTitle=1) {
  int res=1;
  if (sample_labels.size()==0) return 0;
  histosV.reserve(sample_labels.size());
  for (unsigned int i=0; i<sample_labels.size(); ++i) {
    TString histoName= histoNameBase + sample_labels[i];
    TString histoTitle=(setHistoTitle) ? histoName : "";
    eliminateSeparationSigns(histoName);
    TH1D *h=new TH1D(histoName,histoTitle,nBins,xMin,xMax);
    //std::cout << "new histo: " << h->GetName() << "\n";
    h->SetDirectory(0);
    h->SetStats(0);
    h->Sumw2();
    h->GetXaxis()->SetTitle(xAxisLabel);
    h->GetYaxis()->SetTitle(yAxisLabel);
    if (!h) { res=0; break; }
    histosV.push_back(h);
  }
  return res;
}

// -------------------------------------------

inline
int createAnyH1Vec(std::vector<TH1D*> &histosV, const TString &histoNameBase,
		   const std::vector<TString> &sample_labels,
		   int nBins, const double *binEdges,
		   const TString &xAxisLabel="x",
		   const TString &yAxisLabel="counts", int setHistoTitle=1) {
  int res=1;
  if (sample_labels.size()==0) return 0;
  histosV.reserve(sample_labels.size());
  for (unsigned int i=0; i<sample_labels.size(); ++i) {
    TString histoName= histoNameBase + sample_labels[i];
    TString histoTitle=(setHistoTitle) ? histoName : "";
    eliminateSeparationSigns(histoName);
    TH1D *h=new TH1D(histoName,histoTitle,nBins,binEdges);
    //std::cout << "new histo: " << h->GetName() << "\n";
    h->SetDirectory(0);
    h->SetStats(0);
    h->Sumw2();
    h->GetXaxis()->SetTitle(xAxisLabel);
    h->GetYaxis()->SetTitle(yAxisLabel);
    if (!h) { res=0; break; }
    histosV.push_back(h);
  }
  return res;
}

// -------------------------------------------

inline
int createAnyH2Vec(std::vector<TH2D*> &histosV, const TString &histoNameBase,
		   const std::vector<TString> &sample_labels,
		   int nxBins, double xMin, double xMax,
		   int nyBins, double yMin, double yMax,
		   const TString &xAxisLabel="x",
		   const TString &yAxisLabel="y", int setHistoTitle=1) {
  int res=1;
  if (sample_labels.size()==0) return 0;
  histosV.reserve(sample_labels.size());
  for (unsigned int i=0; i<sample_labels.size(); ++i) {
    TString histoName= histoNameBase + sample_labels[i];
    TString histoTitle=(setHistoTitle) ? histoName : "";
    eliminateSeparationSigns(histoName);
    TH2D *h=new TH2D(histoName,histoTitle,
		     nxBins,xMin,xMax,
		     nyBins,yMin,yMax);
    h->SetDirectory(0);
    h->SetStats(0);
    h->Sumw2();
    h->GetXaxis()->SetTitle(xAxisLabel);
    h->GetYaxis()->SetTitle(yAxisLabel);
    if (!h) { res=0; break; }
    histosV.push_back(h);
  }
  return res;
}

//------------------------------------------------------------------------------------------------------------------------

inline
int sameNumBins(const TH1D* h1, const TH1D* h2) {
  return (h1->GetNbinsX() == h2->GetNbinsX()) ? 1:0;
}

//-----------------------------------------------------------

inline
int sameNumBins(const TH2D* h1, const TH2D* h2) {
  return ((h1->GetNbinsX() == h2->GetNbinsX()) &&
	  (h2->GetNbinsY() == h2->GetNbinsY())) ? 1:0;
}

//-----------------------------------------------------------

inline
int sameNumBins(const TMatrixD &a, const TMatrixD &b) {
  return ((a.GetNrows() == b.GetNrows()) &&
	  (a.GetNcols() == b.GetNcols())) ? 1:0;
}


//-----------------------------------------------------------
//-----------------------------------------------------------

template<class histo_t>
inline
int addHistos(histo_t *sum, const std::vector<histo_t*> &vec) {
  for (unsigned int i=0; i<vec.size(); ++i) sum->Add(vec[i]);
  return 1;
}

//-----------------------------------------------------------

template<class histo_t>
inline
histo_t* addHistos(TString newName, const std::vector<histo_t*> &vec) {
  if (vec.size()==0) return NULL;
  histo_t *h=(histo_t*)vec[0]->Clone(newName);
  h->SetDirectory(0);
  h->SetTitle(newName);
  for (unsigned int i=0; i<vec.size(); ++i) h->Add(vec[i]);
  return h;
}

//-----------------------------------------------------------

// error is ignored
int scaleHisto(TH1D *histoNom, const TH1D *histoDenom, int mult=0);
int scaleHisto(TH2D *histoNom, const TH2D *histoDenom, int mult=0);

template<class histo_t>
inline int divide(histo_t *histoNom, const histo_t *histoDenom)
{ return scaleHisto(histoNom,histoDenom,0); }

// error is NOT ignored
int multiplyHisto(TH1D *histoNom, const TH1D *histoDenom, int mult=1);
int multiplyHisto(TH2D *histoNom, const TH2D *histoDenom, int mult=1);

TH1D* convert_TH1F_to_TH1D(const TH1F *h, TString newName);
TH1F* convert_TH1D_to_TH1F(const TH1D *h, TString newName);
TH2D* convert_TH2F_to_TH2D(const TH2F *h, TString newName);
TH2F* convert_TH2D_to_TH2F(const TH2D *h, TString newName);

TH1D* convert_absolute2rshape(const TH1D *h, double norm);
TH1D* convert_rshape2absolute(const TH1D *h, double norm);

//-----------------------------------------------------------

inline void printHisto(const TH1F *h1F) {
  TString hName=TString(h1F->GetName()) + TString("_");
  TH1D *h=convert_TH1F_to_TH1D(h1F,hName);
  printHisto(h);
}

//------------------------------------------------------------------------------------------------------------------------

template<class tObject_t>
inline
int saveVec(TFile &file, const std::vector<tObject_t*> &vec, const TString &subDir="") {
  file.cd();
  if (subDir.Length()) {
    file.mkdir(subDir);
    file.cd(subDir);
  }
  for (unsigned int i=0; i<vec.size(); i++) {
    TString name=vec[i]->GetName();
    eliminateSeparationSigns(name);
    vec[i]->Write(name);
  }
  return 1;
}

//--------------------------------------------------

template<class tObject_t>
inline
int loadVec(TFile &file, std::vector<tObject_t*> &vec, TString subDir="") {
  file.cd();
  if (vec.size()==0) {
    HERE("loadVec(%s,vec,%s): vec.size=0",TString(file.GetName()),subDir);
    return 0;
  }
  int isHisto=(vec[0]->InheritsFrom("TH1")) ? 1:0;
  //if (!isHisto) isHisto=(vec[0]->InheritsFrom("TH2")) ? 1:0;
  if (isHisto && subDir.Length()) {
    if (!file.cd(subDir)) {
      std::cout << "loadVec: failed to change subdir=<" << subDir << ">\n";
      return 0;
    }
  }
  if (!isHisto && subDir.Length() && (subDir[subDir.Length()-1]!='/')) {
    subDir.Append("/");
  }
  for (unsigned int i=0; i<vec.size(); i++) {
    TString name=vec[i]->GetName();
    if (!isHisto && subDir.Length()) name.Prepend(subDir);
    eliminateSeparationSigns(name);
    //std::cout << "reading " << name << std::endl;

    vec[i]->Read(name);
    if (!vec[i]) {
      HERE("loaded null ptr for %s",name);
      return 0;
    }
    vec[i]->SetDirectory(0);
  }
  if (subDir.Length()) file.cd();
  HERE("load ok");
  return 1;
}

//--------------------------------------------------

template<class histo_t>
inline
int saveHisto(TFile &file, histo_t *h, const TString &subDir="", const TString &saveWithName="") {
  file.cd();
  if (subDir.Length() && !file.cd(subDir)) {
    std::cout << " ... creating subdirectory <" << subDir << ">\n";
    file.mkdir(subDir);
    file.cd(subDir);
  }
  TString hname=(saveWithName.Length()) ? saveWithName.Data() : h->GetName();
  h->Write(hname);
  if (subDir.Length()) file.cd();
  return 1;
}

//--------------------------------------------------

template<class histo_t>
inline
int loadHisto(TFile &file, histo_t **h, TString subDir) {
  file.cd();
  if (subDir.Length()) file.cd(subDir);
  TString name=(*h)->GetName();
  eliminateSeparationSigns(name);
  if (!(*h)->Read(name)) {
    std::cout << "failed to load <" << name << ">\n";
    return 0;
  }
  if (*h) {
    (*h)->SetDirectory(0);
  }
  else return 0;
  if (subDir.Length()) file.cd();
  return 1;
}

//--------------------------------------------------

// load from un-opened file
TH1D* LoadHisto1D(TString fieldName, const TString &fname, TString subDir="", int checkBinning=1);

// load from un-opened file
TH2D* LoadHisto2D(TString fieldName, const TString &fname, TString subDir="", int checkBinning=1);
// load from an opened file
TH2D* LoadHisto2D(TFile &fin, TString fieldName, TString subDir="", int checkBinning=1);

//--------------------------------------------------
//--------------------------------------------------

inline
TH1D* Clone(const TH1D* histo, const TString &newName, const TString &newTitle) {
  TH1D *h2=(TH1D*)histo->Clone(newName);
  if (!h2) {
    std::cout << "failed to clone a histo" << std::endl;
  }
  else {
    h2->SetDirectory(0);
    //h2->Sumw2();
    h2->SetTitle(newTitle);
  }
  return h2;
}

//--------------------------------------------------

inline
TH1D* Clone(const TH1D* histo, const TString &newName, int setNameAsTitle=1) {
  return Clone(histo,newName,(setNameAsTitle) ? newName : TString(""));
}

//--------------------------------------------------

inline
TH2D* Clone(const TH2D* histo, const TString &newName, const TString &newTitle) {
  TH2D *h2=(TH2D*)histo->Clone(newName);
  if (!h2) {
    std::cout << "failed to clone a histo" << std::endl;
  }
  else {
    h2->SetDirectory(0);
    //h2->Sumw2();
    h2->SetTitle(newTitle);
  }
  return h2;
}

//--------------------------------------------------

inline
TH2D* Clone(const TH2D* histo, const TString &newName, int setTitle=0) {
  TH2D *h2=(TH2D*)histo->Clone(newName);
  h2->SetDirectory(0);
  //h2->Sumw2();
  if (setTitle) h2->SetTitle(newName);
  return h2;
}

//--------------------------------------------------

inline
int Copy(const TH1D* hSrc, TH1D* hDest) {
  if (!sameNumBins(hSrc,hDest)) {
    std::cout << "Copy(TH1D) : different sizes\n";
    return 0;
  }
  for (int ibin=1; ibin<=hSrc->GetNbinsX(); ++ibin) {
    hDest->SetBinContent(ibin, hSrc->GetBinContent(ibin));
    hDest->SetBinError  (ibin, hSrc->GetBinError(ibin));
  }
  return 1;
}

//--------------------------------------------------

//--------------------------------------------------

inline
int Copy(const TH2D* hSrc, TH2D* hDest) {
  if (!sameNumBins(hSrc,hDest)) {
    std::cout << "Copy(TH2D) : different sizes\n";
    return 0;
  }
  for (int ibin=1; ibin<=hSrc->GetNbinsX(); ++ibin) {
    for (int jbin=1; jbin<=hSrc->GetNbinsY(); ++jbin) {
      hDest->SetBinContent(ibin,jbin, hSrc->GetBinContent(ibin,jbin));
      hDest->SetBinError  (ibin,jbin, hSrc->GetBinError(ibin,jbin));
    }
  }
  return 1;
}

//--------------------------------------------------

// Histo has bins (1,Nx) x (1,Ny).
// Select a histogram (xbin1,xbin2) x (ybin1,ybin2).
//   If the axes are divided at integers (e.g. indices are used), the option
// resetAxis=1 might be useful. It ignores
// the actual values of the original axes.
TH2D* extractSubArea(const TH2D *histoSrc,
		     int xbin1, int xbin2, int ybin1, int ybin2,
		     const TString &newName, int setTitle=0, int resetAxis=0);

//--------------------------------------------------

inline
TH2D* clipToAnalysisUnfBins(const TH2D *histoSrc, TString newName,
			    TString histoTitle, int resetAxis=0) {
  int rangeMin=(DYTools::study2D) ? (DYTools::_nBinsYLowMass+1)  : 1;
  int rangeMax=
    (DYTools::study2D) ? DYTools::nUnfoldingBins : DYTools::nMassBins;
  TH2D* h2=extractSubArea(histoSrc,rangeMin,rangeMax,rangeMin,rangeMax,
			  newName,0,resetAxis);
  if (histoTitle.Length()) h2->SetTitle(histoTitle);
  if (!h2) std::cout << "error in clipToAnalysisUnfBins\n";
  return h2;
}

//--------------------------------------------------
//--------------------------------------------------

inline
double ZpeakCount(const TH2D* h2, double *err=NULL) {
  const int ZpeakIdx1=DYTools::findMassBin(60.) +1;
  const int ZpeakIdx2=DYTools::findMassBin(120.);
  const int yIdx1=1;
  const int yIdx2=DYTools::nYBins[ZpeakIdx1];

  double val=0;
  if (!err) { val=h2->Integral(ZpeakIdx1,ZpeakIdx2,yIdx1,yIdx2); }
  else { val=h2->IntegralAndError(ZpeakIdx1,ZpeakIdx2,yIdx1,yIdx2,*err); }
  return val;
}

//--------------------------------------------------

inline
void printZpeakCount(const TH2D *h2) {
  double err=0;
  std::cout << "ZpeakCount(" << h2->GetName() << ")="
	    << ZpeakCount(h2,&err);
  std::cout << " +/- " << err << "\n";
}

//--------------------------------------------------

inline
double ZpeakCount(const TH1D* h, double *err=NULL) {
  const int ZpeakIdx1=DYTools::findMassBin(60.) +1;
  const int ZpeakIdx2=DYTools::findMassBin(120.);

  double val=0;
  if (!err) { val=h->Integral(ZpeakIdx1,ZpeakIdx2); }
  else { val=h->IntegralAndError(ZpeakIdx1,ZpeakIdx2,*err); }
  return val;
}

//--------------------------------------------------

inline
void printZpeakCount(const TH1D *h1) {
  double err=0;
  std::cout << "ZpeakCount(" << h1->GetName() << ")="
	    << ZpeakCount(h1,&err);
  std::cout << " +/- " << err << "\n";
}

//--------------------------------------------------
//--------------------------------------------------

void writeBinningArrays(TFile &fout, TString producedBy="", int allInfo=1);
int checkBinningArrays(TFile &fin, int printMetaData=0);
int checkBinningRanges(const TVectorD &massBinEdges, const TVectorD &rapidityCounts, const TString &fname);
int checkMatrixSize(const TMatrixD &m, const TString &infoName);

//--------------------------------------------------
//--------------------------------------------------

inline int AsInt(const std::string &s) { return atoi(s.c_str()); }
inline double AsDouble(const std::string &s) { return atof(s.c_str()); }

//--------------------------------------------------
//--------------------------------------------------

TH2D* LoadMatrixFields(TFile &fin, const TString &field, const TString &fieldErr, int loadErr=1, int absoluteRapidity=1);
TH2D* LoadMatrixFields(const TString &fname, int checkBinning, const TString &field, const TString &fieldErr, int loadErr=1, int absoluteRapidity=1);

// LoadThreeMatrices loads (value,valueErr) into h2, and (valueSystErr,value^2) into h2syst

// Load one matrix.
// If reportFieldError=0, no msg is printed, although null ptr is returned
TMatrixD* loadMatrix(const TString &fname, const TString &fieldName, int expect_nRows, int expect_nCols, int reportFieldError=1);

int LoadThreeMatrices(TFile &fin, TH2D **h2, TH2D **h2syst, const TString &field, const TString &fieldErr, const TString &fieldSystErr, int absoluteRapidity=1);

int LoadThreeMatrices(const TString &fileName, TH2D **h2, TH2D **h2syst, const TString &field, const TString &fieldErr, const TString &fieldSystErr, int checkBinning, int absoluteRapidity=1);

//--------------------------------------------------
//--------------------------------------------------

inline
void CreateDir(const TString &fname, int printDir=1) {
  Ssiz_t len=fname.Last('/');
  TString dir=(len==-1) ? fname : fname(0,len);
  if (printDir) std::cout << "dir=" << dir << "\n";
  gSystem->mkdir(dir,kTRUE);
}

//--------------------------------------------------

void prepare(int count,
	     std::vector<TString> &pathV,
	     std::vector<TString> &fnameV,
	     std::vector<TString> &fieldV,
	     std::vector<TString> &labelV,
	     int clear=1,
	     int addEmptyElements=1);

void prepare(int count,
	     std::vector<TString> &pathV,
	     std::vector<TString> &fieldV,
	     std::vector<TString> &labelV,
	     int clear=1,
	     int addEmptyElements=1);


//--------------------------------------------------
//--------------------------------------------------

int saveLatexTable(TString fileTag,
		   const std::vector<TH2D*> &histosV,
		   const std::vector<TString> &labelsV,
		   const char *format,
		   int printErrors,
		   int savePlainFormat);

//--------------------------------------------------

int loadLatexTableTextFile(TString fileTag,
			   std::vector<TH2D*> &histosV,
			   std::vector<TString> &labelsV,
			   int printErrors);

//--------------------------------------------------
//--------------------------------------------------

TCanvas* compareHistos(TH2D* hA, TString labelA,
		       TH2D* hB, TString labelB,
		       int do_removeError=1,
		       TString yAxisLabel="signal_yield");

//--------------------------------------------------

TCanvas* plotProfiles(TString canvName,
		      TH2D* histo, TString label,
		      int do_removeError=1,
		      TString yAxisLabel="signal_yield");

//--------------------------------------------------

TCanvas* plotProfiles(TString canvName,
		      const std::vector<TH2D*> &histosV,
		      const std::vector<TString> &labelsV,
		      std::vector<int> *colorsV=NULL,
		      int do_removeError=1,
		      TString yAxisLabel="signal yield",
		      std::vector<std::vector<TH1D*>*> *hProfV=NULL,
		      std::vector<ComparisonPlot_t*> *cpV=NULL,
		      int delayDraw=0);

//--------------------------------------------------
//--------------------------------------------------


#endif
