#ifndef CrossSection_HH
#define CrossSection_HH

#include <TROOT.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <math.h>

#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
//#include "../Include/UnfoldingMatrix.h"

#ifdef DYee
//#include "../Include/UnfoldingTools.hh"
#include "../Include/ComparisonPlot.hh"
#endif

#include "../Include/BaseClass.hh"

class VXSectD_t;
class MXSectD_t;

//------------------------------------------------------------------------------------------------------------------------

#ifndef DYTools_HH
namespace DYTools {
  const int study2D=0;

  const int _nMassBins2011 = 40;
  const double _massBinLimits2011[_nMassBins2011+1] = 
    {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76, 
     81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141, 
     150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440, 
     510, 600, 1000, 1500}; // 40 bins
  const int _nYBinsMax2011=1; // the largest division into Y bins
  const int _nYBins2011[_nMassBins2011] = { 
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1
  };

  // Constants that define binning in mass and rapidity
  // Note: bin zero is underflow, overflow is neglected
  const int _nMassBins2D = 7;
  const double _massBinLimits2D[_nMassBins2D+1] = 
    {0, // first bin is underflow
     20, 30, 45, 60, 120, 200, 1500
    }; // overflow is very unlikely, do not account for it
  // Rapidity binning is different for different mass bins
  // Note: this implementation neglects underflow and overflow
  // in rapidity.
  const double yRangeMin =  0.0;
  const double yRangeMax =  10.0;
  const int _nYBinsMax2D=25; // the largest division into Y bins
  const int _nYBins2D[_nMassBins2D] = 
    { 25,// underflow, binned like first mass bin 
      25, 25, 25, 25, 25, 10,
    }; // overflow is neglected

  const int nMassBins= (study2D) ? _nMassBins2D : _nMassBins2011;
  const double *massBinLimits= (study2D) ? _massBinLimits2D : _massBinLimits2011;
  const int nYBinsMax= (study2D) ? _nYBinsMax2D : _nYBinsMax2011;
  const int *nYBins= (study2D) ? _nYBins2D : _nYBins2011;
  const int nUnfoldingBins= (study2D) ? 
    ((_nMassBins2D-1) * _nYBinsMax2D + _nYBins2D[_nMassBins2D-1]) : nMassBins;

  // ------

  // This function finds a unique 1D index for (index_m, index_y) pair
  inline 
  int findIndexFlat(int massBin, int yBin){
    int result = -1;
    if( massBin < 0 || massBin > nMassBins || yBin < 0 || yBin > nYBins[massBin] )
      return result;
    
    result = 0;
    for(int i=0; i< massBin; i++) result += nYBins[i];
    result += yBin;
    return result;
  }

  // ------
  /*
  // This function finds a unique 1D index for (index_m, index_y) pair
  inline 
  int findIndexFlat(double mass, double y){
    int massBin=findMassBin(mass);
    int yBin=findAbsYBin(massBin,y);
    return findIndexFlat(massBin,yBin);
  }
  */

};

#endif

//------------------------------------------------------------------------------------------------------------------------


#if !defined(UnfoldingTools_HH) && !defined(MYTOOLS_HH)
namespace unfolding{

  inline
  int unfold(const TVectorD &vin, TVectorD &vout, 
	     const TString &unfoldingConstFileName, 
	     const TString &unfMatrixName="DetInvertedResponse", 
	     int transpose=1) {
    TFile fileConstants(unfoldingConstFileName);
    TMatrixD *DetInvertedResponse_ptr     = (TMatrixD *)fileConstants.FindObjectAny(unfMatrixName);
    fileConstants.Close();
    if (!DetInvertedResponse_ptr) {
      std::cout << "failed to load <" << unfMatrixName << "> from <" << unfoldingConstFileName << ">\n";
      return 0;
    }
    TMatrixD DetInvertedResponse= *DetInvertedResponse_ptr;
    delete DetInvertedResponse_ptr;
    if (transpose) {
      TMatrixD tmp(DetInvertedResponse);
      DetInvertedResponse.Transpose(tmp);
    }

    vout = 0;
    int nBins = vin.GetNoElements();
    for(int i=0; i<nBins; i++){
      for(int j=0; j<nBins; j++){
	vout[i] += DetInvertedResponse(i,j) * vin[j];
      }
    }
    return 1;
  }

  // -------------------------------------------

inline
 int  propagateErrorThroughUnfolding(const TVectorD &errorIn, 
				     TVectorD &errorPropagated,
				     const TString &unfoldingConstFileName,
		     const TString &unfMatrixName="DetInvertedResponse", 
				     int transpose=1) {

    // Read unfolding constants
    std::cout << "propagateErrorThroughUnfolding: Load constants from <" << unfoldingConstFileName << ">" << std::endl;
    
    TFile fileConstants(unfoldingConstFileName);
    TMatrixD DetInvertedResponse     = *(TMatrixD *)fileConstants.FindObjectAny(unfMatrixName);
    fileConstants.Close();
    if (transpose) {
      TMatrixD tmp(DetInvertedResponse);
      DetInvertedResponse.Transpose(tmp);
    }

    // Apply unfolding matrix
    int nBins = errorIn.GetNoElements();
    errorPropagated = 0;
    for(int i=0; i<nBins; i++){
      for(int j=0; j<nBins; j++){
	errorPropagated[i] += pow( DetInvertedResponse   (i,j) * errorIn[j], 2);
      }
      errorPropagated[i] = sqrt(errorPropagated[i]);
    }
    return 1;
  }

  // -------------------------------------------
    //  convert m[nMassBins][ybins] -> v[flat_idx]
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

  int unfold(const TMatrixD &inpM, TMatrixD &outM,
	     const TString &unfConstFName,
	     TVectorD &inpV, TVectorD &outV) {
    int res=
      (flattenMatrix(inpM,inpV) &&
       unfold(inpV,outV,unfConstFName) &&
       deflattenMatrix(outV,outM)) ? 1:0;
    return res;
  }

  // -------------------------------------------  

  int propagateErrorThroughUnfolding(const TMatrixD &inpM, TMatrixD &outM,
	     const TString &unfConstFName,
	     TVectorD &inpV, TVectorD &outV) {
    int res=
      (flattenMatrix(inpM,inpV) &&
       propagateErrorThroughUnfolding(inpV,outV,unfConstFName) &&
       deflattenMatrix(outV,outM)) ? 1:0;
    return res;
  }

  // -------------------------------------------  

  void writeBinningArrays(TFile &fout) {
    using namespace DYTools;
    fout.cd();
    TVectorD mass(nMassBins+1);
    TVectorD rapidityCounts(nMassBins);
    for (int i=0; i<nMassBins+1; i++) mass[i]=massBinLimits[i];
    for (int i=0; i<nMassBins  ; i++) rapidityCounts[i]=nYBins[i];
    mass.Write("massBinLimits"); // was 'massBinning'
    rapidityCounts.Write("rapidityCounts");
  }

  // -------------------------------------------------------

  bool checkBinningArrays(TFile &fin) {
    const char *fncname="unfolding::checkBinningArrays: ";
    TString fileInfo=TString("on file <") + fin.GetName() + ">";
    fin.cd();
    TVectorD* mass= (TVectorD*)fin.FindObjectAny("massBinLimits");
    if (!mass) mass=(TVectorD*)fin.FindObjectAny("massBinning"); // back-port
    //if (!mass) mass=(TVectorD*)fin.FindObjectAny("mass"); // DDBkgr-port
    TVectorD* rapidityCounts= (TVectorD*)fin.FindObjectAny("rapidityCounts");
    if (!mass || !rapidityCounts) {
      std::cout << fncname << "file <" << fin.GetName() << "> does not contain keys massBinning and/or rapidityCounts\n";
      return false;
    }
    bool comparisonOk=true;
    //comparisonOk=checkBinningRanges(*mass,*rapidityCounts, fin.GetName());
    std::cout << "Binning ranges are not checked!\n (this is a restricted version)\n";
    if (mass) delete mass;
    if (rapidityCounts) delete rapidityCounts;
    return comparisonOk;
  }

  // -------------------------------------------------------
};
#endif
//------------------------------------------------------------------------------------------------------------------------


// =====================================

template<class Container_t, class objType_t>
int LoadToContainer_implicitFile(
    Container_t &store, objType_t &obj, int &hasSystErr,
    const TString &objName, const TString &objErrName, 
    const TString &objSystErrName) {
  if (0) obj.Print(); // get rid of compiler complaints
  hasSystErr=0;
  if (objName.Length()) {
    objType_t *ptr= (objType_t*) gFile->Get(objName);
    if (!ptr) { 
      std::cout << "failed to get <" << objName << "> from <" << gFile->GetName() << ">\n";  
      return 0; 
    }
    store.EditXS() = * ptr;
    delete ptr;
  }
  else store.EditXS().Zero();
  if (objErrName.Length()) {
    objType_t *ptr= (objType_t*) gFile->Get(objErrName);
    if (!ptr) { 
      std::cout << "failed to get <" << objErrName << "> from <" << gFile->GetName() << ">\n";  
      return 0; 
    }
    store.EditXSerr() = *ptr;
    delete ptr;
  }
  else store.EditXSerr().Zero();
  if (objSystErrName.Length()) {
    objType_t *ptr= (objType_t*) gFile->Get(objSystErrName);
    if (!ptr) { 
      std::cout << "failed to get <" << objSystErrName << "> from <" << gFile->GetName() << ">\n";  
      return 0; 
    }
    store.EditXSsystErr() = * ptr;
    delete ptr;
    hasSystErr=1;
  }
  else store.EditXSsystErr().Zero();
  return 1;
}

// =====================================

template<class Container_t, class objType_t>
int LoadToContainerFromFile(
    Container_t &store, objType_t &obj, int &hasSystErr,
    TFile &fin, 
    const TString &objName, const TString &objErrName, 
    const TString &objSystErrName) {
  if (0) obj.Print(); // get rid of compiler complaints
  hasSystErr=0;
  if (!fin.IsOpen()) return 0;
  if (objName.Length()) {
    objType_t *ptr= (objType_t*) fin.Get(objName);
    if (!ptr) { 
      std::cout << "failed to get <" << objName << "> from <" << fin.GetName() << ">\n";  
      return 0; 
    }
    store.EditXS() = * ptr;
    delete ptr;
  }
  else store.EditXS().Zero();
  if (objErrName.Length()) {
    objType_t *ptr= (objType_t*) fin.Get(objErrName);
    if (!ptr) { 
      std::cout << "failed to get <" << objErrName << "> from <" << fin.GetName() << ">\n";  
      return 0; 
    }
    store.EditXSerr() = *ptr;
    delete ptr;
  }
  else store.EditXSerr().Zero();
  if (objSystErrName.Length()) {
    objType_t *ptr= (objType_t*) fin.Get(objSystErrName);
    if (!ptr) { 
      std::cout << "failed to get <" << objSystErrName << "> from <" << fin.GetName() << ">\n";  
      return 0; 
    }
    store.EditXSsystErr() = * ptr;
    delete ptr;
    hasSystErr=1;
  }
  else store.EditXSsystErr().Zero();
  return 1;
}


// ----------

template<class Container_t, class objType_t>
int LoadToContainer(Container_t &store, objType_t &obj,
		    int &hasSystErr,
	 const TString &fname, 
	 const TString &objName, const TString &objErrName,
	 const TString &objSystErrName) {
  //HERE("LoadToContainer opening file <%s>",fname);
  TFile fin(fname,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to load file <" << fname << ">\n";
    return 0;
  }
  //HERE("file <%s> opened. Reading it",fname);
  if (!LoadToContainerFromFile(store,obj,hasSystErr,
		       fin,objName,objErrName,objSystErrName)) 
    return 0;
  fin.Close();
  return 1;
}

// ----------

template<class Container_t, class objType_t>
int SaveContainer_implicitFile(
    const Container_t &store, const objType_t &obj,
    const TString &objName, const TString &objErrName, 
    const TString &objSystErrName) {
  if (0) obj.Print(); // get rid of compiler complaints
  if (objName.Length()) store.GetXS().Write(objName);
  if (objErrName.Length()) store.GetXSerr().Write(objErrName);
  if (objSystErrName.Length()) {
    if (store.HasSystErr()==0) {
      std::cout << "object store.FName=<" << store.GetName() << "> has no syst.Err, but the field name is provided (" << objSystErrName << ")\n";
    }
    store.GetXSsystErr().Write(objSystErrName);
  }
  return 1;
}


// ----------

template<class Container_t, class objType_t>
int SaveContainer(
    const Container_t &store, const objType_t &obj,
    TFile &fout, 
    const TString &objName, const TString &objErrName, 
    const TString &objSystErrName) {
  if (0) obj.Print(); // get rid of compiler complaints
  if (!fout.IsOpen()) return 0;
  if (objName.Length()) store.GetXS().Write(objName);
  if (objErrName.Length()) store.GetXSerr().Write(objErrName);
  if (objSystErrName.Length()) {
    if (store.HasSystErr()==0) {
      std::cout << "object store.FName=<" << store.GetName() << "> has no syst.Err, but the field name is provided (" << objSystErrName << ")\n";
    }
    store.GetXSsystErr().Write(objSystErrName);
  }
  return 1;
}


// ----------

template<class Container_t, class objType_t>
int SaveContainer(const Container_t &store, const objType_t &obj,
	 const TString &fname, 
	 const TString &objName, const TString &objErrName,
	 const TString &objSystErrName) {
  TFile fout(fname,"recreate");
  if (!fout.IsOpen()) {
    std::cout << "failed to (re)create file <" << fname << ">\n";
    return 0;
  }
  if (!SaveContainer(store,obj,
		     fout,objName,objErrName,objSystErrName)) 
    return 0;
  fout.Close();
  return 1;
}

// ----------
// ----------

// uses opened file, but only objName is provided

template<class Container_t, class objType_t>
int LoadToContainerFromFile(
   Container_t &store, objType_t &obj, int &hasSystErr,
   TFile &fin, const TString &objName, int loadSystErr) {
  TString systErr;
  if (loadSystErr) systErr=objName+TString("SystErr");
  return LoadToContainerFromFile(store,obj,hasSystErr,
			 fin, 
			 objName, objName+TString("Err"), systErr);
}

// ----------

// opens a file

template<class Container_t, class objType_t>
int LoadToContainer(
	 Container_t &store, objType_t &obj, int &hasSystErr,
	 const TString &fname, 
	 const TString &objName, int loadSystErr=0) {
  TFile fin(fname,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to load file <" << fname << ">\n";
    return 0;
  }
  if (!LoadToContainerFromFile(store,obj,hasSystErr,
			       fin,objName, loadSystErr)) return 0;
  fin.Close();
  return 1;
}


// =====================================
// =====================================

int AddInQuadrature(TVectorD &dest, const TVectorD &e1) {
  for (int i=0; i<e1.GetNoElements(); ++i) {
    double x=dest[i];
    dest[i] = sqrt(x*x + e1[i]*e1[i]);
  }
  return 1;
}

// ----------

int AddInQuadrature(TVectorD &dest, const TVectorD &e1, double c1) {
  double c1sqr= c1*c1;
  for (int i=0; i<e1.GetNoElements(); ++i) {
    double x=dest[i];
    dest[i] = sqrt(x*x + c1sqr *e1[i]*e1[i]);
  }
  return 1;
}

// ----------

int AddInQuadrature(TVectorD &dest, const TVectorD &e1, const TVectorD &e2) {
  for (int i=0; i<e1.GetNoElements(); ++i) {
    double x=dest[i];
    dest[i] = sqrt(x*x + e1[i]*e1[i] + e2[i]*e2[i]);
  }
  return 1;
}

// ----------

//
//  err(a b) = sqrt( (aErr b)^2 + (a bErr)^2 )
//

int SumInRelQuadrature(TVectorD &dest, const TVectorD &err1, const TVectorD &val1, const TVectorD &err2, const TVectorD &val2) {
  for (int i=0; i<err1.GetNoElements(); ++i) {
    dest[i] = sqrt( pow(err1[i]*val2[i],2) + pow(val1[i]*err2[i],2) );
  }
  return 1;
}

// =====================================
// =====================================

class VXSectD_t : public BaseClass_t {
protected:
  TVectorD xSec;
  TVectorD xSecErr;
  TVectorD xSecSystErr;
  TString FName;
  int FHasSystErr;
public:
  VXSectD_t(TString name, int dim=0) : 
    BaseClass_t("VXSectD_t"),
    xSec(dim), xSecErr(dim), xSecSystErr(dim),
    FName(name), FHasSystErr(0)
  {
    this->Zero();
  }

 // ----------

  VXSectD_t(TString name, const VXSectD_t &V) : 
    BaseClass_t("VXSectD_t"),
    xSec(V.xSec), xSecErr(V.xSecErr), xSecSystErr(V.xSecSystErr),
    FName(name), FHasSystErr(V.FHasSystErr)
  {}


 // ----------

  const TString& GetName() const { return FName; }
  void SetName(const TString newName) { FName=newName; }
  int GetDim() const { return xSec.GetNoElements(); }

  const TVectorD& GetCenter() const { return xSec; }
  const TVectorD& GetError() const { return xSecErr; }
  const TVectorD& GetSystError() const { return xSecSystErr; }
  double GetCenter(int i) const { return xSec[i]; }
  double GetError(int i) const { return xSecErr[i]; }
  double GetSystError(int i) const { return xSecSystErr[i]; }
  double GetTotError(int i) const { return sqrt(xSecErr[i]*xSecErr[i] + xSecSystErr[i]*xSecSystErr[i]); }
  TVectorD& EditCenter() { return xSec; }
  TVectorD& EditError() { return xSecErr; }
  TVectorD& EditSystError() { return xSecSystErr; }

  double operator[](int i) const { return xSec[i]; }

  const TVectorD& GetXS() const { return xSec; }
  TVectorD& EditXS() { return xSec; }
  const TVectorD& GetXSerr() const { return xSecErr; }
  TVectorD& EditXSerr() { return xSecErr; }
  const TVectorD& GetXSsystErr() const { return xSecSystErr; }
  TVectorD& EditXSsystErr() { return xSecSystErr; }

  int HasSystErr() const { return FHasSystErr; }
  void HasSystErr(int yes=1) { FHasSystErr=yes; }

  void SetDirectory(void*) const {}  // a dummy
  bool InheritsFrom(const TString &) const { return false; }

 // ----------

  int Assign(const VXSectD_t &v) {
    xSec=v.xSec;
    xSecErr=v.xSecErr;
    xSecSystErr=v.xSecSystErr;
    FHasSystErr=v.FHasSystErr;
    return 1;
  }

 // ----------
 // ----------

  void Zero() { xSec.Zero(); xSecErr.Zero(); xSecSystErr.Zero(); }

 // ----------

  void Scale(double x) {
    xSec *= x;
    xSecErr *= x;
    xSecSystErr *= x;
  }

 // ----------

  // a=x*y
  // (da/a)^2 = (dx/x)^2 + (dy/y)^2
  // (da)^2 = (y * dx)^2 + (a dy/y)^2
  //
  void Scale(double y, double yErr) {
    xSec *= y;
    xSecErr *= y;
    xSecSystErr *= y;
    TVectorD tmp(xSec); // a
    tmp *= (yErr/y);
    xSecSystErr += tmp;
  }

 // ----------

  int Scale(const TVectorD &vec) {
    for (int i=0; i<vec.GetNoElements(); ++i) {
      xSec[i] *= vec[i];
      xSecErr[i] *= vec[i];
      xSecSystErr[i] *= vec[i];
    }
    return 1;
  }

 // ----------

  int Divide(const TVectorD &vec) {
    for (int i=0; i<vec.GetNoElements(); ++i) {
      xSec[i] /= vec[i];
      xSecErr[i] /= vec[i];
      xSecSystErr[i] /= vec[i];
    }
    return 1;
  }

 // ----------

  // a=x/y
  // (da/a)^2 = (dx/x)^2 + (dy/y)^2
  // (da)^2 = (dx/y)^2 + (a dy/y)^2
  //
  void Divide(double y, double yErr) {
    xSec *= (1/y);
    xSecErr *= (1/y);
    xSecSystErr *= (1/y);
    TVectorD tmp(xSec); // a
    tmp *= (yErr/y);
    xSecSystErr += tmp;
  }

 // ----------

  int Divide(const VXSectD_t &V) {
    for (int i=0; i<V.xSec.GetNoElements(); ++i) {
      xSec[i] /= V.xSec[i];
    }
    for (int i=0; i<V.xSecErr.GetNoElements(); ++i) {
      xSecErr[i] /= V.xSec[i];
    }
    for (int i=0; i<V.xSecSystErr.GetNoElements(); ++i) {
      xSecSystErr[i] = sqrt( 
			    pow(xSecSystErr[i]/V.xSec[i],2) +
			    pow(xSec[i] * V.xSecErr[i]/V.xSec[i],2) 
			    );
    }
    return 1;
  }

 // ----------

  int Add(const VXSectD_t &V, double mult=1.) {
    for (int i=0; i<V.xSec.GetNoElements(); ++i) {
      xSec[i] += mult * V.xSec[i];
    }
    int res= 
      AddInQuadrature(xSecErr, V.xSecErr, mult) &&
      AddInQuadrature(xSecSystErr, V.xSecSystErr, mult);
    return res;
  }
 // ----------

  int AddExtraSystError(const TVectorD &V, int relativeErr) {
    int res=1;
    if (relativeErr) {
      TVectorD tmp(V);
      for (int i=0; i<tmp.GetNoElements(); ++i) tmp[i] *= xSec[i];
      res= AddInQuadrature(xSecSystErr, tmp, 1.);
    }
    else {
      res= AddInQuadrature(xSecSystErr, V, 1.);
    }
    if (!res) printError("in AddExtraSystError");
    return res;
  }


 // ----------

  int Accumulate(int idx, double weight) {
    if ((idx<0) || (idx>=xSecErr.GetNoElements())) return 0;
    xSec   [idx] += weight;
    xSecErr[idx] += weight*weight;
    return 1;
  }

 // ----------

  void Accumulate_end(int count) {
    double factor=1/double(count);
    for (int idx=0; idx<xSecErr.GetNoElements(); ++idx) {
      xSec[idx] *= factor;
      double w2= factor*xSecErr[idx];
      xSecErr[idx] = sqrt(w2 - pow(xSec[idx],2));
    }
    return;
  }

 // ----------

  void FinalizeSumErr() {
    for (int idx=0; idx<xSecErr.GetNoElements(); ++idx) {
      double w2=xSecErr[idx];
      xSecErr[idx] = sqrt(w2);
    }
    return;
  }

 // ----------

  int Load(TFile &fin, const TString &objName, int loadSystErr=0) {
    int res=LoadToContainerFromFile(*this,xSec,FHasSystErr, 
				    fin, objName, loadSystErr);
    if (!res) printError("Load");
    return res;
  }

  /*
  int Load(TFile &fin, const TString &objName, int loadSystErr=0) {
    if (!fin.IsOpen()) return 0;
    xSec = * (TVectorD*) fin.Get(objName);
    xSecErr = * (TVectorD*) fin.Get(objName + TString("Err"));
    if (loadSystErr) {
      xSecSystErr = * (TVectorD*) fin.Get(objName + TString("SystErr"));
    }
    else xSecSystErr.Zero();
    FHasSystErr=loadSystErr;
    HERE("HasSystErr=%d",FHasSystErr);
    return 1;
  }
  */

 // ----------

  int Load(TFile &fin, const TString &objName, const TString &objErrName, const TString &objSystErrName) {
    int res=LoadToContainerFromFile(*this,xSec,FHasSystErr,
				    fin, objName,objErrName,objSystErrName);
    if (!res) printError("Load");
    return res;
  }
  
  // ----------

  int Load(const TString &fname, const TString &objName, int loadSystErr=0) {
    int res=LoadToContainer(*this,xSec,FHasSystErr,
			    fname, objName,loadSystErr);
    if (!res) printError("Load(%s)",fname);
    return res;
  }
   
    /*
  int Load(const TString &fname, const TString &objName, int loadSystErr=0) {
    TFile fin(fname,"read");
    if (!fin.IsOpen()) {
      std::cout << "failed to load file <" << fname << ">\n";
      return 0;
    }
    if (!Load(fin,objName, loadSystErr)) return 0;
    fin.Close();
    return 1;
  }
    */
 // ----------

  int Load(const TString &fname, const TString &objName, const TString &objErrName, const TString &objSystErrName) {
    int res=LoadToContainer(*this,xSec,FHasSystErr,
			    fname, objName,objErrName,objSystErrName);
    if (!res) printError("Load");
    return res;
  }

  // ----------

  int Save(TFile &fout, const TString &objName, const TString &objErrName, const TString &objSystErrName) const {
    int res=SaveContainer(*this,xSec,
			  fout, objName,objErrName,objSystErrName);
    if (!res) printError("Save");
    return res;
  }

  // ----------

  int Read(TFile &fin) { 
    TString errName=FName+TString("Err");
    TString errSystName=FName+TString("SystErr");
    int res=LoadToContainerFromFile(*this,xSec,FHasSystErr,fin,
				    FName,errName,errSystName);
    if (!res) printError("Read(fin)");
    return res;
  }

  // ----------

  int Write(TFile &fout) const { 
    TString errName=FName+TString("Err");
    TString errSystName=FName+TString("SystErr");
    int res=SaveContainer(*this,xSec,fout,
			  FName,errName,errSystName);
    if (!res) printError("Write(fout)");
    return res;
  }

  // ----------

  int Read(TString saveName) { 
    TString errName=saveName+TString("Err");
    TString errSystName=saveName+TString("SystErr");
    int res=LoadToContainer_implicitFile(*this,xSec,FHasSystErr,
			  saveName,errName,errSystName);
    if (!res) printError("Read(%s)",saveName);
    return res;
  }

  // ----------

  int Write(TString saveName) const { 
    TString errName=saveName+TString("Err");
    TString errSystName=saveName+TString("SystErr");
    int res=SaveContainer_implicitFile(*this,xSec,
			  saveName,errName,errSystName);
    if (!res) printError("Write(%s)",saveName);
    return res;
  }

  // ----------
  // ----------

  int LoadTH1D(const TString &fname, const TString &hName, int loadError=1) {
    this->Zero();
    TFile fin(fname,"read");
    if (!fin.IsOpen()) {
      return this->reportError("LoadTH1D: failed to open a file <%s>",fname);
    }
    TH1D *h=(TH1D*) fin.Get(hName);
    if (!h) {
      fin.Close();
      return this->reportError("LoadTH1D: failed to get <%s> from <%s>",hName,fname);
    }
    for (int ir=0; ir<h->GetNbinsX(); ++ir) {
      xSec(ir) = h->GetBinContent(ir+1);
      if (loadError) xSecErr(ir)= h->GetBinError(ir+1);
    }
    delete h;
    fin.Close();
    return 1;
  }

  // ----------
  // ----------

  // sqrt(statErr^2 + systErr^2)

  int GetTotalError(TVectorD &totErr) const {
    for (int i=0; i<xSecErr.GetNoElements(); i++) {
      totErr(i) = sqrt( xSecErr[i]*xSecErr[i] + 
			xSecSystErr[i]*xSecSystErr[i] );
    }
    return 1;
  }

  // ----------

  // statErr[i]/xsect[i]

  int GetRelativeStatError(TVectorD &relErr) const {
    for (int i=0; i<xSecErr.GetNoElements(); i++) {
      relErr(i) = xSecErr[i] / xSec[i];
    }
    return 1;
  }

  // ----------

  // totErr[i]/xsect[i]

  int GetRelativeTotError(TVectorD &relErr) const {
    for (int i=0; i<xSecErr.GetNoElements(); i++) {
      double totErr = sqrt( xSecErr[i]*xSecErr[i] + 
			    xSecSystErr[i]*xSecSystErr[i] );
      relErr(i)=totErr/xSec[i];
    }
    return 1;
  }

  // ----------

  int GetRelativeTotalError(TVectorD &relErr) const { return GetRelativeTotError(relErr); }

  // ----------

  /*
  int Unfold(const VXSectD_t &srcVec, const TString &fnameUnf, 
	     int unfoldXS=1, int unfoldXSecErr=1, int unfoldXSecSystErr=1) {
#ifndef MYTOOLS_HH
    using namespace unfolding;
#endif
    xSec.Zero(); xSecErr.Zero(); xSecSystErr.Zero();
    int res=1;
    if (unfoldXS) {
      res= unfold(srcVec.xSec, xSec, fnameUnf);
    }
    if (res && unfoldXSecErr) {
      res= propagateErrorThroughUnfolding(srcVec.xSecErr, xSecErr, fnameUnf);
    }
    if (res && unfoldXSecSystErr) {
      res= propagateErrorThroughUnfolding(srcVec.xSecSystErr, xSecSystErr, fnameUnf);
    }
    if (!res) printError("Unfold");
    return res;
  }

  // ----------

  int UnfoldAndPropagate(const VXSectD_t &srcVec,
			 const TMatrixD &U, int transpose) {
    this->Zero();
    TMatrixD Utrans(TMatrixD::kTransposed,U);
    const TMatrixD *M=(transpose) ? &Utrans : &U;
    xSec = (*M) * srcVec.xSec;
    for (int i=0; i<srcVec.xSecErr.GetNoElements(); ++i) {
      double sum=0;
      for (int j=0; j< M->GetNcols(); ++j) {
	sum += pow ( (*M)(i,j) * srcVec.xSecErr[j], 2);
      }
      xSecErr(i) = sqrt(sum);
    }
    if (srcVec.FHasSystErr) {
      FHasSystErr=1;
      for (int i=0; i<srcVec.xSecSystErr.GetNoElements(); ++i) {
	double sum=0;
	for (int j=0; j< M->GetNcols(); ++j) {
	  sum += pow ( (*M)(i,j) * srcVec.xSecSystErr[j], 2);
	}
	xSecSystErr(i) = sqrt(sum);
      }
    }
    return 1;
  }
  */

  // ----------

  template<class Histo1D_t>
  int FillHisto(Histo1D_t *h, int noError=0) const {
    for (int i=0; i<xSec.GetNoElements(); ++i) {
      h->SetBinContent(i+1,xSec[i]);
      if (noError) h->SetBinError(i+1,0.); else h->SetBinError(i+1,xSecErr[i]);
    }
    return 1;
  }

  // ----------

  template<class Histo1D_t>
  int FillHistoWithError(Histo1D_t *h) const {
    for (int i=0; i<xSec.GetNoElements(); ++i) {
      h->SetBinContent(i+1,xSecErr[i]);
      h->SetBinError(i+1,0.);
    }
    return 1;
  }

  // ----------

  template<class Histo1D_t>
  int FillHistoWithRelError(Histo1D_t *h) const {
    for (int i=0; i<xSec.GetNoElements(); ++i) {
      h->SetBinContent(i+1,xSecErr[i]/xSec[i]);
      h->SetBinError(i+1,0.);
    }
    return 1;
  }

  // ----------

  template<class Histo1D_t>
  int FillHistoWithTotError(Histo1D_t *h) const {
    for (int i=0; i<xSec.GetNoElements(); ++i) {
      double totErr=sqrt(xSecErr[i]*xSecErr[i] + xSecSystErr[i]*xSecSystErr[i]);
      h->SetBinContent(i+1,totErr);
      h->SetBinError(i+1,0.);
    }
    return 1;
  }

  // ----------

  template<class Histo1D_t>
  int FillHistoWithTotRelError(Histo1D_t *h) const {
    for (int i=0; i<xSec.GetNoElements(); ++i) {
      double totErr=sqrt(xSecErr[i]*xSecErr[i] + xSecSystErr[i]*xSecSystErr[i]);
      h->SetBinContent(i+1,totErr/xSec[i]);
      h->SetBinError(i+1,0.);
    }
    return 1;
  }

  // ----------

  // useError=1 (stat) or 2 (syst) or 4 (total)

  int SaveFig(const TString &fname, int useMassBins, int useError, int do_delete=1) const {
    TString hName=TString("histo")+FName;
    TH1F *h=(useMassBins) ? 
      new TH1F(hName,"",DYTools::nMassBins,DYTools::massBinLimits) :
      new TH1F(hName,"",xSec.GetNoElements(),0,xSec.GetNoElements());
    h->SetDirectory(0);
    double err=0;
    for (int i=0; i<xSec.GetNoElements(); ++i) {
      h->SetBinContent(i+1, xSec[i]);
      err=0;
      switch (useError) {
      case 0: break;
      case 1: err=xSecErr[i]; break;
      case 2: err=xSecSystErr[i]; break;
      case 4: err=sqrt(pow(xSecErr[i],2)+pow(xSecSystErr[i],2)); break;
      default: reportError("SaveFig: useError cannot be %d",useError);
      }
      h->SetBinError(i+1, err);
    }
    TString canvName=TString("canv_") + FName;
    TCanvas *cx=new TCanvas(canvName,canvName,600,600);
    if (useMassBins) cx->SetLogx();
    h->Draw((useError) ? "LPE" : "LP");
    cx->SaveAs(fname);
    if (do_delete) {
      delete cx;
      delete h;
    }
    return 1;
  }

  // ----------

  void Print(const char *xsFormat="%6.4lf") const { 
    char format[100];
    if (FHasSystErr) sprintf(format,"%s  %s  %s  %s\n","%2d",xsFormat,xsFormat,xsFormat);
    else sprintf(format,"%s  %s  %s\n","%2d",xsFormat,xsFormat);
    std::cout << "VXSectD_t(" << FName << "):\n";
    HERE("HasSystErr=%d",FHasSystErr);
    for (int i=0; i<xSec.GetNoElements(); ++i) {
	printf(format,i,xSec[i],xSecErr[i],xSecSystErr[i]);
    }
    std::cout << std::endl;
  }


  // ----------

};

// =====================================

class MXSectD_t : public BaseClass_t {
protected:
  TMatrixD xSec;
  TMatrixD xSecErr;
  TMatrixD xSecSystErr;
  TString FName;
  int FHasSystErr;
public:
  MXSectD_t(const TString &name, int nRows, int nCols) : 
    BaseClass_t("MXSectD_t"),
    xSec(nRows,nCols), 
    xSecErr(nRows,nCols), 
    xSecSystErr(nRows,nCols),
    FName(name),
    FHasSystErr(0)
  {
    this->Zero();
  }

 // ----------

  MXSectD_t(const TString &name, const MXSectD_t &M) :
    BaseClass_t("MXSectD_t"),
    xSec(M.xSec),
    xSecErr(M.xSecErr),
    xSecSystErr(M.xSecSystErr),
    FName(name),
    FHasSystErr(M.FHasSystErr) 
  { }

 // ----------

  MXSectD_t(const TString &name, TMatrixD::EMatrixCreatorsOp1 op, const MXSectD_t &M) :
    BaseClass_t("MXSectD_t"),
    xSec(op,M.xSec),
    xSecErr(op,M.xSecErr),
    xSecSystErr(op,M.xSecSystErr),
    FName(name),
    FHasSystErr(M.FHasSystErr) 
  { }

 // ----------

  const TString& GetName() const { return FName; }
  void SetName(const TString &newName) { FName=newName; }

  const TMatrixD& GetCenter() const { return xSec; }
  const TMatrixD& GetError() const { return xSecErr; }
  const TMatrixD& GetSystError() const { return xSecSystErr; }
  TMatrixD &EditCenter() { return xSec; }
  TMatrixD &EditError() { return xSecErr; }
  TMatrixD &EditSystError() { return xSecSystErr; }

  const TMatrixD& GetXS() const { return xSec; }
  TMatrixD& EditXS() { return xSec; }
  const TMatrixD& GetXSerr() const { return xSecErr; }
  TMatrixD& EditXSerr() { return xSecErr; }
  const TMatrixD& GetXSsystErr() const { return xSecSystErr; }
  TMatrixD& EditXSsystErr() { return xSecSystErr; }

  int HasSystErr() const { return FHasSystErr; }
  void HasSystErr(int yes=1) { FHasSystErr=yes; }

  void SetDirectory(void*) const {}  // a dummy
  bool InheritsFrom(const TString &) const { return false; }

 // ----------

  int Assign(const MXSectD_t &m) {
    xSec=m.xSec;
    xSecErr=m.xSecErr;
    xSecSystErr=m.xSecSystErr;
    FHasSystErr=m.FHasSystErr;
    return 1;
  }

 // ----------

  void Zero() {   xSec.Zero();   xSecErr.Zero();   xSecSystErr.Zero();   }

 // ----------

  void Transpose() {
    TMatrixD tmp(xSec);
    xSec.Transpose(tmp);
    tmp=xSecErr;   xSecErr.Transpose(tmp);
    tmp=xSecSystErr;    xSecSystErr.Transpose(tmp);
  }

 // ----------

  int Accumulate(int ir, int ic, double weight) {
    if ((ir<0) || (ir>=xSecErr.GetNrows())) return 0;
    if ((ic<0) || (ic>=xSecErr.GetNcols())) return 0;
    xSec   (ir,ic) += weight;
    xSecErr(ir,ic) += weight*weight;
    return 1;
  }

 // ----------

  int Fill(double mass, double rapidity, double weight) {
    int iM=DYTools::findMassBin(mass);
    int iY=DYTools::findAbsYBin(iM, rapidity);
    if ((iM!=-1) && (iY!=-1)) {
      this->Accumulate(iM,iY,weight);
      return 1;
    }
    return 0;
  }

 // ----------

  void FinalizeSumErr() {
    for (int ir=0; ir<xSecErr.GetNrows(); ++ir) {
      for (int ic=0; ic<xSecErr.GetNcols(); ++ic) {
	double w2=xSecErr(ir,ic);
	xSecErr(ir,ic) = sqrt(w2);
      }
    }
    return;
  }
 // ----------

  // 
  // special use
  //

  int AccumulateForAvg(int ir, int ic, double quantity, double weight) {
    if ((ir<0) || (ir>=xSecErr.GetNrows())) return 0;
    if ((ic<0) || (ic>=xSecErr.GetNcols())) return 0;
    xSec   (ir,ic) += weight;
    xSecErr(ir,ic) += (weight*quantity);
    xSecSystErr(ir,ic) += quantity*(weight*quantity);
    return 1;
  }

 // ----------

  //
  // special use
  //

  int CalcAvgAndRMS(MXSectD_t &res, double badValue=-9999) {
    res.Zero();
    int no_bad_value=1;
    for (int ir=0; ir<xSecErr.GetNrows(); ++ir) {
      for (int ic=0; ic<xSecErr.GetNcols(); ++ic) {
	double w= this->xSec(ir,ic);
	if (w!=0.) {
	  double mean= this->xSecErr(ir,ic)/w;
	  res.xSec(ir,ic) = mean;
	  double errSqr=this->xSecSystErr(ir,ic)/w - mean*mean;
	  if (errSqr<0) std::cout << "negative value!\n";
	  res.xSecErr(ir,ic) = sqrt(fabs(errSqr));
	}
	else {
	  no_bad_value=0;
	  res.xSec(ir,ic) = badValue;
	  res.xSecErr(ir,ic) = badValue;
	}
      }
    }
    return no_bad_value;
  }


 // ----------

  int Load(TFile &fin, const TString &objName, int loadSystErr=0) {
    int res=LoadToContainerFromFile(*this,xSec,FHasSystErr, fin,objName,loadSystErr);
    if (!res) printError("Load(TFile)");
    return res;
  }

  /*
  int Load(TFile &fin, const TString &objName, int loadSystErr=0) {
    if (!fin.IsOpen()) return 0;
    xSec = * (TMatrixD*) fin.Get(objName);
    xSecErr = * (TMatrixD*) fin.Get(objName + TString("Err"));
    if (loadSystErr) {
      xSecSystErr = * (TMatrixD*) fin.Get(objName + TString("SystErr"));
    }
    else xSecSystErr.Zero();
    FHasSystErr=loadSystErr;
    return 1;
  }
  */

 // ----------

  int Load(TFile &fin, const TString &objName, const TString &objErrName, const TString &objSystErrName) {
    int res=LoadToContainerFromFile(*this,xSec,FHasSystErr,
			    fin, objName,objErrName,objSystErrName);
    if (!res) printError("Load(TFile)");
    return res;
  }

  // ----------
  
  int Load(const TString &fname, const TString &objName, int loadSystErr=0) {
    HERE("MXSectD_t::Load  load obj=<%s> from <%s>",objName,fname);
    int res=LoadToContainer(*this,xSec,FHasSystErr, 
			    fname, objName,loadSystErr);
    if (!res) printError("Load(TFile)");
    return res;
  }

  /*
  int Load(const TString &fname, const TString &objName, int loadSystErr=0) {
    TFile fin(fname,"read");
    if (!fin.IsOpen()) {
      return this->reportError("failed to load file <%s>\n",objName.Data());
    }
    if (!Load(fin,objName, loadSystErr)) return 0;
    fin.Close();
    return 1;
  }
  */

 // ----------

  int Load(const TString &fname, const TString &objName, const TString &objErrName, const TString &objSystErrName) {
    int res=LoadToContainer(*this,xSec,FHasSystErr,
			    fname, objName,objErrName,objSystErrName);
    if (!res) printError("Load(fname)");
    return res;
  }

  // ----------

  int Save(TFile &fout, const TString &objName, const TString &objErrName, const TString &objSystErrName) const {
    int res=SaveContainer(*this,xSec,
			  fout, objName,objErrName,objSystErrName);
    if (!res) printError("Save(TFile)");
    return res;
  }

  // ----------

  int Save(const TString &fname, const TString &objName, const TString &objErrName, const TString &objSystErrName) const {
    int res=SaveContainer(*this,xSec,
			  fname, objName,objErrName,objSystErrName);
    if (!res) printError("Save(fname)");
    return res;
  }

  // ----------

  int Read(TFile &fin) { 
    TString errName=FName+TString("Err");
    TString errSystName=FName+TString("SystErr");
    int res=LoadToContainerFromFile(*this,xSec,FHasSystErr,fin,
				    FName,errName,errSystName);
    if (!res) printError("Read(fin)");
    return res;
  }

  // ----------

  int Write(TFile &fout) const { 
    TString errName=FName+TString("Err");
    TString errSystName=FName+TString("SystErr");
    int res=SaveContainer(*this,xSec,fout,
			  FName,errName,errSystName);
    if (!res) printError("Write(fout)");
    return res;
  }

  // ----------

  int Read(TString saveName) { 
    TString errName=saveName+TString("Err");
    TString errSystName=saveName+TString("SystErr");
    int res=LoadToContainer_implicitFile(*this,xSec,FHasSystErr,
				    saveName,errName,errSystName);
    if (!res) printError("Read(%s)",saveName);
    return res;
  }

  // ----------

  int Write(TString saveName) const { 
    TString errName=saveName+TString("Err");
    TString errSystName=saveName+TString("SystErr");
    int res=SaveContainer_implicitFile(*this,xSec,
			  saveName,errName,errSystName);
    if (!res) printError("Write(%s)",saveName);
    return res;
  }

  // ----------
  // ----------

  int LoadTH2D(const TString &fname, const TString &objName, int loadError=1) {
    this->Zero();
    TFile fin(fname,"read");
    if (!fin.IsOpen()) {
      return this->reportError("LoadTH2D: failed to open a file <%s>",fname);
    }
    TH2D *h=(TH2D*) fin.Get(objName);
    if (!h) {
      fin.Close();
      return this->reportError("LoadTH2D: failed to get <%s> from <%s>",objName,fname);
    }
    for (int ir=0; ir<h->GetNbinsX(); ++ir) {
      for (int ic=0; ic<h->GetNbinsY(); ++ic) {
	xSec(ir,ic) = h->GetBinContent(ir+1,ic+1);
	if (loadError) xSecErr(ir,ic)= h->GetBinError(ir+1,ic+1);
      }
    }
    delete h;
    fin.Close();
    return 1;
  }

  // ----------
  // ----------

  int FlattenTo(VXSectD_t &vec) const {
#ifndef MYTOOLS_HH
    using namespace unfolding;
#endif
    vec.HasSystErr( this->FHasSystErr );
    int res=
      flattenMatrix(this->xSec, vec.EditXS()) &&
      flattenMatrix(this->xSecErr, vec.EditXSerr());
    if (res && this->FHasSystErr) {
      res=flattenMatrix(this->xSecSystErr, vec.EditXSsystErr());
    }
    if (!res) {
      this->printError("FlattenTo");
    }
    return res;
  }

  // ----------

  int DeflattenFrom(const VXSectD_t &vec) {
#ifndef MYTOOLS_HH
    using namespace unfolding;
#endif
    this->FHasSystErr = vec.HasSystErr();
    int res=
      deflattenMatrix(vec.GetXS(), this->xSec) &&
      deflattenMatrix(vec.GetXSerr(), this->xSecErr);
    if (res && vec.HasSystErr()) {
      res=deflattenMatrix(vec.GetXSsystErr(),this->xSecSystErr);
    }
    if (!res) {
      this->printError("DeflattenFrom");
    }
    return res;
  }

  // ----------
/*
  int Unfold(const MXSectD_t &yields, const TString &unfConstFName,
	     int unfoldMatrix, int propagateStatErr, int propagateSystErr) {
#ifndef MYTOOLS_HH
    using namespace unfolding;
#endif
    int dim= xSec.GetNrows() * xSec.GetNcols();
    TVectorD vim(dim), vout(dim);
    int res=1;
    FHasSystErr= (yields.FHasSystErr && propagateSystErr) ? 1:0;
    if (unfoldMatrix) { 
      res=unfold(yields.xSec, this->xSec, unfConstFName,
			    vim,vout);
    }
    if (propagateStatErr) {
      res=propagateErrorThroughUnfolding(
		    yields.xSecErr, this->xSecErr,
		    unfConstFName,
		    vim,vout);
    }
    else {
      this->xSecErr.Zero();
    }
    if (propagateSystErr) {
      if (!yields.FHasSystErr) {
	std::cout << "warning: Unfold. propagateSystErr=1, but yields has not syst err\n";
      }
      else {
	res=propagateErrorThroughUnfolding(
		    yields.xSecSystErr, this->xSecSystErr,
		    unfConstFName,
		    vim,vout);
      }
    }
    else {
      this->xSecSystErr.Zero();
    }
    if (!res) this->printError("Unfold(%d,%d,%d)\n",unfoldMatrix,propagateStatErr,propagateSystErr);
    return res;
  }

  // ----------

  int UnfoldAndPropagate(const MXSectD_t &srcM,
			 const TMatrixD &U, int transpose) {
    this->Zero();
    int totDim=U.GetNrows();
    VXSectD_t iniV("iniV",totDim);
    VXSectD_t finV("finV",totDim);
    int res= (srcM.FlattenTo(iniV) &&
	      finV.UnfoldAndPropagate(iniV,U,transpose) &&
	      this->DeflattenFrom(finV)) ? 1:0;
    if (!res) this->reportError("UnfoldAndPropagate");
    return res;
  }

  // ----------

  int UnfoldAndPropagate(const MXSectD_t &srcM,
			 const MXSectD_t &U, int transpose) {
    int res=this->UnfoldAndPropagate(srcM,U.xSec,transpose);
    if (!res) reportError("in UnfoldAndPropagate(MXSectD{%s},MXSectD{%s}",srcM.FName,U.FName);
    return res;
  }
*/

  // ----------

  template<class Histo1D_t>
  int FillHisto(Histo1D_t *h, int ic, int noError=0) const {
    for (int i=0; i<xSec.GetNrows(); ++i) {
      h->SetBinContent(i+1,xSec(i,ic));
      double err=(noError) ? 0. : xSecErr(i,ic);
      h->SetBinError(i+1,err);
    }
    return 1;
  }

  // ----------

  template<class Histo1D_t>
  int FillHistoWithError(Histo1D_t *h, int ic, int dontUseSystErrAsError) const {
    for (int i=0; i<xSecErr.GetNrows(); ++i) {
      h->SetBinContent(i+1,xSecErr(i,ic));
      double err=(dontUseSystErrAsError) ? 0. : xSecSystErr(i,ic);
      h->SetBinError(i+1,err);
    }
    return 1;
  }

  // ----------

  template<class Histo2D_t>
    int FillHisto2D(Histo2D_t *h) const {
    for (int ir=0; ir<xSec.GetNrows(); ++ir) {
      for (int ic=0; ic<xSec.GetNcols(); ++ic) {
	h->SetBinContent(ir+1,ic+1,xSec(ir,ic));
	h->SetBinError(ir+1,ic+1,xSecErr(ir,ic));
      }
    }
    return 1;
  }

  // ----------

};

// =====================================
// =====================================

int LoadAndFlatten(VXSectD_t &xSec,
		   const TString &fname, int xBins, int yBins,
		   const TString &xsName, const TString &xsErrName, 
		   const TString &xsSystErrName) {
  if (xBins*yBins!=xSec.GetDim()) {
    return xSec.reportError("LoadAndFlatten: vec[%d], while xBins=%d, yBins=%d",
			    xSec.GetDim(),xBins,yBins);
  }
  TString tmpName=TString("tmpM_") + xSec.GetName();
  MXSectD_t M(tmpName,xBins,yBins);
  int res=
    (M.Load(fname,xsName,xsErrName,xsSystErrName) &&
     M.FlattenTo(xSec)) ? 1:0;
  if (!res) return xSec.reportError("in LoadAndFlatten");
  return res;
}


// =====================================

#endif
