//
// This file contains implementation of methods related to applying
// electron energy scale corrections.
//

#include "../Include/ElectronEnergyScale.hh"
#include <fstream>
#include <sstream>
#include <algorithm>
#include "MyTools.hh"

#ifdef UseEEM
#include <TTree.h>
#include <TBranch.h>
#endif

#ifdef ZEE_DATA_HH
#include "../Include/PUReweight.hh"
#endif


//------------------------------------------------------

template<class T>
inline void SWAP(T &a, T &b) { T c=a; a=b; b=c; }

template<class T>
inline void SWAP(T &a, T &b,  T &tmp) { tmp=a; a=b; b=tmp; }

//------------------------------------------------------
//------------------------------------------------------

#ifdef ZEE_DATA_HH
template<class Type1_t, class Type2_t>
inline
int QuickSort(std::vector<Type1_t> &arr1, std::vector<Type2_t> &arr2, int calling=0, int print_err=1)
{
  // adapted from Numerical Recipes in C

  unsigned int count=arr1.size();
  if (arr1.size() != arr2.size()) {
    std::cout << "QuickSort: arr1.size=" << arr1.size() << ", arr2.size=" << arr2.size() << "\n";
    return 0;
  }

  const int insertion_sort_limit= 7;
  const int aux_index_size=100;
  if (count==0) return 1;

  long i,ir=count-1, j,k,l=0;
  std::vector<unsigned int> istack(aux_index_size);
  int jstack=0;
  Type1_t a1,tmp1;
  Type2_t a2,tmp2;
  
  for (;;) {
    if (ir-l<insertion_sort_limit) { // insertion sort when subarray is small enough
      for (j=l+1; j<=ir; ++j) {
	a1=arr1[j]; a2=arr2[j];
	for (i=j-1; i>=0; --i) {
	  if (a1>=arr1[i]) break;
	  arr1[i+1]=arr1[i];
	  arr2[i+1]=arr2[i];
	}
	arr1[i+1]=a1; arr2[i+1]=a2;
      }
      if (jstack<1) return 1;
      ir=istack[jstack-1]; // pop stack and begin a new round of partitioning
      l=istack[jstack-2];
      jstack -= 2;
    }
    else {
      // choose median of left, center and right elements as partitioning element a. Also rearrange so that a[l+1]<=a[l]<=a[ir]      
      k=(l+ir)>>1;
      SWAP(arr1[k],arr1[l+1],tmp1);
      SWAP(arr2[k],arr2[l+1],tmp2);
      if (arr1[l+1]>arr1[ir]) {
	SWAP(arr1[l+1],arr1[ir],tmp1);
	SWAP(arr2[l+1],arr2[ir],tmp2);
      }
      if (arr1[l]>arr1[ir]) { 
	SWAP(arr1[l],arr1[ir],tmp1);
	SWAP(arr2[l],arr2[ir],tmp2);
      }
      if (arr1[l+1]>arr1[l]) {
	SWAP(arr1[l+1],arr1[l],tmp1);
	SWAP(arr2[l+1],arr2[l],tmp2);
      }
      // initialize pointers of partitioning
      i=l+1;
      j=ir;
      // partitioning element
      a1=arr1[l]; a2=arr2[l];
      // innermost loop
      for (;;) {
	do ++i; while (a1>arr1[i]); // scan up to find element >a1
	do --j; while (arr1[j]>a1); // scan down to find element <a1
	if (j<i) break; // pointers crossed. Partitioning complete
	SWAP(arr1[i],arr1[j],tmp1);
	SWAP(arr2[i],arr2[j],tmp2);
      }
      arr1[l]=arr1[j]; arr1[j]=a1;
      arr2[l]=arr2[j]; arr2[j]=a2;
      jstack+=2;
      // push pointers to larger subarray on stack, process smaller subarray immediately.
      if (jstack>aux_index_size) {
	if (calling<5) return QuickSort(arr1,arr2,calling+1,print_err);
	if (print_err) printf("QuickSort: aux_index_size is too small in QSort(cnt=%d,type*,type*)",count);
	//return (print_err) ? WAM_infunction("QuickSort") : 0;
	return 0;
      }
      if (ir-l+1>=j-l) {
	istack[jstack-1]=ir;
	istack[jstack-2]=i;
	ir=j-1;
      }
      else {
	istack[jstack-1]=j-1;
	istack[jstack-2]=l;
	l=i;
      }
    }
  }
  return 1;
}

#endif


//------------------------------------------------------
//------------------------------------------------------

ElectronEnergyScale::ElectronEnergyScale(CalibrationSet calibrationSet):
  _calibrationSet(calibrationSet),
  _inpFileName(),
  _isInitialized(false),
  _randomizedStudy(false),
  _nEtaBins(0),
  _etaBinLimits(NULL),
  _dataConst(NULL), _dataConstErr(NULL), _dataConstRandomized(NULL),
  _energyScaleCorrectionRandomizationDone(false),
  _dataSeed(-1),
  _nMCConstants(0),
  _mcConst1Name(), _mcConst2Name(), _mcConst3Name(), _mcConst4Name(),
  _mcConst1(NULL), _mcConst2(NULL), _mcConst3(NULL), _mcConst4(NULL),
  _mcConst1Err(NULL), _mcConst2Err(NULL), _mcConst3Err(NULL), _mcConst4Err(NULL),
  _smearingWidthRandomizationDone(false),
  _mcSeed(-1)
{
  this->init(calibrationSet);
}

//------------------------------------------------------

ElectronEnergyScale::ElectronEnergyScale(const TString &escaleTagName):
  _calibrationSet(UNDEFINED),
  _inpFileName(),
  _isInitialized(false),
  _randomizedStudy(false),
  _nEtaBins(0),
  _etaBinLimits(NULL),
  _dataConst(NULL), _dataConstErr(NULL), _dataConstRandomized(NULL),
  _energyScaleCorrectionRandomizationDone(false),
  _dataSeed(-1),
  _nMCConstants(0),
  _mcConst1Name(), _mcConst2Name(), _mcConst3Name(), _mcConst4Name(),
  _mcConst1(NULL), _mcConst2(NULL), _mcConst3(NULL), _mcConst4(NULL),
  _mcConst1Err(NULL), _mcConst2Err(NULL), _mcConst3Err(NULL), _mcConst4Err(NULL),
  _smearingWidthRandomizationDone(false),
  _mcSeed(-1)
{
  this->init(escaleTagName);
}

//------------------------------------------------------

void ElectronEnergyScale::clear() {
  if (_isInitialized) {
    if (_etaBinLimits) delete[] _etaBinLimits;
    if (_dataConst) delete[] _dataConst;
    if (_dataConstErr) delete[] _dataConstErr;
    if (_dataConstRandomized) delete[] _dataConstRandomized;
    _dataSeed=-1;
    _mcSeed=-1;
    if (_mcConst1) delete[] _mcConst1;
    if (_mcConst2) delete[] _mcConst2;
    if (_mcConst3) delete[] _mcConst3;
    if (_mcConst4) delete[] _mcConst4;
    if (_mcConst1Err) delete[] _mcConst1Err;
    if (_mcConst2Err) delete[] _mcConst2Err;
    if (_mcConst3Err) delete[] _mcConst3Err;
    if (_mcConst4Err) delete[] _mcConst4Err;
  }
}

//------------------------------------------------------

void ElectronEnergyScale::init(const TString &stringWithEScaleTagName, int debug) {
  this->clear();
  //std::cout << "stringWithEScaleTagName=<" << stringWithEScaleTagName << ">\n";
  this->init(ElectronEnergyScale::DetermineCalibrationSet(stringWithEScaleTagName,&_inpFileName),debug);
  if (this->isInitialized()) {
    Ssiz_t pos=stringWithEScaleTagName.Index("RANDOMIZE");
    _randomizedStudy=(pos>=0);
    if (_randomizedStudy) {
      int shift=stringWithEScaleTagName.Contains("RANDOMIZED") ? 10:9;
      //std::cout << "seed from <" << (stringWithEScaleTagName.Data()+pos+shift) << ">\n";
      int seed=atoi(stringWithEScaleTagName.Data()+pos+shift);
      std::cout << "randomization with seed=" << seed << "\n";
      int res=this->randomizeEnergyScaleCorrections(seed);
      if (res) {
	if (stringWithEScaleTagName.Index("INVERT_RAND")>=0) {
	  std::cout << "inverting scale factors\n";
	  res=this->invertRandomizedEnergyScaleCorrections();
	}
	else if (stringWithEScaleTagName.Index("MIRROR_RAND")>=0) {
	  std::cout << "inverting scale factors\n";
	  res=this->mirrorRandomizedEnergyScaleCorrections();
	}
      }
      if (!res) _isInitialized=false;
    }
  }
  return;
}

//------------------------------------------------------

void ElectronEnergyScale::init(CalibrationSet calSet, int debug) {
  _calibrationSet=calSet;
  _isInitialized = false;
   
  if( _calibrationSet == UNDEFINED )
    return;

  std::cout << "_calibrationSet=" << ElectronEnergyScale::CalibrationSetName(_calibrationSet,&_inpFileName) << "\n";

  if( !initializeAllConstants(debug)) {
    std::cout << "failed to initialize\n";
    return;
  }

  if( !initializeExtraSmearingFunction()) {
    std::cout << "failed to prepare extra smearing function\n";
    return;
  }

  _isInitialized = true;
  return;
}

//------------------------------------------------------

bool ElectronEnergyScale::loadInputFile(const TString &fileName, int debug) {
  ifstream fin;
  fin.open(fileName);
  if (!fin) {
    std::cout << "failed to open a file <" << fileName << ">\n";
    return kFALSE;
  }
  std::vector<std::string> lines;
  std::string s;
  while (getline(fin,s)) {
    lines.push_back(s);
  }
  fin.close();
  if (debug) { std::cout << "from <" << fileName << "> loaded "; PrintVec("loaded ",lines,1); }
  // check that the file contains the model we expect
  int ok=0;
  std::string modelFile;
  for (unsigned int i=0; i<lines.size(); ++i) {
    size_t pos=lines[i].find("fit model");
    if (pos != std::string::npos) {
      size_t pos1=lines[i].find(')',pos);
      if (pos1 == std::string::npos) {
	std::cout << "cannot verify model defined in a file: not the right structure\n";
	throw 2;
      }
      std::string model=lines[i].substr(pos+10,pos1-pos-10);
      std::transform(model.begin(), model.end(), model.begin(),
               (int(*)(int)) std::tolower);
      if ((model.find("gauss")!=std::string::npos) && 
	  (_calibrationSet==CalSet_File_Gauss)) ok=1;
      else if ((model.find("voigt")!=std::string::npos) &&
	       (_calibrationSet==CalSet_File_Voigt)) ok=1;
      else if ((model.find("breitwigner")!=std::string::npos) &&
	       (_calibrationSet==CalSet_File_BreitWigner)) ok=1;
      modelFile=model;
      break;
    }
  }
  if (!ok) {
    std::cout << "The object is readied for <" << this->calibrationSetName() << ">, while the file contains model=<" << modelFile << ">\n";
    throw 2;
  }
  // assign constants
  bool res=assignConstants(lines,debug);
  if (!res) std::cout << "error from loadInputFile(" << fileName << ")\n";
  return res;
}

//------------------------------------------------------

bool ElectronEnergyScale::assignConstants(const std::vector<std::string> &lines, int debug) {
  clear();
  int linesHaveNegativeEtas=-1;
  int etaDivCount=0;
  for (unsigned int i=0; i<lines.size(); i++) {
    if ((lines[i].find("g_esfWorkCase")!=std::string::npos) && 
	(lines[i].find("g_esfWorkCaseShortName")==std::string::npos)) {
      if ((lines[i].find("bins on each eta side")!=std::string::npos) ||
	  (lines[i].find("bins for each eta side")!=std::string::npos)) {
	linesHaveNegativeEtas= 1;
      }
      else {
	linesHaveNegativeEtas= 0;
      }
    }
    else if (lines[i].find("EtaDivisionCount")!=std::string::npos) {
      etaDivCount=atoi(lines[i].c_str() + lines[i].find('=') + 1);
    }
  }
  // check that lines have the expected info
  if (debug) std::cout << "linesHaveNegativeEtas=" << linesHaveNegativeEtas << ", etaDivCount=" << etaDivCount << "\n"; 
  if ((linesHaveNegativeEtas==-1) || (etaDivCount==0)) {
    std::cout << "assignConstants: the provided lines do not have a keyword 'g_esfWorkCase' or 'EtaDivisionCount'\n";
    return kFALSE;
  }
  // allocate memory
  if (linesHaveNegativeEtas==0) etaDivCount*=2;
  if (debug) std::cout << "etaDivCount=" << etaDivCount << "\n";
  _nEtaBins = etaDivCount;
  _etaBinLimits = new double[_nEtaBins+1];
  _dataConst = new double[_nEtaBins];
  _dataConstErr = new double[_nEtaBins];
  _dataConstRandomized = new double[_nEtaBins];
  _nMCConstants=1;
  _mcConst1 = new double[_nEtaBins];
  _mcConst2 = 0;
  _mcConst3 = 0;
  _mcConst4 = 0;
  _mcConst1Err = new double[_nEtaBins];
  _mcConst2Err = 0;
  _mcConst3Err = 0;
  _mcConst4Err = 0;
  _mcConst1Name = "smear";

  assert(_etaBinLimits); 
  assert(_dataConst); assert(_dataConstErr);
  assert(_dataConstRandomized);
  assert(_mcConst1); assert(_mcConst1Err);
  int res=ElectronEnergyScale::AssignConstants(lines,_nEtaBins,_etaBinLimits,_dataConst,_dataConstErr,_mcConst1,_mcConst1Err,debug);
  for (int i=0; i<_nEtaBins; ++i) _dataConstRandomized[i]=0;

  if (res) {
    switch(_calibrationSet) {
    case CalSet_File_Voigt: 
      _mcConst2 = new double[_nEtaBins];
      _mcConst2Err = new double[_nEtaBins];
      _mcConst2Name = "gamma";
      res=ElectronEnergyScale::AssignConstants2(lines,_nEtaBins,"gamma_",_mcConst2,_mcConst2Err,debug);
      break;
    default:
      res=res; // dummy to prevent compiler complaints;
    }
  }

  if (!res) std::cout << "failed in this->assignConstants(lines)\n";
  return res;
}

//------------------------------------------------------

bool ElectronEnergyScale::AssignConstants(const std::vector<std::string> &lines, int count, double *eta, double *scale, double *scaleErr, double *smear, double *smearErr, int debug) {
  int etaDivCount=0;
  for (unsigned int i=0; !etaDivCount && (i<lines.size()); ++i) {
    if (lines[i].find("EtaDivisionCount")!=std::string::npos) {
      etaDivCount=atoi(lines[i].c_str()+lines[i].find('=')+1);
    }
  }
  // check etaDivCount count
  if ((etaDivCount!=count) && (etaDivCount*2!=count)) {
    std::cout << "assignConstants: got lines with etaDivCount=" << etaDivCount << ", while allocation states count=" << count << "\n";
    assert(0);
  }

  // Assign etaBins
  for (unsigned int i=0; i<lines.size(); ++i) {
    if (lines[i].find("! bins")!=std::string::npos) {
      std::stringstream s(lines[i].c_str()+7);
      int idx=(etaDivCount==count) ? 0 : etaDivCount;
      for ( ; idx<=count; ++idx) {
	s >> eta[idx];
      }
      if (etaDivCount*2==count) {
	idx=2*etaDivCount;
	for (int ii=0; ii<etaDivCount; ++ii) {
	  eta[ii] = - eta[idx-ii];
	}
      }
      if (fabs(eta[0]+2.5)<1e-6) eta[0]-=1e-5;
      if (fabs(eta[count]-2.5)<1e-6) eta[count]+=1e-5;
      break;
    }
  }

  // Assign scale/smear constants
  double *d,*derr;
  for (unsigned int i=0; i<lines.size(); ++i) {
    d=NULL; derr=NULL;
    if (lines[i].find("scale_")!=std::string::npos) {
      d=scale; derr=scaleErr;
    }
    else if (lines[i].find("smear_")!=std::string::npos) {
      d=smear; derr=smearErr;
    }
    if (d) {
      const char *s=lines[i].c_str();
      int idx=atoi(s + lines[i].find('_') + 1);
      double val=atof(s + lines[i].find(' '));
      double valErr=atof(s + lines[i].find(' ',lines[i].find('.')));
      if  (idx>count) {
	std::cout << "Error:\n";
	std::cout << "line=<" << lines[i] << ">\n";
	std::cout << "contains value for idx=" << idx << ", while count=" << count << "\n";
	std::cout << "in AssignConstants" << std::endl;
	assert(0);
      }
      if (etaDivCount==count) {
	d[idx]=val; derr[idx]=valErr;
      }
      else if (etaDivCount*2==count) {
	d[etaDivCount-idx-1]=val;
	derr[etaDivCount-idx-1]=valErr;
	d[etaDivCount+idx]=val;
	derr[etaDivCount+idx]=valErr;
      }
      else assert(0);
    }
  }
  if (debug) {
    std::cout << "got \n";
    for (unsigned int i=0; i<lines.size(); ++i) std::cout << "--> " << lines[i] << "\n";
    std::cout << "derived \n";
    for (int i=0; i<count; ++i) {
      std::cout << "scale_" << i << "  " << scale[i] << " " << scaleErr[i] << "\n";
    }
    for (int i=0; i<count; ++i) {
      std::cout << "smear_" << i << "  " << smear[i] << " " << smearErr[i] << "\n";
    }
  }
  return true;
}
//------------------------------------------------------

bool ElectronEnergyScale::AssignConstants2(const std::vector<std::string> &lines, int count, const char *parameterNameStart, double *par, double *parErr, int debug) {
  int etaDivCount=0;
  for (unsigned int i=0; !etaDivCount && (i<lines.size()); ++i) {
    if (lines[i].find("EtaDivisionCount")!=std::string::npos) {
      etaDivCount=atoi(lines[i].c_str()+lines[i].find('=')+1);
    }
  }
  // check etaDivCount count
  if ((etaDivCount!=count) && (etaDivCount*2!=count)) {
    std::cout << "assignConstants2: got lines with etaDivCount=" << etaDivCount << ", while allocation states count=" << count << "\n";
    return 0;
  }

  // Assign parameter constants
  double *d,*derr;
  for (unsigned int i=0; i<lines.size(); ++i) {
    d=NULL; derr=NULL;
    if (lines[i].find(parameterNameStart)!=std::string::npos) {
      d=par; derr=parErr;
    }
    if (d) {
      if (lines[i].find('_') == std::string::npos) {
	std::cout << "line <" << lines[i] << "> contains parameterNameStart=<" << parameterNameStart << ">, but has no '_'\n";
	return false;
      }
      std::stringstream ss(lines[i].c_str() + lines[i].find('_') + 1);
      int idx;
      double val,valErr;
      ss >> idx >> val >> valErr;
      if (etaDivCount==count) {
	d[idx]=val; derr[idx]=valErr;
      }
      else if (etaDivCount*2==count) {
	d[etaDivCount-idx-1]=val;
	derr[etaDivCount-idx-1]=valErr;
	d[etaDivCount+idx]=val;
	derr[etaDivCount+idx]=valErr;
      }
      else assert(0);
    }
  }
  if (debug) {
    std::cout << "got \n";
    for (unsigned int i=0; i<lines.size(); ++i) std::cout << "--> " << lines[i] << "\n";
    std::cout << "derived \n";
    for (int i=0; i<count; ++i) {
      std::cout << parameterNameStart << i << "  " << par[i] << " " << parErr[i] << "\n";
    }
  }
  return true;
}

//------------------------------------------------------

bool ElectronEnergyScale::initializeAllConstants(int debug){
  
  bool success = true;
  int nEtaBins1=0;
  //
  // Determine the number of bins
  //
  switch(_calibrationSet) {
  case UNCORRECTED: nEtaBins1=1; break;
  case Date20110901_EPS11_default: nEtaBins1=12; break;
  case Date20120101_default: nEtaBins1=12; break;
  case Date20120802_default: nEtaBins1=12; break;
  case Date20121003FEWZ_default: nEtaBins1=12; break;
  case Date20121025FEWZPU_default: nEtaBins1=12; break;
  case Date20140220_2012_j22_peak_position_IK : nEtaBins1=12; break;
  case Date20140220_2012_j22_peak_position : nEtaBins1=8; break;
  case Date20140220_2012_j22_peak_position_flat : nEtaBins1=1; break;
  case Date20130529_2012_j22_adhoc: nEtaBins1=8; break;
  case CalSet_File_Gauss: 
  case CalSet_File_Voigt:
  case CalSet_File_BreitWigner:
    nEtaBins1=-1; break;
  default:
    std::cout << "ElectronEnergyScale::initializeAllConstants: is not ready for the _calibrationSet= " << ElectronEnergyScale::CalibrationSetName(_calibrationSet,&_inpFileName) << " (" << int(_calibrationSet) << ")\n";
    return false;
  }

  if (nEtaBins1>0) {
    //
    // Allocate memory
    //
    _nEtaBins = nEtaBins1;
    _etaBinLimits = new double[nEtaBins1+1];
    _dataConst    = new double[nEtaBins1];
    _dataConstErr = new double[nEtaBins1];
    _dataConstRandomized = new double[nEtaBins1];
    switch(_calibrationSet) {
    default:
      _nMCConstants = 1;
      _mcConst1Name = "smear";
      _mcConst2Name.Clear();
      _mcConst3Name.Clear();
      _mcConst4Name.Clear();
      _mcConst1 = new double[nEtaBins1];
      _mcConst2 = 0;
      _mcConst3 = 0;
      _mcConst4 = 0;
      _mcConst1Err = new double[nEtaBins1];
      _mcConst2Err = 0;
      _mcConst3Err = 0;
      _mcConst4Err = 0;
    }
  }
  
  // 
  // Assign values
  // 
  switch( _calibrationSet ) { 
  case UNCORRECTED: {    
    if (nEtaBins1!=1) assert(0);
    assert(_etaBinLimits); 
    assert(_dataConst); assert(_dataConstErr);
    assert(_mcConst1); assert(_mcConst1Err);
    int i=0;
    _etaBinLimits[i] = -2.50001; _etaBinLimits[i+1]= 2.50001;
    _dataConst   [i] = 1.0;
    _dataConstErr[i] = 0.0;
    _mcConst1    [i] = 0.;
    _mcConst1Err [i] = 0.;
  }
    break;
  case Date20110901_EPS11_default: {
    //
    // Constants from energy scale calibrations
    // done for Summer 11 result. Note that the 
    // constants are symmetric about eta = 0.
    //
    const int nEtaBins = 12;
    const double etaBinLimits[nEtaBins+1] = 
      {-2.50001, -2.0, -1.5, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.5, 2.0, 2.50001};
    
    const double corrValues[nEtaBins] = 
      {1.04642, 1.00187, 1.01556, 1.00500, 1.00093, 1.00149, 1.00149, 1.00093, 1.00500, 1.01556, 1.00187, 1.04642};
    const double corrErrors[nEtaBins] = 
      {4.28928e-04,3.39718e-04,4.89342e-04,5.80480e-05,1.21192e-05,1.27489e-04,1.27489e-04,1.21192e-05,5.80480e-05,4.89342e-04,3.39718e-04,4.28928e-04};
    
    const double smearValues[nEtaBins] = 
      {2.05888e+00,1.46747e+00,1.14861e+00,7.63770e-01,5.04140e-01,5.27258e-01,5.27258e-01,5.04140e-01,7.63770e-01,1.14861e+00,1.46747e+00,2.05888e+00};
    const double smearErrors[nEtaBins] =
      {2.85889e-02,3.85260e-02,4.26451e-02,3.22979e-02,3.76972e-02,3.32377e-02,3.32377e-02,3.76972e-02,3.22979e-02,4.26451e-02,3.85260e-02,2.85889e-02};

    if (nEtaBins1!=nEtaBins) assert(0);
    assert(_etaBinLimits); 
    assert(_dataConst); assert(_dataConstErr);
    assert(_mcConst1); assert(_mcConst1Err);
    for(int i=0; i<nEtaBins; i++){
      _etaBinLimits[i] = etaBinLimits[i];
      _dataConst   [i] = corrValues[i];
      _dataConstErr[i] = corrErrors[i];
      _mcConst1    [i] = smearValues[i];
      _mcConst1Err [i] = smearErrors[i];
    }
    _etaBinLimits[nEtaBins] = etaBinLimits[nEtaBins];
  }
    break;

  case Date20120101_default: {
    //
    // Data of full 2011. Eta bins like for
    // for Summer 11 result. The
    // constants are NOT symmetric about eta = 0.
    //
    const int nEtaBins = 12;
    //const double etaBinLimits[nEtaBins+1] = 
    //  {-2.50001, -2.0, -1.5, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.5, 2.0, 2.50001};
    std::vector<std::string> lines;
    lines.push_back("! Date 2012 Jan 19\n");
    lines.push_back("! g_esfWorkCase=15 (6 bins on each eta side)");
    lines.push_back("! g_esfWorkCaseShortName= 6binNegs");
    lines.push_back("! g_esfFitModel=1 (fit model Gauss)");
    lines.push_back("scaling sqrt");
    lines.push_back("! bins  -2.50 -2.00 -1.50 -1.20 -0.80 -0.40 0.00 0.40 0.80 1.20 1.50 2.00 2.50");
    lines.push_back("MCOverData=15503.114069");
    lines.push_back("EtaDivisionCount=12");
    lines.push_back("ScalingFactorsCount=12");
    lines.push_back("SmearingFactorsCount=12");
    lines.push_back("scale_0      1.01214 -0.000303694");
    lines.push_back("scale_1     0.989212 -0.000268224");
    lines.push_back("scale_2      1.02211 -0.000280481");
    lines.push_back("scale_3      1.01145 -0.000168768");
    lines.push_back("scale_4      1.00602  -0.00014565");
    lines.push_back("scale_5      1.00573 -0.000144208");
    lines.push_back("scale_6      1.00293   -0.0001399");
    lines.push_back("scale_7      1.00375 -0.000146988");
    lines.push_back("scale_8      1.01098 -0.000170304");
    lines.push_back("scale_9      1.02083 -0.000283934");
    lines.push_back("scale_10     0.989319 -0.000270072");
    lines.push_back("scale_11      1.01284 -0.000304156");
    lines.push_back("smear_0      1.94897   -0.0226282");
    lines.push_back("smear_1      1.56441   -0.0266465");
    lines.push_back("smear_2       1.1956   -0.0284747");
    lines.push_back("smear_3     0.768533   -0.0221135");
    lines.push_back("smear_4     0.432683    -0.028314");
    lines.push_back("smear_5     0.525586   -0.0239693");
    lines.push_back("smear_6     0.455353   -0.0257488");
    lines.push_back("smear_7     0.471485   -0.0277377");
    lines.push_back("smear_8     0.742029   -0.0229296");
    lines.push_back("smear_9      1.19892   -0.0287477");
    lines.push_back("smear_10      1.49028   -0.0275794");
    lines.push_back("smear_11      2.03342   -0.0219498");

    if (nEtaBins1!=nEtaBins) assert(0);

    assert(_etaBinLimits); 
    assert(_dataConst); assert(_dataConstErr);
    assert(_mcConst1); assert(_mcConst1Err);
    //for(int i=0; i<nEtaBins+1; i++) _etaBinLimits[i] = etaBinLimits[i];
    if (!AssignConstants(lines, nEtaBins,_etaBinLimits,_dataConst,_dataConstErr,_mcConst1,_mcConst1Err)) assert(0);
  }
    break;

  case Date20120802_default: {
    //
    // Data of full 2011. Eta bins like for
    // for Summer 11 result. The
    // constants are NOT symmetric about eta = 0.
    // Values derived for the CMSSW44X files with correct isolation variables
    //
    const int nEtaBins = 12;
    //const double etaBinLimits[nEtaBins+1] = 
    //  {-2.50001, -2.0, -1.5, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.5, 2.0, 2.50001};
    std::vector<std::string> lines;
    lines.push_back("! Date 2012 Aug 02\n");
    lines.push_back("! g_esfWorkCase=15 (6 bins on each eta side)\n");
    lines.push_back("! g_esfWorkCaseShortName= 6binNegs\n");
    lines.push_back("! g_esfFitModel=1 (fit model Gauss)\n");
    lines.push_back("scaling sqrt\n");
    lines.push_back("! bins  -2.50 -2.00 -1.50 -1.20 -0.80 -0.40 0.00 0.40 0.80 1.20 1.50 2.00 2.50\n");
    lines.push_back("MCOverData=17365.021297\n");
    lines.push_back("EtaDivisionCount=12\n");
    lines.push_back("ScalingFactorsCount=12\n");
    lines.push_back("SmearingFactorsCount=12\n");
    lines.push_back("scale_0      1.01426 -0.000243511\n");
    lines.push_back("scale_1     0.988547 -0.000260672\n");
    lines.push_back("scale_2       1.0225 -0.000266414\n");
    lines.push_back("scale_3      1.01169 -0.000163012\n");
    lines.push_back("scale_4      1.00637 -0.000165987\n");
    lines.push_back("scale_5        1.006 -0.000153195\n");
    lines.push_back("scale_6      1.00313 -0.000102003\n");
    lines.push_back("scale_7      1.00402 -0.000127814\n");
    lines.push_back("scale_8      1.01126 -0.000140353\n");
    lines.push_back("scale_9      1.02128 -0.000265181\n");
    lines.push_back("scale_10     0.988757 -0.000262493\n");
    lines.push_back("scale_11      1.01436 -0.000265025\n");
    lines.push_back("smear_0      1.89664   -0.0207878\n");
    lines.push_back("smear_1      1.57576   -0.0247795\n");
    lines.push_back("smear_2      1.19665   -0.0271564\n");
    lines.push_back("smear_3     0.780325   -0.0210394\n");
    lines.push_back("smear_4     0.409839   -0.0292847\n");
    lines.push_back("smear_5     0.538319   -0.0227837\n");
    lines.push_back("smear_6     0.450927   -0.0251806\n");
    lines.push_back("smear_7     0.490343   -0.0256544\n");
    lines.push_back("smear_8     0.733717   -0.0220547\n");
    lines.push_back("smear_9      1.19763   -0.0274181\n");
    lines.push_back("smear_10      1.50618    -0.025792\n");
    lines.push_back("smear_11      1.98697   -0.0201719\n");

    if (nEtaBins1!=nEtaBins) assert(0);

    assert(_etaBinLimits); 
    assert(_dataConst); assert(_dataConstErr);
    assert(_mcConst1); assert(_mcConst1Err);
    //for(int i=0; i<nEtaBins+1; i++) _etaBinLimits[i] = etaBinLimits[i];
    if (!AssignConstants(lines, nEtaBins,_etaBinLimits,_dataConst,_dataConstErr,_mcConst1,_mcConst1Err)) assert(0);
  }
    break;

  case Date20121003FEWZ_default: {
    //
    // Data of full 2011. Eta bins like for
    // for Summer 11 result. The
    // constants are NOT symmetric about eta = 0.
    // Values derived for the CMSSW44X files with correct isolation variables.
    // FEWZ correction and MC backgrounds were considered in the derivation
    //
    const int nEtaBins = 12;
    //const double etaBinLimits[nEtaBins+1] = 
    //  {-2.50001, -2.0, -1.5, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.5, 2.0, 2.50001};
    std::vector<std::string> lines;
    lines.push_back("! Date 2012 Oct 03\n");
    lines.push_back("! g_esfWorkCase=15 (6 bins on each eta side)\n");
    lines.push_back("! g_esfWorkCaseShortName= 6binNegs\n");
    lines.push_back("! g_esfFitModel=1 (fit model Gauss)\n");
    lines.push_back("scaling sqrt\n");
    lines.push_back("! bins  -2.50 -2.00 -1.50 -1.20 -0.80 -0.40 0.00 0.40 0.80 1.20 1.50 2.00 2.50\n");
    lines.push_back("MCOverData=17137.515417\n");
    lines.push_back("EtaDivisionCount=12\n");
    lines.push_back("ScalingFactorsCount=12\n");
    lines.push_back("SmearingFactorsCount=12\n");
    lines.push_back("scale_0      1.01437 -0.000273902\n");
    lines.push_back("scale_1     0.988787 -0.000241357\n");
    lines.push_back("scale_2      1.02245 -0.000248856\n");
    lines.push_back("scale_3      1.01171 -0.000198996\n");
    lines.push_back("scale_4      1.00645 -0.000175147\n");
    lines.push_back("scale_5      1.00594 -0.000116214\n");
    lines.push_back("scale_6      1.00319 -0.000128479\n");
    lines.push_back("scale_7      1.00406 -0.000128213\n");
    lines.push_back("scale_8      1.01123  -0.00022103\n");
    lines.push_back("scale_9       1.0212 -0.000330949\n");
    lines.push_back("scale_10     0.988795 -0.000245423\n");
    lines.push_back("scale_11      1.01451  -0.00025763\n");
    lines.push_back("smear_0      1.86203   -0.0207486\n");
    lines.push_back("smear_1       1.4822   -0.0253977\n");
    lines.push_back("smear_2      1.09683    -0.028153\n");
    lines.push_back("smear_3     0.702673   -0.0224599\n");
    lines.push_back("smear_4     0.307007   -0.0341397\n");
    lines.push_back("smear_5     0.460252   -0.0250094\n");
    lines.push_back("smear_6     0.384336   -0.0276896\n");
    lines.push_back("smear_7     0.406731   -0.0291567\n");
    lines.push_back("smear_8     0.650285   -0.0240713\n");
    lines.push_back("smear_9      1.10309   -0.0289236\n");
    lines.push_back("smear_10      1.40315   -0.0266431\n");
    lines.push_back("smear_11      1.95493    -0.020064\n");

    if (nEtaBins1!=nEtaBins) assert(0);

    assert(_etaBinLimits); 
    assert(_dataConst); assert(_dataConstErr);
    assert(_mcConst1); assert(_mcConst1Err);
    //for(int i=0; i<nEtaBins+1; i++) _etaBinLimits[i] = etaBinLimits[i];
    if (!AssignConstants(lines, nEtaBins,_etaBinLimits,_dataConst,_dataConstErr,_mcConst1,_mcConst1Err)) assert(0);
  }
    break;

  case Date20121025FEWZPU_default: {
    //
    // Data of full 2011. Eta bins like for
    // for Summer 11 result. The
    // constants are NOT symmetric about eta = 0.
    // Values derived for the CMSSW44X files with correct isolation variables.
    // FEWZ correction and PU reweighting were considered in the derivation
    // Only MC Zee was used
    //
    const int nEtaBins = 12;
    //const double etaBinLimits[nEtaBins+1] = 
    //  {-2.50001, -2.0, -1.5, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.5, 2.0, 2.50001};
    std::vector<std::string> lines;
    lines.push_back("! Date 2012 Oct 25\n");
    lines.push_back("! g_esfWorkCase=15 (6 bins on each eta side)\n");
    lines.push_back("! g_esfWorkCaseShortName= 6binNegs\n");
    lines.push_back("! g_esfFitModel=1 (fit model Gauss)\n");
    lines.push_back("scaling sqrt\n");
    lines.push_back("! bins  -2.50 -2.00 -1.50 -1.20 -0.80 -0.40 0.00 0.40 0.80 1.20 1.50 2.00 2.50\n");
    lines.push_back("MCOverData=17136.430715");
    lines.push_back("EtaDivisionCount=12");
    lines.push_back("ScalingFactorsCount=12");
    lines.push_back("SmearingFactorsCount=12");
    lines.push_back("scale_0      1.01234 -0.000222138");
    lines.push_back("scale_1     0.987596 -0.000224354");
    lines.push_back("scale_2      1.02165 -0.000250768");
    lines.push_back("scale_3        1.011 -0.000161953");
    lines.push_back("scale_4      1.00578 -0.000135726");
    lines.push_back("scale_5      1.00535 -0.000107557");
    lines.push_back("scale_6      1.00251 -0.000115369");
    lines.push_back("scale_7      1.00341 -0.000133169");
    lines.push_back("scale_8      1.01052 -0.000192806");
    lines.push_back("scale_9      1.02047 -0.000320073");
    lines.push_back("scale_10     0.987861 -0.000233312");
    lines.push_back("scale_11      1.01263 -0.000264971");
    lines.push_back("smear_0      1.95695   -0.0201862");
    lines.push_back("smear_1      1.61288   -0.0242296");
    lines.push_back("smear_2      1.21198   -0.0266632");
    lines.push_back("smear_3     0.805396   -0.0205127");
    lines.push_back("smear_4     0.453941   -0.0268033");
    lines.push_back("smear_5      0.57275   -0.0215568");
    lines.push_back("smear_6     0.485194    -0.023758");
    lines.push_back("smear_7     0.520888   -0.0244349");
    lines.push_back("smear_8     0.761379   -0.0216066");
    lines.push_back("smear_9      1.21117   -0.0272647");
    lines.push_back("smear_10      1.54301   -0.0251034");
    lines.push_back("smear_11      2.04069   -0.0196909");

    if (nEtaBins1!=nEtaBins) assert(0);

    assert(_etaBinLimits); 
    assert(_dataConst); assert(_dataConstErr);
    assert(_mcConst1); assert(_mcConst1Err);
    //for(int i=0; i<nEtaBins+1; i++) _etaBinLimits[i] = etaBinLimits[i];
    if (!AssignConstants(lines, nEtaBins,_etaBinLimits,_dataConst,_dataConstErr,_mcConst1,_mcConst1Err)) assert(0);
  }
    break;

  case Date20130529_2012_j22_adhoc: {
    //
    // Data of full 2012 and Summer12 MC.
    // The constants are taking from:
    // MC:
    //   http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/EgammaAnalysis/ElectronTools/src/ElectronEnergyCalibrator.cc?revision=1.5&view=markup
    //   where the digMC is averaged between R9<0.94 and R9>=0.94 cases
    //   with weights 85% and 15%.
    //   Note: these are not in GeV, but unitless multiplicative constants
    // Data:
    //   Jan 22 rereco of 8 TeV
    //   Constants re kept at 1.00
    //
    const int nEtaBins = 8;
    const double etaBinLimits[nEtaBins+1] = 
      {-2.50001, -2.0, -1.479, -1.0, 0.0, 1.0, 1.479, 2.0, 2.50001};
    
    const double corrValues[nEtaBins] = 
      {1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00};
    // The numbers below, although approximate, should have worked, but they didn't
//       {0.9997, 0.9995, 1.0016, 0.9959, 0.9959, 1.0016, 0.9995, 0.9997};
    const double corrErrors[nEtaBins] = 
      {0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005, 0.005};
    
    const double smearValues[nEtaBins] = 
      {0.0282, 0.0221, 0.0176, 0.0094, 0.0094, 0.0176, 0.0221, 0.0282};
    const double smearErrors[nEtaBins] =
      {0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002, 0.002};

    if (nEtaBins1!=nEtaBins) assert(0);
    assert(_etaBinLimits); 
    assert(_dataConst); assert(_dataConstErr);
    assert(_mcConst1); assert(_mcConst1Err);
    for(int i=0; i<nEtaBins; i++){
      _etaBinLimits[i] = etaBinLimits[i];
      _dataConst   [i] = corrValues[i];
      _dataConstErr[i] = corrErrors[i];
      _mcConst1    [i] = smearValues[i];
      _mcConst1Err [i] = smearErrors[i];
    }
    _etaBinLimits[nEtaBins] = etaBinLimits[nEtaBins];
  }
    break;

  case Date20140220_2012_j22_peak_position_IK: { // another version
    //
    // Constants derived so that the Z peak position in data 
    // matches MC. All constants are the same for data. There is
    // no smearing for MC.
    // Done for 8 TeV data on Jan 22 rereco with momentum regression by Ilya
    //
    const int nEtaBins = 12;
    const double etaBinLimits[nEtaBins+1] = 
      {-2.50001, -2.0, -1.5, -1.2, -0.8, -0.4, 0.0, 0.4, 0.8, 1.2, 1.5, 2.0, 2.50001};

    const double corrValues[nEtaBins] = 
      {0.99841, 0.99841, 0.99841, 0.99841, 0.99841, 0.99841, 
       0.99841, 0.99841, 0.99841, 0.99841, 0.99841, 0.99841};

    const double corrErrors[nEtaBins] = 
      {0.00159, 0.00159, 0.00159, 0.00159, 0.00159, 0.00159,
       0.00159, 0.00159, 0.00159, 0.00159, 0.00159, 0.00159};

    const double smearValues[nEtaBins] = 
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

    const double smearErrors[nEtaBins] =
      {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    
    if (nEtaBins1!=nEtaBins) assert(0);
    assert(_etaBinLimits); 
    assert(_dataConst); assert(_dataConstErr);
    assert(_mcConst1); assert(_mcConst1Err);
    for(int i=0; i<nEtaBins; i++){
      _etaBinLimits[i] = etaBinLimits[i];
      _dataConst   [i] = corrValues[i];
      _dataConstErr[i] = corrErrors[i];
      _mcConst1    [i] = smearValues[i];
      _mcConst1Err [i] = smearErrors[i];
    }
    _etaBinLimits[nEtaBins] = etaBinLimits[nEtaBins];
  }
    break;

  case Date20140220_2012_j22_peak_position: {
    // Similar to Date20140220_2012_j22_peak_position_IK, but fewer
    // eta bins to match EGamma binning. Note: R9-based categories are ignored
    //
    const int nEtaBins = 8;
    const double etaBinLimits[8+1] = { 
      -2.50001, -2.0, -1.4442, -1.0, 0.0, 1.0, 1.4442, 2.0, 2.50001};

    const double corrValues[nEtaBins] = 
      {0.99841, 0.99841, 0.99841, 0.99841, 0.99841, 0.99841, 0.99841, 0.99841};

    const double corrErrors[nEtaBins] = 
      {0.00159, 0.00159, 0.00159, 0.00159, 0.00159, 0.00159, 0.00159, 0.00159};

    const double smearValues[nEtaBins] = {0, 0, 0, 0, 0, 0, 0, 0 };
    const double smearErrors[nEtaBins] = {0, 0, 0, 0, 0, 0, 0, 0 };
    
    if (nEtaBins1!=nEtaBins) assert(0);
    assert(_etaBinLimits); 
    assert(_dataConst); assert(_dataConstErr);
    assert(_mcConst1); assert(_mcConst1Err);
    for(int i=0; i<nEtaBins; i++){
      _etaBinLimits[i] = etaBinLimits[i];
      _dataConst   [i] = corrValues[i];
      _dataConstErr[i] = corrErrors[i];
      _mcConst1    [i] = smearValues[i];
      _mcConst1Err [i] = smearErrors[i];
    }
    _etaBinLimits[nEtaBins] = etaBinLimits[nEtaBins];
  }
    break;

  case Date20140220_2012_j22_peak_position_flat: {
    // Similar to Date20140220_2012_j22_peak_position_IK, but only 1 bin
    //
    const int nEtaBins = 1;
    const double etaBinLimits[1+1] = { -2.50001, 2.50001};
    const double corrValues[nEtaBins] = { 0.99841 };
    const double corrErrors[nEtaBins] = { 0.00159 };
    const double smearValues[nEtaBins] = { 0. };
    const double smearErrors[nEtaBins] = { 0. };

    if (nEtaBins1!=nEtaBins) assert(0);
    assert(_etaBinLimits);
    assert(_dataConst); assert(_dataConstErr);
    assert(_mcConst1); assert(_mcConst1Err);
    for(int i=0; i<nEtaBins; i++){
      _etaBinLimits[i] = etaBinLimits[i];
      _dataConst   [i] = corrValues[i];
      _dataConstErr[i] = corrErrors[i];
      _mcConst1    [i] = smearValues[i];
      _mcConst1Err [i] = smearErrors[i];
    }
    _etaBinLimits[nEtaBins] = etaBinLimits[nEtaBins];
  }
    break;

  case CalSet_File_Gauss: 
  case CalSet_File_Voigt:
  case CalSet_File_BreitWigner: {
    if (_inpFileName.Length()==0) {
      std::cout << "ElectronEnergyScale::initializeAllConstants. Calibration set " << ElectronEnergyScale::CalibrationSetName(_calibrationSet,&_inpFileName) << " requires input file to be set\n";
      assert(0);
    }
    success = loadInputFile(_inpFileName,debug);
  }
    break;

  default:
    std::cout << "ElectronEnergyScale::initializeAllConstants: is not ready for the _calibrationSet= " << ElectronEnergyScale::CalibrationSetName(_calibrationSet,&_inpFileName) << " (" << int(_calibrationSet) << ") [3]\n";
    return false;
  }
  
  return success;
}

//------------------------------------------------------

// The extra smearing function is to provide smearing for
// the mass of an event based on the individual parameters
// of two electrons. Thus the function is actually an 2D
// array of functions, with each pair of (i,j) eta bins for
// a given dielectron candidate corresponding to its unique
// smearing function.
bool ElectronEnergyScale::initializeExtraSmearingFunction(int normalize){

  bool success = true;
  // A sanity check. The function that initializes constants
  // should have been run by now.
  if( _nEtaBins <= 0 && _nEtaBins > nMaxFunctions)
    return false;

  for( int i=0; i<_nEtaBins; i++){
    for( int j=0; j<_nEtaBins; j++){
      TString fname = TString::Format("smearing_function_%03d_%03d", i, j);
      switch(_calibrationSet) {
      case UNCORRECTED: break;
      case Date20110901_EPS11_default:
      case Date20120101_default:
      case Date20120802_default:
      case Date20121003FEWZ_default:
      case Date20121025FEWZPU_default:
      case CalSet_File_Gauss: {
	if(_mcConst1 == 0) continue;
	smearingFunctionGrid[i][j] = new TF1(fname, "gaus(0)", -30, 30);
	smearingFunctionGrid[i][j]->SetNpx(500);
	double si = _mcConst1[i];
	double sj = _mcConst1[j];
	double sij=sqrt(si*si+sj*sj);
	double amp=(normalize) ? 1.0/(sij*sqrt(8*atan(1))) : 1.0;
	smearingFunctionGrid[i][j]->SetParameters(amp,0.0,sij);
      }
	break;
      case Date20140220_2012_j22_peak_position_IK:
      case Date20140220_2012_j22_peak_position:
      case Date20140220_2012_j22_peak_position_flat: {
 	if(_mcConst1 == 0) continue;
	double si = _mcConst1[i];
	double sj = _mcConst1[j];
	double sij=sqrt(si*si+sj*sj);
	if( sij == 0 ){
	  // This means no smearing is needed. However, the Gaussian
	  // generates weird random numbers if its width is zero.
	  // So simply constrain the range in this case.
	  smearingFunctionGrid[i][j] = new TF1(fname, "gaus(0)", -1e-6, 1e-6);
	}else{
	  smearingFunctionGrid[i][j] = new TF1(fname, "gaus(0)", -10, 10);
	}
	smearingFunctionGrid[i][j]->SetNpx(500);
	double amp=(normalize) ? 1.0/(sij*sqrt(8*atan(1))) : 1.0;
	smearingFunctionGrid[i][j]->SetParameters(amp,0.0,sij);
      }
	break;
      case Date20130529_2012_j22_adhoc: {
	if(_mcConst1 == 0) continue;
	// For this calibration set, only the first index out of the two is meaningful.
	smearingFunctionGrid[i][j] = new TF1(fname, "gaus(0)", -10, 10);
	smearingFunctionGrid[i][j]->SetNpx(500);
	double si = _mcConst1[i];
	double amp=(normalize) ? 1.0/(si*sqrt(8*atan(1))) : 1.0;
	smearingFunctionGrid[i][j]->SetParameters(amp,0.0,si);
      }
	break;
      case CalSet_File_BreitWigner: {
	if(_mcConst1 == 0) continue;
	smearingFunctionGrid[i][j] = new TF1(fname, "TMath::BreitWigner(x,[0],[1])", -10, 10);
	smearingFunctionGrid[i][j]->SetNpx(500);
	const double si = _mcConst1[i];
	const double sj = _mcConst1[j];
	const double sij=sqrt(si*si+sj*sj);
	smearingFunctionGrid[i][j]->SetParameters(0.,sij);
      }
	break;
      case CalSet_File_Voigt: {
	if((_mcConst1 == 0) || (_mcConst2 == 0)) continue;
	smearingFunctionGrid[i][j] = new TF1(fname, "TMath::Voigt(x,[0],[1])", -10, 10);
	smearingFunctionGrid[i][j]->SetNpx(500);
	const double si = _mcConst1[i];
	const double sj = _mcConst1[j];
	const double sij=sqrt(si*si+sj*sj);
	const double gammai = _mcConst2[i];
	const double gammaj = _mcConst2[j];
	const double gammaij = sqrt(gammai*gammai + gammaj*gammaj);
	smearingFunctionGrid[i][j]->SetParameters(sij,gammaij);
      }
	break;
      default:
	success = false;
      }
    } // end inner loop over eta bins
  } // end outer loop over eta bins

  return success;
}

//------------------------------------------------------

bool ElectronEnergyScale::isInitialized() const {
  bool yes=_isInitialized;
  if (yes) {
    switch(_calibrationSet) {
    case UNDEFINED: yes=false; break;
    case CalSet_File_Gauss:
    case CalSet_File_BreitWigner:
    case CalSet_File_Voigt:
      if (_inpFileName.Length()==0) yes=false;
      break;
    default:
      yes=true;
    }
  }
  return yes;
}

//------------------------------------------------------

int   ElectronEnergyScale::randomizeEnergyScaleCorrections(int seed){

  if( !_isInitialized ){
    printf("ElectronEnergyScale ERROR: the object is not properly initialized\n");
    return 0;
  }

  TRandom rand;
  rand.SetSeed(seed); std::cout << "\n\n\tSetSeed is called\n\n";
  _dataSeed=seed;
  _energyScaleCorrectionRandomizationDone = true;
  if (_calibrationSet==UNCORRECTED) return 1;

  int specialSeed=0;
  if (abs(seed)==111) {
    std::cout << "ElectronEnergyScale::randomizeEnergyScaleCorrections\n";
    std::cout << " - Randomization superseded by a special |seed|=111\n";
    std::cout << " - The fixed extreme error values will be taken\n";
    specialSeed=(seed>0) ? 1:-1;
  }

  if (_calibrationSet!=Date20140220_2012_j22_peak_position) {
    for(int i=0; i<_nEtaBins; i++){
      double delta=(specialSeed==0) ?
	rand.Gaus(0.0,_dataConstErr[i]) : specialSeed*_dataConstErr[i];
      _dataConstRandomized[i] = _dataConst[i] + delta;
    }
  }
  else {
    std::cout << "symmetrizing the randomized corrections\n";
    for(int i=0; i<_nEtaBins/2; i++){
      double delta=(specialSeed==0) ?
	rand.Gaus(0.0,_dataConstErr[i]) : specialSeed*_dataConstErr[i];
      _dataConstRandomized[i] = _dataConst[i] + delta;
      _dataConstRandomized[_nEtaBins-i-1]= _dataConstRandomized[i];
    }
  }

  return 1;
}

//------------------------------------------------------

// for systematics studies with 2012 data, MC has to be scaled by the
// inverse factors as compared to data

int ElectronEnergyScale::invertRandomizedEnergyScaleCorrections() {
  if (!_energyScaleCorrectionRandomizationDone) {
    std::cout << "invertRandomizedEnergyScaleCorrections: first call randomizeEnergyScaleCorrections\n";
    return 0;
  }
  for (int i=0; i<_nEtaBins; ++i) {
    double s=_dataConstRandomized[i];
    _dataConstRandomized[i] = 1/s;
  }
  return 1;
}

//------------------------------------------------------

// for systematics studies with 2012 data, MC has to be scaled by the
// inverse factors as compared to data

int ElectronEnergyScale::mirrorRandomizedEnergyScaleCorrections() {
  if (!_energyScaleCorrectionRandomizationDone) {
    std::cout << "mirrorRandomizedEnergyScaleCorrections: first call randomizeEnergyScaleCorrections\n";
    return 0;
  }
  for (int i=0; i<_nEtaBins; ++i) {
    double s0=_dataConst[i];
    double s1=_dataConstRandomized[i];
    _dataConstRandomized[i] = 1. + (s1-s0);
  }
  return 1;
}

//------------------------------------------------------

bool ElectronEnergyScale::setCalibrationSet(CalibrationSet calSet) {
  bool ok=kTRUE;
  if (isInitialized() && (calSet==UNCORRECTED)) {
    _calibrationSet = UNCORRECTED;
  }
  else {
    if (calSet==UNCORRECTED) {
      std::cout << "setCalibrationSet(" << ElectronEnergyScale::CalibrationSetName(calSet,NULL) << ") cannot be called for uninitialized object\n";
    }
    else {
      std::cout << "setCalibrationSet(calSet) cannot be called for " << ElectronEnergyScale::CalibrationSetName(calSet,NULL) << ". Use a constructor or init(calSet) instead.\n";
    }
    ok=kFALSE;
    assert(0);
  }
  return ok;
}


//------------------------------------------------------

TH1D* ElectronEnergyScale::createParamHisto(const TString &namebase, const TString &nameTag, const double *params, const double *paramErrs) const {
  int includeGap=0;
  const double gapHi=1.566;
  const double gapLo=1.4442;
  TString name=namebase; name+=nameTag;

  int binCount=_nEtaBins+3+2*includeGap;
  double bins[binCount];
  double vals[binCount];
  double valErrs[binCount];
  for (int i=0; i<binCount; ++i) vals[i]=0;
  for (int i=0; i<binCount; ++i) valErrs[i]=0;

  bins[0]=-3;
  int shift=1;
  for (int i=0; i<_nEtaBins+1; ++i) {
    int idx=i+shift;
    bins[idx]=_etaBinLimits[i];
    vals[idx]=params[i];
    valErrs[idx]=paramErrs[i];
    if (includeGap) {
      if (_etaBinLimits[i]==-1.5) {
	bins[idx]=-gapHi; 
	vals[idx]=0; valErrs[idx]=0;
	bins[idx+1]=-gapLo;
	vals[idx+1]=params[i]; valErrs[idx+1]=paramErrs[i];
	shift++;
      }
      else if (_etaBinLimits[i]==1.5) {
	bins[idx]=gapLo;
	vals[idx]=0; valErrs[idx]=0;
	bins[idx+1]=gapHi;
	vals[idx+1]=params[i]; valErrs[idx+1]=paramErrs[i];
	shift++;
      }
    }
    //std::cout << "idx=" << idx << ", bins[idx]=" << bins[idx] << ", vals[idx]=" << vals[idx] << "\n";
    //std::cout << "i+shift=" << (i+shift) << ", bins[i+shift]=" << bins[i+shift] << ", vals[i+shift]=" << vals[i+shift] << "\n";
  }
  vals[_nEtaBins+shift]=0;
  valErrs[_nEtaBins+shift]=0;
  bins[_nEtaBins+shift+1]=3.;
  TH1D* h=new TH1D(name.Data(),name.Data(), binCount-1, bins);
  h->SetStats(0);
  for (int i=0; i<binCount-2; ++i) {
    h->SetBinContent(i+1,vals[i]);
    h->SetBinError(i+1,valErrs[i]);
  }
  return h;
}

//------------------------------------------------------

TH1D* ElectronEnergyScale::createScaleHisto(const TString &namebase) const {
  TH1D* h= createParamHisto(namebase,"Scale",_dataConst,_dataConstErr);
  if (!h) {
    std::cout << "error in ElectornEnergyScale::createScaleHisto(" << namebase << ")\n";
  }
  return h;
}

//------------------------------------------------------

TH1D* ElectronEnergyScale::createSmearHisto(const TString &namebase, int parameterNo) const {
  char buf[10];
  sprintf(buf,"Smear%d",parameterNo);
  const double *data=(parameterNo==0) ? _mcConst1 : _mcConst2;
  const double *dataErr=(parameterNo==0) ? _mcConst1Err : _mcConst2Err;
  TH1D* h= createParamHisto(namebase,buf,data,dataErr);
  if (!h) {
    std::cout << "error in ElectornEnergyScale::createSmearHisto(" << namebase << ", " << parameterNo << ")\n";
  }
  return h;
}

//------------------------------------------------------

double ElectronEnergyScale::getEnergyScaleCorrection(double eta) const {

  double result = 1.0;
  bool randomize = _randomizedStudy;
  if (_calibrationSet != UNCORRECTED) {
    result = getEnergyScaleCorrectionAny(eta,randomize);
  }

  return result;
}

//------------------------------------------------------

double ElectronEnergyScale::getEnergyScaleCorrectionRandomized(double eta) const {

  double result = 1.0;
  if (_calibrationSet == UNCORRECTED) return result;

  if( !_energyScaleCorrectionRandomizationDone ){
    printf("ElectronEnergyScale ERROR: can not get randomized escale, randomization is not done\n");
    return result;
  }

  bool randomize = true;
  result = getEnergyScaleCorrectionAny(eta,randomize);

  return result;
}

//------------------------------------------------------

double ElectronEnergyScale::getEnergyScaleCorrectionAny(double eta, bool randomize) const {

  double result = 1.0;
  if (_calibrationSet == UNCORRECTED) return result;

  if( !_isInitialized ){
    printf("ElectronEnergyScale ERROR: the object is not properly initialized\n");
    return result;
  }

  for(int i=0; i<_nEtaBins; i++){
    if(eta >= _etaBinLimits[i] && eta < _etaBinLimits[i+1] ){
      if( !randomize )
	result = _dataConst[i];
      else
	result = _dataConstRandomized[i];
      break;
    }
  }

  return result;
}

//------------------------------------------------------

double ElectronEnergyScale::getEnergySmearCorrection(double eta, int the_set) const {
  int ibin=this->getEtaBinIdx(eta)-1;
  std::cout << "ibin=" << ibin << ", the_set=" << the_set << ", _nMCConstants=" << _nMCConstants << "\n";
  double result=0.;
  if ((ibin>=0) && (ibin<_nEtaBins)) {
    switch(the_set) {
    case 0: if (_mcConst1) result=_mcConst1[ibin]; break;
    case 1: if (_mcConst2) result=_mcConst2[ibin]; break;
    case 2: if (_mcConst3) result=_mcConst3[ibin]; break;
    case 3: if (_mcConst4) result=_mcConst4[ibin]; break;
    default: 
      std::cout << "getEnergySmearCorrection: wrong of the_set=" << the_set << "\n";
    }
  }
  std::cout << "returning result=" << result << "\n";
  return result;
}

//------------------------------------------------------

void ElectronEnergyScale::randomizeSmearingWidth(int seed){

  if( !_isInitialized ){
    printf("ElectronEnergyScale ERROR: the object is not properly initialized\n");
    return;
  }

  TRandom rand;
  rand.SetSeed(seed); std::cout << "\n\n\tSetSeed is called\n\n";
  _mcSeed=seed;
  _smearingWidthRandomizationDone = true;
  switch( _calibrationSet ) {
  case UNCORRECTED: break;

  case Date20110901_EPS11_default:
  case Date20120101_default:
  case Date20120802_default:
  case Date20121003FEWZ_default:
  case Date20121025FEWZPU_default:
  case Date20140220_2012_j22_peak_position_IK:
  case Date20140220_2012_j22_peak_position:
  case Date20140220_2012_j22_peak_position_flat:
  case CalSet_File_Gauss: {

    for( int i=0; i<_nEtaBins; i++){
      for( int j=0; j<_nEtaBins; j++){
	TString fname = TString::Format("smearing_function_randomized_%03d_%03d", i, j);
	smearingFunctionGridRandomized[i][j] = new TF1(fname, "gaus(0)", -10, 10);
	smearingFunctionGridRandomized[i][j]->SetNpx(500);
	double si = _mcConst1[i] + rand.Gaus(0.0,_mcConst1Err[i]);
	double sj = _mcConst1[j] + rand.Gaus(0.0,_mcConst1Err[j]);
	double sij= sqrt(si*si+sj*sj);
	smearingFunctionGridRandomized[i][j]->SetParameters(1.0/(sij*sqrt(8*atan(1))),0.0,sij);
	if (i>j) { smearingFunctionGridRandomized[i][j]->SetParameters(smearingFunctionGridRandomized[j][i]->GetParameters()); }
      } // end inner loop over eta bins
    } // end outer loop over eta bins
  }
    break;
  case Date20130529_2012_j22_adhoc: {
    // For this calibration set, only the first index in the functions grid is meaningful
    for( int i=0; i<_nEtaBins; i++){
      for( int j=0; j<_nEtaBins; j++){
	TString fname = TString::Format("smearing_function_randomized_%03d_%03d", i, j);
	smearingFunctionGridRandomized[i][j] = new TF1(fname, "gaus(0)", -10, 10);
	smearingFunctionGridRandomized[i][j]->SetNpx(500);
	double si = _mcConst1[i] + rand.Gaus(0.0,_mcConst1Err[i]);
	smearingFunctionGridRandomized[i][j]->SetParameters(1.0/(si*sqrt(8*atan(1))),0.0,si);
	if (i>j) { smearingFunctionGridRandomized[i][j]->SetParameters(smearingFunctionGridRandomized[j][i]->GetParameters()); }
      } // end inner loop over eta bins
    } // end outer loop over eta bins
  }
    break;
  default:
    // This place should be never reached. This is just a sanity check.
    printf("ElectronEnergyScale ERROR: failed to created randomized smearing function\n");
    throw 1;
  }

  return;
}

//------------------------------------------------------

double ElectronEnergyScale::generateMCSmear(double eta1, double eta2) const {
  
  bool randomize = false; // this is not foreseen to be flagged by _randomizedStudy!!
  if (_calibrationSet == UNCORRECTED) return 0.0;

  double result = generateMCSmearAny(eta1, eta2, randomize);

  return result;
}

//------------------------------------------------------

double ElectronEnergyScale::generateMCSmearRandomized(double eta1, double eta2) const {

  double result = 0.0;
  if (_calibrationSet == UNCORRECTED) return 0.0;

  if( !_smearingWidthRandomizationDone ){
    printf("ElectronEnergyScale ERROR: can not get randomized smear, randomization is not done\n");
    return result;
  }
  
  bool randomize = true;
  result = generateMCSmearAny(eta1, eta2, randomize);

  return result;
}

//------------------------------------------------------

double ElectronEnergyScale::generateMCSmearAny(double eta1, double eta2, bool randomize) const {
  
  double result = 0;
  if (_calibrationSet == UNCORRECTED) return result;

  if( !_isInitialized ){
    printf("ElectronEnergyScale ERROR: the object is not properly initialized\n");
    return result;
  }
  
  int count = 0;
  int ibin = 0;
  int jbin = 0;
  for(int i=0; i<_nEtaBins; i++){
    if(eta1 >= _etaBinLimits[i] && eta1 < _etaBinLimits[i+1] ){
      ibin = i;
      count++;
    }
    if(eta2 >= _etaBinLimits[i] && eta2 < _etaBinLimits[i+1] ){
      jbin = i;
      count++;
    }
  }
  if(count != 2) {
    printf("ElectronEnergyScale: Smear function ERROR\n");
    printf("Failed to obtain index for eta1=%4.2lf, eta2=%4.2lf\n",eta1,eta2);
    throw 1;
  }
 
  if( !randomize)
    result = smearingFunctionGrid[ibin][jbin]->GetRandom();
  else
    result = smearingFunctionGridRandomized[ibin][jbin]->GetRandom();

  return result;
}

//------------------------------------------------------

double ElectronEnergyScale::generateMCSmearSingleEle(double eta) const {
  
  bool randomize = false; // this is not foreseen to be flagged by _randomizedStudy!!
  if (_calibrationSet == UNCORRECTED) return 0.0;

  double result = generateMCSmearAnySingleEle(eta, randomize);

  return result;
}

//------------------------------------------------------

double ElectronEnergyScale::generateMCSmearSingleEleRandomized(double eta) const {

  double result = 0.0;
  if (_calibrationSet == UNCORRECTED) return 0.0;

  if( !_smearingWidthRandomizationDone ){
    printf("ElectronEnergyScale ERROR: can not get randomized smear, randomization is not done\n");
    return result;
  }
  
  bool randomize = true;
  result = generateMCSmearAnySingleEle(eta, randomize);

  return result;
}

//------------------------------------------------------

double ElectronEnergyScale::generateMCSmearAnySingleEle(double eta, bool randomize) const {
  
  double result = 0;
  if (_calibrationSet == UNCORRECTED) return result;

  if( !_isInitialized ){
    printf("ElectronEnergyScale ERROR: the object is not properly initialized\n");
    return result;
  }
  
  int ibin = -1;
  for(int i=0; i<_nEtaBins; i++){
    if(eta >= _etaBinLimits[i] && eta < _etaBinLimits[i+1] ){
      ibin = i;
    }
  }
  if(ibin == -1) {
    printf("ElectronEnergyScale: Smear function ERROR\n");
    printf("Failed to obtain index for eta=%4.2lf\n",eta);
    throw 1;
  }
 
  // The second bin can be any, only the first bin is meaningful
  if( !randomize)
    result = smearingFunctionGrid[ibin][ibin]->GetRandom();
  else
    result = smearingFunctionGridRandomized[ibin][ibin]->GetRandom();

  return result;
}

//------------------------------------------------------

bool ElectronEnergyScale::addSmearedWeightAny(TH1D *hMass, int eta1Bin, int eta2Bin, double mass, double weight, bool randomize) const {
  
  //std::cout << "mass=" << mass << ", weight=" << weight << "\n";
  if( !_isInitialized ){
    printf("ElectronEnergyScale ERROR: the object is not properly initialized\n");
    return kFALSE;
  }
  if (!hMass) {
    std::cout << "this subroutine will do nothing for hMass=NULL\n";
    assert(hMass);
  }
  
  if (_calibrationSet == UNCORRECTED) {
    hMass->Fill(mass,weight);
    return kTRUE;
  }

  if (randomize && !this->isSmearRandomized()) {
    std::cout << "ElectronEnergyScale ERROR: the smearing was not randomized\n";
    return kFALSE;
  }

  eta1Bin--; eta2Bin--;
  assert((eta1Bin>=0)); assert((eta2Bin>=0));
  TF1 *smearFnc = (randomize) ? 
    smearingFunctionGridRandomized[eta1Bin][eta2Bin] :
    smearingFunctionGrid[eta1Bin][eta2Bin];

  TH1D *h=hMass;
  for (int i=1; i<=h->GetNbinsX(); i++) {
    const double xa=h->GetBinLowEdge(i);
    const double xw=h->GetBinWidth(i);
    const double w= smearFnc->Integral( xa-mass, xa-mass+xw );
    h->Fill(xa+0.5*xw, w * weight);
    //std::cout << "adding " << w *weight << " in " << (xa+0.5*xw) << "\n";
  }

  return kTRUE;
}

//------------------------------------------------------

void ElectronEnergyScale::smearDistributionAny(TH1D *destination, int eta1Bin, int eta2Bin, const TH1D *source, bool randomize) const {
  assert(source); assert(destination);
  for (int i=1; i<source->GetNbinsX(); ++i) {
    assert(addSmearedWeightAny(destination,eta1Bin,eta2Bin,source->GetBinCenter(i),source->GetBinContent(i),randomize));
  }
}

//------------------------------------------------------

void ElectronEnergyScale::print() const {

  printf("\nEnergy scale corrections used:\n");
  printf("   Calibration set (%d): %s\n", _calibrationSet, ElectronEnergyScale::CalibrationSetName(this->_calibrationSet,&this->_inpFileName).Data());
  printf("   Smearing function: %s\n",ElectronEnergyScale::CalibrationSetFunctionName(this->_calibrationSet).Data());
  if (_randomizedStudy) printf("  Random data seed=%d\n",_dataSeed);
  if (_smearingWidthRandomizationDone) printf("  Random MC seed=%d\n",_mcSeed);
  printf("   Constants:\n");
  printf("     eta-bin      Escale-const      MC-const-1");
  if ( _mcConst2 ) printf("          MC-const-2");
  if ( _mcConst3 ) printf("          MC-const-3");
  if ( _mcConst4 ) printf("          MC-const-4"); 
  if (_randomizedStudy) printf("          Random Escale-const");
  printf("\n");
  printf("                              %16s    %16s    %16s   %16s\n",
	 _mcConst1Name.Data(), _mcConst2Name.Data(), _mcConst3Name.Data(), _mcConst4Name.Data());
  for(int i=0; i<_nEtaBins; i++){
    printf("   %5.2f- %5.2f  %6.4f+-%5.4f",
    //printf("   %8.5f- %8.5f  %6.4f+-%5.4f",
	   _etaBinLimits[i], _etaBinLimits[i+1],
	   _dataConst[i], _dataConstErr[i]);
    if( _mcConst1 != 0 )
      printf("   %3.1e+- %3.1e", _mcConst1[i], _mcConst1Err[i]);
    if( _mcConst2 != 0 )
      printf("   %3.1e+- %3.1e", _mcConst2[i], _mcConst2Err[i]);
    if( _mcConst3 != 0 )
      printf("   %3.1e+- %3.1e", _mcConst3[i], _mcConst3Err[i]);
    if( _mcConst4 != 0 )
      printf("   %3.1e+- %3.1e", _mcConst4[i], _mcConst4Err[i]);
    if (_randomizedStudy) {
      double rndVal=(_dataConstRandomized) ? _dataConstRandomized[i]  : 0.;
      double mainVal=_dataConst[i];
      if (mainVal==0.0) mainVal=1e-14;
      printf("     %6.4f (% 9.3e%%)", rndVal, (rndVal-mainVal)*100./mainVal);
    }
    printf("\n");
  }
  printf("\n");
  
  return;
}

//------------------------------------------------------

void ElectronEnergyScale::printAsTexTable(const TString &fname) const {
  std::cout << "printAsTexTable(" << fname << ") WARNING: assuming Gaussian distribution\n";
  assert( !_mcConst2 && !_mcConst3 && !_mcConst4 );

  std::ofstream fout(fname.Data());
  assert(fout.is_open());
  fout << "\\begin{table}[htbH]\n";
  fout << "\\begin{center}\n";
  fout << "\\caption{Energy scale correction factors for individual electrons\n";
  fout << "in data as a function of $\\eta$, and extra smearing required for\n";
  fout << "MC to match single electron energy resolution in data.}\n";
  fout << "\\label{tab:energyScaleFactors}\n";
  fout << "\\begin{tabular}{|ccc|c|c|}\n";
  fout << "\\hline\n";
  const char *data_Rnd=(_randomizedStudy) ? "randomized" : "          ";
  const char *mc_Rnd=(_smearingWidthRandomizationDone) ? "randomized" : "";
  fout << "\\multicolumn{3}{|c|}{} &  " << data_Rnd << "  & " << mc_Rnd << " resolution  \\\\\n";
  fout << "\\multicolumn{3}{|c|}{region in $\\eta$} & scale $s_i$  & underestimate $\\sigma_i$, GeV \\\\\n";
  fout << "\\multicolumn{3}{|c|}{} &  (data)     &  (MC)   \\\\\n";
  fout << "\\hline\n";

  printf("   Calibration set (%d): %s\n", _calibrationSet, ElectronEnergyScale::CalibrationSetName(this->_calibrationSet,&this->_inpFileName).Data());
  printf("   Smearing function: %s\n",ElectronEnergyScale::CalibrationSetFunctionName(this->_calibrationSet).Data());
  if (_randomizedStudy) printf("  Random data seed=%d\n",_dataSeed);
  if (_smearingWidthRandomizationDone) printf("  Random MC seed=%d\n",_mcSeed);
  if (_randomizedStudy) printf("          Random Escale-const");
  printf("\n");

  for(int i=0; i<_nEtaBins; i++){
    double eta1=_etaBinLimits[i];
    double eta2=_etaBinLimits[i+1];
    char star1=' ';
    char star2=' ';
    if ((fabs(fabs(eta1)-DYTools::kECAL_GAP_LOW)<1e-4) ||
	(fabs(fabs(eta1)-DYTools::kECAL_GAP_HIGH)<1e-4) ||
	(fabs(fabs(eta1)-1.5)<1e-4)) star1='*';
    if ((fabs(fabs(eta2)-DYTools::kECAL_GAP_LOW)<1e-4) ||
	(fabs(fabs(eta2)-DYTools::kECAL_GAP_HIGH)<1e-4) ||
	(fabs(fabs(eta2)-1.5)<1e-4)) star2='*';
    std::cout << "eta1=" << eta1 << "star1=<" << star1 << ">\n";
    
    fout << Form(" $% 4.1lf$%c &-& $%4.1lf$%c",
		 eta1,star1, eta2,star2);
    fout << Form(" &  $%6.4lf\\pm%6.4lf$",
		 _dataConst[i], _dataConstErr[i]);
    fout << Form("  &  $%4.2lf\\pm%4.2lf$ ",
		 _mcConst1[i], _mcConst1Err[i]);
    fout << "\\\\\n";
  }
  fout << "\\hline\n";
  fout << "\\multicolumn{5}{l}{ * The $|\\eta|$ gap (1.4442-1.566) is excluded}\n";
  fout << "\\end{tabular}\n";
  fout << "\\end{center}\n";
  fout << "\\end{table}\n";
  fout.close();
  std::cout << "file <" << fname << "> created\n";
  return;
}

//------------------------------------------------------
//------------------------------------------------------

TString ElectronEnergyScale::ExtractFileName(const TString &escaleTagName) {
  TString fname;
  if (escaleTagName.Contains("File")) {
    std::string tmp=escaleTagName.Data();
    size_t pos1=tmp.find("File");
    char target[10];
    sprintf(target," \t\n");
    pos1=tmp.find_first_of(' ',pos1+1);
    if (pos1==std::string::npos) { 
      pos1=tmp.find_first_of('_',pos1+6); 
      target[0]='_';
    }
    pos1=tmp.find_first_not_of(" \t\n",pos1+1); // skip File* and spaces
    size_t pos2=tmp.find_first_of(" \n\t",pos1+1);
    if (pos2==std::string::npos) pos2=tmp.size(); // should be -1, but +1 is a correction for the case when the string terminates with the file name only
    if ((pos1==std::string::npos) || (pos2==std::string::npos)) {
      std::cout << "ElectronEnergyScale::ExtractFileName: since escaleTagName contains 'File', the keyword has to be followed by a file name\n";
      std::cout << "  escaleTagName=<" << escaleTagName << ">\n";
      throw 3;
    }
    else {
      fname=tmp.substr(pos1,pos2-pos1);
      //std::cout << "fileName=<" << fileName << ">\n";
    }
  }
  return fname;
}

//------------------------------------------------------

ElectronEnergyScale::CalibrationSet ElectronEnergyScale::DetermineCalibrationSet(const TString &escaleTagName_orig, TString *inputFileName) {
  ElectronEnergyScale::CalibrationSet calibrationSet  = ElectronEnergyScale::UNDEFINED;
  TString fileName;
  TString escaleTagName;
  //std::cout << "entered DetermineCalibrationSet(\"" << escaleTagName_orig << "\", " << ((inputFileName) ? "ptr":"NULL") << ")\n";
  Ssiz_t pos = escaleTagName_orig.Index("#");
  if (pos>0) escaleTagName=escaleTagName_orig(0,pos-1); 
  else escaleTagName=escaleTagName_orig;
  std::cout << "escaleTagName_orig={" << escaleTagName_orig << "}, escaleTagName={" << escaleTagName << "}; pos=" << pos << "\n";
  if ( (pos==0) || escaleTagName.Contains("Date20120802_default") ) {
    calibrationSet = ElectronEnergyScale::Date20120802_default;
  }
  else if ( (pos==0) || escaleTagName.Contains("Date20121025FEWZPU_default") ) {
    calibrationSet = ElectronEnergyScale::Date20121003FEWZ_default;
  }
  else if ( (pos==0) || escaleTagName.Contains("Date20121003FEWZ_default") ) {
    calibrationSet = ElectronEnergyScale::Date20121003FEWZ_default;
  }
  else if ( (pos==0) || escaleTagName.Contains("Date20120101_default") ) {
    calibrationSet = ElectronEnergyScale::Date20120101_default;
  }
  else if ( (pos==0) || escaleTagName.Contains("Date20130529_2012_j22_adhoc") ) {
    calibrationSet = ElectronEnergyScale::Date20130529_2012_j22_adhoc;
  }
  else if ( (pos==0) || escaleTagName.Contains("Date20140220_2012_j22_peak_position_IK") ) {
    calibrationSet = ElectronEnergyScale::Date20140220_2012_j22_peak_position_IK;
  }
  else if ( (pos==0) || escaleTagName.Contains("Date20140220_2012_j22_peak_position_flat") ) {
    calibrationSet = ElectronEnergyScale::Date20140220_2012_j22_peak_position_flat;
  }
  else if ( (pos==0) || escaleTagName.Contains("Date20140220_2012_j22_peak_position") ) {
    calibrationSet = ElectronEnergyScale::Date20140220_2012_j22_peak_position;
  }
  else if ( escaleTagName.Contains("UNCORRECTED")) {
    calibrationSet = ElectronEnergyScale::UNCORRECTED;
  }
  else if ( escaleTagName.Contains("Date20110901_EPS11_default")) {
    calibrationSet = ElectronEnergyScale::Date20110901_EPS11_default;
  }
  else if (escaleTagName.Contains("File")) {
    if (escaleTagName.Contains("File_Gauss") ||
	escaleTagName.Contains("FileGauss")) {
      fileName=ExtractFileName(escaleTagName);
      if (fileName.Length()>0) {
	calibrationSet = ElectronEnergyScale::CalSet_File_Gauss;
      }
    }
    else if (escaleTagName.Contains("File_Voigt") ||
	     escaleTagName.Contains("FileVoigt")) {
      fileName=ExtractFileName(escaleTagName);
      if (fileName.Length()>0) {
	calibrationSet = ElectronEnergyScale::CalSet_File_Voigt;
      }
    }
    else if (escaleTagName.Contains("File_BreitWigner") ||
	     escaleTagName.Contains("FileBreitWigner")) {
      fileName=ExtractFileName(escaleTagName);
      if (fileName.Length()>0) {
	calibrationSet = ElectronEnergyScale::CalSet_File_BreitWigner;
      }
    }
  }
  else{
    //printf("Failed to match escale calibration. Tag: >>%s<<\n", escaleTagName.Data());
    return ElectronEnergyScale::UNDEFINED;
  }
    
  if (fileName.Length()) {
    if (inputFileName) *inputFileName=fileName;
    else {
      std::cout << "warning: ElectronEnergyScale::DetermineCalibrationSet located a keyword meaning the file name will be provided, yet no pointer to the container was supplied\n";
    }
  }
  //std::cout << "DetermineCalibrationSet(" << escaleTagName << ") returns " << ElectronEnergyScale::CalibrationSetName(calibrationSet,&fileName) << " (" << int(calibrationSet) << ")\n";
  return calibrationSet;
}

//------------------------------------------------------

TString ElectronEnergyScale::CalibrationSetName(ElectronEnergyScale::CalibrationSet escaleTag, const TString *fileName) {
  TString name="UNDEFINED";
  switch(escaleTag) {
  case ElectronEnergyScale::UNDEFINED: break;
  case ElectronEnergyScale::UNCORRECTED: name="UNCORRECTED"; break;
  case ElectronEnergyScale::Date20110901_EPS11_default: name="Date20110910_EPS11_default"; break;
  case ElectronEnergyScale::Date20120101_default: name="Date20120101_default"; break;
  case ElectronEnergyScale::Date20120802_default: name="Date20120802_default"; break;
  case ElectronEnergyScale::Date20121003FEWZ_default: name="Date20121003FEWZ_default"; break;
  case ElectronEnergyScale::Date20121025FEWZPU_default: name="Date20121025FEWZPU_default"; break;
  case ElectronEnergyScale::Date20130529_2012_j22_adhoc: name="Date20130529_2012_j22_adhoc"; break;
  case ElectronEnergyScale::Date20140220_2012_j22_peak_position_IK: name="Date20140220_2012_j22_peak_position_IK"; break;
   case ElectronEnergyScale::Date20140220_2012_j22_peak_position: name="Date20140220_2012_j22_peak_position"; break;
   case ElectronEnergyScale::Date20140220_2012_j22_peak_position_flat: name="Date20140220_2012_j22_peak_position_flat"; break;
  case ElectronEnergyScale::CalSet_File_Gauss: 
    name="FileGauss(";
    if (fileName) name+=(*fileName);
    name+=")";
    break;
  case ElectronEnergyScale::CalSet_File_Voigt: 
    name="FileVoigt(";
    if (fileName) name+=(*fileName);
    name+=")";
    break;
  case ElectronEnergyScale::CalSet_File_BreitWigner: 
    name="FileBreitWigner(";
    if (fileName) name+=(*fileName);
    name+=")";
    break;
  default:
    name="CalibrationSetName_undetermined";
  }
  return name;
}

//------------------------------------------------------

TString ElectronEnergyScale::CalibrationSetFunctionName(ElectronEnergyScale::CalibrationSet escaleTag) {
  TString name="Gauss";
  switch (escaleTag) {
  case ElectronEnergyScale::UNDEFINED: name="undefined"; break;
  case ElectronEnergyScale::UNCORRECTED: name="uncorrected"; break;
  case ElectronEnergyScale::Date20110901_EPS11_default: break; // Gauss
  case ElectronEnergyScale::Date20120101_default: break; // Gauss
  case ElectronEnergyScale::Date20120802_default: break; // Gauss
  case ElectronEnergyScale::Date20121003FEWZ_default: break; // Gauss
  case ElectronEnergyScale::Date20121025FEWZPU_default: break; // Gauss
  case ElectronEnergyScale::Date20130529_2012_j22_adhoc: break; // Gauss
  case ElectronEnergyScale::Date20140220_2012_j22_peak_position_IK: break; // Gauss
  case ElectronEnergyScale::Date20140220_2012_j22_peak_position: break; // Gauss
  case ElectronEnergyScale::Date20140220_2012_j22_peak_position_flat: break; // Gauss
  case ElectronEnergyScale::CalSet_File_Gauss: break; // Gauss
  case ElectronEnergyScale::CalSet_File_Voigt: name="Voigt"; break;
  case ElectronEnergyScale::CalSet_File_BreitWigner: name="BreitWigner"; break;
  default:
    name="CalibrationSetFunctionName_undetermined";
  }
  return name;
}

//------------------------------------------------------

TString ElectronEnergyScale::calibrationSetShortName() const {
  TString name;
  switch (this->_calibrationSet) {
  case ElectronEnergyScale::UNDEFINED: name="undefined"; break;
  case ElectronEnergyScale::UNCORRECTED: name="uncorrected"; break;
  case ElectronEnergyScale::Date20110901_EPS11_default: name="default20110901"; break; 
  case ElectronEnergyScale::Date20120101_default: name="default20120101"; break;
  case ElectronEnergyScale::Date20120802_default: name="default20120802"; break;
  case ElectronEnergyScale::Date20121003FEWZ_default: name="default20121003FEWZ"; break;
  case ElectronEnergyScale::Date20121025FEWZPU_default: name="default20121025FEWZPU"; break;
  case ElectronEnergyScale::Date20130529_2012_j22_adhoc: name="adhoc20130529j22"; break;
  case ElectronEnergyScale::Date20140220_2012_j22_peak_position_IK: name="peak20140220IK"; break;
  case ElectronEnergyScale::Date20140220_2012_j22_peak_position: name="peak20140220"; break;
  case ElectronEnergyScale::Date20140220_2012_j22_peak_position_flat: name="peak20140220flat"; break;
  case ElectronEnergyScale::CalSet_File_Gauss: name="Gauss"; break;
  case ElectronEnergyScale::CalSet_File_Voigt: name="Voigt"; break;
  case ElectronEnergyScale::CalSet_File_BreitWigner: name="BreitWigner"; break;
  default:
    name="calibrationSetShortName_undetermined";
  }
  int len=this->_inpFileName.Length();
  if (len) {
    int i=this->_inpFileName.Index("testESF");
    if (i>0) {
      TString temp=this->_inpFileName(i+7, len);
      int idot=temp.Index(".");
      //std::cout << "i=" << i << ", idot=" << idot << ": " << this->_inpFileName(i,100) << "\n";
      TString temp2=temp(0,idot);
      //std::cout << "\ntemp2=<" << temp2 << ">\n";
      if (temp2.Index(name)>0) name=temp2(1,temp2.Length());
      else name+=temp2;
    }
  }
  if (_energyScaleCorrectionRandomizationDone) { name+="_DtRND"; name+=_dataSeed; }
  if (_smearingWidthRandomizationDone) { name+="_McRND"; name+=_mcSeed; }
  return name;
}

//------------------------------------------------------

#ifdef UseEEM
int ElectronEnergyScale::loadEEMFile(const TString &eemFileName, std::vector<std::vector<double>*> &eemData) const {
  EtaEtaMassData_t *eem = new EtaEtaMassData_t();
  int res=1;
  int etaEtaCount = this->numberOfEtaEtaBins();
  eemData.clear(); eemData.reserve(etaEtaCount+1);

  // Read file twice. First time we will determine the number of masses
  // in each (eta1,eta2) bin and allocate memory. The second time we will
  // store the mass value
  const int optimize=1; // optimize memory
  for (int loop=0; res && (loop<1+optimize); ++loop) {
    //std::cout << "loop=" << loop << std::endl;
    std::vector<int> counts(etaEtaCount);
    const char *fname=eemFileName.Data();
    TFile *fin = new TFile(fname);
    assert(fin);
    TTree *tree = (TTree*)fin->Get("Data"); assert(tree);
    tree->SetBranchAddress("Data",&eem);
    TBranch *branch = tree->GetBranch("Data"); assert(branch);
    UInt_t entryMax = tree->GetEntries();
    for (UInt_t ientry=0; ientry < entryMax; ++ientry) {
      branch->GetEntry(ientry);
      int idx=this->getEtaEtaIdx(eem->eta1(),eem->eta2());
      //std::cout << loop << " got " << (*eem) << ", idx=" << idx << "\n";
      if ((idx>=0) && (idx < etaEtaCount)) {
	if (optimize && (loop==0)) {
	  // loops 0 -- count the number of events
	  counts[idx]++;
	}
	else {
	  // loops 1 -- store the memory
	  std::vector<double>* store=eemData[idx];
	  store->push_back(eem->mass());
	}
      }
    }
    //delete branch;
    //delete tree;
    delete fin;
    if (optimize && (loop==0)) {
      // allocate memory
      unsigned int k=0;
      for (int i=0; i<this->numberOfEtaBins(); ++i) {
	for (int j=i; j<this->numberOfEtaBins(); ++j, ++k) {
	  std::vector<double>* tmp=new std::vector<double>();
	  assert(tmp);
	  //std::cout << "counts[k]=" << counts[k] << "\n";
	  tmp->reserve(counts[k]);
	  eemData.push_back(tmp);
	}
      }
    }
  }
  return 1;
}
#endif

//------------------------------------------------------

#ifdef UseEEM
int ElectronEnergyScale::loadEEMFile(const TString &eemFileName, std::vector<std::vector<double>*> &eemData, double massMin, double massMax) const {
  ClearVec(eemData);
  std::vector<std::vector<EtaEtaMassData_t>*> data;
  int res=this->loadEEMFile(eemFileName,data,massMin,massMax);
  if (!res) {
    std::cout << "error in ElectronEnergyScale::loadEEMFile(doubleVV,min,max)\n";
    return 0;
  }
  eemData.reserve(data.size());
  for (unsigned int i=0; i<data.size(); ++i) {
    const std::vector<EtaEtaMassData_t>* src=data[i];
    std::vector<double> *dt=new std::vector<double>();
    dt->reserve(src->size());
    eemData.push_back(dt);
    for (unsigned int ii=0; ii<src->size(); ++ii) {
      dt->push_back((*src)[ii].mass());
    }
  }
  ClearVec(data);
  return 1;
}
#endif

//------------------------------------------------------

#ifdef UseEEM
int ElectronEnergyScale::loadEEMFile(const TString &eemFileName, std::vector<std::vector<EtaEtaMassData_t>*> &eemData, double massMin, double massMax) const {
  EtaEtaMassData_t *eem = new EtaEtaMassData_t();
  int res=1;
  int etaEtaCount = this->numberOfEtaEtaBins();
  eemData.clear(); eemData.reserve(etaEtaCount+1);

  // Read file twice. First time we will determine the number of masses
  // in each (eta1,eta2) bin and allocate memory. The second time we will
  // store the mass value
  const int optimize=1; // optimize memory
  for (int loop=0; res && (loop<1+optimize); ++loop) {
    //std::cout << "loop=" << loop << std::endl;
    std::vector<int> counts(etaEtaCount);
    const char *fname=eemFileName.Data();
    TFile *fin = new TFile(fname);
    assert(fin);
    TTree *tree = (TTree*)fin->Get("Data"); assert(tree);
    tree->SetBranchAddress("Data",&eem);
    TBranch *branch = tree->GetBranch("Data"); assert(branch);
    UInt_t entryMax = tree->GetEntries();
    for (UInt_t ientry=0; ientry < entryMax; ++ientry) {
      branch->GetEntry(ientry);
      int idx=this->getEtaEtaIdx(eem->eta1(),eem->eta2());
      if ((eem->mass()<massMin) || (eem->mass()>massMax)) continue;
      //std::cout << loop << " got " << (*eem) << ", idx=" << idx << "\n";
      if ((idx>=0) && (idx < etaEtaCount)) {
	if (optimize && (loop==0)) {
	  // loops 0 -- count the number of events
	  counts[idx]++;
	}
	else {
	  // loops 1 -- store the memory
	  std::vector<EtaEtaMassData_t>* store=eemData[idx];
	  store->push_back(*eem);
	}
      }
    }
    //delete branch;
    //delete tree;
    delete fin;
    if (optimize && (loop==0)) {
      // allocate memory
      unsigned int k=0;
      for (int i=0; i<this->numberOfEtaBins(); ++i) {
	for (int j=i; j<this->numberOfEtaBins(); ++j, ++k) {
	  std::vector<EtaEtaMassData_t>* tmp=new std::vector<EtaEtaMassData_t>();
	  assert(tmp);
	  //std::cout << "counts[k]=" << counts[k] << "\n";
	  tmp->reserve(counts[k]);
	  eemData.push_back(tmp);
	}
      }
    }
  }
  return 1;
}
#endif
// ------------------------------------------------------------
// ------------------------------------------------------------

#ifdef ZEE_DATA_HH
int ElectronEnergyScale::ProcessZeeDataFile(const char *mc_file_name, const char *data_file_name, std::vector<std::vector<double>*> &mcData, std::vector<std::vector<double>*> &expData, unsigned int debugLoad) {
  ZeeData_t *zee = new ZeeData_t();
  int res=1;
  const int etaEtaCount=this->numberOfEtaEtaBins();

  mcData.clear(); mcData.reserve(etaEtaCount+1);
  expData.clear(); expData.reserve(etaEtaCount+1);
  // Read each file twice. First time we will determine the number of masses
  // in each (eta1,eta2) bin and allocate memory. The second time we will
  // store the mass value
  const int optimize=1; // optimize memory
  const int useScEta=1;
  for (int loop=0; res && (loop<2+2*optimize); ++loop) {
    //std::cout << "loop=" << loop << std::endl;
    std::vector<int> counts(etaEtaCount);
    const int isMC=(loop%2==0) ? 1:0;
    const char *fname=(isMC) ? mc_file_name : data_file_name;
    if (!fname) continue;
    std::vector<std::vector<double>*> *data= (isMC) ? &mcData : &expData;
    TFile *fin = new TFile(fname);
    assert(fin);
    TTree *tree = (TTree*)fin->Get("Events"); assert(tree);
    tree->SetBranchAddress("Events",&zee);
    TBranch *branch = tree->GetBranch("Events");
    UInt_t entryMax = tree->GetEntries();
    for (UInt_t ientry=0; ientry < entryMax; ++ientry) {
      if (debugLoad && (debugLoad<ientry)) continue;
      branch->GetEntry(ientry);
      if (fabs(zee->mass-90)>35) continue;
      int idx=(useScEta) ? 
	this->getEtaEtaIdx(zee->scEta_1,zee->scEta_2) :
	this->getEtaEtaIdx(zee->eta_1,zee->eta_2);
      //std::cout << loop << " got " << (*zee) << ", idx=" << idx << "\n";
      if ((idx>=0) && (idx < etaEtaCount)) {
	if (optimize && (loop/2==0)) {
	  counts[idx]++;
	}
	else {
	  //std::cout << "isMC=" << isMC << ", ientry=" << ientry << ", weights(" << zee->weight << "," << zee->fewz_weight << "," << double(1.) << "), mass=" << zee->mass << "\n";

	  std::vector<double>* store=(*data)[idx];
	  //if (store->capacity()<store->size()+5) 
	  //std::cout << "adding to idx=" << idx << ", count/capacity=" << store->size() << "/" << store->capacity() << std::endl;
	  store->push_back(zee->mass);
	}
      }
    }
    //delete branch;
    //delete tree;
    delete fin;
    if (optimize && (loop/2==0)) {
      // allocate memory
      unsigned int k=0;
      for (int i=0; i<this->numberOfEtaBins(); ++i) {
	for (int j=i; j<this->numberOfEtaBins(); ++j, ++k) {
	  std::vector<double>* tmp=new std::vector<double>();
	  assert(tmp);
	  //std::cout << "counts[k]=" << counts[k] << "\n";
	  tmp->reserve(counts[k]+1);
	  data->push_back(tmp);
	}
      }
    }
  }
  return 1;
}
#endif

// ------------------------------------------------------------

#ifdef ZEE_DATA_HH
int ElectronEnergyScale::ProcessZeeDataFile(const char *mc_file_name, const char *data_file_name, std::vector<std::vector<ZeeData_t*>*> &mcData, std::vector<std::vector<ZeeData_t*>*> &expData, double mass_min, double mass_max, int load_only_idx, unsigned int debugLoad) {
  ZeeData_t *zee = new ZeeData_t();
  int res=1;
  const int etaEtaCount=this->numberOfEtaEtaBins();
  mcData.clear(); mcData.reserve(etaEtaCount+1);
  expData.clear(); expData.reserve(etaEtaCount+1);
  // Read each file twice. First time we will determine the number of masses
  // in each (eta1,eta2) bin and allocate memory. The second time we will
  // store the mass value
  const int optimize=1; // optimize memory
  const int useScEta=1;
  for (int loop=0; res && (loop<2+2*optimize); ++loop) {
    //std::cout << "loop=" << loop << std::endl;
    std::vector<int> counts(etaEtaCount);
    const int isMC=(loop%2==0) ? 1:0;
    const char *fname=(isMC) ? mc_file_name : data_file_name;
    if (!fname) continue;
    std::vector<std::vector<ZeeData_t*>*> *data= (isMC) ? &mcData : &expData;
    TFile *fin = new TFile(fname);
    assert(fin);
    TTree *tree = (TTree*)fin->Get("Events"); assert(tree);
    tree->SetBranchAddress("Events",&zee);
    TBranch *branch = tree->GetBranch("Events");
    UInt_t entryMax = tree->GetEntries();
    for (UInt_t ientry=0; ientry < entryMax; ++ientry) {
      if (debugLoad && (debugLoad<ientry)) continue;
      branch->GetEntry(ientry);
      if ((zee->mass<mass_min) || (zee->mass>mass_max)) continue;
      int idx=(useScEta) ? 
	this->getEtaEtaIdx(zee->scEta_1,zee->scEta_2) :
	this->getEtaEtaIdx(zee->eta_1,zee->eta_2);
      //std::cout << loop << " got " << (*zee) << ", idx=" << idx << "\n";
      if ((load_only_idx!=-1) && (idx!=load_only_idx)) continue;
      if ((idx>=0) && (idx < etaEtaCount)) {
	if (optimize && (loop/2==0)) {
	  counts[idx]++;
	}
	else {
	  //std::cout << "isMC=" << isMC << ", ientry=" << ientry << ", weights(" << zee->weight << "," << zee->fewz_weight << "," << double(1.) << "), mass=" << zee->mass << "\n";

	  std::vector<ZeeData_t*>* store=(*data)[idx];
	  //if (store->capacity()<store->size()+5) 
	  //std::cout << "adding to idx=" << idx << ", count/capacity=" << store->size() << "/" << store->capacity() << std::endl;
	  ZeeData_t *z=new ZeeData_t();
	  z->assign(*zee);
	  store->push_back(z);
	}
      }
    }
    //delete branch;
    //delete tree;
    delete fin;
    if (optimize && (loop/2==0)) {
      // allocate memory
      unsigned int k=0;
      for (int i=0; i<this->numberOfEtaBins(); ++i) {
	for (int j=i; j<this->numberOfEtaBins(); ++j, ++k) {
	  std::vector<ZeeData_t*>* tmp=new std::vector<ZeeData_t*>();
	  assert(tmp);
	  //std::cout << "counts[k]=" << counts[k] << "\n";
	  tmp->reserve(counts[k]+1);
	  data->push_back(tmp);
	}
      }
    }
  }
  return 1;
}
#endif

// ------------------------------------------------------------

#ifdef ZEE_DATA_HH
int ElectronEnergyScale::ProcessZeeDataFiles(const std::vector<std::string> &mc_file_names, const char *data_file_name, std::vector<std::vector<double>*> &mcData, std::vector<std::vector<double>*> &expData, unsigned int debugLoad) {
  mcData.clear(); expData.clear();
  int res=1;
  if (mc_file_names.size()) {
    res=ProcessZeeDataFile(mc_file_names[0].c_str(),data_file_name,mcData,expData,debugLoad);
  }
  else {
    res=ProcessZeeDataFile(NULL,data_file_name,mcData,expData);
  }
  if (res) {
    std::vector<std::vector<double>*> expData_dummy, mcData_tmp;
    for (unsigned int i=1; res && (i<mc_file_names.size()); ++i) {
      res=ProcessZeeDataFile(mc_file_names[i].c_str(),NULL,mcData_tmp,expData_dummy,debugLoad);
      if (res) {
	for (unsigned int k=0; k<mcData_tmp.size(); k++) {
	  std::vector<double> *dt=mcData_tmp[k];
	  std::vector<double> *dest=mcData[k];
	  if (dt->size()) {
	    dest->reserve(dest->size() + dt->size());
	    dest->insert(dest->end(),dt->begin(),dt->end());
	  }
	  delete dt;
	}
      }
    }
  }
  if (!res) std::cout << "error in ProcessZeeDataFiles\n";
  return res;
}
#endif

// ------------------------------------------------------------

#ifdef ZEE_DATA_HH
int ElectronEnergyScale::ProcessZeeDataFileApproximateMCWeight(const char *mc_file_name, const char *data_file_name, std::vector<std::vector<double>*> &mcData, std::vector<std::vector<double>*> &expData, double unitWeight, int include_fewz_weight, int performPUReweight, unsigned int debugLoad) {
  std::cout << "entered ProcessZeeDataFileApproximateMCWeight(" 
	    << ((mc_file_name) ? mc_file_name : "NULL")
	    << "," 
	    << ((data_file_name) ? data_file_name : "NULL")
	    << ")"
	    << "\t includeFEWZweight=" << include_fewz_weight
	    << "\t mcPUReweight=" << performPUReweight
	    << std::endl;

  PUReweight_t puWeight;
  /*
  if (performPUReweight) {
    assert(puWeight.setDefaultFile("DY_m10+pr+a05+o03+pr_4839pb","", 0));
    assert(puWeight.setReference("hNGoodPV_data"));
    assert(puWeight.setActiveSample("hNGoodPV_zee"));
  }
  */

  ZeeData_t *zee = new ZeeData_t();
  int res=1;
  int debug=3; // the distributions
  int debug3testMC=1; 
  int debug2=0; // adding up remaining weights
  TH1D *hraw=0; //(debug==0) ? 0 : new TH1D("hraw","hraw",60,60.,120.);
  TH1D *h1=0; //(debug==0) ? 0 : new TH1D("h1","h1",60,60.,120.);
  TH1D *h2=0; //(debug==0) ? 0 : new TH1D("h2","h2",60,60.,120.);
  if (h1) { h1->SetLineColor(kGreen); h1->SetMarkerColor(kGreen); }
  if (h2) { h2->SetLineColor(kRed);  h2->SetMarkerColor(kRed) ; }
  //if (debug) { assert(hraw); assert(h1); assert(h2); }
  std::vector<TH1D*> hrawV, hfinV;
  const int etaEtaCount=this->numberOfEtaEtaBins();

  if (debug==3) {
    int etaCount=this->numberOfEtaBins();
    hrawV.reserve(etaEtaCount);
    hfinV.reserve(etaEtaCount);
    TH1D *hbase=new TH1D("hbase","hbase",120,60.,120.);
    hbase->SetDirectory(0);
    hbase->GetXaxis()->SetTitle("m_{ee}");
    hbase->GetYaxis()->SetTitle("counts");
    for (int i=0; i<etaCount; ++i) {
      for (int j=i; j<etaCount; ++j) {
	TString hname=Form("hraw_%d_%d",i,j);
	TH1D *h=(TH1D*)hbase->Clone(hname);
	hbase->SetTitle(hname);
	hrawV.push_back(h);
      }
    }
    for (int i=0; i<etaCount; ++i) {
      for (int j=i; j<etaCount; ++j) {
	TString hname=Form("hfin_%d_%d",i,j);
	TH1D *h=(TH1D*)hbase->Clone(hname);
	hbase->SetTitle(hname);
	hfinV.push_back(h);
      }
    }
    delete hbase;
  }

  mcData.clear(); mcData.reserve(etaEtaCount+1);
  expData.clear(); expData.reserve(etaEtaCount+1);
  std::vector<std::vector<double>*> remMCMasses,remMCWeight;
  remMCMasses.reserve(etaEtaCount+1);
  remMCWeight.reserve(etaEtaCount+1);
  double wsum=0., wcount=0.;

  // Read each file twice. First time we will determine the number of masses
  // in each (eta1,eta2) bin and allocate memory. The second time we will
  // store the mass value
  const int optimize=1; // optimize memory
  const int useScEta=1;
  for (int loop=0; res && (loop<2+2*optimize); ++loop) {
    //std::cout << "loop=" << loop << std::endl;
    std::vector<double> counts(etaEtaCount);
    const int isMC=(loop%2==0) ? 1 : 0;
    const int debug3fillHistos=((debug==3) && (isMC == debug3testMC)) ? 1:0;
    const char *fname=(isMC) ? mc_file_name : data_file_name;
    if (!fname) continue;
    std::vector<std::vector<double>*> *data= (isMC) ? &mcData : &expData;
    TFile *fin = new TFile(fname);
    assert(fin);
    TTree *tree = (TTree*)fin->Get("Events"); assert(tree);
    tree->SetBranchAddress("Events",&zee);
    TBranch *branch = tree->GetBranch("Events"); assert(branch);
    UInt_t entryMax = tree->GetEntries();
    for (UInt_t ientry=0; ientry < entryMax; ++ientry) {
      if ((debugLoad!=0) && (ientry>debugLoad)) continue;
      branch->GetEntry(ientry);
      if (fabs(zee->mass-90)>35) continue;
      int idx=(useScEta) ? 
	this->getEtaEtaIdx(zee->scEta_1,zee->scEta_2) :
	this->getEtaEtaIdx(zee->eta_1,zee->eta_2);
      //std::cout << loop << " got " << (*zee) << ", idx=" << idx << "\n";
      if ((idx>=0) && (idx < etaEtaCount)) {
	if (optimize && (loop/2==0)) {
	  double weight=zee->weight;
	  if (isMC && include_fewz_weight) weight*= zee->fewz_weight;
	  double puW=1.0;
	  if (performPUReweight) {
	    //std::cout << "check contents of zee->nPV for getWeightHildreth\n";
	    puW=puWeight.getWeightHildreth( zee->nPV );
	    weight *= puW;
	  }
	  //std::cout << "isMC=" << isMC << ", ientry=" << ientry << ", weights(" << zee->weight << "," << zee->fewz_weight << "," << puW << "), mass=" << zee->mass << "\n";
	  //weight=1.; // debug algo
	  counts[idx]+= weight;
	  if (debug3fillHistos && hrawV.size()) hrawV[idx]->Fill(zee->mass,weight);
	  if ((loop==0) && (fabs(zee->mass-90)<=30)) {
	    wsum+=weight; wcount+=1.;
	    if (fabs(unitWeight)>1e-3) { wcount=1; wsum=unitWeight; }
	  }
	}
	else {
	  std::vector<double>* store=(*data)[idx];
	  //if (store->capacity()<store->size()+5) 
	  //std::cout << "adding to idx=" << idx << ", count/capacity=" << store->size() << "/" << store->capacity() << std::endl;
	  if (!isMC) { 
	    // data - no change
	    store->push_back(zee->mass);
	    if (debug3fillHistos && hfinV.size()) hfinV[idx]->Fill(zee->mass,1.);
	  }
	  else {
	    double w=zee->weight*(wcount/wsum);
	    if (isMC && include_fewz_weight) w *= zee->fewz_weight;
	    if (performPUReweight) w *= puWeight.getWeightHildreth( zee->nPV );
	    if (debug && (ientry<100) && (zee->mass>=60) && (zee->mass<=120)) std::cout << "mass=" << zee->mass << ", w_raw=" << zee->weight << ", w=" << w << "\n";
	    //w=1.; // debug algo
	    if (hraw) hraw->Fill(zee->mass,w);
	    while (w>=1.0) { 
	      store->push_back(zee->mass);
	      if (debug3fillHistos && hfinV.size()) hfinV[idx]->Fill(zee->mass,1.);
	      w-=1.0;
	      if (h1) h1->Fill(zee->mass,1.);
	      if (h2) h2->Fill(zee->mass,1.);
	    }
	    if (fabs(w)>1e-6) {
	      if (!remMCMasses[idx] || !remMCWeight[idx]) { std::cout << "Null ptrs" << std::endl; assert(0); }
	      remMCMasses[idx]->push_back(zee->mass);
	      remMCWeight[idx]->push_back(w);
	    }
	  }
	}
      }
    }
    if (debug) std::cout << "wcount=" << wcount << ", wsum=" << wsum << ", wcount/wsum=" << wcount/wsum << "\n";
   
    //delete branch;
    //delete tree;
    delete fin;
    if (optimize && (loop/2==0)) {
      // allocate memory
      unsigned int k=0;
      for (int i=0; i<this->numberOfEtaBins(); ++i) {
	for (int j=i; j<this->numberOfEtaBins(); ++j, ++k) {
	  std::vector<double>* tmp=new std::vector<double>();
	  assert(tmp);
	  //std::cout << "counts[k]=" << counts[k] << "\n";
	  tmp->reserve(int(counts[k])+1);
	  data->push_back(tmp);
	}
      }
      if (loop==0) {
	k=0;
	for (int i=0; i<this->numberOfEtaBins(); ++i) {
	  for (int j=i; j<this->numberOfEtaBins(); ++j, ++k) {
	    std::vector<double>* tmp=new std::vector<double>();	assert(tmp);
	    tmp->reserve(int(counts[k])+1);
	    remMCMasses.push_back(tmp);
	    tmp=new std::vector<double>();  	assert(tmp);
	    tmp->reserve(int(counts[k])+1);
	    remMCWeight.push_back(tmp);
	  }
	}
      }
    }
  }

  // Add up the remaining MC weights
  int idx=0;
  if (debug2) std::cout << "adding up the remaining MC weights" << std::endl;
  const int debug3fillHistosRem=((debug==3) && debug3testMC) ? 1:0;
  for (int i=0; i<this->numberOfEtaBins(); ++i) {
    for (int j=i; j<this->numberOfEtaBins(); ++j, ++idx) {
      if (remMCWeight[idx]->size()==0) continue;
      std::vector<double>* remMasses=remMCMasses[idx];
      std::vector<double>* remWeight=remMCWeight[idx];
      int do_sort=1;
      if (do_sort && !QuickSort(*remMasses,*remWeight)) {
	std::cout << "QuickSort failed\n";
	return 0;
      }
      //int cx=remMasses->size()-1;
      //int ccc=100;
      if (0) {
	for (unsigned int ii=0; ii<remMasses->size(); ++ii) {
	  std::cout << "ii=" << ii << ", m=" << (*remMasses)[ii] << ", w=" << (*remWeight)[ii] << "\n";
	}
	return 0;
      }
      while (remMasses->size()>2) {
	if (0) {
	  int ii_min=int(remMasses->size())-10;
	  if (ii_min<0) ii_min=0;
	  std::cout << "needs fixing:\n";
	  for (unsigned int ii=(unsigned int)(ii_min); ii<remMasses->size(); ++ii) {
	    std::cout <<"ii=" << ii << ", (" << (*remMasses)[ii] << "," << (*remWeight)[ii] << ")\n";
	  }
	}
	unsigned int ci=remMasses->size()-1, cj=ci-1;
	double m1=(*remMasses)[ci], m2=(*remMasses)[cj];
	for (unsigned int ii=cj; (ii<remMasses->size()); --ii) {
	  if (debug2) std::cout << "ci=" << ci << ", m1=" << m1 << "; cj=" << cj << ", m2=" << m2 << "; ii=" << ii << "; mass=" << (*remMasses)[ii] << "\n";
	  if (fabs(m1-m2)>fabs(m1-(*remMasses)[ii])) { cj=ii; m2=(*remMasses)[ii]; }
	  if (fabs(m1-m2)<0.2) { 
	    if (debug2) std::cout << "break\n"; 
	    break; 
	  }
	  if (do_sort && (fabs(m1-(*remMasses)[ii])>0.5)) { 
	    if (debug2) std::cout << "sorted. break\n"; 
	    break; 
	  }
	}
	if (fabs(m1-m2)>0.2) {
	  if ((*remWeight)[ci]>0.6) {
	    mcData[idx]->push_back(m1);
	    if (debug3fillHistosRem && hfinV.size()) hfinV[idx]->Fill(m1,1.);
	    if (h2) h2->Fill(m1,1.);
	  }
	  //cx=1;
	  if (debug2) std::cout << "erasing ci=" << ci << ", " << (*remMasses)[ci] << "\n";
	  remMasses->erase(remMasses->begin()+ci);
	  remWeight->erase(remWeight->begin()+ci);
	  //std::cout << "modified ci=" << ci << ", " << (*remMasses)[ci] << "\n";
	  //return 0;
	}
	else {
	  double w=(*remWeight)[ci]+(*remWeight)[cj];
	  if (w>=1.) {
	    mcData[idx]->push_back(m1);
	    if (debug3fillHistosRem && hfinV.size()) hfinV[idx]->Fill(m1,1.);
	    w-=1.;
	    if (h2) h2->Fill(m1,1.);
	  }
	  (*remWeight)[cj]=w;
	  if (debug2) std::cout << "branch mi-m2<=0.2 : m1=" << m1 << ", m2=" << m2 << "\n";
	  remMasses->erase(remMasses->begin()+ci);
	  remWeight->erase(remWeight->begin()+ci);
	  if (debug2) std::cout <<"  modified cj is (" << (*remMasses)[cj] << "," << (*remWeight)[cj] << ")\n";
	  //ccc--; if (ccc<=0) return 0;
	}
      }
    }
  }
  if (debug && (hraw || h1 || h2)) {
    TCanvas *ctest = new TCanvas("ctest","ctest",600,600);
    if (hraw) hraw->Draw();
    if (h1) h1->Draw("LP same");
    if (h2) h2->Draw("LP same");
    ctest->Update();
    std::cout << "enter a symbol...";
    char x; cin >> x;
  }
  if (hrawV.size() && hfinV.size()) {
    std::cout << "ProcessZeeDataFileApproximateMCWeight debug print-out\n";
    std::vector<double> diff(hrawV[0]->GetNbinsX());
    std::vector<double> diffratio(hrawV[0]->GetNbinsX());
    for (unsigned int i=0; i<hrawV.size(); ++i) {
      TH1D *h1r=(TH1D*)hrawV[i]->Clone("h1r");
      TH1D *h2f=(TH1D*)hfinV[i]->Clone("h2f");
      h1r->SetDirectory(0);
      h2f->SetDirectory(0);
      std::cout << "h1r->GetName=" << hrawV[i]->GetName() << std::endl;
      h2f->Scale(h1r->Integral()/h2f->Integral());
      double maxdiff=0., maxdifffrac=0.;
      for (int ibin=1; ibin<=h1r->GetNbinsX(); ++ibin) {
	double d=h1r->GetBinContent(ibin);
	diff[ibin-1]= h1r->GetBinContent(ibin) - h2f->GetBinContent(ibin);
	if (d==0.) d=1.;
	if (fabs(diff[ibin-1])>fabs(maxdiff)) {
	  maxdiff=diff[ibin-1];
	  maxdifffrac=maxdiff/d;
	}
	diffratio[ibin-1]= diff[ibin-1]*100./d;
      }
      if (0) {
	std::cout << "raw : ";
	for (int ibin=1; ibin<=h1r->GetNbinsX(); ++ibin) {
	  double d=h1r->GetBinContent(ibin);
	  std::cout << " " << d;
	}
	std::cout << "\n";
	std::cout << "diff: ";
	for (unsigned int ii=0; ii<diff.size(); ++ii) {
	  std::cout << " " << diff[ii];
	}
	std::cout << "\n";
	std::cout <<"frac: ";
	for (unsigned int ii=0; ii<diffratio.size(); ++ii) {
	std::cout << " " << diffratio[ii];
	}
	std::cout << "\n";
      }
      std::cout << "maxdiff=" << maxdiff << ", maxdifffrac=maxdiff/raw=" << maxdifffrac << "\n";
      delete h1r;
      delete h2f;
    }
    ClearVec(hrawV);
    ClearVec(hfinV);
  }
  return 1;
}
#endif

// ------------------------------------------------------------

#ifdef ZEE_DATA_HH
int ElectronEnergyScale::ProcessZeeDataFileApproximateMCWeight(const std::vector<std::string> &mc_file_names, const char *data_file_name, std::vector<std::vector<double>*> &mcData, std::vector<std::vector<double>*> &expData, double unitWeight, int include_fewz_weight, int puReweight, unsigned int debugLoad) {
  std::cout << "entered ProcessZeeDataFileApproximateMCWeight(V)" << std::endl;
  mcData.clear(); expData.clear();
  int res=1;
  if (mc_file_names.size()) {
    res=ProcessZeeDataFileApproximateMCWeight(mc_file_names[0].c_str(),data_file_name,mcData,expData,unitWeight,include_fewz_weight,puReweight,debugLoad);
  }
  else {
    res=ProcessZeeDataFileApproximateMCWeight(NULL,data_file_name,mcData,expData,unitWeight,include_fewz_weight,puReweight,debugLoad);
  }
  if (res) {
    std::vector<std::vector<double>*> expData_dummy, mcData_tmp;
    for (unsigned int i=1; res && (i<mc_file_names.size()); ++i) {
      res=ProcessZeeDataFileApproximateMCWeight(mc_file_names[i].c_str(),NULL,mcData_tmp,expData_dummy,unitWeight,include_fewz_weight,puReweight,debugLoad);
      if (res) {
	for (unsigned int k=0; k<mcData_tmp.size(); k++) {
	  std::vector<double> *dt=mcData_tmp[k];
	  std::vector<double> *dest=mcData[k];
	  if (dt->size()) {
	    dest->reserve(dest->size() + dt->size());
	    dest->insert(dest->end(),dt->begin(),dt->end());
	  }
	  delete dt;
	}
      }
    }
  }
  if (!res) std::cout << "error in ProcessZeeDataFilesApproximateMCWeight(V)\n";
  return res;
}
#endif


// ------------------------------------------------------------
// ------------------------------------------------------------


