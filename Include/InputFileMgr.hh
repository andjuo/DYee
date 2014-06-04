#ifndef InputFileMgr_HH
#define InputFileMgr_HH

#include <TROOT.h>
#include <TString.h>
#include <TObject.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
#include "../Include/BaseClass.hh"
#include "../Include/MyTools.hh"
#include "../Include/DYTools.hh"

// --------------------------------------------------------
// --------------------------------------------------------

typedef enum { _const_none, _const_eff, _const_acceptance, _const_scalefactor } TCorrectionType_t;


// --------------------------------------------------------


class TDescriptiveInfo_t : public TObject {
protected:
  std::vector<std::string> _info;

public:
  TDescriptiveInfo_t() : TObject(), _info() {}
  ~TDescriptiveInfo_t() {}

  unsigned int size() const { return _info.size(); }
  void clear() { _info.clear(); }

  template<class int_type>
  void reserve(int_type reserve_size) { _info.reserve(reserve_size); }

  const std::string get_line(unsigned int idx) const { return _info[idx]; }

  // ----------------

  unsigned int locate(const std::string &key, unsigned int start=0) const {
    const unsigned int m1=(unsigned int)(-1);
    unsigned int idx=m1;
    for (unsigned int i=start; (idx==m1) && (i<_info.size()); ++i) {
      if (PosOk(_info[i],key)) idx=i;
    }
    return idx;
  }

  // ----------------

  void append(const std::string &line) { _info.push_back(line); }

  void print() const {
    std::cout << std::string(80,'-') << "\n";
    std::cout << "File info " << _info.size() << " lines:\n";
    for (unsigned int i=0; i<_info.size(); ++i) {
      std::cout << _info[i];
      if (_info[i].find('\n')==std::string::npos) std::cout << "\n";
    }
    std::cout << std::string(80,'-') << "\n";
  }

  ClassDef(TDescriptiveInfo_t,1)
};


// --------------------------------------------------------

//
// helper class to handle sample inputs
//
class CSample_t
{
public:
  TString          name,label;    // plot item label
  std::vector<Int_t>    colors;    // plot item color
  std::vector<TString>  fnamev;   // ntuple files
  std::vector<Double_t> xsecv;    // per file cross section
  std::vector<TString>  jsonv;    // per file JSON file
  std::vector<Double_t> weightv;  // per file event weight

public:
  CSample_t(){}

  CSample_t(const CSample_t &s) :
    name(s.name), label(s.label),
    colors(s.colors),
    fnamev(s.fnamev),
    xsecv(s.xsecv),
    jsonv(s.jsonv),
    weightv(s.weightv) 
  {}

  ~CSample_t(){ }

  unsigned int size() const { return fnamev.size(); }
  const TString& getName() const { return name; }
  const TString& getLabel() const { return label; }
  template<class int_type>
  TString getFName(int_type idx) const { return fnamev[idx]; }
  template<class int_type>
  TString getJsonFName(int_type idx) const { return jsonv[idx]; }
  template<class int_type>
  double getXsec(int_type idx) const { return xsecv[idx]; }

  TString getFName_nameOnly(int idx) const {
    TString fname=fnamev[idx];
    Ssiz_t pos=fname.Last('/');
    if (pos!=-1) fname=fnamev[idx](pos+1,fname.Length());
    return fname;
  }

  int decodeDefinition(const std::string &line);
  int decodeMCDefinition(const std::string &line);
  int addFile(const std::string &line);

  void clear() {
    name.Clear(); label.Clear(); 
    colors.clear();
    fnamev.clear();
    xsecv.clear();
    jsonv.clear();
    weightv.clear();
  }
};

// --------------------------------------------------------

inline
std::ostream& operator<<(std::ostream& out, const CSample_t &s) {
  out << "   CSample{ name=" << s.name << ", label=" << s.label 
      << ", color[" << s.colors.size() << "]={ ";
  for (unsigned int i=0; i<s.colors.size(); ++i) {
    out << s.colors[i] << " ";
  }
  out << "}, " << s.fnamev.size() << " files }:\n";
  for (unsigned int i=0; i<s.fnamev.size(); i++) {
    out << " " << (i+1) << ". xsec=" << s.xsecv[i] << "  , weight=" << s.weightv[i] << ", fname=<" << s.fnamev[i] << ">, json=<" << s.jsonv[i] << ">\n";
  }
  return out;
}

// --------------------------------------------------------

class InputFileMgr_t : BaseClass_t {
protected:
  TString FLoadedFileName;
  TString FSelectionTag, FConstTag, FXSecTag;
  TString FEnergyScaleTag;
  TString FSpecTagDirUser,FSpecTagFileUser; // auxiliary tags set by user
  TString FSavePlotFormat;
  TString FTriggerTag;
  double FTotLumi;
  int FSelEventsFlag;
  int FFewzCorrFlag;
  int FPUReweightFlag;
  std::vector<TString> FSampleNames;
  std::vector<CSample_t*> FSampleInfos;
  std::vector<TString> FMCSampleNames;
  std::vector<CSample_t*> FMCSignal;
  std::vector<std::string> FUserKeys,FUserValues;

  TString FNtupleDef, FSkimDef;
  TString FRootFileBaseDir; // automatic. Can be changed by user
  TString FNtupleNameExtraTag; // automatic tag
  TString FDirNameExtraTag; // automatic tag
  TDescriptiveInfo_t *FInfoSection; // reserved for TnP section
public:
  InputFileMgr_t() : 
    BaseClass_t("InputFileMgr"),
    FLoadedFileName(),
    FSelectionTag(), FConstTag(), FXSecTag(),
    FEnergyScaleTag(),
    FSpecTagDirUser(),FSpecTagFileUser(),
    FSavePlotFormat(),
    FTriggerTag(),
    FTotLumi(0.), FSelEventsFlag(0),
    FFewzCorrFlag(0),
    FPUReweightFlag(0),
    FSampleNames(), FSampleInfos(),
    FMCSampleNames(), FMCSignal(),
    FUserKeys(), FUserValues(),
    FNtupleDef(), FSkimDef(),
    FRootFileBaseDir(),
    FNtupleNameExtraTag(),
    FDirNameExtraTag(),
    FInfoSection()
  {
    this->autosetRootDir();
  }

  InputFileMgr_t(const InputFileMgr_t &mgr);
  ~InputFileMgr_t();

  void Clear();

  // Access
  const TString& loadedFileName() const { return FLoadedFileName; }
  const TString& selectionTag() const { return FSelectionTag; }
  const TString& constTag() const { return FConstTag; }
  const TString& xsecTag() const { return FXSecTag; }
  const TString& energyScaleTag() const { return FEnergyScaleTag; }
  TString &editEnergyScaleTag() { return FEnergyScaleTag; }
  void clearEnergyScaleTag() { FEnergyScaleTag="UNCORRECTED"; }
  const TString& specTagDirUser() const { return FSpecTagDirUser; }
  const TString& specTagFileUser() const { return FSpecTagFileUser; }
  const TString& savePlotFormat() const { return FSavePlotFormat; }
  const TString& triggerTag() const { return FTriggerTag; }
  double totalLumi() const { return FTotLumi; }
  int selectEventsFlag() const { return FSelEventsFlag; }
  int puReweightFlag() const { return FPUReweightFlag; }
  int fewzFlag() const { return FFewzCorrFlag; }
  const TString& ntupleDef() const { return FNtupleDef; }
  const TString& skimDef() const { return FSkimDef; }
  const TString& ntupleNameExtraTag() const { return FNtupleNameExtraTag; }
  const TString& dirNameExtraTag() const { return FDirNameExtraTag; }
  const TString& rootFileBaseDir() const { return FRootFileBaseDir; }
  // void autosetRootDir(); // see below
  void rootFileBaseDir(TString newdir) { FRootFileBaseDir=newdir; }

  const std::vector<TString>& sampleNames() const { return FSampleNames; }
  const std::vector<CSample_t*>& sampleInfos() const { return FSampleInfos; }
  const std::vector<CSample_t*>& mcSignalInfos() const { return FMCSignal; }

  unsigned int sampleCount() const { return FSampleNames.size(); }

  template<class int_type>
  const TString sampleName(int_type i) const { 
    if (i==(int_type)(-1)) return TString("all");
    return FSampleNames[i]; 
  }

  template<class int_type>
  const CSample_t* sampleInfo(int_type i) const { return FSampleInfos[i]; }
  template<class int_type>
  const CSample_t* mcSampleInfo(int_type i) const { return FMCSignal[i]; }

  int CountDibosonSamples() const;

  unsigned int mcSampleCount() const { return FMCSignal.size(); }
  const std::vector<TString>& mcSampleNames() const { return FMCSampleNames; }
  template<class int_type>
  const TString& mcSampleName(int_type i) const { return FMCSampleNames[i]; }
  template<class int_type>
  const CSample_t* mcInfo(int_type i) const { return FMCSignal[i]; }

  TString dyTag() const {
    int idx=FSelectionTag.Index("/DY");
    TString tag= (idx==-1) ? FSelectionTag : FSelectionTag(idx,FSelectionTag.Length());
    return tag;
  }

  unsigned int userKeyCount() const { return FUserKeys.size(); }
  int userKeyValueExists(const std::string &key) const
  { return (this->userKeyValue(key).size()>0) ? 1:0; }
  std::string userKeyValue(const std::string &key) const;
  TString userKeyValueAsTString(const std::string &key) const { return TString(this->userKeyValue(key).c_str()); }
  int userKeyValueAsInt(const std::string &key) const { return AsInt(this->userKeyValue(key)); }
  double userKeyValueAsDouble(const std::string &key) const { return AsDouble(this->userKeyValue(key)); }

  void addUserKey(const std::string &key, const std::string &value) {
    FUserKeys.push_back(key); FUserValues.push_back(value);
  }

  void addUserKey(const TString &key, const TString &value) {
    FUserKeys.push_back(std::string(key.Data()));
    FUserValues.push_back(std::string(value.Data()));
  }

  void autosetRootDir() {
#ifdef DYee8TeV
    TString dir="../root_files/";
#endif
#ifdef DYee8TeV_reg
    TString dir="../root_files_reg/";
#endif
#ifdef DYee7TeV
    TString dir="../root_files7TeV/";
#endif
    FRootFileBaseDir=dir;
  }

  TString resultBaseDir(const TString &collection,
			//DYTools::TRunMode_t runMode,
			DYTools::TSystematicsStudy_t systMode) const {
    TString dir=FRootFileBaseDir;
    dir.Append(collection);
    dir.Append("/");
    dir.Append(this->dyTag() + 
	       FDirNameExtraTag + 
	       FSpecTagDirUser + 
	       generateSystTag(systMode));
    dir.Append("/");
    return dir;
  }


  TString nTupleDir(//DYTools::TRunMode_t runMode,
		    DYTools::TSystematicsStudy_t systMode) const {
    return resultBaseDir("select",systMode);
  }

  // NtupleNameExtraTag will contain info about E.Scale
  void setNtupleNameExtraTag(const TString &extra) { FNtupleNameExtraTag=extra; }
  void appendNtupleNameExtraTag(const TString &extra) { FNtupleNameExtraTag.Append(extra); }
  void setDirNameExtraTag(const TString &extra) { FDirNameExtraTag=extra; }


  TString resultBaseFileName(const TString &collection, int isNtuple) const {
    TString fname= collection;
    TString extra = DYTools::analysisTag;
    if (!FFewzCorrFlag) extra.Append("_noFEWZ");
    if (!FPUReweightFlag) extra.Append("_noPU");
    if (isNtuple) extra.Append(FNtupleNameExtraTag);
    if (FSpecTagFileUser.Length()) extra.Append(FSpecTagFileUser);
    if (extra.Length()) fname.Append(Form("_%s",extra.Data()));
    fname.Append(".root");
    return fname;
 }
			    

  template<class int_type>
  TString nTupleFullFileName(int_type iSample, DYTools::TSystematicsStudy_t systMode, int unbinned=1) const {
    TString dir=this->nTupleDir(systMode);
    TString ntupleBase=this->sampleName(iSample) + TString("_select");
    TString ntuple= this->resultBaseFileName(ntupleBase,1);
    if (unbinned) {
      TString extra=TString("_") + DYTools::analysisTag;
      ntuple.ReplaceAll(extra,"");
    }
    TString fullName=dir + ntuple;
    return fullName;
  }

  template<class int_type>
  TString yieldFullFileName(int_type iSample, DYTools::TSystematicsStudy_t systMode, int createDir=0) const {
    const int unbinned=0;
    TString fname=this->nTupleFullFileName(iSample,systMode,unbinned);
    fname.ReplaceAll("select","yield");
    if (createDir) CreateDir(fname,1);
    return fname;
  }

  template<class int_type>
  TString yieldFullFileName(int_type iSample, DYTools::TSystematicsStudy_t systMode, const TString &extraTag, int createDir=0) const {
    TString fname=yieldFullFileName(iSample,systMode,createDir);
    int idx=fname.Index(".root");
    fname.Insert(idx,extraTag);
    return fname;
  }

  // signalYieldFullName is based on yield dir
  TString signalYieldFullFileName(DYTools::TSystematicsStudy_t systMode, 
				  int ignoreDebugRunFlag, int createDir=0, 
				  int systematicsFile=0) const {
    TString fname=yieldFullFileName(-1,systMode,createDir);
    if (ignoreDebugRunFlag) fname.ReplaceAll("_DebugRun","");
    const TString target="all_yield";
    Ssiz_t pos=fname.Index(target);
    if (fname.Index(target,pos+1)!=-1) pos=fname.Index(target,pos+1);
    if (pos==-1) { reportError("signalYieldFullName: failed to form a file name"); fname="error.root"; }
    else {
      TString changeTo="bg-subtracted_yield";
      if (systematicsFile) changeTo.Append("Syst");
      fname.Replace(pos,target.Length(),changeTo);
    }
    return fname;
  }

  TString constDir(DYTools::TSystematicsStudy_t systMode, int createDir=0) const {
    TString dir=resultBaseDir("constants",systMode);
    if (createDir) CreateDir(dir,1);
    return dir;
  }

  TString crossSectionDir(DYTools::TSystematicsStudy_t systMode, int createDir=0) const {
    TString dir=resultBaseDir("xsec",systMode);
    if (createDir) CreateDir(dir,1);
    return dir;
  }

  TString crossSectionFullFileName(DYTools::TSystematicsStudy_t systMode,
				   DYTools::TCrossSectionKind_t csKind,
				   int createDir=0,
				   int systematicsFile=0) const {
    if (createDir) TString path= crossSectionDir(systMode,createDir);
    TString fname=this->nTupleFullFileName(-1,systMode,0);
    fname.ReplaceAll("/select/","/xsec/");
    TString target="all_select";
    TString changeTo=TString("xSec_") + CrossSectionKindName(csKind);
    if (systematicsFile) changeTo.Append("Syst");
    Ssiz_t pos=fname.Index(target);
    if (pos!=-1) fname.Replace(pos,target.Length(),changeTo);
    //fname.Prepend(path);
    return fname;
    
  }

  TString theoryDir(DYTools::TSystematicsStudy_t systMode, int createDir=0) const {
    TString dir=resultBaseDir("theory",systMode);
    if (createDir) CreateDir(dir,1);
    return dir;
  }

  TString convertSkim2Ntuple(TString fname) const;

  TString correctionFullFileName(const TString &correctionName, DYTools::TSystematicsStudy_t systMode, int applyNtupleExtraTag) const {
    TString correctionDir=this->constDir(systMode,0);
    //if (systMode!=DYTools::NO_SYST) 
    applyNtupleExtraTag=1;
    TString fname=this->resultBaseFileName(correctionName,applyNtupleExtraTag);
    TString fullName=correctionDir + fname;
    return fullName;
  }

  TString correctionSystFullFileName(const TString &correctionName, DYTools::TSystematicsStudy_t systMode, int applyNtupleExtraTag) const {
    TString corrSystName=correctionName + TString("Syst");
    return correctionFullFileName(corrSystName,systMode,applyNtupleExtraTag);
  }

  // if the user key exist, return the value adjusted for 1D/2D analysis
  // and the number of bins
  int correctionSpecFileName(TString userKey, TString &fileName) const;

  TString theoryFullFileName(const TString &fileNameBase, DYTools::TSystematicsStudy_t systMode, int applyNtupleExtraTag) const {
    TString thDir=this->theoryDir(systMode,0);
    TString fname=this->resultBaseFileName(fileNameBase,applyNtupleExtraTag);
    TString fullName=thDir + fname;
    return fullName;
  }
  

  // TnP section
  const TDescriptiveInfo_t *infoSection() const { return FInfoSection; }
  TString tnpTag() const { return FSelectionTag; }
  //DYTools::TEtBinSet_t etBinsKind() const { return DYTools::ETBINS6; }
  //DYTools::TEtaBinSet_t etaBinsKind() { return DYTools::ETABINS5; }

  TString tnpDir(DYTools::TSystematicsStudy_t systMode, int createDir=0) const { 
    TString dir=resultBaseDir("tag_and_probe",systMode);
    TString userTnpTag=this->userKeyValueAsTString("T&P_PATH_EXTRA");
    if (userTnpTag.Length()) dir.Insert(dir.Length()-1,userTnpTag);
    if (createDir) CreateDir(dir,1);
    return dir;
  }

  TString tnpSelectEventsFName(DYTools::TSystematicsStudy_t systMode, const TString &sampleType, const TString &effType, const TString &triggerName) const {
    const TString u="_";
    const TString b="/";
    TString dir=tnpDir(systMode) + b;
    TString ending=(FPUReweightFlag) ? "_PU.root" : ".root";
    TString fname=Form("selectEvents_%s_%s_%s_%s%s",DYTools::analysisTag.Data(),sampleType.Data(),effType.Data(),triggerName.Data(),ending.Data());
    TString fullName=dir+fname;
    fullName.ReplaceAll("//","/");
    return fullName;
  }

  TString tnpSFSelEventsFullName(DYTools::TSystematicsStudy_t systMode, int createDir=0, int applyNtupleExtraTag=1) const {
    TString dir=this->tnpDir(systMode,createDir);
    applyNtupleExtraTag=1;
    TString fname=this->resultBaseFileName("eventSF_selectedEvents",applyNtupleExtraTag);
    TString rem=TString("_") + DYTools::analysisTag;
    if (fname.Index(rem)!=-1) fname.ReplaceAll(rem,"");
    else fname.ReplaceAll(DYTools::analysisTag,"");
    TString fullName= dir + fname;
    fullName.ReplaceAll("//","/");
    return fullName;
  }

  int tnpEffCalcMethod(DYTools::TDataKind_t, DYTools::TEfficiencyKind_t kind) const; // does not use file info. To be fixed

  TString getTNP_calcMethod(const TDescriptiveInfo_t &info, DYTools::TDataKind_t dataKind, DYTools::TEfficiencyKind_t kind) const;

  std::string getTNP_etetaBinningString(const TDescriptiveInfo_t &info) const;
  TString getTNP_etBinningString(const TDescriptiveInfo_t &info) const;
  TString getTNP_etaBinningString(const TDescriptiveInfo_t &info) const;

  // RECO efficiency needs n-tuples, others can be calculated on skims
  int getTNP_ntuples(const TDescriptiveInfo_t &info, int runOnData,
		     std::vector<TString> &ntupleFileNames,
		     std::vector<TString> &jsonFileNames,
		     int enforceNtupleNames=1) const;

  TString getTNP_calcMethod(DYTools::TDataKind_t dataKind, DYTools::TEfficiencyKind_t kind) const {
    return (FInfoSection) ? this->getTNP_calcMethod(*FInfoSection,dataKind,kind) : TString("getTNP_*: call LoadTnP first");
  }

  std::string getTNP_etetaBinningString() const { return (FInfoSection) ? this->getTNP_etetaBinningString(*FInfoSection) : std::string("getTNP*: call LoadTnP first") ; }
  TString getTNP_etBinningString() const { return (FInfoSection) ? this->getTNP_etBinningString(*FInfoSection) : TString("getTNP*: call LoadTnP first"); }
  TString getTNP_etaBinningString() const {
    return (FInfoSection) ? this->getTNP_etaBinningString(*FInfoSection) : TString("getTNP*: call LoadTnP first");
  }

  int getTNP_ntuples(int runOnData,
		     std::vector<TString> &ntupleFileNames,
		     std::vector<TString> &jsonFileNames) const {
    if (FInfoSection) return getTNP_ntuples(*FInfoSection,runOnData,ntupleFileNames,jsonFileNames);
    return this->reportError("getTNP*: call LoadTnP first");
  }


  /*
  TString constFullName(TCorrectionType_t corr, DYTools::TSystematicsStudy_t systMode, int createDir=0) const {
    TString fname=this->nTupleFullName(-1,systMode);
    fname.ReplaceAll("select","constants");
    fname.ReplaceAll("all_");
    if (createDir) {
      Ssiz_t len=fname.Last('/');
      TString dir=fname(0,len);
      std::cout << "dir=" << dir << "\n";
      gSystem->mkdir(dir,kTRUE);
    }
    return fname;
  }
  */

  TString GetDDBkgFName(const TString fnameBase_inp) const {
    TString tag=this->dyTag();
    TString dir=this->resultBaseDir("ddbkgYield",DYTools::NO_SYST);
    TString fnameBase=fnameBase_inp;
    std::string ver=this->userKeyValue("DDBkgVersion");
    if (ver.size()) fnameBase.Append(TString("_") + TString(ver.c_str()));
    TString fname=dir + fnameBase + 
      TString(Form("_%dD.root",DYTools::study2D+1));
    return fname;
  }

  TString GetTrue2eDDBkgFName() const {
    TString fname=this->GetDDBkgFName("true2eBkgDataPoints");
    return fname;
  }

  TString GetFakeDDBkgFName() const {
    TString fname=this->GetDDBkgFName("fakeBkgDataPoints");
    //TString fname=Form("../root_files/ddbkgYield/%s/fakeBkgDataPoints_%s.root",tag.Data(),DYTools::study2Dstr.Data());
    return fname;
  }


  // Load
  int Load(const TString &inputFile, TDescriptiveInfo_t *TnPSection=NULL);

  // -------------

  int LoadTnP(const TString &inputFile) {
    if (FInfoSection) delete FInfoSection;
    FInfoSection=new TDescriptiveInfo_t();
    if (!FInfoSection) return reportError("LoadTnP: failed to create FInfoSection");
    int res=Load(inputFile,FInfoSection);
    if (!res) printError("in LoadTnP");
    return res;
  }

  // -------------

  // keep only 1st and last sample (data and Zee MC)
  int KeepFirstAndLastSample();
  // if idx=-1, the sample list is cleaned
  int KeepOnlySample(unsigned int idx);
  void removeMCSampleInfo();

  void Print() const { std::cout << *this << std::endl; }
  void PrintSampleFileNames(int nameOnly) const;

  // output
  friend std::ostream& operator<<(std::ostream &out, const InputFileMgr_t &m) {
    out << "InputFileMgr(loadedFile=<" << m.FLoadedFileName << ">, ";
    out << "selectionTag=<" << m.FSelectionTag << ">, ";
    out << "constTag=<" << m.FConstTag << ">, ";
    out << "xsecTag=<" << m.FXSecTag << ">, ";
    out << "energyScaleTag=<" << m.FEnergyScaleTag << ">, ";
    out << "specTagDirUser=<" << m.FSpecTagDirUser << ">, ";
    out << "specTagFileUser=<" << m.FSpecTagFileUser << ">, ";
    out << "NtupleNameExtraTag=<" << m.FNtupleNameExtraTag << ">, ";
    out << "dirNameExtraTag=<" << m.FDirNameExtraTag << ">, ";
    out << "savePlotFormat=<" << m.FSavePlotFormat << ">, ";
    out << "totalLumi=" << m.FTotLumi 
	<< ", weightEvents=" << m.FSelEventsFlag << ", ";
    out << " puReweight=" << m.FPUReweightFlag << ", ";
    out << " fewzCorrection=" << m.FFewzCorrFlag << ", ";
    out << " nTupleDef=<" << m.FNtupleDef 
	<< ">, skimDef=<" << m.FSkimDef << ">";
    out << ";\n";
    out << " rootFileBaseDir=<" << m.FRootFileBaseDir << ">\n";
    out << " " << m.FUserKeys.size() << " user keys: \n";
    for (unsigned int i=0; i<m.FUserKeys.size(); i++) {
      out << "  -- key=<" << m.FUserKeys[i] << ">, value=<" << m.FUserValues[i] << ">\n";
    }
    out << m.FSampleNames.size() << " sample names and " << m.FSampleInfos.size() << " sample infos):\n";
    for (unsigned int i=0; i<m.FSampleNames.size(); ++i) {
      out << " " << (i+1) << " sampleName=" << m.FSampleNames[i] << ":\n";
      out << ( *(m.FSampleInfos[i]) ) << "\n";
    }
    out << m.FMCSampleNames.size() << " mc signal names:\n";
    for (unsigned int i=0; i<m.FMCSampleNames.size(); ++i) {
      out << " " << (i+1) << "  " << m.FMCSampleNames[i] << "\n";
    }
    out << m.FMCSignal.size() << " mc signal file infos:\n";
    for (unsigned int i=0; i<m.FMCSignal.size(); ++i) {
      out << " " << (i+1) << "  " << ( * (m.FMCSignal[i]) ) << "\n";
    }
    return out;
  }
};

// --------------------------------------------------------

// --------------------------------------------------------
// 2011-style input file
/*
class TnPInputFileMgr2011_t { 
protected:
  TString FSampleTypeStr;
  TString FEffTypeStr, FCalcMethodStr;
  TString FEtBinsKindStr, FEtaBinsKindStr;
  TString FDirTag;
  std::vector<TString> FFileNames;
public:
  TnPInputFileMgr2011_t() : FSampleTypeStr(), FEffTypeStr(), FCalcMethodStr(),
			    FEtBinsKindStr(), FEtaBinsKindStr(), FDirTag(),
			    FFileNames() {}

  // cleanup
  void clear() {
    FSampleTypeStr.Clear(); FEffTypeStr.Clear(); FCalcMethodStr.Clear();
    FEtBinsKindStr.Clear(); FEtaBinsKindStr.Clear(); FDirTag.Clear();
    FFileNames.clear();
  }

  // access
  const TString& sampleTypeString() const { return FSampleTypeStr; }
  const TString& effTypeString() const { return FEffTypeStr; }
  const TString& calcMethodString() const { return FCalcMethodStr; }
  const TString& etBinsKindString() const { return FEtBinsKindStr; }
  const TString& etaBinsKindString() const { return FEtaBinsKindStr; }
  const TString& dirTag() const { return FDirTag; }
  const std::vector<TString> fileNames() const { return FFileNames; }
  unsigned int fileCount() const { return FFileNames.size(); }

  const TString& fileName(unsigned int i) const { return FFileNames[i]; }
  const TString& operator[](unsigned int i) const { return FFileNames[i]; }

  // Access with conversion
#ifdef DYToolsUI_HH
  DYTools::TDataKind_t sampleType() const { return DetermineDataKind(FSampleTypeStr); }
  DYTools::TEfficiencyKind_t effType() const { return DetermineEfficiencyKind(FEffTypeStr); }
  DYTools::TTnPMethod_t calcMethod() const { return DetermineTnPMethod(FCalcMethodStr); }
  DYTools::TEtBinSet_t etBinsKind() const { return DetermineEtBinSet(FEtBinsKindStr); }
  DYTools::TEtaBinSet_t etaBinsKind() const { return DetermineEtaBinSet(FEtaBinsKindStr); }
#endif

  // Load 

  int Load(const TString &inputFile);

  friend std::ostream& operator<< (std::ostream& out, const TnPInputFileMgr2011_t &mgr) {
    out << "EffStudyInputMgr2011 (sampleType=" << mgr.FSampleTypeStr << ", effType=" << mgr.FEffTypeStr << ", calcMethod=" << mgr.FCalcMethodStr << ", etBinning=" << mgr.FEtBinsKindStr << ", etaBinning=" << mgr.FEtaBinsKindStr << ", dirTag=<" << mgr.FDirTag << ">, " << mgr.FFileNames.size() << " ntuple files):\n";
    for (unsigned int i=0; i<mgr.FFileNames.size(); ++i) {
      out << (i+1) << " -  " << mgr.FFileNames[i] << "\n";
    }
    return out;
  }
};

// --------------------------------------------------------

class TnPInputFileMgr_t {
protected:
  TString FSampleTypeStr;
  std::vector<TString> FEffTypeStrV, FCalcMethodStrV;
  TString FEtBinsKindStr, FEtaBinsKindStr;
  TString FDirTag;
  std::vector<TString> FFileNames;
public:
  TnPInputFileMgr_t() : FSampleTypeStr(), FEffTypeStrV(), FCalcMethodStrV(),
			 FEtBinsKindStr(), FEtaBinsKindStr(), FDirTag(),
			 FFileNames() {}

  // cleanup
  void clear() {
    FSampleTypeStr.Clear(); FEffTypeStrV.clear(); FCalcMethodStrV.clear();
    FEtBinsKindStr.Clear(); FEtaBinsKindStr.Clear(); FDirTag.Clear();
    FFileNames.clear();
  }

  // access
  const TString& sampleTypeString() const { return FSampleTypeStr; }
  template<class idx_t>
  const TString& effTypeString(idx_t idx) const { return FEffTypeStrV[idx]; }
  template<class idx_t>
  const TString& calcMethodString(idx_t idx) const { return FCalcMethodStrV[idx]; }
  const TString& etBinsKindString() const { return FEtBinsKindStr; }
  const TString& etaBinsKindString() const { return FEtaBinsKindStr; }
  const TString& dirTag() const { return FDirTag; }
  const std::vector<TString> fileNames() const { return FFileNames; }
  unsigned int fileCount() const { return FFileNames.size(); }

  const TString& fileName(unsigned int i) const { return FFileNames[i]; }
  const TString& operator[](unsigned int i) const { return FFileNames[i]; }

  bool hasSameBinCounts(const TnPInputFileMgr_t &mgr) const;

  // Access with conversion
#ifdef DYToolsUI_HH
  DYTools::TDataKind_t sampleType() const { return DetermineDataKind(FSampleTypeStr); }
  template<class idx_t>
  DYTools::TEfficiencyKind_t effType(idx_t idx) const { return DetermineEfficiencyKind(FEffTypeStrV[idx]); }

  template<class idx_t>
  DYTools::TTnPMethod_t effCalcMethod(idx_t idx) const { 
    if (DYTools::efficiencyIsHLT(DYTools::TEfficiencyKind_t(idx))) {
      idx=idx_t(DYTools::HLT);
    }
    return DetermineTnPMethod(FCalcMethodStrV[idx]); 
  }

  int etBinsCount() const { return DYTools::getNEtBins(int(this->etBinsKind())); }
  DYTools::TEtBinSet_t etBinsKind() const { return DetermineEtBinSet(FEtBinsKindStr); }
  DYTools::TEtaBinSet_t etaBinsKind() const { return DetermineEtaBinSet(FEtaBinsKindStr); }
#endif

  // Load 

  int Load(const TString &inputFile);

  friend std::ostream& operator<< (std::ostream& out, const TnPInputFileMgr_t &mgr) {
    out << "EffStudyInputMgr: sampleType=" << mgr.FSampleTypeStr;
    out << ", etBinning=" << mgr.FEtBinsKindStr 
	<< ", etaBinning=" << mgr.FEtaBinsKindStr 
	<< ", dirTag=<" << mgr.FDirTag << ">"; 
    for (unsigned int i=0; i<mgr.FEffTypeStrV.size(); ++i) {
      out << "; effType=" << mgr.FEffTypeStrV[i] 
	  << ", calcMethod=" << mgr.FCalcMethodStrV[i];
    }
    out << "; " << mgr.FFileNames.size() << " ntuple files:\n";
    for (unsigned int i=0; i<mgr.FFileNames.size(); ++i) {
      out << "   " << mgr.FFileNames[i] << "\n";
    }
    return out;
  }
};

// --------------------------------------------------------
// --------------------------------------------------------

class MCInputFileMgr_t {
protected:
  TString FDirTag, FEScaleTag;
  std::vector<TString> FFileNames,FLabels;
  std::vector<Int_t> FColors,FLineStyles;
  std::vector<Double_t> FXSecs,FLumis;
public:
  MCInputFileMgr_t() : FDirTag(),FEScaleTag("20120101_default"),
		       FFileNames(),FLabels(),FColors(),
		       FLineStyles(),FXSecs(),FLumis() {}

  unsigned int size() const { return FFileNames.size(); }
  const TString& dirTag() const { return FDirTag; }
  const TString& escaleTag() const { return FEScaleTag; }
  const std::vector<TString>& fileNames() const { return FFileNames; }
  const std::vector<TString>& labels() const { return FLabels; }
  const std::vector<Int_t>& colors() const { return FColors; }
  const std::vector<Int_t>& lineStyles() const { return FLineStyles; }
  const std::vector<Double_t>& xsecs() const { return FXSecs; }
  const std::vector<Double_t>& lumis() const { return FLumis; }
  template<class idx_t>
  const TString& fileName(idx_t idx) const { return FFileNames[idx]; }
  template<class idx_t>
  const TString& label(idx_t idx) const { return FLabels[idx]; }
  template<class idx_t>
  Int_t color(idx_t idx) const { return FColors[idx]; }
  template<class idx_t>
  Int_t lineStyle(idx_t idx) const { return FLineStyles[idx]; }
  template<class idx_t>
  Double_t xsec(idx_t idx) const { return FXSecs[idx]; }
  template<class idx_t>
  Double_t lumi(idx_t idx) const { return FLumis[idx]; }
  
  int Load(const TString &inputFileName);

  friend std::ostream& operator<<(std::ostream& out, const MCInputFileMgr_t &m) {
    out << "mcInputFileMgr> (" << m.size() << " items):\n";
    for (unsigned int i=0; i<m.size(); ++i) {
      out << "  #" << i << "  <" << m.FFileNames[i] 
	  << "> label=<" << m.FLabels[i]
	  << "> color=" << m.FColors[i]
	  << " lineStyle=" << m.FLineStyles[i]
	  << " xsec=" << m.FXSecs[i] 
	  << " lumi=" << m.FLumis[i]
	  << "\n";
    }
    std::cout << " dirTag=<" << m.FDirTag << ">"
	      << ", escaleTag=<" << m.FEScaleTag 
	      << ">\n";
    return out;
  }
};

// --------------------------------------------------------
// --------------------------------------------------------

class XSecInputFileMgr_t {
protected:
  double FTotLumi;
  TString FName;
  TString FYieldsTag, FConstTag, FEvtEffScaleTag;
  TString FTrigSet;
public:
  XSecInputFileMgr_t() : FTotLumi(0.), FName(),
			 FYieldsTag(), FConstTag(), 
			 FEvtEffScaleTag(),
			 FTrigSet() {}
  void Clear();

  double totLumi() const { return FTotLumi; }
  const TString& name() const { return FName; }
  const TString& yieldsTag() const { return FYieldsTag; }
  const TString& constTag() const { return FConstTag; }
  const TString& evtEffScaleTag() const { return FEvtEffScaleTag; }
  const TString& triggerSetName() const { return FTrigSet; }
  
  int Load(const TString &fname);

  friend std::ostream& operator<<(std::ostream& out, const XSecInputFileMgr_t &mgr) {
    out << "XSecInputFileMgr:\n";
    out << "  -- name=<" << mgr.name() << ">\n";
    out << "  -- yieldsTag=<" << mgr.yieldsTag() << ">\n";
    out << "  -- constTag=<" << mgr.constTag() << ">\n";
    out << "  -- evtEffScaleTag=<" << mgr.evtEffScaleTag() << ">\n";
    out << "  -- triggerSet=<" << mgr.triggerSetName() << ">";
    out << std::endl;
    return out;
  }
};

// --------------------------------------------------------
// --------------------------------------------------------

class EScaleTagFileMgr_t {
protected:
  std::vector<TString> FEScaleTags;
public:
  EScaleTagFileMgr_t() : FEScaleTags() {}
  unsigned int size() const { return FEScaleTags.size(); }
  const TString& escaleTag(int i) const { return FEScaleTags[i]; }
  const std::vector<TString>& getEScaleTags() const { return FEScaleTags; }

  int Load(const TString &inputFileName);
  
  friend std::ostream& operator<<(std::ostream& out, const EScaleTagFileMgr_t &m) {
    out << "eScaleTagFileMgr> (" << m.size() << " items):\n";
    for (unsigned int i=0; i<m.size(); i++) {
      out << "  #" << i << " <" << m.FEScaleTags[i] << ">\n";
    }
    return out;
  }
};
*/

// --------------------------------------------------------

#endif
