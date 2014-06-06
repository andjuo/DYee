#include <assert.h>
#include <fstream>
#include <sstream>
#include "../Include/InputFileMgr.hh"
//#include "../Include/DYTools.hh"

// a file with a function returning the default configuration file names
#include "../Include/InputFileMgr_defaultFiles.hh"
 
// -----------------------------------------------------------
// -----------------------------------------------------------

int GetKeyValue(const std::string &line, std::string &key, std::string &value, int debug=0) {
  if (debug==2) std::cout << "GetKeyValue(" << line << ", debug=" << debug << ")\n";
  if (!line.size() || !PosOk(line,'=')) {
    if (debug==2) std::cout << "keyword not found\n";
    return 0;
  }
  size_t pos=line.find('=');
  size_t p=line.find_last_not_of(" =",pos);
  key=line.substr(0,p+1);
  for (unsigned int i=0; i<key.size(); ++i) key[i]=toupper(key[i]);
  //std::cout << "key=" << key << "\n";
  p=line.find_first_not_of(" =",pos);
  if (!PosOk(p)) p=line.size();
  size_t plast=line.find('#');
  if (!PosOk(plast)) plast=line.size(); else plast--;
  size_t p1=line.find_last_not_of(" \n\t",plast);
  //std::cout << "pos=" << pos << ", p=" << p << ", p1=" << p1 << ", line.size=" << line.size() << "\n";
  value=line.substr(p,p1-p+1);
  if (debug) std::cout << "detected key=<" << key << ">, value=<" << value << ">\n";
  return 1;
}

// -----------------------------------------------------------

int ApplySubstitutions(std::string &line, const std::vector<std::string> &keys, const std::vector<std::string> &values) {
  const int debug=0;
  for (unsigned int i=0; i<keys.size(); ++i) {
    size_t pos=line.find(keys[i]);
    if (PosOk(pos)) {
      if (debug) std::cout << "found key=<" << keys[i] << "> in <" << line << ">\n";
      line.erase(pos,keys[i].size());
      line.insert(pos,values[i]);
      if (debug) std::cout << "insered value=<" << values[i] << "> in <" << line << ">\n";
    }
  }
  return 1;
}

// -----------------------------------------------------------
// -----------------------------------------------------------

int CSample_t::decodeDefinition(const std::string &line) {
  if (line[0]!='$') {
    std::cout << "CSample_t::decodeDefinition: line does not start with '$'\n";
    return 0;
  }
  std::stringstream ss(line);
  char ch;
  int c1,c2,c3,c4;
  ss >> ch >> name >> c1 >> c2 >> c3 >> c4;
  size_t p=line.find('@');
  if (!PosOk(p)) {
    std::cout << "CSample_t::decodeDefinition: line does not contain '@'\n";
    return 0;
  }
  //size_t px=line.find('#');
  //if (PosOk(px)) px--; else px=line.size();
  size_t px=line.size();
  size_t p1=line.find_last_not_of(" \n\t",px);
  label=line.substr(p+1,p1-p+1);
  colors.reserve(4);
  colors.push_back(c1);
  colors.push_back(c2);
  colors.push_back(c3);
  colors.push_back(c4);
  return 1;
}

// -----------------------------------------------------------

int CSample_t::decodeMCDefinition(const std::string &line) {
  if ((line[0]=='$') || (line[0]=='#')) {
    std::cout << "CSample_t::decodeMCDefinition: line cannot start with '" << line[0] << "'\n";
    return 0;
  }
  std::stringstream ss(line);
  std::string fname;
  double xsec;
  int color;
  std::string tmpLabel;
  ss >> fname >> xsec >> color >> tmpLabel;
  size_t p=line.find('@');
  if (!PosOk(p)) {
    std::cout << "CSample_t::decodeDefinition: line does not contain '@'\n";
    return 0;
  }
  //size_t px=tmpLabel.find('#');
  //if (PosOk(px)) px--; else px=tmpLabel.size();
  size_t px=line.size();
  size_t p1=line.find_last_not_of(" \n\t",px);
  label=line.substr(p+1,p1-p+1);
  name=label;
  name.ReplaceAll(" ","_");
  colors.reserve(4);
  colors.push_back(color);
  colors.push_back(color);
  colors.push_back(color);
  colors.push_back(color);
  fnamev.push_back(TString(fname));
  xsecv.push_back(xsec);
  jsonv.push_back(TString());
  weightv.push_back(1.0);
  return 1;
}

// -----------------------------------------------------------

int CSample_t::addFile(const std::string &line_orig) {
  size_t p=line_orig.find('#');
  std::string line=(PosOk(p)) ? line_orig.substr(0,p-1) : line_orig;
  std::stringstream ss(line);
  std::string fname;
  double xsec=0.;
  std::string json;
  ss >> fname >> xsec;
  if (xsec==double(0.)) ss >> json;
  fnamev.push_back(TString(fname));
  xsecv.push_back(xsec);
  jsonv.push_back(TString(json));
  weightv.push_back(1.);
  return 1;
}

// -----------------------------------------------------------
// -----------------------------------------------------------

InputFileMgr_t::InputFileMgr_t (const InputFileMgr_t &mgr) :
  BaseClass_t("InputFileMgr_t"),
  FLoadedFileName(mgr.FLoadedFileName),
  FSelectionTag(mgr.FSelectionTag),
  FConstTag(mgr.FConstTag),
  FXSecTag(mgr.FXSecTag),
  FEnergyScaleTag(mgr.FEnergyScaleTag),
  FSpecTagDirUser(mgr.FSpecTagDirUser),
  FSpecTagFileUser(mgr.FSpecTagFileUser),
  FSavePlotFormat(mgr.FSavePlotFormat),
  FTriggerTag(mgr.FTriggerTag),
  FTotLumi(mgr.FTotLumi),
  FSelEventsFlag(mgr.FSelEventsFlag),
  FFewzCorrFlag(mgr.FFewzCorrFlag),
  FPUReweightFlag(mgr.FPUReweightFlag),
  FSampleNames(),
  FSampleInfos(),
  FMCSampleNames(),
  FMCSignal(),
  FUserKeys(mgr.FUserKeys),
  FUserValues(mgr.FUserValues),
  FNtupleDef(mgr.FNtupleDef),
  FSkimDef(mgr.FSkimDef),
  FRootFileBaseDir(mgr.FRootFileBaseDir),
  FNtupleNameExtraTag(mgr.FNtupleNameExtraTag),
  FDirNameExtraTag(mgr.FDirNameExtraTag),
  FInfoSection(NULL)
{
  if (mgr.FSampleInfos.size()) {
    FSampleNames=mgr.FSampleNames;
    FSampleInfos.reserve(mgr.FSampleInfos.size());
    for (unsigned int i=0; i<mgr.FSampleInfos.size(); ++i) {
      FSampleInfos.push_back(new CSample_t(*mgr.FSampleInfos[i]));
    }
  }
  if (mgr.FMCSignal.size()) {
    FMCSampleNames=mgr.FMCSampleNames;
    FMCSignal.reserve(mgr.FMCSignal.size());
    for (unsigned int i=0; i<mgr.FMCSignal.size(); ++i) {
      FMCSignal.push_back(new CSample_t(*mgr.FMCSignal[i]));
    }
  }
  if (mgr.FInfoSection) {
    FInfoSection=new TDescriptiveInfo_t(*mgr.FInfoSection);
  }
}


// -----------------------------------------------------------

InputFileMgr_t::~InputFileMgr_t() {
  this->Clear();
}

// -----------------------------------------------------------

void InputFileMgr_t::Clear() {
  FLoadedFileName.Clear(); 
  FSelectionTag.Clear(); FConstTag.Clear(); FXSecTag.Clear();
  FEnergyScaleTag.Clear();
  FSpecTagDirUser.Clear(); FSpecTagFileUser.Clear();
  FSavePlotFormat.Clear();
  FTriggerTag.Clear();
  FTotLumi=0.;
  FSelEventsFlag=1;
  FFewzCorrFlag=1;
  FPUReweightFlag=1;
  FSampleNames.clear();
  ClearVec(FSampleInfos);
  FMCSampleNames.clear();
  ClearVec(FMCSignal);
  FUserKeys.clear(); FUserValues.clear();
  FNtupleDef.Clear(); FSkimDef.Clear();
  FNtupleNameExtraTag.Clear();
  FDirNameExtraTag.Clear();
}

// -----------------------------------------------------------

int InputFileMgr_t::CountDibosonSamples() const {
  int count[3];
  count[0]=0; count[1]=0; count[2]=0;
  for (unsigned int i=0; i<FSampleNames.size(); ++i) {
    TString sname=FSampleNames[i];
    if (sname=="ww") count[0]++;
    if (sname=="wz") count[1]++;
    if (sname=="zz") count[2]++;
  }
  int sum=0;
  for (unsigned int i=0; i<3; ++i) {
    sum+=count[i];
    if (count[i]>1) {
      std::cout << "\n\n\tInputFileMgr::CountDibosonSamples detected input error: count[" << i << "]>1\n\n";
      return -1;
    }
  }
  return sum;
}

// -----------------------------------------------------------

std::string InputFileMgr_t::userKeyValue(const std::string &key_input) const {
  std::string key=key_input;
  std::string value;
  // check the supplied version
  for (unsigned int i=0; i<FUserKeys.size(); ++i) {
    if (FUserKeys[i] == key) value=FUserValues[i];
  }
  // check the capitalized version
  if (!value.size()) {
    for (unsigned int i=0; i<key.size(); ++i) key[i]=toupper(key[i]);
    for (unsigned int i=0; i<FUserKeys.size(); ++i) {
      if (FUserKeys[i] == key) value=FUserValues[i];
    }
  }
  return value;
}

// -----------------------------------------------------------

int InputFileMgr_t::Load(const TString &inputfname_inp,
			 TDescriptiveInfo_t *TnPSection) {
  const int c_Start=0;
  const int c_Defs=c_Start+1;
  const int c_Const=c_Defs+1;
  const int c_TnP=c_Const+1;
  const int c_Data=c_TnP+1;
  const int c_Bkg=c_Data+1;
  const int c_MCfiles=c_Bkg+1;

  TString inputfname=inputfname_inp;
#ifdef InputFileMgr_defaultFile_HH
  inputfname=getDefaultFileName(inputfname_inp);
#endif

  std::ifstream ifs;
  ifs.open(inputfname.Data());
  if (!ifs.is_open()) {
    std::cout << "failed to load input file <" << inputfname << ">\n";
    throw 2;
  }

  std::string line;
  int debug=0; // 0, 1, 2
  int state=c_Start;
  std::string key, value;

  // find definitions
  std::vector<std::string> keys;
  std::vector<std::string> values;
  while((state==c_Start) && !ifs.eof() && getline(ifs,line)) {
    if ((line[0]=='#') || !line.size()) {
      if (debug==2) std::cout << "skipping <" << line << ">\n";
      continue;
    }
    if (PosOk(line,"[DEFINITIONS]")) {
      if (debug) std::cout << "found [DEFINITIONS]\n";
      state++;
    }
  }

    if (state!=1) {
    ifs.close();
    return this->reportError("failed to locate [DEFINITIONS] in %s",inputfname);
  }

  if (debug) std::cout << "\n\tload definitions\n";
  while((state==c_Defs) && !ifs.eof() && getline(ifs,line)) {
    if ((line[0]=='#') || !line.size()) {
      if (debug==2) std::cout << "skipping <" << line << ">\n";
      continue;
    }
    if ((line[0]=='%') || (line[0]=='[')) {
      if (debug==2) std::cout << "end of section (" << line[0] << ")\n";
      state++;
    }
    else if (GetKeyValue(line,key,value,debug)) {
      if (debug) std::cout << "adding key=<" << key << ">, value=<" << value << ">\n";
      if (key=="NTUPLE") { 
	FNtupleDef=value; 
	if (debug==2) std::cout << " * nTupleDef=<" << FNtupleDef << ">\n";
      }
      else if (key=="SKIM") {
	FSkimDef=value;
	if (debug==2) std::cout << " * skimDef=<" << FSkimDef << ">\n";
      }
      else if (key=="RESULT_ROOT_DIR") {
	FRootFileBaseDir=value;
	//if (debug==2)
	std::cout << " *!!* setting RESULT_ROOT_DIR=<" << FRootFileBaseDir << ">\n";
      }
      else {
	keys.push_back(key);
	values.push_back(value);
      }
    }
  }

  //
  // Constants
  //

  if (debug) std::cout << "\n\tload constants\n";
  while((line[0]!='[') && !ifs.eof() && getline(ifs,line)) ;
  if (!PosOk(line,"[CONSTANTS]")) {
    ifs.close();
    std::cout << "expected section '[CONSTANTS]' not found\n";
    return 0;
  }
  while((state==c_Const) && !ifs.eof() && getline(ifs,line)) {
    if ((line[0]=='#') || !line.size()) {
      if (debug==2) std::cout << "skipping <" << line << ">\n";
      continue;
    }
    if ((line[0]=='%') || (line[0]=='[')) {
      if (debug==2) std::cout << "end of section (" << line[0] << ")\n";
      state++;
    }
    else if (GetKeyValue(line,key,value,debug)) {
      if (debug) std::cout << "got key=<" << key << ">, value=<" << value << ">\n";
      if (key=="LUMI") {
	FTotLumi=atof(value.c_str());
	if (debug) std::cout << key << "=<" << value << "> (" << FTotLumi << ")\n";
      }
      else if (key=="SELECTEVENTSFLAG") {
	FSelEventsFlag=atoi(value.c_str());
	if (debug) std::cout << key << "=<" << value << "> (" << FSelEventsFlag << ")\n";
      }
      else if (key=="SELECTIONTAG") {
	FSelectionTag=value;
	if (debug) std::cout << key << "=<" << value << ">\n";
      }
      else if (key=="CONSTANTSTAG") {
	FConstTag=value;
	if (debug) std::cout << key << "=<" << value << ">\n";
      }
      else if (key=="XSECTAG") {
	FXSecTag=value;
	if (debug) std::cout << key << "=<" << value << ">\n";
      }
      else if (key=="ESCALESET") {
	FEnergyScaleTag=value;
	if (debug) std::cout << key << "=<" << value << ">\n";
      }
      else if (key=="PLOTEXTENSION") {
	FSavePlotFormat=value;
	if (debug) std::cout << key << "=<" << value << ">\n";
      }
      else if (key=="TRIGGER") {
	FTriggerTag=value;
	if (debug) std::cout << key << "=<" << value << ">\n";
      }
      else if (key=="FEWZ") {
	FFewzCorrFlag=atoi(value.c_str());
	if (debug) std::cout << key << "=<" << value << "> (" << FFewzCorrFlag << ")\n";
	  }
      else if (key=="PUREWEIGHT") {
	FPUReweightFlag=atoi(value.c_str());
	if (debug) std::cout << key << "=<" << value << "> (" << FPUReweightFlag << ")\n";
      }
      else if (key=="SPECTAGDIRUSER") {
	FSpecTagDirUser=value;
	if (debug==2) std::cout << " * specTagDirUser=<" << FSpecTagDirUser << ">\n";
      }
      else if (key=="SPECTAGFILEUSER") {
	FSpecTagFileUser=value;
	if (debug==2) std::cout << " * specTagFileUser=<" << FSpecTagFileUser << ">\n";
      }
      else {
	FUserKeys.push_back(key);
	FUserValues.push_back(value);
	std::cout << " user key=<" << key << ">, with value=<" << value << ">\n";
      }

    }
  }

  //
  // Tag and probe
  //

  if (debug) std::cout << "\n\tload tag&probe section\n";
  while((line[0]!='[') && !ifs.eof() && getline(ifs,line)) ;
  if (!PosOk(line,"[TAG_AND_PROBE]")) {
    ifs.close();
    std::cout << "expected section '[TAG_AND_PROBE]' not found\n";
    return 0;
  }
  if (TnPSection) { 
    TnPSection->clear();
    TnPSection->reserve(10);
  }
  while((state==c_TnP) && !ifs.eof() && getline(ifs,line)) {
    if ((line[0]=='#') || !line.size()) {
      if (debug==2) std::cout << "skipping <" << line << ">\n";
      continue;
    }
    if ((line[0]=='%') || (line[0]=='[')) {
      if (debug==2) std::cout << "end of section (" << line[0] << ")\n";
      state++;
    }
    else if (TnPSection) {
      if (GetKeyValue(line,key,value,debug)) {
	if (debug) std::cout << "got key=<" << key << ">, value=<" << value << ">\n";
	TnPSection->append(line);
      }
    }
  }

  //
  // Data Ntuples
  //

  CSample_t sample;

  if (debug) std::cout << "\n\tload data ntuples section\n";
  while((line[0]!='[') && !ifs.eof() && getline(ifs,line)) ;
  if (!PosOk(line,"[DATA_FILES]")) {
    ifs.close();
    std::cout << "expected section '[DATA_FILES]' not found\n";
    return 0;
  }
  while((state==c_Data) && !ifs.eof() && getline(ifs,line)) {
    if ((line[0]=='#') || !line.size()) {
      if (debug==2) std::cout << "skipping <" << line << ">\n";
      continue;
    }
    if ((line[0]=='%') || (line[0]=='[')) {
      if (debug==2) std::cout << "end of section (" << line[0] << ")\n";
      state++;
    }
    else if (line[0]=='$') {
      sample.clear();
      if (!sample.decodeDefinition(line)) {
	std::cout << "failed to decode sample definition in <" << line << ">\n";
	ifs.close();
	return 0;
      }
    }
    else {
      ApplySubstitutions(line, keys,values);
      if (!sample.addFile(line)) {
	std::cout << "failed to decode file definition in <" << line << ">\n";
	ifs.close();
	return 0;
      }
    }
  }

  if (debug==2) std::cout << "\nData sample:\n" << sample << "\n";
  FSampleNames.push_back(sample.name);
  FSampleInfos.push_back(new CSample_t(sample));


  //
  // Background files
  //

  sample.clear();
  if (debug) std::cout << "\n\tload background definition section\n";
  while((line[0]!='[') && !ifs.eof() && getline(ifs,line)) ;
  if (!PosOk(line,"[BACKGROUND_FILES]")) {
    ifs.close();
    std::cout << "expected section '[BACKGROUND_FILES]' not found\n";
    return 0;
  }
  while((state==c_Bkg) && !ifs.eof() && getline(ifs,line)) {
    if ((line[0]=='#') || !line.size()) {
      if (debug==2) std::cout << "skipping <" << line << ">\n";
      continue;
    }
    if ((line[0]=='%') || (line[0]=='[')) {
      if (debug==2) std::cout << "end of section (" << line[0] << ")\n";
      state++;
    }
    else if (line[0]=='$') {
      if (sample.size()) {
	FSampleNames.push_back(sample.name);
	FSampleInfos.push_back(new CSample_t(sample));
      }
      sample.clear();
      if (!sample.decodeDefinition(line)) {
	std::cout << "failed to decode sample definition in <" << line << ">\n";
	ifs.close();
	return 0;
      }
    }
    else {
      ApplySubstitutions(line, keys,values);
      if (!sample.addFile(line)) {
	std::cout << "failed to decode file definition in <" << line << ">\n";
	ifs.close();
	return 0;
      }
    }
  }

  if (debug==2) std::cout << "\nLast background sample:\n" << sample << "\n";
  FSampleNames.push_back(sample.name);
  FSampleInfos.push_back(new CSample_t(sample));

  //
  // MC signal
  //

  sample.clear();
  if (debug) std::cout << "\n\tload MC signal section\n";
  while((line[0]!='[') && !ifs.eof() && getline(ifs,line)) ;
  if (!PosOk(line,"[MC_SIGNAL]")) {
    ifs.close();
    std::cout << "expected section '[MC_SIGNAL]' not found\n";
    return 0;
  }
  while((state==c_MCfiles) && !ifs.eof() && getline(ifs,line)) {
    if ((line[0]=='#') || !line.size()) {
      if (debug==2) std::cout << "skipping <" << line << ">\n";
      continue;
    }
    if ((line[0]=='%') || (line[0]=='[')) {
      if (debug==2) std::cout << "end of section (" << line[0] << ")\n";
      state++;
    }
    else if (line[0]=='$') {
      if (sample.size()) {
	FSampleNames.push_back(sample.name);
	FSampleInfos.push_back(new CSample_t(sample));
      }
      sample.clear();
      if (!sample.decodeDefinition(line)) {
	std::cout << "failed to decode sample definition in <" << line << ">\n";
	ifs.close();
	return 0;
      }
    }
    else {
      ApplySubstitutions(line, keys,values);
      CSample_t *mcSample=new CSample_t();
      if (!sample.addFile(line)) {
	std::cout << "failed to decode file definition in <" << line << ">\n";
	ifs.close();
	return 0;
      }
      else {
	if (!mcSample->decodeMCDefinition(line)) {
	  std::cout << "failed to decode mc definition in <" << line << ">\n";
	  ifs.close();
	  return 0;
	}
	FMCSampleNames.push_back(mcSample->getLabel());
	FMCSignal.push_back(mcSample);
      }
    }
  }

  FSampleNames.push_back(sample.name);
  FSampleInfos.push_back(new CSample_t(sample));

  ifs.close();
  
  FLoadedFileName=inputfname;
  return DYTools::checkTotalLumi(FTotLumi);
}

// -----------------------------------------------------------

int InputFileMgr_t::KeepFirstAndLastSample() {
  FSampleNames.erase(FSampleNames.begin()+1,FSampleNames.end()-1);
  for (unsigned int i=1; i<FSampleInfos.size()-1; ++i) {
    CSample_t *s=FSampleInfos[i];
    delete s;
  }
  FSampleInfos.erase(FSampleInfos.begin()+1,FSampleInfos.end()-1);
  return 1;
}

// -----------------------------------------------------------

int InputFileMgr_t::KeepOnlySample(unsigned int idx) {
  int isTail=(idx==FSampleNames.size()-1) ? 1:0;
  int isHead=(idx==0) ? 1:0;
  if (idx==(unsigned int)(-1)) {
    // clean entire list
    isTail=1; 
    idx=(unsigned int)(FSampleNames.end()-FSampleNames.begin());
  }
  if (!isTail) FSampleNames.erase(FSampleNames.begin()+idx+1,FSampleNames.end());
  if (!isHead) FSampleNames.erase(FSampleNames.begin(),FSampleNames.begin()+idx);
  for (unsigned int i=0; i<FSampleInfos.size(); ++i) {
    if (i==idx) continue;
    CSample_t *s=FSampleInfos[i];
    delete s;
  }
  if (!isTail) FSampleInfos.erase(FSampleInfos.begin()+idx+1,FSampleInfos.end());
  if (!isHead) FSampleInfos.erase(FSampleInfos.begin(),FSampleInfos.begin()+idx);
  return 1;
}

// -----------------------------------------------------------

void InputFileMgr_t::removeMCSampleInfo() {
  for (unsigned int i=0; i<FMCSignal.size(); ++i) {
    CSample_t *s=FMCSignal[i];
    delete s;
  }
  FMCSampleNames.clear();
  FMCSignal.clear();
}


// -----------------------------------------------------------

void InputFileMgr_t::PrintSampleFileNames(int nameOnly) const {
  std::cout << "\n";
  for (unsigned int i=0; i<FSampleInfos.size(); ++i) {
    const CSample_t *cs=FSampleInfos[i];
    for (unsigned int ii=0; ii<cs->size(); ++ii) {
      TString name=(nameOnly) ? cs->getFName_nameOnly(ii) : cs->getFName(ii);
      std::cout << name << "\n";
    }
  }
  std::cout << "\n";
}

// -----------------------------------------------------------

TString InputFileMgr_t::convertSkim2Ntuple(TString fname) const {
  if (fname.Index(FSkimDef)!=-1) {
    //std::cout << "modifying <" << fname << ">\n";
    fname.ReplaceAll(FSkimDef,FNtupleDef);
  }
  else if (fname.Index(FNtupleDef)!=-1) {
    fname.ReplaceAll(FNtupleDef,FSkimDef);
  }
  else {
    this->reportError("convertSkim2Ntuple: could not convert <%s>",fname);
  }
  //std::cout << "new name=<" << fname << ">\n";
  return fname;
}

// -----------------------------------------------------------

// if the user key exist, return the value adjusted for 1D/2D analysis
// and the number of bins
int InputFileMgr_t::correctionSpecFileName(TString userKey, TString &fileName) const {
  TString specFile=this->userKeyValueAsTString(userKey.Data());
  if (specFile.Length()) {
    std::cout << "user-defined scale factor file <" << specFile << ">\n";
    fileName=specFile;

    if (DYTools::study2D) {
      Ssiz_t idx=specFile.Last('/');
      TString path=specFile(0,idx+1);
      TString fname=specFile(idx+1,specFile.Length());
      fname.ReplaceAll("1D","2D");
      fname.ReplaceAll("nMB41","nMB7");
      fileName= path+fname;
    }
    return 1;
  }
  // not found
  return 0;
}

// -----------------------------------------------------------

int InputFileMgr_t::tnpEffCalcMethod(DYTools::TDataKind_t dataKind, DYTools::TEfficiencyKind_t effKind) const {
  int method=0;
  switch(dataKind) {
  case DYTools::DATA: {
    switch(effKind) {
    case DYTools::RECO:
    case DYTools::ID:
      method=DYTools::FITnFIT;
      break;
    default:
      method=DYTools::COUNTnCOUNT;
    }
  }
    break;
  case DYTools::MC: {
    method=DYTools::COUNTnCOUNT;
  }
    break;
  default:
    reportError("tnpEffCalcMethod: unhandled dataKind=%d",int(dataKind));
  }
  return method;
}

// -----------------------------------------------------------

TString InputFileMgr_t::getTNP_calcMethod(const TDescriptiveInfo_t &info, DYTools::TDataKind_t dataKind, DYTools::TEfficiencyKind_t kind) const {
  std::string dataStr=DataKindName(dataKind).Data();
  std::string kindStr=EfficiencyKindName(kind).Data();
  if (DYTools::efficiencyIsHLT(kind)) kindStr="HLT";
  std::string key=kindStr + std::string("_") + dataStr;
  for (unsigned int i=0; i<key.size(); i++) key[i]=toupper(key[i]);
  unsigned int idx=info.locate(key);
  std::string key1,value;
  if (idx<info.size()) {
    if (!GetKeyValue(info.get_line(idx),key1,value)) {
      reportError("getTNP_calcMethod","code error");
    }
  }
  return TString(value.c_str());
}

// -----------------------------------------------------------

std::string InputFileMgr_t::getTNP_etetaBinningString(const TDescriptiveInfo_t &info) const {
  std::string key="MAP";
  std::string value, key1;
  unsigned int idx=info.locate(key);
  if ((idx<info.size()) && !GetKeyValue(info.get_line(idx),key1,value)) {
    reportError("getTNP_etetaBinningString","code error");
  }
  return value;
}
// -----------------------------------------------------------

TString InputFileMgr_t::getTNP_etBinningString(const TDescriptiveInfo_t &info) const {
  std::string value=getTNP_etetaBinningString(info);
  size_t p=value.find("ETB");
  size_t p1=value.find(" ",p);
  if (!PosOk(p1)) p1=value.size();
  return TString(value.substr(p,p1-p).c_str());
}


// -----------------------------------------------------------

TString InputFileMgr_t::getTNP_etaBinningString(const TDescriptiveInfo_t &info) const {
  std::string value=getTNP_etetaBinningString(info);
  size_t p=value.find("ETAB");
  size_t p1=value.find(" ",p);
  if (!PosOk(p1)) p1=value.size();
  return TString(value.substr(p,p1-p).c_str());
}

 
// -----------------------------------------------------------

int InputFileMgr_t::getTNP_ntuples(const TDescriptiveInfo_t &info,
				   int runOnData,
				   std::vector<TString> &ntupleFileNames,
				   std::vector<TString> &jsonFileNames,
				   int enforceNtupleNames) const
{
  int res=1;
  if (info.size()) res=1; // dummy to prevent compiler complaints
  ntupleFileNames.clear();
  jsonFileNames.clear();
  // obtain file names
  if (!runOnData) {
    for (unsigned int is=0; is<FMCSignal.size(); ++is) {
      CSample_t *sample=FMCSignal[is];
      for (unsigned int i=0; i<sample->size(); ++i) {
	if ((sample->getFName(i).Index("zeem20to500")!=-1) ||
	    (sample->getFName(i).Index("zeem20to200")!=-1)) {
	  ntupleFileNames.push_back(sample->getFName(i));
	  break;
	}
      }
    }
  }
  else {
    const CSample_t *sample= FSampleInfos[0];
    for (unsigned int i=0; i<sample->size(); ++i) {
      ntupleFileNames.push_back(sample->getFName(i));
      jsonFileNames.push_back(sample->getJsonFName(i));
    }
  }

  // change 'tight-loose_skim' to 'ntuple'
  if (enforceNtupleNames) {
    for (unsigned int i=0; i<ntupleFileNames.size(); ++i) {
      TString fname=ntupleFileNames[i];
      if (fname.Index(FSkimDef)!=-1) {
	//std::cout << "modifying i=" << i << " <" << fname << ">\n";
	fname.ReplaceAll(FSkimDef,FNtupleDef);
	//std::cout << "new name=<" << fname << ">\n";
	ntupleFileNames[i]=fname;
      }
    }
  }

  res=(ntupleFileNames.size()>0) ? 1:0;
  if (!res) res=reportError("getTNP_ntuples: no files found (runOnData=%d)",runOnData);
  return res;
}

 
// -----------------------------------------------------------
// -----------------------------------------------------------
/*
int TnPInputFileMgr2011_t::Load(const TString &configFile) {
  this->clear();
  std::ifstream ifs;
  ifs.open(configFile.Data());
  if (!ifs.is_open()) {
    std::cout << "TnPInputFileMgr2011::Load  tried to open the configuration file <" << configFile << ">\n";
    return 0;
  }
  std::string line;
  int state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') break;
    if(state==0){
      // Read 1st line of content: data or MC?
      FSampleTypeStr = TString(line);
      state++;
    }else if(state==1) {
      // Read 2d content line: efficiency type string
      FEffTypeStr = TString(line);
      state++;
    }else if(state==2) {
      // Read 3d content line: fitting mode
      FCalcMethodStr = TString(line);
      state++;
    }else if(state==3) {
      // Read 4th content line: SC ET binning
      FEtBinsKindStr = TString(line);
      state++;
    }else if(state==4) {
      // Read 5th content line: SC eta binning
      FEtaBinsKindStr = TString(line);
      state++;
    }else if(state==5) {
      // Read 5th content line: SC eta binning
      FDirTag = TString(line);
      state++;
    }else if(state==6) {
      FFileNames.push_back(TString(line));
    }
  }
  ifs.close();
  int res=((state==6) && FFileNames.size()) ? 1:0;
  if (!res) {
    std::cout << "Failed to load file <" << configFile << ">\n";
  }
  return res;
}

// -----------------------------------------------------------
// -----------------------------------------------------------

int TnPInputFileMgr_t::Load(const TString &configFile) {
  this->clear();
  std::ifstream ifs;
  ifs.open(configFile.Data());
  if (!ifs.is_open()) {
    std::cout << "TnPInputFileMgr::Load tried to open configFile=<" << configFile << ">\n";
    return 0;
  }
  string line;
  Int_t state=0;
  Int_t subState=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') break;
    if(state==0){
      // Read 1st line of content: data or MC?
      FSampleTypeStr = TString(line);
      state++;
    }else if(state==1) {
      // Read 2d content line: efficiency fitting mode
      size_t pos=line.find(':');
      if (pos==string::npos) {
        std::cout << "expected format is EFFICIENCY:fitting_mode\n";
        return 0;
      }
      subState++;
      int orderErr=0;
      switch(subState) {
      case 1: if (line.find("RECO")==string::npos) orderErr=1; break;
      case 2: if (line.find("ID"  )==string::npos) orderErr=1; break;
      case 3: if (line.find("HLT" )==string::npos) orderErr=1; break;
      default:
	std::cout << "incorrect value of subState=" << subState << "\n";
	return 0;
      }
      if (orderErr) {
	std::cout << "EfficiencyKind:CalculationMethod should be ordered RECO,ID,HLT\n";
	return 0;
      }
      FEffTypeStrV.push_back(TString(line.substr(0,pos)));
      FCalcMethodStrV.push_back(TString(line.c_str()+pos+1));
      if (subState==3) state++;
    }else if(state==2) {
      // Read 3rd content line: SC ET binning
      FEtBinsKindStr = TString(line);
      state++;
    }else if(state==3) {
      // Read 4th content line: SC eta binning
      FEtaBinsKindStr = TString(line);
      state++;
    }else if(state==4) {
      // Read 5th content line: directory tag
      FDirTag = TString(line);
      state++;
    }else if(state==5) {
      FFileNames.push_back(TString(line));
    }
  }
  ifs.close();
  if (state!=5) {
    std::cout << "Load: did not reach state=5\n";
    return 0;
  }

  //std::cout << "Loaded:\n" << *this << "\n";
  return 1;
}

// -----------------------------------------------------------

bool TnPInputFileMgr_t::hasSameBinCounts(const TnPInputFileMgr_t &mgr) const {
  return ((FEtBinsKindStr==mgr.FEtBinsKindStr) &&
	  (FEtaBinsKindStr==mgr.FEtaBinsKindStr) &&
	  (FEffTypeStrV.size() == mgr.FEffTypeStrV.size()) &&
	  (FCalcMethodStrV.size() == mgr.FCalcMethodStrV.size())) ? kTRUE : kFALSE;
}

// -----------------------------------------------------------
// -----------------------------------------------------------

int MCInputFileMgr_t::Load(const TString& inputFileName) {
  std::ifstream ifs;
  ifs.open(inputFileName.Data());
  if (!ifs.is_open()) {
    std::cout << "MCInputFileMgr: failed to open file <" << inputFileName << ">\n";
    return 0;
  }
  std::string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') break;
    if(state == 0){
      FDirTag = TString(line);
      getline(ifs,line);
      stringstream ss3(line); ss3 >> FEScaleTag;
      state++;
      continue;
    }else{
      std::string fname;
      Int_t color1, linesty;
      std::stringstream ss(line);
      Double_t xsec1;
      ss >> fname >> xsec1 >> color1 >> linesty;
      string label1 = line.substr(line.find('@')+1);
      if (label1.size()) {
	FFileNames.push_back(fname);
	FLabels.push_back(label1);
	FColors.push_back(color1);
	FLineStyles.push_back(linesty);
	FXSecs.push_back(xsec1);
	FLumis.push_back(0);
      }
    }
  }
  ifs.close();
  std::cout << "Loaded:\n" << *this << "\n";
  return FFileNames.size();
}

// -----------------------------------------------------------
// -----------------------------------------------------------

void XSecInputFileMgr_t::Clear() {
  FTotLumi=0;
  FName.Clear();
  FYieldsTag.Clear(); FConstTag.Clear(); FEvtEffScaleTag.Clear();
  FTrigSet.Clear();
  return;
}

// -----------------------------------------------------------

int XSecInputFileMgr_t::Load(const TString& inputFileName) {
  Clear();
  std::ifstream ifs;
  ifs.open(inputFileName.Data());
  if (!ifs.is_open()) {
    std::cout << "XSecInputFileMgr: failed to open file <" << inputFileName << ">\n";
    return 0;
  }
  FName=inputFileName;
  std::string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(state==0){
      stringstream ss1(line); ss1 >> FTotLumi;
      state++;
    }else if(state==1){
      FYieldsTag = TString(line);
      state++;
    }else if(state==2){
      FConstTag = TString(line);
      state++;
    } 
    else if (state==3) {
      FEvtEffScaleTag = TString(line);
      if (PosOk(line,"hltEff")) {
	std::cout << "input file probably does not contain a line for tagDir_EventEfficiencyScaleFactorConstants\n";
	assert(0);
      }
      state++;
    }
    else if (state==4) {
      FTrigSet = TString(line);
      break;
    }
  }
  ifs.close();
  if (!DYTools::checkTotalLumi(FTotLumi)) {
    std::cout << "file <" << inputFileName << "> has mismatching total lumi value\n";
    return 0;
  }

  std::cout << "Loaded:\n" << *this << "\n";
  return FTrigSet.Length();
}

// -----------------------------------------------------------
// -----------------------------------------------------------

int EScaleTagFileMgr_t::Load(const TString& inputFileName) {
  std::ifstream ifs;
  ifs.open(inputFileName.Data());
  if (!ifs.is_open()) {
    std::cout << "EScaleTagFileMgr: failed to open file <" << inputFileName << ">\n";
    return 0;
  }
  std::string line;
  TString tag;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') break;
    getline(ifs,line);
    stringstream ss3(line); ss3 >> tag;
    FEScaleTags.push_back(tag);
  }
  ifs.close();
  std::cout << "Loaded:\n" << *this << "\n";
  return FEScaleTags.size();
}

// -----------------------------------------------------------
*/

// -----------------------------------------------------------
// -----------------------------------------------------------

