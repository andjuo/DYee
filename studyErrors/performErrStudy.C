#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "studyHelpers.hh"


// -------------------------------------------------------------
// -------------------------------------------------------------

int performErrStudy() {

  // define info to use
  std::vector<TString> fileTagV;
  std::vector<TString> explainV;
  fileTagV.reserve(20); explainV.reserve(20);
  fileTagV.push_back("_-1_-1__"); explainV.push_back("default");
  //fileTagV.push_back("_203_201__"); explainV.push_back("#Gamma 3, #Gamma 1");

  std::cout << "there are " << fileTagV.size() << " file tags\n";

  // construct file names
  std::vector<TString> fnamesV;
  fnamesV.reserve(fileTagV.size());
  for (unsigned int i=0; i<fileTagV.size(); ++i) {
    TString fname="errStudy" + fileTagV[i] + TString("nExps1000.root");
    fnamesV.push_back(fname);
  }

  // calculate the errors
  std::vector<StudyInfo_t*> infoV;
  infoV.reserve(fnamesV.size());

  for (unsigned int i=0; i<fnamesV.size(); ++i) {
    std::cout << "loading " << fnamesV[i] << "\n";
    StudyInfo_t *info=new StudyInfo_t();
    if (!info->getSidedErrors(fnamesV[i],fileTagV[i])) {
      std::cout << "failed to get info\n";
      return retCodeError;
    }
    infoV.push_back(info);
    info->PrintValues(40, _valTBYieldGen);
    info->PrintValues(41, _valTBYieldGen);
    info->PrintValues(40, _valTBYieldTot);
    info->PrintValues(41, _valTBYieldTot);
  }

  return retCodeOk;
}
// -------------------------------------------------------------
// -------------------------------------------------------------
