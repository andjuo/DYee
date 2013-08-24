#ifndef JSONPARSER
#define JSONPARSER

#include <TString.h>
#include <vector>

class JsonParser {

public:
  
  JsonParser();
  virtual ~JsonParser(){};
  
  int Initialize(const TString &filename);
  bool HasRunLumi(int run, int lumi) const;
  void Reset();
  void Print() const;
  int AssembleNumber(const std::vector<int> &data) const;

  enum Depth {
    d_outside           = 0,
    d_insideOuterBracket = 1,
    d_insideInnerBracketNumberOne = 2,
    d_insideInnerBracketNumberTwo = 3
  };


private:

  std::vector <int>            _runList;
  std::vector < std::vector<int> >  _lumiStatusList;

  bool                    _isInitialized;

  ClassDef(JsonParser,1)

};
#endif
