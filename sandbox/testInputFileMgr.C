#include "../Include/InputFileMgr.hh"

void testInputFileMgr(TString fname="../config_files/data_vilnius.conf") {
  InputFileMgr_t mgr;
  if (!mgr.Load(fname)) {
    std::cout << "failed to load\n";
    return;
  }
  std::cout << "loaded fname=<" << fname << "> ok\n";
  std::cout << mgr;
  return;
}
