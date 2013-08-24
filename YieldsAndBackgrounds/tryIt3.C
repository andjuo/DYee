#include "../Include/MyTools.hh"

void tryIt3() {
  TFile test("test.root","recreate");
  writeIntFlagValues("chk",3,1,2,3);
  test.Close();

  TFile test2("test.root","read");
  TVectorD* values=readFlagValues(test2,"chk",3);
  test2.Close();
  values->Print();
  delete values;
}
