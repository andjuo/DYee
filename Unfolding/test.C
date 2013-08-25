#include "../Include/UnfoldingMatrix.h"

void test() {
  UnfoldingMatrix_t m(UnfoldingMatrix::_cFSR, "test");


  TMatrixD U(3,3);
  for (int i=0, k=0; i<3; ++i) {
    for (int j=0; j<3; ++j, k++) {
      U(i,j)=k;
    }
  }
  TVectorD ini(3), fin(3);
  for (int i=0; i<3; ++i) ini[i]=3-i;

  fin= U*ini;
  std::cout << "ini="; ini.Print();
  std::cout << "U=" ; U.Print();
  std::cout << "fin="; fin.Print();

  TMatrixD tmp=U;
  U.Transpose(tmp);
  int res=unfold(fin,U,ini);
  std::cout << "\"unfolded\" fin="; fin.Print();
}
