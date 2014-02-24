#include <TROOT.h>
#include <iostream>

void addInQuad(double x1, double x2, double x3=0.) {
  double sum2=x1*x1+x2*x2+x3*x3;
  std::cout << "sum(" << x1 << "^2, " << x2 << "^2, " << x3 << "^2 ) = " << sum2 << " = " << sqrt(sum2) << "^2\n";
}
