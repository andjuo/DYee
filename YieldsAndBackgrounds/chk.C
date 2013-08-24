#include "../Include/DYTools.hh"

void chk() {
  double mass=10.;
  mass=0; std::cout << "mass=" << mass << ", validMass=" << DYTools::validMass(mass) << "\n";
  mass=10; std::cout << "mass=" << mass << ", validMass=" << DYTools::validMass(mass) << "\n";
  mass=15; std::cout << "mass=" << mass << ", validMass=" << DYTools::validMass(mass) << "\n";
  mass=1000; std::cout << "mass=" << mass << ", validMass=" << DYTools::validMass(mass) << "\n";
  mass=1500-0.1; std::cout << "mass=" << mass << ", validMass=" << DYTools::validMass(mass) << "\n";
  mass=1500; std::cout << "mass=" << mass << ", validMass=" << DYTools::validMass(mass) << "\n";
  mass=1500.1; std::cout << "mass=" << mass << ", validMass=" << DYTools::validMass(mass) << "\n";
}
