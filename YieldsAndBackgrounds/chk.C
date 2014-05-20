#include "../Include/DYTools.hh"

void chk(int analysisIs2D) {
  if (!DYTools::setup(analysisIs2D)) return;

  double mass=10.;
  mass=0; std::cout << "mass=" << mass << ", validMass=" << DYTools::validMass(mass) << "\n";
  mass=9.99; std::cout << "mass=" << mass << ", validMass=" << DYTools::validMass(mass) << "\n";
  mass=10; std::cout << "mass=" << mass << ", validMass=" << DYTools::validMass(mass) << "\n";
  mass=15; std::cout << "mass=" << mass << ", validMass=" << DYTools::validMass(mass) << "\n";
  mass=1000; std::cout << "mass=" << mass << ", validMass=" << DYTools::validMass(mass) << "\n";
  mass=1500-0.1; std::cout << "mass=" << mass << ", validMass=" << DYTools::validMass(mass) << "\n";
  mass=1500; std::cout << "mass=" << mass << ", validMass=" << DYTools::validMass(mass) << "\n";
  mass=1500.1; std::cout << "mass=" << mass << ", validMass=" << DYTools::validMass(mass) << "\n";
  mass=2000.1; std::cout << "mass=" << mass << ", validMass=" << DYTools::validMass(mass) << "\n";
  mass=3000.1-0.2; std::cout << "mass=" << mass << ", validMass=" << DYTools::validMass(mass) << "\n";
  mass=3000.1; std::cout << "mass=" << mass << ", validMass=" << DYTools::validMass(mass) << "\n";
}
