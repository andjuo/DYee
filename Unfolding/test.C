#include "../Include/UnfoldingMatrix.h"

void test() {

  std::cout << "DYTools::nUnfoldingBins=" << DYTools::nUnfoldingBins << "\n";

  if (0) {
    UnfoldingMatrix_t m(UnfoldingMatrix::_cFSR, "test");
    FlatIndex_t fi;
    fi.setIdx(2350,1.41);
    std::cout << fi << "\n";
    m.fillIni(fi, 1.);
    m.getIniM()->Print();

    std::cout << "yRangeEdge=" << DYTools::yRangeEdge << ", yRangeMax=" << DYTools::yRangeMax << "\n";
  }
  else if (1) {

    FlatIndex_t fi;
    for (int im=0; im<DYTools::nMassBins; ++im) {
      double mass=DYTools::findMassBinCenter(im);
      for (int iy=0; iy<DYTools::nYBins[im]; ++iy) {
	double y=DYTools::findAbsYValue(im,iy);
	fi.setIdx(mass,y);
	std::cout << "mass=" << mass << ", y=" << y << ", fi=" << fi << "\n";
      }
      fi.setIdx(mass,2.4);
      std::cout << "chk mass=" << mass << ", y=" << 2.4 << ", fi=" << fi << "\n";
      fi.setIdx(mass,2.7);
      std::cout << "chk mass=" << mass << ", y=" << 2.7 << ", fi=" << fi << "\n";
    }
    std::cout << "\n";

    fi.setIdx(6000,2.6); 
    std::cout << "chk::  " << fi << "\n";
  }
  else if (0) {
    // test implemented unfolding

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
    std::cout << "res=" << res << "\n";
  }
}
