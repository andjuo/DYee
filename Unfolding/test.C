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
  else if (0) {

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
  else if (1) {
    TMatrixD m1(DYTools::nMassBins,DYTools::nYBinsMax);
    m1.Zero();
    TMatrixD m2(m1);
    TVectorD v(DYTools::nUnfoldingBins), v2(DYTools::nUnfoldingBins);
    
    for (int i=0, fi=0; i<DYTools::nMassBins; ++i) {
      for (int j=0; j<DYTools::nYBins[i]; ++j, ++fi) {
	double val=(i+1) + 0.01*(j+1);
	m1(i,j) = val;
	v2(fi) = val;
	
	int flatIdx=DYTools::findIndexFlat(i,j);
	if (flatIdx!=fi) std::cout << Form("(%d,%d) fi=%d, flatIndex=%d",i,j,fi,flatIdx) << "\n";
      }
    }
    HERE("containers filled");
    flattenMatrix(m1,v);
    deflattenMatrix(v,m2);

    TMatrixD m3(m2);
    m3.Zero();
    deflattenMatrix(v2,m3);

    std::cout << "\n\n m1\n";
    m1.Print();
    std::cout << "\n\n v\n"; v.Print();

    HERE("chk differences\n");
    testMaxDiff("functions: ",m1,m2,0);
    testMaxDiff("auto flat idx; ",m1,m3,0);
    
  }
}
