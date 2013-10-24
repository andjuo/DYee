#include "../DYTools.hh"

int isInt(const TString &str) {
  int idx=str.Index(".");
  std::cout << "search for '.' in <" << str << ">, idx=" << idx << "\n";
  return (idx==-1) ? 1 : 0;
}

int valueAsInt(const TString &str) {
  return atoi(str.Data());
}

double valueAsDouble(const TString &str) {
  return atof(str.Data());
}

void setValue(const TString &str, double &v, int &idx, int isMass, int massIdx) {
  std::cout << "setValue(<" << str << ">)\n";
  if (isInt(str)) {
    std::cout << " .. isInt\n";
    idx=valueAsInt(str);
    v=(isMass) ? DYTools::findMassBinCenter(idx) : DYTools::findAbsYValue(massIdx,idx);
  }
  else {
    std::cout << " .. isNotInt\n";
    v=valueAsDouble(str);
    idx=(isMass) ? DYTools::findMassBin(v) : DYTools::findAbsYBin(massIdx,v);
  }
}


void chkYBins(const TString &massStr, const TString &yStr) {
  int im=-1, iy=-1;
  double m=-99, y=-99;

  setValue(massStr,m,im,1,-1);
  setValue(yStr,y,iy,0, im);

  std::cout << "got massStr=<" << massStr << ">, yStr=<" << yStr << ">\n";
  std::cout << "decoded:\n";
  std::cout << "\tim=" << im << ", mass=" << m << "\n";
  std::cout << "\tiy=" << iy << ", y=" << y << "\n";

  std::cout << "\nlooking for a flat index:\n";
  int idx=DYTools::findIndexFlat(im,iy);
  std::cout << " fi=" << idx << "\n";

  return;
}
