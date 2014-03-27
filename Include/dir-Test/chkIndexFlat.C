#include "../DYTools.hh"

void chkIndexFlat() {

  for (int im=0; im<DYTools::nMassBins; ++im) {
    for (int iy=0; iy<DYTools::nYBins[im]; ++iy) {
      int idx=DYTools::findIndexFlat(im,iy);
      printf("im=%d, iy=%d, flat_idx=%d\n",im,iy,idx);
    }
  }

  for (int im=0; im<DYTools::nMassBins; ++im) {
    for (int iy=0; iy<DYTools::nYBinsMax; ++iy) {
      int idx=DYTools::findIndexFlat(im,iy);
      printf("im=%d, iy=%d, flat_idx=%d\n",im,iy,idx);
    }
  }
}
