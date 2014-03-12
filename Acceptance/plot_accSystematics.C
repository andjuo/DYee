#include "../Efficiency/plot_Systematics.C"

void plot_accSystematics(TString flags="111", int saveCanvas=0) {
  plot_Systematics("acceptance",flags,saveCanvas);
}
