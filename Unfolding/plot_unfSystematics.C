#include "../Efficiency/plot_Systematics.C"

void plot_unfSystematics(TString flags="111", int saveCanvas=0) {
  plot_Systematics("unfYields",flags,saveCanvas);
}
