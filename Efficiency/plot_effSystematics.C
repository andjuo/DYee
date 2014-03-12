#include "../Efficiency/plot_Systematics.C"

void plot_effSystematics(TString flags="111", int saveCanvas=0) {
  plot_Systematics("efficiency",flags,saveCanvas);
}
