#include "../Efficiency/plot_Systematics.C"

void plot_unfSystematics(int analysisIs2D,
			 TString flags="111", int saveCanvas=0) {
  plot_Systematics(analysisIs2D,"unfYields",flags,saveCanvas);
}
