#include "../Efficiency/plot_Systematics.C"

void plot_accSystematics(int analysisIs2D,
			 TString flags="101", int saveCanvas=0) {
  plot_Systematics(analysisIs2D,"acceptance",flags,saveCanvas);
}
