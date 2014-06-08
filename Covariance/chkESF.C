#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"

void chkESF(int analysisIs2D, int plotHistos=-1) {
  if (!DYTools::setup(analysisIs2D)) return;

  TString fname1="../../Results-DYee/root_files_reg/constants/DY_j22_19712pb_egamma_Unregressed_energy/covRhoFileSF_nMB41_asymHLT_Unregressed_energy-allSyst_100_v2.root";
  TString fname2="covRhoFileSF_nMB41_asymHLT_Unregressed_energy-allSyst_100.root";
  TString fname3="covRhoFileSF_nMB41_asymHLT_Unregressed_energy-allSyst_1000.root";

  if (analysisIs2D) {
    std::vector<TString*> tmpV;
    tmpV.push_back(&fname1);
    tmpV.push_back(&fname2);
    tmpV.push_back(&fname3);
    replaceAll(tmpV,"nMB41","nMB7");
    replaceAll(tmpV,"1D","2D");
  }

  TH2D *hRho1=LoadMatrixFields(fname1,1,
			      "scaleFactor","scaleFactorErr",1);
  if (!hRho1) return;
  TH2D *hRho2=LoadMatrixFields(fname2,1,
			      "scaleFactor","scaleFactorErr",1);
  if (!hRho2) return;
  TH2D *hRho3=LoadMatrixFields(fname3,1,
			      "scaleFactor","scaleFactorErr",1);
  if (!hRho3) return;

  std::vector<TH2D*> hV;
  std::vector<TString> labelV;

  int plot1=((plotHistos>0) && ((plotHistos&1)>0)) ? 1:0;
  int plot2=((plotHistos>0) && ((plotHistos&2)>0)) ? 1:0;
  int plot3=((plotHistos>0) && ((plotHistos&4)>0)) ? 1:0;
  if (plotHistos<0) { plot1=1; plot2=1; plot3=1; }

  if (plot1) { hV.push_back(hRho1); labelV.push_back("current"); }
  if (plot2) { hV.push_back(hRho2); labelV.push_back("update 100"); }
  if (plot3) { hV.push_back(hRho3); labelV.push_back("update 1000"); }

  TCanvas *cx=plotProfiles("cx",hV,labelV,NULL,0,
			   "event efficiency scale factor #rho");
  cx->Update();
}
