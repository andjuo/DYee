#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"

//const double reweightCorrection=0.9673619; // count_Z(data)/count_Z(MC)

void plotShapeReweight(int analysisIs2D) {
  if (!DYTools::setup(analysisIs2D)) return;

  TString fname="../../Results-DYee/root_files_reg/yield/DY_j22_19712pb_ApplyEscale/bg-subtracted_yield_1D__peak20140220.root";
  if (analysisIs2D) fname.ReplaceAll("1D","2D");

  TH2D* h2=LoadHisto2D("zeeMCShapeReweight_ddBkg",fname,"ShapeReweight",1);
  if (!h2) return;

  //h2->Scale(1/reweightCorrection);
  removeError(h2);

  std::vector<TH2D*> hV;
  std::vector<TString> labelV;
  hV.push_back(h2);
  labelV.push_back("weight");

  if (!analysisIs2D) {
    if (0) {
      TH1D* h1=createProfileAuto(h2,1,"hShape",1,"escale residual weights");
      TCanvas *cx=new TCanvas("cx","cx",700,700);
      cx->SetLogx();
      h1->Draw();
      cx->Update();
      return;
    }
  }
  TCanvas *cx= plotProfiles("cx",hV,labelV,NULL,1,"escale residual weights");
  cx->Update();
}
