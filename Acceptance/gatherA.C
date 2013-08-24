#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"

void gatherA() {
  TH2D *base=createBaseH2("base","",1);
  
  TString foutName="dyee_Acc_1D_20130823.root";
  TFile fout(foutName,"recreate");

  for (int withFEWZ=0; withFEWZ<2; ++withFEWZ) {
    const TString str=(withFEWZ) ? "withFEWZ" : "noFEWZ";
    const TString fname=(withFEWZ) ? "acceptance_1D.root" : "acceptance_1D_noFEWZ.root";
    const TString fullFName=TString("../root_files/constants/DY_j22_19789pb/") + fname;
    TFile file(fullFName,"read");
    TH2D *hAcc=Clone(base,"hAcc",1);
    TH2D *hPass=Clone(base,"hSumPass",1);
    TH2D* hTotal=Clone(base,"hSumTotal",1);
    TH2D *hSumPassPreFsr=Clone(base,"hSumPassPreFsr",1);
    TH2D *hSumTotalPreFsr=Clone(base,"hSumTotalPreFsr",1);
    TH2D *hAccPreFsr=Clone(base,"hAccPreFsr",1);
    int res=
      loadHisto(file,&hAcc,"") &&
      loadHisto(file,&hPass,"") &&
      loadHisto(file,&hTotal,"") &&
      loadHisto(file,&hSumPassPreFsr,"") &&
      loadHisto(file,&hSumTotalPreFsr,"");
    if (!res) return;
    file.Close();

    hAccPreFsr->Divide(hSumPassPreFsr,hSumTotalPreFsr,1.,1.,"b");
    
    TString uStr=str;
    TH1D *h1Acc=createProfileX(hAcc,1,"hAcc_" + uStr,1,NULL);
    TH1D *h1Pass=createProfileX(hPass,1,"hPass_" + uStr,1,NULL);
    TH1D *h1Total=createProfileX(hTotal,1,"hTotal_" + uStr,1,NULL);
    TH1D *h1AccPreFsr=createProfileX(hAccPreFsr,1,"hAccPreFsr_" + uStr,1,NULL);
    TH1D *h1PassPreFsr=createProfileX(hSumPassPreFsr,1,"hPassPreFsr_" + uStr,1,NULL);
    TH1D *h1TotalPreFsr=createProfileX(hSumTotalPreFsr,1,"hTotalPreFsr_" + uStr,1,NULL);

    printHisto(h1Acc);

    fout.cd();
    //fout.mkdir(str);
    //fout.cd(str);
    h1Acc->Write();
    h1Pass->Write();
    h1Total->Write();
    h1AccPreFsr->Write();
    h1PassPreFsr->Write();
    h1TotalPreFsr->Write();
  }
  fout.Close();
  return;
}
