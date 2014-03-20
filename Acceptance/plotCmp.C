#include "../Include/MyTools.hh"
#include "../Include/ComparisonPlot.hh"

void plotCmp() {
  TFile fin1("Alexey_acc_preFSR.root","read");
  TH2D* h2Alexey=(TH2D*)fin1.Get("hAcc");
  h2Alexey->SetDirectory(0);
  //TMatrixD *Acc=(TMatrixD*)fin.Get("Acc");
  //TMatrixD *AccErr=(TMatrixD*)fin.Get("AccErr");
  fin1.Close();


  //TString path="../root_files_reg/constants/DY_j22_19712pb/";
  //TString fname2= path+ TString("acceptance_2D.root");
  TString fname2="../root_files_reg/constants/DY_j22_19712pb_NoReweight/acceptance_2D___NoReweight.root";
  TString correctionName= "acceptance";

  TFile fin2(fname2,"read");
  TH2D *h2Our=(TH2D*)fin2.Get("hAccPreFsr");
  h2Our->SetDirectory(0);
  fin2.Close();

  std::vector<TH2D*> hIni;
  hIni.push_back(h2Alexey);
  hIni.push_back(h2Our);
  
  std::vector<TString> labelsV;
  labelsV.push_back("Alexey");
  labelsV.push_back("DYee");
  

  std::cout << "Alexey: "; printHisto(h2Alexey);
  std::cout << "Our   : "; printHisto(h2Our);


  std::vector<std::vector<TH1D*>*> hProfV;
  if (!createRapidityProfileVec(hIni,hProfV,labelsV)) return;
  std::cout << "hProfV.size=" << hProfV.size() << "\n";

  for (int im=1; im<7; ++im) {
    TString mStr=Form("M_%2.0lf_%2.0lf",DYTools::massBinLimits[im],DYTools::massBinLimits[im+1]);
    TString cName=TString("canv_") + mStr;
    TCanvas *cx=new TCanvas(cName,cName,500,600);

    ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,"cp",mStr,"|y|",correctionName,"ratio");
    cp.Prepare2Pads(cx);

    //cp.SetLogx();
    TH1D* h=(*hProfV[im])[0];
    h->SetMarkerStyle(20);

    cp.AddHist1D((*hProfV[im])[0],"Alexey","LP",kBlack,1,0);
    cp.AddHist1D((*hProfV[im])[1],"DYee","LP",kBlue,1,0);
  
    //cp.Prepare2Pads(cx);
    cp.Draw(cx);
    cp.TransLegend(-0.35,-0.6);
    if (0) {
      cp.TransLegend(-0.25,0.);
      cp.WidenLegend(0.15,0.);
    }
    cx->Update();
    TString figName=TString("fig-") + mStr;
    SaveCanvas(cx,figName,"dir-accPreFsr");
  }
  return;
}
