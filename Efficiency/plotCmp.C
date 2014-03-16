#include "../Include/MyTools.hh"
#include "../Include/ComparisonPlot.hh"

void plotCmp(int eff) {
  TFile fin1("out1_EE_Jan.root","read");
  TH1D *hEff1=(TH1D*)fin1.Get("eff_postFSRcorr");
  TH1D *hAcc1=(TH1D*)fin1.Get("acc_postFSRcorr");
  hEff1->SetDirectory(0);
  hAcc1->SetDirectory(0);
  fin1.Close();
  TH1D* hAS=(eff) ? hEff1 : hAcc1;

  TString path="../root_files_reg/constants/DY_j22_19712pb/";
  TString fname2= path+ TString((!eff) ? "acceptance_1D.root" : "efficiency_1D.root");
  TString fieldName=(eff) ? "hEfficiency" : "hAcceptance";
  TString correctionName=(eff) ? "efficiency" : "acceptance";

  TFile fin2(fname2,"read");
  TH2D *h2Our=(TH2D*)fin2.Get(fieldName);
  h2Our->SetName("h2Our");
  h2Our->SetDirectory(0);
  fin2.Close();

  TH1D *hOurRaw=createProfileX(h2Our,1,fieldName + TString("Raw"));
  TH1D *hOur=removeLastBin(hOurRaw,fieldName);

  std::cout << "Alexey: "; printHisto(hAS);
  std::cout << "OurRaw: "; printHisto(hOurRaw);
  std::cout << "Our   : "; printHisto(hOur);

  ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,"cp","","#it{M}_{ee}",correctionName,"ratio");
  cp.SetLogx();
  cp.AddHist1D(hAS,"Alexey","LP",kBlack,1,0);
  cp.AddHist1D(hOur,"Our","LP",kBlue,1,0);
  
  TCanvas *cx= new TCanvas("cx","cx",700,700);
  cp.Prepare2Pads(cx);
  cp.Draw(cx);
  cp.TransLegend(0,-0.5);
  cx->Update();
}
