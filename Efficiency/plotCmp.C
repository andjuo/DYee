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
  TString label1="regressed en. (20-500,500-800,800+)";
  TH1D *h2=NULL, *h3=NULL;
  TString label2,label3;

  if (0) {
    TString fname3="../root_files_reg/constants/DY_j22_19712pb_20inf/efficiency_1D.root";
    if (!eff) fname3.ReplaceAll("efficiency","acceptance");
    TFile fin3(fname3,"read");
    TH2D* h2tmp=(TH2D*)fin3.Get(fieldName);
    TH1D* h1tmp=createProfileX(h2tmp,1,"h1tmp");
    label2="regressed en. (20-inf)";
    h2=removeLastBin(h1tmp,fieldName+TString("regEn20inf"));
    delete h1tmp;
    delete h2tmp;
    fin3.Close();
  }

  if (0) {
    TString fname4="../root_files/constants/DY_j22_19789pb/efficiency_1D.root";
    if (!eff) fname4.ReplaceAll("efficiency","acceptance");
    TFile fin4(fname4,"read");
    TH2D* h2tmp=(TH2D*)fin4.Get(fieldName);
    TH1D* h1tmp=createProfileX(h2tmp,1,"h1tmp");
    h3=removeLastBin(h1tmp,fieldName+TString("summer2012"));
    label3="old n-tuples (20-500,500-800,800+)";
    delete h1tmp;
    delete h2tmp;
    fin4.Close();
  }

  if (0 && !eff) {
    TString fname3="../root_files_reg/constants/DY_j22_19712pb/acceptance_1D-PU.root";
    //if (!eff) fname4.ReplaceAll("efficiency","acceptance");
    TFile fin3(fname3,"read");
    TH2D* h2tmp=(TH2D*)fin3.Get(fieldName);
    printHisto(h2tmp);
    TH1D* h1tmp=createProfileX(h2tmp,1,"h1tmp");
    h2=removeLastBin(h1tmp,fieldName+TString("_wPU"));
    label3="regressed (20-500,500-800,800+); wPU";
    delete h1tmp;
    delete h2tmp;
    fin3.Close();
  }

  if (0 && !eff) {
    //TString fname4="../root_files_reg/constants/DY_j22_19712pb/acceptance_1D-PU.root";
    TString fname4="../root_files_reg/constants/DY_j22_19712pb_NoReweight/acceptance_1D___NoReweight.root";
    //if (!eff) fname4.ReplaceAll("efficiency","acceptance");
    TFile fin4(fname4,"read");
    TH2D* h2tmp=(TH2D*)fin4.Get(fieldName);
    printHisto(h2tmp);
    TH1D* h1tmp=createProfileX(h2tmp,1,"h1tmp");
    h3=removeLastBin(h1tmp,fieldName+TString("_noFEWZ"));
    label3="regressed (20-500,500-800,800+); noFEWZ";
    delete h1tmp;
    delete h2tmp;
    fin4.Close();
  }

  std::cout << "Alexey: "; printHisto(hAS);
  //std::cout << "OurRaw: "; printHisto(hOurRaw);
  std::cout << "Our   : "; printHisto(hOur);

  ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,"cp","","#it{M}_{ee}",correctionName,"ratio");
  cp.SetLogx();
  cp.AddHist1D(hAS,"Alexey","LP",kBlack,1,0);
  cp.AddHist1D(hOur,label1,"LP",kBlue,1,0);
  if (h2) cp.AddHist1D(h2,label2,"LP",kGreen+1,1,0);
  if (h3) cp.AddHist1D(h3,label3,"LP",kRed+1,1,0);
  
  TCanvas *cx= new TCanvas("cx","cx",700,850);
  cp.Prepare2Pads(cx);
  cp.Draw(cx);
  cp.TransLegend(0,-0.6);
  if (1 || h2 || h3) {
    cp.TransLegend(-0.15,0.);
    cp.WidenLegend(0.15,0.);
  }
  cx->Update();
}
