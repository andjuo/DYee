#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/ComparisonPlot.hh"

// ------------------------------------

TH1D* loadPUdistr(TString path, TString fname, TString fieldName,
		  int loadTH1F=0);

// ------------------------------------
// ------------------------------------

void plotPUdistr() {
  TH1D* hDataPU= loadPUdistr("../root_files_reg/pileup/8TeV_reg/","dataPileupHildreth_mean_full2012_20131106_repacked.root","pileup_lumibased_data",1);
  hDataPU->Scale(0.83e-3);
  removeError(hDataPU);
  //TString path="dir-1stver/";
  TString path="./";
  TH1D* h_noPU_noFEWZ= loadPUdistr(path,"pu_noPU_noFEWZ.root","npv_noPU_noFewz");
  TH1D* h_noPU_wFEWZ = loadPUdistr(path,"pu_noPU_wFEWZ.root" ,"npv_noPU_wFewz");
  TH1D *h_wPU_noFEWZ = loadPUdistr(path,"pu_wPU_noFEWZ.root" ,"npv_wPU_noFewz");
  TH1D *h_wPU_wFEWZ  = loadPUdistr(path,"pu_wPU_wFEWZ.root"  ,"npv_wPU_wFewz");

  if (0) {
    removeError(h_noPU_noFEWZ);
    removeError(h_noPU_wFEWZ);
    removeError(h_wPU_noFEWZ);
    removeError(h_wPU_wFEWZ);
  }

  if (1) {
    TH1D* h=hDataPU;
    h->SetName("tmp");
    hDataPU=Clone(h_noPU_noFEWZ,"hDataPU","");
    hDataPU->Reset();
    int ibinMax=h->GetNbinsX();
    if (hDataPU->GetNbinsX()<ibinMax) ibinMax=hDataPU->GetNbinsX();
    for (int ibin=1; ibin<=ibinMax; ++ibin) {
      hDataPU->SetBinContent(ibin, h->GetBinContent(ibin));
      hDataPU->SetBinError  (ibin, h->GetBinError(ibin));
    }
  }

  double scale=1/hDataPU->Integral();
  hDataPU->Scale(scale);
  h_noPU_noFEWZ->Scale(scale);
  h_noPU_wFEWZ->Scale(scale);
  h_wPU_noFEWZ->Scale(scale);
  h_wPU_wFEWZ->Scale(scale);

  TString cpName="compPlot";
  TString cpTitle="PU distributions in Zee MC";
  ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,cpName,cpTitle,
		      "#PU_{mean}","a.u.","reweighted/data");
  cp.SetXRange(0.,50.);
  cp.AddHist1D(hDataPU,"data","LP",TAttMarker(kBlack,24,0.8),3,0,1);
  cp.AddHist1D(h_noPU_noFEWZ,"no corr.","LP",TAttMarker(kBlue,20,0.8),3,0,1);
  cp.SkipInRatioPlots(h_noPU_noFEWZ);
  //cp.AddHist1D(h_noPU_wFEWZ, "w/FEWZ","LP",TAttMarker(kOrange,25,0.8),0,0,1);
  cp.AddHist1D(h_wPU_noFEWZ, "w/PU","LP",TAttMarker(kRed+1,5,0.8),2,0,1);
  //cp.AddHist1D(h_wPU_wFEWZ, "w/PU w/FEWZ","LP",TAttMarker(kRed+1,2,0.8),0,0,1);

  TCanvas *cx=new TCanvas("cx","cx",800,800);
  cp.Prepare2Pads(cx);
  cp.Draw(cx);
  cx->Update();
  return;
}

// ------------------------------------
// ------------------------------------

TH1D* loadPUdistr(TString path, TString fname,
		  TString fieldName, int loadTH1F) {
  TString fullFname=path + fname;
  TFile fin(fullFname,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file <" << fin.GetName() << ">\n";
    return NULL;
  }
  TH1D* h=NULL;
  if (loadTH1F) {
    TH1F* hOrig=(TH1F*)fin.Get(fieldName);
    if (!hOrig) {
      std::cout << "failed to get <" << fieldName << "> from <"
		<< fin.GetName() << ">\n";
      return NULL;
    }
    hOrig->SetName("hTmp");
    h= convert_TH1F_to_TH1D(hOrig,fieldName);
    delete hOrig;
  }
  else {
    h=(TH1D*)fin.Get(fieldName);
  }
  h->SetDirectory(0);
  fin.Close();
  return h;
}

// ------------------------------------
