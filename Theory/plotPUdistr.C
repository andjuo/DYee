#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/ComparisonPlot.hh"

// ------------------------------------

TH1D* loadPUdistr(TString fname, TString fieldName, int loadTH1F=0);

// ------------------------------------
// ------------------------------------

void plotPUdistr() {
  TH1D* hDataPU= loadPUdistr("../root_files_reg/pileup/8TeV_reg/dataPileupHildreth_mean_full2012_20131106_repacked.root","pileup_lumibased_data",1);
  hDataPU->Scale(0.83e-3);
  removeError(hDataPU);
  TH1D* h_noPU_noFEWZ= loadPUdistr("pu_noPU_noFEWZ.root","npv_noPU_noFewz");
  TH1D* h_noPU_wFEWZ = loadPUdistr("pu_noPU_wFEWZ.root" ,"npv_noPU_wFewz");
  TH1D *h_wPU_noFEWZ = loadPUdistr("pu_wPU_noFEWZ.root" ,"npv_wPU_noFewz");
  TH1D *h_wPU_wFEWZ  = loadPUdistr("pu_wPU_wFEWZ.root"  ,"npv_wPU_wFewz");

  double scale=1/hDataPU->Integral();
  hDataPU->Scale(scale);
  h_noPU_noFEWZ->Scale(scale);
  h_noPU_wFEWZ->Scale(scale);
  h_wPU_noFEWZ->Scale(scale);
  h_wPU_wFEWZ->Scale(scale);

  TString cpName="compPlot";
  TString cpTitle="PU distributions in Zee MC";
  ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,cpName,cpTitle,
		      "#PU_{mean}","count","ratio");
  cp.AddHist1D(hDataPU,"data","LP",TAttMarker(kBlack,24,0.8),1,1,1);
  cp.AddHist1D(h_noPU_noFEWZ,"no corr.","LP",TAttMarker(kOrange,25,0.8),1,1,1);
  cp.AddHist1D(h_noPU_wFEWZ, "w/FEWZ","LP",TAttMarker(kBlue,20,0.8),2,1,1);
  cp.AddHist1D(h_wPU_noFEWZ, "w/PU","LP",TAttMarker(kGreen+1,5,0.8),3,1,1);
  cp.AddHist1D(h_wPU_wFEWZ, "w/PU w/FEWZ","LP",TAttMarker(kRed+1,2,0.8),4,1,1);

  TCanvas *cx=new TCanvas("cx","cx",800,800);
  cp.Prepare2Pads(cx);
  cp.Draw(cx);
  cx->Update();
  return;
}

// ------------------------------------
// ------------------------------------

TH1D* loadPUdistr(TString fname, TString fieldName, int loadTH1F) {
  TFile fin(fname,"read");
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
