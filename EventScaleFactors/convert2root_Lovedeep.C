#include "../Include/InputFileMgr.hh"
#include "../Include/ComparisonPlot.hh"
#include <sstream>
#include <fstream>
#include <TGraphAsymmErrors.h>


// ---------------------------------------------


struct ThreeM_t {
  TMatrixD M,MerrLo,MerrHi;
public:
  ThreeM_t() : M(6,5), MerrLo(6,5), MerrHi(6,5) {
    M.Zero(), MerrLo.Zero(); MerrHi.Zero();
  }

  void Zero() { M.Zero(); MerrLo.Zero(); MerrHi.Zero(); }


  void Write(TFile &fout, TString field) const {
    fout.cd();
    M.Write(field);
    MerrLo.Write(field + TString("_errLo"));
    MerrHi.Write(field + TString("_errHi"));
  }

  void Write_for_main_code(TFile &fout, int mc) const {
    fout.cd();
    M.Write((mc==1) ? "effArray2DWeighted" : "effArray2D");
    MerrLo.Write((mc==1) ? "effArrayErrLow2DWeighted" : "effArrayErrLow2D");
    MerrHi.Write((mc==1) ? "effArrayErrHigh2DWeighted" : "effArrayErrHigh2D");
  }
};

// ---------------------------------------------

int LoadTable(const std::string &fname, TString label, std::vector<TGraphAsymmErrors*> &grV, std::vector<TH1D*> &hSymmV, ThreeM_t &M, TDescriptiveInfo_t *info) {
  const int debug=1;

  M.Zero();

  std::ifstream fin(fname.c_str());
  if (!fin.is_open()) {
    std::cout << "failed to open a file <" << fname << ">\n";
    return 0;
  }
  if (debug) std::cout << "opened file <" << fname << ">\n";
  
  if (info) {
    info->reserve(info->size()+42);
    info->append(std::string(""));
    info->append(std::string("file: " ) + fname);
    info->append(std::string(""));
  }

  std::string line;
  for (int i=0; i<7; ++i) {
    std::getline(fin,line);
    if (debug) std::cout << line << "\n";
    if (info) info->append(line);
  }

  const int count=6;
  const double ptArr[count+1]={ 10., 15., 20., 30., 40., 50., 500. };
  double *xloArr=new double[count];
  double *xhiArr=new double[count];
  double *xcArr=new double[count];
  double *ycArr =new double[count];
  double *yloArr=new double[count];
  double *yhiArr=new double[count];
  
  for (int i=0; i<count; ++i) {
    xcArr[i]=0.5*(ptArr[i]+ptArr[i+1]);
    if (i==count-1) xcArr[count-1]=100.;
    xloArr[i]=xcArr[i]-ptArr[i];
    xhiArr[i]=ptArr[i+1]-xcArr[i];
  }

  double ptLow,ptHi,etaLow,etaHi,eff,errLo,errHi;
  int lineNo=0;
  int etaRange=0;
  while (!fin.eof() && std::getline(fin,line)) {
    if (info) info->append(line);
    if (debug) std::cout << "got <" << line << ">\n";
    std::stringstream ss(line);
    ss >> ptLow >> ptHi >> etaLow >> etaHi >> eff >> errLo >> errHi;
    ycArr[lineNo]=eff;
    yloArr[lineNo]=errLo;
    yhiArr[lineNo]=errHi;

    M.M(lineNo,etaRange)=eff;
    M.MerrLo(lineNo,etaRange)=errLo;
    M.MerrHi(lineNo,etaRange)=errHi;

    lineNo++;
    if (lineNo%count==0) {
      TString grName=label + TString(Form("_eta_%4.2lf_%4.2lf",etaLow,etaHi));
      TString title=label + TString(Form(" %4.2lf < |#eta| < %4.2lf",etaLow,etaHi));
      TGraphAsymmErrors *gr=new TGraphAsymmErrors(count,xcArr,ycArr,xloArr,xhiArr,yloArr,yhiArr);
      gr->SetName(grName);
      gr->SetTitle(title);
      //gr->SetDirectory(0);
      //gr->SetStats(0);
      grV.push_back(gr);
      lineNo=0;
      etaRange++;

      TString hName=TString("h_") + grName;
      TH1D *h=new TH1D(hName,hName,count,ptArr);
      h->SetDirectory(0);
      h->SetStats(0);
      h->GetXaxis()->SetTitle("{#itp}_{T}");
      h->GetYaxis()->SetTitle(label);
      for (int i=1; i<count; ++i) {
	h->SetBinContent(i, ycArr[i-1]);
	h->SetBinError  (i, 0.5*(yloArr[i-1]+yhiArr[i-1]));
      }
      hSymmV.push_back(h);
    }
  }
  fin.close();
  return 1;
}


// ------------------------------------

void convert2root_Lovedeep () {
  ThreeM_t MMC, MData, MSF;
  
  std::vector<TGraphAsymmErrors*> grV_mc,grV_data,grV_sf;
  std::vector<TH1D*> hV_mc, hV_data,hV_sf;
  if (!LoadTable("dir-Lovedeep/effiGsfIdMediumMC.txt","effMC",grV_mc,hV_mc,MMC,NULL)) {
    std::cout << "failed to load MC\n";
    return;
  }
  if (!LoadTable("dir-Lovedeep/effiGsfIdMediumData.txt","effData",grV_data,hV_data,MData,NULL)) {
    std::cout << "failed to load data\n";
    return;
  }
  if (!LoadTable("dir-Lovedeep/sFGsfIdMedium.txt","sfID",grV_sf,hV_sf,MSF,NULL)) {
    std::cout << "failed to load scale factor\n";
    return;
  }


  for (unsigned int i=0; i<grV_mc.size(); ++i) {
    grV_mc[i]->Print("range");
  }
  for (unsigned int i=0; i<grV_data.size(); ++i) {
    grV_data[i]->Print("range");
  }
  for (unsigned int i=0; i<grV_sf.size(); ++i) {
    grV_sf[i]->Print("range");
  }

  TCanvas *cx = new TCanvas("cx","cx",1200,800);

  std::vector<ComparisonPlot_t*> cpV;
  for (unsigned int i=0; i<5; ++i) {
    TString etaRange=grV_mc[i]->GetTitle();
    etaRange.ReplaceAll("effMC ","");

    ComparisonPlot_t *cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,Form("cp_%d",i),etaRange,"p_{T}","ID efficiency","ID s.f.");
    cp->SetLogx();
    cp->SetYRange(0.2,1.0);
    cp->SetRatioYRange(0.8,1.2);
    cp->AddHist1D(hV_mc[i],"h eff MC","LPE",kOrange,1,0);
    cp->AddHist1D(hV_data[i],"h eff data","LPE",kRed,2,0);
    cp->AddGraph(grV_mc[i],"eff MC","LPE1",kBlack,1,0);
    cp->AddGraph(grV_data[i],"eff data","LPE1",kBlue ,2,0);
    cpV.push_back(cp);
  }

  cpV[0]->Prepare6Pads(cx,1,"comp",0.49);

  for (int i=0; i<5; ++i) {
    cpV[i]->Draw6(cx,1,i+1);

    cx->cd(cpV[i]->getSubPadIdx6(1,i+1,1));
    grV_sf[i]->SetLineColor(kGreen+1);
    grV_sf[i]->SetMarkerColor(kGreen+1);
    grV_sf[i]->Draw("LPE same");
  }
  cx->Update();


  TFile fout("mediumID.root","recreate");
  MMC.Write(fout,"eff_mc");
  MData.Write(fout,"eff_data");
  MSF.Write(fout,"sf");
  saveVec(fout,grV_mc,"effMediumID_MC");
  saveVec(fout,grV_data,"effMediumID_Data");
  saveVec(fout,grV_sf,"sfMediumID");
  saveVec(fout,hV_mc,"effMediumID_MC_histo");
  saveVec(fout,hV_data,"effMediumID_Data_histo");
  fout.Close();
  std::cout << "file <" << fout.GetName() << "> created\n";

  if (1) {
    TString fname="efficiency_TnP_1D_Full2012_dataID_fit-fitEtBins6EtaBins5egamma_PU.root";
    TDescriptiveInfo_t info;
    info.reserve(5);
    info.append("mediumID official scale factors");
    info.append("Web page: https://twiki.cern.ch/twiki/bin/view/Main/EGammaScaleFactors2012");

    TFile foutD(fname,"recreate");
    MData.Write_for_main_code(foutD,0);
    info.Write("info");
    foutD.Close();
    std::cout << "file <" << foutD.GetName() << "> created\n";

    fname.ReplaceAll("dataID_fit-fit","mcID_count-count");
    TFile foutMC(fname,"recreate");
    MMC.Write_for_main_code(foutMC,1);
    info.Write("info");
    foutMC.Close();
    std::cout << "file <" << foutMC.GetName() << "> created\n";
  }
}
