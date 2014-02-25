#include "../Include/InputFileMgr.hh"
#include "../Include/ComparisonPlot.hh"
#include <sstream>
#include <fstream>
#include <TGraphAsymmErrors.h>
#include "convert2root_Lovedeep.C"



// ---------------------------------------------

int LoadTableReco(const std::string &fname, TString label, int isData, std::vector<TGraphAsymmErrors*> &grV, std::vector<TH1D*> &hSymmV, ThreeM_t &M, TDescriptiveInfo_t *info) {
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
  for (int i=0; i<2; ++i) {
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
  char xx;
  while (!fin.eof() && std::getline(fin,line)) {
    if (info) info->append(line);
    if (debug) std::cout << "got <" << line << ">\n";
    while (PosOk(line,'±')) { size_t p=line.find('±'); line[p]=' '; }
    if (debug) std::cout << "mdf <" << line << ">\n";
    std::stringstream ss(line);
    if ((lineNo+1)%6==0) {
      ss >> xx >> ptLow;
      ptHi=500.;
    }
    else ss >> ptLow >> ptHi;
    if (lineNo%6==0) ss >> etaLow >> xx >> etaHi;
    ss >> eff >> xx >> errLo;
    std::cout << "pt range " << ptLow << " - " << ptHi << "; eta range " << etaLow << " - " << etaHi << ", eff= " << eff << " ± " << errLo << "\n";
    if (isData!=1) {
      ss >> eff >> xx >> errLo;
      std::cout << "mc eff= " << eff << " ± " << errLo << "\n";
      if (isData==2) {
	// scale factor
	ss >> eff >> xx >> errLo >> xx >> errHi;
	std::cout << "scale factor " << eff << " ± " << errLo << " ± " << errHi << "\n";
	errLo=sqrt(errLo*errLo+errHi*errHi);
      }
    }
    errHi=errLo;
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

void convert2root_Ilya () {
  ThreeM_t MMC, MData, MSF;
  
  std::vector<TGraphAsymmErrors*> grV_mc,grV_data,grV_sf;
  std::vector<TH1D*> hV_mc, hV_data,hV_sf;
  if (!LoadTableReco("dir-Ilya/effRECO20130822.txt","effMC",0,grV_mc,hV_mc,MMC,NULL)) {
    std::cout << "failed to load MC\n";
    return;
  }
  if (!LoadTableReco("dir-Ilya/effRECO20130822.txt","effData",1,grV_data,hV_data,MData,NULL)) {
    std::cout << "failed to load data\n";
    return;
  }
  if (!LoadTableReco("dir-Ilya/effRECO20130822.txt","sfID",2,grV_sf,hV_sf,MSF,NULL)) {
    std::cout << "failed to load scale factor\n";
    return;
  }

  // From Ilya's presentation at Joint ECAL DPG/EGM POG meeting 
  // on Aug 22, 2013; https://indico.cern.ch/event/268599/
  const double RECOsyst_gap_inPerc   []= { 8.80, 7.30, 2.00, 0.93, 0.45, 0.65 };
  const double RECOsyst_nonGap_inPerc[]= { 5.50, 3.10, 1.20, 0.79, 0.41, 0.45 };
  
  TMatrixD RECOsyst(6,5);
  for (int i=0; i<6; ++i) {
    for (int j=0; j<5; ++j) {
      RECOsyst(i,j) = 0.01*( (j==2) ? RECOsyst_gap_inPerc[i] : RECOsyst_nonGap_inPerc[i]);
    }
  }
  std::cout << "RECO systematics: ";  RECOsyst.Print();

  // values from 2014 Feb 25
  TMatrixD RECOsyst_etaMax25(6,5);
  RECOsyst_etaMax25.Zero();
  RECOsyst_etaMax25(0,3)=0.047;
  RECOsyst_etaMax25(1,3)=0.045;
  RECOsyst_etaMax25(0,4)=0.047;
  RECOsyst_etaMax25(1,4)=0.045;
  RECOsyst_etaMax25(2,4)=0.043;
  RECOsyst_etaMax25(3,4)=0.04;
  RECOsyst_etaMax25(4,4)=0.035;
  RECOsyst_etaMax25(5,4)=0.03;

  std::cout << "RECO systematics due to |eta|<2.5: ";  RECOsyst_etaMax25.Print();

  TMatrixD RECOsyst_tot(6,5);
  for (int iEt=0; iEt<6; ++iEt) {
    for (int iEta=0; iEta<5; ++iEta) {
      RECOsyst_tot(iEt,iEta)= sqrt(pow(RECOsyst(iEt,iEta),2)+pow(RECOsyst_etaMax25(iEt,iEta),2));
    }
  }


  // In Lovedeep&Ilya presentation on Oct 28, 2013
  // the numbers below were named as "absolute errors on efficiencies"
  //const double IDsyst_gap_inPerc   []= { 4.30, 4.60, 2.70, 1.50, 0.28, 0.51 };
  //const double IDsyst_nonGap_inPerc[]= { 4.30, 4.20, 1.40, 0.43, 0.28, 0.45 };

  for (unsigned int i=0; i<grV_mc.size(); ++i) {
    grV_mc[i]->Print("range");
  }
  for (unsigned int i=0; i<grV_data.size(); ++i) {
    grV_data[i]->Print("range");
  }
  for (unsigned int i=0; i<grV_sf.size(); ++i) {
    grV_sf[i]->Print("range");
  }

  TCanvas *cx = new TCanvas("cx","cx",1200,900);

  std::vector<ComparisonPlot_t*> cpV;
  for (unsigned int i=0; i<5; ++i) {
    TString etaRange=grV_mc[i]->GetTitle();
    etaRange.ReplaceAll("effMC ","");

    ComparisonPlot_t *cp=new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,Form("cp_%d",i),etaRange,"p_{T}","RECO efficiency","RECO s.f.");
    cp->SetPrintRatio(1);
    cp->SetLogx();
    cp->SetYAxisTextSizes(0.08, 1., 0.07);
    cp->SetXAxisTextSizes(0.08, 1., 0.07);
    cp->SetYRange(0.65,1.0);
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


  if (1) {

    TFile fout("egammaRECO.root","recreate");
    MMC.Write(fout,"eff_mc");
    MData.Write(fout,"eff_data");
    MSF.Write(fout,"sf");
    RECOsyst.Write("sf_syst_rel_error_egamma");
    RECOsyst_etaMax25.Write("sf_syst_rel_error_maxEta25");
    RECOsyst_tot.Write("sf_syst_rel_error");
    saveVec(fout,grV_mc,"effRECO_MC");
    saveVec(fout,grV_data,"effRECO_Data");
    saveVec(fout,grV_sf,"sfRECO");
    saveVec(fout,hV_mc,"effRECO_MC_histo");
    saveVec(fout,hV_data,"effRECO_Data_histo");
    fout.Close();
    std::cout << "file <" << fout.GetName() << "> created\n";
    
    if (1) {
      TString fname="efficiency_TnP_1D_Full2012_dataRECO_fit-fitEtBins6EtaBins5egamma_PU.root";
      TDescriptiveInfo_t info;
      info.reserve(5);
      info.append("RECO efficiency scale factors");
      info.append("https://twiki.cern.ch/twiki/bin/viewauth/CMS/EgCommissioningAndPhysicsDeliverables");
      
      TFile foutD(fname,"recreate");
      MData.Write_for_main_code(foutD,0);
      RECOsyst.Write("sf_syst_rel_error_egamma");
      RECOsyst_etaMax25.Write("sf_syst_rel_error_maxEta25");
      RECOsyst_tot.Write("sf_syst_rel_error");
      info.Write("info");
      foutD.Close();
      std::cout << "file <" << foutD.GetName() << "> created\n";
      
      fname.ReplaceAll("dataRECO_fit-fit","mcRECO_count-count");
      TFile foutMC(fname,"recreate");
      MMC.Write_for_main_code(foutMC,1);
      RECOsyst.Write("sf_syst_rel_error_egamma");
      RECOsyst_etaMax25.Write("sf_syst_rel_error_maxEta25");
      RECOsyst_tot.Write("sf_syst_rel_error");
      info.Write("info");
      foutMC.Close();
      std::cout << "file <" << foutMC.GetName() << "> created\n";
    }
  }
}
