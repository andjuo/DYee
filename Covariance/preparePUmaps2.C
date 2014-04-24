// Macro assumes that there are available histograms
// at -15%, -10%, -9%, ..., 9%, 10%, 15% cross section weight

#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/PUReweight.hh"
#include "../Include/ComparisonPlot.hh"
#include <TRandom.h>

// --------------------------------------------------------

int printHisto(const TH1F* histo, int exponent=0, int maxLines=-1);
TH1D* convert_TH1F_to_TH1D(const TH1F *h, TString newName);
TH1F* convert_TH1D_to_TH1F(const TH1D *h, TString newName);
TH1D* loadHisto(TString fname, TString fieldName="pileup");

// --------------------------------------------------------

void preparePUmaps2(int nIters, int iSeed=-1) {
  gRandom->SetSeed(iSeed);

  // automatic handling
  PUReweight_t puBase(DYTools::NO_SYST);
  PUReweight_t pu5p(DYTools::PILEUP_5plus);
  PUReweight_t pu5m(DYTools::PILEUP_5minus);

  // manual handling
  TString fnameSrc="../root_files_reg/pileup/8TeV_reg/mcPileupHildreth_mean_full2012_20131106_repacked.root";
  TString fnameTarget="../root_files_reg/pileup/8TeV_reg/dataPileupHildreth_mean_full2012_20131106_repacked.root";
  TString fnameTarget5p=fnameTarget;
  TString fnameTarget5m=fnameTarget;
  fnameTarget5p.ReplaceAll(".root","_plus5percent.root");
  fnameTarget5m.ReplaceAll(".root","_minus5percent.root");

  TH1F *source_f=NULL;
  TH1F *target_f=NULL, *target5p_f=NULL, *target5m_f=NULL;
  TH1D *source=NULL;
  TH1D *target=NULL, *target5p=NULL, *target5m=NULL;

  {
    TFile file(fnameSrc,"read");
    source_f=(TH1F*)file.Get("pileup_simulevel_mc");
    source_f->SetDirectory(0);
    file.Close();
  }
  {
    TFile file(fnameTarget,"read");
    target_f=(TH1F*)file.Get("pileup_lumibased_data");
    target_f->SetDirectory(0);
    file.Close();
  }
  {
    TFile file(fnameTarget5p,"read");
    target5p_f=(TH1F*)file.Get("pileup_lumibased_data");
    target5p_f->SetDirectory(0);
    file.Close();
  }
  {
    TFile file(fnameTarget5m,"read");
    target5m_f=(TH1F*)file.Get("pileup_lumibased_data");
    target5m_f->SetDirectory(0);
    file.Close();
  }

  if (!source_f || !target_f || !target5p_f || !target5m_f) {
    std::cout << "failed to load the starting distributions\n";
    return;
  }

  source  = convert_TH1F_to_TH1D(source_f,"hsource");
  target  = convert_TH1F_to_TH1D(target_f,"hbase");
  target5p= convert_TH1F_to_TH1D(target5p_f,"h5plus");
  target5m= convert_TH1F_to_TH1D(target5m_f,"h5minus");

  if (!source || !target || !target5p || !target5m) {
    std::cout << "failed conversion to TH1D\n";
    return;
  }

  if (0) {
    printHisto(source_f);
    printHisto(target_f);
    printHisto(target5p_f);
    printHisto(target5m_f);
    //return;
  }
  if (0) {
    printHisto(source);
    printHisto(target);
    printHisto(target5p);
    printHisto(target5m);
    //return;
  }

  // ---------------------------------------------------------
  // ---------------------------------------------------------

  removeError1D(source);

  if (1) {
    target->Scale(1e-3);
    target5p->Scale(1e-3);
    target5m->Scale(1e-3);

    TCanvas *cx=new TCanvas("cxStart","cxStart",700,700);
    ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,"cp","",
			"nPU","count","ratio");
    cp.Prepare2Pads(cx);
    cp.AddHist1D(source,"MC","LP",TAttMarker(kRed+1,kDot,1.),1,0,1);
    cp.AddHist1D(target,"base","LP",TAttMarker(kBlack,24,0.8),1,0,1);
    cp.AddHist1D(target5p,"+5%","LP",TAttMarker(kBlue,kFullTriangleUp,0.8),1,0,1);
    cp.AddHist1D(target5m,"-5%","LP",TAttMarker(kGreen+1,kOpenTriangleUp,0.8),1,0,1);
    cp.Draw(cx);
    cx->Update();

    target->Scale(1e3);
    target5p->Scale(1e3);
    target5m->Scale(1e3);
  }

  // ---------------------------------------------------------

  TH1D* h1PU_15minus=loadHisto("MyDataPileupHistogram_true_58990.root");
  TH1D* h1PU_15plus =loadHisto("MyDataPileupHistogram_true_79810.root");
  if (!h1PU_15minus || !h1PU_15plus) return;

  std::vector<TH1D*> h1PUvec;
  h1PUvec.reserve(21);
  double xsecBase=69400;
  int idx5minus=-1, idx5plus=-1, idx0=-1;
  for (int fraction=-10; fraction<=10; fraction++) {
    double xsec=xsecBase*(1+fraction*0.01);
    std::cout << "fraction=" << fraction << ", xsec=" << xsec << "\n";
    TH1D* h1PU=loadHisto(Form("MyDataPileupHistogram_true_%1.0lf.root",xsec));
    if (!h1PU) return;
    if (fraction==-5) idx5minus=h1PUvec.size();
    else if (fraction==0) idx0=h1PUvec.size();
    else if (fraction==5) idx5plus=h1PUvec.size();
    h1PUvec.push_back(h1PU);
  }

  // ---------------------------------------------------------

  if (0) {
    TCanvas *cx=new TCanvas("cxStart","cxStart",700,700);
    ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,"cp","",
			"nPU","count","ratio");
    cp.Prepare2Pads(cx);
    cp.AddHist1D(target,"base","LP",TAttMarker(kBlack,24,0.8),1,0,1);
    cp.AddHist1D(target5p,"+5%","LP",TAttMarker(kBlue,kFullTriangleUp,0.8),1,0,1);
    cp.AddHist1D(target5m,"-5%","LP",TAttMarker(kGreen+1,kOpenTriangleUp,0.8),1,0,1);
    cp.AddHist1D(h1PUvec[idx0],"base","LP",TAttMarker(kBlack,24,0.8),1,0,1);
    cp.AddHist1D(h1PUvec[idx5plus ],"+5%","LP",TAttMarker(kBlue,kFullTriangleUp,0.8),1,0,1);
    cp.AddHist1D(h1PUvec[idx5minus],"-5%","LP",TAttMarker(kGreen+1,kOpenTriangleUp,0.8),1,0,1);
    cp.Draw(cx);
    cx->Update();
  }

  // ---------------------------------------------------------

  std::vector<TH1D*> hRndV;
  hRndV.reserve(nIters);

  for (int i=0; i<nIters; i++) {
    double r=0.05*gRandom->Gaus(0,1.);
    if (fabs(r)>0.15) {
      std::cout << "random r=" << r << "\n";
      return;
    }
    double fr = 1. + r;
    double xsecFracA=0, xsecFracB=0.;

    TH1D *h1a=NULL, *h1b=NULL;
    if ((fr>= 0.85) && (fr<=0.90)) {
      h1a=h1PU_15minus; xsecFracA=0.85;
      h1b=h1PUvec[0];   xsecFracB=0.90;
    }
    else if ((fr>=1.10) && (fr<=1.15)) {
      h1a=h1PUvec.back(); xsecFracA=1.10;
      h1b=h1PU_15plus;    xsecFracB=1.15;
    }
    else {
      for (int fraction=-10; fraction<=10; ++fraction) {
	if ((r>=0.01*fraction) && (r<=0.01*(fraction+1))) {
	  //std::cout << "located fraction=" << fraction << "\n";
	  h1a=h1PUvec[fraction+10]; xsecFracA=1+0.01*fraction;
	  h1b=h1PUvec[fraction+11]; xsecFracB=1+0.01*(fraction+1);
	}
      }
    }
    if (!h1a || !h1b) {
      std::cout << "null histogram\n";
      std::cout << "r=" << r << ", fr=" << fr << "\n";
      if (h1a) std::cout << "h1a->GetName()=" << h1a->GetName() << "\n";
      if (h1b) std::cout << "h1b->GetName()=" << h1b->GetName() << "\n";
      return;
    }

    if (1 && (i<20)) {
      std::cout << "r=" << r << ", xsecFrac range=" << xsecFracA
		<< " .. " << xsecFracB << "\n";
      std::cout << "  " << h1a->GetName() << "\n";
      std::cout << "  " << h1b->GetName() << "\n";
    }

    // weight of the h1b histogram
    double w= (fr - xsecFracA)/(xsecFracB - xsecFracA);

    TString hName=Form("hDRnd_lumibased_data_%d",i);
    TString hTitle=Form("Rnd %d",i);
    std::cout << " - " << hName << "\n";
    TH1D* hRnd=(TH1D*) target->Clone(hName);
    hRnd->SetTitle(hTitle);

    //std::cout << "hRnd->GetNbinsX=" << hRnd->GetNbinsX()
    //<< ", h1a has " << h1a->GetNbinsX() << " bins\n";

    for (int ibin=1; ibin<=hRnd->GetNbinsX(); ++ibin) {
      double va=h1a->GetBinContent(ibin);
      double vb=h1b->GetBinContent(ibin);
      double v= w*vb + (1-w)*va;
      if ((i<5) && (ibin<10)) {
	std::cout << " check: ibin=" << ibin << ", w=" << w
		  << ", va=" << va << ", vb=" << vb
		  << ", v=" << v << "\n";
      }
      hRnd->SetBinContent(ibin, v);
      hRnd->SetBinError(ibin,0.);
    }
    hRndV.push_back(hRnd);
  }

  // ---------------------------------------------------------

  if (1) {
    TCanvas *cx=new TCanvas("cx","cx",700,700);
    ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,"cp","",
			"nPU","count","ratio");
    cp.Prepare2Pads(cx);
    cp.AddHist1D(target,"base","LP",TAttMarker(kBlack,24,0.8),1,0,1);
    cp.AddHist1D(target5p,"+5%","LP",TAttMarker(kBlue,kFullTriangleUp,0.8),1,0,1);
    cp.AddHist1D(target5m,"-5%","LP",TAttMarker(kGreen+1,kOpenTriangleUp,0.8),1,0,1);
    const unsigned int colorCount=5;
    const int autoColor[colorCount] = { kBlue+1, kGreen+1, kOrange+1, kRed+1, kViolet };
    for (unsigned int i=0; i<hRndV.size(); ++i) {
      cp.AddHist1D(hRndV[i],hRndV[i]->GetTitle(),"LP",TAttMarker(autoColor[i%colorCount],kDot,1.),1,0,1);
    }


    TH1D* targetClone=(TH1D*)target->Clone(Form("%s_clone",target->GetName()));
    TH1D* target5pClone=(TH1D*)target5p->Clone(Form("%s_clone",target5p->GetName()));
    TH1D* target5mClone=(TH1D*)target5m->Clone(Form("%s_clone",target5m->GetName()));

    cp.AddHist1D(targetClone,"base","LP",TAttMarker(kWhite,20,0.8),1,0,-1);
    cp.SkipInRatioPlots(targetClone);
    cp.AddHist1D(target5pClone,"+5%","LP",TAttMarker(kWhite,kFullTriangleUp,0.8),1,0,-1);
    cp.SkipInRatioPlots(target5pClone);
    cp.AddHist1D(target5mClone,"-5%","LP",TAttMarker(kWhite,kFullTriangleDown,0.8),1,0,-1);
    cp.SkipInRatioPlots(target5mClone);
    cp.Draw(cx);
    cx->Update();
  }

  // ---------------------------------------------------------

  if (1) {
    std::vector<TH1F*> saveV;
    saveV.reserve(hRndV.size());
    for (unsigned int i=0; i<hRndV.size(); ++i) {
      TString newName=Form("hRnd_lumibased_data_%d",i+1);
      TH1F* h=convert_TH1D_to_TH1F(hRndV[i],newName);
      h->SetDirectory(0);
      saveV.push_back(h);
    }

    TFile fout("randomized_pileup_20140424.root","recreate");
    target_f->Write("pileup_lumibased_data_base");
    target5p_f->Write("pileup_lumibased_data_111");
    target5m_f->Write("pileup_lumibased_data_-111");
    if (!saveVec(fout,saveV,"")) {
      std::cout << " .. failed\n";
    }
    else {
      fout.Close();
      std::cout << "saved to file <" << fout.GetName() << ">\n";
    }
  }

  return;
}

// --------------------------------------------------------
// --------------------------------------------------------

template<class histo_t>
double* getXrange(const histo_t *h) {
  double *arr=new double[h->GetNbinsX()+1];
  for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin) {
    double x1=h->GetBinLowEdge(ibin);
    arr[ibin-1]=x1;
  }
  int i=h->GetNbinsX();
  double x=h->GetBinLowEdge(i);
  double w=h->GetBinWidth(i);
  arr[i]=x+w;
  return arr;
}

template<class histo_t>
double* getXrangeCut(const histo_t *h, int &size, double cut=49.5) {
  size=0;
  for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin) {
    if (h->GetBinLowEdge(ibin) < cut) size++; else break;
  }

  double *arr=new double[size+1];
  for (int ibin=1; ibin<=size; ++ibin) {
    double x1=h->GetBinLowEdge(ibin);
    arr[ibin-1]=x1;
  }
  int i=size;
  double x=h->GetBinLowEdge(i);
  double w=h->GetBinWidth(i);
  arr[i]=x+w;

  if (1) {
    std::cout << "range= ";
    for (int ii=0; ii<=size; ++ii) {
      std::cout << " " << arr[ii];
    }
    std::cout << "\n";
  }
  return arr;
}


TH1D* convert_TH1F_to_TH1D(const TH1F *h, TString newName) {
  int size=0;
  //double *arr=getXrangeCut(h,size,59.5);
  double *arr=getXrangeCut(h,size,100.);
  TH1D *hd=new TH1D(newName,newName,size,arr);
  hd->SetDirectory(0);
  delete arr;
  for (int ibin=1; ibin<=size; ++ibin) {
    hd->SetBinContent(ibin, h->GetBinContent(ibin));
    hd->SetBinError(ibin, h->GetBinError(ibin));
  }
  return hd;
}


TH1F* convert_TH1D_to_TH1F(const TH1D *h, TString newName) {
  double *arr=getXrange(h);
  TH1F *hf=new TH1F(newName,newName,h->GetNbinsX(),arr);
  hf->SetDirectory(0);
  delete arr;
  for (int ibin=1; ibin<=h->GetNbinsX(); ++ibin) {
    hf->SetBinContent(ibin, h->GetBinContent(ibin));
    hf->SetBinError(ibin, h->GetBinError(ibin));
  }
  return hf;
}

// --------------------------------------------------------
// --------------------------------------------------------

TH1D* loadHisto(TString fname, TString fieldName) {
  TFile fin(fname,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file <" << fname << ">\n";
    return NULL;
  }
  TH1D* h1=(TH1D*)fin.Get(fieldName);
  if (!h1) std::cout << "failed to get histo <" << fieldName
		     << "> from file <" << fname << ">\n";
  else h1->SetDirectory(0);
  fin.Close();
  return h1;
}

// --------------------------------------------------------
// --------------------------------------------------------

int printHisto(std::ostream& out, const TH1F* histo, int exponent=0, int maxLines=-1) {
  if (!histo) {
    out << "printHisto: histo is null\n";
    return 0;
  }
  char buf[100];
  const char *format= (exponent) ?
    " %5.2f-%5.2f    %e    %e\n" :
    " %5.2f-%5.2f    %f    %f\n";

  out << "values of " << histo->GetName() << "\n";
  int imax=histo->GetNbinsX();
  int truncated=0;
  if ((maxLines>0) && (imax>maxLines)) { imax=maxLines; truncated=1; }
  for(int i=1; i<=imax; i++) {
    double x=histo->GetBinLowEdge(i);
    double w=histo->GetBinWidth(i);
    sprintf(buf,format,
	    x,x+w,histo->GetBinContent(i),histo->GetBinError(i));
    out << buf;
  }
  if (truncated) out << "... output truncated to " << maxLines << " lines\n";
  return 1;
}

//------------------------------------------------------------------------------------------------------------------------

int printHisto(const TH1F* histo, int exponent, int maxLines) { return printHisto(std::cout, histo, exponent, maxLines); }

//------------------------------------------------------------------------------------------------------------------------
