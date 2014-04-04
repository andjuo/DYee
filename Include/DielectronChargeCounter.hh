#ifndef DielectronChargeCounter_HH
#define DielectronChargeCounter_HH


#include "../Include/DYTools.hh"
#include "../Include/CPlot.hh"
#include "../Include/TDielectron.hh"
#include "../Include/MyTools.hh"
#include <TH2D.h>

// -------------------------------------------

#ifndef MYTOOLS_HH
inline
TH2D* createBaseH2(const TString &histoName, const TString &histoTitle="", int absRapidity=1) {
  double yMin=(absRapidity) ? DYTools::yRangeMin : -DYTools::yRangeMax;
  int yDivCount=(absRapidity) ? DYTools::nYBinsMax : 2*DYTools::nYBinsMax;
  TH2D* h2=new TH2D(histoName,histoTitle,
		    DYTools::nMassBins,DYTools::massBinLimits,
		    yDivCount,yMin, DYTools::yRangeMax);
  h2->SetDirectory(0);
  h2->GetXaxis()->SetTitle("M_{ee} [GeV]");
  h2->GetYaxis()->SetTitle("y");
  return h2;
}

// -------------------------------------------
// -------------------------------------------

inline
TH1D* createProfileX(TH2D *h2, int iyBin, const TString &name, int setTitle=0, const char *title=NULL) {
  /*
  const TArrayD *arr=h2->GetXaxis()->GetXbins();
  arr->Print();
  double bins[h2->GetNbinsX()+1];
  for (int i=0; i<=h2->GetNbinsX(); ++i) {
    bins[i]=(*arr)[i];
  }
  */
  TH1D *h=new TH1D(name,"",h2->GetNbinsX(),h2->GetXaxis()->GetXbins()->GetArray());
  h->SetDirectory(0);
  h->SetStats(0);
  if (setTitle) {
    if (title) h->SetTitle(title);
    else h->SetTitle(name);
  }
  h->GetXaxis()->SetTitle( h2->GetXaxis()->GetTitle() );
  for (int ibin=1; ibin<=h2->GetNbinsX(); ibin++) {
    h->SetBinContent(ibin,h2->GetBinContent(ibin,iyBin));
    h->SetBinError(ibin,h2->GetBinError(ibin,iyBin));
  }
  return h;
}
#endif

// -------------------------------------------
// -------------------------------------------

class DielectronChargeCounter_t {
protected:
public:
  TH2D *Fh2pp,*Fh2pm,*Fh2mm,*Fh2ss;
  TString fNameBase;
  int FAbsRapidity;
public:
  DielectronChargeCounter_t(const TString &nameBase, int use_absRapidity=1) :
    Fh2pp(NULL), Fh2pm(NULL), Fh2mm(NULL), Fh2ss(NULL),
    fNameBase(nameBase),
    FAbsRapidity(use_absRapidity)
  {
    TString namePP= nameBase + TString("_pp");
    TString namePM= nameBase + TString("_pm");
    TString nameMM= nameBase + TString("_mm");
    TString nameSS= nameBase + TString("_ss");
    TH2D *h2=createBaseH2("dcc_tmp",use_absRapidity);
    Fh2pp = (TH2D*)h2->Clone(namePP);
    Fh2pp->SetTitle(namePP);
    Fh2pm = (TH2D*)h2->Clone(namePM);
    Fh2pm->SetTitle(namePM);
    Fh2mm = (TH2D*)h2->Clone(nameMM);
    Fh2mm->SetTitle(nameMM);
    Fh2ss = (TH2D*)h2->Clone(nameSS);
    Fh2ss->SetTitle(nameSS);
  }

  int ptrsOk() const { return ((Fh2pp!=NULL) && (Fh2pm!=NULL) && (Fh2mm!=NULL) && (Fh2ss!=NULL)) ? 1:0; }

  TH2D* getH2PP() { return Fh2pp; }
  TH2D* getH2PM() { return Fh2pm; }
  TH2D* getH2MM() { return Fh2mm; }
  TH2D* getH2SS() { return Fh2ss; }
  const TString &GetName() const { return fNameBase; }
  int absRapidity() const { return FAbsRapidity; }

  void SetDirectory(void*) const {}  // a dummy
  bool InheritsFrom(const TString &testClass) const {
    return (testClass==TString("TH1")) ? true : false;
  }


  double calcIntegral(int iYbin, int the_case) const {
    TH2D* h2=NULL;
    if (the_case==1) h2=Fh2pp;
    else if (the_case==2) h2=Fh2ss;
    else if (the_case==-1) h2=Fh2mm;
    else if (the_case==0) h2=Fh2pm;
    if (!h2) {
      std::cout << "calcIntegral: wrong the_case=" << the_case << "\n";
      return 0;
    }
    double sum=0.;
    for (int ibin=1; ibin<=h2->GetNbinsX(); ++ibin) {
      sum+=h2->GetBinContent(ibin,iYbin);
    }
    return sum;
  }


  TH2D* calcH2Tot(int setTitle=0, const char *title=NULL) const {
    TString name=Fh2pm->GetName();
    name.ReplaceAll("_pm","_tot");
    TH2D* h2=(TH2D*)Fh2pm->Clone(name);
    h2->SetDirectory(0);
    h2->Add(Fh2ss);
    if (setTitle) {
      if (title) h2->SetTitle(name);
      else h2->SetTitle(title);
    }
    return h2;
  }

  TH2D* calcH2SSFrac(TH2D* h2Tot, int setTitle=0, const char *title=NULL) const {
    if (!h2Tot) h2Tot=this->calcH2Tot();
    TString name=Fh2ss->GetName();
    name.ReplaceAll("_ss","_SSfrac");
    TH2D* h2=(TH2D*)Fh2ss->Clone(name);
    h2->SetDirectory(0);
    h2->Divide(h2Tot);
    if (setTitle) {
      if (title) h2->SetTitle(title);
      else h2->SetTitle(name);
    }
    return h2;
  }

  TH1D* calcH2TotProfile(TH2D* h2Tot, int iYbin, int setTitle=0, const char *title=NULL) const {
    if (!h2Tot) h2Tot=this->calcH2Tot(1);
    TH1D *h1tot=NULL;
    h1tot=createProfileX(h2Tot,iYbin,
		     h2Tot->GetTitle() + TString(Form("_profile%d",iYbin)),
		     setTitle,title);
    return h1tot;
  }

  TH1D* calcH2SSFracProfile(TH2D* h2Tot, TH2D* h2ssFrac, int iYbin, int setTitle=0, const char *title=NULL) const {
    if (!h2ssFrac) h2ssFrac=this->calcH2SSFrac(h2Tot,1);
    TH1D *h1ssFrac=NULL;
    h1ssFrac=createProfileX(h2ssFrac,iYbin,
			    h2ssFrac->GetTitle() + TString(Form("_profile%d",iYbin)),
			    setTitle,title);
    return h1ssFrac;
  }

  template<class dielectron_t>
  void Fill(const dielectron_t* d, double weight) {
    int sameSign=0;
    if (d->q_1 == d->q_2) {
      sameSign=d->q_1;
    }
    const double mass=d->mass;
    const double rap=(FAbsRapidity) ? fabs(d->y) : d->y;
    if (sameSign==0) Fh2pm->Fill(mass,rap, weight);
    else {
      if (sameSign>0) Fh2pp->Fill(mass,rap, weight);
      else Fh2mm->Fill(mass,rap, weight);
      Fh2ss->Fill(mass,rap, weight);
    }
  }

  int Write(TFile &file) const {
    file.cd();
    Fh2pp->Write(Fh2pp->GetName());
    Fh2pm->Write(Fh2pm->GetName());
    Fh2mm->Write(Fh2mm->GetName());
    Fh2ss->Write(Fh2ss->GetName());
    return 1;
  }

  int Read(TFile &file) {
    file.cd();
    Fh2pp->Read(Fh2pp->GetName());
    Fh2pm->Read(Fh2pm->GetName());
    Fh2mm->Read(Fh2mm->GetName());
    Fh2ss->Read(Fh2ss->GetName());
    Fh2pp->SetDirectory(0);
    Fh2pm->SetDirectory(0);
    Fh2mm->SetDirectory(0);
    Fh2ss->SetDirectory(0);
    return 1;
  }

  int Write(const TString &dummyStr) const {
    if (0) std::cout << dummyStr << "\n";
    Fh2pp->Write(Fh2pp->GetName());
    Fh2pm->Write(Fh2pm->GetName());
    Fh2mm->Write(Fh2mm->GetName());
    Fh2ss->Write(Fh2ss->GetName());
    return 1;
  }

  int Read(const TString &dummyStr) {
    if (0) std::cout << dummyStr << "\n";
    Fh2pp->Read(Fh2pp->GetName());
    Fh2pm->Read(Fh2pm->GetName());
    Fh2mm->Read(Fh2mm->GetName());
    Fh2ss->Read(Fh2ss->GetName());
    return 1;
  }

  void SetDirectory(TDirectory *dir) {
    Fh2pp->SetDirectory(dir);
    Fh2pm->SetDirectory(dir);
    Fh2mm->SetDirectory(dir);
    Fh2ss->SetDirectory(dir);
  }

};


// -------------------------------------------

#endif
