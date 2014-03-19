#ifndef calcCorrectionSyst_H
#define calcCorrectionSyst_H

#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/EventSelector.hh"


// -------------------------------------------------

struct ApplySystFlags_t {
  int FNoSyst, FPU, FFSR;
  public:

  ApplySystFlags_t() :
    FNoSyst(1), FPU(0), FFSR(0)
  {}

  ApplySystFlags_t(int set_noSyst, int set_pu, int set_fsr) : 
    FNoSyst(set_noSyst), FPU(set_pu), FFSR(set_fsr)
  {}

  ApplySystFlags_t(const ApplySystFlags_t &f) : 
    FNoSyst(f.FNoSyst), FPU(f.FPU), FFSR(f.FFSR)
  {}

  int totSyst() const { return FNoSyst; }
  int noSyst() const { return FNoSyst; }
  int pu() const { return FPU; }
  int fsr() const { return FFSR; }

  TString asTString() const {
    TString flags=Form("%d%d%d",FNoSyst,FPU,FFSR);
    return flags;
  }

  TString asDetailedTString() const {
    TString flags=Form("%d%d%d%d%d",FNoSyst,FPU,FPU,FFSR,FFSR);
    return flags;
  }

  void setNoSyst(int yes) { FNoSyst=yes; }
  void setPU(int yes) { FPU=yes; }
  void setFSR(int yes) { FFSR=yes; }

  void setForAcceptance() {
    FNoSyst=1; FPU=0; FFSR=1;
  }

  void setForEfficiency() {
    FNoSyst=1; FPU=1; FFSR=1;
  }
  
  void adjustForAcceptance() {
    FPU=0;
  }

  void adjustForEfficiency() {
  }

  void adjustForUnfYields() {
  }

  int assignFlags(const TString &str) {
    if (str.Length()!=3) {
      std::cout << "assignFlags(\"" << str << "\") size should be 3\n";
      return 0;
    }
    FNoSyst= (str[0]=='1') ? 1:0;
    FPU    = (str[1]=='1') ? 1:0;
    FFSR   = (str[2]=='1') ? 1:0;
    return 1;
  }
  
  int assignFlagsOdd(const TString &str) {
    if (str.Length()!=5) {
      std::cout << "assignFlags(\"" << str << "\") size should be 5\n";
      return 0;
    }
    FNoSyst= (str[0]=='1') ? 1:0;
    FPU    = (str[1]=='1') ? 1:0;
    FFSR   = (str[3]=='1') ? 1:0;
    return 1;
  }
  
  int assignFlagsEven(const TString &str) {
    if (str.Length()!=5) {
      std::cout << "assignFlags(\"" << str << "\") size should be 5\n";
      return 0;
    }
    FNoSyst= (str[0]=='1') ? 1:0;
    FPU    = (str[2]=='1') ? 1:0;
    FFSR   = (str[4]=='1') ? 1:0;
    return 1;
  }


  void add(const ApplySystFlags_t &f) {
    if (f.FNoSyst) FNoSyst=1;
    if (f.FPU) FPU=1;
    if (f.FFSR) FFSR=1;
  }
  
  static void printFlag(std::ostream &out, const char *name, int value) {
    out << name << ":";
    if (value==1) out << "on"; else out << "off";
  }


  friend std::ostream& operator<<(std::ostream& out, const ApplySystFlags_t &f) {
    out << "useSystFlags(";
    f.printFlag(out,"noSyst",f.FNoSyst); out << ", ";
    f.printFlag(out,"pile-up",f.FPU); out << ", ";
    f.printFlag(out,"FSR",f.FFSR); 
    out << ")";
    return out;
  }
};

// -------------------------------------------------

//void loadCorrectionSystFromFile(int debug, const TString conf,
//				TString correctionKind, 
//				std::vector<TH2D*> &histos,
//				std::vector<std::string> &labels);

//void loadCorrectionSystFromFile(const TString &systFName,
//				TString correctionKind, 
//				std::vector<TH2D*> &histos,
//				std::vector<std::string> &labels);

void printSystTable(std::ostream& out, TString correctionKind, const TH2D* hFsrSyst, const TH2D *hPUSyst, int exponent=1);

inline
void printSystTable(TString correctionKind, const TH2D* hFsrSyst, const TH2D *hPUSyst, int exponent=1) { printSystTable(std::cout,correctionKind,hFsrSyst,hPUSyst,exponent); }

// -------------------------------------------------

int calcCorrectionSyst(int debug, const TString conf, 
		       TString correctionKind, TString fieldName, 
		       const ApplySystFlags_t &flagsPlus,
		       const ApplySystFlags_t &flagsMinus,
		       int printTable,
		       std::vector<TH2D*> *resHistos=NULL,
		       int *save=NULL) {
  DYTools::TRunMode_t runMode=DebugInt2RunMode(debug);

  int fsr5plus= flagsPlus.fsr();
  int fsr5minus=flagsMinus.fsr();
  int pu5plus=  flagsPlus.pu();
  int pu5minus= flagsMinus.pu();
  std::cout << "flagsPlus=" << flagsPlus << ", and flagsMinus=" << flagsMinus << ", decoded as FSR: " << fsr5plus << "," << fsr5minus << ", PU: " << pu5plus << "," << pu5minus << "\n";

  InputFileMgr_t inpMgr0;
  if (!inpMgr0.Load(conf) ||
      !inpMgr0.KeepOnlySample(-1)) return retCodeError;
  // no energy correction for this evaluation
  inpMgr0.removeMCSampleInfo();
  inpMgr0.clearEnergyScaleTag();


  DYTools::TSystematicsStudy_t systModeV[4]={ DYTools::FSR_5plus, DYTools::FSR_5minus, DYTools::PILEUP_5plus, DYTools::PILEUP_5minus };

  // EventSelector modifies info in inpMgr, therefore we need to keep all versions
  std::vector<InputFileMgr_t*> inpMgrV;
  for (int i=0; i<4; ++i) {
    InputFileMgr_t *mgr=new InputFileMgr_t(inpMgr0);
    inpMgrV.push_back(mgr);
    // Construct eventSelector, update mgr and plot directory
    EventSelector_t evtSel(*mgr,runMode,systModeV[i],
			   "","", EventSelector::_selectDefault);
    //mgr->Print();
  }

  // no syst 
  EventSelector_t eventSelector(inpMgr0,runMode,DYTools::NO_SYST,
				"","", EventSelector::_selectDefault);
    
  TString outFName0=inpMgr0.correctionFullFileName(correctionKind,DYTools::NO_SYST,0);
  TH2D *h2Base=LoadHisto2D(fieldName,outFName0);
  int ok=(h2Base!=NULL) ? 1:0;
  if (ok) {
    TString hBaseName=Form("h%s",correctionKind.Data());
    h2Base->SetName(hBaseName);
  }
    
  TH2D *h2FSR5plus=NULL, *h2FSR5minus=NULL;
  TH2D *h2PU5plus=NULL, *h2PU5minus=NULL;
  TH2D *h2FSRsyst=NULL, *h2PUsyst=NULL;

  if (ok && (fsr5plus || fsr5minus)) {
    if (ok && fsr5plus ) { 
      h2FSR5plus= LoadHisto2D(fieldName,inpMgrV[0]->correctionFullFileName(correctionKind,DYTools::FSR_5plus,0));
      if (!h2FSR5plus) ok=0;
    }
    if (ok && fsr5minus) {
      h2FSR5minus= LoadHisto2D(fieldName,inpMgrV[1]->correctionFullFileName(correctionKind,DYTools::FSR_5minus,0));
      if (!h2FSR5minus) ok=0;
    }
    if (ok) {
      if (0) {
	std::cout << "check 1st bin\n";
	std::cout << "base=" << h2Base->GetBinContent(1,1) << ", h2FSR5plus="  << h2FSR5plus->GetBinContent(1,1) << ", h2FSR5minus=" << h2FSR5minus->GetBinContent(1,1) << "\n";
      }
      TString hFSRsystName=Form("h%s_FSRsyst",correctionKind.Data());
      // 1 - take the difference between the variants as well
      h2FSRsyst= getRelDifference(h2Base,hFSRsystName,1,h2FSR5plus,h2FSR5minus);
      if (!h2FSRsyst) ok=0;
      //printHisto(h2FSRsyst);
    }
  }
  
  if (ok && (pu5plus || pu5minus)) {
    if (ok && pu5plus ) { 
      h2PU5plus= LoadHisto2D(fieldName,inpMgrV[2]->correctionFullFileName(correctionKind,DYTools::PILEUP_5plus,0));
      if (!h2PU5plus) ok=0;
    }
    if (ok && pu5minus) {
      h2PU5minus= LoadHisto2D(fieldName,inpMgrV[3]->correctionFullFileName(correctionKind,DYTools::PILEUP_5minus,0));
      if (!h2PU5minus) ok=0;
    }
    if (ok) {
      TString hPUsystName=Form("h%s_PileUpSyst",correctionKind.Data());
      // 1 - take the difference between the variants as well
      h2PUsyst= getRelDifference(h2Base,hPUsystName,1,h2PU5plus,h2PU5minus);
      if (!h2PUsyst) ok=0;
      //printHisto(h2PUsyst);
    }
  }

  if (ok) {
    ok=retCodeOk;
    if (resHistos!=NULL) {
      resHistos->reserve(3);
      resHistos->push_back(h2Base);
      if (h2FSRsyst) resHistos->push_back(h2FSRsyst);
      if (h2PUsyst) resHistos->push_back(h2PUsyst);
    }
  }

  if (ok && printTable) {
    printSystTable(correctionKind,h2FSRsyst,h2PUsyst,(printTable==2) ? 0:1);
  }


  if ((ok==retCodeOk) && resHistos && save && (*save==1)) {
    TString systFName=inpMgr0.correctionSystFullFileName(correctionKind,DYTools::NO_SYST,0);
    TFile fout(systFName,"recreate");
      int res=fout.IsOpen();
      if (res) res=saveVec(fout,*resHistos);
      if (res && h2FSR5plus) res=saveHisto(fout,h2FSR5plus,"details","FSR5plus");
      if (res && h2FSR5minus) res=saveHisto(fout,h2FSR5minus,"details","FSR5minus");
      if (res && h2PU5plus) res=saveHisto(fout,h2PU5plus,"details","PU5plus");
      if (res && h2PU5minus) res=saveHisto(fout,h2PU5minus,"details","PU5minus");
      writeBinningArrays(fout);
      fout.Close();
      std::cout << "systematics saved to file <" << fout.GetName() << "> with res=" << res << "\n";
      if (!res) ok=retCodeError;
  }
  else std::cout << "systematics was not saved\n";
    
  return ok;
}

// -------------------------------------------------
// -------------------------------------------------

int loadCorrectionSystFromFile(const TString &systFName,
			       TString correctionKind,
			       const ApplySystFlags_t &flags,
			       std::vector<TH2D*> &histosV,
			       std::vector<std::string> &labelsV) {

  histosV.clear();
  labelsV.clear();

  TString hBase=Form("h%s",correctionKind.Data());
  TString hFSRsystName=Form("h%s_FSRsyst",correctionKind.Data());
  TString hPUsystName=Form("h%s_PileUpSyst",correctionKind.Data());

  TFile fin(systFName,"read");
  TH2D *h2Base=LoadHisto2D(fin,hBase,"",1);
  if (!h2Base) {
    fin.Close();
    std::cout << "error in plotCorrectionSystFromFile\n";
    return retCodeError;
  }
  TH2D *h2FSR=(flags.fsr()) ? LoadHisto2D(fin,hFSRsystName,"",0) : NULL;
  TH2D *h2PU =(flags.pu() ) ? LoadHisto2D(fin,hPUsystName,"",0) : NULL;
  TH2D *h2FSR5plus=(flags.fsr()) ? LoadHisto2D(fin,"FSR5plus","details",0) : NULL;
  TH2D *h2FSR5minus=(flags.fsr()) ? LoadHisto2D(fin,"FSR5minus","details",0) : NULL;
  TH2D *h2PU5plus =(flags.pu()) ? LoadHisto2D(fin,"PU5plus","details",0) : NULL;
  TH2D *h2PU5minus=(flags.pu()) ? LoadHisto2D(fin,"PU5minus","details",0) : NULL;
  fin.Close();

  histosV.reserve(7);
  histosV.push_back(h2Base);
  histosV.push_back(h2FSR);
  histosV.push_back(h2PU);
  histosV.push_back(h2FSR5plus);
  histosV.push_back(h2FSR5minus);
  histosV.push_back(h2PU5plus);
  histosV.push_back(h2PU5minus);

  labelsV.reserve(7);
  labelsV.push_back(correctionKind.Data());
  labelsV.push_back("FSR syst");
  labelsV.push_back("Pile-up syst");
  labelsV.push_back("FSR 5% plus");
  labelsV.push_back("FSR 5% minus");
  labelsV.push_back("Pile-up 5% plus");
  labelsV.push_back("Pile-up 5% minus");

  return retCodeOk;
}

// -------------------------------------------------


int loadCorrectionSystFromFile(int debug, const TString conf,
			       TString correctionKind,
			       const ApplySystFlags_t &flags,
			       std::vector<TH2D*> &histos,
			       std::vector<std::string> &labels) {
  DYTools::TRunMode_t runMode=DebugInt2RunMode(debug);
  //DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;

  InputFileMgr_t inpMgr0;
  if (!inpMgr0.Load(conf)) return retCodeError;

  // Needed as long as it modifies the inpMgr0 fields
  EventSelector_t eventSelector(inpMgr0,runMode,DYTools::NO_SYST,
				"","", EventSelector::_selectDefault);
  
  TString systFName=inpMgr0.correctionSystFullFileName(correctionKind,DYTools::NO_SYST,0);
  
  if (loadCorrectionSystFromFile(systFName,correctionKind,flags,histos,labels)!=retCodeOk) {
    std::cout << "error in loadCorrectionSystFromFile(conf,correctionKind)\n";
    return retCodeError;
  }
  return retCodeOk;
}


// -------------------------------------------------
// -------------------------------------------------

void printSystTable(std::ostream& out, TString correctionKind, const TH2D* hFsrSyst, const TH2D *hPUSyst, int exponent) {
  if (out != std::cout) HERE("entered printSystTable");
  out << "systematics for " << correctionKind << "\n";

  if (hFsrSyst) printHisto(hFsrSyst);
  if (hPUSyst) printHisto(hPUSyst);
  if (!hFsrSyst && !hPUSyst) {
    const char *msg="\nERROR: both FSR and PU histograms are null\n\n";
    out << msg;
    if (out != std::cout) std::cout << msg;
    return;
  }

  const TH2D *h2=(hFsrSyst!=NULL) ? hFsrSyst : hPUSyst;
  std::string col1Name=(hFsrSyst!=NULL) ? "FSRsyst" : "PileUpSyst";
  const TH2D *h2two=(hFsrSyst!=NULL) ? hPUSyst : NULL;
  std::string col2Name= (h2two) ? "PileUpSyst" : "";

  const char *formatRange= " %6.1f-%6.1f  %5.1f-%5.1f  ";
  const std::string formatRaw   = (exponent) ? "  %7.2e" : " %lf";
  const std::string formatRelRaw= "    %4.2lf";
  std::string format ;
  std::string format2=formatRaw;

  char buf[100];
  
  for (int rel=0; rel<2; rel++) {
    if (rel==0) {
      out << "\n Absolute numbers\n";
      format=formatRaw + formatRaw;
    }
    else {
      out << "\n Relative numbers (%)\n";
      col1Name+="(rel,%)";
      col2Name+="(rel,%)";
      format =std::string("(") + formatRaw + std::string(")*100") + formatRelRaw;
      format2=formatRelRaw;
    }
    out << " mass range   rapidity range   central val.  ";
    out << col1Name;
    if (h2two) out << "    " << col2Name;
    out << "\n";

    for (int ibin=1; ibin<=h2->GetNbinsX(); ++ibin) {
      double x=h2->GetBinLowEdge(ibin);
      double w=h2->GetBinWidth(ibin);
      for (int jbin=1; jbin<=h2->GetNbinsY(); ++jbin) {
	double y=h2->GetYaxis()->GetBinLowEdge(jbin);
	double h=h2->GetYaxis()->GetBinWidth(jbin);
	double val =h2->GetBinContent(ibin,jbin);
	double err1=h2->GetBinError(ibin,jbin);
	double err2=(h2two) ? h2two->GetBinError(ibin,jbin) : 0.0;
	if (rel && (val!=double(0))) { err1/=(0.01*val); err2/=(0.01*val); }
	sprintf(buf,formatRange,x,x+w,y,y+h);
	out << buf;
	sprintf(buf,format.c_str(), val,err1);
	out << buf;
	if (h2two) {
	  sprintf(buf,format2.c_str(), err2);
	  out << buf;
	}
	out << "\n";
      }
    }
  }
  return;
}


// -------------------------------------------------

#endif
