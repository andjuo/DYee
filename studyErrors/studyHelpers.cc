#include "studyHelpers.hh"


// -------------------------------------------------------------
// -------------------------------------------------------------

TString valTypeName(TValueType_t valType) {
  TString name;
  switch(valType) {
  case _valTBYieldGen: name="trueBkgGeneratedRndValue"; break;
  case _valTBYieldTot: name="trueBkgTotalRndValue"; break;
  case _valYieldTot: name="rndSignalYield"; break;
  case _valUnfYield: name="unfYield"; break;
  case _valCS: name="crossSection"; break;
  case _valLast: name="lastIndexValue"; break;
  default:
    std::cout << "valTypeName cannot be determined for valType="
	      << int(valType) << "\n";
    name="__UNKNOWN__";
  }
  return name;
}

// -------------------------------------------------------------

TString sigYieldValTypeName(TSigYieldValueType_t valType) {
  TString name;
  switch(valType) {
  case _syval_sigYieldOrig: name="sigYieldOrigRnd"; break;
  case _syval_sigYieldGen: name="sigYieldGeneratedRnd"; break;
  case _syval_sigYieldActual: name="sigYieldActual"; break;
  case _syval_unfYieldOrig: name="unfYieldOrig"; break;
  case _syval_unfYield: name="unfYieldActual"; break;
  case _syval_csVal: name="crossSection"; break;
  case _syvalLast: name="lastIndexValue"; break;
  default:
    std::cout << "sigYieldValTypeName cannot be determined for "
	      << "sigYieldValType="
	      << int(valType) << "\n";
    name="__UNKNOWN__";
  }
  return name;
}

// -------------------------------------------------------------
// -------------------------------------------------------------

int calculateSidedErrors(const TH1D* hVals, double centralValue,
			 double &errPos, double &errNeg, double &errRMS) {
  errPos=0.; errNeg=0.; errRMS=0.;
  if (!hVals) return 0;
  double sumPosW2=0., sumNegW2=0., sumTotW2=0.;
  double sumPos=0., sumNeg=0., sumTot=0.;
  int countPos=0, countNeg=0;
  for (int ibin=1; ibin<=hVals->GetNbinsX(); ++ibin) {
    double val=hVals->GetBinCenter(ibin);
    double weight=hVals->GetBinContent(ibin);
    double dev= val - centralValue;
    sumTot += weight;
    sumTotW2 += weight*dev*dev;
    if (dev>0) {
      countPos++;
      sumPos += weight;
      sumPosW2 += weight*dev*dev;
    }
    else {
      countNeg++;
      sumNeg += weight;
      sumNegW2 += weight*dev*dev;
    }
  }
  //HERE("countPos=%d",countPos); std::cout << " " << sumPosW2 << "/" << sumPos << "\n";
  //std::cout << "countNeg=" << countNeg << ", " << sumNegW2 << "/" << sumNeg << "\n";
  if (countPos && (sumPos!=double(0))) sumPosW2= sqrt(sumPosW2/sumPos);
  if (countNeg && (sumNeg!=double(0))) sumNegW2= sqrt(sumNegW2/sumNeg);
  if (sumTot!=double(0)) errRMS =sqrt(sumTotW2/sumTot);
  errPos= sumPosW2;
  errNeg= sumNegW2;
  return 1;
}

// -------------------------------------------------------------

void printSidedErrors(const TH1D* hVals, double centralValue) {
  if (centralValue==-9999.) centralValue= hVals->GetMean();
  double errPos=0, errNeg=0, rms=0;
  int ok=calculateSidedErrors(hVals,centralValue,errPos,errNeg,rms);
  std::cout << "sided errors of " << hVals->GetName() << ":\n";
  std::cout << "  - from histo " << hVals->GetMean() << " +- "
	    << hVals->GetRMS() << "\n";
  std::cout << "  - from subroutine ";
  if (!ok) std::cout << " (calculation FAILED) ";
  std::cout << centralValue << " +"
	    << errPos << "-" << errNeg << " (rms=" << rms << ")\n";
}

// -------------------------------------------------------------
// -------------------------------------------------------------

void StudyInfo_t::Clear() {
  fRawHistos40.clear();
  fRawHistos41.clear();
  this->Zero();
}

// -------------------------------------------------------------

void StudyInfo_t::Zero() {
  for (TValueType_t i=_valTBYieldGen; i<_valLast; next(i)) {
    fErr40pos[i]=0.;
    fErr40neg[i]=0.;
    fErr40[i]=0.;
    fErr41pos[i]=0.;
    fErr41neg[i]=0.;
    fErr41[i]=0.;
  }
}

// -------------------------------------------------------------

double StudyInfo_t::getCenter(int ibin, TValueType_t valType,
			      int &ierr) const
{
  double center=0.;
  ierr=0;
  switch(valType) {
  case _valTBYieldGen: {
    const std::vector<TH1D*> *hRawV=
      (ibin==40) ? &fRawHistos40 : &fRawHistos41;
    center= (*hRawV)[valType]->GetMean();
  }
    break;
  case _valTBYieldTot: center= fH2TrueBkg->GetBinContent(ibin,1); break;
  case _valYieldTot: center= fH2SigYield->GetBinContent(ibin,1); break;
  case _valUnfYield: center= fH2UnfYield->GetBinContent(ibin,1); break;
  case _valCS: center= fH2MainCS->GetBinContent(ibin,1); break;
  default:
    std::cout << "getCenter: cannot get center for valType="
	      << int(valType) << "\n";
    ierr=1;
    return 0;
  }
  return center;
}

// -------------------------------------------------------------

int StudyInfo_t::init(TString fname, TString histoTag) {
  TFile fin(fname,"read");
  if (!fin.IsOpen()) {
    std::cout << "StudyInfo::init: failed to open the file <"
	      << fin.GetName() << ">\n";
    return 0;
  }
  std::cout << "file <" << fin.GetName() << "> opened\n";
  fRawHistos40.reserve(_valLast);
  fRawHistos40.clear();
  fRawHistos41.reserve(_valLast);
  fRawHistos41.clear();
  for (int ibin=40; ibin<=41; ++ibin) {
    HERE("ibin=%d",ibin);
    std::vector<TH1D*> *hV= (ibin==40) ? &fRawHistos40 : &fRawHistos41;
    for (TValueType_t i=_valTBYieldGen; i<_valLast; next(i)) {
      std::cout << "get histos for ibin=" << ibin << " valType=" << i << "\n";
      TString fieldBase;
      switch(i) {
      case _valTBYieldGen: fieldBase=Form("h1TrueBkgValAt%d",ibin); break;
      case _valTBYieldTot: fieldBase=Form("h1TrueBkgTotalValAt%d",ibin); break;
      case _valYieldTot: fieldBase=Form("h1SigYieldAt%d",ibin); break;
      case _valUnfYield: fieldBase=Form("h1UnfValAt%d",ibin); break;
      case _valCS: fieldBase=Form("h1CSValAt%d",ibin); break;
      default:
	std::cout << "init: cannot load value of type"
		  << int(i) << "\n";
	return 0;
      }

      TString newName=fieldBase + histoTag;
      TH1D *h1=new TH1D(fieldBase,newName,1,0.,1.);
      if (!loadHisto(fin,&h1,"")) {
	std::cout << "failed to load histo <" << fieldBase << ">\n";
	fin.Close();
	return 0;
      }
      h1->SetName(newName);
      hV->push_back(h1);
    }
  }
  HERE("individual histos loaded");

  TString newName;
  fH2SigYield=new TH2D("signalYieldDDbkg","",1,0.,1.,1,0.,1.);
  fH2UnfYield=new TH2D("unfoldedYieldDDbkg","",1,0.,1.,1,0.,1.);
  fH2TrueBkg =new TH2D("true2eBackgroundFromData","",1,0.,1.,1,0.,1.);
  fH2MainCS  =new TH2D("csYieldDDbkg","",1,0.,1.,1,0.,1.);
  fH2AvgCS=new TH2D("csAvgDistr","",1,0.,1.,1,0.,1.);
  if (!fH2SigYield || !fH2UnfYield || !fH2TrueBkg || !fH2MainCS || !fH2AvgCS) {
    std::cout << "failed to create storage histos\n";
    return 0;
  }
  if (!loadHisto(fin,&fH2SigYield,"") ||
      !loadHisto(fin,&fH2UnfYield,"") ||
      !loadHisto(fin,&fH2TrueBkg,"") ||
      !loadHisto(fin,&fH2MainCS,"") ||
      !loadHisto(fin,&fH2AvgCS,"")) {
    std::cout << "failed to load reference histos\n";
    fin.Close();
    return 0;
  }
  AppendToHistoName(fH2SigYield,histoTag,1);
  AppendToHistoName(fH2UnfYield,histoTag,1);
  AppendToHistoName(fH2TrueBkg,histoTag,1);
  AppendToHistoName(fH2MainCS,histoTag,1);
  AppendToHistoName(fH2AvgCS,histoTag,1);
  fin.Close();
  return 1;
}

// -------------------------------------------------------------

int StudyInfo_t::getSidedErrors(TString fname, TString histoTag) {
  this->Zero();
  if (!this->init(fname,histoTag)) {
    std::cout << "getSidedErrors: failed to load info\n";
    return 0;
  }

  for (TValueType_t i=_valTBYieldGen; i<_valLast; next(i)) {
    //if (i!=_valTBYieldTot) { HERE("\t\tskipping"); continue; }
    for (int ibin=40; ibin<=41; ++ibin) {
      const std::vector<TH1D*> *hRawV=
	(ibin==40) ? &fRawHistos40 : &fRawHistos41;
      double *errPos= (ibin==40) ? fErr40pos : fErr41pos;
      double *errNeg= (ibin==40) ? fErr40neg : fErr41neg;
      double *errRMS= (ibin==40) ? fErr40 : fErr41;
      int ierr=0;
      double center= this->getCenter(ibin, i, ierr);
      if (ierr) {
	std::cout << "error in getSidedErrors at i=" << int(i) << "\n";
	return 0;
      }

      //printHisto((*hRawV)[i]);
      std::cout << "center=" << center << "\n";
      if (!calculateSidedErrors((*hRawV)[i],center,
				errPos[i],errNeg[i],errRMS[i])) {
	std::cout << "error in getSidedErrors at i=" << int(i) << "\n";
	return 0;
      }
    }
  }

  return 1;
}

// -------------------------------------------------------------

void StudyInfo_t::PrintValues(int ibin, TValueType_t valType) const {
  int ierr=0;
  double center= this->getCenter(ibin,valType,ierr);
  double errPos= this->getErrPos(ibin,valType);
  double errNeg= this->getErrNeg(ibin,valType);
  double errRMS= this->getErrRMS(ibin,valType);
  std::cout << " values of " << valType << " at ibin=" << ibin << ":\n";
  std::cout << "  " << center;
  if (ierr) std::cout << "(badVal)";
  std::cout << " +- " << errRMS;
  std::cout << " (asymErr: + " << errPos << " - " << errNeg << ")\n";
}

// -------------------------------------------------------------
// -------------------------------------------------------------

void SigErrStudyInfo_t::Clear() {
  fRawHistos41.clear();
  this->Zero();
}

// -------------------------------------------------------------

void SigErrStudyInfo_t::Zero() {
  for (TSigYieldValueType_t i=_syval_sigYieldOrig; i<_syvalLast; next(i)) {
    fErr41[i]=0.;
  }
}

// -------------------------------------------------------------
/*
double SigErrStudyInfo_t::getCenter(TSigYieldValueType_t valType,
			      int &ierr) const
{
  double center=0.;
  ierr=0;
  switch(valType) {
  case _syval_sigYieldOrig:
  case _syval_sigYieldGen:
  case _syval_sigYieldActual:
  case _syval_unfYieldOrig:
  case _syval_unfYield:
  case _syval_csVal:
    center=fRawHistos41[valType]->GetMean();
    break;
  default:
    std::cout << "getCenter: cannot get center for valType="
	      << int(valType) << "\n";
    ierr=1;
    return 0;
  }
  return center;
}
*/
// -------------------------------------------------------------

int SigErrStudyInfo_t::init(TString fname, TString histoTag) {
  TFile fin(fname,"read");
  if (!fin.IsOpen()) {
    std::cout << "SigErrStudyInfo::init: failed to open the file <"
	      << fin.GetName() << ">\n";
    return 0;
  }
  std::cout << "file <" << fin.GetName() << "> opened\n";
  fRawHistos41.reserve(_syvalLast);
  fRawHistos41.clear();
  std::vector<TH1D*> *hV= &fRawHistos41;
  for (TSigYieldValueType_t i=_syval_sigYieldOrig; i<_syvalLast; next(i)) {
    std::cout << "get histos for valType=" << i << "\n";
    TString fieldBase;
    switch(i) {
    case _syval_sigYieldOrig: fieldBase="h1SigYieldOrigAt41"; break;
    case _syval_sigYieldGen: fieldBase="h1SigYieldGenAt41"; break;
    case _syval_sigYieldActual: fieldBase="h1SigYieldActualAt41"; break;
    case _syval_unfYieldOrig: fieldBase="h1UnfValAt41orig"; break;
    case _syval_unfYield: fieldBase="h1UnfValAt41"; break;
    case _syval_csVal: fieldBase="h1CSValAt41"; break;
    default:
      std::cout << "init: cannot load value of type"
		<< int(i) << "\n";
      return 0;
    }

    TString newName=fieldBase + histoTag;
    TH1D *h1=new TH1D(fieldBase,newName,1,0.,1.);
    if (!loadHisto(fin,&h1,"")) {
      std::cout << "failed to load histo <" << fieldBase << ">\n";
      fin.Close();
      return 0;
    }
    h1->SetName(newName);
    hV->push_back(h1);
  }
  HERE("individual histos loaded");

  TString newName;
  fH2SigYield=new TH2D("signalYieldDDbkg","",1,0.,1.,1,0.,1.);
  fH2UnfYield=new TH2D("unfoldedYieldDDbkg","",1,0.,1.,1,0.,1.);
  fH2MainCS  =new TH2D("csYieldDDbkg","",1,0.,1.,1,0.,1.);
  fH2AvgCS=new TH2D("csAvgDistr","",1,0.,1.,1,0.,1.);
  if (!fH2SigYield || !fH2UnfYield || !fH2MainCS || !fH2AvgCS) {
    std::cout << "failed to create storage histos\n";
    return 0;
  }
  if (!loadHisto(fin,&fH2SigYield,"") ||
      !loadHisto(fin,&fH2UnfYield,"") ||
      !loadHisto(fin,&fH2MainCS,"") ||
      !loadHisto(fin,&fH2AvgCS,"")) {
    std::cout << "failed to load reference histos\n";
    fin.Close();
    return 0;
  }
  AppendToHistoName(fH2SigYield,histoTag,1);
  AppendToHistoName(fH2UnfYield,histoTag,1);
  AppendToHistoName(fH2MainCS,histoTag,1);
  AppendToHistoName(fH2AvgCS,histoTag,1);
  fin.Close();
  return 1;
}

// -------------------------------------------------------------

void SigErrStudyInfo_t::PrintValues(TSigYieldValueType_t valType) const {
  double center= this->getCenter(valType);
  double errRMS= this->getErrRMS(valType);
  std::cout << " values of " << valType << " at ibin=" << 41 << ":";
  std::cout << "  " << center;
  std::cout << " +- " << errRMS;
  std::cout << " (" << errRMS*100/center << "%)";
  std::cout << "\n";
}

// -------------------------------------------------------------

void SigErrStudyInfo_t::PrintValuesAll() const {
  for (TSigYieldValueType_t i=_syval_sigYieldOrig; i<_syvalLast; next(i)) {
    this->PrintValues(i);
  }
}

// -------------------------------------------------------------
// -------------------------------------------------------------
