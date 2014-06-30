#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "studyHelpers.hh"


// -------------------------------------------------------------
// -------------------------------------------------------------

int performSigYieldErrStudy() {

  const int markerCount=5;
  const TAttMarker markerV[markerCount] = { 
    TAttMarker(kBlack,24,1.),
    TAttMarker(kBlue,20,1.),
    TAttMarker(kRed+1,22,1.),
    TAttMarker(kGreen+1,4,1.),
    TAttMarker(kViolet,3,1.)
  };


  // define info to use
  std::vector<TString> fileTagV;
  std::vector<TString> explainV;
  fileTagV.reserve(20); explainV.reserve(20);
  fileTagV.push_back("_-1__"); explainV.push_back("default");
  if (0) {
    fileTagV.push_back("_-1_-mdfLastBinOnly_"); explainV.push_back("def.,mdfLastBin");
    fileTagV.push_back("_-1_-corrSYLastBin-mdfLastBinOnly_"); explainV.push_back("def.,corrSigYield");
    fileTagV.push_back("_-1_-corrSYLastBinUYLastBin-mdfLastBinOnly_"); explainV.push_back("def.,corrSig&UnfYield");
  }
  fileTagV.push_back("_100__"); explainV.push_back("Log-normal #mu, ver.A");
  fileTagV.push_back("_101__"); explainV.push_back("Log-normal #mu, ver.B");
  fileTagV.push_back("_102__"); explainV.push_back("Log-normal #mu, ver.C");
  fileTagV.push_back("_103__"); explainV.push_back("Log-normal #mu,#sigma ver.D");
  fileTagV.push_back("_105__"); explainV.push_back("Log-normal #mu,#sigma ver.E");
  fileTagV.push_back("_200__"); explainV.push_back("Gamma #mu, ver.A");
  fileTagV.push_back("_202__"); explainV.push_back("Gamma #mu, ver.B");
  //fileTagV.push_back("_300__"); explainV.push_back("truncated Gaussian, #mu");

  std::cout << "there are " << fileTagV.size() << " file tags\n";

  // construct file names
  std::vector<TString> fnamesV;
  fnamesV.reserve(fileTagV.size());
  for (unsigned int i=0; i<fileTagV.size(); ++i) {
    TString fname="sigYieldErrStudy" + fileTagV[i] + TString("nExps1000.root");
    fnamesV.push_back(fname);
  }

  // calculate the errors
  std::vector<SigErrStudyInfo_t*> infoV;
  infoV.reserve(fnamesV.size());

  for (unsigned int i=0; i<fnamesV.size(); ++i) {
    std::cout << "loading " << fnamesV[i] << "\n";
    SigErrStudyInfo_t *info=new SigErrStudyInfo_t();
    if (!info->init(fnamesV[i],fileTagV[i])) {
      std::cout << "failed to get info\n";
      return retCodeError;
    }
    infoV.push_back(info);
    info->PrintValuesAll();
  }

  // Make plots
  for (TSigYieldValueType_t iStudy=_syval_sigYieldOrig; 
       iStudy<_syvalLast; next(iStudy)) 
  {
    TString tag=sigYieldValTypeName(iStudy);
    ComparisonPlot_t *cp= new ComparisonPlot_t("cp",tag,
					       "events","count","ratio");
    cp->SetRefIdx(-111);

    for (unsigned int i=0; i<infoV.size(); ++i) {
      cp->AddHist1D(infoV[i]->rawHisto41(iStudy),explainV[i],
		    "LP",markerV[i%markerCount],
		    i/markerCount+1, 0,1);
    }

    TString canvName="cx_" + tag;
    TCanvas *cx= new TCanvas(canvName,canvName,800,800);
    cp->Draw1(cx);
    cx->Update();
  }


  return retCodeOk;
}
// -------------------------------------------------------------
// -------------------------------------------------------------
