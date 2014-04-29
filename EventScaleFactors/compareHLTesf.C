#include <TROOT.h>
#include <TCanvas.h>
#include "../Include/MitStyleRemix.hh"
#include "../Include/ComparisonPlot.hh"
#include "../Include/MyTools.hh"
//#include "../Include/UnfoldingTools.hh"


const int work2D=1;
int iMass=0;

// --------------------------------------------------------------

TH1D* getHistoFromFile(const char *fname, 
		       const char *varName, const char *varErrName,
		       int isMatrix,
		       const char *newHistoName) {
   TFile miofile(fname,"read");
   const int nMBins=DYTools::nMassBins;
   const int nYBins=DYTools::nYBinsMax;
   TMatrixD var(nMBins,nYBins), varErr(nMBins,nYBins);

   if (isMatrix) {
     TMatrixD *x=(TMatrixD*) miofile.Get(varName);
     TMatrixD *xErr=(TMatrixD*) miofile.Get(varErrName);
     if (!x || !xErr) {
       std::cout << "failed to get <" << varName << "> or <" << varErrName << ">\n";
       return NULL;
     }
     var= *x; 
     varErr= *xErr;
     delete x;
     delete xErr;
   }
   else {
     TVectorD *x=(TVectorD*) miofile.Get(varName);
     TVectorD *xErr=(TVectorD*) miofile.Get(varErrName);
     if (!x || !xErr) {
       std::cout << "failed to get <" << varName << "> or <" << varErrName << ">\n";
       return NULL;
     }
     deflattenMatrix(*x, var);
     deflattenMatrix(*xErr, varErr);
     delete x;
     delete xErr;
   }

   TH1D* histo=NULL;

   if (work2D==0) {
     HERE("extracting mass dependence");
     histo=
       extractMassDependence(newHistoName,"",
			     var,varErr,
			     0,0,0);
     HERE("done");
     /*
     std::cout << "variable=\n";
     var.Print();
     std::cout << "varErr=\n";
     varErr.Print();
     std::cout << "histo=\n";
     printHisto(histo,0);
     */
   }
   else {
     HERE("extracting rapidity dependence");
     histo=
       extractRapidityDependence(newHistoName,"",
			     var,varErr,
			     iMass,0);
     HERE("done");
   }


   miofile.Close();
   //histo->Sumw2();
   return histo;
}


// --------------------------------------------------------------

TH1D* getHistoFromFile(const std::string &fname,
		       const std::string &varName,
		       const std::string &varErrName,
		       int isMatrix,
		       const std::string &newHistoName) {
  return getHistoFromFile(fname.c_str(), 
			  varName.c_str(), 
			  varErrName.c_str(),
			  isMatrix,
			  newHistoName.c_str());
}

// --------------------------------------------------------------



// --------------------------------------------------------------

void compareHLTesf(int analysisIs2D) {

  if (!DYTools::setup(analysisIs2D)) {
    std::cout << "failed to initialize the analysis\n";
    return retCodeError;
  }

  if (DYTools::study2D!=work2D) {
    std::cout << "DYTools::study2D=" << DYTools::study2D << ", work2D=" << work2D << "\n";
    return;
  }


  const std::string pathA="../root_files/constants/DY_j22_19789pb/";

  std::string fName0= pathA + std::string("scale_factors_1D.root");
  std::string fName1= pathA + std::string("scale_factors_asymHLT_1D.root");
  if (work2D) {
    TString tmp;
    tmp=fName0.c_str();  tmp.ReplaceAll("1D","2D");  fName0=tmp.Data();
    tmp=fName1.c_str();  tmp.ReplaceAll("1D","2D");  fName1=tmp.Data();
  }

  int performChk=0;
  std::vector<std::string> tags, labels;
  tags.push_back(""); labels.push_back("truncated HLT eff.");
  tags.push_back("_HLTlegs"); labels.push_back("full HLT eff.");

  std::vector<std::string> varNames, varErrNames, newHistoNames;
  std::vector<std::string> yAxisLabel;
  std::vector<int> isMatrix;
  if (!performChk) {
    varNames.push_back("scaleHltFactor_FIArray");
    varErrNames.push_back("scaleHltFactorErr_FIArray");
    isMatrix.push_back(0);
    newHistoNames.push_back("esfHLT");
    yAxisLabel.push_back("HLT scale factor");
  }

  varNames.push_back("scaleFactor");
  varErrNames.push_back("scaleFactorErr");
  isMatrix.push_back(1);
  newHistoNames.push_back("esfTot");
  yAxisLabel.push_back("#rho_{data/MC}");

  /*
  if (0) { // Test flatIndexArr
    varNames.push_back("scaleFactorFlatIdxArray");
    varErrNames.push_back("scaleFactorErrFlatIdxArray");
    isMatrix.push_back(0);
    newHistoNames.push_back("esfTot_chkFI");
  }
  */

  std::vector<ComparisonPlot_t*> compPlotV;
  std::vector<TCanvas*> canvasV;

  for (unsigned int i=0; i<varNames.size(); ++i) {
    //if (i==0) continue;
    std::cout << "processing varName=" << varNames[i] << std::endl;

    int startIMass=(work2D) ? 1:0;
    int endIMass=(work2D) ? DYTools::nMassBins : 1;
    //endIMass=startIMass+1;
    for (int iMassIdx=startIMass; iMassIdx<endIMass; ++iMassIdx) {
      iMass=iMassIdx;
      std::string massStr=(work2D) ? Form("mass_%1.0lf_%1.0lf",DYTools::massBinLimits[iMassIdx],DYTools::massBinLimits[iMassIdx+1]) : "";

      std::string name0=newHistoNames[i] + tags[0] + massStr;
      TH1D *histo0=getHistoFromFile(fName0, varNames[i],varErrNames[i],isMatrix[i],name0);
      return;

      std::string name1=newHistoNames[i] + tags[1] + massStr;
      TH1D *histo1=getHistoFromFile(fName1, varNames[i],varErrNames[i],isMatrix[i],name1);

      TH1D *histoChk0=NULL, *histoChk1=NULL;
      TH1D *histoTest0=NULL;
    /*
      if (performChk) {
      const char* sfName="scaleFactorFlatIdxArray";
      const char* sfErrName="scaleFactorErrFlatIdxArray";
      std::string name3=newHistoNames[i] + tags[3];
      histoChk0=getHistoFromFile(fNameChk0, sfName,sfErrName,!isMatrix[i],name3);
      std::string name4=newHistoNames[i] + tags[4];
      histoChk1=getHistoFromFile(fNameChk1, sfName,sfErrName,!isMatrix[i],name4);
      std::string name5=newHistoNames[i] + tags[5];
      histoTest0=getHistoFromFile(fNameTest0, sfName,sfErrName,!isMatrix[i],name5);
    }
    */
      
    //std::string name2=newHistoNames[i] + tags[2];
      TH1D *histo2=NULL;
    //histo2=getHistoFromFile(fName2, varNames[i],varErrNames[i],isMatrix[i],name2);

      TString canvName=TString("canv_") + TString(varNames[i].c_str());
      if (work2D) canvName.Append(TString("_") + massStr);
      TCanvas *canvas=MakeCanvas(canvName,canvName,600,800);
      canvasV.push_back(canvas);
      
      if (!histo0) std::cout << "histo0 is null\n";
      HERE("creating comparisonPlot obj");

      TString cpName=TString("compPlot_") + TString(varNames[i].c_str());
      ComparisonPlot_t *cp=
	new ComparisonPlot_t(ComparisonPlot_t::_ratioPlain,
			     cpName,
			     (work2D) ? massStr.c_str() : "",
			     (work2D) ? "|y|" : "m(e^{+}e^{-}) [GeV]",
			     yAxisLabel[i].c_str(),
			     "ratio");
      cp->SetRatioTitleOffset(0.5);
      cp->SetRatioYRange(0.96,1.04);
      cp->SetYRange(0.9,1.1);

      if (!work2D) {
	// 1D case
	cp->SetLogx(1);
	if (i==0) {
	  cp->SetRatioYRange(0.99,1.01);
	  cp->SetYRange(0.99,1.01);
	  cp->SetRatioTitleOffset(0.6);
	}
	else if (i==1) {
	  cp->SetRatioYRange(0.96,1.04);
	  cp->SetYRange(0.86,1.0);
	}
      }
      else {
	// 2D case
	if (i==0) {
	  cp->SetYRange(0.98,1.02);
	  cp->SetRatioYRange(0.99,1.01);
	  cp->SetRatioTitleOffset(0.6);
	}
	else {
	  cp->SetYRange(0.75, 1.15);
	}
      }
      
      cp->SetPrintRatios();
      compPlotV.push_back(cp);
      
      cp->Prepare2Pads(canvas,0.01);
      histo1->SetMarkerStyle(24);
      cp->AddHist1D(histo1,labels[1],"LPE1",46,1);
      cp->AddHist1D(histo0,labels[0],"LPE1",9,1); // 38
      if (histo2) cp->AddHist1D(histo2,labels[2],"LPE2",kRed,1);
      if (histoChk0) cp->AddHist1D(histoChk0,labels[3],"LPE1",kGreen+1,1);
      if (histoChk1) cp->AddHist1D(histoChk1,labels[4],"LPE1",kViolet,1);
      if (histoTest0) cp->AddHist1D(histoTest0,labels[5],"LPE1",kBlack,1);
      cp->Draw(canvas,false,"png");
      cp->ChangeLegendPos(-0.05,0.,0.05,0.);

      SaveCanvas(canvas,canvas->GetName(),"dir-compareESF/");
    }
  }

  return;
}
