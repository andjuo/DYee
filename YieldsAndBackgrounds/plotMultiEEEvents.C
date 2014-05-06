#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/ComparisonPlot.hh"
#include "../Include/InputFileMgr.hh"
#include "../Include/EventSelector.hh"
#include "../Include/colorPalettes.hh"

// -------------------------------------

int plotMultiEEEvents(TString confName,
	    DYTools::TRunMode_t runMode=DYTools::NORMAL_RUN,
	    DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST)
{

  if (!DYTools::setup(0)) {
    std::cout << "failed to setup the analysis\n";
    return retCodeError;
  }

  if (DYTools::study2D) {
    std::cout << "Macro is not ready for the 2D case\n";
    return retCodeError;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================

  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(confName)) return retCodeError;
  //mgr.Print();

  // no energy correction for this evaluation
  //inpMgr.clearEnergyScaleTag();

  // Construct eventSelector, update mgr and plot directory
  TString extraTag="";
  TString plotExtraTag;
  EventSelector_t evtSelector(inpMgr,runMode,systMode,
	       extraTag, plotExtraTag, EventSelector::_selectDefault);

  int createDir=0;
  TString yieldFullName= inpMgr.yieldFullFileName(-1,systMode,createDir);
  std::cout << "will work with yieldFile=<" << yieldFullName << ">\n";

  //--------------------------------------------------------------------------------------------------------------
  // Start job
  //==============================================================================================================

  vector<TH2D*> yields; // non-uniform rapidity grid taken into account
  vector<TH2D*> yieldsSameEventLead_BaseH2; // multi-dielectron events
  vector<TH2D*> yieldsSameEventTail_BaseH2; // multi-dielectron events
  vector<TH2D*> yieldsSameEventLead; // multi-dielectron events
  vector<TH2D*> yieldsSameEventTail; // multi-dielectron events

  // the main yield
  createBaseH2Vec(yields,"hYield_",inpMgr.sampleNames());
  // distribution of the same-event candidates
  createBaseH2Vec(yieldsSameEventLead_BaseH2,
		  "hYieldsSameEventLead_BaseH2_",inpMgr.sampleNames());
  createBaseH2Vec(yieldsSameEventTail_BaseH2,
		  "hYieldsSameEventTail_BaseH2_",inpMgr.sampleNames());

  TFile file(yieldFullName,"read");
  int res=1;
  if (res) res=checkBinningArrays(file);
  if (res) res=loadVec(file,yields,"yields");
  if (res) res=loadVec(file,yieldsSameEventLead_BaseH2,"yieldsSameEventLead");
  if (res) res=loadVec(file,yieldsSameEventTail_BaseH2,"yieldsSameEventTail");
  file.Close();
  if (!res) return retCodeError;

  if (!convertBaseH2actualVec(yieldsSameEventLead_BaseH2, yieldsSameEventLead,
			      "hYieldsSameEventLead_",inpMgr.sampleNames(),
			      1)) return retCodeError;
  if (!convertBaseH2actualVec(yieldsSameEventTail_BaseH2, yieldsSameEventTail,
			      "hYieldsSameEventTail_",inpMgr.sampleNames(),
			      1)) return retCodeError;


  std::vector<std::vector<TH1D*>*> hYields;
  if (!createMassProfileVec(yields,hYields,inpMgr.sampleNames()))
    return retCodeError;
  TH1D *hData=(TH1D*)(*hYields[0])[0]->Clone("hDataYield");
  TH1D *hMC  =(TH1D*)(*hYields[0])[1]->Clone("hMCYield");
  for (unsigned int ih=2; ih<hYields[0]->size(); ++ih) {
    hMC->Add((*hYields[0])[ih],1.);
  }
  //printHisto(hData);

  double ymin=0,ymax=1e5;
  TString titleBase="Reg.En., ";
  TString fileBase="fig-RegEn--";
  if (systMode==DYTools::UNREGRESSED_ENERGY) {
    titleBase="UnReg.En, ";  fileBase="fig-UnRegEn--"; 
  }
  else if (systMode==DYTools::APPLY_ESCALE) {
    titleBase="Reg.En.Mdf., "; fileBase="fig-RegEnMdf--";
  }
  for (int plot_case=0; plot_case<1; ++plot_case) {
    TString plotTag;
    TString titleTag;
    switch(plot_case) {
    case 0: plotTag="counts"; titleTag="counts"; break;
    default:
      std::cout << "plot_case=" << plot_case << " is not ready\n";
      return retCodeError;
    }
    vector<TH2D*> *histosLeadAsH2D=NULL;
    vector<TH2D*> *histosTailAsH2D=NULL;
    switch(plot_case) {
    case 0:
      histosLeadAsH2D=&yieldsSameEventLead;
      histosTailAsH2D=&yieldsSameEventTail;
      ymax=600;
      break;
    default:
      std::cout << "plot_case=" << plot_case << " is not ready (2)\n";
      return retCodeError;
    }
    //ymax=1e6;

    std::vector<std::vector<TH1D*>*> histosLead;
    std::vector<std::vector<TH1D*>*> histosTail;
    if (!createMassProfileVec(*histosLeadAsH2D,histosLead,inpMgr.sampleNames())
	||
	!createMassProfileVec(*histosTailAsH2D,histosTail,inpMgr.sampleNames())
	)
      return retCodeError;

    TString cpName=TString("cp_") + plotTag;
    TString cpTitle=" "; //titleBase;
    cpTitle.Append(titleTag);
    ComparisonPlot_t cp(ComparisonPlot_t::_ratioPlain,cpName,cpTitle,
			"#it{M}_{ee} [GeV]","multi #it{ee} counts",
			"ratio to data yield");
    //cp.SetRefIdx(-111); // no comparison
    cp.SetYTitleSize(0.06,1.2);
    cp.SetLogx(1);
    cp.SetLogy(0);
    //cp.SetLogx(0); cp.SetXRange(59.99,120.01);
    cp.SetYRange(ymin,ymax);

    cp.SetLogy(1); cp.SetYRange(1e-3,1e3);

    TString hName=Form("hSumEE_%d",plot_case);
    TH1D *h=(TH1D*)(*histosLead[0])[0]->Clone(hName);
    h->Reset();
    h->SetDirectory(0);
    h->SetTitle(hName);

    cp.AddHist1D(hData,"data","LP skip",kWhite,1,1,-1);
    for (unsigned int ih=0; ih<histosLead[0]->size(); ++ih) {
      TH1D* currHlead=(*histosLead[0])[ih];
      TH1D* currHtail=(*histosTail[0])[ih];
      TString strLead=inpMgr.sampleName(ih) + TString(" lead");
      TString strTail=inpMgr.sampleName(ih) + TString(" tail");
      if (ih==0) {
	cp.AddHist1D(currHlead,strLead,"LP",kBlack,1,0,1);
	cp.AddHist1D(currHtail,strTail,"LP",TAttMarker(kBlack,20,1.),1,0,1);
	//printHisto(currH);
      }
      else {
	h->Add(currHlead,1.);
	h->Add(currHtail,1.);
	int colorLead=inpMgr.sampleInfo(ih)->colors[0];
	int colorTail=kYellow+(ih%2);
	cp.AddToStack(currHlead,strLead,colorLead);
	cp.AddToStack(currHtail,strTail,colorTail);
      }
    }
    h->SetMarkerStyle(20);
    cp.AddHist1D(h,"all samples","LP skip",kBlue,1,1,1);

    TString canvName=Form("cx%d",plot_case);
    TString canvTitle=Form("canv_%d",plot_case);
    TCanvas *cx=new TCanvas(canvName,canvTitle,800,700);
    cp.Prepare2Pads(cx);
    SetSideSpaces(cx,0.0,0.2,0.,0.01);
    cp.Draw(cx,false,"",1,2);
    cp.ChangeLegendPos(0.17,0.,0.1,0.);
    cx->Update();
    
    TString figName=fileBase + plotTag;
    //SaveCanvas(cx,figName);

  }

    return retCodeOk;
}

