#include "../Include/DYTools.hh"
#include "../Include/MyTools.hh"
#include "../Include/ComparisonPlot.hh"
#include <RooRealVar.h>
//#include <RooFormulaVar.h>
#include <RooBreitWigner.h>
#include <RooCBShape.h>
#include <RooFFTConvPdf.h>
#include <RooPlot.h>
#include <RooDataHist.h>
#include <RooFitResult.h>

// sector: -1 (all), 0 (barrel-barrel), 1 (endcap-endcap), 2 (barrel-endcap)
int fitZpeak(int sector=-1,
	     DYTools::TSystematicsStudy_t systMode=DYTools::APPLY_ESCALE) {

  if (!DYTools::setup(0)) {
    std::cout << "failed to initialize\n";
    return retCodeError;
  }

  TString subDirName="mass_Zpeak_1GeV";
  TString histoNameBase="hZpeak";
  TString sectorTag;
  switch(sector) {
  case -1: break;
  case 0: subDirName="mass_ZpeakBB_1GeV"; sectorTag.Append("BB"); break;
  case 1: subDirName="mass_ZpeakEE_1GeV"; sectorTag.Append("EE"); break;
  case 2: subDirName="mass_ZpeakBE_1GeV"; sectorTag.Append("BE"); break;
  default:
    std::cout << "sector=" << sector << " not ready\n";
    return retCodeError;
  }
  histoNameBase.Append(sectorTag);
  histoNameBase.Append("_");

  TString fname="../root_files_reg/yield/DY_j22_19712pb_ApplyEscale/all_yield_1D_mbZpeak__peak20140220.root";
  if (systMode==DYTools::APPLY_ESCALE) {
  }
  else if (systMode==DYTools::ESCALE_DIFF_0000) {
    fname="../root_files_reg/yield/DY_j22_19712pb_EScaleDiff0000/all_yield_1D_mbZpeak.root";
  }
  else {
    std::cout << "Not ready for this systMode\n";
    return retCodeError;
  }

  TFile fin(fname,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open the file <" << fin.GetName() << ">\n";
    return retCodeError;
  }
  TH1D *hData=new TH1D(histoNameBase + TString("data"),"data",1,0.,1.);
  TH1D *hZee= new TH1D(histoNameBase + TString("zee"),"Zee",1,0.,1.);
  if (!loadHisto(fin,&hData,subDirName) ||
      !loadHisto(fin,&hZee, subDirName)) {
    std::cout << "failed to load histos\n";
    fin.Close();
    return retCodeError;
  }
  fin.Close();

  TH1D *hZee_chk=NULL;
  if (0) {
    hZee->SetName("hZee_orig");
    TFile fChk("../root_files_reg/yield/DY_j22_19712pb_ApplyEscale/all_yield_1D_mbZpeak__peak20140220.root","read");
    hZee_chk= new TH1D(histoNameBase + TString("zee"),"Zee escale",1,0.,1.);
    if (!loadHisto(fChk,&hZee_chk,subDirName)) return retCodeError;
    fChk.Close();
  }

  if (1) { // scaling distorts the fit?!
    if (!hZee_chk) {
      std::cout << "hZee->Integral=" << hZee->Integral() << "\n";
      std::cout << "hData->Integral=" << hData->Integral() << "\n";
      hZee->Scale(hData->Integral()/hZee->Integral());
    }
  }


  printHisto(hData);
  printHisto(hZee);
  if (hZee_chk) printHisto(hZee_chk);

  std::cout << "data: " << hData->GetMean() << " +- " << hData->GetRMS() << "\n";
  std::cout << "mc:   " << hZee->GetMean() << " +- " << hZee->GetRMS() << "\n";
  std::cout << "diff: " << (hData->GetMean() - hZee->GetMean()) << "\n";
  if (hZee_chk) {
    std::cout << "mc(chk):   " << hZee_chk->GetMean() << " +- " << hZee_chk->GetRMS() << "\n";
  }

  //return retCodeStop;

  /*** Observables ***/
  RooRealVar mass("mass","#it{M}_{ee}",60.,120.,"GeV");
  RooDataHist dataHist_Data("data","measured data",mass,hData);
  RooDataHist dataHist_MC("mc","simulation",mass,hZee);

  /*** Fit Parameter ***/
  RooRealVar bw_Zmass("bw_Zmass","bw_Zmass",91.1876,87,95);
  RooRealVar bw_Zwidth("bw_Zwidth","bw_Zwidth",2.4952,1,4);

  RooRealVar cb_offset_data("cb_offset_data","cb_offset_data",-10,10);
  RooRealVar cb_sigma_data("cb_sigma_data","cb_sigma_data",0,10);
  RooRealVar cb_offset_mc("cb_offset_mc","cb_offset_mc",-10,10);
  RooRealVar cb_sigma_mc("cb_sigma_mc","cb_sigma_mc",0,10);
  RooRealVar cb_cut("cb_cut","#alpha_mc",0,10);
  RooRealVar cb_power("cb_power","N_mc",0,10);

  //RooRealVar beta("beta","beta",0,-1,1);
//RooFormulaVar massCorr("massCorr","massCorr","@0/(1+@1/2)",RooArgList(mass,beta));

  bw_Zmass.setConstant();
  bw_Zwidth.setConstant();

  /*** Breit Wigner ***/
  RooBreitWigner bw("bw","Breit Wigner",mass,bw_Zmass,bw_Zwidth);

  /*** Crystal Ball for MC ***/
  RooCBShape cb_mc("cb_mc","Crystal Ball (MC)",mass,cb_offset_mc,cb_sigma_mc,
		   cb_cut,cb_power);

  /*** Convolved Model for MC ***/
  RooFFTConvPdf bw_cb_conv_mc("bw_cb_conv_mc","{BWxCB} MC",mass,bw,cb_mc);

  //   Fit MC distribution
  RooFitResult *fitRes_mc=
    bw_cb_conv_mc.fitTo(dataHist_MC,RooFit::Extended(false),
			RooFit::SumW2Error(kTRUE),
			Save(true));

  /*** Crystal Ball for data ***/
  cb_cut.setConstant();
  cb_power.setConstant();
  RooCBShape cb_data("cb_data","Crystal Ball (data)",mass,cb_offset_data,cb_sigma_data,
		   cb_cut,cb_power);

  /*** Convolved Model for data ***/
  RooFFTConvPdf bw_cb_conv_data("bw_cb_conv_data","{BWxCB} data",mass,bw,cb_data);


  //   Fit data distribution
  RooFitResult *fitRes_data=
    bw_cb_conv_data.fitTo(dataHist_Data,RooFit::Extended(false),
			  RooFit::SumW2Error(kTRUE),
			  Save(true));

  //cb_cut.setConstant(false);
  //cb_power.setConstant(false);

  TString explain= sectorTag;
  explain.Append((systMode==DYTools::APPLY_ESCALE) ?
		 " Peak-corrected" : " No correction");
  RooPlot* frame= mass.frame(RooFit::Title(TString("MC.") + explain));
  dataHist_MC.plotOn(frame,RooFit::MarkerColor(kGreen+1)); // Binning(100);
  //dataHist_Data.plotOn(frame,RooFit::MarkerColor(kRed+1)); // Binning(100);
  if (hZee_chk) {
    RooDataHist dataHist_MCchk("mcChk","simulation",mass,hZee_chk);
    dataHist_MCchk.plotOn(frame,RooFit::MarkerColor(kRed+1));
  }
  bw_cb_conv_mc.plotOn(frame,RooFit::LineColor(kBlue),RooFit::LineStyle(2));

  RooPlot* frame_data= mass.frame(RooFit::Title(TString("Data.") + explain));
  dataHist_Data.plotOn(frame_data,RooFit::MarkerColor(kBlack)); // Binning(100);
  bw_cb_conv_data.plotOn(frame_data,RooFit::LineColor(kBlue),RooFit::LineStyle(2));
  //bw_cb_conv_data.paramOn(frame_data, RooFit::Parameters( RooArgSet(cb_offset_data,cb_sigma_data) ), RooFit::Format("NE"), RooFit::Layout(0.25,0.89,0.89));

  //RooPlot* frameMCPar= mass.frame(RooFit::Title("MC fit parameters"));
  bw_cb_conv_mc.paramOn(frame, RooFit::Parameters( RooArgSet(cb_offset_mc,cb_sigma_mc,cb_cut,cb_power) ), RooFit::Format("NE"), RooFit::Layout(0.3,0.9,0.93));
  //RooPlot* frameDataPar= mass.frame(RooFit::Title("Data fit parameters"));
  bw_cb_conv_data.paramOn(frame_data, RooFit::Parameters( RooArgSet(cb_offset_data,cb_sigma_data) ), RooFit::Format("NEU"), RooFit::Layout(0.3,0.89,0.93));

  TCanvas *cx=new TCanvas("cx","cx",1200,600);
  cx->Divide(2,1);
  cx->cd(1);
  frame->Draw();
  //frame->getAttText()->SetTextSize(0.01);
  cx->cd(2);
  frame_data->Draw();
  //cx->Update();

 return retCodeOk;
}
