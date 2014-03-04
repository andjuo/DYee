#include "CovariantEff.h"
//#include <TIter.h>
#include <TKey.h>
#include <TFile.h>
#include <TString.h>

// -------------------------------------------------

TString formNewTitle(TString inpFileName, TString canvName, TString &shortName);

// -------------------------------------------------

void extractFitPlots(TString fname, TString destDir) {
  TString confFileName="../config_files/data_vilnius8TeV_regSSD.conf.py";
  confFileName="../config_files/data_vilnius8TeV_regSSD_etaMax25.conf.py";

  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST;
  systMode=DYTools::UNREGRESSED_ENERGY;
  CovariantEffMgr_t mgr;
  int nExps=0;

  HERE("calling setup");
  assert(mgr.Setup(confFileName,nExps,systMode));
  assert(mgr.initOk());

  std::cout << "\n\nok. Start studies\n";

  TFile fin(fname,"read");
  if (!fin.IsOpen()) {
    std::cout << "failed to open a file <" << fin.GetName() << ">\n";
    return;
  }

  TIter nextkey(fin.GetListOfKeys());
  TKey *key;
  while ( (key = (TKey*)nextkey()) ) {
    printf("key: %s points to an object of class: %s at %lld",
	   key->GetName(),
	   key->GetClassName(),key->GetSeekKey());
    std::cout << "\n";
    if (key->GetClassName()==TString("TCanvas")) {
      TString keyName=key->GetName();
      if ((keyName.Index("canvMassDistr")==-1) && 
	  (keyName.Index("canvDistr")==-1)) {
	
	// modify the key name
	TString shortName;
	TString newTitle=formNewTitle(fname,keyName,shortName);
	if (newTitle.Length()==0) return;

	TCanvas *cx=(TCanvas*)key->ReadObj();
	cx->SetLeftMargin(0.18);
	cx->SetTitle("cx");
	//cx->Draw();
	TObject *obj=NULL;
	TIter canvObjs(cx->GetListOfPrimitives());
	while ((obj=canvObjs())) {
	  //if (obj->InheritsFrom("TNamed")) {
	  //  TNamed *named=(TNamed*)obj;
	  //  std::cout << "the name=" << named->GetTitle() << "\n";
	  //}
	  /*
	  if (obj->InheritsFrom("TPaveText")) {
	    TPaveText* pavetxt=(TPaveText*)obj;
	    std::cout << "pavetxt->text="; pavetxt->Print(); //;txt->GetText() << "\n";
	    TText* txt=pavetxt->GetLine(0);
	    txt->SetText(txt->GetX(),txt->GetY(),newTitle);
	    //txt->AddText(newTitle);
	    std::cout << "modified txt->text="; pavetxt->Print();
	  }
	  */
	  if (obj->InheritsFrom("TH1D")) {
	    TH1D *h=(TH1D*)obj;
	    std::cout << "h->GetTitle=" << h->GetTitle() << std::endl;
	    TString hTitle=h->GetTitle();
	    if (hTitle.Index("RooPlot")!=-1) {
	      h->SetTitle(newTitle);
	      h->GetXaxis()->SetTitle("m_{ee} [GeV]");
	    }
	    else {
	      std::cout << "skipped\n";
	    }
	  }
	}
	cx->Draw();
	cx->Modified();
	cx->Update();
	TString saveName=TString("fig_") + shortName;
	SaveCanvas(cx,saveName,destDir);
	//break;
      }
    }
  }
  fin.Close();

}

      // -------------------------------------

TString formNewTitle(TString inpFName, TString keyName, TString &shortName) {
  TString newTitle;
  shortName.Clear();

  DYTools::TEfficiencyKind_t kind=DYTools::effNONE;
  int idx=-1;

  // determine efficiency
  idx=keyName.Index("-id-");
  if (idx!=-1) kind=DYTools::ID;
  else { 
    idx=keyName.Index("-reco-"); 
    if (idx!=-1) kind=DYTools::RECO;
    else {
      idx=keyName.Index("-hltleg1-");
      if (idx!=-1) kind=DYTools::HLT_leg1;
      else {
	idx=keyName.Index("-hltleg2-");
	if (idx!=-1) kind=DYTools::HLT_leg2;
	else {
	  idx=keyName.Index("-hlt-");
	  if (idx!=-1) kind=DYTools::HLT;
	  else {
	    std::cout << "failed to determine efficiency kind from keyName=<" << keyName << ">\n";
	    return newTitle;
	  }
	}
      }
    }
  }

  // determine et,eta ranges
  int iEt=-1,iEtaMin=-1,iEtaMax=-1;

  for (int i=0; (iEt==-1) && (i<etBinCount); ++i) {
    TString etR=Form("et_%2.1lf_et_%2.1lf",etBinLimits[i],etBinLimits[i+1]);
    etR.ReplaceAll(".","_");
    if (keyName.Index(etR)!=-1) {
      std::cout << "et range is " << etR << "\n";
      iEt=i;
    }
  }
  for (int i=0; (iEtaMin==-1) && (i<etaBinCount); ++i) {
    TString etaR=Form("abs_eta_%4.3lf_abs_eta_%4.3lf",etaBinLimits[i],etaBinLimits[i+1]);
    etaR.ReplaceAll(".","_");
    if (keyName.Index(etaR)!=-1) {
      std::cout << "eta range is " << etaR << "\n";
      iEtaMin=i;
    }
  }
  if (iEtaMin>=0) { 
    iEtaMax=iEtaMin+1;
  }
  else {
    // maybe it was a merged bin
    for (int i=0; (iEtaMin==-1) && (i<etaBinCount-1); ++i) {
      TString etaR=Form("abs_eta_%4.3lf_abs_eta_%4.3lf",etaBinLimits[i],etaBinLimits[i+2]);
      etaR.ReplaceAll(".","_");
      if (keyName.Index(etaR)!=-1) {
	std::cout << "eta range is " << etaR << "\n";
	iEtaMin=i; iEtaMax=i+2;
      }
    }
  }

  int pass=-1;
  if (keyName.Index("-pass")!=-1) pass=1;
  else if (keyName.Index("-fail")!=-1) pass=0;

  if (inpFName.Index("_data")!=-1) {
    newTitle=TString("Data, ");
    shortName="_data";
  }
  else if (inpFName.Index("_mc")!=-1) {
    newTitle=TString("MC, ");
    shortName="_mc";
  }
  else {
    newTitle=TString("data??? ");
    shortName="_xx";
  }

  newTitle.Append(EfficiencyKindName(kind));
  newTitle.Append(Form(", %2.0lf<#it{E}_{T}<%2.0lf GeV, %4.2lf<|#it{#eta}|<%4.2lf",etBinLimits[iEt],etBinLimits[iEt+1],etaBinLimits[iEtaMin],etaBinLimits[iEtaMax]));
  if (pass==1) newTitle.Append(", pass");
  else if (pass==0) newTitle.Append(", fail");
  else newTitle.Append(", fail?");

  shortName.Append(EfficiencyKindName(kind));
  shortName.Append(Form("_Et_%2.0lf_%2.0lf_eta_%4.2lf_%4.2lf",etBinLimits[iEt],etBinLimits[iEt+1],etaBinLimits[iEtaMin],etaBinLimits[iEtaMax]));
  shortName.ReplaceAll(".","_");
  if (pass==1) shortName.Append("_pass");
  else if (pass==0) shortName.Append("_fail");
  else shortName.Append("_xx");

  return newTitle;
}
