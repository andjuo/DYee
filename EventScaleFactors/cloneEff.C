#include "../Include/DYTools.hh"
#include "calcEventEffLink.h"

void cloneEff() {
  

  for (int iBr=0; iBr<6; ++iBr) {
    TString fname="dir-Rami/efficiency_TnP_1D_Full2012_dataRECO_fit-fitEtBins6EtaBins9_PU.root";
    TString outFName="dir-RamiClone/efficiency_TnP_1D_Full2012_dataRECO_fit-fitEtBins7EtaBins9_PU.root";

    if (iBr==0) {
    }
    else if (iBr==1) {
      fname.ReplaceAll("RECO","ID");
      outFName.ReplaceAll("RECO","ID");
    }
    else if (iBr==2) {
      fname.ReplaceAll("RECO_fit-fit","HLT_count-count");
      outFName.ReplaceAll("RECO_fit-fit","HLT_count-count");
    }
    else if (iBr==3) {
      fname.ReplaceAll("dataRECO_fit-fit","mcRECO_count-count");
      outFName.ReplaceAll("dataRECO_fit-fit","mcRECO_count-count");
    }
    else if (iBr==4) {
      fname.ReplaceAll("dataRECO_fit-fit","mcID_count-count");
      outFName.ReplaceAll("dataRECO_fit-fit","mcID_count-count");
    }
    else if (iBr==5) {
      fname.ReplaceAll("dataRECO_fit-fit","mcHLT_count-count");
      outFName.ReplaceAll("dataRECO_fit-fit","mcHLT_count-count");
    }
    else if (iBr==6) {
      fname.ReplaceAll("dataRECO_fit-fit","dataHLTleg1_count-count");
      outFName.ReplaceAll("dataRECO_fit-fit","dataHLTleg1_count-count");
    }
    else if (iBr==7) {
      fname.ReplaceAll("dataRECO_fit-fit","dataHLTleg2_count-count");
      outFName.ReplaceAll("dataRECO_fit-fit","dataHLTleg2_count-count");
    }
    else if (iBr==8) {
      fname.ReplaceAll("dataRECO_fit-fit","mcHLTleg1_count-count");
      outFName.ReplaceAll("dataRECO_fit-fit","mcHLTleg1_count-count");
    }
    else if (iBr==9) {
      fname.ReplaceAll("dataRECO_fit-fit","mcHLTleg2_count-count");
      outFName.ReplaceAll("dataRECO_fit-fit","mcHLTleg2_count-count");
    }
    else {
      std::cout << "iBr error\n";
      return;
    }

    TMatrixD *eff, *effLo, *effHi;
    int weighted=(fname.Index("count-count")>0) ? 1:0;
    if (!loadEff(fname,weighted,&eff,&effLo,&effHi)) {
      return ;
    }
    TMatrixD effNew(eff->GetNrows()+1,eff->GetNcols());
    TMatrixD effLoNew(effNew);
    TMatrixD effHiNew(effNew);

    for (int ir=0; ir<eff->GetNrows(); ++ir) {
      for (int ic=0; ic<eff->GetNcols(); ++ic) {
	effNew(ir,ic) = (*eff)(ir,ic);
	effLoNew(ir,ic) = (*effLo)(ir,ic);
	effHiNew(ir,ic) = (*effHi)(ir,ic);
      }
    }
    // clone the last row
    int ir=eff->GetNrows();
    for (int ic=0; ic<eff->GetNcols(); ++ic) {
      effNew(ir,ic) = (*eff)(ir-1,ic);
      effLoNew(ir,ic) = (*effLo)(ir-1,ic);
      effHiNew(ir,ic) = (*effHi)(ir-1,ic);
    }

    TString field="effArray2D";
    TString fieldErrLo="effArrayErrLow2D";
    TString fieldErrHi="effArrayErrHigh2D";
    if (weighted) {
      field.Append("Weighted");
      fieldErrLo.Append("Weighted");
      fieldErrHi.Append("Weighted");
    }

    TFile outFile(outFName,"recreate");
    effNew.Write(field);
    effLoNew.Write(fieldErrLo);
    effHiNew.Write(fieldErrHi);
    outFile.Close();

    delete eff;
    delete effLo;
    delete effHi;
  }

  return;
}
