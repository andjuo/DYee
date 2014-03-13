// 2014.03.13 Macro to add binningArrays to a file

#include "../Include/DYTools.hh"
#include "../Include/InputFileMgr.hh"

int convertDDBkg2D() {
  if (DYTools::study2D!=1) {
    std::cout << "this macro works on 2D case\n";
    return retCodeError;
  }
  
  for (int ifile=0; ifile<2; ++ifile) {
    if (ifile==1) continue;
    TString fname=(ifile==0) ? "./true2eBkgDataPoints_2D.root" : "./fakeBkgDataPoints_2D.root";
    TString fieldBase=(ifile==0) ? "true2eBackgroundFromData" : "fakeBackgroundFromData";
    Ssiz_t idx=fname.Last('/');
    TString outFname=fname(idx-1,fname.Length());
    outFname(0,1)='.';
    outFname.ReplaceAll("_2D","_20140312_2D");
    std::cout << "outFname=<" << outFname << ">\n";
    TFile fout(outFname,"recreate");
    
    for (int ifield=0; ifield<3; ++ifield) {
      TString field=fieldBase;
      if (ifield==1) field.Append("Error");
      else if (ifield==2) field.Append("ErrorSyst");
      TMatrixD *M=loadMatrix(fname,field,DYTools::nMassBins,DYTools::nYBinsMax, 1);
      fout.cd();
      M->Write(field);
      delete M;
    }
    writeBinningArrays(fout);
    fout.Close();
  }
  return retCodeOk;
}

