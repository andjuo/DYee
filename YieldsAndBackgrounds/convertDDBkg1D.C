#include "../Include/DYTools.hh"
#include "../Include/InputFileMgr.hh"

int convertDDBkg1D() {
  if (!DYTools::setup(0)) {
    return retCodeError;
  }
  
  TString conf="../config_files/data_vilnius8TeV_regSSD.conf.py";
  InputFileMgr_t inpMgr;
  if (!inpMgr.Load(conf)) return retCodeError;

  for (int ifile=0; ifile<2; ++ifile) {
    TString fname=(ifile==0) ? inpMgr.GetTrue2eDDBkgFName() : inpMgr.GetFakeDDBkgFName();
    TString fieldBase=(ifile==0) ? "true2eBackgroundFromData" : "fakeBackgroundFromData";
    Ssiz_t idx=fname.Last('/');
    TString outFname=fname(idx-1,fname.Length());
    outFname(0,1)='.';
    std::cout << "outFname=<" << outFname << ">\n";
    TFile fout(outFname,"recreate");
    
    TMatrixD Mnew(43,1);
    Mnew.Zero();

    for (int ifield=0; ifield<3; ++ifield) {
      TString field=fieldBase;
      if (ifield==1) field.Append("Error");
      else if (ifield==2) field.Append("ErrorSyst");
      TMatrixD *M=loadMatrix(fname,field,41,1, 1);
      for (int ir=0; ir<41; ++ir) Mnew(ir+1,0)=(*M)(ir,0);
      delete M;
      fout.cd();
      Mnew.Write(field);
    }
    writeBinningArrays(fout);
    fout.Close();
  }
  return retCodeOk;
}

