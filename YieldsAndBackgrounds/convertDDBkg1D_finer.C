#include "../Include/DYTools.hh"
#include "../Include/InputFileMgr.hh"

// ----------------------------------------------------------

const int _nMassBins2012 = 41;
const double _massBinLimits2012[_nMassBins2012+1] =
  {15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 64, 68, 72, 76,
   81, 86, 91, 96, 101, 106, 110, 115, 120, 126, 133, 141,
   150, 160, 171, 185, 200, 220, 243, 273, 320, 380, 440,
   510, 600, 1000, 1500, 2000 }; // 41 bin

// ----------------------------------------------------------
// ----------------------------------------------------------

int convertDDBkg1D_finer() {
  if (!DYTools::setup(0)) {
    std::cout << "failed to initialize\n";
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
    
    TMatrixD Mnew(DYTools::nMassBins,1);
    Mnew.Zero();

    for (int ifield=0; ifield<3; ++ifield) {
      TString field=fieldBase;
      if (ifield==1) field.Append("Error");
      else if (ifield==2) field.Append("ErrorSyst");
      TMatrixD *M=loadMatrix(fname,field,41,1, 1);
      for (int ir=0; ir<41; ++ir) {
	std::cout << "using value from M=" << _massBinLimits2012[ir] << " .. "
		  << _massBinLimits2012[ir+1] << " for M="
		  << DYTools::massBinLimits[2*ir] << " .. "
		  << DYTools::massBinLimits[2*ir+1] << " .. "
		  << DYTools::massBinLimits[2*ir+2] << "\n";
	if (ir!=40) {
	  if ((_massBinLimits2012[ir]!=DYTools::massBinLimits[2*ir]) ||
	      (_massBinLimits2012[ir+1]!=DYTools::massBinLimits[2*ir+2])) {
	    std::cout << "range error\n";
	    return retCodeError;
	  }
	}
	double halfValue=0.5*((ifield) ? sqrt(2.) : 1.) * (*M)(ir,0);
	std::cout << "val=" << (*M)(ir,0) << ", halfVal=" << halfValue << "\n";
	if (ir<40) {
	  Mnew(2*ir  ,0)=halfValue;
	  Mnew(2*ir+1,0)=halfValue;
	}
	else Mnew(2*ir,0) = (*M)(ir,0);

	if ((ifile==1) && (ifield==0)) {
	  double val=(*M)(ir,0);
	  double set_last_val=-1;
	  if (ir==39) { set_last_val=2; }
	  if (set_last_val>=0.) {
	    Mnew(2*ir+1, 0)= set_last_val;
	    Mnew(2*ir  , 0)= val - set_last_val;
	  }
	}
      }
      delete M;
      fout.cd();
      Mnew.Print();
      Mnew.Write(field);
    }
    writeBinningArrays(fout,"convertDDBkg1D_finer");
    fout.Close();
  }
  return retCodeOk;
}

