#include "derivePreFsrCS.C"

int run_derivePreFsrCS(int debug, int runCase=3,
		       double user_xmin=0, double user_xmax=-1,
		       double user_ymin=0, double user_ymax=-1,
		       int user_logX=-1, int user_logY=-1)
{
  DYTools::TRunMode_t runMode=DebugInt2RunMode(debug);
  DYTools::TSystematicsStudy_t systMode=DYTools::NO_SYST; 
  TString fname1="../config_files/data_vilnius8TeV_reg_ntuple.conf.py";
  TString fname2="../config_files/data_vilnius8TeV_reg_skim.conf.py";
  
  int runNtuple=((runCase&1)!=0) ? 1:0;
  int runSkim  =((runCase&2)!=0) ? 1:0;

  RangeUser_t rUser(user_xmin,user_xmax,user_ymin,user_ymax,user_logX,user_logY);

  int res=retCodeOk;
  if ((res==retCodeOk) && runNtuple) res=derivePreFsrCS(fname1,runMode,systMode,"ntuple",&rUser);
  if ((res==retCodeOk) && runSkim  ) res=derivePreFsrCS(fname2,runMode,systMode,"skim",&rUser);
  if (res!=retCodeOk) HERE("error in run_derivePreFsrCS\n");
  return res;
}
