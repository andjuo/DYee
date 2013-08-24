{
  gROOT->ProcessLine(".L compareHistos.C+");
  compareHistos("../config_files/data_vilnius7TeV_test.conf.py","../../DrellYanDMDY-20130131/YieldsAndBackgrounds/test.root","test_new.root");
}
