{
  gROOT->ProcessLine(".L compareHistos2.C+");
  compareHistos2("../config_files/data_vilnius7TeV_test.conf.py","test_new.root","test_new_2D.root");
}
