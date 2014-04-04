{
  gROOT->Macro("../DYTools.cc+");
  gROOT->Macro("../DYTools.hh+");

  // to run chkExtractHisto
  gROOT->ProcessLine(".L ../CPlot.cc+");
  gROOT->ProcessLine(".L ../MyTools.cc+");
}
