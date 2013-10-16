{  
  
  gROOT->ProcessLine(".x ../Include/rootlogon.C");

  gROOT->ProcessLine(".L RooCMSShape.cc+");
  gROOT->ProcessLine(".L cutFunctions.cc+");
  gROOT->ProcessLine(".L fitFunctionsCore.cc+");
  gROOT->ProcessLine(".L fitFunctions.cc+");
  gROOT->ProcessLine(".L tnpSelectEvents.hh+");

  //
  // Some compare*C macros need calcEventEffLink.cc
  //
  if (0) {
    std::cout << "\n\n\t\tlinking calcEventEffLink.cc\n\n";
    gROOT->ProcessLine(".L calcEventEffLink.cc+");
  }

}

