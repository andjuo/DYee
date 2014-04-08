{  

  // 
  // Uncomment Include/rootlogon.C, if ComparisonPlot.hh is included in studyCovGen.C
  //
  gROOT->ProcessLine(".x ../Include/rootlogon.C"); 

  // 
  //  These lines needed for DYee CovariantEff module
  //
  gROOT->ProcessLine(".L ../EventScaleFactors/RooCMSShape.cc+");
  gROOT->ProcessLine(".L ../EventScaleFactors/cutFunctions.cc+");
  gROOT->ProcessLine(".L ../EventScaleFactors/fitFunctionsCore.cc+"); // not necessary
  gROOT->ProcessLine(".L ../EventScaleFactors/fitFunctions.cc+"); // not necessary
  gROOT->ProcessLine(".L ../EventScaleFactors/tnpSelectEvents.hh+");
  gROOT->ProcessLine(".L ../EventScaleFactors/calcEventEffLink.cc+");
 
  gROOT->ProcessLine(".L CovariantEff.cc+");


  //gROOT->ProcessLine(".L colorPalettes.hh");
}
