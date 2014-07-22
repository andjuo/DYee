{

  gROOT->ProcessLine(".L ../Include/DYTools.cc+");
  gROOT->ProcessLine(".L ../Include/DYTools.hh+");

  // Load "MIT Style" plotting
  gROOT->Macro("../Include/MitStyleRemix.cc+");
  gROOT->Macro("../Include/CPlot.cc+");

  // Load structures for ntuple analysis
  // The ones below have class definitions and need to be compiled
  gROOT->ProcessLine(".L ../Include/TElectron.hh+");
  gROOT->ProcessLine(".L ../Include/TDielectron.hh+");
  gROOT->ProcessLine(".L ../Include/TMuon.hh+");
  gROOT->ProcessLine(".L ../Include/TEventInfo.hh+");
  gROOT->ProcessLine(".L ../Include/TGenInfo.hh+");
  gROOT->ProcessLine(".L ../Include/TPhoton.hh+");
  gROOT->ProcessLine(".L ../Include/TGenPhoton.hh+");
  gROOT->ProcessLine(".L ../Include/TJet.hh+");
  gROOT->ProcessLine(".L ../Include/TVertex.hh+");
  gROOT->ProcessLine(".L ../Include/EleIDCuts.hh+");
  gROOT->ProcessLine(".L ../Include/JsonParser.cc+");

  gROOT->ProcessLine(".L ../Include/TriggerSelection.hh+");
  gROOT->ProcessLine(".L ../Include/FEWZ.cc+");
  gROOT->ProcessLine(".L ../Include/AccessOrigNtuples.cc+");

  //
  // ComparisonPlot.hh
  gROOT->ProcessLine(".L ../Include/MyTools.cc+");
  gROOT->ProcessLine(".L ../Include/PUReweight.cc+"); // depends on MyTools.hh
  gROOT->ProcessLine(".L ../Include/EventWeight.cc+");
  gROOT->ProcessLine(".L ../Include/ZeeData.hh+"); // depends on EventWeight
  gROOT->ProcessLine(".L ../Include/InputFileMgr.cc+");

  gROOT->ProcessLine(".L ../Include/ElectronEnergyScale.cc+");
  gROOT->ProcessLine(".L ../Include/EventSelector.cc+");

}
