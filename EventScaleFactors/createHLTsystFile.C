#include <TROOT.h>
#include <TFile.h>
#include <TMatrixD.h>
#include <TString.h>
#include <iostream>


void createHLTsystFile() {
  TMatrixD base(6,5);
  base.Zero();

  TMatrixD dataLeg1(base);
  TMatrixD dataLeg2(base);
  TMatrixD mcLeg1(base);
  TMatrixD mcLeg2(base);

  // using non-regressed electron energy
  // numbers from 2014 Feb 26

  // data Leg1
  dataLeg1(0,0)=0.118;
  dataLeg1(1,0)=0.012;

  dataLeg1(0,1)=0.022;
  dataLeg1(1,1)=0.010;

  dataLeg1(0,2)=0.110;
  dataLeg1(1,2)=0.010;
  dataLeg1(2,2)=0.002;
  
  dataLeg1(0,3)=0.085;
  dataLeg1(1,3)=0.010;

  dataLeg1(0,4)=0.680;
  dataLeg1(1,4)=0.010;

  // data Leg2
  dataLeg2(0,0)=0.003;
  dataLeg2(1,0)=0.002;

  dataLeg2(0,1)=0.004;
  dataLeg2(1,1)=0.003;

  dataLeg2(0,2)=0.012;
  dataLeg2(1,2)=0.002;
  dataLeg2(2,2)=0.001;

  dataLeg2(0,3)=0.001;
  dataLeg2(1,3)=0.001;

  dataLeg2(0,4)=0.0065;
  dataLeg2(1,4)=0.0012;

  // mc Leg1
  mcLeg1(0,0)=0.33;
  mcLeg1(1,0)=0.015;
  mcLeg1(2,0)=0.0025;

  mcLeg1(0,1)=0.076;
  mcLeg1(1,1)=0.0046;
  mcLeg1(2,1)=0.001;

  mcLeg1(0,2)=0.132;
  mcLeg1(1,2)=0.0042;
  mcLeg1(2,2)=0.001;

  mcLeg1(0,3)=0.171;
  mcLeg1(1,3)=0.002;

  mcLeg1(0,4)=0.196;
  mcLeg1(1,4)=0.013;

  // mc Leg2
  mcLeg2(0,0)=0.0135;
  mcLeg2(1,0)=0.007;
  mcLeg2(2,0)=0.002;

  mcLeg2(0,1)=0.012;
  mcLeg2(1,1)=0.0045;
  mcLeg2(2,1)=0.001;

  mcLeg2(0,2)=0.0064;
  mcLeg2(1,2)=0.0025;
  mcLeg2(2,2)=0.0015;

  mcLeg2(0,3)=0.003;
  mcLeg2(1,3)=0.0025;
  mcLeg2(2,3)=0.001;
  
  mcLeg2(0,4)=0.003;
  mcLeg2(1,4)=0.0025;


  std::cout << "\ndata HLTleg1: "; dataLeg1.Print();
  std::cout << "\ndata HLTleg2: "; dataLeg2.Print();
  std::cout << "\nmc HLTleg1: "; mcLeg1.Print();
  std::cout << "\nmc HLTleg2: "; mcLeg2.Print();


  TString foutName="unregHLTSystematics20140226.root";
  TFile f(foutName,"recreate");
  dataLeg1.Write("hltLeg1_data_syst_rel_error");
  dataLeg2.Write("hltLeg2_data_syst_rel_error");
  mcLeg1.Write("hltLeg1_mc_syst_rel_error");
  mcLeg2.Write("hltLeg2_mc_syst_rel_error");
  f.Close();
  std::cout << "file <" << f.GetName() << "> created\n";
  return;
}
