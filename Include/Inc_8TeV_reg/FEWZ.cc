#include "../Inc_8TeV_reg/FEWZ.hh"
#include <TFile.h>
#include <TString.h>

// ----------------------------------------------------------

#ifndef _8TeV_analysis_
// Constructor for 7 TeV analysis
FEWZ_t::FEWZ_t(bool loadWeights, bool do_cutZPT100) : fInitialized(kFALSE), fCutZPT100(do_cutZPT100) {
  if (loadWeights) {
    TFile fweights("../root_files/fewz/weights_stepwise_prec10-5_fine12.root");
    if(do_cutZPT100)
      cout << "FEWZ_t NOTE: in MC, for Z/gamma* PT>100 the FEWZ weights for 80<PT<100 GeV are used!" << endl;
    if( !fweights.IsOpen() ) assert(0);
    bool ok=kTRUE;
    for(int i=0; ok && (i<_nMassBinsFEWZ); i++){
      TString hnames = TString::Format("weight_%02d",i+1);
      weights[i] = (TH2D*)fweights.Get(hnames);
      hnames = TString::Format("h_weighterror_%02d",i+1);
      weightErrors[i] = (TH2D*)fweights.Get(hnames);
      if (!weights[i] || !weightErrors[i]) ok=kFALSE;
      else {
	weights[i]->SetDirectory(0);
	weightErrors[i]->SetDirectory(0);
      }
    }
    fInitialized=ok;
  }
}
#endif

// ----------------------------------------------------------

#ifdef _8TeV_analysis_
// Constructor for 8 TeV analysis. Everything is defined in a macro
FEWZ_t::FEWZ_t(bool loadWeights, bool do_cutZPT100)
  : fInitialized(loadWeights), fCutZPT100(do_cutZPT100) {
  std::cout << "FEWZ weights for 8TeV ver.2a (2014.04.23)\n";
  for (int i=0; i<_nMassBinsFEWZ; ++i) {
    weights[i]=NULL;
    weightErrors[i]=NULL;
  }
}
#endif

// ----------------------------------------------------------
// ----------------------------------------------------------

#ifdef _8TeV_analysis_
int FEWZ_t::prepareHisto_8TeV(Int_t mass_bin, int abs_rapidity) {

  // if histo is available, do nothing
  if (weights[mass_bin]) return 1;

  // create the histogram
  TString name=Form("FEWZ_weights_8TeV_%2.0lf_%2.0lf",_massBinLimitsFEWZ[mass_bin],_massBinLimitsFEWZ[mass_bin+1]);
  //std::cout << "preparing histo " << name << ", mass_bin=" << mass_bin << "\n";
  const int nYbins_loc=(abs_rapidity) ? nrapbins_FEWZ8TeV : (2*nrapbins_FEWZ8TeV-1);
  double Ybins[nYbins_loc];
  const double maxRapidity=10.; // artificial limits
  if (abs_rapidity) {
    for (int iY=0; iY<nrapbins_FEWZ8TeV; ++iY) Ybins[iY]=rap_bin_FEWZ8TeV[iY];
    Ybins[nrapbins_FEWZ8TeV-1]=maxRapidity;
  }
  else {
    for (int iY=0; iY<nrapbins_FEWZ8TeV-1; ++iY) {
      Ybins[iY] = -rap_bin_FEWZ8TeV[nrapbins_FEWZ8TeV - 1 - iY];
      Ybins[iY + nrapbins_FEWZ8TeV - 1]=rap_bin_FEWZ8TeV[iY];
    }
    Ybins[0] = -maxRapidity;
    Ybins[nYbins_loc-1] = maxRapidity;
    //for (int iY=0; iY<nYbins_loc; ++iY) std::cout << " " << Ybins[iY]; std::cout << std::endl;
  }

  double ptbins[nptbins_FEWZ8TeV];
  for (int ipt=0; ipt<nptbins_FEWZ8TeV; ++ipt) ptbins[ipt]=pt_bin_FEWZ8TeV[ipt];
  ptbins[nptbins_FEWZ8TeV-1]=800.; // artificial reduction
  TH2D *h2=new TH2D(name,name, nptbins_FEWZ8TeV-1,ptbins, nYbins_loc-1,Ybins);
  h2->SetDirectory(0);
  h2->Reset();
  weights[mass_bin]=h2;
  double mass=this->getMassCenter(mass_bin);
  for (int ipt=0; ipt<nptbins_FEWZ8TeV-1; ++ipt) {
    double pt=0.5*(pt_bin_FEWZ8TeV[ipt] + pt_bin_FEWZ8TeV[ipt+1]);
    for (int iY=0; iY<nYbins_loc-1; ++iY) {
      double y=0.5*(Ybins[iY] + Ybins[iY+1]);
      double w=this->getWeight(mass, pt, y);
      //std::cout << "weight for pt=" << pt << ", y=" << y << " is " << w << "\n";
      h2->SetBinContent(ipt+1,iY+1, w);
      h2->SetBinError(ipt+1,iY+1, 0.);
    }
  }
  return 1;
}
#endif

// ----------------------------------------------------------
// ----------------------------------------------------------

#ifdef _8TeV_analysis_

// A macro from Alexey Svyatkovskiy <asvyatko@purdue.edu> providing
// FEWZ NNLO/NLO correction weight maps (August 08, 2013)

inline
int Find_Index1( double muon_input, bool isPt )
{

   const int nrapbins = 6;
   double rap_bin[nrapbins] = {0.,0.7,1.1,1.9,2.4,1000.0};
   const int nptbins = 11;
   double pt_bin[nptbins] = {0.0,20.0,30.,35.,40.,45.,50.,60.,90.,200.,1000.};
   int Index = -1;

  if( isPt ) {
    for( int i = 0; i < nptbins-1; i++ ) {
      if( muon_input >= pt_bin[i] && muon_input < pt_bin[i+1] ) {
             Index = i;
             break;
      }
      else if( muon_input >= pt_bin[nptbins-1] ) Index = nptbins-2;
    }
    return Index;
  } else {
    for( int i = 0; i < nrapbins-1; i++ ) {
      if( muon_input >= rap_bin[i] && muon_input < rap_bin[i+1] ) {
          Index = i;
          break;
       }
      else if( muon_input >= rap_bin[nrapbins-1] ) Index = nrapbins-2;
    }
    return Index;
  }
}

inline
int Find_Index2( double muon_input, bool isPt )
{

   const int nrapbins = 4;
   double rap_bin[nrapbins] = {0.,0.7,1.9,1000.0};
   const int nptbins = 11;
   double pt_bin[nptbins] = {0.0,20.0,30.,35.,40.,45.,50.,60.,90.,200.,1000.};
   int Index = -1;

  if( isPt ) {
    for( int i = 0; i < nptbins-1; i++ ) {
      if( muon_input >= pt_bin[i] && muon_input < pt_bin[i+1] ) {
             Index = i;
             break;
      }
      else if( muon_input >= pt_bin[nptbins-1] ) Index = nptbins-2;
    }
    return Index;
  } else {
    for( int i = 0; i < nrapbins-1; i++ ) {
      if( muon_input >= rap_bin[i] && muon_input < rap_bin[i+1] ) {
          Index = i;
          break;
       }
      else if( muon_input >= rap_bin[nrapbins-1] ) Index = nrapbins-2; 
    }
    return Index;
  }
}

inline
int Find_Index3( double muon_input, bool isPt )
{

   const int nptbins = 4;
   double pt_bin[nptbins] = {0.0,20.,100.,1000.};
   const int nrapbins = 3;
   double rap_bin[nrapbins] = {0.,0.7,1000.0};
   int Index = -1;

  if( isPt ) {
    for( int i = 0; i < nptbins-1; i++ ) {
      if( muon_input >= pt_bin[i] && muon_input < pt_bin[i+1] ) {
             Index = i;
             break;
      }
      else if( muon_input >= pt_bin[nptbins-1] ) Index = nptbins-2;
    }
    return Index;
  } else {
    for( int i = 0; i < nrapbins-1; i++ ) {
      if( muon_input >= rap_bin[i] && muon_input < rap_bin[i+1] ) {
          Index = i;
          break;
       }
      else if( muon_input >= rap_bin[nrapbins-1] ) Index = nrapbins-2;
    }
    return Index;
  }
}


double weight_FEWZ8TeV( double muon_pt, double muon_rap, double mass) {
 int index_pt = -1;
 int index_rap = -1;

 if (mass >= 15 && mass < 64) {
   index_pt = Find_Index1(muon_pt, true);
   index_rap = Find_Index1(muon_rap, false);
 } else if (mass >= 64 && mass < 120) {
   index_pt = Find_Index2(muon_pt, true);
   index_rap = Find_Index2(muon_rap, false);
 } else if (mass >= 120 && mass < 3000) {
   index_pt = Find_Index3(muon_pt, true);
   index_rap = Find_Index3(muon_rap, false);
 }

//#1
/*
double weights_fewz_1520[5][10] = {{
0.973776, 1.11492, 1.23463, 1.26452, 1.37745, 1.42186, 1.63654, 1.9261, 2.45336, 3.43755}, {
0.982913, 1.11334, 1.21893, 1.27041, 1.31317, 1.44422, 1.61932, 1.87906, 2.42645, 3.14534}, {
0.989479, 1.12629, 1.26198, 1.29244, 1.36864, 1.53303, 1.578, 1.90803, 2.23205, 2.51361}, {
1.00192, 1.14531, 1.27534, 1.31944, 1.38816, 1.42433, 1.61885, 1.93642, 2.46333, 2.78727}, {
1.01619, 1.18475, 1.3078, 1.34552, 1.46668, 1.50327, 1.66022, 1.88626, 2.30752, 2.4581}};
*/

// The updated weights for 15-20 GeV mass bin became known on April 23, 2014
// The changes in the first 3 columns resulted from the correction closure test
double weights_fewz_1520[5][10] = {{
0.90000, 1.23000, 1.25463, 1.26452, 1.37745, 1.42186, 1.63654, 1.9261, 2.45336, 3.43755}, {
0.90000, 1.24000, 1.26893, 1.27041, 1.31317, 1.44422, 1.61932, 1.87906, 2.42645, 3.14534}, {
0.90500, 1.27000, 1.28198, 1.29244, 1.36864, 1.53303, 1.57800, 1.90803, 2.23205, 2.51361}, {
0.90700, 1.29000, 1.30534, 1.31944, 1.38816, 1.42433, 1.61885, 1.93642, 2.46333, 2.78727}, {
0.90800, 1.32000, 1.33780, 1.34552, 1.46668, 1.50327, 1.66022, 1.88626, 2.30752, 2.4581}};

//#2
double weights_fewz_2025[5][10] = {{
0.973776, 1.11492, 1.23463, 1.26452, 1.37745, 1.42186, 1.63654, 1.9261, 2.45336, 3.43755}, {
0.982913, 1.11334, 1.21893, 1.27041, 1.31317, 1.44422, 1.61932, 1.87906, 2.42645, 3.14534}, {
0.989479, 1.12629, 1.26198, 1.29244, 1.36864, 1.53303, 1.578, 1.90803, 2.23205, 2.51361}, {
1.00192, 1.14531, 1.27534, 1.31944, 1.38816, 1.42433, 1.61885, 1.93642, 2.46333, 2.78727}, {
1.01619, 1.18475, 1.3078, 1.34552, 1.46668, 1.50327, 1.66022, 1.88626, 2.30752, 2.4581}};

//#3
double weights_fewz_2530[5][10] = {{
0.968567, 1.07666, 1.18933, 1.2214, 1.26716, 1.35534, 1.4976, 1.83792, 2.27374, 3.26249}, {
0.977465, 1.08471, 1.20161, 1.23233, 1.29367, 1.36965, 1.49842, 1.79827, 2.19182, 3.35555}, {
0.986043, 1.10275, 1.23155, 1.2607, 1.31966, 1.41632, 1.51134, 1.77526, 2.07305, 2.36353}, {
1.0014, 1.13344, 1.25974, 1.30118, 1.36016, 1.41012, 1.60326, 1.85629, 2.31205, 2.51256}, {
1.01962, 1.19647, 1.3294, 1.35609, 1.44448, 1.51468, 1.6651, 1.87865, 2.34528, 2.80925}};

//#4
double weights_fewz_3035[5][10] = {{
0.973858, 1.05056, 1.14885, 1.17712, 1.2215, 1.29861, 1.4307, 1.70399, 2.11095, 3.00655}, {
0.982818, 1.05866, 1.16757, 1.18687, 1.24991, 1.30803, 1.42796, 1.67992, 2.05929, 3.02456}, {
0.991557, 1.07604, 1.18376, 1.2106, 1.26623, 1.34936, 1.44128, 1.66128, 1.95728, 2.23765}, {
1.00573, 1.10398, 1.22029, 1.25223, 1.30114, 1.34194, 1.52048, 1.72962, 2.13181, 2.34171}, {
1.02373, 1.16148, 1.27929, 1.3026, 1.37841, 1.44013, 1.57496, 1.76865, 2.19855, 2.62976}};

//#5
double weights_fewz_3540[5][10] = {{
0.989732, 0.972282, 1.02744, 1.0443, 1.08451, 1.12844, 1.22999, 1.30222, 1.62258, 2.23873}, {
0.998875, 0.980537, 1.06542, 1.05049, 1.11863, 1.12317, 1.21656, 1.32485, 1.66169, 2.0316}, {
1.0081, 0.995935, 1.04037, 1.0603, 1.10594, 1.14847, 1.23108, 1.31935, 1.60998, 1.86}, {
1.01873, 1.01562, 1.10191, 1.10536, 1.12411, 1.1374, 1.27217, 1.3496, 1.59109, 1.82917}, {
1.03605, 1.05651, 1.12896, 1.1421, 1.1802, 1.21647, 1.30452, 1.43867, 1.75835, 2.09128}};

//#6
double weights_fewz_4045[5][10] = {{
0.989732, 0.972282, 1.02744, 1.0443, 1.08451, 1.12844, 1.22999, 1.30222, 1.62258, 2.23873}, {
0.998875, 0.980537, 1.06542, 1.05049, 1.11863, 1.12317, 1.21656, 1.32485, 1.66169, 2.0316}, {
1.0081, 0.995935, 1.04037, 1.0603, 1.10594, 1.14847, 1.23108, 1.31935, 1.60998, 1.86}, {
1.01873, 1.01562, 1.10191, 1.10536, 1.12411, 1.1374, 1.27217, 1.3496, 1.59109, 1.82917}, {
1.03605, 1.05651, 1.12896, 1.1421, 1.1802, 1.21647, 1.30452, 1.43867, 1.75835, 2.09128}};

//#7
double weights_fewz_4550[5][10] = {{
1.00604, 0.956945, 0.994145, 1.00469, 1.03918, 1.05699, 1.13836, 1.21451, 1.4722, 1.88738}, {
1.014, 0.96329, 1.01682, 1.01288, 1.03632, 1.05203, 1.1159, 1.21812, 1.49433, 1.66533}, {
1.02226, 0.976381, 1.01327, 1.01019, 1.04815, 1.0957, 1.13106, 1.2384, 1.45719, 1.56474}, {
1.02934, 0.998431, 1.04915, 1.06283, 1.06603, 1.10242, 1.17645, 1.2618, 1.48553, 1.67461}, {
1.04378, 1.02798, 1.08257, 1.07321, 1.11179, 1.14204, 1.21602, 1.34332, 1.55175, 1.90819}};

//#8
double weights_fewz_5055[5][10] = {{
1.02235, 0.941608, 0.960851, 0.965079, 0.993853, 0.985536, 1.04673, 1.1268, 1.32183, 1.53602}, {
1.02912, 0.946042, 0.968216, 0.975277, 0.954012, 0.980893, 1.01523, 1.11139, 1.32697, 1.29906}, {
1.03643, 0.956827, 0.986172, 0.96008, 0.99036, 1.04292, 1.03104, 1.15745, 1.3044, 1.26948}, {
1.03995, 0.981242, 0.996395, 1.02031, 1.00794, 1.06744, 1.08073, 1.17401, 1.37997, 1.52005}, {
1.05151, 0.99945, 1.03619, 1.00432, 1.04338, 1.06761, 1.12752, 1.24797, 1.34516, 1.72509}};

//#9
double weights_fewz_5560[5][10] = {{
1.02959, 0.941542, 0.95782, 0.952824, 0.984287, 0.987823, 1.02065, 1.10401, 1.27755, 1.4605}, {
1.03567, 0.946896, 0.965881, 0.965349, 0.958662, 0.972948, 1.01625, 1.09353, 1.29119, 1.35351}, {
1.04241, 0.960357, 0.982074, 0.95006, 0.978621, 1.02562, 1.01776, 1.13212, 1.27884, 1.283}, {
1.04794, 0.9791, 0.992944, 1.00913, 1.00984, 1.04553, 1.05514, 1.14065, 1.33224, 1.48418}, {
1.05498, 0.997954, 1.02645, 0.997674, 1.03357, 1.04947, 1.11478, 1.22201, 1.31028, 1.57255}};

//#10
double weights_fewz_6064[5][10] = {{
1.04407, 0.941412, 0.95176, 0.928316, 0.965154, 0.992398, 0.968476, 1.05845, 1.18899, 1.30947}, {
1.04878, 0.948603, 0.961211, 0.945493, 0.967963, 0.957059, 1.01829, 1.05779, 1.21961, 1.46242}, {
1.05437, 0.967419, 0.973877, 0.930019, 0.955141, 0.991, 0.991204, 1.08146, 1.22772, 1.31002}, {
1.06394, 0.974816, 0.986042, 0.98677, 1.01364, 1.00171, 1.00396, 1.07394, 1.23678, 1.41244}, {
1.06194, 0.994962, 1.00698, 0.98438, 1.01396, 1.01317, 1.0893, 1.17009, 1.24054, 1.26749}};

//#11
double weights_fewz_6468[3][10] = {{
1.05132, 0.941347, 0.948729, 0.916061, 0.955588, 0.994685, 0.942392, 1.03567, 1.14471, 1.23395}, {
1.05867, 0.963436, 0.965989, 0.925308, 0.953306, 0.965062, 0.991877, 1.05045, 1.19567, 1.38837}, {
1.06718, 0.986544, 0.992149, 0.97697, 1.00828, 0.989236, 1.03573, 1.09701, 1.19699, 1.28071}};

//#12
double weights_fewz_6872[3][10] = {{
1.05494, 0.958128, 0.95965, 0.927253, 0.959492, 0.989889, 0.955272, 1.03574, 1.14345, 1.22829}, {
1.06191, 0.9751, 0.975414, 0.937342, 0.960194, 0.970575, 0.994642, 1.04966, 1.18674, 1.35406}, {
1.07142, 0.995032, 0.999609, 0.981932, 1.0089, 0.997558, 1.0368, 1.0949, 1.19686, 1.2757}};

//#13
double weights_fewz_7276[3][10] = {{
1.05856, 0.974908, 0.97057, 0.938446, 0.963396, 0.985093, 0.968151, 1.03581, 1.1422, 1.22262}, {
1.06515, 0.986763, 0.984839, 0.949375, 0.967082, 0.976087, 0.997408, 1.04886, 1.1778, 1.31974}, {
1.07565, 1.00352, 1.00707, 0.986895, 1.00952, 1.00588, 1.03786, 1.09279, 1.19673, 1.27068}};

//#14
double weights_fewz_7681[3][10] = {{
1.06217, 0.991689, 0.981491, 0.949639, 0.967299, 0.980297, 0.98103, 1.03587, 1.14095, 1.21695}, {
1.06839, 0.998426, 0.994264, 0.961408, 0.97397, 0.9816, 1.00017, 1.04806, 1.16886, 1.28543}, {
1.07988, 1.01201, 1.01453, 0.991858, 1.01014, 1.01421, 1.03892, 1.09068, 1.1966, 1.26567}};

//#15
double weights_fewz_8186[3][10] = {{
1.06398, 1.00008, 0.986951, 0.955235, 0.969251, 0.977899, 0.98747, 1.03591, 1.14032, 1.21412}, {
1.07001, 1.00426, 0.998977, 0.967425, 0.977414, 0.984357, 1.00155, 1.04766, 1.16439, 1.26827}, {
1.082, 1.01626, 1.01826, 0.994339, 1.01045, 1.01837, 1.03946, 1.08962, 1.19653, 1.26316}};

//#16
double weights_fewz_8691[3][10] = {{
1.06579, 1.00847, 0.992411, 0.960831, 0.971203, 0.9755, 0.99391, 1.03594, 1.13969, 1.21128}, {
1.07164, 1.01009, 1.00369, 0.973442, 0.980858, 0.987113, 1.00294, 1.04726, 1.15992, 1.25111}, {
1.08412, 1.0205, 1.022, 0.996821, 1.01076, 1.02253, 1.03999, 1.08856, 1.19647, 1.26066}};
//#17
double weights_fewz_9196[3][10] = {{
1.06579, 1.00847, 0.992411, 0.960831, 0.971203, 0.9755, 0.99391, 1.03594, 1.13969, 1.21128}, {
1.07164, 1.01009, 1.00369, 0.973442, 0.980858, 0.987113, 1.00294, 1.04726, 1.15992, 1.25111}, {
1.08412, 1.0205, 1.022, 0.996821, 1.01076, 1.02253, 1.03999, 1.08856, 1.19647, 1.26066}};
//#18
double weights_fewz_96101[3][10] = {{
1.06579, 1.00847, 0.992411, 0.960831, 0.971203, 0.9755, 0.99391, 1.03594, 1.13969, 1.21128}, {
1.07164, 1.01009, 1.00369, 0.973442, 0.980858, 0.987113, 1.00294, 1.04726, 1.15992, 1.25111}, {
1.08412, 1.0205, 1.022, 0.996821, 1.01076, 1.02253, 1.03999, 1.08856, 1.19647, 1.26066}};
//#19
double weights_fewz_101106[3][10] = {{
1.06579, 1.00847, 0.992411, 0.960831, 0.971203, 0.9755, 0.99391, 1.03594, 1.13969, 1.21128}, {
1.07164, 1.01009, 1.00369, 0.973442, 0.980858, 0.987113, 1.00294, 1.04726, 1.15992, 1.25111}, {
1.08412, 1.0205, 1.022, 0.996821, 1.01076, 1.02253, 1.03999, 1.08856, 1.19647, 1.26066}};
//#20
double weights_fewz_106110[3][10] = {{
1.07192, 1.01304, 1.00893, 0.979635, 0.99024, 0.993835, 0.960001, 0.987467, 1.05642, 1.15799}, {
1.07959, 1.01509, 0.984092, 0.960777, 0.973696, 0.9924, 0.965069, 0.987102, 1.09771, 1.08016}, {
1.08857, 1.0157, 1.01658, 0.993233, 1.01552, 1.01909, 0.975677, 1.02952, 1.1545, 1.27262}};
//#21
double weights_fewz_110115[3][10] = {{
1.07192, 1.01304, 1.00893, 0.979635, 0.99024, 0.993835, 0.960001, 0.987467, 1.05642, 1.15799}, {
1.07959, 1.01509, 0.984092, 0.960777, 0.973696, 0.9924, 0.965069, 0.987102, 1.09771, 1.08016}, {
1.08857, 1.0157, 1.01658, 0.993233, 1.01552, 1.01909, 0.975677, 1.02952, 1.1545, 1.27262}};
//#22
double weights_fewz_115120[3][10] = {{
1.08058, 1.01161, 0.987128, 0.939578, 0.95139, 0.954759, 0.952376, 0.971983, 1.03031, 1.09207}, {
1.07969, 1.0177, 0.993776, 0.968835, 0.991958, 0.924721, 0.971056, 0.986965, 1.08477, 1.0452}, {
1.08952, 1.01472, 1.00701, 0.95694, 0.991331, 0.988923, 0.956525, 1.01083, 1.13952, 1.13622}};

//#23
double weights_fewz_120126[2][3] = {{
1.08058, 0.982579, 1.05158}, {
1.08458, 0.996471, 1.10972}};

//#24
double weights_fewz_126133[2][3] = {{
1.08975, 0.976689, 1.03077}, {
1.08595, 0.999086, 1.06066}};

//#25
double weights_fewz_133141[2][3] = {{
1.08975, 0.976689, 1.03077}, {
1.08595, 0.999086, 1.06066}};

//#26
double weights_fewz_141150[2][3] = {{
1.09028, 0.986419, 1.059}, {
1.08688, 1.00759, 1.03222}};
//#27
double weights_fewz_150160[2][3] = {{
1.08947, 0.992888, 1.00001}, {
1.09577, 0.99148, 1.01254}};
//#28
double weights_fewz_160171[2][3] = {{
1.09196, 1.00245, 0.943417}, {
1.0938, 0.978614, 0.990971}};
//#29
double weights_fewz_171185[2][3] = {{
1.09484, 1.00878, 0.916319}, {
1.0874, 0.973804, 0.979243}};
//#30
double weights_fewz_185200[2][3] = {{
1.08378, 1.00039, 0.91396}, {
1.0889, 0.969653, 0.95953}};
//#31
double weights_fewz_200220[2][3] = {{
1.06166, 0.983607, 0.909244}, {
1.0919, 0.961353, 0.920104}};
//#32
double weights_fewz_220243[2][3] = {{
1.0506, 0.975215, 0.906886}, {
1.0934, 0.957202, 0.900391}};
//#33
double weights_fewz_243273[2][3] = {{
1.0506, 0.975215, 0.906886}, {
1.0934, 0.957202, 0.900391}};
//#34
double weights_fewz_273320[2][3] = {{
1.0506, 0.975215, 0.906886}, {
1.0934, 0.957202, 0.900391}};
//#35
double weights_fewz_320380[2][3] = {{
1.0506, 0.975215, 0.906886}, {
1.0934, 0.957202, 0.900391}};
//#36
double weights_fewz_380440[2][3] = {{
1.0506, 0.975215, 0.906886}, {
1.0934, 0.957202, 0.900391}};
//#37
double weights_fewz_440510[2][3] = {{
1.05722, 0.979745, 0.867872}, {
1.07688, 0.963899, 0.893238}};
//#38
double weights_fewz_510600[2][3] = {{
1.0273, 0.998866, 0.880407}, {
1.07276, 0.966222, 0.859189}};
//#39
double weights_fewz_6001000[2][3] = {{
1.03718, 1.00361, 0.878818}, {
1.05232, 0.964508, 0.804247}};
//#40
double weights_fewz_10001500[2][3] = {{
0.973859, 1.01897, 0.826874}, {
1.02721, 0.945956, 0.81459}};
//#41
double weights_fewz_15002000[2][3] = {{
1.01542, 1.009672, 0.819218}, {
1.02512, 0.940146, 0.803916}};

    if (mass >= 15 && mass < 20) {
       return weights_fewz_1520[index_rap][index_pt];
    } else if (mass >= 20 && mass < 25) {
      return weights_fewz_2025[index_rap][index_pt];
    } else if (mass >= 25 && mass < 30) {
      return weights_fewz_2530[index_rap][index_pt];
    } else if (mass >= 30 && mass < 35) {
      return weights_fewz_3035[index_rap][index_pt];
    } else if (mass >= 35 && mass < 40) {
      return weights_fewz_3540[index_rap][index_pt];
    } else if (mass >= 40 && mass < 45) {
      return weights_fewz_4045[index_rap][index_pt];
    } else if (mass >= 45 && mass < 50) {
      return weights_fewz_4550[index_rap][index_pt];
    } else if (mass >= 50 && mass < 55) {
      return weights_fewz_5055[index_rap][index_pt];
    } else if (mass >= 55 && mass < 60) {
      return weights_fewz_5560[index_rap][index_pt]; 
    } else if (mass >= 60 && mass < 64) {
      return weights_fewz_6064[index_rap][index_pt];
    } else if (mass >= 64 && mass < 68) {
      return weights_fewz_6468[index_rap][index_pt];
    } else if (mass >= 68 && mass < 72) {
      return weights_fewz_6872[index_rap][index_pt];
    } else if (mass >= 72 && mass < 76) {
      return weights_fewz_7276[index_rap][index_pt];
    } else if (mass >= 76 && mass < 81) {
      return weights_fewz_7681[index_rap][index_pt];
    } else if (mass >= 81 && mass < 86) {
      return weights_fewz_8186[index_rap][index_pt];
    } else if (mass >= 86 && mass < 91) {
      return weights_fewz_8691[index_rap][index_pt];
    } else if (mass >= 91 && mass < 96) {
      return weights_fewz_9196[index_rap][index_pt];
    } else if (mass >= 96 && mass < 101) {
      return weights_fewz_96101[index_rap][index_pt];
    } else if (mass >= 101 && mass < 106) {
      return weights_fewz_101106[index_rap][index_pt];
    } else if (mass >= 106 && mass < 110) {
      return weights_fewz_106110[index_rap][index_pt];
    } else if (mass >= 110 && mass < 115) {
      return weights_fewz_110115[index_rap][index_pt];
    } else if (mass >= 115 && mass < 120) {
        return weights_fewz_115120[index_rap][index_pt];
    } else if (mass >= 120 && mass < 126) {
        return weights_fewz_120126[index_rap][index_pt];
    } else if (mass >= 126 && mass < 133) {
      return weights_fewz_126133[index_rap][index_pt];
    } else if (mass >= 133 && mass < 141) {
      return weights_fewz_133141[index_rap][index_pt];
    } else if (mass >= 141 && mass < 150) {
      return weights_fewz_141150[index_rap][index_pt];
    } else if (mass >= 150 && mass < 160) {
      return weights_fewz_150160[index_rap][index_pt];
    } else if (mass >= 160 && mass < 171) {
      return weights_fewz_160171[index_rap][index_pt];
    } else if (mass >= 171 && mass < 185) {
      return weights_fewz_171185[index_rap][index_pt];
    } else if (mass >= 185 && mass < 200) {
       return weights_fewz_185200[index_rap][index_pt];
    } else if (mass >= 200 && mass < 220) {
       return weights_fewz_200220[index_rap][index_pt];
    } else if (mass >= 220 && mass < 243) {
       return weights_fewz_220243[index_rap][index_pt];
    } else if (mass >= 243 && mass < 273) {
       return weights_fewz_243273[index_rap][index_pt];
    } else if (mass >= 273 && mass < 320) {
       return weights_fewz_273320[index_rap][index_pt];
    } else if (mass >= 320 && mass < 380) {
       return weights_fewz_320380[index_rap][index_pt];
    } else if (mass >= 380 && mass < 440) {
       return weights_fewz_380440[index_rap][index_pt];
    } else if (mass >= 440 && mass < 510) {
       return weights_fewz_440510[index_rap][index_pt];
    } else if (mass >= 510 && mass < 600) {
       return weights_fewz_510600[index_rap][index_pt];
    } else if (mass >= 600 && mass < 1000) {
      return weights_fewz_6001000[index_rap][index_pt];
    } else if (mass >= 1000 && mass < 1500) {
      return weights_fewz_10001500[index_rap][index_pt];
    } else if (mass >= 1500 && mass < 2000) {
      return weights_fewz_15002000[index_rap][index_pt];
    }

    return 1.;
}

// end of A macro from Alexey Svyatkovskiy <asvyatko@purdue.edu> providing


#endif

// ----------------------------------------------------------
// ----------------------------------------------------------

