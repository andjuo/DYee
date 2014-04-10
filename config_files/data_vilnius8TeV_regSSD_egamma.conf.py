# kBlack        1  
# kGray       920
# kWhite        0
# kBlue       600
# kAzure      860
# kCyan       432
# kTeal       840
# kGreen      416
# kSpring     820
# kYellow     400
# kOrange     800
# kRed        632
# kPink       900
# kMagenta    616
# kViolet     880
#
# CAUTION: it is ok to put comments and the end of only some of these lines.
#          Also, there can be no comments between the lines from this point until % sign.

[DEFINITIONS]
PATH_1=/media/sdb3/CMSData2012-regressed
PATH_2=/media/sdb3/CMSData2012-regressed
PATH_3=/media/sdb/CMSData2012-regressed2
PATH_4=
PATH_5=../config_files/JSONs/8TeV/
PATH_SSD=/media/ssd
NTUPLE=ntuple
SKIM=tight-loose_skim
LUMI=19712

# usual RESULT_ROOT_DIR is ../root_files_reg
#RESULT_ROOT_DIR=../../Results-DYee/root_files_reg/


[CONSTANTS]
# luminosity [pb^-1]
Lumi=19712.

# SelectEventsFlag determines how event many events are taken
# 0 => select number of events as expected from luminosity; 
# 1 => weight all events by luminosity     
SelectEventsFlag=1  

# SelectionTag specifies subdirectory for structures
SelectionTag   = DY_j22_19712pb_egamma
ConstantsTag   = DY_j22_19712pb_egamma
XSecTag        = DY_j22_19712pb

# Auxiliary special tags to mark files
SpecTagDirUser=
SpecTagFileUser=

# Electron energy scale
#EScaleSet=Date20130529_2012_j22_adhoc   # Name of energy scale calibrations set. See ElectronEnergyScale.hh.
EScaleSet=UNCORRECTED

# Extension of the saved plots
PlotExtension=png  # png,pdf,C,root, default={png,pdf,root}

# Trigger set
Trigger= Full2012_hltEffOld

# FEWZ correction flag
FEWZ=1

# Pile-up reweight flag
PUReweight=1

# User keys ---------
seedMin=1001
seedMax=1020
Use7TeVMCWeight=1 # Whether renormalize the 1st sample weight to 1 in MC macros
IgnoreDebugRunForYields=1 # prepareYields shouldn't use DebugRun ntuples

# Background type
DDBkg=1
#DDBkgVersion=20131231
#DDBkgVersion=20140301AJ  # AJ
DDBkgVersion=20140312  # from Manny through Alexey

#T&P_ESF_extra=_etaMax24
ScaleFactorTag = DY_j22_19712pb_egamma_Unregressed_energy

[TAG_AND_PROBE]
# Tag and probe options

#T&P_TAG=Date20130529_2012_j22_adhoc // ignored  // SelectionTag is used
#MAP = ETBINS6 ETABINS5
MAP = ETBINS6 ETABINS5egamma


# How to measure efficiencies? COUNTnCOUNT, COUNTnFIT, FITnFIT for pass/fail samples
# Note: for template-based fit of the data for ID efficiency, use FITnFIT only
# Note: for HLT efficiencies both data and MC are sufficiently background free
#   for using count and count method
RECO_DATA=FITnFIT
ID_DATA  =FITnFIT
HLT_DATA =COUNTnCOUNT

RECO_MC  =COUNTnCOUNT
ID_MC    =COUNTnCOUNT
HLT_MC   =COUNTnCOUNT

%

#
#  DATA NTUPLES
#
[DATA_FILES]
$ data 1  1 1 1  @data
PATH_SSD/r12a-del-j22-v1_tight-loose_skim.root  0 PATH_5/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt
PATH_SSD/r12b-del-j22-v1_tight-loose_skim.root  0 PATH_5/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt
PATH_SSD/r12c-del-j22-v1_tight-loose_skim.root  0 PATH_5/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt
PATH_SSD/r12d-del-j22-v1_tight-loose_skim.root  0 PATH_5/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON.txt
%


[BACKGROUND_FILES]

#
# BACKGROUNDS: TTBAR
#
$ ttbar 814 @t#bar{t}
#
# Inclusive NNLO+NNLL x-sec for ttbar at 8 TeV: 245.8 pb
# Reference: http://arxiv.org/pdf/1303.6254v1.pdf
# This is recommended by Rebeca Gonzalez Suarez from top group as a better
# value to use than the official (NLO only) CMS reference
#   https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV
#
# Note: we have powheg and madgraph ntuples. The madgraph sample is the one
# recommended by the top group, the powheg is for (if needed) systematics studies.
#
# If we use a singe sample for ttbar, uncomment the line below.
# Cross section for ttbar->2lX:
#    245.8 * *(0.1080*3)*(0.1080*3) = 25.80
#PATH_2/s12-tt2l-mg-v7c_tight-loose_skim.root 25.80
#
# If we use three samples for mass ranges low-700-1000-inf, use the ntuples
# below.
# Cross section for ttbar->2lX, with restriction of M(tt)<700 GeV:
#     at gen-level, 92.43% of events have M(tt)<700, so
#     xsec = 25.80 * 0.9243 = 23.85 pb
# Cross section for ttbar inclusive, M(tt)>=700 && M(tt)<1000
#     from PREP, the filter efficiency for this sample is 0.074,
#     and before the filter we use the total ttbar inclusive cross section, so
#     xsec = 245.8 * 0.074 = 18.19 pb
#     (note: for madgraph ttbar->2lX sample above, the fraction >700 is a bit different,
#      but the PREP number is used here).
# Cross section for ttbar inclusive, M(tt)>=1000
#     same considerations as above, only PREP filter efficiency is 0.014.
#     xsec = 245.8 * 0.014 = 3.44 pb
# NOTE: the first sample is 2lX, and the other two are inclusive, so cross sections can't be
# directly compared (i.e. if the first sample was inclusive too, its cross section would
# be 245.8 * 0.9243 with truncation at 700 GeV).
PATH_2/s12-tt2l-mg-v7c_massLowTo700_tight-loose_skim.root 23.85
PATH_2/s12-ttj-m700-1000-pwg-v7a_tight-loose_skim.root 18.19
PATH_2/s12-ttj-m1000-pwg-v7a_tight-loose_skim.root 3.44
#
# BACKGROUNDS: W+JETS
#
$ wjets 894 @W+jets
#
# FEWZ NNLO cross section is 37509.0 pb
# taken from the  CMS official reference:
#    https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV
#
PATH_2/s12-wjets-v7a_tight-loose_skim.root 37509.0
#
# BACKGROUNDS: DIBOSONS
#
$ ww 904 @WW
# From MIT's UserCode/MitPhysics/data/xs.dat the cross section is 47.0*(1-0.0305)*(0.1080*3)*(0.1080*3)*1.2151 = 5.812
PATH_2/s12-wwj-v7a_ntuple.root 5.812
$ wz 905 @WZ
#
# In case of inclusive sample:
# From MIT's UserCode/MitPhysics/data/xs.dat the cross section is 18.2*1.23344 = 22.4
# However, this needs to be re-thought (if we really plan to use this number).
# This cross section diverges, and quoted numbers involve cuts. Not clear what it should be.
# PATH_2/s12-wz-v7a_tight-loose_skim.root = 22.4
#
# It is more efficient to use exclusive samples. 
# Cross section for WZ->3ln
#   we start with WZ (m(ll)>12) MCFM NLO number, computed with MSTW PDF from
#     https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV
#   which is 33.85 pb. For exclusive 3ln case:
#   xsec = 33.84 *(0.1080*3)*0.101 = 1.11 pb
# Cross section for WZ->2l2q
#   same idea as above. The exclusive value is:
#   xsec = 33.84 * (0.676)*0.101 = 2.31 pb
PATH_2/s12-wz3ln-v7a_ntuple.root 1.11
PATH_2/s12-wz2q2l-v7a_tight-loose_skim.root 2.31
$ zz 906 @ZZ
# From MIT's UserCode/MitPhysics/data/xs.dat the cross section is 0.364791
PATH_2/s12-zz2l2n-v7a_ntuple.root 0.364791
# From MIT's UserCode/MitPhysics/data/xs.dat the cross section is 2.448690
PATH_2/s12-zz2l2q-v7a_ntuple.root 2.448690
# From MIT's UserCode/MitPhysics/data/xs.dat the cross section is 0.176908
PATH_2/s12-zz4l-v7a_tight-loose_skim.root 0.176908
#
# BACKGROUNDS: Z->TAU TAU
#
$ ztt 855 @Z#rightarrow#tau#tau
# The cross sections for these two samples are exactly the same as for the DYToEE primary signal
# samples. See the origin of the numbers in that section below
PATH_2/s12-ztt1020-v7a_tight-loose_skim.root 3795.4
PATH_2/s12-zttm20-v7a_tight-loose_skim.root 1915.1
#
# BACKGROUNDS: QCD
#
$ qcd 797 @QCD
#
# QCD EMEnriched: six samples are used, skimmed,
# including: 20-30, 30-80, 80-170, 170-250, 250-350, 350+.
# The cross sections are taken from PREP
#           (xsec) * (filter efficiency)
#   20-30   2.886e8 * 0.0101    = 2914000
#   30-80   7.433e7 * 0.0621    = 4616000
#   80-170  1191000.0 * 0.1539  =  183300
#  170-250  30990.0 * 0.148     =    4587
#  250-350  4250.0 * 0.131      =     556.8
#  350-inf  810.0 * 0.11        =      89.10
PATH_2/s12-qcdem2030-v7a_tight-loose_skim.root 2914000
PATH_2/s12-qcdem3080-v7a_tight-loose_skim.root 4616000
PATH_2/s12-qcdem80170-v7a_tight-loose_skim.root 183300
PATH_2/s12-qcdem170250-v7a_tight-loose_skim.root 4587
PATH_2/s12-qcdem250350-v7a_tight-loose_skim.root 556.8
PATH_2/s12-qcdem350-v7a_tight-loose_skim.root 89.10
#
# QCD BCToE: six samples are used, skimmed,
# including: 20-30, 30-80, 80-170, 170-250, 250-350, 350+.
# The cross sections are taken from PREP
#           (xsec) * (filter efficiency)
#   20-30   2.886e8 * 5.8e-4   = 167400
#   30-80   7.424e7 * 0.00225  = 167000
#   80-170  1191000.0 * 0.0109 =  12980
#  170-250  30980.0 * 0.0204   =   632.0
#  250-350  4250.0 * 0.0243    =   103.3
#  350-inf  811.0 * 0.0295     =    23.92
PATH_2/s12-qcdbc2030-v7a_tight-loose_skim.root 167400
PATH_2/s12-qcdbc3080-v7a_tight-loose_skim.root 167000
PATH_2/s12-qcdbc80170-v7a_tight-loose_skim.root  12980
PATH_2/s12-qcdbc170250-v7a_tight-loose_skim.root 632.0
PATH_2/s12-qcdbc250350-v7a_tight-loose_skim.root 103.3
PATH_2/s12-qcdbc350-v7a_tight-loose_skim.root 23.92


[MC_SIGNAL]
#
# SIGNAL MC Z->EE, should be last entry
#
$ zee 426  798 798 803  @Z#rightarrowee

# 
# Explanation of the cross section values used below:
# The DY->ee Powheg NLO cross sections, taken from PREP from
#     http://cms.cern.ch/iCMS/jsp/mcprod/admin/requestmanagement.jsp?dsn=*DYToEE*&campid=Summer12_DR53X_with_regression
#  10-20  : 3708.0 pb
#  20-inf : 1871.0 pb
#  200-inf: 1.483  pb
#  400-inf: 0.1085 pb
#  500-inf: 0.04409 pb
#  700-inf: 0.01025 pb
#  800-inf: 0.005491 pb
# 1000-inf: 0.001796 pb
# 1500-inf: 1.705E-4 pb
# 2000-inf: 2.208E-5 pb
#
# The NNLO for DY->LL (L=e,mu,tau) from generator group is 5745.25 pb
#     https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV
# which means DY->ee is (1/3) of it = 1915.1 pb
#
# We compute the approximation to NNLO for each range as
#    NNLO-approx(range) = (Powheg NLO (range) ) * (NNLO M>20) / (Powheg NLO (range) )
# Thus:
#  10-20  : 3708.0 * (1915.1/1871.0) = 3795.4
#  20-inf : 1871.0 * (1915.1/1871.0) = 1915.1
#  200-inf: 1.483  * (1915.1/1871.0) = 1.5180
#  400-inf: 0.1085 * (1915.1/1871.0) = 0.1111
#  500-inf: 0.04409  * (1915.1/1871.0) = 0.045129
#  700-inf: 0.01025  * (1915.1/1871.0) = 0.01049
#  800-inf: 0.005491 * (1915.1/1871.0) = 0.005620
# 1000-inf: 0.001796 * (1915.1/1871.0) = 0.001838
# 1500-inf: 1.705E-4 * (1915.1/1871.0) = 1.745E-4
# 2000-inf: 2.208E-5 * (1915.1/1871.0) = 2.260E-5
#
# Additionally:
#  20-500: xsec(>20) - xsec(>500)   = 1915.1 - 0.045129 ~= 1915.1
#  500-800: xsec(>500) - xsec(>800) = 0.045129 - 0.005620 ~= 0.039509
#
#    20- 200: xsec(  >20) - xsec( >200) = 1915.1 - 1.5180     = 1913.6
#   200- 400: xsec( >200) - xsec( >400) = 1.5180 - 0.1111     = 1.4069
#   400- 500: xsec( >400) - xsec( >500) = 0.1111 - 0.045129   = 0.06597
#   500- 700: xsec( >500) - xsec( >700) = 0.045129 - 0.01049  = 0.03464
#   700- 800: xsec( >700) - xsec( >800) = 0.01049  - 0.005620 = 0.004870
#   800-1000: xsec( >800) - xsec(>1000) = 0.005620 - 0.001838 = 0.003782
#  1000-1500: xsec(>1000) - xsec(>1500) = 0.001838 - 1.745E-4 = 0.001664
#  1500-2000: xsec(>1500) - xsec(>2000) = 1.745E-4 - 2.260E-5 = 0.0001519
#
# -----
# We will use either 10-20 + 20-inf GeV samples, or 10-20 + 20-500 + 500-800 + 800-inf.
# AJ: Ilya uses scaled CS values, I put here NLO CS values
# -----
# In case we use 20-inf sample:
#PATH_2/s12-zeem20-v7a_tight-loose_skim.root 1915.1
#
# In case we use 20-500-800-inf:
#
#PATH_SSD/s12-zeem20to500-v7a_tight-loose_skim.root 1915.1  600  @Powheg 20-500
#
#PATH_2/s12-zeem500to800-v7a_ntuple.root 0.039509  632 @Powheg 500-800
#
#PATH_2/s12-zeem800-v7a_ntuple.root 0.005620   432 @Powheg 800+
#
# In case we use 20-500-inf:
#/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/ikrav/DrellYan_8TeV/53X_with_regression/s12-zeem500-v7a_ntuple.root 0.045129
#
PATH_SSD/s12-zeem1020-v7a_tight-loose_skim.root 3708.0  800 @Powheg 10-20
#
# new n-tuples with finer binning
#
PATH_SSD/s12-zeem20to200-v7a_tight-loose_skim.root 1869.5   600 1 @Powheg 20-200
PATH_SSD/s12-zeem200to400-v7a_ntuple.root   1.3745     632 1 @Powheg 200-400
PATH_SSD/s12-zeem400to500-v7a_ntuple.root   0.06441    432 1 @Powheg 400-500
PATH_SSD/s12-zeem500to700-v7a_ntuple.root   0.03384    616 1 @Powheg 500-700
PATH_SSD/s12-zeem700to800-v7a_ntuple.root   0.004759   880 1 @Powheg 700-800
PATH_SSD/s12-zeem800to1000-v7a_ntuple.root  0.003695   900 1 @Powheg 800-1000
PATH_SSD/s12-zeem1000to1500-v7a_ntuple.root 0.0016255   625 1 @Powheg 1000-1500
PATH_SSD/s12-zeem1500to2000-v7a_ntuple.root 0.00014842  419 1 @Powheg 1500-2000
PATH_SSD/s12-zeem2000-v7a_ntuple.root 0.00002208  802 1 @Powheg 2000-
%
