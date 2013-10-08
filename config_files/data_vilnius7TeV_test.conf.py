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
PATH_1=/opt/CMSData/2011/2ndProd-ntuple
PATH_2=
PATH_3=
PATH_4=
PATH_5=../config_files/JSONs
NTUPLE=ntuple
SKIM=tight-loose_skim
LUMI=4839


[CONSTANTS]

# luminosity [pb^-1]
Lumi=4839.          

# SelectEventsFlag determines how event many events are taken
# 0 => select number of events as expected from luminosity; 
# 1 => weight all events by luminosity     
SelectEventsFlag=1  

# SelectionTag specifies subdirectory for structures
SelectionTag   =DY_20130801_test
ConstantsTag   =DY_20130801_test
XSecTag        =DY_20130801_test

# Auxiliary special tags to mark files
SpecTagDirUser=
SpecTagFileUser=

# Electron energy scale
EScaleSet=Date20120802_default   # Name of energy scale calibrations set. See ElectronEnergyScale.hh.

# Extension of the saved plots
PlotExtension=png  # png,pdf,C,root, default={png,pdf,root}

# Trigger set
Trigger= Full2011_hltEffOld

# FEWZ correction flag
FEWZ=1

# Pile-up reweight flag
PUReweight=1

# User keys
seedMin=1001
seedMax=1020
Use7TeVMCWeight=1 # Whether renormalize the 1st sample weight to 1 in MC macros
IgnoreDebugRunForYields=1 # prepareYields shouldn't use DebugRun ntuples

[TAG_AND_PROBE]
# Tag and probe options

#T&P_TAG=DY_m10+pr+a05+o03+pr_4839pb # ignored
MAP = ETBINS6 ETABINS5

TargetDistr=pileup_lumibased_data
SourceDistr=pileup_simulevel_mc


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
PATH_1/r11a-del-m10-v1_tight-loose_skim.root  0  PATH_5/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3_in_r11a_del_m10_ntuples.txt
PATH_1/r11a-del-pr-v4_tight-loose_skim.root   0  PATH_5/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON_in_r11a_del_pr_v4_ntuples.txt
PATH_1/r11a-del-a05-v1_tight-loose_skim.root  0  PATH_5/Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v3_in_r11a_del_a05_ntuples_exclJuly11Problem.txt
PATH_1/r11a-del-o03-v1_tight-loose_skim.root  0  PATH_5/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON_in_r11a_del_o03_ntuples.txt
PATH_1/r11b-del-pr-v1_tight-loose_skim.root   0  PATH_5/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON_in_r11b_del_pr_ntuples.txt
%


[BACKGROUND_FILES]

# after '$' the name of the sample is followed by 4 colors:
# 1 general color for plot
# and 3 colors for stack plot (fill, marker and line colors);
# the line ends with the label

#
# BACKGROUNDS: TTBAR
#
$ ttbar 814  634 634 636 @t#bar{t}
# -->inclusive ttbar sample
# This sample is the madgraph sample. While it says TTJets, we are told
# this also include 0-jet bin.
PATH_1/f11-ttj-v14b2-bp_tight-loose_skim.root 157.5

#
# BACKGROUNDS: W+JETS
#
$ wjets 894  810 810 803  @W+jets
# This is a single Madgraph sample, factor: 10438.0*3.0
PATH_1/f11-wjets-v14b-bp_tight-loose_skim.root  31314.0

#
# BACKGROUNDS: tW
#
$ wtop 634  46 46 803  @W+t
PATH_1/f11-wtop-powheg-v14b-bp_ntuple.root  7.87
PATH_1/f11-wtopb-powheg-v14b-bp_ntuple.root 7.87

#
# BACKGROUNDS: DIBOSONS
#
# (this is needed if dibosons all go in one contribution -> )$ diboson 906 @diboson
# Diboson decays specifically to leptons
# For WW->2l2nu, the factors are 43.0*(0.1080*3)*(0.1080*3)

$ ww 904  810 810 803    @WW
##/mnt/unl/home/ikrav/work/ntuples/42X/s11-ww2l-v11-pu_ntuple.root	4.514
##/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/ikrav/Full2011_42X_v2/EE/f11-ww2l-v14b-bp_ntuple.root    4.514
#correction: f11-ww2l-v14b-pu             47.0*(1-0.0305)*(0.1080*3)*(0.1080*3)
PATH_1/f11-ww2l-v14b-bp_ntuple.root    4.783

# For WZ->lnull, the factors are 18.2*(0.1080*3)*0.101
# For inclusive WZ the factor is 18.2
$ wz 905  810 810 803  @WZ
#/mnt/unl/home/ikrav/work/ntuples/42X/s11-wz3l-v11-pu_ntuple.root    0.596  
PATH_1/f11-wz-v14b-bp_ntuple.root  18.2

# The ZZ sample is ZZ->Anything, the factors are 6.77*(1+0.12/1.277)
$ zz 906  810 810 803   @ZZ
#/mnt/unl/home/ikrav/work/ntuples/42X/s11-zz-v11-pu_tight-loose_skim.root  7.406
PATH_1/f11-zz-v14b-bp_tight-loose_skim.root  7.406

#
# BACKGROUNDS: Z->TAU TAU
#
$ ztt 855  807 807 803  @Z#rightarrow#tau#tau
#/mnt/unl/home/ikrav/work/ntuples/42X/s11-zttm20-powheg-v11-pu_tight-loose_skim.root 1666.0
PATH_1/f11-zttm20-powheg-v14b-bp_tight-loose_skim.root 1666.0

#
# BACKGROUNDS: QCD
#
$ qcd 797   885 885 883 @QCD
# EM enriched Individual ntuples
# /mnt/unl/home/ikrav/work/ntuples/42X/s11-qcdem2030-v11-pu_ntuple.root  2454400
# /mnt/unl/home/ikrav/work/ntuples/42X/s11-qcdem3080-v11-pu_ntuple.root  3866200
# /mnt/unl/home/ikrav/work/ntuples/42X/s11-qcdem80170-v11-pu_ntuple.root  139500
# the above combined into skim, factor 2454400 + 3866200 + 139500
#/mnt/unl/home/ikrav/work/ntuples/42X/s11-qcdem-v11-pu_tight-loose_skim.root  6460100
PATH_1/f11-qcdem-v14b-bp_tight-loose_skim.root 6460100
# heavy flavor
# /mnt/unl/home/ikrav/work/ntuples/42X/s11-qcdbc2030-v11-pu_ntuple.root   132160
# /mnt/unl/home/ikrav/work/ntuples/42X/s11-qcdbc3080-v11-pu_ntuple.root   136804
# /mnt/unl/home/ikrav/work/ntuples/42X/s11-qcdbc80170-v11-pu_ntuple.root    9360
# the above combined into skim, factor 132160 + 136804 + 9360
#/mnt/unl/home/ikrav/work/ntuples/42X/s11-qcdbc-v11-pu_tight-loose_skim.root   278324
PATH_1/f11-qcdbce-v14b-bp_tight-loose_skim.root 278324


[MC_SIGNAL]
#
# SIGNAL MC Z->EE, should be last entry
#
$ zee 426  798 798 803  @Z#rightarrowee
# We will use either 10-20 + 20-inf GeV samples, or 10-20 + 20-500 + 500-800 + 800-inf.
# 
# In case we use 20-inf sample:
#/mnt/hadoop/user/uscms01/pnfs/unl.edu/data4/cms/store/user/ikrav/Full2011_42X_v2/EE/f11-zeem20-powheg-v14b-bp_tight-loose_skim.root 1666.0
#
# In case we use 20-500-800-inf:
#
# cross-section of trimmed m20-500 sample is (1666.0 - 0.03335) pb ~= 1666.0 pb
PATH_1/f11-zeem20to500-powheg-v14b-bp_tight-loose_skim.root 1666.0 600  @Powheg 20-500
#
# cross-section of trimmed m500-800 sample is (0.03335 - 0.0037864) pb ~= 0.02956 pb (from PREP)
PATH_1/f11-zeem500to800-powheg-v14b-bp_ntuple.root 0.02956 632 @Powheg 500-800
#
# cross-section of untrimmed m500-inf sample is 0.03335 pb from Prep
#/data/blue/ksung/DYAna2011/s11-zeem500-powheg-v11-pu_ntuple.root 0.03335
#
# cross-section of the 800 sample is unmodified from PREP
PATH_1/f11-zeem800-powheg-v9b-pu_ntuple.root 0.0037864  432 @Powheg 800+
#
# Finally, the 10-20 GeV sample 
# note: instead of xs.dat value 3892.95 we are using
#      powheg_nlo_10-20 * (fewz_nnlo_20-500)/(powheg_nlo_20-500)
#      3216 * 1666/1614 = 3320
PATH_1/f11-zeem1020-powheg-v14b-bp_tight-loose_skim.root 3320  800 @Powheg 10-20
#PATH_1/f11-zeem1020-powheg-v14b-bp_ntuple.root 3320  800 @Powheg 10-20

%
