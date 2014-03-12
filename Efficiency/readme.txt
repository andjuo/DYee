Notes from 2014 March 11

Main code: plotDYEfficiency.C
     The name is kept, since it was used in DrellYanDMDY package

Systematics studies: 
1) run_effSystematics.C
     The code has built-in configuration file name
     Calls plotDYEfficiency.C several times
     Calculation is done by calling functions from Include/calcCorrectionSyst.h


2) plot_effSystematics.C
     Can be called after run_effSystematics.C
     Explicit root-file name is built-in
     Calls Efficiency/plot_Systematics.C


plot_Systematics.C combines histograms and prepares plots
