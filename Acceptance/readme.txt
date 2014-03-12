Notes from 2014 March 11

Main code: plotDYAcceptance.C
     The name is kept, since it was used in DrellYanDMDY package

Systematics studies: 
1) run_accSystematics.C
     The code has built-in configuration file name
     Calls plotDYAcceptance.C several times
     Calculation is done by calling functions from Include/calcCorrectionSyst.h


2) plot_accSystematics.C
     Can be called after run_accSystematics.C
     Calls Efficiency/plot_Systematics.C


plot_Systematics.C combines histograms and prepares plots
     Explicit root-file name is built-in
