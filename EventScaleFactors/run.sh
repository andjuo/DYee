#!/bin/bash

#rootb calcEventEff_new.C+\(0,\"default\",1,DYTools::NORMAL_RUN,DYTools::FSR_5plus\) | tee log-calcEventEff-step1-20140407-FSR105.out
#rootb calcEventEff_new.C+\(0,\"defaultEgamma\",0,DYTools::NORMAL_RUN,DYTools::NO_SYST\) | tee log-calcEventEff-step2-20140407.out
#rootb calcEventEff_new.C+\(0,\"default\",1,DYTools::NORMAL_RUN,DYTools::FSR_5minus\) | tee log-calcEventEff-step1-20140407-FSR095.out
rootb calcEventEff_new.C+\(0,\"default\",1,DYTools::NORMAL_RUN,DYTools::UNREGRESSED_ENERGY\) | tee log-calcEventEff-step1-20140407-UnRegEn.out
#rootb calcEventEff_new.C+\(0,\"default\",1,DYTools::NORMAL_RUN,DYTools::PILEUP_5plus\) | tee log-calcEventEff-step1-20140407-PU5plus.out
#rootb calcEventEff_new.C+\(0,\"default\",1,DYTools::NORMAL_RUN,DYTools::PILEUP_5minus\) | tee log-calcEventEff-step1-20140407-PU5minus.out
rootb calcEventEff_new.C+\(0,\"default\",1,DYTools::NORMAL_RUN,DYTools::UNREG_PU5plus\) | tee log-calcEventEff-step1-20140407-UnRegPU5plus.out
rootb calcEventEff_new.C+\(0,\"default\",1,DYTools::NORMAL_RUN,DYTools::UNREG_PU5minus\) | tee log-calcEventEff-step1-20140407-UnRegPU5minus.out
