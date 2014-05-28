#!/usr/bin/env python

from subprocess import Popen

#data vector 1
name_dv1 = 'ee_results.dat'
#covariance matrix 1
name_cm1 = 'covariance_ee.dat'
#size
length1 = '40'
#data vector 2
name_dv2 = 'mumu_results.dat'
#covariance matrix 2
name_cm2 = 'covariance_mumu.dat'
#length2
length2 = '40' 

#analyzers
Popen('root -b -l -q \'resultCombiner.C(\"'+name_dv1+'\",\"'+name_cm1+'\",\"'+length1+'\",\"'+name_dv2+'\",\"'+name_cm2+'\",\"'+length2+'\")\'',shell=True).wait()
