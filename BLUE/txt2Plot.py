#!/usr/bin/env python

import sys
from ROOT import *
from array import array


binning = [15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,64.,68, 72, 76, 81, 86, 91, 96,
101, 106, 110, 115, 120, 126, 133,
141, 150, 160, 171, 185, 200, 220,
243, 273, 320, 380, 440, 510, 600,1000,1500]

h1 = TH1D('h1','h1',40,array('d',binning))
h2 = TH1D('h2','h2',40,array('d',binning))

#binning file
values = open('real.txt','r')
val = values.readlines()
values.close()

errors = open('andriu.txt','r')
err = errors.readlines()
errors.close()



for i in range(len(val)):
    h1.SetBinContent(i+1,float(val[i]))
    h2.SetBinContent(i+1,float(err[i])) 
    h1.SetBinError(i+1,0.)
    h2.SetBinError(i+1,0.)
 

c = TCanvas()
c.cd()
c.SetLogx()
c.SetLogy()
h1.GetXaxis().SetMoreLogLabels()
h1.GetYaxis().SetRangeUser(0.000000000001,0.03)
h1.SetMarkerStyle(20)
h1.SetMarkerColor(kRed)
h1.SetLineColor(kRed)
h1.GetYaxis().SetTitle("Abs. uncertainty")
h1.GetXaxis().SetTitle("M_{#mu#mu} [GeV]")
h1.Draw("PE")

h2.SetMarkerStyle(20)
h2.SetMarkerColor(kBlue)
h2.SetLineColor(kBlue)
h2.Draw("PEsame")


c.SaveAs("1.png")
