#!/usr/bin/python

import sys
import math
from ellipse import ellipse
import numpy as np
from __autoBook__ import autoBook
import ROOT as r
r.gROOT.SetBatch(1)
r.gStyle.SetOptFit(1)

if __name__ == '__main__':
    scale = 100
    if len(sys.argv) < 2:
        print 'Usage: summarize_ensemble.py <ensemblefile>'
        exit()
    
    infile  = r.TFile.Open(sys.argv[1])
    tree = infile.Get('fitresult')
    outname = sys.argv[1].replace('asymmetry_full','ensemble').replace('.root', '')

    tfile = r.TFile.Open(outname+'.root', 'RECREATE')
    book = autoBook(tfile)
    names = ['delta_Aqq','delta_Aqg','error_Aqq','error_Aqg','pullqq','pullqg']
    fixedLimits = [(-1,1),(-1,1),(0.05,0.45),(0.04,0.14),(-5,5),(-5,5)]
    meanNLL = sum(e.NLL for e in tree) / tree.GetEntries()
    limits = fixedLimits
    wNLL = 40000

    within=0
    for e in tree:
        truth = e.gen_fitX,e.gen_fitY,1.
        mean = e.fitX,e.fitY
        sigmas2 = [[e.fitXX,e.fitXY],
                   [e.fitXY,e.fitYY]]
        el = ellipse(mean=mean,sigmas2=sigmas2)
        if np.dot(truth, el.matrix).dot(truth) < 0: within+=1
        nbins = 30
        book.fill(meanNLL, 'meanNLL', nbins, -6.2e6,-5.8e6)
        book.fill(e.NLL-meanNLL, 'dNLL', nbins, -wNLL,wNLL)
        try:
            values = 100*np.array([e.fitX-e.gen_fitX, e.fitY-e.gen_fitY, e.sigmaX, e.sigmaY])
            values = tuple(values) + (values[0]/values[2], values[1]/values[3])
            for n,v,lim in zip(names,values,limits):
                book.fill(v,n,nbins,*lim)
            book.fill(tuple(values[2:4]), 'errors2D', (nbins,nbins), (limits[2][0],limits[3][0]), (limits[2][1],limits[3][1]))
        except:
            pass
    print within

    c = r.TCanvas()
    c.Divide(2,3)
    for i,k in enumerate(names):
        c.cd(i+1)
        book[k].Fit('gaus','Q')
        book[k].Draw()

    c.Print(outname+'.pdf')
    tfile.Write()
    tfile.Close()
