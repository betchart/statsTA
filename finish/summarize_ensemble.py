#!/usr/bin/env python

import sys
import math
from lib.ellipse import ellipse
import numpy as np
from lib.__autoBook__ import autoBook
import ROOT as r
r.gROOT.SetBatch(1)
r.gStyle.SetOptFit(1)

oneSigmaN2LL = 1.14 * 2

if __name__ == '__main__':
    scale = 100
    if len(sys.argv) < 2:
        print 'Usage: summarize_ensemble.py <ensemblefile>'
        exit()
    
    infile  = r.TFile.Open(sys.argv[1])
    tree = infile.Get('fitresult')
    outname = sys.argv[1].replace('asymmetry_','ensemble').replace('.root', '')
    if outname in sys.argv[1]:
        print "Please rename input file."
        exit()

    mark = len(sys.argv) > 2
    if mark:
        mfile = r.TFile.Open(sys.argv[2])
        mtree = infile.Get('fitresult')
        mtree.GetEntry(0)

    tfile = r.TFile.Open(outname+'.root', 'RECREATE')
    book = autoBook(tfile)
    names = ['delta_Aqq','delta_Aqg','delta_A',
             'error_Aqq','error_Aqg','error_A',
             'pullqq',    'pullqg',  'pull']
    fixedLimits = [(-1.5,1.5),(-1.5,1.5),(-5,5),
                   (0.05,1.05),(0.05,1.05),(0.05,1.05),
                   (-5,5),(-5,5),(-5,5)]
    meanNLL = sum(e.NLL for e in tree) / tree.GetEntries()
    limits = fixedLimits
    wNLL = 40000

    within=0
    tot = 0
    for e in tree:
        truth = e.gen_fitX,e.gen_fitY,1.
        mean = e.fitX,e.fitY
        sigmas2 = [[e.fitXX,e.fitXY],
                   [e.fitXY,e.fitYY]]
        el = ellipse(mean=mean,sigmas2=sigmas2)
        if np.dot(truth, el.matrix).dot(truth) < 0: within+=1
        tot+=1
        nbins = 30
        book.fill(meanNLL, 'meanNLL', nbins, -6.2e6,-5.8e6)
        book.fill(e.NLL-meanNLL, 'dNLL', nbins, -wNLL,wNLL)
        try:
            R = np.array([[1,1],[-1,1]]) # rotate pi/4, except also scale by sqrt(2)
            sigmas = np.array([[e.fitXX,e.fitXY],[e.fitXY,e.fitYY]]) / oneSigmaN2LL
            sigma = math.sqrt(R.dot(sigmas).dot(R.T)[0,0])

            values = 100*np.array([e.fitX-e.gen_fitX, e.fitY-e.gen_fitY, e.fitX+e.fitY+e.correction - e.gen_fitX - e.gen_fitY - e.gen_correction,
                                   e.sigmaX, e.sigmaY, sigma ])
            values = tuple(values) + (values[0]/values[3], values[1]/values[4], values[2]/values[5])
            for n,v,lim in zip(names,values,limits):
                book.fill(v,n,nbins,*lim)
            book.fill(tuple(values[2:4]), 'errors2D', (nbins,nbins), (limits[2][0],limits[3][0]), (limits[2][1],limits[3][1]))
        except:
            pass
    print "%d / %d  :  %.3f" % (within, tot, within / float(tot))

    c = r.TCanvas()
    c.Divide(3,3)
    lines = []
    for i,k in enumerate(names):
        c.cd(i+1)
        book[k].Fit('gaus','Q')
        book[k].Draw()
        if mark and 1 < i < 4:
            x = 100 * (mtree.sigmaX if i==2 else mtree.sigmaY)
            hi = book[k].GetMaximum()
            lines.append(r.TLine(x,0,x,hi))
            lines[-1].SetLineColor(r.kGray)
            lines[-1].Draw()

    c.Print(outname+'.pdf')
    tfile.Write()
    tfile.Close()
