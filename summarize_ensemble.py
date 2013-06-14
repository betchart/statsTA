#!/usr/bin/python

import sys
import math
from ellipse import ellipse
from CLprojection import oneSigmaCLprojection
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
    
    with open(sys.argv[1]) as mFile:
        M = dict([(line.split()[0], tuple(eval(f) for f in line.split()[1:6]))
                  for line in mFile.readlines() if '#' not in line])

    with open(sys.argv[1]) as mFile:
        fA = [eval(i) for i in mFile.readline().split()[1:]]

    def deltas(k): return [ scale*(M[k][i] - fA[i]) for i in [0,1]]
    def errors(k):
        sigmas2 = np.array([[M[k][2],M[k][3]],
                            [M[k][3],M[k][4]]])
        return [scale * oneSigmaCLprojection(sigmas2),
                scale * oneSigmaCLprojection(sigmas2[(1,0),][:,(1,0)])]
    def pulls(k): return [d/e for d,e in zip(deltas(k),errors(k))]

    tfile = r.TFile.Open('.'.join(sys.argv[1].split('.')[:-1]+['root']), 'RECREATE')
    book = autoBook(tfile)
    names = ['delta_Aqq','delta_Aqg','error_Aqq','error_Aqg','pullqq','pullqg']
    limits = [(-1.5,1.5),(-1.5,1.5),(0.34,0.4),(0.09,0.15),(-5,5),(-5,5)]
    truth = tuple(fA + [1.])
    within=0
    for k in M:
        mean = M[k][:2]
        sigmas2 = [[M[k][2], M[k][3]],
                   [M[k][3], M[k][4]]]
        el = ellipse(mean=mean,sigmas2=sigmas2)
        if np.dot(truth, el.matrix).dot(truth) < 0: within+=1
        try:
            values = deltas(k) + errors(k) + pulls(k)
            for n,v,lim in zip(names,values,limits):
                book.fill(v,n,40,*lim)
        except:
            pass
    print within

    c = r.TCanvas()
    c.Divide(2,3)
    for i,k in enumerate(names):
        c.cd(i+1)
        book[k].Fit('gaus','Q')
        book[k].Draw()

    c.Print('.'.join(sys.argv[1].split('.')[:-1]+['pdf']))
    tfile.Write()
    tfile.Close()
