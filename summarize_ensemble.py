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
    if len(sys.argv) < 2:
        print 'Usage: summarize_ensemble.py <ensemblefile>'
        exit()
    
    with open(sys.argv[1]) as mFile:
        M = dict([(line.split()[0], tuple(eval(f) for f in line.split()[1:6]))
                  for line in mFile.readlines() if '#' not in line])

    with open(sys.argv[1]) as mFile:
        fA = [eval(i) for i in mFile.readline().split()[1:]]

    def deltas(k): return [ M[k][i] - fA[i] for i in [0,1]]
    def errors(k):
        sigmas2 = np.array([[M[k][2],M[k][3]],
                            [M[k][3],M[k][4]]])
        return [oneSigmaCLprojection(sigmas2),
                oneSigmaCLprojection(sigmas2[(1,0),][:,(1,0)])]
    def pulls(k): return [d/e for d,e in zip(deltas(k),errors(k))]

    book = autoBook('book')
    names = ['delta_Aqq','delta_Aqg','error_Aqq','error_Aqg','pullqq','pullqg']
    limits = [(-0.013,0.013),(-0.004,0.004),(0.0034,0.004),(0.0010,0.0014),(-4,4),(-4,4)]
    truth = tuple(fA + [1.])
    count=0
    for k in M:
        try:
            values = deltas(k) + errors(k) + pulls(k)
            for n,v,lim in zip(names,values,limits):
                book.fill(v,n,40,*lim)
            mean = M[k][:2]
            sigmas2 = [[M[k][2], M[k][3]],
                       [M[k][3], M[k][4]]]
            el = ellipse(mean=mean,sigmas2=sigmas2)
            level = np.dot(truth, el.matrix).dot(truth)
            print '\t'.join(['%+.5f'%i for i in values+[level]])
            if level<0: count+=1
        except:
            pass
    print count

    c = r.TCanvas()
    c.Divide(2,3)
    for i,k in enumerate(names):
        c.cd(i+1)
        book[k].Fit('gaus')
        book[k].Draw()

    c.Print('.'.join(sys.argv[1].split('.')[:-1]+['pdf']))
