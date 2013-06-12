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
    if len(sys.argv) < 2:
        print 'Usage: summarize_ensemble.py <ensemblefile>'
        exit()
    
    with open(sys.argv[1]) as mFile:
        M = dict([(line.split()[0], tuple(eval(f) for f in line.split()[1:6]))
                  for line in mFile.readlines() if '#' not in line])

    with open(sys.argv[1]) as mFile:
        fA = [eval(i) for i in mFile.readline().split()[1:]]

    def deltas(k): return [ M[k][i] - fA[i] for i in [0,1]]
    def pulls(k): return [d/math.sqrt(M[k][i]) for d,i in zip(deltas(k),[2,4])]

    book = autoBook('book')
    names = ['efAqq','efAqg','pullqq','pullqg']
    limits = [(0,0.01),(0,0.003),(-3,3),(-3,3)]
    truth = tuple(fA + [1.])
    count=0
    for k in M:
        values = [abs(f) for f in deltas(k)] + pulls(k) 
        for n,v,lim in zip(names,values,limits):
            book.fill(v,n,30,*lim)
        mean = M[k][:2]
        sigmas2 = [[M[k][2], M[k][3]],
                   [M[k][3], M[k][4]]]
        el = ellipse(mean=mean,sigmas2=sigmas2)
        level = np.dot(truth, el.matrix).dot(truth)
        print '\t'.join(['%+.5f'%i for i in values+[level]])
        if level<0: count+=1
    print count

    c = r.TCanvas()
    c.Divide(2,2)
    for i,k in enumerate(names):
        c.cd(i+1)
        book[k].Fit('gaus')
        book[k].Draw()

    c.Print('.'.join(sys.argv[1].split('.')[:-1]+['pdf']))
