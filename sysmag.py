#!/usr/bin/python

import sys
import math
from ellipse import ellipse
import numpy as np
from functools import partial

if __name__ == '__main__':

    partitions = ['full','hiM','loM','hiY','loY']
    files = ['data/asymmetry_%s.txt'%s for s in  partitions]
    Ms = []

    for f in files:
        with open(f) as mFile:
            M = dict([(line.split()[0], tuple(eval(f) for f in line.split()[1:6]))
                      for line in mFile.readlines() if '#' not in line])
            Ms.append(M)

    
    form = '\t%.8f'*5
    def deltas(M,k): return [M['central'][i] - M[k][i] for i in [0,1]]
    
    def distance(M,k): return math.sqrt( np.dot(*(2*[deltas(M,k)])) )

    def angle(M,k): return math.atan2(*reversed(deltas(M,k)))

    maxd = [max(distance(M,k) for k in M) for M in Ms]
    def order(k):
        return sum(distance(M,k)/m for M,m in zip(Ms,maxd))

    keys = sorted( set(sum([sorted(M, key=partial(distance,M)) for M in Ms],[])),
                   key=partial(distance,Ms[0]), reverse=True)

    def form(M,key):
        top = sorted(M, key=partial(distance,M))[-5:]
        form = (r"\textbf{%.3f}" if key in top else "%.3f ").rjust(20)
        d = distance(M,key)*100
        return form%d

    headers = {'full':'Full Selection',
               'hiM':r'$m_{\ttbar}>450\GeV$',
               'loM':r'$m_{\ttbar}<450\GeV$',
               'hiY':r'$\tanh\abs{y_{\ttbar}}>0.5$',
               'loY':r'$\tanh\abs{y_{\ttbar}}<0.5$',
               }

    def formkey(key) :
        pairs = [('PD_','pdf'),
                 ('mu0','$\mu$ trig_up'),
                 ('mu1','$\mu$ trig_dn'),
                 ('mu2','$\mu$ id_up'),
                 ('mu3','$\mu$ id_dn'),
                 ('el0','$e$ trig_up'),
                 ('el1','$e$ trig_dn'),
                 ('el2','$e$ id_up'),
                 ('el3','$e$ id_dn'),
                 ('_up','${}^{\uparrow}$'),
                 ('_dn','${}_{\downarrow}$'),
                 ('PU','pileup'),
                 ('RF','RFS'),
                 ('lumi','\lumi'),
                 ('ST','single')
                 ]
        def rep(key,ps):
            if not ps: return key
            return rep(key.replace(*ps[0]), ps[1:])
        return rep(key,pairs)

    vspace=r'''

&&&&&\\

'''
    caption = r'''
\caption{\label{list_systematics} Magnitude of the measurement
         displacement in the plane $(A_c^{y(\QQ)},A_c^{y(\QG)})$ due to
         systematic variations, ordered by decreasing magnitude in the full selection.
%
The five greatest sources of systematic uncertainty in each selection are in bold.}
'''


    print r'\begin{longtable}{lccccc}'
    print r'\hline'
    print r'&\multicolumn{5}{c}{(\%)}'
    print '  &  '.join(headers[p] for p in partitions).join(['  &  ',r'  \\'])
    print r'\hline'
    print r'\hline'
    print r'\endhead'
    print r'\hline'
    print caption
    print r'\endfoot'
    print r'\hline'
    print r'\caption[](continued)'
    print r'\endlastfoot'
    print '\n'.join(
        formkey(key).ljust(15) +' & '+ ' & '.join(form(M,key) for M in Ms) + r'  \\' + (vspace if i%5==4 else '')
        for i,key in enumerate(keys))
    print r'\end{longtable}'
