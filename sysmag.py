#!/usr/bin/python

import sys
import math
from ellipse import ellipse
import numpy as np
from functools import partial

summarize = True

pairs = dict([(item,(item+'_dn',item+'_up')) for item in ['RF','JER','JES','PU','lumi','DY','ST','DY2','ST2'][:-4]])
pairs.update(dict([('PD:%02d-%02d'%(i,i+1),('PD_%02d'%i,'PD_%02d'%(i+1))) for i in range(1,53,2)]))
pairs.update(dict([('muid',('mu2','mu3')),
                   ('mutrig',('mu0','mu1')),
                   ('elid',('el2','el3')),
                   ('eltrig',('el0','el1'))]))

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
    def diff(M,k): return sum(M['central']) - sum(M[k])
    def pairsig(M,p): return math.sqrt( 0.5* (diff(M,pairs[p][0])**2 + diff(M,pairs[p][1])**2) )

    keys = sorted( pairs if summarize else Ms[0], 
                   key=partial(pairsig if summarize else distance,Ms[0]), reverse=True)

    def form(M,key):
        func = pairsig if summarize else distance
        top = sorted(pairs if summarize else M, key=partial(func,M))[-5:]
        form = (r"\textbf{% .3f}" if key in top else "% .3f ").rjust(20)
        d = func(M,key)*100
        return (form%d).ljust(20)

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
                 ('ST','single'),
                 ('mutrig','$\mu$ trig'),
                 ('muid', '$\mu$ id'),
                 ('eltrig','$e$ trig'),
                 ('elid', '$e$ id'),
                 ('PD:','pdf ')
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
        formkey(key).ljust(28) +' & '+ ' & '.join(form(M,key) for M in Ms) + r'  \\' + (vspace if i%5==4 else '')
        for i,key in enumerate(keys))
    print r'\end{longtable}'
