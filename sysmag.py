#!/usr/bin/python

import sys
import math
from ellipse import ellipse
import numpy as np
from functools import partial
import ROOT as r
import os

class fitresult(object):
    pairs = dict([(item,(item+'_dn',item+'_up')) for item in ['Q','JER','JES','PU','lumi','DY','ST','as']])
    pairs.update(dict([('PD:%02d-%02d'%(i,i+1),('PD_%02d'%i,'PD_%02d'%(i+1))) for i in range(1,53,2)]))
    pairs.update(dict([('muid',('mu2','mu3')),
                       ('mutrig',('mu0','mu1')),
                       ('elid',('el2','el3')),
                       ('eltrig',('el0','el1'))]))
    labels = sum(pairs.values(),())

    def __init__(self,partition):
        self.summarize = summarize
        self.tfile = r.TFile.Open('output/%s/asymmetry_%s_sys.root'%(partition,partition))
        self.cfile = r.TFile.Open('output/%s/asymmetry_%s.root'%(partition,partition))
        self.tree = self.tfile.Get('fitresult')
        self.ctree = self.cfile.Get('fitresult')
        self.ctree.GetEntry(0)
        self.extract()

    def extract(self):
        self.values = {}
        for e in self.tree:
            label = max(os.path.commonprefix([e.label,l]) for l in self.labels)
            if label not in self.labels: continue
            dX = e.fitX - self.ctree.fitX
            dY = e.fitY - self.ctree.fitY
            self.values[label] = {'dX':dX,
                                  'dY':dY,
                                  'dA':dX+dY,
                                  'adA':abs(dX+dY),
                                  'd':math.sqrt(dX**2+dY**2)
                              }
        self.order = sorted(self.values, key = lambda k: self.values[k]['d'], reverse=True)
        print self.order[:5]
        
        self.pvalues = dict([(name, math.sqrt(self.values[a]['dA']**2 + self.values[b]['dA']**2)) for name,(a,b) in self.pairs.items()])
        self.porder = sorted(self.pvalues, key = lambda k: self.pvalues[k], reverse=True)
        print self.porder[:5]
        print

    def form(self,key):
        top = key in self.order[:5]+self.porder[:5]
        form = (r"\textbf{% .3f}" if top else "% .3f ").rjust(20)
        d = (self.values[key]['d'] if key in self.order else self.pvalues[key])*100
        return (form%d).ljust(20)
        

if __name__ == '__main__':

    summarize=True
    partitions = [('full','Full Selection'),
                  ('loM',r'$m_{\ttbar}<450\GeV$'),
                  ('hiM',r'$m_{\ttbar}>450\GeV$'),
                  ('loY',r'$\tanh\abs{y_{\ttbar}}<0.5$'),
                  ('hiY',r'$\tanh\abs{y_{\ttbar}}>0.5$')]
    results = [fitresult(p) for p,_ in partitions]

    def formkey(key) :
        subs = [('PD_',r'pdf'),
                ('mu0',r'$\mu$ trig_up'),
                ('mu1',r'$\mu$ trig_dn'),
                ('mu2',r'$\mu$ id_up'),
                ('mu3',r'$\mu$ id_dn'),
                ('el0',r'$e$ trig_up'),
                ('el1',r'$e$ trig_dn'),
                ('el2',r'$e$ id_up'),
                ('el3',r'$e$ id_dn'),
                ('_up',r'${}^{\uparrow}$'),
                ('_dn',r'${}_{\downarrow}$'),
                ('PU',r'pileup'),
                ('RF',r'RFS'),
                ('lumi',r'\lumi'),
                ('ST',r'single'),
                ('mutrig',r'$\mu$ trig'),
                ('muid', r'$\mu$ id'),
                ('eltrig',r'$e$ trig'),
                ('elid',r'$e$ id'),
                ('PD:',r'pdf '),
                ('as',r'\alpha_s')
                 ]
        def rep(key,ps):
            if not ps: return key
            return rep(key.replace(*ps[0]), ps[1:])
        return rep(key,subs)

    vspace=r'''

&&&&&\\

'''
    caption = r'''
\caption{\label{list_systematics} %s 
due to sources of systematic variations, ordered by decreasing
magnitude in the full selection.
%s
The five greatest sources of systematic uncertainty in each selection are in bold.}
'''%('Uncertainty on $A_c^y$' if summarize else 'Magnitude of the measurement displacement in the plane $(A_c^{y(\QQ)},A_c^{y(\QG)})$',
     '%')


    print r'\begin{longtable}{lccccc}'
    print r'\hline'
    print r'&\multicolumn{5}{c}{(\%)}\\'
    print '  &  '.join(zip(*partitions)[1]).join(['  &  ',r'  \\'])
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
        formkey(key).ljust(28) +' & '+ ' & '.join(fr.form(key) for fr in results) + r'  \\' + (vspace if i%5==4 else '')
        for i,key in enumerate(results[0].porder))
    print r'\end{longtable}'
