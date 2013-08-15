import ROOT as r
import numpy as np
import math
from ellipse import ellipse

class fitresult(object):
    def __init__(self,fname, getfirstentry = False):
        self.tfile = r.TFile.Open(fname)
        self.tree = self.tfile.Get('fitresult')
        if getfirstentry: self.tree.GetEntry(0)

class collate(object):
    def __init__(self,partition):
        self.partition = partition
        self.central = fitresult('output/%s/asymmetry_%s.root'%(partition,partition), True)
        self.sys = fitresult('output/%s/asymmetry_%s_sys.root'%(partition,partition))
        self.pois = fitresult('output/%s/asymmetry_%s_t.root'%(partition,partition))
        self.stat = fitresult('output/%s/ens/asymmetry_%sA.root'%(partition,partition))

    def val_err(self, item, scale=1):
        val = getattr(self.central.tree,item) if item not in ['el','mu'] else self.expect(self.central.tree,item)
        def sig2(tree): return sum([(val-(getattr(e,item) if item not in ['el','mu'] else self.expect(e,item)))**2 for e in tree])
        
        stat = sig2(self.stat.tree) / self.stat.tree.GetEntries()
        syst = sig2(self.sys.tree) / 2
        pois = sig2(self.pois.tree) / self.pois.tree.GetEntries()

        return scale*val, scale*math.sqrt(stat+syst+pois)

    def expect(self,tree,lep):
        return sum(getattr(tree,'expect_%s_%s'%(lep,item)) for item in ['tt','wj','mj','st','dy'])

    @property
    def label(self):
        labels = {'full':'Full Selection',
                  'hiM':r'$m_{\ttbar}>450\GeV$',
                  'loM':r'$m_{\ttbar}<450\GeV$',
                  'hiY':r'$\tanh\abs{y_{\ttbar}}>0.5$',
                  'loY':r'$\tanh\abs{y_{\ttbar}}<0.5$',
              }
        return labels[self.partition].ljust(len(labels['loY'])+1)

    def __str__(self):
        return '\n'.join([self.label +' & ' + '  &  '.join('%.2f %.2f'%self.val_err(item,100) for item in ['d_xs_tt','d_xs_wj','f_gg','f_qq','f_qg','f_ag']),
                          self.label +' & ' + '  &  '.join('%.2f %.2f'%self.val_err(item,0.001) for item in ['expect_el_tt','expect_el_wj','expect_el_mj','expect_el_st','expect_el_dy','el']),
                          self.label +' & ' + '  &  '.join('%.2f %.2f'%self.val_err(item,0.001) for item in ['expect_mu_tt','expect_mu_wj','expect_mu_mj','expect_mu_st','expect_mu_dy','mu']),
                      ])

if __name__=='__main__':
    for item in ['full','hiM','loM','hiY','loY']:
        print collate(item)
