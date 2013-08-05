import ROOT as r
import numpy as np
import math
from ellipse import ellipse
from CLprojection import oneSigmaCLprojection

class fitresult(object):
    def __init__(self,fname, getfirstentry = False):
        self.tfile = r.TFile.Open(fname)
        self.tree = self.tfile.Get('fitresult')
        if getfirstentry: self.tree.GetEntry(0)

class lineUD(object):
    def __init__(self,mean,sig):
        self.mean = mean
        self.sig = sig
    def eval(self,T):
        t = 0.5*(1+math.cos(T))
        x = self.mean + self.sig*(2*t-1)
        return x,

class collate(object):
    def __init__(self,partition):
        self.partition = partition
        self.central = fitresult('output/%s/asymmetry_%s.root'%(partition,partition), True)
        self.thr30 = fitresult('output/%s/asymmetry_%s_sys_thr30.root'%(partition,partition), True)
        self.sys = fitresult('output/%s/asymmetry_%s_sys.root'%(partition,partition))
        self.pois = fitresult('output/%s/asymmetry_%s_t.root'%(partition,partition))

        self.mean = np.array([self.central.tree.fitX, self.central.tree.fitY])
        self.thr30_mean = np.array([self.thr30.tree.fitX, self.thr30.tree.fitY])
        self.simmean = self.sim(self.central.tree)
        self.simmean30 = self.sim(self.thr30.tree)

        self.sigmas_stat = np.array([[self.central.tree.fitXX,self.central.tree.fitXY],[self.central.tree.fitXY,self.central.tree.fitYY]])
        self.sigmas_syst = sum([self.sigmas(self.deltas(e)) for e in self.sys.tree],self.sigmas([0,0])) / 2
        self.sigmas_pois = sum([self.sigmas(self.deltas(e)) for e in self.pois.tree],self.sigmas([0,0])) / self.pois.tree.GetEntries()
        self.sigmas_simu = sum([self.sigmas(self.deltas_sim(e)) for e in self.sys.tree],self.sigmas([0,0]))

        self.sigmas_totl = sum([self.sigmas_stat,self.sigmas_syst,self.sigmas_pois],self.sigmas([0,0]))
        self.write('output/asymmetry_%s_points.txt'%partition)
        print self

    def __str__(self):
        return '\n'.join([str(i) for i in [
            '%s % .4f  % .4f  % .4f'%(self.partition.ljust(5), 100*self.mean[0], 100*self.mean[1], 100*sum(self.mean)),
            self.sigmas_stat,
            self.sigmas_pois,
            self.sigmas_syst]])

    def sim(self,e): return np.array([e.Ac_y_ttqq * e.f_qq_hat, e.Ac_y_ttqg * e.f_qg_hat])

    def deltas(self,e): return (e.fitX - self.mean[0], e.fitY - self.mean[1])
    def deltas_sim(self,e): return self.sim(e) - self.simmean

    @staticmethod
    def sigmas(deltas): return np.outer(*(2*[deltas]))


    def write(self,fname):
        R = np.array([[1,1],[-1,1]]) # rotate pi/4, except also scale by sqrt(2)
        d_Aqq = oneSigmaCLprojection(self.sigmas_totl)
        d_Aqg = oneSigmaCLprojection(self.sigmas_totl[(1,0),][:,(1,0)])
        d_Ac = oneSigmaCLprojection(R.dot(self.sigmas_totl.dot(R.T)))
        sim_d_Aqq = oneSigmaCLprojection(self.sigmas_simu)
        sim_d_Aqg = oneSigmaCLprojection(self.sigmas_simu[(1,0),][:,(1,0)])
        sim_d_Ac = oneSigmaCLprojection(R.dot(self.sigmas_simu.dot(R.T)))

        stat = ellipse(mean=list(self.mean), sigmas2=list(self.sigmas_stat))
        syst = ellipse(mean=list(self.mean), sigmas2=list(self.sigmas_syst))
        pois = ellipse(mean=list(self.mean), sigmas2=list(self.sigmas_pois))
        totl = ellipse(mean=list(self.mean), sigmas2=list(self.sigmas_totl))
        simu = ellipse(mean=list(self.simmean), sigmas2=list(self.sigmas_simu))
        poissyst = ellipse(mean=list(self.mean), sigmas2=list(self.sigmas_syst+self.sigmas_pois))

        with open(fname.replace('_points',''), 'w') as wFile:
            print >> wFile, 'central', self.mean[0], self.mean[1]
            print >> wFile, 'thr30',  self.thr30_mean[0], self.thr30_mean[1]
            print >> wFile, 'sim', self.simmean[0], self.simmean[1]
            print >> wFile, 'sim30', self.simmean30[0], self.simmean30[1]

        with open(fname, 'w') as wFile:
            N = 100
            for t in range(N + 1):
                T = t * 2 * math.pi / N
                print >> wFile, '\t'.join(str(f) for f in
                                          sum((item.eval(T)
                                               for item in filter(None, [stat,syst,totl,
                                                                         simu,
                                                                         lineUD(self.mean[0],d_Aqq),
                                                                         lineUD(self.mean[1],d_Aqg),
                                                                         lineUD(sum(self.mean),d_Ac),
                                                                         lineUD(self.simmean[0], sim_d_Aqq),
                                                                         lineUD(self.simmean[1], sim_d_Aqg),
                                                                         lineUD(sum(self.simmean), sim_d_Ac),
                                                                         pois,
                                                                         poissyst
                                                                     ])), ()))


if __name__=='__main__':
    [collate(p) for p in ['full','hiM','loM','hiY','loY']]
