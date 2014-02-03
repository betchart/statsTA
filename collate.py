import ROOT as r
import numpy as np
import math
from utils.ellipse import ellipse

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

oneSigmaN2LL = 1.14 * 2

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

        self.sigmas_stat = np.array([[self.central.tree.fitXX,self.central.tree.fitXY],[self.central.tree.fitXY,self.central.tree.fitYY]]) / oneSigmaN2LL
        self.sigmas_syst = sum([self.sigmas(self.deltas(e)) for e in self.sys.tree],self.sigmas([0,0])) / 2
        self.sigmas_pois = sum([self.sigmas(self.deltas(e)) for e in self.pois.tree],self.sigmas([0,0])) / self.pois.tree.GetEntries()
        self.sigmas_simu = sum([self.sigmas(self.deltas_sim(e)) for e in self.sys.tree],self.sigmas([0,0]))

        self.sigmas_totl = sum([self.sigmas_stat,self.sigmas_syst,self.sigmas_pois],self.sigmas([0,0]))
        self.write('output/asymmetry_%s_points.txt'%partition)

    @property
    def label(self):
        return {'full':'Full Selection',
                'hiM':r'$m_{\ttbar}>450\GeV$',
                'loM':r'$m_{\ttbar}<450\GeV$',
                'hiY':r'$\tanh\abs{y_{\ttbar}}>0.5$',
                'loY':r'$\tanh\abs{y_{\ttbar}}<0.5$',
        }[self.partition]

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
        d_Aqq,d_Aqg = np.sqrt(np.diag(self.sigmas_totl))
        d_Ac = math.sqrt(R.dot(self.sigmas_totl.dot(R.T))[0,0])
        sim_d_Aqq,sim_d_Aqg = np.sqrt(np.diag(self.sigmas_simu))
        sim_d_Ac = math.sqrt(R.dot(self.sigmas_simu.dot(R.T))[0,0])
        syst_d_Ac = math.sqrt(R.dot((self.sigmas_syst+self.sigmas_pois).dot(R.T))[0,0])
        stat_d_Ac = math.sqrt(R.dot(self.sigmas_stat).dot(R.T)[0,0])
        print 100*stat_d_Ac, 100*syst_d_Ac

        stat = ellipse(mean=list(self.mean), sigmas2=list(oneSigmaN2LL * self.sigmas_stat))
        syst = ellipse(mean=list(self.mean), sigmas2=list(oneSigmaN2LL * self.sigmas_syst))
        pois = ellipse(mean=list(self.mean), sigmas2=list(oneSigmaN2LL * self.sigmas_pois))
        totl = ellipse(mean=list(self.mean), sigmas2=list(oneSigmaN2LL * self.sigmas_totl))
        simu = ellipse(mean=list(self.simmean), sigmas2=list(oneSigmaN2LL * self.sigmas_simu))
        poissyst = ellipse(mean=list(self.mean), sigmas2=list(oneSigmaN2LL * (self.sigmas_syst+self.sigmas_pois)))

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
        temp = r'''
        \begin{tabular}{c|r@{.}lr@{.}lr@{.}l}
        \hline
        &\multicolumn{6}{c}{(%s)}\\
        & \multicolumn{2}{c}{$A_c^{y(\QQ)}$} & \multicolumn{2}{c}{$A_c^{y(\QG)}$} & \multicolumn{2}{c}{$A_c^{y}$} \\
        \hline
        \hline
        \POWHEG \sc{ct10} & %s & %s & %s \\
        %s & %s & %s & %s \\
        \hline
        \end{tabular}'''
        def form(n,e):
            n*=100
            e*=100
            err_digit = int(math.floor(math.log(abs(e))/math.log(10)))-1 if e else 0
            scale = float(10**(-err_digit)) if e else 1
            p,d = divmod(round(abs(n*scale))/scale,1)

            return (("%s%d&%0"+str(-err_digit)+"d(%d)")%('-' if n<0 else ' ',p,d*scale,round(e*scale))).ljust(8)

        simcorrection = self.central.tree.f_gg_hat * self.central.tree.Ac_y_ttgg

        print temp%(r'\%',
                    form(self.simmean[0], sim_d_Aqq),
                    form(self.simmean[1], sim_d_Aqg),
                    form(sum(self.simmean)+simcorrection, sim_d_Ac),
                    self.label,
                    form(self.mean[0], d_Aqq),
                    form(self.mean[1], d_Aqg),
                    form(sum(self.mean)+self.central.tree.correction, d_Ac),
                    )



if __name__=='__main__':
    [collate(p) for p in ['full','loM','hiM','loY','hiY']]
