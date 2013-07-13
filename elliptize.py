#!/usr/bin/python

import sys
import math
from ellipse import ellipse
import numpy as np
from CLprojection import oneSigmaCLprojection

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'Usage: elliptize <measurementfile>'
        exit()
    
    with open(sys.argv[1]) as mFile:
        M = dict([(line.split()[0], tuple(eval(f) for f in line.split()[1:6]))
                  for line in mFile.readlines() if '#' not in line])


    with open(sys.argv[1]) as mFile:
        fields = mFile.readline()[1:].split()
        print fields[6:13]
        hats = dict([(line.split()[0], dict((field,eval(f)) for field,f in zip(fields[6:13],line.split()[6:13])))
                     for line in mFile.readlines() if '#' not in line])

    for h,v in hats.items():
        v['f_qq.Ac_y_qq'] = v['fhat_qq'] * v['Ac_y_qq_hat']
        v['f_qg.Ac_y_qg'] = v['fhat_qg'] * v['Ac_y_qg_hat']
    
    form = '\t%.8f'*5
    def deltas(k): return [M['central'][i] - M[k][i] for i in [0,1]]
    def deltasHat(h): return [hats['central'][s] - hats[h][s]
                              for s in ['f_qq.Ac_y_qq','f_qg.Ac_y_qg']]
    
    def distance(k): return math.sqrt( np.dot(*(2*[deltas(k)])) )

    def sigmas(k): return np.outer(*(2*[deltas(k)]))
    def sigmasHat(h): return np.outer(*(2*[deltasHat(h)]))

    def angle(k): return math.atan2(*reversed(deltas(k)))

    print '\n'.join(key.ljust(8) + ": (%.8f)" % distance(key) + form % M[key] 
                    for key in sorted(M, key = distance ))

    
    central = M['central']
    mean = (central[0],central[1])
    simmean = (hats['central']['f_qq.Ac_y_qq'],
               hats['central']['f_qg.Ac_y_qg'])
    stat_sigmas2 = [[central[2],central[3]],
                    [central[3],central[4]]]
    sys_sigmas2 = np.sum(sigmas(k) for k in M) / 2
    sim_sigmas2 = np.sum(sigmasHat(h) for h in hats) / 2
    if not np.all(sim_sigmas2): sim_sigmas2 = np.array([[1e-7,1e-9],[1e-9,1e-7]])
    tot_sigmas2 = stat_sigmas2 + sys_sigmas2

    stat = ellipse(mean=mean, sigmas2=stat_sigmas2)
    syst = ellipse(mean=mean, sigmas2=sys_sigmas2)
    totl = ellipse(mean=mean, sigmas2=tot_sigmas2)
    simu = ellipse(mean=simmean,
                   sigmas2=sim_sigmas2)

    Ac = sum(mean)
    correction = hats['central']['f_gg.Ac_y_gg']
    R = np.array([[1,1],[-1,1]]) # rotate pi/4, except also scale by sqrt(2)
    d_Aqq = oneSigmaCLprojection(tot_sigmas2)
    d_Aqg = oneSigmaCLprojection(tot_sigmas2[(1,0),][:,(1,0)])
    d_Ac = oneSigmaCLprojection(R.dot(tot_sigmas2.dot(R.T)))
    sim_d_Aqq = oneSigmaCLprojection(sim_sigmas2)
    sim_d_Aqg = oneSigmaCLprojection(sim_sigmas2[(1,0),][:,(1,0)])
    sim_d_Ac = oneSigmaCLprojection(R.dot(sim_sigmas2.dot(R.T)))

    print 'xstat', oneSigmaCLprojection(stat_sigmas2)
    print 'ystat', oneSigmaCLprojection(np.array(stat_sigmas2)[(1,0),][:,(1,0)])
    print 'Ac, d_Ac:', Ac, d_Ac
    print 'xsig', d_Aqq
    print 'ysig', d_Aqg
    print 'mean', mean
    print correction

    temp = r'''
            \begin{tabular}{c|r@{.}lr@{.}lr@{.}l}
              \hline
              &\multicolumn{6}{c}{(%s)}\\
              & \multicolumn{2}{c}{$A_c^{y(\QQ)}$} & \multicolumn{2}{c}{$A_c^{y(\QG)}$} & \multicolumn{2}{c}{$A_c^{y}$} \\
              \hline
              \hline
              \POWHEG \sc{ct10} & %s & %s & %s \\
                                & %s & %s & %s \\
              \hline
            \end{tabular}'''
    def form(n,e):
        n*=100
        e*=100
        err_digit = int(math.floor(math.log(abs(e))/math.log(10)))-1 if e else 0
        scale = float(10**(-err_digit)) if e else 1
        p,d = divmod(round(abs(n*scale))/scale,1)
        
        return (("%s%d&%0"+str(-err_digit)+"d(%d)")%('-' if n<0 else ' ',p,d*scale,round(e*scale))).ljust(8)


    print temp%(r'\%',
                form(simmean[0], sim_d_Aqq),
                form(simmean[1], sim_d_Aqg),
                form(sum(simmean)+correction, sim_d_Ac),
                form(mean[0], d_Aqq),
                form(mean[1], d_Aqg),
                form(sum(mean)+correction, d_Ac),
                )


    class line(object):
        def __init__(self,start,end):
            self.start = start
            self.end = end
        def eval(self,T):
            t = 0.5*(1+math.cos(T))
            x = t*self.start[0] + (1-t)*self.end[0]
            y = t*self.start[1] + (1-t)*self.end[1]
            return x,y

    class lineUD(object):
        def __init__(self,mean,sig):
            self.mean = mean
            self.sig = sig
        def eval(self,T):
            t = 0.5*(1+math.cos(T))
            x = self.mean + self.sig*(2*t-1)
            return x,

    with open('_points.'.join(sys.argv[1].split('.')), 'w') as wFile:
        N = 100
        for t in range(N + 1):
            T = t * 2 * math.pi / N
            print >> wFile, '\t'.join(str(f) for f in
                                      sum((getattr(item,'eval')(T)
                                           for item in filter(None, [stat,syst,totl,
                                                                     simu,
                                                                     lineUD(mean[0],d_Aqq),
                                                                     lineUD(mean[1],d_Aqg),
                                                                     lineUD(Ac,d_Ac),
                                                                     lineUD(simmean[0], sim_d_Aqq),
                                                                     lineUD(simmean[1], sim_d_Aqg),
                                                                     lineUD(sum(simmean), sim_d_Ac)
                                                                     ])), ()))
