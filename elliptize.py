#!/usr/bin/python

import sys
import math
from ellipse import ellipse
import numpy as np

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'Usage: elliptize <measurementfile>'
        exit()
    
    with open(sys.argv[1]) as mFile:
        M = dict([(line.split()[0], tuple(eval(f) for f in line.split()[1:6]))
                  for line in mFile.readlines() if '#' not in line])


    with open(sys.argv[1]) as mFile:
        fields = mFile.readline()[1:].split()
        hats = dict([(line.split()[0], dict((field,eval(f)) for field,f in zip(fields[6:],line.split()[6:])))
                     for line in mFile.readlines() if '#' not in line])

    for h,v in hats.items():
        v['f_qq.Ac_y_qq'] = v['fhat_qq'] * v['Ac_y_qq_hat']
        v['f_qg.Ac_y_qq'] = v['fhat_qg'] * v['Ac_y_qg_hat']
    
    form = '\t%.8f'*5
    def deltas(k): return [M['central'][i] - M[k][i] for i in [0,1]]
    def deltasHat(h): return [hats['central'][s] - hats[h][s]
                              for s in ['f_qq.Ac_y_qq','f_qg.Ac_y_qq']]
    
    def distance(k): return math.sqrt( np.dot(*(2*[deltas(k)])) )

    def sigmas(k): return np.outer(*(2*[deltas(k)]))
    def sigmasHat(h): return np.outer(*(2*[deltasHat(h)]))

    def angle(k): return math.atan2(*reversed(deltas(k)))

    print '\n'.join(key.ljust(8) + ": (%.8f)" % distance(key) + form % M[key] 
                    for key in sorted(M, key = distance ))

    
    central = M['central']
    mean = (central[0],central[1])
    stat_sigmas2 = [[central[2],central[3]],
                    [central[3],central[4]]]
    sys_sigmas2 = np.sum(sigmas(k) for k in M) / 2
    sim_sigmas2 = np.sum(sigmasHat(h) for h in hats) / 2
    tot_sigmas2 = stat_sigmas2 + sys_sigmas2
    
    stat = ellipse(mean=mean, sigmas2=stat_sigmas2)
    syst = ellipse(mean=mean, sigmas2=sys_sigmas2)
    totl = ellipse(mean=mean, sigmas2=tot_sigmas2)
    simu = ellipse(mean=(hats['central']['f_qq.Ac_y_qq'],
                         hats['central']['f_qg.Ac_y_qq']),
                   sigmas2=sim_sigmas2) if np.all(sim_sigmas2) else None

    Ac = sum(mean)
    c45 = math.cos(math.pi/4)
    R45 = np.array([[c45,c45],[-c45,c45]])
    print Ac
    print 'xsig', math.sqrt(stat_sigmas2[0][0])
    print 'ysig', math.sqrt(stat_sigmas2[1][1])
    d_Ac = math.sqrt(R45.dot(tot_sigmas2).dot(R45.T)[0,0])
    print d_Ac

    class line(object):
        def __init__(self,start,end):
            self.start = start
            self.end = end
        def eval(self,T):
            t = 0.5*(1+math.cos(T))
            x = t*self.start[0] + (1-t)*self.end[0]
            y = t*self.start[1] + (1-t)*self.end[1]
            return x,y

    with open('_points.'.join(sys.argv[1].split('.')), 'w') as wFile:
        N = 100
        for t in range(N + 1):
            T = t * 2 * math.pi / N
            print >> wFile, '\t'.join(str(f) for f in
                                      sum((getattr(item,'eval')(T)
                                           for item in filter(None, [stat,syst,totl,
                                                                     line(mean,(Ac,0)),
                                                                     line((Ac-d_Ac,0),(Ac+d_Ac,0)),
                                                                     simu
                                                                     ])), ()))
