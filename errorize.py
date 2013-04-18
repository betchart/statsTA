#!/usr/bin/python

import sys
import math
from ellipse import ellipse
import numpy as np

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'Usage: errorize <measurementfile>'
        exit()
    
    with open(sys.argv[1]) as mFile:
        M = dict([(line.split()[0], tuple(eval(f) for f in line.split()[1:3]))
                             for line in mFile.readlines() if '#' not in line])
    
    form = '\t%.8f'*2
    def delta(k): return M['central'][0] - M[k][0]
    
    def distance(k): return abs(delta(k))

    def sigma2(k): return delta(k)**2

    print '\n'.join(key.ljust(8) + ": (%.8f)" % distance(key) + form % M[key] 
                    for key in sorted(M, key = distance ))

    
    central = M['central']
    mean = central[0]
    stat_sigma2 = central[1]**2
    sys_sigma2 = sum(sigma2(k) for k in M) / 2
    
    stat = math.sqrt(stat_sigma2)
    syst = math.sqrt(sys_sigma2)
    totl = math.sqrt(stat_sigma2 + sys_sigma2)

    print "%.8f +/- %.8f stat +/- %.8f sys" % (mean, stat, syst)
