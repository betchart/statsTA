import ROOT as r
r.gROOT.SetBatch(1)
from utils.ellipse import ellipse
import array
import math
import numpy as np

tfile = r.TFile.Open("jiggled_asymmetries.root")
ells = {}

for comp in ['ttqq','ttqg','ttag']:
    for lep in ['el','mu']:
        for epoch in ['orig','extra','total']:
            h = tfile.Get('_'.join([epoch,lep,comp]))
            stats = array.array('d',7*[0.])
            h.GetStats(stats)
            dx2 = stats[3]/stats[0] - stats[2]**2/stats[0]**2
            dy2 = stats[5]/stats[0] - stats[4]**2/stats[0]**2
            dxy = stats[6]/stats[0] - stats[2]*stats[4]/stats[0]**2
            sigmas2 = [[dx2, dxy], [dxy, dy2]]
            mean = ( h.GetMean(1), h.GetMean(2) )
            el = ellipse(mean=mean, sigmas2=sigmas2)
            ells[h.GetName()] = el

order = ['_'.join([epoch,lep,'']) for lep in ['el','mu'] for epoch in ['orig','extra','total']]
for comp in ['ttqq','ttqg','ttag']:
    with open(comp+'_points.txt','w') as out:
        meanel = ells['total_el_'+comp].mean
        meanmu = ells['total_mu_'+comp].mean
        print comp, meanel[0], meanel[1], meanmu[0], meanmu[1]
        print>> out, '# ', '\t'.join(order)
        N = 100
        for t in range(N + 1):
            T = t * 2 * math.pi / N
            print >> out, '\t'.join("%f\t%f\t\t"%ells[name].eval(T) for name in [o+comp for o in order])
