import ROOT as r
import numpy as np

files = ["output/full/ens/asymmetry_full%s.root"%s for s in ['A4','A2','A','2A','4A','8A']]

points = {}
for f in files:
    tfile = r.TFile.Open(f)
    tree = tfile.Get('fitresult')
    sigmaX = sum(e.sigmaX for e in tree) / tree.GetEntries()
    sigmaY = sum(e.sigmaY for e in tree) / tree.GetEntries()
    lumiF = sum(e.lumi_factor for e in tree) / tree.GetEntries()
    points[lumiF] = np.array([sigmaX,sigmaY])

for p in points:
    if p!=1: points[p] /= points[1.0]
points[1] /= points[1]

with open('output/lumi_dep.txt','w') as f:
    for k,v in points.items():
        print>> f, str(k).ljust(5), v[0], v[1]
