#!/usr/bin/python

import sys
import ROOT as r
import systematics
from measurement import measurement
from options import opts

if __name__ == '__main__':
    r.gROOT.SetBatch(1)
    options = opts()

    measure = 0
    if options.partitions==True:
        print ','.join(systematics.partitions)
        exit()
    partitions = options.partitions.split(',')

    allSystematics = [s['label'] for s in systematics.systematics()]
    if options.systematics == True:
        print ','.join(allSystematics)
        exit()
    systs = ([] if not options.systematics else 
             allSystematics if options.systematics=='all' else 
             options.systematics.split(','))

    ensembles = (None if not options.ensembles else tuple(int(i) for i in options.ensembles.split(':')))
    templates = (None if not options.templates else tuple(options.ensembles))
    
    chunk = 10

    if ( len(partitions)==1 and 
         len(systs) <= chunk and 
         ((not ensembles) or ensembles[1]-ensembles[0] <= chunk) and
         ((not templates) or templates[1]-templates[0] == 1) ):
        assert len(filter(None,[systs,ensembles,templates]))<2
        mp = systematics.measurement_pars(partition=partitions[0])
        mp.update({'doVis':options.visualize,
                   'evalSystematics':systs,
                   'ensembles':ensembles})
        measurement(**mp)
    else:
        pass

