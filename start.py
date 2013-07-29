#!/usr/bin/python

import sys
import ROOT as r
import systematics
from ensembles import ensemble_specs
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

    allEnsembles = [e['label'] for e in ensemble_specs()]
    if options.ensembles == True:
        print ','.join(allEnsembles)
        exit()
    ensembles = ([] if not options.ensembles else
                 allEnsembles if options.ensembles=='all' else
                 options.ensembles.split(','))

    ensSlice = ((None,None) if not options.ensSlice else tuple(int(i) for i in options.ensSlice.split(':')))
    templates = (None if not options.templates else tuple(options.ensembles))
    
    if options.batch:
        chunk = 10
        pass
    else:
        mp = systematics.measurement_pars(partition=partitions[0])
        mp.update({'doVis':options.visualize,
                   'evalSystematics':systs,
                   'ensembles':ensembles,
                   'ensSlice':ensSlice})
        print mp
        measurement(**mp)

