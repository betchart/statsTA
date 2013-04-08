import sys
import math
import ROOT as r
import numpy as np

import model
import roo
from falphaLSlice import falphaLSlice
from paraboloid import paraboloid


class topAsymmFit(object):
    @roo.quiet
    def __init__(self, dist, provar, tag):
        self.model = model.topModel(dist=dist, asymmetry=('QueuedBin' in dist), quiet=True)
        self.model.import_data()
        w = self.model.w

        print '\n'.join(str(i) for i in ['', self.model.channels['el'],
                                         '', self.model.channels['mu'], ''])

        self.fitArgs = [r.RooFit.Extended(True), r.RooFit.NumCPU(4),
                        r.RooFit.PrintLevel(-1), r.RooFit.Save()]
        for i in reversed(range(3)):
            central = w.pdf('model').fitTo(w.data('data'), *self.fitArgs[:-1 if i else None])
        central.Print()

        self.print_fracs(w)
        self.print_n(w)
        print

        faL, faLe = w.arg('falphaL').getVal(), w.arg('falphaL').getError()
        faT, faTe = w.arg('falphaT').getVal(), w.arg('falphaT').getError()

        nll = w.pdf('model').createNLL(w.data('data'), *self.fitArgs[:-2])
        pll = nll.createProfile(w.argSet('falphaL,falphaT'))

        def point(faL_, faT_):
            w.arg('falphaL').setVal(faL_)
            w.arg('falphaT').setVal(faT_)
            return faL_, faT_, pll.getVal()

        points = [point(faL + dfaL, faT + dfaT) for dfaL, dfaT in
                  [(0, 0), (-faLe, 0), (faLe, 0), (0, -faTe), (0, faTe), (faLe/2, faTe/2)]]
        parb = paraboloid(points)
        oneSigmaNLL = 1.14
        twoSigmaNLL = 3.0
        oneSigmas = parb.dxy(oneSigmaNLL)
        print oneSigmas
        print math.sqrt(oneSigmas[0,0]), math.sqrt(oneSigmas[1,1])
        param1sigma = parb.parametricEllipse(oneSigmaNLL)
        param2sigma = parb.parametricEllipse(twoSigmaNLL)
        scales = [w.arg(a).getVal() for a in ['Ac_y_ttqq', 'Ac_y_ttqg']]
        with open('points.txt', 'w') as wfile:
            for t in np.arange(0, 2 * math.pi + 0.00001, math.pi / 50):
                point1 = param1sigma.dot([math.cos(t), math.sin(t), 1])
                point2 = param2sigma.dot([math.cos(t), math.sin(t), 1])
                point1 /= point1[2]
                point2 /= point2[2]
                seq = list(parb.xymin) + list(point1[:2]) + list(point2[:2]) + scales
                print>>wfile, '\t'.join(str(f) for f in seq)
        print "Wrote points.txt"

    def print_fracs(self, w):
        for item in \
                ['lumi_mu', 'lumi_el', 'f_gg', 'f_qg', 'f_qq', 'f_ag'] + \
                ['xs_' + i for i in self.model.channels['el'].samples if i != 'data']:
            print "%s: %.04f" % (item, w.arg(item).getVal())
        print

    def print_n(self, w):
        length = 24
        tots = {'el': 0, 'mu': 0}
        print (' ').join(i.rjust(8) for i in [''] + tots.keys())
        for xs in ['tt', 'wj', 'mj', 'st', 'dy']:
            if xs == 'data': continue
            print xs.rjust(length / 3),
            for chan in tots:
                val = w.arg('expect_%s_%s' % (chan, xs)).getVal()
                tots[chan] += val
                print ("%d" % val).rjust(length / 3),
            print
        print '-' * (3 + length)
        print ' '.join(["tot".ljust(length / 3)] +
                       [("%d" % t).rjust(length / 3) for chan, t in tots.items()])
        print


if __name__ == '__main__':
    r.gROOT.SetBatch(1)

    distributions = ['fitTopQueuedBin7TridiscriminantWTopQCD',
                     'fitTopPtOverSumPt_triD',
                     'fitTopTanhRapiditySum_triD',
                     'TridiscriminantQQggQg_triD'
                     ]
    if len(sys.argv) > 1:
        iDist = sys.argv[1]
    else:
        print '\n'.join('%d %s' % iD for iD in enumerate(distributions))
        iDist = raw_input("Which?")
    dist = distributions[int(iDist)] if int(iDist) in range(len(distributions)) else ''
    if not dist:
        print 'not an option'
        exit(0)
    else: print "fitting", dist

    variables = ['alphaL', 'R_ag']
    provar = variables[0]

    topAsymmFit(dist, provar, tag=None)
