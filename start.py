#!/usr/bin/python

import sys
import math
import ROOT as r
import numpy as np

import inputs
import model
import roo
import systematics
from falphaLSlice import falphaLSlice
from paraboloid import paraboloid
from parabola import parabola

class fit(object):
    def __init__(self, label, signal, profileVars, R0_,
                 d_lumi, d_xs_dy, d_xs_st, tag, genPre, sigPre, dirIncrement, quiet = False,
                 hackZeroBins=False):

        self.quiet = quiet
        self.profileVars = profileVars
        self.doAsymm = 'QueuedBin' in signal
        if type(R0_) == tuple:
            diffR0_ = R0_[1]
            R0_ = R0_[0]
        else: diffR0_ = None
        channels = dict([((lep,part),
                          inputs.channel_data(lep, part, tag, signal, sigPre, 
                                              "R%02d" % (R0_ + dirIncrement),
                                              prePre = not dirIncrement, 
                                              hackZeroBins=hackZeroBins and 'QCD'==part))
                         for lep in ['el', 'mu']
                         for part in ['top', 'QCD']
                         ])
        channels['gen'] = inputs.channel_data('mu', 'top', tag,
                                              'genTopDeltaBetazRel; genTopPhiBoost',
                                              sigPrefix = sigPre if not dirIncrement else '',
                                              dirPrefix="R01", getTT=True,
                                              prePre = not dirIncrement)
        
        if diffR0_ :
            for lepPart,chan in channels.items():
                if type(lepPart) != tuple: continue
                lep,part = lepPart
                chan.subtract(inputs.channel_data(lep,part,tag,signal,sigPre,
                                                  "R%02d" % (diffR0_ + dirIncrement)))

        print "###", label
        self.model = model.topModel(channels, asymmetry=self.doAsymm, quiet=True)
        for item in ['d_lumi', 'd_xs_dy', 'd_xs_st']:
            self.model.w.arg(item).setVal(eval(item))
        self.model.import_data()
        self.fitArgs = [r.RooFit.Extended(True), r.RooFit.NumCPU(1),
                        r.RooFit.PrintLevel(-1), r.RooFit.Save()]

        self.doCentral()

    @roo.quiet
    def doCentral(self):
        w = self.model.w
        self.central = w.pdf('model').fitTo(w.data('data'), *self.fitArgs)
        if not self.quiet:
            self.central.Print()
            w.arg('d_qq').Print()
            if self.doAsymm: w.arg('alphaT').Print()

    @roo.quiet
    def doProfile(self):
        w = self.model.w
        meanSigs = [(w.arg(a).getVal(),w.arg(a).getError()) for a in self.profileVars]

        nll = w.pdf('model').createNLL(w.data('data'), *self.fitArgs[:-2])
        pll = nll.createProfile(w.argSet(','.join(self.profileVars)))

        def point(*vals):
            for name,val in zip(self.profileVars,vals) :
                w.arg(name).setVal(val)
            return tuple(vals) + (pll.getVal(),)

        if not self.doAsymm:
            mean,sigma = meanSigs[0]
            points = [point(*p) for p in [(mean,),(mean+sigma,),(mean-sigma,)]]
            parb = parabola(points)
            oneSigmaNLL = 0.5
            oneSigmas = parb.dx(oneSigmaNLL)
            print oneSigmas
            print parb.xmin
            return

        (faL,faLe),(faT,faTe) = meanSigs
        points = [point(faL + dfaL, faT + dfaT) for dfaL, dfaT in
                  [(0, 0), (-faLe, 0), (faLe, 0), (0, -faTe), (0, faTe), (-faLe, faTe)]]

        print '\n'.join(str(p) for p in points)
        parb = paraboloid(points)
        oneSigmaNLL = 1.14
        twoSigmaNLL = 3.0
        oneSigmas = parb.dxy(oneSigmaNLL)
        print oneSigmas
        print math.sqrt(oneSigmas[0,0]), math.sqrt(oneSigmas[1,1])
        print parb.xymin
        param1sigma = parb.parametricEllipse(oneSigmaNLL)
        param2sigma = parb.parametricEllipse(twoSigmaNLL)
        scales = [w.arg(a).getVal() for a in ['Ac_y_ttqq', 'Ac_y_ttqg']]
        with open('stat.txt', 'w') as wfile:
            for t in np.arange(0, 2 * math.pi + 0.00001, math.pi / 50):
                point1 = param1sigma.dot([math.cos(t), math.sin(t), 1])
                point2 = param2sigma.dot([math.cos(t), math.sin(t), 1])
                point1 /= point1[2]
                point2 /= point2[2]
                seq = list(parb.xymin) + list(point1[:2]) + list(point2[:2]) + scales
                print>>wfile, '\t'.join(str(f) for f in seq)
        print "Wrote stat.txt"


class measurement(object):
    def __init__(self, label, signal, profile, R0_, hackZeroBins=False):
        self.central = fit(signal=signal, profileVars=profile, R0_=R0_,
                           hackZeroBins=hackZeroBins, **systematics.central())
        self.central.doProfile()
        return
        syss = []
        for sys in systematics.systematics():
            pars = systematics.central()
            pars.update(sys)
            f = fit(signal=signal, profile=profile, R0_=R0_, quiet=True, **pars)
            try: pass
            except: print 'failed', pars['label']
            else: syss.append(f)


def query(items, default = ()):
    display = '\n'.join('%d %s' % iD for iD in enumerate(items)) + '\nWhich? '
    i = eval(next(default, "input(display)"))
    return items[i] if 0 <= i <= len(items) else None


if __name__ == '__main__':
    r.gROOT.SetBatch(1)

    defaults = iter(sys.argv[1:])
    mp = [query(systematics.measurements, defaults),
          query(systematics.partitions, defaults)]

    if any(n == None for n in mp):
        print 'invalid'
        exit(0)

    print ','.join(mp)
    pars = systematics.measurement_pars(*mp)
    for p,v in pars.items():
        print "   %s: %s" % (p, str(v))
    print

    measurement(','.join(mp), **pars)
