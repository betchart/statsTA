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
from itertools import combinations_with_replacement

class fit(object):
    def __init__(self, label, signal, profileVars, R0_,
                 d_lumi, d_xs_dy, d_xs_st, tag, genPre, sigPre, dirIncrement, quiet = False,
                 hackZeroBins=False, defaults = {}):

        self.label = label
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
        for k,v in defaults.items(): self.model.w.arg(k).setVal(v)
        for item in ['d_lumi', 'd_xs_dy', 'd_xs_st']: self.model.w.arg(item).setVal(eval(item))
        self.model.import_data()
        self.fitArgs = [r.RooFit.Extended(True), r.RooFit.NumCPU(1),
                        r.RooFit.PrintLevel(-1)]

        self.doFit()

    @roo.quiet
    def doFit(self):
        w = self.model.w
        N = len(self.profileVars)

        nll = w.pdf('model').createNLL(w.data('data'), *self.fitArgs[:-1])
        minu = r.RooMinuit(nll)
        minu.setPrintLevel(-1)
        minu.setNoWarn()
        for i in range(100):
            print i
            if not minu.migrad(): break
        print minu.improve()
        print minu.minos(w.argSet(','.join(self.profileVars)))
        pll = nll.createProfile(w.argSet(','.join(self.profileVars)))

        def pllEval(**kwargs) :
            for name,val in kwargs.items(): w.arg(name).setVal(val)
            return pll.getVal()

        def vE(a) : return w.arg(a).getVal(), w.arg(a).getError()

        muSig = [vE(a) for a in self.profileVars]
        points = [p + (pllEval(**dict(zip(self.profileVars, p))),) 
                  for p in [tuple(m + i*s for i, (m, s) in zip(signs, muSig))
                            for signs in combinations_with_replacement([0,-1, 1], N)]]
        print zip(*points)[2]

        oneSigmaNLL = {1: 0.5, 2: 1.14}[N]
        if N == 1:
            parb = parabola(points)
            self.profVal = parb.xmin,
            self.profErr = parb.dx(oneSigmaNLL)
        else:
            parb = paraboloid(points)
            self.profVal = parb.xymin
            self.profErr = parb.dxy(oneSigmaNLL)

        if self.doAsymm: 
            self.scales = np.array([w.arg(a).getVal() for a in ['Ac_y_ttqq', 'Ac_y_ttqg']])
            self.scalesPhi= [w.arg('Ac_phi_%s'%n) for n in ['ttqq','ttgg','ttag','ttqg','tt']]

        self.parb = parb
        self.muSig = muSig
        self.profPLL = pllEval(**dict(zip(self.profileVars,self.profVal)))
        self.fitPLL = pllEval(**dict(zip(self.profileVars,[m for m,_ in muSig])))

        if not self.quiet:
            print zip(*muSig)[0], self.fitPLL
            print self.profVal, self.profPLL
            print self.profErr
            for item in ['d_qq','d_xs_tt','d_xs_wj','factor_elqcd','factor_muqcd']:
                print '\t',
                w.arg(item).Print()
        return

        parb = paraboloid(points)
        oneSigmas = parb.dxy(oneSigmaNLL)
        print oneSigmas
        print math.sqrt(oneSigmas[0,0]), math.sqrt(oneSigmas[1,1])
        print parb.xymin
        param1sigma = parb.parametricEllipse(oneSigmaNLL)
        param2sigma = parb.parametricEllipse(twoSigmaNLL)
        with open('stat.txt', 'w') as wfile:
            for t in np.arange(0, 2 * math.pi + 0.00001, math.pi / 50):
                point1 = param1sigma.dot([math.cos(t), math.sin(t), 1])
                point2 = param2sigma.dot([math.cos(t), math.sin(t), 1])
                point1 /= point1[2]
                point2 /= point2[2]
                seq = list(parb.xymin) + list(point1[:2]) + list(point2[:2]) + scales
                print>>wfile, '\t'.join(str(f) for f in seq)
        print "Wrote stat.txt"

    @classmethod
    def fields(cls): return '#label Ac_y_qq  Ac_y_qg  XX  XY  YY'

    def __str__(self):
        return '\t'.join(str(i) for i in [self.label] +
                         list(self.profVal*self.scales) +
                         [self.profErr[0,0]*self.scales[0]**2,
                          self.profErr[0,1]*self.scales[0]*self.scales[1],
                          self.profErr[1,1]*self.scales[1]**2])

class measurement(object):
    def __init__(self, label, signal, profile, R0_, hackZeroBins=False):
        write = open('data/' + '_'.join(label.split(',')) + '.txt', 'w')
        print >> write, fit.fields()

        self.central = fit(signal=signal, profileVars=profile, R0_=R0_,
                           hackZeroBins=hackZeroBins, **systematics.central())

        print >> write, str(self.central)

        vars = ['alphaL', 'falphaL', 'falphaT','R_ag',
                'd_xs_tt', 'd_xs_wj', 'factor_elqcd', 'factor_muqcd']
        defaults = dict([(v,self.central.model.w.arg(v).getVal()) for v in vars])
        print '\n'.join(str(kv) for kv in defaults.items())

        syss = []
        for sys in systematics.systematics():
            pars = systematics.central()
            pars.update(sys)
            try:
                f = fit(signal=signal, profileVars=profile, R0_=R0_, quiet=True,
                        defaults=defaults, **pars)
                syss.append(f)
                print zip(*f.muSig)[0], f.fitPLL
                print f.profVal, f.profPLL
                print >> write, str(f)
            except:
                print >> write, '#', pars['label'], "FAIL"
            write.flush()

        write.close()

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
