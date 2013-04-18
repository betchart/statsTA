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
from itertools import combinations

class fit(object):
    def __init__(self, label, signal, profileVars, R0_,
                 d_lumi, d_xs_dy, d_xs_st, tag, genPre, sigPre, dirIncrement, quiet = False,
                 hackZeroBins=False, alternateModel=False, defaults = {}, pllPoints=[]):

        self.label = label
        self.quiet = quiet
        self.profileVars = profileVars
        self.pllPoints = pllPoints
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
        self.model = model.topModel(channels, asymmetry=self.doAsymm, quiet=True,
                                    alternate=alternateModel)
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
        for j in range(10):
            minu.setStrategy(2)
            for i in range(10):
                self.status = minu.migrad()
                print i + 10*j,
                sys.stdout.flush()
                if not self.status: break
            if not self.status: break
            minu.setStrategy(1)
            minu.migrad()
        print
        if not self.pllPoints:
            print minu.minos(w.argSet(','.join(self.profileVars)))
        pll = nll.createProfile(w.argSet(','.join(self.profileVars)))
        pll.Print()
        pll.minuit().setStrategy(2)

        def pllEval(**kwargs) :
            for name,val in kwargs.items(): w.arg(name).setVal(val)
            return pll.getVal()

        errMin = {2:0.08, 1:0.005}[N]
        def vE(a) : return w.arg(a).getVal(), max(errMin, w.arg(a).getError())

        muSig = [vE(a) for a in self.profileVars]
        if not self.pllPoints:
            print muSig
            self.pllPoints = [tuple(m + i*s for i, (m, s) in zip(signs, muSig))
                              for signs in sorted(set(combinations([0, 1, 0, -1, 0], N)))]
            if N==1:
                for i in range(len(self.pllPoints)):
                    if self.pllPoints[i][0]<-1 : self.pllPoints[i] = ((muSig[0][0]-1)/2,)
            print self.pllPoints

        points = [p + (pllEval(**dict(zip(self.profileVars, p))),) 
                  for p in self.pllPoints]
        print zip(*points)[-1]

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
        self.fractionHats = [w.arg('f_%s_hat' % a).getVal() for a in ['gg','qg','qq','ag']]

        self.parb = parb
        self.muSig = muSig
        self.profPLL = pllEval(**dict(zip(self.profileVars, self.profVal)))
        self.fitPLL = pllEval(**dict(zip(self.profileVars,[m for m,_ in muSig])))
        choosePLL = ((N==3 and self.profErr[0,0] > 0 and self.profErr[1,1] > 0 
                      and self.profPLL < self.fitPLL) or 
                     (N==1 and self.profErr != None and self.profVal[0]>-1))
        self.best = self.profVal if choosePLL else [m for m,_ in muSig]
        if N==1 and not choosePLL: self.profErr = muSig[0][1]

        if not self.quiet:
            print zip(*muSig)[0], self.fitPLL
            print self.profVal, self.profPLL
            print self.profErr
            for item in ['d_qq','d_xs_tt','d_xs_wj','factor_elqcd','factor_muqcd','alphaL'][:-1 if N==1 else None]:
                print '\t',
                w.arg(item).Print()
        return

    @classmethod
    def fields(cls, doAsymm): 
        if not doAsymm:
            return ('#label d_qq error' +
                    'fhat_gg fhat_qg fhat_qq fhat_ag status')
        return ('#label fqq.Ac_y_qq  fqg.Ac_y_qg  XX  XY  YY fhat_gg ' +
                'fhat_qg fhat_qq fhat_ag Ac_y_gg_hat Ac_y_qq_hat Ac_y_qg_hat fitstatus')

    def __str__(self):
        if not self.doAsymm:
            return '\t'.join(str(i) for i in [self.label] + list(self.best) + [self.profErr] + self.fractionHats + [self.status])

        return '\t'.join(str(i) for i in [self.label] +
                         list(self.best*self.scales) +
                         [self.profErr[0,0]*self.scales[0]**2,
                          self.profErr[0,1]*self.scales[0]*self.scales[1],
                          self.profErr[1,1]*self.scales[1]**2] +
                         self.fractionHats+
                         list(self.scales) + [self.status]
                         )

class measurement(object):
    def __init__(self, label, signal, profile, R0_, hackZeroBins=False, alternateModel=False):
        write = open('data/' + '_'.join(label.split(',')) + '.txt', 'w')
        print >> write, fit.fields('QueuedBin' in signal)

        self.central = fit(signal=signal, profileVars=profile, R0_=R0_,
                           hackZeroBins=hackZeroBins, alternateModel=alternateModel,
                           **systematics.central())

        print >> write, str(self.central)

        if 'QueuedBin' in signal:
            para = self.central.parb.parametricEllipse(0.5)
            angle = self.central.parb.major_angle()
            points = [(x/w,y/w) for x,y,w in [para.dot([math.cos(angle+f*math.pi),
                                                        math.sin(angle+f*math.pi),
                                                        1])
                                              for f in [0, 0.5, 0.75, 1.0, 1.5]]]
            points.insert(0,self.central.parb.xymin)
            self.central.pllPoints = points

            vars_ = ['d_qq' if not alternateModel else 'alphaL_mag',
                     'falphaL', 'falphaT','R_ag',
                     'd_xs_tt', 'd_xs_wj', 'factor_elqcd', 'factor_muqcd']
        else:
            vars_ = ['d_qq','R_ag','d_xs_tt','d_xs_wj','factor_elqcd','factor_muqcd']

        defaults = dict([(v,self.central.model.w.arg(v).getVal()) for v in vars_])
        print '\n'.join(str(kv) for kv in defaults.items())

        syss = []
        for sys in systematics.systematics():
            pars = systematics.central()
            pars.update(sys)
            f = fit(signal=signal, profileVars=profile, R0_=R0_, quiet=False,
                    hackZeroBins=hackZeroBins, alternateModel=alternateModel,
                    defaults=defaults, pllPoints=list(self.central.pllPoints),
                    **pars)
            syss.append(f)
            print >> write, str(f)
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
