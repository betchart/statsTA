import math
import numpy as np
import ROOT as r
import roo
import inputs
import model
from enclosing_ellipse import enclosing_ellipse
from asymmNames import genNameX,genNameY


class fit(object):
    def __init__(self, label, signal, profileVars, R0_,
                 d_lumi, d_xs_dy, d_xs_st, tag, genPre, sigPre, dirIncrement, genDirPre, 
                 quiet = False, hackZeroBins=False, defaults = {},
                 log=None, fixSM=False,altData=None, lumiFactor=1.0):

        self.label = label
        self.quiet = quiet
        self.profileVars = profileVars
        self.doAsymm = 'QueuedBin' in signal
        self.log = log if log else sys.stdout
        if type(R0_) == tuple:
            diffR0_ = R0_[1]
            R0_ = R0_[0]
        else: diffR0_ = None
        prePre = dirIncrement in [0,4,5]
        channels = dict([((lep,part),
                          inputs.channel_data(lep, part, tag, signal, sigPre, 
                                              "R%02d" % (R0_ + dirIncrement),
                                              genDirPre, prePre = prePre,
                                              hackZeroBins=hackZeroBins and 'QCD'==part))
                         for lep in ['el', 'mu']
                         for part in ['top', 'QCD']
                         ])
        channels['gen'] = inputs.channel_data('mu', 'top', tag,
                                              '%s; %s'%(genNameX,genNameY),
                                              sigPrefix = sigPre if dirIncrement in [0,4,5] else '',
                                              dirPrefix=genDirPre, genDirPre=genDirPre, 
                                              getTT=True, prePre = prePre)

        #for n,c in channels.items() :
        #    if n == 'gen': continue
        #    print n
        #    print c.asymmStr()
        #    print

        if diffR0_ :
            for lepPart,chan in channels.items():
                if type(lepPart) != tuple: continue
                lep,part = lepPart
                chan.subtract(inputs.channel_data(lep,part,tag,signal,sigPre,
                                                  "R%02d" % (diffR0_ + dirIncrement),
                                                  genDirPre, prePre = prePre ))

        print "###", label
        print>>self.log, "###", label
        self.model = model.topModel(channels, asymmetry=self.doAsymm, quiet=True)
        self.model.w.arg('lumi_factor').setVal(lumiFactor)
        for k,v in defaults.items(): self.model.w.arg(k).setVal(v)
        for item in ['d_lumi', 'd_xs_dy', 'd_xs_st']: self.model.w.arg(item).setVal(eval(item))

        self.fitArgs = [r.RooFit.Extended(True), r.RooFit.NumCPU(1),
                        r.RooFit.PrintLevel(-1)]
        self.model.import_data(altData)

        if fixSM: self.doSM()
        else:
            self.doFit()

    @roo.quiet
    def doSM(self):
        fixVars = (['R_ag','slosh','falphaL','falphaT'] if self.doAsymm else
                   ['d_qq','R_ag'])
        for item in fixVars: self.model.w.arg(item).setConstant()
        nll = self.model.w.pdf('model').createNLL(self.model.w.data('data'), *self.fitArgs[:-1])
        minu = r.RooMinuit(nll)
        minu.setPrintLevel(-1)
        minu.setNoWarn()
        minu.setStrategy(2)
        minu.migrad()

    @roo.quiet
    def doFit(self):
        w = self.model.w
        nll = w.pdf('model').createNLL(w.data('data'), *self.fitArgs[:-1])
        minu = r.RooMinuit(nll)
        minu.setPrintLevel(-1)
        minu.setNoWarn()
        for j in range(10):
            minu.setStrategy(2)
            for i in range(10):
                self.status = minu.migrad()
                print>>self.log, i + 10*j,
                self.log.flush()
                if not self.status: break
            if not self.status: break
            minu.setStrategy(1)
            minu.migrad()
        print>>self.log
        errMin = 0.12
        oneSigmaNLL = 1.14
        p0 = [w.arg(a).getVal() for a in self.profileVars]
        self.NLL = nll.getVal()
        pll = nll.createProfile(w.argSet(','.join(self.profileVars)))
        print>>self.log, roo.str(nll)
        print>>self.log, roo.str(pll)
        pll.minimizer().setStrategy(2)

        pllCache = {}
        def pllEval(p):
            p = tuple(p)[:2]
            if p not in pllCache:
                for name,val in zip(self.profileVars,p): w.arg(name).setVal(val)
                pllCache[p] = pll.getVal()
            return pllCache[p]

        pllPoints = None
        while not pllPoints:
            cands = [(p0[0]+errMin*math.cos(t),p0[1]+errMin*math.sin(t))
                     for t in np.arange(0,2*math.pi,math.pi/8)]
            minP = min(cands, key=pllEval)
            if pllEval(minP) < pllEval(p0):
                print 'new minimum'
                p0 = minP
                pll.clearAbsMin()
                pll.bestFitObs()
                pllEval(p0)
                print p0
            else: pllPoints = cands
        p0 += pllEval(p0),
        
        targetPLL = oneSigmaNLL + p0[2]
        def contourIntersect(guess,epsilon=0.01,xepsilon=0.001):
            guess += pllEval(guess),
            def point(g): return (p0[0] + g * (guess[0]-p0[0]),
                                  p0[1] + g * (guess[1]-p0[1]))


            def bsearch(lo,hi):
                assert pllEval(lo) < targetPLL
                assert pllEval(hi) > targetPLL
                p = tuple(0.5*(lo[i]+hi[i]) for i in [0,1])
                PLL = pllEval(p)
                if math.sqrt(sum((lo[i]-hi[i])**2 for i in [0,1])) < xepsilon: return p
                if PLL < targetPLL-epsilon: return bsearch(p,hi)
                if PLL > targetPLL+epsilon: return bsearch(lo,p)
                return p

            mlo = guess if pllEval(guess) < targetPLL else None
            mhi = guess if targetPLL < pllEval(guess) else None
            iteration = 0
            while True:
                p = point( math.sqrt((targetPLL-p0[2]) / (guess[2]-p0[2])) )
                if abs(targetPLL-pllEval(p)) < epsilon: return p
                if pllEval(p) < targetPLL and (not mlo or pllEval(mlo) < pllEval(p)): mlo = p
                if pllEval(p) > targetPLL and (not mhi or pllEval(p) < pllEval(mhi)): mhi = p
                if iteration>5 and mlo and mhi: return bsearch(mlo,mhi)
                guess = p + (pllEval(p),)
                iteration += 1

        points = [contourIntersect(p) for p in pllPoints]
        print '\n'.join(str(a) for a in [p0] + points)
        print
        ell = enclosing_ellipse([p[:2] for p in points],p0[:2])
        print '\n'.join([str((x/W,y/W)) for x,y,W in [np.dot( ell.parametric, [math.cos(t),math.sin(t),1]) for t in np.arange(0,2*math.pi,math.pi/8)]])
        
        self.profVal = p0[:2]
        self.profErr = ell.sigmas2
        self.profPLL = p0[2]
        self.best = self.profVal

        self.scales = np.array([w.arg(a).getVal() for a in ['Ac_y_ttqq', 'Ac_y_ttqg']])
        self.scalesPhi= [w.arg('Ac_phi_%s'%n).getVal() for n in ['ttqq','ttgg','ttag','ttqg','tt']]
        self.correction = w.arg('Ac_y_ttgg').getVal() * w.arg('f_gg').getVal()
        self.fractionHats = [w.arg('f_%s_hat' % a).getVal() for a in ['gg','qg','qq','ag']]

        if not self.quiet:
            print>>self.log, self.profVal, self.profPLL
            print>>self.log, self.profErr
            for item in ['d_qq','d_xs_tt','d_xs_wj',
                         'factor_elqcd','factor_muqcd',
                         'f_gg','f_qq','f_qg','f_ag',
                         'R_ag','slosh','alphaL']:
                print>>self.log, '\t', roo.str(w.arg(item))
        self.pll = pll
        return

    @classmethod
    def fields(cls, doAsymm): 
        if not doAsymm:
            return ('#label d_qq error' +
                    'fhat_gg fhat_qg fhat_qq fhat_ag status')
        return ('#label fqq.Ac_y_qq  fqg.Ac_y_qg  XX  XY  YY fhat_gg ' +
                'fhat_qg fhat_qq fhat_ag Ac_y_qq_hat Ac_y_qg_hat f_gg.Ac_y_gg fitstatus NLL Ac_phi_qq_hat Ac_phi_gg_hat Ac_phi_ag_hat Ac_phi_qg_hat')

    def __str__(self):
        if not self.doAsymm:
            return '\t'.join(str(i) for i in [self.label] + list(self.best) +
                             [self.profErr] + self.fractionHats + [self.status])

        return '\t'.join(str(i) for i in [self.label] +
                         list(self.best*self.scales) +
                         [self.profErr[0,0]*self.scales[0]**2,
                          self.profErr[0,1]*self.scales[0]*self.scales[1],
                          self.profErr[1,1]*self.scales[1]**2] +
                         self.fractionHats+
                         list(self.scales) +
                         [self.correction,
                          self.status,
                          self.NLL]+
                         self.scalesPhi[:4]
                         )