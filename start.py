#!/usr/bin/python

import sys
import math
import ROOT as r
import numpy as np

import inputs
import model
import roo
import systematics
import utils
from paraboloid import paraboloid
from parabola import parabola
from itertools import combinations
from asymmNames import genNameX,genNameY

class fit(object):
    def __init__(self, label, signal, profileVars, R0_,
                 d_lumi, d_xs_dy, d_xs_st, tag, genPre, sigPre, dirIncrement, genDirPre, 
                 quiet = False, hackZeroBins=False, defaults = {}, pllPoints=[],
                 log=None, fixSM=False,altData=None, lumiFactor=1.0):

        self.label = label
        self.quiet = quiet
        self.profileVars = profileVars
        self.pllPoints = pllPoints
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

        if fixSM:
            fixVars = (['R_ag','slosh','falphaL','falphaT'] if self.doAsymm else
                       ['d_qq','R_ag'])
            for item in fixVars: self.model.w.arg(item).setConstant()
            nll = self.model.w.pdf('model').createNLL(self.model.w.data('data'), *self.fitArgs[:-1])
            minu = r.RooMinuit(nll)
            minu.setPrintLevel(-1)
            minu.setNoWarn()
            minu.setStrategy(2)
            minu.migrad()
        else:
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
                print>>self.log, i + 10*j,
                self.log.flush()
                if not self.status: break
            if not self.status: break
            minu.setStrategy(1)
            minu.migrad()
        print>>self.log
        if not self.pllPoints:
            print>>self.log, minu.minos(w.argSet(','.join(self.profileVars)))
        self.NLL = nll.getVal()
        pll = nll.createProfile(w.argSet(','.join(self.profileVars)))
        print>>self.log, roo.str(nll)
        print>>self.log, roo.str(pll)
        pll.minuit().setStrategy(2)

        def pllEval(**kwargs) :
            for name,val in kwargs.items(): w.arg(name).setVal(val)
            return pll.getVal()

        errMin = {2:0.08, 1:0.005}[N]
        def vE(a) : return w.arg(a).getVal(), max(errMin, w.arg(a).getError())

        muSig = [vE(a) for a in self.profileVars]
        if not self.pllPoints:
            print>>self.log, muSig
            self.pllPoints = [tuple(m + i*s for i, (m, s) in zip(signs, muSig))
                              for signs in sorted(set(combinations([0, 1, 0, -1, 0], N)))]
            if N==1:
                for i in range(len(self.pllPoints)):
                    if self.pllPoints[i][0]<-1 : self.pllPoints[i] = ((muSig[0][0]-1)/2,)
            print>>self.log, self.pllPoints

        points = [p + (pllEval(**dict(zip(self.profileVars, p))),) 
                  for p in self.pllPoints]
        print>>self.log, zip(*points)[-1]

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
            self.correction = w.arg('Ac_y_ttgg').getVal() * w.arg('f_gg').getVal()
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
            print>>self.log, zip(*muSig)[0], self.fitPLL
            print>>self.log, self.profVal, self.profPLL
            print>>self.log, self.profErr
            for item in ['d_qq','d_xs_tt','d_xs_wj',
                         'factor_elqcd','factor_muqcd',
                         'f_gg','f_qq','f_qg','f_ag',
                         'R_ag','slosh','alphaL'][:-3 if N==1 else None]:
                print>>self.log, '\t', roo.str(w.arg(item))
        self.pll = pll
        return

    @classmethod
    def fields(cls, doAsymm): 
        if not doAsymm:
            return ('#label d_qq error' +
                    'fhat_gg fhat_qg fhat_qq fhat_ag status')
        return ('#label fqq.Ac_y_qq  fqg.Ac_y_qg  XX  XY  YY fhat_gg ' +
                'fhat_qg fhat_qq fhat_ag Ac_y_qq_hat Ac_y_qg_hat f_gg.Ac_y_gg fitstatus')

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
                          self.status]
                         )

class measurement(object):
    def __init__(self, label, signal, profile, R0_, hackZeroBins=False, doVis=False, doSys=False, doEnsembles=True):
        self.isAsymmetry = 'QueuedBin' in signal
        outNameBase = 'data/' + '_'.join(label.split(',')) + ['_nosys',''][int(doSys)]
        write = open(outNameBase + '.txt', 'w')
        log = open(outNameBase + '.log', 'w')
        print >> write, fit.fields('QueuedBin' in signal)


        self.SM = fit(signal=signal, profileVars=profile, R0_=R0_, log=log,
                      hackZeroBins=hackZeroBins, fixSM=True, **systematics.central())
        print >> log, self.SM.model.TTbarComponentsStr()
        #visCanvas = self.SM.model.visualize2D()
        #utils.tCanvasPrintPdf(visCanvas, outNameBase, verbose=False, title='SM', option='(')


        self.central = fit(signal=signal, profileVars=profile, R0_=R0_, log=log,
                           hackZeroBins=hackZeroBins, **systematics.central())

        print >> write, str(self.central)
        self.central.model.print_n()

        if self.isAsymmetry:
            para = self.central.parb.parametricEllipse(0.5)
            angle = self.central.parb.major_angle()
            points = [(x/w,y/w) for x,y,w in [para.dot([math.cos(angle+f*math.pi),
                                                        math.sin(angle+f*math.pi),
                                                        1])
                                              for f in [0, 0.5, 0.75, 1.0, 1.5]]]
            points.insert(0,self.central.parb.xymin)
            self.central.pllPoints = points

            vars_ = ['slosh', 'falphaL', 'falphaT','R_ag',
                     'd_xs_tt', 'd_xs_wj', 'factor_elqcd', 'factor_muqcd']
        else:
            vars_ = ['d_qq','R_ag','d_xs_tt','d_xs_wj','factor_elqcd','factor_muqcd']

        defaults = dict([(v,self.central.model.w.arg(v).getVal()) for v in vars_])
        print>>log, '\n'.join(str(kv) for kv in defaults.items())

        if doVis: self.central.model.visualize2D(printName=outNameBase+'.pdf')
        #utils.tCanvasPrintPdf(visCanvas, outNameBase, verbose=False, title='central')

        if doEnsembles: 
            ensPars = systematics.central()
            ensPars.update({'signal':signal, 'profileVars':profile, 'R0_':R0_, 'log':log, 'hackZeroBins':hackZeroBins})
            self.ensembles(ensPars, lumiFactor=0.5)

        syss = []
        for sys in [[],systematics.systematics()][int(doSys)]:
            pars = systematics.central()
            pars.update(sys)
            f = fit(signal=signal, profileVars=profile, R0_=R0_, quiet=False,
                    hackZeroBins=hackZeroBins,
                    defaults=defaults, pllPoints=list(self.central.pllPoints), log=log,
                    **pars)
            #f.model.visualize(visCanvas)
            #utils.tCanvasPrintPdf(visCanvas, outNameBase, verbose=False, title=pars['label'])
            syss.append(f)
            print >> write, str(f)
            write.flush()

        if False:
            items = 'f_qq,f_gg,f_qg,f_ag'.split(',')
            pars = dict([(p,(self.central.model.w.arg(p))) for p in 'R_ag,slosh'.split(',')])
            values = dict([(i,[]) for i in items])
            for name,p in pars.items():
                p.Print()
                mean = p.getVal()
                err = p.getError()
                probes = [mean+err,mean-err,mean]
                for pro in probes:
                    p.setVal(pro)
                    for i,v in values.items():
                        v.append(self.central.model.w.arg(i).getVal())
            for i,v in values.items():
                print i, '%.4f(%+.4f,%+.4f)'%(v[-1],min(v)-v[-1], max(v)-v[-1])


        #visCanvas.Clear()
        #utils.tCanvasPrintPdf(visCanvas, outNameBase, verbose=True, option=')')
        write.close()
        log.close()

    @roo.quiet
    def ensembles(self, pars, ens='A', lumiFactor=1.0):
        self.central.model.w.arg('lumi_factor').setVal(lumiFactor)
        pars['lumiFactor'] = lumiFactor
        if ens=='B':
            self.central.model.w.arg('falphaL').setVal(0)
            self.central.model.w.arg('falphaT').setVal(0)
            self.central.model.w.arg('slosh').setVal(self.SM.model.w.arg('slosh').getVal())
            self.central.model.w.arg('R_ag').setVal(self.SM.model.w.arg('R_ag').getVal())

        mcstudy = r.RooMCStudy(self.central.model.w.pdf('model'),
                               self.central.model.w.argSet(','.join(self.central.model.observables+['channel'])),
                               r.RooFit.Binned(True),
                               r.RooFit.Extended(True)
                           )
        skip=0
        Nens = 10
        mcstudy.generate(Nens,0,True)
        with open('ensemble_%s_LF%d.txt'%(ens,100*lumiFactor),'w') as ensfile:
            fAqq = self.central.model.w.arg('falphaL').getVal() * self.central.model.w.arg('Ac_y_ttqq').getVal()
            fAqg = self.central.model.w.arg('falphaT').getVal() * self.central.model.w.arg('Ac_y_ttqg').getVal()
            print >> ensfile, '\t'.join(["#truth", '%f'%fAqq, '%f'%fAqg])
            print >> ensfile, fit.fields(self.isAsymmetry), 'NLL'
            for i in range(skip,Nens):
                alt = mcstudy.genData(i)
                pars['label'] = 'ens%d'%i
                f = fit(altData=alt, **pars)
                print >> ensfile, f, f.NLL
                ensfile.flush()


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
