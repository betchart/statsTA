import ROOT as r
from lib import roo
import systematics
from ensembles import ensemble_specs
from fitroutine import fit
import os

class measurement(object):
    def __init__(self, label, signal, profile, R0_, hackZeroBins=False,
                 doVis=False, evalSystematics=[], ensembles=None, ensSlice=(None,None), outDir='output/', templateID=None):
        os.system('mkdir -p %s' % outDir)
        self.outNameBase = (outDir + 
                            '_'.join(label.split(',')) + 
                            ('_t%03d'%templateID if templateID!=None else ''))

        with open(self.outNameBase + 'SM.log', 'w') as log:
            pars = systematics.central()
            pars['label'] += 'SM'
            self.SM = fit(signal=signal, profileVars=profile, R0_=R0_, log=log,
                          hackZeroBins=hackZeroBins, fixSM=True, **pars)
            print >> log, self.SM.model.TTbarComponentsStr()

        with open(self.outNameBase + '.log', 'w') as log:
            pars = systematics.central()
            if templateID!=None: pars['label'] = 'T%03d'%templateID
            self.central = fit(signal=signal, profileVars=profile, R0_=R0_, log=log,
                               hackZeroBins=hackZeroBins, templateID=templateID, **pars)
            self.central.model.print_n(log)
            self.central.ttreeWrite(self.outNameBase+'.root')

        defaults = dict([(v,self.central.model.w.arg(v).getVal()) 
                         for v in ['slosh', 'falphaL', 'falphaT','R_ag',
                                   'd_xs_tt', 'd_xs_wj', 'factor_elqcd', 'factor_muqcd']])

        if doVis: self.central.model.visualize2D(printName=self.outNameBase+'.pdf')

        if ensembles: 
            pars = systematics.central()
            pars.update({'signal':signal, 'profileVars':profile, 'R0_':R0_, 'log':log, 'hackZeroBins':hackZeroBins})
            for ensPars in ensemble_specs():
                if ensPars['label'] not in ensembles: continue
                ensPars.update({'ensSlice':ensSlice})
                self.ensembles(pars, **ensPars)

        for sys in systematics.systematics():
            if sys['label'] not in evalSystematics: continue
            pars = systematics.central()
            pars.update(sys)
            fname = self.outNameBase +'_sys_'+ sys['label'] + '.log'
            with open(fname, 'w') as log:
                f = fit(signal=signal, profileVars=profile, R0_=R0_,
                        quiet=False, hackZeroBins=hackZeroBins,
                        defaults=defaults, log=log, **pars)
                f.ttreeWrite(fname.replace('.log','.root'))

    @roo.quiet
    def ensembles(self, pars, letter='A', lumiFactor=1.0, ensSlice=(None,None), Nens=1000, label=''):
        wGen = self.central.model.w
        wGen.arg('lumi_factor').setVal(lumiFactor)
        pars['lumiFactor'] = lumiFactor

        if letter=='D':
            wGen.arg('falphaL').setVal(0)
            wGen.arg('falphaT').setVal(0)
            wGen.arg('slosh').setVal(self.SM.model.w.arg('slosh').getVal())
            wGen.arg('R_ag').setVal(self.SM.model.w.arg('R_ag').getVal())

        elif letter=='C':
            for item in ['falphaL','falphaT','slosh','R_ag']:
                wGen.arg(item).setVal(self.SM.model.w.arg(item).getVal())

        elif letter=='B':
            for item in ['falphaL','falphaT']:
                wGen.arg(item).setVal(-wGen.arg(item).getVal())

        else: assert letter=='A'

        truth = {'fitX': float(wGen.arg('falphaL').getVal()*self.central.scales[0]),
                 'fitY': float(wGen.arg('falphaT').getVal()*self.central.scales[1])}
        for item in fit.modelItems(): truth[item] = wGen.arg(item).getVal()

        mcstudy = r.RooMCStudy(wGen.pdf('model'),
                               wGen.argSet(','.join(self.central.model.observables+['channel'])),
                               r.RooFit.Binned(True),
                               r.RooFit.Extended(True)
                           )
        mcstudy.generate(Nens,0,True)
        for i in range(Nens)[slice(*ensSlice)]:
            alt = mcstudy.genData(i)
            pars['label'] = '%s_ens%03d'%(label,i)
            with open(self.outNameBase + pars['label'] + '.log', 'w') as log:
                pars['log']=log
                f = fit(altData=alt, **pars)
            f.ttreeWrite(self.outNameBase + pars['label'] + '.root', truth)
