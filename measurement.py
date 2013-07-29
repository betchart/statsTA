import ROOT as r
import roo
import systematics
from fitroutine import fit


class measurement(object):
    def __init__(self, label, signal, profile, R0_, hackZeroBins=False, 
                 doVis=False, evalSystematics=[], ensembles=None):
        outNameBase = 'data/' + '_'.join(label.split(',')) + ['_nosys',''][int(bool(evalSystematics))]
        write = open(outNameBase + '.txt', 'w')
        log = open(outNameBase + '.log', 'w')
        print >> write, fit.fields()


        self.SM = fit(signal=signal, profileVars=profile, R0_=R0_, log=log,
                      hackZeroBins=hackZeroBins, fixSM=True, **systematics.central())
        print >> log, self.SM.model.TTbarComponentsStr()


        self.central = fit(signal=signal, profileVars=profile, R0_=R0_, log=log,
                           hackZeroBins=hackZeroBins, **systematics.central())

        print >> write, str(self.central)
        self.central.model.print_n()

        vars_ = ['slosh', 'falphaL', 'falphaT','R_ag',
                 'd_xs_tt', 'd_xs_wj', 'factor_elqcd', 'factor_muqcd']

        defaults = dict([(v,self.central.model.w.arg(v).getVal()) for v in vars_])
        print>>log, '\n'.join(str(kv) for kv in defaults.items())

        if doVis: self.central.model.visualize2D(printName=outNameBase+'.pdf')

        if ensembles: 
            ensPars = systematics.central()
            ensPars.update({'signal':signal, 'profileVars':profile, 'R0_':R0_, 'log':log, 'hackZeroBins':hackZeroBins})
            self.ensembles(ensPars, lumiFactor=1.0, ens='C', ensSlice=ensembles)

        syss = []
        for sys in systematics.systematics():
            if sys['label'] not in evalSystematics: continue
            pars = systematics.central()
            pars.update(sys)
            f = fit(signal=signal, profileVars=profile, R0_=R0_,
                    quiet=False, hackZeroBins=hackZeroBins,
                    defaults=defaults, log=log, **pars)
            syss.append(f)
            print >> write, str(f)
            write.flush()

        write.close()
        log.close()

    @roo.quiet
    def ensembles(self, pars, ens='A', lumiFactor=1.0, ensSlice=(None,None)):
        self.central.model.w.arg('lumi_factor').setVal(lumiFactor)
        pars['lumiFactor'] = lumiFactor
        if ens=='D':
            self.central.model.w.arg('falphaL').setVal(0)
            self.central.model.w.arg('falphaT').setVal(0)
            self.central.model.w.arg('slosh').setVal(self.SM.model.w.arg('slosh').getVal())
            self.central.model.w.arg('R_ag').setVal(self.SM.model.w.arg('R_ag').getVal())
        if ens=='C':
            for item in ['falphaL','falphaT','slosh','R_ag']:
                self.central.model.w.arg(item).setVal(self.SM.model.w.arg(item).getVal())
        if ens=='B':
            for item in ['falphaL','falphaT']:
                self.central.model.w.arg(item).setVal(-self.central.model.w.arg(item).getVal())

        mcstudy = r.RooMCStudy(self.central.model.w.pdf('model'),
                               self.central.model.w.argSet(','.join(self.central.model.observables+['channel'])),
                               r.RooFit.Binned(True),
                               r.RooFit.Extended(True)
                           )
        Nens = 1000
        mcstudy.generate(Nens,0,True)
        with open('ensemble_%s_LF%d.txt'%(ens,100*lumiFactor),'w') as ensfile:
            with open('ensemble_%s_LF%d.log'%(ens,100*lumiFactor),'w') as enslog:
                pars['log'] = enslog
                fAqq = self.central.model.w.arg('falphaL').getVal() * self.central.model.w.arg('Ac_y_ttqq').getVal()
                fAqg = self.central.model.w.arg('falphaT').getVal() * self.central.model.w.arg('Ac_y_ttqg').getVal()
                print >> ensfile, '\t'.join(["#truth", '%f'%fAqq, '%f'%fAqg])
                print >> ensfile, fit.fields()
                for i in range(Nens)[slice(*ensSlice)]:
                    alt = mcstudy.genData(i)
                    pars['label'] = 'ens%d'%i
                    f = fit(altData=alt, **pars)
                    print >> ensfile, f
                    ensfile.flush()
                    enslog.flush()
