import ROOT as r
import roo
import systematics
from ensembles import ensemble_specs
from fitroutine import fit


class measurement(object):
    def __init__(self, label, signal, profile, R0_, hackZeroBins=False, 
                 doVis=False, evalSystematics=[], ensembles=None, ensSlice=(None,None)):
        self.outNameBase = 'data/' + '_'.join(label.split(',')) + ['_nosys',''][int(bool(evalSystematics))]
        write = open(self.outNameBase + '.txt', 'w')
        log = open(self.outNameBase + '.log', 'w')
        print >> write, fit.fields()


        self.SM = fit(signal=signal, profileVars=profile, R0_=R0_, log=log,
                      hackZeroBins=hackZeroBins, fixSM=True, **systematics.central())
        print >> log, self.SM.model.TTbarComponentsStr()

        self.central = fit(signal=signal, profileVars=profile, R0_=R0_, log=log,
                           hackZeroBins=hackZeroBins, **systematics.central())

        print >> write, str(self.central)
        self.central.model.print_n()
        tfile = r.TFile.Open(self.outNameBase+'.root','RECREATE')
        tree = self.central.ttree()
        tree.Write()
        tfile.Close()

        vars_ = ['slosh', 'falphaL', 'falphaT','R_ag',
                 'd_xs_tt', 'd_xs_wj', 'factor_elqcd', 'factor_muqcd']

        defaults = dict([(v,self.central.model.w.arg(v).getVal()) for v in vars_])
        print>>log, '\n'.join(str(kv) for kv in defaults.items())

        if doVis: self.central.model.visualize2D(printName=self.outNameBase+'.pdf')

        if ensembles: 
            pars = systematics.central()
            pars.update({'signal':signal, 'profileVars':profile, 'R0_':R0_, 'log':log, 'hackZeroBins':hackZeroBins})
            for ensPars in ensemble_specs():
                if ensPars['label'] not in ensembles: continue
                ensPars.update({'ensSlice':ensSlice})
                self.ensembles(pars, **ensPars)

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
    def ensembles(self, pars, letter='A', lumiFactor=1.0, ensSlice=(None,None), Nens=1000, label=''):
        wGen = self.central.model.w
        wGen.arg('lumi_factor').setVal(lumiFactor)
        pars['lumiFactor'] = lumiFactor
        if letter=='D':
            wGen.arg('falphaL').setVal(0)
            wGen.arg('falphaT').setVal(0)
            wGen.arg('slosh').setVal(self.SM.model.w.arg('slosh').getVal())
            wGen.arg('R_ag').setVal(self.SM.model.w.arg('R_ag').getVal())
        if letter=='C':
            for item in ['falphaL','falphaT','slosh','R_ag']:
                wGen.arg(item).setVal(self.SM.model.w.arg(item).getVal())
        if letter=='B':
            for item in ['falphaL','falphaT']:
                wGen.arg(item).setVal(-wGen.arg(item).getVal())


        truth = {'fitX': float(wGen.arg('falphaL')*self.central.scales[0]),
                 'fitY': float(wGen.arg('falphaT')*self.central.scales[1])}
        for item in fit.modelItems(): truth[item] = wGen.arg(item).getVal()

        mcstudy = r.RooMCStudy(wGen.pdf('model'),
                               wGen.argSet(','.join(self.central.model.observables+['channel'])),
                               r.RooFit.Binned(True),
                               r.RooFit.Extended(True)
                           )
        mcstudy.generate(Nens,0,True)
        with open('ensemble_%s_LF%d.txt'%(letter,100*lumiFactor),'w') as ensfile:
            with open('ensemble_%s_LF%d.log'%(letter,100*lumiFactor),'w') as enslog:
                pars['log'] = enslog
                print >> ensfile, '\t'.join(["#truth", '%(fitX)f\t%(fitY)f'%truth])
                print >> ensfile, fit.fields()
                for i in range(Nens)[slice(*ensSlice)]:
                    alt = mcstudy.genData(i)
                    pars['label'] = '%s_ens%d'%(label,i)
                    f = fit(altData=alt, **pars)
                    print >> ensfile, f
                    ensfile.flush()
                    enslog.flush()
                    tfile = r.TFile.Open(self.outNameBase+pars['label']+'.root','RECREATE')
                    tree = f.ttree(truth)
                    tree.Write()
                    tfile.Close()
                    
