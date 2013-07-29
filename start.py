#!/usr/bin/python

import sys
import ROOT as r
import roo
import systematics
from fitroutine import fit
#import utils


class measurement(object):
    def __init__(self, label, signal, profile, R0_, hackZeroBins=False, doVis=False, doSys=False, doEnsembles=False):
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

        vars_ = ['slosh', 'falphaL', 'falphaT','R_ag',
                 'd_xs_tt', 'd_xs_wj', 'factor_elqcd', 'factor_muqcd']

        defaults = dict([(v,self.central.model.w.arg(v).getVal()) for v in vars_])
        print>>log, '\n'.join(str(kv) for kv in defaults.items())

        if doVis: self.central.model.visualize2D(printName=outNameBase+'.pdf')
        #utils.tCanvasPrintPdf(visCanvas, outNameBase, verbose=False, title='central')

        if doEnsembles: 
            ensPars = systematics.central()
            ensPars.update({'signal':signal, 'profileVars':profile, 'R0_':R0_, 'log':log, 'hackZeroBins':hackZeroBins})
            self.ensembles(ensPars, lumiFactor=1.0, ens='C')

        syss = []
        for sys in [[],systematics.systematics()][int(doSys)]:
            pars = systematics.central()
            pars.update(sys)
            f = fit(signal=signal, profileVars=profile, R0_=R0_,
                    quiet=False, hackZeroBins=hackZeroBins,
                    defaults=defaults, log=log, **pars)
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
        skip=0
        Nens = 100
        mcstudy.generate(Nens,0,True)
        with open('ensemble_%s_LF%d.txt'%(ens,100*lumiFactor),'w') as ensfile:
            with open('ensemble_%s_LF%d.log'%(ens,100*lumiFactor),'w') as enslog:
                pars['log'] = enslog
                fAqq = self.central.model.w.arg('falphaL').getVal() * self.central.model.w.arg('Ac_y_ttqq').getVal()
                fAqg = self.central.model.w.arg('falphaT').getVal() * self.central.model.w.arg('Ac_y_ttqg').getVal()
                print >> ensfile, '\t'.join(["#truth", '%f'%fAqq, '%f'%fAqg])
                print >> ensfile, fit.fields(self.isAsymmetry)
                for i in range(skip,Nens):
                    alt = mcstudy.genData(i)
                    pars['label'] = 'ens%d'%i
                    f = fit(altData=alt, **pars)
                    print >> ensfile, f
                    ensfile.flush()
                    enslog.flush()


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
