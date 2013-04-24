import sys
import roo
import utils
import ROOT as r


class topModel(object):
    @roo.quiet
    def __init__(self, channelDict, asymmetry=True, quiet=False, w=None, alternate=False):

        leptons = ['el', 'mu']
        ttcomps = ('qq', 'ag', 'gg', 'qg')
        observables = ['observable', 'tridiscr']

        channels = dict((L, channelDict[(L,'top')]) for L in leptons)
        channels_qcd = dict((L + 'qcd', channelDict[(L, 'QCD')]) for L in leptons)
        gen = channelDict['gen']

        if not w: w = r.RooWorkspace('Workspace')

        for item in ['quiet', 'asymmetry', 'gen', 'channels', 'channels_qcd', 'alternate',
                     'ttcomps', 'observables', 'w']: setattr(self, item, eval(item))

        for item in ['fractions', 'xs_lumi', 'efficiencies', 'shapes', 'qcd', 'asymmetry',
                     'model', 'expressions']: getattr(self, 'import_' + item)(w)

        for item in ['d_lumi', 'd_xs_dy', 'd_xs_st']: w.arg(item).setConstant()

        for v, X in zip(observables, 'XY'):
            w.var(v).setBins(getattr(self.channels['el'].samples['data'].datas[0],
                                     'GetNbins' + X)())

    def import_fractions(self, w):
        [roo.wimport_const(w, "f_%s_hat" % comp, self.gen.samples['tt' + comp].frac)
         for comp in self.ttcomps]

        if self.asymmetry:
            f_qq_max = 0.3
            alphaL_max = 0.99 * min(chan.samples['ttqq'].alphaMax for chan in
                                    self.channels.values() +
                                    self.channels_qcd.values())
            roo.factory(w, "falphaL[0.13, -2, 2]")
            roo.factory(w, "slosh[0.3, 0, 1]")
            e_dqq = '( @0*%f + (1-@0)*abs(@1)/%f) / @2 -1'%(f_qq_max, alphaL_max)
            roo.factory(w, "expr::d_qq('%s',{slosh,falphaL,f_qq_hat})"%e_dqq)
        else:
            roo.factory(w, "d_qq[-0.999999,1]")

        roo.factory(w, "R_ag[%f,0.07,1]" % (self.gen.samples['ttag'].frac /
                                            self.gen.samples['ttqq'].frac))
        roo.factory(w, "expr::f_qq('(1+@0)*@1',{d_qq,f_qq_hat})")
        roo.factory(w, "prod::f_ag(R_ag,f_qq)")
        roo.factory(w, "expr::f_qg('(1-@0-@1)/(1+@2*@3*@4/(@5*@6))'," +\
                        "{f_qq,f_ag,R_ag,f_gg_hat,f_qq_hat,f_ag_hat,f_qg_hat})")
        roo.factory(w, "expr::f_gg('1-@0-@1-@2',{f_qq,f_ag,f_qg})")

    def import_xs_lumi(self, w):
        roo.factory(w, "d_lumi[0,-0.2,0.2]")
        for L, channel in self.channels.items() + self.channels_qcd.items():
            roo.wimport_const(w, 'lumi_%s_hat' % L, channel.lumi)
            roo.factory(w, "expr::lumi_%s('(1+@0)*@1', {d_lumi, lumi_%s_hat})" % (L, L))

        xs_constraints = dict([(samp[:2], (data.xs, data.xs_sigma))
                               for samp, data in self.channels['el'].samples.items()
                               if data.xs > 0])

        for sample, (xs, delta) in xs_constraints.items():
            roo.wimport_const(w, 'xs_%s_hat' % sample, xs)
            roo.factory(w, "d_xs_%s[0,-0.5,1.5]" % sample)
            roo.factory(w, "expr::xs_%s('(1+@0)*@1',{d_xs_%s, xs_%s_hat})" % (3 * (sample,)))

        [roo.factory(w, "prod::xs_tt%s(f_%s,xs_tt)" % (c, c)) for c in self.ttcomps]

    def import_efficiencies(self, w, channels=None):
        if not channels: channels = self.channels
        [roo.wimport_const(w, 'eff_%s_%s' % (lepton, sample), data.eff)
         for lepton, channel in channels.items()
         for sample, data in channel.samples.items()
         if sample != 'data']

    def import_shapes(self, w, channels=None):
        if not channels:
            channels = self.channels
            [roo.factory(w, "%s[-1,1]" % obs) for obs in self.observables]
            roo.factory(w, "channel[%s]" % ','.join("%s=%d" % (s, i)
                                                    for i, s in enumerate(channels)))

        [self.import_shape(w, lepton, sample, data)
         for lepton, channel in channels.items()
         for sample, data in channel.samples.items()
         if sample != 'data']

    @roo.quiet
    def import_shape(self, w, lepton, sample, data):
        name = '_'.join([lepton, sample])
        arglist = r.RooArgList(*[w.var(o) for o in self.observables])
        argset = r.RooArgSet(arglist)

        for i, label in enumerate(['both', 'symm'][
                :None if sample in ['ttag', 'ttqg', 'ttqq', 'dy'] else -1]):
            nL = (name, label)
            roo.wimport(w, r.RooDataHist('_sim_'.join(nL), '', arglist, data.datas[i]))
            roo.wimport(w, r.RooHistPdf('_'.join(nL), '', argset, w.data('_sim_'.join(nL))))

        roo.factory(w, "prod::expect_%s(%s)" %
                    (name, ','.join(['lumi_' + lepton,
                                     'xs_' + sample,
                                     'eff_' + name,
                                     '-1',
                                     'factor_%s' % lepton][:None if 'qcd' in lepton else -2])))

    def import_qcd(self, w):
        [roo.factory(w, "factor_%s[1,0,5]" % lepton) for lepton in self.channels_qcd]

        self.import_efficiencies(w, self.channels_qcd)
        self.import_shapes(w, self.channels_qcd)
        arglist = r.RooArgList(*[w.var(o) for o in self.observables])
        argset = r.RooArgSet(arglist)
        for L, channel in self.channels_qcd.items():
            hist = channel.samples['data'].datas[0]
            dhist = '%s_data_sim_both' % L
            roo.wimport(w, r.RooDataHist(dhist, '', arglist, hist))
            roo.wimport(w, r.RooHistPdf('%s_data_both' % L, '', argset, w.data(dhist)))
            roo.factory(w, "expr::expect_%s_data('@0*%f',{factor_%s})" %
                        (L, hist.Integral(), L))

    def import_asymmetry(self, w):
        if not self.asymmetry: return

        roo.factory(w, "falphaT[0.18, -0.8, 0.8]")
        roo.factory(w, "expr::alphaT('@0/@1',{falphaT,f_qg})")
        roo.factory(w, "expr::alphaL('@0/@1',{falphaL,f_qq})")

        [(roo.factory(w, "SUM::%(n)s( alphaT * %(n)s_both, %(n)s_symm )" % {'n': L + '_ttag'}),
          roo.factory(w, "SUM::%(n)s( alphaT * %(n)s_both, %(n)s_symm )" % {'n': L + '_ttqg'}),
          roo.factory(w, "SUM::%(n)s( alphaL * %(n)s_both, %(n)s_symm )" % {'n': L + '_ttqq'}))
         for L in self.channels.keys() + self.channels_qcd.keys()]

        assert self.gen.samples['tt'].datas[0].GetXaxis().GetTitle() == 'genTopDeltaBetazRel'
        assert self.gen.samples['tt'].datas[0].GetYaxis().GetTitle() == 'genTopPhiBoost'

        for n, d in self.gen.samples.items():
            roo.wimport_const(w, 'Ac_y_' + n, utils.asymmetry(d.datasX[0]))
            roo.wimport_const(w, 'Ac_phi_' + n, utils.asymmetry(d.datasY[0]))
            if not self.quiet:
                w.arg('Ac_y_' + n).Print()
                w.arg('Ac_phi_' + n).Print()

    def import_model(self, w):
        which = dict((i, '_both') for i in ['dy', 'wj', 'st', 'ttgg', 'ttqq', 'ttqg', 'ttag'])
        if self.asymmetry: which.update({'dy': '_symm', 'ttqq': '', 'ttqg': '', 'ttag': ''})

        [roo.factory(w, "SUM::model_%s( expect_%sqcd_data * %sqcd_data_both, %s )" %
                     (lepton, lepton, lepton,
                      ','.join(['expect_%s_%s * %s_%s%s' %
                                (lepton + part, key, lepton + part, key, value)
                                for key, value in which.items()
                                for part in ['', 'qcd']
                                if not (key == 'dy' and 'qcd' in part)])))
         for lepton in self.channels]

        roo.factory(w, "SIMUL::model(channel, %s)" %
                    ', '.join("%(chan)s=model_%(chan)s" %
                              {'chan': lepton} for lepton in self.channels))

    def import_expressions(self, w):
        [roo.factory(w, "sum::expect_%s_tt(%s)" %
                     (L, ','.join(['expect_%s_tt%s' %
                                   (L, f) for f in self.ttcomps])))
         for L in self.channels]

        [roo.factory(w, "sum::expect_%s_notqcd(%s)" %
                     (L, ','.join(['expect_%s_%s' % (L, s)
                                   for s in ['wj', 'st', 'ttgg', 'ttag', 'ttqg', 'ttqq']])))
         for L in self.channels_qcd]

        [roo.factory(w, "sum::expect_%(n)s_mj(expect_%(n)sqcd_data,expect_%(n)sqcd_notqcd)" %
                     {'n': L}) for L in self.channels]

    @roo.quiet
    def import_data(self):
        w = self.w
        obs_ = w.argSet(','.join(self.observables))
        obs = r.RooArgList(obs_)

        datas = [(L, r.RooDataHist('data_' + L, 'N_{obs}^{%s}' % L, obs,
                                   chan.samples['data'].datas[0]))
                 for L, chan in self.channels.items()]

        [roo.wimport(w, dat) for _, dat in datas]
        args = [r.RooFit.Import(*dat) for dat in datas]
        roo.wimport(w, r.RooDataHist('data', 'N_{obs}', obs,
                                     r.RooFit.Index(w.arg('channel')), *args))

    def plot_fracs(self, logy=False, fileName='fractions.pdf'):
        w = self.w
        origLevel = r.gErrorIgnoreLevel
        r.gErrorIgnoreLevel = r.kWarning
        c = r.TCanvas()
        c.Print(fileName + '[')
        c.SetLogy(logy)
        plt = w.var('d_qq').frame()
        plt.SetMaximum(1.0)
        for v in range(0, 94):
            val = v * 0.01 + 0.07
            plt.SetTitle("R_ag = %.3f" % val)
            w.var('R_ag').setVal(val)
            w.arg('f_gg').plotOn(plt, r.RooFit.LineWidth(1), r.RooFit.LineColor(r.kBlue))
            w.arg('f_qg').plotOn(plt, r.RooFit.LineWidth(1), r.RooFit.LineColor(r.kRed))
            w.arg('f_qq').plotOn(plt, r.RooFit.LineWidth(1), r.RooFit.LineColor(r.kGreen))
            w.arg('f_ag').plotOn(plt, r.RooFit.LineWidth(1), r.RooFit.LineColor(r.kViolet))
            plt.Draw()
            c.Print(fileName)
        r.gErrorIgnoreLevel = origLevel
        c.Print(fileName + ']')


    def print_fracs(self):
        w = self.w
        for item in \
                ['lumi_mu', 'lumi_el', 'f_gg', 'f_qg', 'f_qq', 'f_ag'] + \
                ['xs_' + i for i in self.channels['el'].samples if i != 'data']:
            print "%s: %.04f" % (item, w.arg(item).getVal())
        print

    def print_n(self):
        w = self.w
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


    def visualize(self, fileName, canvas=None):

        if not canvas: canvas = r.TCanvas()
        else: canvas.Clear()

        canvas.Divide(5,2,0,0)
        canvas.cd(1)
        
        w = self.w
        
        for j,lep in enumerate(['el','mu']):
            for i in range(5):
                canvas.cd(5*j + i + 1)
                f = w.arg('observable').frame()
                w.arg('tridiscr').setVal(-0.8 + i*0.4)
                tridiscrBinWidth = 0.4

                mod = w.pdf('model_%s'%lep)
                args = [r.RooFit.Project(r.RooArgSet()),
                        r.RooFit.Normalization(tridiscrBinWidth, r.RooAbsReal.Relative),
                        r.RooFit.LineWidth(1),
                        ]
                dargs = [f, 
                         r.RooFit.Cut('%f<=tridiscr && tridiscr<%f'%(-1+i*0.4,-0.6+i*0.4)),
                         r.RooFit.MarkerSize(0.5),
                         r.RooFit.XErrorSize(0)]
                w.data('data_%s'%lep).plotOn(*dargs)

                pf = '' if self.asymmetry else '_both'
                stack = [('_ttqq'+pf,),('_ttqg'+pf,'_ttag'+pf),('_ttgg_both',),('_wj_both',),
                         ('qcd_*',), ('_dy*','_st_both')]
                colors = [r.kViolet,r.kBlue-7,r.kBlue+2,r.kGreen+1,r.kRed,r.kGray]
                for iStack in range(len(stack)):
                    comps = lep+(','+lep).join(sum(stack[iStack:],()))
                    mod.plotOn(f,
                               r.RooFit.Components(comps),
                               r.RooFit.LineColor(colors[iStack]),
                               r.RooFit.FillColor(colors[iStack]),
                               r.RooFit.DrawOption('F'),
                               *args)

                w.data('data_%s'%lep).plotOn(*dargs)
                mod.plotOn(f, r.RooFit.LineColor(r.kOrange), *args)

                f.SetMaximum(4000)
                f.Draw()

        canvas.Print(fileName + '.pdf')
