import sys
import roo
import utils
import math
import ROOT as r
from asymmNames import genNameX,genNameY

def unqueue(h, doIt):
    return utils.unQueuedBins(h,5,[-1,1],[-1,1]) if doIt else h


class topModel(object):
    @roo.quiet
    def __init__(self, channelDict, asymmetry=True, w=None, quiet=True):

        self.fixFractions = True

        leptons = ['el', 'mu']
        ttcomps = ('qq', 'ag', 'gg', 'qg')
        observables = ['XL','XT','tridiscr'] if asymmetry else ['observable', 'tridiscr']

        channels = dict((L, channelDict[(L,'top')]) for L in leptons)
        channels_qcd = dict((L + 'qcd', channelDict[(L, 'QCD')]) for L in leptons)
        gen = channelDict['gen']

        if not w: w = r.RooWorkspace('Workspace')

        for item in ['asymmetry', 'gen', 'channels', 'channels_qcd','quiet',
                     'ttcomps', 'observables', 'w']: setattr(self, item, eval(item))

        for item in ['fractions', 'xs_lumi', 'efficiencies', 'shapes', 'qcd', 'asymmetry',
                     'model', 'expressions']: getattr(self, 'import_' + item)(w)

        for item in ['d_lumi', 'd_xs_dy', 'd_xs_st']: w.arg(item).setConstant()
        if False:
            w.arg('d_xs_tt').setVal(0.16)
            w.arg('d_xs_tt').setConstant()
            w.arg('d_xs_wj').setVal(0.78)
            w.arg('d_xs_wj').setConstant()


    def import_fractions(self, w):
        [roo.wimport_const(w, "f_%s_hat" % comp, self.gen.samples['tt' + comp].frac)
         for comp in self.ttcomps]

        if self.asymmetry:
            f_qq_max = 0.3
            alphaL_max = 0.99 * min(chan.samples['ttqq'].alphaMax for chan in
                                    self.channels.values() +
                                    self.channels_qcd.values())
            roo.factory(w, "falphaL[0.13, -2, 2]")
            roo.factory(w, "slosh[0.44, 0, 1]")
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

        if self.fixFractions:
            w.arg('slosh').setConstant(True)
            w.arg('R_ag').setConstant(True)


    def import_xs_lumi(self, w):
        roo.factory(w, "d_lumi[0,-0.2,0.2]")
        roo.factory(w, 'lumi_factor[1.0]')
        w.arg('lumi_factor').setConstant()
        for L, channel in self.channels.items() + self.channels_qcd.items():
            roo.wimport_const(w, 'lumi_%s_hat' % L, channel.lumi)
            if 'qcd' in L: roo.factory(w, "expr::lumi_%s('(1+@0)*@1', {d_lumi, lumi_%s_hat})" % (L, L))
            else:          roo.factory(w, "expr::lumi_%s('(1+@0)*@1*@2', {d_lumi, lumi_%s_hat, lumi_factor})" % (L, L))

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
            for v, X in zip(self.observables, 'XYZ'[:None if self.asymmetry else -1]):
                w.var(v).setBins(getattr(unqueue(self.channels['el'].samples['data'].datas[0],self.asymmetry),
                                         'GetNbins' + X)())

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
            roo.wimport(w, r.RooDataHist('_sim_'.join(nL), '', arglist, unqueue(data.datas[i], self.asymmetry)))
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
            hist = unqueue(channel.samples['data'].datas[0], self.asymmetry)
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
         for L in self.channels.keys()]

        [(roo.factory(w, "SUM::%(n)s( 0 * %(n)s_both, %(n)s_symm )" % {'n': L + '_ttag'}),
          roo.factory(w, "SUM::%(n)s( 0 * %(n)s_both, %(n)s_symm )" % {'n': L + '_ttqg'}),
          roo.factory(w, "SUM::%(n)s( 0 * %(n)s_both, %(n)s_symm )" % {'n': L + '_ttqq'}))
         for L in self.channels_qcd.keys()]

        assert self.gen.samples['tt'].datas[0].GetXaxis().GetTitle() == genNameX
        assert self.gen.samples['tt'].datas[0].GetYaxis().GetTitle() == genNameY

        for n, d in self.gen.samples.items():
            y,ey = utils.asymmetry(d.datasX[0])
            p,ep = utils.asymmetry(d.datasY[0])
            roo.wimport_const(w, 'Ac_y_' + n, y)
            roo.wimport_const(w, 'Ac_phi_' + n, p)
            roo.wimport_const(w, 'err_Ac_y_' + n, ey)
            roo.wimport_const(w, 'err_Ac_phi_' + n, ep)

    def import_model(self, w):
        which = dict((i, '_both') for i in ['dy', 'wj', 'st', 'ttgg', 'ttqq', 'ttqg', 'ttag'])
        if self.asymmetry: which.update({'dy': '_symm', 'ttqq': '', 'ttqg': '', 'ttag': ''})

        [roo.factory(w, "SUM::%s_mj( expect_%sqcd_data * %sqcd_data_both, %s )" %
                     (lepton, lepton, lepton,
                      ','.join(['expect_%s_%s * %s_%s%s' %
                                (lepton + 'qcd', key, lepton + 'qcd', key, value)
                                for key, value in which.items() if key != 'dy'])))
         for lepton in self.channels]

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
    def import_data(self, altData=None):
        w = self.w
        if altData:
            altData.SetName('data')
            roo.wimport(w, altData)
            return
        obs_ = w.argSet(','.join(self.observables))
        obs = r.RooArgList(obs_)

        dataHists = dict([(L,unqueue(chan.samples['data'].datas[0], self.asymmetry)) for L,chan in self.channels.items()])
        datas = [(L, r.RooDataHist('data_' + L, 'N_{obs}^{%s}' % L, obs, data))  for L, data in dataHists.items()]
        if not self.quiet:
            print self.observables[0]
            print '\n'.join('%s: %.6f (%.6f)'%((L,)+utils.asymmetry(hist.ProjectionX())) for L,hist in dataHists.items())
            print self.observables[1]
            print '\n'.join('%s: %.6f (%.6f)'%((L,)+utils.asymmetry(hist.ProjectionY())) for L,hist in dataHists.items())
            print

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

    def print_n(self,logfile):
        scale = 1000.
        w = self.w
        length = 24

        chans = ['el','mu']
        cross = ['tt', 'wj', 'mj', 'st', 'dy']
        print>>logfile, '\t','&','\t&\t'.join(r'\multicolumn{2}{c}{N_{%s}}'%xs for xs in cross), 'Total', 'Observed'
        for chan in chans:
            tot = 0
            tote2 = 0
            print>>logfile, chan,'&',
            for xs in cross:
                val = w.arg('expect_%s_%s' % (chan, xs)).getVal()
                delta = w.arg('d_xs_%s'%xs)
                factor= w.arg('factor_%sqcd'%chan)
                relerr = max(0.0000001,
                             delta.getError() / (1+delta.getVal()) if delta else
                             factor.getError() / factor.getVal() if xs=='mj' else 0)
                tot += val
                tote2 += (val*relerr)**2
                print>>logfile, utils.roundString(val/scale,relerr*val/scale).rjust(length / 3), '&',
            print>>logfile, utils.roundString(tot/scale,math.sqrt(tote2)/scale).rjust(length / 3), '&',
            print>>logfile, self.channels[chan].samples['data'].datas[0].Integral()/scale, r'\\'
        print>>logfile


    @roo.quiet
    def chi2(self,lep):
        w = self.w
        mod = w.pdf('model_%s'%lep)
        dhist = unqueue(self.channels[lep].samples['data'].datas[0], True)
        exp = mod.generateBinned(w.argSet(','.join(self.observables)), 0, True, False)
        hist = exp.createHistogram(','.join(self.observables), 5, 5, 5)
        chi2 = 5*[0]
        for iX in range(5):
            for iY in range(5):
                for iZ in range(5):
                    ibin = hist.GetBin(iX+1,iY+1,iZ+1)
                    P = hist.GetBinContent(ibin)
                    Q = dhist.GetBinContent(ibin)
                    chi2[iZ] += (P-Q)**2 / P
        del dhist
        del hist
        return sum(chi2)

    @roo.quiet
    def visualize2D(self, printName=''):
        w = self.w
        titles = ['X_{L}','X_{T}','#Delta']
        for v,t in zip(self.observables,titles) :
            w.arg(v).SetTitle(t)

        r.gROOT.ProcessLine(".L tdrstyle.C")
        r.setTDRStyle()
        r.TGaxis.SetMaxDigits(4)
        canvas = r.TCanvas()
        canvas.Print(printName+'[')

        def expected_histogram(pdfname, expect=0):
            print pdfname,expect
            mod = w.pdf(pdfname)
            if not mod:
                w.Print()
                exit()
            args = ','.join(self.observables)
            exp = mod.generateBinned(w.argSet(args), expect, True, False)
            hist = exp.createHistogram(args, 5, 5, 5)
            hist.SetName(pdfname+'genHist')
            return hist

        for j,lep in enumerate(['el','mu']):
            dhist = unqueue(self.channels[lep].samples['data'].datas[0], True)
            hist = expected_histogram('model_'+lep)

            print self.channels[lep].samples.keys()
            hists = dict([(s, expected_histogram(lep+'_'+s+('_both' if s in ['wj','st','ttgg'] else '_symm' if s=='dy' else ''), 
                                                 w.arg('expect_%s_%s'%(lep,s)).getVal()))
                     for s in self.channels[lep].samples if s not in ['data']])
            hists['mj'] = expected_histogram(lep+'_mj', w.arg('expect_%s_mj'%lep).getVal())

            stack = [('ttqq',),('ttqg','ttag'),('ttgg',),('wj',),('mj',), ('dy','st')][::-1]
            colors = [r.kViolet,r.kBlue-7,r.kBlue+2,r.kGreen+1,r.kRed,r.kGray][::-1]
            stackers = []
            for s,c in zip(stack,colors):
                h = hists[s[0]]
                for n in s[1:]: h.Add(hists[n])
                h.SetFillColor(c)
                h.SetLineColor(c)
                stackers.append(h)

            for i in range(3):
                def proj(h):
                    return getattr(h, 'Projection'+'XYZ'[i])(h.GetName()+'xyz'[i], 0, -1, 0 if i==2 else 3, -1 if i==2 else 3)
                model = proj(hist)
                data = proj(dhist)
                comps = [proj(h) for h in stackers]
                hstack = r.THStack('stack','')
                for c in comps: hstack.Add(c)
                model.SetMinimum(0)
                model.SetLineColor(r.kBlue)
                data.SetMinimum(0)
                data.Draw()
                model.Draw('hist same')
                hstack.Draw('hist same')
                data.Draw('same')
                canvas.Print(printName)
                sys.stdout.write(' ')
                if i<2:
                    _,amodel = utils.symmAnti(model)
                    _,adata = utils.symmAnti(data)
                    astackers = [utils.symmAnti(s)[1] for s in stackers]
                    amax = 1.3* max(adata.GetMaximum(),amodel.GetMaximum())
                    adata.SetMaximum(amax)
                    adata.SetMinimum(-amax)
                    amodel.SetMinimum(-amax)
                    amodel.SetMaximum(amax)
                    adata.Draw()
                    amodel.Draw('hist same')
                    for a in astackers:
                        a.Draw('hist same')
                    canvas.Print(printName)
                    sys.stdout.write(' ')
                continue

                canvas.cd(1+j*3+i)
                f = w.arg(self.observables[i]).frame()
                f.SetLineColor(r.kWhite)
                mod = w.pdf('model_%s'%lep)
                toy = mod.generateBinned(w.argSet(','.join(self.observables)), 1e7)
                args = [r.RooFit.ProjWData(toy)]
                dargs = [f,
                         r.RooFit.MarkerSize(1.1),
                         #r.RooFit.XErrorSize(0)
                         ]
                w.data('data_%s'%lep).plotOn(*dargs)

                pf = '' if self.asymmetry else '_both'
                stack = [('_ttqq'+pf,),('_ttqg'+pf,'_ttag'+pf),('_ttgg_both',),('_wj_both',),
                         ('qcd_*',), ('_dy*','_st_both')]
                colors = [r.kViolet,r.kBlue-7,r.kBlue+2,r.kGreen+1,r.kRed,r.kGray]
                for iStack in range(len(stack)):
                    sys.stdout.write(".,;:!|"[iStack])
                    sys.stdout.flush()
                    comps = lep+(','+lep).join(sum(stack[iStack:],()))
                    mod.plotOn(f,
                               r.RooFit.Components(comps),
                               r.RooFit.LineColor(colors[iStack]),
                               r.RooFit.FillColor(colors[iStack]),
                               r.RooFit.DrawOption('F'),
                               *args)

                w.data('data_%s'%lep).plotOn(*dargs)
                #mod.plotOn(f, r.RooFit.LineColor(r.kOrange), *args)

                f.Draw()
                canvas.Print(printName)
                sys.stdout.write(' ')
        print
        canvas.Print(printName+']')

        return



    def TTbarComponentsStr(self):
        '''Print a table with fractions and asymmetries.'''

        format = r'''
        \begin{tabular}{c|r@{.}lr@{.}lr@{.}l}
          \hline
          &\multicolumn{1}{c}{}&\multicolumn{2}{r}{(%s)}&\multicolumn{3}{c}{}\\
          Initial State & \multicolumn{2}{c}{Fraction} & \multicolumn{2}{c}{$\hat{A}_c^T$} & \multicolumn{2}{c}{$\hat{A}_c^L$} \\
          \hline
          \hline
          %s \\
          %s \\
          %s \\
          %s \\
          \hline
          %s \\
          \hline
        \end{tabular}'''

        def form(n,e):
            n*=100
            e*=100
            err_digit = int(math.floor(math.log(abs(e))/math.log(10))) if e else 0
            scale = float(10**(-err_digit)) if e else 1
            p,d = divmod(round(abs(n*scale))/scale,1)

            return (("%s%d&%0"+str(-err_digit)+"d(%d)")%('-' if n<0 else ' ',p,d*scale,round(e*scale))).ljust(8)

        neff = self.gen.samples['tt'].datas[0].GetEffectiveEntries()
        def frac_unc(f) : return math.sqrt(f*(1-f)*neff) / neff

        labels = [r'$\Pq\Paq$', r'$\Paq\Pg$', r'$\Pg\Pg$', '$\Pq\Pg$','$\Pp\Pp$']
        rows = dict((comp,
                     '  &  '.join( [label.rjust(10),
                                    form(self.w.arg('f_%s_hat'%comp).getVal() if comp else 1, frac_unc(self.w.arg('f_%s_hat'%comp).getVal()) if comp else 0),
                                    form(self.w.arg('Ac_phi_tt%s'%comp).getVal(),self.w.arg('err_Ac_phi_tt%s'%comp).getVal()),
                                    form(self.w.arg('Ac_y_tt%s'%comp).getVal(),self.w.arg('err_Ac_y_tt%s'%comp).getVal())
                                    ])) for label,comp in zip(labels,self.ttcomps+('',)))
        
        return format%tuple([r'\%']+[rows[i] for i in ['gg','qq','qg','ag','']])
